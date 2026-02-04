import argparse
import sys
import csv
import multiprocessing
from .gtf_parser import GTFParser
from .bed_parser import BEDParser
from .bam_processor import BAMProcessor
from .quantifier import Quantifier, quantify_gene_core
from .summary import start_summary_worker, start_bed_summary_worker

def parse_args():
    ap = argparse.ArgumentParser(
        description="Calculate read counts/CPM from ChIP-seq BAM and GTF/BED."
    )
    input_group = ap.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--gtf", help="GTF annotation file.")
    input_group.add_argument("--bed", help="BED regions file (0-based, half-open).")
    ap.add_argument("--bam", required=True, help="ChIP-seq BAM file (indexed).")
    ap.add_argument("--out", required=True, help="Output TSV file path.")
    ap.add_argument("--promoter_upstream", type=int, default=2000, 
                    help="Promoter upstream distance from TSS (default: 2000, GTF only)")
    ap.add_argument("--promoter_downstream", type=int, default=2000, 
                    help="Promoter downstream distance from TSS (default: 2000, GTF only)")
    ap.add_argument("--threads", "-t", type=int, default=1,
                    help="Number of threads (processes) to use (default: 1)")
    ap.add_argument("--min_mapq", type=int, default=30,
                    help="Minimum MAPQ to count a read (default: 30)")
    ap.add_argument("--include-secondary", action="store_true",
                    help="Include secondary alignments in counts")
    ap.add_argument("--include-supplementary", action="store_true",
                    help="Include supplementary alignments in counts")
    ap.add_argument("--include-duplicates", action="store_true",
                    help="Include duplicate reads in counts")
    ap.add_argument("--allow-unpaired", action="store_true",
                    help="Allow unpaired reads to be counted")
    ap.add_argument("--allow-improper-pair", action="store_true",
                    help="Allow improperly paired reads to be counted")
    ap.add_argument("--count-all-paired-reads", action="store_true",
                    help="Count both mates instead of read1 only")
    ap.add_argument("--cpm-denominator", choices=["mapped", "filtered"], default="filtered",
                    help="CPM denominator: total mapped reads or reads passing filters")
    ap.add_argument("--summary-out",
                    help="Write summary metrics to this TSV path (default: disabled)")
    return ap.parse_args()

# Global variable for worker processes
_global_bam_processor = None

def init_worker(
    bam_path,
    min_mapq,
    include_secondary,
    include_supplementary,
    include_duplicates,
    require_paired,
    require_proper_pair,
    count_fragments,
    cpm_denominator,
):
    """
    Initialize BAMProcessor in each worker process.
    This runs once per process.
    """
    global _global_bam_processor
    try:
        _global_bam_processor = BAMProcessor(
            bam_path,
            min_mapq=min_mapq,
            include_secondary=include_secondary,
            include_supplementary=include_supplementary,
            include_duplicates=include_duplicates,
            require_paired=require_paired,
            require_proper_pair=require_proper_pair,
            count_fragments=count_fragments,
            cpm_denominator=cpm_denominator,
        )
    except Exception as e:
        print(f"Error initializing worker with BAM {bam_path}: {e}", file=sys.stderr)
        sys.exit(1)

def process_gene_wrapper(args):
    """
    Wrapper function to unpack arguments and call quantify_gene_core
    using the global BAMProcessor instance.
    """
    gene_info, p_up, p_down = args
    return quantify_gene_core(gene_info, _global_bam_processor, p_up, p_down)

def process_bed_wrapper(region):
    """
    Wrapper to process a single BED region using global BAMProcessor.
    """
    chrom = region["chrom"]
    start = region["start"]
    end = region["end"]
    count = _global_bam_processor.count_reads(chrom, start + 1, end)
    total_reads = _global_bam_processor.total_mapped_reads
    cpm = (count / total_reads * 1_000_000) if total_reads > 0 else 0.0
    return {
        "region_id": region["name"],
        "chrom": chrom,
        "start": start,
        "end": end,
        "strand": region.get("strand", "."),
        "count": count,
        "cpm": cpm,
    }

def main():
    args = parse_args()
    
    count = 0
    is_bed = args.bed is not None
    if is_bed:
        headers = [
            "RegionID",
            "Chrom",
            "Start",
            "End",
            "Strand",
            "Count",
            "CPM",
        ]
        parser = BEDParser(args.bed)
        print(f"Processing regions from {args.bed}...", file=sys.stderr)
        regions = list(parser.get_regions())
        summary_conn, summary_proc = start_bed_summary_worker(args, regions)
    else:
        headers = [
            "GeneID",
            "Symbol",
            "GeneType",
            "GeneBody_Count",
            "GeneBody_CPM",
            "Promoter_Count",
            "Promoter_CPM",
        ]
        parser = GTFParser(args.gtf)
        print(f"Processing genes from {args.gtf}...", file=sys.stderr)
        genes = list(parser.get_genes())
        summary_conn, summary_proc = start_summary_worker(args, genes)

    with open(args.out, "w", newline="") as f_out:
        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(headers)

        if args.threads > 1:
            print(f"Running in parallel with {args.threads} threads...", file=sys.stderr)

            with multiprocessing.Pool(
                processes=args.threads,
                initializer=init_worker,
                initargs=(
                    args.bam,
                    args.min_mapq,
                    args.include_secondary,
                    args.include_supplementary,
                    args.include_duplicates,
                    not args.allow_unpaired,
                    not args.allow_improper_pair,
                    not args.count_all_paired_reads,
                    args.cpm_denominator,
                ),
            ) as pool:
                if is_bed:
                    tasks = (region for region in regions)
                    for res in pool.imap(process_bed_wrapper, tasks, chunksize=100):
                        row = [
                            res["region_id"],
                            res["chrom"],
                            res["start"],
                            res["end"],
                            res["strand"],
                            res["count"],
                            f"{res['cpm']:.4f}",
                        ]
                        writer.writerow(row)
                        count += 1
                        if count % 1000 == 0:
                            print(f"Processed {count} regions...", end="\r", file=sys.stderr)
                else:
                    genes_iter = iter(genes)
                    tasks = ((gene, args.promoter_upstream, args.promoter_downstream) for gene in genes_iter)
                    for res in pool.imap(process_gene_wrapper, tasks, chunksize=100):
                        row = [
                            res['gene_id'],
                            res['gene_name'],
                            res['gene_type'],
                            res['gene_body_count'],
                            f"{res['gene_body_cpm']:.4f}",
                            res['promoter_count'],
                            f"{res['promoter_cpm']:.4f}",
                        ]
                        writer.writerow(row)
                        count += 1
                        if count % 1000 == 0:
                            print(f"Processed {count} genes...", end="\r", file=sys.stderr)

        else:
            # Single-threaded mode
            print(f"Loading BAM: {args.bam}", file=sys.stderr)
            try:
                bam_proc = BAMProcessor(
                    args.bam,
                    min_mapq=args.min_mapq,
                    include_secondary=args.include_secondary,
                    include_supplementary=args.include_supplementary,
                    include_duplicates=args.include_duplicates,
                    require_paired=not args.allow_unpaired,
                    require_proper_pair=not args.allow_improper_pair,
                    count_fragments=not args.count_all_paired_reads,
                    cpm_denominator=args.cpm_denominator,
                )
            except Exception as e:
                print(f"Error opening BAM file: {e}", file=sys.stderr)
                sys.exit(1)

            print(f"Total Mapped Reads: {bam_proc.total_mapped_reads}", file=sys.stderr)

            if is_bed:
                for region in regions:
                    count_reads = bam_proc.count_reads(
                        region["chrom"], region["start"] + 1, region["end"]
                    )
                    cpm = (
                        count_reads / bam_proc.total_mapped_reads * 1_000_000
                        if bam_proc.total_mapped_reads > 0
                        else 0.0
                    )
                    row = [
                        region["name"],
                        region["chrom"],
                        region["start"],
                        region["end"],
                        region.get("strand", "."),
                        count_reads,
                        f"{cpm:.4f}",
                    ]
                    writer.writerow(row)
                    count += 1
                    if count % 1000 == 0:
                        print(f"Processed {count} regions...", end="\r", file=sys.stderr)
            else:
                genes_iter = iter(genes)
                quantifier = Quantifier(
                    bam_proc,
                    promoter_upstream=args.promoter_upstream,
                    promoter_downstream=args.promoter_downstream,
                )
                for gene_info in genes_iter:
                    res = quantifier.process_gene(gene_info)

                    row = [
                        res['gene_id'],
                        res['gene_name'],
                        res['gene_type'],
                        res['gene_body_count'],
                        f"{res['gene_body_cpm']:.4f}",
                        res['promoter_count'],
                        f"{res['promoter_cpm']:.4f}",
                    ]
                    writer.writerow(row)
                    count += 1

                    if count % 1000 == 0:
                        print(f"Processed {count} genes...", end="\r", file=sys.stderr)

            bam_proc.close()

    if args.summary_out:
        print("Waiting for summary...", file=sys.stderr)
        if is_bed:
            total_fragments, bed_fragments, err = summary_conn.recv()
            summary_proc.join()
            if err:
                print(f"Summary failed: {err}", file=sys.stderr)
            else:
                with open(args.summary_out, "w", newline="") as f_sum:
                    writer = csv.writer(f_sum, delimiter="\t")
                    writer.writerow(["Metric", "Value"])
                    writer.writerow(["Total_Mapped_Fragments", total_fragments])
                    writer.writerow(["BED_Overlap_Fragments", bed_fragments])
                print(
                    f"\n[Done] Processed {count} regions. Results saved to {args.out}. "
                    f"Summary saved to {args.summary_out}",
                    file=sys.stderr,
                )
        else:
            total_fragments, gene_body_fragments, promoter_fragments, err = summary_conn.recv()
            summary_proc.join()
            if err:
                print(f"Summary failed: {err}", file=sys.stderr)
            else:
                with open(args.summary_out, "w", newline="") as f_sum:
                    writer = csv.writer(f_sum, delimiter="\t")
                    writer.writerow(["Metric", "Value"])
                    writer.writerow(["Total_Mapped_Fragments", total_fragments])
                    writer.writerow(["GeneBody_Overlap_Fragments", gene_body_fragments])
                    writer.writerow(["Promoter_Overlap_Fragments", promoter_fragments])
                print(
                    f"\n[Done] Processed {count} genes. Results saved to {args.out}. "
                    f"Summary saved to {args.summary_out}",
                    file=sys.stderr,
                )
    else:
        if is_bed:
            print(
                f"\n[Done] Processed {count} regions. Results saved to {args.out}",
                file=sys.stderr,
            )
        else:
            print(
                f"\n[Done] Processed {count} genes. Results saved to {args.out}",
                file=sys.stderr,
            )

if __name__ == "__main__":
    main()
