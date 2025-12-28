import argparse
import sys
import csv
import multiprocessing
from .gtf_parser import GTFParser
from .bam_processor import BAMProcessor
from .quantifier import Quantifier
# Need to import quantify_gene_core only inside worker or use it if available
# But it's in .quantifier which we imported Quantifier from. 
# Better to import specifically if needed or just use module path if imported.
# Actually, to use quantify_gene_core in wrapper, we might need to import it inside wrapper or at top level if it's available.
# Let's import it at top level.
from .quantifier import Quantifier, quantify_gene_core

def parse_args():
    ap = argparse.ArgumentParser(
        description="Calculate Gene Body and Promoter Counts/CPM from ChIP-seq BAM and GTF."
    )
    ap.add_argument("--gtf", required=True, help="GTF annotation file.")
    ap.add_argument("--bam", required=True, help="ChIP-seq BAM file (indexed).")
    ap.add_argument("--out", required=True, help="Output TSV file path.")
    ap.add_argument("--promoter_upstream", type=int, default=2000, 
                    help="Promoter upstream distance from TSS (default: 2000)")
    ap.add_argument("--promoter_downstream", type=int, default=2000, 
                    help="Promoter downstream distance from TSS (default: 2000)")
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

def main():
    args = parse_args()
    
    # Headers for output
    headers = [
        "GeneID",
        "Symbol",
        "GeneType",
        "GeneBody_Count",
        "GeneBody_CPM",
        "Promoter_Count",
        "Promoter_CPM",
    ]
    
    count = 0
    parser = GTFParser(args.gtf)
    print(f"Processing genes from {args.gtf}...", file=sys.stderr)

    with open(args.out, "w", newline="") as f_out:
        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(headers)
        
        # Generator for gene info
        genes_iter = parser.get_genes()

        if args.threads > 1:
            print(f"Running in parallel with {args.threads} threads...", file=sys.stderr)
            
            # Prepare arguments generator
            # Each item is (gene_info, p_up, p_down)
            tasks = ((gene, args.promoter_upstream, args.promoter_downstream) for gene in genes_iter)
            
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
                # Use imap to report progress and streaming results
                # chunksize can be tuned, e.g. 100
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
            
            quantifier = Quantifier(
                bam_proc, 
                promoter_upstream=args.promoter_upstream, 
                promoter_downstream=args.promoter_downstream
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

    print(f"\n[Done] Processed {count} genes. Results saved to {args.out}", file=sys.stderr)

if __name__ == "__main__":
    main()
