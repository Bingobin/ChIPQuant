import sys
import multiprocessing

from .bam_processor import BAMProcessor
from .utils import get_promoter_region


def _merge_intervals(intervals):
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = []
    for start, end in intervals:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            if end > merged[-1][1]:
                merged[-1][1] = end
    return [(start, end) for start, end in merged]


def build_intervals_by_chrom(genes, p_up, p_down, bam_references=None):
    ref_set = set(bam_references) if bam_references else None

    def normalize_chrom(chrom):
        if not ref_set:
            return chrom
        if chrom in ref_set:
            return chrom
        if chrom.startswith("chr") and chrom[3:] in ref_set:
            return chrom[3:]
        if "chr" + chrom in ref_set:
            return "chr" + chrom
        return None

    gene_body = {}
    promoter = {}
    for gene in genes:
        chrom = normalize_chrom(gene["chrom"])
        if not chrom:
            continue
        gene_body.setdefault(chrom, []).append((gene["start"] - 1, gene["end"]))
        p_start, p_end = get_promoter_region(
            gene["chrom"],
            gene["start"],
            gene["end"],
            gene["strand"],
            p_up,
            p_down,
        )
        promoter.setdefault(chrom, []).append((p_start - 1, p_end))
    gene_body = {chrom: _merge_intervals(iv) for chrom, iv in gene_body.items()}
    promoter = {chrom: _merge_intervals(iv) for chrom, iv in promoter.items()}
    return gene_body, promoter


def build_bed_intervals_by_chrom(regions, bam_references=None):
    ref_set = set(bam_references) if bam_references else None

    def normalize_chrom(chrom):
        if not ref_set:
            return chrom
        if chrom in ref_set:
            return chrom
        if chrom.startswith("chr") and chrom[3:] in ref_set:
            return chrom[3:]
        if "chr" + chrom in ref_set:
            return "chr" + chrom
        return None

    bed = {}
    for region in regions:
        chrom = normalize_chrom(region["chrom"])
        if not chrom:
            continue
        bed.setdefault(chrom, []).append((region["start"], region["end"]))
    bed = {chrom: _merge_intervals(iv) for chrom, iv in bed.items()}
    return bed


def start_summary_worker(args, genes):
    if not args.summary_out:
        return None, None
    parent_conn, child_conn = multiprocessing.Pipe(duplex=False)
    proc = multiprocessing.Process(
        target=_summary_worker,
        args=(args, genes, child_conn),
    )
    proc.start()
    return parent_conn, proc


def start_bed_summary_worker(args, regions):
    if not args.summary_out:
        return None, None
    parent_conn, child_conn = multiprocessing.Pipe(duplex=False)
    proc = multiprocessing.Process(
        target=_bed_summary_worker,
        args=(args, regions, child_conn),
    )
    proc.start()
    return parent_conn, proc


def _summary_worker(args, genes, conn):
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
        gene_body_intervals, promoter_intervals = build_intervals_by_chrom(
            genes,
            args.promoter_upstream,
            args.promoter_downstream,
            bam_references=bam_proc.bam.references,
        )
        if not promoter_intervals:
            print(
                "Warning: promoter intervals are empty after chromosome normalization.",
                file=sys.stderr,
            )
        total_fragments = bam_proc.total_mapped_reads
        gene_body_fragments = bam_proc.count_fragments_overlapping(gene_body_intervals)
        promoter_fragments = bam_proc.count_fragments_overlapping(promoter_intervals)
        bam_proc.close()
        conn.send((total_fragments, gene_body_fragments, promoter_fragments, None))
    except Exception as e:
        conn.send((None, None, None, str(e)))
    finally:
        conn.close()


def _bed_summary_worker(args, regions, conn):
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
        bed_intervals = build_bed_intervals_by_chrom(
            regions,
            bam_references=bam_proc.bam.references,
        )
        total_fragments = bam_proc.total_mapped_reads
        bed_fragments = bam_proc.count_fragments_overlapping(bed_intervals)
        bam_proc.close()
        conn.send((total_fragments, bed_fragments, None))
    except Exception as e:
        conn.send((None, None, str(e)))
    finally:
        conn.close()
