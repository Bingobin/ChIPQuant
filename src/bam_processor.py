import pysam
import os
import sys

class BAMProcessor:
    def __init__(
        self,
        bam_path,
        min_mapq=30,
        include_secondary=False,
        include_supplementary=False,
        include_duplicates=False,
        require_paired=True,
        require_proper_pair=True,
        count_fragments=True,
        cpm_denominator="filtered",
    ):
        self.bam_path = bam_path
        if not os.path.exists(bam_path):
            raise FileNotFoundError(f"BAM file not found: {bam_path}")
        
        # Check if index exists, warn if not found
        if not os.path.exists(bam_path + ".bai") and not os.path.exists(bam_path[:-4] + ".bai"):
            print(
                f"Warning: BAM index not found for {bam_path}. Random access may fail.",
                file=sys.stderr,
            )

        self.bam = pysam.AlignmentFile(bam_path, "rb")
        self.min_mapq = min_mapq
        self.include_secondary = include_secondary
        self.include_supplementary = include_supplementary
        self.include_duplicates = include_duplicates
        self.require_paired = require_paired
        self.require_proper_pair = require_proper_pair
        self.count_fragments = count_fragments
        self.cpm_denominator = cpm_denominator
        self._warned_missing_chroms = set()
        self.total_mapped_reads = self._get_total_mapped()

    def _get_total_mapped(self):
        """Get total mapped reads in BAM file."""
        if self.cpm_denominator == "filtered":
            return self._count_all_reads_with_filters()
        if self._filters_active():
            print(
                "Warning: CPM denominator uses total mapped reads without filters. "
                "Use --cpm-denominator filtered to match count filters.",
                file=sys.stderr,
            )
        return self.bam.mapped

    def _filters_active(self):
        return (
            self.min_mapq > 0
            or not self.include_secondary
            or not self.include_supplementary
            or not self.include_duplicates
            or self.require_paired
            or self.require_proper_pair
            or self.count_fragments
        )

    def _should_count_read(self, read):
        if read.is_unmapped:
            return False
        if read.mapping_quality < self.min_mapq:
            return False
        if read.is_secondary and not self.include_secondary:
            return False
        if read.is_supplementary and not self.include_supplementary:
            return False
        if read.is_duplicate and not self.include_duplicates:
            return False
        if self.require_paired and not read.is_paired:
            return False
        if self.require_proper_pair and not read.is_proper_pair:
            return False
        return True

    def _count_all_reads_with_filters(self):
        count = 0
        bam = pysam.AlignmentFile(self.bam_path, "rb")
        try:
            for read in bam.fetch(until_eof=True):
                if self.count_fragments:
                    if not self._should_count_read(read):
                        continue
                    if read.is_paired:
                        if (
                            read.is_read1
                            and read.reference_id == read.next_reference_id
                            and read.template_length != 0
                        ):
                            count += 1
                    else:
                        count += 1
                elif self._should_count_read(read):
                    count += 1
        finally:
            bam.close()
        return count

    def count_reads(self, chrom, start, end):
        """
        Count reads within the specified region.
        GTF coordinates are usually 1-based, inclusive.
        pysam/bam usually uses 0-based, half-open.
        
        Args:
            chrom (str): Chromosome name
            start (int): Start position (1-based from GTF)
            end (int): End position (1-based from GTF, inclusive)
        """
        # Convert coordinates: 1-based start -> 0-based start (start - 1)
        # 1-based end (inclusive) -> 0-based end (exclusive) (end)
        # So we use pysam map/fetch(start-1, end)
        
        # Handle chromosome name mismatch (e.g., chr1 vs 1)
        if chrom not in self.bam.references:
            if chrom.startswith("chr") and chrom[3:] in self.bam.references:
                chrom = chrom[3:]
            elif "chr" + chrom in self.bam.references:
                chrom = "chr" + chrom
            else:
                if chrom not in self._warned_missing_chroms:
                    self._warned_missing_chroms.add(chrom)
                    print(
                        f"Warning: chromosome {chrom} not found in BAM references; returning 0.",
                        file=sys.stderr,
                    )
                # Return 0 if chromosome not found
                return 0
                
        try:
            # count() method is more efficient than iterating with fetch()
            if self._filters_active():
                count = 0
                if self.count_fragments:
                    counted_qnames = set()
                    region_start = start - 1
                    region_end = end
                    for read in self.bam.fetch(contig=chrom, start=region_start, stop=region_end):
                        if not self._should_count_read(read):
                            continue
                        if read.is_paired and read.reference_id != read.next_reference_id:
                            continue
                        if read.query_name in counted_qnames:
                            continue
                        if read.is_paired:
                            tlen = abs(read.template_length)
                            if tlen == 0:
                                continue
                            fragment_start = min(read.reference_start, read.next_reference_start)
                            fragment_end = fragment_start + tlen
                        else:
                            fragment_start = read.reference_start
                            fragment_end = read.reference_end
                        if fragment_end > region_start and fragment_start < region_end:
                            counted_qnames.add(read.query_name)
                            count += 1
                else:
                    for read in self.bam.fetch(contig=chrom, start=start-1, stop=end):
                        if self._should_count_read(read):
                            count += 1
                return count
            return self.bam.count(contig=chrom, start=start-1, stop=end)
        except ValueError:
            return 0

    def close(self):
        self.bam.close()

    def iter_fragments(self):
        """
        Yield fragment intervals as (chrom, start, end) in 0-based, half-open coordinates.
        """
        for read in self.bam.fetch(until_eof=True):
            if not self._should_count_read(read):
                continue
            if read.is_paired:
                if not read.is_read1:
                    continue
                if read.reference_id != read.next_reference_id:
                    continue
                tlen = abs(read.template_length)
                if tlen == 0:
                    continue
                start = min(read.reference_start, read.next_reference_start)
                end = start + tlen
            else:
                start = read.reference_start
                end = read.reference_end
            yield read.reference_name, start, end

    def count_fragments_overlapping(self, intervals_by_chrom):
        count = 0
        for chrom, start, end in self.iter_fragments():
            intervals = intervals_by_chrom.get(chrom)
            if not intervals:
                continue
            if _intervals_overlap(intervals, start, end):
                count += 1
        return count

    def count_filtered_fragments(self):
        if not self.count_fragments:
            return self._count_all_reads_with_filters()
        count = 0
        for _ in self.iter_fragments():
            count += 1
        return count

def _intervals_overlap(intervals, start, end):
    """
    intervals: list of (start, end) sorted, non-overlapping, 0-based half-open.
    """
    lo = 0
    hi = len(intervals) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        iv_start, iv_end = intervals[mid]
        if end <= iv_start:
            hi = mid - 1
        elif start >= iv_end:
            lo = mid + 1
        else:
            return True
    return False
