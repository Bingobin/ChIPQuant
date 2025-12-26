import pysam
import os

class BAMProcessor:
    def __init__(self, bam_path):
        self.bam_path = bam_path
        if not os.path.exists(bam_path):
            raise FileNotFoundError(f"BAM file not found: {bam_path}")
        
        # Check if index exists, warn if not found
        if not os.path.exists(bam_path + ".bai") and not os.path.exists(bam_path[:-4] + ".bai"):
             print(f"Warning: BAM index not found for {bam_path}. Random access may fail.", file=sys.stderr)

        self.bam = pysam.AlignmentFile(bam_path, "rb")
        self.total_mapped_reads = self._get_total_mapped()

    def _get_total_mapped(self):
        """Get total mapped reads in BAM file."""
        return self.bam.mapped

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
                # Return 0 if chromosome not found
                return 0
                
        try:
            # count() method is more efficient than iterating with fetch()
            return self.bam.count(contig=chrom, start=start-1, stop=end)
        except ValueError:
            return 0

    def close(self):
        self.bam.close()
