# ChIPQuant

**ChIPQuant** is a simplified, modular Python tool for calculating read counts and CPM (Counts Per Million) from ChIP-seq data (BAM) using either gene annotations (GTF) or custom regions (BED).

## Features

- **Gene Body Quantification**: Counts reads mapped to the full gene span (GTF mode).
- **Promoter Quantification**: Counts reads within a user-defined window around the Transcription Start Site (TSS) (GTF mode).
- **Custom BED Quantification**: Counts reads over user-provided BED regions.
- **CPM Normalization**: Automatically calculates normalized CPM values.
- **Modular Codebase**: Clean, maintainable Python project structure.

## Requirements

- Python 3.6+
- `pysam`

## Installation

1. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

Execute the tool from the project root:

```bash
python chipquant.py --gtf /path/to/genes.gtf \
                    --bam /path/to/chipseq.bam \
                    --out counts_output.tsv \
                    --promoter_upstream 2000 \
                    --promoter_downstream 2000
```

Or, provide custom regions:

```bash
python chipquant.py --bed /path/to/regions.bed \
                    --bam /path/to/chipseq.bam \
                    --out bed_counts.tsv
```

### BED Format Recommendation

Use standard BED (0-based, half-open). Minimal required columns are `chrom`, `start`, `end`.

Example:

```text
chr1	10000	10500	peak_1	0	+
chr1	20000	21000	peak_2	0	-
chr2	5000	8000	peak_3	0	.
```

### Parameters

- `--gtf`: Path to the input GTF annotation file.
- `--bed`: Path to the input BED file (0-based, half-open).
- `--bam`: Path to the indexed BAM file.
- `--out`: Path for the output TSV file.
- `--promoter_upstream`: Base pairs upstream of TSS (default: 2000, GTF only).
- `--promoter_downstream`: Base pairs downstream of TSS (default: 2000, GTF only).
- `--threads` / `-t`: Number of worker processes (default: 1).
- `--min_mapq`: Minimum MAPQ to count a read (default: 30).
- `--include-secondary`: Include secondary alignments (default: false).
- `--include-supplementary`: Include supplementary alignments (default: false).
- `--include-duplicates`: Include duplicate reads (default: false).
- `--allow-unpaired`: Allow unpaired reads to be counted (default: false).
- `--allow-improper-pair`: Allow improperly paired reads to be counted (default: false).
- `--count-all-paired-reads`: Count both mates instead of fragment-level counting (default: false).
- `--cpm-denominator`: CPM denominator strategy: `filtered` (default) or `mapped`.
- `--summary-out`: Write summary metrics to this TSV path (default: disabled).

## Output Format

The output TSV contains (GTF mode):
- Gene Metadata (ID, Symbol, Type)
- **GeneBody_Count**: Raw reads in gene body.
- **GeneBody_CPM**: CPM normalized counts for gene body.
- **Promoter_Count**: Raw reads in promoter region.
- **Promoter_CPM**: CPM normalized counts for promoter region.

The output TSV contains (BED mode):
- Region metadata (RegionID, Chrom, Start, End, Strand)
- **Count**: Raw reads in region.
- **CPM**: CPM normalized counts for region.

## Summary Output

When `--summary-out` is provided, a summary TSV is written containing:
- **Total_Mapped_Fragments**: Total fragments passing filters across the BAM.
- **GeneBody_Overlap_Fragments**: Unique fragments overlapping any gene body region (GTF mode).
- **Promoter_Overlap_Fragments**: Unique fragments overlapping any promoter region (GTF mode).
- **BED_Overlap_Fragments**: Unique fragments overlapping any BED region (BED mode).
