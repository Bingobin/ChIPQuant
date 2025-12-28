# ChIPQuant

**ChIPQuant** is a simplified, modular Python tool for calculating **Gene Body** and **Promoter** read counts and CPM (Counts Per Million) from ChIP-seq data (BAM) and gene annotations (GTF).

## Features

- **Gene Body Quantification**: Counts reads mapped to the full gene span.
- **Promoter Quantification**: Counts reads within a user-defined window around the Transcription Start Site (TSS).
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

Execute the package from the project root:

```bash
python -m src.main --gtf /path/to/genes.gtf \
                   --bam /path/to/chipseq.bam \
                   --out counts_output.tsv \
                   --promoter_upstream 2000 \
                   --promoter_downstream 2000
```

### Parameters

- `--gtf`: Path to the input GTF annotation file.
- `--bam`: Path to the indexed BAM file.
- `--out`: Path for the output TSV file.
- `--promoter_upstream`: Base pairs upstream of TSS (default: 2000).
- `--promoter_downstream`: Base pairs downstream of TSS (default: 2000).
- `--threads` / `-t`: Number of worker processes (default: 1).
- `--min_mapq`: Minimum MAPQ to count a read (default: 30).
- `--include-secondary`: Include secondary alignments (default: false).
- `--include-supplementary`: Include supplementary alignments (default: false).
- `--include-duplicates`: Include duplicate reads (default: false).
- `--allow-unpaired`: Allow unpaired reads to be counted (default: false).
- `--allow-improper-pair`: Allow improperly paired reads to be counted (default: false).
- `--count-all-paired-reads`: Count both mates instead of fragment-level counting (default: false).
- `--cpm-denominator`: CPM denominator strategy: `filtered` (default) or `mapped`.

## Output Format

The output TSV contains:
- Gene Metadata (ID, Symbol, Type, Position)
- **GeneBody_Count**: Raw reads in gene body.
- **GeneBody_CPM**: CPM normalized counts for gene body.
- **Promoter_Count**: Raw reads in promoter region.
- **Promoter_CPM**: CPM normalized counts for promoter region.

## License

MIT
