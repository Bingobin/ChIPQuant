class BEDParser:
    def __init__(self, bed_path):
        self.bed_path = bed_path

    def get_regions(self):
        """
        Generator that parses a BED file and yields region info.
        BED is 0-based, half-open: [start, end)
        Returns dict:
        {
            'chrom': str,
            'start': int,  # 0-based
            'end': int,    # 0-based, exclusive
            'name': str,
            'score': str,
            'strand': str
        }
        """
        with open(self.bed_path, "r") as fh:
            for line in fh:
                if not line.strip():
                    continue
                if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                chrom = parts[0]
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    continue
                if end <= start:
                    continue
                name = parts[3] if len(parts) >= 4 and parts[3] else f"{chrom}:{start}-{end}"
                score = parts[4] if len(parts) >= 5 else ""
                strand = parts[5] if len(parts) >= 6 else "."
                yield {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "name": name,
                    "score": score,
                    "strand": strand,
                }
