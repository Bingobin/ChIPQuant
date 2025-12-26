import sys

class GTFParser:
    def __init__(self, gtf_path):
        self.gtf_path = gtf_path

    def _parse_attributes(self, attr_str):
        """Parse GTF attribute string into a dictionary."""
        d = {}
        for item in attr_str.strip().split(";"):
            item = item.strip()
            if not item or " " not in item:
                continue
            # Compatibility for GTF attributes separated by space
            parts = item.split(" ", 1)
            if len(parts) == 2:
                k, v = parts
                d[k] = v.strip().strip('"')
        return d

    def get_genes(self):
        """
        Generator that parses the GTF file and returns gene information.
        Returns dictionary format:
        {
            'gene_id': str,
            'gene_name': str,
            'gene_type': str,
            'chrom': str,
            'start': int,
            'end': int, # 1-based, inclusive as per GTF usually, but python is 0-based. 
                        # GTF is 1-based. We will keep raw coordinates and handle conversion when using pysam.
            'strand': str
        }
        """
        genes_seen = set()
        
        with open(self.gtf_path, "r") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                
                chrom, source, feature, start, end, score, strand, frame, attrs = parts
                
                # We only care about 'gene' feature lines. If there are multiple transcripts, 
                # GTF usually has a gene line defining the entire gene boundary.
                # If there is no gene line, we would need to aggregate exon/transcripts.
                # For simplicity and adhering to the original script logic, we prioritize feature == 'gene'.
                # Robustness: currently only handling feature == 'gene'.
                if feature != "gene":
                    continue
                
                a = self._parse_attributes(attrs)
                
                gene_id = a.get("gene_id") or a.get("geneId") or a.get("gene")
                if not gene_id:
                    continue
                
                # Deduplicate (Although usually gene lines are unique per gene_id)
                if gene_id in genes_seen:
                    continue
                genes_seen.add(gene_id)
                
                gene_name = a.get("gene_name") or a.get("Name") or gene_id
                gene_type = a.get("gene_type") or a.get("gene_biotype") or "unknown"
                
                yield {
                    'gene_id': gene_id,
                    'gene_name': gene_name,
                    'gene_type': gene_type,
                    'chrom': chrom,
                    'start': int(start),
                    'end': int(end),
                    'strand': strand
                }
