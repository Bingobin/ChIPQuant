from .utils import calculate_cpm, get_promoter_region

class Quantifier:
    def __init__(self, bam_processor, promoter_upstream=2000, promoter_downstream=2000):
        self.bam_processor = bam_processor
        self.p_up = promoter_upstream
        self.p_down = promoter_downstream

    def process_gene(self, gene_info):
        """
        Process a single gene: Calculate Gene Body and Promoter Counts/CPM.
        
        Args:
            gene_info (dict): Dictionary yielded from GTFParser
            
        Returns:
            dict: Original info plus count and cpm results
        """
        return quantify_gene_core(gene_info, self.bam_processor, self.p_up, self.p_down)

def quantify_gene_core(gene_info, bam_processor, p_up, p_down):
    """
    Standalone function for gene quantification logic.
    Useful for multiprocessing where the Quantifier instance cannot be pickled easily
    due to bam_processor (pysam object).
    """
    # 1. Gene Body Quantification
    # Gene Body defined as Start to End
    gb_count = bam_processor.count_reads(
        gene_info['chrom'], 
        gene_info['start'], 
        gene_info['end']
    )
    
    # 2. Promoter Quantification
    p_start, p_end = get_promoter_region(
        gene_info['chrom'],
        gene_info['start'],
        gene_info['end'],
        gene_info['strand'],
        p_up,
        p_down
    )
    
    promoter_count = bam_processor.count_reads(
        gene_info['chrom'],
        p_start,
        p_end
    )
    
    # 3. Calculate CPM
    total_reads = bam_processor.total_mapped_reads
    gb_cpm = calculate_cpm(gb_count, total_reads)
    promoter_cpm = calculate_cpm(promoter_count, total_reads)
    
    # 4. Update result
    result = gene_info.copy()
    result.update({
        'gene_body_count': gb_count,
        'gene_body_cpm': gb_cpm,
        'promoter_start': p_start,
        'promoter_end': p_end,
        'promoter_count': promoter_count,
        'promoter_cpm': promoter_cpm
    })
    
    return result
