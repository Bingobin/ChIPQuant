def calculate_cpm(count, total_mapped_reads):
    """
    Calculate CPM (Counts Per Million).
    CPM = (Count / Total Mapped Reads) * 1,000,000
    """
    if total_mapped_reads <= 0:
        return 0.0
    return (count / total_mapped_reads) * 1_000_000

def get_promoter_region(chrom, start, end, strand, upstream, downstream):
    """
    Calculate promoter region based on gene position and strand.
    
    Args:
        start (int): Gene start (1-based)
        end (int): Gene end (1-based)
        strand (str): '+' or '-'
        upstream (int): Distance upstream of TSS
        downstream (int): Distance downstream of TSS
        
    Returns:
        (promoter_start, promoter_end) 1-based coordinates
    """
    if strand == "+":
        # TSS is start
        # Promoter: (start - upstream) to (start + downstream)
        tss = start
        p_start = tss - upstream
        p_end = tss + downstream
    elif strand == "-":
        # TSS is end
        # Promoter: (end - downstream) to (end + upstream) 
        # Note: 'upstream' for negative strand means larger genomic coordinates (right side).
        # But biologically upstream is 5' end.
        # For negative strand gene, 5' end is at 'end' position.
        # Negative strand upstream (5' side) is coordinate + 
        # Negative strand downstream (3' side) is coordinate -
        # So promoter region should be [end - downstream, end + upstream]
        
        # Example: TSS=1000, upstream=2000, downstream=500.
        # Range should be: from (1000 - 500) to (1000 + 2000) i.e. 500-3000.
        # Verification: 
        # 5' end is at 1000.
        # Upstream 2000bp means going towards 3000.
        # Downstream 500bp means going towards 500.
        # So interval is indeed [TSS - downstream, TSS + upstream].
        
        tss = end
        p_start = tss - downstream
        p_end = tss + upstream
    else:
        # Treat as positive strand by default or ignore
        tss = start
        p_start = tss - upstream
        p_end = tss + downstream
    
    # Ensure coordinates are at least 1
    return max(1, p_start), p_end
