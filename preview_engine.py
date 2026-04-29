# preview_engine.py

def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans('ATCGN', 'TAGCN'))[::-1]

def find_boundaries(seq: str, p1: str, rc_p2: str, min_match: int = 10):
    """Anchor search from 3' end to tolerate 5' truncation."""
    left_end = -1
    match_len_p1 = 0
    for i in range(len(p1), min_match - 1, -1):
        if seq.find(p1[-i:]) != -1:
            left_end = seq.find(p1[-i:]) + i
            match_len_p1 = i
            break
            
    if left_end == -1: return -1, -1, 0, 0
    
    right_start = -1
    match_len_p2 = 0
    for i in range(len(rc_p2), min_match - 1, -1):
        if seq.find(rc_p2[:i], left_end) != -1:
            right_start = seq.find(rc_p2[:i], left_end)
            match_len_p2 = i
            break
            
    if right_start == -1: return -1, -1, 0, 0
    return left_end, right_start, match_len_p1, match_len_p2

def generate_preview(fastq_path: str, p1: str, p2: str, num_reads: int = 5, min_match: int = 10) -> list:
    """Returns a list of (text_string, tag_name) for GUI syntax highlighting."""
    if not fastq_path or not p1 or not p2:
        return [("[Error] Please select a FASTQ file and input both primers.\n", "error")]
        
    p1 = p1.strip().upper()
    p2 = p2.strip().upper()
    rc_p2 = reverse_complement(p2)
    
    output = []
    output.append(("=== Target Alignment Preview (Top 5 Reads) ===\n", "header"))
    output.append((f"Rule: Tolerates 5' truncation. Minimum 3' continuous match: {min_match} bp.\n\n", "normal"))
    
    reads_processed = 0
    with open(fastq_path, 'r', encoding='utf-8') as f:
        while reads_processed < num_reads:
            header = f.readline()
            if not header: break
            seq = f.readline().rstrip()
            f.readline()
            f.readline() # qual
            
            read_id = header.strip().split()[0]
            
            left_end, right_start, m_p1, m_p2 = find_boundaries(seq, p1, rc_p2, min_match)
            target_seq = seq
            strand = "Sense"
            
            if left_end == -1:
                target_seq = reverse_complement(seq)
                left_end, right_start, m_p1, m_p2 = find_boundaries(target_seq, p1, rc_p2, min_match)
                strand = "Antisense (Auto-Flipped)"
                
            if left_end != -1 and right_start != -1:
                var_seq = target_seq[left_end:right_start]
                var_len = len(var_seq)
                
                pad_p1 = " " * (len(p1) - m_p1)
                disp_p1 = p1[-m_p1:]
                disp_p2_rev = p2[-m_p2:][::-1] 
                
                output.append((f"{read_id} [Original: {strand}]\n", "normal"))
                
                output.append((f"5' {pad_p1}", "normal"))
                output.append((disp_p1, "p1"))
                output.append((" 3' (Primer 1)", "normal"))
                if m_p1 < len(p1): output.append((f" [! 5' truncated {len(p1)-m_p1}bp]", "warning"))
                output.append(("\n", "normal"))
                
                output.append((f"5' {pad_p1}", "normal"))
                output.append((disp_p1, "p1"))
                output.append((" -- [", "normal"))
                
                if var_len % 3 != 0:
                    output.append((f"Var: {var_len} bp (Frameshift!)", "error"))
                else:
                    output.append((f"Var: {var_len} bp", "var"))
                    
                output.append(("] -- ", "normal"))
                output.append((rc_p2[:m_p2], "p2"))
                output.append((" 3'\n", "normal"))
                
                spacer = " " * len(f"5' {pad_p1}{disp_p1} -- [Var: {var_len} bp] -- ")
                if var_len % 3 != 0: spacer += " " * len(" (Frameshift!)")
                
                output.append((f"{spacer}3' ", "normal"))
                output.append((disp_p2_rev, "p2"))
                output.append((" 5' (Primer 2, right-to-left)", "normal"))
                if m_p2 < len(p2): output.append((f" [! 5' truncated {len(p2)-m_p2}bp]", "warning"))
                output.append(("\n" + "-" * 75 + "\n", "normal"))
                
            else:
                output.append((f"{read_id} ", "normal"))
                output.append(("[Warning] Dropped. Could not find sufficient primer anchors.\n", "error"))
                output.append(("-" * 75 + "\n", "normal"))
                
            reads_processed += 1
            
    return output