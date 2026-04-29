# core_engine.py
import os
import yaml
import datetime
from collections import Counter

# 标准遗传密码表
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans('ATCGN', 'TAGCN'))[::-1]

class DemultiplexEngine:
    def __init__(self, fastq_paths, round_configs, global_ref, settings, log_callback, output_dir):
        self.fastq_paths = fastq_paths
        self.round_configs = round_configs
        self.settings = settings
        self.log_callback = log_callback
        self.output_dir = output_dir
        
        self.stats = {rnd: {
            "total_assigned": 0,
            "discard_truncated": 0,
            "discard_length_indel": 0,
            "discard_stop_codon": 0,
            "passed_qc": 0
        } for rnd in round_configs.keys()}
        
        self.global_unassigned = 0
        self.dna_counters = {rnd: Counter() for rnd in round_configs.keys()}
        
        self._preprocess_primers()

    def _log(self, text, tag="normal"):
        self.log_callback(text, tag)
        with open(os.path.join(self.output_dir, "terminal_run.log"), 'a', encoding='utf-8') as f:
            f.write(text)

    def _preprocess_primers(self):
        """统一使用 10bp 作为绝对退火锚点"""
        for rnd, cfg in self.round_configs.items():
            fw = cfg['fw_primer']
            rv = cfg['rv_primer']
            
            # Fw: 最后 10bp 作为找坐标的锚点，前面的作为 Barcode 悬臂
            anchor_len = min(10, len(fw))
            cfg['fw_anchor'] = fw[-anchor_len:]
            cfg['fw_barcode'] = fw[:-anchor_len] if len(fw) > anchor_len else ""
            
            # Rv: 反向互补后的前 10bp
            rc_rv = reverse_complement(rv)
            cfg['rv_anchor'] = rc_rv[:10] if len(rc_rv) >= 10 else rc_rv

    def _map_coordinate(self, pos):
        return pos - 1

    def _process_read(self, seq: str):
        """核心处理：带有严格 Delta 校验与单端盲切分拣系统"""
        
        max_barcode_error = 1.0 
        min_victory_margin = 1.0 
        
        score_board = [] 
        
        # --- 1. 容错打分寻找 Fw 锚点 ---
        for rnd, cfg in self.round_configs.items():
            idx = seq.find(cfg['fw_anchor'])
            if idx != -1:
                expected_bc = cfg['fw_barcode']
                read_barcode = seq[max(0, idx - len(expected_bc)) : idx]
                
                score = 0.0
                len_diff = len(expected_bc) - len(read_barcode)
                
                if len_diff > 0:
                    score += len_diff * 1.0 
                    compare_exp = expected_bc[len_diff:]
                else:
                    compare_exp = expected_bc
                    
                for e, r in zip(compare_exp, read_barcode):
                    if r == 'N': score += 0.5   
                    elif e != r: score += 1.0   
                        
                current_fw_end = idx + len(cfg['fw_anchor'])
                score_board.append((score, rnd, current_fw_end))
                
        if not score_board:
            return False 
            
        score_board.sort(key=lambda x: x[0])
        best_score, best_rnd, fw_end_idx = score_board[0]
        
        if best_score > max_barcode_error:
            return False 
            
        if len(score_board) > 1:
            runner_up_score = score_board[1][0]
            if runner_up_score - best_score < min_victory_margin:
                return False 
                
        assigned_rnd = best_rnd
        self.stats[assigned_rnd]["total_assigned"] += 1
        rnd_cfg = self.round_configs[assigned_rnd]
        
        # --- 2. 绝对坐标切片 (彻底修复单端读不穿的问题) ---
        start_offset = self._map_coordinate(rnd_cfg['cds_start'])
        # Python 切片的右边界是不包含的，所以 End 坐标就是偏移量本身
        cds_start_abs = fw_end_idx + start_offset
        cds_end_abs = fw_end_idx + rnd_cfg['cds_end'] 
        
        # 截断检查：如果 Read 的物理长度连用户指定的 CDS End 都覆盖不到，丢弃
        if cds_end_abs > len(seq) or cds_start_abs < 0:
            self.stats[assigned_rnd]["discard_truncated"] += 1
            return True
            
        cds_seq = seq[cds_start_abs : cds_end_abs]
        cds_len = len(cds_seq)
        
        # --- 3. (可选) 大缺失校验 ---
        # 如果 Read 极其碰巧长到了能读出 Rv 锚点，我们才做大缺失校验
        rv_idx = seq.find(rnd_cfg['rv_anchor'], fw_end_idx)
        if rv_idx != -1 and cds_end_abs > rv_idx and not self.settings['adaptive']:
            self.stats[assigned_rnd]["discard_length_indel"] += 1
            return True
        
        # --- 4. 移码校验 ---
        if cds_len % 3 != 0:
            self.stats[assigned_rnd]["discard_length_indel"] += 1
            return True
            
        # --- 5. 翻译与生物学过滤 ---
        aa_seq = ""
        for i in range(0, cds_len, 3):
            codon = cds_seq[i:i+3]
            aa_seq += CODON_TABLE.get(codon, 'X')
            
        if self.settings['filter_stop'] and '*' in aa_seq:
            self.stats[assigned_rnd]["discard_stop_codon"] += 1
            return True
            
        # --- 6. 金库入账 ---
        self.stats[assigned_rnd]["passed_qc"] += 1
        self.dna_counters[assigned_rnd][cds_seq] += 1
        
        return True
    
    def run(self):
        self._log("\n>>> Demultiplexing & Extraction Engine Started...\n", "header")
        
        for file_path in self.fastq_paths:
            self._log(f"[*] Processing file: {os.path.basename(file_path)}...\n", "info")
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    while True:
                        header = f.readline()
                        if not header: break
                        seq = f.readline().rstrip()
                        f.readline()

                        qual = f.readline().rstrip() # 拿到质量分数行
                        
                        # --- 新增：Q-Score 平均质量过滤 (Threshold = 20) ---
                        if sum(ord(c) - 33 for c in qual) / len(qual) < 20.0:
                            self.global_unassigned += 1 # 质量太差，直接按废弃处理
                            continue
                        # 先查正链
                        success = self._process_read(seq)
                        # 翻转查反向互补链
                        if not success:
                            rc_seq = reverse_complement(seq)
                            success = self._process_read(rc_seq)
                            
                        if not success:
                            self.global_unassigned += 1
            except Exception as e:
                self._log(f"[ERROR] Failed reading {file_path}: {e}\n", "error")
                
        self._log(">>> Execution Complete. Compiling Summary...\n", "header")
        return self.dna_counters, self.stats, self.global_unassigned