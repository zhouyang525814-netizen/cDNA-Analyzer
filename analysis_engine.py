# analysis_engine.py
import os
import pandas as pd
import numpy as np

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

def translate_dna(seq):
    return "".join([CODON_TABLE.get(seq[i:i+3], 'X') for i in range(0, len(seq), 3)])

def calculate_gc(seq):
    if not seq: return 0.0
    return ((seq.count('G') + seq.count('C')) / len(seq)) * 100.0

class EnrichmentAnalyzer:
    def __init__(self, dna_counters, stats, round_names, output_dir, log_callback):
        self.dna_counters = dna_counters
        self.stats = stats
        self.round_names = round_names
        self.output_dir = output_dir
        self.log_callback = log_callback

    def _log(self, text, tag="normal"):
        self.log_callback(text, tag)
        with open(os.path.join(self.output_dir, "terminal_run.log"), 'a', encoding='utf-8') as f:
            f.write(text)

    def run(self):
        self._log("\n>>> [Analysis] Booting Enrichment Matrix Engine...\n", "header")
        
        # 1. 降维与聚类
        aa_records = {}
        for rnd in self.round_names:
            for dna, count in self.dna_counters.get(rnd, {}).items():
                aa = translate_dna(dna)
                if aa not in aa_records:
                    aa_records[aa] = {'counts': {r: 0 for r in self.round_names}, 'dna_freq': {}}
                aa_records[aa]['counts'][rnd] += count
                aa_records[aa]['dna_freq'][dna] = aa_records[aa]['dna_freq'].get(dna, 0) + count

        if not aa_records:
            return None

        # 2. 构建 Dataframe
        rows = []
        for aa, data in aa_records.items():
            dom_dna = max(data['dna_freq'].items(), key=lambda x: x[1])[0]
            row = {'Peptide_Seq': aa, 'Dominant_DNA_Seq': dom_dna, 'GC_Percent': calculate_gc(dom_dna)}
            for rnd in self.round_names:
                row[f'Count_{rnd}'] = data['counts'][rnd]
            rows.append(row)
        df = pd.DataFrame(rows)

        # 3. RPM 与 排名
        for rnd in self.round_names:
            total_valid = self.stats[rnd]['passed_qc']
            rpm_col = f'RPM_{rnd}'
            df[rpm_col] = (df[f'Count_{rnd}'] / total_valid) * 1e6 if total_valid > 0 else 0.0
            df[f'Rank_{rnd}'] = df[rpm_col].rank(ascending=False, method='min').astype(int)

        # 4. 计算所有必须的富集组合对
        pseudo = 1.0 
        
        # 4.1 阶梯富集 (Step-wise: 1v0, 2v1, 3v2)
        for i in range(1, len(self.round_names)):
            prev_r, curr_r = self.round_names[i-1], self.round_names[i]
            df[f'Enrich_Step_{curr_r}_vs_{prev_r}'] = np.log2((df[f'RPM_{curr_r}'] + pseudo) / (df[f'RPM_{prev_r}'] + pseudo))
            
        # 4.2 全局富集 (Global: 1v0, 2v0, 3v0)
        first_r = self.round_names[0]
        for i in range(1, len(self.round_names)):
            curr_r = self.round_names[i]
            df[f'Enrich_Global_{curr_r}_vs_{first_r}'] = np.log2((df[f'RPM_{curr_r}'] + pseudo) / (df[f'RPM_{first_r}'] + pseudo))

        df['Present_In_All'] = (df[[f'Count_{r}' for r in self.round_names]] > 0).all(axis=1)

        # 5. 保存
        try:
            # 默认按最后一轮的全局富集排序
            sort_col = f'Enrich_Global_{self.round_names[-1]}_vs_{first_r}' if len(self.round_names) > 1 else f'RPM_{self.round_names[0]}'
            df = df.sort_values(by=sort_col, ascending=False).reset_index(drop=True)
            
            parquet_path = os.path.join(self.output_dir, "Master_Enrichment_Matrix.parquet")
            csv_path = os.path.join(self.output_dir, "Master_Enrichment_Matrix.csv")
            
            df.to_parquet(parquet_path, index=False)
            df.to_csv(csv_path, index=False)
            self._log(f"[*] Exported {len(df):,} peptides to Master Matrix.\n", "success")
            return df
        except Exception as e:
            self._log(f"[ERROR] Failed to save result files: {str(e)}\n", "error")
            return None