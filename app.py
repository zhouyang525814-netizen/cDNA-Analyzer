# app.py
import os
import datetime
import customtkinter as ctk
from tkinter import filedialog, messagebox
from core_engine import DemultiplexEngine
import yaml
from analysis_engine import EnrichmentAnalyzer  # <--- 新增
import traceback
# 设定主题
ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")

def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans('ATCGN', 'TAGCN'))[::-1]

class cDNAAnalyzerApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("cDNA-DISPLAY Analyzer Pro - UI Framework")
        self.geometry("1450x950") 
        
        self.fastq_paths = []
        self.round_data = [] 
        
        self.grid_columnconfigure(0, weight=4) 
        self.grid_columnconfigure(1, weight=6)
        self.grid_rowconfigure(0, weight=1)
        
        self.setup_ui()

    def setup_ui(self):
        # ================= 左侧控制台 (滚动框架) =================
        left_frame = ctk.CTkScrollableFrame(self, fg_color="transparent")
        left_frame.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)
        
        # --- 1. Workspace & Data ---
        frame_data = ctk.CTkFrame(left_frame)
        frame_data.pack(pady=5, fill="x")
        ctk.CTkLabel(frame_data, text="1. Workspace & FASTQ Pool", font=("Arial", 16, "bold")).pack(anchor="w", padx=10, pady=5)
        
        ws_inner = ctk.CTkFrame(frame_data, fg_color="transparent")
        ws_inner.pack(fill="x", padx=10, pady=5)
        self.entry_pin = ctk.CTkEntry(ws_inner, placeholder_text="PIN", width=80)
        self.entry_pin.pack(side="left", padx=(0, 10))
        self.entry_project = ctk.CTkEntry(ws_inner, placeholder_text="Project Name", width=180)
        self.entry_project.pack(side="left", padx=(0, 10))
        
        btn_fastq = ctk.CTkButton(ws_inner, text="Browse FASTQs...", command=self.load_fastq_files, width=140)
        btn_fastq.pack(side="left")
        self.lbl_fastq = ctk.CTkLabel(ws_inner, text="0 files selected", text_color="gray")
        self.lbl_fastq.pack(side="left", padx=10)
        
        btn_clear_fastq = ctk.CTkButton(ws_inner, text="✖", width=30, fg_color="#8B0000", hover_color="#5C0000", command=self.clear_fastq_files)
        btn_clear_fastq.pack(side="left")

        # --- 2. Global Reference ---
        frame_ref = ctk.CTkFrame(left_frame)
        frame_ref.pack(pady=10, fill="x")
        ctk.CTkLabel(frame_ref, text="2. Global Reference Sequence", font=("Arial", 16, "bold")).pack(anchor="w", padx=10, pady=5)
        
        ref_inner = ctk.CTkFrame(frame_ref, fg_color="transparent")
        ref_inner.pack(fill="x", padx=10, pady=5)
        btn_fasta = ctk.CTkButton(ref_inner, text="Load FASTA...", command=self.load_fasta_file, width=120, fg_color="#4A4A4A")
        btn_fasta.pack(side="left", padx=(0, 10))
        
        self.entry_ref = ctk.CTkEntry(ref_inner, placeholder_text="Or paste Full Sequence (5' -> 3') here", width=350)
        self.entry_ref.pack(side="left", fill="x", expand=True)

        # --- 3. Rounds Configuration ---
        self.frame_rounds = ctk.CTkFrame(left_frame)
        self.frame_rounds.pack(pady=5, fill="both", expand=True)
        ctk.CTkLabel(self.frame_rounds, text="3. Round Definitions & Visual CDS", font=("Arial", 16, "bold")).pack(anchor="w", padx=10, pady=5)
        
        self.add_round_ui()
        self.add_round_ui()
        
        btn_add_round = ctk.CTkButton(self.frame_rounds, text="➕ Add Round", fg_color="#4A4A4A", hover_color="#333333", command=self.add_round_ui)
        btn_add_round.pack(pady=10)

        # --- 4. Biological Filters ---
        frame_filters = ctk.CTkFrame(left_frame)
        frame_filters.pack(pady=10, fill="x")
        ctk.CTkLabel(frame_filters, text="4. Quality & Biological Filters", font=("Arial", 16, "bold")).pack(anchor="w", padx=10, pady=5)
        
        self.var_adaptive = ctk.BooleanVar(value=True)
        ctk.CTkCheckBox(frame_filters, text="Adaptive: Allow length variation (In-frame indels only)", variable=self.var_adaptive).pack(anchor="w", padx=10, pady=5)
        
        self.var_stop = ctk.BooleanVar(value=True)
        ctk.CTkCheckBox(frame_filters, text="Discard sequences with premature stop codons (*)", variable=self.var_stop).pack(anchor="w", padx=10, pady=5)

        # --- Action Buttons ---
        frame_actions = ctk.CTkFrame(left_frame, fg_color="transparent")
        frame_actions.pack(pady=5, fill="x")
        
        self.btn_preview = ctk.CTkButton(frame_actions, text="🔍 Preview Alignment & Generate Rulers", fg_color="#1F538D", font=("Arial", 14, "bold"), height=35, command=self.run_preview)
        self.btn_preview.pack(fill="x", padx=10, pady=5)
        
        self.btn_run = ctk.CTkButton(frame_actions, text="▶ Run Pipeline", fg_color="#2A8C55", hover_color="#206A40", font=("Arial", 15, "bold"), height=40, command=self.execute_pipeline)
        self.btn_run.pack(fill="x", padx=10, pady=20) 

        # ================= 右侧控制台 (Terminal) =================
        right_frame = ctk.CTkFrame(self)
        right_frame.grid(row=0, column=1, sticky="nsew", padx=10, pady=10)
        ctk.CTkLabel(right_frame, text="Terminal & Visual Sequence Ruler", font=("Arial", 16, "bold")).pack(anchor="w", padx=10, pady=5)
        
        # Textbox 默认自带鼠标滚轮支持
        self.terminal = ctk.CTkTextbox(right_frame, wrap="none", font=("Courier", 14))
        self.terminal.pack(fill="both", expand=True, padx=10, pady=10)
        
        self.terminal.tag_config("header", foreground="#FFFFFF")
        self.terminal.tag_config("normal", foreground="#CCCCCC")
        self.terminal.tag_config("success", foreground="#00FFFF") 
        self.terminal.tag_config("error", foreground="#FF4444") 
        self.terminal.tag_config("warning", foreground="#FFA500") 
        self.terminal.tag_config("ruler_num", foreground="#FFFF00") 
        self.terminal.tag_config("ruler_tick", foreground="#888888") 
        self.terminal.tag_config("seq", foreground="#00FF00") 
        self.terminal.tag_config("info", foreground="#A9A9A9") # 暗灰色提示
        
        self.log([(">>> System Ready. Enter Data, Reference, and Primers, then click [Preview].\n", "normal")])
        self.log([(">>> (Scroll this terminal using your mouse wheel)\n\n", "info")])

    # --- 动态 UI 交互 ---
    def add_round_ui(self):
        idx = len(self.round_data)
        r_frame = ctk.CTkFrame(self.frame_rounds, fg_color="#333333")
        r_frame.pack(fill="x", pady=5, padx=10)
        
        ctk.CTkLabel(r_frame, text=f"Round {idx}", font=("Arial", 13, "bold")).pack(anchor="w", padx=5, pady=2)
        
        entry_fw = ctk.CTkEntry(r_frame, placeholder_text="Forward Primer (5'->3', incl. Barcode)", width=430)
        entry_fw.pack(padx=5, pady=2)
        entry_rv = ctk.CTkEntry(r_frame, placeholder_text="Reverse Primer (5'->3', anti-sense)", width=430)
        entry_rv.pack(padx=5, pady=2)
        
        cds_frame = ctk.CTkFrame(r_frame, fg_color="transparent")
        cds_frame.pack(fill="x", padx=5, pady=5)
        
        ctk.CTkLabel(cds_frame, text="CDS Start (-/0/+):").pack(side="left")
        entry_start = ctk.CTkEntry(cds_frame, width=55, state="disabled")
        entry_start.pack(side="left", padx=5)
        
        ctk.CTkLabel(cds_frame, text="CDS End:").pack(side="left", padx=(10,0))
        entry_end = ctk.CTkEntry(cds_frame, width=55, state="disabled")
        entry_end.pack(side="left", padx=5)
        
        status_lbl = ctk.CTkLabel(cds_frame, text="(Pending Preview)", text_color="gray", font=("Arial", 11))
        status_lbl.pack(side="left", padx=10)
        
        self.round_data.append({
            "id": idx, "fw_entry": entry_fw, "rv_entry": entry_rv,
            "start_entry": entry_start, "end_entry": entry_end, "status_lbl": status_lbl
        })

    def load_fastq_files(self):
        paths = filedialog.askopenfilenames(filetypes=[("FASTQ Files", "*.fastq *.fq")])
        if paths:
            self.fastq_paths = list(paths)
            self.lbl_fastq.configure(text=f"{len(paths)} files selected", text_color="white")

    def clear_fastq_files(self):
        self.fastq_paths = []
        self.lbl_fastq.configure(text="0 files selected", text_color="gray")

    def load_fasta_file(self):
        path = filedialog.askopenfilename(filetypes=[("FASTA Files", "*.fasta *.fa *.txt")])
        if path:
            seq = ""
            try:
                with open(path, 'r', encoding='utf-8') as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith('>'): continue
                        seq += line
                self.entry_ref.delete(0, "end")
                self.entry_ref.insert(0, seq.upper())
                self.log([(f"[*] Successfully loaded Reference Sequence from FASTA ({len(seq)} bp).\n", "success")])
            except Exception as e:
                messagebox.showerror("Error", f"Failed to read FASTA: {e}")

    def log(self, messages: list):
        for text, tag in messages:
            self.terminal.insert("end", text, tag)
        self.terminal.see("end") # 保证每次打印后都自动滚动到底部

    def _estimate_read_length(self) -> int:
        if not self.fastq_paths: return 150 
        try:
            with open(self.fastq_paths[0], 'r', encoding='utf-8') as f:
                f.readline()
                seq = f.readline().strip()
                if seq: return len(seq)
        except Exception:
            pass
        return 150

    # ================= 核心功能：分块打印序列与标尺 =================
    def _print_wrapped_sequence(self, display_seq: str, chunk_size: int = 40):
        """将序列分行打印，每行附带完整的刻度尺"""
        seq_len = len(display_seq)
        for i in range(0, seq_len, chunk_size):
            chunk = display_seq[i:i+chunk_size]
            nums = ""
            ticks = ""
            
            # 生成该 Chunk 对应的绝对刻度
            for j in range(len(chunk)):
                pos = i + j + 1
                if (pos - 1) % 10 == 0:
                    num_str = str(pos)
                    nums += num_str + " " * (10 - len(num_str))
                    ticks += "|" + "." * 9
            
            # 截断超出的刻度部分
            nums = nums[:len(chunk)]
            ticks = ticks[:len(chunk)]
            
            self.log([(nums + "\n", "ruler_num")])
            self.log([(ticks + "\n", "ruler_tick")])
            self.log([(chunk + "\n\n", "seq")])

    # --- 核心交互：精准生物学悬臂匹配与预览 ---
    def run_preview(self):
        self.terminal.delete("1.0", "end")
        ref_seq = self.entry_ref.get().strip().upper()
        if not ref_seq:
            self.log([("[ERROR] Please input or load the Global Reference Sequence first.\n", "error")])
            return

        avg_read_len = self._estimate_read_length()
        self.log([("=== Generating Visual Sequence Rulers ===\n", "header")])
        self.log([(f"[*] Physical Max Read Length: {avg_read_len} bp\n\n", "normal")])
        
        for r_dict in self.round_data:
            idx = r_dict["id"]
            fw_primer = r_dict["fw_entry"].get().strip().upper()
            rv_primer = r_dict["rv_entry"].get().strip().upper()
            
            if not fw_primer or not rv_primer:
                r_dict["status_lbl"].configure(text="(Skipped)", text_color="gray")
                continue
                
            # 1. 找左锚点：Fw引物的 3' 端最后 10bp
            fw_anchor = fw_primer[-10:] if len(fw_primer) >= 10 else fw_primer
            match_fw = ref_seq.find(fw_anchor)
            
            if match_fw == -1:
                self.log([(f"[!] Round {idx}: Fw Primer 3'-end ({fw_anchor}) not found in Reference.\n", "error")])
                r_dict["status_lbl"].configure(text="Fw Alignment Failed", text_color="#FF4444")
                continue
            start_pos = match_fw + len(fw_anchor)
            
            # 2. 找右锚点：Rv (5'->3') 的前 10bp (退火区)
            rc_rv_full = reverse_complement(rv_primer)
            rc_rv_anchor = rc_rv_full[:10] if len(rc_rv_full) >= 10 else rc_rv_full
            
            match_rv = ref_seq.find(rc_rv_anchor, start_pos)
            
            if match_rv == -1:
                self.log([(f"[!] Round {idx}: Rv Primer matching region ({rc_rv_anchor}) not found in Reference.\n", "error")])
                r_dict["status_lbl"].configure(text="Rv Alignment Failed", text_color="#FF4444")
                continue
                
            # 3. 计算物理读穿限制
            dist_to_rv = match_rv - start_pos
            read_capacity = avg_read_len - len(fw_primer)
            
            self.log([(f"--- Round {idx} Visual Ruler ---\n", "success")])
            self.log([(f"[INFO] Available Read Capacity: {avg_read_len} (Read) - {len(fw_primer)} (Fw Primer) = {read_capacity} bp\n", "info")])
            self.log([(f"[INFO] You can set CDS Start to negative values (up to -{len(fw_primer)}) to include Fw end.\n\n", "info")])
            
            if dist_to_rv <= read_capacity:
                self.log([("Status: ", "normal"), ("Full Read-through\n", "success")])
                display_seq = ref_seq[start_pos : match_rv]
                self._print_wrapped_sequence(display_seq, chunk_size=40)
                r_dict["status_lbl"].configure(text="Ready (Full Read-through)", text_color="#00FFFF")
                
            else:
                self.log([("Status: ", "normal"), ("Truncated Read (Sequence exceeds read capacity)\n", "warning")])
                display_seq = ref_seq[start_pos : start_pos + read_capacity]
                self._print_wrapped_sequence(display_seq, chunk_size=40)
                
                warn_msg = f"[!] WARNING: The actual distance to Rv Anchor is {dist_to_rv} bp, which exceeds your read capacity.\n"
                self.log([(warn_msg, "warning")])
                r_dict["status_lbl"].configure(text="Ready (Truncated Read)", text_color="#FFA500")

            r_dict["start_entry"].configure(state="normal")
            r_dict["end_entry"].configure(state="normal")
            
        self.log([("[System] Alignment complete. Please specify CDS coordinates on the left.\n", "header")])

    def execute_pipeline(self):
        self.terminal.delete("1.0", "end")
        self.log([(">>> Initializing Pipeline Configuration...\n", "header")])
        
        try:
            if not self.fastq_paths: raise ValueError("No FASTQ files selected.")
            if not self.entry_ref.get().strip(): raise ValueError("Reference sequence missing.")
            
            round_configs = {}
            for r_dict in self.round_data:
                fw = r_dict["fw_entry"].get().strip().upper()
                rv = r_dict["rv_entry"].get().strip().upper()
                if not fw or not rv: continue 
                
                start_val = r_dict["start_entry"].get().strip()
                end_val = r_dict["end_entry"].get().strip()
                if not start_val or not end_val:
                    raise ValueError(f"Round {r_dict['id']} is missing CDS Start/End coordinates.")
                
                # --- 支持负数与 0 的坐标校验 ---
                try:
                    start_idx, end_idx = int(start_val), int(end_val)
                except ValueError:
                    raise ValueError(f"Round {r_dict['id']} coordinates must be integers.")
                
                if end_idx < start_idx:
                    raise ValueError(f"Round {r_dict['id']}: CDS End must be greater than or equal to Start.")
                
                # 严格校验：负向扩展不能超过正向引物的总长度
                fw_len = len(fw)
                if start_idx < -fw_len:
                    raise ValueError(f"Round {r_dict['id']}: CDS Start ({start_idx}) cannot exceed Fw Primer length (-{fw_len}).")
                    
                # 计算总长度 (直接相减 + 1)
                length = end_idx - start_idx + 1
                if length % 3 != 0:
                    raise ValueError(f"Round {r_dict['id']}: CDS length ({length} bp) is NOT a multiple of 3! Frameshift detected.")
                
                round_configs[f"Round_{r_dict['id']}"] = {
                    "fw_primer": fw, "rv_primer": rv,
                    "cds_start": start_idx, "cds_end": end_idx
                }
                
            if not round_configs: raise ValueError("No rounds configured properly.")
                
            self.log([("[OK] Inputs & CDS lengths verified.\n", "success")])
            self.log([(f"[*] Adaptive Indel Mode: {self.var_adaptive.get()}\n", "normal")])
            self.log([(f"[*] Filter Premature Stops: {self.var_stop.get()}\n", "normal")])
            self.log([("\n>>> Handoff to Demultiplexing Plugin Interface...\n", "header")])
            
            global_settings = {
                "adaptive": self.var_adaptive.get(),
                "filter_stop": self.var_stop.get()
            }
            self._plugin_demultiplex_interface(self.fastq_paths, round_configs, self.entry_ref.get().strip(), global_settings)
            
        except ValueError as e:
            messagebox.showerror("Input Error", str(e))
            self.log([(f"[ERROR] {str(e)}\n", "error")])

    def _plugin_demultiplex_interface(self, fastq_list, round_configs_dict, reference_seq, settings):
        self.log([("[System] Booting Backend Core Engine...\n", "header")])
        
        # 1. 创建全局 Workspace 和项目文件夹
        pin = self.entry_pin.get().strip() or "0000"
        project_name = self.entry_project.get().strip().replace(" ", "_") or "Unnamed_Project"
        time_str = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        
        base_dir = os.path.join(os.path.expanduser("~"), "Documents", "cDNA_Analyzer_Workspace", pin, f"{project_name}_{time_str}")
        os.makedirs(base_dir, exist_ok=True)
        self.log([(f"[*] Workspace allocated at: {base_dir}\n", "info")])
        
        # 2. 备份用户的所有的输入参数到 config.yaml
        config_data = {
            "timestamp": time_str,
            "project_name": project_name,
            "pin": pin,
            "fastq_files": fastq_list,
            "global_reference": reference_seq,
            "settings": settings,
            "rounds": round_configs_dict
        }
        with open(os.path.join(base_dir, "run_config.yaml"), 'w', encoding='utf-8') as f:
            yaml.dump(config_data, f, default_flow_style=False)
            
        # 3. 启动引擎！
        def log_callback(text, tag):
            self.log([(text, tag)])
            
        engine = DemultiplexEngine(
            fastq_paths=fastq_list,
            round_configs=round_configs_dict,
            global_ref=reference_seq,
            settings=settings,
            log_callback=log_callback,
            output_dir=base_dir
        )
        
        # 接收跑完的数据
        dna_counters, stats, unassigned = engine.run()
        
# 4. 打印最终的漂亮统计报告
        self.log([("\n=== 🏁 Demultiplexing QC Report ===\n", "success")])
        for rnd, stat in stats.items():
            self.log([(f"[{rnd}]\n", "header")])
            self.log([(f"  - Total Demultiplexed:  {stat['total_assigned']:,}\n", "normal")])
            self.log([(f"  - ❌ Drop (Too Short / Truncated): {stat['discard_truncated']:,}\n", "error")])
            # (删除了 Frameshift 这行)
            self.log([(f"  - ❌ Drop (Stop Codon):   {stat['discard_stop_codon']:,}\n", "error")])
            self.log([(f"  - ✅ Passed Valid CDS:    {stat['passed_qc']:,}\n", "success")])
        
        self.log([(f"\n[!] Global Unassigned Reads (Orphan / Garbage): {unassigned:,}\n", "warning")])
        self.log([("\n[System] Core extraction finished! Data is safely stored in memory and log files.\n", "success")])
        
        # ======================== 新增接入 ========================
        # 保证轮次的排序是用户添加的顺序
        ordered_rounds = [f"Round_{r_dict['id']}" for r_dict in self.round_data if r_dict['fw_entry'].get().strip()]
        
        analyzer = EnrichmentAnalyzer(
            dna_counters=dna_counters,
            stats=stats,
            round_names=ordered_rounds,
            output_dir=base_dir,
            log_callback=log_callback
        )
        
        final_df = analyzer.run()
        
        if final_df is not None:
            self.log([("\n[System] Compiling Final Statistical QC Report...\n", "info")])
            
            try:
                import numpy as np
                report_path = os.path.join(base_dir, "QC_Summary_Report.txt")
                
                with open(report_path, 'w', encoding='utf-8') as f:
                    # ==========================================
                    # 0. 基础信息统计
                    # ==========================================
                    f.write("="*85 + "\n")
                    f.write("                cDNA-DISPLAY EXPERIMENT QC & SUMMARY REPORT\n")
                    f.write("="*85 + "\n\n")
                    f.write(f"Project Name    : {project_name}\n")
                    f.write(f"Generation Time : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                    f.write(f"Total Unique Valid Peptides Detected Across All Rounds: {len(final_df):,}\n\n")
                    
                    # ==========================================
                    # 1. 过滤与质量控制 (Demultiplexing & QC Traceability)
                    # ==========================================
                    f.write("--- 1. DEMULTIPLEXING & SEQUENCE FILTERING TRACEABILITY ---\n")
                    f.write(f"[*] Global Unassigned (Orphan/Low Quality) Reads: {unassigned:,}\n\n")
                    
                    header_qc = f"{'Round':<10} | {'Total Assigned':<15} | {'Truncated':<10} | {'Indel/Shift':<12} | {'Stop Codon':<11} | {'Passed QC':<12} | {'Yield (%)':<10}\n"
                    f.write(header_qc)
                    f.write("-" * len(header_qc) + "\n")
                    
                    for rnd in ordered_rounds:
                        s = stats.get(rnd, {})
                        tot = s.get('total_assigned', 0)
                        trunc = s.get('discard_truncated', 0)
                        indel = s.get('discard_length_indel', 0)
                        stop = s.get('discard_stop_codon', 0)
                        passed = s.get('passed_qc', 0)
                        yield_pct = (passed / tot * 100) if tot > 0 else 0.0
                        
                        f.write(f"{rnd:<10} | {tot:<15,} | {trunc:<10,} | {indel:<12,} | {stop:<11,} | {passed:<12,} | {yield_pct:>6.2f}%\n")

                    # ==========================================
                    # 2. 轮次间选择压力与多样性 (Selection Pressure & Complexity)
                    # ==========================================
                    f.write("\n--- 2. LIBRARY DIVERSITY & SELECTION PRESSURE ---\n")
                    header_div = f"{'Round':<10} | {'Passed QC':<12} | {'Unique Seqs':<12} | {'Complexity':<10} | {'Top 100 Takeover':<18} | {'Shannon Entropy'}\n"
                    f.write(header_div)
                    f.write("-" * len(header_div) + "\n")
                    
                    for rnd in ordered_rounds:
                        rpm_col = f"RPM_{rnd}"
                        data = final_df[rpm_col]
                        passed_qc = stats.get(rnd, {}).get('passed_qc', 0)
                        
                        unique_seqs = (data > 0).sum()
                        
                        # Complexity: Unique Seqs / Total Reads (反映测序饱和度与文库复杂性)
                        complexity = (unique_seqs / passed_qc * 100) if passed_qc > 0 else 0
                        
                        total_rpm = data.sum()
                        top100_takeover = (data.nlargest(100).sum() / total_rpm * 100) if total_rpm > 0 else 0
                        
                        p = (data[data > 0] / total_rpm).values
                        entropy = -np.sum(p * np.log2(p)) if total_rpm > 0 else 0
                        
                        f.write(f"{rnd:<10} | {passed_qc:<12,} | {unique_seqs:<12,} | {complexity:>6.2f}%   | {top100_takeover:>6.2f}%            | {entropy:.2f}\n")
                            
                    # ==========================================
                    # 3. Top 20 淘金榜 (移除了 PCR Bias 部分，直接顺延为第3部分)
                    # ==========================================
                    f.write("\n--- 3. TOP 20 ENRICHED CANDIDATES ---\n")
                    global_cols = [c for c in final_df.columns if c.startswith('Enrich_Global_')]
                    
                    if global_cols:
                        sort_col = global_cols[-1]
                        top_20 = final_df.nlargest(20, sort_col)
                        f.write(f"Sorted by: {sort_col}\n\n")
                        
                        header_top = f"{'Rank':<5} | {'Peptide_Seq':<25} | {'GC%':<5} | "
                        header_top += " | ".join([f"{rnd:>8}" for rnd in ordered_rounds])
                        header_top += f" | {'Final_FC':>8}\n"
                        f.write(header_top)
                        f.write("-" * len(header_top) + "\n")
                        
                        for i, (_, row) in enumerate(top_20.iterrows(), 1):
                            seq = str(row['Peptide_Seq'])[:25].ljust(25)
                            gc = f"{row['GC_Percent']:.1f}"
                            rpms = " | ".join([f"{row[f'RPM_{rnd}']:>8.1f}" for rnd in ordered_rounds])
                            fc = f"{row[sort_col]:>8.2f}"
                            
                            f.write(f"{i:<5} | {seq} | {gc:<5} | {rpms} | {fc}\n")
                    else:
                        # 如果没有计算富集倍数 (比如只输入了 1 轮数据)，则按最后一轮的 RPM 排序
                        sort_col = f"RPM_{ordered_rounds[-1]}"
                        top_20 = final_df.nlargest(20, sort_col)
                        f.write(f"Sorted by: {sort_col}\n\n")
                        
                        header_top = f"{'Rank':<5} | {'Peptide_Seq':<25} | {'GC%':<5} | "
                        header_top += " | ".join([f"{rnd:>8}" for rnd in ordered_rounds])
                        f.write(header_top + "\n")
                        f.write("-" * len(header_top) + "\n")
                        
                        for i, (_, row) in enumerate(top_20.iterrows(), 1):
                            seq = str(row['Peptide_Seq'])[:25].ljust(25)
                            gc = f"{row['GC_Percent']:.1f}"
                            rpms = " | ".join([f"{row[f'RPM_{rnd}']:>8.1f}" for rnd in ordered_rounds])
                            f.write(f"{i:<5} | {seq} | {gc:<5} | {rpms}\n")

                # ==========================================
                # 同步输出到 GUI 界面
                # ==========================================
                self.log([("\n📊 Traceable QC Report Generated Successfully!\n", "success")])
                self.log([(f"[*] A detailed text report has been saved to:\n    {report_path}\n", "normal")])
                self.log([("\n🏆 Check the Top Candidates directly in the text file.\n", "info")])
                self.log([("💡 Note: For visualization and bias diagnostics, import 'Master_Enrichment_Matrix.csv' into R, Python, or Prism.\n", "info")])
                
            except Exception as e:
                self.log([(f"\n[ERROR] Failed to generate QC summary: {str(e)}\n", "error")])

if __name__ == "__main__":
    app = cDNAAnalyzerApp()
    app.mainloop()
