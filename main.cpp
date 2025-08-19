import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import math
import random
import heapq
import csv
from collections import defaultdict, Counter

# -----------------------------
# Amino acid properties
# -----------------------------
class AminoAcid:
    def __init__(self, name, type_, charge):
        self.name = name
        self.type = type_
        self.charge = charge

# single-letter to descriptive mapping
AMINO_TABLE = {
    'A': AminoAcid("Alanine","Hydrophobic", "Neutral"),
    'C': AminoAcid("Cysteine","Hydrophobic", "Neutral"),
    'D': AminoAcid("Aspartate","Hydrophilic", "Negative"),
    'E': AminoAcid("Glutamate","Hydrophilic", "Negative"),
    'F': AminoAcid("Phenylalanine","Hydrophobic", "Neutral"),
    'G': AminoAcid("Glycine","Hydrophobic", "Neutral"),
    'H': AminoAcid("Histidine","Hydrophilic", "Positive"),
    'I': AminoAcid("Isoleucine","Hydrophobic", "Neutral"),
    'K': AminoAcid("Lysine","Hydrophilic", "Positive"),
    'L': AminoAcid("Leucine","Hydrophobic", "Neutral"),
    'M': AminoAcid("Methionine","Hydrophobic", "Neutral"),
    'N': AminoAcid("Asparagine","Hydrophilic", "Neutral"),
    'P': AminoAcid("Proline","Hydrophobic", "Neutral"),
    'Q': AminoAcid("Glutamine","Hydrophilic", "Neutral"),
    'R': AminoAcid("Arginine","Hydrophilic", "Positive"),
    'S': AminoAcid("Serine","Hydrophilic", "Neutral"),
    'T': AminoAcid("Threonine","Hydrophilic", "Neutral"),
    'V': AminoAcid("Valine","Hydrophobic", "Neutral"),
    'W': AminoAcid("Tryptophan","Hydrophobic", "Neutral"),
    'Y': AminoAcid("Tyrosine","Hydrophilic", "Neutral"),
}

# -----------------------------
# Core graph model
# -----------------------------
class ProteinGraph:
    def __init__(self):
        self.sequence = []
        self.N = 0
        self.graph = defaultdict(list)   # adjacency list: {u: [(v, weight), ...]}
        self.edge_set = set()            # set of (u,v) pairs for quick membership
        self.positions = {}              # layout positions {index: (x,y)}

    def set_sequence(self, s: str):
        # keep only valid single-letter amino acids
        seq = [c.upper() for c in s.strip() if c.upper() in AMINO_TABLE]
        self.sequence = seq
        self.N = len(seq)
        self.graph.clear()
        self.edge_set.clear()

    def build_backbone(self, energy_min=1, energy_max=10, seed=None):
        if seed is not None:
            random.seed(seed)
        self.graph.clear()
        self.edge_set.clear()
        for i in range(self.N - 1):
            w = random.randint(energy_min, energy_max)
            self.graph[i].append((i + 1, w))
            self.edge_set.add((i, i + 1))

    def insert_edge(self, a, b, w):
        if (a, b) in self.edge_set or a == b or not (0 <= a < self.N and 0 <= b < self.N):
            return False
        self.graph[a].append((b, w))
        self.edge_set.add((a, b))
        return True

    def remove_edge(self, a, b):
        removed_once = False
        new_adj = []
        for v, w in self.graph[a]:
            if not removed_once and v == b:
                removed_once = True
                continue
            new_adj.append((v, w))
        self.graph[a] = new_adj
        if removed_once:
            if not any(v == b for v, _ in self.graph[a]):
                self.edge_set.discard((a, b))
        return removed_once

    def set_edge_weight(self, a, b, new_w):
        for i, (v, w) in enumerate(self.graph[a]):
            if v == b:
                self.graph[a][i] = (v, new_w)
                return True
        return False

    def ensure_edge(self, a, b, w):
        if not (0 <= a < self.N and 0 <= b < self.N) or a == b:
            return False
        if (a, b) in self.edge_set:
            # keep minimum weight
            for i, (v, ww) in enumerate(self.graph[a]):
                if v == b and w < ww:
                    self.graph[a][i] = (v, w)
            return True
        self.graph[a].append((b, w))
        self.edge_set.add((a, b))
        return True

    def insert_random_misfold_edge(self, energy_min=1, energy_max=10):
        if self.N < 3:
            return None
        for _ in range(120):
            a, b = random.randrange(self.N), random.randrange(self.N)
            if a != b and (a, b) not in self.edge_set:
                w = random.randint(energy_min, energy_max)
                self.graph[a].append((b, w))
                self.edge_set.add((a, b))
                misfold_like = (b <= a) or (b - a > 1)
                return (a, b, w, misfold_like)
        return None

    # ---- Analysis helpers ----
    def indegree_counts(self):
        indeg = Counter()
        for u in range(self.N):
            for v, _ in self.graph[u]:
                indeg[v] += 1
        return indeg

    def has_cycle(self):
        visited = [False] * self.N
        recstack = [False] * self.N

        def dfs(u):
            visited[u] = True
            recstack[u] = True
            for v, _ in self.graph[u]:
                if not visited[v]:
                    if dfs(v):
                        return True
                elif recstack[v]:
                    return True
            recstack[u] = False
            return False

        for i in range(self.N):
            if not visited[i] and dfs(i):
                return True
        return False

    def dijkstra(self, source=0, target=None):
        if self.N == 0:
            return [], float("inf")
        if target is None:
            target = self.N - 1
        dist = [float("inf")] * self.N
        prev = [-1] * self.N
        dist[source] = 0
        pq = [(0, source)]
        while pq:
            d, u = heapq.heappop(pq)
            if d > dist[u]:
                continue
            if u == target:
                break
            for v, w in self.graph[u]:
                nd = d + w
                if nd < dist[v]:
                    dist[v] = nd
                    prev[v] = u
                    heapq.heappush(pq, (nd, v))
        if dist[target] == float("inf"):
            return [], float("inf")
        path = []
        at = target
        while at != -1:
            path.append(at)
            at = prev[at]
        path.reverse()
        return path, dist[target]

    def layout_circle(self, cx, cy, radius):
        self.positions.clear()
        if self.N == 0:
            return
        for i in range(self.N):
            angle = 2 * math.pi * i / (self.N if self.N else 1)
            x = cx + radius * math.cos(angle)
            y = cy + radius * math.sin(angle)
            self.positions[i] = (x, y)

    # ---- Error classification ----
    def classify_errors(self):
        issues = []
        seen = set()
        for u in range(self.N):
            for v, w in self.graph[u]:
                key = (u, v)
                if key in seen:
                    issues.append(dict(type="DuplicateEdge", detail=f"{u+1}->{v+1}", u=u, v=v, weight=w))
                seen.add(key)

                if u == v:
                    issues.append(dict(type="SelfLoop", detail=f"{u+1}->{v+1}", u=u, v=v, weight=w))
                elif v == u + 1:
                    pass
                elif v <= u:
                    issues.append(dict(type="BackwardEdge", detail=f"{u+1}->{v+1}", u=u, v=v, weight=w))
                elif v - u > 1:
                    issues.append(dict(type="LongSkipEdge", detail=f"{u+1}->{v+1}", u=u, v=v, weight=w))

                if w >= 9:
                    issues.append(dict(type="HighEnergyEdge", detail=f"{u+1}->{v+1} (E={w})", u=u, v=v, weight=w))

        indeg = self.indegree_counts()
        for i in range(self.N):
            outdeg = len(self.graph[i])
            if i != self.N - 1 and outdeg == 0:
                issues.append(dict(type="DeadEnd", detail=f"node {i+1}", u=i, v=None, weight=None))
            if indeg[i] > 1 and i != 0:
                issues.append(dict(type="MultipleIncoming", detail=f"node {i+1} indeg={indeg[i]}", u=None, v=i, weight=None))

        if self.has_cycle():
            issues.append(dict(type="CyclePresent", detail="Graph contains at least one cycle", u=None, v=None, weight=None))
        return issues

# -----------------------------
# GUI (3-column layout, with file import features)
# -----------------------------
class ProteinGUI(tk.Tk):
    # Colors & fonts
    BG = "#eef7fb"
    PANEL = "#ffffff"
    ACCENT = "#00bcd4"
    BTN = "#02a9b5"
    BTN2 = "#6c63ff"
    WARN = "#ff7043"
    OK = "#43a047"
    TEXT = "#1d2b36"
    NODE_HPHOBIC = "#4aa3ff"
    NODE_HPHILIC = "#ff72b6"
    NODE_PATH = "#31c36a"
    EDGE_PATH = "#ff9800"
    NODE_OUTLINE = "#1E3A8A"

    NODE_R = 22
    FONT = ("Segoe UI", 11)
    FONT_B = ("Segoe UI", 12, "bold")

    def __init__(self):
        super().__init__()
        self.title("Protein Pathway Analyzer")
        self.configure(bg=self.BG)
        self.geometry("1380x820")
        self.minsize(1200, 720)

        # model/state
        self.model = ProteinGraph()
        self.canvas_nodes = {}
        self.canvas_edges = []
        self.current_path_nodes = set()
        self.current_path_edges = set()
        self.selected_idx = None
        self.highlight_hydrophobic = False
        self.highlight_hydrophilic = False

        # status
        self.status_var = tk.StringVar(value="Type a sequence, load a file, then Build Graph.")

        # layout grid columns
        self.grid_columnconfigure(0, weight=0)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=0)
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=0)

        # build UI
        self._build_left_panel()
        self._build_center_canvas()
        self._build_right_panel()
        self._build_statusbar()

    # ---------------- UI parts ----------------
    def _build_left_panel(self):
        left = tk.Frame(self, bg=self.PANEL)
        left.grid(row=0, column=0, sticky="nsw", padx=(12,8), pady=12)

        tk.Label(left, text="Controls", font=self.FONT_B, bg=self.PANEL, fg=self.TEXT).pack(anchor="w", pady=(0,8))

        # manual entry + merge combobox
        row = tk.Frame(left, bg=self.PANEL)
        row.pack(fill="x", pady=(2,6))
        tk.Label(row, text="Sequence", font=self.FONT, bg=self.PANEL).pack(side="left")
        self.entry_seq = tk.Entry(row, width=28, font=self.FONT)
        self.entry_seq.pack(side="left", padx=(8,0))
        self.entry_seq.insert(0, "ACDEFGHIKLMNPQRSTVWY")

        # merge mode and file buttons
        fm = tk.Frame(left, bg=self.PANEL)
        fm.pack(fill="x", pady=(6,6))
        tk.Label(fm, text=" Merge:", bg=self.PANEL).pack(side="left")
        self.merge_mode = tk.StringVar(value="Use Input")
        cmb = ttk.Combobox(fm, textvariable=self.merge_mode, values=["Use Input", "Use File", "Append"], width=12, state="readonly")
        cmb.pack(side="left", padx=(6,4))

        tk.Button(fm, text="Load PDB", command=self.on_load_pdb, bg=self.BTN, fg="white").pack(side="left", padx=(6,0))
        tk.Button(fm, text="Load PDF", command=self.on_load_pdf, bg=self.BTN2, fg="white").pack(side="left", padx=(6,0))

        # build graph
        tk.Button(left, text="Build Graph", command=self.on_build_graph, bg=self.BTN, fg="white", font=self.FONT).pack(fill="x", pady=(10,8))

        # tools area
        box = tk.LabelFrame(left, text="Graph Tools", bg=self.PANEL, labelanchor="n")
        box.pack(fill="x", pady=(4,8))
        tk.Button(box, text="Insert Random Misfold Edge", command=self.on_misfold, bg=self.BTN2, fg="white").pack(fill="x", pady=4)
        tk.Button(box, text="Detect Misfolds (Direct)", command=self.on_detect_misfolds_direct, bg=self.BTN2, fg="white").pack(fill="x", pady=4)
        tk.Button(box, text="Suggest Corrections", command=self.on_suggest_corrections, bg=self.BTN2, fg="white").pack(fill="x", pady=4)

        # paths (ONLY shortest path; removed Start/End Find & Highlight)
        pbox = tk.LabelFrame(left, text="Paths", bg=self.PANEL, labelanchor="n")
        pbox.pack(fill="x", pady=(4,8))
        tk.Button(pbox, text="Shortest Energy Path", command=self.on_shortest_path, bg=self.OK, fg="white").pack(fill="x", pady=4)

        # highlights
        hbox = tk.LabelFrame(left, text="Highlights", bg=self.PANEL, labelanchor="n")
        hbox.pack(fill="x", pady=(4,8))
        tk.Button(hbox, text="Toggle Hydrophobic", command=self.on_toggle_hydrophobic, bg="#2ecc71", fg="white").pack(fill="x", pady=3)
        tk.Button(hbox, text="Toggle Hydrophilic", command=self.on_toggle_hydrophilic, bg="#e67e22", fg="white").pack(fill="x", pady=3)

        # export + reset
        xbox = tk.LabelFrame(left, text="Export", bg=self.PANEL, labelanchor="n")
        xbox.pack(fill="x", pady=(4,8))
        tk.Button(xbox, text="Export PNG", command=self.on_export_png, bg="#4caf50", fg="white").pack(fill="x", pady=3)
        tk.Button(xbox, text="Export CSV", command=self.on_export_csv, bg="#4caf50", fg="white").pack(fill="x", pady=3)
        tk.Button(left, text="Reset", command=self.on_reset, bg=self.WARN, fg="white").pack(fill="x", pady=(6,0))

    def _build_center_canvas(self):
        center = tk.Frame(self, bg=self.BG)
        center.grid(row=0, column=1, sticky="nsew", pady=12)
        center.grid_columnconfigure(0, weight=1)
        center.grid_rowconfigure(0, weight=1)

        self.canvas = tk.Canvas(center, bg=self.BG, highlightthickness=0)
        self.canvas.grid(row=0, column=0, sticky="nsew", padx=4)
        self.canvas.bind("<Configure>", self._redraw_on_resize)
        self.canvas.bind("<Button-1>", self.on_canvas_click)

    def _build_right_panel(self):
        right = tk.Frame(self, bg=self.PANEL)
        right.grid(row=0, column=2, sticky="nse", padx=(8,12), pady=12)

        tk.Label(right, text="Amino Acid Table", font=self.FONT_B, bg=self.PANEL, fg=self.TEXT).pack(anchor="w", pady=(0,8))

        cols = ("pos", "aa", "name", "type", "charge")
        self.tree = ttk.Treeview(right, columns=cols, show="headings", height=1)
        for col, w in [("pos", 52), ("aa", 50), ("name", 130), ("type", 150), ("charge", 90)]:
            self.tree.heading(col, text=col.upper())
            self.tree.column(col, width=w, anchor="center")
        self.tree.pack(fill="x")
        self.tree.bind("<<TreeviewSelect>>", self.on_table_select)

        # detail card
        self.card = tk.Frame(right, bg=self.PANEL)
        self.card.pack(fill="x", pady=(10,0))
        self.card_title = tk.Label(self.card, text="Details", font=self.FONT_B, bg=self.PANEL, fg=self.TEXT)
        self.card_title.pack(anchor="w")
        self.card_body = tk.Label(self.card, text="Select a row to see info.", font=self.FONT, bg=self.PANEL,
                                  justify="left", wraplength=380)
        self.card_body.pack(anchor="w", pady=(6,0))

    def _build_statusbar(self):
        bar = tk.Frame(self, bg=self.PANEL)
        bar.grid(row=1, column=0, columnspan=3, sticky="ew", padx=12, pady=(0,10))
        tk.Label(bar, textvariable=self.status_var, bg=self.PANEL, anchor="w", font=self.FONT, fg=self.TEXT).pack(fill="x")

    # ---------------- File import helpers ----------------
    def on_load_pdb(self):
        """
        Load sequence from PDB using Biopython (PPBuilder).
        Merge according to self.merge_mode.
        """
        path = filedialog.askopenfilename(title="Select PDB file", filetypes=[("PDB files", "*.pdb *.ent"), ("All files", "*.*")])
        if not path:
            return

        try:
            from Bio.PDB import PDBParser, PPBuilder
        except Exception:
            messagebox.showerror("Missing dependency", "Biopython is required for PDB parsing.\nInstall with:\n\npip install biopython")
            return

        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("prot", path)
            ppb = PPBuilder()
            seqs = []
            for pp in ppb.build_peptides(structure):
                seqs.append(str(pp.get_sequence()))
            if not seqs:
                messagebox.showwarning("No sequence", "Could not extract sequence from PDB.")
                return
            pdb_seq = "".join(seqs)
            # keep only valid single-letter amino acids
            pdb_seq = "".join([c for c in pdb_seq if c.upper() in AMINO_TABLE])
            self._merge_sequence_from_file(pdb_seq, source_label=f"PDB ({path.split('/')[-1]})")
        except Exception as e:
            messagebox.showerror("PDB Load Error", f"Failed to parse PDB:\n{e}")

    def on_load_pdf(self):
        """
        Load sequence-like text from a PDF using pypdf.
        Heuristic: find strings of valid single-letter amino acids (length >= 3),
        join them as sequence text. Merge according to self.merge_mode.
        """
        path = filedialog.askopenfilename(title="Select PDF file", filetypes=[("PDF files", "*.pdf"), ("All files", "*.*")])
        if not path:
            return

        try:
            from pypdf import PdfReader
        except Exception:
            # older package name "PyPDF2" may also be installed; try that
            try:
                import PyPDF2 as _p
                PdfReader = _p.PdfReader
            except Exception:
                messagebox.showerror("Missing dependency", "pypdf (or PyPDF2) is required for PDF parsing.\nInstall with:\n\npip install pypdf")
                return

        try:
            reader = PdfReader(path)
            text_parts = []
            for page in reader.pages:
                txt = page.extract_text() or ""
                text_parts.append(txt)
            raw = "\n".join(text_parts)
            # Heuristic: extract only uppercase letters that correspond to amino acids
            letters = [c for c in raw if c.upper() in AMINO_TABLE]
            seq_candidate = "".join(letters)
            if not seq_candidate:
                messagebox.showwarning("No sequence", "Could not find amino-acid single-letter sequence in PDF.")
                return
            # Optionally compress whitespace (we already filtered letters)
            self._merge_sequence_from_file(seq_candidate, source_label=f"PDF ({path.split('/')[-1]})")
        except Exception as e:
            messagebox.showerror("PDF Load Error", f"Failed to parse PDF:\n{e}")

    def _merge_sequence_from_file(self, file_seq: str, source_label: str = "File"):
        """Merge file-provided sequence into the manual entry according to merge_mode."""
        mode = self.merge_mode.get()
        user_seq = self.entry_seq.get().strip().upper()
        user_seq = "".join([c for c in user_seq if c in AMINO_TABLE])
        if mode == "Use File":
            new_seq = file_seq
        elif mode == "Append":
            # append file seq after user seq (use whichever present)
            new_seq = (user_seq or "") + (file_seq or "")
        else:  # Use Input
            new_seq = user_seq or file_seq
        # set entry and notify
        self.entry_seq.delete(0, tk.END)
        self.entry_seq.insert(0, new_seq)
        self.status_var.set(f"{source_label} loaded ({len(file_seq)} residues). Merge mode: {mode}.")

    # ----------------- Actions -----------------
    def on_build_graph(self):
        seq = self.entry_seq.get().strip().upper()
        self.model.set_sequence(seq)
        if self.model.N == 0:
            messagebox.showwarning("Invalid Input", "No valid amino acids found in sequence.")
            return
        self.model.build_backbone()
        self.current_path_nodes.clear()
        self.current_path_edges.clear()
        self.selected_idx = None

        self.populate_table()
        self.layout_and_draw()
        self._resize_table_and_window()
        self.status_var.set(f"Graph built with {self.model.N} residues.")

    def on_misfold(self):
        res = self.model.insert_random_misfold_edge()
        if not res:
            self.status_var.set("Could not insert misfold edge (sequence too short or no free slot).")
        else:
            a, b, w, like = res
            msg = f"Inserted edge {a+1}->{b+1} (E={w})."
            if like:
                msg += " Looks misfold-like."
            self.status_var.set(msg)
        self.layout_and_draw()

    def on_detect_misfolds_direct(self):
        if self.model.N == 0:
            self.status_var.set("Build a graph first.")
            return
        issues = self.model.classify_errors()
        if not issues:
            messagebox.showinfo("Misfold Detection", "No misfolds or anomalies detected.")
            return
        self._show_issues_table(issues)

    def on_suggest_corrections(self):
        if self.model.N == 0:
            self.status_var.set("Build a graph first.")
            return
        issues = self.model.classify_errors()
        if not issues:
            messagebox.showinfo("Suggestions", "No anomalies detected to correct.")
            return
        suggestions = []
        for it in issues:
            t = it["type"]
            u, v, w = it.get("u"), it.get("v"), it.get("weight")
            if t in ("BackwardEdge", "LongSkipEdge"):
                if u is not None and v is not None:
                    suggestions.append(dict(action="remove_edge", detail=f"Remove {u+1}->{v+1} ({t})", u=u, v=v, param=None))
                    if u < self.model.N - 1:
                        suggestions.append(dict(action="ensure_edge", detail=f"Add backbone {u+1}->{u+2}", u=u, v=u+1, param=3))
            elif t == "SelfLoop" and u is not None and v is not None:
                suggestions.append(dict(action="remove_edge", detail=f"Remove self-loop {u+1}->{v+1}", u=u, v=v, param=None))
            elif t == "DuplicateEdge" and u is not None and v is not None:
                suggestions.append(dict(action="remove_edge", detail=f"Remove duplicate {u+1}->{v+1}", u=u, v=v, param=None))
            elif t == "HighEnergyEdge" and u is not None and v is not None:
                suggestions.append(dict(action="set_edge_weight", detail=f"Reduce energy {u+1}->{v+1} → 4", u=u, v=v, param=4))
            elif t == "DeadEnd" and u is not None and u < self.model.N - 1:
                suggestions.append(dict(action="ensure_edge", detail=f"Fix dead-end: {u+1}->{u+2}", u=u, v=u+1, param=3))
            elif t == "MultipleIncoming" and it["v"] is not None:
                target = it["v"]
                incoming = []
                for uu in range(self.model.N):
                    for vv, ww in self.model.graph[uu]:
                        if vv == target:
                            incoming.append((uu, vv, ww))
                preferred = None
                for uu, vv, ww in incoming:
                    if uu == target - 1:
                        preferred = (uu, vv, ww); break
                if preferred is None and incoming:
                    preferred = min(incoming, key=lambda x: x[2])
                for uu, vv, ww in incoming:
                    if preferred and (uu, vv) == (preferred[0], preferred[1]):
                        continue
                    suggestions.append(dict(action="remove_edge", detail=f"Trim extra incoming: remove {uu+1}->{vv+1}", u=uu, v=vv, param=None))
            elif t == "CyclePresent":
                suggestions.append(dict(action="advice", detail="Cycle present: remove a backward/loop edge to break cycle", u=None, v=None, param=None))
        if not suggestions:
            messagebox.showinfo("Suggestions", "No actionable corrections generated.")
            return
        self._show_corrections_table(suggestions)

    def on_shortest_path(self):
        path, cost = self.model.dijkstra()
        if not path:
            self.status_var.set("No path to last residue.")
            return
        self._animate_path(path, cost)

    def on_toggle_hydrophobic(self):
        self.highlight_hydrophobic = not self.highlight_hydrophobic
        self.layout_and_draw()

    def on_toggle_hydrophilic(self):
        self.highlight_hydrophilic = not self.highlight_hydrophilic
        self.layout_and_draw()

    def on_export_png(self):
        try:
            from PIL import ImageGrab
        except Exception:
            messagebox.showerror("Missing dependency", "Pillow is required for PNG export.\nInstall with:\n\npip install pillow")
            return
        self.canvas.update()
        x = self.canvas.winfo_rootx()
        y = self.canvas.winfo_rooty()
        w = x + self.canvas.winfo_width()
        h = y + self.canvas.winfo_height()
        img = ImageGrab.grab(bbox=(x, y, w, h))
        path = filedialog.asksaveasfilename(title="Save graph as PNG", defaultextension=".png", filetypes=[("PNG Image", "*.png")])
        if path:
            img.save(path)
            self.status_var.set(f"Graph exported to {path}")

    def on_export_csv(self):
        if self.model.N == 0:
            self.status_var.set("Nothing to export. Build a graph first.")
            return
        path = filedialog.asksaveasfilename(title="Save graph data as CSV", defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if not path:
            return
        try:
            with open(path, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["# Nodes"])
                writer.writerow(["Index(1-based)", "AA", "Name", "Type", "Charge"])
                for i, aa in enumerate(self.model.sequence, start=1):
                    info = AMINO_TABLE[aa]
                    writer.writerow([i, aa, info.name, info.type, info.charge])
                writer.writerow([])
                writer.writerow(["# Edges"])
                writer.writerow(["From(1-based)", "To(1-based)", "Energy"])
                for u in range(self.model.N):
                    for v, w in self.model.graph[u]:
                        writer.writerow([u+1, v+1, w])
            self.status_var.set(f"Graph exported to {path}")
        except Exception as e:
            messagebox.showerror("Export Error", f"Failed to export CSV:\n{e}")

    def on_reset(self):
        self.model = ProteinGraph()
        self.tree.delete(*self.tree.get_children())
        self.canvas.delete("all")
        self.canvas_nodes.clear()
        self.canvas_edges.clear()
        self.current_path_edges.clear()
        self.current_path_nodes.clear()
        self.selected_idx = None
        self.tree.config(height=1)
        self.card_body.config(text="Select a row to see info.")
        self._resize_table_and_window(force_min=True)
        self.status_var.set("Reset complete. Enter a sequence and click Build Graph.")

    # ----------------- Table & details -----------------
    def populate_table(self):
        self.tree.delete(*self.tree.get_children())
        for i, aa in enumerate(self.model.sequence, start=1):
            info = AMINO_TABLE[aa]
            self.tree.insert("", "end", values=(i, aa, info.name, info.type, info.charge))
        self.tree.config(height=max(1, self.model.N))

    def on_table_select(self, _evt=None):
        sel = self.tree.selection()
        if not sel:
            return
        vals = self.tree.item(sel[0], "values")
        pos = int(vals[0]) - 1
        aa = vals[1]
        info = AMINO_TABLE[aa]
        self.card_title.config(text=f"Residue {pos+1}: {aa}")
        self.card_body.config(text=f"Name: {info.name}\nType: {info.type}\nCharge: {info.charge}")
        self.status_var.set(f"Selected {aa}{pos+1} – {info.type}, {info.charge}")
        self._select_canvas_node(pos)

    # ----------------- Canvas drawing -----------------
    def layout_and_draw(self):
        self.canvas.delete("all")
        self.canvas_nodes.clear()
        self.canvas_edges.clear()

        w = self.canvas.winfo_width() or 900
        h = self.canvas.winfo_height() or 600
        cx, cy = w // 2, h // 2
        radius = min(w, h) * 0.36
        radius = max(180, min(radius, 290))
        self.model.layout_circle(cx, cy, radius)

        # draw edges first so nodes appear on top
        for u in range(self.model.N):
            for v, w_energy in self.model.graph[u]:
                self._draw_edge(u, v, w_energy)

        for i in range(self.model.N):
            self._draw_node(i)

        # re-apply highlights and path styling
        for idx in self.current_path_nodes:
            if idx in self.canvas_nodes:
                circle, _ = self.canvas_nodes[idx]
                self.canvas.itemconfig(circle, fill=self.NODE_PATH)
        for (u, v) in self.current_path_edges:
            for line, _, uu, vv, _ in self.canvas_edges:
                if uu == u and vv == v:
                    self.canvas.itemconfig(line, fill=self.EDGE_PATH, width=4)

        if self.selected_idx is not None:
            self._focus_ring(self.selected_idx)

    def _edge_color_by_energy(self, w_energy):
        t = max(0.0, min(1.0, (w_energy - 1) / 9.0))
        if t < 0.5:
            k = t / 0.5
            r = int(0 + k * (255 - 0))
            g = int(170 + k * (200 - 170))
            b = 0
        else:
            k = (t - 0.5) / 0.5
            r = int(255 + k * (220 - 255))
            g = int(200 + k * (40 - 200))
            b = int(0 + k * (40 - 0))
        return f"#{r:02x}{g:02x}{b:02x}"

    def _draw_node(self, idx):
        x, y = self.model.positions[idx]
        r = self.NODE_R
        aa = self.model.sequence[idx]
        base_fill = self.NODE_HPHOBIC if AMINO_TABLE[aa].type == "Hydrophobic" else self.NODE_HPHILIC

        outline = self.NODE_OUTLINE
        width = 2
        if self.highlight_hydrophobic and AMINO_TABLE[aa].type == "Hydrophobic":
            outline = "#145a32"; width = 3
        if self.highlight_hydrophilic and AMINO_TABLE[aa].type == "Hydrophilic":
            outline = "#7e3517"; width = 3

        circle = self.canvas.create_oval(x-r, y-r, x+r, y+r, fill=base_fill, outline=outline, width=width)

        label = f"{aa}{idx+1}"
        text_ids = []
        halo_color = "#1a1a1a"
        for dx, dy in ((-1,0),(1,0),(0,-1),(0,1)):
            tid = self.canvas.create_text(x + dx, y - 2 + dy, text=label, fill=halo_color, font=("Segoe UI", 10, "bold"))
            text_ids.append(tid)
        main_id = self.canvas.create_text(x, y - 2, text=label, fill="white", font=("Segoe UI", 10, "bold"))
        text_ids.append(main_id)
        self.canvas_nodes[idx] = (circle, text_ids)

    def _draw_edge(self, u, v, w_energy):
        x1, y1 = self.model.positions[u]
        x2, y2 = self.model.positions[v]
        color = self._edge_color_by_energy(w_energy)
        line = self.canvas.create_line(x1, y1, x2, y2, fill=color, width=2, arrow="last")
        mx, my = (x1 + x2)/2, (y1 + y2)/2
        label_id = self.canvas.create_text(mx+12, my+12, text=str(w_energy), fill="#222", font=("Segoe UI", 9))
        self.canvas_edges.append((line, label_id, u, v, w_energy))

    def _redraw_on_resize(self, _evt):
        if self.model.N > 0:
            self.layout_and_draw()

    # click canvas node => select table row and show details
    def on_canvas_click(self, event):
        x, y = event.x, event.y
        for idx, (circle_id, _) in self.canvas_nodes.items():
            x1, y1, x2, y2 = self.canvas.coords(circle_id)
            cx, cy = (x1 + x2) / 2, (y1 + y2) / 2
            if (x - cx)**2 + (y - cy)**2 <= self.NODE_R**2:
                self._select_table_row(idx)
                info = AMINO_TABLE[self.model.sequence[idx]]
                self.status_var.set(f"Residue {idx+1}: {self.model.sequence[idx]} ({info.name}) | {info.type}, {info.charge}")
                return

    def _select_table_row(self, idx):
        # clear selection, then select the corresponding Treeview entry
        for i in self.tree.selection():
            self.tree.selection_remove(i)
        children = self.tree.get_children()
        if 0 <= idx < len(children):
            self.tree.selection_set(children[idx])
            self.tree.see(children[idx])
            self.on_table_select()

    def _select_canvas_node(self, idx):
        self.selected_idx = idx
        self._focus_ring(idx)

    def _focus_ring(self, idx):
        self.canvas.delete("focusring")
        if idx not in self.canvas_nodes:
            return
        circle_id, _ = self.canvas_nodes[idx]
        x1, y1, x2, y2 = self.canvas.coords(circle_id)
        pad = 6
        self.canvas.create_oval(x1-pad, y1-pad, x2+pad, y2+pad, outline="#222", width=2, dash=(3,2), tags="focusring")

    # --------------- Animation ---------------
    def _animate_path(self, path, cost):
        self.current_path_nodes = set(path)
        self.current_path_edges = set((path[i-1], path[i]) for i in range(1, len(path)))
        self.layout_and_draw()
        for idx in path:
            if idx in self.canvas_nodes:
                circle, _ = self.canvas_nodes[idx]
                self.canvas.itemconfig(circle, fill=self.NODE_PATH)
                self.canvas.update()
                self.after(90)
        self.status_var.set(f"Shortest-energy path highlighted (cost={cost}).")

    # --------------- Popups ---------------
    def _show_issues_table(self, issues):
        win = tk.Toplevel(self)
        win.title("Detected Misfolds / Anomalies")
        cols = ("type", "detail", "u", "v", "weight")
        tree = ttk.Treeview(win, columns=cols, show="headings", height=20)
        for c, w in zip(cols, (140, 260, 70, 70, 80)):
            tree.heading(c, text=c.upper()); tree.column(c, width=w, anchor="center")
        tree.pack(fill="both", expand=True, padx=8, pady=8)
        for it in issues:
            u = "" if it["u"] is None else it["u"] + 1
            v = "" if it["v"] is None else it["v"] + 1
            w = "" if it["weight"] is None else it["weight"]
            tree.insert("", "end", values=(it["type"], it["detail"], u, v, w))
        tk.Button(win, text="Close", command=win.destroy, bg="#777", fg="white").pack(pady=6)

    def _show_corrections_table(self, suggestions):
        win = tk.Toplevel(self)
        win.title("Suggested Corrections")
        cols = ("action", "detail", "u", "v", "param")
        tree = ttk.Treeview(win, columns=cols, show="headings", height=22, selectmode="extended")
        widths = (140, 360, 70, 70, 80)
        for c, w in zip(cols, widths):
            tree.heading(c, text=c.upper()); tree.column(c, width=w, anchor="center")
        tree.pack(fill="both", expand=True, padx=8, pady=8)
        for s in suggestions:
            u = "" if s["u"] is None else s["u"] + 1
            v = "" if s["v"] is None else s["v"] + 1
            param = "" if s["param"] is None else s["param"]
            tree.insert("", "end", values=(s["action"], s["detail"], u, v, param))

        def apply_selected():
            sel = tree.selection()
            if not sel:
                messagebox.showwarning("No selection", "Select one or more suggestions.")
                return
            applied = 0
            for iid in sel:
                action, detail, u1, v1, param = tree.item(iid, "values")
                u = None if u1 == "" else int(u1) - 1
                v = None if v1 == "" else int(v1) - 1
                if param == "":
                    p = None
                else:
                    try:
                        p = int(param)
                    except Exception:
                        p = param
                if action == "remove_edge" and u is not None and v is not None:
                    if self.model.remove_edge(u, v):
                        applied += 1
                elif action == "ensure_edge" and u is not None and v is not None:
                    if self.model.ensure_edge(u, v, p if isinstance(p, int) else 3):
                        applied += 1
                elif action == "set_edge_weight" and u is not None and v is not None:
                    if self.model.set_edge_weight(u, v, p if isinstance(p, int) else 4):
                        applied += 1
            if applied > 0:
                self.status_var.set(f"Applied {applied} correction(s).")
                self.layout_and_draw()
            else:
                self.status_var.set("No corrections applied.")

        btns = tk.Frame(win, bg=self.PANEL)
        btns.pack(pady=(0,8))
        tk.Button(btns, text="Apply Selected", command=apply_selected, bg=self.ACCENT, fg="white").pack(side="left", padx=4)
        tk.Button(btns, text="Close", command=win.destroy, bg="#777", fg="white").pack(side="left", padx=4)

    # --------------- Window / sizing ---------------
    def _resize_table_and_window(self, force_min=False):
        rows = max(1, self.model.N)
        self.tree.config(height=rows)
        row_px = 24
        right_base = 180
        need_h = right_base + rows * row_px
        cur_w = self.winfo_width() or 1380
        cur_h = self.winfo_height() or 820
        new_h = max(cur_h, need_h, 720 if force_min else 0)
        self.geometry(f"{cur_w}x{int(new_h)}")


# -----------------------------
# Run
# -----------------------------
if __name__ == "__main__":
    random.seed()
    app = ProteinGUI()
    app.mainloop()
