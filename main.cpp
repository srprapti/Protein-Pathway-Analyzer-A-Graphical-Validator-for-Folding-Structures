import tkinter as tk
from tkinter import ttk, messagebox
import math
import random
import heapq
from collections import defaultdict

# -----------------------------
# Amino acid properties
# -----------------------------
class AminoAcid:
    def __init__(self, name, type_, charge):
        self.name = name
        self.type = type_
        self.charge = charge

AMINO_TABLE = {
    'A': AminoAcid("Ala", "Hydrophobic", "Neutral"),
    'C': AminoAcid("Cys", "Hydrophobic", "Neutral"),
    'D': AminoAcid("Asp", "Hydrophilic", "Negative"),
    'E': AminoAcid("Glu", "Hydrophilic", "Negative"),
    'F': AminoAcid("Phe", "Hydrophobic", "Neutral"),
    'G': AminoAcid("Gly", "Hydrophobic", "Neutral"),
    'H': AminoAcid("His", "Hydrophilic", "Positive"),
    'I': AminoAcid("Ile", "Hydrophobic", "Neutral"),
    'K': AminoAcid("Lys", "Hydrophilic", "Positive"),
    'L': AminoAcid("Leu", "Hydrophobic", "Neutral"),
    'M': AminoAcid("Met", "Hydrophobic", "Neutral"),
    'N': AminoAcid("Asn", "Hydrophilic", "Neutral"),
    'P': AminoAcid("Pro", "Hydrophobic", "Neutral"),
    'Q': AminoAcid("Gln", "Hydrophilic", "Neutral"),
    'R': AminoAcid("Arg", "Hydrophilic", "Positive"),
    'S': AminoAcid("Ser", "Hydrophilic", "Neutral"),
    'T': AminoAcid("Thr", "Hydrophilic", "Neutral"),
    'V': AminoAcid("Val", "Hydrophobic", "Neutral"),
    'W': AminoAcid("Trp", "Hydrophobic", "Neutral"),
    'Y': AminoAcid("Tyr", "Hydrophilic", "Neutral"),
}

# -----------------------------
# Core graph model
# -----------------------------
class ProteinGraph:
    def __init__(self):
        self.sequence = []
        self.N = 0
        self.graph = defaultdict(list)
        self.edge_set = set()
        self.positions = {}
        self.random_seed = None

    def set_sequence(self, s: str):
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

    def insert_random_misfold_edge(self, energy_min=1, energy_max=10):
        if self.N < 3:
            return None
        for _ in range(20):
            a, b = random.randrange(self.N), random.randrange(self.N)
            if a != b and (a, b) not in self.edge_set:
                w = random.randint(energy_min, energy_max)
                self.graph[a].append((b, w))
                self.edge_set.add((a, b))
                misfold_like = (b <= a) or (b - a > 1)
                return (a, b, w, misfold_like)
        return None

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
            angle = 2 * math.pi * i / self.N
            x = cx + radius * math.cos(angle)
            y = cy + radius * math.sin(angle)
            self.positions[i] = (x, y)

# -----------------------------
# GUI with visual enhancements
# -----------------------------
class ProteinGUI(tk.Tk):
    NODE_R = 25
    EDGE_COLOR = "#444444"
    EDGE_HIGHLIGHT = "#FF8C00"
    NODE_FILL_HYDROPHOBIC = "#6495ED"
    NODE_FILL_HYDROPHILIC = "#FF69B4"
    NODE_PATH = "#32CD32"
    NODE_OUTLINE = "#1E3A8A"
    TEXT_COLOR = "#ffffff"

    BG_COLOR = "#e0f7fa"
    PANEL_BG = "#f0f0f5"
    FONT_MAIN = ("Segoe UI", 12)
    FONT_TITLE = ("Segoe UI", 14, "bold")

    def __init__(self):
        super().__init__()
        self.title("Protein Pathway Analyzer - Tkinter GUI")
        self.geometry("1200x750")
        self.minsize(1050, 650)

        self.model = ProteinGraph()
        self.canvas_items_nodes = {}
        self.canvas_items_edges = []
        self.status_var = tk.StringVar()
        self._build_ui()
        self.status_var.set("Enter a sequence and click Build Graph.")

    def _build_ui(self):
        self.configure(bg=self.BG_COLOR)

        # Left Panel
        left = tk.Frame(self, bg=self.PANEL_BG, padx=12, pady=12)
        left.pack(side=tk.LEFT, fill=tk.Y)

        tk.Label(left, text="Protein Pathway Analyzer", font=self.FONT_TITLE, bg=self.PANEL_BG).pack(anchor="w", pady=(0,8))
        tk.Label(left, text="Amino Acid Sequence:", bg=self.PANEL_BG, font=self.FONT_MAIN).pack(anchor="w")
        self.entry_seq = tk.Entry(left, width=36, font=self.FONT_MAIN)
        self.entry_seq.pack(pady=4)
        self.entry_seq.insert(0, "ACDEFGHIKLMNPQRSTVWY")

        buttons = [
            ("Build Graph", self.on_build_graph),
            ("Insert Misfold Edge", self.on_misfold),
            ("Detect Cycles", self.on_detect_cycle),
            ("Animate Shortest Energy Path", self.on_shortest_path),
            ("Reset", self.on_reset)
        ]
        for text, cmd in buttons:
            b = tk.Button(left, text=text, command=cmd, font=self.FONT_MAIN, bg="#00acc1", fg="white", relief="raised", bd=2)
            b.pack(fill=tk.X, pady=4)

        # Table
        tk.Label(left, text="Amino Acid Table:", bg=self.PANEL_BG, font=self.FONT_MAIN).pack(anchor="w", pady=(12,0))
        self.tree = ttk.Treeview(left, columns=("pos","aa","name","type","charge"), show="headings", height=12)
        for col, w in [("pos", 44), ("aa", 44), ("name", 80), ("type", 120), ("charge", 80)]:
            self.tree.heading(col, text=col.upper())
            self.tree.column(col, width=w, anchor="center")
        self.tree.pack(fill=tk.BOTH, pady=8)

        # Status
        self.lbl_status = tk.Label(left, textvariable=self.status_var, bg=self.PANEL_BG, anchor="w", font=self.FONT_MAIN)
        self.lbl_status.pack(fill=tk.X, pady=(4,0))

        # Right Canvas
        right = tk.Frame(self, bg=self.BG_COLOR, padx=8, pady=12)
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.canvas = tk.Canvas(right, bg=self.BG_COLOR, highlightthickness=0)
        self.canvas.pack(fill=tk.BOTH, expand=True)
        self.canvas.bind("<Configure>", self._redraw_on_resize)
        self.canvas.bind("<Button-1>", self.on_canvas_click)

    # ----------- Actions -----------
    def on_build_graph(self):
        seq = self.entry_seq.get().strip()
        self.model.set_sequence(seq)
        if self.model.N==0:
            messagebox.showwarning("Invalid Input", "No valid amino acids found.")
            return
        self.model.build_backbone()
        self.populate_table()
        self.layout_and_draw()
        self.status_var.set(f"Graph built with {self.model.N} residues.")

    def on_misfold(self):
        res = self.model.insert_random_misfold_edge()
        if not res:
            self.status_var.set("Could not insert misfold edge.")
        else:
            a,b,w,misfold_like = res
            msg = f"Inserted edge {a+1}->{b+1} (energy {w})."
            if misfold_like:
                msg += " Looks like misfold edge."
            self.status_var.set(msg)
        self.layout_and_draw()

    def on_detect_cycle(self):
        if self.model.has_cycle():
            messagebox.showinfo("Cycle Detection", "Misfold detected: cycle present.")
            self.status_var.set("Cycle detected in graph.")
        else:
            messagebox.showinfo("Cycle Detection", "No cycles detected.")
            self.status_var.set("No cycles detected.")

    def on_shortest_path(self):
        path, cost = self.model.dijkstra()
        if not path:
            self.status_var.set("No path to last residue.")
            return
        self.status_var.set(f"Animating minimum-energy path (cost={cost})...")
        self._animate_path(path, cost)

    def _animate_path(self, path, cost):
        self.layout_and_draw()
        for i, idx in enumerate(path):
            circle_id, text_id = self.canvas_items_nodes[idx]
            self.canvas.itemconfig(circle_id, fill=self.NODE_PATH)
            self.canvas.update()
            self.after(100)  # small delay for animation effect
        for i in range(1, len(path)):
            u, v = path[i - 1], path[i]
            for line, _, uu, vv, _ in self.canvas_items_edges:
                if uu == u and vv == v:
                    self.canvas.itemconfig(line, fill=self.EDGE_HIGHLIGHT, width=4)
                    break
        self.status_var.set(f"Shortest-energy path highlighted (cost={cost}).")

    def on_reset(self):
        self.model = ProteinGraph()
        self.tree.delete(*self.tree.get_children())
        self.canvas.delete("all")
        self.canvas_items_nodes.clear()
        self.canvas_items_edges.clear()
        self.status_var.set("Reset complete. Enter sequence and click Build Graph.")

    def populate_table(self):
        self.tree.delete(*self.tree.get_children())
        for i, aa in enumerate(self.model.sequence, start=1):
            info = AMINO_TABLE.get(aa)
            self.tree.insert("", "end", values=(i, aa, info.name, info.type, info.charge))

    # ------------- Canvas Drawing -------------
    def layout_and_draw(self):
        self.canvas.delete("all")
        self.canvas_items_nodes.clear()
        self.canvas_items_edges.clear()

        w,h = self.canvas.winfo_width() or 900, self.canvas.winfo_height() or 600
        cx,cy = w//2, h//2
        radius = min(w,h)*0.35
        radius = max(180,min(radius,260))
        self.model.layout_circle(cx,cy,radius)

        for u in range(self.model.N):
            for v,w_energy in self.model.graph[u]:
                self._draw_edge(u,v,w_energy)

        for i in range(self.model.N):
            self._draw_node(i)

    def _draw_node(self, idx):
        x, y = self.model.positions[idx]
        r = self.NODE_R
        aa = self.model.sequence[idx]
        fill_color = self.NODE_FILL_HYDROPHOBIC if AMINO_TABLE[aa].type == "Hydrophobic" else self.NODE_FILL_HYDROPHILIC
        circle = self.canvas.create_oval(x - r, y - r, x + r, y + r, fill=fill_color, outline=self.NODE_OUTLINE, width=2)
        # text slightly above center for full visibility
        text = self.canvas.create_text(x, y-2, text=f"{aa}{idx+1}", fill=self.TEXT_COLOR, font=("Segoe UI", 10, "bold"))
        self.canvas_items_nodes[idx] = (circle, text)

    def _draw_edge(self,u,v,w_energy):
        x1,y1 = self.model.positions[u]
        x2,y2 = self.model.positions[v]
        line = self.canvas.create_line(x1,y1,x2,y2, fill=self.EDGE_COLOR, width=2)
        mx,my = (x1+x2)/2,(y1+y2)/2
        self.canvas.create_text(mx+12,my+12,text=str(w_energy), fill="#222222", font=("Segoe UI",9))
        self.canvas_items_edges.append((line,None,u,v,w_energy))

    def _redraw_on_resize(self,event):
        if self.model.N==0: return
        self.layout_and_draw()

    # Node Click
    def on_canvas_click(self,event):
        x,y = event.x,event.y
        for idx,(circle_id,_) in self.canvas_items_nodes.items():
            x1,y1,x2,y2 = self.canvas.coords(circle_id)
            cx,cy = (x1+x2)/2,(y1+y2)/2
            if (x-cx)**2 + (y-cy)**2 <= self.NODE_R**2:
                aa = self.model.sequence[idx]
                info = AMINO_TABLE.get(aa)
                self.status_var.set(f"Residue {idx+1}: {aa} ({info.name}) | {info.type}, {info.charge}")
                return

if __name__=="__main__":
    random.seed()
    app = ProteinGUI()
    app.mainloop()
