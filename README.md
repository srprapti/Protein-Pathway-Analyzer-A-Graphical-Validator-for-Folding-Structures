 Draft 1 – Protein Structure Validator

 Overview:
- This version allows manual input of amino acid sequences.
- Builds a graph model of the protein.
- Simulates a misfold by inserting a random edge.
- Labels amino acids with type and charge.
- Detects cycles in the graph (misfold detection).
- Calculates the minimum-energy folding path using Dijkstra’s algorithm.

 Features Used:
- Graph data structure
- Random edge insertion
- Cycle detection
- Shortest path (Dijkstra)

 How to Compile:
Use any C++17 compiler. For Code::Blocks:
- Enable `-std=c++17` in build settings.
 Sample Input: MKTAYIA

