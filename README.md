Protein Pathway Analyzer

A Tkinter-based interactive tool for visualizing and analyzing protein structures as graphs.

This project allows users to input an amino acid sequence, visualize it as a circular protein graph, detect misfolds/cycles, insert random misfold edges, and animate the minimum-energy pathway between residues. It also provides detailed amino acid information on click.

Features

Custom Sequence Input: Enter any valid amino acid sequence and generate the corresponding protein graph.

Graph Visualization: Displays residues as circular nodes with edges representing the protein backbone or misfold connections.

Animated Shortest-Energy Path: Highlights the minimum-energy path dynamically and keeps the path permanently highlighted after animation.

Cycle Detection: Detects potential misfolds in the protein graph.

Interactive Node Inspection: Click on any node to view amino acid properties (name, type, charge).

Misfold Simulation: Insert random edges to simulate potential misfolding or unusual connections.

Responsive UI: Graph layout adapts to window resizing.

Color-Coded Nodes and Edges: Easily distinguish normal residues, misfold edges, and highlighted paths.

Installation

Clone the repository:

git clone https://github.com/your-username/protein-pathway-analyzer.git
cd protein-pathway-analyzer


Install Python dependencies:

This project only uses Python’s standard library (tkinter, random, math, heapq, collections), so no additional packages are required.

Make sure you have Python 3.x installed with Tkinter support.

Run the application:

python protein_pathway_analyzer.py

Usage

Launch the application.

Enter your amino acid sequence (e.g., ACDEFGHIKLMNPQRSTVWY) in the entry field.

Click Build Graph to generate the protein structure.

Use the buttons to:

Insert Misfold Edge

Detect Cycles

Animate Shortest Energy Path

Reset the graph

Click any node to view detailed amino acid properties in the status bar.
