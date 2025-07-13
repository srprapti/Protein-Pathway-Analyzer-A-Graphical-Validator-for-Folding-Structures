# Protein Energy Path Validator with GUI

## Overview
This project analyzes an amino acid sequence, simulates protein folding as an energy-weighted graph, detects misfolds (cycles), finds the minimum-energy folding path, and visualizes it using an SFML-based GUI.

## Features
- Input amino acid sequences (single-letter codes)  
- Graph construction with energy weights  
- Misfold simulation (random extra edges)  
- Cycle detection for misfolds  
- Shortest path (minimum-energy) calculation  
- GUI visualization with clickable amino acid nodes showing properties  

## Project Files
- `main.cpp`: Core logic and graph processing  
- `gui.cpp` / `gui.hpp`: SFML GUI rendering and interaction  
- `LiberationSans-Regular.ttf`: Font for GUI text  
- `README.md`: Project documentation  

## Requirements
- C++11+ compiler  
- SFML library installed (https://www.sfml-dev.org/tutorials/2.5/)  
- Font file included  

## Build & Run
```bash
g++ main.cpp gui.cpp -o protein_app -lsfml-graphics -lsfml-window -lsfml-system
./protein_app

Usage
Enter the amino acid sequence when prompted.

View misfold detection and minimum-energy path in the console.

GUI shows protein graph with edges and nodes colored by role.

Click nodes to see amino acid properties in console.
