# README: Monte Carlo Analysis of Statistical Physical Models
This repository contains my work on the numerical analysis of the Ising model, focusing on advanced simulation techniques and Monte Carlo algorithms. The aim is to explore critical phenomena, phase transitions, and other physical properties of the Ising model, with a particular emphasis on the two-dimensional axial next-nearest neighbor Ising (ANNNI) model.

## Features

- **Monte Carlo Simulations**: Includes implementations of the Wolff algorithm and two-replica methods.
- **Population Annealing**: A robust approach to simulate systems at different temperatures efficiently.
- **Statistical Analysis**: Tools to compute observables such as energy, magnetization, and heat capacity.
- **Parallelization**: OpenMP-based parallelized algorithms for high performance on multi-core systems.
- **Checkpointing**: Save and resume simulations for long-running tasks.

## Repository Structure

```plaintext
.
├── src/                  # C++ source code
│   ├── main.cpp          # Main simulation file
│   ├── lattice.cpp       # Lattice model implementation
│   ├── population.cpp    # Population annealing methods
│   └── ...               # Other files
├── include/              # Header files
│   ├── parameters.h      # Simulation parameters
│   └── ...
├── data/                 # Generated simulation data (excluded from version control)
├── scripts/              # Python scripts for data analysis and plotting
├── tests/                # Unit tests for components
└── README.md             # Project overview
## Getting Started

```
This 



