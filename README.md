# readme: Monte Carlo Analysis of Statistical Physical Models
This repository contains my work on the numerical analysis of the Ising model, where I've focused on combining parallelization and Monte Carlo algorithms. The aim is to explore critical phenomena, phase transitions, and other physical properties of the Ising model, with a particular emphasis on the two-dimensional axial next-nearest neighbor Ising (ANNNI) model.

For now, the `population_annealing` folder is up to date, and this contains our latest simulation files that combine Population Annealing and Two-Replica cluster moves for the ANNNI model.

Before going into detail about the repository itself, here is a quick overview of the model we have been investigating:

## The ANNNI Model

The ANNNI model is just like the Ising model, except with added next-nearest neighbor interactions along one axis (regardless of total number of dimensions). This axis is typically denoted as the $z$-axis. The Hamiltonian of this model is

$$ \mathcal{H} = -J \sum_{\langle nn\rangle} \sigma_i \sigma_j - J_1 \sum_{\langle nnn \rangle} \sigma_i \sigma_j \hspace{2pt} , $$

where $J, J_1$ are coupling constants for the different interactions, $\sigma$ denotes individual spins, and the indices of the sums $\langle nn \rangle$ and $\langle nnn \rangle$ indicate that the sums are being taken over nearest neighbors and axial next-nearest neighbors respectively. Typically, $J_1$ takes an opposite sign to $J$, which introduces competing interactions in the system leading to frustration.


## Features

- **Monte Carlo Simulations**: Includes implementations of the Wolff algorithm and two-replica methods.
- **Population Annealing**: A robust approach to simulate systems at different temperatures efficiently, and allows for parallelization of code.
- **Statistical Analysis**: Tools to compute observables such as heat capacity, cluster wrapping (in any direction), fourier analysis and more.
- **Parallelization**: OpenMP-based parallelized algorithms for high performance on multi-core systems.

## Repository Structure

```plaintext
.
├── population_annealing/        # C++ source code
│   ├── makefile                # Compile and link files
│   ├── build.sh                # Runs make and make clean 
│   ├── main.cpp                # Main simulation file
│   ├── functions.cpp           # Includes neighbor table creation
│   ├── spin_class.cpp          # Implementation of spin object
│   ├── lattice_class.cpp       # Implementation of spin-lattice object, with Wolff cluster moves
│   ├── population_class.cpp    # Implementation of lattice population, with population annealing methods and two-replica cluster moves
│   └── include/                # Header files
│       ├── parameters.h        # Simulation global variables
│       ├── functions.h      
│       ├── spin_class.h      
│       ├── lattice_class.h      
│       └── population_class.h      
│   ├── data/                       # Any generated data would go here
│   ├── production-run/             # Large data sets would go here
│   ├── data_analysis/              # Python scripts for data analysis and visualization
│   └── README.md          
```
## Getting Started
### Prerequisites

Before you begin, ensure you have the following tools installed:

- **C++ Compiler**: A compiler with support for C++17 or later (e.g., `g++`, `clang++`).
  - **GNU Scientific Library (GSL)**: [Download and install GSL](https://www.gnu.org/software/gsl/).
  - **FFTW**: C++ Fast Fourier Transform library [Download and install FFTW](https://www.fftw.org).
  - **OpenMP**: Required for parallelized simulations, usually baked into the compiler but can be installed if needed.
- **Python 3.x**: Along with the following Python libraries:
  - `numpy`
  - `matplotlib`
  - `pandas`
  - `glob`
  - `scipy`

### Installation and Usage

Follow these steps to set up and run the project:

1. **Clone the repository**:
   Download the repository to your local machine using `git`:
   ```bash
   git clone https://github.com/oxsilaree/stat-physics.git
   ```
2. **Edit global variables**:
   In `parameters.h`, one might want to change some parameters
   - Side length of lattices: $L$
   - Maximum inverse temperature to anneal to: `MAX_BETA`
   - Minimum distance between paired replicas: `MIN_DISTANCE`
   - Most other parameters are deprecated at this point.
     
   Additionally, at this step you should change the directories in the `collectData` and `loadData` functions for data collection if necessary.

4. **Build program**
   ```bash
   ./build.sh
   ```
5. **Run program**:
   To run this program, supply the 3 following parameters:
   - `kappa`: This is the coupling constant ratio $\kappa = -J_1/J$. Typical ranges are float values from 0 to 2.
   - `method`: "t" for Two Replica, and "p" for Wolff ("s" for simulated annealing with Wolff is depreciated at the moment)
   - `R`: Population size. For the code to run without incident, this should be some size larger than `MIN_DISTANCE` in `parameters.h`.
   ```bash
   ./main <kappa> <method> <R>
   ```



