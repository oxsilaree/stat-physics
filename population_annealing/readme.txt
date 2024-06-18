The contents of this folder are used to make an executable that runs a population annealing simulation. 

At this time, I am working on having my code written and linked properly in several .h and .cpp files. The minimal working example right now can run a simulation of a single lattice, but there seem to be a couple of bugs that make it run imperfectly.

My goal is to also implement a population annealing class, and use OpenMP to allow it to run several lattices at once on different threads (hence the computing is parallelized). It should then also be able to take data at several points during the simulation and write files that I will use Python to analyze.
