This code is designed to calculate the one-point PDF of DM based on the baryonification framework(https://arxiv.org/pdf/2507.07892). The relevant formulas and theoretical background can be found in (https://arxiv.org/pdf/2601.18784).

The repository consists of two files:
1- pdf: the full version of the code optimized for running on the Marvin HPC system(https://wiki.hpc.uni-bonn.de/). It can, of course, be adapted to other computing environments.
2- Test_notebook: a lightweight notebook version of "pdf" that can be run on a personal computer. The numerical resolutions are intentionally reduced, so the results are not intended for precision analyses; however, it provides a useful starting point for becoming familiar with the functions and numerical methods used throughout the code.
To run this code, the BFC code is required as a prerequisite.
