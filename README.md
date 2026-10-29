# structure-presevingCBO


## 1. Description

This repository contains the MATLAB codes used to generate the numerical results presented in the paper:

**Title:** Computation of entanglement for quantum states by a Consensus-Based Optimization method

**Authors:** Michael Herty, Yijia Tang, and Yizhou Zhou


## 2. Code Description

- 'src/'
  
The src files include several methods used to optimize the entanglement of formation for a given quantum state.

The two structure-preserving CBO methods proposed in the paper: Hermitian preserving and Stiefel preserving CBO, see CBO_Hermitian.m and CBO_Unitary.m

The two projection-based CBO methods for comparison, see CBO_Hermitian_projection.m and CBO_Unitary_projection.m

The four CBO methods work for fixed dimensional optimization problems. Simulations for various dimensions is needed to solve the entanglement problem, see M_wrapper.m 

We also give a multi-species CBO method which allows information exchange of consensus between different dimensions, see Multispecies_CBO_Hermitian.m 

- 'test/'
  
The scripts are able to produce the most important numerical results in the paper.

test_conv.m: Figs. 2, 9

example_states.m: Figs. 3, 8

comparison.m: Figs. 4, 5, 6, 7, 10, 11


