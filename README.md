# structure-presevingCBO


## 1. Description

This repository contains the MATLAB codes used to generate the numerical results presented in the paper:

**Title:** Computation of entanglement for quantum states by a Consensus-Based Optimization method

**Authors:** Michael Herty, Yijia Tang, and Yizhou Zhou


## 2. Code Description

- 'src/'
The src files include several methods used to optimize the entanglement of formation for a given quantum state.
The two structure preserving CBO methods proposed in the paper: Hermitian preserving and Stiefel preserving CBO.
The two projection based CBO methods (Hermitian and Stiefel) for comparison.
The four CBO methods work for fixed dimensional optimization problems. 
We also give a multi-species CBO method which allows information exchange of consensus between different dimension.

- 'test/'
The scripts are able to produce the numerical results in the paper.
: Figs. 2, 9
example_states.m: Figs. 3, 8
comparison.m: Figs. 4, 5, 6, 7, 10, 11


