# EigenvaluesMultilinearFormsGraph
### Introduction
This repository contains code, data and additional material associated to the work Eimear Byrne and I have done on the multilinear froms graph eigenvalues and its applications to tensor-codes bounds. The functions are directly related to the work done in the article mentionned below. Namely, there are MAGMA files that compute the spectrum of the multilinear forms graph using its dual characterisation and you will also find functions that accompanies the article heavy computations for the case of 2x3x3 tensors.

### Related article
_Article:_ **The multilinear forms Cayley graph and the eigenvalue method for tensor codes**, by Eimear Byrne & Lucien François. \
_Status of the paper:_ Submitted to SIDMA (SIAM), under review. Open-source arXiv TBD.

### Organisation note
*Note:* The file `MAGMAFUNCTIONS.m` is a file containing a lot of functions (linear algebra, combinatorics, tensor-algebra, finite fields). It is required to run functions in the file `MultilinearFormsSpectra.m`. If you want to use the latter, download the former and put it into the same folder (or copy-paste the required functions). A documentation for the functions in the file `MAGMAFUNCTIONS.m` is avaliable : [`MAGMAFUNCTIONS_documentation.pdf`](MAGMAFUNCTIONS_documentation.pdf). *doc currently in construction...* 


### AI disclamer
Part of the functions present in this repository have been written with the help of the Copilot autocompletion on VsCode. 


## Assistance in thorough computations for 2x3x3 tensors
A section of the paper focuses on the computation of the 2x3x3 trilinear forms graph spectrum. The document contains thorough computations for the stabiliser sizes of the action of GL(2,q) x GL(3,q) on the subspaces of 2x3 matrices over Fq, the orbit sizes of the action of GL(2,q) x GL(3,q) x GL(3,q) on the set of 2x3x3 tensors over Fq, the spectrum computation and the ratio-type bound. 
In the folder [`Orbits233`](Orbits233/) you will find magma programs that give more details on those computations using integer/rational polynomial expressions (of variable q).

## Computation of the multilinear graph spectrum

### Data files organisation
The computed spectrums are stored in the three files below as follows. For any integer `i` in the number of lines in the csv below, the `i`-th line of the `SUMMARY` file provides the parameters q,n1,...,nm of the tensor space *(of size n1 x ... x nm over Fq)* such that the `i`-th line of the `EIGENVALUE` file contains the different eigenvalues of the multilinear forms graph spectrum for such parameters, and whose multiplicities can be found (with the same order) in the `MULTIPLICITIES` file. More precisely, if that spectrum has been computed, there exists a line $i$ such that the spectrum of the multilinear forms graph $Cay(\mathbb F_q^{n_1 \times ... \times n_m}, \mathbb S)$
+ The file [`DATASET_EIGENVALUES_KnownSpectra.csv`](DATASET_EIGENVALUES_KnownSpectra.csv/) contains lines of the form $\theta_0$, $\theta_1$, ... , $\theta_r$
+ The file [`DATASET_MULTIPLICITIES_KnownSpectra.csv`](DATASET_MULTIPLICITIES_KnownSpectra.csv/) contains lines of the form $m_0$, $m_1$, .... , $m_r$
+ The file [`DATASET_SUMMARY_KnownSpectra.csv`](DATASET_SUMMARY_KnownSpectra.csv/) contains lines of the form $q$, $n_1$, ... , $n_m$.

