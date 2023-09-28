# Tumour-immune-interaction-model

## Run simulation of population densities

The code provided in model.py is a simulation of the interaction between the cancer and T-cell population densities within a tumour. In particular, the model focuses on the effect of binding between cancer cells and T-cells (which increases the death rate of cancer cells and the birth rate of T-cells). The code provided implements a discretisation of an adaptation of a set of PDEs proposed in [[1]](#1). 

The function interaction_matrix_model implements the model by taking as inputs the following parameter values:
- $L$: the interval $[-L,L]$ defines the space where we model the densities
- $l$: number of points in the lattice used to discretise $[-L,L]$
- $\tau$: time-step
- $t_f$: final time of simulation
- $a$ and $A$: used to define initial conditions via $n^0_C(u_i) = 10^4(1 + a \cos (Au_i))$ and $n^0_T(v_j) = 10^4(2 + a \cos (Av_j))$
- $\theta_C$: intra-population competition range of cancer cells
- $\theta_T$: intra-population competition ranges of T-cells
- $\eta$: binding distance between cancer cells and T-cells in $[-L,L]$
- $\alpha_C$: cancer cell division rate
- $\mu_C$: intra-population competition death rate of cancer cells
- $\zeta_C$: cancer cell death rate due to binding
- $\alpha_T$: T-cell division rate
- $\mu_T$: intra-population competition death rate of T-cells
- $\zeta_T$: T-cell birth rate due to binding
- $\Gamma$: interaction matrix of binding affinities between cancer cells and T-cells
- $\lambda_C$: relating to how fast cancer cells can mutate across $[-L,L]$

The function returns:
- nC_matrix: list of arrays of the cancer cell densities at each phenotype for each time-step
- nT_matrix: list of arrays of the T-cell densities at each phenotype for each time-step
- u_vector: an array containing $(u_1, u_2, \dotsc, u_l)$ (which are the lattice points of the discretisation of $[-L,L]$)
- time_vector: list containing the different times at which the model is evaluated.

The value of $\eta$ should be kept fixed to $\eta = 2L$ unless using the original version of the model, where not all cancer cells and T-cells will bind.

The function matrix_exponential_dist can be used to generate the elements of the interaction matrix $\Gamma$ via an exponential distribution. A suggestion of parameter for interaction_matrix_model is included to run the model at the end of the file. 

## References
<a id="1">[1]</a> 
Luis Almeida, Chloe Audebert, Emma Leschiera, and Tommaso Lorenzi. Discrete and continuum models for the coevolutionary dynamics between CD8+ cytotoxic T lymphocytes and tumour cells. Mathematical Medicine and Biology: A Journal of the IMA, 40(2):141–174, 01 2023. ISSN 1477-8602. doi: 10.1093/imammb/dqac017. URL https://doi.org/10.1093/imammb/dqac017
