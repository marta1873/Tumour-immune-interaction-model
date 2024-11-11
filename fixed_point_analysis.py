import numpy as np
import matplotlib.pyplot as plt
from model import create_u_vector, matrix_exponential_dist, g


def calculate_interaction_matrix(u_vector, v_vector, mu_c, theta_c, zeta_c, eta, mu_t, theta_t, zeta_t, gamma_mat):
    '''
    Creates interaction matrix M containing intra-population competition death rates and binding growth/death rates
    Args
    u_vector: array containing points in lattice
    v_vector: : array containing points in lattice
    mu_c: float, cancer cell intra-population competition death rate
    theta_c: float, intra-population competition range of cancer cells
    zeta_c: float, cancer cell binding death rate
    eta: float, binding distance between cancer cells and T-cells (set to eta=2L unless implementing original version of model)
    mu_t: float, T-cell intra-population competition death rate
    theta_t: float, intra-population competition range of T-cells
    zeta_t: float, T-cell binding birth rate
    gamma_mat: 2-dimensional array, (i,j)th element contains binding affinity between ith cancer cell and jth T-cell
    '''
    # Get number of points in lattice
    l_u = len(u_vector)
    l_v = len(v_vector)

    # Create M matrix
    M = np.zeros((2*l_u, 2*l_v))
    for i in range(l_u):
        # M_11
        M[i, :l_u] = mu_c * g(u_vector[i], u_vector, theta_c)
        # M_12
        M[i, l_u:] = zeta_c * np.multiply(gamma_mat[i], g(u_vector[i], v_vector, eta))
        # - M_21
        M[l_u + i, :l_u] = - zeta_t * np.multiply(gamma_mat[:, i], g(v_vector[i], u_vector, eta))
        # M_22
        M[l_u + i, l_u:] = mu_t * g(v_vector[i], v_vector, theta_t)

    return M

def calculate_fixed_point(matrix_m, alpha_c, alpha_t, u_vector, v_vector):
    '''
    Calculates non-trivial fixed point in simplified system (beta_c=0) if it exists
    Args
    matrix_m: 2d array, interaction matrix containing intra-population competition death rates and binding growth/death rates
    alpha_c: float, cancer cell division rate
    alpha_t: float, T-cell division rate
    u_vector: array containing points in lattice
    v_vector: : array containing points in lattice
    Returns
    n_c: array, cancer cell density at fixed point
    n_t: array, T-cell density at fixed point
    '''
    # Create vector bar{alpha} = (alpha_c, ..., alpha_c, alpha_t, ..., alpha_t)
    alpha = np.zeros(len(u_vector) + len(v_vector))
    alpha[: len(u_vector)] = alpha_c * np.ones(len(u_vector))
    alpha[len(u_vector):] = alpha_t * np.ones(len(v_vector))

    # Get lattice step
    step = u_vector[1] - u_vector[0]

    # Only returns fixed point if M^{-1} exists
    if np.linalg.det(matrix_m) != 0:
        # Get fixed point values
        n_vec = np.matmul(np.linalg.inv(matrix_m), alpha)/step

        # Split into cancer cells and T-cells
        n_c = n_vec[: len(u_vector)]
        n_t = n_vec[len(u_vector) :]

        return n_c, n_t


def plot_eigenvalues_jacobian(iterations=15, L=1, num_points_lattice=10, mu_c=5 * 10 ** (-6), theta_c=1.8,
                              zeta_c=5 * 10 ** (-6), eta=2, mu_t=5 * 10 ** (-6), theta_t=1.8, zeta_t=3 * 10 ** (-5),
                              alpha_c=1.5, alpha_t=0.05):
    '''
    Generates a plot of the eigenvalues of the Jacobian for several instances of the interaction matrix gamma
    Args
    iterations: int, number of systems plotted
    L: float, (positive) size of interval determined by [-L,L]
    num_points_lattice: int, number of points in lattice
    mu_c: float, cancer cell intra-population competition death rate
    theta_c: float, intra-population competition range of cancer cells
    zeta_c: float, cancer cell binding death rate
    eta: float, binding distance between cancer cells and T-cells (set to eta=2L unless implementing original version of model)
    mu_t: float, T-cell intra-population competition death rate
    theta_t: float, intra-population competition range of T-cells
    zeta_t: float, T-cell binding birth rate
    alpha_c: float, cancer cell division rate
    alpha_t: float, T-cell division rate
    '''
    # Create u and v vectors
    spatial_step = 2 * L / num_points_lattice
    u_vector = create_u_vector(L, num_points_lattice, spatial_step)
    v_vector = u_vector.copy()

    fig, ax = plt.subplots(1)
    fig.suptitle('Eigenvalues of Jacobian at non-trivial fixed point', fontsize=18)

    for i in range(iterations):
        # Generate gamma matrix
        gamma_mat = matrix_exponential_dist(10, 1)
        # Get interaction matrix M
        m = calculate_interaction_matrix(u_vector, v_vector, mu_c, theta_c, zeta_c, eta, mu_t, theta_t, zeta_t, gamma_mat)
        # Get non-trivial fixed point if it exists
        if calculate_fixed_point(m, alpha_c, alpha_t, u_vector, v_vector):
            nC_fixed, nT_fixed = calculate_fixed_point(m, alpha_c, alpha_t, u_vector, v_vector)
            n_vec = np.concatenate((nC_fixed, nT_fixed))

            # Get Jacobian matrix (first term at non-trivial fixed point is 0)
            jacobian_matrix = - np.matmul(np.diag(n_vec), m) * spatial_step

            # Get eigenvalues
            eigenvalues = np.linalg.eig(jacobian_matrix)[0]

            # Plot eigenvalues
            ax.scatter(eigenvalues.real, eigenvalues.imag, s=10)

    ax.set_xlabel('Real axis', fontsize=16)
    ax.set_ylabel('Imaginary axis', fontsize=16)
    plt.grid()
    plt.show()
