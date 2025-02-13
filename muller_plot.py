import numpy as np
import pandas as pd
from pyfish import fish_plot, process_data, setup_figure
import matplotlib.pyplot as plt


def create_population_df(matrix, time_vector):
    '''
    Creates dataframe of population in format to create muller plots
    Args
    matrix: list of arrays, each array contains the cancer cell densities across the lattice at a time-step
    time_vector: array containing all the time points at which the densities are evaluated
    Returns
    populations_df: pandas dataframe, containing columns "Id" (index of lattice point), "Step" (timestep), "Pop"
    (population density at lattice point "Id" for time "Step").
    '''
    # Get number of lattice points
    num_lattice_points = len(matrix[0])

    # Get number of time steps
    num_time_steps = len(matrix)

    # Create array to store the evolution of the population density over time for each lattice point
    pop_list = np.zeros((num_time_steps * num_lattice_points, 3))

    for i in range(num_lattice_points):
        # For each lattice point create an array containing one row per iteration step and three columns
        index_group = np.zeros((num_time_steps, 3))

        # First column contains id of lattice point (ranges from 1 to num_lattice_points)
        index_group[:, 0] = (i + 1) * np.ones(num_time_steps)

        # Second column contains time steps
        index_group[:, 1] = time_vector

        # Third column contains the population density (from amtrix) of that particular lattice point
        index_group[:, 2] = matrix[:, i]

        # Store information for lattice point
        pop_list[num_time_steps * i: num_time_steps * (i + 1), :] = index_group

    # Create dataframe
    populations_df = pd.DataFrame(pop_list, columns=["Id", "Step", "Pop"])
    return populations_df


def create_parent_df(matrix):
    '''
    Creates dataframe of parent-children cells for Muller plot
    Args
    matrix: list of arrays, each array contains the cancer cell densities across the lattice at a time-step
    Returns
    parent_tree_df: pandas dataframe, containing columns "ParentId" (id of parent cell), "ChildId (id of child cells)" -
    all cells in the lattice are defined to have a parent cell with id 0, cells in the lattice are then numbered in
    increasing order starting at 1
    '''
    # Get number of lattice points
    num_lattice_points = len(matrix[0])

    # Create array to store id's (first column will be kept to 0)
    parent_tree = np.zeros((num_lattice_points, 2))

    # Get Id's of lattice points (ranges from 1 to num_lattice_points)
    id = np.arange(1, num_lattice_points + 1)

    # Set Id of cells in lattice
    parent_tree[:, 1] = id

    # Create dataframe
    parent_tree_df = pd.DataFrame(parent_tree, columns=["ParentId", "ChildId"])
    return parent_tree_df


def muller_plot(nC, time_vector):
    '''
    Args
    nC: list of arrays, each array contains the cancer cell densities across the lattice at a time-step
    time_vector: array containing all the time points at which the densities are evaluated
    '''
    # Convert list of arrays to 2-d array
    nC = np.array(nC)

    # Create dataframe with evolution of densities over time
    population_cancer = create_population_df(nC, time_vector)

    # Create dataframe with child-parent id's
    parent_tree_cancer = create_parent_df(nC)

    # Plot figure
    data_cancer = process_data(population_cancer, parent_tree_cancer, absolute=True, smooth=1)
    setup_figure()
    fish_plot(*data_cancer)
    plt.show()
