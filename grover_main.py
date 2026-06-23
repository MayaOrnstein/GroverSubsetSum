from typing import List
from qiskit.quantum_info import Statevector

from grover_code_rca import grover_subsetsum_circuit

import time
import matplotlib.pyplot as plt
import numpy as np


def num_to_indices(n:int):
    """Helper function that converts a number into indices into the subset
    based on its binary representation
    """
    bit_string = f'{int(n):b}'[::-1]
    indices = []
    for j, bit in enumerate(bit_string):
        if bit == "1": indices.append(j)
    return indices


def run_grover_subset(val_set: List[int], target: int, verbose:int=1, plot_probs:bool=False):
    """Run the grover search algorithm on the subset sum problem

    :param val_set: The set we want to use for the subset sum problem
    :type val_set: List of integers

    :param target: The target value for the sum
    :type target: int

    :param verbose: How verbose we want this function to be
    :type verbose: int

    :param plot_probs: If True, then a plot of the output probabilities will be drawn
    :type plot_probs: bool

    ...

    :return: candidates (tuple of indices into the value set) and corresponding probabilitis
    :rtype: dict
    """
    # Size of the subset
    n = len(val_set)
    # Number of subsets (i.e. size of P(val_set))
    N = 2**n
    # Number of iterations Grover's algorithm requires
    K = int(np.rint(np.pi / (4 * np.arcsin(1 / np.sqrt(N))) - 1/2))

    # Evenly distributed if no solution (we double this proability to establish an arbitrary threshold for found solutions
    threshold = min(2/N, 0.5)

    # Create Grover Search circuit
    start_time = time.time()
    gsc = grover_subsetsum_circuit(vals=val_set, target=target, n_iter=K)
    build_time = time.time() - start_time
    if verbose > 1:
        print(gsc)

    # Compute the output on |00...0>
    start_time = time.time()
    psi = Statevector(gsc)
    compute_time = time.time() - start_time

    # Evaluate Probabilities
    # List of probablities arraged by bit string (i.e. index 5, 00101 indicates that the 0th index and 3rd index are used in the corresping sum)
    probabilities = np.round(psi.probabilities(), 3)[0:2**n]
    subset_probs = {}
    for idx, prob in enumerate(probabilities):
        indices = num_to_indices(idx)
        subset_probs[tuple(indices)] = float(prob)
    candidates = {subset: subset_probs[subset] for subset in subset_probs if subset_probs[subset] > threshold}

    # Plot the probabilites
    if plot_probs:
        plt.plot(np.arange(N), probabilities)
        for i in range(N):
            subset = [val_set[idx] for idx in num_to_indices(i)]
            if sum(subset) == target:
                plt.axvline(x=i, color='red', linestyle='--')
                plt.text(i, probabilities[i], subset, color='red')
        plt.xlabel(r"subset index")
        plt.ylabel("probabilities")
        plt.title(f"Prabability distribution of Grover Subset Sum for {val_set} with target {target}")

    # Compute the candidates for subset sum solutions and print them if desired
    if not candidates and verbose:
        print(f"No candidates found for a subset of {val_set} to sum to {target}.")
    else:
        for candidate in candidates:
            candidate_vals = [val_set[i] for i in candidate]
            if verbose:
                print(f"Candidate {candidate_vals} sum to {sum(candidate_vals)} and the target is {target}. p = {subset_probs[candidate]:02f}")

    # Print timing information if sufficiently verbose
    if verbose > 1:
        print(f"Circuit build time {build_time}s")
        print(f"Circuit evaluation time {compute_time}s")
    return candidates
