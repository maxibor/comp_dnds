from itertools import product, permutations
from typing import List, Dict, Tuple

# Define the alphabet and a helper function to compute Hamming distance
alphabet = ["A", "T", "G", "C"]


def hamming_distance(codon1, codon2):
    """Compute the Hamming distance between two sequences"""
    if len(codon1) != len(codon2):
        raise ValueError("Codons must have the same length")
    return sum(a != b for a, b in zip(codon1, codon2))


# Generate all possible codons


# Precompute the paths for all codon pairs with a Hamming distance of up to 3


def get_all_codons_paths(
    alphabet: List[str] = ["A", "T", "G", "C"], codon_size: int = 3
) -> Dict[Tuple[str, str], List[List[str]]]:
    """Get all possible mutation paths between two codons.
    If there is more than one mutation between the two codons, there are multiple possible ways
    to get from one codon to the other. This function computes all possible paths between two codons.

    Args:
        alphabet (list, optional): Letters of DNA alphabet. Defaults to ['A','T','G','C'].
        codon_size (int, optional): Number of DNA base pairs in a codon. Defaults to 3.

    Returns:
        Dict[Tuple[str, str], List[List[str]]]: A dictionary with codon pairs as keys and a
            list of all possible paths between the two codons as values.
    """
    all_codons = ["".join(c) for c in product(alphabet, repeat=codon_size)]
    precomputed_paths = {}
    for codon1 in all_codons:
        for codon2 in all_codons:
            if 1 <= hamming_distance(codon1, codon2) <= codon_size:
                differing_positions = [
                    i for i, (a, b) in enumerate(zip(codon1, codon2)) if a != b
                ]
                paths = []
                for perm in permutations(differing_positions):
                    current_codon = list(codon1)
                    path = [codon1]
                    for position in perm:
                        current_codon[position] = codon2[position]
                        path.append("".join(current_codon))
                    paths.append(path)
                precomputed_paths[(codon1, codon2)] = paths
    return precomputed_paths
