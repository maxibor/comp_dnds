from Bio.Data import CodonTable
from typing import List, Dict, Tuple


def get_expected(
    codon_table: Dict[str, str], alphabet: List[str] = ["A", "T", "G", "C"]
) -> Tuple[Dict[str, float], Dict[str, float]]:
    """Calculate the expected number of synonymous and non synomymous mutations for each codon.

    Args:
        codon_table (Dict[str, str]): genetic code
        alphabet (List[str], optional): Possible DNA letters. Defaults to ['A','T','G','C'].

    Returns:
        Tuple[Dict[str, float], Dict[str, float]] : expected probability of synonymous and non synonymous mutations for each codon.
    """
    non_syn_expect = {k: 0 for k in codon_table.keys()}
    for codon in codon_table:
        aa = codon_table[codon]
        for base in codon:
            if base not in alphabet:
                raise ValueError(f"Invalid base {base} in codon {codon}")
        else:
            for pos, base in enumerate(codon):
                for letter in alphabet:
                    if letter != base:
                        _ = list(codon)
                        _[pos] = letter
                        mut = "".join(_)
                        try:
                            if codon_table[mut] != aa:
                                non_syn_expect[codon] += 1 / 3
                        except KeyError:
                            non_syn_expect[codon] += 1 / 3
    syn_expect = {k: 3 - v for (k, v) in non_syn_expect.items()}
    return (syn_expect, non_syn_expect)
