from Bio.Data import CodonTable
from typing import List, Tuple
import numpy as np
from comp_dnds.paths import get_all_codons_paths, hamming_distance
from comp_dnds.expected import get_expected
from tqdm import tqdm
import random
from scipy.stats import norm


class dnds:
    """Compute dn and ds for two sequences using the Nei-Gojobori(1986) method."""

    def __init__(
        self,
        alphabet: List[str] = ["A", "T", "G", "C"],
        codon_size: int = 3,
        genetic_code: int = 1,
    ) -> None:
        """Constructor for the dnds class.

        Args:
            alphabet (List[str], optional): Letters of DNA alphabet. Defaults to ['A','T','G','C'].
            codon_size (int, optional): Number of DNA base pairs in a codon. Defaults to 3.
            genetic_code (int, optional): Identifier of the genetic code.
                See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for reference.
                Defaults to 1.
        """
        codon_table = CodonTable.generic_by_id[genetic_code].forward_table
        codon_table.update({k: "stop" for k in CodonTable.generic_by_id[1].stop_codons})
        self.codon_table = {
            k: v
            for k, v in codon_table.items()
            if all(letter in alphabet for letter in k)
        }

        self.alphabet = alphabet
        self.codon_size = codon_size

        self.possible_paths = get_all_codons_paths(alphabet, codon_size)
        self.syn_expect, self.non_syn_expect = get_expected(self.codon_table)

    def __is_synonymous(self, ref_codon: str, obs_codon: str) -> bool:
        """Private method to check if the mutation between two codons is synonymous.
        Args:
            ref_codon (str): reference codon
            obs_codon (str): observed codon
        Returns:
            bool: True if the mutation is synonymous, False otherwise.
        """
        if self.codon_table[ref_codon] == self.codon_table[obs_codon]:
            return True
        return False

    def __is_in_frame(self, seq):
        if len(seq) % self.codon_size == 0:
            return True
        return False

    def __get_codons(self, seq: str) -> List[str]:
        """Private method to get codons from a sequence.

        Args:
            seq (str): DNA sequence

        Returns:
            List[str]: list of codons
        """
        codons = [
            seq[i : i + self.codon_size] for i in range(0, len(seq), self.codon_size)
        ]
        for codon in codons:
            for base in codon:
                if base not in self.alphabet:
                    raise ValueError(f"Invalid base {base} in codon {codon}")
        return codons

    def __comp_obs_sd_nd_single_codon(
        self, ref_codon: str, obs_codon: str
    ) -> Tuple[float, float]:
        """Private method to compute the observed number of synonymous and non-synonymous mutations for a single codon pair.

        Args:
            ref_codon (str): reference codon
            obs_codon (str): observed codon
        Returns:
            Tuple[float, float]: Number of observed synonymous and non-synonymous mutations.
        """

        _s = 0
        _n = 0
        if ref_codon == obs_codon:
            return _s, _n

        dist = hamming_distance(ref_codon, obs_codon)
        if dist == 1:
            if self.__is_synonymous(ref_codon, obs_codon):
                _s = 1
            else:
                _n = 1
        else:
            paths = self.possible_paths[(ref_codon, obs_codon)]
            situations = len(paths)
            for path in paths:
                for p in range(len(path) - 1):
                    if self.__is_synonymous(path[p], path[p + 1]):
                        _s += 1
                    else:
                        _n += 1
            _s = _s / situations
            _n = _n / situations
        return _s, _n

    def __comp_obs_sd_nd(
        self, ref_codons: List[str], obs_codons: List[str]
    ) -> Tuple[float, float]:
        """Private method to compute the observed number of synonymous and non-synonymous mutations.

        Args:
            ref_codons (List[str]): codons from the reference sequence
            obs_codons (List[str]): codons from the observed sequence

        Returns:
            Tuple[float, float]: Number of observed synonymous and non-synonymous mutations.
        """
        syn = 0
        non_syn = 0
        for ref_codon, obs_codon in zip(ref_codons, obs_codons):
            _s, _n = self.__comp_obs_sd_nd_single_codon(ref_codon, obs_codon)
            syn += _s
            non_syn += _n
        return syn, non_syn

    def __comp_exp_s_n(self, ref_codons: List[str]) -> Tuple[float, float]:
        """Compute the expected number of synonymous and non-synonymous mutations for the reference codons.\
        
        Args:
            ref_codons (List[str]): codons from the reference sequence
        Returns:
            Tuple[float, float]: Number of expected synonymous and non-synonymous mutations.
        """
        exp_n = 0
        for codon in ref_codons:
            try:
                exp_n += self.non_syn_expect[codon]
            except KeyError:
                raise KeyError(f"Invalid codon: {codon}")
        exp_s = len(ref_codons) * self.codon_size - exp_n

        return exp_s, exp_n

    def _compute_no_bootstrap(
        self, ref_codons: List[str], obs_codons: List[str]
    ) -> Tuple[float, float]:
        """Compute dn and ds for two sequences using the Nei-Gojobori(1986) method.

        Args:
            ref_seq (str): Reference sequence. Typically the ancestral sequence.
            obs_seq (str): Observed sequence. Typically the derived sequence.

        Returns:
            Tuple[float, float]: Respectively dN and dS values.
        """

        exp_s, exp_n = self.__comp_exp_s_n(ref_codons)
        obs_s, obs_n = self.__comp_obs_sd_nd(ref_codons, obs_codons)

        pn = obs_n / exp_n
        ps = obs_s / exp_s

        dn = -(3 / 4) * np.log(1 - (4 / 3) * pn)
        ds = -(3 / 4) * np.log(1 - (4 / 3) * ps)

        return dn, ds

    def _compute_bootstrap(
        self, ref_codons: List[str], obs_codons: List[str], replicates: int = 1000
    ) -> Tuple[float, float]:
        """Compute dn and ds for two sequences using the Nei-Gojobori(1986) method, and uses bootstrap resampling to compute significance.

        Args:
            ref_seq (str): Reference sequence. Typically the ancestral sequence.
            obs_seq (str): Observed sequence. Typically the derived sequence.
            replicates (int, optional): Number of bootstrap replicates. Defaults to 1000.

        Returns:
            Tuple[float, float, float]: Respectively dN and dS, and z values.
        """

        dN, dS = self._compute_no_bootstrap(ref_codons, obs_codons)

        codons_idxs = [
            random.choices(range(len(ref_codons)), k=len(ref_codons))
            for i in range(replicates)
        ]
        rep_ref_codons = [[ref_codons[i] for i in idx] for idx in codons_idxs]
        rep_obs_codons = [[obs_codons[i] for i in idx] for idx in codons_idxs]

        dNs = np.array([])
        dSs = np.array([])

        for r in tqdm(range(replicates)):
            dn, ds = self._compute_no_bootstrap(rep_ref_codons[r], rep_obs_codons[r])
            dNs = np.append(dNs, dn)
            dSs = np.append(dSs, ds)

        z = (dN - dS) / np.sqrt(np.nanvar(dSs) + np.nanvar(dNs))

        return dN, dS, z

    def compute(
        self, ref_seq: str, obs_seq: str, bootstrap: int = 0
    ) -> Tuple[float, float]:
        """Compute dn and ds for two sequences using the Nei-Gojobori(1986) method.

        Args:
            ref_seq (str): Reference sequence. Typically the ancestral sequence.
            obs_seq (str): Observed sequence. Typically the derived sequence.
            bootstrap (int, optional): Number of bootstrap replicates. Defaults to 0.

        Returns:
            Tuple[float, float, float, float]: Respectively dN and dS, and Z-score, and p-value.
            Z-score and p-value are computed only if bootstrap >= 100.
        """
        if len(ref_seq) != len(obs_seq):
            raise ValueError(
                "Reference and observed sequences must be of the same length."
            )
        if not self.__is_in_frame(ref_seq):
            raise ValueError(f"Sequence length must be a multiple {self.codon_size}")

        ref_codons = self.__get_codons(ref_seq)
        obs_codons = self.__get_codons(obs_seq)

        if bootstrap < 100:
            dn, ds = self._compute_no_bootstrap(ref_codons, obs_codons)
            z = None
            pval = None
        else:
            dn, ds, z = self._compute_bootstrap(
                ref_codons, obs_codons, replicates=bootstrap
            )
            pval = norm.sf(abs(z)) * 2

        return dn, ds, z, pval
