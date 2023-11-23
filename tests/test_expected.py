import pytest
from comp_dnds.expected import get_expected
from Bio.Data import CodonTable


class TestGetExpected:
    # Raises a ValueError when given a codon table with an invalid base in a codon.

    def test_return_type(self):
        codon_table = {"AAA": "K", "AAT": "N", "AAG": "K", "AGC": "N"}
        alphabet = ["A", "T", "G", "C"]
        syn_expect, non_syn_expect = get_expected(codon_table, alphabet)
        assert isinstance(syn_expect, dict)
        assert isinstance(non_syn_expect, dict)

    def test_invalid_base_in_codon(self):
        codon_table = {"AAA": "K", "AAT": "N", "AAG": "K", "AXC": "N"}
        alphabet = ["A", "T", "G", "C"]

        with pytest.raises(ValueError):
            get_expected(codon_table, alphabet)

    # Test return values of get_expected function.
    def test_get_expected(self):
        codon_table = CodonTable.generic_by_id[1].forward_table
        codon_table.update({k: "stop" for k in CodonTable.generic_by_id[1].stop_codons})
        alphabet = ["A", "T", "G", "C"]
        codon_table = {
            k: v
            for k, v in codon_table.items()
            if all(letter in alphabet for letter in k)
        }
        syn_expect, non_syn_expect = get_expected(codon_table)
        syn_target = {
            "ATG": 0,
            "AAA": 0.333,
            "CCC": 1,
            "GGG": 1,
            "TTT": 0.333,
            "TAA": 0.666,
        }
        non_syn_target = {
            "ATG": 3,
            "AAA": 2.667,
            "CCC": 2,
            "GGG": 2,
            "TTT": 2.667,
            "TAA": 2.333,
        }

        assert isinstance(syn_expect, dict)
        assert isinstance(non_syn_expect, dict)

        for key, value in syn_target.items():
            assert value == pytest.approx(syn_expect[key], 0.01)

        for key, value in non_syn_target.items():
            assert value == pytest.approx(non_syn_expect[key], 0.01)

        for key in syn_expect:
            assert 3 == pytest.approx(syn_expect[key] + non_syn_expect[key], 0.001)
