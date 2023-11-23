from comp_dnds.paths import hamming_distance, get_all_codons_paths
import pytest


class TestHammingDistance:
    # Compute the Hamming distance between two identical codons, which should return 0.
    def test_identical_codons(self):
        codon1 = "AAA"
        codon2 = "AAA"
        assert hamming_distance(codon1, codon2) == 0

    # Compute the Hamming distance between two codons that differ in only one position, which should return 1.
    def test_differ_in_one_position(self):
        codon1 = "AAA"
        codon2 = "AAT"
        assert hamming_distance(codon1, codon2) == 1

    # Compute the Hamming distance between two codons that differ in all positions, which should return the length of the codons.
    def test_differ_in_all_positions(self):
        codon1 = "AAA"
        codon2 = "TTT"
        assert hamming_distance(codon1, codon2) == len(codon1)

    # Compute the Hamming distance between two codons that are not strings, which should raise a TypeError.
    def test_non_string_codons(self):
        codon1 = 123
        codon2 = "AAA"
        with pytest.raises(TypeError):
            hamming_distance(codon1, codon2)

    # Compute the Hamming distance between two codons with lengths greater than 3, which should raise a ValueError.
    def test_long_codons(self):
        codon1 = "AAA"
        codon2 = "AAAA"
        with pytest.raises(ValueError):
            hamming_distance(codon1, codon2)


class TestGetAllCodonsPaths:
    def test_returns_dictionary_with_paths(self):
        result = get_all_codons_paths()
        assert isinstance(result, dict)
        for key, value in result.items():
            assert isinstance(key, tuple)
            assert isinstance(value, list)
            for path in value:
                assert isinstance(path, list)
                for codon in path:
                    assert isinstance(codon, str)

    def test_paths(self):
        alphabet = ["A", "T"]
        codon_size = 2
        target = {
            ("AA", "AT"): [["AA", "AT"]],
            ("AA", "TA"): [["AA", "TA"]],
            ("AA", "TT"): [["AA", "TA", "TT"], ["AA", "AT", "TT"]],
            ("AT", "AA"): [["AT", "AA"]],
            ("AT", "TA"): [["AT", "TT", "TA"], ["AT", "AA", "TA"]],
            ("AT", "TT"): [["AT", "TT"]],
            ("TA", "AA"): [["TA", "AA"]],
            ("TA", "AT"): [["TA", "AA", "AT"], ["TA", "TT", "AT"]],
            ("TA", "TT"): [["TA", "TT"]],
            ("TT", "AA"): [["TT", "AT", "AA"], ["TT", "TA", "AA"]],
            ("TT", "AT"): [["TT", "AT"]],
            ("TT", "TA"): [["TT", "TA"]],
        }
        assert get_all_codons_paths(alphabet, codon_size) == target
