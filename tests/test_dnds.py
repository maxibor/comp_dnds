from comp_dnds.dnds import dnds
from comp_dnds.expected import get_expected
from comp_dnds.paths import get_all_codons_paths


import pytest


class TestDnds:
    d = dnds()
    ref_seq = "ATG AAA CCC GGG TTT TAA".replace(" ", "")
    obs_seq = "ATG AAA CGC GGC TAC TAA".replace(" ", "")

    def test_class_instantiation(self):
        assert isinstance(self.d, dnds)
        assert self.d.codon_size == 3
        for codon in self.d.codon_table:
            assert "U" not in codon
            assert len(codon) == self.d.codon_size
        assert "stop" in self.d.codon_table.values()
        assert self.d.possible_paths == get_all_codons_paths(
            alphabet=self.d.alphabet, codon_size=self.d.codon_size
        )
        assert self.d.syn_expect == get_expected(self.d.codon_table)[0]

    def test_is_synonymous(self):
        assert self.d._dnds__is_synonymous("TTA", "TTA") is True
        assert self.d._dnds__is_synonymous("TTA", "TTG") is True
        assert self.d._dnds__is_synonymous("TTA", "TTT") is False

    def test_is_in_frame(self):
        assert self.d._dnds__is_in_frame("TTATTT") is True
        assert self.d._dnds__is_in_frame("TTATTTT") is False

    def test_get_codons(self):
        assert self.d._dnds__get_codons("TTATTT") == ["TTA", "TTT"]

    def test_comp_exp_s_n(self):
        seq = "TTATTT"
        codons = self.d._dnds__get_codons(seq)
        syn, non_syn = get_expected(
            codon_table=self.d.codon_table, alphabet=self.d.alphabet
        )
        assert self.d._dnds__comp_exp_s_n(codons)[1] == sum(
            non_syn[codon] for codon in codons
        )

        syn, non_syn = self.d._dnds__comp_exp_s_n(
            self.d._dnds__get_codons(self.ref_seq)
        )

        assert 3.333 == pytest.approx(syn, 0.01)
        assert 14.666 == pytest.approx(non_syn, 0.01)

    def test_comp_obs_sd_nd_single_codon(self):
        assert self.d._dnds__comp_obs_sd_nd_single_codon("ATG", "ATG") == (0, 0)
        assert self.d._dnds__comp_obs_sd_nd_single_codon("CCC", "CGC") == (0, 1)
        assert self.d._dnds__comp_obs_sd_nd_single_codon("TTT", "TAC") == (1, 1)
        assert self.d._dnds__comp_obs_sd_nd_single_codon("GGG", "GGC") == (1, 0)

    def test_comp_obs_sd_nd(self):
        ref_codons = self.d._dnds__get_codons(self.ref_seq)
        obs_codons = self.d._dnds__get_codons(self.obs_seq)
        assert self.d._dnds__comp_obs_sd_nd(ref_codons, obs_codons) == (2, 2)

    def test_compute_different_length(self):
        with pytest.raises(ValueError):
            self.d.compute("ATG", "ATGAA")

    def test_compute_non_string(self):
        with pytest.raises(TypeError):
            self.d.compute(123, "ATG")

    def test_compute_not_alphabet(self):
        with pytest.raises(ValueError):
            self.d.compute("ATG", "AXG")

    def test_compute(self):
        result = self.d.compute(self.ref_seq, self.obs_seq)
        assert isinstance(result, tuple)
        assert all(isinstance(x, float) for x in result[:2])
        assert result[:2] == pytest.approx((0.150, 1.207), 0.01)

    def test_compute_bootstrap(self):
        ref_seq = "AGT CGG AAC ATG AGA AGA AGG AAT CAC ATC ACG TCG TTA GGC AAG".replace(
            " ", ""
        )
        obs_seq = "AGT CGG AAC ATA AGA AGA AGG AAT CAG TTC ACT TCG TTA GGC AAG".replace(
            " ", ""
        )
        result = self.d.compute(ref_seq, obs_seq, bootstrap=1000)
        assert isinstance(result, tuple)
        dn, ds, z, p = result
        assert isinstance(dn, float)
        assert isinstance(ds, float)
        assert isinstance(z, float)
        assert isinstance(p, float)
        assert dn == pytest.approx(0.089, 0.1)
        assert ds == pytest.approx(0.11, 0.1)
        assert z == pytest.approx(-0.19, 0.2)
        assert p == pytest.approx(0.84, 0.1)
