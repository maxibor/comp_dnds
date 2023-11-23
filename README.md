# comp_dnds

Efficiently computing dN/dS according to the Nei-Gojobori(1986)[^1] method


```python
from comp_dnds import dnds

d = dnds()
ref_seq = "ATG AAA CCC GGG TTT TAA".replace(" ", "")
obs_seq = "ATG AAA CGC GGC TAC TAA".replace(" ", "")
dn, ds = d.compute(ref_seq, obs_seq)
ω = dn/ds
print(ω)
# 0.12468371343300778
```

The dN/dS ratio, often represented as ω (omega), is a metric used in molecular biology and evolutionary biology to measure the selective pressure acting on a protein-coding gene. It compares the rate of non-synonymous substitutions (dN) to the rate of synonymous substitutions (dS) in coding sequences.

**Non-synonymous substitutions (dN)**: These are mutations that change the amino acid in a protein. They can potentially affect the protein's structure and function, and therefore they can be subject to natural selection.

**Synonymous substitutions (dS)**: These are mutations that do not change the amino acid in a protein. Because they don't change the protein's structure or function, they are often considered to be selectively neutral, or at least under less selective pressure than nonsynonymous changes.

> The dN/dS ratio provides insight into the evolutionary forces acting on a gene:

- **dN/dS = 1**: The rate of nonsynonymous substitutions is about the same as the rate of synonymous substitutions. This suggests that the gene is evolving under neutral evolution. 

- **dN/dS < 1**: The rate of non-synonymous changes is lower than the rate of synonymous changes. This suggests that the gene is under purifying or negative selection. Negative selection acts to remove deleterious mutations. 

- **dN/dS > 1**: The rate of non-synonimous changes is higher than the rate of synonymous changes. This suggests that the gene is under positive or adaptive selection. This means that nonsynonymous changes provide some advantage and are being selected for.



[^1]: [Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions](https://doi.org/10.1093/oxfordjournals.molbev.a040410)