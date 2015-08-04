NaivePCA performs PCA for population genotype data. It implements the basic
algorithm as is described by [Patterson et al (2006)][1]. More precisely,
suppose we have *m* samples of ploidy *h* and *n* biallic markers. Let
*G*<sub>*ij*</sub> be the number of non-reference alleles for sample *i* at
marker *j*. NaivePCA computes:

```
\mu_j  = \sum_{i=1}^m G_{ij} / m
p_j    = \mu_j / h
M_{ij} = \frac{G_{ij}-\mu_j}{\sqrt{p_j(1-p_j)}}
X_{ij} = \sum_{k=1}^n M_{ik} M_{jk} / n
```
and finds the eigenvectors of matrix (*X*<sub>*ij*</sub>). Notably, if
*G*<sub>*ij*</sub> is missing data, *M*<sub>*ij*</sub> takes zero and the
computation of the mean needs to be adjusted as well.

The input of NaivePCA looks like:
```
sample1  110022110100202021001122*201
sample2  2012201102*221020211*1222001
```
where a number represents a genotype and other characters are treated as
missing data. For now, NaivePCA does not support real matrices. The output is
TAB-delimited. The first column is the sample name. The *i*-th column gives
the eigenvector for the (*i*-1)-th largest eigenvalue.

[1]: http://www.ncbi.nlm.nih.gov/pubmed/17194218
