# BeviMed 3.0

* Re-naming of parameters in `bevimed` function to match the names of variables in the paper (under submission). 
* The allele count matrix `G` should now be supplied as a matrix with rows corresponding to individuals, not variants.
* `expected_explained` and `explaining_variants` functions have been added, respectively computing the expected number of cases with their disease explained by the given variants, and expected number of pathogenic variants present amongst cases.
