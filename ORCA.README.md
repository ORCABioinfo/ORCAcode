ORCAcode
========

Example for ORCA on drug sensitivity data.
- drug.sensitivity : R object (data frame) containing GI50 (growth
inhibition at 50%) values of 116 chemotherapeutic agents on NCI-60 cell
panel
- no.member : R object (vector of number) containing number of drugs in
each drug group
- orca.functions.R : functions necessary for ORCA
- ORCA.perm.test.R : main file for consolidating the actual and
permutation of ORCA on drug sensitivity dataset
- permutation.R : script for performing permutation on drug sensitivity
dataset
- phyperg.actual.R : script for performing calculation of
hypergeometric p-values on drug sensitivity dataset
