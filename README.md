# MrRobin: Two-sample Mendelian Randomization method ROBust to correlated and some INvalid instruments

The goal of `MrRobin` is to provide tools to conduct mendelian randomization analysis
using the MR-Robin algorithm. MR-Robin is a two-sample Mendelian Randomization method ROBust to
correlated and some INvalid instruments that takes summary statistics from complex trait GWAS
and multi-tissue eQTL analyses as input and
uses a reverse regression random slope mixed model to infer whether a gene is
associated with a complex trait.

## Setup

To install and load functions from `MrRobin`, run the following:

  ```R
  devtools::install_github("kjgleason/MrRobin")
  library("MrRobin")
  ```

## Citation

To cite `MrRobin` in publications, please use:

Fan Yang, Kevin J. Gleason, Jiebiao Wang, The GTEx consortium, Jubao Duan, Xin He, Brandon L. Pierce and Lin S. Chen. CCmed: cross-condition mediation analysis for identifying robust trans-eQTLs and assessing their effects on human traits. bioRxiv (2019), doi:10.1101/803106.
