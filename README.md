
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GCLr <img src="man/figures/logo.png" align="right" height="100"/></a>

## Overview

GCLr is an R package designed to streamline genetic data analyses
performed by Gene Conservation Lab (GCL) staff. Many of the functions in
this package require data pulled directly from the GCL Oracle database
Loki. Because of this, only Alaska Department of Fish and Game staff
with credentials for accessing Loki will be able to use this package to
its full extent. However, the package does allow for users without Loki
credentials to read in genetic data contained in `GENEPOP` and `rubias`
formatted files and convert the data into the objects used that can be
used by GCLr (see `GCLr::genepop2gcl()` and `GCLr::base2gcl`).

Here are some examples of things that this package can be used for:

- **Laboratory workflow**

  - pull data reports directly from Loki into R

    - `get_asl_data()` gets age, sex, and length data

    - `get_extraction_info()` gets DNA extraction information

    - `get_geno()` gets raw genotypes

    - `get_gtseq_metadata` gets GTseq metadata

    - `get_tissue_data()` gets genetic tissue data

    - `get_tissue_locations()` gets location of archived tissues

  - create lab project sample sheets

    - `get_gtseq_sample_sheet()` creates a GTseq project sample sheet

  - break up Loki import files

    - `split_gtscore_loki_import()` splits a GTscore import file into
      multiple files

- **Laboratory quality control (QC)**

  - `qc_template` a template .RMD to compare the original (project) and
    reanalyzed (QC) sample genotypes to check for lab errors in data so
    they can be corrected, and calculates failure and error rates of
    genotypic data

- **Quality assurance (baseline and mixed stock analyses)**

  - `remove_ind_miss_loci()` removes samples with missing genotypes

  - `dupcheck_within_sillys()` identifies duplicated samples

  - `remove_dups()` removes one sample from each duplicate set

  - `find_alt_species()` checks for wrong species

  - `read_genepop_hwe()` reads in `GENEPOP` Hardy-Weinberg equilibrium
    test results

- **Baseline analysis**

  - Population structure

    - `collections_map()` creates an interactive map of collection
      locations

    - `fishers_test()` tests for homogeneity of allele frequencies,

    - `pool_collections()` combines collections into populations,

    - `read_genepop_dis()` reads in `GENEPOP` linkage disequilibrium
      test results

    - `summarize_LD()` summarizes `GENEPOP` linkage disequilibrium test
      results

    - `create_pwfst_tree()` creates a phylogenetic tree based on
      pairwise *F<sub>ST</sub>*

    - `create_mds_plot()` creates an interactive multidimensional
      scaling

    - `locus_stats()` calculates observed heterozygosity,
      *F<sub>IS</sub>*, and *F<sub>ST</sub>* by locus

  - Baseline evaluation

    - `create_rubias_base_eval()` creates `rubias` mixture and baseline
      objects/files for evaluating baseline reporting groups

    - `run_rubias_base_eval()` runs tests to evaluate the
      identifiability of baseline reporting units for genetic mixed
      stock analysis

    - `plot_baseline_eval()` plots the results of the baseline
      evaluation tests

- **Genetic mixed stock analysis (MSA)**

  - `create_rubias_base()` creates a `rubias` reference (aka “baseline”)
    object/file

  - `create_rubias_mix()` creates a `rubias` mixture object/file

  - `run_rubias_mix()` runs an MSA in `rubias`

  - `custom_comb_rubias_output()` summarizes the rubias MSA results

  - `stratified_estimator_rubias()` combines estimates for multiple
    strata into a single set of estimates

  - `summarize_rubias_individual_assign()` summarizes the rubias
    individual assignment results

- **Data conversion**

  - `gcl2fstat()` creates a genotypes file in `FSTAT` format

  - `gcl2nexus()` creates a genotypes file in `NEXUS` format

  - `gcl2genepop()` creates a genotypes file in `GENEPOP` format

  - `genepop2colony()` creates a genotypes file in `COLONY` format

  - `genepop2gcl()` reads in a `GENEPOP` file and creates .gcl objects

  - `base2gcl()` takes a `rubias` baseline object and creates .gcl
    objects

## Installation

You can install the package from GitHub using pak.

``` r

install.packages("pak")

pak::pak("commfish/GCLr")
```

## Reporting issues and requests

If you have any issues running the functions in this package, please
file an issue on [GitHub](https://github.com/commfish/GCLr/issues).

Issues can also be filed if you want to request enhancements to
functions or additional functions to be added to the package.
