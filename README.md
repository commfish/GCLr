
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

  - Individual assignment (IA) evaluation

    - `loo_rate_calc()` calculates leave-one-out (LOO) error rates for
      each reporting group

    - `plot_loo_prec_rec()` plots LOO results as a precision-recall
      curve

    - `IA_thresholds()` calculates IA probablity thresholds based on
      precision-recall standards

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

## Package maintenance

This package generally follows the
[git-flow](https://nvie.com/posts/a-successful-git-branching-model/)
branching model, using semantic versioning to document releases. Below
is a quick summary of the different branches.

- `main`

  - stable version of the package, commits/merges to `main` trigger a
    new version number

- `develop`

  - ongoing, general improvements to the package including minor bug
    fixes

- `feature-branches`

  - specific improvements to the package (i.e. creating functions for
    individual assignment)

- `hotfix`

  - for fixing a serious bug found on `main` in the latest version of
    the package

  - references an [issue](https://github.com/commfish/GCLr/issues)

  - merge back into `main` **and** `develop`, triggers a new version
    number

### Update version of GCLr package

Below is a generalized protocol for updating the package version by
merging changes from the `develop` branch into the `main` branch. If
working from a `feature-branch`, follow this same workflow for
`feature-branch` –\> `develop`, then `develop` –\> `main`.

1.  commit changes on `develop` branch

2.  update `NEWS.md`, but without the version header \#1, commit

3.  update the package version with `usethis::use_version()`, choose
    major, minor, or patch, have it commit for you

4.  push to `develop`

5.  create a [pull request](https://github.com/commfish/GCLr/pulls)
    (base: `main` compare: `develop`) with meaningful title (i.e.,
    merging to version 1.X.X) and brief description

6.  pull requests to `main` require review so we can keep `main` stable!

7.  merge after folks approve, confirm the merge, do not delete
    `develop` branch

8.  e-mail all GCL staff to notify everyone about the update

### Hotfix protocol

1.  pull from `main` to make sure you have the latest and greatest
    version

2.  if your serious bug still exists, create an
    [issue](https://github.com/commfish/GCLr/issues)

3.  create `hotfix_issue_XX` branch (referencing the
    [issue](https://github.com/commfish/GCLr/issues) \#) from `main`

4.  working on the `hotfix_issue_XX` branch, commit necessary changes to
    resolve the issue, note that you can include
    [keywords](https://docs.github.com/en/issues/tracking-your-work-with-issues/linking-a-pull-request-to-an-issue)
    in your commit that will auto-magically close the issue once the
    `hotfix_issue_XX` branch is merged back into `main`

5.  update `NEWS.md`, but without the version header \#1, commit

6.  update the package version with `usethis::use_version()`, choose
    patch, have it commit for you

7.  push to `hotfix_issue_XX`

8.  create a [pull request](https://github.com/commfish/GCLr/pulls)
    (base: `main` compare: `hotfix_issue_XX`) with meaningful title
    (i.e., Hotfix issue \#) and brief description

9.  pull requests to `main` require review so we can keep `main` stable!

10. merge after folks approve, confirm the merge, do not delete
    `hotfix_issue_XX` branch yet

11. create another [pull request](#0) (base: `develop` compare:
    `hotfix_issue_XX`) with meaningful title (i.e., Hotfix issue \#
    merging to develop) and brief description

12. merge, confirm the merge, now you can delete the `hotfix_issue_XX`
    branch

13. e-mail all GCL staff to notify everyone about the update

## Other helpful resources

ADF&G Division of Sport Fisheries [Introduction to
Git](https://adfg-dsf.github.io/Git_book/) provides a good overview of
Git, however, note that generally assumes that you will be working off
of the shared network drive, rather than cloning to your local C:/
drive.

ADF&G’s Reproducible Research [R Best
Practices](https://adfg-dsf.github.io/Best_practice_R/Best_practice_R.html)

[GitKraken](https://www.gitkraken.com/) is a nice GUI alternative to
GitHub for visualizing the commit/branch network
