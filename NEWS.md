## Bug fixes

`loki2r()` and `loki2r_gaps()` - the regular expression used to correct the NA issue in version 0.4.1 had a typo in `gsub` such that values for ALLELE_1 overwrote both ALLELE_1 and ALLELE_2.

# GCLr 0.4.1

## Bug fixes

`loki2r()` and `loki2r_gaps()` - the regular expression used to correct the NA issue in version 0.4.0 did not work with na_if() in two of the loki2r functions. After switching to gsub, the issue is fixed.

# GCLr 0.4.0

## Bug fixes

`loki2r()`, `loki2r_gaps()`, and `loki2r_proj()` - the code in these functions to replace "0" allele calls with NA's did not take into account microsatellite alleles that contain zeros (e.g. "105", "109"), and any allele containing a zero became an NA. This issue has beeen corrected so only "0" alleles are converted to NA's.

`get_tissue_data()` - this function was removing workbench IDs that contain 0.  This has been fixed.

`skeleton.Rmd` - splitting the QC template into project and qc data analysis sections led to issues with removing fish prior to calculating conflict rates. After reading sillys from LOKI, each silly object is now duplicated using the name "SILLY_raw.gcl". The raw (i.e., unedited) SILLYS are used to calculate conflict rates in `combine_conflicts()`. 

`combine_conflicts()` - This now calls on the raw silly objects when running `pool_collections()`.

## Enhancements

Split qc_template into project data and QC data analysis sections. Added function to select high-quality samples for the QC.

## New additions

`select_qc_samps()` This function pulls project genotypes from Loki and removes samples with genotypes at less than 80% of loci, then selects QC samples by plate.

# GCLr 0.3.0

## Bug fixes

The output from `get_tissue_data()` previously contained zeros instead of NA's in the follow variables: IS_MISSING_PAIRED_DATA_EXISTS, WELL_HAS_MORE_THAN_ONE_SAMPLE, IS_PRESENT_IN_DATASHEET, IS_PRESENT_BUT_NOT_IN_DS This was messing up the Loki tissue importer because it imports any non NA values, including zeros, as being TRUE. In the case of IS_MISSING_PAIRED_DATA_EXISTS, importing a column of zero's and 1's wuold make Loki mark all tissues as missing. Now there are only 1's (missing tissues) and NA's (not missing) in the output.

`get_tissue_locations` - this function was not exported previously and the example had the wrong function names. Now the function documentation is correct and the function has been exported so it can be used in the package.

`get_gtseq_sample_sheet` - the output from this function had the columns in the wrong order for the GTseq pipeline, they are now in the correct order.

## New additions

`split_combchip_loki_import()` This function takes a combined Biomark import file and splits it up into multiple files so genotypes can be imported into Loki in smaller batches.

`loo_roc_rate_calc()` takes the leave-one-out output from `rubias::self_assign()` and calculates true positive, false negative, false positive, and true negative assignment error rates for each reporting group. The output from this function can then be used to create a Receiver Operator Characteristic curve plot for determining an appropriate individual assignment threshold and whether a reporting group is sufficiently identifiable for producing individual assignment estimates. A ROC curve plotting function will be added to a future version of GCLr.

`proj_geno_stats()` takes the genotypes report object from `get_geno()`, calculates genotypes statistics, and write them to a formatted Excel workbook.

## Enhancements
`get_geno()` has a new argument `proj.stats`. When `project_name` is supplied and `proj.stats = TRUE`, the function calls on `proj_geno_stats()`, which calculates genotypes statistics and writes them to a formatted Excel workbook. 

# GCLr 0.2.3

## Bug fixes

Issues with two functions used in the `qc_template` markdown have been fixed:

1.  The `read_qc_geno()` now adds 'qc' to the end of the silly code in the SILLY_CODE and SillySource variables of the '.gcl' objects it produces. Not having the 'qc' in the silly code was causing an error with `rubias::close_matching_samples()` because there were individuals with the same sample ID being checked(aka SillySource). 

2. `dupcheck_among_sillys()` was missing some namespaces for functions, those have been added.

# GCLr 0.2.2

## Bug fixes

Two issues with the `qc_template` markdown have been fixed:

1.  The `plot_sucessrate` code chunk would give an error if there were no success rates less than 0.8.

2.  The plot produced by the `qc_conflict_rate` code chunk did not have the correct number of bars when plotting less than 3 conflict rates, and when there was only one conflict rate to plot, nothing would show up on the plot.

## Enhancement

The `qc_conflict_rate` code chunk in the `qc_template` now prints a message below the chunk if no conflicts were found. Previously, a blank plot would appear.

# GCLr 0.2.1

## Bug fix

Fixed rounding issue with `split_gtscore_loki_import()` where it wouldn't produce the last file if the number of lines remaining was less than half of the `nlines` argument.

# GCLr 0.2.0

## New additions

`split_gtscore_loki_import()` takes a large GTscore .csv Loki import file and splits it up into smaller import files.

# GCLr 0.1.0

## Comments

This is the first working version of the GCLr package. With this version, the GCLr repository was branched into main and development. Newer versions will have updates on new features, bug fixes, and improvements since this version.
