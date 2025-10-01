# GCLr 0.10.0

## Enhancements
`locus_stats` - added option to supply a hierfstat data object and added gene diversity (Hs) and allelic richness (Ar) to the output. The function now uses hierfstat::wc instead of hierfstat::varcop to calulate Fst and Fis and hierfstat::basic.stats to calculate Ho and Hs. This change simplified the function and it runs much faster than before.

`qc_template` - removed the code chunk that selects QC fish (not used) and added a code chunk to the end of the markdown to produce a QC bounce workbook with alternate species individuals included in the individual conflicts table if any are found; added individual heterozygosity to flag contaminated samples; added chunk to calculate QC power (how many project fish and QC fish have legit genotypes for QC analysis?); added code chunk to remove microsatellite loci from the project data if the project type is TaqMan or GT-seq and the species is chinook; added code chunk to produce genotyping rates plot to end of project analysis to help select samples for QC analysis.

## New additions
`plot_rubias_IA_zscores` - This function reads in a rubias individual posteriors csv file produced by run_rubias_mix() and produces a histogram of the individual assignment (IA) z-scores.

`calc_ind_het` - This function calculates individual heterozygosity; useful for detecting contaminated samples with GT-seq.

## Documentation
`IA_thresholds` - replaced loo_roc_rate_calc() (old function name) in the documentation with loo_rate_calc() (new function name)

`plot_loo_prec_rec` - replaced loo_roc_rate_calc() (old function name) in the documentation with loo_rate_calc() (new function name)

# GCLr 0.9.0

## Enhancements
`run_rubias_base_eval` - updated parallel loop so it can run on all cores regardless of the number of reporting groups being tested and added a progress bar.

## New additions
`get_collection_info` - this function pulls a collection information report from Loki

# GCLr 0.8.1

## Bug fixes
`custom_comb_rubias_output` - fixed issue where an error would occur when resummarizing from the collection_trace output.

# GCLr 0.8.0

## Enhancements
`custom_comb_rubias_output` - added ability to summarize '.fst' rubias output files and a progress bar that appears while building the repunit trace from files.

`summarize_rubias_base_eval` - previously this function called on
custom_comb_rubias_output() to get the eval estimates; however, evaluations with thousands of test mixtures were taking a long time (days) to summarize, and sometimes it would never finish. Now this is a stand alone function that is reduced down to only what is need for summarizing baseline evaluation mixtures. The parallel process for summarizing the evaluation output has been streamlined to work much faster. Now large numbers of mixtures can be summarize in minutes and the function shows a progress bar. Added ability to summarize '.fst' rubias output files.

`create_rubias_base` - added optional '.fst' output file

`create_rubias_mix` - added optional '.fst' output file

`plot_baseline_eval` - added optional 'group_names' argument to make sure groups are plotted in the correct order.

## Bug fixes
`base_eval_sample_sizes` - the sample sizes returned by this function are now rounded to the nearest whole number. 

`create_rubias_base_eval` - the example code has been corrected so it will work properly. When reading in mixture and baseline csv files, the col_types now default to character to avoid incorrect col_types errors in rubias.

# GCLr 0.7.0

## New additions

`inviv_pca` - reads in a GENEPOP file and performs an individual Principal Component Analysis.

`plot_indiv_pca` - takes the raw output from indiv_pca() and plots the coordinates for the first 3 principal components in an interactive 3-dimensional scatter plot. 

ex_genepop.gen - GENEPOP file for examples; this is the same as ex_genepop.txt, but has the .gen extension required by adegenet::read.genepop()

`qc_bounce` - for lab staff; produces an Excel file containing the project genotypes, the proportion of project fish with at least 80% of loci with scores, average success rate by plate, conflicts and success rates for individual samples, and conflicts and success rates by locus.

# GCLr 0.6.2

## Bug fixes

`get_tissue_location()` - added new freezers and removed old; specified location of output files; added timestamp to output.

# GCLr 0.6.1

## Bug fixes

`stratified_estimator_rubias()` - make `repunit` a factor when `bias_corr = TRUE` and rubias output created from files.

# GCLr 0.6.0

## Enhancements

`run_rubias_mix()` - This function is now able to run multiple MCMC chains in rubias.	To run multiple chains, set nchains = (pick a number). For single chain, nchain = 1 (default). The output from the function is saved as .csv files that have a new column "chain".	

`custom_comb_rubias_output()` - updated to summarize multiple chains. Outputs created previously can still be summarized (backwards compatible)

`stratified_estimator_rubias()` - updated to summarize multiple chains. Outputs created previously can still be summarized (backwards compatible)

## New additions

Three functions have been added to the package for testing reporting groups for individual assignment (IA) and determining assignment threshold cut-offs:

`loo_rate_calc()` - takes the leave-one-out output from `rubias::self_assign()` and calculates true positive, false negative, false positive, and true negative assignment error rates for each reporting group. The output from this function can then be used to create a precision-recall curve plot for determining an appropriate individual assignment threshold and whether a reporting group is sufficiently identifiable for producing individual assignment estimates.

`plot_loo_prec_rec()` - This function takes the unmodified output from `rubias::self_assign()`, calls on `loo_rate_calc()` to produce true positive and precision rates, plots the scaled likelihood for each baseline sample (y-axis) by reporting group or population in a box plot (optional), and plots recall (aka *true positive rate*) (x-axis) by precision (y-axis) in an interactive plot that can be used to assess reporting groups and assignment thresholds for individual assignment analyses.

`IA_thresholds()` - This function takes the unmodified output from `rubias::self_assign()`, calculates recall (aka true positive rate) and precision, and returns a tibble of minimum and maximum assignment thresholds for filtering top individual assignment results in order to meet the desired true positive rate and precision.

# GCLr 0.5.2

## Bug fixes

`plot_genepop_hwe()` - update to match tibble output from `read_genepop_hwe()`.

# GCLr 0.5.1

## Bug fixes

`qc_template.Rmd` - update so Concordance files are correctly found for both chip and GT-seq projects.

# GCLr 0.5.0

## Bug fixes

`qc_template.Rmd` - alternate species check was giving an eror if you did not have at least 1 alternate and 1 failure marker, now the function checks to make sure you have at least one of each.

## New additions

`qc_template()` - this function opens qc_template.Rmd and the documentation contains step by step instructions for running the QC script.

`get_fish_inventory()` - this function pulls a fish inventory from Loki and produces three summary plots: 1) a barplot of the number of samples in Loki by species and year, 2) a barplot of the overall number of fish in Loki by species for the years supplied, 3) a heatmap of the number of fish by month and year, with the overall average by month for the years supplied

# GCLr 0.4.4

## Bug fixes

`dupcheck_qc_conflicts()` - update to properly remove conflict fish that were lost in project or QC due to `remove_ind_miss_loci()`.

`qc_template()` - update directions to restart RStudio.

## Enhancements

`loo_roc_rate_calc()` - update to include *precision* - TP / (TP + FP) in the output.

# GCLr 0.4.3

## Bug fixes

`loki2r_gaps()` - when the function updated the `LocusControl$alleles` variable, it made it a list of vectors instead of a list of tibbles. This was causing functions that access the alleles variable to throw and error.

`gcl2nexus()` - the function gave an error if you only had one haplotype locus.

`qc_template()` - fixed QA table so it calculates `final` sample sizes when `wrong_spp` is `NA`.

`create_prior()` - update documentation to clarify `minval` and `sillyvec` arguments.

`dupcheck_within_silly()` - add a **Note** to the documentation clarifying that `$SILLY_CODE` needs to match the name of the `.gcl` object. This can cause issues if you modify `.gcl` objects outside of `pool_collections()`.

# GCLr 0.4.2

## Bug fixes

`loki2r()` and `loki2r_gaps()` - the regular expression used to correct the NA issue in version 0.4.1 had a typo in `gsub` such that values for ALLELE_1 overwrote both ALLELE_1 and ALLELE_2.

# GCLr 0.4.1

## Bug fixes

`loki2r()` and `loki2r_gaps()` - the regular expression used to correct the NA issue in version 0.4.0 did not work with na_if() in two of the loki2r functions. After switching to gsub, the issue is fixed.

# GCLr 0.4.0

## Bug fixes

`loki2r()`, `loki2r_gaps()`, and `loki2r_proj()` - the code in these functions to replace "0" allele calls with NA's did not take into account microsatellite alleles that contain zeros (e.g. "105", "109"), and any allele containing a zero became an NA. This issue has beeen corrected so only "0" alleles are converted to NA's.

`get_tissue_data()` - this function was removing workbench IDs that contain 0. This has been fixed.

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

2.  `dupcheck_among_sillys()` was missing some namespaces for functions, those have been added.

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
