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

Fixed rounding issue with `split_gtscore_loki_import()`where it wouldn't produce the last file if the number of lines remaining was less than half of the `nlines` argument.

# GCLr 0.2.0

## New additions

`split_gtscore_loki_import()` takes a large GTscore .csv Loki import file and splits it up into smaller import files.

# GCLr 0.1.0

## Comments

This is the first working version of the GCLr package. With this version, the GCLr repository was branched into main and development. Newer versions will have updates on new features, bug fixes, and improvements since this version.
