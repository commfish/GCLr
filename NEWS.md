# GCLr 0.2.1

## Bug fix
Fixed rounding issue with `split_gtscore_loki_import()`where it wouldn't produce the last file if the number of lines remaining was less than half of the `nlines` argument.

# GCLr 0.2.0

## New additions
`split_gtscore_loki_import()` takes a large GTscore .csv Loki import file and splits it up into smaller import files.

# GCLr 0.1.0

## Comments
This is the first working version of the GCLr package. With this version, the GCLr repository was branched into main and development. Newer versions will have updates on new features, bug fixes, and improvements since this version.