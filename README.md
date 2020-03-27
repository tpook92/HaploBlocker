# HaploBlocker
R-package: Calculation of haplotype blocks and libraries

This is the github repository for your R-package HaploBlocker. 
A paper on HaploBlocker has been published in Genetics (https://www.genetics.org/content/212/4/1045)

The newest release 1.5.2 is mostly adding an optional modus to identify non-overlapping blocks. Outputs on defaults should be the same. The newest version of HaploBlocker is using RandomFieldsUtils 0.5.9 (This is more than what is available on CRAN!). 1.4.X and below is using RandomFieldsUtils 0.5.3.

For further questions (torsten.pook@uni-goettingen.de)

For explanation on the usage we refer to our user manual. The wiki is currently not equivilant (https://github.com/tpook92/HaploBlocker/wiki) - its just a couple of figure and not top priority as the user manual (Guidelines_to_HaploBlocker) should do.

## Version history

### Version 1.5.2

Added the option to identify non-overlapping haplotype blocks (overlap_remove = TRUE in block_calculation() )

Minor improvments to block_windowdataset()

Hotfixes in plot_block, block_calculation when providing base-pair info

### Version 1.4.8

HaploBlocker has been officially been published in Genetics (https://www.genetics.org/content/early/2019/05/31/genetics.119.302283)

A frozen version of the current package has been added in "Genetics_submission_version"

Min_majorblock is now automatically reduced in case the current parametrization would lead to a empty blocklist // haplotype library

### Version 1.4.7
Removed typo in plot_block

### Version 1.4.5
Support usage of vcf and pedmap files as input for genetic datasets

Creating of a window dataset (block_windowdataset). Block structure more similar to traditional haplotype block methods (set start/end points etc.)

Improvement to documentation (Guidelines_to_HaploBlocker)

### Version 1.4.4
Hotfix in block_dataset_construction

### Version 1.4.3
New features: selection signature detection (bEHH, iHH)

Improved computing time: (parallel computing, skipping unnecessary steps for large scale datasets)

Minor updates to RandomFieldsUtils: 0.4.0+ now required (available at CRAN)

### Version 1.0.0
Initial release

