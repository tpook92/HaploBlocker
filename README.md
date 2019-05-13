# HaploBlocker
R-package: Calculation of haplotype blocks and libraries

This is the github repository for your R-package HaploBlocker. 
For a preprint we refer to biorvix (https://www.biorxiv.org/content/10.1101/339788v2). A new version of the preprint has been uploaded 18/3/2019.

We just uploaded a new version of the package (1.4.3) that includes function for selection signatures detection (bEHH, iHH), parallel computing for larger windows, and other smaller addon to improve computing time for large datasets as the 1000 Genomes Project.
Note that 1.0.0 is not running with the new version of RandomFieldsUtils (0.4.0) - it is highly recommended to use the newest version of both package - RandomFieldsUtils on CRAN is not the newest version!

For further questions (torsten.pook@uni-goettingen.de)

For explanation on the usage we refer to our user manual. The wiki is currently not equivilant (https://github.com/tpook92/HaploBlocker/wiki) - its just a couple of figure and not top priority as the user manual should do.

## Version history

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

