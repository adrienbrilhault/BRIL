# BRIL 1.0.4

* Fixed bug in depth_values() when parameter `u` is provided as dataframe
* Update article reference (Behaviour Research Method, 2022)
* Add install command for OjaNP (after being removed from Cran)

# BRIL 1.0.3

* Fixed bug with MVE-based recursive median which could loop indefinitely
* Fixed bug with Spatial depth median from package "depth"
* Fixed bug with Sample Median and Multivariate Median being inverted in "median_mv()"
* Added references to documentation
* Added benchmarks code, and data labeling procedure to the repository

# BRIL 1.0.2

* Removed Search from Git Page (pkgdown feature not working)
* Removed title metadata from jupyter notebook (causing warning)
* Small fix to R-CMD-Check workflow to specify the path
* Updates to package documentation
* Fix docker error with libstdc++.so.6 library in mybinder RStudio env

# BRIL 1.0.1

* Added R-CMD-Check to Github workflow
* Updates of the functions documentation
* Added examples to plot.BRIL()
* Fixes to code formatting
* Orthographic corrections
* Bug fix in median_mv() examples
* Update of pkgdown website
    * use of Bootstrap 4 and Cosmo Template
    * jupyter notebooks knitted into Articles
    * added link to Keiji Matsuda website
    * functions reference grouped in sections
    * switched figures format to SVG
* Added a `NEWS.md` file to track changes to the package
