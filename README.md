<!-- badges: start -->
[![R-CMD-check](https://github.com/adrienbrilhault/BRIL/workflows/R-CMD-check/badge.svg)](https://github.com/adrienbrilhault/BRIL/actions)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](pkg/LICENSE)
<!-- [![Analytics](https://ga-beacon.appspot.com/UA-200250973-1/Github-Readme?pixel)](https://github.com/igrigorik/ga-beacon) -->
<!-- badges: end -->
# BRIL

[BRIL](https://github.com/adrienbrilhault/BRIL) (Bootstrap and Refine Iterative Locator) is an [R](https://www.r-project.org) package for multivariate mode estimate.

The BRIL algorithm estimate the main mode through three principal steps:

- *BOOTSTRAP:* a recursive depth (or convex body minimizers) trimming to locate a first estimate.
- *REFINE:* a two-pass outliers filtering, the first relying on euclidean distances to the first estimate and unimodality tests, the second on robust distances and multinormality tests.
- *ITERATE:* after removing the samples selected in the REFINE step from the global distribution, the same process is reapplied, the first mode corresponding to the iteration with the largest number of samples.

To browse and test this package in an online RStudio environment, follow the link below: 

[![Binder](https://tinyurl.com/badgeRStudio)](https://mybinder.org/v2/gh/adrienbrilhault/BRIL/HEAD?urlpath=rstudio)


To cite this package, please refer to:

*Adrien Brilhault, Sergio Neuenschwander, and Ricardo Rios - A New Robust Multivariate Mode Estimator for Eye-tracking Calibration - [Behavior Research Methods, 2022](https://rdcu.be/cI9Pf)*
<br><br>


### Installation

To install BRIL from the [GitHub repository](https://github.com/adrienbrilhault/BRIL), add first the [remotes](https://github.com/r-lib/remotes) package if missing.

```r
install.packages("remotes")
```

Then proceed to the installation of the BRIL package and its dependencies
using the `install_github` function from 
[remotes](https://github.com/r-lib/remotes).

```r
remotes::install_github("adrienbrilhault/BRIL", subdir = "pkg")
```

<br>
*Note: As of May 2022, the
[OjaNP](https://cran.r-project.org/web/packages/OjaNP/) package - one of
of BRIC's dependencies - has been removed from the main CRAN repository.
If missing, install it manually with the following command, and then launch the BRIC
package installation again:*

``` r
remotes::install_github("cran/OjaNP", subdir = "pkg")
```
<br>

### Documentation

The package documentation is available in [BRIL-manual.pdf](https://github.com/adrienbrilhault/BRIL/raw/master/BRIL-manual.pdf), at [adrienbrilhault.github.io/BRIL/](https://adrienbrilhault.github.io/BRIL/), and directly in the R environment with the help commands.

```r
# main documentation page
help(package = "BRIL")

# bril() function help
?BRIL::bril()
```
<br>

### Example

If missing, install the [mvtnorm](https://CRAN.R-project.org/package=mvtnorm)
package from CRAN with the command `install.packages("mvtnorm")`. Then try the
following example, which creates a sample multivariate mixture, and estimates
its main mode.  

```r
library(BRIL)
library(mvtnorm)

XY <- rbind(
  rmvnorm(300, c(0,0), diag(2)*3-1),
  rmvnorm(100, c(15,20), diag(2)),
  rmvnorm(150, c(-10,15), diag(2)*2-0.5),
  rmvnorm(200, c(5,5), diag(2)*200)
)

res <- bril(XY)

print(res)
plot(res)
```
<br>

### Complementary code and data

This repository, beyond the R package itself (in the folder [pkg](https://github.com/adrienbrilhault/BRIL/blob/master/pkg)), provides the datasets et complementary materials related to the article *"Brilhault A, Neuenschwander S, Rios R - A New Robust Multivariate Mode Estimator for Eye-tracking Calibration - 2022"*

#### Datasets

The real-world datasets used in the article, consisting of eye-tracking calibrations, are provided in [data/eyeTrackingCalibrations.Rdata](https://github.com/adrienbrilhault/BRIL/raw/master/data/eyeTrackingCalibrations.Rdata), with an illustrative [video](https://www.youtube.com/watch?v=ZQzThht0VMw). They contain the `x` and `y` positions recorded from the eye-tracker in each trial, with their corresponding `metadata`, which includes the following attributes:

- `session`: subject name and session identifier,
- `trial`: trial number within the session,
- `nTrials`: total number of trials in the session,
- `duration`: duration of trials in the session (in ms),
- `target`: id of the calibration target displayed during the trial,
- `nTargets`: number of calibration targets in the session.

The table `groundTruth` provides the manually annotated reference coordinates for each of the calibration targets in every session. The labelling procedure, shown in this [video](https://www.youtube.com/watch?v=10ZapuMvK1s), is described in the article, and its code available in [data/groundTruthLabeling.R](https://github.com/adrienbrilhault/BRIL/blob/master/data/groundTruthLabeling.R).

The three datasets presented in the article as *Set1*, *Set2*, and *Set3* (with low, medium, and high contamination) correspond to the sessions *juj011a00*; *ded00800*, and *ded005a01*, respectively.

#### Additional code

The [benchmarks folder](https://github.com/adrienbrilhault/BRIL/blob/master/benchmarks) of the repository contains a series of additional scripts and materials used in our study, as for instance:

- The code for the benchmarks on artificial datasets, through a Monte-Carlo procedure ([BenchmarkBRIL_Artificial.R](https://github.com/adrienbrilhault/BRIL/blob/master/benchmarks/BenchmarkBRIL_Artificial.R));
- The code for the benchmark on real-world eye-tracking datasets ([BenchmarkBRIL_EyeTracking.R](https://github.com/adrienbrilhault/BRIL/blob/master/benchmarks/BenchmarkBRIL_EyeTracking.R)), and its analysis ([BenchmarkBRIL_EyeTracking_Analysis.R](https://github.com/adrienbrilhault/BRIL/blob/master/benchmarks/BenchmarkBRIL_EyeTracking_Analysis.R));
- The code for the detection of saccades and fixations in raw eye-tracking signals (files [FixationsDetection.R](https://github.com/adrienbrilhault/BRIL/blob/master/benchmarks/FixationsDetection.R) and [utils_fixations.R](https://github.com/adrienbrilhault/BRIL/blob/master/benchmarks/utils_fixations.R)), described in the technical report [FixationsDetection_Report.pdf](https://github.com/adrienbrilhault/BRIL/raw/master/benchmarks/FixationsDetection_Report.pdf);
- The code for the clustering algorithms used as baseline ([utils_clustering.R](https://github.com/adrienbrilhault/BRIL/blob/master/benchmarks/utils_clustering.R)).

#### Illustrative notebooks

- For a demonstration of our method with artificial data, please refer to the 
notebook [demo_ArtificialData.ipynb](demo_ArtificialData.ipynb) 
( **[?????? run binder](https://mybinder.org/v2/gh/adrienbrilhault/BRIL/HEAD?filepath=demo_ArtificialData.ipynb)** )
- For a demonstration of our method with eye-tracking data, please refer to the
notebook [demo_EyeTrackingData.ipynb](demo_EyeTrackingData.ipynb) 
( **[?????? run binder](https://mybinder.org/v2/gh/adrienbrilhault/BRIL/HEAD?filepath=demo_EyeTrackingData.ipynb)** )

*(note that the online viewer from Github sometimes needs one or two reload if the notebook was not rendered)*

To run these notebooks interactively in JupyterLab, follow this link:

[![Binder](https://tinyurl.com/badgeJupyterLab)](https://mybinder.org/v2/gh/adrienbrilhault/BRIL/HEAD?urlpath=lab)
