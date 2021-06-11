# BRIL


[BRIL](https://github.com/adrienbril/BRIL) (Bootstrap and Refine Iterative Locator) is an [R](https://www.r-project.org) 
package for robust mode estimate in multivariate distributions. To cite this 
package, please refer to:

*Adrien Brilhault, Sergio Neuenschwander, and Ricardo Rios - A New Robust Multivariate Mode Estimator for Eye-tracking Calibration - Behavior Research Methods, 2021 (submitted)*
<br><br>

### Installation

To install BRIL from this [GitHub repository](https://github.com/adrienbrilhault/BRIL) 
, add first the [remotes](https://github.com/r-lib/remotes) package if absent.

```r
install.packages("remotes")
```

Then proceed to the installation of the BRIL package with its dependencies
using the `install_github` function from 
[remotes](https://github.com/r-lib/remotes).

```r
remotes::install_github("adrienbril/BRIL", subdir = "pkg")
```
<br>

### Documentation

After loading the package, refer to its documentation directly in the R 
environment, or, alternatively, consulting the 
[pdf manual](/pkg/BRIL-manual.pdf). 

```r
library(BRIL)

# main documantation page
help(package = "BRIL")

# bril() function help
?bril()
```
<br>

### Example use

If absent, install the [mvtnorm](https://CRAN.R-project.org/package=mvtnorm) 
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

This repository, beyond the R package itself (in the folder [pkg](/pkg/)), provides 
additional datasets et complementary materials related to the article *"Brilhault A, Neuenschwander S, Rios R - A New Robust Multivariate Mode Estimator for Eye-tracking Calibration - 2021"*

#### Datasets

The real-world datasets used in the article, consisting of eye-tracking 
calibrations, are provided in [data/eyeTrackingCalibrations.Rdata](data/eyeTrackingCalibrations.Rdata).
They contain the `x` and `y` positions recorded from the eyetracker in each trial,
with their corresponding `metadata`, which includes the following attributes:

- `session`: subject name and session identifier
- `trial`: trial number within the session
- `nTrials`: total number of trials in the session
- `duration`: duration of trials in the session (in ms)
- `target`: calibration target identifier displayed during the trial
- `nTargets`: number of calibration targets in the session
 
The table `goundTruth` finally provide the manually annotated reference 
coordinates for each of the calibration targets in every session.

The three datasets presented in the article as *Set1*, *Set2*, and *Set3* (with 
low, medium, and high contamination) correspond to the sessions *juj011a00*, 
*ded00800*, and *ded005a01*, respectivly.

#### Illustrative notebooks

 - For a demontration of BRIL algorithm on artificial data, please refer to the 
notebook [demo_ArtificialData.ipynb](demo_EyeTrackingData.ipynb)
 - For a demontration of BRIL algorithm on eye-tracking data, please refer to the
notebook [demo_EyeTrackingData.ipynb](demo_EyeTrackingData.ipynb)

To run these notebook interactively online, follow this link: 

*(note that the online jupyter viewer from Github sometimes need one or two reload if the notebook was not rendered)*