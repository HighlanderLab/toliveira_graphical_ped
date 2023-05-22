# Pedigree based mixed-effects models using directed acyclic graphs

This repository presents an evaluation of different MCMC samplers and parameterizations for the animal model, a statistical model commonly used in animal breeding. The focus is on the efficiency and accuracy of the different implementations, with comparisons made between JAGS and NIMBLE packages as well as various parameterizations of the model. The results show that the transformed additive genetic values approach consistently produces chains with a considerably larger effective sample size due to prior independence among additive genetic values. However, the Cholesky decomposition of the NRM using the Gibbs sampling algorithm is the most time-consuming approach, and the implementation with the transformed additive genetic values constantly produced chains with a considerably lower effective sample size when using Gibbs sampling sampler. Furthermore, NIMBLE was found to be more efficient than JAGS for running complex models. These findings have important implications for the choice of MCMC samplers and parameterizations when working with animal models in animal breeding.

# Installing and Specifying Package Versions for R Script

This document provides instructions for installing and specifying the version of the packages used in the `R` script. The `R` script contains the following packages:

* `pacman`
* `MasterBayes`
* `tidyverse`
* `knitr`
* `runjags`
* `MCMCvis`
* `MCMCglmm`
* `gdata`
* `coda`
* `ggrepel`
* `dplyr`
* `kableExtra`
* `MatrixModels`
* `pedigreemm`
* `prettycode`
* `formattable`
* `AlphaSimR`
* `patchwork`
* `animalModels`

To install the required packages and ensure their proper version, follow these steps:

1. Download and install "Jags" from https://sourceforge.net/projects/mcmc-jags/

2. First, check if the "pacman" package is installed. If not, install it using the following command:

```
if (!require("pacman")) {
  install.packages("pacman")
}
```

3. Next, check if the `MasterBayes` package is installed. If not, install it with a specific version (2.58) using the following command:

If you are executing the code on macOS, please follow these steps to address the possible missing dependencies: i) Ensure that you have the Xcode Command Line Tools installed on your system. You can install them by running the following command in the Terminal:

```
xcode-select --install
```

Additionally, you may need to install the GNU Fortran compiler. This is required for certain R packages. You can download the GNU Fortran compiler from the official GNU Fortran website (https://gcc.gnu.org/wiki/GFortranBinaries) and follow the installation instructions specific to your macOS version. First you need to install `brew` (https://brew.sh/) and then use the command below in the terminal:

```
brew install gcc
```

For Linux OS:

If you are executing the code on a Linux system, please follow this step to address any missing dependencies:

```
sudo apt-get install gfortran
```

For Windows OS:

Make sure you have Rtools installed on your system. Rtools provides the necessary tools and compilers for building R packages. You can download Rtools from the official website (https://cran.r-project.org/bin/windows/Rtools/).


After following the appropriate steps based on your operating system, you can proceed to install `MasterBayes` by executing the following code:

```
if (!require("MasterBayes")) {
  # Install devtools package if not already installed
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  # Install remotes package from GitHub source
  devtools::install_github("r-lib/remotes")
  install.packages(c("quadprog", "mvtnorm", "gdata", "kinship2", "genetics"))
  pacman::p_load(
  quadprog, version = "1.5.8",
  mvtnorm, version = "1.1.3",
  gdata, version = "2.19.0",
  kinship2, version = "1.9.6",
  genetics, version = "1.3.8.1.3"
  )
  remotes::install_version("MasterBayes", version = "2.58", repos = "http://cran.us.r-project.org")
}
```

4. Load the required packages using the `pacman` package with the following command:

```
pacman::p_load(
  tidyverse, version = "1.3.2",
  knitr, version = "1.42",
  runjags, version = "2.2.1-7",
  MCMCvis, version = "0.15.5",
  MCMCglmm, version = "2.34",
  gdata, version = "2.18.0.1",
  MasterBayes, version = "2.58",
  coda, version = "0.19-4",
  ggrepel, version = "0.9.2",
  dplyr, version = "1.0.10",
  kableExtra, version = "1.3.4",
  MatrixModels, version = "0.5-1",
  pedigreemm, version = "0.3-3",
  prettycode, version = "1.1.0",
  formattable, version = "0.2.1",
  AlphaSimR, version = "1.3.4",
  patchwork, version = "1.1.2"
)
```

5. Check if the `animalModels` package is installed. If not, install it from a local source file using the following command:

```
if (!require("animalModels")) {
  install.packages("./R/animalModels_0.0.1.tar.gz", repos = NULL,
                   type="source")
}
```

6. Load the `animalModels` package and use the "prettycode" function with the following command:

```
library(animalModels)
prettycode::prettycode()
```

The R script is compatible with the following `R` environment: `R` version: 4.1.3 (2022-03-10)

Please ensure that your `R` environment matches the above specifications before running the script.


In case you encounter issues using pacman for package installation or prefer to use the remotes package, follow this guide to install the specific versions of the packages using the remotes::install_version() function:

1. Install the `remotes` package if it's not already installed:

```
if (!require("remotes")) {
  install.packages("remotes")
}
```

2. Install the specific versions of each package using the `remotes::install_version()` function:

```
# Replace the "package_name" and "package_version" with the desired package and version
if (!require("package_name")) {
  remotes::install_version("package_name", version = "package_version", repos = "http://cran.us.r-project.org")
}
```

3. Replace "package_name" and "package_version" with the names and versions of the required packages:

```
# Example: Installing tidyverse version 1.3.2
if (!require("tidyverse")) {
  remotes::install_version("tidyverse", version = "1.3.2", repos = "http://cran.us.r-project.org")
}
```

4. Repeat step 3 for all the required packages.

Please note that installing specific versions of packages may lead to compatibility issues if a package depends on a different version of another package. Make sure to verify the compatibility between the packages and their versions to avoid any issues.

