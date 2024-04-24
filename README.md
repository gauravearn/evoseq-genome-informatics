# evoseq

[![Build Status](https://travis-ci.org/sablokgaurav/evoseq.png?branch=master)](https://travis-ci.org/sablokgaurav/evoseq)  [![codecov](https://codecov.io/gh/sablokgaurav/evoseq/branch/master/graph/badge.svg)](https://codecov.io/gh/sablokgaurav/evoseq)

## Installing

<!-- If you're putting `evoseq` on CRAN, it can be installed with

    install.packages("evoseq") -->

The pre-release version of the package can be pulled from GitHub using the [devtools](https://github.com/hadley/devtools) package:

    # install.packages("devtools")
    devtools::install_github("sablokgaurav/evoseq", build_vignettes=TRUE)

## For developers

The repository includes a Makefile to facilitate some common tasks.

### Running tests

`$ make test`. Requires the [testthat](https://github.com/hadley/testthat) package. You can also specify a specific test file or files to run by adding a "file=" argument, like `$ make test file=logging`. `test_package` will do a regular-expression pattern match within the file names. See its documentation in the `testthat` package.

### Updating documentation

`$ make doc`. Requires the [roxygen2](https://github.com/klutometis/roxygen) package.

Gaurav Sablok \
Academic Staff Member \
Bioinformatics \
Institute for Biochemistry and Biology \
University of Potsdam \
Potsdam,Germany
