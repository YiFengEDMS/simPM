
<!-- README.md is generated from README.Rmd. Please edit that file -->

# simPM

# Welcome to the homepage of *simPM*\!

### `simPM`: SIMulation-based power analysis for Planned Missing designs.

`simPM` is an R package that automates the Monte Carlo simulation-based
search for optimal planned missing (PM) and ‘post hoc’ planned missing
(PHPM) designs. This R package is developed and maintained by [Yi
Feng](https://terpconnect.umd.edu/~yifeng94/) & [Dr. Gregory R.
Hancock](https://education.umd.edu/directory/gregory-r-hancock) from the
University of Maryland.

## The story behind *simPM*

> ***“Oh %&$\#\!, they cut my funding\!”***

We have heard that too often, unfortunately, from researchers who rely
on external funding to support their longitudinal studies. There are a
variety of reasons that may lead to the adjusted (shrunken) budget
announced by the funding agency. Further, what often adds to the
challenge is that the researchers are asked to provide a revised study
plan to convince that their project can be continued, showing
satisfactory inferential validity and statistical power, albeit with the
reduced budget.

> ***“What can I do except for firing my consultant and (poor) RAs?”***

Planned missing data methods present a very promising solution for such
challenging situations, providing many practical and methodological
advantages. However, careful planning is critical such that the planned
missing design can preserve enough information and statistical power.

`simPM` was created to help researchers survive the unexpected funding
cut in the course of a longitudinal study. It can be used to find an
**optimal** ‘post hoc’ planned missing design that allows the
researchers to complete the study at a reduced cost, while maintaining
satisfactory level of statistical power for testing the focal
parameters.

## What does *simPM* do?

By automizing the simulation-based power analysis for planned missing
designs in longitudinal context, `simPM` can free the researchers from
manually configuring the possible PM designs, determining their
eligibility, setting up the simulations, and summarizing the results
over replications, which can be tedious and time-consuming work
especially when there is a large number of plausible PHPM designs to be
evaluated.

## How to install *simPM*?

The source code of this R package is made public on the author’s [Github
page](https://github.com/YiFengEDMS/simPM). To install `simPM` on your
local machine, please run the following code

``` r
install.packages("devtools")
library(devtools)
intall_github("YiFengEDMS/simPM")
```

## Useful resources

Below are some tutorials demonstrating the usage of `simPM`:

1.  [Installation](articles/installation.html)
2.  [Package
    manual](https://yifengedms.github.io/simPM/manual/simPM_0.0.0.9000.pdf)
3.  [An example of PHPM with an autoregression and cross-lagged
    model](articles/Autoregressive-Cross-Lagged-Model.html)
4.  [An example of PHPM with a conditional linear
    LGM](articles/Conditional-Linear-Latent-Growth-Model.html)
    <!-- 1. [An example of PHPM with a second-order LGM](articles/Second-Order-Latent-Growth-Model.html) -->
5.  [Wave-level PM designs](articles/Wave-Level-PHPM.html)
6.  [Item-level PM designs](articles/Item-Level-PHPM.html)
7.  [Forward assembly](articles/Forward-Assembly-PHPM.html)
    <!-- 1. [Attrition] -->
8.  [Summarize the optimal design](articles/Summary.html)
9.  [Plot the optimal design](articles/PlotPM.html)

## Citation

Feng, Y. & Hancock, G. R. (in press). Oh %&$\#\!, they cut my funding:
Using ‘post hoc’ planned missing data designs to salvage longitudinal
research . *Child Development*.

  - A pre-peer reviewed version of the article is available
    [here](manuscript/Feng-and-Hancock-2020-prereview-main-text.pdf).

Feng, Y. & Hancock, G. R. (2019). simPM: SIMulation-based power analysis
for Planned Missing designs. R package version 0.0.0.9000.
<https://yifengedms.github.io/simPM/>

## Questions or Suggestions?

Send an email to [yifeng94@umd.edu](yifeng94@umd.edu). We are happy to
hear about your thoughts\!
