### R/CausalMST

[Brian Yandell](http://www.stat.wisc.edu/~yandell)

R/CausalMST is a support library for extending [R/qtl2](http://kbroman.org/qtl2), primarily associated with capabilities.reimplementation of [qtl2ggplot](https://github.com/byandell/qtl2ggplot) (for data visualization). See
[R/qtl2](http://kbroman.org/qtl2) (aka qtl2) for the bigger story of the qtl2 suite of routines.

---

### Installation

R/qtl2 is early in development and so is not yet available on
[CRAN](http://cran.r-project.org).

You can install R/qtl2 from [GitHub](https://github.com/rqtl).

You will need the following packages for DOread:

    install.packages(c("tidyverse", "devtools", "mnormt", "corpcor", "qtl"))

Once you have installed these, install CausalMST using `devtools::install_github()` as

    library(devtools)
    install_github("byandell/CausalMST")

Note that the R/qtl package is not required, but is used in the document `cmst.Rmd` to illustrate how QTL data can be used with these routines. A more modern use of R/qtl2 is presented in `mediate1_test.R`. For its use, see `shinyMediate1.R` in package `qtl2shiny`, where you can also learn about package dependencies. More detail is forthcoming. You will need at the least the following:

    install_github("rqtl/qtl2scan")
    install.packages("grid")

---

### Vignettes

coming ...

---

#### License

[Licensed](License.md) under [GPL-3](http://www.r-project.org/Licenses/GPL-3).
