### R/CausalMST

[Brian Yandell](http://www.stat.wisc.edu/~yandell)

R/DOread is a support library for extending [R/qtl2](http://kbroman.org/qtl2), primarily associated with capabilities.reimplementation of [qtl2ggplot](https://github.com/byandell/qtl2ggplot) (for data visualization). See
[R/qtl2](http://kbroman.org/qtl2) (aka qtl2) for the bigger story of the qtl2 suite of routines.

---

### Installation

R/qtl2 is early in development and so is not yet available on
[CRAN](http://cran.r-project.org).

You can install R/qtl2 from [GitHub](https://github.com/rqtl).

You will need the following packages for DOread:

    install.packages(c("broman", "dplyr"))

Then, install R/qtl2 using `devtools::install_github()`.

    library(devtools)
    install_github(paste0("rqtl/qtl2", c("scan")))

Once you have installed these, install qtl2ggplot as

    install_github("byandell/DOread")

---

### Vignettes

coming ...

---

#### License

[Licensed](License.md) under [GPL-3](http://www.r-project.org/Licenses/GPL-3).
