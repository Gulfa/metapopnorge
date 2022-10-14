drat:::add("ncov-ic")
install.packages(c("abind", "tideyverse", "data.table", "odin", "odin.dust", "dust", "mcstate"))
options(
  repos = structure(c(
    SPLVERSE  = "https://docs.sykdomspulsen.no/drat/",
    CRAN      = "https://cran.rstudio.com"
  ))
)

install.packages("spldata")

devtools::install_github(repo="https://github.com/Gulfa/metapop")
devtools::install_github(repo="https://github.com/Gulfa/metapopnorge")
