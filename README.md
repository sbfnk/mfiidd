
# Model fitting and inference for infectious disease dynamics

This repository contains the material to create the mfiidd course page.

All the raw material is in the folder `Rmd/` and is written in
`Rmarkdown`. Any changes to the Rmarkdown files are automatically
updated on the web site once committed to the `main` branch.

## Local testing

The `html` pages can be generated locally using the function

``` r
rmarkdown::render_site("Rmd/")
```

This creates the `html` pages in the `Rmd/_site` directory.
