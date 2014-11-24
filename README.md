# Model fitting and inference for infectious disease dynamics

This repository contains all the material to create the mfiidd webpage with the practicals.

All the raw material is in the folder `/Rmd` and is written in `Rmarkdown`.

The `html` pages are generated using the function `generate_html` in `R/generate_html.R`.

The idea is then to transfer all the `html` pages to the `gh-pages` branch.

Since the `html` pages link to different folder conatining the slides and external references, thos folder need to be transefered too.
Here is a list of the folder to transfer:

* slides
* external_ref

## TODO

* To have all the `.tex` for the presentations in a `tex` folder. At the moment there is only Helen's presentation.
* To have the `Rintro` as a `html` file as well as the `Rmd` source.
* Generate the pdf for all the practicals

