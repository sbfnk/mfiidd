# Model fitting and inference for infectious disease dynamics

This repository contains all the material to create the mfiidd webpage with the practicals.

All the raw material is in the folder `/Rmd` and is written in `Rmarkdown`.

The `html` pages are generated using the function `generate_html` in `R/generate_html.R`.

The `html` pages are created in `/website`, whose content need to be transfered to the `gh-pages` branch.

There are two additional folder in `website` containing the slides and external references, those folders need to be transefered too as the html pages point to them.

## TODO

* To have all the `.tex` for the presentations in a `tex` folder. At the moment there is only Helen's presentation.
* To have the `Rintro` as a `html` file as well as the `Rmd` source.
* Generate the pdf for all the practicals
* I could not manage to recompile the tex file of Helen using Seb's template, looks like the issue is with the `usefonttheme{seb}` command. Not sure whether I have to install the font manually..

