library(usethis)
library(devtools)
library(pkgdown)

## This script is for setting up development dependencies and should be run once during initial setup or when updating dependencies.
#create_package("/omics/groups/OE0219/internal/MJMC/Tools/quothic")

## Optionally set up GitHub PAT for installations from GitHub
# usethis::edit_r_environ()

## Install necessary packages from GitHub if not available on CRAN or specific versions are needed
# devtools::install_github("jlmelville/uwot@v0.1.16")
# devtools::install_github("satijalab/sctransform@v0.3.5")
# devtools::install_github("satijalab/seurat@v4.3.0")
# devtools::install_github('satijalab/seurat-data')
#if (!require(moments)) install.packages("moments")



## Add dependencies to DESCRIPTION
use_package("Seurat", type = "Imports")
use_package("dplyr", type = "Imports")
use_package("tidyr", type = "Imports")
use_package("ggplot2", type = "Imports")
use_package("mclust", type = "Imports")
use_package("moments", type = "Imports")
use_package("SeuratData", type="Imports")
use_package("cluster", type="Imports")
use_package("viridis", type="Imports")

library(devtools)
devtools::document()
devtools::load_all()

build_site()
build_site_github_pages()

remove.packages("quothic")
# Restart R session
setwd("/omics/groups/OE0219/internal/MJMC/Tools/quothic")
devtools::install(upgrade = "never")



library(usethis)
use_testthat()
use_test("shannon_entropy")

use_vignette("introduction_to_quothic")

