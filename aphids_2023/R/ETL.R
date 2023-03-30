#### SET UP ####

# Pull data from OSF repository: https://osf.io/j5vwy/
# osf2R is defunct, use osfr instead: https://cran.r-project.org/web/packages/osfr/vignettes/getting_started.html
library(tidyverse)
library(osfr)

osf_retrieve_node("https://osf.io/j5vwy/") %>%
  osf_ls_files() %>%
  osf_download("./aphids_2023/osf/",recurse=TRUE,conflicts="overwrite")

# Throws some weird UTF encoding error, just dl the zip
