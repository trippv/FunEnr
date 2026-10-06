usethis::use_git()
usethis::use_github()
install()
usethis::install()
usethis::use_mit_license("Miguel Tripp")
usethis::use_r("Function_TopgoEnrich.R")
check()

devtools::check()
devtools::install()
devtools::document()
usethis::use_readme_md()
usethis::use_vignette("FunEnr", "Functional Enrichment tools")
devtools::install(build_vignettes = TRUE)


# generar funcion separator
usethis::use_r("detect_separator.R")




#add dependencies
usethis::use_package("stringr")
usethis::use_package("GOSemSim")
usethis::use_package("org.Hs.eg.db")
usethis::use_package("topGO")
usethis::use_package("rrvgo")

devtools::check()


#### Atualizacion del paquete del 06/10/2026
### cambios:
##### Actualizacion de los encabezados de ambas funciones
#### Actualización del README
#### Actualziacion de DECRIPTION
devtools::document()

devtools::check()
devtools::build_vignettes()

## verificar el remoto
gert::git_remote_list()
gert::git_status()
devtools::document()
devtools::check()
