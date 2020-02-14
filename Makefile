install:
	Rscript -e 'devtools::document()'
	Rscript -e 'devtools::install(upgrade = FALSE)'
pkgdown:
	Rscript -e 'pkgdown::build_site()'
all: install pkgdown
