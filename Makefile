install:
	Rscript -e 'devtools::document()'
	Rscript -e 'devtools::install(upgrade_dependencies = FALSE)'
pkgdown:
	Rscript -e 'pkgdown::build_site()'
all: install pkgdown
