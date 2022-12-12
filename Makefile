help:
	@echo "Target available:"
	@cat Makefile | perl -ne 'if(/^\w+:$$/){s/://; s/^/\t- /; print}'

test:
	@Rscript -e 'library(devtools); devtools::test()'

doc:
	@Rscript -e 'library(roxygen2); roxygen2::roxygenise()'
install:
	@cd ..; R CMD INSTALL ohmiki

