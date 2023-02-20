MAKEFILE=Makefile

.PHONY: help

#------------------------------------------------------------------
# Targets
#------------------------------------------------------------------


help:
	@echo ""
	@echo "- Available targets:"
	@echo ""
	@perl -ne 'if(	/^(\w+):/){print "\t",$$1,"\n"}' $(MAKEFILE)
	@echo ""
	@echo ""

clean:
	@rm -f src/*.o src/*.so; rm -f scomicstools.Rcheck/dbfmcl/libs/dbfmcl.so; rm -rf ./dbfmcl.Rcheck
	@rm -rf /tmp/dbfmcl; rm -rf *dbf_out.txt; rm -rf *mcl_out.txt 

check: clean
	@rm -rf /tmp/scomicstools; mkdir -p /tmp/scomicstools; cp -r ./* /tmp/scomicstools
	@R CMD check /tmp/scomicstools

doc:
	@echo ">>> Creating a package documentation"
	@echo "library(roxygen2); roxygen2::roxygenise()" | R --slave

install:
	@echo ">>> Installing..."
	@rm -f src/*.o src/*.so
	@R CMD INSTALL .

test:
	@echo ">>> Testing package"
	@rm -rf `ls tests/testthat/| grep -v \R$$`
	@echo "devtools::test()" | R --slave

all: doc install check test
