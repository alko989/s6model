PACKAGE=s6model
R=R
VERSION=`grep "Version" s6model/DESCRIPTION | grep -oi '[0-9.]*'`
TARBALL=${PACKAGE}_${VERSION}.tar.gz

.PHONY: all doc-update install vignette test

all:
	make doc-update && make install && make test && make cover

doc-update:
	echo "roxygen2::roxygenize('$(PACKAGE)', roclets=c('rd', 'collate', 'namespace'))" | R --slave

build-package:
	R CMD build --no-build-vignettes --resave-data=no $(PACKAGE)

install:
	make build-package && echo 'remotes::install_deps("s6model")' | R --vanilla && R CMD INSTALL --preclean --no-multiarch $(TARBALL) && date

test: 
	echo "testthat::test_local('${PACKAGE}')" | R --slave

cover: 
	echo "covr::codecov(path = \"${PACKAGE}\", token = "d9ed0064-076f-4726-bbca-abef84e68339")" | R --slave	

vignette:
	{ \
	for f in ${PACKAGE}/vignettes/*.Rmd; do \
		R -e "library(rmarkdown); render('$$f', 'all')" ;\
	done ;\
	mkdir -p ${PACKAGE}/inst/doc ; \
	mv spict/vignettes/*.pdf spict/inst/doc/ ;\
	mv spict/vignettes/*.html spict/inst/doc/ ;\
	cp spict/vignettes/*.Rmd spict/inst/doc/ ;\	
	}