all: doc install

install:
	Rscript -e "require(devtools);require(digest); install('.');"	

doc:
	Rscript -e "require(devtools);require(digest); document('./');"	

test:
	Rscript -e "library('testthat'); require(devtools); test('./');"

example1:
	Rscript -e "source('examples/example-slices.r')"
