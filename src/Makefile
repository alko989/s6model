all: calcFmsy.so s6model.so

calcFmsy.so: calcFmsy.cpp
	Rscript --vanilla -e "TMB::compile('calcFmsy.cpp')"

s6model.so: s6model.cpp
	Rscript --vanilla -e "TMB::compile('s6model.cpp')"

clean:
	rm -rf *o
