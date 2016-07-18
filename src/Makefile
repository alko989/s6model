all: calcFmsy.so s6model.so s6modelts.so

calcFmsy.so: calcFmsy.cpp
	Rscript --vanilla -e "TMB::compile('calcFmsy.cpp')"

s6model.so: s6model.cpp
	Rscript --vanilla -e "TMB::compile('s6model.cpp')"

s6modelts.so: s6modelts.cpp
	Rscript --vanilla -e "TMB::compile('s6modelts.cpp')"

clean:
	rm -rf *o
