all: calcFmsy.so s6model.so

calcFmsy.so:
	Rscript --vanilla -e 'TMB::compile("calcFmsy.cpp")'

s6model.so:
	Rscript --vanilla -e 'TMB::compile("s6model.cpp")'

clean:
	rm -rf *o
