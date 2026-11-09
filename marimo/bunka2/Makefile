all:
	RFILES = debug.o gaus.o inv_mat.o chash-p.o sha3.o ymo.o

mcElies: $(RFILES)
	cc -o $@ $(RFILES)

clang:
	clang -O3 vmo.c
