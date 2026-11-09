all:
	gcc -O3 -mavx2 vmo.c

niederreiter_demo: niederreiter_demo.c vmo.c patterson_q.c
	gcc -O3 -mavx2 -fopenmp -o $@ niederreiter_demo.c
	
mcElies: $(RFILES)
	cc -o $@ $(RFILES)

clang:
	clang -O3 vmo.c
