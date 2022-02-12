all: photodynam

photodynam: deps
	gcc -O3 -Wno-write-strings -shared -o photodynam.so elliptic.o icirc.o scpolyint.o mttr.o source/photodynam.c

deps:
	gcc -c -O3 -Wno-write-strings source/elliptic.c
	gcc -c -O3 -Wno-write-strings source/icirc.c
	gcc -c -O3 -Wno-write-strings source/scpolyint.c
	gcc -c -O3 -Wno-write-strings source/mttr.c
