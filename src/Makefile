all: ellip kepler phot 
	gfortran -Ofast ellip.o phot.o kep.o photlib.o keplib.o -o wrapper.so

ellip: 
	gfortran -Ofast -unroll=12 -c ellip.f90

kepler: kep.o
	gcc -O3 -c -Wall -Werror -fpic keplib.c

phot: phot.o
	gcc -O3 -c -Wall -Werror -fpic photlib.c

kep.o:
	gfortran -Ofast -c -unroll=5 -fmax-stack-var-size=100000 kep.f90

phot.o: 
	gfortran -Ofast -c -unroll=5 -fmax-stack-var-size=100000 phot.f90

clean: 
	rm *.o
	rm *.so
	rm *.mod
