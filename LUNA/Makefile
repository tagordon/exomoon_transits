all: luna.o
	gcc -shared -o wrapper.so luna.o

luna.c:
	gcc -c -Wall -Werror -fpic -o luna.o luna.c
