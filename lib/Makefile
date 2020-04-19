# Ad-hoc makefile to build static library
# Embedding source code here forks library, but removes dependency on
# external files, which is convenient for building
libtlhash.a: libtlhash.c
	${CC} -I../include -O2 -c libtlhash.c
	ar rcs libtlhash.a libtlhash.o
	-rm libtlhash.o
clean:
	-rm libtlhash.a
