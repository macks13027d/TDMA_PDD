# Compiler and complilation flag are included in Makefile.inc
include Makefile.inc

lib:
	cd src; make all

example:
	cd examples; make all

all:
	cd src; make all
	cd examples; make all

clean:
	cd src; make clean
	cd examples; make clean
