# Compiler and complilation flag are included in Makefile.inc
include Makefile.inc

all:
	cd src; mkdir -p obj; make all

clean:
	cd src; make clean
