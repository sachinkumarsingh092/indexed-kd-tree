.DEFAULT_GOAL := compile

CC := gcc
CC := ${CC}
CFLAGS := -Wall -O0 -g -Wno-unused-function
INCLUDES := -I/usr/local/include
LIBS := /usr/local/lib/libgnuastro.a\
	-lgit2 -ltiff -llzma -ljpeg -L/usr/local/lib -lwcs -lcfitsio -lz \
	/usr/local/lib/libgsl.so /usr/local/lib/libgslcblas.so -lm -lc -lpthread -pthread

PROGNAME := match
SHELL := /usr/bin/env bash
REGION :=


compile: ${PROGNAME}.c
	rm -f ./build/*
	${CC} ${CFLAGS} $< ${INCLUDES} -o ${PROGNAME} ${LIBS} && ./${PROGNAME}

ds9:
	ds9 test-pv.fits -zscale -zoom to fit -regions ${REGION}.reg
