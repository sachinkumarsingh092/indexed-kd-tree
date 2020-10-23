.DEFAULT_GOAL := compile

CC := gcc
CC := ${CC}
CFLAGS := -Wall -O0 -g -Wno-unused-function #-Wno-unused-variable -Wno-unused-but-set-variable
INCLUDES := -I/usr/local/include
LIBS := /usr/local/lib/libgnuastro.a\
	-lgit2 -ltiff -llzma -ljpeg -L/usr/local/lib -lwcs -lcfitsio -lz \
	/usr/lib/libgsl.so /usr/lib/libgslcblas.so -lm -lc -lpthread -pthread

PROGNAME := match
SHELL := /usr/bin/env bash
REGION :=


compile: ${PROGNAME}.c
	rm -f ./build/*.txt ./build/reg/* ./build/*.fits
	${CC} ${CFLAGS} $< ${INCLUDES} -o ${PROGNAME} ${LIBS} && ./${PROGNAME}

ds9:
	ds9 test-pv.fits -zscale -zoom to fit -regions ${REGION}.reg
