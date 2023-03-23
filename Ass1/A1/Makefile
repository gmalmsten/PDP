###############################################################################
# Makefile for assignment 1, Parallel and Distributed Computing 2023.
###############################################################################

CC = mpicc
CFLAGS = -std=c99 -g -O3
LIBS = -lm

EXECUTABLE = stencil

all: $(EXECUTABLE)

stencil: stencil.c stencil.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
	
clean:
	$(RM) $(EXECUTABLE)
