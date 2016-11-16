CC=g++
CFLAGS=-g -Wall -std=c++11

all: build run plot

build:
	$(CC) $(CFLAGS) p3.cpp -o p3.o

run:
	./p3.o

plot:
	./plotgrid.m && xdg-open gridpoints.png

output: build
	./p3.o > p3_output.txt

clean:
	rm *.o outfile.bin gridpoints.png
