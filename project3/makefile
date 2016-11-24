CPP=g++
CPPFLAGS=-g -Wall -std=c++11

all: plot

build:
	$(CPP) $(CPPFLAGS) p3.cpp -o p3

run: build
	./p3

plot: build run
	./plotgrid.m && xdg-open gridpoints.png

output: build
	./p3 > p3_output.txt

clean:
	rm -f p3 outfile.bin gridpoints.png
