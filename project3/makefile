CPP=g++
CPPFLAGS=-g -Wall -std=c++11

all: plot

build:
	$(CPP) $(CPPFLAGS) p3.cpp curvExp.cpp curvStraight.cpp Domain.cpp Curvebase.cpp -o p3

run: build
	./p3

plot: build run
	echo 'task3' | ./plotgrid.m && xdg-open task3.pdf;
	echo 'task4' | ./plotgrid.m && xdg-open task4.pdf;
	echo 'task5' | ./plotgrid.m && xdg-open task5.pdf;

output: build
	./p3 > p3_output.txt

clean:
	rm -f p3 *.bin task{3..5}.pdf
