CPP=g++
CPPFLAGS=-g -Wall -std=c++11

all: plot

build:
	$(CPP) $(CPPFLAGS) p3.cpp curvExp.cpp curvStraight.cpp Domain.cpp Curvebase.cpp -o p3

run: build
	./p3

plot: build run
	echo 'task3' | ./plotgrid.m && xdg-open task3.jpg;
	echo 'task4' | ./plotgrid.m && xdg-open task4.jpg;
	echo 'task5' | ./plotgrid.m && xdg-open task5.jpg;

output: build
	./p3 > p3_output.txt

clean:
	rm -f p3 *.bin task3.jpg task4.jpg task5.jpg
