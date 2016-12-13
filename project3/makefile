CPP=g++
CPPFLAGS=-g -Wall -std=c++11

all: show3 show4 show5

build:
	$(CPP) $(CPPFLAGS) p3.cpp curvExp.cpp curvStraight.cpp Domain.cpp Curvebase.cpp -o p3

run: build
	./p3

show%: task%
	xdg-open $^.jpg;

task%: build run
	echo $@ | ./plotgrid.m

output: build
	./p3 > p3_output.txt

clean:
	rm -f p3 *.bin task*.jpg
