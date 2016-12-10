CPP=g++
CPPFLAGS=-g -Wall -std=c++11

all: run

build:
	$(CPP) $(CPPFLAGS) \
	p4.cpp \
	GFkt.cpp \
	MatrixNEW.cpp \
	./../project3/Domain.cpp \
	./../project3/Curvebase.cpp \
	./../project3/curvExp.cpp \
	./../project3/curvStraight.cpp \
	-o p4

run: build
	./p4

clean:
		del -f p4.exe *.bin