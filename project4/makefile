CPP=g++
CPPFLAGS=-g -Wall -std=c++11


ifeq ($(OS),Windows_NT)
    RM=del -f
else
    RM=rm -f
endif

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
	$(RM) p4.exe p4 *.bin
