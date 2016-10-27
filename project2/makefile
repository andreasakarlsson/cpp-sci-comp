CC=g++
CFLAGS=-g -Wall -std=c++11

p2t1:
	$(CC) $(CFLAGS) p2t1.cpp -o p2t1.o && ./p2t1.o

p2t1output:
	./p2t1.o > p2t1_output.txt

p2t2: r8lib.h r8mat_expm1.h matrix.h
	$(CC) $(CFLAGS) p2t2.cpp -o p2t2.o matrix.cpp r8mat_expm1.cpp r8lib.cpp && ./p2t2.o

p2t2output:
	./p2t2.o > p2t2_output.txt

clean:
	rm *.o
