CC       = g++
CXXFLAGS = -Wall -g -std=c++14 -I/opt/LCIO/include -L/opt/LCIO/lib -llcio -lsio

SRC  = mainsterkruskal.cpp mainsterprim.cpp example.cpp
EXEC = amsterk amsterp example

all: $(EXEC)

amsterk: mainsterkruskal.cpp
	$(CC) $(CXXFLAGS) -o amsterk mainsterkruskal.cpp

amsterp: mainsterprim.cpp
	$(CC) $(CXXFLAGS) -o amsterp mainsterprim.cpp

example: example.cpp
	$(CC) $(CXXFLAGS) -o example example.cpp

clean:
	rm -rf figures/ graphtxt/ figuresexample/ $(EXEC)
