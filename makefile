CC       = g++
CXXFLAGS = -Wall -g -std=c++14 -I/opt/LCIO/include -L/opt/LCIO/lib -llcio -lsio

SRC  = mainsteredmonds.cpp mainsterkruskal.cpp mainsterprim.cpp example.cpp
EXEC = amstere amsterk amsterp example

all: $(EXEC)

amstere: mainsteredmonds.cpp
	$(CC) $(CXXFLAGS) -o amstere mainsteredmonds.cpp

amsterk: mainsterkruskal.cpp
	$(CC) $(CXXFLAGS) -o amsterk mainsterkruskal.cpp

amsterp: mainsterprim.cpp
	$(CC) $(CXXFLAGS) -o amsterp mainsterprim.cpp

example: example.cpp
	$(CC) $(CXXFLAGS) -o example example.cpp

clean:
	rm -rf figures/ graphtxt/ figuresexample/ *.root $(EXEC)
