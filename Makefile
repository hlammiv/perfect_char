LDLIBS = -lm
#CFLAGS = -O0 -g -pg -std=gnu99
#CPPFLAGS = -O0 -g -pg 
#CXX=g++-8
CXXFLAGS = -O3 -g -pg -fopenmp -std=c++11
CFLAGS = -O3 -std=gnu99 -Wall -fopenmp 
#CFLAGS = -O3 -std=gnu99 -Wall -fsanitize=address
all: dym-mod-metro

dym-mod-metro: dym-mod-metro.cpp lattice.o group.o timer.o

clean:
	${RM} dym-mod-metro lattice.o group.o timer.o gmon.out
