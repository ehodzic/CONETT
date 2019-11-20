GUROBIROOT=/path_to_gurobi_root_folder

CFLAGS	= -m64 -O4 -std=c++11
INC		= ${GUROBIROOT}/linux64/include/
CPPLIB   = -L${GUROBIROOT}/linux64/lib/ -lgurobi_c++ -lgurobi81

all: conett

objects:
	g++ -c $(CFLAGS) bitmask.cpp -o bitmask.o
	g++ -c $(CFLAGS) subnetwork.cpp -o subnetwork.o
	g++ -c $(CFLAGS) graph.cpp -o graph.o
	g++ -c $(CFLAGS) entry.cpp -o entry.o

conett: objects
	g++ $(CFLAGS) -o $@ conett.cpp bitmask.o subnetwork.o graph.o entry.o heap.h -I$(INC) $(CPPLIB) -lm

.cpp.o:
	g++ -c $(CFLAGS) $^ -o $@ -I$(INC) $(CPPLIB) -lm

clean:
	rm -rf *.o conett
