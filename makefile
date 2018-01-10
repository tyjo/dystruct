CC=g++
CPPFLAGS=-std=c++0x -Wall -Wno-reorder -Ofast -g -fopenmp -isystem./lib/boost_1_62_0/

OBJS=src/main.o src/variational_kalman_smoother.o src/cavi.o src/snp_data.o

main : $(OBJS)
	$(CC) $(CPPFLAGS) -o bin/dystruct $(OBJS)

src/%.o : src/%.cpp
	$(CC) -c $(CPPFLAGS) $< -o $@

.PHONY : clean
clean:
	rm -f $(OBJS)
