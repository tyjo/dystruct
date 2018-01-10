CC=g++
CPPFLAGS=-std=c++0x -Wall -Wno-reorder -Ofast -g -fopenmp -isystem./boost_1_62_0/
LFLAGS=-L./boost_1_62_0/lib/ -lboost_program_options

OBJS=src/main.o src/variational_kalman_smoother.o src/cavi.o src/snp_data.o

main : $(OBJS)
	$(CC) $(CPPFLAGS) $(LFLAGS) -o bin/dystruct $(OBJS)

src/%.o : src/%.cpp
	$(CC) -c $(CPPFLAGS) $< -o $@

.PHONY : clean
clean:
	rm -f $(OBJS)
