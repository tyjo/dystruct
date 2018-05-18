CC=g++
CPPFLAGS=-std=c++11 -Wall -Wno-reorder -Ofast -g -fopenmp -isystem./lib/boost_1_62_0/

OBJS=src/main.o src/variational_kalman_smoother.o src/svi.o src/snp_data.o src/util.o

main : $(OBJS)
	$(CC) $(CPPFLAGS) -o bin/dystruct $(OBJS)

src/%.o : src/%.cpp
	$(CC) -c $(CPPFLAGS) $< -o $@

.PHONY : clean
clean:
	rm -f $(OBJS)
