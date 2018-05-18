CC=g++-7
CPPFLAGS=-std=c++11 -Wall -Wno-reorder -Ofast -g -fopenmp -isystem../lib/boost_1_62_0/

OBJS=test_input.o ../src/snp_data.o ../src/util.o

test_input : $(OBJS)
	$(CC) $(CPPFLAGS) -o test_input $(OBJS)

test_input.o : test_input.cpp
	$(CC) -c $(CPPFLAGS) test_input.cpp -o test_input.o

.PHONY : clean
clean:
	rm -f $(OBJS)