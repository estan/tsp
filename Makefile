all:
	g++ -g -O2 -std=gnu++0x -Wall -pedantic -o tsp tsp.cpp

clean:
	@-rm -rf tsp > /dev/null 2>&1
