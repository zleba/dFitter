SRC=src/dfitter.cpp


dfitter: obj/dfitter.o
	g++  $^ -o $@

obj/dfitter.o: src/dfitter.cpp
	g++ -c $^ -o $@

