SRC=src/dfitter.cpp


dfitter: obj/dfitter.o
	g++  $^ -o $@
#libQCDNUM.so.0

test:  h1FitB/h12006flux.o h1FitB/qcd_2006.o h1FitB/i_2006_fitb.o h1FitB/i_2006_fita.o  obj/exampleCxx.o
	g++   $^ -Lqcdnum/lib -lQCDNUM  -Wl,-rpath=qcdnum/lib  -lgfortran -o $@


obj/dfitter.o: src/dfitter.cpp
	g++ -c $^ -o $@

obj/exampleCxx.o: src/exampleCxx.cc
	g++ -c -Iqcdnum/include  $^  -o $@


