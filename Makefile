SRC = src/dfitter.cpp src/pdf.cpp
INC = -Iinc -Iqcdnum/include

#Convert src to obj
OBJtmp = $(subst src,obj,${SRC})
OBJ   = $(subst cpp,o,${OBJtmp})

fitBobj = h1FitB/h12006flux.o h1FitB/qcd_2006.o h1FitB/i_2006_fitb.o h1FitB/i_2006_fita.o 

LIBS=-Lqcdnum/lib -lQCDNUM  -Wl,-rpath=qcdnum/lib  -lgfortran

dfitter: ${OBJ} ${fitBobj}
	g++  $^ ${LIBS}  -o $@
#libQCDNUM.so.0

test:  h1FitB/h12006flux.o h1FitB/qcd_2006.o h1FitB/i_2006_fitb.o h1FitB/i_2006_fita.o  obj/exampleCxx.o
	g++   $^ ${LIBS} -o $@

obj/%.o: src/%.cpp
	g++ -c ${INC} $^ -o $@




#obj/exampleCxx.o: src/exampleCxx.cc
	#g++ -c -Iqcdnum/include  $^  -o $@

