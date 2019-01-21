SRC = src/dfitter.cpp src/pdf.cpp src/dplotter.cpp

#Convert src to obj
OBJtmp = $(subst src,obj,${SRC})
OBJ   = $(subst cpp,o,${OBJtmp})

H1fitobj = h1pdf2006/h12006flux.o h1pdf2006/qcd_2006.o h1pdf2006/i_2006_fitb.o h1pdf2006/i_2006_fita.o 
H1fitErrobj = h1pdf2006err/h12006flux.o h1pdf2006err/qcd_2006.o h1pdf2006err/i_2006_fitb.o h1pdf2006err/i_2006_fita.o 


#ROOT includes + libraries
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)

INC = ${ROOTCFLAGS} -Iinc -Iqcdnum/include -IPlottingHelper
LIBS=-Lqcdnum/lib -lQCDNUM  -Wl,-rpath=qcdnum/lib  -lgfortran  ${ROOTLIBS}  -Wl,-rpath,PlottingHelper -LPlottingHelper -lPlottingHelper



dfitter: ${OBJ} ${H1fitobj}
	g++  $^ ${LIBS}  -o $@
#libQCDNUM.so.0

test:  h1FitB/h12006flux.o h1FitB/qcd_2006.o h1FitB/i_2006_fitb.o h1FitB/i_2006_fita.o  obj/exampleCxx.o
	g++   $^ ${LIBS} -o $@

hfTest:    obj/exampleHF.o
	g++   $^ ${LIBS} -o $@


obj/%.o: src/%.cpp
	g++ -c ${INC} $^ -o $@

clean:
	rm obj/*.o




#obj/exampleCxx.o: src/exampleCxx.cc
	#g++ -c -Iqcdnum/include  $^  -o $@

