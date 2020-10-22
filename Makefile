
CXXFLAGS = -std=c++14  -pthread -O3
SEQAN = -I ~/work/software/seqan/include/ 

bcmap: barcodeMapper.cpp functions.o
	g++ barcodeMapper.cpp -o bcmap $(CXXFLAGS) $(SEQAN) functions.o
functions.o: src/functions.cpp
	g++ -c src/functions.cpp $(CXXFLAGS) $(SEQAN)
clean:
	rm functions.o bcmap
