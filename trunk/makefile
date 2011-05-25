 SRC_PATH = src
 BUILD_PATH = bin
 VPATH = ../src
 LIB =  -lpthread -openmp -lX11 -lm
 CC = /opt/intel/Compiler/11.1/073/bin/intel64/icpc
 CC = icpc
 #CC = g++
 FLAG = -mia32 
 #FLAG = -m64 -O3 -fprefetch-loop-arrays -funroll-all-loops 
 #FLAG +=   -g  
 #FLAG = -vec_report2 -mcmodel=medium

 all: bSSEI bSSEIOMP bSSEIMI 
 #all: bfilter bfilterOMP bfilterTile bfilterTileOMP bfilterExpTile bfilterExpTileOMP

bSSEI: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter_SSEintrin.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) 

bSSEIMI: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter_SSEintrin.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) -DNO_COMP 

bSSEIOMP: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter_SSEintrin.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) -DUSE_OMP 

clean:	
	rm $(BUILD_PATH)/*
