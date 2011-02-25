 SRC_PATH = src
 BUILD_PATH = bin
 VPATH = ../src
 LIB =  -lpthread -openmp -lX11 -lm
 CC = /opt/intel/Compiler/11.1/073/bin/intel64/icpc
 CC = icpc
 #CC = g++
 FLAG = -xSSE2 -O2
 #FLAG = -m64 -O3 -fprefetch-loop-arrays -funroll-all-loops 
 #FLAG =   -g  
 #FLAG = -vec_report2 -mcmodel=medium

 all: bSSEa bSSEaOMP bSSEaR bSSEaROMP bSSEI bSSEIOMP bSSEIR bSSEIROMP bfilterTile bfilterTileOMP bLlw bL2w bNptD bNptF bNptFOMP 
 #all: bfilter bfilterOMP bfilterTile bfilterTileOMP bfilterExpTile bfilterExpTileOMP

bNptD: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter_noopt_double.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) 

bNptDOMP: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter_noopt_double.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) -DUSE_OMP 

bNptF: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter_noopt_float.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) 

bNptFOMP: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter_noopt_float.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) -DUSE_OMP 

bSSEa: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter_SSEassembly.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) 

bSSEaOMP: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter_SSEassembly.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) -DUSE_OMP

bSSEaR: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter_SSEassembly_reduction.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) 

bSSEaROMP: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter_SSEassembly_reduction.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) -DUSE_OMP 

bSSEI: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter_SSEintrin.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) 

bSSEIMI: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter_SSEintrin.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) -DNO_COMP 

bSSEIOMP: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter_SSEintrin.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) -DUSE_OMP 

bSSEIR: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter_SSEintrin_reduction.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) 

bSSEIROMP: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter_SSEintrin_reduction.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) -DUSE_OMP 

bLlw: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter_lpunroll_linearwindow.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) 

bL2w: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter_lpunroll_2dwindow.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) 

bfilter: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) 

bfilterOMP: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) -DUSE_OMP 
	
bfilterTile: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) -DUSE_TILES 

bfilterTileOMP: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) -DUSE_TILES -DUSE_OMP 

bfilterExpTile: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) -DEXP_TILES 

bfilterExpTileOMP: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) -DEXP_TILES -DUSE_OMP 
clean:	
	rm $(BUILD_PATH)/*
