 SRC_PATH = src
 BUILD_PATH = bin
 VPATH = ../src
 LIB =  -lpthread -openmp -lX11 -lm
 #CC = /opt/intel/Compiler/11.1/073/bin/intel64/icpc
 CC = /opt/intel/Compiler/11.1/072/bin/intel64/icpc
 FLAG = -g

 all: bNptD bNptF bNptFOMP bSSEa bSSEaR bSSEI bSSEIR bLlw bL2w 

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

bSSEaR: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter_SSEassembly_reduction.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) 

bSSEI: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter_SSEintrin.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) 

bSSEIR: $(SRC_PATH)/CImg.h $(SRC_PATH)/def.h
	$(CC)  $(SRC_PATH)/blfilter_SSEintrin_reduction.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) 

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
clean:	
	rm $(BUILD_PATH)/*
