 SRC_PATH = src
 BUILD_PATH = bin
 VPATH = ../src
 LIB =  -lpthread -openmp -lX11 -lm
 CC = /home/abi/intel/Compiler/11.1/073/bin/ia32/icpc

all: bNptD bNptF bNptFOMP bSSEa bSSEaR bSSEI bSSEIR bLlw bL2w 

bNptD: $(SRC_PATH)/CImg.h
	$(CC)  $(SRC_PATH)/blfilter_noopt_double.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) 

bNptF: $(SRC_PATH)/CImg.h
	$(CC)  $(SRC_PATH)/blfilter_noopt_float.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) 

bNptFOMP: $(SRC_PATH)/CImg.h
	$(CC)  $(SRC_PATH)/blfilter_noopt_float.cpp -o $(BUILD_PATH)/$@ $(FLAG) -DUSE_OMP $(LIB)

bSSEa: $(SRC_PATH)/CImg.h
	$(CC)  $(SRC_PATH)/blfilter_SSEassembly.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) 

bSSEaR: $(SRC_PATH)/CImg.h
	$(CC)  $(SRC_PATH)/blfilter_SSEassembly_reduction.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) 

bSSEI: $(SRC_PATH)/CImg.h
	$(CC)  $(SRC_PATH)/blfilter_SSEintrin.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) 

bSSEIR: $(SRC_PATH)/CImg.h
	$(CC)  $(SRC_PATH)/blfilter_SSEintrin_reduction.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) 

bLlw: $(SRC_PATH)/CImg.h
	$(CC)  $(SRC_PATH)/blfilter_lpunroll_linearwindow.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) 

bL2w: $(SRC_PATH)/CImg.h
	$(CC)  $(SRC_PATH)/blfilter_lpunroll_2dwindow.cpp -o $(BUILD_PATH)/$@ $(FLAG) $(LIB) 
	
clean:	
	rm $(BUILD_PATH)/*
