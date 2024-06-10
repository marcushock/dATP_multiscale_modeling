NVCC=nvcc
ARCH=sm_35 #or sm_61
OBJ=$(patsubst src/%,obj/%,$(patsubst %.cu,%.o,$(wildcard src/*.cu)))
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
LIBS=-lcudadevrt -lboost_system-mt -lboost_thread-mt
NVCC_FLAGS=-rdc=true -arch=$(ARCH) --compiler-options -O3
endif
ifneq ($(UNAME), Darwin)
LIBS=-lcudadevrt -lboost_system -lboost_thread
NVCC_FLAGS=-rdc=true -arch=$(ARCH) --compiler-options -O3 -std=c++11
endif
BINNAME=MCMC_CUDA_10States

default: bin/$(BINNAME)
all: default
bin/$(BINNAME): $(OBJ)
	@mkdir -p bin
	$(NVCC) $(NVCC_FLAGS) $(LIBS) -o $@ $^
obj/%.o: src/%.cu
	@mkdir -p obj
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@
clean:
	@mkdir -p obj bin
	rm -r obj bin
