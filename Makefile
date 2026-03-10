NVCC=nvcc
ARCH= sm_86 #or sm_61 _35
OBJ=$(patsubst src/%,obj/%,$(patsubst %.cu,%.o,$(wildcard src/*.cu)))
UNAME := $(shell uname)
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


# New executable for coop testing
COOP_BIN=bin/coop_testing
COOP_SRC=tests/coop_testing.cu src/compute_coop_factor.cu

coop_test: $(COOP_SRC)
	@mkdir -p bin
	$(NVCC) $(NVCC_FLAGS) $(LIBS) -o $(COOP_BIN) $(COOP_SRC)

UPDATE_RUS_TEST_BIN=bin/update_RUs_helper_test
UPDATE_RUS_TEST_SRC=tests/update_RUs_helper_test.cu src/update_RUs.cu

update_rus_test: $(UPDATE_RUS_TEST_SRC)
	@mkdir -p bin
	$(NVCC) $(NVCC_FLAGS) $(LIBS) -o $(UPDATE_RUS_TEST_BIN) $(UPDATE_RUS_TEST_SRC)

run_update_rus_test: update_rus_test
	./$(UPDATE_RUS_TEST_BIN)

# Docker compilation target
dockercompile:
	docker run --user $(id -u):$(id -g) --gpus all -v $(shell pwd):/workspace -it my-cuda-boost:12.2 bash -c "cd /workspace && make clean && make all"
