CC	= g++
CXX	= g++

RPACKDIR=./libRpack
POSITDIR=./SoftPosit

OMP = -fopenmp

MACRO=

CPPFLAGS= -Wall -O3 $(OMP) $(MACRO) -std=gnu++17 -fpermissive \
	-I/usr/local/cuda/include -I$(POSITDIR)/source/include -I$(RPACKDIR)
LFLAGS = $(OMP) -lm -L$(RPACKDIR) -lRpack -lOpenCL -L/usr/local/cuda/lib64 -L/opt/rocm/opencl/lib \
	$(POSITDIR)/build/Linux-x86_64-GCC/softposit.a -llapack -lopenblas

subdirs := $(RPACKDIR)
.PHONY: default $(subdirs)

default: $(subdirs) run_cho run_lu run_gemm run_gemm_trailing run_lu_bench run_lu_power_bench run_cho_bench \
         run_cho_check run_lu_check

$(subdirs):
	$(MAKE) -C $@

run_cho: main_cholesky.o $(RPACKDIR)/libRpack.a
	$(CC) -o $@ $< $(LFLAGS)

run_cho_check: main_cholesky_check.o $(RPACKDIR)/libRpack.a
	$(CC) -o $@ $< $(LFLAGS)

run_lu: main_lu.o $(RPACKDIR)/libRpack.a
	$(CC) -o $@ $< $(LFLAGS)


run_lu_check: main_lu_check.o $(RPACKDIR)/libRpack.a
	$(CC) -o $@ $< $(LFLAGS)

run_lu_bench: main_lu_bench.o $(RPACKDIR)/libRpack.a
	$(CC) -o $@ $< $(LFLAGS)

run_lu_power_bench: main_lu_power_bench.o $(RPACKDIR)/libRpack.a
	$(CC) -o $@ $< $(LFLAGS)

run_cho_bench: main_cholesky_bench.o $(RPACKDIR)/libRpack.a
	$(CC) -o $@ $< $(LFLAGS)

run_gemm: main_gemm.o $(RPACKDIR)/libRpack.a
	$(CC) -o $@ $< $(LFLAGS)

run_gemm_trailing: main_gemm_trailing.o $(RPACKDIR)/libRpack.a
	$(CC) -o $@ $< $(LFLAGS)

clean:;
	rm -rf $(OBJ) *.o run_* *.ptx

