CC	= g++
CXX	= g++

OMP = -fopenmp

MACRO=-D___MPLAPACK_BUILD_WITH__POSIT__

POSITDIR=../SoftPosit
CPPFLAGS= -Wall -O3 $(OMP) $(MACRO) -std=gnu++17 -fpermissive -I. -I/usr/local/cuda/include \
          -I/opt/rocm/opencl/include/CL \
          -I$(POSITDIR)/source/include

LIBRPACK=libRpack.a

OBJ = util.o \
	Rgetrf.o Rgetrf2.o Rgetrs.o Rgetrf_opencl.o Rgetrf2_opencl.o \
	Rpotrf.o Rpotrf2.o Rpotrs.o Rpotrf_opencl.o \
	Mlsame.o iMlaenv.o Rlaswp.o Rscal.o iRamax.o \
	Rsyrk.o Rtrsm.o Mxerbla.o iMieeeck.o iMparmq.o Risnan.o \
	Rgemm.o	Rgemv.o \
	Rgemm_opencl.o

default:$(LIBRPACK)

clean:;
	rm -rf $(OBJ) cpp_kernel.h __tmp__.cl $(LIBRPACK)

Rgemm_opencl.o : cpp_kernel.h

cpp_kernel.h: __tmp__.cl
	./utils/template-converter cpp_kernel_str $< >| cpp_kernel.h

__tmp__.cl: kernel_entry.cl
	cat $^ >| __tmp__.cl

$(LIBRPACK): $(OBJ)
	ar rcs $@ $(OBJ)

