/* LICENSE HEADER MANAGED BY add-license-header
 *
 *    Copyright 2024  Naohito Nakasato
 *    All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AN
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <assert.h>
#include <CL/cl.h>

namespace __CPP {
#include "cpp_kernel.h"

  class __CPP_Helper0 {
  public:
    __CPP_Helper0() { clinfo(); }; 
    ~__CPP_Helper0();

    void clinfo();
    //    void build();
    //    void build(std::string);
    void build(const char []);
    void setup(unsigned int ip, unsigned int id);      
    void setup();
    cl_kernel getkernel(const char *kernel_name);

  public:
    cl_context ctx;
    cl_command_queue q;
    std::vector<std::pair<int,int> > list;

  private:
    cl_program program;
    cl_platform_id platform_id[16];
    cl_kernel ker;
    cl_uint npl;
    cl_device_id device_id[16];
    cl_device_id dev;
    char pname[128];
    char dname[128];
    char pver[128];
    void dumperror();
  };

  __CPP_Helper0::~__CPP_Helper0 () {
    clReleaseKernel(ker);
    clReleaseProgram(program);
    clReleaseCommandQueue(q);
    clReleaseContext(ctx);
  }
  
  cl_kernel __CPP_Helper0::getkernel(const char *kernel_name)
  {
    cl_int status = CL_SUCCESS;
    cl_kernel res = clCreateKernel(program, kernel_name, &status);
    ker = res;

    if (status != CL_SUCCESS) {  
      fprintf(stderr, "Create kernel %s failed %d\n", kernel_name, status);
      exit(-1);
    }
    // https://arrayfire.com/blog/generating-ptx-files-from-opencl-code/
    {
      // Query binary (PTX file) size
      size_t bin_sz;
      cl_int err = clGetProgramInfo(program, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &bin_sz, NULL);

      // Read binary (PTX file) to memory buffer
      unsigned char *bin = (unsigned char *)malloc(bin_sz);
      err = clGetProgramInfo(program, CL_PROGRAM_BINARIES, sizeof(unsigned char *), &bin, NULL);

      // Save PTX
      char buf[1024];
      sprintf(buf, "%s_ocl.ptx", kernel_name);  
      //      FILE *fp = fopen(buf, "wb");
      //      fwrite(bin, sizeof(char), bin_sz, fp);
      //      fclose(fp);
      free(bin);
    }
    return res;
  }

  void __CPP_Helper0::setup()
  {
    char* p_type_str = getenv("OPENCL_PLATFORM");
    char* d_type_str = getenv("OPENCL_DEVICE");

    if (p_type_str == NULL || d_type_str == NULL) {
      setup(0, 0);
    } else {
      setup(atoi(p_type_str), atoi(d_type_str));
    }
  }

  void __CPP_Helper0::setup(unsigned int ip, unsigned int id)
  {
    cl_int result = CL_SUCCESS;
    std::cerr << "Selected: "; 
    if (npl <= ip) {
      fprintf(stderr, "FATAL: the specifed platform does not exist");
      exit(-1);      
    }
    cl_uint ndev = 0;
    if ((result = clGetDeviceIDs(platform_id[ip], CL_DEVICE_TYPE_ALL, 16, device_id, &ndev)) != CL_SUCCESS) {
      fprintf(stderr, "clGetDeviceIDs() failed : %d\n", result);
      exit(-1);
    }

    if (ndev <= id) {
      fprintf(stderr, "FATAL: the specifed device does not exist");
      exit(-1);      
    }

    clGetPlatformInfo(platform_id[ip],  CL_PLATFORM_NAME, sizeof(pname), pname, NULL);
    clGetPlatformInfo(platform_id[ip],  CL_PLATFORM_VERSION, sizeof(pver), pver, NULL);
    clGetDeviceInfo(device_id[id], CL_DEVICE_NAME, sizeof(dname), dname, NULL);

    std::cerr << pname << " " << pver << "::"	<< dname << "\n";

    // Create Context
    if( (ctx = clCreateContext(NULL, 1, &device_id[id], NULL, NULL, &result)) == NULL ) {
      fprintf(stderr, "Create context failed %d : dev %i\n", result, id);
      exit(-1);
    }

    // use the first device only
    cl_int status = CL_SUCCESS;
    q = clCreateCommandQueue(ctx, device_id[id], CL_QUEUE_PROFILING_ENABLE, &status);
    if (status != CL_SUCCESS) {  
      fprintf(stderr, "Create commandq failed %d\n", status);
      exit(-1);
    }
    dev = device_id[id];
  }


  void __CPP_Helper0::build(const char options[] = NULL)
  {
    cl_int status = CL_SUCCESS;
    {
      char *prog;
      prog = (char *)malloc(sizeof(__CPP::cpp_kernel_str));
      size_t ss[1];
      ss[0] = sizeof(__CPP::cpp_kernel_str);
      strcpy(prog, __CPP::cpp_kernel_str);
      program = clCreateProgramWithSource(ctx, 1, (const char **)&prog, ss, &status);
      if (status != CL_SUCCESS) {  
	fprintf(stderr, "cl create program failed %d\n", status);
	exit(-1);
      }
      free(prog);
    }
    status = clBuildProgram(program, 1, &dev, options, NULL, NULL);
    if(status != CL_SUCCESS) {
      fprintf(stderr, "build failed\n");
      dumperror();
      exit(-1);
    }
  }

  void __CPP_Helper0::dumperror()
  {
    cl_int logStatus;
    char * buildLog = NULL;
    size_t buildLogSize = 0;

    logStatus = clGetProgramBuildInfo (program,
                                       dev,
                                       CL_PROGRAM_BUILD_LOG,
                                       buildLogSize,
                                       buildLog,
                                       &buildLogSize);

    buildLog = (char*)malloc(buildLogSize);
    memset(buildLog, 0, buildLogSize);

    logStatus = clGetProgramBuildInfo (program,
                                       dev,
                                       CL_PROGRAM_BUILD_LOG,
                                       buildLogSize,
                                       buildLog,
                                       NULL);

    fprintf(stderr, "%s\n", buildLog);
    std::cout << logStatus << "\n";

    free(buildLog);
  }

  void __CPP_Helper0::clinfo() 
  {
    cl_int result = CL_SUCCESS;

    if ((result = clGetPlatformIDs(16, platform_id, &npl)) != CL_SUCCESS) {
      fprintf(stderr, "clGetPlatformIDs() failed : %dn", result);
      exit(-1);
    }

    for(unsigned int i = 0; i < npl; i++) {
      clGetPlatformInfo(platform_id[i],  CL_PLATFORM_NAME, sizeof(pname), pname, NULL);
      clGetPlatformInfo(platform_id[i],  CL_PLATFORM_VERSION, sizeof(pver), pver, NULL);
      fprintf(stderr, "platform %d %s %s\n", i, pname, pver);

      // Get Device ID
      cl_uint ndev = 0;
      if ((result = clGetDeviceIDs(platform_id[i], CL_DEVICE_TYPE_ALL, 16, device_id, &ndev)) != CL_SUCCESS) {
	fprintf(stderr, "clGetDeviceIDs() failed : %d\n", result);
	exit(-1);
      }

      for(unsigned int j = 0; j < ndev; j++) {
	clGetDeviceInfo(device_id[j], CL_DEVICE_NAME, sizeof(dname), dname, NULL);
	fprintf(stderr, "\tdevice %d %s\n", j, dname);
	list.push_back(std::make_pair(i, j));
      }
    }
  }

  static bool cl_first       = true;
  static bool cl_first_build = true;
  static __CPP_Helper0 *acc0;

  class __CPP_Helper {
  public:
    __CPP_Helper(int ip = 0, int id = 0) {
      if (cl_first) {
	acc0 = new __CPP_Helper0;
	if (ip == 0) acc0->setup();
	else acc0->setup(ip, id);
	cl_first = false;
      }
    }
    ~__CPP_Helper() {
      delete acc0;
    };

    void build(const char options[]) { 
      if (cl_first_build) {
	acc0->build(options);
	cl_first_build = false;
      }
    }
    
    cl_context getctx() {
      return acc0->ctx;
    }

    cl_command_queue getq() {
      return acc0->q;
    }

    cl_kernel getkernel(const char *kernel_name) {
      return acc0->getkernel(kernel_name);
    }
  };
}
