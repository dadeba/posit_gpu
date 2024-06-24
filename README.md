# GEMM Routines in Posit for GPUs

We have ported the addition and multiplication routines from [SoftPosit](https://gitlab.com/cerlane/SoftPosit.git) as OpenCL kernels. We also created GEMM routines in 32-bit Posit arithmetic. These programs were used for benchmarking in our paper presented at HPC Asia 2024.

Part of the code is derived from [MPLAPACK](https://github.com/nakatamaho/mplapack).

## Build Instructions

1. Clone the SoftPosit repository:
    ```bash
    $ git clone https://gitlab.com/cerlane/SoftPosit.git
    ```

2. Apply the patch and build SoftPosit:
    ```bash
    $ cd SoftPosit
    $ patch -p1 < ../SoftPosit.patch
    $ cd build/Linux-x86_64-GCC
    $ make
    $ cd ../../..
    ```

3. Build the project:
    ```bash
    $ make
    ```

## Test Programs

All programs use GEMM routines in 32-bit Posit arithmetic. You can specify the blocking size for the GEMM routines through an environment variable.

Example of setting the block size to 16:
```bash
$ export OPENCL_GEMM_BLOCKSIZE=16
```

The performance of all programs can vary slightly depending on the blocking size.

## GEMM
- run_gemm
- run_gemm_trailing

## LU decomposition
- run_lu
- run_lu_bench
- run_lu_check
- run_lu_power_bench

## Cholesky decomposition
- run_cho 
- run_cho_bench
- run_cho_check

# Reference 
```bibtex
@inproceedings{10.1145/3635035.3635046,
author = {Nakasato, Naohito and Murakami, Yuki and Kono, Fumiya and Nakata, Maho},
title = {Evaluation of POSIT Arithmetic with Accelerators},
year = {2024},
isbn = {9798400708893},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
url = {https://doi.org/10.1145/3635035.3635046},
doi = {10.1145/3635035.3635046},
booktitle = {Proceedings of the International Conference on High Performance Computing in Asia-Pacific Region},
pages = {62–72},
numpages = {11},
location = {Nagoya, Japan},
series = {HPCAsia '24}
}
```
