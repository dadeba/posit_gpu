Posit for GPUs

1. Clone SoftPosit
 git clone https://gitlab.com/cerlane/SoftPosit.git

2. Patch and build SoftPosit
 cd SoftPosit
 patch -p1 < ../SoftPosit.patch
 cd build/Linux-x86_64-GCC
 make

3. Build 
 make

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
