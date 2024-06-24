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

