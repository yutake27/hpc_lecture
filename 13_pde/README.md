# Final Report README

1. implement in **c++**

   ```bash
   g++ 10_cavity.cpp
   ./a.out
   ```

   You will get `u.txt`. `v.txt`, `p.txt` and you can convert those text files into a image as follows.

   ```bash
   python 10_cavity_draw.py
   ```

   You wil get `10_cavity.png`.

   ```bash
   display 10_cavity.png
   ```

   ![](https://github.com/yutake27/hpc_lecture/blob/master/13_pde/10_cavity_cpp.png)

   (uploaded as 10_cavity_cpp.png)

   

2. implement in **openmp**

   ```bash
   g++ 10_cavity_openmp.cpp -fopenmp
   ./a.out
   ```

   You will get a text file and you covert them into a image  in the same way.

   ```bash
   python 10_cavity_draw.py
   display 10_cavity.png
   ```

   ![](https://github.com/yutake27/hpc_lecture/blob/master/13_pde/10_cavity_openmp.png)

   (uploaded as 10_cavity_openmp.png)

3. implement in **cuda**

   ```bash
   nvcc 10_cavity.cu -arch=sm_60 -O3
   ./a.out
   python 10_cavity_draw.py
   display 10_cavity.png
   ```

   ![](https://github.com/yutake27/hpc_lecture/blob/master/13_pde/10_cavity_cuda.png)

   (uploaded as 10_cavity_cuda.png)





