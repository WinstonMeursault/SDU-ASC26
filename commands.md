# Commands

```bash
cd ~
# 注释的星号(*)后为省略的行为，如 vim 里的操作

# 建立工作区，此后所有操作主要发生在该目录下
mkdir workplace

# 安装工具并克隆 HPL HPCG BLAS CBLAS
sudo apt install git vim mpich libmpich-dev
git clone https://github.com/icl-utk-edu/hpl
git clone https://gitub.com/hpcg-benchmark/hpcg
wget http://www.netlib.org/blas/blas-3.8.0.tgz
wget http://www.netlib.org/blas/blast-forum/cblas.tgz
tar -xf blas-3.8.0.tgz
tar -xf cblas.tgz

# 编译 BLAS
cd BLAS-3.8.0
make
ar rv libblas.a *.o

# 编译 CBLAS
cd ../CBLAS
cp ../BLAS-3.8.0/blas_LINUX.a ./
cp ./blas_LINUX.a ./lib
vim Makefile.in

# * 将 BLLIB 设置为 ../blas_LINUX.a
make
./testing/xzcblat1 # 测试 （全 PASS）
cd ..

# 安装 HPL
cd hpl
cp ./setup/Make.Linux_PII_FBLAS Make.test
vim Make.in
# * 将 arch 后的 UNKNOWN 修改为 test
vim Makefile
# * 将 arch 后的 UNKNOWN 修改为 test
vim Make.test
# * 修改 Make.test：
# ARCH = test
# （新建） WORKPLACE = /home/admin1/workplace
# TOPdir = ${WORKPLACE}/hpl
# 注：MPdir 仅用于定义 MPinc 和 MPlib，由于通过 apt 安装，所以可以不用管
# MPinc = -I/usr/include/x8664-linux-gnu/mpich
# MPlib = /usr/lib/x86_64-linux-gnu/libmpi.so
# LAdir = ${WORKPLACE}/CBLAS/lib
# LAlib = ${LAdir}/cblas_LINUX.a ${LAdir}/blas_LINUX.a
# CC = /usr/bin/mpicc
# LINKER = /usr/bin/mpif77
make arch=test

# 安装 HPCG
cd ../hpcg/setup
vim Make.Linux_MPI
# * 修改 Make.Linux_MPI
# 注：同上，MPdir 可以不用管
# MPinc = -I/usr/include/x86_64-linux-gnu/mpich
# MPlib = /usr/lib/x86_64-linux-gnu/libmpi.so
cd ..
mkdir build
cd build
../configure Linux_MPI
make

# 测试 HPL
cd ../../hpl/bin/test
# 调整 HPL.dat，通过 https://www.advancedclustering.com/act_kb/tune-hpl-dat-file/ 获得所需的 HPL.dat
mv HPL.dat HPL.dat.old
vim HPL.dat
# * 粘贴所得 HPL.dat，保存退出
# 运行 HPL
mpirun -np 4 ./xhpl > test.log
# * 下载 test.log

# 测试 HPCG
cd ../../../hpcg/build/bin
vim hpcg.dat
# * 将第四行修改为 1800
# 运行 HPCG
mpirun -np 16 ./xhpcg
# * 下载输出文件
```

## AMSS-NCKU

### 安装依赖包

#### Update the dependency packages

```bash
sudo apt-get update
```

#### Install the Make/Build Tool

```bash
sudo apt-get install make build-essential
```

#### Install the CUDA tool

```bash
sudo apt-get install nvidia-cuda-toolkit
```

#### Install the MPI  tool

```bash
sudo apt install openmpi-bin libopenmpi-dev
```

#### Install Python

```bash
sudo apt-get install python3 python3-pip
```

#### Install the OpenCV Tool

```bash
sudo apt-get install libopencv-dev
```

#### Install the relevant Python packages

```bash
pip install numpy scipy matplotlib SymPy opencv-python-full torch
```

## References

1. [HPL+CBLAS+MPICH的搭建及测试流程-跨栏背心儿-FastEDA](https://www.fasteda.cn/post/163.html)
