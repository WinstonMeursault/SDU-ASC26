# Commands

## 安装Nvidia驱动

### 更新系统

```bash
sudo apt update
sudo apt upgrade -y
```

### 安装必要依赖包

```bash
sudo apt install build-essential dkms
sudo apt install linux-headers-$(uname -r)
```

### 禁用默认Nouveau驱动

编辑文件

```bash
sudo nano /etc/modprobe.d/blacklist-nouveau.conf
```

输入

```text
blacklist nouveau
options nouveau modeset=0 
```

更新initramfs并重启

```bash
sudo update-initramfs -u
sudo reboot
```

### 安装Nvidia驱动及CUDA工具包

```bash
sudo apt install nvidia-driver-525 nvidia-cuda-toolkit
sudo reboot
```

### 测试

```bash
nvidia-smi
```

## 基准测试

注释的星号(*)后为省略的行为，如 vim 里的操作

### 建立工作区

此后所有操作主要发生在该目录下

```bash
cd ~
mkdir workplace
```

### 安装 HPL HPCG BLAS CBLAS

#### 安装工具并克隆 HPL HPCG BLAS CBLAS

```bash
sudo apt install git vim openmpi-bin openmpi-common libopenmpi-dev
git clone https://github.com/icl-utk-edu/hpl
git clone https://github.com/hpcg-benchmark/hpcg
wget http://www.netlib.org/blas/blas-3.8.0.tgz
wget http://www.netlib.org/blas/blast-forum/cblas.tgz
tar -xf blas-3.8.0.tgz
tar -xf cblas.tgz
```

#### 编译 BLAS

```bash
cd BLAS-3.8.0
make -j64
ar rv libblas.a *.o
```

#### 编译 CBLAS

```bash
cd ../CBLAS
cp ../BLAS-3.8.0/blas_LINUX.a ./
cp ./blas_LINUX.a ./lib
vim Makefile.in # 将 BLLIB 设置为 ../blas_LINUX.a
```

```bash
make -j64
./testing/xzcblat1 # 测试  (全 PASS)
cd ..
```

#### 安装 HPL

```bash
cd hpl
cp ./setup/Make.Linux_PII_FBLAS Make.test
sed -i 's/UNKNOWN/test/g' Make.top
sed -i 's/UNKNOWN/test/g' Makefile

vim Make.test
# 修改 Make.test：
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

make -j64 arch=test
```

#### 安装 HPCG

```bash
cd ../hpcg/setup

vim Make.Linux_MPI
# 修改 Make.Linux_MPI
# 注：同上，MPdir 可以不用管
# MPinc = -I/usr/include/x86_64-linux-gnu/mpi
# MPlib = /usr/lib/x86_64-linux-gnu/libmpi.so

cd ..
mkdir build
cd build
../configure Linux_MPI
make -j64
```

### 测试 HPL

```bash
cd ../../hpl/bin/test

# 调整 HPL.dat，通过 https://www.advancedclustering.com/act_kb/tune-hpl-dat-file/ 获得所需的 HPL.dat
mv HPL.dat HPL.dat.old
vim HPL.dat     # 粘贴所得 HPL.dat，将 device out 项改为 7，保存退出

# 运行 HPL
mpirun -np 64 ./xhpl > test.log
```

 随后下载 HPL.out

### 测试 HPCG

```bash
cd ../../../hpcg/build/bin
vim hpcg.dat    # 将第四行修改为 1800

# 运行 HPCG
mpirun -np 64 ./xhpcg
```

随后下载输出文件 `HPL.out`

## AMSS-NCKU

### 配置环境

#### 更新依赖包

```bash
sudo apt-get update
```

### 安装GCC/GFortran 编译器

```bash
sudo apt-get install gcc gfortran
```

#### 安装Make/Build 工具

```bash
sudo apt-get install make build-essential
```

#### 安装CUDA 编译器

```bash
sudo apt-get install nvidia-cuda-toolkit
```

#### 安装MPI工具

```bash
sudo apt install openmpi-bin libopenmpi-dev
```

#### 安装Python

```bash
sudo apt-get install python3 python3-pip
```

#### 安装OpenCV工具

```bash
sudo apt-get install libopencv-dev python3-opencv
```

#### 安装Python 的相关包

```bash
pip install numpy scipy matplotlib SymPy opencv-python torch
```

### AMSS-NCKU测试

```bash
python3 AMSS_NCKU_program.py
```

## 参考文献

1. [HPL+CBLAS+MPICH的搭建及测试流程-跨栏背心儿-FastEDA](https://www.fasteda.cn/post/163.html)
