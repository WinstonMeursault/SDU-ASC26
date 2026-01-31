# 原文链接：[HPL Algorithm](https://www.netlib.org/benchmark/hpl/algorithm.html)

# 不同部分概览
1. Main Algorithm：介绍了如何原始矩阵分块并分配到各个进程中，并简要提及了原始矩阵将要被循环进行 LU 分解。

2. Panel Factorization：介绍了每一次迭代时平面矩阵组的计算。

3. Panel Broadcast：介绍了整个平面矩阵组由多个进程计算完成后，将如何进行进程间通讯以组合运算结果，总共提供了六种算法。
4. Look-ahead：介绍了剩余子矩阵可以提前更新为平面矩阵组以为将来的分解过程做准备。
5. Update：介绍了如何更新剩余子矩阵，并提供了两种算法。
6. Backward Substitution：介绍了在分解完毕后进行的回代过程。
7. Checking the Solution：介绍了如何验证答案的精确性。

# 术语解释
**注：** 术语翻译后标注星号，表示该术语笔者无法找到对应的翻译，本文章中此类术语的翻译由笔者自行翻译，因此不保证翻译质量。

## LU Factorization —— LU 分解
将矩阵 $\textbf{\textit{A}}$ 分解为**单位**下三角矩阵 $\textbf{\textit{L}}$ 和上三角矩阵 $\textbf{\textit{U}}$ 以及置换矩阵 $P$，即：
$$
\textbf{\textit{A}} = \textbf{\textit{PLU}}
$$
比如：
$$
\textbf{\textit{A}}
=
\begin{pmatrix}
1 & 2 & -1\\
2 & 1 & -2\\
-3 & 1 & 1\\
\end{pmatrix}
=
\begin{pmatrix}
1 & 0 & 0\\
2 & 1 & 0\\
-3 & -\frac{7}{3} & 1
\end{pmatrix}
\begin{pmatrix}
1 & 2 & -1\\
0 & -3 & 0\\
0 & 0 & -2
\end{pmatrix}
$$
可通过递归方式“从外到内一层一层地”进行 LU 分解，具体方法可参考：
具体方法参考：
1. [LU Factorization](https://netlib.org/utk/papers/factor/node7.html)
2. [线性方程组-直接法 1：Gauss消去与LU分解 - 知乎](https://zhuanlan.zhihu.com/p/386954541)

## Pivot —— 主元
消去过程中起主导作用的元素，如矩阵化为阶梯形后非零行第一个非零元素，参考[线性代数矩阵中主元和先导的区别是什么，主元的定义是什么？ - 知乎](https://www.zhihu.com/question/374820472)。

## Block-Cyclic Scheme (The Two-dimensional Block-Cyclic Distribution) —— 二维块循环数据分布
将矩阵分块后循环分配给不同的进程进行处理。
这种分配方式先将 $M \times N$ 的全局矩阵划分为大小为 $MB \times NB$ 的子块，然后把所有子块按照行列循环方式分配到 $MP \times NP$ 的进程网格（Process Grid）上分发给 $P$ 个进程。
如：
![[PixPin_2026-01-31_20-14-15.jpg]]
其中，每个数字代表分配的进程，相同的数字代表相同的进程。对于这个示例，$M = N = 16, MB = NB = 2, MP = NP = 2, P = 4$，这种分配方式实现了计算负载和通信开销的平衡。
参考：[The Two-dimensional Block-Cyclic Distribution](https://netlib.org/scalapack/slug/node75.html)

## Right-Looing Variant —— 右视向量法
同下文的 Crout 和 right-looking，都属于 LU 分解算法，具体内容可参考[Right-Looking Algorithm](https://netlib.org/utk/papers/outofcore/node3.html)

## Trailing Submatrix —— 剩余子矩阵\*
指原文左图中灰色的分块矩阵，至于为什么会对其进行更新，可参考[LU Factorization](https://netlib.org/utk/papers/factor/node7.html)

## Panel —— 平面矩阵组\*
指原文左图中蓝色部分，而红色部分可在蓝色部分计算完成后直接完成计算，具体可参考[LU Factorization](https://netlib.org/utk/papers/factor/node7.html)

### Cartesion Property —— 笛卡尔坐标系特性\*
此处指二维循环块分布具有笛卡尔坐标系特性，即进程分配遵循笛卡尔坐标系（平面直角坐标系）。

## Look-ahead Pipe：前瞻管道\*
在当前平面矩阵组尚未完成分解时，提前将下面几个将要计算的平面矩阵组进行缓存并存入该管道中。至于为什么更新剩余子矩阵需要计算，可参考[LU Factorization](https://netlib.org/utk/papers/factor/node7.html)