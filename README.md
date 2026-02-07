# SDU-ASC26 | 山东大学（威海） ASC26 项目

## 项目目录说明

```text
|-- attachments         预赛准备打包提交的文件
|   |-- HPCG
|   |-- HPL
|
|-- ioFilesLog          输入输出文件记录
|   |-- HPCG
|       |-- [N]hpcg.dat     第 N 次测试参数
|       |-- [N]*.txt        第 N 次测试结果
|   |-- HPL
|       |-- output_file_note.md     解释输出文件中 T/V 列字符串含义
|       |-- 测试流程.md     记录参数设置和测试的过程
|       |-- [01]HPL.dat     Auto-Generated, OLD
|       |-- [02]HPL.dat     Auto-Generated, 128 BlockSize
|       |-- [0N]HPL.dat     第 N 次测试参数
|       |-- [03]参数解释.md     解释 [03]HPL.dat 中参数的设定
|       |-- [N]HPL.out      第 N 次测试结果
|
|-- NR-amssncku         AMSS-NCKU 优化项目
|   |-- amssncku        项目代码和说明文件
|   |-- AMSS-NCKU-Python Debug in Ubuntu2204.pdf
|   |-- How to run AMSS-NCKU-Python in Ubuntu2204.pdf
|   |-- README.md
|
|-- unitreeWorldModel           Unitree World Model Action 优化项目
|       |-- ASC26EWMO_source    ASC26测试集
|       |-- unifolmWorldModelAction
|       |-- README.md
|
|-- LICENSE         开源许可证
|-- README.md       项目说明
|-- commands.md     命令执行参考
|-- HPL Algorithm Notes.md      HPL 算法文档注释
|-- ASC26_Preliminary_Round_Notification.pdf    初赛通知
```

## 基准测试成果

### HPL

- 01 OLD 25.699GFLOPS

### HPCG

- 01 OLD 8.51662 GFLOPS
- 01 59GB Used Memory 26.3280 GFLOPS
