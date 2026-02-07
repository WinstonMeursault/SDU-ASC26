# III. Technical Section
## 2. HPL and HPCG Benchmarks
### Software Environment
| Key                    | Value              |
| ---------------------- | ------------------ |
| Operating System       | Ubuntu 20.04.6 LTS |
| Compilers              | GCC 9.4.0          |
| Mathematical Libraries | BLAS 3.8.0, CBLAS  |
| MPI Implementation     | Open MPI 4.0.3     |
### HPL Part
At first, we selected a problem size of 115k, which could occupy 80% of our server's memory, and other parameters were kept at their original values. Thus, the result was not very satisfying. Since a full scale test is time-consuming—about 6 hours for one test—we follow a strategy assuming that some smaller-scale tests with different parameters are generally equal to their huge-scale version, except the scale. Therefore, results from several smaller-scale tests can reveal the comparative relationship between their huge-scale version. To prove the correctness of this method, we constructed two sets of tests: one set is small-scale and another is huge-scale. Each set contained two tests where only one parameter is different. The result showed the comparative relationship between the small-scale tests were equivalent to the huge-scale version. This enabled us to optimize parameters using this method.
By controlling variables, we constructed some sets of small-scale tests, and each set can find a better parameters in different section of the input file. After the tests were all completed, we find some parameter candidates that has the potentially best performance. Testing all of them yielded a far more satisfying parameter set, whose result has a 0.457 improvement in Gflops over the original input file.
### HPCG Part
Employing the same testing method in HPL, we tried different problem size with a smaller running time of 60s. Then we found the best parameter and used it to have a full time test.