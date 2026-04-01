# Numerical-Algebra
数值代数相关的程序实现，仅供理论学习。

本项目为学习用途，仅小范围使用Copilot自动补全。

## 文件说明

`Gauss_Elim.c`：简单高斯消元法化简增广矩阵。在代码中的编辑区输入矩阵行数和矩阵，运行脚本进行化简。矩阵的非增广部分必须满秩且消元过程主元位置不能出现0。

`Gauss_Elim_Col_Pivot.c`：列主元高斯消元法化简增广矩阵。在代码中的编辑区输入矩阵行数和矩阵，运行脚本进行化简。脚本在每次消元时会通过行互换使得主元的绝对值尽可能大并避免主元为零造成停机。矩阵的非增广部分必须满秩。

`Doolittle.c`：矩阵的直接Doolittle分解（直接LU分解）。数学上可以证明能够进行分解的矩阵必须所有顺序主子式非零。

`LU_Decomp_and_solving.c`：在编辑区输入代表线性方程组`AX=b`的增广矩阵，使用LU分解对系数矩阵`A`进行LU分解并用两步求解方程组。输出内容中，`Y`为中间向量，满足`LY=b`，`X`为最终解向量，通过求解`UX=Y`得出. 

`Cholesky.c`：对称正定矩阵的Cholesky分解。

文件夹`LU_decomp`：LU分解的MATLAB实现。这是一项课程作业，否则我不会去写已经内置好的东西。

`Jacobi_Iter.c`：Jacobi迭代法解方程组。

`Gauss_Seidel_Iter.c`：Gauss-Seidel迭代法解方程组。

`Relaxation_Iter.c`：松弛迭代法解方程组。

`Steepest_Descent.c`：最速下降法解方程组。

`spsMatrix.c`：C实现稀疏矩阵存储读取解决方案。`main`函数为一个示例。

## 关于Library文件夹

数值代数功能实现的库函数封装。

使用方法：

1. 将`numalg.h`与`numalg.a`下载并放入选定的文件夹中

2. 打开Dev C++，点击左上角菜单栏File，依次点击New, Project，新建Console Application的C Project，在弹出的路径选择窗口中将.dev文件保存至1.中所述文件夹

3. 保存成功后出现`main.c`脚本标签页及预设模版。从上方菜单栏开始依次点击Project, Project Options, Parameters, Add library or object，在路径选择窗口转到上述文件夹，选择`numalg.h`与`numalg.a`，将它们添加到`Linker`栏中并确认

4. 回到`main.c`的脚本编辑界面，在`main.c`中导入头文件。即将`#include "numalg.h"`添加在脚本的开头

5. 开始编写你的脚本

库函数封装工作仍在进行中，本文件夹会不定时更新。可以打开`numalg.h`文件查看当前可用的函数。
