# CurveFitting using C++ g2o 

Optimization is a very important problem in any fields. To solve this problem using MATLAB [1] is an easy task, furthermore its user-friendly interface and clear documentation attracts many users. However, MATLAB can’t be used in a production code such as C++. In the past, to implement an algorithm to solve non linear least square problem usually adopts the code from the book of Numerical Recipes [2]. That implementation has a huge problem where one needs to provide the derivative of the cost function. If the cost function to optimize has n parameters, one need to derive its n partial derivative equation plus its implementation code, thus it is very error prone. 
Now, g2o [3] open source library solves this issue. However, many curve fitting examples [4,5] found solves only 3 parameters. In this project, an example to solve 5 parameters.

This example here shows fitting a 3D gaussian curve to estimate 5 parameters. MATLAB is used to generate data points and then a nonlinear least square method is used to estimate the parameters. The result matches what g2o produces.

![Model](https://user-images.githubusercontent.com/80547721/153137026-b692768c-087a-4025-8fc7-e8c4f8409337.jpg)

Figure 1 - (a) Model of 3D Gaussian (b) With Noise

![MATLAB_Result](https://user-images.githubusercontent.com/80547721/153137362-08920e7c-8055-40cd-909d-bd6a579a8297.PNG)

Figure 2 - Result using MATLAB 

![ResultCurveFitting](https://user-images.githubusercontent.com/80547721/153137064-e8dbe005-7ef7-4ebb-a8c8-c30215b10e3a.PNG)

Figure 3 - Result from g2o


Prerequisite:
- eigen: 3.1.4(use this version for the g2o below)
- g2o:  git checkout 20170730_git


[1]https://www.mathworks.com/help/optim/ug/lsqnonlin.html?searchHighlight=lsqnonlin&s_tid=srchtitle_lsqnonlin_1

[2] W. H., Teukolsky, S. A., Vetterling, W. T. and Flannery, B. P. [2002]. Numerical Recipes in C++: The Art of Scientific Computing, 2nd edn, Cambridge University
Press.

[3] https://github.com/raulmur/ORB_SLAM2

[4] https://programmerall.com/article/51531612288/

[5] https://github.com/RainerKuemmerle/g2o/blob/master/g2o/examples/data_fitting/curve_fit.cpp

