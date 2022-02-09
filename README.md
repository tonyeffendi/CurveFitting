# CurveFitting

Optimization is a very important problem in any fields. To solve this problem using MATLAB [1] is an easy task, furthermore its user-friendly interface and clear documentation attracts many users. However, MATLAB can’t be used in a production code such as C++. In the past, to implement an algorithm to solve non linear least square problem usually adopts the code from the book of Numerical Recipes [2]. That implementation has a huge problem where one needs to provide the derivative of the cost function. If the cost function to optimize has n parameters, one need to derive its n partial derivative equation plus its implementation code, thus it is very error prone. 
Now, g2o [3] open source library solves this issue. However, many curve fitting examples [4,5] found solves only 3 parameters. In this project, an example to solve 5 parameters.

[1]https://www.mathworks.com/help/optim/ug/lsqnonlin.html?searchHighlight=lsqnonlin&s_tid=srchtitle_lsqnonlin_1

[2] W. H., Teukolsky, S. A., Vetterling, W. T. and Flannery, B. P. [2002]. Numerical Recipes in C++: The Art of Scientific Computing, 2nd edn, Cambridge University
Press.

[3] https://github.com/raulmur/ORB_SLAM2
[4] https://programmerall.com/article/51531612288/
[5] https://github.com/RainerKuemmerle/g2o/blob/master/g2o/examples/data_fitting/curve_fit.cpp
