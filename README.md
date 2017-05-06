# ICP
This is a Matlab implementation of the Iterative Closest Point algorithm described in *Besl, P.J. & McKay, N.D. 1992, 'A Method for Registration  of 3-D Shapes', IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 14, no. 2, IEEE Computer Society*.

It consists of two files, `icp.m` and `icp_driver.m`. `icp_driver.m` is an example of using the `icp` function.

## Use
`[R, t] = icp( srcPoints, targetPoints, tau )` iteratively seeks the rotation R and translation t which align srcPoints most closely with targetPoints. srcPoints and targetPoints are matrices of 3D points with a column for each point. Note that there is no requirement for srcPoints to have the same number of points as targetPoints. tau is a threshold such that when the RMSE between transformed srcPoints and targetPoints does not change by more than this threshold between iterations the algorithm is deemed to have converged.

`[R, t, k] = ICP( srcPoints, targetPoints, tau )` also returns k, the number of iterations required to complete the match.

## More Information
For more information see [The Besl and McKay paper](http://graphics.stanford.edu/courses/cs164-09-spring/Handouts/paper_icp.pdf)

