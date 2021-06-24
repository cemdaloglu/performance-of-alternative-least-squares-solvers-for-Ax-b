# performance-of-alternative-least-squares-solvers-for-Ax-b

# # Task
In this project you will investigate the performance of alternative least squares solvers for Ax = b. The
numerical techniques to be compared are

• Generalized Minimum Residual (GMRES)

• Conjugate Gradients (CG)

• A+: The Pseudoinverse

The first step is to generate a symmetric matrix A according to the example attached in the following pages
(from the book Numerical Linear Algebra by Lloyd N. Trefethern and David Bau). You will use the values
![image](https://user-images.githubusercontent.com/36455629/123243324-93fa7580-d4eb-11eb-88c2-8355974e983c.png) and for each of them you need to create a matrixAwith dimensions 100x100, 500x500
and 10000 x 10000 (if the computations with this size create problems in your computational tools, find the
largest matrix size for which you can still obtain the results) as explained in the example. So, in total there will
be 6 different matrices A.

Now, for eachAwith dimensions mxm, the next step will be to generate 10 samples of vector x_0 by using
randn(m; 1). Note that the vector ![image](https://user-images.githubusercontent.com/36455629/123243534-ca37f500-d4eb-11eb-90cd-5522a5f18166.png) is sampled from a standard normal distribution. Then, you need to
compute b_0 = A x_0 and subsequently create the vector b = b_0 + w, where ![image](https://user-images.githubusercontent.com/36455629/123243662-e8055a00-d4eb-11eb-9418-6e6d91adb776.png) is a vector
sampled from a normal distribution with standard deviation ![image](https://user-images.githubusercontent.com/36455629/123243717-f3588580-d4eb-11eb-8cca-100bdb7a1fb0.png). You will do this for three different values of
![image](https://user-images.githubusercontent.com/36455629/123243746-f9e6fd00-d4eb-11eb-8a20-845a3b290efc.png), namely 0:0001; 0:01; 1. So in total, there will be 30 different realizations of the vector b, corresponding to
the 3 different noise levels. Let us denote them by ![image](https://user-images.githubusercontent.com/36455629/123243828-0b300980-d4ec-11eb-89ad-d2879cf44ceb.png) Note that the value i is
the noise level index and j is the index of the realizations of x_0.

At this point, you will employ the three methods in order to approximate the least squares solution ![image](https://user-images.githubusercontent.com/36455629/123243941-226ef700-d4ec-11eb-939d-8b97eb47a619.png) to
the system ![image](https://user-images.githubusercontent.com/36455629/123243985-2bf85f00-d4ec-11eb-9335-81e02168f559.png) and define the error vector ![image](https://user-images.githubusercontent.com/36455629/123244052-39154e00-d4ec-11eb-9c85-e7ed11a21dd3.png). Now, for ![image](https://user-images.githubusercontent.com/36455629/123244097-43374c80-d4ec-11eb-9213-27f2055aea30.png), we define

![image](https://user-images.githubusercontent.com/36455629/123244146-4d594b00-d4ec-11eb-850f-a00b2eef416c.png)

![image](https://user-images.githubusercontent.com/36455629/123244192-564a1c80-d4ec-11eb-977b-25f1bdd19e30.png) corresponds to the root-mean squared error on the estimated solution, while ![image](https://user-images.githubusercontent.com/36455629/123244226-5ea25780-d4ec-11eb-9746-7ce8f61220b0.png) corresponds to the rootmean
squared error in the fit to the observations. For each different A and the three methods, you will plot the
log values of ![image](https://user-images.githubusercontent.com/36455629/123244268-6b26b000-d4ec-11eb-8b3d-436972f7fc8e.png) with respect to the three values of ![image](https://user-images.githubusercontent.com/36455629/123244298-74178180-d4ec-11eb-8f83-1444eac33d3a.png). What do you observe? Comment on the results. In
the case of the Conjugate Gradients you will also plot the values of ![image](https://user-images.githubusercontent.com/36455629/123244338-81347080-d4ec-11eb-9956-5287ed147def.png) with respect to those of ![image](https://user-images.githubusercontent.com/36455629/123244364-885b7e80-d4ec-11eb-9400-aeba636629f0.png), as the plot
shown in the example below. Moreover, when you are using the GMRES and Conjugate Gradient approaches:

1. Comment on the stopping criteria you used and investigate what happens if you choose less or significantly
more iterations.

2. Check the orthogonality of the bases the algorithms generate for the Krylov subspace along with the
iterations.

For efficient solution to the tridiagonal system generated in GMRES algorithm, you can use efficient algorithms
such as Tridiagonal matrix algorithm.

## Example
![image](https://user-images.githubusercontent.com/36455629/123244560-b5a82c80-d4ec-11eb-9703-2bc6a73163de.png)

![image](https://user-images.githubusercontent.com/36455629/123244581-be98fe00-d4ec-11eb-9115-dad48ce4587e.png)

# # Solution
In this homework report, three different least square solvers are investigated and discussed. A, x and b are described as expected that is told in the homework statement.

# Generalized Minimum Residual (GMRES)

As stated in the lecture notes 2,

![image](https://user-images.githubusercontent.com/36455629/123245004-1b94b400-d4ed-11eb-8664-a24bcf11fd31.png)

Where ![image](https://user-images.githubusercontent.com/36455629/123245083-2c452a00-d4ed-11eb-8d8a-953fdb45c3fa.png). To find ![image](https://user-images.githubusercontent.com/36455629/123245123-336c3800-d4ed-11eb-97c7-fad6b9174933.png)
, QR decomposition can be applied since it is smaller least squares problem. The following figure shows the steps to find ![image](https://user-images.githubusercontent.com/36455629/123245158-39faaf80-d4ed-11eb-8605-f55adcc1e024.png) the which minimizes the given equation.

![image](https://user-images.githubusercontent.com/36455629/123245204-47179e80-d4ed-11eb-8786-906d6687db03.png)

In the Figure 1, step 2 shows that QR decomposition is applied to ![image](https://user-images.githubusercontent.com/36455629/123245328-6282a980-d4ed-11eb-945c-27c61fa3699a.png)
 and obtained Q and R where Q is a unitary matrix and R is an upper triangular matrix. R has form ![image](https://user-images.githubusercontent.com/36455629/123245375-6dd5d500-d4ed-11eb-95e0-51b308d0a2c8.png)
. Step 3 shows the step 1’s representation with QR decomposition. The final equation of the step 3’s right hand side equation is partitioned in the step 4. Then, step 5 is obtained by combining last two steps. Finally,

![image](https://user-images.githubusercontent.com/36455629/123245425-7d551e00-d4ed-11eb-908e-81497612309f.png)

For the residual ![image](https://user-images.githubusercontent.com/36455629/123245487-8e059400-d4ed-11eb-9853-7dc571700c70.png).

Since, Arnoldi method gives the ![image](https://user-images.githubusercontent.com/36455629/123245536-9b228300-d4ed-11eb-8994-0374786fd8c1.png) is found via QR decomposition method ![image](https://user-images.githubusercontent.com/36455629/123245602-aaa1cc00-d4ed-11eb-9559-70cbf612304b.png)
 can be computed. This is the GMRES algorithm, to summarize the steps of the algorithm;

* ![image](https://user-images.githubusercontent.com/36455629/123246378-6531ce80-d4ee-11eb-9de2-6f9cbea351c4.png) are computed by using Arnoldi’s method where ![image](https://user-images.githubusercontent.com/36455629/123246414-7084fa00-d4ee-11eb-8558-a037b3507055.png) are the new columns of ![image](https://user-images.githubusercontent.com/36455629/123246443-78dd3500-d4ee-11eb-937e-27ed2441240c.png).
* QR decomposition is applied to the new column of ![image](https://user-images.githubusercontent.com/36455629/123246509-8a264180-d4ee-11eb-9bb4-5e309eaa45f4.png)
 and does not change the previous columns.
* Partition is applied to the new column of Q since previous columns’ partition is already done.
* Algorithm stops if the final iteration is applied or the error rate is significantly (can be
changed by the user) small.

As a stopping criterion I checked if the k'th least square error ![image](https://user-images.githubusercontent.com/36455629/123246701-c0fc5780-d4ee-11eb-8d29-8c5e472fcbb1.png) is sufficiently small (10^(-15)
for instance) or if the iterations are over. If significantly smaller stop criteria is chosen and number of
iterations are increased, the error will be lessen and the predicted x will be closer to the actual result.

For the Krylov subspace the algorithm generates ![image](https://user-images.githubusercontent.com/36455629/123246815-e0938000-d4ee-11eb-8f0d-9269430dcf67.png)
 vectors in the lecture notes. ![image](https://user-images.githubusercontent.com/36455629/123246850-ebe6ab80-d4ee-11eb-9802-3fa673d3d9ae.png)
 thus, the bases are orthogonal. However, I could not solve the problem in my GMRES function so I used builtin
gmres function in scipy library.

The following Figures 2, 3 are obtained for this method.

![image](https://user-images.githubusercontent.com/36455629/123246904-fe60e500-d4ee-11eb-8f81-1b15e4122afe.png)

![image](https://user-images.githubusercontent.com/36455629/123246927-0587f300-d4ef-11eb-9c62-cc09909b1dc5.png)

For the Figures 2, 3, 4 and 5, ![image](https://user-images.githubusercontent.com/36455629/123247008-16386900-d4ef-11eb-956b-00fb5fb99912.png)
 noise level indexes are 1, 2 and 3 respectively. Therefore, it can be said that smaller noise result as less error for both ![image](https://user-images.githubusercontent.com/36455629/123247056-23edee80-d4ef-11eb-81c7-1180241f2014.png). Moreover, as the size of A increases, both errors decreases and the computed x gets closer to the actual x. The impact of ![image](https://user-images.githubusercontent.com/36455629/123247091-2d775680-d4ef-11eb-8111-0c59143217e0.png)
 can only be seen in Figure 2 and smaller ![image](https://user-images.githubusercontent.com/36455629/123247102-30724700-d4ef-11eb-8827-3ef1a591203d.png)
 has lower error.

# Conjugate Gradient

This method is applied exactly as in the second lecture notes (page 32). As a stopping criterion I checked if ![image](https://user-images.githubusercontent.com/36455629/123247754-e2117800-d4ef-11eb-877a-54226b04cdc1.png)
 is sufficiently small (10^(-5) for instance) or if the iterations are over. If significantly smaller stop criteria is chosen and number of iterations are increased, the error will be lessen and the predicted x will be closer to the actual result.

For the Krylov subspace the algorithm generates ![image](https://user-images.githubusercontent.com/36455629/123247825-f05f9400-d4ef-11eb-9f59-3d074481cd31.png)
 vectors in the lecture notes. ![image](https://user-images.githubusercontent.com/36455629/123247679-d1610200-d4ef-11eb-8f3c-c798bf829324.png)
 thus, the bases are orthogonal.

![image](https://user-images.githubusercontent.com/36455629/123247205-4f70d900-d4ef-11eb-88d9-3420ed2fac95.png)

![image](https://user-images.githubusercontent.com/36455629/123247230-5566ba00-d4ef-11eb-92fa-d8810e63d654.png)

![image](https://user-images.githubusercontent.com/36455629/123247256-5d265e80-d4ef-11eb-9430-75e5198cf6f7.png)

By looking Figures 5 and 6, it can be said that smaller noise result as less error for both ![image](https://user-images.githubusercontent.com/36455629/123247642-c7d79a00-d4ef-11eb-8d7b-7afcb03aaf03.png). Moreover, as the size of A increases, both errors decreases and the computed x gets closer to the actual x, which are the same results and observations of the previous method. Since similar stop criteria are used for GMRES and CG methods their plots are resulted as expected. The impact of ![image](https://user-images.githubusercontent.com/36455629/123247513-a70f4480-d4ef-11eb-9d70-05860f857799.png)
 cannot be seen. GMRES and CG methods gave similar error results for both ![image](https://user-images.githubusercontent.com/36455629/123247477-9f4fa000-d4ef-11eb-8598-c0cd158bc531.png).

# Pseudo-inverse Method

Since this method uses every data in the input matrix A, it takes the longest time so the maximum size A is chosen to be 1000.

![image](https://user-images.githubusercontent.com/36455629/123247340-792a0000-d4ef-11eb-9220-09635d923c79.png)

Results are similar as the previous two results, so it can be said that similar performance values are obtained but GMRES and CG uses less data so their computation cost is less. Computation time for GMRES and CG are approximately 0.015 seconds for input size 500x500 and 0.15 for Pseudo-inverse method. For larger data pseudo-inverse method is problematic because of its time cost.

