# CSE701-Project-Work-

## Summary
This programme implements the EM algorithm and demonstrates how to use the algorithm for clustering problems

## Author
This code was written by Michael Agronah(agronahm@mcmaster.ca) for a course project in Foundations of Modern Scientific Programming (CSE 701) taken at McMaster University.

## The Expectation Maximization Algorithm
The Expectation Maximization Algorithm is an iterative method for estimating parameters of a Gaussian distribution that maximises likelihoods [1]. Each iteration involves two steps: an expectation process (the E- step) and a maximisation step(M-step)

Given data X, and the 


The implementation follows the process: 
1. Estimate likelihoods of each data point belonging to each cluster 
2. Estimate posteriors of each cluster for each data point 
3. Update the mean and variances
5. Compute the weights of each cluster 
6. Compute the sum of likelihood and compare with a threshold value

The process in repeated iteratedly untill the sum of likelihood is threshood contidion is satified. 

Often a minisiting the sum of likelihood instead of the maximisation is better so we can minimise the negative likelihood. 

#Example for 2 Cluster analysis 
Intial mean x1  x2  
Intial variance v1  v2

Compute likelihoods
Class A: P(x_i|a)=Gausian(x1, v2) 
Class B: P(x_i|b) = Gausian(x2, v2) 

Compute posterior
a_i = P(a|x_i) = P(x_i|a)*prior /normaliser
b_i = P(b|x_i) = P(x_i|b)*prior /normaliser

Update mean 
mean_a = sum (a_ix_i)/(sum(a_i))
mean_b = sum (b_ix_i)/(sum(b_i))

Update variance
variance_a = sum(a_i(x_i - mean_a)^2)/ sum(a_i)
variance_b = sum(b_i(x_i - mean_b)^2)/ sum(b_i)

LIKE = sum (weight_A*likelihoods_A + weight_B*likelihoods_B) 

 

## Usage
The header file (project.hpp) contains classes for writing the data used for the implementation, classes to throw execptions and a class for the EM algorithm. 
The project.cpp file contains defined functions and vector operator overloads.
The main.cpp file test the code and returns the outputs 

## Example
#Input
Initial means  = {2000.0, 5000.0, 7000.0, 4000.0};
Initial variances  = {130.0, 700.0, 100.0, 600.0};

#Output
Results from 2 clusters
First_Cluster: Mean =8767.28 Variance = 2.47081e+06
Second_Cluster:Mean = 8153.19 Variance = 2.3604e+06
Minimised Likelihood  =-0.670399
BIC selection criteria =36.8337

Results from 3 clusters
First_Cluster:  Mean =8581.72 Variance = 4.58361e+07
Second_Cluster: Mean = 8024.27 Variance = 1.13866e+07
Third_Cluster: Mean = 8705.04 Variance = 6.65749e+07
Minimum Likelihood  =-0.290342
BIC selection criteria =48.612

Results from 4 clusters
First_Cluster:  Mean =8649.71 Variance = 4.67268e+07
Second_Cluster: Mean = 8052.96   Variance = 1.15905e+07
Third_Cluster: Mean = 8771.25   Variance = 6.76143e+07
Fourth_Cluster: Mean = 8345.14 Variance = 2.13564e+07
Minimum Likelihood  =-0.22545
BIC selection criteria =58.5793

Grouping into 2 clusters gives the best result

##References
1. Fessler, J. A., & Hero, A. O. (1994). Space-alternating generalized expectation-maximization algorithm. IEEE Transactions on signal processing, 42(10), 2664-2677.

2.Dellaert, F. (2002). The expectation maximization algorithm. Georgia Institute of Technology.
