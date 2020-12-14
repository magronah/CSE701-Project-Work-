# CSE701-Project-Work-

## Summary
This programme implements the EM algorithm and demonstrates how to use the algorithm for clustering problems

## Author
This code was written by Michael Agronah(agronahm@mcmaster.ca) for a course project in Foundations of Modern Scientific Programming (CSE 701) taken at McMaster University.

## The Expectation Maximization Algorithm
The Expectation Maximization Algorithm is an iterative method for estimating parameters of a Gaussian distribution that maximises likelihoods [1]. Each iteration involves two steps: the Expectation (or the E-) step and the Maximisation(M-) step. 

Given observations X= {x_i, x_2, ..., x_n} and k_1,k_2,...,k_m  number of cluster. 

The E-step finds the posterior probability of each sample x_i, belonging to each class k_i. The M-step updates the model parameters and estimates the weights for the individual clusters [2]. 

The implementation is summarised as follows: Given intial values for parameters, for each x_i,

#### 1.compute the likelihood of x_i belonging to each cluster

#### 2.estimate posterior probability of each cluster given x_i 

#### 3.update paramters of the gaussian distribution

#### 4.estimate weights of each cluster 

#### 5.compute total likelihood and compare with a threshold condition

## Demonstration for 1d dataset and 2 clusters 
For a univariant gaussian, the parameters to estimate are the means and the variances. 

Given vectors of intial mean = [x_1 , x_2]  and intial variance = [v_1 , v-2],

##### 1.compute likelihoods of cluster A (P(x_i|a)) and  likelihoods of cluster B ( P(x_i|b)) from a Gaussian distributions with parameters (x_1,v_1) and (x_2, v_2) respectively.

##### 2.compute posterior probabilities of clusters A and B

cluster A: a_i= P(a|x_i) = P(x_i|a)*prior /normaliser

cluster B: b_i= P(b|x_i) = P(x_i|b)*prior /normaliser

##### 3.Update mean 

mean_a = sum (a_i*x_i)/(sum(P(a|x_i))), i= 1,...,n

mean_b = sum (b_i*x_i)/(sum(P(b|x_i))), i= 1,...,n

##### 4. Update variance

variance_a = sum(a_i(x_i - mean_a)^2)/ sum(a_i)

variance_b = sum(b_i(x_i - mean_b)^2)/ sum(b_i)

##### 5.Compute weights 

Cluster A: sum a_i, i= 1,...,n

Cluster B: sum b_i,  i= 1,...,n

#### 6.Compute total likelihoods = sum (weight_A*P(x_i|a) + weight_B*P(x_i|b)),  i= 1,...,n

## Usage
The header file (project.hpp) contains classes for writing the data used for the implementation, classes to throw execptions and a class for the EM algorithm. 
The project.cpp file contains defined functions and vector operator overloads.
The main.cpp file test the code and returns the outputs 

## Example
### Inputs
Initial means  = {2000.0, 5000.0, 7000.0, 4000.0};
Initial variances  = {130.0, 700.0, 100.0, 600.0};

### Output
#### Results from 2 clusters

First_Cluster: Mean =8767.28 Variance = 2.47081e+06

Second_Cluster:Mean = 8153.19 Variance = 2.3604e+06

Minimised Likelihood  =-0.670399

BIC selection criteria =36.8337


#### Results from 3 clusters

First_Cluster:  Mean =8581.72 Variance = 4.58361e+07

Second_Cluster: Mean = 8024.27 Variance = 1.13866e+07

Third_Cluster: Mean = 8705.04 Variance = 6.65749e+07

Minimum Likelihood  =-0.290342

BIC selection criteria =48.612


#### Results from 4 clusters

First_Cluster:  Mean =8649.71 Variance = 4.67268e+07

Second_Cluster: Mean = 8052.96   Variance = 1.15905e+07

Third_Cluster: Mean = 8771.25   Variance = 6.76143e+07

Fourth_Cluster: Mean = 8345.14 Variance = 2.13564e+07

Minimum Likelihood  =-0.22545

BIC selection criteria =58.5793

#### Optimal number of cluster
Grouping into 2 clusters gives the best result

## References
1. Fessler, J. A., & Hero, A. O. (1994). Space-alternating generalized expectation-maximization algorithm. IEEE Transactions on signal processing, 42(10), 2664-2677.

2. Dellaert, F. (2002). The expectation maximization algorithm. Georgia Institute of Technology.
