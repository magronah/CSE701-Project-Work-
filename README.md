# CSE701-Project-Work-

## Summary
This programme implements the EM algorithm and demonstrates how to use the algorithm for clustering problems

## Author
This code was written by Michael Agronah(agronahm@mcmaster.ca) for a course project in Foundations of Modern Scientific Programming (CSE 701) taken at McMaster University.

## Usage
The header file (project.hpp) contains classes for writing the data used for the implementation, classes to throw execptions and a class for the EM algorithm
The project.cpp file contains defined functions and vector operator overloads
The main.cpp file test the code and returns the outputs 

The code requires intial inputs of means and variances for 4 clusters and returns the following 
1. Optimal mean and variance estimates for 2, 3 and 4 clusters, 
2. The BIC values of 2, 3 and 4 clusters and
3. The minimal likelihoods for 2, 3 and 4 clusters. 
