#pragma once
/**
 * @file project.hpp
 * @author Michael Agronah (agronahm@mcmaster.ca)
 * @version 0.1
 * @date 2020-12-13
 * @copyright Copyright (c) 2020
 *  
 *@brief Header file with classes for implementing EM algorithm
 *@details File contains classes to throw exceptions, class to write data and class containing functions and variables for EM algorithm
 */

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;
/**
 * @brief Class to throw exception if data file does not open is negative 
 */
class could_not_open_Datafile
{
};

/**
 * @brief Class with function and parameters for writing the data 
 * 
 */

class WriteData_File
{
public:
    //Function to write data into .txt file
    void WriteData(const int &);

    //Means and variances of gausian distributions to simulate data from
    vector<double> data_mean{10000, 7000};
    vector<double> data_standard_deviation{500, 600};

private:
};

/**
 * @brief Class with parameters and functions for the EM algorithm
 * 
 */

class Max_Expectation
{
public:
    //Fucntion to store data into a vector
    vector<double> data();

    //Fucntion to store data into a vector
    vector<double> Unnormalised_posteriors(const vector<double> &mu, const vector<double> &var);

    //Fucntion overloads for estimating the variances and means
    //2 clusters, 3 clusters and 4 clusters

    vector<double> em(const double (&m)[2], const double (&v)[2]);
    vector<double> em(const double (&m)[3], const double (&v)[3]);
    vector<double> em(const double (&m)[4], const double (&v)[4]);

    //Function to combine results from the overloaded em function
    vector<vector<double>> Exp_Max(vector<double> &mu, vector<double> &var);

    // Exception to be thrown if variance is negative
    class variance_must_be_positive
    {
    };

private:
    //Threshold and  Maximum iterations
    double threshold{1e-8};
    int Max_iterations{1000};

    //Priors for 2 cluster, 3 cluster and 4 cluster analysis
    double prior_2clusters{0.5};
    double prior_3clusters{1.0 / 3.0};
    double prior_4clusters{1.0 / 4.0};
};

/**
 * @brief Class to throw exception f vector sizes differ
 * 
 */
class size_must_match
{
};

/**
 * @brief Operetaor overloads for vectors
 * 
 */
ostream &operator<<(ostream &, const vector<double> &);

vector<double> operator+(const vector<double> &, const vector<double> &);

vector<double> operator+=(vector<double> &, const vector<double> &);

vector<double> operator-(const vector<double> &);

vector<double> operator-(const vector<double> &, const vector<double> &);

vector<double> operator-=(vector<double> &, const vector<double> &);

double operator*(const vector<double> &, const vector<double> &);

vector<double> operator*(const double &, const vector<double> &);

vector<double> operator*(const vector<double> &, const double &);
