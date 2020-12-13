/**
 * @file project.cpp
 * @author Michael Agronah (agronahm@mcmaster.ca)
 * @version 0.1
 * @date 2020-12-13
 * @copyright Copyright (c) 2020
 * @brief File with functions and vector operation overloads
 */

#define _USE_MATH_DEFINES
#include "project.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
//#include <algorithm>

using namespace std;

/**
 * @brief  template for overloading the exponential function
 * 
 * @param b vector of type T
 * @return template <typename T> 
 */

template <typename T>
vector<T> exp(const vector<T> &b)
{
    vector<T> c;
    size_t s{b.size()};

    for (size_t i{0}; i < s; i++)
    {
        c.push_back(exp(b[i]));
    }
    return c;
}

/**
 * @brief operator overloads for << 
 * 
 * @param out 
 * @param v 
 * @return ostream& 
 */

ostream &operator<<(ostream &out, const vector<double> &v)
{
    size_t s{v.size() - 1};
    out << "(";
    for (size_t i{0}; i < s; i++)
        out << v[i] << ", ";
    out << v[s] << ")\n";
    return out;
}

vector<double> operator+(const vector<double> &v, const vector<double> &w)
{
    size_t s{v.size()};
    if (s != w.size())
        throw size_must_match{};
    vector<double> u(s);
    for (size_t i{0}; i < s; i++)
        u[i] = v[i] + w[i];
    return u;
}

vector<double> operator+=(vector<double> &v, const vector<double> &w)
{
    v = v + w;
    return v;
}

vector<double> operator-(const vector<double> &v)
{
    size_t s{v.size()};
    vector<double> u(s);
    for (size_t i{0}; i < s; i++)
        u[i] = -v[i];
    return u;
}

vector<double> operator-(const vector<double> &v, const vector<double> &w)
{
    size_t s{v.size()};
    if (s != w.size())
        throw size_must_match{};
    vector<double> u(s);
    for (size_t i{0}; i < s; i++)
        u[i] = v[i] - w[i];
    return u;
}

vector<double> operator-=(vector<double> &v, const vector<double> &w)
{
    v = v - w;
    return v;
}

double operator*(const vector<double> &v, const vector<double> &w)
{
    size_t s{v.size()};
    if (s != w.size())
        throw size_must_match{};
    double p{0};
    for (size_t i{0}; i < s; i++)
        p += v[i] * w[i];
    return p;
}

vector<double> operator*(const double &x, const vector<double> &v)
{
    size_t s{v.size()};
    vector<double> u(s);
    for (size_t i{0}; i < s; i++)
        u[i] = x * v[i];
    return u;
}

vector<double> operator*(const vector<double> &v, const double &x)
{
    return x * v;
}

/**
 * @brief Function to write data
 * 
 * @param dataSize_perGaussian 
 */

void WriteData_File::WriteData(const int &dataSize_perGaussian)
{
    ofstream outFile;
    string outputFilename = "Data.txt";
    outFile.open("Data.txt");
    {
        random_device rd;
        mt19937 mt(rd());
        normal_distribution<double> dis1(data_mean[0], data_standard_deviation[0]);
        normal_distribution<double> dis2(data_mean[1], data_standard_deviation[1]);

        if (outFile.is_open())
        {
            for (int i = 0; i < dataSize_perGaussian; i++)
            {
                outFile << dis1(mt) << " "
                        << " " << dis2(mt) << endl;
            }
            outFile.close();
        }
        else
        {
            throw could_not_open_Datafile{};
        }
    }
}
/**
 * @brief Function to read data
 * 
 * @return vector<double> 
 */
vector<double> Max_Expectation::data()
{
    vector<double> dat;
    ifstream inFile;
    string inFilename = "Data.txt";
    inFile.open("Data.txt");
    if (inFile.is_open())
    {
        double val;
        while (inFile >> val)
        {
            dat.push_back(val);
        }
        inFile.close();
    }
    else
    {
        throw could_not_open_Datafile{};
    }

    return dat;
}

/**
 * @brief Function to compute unmormalised posteriors
 * 
 * @param mu  means of individual classes
 * @param var variances of individual classes
 * @return vector<double> 
 */
vector<double> Max_Expectation::Unnormalised_posteriors(const vector<double> &mu, const vector<double> &var)
{
    size_t s{var.size()};
    for (size_t i{0}; i < s; i++)
    {
        if (var[i] < 0)
        {
            throw variance_must_be_positive{};
        }
    }

    vector<double> dat{data()};
    vector<double> unnorm_poster{0.0};
    const double prior{0.5};

    vector<double> Likelihoods{0.0};
    double normalizer{0.0};
    double post{0.0};

    size_t k{dat.size()};
    // Computing the likelihoods of data points belonging to a given cluster
    for (size_t j{0}; j < s; j++)
    {
        for (size_t i{0}; i < k; i++)
        {
            normalizer = 1.0 / sqrt(2.0 * M_PI * var[j]);
            post = -(0.5 / var[j]) * pow(dat[i] - mu[j], 2);
            Likelihoods.push_back(log(normalizer) + post);
        }
    }

    // Computing posteriors for each data point
    unnorm_poster = prior * Likelihoods;
    return unnorm_poster;
}

/**
 * @brief Function that implements the EM algorithm for 2 clusters
 * 
 * @param m initial class means
 * @param v initial class variances
 * @return vector<double> optimal means and variances 
 */

vector<double> Max_Expectation::em(const double (&m)[2], const double (&v)[2])
{
    vector<double> x_data{data()};

    //Maximuim iteration and Threshold
    int max_ter = Max_iterations;
    double thresh = threshold;

    // Iteration Counts
    int iterCnt = 0;
    double sum_Likelihood{0.1};

    //Intitialising curent means and variances
    vector<double> current_mean{m[0], m[1]};
    vector<double> current_var{v[0], v[1]};

    //Prior and normaliser for posterior
    double prior = prior_2clusters;

    vector<double> Normaliser;
    vector<double> ClassA_unnormalised;
    vector<double> ClassB_unnormalised;

    vector<double> ClassA_posterior;
    vector<double> ClassB_posterior;

    vector<double> unnorm_varA;
    vector<double> unnorm_varB;

    size_t s{x_data.size()};

    // Implementing expectation maximisation
    while ((iterCnt < max_ter) && (sum_Likelihood > thresh))
    {
        vector<double> unnorm_post{Unnormalised_posteriors(current_mean, current_var)};

        for (size_t i{0}; i < s; i++)
        {
            ClassA_unnormalised.push_back(unnorm_post[i]);
            ClassB_unnormalised.push_back(unnorm_post[i + s]);
        }

        // Compute Normaliser for posterior
        Normaliser = ClassA_unnormalised + ClassB_unnormalised;

        // Compute posteriors for each data point
        for (size_t i{0}; i < s; i++)
        {
            double p = ClassA_unnormalised[i] / Normaliser[i];
            ClassA_posterior.push_back(p);
            ClassB_posterior.push_back(1 - p);
        }

        //normalisers for updating means and variances
        double sum_A = accumulate(ClassA_posterior.begin(), ClassA_posterior.end(), 0.0);
        double sum_B = accumulate(ClassB_posterior.begin(), ClassB_posterior.end(), 0.0);

        // Update mean
        current_mean[0] = (x_data * ClassA_posterior) / sum_A;
        current_mean[1] = (x_data * ClassB_posterior) / sum_B;

        // Update variances
        for (size_t i{0}; i < s; i++)
        {
            unnorm_varA.push_back(pow((x_data[i] - current_mean[0]), 2));
            unnorm_varB.push_back(pow((x_data[i] - current_mean[1]), 2));
        }
        current_var[0] = (ClassA_posterior * unnorm_varA) / sum_A;
        current_var[1] = (ClassB_posterior * unnorm_varB) / sum_B;

        // Computing likelihoods for each cluster
        vector<double> pA = (1.0 / prior) * ClassA_unnormalised;
        vector<double> pB = (1.0 / prior) * ClassB_unnormalised;

        // Computing weights for each class
        double w_A = (sum_A / s);
        double w_B = (sum_B / s);

        // Computing total likelihoods
        vector<double> like{w_A * exp(pA) + w_B * exp(pB)};
        sum_Likelihood = -1.0 * accumulate(like.begin(), like.end(), 0.0);

        // Delete entries in the vectors for reintialisation
        unnorm_post.clear();
        ClassA_unnormalised.clear();
        ClassB_unnormalised.clear();
        ClassA_posterior.clear();
        ClassB_posterior.clear();
        unnorm_varA.clear();
        unnorm_varB.clear();
        like.clear();

        iterCnt++;
    }

    //Compute the likelihoods at optimal parameter values
    double norm{0.0};
    double m_post{0.0};
    double m_likel{0.0};
    double max_likelihood{0.0};
    size_t k{current_mean.size()};

    for (size_t j{0}; j < k; j++)
    {
        for (size_t i{0}; i < s; i++)
        {
            norm = 1.0 / sqrt(2.0 * M_PI * current_var[j]);
            m_post = exp(-(0.5 / current_var[j]) * pow(x_data[i] - current_mean[j], 2));
        }
        m_likel = norm * m_post;
        max_likelihood += m_likel;
    }

    // Compute the BIC value for model selection criteria
    double BIC = 4 * log(s) - 2 * log(max_likelihood);
    return {current_mean[0], current_mean[1], current_var[0], current_var[1], BIC, sum_Likelihood};
}
/**
 * @brief Function that implements the EM algorithm for 3 clusters
 * 
 * @param m initial class means
 * @param v initial class variances
 * @return vector<double> optimal means and variances 
 */

vector<double> Max_Expectation::em(const double (&m)[3], const double (&v)[3])
{
    int max_ter = Max_iterations;
    int iterCnt = 0;
    double sum_Likelihood_3clust{0.1};
    double thresh = threshold;

    //intitailising means and variances for EM algorithm
    vector<double> current_mean_3clust{m[0], m[1], m[2]};
    vector<double> current_var_3clust{v[0], v[1], v[2]};

    vector<double> x_data{data()};
    double prior = prior_3clusters;

    vector<double> Normaliser_3clust;

    vector<double> ClassA_unnormalised_3clust;
    vector<double> ClassB_unnormalised_3clust;
    vector<double> ClassC_unnormalised_3clust;

    vector<double> ClassA_posterior_3clust;
    vector<double> ClassB_posterior_3clust;
    vector<double> ClassC_posterior_3clust;

    vector<double> unnorm_varA_3clust;
    vector<double> unnorm_varB_3clust;
    vector<double> unnorm_varC_3clust;

    //normalisers for updating means and variances
    double sum_A_3clust{0.0};
    double sum_B_3clust{0.0};
    double sum_C_3clust{0.0};

    size_t s{x_data.size()};
    while ((iterCnt < max_ter) && (sum_Likelihood_3clust > thresh))
    {
        vector<double> unnorm_post_3clust{Unnormalised_posteriors(current_mean_3clust, current_var_3clust)};

        for (size_t i{0}; i < s; i++)
        {
            ClassA_unnormalised_3clust.push_back(unnorm_post_3clust[i]);
            ClassB_unnormalised_3clust.push_back(unnorm_post_3clust[i + s]);
            ClassC_unnormalised_3clust.push_back(unnorm_post_3clust[i + 2 * s]);
        }

        // Compute Normaliser
        Normaliser_3clust = ClassA_unnormalised_3clust + ClassB_unnormalised_3clust + ClassC_unnormalised_3clust;

        for (size_t i{0}; i < s; i++)
        {
            double p_classA = ClassA_unnormalised_3clust[i] / Normaliser_3clust[i];
            double p_classB = ClassB_unnormalised_3clust[i] / Normaliser_3clust[i];

            ClassA_posterior_3clust.push_back(p_classA);
            ClassB_posterior_3clust.push_back(p_classB);
            ClassC_posterior_3clust.push_back(1 - p_classA - p_classB);
        }

        //Compute the normalisers for updating means and variances
        sum_A_3clust = accumulate(ClassA_posterior_3clust.begin(), ClassA_posterior_3clust.end(), 0.0);
        sum_B_3clust = accumulate(ClassB_posterior_3clust.begin(), ClassB_posterior_3clust.end(), 0.0);
        sum_C_3clust = accumulate(ClassC_posterior_3clust.begin(), ClassC_posterior_3clust.end(), 0.0);

        for (size_t i{0}; i < s; i++)
        {
            unnorm_varA_3clust.push_back(pow((x_data[i] - current_mean_3clust[0]), 2));
            unnorm_varB_3clust.push_back(pow((x_data[i] - current_mean_3clust[1]), 2));
            unnorm_varC_3clust.push_back(pow((x_data[i] - current_mean_3clust[2]), 2));
        }

        // Update mean
        current_mean_3clust[0] = (x_data * ClassA_posterior_3clust) / sum_A_3clust;
        current_mean_3clust[1] = (x_data * ClassB_posterior_3clust) / sum_B_3clust;
        current_mean_3clust[2] = (x_data * ClassC_posterior_3clust) / sum_C_3clust;

        // Update variances
        current_var_3clust[0] = (ClassA_posterior_3clust * unnorm_varA_3clust) / sum_A_3clust;
        current_var_3clust[1] = (ClassB_posterior_3clust * unnorm_varB_3clust) / sum_B_3clust;
        current_var_3clust[2] = (ClassC_posterior_3clust * unnorm_varC_3clust) / sum_C_3clust;

        // Compute weights for each class
        vector<double> pA_3clust = (1.0 / prior) * ClassA_unnormalised_3clust;
        vector<double> pB_3clust = (1.0 / prior) * ClassB_unnormalised_3clust;
        vector<double> pC_3clust = (1.0 / prior) * ClassC_unnormalised_3clust;

        // Computing total likelihoods
        vector<double> like_3clust{(sum_A_3clust / s) * exp(pA_3clust) + (sum_B_3clust / s) * exp(pB_3clust) + (sum_C_3clust / s) * exp(pC_3clust)};
        sum_Likelihood_3clust = -1.0 * accumulate(like_3clust.begin(), like_3clust.end(), 0.0);

        // Delete entries in the vectors for reintialisation
        unnorm_post_3clust.clear();
        ClassA_unnormalised_3clust.clear();
        ClassB_unnormalised_3clust.clear();
        ClassC_unnormalised_3clust.clear();

        ClassA_posterior_3clust.clear();
        ClassB_posterior_3clust.clear();
        ClassC_posterior_3clust.clear();

        unnorm_varA_3clust.clear();
        unnorm_varB_3clust.clear();
        unnorm_varC_3clust.clear();
        like_3clust.clear();

        iterCnt++;
    }

    //Compute the likelihoods at optimal parameter values
    double norm_3clust{0.0};
    double m_post_3clust{0.0};
    double m_likel_3clust{0.0};
    double max_likelihood_3clust{0.0};

    size_t k{current_mean_3clust.size()};
    for (size_t j{0}; j < k; j++)
    {
        for (size_t i{0}; i < s; i++)
        {
            norm_3clust = 1.0 / sqrt(2.0 * M_PI * current_var_3clust[j]);
            m_post_3clust = exp(-(0.5 / current_var_3clust[j]) * pow(x_data[i] - current_mean_3clust[j], 2));
        }
        m_likel_3clust = norm_3clust * m_post_3clust;
        max_likelihood_3clust += m_likel_3clust;
    }

    // Compute the BIC value for model selection criteria
    double BIC_3clust = 6 * log(s) - 2 * log(max_likelihood_3clust);
    return {current_mean_3clust[0], current_mean_3clust[1], current_mean_3clust[2], current_var_3clust[0], current_var_3clust[1], current_var_3clust[2], BIC_3clust, sum_Likelihood_3clust};
}

/**
 * @brief Function that implements the EM algorithm for 4 clusters
 * 
 * @param m initial class means
 * @param v initial class variances
 * @return vector<double> optimal means and variances 
 */

vector<double> Max_Expectation::em(const double (&m)[4], const double (&v)[4])
{
    int max_ter = Max_iterations;
    int iterCnt = 0;
    double sum_Likelihood_4clust{0.1};
    double thresh = threshold;
    vector<double> current_mean_4clust{m[0], m[1], m[2], m[3]};
    vector<double> current_var_4clust{v[0], v[1], v[2], v[3]};

    vector<double> x_data{data()};
    double prior = prior_4clusters;

    // A vectors for the normalising constants
    vector<double> Normaliser_4clust;

    // Vectors of the unnormalised posteriors for each class
    vector<double> ClassA_unnormalised_4clust;
    vector<double> ClassB_unnormalised_4clust;
    vector<double> ClassC_unnormalised_4clust;
    vector<double> ClassD_unnormalised_4clust;

    // Vectors of posteriors for each class
    vector<double> ClassA_posterior_4clust;
    vector<double> ClassB_posterior_4clust;
    vector<double> ClassC_posterior_4clust;
    vector<double> ClassD_posterior_4clust;

    // Vectors of unnormalised variances for each class
    vector<double> unnorm_varA_4clust;
    vector<double> unnorm_varB_4clust;
    vector<double> unnorm_varC_4clust;
    vector<double> unnorm_varD_4clust;

    //Normalisers for each class
    double sum_A_4clust{0.0};
    double sum_B_4clust{0.0};
    double sum_C_4clust{0.0};
    double sum_D_4clust{0.0};

    size_t s{x_data.size()};
    while ((iterCnt < max_ter) && (sum_Likelihood_4clust > thresh))
    {
        vector<double> unnorm_post_4clust{Unnormalised_posteriors(current_mean_4clust, current_var_4clust)};

        for (size_t i{0}; i < s; i++)
        {
            ClassA_unnormalised_4clust.push_back(unnorm_post_4clust[i]);
            ClassB_unnormalised_4clust.push_back(unnorm_post_4clust[i + s]);
            ClassC_unnormalised_4clust.push_back(unnorm_post_4clust[i + 2 * s]);
            ClassD_unnormalised_4clust.push_back(unnorm_post_4clust[i + 3 * s]);
        }

        // Compute vector of posterior Normalisers
        Normaliser_4clust = ClassA_unnormalised_4clust + ClassB_unnormalised_4clust + ClassC_unnormalised_4clust + ClassD_unnormalised_4clust;

        // Compute the posterior for each data point

        for (size_t i{0}; i < s; i++)
        {
            double p_classA = ClassA_unnormalised_4clust[i] / Normaliser_4clust[i];
            double p_classB = ClassB_unnormalised_4clust[i] / Normaliser_4clust[i];
            double p_classC = ClassC_unnormalised_4clust[i] / Normaliser_4clust[i];

            ClassA_posterior_4clust.push_back(p_classA);
            ClassB_posterior_4clust.push_back(p_classB);
            ClassC_posterior_4clust.push_back(p_classC);
            ClassD_posterior_4clust.push_back(1 - (p_classA + p_classB + p_classC));
        }

        // Compute Normalisers for updating variance and means of each class
        sum_A_4clust = accumulate(ClassA_posterior_4clust.begin(), ClassA_posterior_4clust.end(), 0.0);
        sum_B_4clust = accumulate(ClassB_posterior_4clust.begin(), ClassB_posterior_4clust.end(), 0.0);
        sum_C_4clust = accumulate(ClassC_posterior_4clust.begin(), ClassC_posterior_4clust.end(), 0.0);
        sum_D_4clust = accumulate(ClassD_posterior_4clust.begin(), ClassD_posterior_4clust.end(), 0.0);

        for (size_t i{0}; i < s; i++)
        {
            unnorm_varA_4clust.push_back(pow((x_data[i] - current_mean_4clust[0]), 2));
            unnorm_varB_4clust.push_back(pow((x_data[i] - current_mean_4clust[1]), 2));
            unnorm_varC_4clust.push_back(pow((x_data[i] - current_mean_4clust[2]), 2));
            unnorm_varD_4clust.push_back(pow((x_data[i] - current_mean_4clust[3]), 2));
        }

        // Update the means of each class
        current_mean_4clust[0] = (x_data * ClassA_posterior_4clust) / sum_A_4clust;
        current_mean_4clust[1] = (x_data * ClassB_posterior_4clust) / sum_B_4clust;
        current_mean_4clust[2] = (x_data * ClassC_posterior_4clust) / sum_C_4clust;
        current_mean_4clust[3] = (x_data * ClassD_posterior_4clust) / sum_D_4clust;

        // Updating the variances of each class
        current_var_4clust[0] = (ClassA_posterior_4clust * unnorm_varA_4clust) / sum_A_4clust;
        current_var_4clust[1] = (ClassB_posterior_4clust * unnorm_varB_4clust) / sum_B_4clust;
        current_var_4clust[2] = (ClassC_posterior_4clust * unnorm_varC_4clust) / sum_C_4clust;
        current_var_4clust[3] = (ClassD_posterior_4clust * unnorm_varD_4clust) / sum_D_4clust;

        //
        vector<double> pA_4clust = (1.0 / prior) * ClassA_unnormalised_4clust;
        vector<double> pB_4clust = (1.0 / prior) * ClassB_unnormalised_4clust;
        vector<double> pC_4clust = (1.0 / prior) * ClassC_unnormalised_4clust;
        vector<double> pD_4clust = (1.0 / prior) * ClassD_unnormalised_4clust;

        // Computing total likelihoods
        vector<double> like_4clust{(sum_A_4clust / s) * exp(pA_4clust) + (sum_B_4clust / s) * exp(pB_4clust) + (sum_C_4clust / s) * exp(pC_4clust) + (sum_D_4clust / s) * exp(pD_4clust)};
        sum_Likelihood_4clust = -1.0 * accumulate(like_4clust.begin(), like_4clust.end(), 0.0);

        // Delete entries in the vectors for reintialisation
        unnorm_post_4clust.clear();
        ClassA_unnormalised_4clust.clear();
        ClassB_unnormalised_4clust.clear();
        ClassC_unnormalised_4clust.clear();
        ClassD_unnormalised_4clust.clear();

        ClassA_posterior_4clust.clear();
        ClassB_posterior_4clust.clear();
        ClassC_posterior_4clust.clear();
        ClassD_posterior_4clust.clear();

        unnorm_varA_4clust.clear();
        unnorm_varB_4clust.clear();
        unnorm_varC_4clust.clear();
        unnorm_varD_4clust.clear();
        like_4clust.clear();

        iterCnt++;
    }

    //Compute the likelihoods at optimal parameter values
    double norm_4clust{0.0};
    double m_post_4clust{0.0};
    double m_likel_4clust{0.0};
    double max_likelihood_4clust{0.0};

    size_t k{current_mean_4clust.size()};
    for (size_t j{0}; j < k; j++)
    {
        for (size_t i{0}; i < s; i++)
        {
            norm_4clust = 1.0 / sqrt(2.0 * M_PI * current_var_4clust[j]);
            m_post_4clust = exp(-(0.5 / current_var_4clust[j]) * pow(x_data[i] - current_mean_4clust[j], 2));
        }
        m_likel_4clust = norm_4clust * m_post_4clust;
        max_likelihood_4clust += m_likel_4clust;
    }

    // Compute the BIC value for model selection criteria
    double BIC_4clust = 8 * log(s) - 2 * log(max_likelihood_4clust);
    return {current_mean_4clust[0], current_mean_4clust[1], current_mean_4clust[2], current_mean_4clust[3], current_var_4clust[0], current_var_4clust[1], current_var_4clust[2], current_var_4clust[3], BIC_4clust, sum_Likelihood_4clust};
}

/**
 * @brief Function takes results of estimates from the 2 cluster, 3 cluster and 4 cluster analysis
 * 
 * @param mu initial class means
 * @param var initial class variances
 * @return vector<vector<double>> optimal means and variances for 2, 3 and 4 clusters 
 */
vector<vector<double>> Max_Expectation::Exp_Max(vector<double> &mu, vector<double> &var)
{
    double mu_2clust[2] = {mu[0], mu[1]};
    double var_2clust[2] = {var[0], var[1]};

    double mu_3clust[3] = {mu[0], mu[1], mu[2]};
    double var_3clust[3] = {var[0], var[1], var[2]};

    double mu_4clust[4] = {mu[0], mu[1], mu[2], mu[3]};
    double var_4clust[4] = {var[0], var[1], var[2], var[3]};

    vector<double> results_2clust{em(mu_2clust, var_2clust)};
    vector<double> results_3clust{em(mu_3clust, var_3clust)};
    vector<double> results_4clust{em(mu_4clust, var_4clust)};

    return {results_2clust, results_3clust, results_4clust};
}
