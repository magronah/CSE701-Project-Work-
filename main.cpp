/**
 * @file main.cpp
 * @author Michael Agronah (agronahm@mcmaster.ca)
 * @brief A program to cluster data into classes
 * @version 0.1
 * @date 2020-12-12
 * 
 * @copyright Copyright (c) 2020
 * 
 */

#include "project.hpp"
#include <iostream>
#include <vector>

using namespace std;
/**
 * @mainpage This is a program demonstrating how to implement the EM algorithm for clustering  problem. The program was writen as a project for the course (CSE 701) Foundations of Modern Scientific Programming, taken at McMaster Univertsity <span style="color: green;">HTML tags.</span>. 
 */
int main()
{
     try
     {
          WriteData_File wr;

          int dataSize_perGaussian = 100;
          wr.WriteData(dataSize_perGaussian);

          Max_Expectation mn;
          vector<double> var = {130.0, 70.0, 100.0, 60.0};
          vector<double> mu = {2000.0, 5000.0, 700.0, 4000.0};

          vector<vector<double>> res{mn.Exp_Max(mu, var)};
          vector<double> res_2clusters{res[0]};
          vector<double> res_3clusters{res[1]};
          vector<double> res_4clusters{res[2]};

          cout << "  " << endl;
          cout << "Results from 2 clusters" << endl;
          cout << "  " << endl;

          cout << "First_Cluster: "
               << "Mean =" << res_2clusters[0] << " "
               << "Variance = " << res_2clusters[2] << endl;

          cout << "Second_Cluster:"
               << "Mean = " << res_2clusters[1] << " "
               << "Variance = " << res_2clusters[3] << endl;

          cout << "Minimised Likelihood  =" << res_2clusters[5] << endl;
          cout << "BIC selection criteria =" << res_2clusters[4] << endl;

          cout << "  " << endl;
          cout << "  " << endl;

          cout << "Results from 3 clusters" << endl;
          cout << "  " << endl;

          cout << "First_Cluster:  "
               << "Mean =" << res_3clusters[0] << " "
               << "Variance = " << res_3clusters[3] << endl;

          cout << "Second_Cluster: "
               << "Mean = " << res_3clusters[1] << " "
               << "Variance = " << res_3clusters[4] << endl;

          cout << "Third_Cluster: "
               << "Mean = " << res_3clusters[2] << " "
               << "Variance = " << res_3clusters[5] << endl;

          cout << "Minimum Likelihood  =" << res_3clusters[7] << endl;
          cout << "BIC selection criteria =" << res_3clusters[6] << endl;

          cout << "  " << endl;
          cout << "  " << endl;

          cout << "Results from 4 clusters" << endl;
          cout << "  " << endl;

          cout << "First_Cluster:  "
               << "Mean =" << res_4clusters[0] << " "
               << "Variance = " << res_4clusters[4] << endl;

          cout << "Second_Cluster: "
               << "Mean = " << res_4clusters[1] << " "
               << "  Variance = " << res_4clusters[5] << endl;

          cout << "Third_Cluster: "
               << "Mean = " << res_4clusters[2] << " "
               << "  Variance = " << res_4clusters[6] << endl;

          cout << "Fourth_Cluster: "
               << "Mean = " << res_4clusters[3] << " "
               << "Variance = " << res_4clusters[7] << endl;

          cout << "Minimum Likelihood  =" << res_4clusters[9] << endl;
          cout << "BIC selection criteria =" << res_4clusters[8] << endl;

          cout << "  " << endl;
          cout << "  " << endl;

          if ((res_4clusters[8] < res_3clusters[6]) && (res_4clusters[8] < res_2clusters[4]))
          {
               cout << "Grouping into 4 clusters gives the best result" << endl;
          }
          else if ((res_3clusters[6] < res_4clusters[8]) && (res_3clusters[6] < res_2clusters[4]))
          {
               cout << "Grouping into 3 clusters gives the best result" << endl;
          }
          else
          {
               cout << "Grouping into 2 clusters gives the best result" << endl;
          }
     }
     catch (const Max_Expectation::variance_must_be_positive &e)
     {
          cout << "Expectation Maximisation failed: Negative variance encountered.\n";
     }
     catch (const could_not_open_Datafile &e)
     {
          cout << "Error occured: Data File could not be opened\n";
     }
     catch (const size_must_match &e)
     {
          cout << "Error: vector sizes do not much .\n";
     }
}
