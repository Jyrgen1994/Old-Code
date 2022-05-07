#include <iostream>
#include <Eigen/Dense>
#include "Bspline.h"
#include <fstream>
#include <vector>
#include <limits>
#include <algorithm>
#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;


class Functions : public Bspline
{
private:
    double a_0,e;
    int Z;
public:
    Functions(int Z = 1,int order = 4,int n = 22,int resolution = 1000,int get_coeff = 0,int potential = 0,int is_lin_spaced = 0,double grid_min = 0.0,double grid_max =10.5);
        
    double bsplineFunction(int i , int j);
    vector<double> hamiltonianFunction(double r,int i, int j,int l_quantum);
    vector<double> gaussianQuadratureHam(double range_low,double range_high ,int i , int j ,int l_quantum);
    double gaussianQuadratureRho(double range_low,double range_high);
    double gaussianQuadratureV(double range_low,double range_high ,int coeff);
    double gaussianQuadratureB(double range_low,double range_high ,int i , int j ,int l_quantum);
    double hydrogen_s1(double r,double pi);
    double Vee(int coeff,double r);
    //VectorXd chargeDistribution(vector<tuple<int,int>> N_vals,vector<MatrixXd> P_l);
    double exchangePotential(double rho_r);
    double rho(double r ,VectorXd Bspline_final);

} ;