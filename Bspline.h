#ifndef BSPLINE_H
#define BSPLINE_H

#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <limits>
#include <algorithm>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include<Eigen/IterativeLinearSolvers> 	
#include<Eigen/SparseLU>
using namespace std;
using namespace Eigen;

typedef std::numeric_limits< double > dbl;

class Bspline{
private:

	
public:
	int num_knotpoints,order,grid_res,num_Bsplines,num_physical_points,NU,potential_choice,is_lin,get_coeff;
	double grid_max,grid_min,delta_step,R_1,R_2,V_1,V_2;
	VectorXd t,Bspline_final,Bspline_d1,Bspline_d2,t_min,t_max,function;
	MatrixXd Bsplines,A_B,A_Bd1,A_Bd2,M_d2,H,B;
	SparseMatrix<double> M_d2_sparse;
	MatrixXd EigenVectors;
	MatrixXd V_matrix;
	VectorXd grid;
	VectorXd EigenValues;
	vector<vector<int>> N_occ_internal;
	vector<VectorXd> P_nl_coeff;
	vector<VectorXd> Energies; 
	VectorXd rho_vec;
	VectorXd xi;
	vector<double> V_ee_old;
	vector<vector<double>> V_ee_olds;
	bool convergence_help;
	vector<MatrixXd> H_list;
	vector<MatrixXd> Vee_list;

	VectorXd grid_Vdir;
	MatrixXd A_B_Vdir;
	
	double eta;
	int counter_checker,V_ee_cnt;
	ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
	SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;
	
	Bspline(int order = 4,int n = 22,int resolution = 1000,int get_coeff = 0,int potential = 0,int is_lin_spaced = 0,double grid_min = 0.0,double grid_max =10.5);

	void setUpKnotSequence();
	void setUpBSplines(double x,int is_eigen_problem);
	void sortKnot();	
	void sortGrid();
	void createSystem();
	void createEigenSystem(int l_quantum);
	void writeBSplines(MatrixXd result);
	void writeBsplineSecondDerivative(MatrixXd result);
	void writePhi(VectorXd results,string file_name);
	void solveEigenEquation(int l_quantum);
	void writeEigenSolution(VectorXcd eigen_vals,MatrixXd eigen_vecs,int l_quantum);


    virtual double bsplineFunction(int i , int j) = 0;
    virtual vector<double> hamiltonianFunction(double r,int i, int j,int l_quantum) = 0;
    virtual vector<double> gaussianQuadratureHam(double range_low,double range_high ,int i , int j ,int l_quantum) = 0;
    virtual double gaussianQuadratureB(double range_low,double range_high ,int i , int j ,int l_quantum) = 0;
	virtual double gaussianQuadratureV(double range_low,double range_high ,int coeff)=0;
    virtual double gaussianQuadratureRho(double range_low,double range_high) = 0	;
	virtual double hydrogen_s1(double r,double pi) = 0;
	//virtual VectorXd chargeDistribution(vector<tuple<int,int>> N_vals,vector<MatrixXd> P_l) = 0;
	virtual double exchangePotential(double r) = 0;
	virtual double rho(double r ,VectorXd Bspline_final) = 0;
	

	MatrixXd getBSpline();
	MatrixXd getBSplineFirstDeriv();
	MatrixXd getBSplineSecondDeriv();
	VectorXd solveSystemCoefficients(VectorXd charge_density,int Z);
	VectorXd getFunction(VectorXd coefficients);
	VectorXd getEigenFunction(VectorXd coefficients);
	MatrixXd getHamiltonian();
	MatrixXd getBSplineMatrix();


	vector<VectorXd> getEigenCoeffs(int l_max,vector<vector<int>> N_occs,vector<VectorXd> coeffs,int Z);
	vector<vector<int>> getOccNumebers(int nr_electrons);
	VectorXd normalize(VectorXd wave_function);
	int getMaxL(vector<vector<int>> N_occ);

};
#endif