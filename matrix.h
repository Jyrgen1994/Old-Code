#include <Eigen/Dense>
#include <Eigen/Core>
#include <array>
#include <vector> 
#include <complex>
#include "Potential.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <fstream>
#include <ctime>
#include <limits>
#include <iomanip>

using Eigen::MatrixXcd;
using Eigen::internal::BandMatrix;
using namespace std;	
using namespace Eigen;

typedef SparseMatrix<double,ColMajor> SpMat;

class matrix{
/*
	Matrix class, used to intialize, and constuct different kind of matirces,
	such as diagonal, banded-diagonal, or any desired sparse system.
*/
private:
	MatrixXcd MatrixC,eigenVectors;
	MatrixXd MatrixR;
	VectorXcd inverse_power_solution,solution,left_hand_side,eigenValues;
	int size;
public:
	matrix(int size_in)
	{
		MatrixC = MatrixXcd::Constant(size_in,size_in,0.0);
		MatrixR = MatrixXd::Constant(size_in,size_in,0.0);
		solution = VectorXcd::Constant(size_in,0.0);
		left_hand_side  = VectorXcd::Constant(size_in,0.0);
		eigenValues = VectorXcd::Constant(size_in,0.0);
		eigenVectors = MatrixXcd::Constant(size_in,size_in,0.0);
		inverse_power_solution = VectorXcd::Constant(size_in,0.0);
		size = size_in;
	}

	void constructAlphaMatrix(int numOfBoundaries,vector<double> k,vector<double> xBoundaries)
	{
		complex<double> j(0.0,1.0);//define imaginary unit
		int matrix_size = 2*((k.size() - 2) +1) + 1;
		if (numOfBoundaries != (k.size() -1))
		{
			cout<<"Error : Number of boundaries does not match matrix size."<<endl;
		}
		else
		{
			//k contains both boundary ks & all kappas (2(M + 1)) in terms of black board calc.
			// size = 2(M+1) +1 for setting A = 1.
			
			int k_counter = 0; 
			for (int i = 0; i < matrix_size;i++)
			{
				if (i == 0)
				{
					MatrixC(i,i) = 1.0;
				}
				//FIRST BOUNDARY
				else if (0 < i &&  i < 3)
				{
					cout << "check_first : " <<i<< endl;
					if (i % 2 != 0)
					{
						MatrixC(i,i-1) =  exp(j*k[k_counter]*xBoundaries[k_counter]);
						MatrixC(i,i)   =  exp(-j*k[k_counter]*xBoundaries[k_counter]);
						MatrixC(i,i+1) = -exp(k[k_counter + 1]*xBoundaries[k_counter]);
						MatrixC(i,i+2) = -exp(-k[k_counter + 1]*xBoundaries[k_counter]);
					//	cout<< "test1 :"<< i <<endl;	
					}
					else
					{
						MatrixC(i,i-2) =  j*k[k_counter]*exp(j*k[k_counter]*xBoundaries[k_counter]);
						MatrixC(i,i-1) = -j*k[k_counter]*exp(-j*k[k_counter]*xBoundaries[k_counter]);
						MatrixC(i,i)   = -k[k_counter + 1]*exp(k[k_counter + 1]*xBoundaries[k_counter]);
						MatrixC(i,i+1) =  k[k_counter + 1]*exp(-k[k_counter + 1]*xBoundaries[k_counter]);
						k_counter++;	
					//	cout<< "test2 : "<< i << " : " << Matrix(i,i)<<endl;
					}
				}
				//BULK POTENTIALS
				else if(3 <= i && i < (matrix_size - 2) && k.size()>3)
				{
					if ((i % 2) != 0)
					{
						MatrixC(i,i-1) =  exp(k[k_counter]*xBoundaries[k_counter]);
						MatrixC(i,i)   =  exp(-k[k_counter]*xBoundaries[k_counter]);
						MatrixC(i,i+1) = -exp(k[k_counter + 1]*xBoundaries[k_counter]);
						MatrixC(i,i+2) = -exp(-k[k_counter + 1]*xBoundaries[k_counter]);
					}
					else
					{
						MatrixC(i,i-2) =  k[k_counter]*exp(k[k_counter]*xBoundaries[k_counter]);
						MatrixC(i,i-1) = -k[k_counter]*exp(-k[k_counter]*xBoundaries[k_counter]);
						MatrixC(i,i)   = -k[k_counter + 1]*exp(k[k_counter + 1]*xBoundaries[k_counter]);
						MatrixC(i,i+1) =  k[k_counter + 1]*exp(-k[k_counter + 1]*xBoundaries[k_counter]);
						k_counter ++;
					}
				}
				else if (i >= (matrix_size -2))
				{
					cout << "check_end : " <<i<< endl;
					if (i % 2 != 0)
					{
						MatrixC(i,i-1) =  exp(k[k_counter]*xBoundaries[k_counter]);
						MatrixC(i,i)   =  exp(-k[k_counter]*xBoundaries[k_counter]);
						MatrixC(i,i+1) = -exp(j*k[k_counter+1]*xBoundaries[k_counter]);
					}
					else
					{
						MatrixC(i,i-2) = k[k_counter]*exp(k[k_counter]*xBoundaries[k_counter]);
						MatrixC(i,i-1) = -k[k_counter]*exp(-k[k_counter]*xBoundaries[k_counter]);
						MatrixC(i,i)   = -j*k[k_counter + 1]*exp(j*k[k_counter + 1]*xBoundaries[k_counter]);
					}
				}

			}
		}
		
	}

	void constructSEMatrix(VectorXd grid,int band_width)
	{
		Potential pot(1);
		double del_x = grid(1)-grid(0);
		//==============FILLING BOUNDARIES FIRST===========================
		MatrixR.block(0,0,1,3) << -30 - 24*del_x*del_x*pot.PotentialHO2(grid(0)),16,-1;
		MatrixR.block(1,0,1,4) << 16, -30 - 24*del_x*del_x*pot.PotentialHO2(grid(1)),16,-1;
		MatrixR.block(size-2,size-4,1,4) << -1,16,-30 - 24*del_x*del_x*pot.PotentialHO2(grid(grid.size()-2)),16;
		MatrixR.block(size-1,size-3,1,3) << -1,16,-30 - 24*del_x*del_x*pot.PotentialHO2(grid(grid.size()-1)); 

		//==============FILLING BULK OF MATRIX====================
		for (int i = 2;i<size-2;i++)
		{
			MatrixR.block(i,i-2,1,5)<<-1,16,-30 -24*del_x*del_x*pot.PotentialHO2(grid(i)),16, -1;
		}
		MatrixR *= -(1.0/(24*del_x*del_x));
	}

	void solveSystem()
	{
		solution = MatrixC.colPivHouseholderQr().solve(left_hand_side);
	}

	void EigenSolve()
	{
		//SelfAdjoint
		time_t time_req;
		time_req = clock();
		SelfAdjointEigenSolver<MatrixXd> ES(MatrixR);
		eigenVectors = ES.eigenvectors();
		eigenValues = ES.eigenvalues();
		time_req = clock() -time_req;
		cout <<"Time Elapsed For EigenSolve: "<<std::setprecision(15)<<(double)time_req/CLOCKS_PER_SEC<<endl;

	}

	void writeEigenValues()
	{
		ofstream file_EV;
		file_EV.open("eigenValues.txt");
		for (int i = 0; i<eigenValues.size();i++)
		{
			file_EV << setprecision(15)<<real(eigenValues(i))<<"\n";
		}
		file_EV.close();
	}


	void writeEigenVectors()
	{
		ofstream file_EV;
		file_EV.open("EigenVector.txt");
		for (int i = 0 ; i<eigenVectors.row(0).size();i++)
		{
			for (int j = 0; j < eigenVectors.col(0).size();j++)
			{
				file_EV <<setprecision(15)<< real(eigenVectors(j,i))<<","<<imag(eigenVectors(j,i))<<"\n";
			}
			file_EV<<"\n";
		}
		file_EV.close();
	}

	void writeEigenVectorIterative(MatrixXd eigenVector)
	{
		ofstream file_EVI;
		file_EVI.open("EigenVectorsI.txt");
		eigenVector.normalize();
		for(int i = 0;i<eigenVector.row(0).size();i++)
		{
			file_EVI<<setprecision(15)<<eigenVector.col(i)<<"\n\n";
		}
		file_EVI.close();
	}

	void writeEigenValuesIterative(VectorXd eigen_values)
	{
		ofstream file_EVI;
		file_EVI.open("EigenValuesI.txt");
		for (int i = 0; i<eigen_values.size();i++)
		{
			file_EVI<<setprecision(15)<<eigen_values(i)<<"\n";
		}
		file_EVI.close();
	}

	
	double inversePowerIteration(VectorXd initial_guess,double tolerance,int numEigenSol,double del_z)
	{
		/*
		METHOD:
		Inverse Power Iteration method used to calculate lowest Eigenvalue of system.
		Add shift to calculate higher eigenvalues.
		INPUT:
		Initial guess for eigenvector
		Tolerance for convergence
		size of shift
		OUTPUT:
		lowest eigenvalue for the shifted system
		*/
		time_t time_req;
		ofstream file_EVI_vec,file_EVI_val;
		file_EVI_vec.open("EigenVectorsI.txt");
		file_EVI_val.open("EigenValuesI.txt");
		MatrixXd ID =  MatrixXd::Identity(size,size),matrix_shifted = MatrixXd::Constant(size,size,0.0),eigenVectors(size,numEigenSol);
		SpMat spare_shift_matrix;
		VectorXd old_sol,new_sol,eigenValues(numEigenSol);
		SparseLU<SparseMatrix<double, ColMajor>, COLAMDOrdering<Index> > solver;
		double err ,lambda_new ,lambda_old,shift;
		int iter;
		time_req = clock();
		for(int i = 0;i <numEigenSol;i++)
		{	
			
			shift = 0.4 + 0.8*i;
			old_sol = initial_guess;
			new_sol = VectorXd::Constant(size,1.0);
			matrix_shifted = MatrixR - shift*ID;
			//MatrixXd Matrix_inverse = matrix_shifted.inverse();
			spare_shift_matrix = matrix_shifted.sparseView();
			//cg.compute(spare_shift_matrix);
			solver.analyzePattern(spare_shift_matrix);
			solver.factorize(spare_shift_matrix);
			err =1.0;
			lambda_new = 0.0;
			lambda_old =10000.0;	
			iter = 0;
			while(tolerance <err)
			{
				iter++;
				new_sol = solver.solve(old_sol);
				lambda_new = (old_sol.dot(old_sol))/(old_sol.dot(new_sol));
				err = fabs(lambda_new-lambda_old);	
				lambda_old = lambda_new;
				old_sol = new_sol;
			}
			eigenVectors.col(i) = new_sol.normalized();
			cout<<"iterations  "<<iter<<endl;
			eigenValues(i) = lambda_new + shift;
			//cout << "iterations for convergence : "<<iter <<" for "<<i<<"th eigenvector/value"<<endl;
		}
		writeEigenVectorIterative(eigenVectors);
		writeEigenValuesIterative(eigenValues);
		time_req = clock() - time_req;
		cout <<"time Elapsed for inversePowerIteration: "<<std::setprecision(15)<<(double)time_req/CLOCKS_PER_SEC<<endl;
		return time_req;
	}

	VectorXcd getSolution()
	{
		return solution;
	}

	MatrixXcd getMC()
	{
		return MatrixC;
	}

	MatrixXd getMR()
	{
		return MatrixR;
	}

};
