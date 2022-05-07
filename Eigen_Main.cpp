#include <iostream>
#include <Eigen/Dense>
#include "matrix.h"
#include <vector>
#include <Eigen/Sparse>
#include <ctime>
#include <limits>
#include <iomanip>

using namespace Eigen;
using namespace std;
typedef SparseMatrix<double> SpMat;


VectorXd InitialGuess(int N,VectorXd grid)
{
	VectorXd initial_guess = VectorXd::Constant(N,0.0);
	for (int i = 0;i<N;i++)
	{
		initial_guess(i) =1.0;
	}
	return initial_guess;

}


int main()
{	
	SpMat sparse;
	int N = 500;	
	MatrixXd MatrixSE(N,N);
	matrix mold_matrix(N);
	VectorXd initial_guess;

	VectorXd grid(N);
	VectorXd v;
	v = VectorXd::LinSpaced(N,-6,6);
	mold_matrix.constructSEMatrix(v,1);
	MatrixSE = mold_matrix.getMR();
	sparse = MatrixSE.sparseView();

	mold_matrix.EigenSolve();


	mold_matrix.writeEigenVectors();
	mold_matrix.writeEigenValues();

	initial_guess = InitialGuess(N,grid);
	double time_iter = mold_matrix.inversePowerIteration(initial_guess,1e-14,10,v(0)-v(1));	
	//cout << MatrixSE<<endl;
	//cout << print*print.transpose() << endl;
	return 1;
}
