#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Bspline.h"
#include <fstream>
#include <vector>
#include <limits>
#include <algorithm>
#include <Eigen/Eigenvalues>
#include <sstream>
#include<Eigen/IterativeLinearSolvers> 	
#include<Eigen/SparseLU>

using namespace std;
using namespace Eigen;


Bspline::Bspline(int order,int n,int resolution,int get_coeff,int potential,int is_lin_spaced,double grid_min,double grid_max)
{
	/*Consructor for B-splines*/
	//default value for order = 4
	//num_knotpoints = number of knot points including ghost points (2*(k-1))

	/*===========================INITIALIZATION=============================*/
	
	//cout.precision(dbl::max_digits10);
	this->is_lin = is_lin_spaced;
	this->potential_choice = potential;
	this->order               = order;
	this->num_knotpoints      = n;
	this->num_Bsplines = num_knotpoints - order;
	this->counter_checker = 0;
	this->get_coeff = get_coeff;


	this->NU = 1;
	this->V_1 = 0.0;
	this->V_2 = 0.0;
	this->R_1 = 1.0;
	this->R_2 = 2.0;
	this->delta_step = 0.01;

	this->num_physical_points = num_knotpoints - 2*(order-1);
	if(get_coeff != 0)//get_coeff != 0 when calculating coefficients
	{
		resolution = num_physical_points;
	}	
	this->grid_min = grid_min;
	this->grid_max = grid_max;
	this->grid_res = resolution;

	this->t = VectorXd::Zero(num_knotpoints);
	this->t_min = VectorXd::Ones(order - 1)*(this->grid_min );
	this->t_max = VectorXd::Ones(order - 1)*(this->grid_max );
	this->grid     = VectorXd::Zero(grid_res);
	this->Bsplines = MatrixXd::Zero(num_knotpoints ,order);
	this->Bspline_d1 = VectorXd::Zero(num_knotpoints);
	this->Bspline_d2 = VectorXd::Zero(num_knotpoints);
	this->Bspline_final = VectorXd::Zero(num_knotpoints);
	this->A_B   = MatrixXd::Zero(resolution, num_Bsplines);
	this->A_Bd1 = MatrixXd::Zero(resolution, num_Bsplines);
	this->A_Bd2 = MatrixXd::Zero(resolution, num_Bsplines);
	this->M_d2  = MatrixXd::Zero(num_Bsplines,num_Bsplines);
	this->function = VectorXd::Zero(num_physical_points);
	this->H = MatrixXd::Zero(num_Bsplines,num_Bsplines); 
	this->B = MatrixXd::Zero(num_Bsplines,num_Bsplines);
	this->V_matrix = MatrixXd::Zero(num_Bsplines,num_Bsplines);
	this->EigenVectors = MatrixXd::Zero(num_Bsplines-2,num_Bsplines-2);

	this->P_nl_coeff = { };
	this->Energies = { };
	this->rho_vec = VectorXd::Zero(grid_res);
	this->V_ee_old = { };
	this->V_ee_olds = {};
	this->eta = 0.6;
	this->V_ee_cnt =0;
	this->convergence_help = false;
	this->H_list = { };
	this->Vee_list = { };
	/*=====================================================================*/

}

void Bspline::setUpBSplines(double x,int is_Eigen_problem)
{
	/* Method: 
	 * Fills Bspline matrix (grid x order) with Bsplines of order up to
	 * this->order (=4 default)  for a given 
	 * grid point (input).
	 */
	//counter_checker = 0;
	int j;
	double tol = 1e-10;
	double term_1,term_2;


	/*===========RESET BSPLINE FOR EACH GRID POINT=========*/
	Bsplines = MatrixXd::Zero(num_knotpoints ,order);
	Bspline_final = VectorXd::Zero(num_knotpoints);
	Bspline_d1 = VectorXd::Zero(num_knotpoints);
	Bspline_d2 = VectorXd::Zero(num_knotpoints);
	/*=====================================================*/

	for (int k = 0; k<order; k++)
	{
		for (int i = 0; i < num_Bsplines;i++)
		{	
			
			/*for order = 1: */
			if (k == 0)
			{
				if (t(i) <= x && x < t(i+1))
				{
					Bsplines(i,k) = 1.0;
				}
				if (fabs(x- t(num_knotpoints-1)) < tol && i == num_Bsplines-1)
				{
					Bsplines(i,k) = 1.0;
				}
			}
			/*for order  >1: */
			else 
			{
				if (fabs(Bsplines(i,k-1)) > tol)
				{
					Bsplines(i,k) = (x - t(i))/(t(i + k) -t( i ))*Bsplines(i,k-1);
				}

				if (fabs(Bsplines(i+1,k-1)) > tol)
				{
					Bsplines(i,k) += (t(i+k+1) - x)/(t( i + k + 1 ) - t( i + 1))*Bsplines(i+1,k-1);
				}

				if (fabs(x - t(num_knotpoints - 1)) < tol && i == num_Bsplines- 1)
				{
					Bsplines(i,k) = 1.0;
				}
			}
			if (k == order-1)
			{
				/*Get Derivatives of Bslines ORDER 1:*/
				if(fabs(Bsplines(i,k-1)) >tol  )
				{
					Bspline_d1(i) += (order -1)*Bsplines(i,k-1)/(t(i + k) - t(i));
				}
				if(fabs(Bsplines(i+1,k-1)) >tol )
				{
					Bspline_d1(i) -= (order -1)*Bsplines(i+1,k-1)/(t(i + k + 1) - t(i + 1));
				}
				/*ORDER 2:*/
				if(fabs(Bsplines(i,k-2)) >tol )
				{
					Bspline_d2(i) += (order - 1)*(order - 2)*Bsplines(i,k-2)/((t(i + k) -t(i))*(t(i + k - 1) -t(i)));
				}

				if(fabs(Bsplines(i+1,k-2)) >tol)
				{
					term_1 = 1.0/((t(i + k) -t(i))*(t(i + k) - t(i + 1)));
					term_2 = 1.0/((t(i + k + 1) - t(i+1))*(t(i + k) - t(i + 1)));
					Bspline_d2(i) -= (order - 1)*(order - 2)*Bsplines(i+1,k-2)*(term_1 + term_2);
				}

				if(fabs(Bsplines(i+2,k-2)) >tol )
				{
					Bspline_d2(i) += (order -1)*(order - 2)*Bsplines(i+2,k-2)/((t(i + k + 1) -t(i + 1))*(t(i + k + 1) -t(i + 2)));
				}	 
			}
		}
		
		//Fill Matrices holding Bsplies,first & second derivatives:
		if(k == order-1 && !is_Eigen_problem)
		{
			A_B.row(counter_checker) = Bsplines.col(k).segment(0,num_Bsplines);
			A_Bd1.row(counter_checker) = Bspline_d1.segment(0,num_Bsplines);
			A_Bd2.row(counter_checker) = Bspline_d2.segment(0,num_Bsplines);
			counter_checker ++;
		}
		else
		{
			Bspline_final = Bsplines.col(k).segment(0,num_Bsplines);
		}
	}	
}


void Bspline::setUpKnotSequence()
{
	/*Initialize knot sequence*/
	
	int N = num_physical_points;
	double exp_min,exp_max;
	double delta = 1e-8;
	t.segment(0,order-1) = t_min;
	VectorXd extra_points = VectorXd::Zero(NU);
	extra_points << 10000;

	
	if(is_lin == 1)
	{
		/*========================LINSPACED T==================*/
		t.segment(order-1,num_physical_points) = VectorXd::LinSpaced(num_physical_points,grid_min,grid_max);
		t.segment(num_knotpoints-(order-1),order-1) = t_max;
		/*========== CORRESPONDING GRID ===================================*/
		grid =VectorXd::LinSpaced(grid_res,grid_min + delta_step,grid_max - delta_step);
	}
	else if( is_lin == 2 )
	{
		/*================EXPONENTIAL GRID =========================*/
		double min = -7.0;
		double max = log10(grid_max);
		t.segment(0,order-1) = VectorXd::Ones(order-1)*pow(10,min);
		t.segment(num_knotpoints -(order-1),order-1) = VectorXd::Ones(order-1)*pow(10,max);
		//cout<<"TEST ORDER\n"<<VectorXd::Ones(order-1)*exp(min)<<endl;
		VectorXd log_spaced_grid = VectorXd::LinSpaced(grid_res,min,max);
		VectorXd log_spaced_knot = VectorXd::LinSpaced(num_physical_points,min,max);
		for(int i = 0; i < log_spaced_grid.size();i++)
		{
			grid(i) = pow(10,log_spaced_grid(i));

		}		
		for(int i = 0; i < log_spaced_knot.size();i++)
		{
			t( (order-1) +i) = pow(10,log_spaced_knot(i));
		}
	}
	else
	{
		/*=============== NON-LINSPACE T===============================*/
		t.segment(3,N-NU) = VectorXd::LinSpaced(N-NU,grid_min,grid_max);
		t.segment(num_knotpoints-(order -1) -NU , order - 1) = t_max;
		t.segment(num_knotpoints-NU,NU) = extra_points;
		/*========== CORRESPONDING GRID ===================================*/
		grid.segment(0,grid_res-NU) = VectorXd::LinSpaced(grid_res-NU,grid_min + delta_step,grid_max-delta_step);
		grid.segment(grid_res-NU,NU) = extra_points;
	
	}
	
		

	/*EXTRA SHIT*/
	sortGrid();
	sortKnot();
	if ( grid_max > t(num_knotpoints-1) )
	{
		cout<<"in here?"<<(100.0>t(num_knotpoints-1))<<endl;
		t.segment(num_knotpoints-(order-1),(order-1))<<t(num_knotpoints-1),t(num_knotpoints-1),t(num_knotpoints-1);
		cout<<"assign fail?"<<endl;
	}
	//cout << "T !!!! "<<grid_min<<" "<<grid_max<<" "<< num_physical_points<<"\n" <<t<<endl;
}


void Bspline::sortGrid()
{
	vector<double> sort_vec;
	for(int i = 0;i<grid.size();i++)
	{
		sort_vec.push_back(grid(i));
	}
	sort(sort_vec.begin(),sort_vec.end());
	for (int i = 0;i<grid.size();i++)
	{
		grid(i) = sort_vec[i];
	}

}

void Bspline::sortKnot()
{
	vector<double> sort_vec;
	for(int i = 0;i<t.size();i++)
	{
		sort_vec.push_back(t(i));
	}
	sort(sort_vec.begin(),sort_vec.end());
	for (int i = 0;i<t.size();i++)
	{
		t(i) = sort_vec[i];
	}

}


VectorXd Bspline::getFunction(VectorXd coefficients)
{
	function = A_B*coefficients;
	return function;
}

VectorXd Bspline::getEigenFunction(VectorXd coefficients)
{
	MatrixXd A_B_block = A_B.block(0,1,A_B.rows(),num_Bsplines-2);
	function = A_B_block*coefficients;
	return function;
}

vector<vector<int>> Bspline::getOccNumebers(int nr_electrons)
{
	int e_left = nr_electrons;//only for neutral atoms at the moment...
    int full_shell;
    int n_quantum = 1;
    int l_quantum =0;
    vector<vector<int>> N_vals;
    vector<MatrixXd> P_l;
    vector<int> temp = {0 , 0 , 0};
	while (e_left != 0)//find occupation numbers
    {
        full_shell = 2*(2*l_quantum + 1 );
        if ( e_left >= full_shell)//check if enough electrons left to fill shell
        {
            temp[0] = full_shell;
			temp[1] = l_quantum;
			temp[2] = n_quantum;
            N_vals.push_back(temp);  
			N_occ_internal.push_back(temp);
            e_left -= full_shell;
            if( (l_quantum == n_quantum -1) && (e_left != 0) )
            {
                n_quantum++;
                l_quantum = 0;
            }
            else
            {
                l_quantum++;
            }
            
        }
        else //if shell not full, take leftover electrons 
        {
            temp[0] = e_left;
			temp[1] = l_quantum;
			temp[2] = n_quantum;
            N_vals.push_back(temp);
			N_occ_internal.push_back(temp);
            e_left -= e_left;
        }
    }
	vector<int> old = N_vals[0];
	for(int i = 1; i < N_vals.size();i++)
	{
		if(old[1] > N_vals[i][1])
		{
			swap(N_vals[i-1],N_vals[i]);
		}
		old = N_vals[i];
	}
	N_occ_internal = N_vals;
	return N_vals;
}

int Bspline::getMaxL(vector<vector<int>> N_occ)
{
	int l_max = 0;
	for(int i = 0; i <N_occ.size();i++)
	{	
		if(N_occ[i][1] >l_max)
		{
			l_max = N_occ[i][1];
		}
	}
	return l_max;
}

vector<VectorXd> Bspline::getEigenCoeffs(int l_max , vector<vector<int>> N_occs, vector<VectorXd> coeffs,int Z)
{
	vector<VectorXd> coeff_eigen;
	vector<VectorXd> eigen_vals;
	P_nl_coeff = coeffs;
	if(P_nl_coeff.size() != 0)
	{
		
		for(int i = 0; i < grid.size();i++)
		{	
			rho_vec(i) = rho(grid(i),A_B.row(i));
		}
		//rho_vec = normalize(rho_vec);
		//cout<<rho_vec<<endl;
		xi = solveSystemCoefficients(rho_vec,Z);
	}
	for(int i = 0; i < l_max + 1 ;i++)
	{

		createEigenSystem(i);
		//cout<<"aftercreate : "<<V_ee_olds.size()<<endl;
		solveEigenEquation(i);	
		H_list.push_back(H);
		Vee_list.push_back(V_matrix);
		for( int j = 0; j <N_occs.size();j++)
		{
			if (N_occs[j][1] == i)
			{
				coeff_eigen.push_back(EigenVectors.col(N_occs[j][2]-i-1));//store eigenvectors(coefficients) corresponding to different l values}
				//cout<<"here??? CEHCK I "<<i<<"  "<<j<<endl;
			}
		}
		eigen_vals.push_back(EigenValues);
		V_ee_olds[i] = V_ee_old;	
		V_ee_old = { };
		V_ee_cnt = 0;
	}
	Energies = eigen_vals;
	return coeff_eigen;
}


void Bspline::createEigenSystem(int l_quantum)
{
	createSystem();
	M_d2_sparse = M_d2.sparseView();
	//cg.compute(M_d2_sparse);
	solver.analyzePattern(M_d2_sparse);
	solver.factorize(M_d2_sparse);
	H = MatrixXd::Zero(num_Bsplines,num_Bsplines);
	B = MatrixXd::Zero(num_Bsplines,num_Bsplines);
	V_matrix = MatrixXd::Zero(num_Bsplines,num_Bsplines);
	vector<double> res;
	for(int i = 0; i < H.cols(); i++)
	{
		//cout<<"test i "<<i <<endl;
		for(int j = 0; j < H.rows();j++)
		{
			if ((min(i,j) + order) >(max(i,j)))
			{
				for (int k = max(i,j); k <(min(i,j) + order);k++)
				{
					if (k >= (order -1))
					{
						res = gaussianQuadratureHam(t(k),t(k+1) , i , j , l_quantum);
						H(i,j) +=  res[0];
						B(i,j) +=  gaussianQuadratureB(t(k),t(k+1) , i , j ,  l_quantum);
						V_matrix(i,j) += res[1];
					}
				}
			}
			
		}

	}
	
}

void Bspline::solveEigenEquation(int l_quantum)
{	
	MatrixXd H_block(num_Bsplines-2,num_Bsplines-2),B_block(num_Bsplines-2,num_Bsplines-2);
	H_block = H.block(1,1,num_Bsplines-2,num_Bsplines-2);
	B_block = B.block(1,1,num_Bsplines-2,num_Bsplines-2);

	GeneralizedSelfAdjointEigenSolver<MatrixXd> ges;
	ges.compute(H_block,B_block); 
	EigenVectors = ges.eigenvectors().real();
	EigenValues = ges.eigenvalues().real();
	writeEigenSolution(ges.eigenvalues(),EigenVectors,l_quantum);
}

void Bspline::writeEigenSolution(VectorXcd eigen_vals,MatrixXd eigen_vecs,int l_quantum)
{
	ofstream file_vals;
	ofstream file_vecs;
	string file_name_vals;
	string file_name_vecs;

	file_name_vals = "EigenVals" + to_string(l_quantum) +".txt";
	file_name_vecs = "EigenVecs" + to_string(l_quantum) +".txt";

	file_vals.open(file_name_vals);
	file_vecs.open(file_name_vecs);
	for (int i = 0 ; i < eigen_vals.size();i++)
	{
		file_vals<< real(eigen_vals(i))<<"\n";
		for(int j = 0 ; j < eigen_vecs.cols();j++)
		{
			eigen_vecs.col(i).normalize();
			file_vecs<< eigen_vecs.col(i)(j)<<"\n";	
		}
		file_vecs << "\n";
	}
}


void Bspline::createSystem()
{
	MatrixXd block = A_Bd2;
	M_d2.conservativeResize(num_Bsplines,num_Bsplines);
	M_d2.block(1,0,num_physical_points,num_Bsplines) = block;
	M_d2(0,0) = 1.0;
	M_d2(M_d2.rows()-1,M_d2.cols()-1) = 1.0;
}

VectorXd Bspline::solveSystemCoefficients(VectorXd charge_density,int Z)
{
	VectorXd B(num_Bsplines),solution(num_Bsplines),phys_grid(num_physical_points);
	phys_grid = t.segment(3,num_physical_points);
	
	//createSystem();
	//M_d2_sparse = M_d2.sparseView();
	
	double pi = acos(-1),Q =1.0,rho,a_0 = 5.29e-11;
	V_1 = ((4.0/3.0)*pi *R_1*R_1*R_1);
	V_2 = ((4.0/3.0)*pi *R_2*R_2*R_2);
	
	for(int i = 0; i < num_physical_points;i++)
	{

		if (potential_choice == 0)
		{
			/*Homogeneous sphere*/
			if ( phys_grid(i) <= R_1)
			{
				rho = Q/V_1;
			} 
			else
			{	
				rho = 0.0;
			}
		}
		else if (potential_choice == 1)
		{
			/*Homogeneous shell*/
			if(R_1 <= phys_grid(i) && phys_grid(i) <= R_2)
			{
				rho = Q/(V_2-V_1);
			}
			else 
			{
				rho = 0.0;
			}
		}
		else if (potential_choice == 2)
		{
			rho =  hydrogen_s1(phys_grid(i),pi);
		}
		else if (potential_choice == 3)
		{
			rho = charge_density(i);
		}
		else
		{
			cout << "Potential does not exist.";
		}
		
		B(i+1) = -4.0*phys_grid(i)*pi*rho;
	}
	B(0) = 0.0;
	B(num_Bsplines-1) = 1.0*Z;
	cout<<"TETING!!! "<<1.0*Z<<endl;
	
	/*Solve system of coefficients*/
	solution = solver.solve(B);
	//solution = cg.solve(B);
	//solution = M_d2.householderQr().solve(B);
	return solution;
}

VectorXd Bspline::normalize(VectorXd wave_function)
{
	VectorXd normalized_wave_function;
	VectorXd phys_grid = t.segment(order, t.size()-order);
	double norm_factor = 0.0;
	for (int i = 0; i < phys_grid.size()-1;i++)
	{
		norm_factor += gaussianQuadratureRho(phys_grid(i),phys_grid(i+1));
	}
	norm_factor = pow(norm_factor,0.5);
	normalized_wave_function = wave_function*(1.0/norm_factor);
	return normalized_wave_function;
}

MatrixXd Bspline::getBSpline()
{
	return A_B;	
}


MatrixXd Bspline::getBSplineFirstDeriv()
{
	return A_Bd1;
}

MatrixXd Bspline::getBSplineSecondDeriv()
{
	return A_Bd2;
}

MatrixXd Bspline::getHamiltonian()
{
	return H;
}

MatrixXd Bspline::getBSplineMatrix()
{
	return B;
}




void Bspline::writeBSplines(MatrixXd results)
{
	ofstream file;
	file.open("Bsplines.txt");
	for (int i = 0; i <grid_res ;i++)
	{
		file<<results.row(i).transpose()<<endl;
		file<<" "<<endl;
	}
	file.close();
}

void Bspline::writeBsplineSecondDerivative(MatrixXd results)
{
	ofstream file;
	file.open("Bsplines_second_deriv.txt");
	for (int i = 0; i <grid_res ;i++)
	{
		file<<results.row(i).transpose()<<endl;
		file<<" "<<endl;
	}
	file.close();
}

void Bspline::writePhi(VectorXd results,string file_name)
{
	ofstream file;
	file.open(file_name);
	file<<results;
	file.close();
}