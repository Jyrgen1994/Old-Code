#include <iostream>
#include <Eigen/dense>
#include "Bspline.h"
#include "Functions.h"
#include <vector>

using namespace Eigen;
using namespace std;


void fillZlist(vector<vector<int>> *Z_list)
{	
	for(int i = 2;i<19;i++)
	{
		(*Z_list).push_back({i,i-1});
	}
}

VectorXd squareVec(VectorXd vector)
{
	VectorXd res = VectorXd::Zero(vector.size());
	for(int i = 0; i < vector.size();i++)
	{
		res(i) = vector(i)*vector(i);
	}
	return res;
}

int main()
{
	/*==SOME INTIAL PARAMTERS==*/
	int n = 200;
	double Q = 1.0;
	int order = 4;
	double R_min = 0;
	double R_max = 1000.0;
	int Z = 10;
	int type_grid = 2;
	/*======================*/
	vector<double> E_ion = { };
	vector<vector<int>> Z_list = { };
	//fillZlist(&Z_list);
	Z_list.push_back({19,19});
	cout<<Z_list.size()<<endl;

	int resolution = n-2*(order-1);
	MatrixXd result(resolution,n),resultDeriv1(resolution,n),resultDeriv2(resolution,n);
	VectorXd temp(resolution),coefficients(n-order),sol(n-order);
	VectorXd V_dir,V_exch,V_ee;
	
	vector<VectorXd> eigen_coeff = { };
	vector<VectorXd> eigen_vectors = { };
	vector<vector<int>> N_occs = { };
	for(int i = 0; i < Z_list.size();i++)
	{
		for(int j = 0; j <2;j++)
		{
			Bspline *coeff = new Functions(Z_list[i][0],order,n,resolution,1,3,type_grid,R_min,R_max);
			//Bspline *V_dir = new Functions(Z,order,n,resolution,1,3,type_grid,R_min,10.0);
			
			
			coeff->setUpKnotSequence();
			//V_dir->setUpKnotSequence();

			for (int i = 0; i <resolution;i++)
			{
				coeff->setUpBSplines(coeff->grid(i),0);
			}


			//In N_occs we store orbital information
			//N_occs = coeff->getOccNumebers(Z_list[i][j]);
			if(j == 0)
			{
				N_occs = {{2,0,1},{2,0,2},{2,0,3},{6,1,2},{6,1,3},{1,2,3}};
				coeff->N_occ_internal = N_occs;
			}
			else
			{
				N_occs = {{2,0,1},{2,0,2},{2,0,3},{1,0,4},{6,1,2},{6,1,3}};
				coeff->N_occ_internal = N_occs;
			}
			
				
			

			
			
			
			for (int i = 0; i < N_occs.size();i++)
			{
				cout<<N_occs[i][0]<<" "<<N_occs[i][1]<<" "<<N_occs[i][2]<<" "<<endl;	
			}
			eigen_coeff = { };
			eigen_vectors = { };
			double Energy = 0.0;
			double Energy_old = 2000.0;
			double Vee = 0.0 ;
			int l_max;
			l_max = coeff->getMaxL(N_occs);
			for(int i = 0; i <l_max+1;i++)
			{
				coeff->V_ee_olds.push_back({});
			}
			cout<<"LMAX "<<l_max<<endl;
			//fabs(Energy-Energy_old) >1.0;
			int iter =0;
			while( fabs(Energy-Energy_old) >0.1)
			{	
				Energy_old = Energy;
				Energy = 0.0;
				if(iter > 1)
				{
					coeff->convergence_help = true;
				}
				eigen_coeff = coeff->getEigenCoeffs(l_max,N_occs,eigen_coeff,Z_list[i][j]);
				int cnt = 0;
				for(int m = 0; m < l_max+1; m++)
				{
					for(int j = 0; j < N_occs.size();j++)
					{
						if(N_occs[j][1] == m)
						{
							cout<<(0.5*(coeff->Vee_list[N_occs[j][1]]).block(1,1,coeff->num_Bsplines-2,coeff->num_Bsplines-2)*eigen_coeff[j]).dot(eigen_coeff[j])<<endl;
							Vee =(0.5*(coeff->Vee_list[N_occs[j][1]]).block(1,1,coeff->num_Bsplines-2,coeff->num_Bsplines-2)*eigen_coeff[j]).dot(eigen_coeff[j]);
							Energy += N_occs[j][0]*(coeff->Energies[m](N_occs[j][2] -m -1) - Vee);
							//cout<<"TEST cnt : "<<cnt<<" j : "<<j<<endl;
							cnt++;
							//cout<<"GET EIGEN COEFFS : "<<m<<" "<<N_occs[j][2] -m -1<<" "<<N_occs[j][0]<<" and "<<N_occs[j][1]<<" and "<<N_occs[j][2]<<endl;
							cout<<"test orbital energies "<<coeff->Energies[m](N_occs[j][2] -m -1)<<" : "<<N_occs[j][2]<<endl;
						}
					}
				}
				coeff->H_list = { };
				coeff->Vee_list = { };
				cout<<"test E "<<fabs(Energy-Energy_old)<<" "<<Energy<<endl;
				eigen_vectors.push_back(coeff->rho_vec);
				iter++;
				//cout<<"H3r3??"<<endl;
			}
			E_ion.push_back(Energy);
			if (j == 0)
			{
				ofstream file_vectors;
				file_vectors.open("phiEigen.txt");
				for(int i = 0; i < eigen_vectors.size();i++)
				{
					file_vectors<<eigen_vectors[i]<<endl;
					file_vectors<<"\n";
				}
				file_vectors.close();
			}
			
			delete coeff;
		}
	}
	for(int i = 0; i < E_ion.size();i++)
	{
		cout<<"ENERGIES NA "<<E_ion[i]<<endl;
	}
	/*ofstream E_file;
	E_file.open("IonEnergy.txt");
	for (int i = 0; i < E_ion.size();i++)
	{
		E_file<<E_ion[i]<<endl;
		if(((i + 1) % 2) == 0 && i != 0)
		{
			E_file<<"\n";
		}
	}*/

			/*ofstream file_vectors;
			file_vectors.open("phiEigen.txt");
			for(int i = 0; i < eigen_vectors.size();i++)
			{
				file_vectors<<eigen_vectors[i]<<endl;
				file_vectors<<"\n";
			}
			file_vectors.close();*/

			//cout<<eigen_coeff[1].size()<<endl;
			//coeff->writePhi(coeff->getEigenFunction(eigen_coeff[1]),"Phi.txt");
			//coeff->writePhi(coeff->getEigenFunction(eigen_coeff[2]),"Phi2.txt");
			//cout<<N_occs.size()<<endl;	
	
	







	#pragma region Region_coeff//Coeff building part
	//coeff->chargeDistribution();
	/*coeff->createEigenSystem(0);
	coeff->solveEigenEquation(0);
	result = coeff->getBSpline();
	resultDeriv2 = coeff->getBSplineSecondDeriv(); 
	coefficients = coeff->solveSystemCoefficients();
	coeff->writeBSplines(result);
	cout<<"grid\n"<<coeff->grid<<endl;
	//cout<<resultDeriv2<<endl;
	*/
	#pragma endregion Region_coeff
	
	/*======= ACTUAL PLOT =======*/
	
	



	/*
	#pragma region Region_solve//Soloving part
	int resolution_2 = 2000;


	VectorXd sol2(resolution_2);
	VectorXd sol3(resolution_2);
	MatrixXd sol_arr(resolution_2,n-2*(order-1));



	Bspline *solution = new Functions(Z,order,n,resolution_2,0,0,type_grid,R_min,R_max);
	solution->setUpKnotSequence();
	for(int i = 0; i < resolution_2;i++)
	{
		solution->setUpBSplines(solution->grid(i),0);
	}
	cout<<"hej"<<endl;

	VectorXd charge_density = solution->chargeDistribution(N_occs,coeff_eigen);
	coefficients = coeff->solveSystemCoefficients(charge_density,Z);



	ofstream file; 
	file.open("eigenvector.txt");
	for(int i =0; i < sol_arr.cols();i++)
	{
		//cout<<"testing "<<i<<endl;
		sol_arr.col(i) = solution->getEigenFunction(coeff->EigenVectors.col(i));
		file<<solution->getEigenFunction(coeff->EigenVectors.col(i))<<"\n \n";
	}
	file.close();






/*
	sol2 = solution->getFunction(coefficients);
	sol3  = solution->getEigenFunction(coeff->EigenVectors.col(1));
	//cout<< sol3<<endl;
	solution->writeBSplines(solution->getBSpline());
	solution->writePhi(sol3);
	//solution->createEigenSystem();
	
	#pragma endregion Region_solve
	*/
	return 1;
}


