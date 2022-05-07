#include "Functions.h"

Functions::Functions(int Z,int order ,int n ,int resolution ,int get_coeff ,int potential ,int is_lin_spaced,double grid_min,double grid_max ) : 
Bspline( order , n , resolution , get_coeff , potential , is_lin_spaced, grid_min, grid_max )
{
    this->a_0 = 1.0;
    this->e =1.0;
    this->Z = Z;
	cout<<"TEST Z "<<Z<<endl;
}


vector<double> Functions::hamiltonianFunction(double r,int i, int j,int l_quantum)
{
	//i,j is for Bspline i & j 
	double R_0 = 1.5;
	double res= 0.0,res_2 = 0.0;
    double V_ee,V_dir;
    double pi = acos(-1);

    if(P_nl_coeff.size() == 0)
    {
        V_ee = 0.0;
    }
    else
    {
        V_dir = (1.0/r)*xi.dot(Bspline_final);
        V_ee = exchangePotential(rho(r,Bspline_final)) + V_dir;
	    if(convergence_help)
        {
			//cout<<"in the convergence! "<<V_ee_olds.size()<<" "<<V_ee_olds[0].size()<<endl;
            V_ee = (1.0 - eta)*V_ee + eta*V_ee_olds[l_quantum][V_ee_cnt];
			V_ee_cnt++;
			//cout<<"in the convergece 2!"<<V_ee_cnt<<endl;
		}
        V_ee_old.push_back(V_ee);
    }
	res_2 = V_ee*Bspline_final(i)*Bspline_final(j);
	res = (a_0*e/2.0)*Bspline_d1(i)*Bspline_d1(j) + ((a_0*e/2.0)*l_quantum*(l_quantum + 1)/(pow(r,2.0)) - Z/r + e*V_ee)*Bspline_final(i)*Bspline_final(j);
	return {res,res_2};
}


double Functions::bsplineFunction(int i , int j)
{
	return Bspline_final(i)*Bspline_final(j); 
}

vector<double> Functions::gaussianQuadratureHam(double range_low,double range_high ,int i , int j ,int l_quantum)
{

	vector<double> absciccas = { -0.3399810435848563,0.3399810435848563,-0.8611363115940526,0.8611363115940526};
	vector<double> weights   = {0.6521451548625461,0.6521451548625461,0.3478548451374538,0.3478548451374538};
	double integral_result   = 0.0;
	double integral_result_2 = 0.0;
	double r = 0.0;
	vector<double> res;
	for (int m = 0; m < absciccas.size();m++)
	{
		r = absciccas[m]*(range_high - range_low)/2.0 + (range_high + range_low)/2.0;
		setUpBSplines(r,1);//set up Bsplines for the abscissas
		res = hamiltonianFunction( r , i , j , l_quantum );
		integral_result += weights[m]*(res[0]);
		integral_result_2 += weights[m]*(res[1]);
	}
	integral_result *= (range_high - range_low)/2.0;
	integral_result_2 *= (range_high - range_low)/2.0;
	return {integral_result,integral_result_2};
}


double Functions::gaussianQuadratureB(double range_low,double range_high ,int i , int j ,int l_quantum)
{

	vector<double> absciccas = { -0.3399810435848563,0.3399810435848563,-0.8611363115940526,0.8611363115940526};
	vector<double> weights   = {0.6521451548625461,0.6521451548625461,0.3478548451374538,0.3478548451374538};
	double integral_result   = 0.0;
	double r = 0.0;
	for (int m = 0; m < absciccas.size();m++)
	{
		r = absciccas[m]*(range_high - range_low)/2.0 + (range_high + range_low)/2.0;
		setUpBSplines(r,1);//set up Bsplines for the abscissas
		integral_result += weights[m]*(bsplineFunction(  i , j ));
	}
	integral_result *= (range_high - range_low)/2.0;
	return integral_result;
}

double Functions::gaussianQuadratureRho(double range_low,double range_high)
{
    vector<double> absciccas = { -0.3399810435848563,0.3399810435848563,-0.8611363115940526,0.8611363115940526};
	vector<double> weights   = {0.6521451548625461,0.6521451548625461,0.3478548451374538,0.3478548451374538};
	double integral_result   = 0.0;
	double r = 0.0;
	for (int m = 0; m < absciccas.size();m++)
	{
		r = absciccas[m]*(range_high - range_low)/2.0 + (range_high + range_low)/2.0;
		setUpBSplines(r,1);//set up Bsplines for the abscissas
		integral_result += weights[m]*pow(rho( r , Bspline_final),2.0);
    }
	integral_result *= (range_high - range_low)/2.0;

    return integral_result;
}

double Functions::Vee(int coeff,double r)
{
		double V_dir=0.0,V_ee=0.0,P_nl=0.0;
		if(P_nl_coeff.size()!=0)
		{
			V_dir = (1.0/r)*xi.dot(Bspline_final);
			V_ee = exchangePotential(rho(r,Bspline_final)) + V_dir;
		}
		if(P_nl_coeff.size()!=0)
		{
			cout<<"test "<<P_nl_coeff[coeff].size()<<" "<<Bspline_final.size()<<endl;
			P_nl = P_nl_coeff[coeff].dot(Bspline_final.segment(1,num_Bsplines-2));
		}
		return P_nl*V_ee*P_nl;
}

double Functions::gaussianQuadratureV(double range_low,double range_high ,int coeff)
{

	vector<double> absciccas = { -0.3399810435848563,0.3399810435848563,-0.8611363115940526,0.8611363115940526};
	vector<double> weights   = {0.6521451548625461,0.6521451548625461,0.3478548451374538,0.3478548451374538};
	double integral_result   = 0.0;
	double r = 0.0;
	/*for(int i = 0; i < grid.size();i++)
	{	
		rho_vec(i) = rho(grid(i),A_B.row(i));
	}
	xi = solveSystemCoefficients(rho_vec,Z);
	*/
	for (int m = 0; m < absciccas.size();m++)
	{
		r = absciccas[m]*(range_high - range_low)/2.0 + (range_high + range_low)/2.0;
		setUpBSplines(r,1);//set up Bsplines for the abscissas
		integral_result += weights[m]*(Vee( coeff , r));
	}
	integral_result *= (range_high - range_low)/2.0;
	return integral_result;
}

double Functions::hydrogen_s1(double r,double pi)
{
	double a_0 = 1.0;
	double Z = 1.0;
	double phi = 0.0;
	phi = (1/sqrt(pi))*pow(Z/a_0,1.5)*exp(-Z*r/a_0);
	return pow(phi,2);
}


double Functions::exchangePotential(double rho_r)
{   
    double pi = acos(-1);
    double V_exc;
    V_exc = - 3.0*pow((3.0*rho_r)/(e*8.0*pi),0.3333);
    return V_exc;
}

double Functions::rho(double r,VectorXd Bspline_final)
{
    double pi = acos(-1);
    double constant = (1.0/(4.0*pi))*e;
	double rho = 0.0;
    int N_j;
    for(int i = 0; i < N_occ_internal.size();i++)
    {
        N_j = N_occ_internal[i][0];
        rho += N_j* pow(P_nl_coeff[i].dot(Bspline_final.segment(1,P_nl_coeff[i].size()))/r,2.0);
    }
    return constant*rho;
} 