#include <vector>
#include <complex>

class Potential{

private:
	double Pot_HO;
public:
	Potential(double test)
	{
		Pot_HO = test;
	}

	double PotentialHO(double z)
	{
		return (1.0/2)*z*z; 
	}
	double PotentialHO2(double z)
	{
		 return (1.0/2)*z*z + 4*exp(-4*z*z);
	}

};