/*
Author:Markus Tollet
Last changed: 06/03-19
Monte-Carlo simulation 1.0 

Run program using: gcc MC_sim_v1.c -Wall -I...(path to folder)/Gnuwin32/include -L...(path to folder)/GnuWin32/lib -lm -lgsl - lgslcblas -o outputfile_name
Assuming gsl is installed and located in folder:GnuWin32
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

typedef struct positions{
	/*
	Store x,y,z components of r for all particles
	*/
	double *pos_x;
	double *pos_y;
	double *pos_z;
}Pos;

typedef struct  posPar{
	double x;
	double y;
	double z;
}particle;

void Initialize(Pos *positions,int N,double dr,double offset){
	//Initialize positions for IC being lattice like, velocities for given IC
	int N_cube = 2;
	int x_pos=0,y_pos=0,z_pos=0;

	//x_pos = y_pos = z_pos = 0;

	//Creates smallest perfect cube that fit N particles

	for (;N_cube*N_cube*N_cube<N;N_cube++){}


	int i = 0;

	while (i<N){
		//positions ->pos_x[i] = 1; 
		//positions ->pos_y[i] = 2;
		//positions ->pos_z[i] = 3;

		
		(positions ->pos_x)[i] = dr*x_pos + offset; 
		(positions ->pos_y)[i] = dr*y_pos + offset;
		(positions ->pos_z)[i] = dr*z_pos + offset;
		
		//printf("test1: %lf %lf %lf %d\n",positions ->pos_x[i],positions ->pos_y[i],positions ->pos_z[i],i);
		x_pos ++;
		if(x_pos == N_cube){//x row filled
			x_pos =0;
			y_pos ++;
		}
		if(y_pos == N_cube){//y row filled
			y_pos =0;
			z_pos ++;
		}
		i++;
	}	

}

void initialData(int *N,int *sample,int *cycles){
	/*
	Function:
	Read simulation parameters
	*/
	int done = 0;
	int i = 0;
	int res=0;
	char clr;
	printf("Welcome to Monte-Carlo simulation of Ar-fluid!\n");
	while(!done){
		if (i==0){
			printf("Nr of paricles:\n");
			res =fscanf(stdin,"%d",N);
		}
		if(i==1){
			printf("Sample frequency?:\n");
			res = fscanf(stdin,"%d",sample);
		}
		if(i==2){
			printf("Nr of Monte-Carlo cycles:\n");
			res = fscanf(stdin,"%d",cycles);
		}
		if(i==3){
			done = 1;//exits loop
		}
		if (res !=0){
			i++;
		}
		//Loop to clear buffer in case of bad input
		while ((clr = getchar()) != '\n' && clr != EOF) {}
	}
}

void allocateArray(Pos *positions,int N){
	positions ->pos_x = (double *)malloc(N*sizeof(double));
	positions ->pos_y = (double *)malloc(N*sizeof(double));
	positions ->pos_z = (double *)malloc(N*sizeof(double));
}

double energyParticle(particle * part,Pos*positions,int index,int N,double L,double*E_P ){
	double dx=0,dy=0,dz=0;
	double dr_squared=0,dr_6=0,dr_12=0;
	double E_pot = 0;

	for (int i = 0; i<N;i++){
		if (i != index){
			dx = part->x -(positions->pos_x[i]);
			dy = part->y -(positions->pos_y[i]);
			dz = part->z -(positions->pos_z[i]);
			//Minimum Image Convention
			dx -= L*round(dx/L);
			dy -= L*round(dy/L);
			dz -= L*round(dz/L);

			dr_squared = dx*dx + dy*dy + dz*dz;
			dr_6 = 1/(dr_squared*dr_squared*dr_squared);
			dr_12 = dr_6*dr_6;

			E_pot += 4*(dr_12-dr_6);
		}
	}
	return E_pot;
}

void energySample(Pos *positions,int N,double L,double*E_P){
	double dx=0,dy=0,dz=0;
	double dr_squared=0,dr_6=0,dr_12=0;
	double E_pot = 0;


	for(int i = 0;i<(N-1);i++){
		for(int j = i+1;j<N;j++){
			dx = (positions->pos_x[i])-(positions->pos_x[j]);
			dy = (positions->pos_y[i])-(positions->pos_y[j]);
			dz = (positions->pos_z[i])-(positions->pos_z[j]);

			//Minimum Image Convention
			dx -= L*round(dx/L);
			dy -= L*round(dy/L);
			dz -= L*round(dz/L);

			dr_squared = dx*dx + dy*dy + dz*dz;
			dr_6 = 1/(dr_squared*dr_squared*dr_squared);
			dr_12 = dr_6*dr_6;

			E_pot += 4*(dr_12-dr_6);
		}
	}
	(*E_P) = E_pot;
	
}

double Rand(){
	//Returns a uniform random value in 0<x<1
	return (rand() + 1.0)/(RAND_MAX+2.0);
}



void trialMove(Pos*position,int N,double L,double *E_P){
	/*Attempt to perform Trial move of random particle o with random displacement*/
	int o;
	double dr,randNr,beta,enOld,enNew;
	
	dr = 0.3;
	beta = 1;
	particle part;

	o = (int)(Rand()*N) +1;
	
	part.x = (position ->pos_x)[o];
	part.y = (position ->pos_y)[o];
	part.z = (position ->pos_z)[o];

	enOld = energyParticle(&part,position,o,N,L,E_P); 

	randNr = Rand();
	//attempt to translate particle o
	part.x = (position->pos_x)[o]+ (randNr - 0.5)*dr;
	part.y = (position->pos_y)[o]+ (randNr - 0.5)*dr;
	part.z = (position->pos_z)[o]+ (randNr - 0.5)*dr;

	enNew = energyParticle(&part,position,o,N,L,E_P);
	//If condition fullfilled, accept move update particle o
	if (randNr < exp(-beta*(enNew- enOld))){
		(position ->pos_x)[o] = part.x;
		(position ->pos_y)[o] = part.y;
		(position ->pos_z)[o] = part.z;
	}

}

void WriteEnergy(FILE *fp,double E_pot){
	fprintf(fp,"%lf \n",E_pot);
}


int main(){
	/*Main function, initialize and run MC simulation*/
	int N,sample,cycles,done,i;
	double dr,L,offset,E_P,E_P_sample;

	E_P = E_P_sample = 0;
	dr = 1.67;//initial displacement of particles in lattice
	offset = 0.5;//offset 

	
	initialData(&N,&sample,&cycles);//Ask user for simulation parameters

	Pos positions;
	allocateArray(&positions,N);//allocate memory for positions
	Initialize(&positions,N,dr,offset);



	L = dr*pow(N,0.3333) + offset;//Box length

	//================MAIN LOOP================
	FILE *write_path;
	char * file_name_E = "E_MC.txt";
	write_path = fopen(file_name_E,"w");

	i = done = 0;
	while(!done){
		trialMove(&positions,N,L,&E_P); 
		if(!(i%sample)){
			energySample(&positions,N,L,&E_P_sample);
			WriteEnergy(write_path,E_P_sample);
			printf("E_P sample: %lf \n",E_P_sample);
		}
		if (i == cycles-1){
			done = 1;
		}
		i++;
	}
	fclose(write_path);
	return 0;
}