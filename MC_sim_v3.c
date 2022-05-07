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

#define BAR "####################################################################################################"
#define BARWIDTH 100

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

void wrapping(Pos * positions ,int N,double L,int i){
		/*wraps coordinate of particle i */
		if (positions ->pos_x[i]<0){
			positions ->pos_x[i] += L;
		}else if(positions ->pos_x[i]>L){
			positions ->pos_x[i] -= L;
		}

		if (positions ->pos_y[i]<0){
			positions ->pos_y[i] += L;
		}else if(positions ->pos_y[i]>L){
			positions ->pos_y[i] -= L;
		}
		
		if (positions ->pos_z[i]<0){
			positions ->pos_z[i] += L;
		}else if(positions ->pos_z[i]>L){
			positions ->pos_z[i] -= L;
		}	
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

void RDFSample(double *g,int N,double r,double dr,Pos *positions,double L,double V,double delg,double r_in){
	int index;
	if (r_in < 0.5*L){
		index = (int)(r_in/delg);
		g[index] += 2;			
	}
}


void energySample(Pos *positions,int N,double L,double*E_P,double* g,double r_samp,double dr_samp,double V,double delg,int *frames,int cycle){
	double dx=0,dy=0,dz=0;
	double dr_squared=0,dr_6=0,dr_12=0;
	double E_pot = 0;
	
	for(int i = 0;i<(N-1);i++){
		wrapping(positions,N,L,i);
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
			if (cycle>1e6){
				RDFSample(g,N,r_samp,dr_samp,positions,L,V,delg,sqrt(dr_squared));
			}

		}
	}
	if(cycle >1e6){
		(*frames) ++;
	}
	(*E_P) = E_pot;
	
}

double Rand(){
	//Returns a uniform random value in 0<x<1
	return (rand() + 1.0)/(RAND_MAX+2.0);
}



void trialMove(Pos*position,int N,double L,double *E_P,int *accept_cnt,double T,double del){
	/*Attempt to perform Trial move of random particle o with random displacement*/
	int o;
	double randNr,beta,enOld,enNew,ranX,ranY,ranZ;

	beta = 1/T;//1/kB*T
	particle part;

	o = round(Rand()*N);
	
	part.x = (position ->pos_x)[o];
	part.y = (position ->pos_y)[o];
	part.z = (position ->pos_z)[o];

	enOld = energyParticle(&part,position,o,N,L,E_P); 

	ranX = Rand();
	ranY = Rand();
	ranZ = Rand();
	//attempt to translate particle o
	part.x = (position->pos_x)[o]+ (ranX - 0.5)*del;
	part.y = (position->pos_y)[o]+ (ranY - 0.5)*del;
	part.z = (position->pos_z)[o]+ (ranZ - 0.5)*del;

	enNew = energyParticle(&part,position,o,N,L,E_P);
	//If condition fullfilled, accept move update particle o
	randNr = Rand();
	if (randNr < exp(-beta*(enNew- enOld))){
		(position ->pos_x)[o] = part.x;
		(position ->pos_y)[o] = part.y;
		(position ->pos_z)[o] = part.z;
		(*accept_cnt) ++;
	}
}

void WriteEnergy(FILE *fp,double E_pot){
	fprintf(fp,"%lf \n",E_pot);
}

void writeRDF(FILE *fp,double g){
	fprintf(fp,"%lf \n",g);
}


void loadingBar(int i, int cycles){
	double percent = ((double)i)/cycles;
	int val = (int)(percent *100);
	int lpad = (int)(percent*BARWIDTH);
	int rpad = BARWIDTH -lpad;
	printf("\rSimulation Progress: [%.*s%*s](%3d%%)",lpad,BAR,rpad,"",val);
	fflush(stdout);
}

int main(){
	/*Main function, initialize and run MC simulation*/
	int N,sample,cycles=0,done,i,accept_cnt=0;
	double dr,L,offset,E_P,E_P_sample;

	E_P = E_P_sample = 0;
	dr = 1.67;//initial displacement of particles in lattice
	offset = 0.5;//offset 

	
	initialData(&N,&sample,&cycles);//Ask user for simulation parameters

	Pos positions;
	allocateArray(&positions,N);//allocate memory for positions
	Initialize(&positions,N,dr,offset);



	L = dr*pow(N,0.3333) + offset;//Box length

	//============ FILE PATH ETC...==============
	FILE *write_path,*write_path_RDF;
	char * file_name_E = "E_MC.txt";
	char * file_name_RDF = "RDF_MC.txt";
	write_path = fopen(file_name_E,"w");
	


	//======RDF SHIT===============
	double pi = acos(-1);
	double dr_samp = 0.3; 
	double r_samp =0;
	double nid = 0;
	double vb = 0;
	double *g;

	int steps = 400,frames = 0;
	g = (double*)malloc(steps*sizeof(double));
	double delg = L/(2*steps);
	memset(g,0,sizeof(double));
	double R = pow(3/(4*pi),0.333)*L;
	double V = (4/3)*pi*R*R*R;
	double rho = N/V;

	double T = 1.4;
	double ratio_ideal = 0.20,ratio = 0.2;
	double del = 2.5;
	
	//================MAIN LOOP================
	i = done = 0;
	write_path_RDF = fopen(file_name_RDF,"w");
	while(!done){
		trialMove(&positions,N,L,&E_P,&accept_cnt,T,del);
		loadingBar(i,cycles);
		//ratio = (double)accept_cnt/(i+1);
		
		/*
		if (i>10000){
			if (ratio > 0.25){
				del *= 1.05;
			} 
			if( ratio <0.25){
				del *= 0.95;
			}
		} 	
		*/
		if (!(i%4000) && i>1e6){
			for (int i = 0; i<steps;i++){
				r_samp = delg*(0.5 +i);
				vb = ((i+1)*(i+1)*(i+1) -i*i*i)*delg*delg*delg;
				nid = (4/3)*pi*vb*rho;
				g[i] /= (40*N*nid); 
				writeRDF(write_path_RDF,g[i]);
			}
			fprintf(write_path_RDF,"\n");
		}
		
		if(!(i%sample)){
			energySample(&positions,N,L,&E_P_sample,g,r_samp,dr_samp,V,delg,&frames,i);
			WriteEnergy(write_path,E_P_sample);
		}
		if (i == cycles){
			done = 1;
		}
		i++;
	}
	printf("\n Accept ratio: %lf %lf \n", ((double)accept_cnt)/cycles,del);
	printf("delta g: %lf \n",delg);
	fclose(write_path_RDF);
	fclose(write_path);
	return 0;
}