       #include <stdio.h>
       #include <math.h>

      const int n = 37;  // number of knotpoints (including ghost points)
      const int kord = 4;  // order of the B-splines
      const int imax = 500;  // number of points where we want the B-splines
      void setup(int n, int kord, double *t, double xx, double *Bavx, double *dBdx, double *dB2dx);
      int main()
      {
	double h,hh;
	double t[n],x[imax];        
        double Bavx[imax],dBdx[imax],dB2dx[imax];
        double xmax,xx;
	int i,j,k;
        int jchoice;
        FILE *result;
      

        x[0]=0.0;
        x[imax-1]=1.0;
        h=x[imax-1]/(n-2*(kord-1)-1);
        xmax=500.0;
        jchoice=1;  // print out this B-spline
       
/*      t holds the knot sequence; it does not need to be linear (as here)*/
        for (i = 0; i < kord; i = i + 1) {t[i]=0.0;}
        for (i = kord; i <= n-kord; i = i + 1) {t[i]=t[i-1]+h;}
        for (i = n-kord+1; i < n; i = i + 1) {t[i]=t[i-1];} 
        printf( "the knot sequence \n" );
        for(i=0; i < n; i = i + 1) {printf( "   %d  %f \n", i, t[i]);}

/*      x holds the grid on which we want the B-splines - we can chose to generate*/
/*      them in any point between t(first) and t(last) */
        
        
        hh= x[imax-1]/(xmax-1.0);
        for (i = 1; i < imax-1; i = i + 1) {x[i]=x[i-1]+hh;}

/*      write out the selected Bspline jchoice to file result.dat*/
        result = fopen("result.dat", "w");

/*      We loop over all the x-values and generate all Bsplines in this point - i.e. also those that are zero*/ 
        for (i = 0; i < imax; i = i + 1) 
	  {
          setup( n,  kord,  t,   x[i],  Bavx,  dBdx,  dB2dx);
          fprintf(result, "%f %f %f %f \n ",x[i],Bavx[jchoice-1],dBdx[jchoice-1],dB2dx[jchoice-1] );

	  }
         

        return 0;
      }
        void setup(int n, int kord, double *t,  double xx, double *Bavx, double *dBdx, double *dB2dx)
        {	  
        double B[n-kord+2][kord];
        double tol;

      	int i,j,k;  

            
      	/* The output arrays are put to zero*/
              tol=1e-10;
              for (j = 0; j < n-kord; j = j + 1)	 	  
      	    {
                       Bavx[j]=0.0;
                       dBdx[j]=0.0;
                       dB2dx[j]=0.0;
      	    }
               
             
                printf( "   calculate for the grid point %d  %f \n", i, xx);
                for (j = 0; j < n-kord+2; j = j + 1)
      	    {for (k= 0; k< kord; k = k + 1)	  
      	      {B[j][k]=0.0; }}
      	  
               
                   
      	/* k=1 */
              	         
                for (j = 0; j < n-kord; j = j + 1)
                          { 		
                          if( t[j] <= xx && xx < t[j+1] )
              	     {
                             B[j][0]=1.0;
                             /*printf("k=1  spline %d  %f  %f  %f \n", j+1,B[j][0],t[j],t[j+1]);*/
              	     }
              	    if(fabs(xx-t[n-1]) != 0.0 && j==n-kord -1)
              		    { B[j][0]=1.0;
                                    printf(" for k= %d B %d  = %f in last point %f \n",1,j+1, B[j][0],xx);
                                    }
      	    }
      	  /* k>1 */
                for (k = 1; k < kord; k = k + 1)
                  { 
                    for (j = 0; j < n-kord; j = j + 1)
       /* since both k and j starts at 0 we need to add +1 when both k and j is in the index of t */
      		/* this +1 is kept explicit for transparancy */
                     { 
      	        if(fabs(B[j][k-1]) >  tol)
      	        { 
                  B[j][k]=(xx-t[j])/(t[j+k-1+1]-t[j])*B[j][k-1];
      		  /* printf( "   spline j = %d k = %d %f %f %f \n", j+1,k+1,B[j][k],B[j][k-1],t[j+k-1+1]-t[j]);*/}	
                      if(fabs(B[j+1][k-1])>tol)
                      { 
                        B[j][k]= B[j][k]+ (t[j+k+1]-xx)/(t[j+k+1]-t[j+1])*B[j+1][k-1];
                        /*printf( "   spline j = %d k = %d %f %f %f \n", j+1,k+1,B[j][k],B[j+1][k-1],t[j+k+1]-t[j+1]);*/}
                     
                      if(fabs(xx-t[n-1]) < tol && j==n-kord -1)
      		    { B[j][k]=1.0;
                            printf(" for k= %d B %d  = %f in last point %f \n",k+1,j+1, B[j][k],xx);
                            }
      	      
      	      if(k==kord-1)
      		{ 
                         /* the Splines*/
                        Bavx[j]=B[j][k];
      		  /* if(fabs(xx-t[n-1]) < tol && j==n-kord-1)
      		     { Bavx[j]=1.0;}	 */
                        
                        /* the first order derivatives*/    
                        if(fabs(B[j][k-1]) >  tol)
      		    { dBdx[j]=dBdx[j]+(kord-1)*B[j][k-1]/(t[j+k-1+1]-t[j]);
      		    }
                        if(fabs(B[j+1][k-1]) >  tol)
      		    { dBdx[j]=dBdx[j]-(kord-1)*B[j+1][k-1]/(t[j+k+1]-t[j+1]);
      		    }

                         /* the 2nd order derivatives*/ 
                        if(fabs(B[j][k-2]) >  tol)  
      		    {
                dB2dx[j]=dB2dx[j]+(kord-1)*(kord-2)*B[j][k-2]/
      			 ((t[j+k-1+1]-t[j])*(t[j+k-2+1]-t[j]));
      		    } 
                        if(fabs(B[j+1][k-2]) >  tol)  
      		    {
                dB2dx[j]=dB2dx[j]-(kord-1)*(kord-2)*B[j+1][k-2]*(
      		     1.0/((t[j+k-1+1]-t[j])*(t[j+k-1+1]-t[j+1]))+
      		     1.0/((t[j+k+1]-t[j+1])*(t[j+k-1+1]-t[j+1])));

      		    } 
                        if(fabs(B[j+2][k-2]) >  tol)  
      		    {
                dB2dx[j]=dB2dx[j]+(kord-1)*(kord-2)*B[j+2][k-2]/
      			((t[j+k+1]-t[j+1])*(t[j+k+1]-t[j+2]));                   
      		    } 
      		  {printf( "   spline %d in %f  is %f  %f  %f \n", j+1,xx,Bavx[j],dBdx[j],dB2dx[j]);} 
                               
      		}
      	      }
      	    }
	}
        
        

