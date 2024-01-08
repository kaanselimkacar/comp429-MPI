/* 
 * Solves the Panfilov model using an explicit numerical scheme.
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory 
 * and reimplementation by Scott B. Baden, UCSD
 * 
 * Modified and  restructured by Didem Unat, Koc University
 *
 */
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <math.h>
#include <sys/time.h>

//TODO: include mpi.h
#include <mpi.h>
using namespace std;


// Utilities
// 

// Timer
// Make successive calls and take a difference to get the elapsed time.
static const double kMicro = 1.0e-6;
double getTime()
{
    struct timeval TV;
    struct timezone TZ;

    const int RC = gettimeofday(&TV, &TZ);
    if(RC == -1) {
            cerr << "ERROR: Bad call to gettimeofday" << endl;
            return(-1);
    }

    return( ((double)TV.tv_sec) + kMicro * ((double)TV.tv_usec) );

}  // end getTime()

// Allocate a 2D array
double **alloc2D(int m,int n){
   double **E;
   int nx=n, ny=m;
   E = (double**)malloc(sizeof(double*)*ny + sizeof(double)*nx*ny);
   assert(E);
   int j;
   for(j=0;j<ny;j++) 
     E[j] = (double*)(E+ny) + j*nx;
   return(E);
}
    
// Reports statistics about the computation
// These values should not vary (except to within roundoff)
// when we use different numbers of  processes to solve the problem
 double stats(double **E, int m, int n, double *_mx){
     double mx = -1;
     double l2norm = 0;
     int i, j;
     for (j=1; j<=m; j++)
       for (i=1; i<=n; i++) {
	   l2norm += E[j][i]*E[j][i];
	   if (E[j][i] > mx)
	       mx = E[j][i];
      }
     *_mx = mx;
     l2norm /= (double) ((m)*(n));
     l2norm = sqrt(l2norm);
     return l2norm;
 }

// External functions
extern "C" {
    void splot(double **E, double T, int niter, int m, int n);
}
void cmdLine(int argc, char *argv[], double& T, int& n, int& px, int& py, int& plot_freq, int& no_comm, int&num_threads);


void simulate (double** E,  double** E_prev,double** R,
	       const double alpha, const int n, const int m, const double kk,
	       const double dt, const double a, const double epsilon,
	       const double M1,const double  M2, const double b)
{
  int i, j; 
    /* 
     * Copy data from boundary of the computational box 
     * to the padding region, set up for differencing
     * on the boundary of the computational box
     * Using mirror boundaries
     */
  
    for (j=1; j<=m; j++) 
      E_prev[j][0] = E_prev[j][2];
    // no need for a second loop, can combine loops
    for (j=1; j<=m; j++) 
      E_prev[j][n+1] = E_prev[j][n-1];
   
    // TODO: fix this for MPI programming
    
    /* 
    for (i=1; i<=n; i++) 
      E_prev[0][i] = E_prev[2][i];
    // no need for a second loop, can combine loops
    for (i=1; i<=n; i++) 
      E_prev[m+1][i] = E_prev[m-1][i];
    */

    // Solve for the excitation, the PDE
    for (j=1; j<=m; j++){
      for (i=1; i<=n; i++) {
	    E[j][i] = E_prev[j][i]+alpha*(E_prev[j][i+1]+E_prev[j][i-1]-4*E_prev[j][i]+E_prev[j+1][i]+E_prev[j-1][i]);
      }
    }
    
    /* 
     * Solve the ODE, advancing excitation and recovery to the
     *     next timtestep
     */
    for (j=1; j<=m; j++){
      for (i=1; i<=n; i++)
	    E[j][i] = E[j][i] -dt*(kk* E[j][i]*(E[j][i] - a)*(E[j][i]-1)+ E[j][i] *R[j][i]);
    }
    
    for (j=1; j<=m; j++){
      for (i=1; i<=n; i++)
	    R[j][i] = R[j][i] + dt*(epsilon+M1* R[j][i]/( E[j][i]+M2))*(-R[j][i]-kk* E[j][i]*(E[j][i]-b-1));
    }
    
}

// Main program
int main (int argc, char** argv)
{
  /*
   *  Solution arrays
   *   E is the "Excitation" variable, a voltage
   *   R is the "Recovery" variable
   *   E_prev is the Excitation variable for the previous timestep,
   *      and is used in time integration
   */
  double **E, **R, **E_prev;
  
  // Various constants - these definitions shouldn't change
  const double a=0.1, b=0.1, kk=8.0, M1= 0.07, M2=0.3, epsilon=0.01, d=5e-5;
  
  double T=1000.0;
  int m=200,n=200;
  int plot_freq = 0;
  int px = 1, py = 1;
  int no_comm = 0;
  int num_threads=1; 
  
  int P = 1, myrank = 0; 
  int tag1 = 1, tag2 = 2;
     
  /* Initializations */
  MPI_Init(&argc, &argv);
  
  
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  cmdLine( argc, argv, T, n,px, py, plot_freq, no_comm, num_threads);
  m = n;  

  
  MPI_Request reqs[4];
  MPI_Request reqs_2[2];
  MPI_Status mpi_stats[4];
  MPI_Status mpi_stats_2[2];
  // matrix size N*N
  // for part 1 each processor should work on N/P rows 
  // TODO: add local variables myE, myE_prev, my_R
  double **myE, **myR, **myE_prev;
  // Allocation
  int myRowSize = m/P;
  myE = alloc2D(myRowSize+2 , n+2); 
  myE_prev = alloc2D(myRowSize+2 , n+2); 
  myR = alloc2D(myRowSize+2 , n+2); 

  int i,j;
  // Initialization
  for (j=1; j<=myRowSize; j++)
    for (i=1; i<=n; i++)
      myE_prev[j][i] = myR[j][i] = 0;
  
  for (j=1; j<=myRowSize; j++)
    for (i=n/2+1; i<=n; i++)
      myE_prev[j][i] = 1.0;
  
  for (j=myRowSize/2+1; j<=myRowSize; j++)
    for (i=1; i<=n; i++)
      myR[j][i] = 1.0;


  // Allocate contiguous memory for solution arrays
  // The computational box is defined on [1:m+1,1:n+1]
  // We pad the arrays in order to facilitate differencing on the 
  // boundaries of the computation box
  E = alloc2D(m+2,n+2);
  E_prev = alloc2D(m+2,n+2);
  R = alloc2D(m+2,n+2);
  
  // Initialization
  for (j=1; j<=m; j++)
    for (i=1; i<=n; i++)
      E_prev[j][i] = R[j][i] = 0;
  
  for (j=1; j<=m; j++)
    for (i=n/2+1; i<=n; i++)
      E_prev[j][i] = 1.0;
  
  for (j=m/2+1; j<=m; j++)
    for (i=1; i<=n; i++)
      R[j][i] = 1.0;
  
  
  double dx = 1.0/n;

  // For time integration, these values shouldn't change 
  double rp= kk*(b+1)*(b+1)/4;
  double dte=(dx*dx)/(d*4+((dx*dx))*(rp+kk));
  double dtr=1/(epsilon+((M1/M2)*rp));
  double dt = (dte<dtr) ? 0.95*dte : 0.95*dtr;
  double alpha = d*dt/(dx*dx);

  if (myrank == 0)
  {
    cout << "Num Processes   : " << P << endl;
    cout << "Grid Size       : " << n << endl; 
    cout << "Duration of Sim : " << T << endl; 
    cout << "Time step dt    : " << dt << endl; 
    cout << "Process geometry: " << px << " x " << py << endl;
    if (no_comm)
      cout << "Communication   : DISABLED" << endl;
  
    cout << endl;
  }
  // Start the timer
  double t0 = getTime();
  
 
  // Simulated time is different from the integer timestep number
  // Simulated time
  double t = 0.0;
  // Integer timestep number
  int niter=0;

  // if single process
 double **tempE;    
 double **tempE_prev;    
 double **tempR;    
 if (P == 1){
    tempE = E;
    tempE_prev = E_prev;
    tempR = R;
    E = myE;
    E_prev = myE_prev;
    R = myR;
  }
  
  cout << "myRowSize = " << myRowSize << " my n = " << n << endl;
  while (t<T) {
    t += dt;
    niter++;
    
    //TODO: add communication
    // send north -- tag1
    // send south -- tag2
    // receive north -- tag2
    // receive south -- tag1
    if (myrank == 0 && P != 1){
        // send and receive data from south, mirror for north
        
        // receive from south
        MPI_Irecv(&myE_prev[myRowSize+1][1], n, MPI_DOUBLE, myrank+1, tag1, MPI_COMM_WORLD, &reqs[0]);

        // send to south
        MPI_Isend(&myE_prev[myRowSize][1], n, MPI_DOUBLE, myrank+1, tag2, MPI_COMM_WORLD, &reqs[1]);
        // mirror for north    
        for (i=1; i<=n; i++) 
          myE_prev[0][i] = myE_prev[2][i];
    }
    else if (myrank == P-1 && P != 1){
        // send and receive data from north, mirror for south
        
        // receive from north
        MPI_Irecv(&myE_prev[0][1], n, MPI_DOUBLE, myrank-1, tag2, MPI_COMM_WORLD, &reqs[0]);

        // send to north
        MPI_Isend(&myE_prev[1][1], n, MPI_DOUBLE, myrank-1, tag1, MPI_COMM_WORLD, &reqs[1]);
        
        // mirror for south
        for (i=1; i<=n; i++) 
          myE_prev[myRowSize+1][i] = myE_prev[myRowSize-1][i];
    }
    else if (P != 1){
        // send and receive data both from north and south
        
        // receive from north
        MPI_Irecv(&myE_prev[0][1], n, MPI_DOUBLE, myrank-1, tag2, MPI_COMM_WORLD, &reqs[0]);
        // receive from south
        MPI_Irecv(&myE_prev[myRowSize+1][1], n, MPI_DOUBLE, myrank+1, tag1, MPI_COMM_WORLD, &reqs[1]);
    
        // send to north
        MPI_Isend(&myE_prev[1][1], n, MPI_DOUBLE, myrank-1, tag1, MPI_COMM_WORLD, &reqs[2]);
        // send to south
        MPI_Isend(&myE_prev[myRowSize][1], n, MPI_DOUBLE, myrank+1, tag2, MPI_COMM_WORLD, &reqs[3]);
        
    }
    
    if ( (myrank == 0 || myrank == P -1) && P != 1){
      MPI_Waitall(2, reqs_2, mpi_stats_2);    
    }
    else if (P != 1) {
      MPI_Waitall(4, reqs, mpi_stats);    
    }
   
    // for single process 
    if (P == 1){
      for (i=1; i<=n; i++) 
        myE_prev[0][i] = myE_prev[2][i];
      // no need for a second loop, can combine loops
      for (i=1; i<=n; i++) 
        myE_prev[m+1][i] = myE_prev[m-1][i];
    }

    simulate(myE, myE_prev, myR, alpha, n, myRowSize, kk, dt, a, epsilon, M1, M2, b); 
    
    //simulate(E, E_prev, R, alpha, n, m, kk, dt, a, epsilon, M1, M2, b); 
    
    //swap current E with previous E
    //double **tmp = E; E = E_prev; E_prev = tmp;
   
    double **tmp2 = myE; myE = myE_prev; myE_prev = tmp2; 
    if (plot_freq){
      // TODO: Gather data from all processes
      if (P != 1){
        MPI_Gather(&myE[1][0], (n+2) * (myRowSize), MPI_DOUBLE, &E[1][0], n*(n+2), MPI_DOUBLE, 0, MPI_COMM_WORLD ); 
      }
      if (myrank == 0 && P != 1){
        int k = (int)(t/plot_freq);
        if ((t - k * plot_freq) < dt){
	        splot(E,t,niter,m+2,n+2);
        }
      }
      else if (P == 1){
        int k = (int)(t/plot_freq);
        if ((t - k * plot_freq) < dt){
	        splot(myE,t,niter,m+2,n+2);
        }
      }
    }
  }//end of while loop

  // TODO: Gather data from all processes
  if (P != 1){
    MPI_Gather(&myE[1][0], (n+2) * myRowSize, MPI_DOUBLE, &E[1][0], n*(n+2), MPI_DOUBLE, 0, MPI_COMM_WORLD ); 
  }
  if (myrank == 0)
  {  
    double time_elapsed = getTime() - t0;

    double Gflops = (double)(niter * (1E-9 * n * n ) * 28.0) / time_elapsed ;
    double BW = (double)(niter * 1E-9 * (n * n * sizeof(double) * 4.0  ))/time_elapsed;

    cout << "Number of Iterations        : " << niter << endl;
    cout << "Elapsed Time (sec)          : " << time_elapsed << endl;
    cout << "Sustained Gflops Rate       : " << Gflops << endl; 
    cout << "Sustained Bandwidth (GB/sec): " << BW << endl << endl; 

    double mx;
    double l2norm = stats(E_prev,m,n,&mx);
    cout << "Max: " << mx <<  " L2norm: "<< l2norm << endl;

    if (plot_freq){
      cout << "\n\nEnter any input to close the program and the plot..." << endl;
      getchar();
    }
  }
  
  free (E);
  free (E_prev);
  free (R);
  
  if (P == 1){
    free (tempE);
    free (tempE_prev);
    free (tempR);
  }
  else{
    free (myE);
    free (myE_prev);
    free (myR);
  }
  MPI_Finalize(); 
  return 0;
}
