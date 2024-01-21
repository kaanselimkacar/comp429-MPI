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

// include mpi.h
#include <mpi.h>
// include openmp
#include <omp.h>
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
       //if (E[j][i] > 1 || E[j][i] < 0)
        //cout << "value of " << E[j][i] << " at row = " << j << " col = " << i << endl; 
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
    /*
    #pragma omp parallel for
    for (j=1; j<=m; j++){ 
      E_prev[j][0] = E_prev[j][2];
      E_prev[j][n+1] = E_prev[j][n-1];
    }
    */
    /*
    for (j=1; j<=m; j++) 
      E_prev[j][0] = E_prev[j][2];
    // no need for a second loop, can combine loops
    for (j=1; j<=m; j++) 
      E_prev[j][n+1] = E_prev[j][n-1];
    */
    // fix this for MPI programming
    
    /* 
    for (i=1; i<=n; i++) 
      E_prev[0][i] = E_prev[2][i];
    // no need for a second loop, can combine loops
    for (i=1; i<=n; i++) 
      E_prev[m+1][i] = E_prev[m-1][i];
    */
    //#pragma omp parallel
    {

      // Solve for the excitation, the PDE
      //#pragma omp for collapse(2)
      for (j=1; j<=m; j++){
        for (i=1; i<=n; i++) {
	       E[j][i] = E_prev[j][i]+alpha*(E_prev[j][i+1]+E_prev[j][i-1]-4*E_prev[j][i]+E_prev[j+1][i]+E_prev[j-1][i]);
          }
        }
    
        /* 
        * Solve the ODE, advancing excitation and recovery to the
        *     next timtestep
        */
        //#pragma omp for collapse(2)  
        for (j=1; j<=m; j++){
          for (i=1; i<=n; i++){
	          E[j][i] = E[j][i] -dt*(kk* E[j][i]*(E[j][i] - a)*(E[j][i]-1)+ E[j][i] *R[j][i]);
	          R[j][i] = R[j][i] + dt*(epsilon+M1* R[j][i]/( E[j][i]+M2))*(-R[j][i]-kk* E[j][i]*(E[j][i]-b-1));
          }
        }
    
      }// pragma omp parallel
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
  MPI_Init(&argc, &argv);
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
     
  /* Initializations */
  
  
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  cmdLine( argc, argv, T, n,px, py, plot_freq, no_comm, num_threads);
  m = n;  

  //TODO: add support for when the process number does not divide n
  // possible solution: add + 1 to myRowSize, and check for matrix size
  // if rowIndex * myrank * (n+2) < n+2 * n+2
  MPI_Request reqs[8] = {MPI_REQUEST_NULL,
                          MPI_REQUEST_NULL,
                           MPI_REQUEST_NULL,
                            MPI_REQUEST_NULL,
                             MPI_REQUEST_NULL,
                              MPI_REQUEST_NULL,
                               MPI_REQUEST_NULL,
                                MPI_REQUEST_NULL};  
  MPI_Status mpi_stats[8];  

  // matrix size N*N
  // for part 1 each processor should work on N/P rows 
  // add local variables myE, myE_prev, my_R
  double **myE, **myR, **myE_prev;
  


  // Allocate contiguous memory for solution arrays
  // The computational box is defined on [1:m+1,1:n+1]
  // We pad the arrays in order to facilitate differencing on the 
  // boundaries of the computation box
  //cout << "everyday I'm mallocing myrank = " << myrank << endl;
  E = alloc2D(m+2,n+2);
  E_prev = alloc2D(m+2,n+2);
  R = alloc2D(m+2,n+2);
  //cout << "finished I'm mallocing myrank = " << myrank << endl;
 
  // create temp variables to free them later
  double *tempE = &E[0][0];
  double *tempE_prev = &E_prev[0][0];
  double *tempR = &R[0][0];
 
  // Initialization
  int i,j;
  for (j=1; j<=m; j++)
    for (i=1; i<=n; i++)
      E_prev[j][i] = R[j][i] = 0;
  
  for (j=1; j<=m; j++)
    for (i=n/2+1; i<=n; i++)
      E_prev[j][i] = 1.0;
  
  for (j=m/2+1; j<=m; j++)
    for (i=1; i<=n; i++)
      R[j][i] = 1.0;
  
  // Allocation
  int myX = myrank % px;
  int myY = myrank / px;
  

  int usualRowSize = m / py;
  usualRowSize = (py * usualRowSize == m) ? usualRowSize : (usualRowSize + 1);
  int usualColSize = n / px;
  usualColSize = (px * usualColSize == n) ? usualColSize : (usualColSize + 1);
  int smallRowSize = m - (py-1) * usualRowSize;
  int smallColSize = n - (px-1) * usualColSize;
  //cout << "usualColSize = " << usualColSize << " usualRowSize = "<< usualRowSize << endl;  

  int myRowSize = m, myColSize = n;
  if (myY == (py - 1)){
    myRowSize = smallRowSize; 
  }  
  else {
    myRowSize = usualRowSize;
  }
   
  if (myX == (px - 1)){
    myColSize = smallColSize; 
  }  
  else {
    myColSize = usualColSize;
  }
  //cout << "everyday I'm mallocing myrank = " << myrank << endl;
  myE = alloc2D(myRowSize+2 , myColSize+2); 
  myE_prev = alloc2D(myRowSize+2 , myColSize+2); 
  myR = alloc2D(myRowSize+2 , myColSize+2); 
  //cout << "finished I'm mallocing myrank = " << myrank << endl;
  
  /* 6* 6, lets say P = 3, then myrank e {0,1,2}, myRowSize = 2
 *  for 0 it should be j*1 + 0   j + myrank * myRowSize
 *  for 1 it should be j*1 + 2
 *  for 2 it should be j*1 + 4 
  12345678 
  12345678 -- 0
  12345678 -- 0
  12345678 -- 1
  12345678 -- 1
  12345678 -- 2
  12345678 -- 2
  12345678
  */
  // Initialization of myE_prev and myR
  
  for (j=1; j<=myRowSize; j++){
    for (i=1; i<=myColSize; i++){
      myE_prev[j][i] = E_prev[myY * usualRowSize + j][myX * usualColSize + i];
      myR[j][i]      =      R[myY * usualRowSize + j][myX * usualColSize + i];
    }
  }
  //cout << " new myRowSize = " << myRowSize << " myrank= " << myrank << " myY = " << myY << " myX = " << myX << endl; 
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
  
  //cout << "before barrier myrank = " << myrank << endl;
  //MPI_Barrier(MPI_COMM_WORLD);
  //cout << "after barrier myrank = " << myrank << endl;
  // Start the timer
  double t0 = getTime();
  
 
  // Simulated time is different from the integer timestep number
  // Simulated time
  double t = 0.0;
  // Integer timestep number
  int niter=0;
  
  int tag1 = 1, tag2 = 2, tag3 = 3, tag4 = 4;
  //cout << "everyday I'm mallocing myrank = " << myrank << endl;
  double *recvBuffWest = (double *) malloc((sizeof(double) * (myRowSize + 2)));;
  double *recvBuffEast = (double *) malloc((sizeof(double) * (myRowSize + 2)));;
  double *sendBuffWest = (double *) malloc((sizeof(double) * (myRowSize + 2)));;
  double *sendBuffEast = (double *) malloc((sizeof(double) * (myRowSize + 2)));;
  
  // init recv arrays
  for (i = 1 ; i <= myRowSize; i++) {
      recvBuffWest[i] = 0;
      recvBuffEast[i] = 0;
  }
  
  //cout << "usualRowSize = " << usualRowSize << " smallRowSize = " << smallRowSize << "myY = " << myY << endl;
  
  //cout << "on dead homies myX = " << myX << "  myY = " << myY << " myrank = " << myrank << " myRowSize = " << myRowSize << " myColSize = " << myColSize << endl;
  // no errors, but getting false L2 and Linf values
  // probable cause: communication with east and west at the same time 
  //cout << "Still alive !!  myrank = " << myrank << endl;
  while (t<T) {
    t += dt;
    niter++;
    
    // send north    -- tag1
    // send south    -- tag2
    // receive north -- tag2
    // receive south -- tag1
    //
    // send east     -- tag3
    // send west     -- tag4
    // receive east  -- tag4
    // receive west  -- tag3
   
    // ALLAH BU MILLETE BU KODU BI DAHA YAZDIRMASIN

    if (myY > 0){
      // send/recv north
      
      // receive from north
      MPI_Irecv(&myE_prev[0][1], myColSize, MPI_DOUBLE, myrank - px, tag2, MPI_COMM_WORLD, &reqs[0]);

      // send to north
      MPI_Isend(&myE_prev[1][1], myColSize, MPI_DOUBLE, myrank - px, tag1, MPI_COMM_WORLD, &reqs[1]);
    }
    
    if (myY < (py - 1)){
      // send/recv south
      
      // receive from south
      MPI_Irecv(&myE_prev[myRowSize+1][1], myColSize, MPI_DOUBLE, myrank + px, tag1, MPI_COMM_WORLD, &reqs[2]);

      // send to south
      MPI_Isend(&myE_prev[myRowSize][1], myColSize, MPI_DOUBLE, myrank + px, tag2, MPI_COMM_WORLD, &reqs[3]);
    }
    
    if (myX > 0){
      // send/recv from west
      // recv from west
      MPI_Irecv(&recvBuffWest[1], myRowSize, MPI_DOUBLE, myrank - 1, tag3, MPI_COMM_WORLD, &reqs[4]);
      // pack the message first
      for (i = 1; i <= myRowSize; i++)
        sendBuffWest[i] = myE_prev[i][1];
      // send to west
      MPI_Isend(&sendBuffWest[1], myRowSize, MPI_DOUBLE, myrank - 1, tag4, MPI_COMM_WORLD, &reqs[5]);
    }
    
    if (myX < (px - 1)){
      // send/recv from east
      // recv from east
      MPI_Irecv(&recvBuffEast[1], myRowSize, MPI_DOUBLE, myrank + 1, tag4, MPI_COMM_WORLD, &reqs[6]);
      // pack the message first
      for (i = 1; i <= myRowSize; i++)
        sendBuffEast[i] = myE_prev[i][myColSize];
      // send to east
      MPI_Isend(&sendBuffEast[1], myRowSize, MPI_DOUBLE, myrank + 1, tag3, MPI_COMM_WORLD, &reqs[7]);
    }

    // mirror north ----- processes with myY = 0
    // mirror south ----- processes with myY = py-1
    // mirror east  ----- myX = px -1
    // mirror west  ----- myX = 0
    if (myY  == 0){
      // mirror north
      for (i = 1; i <= myColSize; i++)
        myE_prev[0][i] = myE_prev[2][i];
    }
     
    if (myY  == py - 1){
      // mirror south
      for (i = 1; i <= myColSize; i++)
        myE_prev[myRowSize + 1][i] = myE_prev[myRowSize - 1][i];
    }

    if (myX  == 0){
      // mirror west
      for (i = 1; i <= myRowSize; i++)
        myE_prev[i][0] = myE_prev[i][2];
    }
    
    if (myX  == px - 1){
      // mirror east
      for (i = 1; i <= myRowSize; i++)
        myE_prev[i][myColSize + 1] = myE_prev[i][myColSize - 1];
    }
   
    //cout << "Waiting! t = " << t << "  myrank = " << myrank << endl;
    // wait for Isends and Irecvs
    for (i = 0; i < 8; i++){
        if (reqs[i] == MPI_REQUEST_NULL)
            continue;
       
        //cout << "waiting! myrank = " << myrank << endl;
        MPI_Wait(&reqs[i], &mpi_stats[i]);
        reqs[i] = MPI_REQUEST_NULL;
        //cout << "waited! myrank = " << myrank << endl;
        // debug
        //int count;
        //MPI_Get_count(&mpi_stats[i], MPI_DOUBLE, &count);
        //cout << "myrank = " << myrank << " received " << count << " elements with mpi_tag = " << mpi_stats[i].MPI_TAG << endl; 
    }
    //cout << "Waited! t = " << t << "  myrank = " << myrank << endl;

    // unpack the messages
    if (myX > 0){
      // unpack message for west
      for (i = 1; i <= myRowSize; i++)
        myE_prev[i][0] = recvBuffWest[i];
    } 
    if (myX < (px - 1)){
      // unpack message for east
      for (i = 1; i <= myRowSize; i++)
        myE_prev[i][myColSize+1] = recvBuffEast[i];
    }
 
    //cout << "Simulating t = " << t << "  myrank = " << myrank << endl;
    simulate(myE, myE_prev, myR, alpha, myColSize, myRowSize, kk, dt, a, epsilon, M1, M2, b); 
   
    double **tmp2 = myE; myE = myE_prev; myE_prev = tmp2; 
    if (plot_freq){
      // Gather data from all processes
      if (P != 1){
        // TODO: fix for 2d
        //MPI_Gather(&myE[1][0], (n+2) * (myRowSize), MPI_DOUBLE, &E[1][0], myRowSize*(n+2), MPI_DOUBLE, 0, MPI_COMM_WORLD ); 
       
    // copy the data of rank 0 to E
    if (myrank == 0){
      for (i = 1; i <= myRowSize; i++){
        for(j = 1; j <= myColSize; j++){
          E[i][j] = myE[i][j];
        }
      }
    }
    // each process sends their data to rank 0
    if (myrank == 0)
    {
        //double *recvBuffer = (double *) malloc(sizeof(double) * usualRowSize * (usualColSize + 2) + 2);
      double **recvBuffer = alloc2D(usualRowSize + 2, usualColSize + 2);
      for (int rank = 1; rank < P; rank++)
      {
        int senderX = rank % px , senderY = rank / px;
        if (senderX == (px - 1) && senderY == (py - 1)){
          MPI_Recv(&recvBuffer[1][1], (smallColSize+2) * smallRowSize, MPI_DOUBLE, rank, rank, MPI_COMM_WORLD, &mpi_stats[0]); 
          // cook
          for(i = 1; i <= smallRowSize; i++){
            for(j = 1; j <= smallColSize; j++){
              E[senderY * usualRowSize + i][senderX * usualColSize + j] = recvBuffer[i][j];
            }
          }
        }
        else if (senderX == (px - 1)){
          MPI_Recv(&recvBuffer[1][1], (smallColSize+2) * usualRowSize, MPI_DOUBLE, rank, rank, MPI_COMM_WORLD, &mpi_stats[0]); 
          for(i = 1; i <= usualRowSize; i++){
            for(j = 1; j <= smallColSize; j++){
              E[senderY * usualRowSize + i][senderX * usualColSize + j] = recvBuffer[i][j];
            }
          }
            
        }
        else if (senderY == (py - 1)){
          MPI_Recv(&recvBuffer[1][1], (usualColSize+2) * smallRowSize, MPI_DOUBLE, rank, rank, MPI_COMM_WORLD, &mpi_stats[0]); 
          for(i = 1; i <= smallRowSize; i++){
            for(j = 1; j <= usualColSize; j++){
              E[senderY * usualRowSize + i][senderX * usualColSize + j] = recvBuffer[i][j];
            }
          }
           
        }
        else {
          MPI_Recv(&recvBuffer[1][1], (usualColSize+2) * usualRowSize, MPI_DOUBLE, rank, rank, MPI_COMM_WORLD, &mpi_stats[0]); 
          for(i = 1; i <= usualRowSize; i++){
            for(j = 1; j <= usualColSize; j++){
              E[senderY * usualRowSize + i][senderX * usualColSize + j] = recvBuffer[i][j];
            }
          }

        }
      } // end of for
     free(recvBuffer); 
    }// end of if   
    else{
        MPI_Send(&myE[1][1], myRowSize * (myColSize + 2), MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD);
        }
      if (myrank == 0){ 
        int k = (int)(t/plot_freq);
        if ((t - k * plot_freq) < dt){
	      splot(E,t,niter,m+2,n+2);
        }
      }
    } // end of if P!=1
      else{
        int k = (int)(t/plot_freq);
        if ((t - k * plot_freq) < dt){
	      splot(myE,t,niter,m+2,n+2);
        }
        
      }
    MPI_Barrier(MPI_COMM_WORLD);
    }// end of if plot_freq
  
  }//end of while loop
  
  //cout << "Still alive !! Still alive!!  myrank = " << myrank << endl;
  // TODO: fix for 2d
  if (P != 1){
    //MPI_Gather(&myE_prev[1][0], (n+2) * myRowSize, MPI_DOUBLE, &E_prev[1][0], myRowSize*(n+2), MPI_DOUBLE, 0, MPI_COMM_WORLD );
    
    /*
    //debug
    //print of elements of myE_prev of each process
    for (int k = 0; k < P ; k++){
      if (myrank == k){
        for (i = 1 ; i <= myRowSize; i++){
          for (j = 1; j <= myColSize; j++)
            cout << myE_prev[i][j];
          cout << endl;
        }      
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    */
    // copy the data of rank 0 to E
    if (myrank == 0){
      for (i = 1; i <= myRowSize; i++){
        for(j = 1; j <= myColSize; j++){
          E_prev[i][j] = myE_prev[i][j];
        }
      }
    }
    // each process sends their data to rank 0
    if (myrank == 0)
    {
        //double *recvBuffer = (double *) malloc(sizeof(double) * usualRowSize * (usualColSize + 2) + 2);
      double **recvBuffer = alloc2D(usualRowSize + 2, usualColSize + 2);
      for (int rank = 1; rank < P; rank++)
      {
        int senderX = rank % px , senderY = rank / px;
        if (senderX == (px - 1) && senderY == (py - 1)){
          MPI_Recv(&recvBuffer[1][1], (smallColSize+2) * smallRowSize, MPI_DOUBLE, rank, rank, MPI_COMM_WORLD, &mpi_stats[0]); 
          // cook
          for(i = 1; i <= smallRowSize; i++){
            for(j = 1; j <= smallColSize; j++){
              E_prev[senderY * usualRowSize + i][senderX * usualColSize + j] = recvBuffer[i][j];
            }
          }
        }
        else if (senderX == (px - 1)){
          MPI_Recv(&recvBuffer[1][1], (smallColSize+2) * usualRowSize, MPI_DOUBLE, rank, rank, MPI_COMM_WORLD, &mpi_stats[0]); 
          for(i = 1; i <= usualRowSize; i++){
            for(j = 1; j <= smallColSize; j++){
              E_prev[senderY * usualRowSize + i][senderX * usualColSize + j] = recvBuffer[i][j];
            }
          }
            
        }
        else if (senderY == (py - 1)){
          MPI_Recv(&recvBuffer[1][1], (usualColSize+2) * smallRowSize, MPI_DOUBLE, rank, rank, MPI_COMM_WORLD, &mpi_stats[0]); 
          for(i = 1; i <= smallRowSize; i++){
            for(j = 1; j <= usualColSize; j++){
              E_prev[senderY * usualRowSize + i][senderX * usualColSize + j] = recvBuffer[i][j];
            }
          }
           
        }
        else {
          MPI_Recv(&recvBuffer[1][1], (usualColSize+2) * usualRowSize, MPI_DOUBLE, rank, rank, MPI_COMM_WORLD, &mpi_stats[0]); 
          for(i = 1; i <= usualRowSize; i++){
            for(j = 1; j <= usualColSize; j++){
              E_prev[senderY * usualRowSize + i][senderX * usualColSize + j] = recvBuffer[i][j];
            }
          }

        }
      } // end of for
     free(recvBuffer); 
    }// end of if   
    else{
        MPI_Send(&myE_prev[1][1], myRowSize * (myColSize + 2), MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD);
    }
    /*
    int currentY = -1;
    for (i = 1; i <= m; i++){
      if (( (i - 1) % usualRowSize) == 0){
        currentY++;
       }
       //int myIndex = (i % myRowSize == 0) ? myRowSize : i % myRowSize; // myIndex ranges from 1 ,2 , ... , myRowSize
       int myIndex = (i % usualRowSize == 0) ? usualRowSize : i % usualRowSize; // myIndex ranges from 1 ,2 , ... , myRowSize
       if (myrank == 0){
       // receive the messages
         for (j = 1; j <= px; j++){
           if (j == 1 && currentY == 0)
             continue;
           if (j == px){
             MPI_Recv(&E_prev[i][usualColSize * (j-1) + 1], smallColSize, MPI_DOUBLE, currentY * px + j - 1, currentY * px + j - 1, MPI_COMM_WORLD, &mpi_stats[0]);
           }
            else{
             MPI_Recv(&E_prev[i][usualColSize * (j-1) + 1], usualColSize, MPI_DOUBLE, currentY * px + j - 1, currentY * px + j - 1, MPI_COMM_WORLD, &mpi_stats[0]);
            }
            //
            //
            // debug
            //int count;
            //MPI_Get_count(&mpi_stats[0], MPI_DOUBLE, &count);
            //if (i == 0)
            //cout << "myrank = " << myrank << "received " << count << " elements with mpi_tag = " << mpi_stats[0].MPI_TAG << endl; 
          }
        }
        else if (myY == currentY){
          // send a single row of data 
          for (j = 1; j <= px; j++){
            if (j == 1 && currentY == 0)
              continue;
            // TODO: fix index i
            // maybe i % usualRowSize
            if (myrank == (currentY * px + j - 1)){
              //cout << "myrank = " << myrank << " myIndex = " << myIndex << " starting col index = " << usualColSize * (j-1) + 1 << " myColSize = " << myColSize << " smallColSize = " << smallColSize << endl; 
              MPI_Send(&myE_prev[myIndex][usualColSize * (j-1) + 1], myColSize, MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD);
              //if (myIndex == m) cout << "myrank = " << myrank << " currentY = " << currentY << " px = " << px << endl;
              }
          }
        }
        MPI_Barrier(MPI_COMM_WORLD);
      } // end of for
      */
  } // end of if P != 1
  
  //cout << "Still alive !! Still alive!!  Still alive!! myrank = " << myrank << endl;
  if (myrank == 0)
  { 
   /*
    //nan check
   for (i = 1; i <= n; i++)
     for (j = 1; j <= m; j++)
        if (isnan(E_prev[i][j]))
          cout << "NANNNN i = " << i << " j = " << j  << endl;
     */
    double time_elapsed = getTime() - t0;

    double Gflops = (double)(niter * (1E-9 * n * n ) * 28.0) / time_elapsed ;
    double BW = (double)(niter * 1E-9 * (n * n * sizeof(double) * 4.0  ))/time_elapsed;

    cout << "Number of Iterations        : " << niter << endl;
    cout << "Elapsed Time (sec)          : " << time_elapsed << endl;
    cout << "Sustained Gflops Rate       : " << Gflops << endl; 
    cout << "Sustained Bandwidth (GB/sec): " << BW << endl << endl; 
   
    /* 
    // debug
    // print out elements of E_prev
    for (i = 1; i <= m ; i++){
      for (j = 1; j <= n; j++)
        cout << E_prev[i][j];
      cout << endl;
    }
    */
    double mx;
    double l2norm;
    if (P == 1){
        l2norm = stats(myE_prev,m,n,&mx);
    }
    else{
        l2norm = stats(E_prev,m,n,&mx);
    }
    //double l2norm = stats(E,m,n,&mx);
    cout << "Max: " << mx <<  " L2norm: "<< l2norm << endl;

    if (plot_freq){
      cout << "\n\nEnter any input to close the program and the plot..." << endl;
      getchar();
    }
  }
  
  free(sendBuffWest);
  free(sendBuffEast);
  free(recvBuffWest);     
  free(recvBuffEast);     

  free(myE); 
  free(myE_prev); 
  free(myR);
  

  free (E);
  free (E_prev);
  free (R);
   
  MPI_Finalize(); 
  return 0;
}
