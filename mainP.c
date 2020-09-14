#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "mpi.h"

#define DSIZE 2000
#define h     1
#define PI 3.14159265359


int main(int argc, char* argv[])
{
	//- general variables used in program 
	int i,j; 
	unsigned int SIZE = (int)DSIZE*(int)DSIZE ; 
	int myNumberofNonZero = 0,mySIZE=0,NumberofNonZero;
	int BandWidth,IrowIndex,Nine_BandWidth,myNine_BandWidth; 
	double TotalElement, Ratio_NonZero, Ratio_BandWidth; 
	int nProc,myID;

	//- define nine-banded dominant matrix constants 
	double A, B, C;
	A = (-1)/(12*(double)h*(double)h); 
	B = 16/(12*(double)h*(double)h); 
	C = (-5)/((double)h*(double)h); 	
	//- initilize MPI 
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nProc);	
	MPI_Comm_rank(MPI_COMM_WORLD,&myID); 
	//- calculate rowLoad per core and residual if there is  
	int rowLoad = (int)SIZE/nProc;
	int residual = (int)SIZE % nProc;
	//- calculate myNumberofNonZero and mySIZE as distributed
	for(i=myID*rowLoad; i<myID*rowLoad+rowLoad; i++)
	{
		myNumberofNonZero++;
		mySIZE++;

		if((i+1)<SIZE)
			myNumberofNonZero++;
		
		if((i+2)<SIZE)
			myNumberofNonZero++;
		
		if((i+(int)DSIZE)<SIZE)
			myNumberofNonZero++; 
		
		if((i+2*(int)DSIZE)<SIZE)
			myNumberofNonZero++;
	}
	if(myID==nProc-1 && residual!=0)
	{
		for(i=myID*rowLoad+rowLoad; i<myID*rowLoad+rowLoad+residual; i++)
		{
			myNumberofNonZero++;
			mySIZE++; 			

			if((i+1)<SIZE)
			myNumberofNonZero++;
		
			if((i+2)<SIZE)
			myNumberofNonZero++;
		
			if((i+(int)DSIZE)<SIZE)
			myNumberofNonZero++; 
		
			if((i+2*(int)DSIZE)<SIZE)
			myNumberofNonZero++;
		}
	}
	//- sum all myNumberofNonZero values and reduce value to all processes
	MPI_Reduce(&myNumberofNonZero,&NumberofNonZero,1,MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
	
	if(myID==0)
	{
		//- calculate some statistics and print those on screen 		
		TotalElement = (double)SIZE * (double)SIZE;
		Ratio_NonZero = (double)NumberofNonZero / TotalElement;
		printf("Matrix Size %d*%d\n", SIZE, SIZE);
		printf("Degree of Parallelization %d\n", nProc); 
		printf("Number of NonZero stored = %d\n", NumberofNonZero);
		//- because matrix is symmetric, dublicate Ratio_NonZero while print it on the screen 
		printf("NonZero Dominance = %.16f\n", 2*Ratio_NonZero);
		printf("rowLoad per processor %d\n", rowLoad); 
		printf("residual %d rows and it is added to %d. processor\n", residual, nProc-1); 
	}
	//- allocate vectors according to their sizes 
	double* vector=(double*)malloc(SIZE*sizeof(double)); 
	double* myValues=(double*)malloc(myNumberofNonZero*sizeof(double)); 	
	int* myColumns=(int*)malloc(myNumberofNonZero*sizeof(int)); 
	int rowIndexSize=mySIZE+1;
	int* rowIndex=(int*)malloc(rowIndexSize*sizeof(int));
	double* RHS=(double*)calloc(mySIZE,sizeof(double));
	//- assign values to vector
	double deltaH_=1.0/(double)SIZE;  	
	double H_=0;
	for(i=0; i<SIZE; i++)
	{
		vector[i]=cos(H_*PI);
		H_+=deltaH_; 
	} 
 	//- assign values to sparse Matrice's vectors and calculate some values to use in some statistic calculations  
	int cnt=0; 
	IrowIndex=-1;
	rowIndex[0]=0;
	Nine_BandWidth=0;  
	myNine_BandWidth=0; 
	for(i=0; i<mySIZE; i++)
	{
		myValues[cnt]  = C;
		myColumns[cnt] = i+rowLoad*myID;  
		cnt++;  
		IrowIndex++; 
		BandWidth = 1; 
		
		if(((i+rowLoad*myID)-2*(int)DSIZE)>=0)
		{
			BandWidth++;
 		}
		
		if(((i+rowLoad*myID)-(int)DSIZE)>=0)
		{
			BandWidth++; 
		}
		
		if(((i+rowLoad*myID)-2)>=0)
		{
			BandWidth++;
		}
 		
		if(((i+rowLoad*myID)-1)>=0)
		{
			BandWidth++;
		}
		
		if(((i+rowLoad*myID)+1)<SIZE)
		{
			myValues[cnt]  = B     ;
			myColumns[cnt] = (i+rowLoad*myID)+(1) ;  
			cnt++; 
			IrowIndex++; 
			BandWidth++; 
		}
		
		if(((i+rowLoad*myID)+2)<SIZE)
		{
			myValues[cnt]  = A      ; 
			myColumns[cnt] = (i+rowLoad*myID)+(2)  ;
			cnt++; 
			IrowIndex++; 
			BandWidth++;
		}

		if(((i+rowLoad*myID)+(int)DSIZE)<SIZE)
		{
			myValues[cnt]  = B               ; 
			myColumns[cnt] = (i+rowLoad*myID)+((int)DSIZE)  ;
			cnt++; 
			IrowIndex++; 
			BandWidth++; 
		}
		
		if(((i+rowLoad*myID)+2*(int)DSIZE)<SIZE)
		{
			myValues[cnt]  = A                 ;
			myColumns[cnt] = (i+rowLoad*myID)+(2*(int)DSIZE)  ; 
			cnt++;
			IrowIndex++;  
			BandWidth++;
		}
		
		if( BandWidth == 9 )
				myNine_BandWidth++;
	
		if( i != (mySIZE-1) )
			rowIndex[i+1] = IrowIndex + 1 ;
			
	}
	rowIndex[mySIZE] = myNumberofNonZero;
	//- reduce Nine_BandWidth value to zero ID processor to print it and its statistic 
	MPI_Reduce(&myNine_BandWidth,&Nine_BandWidth, 1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD); 
	if(myID==0)
	{
		Ratio_BandWidth = (double)Nine_BandWidth / (double)SIZE ; 	
		printf("Number of Nine BandWidth Row = %d\n", Nine_BandWidth );	
		printf("Nine BandWidth Dominance = %.16f\n", Ratio_BandWidth );	
	}
	//- start measuring wall clock time and matrix-vector multiplication  
	double time1_=MPI_Wtime(); 
	int   bandI_,diagColI_,colI_; //,realRow_;
	double sum;
	cnt=0;
	for(i=0; i<mySIZE; i++)
	{
		bandI_=rowIndex[i+1]-rowIndex[i];
		sum=0.0;
		//realRow_=myID*rowLoad+i;
		//- for stored part
		for(j=0; j<bandI_;j++)
		{
			if(j==0)
			{
				//- for diagonal part
				diagColI_=myColumns[cnt];
				sum+=myValues[cnt]*vector[diagColI_];
				cnt++; 
			}
			else
			{
				//- if element is not on diagonal of matrix  
				colI_=myColumns[cnt];
				sum+=myValues[cnt]*vector[colI_]; 
				
				//- for sym part  
				//RHS[colI_]+=myValues[cnt]*vector[realRow_];
				cnt++;
			}
		}
		//- for unstored part 
		if(((i+rowLoad*myID)-2*(int)DSIZE)>=0)
		{
			colI_=((i+rowLoad*myID)-2*(int)DSIZE);
			sum+=A*vector[colI_];			
 		}
		if(((i+rowLoad*myID)-(int)DSIZE)>=0)
		{
		 	colI_=((i+rowLoad*myID)-(int)DSIZE);
			sum+=B*vector[colI_]; 
		}
		if(((i+rowLoad*myID)-2)>=0)
		{
			colI_=((i+rowLoad*myID)-2);
			sum+=A*vector[colI_];
		}			
		if(((i+rowLoad*myID)-1)>=0)
		{
			colI_=((i+rowLoad*myID)-1);
			sum+=B*vector[colI_]; 
		}
		//RHS[realRow_]+=sum; 
		RHS[i]+=sum; 
	}
	//- wait all processors finish their jobs 
	MPI_Barrier(MPI_COMM_WORLD); 
	//- finish measuring and print Wall Clock Time 
	if(myID==0)
	{ 
		double time2_=MPI_Wtime(); 
		printf("Wall Clock Time = %.16f seconds\n", time2_-time1_); 
	}
	//if(myID==1)
		//for(i=0;i<mySIZE;i++)
			//printf("RHS[%d] = %.16f\n", i,RHS[i]);
	//- deallocate vectors 
	free(vector); 
	free(myValues);
	free(myColumns); 
	free(rowIndex); 
	free(RHS);
	//- finalize MPI 
	MPI_Finalize();
	return 0;
}


	

	 

	
