#include "matrix.h"


extern "C" {

    void los_cp (double* times,        // transition times
		 double* ma,           // transition matrices
		 int* len,             // number of transitions
		 int* rows,            // row number of transition matrice
		 int* cols,            // colum nnumber of transition matrice
		 double* los1,         // LOS given state 1
		 double* los0,         // LOS given state 0
		 double* phi2case,     
		 double* phi2control,
		 double* phi3case,
		 double* phi3control,
		 double* tau         
	) {
        
	Vector Times(times,*len);
    
	Array Ma(ma,*rows,*cols,*len );

	Vector Los1(los1, *len);
	Los1[*len] = *tau;
 
	Vector Los0(los0,*len);
	Los0[*len] = *tau;

	Vector Phi2case(phi2case,*len);

	Vector Phi2control(phi2control,*len);

	Vector Phi3case(phi3case,*len);

	Vector Phi3control(phi3control,*len);


	Matrix Diag(*rows, *cols);    
	Diag.identity();

	Array A;
	A.push_back(Diag);

	Array A2;
	A2.push_back(Diag);

	Vector T;
	T.push_back(*tau);
    
	Vector T2;

	for(int i = (Times.size() - 2); i >= 0; --i) {
	    itVector vpos = T.begin();
	    T.insert(vpos, Times[i+1]);

	    itVector vpos2 = T2.begin();
	    T2.insert(vpos2, Times[i+1]);

	    Vector Diff = T.diff();

	    A = Ma[i+1]*A;
    
	    Vector a11;
	    Vector a00;
	    Vector a01;

	    for(int j = 0; j < A.size(); ++j) {
		a11.push_back( A[j][1][1] );
		a00.push_back( A[j][0][0] );
		a01.push_back( A[j][0][1] );
	    }

	    Los1[i] = Times[i+1] + scalar(Diff, a11);
      
	    Los0[i] = Times[i+1] + scalar(Diff, (a00 + a01));

	    if( i == (Times.size() - 2)) {
		Phi2case[i] = Times[(Times.size()-1)] * A[(A.size()-1)][1][2];
		
		Phi3case[i] = Times[(Times.size()-1)] * A[(A.size()-1)][1][3];	
	    }
	    else {
		Vector Diff2 = T2.diff();
		
		//cout << Diff2 << endl;
		    
		A2 = Ma[i+1]*A2;
		
		Vector a12;
		Vector a13;
		
		for(int l = 0; l < A2.size(); ++l) {
		    a12.push_back( A2[l][1][2] );
		    a13.push_back( A2[l][1][3] );
		}	
		
		Phi2case[i] = (Times[(Times.size()-1)] * A[(A.size()-1)][1][2]) - scalar(Diff2, a12);
		
		Phi3case[i] = (Times[(Times.size()-1)] * A[(A.size()-1)][1][3]) - scalar(Diff2, a13);
		
		
		// stack identity matrix on top for the next loop	           
		itArray apos2 = A2.begin();
		A2.insert(apos2, Diag);
	    }
	    
	    Phi2control[i] = A[(A.size()-1)][1][2] * Los0[i];
	    
	    Phi3control[i] = A[(A.size()-1)][1][3] * Los0[i];
		
	    // stack identity matrix on top for the next loop	           
	    itArray apos = A.begin();
	    A.insert(apos, Diag);
	}
	
	Los1.as_double(los1);

	Los0.as_double(los0);

	Phi2case.as_double(phi2case);

	Phi2control.as_double(phi2control);

	Phi3case.as_double(phi3case);

	Phi3control.as_double(phi3control);
    }
}


/* To be called when there's no competing risks */
extern "C" {

    void los_nocp(double *times,
		  double *ma,
		  int *len,
		  int *rows,
		  int *cols,
		  double *los1,
		  double *los0,
		  double *tau) {


	Vector Times(times,*len);
    
	Array Ma(ma,*rows,*cols,*len );

	Vector Los1(los1, *len);
	Los1[*len] = *tau;
 
	Vector Los0(los0,*len);
	Los0[*len] = *tau;

	Matrix Diag(*rows, *cols);    
	Diag.identity();

	Array A;
	A.push_back(Diag);

	Vector T;
	T.push_back(*tau);
    
	for(int i = (Times.size() - 2); i >= 0; --i) {
	    itVector vpos = T.begin();
	    T.insert(vpos, Times[i+1]);
	 
	    Vector Diff = T.diff();

	    A = Ma[i+1]*A;
    
	    Vector a11;
	    Vector a00;
	    Vector a01;

	    for(int j = 0; j < A.size(); ++j) {
		a11.push_back( A[j][1][1] );
		a00.push_back( A[j][0][0] );
		a01.push_back( A[j][0][1] );
	    }

	    Los1[i] = Times[i+1] + scalar(Diff, a11);
      
	    Los0[i] = Times[i+1] + scalar(Diff, (a00 + a01));

	    // stack identity matrix on top for the next loop	           
	    itArray apos = A.begin();
	    A.insert(apos, Diag);
	}
	
	Los1.as_double(los1);

	Los0.as_double(los0);
    }
}
