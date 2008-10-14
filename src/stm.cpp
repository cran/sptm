#include <cstdlib>
#include <iostream>
#include <math.h>

using namespace std;

extern "C" {  
       
       
void V1 (int * _n, int * _p, double * _u, double * _V) {
	
	int n = *_n;
	int p = *_p;
//		double (*u)[n][n] = (double (*)[n][n]) _u;
//		double (*V)[p] = (double (*)[p]) _V;
   
	
	for (int q1=0; q1<p; q1++)
		for (int q2=0; q2<p; q2++)
//				V[q1][q2] = 0;
			_V[q1*p+q2] = 0;

//		// print u.ij for debugging
//		for (int i=0; i<2; i++)
//			for (int j=0; j<2; j++) {
//				for (int q1=0; q1<p; q1++) {
//					cout << u[i][j][q1] << " ";
//				}
//				cout << endl;
//			}
	
	for (int i=0; i<n; i++)
		for (int j=0; j<n; j++)
			// save (u[i][j][q1] + u[j][i][q1]) in a var here does not make a difference
			for (int k=0; k<n; k++) {
				// enforce k!=j
				if(k==j) continue;
				// strange but true: at least on orca1 that moving loop q1,q2 outside loop i,j,k takes much longer
				for (int q1=0; q1<p; q1++)
					// V is symmetric matrix, so only need to compute half of the matrix
					// let q2 go from 1 to p almost triples the time
					for (int q2=q1; q2<p; q2++) 
//							V[q1][q2] += (u[q1][i][j] + u[q1][j][i]) * (u[q2][i][k] + u[q2][k][i]);
						_V[q1*p+q2] += (_u[q1*n*n+i*n+j] + _u[q1*n*n+j*n+i]) * (_u[q2*n*n+i*n+k] + _u[q2*n*n+k*n+i]);
			}
						
	for (int q1=0; q1<p; q1++) {
		for (int q2=q1; q2<p; q2++) {
//				V[q2][q1] = V[q1][q2];
			_V[q2*p+q1] = _V[q1*p+q2];
		}
	}		
}


void V1_ph1 (int * _n, int * _p, double * _e, double * _V, double * _sampling_p) {
	
	int n = *_n;
	int p = *_p;
//		double (*e)[n][n] = (double (*)[n][n]) _e;
//		double (*V)[p] = (double (*)[p]) _V;
	
	for (int q1=0; q1<p; q1++)
		for (int q2=0; q2<p; q2++)
			_V[q1*p+q2] = 0;

	for (int i=0; i<n; i++)
		for (int j=0; j<n; j++)
			for (int k=0; k<n; k++) {
				// enforce k!=j
				if(k==j) continue;
				// strange but true: at least on orca1 that moving loop q1,q2 outside loop i,j,k takes much longer
				for (int q1=0; q1<p; q1++)
					// V is symmetric matrix, so only need to compute half of the matrix
					// let q2 go from 1 to p almost triples the time
					for (int q2=q1; q2<p; q2++) 
						_V[q1*p+q2] += (_e[q1*n*n+i*n+j] + _e[q1*n*n+j*n+i]) * (_e[q2*n*n+i*n+k] + _e[q2*n*n+k*n+i]) / (_sampling_p[i]*_sampling_p[j]*_sampling_p[k]);
			}
						
	for (int q1=0; q1<p; q1++) {
		for (int q2=q1; q2<p; q2++) {
			_V[q2*p+q1] = _V[q1*p+q2];
		}
	}		
}

void get_phi (int * _n, int * _p, double * _u, double * _phi) {
	
	int n = *_n;
	int p = *_p;
//		double (*u)[n][n] = (double (*)[n][n]) _u;
//		double (*phi)[n] = (double (*)[n]) _phi;
	
	for (int q=0; q<p; q++)
		for (int i=0; i<n; i++) {
			_phi[q*n+i] = 0;
                for (int j=0; j<n; j++) {
			        _phi[q*n+i] += (_u[q*n*n+i*n+j] + _u[q*n*n+j*n+i]);
                }                        
        }
}


void V2 (int * _n, int * _p, double * Z, double * V, 
     double * Xcc, int * N, int * d, double * X ) {
	
	int n = *_n;
	int p = *_p;
	int capN= *N;
	
	for (int q1=0; q1<p; q1++)
		for (int q2=0; q2<p; q2++)
			V[q1*p+q2] = 0;

	// int N1=1; // for debugging
	//double q_t[p];
    // cannot use new b/c it is C not C++
    double * q_t = (double *)malloc(p * sizeof(double));
	for (int i=0; i<capN; i++) {

		double _t = X[i];
		
		for (int q=0; q<p; q++) {
			q_t[q]=0;
			for (int j=0; j<n*n; j++) {
				q_t[q] += Z[j*p+q] * (Xcc[j]>=_t?1:0);
			}
		}			

		int c=0;
		for (int i1=0; i1<capN; i1++)
			c += X[i1]>=_t;
			
		for (int q1=0; q1<p; q1++) {
			for (int q2=q1; q2<p; q2++) {
				V[q1*p+q2] += (1-_t)/(c*c) * q_t[q1]*q_t[q2];
				// for debugging
				//cout << V[q1][q2] << endl;
			}
        }
        		
	}
	free(q_t);
		
	for (int q1=0; q1<p; q1++) {
		for (int q2=q1; q2<p; q2++) {
			V[q2*p+q1] = V[q1*p+q2];
		}
	}		
	
//		delete [] q_t;

}

}

// commented out b/c cran does not like e[n][n][p] or V[p][p]
//int main(int argc, char *argv[])
//{
//	int n=3, p=1;
//	double e[n][n][p];
//	double V[p][p];
//	
//	for (int i=0; i<n; i++)
//		for (int j=0; j<n; j++) 
//			for (int q1=0; q1<p; q1++) {
//				e[i][j][q1]	= (pow(i+1,2)-(j+1))+q1+1;
//				//cout << e[i][j][q1] << endl;
//			}
//		
//	V1(&n, &p, &e[0][0][0], &V[0][0]);
//
////	for (int q1=0; q1<p; q1++) {
////		for (int q2=0; q2<p; q2++)
////			cout << V[q1][q2] << " ";
////		cout << endl;
////	}
//
//    return 1;
//}

