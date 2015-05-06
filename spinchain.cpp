/***
	Finding extreme eigenvalues of quantum spin chains
	The objective is to construct this matrix and store it in sparse matrix form, 
	and then find the ground state and first excited state of the Hamiltonian using sparse matrix algorithms 
	- 	T Mathialakan, PHY 905, ECE, MSU.
***/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>	

/*
	A pair of spin
*/
float pair(bool a, bool b){ 
	//if (a&b) return -1.0; // 11>  -1;
	if (a|b) return 0.0;	// 01>  or 10> 0 // 11>  0;
	return 1.0;	//00>  1
}

/*
	Set k_th bit in a given unsigned integer n
*/
inline void set_bit(unsigned& n, unsigned k)
{
	n |= (1 << k);
}

/*
	Clear k_th bit in a given unsigned integer n
*/
inline void clear_bit(unsigned& n, unsigned k)
{
	n &= ~(1 << k);
}

/*
	Test weather k_th bit in a given unsigned integer n is set
*/
inline bool test_bit(unsigned n, unsigned k)
{
	return (n & (1 << k)) != 0;
}

/*
	Print first k bits of a given unsigned integer n from right to left
*/
inline void print_bit(unsigned n, int k)
{
	for(int i=k-1; i >= 0; i--)
		printf("\t %d", (test_bit(n, i)? 1:0));
	printf("\n");
}

/*
	Convert an integer to bits of a boolean array
*/
bool* int_to_bit(unsigned n, int k, bool* bits)
{
	//bool* bits = (bool*)malloc(sizeof(bool)*k);
	printf("k: %d \n", k);
	for(int i=k-1; i >= 0; i--){
		bits[i] = test_bit(n, i);
		printf("\t %d \n", bits[i]);}
}

/*
	
*/
void print_vect(int* a, int N){
	for(int i=0; i<N; i++)
		printf("\t %d", a[i]);
	printf("\n");
}

/*
	Print sparse matrix
*/
inline void print_matrix(int* start, int* end, int* locations, float* values, int N )
{
	 // printf("values ************ \n");
	// for(int i=0; i<16; i++)
		// printf("\t %f", values[i]);
	// printf("\n");
	int nz_st =0, nz_end=0; 
	for(int i=0; i<N; i++){
		nz_st = locations[start[i]];
		nz_end = locations[end[i]];
		//printf("\t i %d nz_st %d, nz_end %d \n",i, nz_st, nz_end );
		for(int j=0; j<nz_st; j++)	printf("\t 0");
		for(int j=nz_st; j<=nz_end; j++)	printf("\t %f",values[j] ); 
		for(int j=nz_end+1; j<N; j++)	printf("\t 0");
		printf("\n");
	}			
}

/*
	Change  an unsigned integer n by setting the bit values in a given range. It changes the state  
*/
void change_state(unsigned& n, int row, int start[], int end[], int locations[], bool values[] )
{ 
	int s = start[row];
	int e = end[row];
	for(int i=s; i<e; i++)
		if(values[i]) set_bit(n, locations[i]);
}

/*
	Calculate values for Hamiltonian matrix
*/
float hamilt(unsigned i, unsigned j, float J, float lambda, int L){
	float val_diag =0.0, val_flip =0.0;
	bool* bitsi = (bool*)malloc(sizeof(bool)*L);
	int_to_bit( i, L, bitsi);
	printf("bitsi ************ \n");
	for(int i=0; i<L; i++)
		printf("\t %d", bitsi[i]);
	printf("\n");
	// // bitsj = int_to_bit(j, N);
	for (int i=0; i<L-1; i++){
	// printf("\t pp:  %d \n", pair(bitsi[i], bitsi[i+1]));
		val_diag+=pair(bitsi[i], bitsi[i+1]);
		val_flip+=pair(!bitsi[i], bitsi[i+1])+pair(bitsi[i], !bitsi[i+1]);
	}//HL(bitsi[i])
	printf("\t val_flip: %f, val_diag: %f , J: %f\n", val_flip, val_diag, J);
	return ((val_diag));
//return 12;
}

/*
	Create Sparse matrix
*/
void build_sparse( int* start, int* end, int* locations, float* values, float J, float lambda, int L ){
	bool nz;
	int k =0;
	float val =0.0;
	int N = pow(2, L);
	for(unsigned i=0; i<N; i++){
		nz = false;
		for(unsigned j=0; j<N; j++){
			val = hamilt(i,j, J, lambda, L);
			if(nz){
				values[k] = val;
				locations[k] = j;
				end[i] = k; //printf("nz k: %d, i: %d, j : %d end: %d  locations[%d]: %d\n", k, i, j, end[i], k, locations[k]);	
				k++; 				
			}else{
				if(val!=0){
					nz = true;
					values[k] = val;
					locations[k] = j;
					start[i] = k;
					end[i] = k; //printf("k: %d, i: %d, j: %d start: %d end: %d , locations[%d]: %d\n", k, i, j, start[i], end[i], k,locations[k]);	
					k++; 				
				}
			}
		}	
	}
	// printf("values ************ \n");
	// for(int i=0; i<16; i++)
		// printf("\t %f", values[i]);
	// printf("\n");
	// printf("location \n");
	// print_vect(locations, (N*N));
	// printf("end \n"); print_vect(end, N);
	// printf("st \n"); print_vect(start, N);
}

/*
	Lanczos diagonalisation  
*/
void diag( float A, float* beta, float* alpha, float* w, float* v, int N){
	int i=0;
	for( i=0; i<N; i++){
		w[i] = A*v[i];
		alpha[i] = w[i]*v[i];
		w[i] -= (alpha[i]*v[i] - beta[i]*v[i-1]);
		beta[i+1] = abs(w[i]);
		v[i+1] = w[i]/beta[i+1];
 	}
	w[i] = A*v[i];
	alpha[N] = w[N]*v[N];
	
}

/*
	Matrix vector product
*/
void matrix_product(int* start, int* end, int* locations, float* values, float* vector, float* r_vect, int N){
	int nz_st =0, nz_end=0; 
	for(int i=0; i<N; i++){
		nz_st = locations[start[i]];
		nz_end = locations[end[i]];
		r_vect[i] =0.0;
		for(int j=0; j<N; j++)
			if ((nz_st<=j)&&(j<=nz_end))
				r_vect[i] += values[j]*vector[j];
	}		
}


/*
	Total magnetization in a Lattice
*/
int magnetization(bool L[], int n){
	int sum=0;
	 for(int i=0; i<n-1; i++)
		 sum=(int)L[i];
	return 2*sum - n;
}

/*
	The energy by a single spin formed with neighbours
*/
int spotvalue(bool L[], int k, bool spin, int c, int n){
	int val=0;
	if(k%c > 0) val += pair(L[k-1], spin ); // Left neighbour
	if(k%c < c-1) val += pair(L[k+1], spin ); // Right neighbour
	if(k>c) val += pair(L[k-c], spin ); // Above neighbour
	if(k<n-c) val += pair(L[k+c], spin); // Below neighbour
	return val;
}

/*
	The energy by a single spin formed with neighbours in periodic boundary condition
*/
int spotvalue_PBC(bool L[], int k, bool spin, int c, int n){
	int val=0, r=c;
	if(k%c > 0) val += pair(L[k-1], spin );  else val += pair(L[k+c-1], spin ); // Left neighbour
	if(k%c < c-1) val += pair(L[k+1], spin ); else val += pair(L[k-c+1], spin ); // Right neighbour
	if(k>c) val += pair(L[k-c], spin ); else val += pair(L[k+c*(r-1)], spin ); // Above neighbour
	if(k<n-c) val += pair(L[k+c], spin); else val += pair(L[k-c*(r-1)], spin ); // Below neighbour
	return val;
}



int main(int argc, char** arg){
	printf("spin chain  \n");
	// unsigned bits = 9;
	// unsigned un =2;
	// print_bit(bits, 5);
	// printf( "is set   %d , bits %d \n", test_bit(bits, un), bits);
	// set_bit(bits, un);
	// printf( "is set   %d , bits %d \n", test_bit(bits, un), bits);
	// clear_bit(bits, un);
	// printf( "is set   %d , bits %d \n", test_bit(bits, un), bits);
	// print_bit(bits, 5);
	float J = 1.0;
	float lambda = 1.0;
	int L = 6;
	int N = pow(2, L);
	int size = N*N ; // number of elements in sparse matrix 
	float *values = (float *) malloc(sizeof(float)*size); 
	int *locations = (int *) malloc(sizeof(int)*size);
	int *start = (int *) malloc(sizeof(int)*N);
	int *number = (int *) malloc(sizeof(int)*N);	
	//int *state = (int *) malloc(sizeof(int)*N);
	build_sparse( start, number, locations, values, J,lambda, L );
	print_matrix( start, number, locations, values, N );
	return 0;
}
