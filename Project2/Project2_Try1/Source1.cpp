#include <iostream>
#include <cmath>
#include <armadillo>
#include "time.h"


using namespace std;
using namespace arma;


void rotation(mat &A, mat &R, int k, int l, int n) {
	// Take in our matrix to be rotated A, our rotation matrix R, the index values of the largest off diagonal elements k and l, and the x and y dimensions n
	// Rotate A and update R to reflect it, changing their values at their location


	//Declare the values needed for rotation, t, tau, c (cosine), and s (sine)
	double t, tau;
	double c, s;


	//Store the 3 values A(k,k), A(l,l), and A(k,l) so the value isn't lost during rotation
	double a_kk = A(k, k);
	double a_ll = A(l, l);
	double a_kl = A(k, l);

	//Declare the iterated values now. Used in the rotation process
	double a_ik, a_il;
	double r_ik, r_il;
	
	//Define tau
	tau = (A(k, k) - A(l, l)) / (2.0 * A(k, l));

	//If tau is positive t is as specified below
	if (tau > 0.0) {
		t = 1.0 / (tau + sqrt(1.0 + tau*tau));

	}
	
	//If tau is negative (or zero) then t is as specified below
	else {
		t=-1.0 / (-tau + sqrt(1.0 + tau*tau));
	};

	//Define c and s
	c = 1.0 / sqrt(1.0 + t*t);
	s = t*c;

	//Set the values for our new matrix elements A(k,k), A(l,l), A(k,l), and A(l,k)
	A(k, k) = a_kk*c*c - 2.0 * a_kl*c*s + a_ll*s*s;
	A(l, l) = a_ll*s*s - 2.0 * a_kl*c*s + a_kk*c*c;
	A(k, l) = 0;
	A(l, k) = 0;

	//Iterate through all matrix elements in row k, column k, row l, and column l	
	for (int i = 0; i < n; i++) {
		//for all values of i smaller than n

		//Skip values of i that are equal to k or l (we've already defined these)
		if (i != k && i != l) {

			//Store the values of A(i,k) and A(i,l) so we can use them later without them changing.
			a_ik = A(i, k);
			a_il = A(i, l);

			//Set our new values of A(i,k), A(i,l), A(l,i), and (k,i)
			A(i, k) = a_ik*c - a_il*s;
			A(i, l) = a_il*c - a_ik*s;
			A(k, i) = A(i, k);
			A(l, i) = A(i, l);

		};

		//Store our values R(i,k) and R(i,l) so we can use them without them changing
		r_ik = R(i, k);
		r_il = R(i, l);

		//Set our new values of R(i,k) and R(i,l)
		R(i, k) = c*r_ik - s*r_il;
		R(i, l) = c*r_il + s*r_ik;
	};
	return;
}

double max_off_diag(mat A, int &k, int &l,int n) {
	//Take in matrix A, the index values k and l, and the size of our matrix n
	//Return the index values k and l of our max off diag value

	//Initialize our place holder index values i and j, the arma index value ij, and our maximum value max
	int i, j;
	uword ij;
	double max,min;

	//set the # of max iterations to n*n so we don't end up looking past the end of the matrix
	int max_iterations = n+1;

	//Using arma::mat built in functions, find the index value of our maximum off diagonal value
	for (int iteration = 0; iteration < max_iterations; iteration++) {
		max = A.max();
		min = A.min();
		

		if (max == 0.0 && min==0.0) {
			// If the maximum value of the matrix is 0 (i.e.It's a zero matrix) returns the first possible index and max value 0.
			// This is a worst case scenario. Ideally there will be an absolute off diagonal value
			k = 0;
			l = 1;
			max = 1.0e-9;
			return max;
		}

		else if (fabs(max) > fabs(min)) {

			//Using built in functions, find the current max value of the matrix and store its index values in i and j
			ij = A.index_max();
			
		}
		else {
			ij = A.index_min();
		};

		i = ij % n;
		j = ij / n;
		//If the value is on the diagonal set it to zero and start over.
		if (i == j) {
			A(i, j) =0;
		}

		//If the value is not on the diagonal set k and l to i and j, max to A(i,j) and return max
		else {
			k = i;
			l = j;
			max = A(i, j);
			return max;
		};

	};
	//In case of catostrophic failure return 0
	return 0;
}

void jacobi(mat &A, mat &R, int n) {

	// A is our starting matrix, R is our rotation matrix and is empty, n is the x and y dimensions of the matrix A.

	// Initialize R as an n by n identity matrix

	R.eye(n, n);

	//Initialize some constants and variables we'll need

	int k, l;

	double epsilon = 1.0e-8;

	//Our maximum number of iterations is n^3
	double max_num_iterations = (double)n*(double)n* (double)n;

	//We're currently on iteration 0
	int iterations = 0;


	//find the current biggest off diagonal value of matrix A

	double max = max_off_diag(A, k, l, n);

	//While the absolute value of the max off diagonal value of A is bigger than epsiolon (almost 0) and we have not yet hit n^3 number of iterations
	while (fabs(max) > epsilon && (double)iterations < max_num_iterations) {


		//Rotate the matrix using the previously found values of k and l
		rotation(A, R, k, l, n);

		//Find the new max off diagonal value
		max = max_off_diag(A, k, l, n);

		//Increment the number of iterations
		iterations++;
	};
}

void unit_test_max_off_diag() {
	//test max_off_diag. If it works k will be 0 and l will be 1

	int n = 5;
	int k, l;

	mat A(n, n, fill::eye);

	A(0, 0), A(1, 1), A(2, 2), A(3, 3), A(4, 4) = 12.0;
	A(0, 1) = 5.0;

	cout << "Matrix A (max value test):" << endl;
	A.print();

	max_off_diag(A, k, l, n);

	cout << "k:" << k << ", l:" << l << endl;

	cin >> n;

}


void unit_test_jacobi() {
	//Tests jacobi. If the other two functions worked then jacobi should work.
	
	//Set the desired size

	int n = 5;

	//Create 2 n*n matricies, A and R.

	mat A(n, n, fill::zeros);
	mat R;
	

	//Set some random test values of A

	A(0, 0) = 16.0;
	A(0, 1) = 2.0;
	A(1, 2) = 5.0;
	A(1, 1) = 25.0;
	A(1, 3) = 4.0;
	A(2, 2) = 14.0;
	A(2, 3) = 1.0;
	A(3, 4) = 3.0;
	A(2, 4) = 2.0;

	cout << "Matrix A (Jacobi Test):" << endl;
	A.print();

	
	//Run our jacobi function
	jacobi(A,R,n);


	//print our final values of A and R

	cout << "Matrix A (Final):" << endl;
	A.print();

	cout << endl;

	cout << "Matrix R:" << endl;
	R.print();

	cout << endl;


	//Wait for input
	cin >> n;

}

mat matrix_generator_n_i(mat &A, int n, double h, double p = 0) {

	A.zeros(n, n);

	double e = -1 / (h*h);
	
	for (int i = 0; i < n; i++) {

		A(i, i) = 2 / (h*h) + (p + (double)(i + 1)*h)*(p + (double)(i + 1)*h);

		if (i != (n - 1)) {

			A(i, i + 1) = e;
			A(i + 1, i) = e;
		}
	};

	return A;
}

mat matrix_generator_i(mat &A, int n, double h, double w, double p = 0) {

	A.zeros(n, n);

	double e = -1 / (h*h);

	for (int i = 0; i < n; i++) {

		A(i, i) = 2 / (h*h) + (w*(p + (double)(i + 1)*h))*(w*(p + (double)(i + 1)*h))+1/(p + (double)(i + 1)*h);

		if (i != (n - 1)) {

			A(i, i + 1) = e;
			A(i + 1, i) = e;
		}
	};

	return A;
}

void non_interacting(int n, double h) {
	
	//int wait;

	mat A;
	mat R;

	matrix_generator_n_i(A, n, h);

	jacobi(A, R, n);

	cout << "For the non-interacting case" << endl;
	cout << "Matrix A=" << endl;
	A.print();
	cout << "" << endl;
	cout << "Matrix R=" << endl;
	R.print();
	cout << "" << endl;
	//cin >> wait;

}

void interacting_case(int n, double h) {

	clock_t start, finish; // declare start and final time
	//int wait;

	//Create 2 n*n matricies, A and R.

	mat A;
	mat R;

	double omega[4] = { 0.01, 0.5, 1, 5 };

	for (int i = 0; i < 4; i++) {
		
		matrix_generator_i(A, n, h, omega[i]);
		start = clock();
		jacobi(A, R, n);
		finish = clock();

		cout << "For the interacting case where w=" << omega[i] << endl;
		cout << "" << endl;
		cout << "Matrix A=" << endl;
		A.print();
		cout << "" << endl;
		cout << "Matrix R=" << endl;
		R.print();
		cout <<""<< endl;
		cout << "Runtime: " << (finish - start) << endl;

		//cin >> wait;

		A.clear();
		R.clear();

	}
	
}

int main() {
	//Run 2 unit tests to see if my functions work

	//unit_test_max_off_diag();

	//unit_test_jacobi();

	//mat n_i;

	//mat i;

	//matrix_generator_n_i(n_i, 5, 5);

	//cout << "Non-Interacting Matrix:" << endl;
	//n_i.print();
	//cout << "" << endl;

	//matrix_generator_i(i, 5, .5, 5);

	

	//cout << "Interacting Matrix:" << endl;
	//i.print();
	//cout << "" << endl;
	int x;
	//cin >> x;

	non_interacting(5,5);

	interacting_case(5,5);

	cin >> x;
}
