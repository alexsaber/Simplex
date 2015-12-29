// LPSolve(OPS).cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include "Funktion.h"



using namespace std;


double* lpsolve(int n, double *c, int k, double **A, double *b){
	//Schlupfvariablen, befuellen mit '1'
	int *s_vars = new int[k];
	for (int i = 0; i < k; i++) {
		s_vars[i] = 1;
	}

	//create table





	return nullptr;
}



int _tmain(int argc, _TCHAR* argv[]){

	ifstream file("testfile.txt");
	int n, k; // n Zeilen, k Spalten
	file >> n >> k;
	cout << "n: " << n << ", k: " << k << endl;

	double *vector_c = new double[n];
	for (int i = 0; i < n; i++) {
		file >> vector_c[i];
	}
	//testing vector_c
	cout << "vector_c" << endl;
	for (int i = 0; i < n; i++) {
		cout << vector_c[i] << endl;
	}

	double** matrix = new double*[n];
	double *vector_b = new double[k];
	int vector_b_counter = 0;
	for (int i = 0; i < k; i++) {
		matrix[i] = new double[n];
		for (int j = 0; j < n; j++) {
			file >> matrix[i][j];	
		}
		file >> vector_b[vector_b_counter++];
		
	}

	//testing vector_c
	cout << "vector_b" << endl;
	for (int i = 0; i < k; i++) {
		cout << vector_b[i] << endl;
	}


	//results
	double* results = lpsolve(n, vector_c, k, matrix, vector_b);
	for (int i = 0; i < n; i++) {
//		cout << vector_b[i] << endl;
	}


	system("pause");


	return 0;
}

