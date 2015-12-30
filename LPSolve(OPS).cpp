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

	//create table and fill up with the 
	double** table = new double*[k + 1];//Spalten

	int schlupf_vars_counter_i = 0;//fuer Schluepfvariable
	int schlupf_vars_counter_j = n;//fuer Schluepfvariable

	for (int i = 0; i < k; i++){
		table[i] = new double[k + n + 2];//Zeilen
		//kopiere die Werte aus der ausgelesenen Matrix
		for (int j = 0; j < n; j++){
			table[i][j] = A[i][j];
		}
		for (int j = n; j < n + k; j++){
			if (i == schlupf_vars_counter_i && j == schlupf_vars_counter_j){
				table[i][j] = 1;
				schlupf_vars_counter_i++;
				schlupf_vars_counter_j++;
			}	
			else
				table[i][j] = 0;
		}
		
		table[i][n + k] = b[i];
	}

	table[k] = new double[k + n + 1];//Letzte Zeile

	//multiplying Zielfukntionargumente with -1
	for (int i = 0; i < n; i++){
		table[k][i] = c[i] * -1;
	}
	//fillng up Schluepfvariablen for last row with zeros
	for (int i = n; i < n+k+1; i++){
		table[k][i] = 0;
	}

	//filling up Quotient column with zeros
	for (int i = 0; i < k + 1; i++){
		table[i][k + n + 1] = 0;
	}

	//testing tabelle
	cout << "testing table" << endl;
	for (int i = 0; i < k + 1; i++) {
		for (int j = 0; j < k + n + 2; j++) {
			cout << table[i][j] << " ";
		}
		cout << endl;
	}


	//----BEGIN CALCULATION-----
	

	int number_of_rows = k + 1;
	int number_of_columns = k + n + 2;
	int index_quotient_column = number_of_columns - 1;

	cout << "number_of_rows: " << number_of_rows << endl;
	cout << "number_of_columns: " << number_of_columns << endl;

	//find the highest number (without "-") in the last row OR find pivotcolumn

	int index_pivot_column = 0;
	double temp_highest = 0;

	for (int i = 0; i < index_quotient_column; i++) {
		if (i == 0){
			temp_highest = table[number_of_rows - 1][i];
			index_pivot_column = i;
		}
		else 
			if (table[number_of_rows - 1][i] < 0 && table[number_of_rows - 1][i] < temp_highest){
				temp_highest = table[number_of_rows - 1][i];
				index_pivot_column = i;
			}
	}

	cout << "The highest number (without \" - \") is :" << temp_highest << " at index " << index_pivot_column << endl;

	//calculate Quotient
	
		for (int i = 0; i < number_of_rows - 1; i++) {
			if (table[i][index_pivot_column] != 0){
				table[i][index_quotient_column] = table[i][number_of_columns - 2] / table[i][index_pivot_column];
			}
		}
		cout << "Quotient: " << endl;

		for (int i = 0; i < number_of_rows - 1; i++) {

			cout << table[i][index_quotient_column] << " ";
			
		}
		cout << endl;



		//find the smallest number in the Quotient column OR find pivot row

		int index_pivot_row = 0;
		double temp_smallest = 0;

		for (int i = 0; i < number_of_rows; i++) {
			if (i == 0 && table[i][index_quotient_column] != 0){
				temp_smallest = table[i][index_quotient_column];
				index_pivot_row = i;
			}
			else
				if (table[i][index_quotient_column] != 0 && table[i][index_quotient_column] < temp_smallest){
					temp_smallest = table[i][index_quotient_column];
					index_pivot_row = i;
				}
		}

		cout << "The smallest quotient number is " << temp_smallest << " at index " << index_pivot_row << endl;

		//converting the value of the pivot element to 1 by deviding each element in the pivot row by pivot element

		cout << "Deviding each element in the row " << index_pivot_row << " by " << table[index_pivot_row][index_pivot_column] << endl;

		double temp_dev = table[index_pivot_row][index_pivot_column];

		for (int i = 0; i < number_of_columns - 1; i++) {
			if (table[index_pivot_row][index_pivot_column] == 0){
				//!!!!!!!!!!!!!!!!TO THINK ABOUT!!!!!!!!!!!!!!!!
				break;
			}
			table[index_pivot_row][i] /= temp_dev;
		}

		cout << "Results after devision: " << endl;
		
		for (int i = 0; i < number_of_columns - 1; i++) {
			cout << table[index_pivot_row][i] << " ";
		}

		cout << endl;


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

	double** matrix = new double*[k];//Spalten
	double *vector_b = new double[k];
	int vector_b_counter = 0;
	for (int i = 0; i < k; i++) {
		matrix[i] = new double[n];//Zeilen
		for (int j = 0; j < n; j++) {
			file >> matrix[i][j];	
		}
		file >> vector_b[vector_b_counter++];
		
	}
	//testing matrix
	cout << "matrix" << endl;
	for (int i = 0; i < k; i++) {
		for (int j = 0; j < n; j++) {
			cout << matrix[i][j] << endl;
		}
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

