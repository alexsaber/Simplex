// LPSolve(OPS).cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "Funktion.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using namespace std;

bool show_log = false;


void output_basic_vars(double ** table, int k, int n){
//finding basic and non-basic variables
bool *basic_variables = new bool[k + n];
for (int i = 0; i < k + n; i++){
	basic_variables[i] = false;
}


for (int j = 0; j < k + n; j++) {//spalten
	for (int i = 0; i < k + 1; i++) {//zeilen
		if (table[i][j] != 0.0 && table[i][j] != 1.0) {
			break;
		}

		if (i == k && table[i][j] == 0)
			basic_variables[j] = true;
	}
}
cout << "<div class = 'well'>" << endl;

cout << "Basisvariablen: <b>" << endl;

for (int i = 0; i < k + n; i++){
	if (basic_variables[i] == true){
		cout << "X" << i + 1;
		if (i != k + n - 1) cout << ", ";
	}

}

cout << "</b></div>"<< endl;

}
void output_table(double ** table, int k, int n, string table_name){


	//outputting names of columns
	cout.precision(4);

	cout << "<div class='container'>" <<
		"<h4>" << table_name << "</h4>" <<
		"<table class='table table-bordered'>" <<
		" <thead>"
		"<tr>";

	for (int j = 0; j < k + n; j++) {
		cout << "<th>X" << j + 1 << "</th>";
	}

	cout << "<th>Rechte Seite</th>" << "\n" << "<th>Quotient</th>" << endl;

	cout << "\n</tr>"		<< "\n</thead>"
		<< "<tbody>";

	for (int i = 0; i < k + 1; i++) {
		cout << "<tr>";
		for (int j = 0; j < k + n + 2; j++) {
			cout  << "<td>" <<  table[i][j] << "</td>\n";
		}
		cout << " </tr>";
	}



	cout <<
		"</tbody>" <<
		"</table>" <<
		"</div>";


}


bool hasPositiveValues(double **table, int number_of_columns, int number_of_rows){
	for (int i = 0; i < number_of_columns - 2; i++) {
		if (table[number_of_rows-1][i] > 0)
			return true;
	}
	return false;
}

void sensAnalysis(int n, int k,double **basis_matrix, double **optimal_matrix, double *vector_c){
	cout << "<div class = 'well'>" << endl;
	//n = iterator_matrix
	//basis = basis_matrix
	//k = iterator_vektor

	int* bv = new int[];
	int* nbv = new int[];
	int bv_count = 0;
	int nbv_count = 0;
	double* c = new double[];
	double** a = new double*[k];
	for (int i = 0; i < n+k+1; i++){
		a[i] = new double[n+k+1];
	}

	//cout << endl;
	//cout << "Basic & Non-Basic Variables as positions. " << endl;
	//cout << "The order is determent.";
	for (int h=0,i=0,j = 0; i < n+k; i++){
		if (optimal_matrix[k][i] == 0){
			bv[h] = i;
			h++;
			bv_count++;
		}else{
			nbv[j] = i;
			j++;
			nbv_count++;
		}
	}

	/*
	*if u want to know the positions of the basic Variables.

	cout << "Basic Variables as positions:" << endl;
	for (int i = 0; i < bv_count; i++){
		cout << "Basic Variable " << i << ":	" << bv[i] << endl; 
	}

	cout << "Non Basic Variables as positions:" << endl;
	for (int i = 0; i < nbv_count; i++){
		cout << "Non Basic Variable " << i << ":	" << nbv[i] << endl; 
	}
	*/

	MatrixXd cb(bv_count,1);
	MatrixXd cnbv(nbv_count,1);
	MatrixXd B(bv_count,bv_count);
	MatrixXd Bi(bv_count,bv_count);
	MatrixXd NB(bv_count,nbv_count);

	cout << endl;

	cout << "The Vectors a[i] are being computed." << endl;
	cout << "Dependent on the position of the Basic-Variables." << endl;
	cout << "The order is determent." << endl;
	cout << "a[i]s: " << endl;
	for (int i = 0; i < n+k; i++){
		cout << "a" << i << " = [ ";
		for (int j = 0; j < k; j++){
			a[j][i] = basis_matrix[j][i];
			cout << a[j][i] << " ";
		}
		cout << "]" << endl;
	}

	
	cout << endl;

	cout << "The Basic-Variables itself." << endl;
	cout << "B: " << endl;
	for (int i = 0; i < bv_count; i++){
		for (int j = 0; j < bv_count; j++){
			B(j,i) = a[j][bv[i]];
		}
	}

	cout << B << endl;

	cout << endl;

	cout << "The invers of the Basic-Variables  is being computed." << endl;
	Bi = B.inverse();
	
	cout << "The Non-Basic-Variables itself." << endl;
	cout << "NBV: " << endl;
	for (int i = 0; i < nbv_count; i++){
		for (int j = 0; j < bv_count; j++){
			NB(j,i) = a[j][nbv[i]];
		}
	}

	cout << NB << endl;	

	cout << endl;

	cout << "The constants of the basic variables (CB) as in vector_c." << endl;
	cout << "The order is determent." << endl;
	cout << "CB: " << endl;
	for (int i = 0; i < bv_count; i ++){
		if (bv[i] >= n){
			cb(i,0) = 0.0;
		}else{
			cb(i,0) = vector_c[bv[i]];
		}
	}
	cout << " [ ";
	for (int i = 0; i < bv_count; i++){
		cout << cb(i,0) << " " ;
	}
	cout << "]" << endl;	

	cout << endl;

	cout << "The Constants of the non-basic-variables (CNBV) as in vector_c." << endl;
	cout << "The order is determent." << endl;
	cout << "CNBV: " << endl;
	for (int i = 0; i < nbv_count; i ++){
		if (nbv[i] >= n){
			cnbv(i,0)= 0.0;
		}else{
			cnbv(i,0) = vector_c[nbv[i]];
		}
	}

	cout << " [ ";
	for (int i = 0; i < nbv_count; i++){
		cout << cnbv(i,0) << " " ;
	}
	cout << "]" << endl;	

	cout << endl;

	//Changing the Object function coefficiant of a basic variable
	double cbv_dot_Biaj = 0.0;
	double cnbv_dot_Biaj = 0.0;
	double cbv0 = 0.0;
	double cnbv0 = 0.0;
	double a1 = 0.0;
	double b2 = 0.0;
	double* temp = new double[nbv_count];
	for (int i = 0; i < nbv_count; i++){
		temp[i] = 0.0;
	}
	double delta = 0.0;

	cout << endl;
	for (int i = 0; i < bv_count; i ++){
		cbv_dot_Biaj = 0.0;
		cbv0 = 0.0;
		a1 = 0.0;
		b2 = 0.0;
		for (int m = 0; m < bv_count; m++){
			temp[m] = 0.0;
		}
		delta = 0.0;
		for (int j = 0; j < bv_count; j++){
			for (int m = 0; m < bv_count; m++){
				a1 = a[m][i];
				b2 = Bi(j,m);
				temp[j] += a1*b2;
			}
			if (j < nbv_count){
				cbv0 = cnbv(j,0);
			}else{
				cbv0 = 0.0;
			}
			cbv_dot_Biaj += cbv0*temp[j];
		}
		while (cbv_dot_Biaj + delta <= 0){
			delta += 1;
		};
		
		cout << "The Coefficiants of BV" << i << " can be added up by a maximum of  " << delta << "." << endl;
		
	}

	cout << endl;

	//Changing the Object function coefficiant of a non basic variable
	cnbv0 = 0.0;
	a1 = 0.0;
	b2 = 0.0;
	for (int i = 0; i < nbv_count; i++){
		temp[i] = 0.0;
	}
	delta = 0.0;
	cout << endl;
	for (int i = 0; i < nbv_count; i ++){
		cnbv_dot_Biaj = 0.0;
		cnbv0 = 0.0;
		a1 = 0.0;
		b2 = 0.0;
		for (int m = 0; m < nbv_count; m++){
			temp[m] = 0.0;
		}
		delta = 0.0;
		for (int j = 0; j < bv_count; j++){
			for (int m = 0; m < bv_count; m++){
				a1 = a[m][i];
				b2 = Bi(j,m);
				temp[j] += a1*b2;
			}
			if (j < nbv_count){
				cnbv0 = cnbv(j,0);
			}else{
				cnbv0 = 0.0;
			}
			cnbv_dot_Biaj += cnbv0*temp[j];
		}
		while (cnbv_dot_Biaj - cb(i,0) - delta > 0){
			delta += 1;
		};
		
		cout << "The Coefficiants of NBV" << i << " can be added up by a maximum of  " << delta << "." << endl;
		
	}

	cout << endl;
	
	//Adding a new variable or activity (20,[1,2,3,4])
	bool test = false;
	cnbv_dot_Biaj = 0.0;
	cnbv0 = 0.0;
	a1 = 0.0;
	b2 = 0.0;
	for (int i = 0; i < nbv_count; i++){
		temp[i] = 0.0;
	}
	delta = 0.0;

	for (int j = 0; j < bv_count; j++){
		for (int m = 0; m < bv_count; m++){
			a1++;
			b2 = Bi(j,m);
			temp[j] += a1*b2;
		}
		if (j < nbv_count){
		cnbv0 = cnbv(j,0);
		}else{
			cnbv0 = 0.0;
		}
	cnbv_dot_Biaj += cnbv0*temp[j];
	}
	if (cnbv_dot_Biaj - 20 > 0){
		cout << "Adding a new activity <20> with Recources" << endl << " [ ";
		for (int i = 0; i<k;i++){
			cout << i+1 << " ";
		}
		cout << "]" << endl << "is NOT Reasonable." << endl ;
	}
	else{
		cout << "Adding a new activity <20> with Recources" << endl << " [ ";
		for (int i = 0; i<k;i++){
			cout << i+1 << " ";
			test = true;
		}
		cout << "]" << endl << "is REASONABLE." << endl ;
	}

	cout << endl;
	
	//Testen welche Aktivität sinnvoll wäre
if (!test){
	delta = 1.0;
	while (cnbv_dot_Biaj - delta > 0 && delta != 0.0){
		delta += 1;
	}
	cout << "Adding a new activity <" << delta << "> with Recources" << endl << " [ ";
	for (int i = 0; i<k;i++){
		cout << i+1 << " ";
	}
	cout << "]" << endl << "is REASONABLE." << endl ;

	cout << endl;

	//Testen welche Resourcen sinnvoll wären
	int a2 = 0;
	int temp_a = 0;
	double* temp1 = new double[bv_count];
	cnbv_dot_Biaj = 0.0;
	cnbv0 = 0.0;
	b2 = 0.0;
	for (int i = 0; i < nbv_count; i++){
		temp[i] = 0.0;
	}

	for (int i = 0; i < bv_count; i++){
		temp1[i] = 0.0;
	}
	while (cnbv_dot_Biaj - 20 < 0){
		a2 += 1;
		for (int j = 0; j < bv_count; j++){
			switch (a2%k){
				case 0: temp_a = 0; break;
				case 1: temp_a = 1; break;
				case 2: temp_a = 2; break;
				case 3: temp_a = 3; break;
				default: temp_a = a2%k; break;
			}
			for (int m = 0; m < bv_count; m++){
				b2 = Bi(j,m);
				temp[j] += a2*b2;
				temp1[temp_a] = a2;
			}
			if (j < nbv_count){
			cnbv0 = cnbv(j,0);
			}else{
				cnbv0 = 0.0;
			}
			cnbv_dot_Biaj += cnbv0*temp[j];
		}
	}
	cout << "To add the activity <20> one would have to add the resource" << endl << " [ ";
	for (int i = 0; i < bv_count; i++){
		cout << temp1[i] << " ";
	}
	cout << "]." << endl;
}
	cout << "</div>" << endl;
}




double* lpsolve(int n, double *c, int k, double **A, double *b){
	//Schlupfvariablen, befuellen mit '1'
	int *s_vars = new int[k];
	for (int i = 0; i < k; i++) {
		s_vars[i] = 1;
	}

	//create table and fill up with the 
	double** table = new double*[k + 1];//Spalten
	double** basis_table = new double*[k + 1];//Spalten
	for (int i = 0; i < k+1; i++){
		basis_table[i] = new double[k + n + 2];//Zeilen
	}

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
	
	//multiplying Zielfukntionargumente with -1 (was, now just copying)
	for (int i = 0; i < n; i++){
		table[k][i] = c[i];
	}
	
	//fillng up Schluepfvariablen for last row with zeros
	for (int i = n; i < n + k + 1; i++){
		table[k][i] = 0;
	}

	//filling up Quotient column with zeros
	for (int i = 0; i < k + 1; i++){
		table[i][k + n + 1] = 0;
	}

	//filling up basic_table

	for (int i = 0; i < k + 1; i++) {
		for (int j = 0; j < k + n + 2; j++) {
			basis_table[i][j] = table[i][j];
		}
	}

	//continuing with outputting restrictions

	cout << "Die Restriktionen: <b>";
	for (int i = 0; i < k; i++) {//zeilen
		cout << "<br>" ;
		for (int j = 0; j < n + 1; j++) {//spalten
			if (j < n - 1)
				cout << "x" << i + 1 << "*" << table[i][j] << " + ";
			if (j == n - 1)
				cout << "x" << i + 1 << "*" << table[i][j];
			if (j == n)
				cout << " <= " << table[i][j + n];
		}
	}
	cout << "</b>";
	//close div of first well
	cout << "</div>";




	output_table(table, k, n, "Ausgangstableau:");

	output_basic_vars(table, k, n);

	//cout << "RS steht fuer die \"rechte Seite\" des (urspruenglichen) Gleichungssystems." << endl;
	cout << endl;
	//----BEGIN CALCULATION-----

	cout.precision(4);


	int number_of_rows = k + 1;
	int number_of_columns = k + n + 2;
	int index_quotient_column = number_of_columns - 1;
	int index_rechte_seite_column = number_of_columns - 2;

	if (show_log) cout << "number_of_rows: " << number_of_rows << endl;
	if (show_log) cout << "number_of_columns: " << number_of_columns << endl;

	int iteration_number = 1;

	//find the highest number (without "-") in the last row OR find pivotcolumn
do{
	int index_pivot_column = 0;
	double temp_highest = 0;

	for (int i = 0; i < index_quotient_column; i++) {
		if (i == 0){
			temp_highest = table[number_of_rows - 1][i];
			index_pivot_column = i;
		}
		else
			if (table[number_of_rows - 1][i] > 0 && temp_highest < table[number_of_rows - 1][i]){
				temp_highest = table[number_of_rows - 1][i];
				index_pivot_column = i;
			}
	}
	cout << "<div class = 'well'>" << endl;
	cout << "Die groesste positive Zahl in der ZF-Spalte ist " << temp_highest << "<br>" << endl;
	cout << "Spalte Nummer " << index_pivot_column + 1 << " ist die Pivotspalte" << "<br>" << endl;
	cout << endl;
	//calculate Quotient

	cout << "Kalkuliere Quotient" << "<br>" << endl;

	cout << "</div>" << endl;

	for (int i = 0; i < number_of_rows - 1; i++) {
		if (table[i][index_pivot_column] != 0){
			table[i][index_quotient_column] = table[i][number_of_columns - 2] / table[i][index_pivot_column];
		}
		else
			table[i][index_quotient_column] = -1;
	}


	



	//find the smallest number in the Quotient column OR find pivot row

	int index_pivot_row = 0;
	double temp_smallest = 0; //what if it can't be set at all???
	bool first_set = false;

	for (int i = 0; i < number_of_rows; i++) {
		if (first_set == false && table[i][index_quotient_column] > 0){
			temp_smallest = table[i][index_quotient_column];
			index_pivot_row = i;
			first_set = true;
		}
		else
			if (table[i][index_quotient_column] > 0 && table[i][index_quotient_column] < temp_smallest){
				temp_smallest = table[i][index_quotient_column];
				index_pivot_row = i;
			}
	}


	//cout << "Tableau nach Kalkulation (der kleinste Quotien ist mit >*< hervorgehoben):" << endl;

	//old, before html
	/*
	for (int j = 0; j < k + n; j++) {
		cout << "X" << j + 1 << "\t";
	}
	cout << "RS" << "\t" << "Q" << endl;

	for (int i = 0; i < number_of_rows; i++) {
		for (int j = 0; j < number_of_columns; j++) {
			if (j == index_quotient_column && table[i][j] == temp_smallest)
				cout << ">" << table[i][j] << "<" << "\t";
			else
				cout << table[i][j] << "\t";

		}
		cout << endl;
	}
	*/
	//////////////////////////////////////////////////



	//outputting names of columns
	cout.precision(4);

	cout << "<div class='container'>" <<
		"<h4>" << "Tableau nach Kalkulation(der kleinste Quotient ist rot)" << "</h4>" <<
		"<table class='table table-bordered'>" <<
		" <thead>"
		"<tr>";

	for (int j = 0; j < k + n; j++) {
		cout << "<th>X" << j + 1 << "</th>";
	}

	cout << "<th>Rechte Seite</th>" << "\n" << "<th>Quotient</th>" << endl;

	cout << "\n</tr>"		<< "\n</thead>"
		<< "<tbody>";

	for (int i = 0; i < k + 1; i++) {
		cout << "<tr>";
		for (int j = 0; j < k + n + 2; j++) {
			if (j == index_quotient_column && table[i][j] == temp_smallest)
				cout << "<td><font color='red'>" << table[i][j] << "</font></td>\n";
			else
				cout << "<td>" << table[i][j] << "</td>\n";
			
		}
		cout << " </tr>";
	}



	cout <<
		"</tbody>" <<
		"</table>" <<
		"</div>";








	////////////////////////////////////////////////////




	cout << endl;
	/*
	cout << "Jetzt kann man das Pivotelement bestimmen. (ist mit >*< hervorgehoben)" << endl;

	for (int j = 0; j < k + n; j++) {
		cout << "X" << j + 1 << "\t";
	}
	cout << "RS" << "\t" << "Q" << endl;

	for (int i = 0; i < number_of_rows; i++) {
		for (int j = 0; j < number_of_columns; j++) {
			if (index_pivot_row == i && index_pivot_column == j)
				cout << ">" << table[i][j] << "<" << "\t";
			else
				cout << table[i][j] << "\t";

		}
		cout << endl;
	}
	*/
	//////////////////////////////////////////////////



	//outputting names of columns
	cout.precision(4);

	cout << "<div class='container'>" <<
		"<h4>" << "Das Pivotelement. (ist rot)" << "</h4>" <<
		"<table class='table table-bordered'>" <<
		" <thead>"
		"<tr>";

	for (int j = 0; j < k + n; j++) {
		cout << "<th>X" << j + 1 << "</th>";
	}

	cout << "<th>Rechte Seite</th>" << "\n" << "<th>Quotient</th>" << endl;

	cout << "\n</tr>"		<< "\n</thead>"
		<< "<tbody>";

	for (int i = 0; i < k + 1; i++) {
		cout << "<tr>";
		for (int j = 0; j < k + n + 2; j++) {
			if (index_pivot_row == i && index_pivot_column == j)
				cout << "<td><font color='red'>" << table[i][j] << "</font></td>\n";
			else
				cout << "<td>" << table[i][j] << "</td>\n";

		}
		cout << " </tr>";
	}



	cout <<
		"</tbody>" <<
		"</table>" <<
		"</div>";








	////////////////////////////////////////////////////




	cout << endl;
	//converting the value of the pivot element to 1 by deviding each element in the pivot row by pivot element
	cout << "<div class = 'well'>" << endl;
	cout << "Dividiere alle Werte in der Pivotspalte " << index_pivot_row+1 << " durch das Pivotelement mit Wert " << table[index_pivot_row][index_pivot_column] << endl;
	cout << "</div>" << endl;
	double temp_dev = table[index_pivot_row][index_pivot_column];

	for (int i = 0; i < number_of_columns - 1; i++) {
		if (table[index_pivot_row][index_pivot_column] == 0){
			//!!!!!!!!!!!!!!!!TO THINK ABOUT!!!!!!!!!!!!!!!!
			break;
		}
		table[index_pivot_row][i] /= temp_dev;
	}


	output_table(table, k, n, "Die Tabelle nach der Division");

	cout << endl;
	cout << "<div class = 'well'>" << endl;
	cout << "Mit Hilfe vom Gaussschen Eliminationsverfahren alle Werte in der Pivotspalte zu 0 bringen (ausser das Privotelement)" << endl;
	cout << "</div>" << endl;
	for (int i = 0; i < number_of_rows; i++) {
		if (i != index_pivot_row && table[i][index_pivot_column] != 0){
			double temp_multiplyer = table[i][index_pivot_column];
			for (int j = 0; j < number_of_columns - 1; j++) {
				table[i][j] = table[i][j] - (table[index_pivot_row][j] * temp_multiplyer);
			}
		}

	}



	output_table(table, k, n, "Das Tableu nachdem:");

	cout << endl;

	//testing tabelle
	//cout << "Interation Nummer " << iteration_number++ << " is abgeschlossen. Zeige das Tableau" << endl;
	//cout.precision(4);

	string interation_first = "Die Tabelle nach der Iteration Nr. " + std::to_string(iteration_number);
	iteration_number++;
	output_table(table, k, n, interation_first);

	output_basic_vars(table, k, n);

	if (show_log) cout << "hasPositiveValues(table, number_of_columns, number_of_rows) returned: " << hasPositiveValues(table, number_of_columns, number_of_rows) << endl;

}while (hasPositiveValues(table, number_of_columns, number_of_rows));



	double* results_first = new double[n];
{

	
	for (int i = 0; i < n; i++){
		results_first[i] = 0;
	}
	int results_counter = 0;

	for (int i = 0; i < n; i++) {
		int row_index_where_ONE_was_found = -1;
		bool not_only_zeros_and_ONE = false;
		for (int j = 0; j < number_of_rows; j++) {
			if (table[j][i] == 1){
				row_index_where_ONE_was_found = j;
			}
			else if (table[j][i] != 0){
				not_only_zeros_and_ONE = true;
			}
		}
		if (not_only_zeros_and_ONE != true && row_index_where_ONE_was_found >= 0){
			results_first[results_counter++] = table[row_index_where_ONE_was_found][number_of_columns - 2];
		}
		else
			results_counter++;
	}


	//results
	cout << "<div class = 'well'>" << endl;

	cout << endl << "<h3 style='color:red;'>Eine optimale Loesung ist gefunden! </h3>   <h3> <b>" << endl;
	for (int i = 0; i < n; i++) {
		if (i < n - 1)
			cout << "x" << i + 1 << " = " << results_first[i] << ", ";
		else
			cout << "x" << i + 1 << " = " << results_first[i];

	}
	cout << endl << endl;
	cout.precision(4);
	cout << "</b> </h3> <br> <h4> Somit liefert die Zielfunktion den Wert <b>" << table[number_of_rows-1][index_rechte_seite_column] * (-1) << "</b>" << endl;
	
}

cout << "</h4></div>" << endl;
	cout << endl;









	//////////////////////////////////////////////
	//finding basic and non-basic variables
	bool *basic_variables = new bool[k + n];
	for (int i = 0; i < k + n; i++){
		basic_variables[i] = false;
	}


	for (int j = 0; j < k + n; j++) {//spalten
		for (int i = 0; i < k + 1; i++) {//zeilen
			if (table[i][j] != 0.0 && table[i][j] != 1.0) {
				break;
			}
				
			if (i == k && table[i][j] == 0) 
				basic_variables[j] = true;	
		}
	}


	///////////////////////////////////////////////////////////////
	//begin of testing if a non-basic variable has 0 in Zielfunktion
	if (show_log) cout << "Testing if there are multiple solutions" << endl;

	int *index_column_another_solution = nullptr;

	for (int i = 0; i < k + n; i++){
		if (basic_variables[i] == false && (table[k][i] == 0)){
			index_column_another_solution = new int;
			*index_column_another_solution = i;
		}
	}
	cout << endl;
	
	if (index_column_another_solution == nullptr){
		cout << "<div class = 'well'><h4>" << endl;
		cout << "Dieses LP hat nur eine Loesung." << endl;
		cout << "<h4></div>" << endl;
	}
		
	else{
		cout << "<div class = 'well'>" << endl;
		cout << "<h3 style='color:green;'>Dieses LP hat mehrere Loesungen!</h3>" << endl << 
			"<h4><b>x" << *index_column_another_solution+1 << " </b> kann zu Basisbariable umgerechnet werden, ohne den Zielfunktionswert zu benachteiligen" << endl;
		cout << "</h4></div>" << endl;
		iteration_number = 1;

		/////////////////////////////////////////////////////////
		//calculate the coordinates of the new point


		//code from above START

		//calculate Quotient
		cout << "<div class = 'well'>" << endl;
		cout << "Kalkuliere Quotient" << "<br>" << endl;

		cout << "</div>" << endl;

		for (int i = 0; i < number_of_rows - 1; i++) {
			if (table[i][*index_column_another_solution] != 0){
				table[i][index_quotient_column] = table[i][number_of_columns - 2] / table[i][*index_column_another_solution];
			}
			else
				table[i][index_quotient_column] = -1;
		}






		//find the smallest number in the Quotient column OR find pivot row

		int index_pivot_row = 0;
		double temp_smallest = 0; //what if it can't be set at all???
		bool first_set = false;

		for (int i = 0; i < number_of_rows; i++) {
			if (first_set == false && table[i][index_quotient_column] > 0){
				temp_smallest = table[i][index_quotient_column];
				index_pivot_row = i;
				first_set = true;
			}
			else
				if (table[i][index_quotient_column] > 0 && table[i][index_quotient_column] < temp_smallest){
					temp_smallest = table[i][index_quotient_column];
					index_pivot_row = i;
				}
		}


		//cout << "Tableau nach Kalkulation (der kleinste Quotien ist mit >*< hervorgehoben):" << endl;

		//old, before html
		/*
		for (int j = 0; j < k + n; j++) {
		cout << "X" << j + 1 << "\t";
		}
		cout << "RS" << "\t" << "Q" << endl;

		for (int i = 0; i < number_of_rows; i++) {
		for (int j = 0; j < number_of_columns; j++) {
		if (j == index_quotient_column && table[i][j] == temp_smallest)
		cout << ">" << table[i][j] << "<" << "\t";
		else
		cout << table[i][j] << "\t";

		}
		cout << endl;
		}
		*/
		//////////////////////////////////////////////////



		//outputting names of columns
		cout.precision(4);

		cout << "<div class='container'>" <<
			"<h4>" << "Tableau nach Kalkulation(der kleinste Quotient ist rot)" << "</h4>" <<
			"<table class='table table-bordered'>" <<
			" <thead>"
			"<tr>";

		for (int j = 0; j < k + n; j++) {
			cout << "<th>X" << j + 1 << "</th>";
		}

		cout << "<th>Rechte Seite</th>" << "\n" << "<th>Quotient</th>" << endl;

		cout << "\n</tr>"			<< "\n</thead>"
			<< "<tbody>";

		for (int i = 0; i < k + 1; i++) {
			cout << "<tr>";
			for (int j = 0; j < k + n + 2; j++) {
				if (j == index_quotient_column && table[i][j] == temp_smallest)
					cout << "<td><font color='red'>" << table[i][j] << "</font></td>\n";
				else
					cout << "<td>" << table[i][j] << "</td>\n";

			}
			cout << " </tr>";
		}



		cout <<
			"</tbody>" <<
			"</table>" <<
			"</div>";








		////////////////////////////////////////////////////




		cout << endl;
		/*
		cout << "Jetzt kann man das Pivotelement bestimmen. (ist mit >*< hervorgehoben)" << endl;

		for (int j = 0; j < k + n; j++) {
		cout << "X" << j + 1 << "\t";
		}
		cout << "RS" << "\t" << "Q" << endl;

		for (int i = 0; i < number_of_rows; i++) {
		for (int j = 0; j < number_of_columns; j++) {
		if (index_pivot_row == i && index_pivot_column == j)
		cout << ">" << table[i][j] << "<" << "\t";
		else
		cout << table[i][j] << "\t";

		}
		cout << endl;
		}
		*/
		//////////////////////////////////////////////////



		//outputting names of columns
		cout.precision(4);

		cout << "<div class='container'>" <<
			"<h4>" << "Das Pivotelement. (ist rot)" << "</h4>" <<
			"<table class='table table-bordered'>" <<
			" <thead>"
			"<tr>";

		for (int j = 0; j < k + n; j++) {
			cout << "<th>X" << j + 1 << "</th>";
		}

		cout << "<th>Rechte Seite</th>" << "\n" << "<th>Quotient</th>" << endl;

		cout << "\n</tr>"			<< "\n</thead>"
			<< "<tbody>";

		for (int i = 0; i < k + 1; i++) {
			cout << "<tr>";
			for (int j = 0; j < k + n + 2; j++) {
				if (index_pivot_row == i && *index_column_another_solution == j)
					cout << "<td><font color='red'>" << table[i][j] << "</font></td>\n";
				else
					cout << "<td>" << table[i][j] << "</td>\n";

			}
			cout << " </tr>";
		}



		cout <<
			"</tbody>" <<
			"</table>" <<
			"</div>";








		////////////////////////////////////////////////////




		cout << endl;
		//converting the value of the pivot element to 1 by deviding each element in the pivot row by pivot element
		cout << "<div class = 'well'>" << endl;
		cout << "Dividiere alle Werte in der Pivotspalte " << index_pivot_row + 1 << " durch das Pivotelement mit Wert " << table[index_pivot_row][*index_column_another_solution] << endl;
		cout << "</div>" << endl;
		double temp_dev = table[index_pivot_row][*index_column_another_solution];

		for (int i = 0; i < number_of_columns - 1; i++) {
			if (table[index_pivot_row][*index_column_another_solution] == 0){
				//!!!!!!!!!!!!!!!!TO THINK ABOUT!!!!!!!!!!!!!!!!
				break;
			}
			table[index_pivot_row][i] /= temp_dev;
		}


		output_table(table, k, n, "Die Tabelle nach der Division");

		cout << endl;
		cout << "<div class = 'well'>" << endl;
		cout << "Mit Hilfe vom Gaussschen Eliminationsverfahren alle Werte in der Pivotspalte zu 0 bringen (ausser das Privotelement)" << endl;
		cout << "</div>" << endl;
		for (int i = 0; i < number_of_rows; i++) {
			if (i != index_pivot_row && table[i][*index_column_another_solution] != 0){
				double temp_multiplyer = table[i][*index_column_another_solution];
				for (int j = 0; j < number_of_columns - 1; j++) {
					table[i][j] = table[i][j] - (table[index_pivot_row][j] * temp_multiplyer);
				}
			}

		}



		output_table(table, k, n, "Das Tableu nachdem:");

		cout << endl;

		//testing tabelle
		//cout << "Interation Nummer " << iteration_number++ << " is abgeschlossen. Zeige das Tableau" << endl;
		//cout.precision(4);

		output_table(table, k, n, "Die Tabelle nach der Iteration");

		output_basic_vars(table, k, n);







		//code from above END


		double* results_second = new double[n];
		for (int i = 0; i < n; i++){
			results_second[i] = 0;
		}
		int results_counter = 0;

		for (int i = 0; i < n; i++) {
			int row_index_where_ONE_was_found = -1;
			bool not_only_zeros_and_ONE = false;
			for (int j = 0; j < number_of_rows; j++) {
				if (table[j][i] == 1){
					row_index_where_ONE_was_found = j;
				}
				else if (table[j][i] != 0){
					not_only_zeros_and_ONE = true;
				}
			}
			if (not_only_zeros_and_ONE != true && row_index_where_ONE_was_found >= 0){
				results_second[results_counter++] = table[row_index_where_ONE_was_found][number_of_columns - 2];
			}
			else
				results_counter++;
		}
		cout << "<div class = 'well'>" << endl;
		cout << endl << "<h3 style='color:red;'> Die zweite optimale Loesung ist gefunden! </h3> <h3> Folgender Ausdruck zeigt alle moeglichen Loesungen: <br> <b>" << endl;

		cout << "(";

		for (int i = 0; i < n; i++){
			cout << "x" << i + 1;
			if (i != n - 1) cout << ", ";
		}

		cout << ") = delta * ";

		cout << "(";

		for (int i = 0; i < n; i++){
			cout << results_first[i];
			if (i != n - 1) cout << ", ";
		}

		cout << ") + (1 - delta) * (";

		for (int i = 0; i < n; i++){
			cout << results_second[i];
			if (i != n - 1) cout << ", ";
		}

		cout << "), wo delta beliebigen Wert innerhalb des Intervals [0, 1] hat" << endl;
		

		cout << endl;
		cout << "</b> <br> <br> Die Zielfunktion wird immer den Wert " << table[number_of_rows - 1][index_rechte_seite_column] * (-1) << " zurueckliefern.</h3>" << endl;

		cout << "</div>" << endl;

	}

	//sensitivity

	sensAnalysis(n, k, basis_table, table, c);

	return results_first;
}



int main(int argc, char* argv[]){

	cout << "<!DOCTYPE html>\n"
		<< "<html lang='de'>\n"
		<< "<head>\n"
		<< "<title>Bootstrap Example</title>\n"
		<< "<meta charset='utf-8'>\n"
		<< "<meta name='viewport' content='width=device-width, initial-scale=1'>\n"
		<< "<link rel='stylesheet' href='http://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap.min.css'>\n"
		<< "<script src='https://ajax.googleapis.com/ajax/libs/jquery/1.11.3/jquery.min.js'></script>\n"
		<< "<script src='http://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/js/bootstrap.min.js'></script>\n"
		<< "</head>\n"
		<< "<body>\n"

		<< "<div class='container'>\n" ;






	
	if (argc < 2){
		cout << "<p><div class = 'alert alert-danger' role = 'alert'>Man muss eine *.txt Datei als erster Parameter uebergeben!</p>"<< endl;

		cout << "</div>";

		cout << "</body>\n";
		cout << "</html>\n";

		return 1;
	}

	ifstream file(argv[1]);
	if (!file.good()) {
		cout << "<p><div class = 'alert alert-danger' role = 'alert'>Etwas stimmt mit der Datei nicht!</p>" << endl;

		cout << "</div>";

		cout << "</body>\n";
		cout << "</html>\n";

		return 1;
	}

	std::string file_str = argv[1];
	if (file_str.find(".txt") != std::string::npos) {
		cout << "<p><div class = 'alert alert-success' role = 'alert'>Die datei gefunden!</p>" << endl;
	}
	else{
		cout << "<p><div class = 'alert alert-danger' role = 'alert'>Das soll eine txt Datei sein!</p>" << endl;

		cout << "</div>";

		cout << "</body>\n";
		cout << "</html>\n";

		return 1;
	}
	show_log = false;
	//bool min_active = false;
	for (int i = 1; i < argc; i++){
		if (show_log) cout << argv[i] << endl;
		std::string arg_str = argv[i];
		std::string log_str = "show_log";
		if (arg_str == log_str) {
			cout << "Log-Infos sind aktiv" << endl;
			show_log = true;
		}
		/*
		if (arg_str.compare("min")) {
			min_active = true;
		}
		*/
	}
	/*
	if (min_active)
		cout << "Die Zielfunktion wird minimiert" << endl;
	else
		cout << "Die Zielfunktion wird maximiert" << endl;
	*/

	/*
	string testfile;
	ifstream file("testfile.txt");
	
	bool input = false;
	while(!input){
	try{
		cout << "Name of Input file without format (this has to be a .txt file): ";
		cin >> testfile;
		ifstream file(testfile + ".txt");
		if (testfile[0] != 't' && testfile[1] != 'e' && testfile[2] != 's' && testfile[3] != 't'){
			throw 1;
		}
		input = !input;
	}catch(int e){
		cout << "Der Name der Datei muss mit <test> beginnen." << endl;
	}
	}
	*/


	cout << "</div>";

	cout << "<div class = 'well'>" << endl;

	int n, k; // n Zeilen, k Spalten
	file >> n >> k;
	if (show_log) cout << "n: " << n << ", k: " << k << endl;

	cout << "<p>Vorgegeben sind <b>" << n << "</b> Variablen mit <b>" << k << "</b> Restriktionen</p>" << endl;



	double *vector_c = new double[n];
	for (int i = 0; i < n; i++) {
		file >> vector_c[i];
	}
	/*
	if (min_active){
		for (int i = 0; i < n; i++) {
			double value = vector_c[i];
			vector_c[i] = value * (-1.0);
		}
	}
	*/

	//testing vector_c
	if (show_log){
		cout << "vector_c" << endl;
		for (int i = 0; i < n; i++) {
			cout << vector_c[i] << endl;
		}
	}

	cout << "<p>Die Zielfunktion:" << endl
		<< "<br><b>";
		for (int i = 0; i < n; i++) {
			if (i != n - 1)
				cout << "x" << i + 1 << "*" << vector_c[i] << " + ";
			else
				cout << "x" << i + 1 << "*" << vector_c[i];
		}
		cout << "</b></p>";



	

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
	if (show_log){
		cout << "matrix" << endl;
		for (int i = 0; i < k; i++) {
			for (int j = 0; j < n; j++) {
				cout << matrix[i][j] << endl;
			}
		}
	}
	


	//testing vector_b
	if (show_log){
		cout << "vector_b" << endl;
		for (int i = 0; i < k; i++) {
			cout << vector_b[i] << endl;
		}
	}

	//testing if there are negativ numbers

	for (int i = 0; i < k; i++) {
		if (vector_b[i] < 0){

			cout << "<p><div class = 'alert alert-danger' role = 'alert'>Es gilt Nichtnegativitätsbedingung, deswegen dürfen sich auf der rechten Sei der Ungleichungen mit keine negativen Zahlen vorkommen!</p>" << endl;

			cout << "</div>";

			cout << "</body>\n";
			cout << "</html>\n";

			return 1;
		}
	}


	//cout << "</div>";

	//results

	double* results = lpsolve(n, vector_c, k, matrix, vector_b);

	/*
	cout << endl << "Die optimale Loesung ist: " << endl;
	for (int i = 0; i < n; i++) {
		cout << "x" << i+1 << " = " << results[i] << ", ";
	}
	cout << endl << endl;
	*/

	//system("pause");








	cout << "</div> </body>\n";
	cout << "</html>\n";

	return 0;
}
