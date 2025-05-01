#include <iostream>                       
#include <vector>                          
#include <iomanip>   
#include <algorithm> 
#include <chrono>

using namespace std;

void shiftMatrixStep(vector<int>& rowKs, vector<int>& colKs, int step, int n)
{

	if (step == 0){
		for (int i = 0; i < n; i++)
		{
			rowKs[i] = i;
			colKs[i] = i;
		}
	}
	else
	{
		fill(rowKs.begin(), rowKs.end(), step);
		fill(colKs.begin(), colKs.end(), step);
	}
}

void shiftMatrixRowsByKs(vector<vector<int>>& mat,
	const vector<int>& rowKs) {
	int rows = mat.size();
	if (rows == 0) return;

	for (int r = 0; r < rows; ++r) {
		int cols = mat[r].size();
		int k = rowKs[r] % cols;                
		// rotate [begin, begin+k) to the end → left shift by k
		rotate(mat[r].begin(),
			mat[r].begin() + k,
			mat[r].end());
	}

	//// Print result
	//cout << "Row-wise shifts:\n";
	//for (const auto& row : mat) {
	//	for (int x : row)
	//		cout << x << ' ';
	//	cout << '\n';
	//}
	//cout << '\n\n';
}

void shiftMatrixColsByKs(vector<vector<int>>& mat,
	const vector<int>& colKs) {
	int rows = mat.size();
	if (rows == 0) return;
	int cols = mat[0].size();


	for (int c = 0; c < cols; ++c) {
		int k = colKs[c] % rows;             
		if (k == 0) continue;

		vector<int> colData(rows);
		for (int r = 0; r < rows; ++r)
			colData[r] = mat[r][c];

		// Rotate it so that [0..k) moves to the end → upward shift by k
		rotate(colData.begin(),
			colData.begin() + k,
			colData.end());

		for (int r = 0; r < rows; ++r)
			mat[r][c] = colData[r];
	}

	//// Print result
	//cout << "Column-wise shifts:\n";
	//for (const auto& row : mat) {
	//	for (int x : row)
	//		cout << x << ' ';
	//	cout << '\n';
	//}
	//cout << '\n\n';
}


void Matrix_Mul(const vector<vector<int>>& A, const vector<vector<int>>& B, vector<vector<int>>& C, int n)

{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
				C[i][j] += A[i][j] * B[i][j];
		}
	}

}

int main() {

	int n;
	vector<vector<int>> A;
	vector<vector<int>> B;
	vector<vector<int>> C;

	char choice;

	cout << "Enter R to generate matrices randomly or another for manual input: " << endl;

	cin >> choice;

	if (choice == 'R' || choice == 'r')
	{


		srand(time(0));

		cout << "Enter the size of the Matrix: " << endl;

		cin >> n;

		while (n == 1)

		{
			cout << "You must enter size of matrix bigger than 1" << endl;

			cin >> n;

		}

		A.resize(n, vector<int>(n));
		B.resize(n, vector<int>(n));
		C.resize(n, vector<int>(n));

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				A[i][j] = rand() % 100;
				B[i][j] = rand() % 100;
			}
		}

		// cout << "Matrix A is:" << endl;

		// for (int i = 0; i < n; i++)
		// {
		// 	for (int j = 0; j < n; j++)

		// 	{
		// 		cout << setw(8) << A[i][j];
		// 	}

		// 	cout << endl;
		// }

		// cout << "Matrix B is:" << endl;

		// for (int i = 0; i < n; i++)
		// {
		// 	for (int j = 0; j < n; j++)

		// 	{
		// 		cout << setw(8) << B[i][j];
		// 	}

		// 	cout << endl;
		// }
	}

	else
	{

		cout << "Enter the size of the Matrix: " << endl;

		cin >> n;

		while (n == 1)

		{
			cout << "You must enter size of matrix bigger than 1" << endl;

			cin >> n;

		}

		cout << "Enter elements of matrix A: " << endl;

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)

				cin >> A[i][j];
		}

		// cout << "The entered matrix A is: " << endl;

		// for (int i = 0; i < n; i++)
		// {
		// 	for (int j = 0; j < n; j++)

		// 	{
		// 		cout << setw(10) << A[i][j];
		// 	}

		// 	cout << endl;
		// }

		cout << endl << "Enter elements of matrix B: " << endl;

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)

				cin >> B[i][j];
		}

		// cout << "The entered matrix B is:" << endl;

		// for (int i = 0; i < n; i++)
		// {
		// 	for (int j = 0; j < n; j++)

		// 	{
		// 		cout << setw(10) << B[i][j];
		// 	}

		// 	cout << endl;
		// }

		cout << endl;
	}

	A.resize(n, vector<int>(n));
	B.resize(n, vector<int>(n));
	C.resize(n, vector<int>(n));
	vector<vector<int>> A_shifted = A;
	vector<vector<int>> B_shifted = B;
	vector<int> rowKs(n), colKs(n);

    auto start = chrono::high_resolution_clock::now();

	for (int i = 0; i < n; i++) {

		if (i == 0){
			shiftMatrixStep(rowKs, colKs, i, n);
			shiftMatrixRowsByKs(A_shifted, rowKs);
			shiftMatrixColsByKs(B_shifted, colKs);
			Matrix_Mul(A_shifted, B_shifted, C, n);
		}

		else {
			shiftMatrixStep(rowKs, colKs, 1, n);
			shiftMatrixRowsByKs(A, rowKs);
			shiftMatrixColsByKs(B, colKs);
			Matrix_Mul(A, B, C, n);
		}
	}

    auto stop = chrono::high_resolution_clock::now();


    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);

    cout << "Took Serial " << duration.count() << " milliseconds\n";

	// for (const auto& row : C) {
	// 	for (int x : row)
	// 		cout << x << ' ';
	// 	cout << '\n';
	// }

	return 0;          
}