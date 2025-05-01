#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <chrono>
#include <omp.h>

using namespace std;

void shiftMatrixStep(vector<int>& rowKs, vector<int>& colKs, int step, int n) {
	if (step == 0) {
		for (int i = 0; i < n; i++) {
			rowKs[i] = i;
			colKs[i] = i;
		}
	} else {
		fill(rowKs.begin(), rowKs.end(), 1); // shift by 1 for next steps
		fill(colKs.begin(), colKs.end(), 1);
	}
}

void shiftMatrixRowsByKs(vector<vector<int>>& mat, const vector<int>& rowKs) {
	int rows = mat.size();
	if (rows == 0) return;

#pragma omp parallel for
	for (int r = 0; r < rows; ++r) {
		int cols = mat[r].size();
		int k = rowKs[r] % cols;
		rotate(mat[r].begin(), mat[r].begin() + k, mat[r].end());
	}
}

void shiftMatrixColsByKs(vector<vector<int>>& mat, const vector<int>& colKs) {
	int rows = mat.size();
	if (rows == 0) return;
	int cols = mat[0].size();

#pragma omp parallel for
	for (int c = 0; c < cols; ++c) {
		int k = colKs[c] % rows;
		if (k == 0) continue;

		vector<int> colData(rows);
		for (int r = 0; r < rows; ++r)
			colData[r] = mat[r][c];

		rotate(colData.begin(), colData.begin() + k, colData.end());

		for (int r = 0; r < rows; ++r)
			mat[r][c] = colData[r];
	}
}

void Matrix_Mul(const vector<vector<int>>& A, const vector<vector<int>>& B, vector<vector<int>>& C, int n) {
#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			C[i][j] += A[i][j] * B[i][j];
		}
	}
}

int main() {
	int n, numThreads;
	cout << "Enter number of threads: ";
	cin >> numThreads;
	omp_set_num_threads(numThreads);

	vector<vector<int>> A, B, C;
	char choice;

	cout << "Enter R to generate matrices randomly or another for manual input: ";
	cin >> choice;

	cout << "Enter the size of the Matrix: ";
	cin >> n;
	while (n == 1) {
		cout << "You must enter size of matrix bigger than 1" << endl;
		cin >> n;
	}

	A.resize(n, vector<int>(n));
	B.resize(n, vector<int>(n));
	C.assign(n, vector<int>(n, 0));

	if (choice == 'R' || choice == 'r') {
		srand(time(0));
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = rand() % 100;
				B[i][j] = rand() % 100;
			}
		}
	} else {
		cout << "Enter matrix A:\n";
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				cin >> A[i][j];

		cout << "Enter matrix B:\n";
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				cin >> B[i][j];
	}

    // cout << "Matrix A is:" << endl;

    // for (int i = 0; i < n; i++)
    // {
    //     for (int j = 0; j < n; j++)

    //     {
    //         cout << setw(8) << A[i][j];
    //     }

    //     cout << endl;
    // }

    // cout << "Matrix B is:" << endl;

    // for (int i = 0; i < n; i++)
    // {
    //     for (int j = 0; j < n; j++)

    //     {
    //         cout << setw(8) << B[i][j];
    //     }

    //     cout << endl;
    // }

	auto start = chrono::high_resolution_clock::now();

	vector<int> rowKs(n), colKs(n);
	// Initial shift
	shiftMatrixStep(rowKs, colKs, 0, n);
	shiftMatrixRowsByKs(A, rowKs);
	shiftMatrixColsByKs(B, colKs);

	for (int step = 0; step < n; step++) {
		Matrix_Mul(A, B, C, n);

		// Prepare for next shift
		if (step < n - 1) {
			shiftMatrixStep(rowKs, colKs, 1, n);
			shiftMatrixRowsByKs(A, rowKs);
			shiftMatrixColsByKs(B, colKs);
		}
	}

	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);

	cout << "Parallelized Cannon's Algorithm took " << duration.count() << " milliseconds\n";
	// cout << "Resultant Matrix C:\n";
	// for (const auto& row : C) {
	// 	for (int x : row)
	// 		cout << setw(8) << x;
	// 	cout << '\n';
	// }

	return 0;
}
