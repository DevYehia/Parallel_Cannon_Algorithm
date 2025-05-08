#include <iostream>
#include <vector>
#include <algorithm>   // rotate, fill
#include <cmath>       // sqrt, ceil, floor
#include <iomanip>     // setw
//#include <chrono>      // system_clock
//#include <random>      // mt19937, uniform_int_distribution

using namespace std;

// A Block is a blockSize x blockSize matrix of ints
using Block = vector<vector<int>>;
// A Grid is gridSize rows of gridSize Blocks
using Grid = vector<vector<Block>>;

// Print an N x N matrix
void printMatrix(const vector<vector<int>>& matrix) {
    int N = matrix.size();
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < N; ++c) {
            cout << setw(6) << matrix[r][c];
        }
        cout << "\n";
    }
    cout << "\n";
}

// Pad an N x N matrix up to paddedSize x paddedSize with zeros
vector<vector<int>> padMatrix(const vector<vector<int>>& matrix,
    int paddedSize) {
    int N = matrix.size();
    vector<vector<int>> result(paddedSize, vector<int>(paddedSize, 0));
    for (int r = 0; r < N; ++r)
        for (int c = 0; c < N; ++c)
            result[r][c] = matrix[r][c];
    return result;
}

// Break an N x N matrix into gridSize rows of gridSize blocks,
// each block is blockSize x blockSize
Grid makeBlocks(const vector<vector<int>>& matrix,
    int gridSize,
    int blockSize)
{
    int N = matrix.size();
    Grid blocks(gridSize,
        vector<Block>(gridSize,
            Block(blockSize, vector<int>(blockSize, 0))));
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < N; ++c) {
            int blockRow = r / blockSize;
            int blockCol = c / blockSize;
            int inBlockRow = r % blockSize;
            int inBlockCol = c % blockSize;
            blocks[blockRow][blockCol][inBlockRow][inBlockCol]
                = matrix[r][c];
        }
    }
    return blocks;
}

// Reassemble gridSize rows of gridSize blocks,
// each blockSize x blockSize, into one big matrix
vector<vector<int>> assemble(const Grid& blocks) {
    int gridSize = blocks.size();
    int blockSize = blocks[0][0].size();
    int N = gridSize * blockSize;
    vector<vector<int>> matrix(N, vector<int>(N, 0));
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < N; ++c) {
            int blockRow = r / blockSize;
            int blockCol = c / blockSize;
            int inBlockRow = r % blockSize;
            int inBlockCol = c % blockSize;
            matrix[r][c]
                = blocks[blockRow][blockCol][inBlockRow][inBlockCol];
        }
    }
    return matrix;
}

// Rotate each row of the block grid left by the amounts in rowShifts
void shiftBlockRows(Grid& blocks,
    const vector<int>& rowShifts) {
    int gridSize = blocks.size();
    for (int r = 0; r < gridSize; ++r) {
        int shift = rowShifts[r] % gridSize;
        if (shift == 0) continue;
        rotate(blocks[r].begin(),
            blocks[r].begin() + shift,
            blocks[r].end());
    }
}

// Rotate each column of the block grid up by the amounts in colShifts
void shiftBlockCols(Grid& blocks,
    const vector<int>& colShifts) {
    int gridSize = blocks.size();
    for (int c = 0; c < gridSize; ++c) {
        int shift = colShifts[c] % gridSize;
        if (shift == 0) continue;
        vector<Block> column(gridSize);
        for (int r = 0; r < gridSize; ++r)
            column[r] = blocks[r][c];
        rotate(column.begin(),
            column.begin() + shift,
            column.end());
        for (int r = 0; r < gridSize; ++r)
            blocks[r][c] = column[r];
    }
}

// Multiply two blockSize x blockSize blocks A and B into C
void multiplyAcc(const Block& A,
    const Block& B,
    Block& C)
{
    int blockSize = A.size();
    for (int i = 0; i < blockSize; ++i)
        for (int k = 0; k < blockSize; ++k)
            for (int j = 0; j < blockSize; ++j)
                C[i][j] += A[i][k] * B[k][j];
}

// Cannon multiplication emulation: computes A x B = C
// using processCount virtual processes, padding as needed
void cannonMultiply(const vector<vector<int>>& matrixA,
    const vector<vector<int>>& matrixB,
    vector<vector<int>>& matrixC,
    int                        processCount)
{
    int matrixSize = matrixA.size();

    // 1) Determine gridSize: if processCount is a perfect square,
    //    use exact sqrt; otherwise ceil(sqrt)
    double sqrtP = sqrt(double(processCount));
    int    nearestRoot = int(floor(sqrtP + 0.5));
    bool   isSquare = (nearestRoot * nearestRoot == processCount);
    int    gridSize = isSquare ? nearestRoot
        : int(ceil(sqrtP));

    // 2) Determine blockSize so gridSize*blockSize >= matrixSize
    bool dividesEvenly = (matrixSize % gridSize == 0);
    int  blockSize = dividesEvenly
        ? matrixSize / gridSize
        : int(ceil(double(matrixSize) / gridSize));
    int  paddedSize = gridSize * blockSize;

    // 3) Pad A and B if needed
    vector<vector<int>> paddedA = (isSquare && dividesEvenly)
        ? matrixA
        : padMatrix(matrixA, paddedSize);
    vector<vector<int>> paddedB = (isSquare && dividesEvenly)
        ? matrixB
        : padMatrix(matrixB, paddedSize);

    // 4) Partition into blocks
    Grid blockGridA = makeBlocks(paddedA, gridSize, blockSize);
    Grid blockGridB = makeBlocks(paddedB, gridSize, blockSize);

    // 5) Initial skew: row i left by i, column j up by j
    vector<int> rowShifts(gridSize), colShifts(gridSize);
    for (int i = 0; i < gridSize; ++i) {
        rowShifts[i] = i;
        colShifts[i] = i;
    }
    shiftBlockRows(blockGridA, rowShifts);
    shiftBlockCols(blockGridB, colShifts);

    // 6) Allocate zeroed C blocks
    Grid blockGridC(gridSize,
        vector<Block>(gridSize,
            Block(blockSize, vector<int>(blockSize, 0))));

    // 7) gridSize steps of multiply + rotate
    for (int step = 0; step < gridSize; ++step) {
        // local multiply-accumulate
        for (int r = 0; r < gridSize; ++r) {
            for (int c = 0; c < gridSize; ++c) {
                multiplyAcc(blockGridA[r][c],
                    blockGridB[r][c],
                    blockGridC[r][c]);
            }
        }
        // rotate each row/column by 1 for next step
        fill(rowShifts.begin(), rowShifts.end(), 1);
        fill(colShifts.begin(), colShifts.end(), 1);
        shiftBlockRows(blockGridA, rowShifts);
        shiftBlockCols(blockGridB, colShifts);
    }

    // 8) Reassemble and trim to original size
    vector<vector<int>> paddedC = assemble(blockGridC);
    for (int r = 0; r < matrixSize; ++r) {
        for (int c = 0; c < matrixSize; ++c) {
            matrixC[r][c] = paddedC[r][c];
        }
    }
}

int main() {
    int matrixSize;
    cout << "Matrix dimension n: ";
    cin >> matrixSize;
    while (matrixSize <= 1) {
        cout << "Size must be > 1. Enter n again: ";
        cin >> matrixSize;
    }

    cout << "Randomize matrices? (y/n): ";
    char randomizeChoice;
    cin >> randomizeChoice;

    vector<vector<int>> matrixA(matrixSize, vector<int>(matrixSize)),
        matrixB(matrixSize, vector<int>(matrixSize)),
        matrixC(matrixSize, vector<int>(matrixSize, 0));

    if (randomizeChoice == 'y' || randomizeChoice == 'Y') {
        srand(time(0));
        for (int r = 0; r < matrixSize; ++r)
            for (int c = 0; c < matrixSize; ++c) {
                matrixA[r][c] = rand() % 20;
                matrixB[r][c] = rand() % 20;
            }
        cout << "\nMatrix A:\n"; printMatrix(matrixA);
        cout << "Matrix B:\n"; printMatrix(matrixB);
    }
    else {
        cout << "\nEnter A (" << matrixSize << " x " << matrixSize << "):\n";
        for (int r = 0; r < matrixSize; ++r)
            for (int c = 0; c < matrixSize; ++c)
                cin >> matrixA[r][c];
        cout << "\nEnter B (" << matrixSize << " x " << matrixSize << "):\n";
        for (int r = 0; r < matrixSize; ++r)
            for (int c = 0; c < matrixSize; ++c)
                cin >> matrixB[r][c];
    }

    int processCount;
    cout << "Simulate how many processes? ";
    cin >> processCount;

    cout << "\nStarting Cannon emulation with "
        << processCount << " processes.\n\n";

    double sqrtP = sqrt(double(processCount));
    int    nearestRoot = int(floor(sqrtP + 0.5));
    if (nearestRoot * nearestRoot == processCount && matrixSize % nearestRoot == 0) {
        cout << "Using grid " << nearestRoot << " x " << nearestRoot
            << " with block size " << (matrixSize / nearestRoot) << ".\n\n";
    }
    else {
        int nearestRoot = int(ceil(sqrtP));
        int blockSize = (matrixSize % nearestRoot == 0) ? (matrixSize / nearestRoot)
            : int(ceil(double(matrixSize) / nearestRoot));
        cout << "Using grid " << nearestRoot << " x " << nearestRoot
            << ", padded block size " << blockSize << ".\n\n";
    }

    cannonMultiply(matrixA, matrixB, matrixC, processCount);

    cout << "Result C = A x B:\n";
    printMatrix(matrixC);

    return 0;
}
