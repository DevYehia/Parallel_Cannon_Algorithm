#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>      // system_clock

using namespace std;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int P, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // 1) Must have P = q*q
    int q = (int)std::sqrt(P);
    if (q * q != P)
    {
        if (rank == 0)
            std::cerr << "Error: Number of processes must be a perfect square.\n";
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // 2) Build a 2D Cartesian communicator, periodic in both dims
    MPI_Comm comm2d;
    int dims[2] = {q, q};
    int periods[2] = {1, 1}; // wraparound
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm2d);

    // Get my coords in the grid
    int coords[2];
    MPI_Cart_coords(comm2d, rank, 2, coords);
    int myRow = coords[0], myCol = coords[1];

    // 3) Root reads n, broadcasts to all
    int n;
    if (rank == 0)
    {
        std::cout << "Enter matrix dimension n: ";
        std::cin >> n;
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // 4) Compute blockSize and padded size
    int blockSize = (n + q - 1) / q; // = ceil(n / q)
    int nPadded = q * blockSize;     // padded dimension

    // 5) Allocate and (on root) read + pad A and B
    std::vector<int> Aflat, Bflat;
    if (rank == 0)
    {
        Aflat.assign(nPadded * nPadded, 0);
        Bflat.assign(nPadded * nPadded, 0);

        char randChoice;
        std::cout << "Randomize matrices(y/n)\n";
        std::cin >> randChoice;

        if (randChoice == 'y' || randChoice == 'Y')
        {
            srand(time(0));
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    Aflat[i * nPadded + j] = rand() % 20;
                }
            }
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {

                    Bflat[i * nPadded + j] = rand() % 20;
                }
            }
        }
        else
        {
            std::cout << "Enter matrix A (" << n << "x" << n << "):\n";
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    int x;
                    std::cin >> x;
                    Aflat[i * nPadded + j] = x;
                }
            }
            std::cout << "Enter matrix B (" << n << "x" << n << "):\n";
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    int x;
                    std::cin >> x;
                    Bflat[i * nPadded + j] = x;
                }
            }
        }
    }

    // 6) Allocate local blocks and result block
    std::vector<int> Ablock(blockSize * blockSize),
        Bblock(blockSize * blockSize),
        Cblock(blockSize * blockSize, 0);

    // 7) Create MPI datatype for a blockSize x blockSize submatrix
    MPI_Datatype blockType;
    MPI_Type_vector(blockSize, blockSize, nPadded, MPI_INT, &blockType);
    MPI_Type_create_resized(blockType, 0, sizeof(int), &blockType);
    MPI_Type_commit(&blockType);

    // 8) Compute displacements for Scatterv/Gatherv
    std::vector<int> displs(P), counts(P, 1);
    if (rank == 0)
    {
        for (int i = 0; i < q; ++i)
        {
            for (int j = 0; j < q; ++j)
            {
                displs[i * q + j] = i * nPadded * blockSize + j * blockSize;
            }
        }
    }

    auto start = chrono::high_resolution_clock::now();
    // 9) Scatter the blocks of A and B
    MPI_Scatterv(
        Aflat.data(), counts.data(), displs.data(), blockType,
        Ablock.data(), blockSize * blockSize, MPI_INT,
        0, comm2d);
    MPI_Scatterv(
        Bflat.data(), counts.data(), displs.data(), blockType,
        Bblock.data(), blockSize * blockSize, MPI_INT,
        0, comm2d);


    // 10) Initial alignment ("skew")
    MPI_Status status;
    int src, dst;
    // 10a) Shift A left by myRow steps
    MPI_Cart_shift(comm2d, 1, -1, &src, &dst);
    for (int i = 0; i < myRow; ++i)
    {
        MPI_Sendrecv_replace(
            Ablock.data(), blockSize * blockSize, MPI_INT,
            dst, 0, src, 0, comm2d, &status);
    }
    // 10b) Shift B up by myCol steps
    MPI_Cart_shift(comm2d, 0, -1, &src, &dst);
    for (int i = 0; i < myCol; ++i)
    {
        MPI_Sendrecv_replace(
            Bblock.data(), blockSize * blockSize, MPI_INT,
            dst, 0, src, 0, comm2d, &status);
    }

    // 11) The main Cannon loop
    for (int step = 0; step < q; ++step)
    {
        // 11a) Local multiply-accumulate
        for (int i = 0; i < blockSize; ++i)
        {
            for (int k = 0; k < blockSize; ++k)
            {
                int a = Ablock[i * blockSize + k];
                for (int j = 0; j < blockSize; ++j)
                {
                    Cblock[i * blockSize + j] +=
                        a * Bblock[k * blockSize + j];
                }
            }
        }
        // 11b) Shift A one step left
        MPI_Cart_shift(comm2d, 1, -1, &src, &dst);
        MPI_Sendrecv_replace(
            Ablock.data(), blockSize * blockSize, MPI_INT,
            dst, 0, src, 0, comm2d, &status);
        // 11c) Shift B one step up
        MPI_Cart_shift(comm2d, 0, -1, &src, &dst);
        MPI_Sendrecv_replace(
            Bblock.data(), blockSize * blockSize, MPI_INT,
            dst, 0, src, 0, comm2d, &status);
    }


    // 12) Gather Cblocks back to root into paddedCflat
    std::vector<int> paddedCflat;
    if (rank == 0)
    {
        paddedCflat.assign(nPadded * nPadded, 0);
    }
    MPI_Gatherv(
        Cblock.data(), blockSize * blockSize, MPI_INT,
        paddedCflat.data(), counts.data(), displs.data(), blockType,
        0, comm2d);

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    
    if(rank == 0){
        std::cout << "Duration is " << duration.count() << " milliseconds\n";
    }
    // 13) Root prints only the top-left n x n of paddedCflat
    // if (rank == 0)
    // {
    //     std::cout << "Result C = A x B:\n";
    //     for (int i = 0; i < n; ++i)
    //     {
    //         for (int j = 0; j < n; ++j)
    //         {
    //             std::cout << paddedCflat[i * nPadded + j] << ' ';
    //         }
    //         std::cout << '\n';
    //     }
    // }

    MPI_Type_free(&blockType);
    MPI_Comm_free(&comm2d);
    MPI_Finalize();
    return 0;
}