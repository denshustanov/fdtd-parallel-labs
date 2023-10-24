#include "heat_scheme.cpp"
#include <chrono>
#include <iostream>
#include <mpi.h>

int main(int argc, char **argv)
{
    int rank, total;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &total);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int I = 9;
    int J = 99;
    const int iterations = 1;
    ModelParams params = prepareParams(I, J);
    int vecSize = (I + 1) * (J + 1);
    double *m = (double *)calloc(vecSize, sizeof(double));

    /*
        Tasks preparation
    */
    int blocks[total];
    int starts[total];
    int ends[total];

    int start = 0;
    int end = 0;
    int block_size = 0;
    if (rank == 0)
    {
        for (int i = 0; i < total; i++)
        {
            blocks[i] = (I + 1) / total + (i < ((I + 1) % total) ? 1 : 0);
            end += blocks[i];
            starts[i] = start;
            ends[i] = end;

            std::cout << calcIndexCol(start, 0, J) << " " << calcIndexCol(end - 1, J, J) << std::endl;
            start += blocks[i];
        }
    }

    MPI_Scatter(blocks, 1, MPI_INT, &block_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(starts, 1, MPI_INT, &start, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(ends, 1, MPI_INT, &end, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // initial condition
    for (int i = start; i < end; i++)
    {
        m[calcIndexCol(i, 0, J)] = initCondition(i, I);
    }

    double bottom;
    double bottomRight;
    double bottomLeft;

    MPI_Status status;

    // iterating on time layers
    for (int j = 1; j <= J; j++)
    {

        /*
            Step one - communication
            exchanging block boundary values from previous layer with adjacent processes
        */

        // processes that dont ever calculate right boundary
        if (rank < total - 1)
        {

            // std::cout << "rank: " << rank << " sending " << end - 1 << " receiving " << end << std::endl;
            // sending rightmost value
            MPI_Send(
                m + calcIndexCol(end - 1, j - 1, J),
                1,
                MPI_INT,
                rank + 1,
                NULL,
                MPI_COMM_WORLD);
            // receiving value to the right of rightmost
            MPI_Recv(
                m + calcIndexCol(end, j + 1, J),
                1,
                MPI_INT,
                rank + 1,
                NULL,
                MPI_COMM_WORLD,
                &status);

            // left boundary appears here
        }
        // processes that dont ever calculate left boundary
        if (rank > 0)
        {
            // sending leftmost value
            // std::cout << "rank: " << rank << " sending " << start << " receiving " << start - 1 << std::endl;
            MPI_Send(
                m + calcIndexCol(start, j - 1, J),
                1,
                MPI_INT,
                rank - 1,
                NULL,
                MPI_COMM_WORLD);
            // receiving value to the left of leftmost
            MPI_Recv(
                m + calcIndexCol(start - 1, j - 1, J),
                1,
                MPI_INT,
                rank - 1,
                NULL,
                MPI_COMM_WORLD,
                &status);
        }

        /*
            Step two - calculations
        */
        if (rank == 0)
        {
            // left right boundary
            bottom = m[calcIndexCol(0, j - 1, J)];
            bottomLeft = m[calcIndexCol(1, j - 1, J)];
            m[calcIndexCol(0, j, J)] = leftCondition(bottom, bottomLeft, params);
            for (int i = start + 1; i < end; ++i)
            {
                bottom = m[calcIndexCol(i, j - 1, J)];
                bottomRight = m[calcIndexCol(i + 1, j - 1, J)];
                bottomLeft = m[calcIndexCol(i - 1, j - 1, J)];
                m[calcIndexCol(i, j, I)] = innerNode(bottom, bottomLeft, bottomRight, params);
            }
        }
        else if (rank == total - 1)
        {
            // right boundary
            m[calcIndexCol(I, j, J)] = 0;
            // other nodes
            for (int i = start; i < end - 1; ++i)
            {
                bottom = m[calcIndexCol(i, j - 1, J)];
                bottomRight = m[calcIndexCol(i + 1, j - 1, J)];
                bottomLeft = m[calcIndexCol(i - 1, j - 1, J)];
                m[calcIndexCol(i, j, I)] = innerNode(bottom, bottomLeft, bottomRight, params);
            }
        }
        else
        {
            // other nodes
            for (int i = start; i < end; ++i)
            {
                if (rank == 1)
                {
                    std::cout << "rank " << rank << " " << i << " " << j << '\n';
                }
                bottom = m[calcIndexCol(i, j - 1, J)];
                bottomRight = m[calcIndexCol(i + 1, j - 1, J)];
                bottomLeft = m[calcIndexCol(i - 1, j - 1, J)];
                m[calcIndexCol(i, j, I)] = innerNode(bottom, bottomLeft, bottomRight, params);
            }
        }
    }

    /*
        Sending results to root process
    */

    if (rank == 0)
    {
        for (int i = 1; i < total; i++)
        {
            MPI_Recv(
                m + calcIndexCol(starts[i], 0, J),
                (J + 1) * blocks[i],
                MPI_DOUBLE,
                i,
                NULL,
                MPI_COMM_WORLD,
                &status);
        }
    }
    else
    {
        MPI_Send(
            m + calcIndexCol(start, 0, J),
            (J + 1) * block_size,
            MPI_DOUBLE,
            0,
            NULL,
            MPI_COMM_WORLD);
    }

    std::string fname = "output_mpi_" + std::to_string(rank) + ".txt";

    saveMatrix(m, I, J, fname.c_str());

    std::cout << "p " << rank << " of " << total << ", block size: " << block_size << ", start: " << start << ", end: " << end - 1 << std::endl;

    MPI_Finalize();
    return 0;
}