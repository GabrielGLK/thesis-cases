#!/bin/sh

echo "start!"

cd distance/0.55
CC='mpicc -D_MPI=12' make rising.tst

cd ../0.6
CC='mpicc -D_MPI=12' make rising.tst

cd ../0.65
CC='mpicc -D_MPI=12' make rising.tst

cd ../0.7
CC='mpicc -D_MPI=12' make rising.tst

cd ../0.75
CC='mpicc -D_MPI=12' make rising.tst

cd ../0.8
CC='mpicc -D_MPI=12' make rising.tst

cd ../../size/0.55
CC='mpicc -D_MPI=12' make rising.tst

cd ../0.6
CC='mpicc -D_MPI=12' make rising.tst

cd ../0.65
CC='mpicc -D_MPI=12' make rising.tst

cd ../0.7
CC='mpicc -D_MPI=12' make rising.tst

cd ../0.75
CC='mpicc -D_MPI=12' make rising.tst

cd ../0.8
CC='mpicc -D_MPI=12' make rising.tst

echo "finished!!!"