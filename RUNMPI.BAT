
 nprocs=384

## nprocs=48

## nprocs=9

## nprocs=480

## nprocs=60

## nprocs=12

mpiargs="-npernode 16 --bysocket"

echo "RUNNING MPI_TEST"
mpirun $mpiargs -n $nprocs mpi_test

