
nprocs=60

## nprocs=12

mpiargs="-npernode 16 --bysocket"

echo "RUNNING MPI_TEST"
mpirun $mpiargs -n $nprocs mpi_test

