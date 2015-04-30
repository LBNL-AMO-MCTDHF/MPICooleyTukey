
nprocs=60

## nprocs=12

mpiargs="-npernode 16 --bysocket"

echo "RUNNING MPI_TEST"
mpirun $mpiargs -n $nprocs mpi_test > mpi_test.out

echo "RUNNING FT_TEST"
./ft_test > ft_test.out

echo "MAKING README"

rm -f README

cp README.H README

diff ft_test.out mpi_test.out >> README

chmod -w README

ls -tlr README

echo "OK DONE"


