
nprocs=6

echo "RUNNING MPI_TEST"
mpirun -n $nprocs mpi_test > mpi_test.out

echo "RUNNING FT_TEST"
./ft_test > ft_test.out

echo "MAKING README"

cp README.H README

diff ft_test.out mpi_test.out >> README

chmod -w README

ls -tlr README

echo "OK DONE"


