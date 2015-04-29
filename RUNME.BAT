
echo "RUNNING MPI_TEST"
mpirun -n 12 mpi_test > mpi_test.out

echo "RUNNING FT_TEST"
./ft_test > ft_test.out

echo "MAKING README"

cp README.H README

diff ft_test.out mpi_test.out >> README

ls -tlr README

echo "OK DONE"


