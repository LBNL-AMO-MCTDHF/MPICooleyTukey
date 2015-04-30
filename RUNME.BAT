
nprocs=60

## nprocs=12

mpiargs="-npernode 16 --bysocket"

echo -n "RUNNING MPI_TEST   "
date
mpirun $mpiargs -n $nprocs mpi_test > mpi_test.out
echo  -n "DONE MPI_TEST   "
date

echo  -n "RUNNING FT_TEST   "
date
./ft_test > ft_test.out
echo  -n "DONE FT_TEST   "
date

echo "MAKING README"

rm -f README

cp README.H README

diff ft_test.out mpi_test.out >> README

# too big
rm ft_test.out mpi_test.out  

chmod -w README

ls -tlr README

echo "OK DONE"


