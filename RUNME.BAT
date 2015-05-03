
nprocs=12

## nprocs=384

## nprocs=48

## nprocs=9

## nprocs=480

## nprocs=60

## nprocs=12

mpiargs="-npernode 16 --bysocket"

echo "RUNNING MPI_TEST"
mpirun $mpiargs -n $nprocs mpi_test > SCRATCH/mpi_test.out
echo "DONE MPI_TEST"

## broken
## echo "RUNNING FT_TEST (broken)"
## mpirun -n 1 ./ft_test > SCRATCH/ft_test.out
## echo "DONE FT_TEST"

echo "MAKING README"

rm -f README

cp README.H README

### broken  diff SCRATCH/ft_test.out SCRATCH/mpi_test.out |grep -v "MPI STOP" >> README

echo "ft_test is broken" >> README
echo "here is top of mpi_test.out" >> README
echo >> README
head -n 150 SCRATCH/mpi_test.out >> README


# too big
# but on scratch, not in repo 
# rm ft_test.out mpi_test.out  

chmod -w README

ls -tlr README

echo "OK DONE"


