
 ## 
 ## the object is as follows:
 ## 
 ## have OUT-OF-PLACE MPI FFT, and OUT-OF-PLACE INVERSE MPI FFT
 ## 
 ## OUT-OF-PLACE MPI FFT performs FFT but leaves the FT in a different order
 ##   depending upon number of processors.  This is Cooley-Tukey without
 ##   reordering.  We do not need to distribute the FT the same way as we
 ##   are distributing the input in position space (as contiguous blocks)
 ##   so we will not do that.
 ## 
 ## OUT-OF-PLACE INVERSE MPI FFT performs the inverse of this operation
 ##
 ## RUNME.BAT creates README  (This header text is in README.H)
 ##   diff of ft_test.out and mpi_test.out follows
 ##
 ## IT WORKS (v8.00)
 ##
 ##			  
