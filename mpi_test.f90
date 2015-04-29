
program mpi_test
  implicit none
  integer :: myrank,nprocs,mpifileptr

  call mpistart()
  call getmyranknprocs(myrank,nprocs)
  call ctset()

  if (myrank.eq.1) then
     mpifileptr=6
  else
     mpifileptr=987
     open(mpifileptr,file="/dev/null", status="unknown")
  endif

  call mpi_core(myrank,nprocs,mpifileptr)
end program mpi_test



subroutine mpi_core(myrank,nprocs,mpifileptr)
  implicit none
  integer, intent(in) :: myrank,nprocs,mpifileptr
  integer, parameter :: size1=13
  integer :: size
  complex*16 :: input(size1*nprocs), output(size1*nprocs),zoutput1(size1,nprocs),&
       input1(size1,nprocs),output1(size1,nprocs),input0(size1*nprocs)
  real*8 :: realarray(size1*nprocs),randomamount
  integer :: primefactors(7),numfactors,i,proclist(nprocs)

  size=nprocs*size1
  do i=1,size
     input(i)=exp(-249*(i-(size+1)*0.5d0)**2/size**2)  * (-1)**i
  enddo

  randomamount=1d-2

  call RANDOM_NUMBER(realarray)
  input(:)=input(:)+realarray(:) * randomamount

  call RANDOM_NUMBER(realarray)
  input(:)=input(:)+realarray(:)**3 * (0d0,4d0) * randomamount

  write(mpifileptr,*) "Go mpi_test. Dimensions are ",size1,nprocs
!!TEMP
!!TEMP  call getallprimefactors(nprocs,numfactors,primefactors)
!!TEMP
  numfactors=1; primefactors(:)=1; primefactors(1)=nprocs

  write(mpifileptr,*) "     Prime factors of ",nprocs," are"
  write(mpifileptr,*) primefactors(1:numfactors)
  write(mpifileptr,*)

  write(mpifileptr,*) "## call ft. randomamount is",randomamount
  call myzfft1d(input,output,size,1)

  input1(:,:)=RESHAPE(input,(/size1,nprocs/))

  do i=1,nprocs
     proclist(i)=i
  enddo

  write(mpifileptr,*) "## call cooleytukey_outofplace_mpi"
  call cooleytukey_outofplace_mpi(input1(:,myrank),zoutput1(:,myrank),size1,primefactors,proclist,nprocs,myrank,1)

  call mympigather(zoutput1(:,myrank),zoutput1,size1)


  write(mpifileptr,*) "## call cooleytukey_outofplace_inverse_mpi"
  call cooleytukey_outofplace_inverse_mpi(zoutput1(:,myrank),input1(:,myrank),size1,primefactors,proclist,nprocs,myrank,1)

  call mympigather(input1(:,myrank),input0,size1)

!!$  call cooleytukey_replace(1,zoutput1,output1,size/nprocs,nprocs,1)


  if (myrank.eq.1) then
     print *, "## ok"
     print *
     write(*,'(A5,100A10)') " ","input","input" !!,"output","output"
     print *
     
     do i=1,size
        write(*,'(I5,100F10.5)') i, &
             abs(input(i)), &
             abs(input0(i))
!!$             abs(output(i)), &
!!$             abs(output1(i))
        write(*,'(I5,100F10.5)') i, &
             real(input(i)), &
             real(input0(i))
!!$             real(output(i)), &
!!$             real(output1(i))
        print *
     enddo
  endif

  call mpibarrier()
  write(mpifileptr,*) "OK MPI_TEST DONE!"
  call mpistop()

end subroutine mpi_core

