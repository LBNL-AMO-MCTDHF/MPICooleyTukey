
program mpi_test
  implicit none
  integer :: myrank,nprocs,mpifileptr

  call mpistart()
  call getmyranknprocs(myrank,nprocs)

  if (myrank.eq.1) then
     mpifileptr=6
  else
     mpifileptr=987
     open(mpifileptr,file="/dev/null", status="unknown")
  endif

  call ct_init(1, mpifileptr)

  call mpi_core(myrank,nprocs,mpifileptr)
end program mpi_test



subroutine mpi_core(myrank,nprocs,mpifileptr)
  implicit none
  integer, intent(in) :: myrank,nprocs,mpifileptr
!!$  integer, parameter :: size1=70000
  integer, parameter :: size1=7
!!$  integer, parameter :: size1=35
!!$  integer, parameter :: size1=13
!!$  integer, parameter :: size1=26
!!$  integer, parameter :: size1=52
  integer :: size
  complex*16 :: input(size1*nprocs), zoutput1(size1,nprocs),&
       input1(size1,nprocs), & !!TEMP output(size1*nprocs),
       input0(size1*nprocs)  !!,output1(size1,nprocs)
  real*8 :: realarray(size1*nprocs),randomamount
  integer :: i

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

  write(mpifileptr,*) "randomamount is ", randomamount

!!$  write(mpifileptr,*) "## call ft.  "
!!$  call myzfft1d(input,output,size,1)

  input1(:,:)=RESHAPE(input,(/size1,nprocs/))

  if (myrank.eq.1) then
     write(mpifileptr,'(A50,$)') "## call cooleytukey_outofplace_mpi  "
     call system('date')
  else
     write(mpifileptr,*) "## call cooleytukey_outofplace_mpi  "
  endif
  call ctdim(1)
  call cooleytukey_outofplace_forward_mpi(input1(:,myrank),zoutput1(:,myrank),1,1,size1,1)

  call mympigather(zoutput1(:,myrank),zoutput1,size1)

  write(mpifileptr,*) "## call cooleytukey_outofplace_backward_mpi"
  call ctdim(1)
  call cooleytukey_outofplace_backward_mpi(zoutput1(:,myrank),input1(:,myrank),1,1,size1,1)
  input1(:,myrank)=input1(:,myrank)/size1/nprocs

  if (myrank.eq.1) then
     write(mpifileptr,'(A50,$)') "## done cooleytukey_outofplace_backward_mpi  "
     call system('date')
  else
     write(mpifileptr,*) "## done cooleytukey_outofplace_backward_mpi  "
  endif

  call mympigather(input1(:,myrank),input0,size1)

!!$  call cooleytukey_replace(1,zoutput1,output1,size/nprocs,nprocs,1)


  if (myrank.eq.1) then
     print *, "## ok"
     print *
     write(*,'(A5,100A10)') " ","input","input" !!,"output","output"
     print *
     
     do i=1,size
        write(*,'(I20,100F10.5)') i, &
             abs(input(i)), &
             abs(input0(i))
!!$             abs(output(i)), &
!!$             abs(output1(i))
        write(*,'(I20,100F10.5)') i, &
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


subroutine myzfft3d(in,out,dim1,dim2,dim3,howmany)
  implicit none
  integer, intent(in) :: dim1,dim2,dim3,howmany
  complex*16, intent(in) :: in(dim1,dim2,dim3,howmany)
  complex*16, intent(out) :: out(dim1,dim2,dim3,howmany)
  print *, "myzfft3d not programmed mpi_test"; call mpistop()
  out(:,:,:,:)=in(:,:,:,:)  !! avoid warn unused
end subroutine myzfft3d
  
