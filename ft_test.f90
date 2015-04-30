
program ft_test
  implicit none
  integer, parameter :: size1=70000, size2=480
!!$  integer, parameter :: size1=70000, size2=60
!!$  integer, parameter :: size1=7, size2=60
!!$  integer, parameter :: size1=35, size2=12
!!$  integer, parameter :: size1=13, size2=12
!!$  integer, parameter :: size1=26, size2=6
!!$  integer, parameter :: size1=52, size2=3

  integer, parameter :: size=size1*size2
  complex*16 :: input(size), output(size),zoutput1(size),&
       input1(size) !!,output1(size)
  real*8 :: realarray(size),randomamount
  integer :: primefactors(7),numfactors,i

  do i=1,size
     input(i)=exp(-249*(i-(size+1)*0.5d0)**2/size**2)  * (-1)**i
  enddo

  randomamount=1d-2

  call RANDOM_NUMBER(realarray)
  input(:)=input(:)+realarray(:) * randomamount

  call RANDOM_NUMBER(realarray)
  input(:)=input(:)+realarray(:)**3 * (0d0,4d0) * randomamount

  print *, "Go ft_test. Dimensions are ",size1,size2
  call getallprimefactors(size2,numfactors,primefactors)
  print *, "     Prime factors of ",size2," are"
  print *, primefactors(1:numfactors)
  print *

  write(*,*) "randomamount is ", randomamount
  write(*,'(A50,$)') "## call ft.  "
  call system('date')
  call myzfft1d(input,output,size,1)

  write(*,'(A50,$)') "## done ft   "
  call system('date')

  write(*,'(A50,$)') "## call cooleytukey_outofplace  "
  call system('date')
  call cooleytukey_outofplace(input,zoutput1,size1,primefactors,1)

  print *, "## call cooleytukey_outofplace_inverse"
  call cooleytukey_outofplace_inverse(zoutput1,input1,size1,primefactors,1)
  write(*,'(A50,$)') "## done cooleytukey_outofplace_inverse  "
  call system('date')



!!$  call cooleytukey_replace(1,zoutput1,output1,size1,primefactors,1)

  print *, "## ok"
  print *
  write(*,'(A5,100A10)') " ","input","input" !!,"output","output"
  print *

  do i=1,size
     write(*,'(I20,100F10.5)') i, &
          abs(input(i)), &
          abs(input1(i))
!!$          abs(output(i)), &
!!$          abs(output3(i))
     write(*,'(I20,100F10.5)') i, &
          real(input(i)), &
          real(input1(i))
!!$          real(output(i)), &
!!$          real(output3(i))
     print *
  enddo
  
end program ft_test

