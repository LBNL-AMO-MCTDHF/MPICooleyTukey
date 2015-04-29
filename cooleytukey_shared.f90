

subroutine getprimefactor(dim,myfactor)
  implicit none
  integer, intent(in) :: dim
  integer, intent(out) :: myfactor
  integer :: iprime
  integer, parameter :: numprimes=31
  integer, parameter :: primelist(numprimes)=&
       (/  2,  3,  5,  7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,&
          73, 79, 83, 89, 97,101,103,107,109,113,127 /)  ! no need to go remotely this high

  myfactor=dim

  do iprime=1,numprimes
     if (mod(dim,primelist(iprime)).eq.0) then
        myfactor=primelist(iprime)
        return
     endif
  enddo

end subroutine getprimefactor

subroutine getallprimefactors(dim,numfactors,allfactors)
  implicit none
  integer, intent(in) :: dim
  integer, intent(out) :: allfactors(7),numfactors
  integer :: thisdim,flag
  allfactors(:)=1
  numfactors=1
  thisdim=dim
  flag=0

798 if (numfactors.eq.7) then
     allfactors(7)=thisdim
     return
  else
     call getprimefactor(thisdim,allfactors(numfactors))
     if (allfactors(numfactors).eq.thisdim) then
        return
     endif
     thisdim=thisdim/allfactors(numfactors)
     numfactors=numfactors+1
  endif
go to 798
end subroutine getallprimefactors
  


!! haven't done proper version with recursion
!
!subroutine cooleytukey_replace(blocksize,intranspose,out,dim1,dim2,howmany)
!  implicit none
!  integer, intent(in) :: dim1,dim2,howmany,blocksize
!  complex*16, intent(in) :: intranspose(blocksize,dim1,dim2,howmany)
!  complex*16, intent(out) :: out(blocksize,dim1,dim2,howmany)
!  integer :: ii,jj
!
!  do ii=1,howmany
!     do jj=1,blocksize
!        out(jj,:,:,ii)=RESHAPE(TRANSPOSE(intranspose(jj,:,:,ii)),(/dim1,dim2/))
!     enddo
!  enddo
!
!end subroutine cooleytukey_replace



  

subroutine gettwiddlesmall(twiddlefacs,dim1,dim2)
  implicit none
  integer, intent(in) :: dim1,dim2
  complex*16, intent(out) :: twiddlefacs(dim1)
  complex*16 :: phi
  integer :: k1, itwiddle(dim1)
  real*8, parameter :: pi=3.14159265358979323846264338327950d0
  phi=exp((0d0,-2d0) * pi / (dim1*dim2))
  do k1=1,dim1
     itwiddle(k1)=(k1-1)
  enddo
  twiddlefacs(:)=phi**itwiddle(:)
end subroutine gettwiddlesmall



