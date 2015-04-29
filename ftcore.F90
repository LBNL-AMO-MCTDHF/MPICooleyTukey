


!! CORE SUBROUTINES

#ifndef FFTWFLAG

subroutine myzfft1d(in,out,dim,howmany)
  implicit none
  integer, intent(in) :: dim,howmany
  integer :: k
  complex*16, intent(in) :: in(dim,howmany)
  complex*16, intent(out) :: out(dim,howmany)
  complex*16 :: wsave(4*dim+15,howmany)   ! MAKE BIGGER IF SEGFAULT... iffy
  out(:,:)=in(:,:)
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(in,out,dim,howmany,wsave)
!$OMP DO SCHEDULE(STATIC)
  do k=1,howmany
     call zffti(dim,wsave(:,k))
     call zfftf(dim,out(:,k),wsave(:,k))
  enddo
!$OMP END DO
!$OMP END PARALLEL
end subroutine myzfft1d


#else


subroutine myzfft1d(in,out,dim,howmany)
  use, intrinsic :: iso_c_binding
  implicit none
  include "fftw3.f03"
  integer, intent(in) :: dim,howmany
  type(C_PTR),save :: plan,cdata
  complex*16,intent(in) :: in(dim,howmany)
  complex*16, intent(out) :: out(dim,howmany)
  integer :: ostride,istride,onembed(1),inembed(1),idist,odist, dims(1)
!!$  integer,save :: icalled=0
  complex(C_DOUBLE_COMPLEX) :: data(dim,howmany)

  inembed(1)=dim; onembed(1)=dim; idist=dim; odist=dim; istride=1; ostride=1; dims(1)=dim

!!$  if (icalled.eq.0) then

     plan=fftw_plan_many_dft(1,dims,howmany,data,inembed,istride,idist,data,onembed,ostride,odist,FFTW_FORWARD,FFTW_EXHAUSTIVE) 

!!$  endif
!!$  icalled=1

  data(:,:)=in(:,:)
  call fftw_execute_dft(plan, data,data)
  out(:,:)=data(:,:)

  call fftw_destroy_plan(plan)

end subroutine myzfft1d

#endif


subroutine myzfft1d_inverse(in,out,dim,howmany)
  implicit none
  integer, intent(in) :: dim,howmany
  complex*16, intent(in) :: in(dim,howmany)
  complex*16, intent(out) :: out(dim,howmany)
  complex*16 :: &                               !! AUTOMATIC
       inconjg(dim,howmany), &
       outconjg(dim,howmany)

  inconjg(:,:)=CONJG(in(:,:))
  call myzfft1d(inconjg,outconjg,dim,howmany)
  out(:,:)=CONJG(outconjg(:,:))/dim

end subroutine myzfft1d_inverse


