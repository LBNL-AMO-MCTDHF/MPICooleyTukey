


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


subroutine myzfft1d_slowindex_local(in,out,dim1,dim2,howmany)
  implicit none
  integer, intent(in) :: dim1,dim2,howmany
  complex*16, intent(in) :: in(dim1,dim2,howmany)
  complex*16, intent(out) :: out(dim1,dim2,howmany)
  complex*16 :: intrans(dim2,dim1,howmany),outtrans(dim2,dim1,howmany)
  integer :: ii
  do ii=1,howmany
     intrans(:,:,ii)=TRANSPOSE(in(:,:,ii))
  enddo
  call myzfft1d(intrans,outtrans,dim2,dim1*howmany)
  do ii=1,howmany
     out(:,:,ii)=TRANSPOSE(outtrans(:,:,ii))
  enddo
end subroutine myzfft1d_slowindex_local


#else


subroutine myzfft1d(in,out,dim,howmany)
  use, intrinsic :: iso_c_binding
  implicit none
  include "fftw3.f03"
  integer, intent(in) :: dim,howmany
  complex*16,intent(in) :: in(dim,howmany)
  complex*16, intent(out) :: out(dim,howmany)
  call myzfft1d0(1,in,out,dim,howmany)
end subroutine myzfft1d


subroutine myzfft1d_slowindex_local(in,out,dim1,dim2,howmany)
  implicit none
  integer, intent(in) :: dim1,dim2,howmany
  complex*16, intent(in) :: in(dim1,dim2,howmany)
  complex*16, intent(out) :: out(dim1,dim2,howmany)
  call myzfft1d0(dim1,in,out,dim2,howmany)
end subroutine myzfft1d_slowindex_local


recursive subroutine myzfft1d0(blockdim,in,out,dim,howmany)
  use, intrinsic :: iso_c_binding
  implicit none
  include "fftw3.f03"
  integer, intent(in) :: dim,howmany,blockdim
  complex*16 :: in(blockdim,dim,howmany)    !! cannot be declared intent(in)...hmmm...
  complex*16, intent(out) :: out(blockdim,dim,howmany)
  integer, parameter :: maxplans=3
  type(C_PTR),save :: plans(maxplans)
  integer, save :: plandims(maxplans)=-999, planhowmany(maxplans)=-999,&
       planblockdim(maxplans)=-999
  integer,save :: icalleds(maxplans)=0, numplans=0
  integer :: ostride,istride,onembed(1),inembed(1),idist,odist, dims(1),iplan,thisplan
  integer :: myrank,nprocs
  call getmyranknprocs(myrank,nprocs)

!!$  KEEPME           EITHER WORK BLOCKDIM=1            KEEPME
!!$
!!$  inembed(1)=dim; onembed(1)=dim; idist=dim; odist=dim; istride=1; ostride=1; dims(1)=dim
!!$  inembed(1)=dim; onembed(1)=dim; idist=1; odist=1; istride=1; ostride=1; dims(1)=dim
!!$

  inembed(1)=dim; onembed(1)=dim; idist=1; odist=1; istride=blockdim; ostride=blockdim; dims(1)=dim

  if (numplans.eq.0) then
     numplans=1
     thisplan=1
     plandims(thisplan)=dim; planhowmany(thisplan)=howmany;
     planblockdim(thisplan)=blockdim
  else
     thisplan= -99
     do iplan=1,numplans
        if (plandims(iplan).eq.dim.and.planhowmany(iplan).eq.howmany&
             .and.planblockdim(iplan).eq.blockdim) then
           if (icalleds(iplan).eq.0) then
              print *, "ERROR, plan not done ",iplan,dim,howmany; call mpistop()
           endif
           thisplan=iplan
           exit
        endif
     enddo
     if (thisplan.eq.-99) then
        if (numplans.eq.maxplans) then
           print *,  "all plans taken!", maxplans; call mpistop()
        endif
        numplans=numplans+1
        thisplan=numplans
        plandims(thisplan)=dim; planhowmany(thisplan)=howmany;
        planblockdim(thisplan)=blockdim
     endif
  endif
  if (icalleds(thisplan).eq.0) then
     if (myrank.eq.1) then
        print *, "       Making a 1D FFT plan ", thisplan,  howmany, blockdim
        print *, "       ", dims
     endif
     plans(thisplan) = fftw_plan_many_dft(1,dims,howmany*blockdim,in,inembed,istride,idist,out,onembed,ostride,odist,FFTW_FORWARD,FFTW_EXHAUSTIVE) 
  endif
  icalleds(thisplan)=1    

  call fftw_execute_dft(plans(thisplan), in,out)

end subroutine myzfft1d0

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


