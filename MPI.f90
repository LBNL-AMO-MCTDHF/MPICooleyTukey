
!! PARALLEL SUBROUTINES.


module mpimod
  implicit none
  include "mpif.h"
  integer :: myrank = -1
  integer :: nprocs = -1
  integer :: mpifileptr = -1
  integer :: mpiatime,mpibtime
  integer :: mpitime=0, nonmpitime=0
  integer :: MPI_GROUP_WORLD
end module mpimod

subroutine getmyranknprocs(outmyrank,outnprocs)
  use mpimod
  implicit none
  integer, intent(out) :: outmyrank,outnprocs
  if (nprocs.lt.1.or.myrank.gt.nprocs.or.myrank.lt.1) then
     print *, "ACK, has mpistart been called? myrank,nprocs bad values=",myrank,nprocs; stop
  endif
  outmyrank=myrank; outnprocs=nprocs
end subroutine getmyranknprocs


subroutine getWorldCommGroup(outcommunicator,outgroup)
  use mpimod
  implicit none
  integer, intent(out) :: outcommunicator,outgroup
  outcommunicator=MPI_COMM_WORLD
  outgroup=MPI_GROUP_WORLD
end subroutine getWorldCommGroup


subroutine mpistart()
  use mpimod
  implicit none
  integer :: ierr


  call MPI_INIT(ierr)
  if (ierr/=0) then; print *,  "MPI ERR 1";  stop; 
  endif

  call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr)   
  if (ierr/=0) then;   print *,  "MPI ERR 2";   stop; 
  endif

  myrank=myrank+1

  call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr)     
  if (ierr/=0) then;   print *,   "MPI ERR 3";   stop; 
  endif

  call MPI_COMM_GROUP(MPI_COMM_WORLD,MPI_GROUP_WORLD,ierr)
  if (ierr/=0) then;   print *,   "MPI grouperr";   stop; 
  endif

  if (nprocs.eq.1.or.myrank.eq.1) then
     mpifileptr=6
  else
     mpifileptr=987
     open(mpifileptr,file="/dev/null",status="unknown")
  endif

  call system_clock(mpibtime);  call system_clock(mpiatime)

end subroutine mpistart


subroutine mpibarrier()
  use mpimod
  implicit none
  integer :: ierr
  call system_clock(mpiatime)
  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  if (ierr/=0) then
     write(mpifileptr,*) "MPI ERR 3"; call mpistop()
  endif
  call system_clock(mpibtime)
  mpitime=mpitime+mpibtime-mpiatime
end subroutine mpibarrier



subroutine mpistop()
  use mpimod
  implicit none
  integer :: ierr

  write(mpifileptr,*) "MPI STOP.... ",myrank
  if (mpifileptr.ne.6) then
     close(mpifileptr)
  endif

!  ierr=798
!  call mpi_abort(MPI_COMM_WORLD,798,ierr)
!  if (ierr/=0) then
!     write(*,*)  "MPI ABORT ERR ",ierr,myrank
!  endif

  call mpi_finalize(ierr)
  if (ierr/=0) then
     write(*,*)  "MPI FINALIZE ERR ",ierr,myrank
  endif
!!$  call sleep(2)
!!$  write(*,*) "MPI STOP ",myrank
  stop
end subroutine mpistop



subroutine mympicomplexbcast_local(input, source, isize,MPI_COMM_LOCAL)
  use mpimod
  implicit none
  integer,intent(in) :: isize,source,MPI_COMM_LOCAL
  complex*16, intent(inout) :: input(isize)
  integer :: ierr, isource

  isource=source-1
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_bcast(input,isize,MPI_DOUBLE_COMPLEX,isource,MPI_COMM_LOCAL,ierr)
  if (ierr/=0) then
     write(mpifileptr,*) "ERR mympibcast"; call mpistop()
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympicomplexbcast_local


subroutine mympisendrecv_complex_local(sendbuf, recvbuf, dest, source, tag, isize, MPI_COMM_LOCAL)
  use mpimod
  implicit none
  integer, intent(in) :: dest,source,tag,isize,MPI_COMM_LOCAL
  complex*16, intent(in) :: sendbuf(isize)
  complex*16, intent(out) :: recvbuf(isize)
  integer :: ierr, idest,  isource

  idest=dest-1
  isource=source-1
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_sendrecv(sendbuf,isize,MPI_DOUBLE_COMPLEX,idest,tag,&
       recvbuf,isize, MPI_DOUBLE_COMPLEX,isource,tag,MPI_COMM_LOCAL,MPI_STATUS_IGNORE,ierr)
  if (ierr/=0) then
     write(mpifileptr,*) "ERR mympisendrecv",ierr; call mpistop()
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympisendrecv_complex_local




subroutine mympialltoall_complex(input, output, count)
  use mpimod
  implicit none
  integer :: ierr, count
  complex*16 :: input(count,nprocs), output(count,nprocs)

  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  call mpi_alltoall(input(:,:),count,MPI_DOUBLE_COMPLEX,output(:,:),count,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
  if (ierr.ne.0) then
     write(mpifileptr,*) "ERROR ALLTOALL ", ierr; call mpistop()
  endif

  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime

end subroutine mympialltoall_complex


subroutine mympigather(vectorin,vectorsout,insize)
  use mpimod
  implicit none
  integer :: ierr,insize
  complex*16 :: vectorin(insize),vectorsout(insize,nprocs)

  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  call mpi_allgather(vectorin(:),insize,MPI_DOUBLE_COMPLEX,vectorsout(:,:),insize,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
  if (ierr.ne.0) then
     write(mpifileptr,*) "ORBGATHER ERR ", ierr; call mpistop()
  endif

  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympigather



