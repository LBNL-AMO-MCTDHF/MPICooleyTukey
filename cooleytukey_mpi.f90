

#define MAXFACTORS 7

!!$
!!$Apache License
!!$                           Version 2.0, January 2004
!!$                        http://www.apache.org/licenses/
!!$
!!$   TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION
!!$
!!$   1. Definitions.
!!$
!!$      "License" shall mean the terms and conditions for use, reproduction,
!!$      and distribution as defined by Sections 1 through 9 of this document.
!!$
!!$      "Licensor" shall mean the copyright owner or entity authorized by
!!$      the copyright owner that is granting the License.
!!$
!!$      "Legal Entity" shall mean the union of the acting entity and all
!!$      other entities that control, are controlled by, or are under common
!!$      control with that entity. For the purposes of this definition,
!!$      "control" means (i) the power, direct or indirect, to cause the
!!$      direction or management of such entity, whether by contract or
!!$      otherwise, or (ii) ownership of fifty percent (50%) or more of the
!!$      outstanding shares, or (iii) beneficial ownership of such entity.
!!$
!!$      "You" (or "Your") shall mean an individual or Legal Entity
!!$      exercising permissions granted by this License.
!!$
!!$      "Source" form shall mean the preferred form for making modifications,
!!$      including but not limited to software source code, documentation
!!$      source, and configuration files.
!!$
!!$      "Object" form shall mean any form resulting from mechanical
!!$      transformation or translation of a Source form, including but
!!$      not limited to compiled object code, generated documentation,
!!$      and conversions to other media types.
!!$
!!$      "Work" shall mean the work of authorship, whether in Source or
!!$      Object form, made available under the License, as indicated by a
!!$      copyright notice that is included in or attached to the work
!!$      (an example is provided in the Appendix below).
!!$
!!$      "Derivative Works" shall mean any work, whether in Source or Object
!!$      form, that is based on (or derived from) the Work and for which the
!!$      editorial revisions, annotations, elaborations, or other modifications
!!$      represent, as a whole, an original work of authorship. For the purposes
!!$      of this License, Derivative Works shall not include works that remain
!!$      separable from, or merely link (or bind by name) to the interfaces of,
!!$      the Work and Derivative Works thereof.
!!$
!!$      "Contribution" shall mean any work of authorship, including
!!$      the original version of the Work and any modifications or additions
!!$      to that Work or Derivative Works thereof, that is intentionally
!!$      submitted to Licensor for inclusion in the Work by the copyright owner
!!$      or by an individual or Legal Entity authorized to submit on behalf of
!!$      the copyright owner. For the purposes of this definition, "submitted"
!!$      means any form of electronic, verbal, or written communication sent
!!$      to the Licensor or its representatives, including but not limited to
!!$      communication on electronic mailing lists, source code control systems,
!!$      and issue tracking systems that are managed by, or on behalf of, the
!!$      Licensor for the purpose of discussing and improving the Work, but
!!$      excluding communication that is conspicuously marked or otherwise
!!$      designated in writing by the copyright owner as "Not a Contribution."
!!$
!!$      "Contributor" shall mean Licensor and any individual or Legal Entity
!!$      on behalf of whom a Contribution has been received by Licensor and
!!$      subsequently incorporated within the Work.
!!$
!!$   2. Grant of Copyright License. Subject to the terms and conditions of
!!$      this License, each Contributor hereby grants to You a perpetual,
!!$      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
!!$      copyright license to reproduce, prepare Derivative Works of,
!!$      publicly display, publicly perform, sublicense, and distribute the
!!$      Work and such Derivative Works in Source or Object form.
!!$
!!$   3. Grant of Patent License. Subject to the terms and conditions of
!!$      this License, each Contributor hereby grants to You a perpetual,
!!$      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
!!$      (except as stated in this section) patent license to make, have made,
!!$      use, offer to sell, sell, import, and otherwise transfer the Work,
!!$      where such license applies only to those patent claims licensable
!!$      by such Contributor that are necessarily infringed by their
!!$      Contribution(s) alone or by combination of their Contribution(s)
!!$      with the Work to which such Contribution(s) was submitted. If You
!!$      institute patent litigation against any entity (including a
!!$      cross-claim or counterclaim in a lawsuit) alleging that the Work
!!$      or a Contribution incorporated within the Work constitutes direct
!!$      or contributory patent infringement, then any patent licenses
!!$      granted to You under this License for that Work shall terminate
!!$      as of the date such litigation is filed.
!!$
!!$   4. Redistribution. You may reproduce and distribute copies of the
!!$      Work or Derivative Works thereof in any medium, with or without
!!$      modifications, and in Source or Object form, provided that You
!!$      meet the following conditions:
!!$
!!$      (a) You must give any other recipients of the Work or
!!$          Derivative Works a copy of this License; and
!!$
!!$      (b) You must cause any modified files to carry prominent notices
!!$          stating that You changed the files; and
!!$
!!$      (c) You must retain, in the Source form of any Derivative Works
!!$          that You distribute, all copyright, patent, trademark, and
!!$          attribution notices from the Source form of the Work,
!!$          excluding those notices that do not pertain to any part of
!!$          the Derivative Works; and
!!$
!!$      (d) If the Work includes a "NOTICE" text file as part of its
!!$          distribution, then any Derivative Works that You distribute must
!!$          include a readable copy of the attribution notices contained
!!$          within such NOTICE file, excluding those notices that do not
!!$          pertain to any part of the Derivative Works, in at least one
!!$          of the following places: within a NOTICE text file distributed
!!$          as part of the Derivative Works; within the Source form or
!!$          documentation, if provided along with the Derivative Works; or,
!!$          within a display generated by the Derivative Works, if and
!!$          wherever such third-party notices normally appear. The contents
!!$          of the NOTICE file are for informational purposes only and
!!$          do not modify the License. You may add Your own attribution
!!$          notices within Derivative Works that You distribute, alongside
!!$          or as an addendum to the NOTICE text from the Work, provided
!!$          that such additional attribution notices cannot be construed
!!$          as modifying the License.
!!$
!!$      You may add Your own copyright statement to Your modifications and
!!$      may provide additional or different license terms and conditions
!!$      for use, reproduction, or distribution of Your modifications, or
!!$      for any such Derivative Works as a whole, provided Your use,
!!$      reproduction, and distribution of the Work otherwise complies with
!!$      the conditions stated in this License.
!!$
!!$   5. Submission of Contributions. Unless You explicitly state otherwise,
!!$      any Contribution intentionally submitted for inclusion in the Work
!!$      by You to the Licensor shall be under the terms and conditions of
!!$      this License, without any additional terms or conditions.
!!$      Notwithstanding the above, nothing herein shall supersede or modify
!!$      the terms of any separate license agreement you may have executed
!!$      with Licensor regarding such Contributions.
!!$
!!$   6. Trademarks. This License does not grant permission to use the trade
!!$      names, trademarks, service marks, or product names of the Licensor,
!!$      except as required for reasonable and customary use in describing the
!!$      origin of the Work and reproducing the content of the NOTICE file.
!!$
!!$   7. Disclaimer of Warranty. Unless required by applicable law or
!!$      agreed to in writing, Licensor provides the Work (and each
!!$      Contributor provides its Contributions) on an "AS IS" BASIS,
!!$      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
!!$      implied, including, without limitation, any warranties or conditions
!!$      of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A
!!$      PARTICULAR PURPOSE. You are solely responsible for determining the
!!$      appropriateness of using or redistributing the Work and assume any
!!$      risks associated with Your exercise of permissions under this License.
!!$
!!$   8. Limitation of Liability. In no event and under no legal theory,
!!$      whether in tort (including negligence), contract, or otherwise,
!!$      unless required by applicable law (such as deliberate and grossly
!!$      negligent acts) or agreed to in writing, shall any Contributor be
!!$      liable to You for damages, including any direct, indirect, special,
!!$      incidental, or consequential damages of any character arising as a
!!$      result of this License or out of the use or inability to use the
!!$      Work (including but not limited to damages for loss of goodwill,
!!$      work stoppage, computer failure or malfunction, or any and all
!!$      other commercial damages or losses), even if such Contributor
!!$      has been advised of the possibility of such damages.
!!$
!!$   9. Accepting Warranty or Additional Liability. While redistributing
!!$      the Work or Derivative Works thereof, You may choose to offer,
!!$      and charge a fee for, acceptance of support, warranty, indemnity,
!!$      or other liability obligations and/or rights consistent with this
!!$      License. However, in accepting such obligations, You may act only
!!$      on Your own behalf and on Your sole responsibility, not on behalf
!!$      of any other Contributor, and only if You agree to indemnify,
!!$      defend, and hold each Contributor harmless for any liability
!!$      incurred by, or claims asserted against, such Contributor by reason
!!$      of your accepting any such warranty or additional liability.
!!$
!!$   END OF TERMS AND CONDITIONS
!!$
!!$   Copyright 2015 the regents of the University of California
!!$
!!$   Licensed under the Apache License, Version 2.0 (the "License");
!!$   you may not use this file except in compliance with the License.
!!$   You may obtain a copy of the License at
!!$
!!$       http://www.apache.org/licenses/LICENSE-2.0
!!$
!!$   Unless required by applicable law or agreed to in writing, software
!!$   distributed under the License is distributed on an "AS IS" BASIS,
!!$   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!$   See the License for the specific language governing permissions and
!!$   limitations under the License.


module fileptrmod
  implicit none
  integer :: mpifileptr = 6
end module fileptrmod

module ct_mpimod
  implicit none
  integer :: myrank = -1
  integer :: nprocs = -1
  integer :: CT_COMM_WORLD = -1
  integer :: CT_GROUP_WORLD = -1
end module ct_mpimod


subroutine ctset()
  use fileptrmod
  use ct_mpimod
  implicit none
  call getmyranknprocs(myrank,nprocs)
  call getworldcommgroup(CT_COMM_WORLD,CT_GROUP_WORLD)
  if (myrank.eq.1) then
     mpifileptr=6
  else
     mpifileptr=989
     open(mpifileptr,file="/dev/null", status="unknown")
  endif
end subroutine ctset



!! INVERSE OF cooleytukey_outofplace_mpi

subroutine cooleytukey_outofplace_inverse_mpi(intranspose,out,dim1,pf,proclist,localnprocs,localrank,howmany)
  implicit none
  integer, intent(in) :: dim1,localnprocs,localrank,howmany,pf(MAXFACTORS),proclist(localnprocs)
  complex*16, intent(in) :: intranspose(dim1,howmany)
  complex*16, intent(out) :: out(dim1,howmany)
  complex*16 ::  intransconjg(dim1,howmany),  outconjg(dim1,howmany)

  intransconjg(:,:)=CONJG(intranspose(:,:))
  call cooleytukey_outofplaceinput_mpi(intransconjg,outconjg,dim1,pf,proclist,localnprocs,localrank,howmany)
  out(:,:)=CONJG(outconjg(:,:))/dim1/localnprocs

end subroutine cooleytukey_outofplace_inverse_mpi



subroutine twiddlemult_mpi(in,out,dim1,numfactored,factorlist,localnumprocs,ctrank,howmany)
  use fileptrmod
  use ct_mpimod  !!TEMP
  implicit none
  integer, intent(in) :: dim1,howmany,localnumprocs,numfactored,ctrank,factorlist(numfactored)
  complex*16, intent(in) :: in(dim1,howmany)
  complex*16, intent(out) :: out(dim1,howmany)
  complex*16 :: twiddle1(dim1,numfactored),tt1(dim1)
  integer :: ii

  if (localnumprocs.ne.nprocs.or.numfactored.ne.1.or.ctrank.ne.myrank) then
     write(mpifileptr,*) "ACK, twiddlemult_mpi not done",localnumprocs,nprocs,numfactored,&
          ctrank,myrank; call mpistop()
  endif
  call gettwiddlesmall(twiddle1(:,:),dim1*numfactored,localnumprocs)

  tt1(:)=twiddle1(:,1)**(ctrank-1)   !!! ???? crux

  do ii=1,howmany
     out(:,ii) = in(:,ii) * tt1(:)
  enddo

end subroutine twiddlemult_mpi


subroutine checkdivisible(number,factor)
  use fileptrmod
  implicit none
  integer,intent(in) :: number,factor
  if ((number/factor)*factor.ne.number) then
     write(mpifileptr,*) "Divisibility failure", number,factor; call mpistop()
  endif
end subroutine checkdivisible


subroutine checkbetween(number,greaterthan,lessthan)
  use fileptrmod
  implicit none
  integer,intent(in) :: number,greaterthan,lessthan
  if (number.le.greaterthan.or.number.ge.lessthan) then
     write(mpifileptr,*) "Error checkbetween",greaterthan,number,lessthan; call mpistop()
  endif
end subroutine checkbetween


!! fourier transform with OUT-OF-PLACE OUTPUT. 

recursive subroutine cooleytukey_outofplace_mpi(in,outtrans,dim1,pf,proclist,localnprocs,localrank,howmany)
  use fileptrmod
  implicit none
  integer, intent(in) :: dim1,localnprocs,localrank,howmany,pf(MAXFACTORS),proclist(localnprocs/pf(1),pf(1))
  complex*16, intent(in) :: in(dim1,howmany)
  complex*16, intent(out) :: outtrans(dim1,howmany)
  complex*16 ::  tempout(dim1,howmany),  outtemp(dim1,howmany)
  integer :: depth, newrank, dim2, newpf(MAXFACTORS),newproclist(localnprocs/pf(1)),ctrank,ctset(pf(1))

  call checkdivisible(localnprocs,pf(1))
  call checkbetween(localrank,0,localnprocs+1)

  depth=localnprocs/pf(1)
  ctrank=(localrank-1)/depth+1
  newrank=mod(localrank-1,depth)+1
  ctset(:)=proclist(newrank,:)
  newproclist=proclist(:,ctrank)

  if (proclist(newrank,ctrank).ne.localrank) then
     write(mpifileptr,*) "RANK FAIL INVERSE",proclist(newrank,ctrank),localrank,newrank,depth,ctrank,pf(1); call mpistop()
  endif

  call myzfft1d_slowindex_mpi(in,tempout,pf(1),ctrank,ctset,dim1*howmany)
  call twiddlemult_mpi(tempout,outtemp,dim1,depth,newproclist,pf(1),ctrank,howmany)
  if (depth.eq.1) then
     call myzfft1d(outtemp,outtrans,dim1,howmany)
  else
     newpf(1:MAXFACTORS-1)=pf(2:MAXFACTORS); newpf(MAXFACTORS)=1

     call cooleytukey_outofplace_mpi(outtemp,outtrans,dim1,newpf,newproclist,depth,newrank,howmany)
  endif
end subroutine cooleytukey_outofplace_mpi



recursive subroutine cooleytukey_outofplaceinput_mpi(intranspose,out,dim1,pf,proclist,localnprocs,localrank,howmany)
  use fileptrmod
  implicit none
  integer, intent(in) :: dim1,localnprocs,localrank,howmany,pf(MAXFACTORS),proclist(localnprocs/pf(1),pf(1))
  complex*16, intent(in) :: intranspose(dim1,howmany)
  complex*16, intent(out) :: out(dim1,howmany)
  complex*16 ::   temptrans(dim1,howmany),outtrans(dim1,howmany)
  integer :: depth, newrank, newpf(MAXFACTORS),newproclist(localnprocs/pf(1)),ctrank,ctset(pf(1))

  call checkdivisible(localnprocs,pf(1))
  call checkbetween(localrank,0,localnprocs+1)

  depth=localnprocs/pf(1)
  ctrank=(localrank-1)/depth+1
  newrank=mod(localrank-1,depth)+1
  ctset(:)=proclist(newrank,:)
  newproclist=proclist(:,ctrank)


  if (proclist(newrank,ctrank).ne.localrank) then
     write(mpifileptr,*) "RANK FAIL INVERSE",proclist(newrank,ctrank),localrank,newrank,depth,ctrank,pf(1); call mpistop()
  endif


  if (depth.eq.1) then
     call myzfft1d(intranspose,temptrans,dim1,howmany)
  else
     newpf(1:MAXFACTORS-1)=pf(2:MAXFACTORS); newpf(MAXFACTORS)=1
     call cooleytukey_outofplaceinput_mpi(intranspose,temptrans,dim1,newpf,newproclist,depth,newrank,howmany)
  endif

  call twiddlemult_mpi(temptrans,outtrans,dim1,depth,newproclist,pf(1),ctrank,howmany)
  call myzfft1d_slowindex_mpi(outtrans,out,pf(1),ctrank,ctset,dim1*howmany)

end subroutine cooleytukey_outofplaceinput_mpi




subroutine myzfft1d_slowindex_mpi(in,out,localnumprocs,ctrank,proclist,totsize)
  implicit none
  integer, intent(in) :: totsize,localnumprocs,ctrank,proclist(localnumprocs)
  complex*16, intent(in) :: in(totsize)
  complex*16, intent(out) :: out(totsize)
  complex*16 :: fouriermatrix(localnumprocs,localnumprocs),twiddle(localnumprocs)
  integer :: ii

  call gettwiddlesmall(twiddle,localnumprocs,1)

  do ii=1,localnumprocs
     fouriermatrix(:,ii)=twiddle(:)**(ii-1)
  enddo

  call simple_summa(in,out,fouriermatrix,totsize,ctrank,localnumprocs,proclist)

end subroutine myzfft1d_slowindex_mpi




subroutine simple_circ(in, out,mat,howmany,ctrank,localnumprocs,proclist)
  use fileptrmod
  use ct_mpimod
  implicit none
  integer, intent(in) :: howmany,ctrank,localnumprocs,proclist(localnumprocs)
  complex*16, intent(in) :: in(howmany), mat(localnumprocs,localnumprocs)
  complex*16, intent(out) :: out(howmany)
  complex*16 :: work2(howmany),work(howmany)
  integer :: ibox,jbox,deltabox,nnn,CT_GROUP_LOCAL,CT_COMM_LOCAL,ierr,procshift(localnumprocs)

!!TEMP
  if (localnumprocs.ne.nprocs.or.myrank.ne.ctrank) then
     write(mpifileptr,*) "simple_Circ not done",localnumprocs,nprocs,myrank,ctrank; call mpistop()
  endif

  procshift(:)=proclist(:)-1
  call mpi_group_incl(CT_GROUP_WORLD,localnumprocs,procshift,CT_GROUP_LOCAL,ierr)
  if (ierr.ne.0) then
     write(mpifileptr,*) "Error group incl simple_circ",ierr; call mpistop()
  endif
  call mpi_comm_create(CT_COMM_WORLD, CT_GROUP_LOCAL, CT_COMM_LOCAL,ierr)
  if (ierr.ne.0) then
     write(mpifileptr,*) "Error comm create simple_circ",ierr; call mpistop()
  endif

  nnn=1

  out(:)=0

  do deltabox=0,nprocs-1

     ibox=mod(localnumprocs+ctrank-1+deltabox,localnumprocs)+1
     jbox=mod(localnumprocs+ctrank-1-deltabox,localnumprocs)+1

     work(:)=in(:)*mat(ibox,ctrank)

     if (deltabox.ne.0) then
        call mympisendrecv(work(:),work2(:),ibox,jbox,deltabox,howmany,CT_COMM_LOCAL)
        out(:)=out(:)+work2(:)
     else
        out(:)=out(:)+work(:)
     endif
  enddo

end subroutine simple_circ



subroutine simple_summa(in, out,mat,howmany,ctrank,localnumprocs,proclist)
  use fileptrmod
  use ct_mpimod
  implicit none
  integer, intent(in) :: howmany,ctrank,localnumprocs,proclist(localnumprocs)
  complex*16, intent(in) :: in(howmany), mat(localnumprocs,localnumprocs)
  complex*16, intent(out) :: out(howmany)
  complex*16 :: work(howmany)
  integer :: ibox,nnn,CT_GROUP_LOCAL,CT_COMM_LOCAL,ierr,procshift(localnumprocs)

!!TEMP
  if (localnumprocs.ne.nprocs.or.myrank.ne.ctrank) then
     write(mpifileptr,*) "simple_Circ not done",localnumprocs,nprocs,myrank,ctrank; call mpistop()
  endif

  procshift(:)=proclist(:)-1
  call mpi_group_incl(CT_GROUP_WORLD,localnumprocs,procshift,CT_GROUP_LOCAL,ierr)
  if (ierr.ne.0) then
     write(mpifileptr,*) "Error group incl simple_circ",ierr; call mpistop()
  endif
  call mpi_comm_create(CT_COMM_WORLD, CT_GROUP_LOCAL, CT_COMM_LOCAL,ierr)
  if (ierr.ne.0) then
     write(mpifileptr,*) "Error comm create simple_circ",ierr; call mpistop()
  endif

  nnn=1

  out(:)=0d0
  do ibox=1,localnumprocs
     if (ctrank.eq.ibox) then
        work(:)=in(:)
     endif
     call mympibcast(work(:),ibox,howmany,CT_COMM_LOCAL)
     out(:)=out(:)+work(:)*mat(ctrank,ibox)
  enddo

end subroutine simple_summa
