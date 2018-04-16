c$Id:
      subroutine elmt06(d, ul, xl, ix, tl, s, r, ndf, ndm, nst, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c      Modification log                                Date (dd/mm/year)
c      Original version                                    04/08/2015
c      Purpose: 3-D electrophysiological element
c      Remark: This a completely electrical formulation
c-----[--.----+----.----+----.-----------------------------------------]

c     Inputs:
c     d(*)          - Element data parameters (Moduli, body loads, etc)
c     ul(ndf,nen,j) - Element nodal solution parameters (nen is number of nodes on an element)
c                     j = 1: Displacement
c                     j = 2: Increment 
c     xl(ndm,nen)   - Element nodal reference coordinates
c     ix(nen)       - Element global node numbers
c     tl(nen)       - Element nodal temperature numbers
c     s(nst,nst)    - Element matrix (e.g., stiffness, mass)
c     r(nst)        - Element vector (e.g., residual, mass). May also be used as r(nst)
c     ndf           - Number unknowns (max) per node
c     ndm           - Space dimension of mesh
c     nst           - Size of element arrays S and R. N.B. Normally nst = ndf*nen
c     isw           - Task parameter to control computation
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'        ! o, head
      include  'cdata.h'        ! numnp, numel, nummat, nen, neq, ipr
      include  'cdat1.h'        ! ndd, nie, nud
      include  'eldata.h'       ! dm, n, ma, mct, iel, nel, pstyp, eltyp
      include  'hdata.h'        ! nh1, nh2, nh3, ht1, ht2, ht3
      include  'iofile.h'       ! ior, iow, ilg
      include  'pmod2d.h'       ! stype, etype, dtype, plasfl, viscfl, hflag, gflag
      include  'tdata.h'        ! ttim, dt, c1, c2, c3, c4, c5, chi, dtcr
      include  'comblk.h'       ! hr, mr
      include  'strnum.h'       ! iste, istv
      include   'part0.h'       ! npart
      include  'qudshp.h'       ! quad, ttfl, nurbfl, jac, lint, npm, nvm, sg2, el2, shp2, sg3, el3, shp3
      include 'counts.h'
      include  'pointer.h'      !np
      
c     Inputs/outputs
      integer   ix(*), ndf, ndm, nst, isw
      real*8    d(*), ul(ndf,nen,*), xl(ndm,*), tl(*)
      real*8    s(nst,nst), r(nst)
      
c     Dummies
      integer   a, b, i
      integer   l, ord
      integer   n1, n3, nhlo, nhln
      integer   imat
      real*8    sv(5,16),dlon,dtra
      real*8    iionk(nel),iionk_1(nel),diionk_1(nel)
      real*8    xr(3), xh(3),xhn(3),theta
      real*8    scratch(4,64)
      real*8    ua(nel),uan(nel),duan(nel)
      real*8    uak(nel),uak_1(nel)
      real*8    res, xtol        

      save
      
      
c------------------------------------!
c     OUTPUT ELEMENT DESCRIPTION     !
c         (isw = 0)                  !
c------------------------------------!
      if(isw.eq.0 .and.ior.lt.0) then
         write(*,*) 'Elmt 3: Monodomain model 3D'
      endif

c-------------------------------------------------------------!
c     INPUT/OUTPUT OF PROPERTY DATA AFTER COMMAND: 'mate'     !
c                          (isw = 1)                          !
c-------------------------------------------------------------!
      if(isw.eq.1) then
         
         write(*,2000)
         
         call inmate(d, tl, 0, 8)
         
         n1  = nh1
         n3  = nh3
         lint = 4 ! assuming order 2 and tet elements!
         nh1 = lint*n1
         nh3 = lint*n3
         
         stype = nint(d(16))
         etype = nint(d(17))
         dtype = nint(d(18))
         
!     Set plot sequence for 4/10-nod tet 
         if(nen.eq.4) then
            call pltet4(iel)
         elseif(nen.eq.10) then
            call pltet10(iel)
         endif
         
c---------------------------------!
c     CHECK ELEMENT FOR ERRORS    !
c            (isw = 2)            !
c---------------------------------!
      elseif(isw.eq.2) then
         
         if(nel.eq.4 .or. nel.eq.10) then
            call cktets ( n, ix, xl, ndm, nel, scratch )
         else
            call ckbrk8 ( n, ix, xl, ndm, nel, scratch )
         endif         
         
c-------------------------------------------------------------!
c      COMPUTE TANGENT STIFFNESS AND RESIDUAL FORCE VECTOR    !
c                       (isw = 3 , 6)                         !
c      OUTPUT ELEMENT VARIABLES AND ITS NODAL PROJECTIONS     !
c                       (isw = 4 , 8)                         !
c                        J INTEGRALS                          !
c                        (isw = 16)                           !
c-------------------------------------------------------------!
      elseif (isw.eq.3 .or.  isw .eq.4 .or.
     &        isw .eq.8 .or. isw.eq.6 .or. isw.eq.16) then
     
      
        
        dlon = d(ndd-nud+1)*d(ndd-nud+3)/(d(ndd-nud+1)+d(ndd-nud+3))
         dtra = d(ndd-nud+2)*d(ndd-nud+4)/(d(ndd-nud+2)+d(ndd-nud+4))
        
          
      if(npart .eq. 1 .or. npart .eq. 3 ) then ! if partition 1 - diffusive part
      
        
        
       theta = 0.5
c ..... diffusive part
       do a=1,nel
           duan(a) = ul(1,a,2)
           ua(a) = ul(1,a,1)
           uan(a) = ul(1,a,1)-ul(1,a,2)
           
       enddo

C        if (n .eq. 1) then
C            write(*,*) '------PART 1 DIFF-----'
C            write(*,*) 'iteration :', niter
C            do a=1,nel
C            write(*,*) 'iionk: ',tl(a)
C            enddo
C          endif
       
c ..... diffusive part

         call diffusion_os(dlon,dtra,s,r,ua,uan,duan,
     &                 theta,dt*0.5,xl,nel,ndm,nst,isw)
     
         if (isw .eq. 8) then
                  do a=1,nel
                     r(a) = r(a) + 1.0
                     s(a,1) = s(a,1) + tl(a)
                  enddo
                     iste = 1
                     istv = 1
         endif
c ..... iionic part
       
       elseif(npart .eq. 2) then
       
         
         theta = 0.5
         imat = nint(d(ndd-nud))
         iionk_1 = 0.d0
         diionk_1 = 0.d0
         
         do a=1,nel
         
           iionk = tl(a)
           uak_1(a) = ul(1,a,1)
           uak(a) = ul(1,a,1)-ul(1,a,2)
       
         enddo

C          if (n .eq. 1) then
C            write(*,*) '------PART 2 IONC-----'
C            write(*,*) 'iteration :', niter
C            do a=1,nel
C            write(*,*) 'iionk: ',iionk(a)
C            enddo
C          endif

         call  ionic_os2(s,r,uak_1,uak,
     &               d(ndd-nud+1),iionk,
     &               n1,theta,dt,imat,isw)
                  
         endif !npart  
c-------------------------------------------------------------!
c            INITIALIZE ELECTROPHYSIOLOGICAL MODELS           !
c                        (isw = 14)                           !
c-------------------------------------------------------------!
       elseif (isw.eq.14 .or. isw.eq.12) then
         
         if(isw. eq. 14) then
         imat = nint(d(ndd-nud))
         
         do a=1,nel
            nhlo = (a-1)*n1 + nh1
            nhln = (a-1)*n1 + nh2
            
            call calc_material(imat,ul(1,a,1),d(ndd-nud+1), 
     &                   hr(nhlo), hr(nhln)
     &                   ,iionk_1(a),diionk_1(a),n1,isw)
            
            dlon = d(ndd-nud+1)*d(ndd-nud+3)/(d(ndd-nud+1)+d(ndd-nud+3))
            dtra = d(ndd-nud+2)*d(ndd-nud+4)/(d(ndd-nud+2)+d(ndd-nud+4))
            
            hr(np(38)+(ix(a)-1))= 0.d0
           
         enddo

        endif
        if (isw .eq. 12 ) then

         if (npart .eq. 1)then
          do a=1,nel
            nhlo = (a-1)*n1 + nh1
            nhln = (a-1)*n1 + nh2
            
             call calc_material(imat,ul(1,a,1),d(ndd-nud+1), 
     &                   hr(nhlo), hr(nhln)
     &                   ,iionk_1(a),diionk_1(a),n1,isw)
            
             hr(np(38)+(ix(a)-1))= iionk_1(a)
          enddo
         
        endif
C          if (n .eq. 1) then
C            write(*,*) '------TIME-----'
C            write(*,*) 'npart: ',npart
C            do a=1,nel
C            write(*,*) 'iiok: ',iionk_1(a)
C            enddo
C          endif

         endif
      endif
      
c     FORMAT 
      
 2000 format(
     &     /5x,'3 - D   E l e c t r o p h y s i o l o g i c a l',
     &     '   T e t r a h e d r a l    E l e m e n t'/)
      
 2001 format(a1,20a4//5x,'Element Variables'//
     &     '   Elem.   Matl.   1-Coord   2-Coord   3-Coord',
     &     '   fphi   phidot   1-gphi    2-gphi    3-gphi',
     &     '   1-flux    2-flux    3-flux')
      
 2002 format(i8,i8,1p,6e12.4/i8,1p,6e12.4/1x)
 
      end


     


