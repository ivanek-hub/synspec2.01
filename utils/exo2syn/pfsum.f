      program isosum
!     ==============
!
!     evaluation of isotope-averaged partition function for 
!     diatomic molecules from EXOMOL data
!
      implicit none

      integer, parameter :: mt=10000,mato=26,miso=58,mcomp=16,melem=5
      real :: tem(mt,mcomp),pf(mt,mcomp),partf(mt),pf0(mt)
      real :: fra(mato,miso),t(mt),fract(mcomp),tin(mt),pin(mt)
      real :: dtem(mcomp),pfint(mt)
      integer ::  ntem(mcomp)
      character(len=30) :: filepf

      real    :: fr,frac,fratot,framax,tt,pp,tmin,tmax
      integer :: iatm,isot,i,iat1,iat2,iso,iso1,iso2,isomax,ios,it,       &
     &           niso,ntmx,ntmin,ntmax

!
!     read (solar) isotopic fractions
!
      open(unit=7,file='isotops',status='old')
      readfrac: do
         read(7,*,iostat=ios) iatm,isot,fr
         if(ios.ne.0) exit
         fra(iatm,isot)=fr*0.01
      end do readfrac
      close(7)
!
      do i=1,mt
         partf(i)=0.
      end do
      fratot=0.
      read(5,*) iat1,iat2
      iso=0
      isomax=1
      framax=0.
      iso=0
      isotops: do
! 
!        read partition funcion foor all isotopologues
!
         read(5,*,iostat=ios) iso1,iso2,filepf
         if(ios.ne.0) exit isotops
         frac=fra(iat1,iso1)*fra(iat2,iso2)
         fratot=fratot+frac
         iso=iso+1
         fract(iso)=frac
         if(frac.ge.framax) then
            framax=frac
            isomax=iso
         end if
         open(unit=7,file=filepf,status='old')
         it=0
         readpf: do
            read(7,*,iostat=ios) tt,pp
            if(ios.ne.0) exit readpf
            it=it+1
            tem(it,iso)=tt
            pf(it,iso)=pp
         end do readpf
         ntem(iso)=it
         dtem(iso)=tem(2,iso)-tem(1,iso)
         close(7)
         write(6,"('atom1',2i5,'    atom2',i7,2i6,f6.1,4f8.4)")           &
     &         iat1,iso1,iat2,iso2,ntem(iso),dtem(iso),                   &
     &         fra(iat1,iso1),fra(iat2,iso2),frac,fratot
      end do isotops
      niso=iso
      write(*,*) 'niso',niso
!
!     *** the case of just one isotopologe
!
      if(niso.eq.1) then
         iso=1
         do it=1,ntem(iso)
            t(it)=tem(it,iso)
            partf(it)=pf(it,iso)
         end do
         ntmx=ntem(iso)
         write(*,*) 'ntmx',iso,ntmx,t(1),t(ntmx)
       else
!
!     *** more isotopologes
!
         ntmin=ntem(1)
         ntmax=ntem(1)
         tmin=tem(1,1)
         tmax=tem(ntem(1),1)
         do iso=2,niso
            ntmin=min(ntmin,ntem(iso))
            ntmax=max(ntmax,ntem(iso))
            tmin=min(tmin,tem(1,iso))
            tmax=max(tmax,tem(ntem(iso),iso))
         end do
         ntmx=ntem(isomax)
!
!        the case of all the isotop. tables have the same T_min and T_max
!        and also the same numer of temperatures - no interpolation
!
         if((tmax-tmin)/tmin.lt.0.0001.and.ntmax.eq.ntmin) then
            do it=1,ntmax
               t(it)=tem(it,1)
               do iso=1,niso
                  partf(it)=partf(it)+pf(it,iso)*fract(iso)
               end do
            end do
!
!         other cases - some interpolations needed
!
          else
            do it=1,ntmx
               t(it)=tem(it,isomax)
               partf(it)=pf(it,isomax)
               pf0(it)=partf(it)
            end do
            isotp: do iso=1,niso
               if(iso.eq.isomax) exit isotp
               do it=1,ntem(iso)
                  tin(it)=tem(it,iso)
                  pin(it)=pf(it,iso)
               end do
               call interp(tin,pin,t,pfint,ntem(iso),ntem(isomax),        &
     &                     2,0,0)
               do it=1,ntem(isomax)
                  if(t(it).lt.tin(1)) pfint(it)=0.
                  if(t(it).gt.tin(ntem(iso))) pfint(it)=0.
                  partf(it)=partf(it)+pfint(it)*fract(iso)
               end do
            end do isotp
         end if
      end if
!
!     output of isotopologue-averaged partition function on unit 10
!
      do it=1,ntmx
         write(10,"(f8.1,f16.4)") t(it),partf(it)/fratot
      end do
      close(10)
!
      end program isosum
!
!           
!      ******************************************************************
!


      SUBROUTINE INTERP(X,Y,XX,YY,NX,NXX,NPOL,ILOGX,ILOGY)
!     ====================================================
!
!     General interpolation procedure of the (NPOL-1)-th order
!
!     for  ILOGX = 1  logarithmic interpolation in X
!     for  ILOGY = 1  logarithmic interpolation in Y
!
!     Input:
!      X    - array of original x-coordinates
!      Y    - array of corresponding functional values Y=y(X)
!      NX   - number of elements in arrays X or Y
!      XX   - array of new x-coordinates (to which is to be
!             interpolated
!      NXX  - number of elements in array XX
!     Output:
!      YY   - interpolated functional values YY=y(XX)
!
      integer, parameter :: mx=10000
      real :: X(MX),Y(MX),XX(MX),YY(MX)
!
!     no interpolation for NPOL.LE.0 or NX.le.0
!
      IF(NPOL.LE.0.OR.NX.LE.0) THEN
         N=NX
         IF(NXX.GE.NX) N=NXX
         DO I=1,N
            XX(I)=X(I)
            YY(I)=Y(I)
         END DO
         RETURN
      END IF
!
!     interpolation
!
!     if required, compute logarithms to be interpolated
!
      IF(ILOGX.GT.0) THEN
      DO I=1,NX
         X(I)=LOG10(X(I))
      END DO
      DO I=1,NXX
         XX(I)=LOG10(XX(I))
      END DO
      END IF
      IF(ILOGY.GT.0) THEN
      DO I=1,NX
         Y(I)=LOG10(Y(I))
      END DO
      END IF
!
      NM=(NPOL+1)/2
      NM1=NM+1
      NUP=NX+NM1-NPOL
      DO ID=1,NXX
         XXX=XX(ID)
         DO I=NM1,NUP
            IF(XXX.LE.X(I)) EXIT
         END DO
         J=I-NM
         JJ=J+NPOL-1
         YYY=0.
         DO K=J,JJ
            T=1.
            DO M=J,JJ
               IF(K.EQ.M) CYCLE
               T=T*(XXX-X(M))/(X(K)-X(M))
            END DO
            YYY=Y(K)*T+YYY
         END DO
         YY(ID)=YYY
      END DO
!
      IF(ILOGX.GT.0) THEN
      DO I=1,NX
         X(I)=EXP(X(I)*2.30258509299405)
      END DO
      DO I=1,NXX
         XX(I)=EXP(XX(I)*2.30258509299405)
      END DO
      END IF
      IF(ILOGY.GT.0) THEN
      DO I=1,NX
         Y(I)=EXP(Y(I)*2.30258509299405)
      END DO
      DO I=1,NXX
         YY(I)=EXP(YY(I)*2.30258509299405)
      END DO
      END IF
! 
      RETURN
      END SUBROUTINE INTERP

