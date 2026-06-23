      include 'eos_modules.f'

      program eosmol
      use accura
      use params
      use modelp
      use eospar
      implicit real(dp) (a-h,o-z)
      integer :: ielpri(15)
      data ielpri/1,6,7,8,9,11,12,13,14,15,16,20,22,24,26/

      call eos_alloc

      ipri=1
      ifmol=1
      tmolim=10001.
      moltab=1
      irwtab=1
      ieqbc=1
      nd=1
      read(5,*) t,an
      if(t.gt.1.e4) then
        write(*,*) 'T > 10000; stop'
        stop
      end if
     
      write(*,"(/'************************************'/                  & 
     &           'SOLUTION OF THE CHEMICAL EQUILIBRIUM'/                  &
     &           '************************************'/                  &
     &           ' FOR T, N =' f10.1,1pe10.1/)") t,an
      temp(1)=t
      call anesti(t)
      aein=an*anerel
      call state0
      call moleq(t,an,aein,ane,ipri)
!     write(6,"('t,an,ane ',f10.1,1p2e12.3)") t,an,ane
      call eospri

      do i=1,15
         call elebear(ielpri(i))
      end do

      close(33)
      close(34)
      end program eosmol
!
      include 'eos_alloc.f'
!
!     ***********************************************************************
!
      subroutine anesti(t)

      use accura
      use params

      real(dp) :: t

      if(t.ge.9000.) then
         anerel=0.5
       else if(t.ge.4000.) then
         anerel=0.5*10.**(-9.+t/1000.)
       else if(t.ge.1000.) then 
         anerel=1.e-6
       else
         anerel=1.e-11
      end if
!     write(*,*) 'anerel',anerel

      end subroutine anesti
    

!
!     ***********************************************************************
!

      subroutine moleq(tt,an,aein,ane,ipri)
!     ========================================
!
!     calculation of the equilibrium state of atoms and molecules
!
!             tt    - temperature [K]
!             an    - number density
!             aein  - initial estimate of the electron density
!
!     Output: ane    - electron density
!
!     Input data for molecules  iven in the file
!     tsuji.molec
!
      use accura
      use params
      use modelp
      use eospar
      implicit real(dp) (a-h,o-z),logical (l)

      character(len=128) :: MOLEC
      INTEGER        :: NATOMM(5),NELEMM(5)
      REAL(DP), SAVE :: EMASS(100)
      REAL(DP)       :: uelem(100),ull(100),anden(800),aelem(100),        &
     &                  denso(mdepth),eleco(mdepth),wmmo(mdepth)
!
!     data nmetal/92/
!
      data iread/1/
!
      MOLEC ='data/tsuji.molec_bc2'
      if(moltab.eq.0) MOLEC='data/tsuji.molec_orig'
!
        ECONST=4.342945E-1
        AVO=0.602217E+24
        SPA=0.196E-01
        GRA=0.275423E+05
        AHE=0.100E+00
        tk=1./(tt*1.38054e-16)
        pgas=an/tk
        sahcon=1.87840e20*tt*sqrt(tt)
!       nimax=3000
!       nimax=50
        eps=1.e-5
        switer=0.0
        id=1
!
!---- data for atoms  ----------------
!
      if(iread.eq.1) then
!
      do i=1,nmetal
         ia=i
         nelemx(i)=ia
         ccomp(ia)=abndd(ia,id)
         xip(ia)=enev(ia,1)
         xi2(ia)=enev(ia,2)
         emass(ia)=amas(ia)
      end do
!
!---- read molecular data from a table  ----------------------
!
        J=0
        OPEN(UNIT=26,FILE=MOLEC,STATUS='OLD')
        READMOL: DO
           J=J+1
           IF(MOLTAB.GE.1)                                                &
     &     READ (26,"(a8,5e13.5,9i3)",IOSTAT=IOS) CMOL(J),                &
     &        (C(J,K),K=1,5),MMAX(J),(NELEMM(M),NATOMM(M),M=1,4)
           IF(MOLTAB.EQ.0)                                                &
     &     READ (26,"(A8,E11.5,4E12.5,I1,(I2,I3),3(I2,I2))",IOSTAT=IOS)   &
     &        CMOL(J),(C(J,K),K=1,5),MMAX(J),                             &
     &        (NELEMM(M),NATOMM(M),M=1,4)
!
           IF(IOS.NE.0) EXIT READMOL
           MMAXJ=MMAX(J)
           IF(MMAXJ.EQ.0) THEN
              EXIT READMOL
            ELSE
              DO M=1,MMAXJ
                 NELEM(M,J)=NELEMM(M)
                 NATO(M,J)=NATOMM(M)
              END DO
!             write(6,"(i5,a10)") j,cmol(j)
              CYCLE READMOL
           END IF
        END DO READMOL

        NMOLEC=J-1
        close(26)
!
        DO I=1,NMETAL
           NELEMI=NELEMX(I)
           P(NELEMI)=1.e-70_dp
        END DO
        iread=0
      endif
!
!---- end of reading atomic and molecular data  ----------------------
!
        p(99)= aein/tk
        pesave=p(99)
        p(99)=pesave
!
        THETA=5040./tt
        TEM=tt
        PGLOG=log10(Pgas)
        PG=Pgas
!
        CALL RUSSEL(TEM,PG)
!
        PE=P(99)
        ane=pe*tk
        PELOG=log10(PE)
        emass(99)=5.486e-4
        uelem(99)=2.
        aelem(99)=pe*tk/(2.*sahcon*emass(99)**1.5)
        ull(99)=log10(aelem(99))
!
!----atoms-----------------------------------------------------------------
!
        tmass=0.
        atot=0.
        DO I=1,NMETAL
           NELEMI=NELEMX(I)
           FPLOG=log10(FP(NELEMI))
           anden(i)=(p(nelemi)+1.e-70_dp)*tk
           tmass=tmass+anden(i)*emass(nelemi)
           call irwpf(nelemi,1,0,tt,u0)
           uelem(nelemi)=u0
           aelem(nelemi)=anden(i)/(u0*sahcon*emass(nelemi)**1.5)
           ull(nelemi)=log10(aelem(nelemi))
           rrr(id,1,nelemi)=anden(i)/u0
           anato(nelemi,id)=anden(i)
           pfato(nelemi,id)=u0
           atot=atot+anden(i)
        END DO
        an1=anden(1)
!
!---- positive ions ---------------------------------------------------------
!
        DO I=1,NMETAL
           NELEMI=NELEMX(I)
           PLOG= log10(P(NELEMI)+1.0e-70_dp)
           XKPLOG=log10(XKP(NELEMI)+1.0e-70_dp)
           PIONL=PLOG+XKPLOG-PELOG
           anden(i+nmetal)=exp(pionl/econst)*tk
           tmass=tmass+anden(i+nmetal)*emass(nelemi)
           call irwpf(nelemi,2,0,tt,u1)
           anion(nelemi,id)=anden(i+nmetal)
           pfion(nelemi,id)=u1
           rrr(id,2,nelemi)=anden(i+nmetal)/u1
           if(nelemi.ge.2.and.nelemi.le.30) then
              x2log=log10(XK2(NELEMI)+1.0e-70_dp)
              pion2=pionl+x2log-pelog
              anion2(nelemi,id)=exp(pion2/econst)*tk
           end if
           atot=atot+anden(i+nmetal)+anion2(nelemi,id)

           if(id.eq.1) write(78,"(i5,a4,1p2e11.3,0p2f10.3)")              & 
     &        i,typat(i),anato(i,id),anion(i,id),                         &
     &        log10(anato(i,id)),log10(anion(i,id))
!          if(id.eq.1) write(*,"('moleq',i5,a4,5f10.3)")                  & 
!    &      i,typat(i),log10(anion(i,id)),                                &
!    &      plog,xkplog,pelog,pionl
        END DO
        anion2(1,id)=0.
!
!---- molecules-------------------------------------------------------------
!
        DO J=1,NMOLEC
           jm=j+2*nmetal
           PMOLL=log10(PPMOL(J)+1.0e-70_dp)
           anden(jm)=exp(pmoll/econst)*tk
           rrmol(j,id)=0.
           umoll=1.
           if(pmoll.gt.-60.) then
              umoll=log10(anden(jm))+c(j,2)*theta
              amasm=0.
              do jjj=1,mmax(j)
                 i=nelem(jjj,j)
                 amasm=amasm+NATO(jjj,j)*emass(i)
                 umoll=umoll-NATO(jjj,j)*ull(i)
              end do
              ammol(j)=amasm
              tmass=tmass+anden(jm)*amasm
              umoll=exp(umoll/econst)/(sahcon*amasm**1.5)
!
!     replace with EXOMOL data whenever available
!
              um=0.
              if(ipfexo.gt.0.and.ipfbc.le.1.and.tt.le.9000.)              &
     &        call exopf(j,tt,um)
              if(um.gt.0.) then
                 umoll=um
               else
!
!     or with modified Irwin (Barklem & Collet) data whenever available
!
                 call irwpf(0,0,j,tt,um)
                 if(um.gt.0.) umoll=um
              end if
              if(ipfbc.eq.1.and.tt.lt.1000..or.ipfbc.eq.2) then
                 call bcdata(j,1,tt,um)
                 if(um.gt.0.) umoll=um
              end if
!
!     H-
!
              if(j.eq.1) umoll=1.
!
!     set up array RRR = number density/partition function
!
              rrmol(j,id)=anden(jm)/umoll
           end if
!
           anmol(j,id)=anden(jm)
           pfmol(j,id)=umoll
           atot=atot+anden(jm)

         if(id.eq.1) write(79,"(i4,a8,1pe12.3,0pf10.3)")                     &
     &     j,cmol(j),anmol(j,id),log10(anmol(j,id))

!
!         write(6,"(i5,a8,2f10.3,1p2e12.4,0pf12.3,1p2e12.4)")                &
!    &          j,cmol(j),apmlog(j),pmoll,anden(jm),umoll,                   &
!    &          amasm,tmass,atot
        END DO

        jm=2*nmetal
        anhm(id)=anden(1+jm)
        anh2(id)=anden(2+jm)
        anch(id)=anden(5+jm)
        anoh(id)=anden(4+jm)
        ahn(id)=anato(1,id)
        ahp(id)=anion(1,id)
        ahen(id)=anato(2,id)
!
!
!     save new density, molecular weight, and abundances of
!     atomic species
!
      ipri1=ipri
      denso(id)=dens(id)
      eleco(id)=elec(id)
      wmmo(id)=wmm(id)
      dens(id)=tmass*hmass
      elec(id)=pe*tk
      wmm(id)=dens(id)/(an-elec(id))
      ane=elec(id)
!
      do i=1,nmetal
         NELEMI=NELEMX(I)
         ia=iatex(nelemi)
         if(ia.gt.0) then
            attot(ia,id)=(anato(nelemi,id)+anion(nelemi,id))
         end if
      end do

      close(78)
      close(79)
!
      RETURN
      END SUBROUTINE MOLEQ

!
!     ***********************************************************************
!

      SUBROUTINE RUSSEL(TEM,PG)
!     =========================
!
      use accura
      use params
      use modelp
      use eospar
      implicit real(dp) (a-h,o-z),logical (l)

      REAL(DP) ::  FX(100),DFX(100),Z(100),PREV(100),WA(100),             &
     &             UIIDU2(100)
!
      ECONST=4.3426E-1
      XKCON=6.667343E-1
      EPSDIE=5.0E-5
      T=5040.4/TEM
      PGLOG=log10(PG)


      tk=1./(tem*1.38054e-16)
!
!    HEH=helium/hydrogen ratio by number
!
      HEH=CCOMP(2)/CCOMP(1)
!     HEH=YTOT(1)-UN
!
!    evaluation of log XKP(MOL)
!

      call equcon(tem)
!
      apmlog(1)=-log10(1.0353e-16/tem/sqrt(tem)*tk*exp(8762.9/tem))
      DHH=APMLOG(2)
      DHH=EXP(DHH/ECONST)
!
!  evaluation of the ionization constants
!
      TEM25=TEM**2*SQRT(TEM)
      DO I=1,NMETAL
         NELEMI = NELEMX(I)
!
! calculation of the partition functions following Irwin (1981)
!
         call irwpf(nelemi,1,0,tem,g0)
         call irwpf(nelemi,2,0,tem,g1)
         call irwpf(nelemi,3,0,tem,g2)
!        uiidui(nelemi)=g1/g0*0.6665
         uiidui(nelemi)=g1/g0*xkcon
         uiidu2(nelemi)=g2/g1*xkcon
!
         XKP(NELEMI)=UIIDUI(NELEMI)*TEM25*                              &
     &               EXP(-XIP(NELEMI)*T/ECONST)
         XK2(NELEMI)=UIIDU2(NELEMI)*TEM25*                              &
     &               EXP(-XI2(NELEMI)*T/ECONST)
         xk2(nelemi)=max(xk2(nelemi),1.e-70_dp)
      END DO
      HKP=XKP(1)
      XK2(1)=0.
!
!   preliminary value of PH at high temperatures
!
      HKP=XKP(1)
      IF(T.LT.0.6) THEN
         PPH=SQRT(HKP*(PG/(1.0+HEH)+HKP))-HKP
         PH=PPH**2/HKP
      END IF
!
!  evaluation of the fictitious pressures of hydrogen
!     PG=PH+PHH+2.0*PPH+HEH*(PH+2.0*PHH+PPH)
!
      U=(1.0+2.0*HEH)/DHH
      Q=1.0+HEH
      R=(2.0+HEH)*SQRT(HKP)
      S=-1.0*PG
!
!     first estimate of PH at low temperatures
!
      IF(T.GE.0.6) PH=(-Q+SQRT(Q**2-4.*U*S))/(2.*U)
      X=SQRT(PH)
!
!     Russell iterations
!
      ITERAT=0
      ITRAS: DO
         F=((U*X**2+Q)*X+R)*X+S
         DF=2.0*(2.0*U*X**2+Q)*X+R
         XR=X-F/DF
         IF(ABS((X-XR)/XR).GT.EPSDIE) THEN
            ITERAT=ITERAT+1
            IF(ITERAT.GT.50) THEN
               WRITE(6,"(/' NOT CONVERGE IN RUSSEL '                    &
     &         /// 'TEM=',F9.2,5X,'PG=',E12.5,5X,                       &
     &         'X1=',E12.5,5X,'X2=',E12.5,5X,'PH=',E12.5/////)")        &
     &         TEM,PG,X,XR,PH
               EXIT ITRAS
             ELSE
               X=XR
               CYCLE ITRAS
            END IF
          ELSE
            EXIT ITRAS
         END IF
      END DO ITRAS

      PH=XR**2
      PHH=PH**2/DHH
      PPH=SQRT(HKP*PH)
      FPH=PH+2.0*PHH+PPH
      P(100)=PPH
!
!   evaluation of the fictitious pressure of each element
!
      DO I=1,NMETAL
         NELEMI=NELEMX(I)
         FP(NELEMI)=CCOMP(NELEMI)*FPH
      END DO
!
      PE=P(99)
!
!     Russell equations
!
      NITERR = 0
      RUSS: DO
        DO I=1,NMETAL
           NELEMI=NELEMX(I)
!          FX(NELEMI)=-FP(NELEMI)+P(NELEMI)*(1.0+XKP(NELEMI)/PE)
           DFX(NELEMI)=1.0+XKP(NELEMI)/PE*(1.0+XK2(NELEMI)/PE)
           FX(NELEMI)=-FP(NELEMI)+P(NELEMI)*DFX(NELEMI)
        END DO
!
        SPNION=0.0
        spnplu=0.
        DO J=1,NMOLEC
           lpri=tem.lt.4103..and.(j.le.3.or.j.eq.444)
           MMAXJ=MMAX(J)
           PMOLJL=-APMLOG(J)
           natot=0
           DO M=1,MMAXJ
              NELEMJ=NELEM(M,J)
              NATOMJ=NATO(M,J)
              natot=natot+natomj
              PMOLJL=PMOLJL+DFLOAT(NATOMJ)*log10(P(NELEMJ))
              if(nelemj.eq.90) pmoljl=-60.
           END DO
           if(natot.ge.6) pmoljl=-60.
!
           PMOLJ=EXP(PMOLJL/ECONST)
           DO M=1,MMAXJ
              NELEMJ=NELEM(M,J)
              NATOMJ=NATO(M,J)
              ATOMJ=DFLOAT(NATOMJ)
              IF(NELEMJ.EQ.99) then
                 if(natomj.ge.0) then
                    SPNION=SPNION+PMOLJ*NATOMJ
                  else
                    SPNPLU=SPNPLU-PMOLJ*NATOMJ
                 end if
!             if(lpri)
!    &        write(6,"('ion',3i5,a8,1p4e15.5)")                          &
!    &        niterr,j,m,cmol(j),p(99),pmolj,spnplu,spnion                
              end if
              DO I=1,NMETAL
                 NELEMI=NELEMX(I)
                 IF(NELEMJ.EQ.NELEMI) THEN
                    FX(NELEMI)=FX(NELEMI)+ATOMJ*PMOLJ
                    DFX(NELEMI)=DFX(NELEMI)+ATOMJ**2*                     &
     &                          PMOLJ/P(NELEMI)
                 END IF
              END DO
           END DO
           PPMOL(J)=PMOLJ
!          if(tem.lt.3650.)
!    &     write(66,"(2i4,a8,i3,2f10.3,2x,5(2i3,f8.2,2x))")              &
!    &     niterr,j,cmol(j),mmaxj,apmlog(j),pmoljl,                      &
!    &     (nelem(m,j),nato(m,j),log10(p(nelem(m,j))),m=1,mmaxj)

        END DO
!
!   solution of the Russell equations by Newton-Raphson method
!
        DO I=1,NMETAL
           NELEMI=NELEMX(I)
           WA(I)=log10(P(NELEMI)+1.0e-70_dp)
        END DO
        IMAXP1=NMETAL+1
        WA(IMAXP1)=log10(PE+1.0e-70_dp)
        DELTRS = 0.0
        DO I=1,NMETAL
           NELEMI=NELEMX(I)
           PREV(NELEMI)=P(NELEMI)-FX(NELEMI)/DFX(NELEMI)
           PREV(NELEMI)=ABS(PREV(NELEMI))
           IF(PREV(NELEMI).LT.1.0D-70) PREV(NELEMI)=1.0e-70_dp
           Z(NELEMI)=PREV(NELEMI)/P(NELEMI)
           DELTRS=DELTRS+ABS(Z(NELEMI)-1.0)
           IF(SWITER.GT.0.0) THEN
              P(NELEMI)=(PREV(NELEMI)+P(NELEMI))*0.5
            ELSE
              P(NELEMI)=PREV(NELEMI)
           END IF
        END DO
!
!   ionization equilibrium
!
        PEREV = spnplu
        DO I=1,NMETAL
           NELEMI = NELEMX(I)
           PEREV=PEREV+XKP(NELEMI)*P(NELEMI)*(1.+xk2(nelemi)/pe)
        END DO
!
        perev0=perev
!!      write(6,"('niterr',i4,1p4e12.4)")                                 &
!!   &  niterr,perev0,spnion,pe,spnion/pe
        PEREV=SQRT(PEREV/(1.0+SPNION/PE))
        DELTRS=DELTRS+ABS((PE-PEREV)/PE)
!!      write(6,"('niterr',i4,1p7e12.4)")                                 &
!!   &  niterr,tem,pg*tk,fph*tk,pe*tk,perev*tk,                           &
!!   &    (perev+pe)*0.5*tk,deltrs
!       PE=(PEREV+PE)*0.5
        pe=perev
        P(99)=PE
        IF(DELTRS.LE.EPS) THEN
           EXIT RUSS
         ELSE
           NITERR=NITERR+1
           IF(NITERR.LE.300) THEN
              CYCLE RUSS
            ELSE
              WRITE(6,"('*DOES NOT CONVERGE AFTER ',I4,' ITERATIONS')")   &
     &        niterr
            EXIT RUSS
           END IF
        END IF
      END DO RUSS
!
      RETURN
      END SUBROUTINE RUSSEL


!
!     ***********************************************************************
!

      subroutine bcdata(indtsu,ifpf,t,pfeq)
!     =====================================
!
!     evaluates either ooartition function or en equuilibrium constant
!     using the Barklem & Collet (2016 - BC) data
!
!     Input:  indtsu - index in the Tsuji-like indexing
!             ifpf switch for evaluted quantity:
!                  = 0 - evaluted the equilibrium constant
!                 != 0 - partition function
!             t - temperature
!     Output: pfeq - either the partition function or equil.const.
!
!     array  IBC is the BC indec corresponfing to the Tsuji-like index indtsu
!
      use accura
      use modelp,only : tlbc,pfbc,eqbc
      implicit real(dp) (a-h,o-z)

      integer  :: ibc(468)
!
      data ibc/                                                           &
     &   0,   1,   0,  38,  36, 102, 101,   4,   5,   6, 118,  37,        &
     &   0,   0,   0,  45,  43,   0,   0, 106, 104,   0, 123, 121,        &
     & 138, 140,  13, 231, 145, 146, 160,  41,  39,  48, 168,  46,        &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0, 190, 195,   0,        &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 254, 252,        &
     & 257,   0, 258,   0,   0,   0, 235,   0,  33, 131, 185, 249,        &
     &   0,   0,   0, 132, 186, 250,   0,   0,   0,  35,   0, 133,        &
     &   0, 228, 187, 251,   0,   0,   0,   0,  32,  97, 115,   0,        &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,        &
     &   0,   0,   0,   0,   0,   0,  98,   0,   0,   0,   0,   7,        &
     & 134,  40, 135, 188,   0, 136, 229, 189,   0, 253,   0,   0,        &
     &  42, 137,   0, 230,   0,   0,   0,   0,   0,   0,   0, 191,        &
     &   0, 255,   0,  44,   0,   0, 105,   0, 139,   0, 232, 192,        &
     &   0, 256,   0,   0,   0,   0,   0,  14, 107,   0,   0,   0,        &
     & 141,   0,   0,   0,   0,   0,   0,  47, 142, 194, 143, 233,        &
     &   0,   0,   0, 144,   0,   0, 234, 196,   0,   0,   0,   0,        &
     &   0,   0,   0,   0,   0,  50,   0, 147,   0,   0,   0,   0,        &
     &   0,  51, 148,   0,   0,   0, 197,   0, 260,  52,   0, 149,        &
     &   0,   0,   0, 261,   0,   0,  54, 150, 198,   0,  55, 151,        &
     & 237, 199, 262,  62, 206, 269,  63, 158,   0, 241, 207,   0,        &
     & 270,   0,   0,   0,   0, 159,   0,   0, 242, 208,   0, 271,        &
     &   0,   0, 128,   0,   0,   0,   0,   0,   0,  70,  72, 167,        &
     &   0, 245, 216,   0, 279,   0,   0,   0,   0,   0,   0,   0,        &
     &   0, 246, 217,   0,   0,   0,   0,   0,   0,   0,   0,   0,        &
     &   0,   0, 238, 153,  58, 202, 265, 161,   0,   2,   3,   8,        &
     &   9,  10,  11,  12,  15,  16,  17,  18,  19,  20,  21,  22,        &
     &  31,  34,  49,  53,  56,  57,  59,  60,  61,  64,  65,  66,        &
     &  67,  68,  69,  71,  73,  74,  75,  76,  77,  78,  79,  99,        &
     & 100, 103, 108, 109, 110, 111, 112, 116, 117, 119, 120, 124,        &
     & 125, 126, 127, 152, 154, 155, 156, 157, 162, 163, 164, 165,        &
     & 166, 169, 170, 171, 172, 173, 174, 175, 176, 177, 184, 193,        &
     & 200, 201, 203, 204, 205, 209, 210, 211, 212, 213, 214, 215,        &
     & 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 236, 239,        &
     & 240, 243, 244, 247, 248, 259, 263, 264, 266, 267, 268, 272,        &
     & 273, 274, 275, 276, 277, 278, 280, 281, 282, 283, 284, 285,        &
     & 286, 287, 288, 289, 290,   0,  23,  24,  25,  26,  27,  28,        &
     &  29,  30,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,        &
     &  90,  91,  92,  93,  94,  95,  96, 113, 114, 129, 130, 178,        &
     & 179, 180, 181, 182, 183,   0,   0,   0,   0,   0,   0,   0/
!
      data iread/1/
!
      if(iread.eq.1) then
!
!---  read data from BC tables in the first call
!
         mt=1000
         open(67,file='./data/BC.eq',status='old')
         open(68,file='./data/BC.pf',status='old')
         read(67,*) (tlbc(m),m=1,mt)
         read(68,*) (tlbc(m),m=1,mt)
         do i=1,291
            read(67,*)
            read(68,*)
            read(67,*) (eqbc(i,m),m=1,mt)
            read(68,*) (pfbc(i,m),m=1,mt)
         end do
         close(67)
         close(68)
         iread=0
      end if
!
      pfeq=0.
      if(ifpf.eq.0) pfeq=-15000.
      if(t.gt.1.e4) return
      if(indtsu.le.0.or.indtsu.gt.468) return
      in=ibc(indtsu)
      inbc=in
      if(in.eq.0) return

!
!     new table 7
!
      if(ifpf.eq.0) then
         call bceqco(in,t,bcel)
         pfeq=10.**bcel
         return
      end if
!
!     temperature interpolation in the BC tables
!
      iii=int(t)/10
      indt=min(iii,999)
      a1=1.- (log10(t)-tlbc(indt))/(tlbc(indt+1)-tlbc(indt))
!
      if(ifpf.eq.0) then
         pfeq=a1*eqbc(in,indt) + (1.-a1)*eqbc(in,indt+1) + 1.
       else
         pfeq=a1*pfbc(in,indt) + (1.-a1)*pfbc(in,indt+1)
!        pfeq=10.**pfeq
      end if
      pfeq=10.**pfeq
!
      return
      end subroutine bcdata

!
!     ***********************************************************************
!

      subroutine irwpf(jatom,ion,indmol,t,u)
!     ======================================
!
!     partition functions adter Irwin (1981), ApJS. 45, 621.
!     updated with the data of Barklem & Collet (2016)
!     set to the Irwin format by Y. Ossorio
!
!     Input: jatom - atomic number; if =0 - molecules
!            ion - ionization degree
!            indmol - index of a molecule in the new Tsuji-type
!                     indexing (from file tsuji.molec_bc2)
!            t - temperature
!     Output: u - partition function
!
!     array IRWIND(I) - the Irwin index corresponding to Tsuji
!           index I
!           if =0 - molecule I has no data in the Irwin table
!
      use accura
      use params
      implicit real(dp) (a-h,o-z)

      real(dp),save :: a(6,3,92),aa(6),am(6,500)
      real(dp)      :: spec(500)
      integer       :: irwind(478)
!
      data irwind/                                                        &
     &    0,   1,  28,   4,   2,   7,   6,   5,   8,  10,                 &
     &    9,   3,  18,  25,  53,  29,  43,   0,  17, 153,                 &
     &   52,  55, 167,  44,  45, 182,  74,  46,  11, 187,                 &
     &  201,  31,  27,  99, 209,  24,  22,  20,  21,  65,                 &
     &   35,  19,  54,  23,   0,  14,  58,   0,  32,  12,                 &
     &   47,  16,   0,  34,   0,   0,  30,   0,  13,  33,                 &
     &   61,  63, 292,  57,  59,  66, 272,   0,  94, 175,                 &
     &  226, 286,   0,   0,   0, 176, 227, 287,   0,   0,                 &
     &    0,  96,   0, 177,   0, 267, 228, 288,   0,   0,                 &
     &    0,   0,  93, 147, 162, 5*0,                                     &
     &    0,  50,   0,   0,   0,   0,  36,   0,  64,   0,                 &
     &    0,  48,   0,   0, 148,   0,   0,  26,  49,  70,                 &
     &  178,  97, 170, 229,   0, 180, 268, 230,   0, 289,                 &
     &    0,   0,  15, 181,   0, 269, 4*0,                                &
     &    0,   0,   0, 231,   0, 290,   0,  38,   0,   0,                 &
     &  152,  39,  40,   0,  41, 232,   0, 291,   0,   0,                 &
     &    0,   0,   0,  75, 154,   0,   0,   0, 183,   0,                 &
     &    0,   0,   0,   0,   0,  98, 184, 234, 185, 270,                 &
     &    0,   0,   0, 186,   0,   0, 271, 235,   0,   0,                 &
     &   62,   0,   0,   0,   0,   0,   0, 101,   0, 188,                 &
     &    0,   0,   0,   0,   0, 102, 189, 3*0,                           &
     &  236,   0, 294,  67,   0, 190,   0,   0,   0, 295,                 &
     &    0,   0, 104, 191, 237,   0, 105, 192, 274, 238,                 &
     &  296, 112, 245, 303, 113, 199,   0, 278, 246,   0,                 &
     &  304,   0,   0,   0,   0, 200,   0,   0, 279, 247,                 &
     &    0, 305,   0,   0, 172, 5*0,                                     &
     &    0, 120, 122, 208,   0, 282, 255,   0, 312,   0,                 &
     &  7*0,                               283, 256,   0,                 &
     & 10*0,                                                              &
     &  275, 194, 108, 241, 299, 202,   0,  68,  69,  71,                 &
     &   72,  73,  42,  37,  76,  77,  78,  79,  80,  81,                 &
     &   82,  83,  92,  95, 100, 103, 106, 107, 109, 110,                 &
     &  111, 114, 115, 116, 117, 118, 119, 121, 123, 124,                 &
     &  125, 126, 127, 128, 129, 149, 150, 151, 155, 156,                 &
     &  157, 158, 159, 163, 164, 165, 166, 168, 169, 170,                 &
     &  171, 193, 195, 196, 197, 198, 203, 204, 205, 206,                 &
     &  207, 210, 211, 212, 213, 214, 215, 216, 217, 218,                 &
     &  225, 233, 239, 240, 242, 243, 244, 248, 249, 250,                 &
     &  251, 252, 253, 254, 257, 258, 259, 260, 262, 262,                 &
     &  263, 264, 265, 266, 273, 276, 277, 280, 282, 284,                 &
     &  285, 293, 297, 298, 300, 301, 302, 306, 307, 308,                 &
     &  309, 310, 311,  60, 313, 314, 315, 316, 317, 318,                 &
     &  319, 320, 321, 322, 323, 324,  84,  85,  86,  87,                 &
     &   88,  89,  90,  91, 130, 131, 132, 133, 134, 135,                 &
     &  136, 137, 138, 139, 140, 141, 142, 143, 144, 145,                 &
     &  146, 160, 161, 173, 174, 210, 220, 221, 222, 223,                 &
     &  224,16*0,  56/
!
      data iread /0/
!
!     call old Irwin routine MPARTF if desired
!
      if(irwtab.eq.0) then
         call mpartf(jatom,ion,indmol,t,u)
         return
      end if
!
!     read data if first call:
!
      if(iread.ne.1) then
         if(irwtab.eq.0) then
            open(67,file= './data/irwin_orig.dat',status='old')
          else
            open(67,file= './data/irwin_bc.dat',status='old')
         end if
         read(67,*)
         read(67,*)
         atoms: do j=1,92
            ions: do i=1,3
               if(j.eq.1.and.i.eq.3) cycle ions
               sp=float(j)+float(i-1)/100.
               read(67,*) spc,aa
               do k=1,6
                  a(k,i,j)=aa(k)
               end do
            end do ions
         end do atoms
!
         read(67,*)
         read(67,*)
         read(67,*)
         molec: do i=1,324
            read(67,*,iostat=ios) spec(i),aa
            if(ios.ne.0) exit molec
            do j=1,6
               am(j,i)=aa(j)
            end do
         end do molec
         close(67)
         iread=1
      end if
!
!     evaluation of the partition function
!     stop if T is out of limits of Irwin's tables
!
      tl=log(t)
      u=0.
      if(t.lt.1000.) then
         if(indmol.gt.0) then
            call bcdata(indmol,1,t,ulog)
            u=ulog
!           write(*,*) 'pfmol',indmol,t,ulog
            return
          else
            tl=6.9077553
         end if
       else if(t.gt.16000.) then
         write(*,*) 'irwpf, T=',t
         stop 'partf; temp>16000 K'
      end if
!
!     atomic species
!
      if(jatom.gt.0.and.ion.gt.0) then
        ulog=    a(1,ion,jatom)+                                          &
     &       tl*(a(2,ion,jatom)+                                          &
     &       tl*(a(3,ion,jatom)+                                          &
     &       tl*(a(4,ion,jatom)+                                          &
     &       tl*(a(5,ion,jatom)+                                          &
     &       tl*(a(6,ion,jatom))))))
       if(jatom.eq.5.and.ion.eq.3) ulog=1.
        u=exp(ulog)
        return
      end if
!
!     molecular species
!
      if(indmol.gt.0) then
         indm=irwind(indmol)
         if(indm.le.0) return
         ulog=    am(1,indm)+                                             &
     &        tl*(am(2,indm)+                                             &
     &        tl*(am(3,indm)+                                             &
     &        tl*(am(4,indm)+                                             &
     &        tl*(am(5,indm)+                                             &
     &        tl*(am(6,indm))))))
        u=exp(ulog)
      end if
      return
      end subroutine irwpf

!
!     ***********************************************************************
!


      subroutine exopf(indmol,t,u)
!     ============================
!
!     oartition functions from EXOMOL for 32 molewcular species
!
      use accura
      use params
      use modelp, only : pf
      implicit real(dp) (a-h,o-z)

      integer, parameter :: nmol=32
      character(len=4) :: filpf(nmol)
      character(len=7) :: fil
      character(len=6) :: fil1
      character(len=1) :: fil0
      character(len=17):: fil5
      character(len=18):: fil6
      integer :: indtsu(nmol),ntemp(nmol)
!
      data filpf/                                                         &
     &  ' AlO','  C2','  CH','  CN','  CO',                               &
     &  '  CS',' CaH',' CaO',' CrH',' FeH',                               &
     &  '  H2',' HCl','  HF',' MgH',' MgO',                               &
     &  '  N2','  NH','  NO','  NS',' NaH',                               &
     &  '  OH','  PH','  SH',' SiH',' SiO',                               &
     &  ' SiS',' TiH',' TiO','  VO',                                      &
     &  ' H2O',' H2S',' CO2'/
      data ntemp/                                                         &
     &     9,  10,   8,   3,   9,   3,   3,   8,   3,  10,                &
     &    10,   5,   5,   3,   5,   9,   5,   5,   5,   5,                &
     &     5,   4,   5,   5,   9,   5,  48,   8,   8,  10,                &
     &     3,   5/
      data indtsu/                                                        &
     &   134,   8,   5,   7,   6,  20,  34, 179, 198, 214,                &
     &     2,  36,  33,  32, 126,   9,  12,  11,  23, 122,                &
     &     4, 148,  16,  17,  25,  28, 315,  29,  30,   3,                &
     &    57,  44/
      data iread /1/
!
      if(iread.eq.1) then
         do i=1,nmol
            ntemp(i)=ntemp(i)*1000
         end do
         ntemp(27)=ntemp(27)/10
         do i=1,nmol
            fil=filpf(i)//'.pf'
            fil1=fil(2:)
            fil0=fil1(:1)
            if(fil0.eq.' ') then
               fil5='data/EXOMOL/'//fil1(2:)
               open(unit=67,file=fil5,status='old')
             else
               fil6=fil1
               open(unit=67,file='data/EXOMOL/'//fil6,status='old')
            end if
            do j=1,ntemp(i)
               read(67,*) tt,pf(i,j)
            end do
            close(67)
         end do
         iread=0
      end if
!
      ie=0
      u=0.
      do i=1,nmol
         if(indtsu(i).eq.indmol) ie=i
      end do
      if(ie.eq.0) return
!
      tmax=float(ntemp(ie))
      if(t.le.tmax) then
         j=int(t)
         u=pf(ie,j)
       else
         call irwpf(0,0,indmol,tmax,umx)
         call irwpf(0,0,indmol,t,uirw)
         u=pf(ie,ntemp(ie))/umx*uirw
      end if
!
      return
      end subroutine exopf
!
!     ***********************************************************************
!


      subroutine equcon(t)
!     ====================
!
!     correction for equilibrium constant for moleculare ions to express
!     Kp from the B&C tables in a Tsuji-like form

      use accura
      use params, only: nmolec,ieqbc,amas
      use modelp, only: cmol
      use eospar
      implicit real(dp) (a-h,o-z),logical (l)

      real(dp), intent(in) :: t
      real(dp) :: eaf(99)
      integer  :: ineg(17)
      character(len=3) dyp(99)

      integer  :: ibc(468)
!
      data ibc/                                                           &
     &   0,   1,   0,  38,  36, 102, 101,   4,   5,   6, 118,  37,        &
     &   0,   0,   0,  45,  43,   0,   0, 106, 104,   0, 123, 121,        &
     & 138, 140,  13, 231, 145, 146, 160,  41,  39,  48, 168,  46,        &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0, 190, 195,   0,        &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 254, 252,        &
     & 257,   0, 258,   0,   0,   0, 235,   0,  33, 131, 185, 249,        &
     &   0,   0,   0, 132, 186, 250,   0,   0,   0,  35,   0, 133,        &
     &   0, 228, 187, 251,   0,   0,   0,   0,  32,  97, 115,   0,        &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,        &
     &   0,   0,   0,   0,   0,   0,  98,   0,   0,   0,   0,   7,        &
     & 134,  40, 135, 188,   0, 136, 229, 189,   0, 253,   0,   0,        &
     &  42, 137,   0, 230,   0,   0,   0,   0,   0,   0,   0, 191,        &
     &   0, 255,   0,  44,   0,   0, 105,   0, 139,   0, 232, 192,        &
     &   0, 256,   0,   0,   0,   0,   0,  14, 107,   0,   0,   0,        &
     & 141,   0,   0,   0,   0,   0,   0,  47, 142, 194, 143, 233,        &
     &   0,   0,   0, 144,   0,   0, 234, 196,   0,   0,   0,   0,        &
     &   0,   0,   0,   0,   0,  50,   0, 147,   0,   0,   0,   0,        &
     &   0,  51, 148,   0,   0,   0, 197,   0, 260,  52,   0, 149,        &
     &   0,   0,   0, 261,   0,   0,  54, 150, 198,   0,  55, 151,        &
     & 237, 199, 262,  62, 206, 269,  63, 158,   0, 241, 207,   0,        &
     & 270,   0,   0,   0,   0, 159,   0,   0, 242, 208,   0, 271,        &
     &   0,   0, 128,   0,   0,   0,   0,   0,   0,  70,  72, 167,        &
     &   0, 245, 216,   0, 279,   0,   0,   0,   0,   0,   0,   0,        &
     &   0, 246, 217,   0,   0,   0,   0,   0,   0,   0,   0,   0,        &
     &   0,   0, 238, 153,  58, 202, 265, 161,   0,   2,   3,   8,        &
     &   9,  10,  11,  12,  15,  16,  17,  18,  19,  20,  21,  22,        &
     &  31,  34,  49,  53,  56,  57,  59,  60,  61,  64,  65,  66,        &
     &  67,  68,  69,  71,  73,  74,  75,  76,  77,  78,  79,  99,        &
     & 100, 103, 108, 109, 110, 111, 112, 116, 117, 119, 120, 124,        &
     & 125, 126, 127, 152, 154, 155, 156, 157, 162, 163, 164, 165,        &
     & 166, 169, 170, 171, 172, 173, 174, 175, 176, 177, 184, 193,        &
     & 200, 201, 203, 204, 205, 209, 210, 211, 212, 213, 214, 215,        &
     & 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 236, 239,        &
     & 240, 243, 244, 247, 248, 259, 263, 264, 266, 267, 268, 272,        &
     & 273, 274, 275, 276, 277, 278, 280, 281, 282, 283, 284, 285,        &
     & 286, 287, 288, 289, 290,   0,  23,  24,  25,  26,  27,  28,        &
     &  29,  30,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,        &
     &  90,  91,  92,  93,  94,  95,  96, 113, 114, 129, 130, 178,        &
     & 179, 180, 181, 182, 183,   0,   0,   0,   0,   0,   0,   0/
!
      data ineg/1,0,0,0,0,2,0,3,4,0,0,0,0,5,0,6,7/

      data eaf /                                                          &
     &  0.7541950, -0.5000000,  0.6180490, -0.5000000,  0.2797230,        &
     &  1.2621226, -0.0700000,  1.4611136,  3.4011898, -1.2000000,        &
     &  0.5479260, -0.4000000,  0.4328300,  1.3895212,  0.7466070,        &
     &  2.0771042,  3.6127250, -1.0000000,  0.5014590,  0.0245500,        &
     &  0.1880000,  0.0755400,  0.5276600,  0.6758400, -0.5000000,        &
     &  0.1532360,  0.6622600,  1.1571600,  1.2357800, -0.6000000,        &
     &  0.3012000,  1.2326764,  0.8048000,  2.0206047,  3.3635880,        &
     & -1.0000000,  0.4859160,  0.0520600,  0.3070000,  0.4332800,        &
     &  0.9174000,  0.7473000,  0.5500000,  1.0463800,  1.1428900,        &
     &  0.5621400,  1.3044700, -0.7000000,  0.3839200,  1.1120700,        &
     &  1.0474010,  1.9708750,  3.0590520, -0.8000000,  0.4716300,        &
     &  0.1446200,  0.5500000,  0.5700000,  0.9620000,  1.9160000,        &
     &  0.1290000,  0.1620000,  0.1160000,  0.1370000,  1.1650000,        &
     &  0.3520000,  0.3380000,  0.3120000,  1.0290000, -0.0200000,        &
     &  0.2388000,  0.1780000,  0.3230000,  0.8162600,  0.0603960,        &
     &  1.0778000,  1.5643600,  2.1251000,  2.3086100, -0.5000000,        &
     &  0.3770000,  0.3567210,  0.9423620,  1.4000000,  2.4200000,        &
     & -0.7000000,  0.4860000,  0.1000000,  0.3500000,  1.1700000,        &
     &  0.5500000,  0.5300000,  0.4800000, -0.5000000,  0.1000000,        &
     &  0.2800000, -1.7200000, -1.0100000, -0.3000000/

      DATA DYP/'  H',' He',' Li',' Be','  B','  C',
     &         '  N','  O','  F',' Ne',' Na',' Mg',                       &
     &         ' Al',' Si','  P','  S',' Cl',' Ar',                       &
     &         ' K ',' Ca',' Sc',' Ti',' V ',' Cr',                       &
     &         ' Mn',' Fe',' Co',' Ni',' Cu',' Zn',                       &
     &         ' Ga',' Ge',' As',' Se',' Br',' Kr',                       &
     &         ' Rb',' Sr',' Y ',' Zr',' Nb',' Mo',                       &
     &         ' Tc',' Ru',' Rh',' Pd',' Ag',' Cd',                       &
     &         ' In',' Sn',' Sb',' Te',' I ',' Xe',                       &
     &         ' Cs',' Ba',' La',' Ce',' Pr',' Nd',                       &
     &         ' Pm',' Sm',' Eu',' Gd',' Tb',' Dy',                       &
     &         ' Ho',' Er',' Tm',' Yb',' Lu',' Hf',                       &
     &         ' Ta',' W ',' Re',' Os',' Ir',' Pt',                       &
     &         ' Au',' Hg',' Tl',' Pb',' Bi',' Po',                       &
     &         ' At',' Rn',' Fr',' Ra',' Ac',' Th',                       &
     &         ' Pa',' U ',' Np',' Pu',' Am',' Cm',                       &
     &         ' Bk',' Cf',' Es'/
!     

!
!     tk=t*1.38054e-16         ! = k T
!     betae2=4.8298e15*t**1.5  ! = 2 * (2 pi m_e k/h**2)**1.5 * T**1.5
      conl=0.17602759          ! log10(2 * (2 pi m_e k/h**2)**1.5 * k)
      TH=5040./T
 
      pdefl=0.


      write(35,"('   TS   BC              ion  neu           ',           &
     & ' eaf(ion)  eaf(neu)    U(ion+)   U(ion)     corr    ',            &
     & 'K_IH+c   K_BC+c    K_YO   K_YO-(K_BC+c)'/)")
      write(36,"('   TS   BC              ion  neu           ',           &
     & ' eaf(ion)  eaf(neu)    U(ion-)   U(ion)     corr    ',            &
     & 'K_IH+c   K_BC+c    K_YO   K_YO-(K_BC+c)'/)")


      molloop: do j=1,nmolec
         APLOGJ=C(J,5)
         DO K=1,4
            KM5=5-K
            APLOGJ=APLOGJ*TH + C(J,KM5)
         END DO
         apts=aplogj
         apmlog(j)=aplogj
         inbc=ibc(j)
         if(ieqbc.le.0.or.j.gt.468) cycle molloop
         if(ibc(j).eq.0) cycle molloop

         inbc=ibc(j)
         nat=0

         call bceqco(inbc,t,aplogj)

!        write(*,"(2i4,a8,f8.1,f10.3)") j,inbc,cmol(j),t,aplogj

         lneg2=mmax(j).eq.3.and.nelem(3,j).eq.99.and.nato(3,j).eq.1.      &
     &         and.nato(1,j).eq.1.and.nato(2,j).eq.1
         lneg1=mmax(j).eq.2.and.nelem(2,j).eq.99.and.nato(2,j).eq.1.      &
     &      and.nato(1,j).eq.2

         lpos2=mmax(j).eq.3.and.nelem(3,j).eq.99.and.nato(3,j).eq.-1
         lpos1=mmax(j).eq.2.and.nelem(2,j).eq.99.and.nato(2,j).eq.-1
!
!        positive ions
!        -------------
!
         if(lpos2.or.lpos1) then
!
!           positive diatomics  (XY+)
!
            if(lpos2) then
               if(nelem(1,j).eq.1.or.nelem(2,j).eq.1) then
!
!                 ionized hydrides
!
                  if(nelem(1,j).eq.1) nat=nelem(2,j)
                  if(nelem(2,j).eq.1) nat=nelem(1,j)
                  if(nelem(1,j).eq.2.or.nelem(2,j).eq.2) nat=1
                  if(nelem(1,j).eq.9.or.nelem(2,j).eq.9) nat=1
                  if(nelem(1,j).eq.10.or.nelem(2,j).eq.10) nat=1
               end if
!
               if(nelem(1,j).eq.8.or.nelem(2,j).eq.8) then
!
!                 ionized oxides
!
                  if(nelem(1,j).eq.8) nat=nelem(2,j)
                  if(nelem(2,j).eq.8) nat=nelem(1,j)
               end if
!
!              OH+,CN+,NS+,NO+
!
               if(nelem(1,j).eq.8.and.nelem(2,j).eq.1) nat=8
               if(nelem(1,j).eq.6.and.nelem(2,j).eq.7) nat=6
               if(nelem(1,j).eq.7.and.nelem(2,j).eq.16) nat=16
               if(nelem(1,j).eq.7.and.nelem(2,j).eq.8) nat=8

               nio=nat
               if(nio.eq.nelem(1,j)) neu=nelem(2,j)
               if(nio.eq.nelem(2,j)) neu=nelem(1,j)
               amasr=amas(nio)*amas(neu)/(amas(nio)+amas(neu))
            end if
!
!           ionized homonuclear molecules (X2+)
!
            if(lpos1) then
               nat=nelem(1,j)
               nio=nat
               neu=nat
               amasr=amas(nat)*0.5
            end if

            call irwpf(nat,1,0,t,uato)
            call irwpf(nat,2,0,t,uion)
            uato=log10(uato)
            uion=log10(uion)
            corr=conl+uato-uion-2.5*log10(t)+5039.9*xip(nat)/t
            apmlog(j)=aplogj+corr
!           write(*,"(i4,3f11.3)") j,aplogj,corr,apmlog(j)
            write(35,"(2i5,3x,a8,2i5,2x,a3,'+',a3,2x,2f9.3,2x,
     &      2f9.3,2x,f9.3,2x,5f9.3)")                                     &
     &      j,inbc,cmol(j),nio,neu,dyp(nio),dyp(neu),eaf(nio),eaf(neu),   &
     &      uion,uato,                                                    &
     &      corr,pdefl+corr,aplogj+corr,apts,apts-aplogj-corr

            cycle molloop
         end if
!
!        negative ions
!        ------------
!
         if(lneg2.or.lneg1) then
!
!           negative ion of a diatomic molecule (XY-)
!
            if(lneg2) then
               if(eaf(nelem(1,j)).ge.eaf(nelem(2,j))) then
                  nat=nelem(1,j)
                else
                  nat=nelem(2,j)
               end if
               nio=nat
               if(nio.eq.nelem(1,j)) neu=nelem(2,j)
               if(nio.eq.nelem(2,j)) neu=nelem(1,j)
               amasr=amas(nio)*amas(neu)/(amas(nio)+amas(neu))
            end if
!
!           negative ion of a homonuclear molecule (X2-)
!
            if(lneg1) then
               nat=nelem(1,j)
               amasr=0.5*amas(nat)
               nio=nat
               neu=nat
            end if

            uion=1.
            if(ineg(nio).gt.0) call bcnega(ineg(nio),t,uion)
            call irwpf(nat,1,0,t,uato)
            uato=log10(uato)
            uion=log10(uion)
            corr=-conl+uato-uion+2.5*log10(t)-5039.9*eaf(nat)/t
            apmlog(j)=aplogj+corr
!           write(*,"(i4,3f11.3)") j,aplogj,corr,apmlog(j)

            write(36,"(2i5,3x,a8,2i5,2x,a3,'-',a3,2x,2f9.3,2x,
     &      2f9.3,2x,f9.3,2x,5f9.3)")                                     &
     &      j,inbc,cmol(j),nio,neu,dyp(nio),dyp(neu),eaf(nio),eaf(neu),   &
     &      uion,uato,                                                    &
     &      corr,pdefl+corr,aplogj+corr,apts,apts-aplogj-corr

         end if
      end do molloop

      return
      end subroutine equcon
!
!     ***********************************************************************
!

      subroutine eospri
!     =================
!
!     Outprint of Equation of State parameters
!
      use accura
      use params
      use modelp
      use eospar
      implicit real(dp) (a-h,o-z)

      integer  :: neleme(38),insm(20)
      real(dp) :: amh2(5),xml(20)
      data neleme/ 1, 2, 3, 4, 5, 6, 7, 8, 9,                             &
     &            11,12,13,14,15,16,17,19,20,                             &
     &            21,22,23,24,25,26,28,29,32,                             &
     &            35,37,38,39,40,41,53,56,57,58,60/
      data amh2/1.13390E+01,-2.97499E+00,4.10842E-02,-3.58550E-03,        &
     &          1.31844E-04/
      data insm/2,3,4,5,6,7,8,12,17,25,29,30,32,34,122,126,134,           &
     &          179,198,214/
      data init/1/
!
!     id=idstd
      istp=1
      nd=1
      if(ifeos.lt.0) istp=-ifeos
      tmolim=10000.
      ifmol=1
!
      do id=1,nd,istp
         t=temp(id)
         ane=elec(id)
         rho=dens(id)
         ann = dens(id)/wmm(id)+elec(id)
!
!!!!!    if(ifmol.eq.0.or.t.gt.tmolim) then
!!!!!       it=0
!!!!!       itera: do
!!!!!          ann0=ann
!!!!!          it=it+1
!!!!!          call eldens(id,t,ann,ane,anh,anp)
!!!!!          anmol(1,id)=anhmi
!!!!!          anmol(2,id)=ahmol
!!!!!          anato(1,id)=anh
!!!!!          anion(1,id)=anp
!!!!!          hpop=dens(id)/wmy(id)/hmass
!!!!!          do i=1,nmetal
!!!!!             j=neleme(i)
!!!!!             anato(j,id)=anato(j,id)*hpop
!!!!!             anion(j,id)=anion(j,id)*hpop
!!!!!             if(j.ge.2.and.j.le.30) anion2(j,id)=anion2(j,id)*hpop
!!!!!          end do
!!!!!          anato(1,id)=anh
!!!!!          anion(1,id)=anp
!!!!!          wmm(id)=wmy(id)/(ytot(id)-anmol(2,id)/hpop)*hmass
!!!!!          ann=dens(id)/wmm(id)+ane
!!!!!          if((ann-ann0)/ann0.le.1.e-5.or.it.gt.20) exit
!!!!!       end do itera
!!!!!    end if
!
         nmetae=38
!!       write(6,"(' **** DEPTH ID',i4,'  T =',f10.1/)") id,temp(id)
         write(*,*) ''
         write(*,*) 'atomic number densities and partition functions'
         write(*,*) ''
         atot=0.
         do i=1,nmetae
            j=neleme(i)
            if(j.le.28)                                                   &
     &         write(6,"(i4,a3,3x,1p2e12.4)")                             &
     &         j,typat(j),anato(j,id),pfato(j,id)
            atot=atot+anato(j,id)
         end do
         write(*,*) ''
         write(*,*) 'ionic number densities and partition functions'
         write(*,*) ''
         ctot=0.
         do i=1,nmetae
            j=neleme(i)
            if(j.le.28)                                                   &
     &         write(6,"(i4,a3,'+',2x,1p2e12.4)")                         &
     &         j,typat(j),anion(j,id),pfion(j,id)
            atot=atot+anion(j,id)
            ctot=ctot+anion(j,id)
         end do
!
         if(ifmol.gt.0.and.t.le.tmolim) then
             write(6,"(/ 'Molecular number densities ',                   &
     &          'and partition functions'/)")
             do i=1,nmolec
                if(anmol(i,id).gt.ann*1.e-15)                             &
     &             write(6,"(i4,1x,A8,1x,1pe12.4,1x,e12.4)")              &
     &             i, cmol(i), anmol(i,id), pfmol(i,id)
                atot=atot+anmol(i,id)
             end do
         end if
!
         ahmi=1.0353e-16/t/sqrt(t)*exp(8762.9/t)*                         &
     &        anato(1,id)*ane
!
!        original B&C H2+
!
         APLOGJ=amh2(5)
         te=5040./t
         DO K=1,4
            KM5=5-K
            APLOGJ=APLOGJ*TE + amh2(KM5)
         END DO
         tk=1.38054e-16*t
         ph2=-aplogj+log10(anato(1,id)*anion(1,id))+2.*log10(tk)
         anh2b=(10.**ph2)/tk

         htot=anato(1,id)+anion(1,id)+anmol(1,id)+                        &
     &         2.*(anmol(2,id)+anmol(3,id))+anmol(4,id)+anmol(5,id)+      &
     &         anmol(12,id)+2.*anmol(13,id)+anmol(14,id)+                 &
     &         anmol(15,id)+                                              &
     &         anmol(16,id)+anmol(17,id)+anmol(32,id)+anmol(34,id)+       &
     &         4.*anmol(37,id)+2.*anmol(38,id)+3.*anmol(39,id)+           &
     &         2.*anmol(40,id)+3.*anmol(41,id)+2.*anmol(57,id)+           &
     &         anmol(118,id)+anmol(133,id)+                               &
     &         2.*anmol(140,id)+3.*anmol(141,id)+4.*anmol(142,id)+        &
     &         anmol(148,id)+2.*anmol(149,id)+anmol(222,id)
         ahe= (anato(2,id)+anion(2,id)+anion2(2,id))/htot
         aca= (anato(6,id)+anion(6,id)+anion2(6,id))/htot
         acm= (anmol(5,id)+anmol(6,id)+                                   &
     &        anmol(7,id)+2.*(anmol(8,id)+2.*anmol(13,id))+               &
     &        anmol(14,id)+2.*anmol(15,id)+anmol(20,id)+                  &
     &        anmol(37,id)+anmol(38,id)+anmol(39,id)+                     &
     &        anmol(44,id)+anmol(118,id)+anmol(119,id)+                   &
     &        anmol(437,id)+anmol(453,id)                                 &
     &        )/htot
         ana= (anato(7,id)+anion(7,id)+anion2(7,id))/htot
         anm= (anmol(7,id)+2.*anmol(9,id)+anmol(11,id)+                   &
     &        anmol(12,id)+anmol(14,id)+anmol(23,id)+                     &
     &        anmol(24,id)+anmol(40,id)+anmol(41,id)+                     &
     &        anmol(109,id)+anmol(152,id)+anmol(347,id)+                  &
     &        anmol(438,id)+anmol(452,id)+anmol(454,id)                   &
     &        )/htot
         aoa= (anato(8,id)+anion(8,id)+anion2(8,id))/htot
         aom= (anmol(3,id)+anmol(4,id)+                                   &
     &        anmol(6,id)+2.*anmol(10,id)+anmol(11,id)+anmol(25,id)+      &
     &        anmol(26,id)+anmol(29,id)+anmol(30,id)+anmol(31,id)+        &
     &        anmol(35,id)+2.*anmol(44,id)+anmol(49,id)+anmol(51,id)+     &
     &        anmol(54,id)+2.*anmol(56,id)+anmol(65,id)+                  &
     &        2.*anmol(66,id)+anmol(84,id)+anmol(109,id)+                 &
     &        anmol(113,id)+anmol(115,id)+anmol(118,id)+                  &
     &        anmol(119,id)+anmol(126,id)+anmol(134,id)+                  &
     &        anmol(153,id)+anmol(179,id)+anmol(184,id)+                  &
     &        2.*anmol(185,id)+anmol(200,id)+anmol(216,id)+               &
     &        anmol(221,id)+2.*anmol(247,id)+anmol(292,id)+               &
     &        anmol(439,id)+anmol(453,id)+anmol(454,id)                   &
     &        )/htot
         ac=aca+acm
         an=ana+anm
         ao=aoa+aom
         write(6,"(/'EOS useful quantities - summary'//                   &
     &      'T,rho       ',f13.2,1pe13.5/                                 &
     &      'N           ',1p2e13.5/                                      &
     &      'n_e         ',1p2e13.5/                                      &
     &      'H,H+,H-,H2  ',1p4e13.5/                                      &
     &      'H2-,H2+,H2+b',1p3e13.5/                                      &
     &      'Htot        ',1pe13.5/                                       &
     &      'H-          ',1p3e13.5/                                      &
     &      'C,C+,CO,CH4 ',1p4e13.5/                                      &
     &      'N,N+,N2,NH3 ',1p4e13.5/                                      &
     &      'O,O+,H2O,CO ',1p4e13.5/                                      &
     &      'He/H        ',1p2e13.5/                                      &
     &      'C/H         ',1p2e13.5/                                      &
     &      'N/H         ',1p2e13.5/                                      &
     &      'O/H         ',1p2e13.5/)")                                   &
     &       t,dens(id),ann,atot+ane,ane,ctot-anmol(1,id),                &
     &       anato(1,id),anion(1,id),                                     &
     &       anmol(1,id),anmol(2,id),                                     &
     &       anmol(312,id),anmol(426,id),anh2b,                           &
     &       htot,                                                        &
     &       anmol(1,id),ahmi,anmol(1,id)/ahmi,                           &
     &       anato(6,id),anion(6,id),anmol(6,id),anmol(37,id),            &
     &       anato(7,id),anion(7,id),anmol(9,id),anmol(41,id),            &
     &       anato(8,id),anion(8,id),anmol(3,id),anmol(6,id),             &
     &       ahe,ahe/abndd(2,id),                                         &
     &       ac,ac/abndd(6,id),                                           &
     &       an,an/abndd(7,id),                                           &
     &       ao,ao/abndd(8,id)
         act=ac*htot
         ant=an*htot
         aot=ao*htot
!
      if(init.eq.1) then
         write(52,"('    T      rho     w_mol    Ne/Ntot  N(Htot)    '    &
     &     'n(H)   n(H2)',6x,                                             &
     &     'a(He)   a(C)    a(N)    a(O)   molfr(C) molfr(N) ',           &
     &     'molfr(O)'/)")

         write(51,"('    T      rho     w_mol      N        Ne',          &
     &     '     N(Htot)   ',                                             &
     &     'N(H)    N(H+)    N(H-)   N(H2)    N(H2-)    N(H2+)'/)")
         write(53,"(' log10(N/U)'/'   T     rho   ',20a6/)")              &
     &      (cmol(insm(i)),i=1,20)
         write(54,"(' log10[N/n(H)]'/'   T     rho   ',20a6/)")           &
     &      (cmol(insm(i)),i=1,20)
!
          init=0
       end if
!
       write(52,"(f8.1,1pe9.2,0pf8.5,1x,1p4e9.2,1x,0p4f8.5,1x,1p3e9.2,    &
     &     1x,3e9.2,1x,3e9.2)") t,dens(id),wmm(id)/hmass,ane/ann,         &
     &     htot,anato(1,id),2.*anmol(2,id),                               &
     &     ahe/abndd(2,id),ac/abndd(6,id),an/abndd(7,id),ao/abndd(8,id),  &
     &     acm/ac,anm/an,aom/ao
!
       write(51,"(f8.1,1pe9.2,0pf8.5,1x,1p10e9.2)")                       &
     &     t,dens(id),wmm(id)/hmass,ann,ane,htot,                         &
     &     anato(1,id),anion(1,id),anmol(1,id),anmol(2,id),               &
     &     anmol(312,id),anmol(426,id)
!
      if(ifmol.gt.0.and.t.le.tmolim) then
         do i=1,20
            im=insm(i)
            xml(i)=log10(anmol(im,id)/pfmol(im,id))
         end do
         write(53,"(2f6.1,1x,20f6.1)")                                    &
     &      t,log10(dens(id)),(xml(i),i=1,20)
         do i=1,20
            im=insm(i)
            xml(i)=log10(anmol(im,id)/htot)
         end do
         write(54,"(2f6.1,1x,20f6.1)")                                    &
     &      t,log10(dens(id)),(xml(i),i=1,20)
      end if
!
      end do

      return
      end subroutine eospri

!
!     ***********************************************************************
!


      SUBROUTINE STATE0
!     =================
!
!     Initialization of the basic parameters for the Saha equation
!
      use accura
      use params
      implicit real(dp) (a-h,o-z),logical (l)

      real(dp), parameter :: enhe1=24.5799,enhe2=54.3999
      character(len=4)    :: DYP(MATOM)
!     character(len=80)   :: dum
      REAL(DP)            :: D(3,MATOM),XI(8,MATOM),                      &
     &                       ABUN0(MATOM),ABUN1(MATOM)
!
      DATA DYP/' H  ',' He ',' Li ',' Be ',' B  ',' C  ',                 &
     &         ' N  ',' O  ',' F  ',' Ne ',' Na ',' Mg ',                 &
     &         ' Al ',' Si ',' P  ',' S  ',' Cl ',' Ar ',                 &
     &         ' K  ',' Ca ',' Sc ',' Ti ',' V  ',' Cr ',                 &
     &         ' Mn ',' Fe ',' Co ',' Ni ',' Cu ',' Zn ',                 &
     &         ' Ga ',' Ge ',' As ',' Se ',' Br ',' Kr ',                 &
     &         ' Rb ',' Sr ',' Y  ',' Zr ',' Nb ',' Mo ',                 &
     &         ' Tc ',' Ru ',' Rh ',' Pd ',' Ag ',' Cd ',                 &
     &         ' In ',' Sn ',' Sb ',' Te ',' I  ',' Xe ',                 &
     &         ' Cs ',' Ba ',' La ',' Ce ',' Pr ',' Nd ',                 &
     &         ' Pm ',' Sm ',' Eu ',' Gd ',' Tb ',' Dy ',                 &
     &         ' Ho ',' Er ',' Tm ',' Yb ',' Lu ',' Hf ',                 &
     &         ' Ta ',' W  ',' Re ',' Os ',' Ir ',' Pt ',                 &
     &         ' Au ',' Hg ',' Tl ',' Pb ',' Bi ',' Po ',                 &
     &         ' At ',' Rn ',' Fr ',' Ra ',' Ac ',' Th ',                 &
     &         ' Pa ',' U  ',' Np ',' Pu ',' Am ',' Cm ',                 &
     &         ' Bk ',' Cf ',' Es '/
!
!    Standard atomic constants for first 99 species
!      Abundances for the first 30 from Grevesse & Sauval,
!         (1998, Space Sci. Rev. 85, 161)
!
!            Element Atomic  Solar    Std.
!                    weight abundance highest
!
!                                     ionization stage
      DATA D/ 1.008, 1.e0, 2.,                                            &
     &        4.003, 1.00e-1, 3.,                                         &
     &        6.941, 1.26e-11, 3.,                                        &
     &        9.012, 2.51e-11, 3.,                                        &
     &       10.810, 5.0e-10, 4.,                                         &
     &       12.011, 3.31e-4, 5.,                                         &
     &       14.007, 8.32e-5, 5.,                                         &
     &       16.000, 6.76e-4, 5.,                                         &
     &       18.918, 3.16e-8, 4.,                                         &
     &       20.179, 1.20e-4, 4.,                                         &
     &       22.990, 2.14e-6, 4.,                                         &
     &       24.305, 3.80e-5, 4.,                                         &
     &       26.982, 2.95e-6, 4.,                                         &
     &       28.086, 3.55e-5, 5.,                                         &
     &       30.974, 2.82e-7, 5.,                                         &
     &       32.060, 2.14e-5, 5.,                                         &
     &       35.453, 3.16e-7, 5.,                                         &
     &       39.948, 2.52e-6, 5.,                                         &
     &       39.098, 1.32e-7, 5.,                                         &
     &       40.080, 2.29e-6, 5.,                                         &
     &       44.956, 1.48e-9, 5.,                                         &
     &       47.900, 1.05e-7, 5.,                                         &
     &       50.941, 1.00e-8, 5.,                                         &
     &       51.996, 4.68e-7, 5.,                                         &
     &       54.938, 2.45e-7, 5.,                                         &
     &       55.847, 3.16e-5, 5.,                                         &
     &       58.933, 8.32e-8, 5.,                                         &
     &       58.700, 1.78e-6, 5.,                                         &
     &       63.546, 1.62e-8, 5.,                                         &
     &       65.380, 3.98e-8, 5.,                                         &
     &       69.72 ,   1.34896324e-09  ,  3.,                             &
     &       72.60 ,   4.26579633e-09  ,  3.,                             &
     &       74.92 ,   2.34422821e-10  ,  3.,                             &
     &       78.96 ,   2.23872066e-09  ,  3.,                             &
     &       79.91 ,   4.26579633e-10  ,  3.,                             &
     &       83.80 ,   1.69824373e-09  ,  3.,                             &
     &       85.48 ,   2.51188699e-10  ,  3.,                             &
     &       87.63 ,   8.51138173e-10  ,  3.,                             &
     &       88.91 ,   1.65958702e-10  ,  3.,                             &
     &       91.22 ,   4.07380181e-10  ,  3.,                             &
     &       92.91 ,   2.51188630e-11  ,  3.,                             &
     &       95.95 ,   9.12010923e-11  ,  3.,                             &
     &       99.00 ,   1.00000000e-24  ,  3.,                             &
     &       101.1 ,   6.60693531e-11  ,  3.,                             &
     &       102.9 ,   1.23026887e-11  ,  3.,                             &
     &       106.4 ,   5.01187291e-11  ,  3.,                             &
     &       107.9 ,   1.73780087e-11  ,  3.,                             &
     &       112.4 ,   5.75439927e-11  ,  3.,                             &
     &       114.8 ,   6.60693440e-12  ,  3.,                             &
     &       118.7 ,   1.38038460e-10  ,  3.,                             &
     &       121.8 ,   1.09647810e-11  ,  3.,                             &
     &       127.6 ,   1.73780087e-10  ,  3.,                             &
     &       126.9 ,   3.23593651e-11  ,  3.,                             &
     &       131.3 ,   1.69824373e-10  ,  3.,                             &
     &       132.9 ,   1.31825676e-11  ,  3.,                             &
     &       137.4 ,   1.62181025e-10  ,  3.,                             &
     &       138.9 ,   1.58489337e-11  ,  3.,                             &
     &       140.1 ,   4.07380293e-11  ,  3.,                             &
     &       140.9 ,   6.02559549e-12  ,  3.,                             &
     &       144.3 ,   2.95120943e-11  ,  3.,                             &
     &       147.0 ,   1.00000000e-24  ,  3.,                             &
     &       150.4 ,   9.33254366e-12  ,  3.,                             &
     &       152.0 ,   3.46736869e-12  ,  3.,                             &
     &       157.3 ,   1.17489770e-11  ,  3.,                             &
     &       158.9 ,   2.13796216e-12  ,  3.,                             &
     &       162.5 ,   1.41253747e-11  ,  3.,                             &
     &       164.9 ,   3.16227767e-12  ,  3.,                             &
     &       167.3 ,   8.91250917e-12  ,  3.,                             &
     &       168.9 ,   1.34896287e-12  ,  3.,                             &
     &       173.0 ,   8.91250917e-12  ,  3.,                             &
     &       175.0 ,   1.31825674e-12  ,  3.,                             &
     &       178.5 ,   5.37031822e-12  ,  3.,                             &
     &       181.0 ,   1.34896287e-12  ,  3.,                             &
     &       183.9 ,   4.78630102e-12  ,  3.,                             &
     &       186.3 ,   1.86208719e-12  ,  3.,                             &
     &       190.2 ,   2.39883290e-11  ,  3.,                             &
     &       192.2 ,   2.34422885e-11  ,  3.,                             &
     &       195.1 ,   4.78630036e-11  ,  3.,                             &
     &       197.0 ,   6.76082952e-12  ,  3.,                             &
     &       200.6 ,   1.23026887e-11  ,  3.,                             &
     &       204.4 ,   6.60693440e-12  ,  3.,                             &
     &       207.2 ,   1.12201834e-10  ,  3.,                             &
     &       209.0 ,   5.12861361e-12  ,  3.,                             &
     &       210.0 ,   1.00000000e-24  ,  3.,                             &
     &       211.0 ,   1.00000000e-24  ,  3.,                             &
     &       222.0 ,   1.00000000e-24  ,  3.,                             &
     &       223.0 ,   1.00000000e-24  ,  3.,                             &
     &       226.1 ,   1.00000000e-24  ,  3.,                             &
     &       227.1 ,   1.00000000e-24  ,  3.,                             &
     &       232.0 ,   1.20226443e-12  ,  3.,                             &
     &       231.0 ,   1.00000000e-24  ,  3.,                             &
     &       238.0 ,   3.23593651e-13  ,  3.,                             &
     &       237.0 ,   1.00000000e-24  ,  3.,                             &
     &       244.0 ,   1.00000000e-24  ,  3.,                             &
     &       243.0 ,   1.00000000e-24  ,  3.,                             &
     &       247.0 ,   1.00000000e-24  ,  3.,                             &
     &       247.0 ,   1.00000000e-24  ,  3.,                             &
     &       251.0 ,   1.00000000e-24  ,  3.,                             &
     &       254.0 ,   1.00000000e-24  ,  3./
!
      data abun0 /                                                        &
     &  12.00,10.93, 1.05, 1.38, 2.70, 8.39, 7.78, 8.66, 4.56, 7.84,      &
     &   6.17, 7.53, 6.37, 7.51, 5.36, 7.14, 5.50, 6.18, 5.08, 6.31,      &
     &   3.05, 4.90, 4.00, 5.64, 5.39, 7.45, 4.92, 6.23, 4.21, 4.60,      &
     &   2.88, 3.58, 2.29, 3.33, 2.56, 3.28, 2.60, 2.92, 2.21, 2.59,      &
     &   1.42, 1.92,-9.99, 1.84, 1.12, 1.69, 0.94, 1.77, 1.60, 2.00,      &
     &   1.00, 2.19, 1.51, 2.27, 1.07, 2.17, 1.13, 1.58, 0.71, 1.45,      &
     &  -9.99, 1.01, 0.52, 1.12, 0.28, 1.14, 0.51, 0.93, 0.00, 1.08,      &
     &   0.06, 0.88,-0.17, 1.11, 0.23, 1.45, 1.38, 1.64, 1.01, 1.13,      &
     &   0.90, 2.00, 0.65,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99, 0.06,      &
     &  -9.99,-0.52,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99/
!
      data abun1 /                                                        &
     &  12.00,10.93, 3.26, 1.38, 2.79, 8.43, 7.83, 8.69, 4.56, 7.93,      &
     &   6.24, 7.60, 6.45, 7.51, 5.41, 7.12, 5.50, 6.40, 5.08, 6.34,      &
     &   3.15, 4.95, 3.93, 5.64, 5.43, 7.50, 4.99, 6.22, 4.19, 4.56,      &
     &   3.04, 3.65, 2.30, 3.34, 2.54, 3.25, 2.36, 2.87, 2.21, 2.58,      &
     &   1.46, 1.88,-9.99, 1.75, 1.06, 1.65, 1.20, 1.71, 0.76, 2.04,      &
     &   1.01, 2.18, 1.55, 2.24, 1.08, 2.18, 1.10, 1.58, 0.72, 1.42,      &
     &  -9.99, 0.96, 0.52, 1.07, 0.30, 1.10, 0.48, 0.92, 0.10, 0.92,      &
     &   0.10, 0.85,-0.12, 0.65, 0.26, 1.40, 1.38, 1.62, 0.80, 1.17,      &
     &   0.77, 2.04, 0.65,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99, 0.06,      &
     &  -9.99,-0.54,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99/
!
!
!     Ionization potentials for first 99 species:
!
!     Element Ionization potentials (eV)
!              I     II      III     IV       V     VI     VII    VIII
!
      DATA XI/                                                            &
     &       13.595,  0.   ,  0.   ,  0.   ,  0.  ,  0.  ,  0.  ,  0.  ,  &
     &       24.580, 54.400,  0.   ,  0.   ,  0.  ,  0.  ,  0.  ,  0.  ,  &
     &        5.392, 75.619,122.451,  0.   ,  0.  ,  0.  ,  0.  ,  0.  ,  &
     &        9.322, 18.206,153.850,217.713,  0.  ,  0.  ,  0.  ,  0.  ,  &
     &        8.296, 25.149, 37.920,259.298,340.22,  0.  ,  0.  ,  0.  ,  &
     &       11.264, 24.376, 47.864, 64.476,391.99,489.98,  0.  ,  0.  ,  &
     &       14.530, 29.593, 47.426, 77.450, 97.86,551.93,667.03,  0.  ,  &
     &       13.614, 35.108, 54.886, 77.394,113.87,138.08,739.11,871.39,  &
     &       17.418, 34.980, 62.646, 87.140,114.21,157.12,185.14,953.6 ,  &
     &       21.559, 41.070, 63.500, 97.020,126.30,157.91,207.21,239.0 ,  &
     &        5.138, 47.290, 71.650, 98.880,138.37,172.09,208.44,264.16,  &
     &        7.664, 15.030, 80.120,102.290,141.23,186.49,224.9 ,265.96,  &
     &        5.984, 18.823, 28.440,119.960,153.77,190.42,241.38,284.53,  &
     &        8.151, 16.350, 33.460, 45.140,166.73,205.11,246.41,303.07,  &
     &       10.484, 19.720, 30.156, 51.354, 65.01,220.41,263.31,309.26,  &
     &       10.357, 23.400, 35.000, 47.290, 72.50, 88.03,280.99,328.8 ,  &
     &       12.970, 23.800, 39.900, 53.500, 67.80, 96.7 ,114.27,348.3 ,  &
     &       15.755, 27.620, 40.900, 59.790, 75.00, 91.3 ,124.0 ,143.46,  &
     &        4.339, 31.810, 46.000, 60.900, 82.6 , 99.7 ,118.0 ,155.0 ,  &
     &        6.111, 11.870, 51.210, 67.700, 84.39,109.0 ,128.0 ,147.0 ,  &
     &        6.560, 12.890, 24.750, 73.900, 92.0 ,111.1 ,138.0 ,158.7 ,  &
     &        6.830, 13.630, 28.140, 43.240, 99.8 ,120.0 ,140.8 ,168.5 ,  &
     &        6.740, 14.200, 29.700, 48.000, 65.2 ,128.9 ,151.0 ,173.7 ,  &
     &        6.763, 16.490, 30.950, 49.600, 73.0 , 90.6 ,161.1 ,184.7 ,  &
     &        7.432, 15.640, 33.690, 53.000, 76.0 , 97.0 ,119.24,196.46,  &
     &        7.870, 16.183, 30.652, 54.800, 75.0 , 99.1 ,125.0 ,151.06,  &
     &        7.860, 17.060, 33.490, 51.300, 79.5 ,102.0 ,129.0 ,157.0 ,  &
     &        7.635, 18.168, 35.170, 54.900, 75.5 ,108.0 ,133.0 ,162.0 ,  &
     &        7.726, 20.292, 36.830, 55.200, 79.9 ,103.0 ,139.0 ,166.0 ,  &
     &        9.394, 17.964, 39.722, 59.400, 82.6 ,108.0 ,134.0 ,174.0 ,  &
     &        6.000,  20.509,   30.700, 99.99,99.99,99.99,99.99,99.99,    &
     &        7.89944,15.93462, 34.058, 45.715,99.99,99.99,99.99,99.99,   &
     &        9.7887, 18.5892,  28.351, 99.99,99.99,99.99,99.99,99.99,    &
     &        9.750,21.500, 32.000, 99.99,99.99,99.99,99.99,99.99,        &
     &       11.839,21.600, 35.900, 99.99,99.99,99.99,99.99,99.99,        &
     &       13.995,24.559, 36.900, 99.99,99.99,99.99,99.99,99.99,        &
     &        4.175,27.500, 40.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        5.692,11.026, 43.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        6.2171,12.2236, 20.5244,60.607,99.99,99.99,99.99,99.99,     &
     &        6.63390,13.13,23.17,34.418,80.348,99.99,99.99,99.99,        &
     &        6.879,14.319, 25.039, 99.99,99.99,99.99,99.99,99.99,        &
     &        7.099,16.149, 27.149, 99.99,99.99,99.99,99.99,99.99,        &
     &        7.280,15.259, 30.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        7.364,16.759, 28.460, 99.99,99.99,99.99,99.99,99.99,        &
     &        7.460,18.070, 31.049, 99.99,99.99,99.99,99.99,99.99,        &
     &        8.329,19.419, 32.920, 99.99,99.99,99.99,99.99,99.99,        &
     &        7.574,21.480, 34.819, 99.99,99.99,99.99,99.99,99.99,        &
     &        8.990,16.903, 37.470, 99.99,99.99,99.99,99.99,99.99,        &
     &        5.784,18.860, 28.029, 99.99,99.99,99.99,99.99,99.99,        &
     &        7.342,14.627, 30.490,72.3,99.99,99.99,99.99,99.99,          &
     &        8.639,16.500, 25.299,44.2,55.7,99.99,99.99,99.99,           &
     &        9.0096,18.600, 27.96, 37.4,58.7,99.99,99.99,99.99,          &
     &       10.454,19.090, 32.000, 99.99,99.99,99.99,99.99,99.99,        &
     &       12.12984,20.975,31.05,45.,54.14,99.99,99.99,99.99,           &
     &        3.893,25.100, 35.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        5.210,10.000, 37.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        5.580,11.060, 19.169, 99.99,99.99,99.99,99.99,99.99,        &
     &        5.650,10.850, 20.080, 99.99,99.99,99.99,99.99,99.99,        &
     &        5.419,10.550, 23.200, 99.99,99.99,99.99,99.99,99.99,        &
     &        5.490,10.730, 20.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        5.550,10.899, 20.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        5.629,11.069, 20.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        5.680,11.250, 20.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        6.159,12.100, 20.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        5.849,11.519, 20.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        5.930,11.670, 20.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        6.020,11.800, 20.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        6.099,11.930, 20.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        6.180,12.050, 23.700, 99.99,99.99,99.99,99.99,99.99,        &
     &        6.250,12.170, 20.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        6.099,13.899, 19.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        7.000,14.899, 23.299, 99.99,99.99,99.99,99.99,99.99,        &
     &        7.879,16.200, 24.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        7.86404,17.700, 25.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        7.870,16.600, 26.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        8.500,17.000, 27.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        9.100,20.000, 28.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        8.95868,18.563,33.227, 99.99,99.99,99.99,99.99,99.99,       &
     &        9.220,20.500, 30.000, 99.99,99.99,99.99,99.99,99.99,        &
     &       10.430,18.750, 34.200, 99.99,99.99,99.99,99.99,99.99,        &
     &        6.10829,20.4283,29.852,50.72,99.99,99.99,99.99,99.99,       &
     &        7.416684,15.0325,31.9373,42.33,69.,99.99,99.99,99.99,       &
     &        7.285519,16.679, 25.563,45.32,56.0,88.,99.99,99.99,         &
     &        8.430,19.000, 27.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        9.300,20.000, 29.000, 99.99,99.99,99.99,99.99,99.99,        &
     &       10.745,20.000, 30.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        4.000,22.000, 33.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        5.276,10.144, 34.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        6.900,12.100, 20.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,        &
     &        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99/

!
!     An element (hydrogen through zinc) can be considered in one of
!     the three following options:
!     1. explicitly - some of energy levels of some of its ionization
!                     states are considered explicitly, ie. their
!                     populations are determined by solving statistical
!                     equilibrium
!     2. implicitly - the atom is assumed not to contribute to
!                     opacity; but is allowed to contribute to the
!                     total number of particles and to the total charge;
!                     the latter is evaluated assuming LTE ionization
!                     balance, ie. by solving a set of Saha equations
!     3. not considered at all
!
!     Input:
!
!     For each element from 1 (hydrogen) to NATOMS, the following
!     parameters:
!
!     MA     =  0  - if the element is not considered (option 3)
!            =  1  - if the element is non-explicit (option 2)
!            =  2  - if the element is explicit (option 1)
!            =  4  - if the element is semi-explicit (i.e. behaves
!                    like MA=2 for continua and MA=1 for lines
!     NA0,NAK - have the meaning only for MA=2; indicate that the
!               explicit energy levels of the present species have
!               the indices between NA0 and NAK (NAK is thus the index
!               of the highest ionization state, which is represented
!               as one-level ion).
!     ION    -  has the meaning for MA=1 only;
!               if ION=0, standard number of ionization degrees is
!                         considered
!                         (counting the neutral state also; so for
!                         instance to treat all stages of He requires
!                         ION=3, which is a default anyhow).
!               if ION>0, then ION ionization degrees is considered
!     MODPF  -  mode of evaluation of partition functions
!            =  0  -  standard evaluation (see procedure PARTF)
!            >  0  -  partition functions evaluated from the
!                     Opacity Project ionization fraction tables
!            <  0  -  non-standard evaluation, by user supplied
!                     procedure PFSPEC
!     ABN    -  if ABN=0, solar abundance is assumed (given above;
!                         abundance here is assumed as relative
!                         to hydrogen by number
!               if ABN>0, non-solar abundance ABN is assumed; in an
!                         arbitrary scale
!               if ABN<0, non-solar abundance ABN is assumed;
!                         (-ABN times the solar value)
!     PFS    -  see above
!
      iabset=0
      natoms=99
      ND=1
!     read(ibuff,'(a80)') dum
!     read(dum,*,iostat=kstat) natoms,iabset
!     if(kstat.ne.0) READ(dum,*) NATOMS
!!    WRITE(6,"(//' CHEMICAL ELEMENTS INCLUDED'/                          &
!!   &            ' --------------------------'//                         &
!!   & ' NUMBER  ELEMENT           ABUNDANCE'/1H ,16X,                    &
!!   & 'A=N(ELEM)/N(H)  A/A(SOLAR)'/)")
      IAT=0
      IREF=1
!
      DO I=1,MATOM
         if(iabset.eq.1) then
            d(2,i)=10.**(abun1(i)-12.)
          else if(iabset.ne.2) then
            d(2,i)=10.**(abun0(i)-12.)
         end if
      END DO
      DO ID=1,ND
         YTOT(ID)=0.
         WMY(ID)=0.
      END DO
!
      DO I=1,MATOM
         TYPAT(I)=DYP(I)
!        LGR(I)=.TRUE.
!        LRM(I)=.TRUE.
         AMAS(I)=D(1,I)
         ABND(I)=D(2,I)
         IONIZ(I)=int(D(3,I))
         isemex(i)=0
!
         DO J=1,9
            IF(J.LE.8) ENEV(I,J)=xi(J,I)
         END DO
         MA=1
         ion=3
         if(i.eq.1) ion=2
         IF(MA.GT.0) THEN
            LGR(I)=.FALSE.
!           IF(ABN.GT.0) ABND(I)=ABN
!           IF(ABN.LT.0) ABND(I)=ABS(ABN)*D(2,I)
            IF(ION.NE.0) IONIZ(I)=ION
            DO ID=1,ND
               ABNDD(I,ID)=ABND(I)
            END DO
            LRM(I)=.FALSE.
            IATEX(I)=0
            DO ID=1,ND
               YTOT(ID)=YTOT(ID)+ABNDD(I,ID)
               WMY(ID)=WMY(ID)+ABNDD(I,ID)*AMAS(I)
            END DO
            ABN=ABND(I)/D(2,I)
!           WRITE(6,"(I4,3X,A5,1P2E14.2)") I,TYPAT(I),ABND(I),ABN
         END IF
      END DO
      IF(MOD(IMODE,10).LE.1) NATOMS=MATOM
      DO ID=1,ND
         WMM(ID)=WMY(ID)*HMASS/YTOT(ID)
      END DO

      READ(56,*,IOSTAT=IOS) NCHANG
      IF(IOS.NE.0) RETURN
      IF(NCHANG.LE.0) RETURN

      WRITE(6,"(//'CHEMICAL ABUNDANCES  CHANGED'                          &
     &           /'----------------------------'/                         &
     &      ' NUMBER  ELEMENT           ABUNDANCE'/1H ,16X,               &
     &      'A=N(ELEM)/N(H)  A/A(SOLAR)'/)")
      DO II=1,NCHANG
         READ(56,*) I,ABN
         ABND(I)=D(2,I)
         IF(ABN.GT.0) ABND(I)=ABN
         IF(ABN.LT.0) ABND(I)=-ABN*D(2,I)
!        IF(ABN.LT.0) ABND(I)=-ABN*ABNDD(I,1)
         if(abn.gt.1.) abnd(i)=10.**(abn-12.)
         DO ID=1,ND
            ABNDD(I,ID)=ABND(I)
         END DO
         ABNR=ABND(I)/D(2,I)
         WRITE(6,"(1X,I4,3X,A5,1P2E14.2)") I,TYPAT(I),ABND(I),ABNR
      END DO
!
      RETURN
      END SUBROUTINE STATE0

!
!     ***********************************************************************
!


      subroutine mpartf(jatom,ion,indmol,t,u)
!     =======================================
!
!     yields partition functions with polynomial data from
!     ref. Irwin, A.W., 1981, ApJ Suppl. 45, 621.
!     ln u(temp)=sum(a(i)*(ln(temp))**(i-1)) 1<=a<=6
!
!     Input:
!       jatom = element number in periodic table
!       ion   = 1 for neutral, 2 for once ionized and 3 for twice ionized
!       indmol= index of a molecular specie (Tsuji index)
!       temp  = temperature
!     Output:
!       u     = partf.(linear scale) for iat,ion, or indmol, and temperature t
!
!
      use accura
      use params, only : irwtab
      implicit real(dp) (a-h,o-z)

      real(dp), save :: a(6,3,92),aa(6),am(6,500)
      integer        :: indtsu(324)
      integer, save  :: irw(500), nummol

      data indtsu / 2,  5, 12, 4, 8, 7, 6,                                &
     &              9, 11, 10, 29, 50, 59, 46, 133,  52, 19,              &
     &             13, 42, 38, 39, 37, 44, 36,  14, 118, 33,              &
     &              3, 16, 57, 32, 49, 60, 54,  41, 107,304,              &
     &            148, 152, 153, 155, 303, 17,  24,  25, 28, 51,          &
     &            112, 119, 102,   0,  21, 15,  43,  22,478, 64,          &
     &             47,  65, 414,  61, 191, 62 ,109,  40, 66,214,          &
     &             120*0,    30, 136*0/

      data iread /0/
!
!     read data if first call:
!
      if(iread.ne.1) then
        if(irwtab.eq.0) then
           open(67,file= './data/irwin_orig.dat',status='old')
           nummol=66
         else
           open(67,file= './data/irwin_bc.dat',status='old')
           nummol=324
        end if
        read(67,*)
        read(67,*)
        do j=1,92
          icyk: do i=1,3
            if(j.eq.1.and.i.eq.3) cycle icyk
            sp=float(j)+float(i-1)/100.
            read(67,*) spec,aa
            do k=1,6
                a(k,i,j)=aa(k)
            end do
          end do icyk
       end do
!
       read(67,*)
       read(67,*)
       read(67,*)
       do i=1,500
          irw(i)=0
       end do
       do i=1,nummol
          read(67,*) spec,aa
          indm=indtsu(i)
          if(indm.gt.0) then
             irw(indm)=i
             do j=1,6
                am(j,indm)=aa(j)
             end do
          end if
        end do
        close(67)
        iread=1
      end if
!
!     evaluation of the partition function
!     stop if T is out of limits of Irwin's tables
!
        if(t.lt.1000.) then
          stop 'partf; temp<1000 K'
          write(*,*)  'warning: temp<1000 K in PARTF'
        else if(t.gt.16000.) then
          write(*,*) 'mpartf T=',t
          stop 'partf; temp>16000 K'
        endif
        tl=log(t)
        u=0.
!
!     atomic species
!
      if(jatom.gt.0.and.ion.gt.0) then
        ulog=    a(1,ion,jatom)+                                          &
     &       tl*(a(2,ion,jatom)+                                          &
     &       tl*(a(3,ion,jatom)+                                          &
     &       tl*(a(4,ion,jatom)+                                          &
     &       tl*(a(5,ion,jatom)+                                          &
     &       tl*(a(6,ion,jatom))))))
        if(jatom.eq.5.and.ion.eq.3) ulog=1.
        u=exp(ulog)
      end if
!
!     molecular species
!
      if(indmol.gt.0) then
        indm=indmol
        if(irw(indm).gt.0) then
           ulog=    am(1,indm)+                                           &
     &          tl*(am(2,indm)+                                           &
     &          tl*(am(3,indm)+                                           &
     &          tl*(am(4,indm)+                                           &
     &          tl*(am(5,indm)+                                           &
     &          tl*(am(6,indm))))))
           u=exp(ulog)
        end if
      end if
      return
      end subroutine mpartf
!
!     ***********************************************************************
!
      subroutine tsuji_orig(tem,aptsuji)
!     ==================================
!
!     equilibrium constant for original Tsuji data
!
      use accura
!     use params
!     use modelp
!     use eospar
      implicit real(dp) (a-h,o-z),logical (l)

      real(dp) :: aptsuji(600)
      real(dp), save :: c0(600,5)
      character(len=128) :: MOLEC
      character(len=8)    :: cml
!     INTEGER        :: NATOMM(5),NELEMM(5)

      DATA IREAD /1/
!
!---- read molecular data from a table  ----------------------
!       
!     MOLTAB=0
!     if(moltab.eq.0) 
      MOLEC='data/tsuji.molec_orig'

      if(iread.eq.1) then
        J=0
        OPEN(UNIT=26,FILE=MOLEC,STATUS='OLD')
        READMOL2: DO
           J=J+1
           READ (26,"(a8,5e13.5,9i3)",IOSTAT=IOS)                             &
     &        cml,(C0(J,K),K=1,5) !,MMAX(J),(NELEMM(M),NATOMM(M),M=1,4)
!          write(*,*) 'j',j,cmol(j),c0(j,1) 
           IF(IOS.NE.0) EXIT READMOL2
!          MMAXJ=MMAX(J)
!          IF(MMAXJ.EQ.0) THEN
!             EXIT READMOL2
!          END IF
        END DO READMOL2
        
        NMOLEC0=J-1
        close(26)
      end if

      t=5040.4/tem
      aptsuji=0.
      DO J=1,nmolec0
         APLOGJ=C0(J,5)
         DO K=1,4
            KM5=5-K
            APLOGJ=APLOGJ*T + C0(J,KM5)
         END DO
         aptsuji(j)=aplogj
      END DO

      RETURN
      END SUBROUTINE TSUJI_ORIG

!
!     ***********************************************************************
!

      subroutine bceqco(in,t,bcel)
!     ===========================
!
!     new valuation of BC equiloibrium constant
!
!     input:  in   - BC index
!             t    - temperature
!     output: bcel - log(Kp)

      use accura
      implicit none

      integer, intent(in)   :: in
      real(dp), intent(in)  :: t
      real(dp), intent(out) :: bcel
      real(dp)              :: tt,a1

      real(dp),save :: tl(42),bceql(291,42)
      real(dp) :: x(42)
      character(len=5) :: molab(291)
      character(len=7) :: dum
      integer          :: iread,i,j,ios
      data iread /1/

      if(iread.eq.1) then
         open(67,file='data/BCtab7.dat',status='old')
         read(67,*)
         read(67,*)
         read(67,"(a7,42e14.5)") dum,(tl(i),i=1,42)
         read(67,*)
         do j=1,291
            read(67,"(a5,42e14.5)",iostat=ios)                            &
     &      molab(j),(x(i),i=1,42)
!    &      molab,(bceql(in,i),i=1,42)
            do i=1,42
               bceql(j,i)=x(i)
            end do
         end do
         close(67)
         do i=1,42
           tl(i)=log10(tl(i))
         end do
      end if

      if(t.gt.1.e4) return
      tt=log10(t)

      do i=1,42
         if(tt.le.tl(i)) exit
      end do

      a1=(tt-tl(i-1))/(tl(i)-tl(i-1))
      bcel=a1*bceql(in,i)+(1.-a1)*bceql(in,i-1)+1.

      return
      end subroutine bceqco                                                  

!
!     ***********************************************************************
!
      
      subroutine bcnega(in,t,bcnl)
!     ===========================
!     
!     valuation of BC partition functions for some negative ions
!
      use accura
      implicit none
         
      integer, intent(in)   :: in                                         
      real(dp), intent(in)  :: t
      real(dp), intent(out) :: bcnl
      real(dp)              :: tt,a1 
     
      real(dp),save :: tl(42),bcnel(7,42)
      character(len=6) :: molab(7)
      character(len=9) :: dum
      integer          :: iread,i,j
      data iread /1/

      if(iread.eq.1) then
         open(67,file='data/BCpfneg.dat',status='old')
         read(67,*)
         read(67,*)
         read(67,"(a9,42e14.5)") dum,(tl(i),i=1,42)
         read(67,*)
         do j=1,7
            read(67,"(a6,42e14.5)") molab(j),(bcnel(j,i),i=1,42)
         end do
         close(67)
         iread=0
         do i=1,42
           tl(i)=log10(tl(i))
         end do
      end if
      
      if(t.gt.1.e4) return
      tt=log10(t)
      
      do i=1,42
         if(tt.le.tl(i)) exit
      end do
      
      a1=(tt-tl(i-1))/(tl(i)-tl(i-1))
      bcnl=a1*bcnel(in,i)+(1.-a1)*bcnel(in,i-1)

      return
      end subroutine bcnega
      
!
!     ***********************************************************************
!




      subroutine elebear(ielem)
!     =========================
!
      use accura
      use eospar
      use modelp, only : cmol
      use params, only : nmolec
      implicit none

      integer, intent(in) :: ielem
      real(dp)            :: anm(600)
      real(dp)            :: antot
      integer             :: indx(600),ints(600)
      integer             :: id,indb,j,m,ibear,nbear,i,in
      character(len=3)    :: dyp(92)
      character(len=8)    :: labl

      DATA DYP/' H ',' He',' Li',' Be',' B ',' C ',
     &         ' N ',' O ',' F ',' Ne',' Na',' Mg',                       &
     &         ' Al',' Si',' P ',' S ',' Cl',' Ar',                       &
     &         ' K ',' Ca',' Sc',' Ti',' V ',' Cr',                       &
     &         ' Mn',' Fe',' Co',' Ni',' Cu',' Zn',                       &
     &         ' Ga',' Ge',' As',' Se',' Br',' Kr',                       &
     &         ' Rb',' Sr',' Y ',' Zr',' Nb',' Mo',                       &
     &         ' Tc',' Ru',' Rh',' Pd',' Ag',' Cd',                       &
     &         ' In',' Sn',' Sb',' Te',' I ',' Xe',                       &
     &         ' Cs',' Ba',' La',' Ce',' Pr',' Nd',                       &
     &         ' Pm',' Sm',' Eu',' Gd',' Tb',' Dy',                       &
     &         ' Ho',' Er',' Tm',' Yb',' Lu',' Hf',                       &
     &         ' Ta',' W ',' Re',' Os',' Ir',' Pt',                       &
     &         ' Au',' Hg',' Tl',' Pb',' Bi',' Po',                       &
     &         ' At',' Rn',' Fr',' Ra',' Ac',' Th',                       &
     &         ' Pa',' U '/

      id=1
      indb=0
      antot=0.
      molloop: do j=1,nmolec
         ibear=0
         do m=1,mmax(j)
            if(nelem(m,j).eq.ielem) ibear=nato(m,j)
         end do
         if(ibear.eq.0) cycle molloop 
         indb=indb+1
         ints(indb)=j
         anm(indb)=anmol(j,id)*float(ibear) 
         antot=antot+anm(indb) 
!        write(*,"(i4,a8,i4,1p2e11.3)") j,cmol(j),indb,anm(indb),antot
      end do molloop  

      nbear=indb+2
      anm(indb+1)=anato(ielem,id)
      anm(indb+2)=anion(ielem,id)
      antot=antot+anm(indb+1)+anm(indb+2)
      call indexx(nbear,anm,indx)
!     write(*,*)
     
!     do i=1,nbear
!        write(*,"(2i4,1pe11.3)") i,indx(i),anm(indx(i))
!     end do
!     write(*,*)

      write(*,"(/'Distribution of nuclei of',a8
     &          /'================================='/
     &   'species     N_nuc     N/N_tot'/)") dyp(ielem)

      do i=nbear,1,-1
         in=indx(i)
         if(anm(in)/antot.lt.1.e-15) cycle
         if(in.eq.nbear-1) then
           labl=dyp(ielem)
          else if(in.eq.nbear) then
!          labl=dyp(ielem)//'+'
           labl=trim(dyp(ielem))//'+'
          else
           labl=cmol(ints(in))
         end if
 
         write(6,"(a8,1p2e11.3)") labl,anm(in),anm(in)/antot
      end do

      end subroutine elebear
             


C
C           
C ********************************************************************
C           
C           
      SUBROUTINE INDEXX(N,ARRIN,INDX)
C     ===============================
C           
C     Sorting routine
C              
      use accura 
      implicit none

      REAL(DP), INTENT(IN)  :: ARRIN(N)
      INTEGER,INTENT(OUT)   :: INDX(N)
      INTEGER               :: I,J,N,M,IR,INDXT
      REAL(DP)              :: QQ
               
      DO J=1,N 
         INDX(J)=J
      END DO      
      M=N/2+1   
      IR=N        
      OUTER: DO
      IF(M.GT.1)THEN
         M=M-1
         INDXT=INDX(M)
         QQ=ARRIN(INDXT)
       ELSE
         INDXT=INDX(IR)
         QQ=ARRIN(INDXT)
         INDX(IR)=INDX(1)
         IR=IR-1  
         IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN   
         END IF      
      END IF       
      I=M            
      J=M+M          
      DO WHILE(J.LE.IR)
         IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
         END IF
         IF(QQ.LT.ARRIN(INDX(J))) THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
         END IF
      END DO
      INDX(I)=INDXT
      END DO OUTER

      END SUBROUTINE INDEXX
C
