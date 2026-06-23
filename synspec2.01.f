      include 'syn_modules.f' 
 
      PROGRAM SYNSPEC 
! 
! =====================================================================I 
!                                                                      I 
! Program for evaluting synthetic spectra for a given model atmosphere I 
!                                                                      I 
! *******************                                                  I 
! VERSION SYNSPEC2.01                                                  I 
! *******************                                                  I 
!                                                                      I 
! Main input files:                                                    I
!        the same as input to TLUSTY - unit 5                          I
!        keyword parameters                                            I 
!        model atmosphere - unit 8                                     I 
!        line list(s)                                                  I 
!                                                                      I 
! Output: diagnostic outprint - unit 6 (several procedures)            I 
!         synthetic spectrum  - unit 7 (procedure OUTPRI)              I 
!         flux in continuum   - unit 17 (procedure OUTPRI)             I 
!         identification table- unit 12 (procedure INIBLA)             I 
!         partial equiv.widths- unit 16 (procedure OUTPRI)             I 
!         elapsed time        - unit 69 (procedure TIMING)             I 
!                                                                      I 
!      -- if specific intensities are also calculated (set up by the   I 
!         input on unit 55), there are two aditional output files:     I 
!                                                                      I 
!         specific intensities - unit 10                               I 
!         specific intensities in continuum - unit 18                  I 
!                                                                      I 
! Basic options: controlled by switch IMODE                            I 
! IMODE    =  0 - normal synthetic spectrum                            I 
!                 (ie. identification table + emergent flux)           I 
!          =  1 - detailed profiles of a few individual lines          I 
!          =  2 - emergent flux in the continuum (without the          I 
!                 contribution of lines)                               I 
!          = -1 - only identification table, ie. a list of lines which I 
!                 contribute to opacity in a given wavelength          I 
!                 region, together with their approximate equivalent   I 
!                 widths. Synthetic spectrum is not calculated.        I 
!          = -2 - the "iron curtain" option, ie. a monochromatic       I 
!                 opacity for a homogeneous slab of a given T and n_e  I 
!          = -3 - opacity table                                        I
!          = -4 - opacity table eith the continuum only                I      
!                                                                      I 
!                                                                      I 
! ==================================================================== I 
! 
! 
      use accura 
      use params 
      use lindat 
      use modelp 
      use synthp 
      use molist 
 
      implicit real(dp) (a-h,o-z),logical (l) 
! 
      write(*,*) 
      write(*,*) '------------------------------------' 
      write(*,*) 'RUN WITH syn62j GENERATED 05/23/2026' 
      write(*,*) '------------------------------------' 
      write(*,*) 
      OPEN(UNIT=12,STATUS='UNKNOWN') 
      OPEN(UNIT=14,STATUS='UNKNOWN') 
! 
!     INITIALIZATION - INPUT OF BASIC PARAMETERS AND MODEL ATMOSPHERE 
!                      (independently on the basic mode) 
!!!   call syn_alloc 
 
      CALL INITIA 
      if(ifeos.gt.0) imode=-3 
      if(ibfac.gt.1) then 
         LTE0=LTE 
         LTE=.TRUE. 
      END IF 
!     if(imode.le.-3) write(*,*) 'imode,inmod,nd',imode,inmod,nd 
      IF(IMODE.GE.-2.AND.IFEOS.LE.0) THEN 
         IF(INMOD.GT.0) CALL INPMOD 
         IF(INMOD.EQ.0) CALL INKUR 
         IF(ICHANG.NE.0) CALL CHANGE 
         IF(IBFAC.GT.1) THEN 
            CALL INPBF 
            LTE=LTE0 
         END IF 
         IF(IFWIN.GT.1) CALL SETWIN 
       ELSE 
         IDSTD=1 
         CALL INGRID(0,inext,0) 
      END IF 
! 
      CALL INIBL0 
      CALL INIMOD 
      CALL TINT 
! 
      IMODE0=IMODE 
      IF(IMODE0.EQ.-4) IMODE=2 
! 
!     ***** loop over grids (opacity table T,rho pairs) 
! 
      igrd=0 
      GRIDLOOP: DO 
! 
         IF(IMODE0.LE.-3.and.ifeos.le.0) CALL INIBL1(IGRD) 
         IF(IFMOL.GT.0) then 
            CALL MOLINI 
!           write(*,*) 'BEF EOSPRI in MAIN',igrd 
            if(ifeos.ne.0) call eospri 
!           write(*,*) 'AFT EOSPRI in MAIN',igrd 
         end if 
! 
!        zero abundances for selected species (if required) 
! 
         if(imode0.le.-3) call abnchn(1) 
! 
         IBLANK=0 
         NXTSET=0 
         IF(IFMOL.GT.0.AND.IMODE.LT.2) THEN 
            DO ILIST=1,NMLIST 
               NXTSEM(ILIST)=0 
               INACTM(ILIST)=0 
               NLINMT(ILIST)=0 
            END DO 
         END IF 
! 
!        read line list(s) 
! 
         IF(IFEOS.LE.0) THEN 
            IF(IMODE.LT.2) CALL INILIN 
            IF(IFMOL.GT.0.AND.IMODE.LT.2) THEN 
               DO ILIST=1,NMLIST 
                  IF(IMODE.EQ.-3.AND.TEMP(1).LT.TMLIM(ILIST))             & 
     &            CALL INMOLI(ILIST) 
                  IF(IMODE.GE.-2.and.imode.le.1) CALL INMOLI(ILIST) 
!!              write(*,*) 'after INMOLI with ilist',ilist 
               END DO 
            END IF 
         END IF 
! 
!        ACTUAL CALCULATION OF THE SYNTHETIC SPECTRUM 
! 
         EOSIF: IF(IFEOS.LE.0) THEN 
            IEND=0 
            BLANKLOOP: DO 
               IBLANK=IBLANK+1 
!!        write(*,"(//' *** BEF RESOLV - IGRD,IBLANK,NBLANK,INEXT', 
!!   *           5i8/)") igrd,iblank,nblank,inext,nd 
               IF(IFWIN.LE.0) THEN 
                  CALL RESOLV 
                ELSE 
                  CALL RESOLW 
               END IF 
               IF(IMODE0.GE.0) THEN 
                  if(ifreq.le.10.and.inmod.le.1) then 
                     CALL RTECD 
                   else 
                     call RTE 
                  END IF 
                  CALL OUTPRI 
               END IF 
               if((imode.ge.0.and.imode.ne.7.and.iprin.ge.1).or.          & 
     &            (imode.lt.0.and.iprin.ge.2)) then 
                  CALL IDTAB 
                  IF(IFMOL.GT.0) CALL IDMTAB 
               end if 
!!             write(*,"(//' *** AFT RESOLV - IBLANK,NBLANK,IEND', 
!!   *         4i8/)") iblank,nblank,iend 
               IF(IBLANK.LT.NBLANK) THEN 
                  CYCLE BLANKLOOP 
                ELSE 
                  IEND=1 
               END IF 
!               write(*,"(/' *** AFTER BLANKLOOP: NXTSET,IRLIST', 
!    *          4i8//)") nxtset,irlist,nd 
               IF(NXTSET.EQ.1.AND.IRLIST.EQ.0) THEN 
                  IF(IMODE.LT.2) THEN 
                     CALL INILIN 
                     CYCLE BLANKLOOP 
                   ELSE 
                     IEND=1 
                  END IF 
               END IF 
!               write(*,*) 'IRLIST',irlist,nd 
               IF(IFMOL.GT.0.AND.IMODE.LT.2.AND.IRLIST.GT.0) THEN 
                  DO ILIST=1,NMLIST 
!!                   write(*,*) 'ILIST,NXTSEM,INACTM',ilist, 
!!   *               nxtsem(ilist),inactm(ilist) 
                     IF(NXTSEM(ILIST).EQ.1.and.inactm(ilist).eq.0) THEN 
                        CALL INMOLI(ILIST) 
                        IBLANK=0 
!                        write(*,*) 'IBLANK in IF',iblank 
                        CYCLE BLANKLOOP 
                       ELSE 
                        IEND=1 
                      END IF 
                  END DO 
               END IF 
!               write(*,*) 'IEND near end - iend,igrd',iend,igrd 
               IF(IEND.EQ.1) EXIT BLANKLOOP 
            END DO BLANKLOOP 
         END IF EOSIF 
! 
         ifin=1 
         if(imode0.lt.-2) then 
            call ingrid(1,inext,igrd) 
            igrd=igrd+1 
            if(inext.gt.0) cycle gridloop 
         end if 
         if(ifin.eq.1) exit gridloop 
      END DO GRIDLOOP 
! 
!     ***** iend if loop over grids (opacity table T,rho pairs) 
! 
 
      if(imode0.le.-3.and.ifeos.le.0) call fingrd 
      call timing(2,iblank) 
 
      close(6) 
      close(7) 
      close(12) 
      close(14) 
      close(15) 
      close(16) 
      close(17) 
      close(69) 
      close(51) 
      close(52) 
      close(53) 
      close(54) 
      close(66) 
      close(67) 
      close(77) 
 
      END PROGRAM SYNSPEC 
! 
! ******************************************************************** 
! 
!     allocation routines 
! 
      include 'syn_alloc.f' 
! 
! 
!     **************************************************************** 
! 
! 
 
      SUBROUTINE INITIA 
!     ================= 
! 
!     driver for input and initializations 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP),PARAMETER :: WI1=911.753578, WI2=227.837832 
      CHARACTER(LEN=4)   ::TYPION(MIOEX),TYPIOI 
      CHARACTER(LEN=40)  :: FILEI 
      CHARACTER(LEN=20)  :: FINSTD 
      CHARACTER(LEN=1)   :: BLNK  = ' ' 
      INTEGER            :: IATI(MION),IZI(MION) 
      INTEGER            :: IGLE(18),IGMN(25),IGFE(26),IGNI(28) 
      IGLE = (/2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1/) 
      IGMN = (/2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,                       & 
     &          10,21,28,25,6,7,6/) 
      IGFE = (/2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,                       & 
     &          10,21,28,25,6,25,30,25/) 
      IGNI = (/2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,                       & 
     &          10,21,28,25,6,25,28,21,10,21/) 
! 
      CALL READBF 
! 
! ------------------------------------ 
! Basic input parameters - atmospheres 
! ------------------------------------ 
! 
      IF(INMOD.LE.1) THEN 
         READ(IBUFF,*) TEFF,GRAV 
       ELSE IF(INMOD.EQ.2) THEN 
! 
! ------------------------------ 
! Basic input parameters - disks 
! ------------------------------ 
! 
         READ(IBUFF,*) DISPAR 
      END IF 
! 
! ---------------------------- 
! other basic input parameters 
! ---------------------------- 
! 
      READ(IBUFF,*) LTE,LTGREY 
! 
      CALL READKW
!
      READ(IBUFF,*) FINSTD 
      IF(IFKEY.EQ.0) CALL NSTPAR(FINSTD) 
      IF(IMODE.LT.-1) THEN 
         ND=1 
         IDSTD=1 
      END IF 
      IF(LTE) INLTE=0 
      CALL PRKEYW
! 
!     allocate some arrays 
! 
      call syn_alloc 
      if(ihyddk.gt.0.or.ihgom.gt.0) call alloc_hydprf 
! 
! ---------------------------- 
! Frequency points and weights 
! ---------------------------- 
! 
      READ(IBUFF,*) NFREAD 
      NJREAD=NFREAD 
! 
      IF(NJREAD.LT.0) THEN 
         NJREAD=-NJREAD 
         NFREQC=NJREAD 
         DO IJ=1,NJREAD 
            READ(IBUFF,*) FREQEXP 
         END DO 
       ELSE 
         NFREQC=NJREAD 
      END IF 
! 
!     if(imode.ge.0)                                                      & 
!    &WRITE(6,"(31X,'*******************************************'/        & 
!    & 31X,'I',41X,'I'/                                                   & 
!    & 31X,'I   S Y N T H E T I C   S P E C T R U M   I'/                 & 
!    & 31X,'I',41X,'I'/                                                   & 
!    & 31X,'I',8X,'FOR MODEL ATMOSPHERE WITH',8X,'I'/                     & 
!    & 31X,'I',41X,'I'/                                                   & 
!    & 31X,'I',14X,'TEFF  =',F7.0,13X,'I'/                                & 
!    & 31X,'I',14X,'LOG G =',F7.2,13X,'I'/                                & 
!    & 31X,'I',41X,'I'/                                                   & 
!    & 31X,'*******************************************')") TEFF,GRAV 
! 
! ---------------------------------------------------- 
!     turbulent velocities 
! ---------------------------------------------------- 
! 
      IF(VTB.LT.1.E3) VTB=VTB*1.E5 
      DO ID=1,ND 
         VTURB(ID)=VTB 
      END DO 
! 
! ---------------------------------------------------- 
! Input parameters for explicit and non-explicit atoms 
! ---------------------------------------------------- 
! 
!     Input parameters are read by procedure STATE 
!     (see description there) 
! 
      CALL STATE0(1) 
      ID=1 
      IF(IPRIN.GE.1)                                                      & 
     &   WRITE(6,"(//' ----------------------'/                           & 
     & ' YTOT   =',F11.5/' WMY    =',1PE15.5/                             & 
     & ' WMM    =',E15.5)") YTOT(ID),WMY(ID),WMM(ID) 
      DO I=1,MLEVEL 
         ILK(I)=0 
         iexpl(i)=0 
         iltot(i)=0 
      END DO 
! 
! -------------------------------------------------------------- 
! Input of parameters for explicit ions, levels, and transitions 
! -------------------------------------------------------------- 
! 
      ILEV=0 
      IATLST=0 
      ION=0 
      IA=0 
      IUNIT=34 
      NATOM=0 
      WRITE(6,"(//' EXPLICIT IONS INCLUDED'/                              & 
     &            ' ----------------------'//                             & 
     &  ' ION     N0    N1    NK    IZ'/)") 
 
      LEVELS: DO 
         READ(IBUFF,*,IOSTAT=IOS) IATII,IZII,NLEVSI,ILASTI,ILVLIN,        & 
     &              NONSTD,TYPIOI,FILEI 
         IF(IOS.NE.0) EXIT LEVELS 
         IF(ILASTI.EQ.0) THEN 
            ION=ION+1 
            IATI(ION)=IATII 
            IZI(ION)=IZII 
            NLEVS(ION)=NLEVSI 
            TYPION(ION)=TYPIOI 
            FIDATA(ION)=FILEI 
            NLLIM(ION)=ILVLIN 
            ILIMITS(ION)=-1 
            IUPSUM(ION)=0 
            FIBFCS(ION)=BLNK 
            MODEFF=1 
            NFF=0 
            IF(IATI(ION).EQ.1.AND.IZI(ION).EQ.0) THEN 
               IUPSUM(ION)=-100 
               MODEFF=2 
            END IF 
            IF(IATI(ION).EQ.2.AND.IZI(ION).EQ.1) THEN 
               MODEFF=2 
            END IF 
            IF(NONSTD.GE.10) THEN 
              WRITE(*,*)'INITIA: QUANTUM NUMBERS AND ENERGY LIMITS WILL' 
              WRITE(*,*)'        BE IGNORED FOR ION ',IATII,'    ',IZII 
              ILIMITS(ION)=0 
              NONSTD=NONSTD-10 
            END IF 
            IF(NONSTD.GT.0) THEN 
               READ(IBUFF,*) IUPSUM(ION),ICUP,MODEFF,NFF 
             ELSE IF(NONSTD.LT.0) THEN 
               READ(IBUFF,*) ifil1,ifil2,FIODF1(ION),                     & 
     &                    FIODF2(ION),FIBFCS(ION) 
               IF(FIBFCS(ION).NE.' ') THEN 
                  IUNIT=IUNIT+1 
                  INBFCS(ION)=IUNIT 
               END IF 
               IUPSUM(ION)=1 
            END IF 
! 
            IF(IATI(ION).EQ.IATLST) THEN 
               NFIRST(ION)=ILEV 
             ELSE 
               NFIRST(ION)=ILEV+1 
               IATLST=IATI(ION) 
               IA=IATEX(IATLST) 
               N0A(IA)=NFIRST(ION) 
               NATOM=MAX(NATOM,IA) 
            END IF 
            NLAST(ION)=NFIRST(ION)+NLEVS(ION)-1 
            NNEXT(ION)=NLAST(ION)+1 
            ILEV=NNEXT(ION) 
            IZ(ION)=IZI(ION)+1 
            IF(NFF.GT.0) FF(ION)=EH/H*IZ(ION)*IZ(ION)/NFF/NFF 
! 
            N0I=NFIRST(ION) 
            N1I=NLAST(ION) 
            NKI=NNEXT(ION) 
            IFREE(ION)=MODEFF 
            DO II=N0I,N1I 
               IEL(II)=ION 
               IATM(II)=IA 
            END DO 
            ILK(NKI)=ION 
            IATM(NKI)=IA 
! 
            IF(NUMAT(IA).EQ.1) THEN 
               IATH=IA 
               IF(IZ(ION).EQ.1) IELH=ION 
               IF(IZ(ION).EQ.0) IELHM=ION 
            END IF 
            IF(NUMAT(IA).EQ.2) THEN 
               IATHE=IA 
               IF(IZ(ION).EQ.1) IELHE1=ION 
               IF(IZ(ION).EQ.2) IELHE2=ION 
            END IF 
! 
            IF(IPRIN.GE.0)                                                & 
     &       WRITE(6,"(A4,4I6)") TYPION(ION),N0I,N1I,NKI,IZ(ION) 
! 
          ELSE IF(ILASTI.GT.0) THEN 
            ENION(ILEV)=0. 
            G(ILEV)=ILASTI 
            NQUANT(ILEV)=1 
            TYPLEV(ILEV)=TYPIOI 
            IFWOP(ILEV)=0 
            IEL(ILEV)=ION 
            NKA(IA)=NNEXT(ION) 
            IF(ILASTI.EQ.1.AND.IATII.GT.IZII) THEN 
               IF(IATII.LT.25) THEN 
                  G(ILEV)=IGLE(IATII-IZII) 
                ELSE IF(IATII.EQ.25) THEN 
                  G(ILEV)=IGMN(IATII-IZII) 
                ELSE IF(IATII.EQ.26) THEN 
                  G(ILEV)=IGFE(IATII-IZII) 
                ELSE IF(IATII.EQ.28) THEN 
                  G(ILEV)=IGNI(IATII-IZII) 
               ENDIF 
            ENDIF 
          ELSE 
            EXIT LEVELS 
         END IF 
      END DO LEVELS 
 
      NION=ION 
      NLEVEL=NKI 
! 
      if(iath.gt.0) then 
         N0H=N0A(IATH) 
         N1H=NLAST(IELH) 
         NKH=NNEXT(IELH) 
         N0HN=NFIRST(IELH) 
         N0M=0 
         IF(IELHM.GT.0) THEN 
            N0M=NFIRST(IELHM) 
            IOPHMI=0 
         end if 
       else 
         n0h=0 
         n1h=0 
         nkh=0 
         n0hn=0 
      end if 
! 
      IF(IPRIN.GE.1) WRITE(6,"(//' BASIC INPUT PARAMETERS'/               & 
     &            ' ----------------------'/                              & 
     &            ' INMOD  =',I5/                                         & 
     &            ' ND     =',I5/                                         & 
     &            ' IDSTD  =',I5/                                         & 
     &            ' NATOM  =',I5/                                         & 
     &            ' NION   =',I5/                                         & 
     &            ' NLEVEL =',I5/                                         & 
     &            ' IELH   =',I5/                                         & 
     &            ' IELHM  =',I5/)")                                      & 
     &   INMOD,ND,IDSTD,NATOM,NION,NLEVEL,                                & 
     &   IELH,IELHM 
! 
! ----------------------------------------- 
! Parameters for individual explicit levels 
! ----------------------------------------- 
! 
!     IMER=0 
!     ITR=0 
!     IC=0 
!     IL=0 
!     IP=0 
      NMER=0 
! 
      DO ION=1,NION 
         CALL RDATA(ION) 
         NFF=NQUANT(NLAST(ION))+1 
         IF(NFF.GT.0) FF(ION)=EH/H*IZ(ION)*IZ(ION)/NFF/NFF 
      END DO 
! 
      IF(IPRIN.GE.2) WRITE(6,"(//' EXPLICIT ENERGY LEVELS INCLUDED'/      & 
     &            ' -------------------------------'//                    & 
     & ' NO.    LEVEL    ION   ION.EN.(ERG)        G   NQUANT',           & 
     & '  IEL  ILK  IAT'/)") 
      DO I=1,NLEVEL 
         IF(IPRIN.GE.2)                                                   & 
     &   WRITE(6,"(I4,2X,A10,A4,1PE15.7,0PF10.2,4I5)")                    & 
     &   I,TYPLEV(I),TYPION(IEL(I)),ENION(I),G(I),                        & 
     &                NQUANT(I),IEL(I),ILK(I),IATM(I) 
      END DO 
! 
! ----------------------------------------- 
! Input parameters for additional opacities 
! ----------------------------------------- 
! 
      IF(IPRIN.GE.0) WRITE(6,"(//' ADDITIONAL OPACITY SOURCES'/           & 
     & ' --------------------------'/                                     & 
     & ' IOPHMI  (H-  OPACITY IN LTE)       =',I3/                        & 
     & ' IOPH2P  (H2+  OPACITY)             =',I3/                        & 
     & ' IOPHEM  (HE- B-F AND F-F)          =',I3/                        & 
     & ' IOPCH   (CH OPACITY)               =',I3/                        & 
     & ' IOPOH   (OH OPACITY)               =',I3/                        & 
     & ' IOPH2M  (H2- OPACITY)              =',I3/                        & 
     & ' IOH2H2  (CIA H2-H2 OPACITY         =',I3/                        & 
     & ' IOH2HE  (CIA H2-He OPACITY         =',I3/                        & 
     & ' IOH2H1  (CIA H2-H  OPACITY         =',I3/                        & 
     & ' IOHHE   (CIA H-He OPACITY          =',I3/                        & 
     & ' IRSCT   (RAYLEIGH SCAT. ON H I)    =',I3/                        & 
     & ' IRSCH2  (RAYLEIGH SCAT. ON H2      =',I3/                        & 
     & ' IRSCHE  (RAYLEIGH SCAT. ON HE I)   =',I3/                        & 
     & ' IOPHLI  (LYMAN LINES WINGS)        =',I3)")                      & 
     &           IOPHMI,IOPH2P,IOPHEM,IOPCH,IOPOH,                        & 
     &           IOPH2M,IOH2H2,IOH2HE,IOH2H1,IOHHE,                       & 
     &           IRSCT,IRSCH2,IRSCHE,IOPHLI 
! 
! 
      IF(VTB.LT.1.E3) VTB=VTB*1.E5 
      DO ID=1,ND 
         VTURB(ID)=VTB 
      END DO 
      WRITE(6,"(/' TURBULENT VELOCITY  -  DEPTH-INDEPENDENT  VTURB =',    & 
     &  1PE10.3,'  KM/S'/                                                 & 
     &  ' ------------------'/)") VTB*1.E-5 
      DO I=1,ND 
         VTURB(I)=VTURB(I)*VTURB(I) 
      END DO 
! 
      CALL GETLAL 
      write(*,*) 'after getlal'
 
      RETURN 
      END SUBROUTINE INITIA 
! 
! 
! ************************************************************************ 
! 
! 
! 
      SUBROUTINE RDATA(ION) 
!     ===================== 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: WI1=911.753578, WI2=227.837832,              & 
     &                       T15=1.e-15,ECONST= 5.03411142E15 
      INTEGER,PARAMETER   :: MCFIT=10 
      CHARACTER(LEN=1)    :: A 
      CHARACTER(LEN=1000) :: CADENA 
      CHARACTER(LEN=100)  :: DUM 
      REAL(DP)            :: CTEMP(MCFIT),CRATE(MCFIT) 
 
      DATA IEXP0/0/ 
 
      IUNIT=94 
      OPEN(IUNIT,FILE=FIDATA(ION),STATUS='OLD') 
! 
!     read the first record - a label for the energy level input 
! 
      READ(IUNIT,"(A1)") A 
! 
!   ----------------------------------------------------- 
!   input parameters for explicit energy levels 
!   ----------------------------------------------------- 
! 
!   If ILIMITS(ION) < 0, the program finds out whether energy and 
!   quantum numbers are included in the input data files 
 
      IF(ILIMITS(ION).LT.0) THEN 
         READ(IUNIT,'(1000A)')CADENA 
         BACKSPACE(IUNIT) 
         CALL COUNT_WORDS(CADENA,NOW) 
         IF(NOW.LT.14) THEN 
            ILIMITS(ION)=0 
          ELSE 
            ILIMITS(ION)=1 
         END IF 
      END IF 
 
!   Standard format: ENION(I),G(I),NQUANT(I),TYPLEV(I),ifwop(i) 
 
      IF (ILIMITS(ION).EQ.0) THEN 
! 
         DO IL=1,NLEVS(ION) 
            I=IL+NFIRST(ION)-1 
            IE=IEL(I) 
            N0I=NFIRST(IE) 
            NKI=NNEXT(IE) 
            ia=numat(iatm(n0i)) 
            if(isemex(ia).le.1) then 
               iexp0=iexp0+1 
               iexpl(i)=iexp0 
               iltot(iexp0)=i 
               if(il.eq.nlevs(ion)) then 
                  if(nki.eq.nka(iatm(i))) then 
                     iexp0=iexp0+1 
                     iexpl(nki)=iexp0 
                     iltot(iexp0)=nki 
                  end if 
               end if 
            end if 
            IQ=I-N0I+1 
            X=IQ*IQ 
            ifwop(i)=0 
            IZZ=IZ(IE) 
            READ(IUNIT,*)                                                 & 
     &         ENION(I),G(I),NQUANT(I),TYPLEV(I),ifwop(i) 
            if(ifwop(i).lt.0.and.i.ne.nlast(ie))                          & 
     &         call quit('conflict in negative ifwop') 
            if(ifwop(i).ge.2) ifwop(i)=0 
            IF(I.LT.NKI) THEN 
               E=ENION(I) 
               E0=E 
               IF(E.LT.0.) THEN 
                  E=-E 
                  E0=E 
               END IF 
               IF(E.EQ.0.) THEN 
!                 if(izz.le.2) then 
                  if(izz.le.-2) then 
                  w0=wi1 
                  if(izz.eq.2) w0=wi2 
                  WL0=W0*X 
                  IF(WL0.GT.2000.) THEN 
                     ALM=1.E8/(WL0*WL0) 
                     XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM) 
                     WL0=WL0/(XN1*1.e-6+1.e0) 
                  END IF 
                  E0=H*CL*1.e8/WL0 
                  else 
                  E0=EH*IZZ*IZZ/X 
                  end if 
               END IF 
               IF(E.GT.1.D-7.AND.E.LT.100.) E0=1.6018e-12*E 
               IF(E.GT.100..AND.E.LT.1.e7) E0=1.9857e-16*E 
               IF(E.GT.1.e7) E0=H*E 
               IF(ENION(I).GE.0.) THEN 
                 ENION(I)=E0 
               ELSE 
                 ENION(I)=-E0 
               ENDIF 
               IF(G(I).EQ.0.) G(I)=2.*X 
               IF(NQUANT(I).EQ.0) NQUANT(I)=IQ 
             ELSE 
               IF(G(I).EQ.0..AND.NKI.EQ.NKA(IATM(I))) G(I)=1. 
            END IF 
            if(ifwop(i).lt.0) then 
               enion(i)=0. 
               ff(ie)=0. 
               NMER=NMER+1 
               IMER=NMER 
               IMRG(I)=IMER 
               IIMER(IMER)=I 
            endif 
            fropc(i)=0. 
         END DO 
 
!     Upgraded format including limits for energies, and quantum numbers 
 
      ELSE 
 
         DO IL=1,NLEVS(ION) 
            I=IL+NFIRST(ION)-1 
            IE=IEL(I) 
            N0I=NFIRST(IE) 
            NKI=NNEXT(IE) 
            ia=numat(iatm(n0i)) 
            if(isemex(ia).le.1) then 
               iexp0=iexp0+1 
               iexpl(i)=iexp0 
               iltot(iexp0)=i 
               if(il.eq.nlevs(ion)) then 
                  if(nki.eq.nka(iatm(i))) then 
                     iexp0=iexp0+1 
                     iexpl(nki)=iexp0 
                     iltot(iexp0)=nki 
                  end if 
               end if 
            end if 
            IQ=I-N0I+1 
            X=IQ*IQ 
            ifwop(i)=0 
            IZZ=IZ(IE) 
            READ(IUNIT,*)                                                 & 
     &      ENION(I),G(I),NQUANT(I),TYPLEV(I),ifwop(i),frdodf,imodl,      & 
     &      ENION1(I),ENION2(I),                                          & 
     &      SQUANT1(I),SQUANT2(I),                                        & 
     &      LQUANT1(I),LQUANT2(I),                                        & 
     &      PQUANT1(I),PQUANT2(I) 
            if(ifwop(i).lt.0.and.i.ne.nlast(ie))                          & 
     &         call quit('conflict in negative ifwop') 
            if(ifwop(i).ge.2) ifwop(i)=0 
            IF(I.LT.NKI) THEN 
 
!              check and, if necessary, transform ENION(I) 
 
               E=ENION(I) 
               E0=E 
               IF(E.LT.0.) THEN 
                  E=-E 
                  E0=E 
               END IF 
               IF(E.EQ.0.) THEN 
!                 if(izz.le.2) then 
                  if(izz.le.-2) then 
                  w0=wi1 
                  if(izz.eq.2) w0=wi2 
                  WL0=W0*X 
                  IF(WL0.GT.2000.) THEN 
                     ALM=1.E8/(WL0*WL0) 
                     XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM) 
                     WL0=WL0/(XN1*1.e-6+1.e0) 
                  END IF 
                  E0=H*CL*1.e8/WL0 
                  else 
                  E0=EH*IZZ*IZZ/X 
                  end if 
               END IF 
               IF(E.GT.1.e-7.AND.E.LT.100.) E0=1.6018e-12*E 
               IF(E.GT.100..AND.E.LT.1.e7) E0=1.9857e-16*E 
               IF(E.GT.1.e7) E0=H*E 
               IF(ENION(I).GE.0.) THEN 
                 ENION(I)=E0 
               ELSE 
                 ENION(I)=-E0 
               ENDIF 
 
!           check and, if necessary, transform ENION1(I) 
 
               E=ENION1(I) 
               E0=E 
               IF(E.LT.0.) THEN 
                  E=-E 
                  E0=E 
               END IF 
               IF(E.EQ.0.) THEN 
!                 if(izz.le.2) then 
                  if(izz.le.-2) then 
                  w0=wi1 
                  if(izz.eq.2) w0=wi2 
                  WL0=W0*X 
                  IF(WL0.GT.2000.) THEN 
                     ALM=1.E8/(WL0*WL0) 
                     XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM) 
                     WL0=WL0/(XN1*1.e-6+1.e0) 
                  END IF 
                  E0=H*CL*1.e8/WL0 
                  else 
                  E0=EH*IZZ*IZZ/X 
                  end if 
               END IF 
               IF(E.GT.1.D-7.AND.E.LT.100.) E0=1.6018e-12*E 
               IF(E.GT.100..AND.E.LT.1.e7) E0=1.9857e-16*E 
               IF(E.GT.1.e7) E0=H*E 
               IF(ENION1(I).GE.0.) THEN 
                 ENION1(I)=E0 
               ELSE 
                 ENION1(I)=-E0 
               ENDIF 
 
!           check and, if necessary, transform ENION2(I) 
 
               E=ENION2(I) 
               E0=E 
               IF(E.LT.0.) THEN 
                  E=-E 
                  E0=E 
               END IF 
               IF(E.EQ.0.) THEN 
!              if(izz.le.2) then 
                  if(izz.le.-2) then 
                  w0=wi1 
                  if(izz.eq.2) w0=wi2 
                  WL0=W0*X 
                  IF(WL0.GT.2000.) THEN 
                     ALM=1.E8/(WL0*WL0) 
                     XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM) 
                     WL0=WL0/(XN1*1.e-6+1.e0) 
                  END IF 
                  E0=H*CL*1.e8/WL0 
                  else 
                  E0=EH*IZZ*IZZ/X 
                  end if 
               END IF 
               IF(E.GT.1.D-7.AND.E.LT.100.) E0=1.6018e-12*E 
               IF(E.GT.100..AND.E.LT.1.e7) E0=1.9857e-16*E 
               IF(E.GT.1.e7) E0=H*E 
               IF(ENION2(I).GE.0.) THEN 
                 ENION2(I)=E0 
                ELSE 
                 ENION2(I)=-E0 
               ENDIF 
 
! 
!              Enforce an energy tolerance of 10% when the input files 
!              do not have  any (e.g. pure levels in MODION models) 
! 
               IF((ENION1(I)-ENION(I))/ENION(I).LT.1e-6)                  & 
     &         ENION1(I)=ENION(I)*(1.+ERANGE) 
               IF((ENION(I)-ENION2(I))/ENION(I).LT.1e-6)                  & 
     &         ENION2(I)=ENION(I)*(1.-ERANGE) 
 
! 
!              Convert ENION1,ENION2 to cm-1 from the ground level 
!              so they can be directly used in NLTSET 
! 
 
               ENION1(I)=(ENION(N0I)-ENION1(I))*ECONST 
               ENION2(I)=(ENION(N0I)-ENION2(I))*ECONST 
 
 
               IF(G(I).EQ.0.) G(I)=2.*X 
               IF(NQUANT(I).EQ.0) NQUANT(I)=IQ 
             ELSE 
               IF(G(I).EQ.0..AND.NKI.EQ.NKA(IATM(I))) G(I)=1. 
            END IF 
            if(ifwop(i).lt.0) then 
               write(*,*)'RDATA:  IFWOP<0 and ILIMITS is not 0' 
               stop 
               enion(i)=0. 
               ff(ie)=0. 
               IMER=IMER+1 
               IMRG(I)=IMER 
               IIMER(IMER)=I 
            endif 
            fropc(i)=0. 
         END DO 
 
      END IF 
 
! 
! ---------------------------------------------------------------------- 
! 
!   skip lines if more levels than needed, and skip the continuum transition 
!   label 
! 
      DO 
         READ(IUNIT,"(A1)") A 
         IF(A.EQ.'*') EXIT 
      END DO 
      II0=NFIRST(ION)-1 
      ILLIM=NLLIM(ION)+II0 
      JCORR=0 
! 
!   ----------------------------------------------------- 
!   input parameters for continuum transitions 
!   ----------------------------------------------------- 
! 
      CONT: DO 
         READ(IUNIT,'(A100)',IOSTAT=IOS) DUM 
         IF(IOS.NE.0) EXIT CONT 
         READ(DUM,*,IOSTAT=KSTAT) II,JJ,MODE,                             & 
     &      IFANCY,ICOLIS,                                                & 
     &      IFRQ0,IFRQ1,OSC,CPARAM,NCOL 
         IF(KSTAT.NE.0) THEN 
            READ(DUM,*,IOSTAT=IOS2) II,JJ,MODE,IFANCY,ICOLIS,             & 
     &      IFRQ0,IFRQ1,OSC,CPARAM 
            IF(IOS2.NE.0) EXIT CONT 
            NCOL=0 
         END IF 
         IF (NCOL.NE.0) THEN 
            DO IIC=1,NCOL 
               READ(IUNIT,*) ITYPE, NCTEMP 
               READ(IUNIT,*) (CTEMP(IFIT),IFIT=1,NCTEMP) 
               READ(IUNIT,*) (CRATE(IFIT),IFIT=1,NCTEMP) 
            END DO 
         END IF 
! 
         IF(II.EQ.0) THEN 
            IF(JJ.EQ.0) THEN 
               CLOSE(IUNIT) 
               RETURN 
            END IF 
            II0=JJ-1 
            CYCLE CONT 
         END IF 
         IF(IABS(MODE).GT.100) READ(IUNIT,*) FR0INP 
         if(iabs(mode).eq.2) then 
            READ(IUNIT,*) kdo 
            CYCLE CONT 
         end if 
         IF(IFANCY.GT.49.and.ifancy.lt.100) IASV=1 
         if(iabs(mode).eq.3.or.iabs(mode).eq.4) cycle cont 
         IF(IABS(MODE).EQ.5 .OR. IABS(MODE).EQ.15) THEN 
            READ(IUNIT,*) FROPCI 
            if(ion.eq.ielh) then 
               if(ii.eq.1.and.cutlym.ne.0) fropci=-cutlym 
               if(ii.eq.2.and.cutbal.ne.0) fropci=-cutbal 
            end if 
            if(abs(fropci).lt.1.e10) fropci=2.997925e18/fropci 
         END IF 
         IF(II.EQ.1) JCORR=NLEVS(ION)+1-JJ 
         II=II+II0 
         JJ=JJ+II0+JCORR 
         FROPC(II)=FROPCI 
         N0I=NFIRST(IE) 
         NKI=NNEXT(IE) 
         IF(JJ.GE.NKI) THEN 
         LPC=.FALSE. 
         IF(IELHE2.GE.0) THEN 
            IF(II.GE.NFIRST(IELHE2).AND.II.LE.NLAST(IELHE2)               & 
     &     .AND.IFWOP(II).GE.0) LPC=.TRUE. 
         END IF 
         IF(II.GE.N0HN.AND.II.LE.N1H.AND.IFWOP(II).GE.0) LPC=.TRUE. 
         IF(LPC) THEN 
            MODE=5 
            XI=NQUANT(II) 
            X2=XI+3. 
            if(ii.ge.8) x2=xi+2. 
            IF(FROPC(II).GE.0.) THEN 
               FROPC(II)=ENION(II)/6.6256E-27*(1.-XI*XI/(X2*X2)) 
             ELSE 
               FROPC(II)=ABS(FROPC(II)) 
            END IF 
         END IF 
         END IF 
         IF(MODE.EQ.0) THEN 
            IF(II.LT.NLAST(ION)) CYCLE CONT 
            IF(II.EQ.NLAST(ION)) EXIT CONT 
         END IF 
! 
!   ----------------------------------------------------- 
!   Additional input parameters for continuum transitions 
!   ----------------------------------------------------- 
! 
!        Only for IFANCY = 2, 3, or 4 
!        S0BF, ALFBF, BETBF, GAMBF  - parameters for evaluation the 
!        photoionization cross-section 
! 
         IF(IFANCY.GE.2.AND.IFANCY.LE.4)                                  & 
     &   READ(IUNIT,*) S0BF(II),ALFBF(II),BETBF(II),GAMBF(II) 
! 
!   ----------------------------------------------------- 
!   Additional input parameters for continuum transitions -TOPBASE DATA 
!   ----------------------------------------------------- 
! 
!        Only for IFANCY > 100 there are IFANCY-100 fit points 
! 
!        XTOP(MFIT,MCROSS) -  x = alog10(nu/nu0) of a fit point 
!        CTOP(MFIT,MCROSS) -  sigma = alog10(sigma/10^-18) of a fit point 
! 
!        there are IFANCY-100 fit points 
! 
         IF(IFANCY.GT.100) THEN 
            NFIT=IFANCY-100 
            IF(NFIT.GT.MFIT) call quit(' nfit too large (TOPBASE fits)') 
            READ(IUNIT,*) (XTOP(IFIT,II),IFIT=1,NFIT) 
            READ(IUNIT,*) (CTOP(IFIT,II),IFIT=1,NFIT) 
         END IF 
         IBF(II)=IFANCY 
         INDEXP(II)=IABS(MODE) 
         IF(II.LT.NLAST(ION)) THEN 
            CYCLE CONT 
          ELSE 
            EXIT CONT 
         END IF 
      END DO CONT 
 
      DO 
        READ(IUNIT,"(A1)",IOSTAT=IOS) A 
        IF(IOS.NE.0) THEN 
           CLOSE(IUNIT) 
           RETURN 
        END IF 
        IF(A.EQ.'*') EXIT 
      END DO 
! 
!  ----------------------------------------------------------- 
!  Input parameters for line transitions 
!  ----------------------------------------------------------- 
! 
      LINES: DO 
         READ(IUNIT,*,IOSTAT=IOS) II,JJ,MODE,IFANCY,ICOLIS,               & 
     &                            IFRQ0,IFRQ1,OSC,CPARAM 
         IF(IOS.NE.0) EXIT LINES 
         IF(IABS(MODE).GT.100) READ(IUNIT,*) FR0INP 
         IF(JJ.GT.NLEVS(ION)) THEN 
            IF(IABS(MODE).EQ.2) THEN 
               READ(IUNIT,*) K1,K2,K3,X1,X2,X3,K4 
               CYCLE LINES 
            END IF 
            IF(IABS(MODE).EQ.1) READ(IUNIT,*) LCMP 
            IF(IABS(IFANCY).EQ.1) READ(IUNIT,*) GAMR,STARK1,STARK2,       & 
     &         STARK3,VDWH 
            CYCLE LINES 
         END IF 
         if(iabs(mode).eq.2) then 
            READ(IUNIT,*) K1,K2,K3,X1,X2,X3,K4 
            cycle lines 
         end if 
         if(iabs(mode).eq.3.or.iabs(mode).eq.4) cycle lines 
         IF(MODE.EQ.0) CYCLE LINES 
! 
!  ----------------------------------------------------------- 
!  Additional input parameters for "clasical" line transitions 
!   (i.e. those not represented by ODF's - ie ABS(MODE)=1) 
!  ----------------------------------------------------------- 
! 
         READ(IUNIT,*) LCOMP,INTMOD,NF,XMAX,TSTD 
         IF(IABS(IFANCY).EQ.1) READ(IUNIT,*) GAMR,STARK1,STARK2,          & 
     &       STARK3,VDWH 
      END DO LINES 
! 
      CLOSE(IUNIT) 
      RETURN 
      END SUBROUTINE RDATA 
! 
! 
!    ***************************************************************** 
! 
! 
      SUBROUTINE NSTPAR(FINSTD) 
!     ========================== 
! 
!     settiing up the default values of various input flags, and 
!     input of non-standard values of various input flags and parameters 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      INTEGER, PARAMETER :: MVAR=77 
      INTEGER, PARAMETER :: INPFI=4 
      CHARACTER*(*) FINSTD 
      CHARACTER(LEN=80) :: TEXT 
      CHARACTER(LEN=6)  :: PVALUE(MVAR),VARNAM(MVAR) 
      CHARACTER(LEN=20) :: BLNK 
      CHARACTER(LEN=6)  :: BLNK6 
! 
      DATA VARNAM /'IATREF',                                              & 
     &             'BERGFC','IHYDPR','NUNHHE','STHE  ',                   & 
     &             'NFREQS','IBFAC ',                                     & 
     &             'INTRPL','ICHANG','IFEOS ',                            & 
     &             'IOPHMI','IOPH2P','IOPHEM','IOPCH ','IOPOH ',          & 
     &             'IOPH2M','IOH2H2','IOH2HE','IOH2H1','IOHHE ',          & 
     &             'IRSCT ','IRSCH2','IRSCHE',                            & 
     &             'TRAD  ','WDIL  ',                                     & 
     &             'VTB   ','IFMOL ','TMOLIM',                            & 
     &             'MOLTAB','IRWTAB','IIRWIN','IPFEXO',                   & 
     &             'IPFBC ','IEQBC ',                                     & 
     &             'CUTLYM','CUTBAL','IHXENB',                            & 
     &             'GSSTD ','GWSTD ',                                     & 
     &             'IHGOM ','HGLIM ',                                     & 
     &             'ERANGE',                                              & 
     &             'ISPICK','ILPICK','IPPICK',                            & 
     &             'IHYDDK',                                              & 
!-------------------------------------------------------------------------- 
     &             'IMODE ','IDSTD ','IPRIN ',                            & 
     &             'INMOD ','ICHEMC',                                     & 
     &             'IOPHLI','NUNALP','NUNBET','NUNGAM','NUNBAL',          & 
     &             'IFREQ ','INLTE ','ICONTL','INLIST','IFHE2 ',          & 
     &             'IHYDPR','IHE1PR','IHE2PR',                            & 
     &             'ALAM0 ','ALAST ','CUTOF0','CUTOFS','RELOP',           & 
     &             'SPACE',                                               & 
     &             'NMLIST','IUNIM1','IUNIM2',                            & 
     &             'VTB   ',                                              & 
     &             'NMU0  ','ANG0  ','IFLUX '/ 
 
! 
      DATA PVALUE /'     1',                                              & 
     &             '  1.e0','     0','     0',' 1.e19',                   & 
     &             '   120','     0',                                     & 
     &             '     0','     0','     0',                            & 
     &             '     1','     1','     1','     1','     1',          & 
     &             '     1','     1','     1','     1','     1',          & 
     &             '     1','     1','     1',                            & 
     &             '    0.','    0.',                                     & 
     &             '    2.','     1',' 9000.',                            & 
     &             '     1','     1','     1','     1',                   & 
     &             '     1','     1',                                     & 
     &             '    0.','    0.','     0',                            & 
     &             '3.1e-5','1.0e-7',                                     & 
     &             '     0',' 1.e18',                                     & 
     &             '  0.10',                                              & 
     &             '     1','     1','     1',                            & 
     &             '     0',                                              & 
!-------------------------------------------------------------------------- 
     &             '     0','     0','     0',                            & 
     &             '     1','     0',                                     & 
     &             '     0','     0','     0','     0','     0',          & 
     &             '     1','     1','     0','     0','     0',          & 
     &             '     0','     0','     0',                            & 
     &             '     0',' 10000','    20','     0',' 1.e-4',          & 
     &             ' 1.e-2',                                              & 
     &             '     0','    20','    21',                            & 
     &             '    2.',                                              & 
     &             '     0','   0.1','     0'/ 
! 
      DATA BLNK/'                    '/,BLNK6/'      '/ 
! 
      IF(FINSTD.NE.BLNK) THEN 
        OPEN(UNIT=INPFI,FILE=FINSTD,STATUS='UNKNOWN') 
        write(6,"(/'CONTENTS OF KEYWORD FILE',3x,a8/                      & 
     &          '------------------------')") FINSTD 
      END IF 
! 
      INDV=-1 
! 
!    go through the input file line by line 
! 
!     write(6,"(/' INPUT KEYWORD PARAMETERS, FROM FILE:',3x,a8/           & 
!    &        ' -------------------------')") FINSTD 
! 
!     write(6,"(/'CONTENTS OF KEYWORD FILE',3x,a8/                        & 
!    &          '------------------------')") FINSTD 
      READLINES: DO 
         K0=1 
         READ(INPFI,"(A)",IOSTAT=IOS) TEXT 
         IF(IOS.NE.0) EXIT READLINES 
         WRITE(6,*) TEXT 
         IN: DO 
            CALL GETWRD(TEXT,K0,K1,K2) 
            IF(K1.EQ.0) CYCLE READLINES 
            K0=K2+2 
            IF(TEXT(K1:K2).EQ.'=') CYCLE IN 
            INDV=-INDV 
            IF(INDV.EQ.1) THEN 
               DO I=1,MVAR 
                  IF(TEXT(K1:K2).EQ.VARNAM(I)(1:K2-K1+1)) THEN 
                     IVAR=I 
                     PVALUE(IVAR)=BLNK6 
                     PVALUE(IVAR)(6-K2+K1:6)=TEXT(K1:K2) 
                     CYCLE IN 
                  END IF 
               END DO 
               CALL GETWRD(TEXT,K0,K1,K2) 
               IF(K1.EQ.0) THEN 
                  K0=1 
                  IN2: DO 
                     READ(INPFI,"(A)",IOSTAT=IOS) TEXT 
                     IF(IOS.NE.0) EXIT READLINES 
                     CALL GETWRD(TEXT,K0,K1,K2) 
                     IF(K1.EQ.0) CYCLE IN2 
                  END DO IN2 
               END IF 
               K0=K2+2 
               INDV=-INDV 
               CYCLE IN 
               IVAR=I 
             ELSE 
               PVALUE(IVAR)=BLNK6 
               PVALUE(IVAR)(6-K2+K1:6)=TEXT(K1:K2) 
            END IF 
         END DO IN 
      END DO READLINES 
! 
      DO I=1,MVAR 
         WRITE(84,"(1X,A)") PVALUE(I) 
      END DO 
! 
      CLOSE(UNIT=84) 
      REWIND(84) 
      READ(84,*)                                                          & 
     &             IATREF,                                                & 
     &             BERGFC,IHYDPR,NUNHHE,STHE  ,                           & 
     &             NFREQS,IBFAC ,                                         & 
     &             INTRPL,ICHANG,IFEOS ,                                  & 
     &             IOPHMI,IOPH2P,IOPHEM,IOPCH ,IOPOH ,                    & 
     &             IOPH2M,IOH2H2,IOH2HE,IOH2H1,IOHHE ,                    & 
     &             IRSCT ,IRSCH2,IRSCHE,                                  & 
     &             TRAD  ,WDIL  ,                                         & 
     &             VTB   ,IFMOL ,TMOLIM,                                  & 
     &             MOLTAB,IRWTAB,IIRWIN,IPFEXO,                           & 
     &             IPFBC ,IEQBC ,                                         & 
     &             CUTLYM,CUTBAL,IHXENB,                                  & 
     &             GSSTD ,GWSTD ,                                         & 
     &             IHGOM ,HGLIM ,                                         & 
     &             ERANGE,                                                & 
     &             ISPICK,ILPICK,IPPICK,                                  & 
     &             IHYDDK 
!-------------------------------------------------------------------------- 
      IF(IF55.EQ.0)READ(84,*)                                             & 
     &             IMODE ,IDSTD ,IPRIN ,                                  & 
     &             INMOD ,ICHEMC,                                         & 
     &             IOPHLI,NUNALP,NUNBET,NUNGAM,NUNBAL,                    & 
     &             IFREQ ,INLTE ,ICONTL,INLIST,IFHE2 ,                    & 
     &             IHYDPR,IHE1PR,IHE2PR,                                  & 
     &             ALAM0 ,ALAST ,CUTOF0,CUTOFS,RELOP,SPACE,               & 
     &             NMLIST,IUNIM1,IUNIM2,                                  & 
     &             VTB   ,                                                & 
     &             NMU0  ,ANG0  ,IFLUX 
! 
      close(84) 
 
      if(imode.le.-3) then 
         irsct=0 
         irsche=0 
         irsch2=0 
      end if 
 
      if(imode.ge.0.and.teff.gt.12000.) ifmol=0 
 
      if(ifmol.eq.0) then 
         ioh2h2=0 
         ioh2he=0 
         ioh2h1=0 
         iohhe =0 
         ioph2m=0 
         iopoh =0 
         iopch =0 
      end if 
 
      write(6,"(/'IFMOL = ',i2)") ifmol 
 
      NMLIS0=NMLIST 
! 
      RETURN 
      END SUBROUTINE NSTPAR 
! 
! 
!    *************************************************************** 
! 
! 
        subroutine count_words(cadena,n) 
! 
!       Counts the number of words separated by blanks in a string 
! 
        character(len=1000) :: cadena 
        character(len=1   ) :: a,b 
 
        n=0 
        a=cadena(1:1) 
        if (a.ne.' ') n=1 
        do i=2,len(cadena) 
           b=cadena(i:i) 
           if(b.ne.' '.and.a.eq.' ') n=n+1 
           a=b 
        end do 
        end subroutine count_words 
! 
! 
!    *************************************************************** 
! 
! 
      SUBROUTINE GETWRD(TEXT,K0,K1,K2) 
! 
!  FINDS NEXT WORD IN TEXT FROM INDEX K0. NEXT WORD IS TEXT(K1:K2) 
!  THE NEXT WORD STARTS AT THE FIRST ALPHANUMERIC CHARACTER AT K0 
!  OR AFTER. IT ENDS WITH THE LAST ALPHANUMERIC CHARACTER IN A ROW 
!  FROM THE START 
! 
!  TAKEN FROM MULTI - M. CARLSSON (1976) 
! 
      INTEGER, PARAMETER :: MSEPAR=7 
      LOGICAL            :: LSEP,LSEPE 
      CHARACTER*(*) TEXT 
      CHARACTER(LEN=1)   :: SEPAR(MSEPAR) 
      DATA SEPAR/' ','(',')','=','*','/',','/ 
! 
      K1=0 
      ILOOP: DO I=K0,LEN(TEXT) 
         LSEP=.TRUE. 
         IF(K1.EQ.0) THEN 
            JLOOP: DO J=1,MSEPAR 
               IF(TEXT(I:I).EQ.SEPAR(J)) THEN 
                  LSEP=.FALSE. 
                  EXIT JLOOP 
               END IF 
            END DO JLOOP 
 
            IF(LSEP) K1=I 
! 
!           NOT START OF WORD 
! 
         ELSE 
            LSEPE=.TRUE. 
            DO J=1,MSEPAR 
               IF(TEXT(I:I).EQ.SEPAR(J)) THEN 
                  LSEPE=.FALSE. 
                  EXIT ILOOP 
               END IF 
            END DO 
        END IF 
      END DO ILOOP 
! 
!  NO NEW WORD. RETURN K1=K2=0 
! 
      IF(LSEPE) THEN 
         K1=0 
         K2=0 
         RETURN 
      END IF 
! 
!  NEW WORD IN TEXT(K1:I-1) 
! 
      K2=I-1 
! 
      RETURN 
      END 
! 
! 
! ******************************************************************* 
! 
! 
      subroutine readkw 
!     ================= 
! 
!     reads keyword parameters, first from file 'kwords' if it exists 
!     and then form 'fort.55', if it exists 
! 
      use params 
      implicit none 
      integer :: ios,iosk 
 
      open(unit=55,file='fort.55',status='old',iostat=ios) 
      if(ios.eq.0) then 
         if55=1 
         write(6,"('FILE fort.55 EXISTS')") 
       else 
         if55=0 
         write(6,"('FILE fort.55 DOES NOT EXIST')") 
      end if 
      if(if55.eq.1) then 
         read(55,*,iostat=ios) imode 
         if(ios.ne.0) then 
            write(6,"('FILE fort.55 IS HOWEVER EMPTY')") 
            if55=0 
         end if 
         rewind(55) 
      end if 
! 
      open(unit=56,file='kwords',status='old',iostat=iosk) 
      if(iosk.eq.0) then 
         ifkey=1 
         write(6,"('FILE kwords  EXISTS')") 
       else 
         ifkey=0 
         write(6,"('FILE kwords  DOES NOT EXIST')") 
      end if 
! 
      if(ifkey.eq.1) call nstpar('kwords') 
      if(if55.eq.1) call read55 

      call ldimen 
      call readnd 
      return 
      end subroutine readkw 
! 
! 
! 
!     ******************************************************************* 
! 
! 
      subroutine read55 
!     ================= 
! 
!     read (optional) input filw fort.55 
! 
      use accura 
      use params 
      implicit none 
      integer :: ios 
 
      igrdd=0 
      IFWIN=0 
      NDSTEP=0 
      if(ifeos.gt.0) return 
 
!!    open(iunit=55,file='fort.55',status='old',iostat=ios) 
!!    if(ios.ne.0) return 
 
      READ(55,*,iostat=ios) IMODE,IDSTD,IPRIN 
      IF(IOS.NE.0) THEN 
         WRITE(6,"(/'BUT FORT.55 IS EMPTY'/)") 
         IF55=0 
         RETURN 
      END IF 
 
      READ(55,*) INMOD,INTRPL,ICHANG,ICHEMC 
      READ(55,*) IOPHLI,nunalp,nunbet,nungam,nunbal 
      IOPHLI=0 
 
      READ(55,*) IFREQ,INLTE,ICONTL,INLIST,IFHE2 
      READ(55,*) IHYDPR,IHE1PR,IHE2PR 
      READ(55,*) ALAM0,ALAST,CUTOF0,CUTOFS,RELOP,SPACE 
 
      WRITE(6,"(/'SEVERAL KEYWORD PARAMS READ FROM FORT.55:'/              & 
     &           '-----------------------------------------')") 
      WRITE(6,"('IMODE,IDSTD,IPRIN    = ',3i4)") IMODE,IDSTD,IPRIN 
      WRITE(6,"('INMOD,ICHEMC         = ',2i4)") INMOD,ICHEMC 
      WRITE(6,"('IHYDPR,IHE1PR,IHE2PR = ',3i4)") IHYDPR,IHE1PR,IHE2PR 
      WRITE(6,"('IFREQ,INLTE,INLIST   = ',3i4)") IFREQ,INLTE,INLIST 
      WRITE(6,"('ALAM0                = ',f15.5/                           & 
     &          'ALAST                = ',f15.5/                           & 
     &          'CUTOF0               = ',f15.5/                           & 
     &          'RELOP                = ',f15.5/                           & 
     &          'SPACE                = ',f15.5)")                         & 
     &           ALAM0,ALAST,CUTOF0,RELOP,SPACE 
      return 
      end subroutine read55 
!
!
!     ************************************************************
!
!
      subroutine ldimen
!     =================
!
!     setting dimensions for dealing with lines
!
      use params
      use lindat
      use molist

      character(len=80) :: a,b,c

      mlin0  = 1000000
      mlin   =   10000
      mlinm0 = 1000000
      mlinm  =   50000
      msftab = 6000000
      if(imode.ge.-2) msftab=1

      write(*,"(/'DEFAULT DIMENSIONS OF LARGE ARRAYS'/                    &
     &           '----------------------------------'/                    &
     &           'MLIN0   ',i10,'  (atomic lines stored)'/                &
     &           'MLIN    ',i10,'  (atomic lines in in one set)'/         &
     &           'MLINM0  ',i10,'  (molecular lines per list stored)'/    &
     &           'MLINM   ',i10,'  (molecular lines in one set)'/         &
     &           'MSFTAB  ',i10,'  (number of internal freq.points)'      &
     &           /)") mlin0,mlin,mlinm0,mlinm,msftab

      open(unit=57,file='dimens',status='old',iostat=ios)
      if(ios.ne.0) return

      readl: do
         read(57,"(a80)",iostat=ios2) a
         if(ios2.ne.0) exit readl
         i=index(a,'=')
         b=trim(a(1:i-1))
         c=trim(a(i+1:80))
         read(c,*) idim
         if(b.eq.'mlin'.or.b.eq.'MLIN') mlin=idim
         if(b.eq.'mlin0'.or.b.eq.'MLIN0') mlin0=idim
         if(b.eq.'mlinm'.or.b.eq.'MLINM') mlinm=idim
         if(b.eq.'mlinm0'.or.b.eq.'MLINM0') mlinm0=idim
         if(b.eq.'msftab'.or.b.eq.'MSFTAB') msftab=idim
      end do readl


      write(*,"(/'ACTUAL DIMENSIONS OF LARGE ARRAYS'/                     &
     &           '---------------------------------'/                     &
     &           'MLIN0   ',i10,'  (atomic lines stored)'/                &
     &           'MLIN    ',i10,'  (atomic lines in in one set)'/         &
     &           'MLINM   ',i10,'  (molecular lines per list stored)'/    &
     &           'MLINM0  ',i10,'  (molecular lines in one set)'/         &
     &           'MSFTAB  ',i10,'  (number of internal freq.points)'      &
     &           /)") mlin0,mlin,mlinm0,mlinm,msftab

      close(57)
      return
      end subroutine ldimen
! 
! 
! ******************************************************************* 
! 
! 
      subroutine readnd 
!     ================= 
! 
      use params 
      implicit none 
      integer :: npar,nt,ios,ndeq1 
 
      ndeq1=0 
      if(imode.le.-3) then 
         read(2,*) nt 
         rewind(2) 
         if(nt.gt.0) ndeq1=1 
       else if(imode.eq.-1) then 
         ndeq1=1 
      end if 
 
      if(ndeq1.eq.1) then 
         nd=1 
       else 
         read(8,*,iostat=ios) nd,npar 
         if(ios.eq.0) then 
            if(inmod.eq.2) inmod=1 
            write(6,"(/'INPUT MODEL tlusty')") 
          else 
            rewind(8) 
            inmod=0 
            write(6,"(/'INPUT MODEL Kurucz')") 
            READ(8,"(//////////////////////10X,I3/)") ND 
            ND=ND-1 
         end if 
         rewind(8) 
      end if 
 
      mdepth=nd 
      write(*,"(/'ND = MDEPTH =',i4/)") nd 
      return 
      end subroutine readnd 
 
!
!
!     ************************************************************
!
!

      subroutine prkeyw
!     =================

      use params

      write(6,"('IMPORTANT KEYWORDS PARAMETERS'/                            &
     &          '-----------------------------')")
      write(6,"('IMODE, IDSTD, IPRIN, IFEOS,IFMOL  =',5i5)")                &
     &           IMODE, IDSTD, IPRIN, IFEOS,IFMOL
      write(6,"('IRWTAB,IIRWIN,IPFEXO,IPFBC,IEQBC  =',5i5)")                &
     &           IRWTAB,IIRWIN,IPFEXO,IPFBC,IEQBC
      write(6,"('IHYDPR,IHE1PR,IHE2PR,IHGOM,IHYDDK =',5i5)")                &
     &           IHYDPR,IHE1PR,IHE2PR,IHGOM,IHYDDK
      write(6,"('NFREQS,IBFAC ,NUNHHE,IFHE2,INLTE  =',5i5/)")               &
     &           NFREQS,IBFAC ,NUNHHE,IFHE2,INLTE

      end subroutine prkeyw
! 
! 
!     **************************************************************** 
! 
! 
      SUBROUTINE STATE0(MODOLD) 
!     ========================= 
! 
!     Initialization of the basic parameters for the Saha equation 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      real(dp), parameter :: enhe1=24.5799,enhe2=54.3999 
      character(len=4)    :: DYP(MATOM) 
      character(len=80)   :: dum 
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
! 
!      DATA XIFE /8*0.,233.6,262.1/ 
!      DATA NTOTA /99/ 
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
      iabu12=0
      read(ibuff,'(a80)') dum 
      read(dum,*,iostat=kstat1) natoms,iabset i,iabu12
      if(kstat1.ne.0) read(dum,*,iostat=kstat2) natoms,iabset 
      if(kstat2.ne.0) READ(dum,*) NATOMS 
      WRITE(6,"(//' CHEMICAL ELEMENTS INCLUDED'/                          & 
     &            ' --------------------------'//                         & 
     & ' NUMBER  ELEMENT           ABUNDANCE'/1H ,16X,                    & 
     & 'A=N(ELEM)/N(H)  A/A(SOLAR)'/)") 
      IAT=0 
      IREF=0 
      IF(NATOMS.LT.0) NATOMS=-NATOMS 
! 
      DO I=1,MATOM 
         DO J=1,MION0 
            RR(I,J)=0. 
         END DO 
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
         LGR(I)=.TRUE. 
         LRM(I)=.TRUE. 
         IATEX(I)=-1 
         IF(I.LE.NATOMS) THEN 
            IF(MODOLD.EQ.0) THEN 
               READ(IBUFF,*) MA,NA0,NAK,ION,MODPF(I),ABN,                 & 
     &                       (PFSTD(J,I),J=1,5) 
               MA=IABS(MA) 
             ELSE 
               READ(IBUFF,*) MA,ABN,MODPF(I) 
               ION=0 
            END IF 
          ELSE IF(MOD(IMODE,10).LE.1.and.imode.ne.-4) THEN 
            MA=1 
            ABN=0. 
            ION=0 
            MODPF(I)=0 
           ELSE 
            MA=0 
         END IF 
         AMAS(I)=D(1,I) 
         ABND(I)=D(2,I) 
         if(iref.gt.0) abnd(i)=d(2,i)*abnd(iref)/d(2,iref) 
         IONIZ(I)=int(D(3,I)) 
         isemex(i)=0 
! 
!        increase the standard highest ionization for Teff>30,000 K 
! 
         IF(TEFF.GT.3.D4) THEN 
           IF(I.LE.8) IONIZ(I)=I+1 
           IF(I.GT.8.and.i.le.30) IONIZ(I)=9 
         END IF 
! 
         DO J=1,9 
            IF(J.LE.8) ENEV(I,J)=xi(J,I) 
            if(enev(i,j).ge.enhe2) then 
               inpot(i,j)=3 
             else if(enev(i,j).ge.enhe1) then 
               inpot(i,j)=2 
             else 
               inpot(i,j)=1 
            end if 
         END DO 
         IF(MA.GT.0) THEN 
            LGR(I)=.FALSE. 
            IF(ABN.GT.0) ABND(I)=ABN 
            IF(ABN.LT.0) ABND(I)=ABS(ABN)*D(2,I) 
            IF(IABU12.NE.0) ABND(I)=10.**(12.-ABN)
            IF(ION.NE.0) IONIZ(I)=ION 
            IF(ABN.GT.1.E6) THEN 
               READ(IBUFF,*) (ABNDD(I,ID),ID=1,ND) 
             ELSE 
               DO ID=1,ND 
                  ABNDD(I,ID)=ABND(I) 
               END DO 
            END IF 
            IF(MA.EQ.1) THEN 
               LRM(I)=.FALSE. 
               IATEX(I)=0 
             ELSE 
               IAT=IAT+1 
               IATEX(I)=IAT 
               if(ma.eq.4) isemex(i)=1 
               if(ma.eq.5) isemex(i)=2 
               IF(IAT.EQ.IATREF) THEN 
                  IREF=I 
                  DO ID=1,ND 
                     ABNREF(ID)=ABNDD(I,ID) 
                  END DO 
               END IF 
! 
!              store parameters for explicit atoms 
! 
               DO ID=1,ND 
                  ABUND(IAT,ID)=ABNDD(I,ID) 
               END DO 
               AMASS(IAT)=AMAS(I)*HMASS 
               NUMAT(IAT)=I 
               IF(MODOLD.EQ.0) THEN 
                  N0A(IAT)=NA0 
                  NKA(IAT)=NAK 
               END IF 
            END IF 
            DO ID=1,ND 
               YTOT(ID)=YTOT(ID)+ABNDD(I,ID) 
               WMY(ID)=WMY(ID)+ABNDD(I,ID)*AMAS(I) 
            END DO 
            ABN=ABND(I)/D(2,I) 
            IF(MA.EQ.1) WRITE(6,"(I4,3X,A5,1P2E14.2)")                    & 
     &         I,TYPAT(I),ABND(I),ABN 
            IF(MA.EQ.2) WRITE(6,"(I4,3X,A5,1P2E14.2,3X,                   & 
     &       'EXPLICIT: IAT=',I3)")                                       & 
     &       I,TYPAT(I),ABND(I),ABN,IAT
         END IF 
      END DO 
      IF(MOD(IMODE,10).LE.1) NATOMS=MATOM 
      DO ID=1,ND 
         WMM(ID)=WMY(ID)*HMASS/YTOT(ID) 
      END DO 
      DO JJ=1,NATOMS 
         DO ID=1,ND 
            RELAB(JJ,ID)=1. 
         END DO 
      END DO 
! 
      IF(ICHEMC.EQ.1) THEN 
! 
!        abundance change with respect to the model atmosphere input 
!       (unit 5); 
!        this option is switched on by the parameter ICHEMC (read from 
!        unit 55), if it is non-zero, an additional input from 
!        unit 56 is required 
! 
!        unit 56 input: 
! 
!        NCHANG  -  number of chemical elements for which the abundances 
!                are going to be changes; 
! 
!        then there are NCHANG records, each contains: 
! 
!        I       - atomic number 
!        ABN     - new abundance; coded using the same conventions as in 
!               the standard input 
! 
         READ(56,*,IOSTAT=IOS) NCHANG 
         IF(IOS.NE.0) THEN 
            WRITE(6,"(//' CHEMICAL COMPOSITION COULD NOT BE READ FROM ',  & 
     &           'UNIT 56'//' STOP.')") 
            STOP 
         END IF 
 
         WRITE(6,"(//'  CHEMICAL ELEMENTS INCLUDED - CHANGED (unit 56)'   & 
     &        /'    --------------------------'//                         & 
     &         ' NUMBER  ELEMENT           ABUNDANCE'/1H ,16X,            & 
     &         'A=N(ELEM)/N(H)  A/A(SOLAR)'/)") 
         DO II=1,NCHANG 
            READ(56,*) I,ABN 
            ABND(I)=D(2,I) 
            IF(ABN.GT.0) ABND(I)=ABN 
            IF(ABN.LT.0) ABND(I)=-ABN*D(2,I) 
!           IF(ABN.LT.0) ABND(I)=-ABN*ABNDD(I,1) 
            if(iabu12.ne.0) abnd(i)=10.**(abn-12.) 
            DO ID=1,ND 
               ABNDD(I,ID)=ABND(I) 
            END DO 
            LGR(I)=.FALSE. 
            IATX=IATEX(I) 
            IF(IATX.GT.0) THEN 
               DO ID=1,ND 
                  RELAB(IATX,ID)=ABNDD(I,ID)/ABUND(IATX,ID) 
                  ABUND(IATX,ID)=ABNDD(I,ID) 
               END DO 
            END IF 
            ABNR=ABND(I)/D(2,I) 
            WRITE(6,"(1X,I4,3X,A5,1P2E14.2)") I,TYPAT(I),ABND(I),ABNR 
         END DO 
      END IF 
! 
!     renormalize abundances to have the standard element abundance 
!     equal to unity 
! 
      IF(IREF.LE.1) RETURN 
      write(6,"(//'    CHEMICAL ELEMENTS INCLUDED - RENORMALIZATION'/     & 
     &            '    --------------------------'//                      & 
     &          ' NUMBER  ELEMENT           ABUNDANCE'/1H ,16X,           & 
     &          'A=N(ELEM)/N(H)  A/A(SOLAR)'/)") 
      DO I=1,MATOM 
         IAT=IATEX(I) 
         IF(IAT.GE.0) THEN 
            DO ID=1,ND 
               ABNDD(I,ID)=ABNDD(I,ID)/ABNREF(ID) 
               YTOT(ID)=YTOT(ID)+ABNDD(I,ID) 
               WMY(ID)=WMY(ID)+ABNDD(I,ID)*AMAS(I) 
            END DO 
            ABNR=ABND(I)/D(2,I) 
            IF(IAT.EQ.0) THEN 
               WRITE(6,"(1X,I4,3X,A5,1P2E14.2)") I,TYPAT(I),ABND(I),ABNR 
             ELSE 
               DO ID=1,ND 
                  ABUND(IAT,ID)=ABNDD(I,ID) 
               END DO 
               WRITE(6,"(1X,I4,3X,A5,1P2E14.2,3X,                         & 
     &           'EXPLICIT: IAT=',I3,'  N0A=',I3,'  NKA=',I3)")           & 
     &           I,TYPAT(I),ABND(I),ABNR,IAT,N0A(IAT),NKA(IAT) 
            END IF 
         END IF 
      END DO 
      DO ID=1,ND 
         WMM(ID)=WMY(ID)*HMASS/YTOT(ID) 
      END DO 
      RETURN 
      END SUBROUTINE STATE0 
! 
! 
!     **************************************************************** 
! 
! 
      SUBROUTINE INIMOD 
! 
!   VALUES OF  N(ION)/U(ION) FOR ALL THE ATOMS 
!   AND IONS CONSIDERED 
! 
      use accura 
      use params 
      use modelp 
      use lindat 
      implicit real(dp) (a-h,o-z),logical (l) 
! 
!     1. "low-temperature" ionization fractions 
!         (using Hamburg partition functions) 
! 
      DEPTHS: DO ID=1,ND 
         IF(IFMOL.EQ.0.OR.TEMP(ID).GE.TMOLIM) THEN 
            CALL STATE(ID,TEMP(ID),ELEC(ID),S1) 
            HPOP=DENS(ID)/WMM(ID)/YTOT(ID) 
            DO J=1,MION0 
               DO I=1,MATOM 
                  RRR(ID,J,I)=RR(I,J)*HPOP 
               END DO 
            END DO 
            DO IAT=1,NATOM 
               ATTOT(IAT,ID)=HPOP*ABUND(IAT,ID) 
            END DO 
          ELSE 
            HPOP=ATTOT(1,ID) 
         END IF 
         IF(ID.NE.IDSTD) CYCLE DEPTHS 
         TSTD=TEMP(ID) 
         VTS=VTURB(ID) 
         DSTD=SQRT(1.4E7*TSTD+VTS) 
         WRITE(6,"(/' N/U  AT THE STANDARD DEPTH  (ID =',I3,              & 
     &         ' ; T,Ne = ',F8.1,1P2E12.3,' )'/                           & 
     &         ' --------------------------'//)")                         & 
     &         ID,TEMP(ID),ELEC(ID),hpop 
!        DO I=1,MATOM 
         DO I=1,30 
            WRITE(6,"(A4,1P8E9.2)") TYPAT(I),(RRR(ID,J,I),J=1,MION0-1) 
         END DO 
      END DO DEPTHS 
! 
!     2. "high-temperature" ionization fractions 
!         (using the Opacity Project ionization fractions) 
! 
      IF(TEFF.LT.0.) THEN 
         CALL FRAC1 
         ID=IDSTD 
         HPOP=DENS(ID)/WMM(ID)/YTOT(ID) 
         WRITE(6,"(/' N/U  AT THE STANDARD DEPTH  - OP DATA',             & 
     & '  (ID =',I3,' ; T,Ne = ',F8.1,1PE12.3,' )'//)")                   & 
     &    ID,TEMP(ID),ELEC(ID) 
         DO I=1,MATOM 
            WRITE(6,"(A4,(1P8E9.2))")                                     & 
     &      TYPAT(I),(RRR(ID,J,I)/hpop,J=1,MION) 
            ioniz(i)=i+1 
         END DO 
      END IF 
! 
!     3. "high-temperature" ionization fractions 
!         using new Franck Delahaye partition functions 
! 
      if(teff.gt.-2.) return 
      write(6,*) 'teff inimode', id,teff,dens(1) 
      write(6,*) 'iteos=',iteos 
 
      call frac_fd(iteos) 
      ID=IDSTD 
      HPOP=DENS(ID)/WMM(ID)/YTOT(ID) 
      write(6,*) 'inimod',idstd,dens(id),hpop 
      DO I=1,28 
         WRITE(6,"(A4,(1P8E9.2))") TYPAT(I),(RRR(ID,J,I)/hpop,J=1,i+1) 
      end do 
      t=temp(id) 
      ane=elec(id) 
      X=SQRT(T/ANE) 
      XMX=2.145E4*SQRT(X) 
      DO I=1,28 
         do j=1,i+1 
           X=REAL(J) 
           XMAX=XMX*SQRT(X) 
           CALL PARTF(I,J,T,ANE,XMAX,U) 
           PFSTD(J,I)=U 
         end do 
         write(6,"(A4,(1P8E9.2))") TYPAT(I),(PFSTD(J,I),J=1,i+1) 
      END DO 
 
      RETURN 
      END SUBROUTINE INIMOD 
! 
! 
! ******************************************************************** 
! 
      SUBROUTINE STATE(ID,TE,ANE,Q) 
! 
!     modified LTE Saha equations - possibly using 
!     radiation temperatures after 
!     Schaerer and Schmutz AA 288, 321, 1994 
! 
      use accura 
      use params 
      use eospar 
      use wincom 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: FFI(MION0) 
! 
      Q=0. 
      ATOMS: DO I=1,NATOMS 
         IF(LGR(I)) CYCLE ATOMS 
         ION=IONIZ(I) 
         RQ=0. 
         RS=1. 
!        T=TRAD(INPOT(I,1),ID) 
!        if(t.le.0.) t=te 
         T=TE
         X=SQRT(T/ANE) 
         XMX=2.145E4*SQRT(X) 
         CALL PARTF(I,1,T,ANE,XMX,UM) 
         PFSTD(1,I)=UM 
         JMAX=1 
         DO J=2,ION 
            J1=J-1 
!           T=TRAD(INPOT(I,J),ID) 
!           if(t.le.0.) t=te 
            T=TE
            TLN=LOG(T)*1.5 
            TK=BOLK*T 
            THL=11605./T 
            X=SQRT(T/ANE) 
            XMX=2.145E4*SQRT(X) 
            DCH=EH/XMX/XMX/TK 
            DCHT=DCH*J1 
            FI=36.113+TLN-THL*ENEV(I,J1)+DCHT 
            X=J 
            XMAX=XMX*SQRT(X) 
            CALL PARTF(I,J,T,ANE,XMAX,U) 
            PFSTD(J,I)=U 
            FI=EXP(FI)*U/UM/ANE 
            FFI(J)=FI 
            IF(FFI(J).GT.1.) JMAX=J 
            UM=U 
         END DO 
         RQ=JMAX-1 
         IF(JMAX.LT.ION) THEN 
            R=1. 
            RQ=JMAX-1 
            DO J=JMAX+1,ION 
               R=R*FFI(J) 
               RR(I,J)=R/PFSTD(J,I) 
               RS=RS+R 
               RQ=RQ+(J-1)*R 
            END DO 
         END IF 
         IF(JMAX.GT.1) THEN 
            R=1. 
            DO JJ=1,JMAX-1 
               J=JMAX-JJ 
               R=R/FFI(J+1) 
               RR(I,J)=R/PFSTD(J,I) 
               RS=RS+R 
               RQ=RQ+(J-1)*R 
            END DO 
         END IF 
         ABND(I)=ABNDD(I,ID) 
         RR(I,JMAX)=ABND(I)/RS 
         DO J=1,ION 
            IF(J.NE.JMAX) RR(I,J)=RR(I,J)*RR(I,JMAX) 
            if(rr(i,j).lt.1.e-35) rr(i,j)=0. 
         END DO 
         RR(I,JMAX)=RR(I,JMAX)/PFSTD(JMAX,I) 
         X=RQ/RS 
!        IF(LRM(I)) CYCLE ATOMS 
         if(i.gt.1) Q=X*ABND(I)+Q 
         anato(i,id)=rr(i,1)*pfstd(1,i) 
         pfato(i,id)=pfstd(1,i) 
         anion(i,id)=rr(i,2)*pfstd(2,i) 
         pfion(i,id)=pfstd(2,i) 
         ahn(id)=anato(1,id) 
         ahp(id)=anion(1,id) 
         ahen(id)=anato(2,id) 
      END DO ATOMS 
! 
      do i=2,30 
         anion2(i,id)=rr(i,3)*pfstd(3,i) 
      end do 
! 
      do imol=1,600 
         anmol(imol,id)=0. 
         pfmol(imol,id)=0. 
      end do 
! 
      RETURN 
      END SUBROUTINE STATE 
! 
! ******************************************************************** 
! 
      SUBROUTINE TINT 
! 
!     LOGARITHMIC INTERPOLATION COEFFICIENTS FOR INTERPOLATION OF 
!     TEMP(ID) TO THE VALUES  5000,10000,20000,40000 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: TT(4) 
      TT =  (/3.699, 4.000, 4.301, 4.602/) 
! 
      DO ID=1,ND 
         T=LOG10(TEMP(ID)) 
         J=3 
         IF(T.GT.TT(3)) J=4 
         JT(ID)=J 
         X=(TT(J)-TT(J-1))*(TT(J)-TT(J-2))*(TT(J-1)-TT(J-2)) 
         TI0(ID)=(T-TT(J-2))*(T-TT(J-1))*(TT(J-1)-TT(J-2))/X 
         TI1(ID)=(T-TT(J-2))*(TT(J)-T)*(TT(J)-TT(J-2))/X 
         TI2(ID)=(T-TT(J-1))*(T-TT(J))*(TT(J)-TT(J-1))/X 
      END DO 
      RETURN 
      END SUBROUTINE TINT 
! 
! ******************************************************************** 
! 
      SUBROUTINE INIBL0 
! 
!     AUXILIARY INITIALIZATION PROCEDURE 
! 
      use accura 
      use params 
      use modelp 
      use lindat 
      use synthp 
      use wincom 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      real(dp), parameter  :: un=1. 
      REAL(DP)             ::  ABSO(MFREQ),EMIS(MFREQ),SCAT(MFREQ) 
! 
      IF(IDSTD.EQ.0) THEN 
         ID1=5 
         NDSTEP=(ND-2*ID1)/2 
         IDSTD=2*ND/3 
       ELSE IF(IDSTD.LT.0) THEN 
         ID1=1 
         NDSTEP=-IDSTD 
         IDSTD=2*ND/3 
      END IF 
      if(imode.le.-3) ndstep=1 
! 
      alam0s=alam0 
      alasts=alast 
      cutof0s=cutof0 
      cutofss=cutofs 
      relops=relop 
      spaces=space 
! 
!     if ALAST.lt.0 - set up vacuum wavelengths everywhere 
! 
      vaclim=2000. 
      if(alast.lt.0.) then 
         alast=abs(alast) 
         alasts=alast 
         vaclim=1.e18 
      end if 
! 
      if(inlte.lt.10) then 
         lasdel=.true. 
       else if(inlte.le.20) then 
         inlte=inlte-10 
         lasdel=.false. 
       else if(inlte.le.30) then 
         inlte=inlte-20 
         ifreq=11 
         lasdel=.true. 
       else if(inlte.le.40) then 
         inlte=inlte-30 
         ifreq=11 
         lasdel=.false. 
      end if 
! 
      IF(IMODE.NE.2.AND.IMODE.NE.-4) CALL LISTS 
! 
!     VTB    - turbulent velocity (in km/s). In non-negative, this 
!              value overwrites the value given by the standard input 
! 
      read(55,*,iostat=ios) VTB 
      IF(IOS.EQ.0) THEN 
         if(ifwin.le.0.and.vtb.ge.0) then 
            WRITE(6,"(//' TURBULENT VELOCITY  -  CHANGED TO   VTURB =',   & 
     &      1PE10.3,'  KM/S'/' ------------------'/)") VTB 
            do id=1,nd 
               vturb(id)=vtb*vtb*1.e10 
            end do 
         end if 
! 
         TSTD=TEMP(IDSTD) 
         VTS=VTURB(IDSTD) 
         DSTD=SQRT(1.4E7*TSTD+VTS) 
      END IF 
! 
!     angle points (in case the specific intensities are evaluated 
! 
!     NMU0      -    number of angles: 
!               >0 - and if also ANG0>0, angles (mu's) equidistant 
!                    between 1 and ANG0 
!               >0 - and if also ANG0<0, angles (mu's) equidistant 
!                    between 0.7 and ANG0, and sinuses equidistatnt for 
!                    others 
!               <0 - angles read in the next record 
!     ANG0      -    minimum mu (see above) 
!     IFLUX     -    mode for evaluating angle-dependent intensities and 
!                    the corresponding flux: 
!               =0 - no specifiec intensities are evaluated; only usual 
!                    flux is stored (unit 7 and 17) 
!               =1 - specific intensities are evaluated; 
!                    and stored on unit 18 
!               =2 - (interesting only for the case of macroscopic 
!                    velocity field); specific intensities evaluated by 
!                    a simple formal solution (RESOLV) 
! 
      NMU0=1 
      ANG0=1. 
      ANGL(1)=1. 
      WANGL(1)=0. 
      IFLUX=0 

      win: if(ifwin.le.0) then 
         READ(55,*,iostat=ios) NMU0,ANG0,IFLUX 
         IF(IOS.EQ.0) THEN 
! 
!        determinantion of the angle points and weights 
! 
            IF(NMU0.LT.0) THEN 
               NMU0=IABS(NMU0) 
               READ(55,*) (ANGL(IMU),IMU=1,NMU0) 
               DO IMU=2,NMU0-1 
                  WANGL(IMU)=0.5*(ANGL(IMU-1)+ANGL(IMU+1)) 
               END DO 
               WANGL(1)=0.5*(ANGL(1)-ANGL(2)) 
               WANGL(NMU0)=0.5*(ANGL(NMU0-1)-ANGL(NMU0)) 
             ELSE 
               IF(ANG0.GT.0.) THEN 
                  IF(NMU0.GT.1) THEN 
                  DMU=(1.-ANG0)/(NMU0-1) 
                  DO IMU=1,NMU0 
                     ANGL(IMU)=1.-(IMU-1)*DMU 
                     WANGL(IMU)=DMU 
                  END DO 
                  WANGL(1)=0.5*DMU 
                  WANGL(NMU0-1)=0.5*DMU 
                  WANGL(NMU0)=2.*DMU 
                  END IF 
                ELSE 
                  ANGH=0.70710678 
                  DMU=ANGH/(NMU0-1) 
                  DO IMU=1,NMU0 
                     ANGL(IMU)=(IMU-1)*DMU 
                     ANGL(IMU)=SQRT(1.-ANGL(IMU)**2) 
                     IF(IMU.GT.1.AND.IMU.LT.NMU0)                         & 
     &                  WANGL(IMU)=0.5*(ANGL(IMU-1)+ANGL(IMU+1)) 
                  END DO 
                  WANGL(1)=0.5*(ANGL(1)-ANGL(2)) 
                  WANGL(NMU0)=0.5*(ANGL(NMU0-1)-ANGL(NMU0)) 
                  IF(ANG0.LT.0.) DMU=(ANGH+ANG0)/(NMU0-1) 
                  DO IMU=1,NMU0-2 
                     ANGL(IMU+NMU0)=ANGH-IMU*DMU 
                     WANGL(IMU+NMU0)=DMU 
                  END DO 
                  WANGL(NMU0)=WANGL(NMU0)+0.5*DMU 
                  WANGL(2*NMU0-3)=0.5*DMU 
                  WANGL(2*NMU0-2)=2.*DMU 
                  NMU0=2*NMU0-2 
               END IF 
            END IF 
            IF(NMU0.GT.0)                                                 & 
     &      WRITE(6,"(//' SPECIFIC INTENSITIES COMPUTED FOR',I3,          & 
     &            ' ANGLES  mu=cos(theta) ='/                             & 
     &            ' ---------------------------------',                   & 
     &            '------------------------'//                            & 
     &            (10F7.2))") NMU0,(ANGL(I),I=1,NMU0) 
         END IF 
       else 
         call velset 
      end if win 
! 
      IF(IMODE.EQ.-1) THEN 
         INLTE=0 
         CUTOF0=0. 
      END IF 
! 
!     continuum frequencies 
! 
      if(ifwin.le.0) then 
         alam0=alam0s 
         if(alam0s.eq.0.) alam0=5.e7/temp(1)/10. 
         if(alam0s.lt.0.) alam0=-5.e7/temp(1)/alam0s 
         alast=alasts 
         if(alasts.eq.0.) alast=5.e7/temp(1)*20. 
         if(alasts.lt.0.) alast=-5.e7/temp(1)*alasts 
!        if(alast.gt.1.e5) alast=1.e5 
         ALAMC=(ALAM0+ALAST)*0.5 
         if(space.eq.0.) space=4.3e-8*sqrt(temp(idstd))*alamc 
         if(space.lt.0.) space=-5.72e-8*sqrt(temp(idstd))*alamc*space 
         SPACF=2.997925E18/ALAMC/ALAMC*SPACE 
         WRITE(6,"(//'----------------------------------------------'/    & 
     &           ' BASIC INPUT PARAMETERS FOR SYNTHETIC SPECTRA'/         & 
     &           ' ---------------------------------------------'/        & 
     &           ' INITIAL LAMBDA',28X,1H=,F10.3,' ANGSTROMS'/            & 
     &           ' FINAL   LAMBDA',28X,1H=,F10.3,' ANGSTROMS'/            & 
     &           ' CUTOFF PARAMETER',26X,1H=,F10.3,' ANGSTROMS'/          & 
     &           ' MINIMUM VALUE OF (LINE OPAC.)/(CONT.OPAC) =',1PE10.1/  & 
     &           ' MAXIMUM FREQUENCY SPACING',17X,1H=,1PE10.3,'  I.E.',   & 
     &             0PF6.3,'  ANGSTROMS'/                                  & 
     &           ' ---------------------------------------------'/)")     & 
     &           ALAM0,ALAST,CUTOF0,RELOP,SPACF,SPACE 
         CUTOF0=0.1*CUTOF0 
         SPACE0=SPACE*0.1 
         ALAM0=1.e-1*ALAM0 
         ALAST=1.e-1*ALAST 
         ALAMC=ALAMC*0.1 
         ALST00=ALAST 
         FRLAST=2.997925e17/ALAST 
         NFREQ=2 
         FREQ(1)=2.997925e17/ALAM0 
         FREQ(2)=FRLAST 
! 
      else 
! 
         spacon=cutofs 
         IF(SPACON.EQ.0) SPACON=3. 
         XFR=(ALAST-ALAM0)/SPACON 
         NFREQC=int(XFR)+1 
         NFREQC=MIN(NFREQC,MFREQC) 
         NFREQC=MAX(NFREQC,2) 
         DLAMLO=LOG10(ALAST/ALAM0)/(NFREQC-1) 
         AL0L=LOG10(ALAM0) 
         alambe=alam0 
         DO IJ=1,NFREQC 
            AL=AL0L+(IJ-1)*DLAMLO 
            ALAM=EXP(2.3025851*AL) 
            WLAMC(IJ)=ALAM 
            FREQC(IJ)=2.997925E18/ALAM 
         END DO 
         ALAMC=(ALAM0+ALAST)*0.5 
         SPACF=2.997925E18/ALAMC/ALAMC*SPACE 
         WRITE(6,"(//'----------------------------------------------'/    & 
     &           ' BASIC INPUT PARAMETERS FOR SYNTHETIC SPECTRA'/         & 
     &           ' ---------------------------------------------'/        & 
     &           ' INITIAL LAMBDA',28X,1H=,F10.3,' ANGSTROMS'/            & 
     &           ' FINAL   LAMBDA',28X,1H=,F10.3,' ANGSTROMS'/            & 
     &           ' CUTOFF PARAMETER',26X,1H=,F10.3,' ANGSTROMS'/          & 
     &           ' MINIMUM VALUE OF (LINE OPAC.)/(CONT.OPAC) =',1PE10.1/  & 
     &           ' MAXIMUM FREQUENCY SPACING',17X,1H=,1PE10.3,'  I.E.',   & 
     &             0PF6.3,'  ANGSTROMS'/                                  & 
     &           ' ---------------------------------------------'/)")     & 
     &           ALAM0,ALAST,CUTOF0,RELOP,SPACF,SPACE 
         CUTOF0=0.1*CUTOF0 
         SPACE0=SPACE*0.1 
         ALAM0=1.e-1*ALAM0 
         ALAST=1.e-1*ALAST 
         ALAMC=ALAMC*0.1 
         ALST00=ALAST 
         FRLAST=2.997925e17/ALAST 
         NFREQ=2 
         FREQ(1)=2.997925e17/ALAM0 
         FREQ(2)=FRLAST 
! 
      end if 
! 
      CALL SIGAVS 
      ILIN0=0 
      IF(IHYDPR.NE.0) THEN 
         CALL HYDINI 
         CALL XENINI 
      END IF 
      IF(IHE1PR.GT.0.OR.IHE2PR.GT.0) CALL ALLOC_HEPRF
      IF(IHE1PR.GT.0) CALL HE1INI 
      IF(IHE2PR.GT.0) CALL HE2INI 
! 
!     auxiliary quantities for dissolved fractions 
! 
      DO ID=1,ND 
         CALL DWNFR0(ID) 
         CALL WNSTOR(ID) 
      END DO 
! 
!     pretabulate expansion coefficients for the Voigt function 
! 
      CALL PRETAB 
! 
!     calculate the characteristic standard opacity 
! 
      IF(IMODE.LE.2) THEN 
         if(ifwin.le.0.and.ndstep.eq.0) then 
! 
!        old procedure 
! 
            CALL CROSET 
            DO ID=1,ND 
               CALL OPAC(ID,ABSO,EMIS,SCAT) 
               ABSTD(ID)=MIN(ABSO(1),ABSO(2)) 
            END DO 
          else 
! 
!       new procedure 
! 
            if(ifwin.le.0) then 
               nfreqc=ifix(real(cutofs,4)) 
               if(nfreqc.eq.0) nfreqc=mfreq 
               all0=log(alam0) 
               all1=log(alast) 
               dlc=(all1-all0)/(nfreqc-1) 
               do ijc=1,nfreqc 
                  wlamc(ijc)=exp(all0+(ijc-1)*dlc) 
                  freqc(ijc)=2.997925e17/wlamc(ijc) 
               end do 
               CALL CROSEW 
               do id=1,nd 
                  CALL OPACON(ID) 
                  do ijc=1,nfreqc 
                     abstdw(ijc,id)=absoc(ijc) 
                  end do 
               end do 
             else 
               CALL CROSEW 
               DO ID=1,ND 
                  CALL OPACW(ID,ABSO,EMIS,0) 
                  DO IJ=1,NFREQC 
                     ABSTDW(IJ,ID)=ABSOC(IJ)/DENSCON(ID) 
                  END DO 
               END DO 
            end if 
         end if 
      END IF 
! 
      write(*,"(/'IDSTD, NDSTEP = ',2i5/)") IDSTD,NDSTEP 
      RETURN 
      END SUBROUTINE INIBL0 
! 
! *********************************************************************** 
! 
 
      SUBROUTINE INIBL1(IGRD) 
!     ======================= 
! 
!     AUXILIARY INITIALIZATION PROCEDURE 
! 
      use accura 
      use params 
      use modelp 
      use lindat 
      use synthp 
      use wincom 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      real(dp),parameter :: un=1.,bnc=1.4743e-2,hkc=4.79928e4,            & 
     &                      clc=2.997925e17 
      REAL(DP) ::           ABSO(MFREQ),EMIS(MFREQ),SCAT(MFREQ) 
! 
!     auxiliary quantities for dissolved fractions 
! 
      DO ID=1,ND 
         CALL DWNFR0(ID) 
         CALL WNSTOR(ID) 
         anh2(id)=0. 
         anhm(id)=0. 
         anch(id)=0. 
         anoh(id)=0. 
      END DO 
      CALL TINT 
! 
!     reset wavelengths in case of opacity grid calculations 
! 
      if(igrd.ge.0) then 
         alam0=alam0s 
         if(alam0s.eq.0.) alam0=5.e7/temp(1)/10. 
         if(alam0s.lt.0.) alam0=-5.e7/temp(1)/alam0s 
         alast=alasts 
         if(alasts.eq.0.) alast=5.e7/temp(1)*20. 
         if(alasts.lt.0.) alast=-5.e7/temp(1)*alasts 
         cutof0=cutof0s 
         cutofs=cutofss 
         relop=relops 
         if(relops.eq.0) then 
            relop=1.e-15 
            if(temp(1).lt.2.e6) relop=1.e-6 
            if(temp(1).lt.1.e6) relop=1.e-5 
            if(temp(1).lt.1.e5) relop=1.e-4 
         end if 
         space=spaces 
         ALAMC=(ALAM0+ALAST)*0.5 
         if(space.eq.0.) space=4.3e-8*sqrt(temp(idstd))*alamc 
         if(space.lt.0.) space=-5.72e-8*sqrt(temp(idstd))*alamc*space 
         SPACF=2.997925E18/ALAMC/ALAMC*SPACE 
         CUTOF0=0.1*CUTOF0 
         SPACE0=SPACE*0.1 
         ALAM0=1.D-1*ALAM0 
         ALAST=1.D-1*ALAST 
         ALAMC=ALAMC*0.1 
         ALST00=ALAST 
         FRLAST=CLC/ALAST 
! 
         nfreqc=ifix(real(cutofs,4)) 
         if(nfreqc.eq.0) nfreqc=mfreq 
         all0=log(alam0) 
         all1=log(alast) 
         dlc=(all1-all0)/(nfreqc-1) 
         xcc0=hkc/temp(1) 
         do ijc=1,nfreqc 
            wlamc(ijc)=exp(all0+(ijc-1)*dlc) 
            freqc(ijc)=clc/wlamc(ijc) 
         end do 
         id=1 
         CALL CROSEW 
         CALL OPACON(ID) 
         wc0=(freqc(1)-freqc(2))*0.5 
         wc1=(freqc(nfreqc-1)-freqc(nfreqc))*0.5 
         do ijc=2,nfreqc-1 
            absoc(ijc)=min(absoc(ijc),1.e30) 
            write(26,"(f11.3,1p5e13.5)")                                  & 
     &         wlamc(ijc)*10.,log(absoc(ijc)/dens(1)) 
         end do 
! 
         do ijc=1,nfreqc 
            abstdw(ijc,id)=absoc(ijc) 
         end do 
! 
      end if 
! 
!     calculate the characteristic standard opacity 
! 
      IF(IMODE.LE.2.and.imode.ge.-2) THEN 
         if(ifwin.le.0) then 
            CALL CROSET 
            DO ID=1,ND 
               CALL OPAC(ID,ABSO,EMIS,SCAT) 
               ABSTD(ID)=MIN(ABSO(1)+SCAT(1),ABSO(2)+SCAT(2)) 
            END DO 
          else 
            CALL CROSEW 
            DO ID=1,ND 
               CALL OPACW(ID,ABSO,EMIS,0) 
               DO IJ=1,NFREQC 
                  denscon(id)=1. 
                  ABSTDW(IJ,ID)=ABSOC(IJ)/DENSCON(ID) 
               END DO 
            END DO 
         end if 
      END IF 
! 
      RETURN 
      END SUBROUTINE INIBL1 
 
! 
! *********************************************************************** 
! 
      SUBROUTINE RESOLV 
! 
!     driver for evaluating opacities and emissivities which then 
!     enter the solution of the radiative transfer equation 
!     (RTE or RTEDFE) 
! 
      use accura 
      use params 
      use modelp 
      use lindat 
      use synthp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: ABSO(MFREQ),EMIS(MFREQ),SCAT(MFREQ) 
! 
      IHYL=-1 
! 
!     if(imode.le.-3) call abnchn(1) 
! 
!     set up the partial line list for the current interval 
! 
      CALL INISET 
!!    write(*,*) 'after INISET - iblank,nblank',iblank,nblank 
      if(ifmol.gt.0) then 
         do ilist=1,nmlist 
!!          write(*,*) 'before molset ilist',ilist,amlist(ilist) 
            call molset(ilist) 
!!          write(*,*) 'after molset ilist',ilist 
!!          write(*,*) 
         end do 
      end if 
! 
!     select possible hydrogen lines that may contribute to the opacity 
! 
      IF(IMODE.NE.-1) CALL HYLSET 
! 
!     select possible He II lines that may contribute to the opacity 
! 
      IF(IMODE.NE.-1) CALL HE2SET 
! 
!     output of information about selected lines 
! 
      CALL INIBLA 
      if(ifmol.gt.0) call iniblm 
! 
!     photoinization cross-sections 
! 
      CALL CROSET 
! 
!     monochromatic opacity and emissivity including all contributing 
!     lines and continua 
! 
      IF(IMODE.GE.-1) THEN 
         DO ID=1,ND 
!!          write(*,*) 'before opac',id 
            CALL OPAC(ID,ABSO,EMIS,SCAT) 
!!          write(*,*) 'after opac',id,abso(1),abso(2) 
            ABSTD(ID)=0.5*(ABSO(1)+ABSO(2)) 
            DO IJ=1,NFREQ 
               CH(IJ,ID)=ABSO(IJ) 
               ET(IJ,ID)=EMIS(IJ) 
               SC(IJ,ID)=SCAT(IJ) 
            END DO 
            if(imode0.eq.-4) call ougrid(abso) 
         END DO 
!!       write(*,*) 'after all opac' 
! 
!        output of information about selected hydrogen lines 
! 
         CALL INIBLH 
! 
!        the iron curtain or opacity table  option - output of monochromatic opa
! 
       ELSE IF(IMODE.EQ.-2) THEN 
         ID=1 
         write(27,"(1p3e15.4)") temp(id),dens(id),elec(id) 
         CALL OPAC(ID,ABSO,EMIS,SCAT) 
         DO IJ=3,NFREQ-1 
            ABSO(IJ)=(ABSO(IJ)+SCAT(IJ))/HPOP 
            WRITE(27,"(f15.3,1p2e15.5)") WLAM(IJ),ABSO(IJ),scat(ij) 
         END DO 
       else 
         id=1 
         call opac(id,abso,emis,scat) 
         ch(1,id)=abso(1) 
         ch(2,id)=abso(2) 
         call ougrid(abso) 
      END IF 
      RETURN 
      END SUBROUTINE RESOLV 
! 
! ******************************************************************* 
! 
      SUBROUTINE RTE 
! 
!     solution of the radiative transfer equation by Feautrier method 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use lindat 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: D(3,3,MDEPTH),ANU(3,MDEPTH),AANU(MDEPTH),DDD(MDEPTH),   & 
     &            AA(3,3),BB(3,3),CC(3,3),VL(3),AMU(3),WTMU(3),           & 
     &            DT(MDEPTH),TAU(MDEPTH),                                 & 
     &            RDD(MDEPTH),FKK(MDEPTH),ST0(MDEPTH),SS0(MDEPTH),        & 
     &            RINT(MDEPTH,MMU) 
      CHARACTER(LEN=4) :: TYPION(9) 
      REAL(DP) :: CINT1(MDEPTH),CINT2(MDEPTH),                            & 
     &            CTRI(MDEPTH),CTRR(MDEPTH),XKAR(MDEPTH),                 & 
     &            ABXLI(MFREQ),EMXLI(MFREQ) 
      REAL(DP), PARAMETER :: UN=1.e0, HALF=0.5e0 
      REAL(DP), PARAMETER :: THIRD=UN/3., QUART=UN/4., SIXTH=UN/6.,       & 
     &                       TAUREF = 2./3. 
      DATA AMU/.887298334620742e0,.5e0,.112701665379258e0/,               & 
     &     WTMU/.277777777777778e0,.444444444444444e0,.277777777777778e0  & 
     &         / 
      DATA TYPION /' I  ',' II ',' III',' IV ',' V  ',                    & 
     &             ' VI ',' VII','VIII',' IX '/ 
! 
      NMU=3 
      ND1=ND-1 
! 
!     Overall loop over frequencies 
! 
      FRLOOP: DO IJ=1,NFREQ 
      TAUMIN=CH(IJ,1)/DENS(1)*DM(1)*HALF 
      TAU(1)=TAUMIN 
      IREF=1 
      DO I=1,ND1 
         DT(I)=(DM(I+1)-DM(I))*(CH(IJ,I+1)/DENS(I+1)+CH(IJ,I)/DENS(I))*   & 
     &         HALF 
         ST0(I)=ET(IJ,I)/CH(IJ,I) 
         SS0(I)=-SC(IJ,I)/CH(IJ,I) 
         TAU(I+1)=TAU(I)+DT(I) 
         IF(TAU(I).LE.TAUREF.AND.TAU(I+1).GT.TAUREF) IREF=I 
      END DO 
      IREFD(IJ)=IREF 
      ST0(ND)=ET(IJ,ND)/CH(IJ,ND) 
      SS0(ND)=-SC(IJ,ND)/CH(IJ,ND) 
      FR=FREQ(IJ) 
      BNU=BN*(FR*1.E-15)**3 
      PLAND=BNU/(EXP(HK*FR/TEMP(ND  ))-UN) 
      DPLAN=BNU/(EXP(HK*FR/TEMP(ND-1))-UN) 
      DPLAN=(PLAND-DPLAN)/DT(ND1) 
! 
!   +++++++++++++++++++++++++++++++++++++++++ 
!   FIRST PART  -  VARIABLE EDDINGTON FACTORS 
!   +++++++++++++++++++++++++++++++++++++++++ 
! 
      ALB1=0. 
      DO I=1,NMU 
! 
!   ************************ 
!   UPPER BOUNDARY CONDITION 
!   ************************ 
! 
         ID=1 
         DTP1=DT(1) 
         Q0=0. 
         P0=0. 
! 
!       allowance for non-zero optical depth at the first depth point 
! 
         TAMM=TAUMIN/AMU(I) 
         IF(TAMM.GT.0.01) THEN 
            P0=UN-EXP(-TAMM) 
          ELSE 
            P0=TAMM*(UN-HALF*TAMM*(UN-TAMM*THIRD*(UN-QUART*TAMM))) 
         END IF 
         EX=UN-P0 
         Q0=Q0+P0*AMU(I)*WTMU(I) 
! 
         DIV=DTP1/AMU(I)*THIRD 
         VL(I)=DIV*(ST0(ID)+HALF*ST0(ID+1))+ST0(ID)*P0 
         DO J=1,NMU 
            BB(I,J)=SS0(ID)*WTMU(J)*(DIV+P0)-ALB1*WTMU(J) 
            CC(I,J)=-HALF*DIV*SS0(ID+1)*WTMU(J) 
         END DO 
         BB(I,I)=BB(I,I)+AMU(I)/DTP1+UN+DIV 
         CC(I,I)=CC(I,I)+AMU(I)/DTP1-HALF*DIV 
         ANU(I,ID)=0. 
      END DO 
! 
!     Matrix inversion: instead of calling MATINV, a very fast inlined 
!     routine MINV3 for a specific 3 x 3 matrix inversion 
! 
!     CALL MATINV(BB,NMU,3) 
! 
!     ****************************** 
      BB(2,1)=BB(2,1)/BB(1,1) 
      BB(2,2)=BB(2,2)-BB(2,1)*BB(1,2) 
      BB(2,3)=BB(2,3)-BB(2,1)*BB(1,3) 
      BB(3,1)=BB(3,1)/BB(1,1) 
      BB(3,2)=(BB(3,2)-BB(3,1)*BB(1,2))/BB(2,2) 
      BB(3,3)=BB(3,3)-BB(3,1)*BB(1,3)-BB(3,2)*BB(2,3) 
! 
      BB(3,2)=-BB(3,2) 
      BB(3,1)=-BB(3,1)-BB(3,2)*BB(2,1) 
      BB(2,1)=-BB(2,1) 
! 
      BB(3,3)=UN/BB(3,3) 
      BB(2,3)=-BB(2,3)*BB(3,3)/BB(2,2) 
      BB(2,2)=UN/BB(2,2) 
      BB(1,3)=-(BB(1,2)*BB(2,3)+BB(1,3)*BB(3,3))/BB(1,1) 
      BB(1,2)=-BB(1,2)*BB(2,2)/BB(1,1) 
      BB(1,1)=UN/BB(1,1) 
! 
      BB(1,1)=BB(1,1)+BB(1,2)*BB(2,1)+BB(1,3)*BB(3,1) 
      BB(1,2)=BB(1,2)+BB(1,3)*BB(3,2) 
      BB(2,1)=BB(2,2)*BB(2,1)+BB(2,3)*BB(3,1) 
      BB(2,2)=BB(2,2)+BB(2,3)*BB(3,2) 
      BB(3,1)=BB(3,3)*BB(3,1) 
      BB(3,2)=BB(3,3)*BB(3,2) 
!     ****************************** 
! 
      DO I=1,NMU 
         DO J=1,NMU 
            S=0. 
            DO K=1,NMU 
               S=S+BB(I,K)*CC(K,J) 
            END DO 
            D(I,J,ID)=S 
            ANU(I,1)=ANU(I,1)+BB(I,J)*VL(J) 
         END DO 
      END DO 
! 
!   ******************* 
!   NORMAL DEPTH POINTS 
!   ******************* 
! 
      DO ID=2,ND1 
         DTM1=DTP1 
         DTP1=DT(ID) 
         DT0=HALF*(DTM1+DTP1) 
         AL=UN/DTM1/DT0 
         GA=UN/DTP1/DT0 
         BE=AL+GA 
         A=(UN-HALF*AL*DTP1*DTP1)*SIXTH 
         C=(UN-HALF*GA*DTM1*DTM1)*SIXTH 
         B=UN-A-C 
         VL0=A*ST0(ID-1)+B*ST0(ID)+C*ST0(ID+1) 
         DO I=1,NMU 
            DO J=1,NMU 
               AA(I,J)=-A*SS0(ID-1)*WTMU(J) 
               CC(I,J)=-C*SS0(ID+1)*WTMU(J) 
               BB(I,J)=B*SS0(ID)*WTMU(J) 
            END DO 
         END DO 
         DO I=1,NMU 
            DIV=AMU(I)**2 
            VL(I)=VL0 
            AA(I,I)=AA(I,I)+DIV*AL-A 
            CC(I,I)=CC(I,I)+DIV*GA-C 
            BB(I,I)=BB(I,I)+DIV*BE+B 
         END DO 
         DO I=1,NMU 
            S1=0. 
            DO J=1,NMU 
               S=0. 
               S1=S1+AA(I,J)*ANU(J,ID-1) 
               DO K=1,NMU 
                  S=S+AA(I,K)*D(K,J,ID-1) 
               END DO 
               BB(I,J)=BB(I,J)-S 
            END DO 
            VL(I)=VL(I)+S1 
         END DO 
! 
!     Matrix inversion: instead of calling MATINV, a very fast inlined 
!     routine MINV3 for a specific 3 x 3 matrix inversion 
! 
!     CALL MATINV(BB,NMU,3) 
! 
!     ****************************** 
      BB(2,1)=BB(2,1)/BB(1,1) 
      BB(2,2)=BB(2,2)-BB(2,1)*BB(1,2) 
      BB(2,3)=BB(2,3)-BB(2,1)*BB(1,3) 
      BB(3,1)=BB(3,1)/BB(1,1) 
      BB(3,2)=(BB(3,2)-BB(3,1)*BB(1,2))/BB(2,2) 
      BB(3,3)=BB(3,3)-BB(3,1)*BB(1,3)-BB(3,2)*BB(2,3) 
! 
      BB(3,2)=-BB(3,2) 
      BB(3,1)=-BB(3,1)-BB(3,2)*BB(2,1) 
      BB(2,1)=-BB(2,1) 
! 
      BB(3,3)=UN/BB(3,3) 
      BB(2,3)=-BB(2,3)*BB(3,3)/BB(2,2) 
      BB(2,2)=UN/BB(2,2) 
      BB(1,3)=-(BB(1,2)*BB(2,3)+BB(1,3)*BB(3,3))/BB(1,1) 
      BB(1,2)=-BB(1,2)*BB(2,2)/BB(1,1) 
      BB(1,1)=UN/BB(1,1) 
! 
      BB(1,1)=BB(1,1)+BB(1,2)*BB(2,1)+BB(1,3)*BB(3,1) 
      BB(1,2)=BB(1,2)+BB(1,3)*BB(3,2) 
      BB(2,1)=BB(2,2)*BB(2,1)+BB(2,3)*BB(3,1) 
      BB(2,2)=BB(2,2)+BB(2,3)*BB(3,2) 
      BB(3,1)=BB(3,3)*BB(3,1) 
      BB(3,2)=BB(3,3)*BB(3,2) 
!     ****************************** 
! 
         DO I=1,NMU 
            ANU(I,ID)=0. 
            DO J=1,NMU 
               S=0. 
               DO K=1,NMU 
                  S=S+BB(I,K)*CC(K,J) 
               END DO 
               D(I,J,ID)=S 
               ANU(I,ID)=ANU(I,ID)+BB(I,J)*VL(J) 
            END DO 
         END DO 
      END DO 
! 
!   ************ 
!   LOWER BOUNDARY CONDITION 
!   ************ 
! 
      ID=ND 
! 
!     First option: 
!     b.c. is different from stellar atmospheres; expresses symmetry 
!     at the central plane   I(taumax,-mu,nu)=I(taumax,+mu,nu) 
! 
      IF(IFZ0.EQ.0) THEN 
         B=DTP1*HALF 
         A=0. 
         DO I=1,NMU 
            BI=B/AMU(I) 
            AI=A/AMU(I) 
            VL(I)=ST0(ID)*BI+ST0(ID-1)*AI 
            DO J=1,NMU 
               AA(I,J)=-AI*SS0(ID-1)*WTMU(J) 
               BB(I,J)=BI*SS0(ID)*WTMU(J) 
            END DO 
            AA(I,I)=AA(I,I)+AMU(I)/DTP1-AI 
            BB(I,I)=BB(I,I)+AMU(I)/DTP1+BI 
         END DO 
         DO I=1,NMU 
            S1=0. 
            DO J=1,NMU 
               S=0. 
               S1=S1+AA(I,J)*ANU(J,ID-1) 
               DO K=1,NMU 
                  S=S+AA(I,K)*D(K,J,ID-1) 
               END DO 
               BB(I,J)=BB(I,J)-S 
            END DO 
            VL(I)=VL(I)+S1 
         END DO 
! 
!     Second option: 
!     b.c. is the same as in stellar atmospheres - the last depth point 
!     is not at the central plane 
! 
      ELSE 
         DO I=1,NMU 
            AA(I,I)=AMU(I)/DTP1 
            VL(I)=PLAND+AMU(I)*DPLAN+AA(I,I)*ANU(I,ID-1) 
            DO J=1,NMU 
               BB(I,J)=-AA(I,I)*D(I,J,ID-1) 
            END DO 
            BB(I,I)=BB(I,I)+AA(I,I)+UN 
         END DO 
      END IF 
! 
!     Matrix inversion: instead of calling MATINV, a very fast inlined 
!     routine MINV3 for a specific 3 x 3 matrix inversion 
! 
!     CALL MATINV(BB,NMU,3) 
! 
!     ****************************** 
      BB(2,1)=BB(2,1)/BB(1,1) 
      BB(2,2)=BB(2,2)-BB(2,1)*BB(1,2) 
      BB(2,3)=BB(2,3)-BB(2,1)*BB(1,3) 
      BB(3,1)=BB(3,1)/BB(1,1) 
      BB(3,2)=(BB(3,2)-BB(3,1)*BB(1,2))/BB(2,2) 
      BB(3,3)=BB(3,3)-BB(3,1)*BB(1,3)-BB(3,2)*BB(2,3) 
! 
      BB(3,2)=-BB(3,2) 
      BB(3,1)=-BB(3,1)-BB(3,2)*BB(2,1) 
      BB(2,1)=-BB(2,1) 
! 
      BB(3,3)=UN/BB(3,3) 
      BB(2,3)=-BB(2,3)*BB(3,3)/BB(2,2) 
      BB(2,2)=UN/BB(2,2) 
      BB(1,3)=-(BB(1,2)*BB(2,3)+BB(1,3)*BB(3,3))/BB(1,1) 
      BB(1,2)=-BB(1,2)*BB(2,2)/BB(1,1) 
      BB(1,1)=UN/BB(1,1) 
! 
      BB(1,1)=BB(1,1)+BB(1,2)*BB(2,1)+BB(1,3)*BB(3,1) 
      BB(1,2)=BB(1,2)+BB(1,3)*BB(3,2) 
      BB(2,1)=BB(2,2)*BB(2,1)+BB(2,3)*BB(3,1) 
      BB(2,2)=BB(2,2)+BB(2,3)*BB(3,2) 
      BB(3,1)=BB(3,3)*BB(3,1) 
      BB(3,2)=BB(3,3)*BB(3,2) 
!     ****************************** 
! 
      DO I=1,NMU 
         ANU(I,ID)=0. 
         DO J=1,NMU 
            D(I,J,ID)=0. 
            ANU(I,ID)=ANU(I,ID)+BB(I,J)*VL(J) 
         END DO 
      END DO 
! 
!   ************ 
!   BACKSOLUTION 
!   ************ 
! 
      ID=ND 
      FKK(ND)=THIRD 
      AJ=0. 
      AK=0. 
      DO I=1,NMU 
         RMU=WTMU(I)*ANU(I,ID) 
         AJ=AJ+RMU 
         AK=AK+RMU*AMU(I)*AMU(I) 
      END DO 
      RDD(ID)=AJ 
      FKK(ND)=AK/AJ 
      DO ID=ND-1,1,-1 
         DO I=1,NMU 
            DO J=1,NMU 
               ANU(I,ID)=ANU(I,ID)+D(I,J,ID)*ANU(J,ID+1) 
            END DO 
         END DO 
         AJ=0. 
         AK=0. 
         DO I=1,NMU 
            DIV=WTMU(I)*ANU(I,ID) 
            AJ=AJ+DIV 
            AK=AK+DIV*AMU(I)**2 
         END DO 
         FKK(ID)=AK/AJ 
      END DO 
! 
!     surface Eddington actor 
! 
      AH=0. 
      DO I=1,NMU 
         AH=AH+WTMU(I)*AMU(I)*ANU(I,1) 
      END DO 
      FH=AH/AJ-HALF*ALB1 
! 
!     FKK(ND)=THIRD 
! 
! 
!   +++++++++++++++++++++++++++++++++++++++++ 
!   SECOND PART  -  DETERMINATION OF THE MEAN INTENSITIES 
!   RECALCULATION OF THE TRANSFER EQUATION WITH GIVEN EDDINGTON FACTORS 
!   +++++++++++++++++++++++++++++++++++++++++ 
! 
      DTP1=DT(1) 
      DIV=DTP1*THIRD 
      BBB=FKK(1)/DTP1+FH+DIV+SS0(1)*(DIV+Q0) 
      CCC=FKK(2)/DTP1-HALF*DIV*(UN+SS0(2)) 
      VLL=DIV*(ST0(1)+HALF*ST0(2))+ST0(1)*Q0 
      AANU(1)=VLL/BBB 
      DDD(1)=CCC/BBB 
      DO ID=2,ND1 
         DTM1=DTP1 
         DTP1=DT(ID) 
         DT0=HALF*(DTP1+DTM1) 
         AL=UN/DTM1/DT0 
         GA=UN/DTP1/DT0 
         A=(UN-HALF*DTP1*DTP1*AL)*SIXTH 
         C=(UN-HALF*DTM1*DTM1*GA)*SIXTH 
         AAA=AL*FKK(ID-1)-A*(UN+SS0(ID-1)) 
         CCC=GA*FKK(ID+1)-C*(UN+SS0(ID+1)) 
         BBB=(AL+GA)*FKK(ID)+(UN-A-C)*(UN+SS0(ID)) 
         VLL=A*ST0(ID-1)+C*ST0(ID+1)+(UN-A-C)*ST0(ID) 
         BBB=BBB-AAA*DDD(ID-1) 
         DDD(ID)=CCC/BBB 
         AANU(ID)=(VLL+AAA*AANU(ID-1))/BBB 
      END DO 
! 
!     Lower boundary condition 
!     1.option -  different from stellar atmospheres 
! 
      IF(IFZ0.EQ.0) THEN 
         B=DTP1*HALF 
         BBB=FKK(ND)/DTP1+B*(UN+SS0(ND)) 
         AAA=FKK(ND-1)/DTP1 
         VLL=B*ST0(ND) 
       ELSE 
! 
!     Lower boundary condition 
!     2.option - stellar atmospheric 
! 
         BBB=FKK(ND)/DTP1+HALF 
         AAA=FKK(ND1)/DTP1 
         VLL=HALF*PLAND+DPLAN*THIRD 
      END IF 
      BBB=BBB-AAA*DDD(ND1) 
      RDD(ND)=(VLL+AAA*AANU(ND1))/BBB 
      DO IID=1,ND1 
         ID=ND-IID 
         RDD(ID)=AANU(ID)+DDD(ID)*RDD(ID+1) 
      END DO 
      FLUX(IJ)=FH*RDD(1) 
! 
!     if needed (if iprin.ge.3), output of interesting physical 
!     quantities at the monochromatic optical depth  tau(nu)=2/3 
! 
      IF(IPRIN.ge.3) THEN 
      T0=LOG(TAU(IREF+1)/TAU(IREF)) 
      X0=LOG(TAU(IREF+1)/TAUREF)/T0 
      X1=LOG(TAUREF/TAU(IREF))/T0 
      DMREF=EXP(LOG(DM(IREF))*X0+LOG(DM(IREF+1))*X1) 
      TREF=EXP(LOG(TEMP(IREF))*X0+LOG(TEMP(IREF+1))*X1) 
      STREF=EXP(LOG(ST0(IREF))*X0+LOG(ST0(IREF+1))*X1) 
      SCREF=EXP(LOG(-SS0(IREF))*X0+LOG(-SS0(IREF+1))*X1) 
      SSREF=EXP(LOG(-SS0(IREF)*RDD(IREF))*X0+                             & 
     &           LOG(-SS0(IREF+1)*RDD(IREF+1))*X1) 
      SREF=STREF+SSREF 
      ALM=2.997925E18/FREQ(IJ) 
!     WRITE(96,"(I3,F10.3,I4,1PE10.3,0PF10.1,1X,1P3E10.3,E11.3)") 
!    *   IJ,ALM,IREF,DMREF,TREF,SCREF,STREF,SSREF,SREF 
      END IF 
! 
!   THIRD PART  -  DETERMINATION OF THE SPECIFIC INTENSITIES 
!   RECALCULATION OF THE TRANSFER EQUATION WITH GIVEN SOURCE FUNCTION 
! 
      if(iflux.eq.0) return 
      DO IMU=1,NMU0 
      ANX=ANGL(IMU) 
      DTP1=DT(1) 
      DIV=DTP1*THIRD/ANX 
! 
      TAMM=TAUMIN/ANX 
      IF(TAMM.LT.0.01) THEN 
         P0=TAMM*(UN-HALF*TAMM*(UN-TAMM*THIRD*(UN-QUART*TAMM))) 
       ELSE 
         P0=UN-EXP(-TAMM) 
      END IF 
! 
      BBB=ANX/DTP1+UN+DIV 
      CCC=ANX/DTP1-HALF*DIV 
      VLL=(DIV+P0)*(ST0(1)-SS0(1)*RDD(1))                                 & 
     &    +HALF*DIV*(ST0(2)-SS0(2)*RDD(2)) 
      AANU(1)=VLL/BBB 
      DDD(1)=CCC/BBB 
      DIV=ANX*ANX 
      DO ID=2,ND1 
         DTM1=DT(ID-1) 
         DTP1=DT(ID) 
         DT0=HALF*(DTP1+DTM1) 
         AL=UN/DTM1/DT0 
         GA=UN/DTP1/DT0 
         A=(UN-HALF*DTP1*DTP1*AL)*SIXTH 
         C=(UN-HALF*DTM1*DTM1*GA)*SIXTH 
         AAA=DIV*AL-A 
         CCC=DIV*GA-C 
         BBB=DIV*(AL+GA)+UN-A-C 
         VLL=A*(ST0(ID-1)-SS0(ID-1)*RDD(ID-1))+                           & 
     &       C*(ST0(ID+1)-SS0(ID+1)*RDD(ID+1))+                           & 
     &       (UN-A-C)*(ST0(ID)-SS0(ID)*RDD(ID)) 
         BBB=BBB-AAA*DDD(ID-1) 
         DDD(ID)=CCC/BBB 
         AANU(ID)=(VLL+AAA*AANU(ID-1))/BBB 
      END DO 
! 
!     Lower boundary condition 
!     1.option -  different from stellar atmospheres 
! 
      IF(IFZ0.EQ.0) THEN 
         B=DTP1*HALF/ANX 
         BBB=ANX/DTP1+B*(UN+SS0(ND)) 
         AAA=ANX/DTP1 
         VLL=B*ST0(ND) 
       ELSE 
! 
!     Lower boundary condition 
!     2.option - stellar atmospheric 
! 
         AAA=ANX/DTP1 
         BBB=AAA+UN 
         VLL=PLAND+ANX*DPLAN 
      END IF 
! 
      RINT(ND,IMU)=(VLL+AAA*AANU(ND1))/(BBB-AAA*DDD(ND1)) 
      DO IID=1,ND1 
         ID=ND-IID 
         RINT(ID,IMU)=AANU(ID)+DDD(ID)*RINT(ID+1,IMU) 
      END DO 
      END DO 
! 
      FLX=0. 
      DO IMU=1,NMU0 
         RINT(1,IMU)=RINT(1,IMU)/HALF 
         FLX=FLX+ANGL(IMU)*WANGL(IMU)*RINT(1,IMU) 
      END DO 
      FLX=FLX*HALF 
!     FLUX(IJ)=FLX 
! 
!     output of emergent specific intensities to Unit 10 
!     and 18 (continuum) 
! 
      IF(IJ.GT.2) THEN 
         WRITE(10,"(f10.3,1pe15.5/(1P5E15.5))")                           & 
     &      WLAM(IJ),FLX,(RINT(1,IMU),IMU=1,NMU0) 
       ELSE 
         WRITE(18,"(f10.3,1pe15.5/(1P5E15.5))")                           & 
     &      WLAM(IJ),FLX,(RINT(1,IMU),IMU=1,NMU0) 
      END IF 
! 
      if(iprin.ne.4) cycle frloop 
! 
!     compute contribution function C_i (ctri) and C_r (ctrr) 
!     following Magain (1986, A&A 163, 135) 
! 
      if(ijctr(ij).gt.0) then 
         xfr0=(freq(ij)-freq(2))/(freq(1)-freq(2)) 
         tauc=ch(1,1)/dens(1)*dm(1)*half 
         do id=1,nd 
           chc1=ch(1,id) 
           chc2=ch(2,id) 
           chcc=chc2+xfr0*(chc1-chc2) 
           etc1=et(1,id) 
           etc2=et(2,id) 
           etcc=etc2+xfr0*(etc1-etc2) 
           stcc=etcc/chcc 
           cint=cint2(id)+xfr0*(cint1(id)-cint2(id)) 
           avx=(chc1+chc2)*0.5*relop 
           call linop(id,abxli,emxli,avx) 
           sli0=emxli(ij)/abxli(ij) 
           abt0=ch(ij,id) 
           emt0=et(ij,id) 
           stt0=emt0/abt0 
           Xkar(id)=abxli(ij)+chcc*stcc/cint 
           ctri(id)=tauc*abt0/chc1*stt0*exp(-tau(id)) 
           if(tau(id).gt.70.) ctri(id)=0. 
           ctrr(id)=tauc/chc1*abxli(ij)*(un-sli0/cint) 
           if(id.lt.nd) then 
             dtc=(ch(1,id+1)/dens(id+1)+ch(1,id)/dens(id)) 
             tauc=tauc+half*dtc*(dm(id+1)-dm(id)) 
           endif 
         end do 
         taurs=Xkar(1)/dens(1)*dm(1)*half 
         xcti=ctri(1)*half*(dm(2)-dm(1)) 
         xctr=ctrr(1)*half*(dm(2)-dm(1)) 
         do i=1,nd-1 
           ctrr(i)=ctrr(i)*exp(-taurs) 
           if(i.eq.1) xctr=xctr*exp(-taurs) 
           if(i.gt.1) then 
             xcti=xcti+ctri(i)*half*(dm(i+1)-dm(i-1)) 
             xctr=xctr+ctrr(i)*half*(dm(i+1)-dm(i-1)) 
           endif 
           if(taurs.gt.70.) ctrr(i)=0. 
           dtrs=(dm(i+1)-dm(i))*(Xkar(i+1)/dens(i+1)+Xkar(i)/dens(i)) 
           taurs=taurs+half*dtrs 
         end do 
         ctrr(nd)=0. 
         alam=2.997925e18/freq(ij) 
         il0=ijctr(ij) 
         il=indlin(il0) 
         iat=indat(il)/100 
         ion=mod(indat(il),100) 
         write(97,"(i5,f11.4,2x,2a4,i8,1pe12.4,0pf10.1)")                 & 
     &      il,alam,typat(iat),typion(ion),iref,dmref,tref 
         do id=1,nd 
           ctrip=ctri(id)/xcti 
           ctrrp=ctrr(id)/xctr 
           write(97,"(i4,1p4e12.4)") id,dm(id),tau(id),ctrip,ctrrp 
         end do 
       else if(ij.eq.1) then 
         do id=1,nd 
           cint1(id)=rint(id,nmu0) 
         end do 
       else if(ij.eq.2) then 
         do id=1,nd 
           cint2(id)=rint(id,nmu0) 
         end do 
      endif 
! 
!     end of the global loop over frequencies 
! 
      END DO FRLOOP 
      RETURN 
      END SUBROUTINE RTE 
! 
! ******************************************************************** 
! 
      SUBROUTINE OUTPRI 
! 
!     Output of synthetic spectrum 
! 
!     Output onto unit 7 serves as an input to the next program 
!     ROTINS, which performs convolutions for the rotational and 
!     instrumental broadening, and plots the synthetic spectrum 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: UN=1.,CAS=1./2.997925e18,EQWC=1.19917e22 
      REAL(DP), PARAMETER :: PI2=3.141592654/2. 
      REAL(DP)            :: FLX(3),REL(3),ALX(3) 
! 
      if(ifwin.le.0) then 
! 
!     output of synthetic spectrum on unit 7 
! 
         DO IJ=3,NFREQ-1 
            FLAM=FLUX(IJ)*FREQ(IJ)*FREQ(IJ)*CAS 
            WRITE(7,"(F12.5,1PE15.5)") WLAM(IJ),FLAM 
         END DO 
! 
!        output of the continuum flux on unit 17 
! 
         FLAM=FLUX(1)*FREQ(1)*FREQ(1)*CAS 
         WRITE(17,"(F12.5,1PE15.5)") WLAM(1),FLAM 
         IF(IBLANK.EQ.NBLANK) THEN 
            FLAM=FLUX(NFREQ)*FREQ(NFREQ)*FREQ(NFREQ)*CAS 
            WRITE(7,"(F12.5,1PE15.5)") WLAM(NFREQ),FLAM 
            FLAM=FLUX(2)*FREQ(2)*FREQ(2)*CAS 
            WRITE(17,"(F12.5,1PE15.5)") WLAM(2),FLAM 
         END IF 
       else 
         DO IJ=1,NFROBS 
            FLAM=FLUX(IJ)*FRQOBS(IJ)*FRQOBS(IJ)*CAS*0.5 
            flam=max(flam,1.e-40) 
            WRITE(7,"(F12.5,1PE15.5)") WLOBS(IJ),FLAM 
         END DO 
      end if 
! 
!     unit 6 and 16 outputs 
! 
      if(iprin.lt.3) return 
      if(iprin.ge.3) then 
         WRITE(6,"(/' EMERGENT RADIATION'/' ------------------'/)") 
         WRITE(6,"(3('   LAMBDA  LOG HLAM    REL')/)") 
      end if 
      K1=0 
      EQW=0. 
      EQWP=0. 
      IF(IBLANK.EQ.1) EQWT=0. 
      IF(IBLANK.EQ.1) EQWTP=0. 
      XX=UN/(FREQ(2)-FREQ(1)) 
      XXX=UN/(FREQ(1)+FREQ(2))/(FREQ(1)+FREQ(2)) 
      if(ifwin.le.0) then 
         DO IJ=1,NFREQ 
            FLAM=FLUX(IJ)*FREQ(IJ)*FREQ(IJ)*CAS 
            CONT=((FREQ(IJ)-FREQ(1))*FLUX(2)+(FREQ(2)-FREQ(IJ))*          & 
     &           FLUX(1))*XX 
            RE0=FLUX(IJ)/CONT 
            EQW=EQW+(UN-RE0)*W(IJ) 
            REP=RE0 
            IF(REP.GT.UN) REP=UN 
            EQWP=EQWP+(UN-REP)*W(IJ) 
            K1=K1+1 
            FLX(K1)=LOG10(FLAM) 
            ALX(K1)=WLAM(IJ) 
            REL(K1)=RE0 
            IF(K1.EQ.3.OR.IJ.EQ.NFREQ) THEN 
               WRITE(6,"(3(2X,F9.3,F8.4,F7.3))")                          & 
     &         (ALX(I),FLX(I),REL(I),I=1,K1) 
               K1=0 
            END IF 
         END DO 
       else 
         DO IJ=1,NFROBS 
            FLAM=FLUX(IJ)*FREQ(IJ)*FREQ(IJ)*CAS 
            CONT=((FRQOBS(IJ)-FREQ(1))*FLUX(2)+                           & 
     &        (FREQ(2)-FRQOBS(IJ))*FLUX(1))*XX 
            RE0=FLUX(IJ)/CONT 
            EQW=EQW+(UN-RE0)*W(IJ) 
            REP=RE0 
            IF(REP.GT.UN) REP=UN 
            EQWP=EQWP+(UN-REP)*W(IJ) 
            if(iprin.gt.0) then 
               K1=K1+1 
               FLX(K1)=LOG10(FLAM) 
               ALX(K1)=WLAM(IJ) 
               REL(K1)=RE0 
               IF(K1.EQ.3.OR.IJ.EQ.NFREQ) THEN 
                  WRITE(6,"(3(2X,F9.3,F8.4,F7.3))")                       & 
     &            (ALX(I),FLX(I),REL(I),I=1,K1) 
                  K1=0 
               END IF 
            end if 
         END DO 
      end if 
! 
!     output of partial equivalent widths on unit 16 
! 
      EQW=EQW*EQWC*XXX 
      EQWT=EQWT+EQW 
      EQWP=EQWP*EQWC*XXX 
      EQWTP=EQWTP+EQWP 
      if(iprin.gt.2) WRITE(6,"(/,'  EQUIVALENT WIDTH THIS SET  =',2F8.1   & 
     &   ,' mA'/ '  EQUIVALENT WIDTH TOTAL     =',2F8.1,' mA'//)")        & 
     &   EQW,EQWP,EQWT,EQWTP 
      WRITE(16,"(2F12.3,4F12.1)") WLAM(1),WLAM(2),EQW,EQWP,EQWT,EQWTP 
! 
      RETURN 
      END SUBROUTINE OUTPRI 
! 
! ******************************************************************** 
! 
      SUBROUTINE CROSET 
! 
!     SET UP ARRAY CROSS  - PHOTOIONIZATION CROSS-SECTIONS 
! 
      use accura 
      use params 
      use modelp, only : cross,fropc,indexp 
      use synthp 
!     use wincom 
      implicit real(dp) (a-h,o-z),logical (l) 
 
! 
      IJ0=2 
      IF(NFREQ.EQ.1) IJ0=1 
      IF(IMODE.EQ.2) IJ0=NFREQ 
      DO IJ=1,IJ0 
         DO IT=1,MCROSS 
            CROSS(IT,IJ)=0. 
         END DO 
      END DO 
      DO IT=1,NLEVEL 
         IF(INDEXP(IT).NE.5) THEN 
            DO IJ=1,IJ0 
               FR=FREQ(IJ) 
               CROSS(IT,IJ)=SIGK(FR,IT,0) 
            END DO 
         ELSE 
            DO IJ=1,IJ0 
               FR=FREQ(IJ) 
               CROSS(IT,IJ)=SIGK(FR,IT,1) 
               IF(FR.LT.FROPC(IT)) CROSS(IT,IJ)=0. 
            END DO 
         END IF 
      END DO 
! 
      RETURN 
      END SUBROUTINE CROSET 
! 
! ******************************************************************** 
! 
      SUBROUTINE CROSEW 
! 
!     PHOTOIONIZATION CROSS-SECTIONS 
! 
      use accura 
      use params 
      use modelp, only : cross,fropc,indexp 
      use synthp 
      use wincom 
      implicit real(dp) (a-h,o-z),logical (l) 
 
! 
      IJ0=NFREQC 
      DO IJ=1,IJ0 
         DO IT=1,MCROSS 
            CROSS(IT,IJ)=0. 
         END DO 
      END DO 
      DO IT=1,NLEVEL 
         IF(INDEXP(IT).NE.5) THEN 
            DO IJ=1,IJ0 
               FR=FREQC(IJ) 
               CROSS(IT,IJ)=SIGK(FR,IT,0) 
            END DO 
         ELSE 
            DO IJ=1,IJ0 
               FR=FREQC(IJ) 
               CROSS(IT,IJ)=SIGK(FR,IT,1) 
               IF(FR.LT.FROPC(IT)) CROSS(IT,IJ)=0. 
            END DO 
         END IF 
      END DO 
! 
      RETURN 
      END SUBROUTINE CROSEW 
! 
! ******************************************************************** 
! 
! 
 
      FUNCTION SIGK(FR,ITR,MODE) 
!     ========================== 
! 
!     driver for evaluating the photoionization cross-sections 
! 
!     Input: FR  -  frequency 
!            ITR -  index of the transition 
!            mode - =0 cross-section equal to zero longward of edge 
!            mode - >0 cross-section non-zero (extrapolated) longward of edge 
! 
      use accura 
      use params 
      use modelp, only : ctop,xtop,fropc 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: SIH0=2.815e29, E10=2.3025851,                & 
     &                       WI1=911.753878, WI2=227.837832, UN=1.e0 
      REAL(DP) :: XFIT(MFIT) ,                                            & 
     &            SFIT(MFIT)   ! local array containing sigma for OP data 
! 
!     PEACH(X,S,A,B)  =A*X**S*(B+X*(1.-B))*1.E-18 
!     HENRY(X,S,A,B,C)=A*X**S*(C+X*(B-2.*C+X*(1.+C-B)))*1.E-18 
! 
      SIGK=0. 
      II=ITR 
      FR0=ENION(II)/6.6256E-27 
      IF(FR0.LE.0.) RETURN 
      wl0=2.997925e18/fr0 
! 
!     wavelength with an explicit correction to the air wavalength 
! 
      IF(WL0.GT.vaclim) THEN 
         ALM=1.E8/(WL0*WL0) 
         XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM) 
         WL0=WL0/(XN1*1.e-6+UN) 
         fr0=2.997925e18/wl0 
      END IF 
! 
      IF(mode.eq.0  .and. FR.LT.FR0) RETURN 
! 
!     IBF(ITR) is the switch controlling the mode of evaluation of the 
!        cross-section: 
!      = 0  hydrogenic cross-section, with Gaunt factor set to 1 
!      = 1  hydrogenic cross-section with exact Gaunt factor 
!      = 2  Peach-type expression (see function PEACH) 
!      = 3  Henry-type expression (see function HENRY) 
!      = 4  Butler new calculations 
!      = 7  hydrogenic cross-section with Gaunt factor from K. Werner 
!      = 9  Opacity project fits (routine TOPBAS - interpolations) 
!      > 100 - cross-sections extracted form TOPBASE, for several points 
!           In this case, IBF-100 is the number of points 
!      < 0  non-standard, user supplied expression (user should update 
!           subroutine SPSIGK) 
! 
!      for H- : for any IBF > 0  - standard expression 
!      for He I: 
!       for IBF = 11 or = 13  -  Opacity Project cross section 
!                Seaton-Ferney's cubic fits, Hummer's procedure (HEPHOT) 
!           IBF = 11  means that the multiplicity S=1 (singlet) 
!           IBF = 13  means that the multiplicity S=3 (triplet) 
!       for IBF = 10  - cross section, based on Opacity Project, but 
!                       appropriately averaged for an averaged level 
! 
! 
      IB=IBF(ITR) 
      IQ=NQUANT(II) 
      IE=IEL(II) 
      IF(IE.EQ.IELHM) THEN 
         SIGK=SBFHMI(FR) 
         RETURN 
      END IF 
      IF(IE.EQ.IELHE1.AND.IB.GE.10.AND.IB.LE.13) THEN 
         SIGK=SBFHE1(II,IB,FR) 
         RETURN 
      END IF 
! 
      CH=IZ(IE)*IZ(IE) 
      IQ5=IQ*IQ*IQ*IQ*IQ 
! 
      SELECT CASE(IB) 
      CASE(0) 
! 
!        hydrogenic expression (for IBF = 0) 
! 
         SIGK=SIH0/FR/FR/FR*CH*CH/IQ5 
! 
!        exact hydrogenic - with Gaunt factor (for IBF=1) 
! 
       CASE(1) 
         SIGK=SIH0/FR/FR/FR*CH*CH/IQ5 
!        IF(FR.GE.FR0.OR.(IE.EQ.IELH.AND.IQ.LE.3)) 
!    *   SIGK=SIGK*GAUNT(IQ,FR/CH) 
         fr0l=0.95*fr0 
         if(fr.ge.fr0) then 
            sigk=sigk*gaunt(iq,fr/ch) 
          else if(fr.ge.fr0l) then 
            gau0=gaunt(iq,fr0/ch) 
            corg=(fr-fr0l)/(fr0-fr0l)*(gau0-1.)+1. 
            sigk=sigk*corg 
         end if 
       CASE(2) 
! 
!        Peach-type formula (for IBF=2) 
! 
         IF(GAMBF(II).GT.0) THEN 
            IF(GAMBF(II).LT.1.E6) THEN 
              FR0=2.997925E18/GAMBF(II) 
            ELSE 
              FR0=GAMBF(II) 
            END IF 
            IF(FR.LT.FR0) RETURN 
         END IF 
         FREL=FR0/FR 
         SIGK=PEACH(FREL,S0BF(II),ALFBF(II),BETBF(II)) 
       CASE(3) 
! 
!        Henry-type formula (for IBF=3) 
! 
         FREL=FR0/FR 
         SIGK=HENRY(FREL,S0BF(II),ALFBF(II),BETBF(II),GAMBF(II)) 
       CASE(4) 
! 
!     Butler expression 
! 
         FREL=FR0/FR 
         XL=LOG(FREL) 
         SL=S0BF(II)+XL*(ALFBF(II)+XL*BETBF(II)) 
         SIGK=EXP(SL) 
! 
!     exact hydrogenic - with Gaunt factor from K Werner (for IBF=7) 
! 
       CASE(7) 
         IQ5=IQ*IQ*IQ*IQ*IQ 
         SIGK=SIH0/(FR*FR*FR)*CH*CH/IQ5*GNTK(IQ,FR/CH) 
! 
!     selected Opacity Project data (for IBF=9) 
!     (c.-s. evaluated by routine TOPBAS which needs an input file RBF.DAT) 
! 
       CASE(9) 
         SIGK=TOPBAS(FR,FR0,TYPLEV(II)) 
      END SELECT 
! 
!     other Opacity Project data (for IBF>100) 
!     (c.-s. evaluated by interpolating from direct input data) 
! 
      IF(IB.GT.100) THEN 
         NFIT=IB-100 
         X = LOG10(FR/FR0) 
         IF(X.LT.XTOP(1,II)) THEN 
            SIGM=0. 
          ELSE 
            DO IFIT = 1,NFIT 
               XFIT(IFIT) = XTOP(IFIT,II) 
               SFIT(IFIT) = CTOP(IFIT,II) 
            END DO 
            SIGM  = YLINTP (X,XFIT,SFIT,NFIT,MFIT) 
            SIGM  = 1.D-18*EXP(E10*SIGM) 
         END IF 
         SIGK=SIGM 
       ELSE IF(IB.LT.0) THEN 
         CALL SPSIGK(ITR,IB,FR,SIGSP) 
         SIGK=SIGSP 
      END IF 
      if(iatm(ii).eq.iath.and.ii.gt.n0hn+2.and.ib.le.1.and.fr.lt.fr0)     & 
     &   then 
         fr1=fropc(ii) 
         frdec=min(fr1*1.25,fr0) 
         if(fr.gt.fr1.and.fr.lt.frdec)                                    & 
     &      sigk=sigk*(fr-fr1)/(frdec-fr1) 
      end if 
      RETURN 
 
      CONTAINS 
         REAL(DP) FUNCTION PEACH(X,S,A,B) 
            REAL(DP), INTENT(IN) :: X,S,A,B 
            PEACH  =A*X**S*(B+X*(1.-B))*1.e-18 
         END FUNCTION PEACH 
         REAL(DP) FUNCTION HENRY(X,S,A,B,C) 
            REAL(DP), INTENT(IN) :: S,X,A,B,C 
            HENRY=A*X**S*(C+X*(B-2.*C+X*(1.+C-B)))*1.e-18 
         END FUNCTION HENRY 
 
      END FUNCTION SIGK 
! 
! 
!     **************************************************************** 
! 
! 
 
      FUNCTION GAUNT(I,FR) 
!     ==================== 
! 
!     Hydrogenic bound-free Gaunt factor for the principal quantum 
!     number I and frequency FR 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      X=FR/2.99793E14 
      GAUNT=1. 
      IF(I.EQ.1) THEN 
      GAUNT=1.2302628+X*(-2.9094219E-3+X*(7.3993579E-6-8.7356966E-9*X))   & 
     &+(12.803223/X-5.5759888)/X 
       ELSE IF(I.EQ.2) THEN 
      GAUNT=1.1595421+X*(-2.0735860E-3+2.7033384E-6*X)+(-1.2709045+       & 
     &(-2.0244141/X+2.1325684)/X)/X 
       ELSE IF(I.EQ.3) THEN 
      GAUNT=1.1450949+X*(-1.9366592E-3+2.3572356E-6*X)+(-0.55936432+      & 
     &(-0.23387146/X+0.52471924)/X)/X 
       ELSE IF(I.EQ.4) THEN 
      GAUNT=1.1306695+X*(-1.3482273E-3+X*(-4.6949424E-6+2.3548636E-8*X))  & 
     &+(-0.31190730+(0.19683564-5.4418565E-2/X)/X)/X 
       ELSE IF(I.EQ.5) THEN 
      GAUNT=1.1190904+X*(-1.0401085E-3+X*(-6.9943488E-6+2.8496742E-8*X))  & 
     &+(-0.16051018+(5.5545091E-2-8.9182854E-3/X)/X)/X 
       ELSE IF(I.EQ.6) THEN 
      GAUNT=1.1168376+X*(-8.9466573E-4+X*(-8.8393133E-6+3.4696768E-8*X))  & 
     &+(-0.13075417+(4.1921183E-2-5.5303574E-3/X)/X)/X 
       ELSE IF(I.EQ.7) THEN 
      GAUNT=1.1128632+X*(-7.4833260E-4+X*(-1.0244504E-5+3.8595771E-8*X))  & 
     &+(-9.5441161E-2+(2.3350812E-2-2.2752881E-3/X)/X)/X 
       ELSE IF(I.EQ.8) THEN 
      GAUNT=1.1093137+X*(-6.2619148E-4+X*(-1.1342068E-5+4.1477731E-8*X))  & 
     &+(-7.1010560E-2+(1.3298411E-2 -9.7200274E-4/X)/X)/X 
       ELSE IF(I.EQ.9) THEN 
      GAUNT=1.1078717+X*(-5.4837392E-4+X*(-1.2157943E-5+4.3796716E-8*X))  & 
     &+(-5.6046560E-2+(8.5139736E-3-4.9576163E-4/X)/X)/X 
       ELSE IF(I.EQ.10) THEN 
      GAUNT=1.1052734+X*(-4.4341570E-4+X*(-1.3235905E-5+4.7003140E-8*X))  & 
     &+(-4.7326370E-2+(6.1516856E-3-2.9467046E-4/X)/X)/X 
      END IF 
      RETURN 
      END FUNCTION GAUNT 
! 
! 
!     **************************************************************** 
! 
! 
 
      FUNCTION GNTK(I,FR) 
!     =================== 
! 
!     Hydrogenic bound-free Gaunt factor for the principal quantum 
!     number I and frequency FR (from Klaus Werner) 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      Y=1./FR 
      IF(I.EQ.1) THEN 
         GNTK=0.9916+Y*(2.71852e13-Y*2.26846e30) 
       ELSE IF(I.EQ.2) THEN 
         GNTK=1.1050-Y*(2.37490e14-Y*4.07677e28) 
       ELSE IF (I.EQ.3) THEN 
         GNTK=1.1010-Y*(0.98632e14-Y*1.03540e28) 
       ELSE 
         GNTK=1. 
      END IF 
      END FUNCTION GNTK 
! 
! 
!     **************************************************************** 
! 
! 
 
      SUBROUTINE SPSIGK(ITR,IB,FR,SIGSP) 
!     ================================== 
! 
!     Non-standard evaluation of the photoionization cross-sections 
!     Basically user-suppled procedure; here are some examples 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      SIGSP=0. 
      if(itr.le.0) return 
! 
!     Special formula for the He I ground state 
! 
      IF(IB.EQ.-201) SIGSP=7.3E-18*EXP(1.373-2.311E-16*FR) 
! 
!     Special formula for the averaged <n=2> level of He I 
! 
      IF(IB.EQ.-202) SIGSP=SGHE12(FR) 
! 
!     Carbon ground configuration levels 2p2 1D and 1S 
! 
      IF(IB.EQ.-602.OR.IB.EQ.-603) THEN 
         CALL CARBON(IB,FR,SG) 
         SIGSP=SG 
      END IF 
! 
!     Hidalgo (Ap.J. 153, 981, 1968) photoionization data 
! 
      IF(IB.LE.-101.AND.IB.GE.-137) SIGSP=HIDALG(IB,FR) 
! 
!     Reilman and Manson (Ap.J. Suppl. 40, 815, 1979) photoionization data 
! 
      IF(IB.LE.-301.AND.IB.GE.-337) SIGSP=REIMAN(IB,FR) 
      RETURN 
      END SUBROUTINE SPSIGK 
! 
! 
! 
!     **************************************************************** 
! 
! 
 
      SUBROUTINE CARBON(IB,FR,SG) 
!     =========================== 
! 
!     Photoionization cross-section for neutral carbon 2p1D and 2p1S 
!     levels (G.B.Taylor - private communication) 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) ::  FR2(34),SG2(34),FR3(45),SG3(45) 
      DATA FR2/ 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82,     & 
     &                0.83,       0.85, 0.86, 0.87, 0.88, 0.89, 0.90,     & 
     &          0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99,     & 
     &          1.00, 1.10, 1.20, 1.30, 1.45, 1.50, 1.60, 1.80, 2./ 
      DATA SG2/ 12.04, 12.03, 12.09, 12.26, 12.60, 13.24, 14.36, 16.24,   & 
     &          19.28, 23.94, 37.41, 42.88, 44.76, 43.41, 40.46, 37.19,   & 
     &          34.26, 31.82, 29.96, 28.57, 27.68, 27.37, 27.84, 29.69,   & 
     &          34.45, 46.35, 13.80, 11.54, 10.40,  8.96,  8.54,  7.47,   & 
     &           6.53,  5.66/ 
      DATA FR3/ 0.66, 0.68, 0.70, 0.72, 0.74, 0.76, 0.78, 0.80, 0.82,     & 
     &          0.84, 0.86, 0.864,0.866,0.868,0.87, 0.874,0.876,0.88,     & 
     &          0.882,0.884,0.886,0.888,0.89 ,0.894,0.896,0.898,0.90,     & 
     &          0.904,0.908,0.910,0.920,0.94, 0.98, 1.00, 1.10, 1.20,     & 
     &          1.26, 1.34, 1.36, 1.40, 1.46, 1.60, 1.70, 1.80, 2./ 
      DATA SG3/ 13.94, 13.29, 12.56, 11.73, 10.82, 10.18,  8.62,  7.27,   & 
     &           5.74,  4.14,  4.61,  5.92,  6.94,  8.34, 10.21, 16.12,   & 
     &          20.64, 34.56, 44.82, 57.71, 73.09, 89.99,106.38,127.08,   & 
     &         128.38,124.44,117.17, 99.32, 82.95, 76.05, 52.65, 33.23,   & 
     &          21.29, 18.69, 12.62, 11.44,  9.77,  7.53, 10.47,  9.65,   & 
     &          10.19,  7.28,  6.70,  6.11,  4.96/ 
      DATA NC2,NC3/34,45/ 
      DATA FR0/3.28805E15/ 
      F=FR/FR0 
      IF(IB.EQ.-602) THEN 
         J=2 
         IF(F.GT.FR2(1)) THEN 
            DO I=2,NC2 
               J=I 
               IF(F.GT.FR2(I-1).AND.F.LE.FR2(I)) EXIT 
            END DO 
         END IF 
         SG=(F-FR2(J-1))/(FR2(J)-FR2(J-1))*(SG2(J)-SG2(J-1))+SG2(J-1) 
         SG=SG*1.E-18 
       ELSE IF(IB.EQ.-603) THEN 
         J=2 
         IF(F.GT.FR3(1)) THEN 
            DO I=2,NC3 
               J=I 
               IF(F.GT.FR3(I-1).AND.F.LE.FR3(I)) EXIT 
            END DO 
          END IF 
         SG=(F-FR3(J-1))/(FR3(J)-FR3(J-1))*(SG3(J)-SG3(J-1))+SG3(J-1) 
         SG=SG*1.E-18 
      END IF 
 
      RETURN 
      END SUBROUTINE CARBON 
! 
! 
!     **************************************************************** 
! 
 
      FUNCTION SGHE12(FR) 
!     =================== 
! 
!     Special formula for the photoionization cross-section from the 
!     averaged <n=2> level of He I 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      DATA C1/3.E0/,C2/9.E0/,C3/1.6E1/,                                   & 
     & A1/6.45105E-18/,A2/3.02E-19/,A3/9.9847E-18/,A4/1.1763673E-17/,     & 
     & A5/3.63662E-19/,A6/-2.783E2/,A7/1.488E1/,A8/-2.311E-1/,            & 
     & E1/3.5E0/,E2/3.6E0/,E3/1.91E0/,E4/2.9E0/,E5/3.3E0/ 
      X=FR*1.E-15 
      XX=LOG(FR) 
      SGHE12=(C1*(A1/X**E1+A2/X**E2)+A3/X**E3+C2*(A4/X**E4+A5/X**E5)+     & 
     &       C1*EXP(A6+XX*(A7+XX*A8)))/C3 
      RETURN 
      END FUNCTION SGHE12 
! 
! 
!     **************************************************************** 
! 
! 
      FUNCTION HIDALG(IB,FR) 
!     ====================== 
! 
!     Read table of wavelengths and photo-ionization cross-sections 
!     from Hidalgo (1968, Ap. J., 153, 981) for the species indicated by IB 
!     (Hidalgo's number = INDEX = -IB-100). 
!     Compute linearly interpolated value of the cross-section 
!     at the frequency FR. 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) ::  WL1(20),WL2(20),WLI(20),SIG0(20,24),SIGS(20) 
! 
      DATA WL1 /                                                          & 
     &  39.1, 80.9, 97.6,100.1,104.3,107.2,108.7,111.9,113.6,115.4,       & 
     & 117.1,119.0,124.8,126.9,129.1,131.3,133.6,136.0,138.5,141.1/ 
      DATA WL2 /                                                          & 
     &  68.5, 80.9,100.1,120.9,158.8,165.7,177.3,190.6,200.7,206.2,       & 
     & 211.9,218.0,224.5,231.3,246.3,5*0./ 
      DATA SIG0 /                                                         & 
     &120*0.,                                                             & 
     &.0460,.2400,.3500,.3700,.4000,.4300,.4400,.4600,.4700,.4900,        & 
     &.5000,.5200,.5700,.6200, 6*0.,                                      & 
     & 80*0.,                                                             & 
     &.0092,.1000,.1900,.2100,.2300,.2500,.2600,.2900,.3000,.3200,        & 
     &.3400,.3500,.4100,.4300,.4500,.4800,.5000,.5300,.5600,.5900,        & 
     & 20*0.,                                                             & 
     &.3400,.4600,.6300,.7700,.9100,1.080, 14*0.,                         & 
     & 20*0.,                                                             & 
     &.0064,.1100,.2200,.4100,.9400,1.000,1.300,1.600, 12*0.,             & 
     & 80*0.,                                                             & 
     &.0370,.0650,.1300,.2400,.5500,.6300,.7700,.9500,1.100,1.250,        & 
     & 10*0.,                                                             & 
     & 40*0.,                                                             & 
     &.0220,.0390,.0800,.1500,.3500,.4000,.4900,.6200,.7200,.7800,        & 
     &.8500,.9300,1.020,                                                  & 
     & 7*0./ 
! 
      INDEX=-IB-100 
      NUM=20 
      IF(INDEX.GE.13.AND.INDEX.LE.27) NUM=15 
      DO I=1,NUM 
         IF(INDEX.LT.13) WLI(I)=WL1(I) 
         IF(INDEX.GE.13) WLI(I)=WL2(I) 
         SIGS(I)=SIG0(I,INDEX) 
      END DO 
! 
      WLAM=2.997925E18/FR 
      IL=1 
      IR=NUM 
      INTR: DO I=1,NUM-1 
         IF(WLAM.GE.WLI(I).AND.WLAM.LE.WLI(I+1)) THEN 
            IL=I 
            IR=I+1 
            EXIT INTR 
         END IF 
      END DO INTR 
! 
!     LINEAR INTERPOLATION: 
! 
      SIGM=(SIGS(IR)-SIGS(IL))*(WLAM-WLI(IL))/(WLI(IR)-WLI(IL))           & 
     &      + SIGS(IL) 
! 
!     IF OUTSIDE WAVELENGTH RANGE SET TO FIRST(LAST) VALUE: 
! 
       IF(WLAM.LE.WLI(1)) SIGM=SIGS(1) 
       IF(WLAM.GE.WLI(NUM)) SIGM=SIGS(NUM) 
! 
!     IF LAST NON-ZERO SIG VALUES, NO INTERPOLATION: 
! 
!       IF(SIGS(IR).EQ.0.) SIGM=SIGS(IL) 
! 
      HIDALG=SIGM*1.E-18 
      RETURN 
      END FUNCTION HIDALG 
! 
! 
!     **************************************************************** 
! 
! 
      FUNCTION REIMAN(IB,FR) 
!     ====================== 
! 
!     Read table of photon energies and photo-ionization cross-sections 
!     from Reilman & Manson (1979, Ap. J. Suppl., 40, 815) for the species 
!     indicated by IB 
! 
!     Compute linearly interpolated value of the cross-section 
!     at the frequency FR. 
! 
!     (At the moment, only a few transitions are considered) 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: HEV(30),F0(30),SIG0(30,2),SIGS(30) 
! 
      DATA HEV /                                                          & 
     & 130.,160.,190.,210.,240.,270.,300.,330.,360.,390.,                 & 
     & 420.,450.,480.,510.,540.,570.,600.,630.,660.,690.,                 & 
     & 720.,750.,780.,810.,840.,870.,900.,930.,960.,990./ 
      DATA SIG0 /                                                         & 
     & 3*0.,     4.422E-1, 3.478E-1,                                      & 
     & 2.794E-1, 2.286E-1, 1.899E-1, 1.598E-1, 1.360E-1,                  & 
     & 1.169E-1, 1.013E-1, 8.845E-2, 7.776E-2, 6.877E-2,                  & 
     & 6.114E-2, 5.463E-2, 4.904E-2, 4.419E-2, 3.998E-2,                  & 
     & 3.629E-2, 3.305E-2, 3.019E-2, 2.766E-2, 2.540E-2,                  & 
     & 2.339E-2, 2.158E-2, 1.996E-2, 1.850E-2, 1.718E-2,                  & 
     & 4*0.,     1.981E-1, 1.584E-1,                                      & 
     & 1.290E-1, 1.066E-1, 8.932E-2, 7.567E-2, 6.475E-2,                  & 
     & 5.589E-2, 4.862E-2, 4.259E-2, 3.754E-2, 3.329E-2,                  & 
     & 2.966E-2, 2.656E-2, 2.388E-2, 2.157E-2, 1.954E-2,                  & 
     & 1.777E-2, 1.621E-2, 1.484E-2, 1.362E-2, 1.253E-2,                  & 
     & 1.155E-2, 1.067E-2, 9.888E-3, 9.179E-3/ 
! 
      INDEX=-IB-300 
      NUM=30 
      DO I=1,NUM 
         F0(I)=HEV(I)*2.418573E14 
         SIGS(I)=SIG0(I,INDEX) 
      END DO 
! 
      IL=1 
      IR=NUM 
      INTR: DO I=1,NUM-1 
         IF(FR.GE.F0(I).AND.FR.LE.F0(I+1)) THEN 
            IL=I 
            IR=I+1 
            EXIT INTR 
         END IF 
      END DO INTR 
! 
!     LINEAR INTERPOLATION: 
! 
      SIGM=(SIGS(IR)-SIGS(IL))*(FR-F0(IL))/(F0(IR)-F0(IL))                & 
     &      + SIGS(IL) 
! 
!     IF OUTSIDE WAVELENGTH RANGE SET TO FIRST(LAST) VALUE: 
! 
       IF(FR.LE.F0(1)) SIGM=SIGS(1) 
       IF(FR.GE.F0(NUM)) SIGM=SIGS(NUM) 
! 
!     IF LAST NON-ZERO SIG VALUES, NO INTERPOLATION: 
! 
!       IF(SIGS(IR).EQ.0.) SIGM=SIGS(IL) 
! 
      REIMAN=SIGM*1.E-18 
      RETURN 
      END FUNCTION REIMAN 
! 
! 
!     **************************************************************** 
! 
! 
 
      FUNCTION SBFHE1(II,IB,FR) 
!     ========================= 
! 
!     Calculates photoionization cross sections of neutral helium 
!     from states with n = 1, 2, 3, 4. 
! 
!     The levels are either non-averaged (l,s) states, or some 
!     averaged levels. 
!     The program allows only two standard possibilities of 
!     constructing averaged levels: 
!     i)  all states within given principal quantum number n (>1) are 
!         lumped together 
!     ii) all siglet states for given n, and all triplet states for 
!         given n are lumped together separately (there are thus two 
!         explicit levels for a given n) 
! 
!     The cross sections are calculated using appropriate averages 
!     of the Opacity Project cross sections, calculated by procedure 
!     HEPHOT 
! 
!     Input parameters: 
!      II    - index of the lower level (in the numbering of explicit 
!              levels) 
!      IB    - photoionization switch IBF for the given transition 
!            = 10  -  means that the given transition is from an 
!                     averaged level 
!            = 11  -  the given transition is from non-averaged 
!                     singlet state 
!            = 13  -  the given transition is from non-averaged 
!                     triplet state 
!      FR    - frequency 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
! 
      NI=NQUANT(II) 
      IGI=INT(G(II)+0.01) 
      IS=IB-10 
      sbfhe1=0. 
! 
!     ---------------------------------------------------------------- 
!     IB=11 or 13  - photoionization from an non-averaged (l,s) level 
!     ---------------------------------------------------------------- 
! 
      IF(IS.EQ.1.OR.IS.EQ.3) THEN 
         IL=(IGI/IS-1)/2 
         SBFHE1=HEPHOT(IS,IL,NI,FR) 
      END IF 
! 
!     ---------------------------------------------------------------- 
!     IS=10 - photoionization from an averaged level 
!     ---------------------------------------------------------------- 
! 
      IF(IS.EQ.0) THEN 
         IF(NI.EQ.2) THEN 
! 
! ********    photoionization from an averaged level with n=2 
! 
            IF(IGI.EQ.4) THEN 
! 
!      a) lower level is an averaged singlet state 
! 
               SBFHE1=(HEPHOT(1,0,2,FR)+3.*HEPHOT(1,1,2,FR))/9. 
            ELSE IF(IGI.EQ.12) THEN 
! 
!      b) lower level is an averaged triplet state 
! 
               SBFHE1=(HEPHOT(3,0,2,FR)+3.*HEPHOT(3,1,2,FR))/9. 
            ELSE IF(IGI.EQ.16) THEN 
! 
!      c) lower level is an average of both singlet and triplet states 
! 
               SBFHE1=(HEPHOT(1,0,2,FR)+3.*(HEPHOT(1,1,2,FR)+             & 
     &                HEPHOT(3,0,2,FR))+9.*HEPHOT(3,1,2,FR))/1.6e1 
            ELSE 
              call quit('hephot error') 
            END IF 
! 
! 
! ********    photoionization from an averaged level with n=3 
! 
         ELSE IF(NI.EQ.3) THEN 
            IF(IGI.EQ.9) THEN 
! 
!      a) lower level is an averaged singlet state 
! 
               SBFHE1=(HEPHOT(1,0,3,FR)+3.*HEPHOT(1,1,3,FR)+              & 
     &                5.*HEPHOT(1,2,3,FR))/9. 
            ELSE IF(IGI.EQ.27) THEN 
! 
!      b) lower level is an averaged triplet state 
! 
               SBFHE1=(HEPHOT(3,0,3,FR)+3.*HEPHOT(3,1,3,FR)+              & 
     &                5.*HEPHOT(3,2,3,FR))/9. 
            ELSE IF(IGI.EQ.36) THEN 
! 
!      c) lower level is an average of both singlet and triplet states 
! 
               SBFHE1=(HEPHOT(1,0,3,FR)+3.*HEPHOT(1,1,3,FR)+              & 
     &                5.*HEPHOT(1,2,3,FR)+                                & 
     &                3.*HEPHOT(3,0,3,FR)+9.*HEPHOT(3,1,3,FR)+            & 
     &                15.*HEPHOT(3,2,3,FR))/3.6 
            ELSE 
               call quit('hephot error') 
            END IF 
         ELSE IF(NI.EQ.4) THEN 
! 
! ********    photoionization from an averaged level with n=4 
! 
            IF(IGI.EQ.16) THEN 
! 
!      a) lower level is an averaged singlet state 
! 
               SBFHE1=(HEPHOT(1,0,4,FR)+3.*HEPHOT(1,1,4,FR)+              & 
     &                 5.0*HEPHOT(1,2,4,FR)+                              & 
     &                 7.0*HEPHOT(1,3,4,FR))/16. 
            ELSE IF(IGI.EQ.48) THEN 
! 
!      b) lower level is an averaged triplet state 
! 
               SBFHE1=(HEPHOT(3,0,4,FR)+3.*HEPHOT(3,1,4,FR)+              & 
     &                 5.*HEPHOT(3,2,4,FR)+                               & 
     &                 7.*HEPHOT(3,3,4,FR))/16. 
            ELSE IF(IGI.EQ.64) THEN 
! 
!      c) lower level is an average of both singlet and triplet states 
! 
               SBFHE1=(HEPHOT(1,0,4,FR)+3.D0*HEPHOT(1,1,4,FR)+            & 
     &                 5.*HEPHOT(1,2,4,FR)+                               & 
     &                 7.*HEPHOT(1,3,4,FR)+                               & 
     &                 3.*HEPHOT(3,0,4,FR)+                               & 
     &                 9.*HEPHOT(3,1,4,FR)+                               & 
     &                 15.*HEPHOT(3,2,4,FR)+                              & 
     &                 21.*HEPHOT(3,3,4,FR))/64. 
            ELSE 
               call quit('hephot error') 
            END IF 
         ELSE 
            call quit('hephot error') 
         END IF 
      END IF 
      RETURN 
      END FUNCTION SBFHE1 
! 
! 
!     **************************************************************** 
! 
! 
 
      FUNCTION HEPHOT(S,L,N,FREQ) 
!     =========================== 
! 
!           EVALUATES HE I PHOTOIONIZATION CROSS SECTION USING SEATON 
!           FERNLEY'S CUBIC FITS TO THE OPACITY PROJECT CROSS SECTIONS 
!           UP TO SOME ENERGY "EFITM" IN THE RESONANCE-FREE ZONE.  BEYOND 
!           THIS ENERGY LINEAR FITS TO LOG SIGMA IN LOG (E/E0) ARE USED. 
!           THIS EXTRAPOLATION SHOULD BE USED UP TO THE BEGINNING OF THE 
!           RESONANCE ZONE "XMAX", BUT AT PRESENT IT IS USED THROUGH IT. 
!           BY CHANGING A FEW LINES THAT ARE PRESENTLY COMMENTED OUT, 
!           FOR ENERGIES IN THE RESONANCE ZONE A VALUE OF 1/100 OF THE 
!           THRESHOLD CROSS SECTION IS USED -- THIS IS PURELY AD HOC AND 
!           ONLY A TEMPORARY MEASURE.  OBVIOUSLY ANY OTHER VALUE OR FUNCTIONAL 
!           FORM CAN BE INSERTED HERE. 
! 
!           CALLING SEQUENCE INCLUDES: 
!                S = MULTIPLICITY, EITHER 1 OR 3 
!                L = ANGULAR MOMENTUM, 0, 1, OR 2; 
!                    for L > 2 - hydrogenic expresion 
!                FREQ = FREQUENCY 
! 
!           DGH JUNE 1988 JILA, slightly modified by I.H. 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      INTEGER  :: S,L,SS,LL 
      REAL(DP) :: COEF(4,53),FL0(53),A(53),B(53),XFITM(53) 
      INTEGER  :: N0(3,2),IST(3,2) 
! 
      DATA IST/1,36,20,11,45,28/ 
      DATA N0/1,2,3,2,2,3/ 
! 
      DATA FL0/                                                           & 
     & 2.521e-01,-5.381e-01,-9.139e-01,-1.175e+00,-1.375e+00,-1.537e+00,  & 
     &-1.674e+00,-1.792e+00,-1.896e+00,-1.989e+00,-4.555e-01,-8.622e-01,  & 
     &-1.137e+00,-1.345e+00,-1.512e+00,-1.653e+00,-1.774e+00,-1.880e+00,  & 
     &-1.974e+00,-9.538e-01,-1.204e+00,-1.398e+00,-1.556e+00,-1.690e+00,  & 
     &-1.806e+00,-1.909e+00,-2.000e+00,-9.537e-01,-1.204e+00,-1.398e+00,  & 
     &-1.556e+00,-1.690e+00,-1.806e+00,-1.909e+00,-2.000e+00,-6.065e-01,  & 
     &-9.578e-01,-1.207e+00,-1.400e+00,-1.558e+00,-1.692e+00,-1.808e+00,  & 
     &-1.910e+00,-2.002e+00,-5.749e-01,-9.352e-01,-1.190e+00,-1.386e+00,  & 
     &-1.547e+00,-1.682e+00,-1.799e+00,-1.902e+00,-1.995e+00/ 
! 
      DATA XFITM/                                                         & 
     & 3.262e-01, 6.135e-01, 9.233e-01, 8.438e-01, 1.020e+00, 1.169e+00,  & 
     & 1.298e+00, 1.411e+00, 1.512e+00, 1.602e+00, 7.228e-01, 1.076e+00,  & 
     & 1.206e+00, 1.404e+00, 1.481e+00, 1.464e+00, 1.581e+00, 1.685e+00,  & 
     & 1.777e+00, 9.586e-01, 1.187e+00, 1.371e+00, 1.524e+00, 1.740e+00,  & 
     & 1.854e+00, 1.955e+00, 2.046e+00, 9.585e-01, 1.041e+00, 1.371e+00,  & 
     & 1.608e+00, 1.739e+00, 1.768e+00, 1.869e+00, 1.803e+00, 7.360e-01,  & 
     & 1.041e+00, 1.272e+00, 1.457e+00, 1.611e+00, 1.741e+00, 1.855e+00,  & 
     & 1.870e+00, 1.804e+00, 9.302e-01, 1.144e+00, 1.028e+00, 1.210e+00,  & 
     & 1.362e+00, 1.646e+00, 1.761e+00, 1.863e+00, 1.954e+00/ 
! 
      DATA A/                                                             & 
     & 6.95319e-01, 1.13101e+00, 1.36313e+00, 1.51684e+00, 1.64767e+00,   & 
     & 1.75643e+00, 1.84458e+00, 1.87243e+00, 1.85628e+00, 1.90889e+00,   & 
     & 9.01802e-01, 1.25389e+00, 1.39033e+00, 1.55226e+00, 1.60658e+00,   & 
     & 1.65930e+00, 1.68855e+00, 1.62477e+00, 1.66726e+00, 1.83599e+00,   & 
     & 2.50403e+00, 3.08564e+00, 3.56545e+00, 4.25922e+00, 4.61346e+00,   & 
     & 4.91417e+00, 5.19211e+00, 1.74181e+00, 2.25756e+00, 2.95625e+00,   & 
     & 3.65899e+00, 4.04397e+00, 4.13410e+00, 4.43538e+00, 4.19583e+00,   & 
     & 1.79027e+00, 2.23543e+00, 2.63942e+00, 3.02461e+00, 3.35018e+00,   & 
     & 3.62067e+00, 3.85218e+00, 3.76689e+00, 3.49318e+00, 1.16294e+00,   & 
     & 1.86467e+00, 2.02110e+00, 2.24231e+00, 2.44240e+00, 2.76594e+00,   & 
     & 2.93230e+00, 3.08109e+00, 3.21069e+00/ 
! 
      DATA B/                                                             & 
     &-1.29000e+00,-2.15771e+00,-2.13263e+00,-2.10272e+00,-2.10861e+00,   & 
     &-2.11507e+00,-2.11710e+00,-2.08531e+00,-2.03296e+00,-2.03441e+00,   & 
     &-1.85905e+00,-2.04057e+00,-2.02189e+00,-2.05930e+00,-2.03403e+00,   & 
     &-2.02071e+00,-1.99956e+00,-1.92851e+00,-1.92905e+00,-4.58608e+00,   & 
     &-4.40022e+00,-4.39154e+00,-4.39676e+00,-4.57631e+00,-4.57120e+00,   & 
     &-4.56188e+00,-4.55915e+00,-4.41218e+00,-4.12940e+00,-4.24401e+00,   & 
     &-4.40783e+00,-4.39930e+00,-4.25981e+00,-4.26804e+00,-4.00419e+00,   & 
     &-4.47251e+00,-3.87960e+00,-3.71668e+00,-3.68461e+00,-3.67173e+00,   & 
     &-3.65991e+00,-3.64968e+00,-3.48666e+00,-3.23985e+00,-2.95758e+00,   & 
     &-3.07110e+00,-2.87157e+00,-2.83137e+00,-2.82132e+00,-2.91084e+00,   & 
     &-2.91159e+00,-2.91336e+00,-2.91296e+00/ 
! 
      DATA ((COEF(I,J),I=1,4),J=1,10)/                                    & 
     & 8.734e-01,-1.545e+00,-1.093e+00, 5.918e-01, 9.771e-01,-1.567e+00,  & 
     &-4.739e-01,-1.302e-01, 1.174e+00,-1.638e+00,-2.831e-01,-3.281e-02,  & 
     & 1.324e+00,-1.692e+00,-2.916e-01, 9.027e-02, 1.445e+00,-1.761e+00,  & 
     &-1.902e-01, 4.401e-02, 1.546e+00,-1.817e+00,-1.278e-01, 2.293e-02,  & 
     & 1.635e+00,-1.864e+00,-8.252e-02, 9.854e-03, 1.712e+00,-1.903e+00,  & 
     &-5.206e-02, 2.892e-03, 1.782e+00,-1.936e+00,-2.952e-02,-1.405e-03,  & 
     & 1.845e+00,-1.964e+00,-1.152e-02,-4.487e-03/ 
      DATA ((COEF(I,J),I=1,4),J=11,19)/                                   & 
     & 7.377e-01,-9.327e-01,-1.466e+00, 6.891e-01, 9.031e-01,-1.157e+00,  & 
     &-7.151e-01, 1.832e-01, 1.031e+00,-1.313e+00,-4.517e-01, 9.207e-02,  & 
     & 1.135e+00,-1.441e+00,-2.724e-01, 3.105e-02, 1.225e+00,-1.536e+00,  & 
     &-1.725e-01, 7.191e-03, 1.302e+00,-1.602e+00,-1.300e-01, 7.345e-03,  & 
     & 1.372e+00,-1.664e+00,-8.204e-02,-1.643e-03, 1.434e+00,-1.715e+00,  & 
     &-4.646e-02,-7.456e-03, 1.491e+00,-1.760e+00,-1.838e-02,-1.152e-02/ 
      DATA ((COEF(I,J),I=1,4),J=20,27)/                                   & 
     & 1.258e+00,-3.442e+00,-4.731e-01,-9.522e-02, 1.553e+00,-2.781e+00,  & 
     &-6.841e-01,-4.083e-03, 1.727e+00,-2.494e+00,-5.785e-01,-6.015e-02,  & 
     & 1.853e+00,-2.347e+00,-4.611e-01,-9.615e-02, 1.955e+00,-2.273e+00,  & 
     &-3.457e-01,-1.245e-01, 2.041e+00,-2.226e+00,-2.669e-01,-1.344e-01,  & 
     & 2.115e+00,-2.200e+00,-1.999e-01,-1.410e-01, 2.182e+00,-2.188e+00,  & 
     &-1.405e-01,-1.460e-01/ 
      DATA ((COEF(I,J),I=1,4),J=28,35)/                                   & 
     & 1.267e+00,-3.417e+00,-5.038e-01,-1.797e-02, 1.565e+00,-2.781e+00,  & 
     &-6.497e-01,-5.979e-03, 1.741e+00,-2.479e+00,-6.099e-01,-2.227e-02,  & 
     & 1.870e+00,-2.336e+00,-4.899e-01,-6.616e-02, 1.973e+00,-2.253e+00,  & 
     &-3.972e-01,-8.729e-02, 2.061e+00,-2.212e+00,-3.072e-01,-1.060e-01,  & 
     & 2.137e+00,-2.189e+00,-2.352e-01,-1.171e-01, 2.205e+00,-2.186e+00,  & 
     &-1.621e-01,-1.296e-01/ 
      DATA ((COEF(I,J),I=1,4),J=36,44)/                                   & 
     & 1.129e+00,-3.149e+00,-1.910e-01,-5.244e-01, 1.431e+00,-2.511e+00,  & 
     &-3.710e-01,-1.933e-01, 1.620e+00,-2.303e+00,-3.045e-01,-1.391e-01,  & 
     & 1.763e+00,-2.235e+00,-1.829e-01,-1.491e-01, 1.879e+00,-2.215e+00,  & 
     &-9.003e-02,-1.537e-01, 1.978e+00,-2.213e+00,-2.066e-02,-1.541e-01,  & 
     & 2.064e+00,-2.220e+00, 3.258e-02,-1.527e-01, 2.140e+00,-2.225e+00,  & 
     & 6.311e-02,-1.455e-01, 2.208e+00,-2.229e+00, 7.977e-02,-1.357e-01/ 
      DATA ((COEF(I,J),I=1,4),J=45,53)/                                   & 
     & 1.204e+00,-2.809e+00,-3.094e-01, 1.100e-01, 1.455e+00,-2.254e+00,  & 
     &-4.795e-01, 6.872e-02, 1.619e+00,-2.109e+00,-3.357e-01,-2.532e-02,  & 
     & 1.747e+00,-2.065e+00,-2.317e-01,-5.224e-02, 1.853e+00,-2.058e+00,  & 
     &-1.517e-01,-6.647e-02, 1.943e+00,-2.055e+00,-1.158e-01,-6.081e-02,  & 
     & 2.023e+00,-2.070e+00,-6.470e-02,-6.800e-02, 2.095e+00,-2.088e+00,  & 
     &-2.357e-02,-7.250e-02, 2.160e+00,-2.107e+00, 1.065e-02,-7.542e-02/ 
! 
!     Hydrogenic expression for L > 2 
!      [multiplied by relative population of state (s,l,n), ie. 
!       by  stat.weight(s,l)/stat.weight(n)] 
! 
      IF(L.GT.2) THEN 
         GN=2.*N*N 
         HEPHOT=2.815e29/FREQ/FREQ/FREQ/N**5*(2*L+1)*S/GN 
         RETURN 
      END IF 
! 
!     SELECT BEGINNING AND END OF COEFFICIENTS 
! 
      SS=(S+1)/2 
      LL=L+1 
      NSL0=N0(LL,SS) 
      I=IST(LL,SS)+N-NSL0 
! 
!     EVALUATE CROSS SECTION 
! 
      FL=LOG10(FREQ/3.28805E15) 
      X=FL-FL0(I) 
      X=max(x,0.) 
      IF(X.GE.-0.001) THEN 
         IF(X.LT.XFITM(I)) THEN 
            P=COEF(4,I) 
            DO K=1,3 
               P=X*P+COEF(4-K,I) 
            END DO 
            HEPHOT=1.e-18*10.**P 
          ELSE 
            HEPHOT=1.e-18*10.**(A(I)+B(I)*X) 
          END IF 
       ELSE 
          HEPHOT=0. 
      END IF 
      RETURN 
      END FUNCTION HEPHOT 
! 
! 
!     **************************************************************** 
! 
! 
 
      FUNCTION TOPBAS(FREQ,FREQ0,TYPLV) 
!     ================================== 
! 
!     Procedure calculates the photo-ionisation cross section SIGMA in 
!     [cm^2] at frequency FREQ. FREQ0 is the threshold frequency from 
!     level I of ion KI. Threshold cross-sections will be of the order 
!     of the numerical value of 10^-18. 
!     Opacity-Project (OP) interpolation fit formula 
! 
      use accura 
      use params 
      use topdat 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: E10=2.3025851 
      CHARACTER(LEN=10)   ::  TYPLV 
      REAL(DP)            :: XFIT(MOP),SFIT(MOP) 
! 
!     Read OP data if not yet done 
! 
      TOPBAS=0. 
      IF (.NOT.LOPREA) CALL OPDATA 
      X = LOG10(FREQ/FREQ0) 
      TOPD: DO IOP = 1,NTOTOP 
         IF (IDLVOP(IOP).EQ.TYPLV) THEN 
!           Level is not found ,or no data for this level, in RBF.DAT 
            IF (NOP(IOP).LE.0) THEN 
               WRITE(61,"('SIGMA.......: OP DATA NOT AVAILABLE',          & 
     &         ' FOR LEVEL ',A10)") TYPLV 
               RETURN 
            END IF 
!           level has been detected in OP-data file 
 
            DO IFIT = 1,NOP(IOP) 
               XFIT(IFIT) = XOP(IFIT,IOP) 
               SFIT(IFIT) = SOP(IFIT,IOP) 
            END DO 
            SIGM  = YLINTP (X,XFIT,SFIT,NOP(IOP),MOP) 
            SIGM  = 1.e-18*EXP(E10*SIGM) 
            TOPBAS=SIGM 
            EXIT TOPD 
         END IF 
      END DO TOPD 
      RETURN 
      END FUNCTION TOPBAS 
! 
 
!     ****************************************************************** 
! 
! 
      SUBROUTINE OPDATA 
!     ================= 
! 
!     Procedure reads photo-ionization cross sections fit coefficients 
!     based on Opacity-Project (OP) data from file RBF.DAT 
!     Data, as stored, requires linear interpolation. 
! 
!     Meaning of global variables: 
!        NTOTOP    = total number of levels in Opacity Project data 
!        IDLVOP() = level identifyer of current level 
!        NOP()     = number of fit points for current level 
!        XOP(,)    = x     = alog10(nu/nu0)       of fit point 
!        SOP(,)    = sigma = alog10(sigma/10^-18) of fit point 
! 
      use accura 
      use params 
      use topdat 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      CHARACTER(LEN=4) :: IONID 
! 
      CALL ALLOC_TOPDAT

      OPEN (UNIT=40,FILE='RBF.DAT',STATUS='OLD') 
!     Skip header 
      DO IREAD = 1, 21 
         READ (40,*) 
      END DO 
      IOP = 0 
!         = initialize sequential level index op Opacity Project data 
!     Read number of elements in file 
      READ (40,*) NEOP 
      DO IEOP = 1, NEOP 
!        Skip element name header 
         DO IREAD = 1, 3 
            READ (40,*) 
         END DO 
!        Read number of ionization stages of current element in  file 
         READ (40,*) NIOP 
         DO IIOP = 1, NIOP 
!           Read ion identifyer, atomic & electron number, # of levels 
!           for current ion 
            READ (40,*) IONID, IATOM_OP, IELEC_OP, NLEVEL_OP 
            DO ILOP = 1, NLEVEL_OP 
!              Increase sequential level index of Opacity Project data 
               IOP = IOP+1 
!              Read level identifyer and number of sigma fit points 
               READ (40,*) IDLVOP(IOP), NOP(IOP) 
!              Read normalized log10 frequency and log10 cross section values 
               DO IS = 1, NOP(IOP) 
                  READ (40,*) INDEX, XOP(IS,IOP), SOP(IS,IOP) 
               END DO 
            END DO 
         END DO 
      END DO 
      NTOTOP  = IOP 
!             = total number of levels in Opacity Project data 
      LOPREA  = .TRUE. 
!             = set flag as data has been read in 
! 
      RETURN 
      END SUBROUTINE OPDATA 
! 
! 
! 
!     ****************************************************************** 
! 
! 
      FUNCTION YLINTP(XINT,X,Y,N,NTOT) 
!     ================================= 
! 
!     linear interpolation routine. Determines YINT = Y(XINT) from 
!     grid Y(X) with N points and dimension. 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) ::  X(NTOT),Y(NTOT) 
! 
!     bisection (see Numerical Recipes par 3.4 page 90) 
      JL = 0 
      JU = N+1 
      DO 
         IF(JU-JL.GT.1) THEN 
            JM = (JU+JL)/2 
            IF ((X(N).GT.X(1)).EQV.(XINT.GT.X(JM))) THEN 
               JL = JM 
            ELSE 
               JU = JM 
            END IF 
            CYCLE 
          ELSE 
            EXIT 
         END IF 
      END DO 
      J = JL 
      IF (J.EQ.N) J = J-1 
      IF (J.EQ.0) J = J+1 
      RC         = (Y(J+1)-Y(J))/(X(J+1)-X(J)) 
      YLINTP = RC*(XINT-X(J))+Y(J) 
! 
      RETURN 
      END FUNCTION YLINTP 
! 
! 
!     **************************************************************** 
! 
! 
 
      SUBROUTINE OPAC(ID,ABSO,EMIS,SCAT) 
!     ======================================== 
! 
!     Absorption, emission, and scattering coefficients 
!     at depth ID and for several frequencies (some or all) 
! 
!     Input: ID    - depth index 
!     Output: ABSO - array of absorption coefficient 
!             EMIS - array of emission coefficient 
!             SCAT - array of scattering coefficient (all scattering 
!                    mechanisms except electron scattering) 
! 
! 
      use accura 
      use params 
      use modelp 
      use lindat 
      use synthp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) ::  ABSO(MFREQ),EMIS(MFREQ),SCAT(MFREQ),                   & 
     &              ABLIN(MFREQ),EMLIN(MFREQ) 
      character(len=4) :: typion(30) 
      REAL(DP),PARAMETER:: UN=1.,TEN15=1.E-15,CSB=2.0706E-16,CFF=3.694E8 
! 
      typion(1)=' H- ' 
      typion(2)=' HI ' 
      IF(IMODE.EQ.-1.AND.ID.NE.IDSTD) RETURN 
      T=TEMP(ID) 
      ANE=ELEC(ID) 
      T1=UN/T 
      HKT=HK*T1 
      TK=HKT/H 
      SRT=UN/SQRT(T) 
      SGFF=CFF*SRT 
      CON=CSB*T1*SRT 
      conts=1.e-36/con 
      ABLY=0. 
      EMLY=0. 
      SCLY=0. 
      sce=ane*sige 
      IJ0=2 
      IF(NFREQ.EQ.1) IJ0=1 
      IF(IMODE.EQ.2) IJ0=NFREQ 
      M=3 
      IF(ICONTL.EQ.1) M=1 
      mm=80 
      k=45 
      k=1 
      lpri=.false. 
! 
!     Opacity and emissivity in continuum 
!     **** calculated only in the first and the last frequency ***** 
! 
      CONTFRQ: DO IJ=1,IJ0 
         lpri=id.eq.1.and.ij.eq.1 
         lpri=.false. 
         FR=FREQ(IJ) 
         FR15=FR*TEN15 
         BNU=BN*FR15*FR15*FR15 
         HKF=HKT*FR 
         ABF=0. 
         EBF=0. 
         AFF=0. 
         al=2.997925e18/fr 
         IONS: DO IL=1,NION 
            N0I=NFIRST(IL) 
            N1I=NLAST(IL) 
            NKE=NNEXT(IL) 
            XN=POPUL(NKE,ID) 
! 
!           Bound-free contribution + possibly 
!            pseudo-continuum (accounting for dissolved fraction) 
! 
            BFOPAC: DO II=N0I,N1I 
               SG=0. 
               IF(IFWOP(II).LT.0) THEN 
                  SG=SGMERG(II,ID,FR) 
                ELSE 
                  SG=CROSS(II,IJ) 
                  IF(INDEXP(II).EQ.5) THEN 
                     IZZ=IZ(IEL(II)) 
                     FR0=ENION(II)/6.6256E-27 
                     CALL DWNFR1(FR,FR0,ID,IZZ,DW1) 
                     SG=SG*DW1 
                  END IF 
               END IF 
               if(sg.le.0.) CYCLE BFOPAC 
               ABF=ABF+SG*POPUL(II,ID) 
               XX=SG*XN*EXP(ENION(II)*TK)*WOP(II,ID) 
               IF(XX.lt.conts) cycle bfopac 
               EBF=EBF+XX*CON*G(II)/G(NKE) 
               jj=nke 
               if(lpri.and.mod(id,mm).eq.k)                               & 
     &         write(6,"('bf',i4,i6,f10.3,2i5,2x,a4,2x,2i4,1p5e14.7)")    & 
     &           id,ij,al,ii,jj,typion(il),ii-n0i+1,jj-n0i+1,             & 
     &           popul(ii,id),sg,sg*popul(ii,id),abf,ebf 
            END DO BFOPAC 
 
            IT=IFREE(IL) 
            IF(IT.EQ.0) CYCLE IONS 
! 
!           Free-free contribution 
! 
            IE=IL 
            IF(IE.EQ.IELHM) THEN 
               SFF=SFFHMI(XN,FR,T) 
             ELSE 
               CHA=IZ(IL)*IZ(IL) 
               SF1=CHA*XN*SGFF/(FR*FR*FR) 
               HKFM=HKT*MIN(FF(IL),FR) 
               SF2=EXP(HKFM) 
               IF(IT.EQ.2) THEN 
                  SG=GFREE(T,FR/CHA) 
                  SF2=SF2+SG-UN 
               END IF 
               SFF=SF1*SF2 
            END IF 
            AFF=AFF+SFF 
            if(lpri.and.mod(id,mm).eq.k)                                  & 
     &      write(6,"('ffh',i4,i6,i4,2x,a4,2x,1p4e14.6)")                 & 
     &      id,ij,il,typion(il),sff,aff 
 
         END DO IONS 
! 
!        Additional opacities 
! 
         CALL OPADD(0,ID,FR,ABAD,EMAD,SCAD) 
         if(lpri.and.mod(id,mm).eq.k)                                     & 
     &   write(6,"('ad',i4,i6,1p4e14.6)") id,ij,abad,emad,scad 
         IF(IOPHLI.NE.0) CALL LYMLIN(ID,FR,ABLY,EMLY,SCLY) 
! 
!        Total continuum opacity and emissivity 
! 
         X=EXP(-HKF) 
         X1=UN-X 
         BNE=BNU*X*ANE 
!        ABSO(IJ)=ABF+ANE*(X1*AFF-X*EBF)+ABAD+ABLY 
         ABSO(IJ)=ABF+ANE*(X1*AFF-X*EBF)+ABAD 
!        if(imode.ge.0) abso(ij)=abso(ij)+scad 
         EMIS(IJ)=BNE*(AFF+EBF)+EMAD+EMLY 
         SCAT(IJ)=SCAD+SCLY+sce 
         if(lpri.and.mod(id,mm).eq.k)                                     & 
     &   write(6,"('opac1',i8,i4,1p4e14.6)") ij,id,abso(ij),              & 
     &   emis(ij),scat(ij) 
         IF(IJ.EQ.1) THEN 
            ABLY1=ABLY 
            EMLY1=EMLY 
            SCLY1=SCLY 
         END IF 
      END DO CONTFRQ 
!!    write(*,*) 'after cont opac' 
 
      AVAB=(ABSO(1)+ABSO(2)+SCAT(1)+SCAT(2))*0.5*RELOP 
      IF(NFREQ.LE.2.OR.IMODE.EQ.-1) RETURN 
 
      IF(IMODE.NE.2) THEN 
! 
!        interpolated continuum opacity, emissivity, and scattering 
!        for all frequencies 
! 
         DO IJ=3,NFREQ 
            ABSO(IJ)=FRX1(IJ)*ABSO(2)+FRX2(IJ)*ABSO(1) 
            EMIS(IJ)=FRX1(IJ)*EMIS(2)+FRX2(IJ)*EMIS(1) 
            SCAT(IJ)=FRX1(IJ)*SCAT(2)+FRX2(IJ)*SCAT(1)
         END DO 
! 
!        hydrogen lines -- for IHYL = 0 
!        *** calculated only for the first and the last frequency 
!        and interpolated hydrogen line opacity and emissivity 
!        for all frequencies 
! 
         IF(IHYL.EQ.0) THEN 
            CALL HYDLIN(ID,1,2,ABLIN,EMLIN) 
            DO IJ=M,NFREQ 
               ABSO(IJ)=ABSO(IJ)+FRX1(IJ)*ABLIN(2)+FRX2(IJ)*ABLIN(1) 
               EMIS(IJ)=EMIS(IJ)+FRX1(IJ)*EMLIN(2)+FRX2(IJ)*EMLIN(1) 
            END DO 
         END IF 
! 
!        **** Opacity and emissivity in lines **** 
! 
         CALL LINOP(ID,ABLIN,EMLIN,AVAB) 
         DO IJ=3,NFREQ 
            ABSO(IJ)=ABSO(IJ)+ABLIN(IJ) 
            EMIS(IJ)=EMIS(IJ)+EMLIN(IJ) 
!         if(ij.eq.3) write(*,"('alinop',2i4,1p3e12.3)")                  &
!    &    ij,id,abso(ij),emis(ij),scat(ij)
         END DO 
! 
!        **** Opacity and emissivity in molecular lines **** 
! 
         if(ifmol.gt.0) then 
            do ilist=1,nmlist 
               CALL MOLOP(ID,ABLIN,EMLIN,AVAB,ILIST) 
               DO IJ=3,NFREQ 
                     ABSO(IJ)=ABSO(IJ)+ABLIN(IJ) 
                     EMIS(IJ)=EMIS(IJ)+EMLIN(IJ) 
               END DO 
            end do 
         end if 
      END IF 
! 
!     **** Detailed opacity and emissivity in hydrogen lines **** 
!          (for IHYL=1) 
! 
      IF(IHYL.GT.0.OR.IMODE.EQ.2) THEN 
         CALL HYDLIN(ID,M,NFREQ,ABLIN,EMLIN) 
         DO IJ=M,NFREQ 
            a=abso(ij) 
            e=emis(ij) 
            ABSO(IJ)=ABSO(IJ)+ABLIN(IJ) 
            EMIS(IJ)=EMIS(IJ)+EMLIN(IJ) 
            al=2.997925e18/freq(ij) 
            if(lpri.and.mod(id,mm).eq.k)                                  & 
     &       write(6,"('tot.opac',2i5,f12.3,1p7e12.4)") ij,id,al,         & 
     &       abso(ij),emis(ij),scat(ij),a,ablin(ij),e,emlin(ij) 
         END DO 
      END IF 
! 
!     **** Detailed opacity and emissivity in HE II lines **** 
!          (for IHE2L=1) 
! 
      IF(IHE2L.GT.0) THEN 
         CALL HE2LIN(ID,M,NFREQ,ABLIN,EMLIN) 
         DO IJ=M,NFREQ 
            ABSO(IJ)=ABSO(IJ)+ABLIN(IJ) 
            EMIS(IJ)=EMIS(IJ)+EMLIN(IJ) 
         END DO 
      END IF 
! 
!     opacity due to detailed photoinization cross-section 
!     (from tables; including resonance features) 
!     The two routines may be called and correspond to different formats 
!     as well as difference in INPUT! 
! 
      CALL PHTION(ID,ABSO,EMIS,FREQ,NFREQ) 
      CALL PHTX(ID,ABSO,EMIS,FREQ,0) 
! 
       if(imode.ge.0) then 
          do ij=1,nfreq 
             abso(ij)=abso(ij)+scat(ij) 
          end do 
       end if 
! 
      IF(ICONTL.EQ.1) RETURN 
      ABSO(1)=ABSO(1)-ABLY1 
      EMIS(1)=EMIS(1)-EMLY1 
      SCAT(1)=SCAT(1)-SCLY1 
      ABSO(2)=ABSO(2)-ABLY 
      EMIS(2)=EMIS(2)-EMLY 
      SCAT(2)=SCAT(2)-SCLY 
 
      RETURN 
      END SUBROUTINE OPAC 
! 
! 
!     **************************************************************** 
! 
! 
 
      SUBROUTINE OPACW(ID,ABSO,EMIS,MODC) 
!     =================================== 
! 
!     Absorption, emission, and scattering coefficients 
!     at depth ID and for several frequencies (some or all) 
! 
!     Input: ID    - depth index 
!     Output: ABSO - array of absorption coefficient 
!             EMIS - array of emission coefficient 
!             SCAT - array of scattering coefficient (all scattering 
!                    mechanisms except electron scattering) 
! 
! 
      use accura 
      use params 
      use modelp 
      use lindat 
      use synthp 
      use wincom 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: ABSO(MFREQ),EMIS(MFREQ),SCAT(MFREQ),                    & 
     &            ABLIN(MFREQ),EMLIN(MFREQ),                              & 
     &            ABL1(MFREQC),EML1(MFREQC),SCL1(MFREQC) 
      REAL(DP),PARAMETER:: UN=1.,TEN15=1.E-15,CSB=2.0706E-16,CFF=3.694E8 
! 
      IF(IMODE.EQ.-1.AND.ID.NE.IDSTD) RETURN 
      T=TEMP(ID) 
      ANE=ELEC(ID) 
      T1=UN/T 
      HKT=HK*T1 
      TK=HKT/H 
      SRT=UN/SQRT(T) 
      SGFF=CFF*SRT 
      CON=CSB*T1*SRT 
      conts=1.e-36/con 
      ABLY=0. 
      EMLY=0. 
      SCLY=0. 
      IJ0=2 
      IF(NFREQ.EQ.1) IJ0=1 
      IF(IMODE.EQ.2) IJ0=NFREQ 
      M=3 
! 
!     Opacity and emissivity in continuum 
!     **** calculated only for the continuum frequencies ***** 
! 
      CONTFRQ: DO IJ=1,NFREQC 
         FR=FREQC(IJ) 
         FR15=FR*TEN15 
         BNU=BN*FR15*FR15*FR15 
         HKF=HKT*FR 
         ABF=0. 
         EBF=0. 
         AFF=0. 
         IONS: DO IL=1,NION 
            N0I=NFIRST(IL) 
            N1I=NLAST(IL) 
            NKE=NNEXT(IL) 
            XN=POPUL(NKE,ID) 
! 
!           Bound-free contribution + possibly 
!            pseudo-continuum (accounting for dissolved fraction) 
! 
            BFOPAC: DO II=N0I,N1I 
               SG=0. 
               IF(IFWOP(II).LT.0) THEN 
                  SG=SGMERG(II,ID,FR) 
                ELSE 
                  SG=CROSS(II,IJ) 
                  IF(INDEXP(II).EQ.5) THEN 
                     IZZ=IZ(IEL(II)) 
                     FR0=ENION(II)/6.6256E-27 
                     CALL DWNFR1(FR,FR0,ID,IZZ,DW1) 
                     SG=SG*DW1 
                  END IF 
               END IF 
               ABF=ABF+SG*POPUL(II,ID) 
               XX=SG*XN*EXP(ENION(II)*TK)*WOP(II,ID) 
               IF(XX.lt.conts) cycle bfopac 
               EBF=EBF+XX*CON*G(II)/G(NKE) 
            END DO BFOPAC 
            IT=IFREE(IL) 
            IF(IT.EQ.0) CYCLE IONS 
! 
!           Free-free contribution 
! 
            IE=IL 
            IF(IE.EQ.IELHM) THEN 
               SFF=SFFHMI(XN,FR,T) 
             ELSE 
               CHA=IZ(IL)*IZ(IL) 
               SF1=CHA*XN*SGFF/(FR*FR*FR) 
               HKFM=HKT*MIN(FF(IL),FR) 
               SF2=EXP(HKFM) 
               IF(IT.EQ.2) THEN 
                  SG=GFREE(T,FR/CHA) 
                  SF2=SF2+SG-UN 
               END IF 
               SFF=SF1*SF2 
            END IF 
            AFF=AFF+SFF 
         END DO IONS 
! 
!        Additional opacities 
! 
         CALL OPADD(0,ID,FR,ABAD,EMAD,SCAD) 
         IF(IOPHLI.NE.0) CALL LYMLIN(ID,FR,ABLY,EMLY,SCLY) 
! 
!        Total opacity and emissivity 
! 
         X=EXP(-HKF) 
         X1=UN-X 
         BNE=BNU*X*ANE 
         ABSOC(IJ)=ABF+ANE*(X1*AFF-X*EBF)+ANE*SIGE+ABAD+ABLY 
         if(imode.ge.0) absoc(ij)=absoc(ij)+scad 
         EMISC(IJ)=BNE*(AFF+EBF)+EMAD+EMLY 
         SCATC(IJ)=SCAD+SCLY 
         ABL1(IJ)=ABLY 
         EML1(IJ)=EMLY 
         SCL1(IJ)=SCLY 
      END DO CONTFRQ 
! 
      if(modc.eq.0) return 
! 
      IF(NFREQ.LE.2.OR.IMODE.EQ.-1) RETURN 
! 
!     interpolated continuum and hydrogen line opacity and emissivity 
!     for all frequencies 
! 
      DO IJ=1,NFREQ 
         IJC=IJCINT(IJ) 
         ABSO(IJ)=FRX1(IJ)*ABSOC(IJC)+(1.-FRX1(IJ))*ABSOC(IJC+1) 
         EMIS(IJ)=FRX1(IJ)*EMISC(IJC)+(1.-FRX1(IJ))*EMISC(IJC+1) 
         SCAT(IJ)=FRX1(IJ)*SCATC(IJC)+(1.-FRX1(IJ))*SCATC(IJC+1) 
      END DO 
      IF(IMODE.NE.2) THEN 
! 
!     **** Opacity and emissivity in lines **** 
! 
         CALL LINOPW(ID,ABLIN,EMLIN) 
         DO IJ=1,NFREQ 
            ABSO(IJ)=ABSO(IJ)+ABLIN(IJ) 
            EMIS(IJ)=EMIS(IJ)+EMLIN(IJ) 
         END DO 
! 
!     **** Opacity and emissivity in molecular lines **** 
! 
         if(ifmol.gt.0) then 
            do ilist=1,nmlist 
               CALL MOLOP(ID,ABLIN,EMLIN,AVAB,ILIST) 
               DO IJ=1,NFREQ 
                  ABSO(IJ)=ABSO(IJ)+ABLIN(IJ) 
                  EMIS(IJ)=EMIS(IJ)+EMLIN(IJ) 
               END DO 
            end do 
         end if 
      END IF 
! 
!     **** Detailed opacity and emissivity in hydrogen lines **** 
! 
      CALL HYDLIW(ID,ABLIN,EMLIN) 
      DO IJ=1,NFREQ 
         ABSO(IJ)=ABSO(IJ)+ABLIN(IJ) 
         EMIS(IJ)=EMIS(IJ)+EMLIN(IJ) 
      END DO 
! 
!     **** Detailed opacity and emissivity in HE II lines **** 
!          (for IHE2L=1) 
! 
      CALL HE2LIW(ID,ABLIN,EMLIN) 
      DO IJ=1,NFREQ 
         ABSO(IJ)=ABSO(IJ)+ABLIN(IJ) 
         EMIS(IJ)=EMIS(IJ)+EMLIN(IJ) 
      END DO 
! 
!     opacity due to detailed photoinization cross-section 
!     (from tables; including resonance features) 
!     The two routines may be called and correspond to different formats 
!     as well as difference in INPUT! 
! 
      CALL PHTION(ID,ABSO,EMIS,FREQ,NFREQ) 
      CALL PHTX(ID,ABSO,EMIS,FREQ,0) 
! 
      IF(ICONTL.EQ.1) RETURN 
      DO IJ=1,NFREQC 
         ABSOC(IJ)=ABSOC(IJ)-ABL1(IJ) 
         EMISC(IJ)=EMISC(IJ)-EML1(IJ) 
         SCATC(IJ)=SCATC(IJ)-SCL1(IJ) 
      END DO 
      RETURN 
      END SUBROUTINE OPACW 
! 
! 
! ******************************************************************** 
! 
! 
 
      SUBROUTINE OPACON(ID) 
!     ===================== 
! 
!     Absorption, emission, and scattering coefficients 
!     at depth ID and for several frequencies (some or all) 
! 
!     Input: ID    - depth index 
!     Output: ABSO - array of absorption coefficient 
!             EMIS - array of emission coefficient 
!             SCAT - array of scattering coefficient 
! 
! 
      use accura 
      use params 
      use modelp 
      use lindat 
      use synthp 
      use wincom 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP),PARAMETER:: UN=1.,TEN15=1.E-15,CSB=2.0706E-16,CFF=3.694E8 
! 
      T=TEMP(ID) 
      ANE=ELEC(ID) 
      T1=UN/T 
      HKT=HK*T1 
      TK=HKT/H 
      SRT=UN/SQRT(T) 
      SGFF=CFF*SRT 
      CON=CSB*T1*SRT 
      ABLY=0. 
      EMLY=0. 
      SCLY=0. 
      sce=ane*sige 
! 
!     Opacity and emissivity in continuum 
!     **** calculated only for the continuum frequencies ***** 
! 
      FRQ: DO IJ=1,NFREQC 
         FR=FREQC(IJ) 
         FR15=FR*TEN15 
         BNU=BN*FR15*FR15*FR15 
         HKF=HKT*FR 
         ABF=0. 
         EBF=0. 
         AFF=0. 
         IONS: DO IL=1,NION 
            N0I=NFIRST(IL) 
            N1I=NLAST(IL) 
            NKE=NNEXT(IL) 
            XN=POPUL(NKE,ID) 
! 
!           Bound-free contribution + possibly 
!           pseudo-continuum (accounting for dissolved fraction) 
! 
            BFOPAC: DO II=N0I,N1I 
               SG=0. 
               IF(IFWOP(II).LT.0) THEN 
                  SG=SGMERG(II,ID,FR) 
                ELSE 
                  SG=CROSS(II,IJ) 
                  if(sg.le.0.) cycle bfopac 
                  IF(INDEXP(II).EQ.5) THEN 
                     IZZ=IZ(IEL(II)) 
                     FR0=ENION(II)/6.6256E-27 
                     CALL DWNFR1(FR,FR0,ID,IZZ,DW1) 
                     SG=SG*DW1 
                  END IF 
               END IF 
               if(popul(ii,id).lt.1.e-20.or.xn.lt.1.e-20) cycle bfopac 
               ABF=ABF+SG*POPUL(II,ID) 
               XX=SG*XN*EXP(ENION(II)*TK-hkf)*WOP(II,ID) 
               ee=exp(enion(ii)*tk-hkf) 
               EBF=EBF+XX*CON*G(II)/G(NKE) 
            END DO BFOPAC 
            IT=IFREE(IL) 
            IF(IT.EQ.0) CYCLE IONS 
! 
!           Free-free contribution 
! 
            IE=IL 
            IF(IE.EQ.IELHM) THEN 
               SFF=SFFHMI(XN,FR,T) 
             ELSE 
               CHA=IZ(IL)*IZ(IL) 
               SF1=CHA*XN*SGFF/(FR*FR*FR) 
               HKFM=HKT*MIN(FF(IL),FR) 
               SF2=EXP(HKFM) 
               IF(IT.EQ.2) THEN 
                  SG=GFREE(T,FR/CHA) 
                  SF2=SF2+SG-UN 
               END IF 
               SFF=SF1*SF2 
            END IF 
            AFF=AFF+SFF 
         END DO IONS 
! 
!        Additional opacities 
! 
         CALL OPADD(0,ID,FR,ABAD,EMAD,SCAD) 
         IF(IOPHLI.NE.0) CALL LYMLIN(ID,FR,ABLY,EMLY,SCLY) 
! 
!        Total opacity and emissivity 
! 
         X=EXP(-HKF) 
         X1=UN-X 
         BNE=BNU*X*ANE 
         ABSOC(IJ)=ABF+ANE*(X1*AFF-EBF)+ABAD+ABLY 
         EMISC(IJ)=BNE*AFF+BNU*ANE*EBF+EMAD+EMLY 
         SCATC(IJ)=SCAD+SCLY+sce 
 
      END DO FRQ 
! 
      CALL PHTION(ID,ABSOC,EMISC,FREQC,NFREQC) 
      CALL PHTX(ID,ABSOC,EMISC,FREQC,1) 
! 
      RETURN 
      END SUBROUTINE OPACON 
! 
! 
! ******************************************************************** 
! 
! 
      FUNCTION SGMERG(II,ID,FR) 
!     ========================= 
!     formal routine - taken from TLUSTY, but not used here 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: FRH=3.28805E15,PH2=2.815e29*2.,              & 
     &                       EHB=157802.77355 
! 
      sgmerg=0. 
      if(id.gt.0) return 
      IE=IEL(II) 
      CHA=IZ(IE)*IZ(IE) 
      g(ii)=gmer(imrg(ii),id) 
      T1=1./TEMP(ID) 
      EX=EHB*CHA*T1 
      II0=NQUANT(II-1)+1 
      SUM=0. 
      SUD=0. 
      DO I=II0,NLMX 
         X=I 
         XI=1./(X*X) 
         FREDG=FRH*CHA*XI 
         IF(FR.LT.FREDG) CYCLE 
         EXI=EXP(EX*XI) 
          S=EXI*WNHINT(I,ID)*XI/X 
          SUM=SUM+S 
      END DO 
      SG0=PH2/(FR*FR*FR*G(II))*CHA*CHA 
      SGMERG=SUM*SG0 
      RETURN 
      END FUNCTION SGMERG 
! 
! 
!     **************************************************************** 
! 
      FUNCTION GFREE(T,FR) 
!     ==================== 
! 
!     Hydrogenic free-free Gaunt factor, for temperature T and 
!     frequency FR 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      THET=5040.4/T 
      IF(THET.LT.4.E-2) THET=4.E-2 
      X=FR/2.99793E14 
      IF(X.LE.1.) THEN 
         IF(X.LT.0.2) X=0.2 
         GFREE=(1.0823+2.98E-2/THET)+(6.7E-3+1.12E-2/THET)/X 
         RETURN 
       ELSE 
         C1=(3.9999187E-3-7.8622889E-5/THET)/THET+1.070192 
         C2=(6.4628601E-2-6.1953813E-4/THET)/THET+2.6061249E-1 
         C3=(1.3983474E-5/THET+3.7542343E-2)/THET+5.7917786E-1 
         C4=3.4169006E-1+1.1852264E-2/THET 
         GFREE=((C4/X-C3)/X+C2)/X+C1 
      END IF 
      RETURN 
      END FUNCTION GFREE 
! 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE LYMLIN(ID,FREQ,ABLY,EMLY,SCLY) 
!     ========================================= 
! 
!     OPACITY OF THE LYMAN LINES WINGS (ALPHA - DELTA) 
!     WITH APPROXIMATE PARTIAL REDISTRIBUTION 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: SN(4),SR(4),SS(4),GS(4),FRLY(4),BNLY(4),GA(4) 
      DATA FRLY / 2.4660375E15, 2.9227111E15, 3.0825469E15, 3.156528E15/  & 
     &    ,BNLY / 5.527E-2,     4.090E-2,     2.699E-2,     1.855E-2  /,  & 
     &     SN   / 1.308E5,      5.280E3,      5.847E2,      1.078E2   /,  & 
     &     SR   / 1.218E-16,    9.196E-17,    1.058E-16,    1.296E-16 /,  & 
     &     SS   / 9.478E-3,     1.600E-2,     1.441E-2,     1.547E-2  /,  & 
     &     GS   / 7.237E-8,     5.432E-6,     5.821E-5,     4.027E-4  /,  & 
     &     GA   / 1.000,        1.791,        2.362,        2.801 / 
 
      integer, save :: ifstrk,ifnat,ifres,ifprd,ifsti 
! 
      data icomp/0/ 
      if(iath.le.0) return 
      if(icomp.eq.0) then 
         icomp=1 
         read(4,*,iostat=ios) ifstrk,ifnat,ifres,ifprd,ifsti 
         if(ios.ne.0) then 
            ifstrk=0 
            ifnat=1 
            ifres=1 
            ifprd=0 
            ifsti=0 
            if(iophli.lt.0) then 
               ifstrk=1 
               ifprd=1 
            end if 
         end if 
      end if 
! 
      ABLY=0. 
      EMLY=0. 
      SCLY=0. 
 
      if(freq.gt.3.3e15) return 
 
      P=POPUL(N0HN,ID) 
      T=TEMP(ID) 
      ANE=ELEC(ID) 
      DO I=1,4 
         DFR=ABS(FRLY(I)-FREQ) 
         IF(DFR.LE.5.E11) DFR=1.E12 
         DFR2=DFR*DFR 
         DFRS=SQRT(DFR) 
         COR=(2.*FREQ/(FREQ+FRLY(I)))**2 
         F=1. 
         IF(iabs(IOPHLI).EQ.2) F=FEAUTR(FREQ,ID) 
         STARK=SS(I)*ANE*F/DFR2/DFRS 
         if(ifstrk.eq.0) stark=0. 
         if(ifnat.eq.0) sn(i)=0. 
         if(ifres.eq.0) sr(i)=0. 
         SGLY=SN(I)*(1.+SR(I)*P)*COR/DFR2+STARK 
         sgly=sgly*wnhint(i+1,id) 
         GAMA=1./(GA(I)+GS(I)*ANE*F/DFRS) 
         if(ifprd.eq.0) gama=0. 
         ABLY=ABLY+P*SGLY 
         EMLY=EMLY+POPUL(N0HN+I,ID)*SGLY*BNLY(I)*(1.-GAMA) 
         if(ifsti.ne.0) ably=ably-popul(n0hn+i,id)*sgly/(i+1)/(i+1) 
         SCLY=SCLY+P*SGLY*GAMA 
      END DO 
      RETURN 
      END SUBROUTINE LYMLIN 
! 
! ******************************************************************** 
! 
      FUNCTION FEAUTR(FREQ,ID) 
!     ======================== 
! 
!     LYMAN-ALPHA STARK BROADENING AFTER N.FEAUTRIER 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: DL(20),F05(20),F10(20),F20(20),F40(20),X(4) 
      DATA F05 / 0.0537, 0.0964, 0.1330, 0.3105, 0.4585, 0.6772, 0.8229,  & 
     &           0.8556, 0.9250, 0.9618, 0.9733, 1.1076, 1.0644, 1.0525,  & 
     &           0.8841, 0.8282, 0.7541, 0.7091, 0.7164, 0.7672/ 
      DATA F10 / 0.1986, 0.2764, 0.3959, 0.5740, 0.7385, 0.9448, 1.0292,  & 
     &           1.0317, 0.9947, 0.8679, 0.8648, 0.9815, 1.0660, 1.0793,  & 
     &           1.0699, 1.0357, 0.9245, 0.8603, 0.8195, 0.7928/ 
      DATA F20 / 0.4843, 0.5821, 0.7003, 0.8411, 0.9405, 1.0300, 1.0029,  & 
     &           0.9753, 0.8478, 0.6851, 0.6861, 0.8554, 0.9916, 1.0264,  & 
     &           1.0592, 1.0817, 1.0575, 1.0152, 0.9761, 0.9451/ 
      DATA F40 / 0.7862, 0.8566, 0.9290, 0.9915, 1.0066, 0.9878, 0.8983,  & 
     &           0.8513, 0.6881, 0.5277, 0.5302, 0.6920, 0.8607, 0.9111,  & 
     &           0.9651, 1.0793, 1.1108, 1.1156, 1.1003, 1.0839/ 
      DATA DL / -150., -120., -90., -60., -40., -20., -10., -8., -4.,     & 
     &          -2., 2., 4., 8., 10., 20., 40., 60., 90., 120., 150./ 
      DLAM=2.997925E18/FREQ-1215.685 
      DO I=2,20 
         IF(DLAM.LE.DL(I)) EXIT 
      END DO 
      J=I-1 
      C=DL(J)-DL(I) 
      A=(DLAM-DL(I))/C 
      B=(DL(J)-DLAM)/C 
      X(1)=F05(J)*A+F05(I)*B 
      X(2)=F10(J)*A+F10(I)*B 
      X(3)=F20(J)*A+F20(I)*B 
      X(4)=F40(J)*A+F40(I)*B 
      J=JT(ID) 
      Y=TI0(ID)*X(J)+TI1(ID)*X(J-1)+TI2(ID)*X(J-2) 
      FEAUTR=0.5*(Y+1.) 
      RETURN 
      END FUNCTION FEAUTR 
! 
! ******************************************************************** 
! 
      SUBROUTINE HYLSET 
!     ================= 
! 
!     Initialization procedure for treating the hydrogen line opacity 
! 
      use accura 
      use params 
      use synthp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: ALB(15) 
      DATA ALB /656.28,486.13,434.05,410.17,397.01,                       & 
     &          388.91,383.54,379.79,377.06,375.02,                       & 
     &          373.44,372.19,371.20,370.39,369.72/ 
! 
!     IHYL=-1  -  hydrogen lines are excluded a priori 
! 
      IHYL=-1 
      if(iath.le.0) return 
      IF(FREQ(2).GE.3.28805E15) RETURN 
      AL0=2.997925E17/FREQ(1) 
      AL1=2.997925E17/FREQ(2) 
      IF(AL0.GT.364..AND.AL1.LT.364.6) RETURN 
      IF(AL0.GT.560..AND.AL1.LT.580.) RETURN 
      IF(AL0.GT.720..AND.AL1.LT.820.3) RETURN 
! 
!     otherwise, hydrogen lines are included 
! 
      IHYL=0 
      M20=40 
      IF(AL1.LT.364.6) THEN 
         ILOWH=1 
         FRION=3.28805E15 
         M10=int(SQRT(3.28805E15/ABS(FRION-FREQ(2)))) 
         IF(FRION.GT.FREQ(1)) M20=int(SQRT(3.28805E15/(FRION-FREQ(1)))) 
         IHYL=1 
         IF(AL0.GT.123.) IHYL=0 
         IF(AL0.GT.104..AND.AL1.LT.120.) IHYL=0 
         IF(AL0.GT.98.5.AND.AL1.LT.102.) IHYL=0 
         IF(IMODE.EQ.2.OR.IHYDPR.NE.0.OR.GRAV.GE.6.) IHYL=1 
       ELSE IF(AL1.LT.820.) THEN 
         ILOWH=2 
         if(vaclim.lt.3600.) then 
         FRION=8.2225E14 
         M10=int(SQRT(3.289017E15/ABS(FRION-FREQ(2)))) 
         else 
         FRION=8.22013E14 
         M10=int(SQRT(3.28805E15/ABS(FRION-FREQ(2)))) 
         end if 
         IF(FRION.GT.FREQ(1)) M20=int(SQRT(3.289017E15/(FRION-FREQ(1)))) 
         DO I=1,15 
            AL=ALB(I) 
            IF(AL.LT.AL0-1..OR.AL.GT.AL1+1.) CYCLE 
            IHYL=1 
            EXIT 
         END DO 
         IF(IMODE.EQ.2.OR.IHYDPR.NE.0.OR.GRAV.GE.6.) IHYL=1 
       ELSE 
         ILOWH=3 
         IHYL=1 
      END IF 
! 
      ihyl=1 
! 
      RETURN 
      END SUBROUTINE HYLSET 
! 
! ******************************************************************** 
! 
      SUBROUTINE HYLSEW(IJ) 
!     ===================== 
! 
!     Initialization procedure for treating the hydrogen line opacity 
! 
      use accura 
      use params 
      use synthp 
      use wincom

      implicit real(dp) (a-h,o-z),logical (l) 
 
! 
!     IHYL=-1  -  hydrogen lines are excluded a priori 
! 
      IHYLW(IJ)=0 
      if(iath.le.0) return 
      FR=FREQ(IJ) 
      IF(FR.GE.3.28805E15) RETURN 
      AL0=2.997925E17/FR 
      AL1=AL0 
      IF(grav.lt.6.) then 
         IF(AL0.GT.160..AND.AL1.LT.364.6) RETURN 
         IF(AL0.GT.506..AND.AL1.LT.630.) RETURN 
         IF(AL0.GT.680..AND.AL1.LT.820.3) RETURN 
       else 
         IF(AL0.GT.540..AND.AL1.LT.600.) RETURN 
         IF(AL0.GT.720..AND.AL1.LT.820.3) RETURN 
      end if 
! 
!     otherwise, hydrogen lines are included 
! 
      IHYLW(IJ)=1 
      M20W(IJ)=40 
      IF(AL1.LT.364.6) THEN 
         ILOWHW(IJ)=1 
         FRION=3.28805E15 
       ELSE IF(AL1.LT.820.) THEN 
         ILOWHW(IJ)=2 
         FRION=8.2225E14 
       ELSE IF(AL1.LT.1458.) THEN 
         ILOWHW(IJ)=3 
         FRION=3.6544142E14 
       ELSE IF(AL1.LT.2278.) THEN 
         ILOWHW(IJ)=4 
         FRION=2.0555837E14 
       ELSE IF(AL1.LT.3281.) THEN 
         ILOWHW(IJ)=5 
         FRION=1.315589E14 
       ELSE IF(AL1.LT.4466.) THEN 
         ILOWHW(IJ)=6 
         FRION=9.136394E13 
       ELSE 
         ILOWHW(IJ)=7 
         FRION=6.7120228E13 
      END IF 
      IF(FRION.GT.FR) M10W(IJ)=int(SQRT(3.289017E15/ABS(FRION-FR))) 
!    * I3,'  TO ',I3/) 
      RETURN 
      END SUBROUTINE HYLSEW 
! 
! ******************************************************************** 
! 
      SUBROUTINE HYDLIN(ID,I0,I1,ABSOH,EMISH) 
!     ======================================= 
! 
!     opacity and emissivity of hydrogen lines 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use hydxen

      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: FRH1=3.28805E15,FRH2=FRH1/4.,                & 
     &  UN=1.,SIXTH=1./6.,CPP=4.1412E-16,CPJ=157803.,                     & 
     &  C00=1.25E-9,CDOP=1.284523E12,CID=0.02654,TWO=2.,                  & 
     &  CPJ4=CPJ/4.,AL10=2.3025851,CINV=UN/2.997925E18,                   & 
     &  CID1=0.01497 
      REAL(DP) ::  PJ(40),PRF0(54),                                       & 
     &             ABSO(MFREQ),EMIS(MFREQ),ABSOH(MFREQ),EMISH(MFREQ) 
      REAL(DP) :: WLIR(15) 
      INTEGER  :: irlow(15),irupp(15),inij(18) 
!     LOGICAL  :: LPR,LID 
 
      DATA FRH    /3.289017E15/ 
      data wlir/                                                          & 
     &    123680., 75005., 59066., 51273.,190570.,113060.,                & 
     &     87577., 75061.,277960.,162050.,123840.,105010.,                & 
     &    223340.,168760.,141790./ 
      data irlow/4*6, 4*7, 4*8, 3*9/ 
      data irupp/7,8,9,10,8,9,10,11,9,10,11,12,11,12,13/ 
      data nlinir/15/ 
      data inij/7*0,9,11,12,14,16,19,20,20,22,24,26/ 
      data jend/30/ 
 
      DATA INIT /0/ 
! 
      lid=id.eq.45 
      lid=.false. 
      DO IJ=I0,I1 
         ABSOH(IJ)=0. 
         EMISH(IJ)=0. 
      END DO 
! 
      if(iath.le.0.or.rrr(1,1,1).eq.0.) return 
      izz=1 
! 
      IF(INIT.EQ.0) THEN 
         DO I=1,4 
            DO J=I+1,22 
               CALL STARK0(I,J,IZZ,XK,WL0,FIJ,FIJ0) 
               WLINE(I,J)=WL0 
               OSCH(I,J)=FIJ+FIJ0 
            END DO 
         END DO 
         INIT=1 
      END IF 
      DO IJ=I0,I1 
         ABSO(IJ)=0. 
         EMIS(IJ)=0. 
      END DO 
! 
       if(ilowh.le.0) return 
! 
      T=TEMP(ID) 
      T1=UN/T 
      SQT=SQRT(T) 
      ANE=ELEC(ID) 
      ANES=EXP(SIXTH*LOG(ANE)) 
      TL=LOG10(T) 
      ANEL=LOG10(ANE) 
! 
!     populations of the first 40 levels of hydrogen 
! 
      ANP=POPUL(NKH,ID) 
      PP=CPP*ANE*ANP*T1/SQT 
      NLH=N1H-N0HN+1 
!     if(ifwop(n1h).lt.0) nlh=nlh-1 
      nlh=nlh-1 
      DO IL=1,NLH 
         X=IL*IL 
         PJ(IL)=POPUL(N0HN+IL-1,ID) 
      END DO 
      DO IL=NLH+1,40 
         X=IL*IL 
         PJ(IL)=PP*EXP(CPJ/X*T1)*X*wnhint(il,id) 
      END DO 
      p2=pp*exp(cpj4*t1)*4.*wnhint(2,id) 
! 
!     Frequency- and line-independent parameters for evaluating the 
!     asymptotic Stark profile 
! 
      F00=C00*ANES*ANES*ANES*ANES 
      DOP0=1.E8*SQRT(1.65E8*T+VTURB(ID)) 
! 
! ------------------------------------------------------------------- 
!     overall loop over spectral series (only in the infrared region) 
! ------------------------------------------------------------------- 
! 
      ISERL=ILOWH 
      ISERU=ILOWH 
! 
      if(wlam(i0).gt.14000.) iseru=4 
      if(wlam(i0).gt.22700.) iseru=5 
      if(wlam(i0).gt.32800.) iseru=6 
      if(wlam(i0).gt.44660.) iseru=7 
      if(wlam(i0).gt.60000.) iserl=4 
 
      if(ilowh.eq.2.and.wlam(i0).le.4340..and.grav.ge.8.5) iserl=1 
! 
      if(iserl.eq.3.and.iseru.eq.3.and.nunbal.gt.0) iserl=2 
      DO IJ=I0,I1 
         ABSO(IJ)=0. 
         EMIS(IJ)=0. 
      END DO 
! 
!     ======================== 
!     loop over spectral series 
!     ======================== 
! 
      LOWLEV: DO I=ISERL,ISERU 
! 
!        skip the following calculations if one uses the Gomez tables 
! 
         if(ihgom.gt.0.and.elec(id).gt.hglim) then 
            if(i.ge.1.and.i.le.ihgom) then 
               call ghydop(id,i0,i1,pj,absoh,emish) 
               exit lowlev 
            end if 
         end if 
! 
         II=I*I 
         XII=UN/II 
         POPI=PJ(I) 
         IF(I.EQ.1) FRH=3.28805E15 
! 
!        determination of which hydrogen lines contribute in a current 
!        frequency region 
! 
         M1=M10 
         IF(I.LT.ILOWH) M1=ILOWH-1 
         M2=M1+1 
         M1=M1-1 
         M2=M20+3 
         IF(M1.LT.I+1) M1=I+1 
         if(grav.gt.3.) then 
            m2=m2+5 
            m1=m1-3 
            if(m1.gt.i+6) m1=m1-3 
         end if 
!        new! 
         if(i.ge.3) then 
            m1=i+1 
            m2=i+40 
         end if 
         if(i.ge.4) m2=i+20 
         if(i.ge.6) m2=i+10 
! 
!        loop over lines which contribute at given wavelength region 
! 
         m1=min(m1,40) 
         m2=min(m2,40) 
         m1=max(m1,i+1) 
         m2=max(m2,i+2) 
         UPPLEV: DO J=M1,M2 
            ILINE=0 
            JJ=J*J 
            XJJ=UN/JJ 
            ABTRA=PJ(I)*WNHINT(J,ID) 
            EMTRA=PJ(J)*WNHINT(I,ID)*II*XJJ*EXP(CPJ*(XII-XJJ)*T1) 
            if(i.le.2.and.j.le.i+2) then 
               abtra=pj(i) 
               emtra=pj(j)*wnhint(i,id)/wnhint(j,id)*                     & 
     &            ii*xjj*exp(cpj*(xii-xjj)*t1) 
            end if 
            IF(I.LE.4.AND.J.LE.22) ILINE=ILIN0(I,J) 
! 
!           quasi-molecular opacity for Lyman-alpha and beta satellites 
! 
            lquasi=i.eq.1.and.j.eq.2.and.nunalp.gt.0 
            lquasi=lquasi.or.i.eq.1.and.j.eq.3.and.nunbet.gt.0 
            lquasi=lquasi.or.i.eq.1.and.j.eq.4.and.nungam.gt.0 
            lquasi=lquasi.or.i.eq.2.and.j.eq.3.and.nunbal.gt.0 
            lalhhe=i.eq.1.and.j.eq.2.and.nunhhe.gt.0 
            if(lquasi) then 
               DO IJ=I0,I1 
                  call allard(wlam(ij),popi,anp,sg,i,j) 
                  ABSO(IJ)=ABSO(IJ)+SG*ABTRA 
                  EMIS(IJ)=EMIS(IJ)+SG*EMTRA 
               END DO 
            end if 
            ahe=0. 
            if(iathe.gt.0) ahe=popul(n0a(iathe),id) 
            if(lalhhe.and.ahe.gt.0.) then 
               rel=1./6.2831855 
               do ij=i0,i1 
                  call lyahhe(wlam(ij),ahe,sg0) 
                  sg=sg0*rel 
                  abso(ij)=abso(ij)+sg*abtra 
                  emis(ij)=emis(ij)+sg*emtra 
               end do 
            end if 
! 
!           lines with special Stark broadening tables 
! 
            NOSPECI: IF(IHYDDK.EQ.0.OR.I.NE.1.OR.J.GT.3) THEN 
               SPECPROF: IF(ILINE.GT.0) THEN 
                  FID=CID*OSCH(I,J) 
! 
!                 switch to either original Lemke/Tremblay of Xenomorph 
! 
                  LEMXEN: if(ilxen(i,j).eq.0.or.anel.lt.xnemin) then 
! 
!                    original Lemke/Tremblay 
! 
                     NWL=NWLHYD(ILINE) 
                     DO IWL=1,NWL 
                        PRF0(IWL)=PRFHYD(ILINE,ID,IWL) 
                     END DO 
                     DO IJ=I0,I1 
                        AL=ABS(WLAM(IJ)-WLINE(I,J)) 
                        IF(AL.LT.1.E-4) AL=1.E-4 
                        IF(ILEMKE.EQ.1) AL=AL/F00 
                        AL=LOG10(AL) 
                        DO IWL=1,NWL-1 
                           IW0=IWL 
                           IF(AL.LE.WLHYD(ILINE,IWL+1)) EXIT 
                        END DO 
                        IW1=IW0+1 
                        PRFF=(PRF0(IW0)*(WLHYD(ILINE,IW1)-AL)+PRF0(IW1)*  & 
     &                 (AL-WLHYD(ILINE,IW0)))/                            & 
     &                 (WLHYD(ILINE,IW1)-WLHYD(ILINE,IW0)) 
                       SG=EXP(PRFF*AL10)*FID 
                        sg0=EXP(PRFF*AL10) 
                        IF(ILEMKE.EQ.1) SG=SG*WLINE(I,J)**2*CINV/F00 
                        ABSO(IJ)=ABSO(IJ)+SG*ABTRA 
                        EMIS(IJ)=EMIS(IJ)+SG*EMTRA 
                        alm=2.997925e18/freq(ij) 
                        xb=1.62068e-5 
                        if(lid.and.abs(alm-1106.7).lt.0.1)                & 
     &                  write(6,"('lem',f10.3,2i3,1p7e11.3)")             & 
     &                 alm,i,j,sg,abtra,emtra,                            & 
     &                 sg*abtra,sg*emtra*xb,abso(ij),emis(ij)*xb 
                     END DO 
! 
!                  XENOMORPH data for selected lines 
! 
                   else 
                     ixn=ilxen(i,j) 
                     nwl=nwlxen(ixn) 
                     fr0l=2.997925e18/wline(i,j) 
                     do ij=i0,i1 
                        al=(freq(ij)-fr0l)/f00 
                        if(abs(al).lt.1.e-4) al=1.e-4 
                        all=log10(abs(al)) 
                        do iwl=1,nwl-1 
                           iw0=iwl 
                           if(all.le.alxen(ixn,iwl+1)) exit 
                        end do 
                        iw1=iw0+1 
                        if(al.gt.0.) then 
                           prff=(prfb(ixn,id,iw0)*(alxen(ixn,iw1)-all)+   & 
     &                     prfb(ixn,id,iw1)*(all-alxen(ixn,iw0)))/        & 
     &                     (alxen(ixn,iw1)-alxen(ixn,iw0)) 
                         else 
                           prff=(prfr(ixn,id,iw0)*(alxen(ixn,iw1)-all)+   & 
     &                     prfr(ixn,id,iw1)*(all-alxen(ixn,iw0)))/        & 
     &                     (alxen(ixn,iw1)-alxen(ixn,iw0)) 
                        end if 
                        sg=exp(prff*al10)*fid/f00 
                        ABSO(IJ)=ABSO(IJ)+SG*ABTRA 
                        EMIS(IJ)=EMIS(IJ)+SG*EMTRA 
                     end do 
                  end if LEMXEN 
! 
!                lines without special Stark broadening tables 
! 
                ELSE 
                  CALL STARK0(I,J,izz,XKIJ,WL0,FIJ,FIJ0) 
                  if((wl0.le.wlam(i1).and.1.25*wl0.gt.wlam(i0)).or.      & 
     &            (wl0.ge.wlam(i0).and.0.75*wl0.lt.wlam(i1))) then 
                     FXK=F00*XKIJ 
                     FXK1=UN/FXK 
                     DOP=DOP0/WL0 
                     DBETA=WL0*WL0*CINV*FXK1 
                     BETAD=DOP*DBETA 
                     FID=CID*FIJ*DBETA 
!                    FID0=CID1*FIJ0/DOP 
                     CALL DIVSTR(AD,DIV) 
                     fac=two 
                     if(lquasi) fac=un 
                     DO IJ=I0,I1 
                        fr=freq(ij) 
                        BETA=ABS(WLAM(IJ)-WL0)*FXK1 
                        IF(I.LT.5) THEN 
                           SG=STARKA(BETA,AD,DIV,fac)*FID 
                           if(iophli.eq.2.and.i.eq.1.and.j.eq.2)          & 
     &                        sg=sg*feautr(fr,id) 
                         ELSE 
                           SG=STARKIR(II,JJ,T,ANE,BETA)*FID 
                        END IF 
                        ABSO(IJ)=ABSO(IJ)+SG*ABTRA 
                        EMIS(IJ)=EMIS(IJ)+SG*EMTRA 
                        alm=2.997925e18/freq(ij) 
                       xb=1.62068e-5 
                       if(lid.and.abs(alm-1106.7).lt.0.1)                 & 
     &                 write(6,"('hhl',f10.3,2i3,1p7e11.3)")              & 
     &                 alm,i,j,sg,abtra,emtra,                            & 
     &                 sg*abtra,sg*emtra*xb,abso(ij),emis(ij)*xb 
                     END DO 
                  END IF 
               END IF SPECPROF 
             ELSE 
! 
!           Koester tables for Ly-alpha and beta 
! 
                anh2m=anh2(id) 
                if(j.eq.2) then 
                   do ij=i0,i1 
                      wld=wlam(ij) 
                      call lalpdk(wld,t,ane,anp,pj(1),anh2m,pr) 
                      abso(ij)=abso(ij)+pr*abtra 
                      emis(ij)=emis(ij)+pr*emtra 
                       alm=2.997925e18/freq(ij) 
                       xb=1.62068e-5 
                       sg=pr 
                       if(lid.and.abs(alm-1106.7).lt.0.1)                 & 
     &                 write(6,"('dka',f10.3,2i3,1p7e11.3)")              & 
     &                 alm,i,j,sg,abtra,emtra,                            & 
     &                 sg*abtra,sg*emtra*xb,abso(ij),emis(ij)*xb 
                   end do 
                 else if(j.eq.3) then 
                   do ij=i0,i1 
                      wld=wlam(ij) 
                      call lbetdk(wld,t,ane,anp,pj(1),pr) 
                      abso(ij)=abso(ij)+pr*abtra 
                      emis(ij)=emis(ij)+pr*emtra 
                      alm=2.997925e18/freq(ij) 
                      xb=1.62068e-5 
                      sg=pr 
                      if(lid.and.abs(alm-1106.7).lt.0.1)                  & 
     &                write(6,"('dkb',f10.3,2i3,1p7e11.3)")               & 
     &                alm,i,j,sg,abtra,emtra,                             & 
     &                sg*abtra,sg*emtra*xb,abso(ij),emis(ij)*xb 
                   end do 
                end if 
            END IF NOSPECI 
         END DO UPPLEV 
      END DO LOWLEV 
! 
!     far infrared hydrogen lines 
! 
      IF(WLAM(I1).GT.58000.) THEN 
         LOWL: DO I=8,18 
            II=I*I 
            XII=UN/II 
            JINI=INIJ(I) 
            UPPL: DO J=JINI,JEND 
               JJ=J*J 
               XJJ=UN/JJ 
               CALL STARK0(I,J,izz,XKIJ,WL0,FIJ,FIJ0) 
               if((wl0.le.wlam(i1).and.1.5*wl0.gt.wlam(i0)).or.           & 
     &         (wl0.ge.wlam(i0).and.0.5*wl0.lt.wlam(i1))) then 
                  FXK=F00*XKIJ 
                  FXK1=UN/FXK 
                  DOP=DOP0/WL0 
                  DBETA=WL0*WL0*CINV*FXK1 
                  BETAD=DOP*DBETA 
                  FID=CID*FIJ*DBETA 
                  CALL DIVSTR(AD,DIV) 
                  fac=two 
                  DO IJ=I0,I1 
                     fr=freq(ij) 
                     BETA=ABS(WLAM(IJ)-WL0)*FXK1 
                     SG=STARKIR(II,JJ,T,ANE,BETA)*FID 
                     ABSO(IJ)=ABSO(IJ)+SG*ABTRA 
                     EMIS(IJ)=EMIS(IJ)+SG*EMTRA 
                  END DO 
               END IF 
            END DO UPPL 
         END DO LOWL 
      END IF 
! 
      if(wlam(i1).gt.5.e5) then 
         do ij=i0,i1 
            fr=freq(ij) 
            do ilir=1,nlinir 
               if(wlam(ij).gt.wlir(ilir)*0.95.and.                        & 
     &            wlam(ij).lt.wlir(ilir)*1.05) then 
                  j=irupp(ilir) 
                  JJ=J*J 
                  i=irlow(ilir) 
                  II=I*I 
                  XII=UN/II 
                  XJJ=UN/JJ 
                  ABTRA=PJ(I)*WNHINT(J,ID) 
                  EMTRA=PJ(J)*WNHINT(I,ID)*II*XJJ*EXP(CPJ*(XII-XJJ)*T1) 
                  CALL STARK0(I,J,izz,XKIJ,WL0,FIJ,FIJ0) 
                  FXK=F00*XKIJ 
                  FXK1=UN/FXK 
                  DOP=DOP0/WL0 
                  DBETA=WL0*WL0*CINV*FXK1 
                  BETAD=DOP*DBETA 
                  FID=CID*FIJ*DBETA 
                  CALL DIVSTR(AD,DIV) 
                  fac=two 
                  BETA=ABS(WLAM(IJ)-WL0)*FXK1 
                  SG=STARKA(BETA,AD,DIV,fac)*FID 
                  ABSO(IJ)=ABSO(IJ)+SG*ABTRA 
                  EMIS(IJ)=EMIS(IJ)+SG*EMTRA 
               end if 
            end do 
         end do 
      end if 
! 
!     ---------------------------- 
!     total opacity and emissivity 
!     ---------------------------- 
! 
      DO IJ=I0,I1 
         F=FREQ(IJ) 
         F15=F*1.E-15 
         XKF=EXP(-4.79928e-11*F*T1) 
         XKFB=XKF*1.4743E-2*F15*F15*F15 
         ABSOH(IJ)=ABSO(IJ)-XKF*EMIS(IJ) 
         EMISH(IJ)=XKFB*EMIS(IJ) 
      END DO 
      RETURN 
      END SUBROUTINE HYDLIN 
! 
! 
! ******************************************************************** 
! 
      SUBROUTINE HYDLIW(ID,ABSOH,EMISH) 
!     ================================= 
! 
!     opacity and emissivity of hydrogen lines 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use wincom 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: FRH1=3.28805E15,FRH2=FRH1/4.,                & 
     &  UN=1.,SIXTH=1./6.,CPP=4.1412E-16,CPJ=157803.,                     & 
     &  C00=1.25E-9,CDOP=1.284523E12,CID=0.02654,TWO=2.,                  & 
     &  CPJ4=CPJ/4.,AL10=2.3025851,CINV=UN/2.997925E18,                   & 
     &  CID1=0.01497 
      REAL(DP) ::  PJ(40),PRF0(54),                                       & 
     &             ABSO(MFREQ),EMIS(MFREQ),ABSOH(MFREQ),EMISH(MFREQ) 
 
      DATA FRH    /3.289017E15/ 
      DATA INIT /0/ 
 
      if(iath.le.0) return 
      izz=1 
! 
      IF(INIT.EQ.0) THEN 
         DO I=1,4 
            DO J=I+1,22 
               CALL STARK0(I,J,IZZ,XK,WL0,FIJ,FIJ0) 
               WLINE(I,J)=WL0 
               OSCH(I,J)=FIJ+FIJ0 
            END DO 
         END DO 
         INIT=1 
      END IF 
      DO IJ=1,NFREQ 
         ABSO(IJ)=0. 
         EMIS(IJ)=0. 
         ABSOH(IJ)=0. 
         EMISH(IJ)=0. 
      END DO 
      T=TEMP(ID) 
      T1=UN/T 
      SQT=SQRT(T) 
      ANE=ELEC(ID) 
      ANES=EXP(SIXTH*LOG(ANE)) 
! 
!     populations of the first 40 levels of hydrogen 
! 
      ANP=POPUL(NKH,ID) 
      PP=CPP*ANE*ANP*T1/SQT 
      NLH=N1H-N0HN+1 
      if(ifwop(n1h).lt.0) nlh=nlh-1 
      DO IL=1,40 
         X=IL*IL 
         IF(IL.LE.NLH) PJ(IL)=POPUL(N0HN+IL-1,ID) 
         IF(IL.GT.NLH) PJ(IL)=PP*EXP(CPJ/X*T1)*X*wnhint(il,id) 
      end do 
      p2=pp*exp(cpj4*t1)*4.*wnhint(2,id) 
! 
!     Frequency- and line-independent parameters for evaluating the 
!     asymptotic Stark profile 
! 
      F00=C00*ANES*ANES*ANES*ANES 
      DOP0=1.E8*SQRT(1.65E8*T+VTURB(ID)) 
! 
! ------------------------------------------------------------------- 
!     overall loop over spectral series (only in the infrared region) 
! ------------------------------------------------------------------- 
! 
      FREQLOOP: DO IJ=1,NFREQ 
         IF(IHYLW(IJ).LE.0) CYCLE FREQLOOP 
         ISERL=ILOWHW(IJ) 
         ISERU=ILOWHW(IJ) 
         IF(WLAM(IJ).GT.17000..AND.WLAM(IJ).LE.21000.) THEN 
            ISERL=3 
            ISERU=4 
          ELSE IF(WLAM(IJ).GT.22700..AND.WLAM(IJ).LE.29000.) THEN 
            ISERL=4 
            ISERU=5 
          ELSE IF(WLAM(IJ).GT.32800..AND.WLAM(IJ).LE.37000.) THEN 
            ISERL=5 
            ISERU=6 
          ELSE IF(WLAM(IJ).GT.37000..AND.WLAM(IJ).LE.44600.) THEN 
            ISERL=4 
            ISERU=6 
          ELSE IF(WLAM(IJ).GT.44660..AND.WLAM(IJ).LE.58300.) THEN 
            ISERL=5 
            ISERU=7 
          ELSE IF(WLAM(IJ).GT.58300..AND.WLAM(IJ).LE.72000.) THEN 
            ISERL=6 
            ISERU=8 
          ELSE IF(WLAM(IJ).GT.72000..AND.WLAM(IJ).LE.73800.) THEN 
            ISERL=5 
            ISERU=8 
          ELSE IF(WLAM(IJ).GT.73800..AND.WLAM(IJ).LE.77000.) THEN 
            ISERL=5 
            ISERU=9 
          ELSE IF(WLAM(IJ).GT.77000.) THEN 
            ISERL=6 
            ISERU=9 
         END IF 
! 
         if(iserl.eq.3.and.iseru.eq.3.and.nunbal.gt.0) iserl=2 
! 
         ABSO(IJ)=0. 
         EMIS(IJ)=0. 
 
         LOWLEV: DO I=ISERL,ISERU 
            II=I*I 
            XII=UN/II 
            PLTEI=PP*EXP(CPJ*T1*XII)*II 
            POPI=PJ(I) 
            IF(I.EQ.1) FRH=3.28805E15 
! 
!           determination of which hydrogen lines contribute in a current 
!           frequency region 
! 
            M1=M10W(IJ) 
            IF(I.LT.ILOWHW(IJ)) M1=ILOWHW(IJ)-1 
            M2=M1+1 
            M1=M1-1 
            M2=M20W(IJ)+3 
            IF(M1.LT.I+1) M1=I+1 
            if(grav.gt.3.) then 
               m2=m2+5 
               m1=m1-3 
               if(m1.gt.i+6) m1=m1-3 
            end if 
            if(grav.gt.6.) then 
               m2=m2+2 
               m1=m1-1 
               if(m1.gt.i+6) m1=m1-1 
            end if 
            IF(M1.LT.I+1) M1=I+1 
            IF(M2.GT.40) M2=40 
! 
            A=0. 
            E=0. 
! 
!           loop over lines which contribute at given wavelength region 
! 
            UPPLEV: DO J=M1,M2 
               IF(I.EQ.1.AND.J.LE.5.AND.IOPHLI.LT.0) CYCLE UPPLEV 
               ILINE=0 
               JJ=J*J 
               XJJ=UN/JJ 
               ABTRA=PJ(I)*WNHINT(J,ID) 
               EMTRA=PJ(J)*WNHINT(I,ID)*II*XJJ*EXP(CPJ*(XII-XJJ)*T1) 
               if(i.le.2.and.j.le.i+2) then 
                  abtra=pj(i) 
                  emtra=pj(j)*wnhint(i,id)/wnhint(j,id)*                  & 
     &               ii*xjj*exp(cpj*(xii-xjj)*t1) 
               end if 
               IF(I.LE.4.AND.J.LE.22) ILINE=ILIN0(I,J) 
! 
!           quasi-molecular opacity for Lyman-alpha and beta satellites 
! 
               lquasi=i.eq.1.and.j.eq.2.and.nunalp.gt.0 
               lquasi=lquasi.or.i.eq.1.and.j.eq.3.and.nunbet.gt.0 
               lquasi=lquasi.or.i.eq.1.and.j.eq.4.and.nungam.gt.0 
               lquasi=lquasi.or.i.eq.2.and.j.eq.3.and.nunbal.gt.0 
               if(lquasi) then 
                  CALL STARK0(I,J,izz,XKIJ,WL0,FIJ,FIJ0) 
                  FXK=F00*XKIJ 
                  FXK1=UN/FXK 
                  DOP=DOP0/WL0 
                  DBETA=WL0*WL0*CINV*FXK1 
                  BETAD=DOP*DBETA 
                  FID=CID*FIJ*DBETA 
                  CALL DIVSTR(AD,DIV) 
                  fr=freq(ij) 
                  BETA=ABS(WLAM(IJ)-WL0)*FXK1 
                  call allard(wlam(ij),popi,anp,sg,i,j) 
                  sg=sg+STARKA(BETA,AD,DIV,UN)*FID 
                  ABSO(IJ)=ABSO(IJ)+SG*ABTRA 
                  EMIS(IJ)=EMIS(IJ)+SG*EMTRA 
                  CYCLE UPPLEV 
               end if 
! 
!        lines with special Stark broadening tables 
! 
               IF(ILINE.GT.0) THEN 
                  NWL=NWLHYD(ILINE) 
                  DO IWL=1,NWL 
                     PRF0(IWL)=PRFHYD(ILINE,ID,IWL) 
                  END DO 
                  FID=CID*OSCH(I,J) 
                  AL=ABS(WLAM(IJ)-WLINE(I,J)) 
                  IF(AL.LT.1.E-4) AL=1.E-4 
                  IF(ILEMKE.EQ.1) AL=AL/F00 
                  AL=LOG10(AL) 
                  DO IWL=1,NWL-1 
                     IW0=IWL 
                     IF(AL.LE.WLHYD(ILINE,IWL+1)) EXIT 
                  END DO 
                  IW1=IW0+1 
                  PRFF=(PRF0(IW0)*(WLHYD(ILINE,IW1)-AL)+PRF0(IW1)*        & 
     &              (AL-WLHYD(ILINE,IW0)))/                               & 
     &              (WLHYD(ILINE,IW1)-WLHYD(ILINE,IW0)) 
                  SG=EXP(PRFF*AL10)*FID 
                  IF(ILEMKE.EQ.1) SG=SG*WLINE(I,J)**2*CINV/F00 
                  ABSO(IJ)=ABSO(IJ)+SG*ABTRA 
                  EMIS(IJ)=EMIS(IJ)+SG*EMTRA 
! 
!         lines without special Stark broadening tables 
! 
                ELSE 
                  CALL STARK0(I,J,izz,XKIJ,WL0,FIJ,FIJ0) 
                  FXK=F00*XKIJ 
                  FXK1=UN/FXK 
                  DOP=DOP0/WL0 
                  DBETA=WL0*WL0*CINV*FXK1 
                  BETAD=DOP*DBETA 
                  FID=CID*FIJ*DBETA 
                  CALL DIVSTR(AD,DIV) 
                  fr=freq(ij) 
                  BETA=ABS(WLAM(IJ)-WL0)*FXK1 
                  SG=STARKA(BETA,AD,DIV,TWO)*FID 
                  if(iophli.eq.2.and.i.eq.1.and.j.eq.2)                   & 
     &               sg=sg*feautr(fr,id) 
                  ABSO(IJ)=ABSO(IJ)+SG*ABTRA 
                  EMIS(IJ)=EMIS(IJ)+SG*EMTRA 
               END IF 
            END DO UPPLEV 
         END DO LOWLEV 
! 
!     ---------------------------- 
!     total opacity and emissivity 
!     ---------------------------- 
! 
         F=FREQ(IJ) 
         F15=F*1.E-15 
         XKF=EXP(-4.79928e-11*F*T1) 
         XKFB=XKF*1.4743E-2*F15*F15*F15 
         if(abso(ij).le.0. .and. lasdel) then 
            abso(ij)=0. 
            emis(ij)=0. 
         end if 
         ABSOH(IJ)=ABSO(IJ)-XKF*EMIS(IJ) 
         EMISH(IJ)=XKFB*EMIS(IJ) 
      END DO FREQLOOP 
 
      RETURN 
      END SUBROUTINE HYDLIW 
! 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE HE2SET 
!     ================= 
! 
!     Initialization procedure for treating the He II line opacity 
! 
      use accura 
      use params 
      use synthp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) ::  FRHE(12) 
      DATA FRHE /1.3158153e+16, 3.2895381e+15, 1.4624854e+15,             & 
     &           8.2261878e+14, 5.2647201e+14, 3.6560459e+14,             & 
     &           2.6860713e+14, 2.0565220e+14, 1.6249055e+14,             & 
     &           1.3161730e+14, 1.0877460e+14, 9.1400851e+13/ 
! 
!     IHE2L=-1  -  He II lines are excluded a priori 
! 
      IHE2L=-1 
      IF(IFHE2.LE.0) RETURN 
      IF(FREQ(2).GE.1.315812E16) RETURN 
      AL0=2.997925E17/FREQ(1) 
      AL1=2.997925E17/FREQ(2) 
!      IF(AL0.GT.390.) RETURN 
      if(grav.lt.6.) then 
         IF(AL0.GT.31..AND.AL1.LT.91.1) RETURN 
         IF(AL0.GT.26.1.AND.AL1.LT.29.8) RETURN 
         IF(AL0.GT.24.8.AND.AL1.LT.25.1) RETURN 
         IF(AL0.GT.122.1.AND.AL1.LT.162.9) RETURN 
         IF(AL0.GT.165.1.AND.AL1.LT.204.9) RETURN 
         IF(AL0.GT.109..AND.AL1.LT.120.9) RETURN 
         IF(AL0.GT.103..AND.AL1.LT.107.9) RETURN 
         IF(AL0.GT.99.7.AND.AL1.LT.102.) RETURN 
         IF(AL0.GT.320.8.AND.AL1.LT.364.4) RETURN 
         IF(AL0.GT.273.8.AND.AL1.LT.319.8) RETURN 
         IF(AL0.GT.251.6.AND.AL1.LT.272.8) RETURN 
         IF(AL0.GT.239.0.AND.AL1.LT.250.6) RETURN 
         IF(AL0.GT.231.1.AND.AL1.LT.238.0) RETURN 
         IF(AL0.GT.225.8.AND.AL1.LT.230.1) RETURN 
       else if(grav.lt.7.) then 
         IF(AL0.GT.33..AND.AL1.LT.91.1) RETURN 
         IF(AL0.GT.124.1.AND.AL1.LT.160.9) RETURN 
         IF(AL0.GT.167.1.AND.AL1.LT.202.9) RETURN 
         IF(AL0.GT.111..AND.AL1.LT.118.9) RETURN 
         IF(AL0.GT.322.8.AND.AL1.LT.364.4) RETURN 
         IF(AL0.GT.275.8.AND.AL1.LT.317.8) RETURN 
         IF(AL0.GT.253.6.AND.AL1.LT.270.8) RETURN 
         IF(AL0.GT.241.0.AND.AL1.LT.248.6) RETURN 
         IF(AL0.GT.233.1.AND.AL1.LT.236.0) RETURN 
       else 
         IF(AL0.GT.39..AND.AL1.LT.91.1) RETURN 
         IF(AL0.GT.134.1.AND.AL1.LT.150.9) RETURN 
         IF(AL0.GT.177.1.AND.AL1.LT.202.9) RETURN 
      end if 
! 
!     otherwise, He II lines are included 
! 
      IHE2L=1 
      MHE10=60 
      MHE20=60 
      IF(AL1.LT.91.) THEN 
         ILWHE2=1 
       ELSE IF(AL0.LT.204.) THEN 
         ILWHE2=2 
       ELSE IF(AL0.LT.364.) THEN 
         ILWHE2=3 
       ELSE IF(AL0.LT.569.) THEN 
         ILWHE2=4 
       ELSE IF(AL0.LT.819.) THEN 
         ILWHE2=5 
       ELSE IF(AL0.LT.1116.) THEN 
         ILWHE2=6 
       ELSE IF(AL0.LT.1457.) THEN 
         ILWHE2=7 
       ELSE IF(AL0.LT.1844.) THEN 
         ILWHE2=8 
       ELSE IF(AL0.LT.2277.) THEN 
         ILWHE2=9 
       ELSE IF(AL0.LT.2756.) THEN 
         ILWHE2=10 
       ELSE IF(AL0.LT.3279.) THEN 
         ILWHE2=11 
       ELSE 
         ILWHE2=12 
      END IF 
      FRION=FRHE(ILWHE2) 
      FR1=FRION*ILWHE2*ILWHE2 
      IF(FRION.GT.FREQ(2)) MHE10=int(SQRT(FR1/(FRION-FREQ(2)))) 
      IF(FRION.GT.FREQ(1)) MHE20=int(SQRT(FR1/(FRION-FREQ(1))) ) 
      RETURN 
      END SUBROUTINE HE2SET 
! 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE HE2SEW(IJ) 
!     ===================== 
! 
!     Initialization procedure for treating the He II line opacity 
! 
      use accura 
      use params 
      use synthp 
      use wincom

      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: FRHE(12) 
      DATA FRHE /1.3158153e+16, 3.2895381e+15, 1.4624854e+15,             & 
     &           8.2261878e+14, 5.2647201e+14, 3.6560459e+14,             & 
     &           2.6860713e+14, 2.0565220e+14, 1.6249055e+14,             & 
     &           1.3161730e+14, 1.0877460e+14, 9.1400851e+13/ 
! 
!     IHE2L=-1  -  He II lines are excluded a priori 
! 
      IHE2LW(IJ)=-1 
      IF(IFHE2.LE.0) RETURN 
      FR=FREQ(IJ) 
      AL0=2.997925E17/FR 
      AL1=2.997925E17/FR 
      if(grav.lt.6.) then 
         IF(AL0.GT.31..AND.AL1.LT.91.1) RETURN 
         IF(AL0.GT.26.1.AND.AL1.LT.29.8) RETURN 
         IF(AL0.GT.24.8.AND.AL1.LT.25.1) RETURN 
         IF(AL0.GT.122.1.AND.AL1.LT.162.9) RETURN 
         IF(AL0.GT.165.1.AND.AL1.LT.204.9) RETURN 
         IF(AL0.GT.109..AND.AL1.LT.120.9) RETURN 
         IF(AL0.GT.103..AND.AL1.LT.107.9) RETURN 
         IF(AL0.GT.99.7.AND.AL1.LT.102.) RETURN 
         IF(AL0.GT.320.8.AND.AL1.LT.364.4) RETURN 
         IF(AL0.GT.273.8.AND.AL1.LT.319.8) RETURN 
         IF(AL0.GT.251.6.AND.AL1.LT.272.8) RETURN 
         IF(AL0.GT.239.0.AND.AL1.LT.250.6) RETURN 
         IF(AL0.GT.231.1.AND.AL1.LT.238.0) RETURN 
         IF(AL0.GT.225.8.AND.AL1.LT.230.1) RETURN 
       else if(grav.lt.7.) then 
         IF(AL0.GT.33..AND.AL1.LT.91.1) RETURN 
         IF(AL0.GT.124.1.AND.AL1.LT.160.9) RETURN 
         IF(AL0.GT.167.1.AND.AL1.LT.202.9) RETURN 
         IF(AL0.GT.111..AND.AL1.LT.118.9) RETURN 
         IF(AL0.GT.322.8.AND.AL1.LT.364.4) RETURN 
         IF(AL0.GT.275.8.AND.AL1.LT.317.8) RETURN 
         IF(AL0.GT.253.6.AND.AL1.LT.270.8) RETURN 
         IF(AL0.GT.241.0.AND.AL1.LT.248.6) RETURN 
         IF(AL0.GT.233.1.AND.AL1.LT.236.0) RETURN 
       else 
         IF(AL0.GT.39..AND.AL1.LT.91.1) RETURN 
         IF(AL0.GT.134.1.AND.AL1.LT.150.9) RETURN 
         IF(AL0.GT.177.1.AND.AL1.LT.202.9) RETURN 
      end if 
! 
!     otherwise, He II lines are included 
! 
      IHE2LW(IJ)=1 
      MHE10W(IJ)=60 
      MHE20W(IJ)=60 
      IF(AL1.LT.91.) THEN 
         ILWHEW(IJ)=1 
       ELSE IF(AL0.LT.204.) THEN 
         ILWHEW(IJ)=2 
       ELSE IF(AL0.LT.364.) THEN 
         ILWHEW(IJ)=3 
       ELSE IF(AL0.LT.569.) THEN 
         ILWHEW(IJ)=4 
       ELSE IF(AL0.LT.819.) THEN 
         ILWHEW(IJ)=5 
       ELSE IF(AL0.LT.1116.) THEN 
         ILWHEW(IJ)=6 
       ELSE IF(AL0.LT.1457.) THEN 
         ILWHEW(IJ)=7 
       ELSE IF(AL0.LT.1844.) THEN 
         ILWHEW(IJ)=8 
       ELSE IF(AL0.LT.2277.) THEN 
         ILWHEW(IJ)=9 
       ELSE IF(AL0.LT.2756.) THEN 
         ILWHEW(IJ)=10 
       ELSE IF(AL0.LT.3279.) THEN 
         ILWHEW(IJ)=11 
       ELSE 
         ILWHEW(IJ)=12 
      END IF 
      FRION=FRHE(ILWHEW(IJ)) 
      FR1=FRION*ILWHEW(IJ)*ILWHEW(IJ) 
      IF(FRION.GT.FR) MHE10W(IJ)=int(SQRT(FR1/(FRION-FR))) 
      RETURN 
      END SUBROUTINE HE2SEW 
! 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE HE2LIN(ID,I0,I1,ABSOH,EMISH) 
! 
!     opacity and emissivity of He II lines  (these which are not considered 
!     explicitly) 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use heprf 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP),   PARAMETER :: UN=1.,SIXTH=1./6.,                         & 
     &            CPP=4.1412E-16,CPJ=631479.,                             & 
     &            C00=1.25E-9,CDOP=1.284523E12,CID=0.02654,TWO=2.,        & 
     &            CPJ4=CPJ/4.,AL10=2.3025851,CINV=UN/2.997925E18,         & 
     &            CID1=0.01497 
      REAL(DP) :: PJ(80),FRHE(12),OSCHE2(19),PRF0(36),                    & 
     &            ABSO(MFREQ),EMIS(MFREQ),ABSOH(MFREQ),EMISH(MFREQ) 
 
      DATA FRHE /1.3158153e+16, 3.2895381e+15, 1.4624854e+15,             & 
     &           8.2261878e+14, 5.2647201e+14, 3.6560459e+14,             & 
     &           2.6860713e+14, 2.0565220e+14, 1.6249055e+14,             & 
     &           1.3161730e+14, 1.0877460e+14, 9.1400851e+13/ 
      DATA OSCHE2/6.407E-1, 1.506E-1, 5.584E-2, 2.768E-2,                 & 
     &        1.604E-2, 1.023E-2, 6.980E-3,                               & 
     &        8.421E-1, 3.230E-2, 1.870E-2, 1.196E-2, 8.187E-3,           & 
     &        5.886E-3, 4.393E-3, 3.375E-3, 2.656E-3,                     & 
     &        1.038,    1.793E-1, 6.549E-2/ 
! 
      I=ILWHE2 
      izz=2 
      DO IJ=I0,I1 
         ABSO(IJ)=0. 
         EMIS(IJ)=0. 
         ABSOH(IJ)=0. 
         EMISH(IJ)=0. 
      END DO 
      T=TEMP(ID) 
      T1=UN/T 
      SQT=SQRT(T) 
      ANE=ELEC(ID) 
      ANES=EXP(SIXTH*LOG(ANE)) 
! 
!     He III populations (either LTE or NLTE, depending on input model) 
! 
      IF(IELHE2.GT.0) THEN 
         ANP=POPUL(NNEXT(IELHE2),ID) 
         NLHE2=NLAST(IELHE2)-NFIRST(IELHE2)+1 
       ELSE 
         ANP=RRR(ID,3,2) 
         NLHE2=0 
      END IF 
! 
!     populations of the first 60 levels of He II 
! 
      PP=CPP*ANE*ANP*T1/SQT 
      DO IL=1,60 
         X=IL*IL 
         IIL=NFIRST(IELHE2)+IL-1 
         IF(IL.LE.NLHE2) PJ(IL)=POPUL(IIL,ID) 
         IF(IL.GT.NLHE2) PJ(IL)=PP*EXP(CPJ/X*T1)*X*wnhe2(il,id) 
      END DO 
! 
!     Frequency- and line-independent parameters for evaluating the 
!     asymptotic Stark profile 
! 
      F00=3.906e-11*ANES*ANES*ANES*ANES 
      DOP0=1.E8*SQRT(4.12E7*T+VTURB(ID)) 
! 
! ------------------------------------------------------------------- 
!     overall loop over spectral series (only in the infrared region) 
! ------------------------------------------------------------------- 
! 
      ISERU=ILWHE2 
      IF(ILWHE2.LE.3) THEN 
         ISERL=ILWHE2 
       ELSE IF(ILWHE2.LE.5) THEN 
         ISERL=ILWHE2-1 
       ELSE IF(ILWHE2.LE.7) THEN 
         ISERL=ILWHE2-2 
       ELSE IF(ILWHE2.LE.9) THEN 
         ISERL=ILWHE2-3 
       ELSE 
         ISERL=ILWHE2-4 
      END IF 
! 
      DO IJ=I0,I1 
         ABSO(IJ)=0. 
         EMIS(IJ)=0. 
      END DO 
! 
      LOWLEV: DO I=ISERL,ISERU 
         II=I*I 
         XII=UN/II 
         POPI=PJ(I) 
! 
!        determination of which He II lines contribute in a current 
!        frequency region 
! 
         M1=MHE10 
         IF(I.LT.ILWHE2.AND.FRHE(I).GT.FREQ(2)) THEN 
            M1=int(SQRT(FRHE(I)*II/(FRHE(I)-FREQ(2)))) 
         END IF 
         M2=M1+1 
         IF(M1.LT.I+1) M1=I+1 
         M1=M1-1 
         M2=MHE20+3 
         IF(M2.GT.60) M2=60 
         if(grav.gt.6.) then 
            m2=m2+5 
            m1=m1-3 
            if(m1.gt.i+6) m1=m1-3 
         end if 
         IF(M1.LT.I+1) M1=I+1 
         IF(M2.GT.60) M2=60 
! 
!        loop over lines which contribute at given wavelength region 
! 
         UPPLEV: DO J=M1,M2 
            ILINE=0 
            JJ=J*J 
            XJJ=UN/JJ 
            ABTRA=PJ(I)*WNHE2(J,ID) 
            EMTRA=PJ(J)*WNHE2(I,ID)*II*XJJ*EXP(CPJ*(XII-XJJ)*T1) 
!           IF(I.LE.2) THEN 
!           WLIN=227.838/(XII-1./JJ) 
!         ELSE 
!           WLIN=227.7776/(XII-1./JJ) 
!        END IF 
! 
            WL00=227.9384 
            IF(ILHE2(ILINE).GE.3.AND.VACLIM.LE.2001.) WL00=227.7776 
            WLIN=WL00/(XII-XJJ) 
! 
            IF(I.EQ.2) THEN 
               IF(J.EQ.3.AND.IHE2PR.GT.0) ILINE=1 
             ELSE IF(I.EQ.3) THEN 
               IF(J.EQ.4.AND.IHE2PR.GT.0) ILINE=8 
               IF(J.GT.5.AND.J.LE.10.AND.IHE2PR.GT.0) ILINE=J-3 
             ELSE IF(I.EQ.4) THEN 
               IF(J.LE.7.AND.IHE2PR.GT.0) ILINE=J+12 
               IF(J.GE.8.AND.J.LE.15.AND.IHE2PR.GT.0) ILINE=J+1 
            END IF 
            IF(ILINE.GT.0) THEN 
               NWL=NWLHE2(ILINE) 
               DO IWL=1,NWL 
                  PRF0(IWL)=PRFHE2(ILINE,ID,IWL) 
               END DO 
               FID=CID*OSCHE2(ILINE) 
               DO IJ=I0,I1 
                  AL=ABS(WLAM(IJ)-WLIN) 
                  IF(AL.LT.1.E-4) AL=1.E-4 
                  AL=LOG10(AL) 
                  DO IWL=1,NWL-1 
                     IW0=IWL 
                     IF(AL.LE.WLHE2(ILINE,IWL+1)) EXIT 
                  END DO 
                  IW1=IW0+1 
                  PRFF=(PRF0(IW0)*(WLHE2(ILINE,IW1)-AL)+PRF0(IW1)*        & 
     &             (AL-WLHE2(ILINE,IW0)))/                                & 
     &             (WLHE2(ILINE,IW1)-WLHE2(ILINE,IW0)) 
                  SG=EXP(PRFF*AL10)*FID 
                  ABSO(IJ)=ABSO(IJ)+SG*ABTRA 
                  EMIS(IJ)=EMIS(IJ)+SG*EMTRA 
               END DO 
             ELSE 
               CALL STARK0(I,J,izz,XKIJ,WL0,FIJ,FIJ0) 
               FXK=F00*XKIJ 
               FXK1=UN/FXK 
               DOP=DOP0/WL0 
               DBETA=WL0*WL0*CINV*FXK1 
               BETAD=DOP*DBETA 
               FID=CID*FIJ*DBETA 
!              FID0=CID1*FIJ0/DOP 
               CALL DIVHE2(AD,DIV) 
               DO IJ=I0,I1 
                  BETA=ABS(WLAM(IJ)-WL0)*FXK1 
                  SG=STARKA(BETA,AD,DIV,UN)*FID 
!                 if(fid0.gt.0.) then 
!                    xd=beta/betad 
!                    if(xd.lt.5.) sg=sg+exp(-xd*xd)*fid0 
!                 end if 
                  ABSO(IJ)=ABSO(IJ)+SG*ABTRA 
                  EMIS(IJ)=EMIS(IJ)+SG*EMTRA 
               END DO 
            END IF 
         END DO UPPLEV 
      END DO LOWLEV 
! 
!     ---------------------------- 
!     total opacity and emissivity 
!     ---------------------------- 
! 
      DO IJ=I0,I1 
         F=FREQ(IJ) 
         F15=F*1.E-15 
         XKF=EXP(-4.79928e-11*F*T1) 
         XKFB=XKF*1.4743E-2*F15*F15*F15 
         ABSOH(IJ)=ABSO(IJ)-XKF*EMIS(IJ) 
         EMISH(IJ)=XKFB*EMIS(IJ) 
      END DO 
      RETURN 
      END 
! 
! ******************************************************************** 
! 
      SUBROUTINE HE2LIW(ID,ABSOH,EMISH) 
!     ================================= 
! 
!     opacity and emissivity of He II lines  (these which are not considered 
!     explicitly) 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use wincom 
      use heprf 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP),   PARAMETER :: UN=1.,SIXTH=1./6.,                         & 
     &            CPP=4.1412E-16,CPJ=631479.,                             & 
     &            C00=1.25E-9,CDOP=1.284523E12,CID=0.02654,TWO=2.,        & 
     &            CPJ4=CPJ/4.,AL10=2.3025851,CINV=UN/2.997925E18,         & 
     &            CID1=0.01497 
      REAL(DP) :: PJ(80),FRHE(12),OSCHE2(19),PRF0(36),                    & 
     &            ABSO(MFREQ),EMIS(MFREQ),ABSOH(MFREQ),EMISH(MFREQ) 
 
      DATA FRHE /1.3158153e+16, 3.2895381e+15, 1.4624854e+15,             & 
     &           8.2261878e+14, 5.2647201e+14, 3.6560459e+14,             & 
     &           2.6860713e+14, 2.0565220e+14, 1.6249055e+14,             & 
     &           1.3161730e+14, 1.0877460e+14, 9.1400851e+13/ 
      DATA OSCHE2/6.407E-1, 1.506E-1, 5.584E-2, 2.768E-2,                 & 
     &        1.604E-2, 1.023E-2, 6.980E-3,                               & 
     &        8.421E-1, 3.230E-2, 1.870E-2, 1.196E-2, 8.187E-3,           & 
     &        5.886E-3, 4.393E-3, 3.375E-3, 2.656E-3,                     & 
     &        1.038,    1.793E-1, 6.549E-2/ 
! 
      I=ILWHE2 
      izz=2 
      DO IJ=1,NFREQ 
         ABSO(IJ)=0. 
         EMIS(IJ)=0. 
         ABSOH(IJ)=0. 
         EMISH(IJ)=0. 
      END DO 
      IF(IFHE2.LE.0) RETURN 
      T=TEMP(ID) 
      T1=UN/T 
      SQT=SQRT(T) 
      ANE=ELEC(ID) 
      ANES=EXP(SIXTH*LOG(ANE)) 
! 
!     He III populations (either LTE or NLTE, depending on input model) 
! 
      IF(IELHE2.GT.0) THEN 
         ANP=POPUL(NNEXT(IELHE2),ID) 
         NLHE2=NLAST(IELHE2)-NFIRST(IELHE2)+1 
       ELSE 
         ANP=RRR(ID,3,2) 
         NLHE2=0 
      END IF 
! 
!     populations of the first 60 levels of He II 
! 
      PP=CPP*ANE*ANP*T1/SQT 
      DO IL=1,60 
         X=IL*IL 
         IIL=NFIRST(IELHE2)+IL-1 
         IF(IL.LE.NLHE2) PJ(IL)=POPUL(IIL,ID) 
         IF(IL.GT.NLHE2) PJ(IL)=PP*EXP(CPJ/X*T1)*X*wnhe2(il,id) 
      END DO 
! 
!     Frequency- and line-independent parameters for evaluating the 
!     asymptotic Stark profile 
! 
      F00=3.906e-11*ANES*ANES*ANES*ANES 
      DOP0=1.E8*SQRT(4.12E7*T+VTURB(ID)) 
! 
! ------------------------------------------------------------------- 
!     overall loop over spectral series (only in the infrared region) 
! ------------------------------------------------------------------- 
! 
      FREQLOOP: DO IJ=1,NFREQ 
         ABSO(IJ)=0. 
         EMIS(IJ)=0. 
         IF(IHE2LW(IJ).le.0) CYCLE FREQLOOP 
         I=ILWHEW(IJ) 
         FR=FREQ(IJ) 
         ISERU=ILWHEW(IJ) 
         IF(ILWHEW(IJ).LE.3) THEN 
            ISERL=ILWHEW(IJ) 
          ELSE IF(ILWHEW(IJ).LE.5) THEN 
            ISERL=ILWHEW(IJ)-1 
          ELSE IF(ILWHEW(IJ).LE.7) THEN 
            ISERL=ILWHEW(IJ)-2 
          ELSE IF(ILWHEW(IJ).LE.9) THEN 
            ISERL=ILWHEW(IJ)-3 
          ELSE 
            ISERL=ILWHEW(IJ)-4 
         END IF 
! 
! 
         LOWLEV: DO I=ISERL,ISERU 
            II=I*I 
            XII=UN/II 
            PLTEI=PP*EXP(CPJ*T1*XII)*II 
            POPI=PJ(I) 
! 
!           determination of which He II lines contribute in a current 
!           frequency region 
! 
            M1=MHE10W(IJ) 
            IF(I.LT.ILWHEW(IJ).AND.FRHE(I).GT.FR) THEN 
               M1=int(SQRT(FRHE(I)*II/(FRHE(I)-FR))) 
            END IF 
            M2=M1+1 
            IF(M1.LT.I+1) M1=I+1 
            M1=M1-1 
            M2=MHE20W(IJ)+3 
            IF(M2.GT.60) M2=60 
            if(grav.gt.6.) then 
               m2=m2+5 
               m1=m1-3 
               if(m1.gt.i+6) m1=m1-3 
            end if 
            IF(M1.LT.I+1) M1=I+1 
            IF(M2.GT.60) M2=60 
! 
!           loop over lines which contribute at given wavelength region 
! 
            UPPLEV: DO J=M1,M2 
               ILINE=0 
               JJ=J*J 
               XJJ=UN/JJ 
               ABTRA=PJ(I)*WNHE2(J,ID) 
               EMTRA=PJ(J)*WNHE2(I,ID)*II*XJJ*EXP(CPJ*(XII-XJJ)*T1) 
               IF(I.LE.2) THEN 
                  WLIN=227.838/(XII-1./JJ) 
                ELSE 
                  WLIN=227.7776/(XII-1./JJ) 
               END IF 
               IF(I.EQ.2) THEN 
                  IF(J.EQ.3.AND.IHE2PR.GT.0) ILINE=1 
                ELSE IF(I.EQ.3) THEN 
                  IF(J.EQ.4.AND.IHE2PR.GT.0) ILINE=8 
                  IF(J.GT.5.AND.J.LE.10.AND.IHE2PR.GT.0) ILINE=J-3 
                ELSE IF(I.EQ.4) THEN 
                  IF(J.LE.7.AND.IHE2PR.GT.0) ILINE=J+12 
                  IF(J.GE.8.AND.J.LE.15.AND.IHE2PR.GT.0) ILINE=J+1 
               END IF 
               IF(ILINE.GT.0) THEN 
                  NWL=NWLHE2(ILINE) 
                  DO IWL=1,NWL 
                     PRF0(IWL)=PRFHE2(ILINE,ID,IWL) 
                  END DO 
                  FID=CID*OSCHE2(ILINE) 
                  AL=ABS(WLAM(IJ)-WLIN) 
                  IF(AL.LT.1.E-4) AL=1.E-4 
                  AL=LOG10(AL) 
                  DO IWL=1,NWL-1 
                     IW0=IWL 
                     IF(AL.LE.WLHE2(ILINE,IWL+1)) EXIT 
                  END DO 
                  IW1=IW0+1 
                  PRFF=(PRF0(IW0)*(WLHE2(ILINE,IW1)-AL)+PRF0(IW1)*        & 
     &             (AL-WLHE2(ILINE,IW0)))/                                & 
     &             (WLHE2(ILINE,IW1)-WLHE2(ILINE,IW0)) 
                  SG=EXP(PRFF*AL10)*FID 
                  ABSO(IJ)=ABSO(IJ)+SG*ABTRA 
                  EMIS(IJ)=EMIS(IJ)+SG*EMTRA 
                ELSE 
                  CALL STARK0(I,J,izz,XKIJ,WL0,FIJ,FIJ0) 
                  FXK=F00*XKIJ 
                  FXK1=UN/FXK 
                  DOP=DOP0/WL0 
                  DBETA=WL0*WL0*CINV*FXK1 
                  BETAD=DOP*DBETA 
                  FID=CID*FIJ*DBETA 
                  CALL DIVHE2(AD,DIV) 
                  BETA=ABS(WLAM(IJ)-WL0)*FXK1 
                  SG=STARKA(BETA,AD,DIV,UN)*FID 
                  ABSO(IJ)=ABSO(IJ)+SG*ABTRA 
                  EMIS(IJ)=EMIS(IJ)+SG*EMTRA 
               END IF 
            END DO UPPLEV 
         END DO LOWLEV 
! 
!        ---------------------------- 
!        total opacity and emissivity 
!        ---------------------------- 
! 
         F=FREQ(IJ) 
         F15=F*1.E-15 
         XKF=EXP(-4.79928e-11*F*T1) 
         XKFB=XKF*1.4743E-2*F15*F15*F15 
         ABSOH(IJ)=ABSO(IJ)-XKF*EMIS(IJ) 
         EMISH(IJ)=XKFB*EMIS(IJ) 
      END DO FREQLOOP 
 
      RETURN 
      END SUBROUTINE HE2LIW 
! 
! ******************************************************************** 
! 
      SUBROUTINE STARK0(I,J,IZZ,XKIJ,WL0,FIJ,FIJ0) 
! 
!     Auxiliary procedure for evaluating the approximate Stark profile 
!     of hydrogen lines - sets up necessary frequency independent 
!     parameters 
! 
!     Input:  I     - principal quantum number of the lower level 
!             J     - principal quantum number of the upper level 
!             IZZ   - ionic charge (IZZ=1 for hydrogen, etc.) 
!     Output: XKIJ  - coefficients K(i,j) for the Hotzmark profile; 
!                     exact up to j=6, asymptotic for higher j 
!             WL0   - wavelength of the line i-j 
!             FIJ   - Stark f-value for the line i-j 
!             FIJ0  - f-value for the undisplaced component of the line 
! 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP),    PARAMETER :: RYD1=911.763811,RYD2=911.495745,          & 
     &             CXKIJ=5.5E-5,WI1=911.753578, WI2=227.837832,           & 
     &             UN=1.,TEN=10.,TWEN=20.,HUND=100. 
      REAL(DP) ::  FSTARK(10,4),XKIJT(5,4),FOSC0(10,4),FADD(5,5) 
      DATA XKIJT/3.56E-4,5.23E-4,1.09E-3,1.49E-3,2.25E-3,.0125,.0177,     & 
     & .028,.0348,.0493,.124,.171,.223,.261,.342,.683,.866,1.02,1.19,     & 
     & 1.46/ 
      DATA FSTARK/  .1387,    .0791,   .02126,   .01394,   .00642,        & 
     &           4.814E-3, 2.779E-3, 2.216E-3, 1.443E-3, 1.201E-3,        & 
     &              .3921,    .1193,   .03766,   .02209,   .01139,        & 
     &           8.036E-3, 5.007E-3,  3.85E-3, 2.658E-3, 2.151E-3,        & 
     &              .6103,    .1506,   .04931,   .02768,   .01485,        & 
     &             .01023, 6.588E-3, 4.996E-3, 3.524E-3, 2.838E-3,        & 
     &              .8163,    .1788,   .05985,   .03189,   .01762,        & 
     &             .01196, 7.825E-3, 5.882E-3, 4.233E-3, 3.375E-3/ 
      DATA FOSC0 / 0.27746,  0., 0.00773,  0., 0.00134, 0.,               & 
     &             0.000404, 0., 0.000162, 0.,                            & 
     &             0.24869,  0., 0.00701,  0., 0.00131, 0.,               & 
     &             0.000422, 0., 0.000177, 0.,                            & 
     &             0.23175,  0., 0.00653,  0., 0.00118, 0.,               & 
     &             0.000392, 0., 0.000169, 0.,                            & 
     &             0.22148,  0.0005, 0.00563, 0.0004, 0.00108, 0.,        & 
     &             0.000362, 0., 0.000159, 0./ 
      DATA FADD /  1.231, 0.2069, 7.448E-2, 3.645E-2, 2.104E-2,           & 
     &             1.424, 0.2340, 8.315E-2, 4.038E-2, 2.320E-2,           & 
     &             1.616, 0.2609, 9.163E-2, 4.416E-2, 2.525E-2,           & 
     &             1.807, 0.2876, 1.000E-1, 4.787E-2, 2.724E-2,           & 
     &             1.999, 0.3143, 1.083E-1, 5.152E-2, 2.918E-2/ 
! 
      II=I*I 
      JJ=J*J 
      JMIN=J-I 
      IF(JMIN.LE.5.and.i.le.4) THEN 
         XKIJ=XKIJT(JMIN,I) 
       ELSE 
         XKIJ=CXKIJ*(II*JJ)*(II*JJ)/(JJ-II) 
      END IF 
      IF(I.LE.4) THEN 
         IF(JMIN.LE.10) THEN 
            FIJ=FSTARK(JMIN,I) 
            FIJ0=FOSC0(JMIN,I) 
          ELSE 
            CFIJ=((TWEN*I+HUND)*J/(I+TEN)/(JJ-II)) 
            FIJ=FSTARK(10,I)*CFIJ*CFIJ*CFIJ 
            FIJ0=0. 
         END IF 
       ELSE IF(I.LE.9) THEN 
         IF(JMIN.LE.5) THEN 
            FIJ=FADD(JMIN,I-4) 
            FIJ0=0. 
          ELSE 
            CFIJ=((TEN*I+25.)*J/(I+5.)/(JJ-II)) 
            FIJ=FADD(5,I-4)*CFIJ*CFIJ*CFIJ 
            FIJ0=0. 
         END IF 
       ELSE 
         CFIJ=UN*J/(JJ-II) 
         FIJ=1.96*I*CFIJ*CFIJ*CFIJ 
         FIJ0=0. 
      END IF 
! 
!     wavelength with an explicit correction to the air wavalength 
! 
      w0=wi1 
      if(izz.eq.2) w0=wi2 
      WL0=W0/(UN/II-UN/JJ) 
      IF(WL0.GT.vaclim) THEN 
         ALM=1.E8/(WL0*WL0) 
         XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM) 
         WL0=WL0/(XN1*1.D-6+UN) 
      END IF 
      RETURN 
      END SUBROUTINE STARK0 
! 
! ******************************************************************** 
! 
      FUNCTION STARKA(BETA,A,DIV,FAC) 
! 
!     Approximate expressions for the hydrogen Stark profile 
! 
!     Input: BETA  - delta lambda in beta units, 
!            BETAD - Doppler width in beta units 
!            A     - auxiliary parameter 
!                    A=1.5*LOG(BETAD)-1.671 
!            DIV   - only for A > 1; division point between Doppler 
!                    and asymptotic Stark wing, expressed in units 
!                    of betad. 
!                    DIV = solution of equation 
!                    exp(-(beta/betad)**2)/betad/sqrt(pi)= 
!                     = 1.5*FAC*beta**-5/2 
!                    (ie. the point where Doppler profile is equal to 
!                     the asymptotic Holtsmark) 
!                    In order to save computer time, the division point 
!                    DIV is calculated in advance by routine DIVSTR. 
!            FAC   - factor by which the Holtsmark profile is to be 
!                    multiplied to get total Stark Profile 
!                    FAC should be taken to 2 for hydrogen, (and =1 
!                    for He II) 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP),PARAMETER :: F0=-0.5758228,F1=0.4796232,F2=0.0720948/2.,   & 
     &                      AL=1.26,SD=0.5641895,SLO=-2.5,TRHA=1.5,       & 
     &                      BL1=1.52,BL2=8.325,SAC=0.07966/2. 
 
      XD=BETA/BETAD 
! 
!     for a > 1 Doppler core + asymptotic Holtzmark wing with division 
!               point DIV 
! 
      IF(A.GT.AL) THEN 
         IF(XD.LE.DIV) THEN 
            STARKA=SD*EXP(-XD*XD)/BETAD 
          ELSE 
            STARKA=TRHA*FAC*EXP(SLO*LOG(BETA)) 
         END IF 
       ELSE 
! 
!     empirical formula for a < 1 
! 
         IF(BETA.LE.BL1) THEN 
            STARKA=SAC*FAC 
          ELSE IF(BETA.LT.BL2) THEN 
            XL=LOG(BETA) 
            FL=(F0*XL+F1)*XL 
            STARKA=F2*FAC*EXP(FL) 
          ELSE 
            STARKA=TRHA*FAC*EXP(SLO*LOG(BETA)) 
         END IF 
      END IF 
      RETURN 
      END FUNCTION STARKA 
! 
! ******************************************************************* 
! ******************************************************************* 
! 
      FUNCTION STARKIR(II,JJ,T,ANE,BETA) 
!     ================================== 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: PI=3.14159265,PI2=2.*PI,                     & 
     &           OS0=0.026564,RYD=3.28805E15,                             & 
     &           Y2CON=PI*PI*0.5/OS0/CL 
! 
      DEL=BETA/DBETA 
      HKT=HK/T 
      XII=II 
      XJJ=JJ 
      XX=XII/XJJ 
      DD=2.*XJJ*RYD/DEL 
      Y1=XJJ*DEL*0.5*HKT 
      Y2=Y2CON*DEL**2/ANE 
      QSTAT=1.5+.5*(Y1**2-1.384)/(Y1**2+1.384) 
      QIMPA=0. 
      IF(Y1.LE.8..AND.Y1.LT.Y2) THEN 
         EXY2=0. 
         IF(Y2.LE.8.) EXY2=EXPINT(Y2) 
         QIMPA=1.438*SQRT(Y1*(1.-XX))*(.4*EXP(-Y1)+EXPINT(Y1)-.5*EXY2) 
      END IF 
      IF(BETA.LE.20.) THEN 
         PROF=8./(80.+BETA**3) 
         RATIO=QSTAT+QIMPA 
       ELSE 
         PROF=1.5/BETA/BETA/SQRT(BETA) 
         DIOI=PI2*1.48E-25*DD*ANE*(SQRT(DD)*                              & 
     &        (1.3*QSTAT+.3*QIMPA)-3.9*RYD*HKT) 
         RATIO=QSTAT*MIN(1.+DIOI,1.25)+QIMPA 
      END IF 
      STARKIR=PROF*RATIO 
      RETURN 
      END FUNCTION STARKIR 
 
 
! 
! ******************************************************************* 
! ******************************************************************* 
! 
      SUBROUTINE DIVSTR(A,DIV) 
!     ============================== 
! 
!     Auxiliary procedure for STARKA - determination of the division 
!     point between Doppler and asymptotic Stark profiles 
! 
!     Input:  BETAD - Doppler width in beta units 
!     Output: A     - auxiliary parameter 
!                     A=1.5*LOG(BETAD)-1.671 
!             DIV   - only for A > 1; division point between Doppler 
!                     and asymptotic Stark wing, expressed in units 
!                     of betad. 
!                     DIV = solution of equation 
!                     exp(-(beta/betad)**2)/betad/sqrt(pi)=3*beta**-5/2 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: UN=1.,TWO=2.,UNQ=1.25,UNH=1.5,TWH=2.5,       & 
     &                       CA=1.671,BL=5.821,AL=1.26,CX=0.28,DX=0.0001 
! 
      A=UNH*LOG(BETAD)-CA 
      IF(BETAD.LT.BL) RETURN 
      IF(A.GE.AL) THEN 
         X=SQRT(A)*(UN+UNQ*LOG(A)/(4.*A-5.)) 
      ELSE 
         X=SQRT(CX+A) 
      ENDIF 
      DO I=1,5 
         XN=X*(UN-(X*X-TWH*LOG(X)-A)/(TWO*X*X-TWH)) 
         IF(ABS(XN-X).LE.DX) EXIT 
         X=XN 
      END DO 
      DIV=X 
      RETURN 
      END SUBROUTINE DIVSTR 
! 
! ******************************************************************** 
! 
      SUBROUTINE HYDINI 
! 
!     Initializes necessary arrays for evaluating hydrogen line profiles 
!     from the Lemke, Tremblay-Bergeron, or Schoening-Butler tables 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      INTEGER :: IILW(100),IIUP(100) 
      CHARACTER(LEN=1) ::  CHAR 
 
      DATA INIT /0/ 
! 
      IF(INIT.EQ.0) THEN 
         DO I=1,4 
            DO J=I+1,22 
               CALL STARK0(I,J,IZZ,XK,WL0,FIJ,FIJ0) 
               WLINE(I,J)=WL0 
!              OSCH(I,J)=FIJ+FIJ0
            END DO 
         END DO 
         INIT=1 
      END IF 
      DO I=1,4 
         DO J=1,22 
            ILIN0(I,J)=0 
          write(6,"(3i4,f12.3)") i,j,ilin0(i,j),wline(i,j)
         END DO 
      END DO 
! 
! -------------------------------------------- 
!     Schoening-Butler tables - for IHYDPR < 0 
! -------------------------------------------- 
! 
      IF(IHYDPR.LT.0) THEN 
         IHYDPR=67 
         ILEMKE=0 
         NLINE=12 
! 
         OPEN(UNIT=IHYDPR,FILE='./data/hydprf.dat',STATUS='OLD') 
         write(6,*) ' reading Schoening-Butler tables' 
! 
         DO I=1,12 
            READ(IHYDPR,"(1X)") 
         END DO 
         LINES: DO ILINE=1,NLINE 
! 
!        read the tables, which have to be stored in file 
!        unit IHYDPR (which is the input parameter in the progarm) 
! 
            READ(IHYDPR,"(12X,I1,9X,I1)") I,J 
            IF(ILINE.EQ.12) J=10 
            WL0=WLINE(I,J) 
            ILIN0(I,J)=ILINE 
            READ(IHYDPR,*) CHAR,NWL,(WL(I,ILINE),I=1,NWL) 
            READ(IHYDPR,*) CHAR,NT,(XT(I,ILINE),I=1,NT) 
            READ(IHYDPR,*) CHAR,NE,(XNE(I,ILINE),I=1,NE) 
            READ(IHYDPR,"(1X)") 
            NWLH(ILINE)=NWL 
            NWLHYD(ILINE)=NWL 
            NTH(ILINE)=NT 
            NEH(ILINE)=NE 
! 
            DO I=1,NWL 
               IF(WL(I,ILINE).LT.1.E-4) WL(I,ILINE)=1.E-4 
               WLHYD(ILINE,I)=LOG10(WL(I,ILINE)) 
            END DO 
! 
            DO IE=1,NE 
               DO IT=1,NT 
                  READ(IHYDPR,"(1X)") 
                  READ(IHYDPR,*) (PRF(IWL,IT,IE,ILINE),IWL=1,NWL) 
               END DO 
            END DO 
! 
!        coefficient for the asymptotic profile is determined from 
!        the input data 
! 
            XCLOG=PRF(NWL,1,1,ILINE)+2.5*LOG10(WL(NWL,ILINE))+31.5304-    & 
     &         XNE(1,ILINE)-2.*LOG10(WL0) 
            XKLOG=0.6666667*(XCLOG-0.176) 
            XK=EXP(XKLOG*2.3025851) 
! 
            DO ID=1,ND 
! 
!           temperature is modified in order to account for the 
!           effect of turbulent velocity on the Doppler width 
! 
               T=TEMP(ID)+6.06E-9*VTURB(ID) 
               ANE=ELEC(ID) 
               TL=LOG10(T) 
               ANEL=LOG10(ANE) 
               F00=1.25E-9*ANE**0.666666667 
               FXK=F00*XK 
               DOP=1.E8/WL0*SQRT(1.65E8*T) 
               DBETA=WL0*WL0/2.997925E18/FXK 
               BETAD=DBETA*DOP 
! 
!       interpolation to the actual values of temperature and electron 
!       density. The result is stored at array PRFHYD, having indices 
!       ILINE (line number: 1 for L-alpha,..., 4 for H-delta, etc.); 
!                           5 for H-alpha,..., 8 for H-delta, etc.) 
!       ID - depth index 
!       IWL - wavelength index 
! 
               DO IWL=1,NWL 
                  CALL INTHYD(PROF,TL,ANEL,IWL,ILINE) 
                  PRFHYD(ILINE,ID,IWL)=PROF 
               END DO 
            END DO 
         END DO LINES 
         CLOSE(IHYDPR) 
! 
         IHYDPR=-IHYDPR 
         RETURN 
      END IF 
! 
! --------------------------------- 
!     read Lemke or Tremblay tables 
! --------------------------------- 
! 
      if(ihydpr.lt.20) ihydpr=ihydpr+20 
      if(ihydpr.eq.21) then 
         open(unit=ihydpr,file='./data/lemke.dat',status='old') 
         write(6,"(' -----------'/                                        & 
     &       ' reading Lemke tables; ihydpr =',i3,/                       & 
     &       ' -----------')") ihydpr 
       else if(ihydpr.eq.22) then 
         open(unit=ihydpr,file='./data/tremblay.dat',status='old') 
         write(6,"(' -----------'/                                        & 
     &       ' reading Tremblay tables; ihydpr =',i3,/                    & 
     &       ' -----------')") ihydpr 
      end if 
! 
      ILEMKE=1 
      READ(IHYDPR,*) NTAB 
      write(6,"(' ntab',i4)") ntab 
      DO ITAB=1,NTAB 
         ILINEB=ILINE 
         READ(IHYDPR,*) NLLY 
         DO ILI=1,NLLY 
            ILINE=ILINE+1 
            READ(IHYDPR,*) I,J,ALMIN,ANEMIN,TMIN,DLA,DLE,DLT,             & 
     &                     NWL,NE,NT 
            WL0=WLINE(I,J) 
            ILIN0(I,J)=ILINE 
            NWLH(ILINE)=NWL 
            NWLHYD(ILINE)=NWL 
            NTH(ILINE)=NT 
            NEH(ILINE)=NE 
            iilw(iline)=i 
            iiup(iline)=j 
            DO IWL=1,NWL 
               WL(IWL,ILINE)=ALMIN+(IWL-1)*DLA 
               WLHYD(ILINE,IWL)=WL(IWL,ILINE) 
               WL(IWL,ILINE)=EXP(2.3025851*WL(IWL,ILINE)) 
            END DO 
            DO INE=1,NE 
               XNE(INE,ILINE)=ANEMIN+(INE-1)*DLE 
            END DO 
            DO IT=1,NT 
               XT(IT,ILINE)=TMIN+(IT-1)*DLT 
            END DO 
         END DO 
! 
         DO ILI=1,NLLY 
            ILNE=ILINEB+ILI 
            NWL=NWLH(ILNE) 
            READ(IHYDPR,"(1X)") 
            DO INE=1,NEH(ILNE) 
               DO IT=1,NTH(ILNE) 
                  READ(IHYDPR,*) QLT,(PRF(IWL,IT,INE,ILNE),IWL=1,NWL) 
               END DO 
            END DO 
! 
            i=iilw(ilne) 
            j=iiup(ilne) 
            DO ID=1,ND 
               CALL HYDTAB(I,J,ID) 
            END DO 
         END DO 
      END DO 
      NLIHYD=ILNE 
      CLOSE(IHYDPR) 
! 
!     If required, ainitialze Koester tables fro -alpha and beta 
! 
      IF(IHYDDK.GT.0) CALL DKINI 
! 
      RETURN 
      END SUBROUTINE HYDINI 
! 
! 
! ******************************************************************** 
! 
! 
 
      SUBROUTINE HYDTAB(I,J,ID) 
! 
!     interpolated hydrogen line broadening table for line I->J and 
!     for parameters (TEMP, ELEC) at depth ID 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
! 
      ILINE=ILIN0(I,J) 
      IF(ILINE.EQ.0) RETURN 
      WL0=WLINE(I,J) 
      NWL=NWLH(ILINE) 
! 
!     coefficient for the asymptotic profile is determined from 
!     the input data 
! 
      if(id.eq.1) then 
         XCLOG=PRF(NWL,1,1,ILINE)+2.5*WLHYD(ILINE,NWL)-0.477121 
         XKLOG=0.6666667*XCLOG 
         XK=EXP(XKLOG*2.3025851) 
      end if 
! 
!     temperature is modified in order to account for the 
!     effect of turbulent velocity on the Doppler width 
! 
      T=TEMP(ID)+6.06E-9*VTURB(ID) 
      ANE=ELEC(ID) 
      TL=LOG10(T) 
      ANEL=LOG10(ANE) 
      F00=1.25E-9*ANE**0.666666667 
      FXK=F00*XK 
      DOP=1.E8/WL0*SQRT(1.65E8*T) 
      DBETA=WL0*WL0/2.997925E18/FXK 
      BETAD=DBETA*DOP 
! 
!     interpolation to the actual values of temperature and electron 
!     density. The result is stored at array PRFHYD, having indices 
!       ILINE - line number 
!       ID    - depth index 
!       IWL   - wavelength index 
! 
      DO IWL=1,NWL 
         CALL INTHYD(PROF,TL,ANEL,IWL,ILINE) 
         PRFHYD(ILINE,ID,IWL)=PROF 
      END DO 
! 
      RETURN 
      END SUBROUTINE HYDTAB 
! 
! ******************************************************************** 
! 
      SUBROUTINE INTHYD(W0,X0,Z0,IWL,ILINE) 
! 
!     Interpolation in temperature and electron density from the 
!     hydrogen odening tables to the actual valus of 
!     temperature and electron density 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: TWO=2. 
      REAL(DP)            ::  ZZ(3),XX(3),WX(3),WZ(3) 
! 
      NX=3 
      NZ=3 
      NT=NTH(ILINE) 
      NE=NEH(ILINE) 
      BETA=WL(IWL,ILINE)/FXK 
      IF(ILEMKE.EQ.1) THEN 
         BETA=WL(IWL,ILINE)/XK 
         NX=2 
         NZ=2 
      END IF 
! 
!     for values lower than the lowest grid value of electron density 
!     the profiles are determined by the approximate expression 
!     (see STARKA); not by an extrapolation in the HYD tables which may 
!     be very inaccurate 
! 
      IF(Z0.LT.XNE(1,ILINE)*0.99.OR.Z0.GT.XNE(NE,ILINE)*1.01) THEN 
         CALL DIVSTR(A,DIV) 
         W0=STARKA(BETA,A,DIV,TWO)*DBETA 
         W0=LOG10(W0) 
         RETURN 
      END IF 
! 
!     Otherwise, one interpolates (or extrapolates for higher than the 
!     highes grid value of electron density) in the HYD tables 
! 
      DO IZZ=1,NE-1 
         IPZ=IZZ 
         IF(Z0.LE.XNE(IZZ+1,ILINE)) EXIT 
      END DO 
      N0Z=IPZ-NZ/2+1 
      IF(N0Z.LT.1) N0Z=1 
      IF(N0Z.GT.NE-NZ+1) N0Z=NE-NZ+1 
      N1Z=N0Z+NZ-1 
! 
      ZINT: DO IZZ=N0Z,N1Z 
         I0Z=IZZ-N0Z+1 
         ZZ(I0Z)=XNE(IZZ,ILINE) 
! 
!     Likewise, the approximate expression instead of extrapolation 
!     is used for higher that the highest grid value of temperature, 
!     if the Doppler width expressed in beta units (BETAD) is 
!     sufficiently large (> 10) 
! 
         IF(X0.GT.1.01*XT(NT,ILINE).AND.BETAD.GT.10.) THEN 
            CALL DIVSTR(A,DIV) 
            W0=STARKA(BETA,A,DIV,TWO)*DBETA 
            W0=LOG10(W0) 
            RETURN 
         END IF 
! 
!     Otherwise, normal inter- or extrapolation 
! 
!     Both interpolations (in T as well as in electron density) are 
!     by default the quadratic interpolations in logarithms 
! 
         DO IX=1,NT-1 
            IPX=IX 
            IF(X0.LE.XT(IX+1,ILINE)) EXIT 
         END DO 
         N0X=IPX-NX/2+1 
         IF(N0X.LT.1) N0X=1 
         IF(N0X.GT.NT-NX+1) N0X=NT-NX+1 
         N1X=N0X+NX-1 
         DO IX=N0X,N1X 
            I0=IX-N0X+1 
            XX(I0)=XT(IX,ILINE) 
            WX(I0)=PRF(IWL,IX,IZZ,ILINE) 
         END DO 
         IF(WX(1).LT.-99..OR.WX(2).LT.-99..OR.WX(3).LT.-99.) THEN 
            CALL DIVSTR(A,DIV) 
            W0=STARKA(BETA,A,DIV,TWO)*DBETA 
            W0=LOG10(W0) 
            RETURN 
          ELSE 
            WZ(I0Z)=YINT(XX,WX,X0) 
         END IF 
      END DO ZINT 
      W0=YINT(ZZ,WZ,Z0) 
      RETURN 
      END SUBROUTINE INTHYD 
! 
! ******************************************************************** 
! 
      FUNCTION YINT(XL,YL,XL0) 
! 
!     Quadratic interpolation routine 
! 
!     Input:  XL - array of x 
!             YL - array of f(x) 
!             XL0 - the point x(0) to which one interpolates 
! 
      use accura 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: XL(3),YL(3) 
 
      A0=(XL(2)-XL(1))*(XL(3)-XL(2))*(XL(3)-XL(1)) 
      A1=(XL0-XL(2))*(XL0-XL(3))*(XL(3)-XL(2)) 
      A2=(XL0-XL(1))*(XL(3)-XL0)*(XL(3)-XL(1)) 
      A3=(XL0-XL(1))*(XL0-XL(2))*(XL(2)-XL(1)) 
      YINT=(YL(1)*A1+YL(2)*A2+YL(3)*A3)/A0 
      RETURN 
      END FUNCTION YINT 
! 
! ******************************************************************** 
! 
! 
 
      SUBROUTINE HE1INI 
!     ================= 
! 
!     Initializes necessary arrays for evaluating the He I line 
!     absorption profiles using data calculated by Barnard, Cooper 
!     and Smith JQSRT 14, 1025, 1974 (for 4471) 
!     or Shamey, unpublished PhD thesis, 1969 (for other lines) 
! 
!     This procedure is quite analogous to HYDINI for hydrogen lines 
! 
      use accura 
      use params 
      use modelp 
      use heprf 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      DATA NT /4/ 
! 
      IH=67 
      OPEN(UNIT=IH,FILE='./data/he1prf.dat',STATUS='OLD') 
! 
!        read the Barnard, Cooper, Smith tables for He I 4471 line, 
!        which have to be stored in file unit IH 
! 
      NE=7 
      DO IE=1,NE 
         READ(IH,"(/9X,I2,7X,F10.3,13X,I2,6X,E8.1,7X,I3/)")               & 
     &      IL,WL0,IE1,XXNE,NWL 
         NWLAM(IE,1)=NWL 
         XNE447(IE)=LOG10(XXNE) 
         DO I=1,NWL 
            READ(IH,"(5E10.2)") DLM447(I,IE),                             & 
     &                (PRF447(I,IT,IE),IT=1,NT) 
         END DO 
      END DO 
! 
!     read Shamey's tables for He I 4387, 4026, and 4922 lines 
!     which have to be stored in file unit IH 
! 
      NE=8 
      DO ILN=1,3 
         DO IE=1,NE 
            READ(IH,"(/9X,I2,7X,F10.3,13X,I2,6X,E8.1,7X,I3/)")            & 
     &      IL,WL0,IE1,XXNE,NWL 
            NWLAM(IE,ILN+1)=NWL 
            XNEHE1(IE)=LOG10(XXNE) 
            DO I=1,NWL 
               READ(IH,*) DLMHE1(I,IE,ILN),                               & 
     &                    (PRFHE1(I,IT,IE,ILN),IT=1,NT) 
            END DO 
         END DO 
      END DO 
      CLOSE(IH) 
! 
      RETURN 
      END SUBROUTINE HE1INI 
! 
! ******************************************************************** 
! 
 
      FUNCTION WTOT(T,ANE,ID,ILINE) 
!     ============================= 
! 
!     Evaluates the total (electron + ion) impact Stark width 
!     for four HeI lines 
!     After Griem (1974); and Barnard, Cooper, Smith (1974) JQSRT 14, 
!     1025 for the 4471 line 
! 
!     Input: T     - temperature 
!            ANE   - electron density 
!            ID    - depth index 
!            ILINE - index of the line ( = 1  for 4471, 
!                                        = 2  for 4387, 
!                                        = 3  for 4026, 
!                                        = 4  for 4922) 
!     Output: WTOT - Stark width in Angstroms 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) ::  ALPH0(4,4),W0(4,4),ALM0(4) 
      DATA ALPH0 / 0.107, 0.119, 0.134, 0.154,                            & 
     &             0.206, 0.235, 0.272, 0.317,                            & 
     &             0.172, 0.193, 0.218, 0.249,                            & 
     &             0.121, 0.136, 0.157, 0.184/ 
      DATA W0    / 1.460, 1.269, 1.079, 0.898,                            & 
     &             6.130, 5.150, 4.240, 3.450,                            & 
     &             4.040, 3.490, 2.960, 2.470,                            & 
     &             2.312, 1.963, 1.624, 1.315/ 
      DATA ALM0  / 4471.50, 4387.93, 4026.20, 4921.93/ 
      DATA INITI/1/ 
! 
!     change central wavelengths to vacuum, if required 
! 
      if(initi.eq.1.and.vaclim.gt.2000.) then 
         do i=1,4 
            wl0=alm0(i) 
            ALM=1.E8/(WL0*WL0) 
            XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM) 
            WL0=WL0*(XN1*1.D-6+1.) 
            alm0(i)=wl0 
!           write(*,*) 'wtot',i,alm0(i) 
         end do 
         initi=0 
      END IF 
! 
      I=JT(ID) 
      ALPHA=(TI0(ID)*ALPH0(I,ILINE)+TI1(ID)*ALPH0(I-1,ILINE)+             & 
     &      TI2(ID)*ALPH0(I-2,ILINE))*(ANE*1.E-13)**0.25 
      WE=   (TI0(ID)*W0(I,ILINE)+TI1(ID)*W0(I-1,ILINE)+                   & 
     &      TI2(ID)*W0(I-2,ILINE))*ANE*1.E-16 
      F0=1.884E19/ALM0(ILINE)/ALM0(ILINE) 
      SIG=(4.32E-5*WE/SQRT(T)*F0/ANE**0.3333)**0.3333 
      WTOT=WE*(1.+1.36/SIG*ALPHA**0.8889) 
      RETURN 
      END FUNCTION WTOT 
 
! 
! ******************************************************************** 
! 
 
      FUNCTION EXTPRF(DLAM,IT,ILINE,ANEL,DLAST,PLAST) 
!     =============================================== 
! 
!     Extrapolation in wavelengths in Shamey, or Barnard, Cooper, 
!     Smith tables 
!     Special formula suggested by Cooper 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) ::  W0(4,4) 
      DATA W0    / 1.460, 1.269, 1.079, 0.898,                            & 
     &             6.130, 5.150, 4.240, 3.450,                            & 
     &             4.040, 3.490, 2.960, 2.470,                            & 
     &             2.312, 1.963, 1.624, 1.315/ 
! 
      WE=W0(IT,ILINE)*EXP(ANEL*2.3025851)*1.E-16 
      DLASTA=ABS(DLAST) 
      D52=DLASTA*DLASTA*SQRT(DLASTA) 
      F=D52*(PLAST-WE/3.14159/DLAST/DLAST) 
      EXTPRF=(WE/3.14159+F/SQRT(ABS(DLAM)))/DLAM/DLAM 
      RETURN 
      END FUNCTION EXTPRF 
 
! 
! ******************************************************************** 
! 
 
      FUNCTION PHE1(ID,FREQ,ILINE) 
!     ============================ 
! 
!     Absorption profile for four lines of He I, given by 
!     Barnard, Cooper, Smith (1974) JQSRT 14, 1025 for the 4471 line; 
!     Shamey (1969) PhD thesis, for other lines 
! 
!     Input: ID    - depth index 
!            FREQ  - frequency 
!            ILINE - index of the line ( = 1  for 4471, 
!                                        = 2  for 4387, 
!                                        = 3  for 4026, 
!                                        = 4  for 4922) 
! 
!     Output: PHE1 - profile coefficient in frequency units, 
!                    normalized to sqrt(pi) [not unity] 
! 
      use accura 
      use params 
      use modelp 
      use heprf 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      INTEGER, PARAMETER :: NT=4 
      REAL(DP) ::  WLAM0(4),XT0(NT),XX(3),WX(3),YY(2),PP(2),ZZ(3),WZ(3) 
      DATA WLAM0 / 4471.50, 4387.93, 4026.20, 4921.93/ 
      DATA XT0/ 3.699, 4.000, 4.301, 4.602/ 
      DATA INITI/1/ 
 
      SAVE WLAM0 
! 
!     change central wavelengths to vacuum, if required 
! 
      if(initi.eq.1.and.vaclim.gt.2000.) then 
         do i=1,4 
            wl0=wlam0(i) 
            ALM=1.E8/(WL0*WL0) 
            XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM) 
            WL0=WL0*(XN1*1.e-6+1.) 
            wlam0(i)=wl0 
!           write(*,*) 'phe21',i,wlam0(i) 
         end do 
         initi=0 
      END IF 
! 
!     temperature is modified in order to account for the 
!     effect of turbulent velocity on the Doppler width 
! 
      T=TEMP(ID)+2.42E-8*VTURB(ID) 
      TL=LOG10(T) 
      ANE=ELEC(ID) 
      ANEL=LOG10(ANE) 
      ALAM=2.997925E18/FREQ 
      DLAM=ALAM-WLAM0(ILINE) 
      DOPL=SQRT(4.125E7*T)*WLAM0(ILINE)/2.997925E10 
! 
!     isolated line approximation for low electron densities 
! 
      IF(TL.GT.XT0(NT)+0.1) THEN 
         A=WTOT(T,ANE,ID,ILINE)/DOPL 
         V=ABS(DLAM)/DOPL 
         V1=ABS(ALAM-4471.682)/DOPL 
         PHE1=VOIGTK(A,V) 
         IF(ILINE.EQ.1) PHE1=(8.*PHE1+VOIGTK(A,V1))/9. 
         RETURN 
      END IF 
! 
!     otherwise, interpolation (or extrapolation) in tables 
! 
      NX=3 
      NZ=3 
      NY=2 
      NE=8 
      ILNE=ILINE-1 
      IF(ILINE.EQ.1) NE=7 
! 
!     Interpolation in electron density 
! 
      DO JZ=1,NE-1 
         IPZ=JZ 
         IF(ILINE.EQ.1.AND.ANEL.LE.XNE447(JZ+1)) EXIT 
         IF(ILINE.NE.1.AND.ANEL.LE.XNEHE1(JZ+1)) EXIT 
      END DO 
      N0Z=IPZ-NZ/2+1 
      IF(N0Z.LT.1) N0Z=1 
      IF(N0Z.GT.NE-NZ+1) N0Z=NE-NZ+1 
      N1Z=N0Z+NZ-1 
      ZINT: DO JZ=N0Z,N1Z 
         I0Z=JZ-N0Z+1 
         IF(ILINE.EQ.1) ZZ(I0Z)=XNE447(JZ) 
         IF(ILINE.NE.1) ZZ(I0Z)=XNEHE1(JZ) 
! 
!        Interpolation in temperature 
! 
         DO IX=1,NT-1 
            IPX=IX 
            IF(TL.LE.XT0(IX+1)) EXIT 
         END DO 
         N0X=IPX-NX/2+1 
         IF(N0X.LT.1) N0X=1 
         IF(N0X.GT.NT-NX+1) N0X=NT-NX+1 
         N1X=N0X+NX-1 
         XINT: DO IX=N0X,N1X 
            I0X=IX-N0X+1 
            XX(I0X)=XT0(IX) 
! 
!           Interpolation in wavelength 
! 
!           1. For delta lambda beyond tabulated values - special 
!              extrapolation (Cooper's suggestion) 
! 
            NLST=NWLAM(JZ,ILINE) 
            IF(ILINE.EQ.1) THEN 
               D1=DLM447(1,JZ) 
               D2=DLM447(NLST,JZ) 
               IF(DLAM.LT.D1) THEN 
                  PRF0=EXTPRF(DLAM,IX,ILINE,ZZ(I0Z),D1,PRF447(1,IX,JZ)) 
                  WX(I0X)=PRF0 
                  CYCLE ZINT 
                ELSE IF(DLAM.GT.D2) THEN 
                  PRF0=EXTPRF(DLAM,IX,ILINE,ZZ(I0Z),D2,                   & 
     &                PRF447(NLST,IX,JZ)) 
                  WX(I0X)=PRF0 
                  CYCLE XINT 
               END IF 
             ELSE 
               D1=DLMHE1(1,JZ,ILNE) 
               D2=DLMHE1(NLST,JZ,ILNE) 
               IF(DLAM.LT.D1) THEN 
                  PRF0=EXTPRF(DLAM,IX,ILINE,ZZ(I0Z),D1,                   & 
     &                PRFHE1(1,IX,JZ,ILNE)) 
                  WX(I0X)=PRF0 
                  CYCLE XINT 
                ELSE IF(DLAM.GT.D2) THEN 
                  PRF0=EXTPRF(DLAM,IX,ILINE,ZZ(I0Z),D2,                   & 
     &                PRFHE1(NLST,IX,JZ,ILNE)) 
                  WX(I0X)=PRF0 
                  CYCLE XINT 
               END IF 
            END IF 
! 
!           normal linear interpolation in wavelength 
!           (for 4471, linear interpolation in logarithms) 
! 
            DO IY=1,NLST-1 
               IPY=IY 
               IF(ILINE.EQ.1.AND.DLAM.LE.DLM447(IY+1,JZ)) EXIT 
               IF(ILINE.NE.1.AND.DLAM.LE.DLMHE1(IY+1,JZ,ILNE)) EXIT 
            END DO 
            N0Y=IPY-NY/2+1 
            IF(N0Y.LT.1) N0Y=1 
            IF(N0Y.GT.NLST-NY+1) N0Y=NLST-NY+1 
            N1Y=N0Y+NY-1 
            DO IY=N0Y,N1Y 
               I0=IY-N0Y+1 
               IF(ILINE.EQ.1) YY(I0)=DLM447(IY,JZ) 
               IF(ILINE.EQ.1) PP(I0)=LOG(PRF447(IY,IX,JZ)) 
               IF(ILINE.NE.1) YY(I0)=DLMHE1(IY,JZ,ILNE) 
               IF(ILINE.NE.1) PP(I0)=PRFHE1(IY,IX,JZ,ILNE) 
           END DO 
           IF(ILINE.NE.1) THEN 
              WX(I0X)=(PP(2)*(DLAM-YY(1))+PP(1)*(YY(2)-DLAM))/            & 
     &                (YY(2)-YY(1)) 
            ELSE 
             WX(I0X)=(PP(2)*(DLAM-YY(1))+PP(1)*(YY(2)-DLAM))/             & 
     &                (YY(2)-YY(1)) 
             WX(I0X)=EXP(WX(I0X)) 
           END IF 
        END DO XINT 
        WZ(I0Z)=YINT(XX,WX,TL) 
      END DO ZINT 
      W0=YINT(ZZ,WZ,ANEL) 
      PHE1=W0*DOPL*1.772454 
      RETURN 
      END FUNCTION PHE1 
 
! 
! ******************************************************************** 
! 
 
      SUBROUTINE HE2INI 
!     ================= 
! 
!     Initializes necessary arrays for evaluating the He II line 
!     absorption profiles using data calculated by Schoening and 
!     Butler 
! 
!     This procedure is quite analogous to HYDINI for hydrogen lines 
! 
      use accura 
      use params 
      use modelp 
      use heprf 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      DATA NLINE1 /19/ 
! 
      IH=67 
      OPEN(UNIT=IH,FILE='./data/he2prf.dat',STATUS='OLD') 
! 
      DO ILINE=1,NLINE1 
! 
!     read the Schoening and Butler tables, which have to be stored 
!     in file he23prf.dat 
! 
         READ(IH,"(//14X,I2,9X,I2/)") ILHE2(ILINE),IUHE2(ILINE) 
!        IF(ILHE2(ILINE).LE.2) THEN 
!           WL00=227.838 
!         ELSE 
!           WL00=227.7776 
!        END IF 
! 
         WL00=227.9384 
         IF(ILHE2(ILINE).GE.3.AND.VACLIM.LE.2001.) WL00=227.7776 
         WL0=WL00/(1./ILHE2(ILINE)**2-1./IUHE2(ILINE)**2) 
! 
         READ(IH,*) NWL2,(WL2(I,ILINE),I=1,NWL2) 
         READ(IH,"(2X,I4,F10.3,5F12.3)") NT2,(XT2(I),I=1,NT2) 
         READ(IH,"(2X,I4,F10.2,5F12.2/4X,5F12.2)")                        & 
     &      NE2,(XNE2(I,ILINE),I=1,NE2) 
         READ(IH,"(1X)") 
         NWLHE2(ILINE)=NWL2 
! 
         DO I=1,NWL2 
            IF(WL2(I,ILINE).LT.1.E-4) WL2(I,ILINE)=1.E-4 
            WLHE2(ILINE,I)=LOG10(WL2(I,ILINE)) 
         END DO 
! 
         DO IE=1,NE2 
            DO IT=1,NT2 
               READ(IH,"(1X)") 
               READ(IH,"(10F8.3)") (PRF2(IWL,IT,IE),IWL=1,NWL2) 
            END DO 
         END DO 
! 
!        coefficient for the asymptotic profile is determined from 
!        the input data 
! 
         XCLOG=PRF2(NWL2,1,1)+2.5*LOG10(WL2(NWL2,ILINE))+31.831-          & 
     &         XNE2(1,ILINE)-2.*LOG10(WL0) 
         XKLOG=0.6666667*(XCLOG-0.176) 
         XK=EXP(XKLOG*2.3025851) 
         DO ID=1,ND 
            T=TEMP(ID)+2.42E-8*VTURB(ID) 
            ANE=ELEC(ID) 
            TL=LOG10(T) 
            ANEL=LOG10(ANE) 
            F00=1.25E-9*ANE**0.666666667 
            FXK=F00*XK 
            DOP=1.E8/WL0*SQRT(4.12E7*T) 
            DBETA=WL0*WL0/2.997925E18/FXK 
            BETAD=DBETA*DOP 
! 
!     interpolation to the actual values of temperature and electron 
!     density. The result is stored at array PRFHE2, which has indices 
!     ILINE  - index of line 
!     ID     - depth index 
!     IWL    - wavelength index (notice that the wavelength grid may 
!              generally be different for different lines 
! 
            DO IWL=1,NWL2 
               CALL INTHE2(PROF,TL,ANEL,IWL,ILINE) 
               PRFHE2(ILINE,ID,IWL)=PROF 
           END DO 
        END DO 
      END DO 
      CLOSE(IH) 
! 
      RETURN 
      END SUBROUTINE HE2INI 
! 
! ******************************************************************** 
! 
! 
 
      SUBROUTINE INTHE2(W0,X0,Z0,IWL,ILINE) 
!     ===================================== 
! 
!     Interpolation in temperature and electron density from the 
!     Schoening and Butler tables for He II lines to the actual 
!     actual values of temperature and electron density 
! 
!     This procedure is quite analogous to INTHYD for hydrogen lines 
! 
      use accura 
      use params 
      use heprf 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: UN=1. 
      REAL(DP)            :: ZZ(3),XX(3),WX(3),WZ(3) 
! 
      NX=3 
      NZ=3 
! 
!     for values lower than the lowest grid value of electron density 
!     the profiles are determined by the approximate expression 
!     (see STARKA); not by an extrapolation in the tables which may 
!     be very inaccurate 
! 
      IF(Z0.LT.XNE2(1,ILINE)*0.99.OR.Z0.GT.XNE2(NE2,ILINE)*1.01) THEN 
         CALL DIVHE2(A,DIV) 
         W0=STARKA(WL2(IWL,ILINE)/FXK,A,DIV,UN)*DBETA 
         W0=LOG10(W0) 
         RETURN 
      END IF 
! 
!     Otherwise, one interpolates (or extrapolates for higher than the 
!     highes grid value of electron density) in the Schoening and 
!     Butler tables 
! 
      DO IZZ=1,NE2-1 
         IPZ=IZZ 
         IF(Z0.LE.XNE2(IZZ+1,ILINE)) EXIT 
      END DO 
      N0Z=IPZ-NZ/2+1 
      IF(N0Z.LT.1) N0Z=1 
      IF(N0Z.GT.NE2-NZ+1) N0Z=NE2-NZ+1 
      N1Z=N0Z+NZ-1 
! 
      ZINT: DO IZZ=N0Z,N1Z 
         I0Z=IZZ-N0Z+1 
         ZZ(I0Z)=XNE2(IZZ,iline) 
! 
!     Likewise, the approximate expression instead of extrapolation 
!     is used for higher that the highest grid value of temperature, 
!     if the Doppler width expressed in beta units (BETAD) is 
!     sufficiently large (> 10) 
! 
         IF(X0.GT.1.01*XT2(NT2).AND.BETAD.GT.10.) THEN 
            W0=STARKA(WL2(IWL,ILINE)/FXK,A,DIV,UN)*DBETA 
            W0=LOG10(W0) 
            RETURN 
         END IF 
! 
!     Otherwise, normal inter- or extrapolation 
! 
!     Both interpolations (in T as well as in electron density) are 
!     by default the quadratic interpolations in logarithms 
! 
         DO IX=1,NT2-1 
            IPX=IX 
            IF(X0.LE.XT2(IX+1)) EXIT 
         END DO 
         N0X=IPX-NX/2+1 
         IF(N0X.LT.1) N0X=1 
         IF(N0X.GT.NT2-NX+1) N0X=NT2-NX+1 
         N1X=N0X+NX-1 
         DO IX=N0X,N1X 
            I0=IX-N0X+1 
            XX(I0)=XT2(IX) 
            WX(I0)=PRF2(IWL,IX,IZZ) 
         END DO 
         WZ(I0Z)=YINT(XX,WX,X0) 
      END DO ZINT 
      W0=YINT(ZZ,WZ,Z0) 
      RETURN 
      END SUBROUTINE INTHE2 
! 
! ******************************************************************** 
! 
! 
 
      SUBROUTINE DIVHE2(A,DIV) 
!     ======================== 
! 
!     Auxiliary procedure for evaluating approximate Stark profile 
!     for He II lines 
!     This procedure is quite analogous to DIVSTR for hydrogen; 
!     the only difference is a somewhat different definition 
!     of the parameter A ,ie. A for He II is equal to A for hydrogen 
!     minus ln(2) 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP),PARAMETER :: UN=1.,TWO=2.,UNQ=1.25,UNH=1.5,TWH=2.5,        & 
     &                      CA=0.978,BL=5.821,AL=1.26,CX=0.28,DX=0.0001 
! 
      A=UNH*LOG(BETAD)-CA 
      IF(BETAD.LT.BL) RETURN 
      IF(A.GE.AL) THEN 
         X=SQRT(A)*(UN+UNQ*LOG(A)/(4.*A-5.)) 
      ELSE 
         X=SQRT(CX+A) 
      ENDIF 
      DO I=1,5 
         XN=X*(UN-(X*X-TWH*LOG(X)-A)/(TWO*X*X-TWH)) 
         IF(ABS(XN-X).LE.DX) EXIT 
         X=XN 
      END DO 
      DIV=X 
      RETURN 
      END SUBROUTINE DIVHE2 
! 
! ******************************************************************** 
! 
! 
 
      SUBROUTINE PHE2(ISPEC,ID,ABLIN,EMLIN) 
!     ===================================== 
! 
!     Evaluation of the opacity and emissivity in a given He II line, 
!     using profile coefficients calculated by Schoening and Butler. 
! 
!     Input: ISPEC - line index, defined in HE2INI 
!            ID    - depth index 
!     Output: ABLIN - absorption coefficient 
!             EMLIN - emission coefficient 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use heprf 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) ::  ABLIN(MFREQ),EMLIN(MFREQ),OSCHE2(19),                  & 
     &             PRF0(40),WLL(40) 
      DATA OSCHE2/6.407E-1, 1.506E-1, 5.584E-2, 2.768E-2,                 & 
     &        1.604E-2, 1.023E-2, 6.980E-3,                               & 
     &        8.421E-1, 3.230E-2, 1.870E-2, 1.196E-2, 8.187E-3,           & 
     &        5.886E-3, 4.393E-3, 3.375E-3, 2.656E-3,                     & 
     &        1.038,    1.793E-1, 6.549E-2/ 
! 
!     ILINE - line index 
! 
      ILINE=ISPEC-5 
! 
      DO IWL=1,NWLHE2(ILINE) 
         PRF0(IWL)=PRFHE2(ILINE,ID,IWL) 
         WLL(IWL)=WLHE2(ILINE,IWL) 
      END DO 
! 
      I=ILHE2(ILINE) 
      J=IUHE2(ILINE) 
      II=I*I 
      JJ=J*J 
!     IF(I.LE.2) THEN 
!        WLIN=227.838/(1./II-1./JJ) 
!      ELSE 
!        WLIN=227.7776/(1./II-1./JJ) 
!     END IF 
! 
      WL00=227.9384 
      IF(ILHE2(ILINE).GE.3.AND.VACLIM.LE.2001.) WL00=227.7776 
      WLIN=WL00/(1./ILHE2(ILINE)**2-1./IUHE2(ILINE)**2) 
! 
      T=TEMP(ID) 
! 
!     He III population (either LTE or NLTE, depending on input model) 
! 
      IF(IELHE2.GT.0.and.inlte.gt.0) THEN 
         PP=POPUL(NNEXT(IELHE2),ID) 
         NLHE2=NLAST(IELHE2)-NFIRST(IELHE2)+1 
       ELSE 
         PP=RRR(ID,3,2) 
         NLHE2=0 
      END IF 
! 
!     population of the lower level of the given transition 
!     (again either LTE or NLTE) 
! 
      PP=PP*ELEC(ID)*4.1412E-16/T/SQRT(T)*II 
      IF(I.LE.NLHE2.and.inlte.gt.0) THEN 
         POPI=POPUL(NFIRST(IELHE2)+I-1,ID) 
       ELSE 
         POPI=PP*EXP(631479./T/II) 
      END IF 
! 
!     population of the upper level of the given transition 
!     (again either LTE or NLTE) 
! 
      IF(J.LE.NLHE2) THEN 
         POPJ=POPUL(NFIRST(IELHE2)+J-1,ID)*II/JJ 
       ELSE 
         POPJ=PP*EXP(631479./T/JJ) 
      END IF 
 
! 
!     loop over frequency points - opacity and emissivity in the given line 
!     absorption coefficent is found by interpolating in previously 
!     calculated tables, based on calculations of Schoening and Butler 
!     (see procedure HE2INI) 
! 
      FID=0.02654*OSCHE2(ILINE) 
      DO IJ=3,NFREQ 
         AL=ABS(WLAM(IJ)-WLIN) 
         IF(AL.LT.1.E-4) AL=1.E-4 
         AL=LOG10(AL) 
         DO IWL=1,NWLHE2(ILINE)-1 
            IW0=IWL 
            IF(AL.LE.WLL(IWL+1)) EXIT 
         END DO 
         IW1=IW0+1 
         PRH=(PRF0(IW0)*(WLL(IW1)-AL)+PRF0(IW1)*(AL-WLL(IW0)))/           & 
     &       (WLL(IW1)-WLL(IW0)) 
         SG=EXP(PRH*2.3025851)*FID 
         if((popi-popj).le.0. .and. lasdel) cycle 
         ABLIN(IJ)=ABLIN(IJ)+SG*(POPI-POPJ) 
         EMLIN(IJ)=EMLIN(IJ)+SG*POPJ*1.4747E-2*(FREQ(IJ)*1.E-15)**3 
      END DO 
      RETURN 
      END SUBROUTINE PHE2 
! 
! ******************************************************************** 
! 
! 
 
      FUNCTION ISPEC(IAT,ION,ALAM) 
!     ============================ 
! 
!     Auxiliary procedure for INISET 
! 
!     Input:  IAT  - atomic number 
!             ION  - ion (=1 for neutrals, =2 for once ionized, etc.) 
!             ALAM - wavelength in nanometers 
!     Output: ISPEC - parameter specifying whether the given line 
!                     is taken with a special (pretabulated) absorption 
!                     profile - only for hydrogen and helium 
!                   = 0  - profile is taken as an ordinary Voigt profile 
!                   > 0  - special profile 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
! 
      ISPEC=0 
      IF(IAT.GT.2) RETURN 
! 
      IF(IAT.EQ.1) THEN 
         ISPEC=1 
         RETURN 
       ELSE 
         IF(ION.EQ.1) THEN 
            IF(ABS(ALAM-447.1).LT.0.5.AND.IHE1PR.GT.0) ISPEC=2 
            IF(ABS(ALAM-438.8).LT.0.2.AND.IHE1PR.GT.0) ISPEC=3 
            IF(ABS(ALAM-402.6).LT.0.2.AND.IHE1PR.GT.0) ISPEC=4 
            IF(ABS(ALAM-492.2).LT.0.2.AND.IHE1PR.GT.0) ISPEC=5 
          ELSE 
! 
            IF(ALAM.LT.163..OR.ALAM.GT.1012.7) RETURN 
            IF(ALAM.LT.321.) THEN 
               IF(ABS(ALAM-164.0).LT.0.2.AND.IHE2PR.GT.0) ISPEC=6 
               IF(ABS(ALAM-320.3).LT.0.2.AND.IHE2PR.GT.0) ISPEC=7 
               IF(ABS(ALAM-273.3).LT.0.2.AND.IHE2PR.GT.0) ISPEC=8 
               IF(ABS(ALAM-251.1).LT.0.2.AND.IHE2PR.GT.0) ISPEC=9 
               IF(ABS(ALAM-238.5).LT.0.2.AND.IHE2PR.GT.0) ISPEC=10 
               IF(ABS(ALAM-230.6).LT.0.2.AND.IHE2PR.GT.0) ISPEC=11 
               IF(ABS(ALAM-225.3).LT.0.2.AND.IHE2PR.GT.0) ISPEC=12 
             ELSE IF(ALAM.LT.541.) THEN 
               IF(ALAM.LT.392.3) RETURN 
               IF(ABS(ALAM-468.6).LT.0.2.AND.IHE2PR.GT.0) ISPEC=13 
               IF(ABS(ALAM-485.9).LT.0.2.AND.IHE2PR.GT.0) ISPEC=14 
               IF(ABS(ALAM-454.2).LT.0.2.AND.IHE2PR.GT.0) ISPEC=15 
               IF(ABS(ALAM-433.9).LT.0.2.AND.IHE2PR.GT.0) ISPEC=16 
               IF(ABS(ALAM-420.0).LT.0.2.AND.IHE2PR.GT.0) ISPEC=17 
               IF(ABS(ALAM-410.0).LT.0.2.AND.IHE2PR.GT.0) ISPEC=18 
               IF(ABS(ALAM-402.6).LT.0.2.AND.IHE2PR.GT.0) ISPEC=19 
               IF(ABS(ALAM-396.8).LT.0.2.AND.IHE2PR.GT.0) ISPEC=20 
               IF(ABS(ALAM-392.3).LT.0.2.AND.IHE2PR.GT.0) ISPEC=21 
             ELSE 
               IF(ABS(ALAM-1012.4).LT.0.2.AND.IHE2PR.GT.0) ISPEC=22 
               IF(ABS(ALAM-656.0).LT.0.2.AND.IHE2PR.GT.0) ISPEC=23 
               IF(ABS(ALAM-541.2).LT.0.2.AND.IHE2PR.GT.0) ISPEC=24 
            END IF 
         END IF 
      END IF 
      RETURN 
      END FUNCTION ISPEC 
! 
! 
!     ****************************************************************** 
! 
! 
 
      SUBROUTINE HESET(IL,ALM,EXCL,EXCU,ION,IPRF0,ILWN,IUPN) 
!     ====================================================== 
! 
!     Auxiliary procedure for INISET - set up quantities: 
!     IPRF0      - index for the procedure evaluating standard absorption 
!                  profile coefficient for He I lines - see GAMHE 
!     ILWN,IUPN  - only in NLTE option is switched on; 
!                  indices of the lower and upper level associated with 
!                  the given line 
! 
!     Input: IL - line index 
!            ALM - line wavelength in nm 
!            EXCL - excitation potential of the lower level (in cm**-1) 
!            EXCU - excitation potential of the upper level (in cm**-1) 
!            ION  - ionisation degree (1=neutrals, 2=once ionized, etc.) 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      INTEGER ::  JU(24),NU(24),IT(24) 
      DATA IT/1,1,0,1,0,0,0,1,0,0,0,1,1,0,0,0,1,0,1,0,0,0,0,0/ 
      DATA NU/6,6,9,3,8,4,7,5,6,6,5,4,4,4,3,4,3,3,5,5,7,8,10,2/ 
      DATA JU/15,3,5,9,5,3,5,3,5,1,1,15,3,5,3,1,15,5,15,5,1,1,1,9/ 
! 
! ******* He II  *********** 
! 
      IF(ION.EQ.2) THEN 
         IF(IELHE2.LE.0) RETURN 
         N0I=NFIRST(IELHE2) 
         NLHE2=NLAST(IELHE2)-N0I+1 
         XL=SQRT(1./(1.-EXCL/438916.146)) 
         ILW=INT(XL) 
         IF((FLOAT(ILW)-XL).LT.0.) ILW=ILW+1 
         XU=SQRT(1./(1.-EXCU/438916.146)) 
         IUN=INT(XU) 
         IF((FLOAT(IUN)-XU).LT.0.) IUN=IUN+1 
         IF(ILW.LE.NLHE2) ILWN=ILW+N0I-1 
         IF(IUN.LE.NLHE2) IUPN=IUN+N0I-1 
         RETURN 
      END IF 
! 
! ******* He I  *********** 
! 
!     switch IPRF0 - see GAMHE 
! 
      IL1=IL 
      ALAM=ALM*10. 
      IPRF=0 
      IF(ABS(ALAM-3819.60).LT.1.) IPRF=1 
      IF(ABS(ALAM-3867.50).LT.1.) IPRF=2 
      IF(ABS(ALAM-3871.79).LT.1.) IPRF=3 
      IF(ABS(ALAM-3888.65).LT.1.) IPRF=4 
      IF(ABS(ALAM-3926.53).LT.1.) IPRF=5 
      IF(ABS(ALAM-3964.73).LT.1.) IPRF=6 
      IF(ABS(ALAM-4009.27).LT.1.) IPRF=7 
      IF(ABS(ALAM-4120.80).LT.1.) IPRF=8 
      IF(ABS(ALAM-4143.76).LT.1.) IPRF=9 
      IF(ABS(ALAM-4168.97).LT.1.) IPRF=10 
      IF(ABS(ALAM-4437.55).LT.1.) IPRF=11 
      IF(ABS(ALAM-4471.50).LT.1.) IPRF=12 
      IF(ABS(ALAM-4713.20).LT.1.) IPRF=13 
      IF(ABS(ALAM-4921.93).LT.1.) IPRF=14 
      IF(ABS(ALAM-5015.68).LT.1.) IPRF=15 
      IF(ABS(ALAM-5047.74).LT.1.) IPRF=16 
      IF(ABS(ALAM-5875.70).LT.1.) IPRF=17 
      IF(ABS(ALAM-6678.15).LT.1.) IPRF=18 
      IF(ABS(ALAM-4026.20).LT.1.) IPRF=19 
      IF(ABS(ALAM-4387.93).LT.1.) IPRF=20 
      IF(ABS(ALAM-4023.97).LT.1.) IPRF=21 
      IF(ABS(ALAM-3935.91).LT.1.) IPRF=22 
      IF(ABS(ALAM-3833.55).LT.1.) IPRF=23 
      IF(ABS(ALAM-10830.0).LT.1.) IPRF=24 
      IF(IPRF.GT.0.AND.IPRF.LE.20) IPRF0=IPRF 
! 
!     Indices of NLTE levels associated with the given line 
! 
      IF(INLTE.gt.5.OR.IELHE1.EQ.0) RETURN 
      N0I=NFIRST(IELHE1) 
      N1I=NLAST(IELHE1) 
      HC=CL*H 
      EION=ENION(N0I)/HC 
      ILW=0 
      IUN=0 
      NQL=0 
 
      IF(IPRF.GT.0) NQL=NU(IPRF) 
      DO I=N0I,N1I 
         NQ=NQUANT(I) 
         EX=EION-ENION(I)/HC 
         IF(ABS(EXCL-EX).LT.100.) THEN 
            ILW=I 
            IGL=INT(G(I)+0.001) 
         END IF 
         IF(NQ.EQ.NQL) THEN 
            IG=INT(G(I)+0.001) 
            IF(IT(IPRF).EQ.0) THEN 
               IF(NQ.EQ.2.AND.IG.EQ.JU(IPRF)) IUN=I 
               IF(NQ.EQ.3) THEN 
                  IF(IG.EQ.JU(IPRF)) THEN 
                     IF(IG.EQ.1.OR.IG.EQ.5) IUN=I 
                     IF(IG.EQ.3.AND.IGL.EQ.1) IUN=I 
                   ELSE 
                     IF(IG.EQ.9) IUN=I 
                  END IF 
               END IF 
               IF(NQ.EQ.4) THEN 
                  IF(IG.EQ.JU(IPRF)) THEN 
                     IF(IG.EQ.1.OR.IG.EQ.5.OR.IG.EQ.7) IUN=I 
                     IF(IG.EQ.3.AND.IGL.EQ.1) IUN=I 
                   ELSE 
                     IF(IG.EQ.16) IUN=I 
                  END IF 
               END IF 
               IF(IG.EQ.25.OR.IG.EQ.36) IUN=I 
               IF(IG.EQ.49.OR.IG.EQ.64.OR.IG.EQ.81) IUN=I 
               IF(IG.EQ.100.OR.IG.EQ.121.OR.IG.EQ.144) IUN=I 
             ELSE 
               IF(NQ.EQ.3) THEN 
                  IF(IG.EQ.JU(IPRF)) THEN 
                     IF(IG.EQ.9.OR.IG.EQ.15) IUN=I 
                     IF(IG.EQ.3.AND.IGL.EQ.9) IUN=I 
                   ELSE 
                     IF(IG.EQ.27) IUN=I 
                  END IF 
               END IF 
               IF(NQ.EQ.4) THEN 
                  IF(IG.EQ.JU(IPRF)) THEN 
                     IF(IG.EQ.9.OR.IG.EQ.15.OR.IG.EQ.21) IUN=I 
                     IF(IG.EQ.3.AND.IGL.EQ.9) IUN=I 
                   ELSE 
                     IF(IG.EQ.48) IUN=I 
                  END IF 
               END IF 
               IF(IG.EQ.75) IUN=I 
               IF(IG.EQ.108.OR.IG.EQ.147.OR.IG.EQ.192) IUN=I 
               IF(IG.EQ.243.OR.IG.EQ.300.OR.IG.EQ.363) IUN=I 
            END IF 
            IF(NQ.EQ.2.AND.IG.EQ.16) IUN=I 
            IF(NQ.EQ.3.AND.IG.EQ.36) IUN=I 
            IF(NQ.EQ.4.AND.IG.EQ.64) IUN=I 
            IF(NQ.EQ.5.AND.IG.EQ.100) IUN=I 
            IF(NQ.EQ.6.AND.IG.EQ.144) IUN=I 
            IF(NQ.EQ.7.AND.IG.EQ.196) IUN=I 
            IF(NQ.EQ.8.AND.IG.EQ.256) IUN=I 
            IF(NQ.EQ.9.AND.IG.EQ.324) IUN=I 
            IF(NQ.EQ.10.AND.IG.EQ.400) IUN=I 
         END IF 
      END DO 
      ILWN=ILW 
      IUPN=IUN 
! 
      RETURN 
      END SUBROUTINE HESET 
! 
! 
! ******************************************************************** 
! 
      SUBROUTINE INISET 
!     ================= 
! 
!     SELECTION OF LINES THAT MAY CONTRIBUTE, 
!     SET UP AUXILIARY FIELDS CONTAINING LINE PARAMETERS, 
!     SET UP THE SET OF FREQUENCY POINTS 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use lindat 
      use wincom 
      use molist 
 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      INTEGER, SAVE :: ILLAST 
! 
      DATA CNM,CAS /2.997925e17,2.997925e18/ 
! 
      DO I=1,MFRQ 
         W(I)=0. 
         IJCTR(I)=0 
      END DO 
! 
      IL0=0 
      IPRSET=0 
      NLIN=0 
      IREADP=1 
      IRLIST=0 
      IF(IBLANK.LE.1.OR.IMODE.EQ.1.OR.IMODE.EQ.-1) IREADP=0 
      IF(IBLANK.LE.1) APREV=0. 
      FRMIN=CNM/ALAM0 
      FRM=FRMIN 
      if(ifwin.le.0) then 
         ij0=3 
       else 
         ij0=1 
      end if 
      IJ=IJ0 
      FREQ(IJ0)=FRM 
      SPACE=SPACE0 
!     IF(ALAMC.GT.0.) SPACE=SPACE0*ALAM0/ALAMC 
      IF(SPACE0.LT.0.) SPACE=-SPACE0 
! 
!     ==== for IMODE=2, i.,e. the continuum-only mode, selecting frequencies 
! 
      IMOD2: IF(IMODE.EQ.2) THEN 
         NFRP=NFREQS+1 
         W0=SPACE 
         FRACT=FREQ(IJ) 
         ALACT=CNM/FRACT 
! 
         DO K=1,NFRP 
            FRACT=FRACT-W0 
            ALACT=ALACT+W0 
            IJ=IJ+1 
            IF(IJ.GT.NFREQS) EXIT 
            FREQ(IJ)=CNM/ALACT 
            W(IJ)=W(IJ)+(FREQ(IJ-1)-FREQ(IJ))*0.5 
            W(IJ-1)=W(IJ-1)+(FREQ(IJ-1)-FREQ(IJ))*0.5 
         END DO 
 
         FRMAX=FREQ(NFREQS) 
         ALAM1=CNM/FRMAX 
         NFREQ=NFREQS 
         IJMAX=NFREQS 
         NBLANK=IBLANK+1 
 
       ELSE 
! 
!        ==== for IMODE.ne.2 - standard synthetic spectrum mode; 
!             selection of lines and corresponding frequencies 
! 
         ISTR=0 
         IJMAX=0 
         if(ifwin.le.0) then 
            CUTOFF=CUTOF0 
            DOPSTD=1.E7/ALAM0*DSTD 
            DISTAN=0.15*DOPSTD 
            SPAC=3.E16/ALAM0/ALAM0*SPACE 
            DISTA0=0.14*SPAC 
            ASTD=1.0 
            AVAB=ABSTD(IDSTD)*RELOP 
         end if 
         FRLI0=FRMIN 
         IF(IBLANK.GE.2.AND.IMODE.EQ.-1) IL0=ILLAST 
! 
!  ****  loop over contributing lines 
! 
         LINES: DO 
! 
!           set up indices of lines 
!           IL0 - the current index of line in the numbering of all lines 
! 
            IF(IREADP.EQ.1) THEN 
               IPRSET=IPRSET+1 
               IL0=INDLIP(IPRSET) 
               IF(FREQ0(IL0).LT.FRMIN) THEN 
                  IREADP=0 
                  IL0=INDLIP(IPRSET-1)+1 
               END IF 
             ELSE 
               IL0=IL0+1 
            END IF 
            IF(IL0.GT.NLIN0) EXIT LINES 
            FRLIM=FRLI0 
            FR0=FREQ0(IL0) 
            ALAM=CNM/FR0 
!!          write(*,*) 'il0,alam',il0,alam,alam0-cutoff,alam1+cutoff 
 
            if(ifwin.gt.0) then 
               IF(ALAMC.GT.0.) SPACE=SPACE0*ALAM/ALAMC 
               IF(SPACE0.LT.0.) SPACE=-SPACE0 
               CUTOFF=CUTOF0*ALAM/ALAMC 
               DOPSTD=1.E7/ALAM*DSTD 
               DISTAN=0.15*DOPSTD 
               SPAC=SPACE 
               IF(MOD(IFREQ,10).GT.0) SPAC=3.E16/ALAM/ALAM*SPACE 
               DISTA0=0.14*SPAC 
            end if 
! 
!           set up a different starting wavelength for IMODE=1 
! 
            IF(IMODE.EQ.1.AND.ISTR.NE.1.AND.IJ.EQ.3.AND.                  & 
     &         ALAM.GE.ALAM0+2.*CUTOFF) THEN 
               ALAM0=ALAM-CUTOFF+0.0001 
               FRMIN=CNM/ALAM0 
               FRM=FRMIN 
               IJ=IJ0 
               FREQ(IJ0)=FRM 
            END IF 
            IF(ALAM.LT.ALAM0-CUTOFF) CYCLE LINES 
            IF(IJ.GT.NFREQS.AND.ALAM.GT.ALAM1+CUTOFF)  THEN 
               NBLANK=IBLANK+1 
!!          write(*,"('ij,nfr,al,al+c,iblank,nblank',2i6,2f12.3,2i8)") 
!!   *         ij,nfreqs,alam*10.,(alam1+cutoff)*10.,iblank,nblank 
               EXIT LINES 
            END IF 
! 
!           SECOND SELECTION : FOR LINE STRENGHTS 
! 
            ISTR=0 
            IF(IMODE.GE.1) THEN 
               ISTR=1 
             ELSE 
               EXT=EXTIN(IL0) 
               FRLI0=FR0-EXT-SPAC 
               IF(FRLI0.GT.FRLIM) FRLI0=FRLIM 
               FRMIV=FRMIN 
               if(ifwin.gt.0) frmiv=frmiv*(1.+vinf/2.997925e10) 
               IF(ALAM.LT.ALAM0.AND.FR0-FRMIV.GT.EXT+SPAC) CYCLE LINES 
               ISTR=1 
               FRMAV=FRMAX 
               if(ifwin.gt.0) frmav=frmav*(1.-vinf/2.997925e10) 
               IF(IJ.GE.NFREQS+1.AND.FRMAv-FR0.GT.EXT+SPAC) CYCLE LINES 
            END IF 
 
! 
            NLIN=NLIN+1 
            if(nlin.gt.mlin) call quit(' too many lines in a set') 
            INDLIN(NLIN)=IL0 
            ALAMCU=ALAM+CUTOFF 
            IF(IJ.GE.NFREQS+1) CYCLE LINES 
            IF(FR0.GT.FRMIN) CYCLE LINES 
! 
!           FREQUENCY POINTS AND WEIGHTS 
! 
            DELT=ABS(FRM-FR0) 
            IF(DELT.LT.DISTA0.AND.IMODE.NE.1) CYCLE LINES 
            DFREL=CNM*(1./FR0-1./FRM)/SPACE 
            NFRP=int(DFREL)+1 
            IF(NFRP.LE.2) NFRP=2 
            W0=CNM*(1./FR0-1.D0/FRM)/NFRP 
            FRM=FR0 
            FRACT=FREQ(IJ) 
            ALACT=CNM/FRACT 
! 
            FRQ: DO K=1,NFRP 
               FRACT=FRACT-W0 
               ALACT=ALACT+W0 
               IF(IMODE.LT.1.AND.NFRP.NE.2) THEN 
                  IF(FRACT.LT.FRLIM.AND.FRACT.GT.FR0+EXT+SPAC) THEN 
                     CYCLE FRQ 
                  END IF 
               END IF 
               IJ=IJ+1 
!              write(*,*) 'frq',k,ij,alact 
               IF(IJ.GT.NFREQS) EXIT FRQ 
               FREQ(IJ)=CNM/ALACT 
               W(IJ)=W(IJ)+(FREQ(IJ-1)-FREQ(IJ))*0.5 
               W(IJ-1)=W(IJ-1)+(FREQ(IJ-1)-FREQ(IJ))*0.5 
            END DO FRQ 
!           write(*,"('EFR',3i8,1pe10.3,0pf12.3)") 
!    *      ij,nfreqs,nfrp,alact*10.,alamcu*10. 
 
            IF(IJ.LE.NFREQS) THEN 
               IJCTR(IJ)=IL0 
!              write(*,*) '*** ij,ijctr',ij,ijctr(ij) 
               DISTA0=DISTAN 
               CYCLE LINES 
             ELSE 
               FRMAX=FREQ(NFREQS) 
               ALAM1=CNM/FRMAX 
               NFREQ=NFREQS 
!              write(*,*) '++++ ij,nfreqs,alam1',ij,nfreqs,alam1 
               IF(IMODE.EQ.2) EXIT LINES 
               CYCLE LINES 
            END IF 
         END DO LINES 
! 
         NBLANK=IBLANK+1 
         IF(IJ.GE.NFREQS+1) THEN 
            IJMAX=NFREQS 
            NFREQ=NFREQS 
          ELSE 
            IJMAX=IJ 
            IJMAX=MIN(IJMAX,NFREQS) 
            NFREQ=IJMAX 
         END IF 
 
!!       write(*,"('end loop lines',6i8,1p2e13.5)") 
!!   *   ij,ijmax,nfreq,nfreqs,iblank,nblank,freq(ijmax),frlast 
 
         IF(FREQ(IJMAX).LE.FRLAST) NBLANK=IBLANK 
!!       write(*,*) 'alm00',alm00 
!!       write(*,*) 'alm01',freq(ijmax),0.999999*cnm/alm00 
         if(alm00.gt.0.) then 
            if(freq(ijmax).ge.0.999999*cnm/alm00.and.iblank.gt.1)         & 
     &         nblank=iblank 
         end if 
 
!!       write(*,"(' iblank,nblank,alm00',2i8,f12.3)") 
!!   *     iblank,nblank,alm00*10. 
! 
!        correction for molecular lines 
! 
!!       write(*,*) 'nmlist,ifmol',nmlist,ifmol 
         if(nmlist.gt.0.and.ifmol.gt.0) then 
            do ilist=1,nmlist 
               if(alastm(ilist).gt.0..and.alastm(ilist).le.alact) then 
                  nblank=iblank 
                  irlist=1 
!!              write(*,"(' MOL',i4,2f12.3,3i8/)") 
!!   *          ilist,10.*alastm(ilist),10.*alact,iblank,nblank,irlist 
               end if 
            end do 
         end if 
 
      END IF IMOD2 
! 
      IF(IFWIN.LE.0) THEN 
         FREQ(1)=FREQ(3) 
         FREQ(2)=FREQ(IJMAX) 
         W(1)=0.5*(FREQ(1)-FREQ(2)) 
         W(2)=W(1) 
      END IF 
! 
!     truncate the interval if the required end is reached 
! 
      ijmx=2 
      if(ifwin.gt.0) ijmx=ijmax 
      IF(FREQ(ijmx).LT.FRLAST) THEN 
         FREQ(ijmx)=FRLAST 
         if(ifwin.le.0) then 
         W(1)=0.5*(FREQ(1)-FREQ(2)) 
         W(2)=W(1) 
         end if 
         DO IJ=IJ0,NFREQ 
            IF(FREQ(IJ).LT.FRLAST) EXIT 
            IJMAX=IJ 
         END DO 
         NFREQ=IJMAX+1 
         FREQ(NFREQ)=FRLAST 
         W(NFREQ)=0.5*(FREQ(NFREQ-1)-FREQ(NFREQ)) 
         W(NFREQ-1)=W(NFREQ)+0.5*(FREQ(NFREQ-2)-FREQ(NFREQ-1)) 
      END IF 
! 
!     frequency interpolation coefficients 
! 
      IMODM1: IF(IMODE.NE.-1) THEN 
         if(ifwin.le.0) then 
            XX=FREQ(2)-FREQ(1) 
            DO IJ=1,NFREQ 
               WLAM(IJ)=2.997925E18/FREQ(IJ) 
               FRX1(IJ)=(FREQ(IJ)-FREQ(1))/XX 
               FRX2(IJ)=(FREQ(2)-FREQ(IJ))/XX 
            END DO 
          else 
            DO IJ=1,NFREQ 
               WLAM(IJ)=CAS/FREQ(IJ) 
               frqobs(ij)=freq(ij) 
               wlobs(ij)=wlam(ij) 
               fr=freq(ij) 
               BNUE(IJ)=BN*fr*fr*fr 
               DO IJCI=1,NFREQC-1 
                  IF(WLAM(IJ).LE.WLAMC(IJCI)) EXIT 
               END DO 
               IJC=IJCI 
               IJCINT(IJ)=MAX(IJC-1,1) 
               IJCI=IJCINT(IJ) 
               FRX1(IJ)=(FREQ(IJ)-FREQC(IJCI+1))/                         & 
     &            (FREQC(IJCI)-FREQC(IJCI+1)) 
            END DO 
            nfrobs=nfreq 
            xx=freq(nfreq)-freq(1) 
         end if 
! 
!        frequency indices of the line centers 
! 
         DFRCON=real(NFREQ-ij0) 
         DFRCON=-DFRCON/XX 
         IFRCON=INT(DFRCON) 
         LINC: DO IL=1,NLIN 
            fr0=freq0(indlin(il)) 
            XJC=3.+DFRCON*(FREQ(1)-FR0) 
            IJC=INT(XJC) 
            IJCNTR(IL)=IJC 
            if(ijc.le.ij0.or.ijc.ge.nfreq) cycle linc 
            if(fr0.lt.freq(ijc)) then 
               ijc0=ijc 
               dfr0=freq(ijc0)-fr0 
               hig: do 
                  ijc0=ijc0+1 
                  dfr=abs(freq(ijc0)-fr0) 
                  if(dfr.lt.dfr0) then 
                     ijc=ijc0 
                     ijc0=ijc0+1 
                     dfr0=dfr 
                     cycle hig 
                   else 
                     exit hig 
                  end if 
               end do hig 
             else if(fr0.gt.freq(ijc)) then 
               ijc0=ijc 
               dfr0=fr0-freq(ijc0) 
               low: do 
                  ijc0=ijc0-1 
                  dfr=abs(freq(ijc0)-fr0) 
                  if(dfr.lt.dfr0) then 
                     ijc=ijc0 
                     ijc0=ijc0-1 
                     dfr0=dfr 
                     cycle low 
                   else 
                     exit low 
                  end if 
              end do low 
            end if 
            IJCNTR(IL)=IJC 
         END DO LINC 
      END IF IMODM1 
!!    write(*,*) 'after linc' 
! 
      if(ifwin.gt.0) then 
! 
!     set up switches for hydrogen and He II line opacity 
! 
         DO  IJ=1,NFREQ 
            call hylsew(ij) 
            call he2sew(ij) 
         END DO 
      end if 
! 
      NSP=0 
      DO IL=1,NLIN 
         IL0=INDLIN(IL) 
         ISP=ISPRF(IL0) 
         IF(ISP.GT.5) THEN 
            NSP=NSP+1 
            ISP0(NSP)=ISP 
         END IF 
         INDLIP(IL)=INDLIN(IL) 
      END DO 
      if(ifwin.le.0) then 
         ILLAST=INDLIN(NLIN) 
       else 
         ILLAST=0 
         IF(NLIN.GT.0) ILLAST=INDLIN(NLIN) 
      end if 
!!       write(*,"('8iblank,nblank,alm00',2i8,f12.3)") 
!!   *     iblank,nblank,alm00*10. 
! 
      CALL READPH 
! 
      IF(ALAM0.LE.APREV+0.001) NBLANK=IBLANK 
!!       write(*,"('9iblank,nblank,alm00',2i8,f12.3)") 
!!   *     iblank,nblank,alm00*10. 
      APREV=ALAM0 
      ALAM0=ALAM1 
      ALM00=CNM/FREQ(NFREQ) 
!!    write(66,"('===== alam0,aprev,alm00,iblank,nblank',3f12.3,2i8)") 
!!   *   alam0,aprev,alm00,iblank,nblank 
! 
      RETURN 
      END SUBROUTINE INISET 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE READPH 
!     ================= 
! 
!     Auxiliary routine for LINSET - read table of detailed 
!     photoinization cross-section from unit IPHT1, 
!     and interpolate to the set of current wavelengths (WLAM) 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use lindat 
      use photcs 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: PHT0(MPHOT),PHT1(MPHOT) 
      INTEGER, SAVE  :: IPHT(MPHOT),IEND(MPHOT),                          & 
     &            IFILE(MPHOT),NELEM(MPHOT),INDEX(MPHOT,MPHOT) 
      INTEGER, PARAMETER :: IPHT0=57 
      INTEGER            :: IALLOC
      DATA IALLOC /1/ 
! 
!     initialization - read basic information about files where the 
!                      cross-sections are stored, 
!                      and basic parameters for starting levels 
! 
      IF(IALLOC.EQ.1) THEN
         CALL ALLOC_PHOTCS
         IALLOC=0
      END IF

      IF(IBLANK.LE.1) THEN 
         NPHT=0 
         IPHT1=0 
         NUMFIL=0 
!        DO IJ=1,MFRQ 
!           DO I=1,MPHOT 
!              PHOT(IJ,I)=0. 
!           END DO 
!        END DO 
         READ(IPHT0,*,IOSTAT=IOS) NPHT 
         IF(IOS.NE.0) RETURN 
         IF(NPHT.LE.0) RETURN 
         CALL ALLOC_PHOTCS 
         DO IJ=1,MFRQ 
            DO I=1,MPHOT 
               PHOT(IJ,I)=0. 
            END DO 
         END DO 
         npht1=npht 
         READ(IPHT0,*) (IPHT(I),I=1,NPHT) 
         READ(IPHT0,*) (APHT(I),I=1,NPHT) 
         READ(IPHT0,*) (EPHT(I),I=1,NPHT) 
         READ(IPHT0,*) (GPHT(I),I=1,NPHT) 
         READ(IPHT0,*) (JPHT(I),I=1,NPHT) 
! 
!     determination of the number of files (NFILE) and the 
!     partitioning of the individual cross-section to the corresponding 
!     files 
! 
         NUMFIL=1 
         IFILE(1)=1 
         NELEM(1)=1 
         INDEX(1,1)=1 
         IF(NPHT.GT.1) THEN 
            ILOOP: DO I=2,NPHT 
               JLOOP: DO J=1,I-1 
                  IF(IPHT(I).EQ.IPHT(J)) THEN 
                     IFILE(I)=IFILE(J) 
                     NELEM(IFILE(I))=NELEM(IFILE(I))+1 
                     INDEX(IFILE(I),NELEM(IFILE(I)))=I 
                     CYCLE ILOOP 
                  END IF 
               END DO JLOOP 
               NUMFIL=NUMFIL+1 
               IFILE(I)=NUMFIL 
               NELEM(NUMFIL)=1 
               INDEX(NUMFIL,1)=I 
            END DO ILOOP 
         END IF 
         DO IFIL=1,NUMFIL 
            IEND(IFIL)=0 
         END DO 
      END IF 
      IF(NUMFIL.LE.0) RETURN 
! 
!     loop over individual files containing the photoionization data 
! 
      READFILES: DO IFIL=1,NUMFIL 
         IF(IEND(IFIL).EQ.1) CYCLE READFILES 
         IF(IEND(IFIL).EQ.2) CYCLE READFILES 
         NPHT1=NELEM(IFIL) 
         IPHT1=IPHT(INDEX(IFIL,1)) 
         IF(IBLANK.LE.1) THEN 
            IPHT1LOOP: DO 
               READ(IPHT1,*,IOSTAT=IOS2) WPHT1,(PHT1(I),I=1,NPHT1) 
               IF(IOS2.EQ.0) THEN 
                  IF(WPHT1.LT.WLAM(1)) THEN 
                     CYCLE IPHT1LOOP 
                   ELSE 
                     EXIT IPHT1LOOP 
                  END IF 
               END IF 
            END DO IPHT1LOOP 
            BACKSPACE(IPHT1) 
            BACKSPACE(IPHT1) 
            READ(IPHT1,*,IOSTAT=IOS2) WPHT0,(PHT0(I),I=1,NPHT1) 
          ELSE 
            BACKSPACE(IPHT1) 
            BACKSPACE(IPHT1) 
            READ(IPHT1,*,IOSTAT=IOS3) WPHT0,(PHT0(I),I=1,NPHT1) 
            IF(IOS3.EQ.0) THEN 
               READ(IPHT1,*,IOSTAT=IOS4) WPHT1,(PHT1(I),I=1,NPHT1) 
            END IF 
         END IF 
         DW=WPHT1-WPHT0 
         A1=(WPHT1-WLAM(3))/DW 
         A2=(WLAM(3)-WPHT0)/DW 
         DO I=1,NPHT1 
            INDX=INDEX(IFIL,I) 
            PHOT(1,INDX)=0. 
            PHOT(2,INDX)=0. 
            PHOT(3,INDX)=(A1*PHT0(I)+A2*PHT1(I))*1.E-18 
            DO IJ=4,MFRQ 
               PHOT(IJ,INDX)=0. 
            END DO 
         END DO 
         FREQS: DO IJ=4,MFRQ 
            IF(WLAM(IJ).LE.WPHT1) THEN 
               A1=(WPHT1-WLAM(IJ))/DW 
               A2=(WLAM(IJ)-WPHT0)/DW 
               DO I=1,NPHT1 
                  INDX=INDEX(IFIL,I) 
                  PHOT(IJ,INDX)=(A1*PHT0(I)+A2*PHT1(I))*1.E-18 
               END DO 
             ELSE 
               WPHT0=WPHT1 
               DO I=1,NPHT1 
                  PHT0(I)=PHT1(I) 
               END DO 
               IFSML=0 
               READPHT: DO 
                  READ(IPHT1,*,IOSTAT=IOS2) WPHT1,(PHT1(I),I=1,NPHT1) 
                  IF(IOS2.EQ.0) EXIT READPHT 
                  IF(WPHT1.LT.WLAM(IJ)) THEN 
                     IFSML=1 
                     CYCLE READPHT 
                  END IF 
                  IF(IFSML.EQ.1) THEN 
                     BACKSPACE(IPHT1) 
                     BACKSPACE(IPHT1) 
                     READ(IPHT1,*,IOSTAT=IOS5) WPHT0,(PHT0(I),I=1,NPHT1) 
                     IF(IOS5.EQ.0) THEN 
                        READ(IPHT1,*,IOSTAT=IOS6)                         & 
     &                       WPHT1,(PHT1(I),I=1,NPHT1) 
                     END IF 
                  END IF 
               END DO READPHT 
               DW=WPHT1-WPHT0 
               A1=(WPHT1-WLAM(IJ))/DW 
               A2=(WLAM(IJ)-WPHT0)/DW 
               DO I=1,NPHT1 
                  INDX=INDEX(IFIL,I) 
                  PHOT(IJ,INDX)=(A1*PHT0(I)+A2*PHT1(I))*1.E-18 
               END DO 
            END IF 
            CYCLE FREQS 
            IEND(IFIL)=1 
         END DO FREQS 
         IF(IEND(IFIL).EQ.2) THEN 
            DO I=1,NPHT1 
               INDX=INDEX(IFIL,I) 
               PHOT(IJ,INDX)=0. 
            END DO 
         END IF 
         PHOT(1,INDX)=PHOT(3,INDX) 
         PHOT(2,INDX)=PHOT(MFRQ,INDX) 
         CYCLE READFILES 
         IEND(IFIL)=2 
         DO IJ=1,MFREQ 
            DO I=1,NELEM(IFIL) 
               INDX=INDEX(IFIL,I) 
               PHOT(IJ,INDX)=0. 
            END DO 
         END DO 
      END DO READFILES 
      RETURN 
      END SUBROUTINE READPH 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE INILIN 
!     ================= 
! 
!     read in the input line list, 
!     selection of lines that may contribute, 
!     set up auxiliary fields containing line parameters, 
! 
!     Input of line data - unit 19: 
! 
!     For each line, one (or two) records, containing: 
! 
!    ALAM    - wavelength (in nm) 
!    ANUM    - code of the element and ion (as in Kurucz-Peytremann) 
!              (eg. 2.00 = HeI; 26.00 = FeI; 26.01 = FeII; 6.03 = C IV) 
!    GF      - log gf 
!    EXCL    - excitation potential of the lower level (in cm*-1) 
!    QL      - the J quantum number of the lower level 
!    EXCU    - excitation potential of the upper level (in cm*-1) 
!    QU      - the J quantum number of the upper level 
!    AGAM    = 0. - radiation damping taken classical 
!            > 0. - the value of Gamma(rad) 
! 
!     There are now two possibilities, called NEW and OLD, of the next 
!     parameters: 
!     a) NEW, next parameters are: 
!    GS      = 0. - Stark broadening taken classical 
!            > 0. - value of log gamma(Stark) 
!    GW      = 0. - Van der Waals broadening taken classical 
!            > 0. - value of log gamma(VdW) 
!    INEXT   = 0  - no other record necessary for a given line 
!            > 0  - a second record is present, see below 
! 
!    The following parameters may or may not be present, 
!    in the same line, next to INEXT: 
!    ISQL   >= 0  - value for the spin quantum number (2S+1) of lower level 
!            < 0  - value for the spin number of the lower level unknown 
!    ILQL   >= 0  - value for the L quantum number of lower level 
!            < 0  - value for L of the lower level unknown 
!    IPQL   >= 0  - value for the parity of lower level 
!            < 0  - value for the parity of the lower level unknown 
!    ISQU   >= 0  - value for the spin quantum number (2S+1) of upper level 
!            < 0  - value for the spin number of the upper level unknown 
!    ILQU   >= 0  - value for the L quantum number of upper level 
!            < 0  - value for L of the upper level unknown 
!    IPQU   >= 0  - value for the parity of upper level 
!            < 0  - value for the parity of the upper level unknown 
!    (by default, the program finds out whether these quantum numbers 
!     are included, but the user can force the program to ignore them 
!     if present by setting INLIST=10 or larger 
! 
!    If INEXT was set to >0 then the following record includes: 
!    WGR1,WGR2,WGR3,WGR4 - Stark broadening values from Griem (in Angst) 
!                   for T=5000,10000,20000,40000 K, respectively; 
!                   and n(el)=1e16 for neutrals, =1e17 for ions. 
!    ILWN    = 0  - line taken in LTE (default) 
!            > 0  - line taken in NLTE, ILWN is then index of the 
!                   lower level 
!            =-1  - line taken in approx. NLTE, with Doppler K2 function 
!            =-2  - line taken in approx. NLTE, with Lorentz K2 function 
!    IUN     = 0  - population of the upper level in LTE (default) 
!            > 0  - index of the lower level 
!    IPRF    = 0  - Stark broadening determined by GS 
!            < 0  - Stark broadening determined by WGR1 - WGR4 
!            > 0  - index for a special evaluation of the Stark 
!                   broadening (in the present version inly for He I - 
!                   see procedure GAMHE) 
!      b) OLD, next parameters are 
!     IPRF,ILWN,IUN - the same meaning as above 
!     next record with WGR1-WGR4 - again the same meaning as above 
!     (this record is automatically read if IPRF<0 
! 
!     The only differences between NEW and OLD is the occurence of 
!     GS and GW in NEW, and slightly different format of reading. 
! 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use lindat 
      implicit real(dp) (a-h,o-z),logical (l) 
! 
      REAL(DP), PARAMETER  ::                                             & 
     &           C1     = 2.3025851,                                      & 
     &           C2     = 4.2014672,                                      & 
     &           C3     = 1.4387886,                                      & 
     &           CNM    = 2.997925e17,                                    & 
     &           ANUMIN = 1.9,                                            & 
     &           ANUMAX = 99.31,                                          & 
     &           AHE2   = 2.01,                                           & 
     &           EXT0   = 3.17,                                           & 
     &           UN     = 1.0,                                            & 
     &           TEN    = 10.,                                            & 
     &           HUND   = 1.e2,                                           & 
     &           TENM4  = 1.e-4,                                          & 
     &           TENM8  = 1.e-8,                                          & 
     &           OP4    = 0.4,                                            & 
     &           AGR0=2.4734E-22,                                         & 
     &           XEH=13.595, XET=8067.6, XNF=25.,                         & 
     &           R02=2.5, R12=45., VW0=4.5E-9,                            & 
     &           ENHE1=198310.76, ENHE2=438908.85 
      CHARACTER(LEN=1000) :: CADENA 
 
      DATA INLSET /0/ 
! 
      if(ibin(0).eq.0) then 
         open(unit=19,file=amlist(0),status='old') 
       else 
         open(unit=19,file=amlist(0),form='unformatted',status='old') 
      end if 
      if(imode.lt.-2) then 
         call inilin_grid 
         return 
      end if 
! 
      if(ndstep.eq.0) then 
         write(6,"(/' lines are rejected based on opacities at the',      & 
     & ' standard depth:'/                                                & 
     & ' ID =',i4,'  T = ',f10.1,',   DENS = ',1pe10.3/)")                & 
     &   idstd,temp(idstd),dens(idstd) 
       else 
         write(6,"(/' lines are rejected based on opacities',             & 
     &   ' at depths:'/)") 
         do id=1,nd,ndstep 
            write(6,"(' ID =',i4,'  T = ',f10.1,',  DENS = ',1pe10.3/)")  & 
     &      id,temp(id),dens(id) 
         end do 
      end if 
! 
      IL=0 
      INNLT0=0 
      IGRIE0=0 
      IF(NXTSET.EQ.1) THEN 
          ALAM0=ALM00 
          ALAST=ALST00 
          FRLAST=CNM/ALAST 
          NXTSET=0 
          REWIND 19 
      END IF 
      ALAM00=ALAM0 
      ALAST=CNM/FRLAST 
      ALAST0=ALAST 
      DOPSTD=1.E7/ALAM0*DSTD 
      DOPLAM=ALAM0*ALAM0/CNM*DOPSTD 
      AVAB=ABSTD(IDSTD)*RELOP 
      ASTD=1.0 
!     IF(GRAV.GT.6.) ASTD=0.1 
      CUTOFF=CUTOF0 
      ALAST=CNM/FRLAST 
      IF(INLTE.GE.1.AND.INLSET.EQ.0) THEN 
         CALL NLTSET(0,IL,IAT,ION,ALAM0,EXCL,EXCU,QL,QU,                  & 
     &         ISQL,ILQL,IPQL,ISQU,ILQU,IPQU,IEVEN,INNLT0,ILMATCH) 
         INLSET=1 
         ILMATCH=0 
         ILSEARCH=0 
         ILFOUND=0 
         ILFAIL=0 
         ILMULT=0 
      END IF 
! 
! 
!     Check whether any ion needs to compare quantum number limits 
! 
      MAXILIMITS=0 
      DO I=1,NION 
        IF (ILIMITS(I).EQ.1) MAXILIMITS=1 
      END DO 
      IF (MAXILIMITS.EQ.0.and.inlist.gt.0) INLIST=20 
! 
!     If INLIST=0 or 10, the program checks for the number of words 
!     present in the first line of the file to determine if quantum 
!     numbers are included. If  INLINST=11, they will be ignored anyway 
 
      IADQN=0 
      IF(ibin(0).eq.0) then 
        CADENA=' ' 
        READ(19,'(1000a)') CADENA 
        BACKSPACE(19) 
        CALL COUNT_WORDS(CADENA,NOW) 
        IF(NOW.LT.12) THEN 
           WRITE(11,*) 'INILIN: NO quantum numbers given in linelist' 
         ELSE 
           IADQN=1 
        END IF 
        if(inlist.ge.10)                                                  & 
     &  write(11,*) 'INILIN: if present, quant. num. limits are ignored' 
      ELSE 
        read(19,iostat=ios) ALAM,ANUM,GF,EXCL,QL,EXCU,QU,AGAM,            & 
     &                 GS,GW,INEXT,ISQL,ILQL,IPQL,ISQU,ILQU,IPQU 
        if(ios.eq.0) then 
!          BACKSPACE(19) 
           IADQN=1 
         else 
           backspace(19) 
           read(19) ALAM,ANUM,GF,EXCL,QL,EXCU,QU,AGAM,                    & 
     &              GS,GW,INEXT 
           backspace(19) 
        end if 
 
        if(iadqn.eq.0)                                                    & 
     &     write(11,*) 'INILIN: no quantum numbers in binary linelist' 
        IF(INLIST.GE.10) THEN 
           write(11,*)                                                    & 
     &     'INILIN: if present, quant. num. limits are ignored' 
        END IF 
      END IF 
 
      rstd=1.e4 
      if(relop.gt.0.) rstd=1./relop 
      afac=10. 
      if(iat.gt.15.and.iat.ne.26) afac=1. 
      afac=afac*rstd*astd 
! 
!     first part of reading line list - read only lambda, and 
!     skip all lines with wavelength below ALAM0-CUTOFF 
! 
      ALAM=0. 
      IJC=2 
      DO 
         if(ibin(0).eq.0) then 
            READ(19,"(F10.4)") ALAM 
          ELSE 
            READ(19) ALAM 
         END IF 
         IF(ALAM.GE.ALAM0-CUTOFF) EXIT 
      END DO 
 
      BACKSPACE(19) 
! 
!     ================== 
!     read the line list 
!     ================== 
! 
      READLINES: DO 
         ILWN=0 
         IUN=0 
         IPRF=0 
         GS=0. 
         GW=0. 
         IF(IBIN(0).EQ.0) THEN 
            IF(IADQN.EQ.0) THEN 
               READ(19,*,IOSTAT=IOS) ALAM,ANUM,GF,EXCL,QL,EXCU,QU,        & 
     &            AGAM,GS,GW 
             ELSE 
               READ(19,*,IOSTAT=IOS) ALAM,ANUM,GF,EXCL,QL,EXCU,QU,        & 
     &            AGAM,GS,GW,INEXT,ISQL,ILQL,IPQL,ISQU,ILQU,IPQU 
            END IF 
          ELSE 
            IF(IADQN.EQ.0) THEN 
               READ(19,IOSTAT=IOS) ALAM,ANUM,GF,EXCL,QL,EXCU,QU,          & 
     &            AGAM,GS,GW 
             ELSE 
               READ(19,IOSTAT=IOS) ALAM,ANUM,GF,EXCL,QL,EXCU,QU,          & 
     &            AGAM,GS,GW,INEXT,ISQL,ILQL,IPQL,ISQU,ILQU,IPQU 
            END IF 
         END IF 
!         write(*,*) 'alam,anum,gf,ios',alam,anum,gf,ios 
! 
!        test for the error in or end of line list 
! 
         IF(IOS.GT.0) THEN 
            CYCLE READLINES 
          ELSE IF(IOS.LT.0) THEN 
            EXIT READLINES 
          ELSE 
! 
!           line is in inside the computational range 
! 
            IF(INLIST.GE.10) THEN 
               IF(ISPICK.EQ.0) THEN 
                 ISQL=-1 
                 ISQU=-1 
               END IF 
               IF(ILPICK.EQ.0) THEN 
                 ILQL=-1 
                 ILQU=-1 
               END IF 
               IF(IPPICK.EQ.0) THEN 
                 IPQL=-1 
                 IPQU=-1 
               END IF 
 
            END IF 
!            write(*,*) 'inlist',inlist,ispick,ilpick,ippick 
! 
!           change wavelength to vacuum for lambda > 2000 
! 
            if(alam.gt.200..and.vaclim.gt.2000.) then 
               wl0=alam*10. 
               ALM=1.E8/(WL0*WL0) 
               XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM) 
               WL0=WL0*(XN1*1.e-6+UN) 
               alam=wl0*0.1 
            END IF 
! 
!           first selection : for a given interval a atomic number 
! 
            IF(ALAM.GT.ALAST+CUTOFF) EXIT READLINES 
            IF(ANUM.LT.ANUMIN.OR.ANUM.GT.ANUMAX) CYCLE READLINES 
            IF(ABS(ANUM-AHE2).LT.TENM4.AND.IFHE2.GT.0) CYCLE READLINES 
! 
!           second selection : for line strenghts 
! 
            FR0=CNM/ALAM 
            IAT=INT(ANUM) 
            FRA=(ANUM-FLOAT(IAT)+TENM4)*HUND 
            ION=INT(FRA)+1 
            IF(ION.GT.IONIZ(IAT)) CYCLE READLINES 
            IEVEN=1 
            EXCL=ABS(EXCL) 
            EXCU=ABS(EXCU) 
            IF(EXCL.GT.EXCU) THEN 
               FRA=EXCL 
               EXCL=EXCU 
               EXCU=FRA 
               FRA=QL 
               QL=QU 
               QU=FRA 
               IEVEN=0 
               IFRA=ISQL 
               ISQL=ISQU 
               ISQU=IFRA 
 
               IFRA=ILQL 
               ILQL=ILQU 
               ILQU=IFRA 
 
               IFRA=IPQL 
               IPQL=IPQU 
               IPQU=IFRA 
            END IF 
            GFP=C1*GF-C2 
            EPP=C3*EXCL 
            DOPSTD=1.E7/ALAM*DSTD 
            DOPLAM=ALAM*ALAM/CNM*DOPSTD 
! 
!     ************************************** 
!     rejecting weak lines 
!     ************************************* 
! 
            REJ: if(ndstep.eq.0.and.ifwin.eq.0) then 
! 
!              old procedure for rejecting lines 
! 
               GX=GFP-EPP/TSTD 
               AB0=0. 
!              write(*,*) 'rej',iat,ion,rrr(idstd,ion,iat) 
               if(gx.gt.-50)                                              & 
     &         AB0=EXP(GFP-EPP/TSTD)*RRR(IDSTD,ION,IAT)/DOPSTD/AVAB 
               IF(AB0.LT.UN) CYCLE READLINES 
! 
             else 
! 
!              new procedure for rejecting lines 
! 
               do ijcn=ijc,nfreqc 
                  if(fr0.ge.freqc(ijcn)) exit 
               end do 
               ijc=ijcn 
               if(ijc.gt.nfreqc) ijc=nfreqc 
               tkm=1.65e8/amas(iat) 
               DP0=3.33564E-11*FR0 
               do id=1,nd,ndstep 
                  td=temp(id) 
                  gx=gfp-epp/td 
                  ab0=0. 
                  if(gx.gt.-50) then 
                     dops=dp0*sqrt(tkm/td+vturb(id)) 
                     AB0=EXP(gx)*RRR(ID,ION,IAT)/(DOPS*abstdw(ijc,id)*    & 
     &               relop) 
                  end if 
                  if(ab0.ge.un) exit REJ 
               end do 
               cycle readlines 
            end if REJ 
! 
!           truncate line list if there are more lines than maximum 
!           alowable (given by MLIN0) 
! 
            IL=IL+1 
            IF(IL.GT.MLIN0) THEN 
               WRITE(6,"(' **** MORE LINES THAN MLIN0;  ',                & 
     &         'LINE LIST TRUNCATED '/                                    & 
     &'        AT LAMBDA',F15.4,'  NM'/)") ALAM 
               IL=MLIN0 
               ALAST=CNM/FREQ0(IL)-CUTOFF 
               FRLAST=CNM/ALAST 
               NXTSET=1 
               EXIT READLINES 
            END IF 
! 
!     ============================================= 
!     line is selected, set up necessary parameters 
!     ============================================= 
! 
!            store parameters for selected lines 
! 
            FREQ0(IL)=FR0 
            EXCL0(IL)=real(EPP) 
            EXCU0(IL)=real(EXCU*C3) 
            GF0(IL)=real(GFP) 
            INDAT(IL)=100*IAT+ION 
! 
!           indices for corresponding excitation temperatures 
!           of the lower and upper levels 
!           (for winds) 
! 
            if(ifwin.gt.0) then 
               IJCONT(IL)=IJC 
!C             if(excl.ge.enhe2) then 
!C                ipotl(il)=3 
!C              else if(excl.ge.enhe1) then 
!C                ipotl(il)=2 
!C              else 
!C                ipotl(il)=1 
!C             end if 
            end if 
! 
!     ****** line broadening parameters ***** 
! 
!            1) natural broadening 
! 
            IF(AGAM.GT.0.) THEN 
               GAMR0(IL)=real(EXP(C1*AGAM)) 
             ELSE 
               GAMR0(IL)=real(AGR0*FR0*FR0) 
            END IF 
! 
!           if Stark or Van der Waals broadenig assumed classical, 
!           evaluate the effective quantum number 
! 
            IF(GS.EQ.0..OR.GW.EQ.0) THEN 
               Z=FLOAT(ION) 
               XNEFF2=Z**2*(XEH/(ENEV(IAT,ION)-EXCU/XET)) 
               IF(XNEFF2.LE.0..OR.XNEFF2.GT.XNF) XNEFF2=XNF 
            END IF 
! 
!           2) Stark broadening 
! 
            IF(GS.NE.0.) THEN 
               GS0(IL)=real(EXP(C1*GS)) 
             ELSE 
               GS0(IL)=real(TENM8*XNEFF2*XNEFF2*SQRT(XNEFF2)) 
            END IF 
! 
!           3) Van der Waals broadening 
! 
            IF(GW.NE.0.) THEN 
               GW0(IL)=real(EXP(C1*GW)) 
             ELSE 
               IF(IAT.EQ.2) THEN 
                  R2=0. 
                ELSE IF(IAT.LT.21) THEN 
                  R2=R02*(XNEFF2/Z)**2 
                ELSE IF(IAT.LT.45) then 
                  R2=(R12-FLOAT(IAT))/Z 
                ELSE 
                  R2=0.5 
               END IF 
               GW0(IL)=real(VW0*R2**OP4) 
            END IF 
! 
!           evaluation of EXTIN0 - the distance (in delta freq.) where 
!           the line is supposed to contribute to the total opacity 
! 
            call profil(il,iat,idstd,agam) 
            IF(IAT.LE.2) THEN 
               EXT=SQRT(10.*AB0) 
             ELSE IF(IAT.LE.14) THEN 
               EX0=AB0*ASTD*10. 
               EXT=EXT0 
               IF(EX0.GT.TEN) EXT=SQRT(EX0) 
             ELSE 
               EX0=AB0*ASTD 
               EXT=EXT0 
               IF(EX0.GT.TEN) EXT=SQRT(EX0) 
            END IF 
            EXTIN0=EXT*DOPSTD 
            EXTIN(IL)=real(EXTIN0) 
! 
!           4) parameters for a special profile evaluation: 
! 
!           a) special He I and He II line broadening parameters 
! 
            ISPRFF=0 
            IF(IAT.LE.2) ISPRFF=ISPEC(IAT,ION,ALAM) 
            IF(IAT.EQ.2) CALL HESET(IL,ALAM,EXCL,EXCU,ION,IPRF,ILWN,IUN) 
            ISPRF(IL)=ISPRFF 
            IPRF0(IL)=IPRF 
! 
!           implied NLTE option 
! 
            if(inlte.eq.-2.or.inlte.eq.12) then 
                if(iat.le.20.and.excl.le.1000.) qu=-abs(qu) 
             else if(inlte.eq.-3) then 
                if(excl.le.1000.) qu=-abs(qu) 
             else if(inlte.eq.-4) then 
                qu=-abs(qu) 
            end if 
! 
!     ========================= 
!     NLTE lines initialization 
!     ======================== 
! 
            INDNLT(IL)=0 
            IF(QU.LT.0..OR.QL.LT.0.) THEN 
               ILWN=-1 
               QU=ABS(QU) 
               QL=ABS(QL) 
            END IF 
            IF(ILWN.LT.0.AND.INLTE.NE.0) THEN 
               INNLT0=INNLT0+1 
               INDNLT(IL)=INNLT0 
               IF(INNLT0.GT.MNLT) THEN 
                  WRITE(6,"(' **** MORE LINES IN NLTE OPTION THAN MNLT'/  & 
     &'           FOR LINES WITH LAMBDA GREATER THAN',F15.4,'  NM'/)")    & 
     &            ALAM 
                  EXIT READLINES 
               END IF 
               GI=2.*QL+UN 
               GJ=2.*QU+UN 
               CALL NLTE(IL,ILWN,IUN,GI,GJ) 
               ILOWN(IL)=ILWN 
               IUPN(IL)=IUN 
            END IF 
            IF(ILWN.GT.0.AND.INLTE.NE.0) THEN 
               INNLT0=INNLT0+1 
               INDNLT(IL)=INNLT0 
               IF(INNLT0.GT.MNLT) THEN 
                  WRITE(6,"(' **** MORE LINES IN NLTE OPTION THAN MNLT'/  & 
     &'           FOR LINES WITH LAMBDA GREATER THAN',F15.4,'  NM'/)")    & 
     &            ALAM 
                  EXIT READLINES 
               END IF 
               GI=2.*QL+UN 
               GJ=2.*QU+UN 
               CALL NLTE(IL,ILWN,IUN,GI,GJ) 
               ILOWN(IL)=ILWN 
               IUPN(IL)=IUN 
            END IF 
            IF(ILWN.EQ.0.AND.INLTE.GE.1) THEN 
               ILMATCH=-1 
               CALL NLTSET(1,IL,IAT,ION,ALAM,EXCL,EXCU,QL,QU,             & 
     &         ISQL,ILQL,IPQL,ISQU,ILQU,IPQU,IEVEN,INNLT0,ILMATCH) 
! 
!              Success accounting for nlte lines matched with 
!              quantum numbers and energy limits 
! 
!              nlte lines searched  matching energies and quantum numbers 
! 
               IF(ILMATCH.GE.0) THEN 
                   ILSEARC=ILSEARCH+1 
!                  nlte lines not found matching 
                   IF (ILMATCH.EQ.0) THEN 
                      ILFAIL=ILFAIL+1 
!                 nlte lines with multiple matches 
                    ELSE IF (ILMATCH.EQ.2) THEN 
                       ILMULT=ILMULT+1 
!                  nlte lines uniquely matched 
                     ELSE IF (ILMATCH.EQ.1) THEN 
                        ILFOUND=ILFOUND+1 
                   END IF 
               END IF 
 
               IF(INDNLT(IL).GT.0) THEN 
                  IF(INDNLT(IL).GT.MNLT) THEN 
                     WRITE(6,"(' **** MORE LINES IN NLTE OPTION',         & 
     &               'THAN MNLT'/                                         & 
     &               '   FOR LINES WITH LAMBDA GREATER THAN',             & 
     &               F15.4,'  NM'/)") ALAM 
                     EXIT READLINES 
                  END IF 
                  GI=2.*QL+UN 
                  GJ=2.*QU+UN 
                  ILWN=ILOWN(IL) 
                  IUN=IUPN(IL) 
                  IF(ILWN.EQ.IUN.AND.GI.EQ.GJ) THEN 
                     INDNLT(IL)=0 
                     ILOWN(IL)=0 
                     IUPN(IL)=0 
                   ELSE 
                     CALL NLTE(IL,ILWN,IUN,GI,GJ) 
                  END IF 
               END IF 
            END IF 
         END IF 
      END DO READLINES 
! 
      NLIN0=IL 
      NNLT=INNLT0 
      ALM1=CNM/FREQ0(1) 
      IF(ALAM0.LT.ALM1.AND.IMODE.NE.1) THEN 
         ALAM0=ALM1-4.*DOPLAM 
         IF(ALAM0.LT.ALAM00) ALAM0=ALAM00 
      END IF 
      ALM2=CNM/FREQ0(NLIN0) 
      IF(NLIN0.GT.1) ALM2=CNM/FREQ0(NLIN0-1) 
      IF(ALAST.GT.ALM2.AND.IMODE.NE.1) THEN 
         ALAST=ALM2-4.*DOPLAM 
         IF(ALAST.GT.ALAST0) ALAST=ALAST0 
         FRLAST=CNM/ALAST 
      END IF 
      IBLANK=0 
! 
      WRITE(11,*)'INILIN: NLTE matches using Energies and SLP limits --' 
      WRITE(11,*)ILSEARCH,' lines searched' 
      WRITE(11,*)ILFAIL,' lines unmatched -- set to LTE' 
      WRITE(11,*)ILMULT,' lines with multiple matches' 
      WRITE(11,*)ILFOUND,' lines uniquely matched' 
      WRITE(11,*) '--------------------------------------------------' 
! 
      WRITE(*,*)' ---------------------------------------------------' 
      WRITE(6,"(/' LINES - TOTAL        :',I10                            & 
     &       /' LINES - NLTE         :',I10/)") NLIN0,NNLT 
      RETURN 
      END SUBROUTINE INILIN 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE INILIN_grid 
!     ====================== 
! 
!     read in the input line list, 
!     selection of lines that may contribute, 
!     set up auxiliary fields containing line parameters, 
! 
!     Input of line data - unit 19: 
! 
!     For each line, one (or two) records, containing: 
! 
!    ALAM    - wavelength (in nm) 
!    ANUM    - code of the element and ion (as in Kurucz-Peytremann) 
!              (eg. 2.00 = HeI; 26.00 = FeI; 26.01 = FeII; 6.03 = C IV) 
!    GF      - log gf 
!    EXCL    - excitation potential of the lower level (in cm*-1) 
!    QL      - the J quantum number of the lower level 
!    EXCU    - excitation potential of the upper level (in cm*-1) 
!    QU      - the J quantum number of the upper level 
!    AGAM    = 0. - radiation damping taken classical 
!            > 0. - the value of Gamma(rad) 
! 
!     There are now two possibilities, called NEW and OLD, of the next 
!     parameters: 
!     a) NEW, next parameters are: 
!    GS      = 0. - Stark broadening taken classical 
!            > 0. - value of log gamma(Stark) 
!    GW      = 0. - Van der Waals broadening taken classical 
!            > 0. - value of log gamma(VdW) 
!    INEXT   = 0  - no other record necessary for a given line 
!            > 0  - next record is read, which contains: 
!    WGR1,WGR2,WGR3,WGR4 - Stark broadening values from Griem (in Angst) 
!                   for T=5000,10000,20000,40000 K, respectively; 
!                   and n(el)=1e16 for neutrals, =1e17 for ions. 
!    ILWN    = 0  - line taken in LTE (default) 
!            > 0  - line taken in NLTE, ILWN is then index of the 
!                   lower level 
!            =-1  - line taken in approx. NLTE, with Doppler K2 function 
!            =-2  - line taken in approx. NLTE, with Lorentz K2 function 
!    IUN     = 0  - population of the upper level in LTE (default) 
!            > 0  - index of the lower level 
!    IPRF    = 0  - Stark broadening determined by GS 
!            < 0  - Stark broadening determined by WGR1 - WGR4 
!            > 0  - index for a special evaluation of the Stark 
!                   broadening (in the present version inly for He I - 
!                   see procedure GAMHE) 
!      b) OLD, next parameters are 
!     IPRF,ILWN,IUN - the same meaning as above 
!     next record with WGR1-WGR4 - again the same meaning as above 
!     (this record is automatically read if IPRF<0 
! 
!     The only differences between NEW and OLD is the occurence of 
!     GS and GW in NEW, and slightly different format of reading. 
! 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use lindat 
      use optabl 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER ::                                              & 
     &           C1     = 2.3025851,                                      & 
     &           C2     = 4.2014672,                                      & 
     &           C3     = 1.4387886,                                      & 
     &           CNM    = 2.997925e17,                                    & 
     &           ANUMIN = 1.9,                                            & 
     &           ANUMAX = 99.31,                                          & 
     &           AHE2   = 2.01,                                           & 
     &           EXT0   = 3.17,                                           & 
     &           UN     = 1.0,                                            & 
     &           TEN    = 10.,                                            & 
     &           HUND   = 1.e2,                                           & 
     &           TENM4  = 1.e-4,                                          & 
     &           TENM8  = 1.e-8,                                          & 
     &           OP4    = 0.4,                                            & 
     &           AGR0=2.4734E-22,                                         & 
     &           XEH=13.595, XET=8067.6, XNF=25.,                         & 
     &           R02=2.5, R12=45., VW0=4.5E-9,                            & 
     &           bnc=1.4743e-2,hkc=4.79928e-11,                           & 
     &           ENHE1=198310.76, ENHE2=438908.85 
 
      DATA INLSET /0/ 
! 
      if(irelin.eq.0) return 
! 
      relop0=relop 
      relop=1.e-3*relop 
      if(relop.gt.1.e-4) relop=1.e-4 
      if(relop.lt.1.e-5) relop=1.e-5 
      ijcon=2 
      IL=0 
      INNLT0=0 
      IGRIE0=0 
      IF(NXTSET.EQ.1) THEN 
          ALAM0=ALM00 
          ALAST=ALST00 
          FRLAST=CNM/ALAST 
          NXTSET=0 
          REWIND 19 
      END IF 
      ALAM00=ALAM0 
      ALAST=CNM/FRLAST 
      ALAST0=ALAST 
      DOPSTD=1.E7/ALAM0*DSTD 
      DOPLAM=ALAM0*ALAM0/CNM*DOPSTD 
      AVAB=ABSTD(IDSTD)*RELOP 
      id=idstd 
      dstdid=sqrt(1.4e7*temp(idstd)) 
      ASTD=1.0 
!     IF(GRAV.GT.6.) ASTD=0.1 
      CUTOFF=CUTOF0 
      ALAST=CNM/FRLAST 
      absta=absoc(1) 
!     write(6,"(/' read line list with alam0, alast',2f10.3,1p3e11.3/)") 
!    *   alam0,alast,abstd(idstd),absta 
! 
      rstd=1.e4 
      if(relop.gt.0.) rstd=1./relop 
      afac=10. 
      if(iat.gt.15.and.iat.ne.26) afac=1. 
      afac=afac*rstd*astd 
! 
      afac=afac*rstd*astd 
      afilin=alast 
! 
!     first part of reading line list - read only lambda, and 
!     skip all lines with wavelength below ALAM0-CUTOFF 
! 
      ALAM=0. 
      DO 
         if(ibin(0).eq.0) then 
            read(19,"(F10.4)") alam 
          else 
            read(19) alam 
         end if 
         IF(ALAM.GT.ALAM0-CUTOFF) EXIT 
      END DO 
      BACKSPACE(19) 
! 
      READLINES: DO 
         IUN=0 
         IPRF=0 
         GS=0. 
         GW=0. 
         IF(IBIN(0).EQ.0) THEN 
            READ(19,*,IOSTAT=IOS) ALAM,ANUM,GF,EXCL,QL,EXCU,QU,AGAM,      & 
     &                        GS,GW 
          ELSE 
            read(19,IOSTAT=IOS) ALAM,ANUM,GF,EXCL,QL,EXCU,QU,AGAM,        & 
     &                        GS,GW 
         END IF 
         IF(IOS.GT.0) CYCLE READLINES 
         IF(IOS.LT.0) EXIT READLINES 
! 
!        change wavelength to vacuum for lambda > 2000 
! 
         if(alam.gt.200..and.vaclim.gt.2000.) then 
            wl0=alam*10. 
            ALM=1.E8/(WL0*WL0) 
            XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM) 
            WL0=WL0*(XN1*1.D-6+UN) 
                alam=wl0*0.1 
         END IF 
! 
!        first selection : for a given interval a atomic number 
! 
         IF(ALAM.GT.ALAST+CUTOFF) EXIT READLINES 
! 
!        second selection : for line strengths 
! 
         FR0=CNM/ALAM 
         if(inlist.ge.0) then 
            IAT=ifix(real(ANUM,4)) 
            FRA=(ANUM-FLOAT(IAT)+TENM4)*HUND 
            ION=INT(FRA)+1 
            IF(ION.GT.IONIZ(IAT)) CYCLE READLINES 
            IEVEN=1 
            EXCL=ABS(EXCL) 
            EXCU=ABS(EXCU) 
            IF(EXCL.GT.EXCU) THEN 
               FRA=EXCL 
               EXCL=EXCU 
               EXCU=FRA 
               FRA=QL 
               QL=QU 
               QU=FRA 
               IEVEN=0 
            END IF 
            GFP=C1*GF-C2 
            EPP=C3*EXCL 
            else 
            IF(ION.GT.IONIZ(IAT)) CYCLE READLINES 
         end if 
! 
         if(fr0.lt.freqc(ijcon)) then 
            ijcon=ijcon+1 
            absta=0.5*(absoc(ijcon)+scatc(ijcon)+                         & 
     &             absoc(ijcon-1)+scatc(ijcon-1)) 
         end if 
         abstd(id)=absta 
         dop=1.e7/alam*dstdid 
         abct=exp(gfp-epp/temp(id))*rrr(id,ion,iat) 
         abid=abct/dop/absta 
         ext=sqrt(abid*afac)*dop 
! 
         ALAX0=12. 
         if(imode.eq.-6) cycle readlines 
         if(alam.lt.afilin) then 
            if(abid.ge.relop) then 
               afilin=alam 
             else 
               if(abid.lt.relop*1.e-6) cycle readlines 
            end if 
          else if(alam.lt.9500.) then 
            if(abid.lt.relop) cycle readlines 
          else if(alam.lt.9950.) then 
            if(abid.lt.relop*1.e-9) cycle readlines 
          else 
            if(abid.lt.relop*1.e-19) cycle readlines 
         end if 
! 
         IF(ANUM.LT.ANUMIN.OR.ANUM.GT.ANUMAX) CYCLE READLINES 
         IF(ANUM.GT.ANUMAX) CYCLE READLINES 
         IF(ABS(ANUM-AHE2).LT.TENM4.AND.IFHE2.GT.0) CYCLE READLINES 
! 
         extin0=ext 
! 
!        truncate line list if there are more lines than maximum allowable 
!        (given by MLIN0) 
! 
         IL=IL+1 
         IF(IL.GT.MLIN0) THEN 
            WRITE(6,"(' ** MORE LINES THAN MLIN0, LINE LIST TRUNCATED '/  & 
     &'          AT LAMBDA',F15.4,'  NM'/)") ALAM 
            IL=MLIN0 
            ALAST=CNM/FREQ0(IL)-CUTOFF 
            FRLAST=CNM/ALAST 
            NXTSET=1 
            EXIT READLINES 
         END IF 
! 
!     ============================================= 
!     line is selected, set up necessary parameters 
!     ============================================= 
! 
!        evaluation of EXTIN0 - the distance (in delta frequency) where 
!        the line is supposed to contribute to the total opacity 
! 
!        store parameters for selected lines 
! 
         FREQ0(IL)=FR0 
         EXCL0(IL)=real(EPP,4) 
         EXCU0(IL)=real(EXCU*C3,4) 
         GF0(IL)=real(GFP,4) 
         EXTIN(IL)=real(EXTIN0,4) 
         INDAT(IL)=100*IAT+ION 
! 
!     ****** line broadening parameters ***** 
! 
!        1) natural broadening 
! 
         IF(AGAM.GT.0.) THEN 
            GAMR0(IL)=real(EXP(C1*AGAM),4) 
          ELSE 
            GAMR0(IL)=real(AGR0*FR0*FR0,4) 
         END IF 
! 
!        if Stark or Van der Waals broadening assumed classical, 
!        evaluate the effective quantum number 
! 
         IF(GS.EQ.0..OR.GW.EQ.0) THEN 
            Z=FLOAT(ION) 
            XNEFF2=Z**2*(XEH/(ENEV(IAT,ION)-EXCU/XET)) 
            IF(XNEFF2.LE.0..OR.XNEFF2.GT.XNF) XNEFF2=XNF 
         END IF 
! 
!        2) Stark broadening 
! 
         IF(GS.NE.0.) THEN 
            GS0(IL)=real(EXP(C1*GS),4) 
          ELSE 
            GS0(IL)=real(TENM8*XNEFF2*XNEFF2*SQRT(XNEFF2),4) 
         END IF 
! 
!        3) Van der Waals broadening 
! 
         IF(GW.NE.0.) THEN 
            GW0(IL)=real(EXP(C1*GW),4) 
          ELSE 
            IF(IAT.LT.21) THEN 
               R2=R02*(XNEFF2/Z)**2 
             ELSE IF(IAT.LT.45) then 
               R2=(R12-FLOAT(IAT))/Z 
             ELSE 
               R2=0.5 
            END IF 
            GW0(IL)=real(VW0*R2**OP4,4) 
         END IF 
! 
!        4) parameters for a special profile evaluation: 
! 
!        a) special He I and He II line broadening parameters 
! 
         ISPRFF=0 
         IF(IAT.LE.2) ISPRFF=ISPEC(IAT,ION,ALAM) 
         IF(IAT.EQ.2) CALL HESET(IL,ALAM,EXCL,EXCU,ION,IPRF,ILWN,IUN) 
         ISPRF(IL)=ISPRFF 
         IPRF0(IL)=IPRF 
! 
      END DO READLINES 
! 
      NLIN0=IL 
      NNLT=INNLT0 
      NGRIEM=IGRIE0 
      ALM1=CNM/FREQ0(1) 
      IF(ALAM0.LT.ALM1.AND.IMODE.NE.1) THEN 
         ALAM0=ALM1-4.*DOPLAM 
         IF(ALAM0.LT.ALAM00) ALAM0=ALAM00 
      END IF 
      ALM2=CNM/FREQ0(NLIN0) 
      IF(NLIN0.GT.1) ALM2=CNM/FREQ0(NLIN0-1) 
      IF(ALAST.GT.ALM2.AND.IMODE.NE.1) THEN 
         ALAST=ALM2-4.*DOPLAM 
         IF(ALAST.GT.ALAST0) ALAST=ALAST0 
         FRLAST=CNM/ALAST 
      END IF 
      IBLANK=0 
      relop=relop0 
! 
      WRITE(6,"(/' ATOMIC LINES        :',I10/)") NLIN0 
      RETURN 
      END SUBROUTINE INILIN_grid 
! 
! 
! ******************************************************************** 
! 
! 
 
      SUBROUTINE INIBLA 
!     ================= 
! 
!     driving procedure for treating a partial line list for the 
!     current wavelength region 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use lindat 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: DP0=3.33564E-11, DP1=1.651E8,                & 
     &                       VW1=0.42, VW2=0.45,TENM4=1.E-4,UN=1. 
! 
      IF(NLIN.EQ.0) RETURN 
      XX=FREQ(1) 
      IF(NFREQ.GE.2) XX=0.5*(FREQ(1)+FREQ(2)) 
      if(ifwin.gt.0) XX=0.5*(FREQC(1)+FREQC(NFREQC)) 
      BNU=BN*(XX*1.E-15)**3 
      HKF=HK*XX 
      if(ifwin.gt.0) XX=un 
      DO ID=1,ND 
         T=TEMP(ID) 
         ANE=ELEC(ID) 
         EXH=EXP(HKF/T) 
         EXHK(ID)=UN/EXH 
         PLAN(ID)=BNU/(EXH-UN) 
         STIM(ID)=UN-EXHK(ID) 
         if(iath.gt.0) then 
            ANP=POPUL(NKH,ID) 
            AH=DENS(ID)/WMM(ID)/YTOT(ID)-ANP 
          else 
            ah=rrr(id,1,1) 
         end if 
         AHE=RRR(ID,1,2)
         AH2=0.
         IF(IFMOL.EQ.1.AND.T.LT.TMOLIM) AH2=ANH2(ID) 
         VDWC(ID)=(AH+VW1*AHE+0.85*AH2)*(T*TENM4)**VW2 
         DO IAT=1,MATOM 
            IF(AMAS(IAT).GT.0.)                                           & 
     &      DOPA1(IAT,ID)=UN/(XX*DP0*SQRT(DP1*T/AMAS(IAT)+VTURB(ID))) 
         END DO 
      END DO 
      RETURN 
      END SUBROUTINE INIBLA 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE IDTAB 
!     ================ 
! 
!     output of selected line parameters (identification table) 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use lindat 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      CHARACTER(LEN=4) :: TYPION(30) 
      CHARACTER(LEN=4) :: APB,AP0,AP1,AP2,AP3,AP4,APR 
! 
      REAL(DP), PARAMETER :: C1=2.3025851, C2=4.2014672, C3=1.4387886 
      DATA TYPION /' I  ',' II ',' III',' IV ',' V  ',                    & 
     &             ' VI ',' VII','VIII',' IX ',' X  ',                    & 
     &             ' XI ',' XII','XIII',' XIV',' XV ',                    & 
     &             ' XVI','XVII',' 18 ',' XIX',' XX ',                    & 
     &             ' XXI','XXII',' 23 ','XXIV','XXV ',                    & 
     &             'XXVI',' 27 ',' 28 ','XXIX',' XXX'/ 
      DATA APB,AP0,AP1,AP2,AP3,AP4 /'    ','   .','   *','  **',' ***',   & 
     &                              '****'/ 
! 
      write(12,*) 'nlin',nlin 
      IF(NLIN.EQ.0.OR.IPRIN.LE.-2) RETURN 
! 
      ALM0=2.997925e18/FREQ(1) 
      ALM1=2.997925e18/FREQ(2) 
      if(ifwin.gt.0) ALM0=2.997925e18/FRQOBS(1) 
      if(ifwin.gt.0) ALM1=2.997925e18/FRQOBS(NFREQ) 
!     if(iprin.ge.2) then 
!     IF(IMODE.GE.0.OR.(IMODE.EQ.-1.AND.IBLANK.EQ.1)) 
!        WRITE(6,"(/1H ,13X, 
!    * 'LAMBDA    ATOM    LOG GF       ELO    LINE/CONT',2X, 
!    * 'EQ.WIDTH'/)") 
!     end if 
! 
      DO IL0=1,NLIN 
         IL=INDLIN(IL0) 
         ALAM=2.997925e18/FREQ0(IL) 
         ID=IDSTD 
         IJCN=IJCNTR(IL0) 
         ID0=0 
         IF(IJCN.GE.1.AND.IJCN.LE.NFREQS) ID0=IREFD(IJCN) 
         IF(ID0.GT.0.and.id0.lt.nd) ID=ID0 
         IAT=INDAT(IL)/100 
         ION=MOD(INDAT(IL),100) 
         CALL PROFIL(IL,IAT,ID,AGAM) 
         ABCNT=EXP(GF0(IL)-EXCL0(IL)/TEMP(ID))*RRR(ID,ION,IAT)*           & 
     &         STIM(ID) 
         absta=min(ch(1,idstd),ch(2,idstd)) 
         if(ifwin.le.0) then 
            DOP1=DOPA1(IAT,ID) 
            str0=abcnt*dop1/absta 
          else 
            DOP1=DOPA1(IAT,ID)/FREQ0(IL) 
            STR0=ABCNT*DOP1/ABSTDW(IJCONT(IL),ID) 
         end if 
         GF=(GF0(IL)+C2)/C1 
         EXCL=EXCL0(IL)/C3 
         IF(STR0.LE.1.2) THEN 
            WW1=0.886*STR0*(1.-STR0*(0.707-STR0*0.577)) 
          ELSE 
            WW1=SQRT(LOG(STR0)) 
         END IF 
         IF(STR0.GT.55.) THEN 
            WW2=0.5*SQRT(3.14*AGAM*STR0) 
            IF(WW2.GT.WW1) WW1=WW2 
         END IF 
         EQW=ALAM/FREQ0(IL)*1.E3/DOP1*WW1 
         STR=EQW*10. 
         APR=APB 
         IF(STR.GE.1.E0.AND.STR.LT.1.E1) APR=AP0 
         IF(STR.GE.1.E1.AND.STR.LT.1.E2) APR=AP1 
         IF(STR.GE.1.E2.AND.STR.LT.1.E3) APR=AP2 
         IF(STR.GE.1.E3.AND.STR.LT.1.E4) APR=AP3 
         IF(STR.GE.1.E4) APR=AP4 
         if(alam.ge.alm0.and.alam.lt.alm1) then 
         ill=ilown(il) 
         ilu=iupn(il) 
         if(ill.gt.0) ill=ill-nfirst(iel(ill))+1 
         if(ilu.gt.0) ilu=ilu-nfirst(iel(ilu))+1 
 
         WRITE(12,"(F11.3,2X,A4,A4,F7.2,F12.3,1PE11.2,0PF8.1,1X,          & 
     &       A4,3i4)") ALAM,TYPAT(IAT),TYPION(ION),GF,EXCL,               & 
     &                 STR0,EQW,APR,ill,ilu,id 
         end if 
      END DO 
! 
      RETURN 
      END SUBROUTINE IDTAB 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE INIBLH 
!     ================= 
! 
!     output information about hydrogen lines 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use lindat 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      CHARACTER(LEN=4) :: TYPION(30) 
      CHARACTER(LEN=4) ::  APB,AP0,AP1,AP2,AP3,AP4,APR 
! 
      REAL(DP), PARAMETER :: C1=2.3025851, C2=4.2014672, C3=1.4387886,    & 
     &                       DP0=3.33564E-11, DP1=1.651E8,                & 
     &                       VW1=0.42, VW2=0.45,TENM4=1.E-4,UN=1. 
 
      DATA TYPION /' I  ',' II ',' III',' IV ',' V  ',                    & 
     &             ' VI ',' VII','VIII',' IX ',' X  ',                    & 
     &             ' XI ',' XII','XIII',' XIV',' XV ',                    & 
     &             ' XVI','XVII',' 18 ',' XIX',' XX ',                    & 
     &             ' XXI','XXII',' 23 ','XXIV','XXV ',                    & 
     &             'XXVI',' 27 ',' 28 ','XXIX',' XXX'/ 
      DATA APB,AP0,AP1,AP2,AP3,AP4 /'    ','   .','   *','  **',' ***',   & 
     &                              '****'/ 
! 
      IF(IPRIN.LE.-2.OR.IHYL.LT.0) RETURN 
      ALM0=2.997925e18/FREQ(1) 
      ALM1=2.997925e18/FREQ(2) 
      XX=FREQ(1) 
      IF(NFREQ.GE.2) XX=0.5*(FREQ(1)+FREQ(2)) 
      BNU=BN*(XX*1.E-15)**3 
      HKF=HK*XX 
! 
      IAT=1 
      ION=1 
      IZZ=1 
      ID=IDSTD 
      T=TEMP(ID) 
      ANE=ELEC(ID) 
      EXH=EXP(HKF/T) 
      EXHK(ID)=UN/EXH 
      PLAN(ID)=BNU/(EXH-UN) 
      STIM(ID)=UN-EXHK(ID) 
      DOPA1(IAT,ID)=UN/(XX*DP0*SQRT(DP1*T/AMAS(IAT)+VTURB(ID))) 
      ISERL=ILOWH 
      ISERU=ILOWH 
      IF(alm0.GT.17000..AND.alm1.LT.21000.) THEN 
         ISERL=3 
         ISERU=4 
       ELSE IF(alm0.GT.22700.) THEN 
         ISERL=4 
         ISERU=5 
         IF(alm0.GT.32800.) ISERU=6 
         IF(alm0.GT.44660.) ISERU=7 
      END IF 
! 
      DO I=ISERL,ISERU 
         II=I*I 
         XII=UN/II 
         M1=M10 
         IF(I.LT.ILOWH) M1=ILOWH-1 
         M2=M1+1 
         IF(M1.LT.I+1) M1=I+1 
         M1=M1-1 
         M2=M20+3 
         IF(M1.LT.I+1) M1=I+1 
         if(grav.gt.3.) then 
            m2=m2+5 
            m1=m1-3 
            if(m1.gt.i+6) m1=m1-3 
         end if 
         if(grav.gt.6.) then 
            m2=m2+2 
            m1=m1-1 
            if(m1.gt.i+6) m1=m1-1 
         end if 
         IF(M1.LT.I+1) M1=I+1 
         IF(M2.GT.20) M2=20 
         ILINH=0 
         DO J=M2,M1,-1 
            CALL STARK0(I,J,izz,XKIJ,WL0,FIJ,FIJ0) 
            ALAM=WL0 
            if(alam.ge.alm0.and.alam.lt.alm1) then 
               ILINH=ILINH+1 
               GH=2.*II 
               GF=LOG10(FIJ*GH) 
               EXCL=109679.*(1.-XII) 
               EXCL0H=EXCL*C3 
               GF0H=GF*C1-C2 
               ABCNT=EXP(GF0H-EXCL0H/TEMP(ID))*RRR(ID,ION,IAT)*           & 
     &                DOPA1(IAT,ID)*STIM(ID) 
               STR0=ABCNT/ABSTD(ID) 
               IF(STR0.LE.1.2) THEN 
                  WW1=0.886*STR0*(1.-STR0*(0.707-STR0*0.577)) 
                ELSE 
                  WW1=SQRT(LOG(STR0)) 
               END IF 
               IF(STR0.GT.55.) THEN 
                  agam=0.01 
                  WW2=0.5*SQRT(3.14*AGAM*STR0) 
                  IF(WW2.GT.WW1) WW1=WW2 
               END IF 
               EQW=ALAM*ALAM/3.E18*1.E3/DOPA1(IAT,ID)*WW1 
               STR=EQW*10. 
               APR=APB 
               IF(STR.GE.1.E0.AND.STR.LT.1.E1) APR=AP0 
               IF(STR.GE.1.E1.AND.STR.LT.1.E2) APR=AP1 
               IF(STR.GE.1.E2.AND.STR.LT.1.E3) APR=AP2 
               IF(STR.GE.1.E3.AND.STR.LT.1.E4) APR=AP3 
               IF(STR.GE.1.E4) APR=AP4 
!              if(iprin.ge.2) 
!    *         WRITE(6,"(F10.3,2X,2A4,F7.2,F12.3,1PE11.2,0PF8.1, 
!    *           1X,A4,2i3)") 
!    *           ALAM,TYPAT(IAT),TYPION(ION),GF,EXCL, STR0,EQW,APR,i,j 
               WRITE(14,"(F10.3,2X,2A4,F7.2,F12.3,1PE11.2,0PF8.1,         & 
     &           1X,A4,2i3)")                                             & 
     &           ALAM,TYPAT(IAT),TYPION(ION),GF,EXCL, STR0,EQW,APR,i,j 
            end if 
         END DO 
      END DO 
! 
      RETURN 
      END SUBROUTINE INIBLH 
! 
! ******************************************************************** 
! 
! 
 
      SUBROUTINE NLTSET(MODE,IL,IAT,ION,ALAMO,EXCL,EXCU,QL,QU,            & 
     &         ISQL,ILQL,IPQL,ISQU,ILQU,IPQU,IEVEN,INNLT0,ILMATCH) 
!     =============================================================== 
! 
!     NLTE option -  automatic assignement of level indices 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use lindat 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL*8, ALLOCATABLE, SAVE :: ELIMEV(:,:),ELIMOD(:,:),               & 
     &                             ELIML(:,:),                            & 
     &                             ENREV(:,:),ENROD(:,:) 
      INTEGER,ALLOCATABLE, SAVE :: INDEV(:,:),INDOD(:,:),                 & 
     &                             INDLV(:,:),INDIO(:),                   & 
     &                             NEVEN(:),NODD(:),NLEVE(:),             & 
     &                             IATN(:),IONN(:) 
      INTEGER, SAVE             :: NODD0,NNION 
 
      INTEGER, PARAMETER  :: MNION = MIOEX,MNLEV = MLEVEL,INLLEV=13 
      REAL(DP), PARAMETER :: ECONST= 5.03411142E15 
      CHARACTER(LEN=10) :: typ 
      character(len=4)  :: typ1 
      character(len=2)  :: typ2 
      character(len=2)  :: typin(60) 
      data typin /' 1',' 2',' 3',' 4',' 5',' 6',' 7',' 8',' 9','10',      & 
     &            '11','12','13','14','15','16','17','18','19','20',      & 
     &            '21','22','23','24','25','26','27','28','29','30',      & 
     &            '31','32','33','34','35','36','37','38','39','40',      & 
     &            '41','42','43','44','45','46','47','48','49','50',      & 
     &            '51','52','53','54','55','56','57','58','59','60'/ 
! 
!     +++++++++++++++++++++++++++ 
!     MODE = 0  -  initialization 
!     +++++++++++++++++++++++++++ 
! 
      IF(MODE.EQ.0) THEN 
         NNION=0 
         READ(INLLEV,*,IOSTAT=IOS) NNION 
         IF(IOS.NE.0) RETURN 
 
         ALLOCATE (ELIMEV(MNION,MNLEV),ELIMOD(MNION,MNLEV),               & 
     &             ELIML(MNION,MNLEV),                                    & 
     &             ENREV(MNION,MNLEV),ENROD(MNION,MNLEV)) 
         ALLOCATE (INDEV(MNION,MNLEV),INDOD(MNION,MNLEV),                 & 
     &             INDLV(MNION,MNLEV),INDIO(MNION),                       & 
     &             NEVEN(MNION),NODD(MNION),NLEVE(MNION),                 & 
     &             IATN(MNION),IONN(MNION)) 
 
         IONS: DO I=1,NNION 
            READ(INLLEV,*) IATN(I),IONN(I) 
            READ(INLLEV,*) NEVEN(I) 
            IF(NEVEN(I).GT.0) THEN 
               DO J=1,NEVEN(I) 
                  READ(INLLEV,*) ELIMEV(I,J) 
               END DO 
               READ(INLLEV,*) NODD(I) 
               NODD0=NODD(I) 
               IF(NODD(I).GT.0) THEN 
                  DO J=1,NODD(I) 
                     READ(INLLEV,*) ELIMOD(I,J) 
                  END DO 
                ELSE 
                  NODD(I)=NEVEN(I) 
                  DO J=1,NODD(I) 
                     ELIMOD(I,J)=ELIMEV(I,J) 
                  END DO 
               END IF 
               INDION=0 
               DO IONEX=1,NION 
                  N0I=NFIRST(IONEX) 
                  IA=NUMAT(IATM(N0I)) 
                  IF(IA.EQ.IATN(I).AND.IZ(IONEX)-1.EQ.IONN(I))            & 
     &            INDION=IONEX 
               END DO 
               IF(INDION.LE.0) THEN 
                  call quit(' INCONSISTENCY IN UNIT 13 INPUT - NLTE') 
               END IF 
               NOFF=NFIRST(INDION)-1 
! 
               ine=1 
               ino=1 
               do ii=nfirst(indion),nlast(indion) 
                  TYP=TYPLEV(II) 
                  typ1=typ(2:5) 
                  typ2=typ(8:9) 
                  iev=0 
                  if(typ1.eq.'even') iev=1 
                  do k=1,60 
                     if(typin(k).eq.typ2) ind=k 
                  end do 
                  if(iev.eq.1) then 
                     indev(i,ine)=ii 
                     write(11,*) 'super-e ',i,ii,ine,elimev(i,ine) 
                     ine=ine+1 
                   else 
                     indod(i,ino)=ii 
                     write(11,*) 'super-o ',i,ii,ino,elimod(i,ino) 
                     ino=ino+1 
                  end if 
               end do 
            END IF 
         END DO IONS 
! 
         INDION=NNION 
         IONX: DO IONEX=1,NION 
            N0I=NFIRST(IONEX) 
            IA=NUMAT(IATM(N0I)) 
            if(isemex(ia).ge.1) CYCLE IONX 
            IONM1=IZ(IONEX)-1 
            IF(IA.EQ.1.OR.IA.EQ.2) CYCLE IONX 
            DO I=1,NNION 
               IF(IA.EQ.IATN(I).AND.IONM1.EQ.IONN(I)) CYCLE IONX 
            END DO 
            IF(NFIRST(IONEX).EQ.NLAST(IONEX)) CYCLE IONX 
            INDION=INDION+1 
            EION=ENION(NFIRST(IONEX)) 
            NLEVE(INDION)=NLAST(IONEX)-NFIRST(IONEX)+1 
            INDIO(INDION)=IONEX 
            NEVEN(INDION)=0 
            IATN(INDION)=IA 
            IONN(INDION)=IONM1 
            DELE=0. 
            DO II=NFIRST(IONEX),NLAST(IONEX) 
               I=II-NFIRST(IONEX)+1 
               E=(EION-ENION(II))*ECONST 
               IF(II.LT.NLAST(IONEX)) THEN 
                  E1=(EION-ENION(II+1))*ECONST 
                  DELE=0.5*(E1-E) 
                  ELIML(INDION,I)=E+DELE 
                ELSE 
                  IF(INLTE.GE.2) THEN 
                     ELIML(INDION,I)=E+DELE 
                   ELSE 
                     ELIML(INDION,I)=EION*ECONST 
                  END IF 
               END IF 
               INDLV(INDION,I)=II 
            END DO 
         END DO IONX 
         NNION=INDION 
 
!        Header for the table with the level assignments 
! 
         if(inlte.gt.0.and.iprin.ge.1)                                    & 
     &      WRITE(11,*)'NLTSET: IAT ION      LAMBDA     EXCL '//          & 
     &      '      EXCU    ILWN  IUN' 
 
         RETURN 
      END IF 
! 
! 
!     ++++++++++++++++++++++++++++++++++++++++++ 
!     MODE > 0  -  level indices for the line IL 
!     ++++++++++++++++++++++++++++++++++++++++++ 
! 
      IF(NNION.LE.0) RETURN 
      INION=0 
      IONM1=ION-1 
      DO I=1,NNION 
         IF(IAT.EQ.IATN(I).AND.IONM1.EQ.IONN(I)) INION=I 
      END DO 
      if(isemex(iat).ge.1) RETURN 
      IF(INION.LE.0) RETURN 
      IF(NEVEN(INION).EQ.0) IEVEN=2 
      IF(NEVEN(INION).LT.0) THEN 
         NEV1=-NEVEN(INION) 
         IF(IEVEN.EQ.1) THEN 
            ILWN=0 
            J=1 
            DO WHILE (ILWN.EQ.0 .AND. J.LE.NEV1) 
               IF(QL.EQ.ELIMEV(INION,J)) THEN 
                  DE=ENREV(INION,J) 
                  IF(EXCL.NE.0.) DE=(EXCL-DE)/EXCL 
                  IF(ABS(DE).LT.1.D-5) ILWN=INDEV(INION,J) 
               END IF 
               J=J+1 
            END DO 
            IUN=0 
            J=1 
            DO WHILE (IUN.EQ.0 .AND. J.LE.NODD(INION)) 
               IF(QU.EQ.ELIMOD(INION,J)) THEN 
                  DE=(EXCU-ENROD(INION,J))/EXCU 
                  IF(ABS(DE).LT.1.D-5) IUN=INDOD(INION,J) 
               END IF 
               J=J+1 
            END DO 
         ELSE IF(IEVEN.EQ.0) THEN 
            ILWN=0 
            J=1 
            DO WHILE (ILWN.EQ.0 .AND. J.LE.NODD(INION)) 
               IF(QL.EQ.ELIMOD(INION,J)) THEN 
                  DE=ENROD(INION,J) 
                  IF(EXCL.NE.0.) DE=(EXCL-DE)/EXCL 
                  IF(ABS(DE).LT.1.D-5) ILWN=INDOD(INION,J) 
               END IF 
               J=J+1 
            END DO 
            IUN=0 
            J=1 
            DO WHILE (IUN.EQ.0 .AND. J.LE.NEV1) 
               IF(QU.EQ.ELIMEV(INION,J)) THEN 
                  DE=(EXCU-ENREV(INION,J))/EXCU 
                  IF(ABS(DE).LT.1.D-5) IUN=INDEV(INION,J) 
               END IF 
               J=J+1 
            END DO 
         END IF 
         RETURN 
      END IF 
! 
! 
      ILWN=0 
      IUN=0 
      PARITY: IF(IEVEN.EQ.1) THEN 
         IND=0 
         JLOOP1: DO J=1,NEVEN(INION) 
            IF(EXCL.LE.ELIMEV(INION,J)) THEN 
               IND=J 
               EXIT JLOOP1 
            END IF 
         END DO JLOOP1 
         IF(IND.GT.0) THEN 
            ILWN=INDEV(INION,IND) 
! 
            IND=0 
            JLOOP2: DO J=1,NODD(INION) 
               IF(EXCU.LE.ELIMOD(INION,J)) THEN 
                  IND=J 
                  EXIT JLOOP2 
               END IF 
            END DO JLOOP2 
            IF(IND.GT.0) THEN 
               IUN=INDOD(INION,IND) 
            END IF 
         END IF 
! 
       ELSE IF(IEVEN.EQ.0) THEN 
         IND=0 
         JLOOP3: DO J=1,NODD(INION) 
            IF(EXCL.LE.ELIMOD(INION,J)) THEN 
               IND=J 
               EXIT JLOOP3 
            END IF 
         END DO JLOOP3 
         IF(IND.GT.0) THEN 
            ILWN=INDOD(INION,IND) 
! 
            IND=0 
            JLOOP4: DO J=1,NEVEN(INION) 
               IF(EXCU.LE.ELIMEV(INION,J)) THEN 
                  IND=J 
                  EXIT JLOOP4 
               END IF 
            END DO JLOOP4 
            IF(IND.GT.0) IUN=INDEV(INION,IND) 
         END IF 
! 
!        transition between levels without a distinction in parity 
! 
       ELSE 
 
        MATCH: IF (ILIMITS(INDIO(INION)).EQ.0.OR.INLIST.GE.10) THEN 
! 
!        level identification: using only energy limits 
! 
         IND=0 
         JLOOP5: DO J=1,NLEVE(INION) 
            IF(EXCL.LE.ELIML(INION,J)) THEN 
               IND=J 
               EXIT JLOOP5 
            END IF 
         END DO JLOOP5 
         IF(IND.GT.0) THEN 
            ILWN=INDLV(INION,IND) 
! 
            IND=0 
            JLOOP6: DO J=1,NLEVE(INION) 
               IF(EXCU.LE.ELIML(INION,J)) THEN 
                  IND=J 
                  EXIT JLOOP6 
               END IF 
            END DO JLOOP6 
            IUN=INDLV(INION,IND) 
         END IF 
 
        ELSE IF (ILIMITS(INDIO(INION)).EQ.1.and.inlist.lt.10) THEN 
! 
!        level identification: using energy limits and quantum numbers 
! 
 
         IND=0 
         INMATCHL=0 
         DO J=1,NLEVE(INION) 
            IF(EXCL.GE.ENION1(INDLV(INION,J)).AND.                        & 
     &         EXCL.LE.ENION2(INDLV(INION,J)).AND.                        & 
     &         ((IPQL.GE.PQUANT1(INDLV(INION,J)).AND.                     & 
     &         IPQL.LE.PQUANT2(INDLV(INION,J))).OR.                       & 
     &         (IPQL.EQ.-1))                    .AND.                     & 
     &         ((ISQL.GE.SQUANT1(INDLV(INION,J)).AND.                     & 
     &         ISQL.LE.SQUANT2(INDLV(INION,J))).OR.                       & 
     &         (ISQL.EQ.-1))                    .AND.                     & 
     &         ((ILQL.GE.LQUANT1(INDLV(INION,J)).AND.                     & 
     &         ILQL.LE.LQUANT2(INDLV(INION,J))).OR.                       & 
     &         (ILQL.EQ.-1))                                              & 
     &                                             ) THEN 
               IND=J 
               INMATCHL=INMATCHL+1 
            END IF 
         END DO 
         IF (INMATCHL.GT.1)                                               & 
     &       WRITE(11,'(A55,1X,F12.4)')                                   & 
     &       ' NLTSET: WARNING-- multiple matches for lower level of ',   & 
     &       ALAMO 
         IF (INMATCHL.LE.0) THEN 
            ILWN=0 
            IUN=0 
          ELSE 
            ILWN=INDLV(INION,IND) 
         END IF 
! 
! 
         IND=0 
         INMATCHU=0 
         DO J=1,NLEVE(INION) 
            IF(EXCU.GE.ENION1(INDLV(INION,J))   .AND.                     & 
     &         EXCU.LE.ENION2(INDLV(INION,J))   .AND.                     & 
     &         ((IPQU.GE.PQUANT1(INDLV(INION,J)).AND.                     & 
     &         IPQU.LE.PQUANT2(INDLV(INION,J))).OR.                       & 
     &         (IPQU.EQ.-1))                    .AND.                     & 
     &         ((ISQU.GE.SQUANT1(INDLV(INION,J)).AND.                     & 
     &         ISQU.LE.SQUANT2(INDLV(INION,J))).OR.                       & 
     &         (ISQU.EQ.-1))                    .AND.                     & 
     &         ((ILQU.GE.LQUANT1(INDLV(INION,J)).AND.                     & 
     &         ILQU.LE.LQUANT2(INDLV(INION,J))).OR.                       & 
     &         (ILQU.EQ.-1))                                              & 
     &                                             ) THEN 
 
               IND=J 
               INMATCHU=INMATCHU+1 
            END IF 
         END DO 
         IF (INMATCHU.GT.1)                                               & 
     &       WRITE(11,'(A55,1X,F12.4)')                                   & 
     &       ' NLTSET: WARNING-- multiple matches for upper level of ',   & 
     &       ALAMO 
         IF (INMATCHU.LE.0) THEN 
            IUN=0 
          ELSE 
            IUN=INDLV(INION,IND) 
         END IF 
 
         IF (INMATCHL.EQ.0.or.INMATCHU.EQ.0) THEN 
                ILMATCH=0 
           ELSE IF (INMATCHL.GT.1.or.INMATCHU.GT.1) THEN 
                ILMATCH=2 
           ELSE 
                ILMATCH=1 
         ENDIF 
 
        ELSE 
 
         write(11,*)('ILIMITS is neither 0 or 1') 
 
        END IF MATCH 
 
        if(inlte.gt.0.and.iprin.ge.1)                                     & 
     &     WRITE(11,'(10x,2(i2,1x),3x,3(F10.3,1x),2(i4,1x))')IAT,ION,     & 
     &              ALAMO,EXCL,EXCU,ILWN,IUN 
 
      END IF PARITY 
! 
! 
      IF(INLTE.EQ.5) THEN 
         INNLT0=INNLT0+1 
         INDNLT(IL)=INNLT0 
        ELSE IF(INLTE.EQ.4) THEN 
         IF(ILWN.GT.0.AND.IUN.GT.0) THEN 
            INDNLT(IL)=-1 
         END IF 
        ELSE IF(INLTE.EQ.3) THEN 
         IF(ILWN.GT.0) THEN 
            INDNLT(IL)=-1 
         END IF 
        ELSE 
         INDNLT(IL)=-1 
      END IF 
      BNUL(IL)=real(BN*(FREQ0(IL)*1.E-15)**3) 
      ILOWN(IL)=ILWN 
      IUPN(IL)=IUN 
      RETURN 
      END SUBROUTINE NLTSET 
! 
! ******************************************************************** 
! ******************************************************************** 
! 
! 
 
      SUBROUTINE PHTION(ID,ABSO,EMIS,FRE,NFRE) 
!     ======================================== 
! 
!     Opacity due to detailed photoionization (read from tables by 
!     routine READPH) 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use lindat 
      use photcs 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: ABSO(MFRQ),EMIS(MFRQ),PLANF(MFRQ),STIMU(MFRQ),          & 
     &            FRE(MFRQ) 
      REAL(DP), PARAMETER :: C3=1.4387886 
! 
      IF(NPHT.LE.0) RETURN 
      T=TEMP(ID) 
      DO IJ=1,NFRE 
         XX=FRE(IJ) 
         X15=XX*1.E-15 
         BNU=BN*X15*X15*X15 
         HKF=HK*XX 
         EXH=EXP(HKF/T) 
         PLANF(IJ)=BNU/(EXH-1.) 
         STIMU(IJ)=1.-1./EXH 
      END DO 
      DO I=1,NPHT 
         IF(JPHT(I).LE.0) THEN 
            IAT=int(APHT(I)) 
            X=(APHT(I)-FLOAT(IAT)+1.E-4)*1.E2 
            ION=INT(X)+1 
            POP=RRR(ID,ION,IAT)*GPHT(I)*EXP(-EPHT(I)*C3/T) 
          ELSE 
            JJ=JPHT(I) 
            POP=POPUL(JJ,ID) 
         END IF 
         DO IJ=1,NFRE 
            ABP=PHOT(IJ,I)*POP*STIMU(IJ) 
            ABSO(IJ)=ABSO(IJ)+ABP 
            EMIS(IJ)=EMIS(IJ)+ABP*PLANF(IJ) 
         END DO 
      END DO 
      RETURN 
      END SUBROUTINE PHTION 
 
! 
! ******************************************************************** 
! ******************************************************************** 
! 
      SUBROUTINE NLTE(IL,ILW,IUN,GI,GJ) 
!     =========================================== 
! 
!     Control procedure for the NLTE option 
! 
      use accura 
      use params 
      use modelp 
      use lindat 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: UN   = 1.,                                   & 
     &                       C3   = 1.4387886,                            & 
     &                       XET  = 8067.6,                               & 
     &                       XET3 = XET*C3 
! 
!     CALCULATION OF THE 
!     CENTRAL OPACITY (ABCENT) AND THE LINE SOURCE FUNCTION  (SLIN) 
! 
      if(gi.le.0..or.gj.le.0.) return 
      ILNLT=INDNLT(IL) 
      IF(ILNLT.LE.0) RETURN 
      IAT=INDAT(IL)/100 
      ION=MOD(INDAT(IL),100) 
      EGF=EXP(GF0(IL)) 
      BNU=BN*(FREQ0(IL)*1.E-15)**3 
      DP0=3.33564E-11*FREQ0(IL) 
      DP1=1.651E8/AMAS(IAT) 
      IF(ILW.GT.0) THEN 
! 
!     line is a transition between explicit levels of the 
!     input model 
! 
         NKI=NNEXT(IEL(ILW)) 
         DO ID=1,ND 
            T=TEMP(ID) 
            COR=1. 
            PP=PNLT(IAT,ION,ID) 
            IF(ILW.GT.0) THEN 
               PI=POPUL(ILW,ID)/G(ILW) 
             ELSE 
               PI=PP*EXP((ENEV(IAT,ION)*XET3-EXCL0(IL))/T) 
            END IF 
            IF(IUN.GT.0) THEN 
               PJ=POPUL(IUN,ID)/G(IUN) 
               cor=(excu0(il)-excl0(il)+                                  & 
     &            (enion(iun)-enion(ilw))/1.38054e-16)/t 
               cor=exp(cor) 
             ELSE 
               PJ=PP*EXP((ENEV(IAT,ION)*XET3-EXCU0(IL))/T) 
            END IF 
            if(pj.gt.0.) then 
               X=PI/PJ*cor 
             else 
               x=un 
            end if 
            IF(X.EQ.UN) X=EXP(4.79928E-11*FREQ0(IL)/T) 
            DOP=DP0*SQRT(DP1*T+VTURB(ID)) 
            SLIN(ILNLT,ID)=BNU/(X-UN) 
            if(pi.gt.0.) ABCENT(ILNLT,ID)=PI*(UN-UN/X)*EGF/DOP 
         END DO 
         RETURN 
      END IF 
! 
!     Approximate NLTE for resonance lines - second order escape 
!     probablity theory form of the source function 
! 
!     Optical depth scale 
! 
      ALMIL=2.997925E17/FREQ0(IL) 
      HKF=HK*FREQ0(IL) 
      DO ID=1,ND 
         T=TEMP(ID) 
         DOP=DP0*SQRT(DP1*T+VTURB(ID)) 
         X=EXP(HKF/T) 
         ABCENT(ILNLT,ID)=EGF*EXP(-EXCL0(IL)/T)*RRR(ID,ION,IAT)/          & 
     &                    DOP*(1.-1./X) 
         ABA=ABSTD(ID)+ABCENT(ILNLT,ID)*1.77245 
         if(ifwin.gt.0)                                                   & 
     &   ABA=ABSTDW(IJCONT(IL),ID)+ABCENT(ILNLT,ID)*1.77245 
         IF(ID.EQ.1) THEN 
            ABM=ABA/DENS(1) 
            TAU=0.5*DM(1)*ABM 
          ELSE 
            AB0=ABA/DENS(ID) 
            TAU=TAU+0.5*(DM(ID)-DM(ID-1))*(AB0+ABM) 
            ABM=AB0 
         END IF 
! 
!        approximate epsilon after Kastner 
! 
         E=EPS(T,ELEC(ID),ALMIL,ION,IUN) 
         XK2=XK2DOP(TAU) 
         SLIN(ILNLT,ID)=SQRT(E/(E+(1.-E)*XK2))*BNU/(X-1.) 
      END DO 
      RETURN 
      END SUBROUTINE NLTE 
! 
! ******************************************************************** 
! 
! 
 
      SUBROUTINE LINOP(ID,ABLIN,EMLIN,AVAB) 
!     ===================================== 
! 
!     TOTAL LINE OPACITY (ABLIN)  AND EMISSIVITY (EMLIN) 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use lindat 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER ::                                              & 
     &           UN     = 1.,                                             & 
     &           EXT0   = 3.17,                                           & 
     &           TEN    = 10.,                                            & 
     &           C3     = 1.4387886,                                      & 
     &           XET    = 8067.6,                                         & 
     &           XET3   = XET*C3 
      REAL(DP):: ABLIN(MFREQ),EMLIN(MFREQ),ABLINN(MFREQ) 
! 
      DO IJ=1,NFREQ 
         ABLIN(IJ)=0. 
         ABLINN(IJ)=0. 
         EMLIN(IJ)=0. 
      END DO 
! 
      IF(NLIN.EQ.0) RETURN 
! 
!     overall loop over contributing lines 
! 
      TEM1=UN/TEMP(ID) 
      LINELOOP: DO I=1,NLIN 
         IL=INDLIN(I) 
         INNLT=INDNLT(IL) 
         IAT=INDAT(IL)/100 
         ION=MOD(INDAT(IL),100) 
         LPR=.TRUE. 
         ISP=ISPRF(IL) 
         IF(ISP.GT.1.AND.ISP.LE.5) LPR=.FALSE. 
         IF (ISP.GE.6) CYCLE LINELOOP 
         CALL PROFIL(IL,IAT,ID,AGAM) 
         DOP1=DOPA1(IAT,ID) 
         FR0=FREQ0(IL) 
         IF(INNLT.EQ.0) THEN 
            AB0=EXP(GF0(IL)-EXCL0(IL)*TEM1)*RRR(ID,ION,IAT)*              & 
     &          DOP1*STIM(ID) 
          ELSE IF(INNLT.GT.0) THEN 
            AB0=ABCENT(INNLT,ID) 
            SL0=SLIN(INNLT,ID) 
          ELSE 
            ILW=ILOWN(IL) 
            IUN=IUPN(IL) 
            COR=1. 
            PP=PNLT(IAT,ION,ID) 
            IF(ILW.GT.0) THEN 
               PI=POPUL(ILW,ID)/G(ILW) 
             ELSE 
               PI=PP*EXP((ENEV(IAT,ION)*XET3-EXCL0(IL))*TEM1) 
            END IF 
            IF(IUN.GT.0) THEN 
               PJ=POPUL(IUN,ID)/G(IUN) 
               cor=(excu0(il)-excl0(il)+                                  & 
     &             (enion(iun)-enion(ilw))/1.38054e-16)*tem1 
               cor=exp(cor) 
             ELSE 
               PJ=PP*EXP((ENEV(IAT,ION)*XET3-EXCU0(IL))*TEM1) 
            END IF 
            if(pj.gt.0.) then 
               X=PI/PJ*cor 
             else 
               x=un 
            end if 
            IF(X.EQ.UN) X=EXP(4.79928E-11*FREQ0(IL)*TEM1) 
            SL0=BNUL(IL)/(X-UN) 
            ab0=0. 
            if(pi.gt.0.) AB0=PI*(UN-UN/X)*EXP(GF0(IL))*DOP1 
         END IF 
!        if(id.le.10) write(*,"(6i5,1p4e12.3)")                           &
!    &   id,i,il,innlt,iat,ion,rrr(id,ion,iat),dop1,stim(id),ab0
         if(ab0.le.0.and.lasdel) CYCLE LINELOOP 
! 
!        set up limiting frequencies where the line I is supposed to 
!        contribute to the opacity 
! 
         EX0=AB0/AVAB*AGAM 
         EXT=EXT0 
         IF(EX0.GT.TEN) EXT=SQRT(EX0) 
         EXT=EXT/DOP1 
         XIJEXT=DFRCON*EXT+1.5 
!        IJ1=MAX(IJCNTR(I)-IJEXT,3) 
!        IJ2=MIN(IJCNTR(I)+IJEXT,NFREQS) 
         IJ1=int(MAX(float(IJCNTR(I))-XIJEXT,3.)) 
         IJ2=int(MIN(float(IJCNTR(I))+XIJEXT,float(NFREQS))) 
         IF(IJ1.GE.NFREQ.OR.IJ2.LE.2) CYCLE LINELOOP 
! 
         LTELINE: IF(INNLT.EQ.0) THEN 
! 
!        ********* 
!        LTE lines 
!        ********* 
! 
!        first, He I lines with occuopation probabilities 
! 
            if(iat.eq.2.and.ion.eq.1) then 
               ELO=EXCL0(IL)/C3 
               EUP=EXCU0(IL)/C3 
               CALL HE1LIN(ID,FR0,ELO,EUP,WLO,WUP) 
               WLU=WLO*WUP 
               DO IJ=IJ1,IJ2 
                  FR=FREQ(IJ) 
                  XF=ABS(FR-FR0)*DOP1 
                  IF(.NOT.LPR) THEN 
                     ABL=AB0*WLU*PHE1(ID,FR,ISP-1) 
                   ELSE 
                     ABL=AB0*WLU*VOIGTK(AGAM,XF) 
                  END IF 
                  ABLINN(IJ)=ABLINN(IJ)+ABL 
                  EMLIN(IJ)=EMLIN(IJ)+ABL*PLAN(ID) 
               END DO 
             ELSE 
! 
!            other species 
! 
               DO IJ=IJ1,IJ2 
                  XF=ABS(FREQ(IJ)-FR0)*DOP1 
                  ABLIN(IJ)=ABLIN(IJ)+AB0*VOIGTK(AGAM,XF) 
!            if(id.le.10.and.ij.eq.3)                                     &
!    &       write(6,"(6i5,1p7e12.3)") id,i,il,innlt,iat,ion,             &
!    &       freq(ij),fr0,dop1,xf,agam,voigtk(agam,xf),ab0
               END DO 
            END IF 
! 
!        ********** 
!        NLTE LINES 
!        ********** 
! 
          ELSE 
            IF(LPR) THEN 
! 
               DO IJ=IJ1,IJ2 
                  XF=ABS(FREQ(IJ)-FR0)*DOP1 
                  ABL=AB0*VOIGTK(AGAM,XF) 
                  ABLINN(IJ)=ABLINN(IJ)+ABL 
                  EMLIN(IJ)=EMLIN(IJ)+ABL*SL0 
               END DO 
! 
!        again, special expressions for 4 selected He I lines 
! 
             ELSE 
               DO IJ=3,NFREQ 
                  FR=FREQ(IJ) 
                  ABL=AB0*PHE1(ID,FR,ISP-1) 
                  ABLINN(IJ)=ABLINN(IJ)+ABL 
                  EMLIN(IJ)=EMLIN(IJ)+ABL*SL0 
               END DO 
            END IF 
         END IF LTELINE 
      END DO LINELOOP 
! 
      DO IJ=3,NFREQ 
         EMLIN(IJ)=EMLIN(IJ)+ABLIN(IJ)*PLAN(ID) 
         ABLIN(IJ)=ABLIN(IJ)+ABLINN(IJ) 
      END DO 
! 
!     special routine for selected He II lines 
! 
      IF(NSP.EQ.0) RETURN 
      DO IS=1,NSP 
         ISP=ISP0(IS) 
         IF(ISP.GE.6.AND.ISP.LE.24) CALL PHE2(ISP,ID,ABLIN,EMLIN) 
      END DO 
! 
      RETURN 
      END SUBROUTINE LINOP 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE HE1LIN(ID,FR0,ELO,EUP,WLO,WUP) 
!     ========================================= 
! 
!     compute the occupation probabilities for He I levels 
!     corresponding to a line from the atomic line list. 
!     the occupation prtobabilities WLO and WUP are assumed hydrogenic 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use lindat 
      implicit real(dp) (a-h,o-z),logical (l) 
 
!     real*4 elo,eup 
      REAL(DP), PARAMETER :: ENH=109678.7568, ENHE=198305.,               & 
     &                       C3= 1.4387886 
! 
      ALM=2.997925E17/FR0 
! 
!     first, set the vacuum wavelength for the line, if needed 
! 
      if(alm.gt.200..and.vaclim.gt.2000.) then 
         wl0=alm*10. 
         ALA=1.E8/(WL0*WL0) 
         XN1=64.328+29498.1/(146.-ALA)+255.4/(41.-ALA) 
         WL0=WL0*(XN1*1.e-6+1.) 
         ala=wl0*0.1 
         alm=ala 
      END IF 
!     eup=min(1.e7/alm-elo,198304.) 
!     EUP=FR0/2.997925E10+ELO 
      xnnl=max(sqrt(enh/(enhe-elo)),1.) 
      xnnu=sqrt(enh/(enhe-eup)) 
      nlo=int(xnnl+0.1) 
      nup=int(xnnu+0.1) 
      wlo=wnhint(nlo,id) 
      if(nup.gt.nlmx) then 
         wup=0.0 
         return 
      end if 
      wup=wnhint(nup,id) 
      return 
      end SUBROUTINE HE1LIN 
 
! 
! ******************************************************************** 
! 
! 
 
      SUBROUTINE LINOPW(ID,ABLIN,EMLIN) 
!     ================================= 
! 
!     TOTAL LINE OPACITY (ABLIN)  AND EMISSIVITY (EMLIN) 
!     (a variant for winds) 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use lindat 
      use wincom 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER ::                                              & 
     &           UN     = 1.,                                             & 
     &           EXT0   = 3.17,                                           & 
     &           TEN    = 10.,                                            & 
     &           C3     = 1.4387886,                                      & 
     &           XET    = 8067.6,                                         & 
     &           XET3   = XET*C3 
      REAL(DP):: ABLIN(MFREQ),EMLIN(MFREQ),ABLINN(MFREQ) 
! 
      DO IJ=1,NFREQ 
         ABLIN(IJ)=0. 
         ABLINN(IJ)=0. 
         EMLIN(IJ)=0. 
      END DO 
      wdil(id)=1. 
      plw=plan(id)*wdil(id) 
! 
      IF(NLIN.EQ.0) RETURN 
! 
!     overall loop over contributing lines 
! 
      TEM1=UN/TEMP(ID) 
      HKT=HK*TEM1 
      xx=freq(nopac)-freq(1) 
      DFRCON=NOPAC-1 
      DFRCON=-DFRCON/XX 
      IFRCON=int(DFRCON) 
      LINELOOP: DO I=1,NLIN 
         IL=INDLIN(I) 
         INNLT=INDNLT(IL) 
! 
!        rejecting lines for v > velmax 
! 
         if(ilvi(id).gt.0) then 
            if(innlt.eq.0) then 
               cycle lineloop 
             else 
               if(nltoff.ne.0) cycle lineloop 
            end if 
         end if 
! 
! 
!     frequency indices of the line centers 
! 
         if (id.eq.1) then 
 
            fr0=freq0(il) 
            XJC=3.+DFRCON*(FREQ(1)-FR0) 
            IJC=int(XJC) 
            IJCNTR(I)=IJC 
            if(ijc.gt.1.or.ijc.lt.nopac) then 
               if(fr0.lt.freq(ijc)) then 
                  ijc0=ijc 
                  dfr0=freq(ijc0)-fr0 
                  low: do 
                     ijc0=ijc0+1 
                     dfr=abs(freq(ijc0)-fr0) 
                     if(dfr.lt.dfr0) then 
                        ijc=ijc0 
                        ijc0=ijc0+1 
                        dfr0=dfr 
                        cycle low 
                      else 
                        exit low 
                     end if 
                  end do low 
                else if(fr0.gt.freq(ijc)) then 
                  ijc0=ijc 
                  dfr0=fr0-freq(ijc0) 
                  high: do 
                     ijc0=ijc0-1 
                     dfr=abs(freq(ijc0)-fr0) 
                     if(dfr.lt.dfr0) then 
                        ijc=ijc0 
                        ijc0=ijc0-1 
                        dfr0=dfr 
                        cycle high 
                      else 
                        exit high 
                     end if 
                  end do high 
               end if 
               IJCNTR(I)=IJC 
            end if 
         end if 
! 
         IAT=INDAT(IL)/100 
         ION=MOD(INDAT(IL),100) 
         FR0=FREQ0(IL) 
         LPR=.TRUE. 
         ISP=ISPRF(IL) 
         IF(ISP.GT.1.AND.ISP.LE.5) LPR=.FALSE. 
         IF (ISP.GE.6) CYCLE LINELOOP 
         CALL PROFIL(IL,IAT,ID,AGAM) 
         DOP1=DOPA1(IAT,ID)/FR0 
         FR0=FREQ0(IL) 
         IF(INNLT.EQ.0) THEN 
            if(itrad.le.0) then 
            AB0=EXP(GF0(IL)-EXCL0(IL)*TEM1)*RRR(ID,ION,IAT)*              & 
     &          DOP1*(1.-exp(-hkt*fr0)) 
            else 
!C          trl=trad(ipotl(il),id) 
            trl=temp(id) 
            xx=exp(-hkt*fr0) 
            AB0=EXP(GF0(IL)-EXCL0(IL)/trl)*RRR(ID,ION,IAT)*               & 
     &          DOP1*(1.-xx) 
            if(excl0(il).gt.2000.) ab0=ab0*wdil(id) 
            pla=1.4743e-2*(fr0*1.e-15)**3*xx/(1.-xx) 
            sl0=pla*wdil(id) 
            end if 
          ELSE IF(INNLT.GT.0) THEN 
            AB0=ABCENT(INNLT,ID) 
            SL0=SLIN(INNLT,ID) 
          ELSE 
            ILW=ILOWN(IL) 
            IUN=IUPN(IL) 
            COR=1. 
            PP=PNLT(IAT,ION,ID) 
            IF(ILW.GT.0) THEN 
               PI=POPUL(ILW,ID)/G(ILW) 
             ELSE 
               PI=PP*EXP((ENEV(IAT,ION)*XET3-EXCL0(IL))*TEM1) 
            END IF 
            IF(IUN.GT.0) THEN 
               PJ=POPUL(IUN,ID)/G(IUN) 
               cor=(excu0(il)-excl0(il)+                                  & 
     &             (enion(iun)-enion(ilw))/1.38054e-16)*tem1 
               cor=exp(cor) 
             ELSE 
               PJ=PP*EXP((ENEV(IAT,ION)*XET3-EXCU0(IL))*TEM1) 
            END IF 
            if(pj.gt.0.) then 
               X=PI/PJ*cor 
             else 
               x=un 
            end if 
            IF(X.EQ.UN) X=EXP(4.79928E-11*FREQ0(IL)*TEM1) 
            SL0=BNUL(IL)/(X-UN) 
            ab0=0. 
            if(pi.gt.0.) AB0=PI*(UN-UN/X)*EXP(GF0(IL))*DOP1 
         END IF 
         if(ab0.le.0.and.lasdel) cycle lineloop 
! 
!        set up limiting frequencies where the line I is supposed to 
!        contribute to the opacity 
! 
!        if(ifwin.le.0) then 
            avabw=abstdw(ijcont(il),id)*relop 
            EX0=AB0/AVABw*AGAM 
            EXT=EXT0 
            IF(EX0.GT.TEN) EXT=SQRT(EX0) 
            EXT=EXT/DOP1 
            IJEXT=int((DFRCON*EXT)+1.5) 
            IJ1=MAX(IJCNTR(I)-IJEXT,1) 
            IJ2=MIN(IJCNTR(I)+IJEXT,NFREQ) 
            IF(IJ1.GE.NFREQ.OR.IJ2.LE.2) CYCLE LINELOOP 
!        end if 
! 
         IF(INNLT.EQ.0.and.itrad.le.0) THEN 
! 
!        ********* 
!        LTE lines 
!        ********* 
! 
         IF(LPR) THEN 
! 
            DO IJ=IJ1,IJ2 
               XF=ABS(FREQ(IJ)-FR0)*DOP1 
               ABLIN(IJ)=ABLIN(IJ)+AB0*VOIGTK(AGAM,XF) 
            END DO 
! 
!        special expressions for 4 selected He I lines 
! 
          ELSE 
            DO IJ=1,NFREQ 
               FR=FREQ(IJ) 
               ABL=AB0*PHE1(ID,FR,ISP-1) 
               ABLIN(IJ)=ABLIN(IJ)+ABL 
            END DO 
         END IF 
! 
!        ********** 
!        NLTE LINES 
!        ********** 
! 
       ELSE 
         IF(LPR) THEN 
! 
            DO IJ=IJ1,IJ2 
               XF=ABS(FREQ(IJ)-FR0)*DOP1 
               ABL=AB0*VOIGTK(AGAM,XF) 
               ABLINN(IJ)=ABLINN(IJ)+ABL 
               if(ilne(id).le.0) EMLIN(IJ)=EMLIN(IJ)+ABL*SL0 
           END DO 
! 
!        again, special expressions for 4 selected He I lines 
! 
         ELSE 
            DO IJ=1,NFREQ 
               FR=FREQ(IJ) 
               ABL=AB0*PHE1(ID,FR,ISP-1) 
               ABLINN(IJ)=ABLINN(IJ)+ABL 
               if(ilne(id).le.0) EMLIN(IJ)=EMLIN(IJ)+ABL*SL0 
            END DO 
         END IF 
      END IF 
      END DO LINELOOP 
! 
      if(vel(id).le.velmax) then 
      DO IJ=1,NFREQ 
         PLA=BNUE(IJ)/(EXP(HKT*FREQ(IJ))-1.) 
         EMLIN(IJ)=EMLIN(IJ)+ABLIN(IJ)*pla*wdil(id) 
         ABLIN(IJ)=ABLIN(IJ)+ABLINN(IJ) 
      END DO 
      end if 
! 
!     special routine for selected He II lines 
! 
      IF(NSP.EQ.0) RETURN 
      DO IS=1,NSP 
         ISP=ISP0(IS) 
         IF(ISP.GE.6.AND.ISP.LE.24) CALL PHE2(ISP,ID,ABLIN,EMLIN) 
      END DO 
! 
      RETURN 
      END SUBROUTINE LINOPW 
! 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE PROFIL(IL,IAT,ID,AGAM) 
!     ================================= 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use lindat 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: WGR(4) 
      REAL(DP), PARAMETER :: PI4=7.95774715E-2 
! 
      IPRF=IPRF0(IL) 
      T=TEMP(ID) 
      ANE=ELEC(ID) 
! 
!     radiative broadening (classical) 
! 
      AGAM=GAMR0(IL) 
! 
!     Stark broadening - standard (given in the line list or classical) 
! 
      IF(IPRF.EQ.0) THEN 
         AGAM=AGAM+GS0(IL)*ANE 
! 
!     Stark broadening - special expressions for He I 
! 
       ELSE IF(IPRF.GT.0) THEN 
         ANP=POPUL(NKH,ID) 
         CALL GAMHE(IPRF,T,ANE,ANP,ID,GAM) 
         AGAM=AGAM+GAM 
! 
!     Stark broadening - Griem 
! 
       ELSE 
         DO I=1,4 
            WGR(I)=WGR0(I,IGRIEM(IL)) 
         END DO 
         FR=FREQ0(IL) 
         ION=MOD(INDAT(IL),100) 
         CALL GRIEM(ID,T,ANE,ION,FR,WGR,GAM) 
         AGAM=AGAM+GAM 
      END IF 
! 
!     Van Der Waals broadening 
! 
      AGAM=AGAM+GW0(IL)*VDWC(ID) 
! 
!     final Voigt parameter a 
! 
      DOP1=DOPA1(IAT,ID) 
      if(ifwin.gt.0) DOP1=DOP1/FREQ0(IL) 
      AGAM=AGAM*DOP1*PI4 
! 
      RETURN 
      END SUBROUTINE PROFIL 
! 
! ******************************************************************** 
! 
      SUBROUTINE GRIEM(ID,T,ANE,ION,FR,WGR,GAM) 
!     ========================================= 
! 
!     STARK DAMPING PARAMETER (GAM) CALCULATED FROM INPUT VALUES 
!     OF STARK WIDTHS FOR  T=5000, 10000, 20000, 40000 K, 
!     AND FOR  NE=1.E16 (FOR NEUTRALS)  OR  NE = 1.E17 (FOR IONS) 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: WGR(4) 
      if(t.le.0.) return 
      J=JT(ID) 
      GAM=(TI0(ID)*WGR(J)+TI1(ID)*WGR(J-1)+TI2(ID)*WGR(J-2))              & 
     &    *ANE*1.E-10*FR*1.E-10*FR*4.2E-14 
      IF(ION.GT.1) GAM=GAM*0.1 
      IF(GAM.LT.0.) GAM=0. 
      RETURN 
      END SUBROUTINE GRIEM 
! 
! ******************************************************************** 
! 
      SUBROUTINE GAMHE(IND,T,ANE,ANP,ID,GAM) 
!     ====================================== 
! 
!     NEUTRAL HELIUM STARK BROADENING PARAMETERS 
!     AFTER DIMITRIJEVIC AND SAHAL-BRECHOT, 1984, J.Q.S.R.T. 31, 301 
!     OR FREUDENSTEIN AND COOPER, 1978, AP.J. 224, 1079  (FOR C(IND).GT.0) 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: W(5,20),V(4,20),C(20) 
! 
!   ELECTRONS T= 5000   10000   20000   40000     LAMBDA 
! 
      DATA W /  5.990,  6.650,  6.610,  6.210,    3819.60,                & 
     &          2.950,  3.130,  3.230,  3.300,    3867.50,                & 
     &          0.000,  0.000,  0.000,  0.000,    3871.79,                & 
     &          0.142,  0.166,  0.182,  0.190,    3888.65,                & 
     &          0.000,  0.000,  0.000,  0.000,    3926.53,                & 
     &          1.540,  1.480,  1.400,  1.290,    3964.73,                & 
     &         41.600, 50.500, 57.400, 65.800,    4009.27,                & 
     &          1.320,  1.350,  1.380,  1.460,    4120.80,                & 
     &          7.830,  8.750,  8.690,  8.040,    4143.76,                & 
     &          5.830,  6.370,  6.820,  6.990,    4168.97,                & 
     &          0.000,  0.000,  0.000,  0.000,    4437.55,                & 
     &          1.630,  1.610,  1.490,  1.350,    4471.50,                & 
     &          0.588,  0.620,  0.641,  0.659,    4713.20,                & 
     &          2.600,  2.480,  2.240,  1.960,    4921.93,                & 
     &          0.627,  0.597,  0.568,  0.532,    5015.68,                & 
     &          1.050,  1.090,  1.110,  1.140,    5047.74,                & 
     &          0.277,  0.298,  0.296,  0.293,    5875.70,                & 
     &          0.714,  0.666,  0.602,  0.538,    6678.15,                & 
     &          3.490,  3.630,  3.470,  3.190,    4026.20,                & 
     &          4.970,  5.100,  4.810,  4.310,    4387.93/ 
! 
!   PROTONS   T= 5000   10000   20000   40000 
! 
      DATA V /  1.520,  4.540,  9.140, 10.200,                            & 
     &          0.607,  0.710,  0.802,  0.901,                            & 
     &          0.000,  0.000,  0.000,  0.000,                            & 
     &          0.0396, 0.0434, 0.0476, 0.0526,                           & 
     &          0.000,  0.000,  0.000,  0.000,                            & 
     &          0.507,  0.585,  0.665,  0.762,                            & 
     &          0.930,  1.710, 13.600, 27.200,                            & 
     &          0.288,  0.325,  0.365,  0.410,                            & 
     &          1.330,  6.800, 12.900, 14.300,                            & 
     &          1.100,  1.370,  1.560,  1.760,                            & 
     &          0.000,  0.000,  0.000,  0.000,                            & 
     &          1.340,  1.690,  1.820,  1.630,                            & 
     &          0.128,  0.143,  0.161,  0.181,                            & 
     &          2.040,  2.740,  2.950,  2.740,                            & 
     &          0.187,  0.210,  0.237,  0.270,                            & 
     &          0.231,  0.260,  0.291,  0.327,                            & 
     &          0.0591, 0.0650, 0.0719, 0.0799,                           & 
     &          0.231,  0.260,  0.295,  0.339,                            & 
     &          2.180,  3.760,  4.790,  4.560,                            & 
     &          1.860,  5.320,  7.070,  7.150/ 
      DATA C /2*0.,1.83E-4,0.,1.13E-4,5*0.,1.6E-4,9*0./ 
! 
      IF(W(1,IND).EQ.0.) THEN 
         GAM=C(IND)*T**0.16667*ANE 
       ELSE 
         J=JT(ID) 
         GAM=((TI0(ID)*W(J,IND)+TI1(ID)*W(J-1,IND)+TI2(ID)*W(J-2,IND))    & 
     &        *ANE                                                        & 
     &       +(TI0(ID)*V(J,IND)+TI1(ID)*V(J-1,IND)+TI2(ID)*V(J-2,IND))    & 
     &        *ANP)*1.884E3/W(5,IND)**2 
         IF(GAM.LT.0.) GAM=0. 
      END IF 
      END SUBROUTINE GAMHE 
! 
! ******************************************************************** 
! 
      FUNCTION EPS(T,ANE,ALAM,ION,N) 
!     ============================== 
! 
!     NLTE PARAMETER EPSILON (COLLISIONAL/SPONTANEOUS DEEXCITATION) 
!     AFTER  KASTNER, 1981, J.Q.S.R.T. 26, 377 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      DATA CK0,CK1 /7.75E-8, 2.58E-8/ 
 
      X=1.438E8/ALAM/T 
      XKT=12390./ALAM 
      TT=0.75*X 
      T1=TT+1. 
      A=4.36E7*XKT*XKT/(1.-EXP(-X)) 
      IF(ION.NE.1) THEN 
         B=1.1+LOG(T1/TT)-0.4/T1/T1 
         C=X*B*SQRT(T)/XKT/XKT*ANE 
         IF(N.EQ.0) C=CK0*C 
         IF(N.NE.0) C=CK1*C 
       ELSE 
         C=2.16/T/SQRT(T)/X**1.68*ANE 
      END IF 
      EPS=C/(C+A) 
      RETURN 
      END FUNCTION EPS 
! 
! ******************************************************************** 
! 
      FUNCTION XK2DOP(TAU) 
!     ==================== 
! 
!     KERNEL FUNCTION K2  (AUXILIARY PROCEDURE TO NLTE) 
!     AFTER  HUMMER,  1981, J.Q.S.R.T. 26, 187 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      DATA PI2SQ,PISQ /2.506628275,  1.772453851/ 
      DATA A0,A1,A2,A3,A4 /                                               & 
     &  1.e0,  -1.117897000e-1,  -1.249099917e-1,  -9.136358767e-3,       & 
     &         -3.370280896e-4/ 
      DATA B0,B1,B2,B3,B4,B5 /                                            & 
     &  1.e0,   1.566124168e-1,   9.013261660e-3,   1.908481163e-4,       & 
     &         -1.547417750e-7,  -6.657439727e-9/ 
      DATA C0,C1,C2,C3,C4 /                                               & 
     &  1.0e0,   1.915049608e01,   1.007986843e02,   1.295307533e02,      & 
     &         -3.143372468e01/ 
      DATA D0,D1,D2,D3,D4,D5/                                             & 
     &  1.e0,   1.968910391e01,   1.102576321e02,   1.694911399e02,       & 
     &         -1.669969409e01,  -3.666448000e01/ 
      XK2DOP=1.e0 
      IF(TAU.LE.0.) RETURN 
      IF(TAU.LE.11.) THEN 
         P=A0+TAU*(A1+TAU*(A2+TAU*(A3+TAU*A4))) 
         Q=B0+TAU*(B1+TAU*(B2+TAU*(B3+TAU*(B4+TAU*B5)))) 
         XK2DOP=TAU/PI2SQ*LOG(TAU/PISQ)+P/Q 
       ELSE 
         X=1.e0/LOG(TAU/PISQ) 
         P=C0+X*(C1+X*(C2+X*(C3+X*C4))) 
         Q=D0+X*(D1+X*(D2+X*(D3+X*(D4+X*D5)))) 
         XK2DOP=P/Q/2.D0/TAU/SQRT(LOG(TAU/PISQ)) 
      END IF 
      RETURN 
      END FUNCTION XK2DOP 
! 
! ******************************************************************** 
! 
      SUBROUTINE INKUR 
!     ================ 
! 
!     Input of a Kurucz model atmosphere 
! 
!     Input values (extracted from the Kurucz files): 
!      TEF, G  - effective temperature, log g (appears only in output) 
!      ND      - number of depth points 
!     and for each depth: 
!      DM      - m, m is the mass depth coordinate 
!      T       - temperature 
!      P       - gass pressure 
!      ANE     - electron density 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      READ(8,"(4X,F8.0,9X,F8.5)") TEF,GRAV 
      READ(8,"(/////////////////////10X,I3/)") ND 
      ND=ND-1 
      WRITE(6,"(' INPUT KURUCZ MODEL FOR TEFF=',F7.0,'   LOG G =',        & 
     &       F7.2//1H ,7X,'MASS',9X,'T',9X,'NE',9X,'DENS'/                & 
     &       '-----------------------------------------------'/)")        & 
     &       TEF,GRAV 
      WRITE(77,"(' INPUT KURUCZ MODEL FOR TEFF=',F7.0,'   LOG G =',       & 
     &       F7.2//1H ,7X,'MASS',9X,'T',9X,'NE',9X,'DENS'/                & 
     &       '-----------------------------------------------'/)")        & 
     &       TEF,GRAV 
      DO ID=1,ND 
         READ(8,*) DM(ID),TEMP(ID),P,ELEC(ID) 
         AN=P/TEMP(ID)/BOLK 
         DENS(ID)=WMM(ID)*(AN-ELEC(ID)) 
         WRITE(6,"(I5,1PE10.3,0PF10.1,1P2E12.3)")                         & 
     &      ID,DM(ID),TEMP(ID),ELEC(ID),DENS(ID) 
         WRITE(77,"(I5,1PE10.3,0PF10.1,1P2E12.3)")                        & 
     &      ID,DM(ID),TEMP(ID),ELEC(ID),DENS(ID) 
         T=TEMP(ID) 
         IF(IFMOL.GT.0.AND.T.LT.TMOLIM) THEN 
            AEIN=ELEC(ID) 
            CALL MOLEQ(ID,T,AN,AEIN,ANE,1) 
          ELSE 
            DO IAT=1,NATOM 
               ATTOT(IAT,ID)=DENS(ID)/WMM(ID)/YTOT(ID)*ABUND(IAT,ID) 
            END DO 
         END IF 
         CALL WNSTOR(ID) 
         CALL SABOLF(ID) 
         CALL RATMAT(ID,ESEMAT,BESE) 
         CALL LEVSOL(ESEMAT,BESE,POPLTE,NLEVEL) 
         DO J=1,NLEVEL 
            POPUL(J,ID)=POPLTE(J) 
         END DO 
      END DO 
! 
      CLOSE(8) 
      RETURN 
      END SUBROUTINE INKUR 
! 
! ******************************************************************** 
! 
! 
! 
      SUBROUTINE INPMOD 
!     ================= 
! 
!     Read an initial model atmosphere from unit 8 
!     File 8 contains: 
!      1. NDPTH -  number of depth points in which the initial model is 
!                  given (if not equal to ND, routine interpolates 
!                  automatically to the set DM by linear interpolation 
!                  in log(DM) 
!         NUMPAR - number of input model parameters in each depth 
!                  = 3 for LTE model - ie. N, T, N(electron); 
!                  > 3 for NLTE model) 
!      2. DEPTH(ID),ID=1,NDPTH - mass-depth points for the input model 
!      3. for each depth: 
!                 T   - temperature 
!                 ANE - electron density 
!                 RHO - mass density 
!                 level populations - only for NLTE input model 
!                       Number of input level populations need not be 
!                       equal to NLEVEL; in that case the procedure 
!                       CHANGE is called from START to calculate the 
!                       remaining level populations 
! 
!     Note: The output file 7, which is created by this program 
!           (procedure OUTPUT) has the same structure as file 8 
!           and may thus be used as input to another run of the 
!           program 
!     INTRPL - switch indicating whether (and, if so, how) interpolate 
!              the initial model if the depth scales for the input model 
!              and the present depth scale are different 
!            = 0  -  no interpolation, i.e. scale DEPTH coincides with DM 
!            > 0  -  polynomial interpolation of the (INTRPL-1)th order 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      INTEGER, PARAMETER :: MINPUT=MLEVEL+4 
      REAL(DP) :: TOTN(MDEPTH),PLTE(MLEVEL),X(MINPUT),                    & 
     &          DEPTH(MDEPTH) 
! 
      NUMLT=3 
      IF(INMOD.EQ.2) NUMLT=4 
      READ(8,*) NDPTH,NUMPAR 
      READ(8,*) (DEPTH(I),I=1,NDPTH) 
      notnd=0 
      if(ndpth.ne.nd) notnd=1 
!     write(*,*) 'nd,ndpth,notnd',nd,ndpth,notnd 
      ND=NDPTH 
      NUMP=ABS(NUMPAR) 
      DO ID=1,NDPTH 
         READ(8,*) (X(I),I=1,NUMP) 
         TEMP(ID)=X(1) 
         ELEC(ID)=X(2) 
         DENS(ID)=X(3) 
         if(notnd.eq.1) then 
            wmm(id)=wmm(1) 
            ytot(id)=ytot(1) 
            wmy(id)=wmy(1) 
         end if 
         TOTN(ID)=DENS(ID)/WMM(ID)+ELEC(ID) 
!        write(*,"('inp',i4,1p5e11.3)") id,temp(id),elec(id),dens(id),    & 
!    &                                  totn(id),wmm(id) 
         CALL WNSTOR(ID) 
         CALL SABOLF(ID) 
         IP=NUMLT 
         IF(NUMPAR.LT.0) THEN 
            IP=IP+1 
            TOTN(ID)=X(IP) 
         END IF 
         IF(INMOD.EQ.2) IP=IP+1 
! 
!        first compute LTE level populations for all levels, 
!        i.e. explicit, semi-explisit, and quasi-explicit 
! 
         NLEV0=NLEVEL 
         TEMP(ID)=X(1) 
         ELEC(ID)=X(2) 
         DENS(ID)=X(3) 
         t=temp(id) 
         if(ifmol.gt.0.and.t.lt.tmolim) then 
            ipri=1 
            aein=elec(id) 
            an=totn(id) 
!        write(*,*) 'before moleq',t,an,aein 
            call moleq(id,t,an,aein,ane,ipri) 
!        write(*,*) 'after moleq',t,an,aein,ane 
          else 
            if(imode.gt.-2) then 
               DO IAT=1,NATOM 
                  ATTOT(IAT,ID)=DENS(ID)/WMM(ID)/YTOT(ID)*ABUND(IAT,ID) 
!                 if(id.eq.3) write(*,*) 'attot',iat,dens(id),ytot(id), 
!    *            abund(iat,id),attot(iat,id) 
               END DO 
             else 
               DO IAT=1,NATOM 
                  ATTOT(IAT,ID)=DENS(ID)/WMM(1)/YTOT(1)*ABUND(IAT,1) 
               END DO 
            end if 
         end if 
         CALL WNSTOR(ID) 
         CALL SABOLF(ID) 
         CALL RATMAT(ID,ESEMAT,BESE) 
         CALL LEVSOL(ESEMAT,BESE,POPLTE,NLEV0) 
         DO I=1,NLEV0 
            POPUL(I,ID)=POPLTE(I) 
            PLTE(I)=POPLTE(I) 
         END DO 
! 
!        if the input file fort.8 contains also NLTE level populations 
!        of b-factors, replace the LTE populations by those 
! 
         IF(NUMP.GT.IP) THEN 
            NLEV0=NUMP-IP 
            DO I=1,NLEV0 
               j=iltot(i) 
               POPUL(J,ID)=X(IP+I)*RELAB(IATM(I),ID) 
            END DO 
! 
!           in the case the input "NLTE populations are in fact b-factors, 
!           compute the real populations 
! 
            if(ibfac.eq.1) then 
               do i=1,nlev0 
                  j=iltot(i) 
                  popul(j,id)=popul(j,id)*plte(j) 
               end do 
            end if 
         END IF 
      END DO 
! 
      close(8) 
! 
      write(6,"(/' INPUT TLUSTY MODEL'/                                   & 
     &        ' ------------------'/                                      & 
     &         8X,'MASS',9X,'T',9X,'NE',9X,'DENS',9X,'NH_1'//)") 
         nd=ndpth 
         DO ID=1,ND 
            DM(ID)=DEPTH(ID) 
            write(6,"(i6,1pe10.3,0pf10.1,1p4e12.3)")                      & 
     &       id,dm(id),temp(id),elec(id),dens(id),popul(1,id) 
         END DO 
! 
      DO ID=1,ND 
         BCON=ELEC(ID)/TEMP(ID)/SQRT(TEMP(ID))*2.0706E-16 
         DO IONE=1,NION 
            ION=IZ(IONE) 
            IAT=NUMAT(IATM(NFIRST(IONE))) 
            NKI=NNEXT(IONE) 
            IF(ION.GT.0) PNLT(IAT,ION,ID)=POPUL(NKI,ID)/G(NKI)*BCON 
         END DO 
      END DO 
! 
!     check abundances 
! 
!     CALL CHCKAB 
      RETURN 
      END SUBROUTINE INPMOD 
 
! 
! ******************************************************************** 
! 
! 
 
      SUBROUTINE INPBF 
!     ================ 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      INTEGER,PARAMETER :: MINPUT=MLEVEL+4 
      REAL(DP) ::  DEPTH(MDEPTH),X(MINPUT),XX(MDEPTH),BF(MDEPTH) 
! 
      OPEN(8,FILE='bfactors',STATUS='OLD') 
      NUMLT=3 
      IF(INMOD.EQ.2) NUMLT=4 
      READ(8,*) NDPTH,NUMPAR 
      READ(8,*) (DEPTH(I),I=1,NDPTH) 
      IF(NUMPAR.LT.0) NUMLT=NUMLT+1 
      NUMP=ABS(NUMPAR) 
      NLEV=NUMP-NUMLT 
      DO ID=1,NDPTH 
         READ(8,*) (X(I),I=1,NUMP) 
         DO I=1,NLEV 
            POPUL0(I,ID)=X(I+NUMLT) 
         END DO 
      END DO 
      CLOSE(8) 
! 
!     interpolate the input b-factors to the original DM-scale; 
!     compute new NLTE populations 
! 
      DO I=1,NLEV 
         DO ID=1,NDPTH 
            XX(ID)=POPUL0(I,ID) 
         END DO 
         CALL INTERP(DEPTH,XX,DM,BF,NDPTH,ND,2,1,1) 
         DO ID=1,ND 
            POPUL(I,ID)=POPUL(I,ID)*BF(ID) 
         END DO 
      END DO 
! 
      RETURN 
      END SUBROUTINE INPBF 
 
! 
! 
!     **************************************************************** 
! 
! 
 
      SUBROUTINE LEVSOL(A,B,POPP,NLVCAL) 
!     ================================== 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: A(MLEVEL,MLEVEL),B(MLEVEL),POPP(MLEVEL),                & 
     &          POPP1(MLEVEL) 
! 
!     new populations by inverting several partial rate matrices for the 
!     individual chemical species 
! 
      if(nlvcal.le.0) return 
         ATOMS: DO IAT=1,NATOM 
            N1=N0A(IAT) 
            NK=NKA(IAT) 
            IF(N1.LE.0) THEN 
               DO I=N0A(IAT),NKA(IAT) 
                  N1=I 
                  IF(I.GT.0) EXIT 
               END DO 
            END IF 
            IF(N1.LE.0) CYCle ATOMS 
            NLP=NK-N1+1 
            DO I=N1,NK 
               DO J=N1,NK 
                  ESEMAT(I-N1+1,J-N1+1)=A(I,J) 
               END DO 
               BESE(I-N1+1)=B(I) 
            END DO 
            CALL LINEQS(ESEMAT,BESE,POPP1,NLP,MLEVEL) 
            DO I=N1,NK 
                POPP(I)=POPP1(I-N1+1) 
            END DO 
         END DO ATOMS 
      RETURN 
      END SUBROUTINE LEVSOL 
 
! 
! 
!     **************************************************************** 
! 
 
      SUBROUTINE CHANGE 
!     ================= 
! 
!     This procedure controls an evaluation of initial level 
!     populations in case where the system of explicit levels 
!     (ie. the choice of explicit level, their numbering, or their 
!     total number) is not consistent with that for the input level 
!     populations read by procedure INPMOD. 
!     Obviously, this procedure need be used only for NLTE input models. 
! 
!     Input from unit 5: 
!     For each explicit level, II=1,NLEVEL, the following parameters: 
!      IOLD   -  NE.0 - means that population of this level is 
!                       contained in the set of input populations; 
!                       IOLD is then its index in the "old" (i.e. input) 
!                       numbering. 
!                       All the subsequent parameters have no meaning 
!                       in this case. 
!             -  EQ.0 - means that this level has no equivalent in the 
!                       set of "old" levels. Population of this level 
!                       has thus to be evaluated. 
!      MODE   -    indicates how the population is evaluated: 
!             = 0  - population is equal to the population of the "old" 
!                    level with index ISIOLD, multiplied by REL; 
!             = 1  - population assumed to be LTE, with respect to the 
!                    first state of the next ionization degree whose 
!                    population must be contained in the set of "old" 
!                    (ie. input) populations, with index NXTOLD in the 
!                    "old" numbering. 
!                    The population determined of this way may further 
!                    be multiplied by REL. 
!             = 2  - population determined assuming that the b-factor 
!                    (defined as the ratio between the NLTE and 
!                    LTE population) is the same as the b-factor of 
!                    the level ISINEW (in the present numbering). The 
!                    level ISINEW must have the equivalent in the "old" 
!                    set; its index in the "old" set is ISIOLD, and the 
!                    index of the first state of the next ionization 
!                    degree, in the "old" numbering, is NXTSIO. 
!                    The population determined of this way may further 
!                    be multiplied by REL. 
!             = 3  - level corresponds to an ion or atom which was not 
!                    explicit in the old system; population is assumed 
!                    to be LTE. 
!      NXTOLD  -  see above 
!      ISINEW  -  see above 
!      ISIOLD  -  see above 
!      NXTSIO  -  see above 
!      REL     -  population multiplier - see above 
!                 if REL=0, the program sets up REL=1 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: S = 2.0706E-16 
      IFESE=0 
      LEVELS: DO II=1,NLEVEL 
         READ(ICHANG,*) IOLD,MODE,NXTOLD,ISINEW,ISIOLD,NXTSIO,REL 
         IF(MODE.GE.3) IFESE=IFESE+1 
         IF(REL.EQ.0.) REL=1. 
         DEPTHS: DO ID=1,ND 
            IF(IOLD.NE.0) THEN 
               POPUL0(II,ID)=POPUL(IOLD,ID) 
               CYCLE DEPTHS 
            END IF 
            IF(MODE.EQ.0) THEN 
               POPUL0(II,ID)=POPUL(ISIOLD,ID)*REL 
               CYCLE DEPTHS 
             ELSE 
               T=TEMP(ID) 
               ANE=ELEC(ID) 
               IF(MODE.LE.2) THEN 
                  NXTNEW=NNEXT(IEL(II)) 
                  SB=S/T/SQRT(T)*G(II)/G(NXTNEW)*EXP(ENION(II)/T/BOLK) 
                  IF(MODE.LE.1) THEN 
                     POPUL0(II,ID)=SB*ANE*POPUL(NXTOLD,ID)*REL 
                     CYCLE DEPTHS 
                   ELSE 
                     KK=ISINEW 
                     KNEXT=NNEXT(IEL(KK)) 
                     SBK=S/T/SQRT(T)*G(KK)/G(KNEXT)*                      & 
     &                   EXP(ENION(KK)/T/BOLK) 
                     POPUL0(II,ID)=SB/SBK*POPUL(NXTOLD,ID)/               & 
     &                   POPUL(NXTSIO,ID)*POPUL(ISIOLD,ID)*REL 
                     CYCLE DEPTHS 
                  END IF 
                ELSE 
                  IF(IFESE.EQ.1) THEN 
                     CALL SABOLF(ID) 
                     CALL RATMAT(ID,ESEMAT,BESE) 
                     CALL LINEQS(ESEMAT,BESE,POPLTE,NLEVEL,MLEVEL) 
                  END IF 
                  POPUL0(II,ID)=POPLTE(II) 
               END IF 
            END IF 
         END DO DEPTHS 
      END DO LEVELS 
      DO I=1,NLEVEL 
         DO ID=1,ND 
            POPUL(I,ID)=POPUL0(I,ID) 
         END DO 
      END DO 
      RETURN 
      END SUBROUTINE CHANGE 
! 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE RATMAT(ID,A,B) 
! 
!     LTE RATE MATRIX  (SAHA-BOLTZMANN EQS. + CHARGE CONSERVATION EQ.) 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      real(dp), parameter :: un=1. 
      REAL(DP) :: A(MLEVEL,MLEVEL),B(MLEVEL) 
! 
      ANE=ELEC(ID) 
      DO I=1,NLEVEL 
         B(I)=0. 
         DO J=1,NLEVEL 
            A(J,I)=0. 
         END DO 
      END DO 
! 
      DO IAT=1,NATOM 
         N0I=N0A(IAT) 
         NKI=NKA(IAT) 
         N1I=NKI-1 
         NREFI=NKI 
         DO I=N0I,N1I 
            A(I,I)=1. 
            N=NNEXT(IEL(I)) 
            A(I,N)=-ANE*SBF(I)*WOP(I,ID) 
         END DO 
         DO I=N0I,NKI 
            IL=ILK(I) 
            A(NREFI,I)=UN 
            IF(IL.NE.0) A(NREFI,I)=1.+ANE*USUM(IL) 
         END DO 
         B(NREFI)=ATTOT(IAT,ID) 
      END DO 
! 
      RETURN 
      END SUBROUTINE RATMAT 
! 
!     **************************************************************** 
! 
      SUBROUTINE SABOLF(ID) 
!     ===================== 
! 
!     Saha-Boltzmann factors (SBF) 
!     and "upper sums" - sum of Saha-Boltzmann factors for upper, LTE, 
!     levels which are not included explicitly (USUM), and derivatives 
!     wrt. temperature (T) and electron density (DUSUMN) 
! 
!     Input: ID  - depth index 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: UH=1.5,                                      & 
     &                       CMAX=2.154e4,CCON=2.0706e-16,TWO=2. 
! 
!     DCHI - approximate lowering of ionization potential for neutrals 
!      Actual lowering is DCHI*effective charge, and is considered only 
!      if IUPSUM(ION).GT.0 
! 
      T=TEMP(ID) 
      SQT=SQRT(T) 
      ANE=ELEC(ID) 
      STANE=SQRT(T/ANE) 
      XMAX=CMAX*SQRT(STANE) 
      TK=BOLK*T 
      CON=CCON/T/SQT 
! 
!     Saha-Boltzmann factors 
! 
      IONS: DO ION=1,NION 
         QZ=IZ(ION) 
         CFN=CON/G(NNEXT(ION)) 
         DCH=0. 
         IUPS=IUPSUM(ION) 
         SSBF=0. 
         USUM(ION)=0. 
         nlst=nlast(ion) 
         if(ifwop(nlst).ge.0) then 
             nl1up=nquant(nlst)+1 
          else 
             nl1up=nquant(nlst) 
         end if 
         DO II=NFIRST(ION),NLAST(ION) 
            if(ifwop(ii).lt.0) then 
               E=EH*QZ*QZ/TK 
               SUM=0. 
               DO J=nl1up,NLMX 
                  XJ=J 
                  XI=J*J 
                  X=E/XI 
                  FI=XI*EXP(X)*WNHINT(J,ID) 
                  SUM=SUM+FI 
               END DO 
               g(ii)=sum*two 
               gmer(imrg(ii),id)=g(ii) 
            end if 
            X=ENION(II)/TK 
            if(x.gt.110.) x=110. 
            SB=CFN*G(II)*EXP(X) 
            SBF(II)=SB 
            SSBF=SSBF+SB 
         END DO 
! 
!     Upper sums 
! 
         if(ifwop(nlst).lt.0) cycle ions 
         if(iups.eq.0) then 
! 
!     1. More exact approach - using (exact) partition functions 
! 
            IAT=NUMAT(IATM(NFIRST(ION))) 
            XMX=XMAX*SQRT(QZ) 
            CALL PARTF(IAT,IZ(ION),T,ANE,XMX,U) 
            EE=ENION(NFIRST(ION))/TK 
            if(ee.gt.110.) ee=110. 
            CFE=CFN*EXP(EE) 
            USUM(ION)=CFE*U-SSBF 
            xx=(ssbf-sbf(nfirst(ion)))/sbf(nfirst(ion)) 
            IF(USUM(ION).LT.0.or.ee.ge.109.or.xx.lt.1.e-7) USUM(ION)=0. 
            IF(USUM(ION).LT.0.) USUM(ION)=0. 
! 
!     2. Approximate approach - summation over fixed number of upper 
!        levels, assumed hydrogenic (ie. their ionization energy and 
!        statistical weight hydrogenic) 
! 
          else if(iups.gt.0) then 
            SUM=0. 
            DSUM=0. 
            E=EH*QZ*QZ/TK 
            DO J=NQUANT(NLAST(ION))+1,IUPS 
               XI=J*J 
               X=E/XI 
               FI=XI*EXP(X) 
               SUM=SUM+FI 
            END DO 
            USUM(ION)=SUM*CON*TWO 
! 
!        3. occupation probability form 
! 
         else 
            SUM=0. 
            DSUM=0. 
            E=EH*QZ*QZ/TK 
            DO J=NQUANT(NLAST(ION))+1,NLMX 
               XJ=J 
               XI=J*J 
               X=E/XI 
               FI=XI*EXP(X)*WNHINT(J,ID) 
               SUM=SUM+FI 
            END DO 
            USUM(ION)=SUM*CON*TWO 
         end if 
      END DO IONS 
      RETURN 
      END SUBROUTINE SABOLF 
! 
! ******************************************************************** 
! 
! 
 
      SUBROUTINE OPADD(MODE,ID,FR,ABAD,EMAD,SCAD) 
!     =========================================== 
! 
!     Additional opacities 
!     This is basically user-supplied procedure; here are some more 
!     important non-standard opacity sources, namely 
!     Rayleigh scattering, H- opacity, H2+ opacity, and additional 
!     opacity of He I and He II. 
! 
!     Input parameters: 
!     MODE  - controls the nature and the amount of calculations 
!           = -1 - (OPADD called from START) evaluation of relevant 
!                  depth-dependent quantities (usually photoionization 
!                  cross-sections, but also possibly other), which are 
!                  stored in array CROS 
!           = 0  - evaluation of an additional opacity, emissivity, and 
!                  scattering - for procedure OPAC0 
!     ID    - depth index 
!     FR    - frequency 
! 
!     Output: 
! 
!     ABAD  - absorption coefficient (at frequency FR and depth ID) 
!     EMAD  - emission coefficient (at frequency FR and depth ID) 
!     SCAD  - scattering coefficient (at frequency FR and depth ID) 
! 
! 
      use accura 
      use params 
      use modelp 
      use eospar 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER ::                                              & 
     &           FRAYH  =  2.463E15,                                      & 
     &           FRAYHe =  5.150E15,                                      & 
     &           FRAYH2 =  2.922E15,                                      & 
     &           CLS    =  2.997925e18 
! 
      AB0=0. 
      AB1=0. 
      ABAD=0. 
      EMAD=0. 
      SCAD=0. 
      lpri=id.eq.1.and.cls/fr.le.3501. 
      lpri=id.eq.1 
      al=2.997925e18/fr 
      lpri=.false. 
! 
      if(iath.gt.0) then 
         N0HN=NFIRST(IELH) 
         NKH=NKA(IATH) 
         ANH=POPUL(N0HN,ID) 
         ANP=POPUL(NKH,ID) 
       ELSE 
         anh=ahn(id) 
         anp=ahp(id) 
      END IF 
! 
      if(ielhe1.gt.0) then 
         anhe=popul(NFIRST(IELHE1),ID) 
       else 
         anhe=ahen(id) 
      end if 
! 
      IF(MODE.GE.0) THEN 
         T=TEMP(ID) 
         ANE=ELEC(ID) 
         HKT=HK/T 
         T32=1./T/SQRT(T) 
      END IF 
! 
      IT=NLEVEL 
! 
!   ----------------------- 
!   HI  Rayleigh scattering 
!   ----------------------- 
! 
      IF(IRSCT.NE.0.AND.IOPHLI.NE.1.AND.IOPHLI.NE.2) THEN 
         X=1.e0/(CLS/MIN(FR,FRAYH))**2 
         SG=(5.799E-13+(1.422E-6+2.784*X)*X)*X*X 
!        ABAD=POPUL(N0HN,ID)*SG 
         scad=anh*sg 
         if(lpri) write(*,"('H-Ray ',1p5e11.3)") al,anh,sg,scad 
      END IF 
      IF(IOPHMI.NE.0) THEN 
! 
!   ---------------------------- 
!   H-  bound-free and free-free 
!   ---------------------------- 
!     Note: IOPHMI must not by taken non-zero if H- is considered 
!           explicitly, because H- opacity would be taken twice 
! 
          SG=SBFHMI(FR) 
         if(lpri) write(*,"('H-bf ',1p5e11.3)") al,sg 
          XHM=8762.9/T 
          SB=1.0353E-16*T32*EXP(XHM)*POPUL(N0HN,ID)*ANE*SG 
          SF=SFFHMI(ANH,FR,T)*ANE 
          AB0=SB+SF 
         if(lpri) write(*,"('H-ff ',1p5e11.3)") al,sb,sf,ab0 
      END IF 
! 
!   ----------------------- 
!   He I  Rayleigh scattering 
!   ----------------------- 
! 
      IF(IRSCHE.NE.0.AND.MODE.GE.0) THEN 
         X=(CLS/MIN(FR,FRAYHe))**2 
         CS=5.484E-14/X/X*(1.+(2.44E5+5.94E10/(X-2.90E5))/X)**2 
         sg=anhe*cs 
!        abad=abad+sg 
         scad=scad+sg 
         if(lpri) write(*,"('HeRay',1p5e11.3)") al,cs,anhe,sg,scad 
      END IF 
! 
!   ----------------------- 
!   H2  Rayleigh scattering 
!   ----------------------- 
! 
      IF(IRSCH2.NE.0.AND.MODE.GE.0.AND.IFMOL.GT.0) THEN 
         X=(CLS/MIN(FR,FRAYH2))**2 
           X2=1./X/X 
           CS=(8.14E-13+1.28E-6/X+1.61*X2)*X2 
           sg=cs*anh2(id) 
!          abad=abad+sg 
           scad=scad+sg 
         if(lpri) write(*,"('H2Ray',1p5e11.3)") al,cs,anhe,sg,ab0 
        END IF 
! 
      IF(IOPH2P.GT.0.AND.IFMOL.GT.0.and.                                  & 
     &   t.lt.tmolim.and.fr.lt.3.28e15) THEN 
! 
!   ----------------------------- 
!   H2+  bound-free and free-free 
!   ----------------------------- 
! 
         X=FR*1.E-15 
         SG1=(-7.342E-3+(-2.409+(1.028+(-4.23E-1+                         & 
     &       (1.224E-1-1.351E-2*X)*X)*X)*X)*X)*1.602E-12/BOLK 
         IT=IT+1 
         X=LOG(FR) 
         SG2=-3.0233E3+(3.7797E2+(-1.82496E1+(3.9207E-1-                  & 
     &        3.1672E-3*X)*X)*X)*X 
         X2=-SG1/T+SG2 
         SB=0. 
         IF(X2.GT.-150.) SB=ANH*ANP*EXP(X2) 
         AB0=AB0+SB 
         if(lpri) write(*,"('H2+  ',1p5e11.3)") al,sg1,sg2,sb,ab0 
      END IF 
!     end if 
! 
!   ----------------------------- 
!   He-  free-free 
!   ----------------------------- 
! 
      if(mode.ge.0.and.iophem.gt.0) then 
         C1=3.397e-46+(-5.216e-31+7.039e-15/FR)/FR 
         C2=-4.116e-42+(1.067e-26+8.135e-11/FR)/FR 
         C3=5.081e-37+(-8.724e-23-5.659e-8/FR)/FR 
         cs=c1*t+c2+c3/t 
         sg=anhe*ane*cs 
         ab0=ab0+sg 
         if(lpri) write(*,"('He-ff',1p5e11.3)") al,cs,anhe,sg,ab0 
      end if 
! 
!   ----------------------------- 
!   H2-  free-free 
!   ----------------------------- 
! 
      IF(IOPH2M.NE.0.AND.MODE.GE.0.AND.IFMOL.GT.0.AND.T.LT.TMOLIM) THEN 
         call h2minus(t,anh2(id),ane,fr,oph2) 
         ab1=ab1+oph2 
         if(lpri) write(*,*) 'H2-' 
      END IF 
! 
!   ----------------------------- 
!     CH and OH continuuum opacity 
!   ----------------------------- 
! 
      if(mode.ge.0.and.ifmol.gt.0.and.t.lt.tmolim) then 
         if(iopch.gt.0) ab0=ab0+sbfch(fr,t)*anch(id) 
         if(iopoh.gt.0) ab0=ab0+sbfoh(fr,t)*anoh(id) 
         if(lpri) write(*,*) 'CH,OH' 
! 
!     --------------------------- 
!     CIA H2-H2 opacity 
!     --------------------------- 
! 
         if(ioh2h2.gt.0) then 
            call cia_h2h2(t,anh2(id),fr,oph2) 
            ab1=ab1+oph2 
         if(lpri) write(*,*) 'H2H2 CIA' 
         end if 
! 
!     --------------------------- 
!     CIA H2-He opacity 
!     --------------------------- 
! 
         if(ioh2he.gt.0) then 
            call cia_h2he(t,anh2(id),anhe,fr,oph2) 
            ab1=ab1+oph2 
         if(lpri) write(*,*) 'H2He CIA' 
         end if 
! 
!     --------------------------- 
!     CIA H2-H opacity 
!     --------------------------- 
! 
         if(ioh2h1.gt.0) then 
            call cia_h2h(t,anh2(id),anh,fr,oph2) 
            ab1=ab1+oph2 
         if(lpri) write(*,*) 'H2H CIA' 
         end if 
! 
!     --------------------------- 
!     CIA H-He opacity 
!     --------------------------- 
! 
         if(iohhe.gt.0) then 
            call cia_hhe(t,anh,anhe,fr,oph2) 
            ab1=ab1+oph2 
         if(lpri) write(*,*) 'HHe CIA' 
         end if 
      end if 
! 
!     ---------------------------------------------- 
!     The user may supply more opacity sources here: 
!     ---------------------------------------------- 
! 
!     Finally, actual absorption and emission coefficients 
 
      IF(MODE.LT.0) RETURN 
!      write(*,*) 'before hydhop',id,fr,abad 
      call hydhop(id,fr,abh,emh) 
!     if(lpri) write(*,*) 'after hydhop',id,fr,abh,emh 
      X=EXP(-HKT*FR) 
      X1=1.-X 
      BNX=BN*(FR*1.E-15)**3*X 
      ABAD=ABAD+X1*AB0+AB1+abh 
      EMAD=EMAD+BNX*(AB0+AB1/X1)+emh 
      if(lpri) write(*,"('adtot',1p7e11.3)") al,x,x1,ab0,ab1,abad,emad 
      RETURN 
      END SUBROUTINE OPADD 
! 
! 
!     **************************************************************** 
! 
! 
      subroutine hydhop(id,fr,abad,emad) 
!     ================================== 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
 
      REAL(DP), PARAMETER :: SIH0=2.815e29 
      REAL(DP), PARAMETER :: CPP=4.1412E-16,CPJ=157803. 
      REAL(DP), PARAMETER :: FRH=3.28805E15 
      cha=1. 
! 
!     populations of the first 40 levels of hydrogen 
! 
      t=temp(id) 
      hkt=hk/t 
      ane=elec(id) 
      ANP=POPUL(NKH,ID) 
      PP=CPP*ANE*ANP/T**1.5 
! 
      nlas=n1h 
      if(ifwop(n1h).lt.0) nlas=n1h-1 
      nlash=nlas-n0hn+1 
      sigk=0. 
      abh=0. 
      do ii=nlash+1,40 
         xi=ii 
         xi2=xi*xi 
         fr0=frh/xi2 
         fr1=fr0*(1./xi2-1./(xi+5.)**2) 
         if(fr.ge.fr1) then 
            SIGK=SIH0/FR/FR/FR*CHA*CHA/(xi2**2*xi) 
            frdec=min(fr1*1.25,fr0) 
            if(fr.gt.fr1.and.fr.lt.frdec)                                 & 
     &         sigk=sigk*(fr-fr1)/(frdec-fr1) 
            PJ=PP*EXP(CPJ/XI2/T)*XI2*wnhint(ii,id) 
            abh=abh+pj*sigk 
         end if 
      end do 
! 
      X=EXP(-HKT*FR) 
      X1=1.-X 
      BNX=BN*(FR*1.E-15)**3*X 
      ABAD=X1*ABH 
      EMAD=BNX*ABH 
      RETURN 
      END subroutine hydhop 
 
! 
! 
!     **************************************************************** 
! 
! 
 
      function wn(xn,a,id,z) 
!     ====================== 
! 
!     evaluation of the occupation probablities for a hydrogenic ion 
!     using eqs (4.26), and (4.39) of Hummer,Mihalas Ap.J. 331, 794, 1988. 
!     approximate evaluation of Q(beta) - Hummer 
! 
!     Input: xn  - real number corresponding to quantum number n 
!            a   - correlation parameter 
!            id  - depth index 
!            z   - ionic charge 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      real(dp), parameter :: p1=0.1402,p2=0.1285,p3=1.,p4=3.15,p5=4.,     & 
     &                       tkn=3.01,ckn=5.33333333,cb=8.59e14,          & 
     &                       f23=-2./3.,un=1.,a0=5.29177e-9 
! 
!     evaluation of k(n) 
! 
      if(xn.le.tkn) then 
         xkn=un 
       else 
         xkn=ckn*xn/(xn+un)/(xn+un) 
      end if 
! 
!     evaluation of beta 
! 
!     beta=cb*bergfc*z*z*z*xkn/(xn*xn*xn*xn)*exp(f23*log(elec(id))) 
      beta=cb*z*z*z*xkn/(xn*xn*xn*xn)*exp(f23*log(elec(id))) 
! 
!     approximate expression for Q(beta) 
! 
      x=exp(p4*log(un+p3*a)) 
!     c1=p1*(x+p5*z*a*a*a)    ! previous expression -ERROR !!!!!! 
      c1=p1*(x+p5*(z-un)*a*a*a) 
      c2=p2*x 
      f=(c1*beta*beta*beta)/(un+c2*beta*sqrt(beta)) 
      wp=f/(un+f) 
! 
!     contribution from neutral particles 
! 
!     xn2=xn*xn+un 
      xnh=0. 
      xnhe1=0. 
      if(ielh.gt.0) xnh=popul(nfirst(ielh),id) 
      if(ielhe1.gt.0) xnhe1=popul(nfirst(ielhe1),id) 
      aneut=xnh+xnhe1 
      w0=exp(-4.189*aneut*(a0*(un+xn*xn))**3) 
!     w0=exp(wa0*xn2*xn2*xn2*(xnh+xnhe1)) 
!     W0=1. 
      wn=wp*w0 
      return 
      end function wn 
! 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE WNSTOR(ID) 
!     ===================== 
! 
!     Stores occupation probabilities for hydrogen levels 
!     in common WNCOM for further use 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      real(dp), parameter :: p1=0.1402,p2=0.1285,p3=1.,p4=3.15,p5=4.,     & 
     &                       tkn=3.01,ckn=5.33333333,cb=8.59e14,un=1.,    & 
     &                       f23=-2./3.,ccor=0.09,two=2.,sixth=1./6. 
! 
      ANE=ELEC(ID) 
      A=CCOR*EXP(SIXTH*LOG(ANE))/SQRT(TEMP(ID)) 
      DO I=1,NLMX 
         XN=I 
         WNHINT(I,ID)=wn(xn,a,id,un) 
         WNHE2(I,ID)=wn(xn,a,id,two) 
      END DO 
! 
!     array WOP - occupation probabilities for explicit levels 
! 
      do ii=1,nlevel 
         wop(ii,id)=un 
         if(ifwop(ii).le.0) cycle 
         ie=iel(ii) 
         nq=nquant(ii) 
         if(iz(ie).eq.1) then 
            wop(ii,id)=wnhint(nq,id) 
          else if(iz(ie).eq.2) then 
            wop(ii,id)=wnhe2(nq,id) 
          else 
            z=iz(ie) 
            xn=nq 
            wop(ii,id)=wn(xn,a,id,z) 
         end if 
      END DO 
      RETURN 
      END SUBROUTINE WNSTOR 
! 
! 
! ******************************************************************** 
! 
! 
      subroutine quit(text) 
!     ===================== 
! 
!     stops the program and writes a text 
! 
      use accura 
      use params 
      character*(*) text 
      write(6,"(1x,a)") text 
      stop 
      end subroutine quit 
! 
! 
 
! 
! ******************************************************************* 
! 
 
 
      function voigte(a,vs) 
!     ===================== 
! 
!     computes a voigt function  h = h(a,v) 
!     a=gamma/(4*pi*dnud)   and  v=(nu-nu0)/dnud.  this  is  done after 
!     traving (landolt-bernstein, p. 449). 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      real(dp), parameter :: un=1., two=2.,                               & 
     &                       sqp=1.772453851,sq2=1.414213562 
      real(dp) ::  ak(19),a1(5) 
      data ak      /-1.12470432, -0.15516677,  3.28867591, -2.34357915,   & 
     &  0.42139162, -4.48480194,  9.39456063, -6.61487486,  1.98919585,   & 
     & -0.22041650, 0.554153432, 0.278711796,-0.188325687, 0.042991293,   & 
     &-0.003278278, 0.979895023,-0.962846325, 0.532770573,-0.122727278/ 
! 
      v = abs(vs) 
      u = a + v 
      v2 = v*v 
      ex=0. 
      if(v2.lt.100.) ex = exp(-v2) 
      k = 1 
 
      if(a.eq.0.0) then 
         voigte=0. 
         if(v2.lt.100.) voigte=exp(-v2) 
         return 
       else if(a.gt.0.2) then 
         if(a.gt.1.4.or.u.gt.3.2) then 
            a2 = a*a 
            u = sq2*(a2 + v2) 
            u2 = un/(u*u) 
            voigte = sq2/sqp*a/u*(1. + u2*(3.*v2 - a2) +                  & 
     &               u2*u2*(15.*v2*v2 - 30.*v2*a2 + 3.*a2*a2)) 
            return 
          else 
            ex=0. 
            if(v2.lt.100.)ex = exp(-v2) 
            k = 2 
         end if 
       else if(v.ge.5) then 
         voigte=a*(15. + 6.*v2 + 4.*v2*v2)/(4.*v2*v2*v2*sqp) 
         return 
      end if 
 
      quo=un 
      if(v.lt.2.4) then 
         m = 6 
         if (v.lt.1.3) m = 1 
       else 
         quo=un/(v2-1.5) 
         m=11 
      end if 
      do i=1,5 
         a1(i) = ak(m) 
         m = m + 1 
      end do 
      h1 = quo*(a1(1) + v*(a1(2) + v*(a1(3) + v*(a1(4) + v*a1(5))))) 
 
      if(k.le.1) then 
         voigte=h1*a + ex*(un + a*a*(un - two*v2)) 
         return 
      end if 
 
      pqs=two/sqp 
      h1p = h1 + pqs*ex 
      h2p = pqs*h1p - two*v2*ex 
      h3p = (pqs*(un - ex*(un - two*v2)) - two*v2*h1p)/3. + pqs*h2p 
      h4p = (two*v2*v2*ex - pqs*h1p)/3. + pqs*h3p 
      psi = ak(16) + a*(ak(17) + a*(ak(18) + a*ak(19))) 
      voigte=psi*(ex + a*(h1p + a*(h2p + a*(h3p + a*h4p)))) 
 
      return 
      end function voigte 
! 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE SIGAVS 
!     ================= 
! 
!     Read bound-free cross-sections for averaged levels 
!     from the unit INSA (given by IFANCY), with increasing frequencies 
!     It assumes that all continuum transitions for a given ion are 
!     given in a successive order in the data (i.e. as in TLUSTY for 
!     explicit levels. For other levels, additional input data in 
!     unit 54 !! 
! 
      use accura 
      use params 
      use synthp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: HCCM=H*2.997925e10,BAM=1.e-18 
      REAL(DP) ::  CRD(MFCRA),XIFE(8),FRD(MFCRA) 
! 
      DATA XIFE/63480.,130563.,247220.,442000.,605000.,799000.,           & 
     &          1008000.,1218380./ 
! 
      FR1=FREQ(1) 
      FR2=FREQ(2) 
      NUNIT=0 
      NQHT=0 
      IONS: DO I=1,NION 
        IF(IASV.EQ.0) EXIT IONS 
        N1=NFIRST(I) 
        N2=NLAST(I) 
        INSA=0 
 
        DO II=N1,N2 
          NFCR(II)=2 
          FRECR(II,1)=FR1 
          FRECR(II,2)=FR2 
          CROSR(II,1)=0. 
          CROSR(II,2)=0. 
          INSB=IBF(II) 
          IF(INSB.LT.50.OR.INSB.GT.100) CYCLE 
          IF(INSA.EQ.0) INSA=INSB 
          IF(INSA.NE.INSB)                                                & 
     &     call quit(' Incoherent file units in SIGAVS') 
        END DO 
 
        IF(INSA.EQ.0) CYCLE IONS 
        IF(FIBFCS(I).NE.' ') THEN 
           INSA=INBFCS(I) 
           OPEN(INSA,FILE=FIBFCS(I),STATUS='OLD') 
        END IF 
        READ(INSA,*,IOSTAT=IOS) IIAT,IIZ,NSUP 
        IF(IOS.NE.0) THEN 
           call quit('ERROR IN FILE FOR BF SIG OF AVER.LEVELS (1)') 
        END IF 
 
        ATI=IIAT+0.01*(IIZ-1) 
        NBFI=NSUP 
        IF(NSUP.GT.(N2-N1+1)) NBFI=(N2-N1+1) 
 
        LEVELS: DO II=1,NBFI 
          READ(INSA,*,IOSTAT=IOS) IILO,EELO,GGLO,NFCRR 
          IF(IOS.NE.0) THEN 
            call quit('ERROR IN FILE FOR BF SIG OF AVER.LEVELS (2)') 
          END IF 
          IK=N1+IILO-1 
          IF (IK.GT.N2 .OR. IK.LT.N1)                                     & 
     &      call quit(' Inconsistent level numbering in SIGAVS') 
          IF(IIAT.EQ.26) ECMR=XIFE(IIZ)-EELO 
          READ(INSA,*,IOSTAT=IOS) FR0,CR0 
          IF(IOS.NE.0) THEN 
             call quit('ERROR IN FILE FOR BF SIG OF AVER.LEVELS (1)') 
          END IF 
          NFD=1 
          FRD(NFD)=FR0 
          CRD(NFD)=CR0 
          LUV=.FALSE. 
 
          FRQ: DO IJ=1,NFCRR-1 
            READ(INSA,*,IOSTAT=IOS2) FRIN,CRIN 
            IF(IOS2.NE.0) THEN 
              call quit('ERROR IN FILE FOR BF SIG OF AVER.LEVELS (1)') 
            END IF 
            IF(LUV) CYCLE FRQ 
            IF(FRIN.GT.FR1) THEN 
              IF(FR0.LE.FR2.AND.IJ.GT.1) THEN 
                NFD=NFD+1 
                FRD(NFD)=FR0 
                CRD(NFD)=CR0 
              ENDIF 
              NFD=NFD+1 
              FRD(NFD)=FRIN 
              CRD(NFD)=CRIN 
              LUV=.TRUE. 
            ELSE IF(FRIN.GT.FR2) THEN 
              IF(FR0.LE.FR2.AND.IJ.GT.1) THEN 
                NFD=NFD+1 
                FRD(NFD)=FR0 
                CRD(NFD)=CR0 
              ENDIF 
              NFD=NFD+1 
              FRD(NFD)=FRIN 
              CRD(NFD)=CRIN 
              FR0=FRIN 
              CR0=CRIN 
            ELSE 
              FR0=FRIN 
              CR0=CRIN 
            ENDIF 
            IF(NFD.GT.MFCRA)                                              & 
     &        call quit(' Too many frequencies in SIGAVS') 
          END DO FRQ 
          CRMX(IK)=0. 
          DO IJ=1,NFD 
            CRMX(IK)=MAX(CRMX(IK),CRD(IJ)) 
          END DO 
          IF(CRMX(IK).GT.0.) THEN 
            NFCR(IK)=NFD 
            DO IJ=1,NFD 
              FRECR(IK,IJ)=FRD(NFD-IJ+1) 
              CROSR(IK,IJ)=CRD(NFD-IJ+1)*BAM 
            END DO 
          ENDIF 
         END DO LEVELS 
      END DO IONS 
! 
      READ(50,*,IOSTAT=IOS3) NUNIT 
      IF(IOS3.NE.0) RETURN 
      IF(NUNIT.LE.0) RETURN 
      WRITE(6,"(///,' DETAILED PHOTOIONIZATION CROSS-SECTIONS',           & 
     & ' (NON-EXPLICIT LEVELS)',/,                                        & 
     & ' ---------------------------------------',/)") 
      FILES: DO IN=1,NUNIT 
        READ(50,*,IOSTAT=IOS4) ATIR,INSA,NQHTR 
        IF(IOS4.NE.0) RETURN 
        NQHT=NQHT+NQHTR 
        IF(NQHT.GT.MPHOT)                                                 & 
     &    call quit(' Too many BF cross-sections in SIGAVS') 
        READ(INSA,*,IOSTAT=IOS) IIAT,IIZ,NSUP 
        IF(IOS2.NE.0) THEN 
           call quit('ERROR IN FILE FOR BF SIG OF AVER.LEVELS (2)') 
        END IF 
! 
!       check the total number of superlevels 
! 
        IF(NQHTR.GT.NSUP) THEN 
           WRITE(6,"('NQHTR=',i4,' in Unit 50 input greater than NSUP=',  & 
     &              i4,/' program resets NQHTR to NSUP'/)") NQHTR,NSUP 
           NQHTR=NSUP 
        END IF 
! 
!       loop over superlevels - read cross-sections 
! 
        SUPLEVELS: DO I=1,NQHTR 
          IK=NQHT-NQHTR+I 
          READ(INSA,*,IOSTAT=IOS5) IILO,EELO,GGLO,NFCRR 
          IF(IOS5.NE.0) THEN 
             call quit('ERROR IN FILE FOR BF SIG OF AVER.LEVELS (2)') 
          END IF 
          AQHT(IK)=ATIR 
          EQHT(IK)=EELO 
          GQHT(IK)=GGLO 
          READ(INSA,*) FR0,CR0 
          NFD=1 
          FRD(NFD)=FR0 
          CRD(NFD)=CR0 
          LUV=.FALSE. 
          FRQ2: DO IJ=1,NFCRR-1 
            READ(INSA,*) FRIN,CRIN 
            IF(LUV) CYCLE FRQ2 
            IF(FRIN.GT.FR1) THEN 
              IF(FR0.LE.FR2.AND.IJ.GT.1) THEN 
                NFD=NFD+1 
                FRD(NFD)=FR0 
                CRD(NFD)=CR0 
              ENDIF 
              NFD=NFD+1 
              FRD(NFD)=FRIN 
              CRD(NFD)=CRIN 
              LUV=.TRUE. 
            ELSE IF(FRIN.GT.FR2) THEN 
              IF(FR0.LE.FR2.AND.IJ.GT.1) THEN 
                NFD=NFD+1 
                FRD(NFD)=FR0 
                CRD(NFD)=CR0 
              ENDIF 
              NFD=NFD+1 
              FRD(NFD)=FRIN 
              CRD(NFD)=CRIN 
              FR0=FRIN 
              CR0=CRIN 
            ELSE 
              FR0=FRIN 
              CR0=CRIN 
            ENDIF 
          END DO FRQ2 
          CRMY(IK)=0. 
          DO IJ=1,NFD 
            CRMY(IK)=MAX(CRMY(IK),CRD(IJ)) 
          END DO 
          IF(CRMY(IK).GT.0.) THEN 
            WRITE(6,"(F7.2,I6,F13.3,I8)") ATIR,IILO,EELO,NFD 
            NFQHT(IK)=NFD 
            DO IJ=1,NFD 
              FRECQ(IK,IJ)=FRD(NFD-IJ+1) 
              QHOT(IK,IJ)=CRD(NFD-IJ+1)*BAM 
            END DO 
          ENDIF 
        END DO SUPLEVELS 
      END DO FILES 
 
      RETURN 
      END SUBROUTINE SIGAVS 
! 
! 
! ******************************************************************** 
! 
! 
 
      SUBROUTINE PHTX(ID,ABSO,EMIS,fre,icon) 
!     ====================================== 
! 
!     Opacity due to detailed photoionization (read from tables by 
!     routine SIGAVS) 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use lindat 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) ::           ABSO(MFREQ),EMIS(MFREQ),PLANF(MFREQ),         & 
     &                      STIMU(MFREQ),FRE(MFREQ) 
      INTEGER, SAVE      :: IJP(MLEVEL),IJQ(MPHOT) 
      REAL(DP),PARAMETER :: C3=1.4387886 
! 
      IF(IASV.EQ.0 .AND. NQHT.EQ.0) RETURN 
      T=TEMP(ID) 
      nfre=nfreq 
      ij0=3 
      if(icon.eq.1) then 
         ij0=1 
         nfre=nfreqc 
      end if 
! 
      DO IJ=1,NFRE 
         XX=FRE(IJ) 
         X15=XX*1.E-15 
         BNU=BN*X15*X15*X15 
         HKF=HK*XX 
         EXH=EXP(HKF/T) 
         PLANF(IJ)=BNU/(EXH-1.) 
         STIMU(IJ)=1.-1./EXH 
      END DO 
! 
      IF(IASV.NE.0) THEN 
         IF(ID.EQ.1) THEN 
           LEVELS: DO I=1,NLEVEL 
             IF(CRMX(I).EQ.0.) CYCLE LEVELS 
             IK1=MAX0(2,IJP(I)) 
             DO IJ=3,NFRE 
               IKLOOP: DO IK=IK1,NFCR(I) 
                 IF(FRECR(I,IK).LT.FRE(IJ)) THEN 
                   IK2=IK 
                   EXIT IKLOOP 
                 ENDIF 
               END DO IKLOOP 
               IK1=IK2 
               IF(IJ.EQ.3) IJP(I)=IK1 
               DFR=(FRE(IJ)-FRECR(I,IK1))/(FRECR(I,IK1-1)-FRECR(I,IK1)) 
               PHOTI(I,IJ)=CROSR(I,IK1)+DFR*                              & 
     &                     (CROSR(I,IK1-1)-CROSR(I,IK1)) 
             END DO 
             PHOTI(I,1)=PHOTI(I,3) 
             PHOTI(I,2)=PHOTI(I,NFREQ) 
           END DO LEVELS 
         END IF 
 
         LEVELS2: DO I=1,NLEVEL 
           IF(CRMX(I).EQ.0.) CYCLE LEVELS2 
           POP=POPUL(I,ID) 
           DO IJ=1,NFRE 
             ABA=PHOTI(I,IJ)*POP*STIMU(IJ) 
             ABSO(IJ)=ABSO(IJ)+ABA 
             EMIS(IJ)=EMIS(IJ)+ABA*PLANF(IJ) 
           END DO 
         END DO LEVELS2 
 
      END IF 
! 
      IF(NQHT.EQ.0) RETURN 
      IF(ID.EQ.1) THEN 
        NQLOOP: DO I=1,NQHT 
          IF(CRMY(I).EQ.0.) CYCLE NQLOOP 
          IK1=MAX0(2,IJQ(I)) 
          FRQLOOP: DO IJ=3,NFRE 
            IK2LOOP: DO IK=IK1,NFQHT(I) 
              IF(FRECQ(I,IK).LT.FRE(IJ)) THEN 
                IK2=IK 
                EXIT IK2LOOP 
              ENDIF 
            END DO IK2LOOP 
            IK1=IK2 
            IF(IJ.EQ.3) IJQ(I)=IK1 
            DFR=(FRE(IJ)-FRECQ(I,IK1))/(FRECQ(I,IK1-1)-FRECQ(I,IK1)) 
            PHOTI(I,IJ)=QHOT(I,IK1)+DFR*(QHOT(I,IK1-1)-QHOT(I,IK1)) 
          END DO FRQLOOP 
        END DO NQLOOP 
      END IF 
 
      NQ2LOOP: DO I=1,NQHT 
        IF(CRMY(I).EQ.0.) CYCLE NQ2LOOP 
        IAT=int(AQHT(I)) 
        X=(AQHT(I)-FLOAT(IAT)+1.E-4)*100. 
        ION=INT(X)+1 
        POP=RRR(ID,ION,IAT)*GQHT(I)*EXP(-EQHT(I)*C3/T) 
        DO IJ=3,NFRE 
          ABA=PHOTI(I,IJ)*POP*STIMU(IJ) 
          ABSO(IJ)=ABSO(IJ)+ABA 
          EMIS(IJ)=EMIS(IJ)+ABA*PLANF(IJ) 
        END DO 
      END DO NQ2LOOP 
! 
      RETURN 
      END SUBROUTINE PHTX 
! 
! 
! ******************************************************************** 
! 
      subroutine getlal 
!     ================= 
! 
!     getlal reads in the profile functions for Lyman alpha, beta, gamma, 
!     and Balmer alpha, including the quasi-molecular satellites; 
!     valid for first and second order in neutral and ionized H density 
!     modified routine provided originally by D. Koester 
! 
! 
      use accura 
      use params 
      use allarn 
      implicit real(dp) (a-h,o-z),logical (l) 
! 
!     Lyman alpha 
! 
      if(nunalp.gt.0.or.nunbet.gt.0.or.nungam.gt.0.or.nunbal.gt.0) then 
         call alloc_allarn 
      end if 
 
      write(*,*) 'getlal',nunalp,nunbet,nungam,nunbal
      nxalp=0 
      if(nunalp.gt.0) then 
         nunalp=67 
         open(unit=nunalp,file='./data/laquasi.dat',status='old') 
         read(nunalp,*) nxalp,stnnea,stncha,vneua,vchaa 
         do i=1,nxalp 
            read(nunalp,*) xlalp(i),(plalp(i,j),j=1,NNMAX) 
         end do 
         close(nunalp) 
         stnnea=10.0**stnnea 
         stncha=10.0**stncha 
         iwarna=0 
         close(nunalp) 
         write(*,*) 
         write(*,*) ' read quasi-molecular data for L alpha' 
      end if 
! 
!     Lyman beta 
! 
      nxbet=0 
      if(nunbet.gt.0) then 
         nunbet=67 
         open(unit=nunbet,file='./data/lbquasi.dat',status='old') 
         read(nunbet,*) nxbet,stnneb,stnchb,vneub,vchab 
         do i=1,nxbet 
            read(nunbet,*) xlbet(i),(plbet(i,j),j=1,NNMAX) 
         end do 
         close(nunbet) 
         stnneb=10.0**stnneb 
         stnchb=10.0**stnchb 
         iwarnb=0 
         write(*,*) ' read quasi-molecular data for L beta' 
      end if 
! 
!     Lyman gamma 
! 
      nxgam=0 
      if(nungam.gt.0) then 
         nungam=67 
         open(unit=nunalp,file='./data/lgquasi.dat',status='old') 
         read(nungam,*) nxgam,stnneg,stnchg,vneug,vchag 
         do i=1,nxgam 
            read(nungam,*) xlgam(i),(plgam(i,j),j=1,NNMAX) 
         end do 
         close(nungam) 
         stnneg=10.0**stnneg 
         stnchg=10.0**stnchg 
         iwarng=0 
         write(*,*) ' read quasi-molecular data for L gamma' 
      end if 
! 
!     Balmer alpha 
! 
      nxbal=0 
      if(nunbal.gt.0) then 
         nunbal=67 
         open(unit=nunalp,file='./data/lhquasi.dat',status='old') 
         read(nunbal,*) nxbal,stnnec,stnchc,vneuc,vchac 
         do i=1,nxbal 
            read(nunbal,*) xlbal(i),(plbal(i,j),j=1,NNMAX) 
         end do 
         close(nunbal) 
         stnnec=10.0**stnnec 
         stnchc=10.0**stnchc 
         iwarnc=0 
         write(*,*) ' read quasi-molecular data for H alpha' 
      end if 
      write(*,*) 
      return 
      end subroutine getlal 
! 
! 
! ******************************************************************** 
! 
      subroutine allard(xl,hneutr,hcharg,prof,iq,jq) 
!     ============================================== 
! 
!     quasi-molecular opacity for Lyman alpha, beta, and Balmer alpha 
!     modified routine provided originally by D. Koester 
! 
!     Input:  xl:  wavelength in [A] 
!             hneutr:  neutral H particle density [cm-3] 
!             hcharg: ionized H particle density [cm-3] 
!             iq:   quantum number of the lower level 
!             jq:   quantum number of the upper level; 
!                   =2  -  Lyman alpha 
!                   =3  -  Lyman beta 
!     Output: prof:  Lyman alpha line profile, normalized to 1.0e8 
!             if integrated over A; 
!             It then renormalized by multiplying by 
!             8.853e-29*lambda_0^2*f_ij 
! 
      use accura 
      use params 
      use allarn 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      real(dp), parameter ::                                              & 
     &           xnorma=8.8528e-29*1215.6*1215.6*0.41618,                 & 
     &           xnormb=8.8528e-29*1025.73*1025.7*0.0791,                 & 
     &           xnormg=8.8528e-29*972.53*972.53*0.0290,                  & 
     &           xnormc=8.8528e-29*6562.*6562.*0.6407 
 
      prof=0. 
! 
!     Lyman alpha 
! 
      if(iq.eq.1.and.jq.eq.2) then 
         if(xl.lt.xlalp(1)) return 
         vn1=hneutr/stnnea 
         vn2=hcharg/stncha 
         vns=vn1*vneua+vn2*vchaa 
         write(*,*) 'vn',vn1,vneua,vn1*vneua,vn2,vchaa,vn2*vchaa
         write(*,*) 'iwarna',iwarna
         if(iwarna.eq.0) then 
            if(vn1*vneua.gt.0.3.or.vn2*vchaa.gt.0.3) then 
               write(*,*) '          warning: density too high for',      & 
     &           ' Lyman alpha expansion' 
               iwarna=1 
            endif 
         endif 
         vn11=vn1*vn1 
         vn22=vn2*vn2 
         vn12=vn1*vn2 
         xnorm=1.0/(1.0+vns+0.5*vns*vns) 
! 
         if(xl.le.xlalp(nxalp)) then 
            jl=0 
            ju=nxalp+1 
            loca: do 
               if(ju-jl.gt.1) then 
                  jm=(ju+jl)/2 
                  if((xlalp(nxalp).gt.xlalp(1)).eqv.(xl.gt.xlalp(jm)))    & 
     &             then 
                     jl=jm 
                   else 
                     ju=jm 
                  end if 
                  cycle loca 
                else 
                  exit loca 
               end if 
            end do loca 
            j=jl 
! 
            if(j.eq.0) j=1 
            if(j.eq.nxalp) j=j-1 
            a1=(xl-xlalp(j))/(xlalp(j+1)-xlalp(j)) 
            p1=  vn1*((1.0-a1)*plalp(j,1)+a1*plalp(j+1,1)) 
            p11=vn11*((1.0-a1)*plalp(j,2)+a1*plalp(j+1,2)) 
            p2=  vn2*((1.0-a1)*plalp(j,3)+a1*plalp(j+1,3)) 
            p22=vn22*((1.0-a1)*plalp(j,4)+a1*plalp(j+1,4)) 
            p12=vn12*((1.0-a1)*plalp(j,5)+a1*plalp(j+1,5)) 
            prof=(p1+p2+p11+p22+p12)*xnorm*xnorma 
! 
          else 
            j=nxalp-1 
            a1=1. 
            p1=  vn1*((1.0-a1)*plalp(j,1)+a1*plalp(j+1,1)) 
            p11=vn11*((1.0-a1)*plalp(j,2)+a1*plalp(j+1,2)) 
            p2=  vn2*((1.0-a1)*plalp(j,3)+a1*plalp(j+1,3)) 
            p22=vn22*((1.0-a1)*plalp(j,4)+a1*plalp(j+1,4)) 
            p12=vn12*((1.0-a1)*plalp(j,5)+a1*plalp(j+1,5)) 
            pro0=(p1+p2+p11+p22+p12)*xnorm*xnorma 
            xlas=xlalp(nxalp) 
            x0=1215.67 
            dxlas=xlalp(nxalp)-x0 
            dx=xl-x0 
            prof=pro0/(dx/dxlas)**2.5 
! 
         end if 
         return 
      end if 
! 
!     Lyman beta 
! 
      if(iq.eq.1.and.jq.eq.3) then 
         if(nxbet.eq.0) return 
         if(xl.lt.xlbet(1).or.xl.gt.xlbet(nxbet)) return 
         vn1=hneutr/stnneb 
         vn2=hcharg/stnchb 
         vns=vn1*vneub+vn2*vchab 
         if(iwarnb.eq.0) then 
            if(vn1*vneub.gt.0.3.or.vn2*vchab.gt.0.3) then 
               write(*,*) '          warning: density too high for',      & 
     &           ' Lyman beta expansion' 
               iwarnb=1 
            endif 
         endif 
         vn11=vn1*vn1 
         vn22=vn2*vn2 
         vn12=vn1*vn2 
         xnorm=1.0/(1.0+vns+0.5*vns*vns) 
! 
         jl=0 
         ju=nxbet+1 
         locb: do 
            if(ju-jl.gt.1) then 
               jm=(ju+jl)/2 
               if((xlbet(nxbet).gt.xlbet(1)).eqv.(xl.gt.xlbet(jm))) then 
                  jl=jm 
                else 
                  ju=jm 
               endif 
               cycle locb 
             else 
               exit locb 
            end if 
         end do locb 
         j=jl 
! 
         if(j.eq.0) j=1 
         if(j.eq.nxbet) j=j-1 
         a1=(xl-xlbet(j))/(xlbet(j+1)-xlbet(j)) 
         p1=  vn1*((1.0-a1)*plbet(j,1)+a1*plbet(j+1,1)) 
         p11=vn11*((1.0-a1)*plbet(j,2)+a1*plbet(j+1,2)) 
         p2=  vn2*((1.0-a1)*plbet(j,3)+a1*plbet(j+1,3)) 
         p22=vn22*((1.0-a1)*plbet(j,4)+a1*plbet(j+1,4)) 
         p12=vn12*((1.0-a1)*plbet(j,5)+a1*plbet(j+1,5)) 
         prof=(p1+p2+p11+p22+p12)*xnorm*xnormb 
         return 
      end if 
! 
!     Lyman gamma 
! 
      if(iq.eq.1.and.jq.eq.4) then 
         if(nxgam.eq.0) return 
         if(xl.lt.xlgam(1).or.xl.gt.xlgam(nxgam)) return 
         vn1=hneutr/stnneg 
         vn2=hcharg/stnchg 
         vns=vn1*vneug+vn2*vchag 
         if(iwarng.eq.0) then 
            if(vn1*vneug.gt.0.3.or.vn2*vchag.gt.0.3) then 
               write(*,*) '          warning: density too high for',      & 
     &           ' Lyman gamma expansion' 
               iwarng=1 
            endif 
         endif 
         vn11=vn1*vn1 
         vn22=vn2*vn2 
         vn12=vn1*vn2 
         xnorm=1.0/(1.0+vns+0.5*vns*vns) 
! 
         jl=0 
         ju=nxgam+1 
         locg: do 
            if(ju-jl.gt.1) then 
               jm=(ju+jl)/2 
               if((xlgam(nxgam).gt.xlgam(1)).eqv.(xl.gt.xlgam(jm))) then 
                  jl=jm 
               else 
                  ju=jm 
               end if 
               cycle locg 
             else 
               exit locg 
            end if 
         end do locg 
         j=jl 
! 
         if(j.eq.0) j=1 
         if(j.eq.nxgam) j=j-1 
         a1=(xl-xlgam(j))/(xlgam(j+1)-xlgam(j)) 
         p1=  vn1*((1.0-a1)*plgam(j,1)+a1*plgam(j+1,1)) 
         p11=vn11*((1.0-a1)*plgam(j,2)+a1*plgam(j+1,2)) 
         p2=  vn2*((1.0-a1)*plgam(j,3)+a1*plgam(j+1,3)) 
         p22=vn22*((1.0-a1)*plgam(j,4)+a1*plgam(j+1,4)) 
         p12=vn12*((1.0-a1)*plgam(j,5)+a1*plgam(j+1,5)) 
         prof=(p1+p2+p11+p22+p12)*xnorm*xnormg 
         return 
      end if 
! 
!     Balmer alpha 
! 
      if(iq.eq.2.and.jq.eq.3) then 
         if(xl.lt.xlbal(1).or.xl.gt.xlbal(nxbal)) return 
         vn1=0. 
         vn2=hcharg/stnchc 
         vns=vn1*vneuc+vn2*vchac 
         vn11=vn1*vn1 
         vn22=vn2*vn2 
         vn12=vn1*vn2 
         xnorm=1.0/(1.0+vns+0.5*vns*vns) 
! 
         jl=0 
         ju=nxbal+1 
         loch: do 
            if(ju-jl.gt.1) then 
               jm=(ju+jl)/2 
               if((xlbal(nxbal).gt.xlbal(1)).eqv.(xl.gt.xlbal(jm))) then 
                  jl=jm 
               else 
                  ju=jm 
               endif 
               cycle loch 
             else 
               exit loch 
            end if 
         end do loch 
         j=jl 
! 
         if(j.eq.0) j=1 
         if(j.eq.nxbal) j=j-1 
         a1=(xl-xlbal(j))/(xlbal(j+1)-xlbal(j)) 
         p1=  vn1*((1.0-a1)*plbal(j,1)+a1*plbal(j+1,1)) 
         p11=vn11*((1.0-a1)*plbal(j,2)+a1*plbal(j+1,2)) 
         p2=  vn2*((1.0-a1)*plbal(j,3)+a1*plbal(j+1,3)) 
         p22=vn22*((1.0-a1)*plbal(j,4)+a1*plbal(j+1,4)) 
         p12=vn12*((1.0-a1)*plbal(j,5)+a1*plbal(j+1,5)) 
         prof=(p1+p2+p11+p22+p12)*xnorm*xnormc 
      end if 
! 
      return 
      end subroutine allard 
! 
! 
!     ************************************************************** 
! 
! 
      subroutine lyahhe(xl,ahe,prof) 
!     ============================== 
! 
!     Lyman alpha broadening by helium - after N. Allard 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      integer, parameter  :: nxmax=1000 
!     real(dp), parameter :: sthe=1.e21 
      real(dp), save      ::  xlhhe(nxmax),sighhe(nxmax) 
      real(dp)            ::  xlhh0(nxmax),sighh0(nxmax) 
      integer, save       :: nxhhe 
      data iread/0/ 
! 
      if(iread.eq.0) then 
         it=0 
         do i=1,nxmax 
            read(67,*,iostat=ios) xl,sig 
            if(ios.ne.0) exit 
            it=it+1 
            if(nunhhe.eq.1) xl=1./(1.e-8*xl+1./1215.67) 
            xlhh0(it)=xl 
            sighh0(it)=sig 
         end do 
         nxhhe=it 
         do i=1,nxhhe 
            xlhhe(i)=xlhh0(nxhhe-i+1) 
            sighhe(i)=sighh0(nxhhe-i+1) 
         end do 
         close(67) 
         iread=1 
      end if 
! 
      prof=0. 
      if(xl.gt.xlhhe(nxhhe)) return 
      jl=0 
      ju=nxhhe+1 
      loc: do 
         if(ju-jl.gt.1) then 
            jm=(ju+jl)/2 
            if((xlhhe(nxhhe).gt.xlhhe(1)).eqv.(xl.gt.xlhhe(jm))) then 
               jl=jm 
             else 
               ju=jm 
            end if 
            cycle loc 
          else 
            exit loc 
         end if 
      end do loc 
      j=jl 
! 
      if(j.eq.0) j=1 
      if(j.eq.nxhhe) j=j-1 
      a1=(xl-xlhhe(j))/(xlhhe(j+1)-xlhhe(j)) 
      s1=(1.0-a1)*sighhe(j)+a1*sighhe(j+1) 
      prof=s1*ahe/sthe*6.2831855 
      return 
      end subroutine lyahhe 
! 
! 
!     ************************************************************** 
! 
! 
      subroutine readbf 
!     ================= 
! 
!     auxiliary subroutine for enabling reading of input data with 
!     comments 
! 
!     lines beginning with ! or * are understood as comments 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      character(len=80) buff 
 
      readin: do 
         read(5,"(a)",iostat=ios) buff 
         if(ios.ne.0) exit readin 
         if(buff(1:1).eq.'!'.or.buff(1:1).eq.'*') cycle readin 
         write(ibuff,"(a)") buff 
      end do readin 
      rewind ibuff 
      return 
      end subroutine readbf 
! 
! 
!     ******************************************************************* 
! 
! 
 
      SUBROUTINE PRETAB 
!     ================= 
! 
!     pretabulate expansion coefficients for the Voigt function 
!     200 steps per doppler width - up to 10 Doppler widths 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) ::  TABVI(81),TABH1(81) 
      DATA TABVI/0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.,1.1,1.2,1.3,1.4,1.5,    & 
     &1.6,1.7,1.8,1.9,2.,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.,3.1,3.2,  & 
     & 3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,   & 
     & 5.8,6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4,7.6,7.8,8.0,8.2,8.4,8.6,8.8,   & 
     & 9.0,9.2,9.4,9.6,9.8,10.0,10.2,10.4,10.6,10.8,11.0,11.2,11.4,11.6,  & 
     & 11.8,12.0/ 
      DATA TABH1/-1.12838,-1.10596,-1.04048,-.93703,-.80346,-.64945,      & 
     &-.48552,-.32192,-.16772,-.03012,.08594,.17789,.24537,.28981,        & 
     &.31394,.32130,.31573,.30094,.28027,.25648,.231726,.207528,.184882,  & 
     &.164341,.146128,.130236,.116515,.104739,.094653,.086005,.078565,    & 
     & .072129,.066526,.061615,.057281,.053430,.049988,.046894,.044098,   & 
     & .041561,.039250,.035195,.031762,.028824,.026288,.024081,.022146,   & 
     & .020441,.018929,.017582,.016375,.015291,.014312,.013426,.012620,   & 
     & .0118860,.0112145,.0105990,.0100332,.0095119,.0090306,.0085852,    & 
     & .0081722,.0077885,.0074314,.0070985,.0067875,.0064967,.0062243,    & 
     & .0059688,.0057287,.0055030,.0052903,.0050898,.0049006,.0047217,    & 
     & .0045526,.0043924,.0042405,.0040964,.0039595/ 
! 
      N=MVOI 
      VSTEPS=200. 
      DO I=1,N 
         H0TAB(I)=FLOAT(I-1)/VSTEPS 
      END DO 
      CALL INTERP(TABVI,TABH1,H0TAB,H1TAB,81,N,2,0,0) 
      DO I=1,N 
         VV=(FLOAT(I-1)/VSTEPS)**2 
         H0TAB(I)=EXP(-VV) 
         H2TAB(I)=H0TAB(I)-(VV+VV)*H0TAB(I) 
      END DO 
      RETURN 
      END SUBROUTINE PRETAB 
! 
! 
!     ******************************************************************* 
! 
! 
 
      FUNCTION VOIGTK(A,V) 
!     ==================== 
! 
!     Voigt function after Kurucz (in Computational Astrophysics) 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER ::                                              & 
     &           ONE=1., THREE=3., TEN=10., FIFTN=15., TWOH=200.,         & 
     &           C14142=1.4142, C11283=1.12838, C15=1.5,C32=3.2,          & 
     &           C05642=0.5642,C79788=0.79788,C02=0.2,C14=1.4,            & 
     &           C37613=0.37613,C23=2./3.,                                & 
     &           CV1=-.122727278,CV2=.532770573,CV3=-.96284325,           & 
     &           CV4=.979895032 
      IV=int(V*TWOH+C15) 
      IF(A.LT.C02) THEN 
         IF(V.LE.TEN) THEN 
            VOIGTK=(H2TAB(IV)*A+H1TAB(IV))*A+H0TAB(IV) 
          ELSE 
            VOIGTK=C05642*A/(V*V) 
         END IF 
         RETURN 
      END IF 
      IF(A.LE.C14.AND.A+V.LE.C32) THEN 
         VV=V*V 
         HH1=H1TAB(IV)+H0TAB(IV)*C11283 
         HH2=H2TAB(IV)+HH1*C11283-H0TAB(IV) 
         HH3=(ONE-H2TAB(IV))*C37613-HH1*C23*VV+HH2*C11283 
         HH4=(THREE*HH3-HH1)*C37613+H0TAB(IV)*C23*VV*VV 
         VOIGTK=((((HH4*A+HH3)*A+HH2)*A+HH1)*A+H0TAB(IV))*                & 
     &          (((CV1*A+CV2)*A+CV3)*A+CV4) 
       ELSE 
         AA=A*A 
         VV=V*V 
         U=(AA+VV)*C14142 
         UU=U*U 
         VOIGTK=((((AA-TEN*VV)*AA*THREE+FIFTN*VV*VV)/UU+                  & 
     &          THREE*VV-AA)/UU+ONE)*A*C79788/U 
      END IF 
      RETURN 
      END FUNCTION VOIGTK 
! 
! 
!     ******************************************************************* 
! 
! 
 
      SUBROUTINE RTECD 
!     ================ 
! 
!     solution of the radiative transfer equation by Feautrier method 
!     for two continuum points 
!     used when one employs RTEDFE, ie. the DFE method for the 
!     transfer equation for the inner frequency points 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: D(3,3,MDEPTH),ANU(3,MDEPTH),AANU(MDEPTH),DDD(MDEPTH),   & 
     &            AA(3,3),BB(3,3),CC(3,3),VL(3),AMU(3),WTMU(3),           & 
     &            DT(MDEPTH),TAU(MDEPTH),                                 & 
     &            RDD(MDEPTH),FKK(MDEPTH),ST0(MDEPTH),SS0(MDEPTH),        & 
     &            RINT(MDEPTH,MMU) 
      REAL(DP), PARAMETER :: UN=1.,HALF=0.5,THIRD=UN/3.,QUART=UN/4.,      & 
     &                       SIXTH=UN/6.,TAUREF=21./3. 
      DATA AMU  /.887298334620742,.5,.112701665379258/ 
      DATA WTMU /.277777777777778,.444444444444444,.277777777777778/ 
! 
      NMU=3 
      ND1=ND-1 
! 
!     loop over two continuum frequencies 
! 
      FREQS: DO IJ=1,2 
      TAUMIN=CH(IJ,1)/DENS(1)*DM(1)*HALF 
      TAU(1)=TAUMIN 
      DO I=1,ND1 
         DT(I)=(DM(I+1)-DM(I))*(CH(IJ,I+1)/DENS(I+1)+CH(IJ,I)/DENS(I))*   & 
     &         HALF 
         ST0(I)=ET(IJ,I)/CH(IJ,I) 
         SS0(I)=-SC(IJ,I)/CH(IJ,I) 
         TAU(I+1)=TAU(I)+DT(I) 
         IF(TAU(I).LE.TAUREF.AND.TAU(I+1).GT.TAUREF) IREF=I 
!      if(ij.eq.1.and.i.eq.1) write(*,"('rtecd',i4,f8.2,i4,1p4e10.2)") 
!    * ij,3.e18/freq(ij),i,ch(ij,i),et(ij,i),sc(ij,i),dt(i) 
      END DO 
      ST0(ND)=ET(IJ,ND)/CH(IJ,ND) 
      SS0(ND)=-SC(IJ,ND)/CH(IJ,ND) 
      FR=FREQ(IJ) 
      BNU=BN*(FR*1.E-15)**3 
      PLAND=BNU/(EXP(HK*FR/TEMP(ND  ))-UN) 
      DPLAN=BNU/(EXP(HK*FR/TEMP(ND-1))-UN) 
      DPLAN=(PLAND-DPLAN)/DT(ND1) 
! 
!   +++++++++++++++++++++++++++++++++++++++++ 
!   FIRST PART  -  VARIABLE EDDINGTON FACTORS 
!   +++++++++++++++++++++++++++++++++++++++++ 
! 
!   Allowance for wind blanketing 
! 
      ALB1=0. 
      DO I=1,NMU 
! 
!   ************************ 
!   UPPER BOUNDARY CONDITION 
!   ************************ 
! 
         ID=1 
         DTP1=DT(1) 
         Q0=0. 
         P0=0. 
! 
!        allowance for non-zero optical depth at the first depth point 
! 
         TAMM=TAUMIN/AMU(I) 
         IF(TAMM.GT.0.01) THEN 
            P0=UN-EXP(-TAMM) 
          ELSE 
            P0=TAMM*(UN-HALF*TAMM*(UN-TAMM*THIRD*(UN-QUART*TAMM))) 
         END IF 
         EX=UN-P0 
         Q0=Q0+P0*AMU(I)*WTMU(I) 
! 
         DIV=DTP1/AMU(I)*THIRD 
         VL(I)=DIV*(ST0(ID)+HALF*ST0(ID+1))+ST0(ID)*P0 
         DO J=1,NMU 
            BB(I,J)=SS0(ID)*WTMU(J)*(DIV+P0)-ALB1*WTMU(J) 
            CC(I,J)=-HALF*DIV*SS0(ID+1)*WTMU(J) 
         END DO 
         BB(I,I)=BB(I,I)+AMU(I)/DTP1+UN+DIV 
         CC(I,I)=CC(I,I)+AMU(I)/DTP1-HALF*DIV 
         ANU(I,ID)=0. 
      END DO 
! 
!     Matrix inversion: instead of calling MATINV, a very fast inlined 
!     routine MINV3 for a specific 3 x 3 matrix inversion 
! 
!     CALL MATINV(BB,NMU,3) 
! 
!     ****************************** 
      BB(2,1)=BB(2,1)/BB(1,1) 
      BB(2,2)=BB(2,2)-BB(2,1)*BB(1,2) 
      BB(2,3)=BB(2,3)-BB(2,1)*BB(1,3) 
      BB(3,1)=BB(3,1)/BB(1,1) 
      BB(3,2)=(BB(3,2)-BB(3,1)*BB(1,2))/BB(2,2) 
      BB(3,3)=BB(3,3)-BB(3,1)*BB(1,3)-BB(3,2)*BB(2,3) 
! 
      BB(3,2)=-BB(3,2) 
      BB(3,1)=-BB(3,1)-BB(3,2)*BB(2,1) 
      BB(2,1)=-BB(2,1) 
! 
      BB(3,3)=UN/BB(3,3) 
      BB(2,3)=-BB(2,3)*BB(3,3)/BB(2,2) 
      BB(2,2)=UN/BB(2,2) 
      BB(1,3)=-(BB(1,2)*BB(2,3)+BB(1,3)*BB(3,3))/BB(1,1) 
      BB(1,2)=-BB(1,2)*BB(2,2)/BB(1,1) 
      BB(1,1)=UN/BB(1,1) 
! 
      BB(1,1)=BB(1,1)+BB(1,2)*BB(2,1)+BB(1,3)*BB(3,1) 
      BB(1,2)=BB(1,2)+BB(1,3)*BB(3,2) 
      BB(2,1)=BB(2,2)*BB(2,1)+BB(2,3)*BB(3,1) 
      BB(2,2)=BB(2,2)+BB(2,3)*BB(3,2) 
      BB(3,1)=BB(3,3)*BB(3,1) 
      BB(3,2)=BB(3,3)*BB(3,2) 
!     ****************************** 
! 
      DO I=1,NMU 
         DO J=1,NMU 
            S=0. 
            DO K=1,NMU 
               S=S+BB(I,K)*CC(K,J) 
            END DO 
            D(I,J,ID)=S 
            ANU(I,1)=ANU(I,1)+BB(I,J)*VL(J) 
         END DO 
      END DO 
! 
!   ******************* 
!   NORMAL DEPTH POINTS 
!   ******************* 
! 
      DO ID=2,ND1 
         DTM1=DTP1 
         DTP1=DT(ID) 
         DT0=HALF*(DTM1+DTP1) 
         AL=UN/DTM1/DT0 
         GA=UN/DTP1/DT0 
         BE=AL+GA 
         A=(UN-HALF*AL*DTP1*DTP1)*SIXTH 
         C=(UN-HALF*GA*DTM1*DTM1)*SIXTH 
         B=UN-A-C 
         VL0=A*ST0(ID-1)+B*ST0(ID)+C*ST0(ID+1) 
         DO I=1,NMU 
            DO J=1,NMU 
               AA(I,J)=-A*SS0(ID-1)*WTMU(J) 
               CC(I,J)=-C*SS0(ID+1)*WTMU(J) 
               BB(I,J)=B*SS0(ID)*WTMU(J) 
            END DO 
         END DO 
         DO I=1,NMU 
            DIV=AMU(I)**2 
            VL(I)=VL0 
            AA(I,I)=AA(I,I)+DIV*AL-A 
            CC(I,I)=CC(I,I)+DIV*GA-C 
            BB(I,I)=BB(I,I)+DIV*BE+B 
         END DO 
         DO I=1,NMU 
            S1=0. 
            DO J=1,NMU 
               S=0. 
               S1=S1+AA(I,J)*ANU(J,ID-1) 
               DO K=1,NMU 
                  S=S+AA(I,K)*D(K,J,ID-1) 
               END DO 
               BB(I,J)=BB(I,J)-S 
            END DO 
            VL(I)=VL(I)+S1 
         END DO 
! 
!     Matrix inversion: instead of calling MATINV, a very fast inlined 
!     routine MINV3 for a specific 3 x 3 matrix inversion 
! 
!     CALL MATINV(BB,NMU,3) 
! 
!     ****************************** 
      BB(2,1)=BB(2,1)/BB(1,1) 
      BB(2,2)=BB(2,2)-BB(2,1)*BB(1,2) 
      BB(2,3)=BB(2,3)-BB(2,1)*BB(1,3) 
      BB(3,1)=BB(3,1)/BB(1,1) 
      BB(3,2)=(BB(3,2)-BB(3,1)*BB(1,2))/BB(2,2) 
      BB(3,3)=BB(3,3)-BB(3,1)*BB(1,3)-BB(3,2)*BB(2,3) 
! 
      BB(3,2)=-BB(3,2) 
      BB(3,1)=-BB(3,1)-BB(3,2)*BB(2,1) 
      BB(2,1)=-BB(2,1) 
! 
      BB(3,3)=UN/BB(3,3) 
      BB(2,3)=-BB(2,3)*BB(3,3)/BB(2,2) 
      BB(2,2)=UN/BB(2,2) 
      BB(1,3)=-(BB(1,2)*BB(2,3)+BB(1,3)*BB(3,3))/BB(1,1) 
      BB(1,2)=-BB(1,2)*BB(2,2)/BB(1,1) 
      BB(1,1)=UN/BB(1,1) 
! 
      BB(1,1)=BB(1,1)+BB(1,2)*BB(2,1)+BB(1,3)*BB(3,1) 
      BB(1,2)=BB(1,2)+BB(1,3)*BB(3,2) 
      BB(2,1)=BB(2,2)*BB(2,1)+BB(2,3)*BB(3,1) 
      BB(2,2)=BB(2,2)+BB(2,3)*BB(3,2) 
      BB(3,1)=BB(3,3)*BB(3,1) 
      BB(3,2)=BB(3,3)*BB(3,2) 
!     ****************************** 
! 
         DO I=1,NMU 
            ANU(I,ID)=0. 
            DO J=1,NMU 
               S=0. 
               DO K=1,NMU 
                  S=S+BB(I,K)*CC(K,J) 
               END DO 
               D(I,J,ID)=S 
               ANU(I,ID)=ANU(I,ID)+BB(I,J)*VL(J) 
            END DO 
         END DO 
      END DO 
! 
!   ************ 
!   LOWER BOUNDARY CONDITION (SA) 
!   ************ 
! 
      ID=ND 
      DO I=1,NMU 
         AA(I,I)=AMU(I)/DTP1 
         VL(I)=PLAND+AMU(I)*DPLAN+AA(I,I)*ANU(I,ID-1) 
         DO J=1,NMU 
            BB(I,J)=-AA(I,I)*D(I,J,ID-1) 
         END DO 
         BB(I,I)=BB(I,I)+AA(I,I)+UN 
      END DO 
!
!     DISKS
!
      IF(INMOD.EQ.2) THEN
         B=DTP1*HALF
         DO I=1,NMU
            BI=B/AMU(I)
            VL(I)=ST0(ID)*BI
            DO J=1,NMU
               BB(I,J)=BI*SS0(ID)*WTMU(J)
            END DO
            AA(I,I)=AMU(I)/DTP1
            BB(I,I)=BB(I,I)+AMU(I)/DTP1+BI
         END DO
         DO I=1,NMU
            S1=0.
            DO J=1,NMU
               S=0.
               S1=S1+AA(I,J)*ANU(J,ID-1)
               DO K=1,NMU
                  S=S+AA(I,K)*D(K,J,ID-1)
               END DO
               BB(I,J)=BB(I,J)-S
            END DO
            VL(I)=VL(I)+S1
         END DO
      END IF

! 
!     Matrix inversion: instead of calling MATINV, a very fast inlined 
!     routine MINV3 for a specific 3 x 3 matrix inversion 
! 
!     CALL MATINV(BB,NMU,3) 
! 
!     ****************************** 
      BB(2,1)=BB(2,1)/BB(1,1) 
      BB(2,2)=BB(2,2)-BB(2,1)*BB(1,2) 
      BB(2,3)=BB(2,3)-BB(2,1)*BB(1,3) 
      BB(3,1)=BB(3,1)/BB(1,1) 
      BB(3,2)=(BB(3,2)-BB(3,1)*BB(1,2))/BB(2,2) 
      BB(3,3)=BB(3,3)-BB(3,1)*BB(1,3)-BB(3,2)*BB(2,3) 
! 
      BB(3,2)=-BB(3,2) 
      BB(3,1)=-BB(3,1)-BB(3,2)*BB(2,1) 
      BB(2,1)=-BB(2,1) 
! 
      BB(3,3)=UN/BB(3,3) 
      BB(2,3)=-BB(2,3)*BB(3,3)/BB(2,2) 
      BB(2,2)=UN/BB(2,2) 
      BB(1,3)=-(BB(1,2)*BB(2,3)+BB(1,3)*BB(3,3))/BB(1,1) 
      BB(1,2)=-BB(1,2)*BB(2,2)/BB(1,1) 
      BB(1,1)=UN/BB(1,1) 
! 
      BB(1,1)=BB(1,1)+BB(1,2)*BB(2,1)+BB(1,3)*BB(3,1) 
      BB(1,2)=BB(1,2)+BB(1,3)*BB(3,2) 
      BB(2,1)=BB(2,2)*BB(2,1)+BB(2,3)*BB(3,1) 
      BB(2,2)=BB(2,2)+BB(2,3)*BB(3,2) 
      BB(3,1)=BB(3,3)*BB(3,1) 
      BB(3,2)=BB(3,3)*BB(3,2) 
!     ****************************** 
! 
      DO I=1,NMU 
         ANU(I,ID)=0. 
         DO J=1,NMU 
            D(I,J,ID)=0. 
            ANU(I,ID)=ANU(I,ID)+BB(I,J)*VL(J) 
         END DO 
      END DO 
! 
!   ************ 
!   BACKSOLUTION 
!   ************ 
! 
      DO ID=ND-1,1,-1 
         DO I=1,NMU 
            DO J=1,NMU 
               ANU(I,ID)=ANU(I,ID)+D(I,J,ID)*ANU(J,ID+1) 
            END DO 
         END DO 
         AJ=0. 
         AK=0. 
         DO I=1,NMU 
            DIV=WTMU(I)*ANU(I,ID) 
            AJ=AJ+DIV 
            AK=AK+DIV*AMU(I)**2 
         END DO 
         FKK(ID)=AK/AJ 
      END DO 
! 
!     surface Eddington actor 
! 
      AH=0. 
      DO I=1,NMU 
         AH=AH+WTMU(I)*AMU(I)*ANU(I,1) 
      END DO 
      FH=AH/AJ-HALF*ALB1 
! 
      FKK(ND)=THIRD 
! 
! 
!   +++++++++++++++++++++++++++++++++++++++++ 
!   SECOND PART  -  DETERMINATION OF THE MEAN INTENSITIES 
!   RECALCULATION OF THE TRANSFER EQUATION WITH GIVEN EDDINGTON FACTORS 
!   +++++++++++++++++++++++++++++++++++++++++ 
! 
      DTP1=DT(1) 
      DIV=DTP1*THIRD 
      BBB=FKK(1)/DTP1+FH+DIV+SS0(1)*(DIV+Q0) 
      CCC=FKK(2)/DTP1-HALF*DIV*(UN+SS0(2)) 
      VLL=DIV*(ST0(1)+HALF*ST0(2))+ST0(1)*Q0 
      AANU(1)=VLL/BBB 
      DDD(1)=CCC/BBB 
      DO ID=2,ND1 
         DTM1=DTP1 
         DTP1=DT(ID) 
         DT0=HALF*(DTP1+DTM1) 
         AL=UN/DTM1/DT0 
         GA=UN/DTP1/DT0 
         A=(UN-HALF*DTP1*DTP1*AL)*SIXTH 
         C=(UN-HALF*DTM1*DTM1*GA)*SIXTH 
         AAA=AL*FKK(ID-1)-A*(UN+SS0(ID-1)) 
         CCC=GA*FKK(ID+1)-C*(UN+SS0(ID+1)) 
         BBB=(AL+GA)*FKK(ID)+(UN-A-C)*(UN+SS0(ID)) 
         VLL=A*ST0(ID-1)+C*ST0(ID+1)+(UN-A-C)*ST0(ID) 
         BBB=BBB-AAA*DDD(ID-1) 
         DDD(ID)=CCC/BBB 
         AANU(ID)=(VLL+AAA*AANU(ID-1))/BBB 
      END DO 
      BBB=FKK(ND)/DTP1+HALF 
      AAA=FKK(ND1)/DTP1 
      BBB=BBB-AAA*DDD(ND1) 
      VLL=HALF*PLAND+DPLAN*THIRD 
      RDD(ND)=(VLL+AAA*AANU(ND1))/BBB 
      DO IID=1,ND1 
         ID=ND-IID 
         RDD(ID)=AANU(ID)+DDD(ID)*RDD(ID+1) 
      END DO 
      FLUX(IJ)=FH*RDD(1) 
! 
      if(ij.eq.1) then 
         do id=1,nd 
            scc1(id)=-rdd(id)*ss0(id)*ch(1,id) 
         end do 
       else 
         do id=1,nd 
            scc2(id)=-rdd(id)*ss0(id)*ch(2,id) 
         end do 
      end if 
! 
!     if needed (if iprin.ge.3), output of interesting physical 
!     quantities at the monochromatic optical depth  tau(nu)=2/3 
! 
      IF(IPRIN.ge.3) THEN 
      T0=LOG(TAU(IREF+1)/TAU(IREF)) 
      X0=LOG(TAU(IREF+1)/TAUREF)/T0 
      X1=LOG(TAUREF/TAU(IREF))/T0 
      DMREF=EXP(LOG(DM(IREF))*X0+LOG(DM(IREF+1))*X1) 
      TREF=EXP(LOG(TEMP(IREF))*X0+LOG(TEMP(IREF+1))*X1) 
      STREF=EXP(LOG(ST0(IREF))*X0+LOG(ST0(IREF+1))*X1) 
      SCREF=EXP(LOG(-SS0(IREF))*X0+LOG(-SS0(IREF+1))*X1) 
      SSREF=EXP(LOG(-SS0(IREF)*RDD(IREF))*X0+                             & 
     &           LOG(-SS0(IREF+1)*RDD(IREF+1))*X1) 
      SREF=STREF+SSREF 
      ALM=2.997925E18/FREQ(IJ) 
      WRITE(96,"(I3,F10.3,I4,1PE10.3,0PF10.1,1X,1P3E10.3,E11.3)")         & 
     &   IJ,ALM,IREF,DMREF,TREF,SCREF,STREF,SSREF,SREF 
      END IF 
! 
!   ******************************************************************** 
! 
!   THIRD PART  -  DETERMINATION OF THE SPECIFIC INTENSITIES 
!   RECALCULATION OF THE TRANSFER EQUATION WITH GIVEN SOURCE FUNCTION 
! 
      IF(IFLUX.EQ.0) CYCLE FREQS 
      ANGLES: DO IMU=1,NMU0 
         IF(IFLUX.EQ.0) EXIT ANGLES 
         ANX=ANGL(IMU) 
         DTP1=DT(1) 
         DIV=DTP1*THIRD/ANX 
! 
         TAMM=TAUMIN/ANX 
         IF(TAMM.LT.0.01) THEN 
            P0=TAMM*(UN-HALF*TAMM*(UN-TAMM*THIRD*(UN-QUART*TAMM))) 
          ELSE 
            P0=UN-EXP(-TAMM) 
         END IF 
! 
         BBB=ANX/DTP1+UN+DIV 
         CCC=ANX/DTP1-HALF*DIV 
         VLL=(DIV+P0)*(ST0(1)-SS0(1)*RDD(1))                              & 
     &       +HALF*DIV*(ST0(2)-SS0(2)*RDD(2)) 
         AANU(1)=VLL/BBB 
         DDD(1)=CCC/BBB 
         DIV=ANX*ANX 
         DO ID=2,ND1 
            DTM1=DT(ID-1) 
            DTP1=DT(ID) 
            DT0=HALF*(DTP1+DTM1) 
            AL=UN/DTM1/DT0 
            GA=UN/DTP1/DT0 
            A=(UN-HALF*DTP1*DTP1*AL)*SIXTH 
            C=(UN-HALF*DTM1*DTM1*GA)*SIXTH 
            AAA=DIV*AL-A 
            CCC=DIV*GA-C 
            BBB=DIV*(AL+GA)+UN-A-C 
            VLL=A*(ST0(ID-1)-SS0(ID-1)*RDD(ID-1))+                        & 
     &          C*(ST0(ID+1)-SS0(ID+1)*RDD(ID+1))+                        & 
     &          (UN-A-C)*(ST0(ID)-SS0(ID)*RDD(ID)) 
            BBB=BBB-AAA*DDD(ID-1) 
            DDD(ID)=CCC/BBB 
            AANU(ID)=(VLL+AAA*AANU(ID-1))/BBB 
         END DO 
! 
!        Lower boundary condition 
! 
         AAA=ANX/DTP1 
         BBB=AAA+UN 
         VLL=PLAND+ANX*DPLAN 
! 
         RINT(ND,IMU)=(VLL+AAA*AANU(ND1))/(BBB-AAA*DDD(ND1)) 
         DO IID=1,ND1 
            ID=ND-IID 
            RINT(ID,IMU)=AANU(ID)+DDD(ID)*RINT(ID+1,IMU) 
         END DO 
      END DO ANGLES 
! 
      IF(IFLUX.GE.1) THEN 
         FLX=0. 
         DO IMU=1,NMU0 
            RINT(1,IMU)=RINT(1,IMU)/HALF 
            FLX=FLX+ANGL(IMU)*WANGL(IMU)*RINT(1,IMU) 
         END DO 
         FLX=FLX*HALF 
! 
!        output of emergent specific intensities in continuum to Unit 18 
! 
         WRITE(18,"(f10.3,1pe15.5/(1P5E15.5))")                           & 
     &   WLAM(IJ),FLX,(RINT(1,IMU),IMU=1,NMU0) 
      END IF 
      END DO FREQS 
! 
!     call rtedfe for the internal points 
! 
      CALL RTEDFE 
! 
      RETURN 
      END SUBROUTINE RTECD 
! 
! 
!     ******************************************************************* 
! 
! 
 
      SUBROUTINE RTEDFE 
!     ================= 
! 
!     Solution of the radiative transfer equation - frequency by 
!     frequency - for the known source function. 
! 
!     The numerical method used: 
!     Discontinuous Finite Element (DFE) method 
!     Castor, Dykema, Klein, 1992, ApJ 387, 561. 
! 
!     Input through blank COMMON block: 
!      CH     - two-dimensional array  absorption coefficient (frequency, 
!               depth) 
!      ET     - emission coefficient (frequency, depth) 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: ONE=1.,TWO=2.,HALF=0.5,TAUREF=2./3. 
      REAL(DP) :: DT(MDEPTH),ST0(MDEPTH),AB0(MDEPTH),DELDM(MDEPTH),       & 
     &            dtau(mdepth),rip(mdepth),rim(mdepth),riup(mdepth),      & 
     &            AMU(3),WTMU(3),RINT1(MMU),                              & 
     &            AMUI(MMU),AMUW(MMU),TAU(MDEPTH),SS0(MDEPTH) 
! 
!     angle points (AMU) and angular integration weights (WTMU) 
! 
      DATA AMU  /.887298334620742,.5,.112701665379258/ 
      DATA WTMU /.277777777777778,.444444444444444,.277777777777778/ 
! 
      DO I=1,ND-1 
         DELDM(I)=HALF*(DM(I+1)-DM(I)) 
      END DO 
! 
!     angle points 
! 
      IF(IFLUX.EQ.0) THEN 
         NMUS=NMU 
         do i=1,nmu 
            amui(i)=amu(i) 
            amuw(i)=amu(i)*wtmu(i) 
         end do 
       ELSE IF(IFLUX.EQ.1) THEN 
         NMUS=NMU0 
         do i=1,nmus 
            amui(i)=angl(i) 
            amuw(i)=angl(i)*wangl(i) 
         end do 
      END IF 
! 
!     overall loop over frequencies 
! 
      DO IJ=1,NFREQ 
      FR=FREQ(IJ) 
! 
!     total source function 
! 
      DO ID=1,ND 
         AB0(ID)=CH(IJ,ID) 
         SCT=FRX1(IJ)*SCC2(ID)+FRX2(IJ)*SCC1(ID) 
         ST0(ID)=(ET(IJ,ID)+SCT)/AB0(ID) 
         SS0(ID)=-SCT/AB0(ID) 
      END DO 
      AH=0. 
! 
!     optical depth scale 
! 
      TAU(1)=0. 
      IREF=1 
      DO ID=1,ND-1 
         DT(ID)=DELDM(ID)*(AB0(ID+1)/DENS(ID+1)+AB0(ID)/DENS(ID)) 
         TAU(ID+1)=TAU(ID)+DT(ID) 
         IF(TAU(ID).LE.TAUREF.AND.TAU(ID+1).GT.TAUREF) IREF=ID 
      END DO 
      IREFD(IJ)=IREF 
! 
!     quantities for the lower boundary condition 
! 
      FR15=FR*1.D-15 
      BNU=BN*FR15*FR15*FR15 
      PLAND=BNU/(EXP(HK*FR/TEMP(ND))-ONE) 
      DPLAN=BNU/(EXP(HK*FR/TEMP(ND-1))-ONE) 
      DPLAN=(PLAND-DPLAN)/DT(ND-1) 
! 
!     loop over angle poits 
! 
      DO I=1,NMUS 
         do id=1,nd-1 
            dtau(id)=dt(id)/amui(i) 
         enddo 
! 
!           outgoing intensity 
! 
            if(inmod.ne.2) rip(nd)=PLAND+AMUI(I)*DPLAN 
            id=nd-1 
            dt0=dtau(id) 
            dtaup1=dt0+one 
            dtau2=dt0*dt0 
            bb=two*dtaup1 
            cc=dt0*dtaup1 
            aa=dtau2+bb 
            rim(id+1)=(aa*rip(id+1)-cc*st0(id+1)+dt0*st0(id))/bb 
            do id=nd-1,1,-1 
               dt0=dtau(id) 
               dtaup1=dt0+one 
               dtau2=dt0*dt0 
               bb=two*dtaup1 
               cc=dt0*dtaup1 
               aa=one/(dtau2+bb) 
               rim(id)=(two*rim(id+1)+dt0*st0(id+1)+cc*st0(id))*aa 
               rip(id+1)=(bb*rim(id+1)+cc*st0(id+1)-dt0*st0(id))*aa 
            enddo 
            do id=2,nd-1 
               riup(id)=(rim(id)*dtau(id-1)+rip(id)*dtau(id))/            & 
     &                  (dtau(id-1)+dtau(id)) 
            enddo 
            riup(1)=rim(1) 
            riup(nd)=rip(nd) 
! 
         AH=AH+AMUW(I)*RIUP(1) 
         RINT1(I)=RIUP(1) 
         rint1(i)=max(rint1(i),1.e-40) 
! 
!     end of the loop over angle points 
! 
      END DO 
! 
      FLUX(IJ)=AH*HALF 
      if(iflux.ge.1) then 
! 
!     output of emergent specific intensities to Unit 10 (line points) 
!     or 18 (two continuum points) 
! 
      IF(IJ.GT.2) THEN 
         WRITE(10,"(1H ,f10.3,1pe15.5/(1P5E15.5))")                       & 
     &   WLAM(IJ),FLUX(IJ),(RINT1(IMU),IMU=1,NMUS) 
       ELSE 
         WRITE(18,"(1H ,f10.3,1pe15.5/(1P5E15.5))")                       & 
     &   WLAM(IJ),FLUX(IJ),(RINT1(IMU),IMU=1,NMUS) 
      END IF 
      end if 
! 
!     if needed (if iprin.ge.3), output of interesting physical 
!     quantities at the monochromatic optical depth  tau(nu)=2/3 
! 
      IF(IPRIN.GE.3) THEN 
         T0=LOG(TAU(IREF+1)/TAU(IREF)) 
         X0=LOG(TAU(IREF+1)/TAUREF)/T0 
         X1=LOG(TAUREF/TAU(IREF))/T0 
         DMREF=EXP(LOG(DM(IREF))*X0+LOG(DM(IREF+1))*X1) 
         TREF=EXP(LOG(TEMP(IREF))*X0+LOG(TEMP(IREF+1))*X1) 
         STREF=EXP(LOG(ST0(IREF))*X0+LOG(ST0(IREF+1))*X1) 
         SSREF=EXP(LOG(-SS0(IREF))*X0+LOG(-SS0(IREF+1))*X1) 
         SREF=STREF+SSREF 
         ALM=2.997925E18/FREQ(IJ) 
         WRITE(96,"(I3,F10.3,I4,1PE10.3,0PF10.1,1X,1P3E10.3)")            & 
     &      IJ,ALM,IREF,DMREF,TREF,STREF,SSREF,SREF 
      END IF 
! 
!     end of the loop over frequencies 
! 
      END DO 
      RETURN 
      END SUBROUTINE RTEDFE 
! 
! 
!    ******************************************************************* 
! 
! 
      SUBROUTINE PARTF(IAT,IZI,T,ANE,XMAXN,U) 
!     ======================================= 
! 
!     Partition functions 
!     The standard evaluation is for hydrogen through zinc, for 
!     neutrals and first four ionization degrees. 
!     Basically after Traving, Baschek, and Holweger, Abhand. Hamburg. 
!     Sternwarte. Band VIII, Nr. 1 (1966) 
! 
!     For higher atomic numbers  modified Kurucz routine PFSAHA, 
!     called PFHEAV here is used. The routine was provided by 
!     Charles Proffitt. 
! 
!     The routine calls special procedures for Fe and Ni; or 
!     the values based on the tabulated Opacity Project ionization 
!     fractions 
! 
!     Input: 
!      IAT   - atomic number 
!      IZI   - ionic charge (=1 for neutrals, =2 for once ionized, etc) 
!      T     - temperature 
!      ANE   - electron density 
!      XMAXN - principal quantum number of the last bound level 
! 
!     Output: 
!      U     - partition function 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      INTEGER, PARAMETER :: NIONS=123, NSS=222 
      REAL(DP), PARAMETER :: UN=1., HALF=0.5, TWO=2., TRHA=1.5,           & 
     &                       THIRD=UN/3., SIXTH=UN/6. 
      REAL(DP) ::                                                         & 
     &       AHH( 6),  ALB(12),  AB (11),  AC (19),  AN (30),  AO (49),   & 
     &       AF (34),  ANN(23),  ANA(19),  AMG(15),  AAL(17),  ASI(23),   & 
     &       AP (19),  AS (29),  ACL(28),  AAR(25),  AK (30),  ACA(17),   & 
     &       ASC(24),  ATI(33),  AV (33),  ACR(29),  AMN(28),  AFE(35),   & 
     &       ACO(29),  ANI(23),  ACU(20),  AZN(18),                       & 
     &       GHH( 6),  GLB(12),  GB (11),  GC (19),  GN (30),  GO (49),   & 
     &       GF (34),  GNN(23),  GNA(19),  GMG(15),  GAL(17),  GSI(23),   & 
     &       GP (19),  GS (29),  GCL(28),  GAR(25),  GK (30),  GCA(17),   & 
     &       GSC(24),  GTI(33),  GV (33),  GCR(29),  GMN(28),  GFE(35),   & 
     &       GCO(29),  GNI(23),  GCU(20),  GZN(18),                       & 
     &       XL(222),                                                     & 
     &       CH1(66),  CH2(72),  CH3(55),  CH4(29),  CHION(222),          & 
     &       ALF(678), GAM(678) 
      INTEGER :: INDEX0(5,30),                                            & 
     *           IS(123),                                                 & 
     *           IM(222),                                                 & 
     *           IGPR(222),                                               & 
     *           IG0(123),                                                & 
     *           IGLE(28) 
      INTEGER, SAVE :: INDEXS(123),INDEXM(222) 
! 
      EQUIVALENCE   ( AHH(1), ALF(  1)),( ALB(1), ALF(  7)),              & 
     &              ( AB (1), ALF( 19)),                                  & 
     &              ( AC (1), ALF( 30)),( AN (1), ALF( 49)),              & 
     &              ( AO (1), ALF( 79)),( AF (1), ALF(128)),              & 
     &              ( ANN(1), ALF(162)),( ANA(1), ALF(185)),              & 
     &              ( AMG(1), ALF(204)),( AAL(1), ALF(219)),              & 
     &              ( ASI(1), ALF(236)),( AP (1), ALF(259)),              & 
     &              ( AS (1), ALF(278)),( ACL(1), ALF(307)),              & 
     &              ( AAR(1), ALF(335)),( AK (1), ALF(360)),              & 
     &              ( ACA(1), ALF(390)),( ASC(1), ALF(407)),              & 
     &              ( ATI(1), ALF(431)),( AV (1), ALF(464)),              & 
     &              ( ACR(1), ALF(497)),( AMN(1), ALF(526)),              & 
     &              ( AFE(1), ALF(554)),( ACO(1), ALF(589)),              & 
     &              ( ANI(1), ALF(618)),( ACU(1), ALF(641)),              & 
     &              ( AZN(1), ALF(661)) 
      EQUIVALENCE   ( GHH(1), GAM(  1)),( GLB(1), GAM(  7)),              & 
     &              ( GB (1), GAM( 19)),                                  & 
     &              ( GC (1), GAM( 30)),( GN (1), GAM( 49)),              & 
     &              ( GO (1), GAM( 79)),( GF (1), GAM(128)),              & 
     &              ( GNN(1), GAM(162)),( GNA(1), GAM(185)),              & 
     &              ( GMG(1), GAM(204)),( GAL(1), GAM(219)),              & 
     &              ( GSI(1), GAM(236)),( GP (1), GAM(259)),              & 
     &              ( GS (1), GAM(278)),( GCL(1), GAM(307)),              & 
     &              ( GAR(1), GAM(335)),( GK (1), GAM(360)),              & 
     &              ( GCA(1), GAM(390)),( GSC(1), GAM(407)),              & 
     &              ( GTI(1), GAM(431)),( GV (1), GAM(464)),              & 
     &              ( GCR(1), GAM(497)),( GMN(1), GAM(526)),              & 
     &              ( GFE(1), GAM(554)),( GCO(1), GAM(589)),              & 
     &              ( GNI(1), GAM(618)),( GCU(1), GAM(641)),              & 
     &              ( GZN(1), GAM(661)) 
      EQUIVALENCE   ( CH1(1), CHION(  1)),                                & 
     &              ( CH2(1), CHION( 67)),                                & 
     &              ( CH3(1), CHION(139)),                                & 
     &              ( CH4(1), CHION(194))                                 & 
! 
      DATA IGLE/2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,                      & 
     &          10,21,28,25,6,25,28,21,10,21/ 
! 
      DATA INDEX0   /   1,  -1,   0,   0,   0,                            & 
     &                  2,   3,  -1,   0,   0,                            & 
     &                  4,   5,  -2,  -1,   0,                            & 
     &                  6,   7,  -1,  -2,  -1,                            & 
     &                  8,   9,  10,  -1,  -2,                            & 
     &                 11,  12,  13,  14,  -1,                            & 
     &                 15,  16,  17,  18,  19,                            & 
     &                 20,  21,  22,  23,  24,                            & 
     &                 25,  26,  27,  28,  -6,                            & 
     &                 29,  30,  31,  32,  -9,                            & 
     &                 33,  34,  35,  36,  -4,                            & 
     &                 37,  38,  39,  40,  -9,                            & 
     &                 41,  42,  43,  44,  -6,                            & 
     &                 45,  46,  47,  48,  -1,                            & 
     &                 49,  50,  51,  52,  53,                            & 
     &                 54,  55,  56,  57,  58,                            & 
     &                 59,  60,  61,  62,  63,                            & 
     &                 64,  65,  66,  67,  68,                            & 
     &                 69,  70,  71,  72,  73,                            & 
     &                 74,  75,  76,  77,  -9,                            & 
     &                 78,  76,  80,  81,  82,                            & 
     &                 83,  84,  85,  86,  87,                            & 
     &                 88,  89,  90,  91,  92,                            & 
     &                 93,  94,  95,  96,  97,                            & 
     &                 98,  99, 100, 101, 102,                            & 
     &                103, 104, 105, 106, 107,                            & 
     &                108, 109, 110, 111, -25,                            & 
     &                112, 113, 114, 115,  -1,                            & 
     &                116, 117, 118, 119,  -1,                            & 
     &                120, 121, 122, 123,  -1                         / 
! 
      DATA IG0      /   2,                                                & 
     &                  1,   2,                                           & 
     &                  2,   1,                                           & 
     &                  1,   2,                                           & 
     &                  2,   1,   2,                                      & 
     &                  1,   2,   1,   2,                                 & 
     &                  4,   1,   2,   1,   2,                            & 
     &                  5,   4,   1,   2,   1,                            & 
     &                  4,   5,   4,   1,                                 & 
     &                  1,   4,   5,   4,                                 & 
     &                  2,   1,   4,   5,                                 & 
     &                  1,   2,   1,   4,                                 & 
     &                  2,   1,   2,   1,                                 & 
     &                  1,   2,   1,   2,                                 & 
     &                  4,   1,   2,   1,   2,                            & 
     &                  5,   4,   1,   2,   1,                            & 
     &                  4,   5,   4,   1,   2,                            & 
     &                  1,   4,   5,   4,   1,                            & 
     &                  2,   1,   4,   5,   4,                            & 
     &                  1,   2,   1,   4,                                 & 
     &                  4,   3,   4,   1,   4,                            & 
     &                  5,   4,   5,   4,   1,                            & 
     &                  4,   1,   4,   5,   4,                            & 
     &                  7,   6,   1,   4,   5,                            & 
     &                  6,   7,   6,   1,   4,                            & 
     &                  9,  10,   9,   6,   1,                            & 
     &                 10,   9,  10,  20,                                 & 
     &                  9,   6,   9,  28,                                 & 
     &                  2,   1,   6,  21,                                 & 
     &                  1,   2,   1,  10                              / 
! 
      DATA IS       /   1,                                                & 
     &                  1,   1,                                           & 
     &                  1,   1,                                           & 
     &                  2,   1,                                           & 
     &                  1,   2,   1,                                      & 
     &                  1,   2,   2,   1,                                 & 
     &                  2,   2,   3,   2,   1,                            & 
     &                  3,   4,   3,   5,   2,                            & 
     &                  2,   3,   4,   3,                                 & 
     &                  2,   2,   3,   2,                                 & 
     &                  1,   2,   2,   3,                                 & 
     &                  1,   1,   2,   2,                                 & 
     &                  2,   2,   1,   2,                                 & 
     &                  1,   2,   2,   1,                                 & 
     &                  2,   1,   1,   1,   1,                            & 
     &                  3,   2,   1,   2,   2,                            & 
     &                  2,   3,   2,   1,   1,                            & 
     &                  2,   2,   3,   1,   1,                            & 
     &                  1,   2,   3,   3,   2,                            & 
     &                  2,   1,   2,   2,                                 & 
     &                  3,   1,   1,   1,   1,                            & 
     &                  3,   2,   1,   1,   1,                            & 
     &                  2,   3,   1,   1,   1,                            & 
     &                  3,   2,   1,   1,   1,                            & 
     &                  3,   2,   1,   1,   1,                            & 
     &                  3,   2,   2,   1,   1,                            & 
     &                  4,   2,   1,   1,                                 & 
     &                  2,   2,   1,   1,                                 & 
     &                  3,   2,   1,   1,                                 & 
     &                  3,   3,   1,   1                              / 
! 
      DATA IM       /   2,                                                & 
     &                  2,   2,                                           & 
     &                  2,   2,                                           & 
     &                  3,   2,   3,                                      & 
     &                  3,   3,   2,   3,                                 & 
     &                  4,   3,   3,   3,   3,   3,                       & 
     &                  3,   3,   4,   3,   3,   4,   2,   3,   2,   3,   & 
     &                  4,   2,   2,   4,   2,   3,   3,   4,   4,   2,   & 
     &                  3,   4,   2,   2,   2,   3,   3,                  & 
     &                  3,   3,   4,   2,   2,                            & 
     &                  4,   2,   3,   2,   5,   2,   2,                  & 
     &                  2,   2,   3,   2,   4,   2,   2,   4,   2,        & 
     &                  2,   2,   2,   3,   2,   4,   2,   2,             & 
     &                  3,   3,   2,   2,   3,   2,                       & 
     &                  3,   2,   3,   2,   3,   2,   2,                  & 
     &                  5,   4,   4,   4,   3,   3,                       & 
     &                  3,   2,   4,   4,   3,   3,                       & 
     &                  4,   2,   2,   4,   2,   5,   4,   2,   3,   1,   & 
     &                  3,   2,   5,   2,   2,   4,   2,   4,   4,        & 
     &                  2,   2,   3,   2,   4,   2,   2,   4,   4,        & 
     &                  3,   2,   3,   3,   2,   3,                       & 
     &                  4,   2,   2,   4,   2,                            & 
     &                  3,   2,   3,   2,   2,   3,   2,                  & 
     &                  4,   3,   3,   5,   4,   2,   3,                  & 
     &                  6,   4,   3,   6,   3,   5,   4,   2,             & 
     &                  5,   3,   5,   4,   4,   4,   4,   4,             & 
     &                  3,   3,   3,   4,   4,   4,   4,   4,             & 
     &                  3,   2,   3,   4,   4,   4,   4,   4,             & 
     &                  4,   4,   3,   5,   3,   4,   4,   4,   4,        & 
     &                  5,   3,   3,   3,   5,   4,   5,   1,             & 
     &                  6,   3,   5,   3,   5,   1,                       & 
     &                  2,   3,   3,   4,   3,   4,   1,                  & 
     &                  2,   2,   2,   3,   3,   2,   3,   1          / 
! 
      DATA IGPR     /   2,                                                & 
     &                  4,   2,                                           & 
     &                  2,   4,                                           & 
     &                  4,  12,   2,                                      & 
     &                  2,   4,  12,   2,                                 & 
     &                 12,   2,  18,   4,  12,   2,                       & 
     &                 18,  10,  12,  24,   2,  18,   6,   4,  12,   2,   & 
     &                  8,  20,  12,  18,  10,   2,  10,  12,  24,  20,   & 
     &                  2,  18,   6,  18,  10,   4,  12,                  & 
     &                 18,  10,   8,  20,  12,                            & 
     &                 18,  10,   2,  10,  12,  24,  20,                  & 
     &                  8,   4,  18,  10,   8,  20,  12,  18,  10,        & 
     &                  2,   8,   4,  18,  10,   8,  20,  12,             & 
     &                  4,   2,   8,   4,  18,  10,                       & 
     &                  2,  18,   4,  12,   2,   8,   4,                  & 
     &                 12,   2,  18,   4,  12,   2,                       & 
     &                 18,  10,  12,   2,   4,   2,                       & 
     &                  8,  20,  12,  18,  10,  12,   2,  18,   4,  12,   & 
     &                 18,  10,   8,  20,  12,  18,  10,  12,   2,        & 
     &                  8,   4,  18,  10,   8,  20,  12,  18,  12,        & 
     &                  2,   8,   4,  18,  10,   2,                       & 
     &                  8,  20,  12,  18,  10,                            & 
     &                  4,  20,   2,   8,   4,  18,  10,                  & 
     &                 30,  42,  18,  20,   2,  12,  18,                  & 
     &                 56,  56,  28,  42,  10,  20,   2,  12,             & 
     &                 50,  70,  56,  72,  64,  42,  20,   2,             & 
     &                 12,  60,  40,  50,  18,  56,  42,  20,             & 
     &                 14,  10,  50,  12,  72,  50,  56,  42,             & 
     &                 60,  56,  40,  50,  18,  12,  72,  50,  56,        & 
     &                 42,  70,  42,  18,  56,  24,  50,  12,             & 
     &                 20,  56,  42,  18,  56,  50,                       & 
     &                  2,  30,  10,  20,  56,  42,  56,                  & 
     &                  4,   8,  12,   2,  30,  10,  20,  42          / 
! 
      DATA XL       /11.0,                                                & 
     &                8.0,12.0,                                           & 
     &                6.0, 6.0,                                           & 
     &                6.0, 4.0, 8.0,                                      & 
     &                9.0, 6.0, 4.0, 6.0,                                 & 
     &                6.0, 6.0, 5.0, 6.1, 5.0, 6.0,                       & 
     &                6.1, 4.0, 5.0, 3.9, 6.0, 5.0, 4.0, 6.0, 6.3, 6.0,   & 
     &                8.0, 6.0, 3.4, 6.0, 5.0, 3.9, 3.9, 6.0, 4.9, 4.0,   & 
     &                5.9, 5.0, 4.9, 4.0, 4.0, 6.0, 6.0,                  & 
     &                4.0, 4.0, 5.0, 4.0, 4.0,                            & 
     &                5.0, 4.0, 3.9, 4.0, 5.0, 5.0, 4.0,                  & 
     &                6.0, 6.0, 5.0, 4.0, 3.9, 4.0, 4.0, 5.0, 5.0,        & 
     &                7.0, 4.0, 4.0, 4.0, 4.0, 5.0, 5.0, 5.0,             & 
     &                7.0, 7.0, 5.0, 5.0, 5.0, 5.0,                       & 
     &                7.0, 4.0, 7.0, 4.0, 7.0, 5.0, 5.0,                  & 
     &                6.1, 5.9, 5.0, 5.0, 5.0, 7.0,                       & 
     &                5.0, 5.0, 5.0, 7.0, 8.6, 8.0,                       & 
     &                6.0, 5.0, 5.0, 5.0, 5.0, 3.5, 5.0,14.4, 5.0, 4.0,   & 
     &                6.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.2,        & 
     &                6.0, 6.0, 5.1, 5.0, 5.0, 5.0, 5.0, 5.0, 4.0,        & 
     &                7.0, 5.0, 5.0, 6.0, 6.0, 5.0,                       & 
     &                6.0, 5.0, 5.0, 3.6, 4.0,                            & 
     &                5.9, 6.0, 7.0, 5.0, 4.9, 5.0, 4.3,                  & 
     &                4.9, 4.9, 5.0, 5.0, 6.0, 4.6, 3.8,                  & 
     &                5.0, 4.7, 5.0, 5.0, 5.0, 5.0, 6.0, 4.8,             & 
     &                5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,11.2,             & 
     &                5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.2,             & 
     &                6.0, 5.0, 6.0, 7.0, 5.0, 5.0, 5.0, 5.0,             & 
     &                5.0, 5.0, 5.0, 5.0, 5.0, 6.0, 5.0, 3.6, 3.8,        & 
     &                5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 3.0,             & 
     &                5.4, 5.0, 9.0, 5.0, 5.0, 3.0,                       & 
     &                8.0, 6.0, 5.0, 7.0, 5.0, 5.0, 2.9,                  & 
     &                8.0, 5.0, 5.0, 8.0, 5.0, 5.0, 5.0, 2.8          / 
! 
! 
      DATA CH1      /  13.595 ,                                           & 
     &                 24.580 ,  54.403 ,                                 & 
     &                  5.390 , 75.619 ,                                  & 
     &                  9.320 ,  13.278 ,  18.206 ,                       & 
     &                  8.296 ,  25.149 ,  31.146 ,  37.920 ,             & 
     &                 11.256 ,  24.376 ,  30.868 ,  47.871 ,  55.873 ,   & 
     &                 64.476 ,                                           & 
     &                 14.529 ,  16.428 ,  29.593 ,  36.693 ,             & 
     &                 47.426 ,  55.765 ,  63.626 ,  77.450 ,  87.445 ,   & 
     &                 97.863 ,                                           & 
     &                 13.614 ,  16.938 ,  18.630 ,                       & 
     &                 35.108 ,  37.621 ,  40.461 ,  42.584 ,             & 
     &                 54.886 ,  63.733 ,  70.556 ,                       & 
     &                 77.394 ,  87.609 ,  97.077 , 103.911 , 106.116 ,   & 
     &                113.873 , 125.863 ,                                 & 
     &                 17.418 ,  20.009 ,  34.977 ,  39.204 ,  41.368 ,   & 
     &                 62.646 ,  65.774 ,  69.282 ,  71.882 ,             & 
     &                 87.139 ,  97.852 , 106.089 ,                       & 
     &                 21.559 ,  21.656 ,  41.071 ,  44.274 ,             & 
     &                 63.729 ,  68.806 ,  71.434 ,  97.162 , 100.917 / 
      DATA CH2      /   5.138 ,  47.290 ,  47.459 ,  71.647 ,  75.504 ,   & 
     &                 98.880 , 104.778 , 107.864 ,                       & 
     &                  7.644 ,  15.031 ,  80.117 ,  80.393 ,             & 
     &                109.294 , 113.799 ,                                 & 
     &                  5.984 ,  10.634 ,  18.823 ,  25.496 ,             & 
     &                 28.441 , 119.957 , 120.383 ,                       & 
     &                  8.149 ,  16.339 ,  22.894 ,                       & 
     &                 33.459 ,  42.333 ,  45.130 ,                       & 
     &                 10.474 ,  11.585 ,  19.720 ,                       & 
     &                 30.156 ,  51.354 ,  65.007 ,                       & 
     &                 10.357 ,  12.200 ,  13.401 ,  23.405 ,  24.807 ,   & 
     &                 35.047 ,  47.292 ,  57.681 ,  72.474 ,  85.701 ,   & 
     &                 13.014 ,  14.458 ,  23.798 ,  26.041 ,  27.501 ,   & 
     &                 39.904 ,  41.610 ,  53.450 ,  67.801 ,             & 
     &                 15.755 ,  15.933 ,  27.619 ,  29.355 ,             & 
     &                 40.899 ,  42.407 ,  45.234 ,  59.793 ,  75.002 ,   & 
     &                  4.339 ,  31.810 ,  32.079 ,                       & 
     &                 45.738 ,  47.768 ,  50.515 ,                       & 
     &                 60.897 ,  63.890 ,  65.849 ,  82.799 ,  85.150 / 
      DATA CH3      /   6.111 ,   7.808 ,  11.868 ,                       & 
     &                 51.207 ,  51.596 ,  67.181 ,  69.536 ,             & 
     &                  6.538 ,   7.147 ,   8.042 ,                       & 
     &                 12.891 ,  24.752 ,  74.090 ,  91.847 ,             & 
     &                  6.818 ,  6.953 ,  7.411 ,                         & 
     &                 13.635 ,  14.685 ,  28.137 ,  43.236 , 100.083 ,   & 
     &                  6.738 ,   7.101 ,  14.205 ,  15.670 ,  16.277 ,   & 
     &                 29.748 ,  48.464 ,  65.198 ,                       & 
     &                  6.763 ,   8.285 ,   9.221 ,                       & 
     &                 16.493 ,  18.662 ,  30.950 ,  49.580 ,  73.093 ,   & 
     &                  7.432 ,   8.606 ,   9.240 ,  15.636 ,  18.963 ,   & 
     &                 33.690 ,  53.001 ,  76.006 ,                       & 
     &                  7.896 ,   8.195 ,   8.927 ,  16.178 ,  18.662 ,   & 
     &                 30.640 ,  34.607 ,  56.001 ,  79.001           / 
      DATA CH4      /   7.863 ,   8.378 ,   9.160 ,   9.519 ,             & 
     &                 17.052 ,  18.958 ,  33.491 ,  53.001 ,             & 
     &                  7.633 ,   8.793 ,  18.147 ,  20.233 ,  35.161 ,   & 
     &                 56.025 ,                                           & 
     &                  7.724 ,  10.532 ,  10.980 ,                       & 
     &                 20.286 ,  27.985 ,  36.826 ,  61.975 ,             & 
     &                  9.391 ,  17.503 ,  17.166 ,                       & 
     &                 17.959 ,  27.757 ,  28.310 ,  39.701 ,  65.074 / 
! 
      DATA AHH      /  20.4976, 747.5023,                                 & 
     &                 28.1703, 527.8296,  22.2809, 987.7189          / 
      DATA GHH      /  10.853 ,  13.342 ,                                 & 
     &                 21.170 ,  24.125 ,  43.708 ,  53.542           / 
! 
      DATA ALB      /   8.4915,  97.5015,  23.3299, 192.6701,             & 
     &                  9.1849,  32.9263, 183.8887,  19.9563,  88.0437,   & 
     &                  6.0478,  35.9723, 233.9798                    / 
      DATA GLB      /   2.022 ,   4.604 ,  62.032 ,  72.624 ,             & 
     &                  2.735 ,   6.774 ,   8.569 ,  10.750 ,  11.672 ,   & 
     &                  3.967 ,  12.758 ,  16.692                     / 
! 
      DATA AB       /   4.0086,  19.6741, 402.3110,                       & 
     &                  9.7257,  30.9262, 186.3466,  44.1629,  60.8371,   & 
     &                  6.0084,  23.5767,  76.4149                    / 
      DATA GB       /   0.002 ,   3.971 ,   7.882 ,                       & 
     &                  4.720 ,  13.477 ,  22.103 ,  23.056 ,  24.734 ,   & 
     &                  6.000 ,  24.540 ,  32.300                     / 
! 
      DATA AC       /   8.0158,   5.8833,  33.7521, 595.3432,             & 
     &                  4.0003,  17.0841,  82.9154,                       & 
     &                 15.9808,  48.2044, 435.8093,                       & 
     &                 10.0281,  15.7574, 186.2109,                       & 
     &                 15.4127,  55.9559, 243.6311,                       & 
     &                  6.0057,  23.5757,  76.4185                    / 
      DATA GC       /   0.004 ,   1.359 ,   6.454 ,  10.376 ,             & 
     &                  0.008 ,  16.546 ,  21.614 ,                       & 
     &                  5.688 ,  15.801 ,  26.269 ,                       & 
     &                  6.691 ,  25.034 ,  40.975 ,                       & 
     &                 17.604 ,  36.180 ,  47.133 ,                       & 
     &                  8.005 ,  40.804 ,  54.492                     / 
! 
      DATA AN       /  14.0499,  30.8008, 883.1443,                       & 
     &                 10.0000,  16.0000,  64.0000,                       & 
     &                  8.0462,   6.2669,  17.8696, 282.8084,             & 
     &                  7.3751,  33.1390, 215.4829,                       & 
     &                  4.0003,  19.3533,  80.6462,                       & 
     &                 13.0998,  19.6425,  94.3035, 370.9539,             & 
     &                 16.0000,  38.0000,                                 & 
     &                 10.3289,  14.5021, 187.1624, 108.1615, 191.8383,   & 
     &                  6.0044,  23.5612,  76.4344                    / 
      DATA GN       /   2.554 ,   9.169 ,  13.651 ,                       & 
     &                 12.353 ,  13.784 ,  14.874 ,                       & 
     &                  0.014 ,   2.131 ,  15.745 ,  24.949 ,             & 
     &                  6.376 ,  14.246 ,  29.465 ,                       & 
     &                  0.022 ,  31.259 ,  41.428 ,                       & 
     &                  7.212 ,  15.228 ,  34.387 ,  46.708 ,             & 
     &                 46.475 ,  49.468 ,                                 & 
     &                  8.693 ,  37.650 ,  65.479 ,  61.155 ,  79.196 ,   & 
     &                  9.999 ,  60.991 ,  82.262                     / 
! 
      DATA AO       /   4.0029,   5.3656,  36.2853,1044.3447,             & 
     &                131.0217, 868.9779,  14.8533,  93.1466,             & 
     &                 12.7843,   5.6828,  98.0919, 829.4396,             & 
     &                 50.9878, 199.0120,   2.0000,   6.0000,  10.0000,   & 
     &                 10.0000,  30.0000,  50.0000,                       & 
     &                  8.0703,   5.7144,  84.1156, 529.0927,             & 
     &                  5.6609,  28.9355, 111.3620, 494.0413,             & 
     &                 45.5249, 134.4751,                                 & 
     &                  4.0003,  21.2937,  78.7058,                       & 
     &                 12.8293,  16.2730, 123.6578, 327.2396,             & 
     &                 48.7883, 102.2117,  20.0060, 161.9903,             & 
     &                 28.4184,  61.5816,                                 & 
     &                 10.5563,  13.2950, 188.1390,                       & 
     &                 14.6560, 129.4922, 470.8512                    / 
      DATA GO       /   0.022 ,   2.019 ,   9.812 ,  13.087 ,             & 
     &                 13.804 ,  16.061 ,  14.293 ,  16.114 ,             & 
     &                  3.472 ,   7.437 ,  22.579 ,  32.035 ,             & 
     &                 27.774 ,  33.678 ,  28.118 ,  31.019 ,  34.204 ,   & 
     &                 30.892 ,  33.189 ,  36.181 ,                       & 
     &                  0.032 ,   2.760 ,  35.328 ,  48.277 ,             & 
     &                  7.662 ,  16.786 ,  42.657 ,  54.522 ,             & 
     &                 50.204 ,  56.044 ,                                 & 
     &                  0.048 ,  50.089 ,  66.604 ,                       & 
     &                  8.954 ,  18.031 ,  57.755 ,  72.594 ,             & 
     &                 68.388 ,  82.397 ,  31.960 ,  76.876 ,             & 
     &                 75.686 ,  80.388 ,                                 & 
     &                 10.747 ,  52.323 ,  94.976 ,                       & 
     &                 27.405 ,  86.350 , 109.917                     / 
! 
      DATA AF       /   2.0001,  39.9012, 122.0986,                       & 
     &                 10.0000,  30.0000,  50.0000,                       & 
     &                  4.0199,   5.5741,  22.1839, 190.2179,             & 
     &                 53.0383, 126.9616,  31.6894,  75.3105,             & 
     &                 13.5014,   7.9936,  55.7981, 298.7039,             & 
     &                 26.2496,  63.7503,   2.0000,   6.0000,  10.0000,   & 
     &                 28.7150,  71.2850,                                 & 
     &                  8.0153,   6.1931,  21.7287,  48.7780, 278.2782,   & 
     &                178.5560, 421.4435,  51.7632,  95.2368          / 
      DATA GF       /   0.050 ,  13.317 ,  15.692 ,                       & 
     &                 15.361 ,  17.128 ,  18.498 ,                       & 
     &                  0.048 ,   2.735 ,  20.079 ,  30.277 ,             & 
     &                 27.548 ,  32.532 ,  30.391 ,  34.707 ,             & 
     &                  4.479 ,  12.072 ,  31.662 ,  51.432 ,             & 
     &                 44.283 ,  50.964 ,  46.193 ,  50.436 ,  54.880 ,   & 
     &                 50.816 ,  57.479 ,                                 & 
     &                  0.058 ,   3.434 ,  14.892 ,  37.472 ,  69.883 ,   & 
     &                 67.810 ,  83.105 ,  72.435 ,  79.747           / 
! 
      DATA ANN      /  34.5080, 365.4919,  16.5768, 183.4231,             & 
     &                  2.0007,  89.5607, 380.4381,  26.4473,  63.5527,   & 
     &                  4.0342,   5.6162,  11.5176,  72.8273,             & 
     &                 48.5684, 131.4315,  31.1710,  76.8290,             & 
     &                 14.0482,  13.3077,  52.7897, 467.8487,             & 
     &                 54.2196, 195.7800                              / 
      DATA GNN      /  17.796 ,  20.730 ,  17.879 ,  20.855 ,             & 
     &                  0.097 ,  29.878 ,  37.221 ,  31.913 ,  37.551 ,   & 
     &                  0.092 ,   3.424 ,  24.806 ,  46.616 ,             & 
     &                 45.643 ,  54.147 ,  48.359 ,  57.420 ,             & 
     &                  5.453 ,  18.560 ,  46.583 ,  80.101 ,             & 
     &                 70.337 ,  85.789                               / 
! 
      DATA ANA      /  11.6348, 158.3593,                                 & 
     &                 21.0453,  50.9546,  10.1389,  25.8611,             & 
     &                  2.0019,  38.0569, 137.9398,  28.3106,  61.6893,   & 
     &                  4.0334,   5.8560,  18.1786, 208.9142,             & 
     &                 93.6895, 406.3095,  60.4276, 239.5719          / 
      DATA GNA      /   2.400 ,   4.552 ,                                 & 
     &                 34.367 ,  40.566 ,  34.676 ,  40.764 ,             & 
     &                  0.170 ,  44.554 ,  57.142 ,  51.689 ,  60.576 ,   & 
     &                  0.152 ,   4.260 ,  36.635 ,  83.254 ,             & 
     &                 72.561 ,  89.475 ,  75.839 ,  92.582           / 
! 
      DATA AMG      /  10.7445, 291.5057,  53.7488,                       & 
     &                  6.2270,  31.1291, 132.6438,                       & 
     &                 40.4379, 159.5618,  20.3845,  79.6154,             & 
     &                  2.0007, 106.8977, 343.1010,  10.1326, 237.8581/ 
      DATA GMG      /   2.805 ,   6.777 ,   9.254 ,                       & 
     &                  4.459 ,   9.789 ,  13.137 ,                       & 
     &                 57.413 ,  71.252 ,  58.010 ,  71.660 ,             & 
     &                  0.276 ,  74.440 ,  94.447 ,  54.472 ,  95.858 / 
! 
      DATA AAL      /   4.0009,  11.7804, 142.2179,  13.6585,  96.3371,   & 
     &                 10.0807,  49.5843, 285.3343,  14.6872,  59.3122,   & 
     &                  6.3277,  29.5086, 134.1634,                       & 
     &                 46.3164, 153.6833,  22.9896,  77.0103          / 
      DATA GAL      /   0.014 ,   3.841 ,   5.420 ,   3.727 ,   8.833 ,   & 
     &                  4.749 ,  11.902 ,  16.719 ,  11.310 ,  18.268 ,   & 
     &                  6.751 ,  16.681 ,  24.151 ,                       & 
     &                 83.551 , 104.787 ,  84.293 , 105.171           / 
! 
      DATA ASI      /   7.9658,   4.6762,   1.3512, 123.2267, 443.7797,   & 
     &                  4.0000,   7.4186,  24.1754, 60.4060,              & 
     &                 14.4695,  11.9721,  26.5062, 269.0521,             & 
     &                  9.1793,   4.8766,  29.1442,  52.7998,             & 
     &                 13.2674,  36.0417, 180.6910,                       & 
     &                  6.4839,  27.6851, 135.8301                    / 
      DATA GSI      /   0.020 ,   0.752 ,   1.614 ,   5.831 ,   7.431 ,   & 
     &                  0.036 ,   8.795 ,  11.208 ,  13.835 ,             & 
     &                  5.418 ,   7.825 ,  14.440 ,  19.412 ,             & 
     &                  6.572 ,  11.449 ,  18.424 ,  25.457 ,             & 
     &                 15.682 ,  27.010 ,  34.599 ,                       & 
     &                  9.042 ,  24.101 ,  37.445                     / 
! 
      DATA AP       /  13.5211,  22.2130, 353.2583,  10.0000, 150.0000,   & 
     &                  8.0241,   5.8085,  51.7542, 252.4002,             & 
     &                  4.0021,  20.7985,  62.4194, 200.7786,             & 
     &                 11.7414,  63.5124, 179.7420,                       & 
     &                  6.8835,  32.7777, 228.3366                    / 
      DATA GP       /   1.514 ,   5.575 ,   9.247 ,   8.076 ,  10.735 ,   & 
     &                  0.043 ,   1.212 ,   8.545 ,  15.525 ,             & 
     &                  0.074 ,   7.674 ,  16.639 ,  25.118 ,             & 
     &                  8.992 ,  24.473 ,  40.704 ,                       & 
     &                 11.464 ,  33.732 ,  55.455                     / 
! 
      DATA AS       /   3.9615,   5.0780,  15.0944, 362.8588,             & 
     &                 51.5995, 268.4002,  12.0000, 276.0000,             & 
     &                 11.4377,   5.5126, 141.0009, 254.0478,             & 
     &                 33.0518, 126.9479,                                 & 
     &                  4.0707,   4.0637,   5.7245, 144.6376, 106.4909,   & 
     &                  4.0011,  19.2813,  27.5990,  35.1179,             & 
     &                 94.7454, 283.2486,                                 & 
     &                 10.5474,  28.7137,  65.7378,  24.0000          / 
      DATA GS       /   0.053 ,   1.121 ,   5.812 ,   9.425 ,             & 
     &                  8.936 ,  11.277 ,   9.600 ,  12.551 ,             & 
     &                  1.892 ,   3.646 ,  13.550 ,  19.376 ,             & 
     &                 16.253 ,  21.062 ,                                 & 
     &                  0.043 ,   0.123 ,   1.590 ,  13.712 ,  22.050 ,   & 
     &                  0.118 ,   9.545 ,  18.179 ,  31.441 ,             & 
     &                 30.664 ,  56.150 ,                                 & 
     &                 10.704 ,  27.075 ,  50.599 ,  43.034           / 
! 
      DATA ACL      /   2.0007,  62.5048, 669.4942,  29.0259, 130.9740,   & 
     &                  3.9064,   0.3993,   5.3570,  60.3424, 119.9913,   & 
     &                138.1567, 278.8418, 102.3681, 158.6314,             & 
     &                 12.6089,   5.9527, 110.5635, 262.8715,             & 
     &                 69.2035, 100.7960,                                 & 
     &                  7.3458,   5.6638,  44.1256, 202.7846,             & 
     &                  4.0037,  21.8663,  40.5363,  57.5919          / 
      DATA GCL      /   0.110 ,   9.919 ,  12.280 ,  11.017 ,  13.532 ,   & 
     &                  0.092 ,   0.581 ,   1.620 ,  13.121 ,  19.787 ,   & 
     &                 16.365 ,  21.988 ,  18.065 ,  23.594 ,             & 
     &                  2.358 ,   5.708 ,  19.084 ,  30.683 ,             & 
     &                 24.880 ,  33.229 ,                                 & 
     &                  0.102 ,   1.391 ,  14.709 ,  36.968 ,             & 
     &                  0.185 ,  11.783 ,  25.653 ,  44.698           / 
! 
      DATA AAR      /  43.6623, 324.3375,  20.8298, 163.1701,             & 
     &                  2.0026, 137.4515, 258.5445,  62.8129, 149.1867,   & 
     &                  4.0495,  14.4466,  46.8234, 124.6651,             & 
     &                151.9828, 268.0157, 101.1302, 150.8691,             & 
     &                 13.3718,   8.6528,  60.4614, 285.5072,             & 
     &                  6.7655,   4.7684,  12.8631,  54.5260          / 
      DATA GAR      /  12.638 ,  14.958 ,  12.833 ,  15.139 ,             & 
     &                  0.178 ,  17.522 ,  23.584 ,  20.464 ,  25.150 ,   & 
     &                  0.151 ,   1.561 ,  17.399 ,  30.871 ,             & 
     &                 24.684 ,  33.978 ,  27.091 ,  36.481 ,             & 
     &                  2.810 ,   8.877 ,  24.351 ,  44.489 ,             & 
     &                  0.144 ,   1.160 ,  10.210 ,  27.178           / 
! 
      DATA AK       /  12.9782, 148.6673,   6.3493,                       & 
     &                 66.3444, 101.6553,   4.0001,  13.4465,  46.5534,   & 
     &                  2.0171, 116.4767, 713.4965,  63.5907, 396.4079,   & 
     &                  2.0000,  10.0000,  30.0000,                       & 
     &                  4.0702,   5.7791,  52.6795, 327.4539,             & 
     &                 62.8604, 357.1331,  55.9337, 196.0646,             & 
     &                 10.9275,   5.5398,  43.2761,  76.2560,             & 
     &                 42.0000,  18.0000                              / 
      DATA GK       /   1.871 ,   3.713 ,  18.172 ,                       & 
     &                 21.185 ,  27.705 ,   2.059 ,  23.709 ,  28.542 ,   & 
     &                  0.273 ,  26.709 ,  39.640 ,  31.220 ,  41.865 ,   & 
     &                 29.955 ,  37.557 ,  42.862 ,                       & 
     &                  0.228 ,   2.274 ,  21.703 ,  50.191 ,             & 
     &                 32.145 ,  49.262 ,  34.155 ,  51.718 ,             & 
     &                  3.043 ,   5.479 ,  20.547 ,  30.680 ,             & 
     &                 36.275 ,  47.345                               / 
! 
      DATA ACA      /  18.2366,  27.5012, 149.2617,  94.5242, 705.4711,   & 
     &                 11.8706,  14.0710, 106.0547,                       & 
     &                 57.2414, 110.7567,  29.8121,  54.1874,             & 
     &                  2.0184,  97.5784, 282.3939, 209.1871, 252.8129/ 
      DATA GCA      /   2.050 ,   3.349 ,   5.321 ,   4.873 ,   7.017 ,   & 
     &                  1.769 ,   5.109 ,   9.524 ,                       & 
     &                 27.271 ,  41.561 ,  29.172 ,  42.140 ,             & 
     &                  0.394 ,  28.930 ,  52.618 ,  38.593 ,  49.646 / 
! 
      DATA ASC      /   6.0014,  83.1958,  67.3666, 329.4354,             & 
     &                 44.0793, 169.9969, 533.9195,                       & 
     &                 34.1642, 124.8475, 228.9879,                       & 
     &                 11.9979,  16.9280,  28.4778,  82.0418, 234.5360,   & 
     &                  6.0042,   2.7101,  13.9801,  65.3039,             & 
     &                 12.0000,  12.0000,                                 & 
     &                  2.0051,   2.9621,  29.0306                    / 
      DATA GSC      /   0.021 ,   2.056 ,   3.551 ,   5.465 ,             & 
     &                  1.535 ,   3.797 ,   6.203 ,                       & 
     &                  2.389 ,   4.858 ,   7.141 ,                       & 
     &                  0.011 ,   0.430 ,   1.156 ,   3.711 ,   8.863 ,   & 
     &                  0.025 ,   3.499 ,  10.463 ,  18.606 ,             & 
     &                 41.779 ,  57.217 ,                                 & 
     &                  0.539 ,  24.442 ,  51.079                     / 
! 
      DATA ATI      /   7.0887,   8.9186,  17.5633, 206.6832, 438.5735,   & 
     &                654.1721,                                           & 
     &                 38.0462,  69.6271, 364.2845, 832.0408,             & 
     &                 98.8562,  57.9934, 442.1498,                       & 
     &                 19.7843,  32.0637,  37.0895, 110.6682, 288.4946,   & 
     &                521.8837,                                           & 
     &                 10.0000,  34.0000, 120.0000,                       & 
     &                 16.1691,  22.3550,  24.1646,  83.5128, 222.7963,   & 
     &                  6.0020,   4.6177,  25.2636,  52.1162,             & 
     &                 12.0000,   8.0000                              / 
      DATA GTI      /   0.021 ,   0.048 ,   1.029 ,   2.183 ,   4.109 ,   & 
     &                  5.785 ,                                           & 
     &                  0.846 ,   1.792 ,   3.836 ,   5.787 ,             & 
     &                  2.561 ,   4.869 ,   6.340 ,                       & 
     &                  0.023 ,   0.124 ,   0.774 ,   1.810 ,   4.980 ,   & 
     &                  9.585 ,                                           & 
     &                  1.082 ,   4.928 ,  11.279 ,                       & 
     &                  0.041 ,   1.375 ,   4.768 ,  10.985 ,  19.769 ,   & 
     &                  0.048 ,  11.577 ,  24.531 ,  36.489 ,             & 
     &                 54.436 ,  75.373                               / 
! 
      DATA AV       /  15.2627,  23.9869,  51.3053, 570.3384,1650.9417,   & 
     &                162.2829, 298.8303, 908.8852,                       & 
     &                 23.6736,  37.1624,  86.8011, 300.7440, 864.5880,   & 
     &                 57.8961,  79.4605, 214.9007, 864.7425,             & 
     &                 61.8508,  64.0845, 192.8298, 718.2349,             & 
     &                 23.8116,  68.2495, 135.0613, 536.7632,             & 
     &                 15.9543,  22.5542,  71.4921, 248.9544,             & 
     &                  6.0006,   5.8785,  50.5077,  97.6129          / 
      DATA GV       /   0.026 ,   0.145 ,   0.718 ,   2.586 ,   5.458 ,   & 
     &                  2.171 ,   4.153 ,   6.097 ,                       & 
     &                  0.009 ,   0.366 ,   1.504 ,   5.294 ,  10.126 ,   & 
     &                  1.796 ,   2.353 ,   6.068 ,  12.269 ,             & 
     &                  2.560 ,   3.674 ,   6.593 ,  12.880 ,             & 
     &                  0.045 ,   1.684 ,   8.162 ,  21.262 ,             & 
     &                  0.065 ,   1.746 ,  15.158 ,  33.141 ,             & 
     &                  0.077 ,  21.229 ,  44.134 ,  60.203           / 
! 
      DATA ACR      /  30.1842,  79.2847, 149.5293,                       & 
     &                215.3696, 119.1974, 741.4321,                       & 
     &                184.9946,1352.5038, 784.4937,                       & 
     &                 46.6191, 160.1361, 488.0449, 657.1928,             & 
     &                 47.1742, 267.0275, 441.1324, 150.6650,             & 
     &                 24.3768, 122.8359, 285.5092, 794.1654,             & 
     &                 24.2296,  75.0258, 172.9452, 543.6511,             & 
     &                 15.9819,  17.6800,  95.2003, 225.0947          / 
      DATA GCR      /   0.993 ,   3.070 ,   5.673 ,                       & 
     &                  3.339 ,   4.801 ,   7.198 ,                       & 
     &                  2.829 ,   4.990 ,   7.643 ,                       & 
     &                  1.645 ,   3.727 ,   7.181 ,  12.299 ,             & 
     &                  2.902 ,   4.273 ,   8.569 ,  14.912 ,             & 
     &                  0.047 ,   2.566 ,   9.441 ,  21.198 ,             & 
     &                  0.078 ,   2.242 ,  15.638 ,  32.725 ,             & 
     &                  0.103 ,   2.146 ,  26.153 ,  49.381           / 
! 
      DATA AMN      /  53.9107,  81.3931, 546.6945 ,                      & 
     &                144.1893, 407.8029,  45.6177, 298.4423,2410.9335,   & 
     &                 22.6382,  93.8419, 183.9367, 907.5765,             & 
     &                137.0409, 168.6783, 329.0287, 773.2513,             & 
     &                 70.1925,  72.3372, 213.9512, 539.5165,             & 
     &                 24.2373,  93.5415, 456.6167, 506.5484,             & 
     &                 24.7687,  66.9896, 264.1853, 484.0161          / 
      DATA GMN      /   2.527 ,   4.204 ,   6.602 ,                       & 
     &                  4.155 ,   7.321 ,   2.285 ,   5.631 ,   8.448 ,   & 
     &                  1.496 ,   3.839 ,   7.751 ,  13.484 ,             & 
     &                  3.681 ,   6.054 ,   9.934 ,  14.936 ,             & 
     &                  3.531 ,   6.967 ,  15.222 ,  25.069 ,             & 
     &                  0.071 ,   2.896 ,  20.725 ,  37.383 ,             & 
     &                  0.126 ,   2.660 ,  28.528 ,  53.413           / 
! 
      DATA AFE      /  14.4102,   2.7050, 421.6612, 940.1484,             & 
     &                 36.2187,  22.8883, 239.5997, 825.2919,             & 
     &                110.0242, 992.3040, 640.6715,                       & 
     &                 17.0494,  32.3783,  34.3184, 420.9626,1067.2064,   & 
     &                154.0059, 462.1117, 329.8618,                       & 
     &                 15.7906,  47.1186, 279.9292, 692.1005,             & 
     &                 91.0206, 206.3082, 706.9927, 836.6689,             & 
     &                 40.0790,  27.6965,  28.2243,  18.0001,             & 
     &                 24.0899,  89.6340,  51.5756, 241.6980          / 
      DATA GFE      /   0.066 ,   0.339 ,   2.897 ,   6.585 ,             & 
     &                  0.923 ,   1.679 ,   4.620 ,   7.053 ,             & 
     &                  4.249 ,   5.875 ,   7.781 ,                       & 
     &                  0.062 ,   0.283 ,   1.504 ,   5.430 ,  11.210 ,   & 
     &                  2.792 ,   7.627 ,  13.623 ,                       & 
     &                  0.077 ,   3.723 ,  12.137 ,  23.700 ,             & 
     &                  2.688 ,   7.595 ,  15.444 ,  25.587 ,             & 
     &                  3.982 ,   4.677 ,   6.453 ,  23.561 ,             & 
     &                  0.102 ,   3.354 ,  22.954 ,  33.796           / 
! 
      DATA ACO      /  11.9120,  20.4424,  28.3863, 132.5038, 600.7461,   & 
     &                 33.3092, 237.4331, 977.2502,                       & 
     &                 55.5396, 318.8169, 619.6366,                       & 
     &                 32.6900,  83.8694, 107.4378,                       & 
     &                 11.2593,  38.2239,  22.9964, 261.3486, 637.1485,   & 
     &                 23.0233,  41.6599, 264.6460, 181.6699,             & 
     &                 16.0356,   7.8633,  70.3158, 423.3512, 742.3553,   & 
     &                  0.                                            / 
      DATA GCO      /   0.112 ,   0.341 ,   0.809 ,   3.808 ,   6.723 ,   & 
     &                  2.057 ,   3.484 ,   7.210 ,                       & 
     &                  2.405 ,   5.133 ,   8.097 ,                       & 
     &                  2.084 ,   5.291 ,   8.426 ,                       & 
     &                  0.135 ,   0.517 ,   1.606 ,   6.772 ,  12.622 ,   & 
     &                  2.512 ,   4.348 ,   8.253 ,  15.377 ,             & 
     &                  0.132 ,   0.863 ,   3.086 ,  11.789 ,  23.263 ,   & 
     &                  0.                                            / 
! 
      DATA ANI      /   7.1268,  12.4486,  11.9953,  10.0546, 114.1658,   & 
     &                391.2064,                                           & 
     &                 26.3908, 213.8081, 938.7927,                       & 
     &                  4.1421,  37.3781,  25.9712, 333.3397, 311.1633,   & 
     &                 33.1031, 184.1854, 136.7072,                       & 
     &                 11.1915,   5.4174,  53.6793, 460.6781, 380.0056,   & 
     &                  0.                                            / 
      DATA GNI      /   0.026 ,   0.137 ,   0.315 ,   1.778 ,   4.029 ,   & 
     &                  6.621 ,                                           & 
     &                  2.249 ,   4.042 ,   7.621 ,                       & 
     &                  0.191 ,   1.235 ,   3.358 ,   8.429 ,  17.096 ,   & 
     &                  3.472 ,   9.065 ,  16.556 ,                       & 
     &                  0.194 ,   1.305 ,   5.813 ,  14.172 ,  26.169 ,   & 
     &                  0.                                            / 
! 
      DATA ACU      /  11.0549, 238.9423,  10.3077, 126.2990,1073.3876,   & 
     &                 30.0000,  50.0000,  60.0000,                       & 
     &                 19.2984,  50.5974, 240.2021,1216.9016,             & 
     &                 48.3048, 583.2011, 320.4931,                       & 
     &                  4.0155,  70.3264, 313.1213, 536.5331,             & 
     &                  0.                                            / 
      DATA GCU      /   4.212 ,   7.227 ,   1.493 ,   5.859 ,   9.709 ,   & 
     &                  7.081 ,   9.362 ,  10.130 ,                       & 
     &                  2.865 ,   8.260 ,  14.431 ,  18.292 ,             & 
     &                  9.650 ,  14.640 ,  24.320 ,                       & 
     &                  0.337 ,   8.520 ,  16.925 ,  28.342 ,             & 
     &                  0.                                            / 
! 
      DATA AZN      /  15.9880, 484.0042,  18.5863, 123.4134,             & 
     &                  3.0000, 189.0000,                                 & 
     &                  6.1902,  38.9317, 204.8780,                       & 
     &                 10.2588,  89.3771, 370.3640,  30.0000, 128.0000,   & 
     &                 24.6904, 106.7491, 439.5586,                       & 
     &                  0.                                            / 
      DATA GZN      /   4.546 ,   8.840 ,  10.247 ,  16.620 ,             & 
     &                 11.175 ,  16.321 ,                                 & 
     &                  6.113 ,  12.964 ,  16.444 ,                       & 
     &                  7.926 ,  13.633 ,  24.353 ,  16.286 ,  24.910 ,   & 
     &                 10.291 ,  20.689 ,  32.077 ,                       & 
     &                  0.                                            / 
! 
      DATA ICOMP /0/ 
! 
!     save indexs,indexm,index0,is,im,ig0,igpr,                           & 
!    &     xl, chion, alf, gam 
! 
      IF(ICOMP.EQ.0) THEN 
         IND=1 
         DO K=1,NIONS 
            INDEXS(K)=IND 
            IND=IND+IS(K) 
         END DO 
         IND=1 
         DO K=1,NSS 
            INDEXM(K)=IND 
            IND=IND+IM(K) 
         END DO 
         ICOMP=1 
      END IF 
 
!C    IF(PFUNC(IZI,IAT,ID).GT.0.) THEN 
!C       U=PFUNC(IZI,IAT,ID) 
!C       RETURN 
!C    END IF 
! 
      IF((IAT.EQ.26.or.iat.eq.28)                                         & 
     &  .AND.IZI.GE.4.AND.IZI.LE.9) THEN 
         if(iat.eq.26) call pffe(IZI,T,ANE,U) 
         if(iat.eq.28) call pfni(izi,t,u,dut,dun) 
         RETURN 
      END IF 
! 
      IF(IAT.GT.30.AND.IZI.LE.3) THEN 
         CALL PFHEAV(IAT,IZI,3,T,ANE,U) 
         RETURN 
      END IF 
      IF(IAT.GT.8 .AND. IZI.GT.5) then 
         u=igle(iat-izi+1) 
         return 
      end if 
! 
!     Irwin partition functions by default 
! 
      if(iirwin.gt.0.and.t.lt.16000.) then 
         if(izi.le.2) then 
            call irwpf(iat,izi,0,t,u0) 
            u=u0 
            return 
         end if 
       else if(iat.gt.30.and.izi.le.3) then 
         CALL PFHEAV(IAT,IZI,3,T,ANE,U) 
      end if 
! 
      IF(IZI.LE.0.OR.IZI.GT.9.OR.IAT.LE.0.OR.IAT.GT.30) THEN 
         CALL PFSPEC(IAT,IZI,T,ANE,U) 
         RETURN 
      END IF 
! 
      MODE=MODPF(IAT) 
      IF(MODE.LT.0) THEN 
         CALL PFSPEC(IAT,IZI,T,ANE,U) 
       ELSE IF(MODE.GT.0) THEN 
         U=IGLE(IAT-IZI+1) 
       ELSE 
         I0=INDEX0(IZI,IAT) 
         IF(I0.GT.0) THEN 
            QZ=IZI 
            XMAX=XMAXN 
            THET=5040.4/T 
            A=31.321*QZ*QZ*THET 
            XMAX2=XMAX*XMAX 
            QAS1=XMAX*THIRD*(XMAX2+TRHA*XMAX+HALF) 
            IS0=INDEXS(I0) 
            ISS=IS0+IS(I0)-1 
            SU1=0. 
            SQA=0. 
            DO K=IS0,ISS 
               XXL=XL(K) 
               GPR=IGPR(K) 
               X=CHION(K)*THET 
               EX=0. 
               IF(X.LT.30) EX=EXP(-X*2.30258029299405) 
               QAS=(QAS1-XXL*THIRD*(XXL*XXL+TRHA*XXL+HALF)+(XMAX-XXL)*    & 
     &             (UN+A*HALF/XXL/XMAX)*A)*GPR*EX 
               SQA=SQA+QAS 
               M0=INDEXM(K) 
               M1=M0+IM(K)-1 
               AL1=0. 
               DO M=M0,M1 
                  XG=GAM(M)*THET 
                  IF(XG.LE.20.) THEN 
                     XM=EXP(-XG*2.30258029299405)*ALF(M) 
                     AL1=AL1+XM 
                  END IF 
               END DO 
               SU1=SU1+AL1 
            END DO 
            U=IG0(I0) 
            U=U+SU1+SQA 
            IF(U.LT.0.) U=IG0(I0) 
          ELSE IF(I0.LT.0) THEN 
            U=FLOAT(-I0) 
          ELSE 
            CALL PFSPEC(IAT,IZI,T,ANE,U) 
         END IF 
      END IF 
      RETURN 
      END SUBROUTINE PARTF 
! 
! ******************************************************************** 
! 
! 
      subroutine pffe(ion,t,ane,pf) 
!     ============================= 
! 
!     partition functions for Fe IV to Fe IX 
!     after Fischel and Sparks, 1971, NASA SP-3066 
! 
!     Output:  PF   partition function 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      real(dp) :: tt(50),pn(10) 
      integer  :: nca(9) 
      real(dp) :: p4a(22),p4b(10,28),                                     & 
     &            p5a(30),p5b(10,20),                                     & 
     &            p6a(37),p6b(10,13),                                     & 
     &            p7a(40),p7b(10,10),                                     & 
     &            p8a(41),p8b(10,9),                                      & 
     &            p9a(45),p9b(10,5) 
! 
      real(dp), parameter :: xen=2.302585093,xmil=0.001,xmilen=xmil*xen   & 
      real(dp), parameter :: xbtz=1.38054e-16 
! 
      data nca /3*0,22,30,37,40,41,45/                                    & 
     &     nne /10/ 
! 
      data tt /                                                           & 
     & 3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,      & 
     & 20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,                       & 
     & 32.,34.,36.,38.,40.,42.,44.,46.,48.,                               & 
     & 50.,55.,60.,65.,70.,75.,80.,85.,90.,95.,100.,125.,150./ 
! 
      data pn /-2.,-1.,0.,1.,2.,3.,4.,5.,6.,7./ 
! 
      data p4a /                                                          & 
     & 0.778, 0.778, 0.778, 0.779, 0.783, 0.789, 0.801, 0.818,            & 
     & 0.842, 0.871, 0.906, 0.945, 0.987, 1.030, 1.074, 1.117,            & 
     & 1.160, 1.201, 1.242, 1.280, 1.317, 1.353/ 
! 
      data p4b /                                                          & 
     & 1.406,1.393,1.389,7*1.387,                                         & 
     & 1.464,1.434,1.424,1.421,1.420,5*1.419,                             & 
     & 1.546,1.483,1.461,1.454,1.451,1.451,4*1.450,                       & 
     & 1.665,1.547,1.503,1.488,1.482,1.481,4*1.480,                       & 
     & 1.826,1.636,1.553,1.524,1.514,1.510,4*1.509,                       & 
     & 2.024,1.755,1.618,1.564,1.546,1.540,1.538,3*1.537,                 & 
     & 2.480,2.087,1.814,1.674,1.619,1.599,1.593,1.591,1.590,1.590,       & 
     & 2.945,2.489,2.105,1.846,1.717,1.667,1.649,1.643,1.641,1.640,       & 
     & 3.379,2.897,2.452,2.089,1.859,1.751,1.710,1.696,1.691,1.689,       & 
     & 3.774,3.283,2.808,2.381,2.054,1.864,1.782,1.751,1.741,1.738,       & 
     & 4.133,3.637,3.150,2.688,2.292,2.015,1.871,1.814,1.793,1.786,       & 
     & 4.460,3.962,3.468,2.989,2.549,2.199,1.984,1.886,1.848,1.835,       & 
     & 4.757,4.258,3.762,3.274,2.809,2.406,2.121,1.972,1.908,1.886,       & 
     & 5.029,4.530,4.032,3.539,3.061,2.624,2.279,2.073,1.976,1.939,       & 
     & 5.279,4.780,4.281,3.785,3.299,2.840,2.450,2.189,2.051,1.996,       & 
     & 5.510,5.010,4.511,4.013,3.522,3.050,2.628,2.318,2.136,2.057,       & 
     & 6.014,5.514,5.014,4.515,4.018,3.530,3.065,2.666,2.381,2.228,       & 
     & 6.435,5.935,5.435,4.936,4.437,3.943,3.460,3.022,2.658,2.422,       & 
     & 6.794,6.294,5.794,5.294,4.794,4.297,3.807,3.343,2.939,2.631,       & 
     & 7.102,6.602,6.102,5.602,5.102,4.604,4.110,3.638,3.194,2.845,       & 
     & 7.370,6.870,6.370,5.870,5.370,4.871,4.375,3.892,3.439,3.052,       & 
     & 7.606,7.106,6.606,6.106,5.605,5.106,4.608,4.125,3.661,3.249,       & 
     & 7.815,7.315,6.814,6.314,5.814,5.314,4.816,4.333,3.851,3.418,       & 
     & 8.001,7.501,7.001,6.500,6.000,5.500,5.001,4.511,4.032,3.586,       & 
     & 8.168,7.668,7.168,6.668,6.168,5.667,5.168,4.680,4.197,3.741,       & 
     & 8.319,7.819,7.319,6.819,6.319,5.818,5.319,4.832,4.347,3.884,       & 
     & 8.900,8.399,7.899,7.399,6.899,6.398,5.898,5.405,4.917,4.431,       & 
     & 9.294,8.794,8.294,7.793,7.293,6.793,6.292,5.799,5.306,4.824/ 
! 
      data p5a /                                                          & 
     & 1.235, 1.276, 1.301, 1.321, 1.339, 1.359, 1.381, 1.405,            & 
     & 1.432, 1.460, 1.489, 1.518, 1.546, 1.574, 1.601, 1.627,            & 
     & 1.652, 1.675, 1.697, 1.718, 1.738, 1.757, 1.775, 1.792,            & 
     & 1.808, 1.823, 1.838, 1.851, 1.877, 1.900/ 
! 
      data p5b /                                                          & 
     & 1.943,1.928,1.923,7*1.921,                                         & 
     & 2.011,1.964,1.947,1.942,1.941,5*1.940,                             & 
     & 2.144,2.025,1.980,1.965,1.960,1.958,4*1.957,                       & 
     & 2.361,2.137,2.032,1.993,1.980,1.976,1.975,3*1.974,                 & 
     & 2.646,2.315,2.121,2.035,2.004,1.994,1.991,1.990,1.989,1.989,       & 
     & 2.960,2.553,2.260,2.102,2.037,2.015,2.007,2.005,2.004,2.004,       & 
     & 3.274,2.823,2.450,2.205,2.086,2.040,2.025,2.020,2.018,2.018,       & 
     & 3.575,3.101,2.674,2.348,2.158,2.075,2.045,2.036,2.032,2.031,       & 
     & 4.251,3.757,3.275,2.829,2.466,2.234,2.124,2.083,2.069,2.064,       & 
     & 4.822,4.324,3.829,3.346,2.895,2.522,2.278,2.161,2.116,2.100,       & 
     & 5.308,4.808,4.310,3.816,3.334,2.888,2.525,2.297,2.187,2.145,       & 
     & 5.725,5.225,4.726,4.228,3.736,3.260,2.828,2.496,2.294,2.206,       & 
     & 6.088,5.589,5.089,4.590,4.093,3.604,3.139,2.733,2.447,2.291,       & 
     & 6.407,5.907,5.407,4.908,4.409,3.915,3.433,2.988,2.629,2.399,       & 
     & 6.689,6.189,5.689,5.189,4.690,4.193,3.704,3.236,2.832,2.535,       & 
     & 6.940,6.440,5.940,5.440,4.941,4.443,3.949,3.469,3.038,2.687,       & 
     & 7.166,6.666,6.166,5.666,5.166,4.667,4.171,3.684,3.237,2.847,       & 
     & 7.370,6.870,6.369,5.869,5.369,4.870,4.373,3.882,3.417,3.008,       & 
     & 8.150,7.649,7.149,6.649,6.149,5.649,5.149,4.651,4.167,3.700,       & 
     & 8.677,8.177,7.676,7.176,6.676,6.176,5.676,5.176,4.687,4.203/ 
! 
      data p6a /                                                          & 
     & 1.218, 1.273, 1.309, 1.335, 1.358, 1.379, 1.400, 1.421,            & 
     & 1.442, 1.463, 1.484, 1.504, 1.523, 1.542, 1.560, 1.577,            & 
     & 1.594, 1.609, 1.624, 1.638, 1.652, 1.664, 1.677, 1.688,            & 
     & 1.699, 1.709, 1.719, 1.729, 1.746, 1.762, 1.777, 1.790,            & 
     & 1.803, 1.814, 1.825, 1.834, 1.843/ 
! 
      data p6b /                                                          & 
     & 1.862,1.855,1.853,7*1.852,                                         & 
     & 1.958,1.900,1.880,1.874,1.872,5*1.871,                             & 
     & 2.264,2.045,1.944,1.906,1.894,1.890,4*1.888,                       & 
     & 2.776,2.386,2.119,1.984,1.930,1.912,1.906,1.904,2*1.903,           & 
     & 3.321,2.856,2.453,2.165,2.012,1.949,1.927,1.920,1.918,1.917,       & 
     & 3.821,3.333,2.868,2.465,2.178,2.025,1.963,1.941,1.934,1.932,       & 
     & 4.266,3.771,3.285,2.825,2.434,2.164,2.027,1.972,1.953,1.947,       & 
     & 4.662,4.164,3.670,3.187,2.739,2.372,2.135,2.022,1.980,1.965,       & 
     & 5.015,4.516,4.019,3.527,3.052,2.624,2.295,2.102,2.019,1.988,       & 
     & 5.332,4.832,4.344,3.838,3.351,2.889,2.493,2.217,2.075,2.017,       & 
     & 5.618,5.118,4.619,4.121,3.628,3.149,2.711,2.364,2.155,2.058,       & 
     & 6.710,6.210,5.710,5.210,4.711,4.213,3.719,3.241,2.807,2.462,       & 
     & 7.446,6.946,6.446,5.946,5.446,4.946,4.447,3.952,3.474,3.022/ 
! 
      data p7a /                                                          & 
     & 1.074,1.130,1.167,1.194,1.215,1.234,1.250,1.266,1.280,1.293,       & 
     & 1.306,1.318,1.329,1.340,1.350,1.360,1.369,1.378,1.386,1.394,       & 
     & 1.401,1.408,1.415,1.421,1.427,1.433,1.439,1.444,1.454,1.463,       & 
     & 1.471,1.479,1.486,1.492,1.498,1.504,1.509,1.514,1.525,1.534/ 
! 
      data p7b /                                                          & 
     & 1.555,1.546,1.544,1.543,6*1.542,                                   & 
     & 1.617,1.572,1.557,1.552,1.550,1.550,4*1.549,                       & 
     & 1.798,1.648,1.587,1.566,1.559,1.557,4*1.556,                       & 
     & 2.134,1.832,1.666,1.597,1.573,1.565,1.563,1.562,2*1.561,           & 
     & 2.550,2.138,1.836,1.671,1.602,1.578,1.570,1.568,2*1.567,           & 
     & 2.968,2.504,2.102,1.816,1.665,1.603,1.582,1.575,2*1.572,           & 
     & 3.359,2.875,2.419,2.037,1.779,1.651,1.601,1.584,1.579,1.577,       & 
     & 3.718,3.224,2.745,2.305,1.953,1.736,1.636,1.599,1.586,1.582,       & 
     & 5.097,4.598,4.098,3.601,3.110,2.638,2.217,1.899,1.719,1.643,       & 
     & 6.026,5.526,5.026,4.527,4.028,3.531,3.042,2.576,2.170,1.885/ 
! 
      data p8a /                                                          & 
     & 0.809,0.849,0.875,0.894,0.908,0.918,0.927,0.934,0.939,0.944,       & 
     & 0.948,0.952,0.955,0.958,0.960,0.962,0.964,0.966,0.967,0.969,       & 
     & 0.970,0.971,0.973,0.974,0.975,0.975,0.976,0.977,0.978,0.980,       & 
     & 0.981,0.982,0.983,0.984,0.984,0.985,0.986,0.986,0.987,0.988,       & 
     & 0.989/ 
! 
      data p8b /                                                          & 
     & 0.992,0.991,8*0.990,                                               & 
     & 1.000,0.994,0.992,7*0.991,                                         & 
     & 1.032,1.005,0.996,0.993,0.992,5*0.991,                             & 
     & 1.129,1.040,1.008,0.997,0.993,5*0.992,                             & 
     & 1.335,1.132,1.042,1.009,0.998,0.994,0.993,0.993,2*0.992,           & 
     & 1.640,1.312,1.121,1.038,1.007,0.998,0.994,3*0.993,                 & 
     & 1.987,1.573,1.269,1.101,1.030,1.005,0.997,2*0.994,0.993,           & 
     & 3.514,3.017,2.526,2.053,1.628,1.305,1.119,1.039,1.010,1.000,       & 
     & 4.569,4.069,3.569,3.072,2.580,2.103,1.671,1.336,1.136,1.048/ 
! 
      data p9a /39*0.000,0.001,0.002,0.005,0.008,0.014,0.021/ 
! 
      data p9b /                                                          & 
     & 2*0.032,8*0.031,                                                   & 
     & 0.048,0.045,8*0.044,                                               & 
     & 0.076,0.065,0.061,0.060,6*0.059,                                   & 
     & 1.128,0.722,0.429,0.271,0.207,0.184,0.177,0.174,2*0.173,           & 
     & 2.696,2.200,1.712,1.249,0.848,0.564,0.415,0.354,0.333,0.327/ 
! 
      na=nca(ion) 
      nb=50-na 
      pne=log10(ane*xbtz*t) 
      t0=xmil*t 
      j=1 
      if(pne.ge.pn(1)) then 
         if(pne.gt.pn(nne)) then 
           j1=nne 
           j2=nne 
          else 
           j1=nne 
           j2=nne 
           do j=1,nne-1 
              if(pne.ge.pn(j).and.pne.lt.pn(j+1)) exit 
            end do 
          end if 
       else 
         j1=j 
         j2=j1+1 
         if(pne.lt.pn(1)) j2=1 
      end if 
 
      do i=1,49 
         if(t0.ge.tt(i).and.t0.lt.tt(i+1)) exit 
      end do 
 
      i1=i 
      i2=i+1 
      if(t0.gt.tt(50)) then 
        i1=50 
        i2=50 
      endif 
      if(i2.le.na) then 
         select case(ion) 
         case(4) 
           px1=p4a(i1) 
           px2=p4a(i1) 
           py1=p4a(i2) 
           py2=p4a(i2) 
         case(5) 
           px1=p5a(i1) 
           px2=p5a(i1) 
           py1=p5a(i2) 
           py2=p5a(i2) 
         case(6) 
           px1=p6a(i1) 
           px2=p6a(i1) 
           py1=p6a(i2) 
           py2=p6a(i2) 
         case(7) 
           px1=p7a(i1) 
           px2=p7a(i1) 
           py1=p7a(i2) 
           py2=p7a(i2) 
         case(8) 
           px1=p8a(i1) 
           px2=p8a(i1) 
           py1=p8a(i2) 
           py2=p8a(i2) 
         case(9) 
           px1=p9a(i1) 
           px2=p9a(i1) 
           py1=p9a(i2) 
           py2=p9a(i2) 
         end select 
       else if(i1.eq.na) then 
         select case(ion) 
         case(4) 
           px1=p4a(i1) 
           px2=p4a(i1) 
           py1=p4b(j1,i2-na) 
           py2=p4b(j2,i2-na) 
         case(5) 
           px1=p5a(i1) 
           px2=p5a(i1) 
           py1=p5b(j1,i2-na) 
           py2=p5b(j2,i2-na) 
         case(6) 
           px1=p6a(i1) 
           px2=p6a(i1) 
           py1=p6b(j1,i2-na) 
           py2=p6b(j2,i2-na) 
         case(7) 
           px1=p7a(i1) 
           px2=p7a(i1) 
           py1=p7b(j1,i2-na) 
           py2=p7b(j2,i2-na) 
         case(8) 
           px1=p8a(i1) 
           px2=p8a(i1) 
           py1=p8b(j1,i2-na) 
           py2=p8b(j2,i2-na) 
         case(9) 
           px1=p9a(i1) 
           px2=p9a(i1) 
           py1=p9b(j1,i2-na) 
           py2=p9b(j2,i2-na) 
         end select 
      else 
         select case(ion) 
         case(4) 
           px1=p4b(j1,i1-na) 
           px2=p4b(j2,i1-na) 
           py1=p4b(j1,i2-na) 
           py2=p4b(j2,i2-na) 
         case(5) 
           px1=p5b(j1,i1-na) 
           px2=p5b(j2,i1-na) 
           py1=p5b(j1,i2-na) 
           py2=p5b(j2,i2-na) 
         case(6) 
           px1=p6b(j1,i1-na) 
           px2=p6b(j2,i1-na) 
           py1=p6b(j1,i2-na) 
           py2=p6b(j2,i2-na) 
         case(7) 
           px1=p7b(j1,i1-na) 
           px2=p7b(j2,i1-na) 
           py1=p7b(j1,i2-na) 
           py2=p7b(j2,i2-na) 
         case(8) 
           px1=p8b(j1,i1-na) 
           px2=p8b(j2,i1-na) 
           py1=p8b(j1,i2-na) 
           py2=p8b(j2,i2-na) 
         case(9) 
           px1=p9b(j1,i1-na) 
           px2=p9b(j2,i1-na) 
           py1=p9b(j1,i2-na) 
           py2=p9b(j2,i2-na) 
         end select 
      end if 
      dlgunx=px2-px1 
      px=px1+(pne-pn(j1))*dlgunx 
      dlguny=py2-py1 
      py=py1+(pne-pn(j1))*dlguny 
      delt=tt(i2)-tt(i1) 
      if(delt.ne.0.) then 
         dlgut=(py-px)/delt 
         pf=px+(t0-tt(i1))*dlgut 
       else 
         pf=px 
      end if 
      pf=exp(xen*pf) 
      return 
      end subroutine pffe 
 
! 
! ******************************************************************** 
! ******************************************************************** 
! 
 
      SUBROUTINE MATINV(A,N,NR) 
!     ========================= 
! 
!     Matrix inversion 
!      by LU decomposition 
! 
!      A  -  matrix of actual size (N x N) and maximum size (NR x NR) 
!            to be inverted; 
!      Inversion is accomplished in place and the original matrix is 
!      replaced by its inverse 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: A(NR,NR) 
! 
      IF(N.EQ.1) THEN 
         A(1,1)=1.0/A(1,1) 
         RETURN 
      END IF 
! 
      DO I=2,N 
         IM1=I-1 
         DO J=1,IM1 
            JM1=J-1 
            DIV=A(J,J) 
            SUM=0. 
            IF(JM1.GE.1) THEN 
               DO K=1,JM1 
                  SUM=SUM+A(I,K)*A(K,J) 
               END DO 
            END IF 
            A(I,J)=(A(I,J)-SUM)/DIV 
         END DO 
         DO J=I,N 
            SUM=0. 
            DO K=1,IM1 
               SUM=SUM+A(I,K)*A(K,J) 
            END DO 
            A(I,J)=A(I,J)-SUM 
         END DO 
      END DO 
      DO II=2,N 
         I=N+2-II 
         IM1=I-1 
         IF(IM1.GE.1) THEN 
            DO JJ=1,IM1 
               J=I-JJ 
               JP1=J+1 
               SUM=0. 
               IF(JP1.LE.IM1) THEN 
                  DO K=JP1,IM1 
                     SUM=SUM+A(I,K)*A(K,J) 
                  END DO 
               END IF 
               A(I,J)=-A(I,J)-SUM 
            END DO 
         END IF 
      END DO 
      DO II=1,N 
         I=N+1-II 
         DIV=A(I,I) 
         IP1=I+1 
         IF(IP1.LE.N) THEN 
            DO JJ=IP1,N 
               J=N+IP1-JJ 
               SUM=0. 
               DO K=IP1,J 
                  SUM=SUM+A(I,K)*A(K,J) 
               END DO 
               A(I,J)=-SUM/DIV 
            END DO 
         END IF 
         A(I,I)=1.0e0/A(I,I) 
      END DO 
! 
      DO I=1,N 
         DO J=1,I-1 
            K0=I 
            SUM=0. 
            DO K=K0,N 
               SUM=SUM+A(I,K)*A(K,J) 
            END DO 
            A(I,J)=SUM 
         END DO 
         DO J=1,N 
            K0=J 
            SUM=A(I,K0) 
            IF(K0.LT.N) THEN 
               K0=K0+1 
               DO K=K0,N 
                  SUM=SUM+A(I,K)*A(K,J) 
               END DO 
            END IF 
            A(I,J)=SUM 
         END DO 
      END DO 
 
      RETURN 
      END SUBROUTINE MATINV 
! 
! 
!     **************************************************************** 
! 
! 
 
      SUBROUTINE LINEQS(A,B,X,N,NR) 
!     ============================= 
! 
!     Solution of the linear system A*X=B 
!     by Gaussian elimination with partial pivoting 
! 
!     Input: A  - matrix of the linear system, with actual size (N x N), 
!                and maximum size (NR x NR) 
!            B  - the rhs vector 
!     Output: X - solution vector 
!     Note that matrix A and vector B are destroyed here ! 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: A(NR,NR),B(NR),X(NR),D(MLEVEL) 
      INTEGER  :: IP(MLEVEL) 
 
      DO I=1,N 
         DO J=1,N 
            D(J)=A(J,I) 
         END DO 
         IM1=I-1 
         IF(IM1.GE.1) THEN 
            DO J=1,IM1 
               IT=IP(J) 
               A(J,I)=D(IT) 
               D(IT)=D(J) 
               JP1=J+1 
               DO K=JP1,N 
                  D(K)=D(K)-A(K,J)*A(J,I) 
               END DO 
            END DO 
         END IF 
         AM=ABS(D(I)) 
         IP(I)=I 
         DO K=I,N 
            IF(AM.LT.ABS(D(K))) THEN 
               IP(I)=K 
               AM=ABS(D(K)) 
            END IF 
         END DO 
         IT=IP(I) 
         A(I,I)=D(IT) 
         D(IT)=D(I) 
         IP1=I+1 
         IF(IP1.GT.N) EXIT 
         DO K=IP1,N 
            A(K,I)=D(K)/A(I,I) 
         END DO 
      END DO 
 
      DO I=1,N 
         IT=IP(I) 
         X(I)=B(IT) 
         B(IT)=B(I) 
         IP1=I+1 
         IF(IP1.GT.N) EXIT 
         DO J=IP1,N 
            B(J)=B(J)-A(J,I)*X(I) 
         END DO 
      END DO 
 
      DO I=1,N 
         K=N-I+1 
         SUM=0. 
         KP1=K+1 
         IF(KP1.LE.N) THEN 
            DO J=KP1,N 
               SUM=SUM+A(K,J)*X(J) 
            END DO 
         END IF 
         X(K)=(X(K)-SUM)/A(K,K) 
      END DO 
      RETURN 
      END SUBROUTINE LINEQS 
! 
! 
!     **************************************************************** 
! 
! 
 
      FUNCTION EXPINT(X) 
!     ================== 
! 
!     First exponential integral function E1(X) 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
! 
      IF(X.LE.1.0) THEN 
         EXPINT=-LOG(X)-0.57721566+X*(0.99999193+X*(-0.24991055           & 
     &          +X*(0.05519968+X*(-0.00976004+X*0.00107857)))) 
       ELSE 
         EXPINT=EXP(-X)*((0.2677734343+X*(8.6347608925+X*                 & 
     &          (18.059016973+X*(8.5733287401+X))))/                      & 
     &          (3.9584969228+X*(21.0996530827+X*                         & 
     &          (25.6329561486+X*(9.5733223454+X)))))/X 
      END IF 
      RETURN 
      END FUNCTION EXPINT 
! 
! 
!     **************************************************************** 
! 
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
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) ::  X(NX),Y(NX),XX(NXX),YY(NXX) 
      EXP10(X0)=EXP(X0*2.30258509299405) 
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
      IF(ILOGX.NE.0) THEN 
         DO I=1,NX 
            X(I)=LOG10(X(I)) 
         END DO 
         DO I=1,NXX 
            XX(I)=LOG10(XX(I)) 
         END DO 
      END IF 
      IF(ILOGY.NE.0) THEN 
         DO I=1,NX 
            Y(I)=LOG10(Y(I)) 
         END DO 
      END IF 
 
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
               IF(K.NE.M) T=T*(XXX-X(M))/(X(K)-X(M)) 
            END DO 
            YYY=Y(K)*T+YYY 
         END DO 
         YY(ID)=YYY 
      END DO 
      IF(ILOGX.NE.0) THEN 
         DO I=1,NX 
            X(I)=EXP10(X(I)) 
         END DO 
         DO I=1,NXX 
            XX(I)=EXP10(XX(I)) 
         END DO 
      END IF 
      IF(ILOGY.NE.0) THEN 
         DO I=1,NX 
            Y(I)=EXP10(Y(I)) 
         END DO 
         DO I=1,NXX 
           YY(I)=EXP10(YY(I)) 
         END DO 
      END IF 
 
      RETURN 
      END SUBROUTINE INTERP 
! 
! ******************************************************************** 
! 
 
      subroutine intrp(wltab,absop,wlgrid,abgrd,nfr,nfgrid) 
!     ===================================================== 
! 
!     a more efficient interpolation routine - using bisection 
! 
      use accura 
      use params 
      use optabl, only : yint,jint 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      real(dp) :: wltab(nfr),absop(nfr),                                  & 
     &            wlgrid(nfgrid),abgrd(nfgrid) 
! 
!      set up interpolation coefficients for an interpolation 
!      by bisection 
! 
         fr1=wltab(1) 
         fr2=wltab(nfr) 
         do ij=1,nfgrid 
            xint=wlgrid(ij) 
            jl=0 
            ju=nfr+1 
            sele: do 
               if(ju-jl.gt.1) then 
                  jm=(ju+jl)/2 
                  if((fr2.gt.fr1).eqv.(xint.gt.wltab(jm))) then 
                     jl=jm 
                   else 
                     ju=jm 
                  end if 
                  cycle sele 
                else 
                  exit sele 
               end if 
            end do sele 
            j=jl 
            if(j.eq.nfr) j=j-1 
            if(j.eq.0) j=j+1 
            jint(ij)=j 
!           yint(ij)=un/log10(wltab(j+1)/wltab(j)) 
            yint(ij)=1./(wltab(j+1)-wltab(j)) 
         end do 
! 
         do ij=1,nfgrid 
            j=jint(ij) 
            rc=(absop(j+1)-absop(j))*yint(ij) 
!           abgrd(ij)=rc*log10(wlgrid(ij)/wltab(j))+absop(j) 
            abgrd(ij)=rc*(wlgrid(ij)-wltab(j))+absop(j) 
         end do 
      return 
      end subroutine intrp 
 
! 
! ******************************************************************** 
! 
      SUBROUTINE PFSPEC(IAT,IZI,T,ANE,U) 
!     ================================== 
 
!     Non-standard evaluation of the partition function 
!     user supplied procedure 
! 
!     Input: 
!      IAT   - atomic number 
!      IZI   - ionic charge (=1 for neutrals, =1 for once ionized, etc) 
!      T     - temperature 
!      ANE   - electron density 
!      XMAX  - principal quantum number of the last bound level 
! 
!     Output: 
!      U     - partition function 
! 
! 
! Modified from the ATMOS related programme 5-April-1990 
! as in addition to TLUSTY to allow high ionisation states 
! of C, N and O 
! 
! M.A Barstow - University of Leicester, Dept of Physics & Astronomy 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      INTEGER, PARAMETER :: 
     &     MH=100,MHEI=100,MHEII=100,MCI=135,                             & 
     &     MCII=157,MCIII=156,MCIV=55,MCV=15,MCVI=100,MNI=228,MNII=122,   & 
     &     MNIII=133,MNIV=73,MNV=51,MNVI=8,MNVII=100,MOI=174,MOII=191,    & 
     &     MOIII=168,MOIV=166,MOV=115,MOVI=52,MOVII=16,MOVIII=100 
      REAL(DP) :: GHYD(MH),SHYD(MH),ENHYD(MH),                            & 
     &     GHEL(MH),ENHEL(MH),SHEL(MH),                                   & 
     &     GCI(MCI),ENCI(MCI),SCI(MCI),                                   & 
     &     GCII(MCII),ENCII(MCII),SCII(MCII),                             & 
     &     GCIII(MCIII),ENCIII(MCIII),SCIII(MCIII),                       & 
     &     GCIV(MCIV),ENCIV(MCIV),SCIV(MCIV),                             & 
     &     GCV(MCV),ENCV(MCV),SCV(MCV),                                   & 
     &     GNI(MNI),ENNI(MNI),SNI(MNI),                                   & 
     &     GNII(MNII),ENNII(MNII),SNII(MNII),                             & 
     &     GNIII(MNIII),ENNIII(MNIII),SNIII(MNIII),                       & 
     &     GNIV(MNIV),ENNIV(MNIV),SNIV(MNIV),                             & 
     &     GNV(MNV),ENNV(MNV),SNV(MNV),                                   & 
     &     GNVI(MNVI),ENNVI(MNVI),SNVI(MNVI),                             & 
     &     GOI(MOI),ENOI(MOI),SOI(MOI),                                   & 
     &     GOII(MOII),ENOII(MOII),SOII(MOII),                             & 
     &     GOIII(MOIII),ENOIII(MOIII),SOIII(MOIII),                       & 
     &     GOIV(MOIV),ENOIV(MOIV),SOIV(MOIV),                             & 
     &     GOV(MOV),ENOV(MOV),SOV(MOV),                                   & 
     &     GOVI(MOVI),ENOVI(MOVI),SOVI(MOVI),                             & 
     &     GOVII(MOVII),ENOVII(MOVII),SOVII(MOVII) 
      INTEGER :: NHYD(MH),NHEL(MHEI),NCI(MCI),NCII(MCII),                 & 
     &     NCIII(MCIII),NCIV(MCIV),NCV(MCV),NNI(MNI),NNII(MNII),          & 
     &     NNIII(MNIII),NNIV(MNIV),NNV(MNV),NNVI(MNVI),NOI(MOI),          & 
     &     NOII(MOII),NOIII(MOIII),NOIV(MOIV),NOV(MOV),NOVI(MOVI),        & 
     &     NOVII(MOVII) 
      REAL(DP), PARAMETER :: HI=13.5878,HEI=24.587,HEII=54.416,           & 
     &     CVI=489.84,NVII=666.83,OVIII=871.12 
      REAL(DP), PARAMETER :: ZH=1.0,ZHE=2.0,ZC=6.0,ZN=7.0,ZO=8.0 
 
!                           N***=QUANTUM NO. OF LEVEL 
!      DATA FOR IONS        G***=STATISTICAL WEIGHT OF LEVEL 
!                           EN***=ENERGY OF LEVEL 
!                           S*=SCREENING NO. OF LEVEL 
 
        DATA NHYD/ 1, 2, 3, 4, 5, 6,                                      & 
     &          7, 8, 9,10,11,12,                                         & 
     &          13,14,15,16,17,18,                                        & 
     &          19,20,21,22,23,24,                                        & 
     &          25,26,27,28,29,30,                                        & 
     &          31,32,33,34,35,36,                                        & 
     &          37,38,39,40,41,42,                                        & 
     &          43,44,45,46,47,48,                                        & 
     &          49,50,51,52,53,54,                                        & 
     &          55,56,57,58,59,60,                                        & 
     &          61,62,63,64,65,66,                                        & 
     &          67,68,69,70,71,72,                                        & 
     &          73,74,75,76,77,78,                                        & 
     &          79,80,81,82,83,84,                                        & 
     &          85,86,87,88,89,90,                                        & 
     &          91,92,93,94,95,96,                                        & 
     &          97,98,99, 100 / 
        DATA GHYD/ 2.000000, 8.000000, 18.00000,                          & 
     &           32.00000, 50.00000, 72.00000,                            & 
     &           98.00000, 128.0000, 162.0000,                            & 
     &           200.0000, 242.0000, 288.0000,                            & 
     &           338.0000, 392.0000, 450.0000,                            & 
     &           512.0000, 578.0000, 648.0000,                            & 
     &           722.0000, 800.0000, 882.0000,                            & 
     &           968.0000, 1058.000, 1152.000,                            & 
     &           1250.000, 1352.000, 1458.000,                            & 
     &           1568.000, 1682.000, 1800.000,                            & 
     &           1922.000, 2048.000, 2178.000,                            & 
     &           2312.000, 2450.000, 2592.000,                            & 
     &           2738.000, 2888.000, 3042.000,                            & 
     &           3200.000, 3362.000, 3528.000,                            & 
     &           3698.000, 3872.000, 4050.000,                            & 
     &           4232.000, 4418.000, 4608.000,                            & 
     &           4802.000, 5000.000, 5202.000,                            & 
     &           5408.000, 5618.000, 5832.000,                            & 
     &           6050.000, 6272.000, 6498.000,                            & 
     &           6728.000, 6962.000, 7200.000,                            & 
     &           7442.000, 7688.000, 7938.000,                            & 
     &           8192.000, 8450.000, 8712.000,                            & 
     &           8978.000, 9248.000, 9522.000,                            & 
     &           9800.000, 10082.00, 10368.00,                            & 
     &           10658.00, 10952.00, 11250.00,                            & 
     &           11552.00, 11858.00, 12168.00,                            & 
     &           12482.00, 12800.00, 13122.00,                            & 
     &           13448.00, 13778.00, 14112.00,                            & 
     &           14450.00, 14792.00, 15138.00,                            & 
     &           15488.00, 15842.00, 16200.00,                            & 
     &           16562.00, 16928.00, 17298.00,                            & 
     &           17672.00, 18050.00, 18432.00,                            & 
     &           18818.00, 19208.00, 19602.00,                            & 
     &           20000.00/ 
        DATA ENHYD /0.0000000E+00,10.19085000000000,12.07804444444444,    & 
     &           12.73856250000000,13.04428800000000,13.21036111111111,   & 
     &           13.31049795918367,13.37549062500000,13.42004938271605,   & 
     &           13.45192200000000,13.47550413223140,13.49344027777778,   & 
     &           13.50739881656805,13.51847448979592,13.52740977777778,   & 
     &           13.53472265625000,13.54078339100346,13.54586234567901,   & 
     &           13.55016066481994,13.55383050000000,13.55698866213152,   & 
     &           13.55972603305785,13.56211417769376,13.56421006944444,   & 
     &           13.56605952000000,13.56769970414201,13.56916104252401,   & 
     &           13.57046862244898,13.57164328180737,13.57270244444444,   & 
     &           13.57366077003122,13.57453066406250,13.57532268135905,   & 
     &           13.57604584775087,13.57670791836735,13.57731558641975,   & 
     &           13.57787465303141,13.57839016620499,13.57886653517423,   & 
     &           13.57930762500000,13.57971683521713,13.58009716553288,   & 
     &           13.58045127095727,13.58078150826446,13.58108997530864,   & 
     &           13.58137854442344,13.58164889090086,13.58190251736111,   & 
     &           13.58214077467722,13.58236488000000,13.58257593233372,   & 
     &           13.58277492603550,13.58296276254895,13.58314026063100,   & 
     &           13.58330816528926,13.58346715561225,13.58361785164666,   & 
     &           13.58376082045184,13.58389658144211,13.58402561111111,   & 
     &           13.58414834721849,13.58426519250780,13.58437651801461,   & 
     &           13.58448266601563,13.58458395266272,13.58468067033976,   & 
     &           13.58477308977501,13.58486146193772,13.58494601974375,   & 
     &           13.58502697959184,13.58510454274945,13.58517889660494,   & 
     &           13.58525021580034,13.58531866325785,13.58538439111111,   & 
     &           13.58544754155125,13.58550824759656,13.58556663379356,   & 
     &           13.58562281685627,13.58567690625000,13.58572900472489,   & 
     &           13.58577920880428,13.58582760923211,13.58587429138322,   & 
     &           13.58591933564014,13.58596281773932,13.58600480908971,   & 
     &           13.58604537706612,13.58608458527964,13.58612249382716,   & 
     &           13.58615915952180,13.58619463610586,13.58622897444791,   & 
     &           13.58626222272522,13.58629442659280,13.58632562934028,   & 
     &           13.58635587203741,13.58638519366930,13.58641363126212,   & 
     &           13.58644122000000/ 
        DATA SHYD/100*0.0/ 
      DATA NHEL/1,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,            & 
     &        5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,           & 
     &        23,24,25,26,27,                                             & 
     &          28,29,30,31,32,33,                                        & 
     &          34,35,36,37,38,39,                                        & 
     &          40,41,42,43,44,45,                                        & 
     &          46,47,48,49,50,51,                                        & 
     &          52,53,54,55,56,57,                                        & 
     &          58,59,60,61,62,63,                                        & 
     &          64,65,66,67,68,69,                                        & 
     &          70,71,72,73,74,75,                                        & 
     &          76,77,78,79,80,81/ 
      DATA GHEL/1.0  ,3.0  ,1.0  ,5.0  ,3.0  ,1.0  ,3.0  ,                & 
     &          3.0  ,1.0  ,5.0  ,3.0  ,                                  & 
     &         1.0  ,15.0  ,5.0  ,3.0  ,3.0  ,1.0  ,9.0  ,                & 
     &         15.0  ,5.0  ,21.0  ,7.0  ,                                 & 
     &         3.0  ,100.0  ,144.0  ,196.0  ,256.0  ,324.0  ,             & 
     &         400.0  ,484.0  ,                                           & 
     &         576.0  ,676.0  ,784.0  ,900.0  ,1024.0  ,1156.0  ,         & 
     &         1296.0  ,1444.0  ,1600.0  ,1764.0  ,1936.0  ,              & 
     &         2116.0  ,2304.0  ,2500.0  ,2704.0  ,3136.0  ,              & 
     &          3136.000000000000,3364.000000000000,3600.000000000000,    & 
     &          3844.000000000000,4096.000000000000,4356.000000000000,    & 
     &          4624.000000000000,4900.000000000000,5184.000000000000,    & 
     &          5476.000000000000,5776.000000000000,6084.000000000000,    & 
     &          6400.000000000000,6724.000000000000,7056.000000000000,    & 
     &          7396.000000000000,7744.000000000000,8100.000000000000,    & 
     &          8464.000000000000,8836.000000000000,9216.000000000000,    & 
     &          9604.000000000000,10000.00000000000,10404.00000000000,    & 
     &          10816.00000000000,11236.00000000000,11664.00000000000,    & 
     &          12100.00000000000,12544.00000000000,12996.00000000000,    & 
     &          13456.00000000000,13924.00000000000,14400.00000000000,    & 
     &          14884.00000000000,15376.00000000000,15876.00000000000,    & 
     &          16384.00000000000,16900.00000000000,17424.00000000000,    & 
     &          17956.00000000000,18496.00000000000,19044.00000000000,    & 
     &          19600.00000000000,20164.00000000000,20736.00000000000,    & 
     &          21316.00000000000,21904.00000000000,22500.00000000000,    & 
     &          23104.00000000000,23716.00000000000,24336.00000000000,    & 
     &          24964.00000000000,25600.00000000000,26244.00000000000/ 
      DATA ENHEL/0.0  ,19.819  ,20.615  ,20.964  ,                        & 
     &           20.964  ,20.964  ,21.218  ,                              & 
     &           22.718  ,22.920  ,23.007  ,23.007  ,                     & 
     &           23.007  ,23.073  ,23.074  ,                              & 
     &           23.087  ,23.593  ,23.673  ,23.707  ,                     & 
     &           23.736  ,23.736  ,23.737  ,                              & 
     &           23.737  ,23.742  ,24.028  ,24.201  ,                     & 
     &           24.304  ,24.371  ,24.417  ,                              & 
     &           24.449  ,24.473  ,24.491  ,24.506  ,                     & 
     &           24.517  ,24.526  ,24.534  ,                              & 
     &           24.540  ,24.545  ,24.549  ,24.553  ,                     & 
     &           24.556  ,24.559  ,24.562  ,                              & 
     &           24.564  ,24.566  ,24.568  ,24.570  ,                     & 
     &          24.57131951530612,24.57238228299643,24.57334055555556,    & 
     &          24.57420759625390,24.57499462890625,24.57571120293848,    & 
     &          24.57636548442907,24.57696448979592,24.57751427469136,    & 
     &          24.57802008765522,24.57848649584488,24.57891748849441,    & 
     &          24.57931656250000,24.57968679357525,24.58003089569161,    & 
     &          24.58035127095727,24.58065005165289,24.58092913580247,    & 
     &          24.58119021739130,24.58143481213219,24.58166427951389,    & 
     &          24.58187984173261,24.58208260000000,24.58227354863514,    & 
     &          24.58245358727811,24.58262353150587,24.58278412208505,    & 
     &          24.58293603305785,24.58307987882653,24.58321622037550,    & 
     &          24.58334557074911,24.58346839988509,24.58358513888889,    & 
     &          24.58369618382155,24.58380189906348,24.58390262030738,    & 
     &          24.58399865722656,24.58409029585799,24.58417780073462,    & 
     &          24.58426141679661,24.58434137110727,24.58441787439614,    & 
     &          24.58449112244898,24.58456129736163,24.58462856867284,    & 
     &          24.58469309438919,24.58475502191381,24.58481448888889,    & 
     &          24.58487162396122,24.58492654747850,24.58497937212360,    & 
     &          24.58503020349303,24.58507914062500,24.58512627648224/ 
      DATA SHEL/0.375  ,0.622  ,0.622  ,0.842  ,                          & 
     &          0.842  ,0.842  ,0.842  ,0.747  ,                          & 
     &          0.747  ,0.912  ,0.912  ,0.912  ,                          & 
     &          0.993  ,0.993  ,0.912  ,0.810  ,                          & 
     &          0.810  ,0.937  ,0.995  ,0.995  ,                          & 
     &          1.000  ,1.000  ,0.937  ,0.949  ,                          & 
     &          0.958  ,75*1.000  / 
        DATA NCI/2,2,2,2,2,2,3,3,3,3,2,2,2,3,3,3,3,3,                     & 
     &          3,3,3,3,3,2,3,4,4,4,3,3,3,3,3,3,4,3,                      & 
     &          3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,                      & 
     &          4,4,4,5,4,4,4,4,4,5,5,5,5,5,5,5,5,5,                      & 
     &          5,5,5,5,6,5,5,5,5,5,6,6,6,6,6,6,6,7,                      & 
     &          6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,8,8,                      & 
     &          8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,10,10,                      & 
     &          10,11,11,11,2,3,3,3,2,2/ 
        DATA GCI/1.0  ,3.0  ,5.0  ,5.0  ,1.0  ,                           & 
     &          5.0  ,1.0  ,3.0  ,5.0  ,3.0  ,                            & 
     &          7.0  ,5.0  ,3.0  ,3.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,3.0  ,1.0  ,3.0  ,                            & 
     &          5.0  ,5.0  ,1.0  ,9.0  ,5.0  ,                            & 
     &          1.0  ,3.0  ,5.0  ,5.0  ,7.0  ,                            & 
     &          9.0  ,3.0  ,5.0  ,7.0  ,3.0  ,                            & 
     &          3.0  ,3.0  ,5.0  ,3.0  ,1.0  ,                            & 
     &          3.0  ,5.0  ,7.0  ,3.0  ,3.0  ,                            & 
     &          1.0  ,3.0  ,5.0  ,5.0  ,1.0  ,                            & 
     &          5.0  ,5.0  ,7.0  ,9.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,3.0  ,7.0  ,3.0  ,                            & 
     &          5.0  ,3.0  ,1.0  ,3.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,5.0  ,1.0  ,5.0  ,                            & 
     &          5.0  ,7.0  ,9.0  ,3.0  ,5.0  ,                            & 
     &          7.0  ,3.0  ,7.0  ,3.0  ,5.0  ,                            & 
     &          3.0  ,1.0  ,5.0  ,5.0  ,7.0  ,                            & 
     &          9.0  ,3.0  ,5.0  ,7.0  ,3.0  ,                            & 
     &          7.0  ,5.0  ,3.0  ,1.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,9.0  ,3.0  ,5.0  ,                            & 
     &          7.0  ,7.0  ,3.0  ,5.0  ,3.0  ,                            & 
     &          1.0  ,9.0  ,7.0  ,5.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,7.0  ,5.0  ,3.0  ,                            & 
     &          1.0  ,9.0  ,7.0  ,5.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,7.0  ,3.0  ,5.0  ,                            & 
     &          7.0  ,3.0  ,5.0  ,7.0  ,5.0  ,                            & 
     &          3.0  ,5.0  ,7.0  ,3.0  ,3.0  / 
        DATA ENCI/0.0  ,2.0333605e-03,5.3933649e-03,1.263870,2.684086,    & 
     &          4.182672,7.480511,7.482891,7.487915,7.684888,             & 
     &          7.946046,7.946620,7.946474,8.537387,8.640516,             & 
     &          8.643146,8.647287,8.771255,8.846707,8.848247,             & 
     &          8.850785,9.002712,9.171972,9.330682,9.631248,             & 
     &          9.683908,9.685375,9.689256,9.695577,9.697620,             & 
     &          9.701885,9.708156,9.708925,9.710041,9.712769,             & 
     &          9.714380,9.761111,9.833419,9.834406,9.834934,             & 
     &          9.940317,9.942698,9.946449,9.988707,10.05592,             & 
     &          10.08144,10.08328,10.08553,10.13833,10.19809,             & 
     &          10.35278,10.38514,10.38514,10.38514,10.39370,             & 
     &          10.39456,10.39580,10.40021,10.40845,10.41874,             & 
     &          10.42750,10.42990,10.42990,10.52043,10.52041,             & 
     &          10.52041,10.53705,10.58840,10.61635,10.67973,             & 
     &          10.70230,10.70328,10.70328,10.70878,10.70878,             & 
     &          10.71184,10.71407,10.71854,10.72362,10.72523,             & 
     &          10.72684,10.72684,10.86509,10.87426,10.87513,             & 
     &          10.87513,10.87997,10.87997,10.88257,10.88533,             & 
     &          10.88679,10.88964,10.89075,10.89075,10.88980,             & 
     &          10.97789,10.97854,10.97854,10.98597,10.98597,             & 
     &          10.98597,10.98808,10.98913,10.98994,10.98994,             & 
     &          10.98994,11.04474,11.04474,11.04487,11.05280,             & 
     &          11.05280,11.05280,11.05392,11.05429,11.05429,             & 
     &          11.05429,11.09049,11.09049,11.09049,11.09843,             & 
     &          11.09843,11.09843,11.09880,11.13129,11.13129,             & 
     &          11.13129,11.15477,11.15477,11.15477,12.13544,             & 
     &          12.83767,12.84024,12.84331,13.11772,14.86312/ 
        DATA NCII/2,2,2,2,2,2,2,2,2,2,3,3,3,2,3,3,2,2,4,4,4,3,3,3,        & 
     &          4,4,2,2,4,4,5,5,5,3,3,5,5,5,5,6,3,3,3,3,3,3,6,6,          & 
     &          6,6,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,          & 
     &          3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,          & 
     &          4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,          & 
     &          4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,          & 
     &          5,5,5,5,5,5,6,6,6,6,6,6,6/ 
        DATA GCII/2.0  ,4.0  ,2.0  ,4.0  ,6.0  ,                          & 
     &          6.0  ,4.0  ,2.0  ,2.0  ,4.0  ,                            & 
     &          2.0  ,2.0  ,4.0  ,4.0  ,4.0  ,                            & 
     &          6.0  ,6.0  ,4.0  ,2.0  ,2.0  ,                            & 
     &          4.0  ,2.0  ,4.0  ,6.0  ,4.0  ,                            & 
     &          6.0  ,2.0  ,4.0  ,6.0  ,8.0  ,                            & 
     &          2.0  ,2.0  ,4.0  ,2.0  ,4.0  ,                            & 
     &          4.0  ,6.0  ,6.0  ,8.0  ,2.0  ,                            & 
     &          2.0  ,4.0  ,6.0  ,8.0  ,2.0  ,                            & 
     &          4.0  ,4.0  ,6.0  ,6.0  ,8.0  ,                            & 
     &          4.0  ,2.0  ,4.0  ,6.0  ,4.0  ,                            & 
     &          6.0  ,2.0  ,4.0  ,6.0  ,8.0  ,                            & 
     &          10.0  ,2.0  ,4.0  ,6.0  ,8.0  ,                           & 
     &          4.0  ,6.0  ,6.0  ,4.0  ,2.0  ,                            & 
     &          6.0  ,8.0  ,4.0  ,2.0  ,2.0  ,                            & 
     &          4.0  ,6.0  ,2.0  ,4.0  ,2.0  ,                            & 
     &          4.0  ,6.0  ,8.0  ,4.0  ,2.0  ,                            & 
     &          4.0  ,6.0  ,4.0  ,6.0  ,4.0  ,                            & 
     &          6.0  ,8.0  ,10.0  ,2.0  ,4.0  ,                           & 
     &          6.0  ,8.0  ,4.0  ,6.0  ,6.0  ,                            & 
     &          4.0  ,2.0  ,6.0  ,8.0  ,4.0  ,                            & 
     &          6.0  ,8.0  ,10.0  ,6.0  ,8.0  ,                           & 
     &          6.0  ,8.0  ,10.0  ,12.0  ,8.0  ,                          & 
     &          10.0  ,8.0  ,6.0  ,4.0  ,2.0  ,                           & 
     &          6.0  ,4.0  ,4.0  ,2.0  ,2.0  ,                            & 
     &          4.0  ,6.0  ,2.0  ,4.0  ,2.0  ,                            & 
     &          4.0  ,6.0  ,8.0  ,6.0  ,4.0  ,                            & 
     &          2.0  ,6.0  ,8.0  ,4.0  ,6.0  ,                            & 
     &          8.0  ,10.0  ,6.0  ,8.0  ,10.0  ,                          & 
     &          12.0  ,8.0  ,6.0  ,4.0  ,2.0  ,                           & 
     &          2.0  ,4.0  ,6.0  ,8.0  ,6.0  ,                            & 
     &          4.0  ,2.0  / 
        DATA ENCII/0.0  ,7.9350658e-03,5.331397,5.334075,5.337658,        & 
     &          9.290338,9.290624,11.96386,13.71590,13.72101,             & 
     &          14.44900,16.33194,16.33332,17.60895,18.04607,             & 
     &          18.04625,18.65519,18.65582,19.49478,20.14995,             & 
     &          20.15068,20.70119,20.70413,20.70971,20.84491,             & 
     &          20.84496,20.92025,20.92256,20.95094,20.95094,             & 
     &          21.49265,21.73314,21.73405,22.09347,22.13075,             & 
     &          22.13075,22.13075,22.18799,22.18799,22.47211,             & 
     &          22.52747,22.52929,22.53239,22.53689,22.56844,             & 
     &          22.57086,22.82136,22.82136,22.85996,22.85996,             & 
     &          22.89870,23.11398,23.11600,23.11878,23.38108,             & 
     &          23.38522,24.12408,24.27024,24.27201,24.27444,             & 
     &          24.27787,24.37010,24.37079,24.37187,24.37315,             & 
     &          24.60198,24.60332,24.65351,24.65617,24.65793,             & 
     &          24.78982,24.79512,25.06741,25.07039,25.98117,             & 
     &          25.98415,25.98986,26.58329,26.58615,26.62689,             & 
     &          26.62867,26.63139,26.63554,26.75178,26.82771,             & 
     &          26.82771,26.83016,26.89454,26.89578,27.22147,             & 
     &          27.22329,27.22585,27.22930,27.29263,27.29263,             & 
     &          27.29378,27.29509,27.35131,27.35294,27.37703,             & 
     &          27.37957,27.38104,27.41188,27.41302,27.41395,             & 
     &          27.41395,27.41395,27.41409,27.46301,27.46301,             & 
     &          27.46810,27.46936,27.47200,27.47561,27.47330,             & 
     &          27.47864,27.48713,27.49096,27.49330,27.49330,             & 
     &          27.48854,27.49412,27.55688,27.56022,27.99752,             & 
     &          27.99752,27.99752,28.25640,28.25640,28.61124,             & 
     &          28.61124,28.61124,28.61124,28.64683,28.64683,             & 
     &          28.64683,28.66803,26.43629,28.66875,28.66875,             & 
     &          28.66875,28.66875,28.70253,28.70253,28.70253,             & 
     &          28.70253,28.70515,28.70515,28.70515,28.70515,             & 
     &          29.31561,29.31561,29.31561,29.31561,29.33557,             & 
     &          29.33557,29.33557/ 
        DATA NCIII/2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,4,       & 
     &          3,4,4,4,4,3,4,4,4,4,4,4,4,4,3,3,3,4,3,3,3,3,3,3,          & 
     &          3,3,3,3,3,3,5,3,3,3,3,5,5,5,5,3,5,5,5,5,5,5,5,5,          & 
     &          5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,          & 
     &          7,8,8,8,8,9,9,9,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,          & 
     &          4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,          & 
     &          6,6,6,6,6,6,7,7,7,7,7,7/ 
        DATA GCIII/1.0  ,1.0  ,3.0  ,5.0  ,3.0  ,                         & 
     &          1.0  ,3.0  ,5.0  ,5.0  ,1.0  ,                            & 
     &          3.0  ,1.0  ,3.0  ,1.0  ,3.0  ,                            & 
     &          5.0  ,3.0  ,5.0  ,7.0  ,5.0  ,                            & 
     &          1.0  ,3.0  ,5.0  ,3.0  ,3.0  ,                            & 
     &          1.0  ,1.0  ,3.0  ,5.0  ,3.0  ,                            & 
     &          3.0  ,5.0  ,7.0  ,5.0  ,7.0  ,                            & 
     &          9.0  ,3.0  ,7.0  ,3.0  ,5.0  ,                            & 
     &          7.0  ,5.0  ,3.0  ,1.0  ,3.0  ,                            & 
     &          5.0  ,5.0  ,5.0  ,5.0  ,7.0  ,                            & 
     &          9.0  ,3.0  ,5.0  ,7.0  ,3.0  ,                            & 
     &          5.0  ,3.0  ,1.0  ,7.0  ,3.0  ,                            & 
     &          1.0  ,3.0  ,5.0  ,1.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,7.0  ,9.0  ,11.0  ,                           & 
     &          9.0  ,5.0  ,3.0  ,5.0  ,7.0  ,                            & 
     &          9.0  ,7.0  ,3.0  ,3.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,7.0  ,9.0  ,11.0  ,                           & 
     &          9.0  ,5.0  ,5.0  ,7.0  ,9.0  ,                            & 
     &          7.0  ,3.0  ,3.0  ,3.0  ,5.0  ,                            & 
     &          7.0  ,5.0  ,3.0  ,3.0  ,5.0  ,                            & 
     &          7.0  ,3.0  ,5.0  ,7.0  ,1.0  ,                            & 
     &          3.0  ,5.0  ,3.0  ,3.0  ,5.0  ,                            & 
     &          7.0  ,1.0  ,3.0  ,5.0  ,5.0  ,                            & 
     &          5.0  ,3.0  ,5.0  ,7.0  ,5.0  ,                            & 
     &          3.0  ,1.0  ,7.0  ,3.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,1.0  ,3.0  ,5.0  ,                            & 
     &          5.0  ,5.0  ,3.0  ,5.0  ,7.0  ,                            & 
     &          5.0  ,3.0  ,1.0  ,3.0  ,5.0  ,                            & 
     &          7.0  ,1.0  ,3.0  ,5.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,5.0  ,3.0  ,1.0  ,                            & 
     &          3.0  ,5.0  ,7.0  ,1.0  ,3.0  ,5.0  / 
        DATA ENCIII/0.0  ,6.486296,6.489148,6.496191,12.69008,            & 
     &          17.03237,17.03602,17.04185,18.08638,22.62984,             & 
     &          29.52845,30.64541,32.10371,32.19328,32.19396,             & 
     &          32.19555,33.47080,33.45866,33.47146,34.27982,             & 
     &          38.20770,38.21183,38.22034,38.36164,38.43612,             & 
     &          38.64882,39.39549,39.39549,39.39611,39.64054,             & 
     &          39.84380,39.84582,39.84874,39.91699,39.91782,             & 
     &          39.91892,39.97328,40.01022,40.05026,40.05341,             & 
     &          40.05822,40.19756,40.57121,40.86969,40.87231,             & 
     &          40.87686,41.24874,41.30157,41.32848,41.33158,             & 
     &          41.33611,41.85783,41.80309,41.86202,42.14028,             & 
     &          42.16117,42.16444,42.16623,42.32471,42.55869,             & 
     &          42.67342,42.67342,42.67342,42.78661,42.83001,             & 
     &          42.83001,42.83001,42.96405,42.96405,42.96416,             & 
     &          42.96405,42.98029,42.98736,43.03527,43.03550,             & 
     &          43.03579,43.25349,43.98952,44.27370,44.39248,             & 
     &          44.39248,44.39248,44.46592,44.46592,44.46600,             & 
     &          44.47219,44.47673,44.48596,44.48596,44.48596,             & 
     &          44.52591,45.07626,45.24178,45.32720,45.32720,             & 
     &          45.32720,45.38200,45.86543,45.92891,45.92891,             & 
     &          45.92891,46.33929,46.33929,46.33929,46.69749,             & 
     &          46.69749,46.69749,47.25143,47.35238,47.35238,             & 
     &          47.35722,47.64920,47.64920,47.65379,47.81342,             & 
     &          47.83558,48.06245,48.06245,48.06245,48.16114,             & 
     &          48.16114,48.16114,48.20208,50.51542,50.55803,             & 
     &          50.55803,50.55803,50.69428,50.69428,50.69428,             & 
     &          50.77264,50.79460,50.90022,50.90022,50.90022,             & 
     &          50.93829,50.93829,50.93829,52.24497,52.24497,             & 
     &          52.24497,52.31775,52.31775,52.31775,52.43107,             & 
     &          52.43107,52.43107,52.45302,52.45302,52.45302,             & 
     &          53.23251,53.23251,53.23251,53.27802,53.27802,             & 
     &          53.27802/ 
        DATA NCIV/2,2,2,3,3,3,3,3,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,        & 
     &          6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,8,8,          & 
     &          8,8,8,8,8,8,8/ 
        DATA GCIV/2.0  ,2.0  ,4.0  ,2.0  ,2.0  ,                          & 
     &          4.0  ,4.0  ,6.0  ,2.0  ,2.0  ,                            & 
     &          4.0  ,4.0  ,6.0  ,6.0  ,8.0  ,                            & 
     &          2.0  ,2.0  ,4.0  ,4.0  ,6.0  ,                            & 
     &          6.0  ,8.0  ,8.0  ,10.0  ,2.0  ,                           & 
     &          2.0  ,4.0  ,4.0  ,6.0  ,6.0  ,                            & 
     &          8.0  ,8.0  ,10.0  ,10.0  ,12.0  ,                         & 
     &          2.0  ,2.0  ,4.0  ,4.0  ,6.0  ,                            & 
     &          6.0  ,8.0  ,8.0  ,10.0  ,10.0  ,                          & 
     &          12.0  ,2.0  ,4.0  ,6.0  ,8.0  ,                           & 
     &          8.0  ,10.0  ,12.0  ,14.0  ,16.0  / 
        DATA ENCIV/0.0  ,7.995100,8.008378,37.54872,39.68134,             & 
     &          39.68525,40.28040,40.28173,49.76113,50.62434,             & 
     &          50.62599,50.87540,50.87595,50.88784,50.88784,             & 
     &          55.21889,55.65134,55.65221,55.77947,55.77947,             & 
     &          55.78577,55.78578,55.78703,55.78703,58.12002,             & 
     &          58.36774,58.36774,58.44275,58.44275,58.44709,             & 
     &          58.44709,58.44764,58.44764,58.44770,58.44770,             & 
     &          59.84267,60.00038,60.00038,60.04725,60.04725,             & 
     &          60.05156,60.05156,60.05191,60.05191,60.05194,             & 
     &          60.05194,61.05946,61.05946,61.09294,61.09294,             & 
     &          61.09319,61.09319,61.09319,61.09319,61.09319/ 
        DATA NCV/1,2,2,2,2,2,3,3,3,3,4,5,6,7,8/ 
        DATA GCV/1.0  ,3.0  ,1.0  ,3.0  ,5.0  ,                           & 
     &          3.0  ,3.0  ,5.0  ,7.0  ,3.0  ,                            & 
     &          3.0  ,3.0  ,3.0  ,3.0  ,3.0  / 
        DATA ENCV/0.0  ,298.9618,304.4046,304.4030,304.4199,              & 
     &          307.8855,354.2645,354.2645,354.2645,354.5177,             & 
     &          370.9247,378.5349,382.6710,385.1917,386.6807/ 
        DATA NNI/2,2,2,2,2,3,3,3,3,3,2,2,2,3,3,3,3,3,3,3,3,3,3,3,         & 
     &          3,3,3,3,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,          & 
     &          3,3,4,4,4,4,4,4,4,4,4,5,5,5,5,5,4,4,4,4,4,4,4,4,          & 
     &          4,4,4,4,4,4,4,4,4,3,3,3,3,6,6,6,6,6,5,5,5,5,5,5,          & 
     &          5,5,5,5,5,5,5,5,5,5,5,7,7,7,7,7,6,6,6,6,6,6,6,6,          & 
     &          6,6,6,6,6,6,6,6,6,8,8,8,8,8,7,7,7,7,7,7,7,7,7,7,          & 
     &          7,7,7,9,9,9,9,9,8,8,8,8,8,8,8,8,8,8,8,8,8,10,10,10,       & 
     &          10,10,9,9,9,9,9,9,9,9,9,9,9,9,9,11,11,11,11,11,10,        & 
     &          10,10,10,10,10,10,10,10,10,10,10,10,12,12,12,12,12,       & 
     &          11,11,11,11,11,11,11,11,11,11,11,11,11,13,13,12,12,       & 
     &          12,12,12,12,12/ 
        DATA GNI/4.0  ,6.0  ,4.0  ,4.0  ,2.0  ,                           & 
     &          2.0  ,4.0  ,6.0  ,2.0  ,4.0  ,                            & 
     &          6.0  ,4.0  ,2.0  ,2.0  ,2.0  ,                            & 
     &          4.0  ,6.0  ,8.0  ,2.0  ,4.0  ,                            & 
     &          6.0  ,4.0  ,4.0  ,6.0  ,2.0  ,                            & 
     &          4.0  ,6.0  ,4.0  ,2.0  ,4.0  ,                            & 
     &          6.0  ,2.0  ,4.0  ,4.0  ,2.0  ,                            & 
     &          4.0  ,6.0  ,8.0  ,10.0  ,6.0  ,                           & 
     &          8.0  ,2.0  ,4.0  ,6.0  ,2.0  ,                            & 
     &          4.0  ,6.0  ,8.0  ,4.0  ,6.0  ,                            & 
     &          2.0  ,2.0  ,4.0  ,6.0  ,8.0  ,                            & 
     &          2.0  ,4.0  ,6.0  ,4.0  ,2.0  ,                            & 
     &          4.0  ,6.0  ,2.0  ,4.0  ,4.0  ,                            & 
     &          6.0  ,8.0  ,10.0  ,2.0  ,4.0  ,                           & 
     &          6.0  ,8.0  ,4.0  ,2.0  ,6.0  ,                            & 
     &          8.0  ,2.0  ,4.0  ,6.0  ,4.0  ,                            & 
     &          6.0  ,4.0  ,6.0  ,2.0  ,4.0  ,                            & 
     &          2.0  ,4.0  ,6.0  ,2.0  ,4.0  ,                            & 
     &          4.0  ,6.0  ,8.0  ,10.0  ,4.0  ,                           & 
     &          2.0  ,6.0  ,8.0  ,2.0  ,4.0  ,                            & 
     &          6.0  ,8.0  ,2.0  ,4.0  ,6.0  ,                            & 
     &          4.0  ,6.0  ,2.0  ,4.0  ,6.0  ,                            & 
     &          2.0  ,4.0  ,4.0  ,6.0  ,8.0  ,                            & 
     &          10.0  ,2.0  ,4.0  ,6.0  ,8.0  ,                           & 
     &          4.0  ,2.0  ,6.0  ,8.0  ,4.0  ,                            & 
     &          6.0  ,2.0  ,4.0  ,6.0  ,2.0  ,                            & 
     &          4.0  ,6.0  ,2.0  ,4.0  ,2.0  ,                            & 
     &          4.0  ,6.0  ,8.0  ,6.0  ,8.0  ,                            & 
     &          4.0  ,2.0  ,4.0  ,6.0  ,2.0  ,                            & 
     &          4.0  ,6.0  ,2.0  ,4.0  ,2.0  ,                            & 
     &          4.0  ,6.0  ,2.0  ,4.0  ,6.0  ,                            & 
     &          8.0  ,4.0  ,2.0  ,6.0  ,8.0  ,                            & 
     &          4.0  ,6.0  ,2.0  ,4.0  ,6.0  ,                            & 
     &          2.0  ,4.0  ,2.0  ,4.0  ,6.0  ,                            & 
     &          2.0  ,4.0  ,6.0  ,8.0  ,4.0  ,                            & 
     &          2.0  ,6.0  ,8.0  ,4.0  ,6.0  ,                            & 
     &          2.0  ,4.0  ,6.0  ,2.0  ,4.0  ,                            & 
     &          2.0  ,4.0  ,6.0  ,4.0  ,2.0  ,                            & 
     &          6.0  ,8.0  ,2.0  ,4.0  ,6.0  ,                            & 
     &          8.0  ,4.0  ,6.0  ,2.0  ,4.0  ,                            & 
     &          6.0  ,2.0  ,4.0  ,2.0  ,4.0  ,                            & 
     &          6.0  ,4.0  ,2.0  ,6.0  ,8.0  ,                            & 
     &          2.0  ,4.0  ,6.0  ,8.0  ,4.0  ,                            & 
     &          6.0  ,2.0  ,4.0  ,6.0  ,2.0  ,                            & 
     &          4.0  ,4.0  ,2.0  ,2.0  ,4.0  ,                            & 
     &          6.0  ,4.0  ,6.0  / 
        DATA ENNI/0.0  ,2.383371,2.384363,3.575739,3.575739,              & 
     &          10.32619,10.33038,10.33617,10.67904,10.69042,             & 
     &          10.92429,10.92973,10.93217,11.60284,11.75037,             & 
     &          11.75317,11.75780,11.76412,11.83769,11.83997,             & 
     &          11.84472,11.99580,12.00032,12.00975,12.12207,             & 
     &          12.12649,12.35701,12.35614,12.84713,12.85333,             & 
     &          12.86185,12.91211,12.92268,12.97078,12.97568,             & 
     &          12.97693,12.97929,12.98350,12.98958,12.99502,             & 
     &          13.00392,13.00161,13.00483,13.00074,13.01686,             & 
     &          13.01822,13.01983,13.02095,13.03344,13.03636,             & 
     &          13.20179,13.23674,13.23917,13.24364,13.25041,             & 
     &          13.26429,13.26623,13.27127,13.32189,13.61527,             & 
     &          13.62076,13.62945,13.64202,13.65185,13.66270,             & 
     &          13.66493,13.66914,13.67609,13.66580,13.67249,             & 
     &          13.67410,13.68043,13.66588,13.66872,13.67695,             & 
     &          13.68464,13.67869,13.68191,13.68836,13.69398,             & 
     &          13.69673,13.70310,13.70607,13.92292,13.92614,             & 
     &          13.95653,13.96207,13.97100,13.97749,13.98841,             & 
     &          13.97948,13.98097,13.98543,13.99324,13.98568,             & 
     &          13.98754,13.98803,13.99674,13.98865,13.98865,             & 
     &          13.98865,13.99696,13.99237,13.99473,13.99944,             & 
     &          14.00155,14.00384,14.13620,14.14326,14.15244,             & 
     &          14.15045,14.15455,14.15417,14.15417,14.15417,             & 
     &          14.15417,14.15690,14.15690,14.15690,14.16508,             & 
     &          14.15827,14.16025,14.15864,14.16843,14.16313,             & 
     &          14.17035,14.16645,14.16645,14.16831,14.23464,             & 
     &          14.24468,14.25113,14.25212,14.25212,14.25683,             & 
     &          14.25683,14.25683,14.25683,14.25882,14.25882,             & 
     &          14.26043,14.26043,14.26545,14.27073,14.27109,             & 
     &          14.27109,14.27109,14.36247,14.36247,14.31821,             & 
     &          14.31821,14.31821,14.32329,14.32329,14.32329,             & 
     &          14.32329,14.32403,14.32403,14.32465,14.32465,             & 
     &          14.33234,14.33544,14.33494,14.33494,14.33494,             & 
     &          14.36272,14.36272,14.36433,14.36433,14.36433,             & 
     &          14.36830,14.36830,14.36830,14.36830,14.36854,             & 
     &          14.36854,14.37016,14.37016,14.37896,14.38119,             & 
     &          14.38107,14.38107,14.38107,14.39557,14.39557,             & 
     &          14.39768,14.39768,14.39768,14.40152,14.40152,             & 
     &          14.40202,14.40202,14.40264,14.40264,14.40264,             & 
     &          14.40264,14.41206,14.41206,14.41442,14.41442,             & 
     &          14.41442,14.42012,14.42012,14.42099,14.42099,             & 
     &          14.42099,14.42583,14.42583,14.42682,14.42682,             & 
     &          14.42781,14.42781,14.42781,14.42781,14.43636,             & 
     &          14.43636,14.43698,14.43698,14.43698,14.46253,             & 
     &          14.44021,14.44455,14.44455,14.45434,14.45434,             & 
     &          14.45434,14.45980,14.45980/ 
        DATA NNII/2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,2,3,3,3,3,2,3,        & 
     &          3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,3,4,4,4,          & 
     &          4,4,4,4,4,4,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,          & 
     &          4,4,4,4,4,4,4,4,4,4,3,3,3,5,5,5,5,5,5,5,5,5,5,5,          & 
     &          5,5,5,5,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3/ 
        DATA GNII/1.0  ,3.0  ,5.0  ,5.0  ,1.0  ,                          & 
     &          5.0  ,7.0  ,5.0  ,3.0  ,5.0  ,                            & 
     &          3.0  ,1.0  ,5.0  ,1.0  ,3.0  ,                            & 
     &          5.0  ,3.0  ,3.0  ,3.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,3.0  ,3.0  ,1.0  ,                            & 
     &          3.0  ,5.0  ,5.0  ,1.0  ,5.0  ,                            & 
     &          7.0  ,9.0  ,5.0  ,3.0  ,5.0  ,                            & 
     &          7.0  ,5.0  ,3.0  ,1.0  ,7.0  ,                            & 
     &          3.0  ,1.0  ,3.0  ,5.0  ,3.0  ,                            & 
     &          3.0  ,3.0  ,5.0  ,7.0  ,1.0  ,                            & 
     &          3.0  ,5.0  ,3.0  ,5.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,1.0  ,5.0  ,7.0  ,                            & 
     &          9.0  ,5.0  ,3.0  ,5.0  ,7.0  ,                            & 
     &          5.0  ,3.0  ,1.0  ,7.0  ,5.0  ,                            & 
     &          7.0  ,9.0  ,7.0  ,7.0  ,9.0  ,                            & 
     &          11.0  ,3.0  ,9.0  ,7.0  ,5.0  ,                           & 
     &          3.0  ,5.0  ,1.0  ,3.0  ,5.0  ,                            & 
     &          1.0  ,3.0  ,5.0  ,3.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,5.0  ,7.0  ,9.0  ,                            & 
     &          7.0  ,7.0  ,9.0  ,11.0  ,9.0  ,                           & 
     &          1.0  ,3.0  ,5.0  ,7.0  ,9.0  ,                            & 
     &          3.0  ,5.0  ,7.0  ,5.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,9.0  ,11.0  ,7.0  ,                           & 
     &          5.0  ,3.0  ,1.0  ,3.0  ,5.0  ,                            & 
     &          7.0  ,9.0  / 
        DATA ENNII/0.0  ,6.0876831e-03,1.6279284e-02,1.898923,4.052723,   & 
     &          5.848106,11.43604,11.43781,11.43801,13.54146,             & 
     &          13.54146,13.54228,17.87734,18.46259,18.46651,             & 
     &          18.48341,18.49722,19.23384,20.40944,20.64636,             & 
     &          20.65389,20.66582,20.67651,20.94027,21.14861,             & 
     &          21.15298,21.16022,21.59986,22.10340,23.12481,             & 
     &          23.13218,23.14229,23.19670,23.23962,23.24260,             & 
     &          23.24636,23.41565,23.42207,23.42555,23.47490,             & 
     &          23.57225,24.36823,24.37465,24.38944,24.53166,             & 
     &          25.06612,25.13369,25.14001,25.15193,25.18946,             & 
     &          25.19245,25.20124,25.23510,25.46049,25.53877,             & 
     &          25.54572,25.55447,25.58160,25.99668,26.00464,             & 
     &          26.01527,26.02787,26.06667,26.06994,26.07548,             & 
     &          26.12440,26.13011,26.13327,26.16475,26.16510,             & 
     &          26.16800,26.16849,26.17391,26.19663,26.19758,             & 
     &          26.20937,26.20252,26.21087,26.21191,26.21252,             & 
     &          26.22134,26.22182,26.25393,26.25770,26.26368,             & 
     &          26.55921,26.56489,26.58065,26.63554,27.36569,             & 
     &          27.36569,27.36569,27.40948,27.40948,27.40999,             & 
     &          27.41783,27.42901,27.42963,27.43824,27.43947,             & 
     &          27.77609,27.77805,27.78169,27.78704,27.79372,             & 
     &          28.01910,28.02209,28.02755,28.54429,30.17253,             & 
     &          30.17448,30.17763,30.18179,30.18682,30.34387,             & 
     &          30.34864,30.35188,30.41607,30.41652,30.41750,             & 
     &          30.41894,30.42068/ 
        DATA NNIII/2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,3,3,3,3,3,3,3,3,       & 
     &          3,4,3,3,3,3,3,3,4,4,3,3,3,3,4,4,4,4,3,3,3,3,3,3,          & 
     &          3,3,3,3,3,5,3,3,3,3,3,3,3,5,5,3,3,5,5,5,5,6,6,6,          & 
     &          6,6,6,4,4,4,3,3,4,4,4,4,4,4,3,3,4,4,4,4,4,4,4,4,          & 
     &          4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,          & 
     &          4,4,4,4,4,4,4,3,3,5,5,5,5/ 
        DATA GNIII/2.0  ,4.0  ,2.0  ,4.0  ,6.0  ,                         & 
     &          6.0  ,4.0  ,2.0  ,2.0  ,4.0  ,                            & 
     &          4.0  ,6.0  ,4.0  ,2.0  ,2.0  ,                            & 
     &          4.0  ,2.0  ,4.0  ,4.0  ,6.0  ,                            & 
     &          2.0  ,4.0  ,6.0  ,2.0  ,4.0  ,                            & 
     &          2.0  ,2.0  ,4.0  ,2.0  ,4.0  ,                            & 
     &          6.0  ,8.0  ,2.0  ,4.0  ,4.0  ,                            & 
     &          2.0  ,4.0  ,6.0  ,4.0  ,6.0  ,                            & 
     &          6.0  ,8.0  ,4.0  ,6.0  ,2.0  ,                            & 
     &          4.0  ,6.0  ,8.0  ,10.0  ,2.0  ,                           & 
     &          4.0  ,6.0  ,8.0  ,2.0  ,4.0  ,                            & 
     &          6.0  ,6.0  ,4.0  ,2.0  ,6.0  ,                            & 
     &          8.0  ,4.0  ,6.0  ,4.0  ,2.0  ,                            & 
     &          6.0  ,8.0  ,8.0  ,10.0  ,4.0  ,                           & 
     &          6.0  ,6.0  ,8.0  ,8.0  ,10.0  ,                           & 
     &          2.0  ,4.0  ,6.0  ,4.0  ,6.0  ,                            & 
     &          2.0  ,4.0  ,2.0  ,4.0  ,6.0  ,                            & 
     &          8.0  ,2.0  ,4.0  ,4.0  ,6.0  ,                            & 
     &          4.0  ,2.0  ,4.0  ,6.0  ,4.0  ,                            & 
     &          6.0  ,8.0  ,10.0  ,4.0  ,6.0  ,                           & 
     &          2.0  ,4.0  ,6.0  ,8.0  ,6.0  ,                            & 
     &          4.0  ,2.0  ,6.0  ,8.0  ,4.0  ,                            & 
     &          6.0  ,8.0  ,10.0  ,6.0  ,8.0  ,                           & 
     &          6.0  ,8.0  ,10.0  ,12.0  ,8.0  ,                          & 
     &          10.0  ,8.0  ,6.0  ,4.0  ,2.0  ,                           & 
     &          6.0  ,4.0  ,4.0  ,6.0  ,2.0  ,                            & 
     &          4.0  ,6.0  ,8.0  / 
        DATA ENNIII/0.0  ,2.1635452e-02,7.180255,7.098413,7.108480,       & 
     &          12.52548,12.52643,16.24252,18.08651,18.10019,             & 
     &          23.16076,25.17799,25.18006,27.43827,28.56680,             & 
     &          28.56730,30.45896,30.46342,33.13367,33.13441,             & 
     &          35.65022,35.65797,35.67233,36.84229,36.85629,             & 
     &          38.44641,38.32793,38.33453,38.39367,38.39807,             & 
     &          38.40689,38.41771,38.64517,38.64825,38.95919,             & 
     &          39.34056,39.34595,39.35325,39.39646,39.40031,             & 
     &          39.71098,39.71098,39.79651,39.80747,40.55027,             & 
     &          40.94474,40.94909,40.95552,40.96437,41.26192,             & 
     &          41.26358,41.26631,41.26982,41.37555,41.47835,             & 
     &          41.48166,41.68555,41.69232,41.69667,42.12335,             & 
     &          42.13715,42.39634,42.39655,42.48893,42.49769,             & 
     &          42.49625,42.49625,42.54757,42.54757,43.95493,             & 
     &          43.95493,44.00932,44.00932,44.04135,44.04135,             & 
     &          45.69180,45.69957,45.71402,46.28896,46.29317,             & 
     &          46.46321,46.47039,46.71232,46.71811,46.72555,             & 
     &          46.73671,46.81577,46.81788,46.85206,46.86286,             & 
     &          46.92110,47.02857,47.03412,47.04068,47.61238,             & 
     &          47.61238,47.61845,47.62763,47.75000,47.75000,             & 
     &          47.77108,47.77108,49.01428,47.77802,47.88887,             & 
     &          47.88887,47.88887,47.97657,47.97913,47.98245,             & 
     &          47.98245,47.98363,47.98760,48.07270,48.08297,             & 
     &          48.11119,48.11662,48.12305,48.13089,48.12993,             & 
     &          48.14229,48.14024,48.14488,48.15087,48.15427,             & 
     &          48.15307,48.16119,49.16950,49.17073,50.71214,             & 
     &          50.71214,50.71214,50.71214/ 
        DATA NNIV/2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,        & 
     &          3,3,3,3,3,3,4,3,3,3,3,3,4,4,4,3,3,3,3,4,4,4,4,3,          & 
     &          3,3,4,4,4,4,3,4,5,5,5,5,5,5,5,6,6,6,4,4,4,4,5,5,4/ 
        DATA GNIV/1.0  ,1.0  ,3.0  ,7.0  ,3.0  ,                          & 
     &          1.0  ,3.0  ,5.0  ,5.0  ,1.0  ,                            & 
     &          3.0  ,1.0  ,1.0  ,3.0  ,5.0  ,                            & 
     &          3.0  ,5.0  ,7.0  ,5.0  ,1.0  ,                            & 
     &          3.0  ,5.0  ,3.0  ,3.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,3.0  ,1.0  ,3.0  ,                            & 
     &          5.0  ,5.0  ,5.0  ,5.0  ,7.0  ,                            & 
     &          9.0  ,1.0  ,3.0  ,5.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,7.0  ,3.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,5.0  ,3.0  ,1.0  ,                            & 
     &          5.0  ,5.0  ,7.0  ,9.0  ,3.0  ,                            & 
     &          7.0  ,3.0  ,3.0  ,5.0  ,7.0  ,                            & 
     &          7.0  ,9.0  ,11.0  ,3.0  ,5.0  ,                           & 
     &          7.0  ,5.0  ,3.0  ,5.0  ,7.0  ,                            & 
     &          3.0  ,5.0  ,7.0  / 
        DATA ENNIV/0.0  ,8.323934,8.331770,8.349648,16.20427,             & 
     &          21.75491,21.76399,21.77946,23.41898,29.18244,             & 
     &          46.76804,50.15470,50.32483,50.32679,50.33118,             & 
     &          52.06988,52.07031,52.07132,53.20933,57.68086,             & 
     &          57.69048,57.71067,58.64906,59.62210,60.05779,             & 
     &          60.05779,60.07403,60.44809,61.27855,61.27855,             & 
     &          61.29070,61.78379,61.95650,61.97423,61.97423,             & 
     &          61.97423,62.44215,62.44215,62.44215,62.67301,             & 
     &          62.67685,62.68218,62.77282,62.86333,63.40415,             & 
     &          63.40415,63.40415,63.41109,63.41767,63.41767,             & 
     &          63.80760,64.05482,64.05569,64.05706,64.39976,             & 
     &          64.70402,68.21900,68.53058,68.53058,68.53058,             & 
     &          68.73986,68.73986,68.73986,71.28416,71.28416,             & 
     &          71.28416,73.28070,73.60580,73.60580,73.61063,             & 
     &          78.63129,78.63129,78.63129/ 
        DATA NNV/2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,6,         & 
     &          6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8/ 
        DATA GNV/2.0  ,2.0  ,4.0  ,2.0  ,2.0  ,                           & 
     &          4.0  ,4.0  ,6.0  ,2.0  ,2.0  ,                            & 
     &          4.0  ,4.0  ,6.0  ,2.0  ,2.0  ,                            & 
     &          4.0  ,4.0  ,6.0  ,2.0  ,2.0  ,                            & 
     &          4.0  ,4.0  ,6.0  ,6.0  ,8.0  ,                            & 
     &          8.0  ,10.0  ,12.0  ,2.0  ,2.0  ,                          & 
     &          4.0  ,4.0  ,6.0  ,6.0  ,8.0  ,                            & 
     &          8.0  ,10.0  ,12.0  ,14.0  ,2.0  ,                         & 
     &          2.0  ,4.0  ,4.0  ,6.0  ,6.0  ,                            & 
     &          8.0  ,8.0  ,10.0  ,12.0  ,14.0  ,16.0  / 
        DATA ENNV/0.0  ,9.976473,10.00851,56.55396,59.23740,              & 
     &          59.24660,60.05890,60.06188,75.17694,76.26962,             & 
     &          76.26962,76.61120,76.61120,83.55153,84.09893,             & 
     &          84.09893,84.27598,84.27598,88.02306,88.33514,             & 
     &          88.33514,88.43854,88.43742,88.44214,88.44214,             & 
     &          88.44313,88.44313,88.44313,90.68689,90.88043,             & 
     &          90.88043,90.94527,90.94527,90.94912,90.94912,             & 
     &          90.94974,90.94974,90.94974,90.94974,92.40136,             & 
     &          92.53167,92.53167,92.57358,92.57358,92.57618,             & 
     &          92.57618,92.57668,92.57668,92.57668,92.57668,             & 
     &          92.57668/ 
        DATA NNVI/1,2,2,2,2,2,3,4/ 
        DATA GNVI/1.0  ,3.0  ,1.0  ,3.0  ,5.0  ,3.0  ,3.0  ,3.0  / 
        DATA ENNVI/0.0  ,419.8009,426.2953,426.2965,426.3325,             & 
     &          425.7398,497.9737,521.5830/ 
        DATA NOI/2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,3,3,3,3,3,3,3,3,3,         & 
     &          4,4,4,4,4,4,3,3,3,5,5,3,4,4,4,4,4,4,4,4,5,5,5,6,          & 
     &          6,5,5,5,5,5,5,5,5,6,6,6,7,7,6,6,6,6,6,6,6,6,8,8,          & 
     &          7,7,7,7,7,7,7,7,9,9,8,8,8,8,8,8,8,8,10,10,9,9,9,9,        & 
     &          9,9,9,9,11,11,10,10,10,10,10,10,10,10,3,3,3,3,3,3,        & 
     &          3,3,3,3,3,3,4,3,3,3,3,3,3,3,3,3,3,3,4,4,4,2,2,2,3,        & 
     &          3,3,3,3,5,4,4,4,4,4,4,4,4,4,4,4,3,6,5,5,5,5,5,5,5,        & 
     &          5,5,5,7,6,6,6,2/ 
        DATA GOI/5.0  ,3.0  ,1.0  ,5.0  ,1.0  ,                           & 
     &          5.0  ,3.0  ,3.0  ,5.0  ,7.0  ,                            & 
     &          5.0  ,3.0  ,1.0  ,5.0  ,3.0  ,                            & 
     &          9.0  ,7.0  ,5.0  ,5.0  ,3.0  ,                            & 
     &          1.0  ,7.0  ,5.0  ,3.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,5.0  ,3.0  ,1.0  ,                            & 
     &          7.0  ,5.0  ,3.0  ,5.0  ,3.0  ,                            & 
     &          5.0  ,9.0  ,7.0  ,5.0  ,3.0  ,                            & 
     &          1.0  ,7.0  ,5.0  ,3.0  ,5.0  ,                            & 
     &          3.0  ,1.0  ,5.0  ,3.0  ,9.0  ,                            & 
     &          7.0  ,5.0  ,3.0  ,1.0  ,7.0  ,                            & 
     &          5.0  ,3.0  ,5.0  ,3.0  ,1.0  ,                            & 
     &          5.0  ,3.0  ,9.0  ,7.0  ,5.0  ,                            & 
     &          3.0  ,1.0  ,7.0  ,5.0  ,3.0  ,                            & 
     &          5.0  ,3.0  ,9.0  ,7.0  ,5.0  ,                            & 
     &          3.0  ,1.0  ,7.0  ,5.0  ,3.0  ,                            & 
     &          5.0  ,3.0  ,9.0  ,7.0  ,5.0  ,                            & 
     &          3.0  ,1.0  ,7.0  ,5.0  ,3.0  ,                            & 
     &          5.0  ,3.0  ,9.0  ,7.0  ,5.0  ,                            & 
     &          3.0  ,1.0  ,7.0  ,5.0  ,3.0  ,                            & 
     &          5.0  ,3.0  ,9.0  ,7.0  ,5.0  ,                            & 
     &          3.0  ,1.0  ,7.0  ,5.0  ,3.0  ,                            & 
     &          7.0  ,5.0  ,3.0  ,9.0  ,7.0  ,                            & 
     &          5.0  ,5.0  ,3.0  ,1.0  ,7.0  ,                            & 
     &          3.0  ,5.0  ,5.0  ,5.0  ,3.0  ,                            & 
     &          1.0  ,9.0  ,7.0  ,5.0  ,9.0  ,                            & 
     &          11.0  ,9.0  ,7.0  ,7.0  ,7.0  ,                           & 
     &          5.0  ,3.0  ,5.0  ,3.0  ,1.0  ,                            & 
     &          7.0  ,5.0  ,3.0  ,3.0  ,5.0  ,                            & 
     &          5.0  ,9.0  ,7.0  ,5.0  ,9.0  ,                            & 
     &          11.0  ,9.0  ,7.0  ,7.0  ,5.0  ,                           & 
     &          3.0  ,1.0  ,1.0  ,5.0  ,9.0  ,                            & 
     &          7.0  ,5.0  ,9.0  ,11.0  ,9.0  ,                           & 
     &          7.0  ,5.0  ,3.0  ,1.0  ,5.0  ,                            & 
     &          5.0  ,3.0  ,1.0  ,3.0  / 
        DATA ENOI/0.0  ,01.9651687e-02,2.8082693e-02,1.967363,0.4206081,  & 
     &          9.146132,9.521420,10.74028,10.74053,10.74098,             & 
     &          10.98893,10.98886,10.98895,11.83768,11.93056,             & 
     &          12.07869,12.07870,12.07870,12.07872,12.07872,             & 
     &          12.07872,12.08711,12.08711,12.08711,12.28604,             & 
     &          12.28612,12.28627,12.35891,12.35891,12.35891,             & 
     &          12.53927,12.54078,12.54176,12.66092,12.69755,             & 
     &          12.72854,12.75377,12.75377,12.75377,12.75377,             & 
     &          12.75377,12.75911,12.75911,12.75911,12.87829,             & 
     &          12.87829,12.87829,13.02082,13.03891,13.06624,             & 
     &          13.06624,13.06624,13.06624,13.06624,13.06913,             & 
     &          13.06913,13.06913,13.13145,13.13145,13.13145,             & 
     &          13.21004,13.22030,13.23559,13.23559,13.23559,             & 
     &          13.23559,13.23559,13.23740,13.23740,13.23740,             & 
     &          13.32166,13.32807,13.33749,13.33749,13.33749,             & 
     &          13.33749,13.33749,13.33869,13.33869,13.33869,             & 
     &          13.39308,13.39756,13.40353,385.3597,13.40353,             & 
     &          13.40353,13.40353,13.40488,13.40488,13.40488,             & 
     &          13.44262,13.44449,13.44872,13.44872,13.44872,             & 
     &          13.44872,13.44872,13.44966,13.44966,13.44966,             & 
     &          13.47577,13.47812,13.48112,13.48112,13.48112,             & 
     &          13.48112,13.48112,13.48148,13.48148,13.48148,             & 
     &          14.04685,14.04687,14.04730,14.09888,14.09975,             & 
     &          14.10046,14.12320,14.12450,14.12526,14.13382,             & 
     &          14.37218,14.46048,15.22525,15.28698,15.29424,             & 
     &          15.29817,15.40062,15.40062,15.40062,15.40372,             & 
     &          15.40390,15.40622,15.40550,15.41465,15.59420,             & 
     &          15.59514,15.59577,15.65520,15.66431,15.66970,             & 
     &          15.78109,15.78181,15.78222,15.82895,15.94391,             & 
     &          16.01073,16.07676,16.07676,16.07676,16.07836,             & 
     &          16.07844,16.08080,16.08005,16.08545,16.11433,             & 
     &          16.11550,16.11614,16.23505,16.35702,16.35702,             & 
     &          16.35702,16.35702,16.39057,16.39063,16.39308,             & 
     &          16.39308,16.40451,16.40451,16.40451,16.54127,             & 
     &          16.56668,16.56668,16.56668,23.53702/ 
        DATA NOII/2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,2,3,3,3,3,3,3,3,3,        & 
     &          3,3,3,3,3,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,          & 
     &          3,3,3,3,3,3,3,3,3,4,4,4,4,4,3,4,4,4,4,4,4,4,4,3,          & 
     &          3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,          & 
     &          4,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,          & 
     &          5,5,5,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,          & 
     &          5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,3,3,3,4,4,4,4,          & 
     &          4,4,4,4,4,4,3,3,4,4,4,4,4,4,4,5,5,3,3,3,3,3,4/ 
        DATA GOII/4.0  ,6.0  ,4.0  ,4.0  ,2.0  ,                          & 
     &          6.0  ,4.0  ,2.0  ,6.0  ,4.0  ,                            & 
     &          2.0  ,4.0  ,6.0  ,2.0  ,4.0  ,                            & 
     &          2.0  ,2.0  ,2.0  ,4.0  ,6.0  ,                            & 
     &          8.0  ,6.0  ,4.0  ,2.0  ,4.0  ,                            & 
     &          6.0  ,4.0  ,6.0  ,4.0  ,4.0  ,                            & 
     &          2.0  ,2.0  ,4.0  ,2.0  ,6.0  ,                            & 
     &          8.0  ,6.0  ,4.0  ,4.0  ,6.0  ,                            & 
     &          8.0  ,10.0  ,6.0  ,4.0  ,2.0  ,                           & 
     &          2.0  ,4.0  ,2.0  ,4.0  ,6.0  ,                            & 
     &          8.0  ,6.0  ,8.0  ,4.0  ,2.0  ,                            & 
     &          4.0  ,6.0  ,2.0  ,4.0  ,6.0  ,                            & 
     &          2.0  ,4.0  ,6.0  ,2.0  ,4.0  ,                            & 
     &          6.0  ,8.0  ,4.0  ,6.0  ,2.0  ,                            & 
     &          4.0  ,2.0  ,4.0  ,8.0  ,6.0  ,                            & 
     &          10.0  ,8.0  ,4.0  ,6.0  ,2.0  ,                           & 
     &          4.0  ,4.0  ,6.0  ,8.0  ,10.0  ,                           & 
     &          2.0  ,4.0  ,6.0  ,8.0  ,4.0  ,                            & 
     &          6.0  ,4.0  ,2.0  ,4.0  ,2.0  ,                            & 
     &          6.0  ,8.0  ,2.0  ,6.0  ,4.0  ,                            & 
     &          8.0  ,5.80  ,4.0  ,2.0  ,6.0  ,                           & 
     &          8.0  ,10.0  ,12.0  ,8.0  ,10.0  ,                         & 
     &          4.0  ,6.0  ,4.0  ,6.0  ,8.0  ,                            & 
     &          10.0  ,6.0  ,8.0  ,2.0  ,4.0  ,                           & 
     &          6.0  ,2.0  ,4.0  ,6.0  ,4.0  ,                            & 
     &          2.0  ,4.0  ,6.0  ,8.0  ,2.0  ,                            & 
     &          4.0  ,6.0  ,4.0  ,6.0  ,2.0  ,                            & 
     &          4.0  ,6.0  ,8.0  ,6.0  ,4.0  ,                            & 
     &          2.0  ,6.0  ,8.0  ,8.0  ,6.0  ,                            & 
     &          4.0  ,2.0  ,6.0  ,8.0  ,10.0  ,                           & 
     &          12.0  ,8.0  ,10.0  ,4.0  ,6.0  ,                          & 
     &          4.0  ,6.0  ,8.0  ,10.0  ,6.0  ,                           & 
     &          8.0  ,4.0  ,6.0  ,8.0  ,6.0  ,                            & 
     &          8.0  ,4.0  ,6.0  ,2.0  ,4.0  ,                            & 
     &          8.0  ,10.0  ,6.0  ,8.0  ,6.0  ,                           & 
     &          4.0  ,2.0  ,4.0  ,6.0  ,10.0  ,                           & 
     &          12.0  ,2.0  ,4.0  ,4.0  ,6.20  ,                          & 
     &          10.0  ,8.0  ,6.0  ,4.0  ,2.0  ,                           & 
     &          6.0  / 
        DATA ENOII/0.0  ,3.323850,3.326454,5.017305,5.017491,             & 
     &          14.85813,14.87838,14.88860,20.58005,20.57736,             & 
     &          22.96648,22.97954,23.001876,23.41940,23.44172,            & 
     &          24.26523,25.28586,25.63160,25.63849,25.64984,             & 
     &          25.66529,25.66142,25.66154,25.83188,25.83760,             & 
     &          25.84900,26.22564,26.24928,26.30498,26.35845,             & 
     &          26.37943,26.55392,26.56133,28.12621,28.35835,             & 
     &          28.36128,28.51330,28.51270,28.67733,28.68403,             & 
     &          28.69369,28.70637,28.82200,28.83108,28.83932,             & 
     &          28.82414,28.82992,28.85285,28.85711,28.85729,             & 
     &          28.85808,28.86334,28.88355,28.94193,28.95606,             & 
     &          29.06249,29.06893,29.58618,29.59923,29.61924,             & 
     &          29.79726,29.82051,30.42546,30.47162,30.47763,             & 
     &          30.48836,30.50400,30.74951,30.77135,30.80112,             & 
     &          30.81214,31.02747,31.02747,31.14773,31.14812,             & 
     &          31.31967,31.31982,31.37404,31.37430,31.46620,             & 
     &          31.46649,31.55199,31.55199,31.55199,31.56553,             & 
     &          31.61407,31.61407,31.61407,31.61407,31.61407,             & 
     &          31.62925,31.63375,31.63644,31.63766,31.65117,             & 
     &          31.65364,31.67396,31.69345,31.70178,31.71699,             & 
     &          31.70200,31.71709,31.72948,31.72935,31.70999,             & 
     &          31.71043,31.71889,31.73747,31.71911,31.73823,             & 
     &          31.72081,31.72752,31.75062,31.75112,31.75553,             & 
     &          31.75715,31.75586,31.75803,31.95026,31.96318,             & 
     &          31.98375,32.03889,32.06284,32.14771,32.14780,             & 
     &          32.35511,32.35511,32.36540,32.38251,32.39264,             & 
     &          32.39264,32.40412,32.44667,32.46798,32.88345,             & 
     &          32.88345,32.88345,32.88345,32.90963,32.91418,             & 
     &          32.91418,32.92780,32.92780,32.93536,32.94354,             & 
     &          32.95061,32.96264,32.93858,32.94181,32.95049,             & 
     &          32.97082,32.95073,32.97146,32.96227,32.96227,             & 
     &          32.97119,32.97528,32.97826,32.97999,32.97863,             & 
     &          32.97999,33.19875,33.19968,33.20123,34.06365,             & 
     &          34.06901,34.08607,34.08607,34.17174,34.17174,             & 
     &          34.20029,34.20029,34.20504,34.20504,34.21390,             & 
     &          34.21390,34.21960,34.22819,34.22819,34.23350,             & 
     &          34.23350,34.25269,34.25269,34.48530,34.48530,             & 
     &          36.19083,36.18759,36.19109,36.19123,36.19131,37.05294/ 
        DATA NOIII/2,2,2,2,2,2,2,3,3,2,2,2,2,2,2,3,3,3,3,2,2,2,3,3,       & 
     &          3,3,3,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,          & 
     &          2,3,3,3,4,4,4,4,3,3,3,3,3,3,4,4,4,4,4,3,3,3,4,4,          & 
     &          4,4,4,3,3,3,3,4,4,4,4,3,3,3,4,4,4,4,4,4,4,4,5,5,          & 
     &          5,5,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,          & 
     &          5,5,5,5,5,5,5,5,5,3,3,3,6,6,6,6,7,3,3,4,4,4,3,4,          & 
     &          4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,5/ 
        DATA GOIII/1.0  ,3.0  ,5.0  ,5.0  ,1.0  ,                         & 
     &          5.0  ,7.0  ,5.0  ,3.0  ,5.0  ,                            & 
     &          3.0  ,1.0  ,5.0  ,3.0  ,3.0  ,                            & 
     &          1.0  ,3.0  ,5.0  ,3.0  ,5.0  ,                            & 
     &          3.0  ,1.0  ,3.0  ,3.0  ,5.0  ,                            & 
     &          7.0  ,3.0  ,5.0  ,1.0  ,3.0  ,                            & 
     &          5.0  ,5.0  ,1.0  ,5.0  ,7.0  ,                            & 
     &          9.0  ,5.0  ,3.0  ,5.0  ,7.0  ,                            & 
     &          5.0  ,3.0  ,1.0  ,7.0  ,3.0  ,                            & 
     &          3.0  ,5.0  ,7.0  ,1.0  ,1.0  ,                            & 
     &          3.0  ,5.0  ,1.0  ,3.0  ,5.0  ,                            & 
     &          3.0  ,3.0  ,1.0  ,3.0  ,5.0  ,                            & 
     &          7.0  ,9.0  ,3.0  ,3.0  ,5.0  ,                            & 
     &          7.0  ,3.0  ,3.0  ,5.0  ,7.0  ,                            & 
     &          1.0  ,3.0  ,5.0  ,5.0  ,1.0  ,                            & 
     &          3.0  ,5.0  ,7.0  ,5.0  ,5.0  ,                            & 
     &          7.0  ,9.0  ,5.0  ,5.0  ,3.0  ,                            & 
     &          1.0  ,3.0  ,5.0  ,7.0  ,5.0  ,                            & 
     &          3.0  ,1.0  ,7.0  ,3.0  ,1.0  ,                            & 
     &          3.0  ,5.0  ,3.0  ,3.0  ,5.0  ,                            & 
     &          7.0  ,3.0  ,5.0  ,7.0  ,9.0  ,                            & 
     &          11.0  ,1.0  ,3.0  ,5.0  ,7.0  ,                           & 
     &          9.0  ,7.0  ,5.0  ,3.0  ,5.0  ,                            & 
     &          3.0  ,1.0  ,5.0  ,7.0  ,9.0  ,                            & 
     &          5.0  ,7.0  ,9.0  ,5.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,7.0  ,3.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,5.0  ,3.0  ,5.0  ,                            & 
     &          7.0  ,7.0  ,7.0  ,5.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,3.0  ,3.0  ,1.0  ,                            & 
     &          3.0  ,5.0  ,7.0  ,9.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,3.0  ,5.0  ,7.0  ,                            & 
     &          7.0  ,5.0  ,3.0  ,5.0  ,7.0  ,                            & 
     &          9.0  ,3.0  ,5.0  ,7.0  ,1.0  ,                            & 
     &          3.0  ,5.0  ,3.0  / 
        DATA ENOIII/0.0  ,1.4059945e-02,3.8038719e-02,2.513308,5.354124,  & 
     &          7.477820,14.88140,14.88477,14.88550,17.65325,             & 
     &          17.65339,17.65514,23.19140,24.43587,26.09378,             & 
     &          33.13600,33.15068,33.18253,33.85794,35.18196,             & 
     &          35.20895,35.22094,36.07438,36.43500,36.45190,             & 
     &          36.47919,36.89279,36.98353,37.22392,37.23410,             & 
     &          37.25028,38.01204,38.90675,40.22861,40.25288,             & 
     &          40.27497,40.26230,40.57149,40.57759,40.58673,             & 
     &          40.84922,40.86335,40.87098,41.14086,41.25951,             & 
     &          41.97723,41.99266,42.14902,42.56451,43.39812,             & 
     &          43.41013,43.43237,44.22956,44.24270,44.27655,             & 
     &          44.46952,45.03978,45.31862,45.32294,45.33144,             & 
     &          45.34384,45.35962,45.34443,45.43903,45.45230,             & 
     &          45.47797,45.62070,45.69189,45.69899,45.71153,             & 
     &          45.91510,45.92614,45.93959,45.98626,46.25228,             & 
     &          46.44183,45.21283,46.46955,46.62690,46.78899,             & 
     &          46.78899,46.78899,46.82767,46.91713,46.91867,             & 
     &          46.92080,47.01923,47.02679,47.03461,47.20199,             & 
     &          47.20199,47.20199,47.21141,47.24910,48.62968,             & 
     &          48.62968,48.62968,48.69874,48.86141,48.86587,             & 
     &          48.87442,48.91428,48.91908,48.92621,48.93560,             & 
     &          48.94701,49.36293,49.36248,49.36198,49.36323,             & 
     &          49.37332,49.40500,49.41368,49.41845,49.63815,             & 
     &          49.65178,49.65844,49.76514,49.77709,49.79367,             & 
     &          49.78386,49.78386,49.78386,49.81572,49.78386,             & 
     &          49.78386,49.78386,50.01249,50.03133,50.31391,             & 
     &          50.31750,50.32357,51.41365,51.47638,51.47638,             & 
     &          51.47638,52.44297,52.69355,52.85969,53.12613,             & 
     &          53.14089,53.16110,53.31682,54.18348,54.33549,             & 
     &          54.33549,54.34320,54.35460,54.36977,54.46407,             & 
     &          54.47044,54.48261,54.88958,54.88958,54.88958,             & 
     &          55.81414,55.82281,55.82951,56.14741,56.14741,             & 
     &          56.14741,56.31095,56.31095,56.31095,56.73994,             & 
     &          56.73994,56.73994,58.73808/ 
        DATA NOIV/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,        & 
     &          3,3,3,3,3,3,3,3,3,3,3,4,3,3,3,3,3,3,3,3,3,3,3,3,          & 
     &          3,3,4,4,3,3,3,3,3,3,5,3,3,3,3,5,5,5,5,3,4,4,4,3,          & 
     &          3,4,4,6,6,4,4,3,3,3,3,3,3,3,4,4,7,7,4,4,4,4,4,4,          & 
     &          4,4,4,4,4,4,4,4,4,4,3,8,8,4,4,3,3,3,3,3,3,3,3,3,          & 
     &          3,3,3,3,3,5,5,3,3,5,5,3,3,5,5,5,5,5,5,5,5,5,5,5,          & 
     &          3,3,3,3,3,3,3,3,3,6,6,6,6,4,4,3,4,4,7,7,7,7/ 
        DATA GOIV/2.0  ,4.0  ,2.0  ,4.0  ,6.0  ,                          & 
     &          6.0  ,4.0  ,2.0  ,2.0  ,6.0  ,                            & 
     &          4.0  ,6.0  ,4.0  ,2.0  ,4.0  ,                            & 
     &          2.0  ,2.0  ,4.0  ,2.0  ,4.0  ,                            & 
     &          6.0  ,2.0  ,4.0  ,2.0  ,4.0  ,                            & 
     &          2.0  ,4.0  ,6.0  ,8.0  ,4.0  ,                            & 
     &          2.0  ,4.0  ,6.0  ,4.0  ,6.0  ,                            & 
     &          2.0  ,2.0  ,4.0  ,6.0  ,8.0  ,                            & 
     &          10.0  ,2.0  ,4.0  ,6.0  ,8.0  ,                           & 
     &          4.0  ,6.0  ,6.0  ,4.0  ,2.0  ,                            & 
     &          4.0  ,6.0  ,6.0  ,8.0  ,4.0  ,                            & 
     &          2.0  ,2.0  ,4.0  ,2.0  ,4.0  ,                            & 
     &          6.0  ,2.0  ,4.0  ,4.0  ,6.0  ,                            & 
     &          6.0  ,8.0  ,2.0  ,2.0  ,4.0  ,                            & 
     &          6.0  ,6.0  ,8.0  ,2.0  ,4.0  ,                            & 
     &          4.0  ,6.0  ,2.0  ,4.0  ,4.0  ,                            & 
     &          6.0  ,2.0  ,4.0  ,6.0  ,2.0  ,                            & 
     &          4.0  ,4.0  ,6.0  ,6.0  ,8.0  ,                            & 
     &          2.0  ,2.0  ,4.0  ,6.0  ,8.0  ,                            & 
     &          6.0  ,4.0  ,2.0  ,4.0  ,6.0  ,                            & 
     &          6.0  ,8.0  ,4.0  ,6.0  ,6.0  ,                            & 
     &          8.0  ,2.0  ,6.0  ,8.0  ,4.0  ,                            & 
     &          2.0  ,4.0  ,6.0  ,2.0  ,4.0  ,                            & 
     &          6.0  ,8.0  ,2.0  ,4.0  ,6.0  ,                            & 
     &          6.0  ,4.0  ,4.0  ,6.0  ,8.0  ,                            & 
     &          2.0  ,4.0  ,6.0  ,8.0  ,4.0  ,                            & 
     &          6.0  ,6.0  ,4.0  ,2.0  ,4.0  ,                            & 
     &          6.0  ,8.0  ,6.0  ,4.0  ,2.0  ,                            & 
     &          6.0  ,8.0  ,4.0  ,2.0  ,6.0  ,                            & 
     &          4.0  ,2.0  ,4.0  ,6.0  ,6.0  ,                            & 
     &          8.0  ,4.0  ,2.0  ,2.0  ,4.0  ,                            & 
     &          6.0  ,8.0  ,4.0  ,6.0  ,2.0  ,                            & 
     &          4.0  ,6.0  ,2.0  ,4.0  ,6.0  ,8.0  / 
        DATA ENOIV/0.0  ,4.7920357e-02,8.824909,8.841201,8.864076,        & 
     &          15.73825,15.73998,20.37910,22.37705,22.40721,             & 
     &          28.67474,31.63571,31.63934,35.83378,35.83476,             & 
     &          44.33902,48.37428,48.38508,54.37857,54.39532,             & 
     &          54.42593,56.14158,56.17444,57.92984,57.94415,             & 
     &          58.03452,58.04428,58.06108,58.08709,58.79609,             & 
     &          59.33789,59.34961,59.36561,59.84372,59.87542,             & 
     &          60.23497,61.10992,61.36131,61.37108,61.38501,             & 
     &          61.40412,61.93150,61.93509,61.94088,61.94888,             & 
     &          62.18008,62.18691,62.46812,62.48219,62.49133,             & 
     &          63.30199,63.30286,63.32506,63.35387,63.75540,             & 
     &          63.77412,64.30924,64.30999,66.87376,67.85857,             & 
     &          67.86167,68.16618,68.17400,68.44416,68.44416,             & 
     &          68.50069,68.50069,68.74507,70.50282,70.51955,             & 
     &          70.55017,70.76975,70.76975,71.12993,71.15609,             & 
     &          71.21387,71.21387,71.31690,71.33785,71.39315,             & 
     &          71.39737,71.48887,71.50672,71.53300,72.12492,             & 
     &          72.12764,72.47591,72.50269,72.88482,72.88482,             & 
     &          73.16019,73.37047,73.37047,73.37047,73.37047,             & 
     &          73.52322,73.52322,73.52322,73.60108,73.61112,             & 
     &          73.64819,73.65725,73.68911,73.71453,73.93237,             & 
     &          73.95444,74.05078,74.06293,74.06293,74.10930,             & 
     &          74.12628,74.40265,74.40438,74.76035,74.76035,             & 
     &          74.76035,74.76035,75.18896,75.18896,75.18896,             & 
     &          76.30446,76.30806,76.44791,77.47625,77.47625,             & 
     &          77.92433,77.92433,78.12258,78.12258,78.19797,             & 
     &          78.21979,78.41159,78.43242,78.59385,78.59385,             & 
     &          78.59385,78.59385,78.63718,78.63718,78.63718,             & 
     &          78.85769,78.88398,78.91572,78.91572,78.96023,             & 
     &          78.97250,78.98019,80.20107,80.20107,80.72665,             & 
     &          80.72900,81.00314,81.01343,81.37509,81.37509,             & 
     &          81.37509,81.37509,81.42716,81.42716,81.83012,             & 
     &          82.88895,82.88895,83.03365,83.03365,83.03365,83.03365/ 
        DATA NOV/2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,         & 
     &          3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,          & 
     &          4,4,4,4,4,4,5,5,5,5,5,5,4,4,4,4,4,4,4,4,4,4,4,6,          & 
     &          6,6,6,6,4,4,4,6,4,4,4,4,4,7,7,7,7,7,8,8,8,8,5,5,          & 
     &          5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6/ 
        DATA GOV/1.0  ,1.0  ,3.0  ,5.0  ,3.0  ,                           & 
     &          1.0  ,3.0  ,5.0  ,5.0  ,1.0  ,                            & 
     &          3.0  ,1.0  ,3.0  ,1.0  ,3.0  ,                            & 
     &          5.0  ,3.0  ,5.0  ,7.0  ,5.0  ,                            & 
     &          1.0  ,3.0  ,5.0  ,3.0  ,3.0  ,                            & 
     &          3.0  ,5.0  ,7.0  ,3.0  ,1.0  ,                            & 
     &          3.0  ,5.0  ,5.0  ,5.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,1.0  ,5.0  ,3.0  ,                            & 
     &          1.0  ,7.0  ,3.0  ,3.0  ,1.0  ,                            & 
     &          1.0  ,3.0  ,5.0  ,3.0  ,3.0  ,                            & 
     &          5.0  ,7.0  ,5.0  ,7.0  ,3.0  ,                            & 
     &          3.0  ,3.0  ,5.0  ,7.0  ,5.0  ,                            & 
     &          3.0  ,3.0  ,3.0  ,5.0  ,7.0  ,                            & 
     &          3.0  ,1.0  ,3.0  ,5.0  ,5.0  ,                            & 
     &          5.0  ,3.0  ,7.0  ,3.0  ,5.0  ,                            & 
     &          7.0  ,3.0  ,5.0  ,7.0  ,5.0  ,                            & 
     &          5.0  ,3.0  ,1.0  ,7.0  ,3.0  ,                            & 
     &          3.0  ,3.0  ,5.0  ,7.0  ,5.0  ,                            & 
     &          3.0  ,3.0  ,5.0  ,7.0  ,3.0  ,                            & 
     &          3.0  ,5.0  ,7.0  ,1.0  ,3.0  ,                            & 
     &          5.0  ,5.0  ,5.0  ,3.0  ,5.0  ,                            & 
     &          7.0  ,7.0  ,3.0  ,3.0  ,5.0  ,                            & 
     &          7.0  ,1.0  ,3.0  ,5.0  ,5.0  / 
        DATA ENOV/0.0  ,10.18183,10.19878,10.23674,19.68863,              & 
     &          26.48845,26.50776,26.54108,28.73015,35.69651,             & 
     &          67.83862,69.59028,72.01395,72.28146,72.28596,             & 
     &          72.29554,74.50599,74.50733,74.50979,75.95557,             & 
     &          80.97483,80.99497,81.03748,82.38657,83.40436,             & 
     &          83.97941,84.00407,84.04314,84.82139,85.49855,             & 
     &          85.51269,85.53633,86.12596,86.43890,87.33036,             & 
     &          87.33829,87.35107,87.73579,87.80076,87.81837,             & 
     &          87.82866,88.39750,89.17985,89.60004,90.71603,             & 
     &          91.26665,91.26665,91.26888,91.48672,92.04689,             & 
     &          92.04763,92.04937,92.59937,92.97132,98.72523,             & 
     &          99.49233,106.7049,106.7049,106.7049,100.2237,             & 
     &          102.1987,102.8568,103.0377,103.0583,103.0944,             & 
     &          103.1870,103.5465,103.5465,103.5676,103.8792,             & 
     &          103.8829,104.1001,104.2509,104.2990,104.2990,             & 
     &          104.2990,104.3064,104.3181,104.3333,104.4087,             & 
     &          104.5556,104.5689,104.5754,105.0316,105.0733,             & 
     &          106.7358,106.9963,106.9963,106.9963,106.9274,             & 
     &          108.4187,108.5325,108.5325,108.5325,111.4108,             & 
     &          111.5461,111.5461,111.5461,111.7535,111.7535,             & 
     &          111.7535,111.8898,111.9082,112.1444,112.1444,             & 
     &          112.1444,112.3809,115.9379,116.0435,116.0435,             & 
     &          116.0435,116.1501,116.1501,116.1501,116.2166/ 
        DATA NOVI/2,2,2,3,3,3,3,3,4,4,4,4,4,4,4,5,5,5,5,5,6,6,6,6,        & 
     &          6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8/ 
        DATA GOVI/2.0  ,2.0  ,4.0  ,2.0  ,2.0  ,                          & 
     &          4.0  ,4.0  ,6.0  ,2.0  ,2.0  ,                            & 
     &          4.0  ,4.0  ,6.0  ,6.0  ,8.0  ,                            & 
     &          2.0  ,2.0  ,4.0  ,4.0  ,6.0  ,                            & 
     &          2.0  ,4.0  ,4.0  ,6.0  ,6.0  ,                            & 
     &          8.0  ,8.0  ,10.0  ,12.0  ,2.0  ,                          & 
     &          2.0  ,4.0  ,4.0  ,6.0  ,6.0  ,                            & 
     &          8.0  ,8.0  ,10.0  ,12.0  ,14.0  ,                         & 
     &          2.0  ,2.0  ,4.0  ,6.0  ,8.0  ,                            & 
     &          8.0  ,10.0  ,12.0  ,14.0  ,16.0  ,4.0  ,6.0  / 
        DATA ENOVI/0.0  ,11.94909,12.01505,79.35559,82.58831,             & 
     &          82.60773,83.64374,83.65008,105.7219,107.0408,             & 
     &          107.0487,107.4805,107.4831,107.5050,107.5062,             & 
     &          117.6237,118.2920,118.2920,118.5122,118.5122,             & 
     &          124.3735,124.3735,124.5034,124.5034,124.5142,             & 
     &          124.5142,124.5156,124.5156,124.5156,127.8017,             & 
     &          128.0311,128.0311,128.1171,128.1171,128.1243,             & 
     &          128.1243,128.1252,128.1252,128.1252,128.1252,             & 
     &          130.2520,130.3984,130.3984,130.4674,130.4674,             & 
     &          130.4680,130.4680,130.4680,130.4680,130.4680,             & 
     &          130.4693,130.4693/ 
        DATA NOVII/1,2,2,2,2,2,3,3,3,3,3,3,3,4,5,6/ 
        DATA GOVII/1.0  ,3.0  ,1.0  ,3.0  ,5.0  ,                         & 
     &          3.0  ,1.0  ,3.0  ,5.0  ,7.0  ,                            & 
     &          5.0  ,3.0  ,3.0  ,3.0  ,3.0  ,3.0  / 
        DATA ENOVII/0.0  ,561.0761,568.6182,568.6255,568.6938,            & 
     &          573.9532,664.1129,664.1129,664.1129,665.1804,             & 
     &          665.1804,665.1804,665.6218,697.8022,712.7239,720.8449/ 
        DATA SCI/4.179704,4.179868,4.180140,4.284864,4.411317,            & 
     &          4.556712,4.417538,4.418036,4.419087,4.460873,             & 
     &          5.012059,5.012145,5.012123,4.656621,4.682271,             & 
     &          4.682931,4.683973,4.715525,4.735114,4.735517,             & 
     &          4.736181,4.776610,4.823287,5.245868,4.960446,             & 
     &          4.636463,4.637096,4.638772,4.981131,4.981795,             & 
     &          4.983181,4.985224,4.985476,4.985839,4.648973,             & 
     &          4.987257,5.002644,5.026932,5.027268,5.027448,             & 
     &          4.751991,4.753114,4.754886,4.775015,4.807735,             & 
     &          4.820394,4.821310,4.822433,4.849118,4.880081,             & 
     &          4.964531,4.983084,4.983084,4.983084,4.988046,             & 
     &          4.988550,4.989272,4.739796,4.996660,5.002712,             & 
     &          5.007890,5.009317,5.009317,4.830776,4.830763,             & 
     &          4.830763,4.843919,4.885500,4.908801,4.963564,             & 
     &          4.983777,4.984663,4.984663,4.989659,4.989659,             & 
     &          4.992449,4.793379,4.998577,5.003253,5.004741,             & 
     &          5.006231,5.006231,4.972326,4.984214,4.985345,             & 
     &          4.985345,4.991673,4.991673,4.995098,4.831870,             & 
     &          5.000666,5.004451,5.005935,5.005935,5.004665,             & 
     &          4.984616,4.985763,4.985763,4.999064,4.999064,             & 
     &          4.999064,5.002865,5.004758,5.006233,5.006233,             & 
     &          5.006233,4.984145,4.984145,4.984432,5.002990,             & 
     &          5.002990,5.002990,5.005628,5.006508,5.006508,             & 
     &          5.006508,4.983365,4.983365,4.983365,5.006886,             & 
     &          5.006886,5.006886,5.008002,5.012074,5.012074,             & 
     &          5.012074,5.014099,5.014099,5.014099,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000/ 
        DATA SCII/3.322208,3.322644,3.633090,3.633257,3.633480,           & 
     &          3.893420,3.893440,4.089184,4.229172,4.229597,             & 
     &          3.436720,3.692591,3.692789,4.589103,3.953149,             & 
     &          3.953178,4.702748,4.702820,3.603431,3.770060,             & 
     &          3.770254,4.440431,4.441056,4.442240,3.961645,             & 
     &          3.961659,4.991754,4.992091,3.992479,3.992479,             & 
     &          3.697579,3.795690,3.796068,4.770876,4.780956,             & 
     &          3.968260,3.968260,3.994325,3.994325,3.754887,             & 
     &          4.893884,4.894430,4.895358,4.896706,4.906213,             & 
     &          4.906944,3.971236,3.971236,3.996578,3.996578,             & 
     &          5.011171,5.086056,5.086787,5.087796,5.188515,             & 
     &          5.190205,5.591661,5.735427,5.737656,5.740739,             & 
     &          5.745144,5.937513,5.941305,5.947732,5.956571,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000/ 
        DATA SCIII/2.247678,2.511178,2.511298,2.511595,2.783334,          & 
     &          2.988424,2.988601,2.988887,3.040348,3.275480,             & 
     &          2.516355,2.624130,2.770250,2.779440,2.779510,             & 
     &          2.779673,2.913505,2.912204,2.913576,3.001503,             & 
     &          3.471912,3.472452,3.473566,2.656192,3.501991,             & 
     &          2.707107,2.843331,2.843331,2.843447,3.667003,             & 
     &          2.928022,2.928408,2.928967,2.942070,2.942230,             & 
     &          2.942442,2.952919,2.960061,3.725865,3.726323,             & 
     &          3.727023,2.996535,3.802981,3.848412,3.848814,             & 
     &          3.849514,3.907525,3.915897,3.920174,3.920667,             & 
     &          3.921389,4.006182,3.997117,4.006877,2.756046,             & 
     &          4.057183,4.057739,4.058045,4.085242,2.876864,             & 
     &          2.910818,2.910818,2.910818,4.166811,2.957773,             & 
     &          2.957773,2.957773,2.998549,2.998549,2.998583,             & 
     &          2.998549,3.003525,3.005696,3.020441,3.020510,             & 
     &          3.020601,3.088544,2.797248,2.916939,2.968366,             & 
     &          2.968366,2.968366,3.000603,3.000603,3.000640,             & 
     &          3.003372,3.005378,3.009464,3.009464,3.009464,             & 
     &          3.027199,2.830505,2.926040,2.976524,2.976524,             & 
     &          2.976524,3.009360,2.932986,2.982088,2.982088,             & 
     &          2.982088,2.986294,2.986294,2.986294,4.828427,             & 
     &          4.828427,4.828427,5.151012,5.224114,5.224114,             & 
     &          5.227788,5.497263,5.497263,5.502663,5.756044,             & 
     &          5.817122,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000,6.000000,6.000000,6.000000,6.000000,             & 
     &          6.000000/ 
        DATA SCIV/1.644934,1.923884,1.924364,1.778341,1.948965,           & 
     &          1.949284,1.998202,1.998312,1.838941,1.962835,             & 
     &          1.963075,1.999589,1.999669,2.001418,2.001418,             & 
     &          1.874532,1.972046,1.972243,2.001394,2.001394,             & 
     &          2.002842,2.002845,2.003133,2.003133,1.897882,             & 
     &          1.978616,1.978616,2.003384,2.003384,2.004821,             & 
     &          2.004821,2.005002,2.005002,2.005023,2.005023,             & 
     &          1.913888,1.984033,1.984033,2.005113,2.005113,             & 
     &          2.007061,2.007061,2.007217,2.007217,2.007234,             & 
     &          2.007234,1.989961,1.989961,2.009654,2.009654,             & 
     &          2.009800,2.009800,2.009800,2.009800,2.009800/ 
        DATA SCV/0.6309066,0.7688928,0.9242349,0.9241881,0.9246764,       & 
     &          1.026124,1.003322,1.003322,1.003322,1.020118,             & 
     &          1.021842,1.027042,1.033988,1.051912,1.003001/ 
        DATA SNI/4.931870,5.108953,5.109031,5.204087,5.204087,            & 
     &          5.329969,5.330800,5.331948,5.401419,5.403777,             & 
     &          5.968683,5.969459,5.969806,5.605717,5.641185,             & 
     &          5.641868,5.642995,5.644538,5.662621,5.663186,             & 
     &          5.664362,5.702335,5.703489,5.705897,5.734946,             & 
     &          5.736104,5.797976,5.797737,5.588642,5.591228,             & 
     &          5.594790,5.615995,5.620492,5.980873,5.982464,             & 
     &          5.982872,5.983638,5.985012,5.986994,5.988774,             & 
     &          5.991692,5.990931,5.991989,5.990646,5.995945,             & 
     &          5.996395,5.996926,5.997295,6.001428,6.002394,             & 
     &          5.745163,5.761658,5.762813,5.764937,5.768167,             & 
     &          5.774817,5.775745,5.778173,5.802792,5.696104,             & 
     &          5.699982,5.706142,5.715098,5.722151,5.983984,             & 
     &          5.985277,5.987724,5.991767,5.985780,5.989671,             & 
     &          5.990610,5.994303,5.985830,5.987479,5.992273,             & 
     &          5.996771,5.993288,5.995173,5.998955,6.002261,             & 
     &          6.003886,6.255743,6.257061,6.360914,6.362586,             & 
     &          5.757127,5.763044,5.772635,5.779661,5.791555,             & 
     &          5.984846,5.986194,5.990250,5.997386,5.990476,             & 
     &          5.992170,5.992623,6.000597,5.993189,5.993189,             & 
     &          5.993189,6.000802,5.996590,5.998751,6.003087,             & 
     &          6.005032,6.007153,5.793718,5.804320,5.818227,             & 
     &          5.815205,5.821445,5.989322,5.989322,5.989322,             & 
     &          5.989322,5.992901,5.992901,5.992901,6.003714,             & 
     &          5.994695,5.997311,5.995185,6.008173,6.001116,             & 
     &          6.010741,6.005529,6.005529,6.008008,5.801159,             & 
     &          5.821038,5.833978,5.835981,5.835981,5.989852,             & 
     &          5.989852,5.989852,5.989852,5.993398,5.993398,             & 
     &          5.996287,5.996287,6.005342,6.014956,6.015613,             & 
     &          6.015613,6.015613,5.971639,5.971639,5.850568,             & 
     &          5.850568,5.850568,5.990061,5.990061,5.990061,             & 
     &          5.990061,5.991796,5.991796,5.993244,5.993244,             & 
     &          6.011374,6.018783,6.017594,6.017594,6.017594,             & 
     &          5.858175,5.858175,5.863377,5.863377,5.863377,             & 
     &          5.988659,5.988659,5.988659,5.988659,5.989390,             & 
     &          5.989390,5.994151,5.994151,6.020563,6.027375,             & 
     &          6.026996,6.026996,6.026996,5.866343,5.866343,             & 
     &          5.874646,5.874646,5.874646,5.990859,5.990859,             & 
     &          5.992667,5.992667,5.994933,5.994933,5.994933,             & 
     &          5.994933,6.030020,6.030020,6.038991,6.038991,             & 
     &          6.038991,5.873279,5.873279,5.877365,5.877365,             & 
     &          5.877365,5.992039,5.992039,5.996431,5.996431,             & 
     &          6.000838,6.000838,6.000838,6.000838,6.039686,             & 
     &          6.039686,6.042562,6.042562,6.042562,6.018730,             & 
     &          5.886326,5.994597,5.994597,6.047575,6.047575,             & 
     &          6.047575,6.078406,6.078406/ 
        DATA SNII/4.048939,4.049242,4.049750,4.145151,4.258360,           & 
     &          4.356432,4.688145,4.688257,4.688270,4.826217,             & 
     &          4.826217,4.826272,5.142618,4.284333,4.284811,             & 
     &          4.286872,4.288557,5.253337,4.532960,4.564950,             & 
     &          4.565974,4.567596,5.379367,4.605227,4.634192,             & 
     &          4.634804,4.635817,4.698180,4.771750,4.928997,             & 
     &          4.930175,4.931792,4.940517,4.947424,4.947905,             & 
     &          4.948512,4.976005,4.977055,4.977624,4.985716,             & 
     &          5.001774,4.517682,4.519204,4.522714,5.167543,             & 
     &          4.688999,4.706267,4.707887,4.710950,4.720615,             & 
     &          4.721386,4.723658,4.732426,4.791679,5.359476,             & 
     &          5.360877,5.362645,4.824183,4.939473,4.941747,             & 
     &          4.944790,4.948400,4.959554,4.960499,4.962098,             & 
     &          4.976268,4.977930,4.978850,4.988034,4.988136,             & 
     &          4.988983,4.989127,4.990716,4.997378,4.997656,             & 
     &          5.001124,4.999108,5.001566,5.001872,5.002052,             & 
     &          5.004650,5.004791,5.510713,5.511550,5.512880,             & 
     &          4.633615,4.635821,4.641956,4.663452,4.970949,             & 
     &          4.970949,4.970949,4.990886,4.990886,4.991119,             & 
     &          4.994713,4.999842,5.000126,5.004091,5.004655,             & 
     &          5.899771,5.900361,5.901458,5.903069,5.905087,             & 
     &          5.975470,5.976437,5.978201,6.162114,7.000000,             & 
     &          7.000000,7.000000,7.000000,7.000000,7.000000,             & 
     &          7.000000,7.000000,7.000000,7.000000,7.000000,             & 
     &          7.000000,7.000000/ 
        DATA SNIII/3.264886,3.265738,3.559230,3.555733,3.556163,          & 
     &          3.795859,3.795903,3.971288,4.062202,4.062887,             & 
     &          4.328298,4.441761,4.441879,3.362787,4.644640,             & 
     &          4.644671,3.648880,3.649321,3.924340,3.924418,             & 
     &          4.208216,4.209135,4.210838,4.353292,4.355043,             & 
     &          3.749472,4.546074,4.546964,4.554955,4.555551,             & 
     &          4.556746,4.558211,3.785649,3.786211,4.632736,             & 
     &          4.686664,4.687436,4.688480,3.926233,3.926969,             & 
     &          3.987034,3.987034,4.752837,4.754452,4.866728,             & 
     &          4.928826,4.929522,4.930549,4.931964,4.980142,             & 
     &          4.980413,4.980860,4.981436,3.664742,5.015918,             & 
     &          5.016470,5.050785,5.051935,5.052674,5.126587,             & 
     &          5.129026,3.959079,3.959144,5.192320,5.193925,             & 
     &          3.989434,3.989434,4.005147,4.005147,3.968565,             & 
     &          3.968565,3.992409,3.992409,4.006539,4.006539,             & 
     &          5.571515,5.574721,5.580696,6.132490,6.134099,             & 
     &          5.935632,5.939608,6.083616,6.087341,6.092149,             & 
     &          6.099411,6.364476,6.365573,6.178216,6.185985,             & 
     &          6.229222,6.316157,6.320952,6.326658,7.000000,             & 
     &          7.000000,7.000000,7.000000,7.000000,7.000000,             & 
     &          7.000000,7.000000,7.000000,7.000000,7.000000,             & 
     &          7.000000,7.000000,7.000000,7.000000,7.000000,             & 
     &          7.000000,7.000000,7.000000,7.000000,7.000000,             & 
     &          7.000000,7.000000,7.000000,7.000000,7.000000,             & 
     &          7.000000,7.000000,7.000000,7.000000,7.000000,             & 
     &          7.000000,7.000000,7.000000,7.000000,7.000000,             & 
     &          7.000000,7.000000,7.000000/ 
        DATA SNIV/2.226836,2.490623,2.490878,2.491461,2.755431,           & 
     &          2.952339,2.952669,2.953231,3.013266,3.231892,             & 
     &          2.493613,2.749589,2.762857,2.763010,2.763353,             & 
     &          2.901417,2.901452,2.901533,2.994477,3.382731,             & 
     &          3.383611,3.385459,3.472422,3.564919,3.607152,             & 
     &          3.607152,3.608737,3.645438,3.728391,3.728391,             & 
     &          2.639493,3.779903,3.797702,3.799535,3.799535,             & 
     &          3.799535,2.797721,2.797721,2.797721,3.872625,             & 
     &          3.873032,3.873596,3.883205,2.857107,2.934635,             & 
     &          2.934635,2.934635,3.951730,3.952443,3.952443,             & 
     &          2.993446,3.029914,3.030044,3.030246,4.061023,             & 
     &          3.127314,2.880353,2.950477,2.950477,2.950477,             & 
     &          2.998267,2.998267,2.998267,2.959708,2.959708,             & 
     &          2.959708,4.785085,4.873190,4.873190,4.874527,             & 
     &          7.000000,7.000000,7.000000/ 
        DATA SNV/1.634565,1.915400,1.916327,1.771111,1.943796,            & 
     &          1.944398,1.997854,1.998051,1.833396,1.959357,             & 
     &          1.959357,1.999384,1.999384,1.870467,1.969524,             & 
     &          1.969524,2.001982,2.001982,1.895972,1.977561,             & 
     &          1.977561,2.004889,2.004593,2.005843,2.005843,             & 
     &          2.006106,2.006106,2.006106,1.914798,1.983841,             & 
     &          1.983841,2.007186,2.007186,2.008574,2.008574,             & 
     &          2.008797,2.008797,2.008797,2.008797,1.929893,             & 
     &          1.990742,1.990742,2.010469,2.010469,2.011696,             & 
     &          2.011696,2.011930,2.011930,2.011930,2.011930,             & 
     &          2.011930/ 
        DATA SNVI/0.6290283,0.7657156,0.9208641,0.9208946,0.9217644,      & 
     &          0.9074407,1.024314,1.024867/ 
        DATA SOI/5.998809,6.000254,6.000874,6.149045,6.029965,            & 
     &          6.280362,6.354168,6.620857,6.620917,6.621027,             & 
     &          6.681873,6.681856,6.681878,6.554275,6.592577,             & 
     &          6.991942,6.991948,6.991948,6.991953,6.991953,             & 
     &          6.991953,6.994710,6.994710,6.994710,6.749977,             & 
     &          6.750016,6.750087,6.784759,6.784759,6.784759,             & 
     &          7.156593,7.157186,7.157570,6.676266,6.701954,             & 
     &          7.234456,6.993919,6.993919,6.993919,6.993919,             & 
     &          6.993919,6.997045,6.997045,6.997045,6.836974,             & 
     &          6.836974,6.836974,6.746833,6.766088,6.996467,             & 
     &          6.996467,6.996467,6.996467,6.996467,6.999115,             & 
     &          6.999115,6.999115,6.869720,6.869720,6.869720,             & 
     &          6.793480,6.808908,6.999084,6.999084,6.999084,             & 
     &          6.999084,6.999084,7.001479,7.001479,7.001479,             & 
     &          6.826997,6.839929,7.001803,7.001803,7.001803,             & 
     &          7.001803,7.001803,7.003955,7.003955,7.003955,             & 
     &          6.852827,6.864539,7.004705,8.000000,7.004705,             & 
     &          7.004705,7.004705,7.007905,7.007905,7.007905,             & 
     &          6.877349,6.883498,7.007761,7.007761,7.007761,             & 
     &          7.007761,7.007761,7.010591,7.010591,7.010591,             & 
     &          6.890946,6.900387,7.011453,7.011453,7.011453,             & 
     &          7.011453,7.011453,7.012788,7.012788,7.012788,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000/ 
        DATA SOII/4.784610,4.940430,4.940555,5.022952,5.022962,           & 
     &          5.557054,5.558274,5.558889,5.930025,5.929834,             & 
     &          5.160761,5.162283,5.164577,5.214052,5.216705,             & 
     &          6.210938,5.445367,5.490556,5.491464,5.492962,             & 
     &          5.495002,5.494491,5.494507,5.517109,5.517870,             & 
     &          5.519391,5.570158,5.573380,5.580988,6.392210,             & 
     &          6.394130,5.615287,5.616316,5.844496,5.880436,             & 
     &          5.880893,5.904250,5.904674,5.930839,5.931911,             & 
     &          5.933457,5.935489,5.954107,5.955576,5.956912,             & 
     &          5.954453,5.955388,5.959104,5.959794,5.959825,             & 
     &          5.959952,5.960805,5.964087,5.973599,5.975908,             & 
     &          5.993385,5.994448,5.442262,5.445265,5.449878,             & 
     &          5.491284,5.496742,6.232405,5.654758,5.656267,             & 
     &          5.658962,5.662895,5.725537,5.731195,5.738926,             & 
     &          5.741796,6.348958,6.348958,6.373241,6.373322,             & 
     &          6.408604,6.408635,6.419951,6.420006,6.439373,             & 
     &          6.439435,5.943564,5.943564,5.943564,5.947441,             & 
     &          5.961402,5.961402,5.961402,5.961402,6.471051,             & 
     &          5.965786,5.967088,5.967867,5.968222,5.972136,             & 
     &          5.972852,5.978758,6.488330,5.986873,5.991322,             & 
     &          5.986938,5.991353,5.994984,5.994948,5.989273,             & 
     &          5.989404,5.991879,5.997332,5.991945,5.997554,             & 
     &          5.992443,5.994410,6.001196,6.001346,6.002642,             & 
     &          6.003121,6.002740,6.003380,5.576062,5.580966,             & 
     &          5.588796,5.609913,5.619139,6.121709,6.121740,             & 
     &          5.734796,5.734796,5.738977,5.745944,5.750079,             & 
     &          5.750079,5.754774,5.772264,5.781077,5.960446,             & 
     &          5.960446,5.960446,5.960446,5.972282,5.974347,             & 
     &          5.974347,5.980535,5.980535,5.983981,5.987715,             & 
     &          5.990945,5.996456,5.985451,5.986923,5.990890,             & 
     &          6.000214,5.991003,6.000511,5.996286,5.996286,             & 
     &          6.000386,6.002267,6.003636,6.004436,6.003808,             & 
     &          6.004436,6.864734,6.865004,6.865458,6.871479,             & 
     &          6.874276,6.883227,6.883227,6.929313,6.929313,             & 
     &          6.945119,6.945119,6.947771,6.947771,7.214550,             & 
     &          7.214550,6.955942,6.960794,6.960794,6.963802,             & 
     &          6.963802,6.974759,6.974759,6.897856,6.897856,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000/ 
        DATA SOIII/3.980091,3.980605,3.981483,4.073126,4.181011,          & 
     &          4.263698,4.567496,2.851461,2.851508,4.688399,             & 
     &          4.688406,4.688483,4.944256,5.004756,5.087306,             & 
     &          4.201648,4.202927,4.205704,4.265077,5.589531,             & 
     &          5.591178,5.591910,4.466920,4.500863,4.502461,             & 
     &          4.505044,4.544430,5.702087,4.576288,4.577272,             & 
     &          4.578837,4.653335,4.743010,4.880211,4.882788,             & 
     &          4.885133,4.883788,4.916797,4.917453,4.918434,             & 
     &          4.946753,4.948285,4.949112,4.978528,4.991553,             & 
     &          5.071567,5.073311,5.091046,6.092469,5.236801,             & 
     &          5.238239,5.240906,4.450987,4.453166,4.458786,             & 
     &          4.490992,5.440957,5.477274,5.477840,5.478956,             & 
     &          5.480584,5.482658,4.640882,4.657492,4.659830,             & 
     &          4.664354,4.689624,5.526725,5.527675,5.529354,             & 
     &          4.742366,4.744359,4.746791,4.755242,4.803841,             & 
     &          5.629193,5.463434,5.633066,5.655169,4.904211,             & 
     &          4.904211,4.904211,4.911571,5.696494,5.696715,             & 
     &          5.697021,4.948280,4.949739,4.951245,4.983719,             & 
     &          4.983719,4.983719,4.985557,4.992922,4.595488,             & 
     &          4.595488,4.595488,4.614187,5.995187,5.995924,             & 
     &          5.997336,6.003933,6.004729,6.005913,6.007472,             & 
     &          6.009368,6.079757,6.079680,6.079593,6.079808,             & 
     &          6.081549,6.087021,6.088523,6.089349,6.127790,             & 
     &          6.130199,6.131379,6.150372,6.152512,6.155484,             & 
     &          4.922875,4.922875,4.922875,4.932409,4.922875,             & 
     &          4.922875,4.922875,4.991952,4.997716,6.251313,             & 
     &          6.251994,6.253141,4.947119,4.974443,4.974443,             & 
     &          4.974443,5.003925,6.782259,6.828280,6.541492,             & 
     &          6.547457,6.555665,6.965416,7.060263,7.160808,             & 
     &          7.160808,7.166230,7.174317,7.185196,7.256399,             & 
     &          7.261456,7.271212,7.771374,7.771374,7.771374,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000/ 
        DATA SOIV/3.228562,3.230039,3.508826,3.509360,3.510109,           & 
     &          3.741247,3.741307,3.904661,3.977057,3.978160,             & 
     &          4.214302,4.331145,4.331291,4.503491,4.503533,             & 
     &          3.322590,3.617383,3.618198,4.097019,4.098440,             & 
     &          4.101037,4.249484,4.252384,4.410741,4.412061,             & 
     &          4.420406,4.421309,4.422863,4.425270,4.491520,             & 
     &          4.543003,4.544125,4.545658,4.591770,4.594849,             & 
     &          3.506632,4.717019,4.742457,4.743450,4.744866,             & 
     &          4.746809,4.800908,4.801279,4.801878,4.802707,             & 
     &          4.826727,4.827440,4.856910,4.858391,4.859354,             & 
     &          3.927959,3.928085,4.948471,4.951597,4.995503,             & 
     &          4.997566,5.057139,5.057223,3.602068,5.487784,             & 
     &          5.488193,5.528638,5.529685,3.943577,3.943577,             & 
     &          3.956409,3.956409,5.607411,5.152442,5.155901,             & 
     &          5.162243,5.906104,5.906104,5.285100,5.290775,             & 
     &          3.955026,3.955026,5.325924,5.330538,6.007065,             & 
     &          6.007765,6.023023,6.026014,6.030425,6.132527,             & 
     &          6.133009,5.594399,5.600957,3.969004,3.969004,             & 
     &          5.768017,5.824149,5.824149,5.824149,5.824149,             & 
     &          5.865850,5.865850,5.865850,5.887424,5.890223,             & 
     &          5.900586,5.903125,5.912084,5.919259,5.981793,             & 
     &          5.988238,6.512458,4.040435,4.040435,6.034046,             & 
     &          6.039135,6.592915,6.593322,6.679720,6.679720,             & 
     &          6.679720,6.679720,6.791924,6.791924,6.791924,             & 
     &          7.150804,7.152208,7.208681,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000,             & 
     &          8.000000/ 
        DATA SOV/2.212299,2.477108,2.477559,2.478570,2.736373,            & 
     &          2.929941,2.930501,2.931468,2.995396,3.204502,             & 
     &          2.480139,2.586177,2.736415,2.753261,2.753545,             & 
     &          2.754149,2.895501,2.895588,2.895747,2.990361,             & 
     &          3.333698,3.335127,3.338143,3.434917,3.509306,             & 
     &          3.551885,3.553720,3.556629,3.614975,3.666382,             & 
     &          3.667461,3.669268,3.714561,3.738797,3.808601,             & 
     &          3.809227,3.810236,3.840735,3.845908,3.847311,             & 
     &          3.848131,3.893723,3.957266,2.655746,2.780047,             & 
     &          2.842480,2.842480,2.842734,2.867645,2.932266,             & 
     &          2.932352,2.932553,2.996815,3.040747,2.722722,             & 
     &          2.858081,4.369750,4.369750,4.369750,2.990545,             & 
     &          4.293693,4.399677,4.429359,4.432752,4.438707,             & 
     &          4.454041,4.514207,4.514207,4.517767,4.570812,             & 
     &          4.571450,2.913395,2.952782,2.965416,2.965416,             & 
     &          2.965416,4.644915,4.646958,4.649636,2.994349,             & 
     &          4.688903,4.691261,4.692408,4.774586,4.782194,             & 
     &          2.928606,3.022013,3.022013,3.022013,2.997125,             & 
     &          2.933282,2.986425,2.986425,2.986425,5.872365,             & 
     &          5.931636,5.931636,5.931636,6.025977,6.025977,             & 
     &          6.025977,6.090481,6.099398,6.217293,6.217293,             & 
     &          6.217293,6.343698,8.000000,8.000000,8.000000,             & 
     &          8.000000,8.000000,8.000000,8.000000,8.000000/ 
        DATA SOVI/1.626749,1.908750,1.910343,1.765577,1.939605,           & 
     &          1.940666,1.997515,1.997865,1.829541,1.956604,             & 
     &          1.957376,1.999562,1.999822,2.001964,2.002083,             & 
     &          1.867336,1.968342,1.968342,2.001995,2.001995,             & 
     &          1.976059,1.976059,2.004681,2.004681,2.007063,             & 
     &          2.007063,2.007365,2.007365,2.007365,1.914098,             & 
     &          1.982390,1.982390,2.008208,2.008208,2.010369,             & 
     &          2.010369,2.010631,2.010631,2.010631,2.010631,             & 
     &          1.930105,1.987143,1.987143,2.014184,2.014184,             & 
     &          2.014431,2.014431,2.014431,2.014431,2.014431,             & 
     &          2.014965,2.014965/ 
        DATA SOVII/0.6273875,0.7631111,0.9180546,0.9182081,0.9196253,     & 
     &          1.029738,0.9543552,0.9543552,0.9543552,1.004676,          & 
     &          1.004676,1.004676,1.025589,1.027917,1.034432,             & 
     &          1.045346/ 
! 
!       Find index for atom and ion, 10*IAT+IZI 
! 
        IF(IAT.GT.2.AND.IAT.LT.6.OR.                                      & 
     &    (IAT.LT.1.OR.IAT.GT.8)) THEN 
          U=0 
          WRITE(*,*) 'INVALID ATOM IN USER SUPPLIED ROUTINE PARTFUN' 
          STOP 
        END IF 
 
        IND=10*IAT+IZI 
 
        SELECT CASE(IND) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR HYDROGEN 
        CASE(11) 
           CALL PARTDV(T,ANE,ZH,MH,NHYD,GHYD,ENHYD,SHYD,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR HEI 
        CASE(21) 
           CALL PARTDV(T,ANE,ZHE,MHEI,NHEL,GHEL,ENHEL,SHEL,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR HEII 
        CASE(22) 
           CALL PARTDV(T,ANE,ZHE,MHEII,NHYD,GHYD,ENHYD,SHYD,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR CI 
        CASE(61) 
           CALL PARTDV(T,ANE,ZC,MCI,NCI,GCI,ENCI,SCI,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR CII 
        CASE(62) 
           CALL PARTDV(T,ANE,ZC,MCII,NCII,GCII,ENCII,SCII,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR CIII 
        CASE(63) 
           CALL PARTDV(T,ANE,ZC,MCIII,NCIII,GCIII,ENCIII,SCIII,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR CIV 
        CASE(64) 
           CALL PARTDV(T,ANE,ZC,MCIV,NCIV,GCIV,ENCIV,SCIV,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR CV 
        CASE(65) 
           CALL PARTDV(T,ANE,ZC,MCV,NCV,GCV,ENCV,SCV,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR CVI 
        CASE(66) 
           CALL PARTDV(T,ANE,ZC,MH,NHYD,GHYD,ENHYD,SHYD,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR NI 
        CASE(71) 
           CALL PARTDV(T,ANE,ZN,MNI,NNI,GNI,ENNI,SNI,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR NII 
        CASE(72) 
           CALL PARTDV(T,ANE,ZN,MNII,NNII,GNII,ENNII,SNII,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR NIII 
        CASE(73) 
           CALL PARTDV(T,ANE,ZN,MNIII,NNIII,GNIII,ENNIII,SNIII,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR NIV 
        CASE(74) 
           CALL PARTDV(T,ANE,ZN,MNIV,NNIV,GNIV,ENNIV,SNIV,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR NV 
        CASE(75) 
           CALL PARTDV(T,ANE,ZN,MNV,NNV,GNV,ENNV,SNV,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR NVI 
        CASE(76) 
           CALL PARTDV(T,ANE,ZN,MNVI,NNVI,GNVI,ENNVI,SNVI,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR NVII 
        CASE(77) 
           CALL PARTDV(T,ANE,ZN,MH,NHYD,GHYD,ENHYD,SHYD,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR OI 
        CASE(81) 
           CALL PARTDV(T,ANE,ZO,MOI,NOI,GOI,ENOI,SOI,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR OII 
        CASE(82) 
           CALL PARTDV(T,ANE,ZO,MOII,NOII,GOII,ENOII,SOII,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR OIII 
        CASE(83) 
           CALL PARTDV(T,ANE,ZO,MOIII,NOIII,GOIII,ENOIII,SOIII,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR OIV 
        CASE(84) 
           CALL PARTDV(T,ANE,ZO,MOIV,NOIV,GOIV,ENOIV,SOIV,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR OV 
        CASE(85) 
           CALL PARTDV(T,ANE,ZO,MOV,NOV,GOV,ENOV,SOV,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR OVI 
        CASE(86) 
           CALL PARTDV(T,ANE,ZO,MOVI,NOVI,GOVI,ENOVI,SOVI,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR OVII 
        CASE(87) 
           CALL PARTDV(T,ANE,ZO,MOVII,NOVII,GOVII,ENOVII,SOVII,U) 
! 
!       CALCULATING PARTITION FUNCTIONS FOR OVIII 
        CASE(88) 
           CALL PARTDV(T,ANE,ZO,MH,NHYD,GHYD,ENHYD,SHYD,U) 
        END SELECT 
! 
        RETURN 
        END SUBROUTINE PFSPEC 
! 
!     ************************************************************** 
! 
! 
        SUBROUTINE PARTDV(TEMP,DNE,Z,NLEV,NE,GEE,ENRGY,S,U) 
!       =================================================== 
! 
        use accura 
        use params 
        implicit real(dp) (a-h,o-z),logical (l) 
 
        REAL(DP) :: GEE(*),ENRGY(*),S(*) 
        INTEGER  :: NE(*) 
 
        U=0.0 
        ET=TEMP/11604.8 
        P=(14.69e0-0.20-0.6667*LOG10(DNE)) 
! 
        DO I=1,NLEV 
           U1=FLOAT(NE(I)) 
           ZSTAR=Z-S(I) 
           IF (ZSTAR.GT.0)THEN 
               W=P+4.*LOG10(ZSTAR)-4.*LOG10(U1) 
            ELSE 
               W=0.0 
           ENDIF 
           IF (W.GT.1.) W=1. 
! 
           IF ((ENRGY(I)/ET).LT.65.0) THEN 
                U1=GEE(I)*W*EXP(-ENRGY(I)/ET) 
            ELSE 
                U1=0.0 
           ENDIF 
           U=U+U1 
        END DO 
        RETURN 
        END SUBROUTINE PARTDV 
! 
!     ************************************************************** 
! 
      subroutine pfni(ion,t,pf,dut,dun) 
!     ================================= 
! 
!     partition functions for Ni IV to Ni IX 
! 
!     this routine interpolates within a grid 
!     calculated from all levels predicted by 
!     Kurucz (1992), i.e. over 12,000 levels per ion. 
!     the partition functions depend only on T ! 
!     (i.e. no level dissolution with increasing density) 
!     TL  27-DEC-1994, 23-JAN-1995 
! 
!     Output:  PF   partition function 
!              DUT  d(PF)/dT 
!              DUN  d(PF)/d(ANE) (=0 in this version) 
! 
      use accura 
      implicit real(dp) (a-h,o-z) 
! 
      real(dp) :: g0(6) 
      real(dp) :: p4a(190),p4b(170) 
      real(dp) :: p5a(190),p5b(170) 
      real(dp) :: p6a(190),p6b(170) 
      real(dp) :: p7a(190),p7b(170) 
      real(dp) :: p8a(190),p8b(170) 
      real(dp) :: p9a(190),p9b(170) 
      real(dp), parameter :: xen=2.302585093,xmil=0.001 
! 
      data g0/28.,25.,6.,25.,28.,21./ 
! 
      data p4a/                                                           & 
     &    1.447,1.464,1.482,1.501,1.518,1.535,1.551,1.567,1.582,1.596,    & 
     &    1.610,1.623,1.636,1.648,1.659,1.671,1.681,1.692,1.702,1.711,    & 
     &    1.721,1.730,1.739,1.748,1.757,1.765,1.774,1.782,1.791,1.799,    & 
     &    1.808,1.816,1.824,1.833,1.841,1.850,1.859,1.868,1.877,1.886,    & 
     &    1.895,1.905,1.914,1.924,1.934,1.945,1.955,1.966,1.977,1.989,    & 
     &    2.000,2.012,2.025,2.037,2.050,2.063,2.077,2.091,2.105,2.119,    & 
     &    2.134,2.149,2.164,2.179,2.195,2.211,2.227,2.243,2.260,2.276,    & 
     &    2.293,2.310,2.327,2.344,2.362,2.379,2.397,2.414,2.432,2.449,    & 
     &    2.467,2.484,2.502,2.519,2.537,2.554,2.571,2.588,2.606,2.623,    & 
     &    2.640,2.657,2.674,2.690,2.707,2.723,2.740,2.756,2.772,2.788,    & 
     &    2.804,2.819,2.835,2.850,2.866,2.881,2.896,2.911,2.925,2.940,    & 
     &    2.954,2.969,2.983,2.997,3.010,3.024,3.038,3.051,3.064,3.077,    & 
     &    3.090,3.103,3.116,3.128,3.141,3.153,3.165,3.177,3.189,3.201,    & 
     &    3.213,3.224,3.235,3.247,3.258,3.269,3.280,3.291,3.301,3.312,    & 
     &    3.322,3.332,3.343,3.353,3.363,3.373,3.382,3.392,3.402,3.411,    & 
     &    3.421,3.430,3.439,3.448,3.457,3.466,3.475,3.484,3.492,3.501,    & 
     &    3.509,3.518,3.526,3.534,3.542,3.550,3.558,3.566,3.574,3.582,    & 
     &    3.589,3.597,3.604,3.612,3.619,3.626,3.634,3.641,3.648,3.655,    & 
     &    3.662,3.669,3.676,3.682,3.689,3.696,3.702,3.709,3.715,3.722/ 
      data p4b/                                                           & 
     &    3.589,3.597,3.604,3.612,3.619,3.626,3.634,3.641,3.648,3.655,    & 
     &    3.662,3.669,3.676,3.682,3.689,3.696,3.702,3.709,3.715,3.722,    & 
     &    3.728,3.734,3.740,3.747,3.753,3.759,3.765,3.771,3.777,3.782,    & 
     &    3.788,3.794,3.800,3.805,3.811,3.816,3.822,3.827,3.833,3.838,    & 
     &    3.843,3.849,3.854,3.859,3.864,3.869,3.874,3.879,3.884,3.889,    & 
     &    3.894,3.899,3.904,3.909,3.913,3.918,3.923,3.927,3.932,3.936,    & 
     &    3.941,3.945,3.950,3.954,3.959,3.963,3.967,3.972,3.976,3.980,    & 
     &    3.984,3.988,3.993,3.997,4.001,4.005,4.009,4.013,4.017,4.021,    & 
     &    4.024,4.028,4.032,4.036,4.040,4.043,4.047,4.051,4.055,4.058,    & 
     &    4.062,4.065,4.069,4.072,4.076,4.079,4.083,4.086,4.090,4.093,    & 
     &    4.097,4.100,4.103,4.107,4.110,4.113,4.116,4.120,4.123,4.126,    & 
     &    4.129,4.132,4.135,4.138,4.141,4.144,4.148,4.151,4.154,4.157,    & 
     &    4.159,4.162,4.165,4.168,4.171,4.174,4.177,4.180,4.182,4.185,    & 
     &    4.188,4.191,4.193,4.196,4.199,4.202,4.204,4.207,4.210,4.212,    & 
     &    4.215,4.217,4.220,4.223,4.225,4.228,4.230,4.233,4.235,4.238,    & 
     &    4.240,4.243,4.245,4.247,4.250,4.252,4.255,4.257,4.259,4.262,    & 
     &    4.264,4.266,4.268,4.271,4.273,4.275,4.278,4.280,4.282,4.284/ 
      data p5a/                                                           & 
     &    1.398,1.408,1.427,1.446,1.466,1.486,1.506,1.526,1.545,1.564,    & 
     &    1.583,1.601,1.619,1.636,1.652,1.668,1.683,1.698,1.712,1.725,    & 
     &    1.738,1.751,1.763,1.775,1.786,1.797,1.808,1.818,1.828,1.837,    & 
     &    1.846,1.855,1.864,1.873,1.881,1.889,1.897,1.904,1.912,1.919,    & 
     &    1.926,1.933,1.940,1.946,1.953,1.960,1.966,1.972,1.979,1.985,    & 
     &    1.991,1.997,2.003,2.009,2.016,2.022,2.028,2.034,2.040,2.046,    & 
     &    2.052,2.058,2.065,2.071,2.077,2.084,2.090,2.097,2.103,2.110,    & 
     &    2.117,2.124,2.131,2.138,2.145,2.152,2.160,2.167,2.175,2.183,    & 
     &    2.191,2.199,2.207,2.216,2.224,2.233,2.241,2.250,2.259,2.268,    & 
     &    2.278,2.287,2.297,2.306,2.316,2.326,2.336,2.346,2.356,2.367,    & 
     &    2.377,2.387,2.398,2.409,2.419,2.430,2.441,2.452,2.463,2.474,    & 
     &    2.485,2.497,2.508,2.519,2.530,2.542,2.553,2.564,2.576,2.587,    & 
     &    2.599,2.610,2.621,2.633,2.644,2.655,2.667,2.678,2.689,2.701,    & 
     &    2.712,2.723,2.734,2.745,2.757,2.768,2.779,2.790,2.801,2.812,    & 
     &    2.822,2.833,2.844,2.855,2.865,2.876,2.886,2.897,2.907,2.918,    & 
     &    2.928,2.938,2.948,2.958,2.968,2.978,2.988,2.998,3.008,3.018,    & 
     &    3.027,3.037,3.046,3.056,3.065,3.075,3.084,3.093,3.102,3.111,    & 
     &    3.120,3.129,3.138,3.147,3.156,3.164,3.173,3.182,3.190,3.198,    & 
     &    3.207,3.215,3.223,3.232,3.240,3.248,3.256,3.264,3.272,3.279/ 
      data p5b/                                                           & 
     &    3.120,3.129,3.138,3.147,3.156,3.164,3.173,3.182,3.190,3.198,    & 
     &    3.207,3.215,3.223,3.232,3.240,3.248,3.256,3.264,3.272,3.279,    & 
     &    3.287,3.295,3.303,3.310,3.318,3.325,3.333,3.340,3.347,3.355,    & 
     &    3.362,3.369,3.376,3.383,3.390,3.397,3.404,3.411,3.417,3.424,    & 
     &    3.431,3.438,3.444,3.451,3.457,3.464,3.470,3.476,3.483,3.489,    & 
     &    3.495,3.501,3.507,3.514,3.520,3.526,3.531,3.537,3.543,3.549,    & 
     &    3.555,3.561,3.566,3.572,3.578,3.583,3.589,3.594,3.600,3.605,    & 
     &    3.610,3.616,3.621,3.626,3.632,3.637,3.642,3.647,3.652,3.657,    & 
     &    3.662,3.667,3.672,3.677,3.682,3.687,3.692,3.697,3.701,3.706,    & 
     &    3.711,3.716,3.720,3.725,3.729,3.734,3.738,3.743,3.747,3.752,    & 
     &    3.756,3.761,3.765,3.769,3.774,3.778,3.782,3.786,3.790,3.795,    & 
     &    3.799,3.803,3.807,3.811,3.815,3.819,3.823,3.827,3.831,3.835,    & 
     &    3.839,3.843,3.846,3.850,3.854,3.858,3.862,3.865,3.869,3.873,    & 
     &    3.876,3.880,3.884,3.887,3.891,3.894,3.898,3.901,3.905,3.908,    & 
     &    3.912,3.915,3.918,3.922,3.925,3.929,3.932,3.935,3.939,3.942,    & 
     &    3.945,3.948,3.951,3.955,3.958,3.961,3.964,3.967,3.970,3.974,    & 
     &    3.977,3.980,3.983,3.986,3.989,3.992,3.995,3.998,4.001,4.004/ 
      data p6a/                                                           & 
     &    0.778,0.804,0.817,0.834,0.854,0.876,0.901,0.928,0.957,0.987,    & 
     &    1.017,1.048,1.079,1.109,1.139,1.169,1.197,1.225,1.253,1.279,    & 
     &    1.304,1.329,1.353,1.376,1.398,1.419,1.440,1.459,1.478,1.497,    & 
     &    1.515,1.532,1.548,1.564,1.580,1.594,1.609,1.623,1.636,1.649,    & 
     &    1.662,1.674,1.686,1.698,1.709,1.720,1.730,1.740,1.750,1.760,    & 
     &    1.769,1.779,1.788,1.796,1.805,1.813,1.821,1.829,1.837,1.845,    & 
     &    1.852,1.860,1.867,1.874,1.881,1.888,1.894,1.901,1.907,1.914,    & 
     &    1.920,1.926,1.932,1.938,1.944,1.950,1.956,1.962,1.968,1.974,    & 
     &    1.979,1.985,1.991,1.996,2.002,2.007,2.013,2.018,2.024,2.029,    & 
     &    2.035,2.041,2.046,2.052,2.057,2.063,2.068,2.074,2.080,2.086,    & 
     &    2.091,2.097,2.103,2.109,2.115,2.121,2.127,2.133,2.139,2.145,    & 
     &    2.152,2.158,2.164,2.171,2.177,2.184,2.190,2.197,2.204,2.211,    & 
     &    2.218,2.225,2.232,2.239,2.246,2.253,2.261,2.268,2.276,2.283,    & 
     &    2.291,2.298,2.306,2.314,2.322,2.330,2.338,2.346,2.354,2.362,    & 
     &    2.370,2.379,2.387,2.395,2.404,2.412,2.420,2.429,2.438,2.446,    & 
     &    2.455,2.463,2.472,2.481,2.489,2.498,2.507,2.516,2.524,2.533,    & 
     &    2.542,2.551,2.560,2.569,2.577,2.586,2.595,2.604,2.613,2.622,    & 
     &    2.631,2.639,2.648,2.657,2.666,2.675,2.683,2.692,2.701,2.710,    & 
     &    2.718,2.727,2.736,2.744,2.753,2.761,2.770,2.779,2.787,2.796/ 
      data p6b/                                                           & 
     &    2.631,2.639,2.648,2.657,2.666,2.675,2.683,2.692,2.701,2.710,    & 
     &    2.718,2.727,2.736,2.744,2.753,2.761,2.770,2.779,2.787,2.796,    & 
     &    2.804,2.812,2.821,2.829,2.838,2.846,2.854,2.862,2.871,2.879,    & 
     &    2.887,2.895,2.903,2.911,2.919,2.927,2.935,2.943,2.951,2.958,    & 
     &    2.966,2.974,2.982,2.989,2.997,3.005,3.012,3.020,3.027,3.035,    & 
     &    3.042,3.049,3.057,3.064,3.071,3.078,3.086,3.093,3.100,3.107,    & 
     &    3.114,3.121,3.128,3.135,3.141,3.148,3.155,3.162,3.169,3.175,    & 
     &    3.182,3.188,3.195,3.202,3.208,3.214,3.221,3.227,3.234,3.240,    & 
     &    3.246,3.252,3.259,3.265,3.271,3.277,3.283,3.289,3.295,3.301,    & 
     &    3.307,3.313,3.319,3.325,3.330,3.336,3.342,3.348,3.353,3.359,    & 
     &    3.364,3.370,3.376,3.381,3.386,3.392,3.397,3.403,3.408,3.413,    & 
     &    3.419,3.424,3.429,3.434,3.440,3.445,3.450,3.455,3.460,3.465,    & 
     &    3.470,3.475,3.480,3.485,3.490,3.495,3.499,3.504,3.509,3.514,    & 
     &    3.518,3.523,3.528,3.533,3.537,3.542,3.546,3.551,3.555,3.560,    & 
     &    3.564,3.569,3.573,3.578,3.582,3.586,3.591,3.595,3.599,3.604,    & 
     &    3.608,3.612,3.616,3.621,3.625,3.629,3.633,3.637,3.641,3.645,    & 
     &    3.649,3.653,3.657,3.661,3.665,3.669,3.673,3.677,3.681,3.685/ 
      data p7a/                                                           & 
     &    1.398,1.398,1.398,1.398,1.406,1.425,1.443,1.461,1.480,1.498,    & 
     &    1.516,1.534,1.551,1.568,1.585,1.601,1.616,1.631,1.646,1.660,    & 
     &    1.674,1.687,1.700,1.712,1.724,1.736,1.747,1.758,1.768,1.778,    & 
     &    1.788,1.797,1.806,1.815,1.824,1.832,1.840,1.848,1.855,1.863,    & 
     &    1.870,1.877,1.883,1.890,1.896,1.902,1.908,1.914,1.920,1.925,    & 
     &    1.931,1.936,1.941,1.946,1.951,1.956,1.960,1.965,1.969,1.974,    & 
     &    1.978,1.982,1.986,1.990,1.994,1.998,2.001,2.005,2.009,2.012,    & 
     &    2.016,2.019,2.022,2.026,2.029,2.032,2.035,2.038,2.041,2.044,    & 
     &    2.047,2.050,2.053,2.056,2.059,2.061,2.064,2.067,2.069,2.072,    & 
     &    2.075,2.077,2.080,2.082,2.085,2.088,2.090,2.093,2.095,2.098,    & 
     &    2.100,2.103,2.105,2.107,2.110,2.112,2.115,2.117,2.120,2.122,    & 
     &    2.125,2.127,2.130,2.132,2.135,2.137,2.140,2.142,2.145,2.148,    & 
     &    2.150,2.153,2.155,2.158,2.161,2.163,2.166,2.169,2.172,2.175,    & 
     &    2.178,2.180,2.183,2.186,2.189,2.192,2.195,2.198,2.202,2.205,    & 
     &    2.208,2.211,2.215,2.218,2.221,2.225,2.228,2.232,2.235,2.239,    & 
     &    2.243,2.246,2.250,2.254,2.258,2.261,2.265,2.269,2.273,2.277,    & 
     &    2.282,2.286,2.290,2.294,2.299,2.303,2.307,2.312,2.316,2.321,    & 
     &    2.325,2.330,2.335,2.339,2.344,2.349,2.354,2.359,2.364,2.369,    & 
     &    2.374,2.379,2.384,2.389,2.394,2.399,2.405,2.410,2.415,2.420/ 
      data p7b/                                                           & 
     &    2.325,2.330,2.335,2.339,2.344,2.349,2.354,2.359,2.364,2.369,    & 
     &    2.374,2.379,2.384,2.389,2.394,2.399,2.405,2.410,2.415,2.420,    & 
     &    2.426,2.431,2.437,2.442,2.448,2.453,2.459,2.464,2.470,2.476,    & 
     &    2.481,2.487,2.493,2.498,2.504,2.510,2.516,2.521,2.527,2.533,    & 
     &    2.539,2.545,2.551,2.556,2.562,2.568,2.574,2.580,2.586,2.592,    & 
     &    2.598,2.604,2.610,2.616,2.622,2.628,2.634,2.640,2.646,2.652,    & 
     &    2.658,2.664,2.670,2.676,2.682,2.687,2.693,2.699,2.705,2.711,    & 
     &    2.717,2.723,2.729,2.735,2.741,2.747,2.753,2.759,2.764,2.770,    & 
     &    2.776,2.782,2.788,2.794,2.799,2.805,2.811,2.817,2.823,2.828,    & 
     &    2.834,2.840,2.846,2.851,2.857,2.863,2.868,2.874,2.879,2.885,    & 
     &    2.891,2.896,2.902,2.907,2.913,2.918,2.924,2.929,2.935,2.940,    & 
     &    2.945,2.951,2.956,2.962,2.967,2.972,2.978,2.983,2.988,2.993,    & 
     &    2.999,3.004,3.009,3.014,3.019,3.025,3.030,3.035,3.040,3.045,    & 
     &    3.050,3.055,3.060,3.065,3.070,3.075,3.080,3.085,3.090,3.095,    & 
     &    3.099,3.104,3.109,3.114,3.119,3.123,3.128,3.133,3.138,3.142,    & 
     &    3.147,3.152,3.156,3.161,3.165,3.170,3.175,3.179,3.184,3.188,    & 
     &    3.193,3.197,3.202,3.206,3.210,3.215,3.219,3.224,3.228,3.232/ 
      data p8a/                                                           & 
     &    1.447,1.447,1.447,1.447,1.447,1.447,1.459,1.475,1.489,1.504,    & 
     &    1.518,1.531,1.544,1.556,1.568,1.580,1.591,1.602,1.612,1.622,    & 
     &    1.631,1.640,1.649,1.658,1.666,1.674,1.682,1.689,1.696,1.703,    & 
     &    1.710,1.716,1.722,1.728,1.734,1.740,1.745,1.751,1.756,1.761,    & 
     &    1.766,1.770,1.775,1.779,1.784,1.788,1.792,1.796,1.800,1.804,    & 
     &    1.807,1.811,1.814,1.818,1.821,1.824,1.827,1.831,1.834,1.836,    & 
     &    1.839,1.842,1.845,1.848,1.850,1.853,1.855,1.858,1.860,1.863,    & 
     &    1.865,1.867,1.870,1.872,1.874,1.876,1.878,1.880,1.882,1.884,    & 
     &    1.886,1.888,1.890,1.892,1.894,1.896,1.898,1.900,1.902,1.903,    & 
     &    1.905,1.907,1.909,1.911,1.912,1.914,1.916,1.917,1.919,1.921,    & 
     &    1.923,1.924,1.926,1.928,1.929,1.931,1.933,1.934,1.936,1.938,    & 
     &    1.939,1.941,1.943,1.945,1.946,1.948,1.950,1.951,1.953,1.955,    & 
     &    1.957,1.959,1.960,1.962,1.964,1.966,1.968,1.970,1.971,1.973,    & 
     &    1.975,1.977,1.979,1.981,1.983,1.985,1.987,1.989,1.991,1.993,    & 
     &    1.995,1.998,2.000,2.002,2.004,2.006,2.009,2.011,2.013,2.015,    & 
     &    2.018,2.020,2.023,2.025,2.027,2.030,2.032,2.035,2.037,2.040,    & 
     &    2.043,2.045,2.048,2.051,2.053,2.056,2.059,2.062,2.064,2.067,    & 
     &    2.070,2.073,2.076,2.079,2.082,2.085,2.088,2.091,2.094,2.097,    & 
     &    2.100,2.103,2.107,2.110,2.113,2.116,2.120,2.123,2.126,2.130/ 
      data p8b/                                                           & 
     &    2.070,2.073,2.076,2.079,2.082,2.085,2.088,2.091,2.094,2.097,    & 
     &    2.100,2.103,2.107,2.110,2.113,2.116,2.120,2.123,2.126,2.130,    & 
     &    2.133,2.137,2.140,2.143,2.147,2.151,2.154,2.158,2.161,2.165,    & 
     &    2.168,2.172,2.176,2.180,2.183,2.187,2.191,2.195,2.198,2.202,    & 
     &    2.206,2.210,2.214,2.218,2.222,2.226,2.230,2.233,2.237,2.241,    & 
     &    2.245,2.250,2.254,2.258,2.262,2.266,2.270,2.274,2.278,2.282,    & 
     &    2.286,2.291,2.295,2.299,2.303,2.307,2.312,2.316,2.320,2.324,    & 
     &    2.329,2.333,2.337,2.341,2.346,2.350,2.354,2.359,2.363,2.367,    & 
     &    2.371,2.376,2.380,2.384,2.389,2.393,2.397,2.402,2.406,2.410,    & 
     &    2.415,2.419,2.423,2.428,2.432,2.436,2.441,2.445,2.449,2.454,    & 
     &    2.458,2.462,2.467,2.471,2.475,2.480,2.484,2.488,2.493,2.497,    & 
     &    2.501,2.506,2.510,2.514,2.519,2.523,2.527,2.531,2.536,2.540,    & 
     &    2.544,2.548,2.553,2.557,2.561,2.565,2.570,2.574,2.578,2.582,    & 
     &    2.586,2.591,2.595,2.599,2.603,2.607,2.611,2.616,2.620,2.624,    & 
     &    2.628,2.632,2.636,2.640,2.644,2.648,2.652,2.656,2.661,2.665,    & 
     &    2.669,2.673,2.677,2.681,2.685,2.689,2.693,2.696,2.700,2.704,    & 
     &    2.708,2.712,2.716,2.720,2.724,2.728,2.732,2.736,2.739,2.743/ 
      data p9a/                                                           & 
     &    1.322,1.322,1.322,1.322,1.322,1.322,1.322,1.322,1.322,1.325,    & 
     &    1.334,1.342,1.351,1.358,1.366,1.373,1.380,1.386,1.392,1.398,    & 
     &    1.404,1.409,1.415,1.420,1.425,1.429,1.434,1.438,1.442,1.446,    & 
     &    1.450,1.454,1.457,1.461,1.464,1.467,1.470,1.473,1.476,1.479,    & 
     &    1.482,1.485,1.487,1.490,1.492,1.495,1.497,1.499,1.501,1.503,    & 
     &    1.505,1.507,1.509,1.511,1.513,1.515,1.517,1.519,1.520,1.522,    & 
     &    1.524,1.525,1.527,1.528,1.530,1.531,1.533,1.534,1.535,1.537,    & 
     &    1.538,1.539,1.541,1.542,1.543,1.545,1.546,1.547,1.548,1.549,    & 
     &    1.551,1.552,1.553,1.554,1.555,1.556,1.558,1.559,1.560,1.561,    & 
     &    1.562,1.563,1.565,1.566,1.567,1.568,1.569,1.570,1.571,1.573,    & 
     &    1.574,1.575,1.576,1.577,1.579,1.580,1.581,1.582,1.584,1.585,    & 
     &    1.586,1.588,1.589,1.590,1.592,1.593,1.594,1.596,1.597,1.599,    & 
     &    1.600,1.602,1.603,1.605,1.606,1.608,1.609,1.611,1.612,1.614,    & 
     &    1.616,1.617,1.619,1.621,1.622,1.624,1.626,1.628,1.630,1.631,    & 
     &    1.633,1.635,1.637,1.639,1.641,1.643,1.645,1.647,1.649,1.651,    & 
     &    1.653,1.655,1.657,1.659,1.661,1.664,1.666,1.668,1.670,1.673,    & 
     &    1.675,1.677,1.679,1.682,1.684,1.686,1.689,1.691,1.694,1.696,    & 
     &    1.699,1.701,1.704,1.706,1.709,1.711,1.714,1.716,1.719,1.722,    & 
     &    1.724,1.727,1.729,1.732,1.735,1.738,1.740,1.743,1.746,1.749/ 
      data p9b/                                                           & 
     &    1.699,1.701,1.704,1.706,1.709,1.711,1.714,1.716,1.719,1.722,    & 
     &    1.724,1.727,1.729,1.732,1.735,1.738,1.740,1.743,1.746,1.749,    & 
     &    1.751,1.754,1.757,1.760,1.763,1.765,1.768,1.771,1.774,1.777,    & 
     &    1.780,1.783,1.786,1.789,1.792,1.795,1.798,1.801,1.804,1.807,    & 
     &    1.810,1.813,1.816,1.819,1.822,1.825,1.828,1.831,1.834,1.837,    & 
     &    1.840,1.843,1.847,1.850,1.853,1.856,1.859,1.862,1.865,1.869,    & 
     &    1.872,1.875,1.878,1.881,1.884,1.888,1.891,1.894,1.897,1.901,    & 
     &    1.904,1.907,1.910,1.913,1.917,1.920,1.923,1.926,1.930,1.933,    & 
     &    1.936,1.939,1.943,1.946,1.949,1.952,1.956,1.959,1.962,1.965,    & 
     &    1.969,1.972,1.975,1.978,1.982,1.985,1.988,1.992,1.995,1.998,    & 
     &    2.001,2.005,2.008,2.011,2.014,2.018,2.021,2.024,2.027,2.031,    & 
     &    2.034,2.037,2.040,2.044,2.047,2.050,2.053,2.057,2.060,2.063,    & 
     &    2.066,2.070,2.073,2.076,2.079,2.083,2.086,2.089,2.092,2.095,    & 
     &    2.099,2.102,2.105,2.108,2.111,2.115,2.118,2.121,2.124,2.127,    & 
     &    2.131,2.134,2.137,2.140,2.143,2.146,2.149,2.153,2.156,2.159,    & 
     &    2.162,2.165,2.168,2.171,2.175,2.178,2.181,2.184,2.187,2.190,    & 
     &    2.193,2.196,2.199,2.202,2.205,2.208,2.212,2.215,2.218,2.221/ 
! 
      if(t.lt.12000.) then 
        pf=g0(ion-3) 
        dut=0. 
        dun=0. 
        return 
      endif 
! 
      it=int(t/1000) 
      if(it.ge.350) it=349 
      t1=1000.*it 
      t2=t1+1000. 
      select case(ion) 
        case(4) 
        if(t.le.200000.) then 
          xu1=p4a(it-10) 
          xu2=p4a(it-9) 
        else 
          xu1=p4b(it-180) 
          xu2=p4b(it-179) 
        endif 
      case(5) 
        if(t.le.200000.) then 
          xu1=p5a(it-10) 
          xu2=p5a(it-9) 
        else 
          xu1=p5b(it-180) 
          xu2=p5b(it-179) 
        endif 
      case(6) 
        if(t.le.200000.) then 
          xu1=p6a(it-10) 
          xu2=p6a(it-9) 
        else 
          xu1=p6b(it-180) 
          xu2=p6b(it-179) 
        endif 
      case(7) 
        if(t.le.200000.) then 
          xu1=p7a(it-10) 
          xu2=p7a(it-9) 
        else 
          xu1=p7b(it-180) 
          xu2=p7b(it-179) 
        endif 
      case(8) 
        if(t.le.200000.) then 
          xu1=p8a(it-10) 
          xu2=p8a(it-9) 
        else 
          xu1=p8b(it-180) 
          xu2=p8b(it-179) 
        endif 
      case(9) 
        if(t.le.200000.) then 
          xu1=p9a(it-10) 
          xu2=p9a(it-9) 
        else 
          xu1=p9b(it-180) 
          xu2=p9b(it-179) 
        endif 
      end select 
! 
      dxt=xmil*(xu2-xu1) 
      xu=xu1+(t-t1)*dxt 
      pf=exp(xen*xu) 
      dut=xen*pf*dxt 
      dun=0. 
      return 
      end subroutine pfni 
! 
! ****************************************************************** 
! 
! 
      SUBROUTINE PFHEAV(IIZ,JNION,MODE,t,ane,u) 
!     ========================================= 
! 
!     subset of kurucz's pfsaha for Z>28. 
!     removed code for Z<28; crp- 28 aug, 1995 
!     EDITED 27 JULY 1994 BY GMW - REPLACED PT III PF COEFF. AND IP 
!     MODE 3 RETURNS PARTITION FUNCTION 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: DEBCON=1./2.8965E-18,                        & 
     &           TVCON=8.6171E-5,                                         & 
     &           HIONEV=13.595,                                           & 
     &           ONE=1.,                                                  & 
     &           HALF=0.5,                                                & 
     &           THIRD=1./3.,                                             & 
     &           X18=1./18.,                                              & 
     &           X120=1./120.,                                            & 
     &           T211=2000./11. 
! 
!     REAL(DP) :: F(6), 
      REAL(DP) :: IP(6),PART(6),POTLO(6) 
!     REAL(DP) :: FSAVE(6) 
      REAL(DP) :: SCALE(4) 
      INTEGER :: NNN(6*218) 
      INTEGER :: NNN16(54),NNN17(54),NNN18(54),NNN19(54),NNN20(54) 
      INTEGER :: NNN21(54),NNN22(54),NNN23(54),NNN24(54),NNN25(54) 
      INTEGER :: NNN26(54),NNN27(54),NNN28(54),NNN29(54),NNN30(54) 
      INTEGER :: NNN31(54),NNN32(54),NNN33(54),NNN34(54),NNN35(54) 
      INTEGER :: NNN36(54),NNN37(54),NNN38(54),NNN39(54),NNN40(12) 
 
      EQUIVALENCE (NNN( 811-810),NNN16(1)) 
      EQUIVALENCE (NNN( 865-810),NNN17(1)),(NNN( 919-810),NNN18(1)) 
      EQUIVALENCE (NNN( 973-810),NNN19(1)),(NNN(1027-810),NNN20(1)) 
      EQUIVALENCE (NNN(1081-810),NNN21(1)),(NNN(1135-810),NNN22(1)) 
      EQUIVALENCE (NNN(1189-810),NNN23(1)),(NNN(1243-810),NNN24(1)) 
      EQUIVALENCE (NNN(1297-810),NNN25(1)),(NNN(1351-810),NNN26(1)) 
      EQUIVALENCE (NNN(1405-810),NNN27(1)),(NNN(1459-810),NNN28(1)) 
      EQUIVALENCE (NNN(1513-810),NNN29(1)),(NNN(1567-810),NNN30(1)) 
      EQUIVALENCE (NNN(1621-810),NNN31(1)),(NNN(1675-810),NNN32(1)) 
      EQUIVALENCE (NNN(1729-810),NNN33(1)),(NNN(1783-810),NNN34(1)) 
      EQUIVALENCE (NNN(1837-810),NNN35(1)),(NNN(1891-810),NNN36(1)) 
      EQUIVALENCE (NNN(1945-810),NNN37(1)),(NNN(1999-810),NNN38(1)) 
      EQUIVALENCE (NNN(2053-810),NNN39(1)),(NNN(2107-810),NNN40(1)) 
!      ( 1)( 2)   ( 3)( 4)   ( 5)( 6)   ( 7)( 8)   ( 9)(10)   ( IP ) G 
      DATA NNN16/                                                         & 
     & 227027622, 306233052, 356839222, 446052912, 652382292,   763314,   & 
     & 108416342, 222428472, 353944332, 577378932, 110314303,  1814900,   & 
     & 198724282, 293236452, 468362702,  86511123, 136016073,  3516000,   & 
     & 279836622, 461857562, 720693022, 124915873, 192522633,  5600000,   & 
     & 262136422, 501167232,  87911303, 138916483, 190721673,  7900000,   & 
     & 201620781, 231026761, 314737361, 450555381, 692386911,   772301,   & 
     & 109415761, 247938311,  58910042, 190937022,  68311693,  2028903,   & 
     & 897195961, 107212972, 165021182, 260230862, 356940532,  3682900,   & 
     & 100010001, 100410231, 108712611, 167124841, 388460411,   939102/ 
      DATA NNN17/                                                         & 
     & 200020021, 201620761, 223726341, 351352061,  80812472,  1796001,   & 
     & 100610471, 122617301, 300566361, 149924112, 332342352,  3970000,   & 
     & 403245601, 493151431, 529654331, 559358091, 611065171,   600000,   & 
     &  99710051, 104511541, 135016501, 208226431, 321837921,  2050900,   & 
     & 199820071, 204521391, 229124761, 266028451, 302932131,  3070000,   & 
     & 502665261, 755183501, 901496201, 102410942, 117912812,   787900,   & 
     & 422848161, 512153401, 557458941, 636270361, 794489061,  1593000,   & 
     & 100010261, 114613921, 175221251, 249828711, 324436181,  3421000,   & 
     & 403143241, 491856701, 649173781, 840396751, 113013392,   981000/ 
      DATA NNN18/                                                         & 
     & 593676641, 884697521, 105911572, 129515012, 180322212,  1858700,   & 
     & 484470541,  91510972, 125614082, 157017612, 199722912,  2829900,   & 
     & 630172361, 799686381, 919797221, 102810942, 117712832,   975000,   & 
     & 438055511, 691582151,  94510732, 121413672, 152016732,  2150000,   & 
     & 651982921,  94610382, 113212492, 139515462, 169718482,  3200000,   & 
     & 437347431, 498951671, 538559501,  74710812, 169126672,  1183910,   & 
     & 705183611,  93510092, 111614162, 222932532, 427652992,  2160000,   & 
     & 510869921,  87410312, 123116552, 236530712, 377744832,  3590000,   & 
     & 100010001, 100010051, 105012781, 198535971,  65911422,  1399507/ 
      DATA NNN19/                                                         & 
     & 461049811, 522254261, 609088131, 168935052,  68612253,  2455908,   & 
     & 759990901, 101911142, 129017782, 302856642,  99414333,  3690000,   & 
     & 200020011, 200720361, 211523021, 269434141, 459163351,   417502,   & 
     & 100010001, 100110321, 129524961,  61014202, 291753192,  2750004,   & 
     & 473650891, 533156051,  66810932, 232950852,  99915303,  4000000,   & 
     & 100110041, 104111741, 146019721, 281941411, 607785251,   569202,   & 
     & 202621931, 255331271, 384347931, 624085761, 122417632,  1102600,   & 
     & 100010001, 100110321, 129524961,  61014202, 291753192,  4300000,   & 
     & 791587851, 100012192, 155119942, 254031782, 389946932,   637900/ 
      DATA NNN20/                                                         & 
     & 118217102, 220827002, 319036792, 416646512, 513256072,  1223000,   & 
     &  92510012, 104710862, 112311612, 120212472, 132814282,  2050000,   & 
     & 141320802, 291439702, 531170262,  92712273, 162521053,   684000,   & 
     & 354454352, 724689652, 107212643, 148517093, 193321573,  1312900,   & 
     & 209727032, 324537052, 415446282, 510255752, 604965222,  2298000,   & 
     & 256636022, 465759302, 749693962, 116514243, 171520333,   687900,   & 
     & 335157222,  84511463, 147718363, 221826083, 299933893,  1431900,   & 
     & 223725352, 280830972, 340937362, 406844002, 473150632,  2503900,   & 
     & 703972941,  82610822, 154822682, 327244912, 571469372,   709900/ 
      DATA NNN21/                                                         & 
     &  75714552, 274347322, 718897632, 123414913, 174920063,  1614900,   & 
     & 267645462, 669890262, 115514323, 173620673, 242528083,  2714900,   & 
     &  90613732, 184823562, 291735332, 419949102, 565764332,   728000,   & 
     & 131318312, 227126932, 311735452, 397644072, 483852692,  1525900,   & 
     & 204721673, 234725733, 284031463, 348738613, 426546943,  3000000,   & 
     & 176824122, 318941082, 515263202, 761790472, 106112303,   736400,   & 
     & 221934642, 501968372,  88911173, 136316243, 189221613,  1675900,   & 
     & 210622722, 241025422, 267928262, 297731272, 327834282,  2846000,   & 
     & 148520202, 255230902, 364942462, 489656082, 638872352,   746000/ 
      DATA NNN22/                                                         & 
     & 153421292, 288137912, 484660322, 720187062, 101011483,  1807000,   & 
     & 254537212, 492362292, 770592182, 107312243, 137615273,  3104900,   & 
     & 115919651, 320746011, 607576761,  95011642, 141817172,   832900,   & 
     & 755087211, 105913442, 173122222, 282034722, 412247732,  1941900,   & 
     & 180223462, 289735212, 414247632, 538460052, 662672472,  3292000,   & 
     & 200020001, 200220141, 206422141, 257633021, 455164681,   757403,   & 
     & 100810581, 125817401, 260641031,  66210072, 135316982,  2148000,   & 
     & 795887491,  97711762, 156620252, 248329422, 340038582,  3481900,   & 
     & 100010001, 100410241, 109212891, 176827421, 444268771,   899003/ 
      DATA NNN23/                                                         & 
     & 200020021, 201720921, 233329881, 451475371, 127520782,  1690301,   & 
     & 100310281, 114815371, 246138311, 519265531, 791492761,  3747000,   & 
     & 252431921, 368440461, 433746521, 512259221, 723389021,   578400,   & 
     & 100110071, 104611651, 146118581, 225426511, 304734431,  1886000,   & 
     & 200120111, 205021611, 243628031, 317035371, 390442701,  2802900,   & 
     & 232637101, 488058571, 669074381, 816189091,  97210632,   734200,   & 
     & 286335941, 408144471, 479351961, 571862901, 686274341,  1462700,   & 
     & 100010251, 114013811, 175321601, 256829751, 338337901,  3049000,   & 
     & 404043481, 494656811, 646772781, 813490751, 101411372,   863900/ 
      DATA NNN24/                                                         & 
     & 303147981, 618472951, 827392621, 103711702, 131214532,  1650000,   & 
     & 313037601, 429347901, 536260591, 689477591, 862494881,  2529900,   & 
     & 526258801, 657372351, 784284071, 897095741, 102711082,   900900,   & 
     & 440855541, 686481251,  93810792, 125414792, 176321132,  1860000,   & 
     & 349054751, 699883081,  96611302, 134216202, 197724212,  2800000,   & 
     & 405342041, 438645621, 475751071, 587974491, 102214572,  1045404,   & 
     & 568567471, 773485861,  94510362, 112712182, 130914002,  1909000,   & 
     & 514269581,  86910562, 130716652, 215327742, 351843662,  3200000,   & 
     & 100010001, 100010091, 109515351, 291060661, 119621482,  1212716/ 
      DATA NNN25/                                                         & 
     & 414844131, 465649111, 538464651,  87112232, 158019362,  2120000,   & 
     & 615475101, 867797531, 112213462, 157618062, 203622662,  3209900,   & 
     & 200020001, 201020501, 215623871, 283536181, 462756261,   389300,   & 
     & 100010001, 100310371, 119016501, 269146361,  77912412,  2510000,   & 
     & 424445601, 481750061, 516953311, 549356551, 581759791,  3500000,   & 
     & 101210791, 135119351, 282340571, 574580391, 111015062,   521002,   & 
     & 262638611, 504160621, 698579371,  91010692, 129115952,  1000000,   & 
     & 100010001, 100310351, 118416321, 264945521,  76512182,  3700000,   & 
     &  71111992, 172323592, 312540402, 510763182, 765791012,   558000/ 
      DATA NNN26/                                                         & 
     & 204529582, 383647882, 582469262, 807992692, 104911723,  1106000,   & 
     &  94712552, 148416582, 179819212, 203621522, 227424042,  1916900,   & 
     & 295959132, 103515693, 215527593, 335939413, 449650223,   565000,   & 
     &  79718153, 289639443, 495159253, 686877533, 863794813,  1085000,   & 
     & 298640242, 475053692, 596965912, 725379692, 872094692,  2008000,   & 
     & 460693672, 158523823, 327242303, 519661563, 709379783,   541900,   & 
     & 455480232, 114014653, 178521013, 240927073, 299232633,  1055000,   & 
     &  46410533, 183826893, 354443773, 518459633, 674375243,  2320000,   & 
     & 139623042, 364860002,  96114603, 209828633, 373446973,   549000/ 
      DATA NNN27/                                                         & 
     & 460493692, 158523823, 327142303, 519661563, 709279783,  1073000,   & 
     & 455480232, 114014653, 178521013, 240927073, 299232633,  2000000,   & 
     & 131720482, 280535692, 441254492, 676583972, 103412583,   555000,   & 
     & 139623042, 364860002,  96114603, 209828633, 373446973,  1089900,   & 
     & 460493682, 158523823, 327142303, 519661563, 709279783,  2000000,   & 
     &  92915672, 222431062, 444763802,  89612173, 159520253,   562900,   & 
     & 315059662,  97114563, 204627093, 342541693, 490556383,  1106900,   & 
     & 269037812, 520270372,  91111273, 133915483, 172719093,  2000000,   & 
     & 800080571, 851699301, 127617362, 240433032, 444958442,   568000/ 
      DATA NNN28/                                                         & 
     & 125416052, 211828182, 375549622, 644381732, 101112213,  1125000,   & 
     & 800080571, 851699301, 127617362, 240433032, 444958442,  2000000,   & 
     & 240432982, 427555202, 708489962, 112613853, 167319843,   615900,   & 
     & 534793262, 139219123, 247730843, 371043333, 495055893,  1210000,   & 
     & 364145232, 514756362, 604864112, 673870372, 732276072,  2000000,   & 
     & 480767202,  89011393, 144118243, 230028753, 354142883,   584900,   & 
     & 480767192,  89011393, 144118243, 230028753, 354142883,  1151900,   & 
     & 480767202,  89011393, 144118243, 230028753, 354142883,  2000000,   & 
     & 343147532, 645887152, 115314793, 183322063, 257729373,   593000/ 
      DATA NNN29/                                                         & 
     & 343147532, 645887142, 115314793, 183322063, 257729373,  1167000,   & 
     & 343147532, 645887142, 115314793, 183322063, 257729373,  2000000,   & 
     & 222635002, 542276772, 100312353, 145716713, 187020703,   602000,   & 
     & 222635002, 542276772, 100312353, 145716713, 187020703,  1180000,   & 
     & 222635002, 542276772, 100312353, 145716713, 187020703,  2000000,   & 
     & 133715382, 209130152, 429859382,  79410293, 129815983,   609900,   & 
     & 265934782, 497877532, 120517733, 245032063, 400448073,  1193000,   & 
     & 265934782, 497877532, 120517733, 245032063, 400448073,  2000000,   & 
     & 800381111,  87510702, 147621462, 310343462, 585475982,   618000/ 
      DATA NNN30/                                                         & 
     & 156718872, 279244452, 678196342, 128316243, 197823443,  1205000,   & 
     &  93517192, 364666132, 103414613, 192624193, 293334613,  2370000,   & 
     & 100010011, 101310651, 118613951, 169120661, 250629971,   625000,   & 
     & 200120901, 270345231,  81714042, 223533112, 461959862,  1217000,   & 
     & 100312561, 250851931,  91914182, 198626022, 323638692,  2000000,   & 
     & 514664441, 759086851,  99211442, 133315612, 182721252,   609900,   & 
     & 125924831, 438667801,  98714112, 199727872, 380850742,  1389900,   & 
     & 323948621, 661297271, 158626482, 426865032,  93712843,  1900000,   & 
     & 659294081, 128016962, 222528952, 372047062, 585171462,   700000/ 
      DATA NNN31/                                                         & 
     &  99117882, 274638812, 520867322,  84410313, 123314453,  1489900,   & 
     & 187427702, 343739872, 448049452, 539358282, 625266642,  2329900,   & 
     &  65210892, 171325762, 373552252, 705192012, 116414343,   787900,   & 
     & 192837842, 600784802, 111113823, 165419233, 218524383,  1620000,   & 
     &  99117872, 274638812, 520867312,  84410313, 123314453,  2400000,   & 
     & 398981651, 130019172, 273438022, 516168382,  88411163,   797900,   & 
     & 131429482, 523279952, 111414623, 183422233, 262130233,  1770000,   & 
     & 192837842, 600784792, 111113823, 165419233, 218524383,  2500000,   & 
     & 600963001,  75910412, 150121572, 301940972, 539168952,   787000/ 
      DATA NNN32/                                                         & 
     &  73710852, 190731262, 464964142,  83810503, 127315053,  1660000,   & 
     & 131429482, 523279952, 111414623, 183422233, 262130233,  2600000,   & 
     & 110815502, 216829732, 398752322, 672484682, 104612673,   850000,   & 
     & 168225972, 362046562, 566766422, 757484612,  93010103,  1700000,   & 
     &  73710852, 190731262, 464964142,  83810503, 127315053,  2700000,   & 
     & 129117892, 239430882, 388748292, 596173252,  89510843,   910000,   & 
     & 110815502, 216829732, 398752322, 672484682, 104612673,  2000000,   & 
     & 168225972, 362046562, 566766422, 757484612,  93010103,  2800000,   & 
     & 158918512, 207523002, 254328242, 316335762, 407246582,   900000/ 
      DATA NNN33/                                                         & 
     &  98115462, 224930742, 401150612, 623475412,  89910583,  1855900,   & 
     & 146323292, 354651802,  74810923, 161723953, 348749363,  3322700,   & 
     & 203222611, 265731251, 364042301, 494958601, 702084731,   922000,   & 
     & 120521331, 357753801,  75310062, 130516572, 206925452,  2050000,   & 
     & 651780821, 108814772, 195925252, 316338622, 460853882,  3000000,   & 
     & 100010001, 100110111, 105211851, 152122101, 341552811,  1043002,   & 
     & 200320211, 210023021, 268834231, 480472341, 111416912,  1875000,   & 
     & 104012871, 186129471, 458664151,  82410072, 119013732,  3420000,   & 
     & 200420711, 222424271, 265429161, 325637371, 442853911,   610500/ 
      DATA NNN34/                                                         & 
     & 100010021, 101910801, 121414641, 189525811, 358949721,  2041900,   & 
     & 200020311, 216624611, 296337451, 489064791,  85711212,  2979900,   & 
     & 103411711, 147819101, 244331781, 434862751,  93113762,   741404,   & 
     & 204122231, 248227841, 311535621, 429153941, 651976431,  1502800,   & 
     & 100210131, 106812201, 154522671, 381665951,  95512512,  3192900,   & 
     & 400140351, 416944121, 474851591, 564362181, 690477231,   728700,   & 
     & 106814451, 204427341, 350744811, 586879131, 108314772,  1667900,   & 
     & 205523051, 264830231, 345439921, 469156001, 675281671,  2555900,   & 
     & 500950661, 518153561, 559058941, 628968071, 748483501,   843000/ 
      DATA NNN35/                                                         & 
     & 443756241, 696282451,  95411012, 128615262, 182922012,  1900000,   & 
     & 336953201, 682481011,  93810882, 127915272, 184622442,  2700000,   & 
     & 402841621, 431544771, 463148311, 520059491, 734896851,   930000,   & 
     & 576168741, 788387631,  96910642, 116012552, 135014462,  2000000,   & 
     & 490265341, 812797201, 116614322, 179622692, 285035302,  2900000,   & 
     & 100010001, 100010031, 102311051, 133018071, 264539391,  1074500,   & 
     & 402841621, 431544771, 463148311, 520059491, 734996851,  2000000,   & 
     & 576168741, 788387631,  96910642, 116012552, 135014462,  3000000,   & 
     & 200020011, 201220591, 218124481, 296538611, 488859141,   400000/ 
      DATA NNN36/                                                         & 
     & 100010001, 100010031, 102311051, 133018071, 264539401,  2200000,   & 
     & 421645151, 477449611, 511852711, 542455761, 572958821,  3300000,   & 
     & 100010041, 105212131, 153220271, 270435641, 460258111,   527600,   & 
     & 201221791, 258131471, 381645781, 546365131, 777592781,  1014400,   & 
     & 100010001, 100010031, 102311051, 133018071, 264539391,  3400000,   & 
     & 510064491,  82710872, 142718412, 232328712, 348341572,   690000,   & 
     & 228951571,  88513232, 183324132, 305537492, 448152402,  1210000,   & 
     & 723989131, 103511752, 130814352, 155416652, 177018682,  2000000,   & 
     & 620099241, 162725772, 391457072,  80110833, 141818023,   600000/ 
      DATA NNN37/                                                         & 
     & 620099241, 162725772, 391457072,  80110833, 141818023,  1200000,   & 
     & 620099251, 162725772, 391457072,  80110833, 141818023,  2000000,   & 
     & 347877992, 129318323, 240730533, 380546863, 570368573,   600000,   & 
     & 347877992, 129318323, 240730533, 380546863, 570368573,  1200000,   & 
     & 347777992, 129318323, 240730533, 380546863, 570368573,  2000000,   & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000,   & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,   & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000,   & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000/ 
      DATA NNN38/                                                         & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,   & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000,   & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000,   & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,   & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000,   & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000,   & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,   & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000,   & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000/ 
      DATA NNN39/                                                         & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,   & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000,   & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000,   & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,   & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000,   & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000,   & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,   & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000,   & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000/ 
      DATA NNN40/                                                         & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,   & 
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000/ 
      DATA SCALE/.001,.01,.1,1./ 
! 
      if(mode.lt.0) return 
      tk=1.38054e-16*t 
      tv=8.6171e-5*t 
!     LOWERING OF THE IONIZATION POTENTIAL IN VOLTS FOR UNIT ZEFF 
      CHARGE=ANE*2. 
      DEBYE=SQRT(TK*DEBCON/CHARGE) 
!     DEBYE=SQRT(TK/12.5664/4.801E-10**2/CHARGE) 
      POTLOW=MIN(1.,1.44E-7/DEBYE) 
      IF(IIZ.LE.28)then 
         write(6,*) 'Error, routine PFHEAV for Z.GE.28 only' 
         stop 
       endif 
!    removed elements with z<28 
      if(iiz.eq.28) n=1 
      IF(IIZ.GT.28) N=3*IIZ+54-135 
      IF(IIZ.eq.28) NIONS=4 
      IF(IIZ.GT.28) NIONS=3 
      NION2=MIN0(JNION+2,NIONS) 
      N=N-1 
! 
      DO ION=1,NION2 
         Z=ION 
         POTLO(ION)=POTLOW*Z 
         N=N+1 
         nnn6n=nnn(6+6*(N-1)) 
         NNN100=NNN6N/100 
         XN1= NNN100 
         IP(ION)=XN1*1.e-3 
         IG=NNN6N-NNN100*100 
         GGG=IG 
         T2000=IP(ION)*T211 
         IT=MAX0(1,MIN0(9, INT(T/T2000-HALF))) 
         XIT=IT 
         DT=T/T2000-XIT-HALF 
         PMIN=ONE 
         I=(IT+1)/2 
         nnnin=nnn(i+6*(N-1)) 
         K1=NNNIN/100000 
         K2=NNNIN-K1*100000 
         K3=K2/10 
         xk1=k1 
         xk3=k3 
         KSCALE=K2-K3*10 
         P1=XK1*SCALE(KSCALE) 
         P2=XK3*SCALE(KSCALE) 
         KP1=int(P1) 
         PMIN=KP1 
         IF(MOD(IT,2).EQ.0.OR.DT.LT.0..AND.KSCALE.LE.1.AND.               & 
     &      KP1.EQ.INT(P2+.5)) THEN 
            PMIN=KP1 
            xk3=k3 
            P1=XK3*SCALE(KSCALE) 
            nnni1n=nnn(i+1+6*(N-1)) 
            K1=NNNI1N/100000 
            KSCALE=MOD(NNNI1N,10) 
            xk1=k1 
            P2=XK1*SCALE(KSCALE) 
         END IF 
         PART(ION)= MAX (PMIN,P1+(P2-P1)*DT) 
         IF(GGG.NE.0..AND.POTLO(ION).GE..1.AND.T.GE.T2000*4.) THEN 
            IF(T.GT.(T2000*11.)) TV=(T2000*11.)*TVCON 
            D1=.1/TV 
            D2=POTLO(ION)/TV 
            DX=SQRT(HIONEV*Z*Z/TV/D2)**3 
            PART(ION)=PART(ION)+GGG*EXP(-IP(ION)/TV)*                     & 
     &              (DX*(THIRD+(ONE-(HALF+(X18+D2*X120)*D2)*D2)*D2)-      & 
     &              DX*(THIRD+(ONE-(HALF+(X18+D1*X120)*D1)*D1)*D1)) 
         END IF 
      END DO 
      u=part(jnion) 
      RETURN 
      END SUBROUTINE PFHEAV 
! 
! ****************************************************************** 
! 
      subroutine frac1 
!     ================ 
! 
      use accura 
      use params 
      use modelp 
      use opadat 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      real(dp) :: xxt(mdepth),xxe(mdepth) 
      integer  :: kt0(mdepth),kn0(mdepth) 
! 
      do id=1,nd 
         xxt(id)=dlog10(temp(id)) 
         kt0(id)=2*int(20.*xxt(id)) 
         xxe(id)=dlog10(elec(id)) 
         kn0(id)=int(2.*xxe(id)) 
      end do 
! 
      ATOMS: DO IAT=1,30 
         iatnum=iat 
         call fractn(iatnum) 
         if(iatnum.le.0) cycle atoms 
         depths: do id=1,nd 
           if(kt0(id).lt.itemp(1)) then 
              kt1=1 
              write(6,"(' (FRACOP) Extrapol. in T (low)',i4,f7.0)")       & 
     &          id,temp(id) 
            else if(kt0(id).ge.itemp(ntt)) then 
              kt1=ntt-1 
              write(6,"(' (FRACOP) Extrapol. in T (high)',i4,f12.0)")     & 
     &          id,temp(id) 
            else 
             do it=1,ntt 
                if(kt0(id).eq.itemp(it)) then 
                   kt1=it 
                   exit 
                end if 
             end do 
           end if 
 
           if(kn0(id).lt.1) then 
              kn1=1 
            else if(kn0(id).ge.60) then 
              kn1=59 
              write(6,"(' (FRACOP) Extrapol. in Ne (high)',i4,f9.4)")     & 
     &           id,xxe(id) 
            else 
              kn1=kn0(id) 
           end if 
 
           xt1=0.025*itemp(kt1) 
           dxt=0.05 
           at1=(xxt(id)-xt1)/dxt 
           xn1=0.5*kn1 
           dxn=0.5 
           an1=(xxe(id)-xn1)/dxn 
           do ion=1,mion1 
              x11=frac(kt1,kn1,ion) 
              x21=frac(kt1+1,kn1,ion) 
              x12=frac(kt1,kn1+1,ion) 
              x22=frac(kt1+1,kn1+1,ion) 
              x1221=x11*x21*x12*x22 
              if(x1221.eq.0.) then 
                  xx1=x11+at1*(x21-x11) 
                  xx2=x12+at1*(x22-x12) 
                  rrx=xx1+an1*(xx2-xx1) 
              else 
                  x11=dlog10(x11) 
                  x21=dlog10(x21) 
                  x12=dlog10(x12) 
                  x22=dlog10(x22) 
                  xx1=x11+at1*(x21-x11) 
                  xx2=x12+at1*(x22-x12) 
                  rrx=xx1+an1*(xx2-xx1) 
                  rrx=exp(2.3025851*rrx) 
              endif 
              rrr(id,ion,iat)=rrx*abndd(iat,id)*                          & 
     &                        dens(id)/wmm(id)/ytot(id) 
           end do 
         end do depths 
      END DO ATOMS 
! 
      return 
      end subroutine frac1 
! 
! ****************************************************************** 
! 
      subroutine fractn(iatnum) 
!     ========================= 
! 
      use opadat 
      implicit real(dp) (a-h,o-z) 
 
      INTEGER, PARAMETER :: MDAT = 17, INP=71 
      integer  :: ioo(-1:mion1),idat(mion1) 
      real(dp) :: frac0(-1:mion1) 
      real(dp) :: gg(mion1,mdat),g0(mion1),z0(-1:mion1) 
      real(dp) :: uu(mion1,mdat),u0(mion1) 
      real(dp) :: u6(6),u7(7),u8(8),u10(10),u11(11) 
      real(dp) :: u12(12),u13(13),u14(14),u16(16),u18(18),u20(20) 
      real(dp) :: u24(24),u25(25),u26(26),u28(28) 
 
      equivalence (u6(1),uu(1,3)),(u7(1),uu(1,4)),(u8(1),uu(1,5)) 
      equivalence (u10(1),uu(1,6)),(u11(1),uu(1,7)),(u12(1),uu(1,8)) 
      equivalence (u13(1),uu(1,9)),(u14(1),uu(1,10)),(u16(1),uu(1,11)) 
      equivalence (u18(1),uu(1,12)),(u20(1),uu(1,13)),(u24(1),uu(1,14)) 
      equivalence (u25(1),uu(1,15)),(u26(1),uu(1,16)),(u28(1),uu(1,17)) 
      data idat   / 1, 2, 0, 0, 0, 3, 4, 5, 0, 6,                         & 
     &              7, 8, 9,10, 0,11, 0,12, 0,13,                         & 
     &              0, 0, 0,14,15,16, 0,17, 0, 0/ 
      data gg/2.,29*0.,2.,1.,28*0.,                                       & 
     &        2.,1.,2.,1.,6.,9.,24*0.,2.,1.,2.,1.,6.,9.,4.,23*0.,         & 
     &        2.,1.,2.,1.,6.,9.,4.,9.,22*0.,                              & 
     &        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,20*0.,                        & 
     &        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,19*0.,                     & 
     &        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,18*0.,                  & 
     &        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,17*0.,               & 
     &        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,16*0.,            & 
     &        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,14*0.,      & 
     &        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,6.,1.,      & 
     &        12*0.,2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,      & 
     &        6.,1.,2.,1.,10*0.,2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,      & 
     &        6.,9.,4.,9.,6.,1.,10.,21.,28.,25.,6.,7.,6*0.,               & 
     &        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,            & 
     &           6.,1.,10.,21.,28.,25.,6.,7.,6.,5*0.,                     & 
     &        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,            & 
     &           6.,1.,10.,21.,28.,25.,6.,25.,30.,25.,4*0.,               & 
     &        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,            & 
     &           6.,1.,10.,21.,28.,25.,6.,25.,28.,21.,10.,21.,0.,0./ 
      data uu(1,1)/109.6787/ 
      data uu(1,2)/198.3108/ 
      data uu(2,2)/438.9089/ 
      data u6/90.82,196.665,386.241,520.178,3162.395,3952.061/ 
      data u7/117.225,238.751,382.704,624.866,789.537,4452.758,5380.089/ 
      data u8/109.837,283.24,443.086,624.384,918.657,1114.008,5963.135,   & 
     &        7028.393/ 
      data u10/173.93,330.391,511.8,783.3,1018.,1273.8,1671.792,          & 
     &         1928.462,9645.005,10986.876/ 
      data u11/41.449,381.395,577.8,797.8,1116.2,1388.5,1681.5,2130.8,    & 
     &         2418.7,11817.061,13297.676/ 
      data u12/61.671,121.268,646.41,881.1,1139.4,1504.3,1814.3,2144.7,   & 
     &         2645.2,2964.4,14210.261,15829.951/ 
      data u13/48.278,151.86,229.446,967.8,1239.8,1536.3,1947.3,2295.4,   & 
     &         2663.4,3214.8,3565.6,16825.022,18584.138/ 
      data u14/65.748,131.838,270.139,364.093,1345.1,1653.9,1988.4,       & 
     &         2445.3,2831.9,3237.8,3839.8,4222.4,19661.693,21560.63/ 
      data u16/83.558,188.2,280.9,381.541,586.2,710.184,2265.9,2647.4,    & 
     &         3057.7,3606.1,4071.4,4554.3,5255.9,5703.6,26002.663,       & 
     &         28182.535/ 
      data u18/127.11,222.848,328.6,482.4,605.1,734.04,1002.73,1157.08,   & 
     &         3407.3,3860.9,4347.,4986.6,5533.8,6095.5,6894.2,7404.4,    & 
     &         33237.173,35699.936/ 
      data u20/49.306,95.752,410.642,542.6,681.6,877.4,1026.,1187.6,      & 
     &         1520.64,1704.047,4774.,5301.,5861.,6595.,7215.,7860.,      & 
     &         8770.,9338.,41366.,44177.41/ 
      data u24/54.576,132.966,249.7,396.5,560.2,731.02,1291.9,1490.,      & 
     &         1688.,1971.,2184.,2404.,2862.,3098.52,8151.,8850.,         & 
     &         9560.,10480.,11260.,12070.,13180.,13882.,60344.,63675.9/ 
      data u25/59.959,126.145,271.55,413.,584.,771.1,961.44,1569.,        & 
     &         1789.,2003.,2307.,2536.,2771.,3250.,3509.82,9152.,         & 
     &         9872.,10620.,11590.,12410.,13260.,14420.,15162.,           & 
     &         65660.,69137.4/ 
      data u26/63.737,130.563,247.22,442.,605.,799.,1008.,1218.38,        & 
     &         1884.,2114.,2341.,2668.,2912.,3163.,3686.,3946.82,         & 
     &         10180.,10985.,11850.,12708.,13620.,14510.,15797.,          & 
     &         16500.,71203.,74829.6/ 
      data u28/61.6,146.542,283.8,443.,613.5,870.,1070.,1310.,1560.,      & 
     &         1812.,2589.,2840.,3100.,3470.,3740.,4020.,4606.,           & 
     &         4896.2,12430.,13290.,14160.,15280.,16220.,17190.,          & 
     &         18510.,19351.,82984.,86909.4/ 
! 
      if(idat(iatnum).eq.0) then 
         write(6,"(' OP data for element no. ',i3,' do not exist')")      & 
     &   iatnum 
         iatnum=-1 
         return 
      end if 
! 
      g0(iatnum+1)=1. 
      do i=1,iatnum 
        ig0=iatnum-i+1 
        g0(ig0)=gg(i,idat(iatnum)) 
        u0(i)=uu(i,idat(iatnum))*1000. 
      enddo 
! 
      if(iatnum.eq.1) open(inp,file='ioniz.dat',status='old') 
      do it=1,mtemp 
         do ie=1,melec 
            fracm(it,ie)=0. 
            do ion=1,mion1 
               frac(it,ie,ion)=0. 
            end do 
         end do 
      end do 
! 
      read(inp,*) 
      read(inp,*) it0,it1,itstp 
      ntt=(it1-it0)/itstp+1 
! 
      do it=1,ntt 
         read(inp,*) itt,ie0,ie1,iestp 
         itemp(it)=itt 
         net=(ie1-ie0)/iestp+1 
         t=exp(2.3025851*0.025*itt) 
         safac0=sqrt(t)*t/2.07e-16 
         tkcm=0.69496*t 
         do ie=1,net 
            read(inp,"(3i4,2x,4(i4,1x,e9.3))") iee,ion0,ion1,             & 
     &                (ioo(i),frac0(i),i=ion0,min(ion1,ion0+3)) 
            ane=exp(2.3025851*0.25*iee) 
            safac=safac0/ane 
            nio=ion1-ion0 
            if(nio.ge.3) then 
               nlin=nio/4 
               do ilin=1,nlin 
                  read(inp,"(14x,4(i4,1x,e9.3))") (ioo(i),frac0(i),       & 
     &                 i=ion0+4*ilin,min(ion1,ion0+4*ilin+3)) 
               end do 
            end if 
            ieind=iee/2 
            do ion=ion0,ion1 
              if(ion.lt.iatnum) then 
               if(ion.eq.ion0) then 
                  z0(ion)=g0(iatnum-ion) 
               else 
                  z0(ion)=frac0(ion)/frac0(ion-1)*safac*z0(ion-1) 
                  z0(ion)=z0(ion)*exp(-u0(iatnum-ion)/tkcm) 
               endif 
                  frac(it,ieind,iatnum-ion)=frac0(ion)/z0(ion) 
              else 
                  u0hm=6090.5 
                  z0hm=frac0(ion)/frac0(ion-1)*safac 
                  z0hm=z0hm*exp(-u0hm/tkcm) 
                  fracm(it,ieind)=frac0(ion)/z0hm 
              end if 
            end do 
         end do 
      end do 
      return 
      end subroutine fractn 
! 
! 
!     ******************************************************************* 
! 
! 
 
      SUBROUTINE DWNFR0(ID) 
!     ===================== 
! 
!     Auxiliary quantities for dissolved fractions 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: UN=1.,SIXTH=UN/6.,CCOR=0.0 
      real(dp), parameter :: p1=0.1402,p2=0.1285,p3=un,p4=3.15,p5=4. 
      real(dp), parameter :: f23=-2./3. 
! 
      ANE=ELEC(ID) 
      ELEC23(ID)=EXP(F23*LOG(ANE)) 
      ANES=EXP(SIXTH*LOG(ANE)) 
      ACOR=CCOR*ANES/SQRT(TEMP(ID)) 
      X=EXP(P4*LOG(UN+P3*ACOR)) 
      DWC2(ID)=P2*X 
      A3=ACOR*ACOR*ACOR 
      DO IZZ=1,MZZ 
         Z3(IZZ)=IZZ*IZZ*IZZ 
         DWC1(IZZ,ID)=P1*(X+P5*(IZZ-1.)*A3) 
      END DO 
      RETURN 
      END SUBROUTINE DWNFR0 
! 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE DWNFR1(FR,FR0,ID,IZZ,DW1) 
!     ==================================== 
! 
!     dissolved fraction for frequency FR 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: UN=1.,TKN=3.01,CKN=5.33333333,CB=8.59e14 
      REAL(DP), PARAMETER :: SQFRH=5.734152e7,a0=0.529177e-8 
      real(dp), parameter :: wa0=-3.1415926538/6.*a0*a0*a0 
! 
      IF(FR.LT.FR0) THEN 
         XN=SQFRH*IZZ/SQRT(FR0-FR) 
         if(xn.le.tkn) then 
            xkn=un 
          else 
            xn1=un/(xn+un) 
            xkn=ckn*xn*xn1*xn1 
         end if 
         BETA=CB*Z3(IZZ)*XKN/(XN*XN*XN*XN)*ELEC23(ID) 
         beta=beta*bergfc 
         BETA3=BETA*BETA*BETA 
         BETA32=SQRT(BETA3) 
         F=(DWC1(IZZ,ID)*BETA3)/(UN+DWC2(ID)*BETA32) 
! 
!     contribution from neutral particles 
! 
         xn2=xn*xn+un 
         xnh=0. 
         xnhe1=0. 
         if(ielh.gt.0) xnh=popul(nfirst(ielh),id) 
         if(ielhe1.gt.0) xnhe1=popul(nfirst(ielhe1),id) 
         w0=exp(wa0*xn2*xn2*xn2*(xnh+xnhe1)) 
         W0=1. 
! 
         DW1=UN-F/(UN+F)*w0 
       ELSE 
         DW1=UN 
      END IF 
      RETURN 
      END SUBROUTINE DWNFR1 
! 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE CHCKAB 
!     ================= 
! 
!     check input abumdances of explicit atoms (unit 5) and those 
!     which follow from the models atmosphere (unit 7) obtained by 
!     summing all populations and upper sums 
!     The program stops if it finds discrepancy more than 10 % 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      real(dp) :: sumpop(matom),sumiat(matom) 
! 
      IST=0 
      DO ID1=1,3 
         IF(ID1.EQ.1) ID=1 
         IF(ID1.EQ.2) ID=46 
         IF(ID1.EQ.3) ID=ND 
         CALL WNSTOR(ID) 
         ANE=ELEC(ID) 
         CALL SABOLF(ID) 
         DO IAT=1,NATOM 
            SUM=0. 
            sump=0. 
            DO I=N0A(IAT),NKA(IAT) 
               IL=ILK(I) 
               A=1. 
               IF(IL.GT.0) A=1.+ANE*USUM(IL) 
               SUM=SUM+A*POPUL(I,ID) 
               SUMP=SUMP+POPUL(I,ID) 
            END DO 
            SUMIAT(IAT)=SUM 
            SUMPOP(IAT)=SUMP 
         END DO 
         WRITE(6,"(' check of abundances (id =',i3/                       & 
     &   ' computed from model atmosphere  - input abundances'/)") ID 
         DO IAT=1,NATOM 
            X=SUMIAT(IAT)/SUMIAT(IATREF) 
            WRITE(6,"(i5,1p3e20.3)")                                      & 
     &      IAT,X,abund(iat,id),SUMPOP(IAT)/SUMPOP(IATREF) 
            IF(X/abund(iat,id).GT.1.1.OR.X/abund(iat,id).LT.0.9)          & 
     &         ist=ist+1 
         END DO 
      END DO 
      IF(IST.GT.0) THEN 
         WRITE(6,"(' ERROR !!! - inconsistent abundances'/)") 
         STOP 
      END IF 
      RETURN 
      END SUBROUTINE CHCKAB 
! 
! 
! ******************************************************************** 
! 
! 
      subroutine molini 
!     ================= 
! 
!     Initialization of the molecular equilibrium 
! 
      use accura 
      use params 
      use modelp 
      use eospar 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      real(dp) :: hpo(mdepth) 
! 
      aeinit=1.0 
! 
      do id=1,nd 
         t=temp(id) 
         tln=log(t)*1.5 
         thl=11605./t 
         t32=exp(tln) 
 
         do i=1,MMOLEC 
            rrmol(i,id)=0. 
         end do 
 
         hpo(id)=DENS(ID)/WMM(ID)/YTOT(ID) 
 
         if(t.gt.tmolim) cycle 
         HPOP=DENS(ID)/WMM(ID)/YTOT(ID) 
 
         an=dens(id)/wmm(id)+elec(id) 
         aeinit=0.1*an 
         if(t.lt.4000.) aeinit=0.01*an 
         call moleq(id,t,an,aeinit,ane,0) 
! 
!     next initial guess will be the last ane determined for 
!     previous depth point 
 
         aeinit=ane 
! 
         if (id.eq.idstd.and.ifeos.eq.0) then 
            write(6,                                                      & 
     &      "(/ 'Molecular number densities at the standard depth'/)") 
            nmol=nmolec 
            if(id.eq.1) nmol=32 
            do i=1,nmol 
               write(6,"(i4,1x,A8,1x,1pe12.2,1x,e12.2)")                  & 
     &         i, cmol(i), rrmol(i,id), pfmol(i,id)
            end do 
         end if 
      end do 
 
!     update atomic populations once molecular densities are calculated 
 
      if(imode.lt.-4) then 
      do i=1,nlevel 
         iat=numat(iatm(i)) 
         ion=iz(iel(i)) 
         ii=nfirst(iel(i)) 
         ener=(enion(ii)-enion(i))/bolk 
         if((enion(i).eq.0).and.(ilk(i).gt.0)) then 
            ener=0. 
            ion=ion+1 
         end if 
         if(ifwop(i).ge.0) then 
            do id=1,nd 
               popul(i,id)=rrr(id,ion,iat)*g(i)                           & 
     &              *exp(-ener/temp(id)) 
               if(iat.eq.1.and.ion.eq.0) popul(i,id)=anhm(id) 
            end do 
         endif 
      end do 
      end if 
! 
      return 
      end subroutine molini 
! 
! 
! ******************************************************************** 
! 
! 
      subroutine lists 
!     ================ 
! 
      use accura 
      use molist 
      use lindat 
 
      implicit real(dp) (a-h,o-z) 
      character(len=80) :: dum 
      character(len=40) :: am 
      character(len=6)  :: ilab 
      character(len=2)  :: iu 
      real(dp)          :: x(9) 
 
      ilist=0 
 
      iunit=20 
      istp=10 
      ibin(0)=mod(inlist,istp) 
! 
!     ---------------------------------- 
!     old reading of line lists from fort.55; 
!     but if fort.3 is present, it rewrites information from fort.55 
!     ---------------------------------- 
! 
      read(55,*,iostat=ios)                                               & 
     &   nmlist,(iunitm(ilist),ilist=1,nmlist) 
      if(ios.eq.0) then 
         do ilist=1,nmlist 
            write(iu,"(i2)") iunitm(ilist) 
            amlist(ilist) ='fort.' // iu 
            ibin(ilist)=ibin(0) 
            tmlim(ilist)=tmolim 
!           write(*,*) 'ilist',ilist,amlist(ilist),ibin(ilist) 
         end do 
       else 
         nmlist=nmlis0 
         iunitm(1)=iunim1 
         iunitm(2)=iunim2 
      end if 
 
!     ---------------------------------- 
!     reading line list data from unit 3 
!     ---------------------------------- 
! 
      ilist=0 
      amlist(0)='fort.19' 
      read(3,*,iostat=ios) am,ibi 
      if(ios.eq.0) then 
         amlist(0)=am 
         ibin(0)=ibi 
         do 
            read(3,*,iostat=ios) am,ib,tml 
            if(ios.ne.0) exit 
            ilist=ilist+1 
            amlist(ilist)=am 
            ibin(ilist)=ib 
            tmlim(ilist)=tml 
         end do 
         nmlist=ilist 
         if(nmlist.gt.0.and.ifmol.eq.0) then 
            write(*,*) 'NEEDS TO SET IFMOL > 0 with NMLIST>0' 
            stop 
         end if 
      end if 
! 
      ilist=0 
      ilab='ATOMIC' 
      write(6,"((/'************************'/                             & 
     &        ' LINE LISTS:'/                                             & 
     &       /' ILIST  TYPE',9x,'FILENAME     IBIN    TMLIM'/             & 
     &        ' ----------------------------------------------'/          & 
     &        i4,4x,a6,2x,a15,2x,i4,f11.1))")                             & 
     &ilist,ilab,trim(amlist(ilist)),ibin(ilist) 
      ilab='MOLEC ' 
      do ilist=1,nmlist 
         write(6,"(i4,4x,a6,2x,a15,2x,i4,f11.1)")                         & 
     &   ilist,ilab,trim(amlist(ilist)),ibin(ilist),tmlim(ilist) 
      end do 
! 
!     ------------------------------- 
!     detect the type of the line list 
! 
      nvdwli=0 
      nbroad=0 
      do ilist=1,nmlist 
         ivdwli(ilist)=0 
         ibroli(ilist)=0 
 
         if(ibin(ilist).eq.0) then 
            open(unit=iunit,file=amlist(ilist),status='old') 
          else 
            open(unit=iunit,file=amlist(ilist),form='unformatted',        & 
     &        status='old') 
         end if 
! 
!        text list 
! 
         if(ibin(ilist).eq.0) then 
            read(iunit,'(a80)') dum 
            read(dum,*,iostat=kst1) (x(i),i=1,9) 
            if(kst1.eq.0) then 
               np=9 
               nvdwli=nvdwli+1 
               ivdwli(ilist)=nvdwli 
             else 
               read(dum,*,iostat=kst2) (x(i),i=1,7) 
               if(kst2.eq.0) then 
                  np=7 
                  nbroad=nbroad+1 
                  ibroli(ilist)=nbroad 
                else 
                  read(dum,*,iostat=kst3) (x(i),i=1,4) 
                  if(kst3.eq.0) then 
                     np=4 
                   else 
                     write(*,*) 'no applicable format of line list',      & 
     &                          ilist 
                  end if 
               end if 
            end if 
          else 
! 
!         binary list 
! 
            read(iunit,iostat=ios) (x(i),i=1,9) 
            if(ios.eq.0) then 
               np=9 
               nvdwli=nvdwli+1 
               ivdwli(ilist)=nvdwli 
             else 
               read(iunit,iostat=ios1) (x(i),i=1,7) 
               if(ios1.eq.0) then 
                  np=7 
                  nbroad=nbroad+1 
                  ibroli(ilist)=nbroad 
                else 
                  read(iunit,iostat=ios2) (x(i),i=1,4) 
                  if(ios2.eq.0) then 
                     np=4 
                  end if 
               end if 
            end if 
         end if 
         nmpar(ilist)=np 
         rewind(iunit) 
      end do 
 
      MMLIST=NMLIST 
      MBROAD=NBROAD 
      MVDWLI=NVDWlI 
      if(nmlist.gt.0) write(6,"(/                                         & 
     & 'TOTAL NUMBER OF MOLECULAR LINE LISTS        :',i4/                & 
     & 'NUMBER OF LISTS WITH 3 BROADENING PARAMETERS:',I4/                & 
     & 'NUMBER OF LISTS WITH EXOMOL VDW BROADENING  :',I4/)")             & 
     &  nmlist,nbroad,nvdwli 
 
      call alloc_molist 
 
      end subroutine lists 
! 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE INMOLI(ILIST) 
!     ======================== 
! 
!     read in the input molecular line list, 
!     selection of lines that may contribute, 
!     set up auxiliary fields containing line parameters, 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use lindat 
      use molist 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER ::                                              & 
     &           PI4=7.95774715E-2,                                       & 
     &           C1     = 2.3025851,                                      & 
     &           C2     = 4.2014672,                                      & 
     &           C3     = 1.4387886,                                      & 
     &           CNM    = 2.997925e17,                                    & 
     &           EXT0   = 3.17,                                           & 
     &           UN     = 1.0,                                            & 
     &           TEN    = 10.,                                            & 
     &           HUND   = 1.e2,                                           & 
     &           TENM4  = 1.e-4,                                          & 
     &           TENM8  = 1.e-8,                                          & 
     &           OP4    = 0.4,                                            & 
     &           AGR0=2.4734E-22,                                         & 
     &           XEH=13.595, XET=8067.6, XNF=25.,                         & 
     &           R02=2.5, R12=45., VW0=4.5E-9 
 
      save almm00 
 
      if(imode.ne.-3.and.temp(idstd).gt.tmolim) return 
      IUNIT=20 
 
      if(ibin(ilist).eq.0) then 
         open(unit=iunit,file=amlist(ilist),status='old') 
       else 
         open(unit=iunit,file=amlist(ilist),form='unformatted',           & 
     &        status='old') 
      end if 
! 
!     define a conversion table between Kurucz notation and Tsuji table 
!     through array MOLIND 
! 
      do i=1,11000 
         molind(i)=0 
      end do 
      molind(101)=2 
      molind(106)=5 
      molind(107)=12 
      molind(108)=4 
      molind(109)=33 
      molind(111)=122 
      molind(112)=32 
      molind(113)=133 
      molind(114)=17 
      molind(116)=16 
      molind(117)=36 
      molind(120)=34 
      molind(124)=198 
      molind(126)=214 
      molind(606)=8 
      molind(607)=7 
      molind(608)=6 
      molind(614)=21 
      molind(616)=20 
      molind(707)=9 
      molind(708)=11 
      molind(714)=24 
      molind(716)=23 
      molind(808)=10 
      molind(812)=126 
      molind(813)=134 
      molind(814)=25 
      molind(816)=26 
      molind(820)=179 
      molind(822)=29 
      molind(823)=30 
      molind(839)=246 
      molind(840)=31 
      molind(10108)=3 
! 
!     analogous indices for positive ions 
! 
      do i=1,11000 
         ionind(i)=0 
      end do 
      ionind(101)=427 
      ionind(106)=437 
      ionind(107)=438 
      ionind(108)=439 
      ionind(109)=440 
      ionind(111)=0 
      ionind(112)=442 
      ionind(113)=443 
      ionind(114)=444 
      ionind(116)=446 
      ionind(117)=447 
      ionind(120)=0 
      ionind(124)=0 
      ionind(126)=0 
      ionind(606)=429 
      ionind(607)=452 
      ionind(608)=453 
      ionind(614)=0 
      ionind(616)=0 
      ionind(707)=430 
      ionind(708)=454 
      ionind(714)=0 
      ionind(716)=455 
      ionind(808)=431 
      ionind(812)=0 
      ionind(813)=0 
      ionind(814)=457 
      ionind(816)=459 
      ionind(820)=0 
      ionind(822)=0 
      ionind(823)=0 
      ionind(839)=0 
      ionind(840)=0 
      ionind(10108)=0 
! 
!     improved association of Kurucz lable to Tsuji indices 
! 
      call molindx 
! 
!     setting auxiliary parameters for handling the line lists 
! 
 
      ALAST=CNM/FRLAST 
      ALASTM(ILIST)=ALAST 
      IL=0 
!!    write(6,"(2i5,f12.3,1pe12.3)") 
!!   * ilist,nxtsem(ilist),alastm(ILIST),frlasm(ilist) 
      IF(NXTSEM(ILIST).EQ.1) THEN 
          ALAM0=ALM00 
          ALASTM(ILIST)=ALST00 
          FRLASM(ILIST)=CNM/ALASTM(ILIST) 
          NXTSEM(ILIST)=0 
          REWIND IUNIT 
      END IF 
      ALMM00=ALAM0 
      DOPSTD=1.E7/ALAM0*DSTD 
      DOPLAM=ALAM0*ALAM0/CNM*DOPSTD 
      AVAB=ABSTD(IDSTD)*RELOP 
      ASTD=1.0 
      CUTOFF=CUTOF0 
      ALAST=CNM/FRLAST 
!!    write(6,"(2i5,f12.3,1pe12.3)") 
!!   * ilist,nxtsem(ilist),alastm(ILIST),frlasm(ilist) 
! 
!     *********************************** 
! 
!     first part of reading line list - read only lambda, and 
!     skip all lines with wavelength below ALAM0-CUTOFF 
! 
      REWIND IUNIT 
      ALAM=0. 
      IJC=2 
 
!     write(6,*) 
      write(*,*) '----------------------------------------------------' 
      write(*,*) 
! 
      READ1: DO 
         IF(IBIN(ILIST).EQ.0) THEN 
            READ(IUNIT,*,IOSTAT=IOS) ALAM 
          ELSE 
            READ(IUNIT,IOSTAT=IOS) ALAM 
         END IF 
         IF(IOS.NE.0.OR.ALAM.GT.ALAST+CUTOFF) THEN 
            INACTM(ILIST)=1 
            WRITE(6,"(25x,'LIST',I3,'  FILE: ',A,                         & 
     &      '  DECLARED INACTIVE')") ILIST,TRIM(AMLIST(ILIST)) 
            WRITE(6,"('because its first wavelength',f10.3,               & 
     &      '  is larger than the endpoint',f10.3,' A'/)")                & 
     &      alam*10.,(alast+cutoff)*10. 
 
            RETURN 
          ELSE IF(ALAM.LT.ALAM0-CUTOFF) THEN 
            CYCLE READ1 
          ELSE 
            EXIT READ1 
         END IF 
      END DO READ1 
      BACKSPACE(IUNIT) 
! 
!     *********************************** 
! 
!     read the line list 
! 
      ill=0 
      READLINES: DO 
         ill=ill+1 
         np=nmpar(ilist) 
         if(ibin(ilist).eq.0) then 
            if(np.eq.9) then 
               read(iunit,*,iostat=ios1) alam,anum,gf,excl,gr,gh2,xnh2,   & 
     &                                ghe,xnhe 
               if(ios1.ne.0) exit readlines 
             else if(np.eq.7) then 
               READ(IUNIT,*,iostat=ios2) ALAM,ANUM,GF,EXCL,GR,GS,GW 
               if(ios2.ne.0) exit readlines 
             else 
               read(iunit,*,iostat=ios3) alam,anum,gf,excl 
               if(ios3.ne.0) exit readlines 
               gr=2.4e13/alam**2 
               gs=gsstd 
               gw=gwstd 
            end if 
          else 
            if(np.eq.9) then 
               read(iunit,*,iostat=ios4) alam,anum,gf,excl,gr,gh2,xnh2,   & 
     &                              ghe,xnhe 
               if(ios4.ne.0) exit readlines 
             else if(np.eq.7) then 
               READ(IUNIT,*,iostat=ios5) ALAM,ANUM,GF,EXCL,GR,GS,GW 
               if(ios5.ne.0) exit readlines 
             else 
               read(iunit,iostat=ios6) alam,anum,gf,excl 
               if(ios.ne.0) exit readlines 
               gr=2.4e13/alam**2 
               gs=gsstd 
               gw=gwstd 
            end if 
         end if 
! 
!     change wavelength to vacuum for lambda > 2000 
! 
         if(alam.gt.200..and.vaclim.gt.2000.) then 
            wl0=alam*10. 
            ALM=1.E8/(WL0*WL0) 
            XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM) 
            WL0=WL0*(XN1*1.e-6+UN) 
            alam=wl0*0.1 
         END IF 
! 
!        first selection : for a given interval 
! 
         IF(ALAM.GT.ALASTM(ILIST)+CUTOFF) EXIT READLINES 
! 
!        second selection : for line strengths 
! 
         FR0=CNM/ALAM 
         icod=int(anum+tenm4) 
         ion=int((anum-float(icod)+tenm4)*100.) 
 
         if(ion.eq.0) then 
            imol=molind(icod) 
          else 
            imol=ionind(icod) 
         end if 
 
         if(imol.le.0.or.imol.gt.nmolec) cycle readlines 
         EXCL=ABS(EXCL) 
         GFP=C1*GF-C2 
         EPP=C3*EXCL 
         gx=gfp-epp/tstd 
         ab0=0. 
! 
!     ************************************** 
!     rejecting weak lines 
!     ************************************* 
! 
         REJ: if(ndstep.eq.0.and.ifwin.eq.0) then 
! 
!           old procedure for line rejection 
! 
            if(gx.gt.-30)                                                 & 
     &      AB0=EXP(GFP-EPP/TSTD)*RRMOL(IMOL,IDSTD)/DOPSTD/AVAB 
            IF(AB0.LT.UN) CYCLE READLINES 
          else 
! 
!      new procedure for line rejection 
! 
            do ijcn=ijc,nfreqc 
               if(fr0.ge.freqc(ijcn)) exit 
            end do 
            ijc=ijcn 
            if(ijc.gt.nfreqc) ijc=nfreqc 
! 
            tkm=1.65e8/ammol(imol) 
            DP0=3.33564E-11*FR0 
            do id=1,nd,ndstep 
               td=temp(id) 
               gx=gfp-epp/td 
               ab0=0. 
               if(gx.gt.-30) then 
                  dops=dp0*sqrt(tkm*td+vturb(id)) 
                  AB0=EXP(gx)*RRMOL(IMOL,ID)/(DOPS*abstdw(ijc,id)*relop) 
               end if 
               if(ab0.ge.un) exit REJ 
            end do 
            cycle readlines 
         end if REJ 
! 
!        truncate line list if there are more lines than maximum allowable 
!        (given by MLIN0) 
! 
         IL=IL+1 
         IF(IL.GT.MLINM0) THEN 
            WRITE(6,                                                      & 
     &      "(' **** MORE LINES THAN MLINM0, LINE LIST TRUNCATED '/       & 
     &      '       AT LAMBDA',F15.4,'  NM'/)") ALAM 
            IL=MLINM0 
            ALASTM(ILIST)=CNM/FREQM(IL,ILIST)-CUTOFF 
            FRLASM(ILIST)=CNM/ALASTM(ILIST) 
            NXTSEM(ILIST)=1 
            EXIT READLINES 
         END IF 
! 
!        ============================================= 
!        line is selected, set up necessary parameters 
!        ============================================= 
! 
!        evaluation of EXTIN0 - the distance (in delta frequency) where 
!        the line is supposed to contribute to the total opacity 
! 
         EX0=AB0*ASTD*10. 
         EXT=EXT0 
         IF(EX0.GT.TEN) EXT=SQRT(EX0) 
         EXTIN0=EXT*DOPSTD 
! 
!        store parameters for selected lines 
! 
         FREQM(IL,ILIST)=FR0 
         EXCLM(IL,ILIST)=real(EPP) 
         GFM(IL,ILIST)=real(GFP) 
         EXTINM(IL,ILIST)=real(EXTIN0) 
         INDATM(IL,ILIST)=imol 
!         if(ilist.eq.3.or.ilist.eq.14) write(*,"(i5,2f12.4)")            & 
!    &       il,alam,cnm/freqm(il,ilist)
! 
!        ****** non-standard line broadening parameters ***** 
! 
         ibro=ibroli(ilist) 
         if(ibro.ne.0) then 
            GRM(IL,ibro)=real(GR*PI4) 
            GSM(IL,ibro)=real(GS*PI4*3.125e-5) 
            GWM(IL,ibro)=real(GW*PI4) 
         end if 
 
         ivdw=ivdwli(ilist) 
         if(ivdw.ne.0) then 
             gvdwh2(il,ivdw)=real(gh2) 
             gexph2(il,ivdw)=real(xnh2) 
             gvdwhe(il,ivdw)=real(ghe) 
             gexphe(il,ivdw)=real(xnhe) 
         end if 
! 
      END DO READLINES 
! 
!     ================================== 
!     end of reading 
!     ================================== 
 
      NLINM0(ILIST)=IL 
      nlinmt(ilist)=nlinmt(ilist)+nlinm0(ilist) 
      alend(ilist)=cnm/fr0 
      if(il.gt.0) alend(ilist)=cnm/freqm(il,ilist)
! 
      xln=float(il)/1000 
 
      if(il.lt.100000) then 
       write(6,"(' MOLECULAR LINES - FROM LIST:',i4,'   FILE:  ',a,':'/   & 
     & 36x,'WITH SELECTED ',i10,'   LINES')")                             & 
     & ilist,trim(amlist(ilist)),il 
      else if(il.lt.1000000) then 
       write(6,"(' MOLECULAR LINES - FROM LIST:',i4,'   FILE:  ',a,':'/   & 
     & 36x,'WITH SELECTED ',f10.3,' K  LINES')")                          & 
     & ilist,trim(amlist(ilist)),xln 
      else 
       write(6,"(' MOLECULAR LINES - FROM LIST:',i4,'   FILE:  ',a,':'/   & 
     & 36x,'WITH SELECTED ',f10.3,' M  LINES')")                          & 
     & ilist,trim(amlist(ilist)),xln/1.e3 
       end if 
!      write(6,*) 
       if(nlinm0(ilist).gt.0) then           
!         write(6,"('LINES FOR LIST',i4,'  BETWEEN ',2f12.3)")            &
!    &    ilist,cnm/freqm(1,ilist),alend(ilist)
          write(6,"(19x,'BETWEEN ',2f12.3)")                              &
     &    cnm/freqm(1,ilist),alend(ilist)
        else
          inactm(ilist)=1
          write(6,"(' NO SELECTED LINES - SET INACTIVE')")
       end if
!      write(*,*) '----------------------------------------------------' 
!      write(*,*) 
 
      RETURN 
      END SUBROUTINE INMOLI 
! 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE MOLSET(ILIST) 
!     ======================== 
! 
!     Selection of molecular lines that may contribute, 
!     set up auxiliary fields containing line parameters. 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use lindat 
      use molist 
 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      INTEGER, SAVE :: IMLAST 
! 
      REAL(DP), PARAMETER :: CNM=2.997925e17 
! 
!     write(*,*) 'molset-ilist,inactm,alend', 
!    *   ilist,inactm(ilist),alend(ilist) 
      if(inactm(ilist).ne.0) return 
      IL0=0 
      IPRSEM(ILIST)=0 
      NLINM=0 
      IREADM(ILIST)=1 
      IF(IBLANK.LE.1.OR.IMODE.EQ.1.OR.IMODE.EQ.-1) IREADM(ILIST)=0 
      IF(IBLANK.LE.1) APREV=0. 
      ALA0=CNM/FREQ(1) 
      ALA1=CNM/FREQ(2) 
! 
!     skip if current wavelength larger than the largest wavelngth in the 
!     line list 
! 
!!    write(*,*)
!!    write(*,"('MOLSET-BEG',5i6,2f10.3)")                                &
!!   & ilist,iblank,ireadm(ilist),nlinml(ilist),ilastm(ilist),ala0,ala1 
      if(ala0.gt.alend(ilist)) then 
         inactm(ilist)=1 
         write(6,"(/'*** LIST ',i4,' SET INACTIVE AT IBLANK ',i5)")       &
     &         ilist,iblank
         write(6,"(15x,'BECAUSE ALA0 ',f12.3,' > ALEND',f12.3)")          &
     &         ala0,alend(ilist)
         return 
      end if 
! 
      FRMINM=CNM/ALA0 
      FRM=FRMINM 
      SPACE=SPACE0 
      IF(ALAMC.GT.0.) SPACE=SPACE0*ALA0/ALAMC 
      IF(SPACE0.LT.0.) SPACE=-SPACE0 
 
      CUTOFF=CUTOF0*0.2 
      DOPSTD=1.E7/ALA0*DSTD 
      DISTAN=0.15*DOPSTD 
      SPAC=3.E16/ALA0/ALA0*SPACE 
      DISTA0=0.14*SPAC 
      IF(IBLANK.GE.2.AND.IMODE.EQ.-1) IL0=IMLAST 
      FRLI0=FRMINM 
      ASTD=1.0 
      AVAB=ABSTD(IDSTD)*RELOP 

!     if(ilist.eq.1) write(*,"(//'******* MgH',2i5)")                     &
!    & iblank,nlinml(ilist)
! 
!     set up indices of lines 
!     IL0 - is the current index of line in the numbering of all lines 
! 
      mm=3
      SELLINES: DO 
!!       if(ilist.eq.mm) write(*,"(/'sel1',6i10)")                        &
!!   &      ilist,ireadm(ilist),iprsem(ilist),il0,                        &
!!   &      inmlip(iprsem(ilist),ilist),nlinml(ilist)
!!       IF(ilist.eq.mm.and.IPRSEM(ILIST).GT.10) STOP
         IF(IREADM(ILIST).EQ.1) THEN 
            IPRSEM(ILIST)=IPRSEM(ILIST)+1
            if((iprsem(ilist).le.nlinml(ilist).and.nlinml(ilist).gt.0)    &
     &         .or.nlinml(ilist).eq.0) then 
               IL0=INMLIP(IPRSEM(ILIST),ILIST) 
!!             if(ilist.eq.mm) write(*,"('sel2',5i10,2f10.3)")            &
!!   &         ilist,ireadm(ilist),iprsem(ilist),il0
            end if
            if(nlinml(ilist).eq.0) il0=ilastm(ilist)
            IF(FREQM(IL0,ILIST).LT.FRMINM.                                &
     &      or.iprsem(ilist).gt.nlinml(ilist).and.nlinml(ilist).gt.0)     &
     &      THEN 
!!             if(ilist.eq.mm)                                            &
!!   &         write(*,"('sel4',3i8)") ilist,                             &
!!   &         INMLIP(IPRSEM(ILIST)-1,ILIST),il0
               IREADM(ILIST)=0 
!!             if(ilist.eq.mm) write(*,"('sel5',i8)")
!!   &         IPRSEM(ILIST)-1
               IL0=INMLIP(IPRSEM(ILIST)-1,ILIST)+1 
!!             if(ilist.eq.mm) write(*,"('sel6',2i8)") ilist,il0
            END IF 
          ELSE 
            IL0=IL0+1 
!!          if(ilist.eq.mm) write(*,"('sel7',2i8)") ilist,il0
         END IF 
!!       if(ilist.eq.mm) write(*,"('SELF',5i10)")                         &
!!   &    ilist,ireadm(ilist),iprsem(ilist),il0,nlinm0(ilist)
         IF(IL0.GT.NLINM0(ILIST)) EXIT SELLINES 
         FRLIM=FRLI0 
         FR0=FREQM(IL0,ILIST) 
         ALAM=CNM/FR0 
         if(il0.eq.nlinm0(ilist).and.alam.lt.ala0) exit sellines
         if(il0.gt.ilastm(ilist)) ilastm(ilist)=il0
!!       if(ilist.eq.mm)                                                  &
!!   &    write(*,"('==before sel',i5,f12.3)") il0,alam
!        if(ilist.eq.mm) write(*,*) 'SiH',il0
! 
         IF(ALAM.LT.ALA0-CUTOFF) CYCLE SELLINES 
         IF(ALAM.GT.ALA1+CUTOFF) EXIT SELLINES 
! 
!        SECOND SELECTION : FOR LINE STRENGHTS 
! 
         EXT=EXTINM(IL0,ILIST) 
         FRLI0=FR0-EXT-SPAC 
         IF(FRLI0.GT.FRLIM) FRLI0=FRLIM 
!!       if(ilist.eq.mm) write(*,"(i6,1p8e12.4)")                         &
!!   &   il0,fr0,ext,alam,ala0,fr0-frminm,ext+spac
         IF(ALAM.LT.ALA0.AND.FR0-FRMINM.GT.EXT+SPAC) CYCLE SELLINES 
         IF(FREQ(NFREQS)-FR0.GT.EXT+SPAC) CYCLE SELLINES 
! 
         NLINM=NLINM+1 
         if(nlinm.gt.mlinm) then 
            write(*,*) 'nlinm,mlinm',nlinm,mlinm 
            call quit('too many molecular lines in a set') 
         end if 
         INMLIN(NLINM,ILIST)=IL0 
!!       if(ilist.eq.mm) write(*,*) '** after sel.',il0,nlinm,alam 
      END DO SELLINES 
! 
!     frequency indices of the line centers 
! 
      XX=FREQ(2)-FREQ(1) 
      DFRCON=NFREQ-3 
      DFRCON=-DFRCON/XX 
      IFRCON=INT(DFRCON) 
      CENTR: DO IL=1,NLINM 
         fr0=freqm(inmlin(il,ilist),ILIST) 
         XJC=3.+DFRCON*(FREQ(1)-FR0) 
         IJC=INT(XJC) 
         IJCMTR(IL,ILIST)=IJC 
         if(ijc.le.3.or.ijc.ge.nfreq) cycle centr 
         if(fr0.lt.freq(ijc)) then 
            ijc0=ijc 
            dfr0=freq(ijc0)-fr0 
            finbe: do 
               ijc0=ijc0+1 
               dfr=abs(freq(ijc0)-fr0) 
               if(dfr.lt.dfr0) then 
                  ijc=ijc0 
                  ijc0=ijc0+1 
                  dfr0=dfr 
                  cycle finbe 
                else 
                  exit finbe 
               end if 
            end do finbe 
          else if(fr0.gt.freq(ijc)) then 
            ijc0=ijc 
            dfr0=fr0-freq(ijc0) 
            finen: do 
               ijc0=ijc0-1 
               dfr=abs(freq(ijc0)-fr0) 
               if(dfr.lt.dfr0) then 
                  ijc=ijc0 
                  ijc0=ijc0-1 
                  dfr0=dfr 
                  cycle finen 
                else 
                  exit finen 
               end if 
            end do finen 
         end if 
         IJCMTR(IL,ILIST)=IJC 
      END DO CENTR 
! 
      DO  IL=1,NLINM 
         INMLIP(IL,ILIST)=INMLIN(IL,ILIST) 
      END DO 
      NLINML(ILIST)=NLINM 
      IMLAST=INMLIN(NLINML(ILIST),ILIST) 
!!     write(6,"(70x,'ILIST,NLINM,IMLAST',3i7)") ilist,nlinm,imlast

      if(nlinm.eq.0.and.alend(ilist).lt.ala0) then
         inactm(ilist)=1
         write(6,"(/'*** LIST ',i3,' SET INACTIVE at IBLANK ',i5)")       &  
     &         ILIST,IBLANK
         write(6,"(15x,'BECAUSE ALAST',f12.3,' >  AL0 ',f12.3)")          &
     &         alend(ilist),ala0
      end if
! 
      CALL INIBLM 
!      write(*,*) 'after iniblm'
! 
      RETURN 
      END SUBROUTINE MOLSET 
! 
! 
! ******************************************************************** 
! 
! 
 
      SUBROUTINE INIBLM 
!     ================= 
! 
!     driving procedure for treating a partial molecular line list for the 
!     current wavelength region 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use lindat 
      use molist 
 
      implicit real(dp) (a-h,o-z),logical (l) 
! 
      REAL(DP), PARAMETER :: DP0=3.33564E-11, DP1=1.651E8, UN=1. 
! 
      XX=FREQ(1) 
      IF(NFREQ.GE.2) XX=0.5*(FREQ(1)+FREQ(2)) 
      BNU=BN*(XX*1.E-15)**3 
      HKF=HK*XX 
      DO ID=1,ND 
         T=TEMP(ID) 
         EXH=EXP(HKF/T) 
         EXHK(ID)=UN/EXH 
         PLAN(ID)=BNU/(EXH-UN) 
         STIM(ID)=UN-EXHK(ID) 
         DO IMOL=1,NMOLEC 
            IF(AMMOL(IMOL).GT.0.)                                         & 
     &      DOPMOL(IMOL,ID)=UN/(XX*DP0*SQRT(DP1*T/AMMOL(IMOL)+            & 
     &                      VTURB(ID))) 
         END DO 
      END DO 
      RETURN 
      END SUBROUTINE INIBLM 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE IDMTAB 
!     ================= 
! 
!     output of selected molecular line parameters (identification table) 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use lindat 
      use molist 
 
      implicit real(dp) (a-h,o-z),logical (l) 
 
 
      CHARACTER(LEN=4) :: APB,AP0,AP1,AP2,AP3,AP4,APR 
! 
 
      REAL(DP), PARAMETER :: C1=2.3025851, C2=4.2014672, C3=1.4387886,    & 
     &                       PI4=7.95774715E-2 
 
      DATA APB,AP0,AP1,AP2,AP3,AP4 /'    ','   .','   *','  **',' ***',   & 
     &                              '****'/ 
! 
      IF(IMODE.LE.-3.AND.TEMP(1).GT.TMOLIM) RETURN 
      ALM0=2.997925e18/FREQ(1) 
      ALM1=2.997925e18/FREQ(2) 
      if(ifwin.gt.0) ALM1=2.997925e18/FREQ(NFREQ) 
      IF(IPRIN.LE.-2) RETURN 
      if(iprin.ge.3) then 
         IF(IMODE.GE.0) WRITE(6,"(/' ',I4,'. SET (MOLECULAR LINES):',     & 
     & ' INTERVAL  ',F9.3,' -',F9.3,' ANGSTROMS'/                         & 
     &        ' ------------')") IBLANK,ALM0,ALM1 
         IF(IMODE.GE.0.OR.(IMODE.EQ.-1.AND.IBLANK.EQ.1))                  & 
     &   WRITE(6,"(/14X,                                                  & 
     & 'LAMBDA  MOLECULE  LOG GF       ELO    LINE/CONT',2X,              & 
     & 'EQ.WIDTH',8x,'AGAM'/)") 
      end if 
! 
      ID=IDSTD 
      LISTLOOP: DO ILIST=1,NMLIST 
         IF(NLINML(ILIST).EQ.0) CYCLE LISTLOOP 
         ivdw=ivdwli(ilist) 
         ibro=ibroli(ilist) 
         lbstd=ivdw.eq.0.and.ibro.eq.0 
         if(lbstd) then 
            gs=gsstd*3.125e-5 
            gw=gwstd 
         end if 
 
         LINELOOP: DO IL0=1,NLINML(ILIST) 
            IL=INMLIN(IL0,ILIST) 
            ALAM=2.997925D18/FREQM(IL,ILIST) 
            IJCN=IJCMTR(IL0,ILIST) 
            IMOL=INDATM(IL,ILIST) 
            DOP1=DOPMOL(IMOL,ID) 
            ANE=ELEC(ID) 
            IF(LBSTD) THEN 
               AGAM=(2.4E15/ALAM**2+GS*ANE+GW*VDWC(ID))*PI4*DOP1 
             ELSE IF(IBRO.GT.0) THEN 
               AGAM=(GRM(IL,IBRO)+GSM(IL,IBRO)*ANE+GWM(IL,IBRO))*DOP1 
             ELSE IF(IVDW.GT.0) THEN 
               AGAM=GVDW(IL,IVDW,ID)*DOP1 
            END IF 
            ABCNT=EXP(GFM(IL,ILIST)-EXCLM(IL,ILIST)/TEMP(ID))*            & 
     &         RRMOL(IMOL,ID)*DOP1*STIM(ID) 
            absta=min(ch(1,id),ch(2,id)) 
            str0=abcnt/absta 
            if(ifwin.gt.0) STR0=ABCNT/ABSTDW(IJCONT(IL),ID) 
            GF=(GFM(IL,ILIST)+C2)/C1 
            EXCL=EXCLM(IL,ILIST)/C3 
            IF(STR0.LE.1.2) THEN 
               WW1=0.886*STR0*(1.-STR0*(0.707-STR0*0.577)) 
             ELSE 
               WW1=SQRT(LOG(STR0)) 
            END IF 
            IF(STR0.GT.55.) THEN 
               WW2=0.5*SQRT(3.14*AGAM*STR0) 
               IF(WW2.GT.WW1) WW1=WW2 
            END IF 
            EQW=ALAM/FREQM(IL,ILIST)*1.E3/DOP1*WW1 
            STR=EQW*10. 
            APR=APB 
            IF(STR.GE.1.E0.AND.STR.LT.1.E1) APR=AP0 
            IF(STR.GE.1.E1.AND.STR.LT.1.E2) APR=AP1 
            IF(STR.GE.1.E2.AND.STR.LT.1.E3) APR=AP2 
            IF(STR.GE.1.E3.AND.STR.LT.1.E4) APR=AP3 
            IF(STR.GE.1.E4) APR=AP4 
            if(alam.ge.alm0.and.alam.lt.alm1) then 
               WRITE(15,"(F11.3,2X,A4,4X,F7.2,F12.3,1PE11.2,0PF8.1,       & 
     &         1X,A4,i4,1PE10.2)") ALAM,CMOL(IMOL),GF,EXCL,               & 
     &         STR0,EQW,APR,id,AGAM 
!              WRITE(15,"(F11.3,2X,A4,4X,F7.2,F12.3,1PE11.2,0PF8.1, 
!    *         1X,A4,1P4E10.2)") ALAM,CMOL(IMOL),GF,EXCL, 
!    *         STR0,EQW,APR,abcnt,rrmol(imol,id),dop1,stim(id) 
            end if 
         END DO LINELOOP 
      END DO LISTLOOP 
 
      RETURN 
      END SUBROUTINE IDMTAB 
! 
! 
! ******************************************************************** 
! 
! 
 
      SUBROUTINE MOLOP(ID,ABLIN,EMLIN,AVAB,ILIST) 
!     =========================================== 
! 
!     Total molecular line opacity (ABLIN) and emissivity (EMLIN) 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use lindat 
      use molist 
 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: UN     = 1.,                                 & 
     &                       EXT0   = 3.17,                               & 
     &                       TEN    = 10.,                                & 
     &                       PI4    = 7.95774715E-2 
      REAL(DP)            :: ABLIN(MFREQ),EMLIN(MFREQ) 
! 
      DO IJ=1,NFREQ 
         ABLIN(IJ)=0. 
         EMLIN(IJ)=0. 
      END DO 
! 
!     write(*,*) 'molop-ilist,nlinml,inactm',ilist,nlinml(ilist), 
!    *            inactm(ilist) 
      if(temp(id).gt.tmolim) return 
      IF(NLINML(ILIST).EQ.0) RETURN 
      if(inactm(ilist).ne.0) return 
 
      ivdw=ivdwli(ilist) 
      ibro=ibroli(ilist) 
      lbstd=ivdw.eq.0.and.ibro.eq.0 
      if(lbstd) then 
         gs=gsstd*3.125e-5 
         gw=gwstd 
      end if 
! 
!     overall loop over contributing lines 
! 
      TEM1=UN/TEMP(ID) 
      ANE=ELEC(ID) 
      DO I=1,NLINML(ILIST) 
         IL=INMLIN(I,ILIST) 
         IMOL=INDATM(IL,ILIST) 
         FR0=FREQM(IL,ILIST) 
         DOP1=DOPMOL(IMOL,ID) 
         IF(LBSTD) THEN 
            ALAM=2.997925E17/FR0 
            AGAM=(2.4E13/ALAM**2+GS*ANE+GWSTD*VDWC(ID))*PI4*DOP1 
          ELSE IF(IBRO.GT.0) THEN 
            AGAM=(GRM(IL,IBRO)+GSM(IL,IBRO)*ANE+GWM(IL,IBRO))*DOP1 
          ELSE IF(IVDW.GT.0) THEN 
            AGAM=GVDW(IL,IVDW,ID)*DOP1 
         END IF 
!        AGAM=(GRM(IL,ILIST)+GSM(IL,ILIST)*ANE+ 
!    *         GVDW(IL,ILIST,ID))*DOP1 
         AB0=EXP(GFM(IL,ILIST)-EXCLM(IL,ILIST)*TEM1)*RRMOL(IMOL,ID)*      & 
     &           DOP1*STIM(ID) 
! 
!        set up limiting frequencies where the line I is supposed to 
!        contribute to the opacity 
! 
         EX0=AB0/AVAB*AGAM 
         EXT=EXT0 
         IF(EX0.GT.TEN) EXT=SQRT(EX0) 
         EXT=EXT/DOP1 
         XIJEXT=DFRCON*EXT+1.5 
         IJ1=int(MAX(float(IJCMTR(I,ILIST))-XIJEXT,3.)) 
         IJ2=int(MIN(float(IJCMTR(I,ILIST))+XIJEXT,float(NFREQS))) 
         IF(IJ1.LT.NFREQ.AND.IJ2.GT.2) THEN 
            DO IJ=IJ1,IJ2 
               XF=ABS(FREQ(IJ)-FR0)*DOP1 
               ABLIN(IJ)=ABLIN(IJ)+AB0*VOIGTK(AGAM,XF) 
            END DO 
         END IF 
      END DO 
! 
      DO IJ=3,NFREQ 
         EMLIN(IJ)=EMLIN(IJ)+ABLIN(IJ)*PLAN(ID) 
      END DO 
! 
      RETURN 
      END SUBROUTINE MOLOP 
! 
! 
! ******************************************************************** 
! 
! 
      FUNCTION SBFHMI(FR) 
!     =================== 
! 
!     Bound-free cross-section for H- (negative hydrogen ion) 
!     Taken from Kurucz ATLAS9 
! 
!     FROM MATHISEN (1984), AFTER WISHART(1979) AND BROAD AND REINHARDT (1976) 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: WBF(85),BF(85) 
      DATA WBF/  18.00,  19.60,  21.40,  23.60,  26.40,  29.80,  34.30,   & 
     &   40.40,  49.10,  62.60, 111.30, 112.10, 112.67, 112.95, 113.05,   & 
     &  113.10, 113.20, 113.23, 113.50, 114.40, 121.00, 139.00, 164.00,   & 
     &  175.00, 200.00, 225.00, 250.00, 275.00, 300.00, 325.00, 350.00,   & 
     &  375.00, 400.00, 425.00, 450.00, 475.00, 500.00, 525.00, 550.00,   & 
     &  575.00, 600.00, 625.00, 650.00, 675.00, 700.00, 725.00, 750.00,   & 
     &  775.00, 800.00, 825.00, 850.00, 875.00, 900.00, 925.00, 950.00,   & 
     &  975.00,1000.00,1025.00,1050.00,1075.00,1100.00,1125.00,1150.00,   & 
     & 1175.00,1200.00,1225.00,1250.00,1275.00,1300.00,1325.00,1350.00,   & 
     & 1375.00,1400.00,1425.00,1450.00,1475.00,1500.00,1525.00,1550.00,   & 
     & 1575.00,1600.00,1610.00,1620.00,1630.00,1643.91/ 
      DATA BF/   0.067,  0.088,  0.117,  0.155,  0.206,  0.283,  0.414,   & 
     &   0.703,   1.24,   2.33,  11.60,  13.90,  24.30,  66.70,  95.00,   & 
     &   56.60,  20.00,  14.60,   8.50,   7.10,   5.43,   5.91,   7.29,   & 
     &   7.918,  9.453,  11.08,  12.75,  14.46,  16.19,  17.92,  19.65,   & 
     &   21.35,  23.02,  24.65,  26.24,  27.77,  29.23,  30.62,  31.94,   & 
     &   33.17,  34.32,  35.37,  36.32,  37.17,  37.91,  38.54,  39.07,   & 
     &   39.48,  39.77,  39.95,  40.01,  39.95,  39.77,  39.48,  39.06,   & 
     &   38.53,  37.89,  37.13,  36.25,  35.28,  34.19,  33.01,  31.72,   & 
     &   30.34,  28.87,  27.33,  25.71,  24.02,  22.26,  20.46,  18.62,   & 
     &   16.74,  14.85,  12.95,  11.07,  9.211,  7.407,  5.677,  4.052,   & 
     &   2.575,  1.302, 0.8697, 0.4974, 0.1989,    0. / 
!    Bell and Berrington J.Phys.B,vol. 20, 801-806,1987. 
! 
      HMINBF=0. 
      IF(FR.GT.1.82365E14) THEN 
         WAVE=2.99792458E17/FR 
         HMINBF=YLINTP(WAVE,WBF,BF,85,85)*1.E-18 
      END IF 
      SBFHMI=HMINBF 
      RETURN 
      END FUNCTION SBFHMI 
! 
! 
! ******************************************************************** 
! 
! 
 
      FUNCTION SFFHMI(POPI,FR,T) 
!     ========================== 
! 
!     Free-free cross-section for H- (negative hydrogen ion) 
!     Taken from Kurucz ATLAS9 
! 
!     From Bell and Berrington J.Phys.B,vol. 20, 801-806,1987. 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: CONFF=5040.*1.380658E-16, CONTH=5040. 
      REAL(DP), SAVE :: WFFLOG(22),FFLOG(22,11) 
      REAL(DP) :: FFCS(11,22),FFLOG2(22) 
      REAL(DP) :: FFBEG(11,11),FFEND(11,11),FFTT(11) 
      REAL(DP) :: THETAFF(11),WAVEK(22) 
      EQUIVALENCE (FFCS(1,1),FFBEG(1,1)),(FFCS(1,12),FFEND(1,1)) 
! 
      DATA WAVEK/.50,.40,.35,.30,.25,.20,.18,.16,.14,.12,.10,.09,.08,     & 
     & .07,.06,.05,.04,.03,.02,.01,.008,.006/ 
      DATA THETAFF/                                                       & 
     &  0.5,  0.6, 0.8,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.8,  3.6/ 
      DATA FFBEG/                                                         & 
     &.0178,.0222,.0308,.0402,.0498,.0596,.0695,.0795,.0896, .131, .172,  & 
     &.0228,.0280,.0388,.0499,.0614,.0732,.0851,.0972, .110, .160, .211,  & 
     &.0277,.0342,.0476,.0615,.0760,.0908, .105, .121, .136, .199, .262,  & 
     &.0364,.0447,.0616,.0789,.0966, .114, .132, .150, .169, .243, .318,  & 
     &.0520,.0633,.0859, .108, .131, .154, .178, .201, .225, .321, .418,  & 
     &.0791,.0959, .129, .161, .194, .227, .260, .293, .327, .463, .602,  & 
     &.0965, .117, .157, .195, .234, .272, .311, .351, .390, .549, .711,  & 
     & .121, .146, .195, .241, .288, .334, .381, .428, .475, .667, .861,  & 
     & .154, .188, .249, .309, .367, .424, .482, .539, .597, .830, 1.07,  & 
     & .208, .250, .332, .409, .484, .557, .630, .702, .774, 1.06, 1.36,  & 
     & .293, .354, .468, .576, .677, .777, .874, .969, 1.06, 1.45, 1.83/ 
      DATA FFEND/                                                         & 
     & .358, .432, .572, .702, .825, .943, 1.06, 1.17, 1.28, 1.73, 2.17,  & 
     & .448, .539, .711, .871, 1.02, 1.16, 1.29, 1.43, 1.57, 2.09, 2.60,  & 
     & .579, .699, .924, 1.13, 1.33, 1.51, 1.69, 1.86, 2.02, 2.67, 3.31,  & 
     & .781, .940, 1.24, 1.52, 1.78, 2.02, 2.26, 2.48, 2.69, 3.52, 4.31,  & 
     & 1.11, 1.34, 1.77, 2.17, 2.53, 2.87, 3.20, 3.51, 3.80, 4.92, 5.97,  & 
     & 1.73, 2.08, 2.74, 3.37, 3.90, 4.50, 5.01, 5.50, 5.95, 7.59, 9.06,  & 
     & 3.04, 3.65, 4.80, 5.86, 6.86, 7.79, 8.67, 9.50, 10.3, 13.2, 15.6,  & 
     & 6.79, 8.16, 10.7, 13.1, 15.3, 17.4, 19.4, 21.2, 23.0, 29.5, 35.0,  & 
     & 27.0, 32.4, 42.6, 51.9, 60.7, 68.9, 76.8, 84.2, 91.4, 117., 140.,  & 
     & 42.3, 50.6, 66.4, 80.8, 94.5, 107., 120., 131., 142., 183., 219.,  & 
     & 75.1, 90.0, 118., 144., 168., 191., 212., 234., 253., 325., 388./ 
      DATA ISTART/0/ 
! 
      IF(ISTART.EQ.0) THEN 
      ISTART=1 
      DO IWAVE=1,22 
         WFFLOG(IWAVE)=LOG(91.134/WAVEK(IWAVE)) 
         DO ITHETA=1,11 
            FFLOG(IWAVE,ITHETA)=LOG(FFCS(ITHETA,IWAVE)*1.E-26) 
         END DO 
      END DO 
      ENDIF 
! 
      WAVE=2.99792458E17/FR 
      WAVELOG=LOG(WAVE) 
! 
      DO ITHETA=1,11 
         DO IWAVE=1,22 
            FFLOG2(IWAVE)=FFLOG(IWAVE,ITHETA) 
         END DO 
         FFTLOG=YLINTP(WAVELOG,WFFLOG,FFLOG2,22,22) 
         FFTT(ITHETA)=EXP(FFTLOG)/THETAFF(ITHETA)*CONFF 
      END DO 
! 
      THETA=CONTH/T 
      FFTH=YLINTP(THETA,THETAFF,FFTT,11,11) 
      SFFHMI=FFTH*POPI/(1.-exp(-hk*fr/t)) 
      RETURN 
      END FUNCTION SFFHMI 
! 
! 
! 
!     ****************************************************************** 
! 
! 
! ========================================================================= 
! ************************************************************************* 
! ************************************************************************* 
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
 
!     data indtsu / 2,  5, 12, 4, 8, 7, 6, 
!    *              9, 11, 10, 29, 50, 59, 46, 132, 52, 19, 
!    *             13, 42, 38, 39, 37, 44, 36, 14, 118, 33, 
!    *              3, 16, 57, 32, 49, 60, 54, 41, 107,  0, 
!    *            148, 152, 153, 155, 0, 17, 24, 25, 28, 51, 
!    *            112, 119,   0,   0,21, 15, 43, 56,  0, 64, 
!    *             47,  65,   0,  61, 0, 62,118, 40, 66/ 
!     data indtsu / 2,  5, 12, 4, 8, 7, 6, 
!    *              9, 11, 10, 29, 50, 59, 46, 132,  52, 19, 
!    *             13, 42, 38, 39, 37, 44, 36,  14, 117, 33, 
!    *              3, 16, 57, 32, 49, 60, 54,  41, 106,303, 
!    *            147, 151, 152, 154, 302, 17,  24,  25, 28, 51, 
!    *            111, 118, 102,   0,  21, 15,  43,  56,478, 64, 
!    *             47,  65, 413,  61, 190, 62 ,108,  40, 66,214, 
!    *             257*0./ 
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
! ========================================================================= 
! ************************************************************************* 
! ************************************************************************* 
! 
! 
      subroutine moleq(id,tt,an,aein,ane,ipri) 
!     ======================================== 
! 
!     calculation of the equilibrium state of atoms and molecules 
! 
!     Input:  id    - depth point 
!             tt    - temperature [K] 
!             an    - number density 
!             aein  - initial estimate of the electron density 
! 
!     Output: ane    - electron density 
! 
!             rrr(id,j,i) - N/U for the atom with atomic number i and 
!                     ion j (j=1 for neutral, and j=2 for 1st ions) 
!             rrmol(imol,id) - N/U for the molecule with index imol 
!                     (the index is given by the ordering of 
!                     in the input file  tsuji.molec 
! 
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
!          if(id.eq.1) write(*,"('moleq',i5,a4,5f10.3)") 
!    &      i,typat(i),log10(anion(i,id)), 
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
! 
!         if(id.eq.1)                                                     & 
!    &    write(6,"(i5,a8,2f10.3,1p2e12.4,0pf12.3,1p2e12.4)")             & 
!    &          j,cmol(j),apmlog(j),pmoll,anden(jm),umoll,                & 
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
! 
      RETURN 
      END SUBROUTINE MOLEQ 
 
 
! 
! ========================================================================= 
! ************************************************************************* 
! ************************************************************************* 
! 
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
!      ELSE 
!        IF(PG/DHH.LE.0.1) THEN 
!           PH=PG/(1.0+HEH) 
!         ELSE 
!           PH=0.5 * (SQRT(DHH*(DHH+4.0 *PG/(1.0+HEH)))-DHH) 
!        END IF 
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
            IF(ITERAT.GT.150) THEN 
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
!    *        write(6,"('ion',3i5,a8,1p4e15.5)") 
!    *        niterr,j,m,cmol(j),p(99),pmolj,spnplu,spnion 
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
!    *     write(66,"(2i4,a8,i3,2f10.3,2x,5(2i3,f8.2,2x))") 
!    *     niterr,j,cmol(j),mmaxj,apmlog(j),pmoljl, 
!    *     (nelem(m,j),nato(m,j),log10(p(nelem(m,j))),m=1,mmaxj) 
 
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
!       write(6,601) niterr,perev0,spnion,pe,spnion/pe 
        PEREV=SQRT(PEREV/(1.0+SPNION/PE)) 
        DELTRS=DELTRS+ABS((PE-PEREV)/PE) 
!       if(iprin.gt.4) 
!    *     write(6,601) niterr,tem,pg*tk,fph*tk,pe*tk,perev*tk, 
!    *    (perev+pe)*0.5*tk,deltrs 
!       PE=(PEREV+PE)*0.5 
        pe=perev 
        P(99)=PE 
        IF(DELTRS.LE.EPS) THEN 
           EXIT RUSS 
         ELSE 
           NITERR=NITERR+1 
           IF(NITERR.LE.NIMAX) THEN 
              CYCLE RUSS 
            ELSE 
              WRITE(6,"('*DOES NOT CONVERGE AFTER ',I4,' ITERATIONS')")   & 
     &        NIMAX 
            EXIT RUSS 
           END IF 
        END IF 
      END DO RUSS 
! 
      if(iprin.ge.4) then 
         write(6,"('russel iterations ',i4,1p7e13.4)")                    & 
     &      niterr,tem,pg*tk,fph*tk,pe*tk,perev*tk,                       & 
     &      (perev+pe)*0.5*tk,deltrs 
         write(*,*) ' ' 
      end if 
! 
      RETURN 
      END SUBROUTINE RUSSEL 
 
! 
! 
!     *********************************************************************** 
! 
 
      subroutine equcon(t) 
!     ==================== 
! 
!     correction for equilibrium constant for molecular ions to express 
!     Kp from the B&C tables in a Tsuji-like form 
 
      use accura 
      use params, only: nmolec,ieqbc,amas 
      use eospar 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      real(dp), intent(in) :: t 
      real(dp) :: eaf(99) 
      integer  :: ineg(17) 
 
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
! 
!     tk=t*1.38054e-16         ! = k T 
!     betae2=4.8298e15*t**1.5  ! = 2 * (2 pi m_e k/h**2)**1.5 * T**1.5 
      conl=0.17602759          ! log10(2 * (2 pi m_e k/h**2)**1.5 * k) 
      TH=5040./T 
 
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
!          write(*,"(i4,3f11.3)") j,aplogj,corr,apmlog(j) 
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
!          write(*,"(i4,3f11.3)") j,aplogj,corr,apmlog(j) 
         end if 
      end do molloop 
 
      return 
      end subroutine equcon 
 
! 
!     *********************************************************************** 
! 
 
      subroutine bceqco(in,t,bcel) 
!     ============================ 
! 
!     new valuation of BC equilibrium constant from their Table 7 
! 
!     input:  in   - BC index 
!             t    - temperature 
!     output: bcel -log(Kp) 
 
      use accura 
      implicit none 
!     implicit real(dp) (a-h,o-z) 
 
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
!     implicit real(dp) (a-h,o-z) 
 
      integer, intent(in)   :: in 
      real(dp), intent(in)  :: t 
      real(dp), intent(out) :: bcnl 
      real(dp)              :: tt,a1 
 
      real(dp),save    :: tl(42),bcnel(7,42) 
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
! ******************************************************************** 
! 
! 
      SUBROUTINE SETWIN 
!     ================= 
! 
!     Initialisation of an extended radial structure 
!      (spherical symmetry is assumed) 
!     with a continuous connection between the lower quasi-hydrostatic 
!     layers and the upper, supersonic layers. The velocity structure 
!     in the upper layers is a beta-type law (v=vinf*(1-r0/r)^beta). 
! 
!     Additional input are read at the end of Unit 8: 
!      RCORE : Core radius (deepest layer, in solar radii or in cm) 
!      NDRAD : Number of layers 
!      NRCORE: Number of core rays 
!      INRV  : Switch indicating the data to be read: 
!           = 0 : Read an hydrostatic, plane-parallel model only; the 
!                   routine builds the radial points, density and 
!                   velocity structure; 
!           < 0 : Read also an hydrostatic, plane-parallel model, but 
!                   an empirical velocity law V(r) is read at each 
!                   radial point (r(id) is read); 
!           > 0 : Input from an extended model atmosphere; the velocity 
!                   law is read; the density structure is recomputed for 
!                   a possibly different mass-loss rate. 
!      XMDOT  : Mass loss rate (in solar mass/yr) 
!      BETAV, VINF : Parameters of the velocity law (VINF in km/s) 
!      RD, VEL: Radial points, expansion velocity 
! 
!     Synspec version 
! 
      use accura 
      use params 
      use modelp 
      use wincom 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: RSUN=6.96e10 
 
      call alloc_wincom 
! 
!     Read data for spherical atmosphere and velocity law 
! 
      READ(8,*,IOSTAT=IOS) RCORE,NDRAD,NRCORE,INRV,NFIRY,NDF 
      If(IOS.NE.0) RETURN 
      IF(RCORE.LT.1.E5) RCORE=RCORE*RSUN 
      IF(NDRAD.GT.MDEPTH) CALL quit('NDRAD too large') 
      READ(8,*) XMDOT,BETAV,VINF 
      XMDOT=6.30289e25*XMDOT 
      VINF=1.D5*VINF 
      ND=NDRAD 
      DO ID=1,ND 
        READ(8,*) RD(ID),VEL(ID),VTURB(ID),DENSCON(ID) 
        if(denscon(id).eq.0.) denscon(id)=1. 
        vturb(id)=vturb(id)*vturb(id) 
      END DO 
! 
!   Apply density contrast for clumping 
! 
      DO ID=1,ND 
        ELEC(ID) = ELEC(ID) * DENSCON(ID) 
        DENS(ID) = DENS(ID) * DENSCON(ID) 
        DO I=1,NLEVEL 
           POPUL(I,ID) = POPUL(I,ID) * DENSCON(ID) 
        END DO 
      END DO 
! 
!  Set up rays and weights 
! 
      itrad=1 
      call radtem 
      CALL SETRAY 
      CALL WGTJH1 
! 
      RETURN 
      END SUBROUTINE SETWIN 
! 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE SETRAY 
!     ================= 
! 
!     Setup impact rays and angles 
!      (assumes one impact ray tangent to every depth layer) 
! 
      use accura 
      use params 
      use modelp 
      use wincom 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: PI4=4.*3.141592654 
      REAL(DP), PARAMETER :: UN=1., TWO=2., HALF=0.5 
      REAL(DP) :: RS(MDEPF),RDX(MDEPF) 
      REAL(DP) ::  ZIU(MDEPTH),VIU(MDEPTH),ZIUF(MDEPF),VIUF(MDEPF) 
! 
!     Fine radial grid 
! 
      if(ndf.eq.0.or.ndf.eq.nd) then 
         ndf=nd 
         DO ID=1,NDF 
            DENSF(ID)=DENS(ID) 
         END DO 
      else 
         XR1=LOG(DENS(1)) 
         XR2=LOG(DENS(ND)) 
         DXR=(XR2-XR1)/FLOAT(NDF-1) 
         DO ID=1,NDF 
            DENSF(ID)=EXP(XR1+FLOAT(ID-1)*DXR) 
         END DO 
      end if 
! 
! 
!     Impact rays 
! 
      NREXT=ND 
      DO ID=1,NREXT 
        PIM(ID)=RD(ID) 
        NUD(ID)=ID 
      END DO 
      DO IU=1,NRCORE 
        PIM(NREXT+IU)=FLOAT(NRCORE-IU)/FLOAT(NRCORE)*RCORE 
        NUD(NREXT+IU)=ND 
      END DO 
      KMU=NREXT+NRCORE 
! 
!     Angles 
! 
      DO ID=1,ND 
        RD1=UN/RD(ID) 
        DO IU=ID,KMU 
          PRR=PIM(IU)*RD1 
          BMU(IU,ID)=SQRT(UN-PRR*PRR) 
        END DO 
      END DO 
! 
!     Depth increments along each ray 
! 
      DELZ(1,1)=0. 
      DFRQ(1,1)=0. 
      DO IU=2,KMU 
        NUDF(IU)=NUD(IU) 
        IU1=IU 
        IF(IU.GT.ND) IU1=ND 
        DO ID=1,IU1-1 
          DELZ(IU,ID)=BMU(IU,ID)*RD(ID)-BMU(IU,ID+1)*RD(ID+1) 
          DFRQ(IU,ID)=BMU(IU,ID)*VEL(ID)/CL 
          JD=2*NUD(IU)-ID 
          DFRQ(IU,JD)=-DFRQ(IU,ID) 
        END DO 
        DELZ(IU,IU1)=DELZ(IU,IU1-1) 
        DFRQ(IU,IU1)=0. 
        IF(IU.GT.NREXT) DFRQ(IU,ND)=BMU(IU,ND)*VEL(ND)/CL 
      END DO 
! 
! Finer grid along the NFIRY most external rays 
!   velocity steps DVD(ID) 
! 
      XMD4=XMDOT/PI4 
      CLV=UN/CL 
      DO ID=1,ND 
        DVD(ID)=SQRT(1.6e7*TEMP(ID)+VTURB(ID)) * 0.3 
!        DVD(ID)=SQRT(1.6e7*TEMP(ID)) 
      END DO 
      NUDX=ND 
      DO IU=2,NFIRY 
        IF(PIM(IU).GT.0.) THEN 
          DO ID=1,NUD(IU) 
            IID=NUD(IU)-ID+1 
            ZIU(ID)=VEL(IID) 
            VIU(ID)=DFRQ(IU,IID)*CL 
          ENDDO 
        ELSE 
          DO ID=1,NUD(IU) 
            IID=NUD(IU)-ID+1 
            ZIU(ID)=RD(IID) 
            VIU(ID)=DFRQ(IU,IID)*CL 
          ENDDO 
        ENDIF 
        NUDF(IU)=1 
        VIUF(1)=DFRQ(IU,1)*CL 
        DO ID=1,NUD(IU)-1 
          VZ1=DFRQ(IU,ID)*CL 
          VZ2=DFRQ(IU,ID+1)*CL 
          NFG=int((VZ1-VZ2)/DVD(ID))+1 
          XFG=(VZ1-VZ2)/DFLOAT(NFG) 
          IV0=NUDF(IU) 
          DO IV=1,NFG 
            VIUF(IV0+IV)=VZ1-DFLOAT(IV)*XFG 
          ENDDO 
          NUDF(IU)=NUDF(IU)+NFG 
          IF(NUDF(IU).GT.MDEPF )                                          & 
     &      CALL quit('Too many points in fine grid - SETRAY') 
        END DO 
        IF(NUDF(IU).GT.NUDX) NUDX=NUDF(IU) 
        INRP=2 
        IF(IU.GT.8) INRP=4 
        CALL INTERP(VIU,ZIU,VIUF,ZIUF,NUD(IU),NUDF(IU),INRP,0,0) 
        IF(PIM(IU).GT.0.) THEN 
          DO ID=1,NUDF(IU) 
            DMU=VIUF(ID)/ZIUF(ID) 
            RS(ID)=PIM(IU)/SQRT(UN-DMU*DMU) 
            DFRQF(IU,ID)=VIUF(ID)*CLV 
            VELF(IU,ID)=ZIUF(ID) 
            RDX(ID)=XMD4/(RS(ID)*RS(ID)*VELF(IU,ID)) 
            ZIUF(ID)=DMU*RS(ID) 
          END DO 
        ELSE 
          DO ID=1,NUDF(IU) 
            RS(ID)=ZIUF(ID) 
            DFRQF(IU,ID)=VIUF(ID)*CLV 
            VELF(IU,ID)=VIUF(ID) 
            RDX(ID)=XMD4/(RS(ID)*RS(ID)*VELF(IU,ID)) 
          END DO 
        END IF 
        IF(IU.LE.NREXT) THEN 
          DO ID=1,NUDF(IU) 
            JD=2*NUDF(IU)-ID 
            DFRQF(IU,JD)=-DFRQF(IU,ID) 
          END DO 
        END IF 
        DO ID=1,NUDF(IU)-1 
          DELZF(IU,ID)=ZIUF(ID)-ZIUF(ID+1) 
        END DO 
        DELZF(IU,NUDF(IU))=DELZF(IU,NUDF(IU)-1) 
! 
!   Assign depth index 
! 
        KRAY(IU,1)=2 
        DRAY(IU,1)=0. 
        IDK=1 
        DO ID=2,NUDF(IU) 
          DO WHILE (RDX(ID).GE.DENSF(IDK).and.idk.le.ndf) 
            IDK=IDK+1 
          END DO 
!          IDK=IDK+1 
          IF(IDK.GT.NDF) IDK=NDF 
          KRAY(IU,ID)=IDK 
          DRAY(IU,ID)=(RDX(ID)-DENSF(IDK-1))/(DENSF(IDK)-DENSF(IDK-1)) 
        END DO 
        IF(IU.LE.NREXT) THEN 
          DO ID=1,NUDF(IU) 
            JD=2*NUDF(IU)-ID 
            KRAY(IU,JD)=KRAY(IU,ID) 
            DRAY(IU,JD)=DRAY(IU,ID) 
          END DO 
        END IF 
      END DO 
! 
!    remaining rays (without finer grid) 
! 
      IF(NFIRY.LT.KMU) THEN 
        IU=KMU 
        KRAY(IU,1)=2 
        DRAY(IU,1)=0. 
        IDK=1 
        DO ID=2,NUDF(IU) 
          DO WHILE (DENS(ID).GE.DENSF(IDK).and.idk.le.ndf) 
            IDK=IDK+1 
          END DO 
!          IDK=IDK+1 
          IF(IDK.GT.NDF) IDK=NDF 
          KRAY(IU,ID)=IDK 
          DRAY(IU,ID)=(DENS(ID)-DENSF(IDK-1))/(DENSF(IDK)-DENSF(IDK-1)) 
        END DO 
        DO IU=NFIRY+1,KMU 
          DO ID=1,NUDF(IU) 
            KRAY(IU,ID)=KRAY(KMU,ID) 
            DRAY(IU,ID)=DRAY(KMU,ID) 
            DFRQF(IU,ID)=DFRQ(IU,ID) 
            DELZF(IU,ID)=DELZ(IU,ID) 
          ENDDO 
          IF(IU.LE.NREXT) THEN 
            DO ID=1,NUDF(IU) 
              JD=2*NUDF(IU)-ID 
              KRAY(IU,JD)=KRAY(IU,ID) 
              DRAY(IU,JD)=DRAY(IU,ID) 
              DFRQF(IU,JD)=-DFRQF(IU,ID) 
            END DO 
          END IF 
        END DO 
      END IF 
! 
      NFTOT=0 
      DO IU=2,KMU 
        IUD=NUDF(IU) 
        IF(IU.LE.NREXT) IUD=2*NUDF(IU)-1 
        NFTOT=NFTOT+IUD 
      ENDDO 
      write(10,*) 'NFTOT=',NFTOT 
! 
      RETURN 
      END SUBROUTINE SETRAY 
! 
! 
!     **************************************************************** 
! 
! 
      SUBROUTINE WGTJH1 
!     ================= 
! 
!     Angle quadrature weights 
!      from Hummer, Kunasz, & Kunasz, 1973, Comp. Phys. Comm. 6, 38 
! 
!     The present version of this routine assumes that there are 
!      impact rays tangent to every depth layers (i.e. NREXT=ND) 
! 
      use accura 
      use params 
      use wincom 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: UN=1., TWO=2.,HALF=0.5,SIX=6. 
      REAL(DP), PARAMETER :: C03=UN/3.,D03=2./3.,C04=UN/4.,C06=UN/6. 
      REAL(DP), PARAMETER :: C24=UN/24.,C45=UN/45.,D45=2./45.,C72=UN/72. 
      REAL(DP) :: WAJ(MKU),WBJ(MKU),AHH(MKU,4) 
      REAL(DP) :: BMUH(MKU),BMUHP(MKU),WAH(MKU),WBH(MKU) 
      REAL(DP) :: WSD(MKU),WSU(MKU),WSL(MKU),WUU(MKU) 
      REAL(DP) :: WTD(MKU),WTU(MKU),WTL(MKU) 
! 
      DEPTHS: DO ID=1,ND 
        DO IU=ID+1,KMU 
          AHH(IU,1)=BMU(IU,ID)-BMU(IU-1,ID) 
          AHH(IU,2)=AHH(IU,1)*AHH(IU,1) 
          AHH(IU,3)=AHH(IU,2)*AHH(IU,1) 
          AHH(IU,4)=AHH(IU,3)*AHH(IU,1) 
          BMUH(IU)=BMU(IU,ID)*AHH(IU,1) 
          BMUHP(IU)=BMU(IU-1,ID)*AHH(IU,1) 
        END DO 
! 
!     Weights for J 
! 
        WAJ(ID)=HALF*AHH(ID+1,1) 
        WAJ(KMU)=HALF*AHH(KMU,1) 
        WBJ(ID)=-C24*AHH(ID+1,3) 
        WBJ(KMU)=-C24*AHH(KMU,3) 
        WSL(ID+1)=C06*AHH(ID+1,1) 
        WSU(KMU-1)=0. 
        WSD(ID)=C03*AHH(ID+1,1) 
        WSD(KMU)=UN 
        WTL(ID+1)=UN/AHH(ID+1,1) 
        WTU(KMU-1)=0. 
        WTD(ID)=-WTL(ID+1) 
        WTD(KMU)=0. 
        DO IU=ID+1,KMU-1 
          WAJ(IU)=HALF*(AHH(IU,1)+AHH(IU+1,1)) 
          WBJ(IU)=-C24*(AHH(IU+1,3)+AHH(IU,3)) 
          AH1=SIX/(AHH(IU,1)+AHH(IU+1,1)) 
          WSL(IU+1)=C06*AH1*AHH(IU+1,1) 
          WSU(IU-1)=UN-WSL(IU+1) 
          WSD(IU)=TWO 
          WTL(IU+1)=AH1/AHH(IU+1,1) 
          WTU(IU-1)=AH1/AHH(IU,1) 
          WTD(IU)=-SIX/AHH(IU,1)/AHH(IU+1,1) 
        END DO 
        NMUD=KMU-ID+1 
        CALL TRIDAG(WSL,WSD,WSU,WBJ,WUU,NMUD) 
        WMUJ(ID,ID)=WAJ(ID)+WTD(ID)*WUU(ID)+WTU(ID)*WUU(ID+1) 
        WMUJ(KMU,ID)=WAJ(KMU)+WTL(KMU)*WUU(KMU-1)+WTD(KMU)*WUU(KMU) 
        DO IU=ID+1,KMU-1 
          WMUJ(IU,ID)=WAJ(IU)+WTL(IU)*WUU(IU-1)+                          & 
     &                WTD(IU)*WUU(IU)+WTU(IU)*WUU(IU+1) 
        END DO 
! 
!       Weights for emergent flux H 
! 
        IF(ID.GT.1) CYCLE DEPTHS 
        WAH(ID)=HALF*BMUH(ID+1)-C03*AHH(ID+1,2) 
        WAH(KMU)=HALF*BMUHP(KMU)+C03*AHH(KMU,2) 
        WBH(ID)=AHH(ID+1,3)*(C45*AHH(ID+1,1)-C24*BMU(ID+1,ID)) 
        WBH(KMU)=-AHH(KMU,3)*(C45*AHH(KMU,1)+C24*BMU(KMU-1,ID)) 
        WSL(ID+1)=0. 
        WSD(ID)=UN 
        WTL(ID+1)=0. 
        WTD(ID)=0. 
        DO IU=ID+1,KMU-1 
            WAH(IU)=HALF*(BMUH(IU+1)+BMUHP(IU))-                          & 
     &             C03*(AHH(IU+1,2)-AHH(IU,2)) 
            WBH(IU)=-C24*(BMUH(IU+1)*AHH(IU+1,2)+BMUHP(IU)*AHH(IU,2))+    & 
     &             C45*(AHH(IU+1,4)-AHH(IU,4)) 
        END DO 
        CALL TRIDAG(WSL,WSD,WSU,WBH,WUU,NMUD) 
        WMUH(ID)=WAH(ID)+WTD(ID)*WUU(ID)+WTU(ID)*WUU(ID+1) 
        WMUH(KMU)=WAH(KMU)+WTL(KMU)*WUU(KMU-1)+WTD(KMU)*WUU(KMU) 
        DO IU=ID+1,KMU-1 
          WMUH(IU)=WAH(IU)+WTL(IU)*WUU(IU-1)+                             & 
     &             WTD(IU)*WUU(IU)+WTU(IU)*WUU(IU+1) 
        END DO 
! 
      END DO DEPTHS 
! 
!     Weights for H are overwritten by trapezoidal weigths 
! 
      id=1 
      wmuh(1)=bmu(1,id)*(bmu(2,id)-bmu(1,id))*half 
      wmuh(kmu)=bmu(kmu,id)*(bmu(kmu,id)-bmu(kmu-1,id))*half 
      do iu=2,kmu-1 
        wmuh(iu)=bmu(iu,id)*(bmu(iu+1,id)-bmu(iu-1,id))*half 
      end do 
! 
      RETURN 
      END SUBROUTINE WGTJH1 
! 
! 
!     **************************************************************** 
! 
! 
      SUBROUTINE TRIDAG(A,B,C,R,U,N) 
!     ============================== 
! 
!     Solve tridiagonal system of equations 
!      from Numerical Recipes (standard Gaussian elimination) 
! 
      use accura 
      use params 
      use wincom 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: A(N),B(N),C(N),R(N),U(N) 
      REAL(DP) ::  GTRID(MKU) 
! 
      BTRID=B(1) 
      U(1)=R(1)/BTRID 
      DO J=2,N 
        GTRID(J)=C(J-1)/BTRID 
        BTRID=B(J)-A(J)*GTRID(J) 
        U(J)=(R(J)-A(J)*U(J-1))/BTRID 
      END DO 
      DO J=N-1,1,-1 
        U(J)=U(J)-GTRID(J+1)*U(J+1) 
      END DO 
! 
      RETURN 
      END 
! 
! 
!     **************************************************************** 
! 
! 
      SUBROUTINE RESOLW 
!     ================= 
! 
!     driver for evaluating opacities and emissivities which then 
!     enter the solution of the radiative transfer equation (RTEWIN) 
!     Setup opacities for a given frequency set 
!     Oversample in radial and frequency space for later interpolation 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use wincom 
      use lindat 
 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: UN=1., TWO=2., HALF=0.5 
      REAL(DP) :: ABSO(MOPAC),EMIS(MOPAC) 
      REAL(DP) :: ABSD(MDEPTH),ASF(MDEPF),XDS(MDEPTH),XDSF(MDEPF) 
! 
!     set up the partial line list for the current interval 
! 
      CALL INISET 
! 
!     output of information about selected lines 
! 
      IF(IMODE.LT.2) CALL INIBLA 
! 
!  Setup fine grid of frequencies 
! 
      CLV=UN/2.997925E10 
      FQ1=FREQ(1)*(UN+VINF*CLV) 
      FQ2=FREQ(NFREQ)*(UN-VINF*CLV) 
      VXD=SQRT(0.3e7*TSTD)*FREQ(1)*CLV 
      VXS=SPACE0*FREQ(1)*FREQ(1)*CLV*1.e-7 
!     DVX=MAX(VXD,VXS) 
      DVX=VXS 
      NOPAC=int((FQ1-FQ2)/DVX)+1 
      DVX=(FQ1-FQ2)/DFLOAT(NOPAC) 
      NOPAC=NOPAC+3 
      nopac=nfreq 
      WRITE(6,"(/,' Opacity table for',i5,' frequencies and',/,           & 
     &   '                  ',i5,' radial (density) points')")            & 
     &   NOPAC,NDF 
      IF(NOPAC.GT.MOPAC) CALL quit('Too many freqs in fine grid') 
      DO IJ=1,NOPAC 
         FFQ(ij)=FQ1-DFLOAT(ij-1)*DVX 
         fr=freq(ij)*1.d-15 
         BNUE(IJ)=BN*fr*fr*fr 
         DO IJCI=IJC,NFREQC-1 
            IF(WLAM(IJ).LE.WLAMC(IJCI)) EXIT 
         END DO 
         IJC=IJCI 
         IJCINT(IJ)=MAX(IJC-1,1) 
         IJCI=IJCINT(IJ) 
         FRX1(IJ)=(FREQ(IJ)-FREQC(IJCI+1))/                               & 
     &            (FREQC(IJCI)-FREQC(IJCI+1)) 
!         write(80,"(2i5,2f10.3,1p2e11.3)") 
!    *    ij,ijci,wlam(ij),wlamc(ijci),freq(ij),frx1(ij) 
      END DO 
      nfreq=nopac 
      DO JI=1,NOPAC-1 
        FFQV(JI)=UN/(FFQ(JI)-FFQ(JI+1)) 
      END DO 
      FFQV(NOPAC)=UN 
! 
!     the continuum opacities and radiation field - done only once 
! 
!     ----------------------------------- 
      if(iblank.le.1) then 
! 
!     determine the "core" radius and the factor that multiplies 
!     H_nu at ID=1 to get physical flux there (R2F) 
! 
      ID0=ND 
      DO WHILE(TEMP(ID0).GT.TEFF .AND. ID0.GT.1) 
        ID0=ID0-1 
      END DO 
      ID0=ID0+1 
      R2F=RD(1)*RD(1)/RD(ID0)/RD(ID0) 
! 
!     photoinization cross-sections 
! 
      CALL CROSEW 
! 
!     store opacity and emissivity in continuum 
! 
      DO ID=1,ND 
         CALL OPACW(ID,ABSO,EMIS,0) 
         DO IJ=1,NFREQC 
            CHC(IJ,ID)=ABSOC(IJ) / DENSCON(ID) 
            ETC(IJ,ID)=EMISC(IJ) / DENSCON(ID) 
            SCC(IJ,ID)=(SCATC(IJ)+ELEC(ID)*SIGE) / DENSCON(ID) 
         END DO 
      END DO 
! 
!     radiation field in the continuum 
! 
      call rtesca 
      do ij=1,nfreqc 
         write(17,"(F10.4,1PE15.5)") wlamc(ij),fluxc(ij)*r2f 
      end do 
! 
      end if 
!     ----------------------------------- 
! 
!     Store opacity and thermal source function in all frequencies 
!     and depths 
! 
      DO ID=1,ND 
         CALL OPACW(ID,ABSO,EMIS,1) 
         DO IJ=1,NOPAC 
            AB(IJ,ID)=ABSO(IJ) / DENSCON(ID) 
            STH(IJ,ID)=EMIS(IJ)/ABSO(IJ) 
         END DO 
      END DO 
! 
!  Interpolate to a finer radial (density) grid 
! 
      if(ndf.ne.nd)  then 
      DO ID=1,ND 
        XDS(ID)=LOG10(DENS(ID)) 
      END DO 
      DO ID=1,NDF 
        XDSF(ID)=LOG10(DENSF(ID)) 
      END DO 
      DO IJ=1,NOPAC 
        DO ID=1,ND 
          ABSD(ID)=AB(IJ,ID) 
        END DO 
        CALL INTERP(XDS ,ABSD,XDSF ,ASF,ND,NDF,2,0,1) 
        DO ID=1,NDF 
          AB(IJ,ID)=ASF(ID) 
        END DO 
        DO ID=1,ND 
          ABSD(ID)=STH(IJ,ID) 
        END DO 
        CALL INTERP(XDS ,ABSD,XDSF ,ASF,ND,NDF,2,0,1) 
        DO ID=1,NDF 
          STH(IJ,ID)=ASF(ID) 
        END DO 
      END DO 
      DO IJ=1,NFREQC 
        DO ID=1,ND 
          ABSD(ID)=SCC(IJ,ID) 
        END DO 
        CALL INTERP(XDS ,ABSD,XDSF ,ASF,ND,NDF,2,0,1) 
        DO ID=1,NDF 
          SCH(IJ,ID)=ASF(ID) 
        END DO 
      END DO 
      end if 
      write(6,*) ' Done' 
      write(6,*) 
! 
! 
!     Loop on rays, solving radiative transfer equation 
! 
      DO IJ=1,NFREQ 
        FLUX(IJ)=0. 
      END DO 
      DO IU=2,KMU 
        CALL RTEWIN(IU) 
      END DO 
      DO IJ=1,NFREQ 
        FLUX(IJ)=FLUX(IJ)*R2F 
      END DO 
! 
      RETURN 
      END SUBROUTINE RESOLW 
! 
! 
!     **************************************************************** 
! 
! 
      SUBROUTINE RTESCA 
!     ================= 
! 
!     Solution of the radiative transfer equation 
!      for deriving the scattering in continuum 
! 
!     Solution along every rays, for the spherically-symmetric case 
! 
!     Solution in the optical depth scale 
! 
!     The numerical method used: 
!     Discontinuous Finite Element method 
!     Castor, Dykema, Klein, 1992, ApJ 387, 561. 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use wincom 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: UN=1., TWO=2., HALF=0.5,DJMAX=1.e-3 
      INTEGER, PARAMETER  :: NTRALI=10 
      REAL(DP) :: ST0(mdepf ),RAD00(mdepf ),AB0(mdepf ),ALI1(mdepf ),     & 
     &            rip(mdepf ),rim(mdepf ),riin(mdepf ),riup(mdepf ),      & 
     &            aip(mdepf ),aim(mdepf ),aiin(mdepf ),aiup(mdepf ),      & 
     &            dt(mdepf ),dtau(mdepf ),RDX(mdepf ),PTX(mdepf ),        & 
     &            uf(mdepf ),af(mdepf ),ss0(mdepf ),scx(mdepth) 
      REAL(DP) :: densr(mdepf),rdy(mdepf),                                & 
     &            abc0(mdepf),abc1(mdepf),stc0(mdepf),stc1(mdepf),        & 
     &            scc0(mdepf),scc01(mdepf) 
! 
!     overall loop over continuum frequencies 
! 
      FREQLOOP: DO IJ=1,NFREQC 
         FR=FREQC(IJ) 
! 
!        Initialisation of J=B 
! 
         if(ij.eq.1) then 
            FR15=FR*1.e-15 
            BNU=BN*FR15*FR15*FR15 
            HKFR=HK*FR 
            DO ID=1,ND 
              RAD00(ID)=BNU/(EXP(HKFR/TEMP(ID))-UN) 
            END DO 
         end if 
! 
!        Loop over electron scattering 
! 
         ELSCAT: DO ITRALI=1,NTRALI 
            fluxc(ij)=0. 
! 
            DO ID=1,ND 
              RAD1(ID)=0. 
              ALI1(ID)=0. 
            END DO 
! 
!           Loop over impact rays 
! 
            if(nd.eq.ndf) then 
               do id=1,nd 
                  densf(id)=dens(id) 
                  rdx(id)=rad00(id) 
                  abc0(id)=chc(ij,id) 
                  stc0(id)=etc(ij,id)/chc(ij,id) 
                  scc0(id)=scc(ij,id) 
               end do 
             else 
               CALL INTERP(DENS,RAD00,DENSF,RDX,ND,NDF,4,1,0) 
               do id=1,nd 
                  abc1(id)=chc(ij,id) 
                  stc1(id)=etc(ij,id)/chc(ij,id) 
                  scc01(ij)=scc(ij,id) 
               end do 
               CALL INTERP(DENS,abc1,DENSF,abc0,ND,NDF,4,1,0) 
               CALL INTERP(DENS,stc1,DENSF,stc0,ND,NDF,4,1,0) 
               CALL INTERP(DENS,scc01,DENSF,scc0,ND,NDF,4,1,0) 
            end if 
 
            IMPACT: DO IU=1,KMU 
              iud=nud(iu) 
              IF(IU.LE.NFIRY) IUD=NUDF(IU) 
              if(iud.le.1) cycle impact 
              DO ID=1,IUD 
                KY=KRAY(IU,ID) 
                YDR=DRAY(IU,ID) 
                YDR1=UN-DRAY(IU,ID) 
                DENSR(ID)=YDR1*DENSF(KY-1)+YDR*DENSF(KY) 
                AB0(ID)=YDR1*abc0(KY-1)+YDR*abc0(KY) 
                ST0(ID)=YDR1*stc0(KY-1)+YDR*stc0(KY) 
                SC0=YDR1*scc0(KY-1)+YDR*scc0(KY) 
                RDY(id)=YDR1*RDX(KY-1)+YDR*RDX(KY) 
                SS0(ID)=SC0/AB0(ID) 
                ST0(ID)=ST0(ID)+SS0(ID)*RDY(ID) 
              END DO 
              IF(IU.LE.NFIRY) THEN 
                DO ID=1,IUD-1 
                  DTAU(ID)=HALF*(AB0(ID)+AB0(ID+1))*DELZF(IU,ID) 
                END DO 
               ELSE 
                DO ID=1,IUD-1 
                  DT(ID)=HALF*(AB0(ID)+AB0(ID+1)) 
                  DTAU(ID)=DT(ID)*DELZ(IU,ID) 
                END DO 
              END IF 
! 
!             incoming intensity   (TAUMIN=0.) 
! 
              rim(1)=0. 
              aim(1)=0. 
              do id=1,iud-1 
                dt0=dtau(id) 
                dtaup1=dt0+un 
                dtau2=dt0*dt0 
                bb=two*dtaup1 
                cc=dt0*dtaup1 
                aa=un/(dtau2+bb) 
                rip(id)=(bb*rim(id)+cc*st0(id)-dt0*st0(id+1))*aa 
                rim(id+1)=(two*rim(id)+dt0*st0(id)+cc*st0(id+1))*aa 
                aip(id)=(cc+bb*aim(id))*aa 
                aim(id+1)=cc*aa 
              end do 
              do id=2,iud-1 
                dtt=un/(dtau(id-1)+dtau(id)) 
                riin(id)=(rim(id)*dtau(id)+rip(id)*dtau(id-1))*dtt 
                aiin(id)=(aim(id)*dtau(id)+aip(id)*dtau(id-1))*dtt 
              end do 
              riin(1)=rim(1) 
              riin(iud)=rim(iud) 
              aiin(1)=aim(1) 
              aiin(iud)=aim(iud) 
              rip(iud)=rim(iud) 
! 
!             Outgoing intensity 
!             symmetric boundary condition (rim(iud)=riin(iud)) 
!             or diffusion approx. for core rays 
! 
              IF(IU.GT.NREXT) THEN 
                PLAND=BNU/(EXP(HK*FR/TEMP(ND))-UN) 
                DPLAN=PLAND-BNU/(EXP(HK*FR/TEMP(ND-1))-UN) 
                rip(iud)=PLAND+dplan/dtau(iud-1) 
                dt0=dtau(iud-1) 
                dtaup1=dt0+un 
                dtau2=dt0*dt0 
                bb=two*dtaup1 
                cc=dt0*dtaup1 
                aa=dtau2+bb 
                rim(iud)=(aa*rip(iud)-cc*st0(iud)+dt0*st0(iud-1))/bb 
              END IF 
              do id=iud-1,1,-1 
                dt0=dtau(id) 
                dtaup1=dt0+un 
                dtau2=dt0*dt0 
                bb=two*dtaup1 
                cc=dt0*dtaup1 
                aa=un/(dtau2+bb) 
                rip(id+1)=(bb*rim(id+1)+cc*st0(id+1)-dt0*st0(id))*aa 
                rim(id)=(two*rim(id+1)+dt0*st0(id+1)+cc*st0(id))*aa 
                aip(id+1)=(cc+bb*aim(id+1))*aa 
                aim(id)=cc*aa 
              end do 
              do id=2,iud-1 
                dtt=un/(dtau(id-1)+dtau(id)) 
                riup(id)=(rim(id)*dtau(id-1)+rip(id)*dtau(id))*dtt 
                aiup(id)=(aim(id)*dtau(id-1)+aip(id)*dtau(id))*dtt 
              end do 
              riup(1)=rim(1) 
              riup(iud)=rim(iud) 
              aiup(1)=aim(1) 
              aiup(iud)=aim(iud) 
! 
!             symmetrized (Feautrier) intensity  -- (riin+riup)/2 -- 
!             and interpolation in original radial grid 
! 
              do id=1,iud 
                uf(id)=(riup(id)+riin(id)) 
                af(id)=(aiup(id)+aiin(id)) 
              end do 
              if(iu.le.nfiry) then 
                inrp=min(nud(iu),4) 
                call interp(densr,uf,dens,ptx,iud,nud(iu),inrp,1,0) 
                do id=1,nud(iu) 
                  uf(id)=ptx(id) 
                end do 
                call interp(densr,af,dens,ptx,iud,nud(iu),inrp,1,0) 
                do id=1,nud(iu) 
                  af(id)=ptx(id) 
                end do 
                iud=nud(iu) 
              end if 
! 
!             Contribution to J 
! 
              do id=1,nud(iu) 
                rad1(id)=rad1(id)+wmuj(iu,id)*uf(id) 
                ali1(id)=ali1(id)+wmuj(iu,id)*af(id) 
              end do 
              FLUXc(IJ)=FLUXc(IJ)+WMUH(IU)*RIM(1) 
! 
!           End loop over impact rays 
! 
            END DO IMPACT 
! 
!           solution of the transfer equation 
!           Variables: 
!           RAD1    - mean intensity 
! 
            NDX=NUDF(KMU) 
            CALL INTERP(DENSR,SS0,DENS,SCX,NDX,ND,4,1,1) 
            DJTOT=0. 
            DO ID=1,ND 
              RAD1(ID)=RAD1(ID)*HALF 
              ALI1(ID)=ALI1(ID)*HALF 
              SSS=SCX(ID) 
              DELTAJ=(RAD1(ID)-RAD00(ID))/(UN-SSS*ALI1(ID)) 
              RAD00(ID)=RAD00(ID)+DELTAJ 
              DJTOT=MAX(DJTOT,ABS(DELTAJ/RAD00(ID))) 
            END DO 
            write(6,"(' IJ,LAM,ITRALI,DJ',i5,f10.2,i5,1p2e12.3)")         & 
     &          ij,2.997925e18/fr,itrali,djtot,djmax 
            IF(DJTOT.GT.DJMAX.AND.ITRALI.LE.NTRALI) THEN 
               CYCLE ELSCAT 
              ELSE 
               EXIT ELSCAT 
            END IF 
! 
!        end loop for electron scattering 
! 
         END DO ELSCAT 
 
         CALL INTERP(DENS,RAD00,DENSF,RDX,ND,NDF,4,1,0) 
         do id=1,ndf 
           sccf(ij,id)=scc0(ID)*RDX(ID) 
         enddo 
         fluxc(ij)=fluxc(ij)*2.997925e18/wlamc(ij)**2*0.5 
! 
      END DO FREQLOOP 
      RETURN 
      END SUBROUTINE RTESCA 
! 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE RTEWIN(IU) 
!     ===================== 
! 
!     Solution of the radiative transfer equation - frequency by 
!     frequency - for the known source function. 
! 
!     The numerical method used: 
!     Discontinuous Finite Element (DFE) method 
!     Castor, Dykema, Klein, 1992, ApJ 387, 561. 
! 
!     Input through blank COMMON block: 
!      AB     - two-dimensional array  absorption coefficient (frequency, 
!               depth) 
!      STH     - Thermal source function 
! 
!     Version including velocity field and extension 
!      radiative transfer along ray IU 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use wincom 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: UN=1., TWO=2., HALF=0.5 
      REAL(DP), PARAMETER :: TAUREF = 0.6666666666667 
      REAL(DP) :: ST0(2*MDEPF ),TAU(2*MDEPF ),AB0(2*MDEPF ),              & 
     &            rip(2*MDEPF ),rim(2*MDEPF ),                            & 
     &            sctd(2*mdepf) 
! 
      IUD=NUDF(IU) 
      IF(IU.LE.NREXT) IUD=2*NUDF(IU)-1 
      IF(IUD.EQ.1) RETURN 
      IF(NFREQ.GT.1) dlama0=(wlobs(nfrobs)-wlobs(1))/(nfrobs-1) 
! 
!     overall loop over frequencies (observer's frame) 
! 
      FREQLOOP: DO IJ=1,NFROBS 
         FR=FRQOBS(IJ) 
         wl0=wlobs(ij) 
! 
!       Opacity and total source function 
!       interpolation in opacity table 
! 
         IVK=NOPAC-2 
         DO ID=1,IUD 
            KY=KRAY(IU,ID) 
            YDR=DRAY(IU,ID) 
            YDR1=UN-YDR 
            dwlcom=wl0*DFRQF(IU,ID) 
            wlcom=wl0+dwlcom 
            if(wlcom.le.wlam(3)) then 
               abd1=ab(1,ky-1) 
               std1=sth(1,ky-1) 
               abd0=ab(1,ky) 
               std0=sth(1,ky) 
               ij1=1 
             else if(wlcom.ge.wlam(nfreq)) then 
               abd1=ab(nfreq,ky-1) 
               std1=sth(nfreq,ky-1) 
               abd0=ab(nfreq,ky) 
               std0=sth(nfreq,ky) 
               ij1=nfreq 
             else 
               xijap=(wlcom-wlam(3))/dlama0 
               ijap=int(xijap) 
               ijap=max(ijap,1) 
               ijap=min(ijap,nfreq) 
               wlap=wlam(ijap) 
               if(wlcom.lt.wlap) then 
                  ij1=ijap-1 
                  do iji=ijap-1,1,-1 
                     if(wlcom.ge.wlam(iji)) exit 
                  end do 
                  ij1=iji 
                else 
                  ij1=ijap+1 
                  do iji=ijap+1,nfreq 
                     if(wlcom.lt.wlam(iji)) exit 
                  end do 
                  ij1=iji-1 
               end if 
               xfa=(wlam(ij1+1)-wlcom)/(wlam(ij1+1)-wlam(ij1)) 
               abd1=xfa*ab(ij1,ky-1)+(1.-xfa)*ab(ij1+1,ky-1) 
               std1=xfa*sth(ij1,ky-1)+(1.-xfa)*sth(ij1+1,ky-1) 
               abd0=xfa*ab(ij1,ky)+(1.-xfa)*ab(ij1+1,ky) 
               std0=xfa*sth(ij1,ky)+(1.-xfa)*sth(ij1+1,ky) 
            end if 
            AB0(ID)=YDR1*Abd1+YDR*abd0 
            ST0(ID)=YDR1*Std1+YDR*Std0 
! 
!          Add scattering 
! 
            IJC=IJCINT(IJ1) 
            IF(IFREQ.NE.17) THEN 
               SC1=YDR1*SCCF(ijc,KY-1)+YDR*SCCF(ijc,KY) 
               SC2=YDR1*SCCF(ijc+1,KY-1)+YDR*SCCF(ijc+1,KY) 
               SCT=FRX1(ij1)*SC1+(1.-FRX1(ij1))*SC2 
               sctd(id)=sct/ab0(id) 
               ST0(ID)=ST0(ID)+SCT/AB0(ID) 
            END IF 
         END DO 
! 
!       Optical depth scale 
! 
         TAU(1)=0. 
         IREF=1 
         IF(IU.LE.NFIRY) THEN 
            DO ID=1,IUD-1 
               JD=ID 
               IF(ID.GT.NUDF(IU)) JD=2*NUDF(IU)-ID-1 
               DT=HALF*(AB0(ID+1)+AB0(ID))*DELZF(IU,JD) 
               TAU(ID+1)=TAU(ID)+DT 
            END DO 
          ELSE 
            DO ID=1,IUD-1 
               JD=ID 
               IF(ID.GT.NUD(IU)) JD=2*NUD(IU)-ID-1 
               DT=HALF*(AB0(ID+1)+AB0(ID))*DELZ(IU,JD) 
               TAU(ID+1)=TAU(ID)+DT 
            END DO 
         END IF 
         if(iu.eq.kmu) then 
            DO ID=1,IUD-1 
               IF(TAU(ID).LE.TAUREF.AND.TAU(ID+1).GT.TAUREF) IREF=ID 
            END DO 
            irefd(ij)=iref 
         end if 
! 
!        Outgoing intensity 
! 
         IF(IU.LE.NREXT) THEN 
! 
!       1. External rays 
! 
           ndt=iud 
           rip(ndt)=0. 
           dt0=tau(ndt)-tau(ndt-1) 
           dtaup1=dt0+un 
           dtau2=dt0*dt0 
           bb=two*dtaup1 
           cc=dt0*dtaup1 
           aa=dtau2+bb 
           rim(ndt)=(aa*rip(ndt)-cc*st0(ndt)+dt0*st0(ndt-1))/bb 
           do id=1,iud-1 
             jd=iud-id 
             dt0=tau(jd+1)-tau(jd) 
             dtaup1=dt0+un 
             dtau2=dt0*dt0 
             bb=two*dtaup1 
             cc=dt0*dtaup1 
             aa=un/(dtau2+bb) 
             rim(jd)=(two*rim(jd+1)+dt0*st0(jd+1)+cc*st0(jd))*aa 
           enddo 
         ELSE 
! 
!        2. core rays 
! 
           NDT=IUD 
           FR15=FR*1.D-15 
           BNU=BN*FR15*FR15*FR15 
           PLAND=BNU/(EXP(HK*FR/TEMP(ND))-UN) 
           DPLAN=BNU/(EXP(HK*FR/TEMP(ND-1))-UN) 
           DPLAN=(PLAND-DPLAN)/(TAU(IUD)-TAU(IUD-1)) 
           RIP(NDT)=PLAND+DPLAN 
           dt0=tau(ndt)-tau(ndt-1) 
           dtaup1=dt0+un 
           dtau2=dt0*dt0 
           bb=two*dtaup1 
           cc=dt0*dtaup1 
           aa=dtau2+bb 
           rim(ndt)=(aa*rip(ndt)-cc*st0(ndt)+dt0*st0(ndt-1))/bb 
           do id=iud-1,1,-1 
              dt0=tau(id+1)-tau(id) 
              dtaup1=dt0+un 
              dtau2=dt0*dt0 
              bb=two*dtaup1 
              cc=dt0*dtaup1 
              aa=un/(dtau2+bb) 
              rim(id)=(two*rim(id+1)+dt0*st0(id+1)+cc*st0(id))*aa 
           enddo 
         ENDIF 
         FLUX(IJ)=FLUX(IJ)+WMUH(IU)*RIM(1) 
! 
      END DO FREQLOOP 
      RETURN 
      END SUBROUTINE RTEWIN 
! 
! 
! *********************************************************************** 
! 
! 
      SUBROUTINE VELSET 
!     ================= 
! 
!     Determination of the macroscopic velocity as a function of depth 
! 
!     Input: 
! 
!     RSTAR   - stellar radius (in solar radii or in cm) 
!     RMAX    - maximum radial extent (in stellar radii) 
!     AMLOSS  - mass loss rate ( in solar masses per year) 
!     VELMAX  - maximum velocity (= V_infinity) - in km/s 
!     BETA    - beta exponent in the beta-law for velocity 
!     NDRAD   - Number of layers 
!     NRCORE  - Number of core rays 
! 
      use accura 
      use params 
      use modelp 
      use wincom 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      real(dp) :: zz(mdepth),vel0(mdepth),rrel(mdepth),                   & 
     &            den0(mdepth),vel00(mdepth),ind(mdepth),                 & 
     &            densa(mdepth),eleca(mdepth),tempa(mdepth),              & 
     &            rda(mdepth),rrela(mdepth),vel0a(mdepth) 
! 
      un=1 
      two=2. 

      velmax=3.e5
      nltoff=0
      iemoff=0
      itrad=0
      do id=1,nd
        wdil(id)=un
      end do
         read(55,*,iostat=ios2) velmax,ITRAD,nltoff,iemoff
         if(ios2.ne.0) then
            write(6,"(//' velmax (velocity for line rejection)',          &
     &       ' itrad,nltoff,iemoff',f10.1,2i3)")                          &
     &       velmax,itrad,nltoff,iemoff
         end if
      velmax=velmax*1.e5
      do id=1,nd
         ilvi(id)=0
         ilne(id)=0
         if(vel(id).gt.velmax.and.iemoff.eq.0) ilvi(id)=1
         if(vel(id).gt.velmax.and.nltoff.gt.0.and.iemoff.gt.0)            &
     &      ilne(id)=1
      end do


      read(55,*,iostat=ios) rstar,rmax,amloss,vinf,beta,                  & 
     &           ndrad,nrcore,nfiry,ndf,nda 
      if(ios.ne.0) return 
      rstr=rstar 
      if(rstar.lt.1.e5) rstr=rstar*6.9598e10 
      amdot=amloss*6.3029e25 
      RCORE=RSTR 
      XMDOT=amdot 
      BETAV=beta 
      con=amdot/12.566e5 
      conr=con/rstr/rstr 
      nrext0=ndrad-nd 
      zz(nd+nrext0)=0. 
      rd(nd+nrext0)=rstr 
      rrel(nd+nrext0)=1. 
      do iid=1,nd-1 
         id=nd-iid 
         zz(id+nrext0)=zz(id+1+nrext0)+2.*(dm(id+1)-dm(id))/              & 
     &                 (dens(id+1)+dens(id)) 
         rd(id+nrext0)=rstr+zz(id+nrext0) 
         rrel(id+nrext0)=rd(id+nrext0)/rstr 
      end do 
! 
      do id=1+nrext0,nd+nrext0 
         vel0(id)=con/rd(id)**2/dens(id-nrext0) 
         vel00(id)=vel0(id) 
         if(vel00(id).gt.vinf) vel00(id)=vinf 
      end do 
      vin=vel0(nrext0+1) 
      r1=rrel(nrext0+1) 
! 
      if(rrel(1+nrext0).lt.rmax.and.nd.lt.ndrad) then 
      rl1=1.-1./rrel(1+nrext0) 
      rl2=1.-1./rmax 
      drl=(rl2-rl1)/nrext0 
      do id=1,nrext0 
         rlo=rl2-(id-1)*drl 
         rrel(id)=1./(1.-rlo) 
         rd(id)=rrel(id)*rstr 
      end do 
      end if 
! 
      depth: do id=nd+nrext0-1,nrext0+1,-1 
         r0=rrel(id) 
         numid=0 
         do id1=nd+nrext0-1,nrext0+1,-1 
            x=un-r0/rrel(id1) 
            if(x.lt.1.e-6) x=1.e-6 
            v2=vinf*x**beta 
            ind(id1)=0 
            if(v2.ge.vel0(id1)) then 
              ind(id1)=id1 
              numid=numid+1 
            end if 
         end do 
         if(numid.eq.0) exit depth 
         rsum=0. 
         isum=0 
         do id1=nd+nrext0-1,nrext0+1,-1 
            if(ind(id1).gt.0) then 
               rsum=rsum+rrel(id1) 
               isum=isum+id1 
            endif 
         end do 
         rc=rsum/numid 
         idc=isum/numid 
         numid0=numid 
         r00=r0 
      end do depth 
      v1=vel0(idc) 
      r0=(r0+r00)*0.5 
      if(r0.lt.rc) v2=vinf*(un-r0/rc)**beta 
      write(6,"('numid,idc,rc,r0,v1,v2 ',2i4,4f10.5)")                    & 
     & numid0,idc,rc,r0,v1,v2 
! 
      do id=nd+nrext0-1,1,-1 
         if(rrel(id).gt.rc.and.rrel(id).gt.r0)                            & 
     &      vel0(id)=vinf*(1.-r0/rrel(id))**beta 
      end do 
! 
      t1=temp(1) 
      erel=elec(1)/dens(1) 
      do id=nd,1,-1 
         temp(id+nrext0)=temp(id) 
         den0(id+nrext0)=dens(id) 
         elec(id+nrext0)=elec(id) 
         do i=1,nlevel 
            popul(i,id+nrext0)=popul(i,id) 
         end do 
         WMM(ID+nrext0)=WMM(id) 
         WMY(ID+nrext0)=WMY(id) 
         YTOT(ID+nrext0)=YTOT(id) 
         do i=1,natom 
            relab(i,id+nrext0)=relab(i,id) 
            abund(i,id+nrext0)=abund(i,id) 
         end do 
         do i=1,matom 
            abndd(i,id+nrext0)=abndd(i,id) 
         end do 
      end do 
! 
      do id=1,nrext0 
         TEMP(ID)=T1 
         WMM(ID)=WMM(NREXT0+1) 
         WMY(ID)=WMY(NREXT0+1) 
         YTOT(ID)=YTOT(NREXT0+1) 
         do i=1,natom 
            relab(i,id)=relab(i,nrext0+1) 
            abund(i,id)=abund(i,nrext0+1) 
         end do 
         do i=1,matom 
            abndd(i,id)=abndd(i,nrext0+1) 
         end do 
      end do 
      idstd=idstd+nrext0 
! 
      VINF=vinf*1.e5 
      write(6,"('    ID    M      TEMP       ELEC      DENS       ',      & 
     &       '      R       Rrel     VEL'/)") 
      do id=1,nd+nrext0 
         if(vel0(id).gt.0.) dens(id)=con/rd(id)**2/vel0(id) 
         VEL(ID)=vel0(id)*1.e5 
      end do 
! 
      do id=nd,1,-1 
         id1=id+nrext0 
         elec(id1)=elec(id1)*dens(id1)/den0(id1) 
         do i=1,nlevel 
            popul(i,id1)=popul(i,id1)*dens(id1)/den0(id1) 
         end do 
      end do 
! 
      do id=1,nrext0 
         elec(id)=elec(nrext0+1)*dens(id)/dens(nrext0+1) 
         do i=1,nlevel 
            popul(i,id)=popul(i,nrext0+1)*dens(id)/dens(nrext0+1) 
         end do 
      end do 
! 
      ND=NDRAD 
      if(ndf.eq.0) ndf=nd 
      do id=1,nd 
         write(6,"(i3,1pe10.3,0pf8.0,1p3e12.3,0pf10.4,0p2f8.2)")          & 
     &    id,dm(id),temp(id),elec(id),dens(id),rd(id),                    & 
     &    rrel(id),vel0(id) 
         write(96,"(i3,1pe10.3,0pf8.0,1p3e12.3,0pf10.4,0p2f8.2)")         & 
     &    id,dm(id),temp(id),elec(id),dens(id),rd(id),rrel(id),           & 
     &    vel0(id),vel00(id) 
      end do 
! 
      if(nda.gt.0) then 
         XR1=LOG(DENS(1)) 
         XR2=LOG(DENS(ND)) 
         DXR=(XR2-XR1)/FLOAT(NDA-1) 
         DO ID=1,NDA 
            DENSA(ID)=EXP(XR1+FLOAT(ID-1)*DXR) 
         END DO 
         CALL INTERP(DENS,TEMP,DENSA,TEMPA,ND,NDA,3,1,1) 
         CALL INTERP(DENS,ELEC,DENSA,ELECA,ND,NDA,3,1,1) 
         CALL INTERP(DENS,RD,DENSA,RDA,ND,NDA,3,1,1) 
         CALL INTERP(DENS,RREl,DENSA,RRELA,ND,NDA,3,1,1) 
         CALL INTERP(DENS,VEL0,DENSA,VEL0A,ND,NDA,3,1,1) 
         do id=1,nda 
         write(6,"(i3,0pf8.0,1p3e12.3,0pf10.4,0p2f8.2)")                  & 
     &   id,tempa(id),eleca(id),densa(id),rda(id),rrela(id),vel0a(id) 
         write(96,"(i3,0pf8.0,1p3e12.3,0pf10.4,0p2f8.2)")                 & 
     &   id,tempa(id),eleca(id),densa(id),rda(id),rrela(id),vel0a(id) 
         end do 
      end if 
! 
      return 
      end SUBROUTINE VELSET 
! 
! 
! *********************************************************************** 
! 
! 
 
      SUBROUTINE RADTEM 
!     ================= 
! 
!     determination of the radiation temperatures 
!     after Schmutz (1991); inversion done by Newton-Raphson 
! 
      use accura 
      use params 
      use modelp 
      use wincom 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP), PARAMETER :: CON=2.0706e-16, un=1. 
      integer, parameter  :: nterad=3 
! 
      DO ID=1,ND 
         rx=RD(ND)/RD(ID) 
!        WDIL(ID)=0.5*(1.-sqrt(1.-rx*rx)) 
         wdil(id)=un-sqrt(un-rx*rx) 
      END DO 
      DO ITRD=1,NTERAD 
         if(itrad.eq.0) then 
            do id=1,nd 
               trad(itrd,id)=temp(id) 
            end do 
          else 
            II=0 
            JJ=0 
            IF(ITRD.LE.NION) II=NFIRST(ITRD) 
            IF(ITRD.LE.NION) JJ=NNEXT(ITRD) 
            DO ID=1,ND 
               TRAD(ITRD,ID)=TEMP(ID) 
               IF(II.GT.0) THEN 
                  AA=POPUL(JJ,ID)/POPUL(II,ID)*ELEC(ID)*CON 
                  AA=AA*G(II)/G(JJ)/WDIL(ID)/SQRT(TEMP(ID)) 
                  TR=TEMP(ID) 
                  ITER=0 
                  ITERATE: DO 
                     ITER=ITER+1 
                     XX=ENION(II)/BOLK/TR 
                     DTR=(AA*EXP(XX)-TR)/(1.+XX) 
                     DTRR=DTR/TR 
                     TR=TR+DTR 
                     IF(ABS(DTRR).GT.1.E-3.AND.ITER.LT.100) THEN 
                        CYCLE ITERATE 
                      ELSE 
                        EXIT ITERATE 
                     END IF 
                  END DO ITERATE 
                  TRAD(ITRD,ID)=TR 
               END IF 
            END DO 
         end if 
      END DO 
      write(6,"(/' radiation temperatures/')") 
      do id=1,nd 
         write(6,"(i5,4f10.1)") id,temp(id),trad(1,id),trad(2,id),        & 
     &                          trad(3,id) 
      end do 
      RETURN 
      END SUBROUTINE RADTEM 
! 
! 
! *********************************************************************** 
! 
! 
      FUNCTION SBFCH(FR,T) 
!     ==================== 
! 
!     cross-section times partition function for CH 
! 
!     from Kurucz ATLAS9 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      real(dp), parameter :: fihu=500.,fihui=1./fihu,                     & 
     &                       twhu=200.,twhui=1./twhu,                     & 
     &                       tenl=2.30258509299405E0 
! 
      REAL(DP) :: CROSSCH(15,105),PARTCH(41),CROSSCHT(15) 
      REAL(DP) :: C1(150),C2(150),C3(150),C4(150),C5(150) 
      REAL(DP) :: C6(150),C7(150),C8(150),C9(150),C10(150) 
      REAL(DP) :: C11(75) 
! 
      EQUIVALENCE (CROSSCH(1, 1),C1(1)),(CROSSCH(1,11),C2(1)) 
      EQUIVALENCE (CROSSCH(1,21),C3(1)),(CROSSCH(1,31),C4(1)) 
      EQUIVALENCE (CROSSCH(1,41),C5(1)),(CROSSCH(1,51),C6(1)) 
      EQUIVALENCE (CROSSCH(1,61),C7(1)),(CROSSCH(1,71),C8(1)) 
      EQUIVALENCE (CROSSCH(1,81),C9(1)),(CROSSCH(1,91),C10(1)) 
      EQUIVALENCE (CROSSCH(1,101),C11(1)) 
! 
      DATA C1/-38.000,-38.000,-38.000,-38.000,-38.000,-38.000,-38.000,    & 
     &-38.000,-38.000,-38.000,-38.000,-38.000,-38.000,-38.000,-38.000,    & 
     &        -32.727,-31.151,-30.133,-29.432,-28.925,-28.547,-28.257,    & 
     &-28.030,-27.848,-27.701,-27.580,-27.479,-27.395,-27.322,-27.261,    & 
     &        -31.588,-30.011,-28.993,-28.290,-27.784,-27.405,-27.115,    & 
     &-26.887,-26.705,-26.558,-26.437,-26.336,-26.251,-26.179,-26.117,    & 
     &        -30.407,-28.830,-27.811,-27.108,-26.601,-26.223,-25.932,    & 
     &-25.705,-25.523,-25.376,-25.255,-25.154,-25.069,-24.997,-24.935,    & 
     &        -29.513,-27.937,-26.920,-26.218,-25.712,-25.334,-25.043,    & 
     &-24.816,-24.635,-24.487,-24.366,-24.266,-24.181,-24.109,-24.047,    & 
     &        -28.910,-27.341,-26.327,-25.628,-25.123,-24.746,-24.457,    & 
     &-24.230,-24.049,-23.902,-23.782,-23.681,-23.597,-23.525,-23.464,    & 
     &        -28.517,-26.961,-25.955,-25.261,-24.760,-24.385,-24.098,    & 
     &-23.873,-23.694,-23.548,-23.429,-23.329,-23.245,-23.174,-23.113,    & 
     &        -28.213,-26.675,-25.680,-24.993,-24.497,-24.127,-23.843,    & 
     &-23.620,-23.443,-23.299,-23.181,-23.082,-22.999,-22.929,-22.869,    & 
     &        -27.942,-26.427,-25.446,-24.769,-24.280,-23.915,-23.635,    & 
     &-23.416,-23.241,-23.100,-22.983,-22.887,-22.805,-22.736,-22.677,    & 
     &        -27.706,-26.210,-25.241,-24.572,-24.088,-23.728,-23.451,    & 
     &-23.235,-23.063,-22.923,-22.808,-22.713,-22.633,-22.565,-22.507/ 
      DATA C2/-27.475,-26.000,-25.043,-24.382,-23.905,-23.548,-23.275,    & 
     &-23.062,-22.891,-22.753,-22.640,-22.546,-22.467,-22.400,-22.343,    & 
     &        -27.221,-25.783,-24.844,-24.193,-23.723,-23.372,-23.102,    & 
     &-22.892,-22.724,-22.588,-22.476,-22.384,-22.306,-22.240,-22.184,    & 
     &        -26.863,-25.506,-24.607,-23.979,-23.523,-23.182,-22.919,    & 
     &-22.714,-22.550,-22.417,-22.309,-22.218,-22.142,-22.078,-22.023,    & 
     &        -26.685,-25.347,-24.457,-23.835,-23.382,-23.044,-22.784,    & 
     &-22.580,-22.418,-22.286,-22.178,-22.089,-22.014,-21.950,-21.896,    & 
     &        -26.085,-24.903,-24.105,-23.538,-23.120,-22.805,-22.561,    & 
     &-22.370,-22.217,-22.093,-21.991,-21.906,-21.835,-21.775,-21.723,    & 
     &        -25.902,-24.727,-23.936,-23.376,-22.964,-22.654,-22.415,    & 
     &-22.227,-22.076,-21.955,-21.855,-21.772,-21.702,-21.644,-21.593,    & 
     &        -25.215,-24.196,-23.510,-23.019,-22.655,-22.378,-22.163,    & 
     &-21.992,-21.855,-21.744,-21.653,-21.577,-21.513,-21.459,-21.412,    & 
     &        -24.914,-23.937,-23.284,-22.820,-22.475,-22.212,-22.007,    & 
     &-21.845,-21.715,-21.609,-21.522,-21.449,-21.388,-21.336,-21.292,    & 
     &        -24.519,-23.637,-23.039,-22.606,-22.281,-22.030,-21.834,    & 
     &-21.678,-21.552,-21.450,-21.365,-21.295,-21.236,-21.185,-21.142,    & 
     &        -24.086,-23.222,-22.650,-22.246,-21.948,-21.722,-21.546,    & 
     &-21.407,-21.296,-21.205,-21.131,-21.070,-21.018,-20.974,-20.937/ 
      DATA C3/-23.850,-23.018,-22.472,-22.088,-21.805,-21.590,-21.422,    & 
     &-21.289,-21.182,-21.095,-21.024,-20.964,-20.914,-20.872,-20.835,    & 
     &        -23.136,-22.445,-21.994,-21.676,-21.440,-21.259,-21.117,    & 
     &-21.004,-20.912,-20.837,-20.775,-20.723,-20.679,-20.642,-20.611,    & 
     &        -23.199,-22.433,-21.927,-21.573,-21.314,-21.119,-20.969,    & 
     &-20.851,-20.758,-20.682,-20.621,-20.571,-20.529,-20.493,-20.463,    & 
     &        -22.696,-22.020,-21.585,-21.286,-21.071,-20.912,-20.791,    & 
     &-20.697,-20.622,-20.563,-20.514,-20.475,-20.442,-20.414,-20.391,    & 
     &        -22.119,-21.557,-21.194,-20.943,-20.761,-20.624,-20.518,    & 
     &-20.434,-20.367,-20.313,-20.268,-20.231,-20.201,-20.175,-20.153,    & 
     &        -21.855,-21.300,-20.931,-20.673,-20.485,-20.344,-20.235,    & 
     &-20.151,-20.084,-20.031,-19.988,-19.953,-19.924,-19.900,-19.880,    & 
     &        -21.126,-20.673,-20.382,-20.184,-20.044,-19.943,-19.868,    & 
     &-19.811,-19.769,-19.736,-19.710,-19.690,-19.674,-19.662,-19.652,    & 
     &        -20.502,-20.150,-19.922,-19.766,-19.657,-19.578,-19.520,    & 
     &-19.478,-19.446,-19.422,-19.404,-19.390,-19.379,-19.371,-19.365,    & 
     &        -20.030,-19.724,-19.530,-19.399,-19.309,-19.245,-19.199,    & 
     &-19.166,-19.142,-19.125,-19.112,-19.103,-19.096,-19.091,-19.088,    & 
     &        -19.640,-19.364,-19.189,-19.074,-18.996,-18.943,-18.906,    & 
     &-18.881,-18.863,-18.852,-18.844,-18.839,-18.837,-18.836,-18.836/ 
      DATA C4/-19.333,-19.092,-18.939,-18.838,-18.770,-18.725,-18.695,    & 
     &-18.675,-18.662,-18.655,-18.651,-18.649,-18.649,-18.651,-18.653,    & 
     &        -19.070,-18.880,-18.756,-18.674,-18.621,-18.585,-18.562,    & 
     &-18.548,-18.540,-18.536,-18.536,-18.537,-18.539,-18.542,-18.546,    & 
     &        -18.851,-18.708,-18.617,-18.558,-18.521,-18.498,-18.484,    & 
     &-18.477,-18.475,-18.476,-18.478,-18.482,-18.487,-18.493,-18.498,    & 
     &        -18.709,-18.599,-18.533,-18.494,-18.471,-18.459,-18.454,    & 
     &-18.454,-18.457,-18.462,-18.469,-18.476,-18.483,-18.490,-18.498,    & 
     &        -18.656,-18.572,-18.524,-18.497,-18.485,-18.480,-18.482,    & 
     &-18.486,-18.493,-18.501,-18.510,-18.519,-18.527,-18.536,-18.544,    & 
     &        -18.670,-18.613,-18.582,-18.566,-18.561,-18.562,-18.568,    & 
     &-18.575,-18.583,-18.592,-18.601,-18.610,-18.619,-18.627,-18.635,    & 
     &        -18.728,-18.700,-18.687,-18.683,-18.685,-18.691,-18.698,    & 
     &-18.706,-18.715,-18.723,-18.731,-18.739,-18.745,-18.752,-18.758,    & 
     &        -18.839,-18.835,-18.836,-18.842,-18.849,-18.857,-18.865,    & 
     &-18.872,-18.878,-18.883,-18.888,-18.892,-18.895,-18.898,-18.900,    & 
     &        -19.034,-19.041,-19.049,-19.057,-19.064,-19.069,-19.071,    & 
     &-19.071,-19.070,-19.068,-19.065,-19.061,-19.058,-19.054,-19.051,    & 
     &        -19.372,-19.378,-19.382,-19.380,-19.372,-19.359,-19.341,    & 
     &-19.321,-19.300,-19.280,-19.261,-19.243,-19.227,-19.212,-19.199/ 
      DATA C5/-19.780,-19.777,-19.763,-19.732,-19.686,-19.631,-19.573,    & 
     &-19.517,-19.465,-19.419,-19.379,-19.344,-19.314,-19.288,-19.265,    & 
     &        -20.151,-20.133,-20.087,-20.009,-19.911,-19.810,-19.715,    & 
     &-19.631,-19.559,-19.497,-19.446,-19.402,-19.365,-19.333,-19.306,    & 
     &        -20.525,-20.454,-20.312,-20.138,-19.970,-19.825,-19.705,    & 
     &-19.607,-19.528,-19.464,-19.411,-19.367,-19.330,-19.300,-19.274,    & 
     &        -20.869,-20.655,-20.366,-20.104,-19.894,-19.731,-19.604,    & 
     &-19.505,-19.426,-19.363,-19.312,-19.271,-19.236,-19.208,-19.184,    & 
     &        -21.179,-20.768,-20.380,-20.081,-19.856,-19.686,-19.556,    & 
     &-19.454,-19.375,-19.311,-19.260,-19.218,-19.184,-19.155,-19.131,    & 
     &        -21.167,-20.601,-20.206,-19.925,-19.719,-19.565,-19.447,    & 
     &-19.355,-19.283,-19.226,-19.180,-19.143,-19.112,-19.087,-19.066,    & 
     &        -20.918,-20.348,-19.976,-19.720,-19.536,-19.401,-19.299,    & 
     &-19.220,-19.159,-19.112,-19.073,-19.043,-19.018,-18.998,-18.981,    & 
     &        -20.753,-20.204,-19.847,-19.602,-19.427,-19.299,-19.203,    & 
     &-19.129,-19.072,-19.028,-18.993,-18.965,-18.942,-18.924,-18.909,    & 
     &        -20.456,-19.987,-19.677,-19.460,-19.302,-19.186,-19.098,    & 
     &-19.030,-18.978,-18.937,-18.904,-18.878,-18.857,-18.841,-18.827,    & 
     &        -20.154,-19.734,-19.461,-19.272,-19.136,-19.035,-18.960,    & 
     &-18.902,-18.858,-18.824,-18.797,-18.775,-18.759,-18.745,-18.735/ 
      DATA C6/-19.941,-19.544,-19.288,-19.114,-18.992,-18.903,-18.837,    & 
     &-18.788,-18.751,-18.723,-18.701,-18.684,-18.671,-18.661,-18.654,    & 
     &        -19.657,-19.321,-19.104,-18.956,-18.853,-18.779,-18.724,    & 
     &-18.684,-18.655,-18.632,-18.615,-18.602,-18.592,-18.585,-18.579,    & 
     &        -19.388,-19.109,-18.930,-18.810,-18.725,-18.664,-18.620,    & 
     &-18.586,-18.562,-18.543,-18.529,-18.518,-18.510,-18.503,-18.498,    & 
     &        -19.201,-18.953,-18.794,-18.686,-18.611,-18.556,-18.515,    & 
     &-18.485,-18.462,-18.446,-18.433,-18.423,-18.416,-18.410,-18.406,    & 
     &        -18.923,-18.719,-18.588,-18.500,-18.439,-18.396,-18.365,    & 
     &-18.344,-18.328,-18.318,-18.311,-18.307,-18.304,-18.303,-18.302,    & 
     &        -18.614,-18.458,-18.361,-18.298,-18.258,-18.232,-18.216,    & 
     &-18.206,-18.202,-18.201,-18.202,-18.205,-18.208,-18.213,-18.218,    & 
     &        -18.419,-18.295,-18.222,-18.178,-18.153,-18.139,-18.132,    & 
     &-18.131,-18.133,-18.138,-18.143,-18.150,-18.157,-18.164,-18.172,    & 
     &        -18.296,-18.201,-18.148,-18.118,-18.101,-18.094,-18.091,    & 
     &-18.093,-18.096,-18.101,-18.107,-18.113,-18.120,-18.126,-18.132,    & 
     &        -18.021,-17.992,-17.977,-17.970,-17.967,-17.968,-17.970,    & 
     &-17.974,-17.978,-17.983,-17.989,-17.994,-18.000,-18.005,-18.011,    & 
     &        -17.694,-17.686,-17.686,-17.691,-17.698,-17.708,-17.718,    & 
     &-17.729,-17.740,-17.750,-17.761,-17.771,-17.781,-17.790,-17.798/ 
      DATA C7/-17.374,-17.384,-17.400,-17.420,-17.440,-17.462,-17.483,    & 
     &-17.503,-17.523,-17.541,-17.558,-17.575,-17.590,-17.603,-17.616,    & 
     &        -17.169,-17.199,-17.230,-17.262,-17.293,-17.323,-17.351,    & 
     &-17.378,-17.404,-17.427,-17.449,-17.469,-17.488,-17.505,-17.520,    & 
     &        -17.151,-17.184,-17.217,-17.250,-17.282,-17.313,-17.342,    & 
     &-17.369,-17.395,-17.418,-17.440,-17.461,-17.480,-17.497,-17.513,    & 
     &        -17.230,-17.260,-17.290,-17.320,-17.348,-17.375,-17.401,    & 
     &-17.425,-17.448,-17.469,-17.489,-17.508,-17.525,-17.541,-17.556,    & 
     &        -17.379,-17.403,-17.425,-17.446,-17.467,-17.486,-17.505,    & 
     &-17.524,-17.541,-17.558,-17.574,-17.588,-17.602,-17.615,-17.627,    & 
     &        -17.596,-17.604,-17.609,-17.612,-17.616,-17.622,-17.628,    & 
     &-17.636,-17.644,-17.652,-17.661,-17.670,-17.679,-17.687,-17.695,    & 
     &        -17.846,-17.823,-17.795,-17.770,-17.750,-17.735,-17.725,    & 
     &-17.719,-17.716,-17.715,-17.716,-17.719,-17.722,-17.726,-17.730,    & 
     &        -18.089,-18.015,-17.942,-17.882,-17.836,-17.802,-17.777,    & 
     &-17.760,-17.748,-17.740,-17.736,-17.734,-17.733,-17.734,-17.736,    & 
     &        -18.299,-18.156,-18.038,-17.947,-17.881,-17.833,-17.798,    & 
     &-17.774,-17.757,-17.745,-17.738,-17.733,-17.730,-17.729,-17.729,    & 
     &        -18.441,-18.243,-18.096,-17.991,-17.915,-17.860,-17.821,    & 
     &-17.792,-17.772,-17.757,-17.746,-17.738,-17.733,-17.730,-17.728/ 
      DATA C8/-18.474,-18.262,-18.111,-18.004,-17.926,-17.869,-17.826,    & 
     &-17.795,-17.771,-17.753,-17.740,-17.730,-17.722,-17.717,-17.713,    & 
     &        -18.387,-18.191,-18.053,-17.952,-17.878,-17.823,-17.782,    & 
     &-17.752,-17.729,-17.711,-17.698,-17.689,-17.681,-17.676,-17.672,    & 
     &        -18.161,-17.990,-17.874,-17.793,-17.736,-17.696,-17.668,    & 
     &-17.648,-17.634,-17.625,-17.619,-17.616,-17.614,-17.614,-17.615,    & 
     &        -17.908,-17.774,-17.690,-17.637,-17.604,-17.583,-17.572,    & 
     &-17.567,-17.566,-17.568,-17.571,-17.576,-17.581,-17.587,-17.593,    & 
     &        -17.681,-17.589,-17.540,-17.515,-17.506,-17.505,-17.511,    & 
     &-17.520,-17.530,-17.542,-17.554,-17.566,-17.578,-17.589,-17.600,    & 
     &        -17.647,-17.606,-17.584,-17.575,-17.573,-17.576,-17.582,    & 
     &-17.589,-17.597,-17.605,-17.614,-17.623,-17.631,-17.639,-17.646,    & 
     &        -17.300,-17.291,-17.291,-17.297,-17.307,-17.319,-17.333,    & 
     &-17.347,-17.361,-17.375,-17.389,-17.402,-17.415,-17.427,-17.438,    & 
     &        -16.786,-16.802,-16.825,-16.853,-16.883,-16.914,-16.944,    & 
     &-16.974,-17.003,-17.030,-17.055,-17.079,-17.101,-17.122,-17.141,    & 
     &        -16.489,-16.533,-16.579,-16.625,-16.670,-16.713,-16.754,    & 
     &-16.793,-16.830,-16.864,-16.896,-16.925,-16.952,-16.977,-17.000,    & 
     &        -16.694,-16.724,-16.756,-16.789,-16.823,-16.856,-16.888,    & 
     &-16.919,-16.949,-16.976,-17.002,-17.026,-17.048,-17.069,-17.088/ 
      DATA C9/-16.935,-16.951,-16.971,-16.993,-17.016,-17.040,-17.064,    & 
     &-17.088,-17.111,-17.132,-17.153,-17.172,-17.190,-17.206,-17.222,    & 
     &        -17.200,-17.208,-17.220,-17.235,-17.251,-17.269,-17.286,    & 
     &-17.304,-17.322,-17.338,-17.354,-17.369,-17.384,-17.397,-17.409,    & 
     &        -17.597,-17.591,-17.589,-17.590,-17.594,-17.600,-17.608,    & 
     &-17.617,-17.626,-17.635,-17.645,-17.654,-17.662,-17.671,-17.679,    & 
     &        -18.166,-18.134,-18.107,-18.085,-18.068,-18.056,-18.047,    & 
     &-18.041,-18.038,-18.036,-18.035,-18.035,-18.036,-18.038,-18.039,    & 
     &        -19.000,-18.917,-18.838,-18.770,-18.714,-18.669,-18.632,    & 
     &-18.603,-18.579,-18.560,-18.545,-18.532,-18.522,-18.514,-18.507,    & 
     &        -20.313,-19.982,-19.754,-19.592,-19.472,-19.380,-19.309,    & 
     &-19.253,-19.208,-19.172,-19.143,-19.119,-19.099,-19.083,-19.069,    & 
     &        -19.751,-19.611,-19.520,-19.461,-19.423,-19.398,-19.382,    & 
     &-19.372,-19.366,-19.364,-19.363,-19.364,-19.366,-19.368,-19.371,    & 
     &        -19.581,-19.431,-19.337,-19.277,-19.240,-19.218,-19.207,    & 
     &-19.202,-19.203,-19.207,-19.212,-19.220,-19.228,-19.236,-19.245,    & 
     &        -19.685,-19.506,-19.389,-19.311,-19.258,-19.222,-19.199,    & 
     &-19.184,-19.175,-19.170,-19.168,-19.169,-19.171,-19.174,-19.177,    & 
     &        -19.977,-19.756,-19.606,-19.501,-19.425,-19.370,-19.330,    & 
     &-19.300,-19.278,-19.262,-19.250,-19.241,-19.235,-19.230,-19.227/ 
      DATA C10/-20.445,-20.158,-19.958,-19.815,-19.711,-19.633,-19.574,   & 
     &-19.528,-19.493,-19.465,-19.442,-19.425,-19.410,-19.398,-19.389,    & 
     &        -20.980,-20.625,-20.391,-20.229,-20.110,-20.020,-19.949,    & 
     &-19.892,-19.846,-19.807,-19.775,-19.748,-19.724,-19.704,-19.687,    & 
     &        -21.404,-21.023,-20.771,-20.594,-20.461,-20.358,-20.274,    & 
     &-20.205,-20.148,-20.099,-20.058,-20.022,-19.991,-19.965,-19.942,    & 
     &        -21.309,-20.970,-20.753,-20.603,-20.495,-20.412,-20.348,    & 
     &-20.295,-20.252,-20.215,-20.185,-20.158,-20.135,-20.115,-20.098,    & 
     &        -21.221,-20.906,-20.707,-20.574,-20.480,-20.412,-20.361,    & 
     &-20.322,-20.292,-20.268,-20.249,-20.233,-20.221,-20.210,-20.201,    & 
     &        -21.441,-21.097,-20.878,-20.728,-20.623,-20.546,-20.489,    & 
     &-20.446,-20.413,-20.387,-20.368,-20.352,-20.340,-20.330,-20.322,    & 
     &        -21.668,-21.305,-21.071,-20.911,-20.797,-20.713,-20.650,    & 
     &-20.602,-20.565,-20.536,-20.514,-20.496,-20.481,-20.470,-20.460,    & 
     &        -21.926,-21.556,-21.316,-21.150,-21.031,-20.942,-20.874,    & 
     &-20.822,-20.782,-20.750,-20.724,-20.704,-20.687,-20.674,-20.663,    & 
     &        -22.319,-21.937,-21.686,-21.510,-21.380,-21.282,-21.206,    & 
     &-21.147,-21.099,-21.061,-21.031,-21.006,-20.985,-20.968,-20.954,    & 
     &        -22.969,-22.561,-22.288,-22.092,-21.945,-21.832,-21.743,    & 
     &-21.672,-21.616,-21.570,-21.533,-21.503,-21.477,-21.457,-21.439/ 
      DATA C11/-24.001,-23.527,-23.199,-22.957,-22.772,-22.629,-22.516,   & 
     &-22.427,-22.355,-22.297,-22.250,-22.212,-22.180,-22.153,-22.131,    & 
     &        -24.233,-23.774,-23.477,-23.273,-23.128,-23.022,-22.943,    & 
     &-22.883,-22.837,-22.802,-22.774,-22.752,-22.735,-22.721,-22.710,    & 
     &        -24.550,-23.913,-23.521,-23.266,-23.094,-22.976,-22.893,    & 
     &-22.836,-22.796,-22.768,-22.750,-22.737,-22.730,-22.726,-22.725,    & 
     &        -24.301,-23.665,-23.274,-23.019,-22.848,-22.730,-22.648,    & 
     &-22.591,-22.552,-22.525,-22.507,-22.495,-22.489,-22.485,-22.485,    & 
     &        -24.519,-23.883,-23.491,-23.237,-23.065,-22.948,-22.866,    & 
     &-22.809,-22.770,-22.743,-22.724,-22.713,-22.706,-22.703,-22.702/ 
      DATA PARTCH/                                                        & 
     &     203.741,  249.643,  299.341,  353.477,  412.607,  477.237,     & 
     &     547.817,  624.786,  708.543,  799.463,  897.912, 1004.227,     & 
     &    1118.738, 1241.761, 1373.588, 1514.481, 1664.677, 1824.394,     & 
     &    1993.801, 2173.050, 2362.234, 2561.424, 2770.674, 2989.930,     & 
     &    3219.204, 3458.378, 3707.355, 3966.005, 4234.155, 4511.604,     & 
     &    4798.135, 5093.554, 5397.593, 5709.948, 6030.401, 6358.646,     & 
     &    6694.379, 7037.313, 7387.147, 7743.579, 8106.313/ 
      DATA FREQ1/0./ 
! 
      SBFCH=0. 
      IF(FR.NE.FREQ1) THEN 
         FREQ1=FR 
         WAVENO=FR/2.99792458E10 
         EVOLT=WAVENO/8065.479 
         N=int(EVOLT*10.) 
         EN=FLOAT(N)*.1 
         IF(N.LT.20) RETURN 
         IF(N.GE.105) RETURN 
! 
         DO IT=1,15 
            CROSSCHT(IT)=(CROSSCH(IT,N)+(CROSSCH(IT,N+1)-CROSSCH(IT,N))*  & 
     &                (EVOLT-EN)*10.) 
         END DO 
      END IF 
! 
!     interpolate to obtain partition function 
! 
      IF(T.GE.9000.) RETURN 
      IF(N.LT.20) RETURN 
      IF(N.GE.105) RETURN 
      IT=int((T-1000.)*twhui+1.) 
      IF(IT.LT.1) IT=1 
      TN=FLOAT(IT)*twhu+800. 
      PART=PARTCH(IT)+(PARTCH(IT+1)-PARTCH(IT))*(T-TN)*twhui 
! 
!     interpolate to obtain cross-section 
! 
      IT=int((T-2000.)*fihui+1.) 
      IF(IT.LT.1) IT=1 
      TN=FLOAT(IT)*fihu+1500. 
      SBFCH=EXP((CROSSCHT(IT)+(CROSSCHT(IT+1)-CROSSCHT(IT))*              & 
     &     (T-TN)*fihui)*tenl) 
      RETURN 
      END FUNCTION SBFCH 
! 
! 
! *********************************************************************** 
! 
! 
 
      FUNCTION SBFOH(FR,T) 
!     ==================== 
! 
!     cross-section times partition function for OH 
! 
!     from Kurucz ATLAS9 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      real(dp), parameter :: fihu=500.,fihui=1./fihu,                     & 
     &                       twhu=200.,twhui=1./twhu,                     & 
     &                       tenl=2.30258509299405 
      REAL(DP) :: CROSSOH(15,130),PARTOH(41),CROSSOHT(15) 
      REAL(DP) :: C1(150),C2(150),C3(150),C4(150),C5(150) 
      REAL(DP) :: C6(150),C7(150),C8(150),C9(150),C10(150) 
      REAL(DP) :: C11(150),C12(150),C13(150) 
      EQUIVALENCE (CROSSOH(1, 1),C1(1)),(CROSSOH(1,11),C2(1)) 
      EQUIVALENCE (CROSSOH(1,21),C3(1)),(CROSSOH(1,31),C4(1)) 
      EQUIVALENCE (CROSSOH(1,41),C5(1)),(CROSSOH(1,51),C6(1)) 
      EQUIVALENCE (CROSSOH(1,61),C7(1)),(CROSSOH(1,71),C8(1)) 
      EQUIVALENCE (CROSSOH(1,81),C9(1)),(CROSSOH(1,91),C10(1)) 
      EQUIVALENCE (CROSSOH(1,101),C11(1)) 
      EQUIVALENCE (CROSSOH(1,111),C12(1)) 
      EQUIVALENCE (CROSSOH(1,121),C13(1)) 
! 
      DATA C1/-30.855,-29.121,-27.976,-27.166,-26.566,-26.106,-25.742,    & 
     &-25.448,-25.207,-25.006,-24.836,-24.691,-24.566,-24.457,-24.363,    & 
     &        -30.494,-28.760,-27.615,-26.806,-26.206,-25.745,-25.381,    & 
     &-25.088,-24.846,-24.645,-24.475,-24.330,-24.205,-24.097,-24.002,    & 
     &        -30.157,-28.425,-27.280,-26.472,-25.872,-25.411,-25.048,    & 
     &-24.754,-24.513,-24.312,-24.142,-23.997,-23.872,-23.764,-23.669,    & 
     &        -29.848,-28.117,-26.974,-26.165,-25.566,-25.105,-24.742,    & 
     &-24.448,-24.207,-24.006,-23.836,-23.692,-23.567,-23.458,-23.364,    & 
     &        -29.567,-27.837,-26.693,-25.885,-25.286,-24.826,-24.462,    & 
     &-24.169,-23.928,-23.727,-23.557,-23.412,-23.287,-23.179,-23.084,    & 
     &        -29.307,-27.578,-26.436,-25.628,-25.029,-24.569,-24.205,    & 
     &-23.912,-23.671,-23.470,-23.300,-23.155,-23.031,-22.922,-22.828,    & 
     &        -29.068,-27.341,-26.199,-25.391,-24.792,-24.332,-23.969,    & 
     &-23.676,-23.435,-23.234,-23.064,-22.920,-22.795,-22.687,-22.592,    & 
     &        -28.820,-27.115,-25.978,-25.172,-24.574,-24.115,-23.752,    & 
     &-23.459,-23.218,-23.017,-22.848,-22.703,-22.579,-22.470,-22.376,    & 
     &        -28.540,-26.891,-25.768,-24.968,-24.372,-23.914,-23.552,    & 
     &-23.259,-23.019,-22.818,-22.649,-22.504,-22.380,-22.272,-22.177,    & 
     &        -28.275,-26.681,-25.574,-24.779,-24.186,-23.729,-23.368,    & 
     &-23.076,-22.836,-22.636,-22.467,-22.322,-22.198,-22.090,-21.996/ 
      DATA C2/-27.993,-26.470,-25.388,-24.602,-24.014,-23.560,-23.200,    & 
     &-22.909,-22.669,-22.470,-22.301,-22.157,-22.033,-21.925,-21.831,    & 
     &        -27.698,-26.252,-25.204,-24.433,-23.851,-23.401,-23.043,    & 
     &-22.754,-22.515,-22.316,-22.148,-22.005,-21.881,-21.773,-21.679,    & 
     &        -27.398,-26.026,-25.019,-24.267,-23.696,-23.251,-22.896,    & 
     &-22.609,-22.372,-22.174,-22.007,-21.864,-21.741,-21.634,-21.540,    & 
     &        -27.100,-25.791,-24.828,-24.102,-23.543,-23.106,-22.756,    & 
     &-22.472,-22.238,-22.041,-21.875,-21.733,-21.611,-21.504,-21.411,    & 
     &        -26.807,-25.549,-24.631,-23.933,-23.391,-22.964,-22.621,    & 
     &-22.341,-22.109,-21.915,-21.751,-21.610,-21.488,-21.383,-21.290,    & 
     &        -26.531,-25.310,-24.431,-23.761,-23.238,-22.823,-22.488,    & 
     &-22.214,-21.986,-21.795,-21.633,-21.494,-21.374,-21.269,-21.178,    & 
     &        -26.239,-25.066,-24.225,-23.585,-23.082,-22.681,-22.356,    & 
     &-22.089,-21.866,-21.679,-21.520,-21.383,-21.265,-21.162,-21.072,    & 
     &        -25.945,-24.824,-24.017,-23.405,-22.923,-22.538,-22.223,    & 
     &-21.964,-21.748,-21.565,-21.410,-21.276,-21.160,-21.059,-20.970,    & 
     &        -25.663,-24.587,-23.810,-23.222,-22.761,-22.391,-22.088,    & 
     &-21.838,-21.629,-21.452,-21.300,-21.170,-21.057,-20.958,-20.872,    & 
     &        -25.372,-24.350,-23.603,-23.038,-22.596,-22.241,-21.950,    & 
     &-21.710,-21.508,-21.337,-21.190,-21.064,-20.954,-20.858,-20.774/ 
      DATA C3/-25.076,-24.111,-23.396,-22.853,-22.429,-22.088,-21.809,    & 
     &-21.578,-21.384,-21.220,-21.078,-20.957,-20.851,-20.758,-20.676,    & 
     &        -24.779,-23.870,-23.189,-22.669,-22.261,-21.934,-21.667,    & 
     &-21.445,-21.259,-21.101,-20.965,-20.848,-20.746,-20.656,-20.578,    & 
     &        -24.486,-23.629,-22.983,-22.486,-22.095,-21.781,-21.524,    & 
     &-21.311,-21.132,-20.980,-20.850,-20.737,-20.639,-20.553,-20.478,    & 
     &        -24.183,-23.382,-22.774,-22.302,-21.928,-21.627,-21.381,    & 
     &-21.177,-21.005,-20.859,-20.734,-20.625,-20.531,-20.449,-20.376,    & 
     &        -23.867,-23.127,-22.561,-22.116,-21.761,-21.474,-21.238,    & 
     &-21.043,-20.878,-20.738,-20.617,-20.513,-20.423,-20.344,-20.274,    & 
     &        -23.538,-22.862,-22.340,-21.926,-21.592,-21.320,-21.096,    & 
     &-20.909,-20.751,-20.617,-20.502,-20.402,-20.315,-20.239,-20.172,    & 
     &        -23.234,-22.604,-22.120,-21.734,-21.422,-21.166,-20.953,    & 
     &-20.776,-20.625,-20.497,-20.387,-20.291,-20.208,-20.135,-20.071,    & 
     &        -22.934,-22.347,-21.898,-21.541,-21.250,-21.010,-20.811,    & 
     &-20.643,-20.500,-20.378,-20.273,-20.182,-20.102,-20.033,-19.971,    & 
     &        -22.637,-22.092,-21.676,-21.345,-21.075,-20.853,-20.666,    & 
     &-20.508,-20.374,-20.259,-20.159,-20.073,-19.997,-19.931,-19.872,    & 
     &        -22.337,-21.835,-21.452,-21.147,-20.899,-20.693,-20.520,    & 
     &-20.373,-20.247,-20.139,-20.046,-19.964,-19.892,-19.830,-19.774/ 
      DATA C4/-22.049,-21.584,-21.230,-20.950,-20.721,-20.531,-20.372,    & 
     &-20.236,-20.119,-20.019,-19.931,-19.855,-19.788,-19.729,-19.676,    & 
     &        -21.768,-21.337,-21.011,-20.754,-20.544,-20.370,-20.223,    & 
     &-20.098,-19.991,-19.898,-19.817,-19.746,-19.683,-19.628,-19.579,    & 
     &        -21.494,-21.096,-20.796,-20.559,-20.367,-20.208,-20.074,    & 
     &-19.960,-19.861,-19.776,-19.701,-19.636,-19.578,-19.527,-19.482,    & 
     &        -21.233,-20.861,-20.585,-20.368,-20.193,-20.048,-19.926,    & 
     &-19.821,-19.732,-19.654,-19.586,-19.526,-19.473,-19.426,-19.384,    & 
     &        -20.983,-20.635,-20.380,-20.181,-20.021,-19.889,-19.778,    & 
     &-19.683,-19.602,-19.531,-19.469,-19.415,-19.367,-19.324,-19.286,    & 
     &        -20.743,-20.418,-20.182,-19.999,-19.853,-19.733,-19.633,    & 
     &-19.547,-19.474,-19.410,-19.354,-19.305,-19.261,-19.223,-19.189,    & 
     &        -20.515,-20.210,-19.991,-19.824,-19.690,-19.581,-19.490,    & 
     &-19.413,-19.347,-19.290,-19.240,-19.196,-19.157,-19.122,-19.092,    & 
     &        -20.297,-20.011,-19.808,-19.654,-19.532,-19.434,-19.352,    & 
     &-19.282,-19.223,-19.172,-19.127,-19.088,-19.054,-19.023,-18.996,    & 
     &        -20.090,-19.822,-19.633,-19.491,-19.381,-19.291,-19.218,    & 
     &-19.156,-19.103,-19.057,-19.018,-18.983,-18.952,-18.925,-18.901,    & 
     &        -19.893,-19.642,-19.467,-19.337,-19.236,-19.155,-19.089,    & 
     &-19.034,-18.987,-18.946,-18.912,-18.881,-18.854,-18.831,-18.810/ 
      DATA C5/-19.705,-19.472,-19.309,-19.190,-19.098,-19.025,-18.966,    & 
     &-18.917,-18.876,-18.840,-18.810,-18.783,-18.760,-18.739,-18.721,    & 
     &        -19.527,-19.310,-19.161,-19.051,-18.968,-18.903,-18.851,    & 
     &-18.807,-18.771,-18.740,-18.713,-18.690,-18.670,-18.653,-18.637,    & 
     &        -19.357,-19.159,-19.022,-18.922,-18.847,-18.789,-18.743,    & 
     &-18.704,-18.673,-18.646,-18.623,-18.603,-18.586,-18.571,-18.558,    & 
     &        -19.195,-19.016,-18.892,-18.803,-18.736,-18.684,-18.643,    & 
     &-18.610,-18.583,-18.560,-18.540,-18.523,-18.509,-18.496,-18.485,    & 
     &        -19.042,-18.883,-18.772,-18.693,-18.634,-18.589,-18.553,    & 
     &-18.525,-18.501,-18.481,-18.465,-18.451,-18.438,-18.428,-18.419,    & 
     &        -18.894,-18.758,-18.662,-18.593,-18.542,-18.503,-18.473,    & 
     &-18.448,-18.428,-18.412,-18.398,-18.386,-18.376,-18.367,-18.359,    & 
     &        -18.752,-18.639,-18.559,-18.501,-18.458,-18.426,-18.400,    & 
     &-18.380,-18.363,-18.350,-18.338,-18.328,-18.320,-18.313,-18.306,    & 
     &        -18.611,-18.523,-18.460,-18.415,-18.381,-18.355,-18.334,    & 
     &-18.318,-18.304,-18.293,-18.284,-18.276,-18.269,-18.263,-18.258,    & 
     &        -18.471,-18.408,-18.362,-18.329,-18.304,-18.285,-18.269,    & 
     &-18.257,-18.247,-18.238,-18.231,-18.224,-18.219,-18.214,-18.210,    & 
     &        -18.330,-18.290,-18.261,-18.239,-18.223,-18.211,-18.201,    & 
     &-18.192,-18.185,-18.179,-18.174,-18.169,-18.165,-18.162,-18.159/ 
      DATA C6/-18.190,-18.168,-18.154,-18.143,-18.135,-18.129,-18.124,    & 
     &-18.120,-18.116,-18.112,-18.109,-18.106,-18.104,-18.102,-18.100,    & 
     &        -18.055,-18.047,-18.043,-18.042,-18.040,-18.039,-18.039,    & 
     &-18.038,-18.037,-18.036,-18.035,-18.034,-18.033,-18.033,-18.032,    & 
     &        -17.929,-17.931,-17.935,-17.939,-17.943,-17.946,-17.948,    & 
     &-17.950,-17.952,-17.953,-17.955,-17.956,-17.957,-17.958,-17.959,    & 
     &        -17.818,-17.826,-17.834,-17.842,-17.849,-17.855,-17.860,    & 
     &-17.865,-17.869,-17.872,-17.875,-17.878,-17.881,-17.883,-17.886,    & 
     &        -17.724,-17.736,-17.747,-17.758,-17.767,-17.775,-17.782,    & 
     &-17.788,-17.793,-17.798,-17.803,-17.807,-17.811,-17.815,-17.819,    & 
     &        -17.651,-17.665,-17.678,-17.690,-17.701,-17.710,-17.718,    & 
     &-17.725,-17.732,-17.738,-17.744,-17.749,-17.755,-17.760,-17.765,    & 
     &        -17.601,-17.615,-17.629,-17.642,-17.653,-17.663,-17.672,    & 
     &-17.680,-17.688,-17.695,-17.701,-17.708,-17.714,-17.720,-17.726,    & 
     &        -17.572,-17.587,-17.602,-17.614,-17.626,-17.636,-17.645,    & 
     &-17.654,-17.662,-17.670,-17.677,-17.684,-17.691,-17.698,-17.704,    & 
     &        -17.565,-17.581,-17.595,-17.607,-17.619,-17.629,-17.638,    & 
     &-17.647,-17.656,-17.664,-17.671,-17.679,-17.686,-17.693,-17.700,    & 
     &        -17.580,-17.594,-17.608,-17.620,-17.630,-17.640,-17.650,    & 
     &-17.658,-17.667,-17.675,-17.682,-17.690,-17.697,-17.704,-17.711/ 
      DATA C7/-17.613,-17.626,-17.639,-17.649,-17.659,-17.669,-17.677,    & 
     &-17.686,-17.694,-17.701,-17.709,-17.716,-17.723,-17.730,-17.737,    & 
     &        -17.663,-17.675,-17.685,-17.695,-17.703,-17.711,-17.719,    & 
     &-17.727,-17.734,-17.741,-17.748,-17.755,-17.761,-17.768,-17.774,    & 
     &        -17.728,-17.737,-17.745,-17.752,-17.759,-17.766,-17.772,    & 
     &-17.778,-17.785,-17.791,-17.797,-17.803,-17.808,-17.814,-17.820,    & 
     &        -17.803,-17.809,-17.814,-17.818,-17.823,-17.828,-17.832,    & 
     &-17.837,-17.842,-17.847,-17.852,-17.856,-17.861,-17.866,-17.871,    & 
     &        -17.884,-17.886,-17.888,-17.889,-17.891,-17.893,-17.896,    & 
     &-17.899,-17.902,-17.905,-17.908,-17.912,-17.915,-17.919,-17.922,    & 
     &        -17.966,-17.964,-17.961,-17.959,-17.958,-17.958,-17.958,    & 
     &-17.959,-17.960,-17.961,-17.963,-17.964,-17.966,-17.968,-17.970,    & 
     &        -18.040,-18.034,-18.028,-18.023,-18.019,-18.016,-18.013,    & 
     &-18.012,-18.010,-18.010,-18.009,-18.009,-18.009,-18.009,-18.010,    & 
     &        -18.096,-18.087,-18.078,-18.071,-18.065,-18.059,-18.055,    & 
     &-18.051,-18.047,-18.045,-18.042,-18.040,-18.039,-18.037,-18.036,    & 
     &        -18.125,-18.115,-18.105,-18.097,-18.089,-18.082,-18.076,    & 
     &-18.070,-18.065,-18.061,-18.057,-18.053,-18.051,-18.048,-18.046,    & 
     &        -18.120,-18.112,-18.103,-18.095,-18.087,-18.079,-18.072,    & 
     &-18.066,-18.060,-18.055,-18.050,-18.046,-18.042,-18.039,-18.036/ 
      DATA C8/-18.083,-18.078,-18.071,-18.064,-18.057,-18.050,-18.044,    & 
     &-18.037,-18.032,-18.026,-18.022,-18.017,-18.014,-18.010,-18.007,    & 
     &        -18.025,-18.022,-18.017,-18.012,-18.006,-18.000,-17.994,    & 
     &-17.989,-17.984,-17.979,-17.975,-17.971,-17.968,-17.965,-17.963,    & 
     &        -17.957,-17.955,-17.952,-17.948,-17.943,-17.938,-17.934,    & 
     &-17.929,-17.925,-17.922,-17.918,-17.916,-17.913,-17.911,-17.910,    & 
     &        -17.890,-17.889,-17.886,-17.882,-17.879,-17.875,-17.871,    & 
     &-17.867,-17.864,-17.862,-17.860,-17.858,-17.857,-17.856,-17.855,    & 
     &        -17.831,-17.829,-17.826,-17.822,-17.819,-17.815,-17.812,    & 
     &-17.810,-17.807,-17.806,-17.804,-17.803,-17.803,-17.803,-17.803,    & 
     &        -17.786,-17.782,-17.777,-17.773,-17.769,-17.766,-17.763,    & 
     &-17.761,-17.759,-17.758,-17.757,-17.757,-17.757,-17.758,-17.759,    & 
     &        -17.753,-17.747,-17.741,-17.735,-17.731,-17.727,-17.724,    & 
     &-17.722,-17.721,-17.720,-17.720,-17.720,-17.721,-17.722,-17.724,    & 
     &        -17.733,-17.724,-17.716,-17.709,-17.703,-17.699,-17.696,    & 
     &-17.694,-17.693,-17.692,-17.692,-17.693,-17.694,-17.695,-17.697,    & 
     &        -17.723,-17.711,-17.700,-17.691,-17.685,-17.680,-17.676,    & 
     &-17.674,-17.673,-17.672,-17.673,-17.673,-17.675,-17.676,-17.678,    & 
     &        -17.718,-17.702,-17.689,-17.679,-17.672,-17.667,-17.663,    & 
     &-17.660,-17.659,-17.659,-17.659,-17.660,-17.661,-17.663,-17.665/ 
      DATA C9/-17.713,-17.695,-17.681,-17.670,-17.662,-17.656,-17.653,    & 
     &-17.650,-17.649,-17.649,-17.649,-17.650,-17.651,-17.653,-17.655,    & 
     &        -17.705,-17.686,-17.671,-17.660,-17.652,-17.647,-17.643,    & 
     &-17.641,-17.640,-17.640,-17.640,-17.641,-17.643,-17.645,-17.647,    & 
     &        -17.690,-17.671,-17.657,-17.647,-17.640,-17.635,-17.632,    & 
     &-17.630,-17.630,-17.630,-17.631,-17.632,-17.634,-17.636,-17.639,    & 
     &        -17.667,-17.649,-17.637,-17.629,-17.623,-17.619,-17.618,    & 
     &-17.617,-17.617,-17.618,-17.619,-17.621,-17.623,-17.626,-17.628,    & 
     &        -17.635,-17.621,-17.611,-17.605,-17.601,-17.600,-17.599,    & 
     &-17.599,-17.601,-17.602,-17.604,-17.607,-17.609,-17.612,-17.615,    & 
     &        -17.596,-17.585,-17.579,-17.576,-17.575,-17.575,-17.576,    & 
     &-17.578,-17.580,-17.582,-17.585,-17.588,-17.591,-17.595,-17.598,    & 
     &        -17.550,-17.544,-17.542,-17.542,-17.544,-17.546,-17.548,    & 
     &-17.552,-17.555,-17.558,-17.562,-17.566,-17.570,-17.573,-17.577,    & 
     &        -17.501,-17.500,-17.501,-17.504,-17.508,-17.513,-17.517,    & 
     &-17.521,-17.526,-17.530,-17.535,-17.539,-17.544,-17.548,-17.553,    & 
     &        -17.449,-17.452,-17.457,-17.463,-17.470,-17.476,-17.482,    & 
     &-17.488,-17.493,-17.499,-17.504,-17.509,-17.514,-17.519,-17.524,    & 
     &        -17.396,-17.403,-17.412,-17.420,-17.429,-17.437,-17.444,    & 
     &-17.451,-17.458,-17.464,-17.470,-17.476,-17.481,-17.487,-17.492/ 
      DATA C10/-17.344,-17.355,-17.366,-17.377,-17.387,-17.396,-17.405,   & 
     &-17.413,-17.420,-17.427,-17.434,-17.440,-17.446,-17.452,-17.458,    & 
     &        -17.295,-17.307,-17.321,-17.333,-17.345,-17.355,-17.365,    & 
     &-17.373,-17.382,-17.389,-17.397,-17.404,-17.410,-17.417,-17.423,    & 
     &        -17.249,-17.264,-17.278,-17.292,-17.304,-17.316,-17.326,    & 
     &-17.335,-17.344,-17.352,-17.360,-17.368,-17.375,-17.382,-17.389,    & 
     &        -17.209,-17.225,-17.241,-17.255,-17.268,-17.280,-17.291,    & 
     &-17.301,-17.310,-17.319,-17.327,-17.335,-17.343,-17.350,-17.357,    & 
     &        -17.177,-17.194,-17.210,-17.225,-17.239,-17.251,-17.262,    & 
     &-17.272,-17.282,-17.291,-17.300,-17.308,-17.316,-17.324,-17.331,    & 
     &        -17.154,-17.172,-17.189,-17.204,-17.218,-17.230,-17.242,    & 
     &-17.252,-17.262,-17.272,-17.280,-17.289,-17.298,-17.306,-17.314,    & 
     &        -17.144,-17.162,-17.179,-17.194,-17.208,-17.220,-17.232,    & 
     &-17.242,-17.253,-17.262,-17.271,-17.280,-17.289,-17.297,-17.306,    & 
     &        -17.146,-17.164,-17.181,-17.196,-17.210,-17.222,-17.234,    & 
     &-17.245,-17.255,-17.265,-17.274,-17.283,-17.292,-17.301,-17.309,    & 
     &        -17.163,-17.180,-17.197,-17.212,-17.225,-17.237,-17.249,    & 
     &-17.260,-17.270,-17.280,-17.289,-17.298,-17.307,-17.316,-17.325,    & 
     &        -17.193,-17.211,-17.227,-17.241,-17.254,-17.266,-17.277,    & 
     &-17.288,-17.298,-17.308,-17.317,-17.327,-17.336,-17.345,-17.353/ 
      DATA C11/-17.239,-17.256,-17.271,-17.284,-17.297,-17.309,-17.320,   & 
     &-17.330,-17.340,-17.350,-17.359,-17.369,-17.378,-17.387,-17.395,    & 
     &        -17.299,-17.315,-17.329,-17.342,-17.354,-17.365,-17.376,    & 
     &-17.386,-17.396,-17.405,-17.415,-17.424,-17.433,-17.442,-17.451,    & 
     &        -17.373,-17.388,-17.402,-17.414,-17.425,-17.436,-17.446,    & 
     &-17.456,-17.466,-17.475,-17.484,-17.493,-17.502,-17.511,-17.520,    & 
     &        -17.462,-17.476,-17.489,-17.500,-17.511,-17.521,-17.531,    & 
     &-17.541,-17.550,-17.559,-17.569,-17.578,-17.587,-17.595,-17.604,    & 
     &        -17.567,-17.581,-17.592,-17.603,-17.613,-17.623,-17.632,    & 
     &-17.641,-17.651,-17.660,-17.669,-17.678,-17.686,-17.695,-17.704,    & 
     &        -17.689,-17.701,-17.712,-17.722,-17.732,-17.741,-17.750,    & 
     &-17.759,-17.768,-17.777,-17.786,-17.795,-17.803,-17.812,-17.821,    & 
     &        -17.829,-17.840,-17.851,-17.860,-17.869,-17.878,-17.887,    & 
     &-17.896,-17.904,-17.913,-17.922,-17.930,-17.939,-17.948,-17.956,    & 
     &        -17.988,-18.000,-18.010,-18.019,-18.028,-18.036,-18.045,    & 
     &-18.053,-18.062,-18.070,-18.079,-18.087,-18.096,-18.104,-18.112,    & 
     &        -18.171,-18.183,-18.192,-18.201,-18.210,-18.218,-18.227,    & 
     &-18.235,-18.243,-18.252,-18.260,-18.268,-18.277,-18.285,-18.293,    & 
     &        -18.381,-18.393,-18.403,-18.413,-18.422,-18.430,-18.438,    & 
     &-18.447,-18.455,-18.463,-18.471,-18.479,-18.487,-18.495,-18.503/ 
      DATA C12/-18.625,-18.638,-18.650,-18.660,-18.669,-18.678,-18.687,   & 
     &-18.695,-18.703,-18.711,-18.719,-18.726,-18.734,-18.742,-18.750,    & 
     &        -18.912,-18.929,-18.943,-18.955,-18.966,-18.975,-18.984,    & 
     &-18.993,-19.001,-19.008,-19.016,-19.023,-19.031,-19.038,-19.045,    & 
     &        -19.260,-19.283,-19.303,-19.320,-19.333,-19.345,-19.355,    & 
     &-19.364,-19.372,-19.380,-19.387,-19.394,-19.400,-19.407,-19.413,    & 
     &        -19.704,-19.740,-19.771,-19.796,-19.816,-19.832,-19.845,    & 
     &-19.855,-19.863,-19.870,-19.876,-19.882,-19.887,-19.892,-19.897,    & 
     &        -20.339,-20.386,-20.424,-20.454,-20.476,-20.492,-20.502,    & 
     &-20.509,-20.513,-20.516,-20.518,-20.520,-20.521,-20.523,-20.524,    & 
     &        -21.052,-21.075,-21.093,-21.105,-21.114,-21.120,-21.123,    & 
     &-21.125,-21.126,-21.127,-21.128,-21.130,-21.131,-21.133,-21.135,    & 
     &        -21.174,-21.203,-21.230,-21.255,-21.278,-21.299,-21.320,    & 
     &-21.339,-21.357,-21.375,-21.392,-21.408,-21.424,-21.439,-21.454,    & 
     &        -21.285,-21.317,-21.346,-21.372,-21.395,-21.416,-21.435,    & 
     &-21.452,-21.468,-21.483,-21.497,-21.511,-21.524,-21.536,-21.548,    & 
     &        -21.396,-21.429,-21.459,-21.486,-21.511,-21.532,-21.551,    & 
     &-21.569,-21.585,-21.600,-21.614,-21.627,-21.640,-21.652,-21.663,    & 
     &        -21.516,-21.549,-21.580,-21.609,-21.635,-21.658,-21.678,    & 
     &-21.696,-21.713,-21.728,-21.742,-21.755,-21.767,-21.779,-21.790/ 
      DATA C13/-21.651,-21.681,-21.711,-21.738,-21.763,-21.785,-21.804,   & 
     &-21.821,-21.837,-21.851,-21.864,-21.876,-21.887,-21.898,-21.908,    & 
     &        -21.810,-21.831,-21.853,-21.874,-21.893,-21.910,-21.925,    & 
     &-21.938,-21.950,-21.961,-21.971,-21.980,-21.989,-21.998,-22.006,    & 
     &        -22.009,-22.016,-22.026,-22.037,-22.048,-22.058,-22.066,    & 
     &-22.074,-22.081,-22.088,-22.094,-22.099,-22.105,-22.111,-22.117,    & 
     &        -22.353,-22.317,-22.296,-22.284,-22.276,-22.270,-22.266,    & 
     &-22.262,-22.260,-22.258,-22.257,-22.257,-22.257,-22.258,-22.259,    & 
     &        -22.705,-22.609,-22.552,-22.515,-22.488,-22.468,-22.451,    & 
     &-22.438,-22.427,-22.418,-22.410,-22.405,-22.400,-22.397,-22.395,    & 
     &        -22.889,-22.791,-22.731,-22.690,-22.659,-22.634,-22.612,    & 
     &-22.594,-22.579,-22.566,-22.555,-22.546,-22.539,-22.533,-22.528,    & 
     &        -23.211,-23.109,-23.041,-22.989,-22.945,-22.906,-22.872,    & 
     &-22.842,-22.816,-22.793,-22.774,-22.757,-22.743,-22.732,-22.722,    & 
     &        -25.312,-24.669,-24.250,-23.959,-23.746,-23.587,-23.463,    & 
     &-23.366,-23.288,-23.225,-23.173,-23.131,-23.095,-23.066,-23.041,    & 
     &        -25.394,-24.752,-24.333,-24.041,-23.829,-23.669,-23.546,    & 
     &-23.449,-23.371,-23.308,-23.256,-23.214,-23.178,-23.149,-23.124,    & 
     &        -25.430,-24.787,-24.369,-24.077,-23.865,-23.705,-23.582,    & 
     &-23.484,-23.407,-23.344,-23.292,-23.249,-23.214,-23.185,-23.160/ 
      DATA PARTOH/                                                        & 
     &   145.979,  178.033,  211.618,  247.053,  284.584,  324.398,       & 
     &   366.639,  411.425,  458.854,  509.012,  561.976,  617.823,       & 
     &   676.626,  738.448,  803.363,  871.437,  942.735, 1017.330,       & 
     &  1095.284, 1176.654, 1261.510, 1349.898, 1441.875, 1537.483,       & 
     &  1636.753, 1739.733, 1846.434, 1956.883, 2071.080, 2189.029,       & 
     &  2310.724, 2436.155, 2565.283, 2698.103, 2834.571, 2974.627,       & 
     &  3118.242, 3265.366, 3415.912, 3569.837, 3727.077/ 
      DATA FREQ1/0./ 
! 
      SBFOH=0. 
      IF(FR.NE.FREQ1) THEN 
         FREQ1=FR 
         WAVENO=FR/2.99792458E10 
         EVOLT=WAVENO/8065.479 
         N=int(EVOLT*10.-20.) 
         EN=FLOAT(N)*.1+2. 
         IF(N.LE.0) RETURN 
         IF(N.GE.130) RETURN 
         DO IT=1,15 
            CROSSOHT(IT)=(CROSSOH(IT,N)+(CROSSOH(IT,N+1)-CROSSOH(IT,N))*  & 
     &                   (EVOLT-EN)*10.) 
         END DO 
      END IF 
! 
!     interpolate to obtain partition function 
! 
      IF(T.GE.9000.) RETURN 
      IF(N.LE.0) RETURN 
      IF(N.GE.130) RETURN 
      IT=int((T-1000.)*twhui+1.) 
      IF(IT.LT.1) IT=1 
      TN=FLOAT(IT)*twhu+800. 
      PART=PARTOH(IT)+(PARTOH(IT+1)-PARTOH(IT))*(T-TN)*twhui 
! 
!     interpolate to obtain cross-section 
! 
      IT=int((T-2000.)*fihui+1.) 
      IF(IT.LT.1) IT=1 
      TN=FLOAT(IT)*fihu+1500. 
      SBFOH=EXP((CROSSOHT(IT)+(CROSSOHT(IT+1)-CROSSOHT(IT))*              & 
     &      (T-TN)*fihui)*tenl) 
      RETURN 
      END FUNCTION SBFOH 
! 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE XENINI 
!     ================= 
! 
!     Initializes necessary arrays for evaluating hydrogen line profiles 
!     from the XENOMORPH tables 
! 
      use accura 
      use params 
      use modelp 
      use hydxen

      implicit real(dp) (a-h,o-z),logical (l) 
 
! 
      DO I=1,4 
         DO J=1,22 
            ILXEN(I,J)=0 
         END DO 
      END DO 
      if(ihxenb.gt.0) then 
         call alloc_hydxen
         ihxenb=23 
         ihxenr=ihxenb+1 
         open(unit=ihxenb,file='xenomorph.blue.dat',status='old') 
         open(unit=ihxenr,file='xenomorph.red.dat',status='old') 
         write(6,"(' -----------'/                                        & 
     &       ' reading XENOMORPH tables; ihxen =',2i3,/                   & 
     &       ' -----------')") ihxenb,ihxenr 
       else 
         return 
      end if 
! 
! --------------------------------- 
!     read  tables - blue wing 
! --------------------------------- 
! 
      ILINE=0 
      READ(IHXENB,*) NTAB 
      DO ITAB=1,NTAB 
         ILINEB=ILINE 
         READ(IHXENB,*) NLXEN 
         DO ILI=1,NLXEN 
            ILINE=ILINE+1 
            READ(IHXENB,*) I,J,ALMIN,ANEMIN,TMIN,DLA,DLE,DLT,             & 
     &                  NWL,NE,NT 
            XNEMIN=ANEMIN 
            ILXEN(I,J)=ILINE 
            NWLXEN(ILINE)=NWL 
            NTHXEN(ILINE)=NT 
            NEHXEN(ILINE)=NE 
            DO IWL=1,NWL 
               ALXEN(ILINE,IWL)=ALMIN+(IWL-1)*DLA 
            END DO 
            DO INE=1,NE 
               XNEXEN(INE,ILINE)=ANEMIN+(INE-1)*DLE 
            END DO 
            DO IT=1,NT 
               XTXEN(IT,ILINE)=TMIN+(IT-1)*DLT 
            END DO 
         END DO 
! 
         DO ILI=1,NLXEN 
            ILNE=ILINEB+ILI 
            NWL=NWLXEN(ILNE) 
            READ(IHXENB,"(1X)") 
            DO INE=1,NEHXEN(ILNE) 
               DO IT=1,NTHXEN(ILNE) 
                  READ(IHXENB,*) QLT,(PRFXB(ILNE,IWL,IT,INE),IWL=1,NWL) 
               END DO 
            END DO 
! 
         END DO 
      END DO 
      CLOSE(IHXENB) 
! 
! --------------------------------- 
!     read  tables - red wing 
! --------------------------------- 
! 
      ILINE=0 
      READ(IHXENR,*) NTAB 
      DO ITAB=1,NTAB 
         ILINEB=ILINE 
         READ(IHXENR,*) NLXEN 
         DO ILI=1,NLXEN 
            ILINE=ILINE+1 
            READ(IHXENR,*) I,J,ALMIN,ANEMIN,TMIN,DLA,DLE,DLT,             & 
     &                  NWL,NE,NT 
         END DO 
! 
         DO ILI=1,NLXEN 
            ILNE=ILINEB+ILI 
            NWL=NWLXEN(ILNE) 
            READ(IHXENR,"(1X)") 
            DO INE=1,NEHXEN(ILNE) 
               DO IT=1,NTHXEN(ILNE) 
                  READ(IHXENR,*) QLT,(PRFXR(ILNE,IWL,IT,INE),IWL=1,NWL) 
               END DO 
            END DO 
! 
         END DO 
      END DO 
! 
!     interpolation to the actual values of temperature and electron 
!     density 
! 
      do id =1,nd 
         tl=log10(temp(id)) 
         anel=log10(elec(id)) 
         do ili=1,nlxen 
            iline=ilineb+ili 
            nwl=nwlxen(iline) 
            do iwl=1,nwl 
               call intxen(prfb0,prfr0,tl,anel,iwl,iline) 
               prfb(iline,id,iwl)=prfb0 
               prfr(iline,id,iwl)=prfb0 
            end do 
         end do 
      end do 
      CLOSE(IHXENR) 
! 
      RETURN 
      END SUBROUTINE XENINI 
! 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE INTXEN(W0B,W0R,X0,Z0,IWL,ILINE) 
!     ========================================== 
! 
!     Interpolation in temperature and electron density from the 
!     Xenomorph tables for hydrogen lines to the actual valus of 
!     temperature and electron density 
! 
      use accura 
      use params 
      use modelp 
      use hydxen

      implicit real(dp) (a-h,o-z),logical (l) 
 
      REAL(DP) :: ZZ(3),XX(3),WXB(3),WZB(3),WXR(3),WZR(3) 
! 
      NX=2 
      NZ=2 
      NT=NTHXEN(ILINE) 
      NE=NEHXEN(ILINE) 
! 
      DO IZZ=1,NE-1 
         IPZ=IZZ 
         IF(Z0.LE.XNEXEN(IZZ+1,ILINE)) EXIT 
      END DO 
      N0Z=IPZ-NZ/2+1 
      IF(N0Z.LT.1) N0Z=1 
      IF(N0Z.GT.NE-NZ+1) N0Z=NE-NZ+1 
      N1Z=N0Z+NZ-1 
! 
      DO IZZ=N0Z,N1Z 
         I0Z=IZZ-N0Z+1 
         ZZ(I0Z)=XNEXEN(IZZ,ILINE) 
         DO IX=1,NT-1 
            IPX=IX 
            IF(X0.LE.XTXEN(IX+1,ILINE)) EXIT 
         END DO 
         N0X=IPX-NX/2+1 
         IF(N0X.LT.1) N0X=1 
         IF(N0X.GT.NT-NX+1) N0X=NT-NX+1 
         N1X=N0X+NX-1 
         DO IX=N0X,N1X 
            I0=IX-N0X+1 
            XX(I0)=XTXEN(IX,ILINE) 
            WXB(I0)=PRFXB(ILINE,IWL,IX,IZZ) 
            WXR(I0)=PRFXR(ILINE,IWL,IX,IZZ) 
         END DO 
         WZB(I0Z)=YINT(XX,WXB,X0) 
         WZR(I0Z)=YINT(XX,WXR,X0) 
      END DO 
      W0B=YINT(ZZ,WZB,Z0) 
      W0R=YINT(ZZ,WZR,Z0) 
      RETURN 
      END SUBROUTINE INTXEN 
! 
! 
!     ****************************************************************** 
! 
! 
      SUBROUTINE GOMINI 
!     ================= 
! 
      use accura 
      use params 
      use modelp 
      use hydprf 
      implicit real(dp) (a-h,o-z),logical (l) 
      real(dp), allocatable :: hydcrs(:,:,:) 
 
      real(dp) :: temvec(mtabth),elevec(mtabeh) 
 
      allocate (hydcrs(mtabth,mtabeh,mfhtab)) 
! 
      if(ihgom.eq.0) return 
! 
      open(53,file='gomhyd.dat',status='old') 
! 
      read(53,*) nugfreq,nugtemp,nugele 
      read(53,*) 
      read(53,*) (temvec(i),i=1,nugtemp) 
      read(53,*) 
      read(53,*) (elevec(j),j=1,nugele) 
      do it=1,nugtemp 
         temvec(it)=log(temvec(it)*1.161e4) 
      end do 
! 
      EGTAB1 = elevec(1) 
      EGTAB2 = elevec(nugele) 
      TGTAB1 = temvec(1) 
      TGTAB2 = temvec(nugtemp) 
! 
      do k = 1, nugfreq 
         read(53,"(40x,f17.14)") eneev 
         frgtab(k)=3.28805e15/13.595*eneev 
         wlgtab(k)=2.997925e18/frgtab(k) 
         do i = 1, nugtemp 
            read(53,*) (hydcrs(i,j,k),j=1,nugele) 
         end do 
      end do 
      frg1=frgtab(1) 
      frg2=frgtab(nugfreq) 
      close(53) 
! 
!     Interpolate to the actual temperature and electron density 
!     at the individual depth points 
! 
      depths: do id=1,nd 
         if(elec(id).lt.HGLIM) cycle depths 
         rl=log(elec(id)) 
         tl=log(temp(id)) 
! 
         DELTAR=(RL-EGTAB1)/(EGTAB2-EGTAB1)*FLOAT(nugele-1) 
         JR = 1 + IDINT(DELTAR) 
         IF(JR.LT.1) JR = 1 
         IF(JR.GT.(nugele-1)) JR = nugele-1 
         r1i=elevec(jr) 
         r2i=elevec(jr+1) 
         dri=(RL-R1i)/(R2i-R1i) 
         if(JR .eq. 1) dri = 0. 
! 
         DELTAT=(TL-TGTAB1)/(TGTAB2-TGTAB1)*FLOAT(nugtemp-1) 
         JP = 1 + IDINT(DELTAT) 
         IF(JP.LT.1) JP = 1 
         IF(JP.GT.nugtemp-1) JP = nugtemp-1 
         t1i=temvec(jp) 
         t2i=temvec(jp+1) 
         dti=(TL-T1i)/(T2i-T1i) 
         if(JP .eq. 1) dti = 0. 
! 
!        loop over tabular frequencies 
! 
         do jf=1,nugfreq 
            opr1=hydcrs(jp,jr,jf)+dti*                                    & 
     &           (hydcrs(jp+1,jr,jf)-hydcrs(jp,jr,jf)) 
            opr2=hydcrs(jp,jr+1,jf)+dti*                                  & 
     &           (hydcrs(jp+1,jr+1,jf)-hydcrs(jp,jr+1,jf)) 
            opac=opr1+dri*(opr2-opr1) 
            hydopg(jf,id)=opac+log(0.02654*4.1347e-15) 
         end do 
      end do depths 
      deallocate (hydcrs) 
      return 
      end subroutine gomini 
! 
!     ****************************************************\ 
! 
! 
      subroutine ghydop(id,i0,i1,pj,absoh,emish) 
!     ========================================== 
! 
!     hydrogen opacity -- lines + pseudocontinuum from Gomez tables 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use hydprf 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      real(dp) :: absoh(mfreq),emish(mfreq),pj(40) 
! 
      frg1=frgtab(1) 
      frg2=frgtab(nugfreq) 
      frq: do ij=i0,i1 
         fr=freq(ij) 
         if(fr.lt.frg1.or.fr.gt.frg2) cycle frq 
         wla=2.997925e18/fr 
         frl=log10(fr) 
! 
         if(ij.eq.i0) igf=nugfreq 
         ig: do 
            if(wla.gt.wlgtab(igf)) then 
               igf=igf-1 
               cycle ig 
             else 
               exit ig 
            end if 
         end do ig 
         ig0=igf 
         if(ig0.le.2) ig0=2 
         ig1=igf-1 
         abl=(hydopg(ig1,id)-hydopg(ig0,id))*(wla-wlgtab(ig0))/           & 
     &       (wlgtab(ig1)-wlgtab(ig0))+hydopg(ig0,id) 
! 
         ii=1 
         if(freq(ij).gt.8.22013e14) then 
            pp=pj(1)*2. 
          else 
            pp=pj(2)*8. 
         end if 
! 
         F15=FR*1.E-15 
         XKF=EXP(-4.79928e-11*FR/TEMP(ID)) 
         XKFB=XKF*1.4743E-2*F15*F15*F15 
 
         oph=exp(abl)*pp 
         absoh(ij)=absoh(ij)+oph 
         emish(ij)=emish(ij)+oph*xkfb/(1.-xkf) 
      end do frq 
! 
      return 
      end subroutine ghydop 
 
! 
! ******************************************************************** 
! 
      subroutine ingrid(mode,inext,igrd) 
!     ================================== 
! 
!     setting state parameters for the opacity grid calculations 
! 
!     input: 
!           temp1 - lowest value of T 
!           temp2 - largest value of T 
!           ntemp - number of temperature values 
!           dens1 - lowest value of the density parameter 
!           dens2 - largest value of the density parameter 
!           ndens - number of the density parameter values 
! 
!           isdens = 0 - density parameter is electron density 
!                  > 0 - density parameter is mass density 
!                  < 0 - density parameter is gas pressure 
! 
! 
      use accura 
      use params 
      use modelp 
      use lindat 
      use optabl 
      use molist 
 
      implicit real(dp) (a-h,o-z),logical (l) 
 
      real(dp), parameter :: un=1.,ten15=1.e-15,c18=2.997925e18 
      real(dp) :: templ(mdepth)
! 
!     -------------- 
!     initialization 
!     -------------- 
! 
      igrdd=igrd 
      if(mode.eq.0) then 
! 
         read(2,*) ntemp,temp1,temp2 
         if(ntemp.eq.0) temp1=0.
         if(ntemp.eq.0.or.temp1.eq.0.) then
            mttab=mdepth
          else
            mttab=ntemp
         end if
!        write(*,*) 'ntemp,temp1,mttab',ntemp,temp1,mttab
         allocate (densg(mttab,mrtab),nden(mttab))

         read(2,*) idens 
         if(idens.lt.10) then 
            read(2,*) ndens,dens1,dens2 
          else if(idens.lt.20) then 
            read(2,*) ndens,densl1,densl2,densu1,densu2 
          else 
            do it=1,ntemp 
               read(2,*) ndens,densl,densu 
               densg(it,1)=densl 
               densg(it,ndens)=densu 
               nden(it)=ndens 
            end do 
         end if 
         if(idens.lt.20) then 
            do it=1,ntemp 
               nden(it)=ndens 
            end do 
         end if 
         if(ifeos.le.0) then 
            read(2,*) nfgrid,inttab,wlam1,wlam2 
            read(2,*) tabname,ibingr 
         end if
!
         write(6,"(/'INPUT PATRAMETERS FO R THE OPACITY TABLE',           &
     &              ' (from file fort.2)'/                                &
     &              '----------------------------------------'/           &
     &   'nfgrid  =',i8/                                                  &
     &   'inttab  =',i8/                                                  &
     &   'wlam1,2 =',2f10.1/                                              &
     &   'tabname =  ',a/)") nfgrid,inttab,wlam1,wlam2,tabname
      if(ntemp.eq.0) write(6,*) 'table = a model-tailored opacity table'
         if(ibingr.ne.0) write(6,*) 'table stored in the binary format'
! 
!        set mfgrid for the allocation routine to current nfgrid
!
         mfgrid=nfgrid
         call alloc_optabl 
! 
         irsct=0 
         irsche=0 
         irsch2=0 
! 
         wl1=log(wlam1) 
         wl2=log(wlam2) 
         dwl=(wl2-wl1)/(nfgrid-1) 
         do i=1,nfgrid 
            wlgrid(i)=exp(wl1+(i-1)*dwl) 
         end do 
! 
         if(temp1.gt.0.) then 
            at1=log(temp1) 
            at2=log(temp2) 
            dt=0. 
            if(ntemp.gt.1) dt=(at2-at1)/(ntemp-1) 
            do i=1,ntemp 
               templ(i)=at1+(i-1)*dt 
               tempg(i)=exp(templ(i)) 
            end do 
            if(idens.lt.10) then 
               at1=log(dens1) 
               at2=log(dens2) 
               dr=0. 
               ndens=nden(1) 
               if(ndens.gt.1) dr=(at2-at1)/(ndens-1) 
               do i=1,ntemp 
                  do j=1,ndens 
                     densg(i,j)=exp(at1+(j-1)*dr) 
                  end do 
               end do 
             else if(idens.lt.20) then 
               rhol1=log(densl1) 
               rhol2=log(densl2) 
               rhou1=log(densu1) 
               rhou2=log(densu2) 
               do i=1,ntemp 
                  ndens=nden(i) 
                  dens1=rhol1+(rhou1-rhol1)/(at2-at1)*(templ(i)-at1) 
                  dens2=rhol2+(rhou2-rhol2)/(at2-at1)*(templ(i)-at1) 
                  dr=0. 
                  if(ndens.gt.1) dr=(dens2-dens1)/(ndens-1) 
                  do j=1,ndens 
                     densg(i,j)=exp(dens1+(j-1)*dr) 
                  end do 
               end do 
             else 
               do i=1,ntemp 
                  ndens=nden(i) 
                  at1=log(densg(i,1)) 
                  at2=log(densg(i,ndens)) 
                  dr=0. 
                  if(ndens.gt.1) dr=(at2-at1)/(ndens-1) 
                  do j=2,ndens-1 
                     densg(i,j)=exp(at1+(j-1)*dr) 
                  end do 
               end do 
            end if 
! 
            write(6,"(/' COMPUTING AN OPACITY TABLE WITH GRID ',          & 
     &      'PARAMETERS:'/' ===== ntemp, ndens ',2i4)") ntemp,nden(1) 
            write(*,*) 'ntemp,ndens',ntemp,nden(1) 
            do i=1,ntemp 
               ndens=nden(i) 
               write(6,"(f10.1,20f8.2)")                                  & 
     &         tempg(i),(log10(densg(i,j)),j=1,ndens) 
            end do 
          else 
             if(inmod.eq.0) then 
                 call inkur 
               else 
                 call inpmod 
             end if 
 
            ntemp=nd 
            ndens=1 
!           write(*,*) 'ntemp,ndens',ntemp,ndens 
            do it=1,ntemp 
               tempg(it)=temp(it) 
               densg0(it)=dens(it) 
               densg(it,1)=dens(it) 
               elecm(it)=elec(it) 
               nden(it)=ndens 
!            write(*,*) it,tempg(it),densg(it,1),elecm(it),nden(it) 
            end do 
            if(ifeos.le.0) then 
               write(6,"(/' COMPUTING AN OPACITY TABLE WITH ',            & 
     &        'GRID PARAMETERS:'/' ===== ntemp, ndens ',2i4)")            & 
     &         ntemp,ndens 
               do it=1,ntemp 
                  write(6,"(f10.1,1pe12.3)") tempg(it),densg(it,1) 
               end do 
            end if 
            ndens=1 
            idens=2 
         end if 
! 
         nd=1 
         idstd=1 
         inext=1 
         frmx=0. 
         frmn=1.e20 
         idens0=mod(idens,10) 
! 
         indext=1 
         indexn=1 
         ipfreq=0 
         irelin=1 
         temp(1)=tempg(indext) 
! 
         write(6,"(/' ************************************',              & 
     &       /' GRID POINT OF THE OPACITY TABLE WITH:'/                   & 
     &        ' INDEX TEMP, T   ',i4,f10.1/                               & 
     &        ' INDEX DENS, DENS',I4,1PE10.1,                             & 
     &       /' ************************************'/)")                 & 
     &        indext,temp(1), indexn,densg(indext,indexn) 
! 
         if(temp1.le.0.) elec(1)=elecm(indext) 
         call densit(densg(indext,indexn),idens0) 
         if(ntemp.eq.1.and.ndens.eq.1) inext=0 
 
         elecgr(indext,indexn)=elec(1) 
 
         call abnchn(0) 
 
         return 
! 
!     --------------------------------------------- 
!     after computing the table for one T-rho pair: 
!     --------------------------------------------- 
! 
       else if(mode.eq.1) then 
         if(ifeos.le.0) then 
! 
         call timing(1,igrd+1) 
! 
         do i=1,3 
            xli(i)=0. 
         end do 
         do i=1,nmlist 
            xli(i)=float(nlinmt(i))*1.e-3 
         end do 
! 
         if(imode.ge.-5) then 
            if(indext.eq.1.and.indexn.eq.1)                               & 
     &      write(29,"('  it   ir     t       rho      elec',6x,          & 
     &      ' atomic   molec1  molec2  molec3      time'/)") 
            write(29,"(2i4,f9.2,1p2e10.2,2x,0pf8.1,2x,3f8.1,2x,f8.2)")    & 
     &              indext,indexn,temp(1),dens(1),elec(1),                & 
     &              float(nlin0)*1.e-3,                                   & 
     &              (xli(i),i=1,3),dtim 
          else 
            alam0=alam0s 
            if(alam0s.eq.0.) alam0=5.e7/temp(1)/10. 
            if(alam0s.lt.0.) alam0=-5.e7/temp(1)/alam0s 
            alast=alasts 
            if(alasts.eq.0.) alast=5.e7/temp(1)*20. 
            if(alasts.lt.0.) alast=-5.e7/temp(1)*alasts 
            if(alast.gt.1.e5) alast=1.e5 
            write(29,"(1p3e11.3,0pf9.3,0pf12.3)") temp(1),elec(1),        & 
     &         dens(1),alam0,alast 
         end if 
! 
!     ------------------------------------------------ 
!     interpolate and store previously computed table 
!     ------------------------------------------------ 
! 
         nfr=ipfreq 
         nfrtab(indext,indexn)=ipfreq 
!        write(*,*) 'indext,indexn,nfreq',indext,indexn,ipfreq 
!        write(*,"(/'NFR,NFGRID:',2i10)") nfr,nfgrid,inttab 
         if(inttab.eq.1) then 
!        call interp(wltab,absop,wlgrid,abgrd,nfr,nfgrid,2,0,0) 
            call intrp(wltab,absop,wlgrid,abgrd,nfr,nfgrid) 
          else 
            ij=0 
            ijgrd=0 
            grid: do 
               ijgrd=ijgrd+1 
               if(ijgrd.lt.nfgrid) then 
                  wlgr=0.5*(wlgrid(ijgrd)+wlgrid(ijgrd+1)) 
                else 
                  wlgr=wlgrid(nfgrid) 
               end if 
               isum=0 
               sum=0. 
               inn: do 
                  ij=ij+1 
                  if(ij.gt.nfr) exit grid 
                  wlt=wltab(ij) 
                  abl=absop(ij) 
                  if(wlt.le.wlgr) then 
                     sum=sum+exp(abl) 
                     isum=isum+1 
                     cycle inn 
                  end if 
                  if(isum.gt.0) then 
                     abgrd(ijgrd)=log(sum/float(isum)) 
                   else 
                     abg=abl+(absop(ij+1)-abl)/(wltab(ij+1)-wlt)*         & 
     &                   (wlgr-wlt) 
                     abgrd(ijgrd)=abg 
                  end if 
                  if(ijgrd.lt.nfgrid) then 
                     ij=ij-1 
                     cycle grid 
                   else 
                     exit grid 
!                   else if(ijgrd.eq.nfgrid) then 
!                    wlgr=wlgrid(nfgrid) 
!                    sum=0. 
!                    isum=0 
!                    if(ij.lt.nfr) ij=ij-1 
!                    cycle inn 
                  end if 
               end do inn 
            end do grid 
         end if 
! 
         do ij=1,nfgrid 
            absgrd(indext,indexn,ij)=real(abgrd(ij)) 
         end do 
         absgrd(indext,indexn,nfgrid)=absgrd(indext,indexn,nfgrid-1) 
      end if 
! 
!     ------------------------------ 
!     prepare values for a new table 
!     ------------------------------ 
! 
      ipfreq=0 
      ndens=nden(indext) 
      if(indexn.lt.ndens) then 
         indexn=indexn+1 
         rho=densg(indext,indexn) 
         write(6,"(/' ************************************'/              & 
     &        ' GRID POINT OF THE OPACITY TABLE WITH:'/                   & 
     &        ' INDEX TEMP, T   ',i4,f10.1/                               & 
     &        ' INDEX DENS, DENS',I4,1PE10.1/                             & 
     &        ' NFR, NFGRID     ',2i10/                                   & 
     &        ' ************************************'/)")                 & 
     &        indext,tempg(indext),indexn,densg(indext,indexn),           & 
     &        nfr,nfgrid 
         call densit(rho,idens0) 
         inext=1 
       else 
         indexn=1 
         irelin=1 
         if(indext.lt.ntemp) then 
            indext=indext+1 
            temp(1)=tempg(indext) 
            if(temp1.le.0.) then 
               densg(indext,indexn)=densg0(indext) 
               elec(1)=elecm(indext) 
            end if 
            rho=densg(indext,indexn) 
            write(6,"(/' ************************************'/           & 
     &        ' GRID POINT OF THE OPACITY TABLE WITH:'/                   & 
     &        ' INDEX TEMP, T   ',i4,f10.1/                               & 
     &        ' INDEX DENS, DENS',I4,1PE10.1/                             & 
     &        ' NFR, NFGRID     ',2i10/                                   & 
     &        ' ************************************'/)")                 & 
     &       indext,tempg(indext),indexn,densg(indext,indexn),            & 
     &       nfr,nfgrid 
            call densit(rho,idens0) 
            inext=1 
          else 
            inext=0 
         end if 
      end if 
      if(inext.eq.1) then 
         rewind(19) 
         if(inlist.lt.0) rewind(19) 
      end if 
! 
      elecgr(indext,indexn)=elec(1) 
! 
      call abnchn(0) 
      id=1 
         do i=1,4 
            do j=i+1,22 
               call hydtab(i,j,id) 
            end do 
         end do 
      end if 
! 
      return 
      end subroutine ingrid 
! 
! 
! ******************************************************************** 
! 
! 
      subroutine ougrid(abso) 
!     ======================= 
! 
!     output of grid opacities 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use optabl 
      implicit real(dp) (a-h,o-z) 
 
      real(dp), parameter :: un=1.,ten15=1.e-15,c18=2.997925e18 
      REAL(DP) :: ABSO(MFREQ) 
! 
      d1=un/dens(1) 
      if (nfreq.le.3) return 
! 
      if(iprin.lt.4) then 
         do ij=3,nfreq-1 
            abl=log(abso(ij)*d1) 
            ipfreq=ipfreq+1 
            absop(ipfreq)=abl 
            wltab(ipfreq)=2.997925e18/freq(ij) 
         end do 
       else 
         do ij=3,nfreq-1 
            abl=log(abso(ij)*d1) 
            ipfreq=ipfreq+1 
            write(27,"(i10,f14.5,0pf12.5)") ipfreq,c18/freq(ij),abl 
            absop(ipfreq)=abl 
            wltab(ipfreq)=2.997925e18/freq(ij) 
         end do 
      end if 
! 
      return 
      end subroutine ougrid 
! 
! 
! ******************************************************************** 
! 
! 
      subroutine fingrd 
!     ================= 
! 
!     storing the complete, interpolated, opacity table 
! 
      use accura 
      use params 
      use modelp 
      use synthp 
      use optabl 
      implicit real(dp) (a-h,o-z) 
 
! 
      if(ifeos.gt.0) return 
! 
      close(53) 
      iophmp=iophmi 
      if(ielhm.gt.0.and.relabn(1).gt.0.) iophmp=1 
      if(ibingr.eq.0) then 
         open(53,file=tabname,status='unknown') 
         write(53,"('opacity table with element abundances:'/             & 
     &          'element   for EOS   for opacities')") 
         do iat=1,92 
            write(53,"('  ',a4,1p2e12.3)")                                & 
     &         typat(iat),abnd(iat),abnd(iat)*relabn(iat) 
         end do 
         write(53,"(/'molecules - ifmol,tmolim:'/,i4,f10.1)")             & 
     &      ifmol,tmolim 
         write(53,"('additional opacities'/                               & 
     &   '  H-  H2+ He- CH  OH  H2- CIA: H2H2 H2He H2H  HHe'/             & 
     &          6i4,4x,4i4)") iophmp,ioph2p,iophem,iopch,iopoh,           & 
     &          ioph2m,ioh2h2,ioh2he,ioh2h1,iohhe 
         if(idens.lt.10) then 
            ndens=nden(1) 
            write(53,"(/'number of frequencies, temperatures, ',          & 
     &        'densities:'/10x,3i10)") nfgrid,ntemp,nden(1) 
            write(53,"('log temperatures'/(6F11.6))")                     & 
     &         (log(tempg(i)),i=1,ntemp) 
            write(53,"('log densities'/(6F11.6))")                        & 
     &         (log(densg(1,j)),j=1,nden(1)) 
            write(53,"('log electron densities from EOS'/(6f11.6))")      & 
     &      ((log(elecgr(i,j)),j=1,nden(1)),i=1,ntemp) 
            do k = 1, nfgrid 
               write(53,"(/' *** frequency # : ',i8,f15.5/1pe20.8)")      & 
     &         k,wlgrid(k),2.997925e18/wlgrid(k) 
               do j = 1,ndens 
                  write(53,"((1p6e14.6))") (absgrd(i,j,k),i=1,ntemp) 
               end do 
            end do 
          else 
            write(53,"(/'number of frequencies, temperatures, ',          & 
     &         'densities:'/10x,3i10)") nfgrid,ntemp,-nden(1) 
            write(53,"(30i3)") (nden(i),i=1,ntemp) 
            write(53,"('log temperatures'/(6F11.6))") (log(tempg(i)),     & 
     &         i=1,ntemp) 
            write(53,"('log densities')") 
            do i=1,ntemp 
               ndens=nden(i) 
               write(53,"(6f14.6)") (log(densg(i,j)),j=1,ndens) 
            end do 
            write(53,"('log electron densities from EOS')") 
            do i=1,ntemp 
               ndens=nden(i) 
               write(53,"(6f14.6)") (log(elecgr(i,j)),j=1,ndens) 
            end do 
            do k = 1,nfgrid 
               write(53,"(/' *** frequency # : ',i8,f15.5/1pe20.8)")      & 
     &         k,wlgrid(k),2.997925e18/wlgrid(k) 
               do i=1,ntemp 
                  ndens=nden(i) 
                  write(53,"((1p6e14.6))") (absgrd(i,j,k),j=1,ndens) 
               end do 
            end do 
         end if 
       end if 
         do iat=1,92 
            write(63) typat(iat),abnd(iat),abnd(iat)*relabn(iat) 
         end do 
         write(63) ifmol,tmolim 
         write(63) iophmp,ioph2p,iophem,iopch,iopoh,ioph2m,               & 
     &             ioh2h2,ioh2he,ioh2h1,iohhe 
 
         if(idens.lt.10) then 
            ndens=nden(1) 
            write(63) nfgrid,ntemp,nden(1) 
            write(63) (log(tempg(i)),i=1,ntemp) 
            write(63) (log(densg(1,j)),j=1,nden(1)) 
            write(63) ((log(elecgr(i,j)),j=1,nden(1)),i=1,ntemp) 
            do k = 1, nfgrid 
               write(63) 2.997925e18/wlgrid(k) 
               do j = 1,ndens 
                  write(63) (absgrd(i,j,k),i=1,ntemp) 
               end do 
            end do 
          else 
            write(63) nfgrid,ntemp,-nden(1) 
            write(63) (nden(i),i=1,ntemp) 
            write(63) (log(tempg(i)),i=1,ntemp) 
            do i=1,ntemp 
               ndens=nden(i) 
               write(63) (log(densg(i,j)),j=1,ndens) 
            end do 
            do i=1,ntemp 
               ndens=nden(i) 
               write(63) (log(elecgr(i,j)),j=1,ndens) 
            end do 
            do k = 1,nfgrid 
               write(63) 2.997925e18/wlgrid(k) 
               do i=1,ntemp 
                  ndens=nden(i) 
                  write(63) (absgrd(i,j,k),j=1,ndens) 
                  if(k.le.100) write(*,*) 'abs(1)',i,ndens,               & 
     &              (absgrd(i,j,k),j=1,ndens) 
               end do 
            end do 
         end if 
!     end if 
! 
      close(53) 
      close(63) 
      return 
      end subroutine fingrd 
! 
! 
!     ************************************************************* 
! 
! 
      subroutine abnchn(mode) 
!     ======================= 
! 
!     changing abundances (eliminating) species for an 
!     evaluating an opacity table 
! 
      use accura 
      use params 
      use modelp 
      use optabl, only : relabn 
      implicit real(dp) (a-h,o-z) 
 
      data iread/1/ 
! 
      if(iread.eq.1) then 
         do ia=1,matom 
            relabn(ia)=1. 
         end do 
         do 
            read(2,*,iostat=ios) iatom,rela 
            if(ios.ne.0) exit 
            relabn(iatom)=rela 
            write(*,*) 'ABUNDANCES CHANGED (AT.NUMBER, ABUND):',          & 
     &         iatom,rela 
         end do 
         if(relabn(1).eq.0.) then 
            iophmi=0 
            ioph2p=0 
         end if 
         iread=0 
      end if 
! 
      if(mode.eq.0) then 
         do iat=1,natom 
            do ii=n0a(iat),nka(iat) 
               popul0(ii,1)=popul(ii,1) 
            end do 
         end do 
         return 
      end if 
! 
      do iat=1,natom 
         ia=numat(iat) 
         do ii=n0a(iat),nka(iat) 
            popul(ii,1)=popul0(ii,1)*relabn(ia) 
         end do 
      end do 
! 
      do ia=1,matom 
         do io=1,mion0 
            rrr(1,io,ia)=rrr(1,io,ia)*relabn(ia) 
         end do 
      end do 
! 
      return 
      end subroutine abnchn 
! 
! 
!     ************************************************************* 
! 
! 
 
      subroutine densit(rho,idens) 
!     ============================ 
! 
!     determining the state parameters for the opacity grid 
!     calculations 
! 
      use accura 
      use params 
      use modelp 
      implicit real(dp) (a-h,o-z) 
 
      id=1 
      dm(id)=0. 
      IF(IFMOL.EQ.0.OR.TEMP(ID).GT.TMOLIM)                                & 
     &   WMM(ID)=WMY(ID)*HMASS/YTOT(ID) 
         if(idens.eq.0) then 
            ELEC(ID)=rho 
            ane=elec(id) 
            call todens(id,temp(id),an,ane) 
            DENS(ID)=(an-ane)*wmm(id) 
            p=an*bolk*temp(id) 
            WRITE(6,"(' **densit** t,rho,ne',I3,0PF10.1,1P5e11.3)")       & 
     &         ID,TEMP(ID),DENS(ID),ELEC(ID) 
          else if(idens.lt.0) then 
            AN=rho/TEMP(ID)/BOLK 
            CALL ELDENS(ID,TEMP(ID),AN,ANE,ANH,ANP) 
            ELEC(ID)=ANE 
            DENS(ID)=WMM(ID)*(AN-ELEC(ID)) 
            WRITE(6,"(' **densit** t,rho,ne,rho0,an',I3,0PF10.1,          & 
     &         1P5e11.3)") ID,TEMP(ID),DENS(ID),ELEC(ID),ane0,an 
          else if(idens.eq.1) then 
            DENS(ID)=RHO 
            CALL RHONEN(ID,TEMP(ID),RHO,AN,ANE) 
            ELEC(ID)=ANE 
            DENS(ID)=RHO 
            rho0=WMM(ID)*(AN-ANE) 
            WRITE(6,"(' **densit** t,rho,ne,rho0,an',I3,0PF10.1,          & 
     &         1P5e11.3)") IDens,TEMP(ID),DENS(ID),ane,rho0,an 
          else if(idens.eq.2) then 
            CALL RHONEN(ID,TEMP(ID),RHO,AN,ANE) 
            DENS(ID)=RHO 
            ANE=ELEC(ID) 
            rho0=WMM(ID)*(AN-ANE) 
            WRITE(6,"(' **densit** t,rho,ne,rho0,an',I3,0PF10.1,          & 
     &         1P5e11.3)") idens,TEMP(ID),DENS(ID),ane,rho0,an 
         end if 
      CALL INIMOD 
! 
      CALL WNSTOR(ID) 
      CALL SABOLF(ID) 
      CALL RATMAT(ID,ESEMAT,BESE) 
      CALL LEVSOL(ESEMAT,BESE,POPLTE,NLEVEL) 
      DO J=1,NLEVEL 
         POPUL(J,ID)=POPLTE(J) 
      END DO 
! 
      return 
      end subroutine densit 
 
 
! 
! ******************************************************************** 
! 
 
      SUBROUTINE TODENS(ID,T,AN,ANE) 
!     ============================== 
! 
!     determines AN (and ANP, AHTOT, and AHMOL) from T and ANE 
! 
!     Input parameters: 
!     T    - temperature 
!     ANE  - electron number density 
! 
!     Output: 
!     AN    - total particle density 
!     ANP   - proton number density 
!     AHTOT - total hydrogen number density 
!     AHMOL - relative number of hydrogen molecules with respect to the 
!             total number of hydrogens 
! 
      use accura 
      use params 
      use modelp 
      use eospar 
      implicit real(dp) (a-h,o-z) 
 
      real(dp), parameter :: un=1.,two=2.,half=0.5 
! 
      QM=0. 
      Q2=0. 
      QP=0. 
      Q=0. 
      DQN=0. 
      TK=BOLK*T 
      THET=5.0404e3/T 
! 
!     Coefficients entering ionization (dissociation) balance of: 
!     atomic hydrogen          - QH; 
!     negative hydrogen ion    - QM 
!     hydrogen molecule        - QP 
!     ion of hydrogen molecule - Q2 
! 
      QM=1.0353e-16/T/SQRT(T)*EXP(8762.9/T) 
      QH=EXP((15.38287+1.5*LOG10(T)-13.595*THET)*2.30258509299405) 
! 
      if(t.gt.16000.) then 
         ih2=0 
         ih2p=0 
       else 
         QP=TK*EXP((-11.206998+THET*(2.7942767+THET*                      & 
     &      (0.079196803-0.024790744*THET)))*2.30258509299405) 
         Q2=TK*EXP((-12.533505+THET*(4.9251644+THET*                      & 
     &      (-0.056191273+0.0032687661*THET)))*2.30258509299405) 
         ih2=1 
      end if 
 
! 
!     procedure STATE determines Q (and DQN) - the total charge (and its 
!     derivative wrt temperature) due to ionization of all atoms which 
!     are considered (both explicit and non-explicit), by solving the set 
!     of Saha equations for the current values of T and ANE 
! 
      CALL STATE(ID,T,ANE,Q) 
! 
!     Auxiliary parameters for evaluating the elements of matrix of 
!     linearized equations. 
!     Note that complexity of the matrix depends on whether the hydrogen 
!     molecule is taken into account 
!     Treatment of hydrogen ionization-dissociation is based on 
!     Mihalas, in Methods in Comput. Phys. 7, p.10 (1967) 
! 
      G2=QH/ANE 
      G3=0. 
      G4=0. 
      G5=0. 
      D=0. 
      E=0. 
      G3=QM*ANE 
      A=UN+G2+G3 
      D=G2-G3 
         IF(IH2.EQ.0) THEN 
            F1=UN/A 
            FE=D/A+Q 
          ELSE 
            E=G2*QP/Q2 
            B=TWO*(UN+E) 
            GG=ANE*Q2 
            C1=B*(GG*B+A*D)-E*A*A 
            C2=A*(TWO*E+B*Q)-D*B 
            C3=-E-B*Q 
            F1=(SQRT(C2*C2-4.*C1*C3)-C2)*HALF/C1 
            FE=F1*D+E*(UN-A*F1)/B+Q 
         END IF 
         AH=ANE/FE 
         ANH=AH*F1 
      AE=ANH/ANE 
      GG=AE*QP 
      E=ANH*Q2 
      B=ANH*QM 
! 
!      S(1)=AN-ANE-YTOT(ID)*AH 
!      S(2)=ANH*(D+GG)+Q*AH-ANE 
!      S(3)=AH-ANH*(A+TWO*(E+GG)) 
! 
      hhn=A+TWO*(E+GG) 
      anh=ane/(d+gg+q*hhn) 
      ah=anh*hhn 
      an=ane+ytot(id)*ah 
! 
      AHTOT=AH 
      AHMOL=TWO*ANH*(ANH*Q2+ANH/ANE*QP)/AH 
      ANP=ANH/ANE*QH 
      RETURN 
      END SUBROUTINE TODENS 
! 
! 
! *********************************************************************** 
! 
 
      subroutine rhonen(id,t,rho,an,ane) 
!     ================================== 
! 
!     iterative determination of N and Ne from given T and RHO 
! 
! 
!     Input:  T   - temperature 
!             RHO - mass density 
!     Output: AN  - total particle density 
!             ANE - elctron density 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z) 
 
      it=0 
      if(id.eq.1.and.anerel.eq.0.) then 
         anerel=0.5 
         if(t.lt.9000.) anerel=0.4 
         if(t.lt.8000.) anerel=0.1 
         if(t.lt.7000.) anerel=0.01 
         if(t.lt.6000.) anerel=0.001 
         if(t.lt.5500.) anerel=0.0001 
         if(t.lt.5000.) anerel=1.e-5 
         if(t.lt.4000.) anerel=1.e-6 
         if(it.lt.1000.) anerel=1.e-10 
      end if 
      do 
         it=it+1 
         an=rho/wmm(id)/(1.-anerel) 
         ane0=anerel*an 
         call eldens(id,t,an,ane,anh,anp) 
         anerel=ane/an 
!!       write(6,"(/' **** rhonen it,id,t,r,N,Ne,wmm,ner',2i4,f7.0, 
!!   *    1p5e11.4)") it,id,t,rho,an,ane,wmm(id),anerel 
         if(abs((ane-ane0)/ane0).lt.1.e-5) exit 
         if(it.lt.50) cycle 
      end do 
! 
      return 
      end subroutine rhonen 
! 
! ******************************************************************** 
! 
 
      SUBROUTINE ELDENS(ID,T,AN,ANE,ANH,ANP) 
!     ====================================== 
! 
!     Evaluation of the electron density and the total hydrogen 
!     number density for a given total particle number density 
!     and temperature; 
!     by solving the set of Saha equations, charge conservation and 
!     particle conservation equations (by a Newton-Raphson method) 
! 
!     Input parameters: 
!     T    - temperature 
!     AN   - total particle number density 
! 
!     Output: 
!     ANE   - electron density 
!     ANP   - proton number density 
!     AHTOT - total hydrogen number density 
!     AHMOL - relativer number of hydrogen molecules with respect to the 
!             total number of hydrogens 
!     ENERG - part of the internal energy: excitation and ionization 
! 
      use accura 
      use params 
      use modelp 
      use eospar 
      implicit real(dp) (a-h,o-z) 
 
      real(dp), parameter :: un=1.,two=2.,half=0.5 
      REAL(DP) :: R(3,3),S(3),SOL(3) 
! 
      TK=BOLK*T 
      if(ifmol.gt.0.and.t.lt.tmolim) then 
         aein=an*anerel 
         call moleq(id,t,an,aein,ane,0) 
         anerel=ane/an 
         return 
      end if 
! 
      QM=0. 
      Q2=0. 
      QP=0. 
      Q=0. 
      DQN=0. 
      TK=BOLK*T 
      THET=5.0404e3/T 
! 
!     Coefficients entering ionization (dissociation) balance of: 
!     atomic hydrogen          - QH; 
!     negative hydrogen ion    - QM 
!     hydrogen molecule        - Q2 
!     ion of hydrogen molecule - QP 
! 
      IF(IATREF.EQ.IATH) THEN 
         QM=1.0353e-16/T/SQRT(T)*EXP(8762.9/T) 
         QH0=EXP((15.38287+1.5*LOG10(T)-13.595*THET)*2.30258509299405) 
! 
         if(t.gt.16000.) then 
            ih2=0 
          else 
            ih2=1 
            QP=TK*EXP((-11.206998+THET*(2.7942767+THET*                   & 
     &      (0.079196803-0.024790744*THET)))*2.30258509299405) 
            Q2=TK*EXP((-12.533505+THET*(4.9251644+THET*                   & 
     &      (-0.056191273+0.0032687661*THET)))*2.30258509299405) 
         end if 
      END IF 
! 
!     Initial estimate of the electron density 
! 
      if(anerel.le.0.) then 
         if(t.gt.1.e4) then 
            anerel=0.5 
          else 
            if(elec(id).gt.0..and.dens(id).gt.0.) then 
               anerel=elec(id)/(elec(id)+dens(id)/wmm(id)) 
             else 
               anerel=0.1 
            end if 
         end if 
      end if 
! 
      ANE=AN*ANEREL 
      IT=0 
! 
!     Basic Newton-Raphson loop - solution of the non-linear set 
!     for the unknown vector P, consistiong of AH, ANH (neutral 
!     hydrogen number density) and ANE. 
! 
      DO 
         IT=IT+1 
! 
!        procedure STATE determines Q (and DQN) - the total charge (and its 
!        derivative wrt temperature) due to ionization of all atoms which 
!        are considered (both explicit and non-explicit), by solving the set 
!        of Saha equations for the current values of T and ANE 
! 
         CALL STATE(ID,T,ANE,Q) 
         QH=QH0*2./PFSTD(1,1) 
! 
!        Auxiliary parameters for evaluating the elements of matrix of 
!        linearized equations. 
!        Note that complexity of the matrix depends on whether the hydrogen 
!        molecule is taken into account 
!        Treatment of hydrogen ionization-dissociation is based on 
!        Mihalas, in Methods in Comput. Phys. 7, p.10 (1967) 
! 
         IF(IATREF.EQ.IATH) THEN 
            G2=QH/ANE 
            G3=0. 
            G4=0. 
            G5=0. 
            D=0. 
            E=0. 
            G3=QM*ANE 
            A=UN+G2+G3 
            D=G2-G3 
            IF(IT.LE.1) THEN 
               IF(IH2.EQ.0) THEN 
                  F1=UN/A 
                  FE=D/A+Q 
                ELSE 
                  E=G2*QP/Q2 
                  B=TWO*(UN+E) 
                  GG=ANE*Q2 
                  C1=B*(GG*B+A*D)-E*A*A 
                  C2=A*(TWO*E+B*Q)-D*B 
                  C3=-E-B*Q 
                  F1=(SQRT(C2*C2-4.*C1*C3)-C2)*HALF/C1 
                  FE=F1*D+E*(UN-A*F1)/B+Q 
               END IF 
               AH=ANE/FE 
               ANH=AH*F1 
            END IF 
            AE=ANH/ANE 
            GG=AE*QP 
            E=ANH*Q2 
            B=ANH*QM 
! 
!           Matrix of the linearized system R, and the rhs vector S 
! 
            R(1,1)=YTOT(ID) 
            r(1,2)=-two*(anh*q2+gg) 
            R(1,3)=UN 
            R(2,1)=-Q 
            R(2,2)=-D-TWO*GG 
            R(2,3)=UN+B+AE*(G2+GG)-DQN*AH 
            R(3,1)=-UN 
            R(3,2)=A+4.*(anh*q2+GG) 
            R(3,3)=B-AE*(G2+TWO*GG) 
            S(1)=AN-ANE-YTOT(ID)*AH+anh*(anh*q2+gg) 
            S(2)=ANH*(D+GG)+Q*AH-ANE 
            S(3)=AH-ANH*(A+TWO*(anh*q2+GG)) 
! 
!           Solution of the linearized equations for the correction vector P 
! 
            CALL LINEQS(R,S,SOL,3,3) 
! 
!           New values of AH, ANH, and ANE 
! 
            AH=AH+SOL(1) 
            ANH=ANH+SOL(2) 
            DELNE=SOL(3) 
            ANE=ANE+DELNE 
! 
!        hydrogen is not the reference atom 
! 
         ELSE 
! 
!           Matrix of the linearized system R, and the rhs vector S 
! 
            IF(IT.EQ.1) THEN 
               ANE=AN*HALF 
               AH=ANE/YTOT(ID) 
            END IF 
 
!C TEMPORARAY 
            qref=0. 
            dqn=0. 
            dqnr=0. 
!C END TEMPORARY 
 
            R(1,1)=YTOT(ID) 
            R(1,2)=UN 
            R(2,1)=-Q-QREF 
            R(2,2)=UN-(DQN+DQNR)*AH 
            S(1)=AN-ANE-YTOT(ID)*AH 
            S(2)=(Q+QREF)*AH-ANE 
! 
!        Solution of the linearized equations for the correction vector P 
! 
            CALL LINEQS(R,S,SOL,2,3) 
            AH=AH+SOL(1) 
            DELNE=SOL(2) 
            ANE=ANE+DELNE 
         END IF 
! 
!        Convergence criterion 
! 
         IF(ANE.LE.0.) ANE=1.e-7*AN 
         IF(ABS(DELNE/ANE).LT.1.e-6.OR.IT.GE.20) EXIT 
      END DO 
! 
!     ANEREL is the exact ratio betwen electron density and total 
 
!     particle density, which is going to be used in the subseguent 
!     call of ELDENS 
! 
      ANEREL=ANE/AN 
      AHTOT=AH 
      IF(IATREF.EQ.IATH) THEN 
         AHMOL=ANH*ANH*Q2 
         ANP=ANH/ANE*QH 
         ANHMI=ANH*ANE*QM 
         anhn=anh+anp+anhmi+2.*ahmol 
         wmm(id)=wmy(id)/(ytot(id)-ahmol/anhn)*hmass 
         ahn(id)=anh 
         ahp(id)=anp 
      END IF 
! 
      RETURN 
      END SUBROUTINE ELDENS 
! 
 
 
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE TIMING(MOD,ITER) 
!     =========================== 
! 
!     Timing procedure (call machine dependent routine!!) 
! 
      use params, only : dtim 
 
      CHARACTER(LEN=6)    :: ROUT 
      REAL                :: DUMMY(2),TIME,DT 
      REAL, SAVE          :: T0=0. 
!!    REAL, SAVE          :: DTIM 
! 
!     TIME=etime(dummy) 
      call etime(dummy,time) 
 
      DT=TIME-T0 
      T0=TIME 
      IP=ITER 
      IF(MOD.EQ.1) THEN 
         ROUT=' TABLE' 
      ELSE IF(MOD.EQ.2) THEN 
         ROUT=' FINAL' 
      ENDIF 
      WRITE(69,"(I6,2F11.2,2X,A6)") IP,TIME,DT,ROUT 
      dtim=dt 
      RETURN 
      END SUBROUTINE TIMING 
! 
! 
! ******************************************************************** 
! 
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
      if(ifeos.lt.0) istp=-ifeos 
! 
      write(*,*) 'id,nd,istp,idstd',id,nd,istp,idstd 
      do id=1,nd,istp 
         t=temp(id) 
         ane=elec(id) 
         rho=dens(id) 
         ann = dens(id)/wmm(id)+elec(id) 
! 
         if(ifmol.eq.0.or.t.gt.tmolim*1.e5) then 
            it=0 
            itera: do 
               ann0=ann 
               it=it+1 
            write(*,*) 'BE:id,it,t,ann,ane',id,it,t,ann,ane 
               call eldens(id,t,ann,ane,anh,anp) 
            write(*,*) 'AF:id,it,t,ann,ane',id,it,t,ann,ane 
            write(*,*) 'id,it,anh,anp',id,it,anh,anp 
               anmol(1,id)=anhmi 
               anmol(2,id)=ahmol 
               anato(1,id)=anh 
               anion(1,id)=anp 
               hpop=dens(id)/wmy(id)/hmass 
           write(*,*) 'anh,hpop',anh,hpop,nmetal 
           write(*,*) 'neleme',neleme(1) 
               do i=1,nmetal 
                  j=neleme(i) 
 
                  anato(j,id)=anato(j,id)*hpop 
                  anion(j,id)=anion(j,id)*hpop 
             write(*,*) i,j,id,anato(j,id),anion(j,id) 
                  if(j.ge.2.and.j.le.30) anion2(j,id)=anion2(j,id)*hpop 
               end do 
               anato(1,id)=anh 
               anion(1,id)=anp 
               wmm(id)=wmy(id)/(ytot(id)-anmol(2,id)/hpop)*hmass 
               ann=dens(id)/wmm(id)+ane 
               if((ann-ann0)/ann0.le.1.e-5.or.it.gt.20) exit 
            end do itera 
         end if 
! 
         nmetae=38 
         write(6,"(' **** DEPTH ID',i4,'  T =',f10.1/)") id,temp(id) 
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
! 
! ******************************************************************* 
! 
! 
 
      subroutine cia_h2h2(t,ah2,ff,opac) 
!     ===================--============= 
! 
!     CIA H2-H2 opacity 
!     data from Borysow A., Jorgensen U.G., Fu Y. 2001, JQSRT 68, 235 
! 
      use accura 
      implicit real(dp)(A-H,O-Z) 
      integer, parameter :: nlines=1000 
      real(dp) :: freq(nlines),temp(7),alpha(nlines,7) 
      real(dp), parameter :: amagat=2.6867774e+19,fac=1./amagat**2 
      data temp / 1000. , 2000. , 3000. , 4000. , 5000. , 6000. ,         & 
     &            7000. / 
      data ntemp /7/ 
      data ifirst /0/ 
      PARAMETER (CAS=2.997925e10) 
      SAVE FREQ,ALPHA 
 
!     input frequency in Hz but needed wave numbers in cm^-1 
      f=ff/cas 
!     read in CIA tables if this is the first call 
      if (ifirst.eq.0) then 
         write(*,'(a)') 'Reading in H2-H2 CIA opacity tables...' 
         open(10,file="./data/CIA_H2H2.dat",status='old') 
         do i=1,3 
            read (10,*) 
         enddo 
         do i=1,nlines 
            read (10,*) freq(i),(alpha(i,j),j=1,ntemp) 
         enddo 
         close(10) 
 
!     take logarithm of tables prior to doing linear interpolations 
 
         do i=1,nlines 
            do j=1,ntemp 
               alpha(i,j)=log(alpha(i,j)) 
            enddo 
         enddo 
 
         ifirst=1 
      endif 
 
!     locate position in temperature array 
      call locate(temp,ntemp,t,j,ntemp) 
 
      if (j.eq.0) then 
         write(66,'(a,f6.0,a)')                                           & 
     &   'Warning: requested temperature is below',temp(1),' K' 
         write(66,'(a)') 'CIA H2-H2 opacity is extreapolated down' 
         j=1 
!        opac=0. 
!        return 
      endif 
 
!     locate position in frequency array 
      call locate(freq,nlines,f,i,nlines) 
 
!     linearly interpolate in frequency and temperature 
 
      if (j.eq.ntemp) then 
!     hold values constant if off high temperature end of table 
         y1=alpha(i,j) 
         y2=alpha(i+1,j) 
         tt=(f-freq(i))/(freq(i+1)-freq(i)) 
         alp=(1.-tt)*y1 + tt*y2 
      else if (i.eq.0 .or. i.eq.nlines) then 
!     set values to a very small number if off frequency table 
         alp=-50. 
      else 
!     interpolate linearly within table 
         y1=alpha(i,j) 
         y2=alpha(i+1,j) 
         y3=alpha(i+1,j+1) 
         y4=alpha(i,j+1) 
 
         tt=(f-freq(i))/(freq(i+1)-freq(i)) 
         uu=(t-temp(j))/(temp(j+1)-temp(j)) 
 
         alp=(1.-tt)*(1.-uu)*y1 + tt*(1.-uu)*y2 + tt*uu*y3 +              & 
     &       (1.-tt)*uu*y4 
      endif 
 
      alp=exp(alp) 
 
!     final opacity 
 
      opac=fac*ah2*ah2*alp 
! 
      return 
      end subroutine cia_h2h2 
 
! 
! 
! 
! ******************************************************************** 
! 
! 
 
      SUBROUTINE locate(xx,n,x,j,nxdim) 
!     ================================= 
! 
      use accura 
      implicit real(dp)(A-H,O-Z) 
      real(dp) :: xx(nxdim) 
! 
      jl=0 
      ju=n+1 
      loca: do 
         if(ju-jl.gt.1) then 
           jm=(ju+jl)/2 
           if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then 
             jl=jm 
            else 
             ju=jm 
           end if 
           cycle loca 
          else 
           exit loca 
         end if 
      end do loca 
      if(x.eq.xx(1)) then 
         j=1 
       else if(x.eq.xx(n)) then 
         j=n-1 
       else 
         j=jl 
      end if 
      return 
      END SUBROUTINE locate 
 
! 
! 
! ******************************************************************** 
! 
! 
 
      subroutine cia_h2he(t,ah2,ahe,ff,opac) 
!     ====================================== 
! 
!     CIA H2-He opacity 
!     data from Jorgensen U.G., Hammer D., Borysow A., Falkesgaard J., 2000, 
!     Astronomy & Astrophysics 361, 283 
! 
      use accura 
      implicit real(dp)(A-H,O-Z) 
      integer, parameter :: nlines=242 
      real(dp), save :: freq(nlines),temp(7),alpha(nlines,7) 
      real(dp), parameter :: amagat=2.6867774e+19,fac=1./amagat**2 
      data temp / 1000. , 2000. , 3000. , 4000. , 5000. , 6000. ,         & 
     &            7000. / 
      data ntemp /7/ 
      data ifirst /0/ 
      REAL(DP), PARAMETER :: CAS=2.997925e10 
 
!     input frequency in Hz but needed wave numbers in cm^-1 
      f=ff/cas 
!     read in CIA tables if this is the first call 
      if (ifirst.eq.0) then 
         write(*,'(a)') 'Reading in H2-He CIA opacity tables...' 
         open(10,file="./data/CIA_H2He.dat",status='old') 
         do i=1,3 
            read (10,*) 
         end do 
         do i=1,nlines 
            read (10,*) freq(i),(alpha(i,j),j=1,ntemp) 
         end do 
         close(10) 
 
!     take logarithm of tables prior to doing linear interpolations 
 
         do i=1,nlines 
            do j=1,ntemp 
               alpha(i,j)=log(alpha(i,j)) 
            end do 
         end do 
 
         ifirst=1 
      end if 
 
!     locate position in temperature array 
      call locate(temp,ntemp,t,j,ntemp) 
 
      if(j.eq.0) then 
         write(66,'(a,f6.0,a)')                                           & 
     &   'Warning: requested temperature is below',temp(1),' K' 
         write(66,'(a)') 'CIA H2-He opacity is extrapolated down' 
!        opac=0. 
!        return 
         j=1 
      end if 
 
!     locate position in frequency array 
      call locate(freq,nlines,f,i,nlines) 
 
!     linearly interpolate in frequency and temperature 
 
      if(j.eq.ntemp) then 
!        hold values constant if off high temperature end of table 
         y1=alpha(i,j) 
         y2=alpha(i+1,j) 
         tt=(f-freq(i))/(freq(i+1)-freq(i)) 
         alp=(1.-tt)*y1 + tt*y2 
       else if (i.eq.0 .or. i.eq.nlines) then 
!        set values to a very small number if off frequency table 
         alp=-50. 
       else 
!        interpolate linearly within table 
         y1=alpha(i,j) 
         y2=alpha(i+1,j) 
         y3=alpha(i+1,j+1) 
         y4=alpha(i,j+1) 
 
         tt=(f-freq(i))/(freq(i+1)-freq(i)) 
         uu=(t-temp(j))/(temp(j+1)-temp(j)) 
 
         alp=(1.-tt)*(1.-uu)*y1 + tt*(1.-uu)*y2 + tt*uu*y3 +              & 
     &       (1.-tt)*uu*y4 
      end if 
 
      alp=exp(alp) 
 
!     final opacity 
 
      opac=fac*ah2*ahe*alp 
! 
      return 
      end subroutine cia_h2he 
! 
! 
! ******************************************************************* 
! 
! 
 
      subroutine cia_h2h(t,ah2,ah,ff,opac) 
!     ==================================== 
! 
!     CIA H2-H opacity - data taken from TURBOSPEC 
! 
      use accura 
      implicit real(dp)(A-H,O-Z) 
      integer, parameter :: nlines=67 
      real(dp), save :: freq(nlines),alpha(nlines,4) 
      real(dp)       :: temp(4) 
      real(dp), parameter :: amagat=2.6867774e+19,fac=1./amagat**2 
      data temp / 1000. , 1500., 2000. , 2500. / 
      data ntemp /4/ 
      data ifirst /0/ 
      REAL(DP), PARAMETER :: CAS=2.997925e10 
 
!     input frequency in Hz but needed wave numbers in cm^-1 
      f=ff/cas 
 
!     read in CIA tables if this is the first call 
      if(ifirst.eq.0) then 
         write(*,'(a)') 'Reading in H2-H CIA opacity tables...' 
         open(10,file="./data/CIA_H2H.dat",status='old') 
         do i=1,3 
            read (10,*) 
         end do 
         do i=1,nlines 
            read (10,*) freq(i),(alpha(i,j),j=1,ntemp) 
         end do 
         close(10) 
 
!     take logarithm of tables prior to doing linear interpolations 
 
         do i=1,nlines 
            do j=1,ntemp 
               alpha(i,j)=log(alpha(i,j)) 
            enddo 
         enddo 
 
         ifirst=1 
      end if 
 
!     locate position in temperature array 
      call locate(temp,ntemp,t,j,ntemp) 
 
      if(j.eq.0) then 
         write(66,'(a,f6.0,a)')                                           & 
     &   'Warning: requested temperature is below',temp(1),' K' 
         write(66,'(a)') 'CIA H2-H opacity id extraspolsated down' 
!        opac=0. 
!        return 
         j=1 
      end if 
 
!     locate position in frequency array 
      call locate(freq,nlines,f,i,nlines) 
 
!     linearly interpolate in frequency and temperature 
 
      if (j.eq.ntemp) then 
!        hold values constant if off high temperature end of table 
         y1=alpha(i,j) 
         y2=alpha(i+1,j) 
         tt=(f-freq(i))/(freq(i+1)-freq(i)) 
         alp=(1.-tt)*y1 + tt*y2 
       else if (i.eq.0 .or. i.eq.nlines) then 
!        set values to a very small number if off frequency table 
         alp=-50. 
       else 
!        interpolate linearly within table 
         y1=alpha(i,j) 
         y2=alpha(i+1,j) 
         y3=alpha(i+1,j+1) 
         y4=alpha(i,j+1) 
 
         tt=(f-freq(i))/(freq(i+1)-freq(i)) 
         uu=(t-temp(j))/(temp(j+1)-temp(j)) 
 
         alp=(1.-tt)*(1.-uu)*y1 + tt*(1.-uu)*y2 + tt*uu*y3 +              & 
     &       (1.-tt)*uu*y4 
      end if 
 
      alp=exp(alp) 
 
!     final opacity 
 
      opac=fac*ah2*ah*alp 
! 
      return 
      end subroutine cia_h2h 
! 
! 
! ******************************************************************* 
! 
! 
 
      subroutine cia_hhe(t,ah,ahe,ff,opac) 
!     ==================================== 
! 
!     CIA H-He opacity 
!     data from Gustafsson M., Frommhold, L. 2001, ApJ 546, 1168 
! 
      use accura 
      implicit real(dp)(A-H,O-Z) 
      integer, parameter :: nlines=43 
      real(dp), save :: freq(nlines),alpha(nlines,11) 
      real(dp)       :: temp(11) 
      real(dp), parameter :: amagat=2.6867774e+19,fac=1./amagat**2 
      data temp / 1000.,  1500.,  2250., 3000.,  4000.,  5000.,           & 
     &            6000.,  7000., 8000.,  9000., 10000./ 
      data ntemp /11/ 
      data ifirst /0/ 
      REAL(DP), PARAMETER :: CAS=2.997925e10 
 
!     input frequency in Hz but needed wave numbers in cm^-1 
      f=ff/cas 
!     read in CIA tables if this is the first call 
      if (ifirst.eq.0) then 
         write(*,'(a)') 'Reading in H-He CIA opacity tables...' 
         open(10,file="./data/CIA_HHe.dat",status='old') 
         do i=1,3 
            read (10,*) 
         end do 
         do i=1,nlines 
            read (10,*) freq(i),(alpha(i,j),j=1,ntemp) 
         end do 
         close(10) 
 
!     take logarithm of tables prior to doing linear interpolations 
 
         do i=1,nlines 
            do j=1,ntemp 
               alpha(i,j)=log(alpha(i,j)) 
            enddo 
         end do 
 
         ifirst=1 
      end if 
 
!     locate position in temperature array 
      call locate(temp,ntemp,t,j,ntemp) 
 
      if(j.eq.0) then 
         write(66,'(a,f6.0,a)')                                           & 
     &   'Warning: requested temperature is below',temp(1),' K' 
         write(66,'(a)') 'CIA H-He opacity is extrapolated down' 
!        opac=0. 
!        return 
         j=1 
      end if 
 
!     locate position in frequency array 
      call locate(freq,nlines,f,i,nlines) 
 
!     linearly interpolate in frequency and temperature 
 
      if(j.eq.ntemp) then 
!        hold values constant if off high temperature end of table 
         y1=alpha(i,j) 
         y2=alpha(i+1,j) 
         tt=(f-freq(i))/(freq(i+1)-freq(i)) 
         alp=(1.-tt)*y1 + tt*y2 
       else if (i.eq.0 .or. i.eq.nlines) then 
!        set values to a very small number if off frequency table 
         alp=-50. 
       else 
!        interpolate linearly within table 
         y1=alpha(i,j) 
         y2=alpha(i+1,j) 
         y3=alpha(i+1,j+1) 
         y4=alpha(i,j+1) 
 
         tt=(f-freq(i))/(freq(i+1)-freq(i)) 
         uu=(t-temp(j))/(temp(j+1)-temp(j)) 
 
         alp=(1.-tt)*(1.-uu)*y1 + tt*(1.-uu)*y2 + tt*uu*y3 +              & 
     &       (1.-tt)*uu*y4 
      end if 
 
      alp=exp(alp) 
 
!     final opacity 
 
      opac=fac*ah*ahe*alp 
! 
      return 
      end subroutine cia_hhe 
! 
! 
! ******************************************************************* 
! 
! 
      subroutine h2minus(t,anh2,ane,fr,oph2m) 
!     ======================================= 
! 
!     H2- free-free opacity 
! 
!     data from K L Bell 1980 J. Phys. B: At. Mol. Phys. 13 1859, Table 1 
!     The first column is theta=5040/T(K) 
!     The first row are names for each row corresponding to lambda (angstroms) 
!     The last row for 10.0 is linearly extrapolated 
!     The units of everything else is 10^26 cm4/dyn-1 
! 
      use accura 
      use params 
      implicit real(dp) (a-h,o-z) 
 
      real(dp) :: FFthet(9),FFlamb(18),FFkapp(18,9) 
      data FFthet / 0.5, 0.8, 1.0, 1.2, 1.6, 2.0,                         & 
     &     2.8, 3.6, 10.0 / 
      data nthet /9/ 
      data FFlamb /151883., 113913.,  91130.,  60753.,                    & 
     &     45565.,  36452.,  30377.,  22783.,                             & 
     &     18226.,  15188.,  11391.,  9113.,  7594.,                      & 
     &     6509.,  5696.,  5063.,  4142.,  3505./ 
      data nlamb /18/ 
      data FFkapp /                                                       & 
     &     7.16e+01,4.03e+01,2.58e+01,1.15e+01,6.47e+00,                  & 
     &     4.15e+00,2.89e+00,1.63e+00,1.05e+00,7.36e-01,                  & 
     &     4.20e-01,2.73e-01,1.92e-01,1.43e-01,1.10e-01,                  & 
     &     8.70e-02,5.84e-02,4.17e-02,9.23e+01,5.20e+01,                  & 
     &     3.33e+01,1.48e+01,8.37e+00,5.38e+00,3.76e+00,                  & 
     &     2.14e+00,1.39e+00,9.75e-01,5.64e-01,3.71e-01,                  & 
     &     2.64e-01,1.98e-01,1.54e-01,1.24e-01,8.43e-02,                  & 
     &     6.10e-02,1.01e+02,5.70e+01,3.65e+01,1.63e+01,                  & 
     &     9.20e+00,5.92e+00,4.14e+00,2.36e+00,1.54e+00,                  & 
     &     1.09e+00,6.35e-01,4.22e-01,3.03e-01,2.30e-01,                  & 
     &     1.80e-01,1.46e-01,1.01e-01,7.34e-02,1.08e+02,                  & 
     &     6.08e+01,3.90e+01,1.74e+01,9.84e+00,6.35e+00,                  & 
     &     4.44e+00,2.55e+00,1.66e+00,1.18e+00,6.97e-01,                  & 
     &     4.67e-01,3.39e-01,2.59e-01,2.06e-01,1.67e-01,                  & 
     &     1.17e-01,8.59e-02,1.18e+02,6.65e+01,4.27e+01,                  & 
     &     1.91e+01,1.08e+01,6.99e+00,4.91e+00,2.84e+00,                  & 
     &     1.87e+00,1.34e+00,8.06e-01,5.52e-01,4.08e-01,                  & 
     &     3.17e-01,2.55e-01,2.10e-01,1.49e-01,1.11e-01,                  & 
     &     1.26e+02,7.08e+01,4.54e+01,2.04e+01,1.16e+01,                  & 
     &     7.50e+00,5.28e+00,3.07e+00,2.04e+00,1.48e+00,                  & 
     &     9.09e-01,6.33e-01,4.76e-01,3.75e-01,3.05e-01,                  & 
     &     2.53e-01,1.82e-01,1.37e-01,1.38e+02,7.76e+01,                  & 
     &     4.98e+01,2.24e+01,1.28e+01,8.32e+00,5.90e+00,                  & 
     &     3.49e+00,2.36e+00,1.74e+00,1.11e+00,7.97e-01,                  & 
     &     6.13e-01,4.92e-01,4.06e-01,3.39e-01,2.49e-01,                  & 
     &     1.87e-01,1.47e+02,8.30e+01,5.33e+01,2.40e+01,                  & 
     &     1.38e+01,9.02e+00,6.44e+00,3.90e+00,2.68e+00,                  & 
     &     2.01e+00,1.32e+00,9.63e-01,7.51e-01,6.09e-01,                  & 
     &     5.07e-01,4.27e-01,3.16e-01,2.40e-01,2.19e+02,                  & 
     &     1.26e+02,8.13e+01,3.68e+01,2.18e+01,1.46e+01,                  & 
     &     1.08e+01,7.18e+00,5.24e+00,4.17e+00,3.00e+00,                  & 
     &     2.29e+00,1.86e+00,1.55e+00,1.32e+00,1.13e+00,                  & 
     &     8.52e-01,6.64e-01/ 
 
!     locate position in temperature array 
      theta=5040./t 
 
      call locate(FFthet,nthet,theta,j,nthet) 
      if(j.eq.0) then 
         write(66,*) 
         write(66,'(a,f6.0,a)')                                           & 
     &   'warning: requested temperature is outside the ranges' 
         write(66,'(a)') 'h2minus cross-section is set to 0' 
         write(66,*) 
!        j=1 
!        stop 
         oph2m=0. 
         return 
      end if 
      flamb=CL*1.e8/fr 
!     locate position in wavelength array 
      call locate(FFlamb,nlamb,flamb,i,nlamb) 
 
!     linearly interpolate in frequency and temperature 
      if(j.eq.nthet) then 
!        hold values constant if off high temperature end of table 
         y1=FFkapp(i,j) 
         y2=FFkapp(i+1,j) 
         tt=(flamb-FFlamb(i))/(FFlamb(i+1)-FFlamb(i)) 
         Fkappa=(1.-tt)*y1 + tt*y2 
       else if (i.eq.0 .or. i.eq.nthet) then 
!        set values to 0 if off frequency table 
         Fkappa=0.0 
       else 
!        interpolate linearly within table 
         y1=FFkapp(i,j) 
         y2=FFkapp(i+1,j) 
         y3=FFkapp(i+1,j+1) 
         y4=FFkapp(i,j+1) 
 
         tt=(flamb-FFlamb(i))/(FFlamb(i+1)-FFlamb(i)) 
         uu=(theta-FFthet(j))/(FFthet(j+1)-FFthet(j)) 
 
         Fkappa=(1.-tt)*(1.-uu)*y1 + tt*(1.-uu)*y2 + tt*uu*y3 +           & 
     &       (1.-tt)*uu*y4 
      end if 
      pe=ane*BOLK*t 
      oph2m= anh2 * 1.0E-26 *pe * Fkappa 
      return 
      end subroutine h2minus 
! 
! 
!   ********************************************************************** 
! 
! 
      subroutine h2opf(t,pf) 
! 
!     partition function for H2O from EXOMOIL data 
! 
      use accura 
      use params 
      use modelp,only : ttab,pftab 
 
      implicit real(dp) (a-h,o-z) 
! 
      data init /1/ 
! 
      if(init.eq.1) then 
         open(67,file='./data/h2o_exomol.pf',status='old') 
         do i=1,10000 
            read(67,*) ttab(i),pftab(i) 
         end do 
         close(67) 
         init=0 
      end if 
! 
      itab=ifix(real(t)) 
      pf=pftab(itab)+(t-ttab(itab))*(pftab(itab+1)-pftab(itab)) 
      return 
      end subroutine h2opf 
 
! 
! 
!   ********************************************************************** 
! 
! 
      subroutine vopf(t,pf) 
! 
!     partition function for VO from EXOMOILA data 
! 
      use accura 
      use params 
      use modelp, only : ttab,pftab 
 
      implicit real(dp) (a-h,o-z) 
! 
      data init /1/ 
! 
      if(init.eq.1) then 
         open(67,file='./data/vo_exomol.pf',status='old') 
         do i=1,8000 
            read(67,*) ttab(i),pftab(i) 
         end do 
         close(67) 
         init=0 
      end if 
! 
      itab=ifix(real(t)) 
      pf=pftab(itab)+(t-ttab(itab))*(pftab(itab+1)-pftab(itab)) 
      return 
      end subroutine vopf 
 
! 
 
! 
! 
! ******************************************************************* 
! 
! 
      function gvdw(il,ilist,id) 
!     ========================== 
! 
!     evaluation of the Van der Waals broadening parameter 
! 
!     currently, two possibilities, determined by the value of the parameter 
!     ivdwli(ilist) - the mode of evaluation is the same for the whole line list
!       = 0 - standard expression 
!       > 0 - evaluation using EXOMOL data, assuming breadening by H2 and He 
! 
      use accura 
      use params 
      use modelp 
      use lindat 
      use molist 
 
      implicit real(dp) (a-h,o-z) 
! 
!     clasical, original expression 
! 
!C    if(ivdwli(ilist).eq.0) then 
!C       gvdw=gwm(il,ilist)*vdwc(id) 
!C       return 
!C    end if 
! 
!     if(ibroli(ilist).eq.0.and ivdwli(ilist).eq.0) then 
!        gvdw=gwstd*vdwc(id) 
!        return 
!     end if 
! 
!     EXOMOL form - broadening by H2 and He 
! 
!     con= 1.e-6*c*k 
      con=4.1388e-12 
      t=temp(id) 
      anhe=rrr(id,1,2) 
      gvdw=con*t*((296./t)**gexph2(il,ilist)*gvdwh2(il,ilist)*anh2(id)+   & 
     &            (296./t)**gexphe(il,ilist)*gvdwhe(il,ilist)*anhe) 
      return 
      end function gvdw 
! 
! 
! ******************************************************************* 
! 
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
! 
! ******************************************************************* 
! 
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
!       write(*,*) 'uato',jatom,ion,tl,ulog 
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
! ========================================================================= 
 
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
!     Output: pfeq - either the partition function of equil.const. 
! 
!     array  IBC is the BC indec corresponfing to the Tsuji-like index Iindtsu 
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
      if(in.eq.0) return 
! 
!     temperature interpolation in the BC tables 
! 
      iii=int(t)/10 
      indt=min(iii,999) 
      a1=1.- (real(log10(t))-tlbc(indt))/(tlbc(indt+1)-tlbc(indt)) 
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
 
! ============================================================ 
! 
 
! 
! 
! *********************************************************************** 
! 
! 
      subroutine frac_fd(iteos) 
!     ========================= 
! 
      use accura 
      use params 
      use modelp 
      use opadat 
      implicit real(dp) (a-h,o-z) 
 
      integer  :: jdatt(30) 
      real(dp) :: frad(31,30),zc(31,30),felec(30) 
      data jdatt   / 1, 2, 0, 0, 0, 3, 4, 5, 0, 6,                        & 
     &              7, 8, 9,10, 0,11, 0,12,13,14,                         & 
     &              0, 0, 0, 0, 0,15, 0, 0, 0, 0/ 
 
      call alloc_opadat 
 
      if(iteos.le.1.and.igrdd.eq.0) call read_fd 
 
! 
!     zero partition functions 
! 
      do id=1,nd 
         do j=1,30 
            do i=1,31 
               pfunc(i,j,id)=0. 
            end do 
         end do 
      end do 
! 
      do id=1,nd 
         t=temp(id) 
         ane=elec(id) 
         rho=dens(id) 
 
         call ion_fd(t,ane,frad,zc,felec) 
 
         do iat=1,30 
            if(jdatt(iat).gt.0) then 
               do ion=1,iat+1 
                  if(zc(ion,iat).gt.0.)                                   & 
     &            rrr(id,ion,iat)=frad(ion,iat)/zc(ion,iat)*              & 
     &                            abndd(iat,id)*                          & 
     &                            dens(id)/wmm(id)/ytot(id) 
                  pfunc(ion,iat,id)=zc(ion,iat) 
                  frpf0(ion,iat)=frad(ion,iat) 
               end do 
            end if 
         end do 
      end do 
! 
      return 
      end subroutine frac_fd 
 
 
 
      subroutine ion_fd(T,ane,frad,zc,felec) 
!     ====================================== 
! 
!     calculates partition functions, and possibly 
!     ionization fractions by interpolation from 
!     tabulated values provided by Franck Delahaye 
!     (calculations similar to the Opacity Project) 
! 
!     routine written by C. Stehle, modified  by I.H. on May 2011 
! 
      use accura 
      use opadat 
      implicit real(dp) (a-h,o-z) 
 
! Tabulated ionic fractions fract, partition functions Zct 
! mean ion charge felect  for the grid values of MHD data 
!  97 is the number of temperatures, 31 or 29 the number of rho values 
 
      real(dp) :: frad(31,30),zc(31,30),felec(30) 
      real(dp) :: zc_Ne_it0(31,30),zc_Ne_it1(31,30) 
      real(dp) :: frac_Ne_it0(31,30),frac_Ne_it1(31,30) 
 
      integer  :: jdatt(30) 
      data jdatt   / 1, 2, 0, 0, 0, 3, 4, 5, 0, 6,                        & 
     &              7, 8, 9,10, 0,11, 0,12,13,14,                         & 
     &              0, 0, 0, 0, 0,15, 0, 0, 0, 0/ 
 
 
      altref=log10(T) 
      alrref=log10(ane) 
 
! ........ bracketing in T trough the grid (data_ion_read) 
!          ( -> it0,it1) 
 
      if(altref.le.alogtt(1)) altref=alogtt(1) 
      itloop: do it=1,96 
         testT=(altref-alogtt(it))*(altref-alogtt(it+1)) 
         if(testT.le.0.)  then 
           it0=it 
           it1=it+1 
           exit itloop 
         endif 
      end do itloop 
 
      dtemp = alogtt(it1)-alogtt(it0) 
      prop_T= (altref-alogtt(it0))/dtemp 
 
 
! ........ end of bracketing in T trough the grid (data_ion_read) 
!          ( -> it0,it1) 
 
! ....... LOOP on atomic species IAT 
 
       do iat=1,30 
         if (jdatt(iat).gt.0) then 
 
! ........ bracketing in Ne first for it0 
 
          irhloop: do ir=1,nrho(iat)-1 
             if(alrref.le.alognet(it0,1,iat) )                            & 
     &       alrref = alognet(it0,1,iat) 
             TestR=(alrref-alognet(it0,ir,iat))*                          & 
     &             (alrref-alognet(it0,ir+1,iat)) 
             if (TestR.le.0.) then 
                ir0=ir 
                ir1=ir+1 
                exit irhloop 
             end if 
          end do irhloop 
 
! ........ end of bracketing in Ne first for it0 
 
 
! .....  interpolation in Ne of ionic fractions for it0 
 
          Prop_Ne=(alrref-alognet(it0,ir0,iat))/                          & 
     &          (alognet(it0,ir1,iat)-alognet(it0,ir0,iat)) 
 
          do ion = 1, iat+1 
            frac_Ne_it0(ion,iat)= fract(it0,ir0,ion,iat)+                 & 
     &      Prop_Ne*(fract(it0,ir1,ion,iat)-fract(it0,ir0,ion,iat)) 
            zc_Ne_it0(ion,iat)=zct(it0,ir0,ion,iat)+                      & 
     &      Prop_Ne*(zct(it0,ir1,ion,iat)-zct(it0,ir0,ion,iat)) 
         end do 
 
! .....  end of interpolation in Ne of ionic fractions for it0 
 
 
! ........ bracketing in Ne  for it1 
 
          irloop: do ir=1,nrho(iat)-1 
             if(alrref.le.alognet(it1,1,iat) )                            & 
     &         alrref = alognet(it1,1,iat) 
               TestR=(alrref-alognet(it1,ir,iat))*                        & 
     &           (alrref-alognet(it1,ir+1,iat)) 
               if (TestR.le.0.) then 
                 ir0=ir 
                 ir1=ir+1 
                 exit irloop 
             end if 
          end do irloop 
 
! ........ end of  bracketing in Ne  for it1 
 
! .....  interpolation in Ne of ionic fractions for it1 
 
          Prop_Ne=(alrref-alognet(it1,ir0,iat))/                          & 
     &          (alognet(it1,ir1,iat)-alognet(it1,ir0,iat)) 
 
!            if(Prop_Ne .lt.0.) then 
!            write(6,*) ir0,ir1,Prop_Ne,alrref,alognet(it1,ir0,iat), 
!     *      alognet(it1,ir1,iat) 
!            endif 
 
            if(Prop_Ne .lt.0.)                                            & 
     &      write(6,"('N-e interpolation with negative prop_ne, iat= ',   & 
     &      i3)") iat 
 
        do ion = 1, iat+1 
            frac_Ne_it1(ion,iat)= fract(it1,ir0,ion,iat)+                 & 
     &      Prop_Ne*(fract(it1,ir1,ion,iat)-fract(it1,ir0,ion,iat)) 
            zc_Ne_it1(ion,iat)=zct(it1,ir0,ion,iat)+                      & 
     &      Prop_Ne*(zct(it1,ir1,ion,iat)-zct(it1,ir0,ion,iat)) 
        end do 
 
! .....  end of interpolation in Ne of ionic fractions for it1 
 
!  :::::::::: FINAL INTERPOLATION IN T 
!            returns ionic fractions frac, partition function zc 
!            and mean ion charge felec 
 
         felec (iat) =0. 
 
        do ion=1,iat+1 
          frad(ion,iat) =   frac_Ne_it0(ion,iat) +                        & 
     &    prop_T * (frac_Ne_it1(ion,iat)-frac_Ne_it0(ion,iat)) 
          zc(ion,iat) =    zc_Ne_it0(ion,iat)+                            & 
     &     prop_T * (zc_Ne_it1(ion,iat) - zc_Ne_it0(ion,iat)) 
 
 
          charge=ion-1. 
          felec(iat)=  felec(iat)+charge*frad(ion,iat) 
        end do 
        end if 
      end do 
 
! .......end of LOOP on atomic species IAT 
 
      return 
      end subroutine ion_fd 
                   	 
! ****************************************************************** 
 
      subroutine read_fd 
!     ================== 
! 
!     reading Franck Delahaye data for ionization fractions and 
!     partition functions 
! 
!     data tabulated for 97 temperatures and 31 or 29 densities 
! 
!     routine written by C. Stehle, modified  by I.H. on May 2011 
! 
      use accura 
      use opadat 
      implicit real(dp) (a-h,o-z) 
! 
      real(dp) :: am(30) 
      integer  :: jdatt(30),nrh(30),nrt(30) 
      character(len=9) ::  filnam(15) 
! 
!     jdatt - an internal index of the species in the tables 
! 
      data jdatt   / 1, 2, 0, 0, 0, 3, 4, 5, 0, 6,                        & 
     &              7, 8, 9,10, 0,11, 0,12,13,14,                         & 
     &              0, 0, 0, 0, 0,15, 0, 0, 0, 0/ 
       data filnam  /'H_cs.dat','He_fd.dat','C_fd.dat',                   & 
     &               'N_fd.dat','O_fd.dat','Ne_fd.dat',                   & 
     &               'Na_fd.dat','Mg_fd.dat','Al_fd.dat',                 & 
     &               'Si_fd.dat','S_fd.dat','Ar_fd.dat','K_fd.dat',       & 
     &               'Ca_fd.dat','Fe_fd.dat'/ 
      data nrh / 28, 30,  0,  0,  0, 28, 28, 28,  0, 31,                  & 
     &           29, 31, 30, 31,  0, 31,  0, 31, 31, 31,                  & 
     &            0,  0,  0,  0,  0, 31,  0,  0,  0,  0/ 
      data nrt / 97, 97,  0,  0,  0, 97, 97, 97,  0, 97,                  & 
     &           97, 97, 97, 97,  0, 97,  0, 97, 97, 97,                  & 
     &            0,  0,  0,  0,  0, 97,  0,  0,  0,  0/ 
      data am/1.008,  4.003,  6.941,  9.012, 10.810, 12.011,              & 
     &       14.007, 16.000, 18.918, 20.179, 22.990, 24.305,              & 
     &       26.982, 28.086, 30.974, 32.060, 35.453, 39.948,              & 
     &       39.098, 40.080, 44.956, 47.900, 50.941, 51.996,              & 
     &       54.938, 55.847, 58.933, 58.700, 63.546, 65.380/ 
! 
! 
!       loop over species 
! 
 
                    j=jdatt(1) 
 
          nrho(1)=nrh(1) 
          open(66, file=filnam(1),status='unknown') 
 
!        hydrogen 
 
         iat=1 
          iMolecule=1 
          tloop: do it=1,nrt(1) ! loop over T 
            rloop: do ie=1,nrho(1)  ! loop over rho 
               read(66,*,iostat=ios) 
               if(ios.ne.0) exit tloop 
               read(66,"(10X,F13.8,10X,F13.8)",iostat=ios)                & 
     &         alogtt(it),alogrhot(ie,iat) 
               if(ios.ne.0) exit tloop 
               read(66,*) 
               read(66,*,iostat=ios) 
               if(ios.ne.0) exit tloop 
               read(66,"(4(9X,1pe11.4))")                                 & 
     &          fractH2(it,ie),fractH2p(it,ie),fracHm(it,ie),             & 
     &          fract(it,ie,1,iat) 
               read(66,"(4(9X,1pe11.4))")  fract(it,ie,2,iat) 
               read(66,"(7X,1pe19.12,7x,1pe19.12,6x,1pe19.12)")           & 
     &          zctH2(it,ie),zctH2p(it,ie), zctHm(it,ie) 
               read(66,"(8X,1pe19.12,8x,1pe19.12)")                       & 
     &          zct(it,ie,1,iat),zct(it,ie,2,iat) 
               read(66,*) 
               read(66,"(9X,1pe11.4)") felect(it,ie,iat) 
! 
               if(iMolecule.eq.0 ) then 
                   fract(it,ie,1,iat) = fract(it,ie,1,iat)+               & 
     &             2.*fractH2(it,ie) 
                   fract(it,ie,2,iat)= fract(it,ie,2,iat) +               & 
     &             fractH2p(it,ie)-fracHm(it,ie) 
               end if 
 
             end do rloop 
          end do tloop 
        close(66) 
 
        atoms: do iat=2,30 
          j=jdatt(iat) 
          if (j.gt.0) then 
          nrho(iat)=nrh(iat) 
          ntempe(iat)= nrt(iat) 
          write(6,*) 
          open(66, file=filnam(j),status='unknown') 
          write(6,*) iat,filnam(j), nrt(iat),nrho(iat) 
! 
          tloop2: do it=1,nrt(iat) ! loop over T 
            do ie=1,nrho(iat)  ! loop over rho 
               read(66,*,iostat=ios) 
               if(ios.ne.0) exit tloop2 
               read(66,"(10X,F13.8,10X,F13.8)",iostat=ios)                & 
     &          alogtt(it),alogrhot(ie,iat) 
               if(ios.ne.0) exit tloop2 
               read(66,*) 
               read(66,*,iostat=ios) 
               if(ios.ne.0) exit tloop2 
               read(66,"(4(9X,1pe11.4))")                                 & 
     &          (fract(it,ie,ion,iat),ion=1,iat+1) 
               read(66,"(3(9X,1pe19.12))")                                & 
     &          (zct(it,ie,ion,iat),ion=1,iat+1) 
               read(66,*) 
               read(66,"(9X,1pe11.4)") felect(it,ie,iat) 
            end do 
          end do tloop2 
 
          close(66) 
          end if 
 
        end do atoms 
 
!       end of loop over species 
! 
!      electron density for the mass density values from the 
!      tables (different for each species) 
 
       amh= 1.673e-24  ! mass hydrogen atom in g 
 
       do iat=1,30 
          if(jdatt(iat).gt.0) then 
             do it=1,nrt(iat) 
               do ir=1,nrho(iat) 
                  rhoir=10.**alogrhot(ir,iat) 
                  dens_ion= felect(it,ir,iat)*rhoir/(am(iat)*amh) 
                  alognet(it,ir,iat)=log10(dens_ion) 
               end do 
             end do 
          end if 
       end do 
 
       return 
       end subroutine read_fd 
! 
! 
!    ************************************************************************** 
! 
! 
 
! 
!     ********************************************************** 
! 
 
      subroutine dkini 
!     ================ 
! 
!     Initializes necessary arrates for evaluating hydrogen Lyman alpha 
!     and beta line profiles from the tables adopted from Detlev Koester's 
!     program ATMDK 
! 
      use accura 
      implicit real(dp) (a-h,o-z) 
      character(len=45) :: dkfile(5) 
      data dkfile/'/Users/ivan/atmdk/data/ly_alpha_hh.tab',               & 
     &            '/Users/ivan/atmdk/data/ly_alpha_hp2026.tab',           & 
     &            '/Users/ivan/atmdk/data/ly_alpha_h2.tab',               & 
     &            '/Users/ivan/atmdk/data/ly_beta_hh2026.tab',            & 
     &            '/Users/ivan/atmdk/data/ly_beta_hp2026.tab'/ 
! 
      ntab=5 
      write(*,*) 
      do itb=1,ntab 
         call initabdk(itb,dkfile(itb)) 
!        write(*,*) 'reading Koester table:',dkfile(itab) 
      end do 
      return 
      end subroutine dkini 
 
 
! 
!     ********************************************************** 
! 
 
      subroutine initabdk(itab,tabfile) 
!     ================================= 
! 
!     initoializes one particular Koester's table for hydrogen 
 
      use accura 
      use hydprf 
      implicit real(dp) (a-h,o-z) 
 
      character(len=72) :: inpline 
      character(len=40) :: tabfile 
! 
      iun=61 
      open(iun,file=tabfile) 
      write(*,*) 'reading DK table ',tabfile 
      do i=1,100 
         read(iun,'(a)') inpline 
         if (inpline(1:3).eq.'END'.or.inpline(5:7).eq.'---') exit 
      end do 
 
      read(iun,*) nden(itab) 
!        write(*,*) 'itab,nden(itab)',itab,nden(itab) 
      do iden=1,nden(itab) 
         read(iun,*) tabden(itab,iden),numlam(itab,iden),                 & 
     &               numtem(itab,iden),                                   & 
     &               (tabtem(itab,item,iden),item=1,numtem(itab,iden)) 
!        write(*,*) 'iden,tabden(itab,iden)',iden,tabden(itab,iden) 
!        write(*,*) '   numlam',numlam(itab,iden) 
!        write(*,*) 'numtem',numtem(itab,iden) 
        do ilam=1,numlam(itab,iden) 
            read(iun,*) tablam(itab,ilam,iden),                           & 
     &      (prftab(itab,ilam,item,iden),item=1,numtem(itab,iden)) 
        end do 
!        write(*,*) 'numlam',numlam(itab,iden) 
      end do 
      close(iun) 
      return 
      end subroutine initabdk 
! 
!     ********************************************************** 
! 
      subroutine dkprof(itab,wl0,tt0,dd0,prof) 
!     ======================================== 
! 
      use accura 
      use hydprf 
      implicit real(dp) (a-h,o-z) 
 
      real(dp) :: tabdei(mden),wl1(mlam),tl1(mtem) 
! 
      denlog=log10(dd0) 
      ndei=nden(itab) 
      do iden=1,ndei 
         tabdei(iden)=tabden(itab,iden) 
      end do 
      call locat(tabdei,ndei,denlog,jd,fd) 
      jd1=jd 
      jd2=jd1+1 
      fd2=fd 
      fd1=1.-fd2 
! 
!     lower density point 
! 
      prof=1.e-40 
      do ilam=1,numlam(itab,jd1) 
         wl1(ilam)=tablam(itab,ilam,jd1) 
      end do 
      if(wl0.lt.wl1(1).or.wl0.gt.wl1(numlam(itab,jd1))) return 
      call locat(wl1,numlam(itab,jd1),wl0,jl,fl) 
      jl1=jl 
      jl2=jl1+1 
      fl2=fl 
      fl1=1.-fl2 
!      write(*,*) 'jl1,jl2,fl2',jl1,jl2,fl2 
!      write(*,*) 'wl1(jl1),wl1(jl2)',wl1(jl1),wl1(jl2),wl0 
      do item=1,numtem(itab,jd1) 
         tl1(item)=tabtem(itab,item,jd1) 
      end do 
      call locat(tl1,numtem(itab,jd1),tt0,jt,ft) 
      jt1=jt 
      jt2=jt1+1 
!      write(*,*) 'jt1,jt2,ft',jt1,jt2,ft 
      ft2=ft 
      ft1=1.-ft2 
      pr11=fl1*prftab(itab,jl1,jt1,jd1)+fl2*prftab(itab,jl2,jt1,jd1) 
      pr12=fl1*prftab(itab,jl1,jt2,jd1)+fl2*prftab(itab,jl2,jt2,jd1) 
      pr1=ft1*pr11+ft2*pr12 
! 
!     upper density point 
! 
      do ilam=1,numlam(itab,jd2) 
         wl1(ilam)=tablam(itab,ilam,jd2) 
      end do 
      call locat(wl1,numlam(itab,jd2),wl0,jl,fl) 
      jl1=jl 
      jl2=jl1+1 
      fl2=fl 
      fl1=1.-fl2 
      do item=1,numtem(itab,jd2) 
         tl1(item)=tabtem(itab,item,jd2) 
      end do 
      call locat(tl1,numtem(itab,jd2),tt0,jt,ft) 
      jt1=jt 
      jt2=jt1+1 
      ft2=ft 
      ft1=1.-ft2 
      pr21=fl1*prftab(itab,jl1,jt1,jd2)+fl2*prftab(itab,jl2,jt1,jd2) 
      pr22=fl1*prftab(itab,jl1,jt2,jd2)+fl2*prftab(itab,jl2,jt2,jd2) 
      pr2=ft1*pr21+ft2*pr22 
! 
      prof=fd1*pr1+fd2*pr2+denlog-18. 
      prof=10.**prof 
      return 
      end subroutine dkprof 
! 
!     ***************************************************************** 
! 
 
      SUBROUTINE locat(xx,n,x,j,a) 
!     ============================ 
! 
      use accura 
      implicit real(dp) (a-h,o-z) 
      real(dp) :: xx(n) 
! 
      jl=0 
      ju=n+1 
      loc: do 
         if(ju-jl.gt.1) then 
           jm=(ju+jl)/2 
           if((xx(n).gt.xx(1)).eqv.(x.ge.xx(jm)))then 
              jl=jm 
           else 
              ju=jm 
           endif 
           cycle loc 
          else 
           exit loc 
         end if 
      end do loc 
      if(x.le.xx(1)) then 
        j=1 
        a=0. 
       else if(x.ge.xx(n)) then 
        j=n-1 
        a=1. 
       else 
        j=jl 
        a=(x-xx(j))/(xx(j+1)-xx(j)) 
      end if 
      return 
      END SUBROUTINE LOCAT 
 
! 
!     ****************************************************** 
! 
 
      subroutine lalpdk(wl,t,ane,anp,anh,anh2m,prof) 
!     ============================================== 
! 
!     Lyman alpha profile using DK tables 
! 
      use accura 
      implicit real(dp) (a-h,o-z) 
      logical :: lpr 
! 
      lpr=abs(t-13200).lt.100. 
      lpr=.false. 
!     broadening by neutral hydrogen 
      call dkprof(1,wl,t,anh,profhh) 
!     broadening by protons 
      call dkprof(2,wl,t,anp,profhp) 
!     broadening by H2 
      if(anh2m.gt.1.e-10*anh) call dkprof(3,wl,t,anh2m,profh2) 
!     broadening by electrons (after DK) 
      dl=abs(wl-1215.67) 
      zz=7862.67*dl/t 
      zz=0.5+0.5*(zz-1.)/(zz+1.) 
      ff=1.254e-9*ane**0.666667 
      alp=log10(max(dl/ff,1.04194e-3)) 
      profe=1.34896/3.*exp(2.3025851*3.0039-5.7565*alp)/ff 
! 
!     total profile 
! 
      prof=profhh+profhp+profh2+profe 
      prof=prof*1215.67**2*88.53*0.4162*1.e-30 
      if(lpr) write(6,"(f10.3,f8.1,1p3e9.1,2x,5e10.2)")                   & 
     &        wl,t,ane,anh,anp,profhh,profhp,profh2,profe,prof 
      return 
      end subroutine lalpdk 
 
! 
!     ****************************************************** 
! 
 
      subroutine lbetdk(wl,t,ane,anp,anh,prof) 
!     ======================================== 
! 
!     Lyman beta profile using DK tables 
! 
      use accura 
      implicit real(dp) (a-h,o-z) 
      logical ::  lpr 
! 
      lpr=abs(t-13200).lt.100. 
      lpr=.false. 
!     broadening by neutral hydrogen 
      call dkprof(4,wl,t,anh,profhh) 
!     broadening by protons 
      call dkprof(5,wl,t,anp,profhp) 
!     broadening by electrons (after DK) 
      dl=abs(wl-1025.73) 
      ff=1.254e-9*ane**0.666667 
      alp=log10(max(dl/ff,1.5258e-3)) 
      profe=exp(2.3025851*3.25358-5.7565*alp)/ff 
! 
!     total profile 
! 
      prof=profhh+profhp+profe 
      prof=prof*1025.73**2*88.53*0.0791*1.e-30 
      if(lpr) write(6,"(f10.3,f8.1,1p3e9.1,2x,5e10.2)")                   & 
     &        wl,t,ane,anh,anp,profhh,profhp,profe,prof 
      return 
      end subroutine lbetdk 
 
! 
!     ****************************************************** 
! 
 
      subroutine molindx 
!     ================== 
! 
!     a systematic translation of the Kurucz labels to the Tsuji indices 
!     set by Carlos Allende Prieto, Jan 2026 
! 
      use accura 
      use molist 
 
      molind=0 
 
      molind(199)=1 
      molind(101)=2 
      molind(10108)=3 
      molind(108)=4 
      molind(106)=5 
      molind(608)=6 
      molind(607)=7 
      molind(606)=8 
      molind(707)=9 
      molind(808)=10 
      molind(708)=11 
      molind(107)=12 
      molind(10607)=14 
      molind(10606)=15 
      molind(116)=16 
      molind(114)=17 
      molind(616)=20 
      molind(614)=21 
      molind(716)=23 
      molind(714)=24 
      molind(814)=25 
      molind(816)=26 
      molind(1616)=27 
      molind(1416)=28 
      molind(822)=29 
      molind(823)=30 
      molind(840)=31 
      molind(112)=32 
      molind(109)=33 
      molind(120)=34 
      molind(857)=35 
      molind(117)=36 
      molind(10106)=38 
      molind(10107)=40 
      molind(999)=45 
      molind(913)=46 
      molind(920)=47 
      molind(10812)=49 
      molind(10813)=51 
      molind(10811)=54 
      molind(10116)=57 
      molind(1317)=59 
      molind(1117)=60 
      molind(1719)=61 
      molind(10819)=62 
      molind(1720)=63 
      molind(10820)=65 
      molind(1622)=67 
      molind(103)=69 
      molind(308)=70 
      molind(309)=71 
      molind(317)=72 
      molind(10308)=73 
      molind(10104)=75 
      molind(408)=76 
      molind(409)=77 
      molind(417)=78 
      molind(10408)=80 
      molind(105)=82 
      molind(10105)=83 
      molind(508)=84 
      molind(516)=86 
      molind(509)=87 
      molind(517)=88 
      molind(10508)=89 
      molind(699)=92 
      molind(10699)=94 
      molind(10708)=109 
      molind(899)=113 
      molind(10899)=115 
      molind(10608)=118 
      molind(909)=120 
      molind(809)=121 
      molind(111)=122 
      molind(811)=123 
      molind(911)=124 
      molind(812)=126 
      molind(1216)=127 
      molind(912)=128 
      molind(1217)=130 
      molind(113)=133 
      molind(813)=134 
      molind(1316)=136 
      molind(1499)=139 
      molind(10114)=140 
      molind(914)=144 
      molind(1417)=146 
      molind(115)=148 
      molind(10115)=149 
      molind(615)=151 
      molind(715)=152 
      molind(815)=153 
      molind(1516)=155 
      molind(915)=156 
      molind(1517)=158 
      molind(1799)=163 
      molind(1717)=164 
      molind(617)=165 
      molind(817)=169 
      molind(1617)=172 
      molind(10817)=174 
      molind(119)=176 
      molind(819)=177 
      molind(919)=178 
      molind(820)=179 
      molind(1620)=180 
      molind(2121)=182 
      molind(721)=183 
      molind(821)=184 
      molind(1621)=187 
      molind(921)=188 
      molind(922)=190 
      molind(1722)=192 
      molind(723)=195 
      molind(124)=198 
      molind(724)=199 
      molind(824)=200 
      molind(924)=202 
      molind(1724)=204 
      molind(125)=206 
      molind(825)=207 
      molind(1625)=209 
      molind(925)=211 
      molind(1725)=213 
      molind(126)=214 
      molind(826)=216 
      molind(1626)=217 
      molind(926)=218 
      molind(1726)=220 
      molind(128)=223 
      molind(828)=224 
      molind(928)=225 
      molind(1728)=226 
      molind(129)=227 
      molind(829)=228 
      molind(1629)=229 
      molind(929)=230 
      molind(1729)=231 
      molind(137)=232 
      molind(937)=233 
      molind(1737)=234 
      molind(138)=235 
      molind(838)=236 
      molind(1638)=238 
      molind(938)=239 
      molind(1738)=241 
      molind(10838)=243 
      molind(839)=246 
      molind(1639)=249 
      molind(939)=250 
      molind(1739)=252 
      molind(140)=253 
      molind(740)=255 
      molind(1640)=257 
      molind(940)=258 
      molind(1740)=260 
      molind(153)=262 
      molind(156)=263 
      molind(856)=264 
      molind(1656)=266 
      molind(956)=267 
      molind(1756)=269 
      molind(10856)=271 
      molind(757)=274 
      molind(1657)=278 
      molind(957)=279 
      molind(858)=282 
      molind(1658)=284 
      molind(958)=285 
      molind(860)=287 
      molind(960)=289 
      molind(1632)=291 
      molind(832)=292 
      molind(132)=293 
      molind(932)=294 
      molind(1732)=295 
      molind(841)=296 
      molind(303)=298 
      molind(505)=299 
      molind(1111)=300 
      molind(1212)=301 
      molind(1313)=302 
      molind(1414)=303 
      molind(1515)=304 
      molind(1919)=305 
      molind(2929)=306 
      molind(3333)=307 
      molind(3434)=308 
      molind(5151)=309 
      molind(5252)=310 
      molind(5353)=311 
      molind(5555)=312 
      molind(10199)=313 
      molind(104)=314 
      molind(122)=315 
      molind(127)=316 
      molind(130)=317 
      molind(131)=318 
      molind(133)=319 
      molind(134)=320 
      molind(135)=321 
      molind(147)=322 
      molind(148)=323 
      molind(149)=324 
      molind(150)=325 
      molind(151)=326 
      molind(152)=327 
      molind(155)=328 
      molind(170)=329 
      molind(178)=330 
      molind(179)=331 
      molind(180)=332 
      molind(181)=333 
      molind(182)=334 
      molind(183)=335 
      molind(609)=338 
      molind(634)=339 
      molind(635)=340 
      molind(645)=341 
      molind(677)=342 
      molind(678)=343 
      molind(507)=345 
      molind(709)=346 
      molind(713)=347 
      molind(717)=348 
      molind(722)=349 
      molind(733)=350 
      molind(734)=351 
      molind(831)=352 
      molind(833)=353 
      molind(834)=354 
      molind(835)=355 
      molind(837)=356 
      molind(849)=357 
      molind(850)=358 
      molind(851)=359 
      molind(852)=360 
      molind(853)=361 
      molind(865)=362 
      molind(871)=363 
      molind(872)=364 
      molind(873)=365 
      molind(874)=366 
      molind(878)=367 
      molind(882)=368 
      molind(883)=369 
      molind(890)=370 
      molind(916)=372 
      molind(930)=373 
      molind(931)=374 
      molind(933)=375 
      molind(934)=376 
      molind(935)=377 
      molind(947)=378 
      molind(948)=379 
      molind(949)=380 
      molind(950)=381 
      molind(951)=382 
      molind(953)=383 
      molind(955)=384 
      molind(967)=385 
      molind(970)=386 
      molind(971)=387 
      molind(980)=388 
      molind(981)=389 
      molind(982)=390 
      molind(311)=391 
      molind(1533)=392 
      molind(1551)=393 
      molind(416)=394 
      molind(1624)=395 
      molind(1633)=396 
      molind(1634)=397 
      molind(1650)=398 
      molind(1652)=399 
      molind(1682)=400 
      molind(1683)=401 
      molind(1721)=402 
      molind(1730)=403 
      molind(1731)=404 
      molind(1733)=405 
      molind(1734)=406 
      molind(1735)=407 
      molind(1747)=408 
      molind(1748)=409 
      molind(1749)=410 
      molind(1750)=411 
      molind(1751)=412 
      molind(1753)=413 
      molind(1755)=414 
      molind(1770)=415 
      molind(1779)=416 
      molind(1780)=417 
      molind(1781)=418 
      molind(1782)=419 
      molind(1334)=420 
      molind(1434)=421 
      molind(3234)=422 
      molind(1935)=423 
      molind(1452)=424 
      molind(3252)=425 
      molind(1953)=426 
 
      end subroutine molindx 
 
! 
!     ****************************************************** 
! 
 
