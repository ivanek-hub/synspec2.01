      module accura
      integer,parameter::dp=selected_real_kind(15,300)
      end module accura

!     *******************************************************
!
      module params 
!                                                                               
!     Parameters that specify dimensions of arrays                              
!                                                                               
                                                                                
      use accura 

      INTEGER,PARAMETER ::                                                &     
     &           MATEX = 30,                                              &     
     &           MIOEX = 90,                                              &     
     &           MLEVEL= 1650  
      INTEGER :: MDEPTH  
      INTEGER,PARAMETER ::                                                &
     &           MDEPF = 500,                                             &     
     &           MFREQ = 2000,                                            &     
!    &           MFREQ =  120,                                            &     
     &           MFREQC= 2000,                                            &     
     &           MFRQ  = 2000,                                            &     
     &           MOPAC = MFRQ,                                            &     
     &           MMU   = 20,                                              &     
     &           MCROSS= MLEVEL,                                          &     
     &           MFIT  = 1650,                                            &     
     &           MFCRA =  1200,                                           &     
     &           MTRAD =  3,                                              &     
     &           MATOM =  99,                                             &     
     &           MATOMBIG = 99,                                           &     
     &           MION  = 90,                                              &     
     &           MION0 =  9,                                              &     
     &           MMOLEC=500,                                              &     
     &           MPHOT = 10,                                              &     
     &           MZZ   =  2,                                              &     
     &           MMER  =  2,                                              &     
     &           NLMX  = 80,                                              &     
     &           MI1   = MION0-1                                                
      INTEGER,PARAMETER ::                                                &     
     &           MLINH = 78,                                              &     
     &           MHT   = 7,                                               &     
     &           MHE   = 20,                                              &     
     &           MHWL  = 55,                                              &     
     &           MVOI  = 2001                                                   
                                                                                
      INTEGER :: MFGRID
      INTEGER :: MTTAB
      INTEGER :: MSFTAB
      INTEGER,PARAMETER ::                                                &     
     &           MRTAB   =      20,                                       &     
     &           mfhtab  =    1000,                                       &     
     &           mtabth  =      10,                                       &     
     &           mtabeh  =      10                                              
!                                                                               
!     Basic physical constants                                                  
!                                                                               
      REAL(DP),PARAMETER ::                                               &     
     &           H     = 6.6256e-27,                                      &     
     &           CL    = 2.997925e10,                                     &     
     &           BOLK  = 1.38054e-16,                                     &     
     &           HK    = 4.79928144e-11,                                  &     
     &           EH    = 2.17853041e-11,                                  &     
     &           BN    = 1.4743e-2,                                       &     
     &           SIGE  = 6.6516e-25,                                      &     
     &           PI4H  = 1.8966e27,                                       &     
     &           HMASS = 1.67333e-24                                            
!                                                                               
!     Unit number                                                               
!                                                                               
      INTEGER,PARAMETER :: IBUFF=95
      INTEGER ::           IF55,IFKEY                                             
!                                                                               
!     Basic parameters                                                          
!                                                                               
      INTEGER ::    NATOM,                                                &     
     &              NION,                                                 &     
     &              NLEVEL,                                               &     
     &              ND,NDSTEP,                                            &     
     &              NFREQ,NFROBS,NFREQC,NFREQS,                           &     
     &              NMU                                                         
      LOGICAL ::    LTE,LTEGR                                                   
      REAL(DP)  ::  TEFF,                                                 &     
     &              GRAV,                                                 &     
     &              vaclim                                                      
      REAL(DP),ALLOCATABLE ::                                             &     
     &              YTOT(:),                                              &     
     &              WMM(:),                                               &     
     &              WMY(:),                                               &     
     &              ATTOT(:,:)                                                  
      INTEGER ::    IMODE,                                                &     
     &              IMODE0,                                               &     
     &              IFREQ,                                                &     
     &              INLTE,                                                &     
     &              IDSTD,                                                &     
     &              IFWIN,                                                &     
     &              IFEOS,                                                &     
     &              IBFAC,                                                &     
     &              INMOD,INTRPL,ICHANG,ICHEMC,IATREF,ICONTL,             &     
     &              IBLANK,NBLANK,                                        &     
     &              IPRIN,                                                &     
     &              IGRDD,IRELIN                                                
      REAL(DP)  ::  ALAM0,ALAST,CUTOF0,CUTOFS,RELOP,SPACE
      REAL(DP)  ::  ALM00,ALST00,ALAMBE,DLAMLO                                
      INTEGER   ::  NXTSET,INLIST,NMLIS0                   
      REAL      ::  DTIM
!                                                                               
!     Parameters for explicit atoms                                             
!                                                                               
      REAL(DP),ALLOCATABLE ::                                             &     
     &              AMASS(:),                                             &     
     &              ABUND(:,:),                                           &     
     &              RELAB(:,:),                                           &     
     &              SABND(:),                                             &     
     &              FF(:)                                                       
      INTEGER,ALLOCATABLE ::                                              &     
     &              NUMAT(:),                                             &     
     &             N0A(:),                                               &     
     &              NKA(:)                                                      
!                                                                               
!     Parameters for explicit ions                                              
!                                                                               
      INTEGER,ALLOCATABLE ::                                              &     
     &              NFIRST(:),                                            &     
     &              NLAST(:),                                             &     
     &              NNEXT(:),                                             &     
     &              IUPSUM(:),                                            &     
     &              IZ(:),                                                &     
     &              IFREE(:),                                             &     
     &              INBFCS(:),                                            &     
     &              ILIMITS(:)                                                  
                                                                                
      CHARACTER(LEN=10),ALLOCATABLE :: TYPLEV(:)                                
!                                                                               
!     Parameters for explicit levels                                            
!                                                                               
      REAL(DP),ALLOCATABLE ::                                             &     
     &              ENION(:),                                             &     
     &              G(:)                                                        
      INTEGER,ALLOCATABLE ::                                              &     
     &              NQUANT(:),                                            &     
     &              IATM(:),                                              &     
     &              IEL(:),                                               &     
     &              ILK(:),                                               &     
     &              ifwop(:),                                             &     
     &              isemex(:)                                                   
      INTEGER    :: NMER                                                        
                                                                                
      INTEGER,ALLOCATABLE :: IEXPL(:),ILTOT(:)                                  
!                                                                               
!     Limits for explicit levels                                                
!                                                                               
      REAL(DP),ALLOCATABLE ::                                             &     
     &              ENION1(:),                                            &     
     &              ENION2(:)                                                   
      INTEGER*4,ALLOCATABLE ::                                            &     
     &              SQUANT1(:),SQUANT2(:),                                &     
     &              LQUANT1(:),LQUANT2(:),                                &     
     &              PQUANT1(:),PQUANT2(:)                                       
!                                                                               
!     Parameters for all considered transitions                                 
!                                                                               
      INTEGER,ALLOCATABLE ::                                              &     
     &              IBF(:)                                                      
      REAL(DP),ALLOCATABLE ::                                             &     
     &              S0BF(:),                                              &     
     &              ALFBF(:),                                             &     
     &              BETBF(:),                                             &     
     &              GAMBF(:)                                                    
                                                                                
      REAL(DP),ALLOCATABLE ::                                             &     
     &              SGM0(:),                                              &     
     &              FRCH(:),                                              &     
     &              SGEXT1(:,:),                                          &     
     &              GMER(:,:),                                            &     
     &              SGMSUM(:,:,:),                                        &     
     &              SGMG(:,:)                                                   
      INTEGER,ALLOCATABLE ::                                              &     
     &              IMRG(:),                                              &     
     &              IIMER(:)                                                    
!                                                                               
      REAL(DP),ALLOCATABLE ::                                             &     
     &              ELEC23(:),                                            &     
     &              Z3(:),                                                &     
     &              DWC1(:,:),                                            &     
     &              DWC2(:)                                                     
!                                                                               
!     additional opacities                                                      
!                                                                               
      INTEGER    :: IOPADD,                                               &     
     &              IOPHMI,                                               &     
     &              IOPH2P,                                               &     
     &              IOPHEM,                                               &     
     &              IOPCH,                                                &     
     &              IOPOH,                                                &     
     &              IOPH2M,                                               &     
     &              IOH2H2,IOH2HE,IOH2H1,IOHHE,                           &     
     &              IOPHLI,                                               &     
     &              IRSCT,                                                &     
     &              IRSCHE,                                               &     
     &              IRSCH2                                                      
!                                                                               
!     Auxiliary parameters                                                      
!                                                                               
      INTEGER   ::  IATH,IELH,IELHM,N0H,N1H,NKH,N0HN,N0M,                 &     
     &              IATHE,IELHE1,IELHE2                                         
      REAL(DP)  ::  TMOLIM,ERANGE                                             
      REAL(DP)  ::  ANEREL                                                    
      LOGICAL   ::  LASDEL                                                      
                                                                                
      INTEGER   ::  NMOLEC,IFMOL,IEQBC,                                   &     
     &              MOLTAB,IRWTAB,IIRWIN,IPFEXO,IPFBC,IPFEQ,              &     
     &              ISPICK,ILPICK,IPPICK                                        
      INTEGER   ::  NMLI0,IUNIM1,IUNIM2
                                                                                
      CHARACTER(LEN=40),ALLOCATABLE :: FIDATA(:),                         &     
     &              FIODF1(:),FIODF2(:),FIBFCS(:)                               
!                                                                               
!     Parameters for atoms considered in line blanketing opacity                
!                                                                               
      LOGICAL,ALLOCATABLE ::  LGR(:),LRM(:)                                     
      REAL(DP),ALLOCATABLE  ::                                            &     
     &            PFSTD(:,:),RR(:,:),                                     &     
     &            ENEV(:,:),AMAS(:),ABND(:),                              &     
     &            ABNDD(:,:),ABNREF(:)                                          
      CHARACTER(LEN=4),ALLOCATABLE :: TYPAT(:)                                  
      INTEGER,ALLOCATABLE ::                                              &     
     &            MODPF(:),IATEX(:),INPOT(:,:),IONIZ(:),                  &     
     &            NLEVS(:),NLLIM(:)                                             
      INTEGER ::  NATOMS                                                        
!                                                                               
!     parameters for hydrogen Stark broadening tables                           
!                                                                               
      REAL(DP),ALLOCATABLE ::                                             &     
     &              PRFHYD(:,:,:),                                        &     
     &              WLHYD(:,:),                                           &     
     &              WL(:,:),                                              &     
     &              XT(:,:),                                              &     
     &              XNE(:,:),                                             &     
     &              PRF(:,:,:,:),                                         &     
     &              WLINE(:,:),                                           &     
     &              OSCH(:,:)                                                   
      INTEGER,ALLOCATABLE ::                                              &     
     &              NWLHYD(:),                                            &     
     &              NWLH(:),                                              &     
     &              NTH(:),                                               &     
     &              NEH(:),                                               &     
     &              ILIN0(:,:)                                                  
      INTEGER    :: ILEMKE,                                               &     
     &              NLIHYD                                                      
      REAL(DP) ::   XK,FXK,BETAD,DBETA,BERGFC,CUTLYM,CUTBAL                   
      REAL(DP) ::   STHE,HGLIM                                                
      INTEGER ::    IHYDPR,IHE1PR,IHE2PR,                                 &     
     &              IHYL,ILOWH,M10,M20                                          
      INTEGER ::    NUNALP,NUNBET,NUNGAM,NUNBAL,                          &     
     &              nunhhe,ihgom,ihyddk                                         
!     INTEGER,ALLOCATABLE ::                                              &     
!    &              IHYLW(:),ILOWHW(:),                                   &     
!    &              M10W(:),M20W(:)                                             
      INTEGER ::    IFHE2,IHE2L,ILWHE2,MHE10,MHE20                              
!     INTEGER,ALLOCATABLE ::                                              &     
!    &              IHE2LW(:),ILWHEW(:),                                  &     
!    &              MHE10W(:),MHE20W(:)                                         
!                                                                               
!     parameters for the macroscopic velocity field and angles                  
!                                                                               
      REAL(DP),ALLOCATABLE ::                                             &     
     &              ANGL(:),WANGL(:),VELC(:)                                    
      INTEGER ::    NMU0,IFLUX                                                  
                                                                                
                                                                                
      REAL(DP),ALLOCATABLE ::                                             &     
     &              ESEMAT(:,:),BESE(:),POPLTE(:),POPUL0(:,:),            &     
     &              PNLT(:,:,:)                                                 
!                                                                               
!     parameters for the "standard" molecula Strak and vdWaals broadening       
!                                                                               
      REAL(DP)  ::  GSSTD,GWSTD                                               
      REAL(DP)  ::  ZND                                                       
      INTEGER ::    IFZ0                                                        
                                                                                
      REAL(DP),ALLOCATABLE  :: H0TAB(:),H1TAB(:),H2TAB(:)                       
                                                                                
      end module params                                                         
!                                                                         &     
!*****************************************************************              
!                                                                               
      module modelp                                                             

      use accura 
      use params                                                                
!                                                                               
!     Basic parameters of the model atmosphere                                  
!                                                                               
      REAL(DP),ALLOCATABLE ::                                             &     
     &              DM(:),                                                &     
     &              TEMP(:),                                              &     
     &              ELEC(:),                                              &     
     &              DENS(:),                                              &     
     &              ZD(:),                                                &     
     &              VTURB(:),                                             &     
     &              ABSTD(:),                                             &     
     &              ABSTDW(:,:),                                          &     
     &              POPUL(:,:),                                           &     
     &              POPREL(:,:),                                          &     
     &              DMR0(:),                                              &     
     &              DMRP(:),                                              &     
     &              SBF(:),                                               &     
     &              USUM(:),                                              &     
     &              WOP(:,:),                                             &     
     &              WNHINT(:,:),                                          &     
     &              WNHE2(:,:),                                           &     
     &              RRR(:,:,:),                                           &     
     &              TI0(:),                                               &     
     &              TI1(:),                                               &     
     &              TI2(:),                                               &     
     &              VDWC(:),                                              &     
     &              DOPA1(:,:)                                                  
      REAL(DP)     :: VTB                                                       
      INTEGER,ALLOCATABLE :: JT(:)                                              
      CHARACTER(LEN=8),ALLOCATABLE :: CMOL(:)                                   
      REAL(DP),ALLOCATABLE ::                                             &     
     &              RRMOL(:,:),                                           &     
     &              DOPMOL(:,:),                                          &     
     &              AMMOL(:)                                                    
                                                                                
      REAL(DP),ALLOCATABLE ::                                             &     
     &              anh2(:),anch(:),anoh(:),                              &     
     &              anhm(:)                                                     
      REAL(DP)     :: HPOP                                                      
!                                                                               
      REAL(DP),ALLOCATABLE ::                                             &     
     &              OPATM(:,:,:),                                         &     
     &              EMATM(:,:,:),                                         &     
     &              OPATML(:,:),                                          &     
     &              GRADAT(:,:),                                          &     
     &              GRADFA(:,:),                                          &     
     &              POPAT(:,:),                                           &     
     &              DGRAD0(:,:,:),                                        &     
     &              DGRADP(:,:,:)                                               
!                                                                               
      REAL(DP),ALLOCATABLE ::                                             &     
     &              RAD(:,:),                                             &     
     &              RAD0(:,:),                                            &     
     &              FLX0(:,:),                                            &     
     &              flxt(:),                                              &     
     &              flxi(:)                                                     
!                                                                               
      REAL(DP)    ::  XNEMIN                                                    
      INTEGER,ALLOCATABLE ::                                              &     
     &              NWLXEN(:),                                            &     
     &              ILXEN(:,:)                                                  
      INTEGER   ::  IHXENB                                                      
                                                                                
      REAL(DP),ALLOCATABLE :: CROSS(:,:),DENSCON(:)
      REAL(DP),ALLOCATABLE :: CTOP(:,:),XTOP(:,:)                               
                                                                                
      REAL(DP),ALLOCATABLE :: CH(:,:),ET(:,:),SC(:,:)                           
      REAL(DP),ALLOCATABLE :: CHC(:,:),ETC(:,:),SCC(:,:)                        
                                                                                
      REAL(DP),ALLOCATABLE :: AB(:,:),STH(:,:),SCH(:,:)                         
      REAL(DP),ALLOCATABLE :: SCCF(:,:)                                         
                                                                                
      REAL(DP),ALLOCATABLE :: TLBC(:),PFBC(:,:),EQBC(:,:)                       
                                                                                
      REAL(DP),ALLOCATABLE :: TTAB(:),PFTAB(:)                                  
                                                                                
      REAL(DP),ALLOCATABLE :: PF(:,:)                                           
      REAL(DP),ALLOCATABLE :: FROPC(:)                                          
      INTEGER,ALLOCATABLE :: INDEXP(:)                                          
                                                                                
      REAL(DP),ALLOCATABLE :: FLUX(:),FLUXC(:)                                  
      REAL(DP),ALLOCATABLE :: SCC1(:),SCC2(:)                                   
                                                                                
      INTEGER,ALLOCATABLE :: IREFD(:)                                           
                                                                                
      end module modelp                                                         
                                                                                
!                                                                         &     
!*****************************************************************              
!                                                                               
      module lindat                                                             

      use accura
      use params                                                                
                                                                                
      INTEGER,PARAMETER ::                                                &
     &                     MGRIEM  =     10,                              &     
     &                     MNLT    =   2000,                              &     
     &                     MSPHE2  =     20,                              &     
     &                     MMLIN0  =     40                                     
!                                                                               
      INTEGER            ::  MMLIST,NMLIST                                      
      INTEGER            ::  MLIN0,MLIN
                                                                                
       REAL*4,ALLOCATABLE ::                                              &     
     &              EXCL0(:),                                             &     
     &              EXCU0(:),                                             &     
     &              GF0(:),                                               &     
     &              EXTIN(:),                                             &     
     &              BNUL(:),                                              &     
     &              GAMR0(:),                                             &     
     &              GS0(:),                                               &     
     &              GW0(:),                                               &     
     &              WGR0(:,:)                                                   
!                                                                               
      REAL(DP),ALLOCATABLE ::                                             &     
     &              FREQ0(:)                                                    
      INTEGER,ALLOCATABLE ::                                              &     
     &              INDAT(:),                                             &     
     &              INDNLT(:),                                            &     
     &              ILOWN(:),                                             &     
     &              IUPN(:),                                              &     
     &              IJCONT(:),                                            &     
     &              IJCNTR(:),                                            &     
     &              INDLIN(:),                                            &     
     &              INDLIP(:),                                            &     
     &              IJCTR(:)                                                    
                                                                                
      INTEGER ::    NLIN0,NLIN,IRLIST,                                    &     
     &              NNLT,NGRIEM                                                 
!                                                                               
!     INTEGER,ALLOCATABLE :: IBIN(:)                                            
!                                                                               
      INTEGER,ALLOCATABLE ::                                              &     
     &              IPRF0(:),                                             &     
     &              ISPRF(:),                                             &     
     &              IGRIEM(:),                                            &     
     &              ISP0(:)                                                     
      INTEGER ::    NSP                                                         
!                                                                               
      REAL(DP),ALLOCATABLE ::                                             &     
     &              ABCENT(:,:),                                          &     
     &              SLIN(:,:),                                            &     
     &              PLAN(:),                                              &     
     &              STIM(:),                                              &     
     &              EXHK(:)                                                     
!                                                                               
      REAL(DP)  ::  ALAM1,FRMIN,FRLAST,FRLI0,FRLIM,                       &     
     &              SPACE0,TSTD,DSTD,ALAMC,APREV                   
      REAL(DP)  ::  ALAM0s,ALASTs,CUTOF0s,CUTOFSs,RELOPs,SPACEs               
!                                                                               
      REAL(DP)  ::  DFRCON                                                    
                                                                                
      CHARACTER(LEN=40)  ::  AMLIST(0:MMLIN0)                                   
      INTEGER            ::  IBIN(0:MMLIN0)                                     
                                                                                
      end module lindat                                                         
!                                                                         &     
!*****************************************************************              
!                                                                               
      module synthp                                                             

      use accura
      use params                                                                
                                                                                
      REAL(DP),ALLOCATABLE ::                                             &     
     &              FREQ(:),W(:),WLAM(:),                                 &     
     &              FRX1(:),FRX2(:),BNUE(:),                              &     
     &              FRQOBS(:),WLOBS(:),                                   &     
     &              FREQC(:),WLAMC(:),                                    &     
     &              ABSOC(:),EMISC(:),SCATC(:),PLAC(:)                          
                                                                                
                                                                                
      INTEGER,ALLOCATABLE :: IJCINT(:)                                          
                                                                                
      REAL(DP),ALLOCATABLE ::                                             &     
     &              FRECR(:,:),CROSR(:,:),CRMX(:),                        &     
     &              PHOTI(:,:)                                                  
      INTEGER,ALLOCATABLE :: NFCR(:)                                            
      INTEGER ::    IASV                                                        
                                                                                
      REAL(DP),ALLOCATABLE ::                                             &     
     &              FRECQ(:,:),QHOT(:,:),                                 &     
     &              AQHT(:),EQHT(:),GQHT(:),                              &     
     &              CRMY(:)                                                     
      INTEGER,ALLOCATABLE :: NFQHT(:)                                           
      INTEGER             :: NQHT                                               
                                                                                
      end module synthp                                                         
!                                                                         &     
!*****************************************************************              
!                                                                               
      module wincom                                                             

      use accura
      use params                                                                
                                                                                
      INTEGER,PARAMETER  :: MRCORE=20
      INTEGER            :: MKU,                                          &     
     &                      MEXT                                                
!!                          MKU=MDEPTH=MCORE; MEXT=MKU
                                                 
      REAL(DP),ALLOCATABLE :: BMU(:,:),WMUJ(:,:),WMUH(:),                 &     
     &                      RD(:),PIM(:),RAD1(:),                         &     
     &                      DELZ(:,:)                                           
      REAL(DP)             :: RCORE,RFNORM                                      
      INTEGER,ALLOCATABLE:: NUD(:),NUDF(:)                                      
      INTEGER            :: KMU,NREXT,NRCORE,NFIRY,NDF                          
      REAL(DP),ALLOCATABLE :: DELZF(:,:),DFRQF(:,:),                      &     
     &                      VEL(:),DFRQ(:,:),DVD(:)                             
      REAL(DP)    ::          XMDOT,XMD4,BETAV,VINF                             
      REAL(DP),ALLOCATABLE :: FFQ(:),FFQV(:),RDF(:),DENSF(:),             &     
     &                      VELF(:,:),DRAY(:,:)                                 
      INTEGER,ALLOCATABLE:: KRAY(:,:)                                           
      REAL(DP),ALLOCATABLE :: WDIL(:),PLANW(:),TRAD(:,:)  
      INTEGER,ALLOCATABLE:: ILNE(:),ILVI(:)                                     
      INTEGER            :: NOPAC                                               
      INTEGER            :: IEMOFF,NLTOFF,ITRAD                                 
      REAL(DP)             :: VELMAX                                            
                                                                                
      INTEGER,ALLOCATABLE ::                                              &     
     &              IHYLW(:),ILOWHW(:),                                   &     
     &              M10W(:),M20W(:)                                             
      INTEGER,ALLOCATABLE ::                                              &     
     &              IHE2LW(:),ILWHEW(:),                                  &     
     &              MHE10W(:),MHE20W(:)                                         

      end module wincom                                                         
!                                                                         &     
!*****************************************************************              
!                                                                               
      module optabl                                                             

      use accura

      REAL, ALLOCATABLE    :: ABSGRD(:,:,:)                                       

      REAL(DP),ALLOCATABLE :: WLGRID(:),                                  &     
     &                      ABSOP(:),WLTAB(:),ABGRD(:),XLI(:)                   
      REAL(DP),ALLOCATABLE :: TEMPG(:),DENSG(:,:),ELECGR(:,:),            &     
     &                      DENSG0(:),ELECM(:),RELABN(:),                 &     
     &                      YINT(:)                                             
      REAL(DP)             :: TEMP1                                             
      INTEGER,ALLOCATABLE:: NDEN(:),NFRTAB(:,:),JINT(:)                         
      INTEGER            :: NTEMP,NDENS,                                  &     
     &                      NFGRID,IPFREQ,INDEXT,INDEXN,                  &     
     &                      INTTAB,IBINGR,IDENS,IDENS0                          
      CHARACTER(LEN=80)     TABNAME                                             
                                                                                
      end module optabl                                                         
!                                                                         &     
!****&************************************************************              
!                                                                               
      module molist                                                             

      use accura
                                                                                
      INTEGER, PARAMETER ::                                               &  
     &                       MMLINI =      40                                   
                                                                                
      INTEGER             :: MLINM0,MLINM                                                                          
      INTEGER             :: MBROAD,NBROAD,                               &     
     &                       MVDWLI,NVDWLI                                      
                                                                                
      REAL(DP), ALLOCATABLE :: FREQM(:,:)                                       
      REAL*4, ALLOCATABLE :: EXCLM(:,:),                                  &     
     &                       GFM(:,:),                                    &     
     &                       EXTINM(:,:),                                 &     
     &                       GRM(:,:),                                    &     
     &                       GSM(:,:),                                    &     
     &                       GWM(:,:),                                    &     
     &                       GVDWH2(:,:),                                 &     
     &                       GEXPH2(:,:),                                 &     
     &                       GVDWHE(:,:),                                 &     
     &                       GEXPHE(:,:)                                        
                                                                                
      INTEGER,ALLOCATABLE :: INDATM(:,:),                                 &     
     &                       INMLIN(:,:),                                 &     
     &                       INMLIP(:,:),                                 &     
     &                       IJCMTR(:,:),                                 &     
     &                       NLINM0(:),                                   &     
     &                       NLINML(:),                                   &     
     &                       NLINMT(:),                                   &     
     &                       INACTM(:),                                   &     
     &                       NXTSEM(:),                                   &     
     &                       IPRSEM(:),                                   &     
     &                       IREADM(:),                                   &
     &                       ILASTM(:)                                          
      INTEGER,ALLOCATABLE :: MOLIND(:),                                   &     
     &                       IONIND(:)                                          
                                                                                
      REAL(DP),ALLOCATABLE ::  FRLASM(:),                                 &     
     &                       ALASTM(:),                                   &     
     &                       ALEND(:)                                           
                                                                                
!     CHARACTER(LEN=40)  ::  AMLIST(0:MMLINI)                                   
!     INTEGER            ::  IBIN(0:MMLINI),                                    
      INTEGER                IVDWLI(MMLINI),                              &     
     &                       IBROLI(MMLINI),                              &     
     &                       NMPAR(MMLINI),                               &     
     &                       IUNITM(MMLINI)                                     
      REAL(DP)             ::  TMLIM(MMLINI)                                    
                                                                                
      end module molist                                                         
!                                                                         &     
!*****************************************************************              
!                                                                               
      module eospar                                                             

      use accura
                                                                                
      INTEGER,PARAMETER   :: NMETAL=92,NIMAX=200                               
                                                                                
      REAL(DP), ALLOCATABLE :: pfmol(:,:),anmol(:,:),                     &     
     &                       pfato(:,:),anato(:,:),                       &     
     &                       pfion(:,:),anion(:,:),                       &     
     &                       anion2(:,:)                                        
      REAL(DP), ALLOCATABLE :: ahn(:),ahp(:),ahen(:)                            
      REAL(DP)              :: anhmi,ahmol                                      
                                                                                
      REAL(DP), ALLOCATABLE :: C(:,:),PPMOL(:),APMLOG(:),P(:),            &     
     &                       XIP(:),XI2(:),CCOMP(:),UIIDUI(:),            &     
     &                       FP(:),XKP(:),XK2(:)                                
      REAL(DP)              :: EPS,SWITER                                       
      INTEGER,ALLOCATABLE :: NELEM(:,:),NATO(:,:),MMAX(:),                &     
     &                       NELEMX(:)                                          
                                                                                
      end module eospar                                                         
!                                                                         &     
!*****************************************************************              
!                                                                               
      module heprf                                                              

      use accura
                                                                                
      REAL(DP), ALLOCATABLE :: PRFHE2(:,:,:),WLHE2(:,:)                         
      REAL(DP), ALLOCATABLE :: WL2(:,:),XT2(:),XNE2(:,:),PRF2(:,:,:)            
                                                                                
      INTEGER,ALLOCATABLE :: NWLHE2(:),ILHE2(:),IUHE2(:)                        
      INTEGER             :: NWL2,NT2,NE2                                       
                                                                                
      REAL(DP), ALLOCATABLE :: PRFHE1(:,:,:,:),DLMHE1(:,:,:),XNEHE1(:),   &     
     &                       PRF447(:,:,:),DLM447(:,:),XNE447(:)                
      INTEGER,ALLOCATABLE :: NWLAM(:,:)                                         
                                                                                
      end module heprf                                                          
                                                                                
!                                                                         &     
!*****************************************************************              
!                                                                               
                                                                                
      module opadat                                                             

      use accura
                                                                                
      integer, parameter  :: mtemp=100,melec=60,mion1=30                        
      REAL(DP), allocatable :: frac(:,:,:),fracm(:,:)                           
      integer,allocatable :: itemp(:)                                           
      integer             :: ntt                                                
                                                                                
      REAL(DP), allocatable :: pfunc(:,:,:),frpf0(:,:)                          
                                                                                
      REAL(DP), allocatable :: alogtt(:), alogrhot(:,:), fract(:,:,:,:),  &     
     &    zct(:,:,:,:), felect(:,:,:), alognet(:,:,:),                    &     
     &    fractH2(:,:),fractH2p(:,:),fracHm(:,:),                         &     
     &    zctH2(:,:),zctH2p(:,:),zctHm(:,:),                              &     
     &    zcthyd(:,:),zctHp(:,:)                                                
                                                                                
      integer,allocatable :: nrho(:),ntempe(:)                                  
                                                                                
      end module opadat                                                         
                                                                                
!                                                                         &     
!*****************************************************************              
!                                                                               
                                                                                
      module photcs                                                             

      use accura
                                                                                
      REAL(DP), allocatable :: phot(:,:),apht(:),epht(:),gpht(:)                
      REAL(DP)              :: wpht0,wpht1                                      
                                                                                
      integer,allocatable :: jpht(:)                                            
      integer             :: npht                                               
                                                                                
      end module photcs                                                         
                                                                                
!                                                                         &     
!*****************************************************************              
!                                                                               
                                                                                
      module hydprf                                                             

      use accura
                                                                                
      integer, parameter :: mtab=5,                                       &     
     &                      mden=10,mtem=10,mlam=1000                           
                                                                                
      REAL(DP),allocatable :: tabden(:,:),tabtem(:,:,:),                  &     
     &                      tablam(:,:,:),                                &     
     &                      prftab(:,:,:,:)                                     
                                                                                
      integer,allocatable:: nden(:),numtem(:,:),numlam(:,:)                     
                                                                                
      REAL(DP), allocatable :: frgtab(:),wlgtab(:),hydopg(:,:)                  
      integer             :: nugfreq                                            
                                                                                
      end module hydprf                                                         
                                                                                
!                                                                         &     
!*****************************************************************              
!                                                                               
                                                                                
      module topdat                                                             

      use accura
                                                                                
      integer, parameter  :: mmaxop=200,                                   &    
     &                       mop   = 15                                         
                                                                                
      REAL(DP), allocatable :: sop(:,:),xop(:,:)                                
      integer, allocatable:: nop(:)  ! number of fit points for current lev.    
      integer             :: ntotop  ! total number of levels in OP data        
      logical             :: loprea  ! .T. OP data read in; .F. not yet         
      character(len=10),allocatable :: idlvop(:)  ! level identifyer            
                                                                                
      end module topdat 

!                                                                         &
!*****************************************************************
!

      module allarn 

      use accura

      integer, parameter :: nxmax=1400,                                   &
     &                      nnmax=   5

       real(dp), allocatable :: xlalp(:),plalp(:,:)
       real(dp)              :: stnnea,stncha,vneua,vchaa
       integer               :: nxalp,iwarna 

       real(dp), allocatable :: xlbet(:),plbet(:,:)
       real(dp)              :: stnneb,stnchb,vneub,vchab
       integer               :: nxbet,iwarnb 

       real(dp), allocatable :: xlgam(:),plgam(:,:)
       real(dp)              :: stnneg,stnchg,vneug,vchag
       integer               :: nxgam,iwarng 

       real(dp), allocatable :: xlbal(:),plbal(:,:)
       real(dp)              :: stnnec,stnchc,vneuc,vchac
       integer               :: nxbal,iwarnc 

       end module allarn
!                                                                         &
!*****************************************************************        
!    

      module hydxen

      use accura

      real(dp), allocatable ::  prfxb(:,:,:,:),                           &
     &                          prfxr(:,:,:,:),                           &
     &                          PRFB(:,:,:),                              &     
     &                          PRFR(:,:,:),                              &     
     &                          ALXEN(:,:),                               &     
     &                          xtxen(:,:),xnexen(:,:)
       integer, allocatable ::  nthxen(:),nehxen(:)

       end module hydxen













                                                       
