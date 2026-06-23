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
     &           MLEVEL= 1650,                                            &     
     &           MDEPTH= 100,                                             &     
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
!C   &           MHWL  = 55                                               &     
     &           MHWL  = 55,                                              &     
     &           MVOI  = 2001                                                   
                                                                                
      INTEGER,PARAMETER ::                                                &     
     &           MFGRID  =  100000,                                       &     
     &           MTTAB   =      71,                                       &     
     &           MRTAB   =      20,                                       &     
     &           MSFTAB  = 6000000,                                       &     
     &           mfhtab  =    1000,                                       &     
     &           mtabth  =      10,                                       &     
     &           mtabeh  =      10                                              
!                                                                               
!     Basic physical constants                                                  
!                                                                               
      REAL(DP),PARAMETER ::                                               &     
     &           H     = 6.6256D-27,                                      &     
     &           CL    = 2.997925D10,                                     &     
     &           BOLK  = 1.38054D-16,                                     &     
     &           HK    = 4.79928144D-11,                                  &     
     &           EH    = 2.17853041D-11,                                  &     
     &           BN    = 1.4743D-2,                                       &     
     &           SIGE  = 6.6516D-25,                                      &     
     &           PI4H  = 1.8966D27,                                       &     
     &           HMASS = 1.67333D-24                                            
!                                                                               
!     Unit number                                                               
!                                                                               
      INTEGER,PARAMETER :: IBUFF=95                                             
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
      REAL(DP)  ::    TEFF,                                               &     
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
      REAL(DP)    ::  ALM00,ALST00,ALAMBE,DLAMLO                                
      INTEGER   ::  NXTSET,INLIST                                               
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
     &              N0A(:),                                               &     
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
      REAL(DP)    ::  TMOLIM,ERANGE                                             
      REAL(DP)    ::  ANEREL                                                    
      LOGICAL   ::  LASDEL                                                      
                                                                                
      INTEGER   ::  NMOLEC,IFMOL,IEQBC,                                   &     
     &              MOLTAB,IRWTAB,IIRWIN,IPFEXO,IPFBC,IPFEQ,              &     
     &              ISPICK,ILPICK,IPPICK                                        
                                                                                
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
      REAL(DP) ::     XK,FXK,BETAD,DBETA,BERGFC,CUTLYM,CUTBAL                   
      REAL(DP) ::     STHE,HGLIM                                                
      INTEGER ::    IHYDPR,IHE1PR,IHE2PR,                                 &     
     &              IHYL,ILOWH,M10,M20                                          
      INTEGER ::    NUNALP,NUNBET,NUNGAM,NUNBAL,                          &     
     &              nunhhe,ihgom,ihyddk                                         
      INTEGER,ALLOCATABLE ::                                              &     
     &              IHYLW(:),ILOWHW(:),                                   &     
     &              M10W(:),M20W(:)                                             
      INTEGER ::    IFHE2,IHE2L,ILWHE2,MHE10,MHE20                              
      INTEGER,ALLOCATABLE ::                                              &     
     &              IHE2LW(:),ILWHEW(:),                                  &     
     &              MHE10W(:),MHE20W(:)                                         
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
      REAL(DP)  ::    GSSTD,GWSTD                                               
      REAL(DP)  ::    ZND                                                       
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
      REAL(DP),ALLOCATABLE ::                                             &     
     &              PRFXB(:,:,:,:),                                       &     
     &              PRFXR(:,:,:,:),                                       &     
     &              PRFB(:,:,:),                                          &     
     &              PRFR(:,:,:),                                          &     
     &              ALXEN(:,:),                                           &     
     &              XTXEN(:,:),                                           &     
     &              XNEXEN(:,:)                                                 
      REAL(DP)    ::  XNEMIN                                                    
      INTEGER,ALLOCATABLE ::                                              &     
     &              NWLXEN(:),                                            &     
     &              NTHXEN(:),                                            &     
     &              NEHXEN(:),                                            &     
     &              ILXEN(:,:)                                                  
      INTEGER   ::  IHXENB                                                      
                                                                                
      REAL(DP),ALLOCATABLE :: CROSS(:,:)                                        
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
      module eospar                                                             

      use accura
                                                                                
      INTEGER,PARAMETER   :: NMETAL=92,NIMAX=100                               
                                                                                
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














                                                       
