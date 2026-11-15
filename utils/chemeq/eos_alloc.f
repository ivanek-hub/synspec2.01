      subroutine eos_alloc                                                      
                                                                                
      call alloc_params                                                         
      call alloc_modelp                                                         
      call alloc_eospar                                                         
                                                                                
      end subroutine eos_alloc                                       
!                                                                               
!******************************************************************             
!                                                                               
                                                                                
      subroutine alloc_params                                                   
      use params                                                                
!                                                                               
      ALLOCATE      (YTOT(MDEPTH),                                        &     
     &              WMM(MDEPTH),                                          &     
     &              WMY(MDEPTH),                                          &     
     &              ATTOT(MATOM,MDEPTH))                                        
!                                                                               
!     Parameters for explicit atoms                                             
!                                                                               
      ALLOCATE      (AMASS(MATEX),                                        &     
     &              ABUND(MATEX,MDEPTH),                                  &     
     &              RELAB(MATEX,MDEPTH),                                  &     
     &              NUMAT(MATEX),                                         &     
     &              N0A(MATEX),                                           &     
     &              NKA(MATEX),                                           &     
     &              SABND(MATEX),                                         &     
     &              FF(MIOEX))                                                  
!                                                                               
!     Parameters for explicit ions                                              
!                                                                               
      ALLOCATE      (NFIRST(MIOEX),                                       &     
     &              NLAST(MIOEX),                                         &     
     &              NNEXT(MIOEX),                                         &     
     &              IUPSUM(MIOEX),                                        &     
     &              IZ(MIOEX),                                            &     
     &              IFREE(MIOEX),                                         &     
     &              INBFCS(MIOEX),                                        &     
     &              ILIMITS(MIOEX))                                             
!                                                                               
!     Parameters for explicit levels                                            
!                                                                               
      ALLOCATE      (ENION(MLEVEL),                                       &     
     &              G(MLEVEL))                                                  
                                                                                
      ALLOCATE      (NQUANT(MLEVEL),                                      &     
     &              IATM(MLEVEL),                                         &     
     &              IEL(MLEVEL),                                          &     
     &              ILK(MLEVEL),                                          &     
     &              ifwop(mlevel),                                        &     
     &              isemex(matom))                                              
                                                                                
      ALLOCATE      (IEXPL(MLEVEL),ILTOT(MLEVEL))                               
!                                                                               
!     Limits for explicit levels                                                
!                                                                               
      ALLOCATE      (ENION1(MLEVEL),                                      &     
     &              ENION2(MLEVEL))                                             
      ALLOCATE      (SQUANT1(MLEVEL),SQUANT2(MLEVEL),                     &     
     &              LQUANT1(MLEVEL),LQUANT2(MLEVEL),                      &     
     &              PQUANT1(MLEVEL),PQUANT2(MLEVEL))                            
!                                                                               
!     Parameters for all considered transitions                                 
!                                                                               
      ALLOCATE      (IBF(MLEVEL),                                         &     
     &              S0BF(MLEVEL),                                         &     
     &              ALFBF(MLEVEL),                                        &     
     &              BETBF(MLEVEL),                                        &     
     &              GAMBF(MLEVEL))                                              
!                                                                               
      ALLOCATE      (SGM0(MMER),                                          &     
     &              FRCH(MMER),                                           &     
     &              SGEXT1(MMER,MDEPTH),                                  &     
     &              GMER(MMER,MDEPTH),                                    &     
     &              SGMSUM(NLMX,MMER,MDEPTH),                             &     
     &              SGMG(MMER,MDEPTH))                                          
      ALLOCATE      (IMRG(MLEVEL),                                        &     
     &              IIMER(MMER))                                                
!                                                                               
      ALLOCATE      (ELEC23(MDEPTH),                                      &     
     &              Z3(MZZ),                                              &     
     &              DWC1(MZZ,MDEPTH),                                     &     
     &              DWC2(MDEPTH))                                               
!                                                                               
!     Parameters for atoms considered in line blanketing opacity                
!                                                                               
      ALLOCATE      (LGR(MATOM),LRM(MATOM))                                     
      ALLOCATE      (PFSTD(MION,MATOM),MODPF(MATOM),                      &     
     &              RR(MATOM,MION),                                       &     
     &              ENEV(MATOM,MI1),AMAS(MATOM),ABND(MATOM),              &     
     &              ABNDD(MATOM,MDEPTH),ABNREF(MDEPTH),TYPAT(MATOM))            
      ALLOCATE      (IATEX(MATOM),INPOT(MATOM,MION0),                     &     
     &              IONIZ(MATOM))                                               
      ALLOCATE      (NLEVS(MION),NLLIM(MION))                                   
                                                                                
      ALLOCATE     (FIDATA(MION),FIODF1(MION),FIODF2(MION),FIBFCS(MION))        
      ALLOCATE     (TYPLEV(MIOEX))                                              
                                                                                
!                                                                               
!     parameters for hydrogen Stark broadening tables                           
!                                                                               
      ALLOCATE      (PRFHYD(MLINH,MDEPTH,MHWL),                           &     
     &              WLHYD(MLINH,MHWL),                                    &     
     &              NWLHYD(MLINH),                                        &     
     &              WL(MHWL,MLINH),                                       &     
     &              XT(MHT,MLINH),                                        &     
     &              XNE(MHE,MLINH),                                       &     
     &              PRF(MHWL,MHT,MHE,MLINH),                              &     
     &              WLINE(4,22),                                          &     
     &              OSCH(4,22))                                                 
      ALLOCATE      (NWLH(MLINH),                                         &     
     &              NTH(MLINH),                                           &     
     &              NEH(MLINH),                                           &     
     &              ILIN0(4,22),                                          &     
     &              IHYLW(MFREQ),ILOWHW(MFREQ),                           &     
     &              M10W(MFREQ),M20W(MFREQ),                              &     
     &              IHE2LW(MFREQ),ILWHEW(MFREQ),                          &     
     &              MHE10W(MFREQ),MHE20W(MFREQ))                                
!                                                                               
!     parameters for the macroscopic velocity field and angles                  
!                                                                               
      ALLOCATE      (ANGL(MMU),WANGL(MMU),VELC(MDEPTH))                         
                                                                                
      ALLOCATE      (H0TAB(MVOI),H1TAB(MVOI),H2TAB(MVOI))                       
                                                                                
      end subroutine alloc_params                                               
                                                                                
!                                                                               
!*****************************************************************              
!                                                                               
      subroutine alloc_modelp                                                   
      use params                                                                
      use modelp                                                                
!                                                                               
!     Basic parameters of the model atmosphere                                  
!                                                                               
      ALLOCATE      (DM(MDEPTH),                                          &     
     &              TEMP(MDEPTH),                                         &     
     &              ELEC(MDEPTH),                                         &     
     &              DENS(MDEPTH),                                         &     
     &              ZD(MDEPTH),                                           &     
     &              VTURB(MDEPTH),                                        &     
     &              ABSTD(MDEPTH),                                        &     
     &              ABSTDW(MFREQC,MDEPTH),                                &     
     &              POPUL(MLEVEL,MDEPTH),                                 &     
     &              POPREL(MLEVEL,MDEPTH),                                &     
     &              DMR0(MDEPTH),                                         &     
     &              DMRP(MDEPTH),                                         &     
     &              SBF(MLEVEL),                                          &     
     &              USUM(MIOEX),                                          &     
     &              WOP(MLEVEL,MDEPTH),                                   &     
     &              WNHINT(NLMX,MDEPTH),                                  &     
     &              WNHE2(NLMX,MDEPTH),                                   &     
     &              RRR(MDEPTH,MION,MATOM),                               &     
     &              JT(MDEPTH),                                           &     
     &              TI0(MDEPTH),                                          &     
     &              TI1(MDEPTH),                                          &     
     &              TI2(MDEPTH),                                          &     
     &              VDWC(MDEPTH),                                         &     
     &              DOPA1(MATOM,MDEPTH))                                        
      ALLOCATE      (CMOL(MMOLEC))                                              
      ALLOCATE      (RRMOL(MMOLEC,MDEPTH),                                &     
     &              DOPMOL(MMOLEC,MDEPTH),                                &     
     &              AMMOL(MMOLEC),                                        &     
     &              anh2(mdepth),anch(mdepth),anoh(mdepth),               &     
     &              anhm(mdepth))                                               
!                                                                               
      ALLOCATE      (OPATM(MATOM,MFREQ,MDEPTH),                           &     
     &              EMATM(MATOM,MFREQ,MDEPTH),                            &     
     &              OPATML(MATOM,MFREQ),                                  &     
     &              GRADAT(MATOM,MDEPTH),                                 &     
     &              GRADFA(MATOM,MDEPTH),                                 &     
     &              POPAT(MATOM,MDEPTH),                                  &     
     &              DGRAD0(MATOM,MATOM,MDEPTH),                           &     
     &              DGRADP(MATOM,MATOM,MDEPTH))                                 
!                                                                               
      ALLOCATE      (RAD(MFREQ,MDEPTH),                                   &     
     &              RAD0(MFREQ,MDEPTH),                                   &     
     &              FLX0(MFREQ,MDEPTH),                                   &     
     &              flxt(mdepth),                                         &     
     &              flxi(mdepth))                                               
!                                                                               
      ALLOCATE      (PRFXB(MLINH,MHWL,MHT,MHE),                           &     
     &              PRFXR(MLINH,MHWL,MHT,MHE),                            &     
     &              PRFB(MLINH,MDEPTH,MHWL),                              &     
     &              PRFR(MLINH,MDEPTH,MHWL),                              &     
     &              ALXEN(MLINH,MHWL),                                    &     
     &              XTXEN(MHT,MLINH),                                     &     
     &              XNEXEN(MHE,MLINH))                                          
      ALLOCATE      (NWLXEN(MLINH),                                       &     
     &              NTHXEN(MLINH),                                        &     
     &              NEHXEN(MLINH),                                        &     
     &              ILXEN(4,22))                                                
                                                                                
      ALLOCATE      (CROSS(MCROSS,MFRQ))                                        
                                                                                
      ALLOCATE      (CTOP(MFIT,MCROSS),                                   &     
     &              XTOP(MFIT,MCROSS))       
                                                                                
      ALLOCATE      (CH(MFREQ,MDEPTH),ET(MFREQ,MDEPTH),                   &     
     &              SC(MFREQ,MDEPTH))                                           
                                                                                
      ALLOCATE      (CHC(MFREQC,MDEPTH),ETC(MFREQC,MDEPTH),               &     
     &              SCC(MFREQC,MDEPTH))                                         
                                                                                
      ALLOCATE      (AB(MOPAC,MDEPF),STH(MOPAC,MDEPF),SCH(MFREQC,MDEPF))        
                                                                                
      ALLOCATE      (SCCF(MFREQC,mdepf))                                        
                                                                                
      ALLOCATE      (tlbc(1000),pfbc(291,1000),eqbc(291,1000))                  
      ALLOCATE      (ttab(10000),pftab(10000))                                  
                                                                                
      ALLOCATE      (pf(32,10000))                                              
      ALLOCATE      (FROPC(MLEVEL))                                             
      ALLOCATE      (INDEXP(MLEVEL))                                            
                                                                                
      ALLOCATE      (FLUX(MFREQ),FLUXC(MFREQ))                                  
      ALLOCATE      (SCC1(MDEPTH),SCC2(MDEPTH))                                 
                                                                                
      ALLOCATE      (IREFD(MFREQ))                                              
                                                                                
      ALLOCATE      (ESEMAT(MLEVEL,MLEVEL),BESE(MLEVEL),POPLTE(MLEVEL),   &     
     &              POPUL0(MLEVEL,MDEPTH),PNLT(MATOM,MION,MDEPTH))              
                                                                                
      end subroutine alloc_modelp                                               
                                                                                
!                                                                               
!*****************************************************************              
!                                                                               
!                                                                               
                                                                                
      subroutine alloc_eospar                                                   
      use params                                                                
      use eospar                                                                
                                                                                
      ALLOCATE     (pfmol(600,mdepth),anmol(600,mdepth),                  &     
     &              pfato(100,mdepth),anato(100,mdepth),                  &     
     &              pfion(100,mdepth),anion(100,mdepth))                        
      ALLOCATE     (anion2(30,mdepth))                                          
      ALLOCATE     (ahn(mdepth),ahp(mdepth),ahen(mdepth))                       
                                                                                
      ALLOCATE     (C(600,5),PPMOL(600),APMLOG(600),P(100),               &     
     &              XIP(100),XI2(100),CCOMP(100),UIIDUI(100),             &     
     &              FP(100),XKP(100),XK2(100))                                  
      ALLOCATE     (NELEM(5,600),NATO(5,600),MMAX(600),                   &     
     &              NELEMX(100))                                                
                                                                                
      end subroutine alloc_eospar                                               
                                                                                
!                                                                               
!*****************************************************************              
!                                                                               
