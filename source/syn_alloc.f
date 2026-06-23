      subroutine syn_alloc                                                      
                                                                                
      call alloc_params                                                         
      call alloc_modelp                                                         
      call alloc_lindat                                                         
      call alloc_synthp                                                         
      call alloc_eospar                                                         
 
      end subroutine syn_alloc                                                  
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
!     ALLOCATE      (SGM0(MMER),                                          &     
!    &              FRCH(MMER),                                           &     
!    &              SGEXT1(MMER,MDEPTH),                                  &     
      ALLOCATE     (GMER(MMER,MDEPTH))                                    &     
!    &              SGMSUM(NLMX,MMER,MDEPTH),                             &     
!    &              SGMG(MMER,MDEPTH))                                          
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
     &              ILIN0(4,22))                                                
!    &              IHYLW(MFREQ),ILOWHW(MFREQ),                           &     
!    &              M10W(MFREQ),M20W(MFREQ),                              &     
!    &              IHE2LW(MFREQ),ILWHEW(MFREQ),                          &     
!    &              MHE10W(MFREQ),MHE20W(MFREQ))                                
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
!    &              POPREL(MLEVEL,MDEPTH),                                &     
!    &              DMR0(MDEPTH),                                         &     
!    &              DMRP(MDEPTH),                                         &     
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
!     ALLOCATE      (OPATM(MATOM,MFREQ,MDEPTH),                           &     
!    &              EMATM(MATOM,MFREQ,MDEPTH),                            &     
!    &              OPATML(MATOM,MFREQ),                                  &     
!    &              GRADAT(MATOM,MDEPTH),                                 &     
!    &              GRADFA(MATOM,MDEPTH),                                 &     
!    &              POPAT(MATOM,MDEPTH),                                  &     
!    &              DGRAD0(MATOM,MATOM,MDEPTH),                           &     
!    &              DGRADP(MATOM,MATOM,MDEPTH))                                 
!                                                                               
!     ALLOCATE      (RAD(MFREQ,MDEPTH),                                   &     
!    &              RAD0(MFREQ,MDEPTH),                                   &     
!    &              FLX0(MFREQ,MDEPTH),                                   &     
!    &              flxt(mdepth),                                         &     
!    &              flxi(mdepth))                                               
!                                                                               
!     ALLOCATE     (PRFB(MLINH,MDEPTH,MHWL),                              &     
!    &              PRFR(MLINH,MDEPTH,MHWL),                              &     
!    &              ALXEN(MLINH,MHWL))                                          
      ALLOCATE      (NWLXEN(MLINH),                                       &     
     &              ILXEN(4,22))                                                
                                                                                
      ALLOCATE      (CROSS(MCROSS,MFRQ))                                        
      ALLOCATE      (DENSCON(MDEPTH))
                                                                                
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
      subroutine alloc_lindat                                                   
      use params                                                                
      use lindat                                                                
                                                                                
      MPRF=MLIN0

       ALLOCATE    (EXCL0(MLIN0),                                         &     
     &              EXCU0(MLIN0),                                         &     
     &              GF0(MLIN0),                                           &     
     &              EXTIN(MLIN0),                                         &     
     &              BNUL(MLIN0),                                          &     
     &              GAMR0(MPRF),                                          &     
     &              GS0(MPRF),                                            &     
     &              GW0(MPRF),                                            &     
     &              WGR0(4,MGRIEM))                                             
!                                                                               
      ALLOCATE     (FREQ0(MLIN0),                                         &     
     &              INDAT(MLIN0),                                         &     
     &              INDNLT(MLIN0),                                        &     
     &              ILOWN(MLIN0),                                         &     
     &              IUPN(MLIN0),                                          &     
     &              IJCONT(MLIN0),                                        &     
     &              IJCNTR(MLIN),                                         &     
     &              INDLIN(MLIN),                                         &     
     &              INDLIP(MLIN),                                         &     
     &              IJCTR(MFREQ))                                               
!                                                                               
      ALLOCATE     (IPRF0(MPRF),                                          &     
     &              ISPRF(MPRF),                                          &     
     &              IGRIEM(MPRF),                                         &     
     &              ISP0(MSPHE2))                                               
!                                                                               
      ALLOCATE     (ABCENT(MNLT,MDEPTH),                                  &     
     &              SLIN(MNLT,MDEPTH))                                          
!                                                                               
      ALLOCATE     (PLAN(MDEPTH),                                         &     
     &              STIM(MDEPTH),                                         &     
     &              EXHK(MDEPTH))                                               
!                                                                               
      end subroutine alloc_lindat                                               
!                                                                               
!*****************************************************************              
!                                                                               
      subroutine alloc_synthp                                                   
      use params                                                                
      use synthp                                                                
                                                                                
      ALLOCATE      (FREQ(MFREQ),W(MFREQ),WLAM(MFREQ),                    &     
     &              FRX1(MFREQ),FRX2(MFREQ),BNUE(MFREQ),                  &     
     &              FRQOBS(MFREQ),WLOBS(MFREQ),                           &     
     &              FREQC(MFREQC),WLAMC(MFREQC),                          &     
     &              ABSOC(MFREQC),EMISC(MFREQC),SCATC(MFREQC),            &     
     &              PLAC(MFREQC),                                         &     
     &              IJCINT(MFREQ))                                              
      ALLOCATE      (FRECR(MCROSS,MFCRA),CROSR(MCROSS,MFCRA),             &     
     &              CRMX(MCROSS),NFCR(MCROSS))                                  
      ALLOCATE      (PHOTI(MCROSS,MFREQ))                                       
      ALLOCATE      (FRECQ(MPHOT,MFCRA),QHOT(MPHOT,MFCRA),                &     
     &              AQHT(MPHOT),EQHT(MPHOT),GQHT(MPHOT),                  &     
     &              CRMY(MPHOT),NFQHT(MPHOT))                                   
                                                                                
      end subroutine alloc_synthp                                               
!                                                                               
!*****************************************************************              
!                                                                               
      subroutine alloc_wincom                                                   
      use params                                                                
      use wincom                                                                
                                                                                
      ALLOCATE      (BMU(MKU,MDEPTH),WMUJ(MKU,MDEPTH),WMUH(MKU),          &     
     &              RD(MDEPTH),PIM(MKU),RAD1(MDEPTH),                     &     
     &              DELZ(MKU,MDEPTH),NUD(MKU),NUDF(MKU),                  &     
     &              DELZF(MEXT,MDEPF ),DFRQF(MEXT,2*MDEPF ),              &     
     &              VEL(MDEPTH),DFRQ(MKU,2*MDEPTH),DVD(MDEPTH))                 
      ALLOCATE      (FFQ(MOPAC),FFQV(MOPAC),RDF(MDEPF),DENSF(MDEPF),      &     
     &              VELF(MEXT,MDEPF),DRAY(MEXT,2*MDEPF),                  &     
     &              KRAY(MEXT,2*MDEPF),                                   &     
     &              WDIL(MDEPTH),PLANW(MDEPTH),TRAD(MTRAD,MDEPTH))              
      ALLOCATE      (ILNE(MDEPTH),ILVI(MDEPTH))                                 
                                                                                
      ALLOCATE     (IHYLW(MFREQ),ILOWHW(MFREQ),                           &     
     &              M10W(MFREQ),M20W(MFREQ),                              &     
     &              IHE2LW(MFREQ),ILWHEW(MFREQ),                          &     
     &              MHE10W(MFREQ),MHE20W(MFREQ))                                

      end subroutine alloc_wincom                                               
!                                                                               
!*****************************************************************              
!                                                                               
      subroutine alloc_optabl                                                   
      use params                                                                
      use optabl                                                                
                                                                                
!     ALLOCATE (tempg(mttab),densg(mttab,mrtab),elecgr(mttab,mrtab),      &     
      ALLOCATE (tempg(mttab),elecgr(mttab,mrtab),                         &     
     &         densg0(mttab),                                             &     
     &         wlgrid(mfgrid))                                                  
      ALLOCATE (absgrd(mttab,mrtab,mfgrid))                                     
      ALLOCATE (absop(msftab),wltab(msftab))                                    
!     ALLOCATE (nden(mttab),nfrtab(mttab,mrtab))                                
      ALLOCATE (nfrtab(mttab,mrtab))                                
      ALLOCATE (elecm(mdepth))                                                  
      ALLOCATE (relabn(matom))                                                  
      ALLOCATE (abgrd(mfgrid),xli(3))                                           
      ALLOCATE (yint(mfgrid),jint(mfgrid))                                      
                                                                                
                                                                                
      end subroutine alloc_optabl                                               
!                                                                               
!*****************************************************************              
!                                                                               
                                                                                
      subroutine alloc_molist                                                   
                                                                                
      use molist                                                                
      use lindat                                                                
                                                                                
      ALLOCATE     (NLINM0(MMLINI),                                       &     
     &              NLINML(MMLINI),                                       &     
     &              NLINMT(MMLINI),                                       &     
     &              INACTM(MMLINI))                                             
                                                                                
      ALLOCATE     (FREQM(MLINM0,MMLIST),                                 &     
     &              EXCLM(MLINM0,MMLIST),                                 &     
     &              GFM(MLINM0,MMLIST),                                   &     
     &              EXTINM(MLINM0,MMLIST))                                      
                                                                                
      ALLOCATE     (INDATM(MLINM0,MMLIST),                                &     
     &              INMLIN(MLINM,MMLIST),                                 &     
     &              INMLIP(MLINM,MMLIST))                                       
                                                                                
      IF(MBROAD.GT.0) THEN                                                      
         ALLOCATE  (GRM(MLINM0,MBROAD),                                   &     
     &              GSM(MLINM0,MBROAD),                                   &     
     &              GWM(MLINM0,MBROAD))                                         
      END IF                                                                    
                                                                                
                                                                                
      IF(MVDWLI.GT.0) THEN                                                      
         ALLOCATE  (GVDWH2(MLINM0,MVDWLI),                                &     
     &              GEXPH2(MLINM0,MVDWLI),                                &     
     &              GVDWHE(MLINM0,MVDWLI),                                &     
     &              GEXPHE(MLINM0,MVDWLI))                                      
      END IF                                                                    
                                                                                
      ALLOCATE     (IJCMTR(MLINM,MMLIST),                                 &     
     &              NXTSEM(MMLIST),                                       &     
     &              IPRSEM(MMLIST),                                       &     
     &              IREADM(MMLIST),                                       &
     &              ILASTM(MMLIST))      
      ALLOCATE     (FRLASM(MMLIST),                                       &     
     &              ALASTM(MMLIST),                                       &     
     &              ALEND(MMLIST))                                              
                                                                                
      ALLOCATE     (MOLIND(11000),                                        &     
     &              IONIND(11000))                                              
                                                                                
      end subroutine alloc_molist                                               
!                                                                               
!*****************************************************************              
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
      subroutine alloc_heprf                                                    
      use params                                                                
      use heprf                                                                 
                                                                                
      ALLOCATE (PRFHE2(19,MDEPTH,36),WLHE2(19,36),                        &     
     &         NWLHE2(19),ILHE2(19),IUHE2(19))                                  
      ALLOCATE (WL2(36,19),XT2(6),XNE2(11,19),PRF2(36,6,11))                    
      ALLOCATE (PRFHE1(50,4,8,3),DLMHE1(50,8,3),XNEHE1(8),                &     
     &         NWLAM(8,4),                                                &     
     &         PRF447(80,4,7),DLM447(80,7),XNE447(7))                           
                                                                                
      end subroutine alloc_heprf                                                
                                                                                
                                                                                
!                                                                               
!*****************************************************************              
!                                                                               
                                                                                
      subroutine alloc_opadat                                                   
                                                                                
      use params                                                                
      use opadat                                                                
                                                                                
      allocate   (frac(mtemp,melec,mion1),fracm(mtemp,melec))                   
      allocate   (itemp(mtemp))                                                 
                                                                                
      allocate   (pfunc(31,30,mdepth),frpf0(31,30))                             
      allocate   (alogtt(97), alogrhot(31,30), fract(97,31,31,30),        &     
     &           zct(97,31,31,30), felect(97,31,30), alognet(97,31,30),   &     
     &           fractH2(97,31),fractH2p(97,31),fracHm(97,31),            &     
     &           zctH2(97,31),zctH2p(97,31),zctHm(97,31),                 &     
     &            zcthyd(97,31),zctHp(97,31))                                   
      allocate   (nrho(30),ntempe(30))                                          
                                                                                
      end subroutine alloc_opadat                                               
                                                                                
!                                                                               
!*****************************************************************              
!                                                                               
                                                                                
      subroutine alloc_photcs                                                   
                                                                                
      use params                                                                
      use photcs                                                                
                                                                                
      ALLOCATE   (PHOT(MFRQ,MPHOT),                                       &     
     &           APHT(MPHOT),EPHT(MPHOT),GPHT(MPHOT))                           
      ALLOCATE   (JPHT(MPHOT))                                                  
                                                                                
      end subroutine alloc_photcs                                               
!                                                                               
!*****************************************************************              
!                                                                               
                                                                                
      subroutine alloc_hydprf                                                   
                                                                                
      use params                                                                
      use hydprf                                                                
                                                                                
      allocate   (tabden(mtab,mden),tabtem(mtab,mtem,mden),               &     
     &           tablam(mtab,mlam,mden),prftab(mtab,mlam,mtem,mden))            
      allocate   (nden(mtab),numtem(mtab,mden),numlam(mtab,mden))               
                                                                                
      allocate   (frgtab(mfhtab),wlgtab(mfhtab),hydopg(mfhtab,mdepth))          
                                                                                
      end subroutine alloc_hydprf
!
!*****************************************************************
!
                                
       subroutine alloc_topdat

       use topdat

       allocate   (sop(mop,mmaxop),xop(mop,mmaxop))
       allocate   (nop(mmaxop))
       allocate   (idlvop(mmaxop))

       end subroutine alloc_topdat
!
!*****************************************************************
!

      subroutine alloc_hydxen

      use params
      use hydxen

      ALLOCATE      (PRFXB(MLINH,MHWL,MHT,MHE),                           &
     &              PRFXR(MLINH,MHWL,MHT,MHE),                            &
     &              PRFB(MLINH,MDEPTH,MHWL),                              &
     &              PRFR(MLINH,MDEPTH,MHWL),                              &
     &              ALXEN(MLINH,MHWL),                                    &
     &              XTXEN(MHT,MLINH),                                     &
     &              XNEXEN(MHE,MLINH),                                    &
     &              NTHXEN(MLINH),                                        &
     &              NEHXEN(MLINH))                                        &
   
      end subroutine alloc_hydxen
!
!*****************************************************************
!

      subroutine alloc_allarn

      use params
      use allarn

      allocate     (xlalp(nxmax),plalp(nxmax,nnmax),                      &
     &              xlbet(nxmax),plbet(nxmax,nnmax),                      &
     &              xlgam(nxmax),plgam(nxmax,nnmax),                      &
     &              xlbal(nxmax),plbal(nxmax,nnmax))

      end subroutine alloc_allarn
