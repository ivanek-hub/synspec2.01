      program redulis
!     ===============
!
!     reducing a molecular lione list
!
!     INPUT:  fort.20 - original list
!             std.input: mandataory:   gfelim - rejection parameter
!                        optional:     tstd - characteristic T [K]
!                                      wlchar - characteristic wavlength [A]
!                                      rnu0 - characteristic log(N/U)_char
!                                             see Eq. (54) of Paper IV
!            if the least three are not set, 
!            then TSTD=3000, WLCHAR=10000., RNU0=10.
!
!     OUTPUT: fort.12 - reduced line list
!
      implicit real*8(a-h,o-z)
      character(len=4)  :: alab(20)
      character(len=80) :: rec
      dimension rnu(20),cod(20),insm(20),nlito(20),nliin(20)
      data rnu/16.2,10.0,10.8,7.9,11.2,6.5,5.1,8.3,7.5,9.6,
     *          6.1, 8.1, 6.6,8.4, 6.7,2.7,5.1,4.3,5.8,5.4/
      data nlito/20*0/
      data nliin/20*0/
      data alab /'H2  ','H2O ','OH  ','CH  ','CO  ','CN  ',
     *           'C2  ','NH  ','SiH ','SiO ','TiO ','VO  ',
     *           'MgH ','CaH ','NaH ','MgO ','AlO ','CaO ',
     *           'CrH ','FeH '/
      data insm/2,3,4,5,6,7,8,12,17,25,29,30,32,34,122,126,134,
     *          179,198,214/

      iunit=5
      gfelim=-9.
      tstd=3000.
      wlchar=10000.
      rnu0=10.
      read(iunit,'(a80)') rec
      read(rec,*,iostat=kst) gfelim,tstd,wlchar,rnu0
      if(kst.ne.0) then
         read(rec,*,iostat=ios) gfel
         if(ios.eq.0) gflim=gfel
      end if
      wlchar=wlchar*0.1
!
      readmol: do
         ts=tstd
         read(20,'(a80)',iostat=ios2) rec
         if(ios2.ne.0) exit readmol
         read(rec,*) alam,code,gf,el
         if(code.eq.101.00) then
            inm=1
          else if(code.eq.10108.00) then
            inm=2
          else if(code.eq.108.00) then
            inm=3
          else if(code.eq.106.00) then
            inm=4
            ts=5000.
          else if(code.eq.608.00) then
            inm=5
          else if(code.eq.607.00) then
            inm=6
            ts=5000.
          else if(code.eq.606.00) then
            inm=7
          else if(code.eq.107.00) then
            inm=8
            ts=4000.
          else if(code.eq.114.00) then
            inm=9
            ts=4000
             else if(code.eq.814.00) then
            inm=10
          else if(code.eq.822.00) then
            inm=11
          else if(code.eq.823.00) then
            inm=12
          else if(code.eq.112.00) then
            inm=13
          else if(code.eq.120.00) then
            inm=14
          else if(code.eq.111.00) then
            inm=15
          else if(code.eq.812.00) then
            inm=16
          else if(code.eq.813.00) then
            inm=17
          else if(code.eq.820.00) then
            inm=18
          else if(code.eq.124.00) then
            inm=19
          else if(code.eq.126.00) then
            inm=20
         end if
         nlito(inm)=nlito(inm)+1
         cod(inm)=code
      
         gfe=gf-0.625*el/ts+log10(alam/wlchar)+rnu(inm)-rnu0
         if(gfe.gt.gfelim) then
            nliin(inm)=nliin(inm)+1 
            write(12,"(a80)") rec
         end if
!
      end do readmol
!
      write(6,"('gfelim,tstd,wlchar,rnu0 = ',4f10.1/)")                    &
     &  gfelim,tstd,wlchar,rnu0
!
      nto=0
      nin=0
      do inm=1,20
         nto=nto+nlito(inm)
         nin=nin+nliin(inm)
!        if(inm.eq.7) write(*,*) inm,cod(inm),alab(inm)
         write(6,"(i4,f10.2,2x,a4,2i10)")                                 & 
     &      inm,cod(inm),alab(inm),nlito(inm),nliin(inm)
      end do
      write(*,*)
      write(*,*) 'total no',nto
      write(*,*) 'included',nin
!
      end     program redulis
