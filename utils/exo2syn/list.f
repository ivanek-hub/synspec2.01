      program list
!     ============
!
!     program to extract line data for a particular isotopologue
!     of a molecule from the EXOMOL data and output them in the
!     Synspec format
!
!     Standard Input: 
!      type   - a Kurucz-type label of a molecule, e.g 813.00 for AlO
!      molname - a character variable with the chemical label, e.g. 'AlO'
!      iat1,iso1 - atomic number and isotope of the first atom.
!                  e.g  13, 27 for 27Al
!      iat2, iso2 - the same for the second atom, e.g.  8  16 for 16O
!      core - character variable with the core name of the EXOMOL files,
!             e.g. '27Al-16O__ATP'
!      wbeg - lowest wavelength [nm] for a line to be included in the list
!      wend - the same for the longest wavelength
!      gflim - log(gf) limit- lines with log(gf) below gflim are excluded
!     
!     Input files:
!      - the program accumes the the three EXOMOL files are present 
!        in the current directory, e.g.
!        27Al-16O__ATP.states   - energy states
!        27Al-16O__ATP.trans    - transitions
!        27Al-16O__ATP.pf       - partition functions
!      - file "isotops" - solar isotoic fractions of the individual atoms
!      - file with the name [molname].pf, e.g. AlO.pf, located in the
!        directory ../../EXOMOL_pfs (may be changed to suit the user's needs)  
!
!     Output:
!      - standard ouput - a short log file
!      - file fort.10 - the resulting line list in the Synspec format,
!        but in reverse wavelength order  - from the long to short.
!        An accompanied short program "reverse.f" reverses it to the
!        ordering needed by Synspec - short to long wavelengths
!     
!        Notice that wavelengths above 2000A are transformed from the
!        vacuum to air wavelengths.
!        an accompanied short program "reverse.f" reverses it to the
!        ordering needed by Synspec - short to long wavelengths
!
!
      implicit none

      integer, parameter :: dp=selected_real_kind(15,300)
      integer, parameter :: mt=10000,mst=1000000,mato=26,miso=58
      integer            :: ml 

      real(dp), allocatable :: al(:),gl(:),el(:)

      character(len=50) :: core,molname,sfile,tfile,pfile,pafile
      character(len=80) :: cadena
      real(dp) :: fra(mato,miso)
      real(dp) :: e(mst),g(mst)

      real(dp) :: con,wbeg,wend,gflim,fr,sfr,wtbeg,wtend,                 &
     &            t2,pf2,ta2,pfa2,relu,relu1,relu2,fis,fiso,a,wtr,        &
     &            wlam,wlam0,vac,wl0,alm,xn1,gf
      real(dp) :: type
      integer  :: ia,is,iat1,iat2,iso1,iso2,ios,iatm,isot,is1,is2,        &
     &            i,im,n,ninc,nu,nl,nw

!
      con=1.e-14/8./3.1415926/0.02654
!
      do ia=1,mato
         do is=1,miso
            fra(ia,is)=0.
         end do
      end do
!
      read(5,*) type,molname
      read(5,*) iat1,iso1,iat2,iso2
      read(5,*) core
      read(5,*) wbeg,wend,gflim
      write(*,*) 'species:       ',molname
      write(*,*) 'EXOMOL files  :',core
      write(6,"(' IAT1,ISO1,IAT2,ISO2',2i4,2x,2i4)") iat1,iso1,iat2,iso2
      write(6,"(' Wbeg,Wend,gflim    ',3f11.2)") wbeg,wend,gflim
!
!     determine the approprite isotopic fraction of the given isotopologue
!
      open(unit=7,file='isotops',status='old')
      readfrac: do
         read(7,*,iostat=ios) iatm,isot,fr
         if(ios.ne.0) exit
         fra(iatm,isot)=fr*0.01
      end do readfrac
      close(7)
!
!     normalize the isotopic fraction
!
      sfr=0.
      do is1=1,miso
         do is2=1,miso
            sfr=sfr+fra(iat1,is1)*fra(iat2,is2)
         end do
      end do
      write(*,*) 'frtot',sfr
!
      fiso=fra(iat1,iso1)*fra(iat2,iso2)
!
      wtbeg=1.d7/max(wbeg,100.)
      wtend=1.d7/wend
!
      do i=len(core),1,-1
         if(core(i:i).ne.' ') exit
      end do

      sfile= core(1:i) // '.states'
      tfile= core(1:i) // '.trans'
      pfile= core(1:i) // '.pf'
      write(*,*) 'sfile  ',sfile
      write(*,*) 'tfile  ',tfile
      write(*,*) 'pfile  ',pfile
!
!     read file with partition functions for the given isotopologue
!
      open(unit=7,file=pfile,status='old')
      do
         read(7,*) t2,pf2
         if(t2.gt.999.9.and.t2.lt.1000.1) exit
      end do
!
!     read avereaged partition function
!
      do i=len(molname),1,-1
         if(molname(i:i).ne.' ') exit
      end do

      im=i
      pafile= molname(1:im) // '.pf'
      open(unit=7,file=pafile,status='old')
      do
         read(7,*) ta2,pfa2
         if(ta2.gt.999.9.and.ta2.lt.1000.1) exit
      end do
      close(7)
!
!     determine ratio of the partition functions at 1000 K,
!     which is constant with T
!
      relu=1.
      relu1=1.
      relu2=1.
      relu2=pf2/pfa2
      relu=relu2
      fis=fiso
      fiso=log10(fiso/relu)
!
      write(6,"(' T = ',0pf10.1,1pe11.3,0pf10.1,1pe11.3,0pf10.3/)")       & 
     &   t2,pf2,ta2,pfa2,relu2
      write(6,"(' isot.fract. ',f11.4/                                    &
     &       ' U(iso)/U(av)',f11.4/                                       &
     &       ' log(cor.fac)',f11.4)") fis,relu,fiso
!           
      open(unit=11,file=sfile,status='old')
      open(unit=12,file=tfile,status='old')
!
!     read data for energy states
!
      i=0
      do
         i=i+1
         read(11,*,iostat=ios) i,e(i),g(i)
         if(ios.ne.0) exit
      end do
      i=i-1
      write(*,*) 'total number of states',i
      close(11)
!
!     read data for transitions and write line data in the Synspec format
!
      read(12,"(a80)") cadena
      call count_words(cadena,nw)
      backspace(12)
      
      write(*,"(/'number of columns in the .trans file:',i4/)") nw
      if(nw.eq.4) then

      n=0
      ninc=0
      readtr: do
         n=n+1
         read(12,*,iostat=ios) nu,nl,a,wtr
         if(ios.ne.0) exit readtr
         if(e(nu)-e(nl).le.0.) cycle readtr
         if(e(nu).gt.0..and.e(nu).lt.2.e5.and.                            &
     &      e(nl).gt.0..and.e(nl).lt.2.e5) then
            wlam=1.e7_dp/wtr
!           wlam=1.e7_dp/(e(nu)-e(nl))
            if(wlam.ge.wend*0.999999) cycle readtr
            if(wlam.lt.wbeg) cycle readtr
!
!           transform to air wavelengths above 2000 A
!
            vac=200.
            wlam0=wlam
            if(wlam.gt.vac) then
               wl0=wlam*10.
               ALM=1.E8/(WL0*WL0)
               XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM)
               WL0=WL0/(XN1*1.E-6+1.)
               wlam=wl0*0.1
            END IF
!
            gf=log10(a*wlam0*wlam0*g(nu)*con)+fiso
            if(gf.lt.gflim) cycle readtr
            ninc=ninc+1
            write(10,"(f10.4,f10.2,f8.3,f12.3)") wlam,type,gf,e(nl)
         end if
      end do readtr

      else if(nw.eq.3) then

      ninc=0
      readtra: do
         read(12,*,iostat=ios) nu,nl,a
         if(ios.ne.0) exit readtra
         if(e(nu)-e(nl).le.0.) cycle readtra
         if(e(nu).gt.0..and.e(nu).lt.2.e5.and.                            &
     &      e(nl).gt.0..and.e(nl).lt.2.e5) then
            wlam=1.e7_dp/(e(nu)-e(nl))
            if(wlam.ge.wend*0.999999) cycle readtra
            if(wlam.lt.wbeg) cycle readtra
            gf=log10(a*wlam*wlam*g(nu)*con)+fiso
            if(gf.lt.gflim) cycle readtra
            ninc=ninc+1
         end if
      end do readtra
      rewind(12)

      ml=ninc+10
      write(*,*) 'ninc,ml',ninc,ml

      allocate (al(ml),gl(ml),el(ml))

      n=0
      ninc=0
      readtr3: do
         n=n+1
         read(12,*,iostat=ios) nu,nl,a
         if(ios.ne.0) exit readtr3
         if(e(nu)-e(nl).le.0.) cycle readtr3
         if(e(nu).gt.0..and.e(nu).lt.2.e5.and.                            &
     &      e(nl).gt.0..and.e(nl).lt.2.e5) then
            wlam=1.e7_dp/(e(nu)-e(nl))
            if(wlam.ge.wend*0.999999) cycle readtr3
            if(wlam.lt.wbeg) cycle readtr3
!
!           transform to air wavelengths above 2000 A
!
            vac=200.
            wlam0=wlam
            if(wlam.gt.vac) then
               wl0=wlam*10.
               ALM=1.E8/(WL0*WL0)
               XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM)
               WL0=WL0/(XN1*1.E-6+1.)
               wlam=wl0*0.1
            END IF
!
            gf=log10(a*wlam0*wlam0*g(nu)*con)+fiso
            if(gf.lt.gflim) cycle readtr3
            ninc=ninc+1
            al(ninc)=real(wlam)
            gl(ninc)=real(gf)
            el(ninc)=real(e(nl))
            write(11,"(f10.4,f10.2,f8.3,f12.3)") wlam,type,gf,e(nl)
         end if
      end do readtr3

      call sorl(al,gl,el,type,ml,ninc)

      end if

      close(12)
      close(10)
!
      write(*,*) 'read',n,' transitions; included',ninc
      close(6)
      end
!                                        
!           
!    ***************************************************************
!
!                                                                  
        subroutine count_words(cadena,n)                           
!                                                                      
!       Counts the number of words separated by blanks in a string
!
        implicit none
            
        character(len=80)   :: cadena
        character(len=1)    :: a,b

        integer :: n,i
                                 
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

      subroutine sorl(al,gl,el,type,ml,nline)
!     =======================================
!
!     sorting lines by increasing wavelength

      implicit none

      integer, parameter :: dp=selected_real_kind(15,300)
      real(dp), intent(inout)  :: al(ml),gl(ml),el(ml)
      real(dp), intent(in)    :: type
      integer, intent(in) :: ml,nline
      integer             :: iline,j
      real(dp)            :: ali

      alllines: do iline=2,nline
         j=iline
         ali=al(iline)
         lower: do
            if(ali.gt.al(j-1)) then
               al(j)=al(j-1)
               al(j-1)=ali
               j=j-1
               if(j.le.1) exit lower
               cycle lower
             else
               exit lower
            end if
         end do lower
      end do alllines

      do j=1,nline
         write(10,"(f10.4,f10.2,f8.3,f12.3)") al(j),type,gl(j),el(j)
      end do

      end subroutine sorl
