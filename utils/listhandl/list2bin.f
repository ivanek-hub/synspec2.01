      program list2bin
!     ================
!
!     transforms an ASCII line list to a binary
!
!     INPUT:  std input - ASCII list
!
!     OUTPUT: fort.12 - binary list
!             std. output gives an info on number of lines
!
      implicit none

      character(len=102) :: rec
      real   ::  x(10)
      integer :: j(7)
      integer :: n,npr,npi,i,kst,kst1,kst2,kst3,kst4,ios
!
!     first, find the format of the line list by examining the 1st record
!
      n=1
      npr=10
      npi=7
      read(5,'(a102)') rec
      write(*,*) rec
      read(rec,*,iostat=kst) (x(i),i=1,10),(j(i),i=1,7)
      if(kst.ne.0) then
         read(rec,*,iostat=kst1) (x(i),i=1,10),j(1)
         npr=10
         npi=1
         if(kst1.ne.0) then
            read(rec,*,iostat=kst2) (x(i),i=1,9)
            npr=9
            npi=0
            if(kst2.ne.0) then
               read(rec,*,iostat=kst3) (x(i),i=1,7)
               npr=7
               npi=0
               if(kst3.ne.0) then
                  read(rec,*,iostat=kst4) (x(i),i=1,4)
                  npr=4
                  npi=0
                  if(kst4.ne.0) then
                     write(*,*) 
     *               'no applicable format of the line list detected'
                     stop
                  end if
               end if
            end if
         end if
      end if
!
      write(6,"('number of parameters in the list; real,int:',2i4)")      &
     &   npr,npi
!
!     read the list and transport to binary
!
      iolist: do
         read(5,*,iostat=ios) (x(i),i=1,npr),(j(i),i=1,npi)
         if(ios.gt.0) cycle iolist
         if(ios.le.0) exit iolist
         n=n+1
         write(12) (x(i),i=1,npr),(j(i),i=1,npi)
      end do iolist
!
      write(6,"(i10,'  lines included in the binary file fort.12')") n
 
      end program list2bin

