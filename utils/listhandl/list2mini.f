      program mini
c
      implicit none

      character(len=80) :: dum80
      character(len=40) :: dum40
      integer :: ios
c
      do
         read(5,"(a80)",iostat=ios) dum80
         if(ios.ne.0) exit
         dum40=dum80(1:40)
         write(6,"(a40)") dum40
      end do
      end program mini
