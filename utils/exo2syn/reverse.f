      program reverse
 
      character(len=80),allocatable :: a(:)
      allocate (a(40000000000))
 
      i=0
      readloop: do
         i=i+1
         read(5,"(a80)",iostat=ios) a(i)
         if(ios.ne.0) exit readloop
      end do readloop
      n=i-1
 
      do i=n,1,-1
         write(6,"(a80)") a(i)
      end do
      close(6)
 
      end program reverse
