      program mergelist
 
      implicit none 

      real,allocatable ::  alam(:),anum(:),gf(:),excl(:)!,gr(:),gs(:),gw(:)     
      integer, parameter :: nl=50000000
      real    :: al,anu,a,ex
      integer :: i,n,ios,j,k,inewend,ioldend

      allocate (alam(nl),anum(nl),gf(nl),excl(nl))!,gr(nl),gs(nl),gw(nl))
 
      i=1
      do
         read(5,*,iostat=ios) alam(i),anum(i),gf(i),excl(i)
!          anum(i)=0.; gf(i)=0.; excl(i)=0.
!          write(*,*) 'i,al,ios',i,alam(i),ios
         if(ios.ne.0) exit
         i=i+1
      end do
      n=i-1
      write(*,*) 'number of lines in the first list ',n
!
      j=0
      k=1
      inewend=0
      ioldend=0

      anu=0.; a=0.; ex=0.

      read1: do
         read(10,*,iostat=ios) al,anu,a,ex
         if(ios.ne.0) exit read1
         if(ioldend.eq.1) then
            write(11,"(f10.4,f10.2,f8.3,f12.3)") al,anu,a,ex
            j=j+1
            cycle read1
         end if
         if(al.le.alam(k)) then
            write(11,"(f10.4,f10.2,f8.3,f12.3)") al,anu,a,ex
            j=j+1
          else
            mer1: do
            if(alam(k).lt.al) then 
!              write(6,"(2i10,2f15.4,3i10)")                              &
!    &         k,n,alam(k),al,j,ioldend,inewend
               if(k.gt.n) then
                  ioldend=1
                  if(inewend.eq.0) then
                    write(11,"(f10.4,f10.2,f8.3,f12.3)") al,anu,a,ex
                    j=j+1
                    cycle read1
                  end if
                else
                  write(11,"(f10.4,f10.2,f8.3,f12.3)")                    &
     &            alam(k),anum(k),gf(k),excl(k)
                  k=k+1
                  cycle mer1
               end if
             else
               write(11,"(f10.4,f10.2,f8.3,f12.3)") al,anu,a,ex
               j=j+1
               cycle read1
            end if
            end do mer1
         end if
      end do read1

      inewend=1
      mer2: do
         if(k.gt.n) exit mer2
         write(11,"(f10.4,f10.2,f8.3,f12.3)")                             &
     &   alam(k),anum(k),gf(k),excl(k)
         k=k+1
      end do mer2
      ioldend=1

      write(*,*) 'number of lines in the second list',j
 
      close(11)
      end program mergelist
 
