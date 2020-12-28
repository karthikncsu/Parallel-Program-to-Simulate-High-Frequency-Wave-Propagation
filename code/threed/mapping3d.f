       subroutine mapping3d(elem,nelem,locelem,locnelem,mapold2new
     &   ,N,locn)
! subroutine which maps old node to new nodes
       implicit NONE
       INTEGER N,locnelem,locn,nelem,i
       INTEGER elem(nelem,8),locelem(locnelem),mapold2new(N)
       INTEGER,dimension (:), allocatable :: locnod
       INTEGER flag,locnnod
       
       allocate(locnod(1:8*locnelem))

       do i=1,locnelem
       locnod((i-1)*8+1)=elem(locelem(i),1)
       locnod((i-1)*8+2)=elem(locelem(i),2)
       locnod((i-1)*8+3)=elem(locelem(i),3)
       locnod((i-1)*8+4)=elem(locelem(i),4)
       locnod((i-1)*8+5)=elem(locelem(i),5)
       locnod((i-1)*8+6)=elem(locelem(i),6)
       locnod((i-1)*8+7)=elem(locelem(i),7)
       locnod((i-1)*8+8)=elem(locelem(i),8)
       enddo
       
       call intrmvdup(locnod,8*locnelem,locn)
       mapold2new(:)=0

!       write(*,*) "done with intrmbdup"       
!       write(*,*) locnod
       do i=1,locn
       mapold2new(locnod(i))=i
       enddo
!       write(*,*) "mapold2new"
!       write(*,*) mapold2new
      
       end
