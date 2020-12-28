       subroutine mapping(elem,nelem,locelem,locnelem,mapold2new,N,locn)
! subroutine which maps old node to new nodes
       implicit NONE
       INTEGER N,locnelem,locn,nelem,i
       INTEGER elem(nelem,4),locelem(locnelem),mapold2new(N)
       INTEGER,dimension (:), allocatable :: locnod
       INTEGER flag,locnnod
       
       allocate(locnod(1:4*locnelem))

       do i=1,locnelem
       locnod((i-1)*4+1)=elem(locelem(i),1)
       locnod((i-1)*4+2)=elem(locelem(i),2)
       locnod((i-1)*4+3)=elem(locelem(i),3)
       locnod((i-1)*4+4)=elem(locelem(i),4)
       enddo
       
       call intrmvdup(locnod,4*locnelem,locn)
       mapold2new(:)=0
       
!       write(*,*) locnod
       do i=1,locn
       mapold2new(locnod(i))=i
       enddo
!       write(*,*) "mapold2new"
!       write(*,*) mapold2new
      
       end
