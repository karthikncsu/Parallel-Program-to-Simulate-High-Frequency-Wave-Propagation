       subroutine newelem(elem,nelem,lelem,lnelem,mapold2new,N,dumelem)
       implicit NONE
       INTEGER N,nelem,lnelem
       INTEGER elem(nelem,4),lelem(lnelem),mapold2new(N)
       INTEGER dumelem(lnelem,4)
       
       do i=1,lnelem
       dumelem(i,:)=elem(lelem(i),:)
       enddo

      
       end
