        subroutine psphexcnt(dumpo,elem,nelem,nnod,pocount,pomaxcnt)
        implicit NONE
!-------------------------------------------------------------------
! This subroutine finds the points surrounding points and saves in 
! nnod x 8 array(dumpo). It also finds the number of points surro-
! -nding each node and saves in pocount(nnod) array. pomaxcnt is
! the sum of the array pocount.
!
! elem:     elements array          
! nelem:    Number of elements
! nnod:     Number of nodes
!-------------------------------------------------------------------
! Subroutine Variables        
        INTEGER dumpo(nnod,27),elem(nelem,8),nnod,nelem,pocount(nnod)
        INTEGER pomaxcnt
        
! Variables of the subroutine
        integer ipnt1,ipnt2,ipnt3,ipnt4,ipnt5,ipnt6,ipnt7,ipnt8
        integer ielem,inod,ipnt

! Finding the nodes surrounding a node        
        dumpo(:,:)=0

        do ielem=1,nelem
        
        ipnt1=elem(ielem,1)
        ipnt2=elem(ielem,2)
        ipnt3=elem(ielem,3)
        ipnt4=elem(ielem,4)
        ipnt5=elem(ielem,5)
        ipnt6=elem(ielem,6)
        ipnt7=elem(ielem,7)
        ipnt8=elem(ielem,8)

        dumpo(ipnt1,1)=ipnt1
        dumpo(ipnt1,2)=ipnt2
        dumpo(ipnt1,3)=ipnt3
        dumpo(ipnt1,4)=ipnt4
        dumpo(ipnt1,1+9)=ipnt5
        dumpo(ipnt1,2+9)=ipnt6
        dumpo(ipnt1,3+9)=ipnt7
        dumpo(ipnt1,4+9)=ipnt8
        
        dumpo(ipnt2,1)=ipnt2
        dumpo(ipnt2,4)=ipnt3
        dumpo(ipnt2,5)=ipnt4
        dumpo(ipnt2,6)=ipnt1
        dumpo(ipnt2,1+9)=ipnt6
        dumpo(ipnt2,4+9)=ipnt7
        dumpo(ipnt2,5+9)=ipnt8
        dumpo(ipnt2,6+9)=ipnt5
        
        dumpo(ipnt3,1)=ipnt3
        dumpo(ipnt3,6)=ipnt4
        dumpo(ipnt3,7)=ipnt1
        dumpo(ipnt3,8)=ipnt2
        dumpo(ipnt3,1+9)=ipnt7
        dumpo(ipnt3,6+9)=ipnt8
        dumpo(ipnt3,7+9)=ipnt5
        dumpo(ipnt3,8+9)=ipnt6

        dumpo(ipnt4,1)=ipnt4
        dumpo(ipnt4,2)=ipnt3
        dumpo(ipnt4,8)=ipnt1
        dumpo(ipnt4,9)=ipnt2
        dumpo(ipnt4,1+9)=ipnt8
        dumpo(ipnt4,2+9)=ipnt7
        dumpo(ipnt4,8+9)=ipnt5
        dumpo(ipnt4,9+9)=ipnt6
        
        dumpo(ipnt5,1)=ipnt5
        dumpo(ipnt5,2)=ipnt6
        dumpo(ipnt5,3)=ipnt7
        dumpo(ipnt5,4)=ipnt8
        dumpo(ipnt5,1+9+9)=ipnt1
        dumpo(ipnt5,2+9+9)=ipnt2
        dumpo(ipnt5,3+9+9)=ipnt3
        dumpo(ipnt5,4+9+9)=ipnt4

        dumpo(ipnt6,1)=ipnt6
        dumpo(ipnt6,4)=ipnt7
        dumpo(ipnt6,5)=ipnt8
        dumpo(ipnt6,6)=ipnt5
        dumpo(ipnt6,1+9+9)=ipnt2
        dumpo(ipnt6,4+9+9)=ipnt3
        dumpo(ipnt6,5+9+9)=ipnt4
        dumpo(ipnt6,6+9+9)=ipnt1
        
        dumpo(ipnt7,1)=ipnt7
        dumpo(ipnt7,6)=ipnt8
        dumpo(ipnt7,7)=ipnt5
        dumpo(ipnt7,8)=ipnt6
        dumpo(ipnt7,1+9+9)=ipnt3
        dumpo(ipnt7,6+9+9)=ipnt4
        dumpo(ipnt7,7+9+9)=ipnt1
        dumpo(ipnt7,8+9+9)=ipnt2

        dumpo(ipnt8,1)=ipnt8
        dumpo(ipnt8,2)=ipnt7
        dumpo(ipnt8,8)=ipnt5
        dumpo(ipnt8,9)=ipnt6
        dumpo(ipnt8,1+9+9)=ipnt4
        dumpo(ipnt8,2+9+9)=ipnt3
        dumpo(ipnt8,8+9+9)=ipnt1
        dumpo(ipnt8,9+9+9)=ipnt2

        enddo
       
        ! Counting the number of nodes surrounding each node
        pomaxcnt=0
        pocount(:)=0
        do inod=1,nnod
            do ipnt=1,27
            if(dumpo(inod,ipnt).ne.0) then
            pocount(inod)=pocount(inod)+1
            endif
            enddo
        pomaxcnt=pomaxcnt+pocount(inod)
        enddo

      
        end
