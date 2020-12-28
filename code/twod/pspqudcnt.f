        subroutine pspqudcnt(dumpo,elem,nelem,nnod,pocount,pomaxcnt)
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
        INTEGER dumpo(nnod,9),elem(nelem,4),nnod,nelem,pocount(nnod)
        INTEGER pomaxcnt
        
! Variables of the subroutine
        integer ipnt1,ipnt2,ipnt3,ipnt4,ielem,inod,ipnt

! Finding the nodes surrounding a node        
        dumpo(:,:)=0

        do ielem=1,nelem
        
        ipnt1=elem(ielem,1)
        ipnt2=elem(ielem,2)
        ipnt3=elem(ielem,3)
        ipnt4=elem(ielem,4)

        dumpo(ipnt1,1)=ipnt2
        dumpo(ipnt1,2)=ipnt3
        dumpo(ipnt1,3)=ipnt4

        dumpo(ipnt2,3)=ipnt3
        dumpo(ipnt2,4)=ipnt4
        dumpo(ipnt2,5)=ipnt1

        dumpo(ipnt3,5)=ipnt4
        dumpo(ipnt3,6)=ipnt1
        dumpo(ipnt3,7)=ipnt2

        dumpo(ipnt4,1)=ipnt3
        dumpo(ipnt4,7)=ipnt1
        dumpo(ipnt4,8)=ipnt2
        enddo
       
        ! Counting the number of nodes surrounding each node
        pomaxcnt=0
        pocount(:)=0
        do inod=1,nnod
            dumpo(inod,9)=inod
            do ipnt=1,9
            if(dumpo(inod,ipnt).ne.0) then
            pocount(inod)=pocount(inod)+1
            endif
            enddo
        pomaxcnt=pomaxcnt+pocount(inod)
        enddo

      
        end
