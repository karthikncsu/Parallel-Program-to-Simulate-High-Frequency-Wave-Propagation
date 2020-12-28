       subroutine pspqudll(posupo,M,posuloc,nnod,dumpo)
        implicit NONE
!-------------------------------------------------------------------
! This subroutine saves (nnod x 8 array)dumpo array into linked list
! The ending and starting of the linked list is stored in psuloc.
!
! elem:     elements array          
! nelem:    Number of elements
! nnod:     Number of nodes
!-------------------------------------------------------------------
        
! Subroutine Variables        
        INTEGER dumpo(nnod,9),nnod
        INTEGER M,posupo(M),posuloc(nnod+1)
        
! Variables of the subroutine
        integer pnt1,pnt2,inod,ielem,flag,icount,N
       
        flag=1
        posuloc(1)=1
        do inod=1,nnod
            do icount=1,9
            if(dumpo(inod,icount).ne.0) then
                posupo(flag)=dumpo(inod,icount)
                flag=flag+1
            endif
            enddo
        posuloc(inod+1)=flag

        pnt1=posuloc(inod)
        pnt2=posuloc(inod+1)-1
        N=posuloc(inod+1)-posuloc(inod)
        call intsort(posupo(pnt1:pnt2),N)

        enddo
        
        
        end
