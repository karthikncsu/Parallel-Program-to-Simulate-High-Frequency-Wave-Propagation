       subroutine intsort(A,M)
        implicit NONE
       
        INTEGER M,A(M),i,j,swap

        do i=M-1,1,-1
        do j=1,i
            if(A(j)>A(j+1)) then
            swap=A(j)
            A(j)=A(j+1)
            A(j+1)=swap
            endif
        enddo
        enddo
      
        end
