       subroutine intrmvdup(A,M,bndcnt)
! This function sorts the given array in increasing order removes the
! duplicate elements
        implicit NONE
       
        INTEGER M,A(M),i,bndcnt
        INTEGER,dimension (:), allocatable :: Adum
        
        allocate(Adum(1:M))
        Adum(:)=0
        call intsort(A,M)
!        write(*,*) "intsort is good"
        
        bndcnt=1
        Adum(bndcnt)=A(1)

        do i=2,M
        if(A(i).ne.A(i-1)) then
        bndcnt=bndcnt+1
        Adum(bndcnt)=A(i)
        endif
        enddo

        A=Adum        
        deallocate(Adum)
      
        end
