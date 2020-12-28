        subroutine findpos(JAval,JApos,A,M)
        
        INTEGER M
        INTEGER JAval,JApos,A(M),i
        
        Japos=0        
        do i=1,M
            if(A(i)==JAval) then
            JApos=i
            return
            endif
        enddo
        
        write(*,*) "Did not find the Japos value in findpos subroutine"
        call exit(1)
       
        end
