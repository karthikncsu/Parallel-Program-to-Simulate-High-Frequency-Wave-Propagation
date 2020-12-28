        subroutine dotproduct(x,y,M,res)
        INTEGER M
        DOUBLE PRECISION x(M),y(M),res
        integer i
        res=0.d0
        do i=1,M
            res=res+x(i)*y(i)
        enddo
        
        end
