      subroutine invmat3(Ja,invJa,detJa)

        implicit NONE
        DOUBLE PRECISION Ja(3,3),invJa(3,3),detJa
        
        detJa=0
        detJa=Ja(1,1)*(Ja(2,2)*Ja(3,3)-Ja(2,3)*Ja(3,2))
        detJa=detJa-Ja(1,2)*(Ja(2,1)*Ja(3,3)-Ja(3,1)*Ja(2,3))
        detJa=detJa+Ja(1,3)*(Ja(2,1)*Ja(3,2)-Ja(2,2)*Ja(3,1))
        
        if(detJa==0) then
        write(*,*) "Jacobi is singular"
        call exit(1)
        endif
        
        invJa(1,1)=(Ja(2,2)*Ja(3,3)-Ja(2,3)*Ja(3,2))/detJa
        invJa(2,1)=-(Ja(2,1)*Ja(3,3)-Ja(3,1)*Ja(2,3))/detJa
        invJa(3,1)=(Ja(2,1)*Ja(3,2)-Ja(2,2)*Ja(3,1))/detJa
        
        invJa(1,2)=-(Ja(1,2)*Ja(3,3)-Ja(1,3)*Ja(3,2))/detJa
        invJa(2,2)=(Ja(1,1)*Ja(3,3)-Ja(1,3)*Ja(3,1))/detJa
        invJa(3,2)=-(Ja(1,1)*Ja(3,2)-Ja(1,2)*Ja(3,1))/detJa
        
        invJa(1,3)=(Ja(1,2)*Ja(2,3)-Ja(2,2)*Ja(1,3))/detJa
        invJa(2,3)=-(Ja(1,1)*Ja(2,3)-Ja(2,1)*Ja(1,3))/detJa
        invJa(3,3)=(Ja(1,1)*Ja(2,2)-Ja(1,2)*Ja(2,1))/detJa
        

        
        
        end
