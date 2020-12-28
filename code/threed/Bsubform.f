      subroutine Bsubform(dN,Bsub,invJa)
        implicit NONE
        
        DOUBLE PRECISION dN(3),Bsub(6,3),invJa(3,3)
        DOUBLE PRECISION dNdx,DNdy,Dndz
        
        dNdx=invJa(1,1)*dN(1)+invJa(1,2)*dN(2)+invJa(1,3)*dN(3)
        dNdy=invJa(2,1)*dN(1)+invJa(2,2)*dN(2)+invJa(2,3)*dN(3)
        dNdz=invJa(3,1)*dN(1)+invJa(3,2)*dN(2)+invJa(3,3)*dN(3)

        Bsub(:,:)=0.d0

        Bsub(1,1)=dNdx;   Bsub(2,2)=dNdy;    Bsub(3,3)=dNdz
        Bsub(4,2)=dNdz;   Bsub(4,3)=dNdy
        Bsub(5,1)=dNdz;   Bsub(5,3)=dNdx
        Bsub(6,1)=dNdy;   Bsub(6,2)=dNdx
        

        
        
        end
