      subroutine qudasmbly(Kst,Ma,coord,elem,intval,prop,psp,pspl,IA,JA)

        implicit NONE
        integer intval(3)
        integer pspl(intval(2)+1)
        integer psp(pspl(intval(2)+1)-1)
        DOUBLE PRECISION Kst(intval(1)),Ma(intval(1))
        DOUBLE PRECISION coord(intval(2),2),prop(3)
        integer elem(intval(3),4)
        integer nnod,nelem,NZ,IA(intval(2)+1),JA(intval(1))
        
! Internal Variables        

        DOUBLE PRECISION nu,E0,rho0
        DOUBLE PRECISION,dimension (:,:),allocatable :: dN1,dN2,dN3,dN4
        DOUBLE PRECISION,dimension (:,:),allocatable :: B,D,kelem,melem
        DOUBLE PRECISION,dimension (:,:),allocatable :: melemd
        DOUBLE PRECISION,dimension (:,:),allocatable :: Cdum,C,Nval
        DOUBLE PRECISION,dimension (:),allocatable :: N1,N2,N3,N4
        DOUBLE PRECISION Jac,fact,dxdz,dxde,dydz,dyde,zeta,eta,alp,bet
        DOUBLE PRECISION x1,x2,x3,x4,y1,y2,y3,y4
        DOUBLE PRECISION dN1dx,dN2dx,dN3dx,dN4dx,dN1dy,dN2dy,dN3dy,dN4dy
        DOUBLE PRECISION dN1dz,dN2dz,dN3dz,dN4dz,dN1de,dN2de,dN3de,dN4de
        DOUBLE PRECISION N1val,N2val,N3val,N4val
        integer i,j,innod1,innod2,innod3,innod4,igp,ngp,ielem
        integer M,N,K,LDA,LDB,LDC,inod,flag,flag1
        integer ik(4),JAval,JApos

        dN1dz(zeta,eta)=-(1.d0-eta)/4.d0
        dN2dz(zeta,eta)=(1.d0-eta)/4.d0
        dN3dz(zeta,eta)=(1.d0+eta)/4.d0
        dN4dz(zeta,eta)=-(1.d0+eta)/4.d0

        dN1de(zeta,eta)=-(1.d0-zeta)/4.d0
        dN2de(zeta,eta)=-(1.d0+zeta)/4.d0
        dN3de(zeta,eta)=(1.d0+zeta)/4.d0
        dN4de(zeta,eta)=(1.d0-zeta)/4.d0

        N1val(zeta,eta)=(1.d0-zeta)*(1.d0-eta)/4.d0
        N2val(zeta,eta)=(1.d0+zeta)*(1.d0-eta)/4.d0
        N3val(zeta,eta)=(1.d0+zeta)*(1.d0+eta)/4.d0
        N4val(zeta,eta)=(1.d0-zeta)*(1.d0+eta)/4.d0
      
        allocate(dN1(1:4,1:2),dN2(1:4,1:2),dN3(1:4,1:2),dN4(1:4,1:2))
        allocate(B(1:3,1:8),D(1:3,1:3),kelem(1:8,1:8),C(1:2,1:2))
        allocate(Cdum(1:3,1:8),melem(1:8,1:8),melemd(1:8,1:8))
        allocate(N1(1:4),N2(1:4),N3(1:4),N4(1:4),Nval(1:2,1:8))
        
        nu=prop(1);  E0=prop(2);  rho0=prop(3)
        NZ=intval(1); nnod=intval(2); nelem=intval(3)
        
c! Shape function for stiffness matrix
!        write(*,*) "entered Quad assembly"

        fact=sqrt(3.d0)
        zeta=-1.d0/fact; eta=-1.d0/fact        
        dN1(1,1)=dN1dz(zeta,eta);  dN1(1,2)=dN1de(zeta,eta)
        dN2(1,1)=dN2dz(zeta,eta);  dN2(1,2)=dN2de(zeta,eta)
        dN3(1,1)=dN3dz(zeta,eta);  dN3(1,2)=dN3de(zeta,eta)
        dN4(1,1)=dN4dz(zeta,eta);  dN4(1,2)=dN4de(zeta,eta)

        zeta=1.d0/fact; eta=-1.d0/fact       
        dN1(2,1)=dN1dz(zeta,eta);  dN1(2,2)=dN1de(zeta,eta)
        dN2(2,1)=dN2dz(zeta,eta);  dN2(2,2)=dN2de(zeta,eta)
        dN3(2,1)=dN3dz(zeta,eta);  dN3(2,2)=dN3de(zeta,eta)
        dN4(2,1)=dN4dz(zeta,eta);  dN4(2,2)=dN4de(zeta,eta)

        zeta=1.d0/fact; eta=1.d0/fact      
        dN1(3,1)=dN1dz(zeta,eta);  dN1(3,2)=dN1de(zeta,eta)
        dN2(3,1)=dN2dz(zeta,eta);  dN2(3,2)=dN2de(zeta,eta)
        dN3(3,1)=dN3dz(zeta,eta);  dN3(3,2)=dN3de(zeta,eta)
        dN4(3,1)=dN4dz(zeta,eta);  dN4(3,2)=dN4de(zeta,eta)

        zeta=-1.d0/fact; eta=1.d0/fact        
        dN1(4,1)=dN1dz(zeta,eta);  dN1(4,2)=dN1de(zeta,eta)
        dN2(4,1)=dN2dz(zeta,eta);  dN2(4,2)=dN2de(zeta,eta)
        dN3(4,1)=dN3dz(zeta,eta);  dN3(4,2)=dN3de(zeta,eta)
        dN4(4,1)=dN4dz(zeta,eta);  dN4(4,2)=dN4de(zeta,eta)

! Shape function evaluation for Ma matrix
        fact=sqrt(3.d0)
        zeta=-1.d0/fact; eta=-1.d0/fact
        N1(1)=N1val(zeta,eta); N2(1)=N2val(zeta,eta)
        N3(1)=N3val(zeta,eta); N4(1)=N4val(zeta,eta)
        
        zeta=1.d0/fact; eta=-1.d0/fact
        N1(2)=N1val(zeta,eta); N2(2)=N2val(zeta,eta)
        N3(2)=N3val(zeta,eta); N4(2)=N4val(zeta,eta)
        
        zeta=1.d0/fact; eta=1.d0/fact
        N1(3)=N1val(zeta,eta); N2(3)=N2val(zeta,eta)
        N3(3)=N3val(zeta,eta); N4(3)=N4val(zeta,eta)
        
        zeta=-1.d0/fact; eta=1.d0/fact
        N1(4)=N1val(zeta,eta); N2(4)=N2val(zeta,eta)
        N3(4)=N3val(zeta,eta); N4(4)=N4val(zeta,eta)

! Arranign the stiffness and mass matrix in the compressed sparse row 
! storage form
        
        flag=1
        IA(:)=0
        JA(:)=0
        IA(1)=1
        do inod=1,nnod
        innod1=pspl(inod)
        innod2=pspl(inod+1)-1

            do i=innod1,innod2
            JA(flag)=psp(i)*2-1
            flag=flag+1
            JA(flag)=psp(i)*2
            flag=flag+1
            enddo
        IA(2*inod)=flag

            do i=innod1,innod2
            JA(flag)=psp(i)*2-1
            flag=flag+1
            JA(flag)=psp(i)*2
            flag=flag+1
            enddo
        IA(2*inod+1)=flag
        enddo

! Calcaulting the elasticity matrix 
        B(:,:)=0.d0;  Kst(:)=0.d0; Ma(:)=0.d0       
        fact=E0/((1.d0+nu)*(1.d0-2.d0*nu))
        D(1,1)=fact*(1.d0-nu);  D(1,2)=fact*nu;      D(1,3)=0.d0
        D(2,1)=fact*nu;         D(2,2)=fact*(1-nu);  D(2,3)=0.d0
        D(3,1)=0.d0;         D(3,2)=0.d0;         D(3,3)=fact*(0.5-nu)
        
!Assemblying the stiffness and mass matrices
        ngp=4
        do ielem=1,nelem
        innod1=elem(ielem,1)
        innod2=elem(ielem,2)
        innod3=elem(ielem,3)
        innod4=elem(ielem,4)
        
        x1=coord(innod1,1);  y1=coord(innod1,2)   
        x2=coord(innod2,1);  y2=coord(innod2,2)   
        x3=coord(innod3,1);  y3=coord(innod3,2)   
        x4=coord(innod4,1);  y4=coord(innod4,2)
        
        kelem(:,:)=0.d0
        melem(:,:)=0.d0
      
        do igp=1,ngp
        B(:,:)=0.d0
        Cdum(:,:)=0.d0
        Nval(:,:)=0.d0
        
        dxdz=dN1(igp,1)*x1+dN2(igp,1)*x2+dN3(igp,1)*x3+dN4(igp,1)*x4
        dxde=dN1(igp,2)*x1+dN2(igp,2)*x2+dN3(igp,2)*x3+dN4(igp,2)*x4

        dydz=dN1(igp,1)*y1+dN2(igp,1)*y2+dN3(igp,1)*y3+dN4(igp,1)*y4
        dyde=dN1(igp,2)*y1+dN2(igp,2)*y2+dN3(igp,2)*y3+dN4(igp,2)*y4
        
        Jac=dxdz*dyde-dydz*dxde

! Calculating elemental stiffness matrix
        B(1,1)=dyde*dN1(igp,1)-dxde*dN1(igp,2)
        B(1,3)=dyde*dN2(igp,1)-dxde*dN2(igp,2)
        B(1,5)=dyde*dN3(igp,1)-dxde*dN3(igp,2)
        B(1,7)=dyde*dN4(igp,1)-dxde*dN4(igp,2)

        B(2,2)=-dydz*dN1(igp,1)+dxdz*dN1(igp,2)
        B(2,4)=-dydz*dN2(igp,1)+dxdz*dN2(igp,2)
        B(2,6)=-dydz*dN3(igp,1)+dxdz*dN3(igp,2)
        B(2,8)=-dydz*dN4(igp,1)+dxdz*dN4(igp,2)

        B(3,1)=-dydz*dN1(igp,1)+dxdz*dN1(igp,2)
        B(3,2)=dyde*dN1(igp,1)-dxde*dN1(igp,2)
        B(3,3)=-dydz*dN2(igp,1)+dxdz*dN2(igp,2)
        B(3,4)=dyde*dN2(igp,1)-dxde*dN2(igp,2)
        B(3,5)=-dydz*dN3(igp,1)+dxdz*dN3(igp,2)
        B(3,6)=dyde*dN3(igp,1)-dxde*dN3(igp,2)
        B(3,7)=-dydz*dN4(igp,1)+dxdz*dN4(igp,2)
        B(3,8)=dyde*dN4(igp,1)-dxde*dN4(igp,2)
        
        B=B/Jac
    
        M=3; N=8; K=3; alp=Jac; bet=0.d0; LDA=3; LDB=3; LDC=3
        call DGEMM('N','N',M,N,K,alp,D,LDA,B,LDB,bet,Cdum,LDC)

        M=8; N=8; K=3; alp=1.d0; bet=1.d0; LDA=3; LDB=3; LDC=8
        call DGEMM('T','N',M,N,K,alp,B,LDA,Cdum,LDB,bet,kelem,LDC)

! Calculating elemental Ma matrix
        Nval(1,1)=N1(igp); Nval(1,3)=N2(igp)
        Nval(1,5)=N3(igp); Nval(1,7)=N4(igp)
        Nval(2,2)=N1(igp); Nval(2,4)=N2(igp)
        Nval(2,6)=N3(igp); Nval(2,8)=N4(igp)

        M=8; N=8; K=2; alp=Jac*rho0; bet=1.d0; LDA=2; LDB=2; LDC=8
        call DGEMM('T','N',M,N,K,alp,Nval,LDA,Nval,LDB,bet,melem,LDC)
        
        end do
        
        melemd=melem;
        melem=0.d0;
        
        do i=1,8
        melem(i,i)=sum(melemd(i,:))
        enddo
        
        ik(1)=2*(innod1-1); ik(2)=2*(innod2-1)
        ik(3)=2*(innod3-1); ik(4)=2*(innod4-1)

        do j=1,4
        do i=1,4
        inod=ik(i)+1; JAval=ik(j)+1; JApos=0
        M=IA(inod+1)-IA(inod)
        call findpos(JAval,JApos,JA(IA(inod):IA(inod+1)-1),M)
        JApos=IA(inod)-1+JApos
        Kst(JApos)=Kst(JApos)+kelem(2*i-1,2*j-1)
        Ma(JApos)=Ma(JApos)+melem(2*i-1,2*j-1)
        
        inod=ik(i)+1; JAval=ik(j)+2; JApos=0
        M=IA(inod+1)-IA(inod)
        call findpos(JAval,JApos,JA(IA(inod):IA(inod+1)-1),M)
        JApos=IA(inod)-1+JApos
        Kst(JApos)=Kst(JApos)+kelem(2*i-1,2*j)
        Ma(JApos)=Ma(JApos)+melem(2*i-1,2*j)

        inod=ik(i)+2; JAval=ik(j)+1; JApos=0
        M=IA(inod+1)-IA(inod)
        call findpos(JAval,JApos,JA(IA(inod):IA(inod+1)-1),M)
        JApos=IA(inod)-1+JApos
        Kst(JApos)=Kst(JApos)+kelem(2*i,2*j-1)
        Ma(JApos)=Ma(JApos)+melem(2*i,2*j-1)

        inod=ik(i)+2; JAval=ik(j)+2; JApos=0
        M=IA(inod+1)-IA(inod)
        call findpos(JAval,JApos,JA(IA(inod):IA(inod+1)-1),M)
        JApos=IA(inod)-1+JApos
        Kst(JApos)=Kst(JApos)+kelem(2*i,2*j)
        Ma(JApos)=Ma(JApos)+melem(2*i,2*j)

        end do
        end do
        
        end do
        deallocate(dN1,dN2,dN3,dN4,B,D,kelem,melem,Cdum,C,melemd)
        deallocate(N1,N2,N3,N4,Nval)        
        
!        allocate(kelem(1:nnod,1:nnod))
!        kelem(:,:)=0.d0
!        write(*,*) JA
c        write(*,*) IA

c        do i=1,nnod
c            do j=IA(i),IA(i+1)-1
c            kelem(i,JA(j))=Kst(j)
c            enddo        
c        enddo        
c        do i=1,nnod
c        do j=1,nnod
c        write(*,*) i,j,kelem(i,j)
c        enddo
c        enddo
!        deallocate(kelem)
        
        write(*,*) "absKst",sum(abs(Kst)),"Kst:",sum(Kst),"Mass",sum(Ma)

        
  101   format(E13.6,X,E13.6,X,E13.6)
        
        
        end
