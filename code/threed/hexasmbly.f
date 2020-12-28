      subroutine hexasmbly(Kst,Ma,coord,elem,intval,prop,psp,pspl,IA,JA)

        implicit NONE
        integer intval(3)
        integer pspl(intval(2)+1)
        integer psp(pspl(intval(2)+1)-1)
        DOUBLE PRECISION Kst(intval(1)),Ma(3*intval(2))
        DOUBLE PRECISION coord(intval(2),3),prop(3)
        integer elem(intval(3),8)
        integer nnod,nelem,NZ,IA(intval(2)+1),JA(intval(1))
        
! Internal Variables        

        DOUBLE PRECISION nu,E0,rho0
        DOUBLE PRECISION,dimension (:,:),allocatable :: dN1k,dN2k,dN3k
        DOUBLE PRECISION,dimension (:,:),allocatable :: dN4k,dN5k,dN6k
        DOUBLE PRECISION,dimension (:,:),allocatable :: dN7k,dN8k,dN
        DOUBLE PRECISION,dimension (:,:),allocatable :: dN1ma,dN2ma
        DOUBLE PRECISION,dimension (:,:),allocatable :: dN3ma,dN4ma
        DOUBLE PRECISION,dimension (:,:),allocatable :: dN5ma,dN6ma
        DOUBLE PRECISION,dimension (:,:),allocatable :: dN7ma,dN8ma,Nval
        DOUBLE PRECISION,dimension (:,:),allocatable :: Xcoord,intp
        DOUBLE PRECISION,dimension (:,:),allocatable :: D,kelem,melem
        DOUBLE PRECISION,dimension (:,:),allocatable :: Bfull,Bsub
        DOUBLE PRECISION,dimension (:,:),allocatable :: Cdum,Jac,invJac
        DOUBLE PRECISION,dimension (:),allocatable :: N1,N2,N3,N4
        DOUBLE PRECISION,dimension (:),allocatable :: N5,N6,N7,N8
        DOUBLE PRECISION fact,zeta,eta,mu,alp,bet,detJac
        DOUBLE PRECISION x1,x2,x3,x4,x5,x6,x7,x8
        DOUBLE PRECISION y1,y2,y3,y4,y5,y6,y7,y8
        DOUBLE PRECISION z1,z2,z3,z4,z5,z6,z7,z8
        DOUBLE PRECISION dN1dz,dN2dz,dN3dz,dN4dz,dN5dz,dN6dz,dN7dz,dN8dz
        DOUBLE PRECISION dN1de,dN2de,dN3de,dN4de,dN5de,dN6de,dN7de,dN8de
        DOUBLE PRECISION dN1dm,dN2dm,dN3dm,dN4dm,dN5dm,dN6dm,dN7dm,dN8dm
        DOUBLE PRECISION N1val,N2val,N3val,N4val,N5val,N6val,N7val,N8val
        integer i,j,igp,ngp,ielem
        integer innod1,innod2,innod3,innod4
        integer innod5,innod6,innod7,innod8
        integer M,N,K,LDA,LDB,LDC,inod,flag,flag1
        integer ik(8),JAval,JApos

        dN1dz(zeta,eta,mu)=-(1.d0-eta)*(1.d0-mu)/8.d0
        dN2dz(zeta,eta,mu)=(1.d0-eta)*(1.d0-mu)/8.d0
        dN3dz(zeta,eta,mu)=(1.d0+eta)*(1.d0-mu)/8.d0
        dN4dz(zeta,eta,mu)=-(1.d0+eta)*(1.d0-mu)/8.d0
        dN5dz(zeta,eta,mu)=-(1.d0-eta)*(1.d0+mu)/8.d0
        dN6dz(zeta,eta,mu)=(1.d0-eta)*(1.d0+mu)/8.d0
        dN7dz(zeta,eta,mu)=(1.d0+eta)*(1.d0+mu)/8.d0
        dN8dz(zeta,eta,mu)=-(1.d0+eta)*(1.d0+mu)/8.d0

        dN1de(zeta,eta,mu)=-(1.d0-zeta)*(1.d0-mu)/8.d0
        dN2de(zeta,eta,mu)=-(1.d0+zeta)*(1.d0-mu)/8.d0
        dN3de(zeta,eta,mu)=(1.d0+zeta)*(1.d0-mu)/8.d0
        dN4de(zeta,eta,mu)=(1.d0-zeta)*(1.d0-mu)/8.d0
        dN5de(zeta,eta,mu)=-(1.d0-zeta)*(1.d0+mu)/8.d0
        dN6de(zeta,eta,mu)=-(1.d0+zeta)*(1.d0+mu)/8.d0
        dN7de(zeta,eta,mu)=(1.d0+zeta)*(1.d0+mu)/8.d0
        dN8de(zeta,eta,mu)=(1.d0-zeta)*(1.d0+mu)/8.d0

        dN1dm(zeta,eta,mu)=-(1.d0-zeta)*(1.d0-eta)/8.d0
        dN2dm(zeta,eta,mu)=-(1.d0+zeta)*(1.d0-eta)/8.d0
        dN3dm(zeta,eta,mu)=-(1.d0+zeta)*(1.d0+eta)/8.d0
        dN4dm(zeta,eta,mu)=-(1.d0-zeta)*(1.d0+eta)/8.d0
        dN5dm(zeta,eta,mu)=(1.d0-zeta)*(1.d0-eta)/8.d0
        dN6dm(zeta,eta,mu)=(1.d0+zeta)*(1.d0-eta)/8.d0
        dN7dm(zeta,eta,mu)=(1.d0+zeta)*(1.d0+eta)/8.d0
        dN8dm(zeta,eta,mu)=(1.d0-zeta)*(1.d0+eta)/8.d0

        N1val(zeta,eta,mu)=(1.d0-zeta)*(1.d0-eta)*(1.d0-mu)/8.d0
        N2val(zeta,eta,mu)=(1.d0+zeta)*(1.d0-eta)*(1.d0-mu)/8.d0
        N3val(zeta,eta,mu)=(1.d0+zeta)*(1.d0+eta)*(1.d0-mu)/8.d0
        N4val(zeta,eta,mu)=(1.d0-zeta)*(1.d0+eta)*(1.d0-mu)/8.d0
        N5val(zeta,eta,mu)=(1.d0-zeta)*(1.d0-eta)*(1.d0+mu)/8.d0
        N6val(zeta,eta,mu)=(1.d0+zeta)*(1.d0-eta)*(1.d0+mu)/8.d0
        N7val(zeta,eta,mu)=(1.d0+zeta)*(1.d0+eta)*(1.d0+mu)/8.d0
        N8val(zeta,eta,mu)=(1.d0-zeta)*(1.d0+eta)*(1.d0+mu)/8.d0
      
        allocate(dN1k(1:8,1:3),dN2k(1:8,1:3))
        allocate(dN3k(1:8,1:3),dN4k(1:8,1:3))
        allocate(dN5k(1:8,1:3),dN6k(1:8,1:3))
        allocate(dN7k(1:8,1:3),dN8k(1:8,1:3),dN(1:3,1:8))
        
        allocate(dN1ma(1:8,1:3),dN2ma(1:8,1:3))
        allocate(dN3ma(1:8,1:3),dN4ma(1:8,1:3))
        allocate(dN5ma(1:8,1:3),dN6ma(1:8,1:3))
        allocate(dN7ma(1:8,1:3),dN8ma(1:8,1:3))

        allocate(Bfull(1:6,1:24),D(1:6,1:6),Bsub(1:6,1:3))
        allocate(kelem(1:24,1:24),melem(1:24,1:24))
        
        allocate(Cdum(1:6,1:24))
        allocate(N1(1:8),N2(1:8),N3(1:8),N4(1:8))
        allocate(N5(1:8),N6(1:8),N7(1:8),N8(1:8))
        allocate(Xcoord(1:8,1:3),Jac(1:3,1:3),invJac(1:3,1:3))
        allocate(intp(1:8,1:3),Nval(1:3,1:24))

        
        nu=prop(1);  E0=prop(2);  rho0=prop(3)
        NZ=intval(1); nnod=intval(2); nelem=intval(3)
        
c! Shape function for stiffness matrix

        fact=sqrt(3.d0)
        intp(1,1)=-1.d0/fact; intp(1,2)=-1.d0/fact; intp(1,3)=-1.d0/fact         
        intp(2,1)=1.d0/fact;  intp(2,2)=-1.d0/fact; intp(2,3)=-1.d0/fact         
        intp(3,1)=1.d0/fact;  intp(3,2)=1.d0/fact;  intp(3,3)=-1.d0/fact         
        intp(4,1)=-1.d0/fact; intp(4,2)=1.d0/fact;  intp(4,3)=-1.d0/fact         
        intp(5,1)=-1.d0/fact; intp(5,2)=-1.d0/fact; intp(5,3)=1.d0/fact         
        intp(6,1)=1.d0/fact;  intp(6,2)=-1.d0/fact; intp(6,3)=1.d0/fact         
        intp(7,1)=1.d0/fact;  intp(7,2)=1.d0/fact;  intp(7,3)=1.d0/fact         
        intp(8,1)=-1.d0/fact; intp(8,2)=1.d0/fact;  intp(8,3)=1.d0/fact
       
        do igp=1,8

        zeta=intp(igp,1); eta=intp(igp,2); mu=intp(igp,3);

        dN1k(igp,1)=dN1dz(zeta,eta,mu); dN1k(igp,2)=dN1de(zeta,eta,mu)
        dN2k(igp,1)=dN2dz(zeta,eta,mu); dN2k(igp,2)=dN2de(zeta,eta,mu)
        dN3k(igp,1)=dN3dz(zeta,eta,mu); dN3k(igp,2)=dN3de(zeta,eta,mu)
        dN4k(igp,1)=dN4dz(zeta,eta,mu); dN4k(igp,2)=dN4de(zeta,eta,mu)
        dN5k(igp,1)=dN5dz(zeta,eta,mu); dN5k(igp,2)=dN5de(zeta,eta,mu)
        dN6k(igp,1)=dN6dz(zeta,eta,mu); dN6k(igp,2)=dN6de(zeta,eta,mu)
        dN7k(igp,1)=dN7dz(zeta,eta,mu); dN7k(igp,2)=dN7de(zeta,eta,mu)
        dN8k(igp,1)=dN8dz(zeta,eta,mu); dN8k(igp,2)=dN8de(zeta,eta,mu)

        dN1k(igp,3)=dN1dm(zeta,eta,mu); dN5k(igp,3)=dN5dm(zeta,eta,mu)
        dN2k(igp,3)=dN2dm(zeta,eta,mu); dN6k(igp,3)=dN6dm(zeta,eta,mu)
        dN3k(igp,3)=dN3dm(zeta,eta,mu); dN7k(igp,3)=dN7dm(zeta,eta,mu)
        dN4k(igp,3)=dN4dm(zeta,eta,mu); dN8k(igp,3)=dN8dm(zeta,eta,mu)
        
        enddo

c! Shape function values for mass matrix

        fact=sqrt(3.d0)
        intp(1,1)=-1.d0/fact; intp(1,2)=-1.d0/fact; intp(1,3)=-1.d0/fact         
        intp(2,1)=1.d0/fact;  intp(2,2)=-1.d0/fact; intp(2,3)=-1.d0/fact         
        intp(3,1)=1.d0/fact;  intp(3,2)=1.d0/fact;  intp(3,3)=-1.d0/fact         
        intp(4,1)=-1.d0/fact; intp(4,2)=1.d0/fact;  intp(4,3)=-1.d0/fact         
        intp(5,1)=-1.d0/fact; intp(5,2)=-1.d0/fact; intp(5,3)=1.d0/fact         
        intp(6,1)=1.d0/fact;  intp(6,2)=-1.d0/fact; intp(6,3)=1.d0/fact         
        intp(7,1)=1.d0/fact;  intp(7,2)=1.d0/fact;  intp(7,3)=1.d0/fact         
        intp(8,1)=-1.d0/fact; intp(8,2)=1.d0/fact;  intp(8,3)=1.d0/fact
        
        do igp=1,8

        zeta=intp(igp,1); eta=intp(igp,2); mu=intp(igp,3);

        dN1ma(igp,1)=dN1dz(zeta,eta,mu); dN1ma(igp,2)=dN1de(zeta,eta,mu)
        dN2ma(igp,1)=dN2dz(zeta,eta,mu); dN2ma(igp,2)=dN2de(zeta,eta,mu)
        dN3ma(igp,1)=dN3dz(zeta,eta,mu); dN3ma(igp,2)=dN3de(zeta,eta,mu)
        dN4ma(igp,1)=dN4dz(zeta,eta,mu); dN4ma(igp,2)=dN4de(zeta,eta,mu)
        dN5ma(igp,1)=dN5dz(zeta,eta,mu); dN5ma(igp,2)=dN5de(zeta,eta,mu)
        dN6ma(igp,1)=dN6dz(zeta,eta,mu); dN6ma(igp,2)=dN6de(zeta,eta,mu)
        dN7ma(igp,1)=dN7dz(zeta,eta,mu); dN7ma(igp,2)=dN7de(zeta,eta,mu)
        dN8ma(igp,1)=dN8dz(zeta,eta,mu); dN8ma(igp,2)=dN8de(zeta,eta,mu)

        dN1ma(igp,3)=dN1dm(zeta,eta,mu); dN5ma(igp,3)=dN5dm(zeta,eta,mu)
        dN2ma(igp,3)=dN2dm(zeta,eta,mu); dN6ma(igp,3)=dN6dm(zeta,eta,mu)
        dN3ma(igp,3)=dN3dm(zeta,eta,mu); dN7ma(igp,3)=dN7dm(zeta,eta,mu)
        dN4ma(igp,3)=dN4dm(zeta,eta,mu); dN8ma(igp,3)=dN8dm(zeta,eta,mu)

        N1(igp)=N1val(zeta,eta,mu); N5(igp)=N5val(zeta,eta,mu)
        N2(igp)=N2val(zeta,eta,mu); N6(igp)=N6val(zeta,eta,mu)
        N3(igp)=N3val(zeta,eta,mu); N7(igp)=N7val(zeta,eta,mu)
        N4(igp)=N4val(zeta,eta,mu); N8(igp)=N8val(zeta,eta,mu)
        
        enddo
! Arraniging the stiffness in the compressed sparse row storage form
        
        flag=1
        IA(:)=0
        JA(:)=0
        IA(1)=1
        do inod=1,nnod
        innod1=pspl(inod)
        innod2=pspl(inod+1)-1

            do i=innod1,innod2
            JA(flag)=(psp(i)-1)*3+1
            flag=flag+1
            JA(flag)=(psp(i)-1)*3+2
            flag=flag+1
            JA(flag)=(psp(i)-1)*3+3
            flag=flag+1
            enddo
        IA(3*(inod-1)+2)=flag

            do i=innod1,innod2
            JA(flag)=(psp(i)-1)*3+1
            flag=flag+1
            JA(flag)=(psp(i)-1)*3+2
            flag=flag+1
            JA(flag)=(psp(i)-1)*3+3
            flag=flag+1
            enddo
        IA(3*(inod-1)+3)=flag

            do i=innod1,innod2
            JA(flag)=(psp(i)-1)*3+1
            flag=flag+1
            JA(flag)=(psp(i)-1)*3+2
            flag=flag+1
            JA(flag)=(psp(i)-1)*3+3
            flag=flag+1
            enddo
        IA(3*(inod-1)+4)=flag
        enddo
        
! Calcaulting the elasticity matrix 
        D(:,:)=0.d0       
        fact=E0/((1.d0+nu)*(1.d0-2.d0*nu))
        D(1,1)=fact*(1.d0-nu);  D(1,2)=fact*nu;      D(1,3)=fact*nu
        D(2,1)=fact*nu;         D(2,2)=fact*(1-nu);  D(2,3)=fact*nu
        D(3,1)=fact*nu;         D(3,2)=fact*nu;      D(3,3)=fact*(1-nu)
        D(4,4)=fact*(0.5-nu); D(5,5)=fact*(0.5-nu); D(6,6)=fact*(0.5-nu)
        
!Assemblying the stiffness and mass matrices
        ngp=8
        Kst=0.d0;           Ma=0.d0;
        
        write(*,*) "rho:",rho0,"E:",E0

        do ielem=1,nelem
        innod1=elem(ielem,1)
        innod2=elem(ielem,2)
        innod3=elem(ielem,3)
        innod4=elem(ielem,4)
        innod5=elem(ielem,5)
        innod6=elem(ielem,6)
        innod7=elem(ielem,7)
        innod8=elem(ielem,8)
        
        x1=coord(innod1,1);  y1=coord(innod1,2);   z1=coord(innod1,3)   
        x2=coord(innod2,1);  y2=coord(innod2,2);   z2=coord(innod2,3)   
        x3=coord(innod3,1);  y3=coord(innod3,2);   z3=coord(innod3,3)   
        x4=coord(innod4,1);  y4=coord(innod4,2);   z4=coord(innod4,3)
        x5=coord(innod5,1);  y5=coord(innod5,2);   z5=coord(innod5,3)   
        x6=coord(innod6,1);  y6=coord(innod6,2);   z6=coord(innod6,3)   
        x7=coord(innod7,1);  y7=coord(innod7,2);   z7=coord(innod7,3)   
        x8=coord(innod8,1);  y8=coord(innod8,2);   z8=coord(innod8,3)
        
        kelem=0.d0;        melem=0.d0;        Nval=0.d0
        
        do igp=1,8
        Xcoord=0.d0
        Xcoord(1,1)=x1;     Xcoord(1,2)=y1;     Xcoord(1,3)=z1;
        Xcoord(2,1)=x2;     Xcoord(2,2)=y2;     Xcoord(2,3)=z2;
        Xcoord(3,1)=x3;     Xcoord(3,2)=y3;     Xcoord(3,3)=z3;
        Xcoord(4,1)=x4;     Xcoord(4,2)=y4;     Xcoord(4,3)=z4;
        Xcoord(5,1)=x5;     Xcoord(5,2)=y5;     Xcoord(5,3)=z5;
        Xcoord(6,1)=x6;     Xcoord(6,2)=y6;     Xcoord(6,3)=z6;
        Xcoord(7,1)=x7;     Xcoord(7,2)=y7;     Xcoord(7,3)=z7;
        Xcoord(8,1)=x8;     Xcoord(8,2)=y8;     Xcoord(8,3)=z8;
        
        dN=0.d0
        dN(1,1)=dN1k(igp,1); dN(2,1)=dN1k(igp,2); dN(3,1)=dN1k(igp,3)
        dN(1,2)=dN2k(igp,1); dN(2,2)=dN2k(igp,2); dN(3,2)=dN2k(igp,3)
        dN(1,3)=dN3k(igp,1); dN(2,3)=dN3k(igp,2); dN(3,3)=dN3k(igp,3)
        dN(1,4)=dN4k(igp,1); dN(2,4)=dN4k(igp,2); dN(3,4)=dN4k(igp,3)
        dN(1,5)=dN5k(igp,1); dN(2,5)=dN5k(igp,2); dN(3,5)=dN5k(igp,3)
        dN(1,6)=dN6k(igp,1); dN(2,6)=dN6k(igp,2); dN(3,6)=dN6k(igp,3)
        dN(1,7)=dN7k(igp,1); dN(2,7)=dN7k(igp,2); dN(3,7)=dN7k(igp,3)
        dN(1,8)=dN8k(igp,1); dN(2,8)=dN8k(igp,2); dN(3,8)=dN8k(igp,3)

        M=3; N=3; K=8; alp=1.d0; bet=0.d0; LDA=3; LDB=8; LDC=3;
        call DGEMM('N','N',M,N,K,alp,dN,LDA,Xcoord,LDB,bet,Jac,LDC)
        
        call invmat3(Jac,invJac,detJac)
        
        Bfull=0.d0; Bsub(:,:)=0.d0;  Cdum(:,:)=0.d0
        call Bsubform(dN1k(igp,1:3),Bsub,invJac)
        Bfull(1:6,1:3)=Bsub(1:6,1:3)

        call Bsubform(dN2k(igp,1:3),Bsub,invJac)
        Bfull(1:6,4:6)=Bsub(1:6,1:3)

        call Bsubform(dN3k(igp,1:3),Bsub,invJac)
        Bfull(1:6,7:9)=Bsub(1:6,1:3)

        call Bsubform(dN4k(igp,1:3),Bsub,invJac)
        Bfull(1:6,10:12)=Bsub(1:6,1:3)

        call Bsubform(dN5k(igp,1:3),Bsub,invJac)
        Bfull(1:6,13:15)=Bsub(1:6,1:3)

        call Bsubform(dN6k(igp,1:3),Bsub,invJac)
        Bfull(1:6,16:18)=Bsub(1:6,1:3)
        
        call Bsubform(dN7k(igp,1:3),Bsub,invJac)
        Bfull(1:6,19:21)=Bsub(1:6,1:3)
        
        call Bsubform(dN8k(igp,1:3),Bsub,invJac)
        Bfull(1:6,22:24)=Bsub(1:6,1:3)
                
        M=6; N=24; K=6; alp=detJac; bet=0.d0; LDA=6; LDB=6; LDC=6
        call DGEMM('N','N',M,N,K,alp,D,LDA,Bfull,LDB,bet,Cdum,LDC)

        M=24; N=24; K=6; alp=1.d0; bet=1.d0; LDA=6; LDB=6; LDC=24
        call DGEMM('T','N',M,N,K,alp,Bfull,LDA,Cdum,LDB,bet,kelem,LDC)

! Calculating elemental Ma matrix
        dN=0.d0
        dN(1,1)=dN1ma(igp,1); dN(2,1)=dN1ma(igp,2); dN(3,1)=dN1ma(igp,3)
        dN(1,2)=dN2ma(igp,1); dN(2,2)=dN2ma(igp,2); dN(3,2)=dN2ma(igp,3)
        dN(1,3)=dN3ma(igp,1); dN(2,3)=dN3ma(igp,2); dN(3,3)=dN3ma(igp,3)
        dN(1,4)=dN4ma(igp,1); dN(2,4)=dN4ma(igp,2); dN(3,4)=dN4ma(igp,3)
        dN(1,5)=dN5ma(igp,1); dN(2,5)=dN5ma(igp,2); dN(3,5)=dN5ma(igp,3)
        dN(1,6)=dN6ma(igp,1); dN(2,6)=dN6ma(igp,2); dN(3,6)=dN6ma(igp,3)
        dN(1,7)=dN7ma(igp,1); dN(2,7)=dN7ma(igp,2); dN(3,7)=dN7ma(igp,3)
        dN(1,8)=dN8ma(igp,1); dN(2,8)=dN8ma(igp,2); dN(3,8)=dN8ma(igp,3)

        M=3; N=3; K=8; alp=1.d0; bet=0.d0; LDA=3; LDB=8; LDC=3; Jac=0.d0
        call DGEMM('N','N',M,N,K,alp,dN,LDA,Xcoord,LDB,bet,Jac,LDC)
        
        call invmat3(Jac,invJac,detJac)
        
        Nval(:,:)=0.d0
        Nval(1,1)=N1(igp);   Nval(1,4)=N2(igp)
        Nval(1,7)=N3(igp);   Nval(1,10)=N4(igp)
        Nval(1,13)=N5(igp);  Nval(1,16)=N6(igp)
        Nval(1,19)=N7(igp);  Nval(1,22)=N8(igp)
        
        Nval(2,2)=N1(igp);   Nval(2,5)=N2(igp)
        Nval(2,8)=N3(igp);   Nval(2,11)=N4(igp)
        Nval(2,14)=N5(igp);  Nval(2,17)=N6(igp)
        Nval(2,20)=N7(igp);  Nval(2,23)=N8(igp)
        
        Nval(3,3)=N1(igp);   Nval(3,6)=N2(igp)
        Nval(3,9)=N3(igp);   Nval(3,12)=N4(igp)
        Nval(3,15)=N5(igp);  Nval(3,18)=N6(igp)
        Nval(3,21)=N7(igp);  Nval(3,24)=N8(igp)
        
        M=24; N=24; K=3; alp=detJac*rho0; bet=1.d0; LDA=3; LDB=3;LDC=24
        call DGEMM('T','N',M,N,K,alp,Nval,LDA,Nval,LDB,bet,melem,LDC)
        end do
        
        ik(1)=3*(innod1-1); ik(2)=3*(innod2-1)
        ik(3)=3*(innod3-1); ik(4)=3*(innod4-1)
        ik(5)=3*(innod5-1); ik(6)=3*(innod6-1)
        ik(7)=3*(innod7-1); ik(8)=3*(innod8-1)

! Assmeblying elemental stiffness matrix to final stiffness matrix
        do j=1,8
        do i=1,8
        
        inod=ik(i)+1; JAval=ik(j)+1; JApos=0
        M=IA(inod+1)-IA(inod)
        call findpos(JAval,JApos,JA(IA(inod):IA(inod+1)-1),M)
        JApos=IA(inod)-1+JApos
        Kst(JApos)=Kst(JApos)+kelem(3*(i-1)+1,3*(j-1)+1)
        
        inod=ik(i)+1; JAval=ik(j)+2; JApos=0
        M=IA(inod+1)-IA(inod)
        call findpos(JAval,JApos,JA(IA(inod):IA(inod+1)-1),M)
        JApos=IA(inod)-1+JApos
        Kst(JApos)=Kst(JApos)+kelem(3*(i-1)+1,3*(j-1)+2)

        inod=ik(i)+1; JAval=ik(j)+3; JApos=0
        M=IA(inod+1)-IA(inod)
        call findpos(JAval,JApos,JA(IA(inod):IA(inod+1)-1),M)
        JApos=IA(inod)-1+JApos
        Kst(JApos)=Kst(JApos)+kelem(3*(i-1)+1,3*(j-1)+3)

        inod=ik(i)+2; JAval=ik(j)+1; JApos=0
        M=IA(inod+1)-IA(inod)
        call findpos(JAval,JApos,JA(IA(inod):IA(inod+1)-1),M)
        JApos=IA(inod)-1+JApos
        Kst(JApos)=Kst(JApos)+kelem(3*(i-1)+2,3*(j-1)+1)

        inod=ik(i)+2; JAval=ik(j)+2; JApos=0
        M=IA(inod+1)-IA(inod)
        call findpos(JAval,JApos,JA(IA(inod):IA(inod+1)-1),M)
        JApos=IA(inod)-1+JApos
        Kst(JApos)=Kst(JApos)+kelem(3*(i-1)+2,3*(j-1)+2)
        
        inod=ik(i)+2; JAval=ik(j)+3; JApos=0
        M=IA(inod+1)-IA(inod)
        call findpos(JAval,JApos,JA(IA(inod):IA(inod+1)-1),M)
        JApos=IA(inod)-1+JApos
        Kst(JApos)=Kst(JApos)+kelem(3*(i-1)+2,3*(j-1)+3)
        
        inod=ik(i)+3; JAval=ik(j)+1; JApos=0
        M=IA(inod+1)-IA(inod)
        call findpos(JAval,JApos,JA(IA(inod):IA(inod+1)-1),M)
        JApos=IA(inod)-1+JApos
        Kst(JApos)=Kst(JApos)+kelem(3*(i-1)+3,3*(j-1)+1)

        inod=ik(i)+3; JAval=ik(j)+2; JApos=0
        M=IA(inod+1)-IA(inod)
        call findpos(JAval,JApos,JA(IA(inod):IA(inod+1)-1),M)
        JApos=IA(inod)-1+JApos
        Kst(JApos)=Kst(JApos)+kelem(3*(i-1)+3,3*(j-1)+2)
        
        inod=ik(i)+3; JAval=ik(j)+3; JApos=0
        M=IA(inod+1)-IA(inod)
        call findpos(JAval,JApos,JA(IA(inod):IA(inod+1)-1),M)
        JApos=IA(inod)-1+JApos
        Kst(JApos)=Kst(JApos)+kelem(3*(i-1)+3,3*(j-1)+3)

        end do
        end do

! Assmeblying elemental mass matrix to global mass matrix
        do i=1,8
        Ma(ik(i)+1)=Ma(ik(i)+1)+sum(melem(3*(i-1)+1,:))
        Ma(ik(i)+2)=Ma(ik(i)+2)+sum(melem(3*(i-1)+2,:));
        Ma(ik(i)+3)=Ma(ik(i)+3)+sum(melem(3*(i-1)+3,:));
        enddo
!        if(ielem==2) write(*,*) ielem,melem(ielem,:)
        end do

        deallocate(dN1k,dN2k,dN3k,dN4k,dN5k,dN6k,dN7k,dN8k,dN)
        deallocate(dN1ma,dN2ma,dN3ma,dN4ma,dN5ma,dN6ma,dN7ma,dN8ma)
        deallocate(Bfull,D,Bsub,kelem,melem,Cdum)
        deallocate(N1,N2,N3,N4,N5,N6,N7,N8,Xcoord,Jac,invJac,intp,Nval)
        
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
        
        write(*,*) "sum(absKst)",sum(abs(Kst)),"sum(Kst):",sum(Kst)
        write(*,*) "sum(Mass)",sum(Ma)
!        write(*,*) elem
!        write(*,*) Ma

        
  101   format(E13.6,X,E13.6,X,E13.6)
        
        
        end
