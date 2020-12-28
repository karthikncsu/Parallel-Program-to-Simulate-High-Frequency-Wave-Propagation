        subroutine FEM_expl()
! Subroutine to solve using explicit method and lumped mass method

        use variables
        implicit NONE
        DOUBLE PRECISION unew(3*nnod),uold(3*nnod),uold0(3*nnod),prop(3)
        DOUBLE PRECISION f(3*nnod),ctime,etime,L2norm
        DOUBLE PRECISION alpha,c0,time,bclv_ti,bcrv_ti,Fx,Fy,Fz
        DOUBLE PRECISION omega,Ac,Bc,Cc,excit,alp,bet,K1x,K1y,MFLOPS
        integer inod,nsol,pnt1,pnt2,pnt3,pnt4,pnt5,pnt6,pnt7,pnt8,iface
        integer M,N,LDA,INCX,pomaxcnt,intvals(3),JAval,JApos
        integer, dimension (:,:),allocatable :: dumpo
        integer, dimension (:),allocatable ::pocount,IA,JA
        integer*8 totop
!--------------------------------------------------------------------
! Points surrounding points
!--------------------------------------------------------------------
        write(*,*) "Finding point surrounding points"

        allocate(dumpo(1:nnod,1:27),pocount(1:nnod))
        call psphexcnt(dumpo,elem,nelem,nnod,pocount,pomaxcnt)

        allocate(psp(1:pomaxcnt),pspl(nnod+1))
        call psphexll(psp,pomaxcnt,pspl,nnod,dumpo)
        
        deallocate(dumpo,pocount)
!--------------------------------------------------------------------
! Assemblying the matrices
!--------------------------------------------------------------------
        M=9*pomaxcnt
        Write(*,*) "CSR Matrix Size:", M
        allocate (Kstif(1:M),Ma3d(1:3*nnod),IA(1:nnod*3+1),JA(1:M))
        write(*,*) "Assemblying matrices"
        prop(1)=nu;  prop(2)=E0;  prop(3)=rho0
        intvals(1)=M; intvals(2)=nnod; intvals(3)=nelem;
       call hexasmbly(Kstif,Ma3d,coord,elem,intvals,prop,psp,pspl,IA,JA)
       
c!--------------------------------------------------------------------
c! Post processing to an animation file
c!--------------------------------------------------------------------
        nsol=1
        time=t0
        omega=2.0*PI*freq

        unew(:)=0.d0
        uold(:)=0.d0
        uold0(:)=0.d0
        
        open(25,file='animation3d_wave.dat')
        write(25,*) "npoints",pptot

        do inod = 1,pptot
        pnt1=ppnod(inod)
        K1x=coord(pnt1,1)
        K1y=coord(pnt1,2)
        write (25,*) K1x,K1y
        end do
        write(25,*) time

        do inod = 1,pptot
        pnt1=ppnod(inod)
        write (25,*) unew(3*(pnt1-1)+3)
        end do
        
!--------------------------------------------------------------------
! Time Stepping started
!--------------------------------------------------------------------
        write(*,*) "Start of time stepping"
        ctime=MPI_WTIME()
        totop=0
        
        do while(time<=tf)

        time=time+dt
        if(time<10.0*PI/omega) then
           excit=sin(omega*time)*(sin(omega*time/10.0))**2
        else
           excit=0.d0
        end if
        
        f(:)=0.d0
        do iface=1,nface
!        if(int(face(iface,5))==1.0) then
        pnt1=int(face(iface,1))
        pnt2=int(face(iface,2))
        pnt3=int(face(iface,3))
        pnt4=int(face(iface,4))

        Fx=face(iface,6)*face(iface,9)*excit/4.d0
        Fy=face(iface,7)*face(iface,9)*excit/4.d0
        Fz=face(iface,8)*face(iface,9)*excit/4.d0
        
        f(3*(pnt1-1)+1)=f(3*(pnt1-1)+1)+Fx
        f(3*(pnt1-1)+2)=f(3*(pnt1-1)+2)+Fy
        f(3*(pnt1-1)+3)=f(3*(pnt1-1)+3)+Fz

        f(3*(pnt2-1)+1)=f(3*(pnt2-1)+1)+Fx
        f(3*(pnt2-1)+2)=f(3*(pnt2-1)+2)+Fy
        f(3*(pnt2-1)+3)=f(3*(pnt2-1)+3)+Fz

        f(3*(pnt3-1)+1)=f(3*(pnt3-1)+1)+Fx
        f(3*(pnt3-1)+2)=f(3*(pnt3-1)+2)+Fy
        f(3*(pnt3-1)+3)=f(3*(pnt3-1)+3)+Fz

        f(3*(pnt4-1)+1)=f(3*(pnt4-1)+1)+Fx
        f(3*(pnt4-1)+2)=f(3*(pnt4-1)+2)+Fy
        f(3*(pnt4-1)+3)=f(3*(pnt4-1)+3)+Fz
!        endif
        end do
!        Write(*,*) time,sum(f)
!        totop=totop+21*nface

        do inod=1,3*nnod
        pnt1=IA(inod)
        pnt2=IA(inod+1)-1
        M=pnt2-pnt1+1
        call dotproduct(Kstif(pnt1:pnt2),uold(JA(pnt1:pnt2)),M,alpha)
        f(inod)=dt*dt*(f(inod)-alpha)
        unew(inod)=f(inod)/Ma3d(inod)+2.d0*uold(inod)-uold0(inod)
        end do
        
        if(nsol<=int(time/tsave)) then
        write(25,*) time
 
        do inod = 1,pptot
        pnt1=ppnod(inod)
        write (25,*) unew(3*(pnt1-1)+3)
        end do

        L2norm=sqrt(sum(unew(:)**2))  
        write(*,*) "t=",time,"L2norm",L2norm,int(time/dt)
        nsol=nsol+1
        endif
        
        uold0(:)=uold(:)
        uold(:)=unew(:)
             
        end do
!--------------------------------------------------------------------
! End of time stepping
!--------------------------------------------------------------------
        
        etime=MPI_WTIME()
        etime=etime-ctime
        MFLOPS=(sum(IA(2:3*nnod+1)-IA(1:3*nnod))*2+6*3*nnod)
        MFLOPS=(tf/dt)*MFLOPS*1e-6/etime
        
        write(*,*) "Comp. Time:",etime,"MFLOPS:",MFLOPS
       
        CLOSE(25)
  101   format(E13.6,X,E13.6,X,E13.6)
        
        
        end
