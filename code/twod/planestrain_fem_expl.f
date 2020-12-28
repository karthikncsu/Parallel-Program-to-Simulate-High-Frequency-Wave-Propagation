        subroutine planestrain_fem_expl()
! Subroutine to solve using explicit method and lumped mass method

        use variables
        implicit NONE
        DOUBLE PRECISION unew(2*nnod),uold(2*nnod),uold0(2*nnod),prop(3)
        DOUBLE PRECISION f(2*nnod),A(2,2),B(2,2),C(2,2),ctime,etime
        DOUBLE PRECISION alpha,c0,time,bclv_ti,bcrv_ti,Fx,Fy,MFLOPS
        DOUBLE PRECISION omega,Ac,Bc,Cc,excit,alp,bet,K1x,K1y
        DOUBLE PRECISION Mlump(2*nnod)
        integer inod,nsol,pnt1,pnt2,pnt3,pnt4,pnt5,pnt6,pnt7,pnt8,iface
        integer M,N,LDA,INCX,pomaxcnt,intvals(3),JAval,JApos
        integer, dimension (:,:),allocatable :: dumpo
        integer, dimension (:),allocatable ::pocount,IA,JA
        
        f(:)=0.d0

        do iface=1,nface
!        if(int(face(iface,3))==1.0) then
        pnt1=int(face(iface,1))
        pnt2=int(face(iface,2))
        
        Fx=face(iface,4)*face(iface,6)/2.d0
        Fy=face(iface,5)*face(iface,6)/2.d0
        
        f(2*(pnt1-1)+1)=f(2*(pnt1-1)+1)+Fx
        f(2*(pnt1-1)+2)=f(2*(pnt1-1)+2)+Fy

        f(2*(pnt2-1)+1)=f(2*(pnt2-1)+1)+Fx
        f(2*(pnt2-1)+2)=f(2*(pnt2-1)+2)+Fy
!        write(*,*) iface,Fx,Fy
!        endif
        end do
!--------------------------------------------------------------------
! Points surrounding points
!--------------------------------------------------------------------
        write(*,*) "Finding point surrounding points"

        allocate(dumpo(1:nnod,1:9),pocount(1:nnod))
        call pspqudcnt(dumpo,elem,nelem,nnod,pocount,pomaxcnt)

        allocate(psp(1:pomaxcnt),pspl(nnod+1))
        call pspqudll(psp,pomaxcnt,pspl,nnod,dumpo)
    
        deallocate(dumpo,pocount)
        Write(*,*) "CSR Matrix Size:", pomaxcnt
!--------------------------------------------------------------------
! Assemblying the matrices
!--------------------------------------------------------------------
        M=4*pomaxcnt
        allocate (Kstif(1:M),Mass(1:M),IA(1:nnod*2+1),JA(1:M))
        write(*,*) "Assemblying matrices"
        prop(1)=nu;  prop(2)=E0;  prop(3)=rho0
        intvals(1)=M; intvals(2)=nnod; intvals(3)=nelem;
       call qudasmbly(Kstif,Mass,coord,elem,intvals,prop,psp,pspl,IA,JA)
       
! if lumped mass, only diagonal terms are needed. Currently only lumped
! mass is used. In future, if implict method is implemented, below code
! need to be conditioned properly
        do inod=1,nnod*2
        JAval=inod; JApos=0
        M=IA(inod+1)-IA(inod)
        call findpos(JAval,JApos,JA(IA(inod):IA(inod+1)-1),M)
        JApos=IA(inod)-1+JApos
        Mlump(inod)=Mass(JApos)
        if(Mlump(inod)==0) then
        write(*,*) "singular values in mass matrix"
        call exit(1)
        endif
        end do
        
        deallocate(Mass)
!--------------------------------------------------------------------
! Post processing to an animation file
!--------------------------------------------------------------------
        nsol=1
        time=t0
        omega=2.0*PI*freq

        unew(:)=0.d0
        uold(:)=0.d0
        uold0(:)=0.d0
        
        open(25,file='animation_wave.dat')
        write(25,*) "npoints",pptot,"Tpoints",int(tf/tsave)+1 

        do inod = 1,pptot
        pnt1=ppnod(inod)
        K1x=coord(pnt1,1)
        write (25,*) K1x
        end do
        
        write(25,*) time

        do inod = 1,pptot
        pnt1=ppnod(inod)
        write (25,*) unew(2*pnt1)
        end do
        
!--------------------------------------------------------------------
! Time Stepping started
!--------------------------------------------------------------------
        write(*,*) "Start of time stepping"
        call cpu_time(ctime)
        do while(time<=tf)

        time=time+dt
        if(time<10.0*PI/omega) then
           excit=sin(omega*time)*(sin(omega*time/10.0))**2
        else
           excit=0.d0
        end if
        
        f(:)=0.d0
        do iface=1,nface
!        if(int(face(iface,3))==1.0) then
        pnt1=int(face(iface,1))
        pnt2=int(face(iface,2))

        Fx=face(iface,4)*face(iface,6)/2.d0
        Fy=face(iface,5)*face(iface,6)/2.d0
        
        f(2*(pnt1-1)+1)=f(2*(pnt1-1)+1)+Fx*excit
        f(2*(pnt1-1)+2)=f(2*(pnt1-1)+2)+Fy*excit

        f(2*(pnt2-1)+1)=f(2*(pnt2-1)+1)+Fx*excit
        f(2*(pnt2-1)+2)=f(2*(pnt2-1)+2)+Fy*excit
!        endif
        end do

       
        do inod=1,2*nnod
        pnt1=IA(inod)
        pnt2=IA(inod+1)-1
        M=pnt2-pnt1+1
        call dotproduct(Kstif(pnt1:pnt2),uold(JA(pnt1:pnt2)),M,alpha)
        f(inod)=dt*dt*(f(inod)-alpha)
        unew(inod)=f(inod)/Mlump(inod)+2.d0*uold(inod)-uold0(inod)
        end do
        
        if(nsol<=int(time/tsave)) then
        write(25,*) time
        do inod = 1,pptot
        pnt1=ppnod(inod)
        write (25,*) unew(2*pnt1)
        end do
        
        write(*,*) "t=",time,"L2norm",sqrt(sum(unew(:)**2)),int(time/dt)
        nsol=nsol+1
        endif
        
        uold0(:)=uold(:)
        uold(:)=unew(:)
             
        end do
!--------------------------------------------------------------------
! End of time stepping
!--------------------------------------------------------------------
        
        call cpu_time(etime)
        etime=etime-ctime
        MFLOPS=8*pomaxcnt+16*nnod+nface*12
        MFLOPS=tf/dt*MFLOPS*1e-6/etime
        
        write(*,*) "comp. Time:",etime,"MFLOPS:",MFLOPS
       
        CLOSE(25)
  101   format(E13.6,X,E13.6,X,E13.6)
        
        
        end
