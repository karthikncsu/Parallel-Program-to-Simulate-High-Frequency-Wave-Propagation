        subroutine planestrain_fem_expl_parallel()
! Subroutine to solve using explicit method and lumped mass method

        use variables
        IMPLICIT NONE

        DOUBLE PRECISION unew(2*nnod),uold(2*nnod),uold0(2*nnod),prop(3)
        DOUBLE PRECISION f(2*nnod),A(2,2),B(2,2),C(2,2),ctime,etime
        DOUBLE PRECISION alpha,c0,time,bclv_ti,bcrv_ti,Fx,Fy,Ux,Uy
        DOUBLE PRECISION omega,Ac,Bc,Cc,excit,alp,bet,K1x,K1y
        DOUBLE PRECISION Mlump(2*nnod),L2,L2norm,MFLOPS,MFLOPS_avg
        integer inod,nsol,pnt1,pnt2,pnt3,pnt4,pnt5,pnt6,pnt7,pnt8,iface
        integer M,N,LDA,INCX,pomaxcnt,intvals(3),JAval,JApos,iproc,dum
        integer, dimension (:,:),allocatable :: dumpo
        integer, dimension (:),allocatable :: pocount,IA,JA
        INTEGER, dimension (:), allocatable :: rdofr,rdofs,ldofr,ldofs
        INTEGER rdofrcnt,rdofscnt,ldofrcnt,ldofscnt,tag1,tag2
        DOUBLE PRECISION,dimension (:),allocatable :: rdofrval,rdofsval
        DOUBLE PRECISION,dimension (:),allocatable :: ldofrval,ldofsval
        DOUBLE PRECISION,dimension (:),allocatable :: u1save
        DOUBLE PRECISION,dimension (:,:),allocatable ::coordsave
        character(len=100) :: outname,test        

        
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
!        write(*,*) "Finding point surrounding points"

        allocate(dumpo(1:nnod,1:9),pocount(1:nnod))
        call pspqudcnt(dumpo,elem,nelem,nnod,pocount,pomaxcnt)

        allocate(psp(1:pomaxcnt),pspl(nnod+1))
        call pspqudll(psp,pomaxcnt,pspl,nnod,dumpo)
    
        deallocate(dumpo,pocount)
        Write(*,*) myid,"CSR Matrix Size:", pomaxcnt
!--------------------------------------------------------------------
! Assemblying the matrices
!--------------------------------------------------------------------
        M=4*pomaxcnt
        allocate (Kstif(1:M),Mass(1:M),IA(1:nnod*2+1),JA(1:M))
!        write(*,*) "Assemblying matrices"
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
! Book keeping the DOFs need to be exchanged between the processors
!--------------------------------------------------------------------
        rdofrcnt=2*rnodrcnt;  rdofscnt=2*rnodscnt
        ldofrcnt=2*lnodrcnt;  ldofscnt=2*lnodscnt
        allocate(rdofr(1:rdofrcnt),rdofs(1:rdofscnt))
        allocate(ldofr(1:ldofrcnt),ldofs(1:ldofscnt))

        allocate(rdofrval(1:rdofrcnt),rdofsval(1:rdofscnt))
        allocate(ldofrval(1:ldofrcnt),ldofsval(1:ldofscnt))
        
        write(*,*) myid,ldofrcnt,ldofscnt,rdofrcnt,rdofscnt

        do inod=1,rnodrcnt
        rdofr(inod*2-1)=rnodr(inod)*2-1
        rdofr(inod*2)=rnodr(inod)*2
!        write(*,*)myid,"rnodr",coord(rnodr(inod),1),coord(rnodr(inod),2)
        enddo

        do inod=1,rnodscnt
        rdofs(inod*2-1)=rnods(inod)*2-1
        rdofs(inod*2)=rnods(inod)*2
!        write(*,*)myid,"rnods",coord(rnods(inod),1),coord(rnods(inod),2)
        enddo

        do inod=1,lnodrcnt
        ldofr(inod*2-1)=lnodr(inod)*2-1
        ldofr(inod*2)=lnodr(inod)*2
!        write(*,*)myid,"lnodr",coord(lnodr(inod),1),coord(lnodr(inod),2)
        enddo

        do inod=1,lnodscnt
        ldofs(inod*2-1)=lnods(inod)*2-1
        ldofs(inod*2)=lnods(inod)*2
!        write(*,*)myid,"lnods",coord(lnods(inod),1),coord(lnods(inod),2)
        enddo
        tag1=1; tag2=2;
!--------------------------------------------------------------------
! Post processing to an animation file
!--------------------------------------------------------------------
        nsol=1
        time=t0
        omega=2.0*PI*freq

        unew(:)=0.d0
        uold(:)=0.d0
        uold0(:)=0.d0
        

        if(myid.ne.0) then
! All processors sending the data to processor 0
        call MPI_Send(coord(ppnod,1:2),2*pptot,MPI_DOUBLE_PRECISION,
     &   0,100+(myid+1),comm1d,ierr)

        else
! Processor 0 receving the data and writing to a file
        write (outname, "('animation2d_parallel.dat')")
        
        open(25,file=outname)
        
        write(25,*) "npoints",sum(locpptot), "ntime",int(tf/tsave)+1
        do inod = 1,pptot 
        pnt1=ppnod(inod)
        write (25,*) coord(pnt1,1)
        end do

       
        do iproc=1,nprocs-1 
        allocate(coordsave(1:locpptot(iproc+1),1:2))

        call MPI_Recv(coordsave,2*locpptot(iproc+1),MPI_DOUBLE_PRECISION
     &  ,iproc,100+(iproc+1),comm1d,status,ierr)
     
        do inod = 1,locpptot(iproc+1)
        write (25,*) coordsave(inod,1)
        enddo

        
        deallocate(coordsave)
        enddo

        write(25,*) time
        do inod=1,sum(locpptot)
        write(25,*) 0.0
        enddo

        endif
       
        call MPI_cart_shift(comm1d,0,1,left,right,ierr)
!--------------------------------------------------------------------
! Time Stepping started
!--------------------------------------------------------------------
!        write(*,*) "Start of time stepping"
        allocate(u1save(1:maxval(locpptot)))
        
        ctime=MPI_Wtime()
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
        
        if(myid.ne.nprocs-1) rdofsval=uold(rdofs);
        if(myid.ne.0) ldofsval=uold(ldofs)
     
       call MPI_Sendrecv(rdofsval,rdofscnt,MPI_DOUBLE_PRECISION
     &        ,right,1,ldofrval,ldofrcnt,MPI_DOUBLE_PRECISION,left
     &              ,1,comm1d,status,ierr)

        call MPI_Sendrecv(ldofsval,ldofscnt,MPI_DOUBLE_PRECISION
     &        ,left,2,rdofrval,rdofrcnt,MPI_DOUBLE_PRECISION,right
     &              ,2,comm1d,status,ierr)
   
        if(myid.ne.nprocs-1) uold(rdofr)=rdofrval;        
        if(myid.ne.0) uold(ldofr)=ldofrval
        
        do inod=1,2*nnod
        pnt1=IA(inod)
        pnt2=IA(inod+1)-1
        M=pnt2-pnt1+1
        call dotproduct(Kstif(pnt1:pnt2),uold(JA(pnt1:pnt2)),M,alpha)
        f(inod)=dt*dt*(f(inod)-alpha)
        unew(inod)=f(inod)/Mlump(inod)+2.d0*uold(inod)-uold0(inod)
        end do

! Post processing section        
        if(nsol<=int(time/tsave)) then

        L2=sum(unew(:)**2)-sum(unew(ldofr)**2)-sum(unew(rdofr)**2)

        call MPI_REDUCE(L2,L2norm,1,MPI_DOUBLE_PRECISION
     &  ,MPI_SUM,0,comm1d,ierr)

        if(myid==0) then
        write(*,*) myid,"t=",time,"L2norm",sqrt(L2norm),int(time/dt)
        endif
       
        if(myid.ne.0) then
! All processors sending the data to processor 0
        call MPI_Send(unew(2*(ppnod-1)+2),pptot,MPI_DOUBLE_PRECISION,
     &   0,1000+(myid+1),comm1d,ierr)
        else
! Processor 0 receving the data and writing to a file
        write(25,*) time
        
        do inod = 1,pptot
        pnt1=ppnod(inod)
        write (25,*) unew(2*(pnt1-1)+2)
        end do
        
        do iproc=1,nprocs-1 
        
        call MPI_Recv(u1save(1:locpptot(iproc+1)),locpptot(iproc+1),
     &  MPI_DOUBLE_PRECISION,iproc,1000+(iproc+1),comm1d,status,ierr)
     
        do inod = 1,locpptot(iproc+1)
        write (25,*) u1save(inod)
        enddo

        enddo
        endif !(myid.ne.0)
        
        
        nsol=nsol+1
        endif
! End of Post processing section        
        
        uold0(:)=uold(:)
        uold(:)=unew(:)
             
        end do
!--------------------------------------------------------------------
! End of time stepping
!--------------------------------------------------------------------
        
        etime=MPI_wtime()
        etime=etime-ctime
        MFLOPS=8*pomaxcnt+16*nnod+nface*12
        MFLOPS=tf/dt*MFLOPS*1e-6

        call MPI_REDUCE(MFLOPS,MFLOPS_avg,1,MPI_DOUBLE_PRECISION
     &  ,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if(myid==0) then
        write(*,*) myid,"comp. Time:",etime,"MFLOPS:",MFLOPS_avg/etime
        endif
        
        CLOSE(25)
        
  101   format(E13.6,X,E13.6,X,E13.6)
        
        end

