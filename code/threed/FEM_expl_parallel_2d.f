        subroutine FEM_expl_parallel_2d()
! Subroutine to solve using explicit method and lumped mass method

        use variables
        implicit NONE
        DOUBLE PRECISION unew(3*nnod),uold(3*nnod),uold0(3*nnod),prop(3)
        DOUBLE PRECISION f(3*nnod),ctime,etime,L2,L2norm
        DOUBLE PRECISION alpha,c0,time,bclv_ti,bcrv_ti,Fx,Fy,Fz
        DOUBLE PRECISION omega,excit,alp,bet,K1x,K1y,MFLOPS,MFLOPS_avg
        integer inod,nsol,pnt1,pnt2,pnt3,pnt4,pnt5,pnt6,pnt7,pnt8,iface
        integer M,N,LDA,INCX,pomaxcnt,intvals(3),JAval,JApos,totop,iproc
        integer, dimension (:,:),allocatable :: dumpo
        integer, dimension (:),allocatable ::pocount,IA,JA
! Parallel Variables
        INTEGER rdofrcnt,rdofscnt,ldofrcnt,ldofscnt
        INTEGER rdofr(3*rnodrcnt),rdofs(3*rnodscnt)
        INTEGER ldofr(3*lnodrcnt),ldofs(3*lnodscnt)
        DOUBLE PRECISION rdofrval(3*rnodrcnt),rdofsval(3*rnodscnt)
        DOUBLE PRECISION ldofrval(3*lnodrcnt),ldofsval(3*lnodscnt)

        INTEGER tdofrcnt,tdofscnt,bdofrcnt,bdofscnt
        INTEGER tdofr(3*tnodrcnt),tdofs(3*tnodscnt)
        INTEGER bdofr(3*bnodrcnt),bdofs(3*bnodscnt)
        DOUBLE PRECISION tdofrval(3*tnodrcnt),tdofsval(3*tnodscnt)
        DOUBLE PRECISION bdofrval(3*bnodrcnt),bdofsval(3*bnodscnt)

        INTEGER ltdofrcnt,ltdofscnt,lbdofrcnt,lbdofscnt
        INTEGER ltdofr(3*ltnodrcnt),ltdofs(3*ltnodscnt)
        INTEGER lbdofr(3*lbnodrcnt),lbdofs(3*lbnodscnt)
        DOUBLE PRECISION ltdofrval(3*ltnodrcnt),ltdofsval(3*ltnodscnt)
        DOUBLE PRECISION lbdofrval(3*lbnodrcnt),lbdofsval(3*lbnodscnt)

        INTEGER rtdofrcnt,rtdofscnt,rbdofrcnt,rbdofscnt
        INTEGER rtdofr(3*rtnodrcnt),rtdofs(3*rtnodscnt)
        INTEGER rbdofr(3*rbnodrcnt),rbdofs(3*rbnodscnt)
        DOUBLE PRECISION rtdofrval(3*rtnodrcnt),rtdofsval(3*rtnodscnt)
        DOUBLE PRECISION rbdofrval(3*rbnodrcnt),rbdofsval(3*rbnodscnt)
        
        DOUBLE PRECISION,dimension (:,:),allocatable ::coordsave
        DOUBLE PRECISION,dimension (:),allocatable ::u1save
        
        character(len=100) :: outname,test        
              
!--------------------------------------------------------------------
! Points surrounding points
!--------------------------------------------------------------------
       write(*,*) "Solving the wave propogation problem"
        
        nnod=locnnod(myid+1)
        
        allocate(dumpo(1:nnod,1:27),pocount(1:nnod))
        call psphexcnt(dumpo,elem,nelem,nnod,pocount,pomaxcnt)
        
        if(pomaxcnt<0) then
        write(*,*) "pomaxcnt in fem_exmpl_parallel_2d.f is less than 0"
        write(*,*) "exiting the program, pomaxcnt:",pomaxcnt    
        endif

        allocate(psp(1:pomaxcnt),pspl(nnod+1))
        call psphexll(psp,pomaxcnt,pspl,nnod,dumpo)
        
        deallocate(dumpo,pocount)
!--------------------------------------------------------------------
! Assemblying the matrices
!--------------------------------------------------------------------
        if(myid==0) write(*,*) "Assemblying matrices"
        
        M=9*pomaxcnt
        Write(*,*) myid,"CSR Matrix Size:", M
        
        allocate (Kstif(1:M),Ma3d(1:3*nnod),IA(1:nnod*3+1),JA(1:M))
        prop(1)=nu;  prop(2)=E0;  prop(3)=rho0
        intvals(1)=M; intvals(2)=nnod; intvals(3)=nelem;
       call hexasmbly(Kstif,Ma3d,coord,elem,intvals,prop,psp,pspl,IA,JA)

!--------------------------------------------------------------------
! Book keeping the DOFs need to be exchanged between the processors
!--------------------------------------------------------------------

! left and right dofs

        rdofrcnt=3*rnodrcnt;  rdofscnt=3*rnodscnt
        ldofrcnt=3*lnodrcnt;  ldofscnt=3*lnodscnt

        do inod=1,rnodrcnt
        rdofr((inod-1)*3+1)=(rnodr(inod)-1)*3+1
        rdofr((inod-1)*3+2)=(rnodr(inod)-1)*3+2
        rdofr((inod-1)*3+3)=(rnodr(inod)-1)*3+3
        enddo

        do inod=1,rnodscnt
        rdofs((inod-1)*3+1)=(rnods(inod)-1)*3+1
        rdofs((inod-1)*3+2)=(rnods(inod)-1)*3+2
        rdofs((inod-1)*3+3)=(rnods(inod)-1)*3+3
        enddo

        do inod=1,lnodrcnt
        ldofr((inod-1)*3+1)=(lnodr(inod)-1)*3+1
        ldofr((inod-1)*3+2)=(lnodr(inod)-1)*3+2
        ldofr((inod-1)*3+3)=(lnodr(inod)-1)*3+3
        enddo

        do inod=1,lnodscnt
        ldofs((inod-1)*3+1)=(lnods(inod)-1)*3+1
        ldofs((inod-1)*3+2)=(lnods(inod)-1)*3+2
        ldofs((inod-1)*3+3)=(lnods(inod)-1)*3+3
        enddo
        
! Bottom and top dofs

        tdofrcnt=3*tnodrcnt;  tdofscnt=3*tnodscnt
        bdofrcnt=3*bnodrcnt;  bdofscnt=3*bnodscnt

        do inod=1,tnodrcnt
        tdofr((inod-1)*3+1)=(tnodr(inod)-1)*3+1
        tdofr((inod-1)*3+2)=(tnodr(inod)-1)*3+2
        tdofr((inod-1)*3+3)=(tnodr(inod)-1)*3+3
        enddo

        do inod=1,tnodscnt
        tdofs((inod-1)*3+1)=(tnods(inod)-1)*3+1
        tdofs((inod-1)*3+2)=(tnods(inod)-1)*3+2
        tdofs((inod-1)*3+3)=(tnods(inod)-1)*3+3
        enddo

        do inod=1,bnodrcnt
        bdofr((inod-1)*3+1)=(bnodr(inod)-1)*3+1
        bdofr((inod-1)*3+2)=(bnodr(inod)-1)*3+2
        bdofr((inod-1)*3+3)=(bnodr(inod)-1)*3+3
        enddo

        do inod=1,bnodscnt
        bdofs((inod-1)*3+1)=(bnods(inod)-1)*3+1
        bdofs((inod-1)*3+2)=(bnods(inod)-1)*3+2
        bdofs((inod-1)*3+3)=(bnods(inod)-1)*3+3
        enddo

! left top and right top dofs

        rtdofrcnt=3*rtnodrcnt;  rtdofscnt=3*rtnodscnt
        ltdofrcnt=3*ltnodrcnt;  ltdofscnt=3*ltnodscnt

        do inod=1,rtnodrcnt
        rtdofr((inod-1)*3+1)=(rtnodr(inod)-1)*3+1
        rtdofr((inod-1)*3+2)=(rtnodr(inod)-1)*3+2
        rtdofr((inod-1)*3+3)=(rtnodr(inod)-1)*3+3
        enddo

        do inod=1,rtnodscnt
        rtdofs((inod-1)*3+1)=(rtnods(inod)-1)*3+1
        rtdofs((inod-1)*3+2)=(rtnods(inod)-1)*3+2
        rtdofs((inod-1)*3+3)=(rtnods(inod)-1)*3+3
        enddo

        do inod=1,ltnodrcnt
        ltdofr((inod-1)*3+1)=(ltnodr(inod)-1)*3+1
        ltdofr((inod-1)*3+2)=(ltnodr(inod)-1)*3+2
        ltdofr((inod-1)*3+3)=(ltnodr(inod)-1)*3+3
        enddo

        do inod=1,ltnodscnt
        ltdofs((inod-1)*3+1)=(ltnods(inod)-1)*3+1
        ltdofs((inod-1)*3+2)=(ltnods(inod)-1)*3+2
        ltdofs((inod-1)*3+3)=(ltnods(inod)-1)*3+3
        enddo

! left and right bottom dofs

        rbdofrcnt=3*rbnodrcnt;  rbdofscnt=3*rbnodscnt
        lbdofrcnt=3*lbnodrcnt;  lbdofscnt=3*lbnodscnt

        do inod=1,rbnodrcnt
        rbdofr((inod-1)*3+1)=(rbnodr(inod)-1)*3+1
        rbdofr((inod-1)*3+2)=(rbnodr(inod)-1)*3+2
        rbdofr((inod-1)*3+3)=(rbnodr(inod)-1)*3+3
        enddo

        do inod=1,rbnodscnt
        rbdofs((inod-1)*3+1)=(rbnods(inod)-1)*3+1
        rbdofs((inod-1)*3+2)=(rbnods(inod)-1)*3+2
        rbdofs((inod-1)*3+3)=(rbnods(inod)-1)*3+3
        enddo

        do inod=1,lbnodrcnt
        lbdofr((inod-1)*3+1)=(lbnodr(inod)-1)*3+1
        lbdofr((inod-1)*3+2)=(lbnodr(inod)-1)*3+2
        lbdofr((inod-1)*3+3)=(lbnodr(inod)-1)*3+3
        enddo

        do inod=1,lbnodscnt
        lbdofs((inod-1)*3+1)=(lbnods(inod)-1)*3+1
        lbdofs((inod-1)*3+2)=(lbnods(inod)-1)*3+2
        lbdofs((inod-1)*3+3)=(lbnods(inod)-1)*3+3
        enddo

!--------------------------------------------------------------------
! End of book keeping the DOFs need to be exchanged between the processors
!--------------------------------------------------------------------

        call MPI_cart_shift(comm1d,0,1,left,right,ierr)
        call MPI_cart_shift(comm1d,1,1,bottom,top,ierr)

        if(proccoord(1).ne.nprocx-1.and.proccoord(2).ne.nprocy-1)then
          righttop=proccordtomyid(proccoord(1)+1,proccoord(2)+1)
        else
          righttop=MPI_PROC_NULL
        endif
     
        if(proccoord(1).ne.0.and.proccoord(2).ne.nprocy-1) then
            lefttop=proccordtomyid(proccoord(1)-1,proccoord(2)+1)
        else
            lefttop=MPI_PROC_NULL
        endif
        
        if(proccoord(1).ne.nprocx-1.and.proccoord(2).ne.0) then
          rightbottom=proccordtomyid(proccoord(1)+1,proccoord(2)-1)
        else
          rightbottom=MPI_PROC_NULL
        endif
        
        if(proccoord(1).ne.0.and.proccoord(2).ne.0) then
          leftbottom=proccordtomyid(proccoord(1)-1,proccoord(2)-1)
        else
          leftbottom=MPI_PROC_NULL
        endif 
!--------------------------------------------------------------------
! Post processing to an animation file
!--------------------------------------------------------------------
        nsol=1
        time=t0
        omega=2.d0*PI*freq

       
        if(myid.ne.0) then
! All processors sending the data to processor 0
        call MPI_Send(coord(ppnod,1:2),2*pptot,MPI_DOUBLE_PRECISION,
     &   0,100+(myid+1),comm1d,ierr)

        else
! Processor 0 receving the data and writing to a file
        write (outname, "('animation3d_parallel.dat')")
        
        open(25,file=outname)
        write(25,*) "npoints",sum(locpptot), "ntime",int(tf/tsave)+1
        do inod = 1,pptot
        pnt1=ppnod(inod)
        write (25,*) coord(pnt1,1),coord(pnt1,2)
        end do
        
        do iproc=1,nprocs-1 
        allocate(coordsave(1:locpptot(iproc+1),1:2))

        call MPI_Recv(coordsave,2*locpptot(iproc+1),MPI_DOUBLE_PRECISION
     &  ,iproc,100+(iproc+1),comm1d,status,ierr)
        do inod = 1,locpptot(iproc+1)
        write (25,*) coordsave(inod,1),coordsave(inod,2)
        enddo
        
        deallocate(coordsave)
        enddo

        write(25,*) time
        do inod=1,sum(locpptot)
        write(25,*) 0.0
        enddo

        endif

        if(myid==0) then
! BCAST is to make remaining processor wait untill the 1st processor
! writes the data to a file
            call MPI_BCAST(1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        else
            call MPI_BCAST(iproc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        endif        

!--------------------------------------------------------------------
! Time Stepping started
!--------------------------------------------------------------------
        if(myid==0) write(*,*) "Start of time stepping"
        ctime=MPI_Wtime()
        totop=0

        allocate(u1save(1:maxval(locpptot)))
        u1save=0.d0

        unew(:)=0.d0
        uold(:)=0.d0
        uold0(:)=0.d0
                
        do while(time<=tf)

        time=time+dt
        if(exc_typ==3) then
            if(time<10.0*PI/omega) then
            excit=sin(omega*time)*(sin(omega*time/10.0))**2
            else
            excit=0.d0
            end if
        else if(exc_typ==2) then
            excit=sin(omega*time)
        endif
        
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
        
!        totop=totop+21*nface
!        write(*,*) myid,time,nnod

!--------------------------------------------------------
!   Start of parallel Processing part inside the time stepping
!--------------------------------------------------------
! Right and left boundary nodes
        if(proccoord(1).ne.nprocx-1)   rdofsval=uold(rdofs)
        if(proccoord(1).ne.0)          ldofsval=uold(ldofs)

        call MPI_Sendrecv(rdofsval,rdofscnt,MPI_DOUBLE_PRECISION
     &  ,right,501,ldofrval,ldofrcnt,MPI_DOUBLE_PRECISION,left,501
     &  ,comm1d,status,ierr)

        call MPI_Sendrecv(ldofsval,ldofscnt,MPI_DOUBLE_PRECISION
     &  ,left,502,rdofrval,rdofrcnt,MPI_DOUBLE_PRECISION,right,502
     &  ,comm1d,status,ierr)

        if(proccoord(1).ne.nprocx-1)   uold(rdofr)=rdofrval
        if(proccoord(1).ne.0)          uold(ldofr)=ldofrval

! Top and bottom boundary nodes
        if(proccoord(2).ne.nprocy-1)   tdofsval=uold(tdofs)
        if(proccoord(2).ne.0)          bdofsval=uold(bdofs)

        call MPI_Sendrecv(tdofsval,tdofscnt,MPI_DOUBLE_PRECISION
     &  ,top,503,bdofrval,bdofrcnt,MPI_DOUBLE_PRECISION,bottom,503
     &  ,comm1d,status,ierr)

        call MPI_Sendrecv(bdofsval,bdofscnt,MPI_DOUBLE_PRECISION
     &  ,bottom,504,tdofrval,tdofrcnt,MPI_DOUBLE_PRECISION,top,504
     &  ,comm1d,status,ierr)

        if(proccoord(2).ne.nprocy-1)   uold(tdofr)=tdofrval
        if(proccoord(2).ne.0)          uold(bdofr)=bdofrval

! right top corner nodes
        if(proccoord(1).ne.nprocx-1.and.proccoord(2).ne.nprocy-1) 
     &  rtdofsval=uold(rtdofs)

        call MPI_Sendrecv(rtdofsval,rtdofscnt,MPI_DOUBLE_PRECISION
     &  ,righttop,505,lbdofrval,lbdofrcnt,MPI_DOUBLE_PRECISION
     &  ,leftbottom,505,comm1d,status,ierr)

        if(proccoord(1).ne.0.and.proccoord(2).ne.0)                
     &   uold(lbdofr)=lbdofrval

! Left top corner nodes             
        if(proccoord(1).ne.0.and.proccoord(2).ne.nprocy-1) 
     &   ltdofsval=uold(ltdofs)
        
        call MPI_Sendrecv(ltdofsval,ltdofscnt,MPI_DOUBLE_PRECISION
     &  ,lefttop,506,rbdofrval,rbdofrcnt,MPI_DOUBLE_PRECISION
     &  ,rightbottom,506,comm1d,status,ierr)
     
        if(proccoord(1).ne.nprocx-1.and.proccoord(2).ne.0) 
     &   uold(rbdofr)=rbdofrval

! right bottom corner nodes
        if(proccoord(1).ne.nprocx-1.and.proccoord(2).ne.0)
     &   rbdofsval=uold(rbdofs)
        
        call MPI_Sendrecv(rbdofsval,rbdofscnt,MPI_DOUBLE_PRECISION
     &  ,rightbottom,507,ltdofrval,ltdofrcnt,MPI_DOUBLE_PRECISION
     &  ,lefttop,507,comm1d,status,ierr)
     
        if(proccoord(1).ne.0.and.proccoord(2).ne.nprocy-1)
     &   uold(ltdofr)=ltdofrval

! left bottom corner nodes            
        if(proccoord(1).ne.0.and.proccoord(2).ne.0)
     &   lbdofsval=uold(lbdofs)

        call MPI_Sendrecv(lbdofsval,lbdofscnt,MPI_DOUBLE_PRECISION
     &  ,leftbottom,508,rtdofrval,rtdofrcnt,MPI_DOUBLE_PRECISION
     &  ,righttop,508,comm1d,status,ierr)
       
        if(proccoord(1).ne.nprocx-1.and.proccoord(2).ne.nprocy-1)
     &   uold(rtdofr)=rtdofrval
        
!--------------------------------------------------------
!   End of parallel Processing part
!--------------------------------------------------------
!        nnod=locnnod(myid+1)
        
        do inod=1,3*nnod
        pnt1=IA(inod)
        pnt2=IA(inod+1)-1
        M=pnt2-pnt1+1
        call dotproduct(Kstif(pnt1:pnt2),uold(JA(pnt1:pnt2)),M,alpha)
        f(inod)=dt*dt*(f(inod)-alpha)
        unew(inod)=f(inod)/Ma3d(inod)+2.d0*uold(inod)-uold0(inod)
        
!        totop=totop+(2*M-1)+7
        end do
        
!----------------------------------------------------------------- 
! Writing data at each time step
!-----------------------------------------------------------------
        if(nsol<=int(time/tsave)) then

        L2=sum(unew(:)**2)!-sum(unew(ldofr)**2)-sum(unew(rdofr)**2)      
        
        call MPI_REDUCE(L2,L2norm,1,MPI_DOUBLE_PRECISION
     &  ,MPI_SUM,0,comm1d,ierr)
              
        if(myid.ne.0) then
! All processors sending the data to processor 0
        call MPI_Send(unew(3*(ppnod-1)+3),pptot,MPI_DOUBLE_PRECISION,
     &   0,1000+(myid+1),comm1d,ierr)
        else
! Processor 0 receving the data and writing to a file
        write(25,*) time
        
        do inod = 1,pptot
        pnt1=ppnod(inod)
        write (25,*) unew(3*(pnt1-1)+3)
        end do
        
        do iproc=1,nprocs-1 

        call MPI_Recv(u1save(1:locpptot(iproc+1)),locpptot(iproc+1),
     &  MPI_DOUBLE_PRECISION,iproc,1000+(iproc+1),comm1d,status,ierr)
     
        do inod = 1,locpptot(iproc+1)
        write (25,*) u1save(inod)
        enddo

        enddo
        
        write(*,*) "t=",time,"L2norm",sqrt(L2norm),int(time/dt)        
        
        endif !(myid.ne.0)
        
        nsol=nsol+1
        endif !nsol<=int(time/tsave)
!-----------------------------------------------------------------
! End of Writing data at each time step
!-----------------------------------------------------------------        

        uold0=uold
        uold=unew
             
        end do
!--------------------------------------------------------------------
! End of time stepping
!--------------------------------------------------------------------
        
        etime=MPI_wtime()
        etime=etime-ctime

        MFLOPS=(sum(IA(2:3*nnod+1)-IA(1:3*nnod))*2+6*3*nnod)
        MFLOPS=(tf/dt)*MFLOPS*1e-6/etime
                
        call MPI_REDUCE(MFLOPS,MFLOPS_avg,1,MPI_DOUBLE_PRECISION
     &  ,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if(myid==0) then
        write(*,*) myid,"comp. Time:",etime,"MFLOPS:",MFLOPS_avg
        write(*,*) "Results file saved to the file:",outname
        write(*,*) "Number of postprocessing nodes:",sum(locpptot)
        CLOSE(25)
        
        open(100,file="comp_time.dat",ACCESS="Append")        
        write(100,*) sum(locnnod),int(time),nprocs,etime,int(MFLOPS_avg)
        close(100)
        endif !(myid==0)
        

  101   format(E13.6,X,E13.6,X,E13.6)
        
        
        end
