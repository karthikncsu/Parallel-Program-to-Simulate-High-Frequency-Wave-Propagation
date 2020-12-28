        subroutine divisionofgrid3d_1dx_parallel()
! Subroutine to division the mesh into sub-meshes for parallel code 
        use variables
        
        IMPLICIT NONE
        DOUBLE PRECISION coordmin(spcdim),coordmax(spcdim)
        DOUBLE PRECISION dlenx,dminx,dmaxx,coordlen(spcdim),fact
        DOUBLE PRECISION dleny,dminy,dmaxy,xleft,xright,ytop,ybottom
        INTEGER axis,ispc,ielem,maxb,bndcnt,inod,maxf,maxe
        INTEGER genflag,iproc,allstuff(4)
        INTEGER ipnt1,ipnt2,ipnt3,ipnt4,ipnt5,ipnt6,ipnt7,ipnt8,iface

        INTEGER,dimension (:), allocatable :: dumlnodr,dumlnods
        INTEGER,dimension (:), allocatable :: dumrnodr,dumrnods
        INTEGER,dimension (:), allocatable :: locelem
        INTEGER dumnnodls,dumnnodlr,dumnnodrs,dumnnodrr

        INTEGER,dimension (:), allocatable :: dumppnod
        DOUBLE PRECISION, dimension (:,:), allocatable :: dumcoord
        DOUBLE PRECISION, dimension (:,:), allocatable :: dumface
        INTEGER,dimension (:,:), allocatable :: dumelem
        INTEGER dumnelem,dumnnod,dumnface,dumpptot
        DOUBLE PRECISION, dimension (:), allocatable :: xc,yc,zc

        INTEGER,dimension (:), allocatable :: mapold2new

        INTEGER iprocx,iprocy
        INTEGER proccoordall(2,nprocs)     
        DOUBLE PRECISION ix,iy,iz
        DOUBLE PRECISION ipmaxx,ipminx,ipmaxy,ipminy
        DOUBLE PRECISION xmin,xmax,ymin,ymax,stime,etime
        INTEGER dumimin,dumimax

!-----------------------------------------------------------------
!Start of the program
!-----------------------------------------------------------------
        stime=MPI_WTIME()
        call MPI_CART_COORDS(comm1d,myid,2,proccoord,ierr)
        
        call MPI_ALLGATHER(proccoord,2,MPI_INTEGER,proccoordall,2,
     &   MPI_INTEGER,comm1d,ierr)
     
        allocate(proccordtomyid(0:nprocx-1,0:nprocy-1))
        proccordtomyid=0

         do iproc=1,nprocs
         ipnt1=proccoordall(1,iproc)
         ipnt2=proccoordall(2,iproc)

         proccordtomyid(ipnt1,ipnt2)=iproc-1
         enddo

        if(myid==0) then
            call readgrid3d()
            allstuff(1)=nelem; allstuff(2)=nnod;
            allstuff(3)=nface; allstuff(4)=pptot
            call MPI_BCAST(allstuff,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

          call MPI_BCAST(elem,8*nelem,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

          call MPI_BCAST(coord,3*nnod,MPI_DOUBLE_PRECISION
     &     ,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(face,9*nface,MPI_DOUBLE_PRECISION
     &     ,0,MPI_COMM_WORLD,ierr)
     
          call MPI_BCAST(ppnod,pptot,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        else
            call MPI_BCAST(allstuff,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

            nelem=allstuff(1); nnod=allstuff(2);
            nface=allstuff(3); pptot=allstuff(4)
            
            allocate(elem(1:nelem,1:8),coord(1:nnod,1:3))
            allocate(face(1:nface,1:9),ppnod(1:pptot))

          call MPI_BCAST(elem,8*nelem,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

          call MPI_BCAST(coord,3*nnod,MPI_DOUBLE_PRECISION
     &     ,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(face,9*nface,MPI_DOUBLE_PRECISION
     &     ,0,MPI_COMM_WORLD,ierr)
     
          call MPI_BCAST(ppnod,pptot,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        endif
        
        if(myid==0) write(*,*) "sent all the elements to
     &    different processors"

! aldumating the memory to variables
        maxe=nelem*2/(nprocx*nprocy)
        maxf=nface*2/(nprocx*nprocy)
        maxb=nnod*2/(nprocx*nprocy)
        
        allocate(dumlnodr(1:maxb),dumlnods(1:maxb))
        allocate(dumrnodr(1:maxb),dumrnods(1:maxb))
        allocate(locelem(1:maxe))
        allocate(dumface(1:maxf,1:9),xc(1:maxe),yc(1:maxe),zc(1:maxe))

        dumlnodr=0;   dumlnods=0;    dumrnodr=0;     dumrnods=0; 
        dumnnodls=0;  dumnnodlr=0;   dumnnodrs=0;    dumnnodrr=0;
        dumnelem=0;   locelem=0

       
        dlenx=maxval(coord(:,1))-minval(coord(:,1))
        dminx=minval(coord(:,1))-0.01*dlenx;
        dmaxx=maxval(coord(:,1))+0.01*dlenx

        dleny=maxval(coord(:,2))-minval(coord(:,2))
        dminy=minval(coord(:,2))-0.01*dleny;
        dmaxy=maxval(coord(:,2))+0.01*dleny

        fact=1.023567
        dlenx=fact*dlenx/dble(nprocx);
        dleny=fact*dleny/dble(nprocy);

        if(myid==0) then
        write(*,*) "min x-coord(actual,used)",minval(coord(:,1)),dminx
        write(*,*) "min y-coord(actual,used)",minval(coord(:,2)),dminy
        write(*,*) "Length of each domain along x-dir:",dlenx
        write(*,*) "Length of each domain along y-dir:",dleny
        endif
        allocate(label(1:maxe))
        label=1000

        call MPI_cart_shift(comm1d,0,1,left,right,ierr)        
!-----------------------------------------------------
! finding the elements of each processors
!-----------------------------------------------------

        iprocx=proccoord(1); iprocy=proccoord(2);
        xleft=dminx+iprocx*dlenx;   xright=dminx+(iprocx+1)*dlenx
        ybottom=dminy+iprocy*dleny;   ytop=dminy+(iprocy+1)*dleny

       do ielem=1,nelem
        ipminx=minval(coord(elem(ielem,:),1))
        ipmaxx=maxval(coord(elem(ielem,:),1))

        ipnt1=elem(ielem,1); ipnt2=elem(ielem,2)
        ipnt3=elem(ielem,3); ipnt4=elem(ielem,4)
        ipnt5=elem(ielem,5); ipnt6=elem(ielem,6)
        ipnt7=elem(ielem,7); ipnt8=elem(ielem,8) 
        
              
        if(ipminx>xleft.and.ipmaxx<xright) then 
     
         dumnelem=dumnelem+1;     locelem(dumnelem)=ielem

         xc(dumnelem)=sum(coord(elem(ielem,:),1))/8.0
         yc(dumnelem)=sum(coord(elem(ielem,:),2))/8.0
         zc(dumnelem)=maxval(coord(elem(ielem,:),3))

         label(dumnelem)=myid
        
        elseif(ipminx<xleft.and.ipmaxx>xleft) then

         dumnelem=dumnelem+1;     locelem(dumnelem)=ielem

! Left side receviing nodes
         dumnnodlr=dumnnodlr+1;     dumlnodr(dumnnodlr)=ipnt1         
         dumnnodlr=dumnnodlr+1;     dumlnodr(dumnnodlr)=ipnt4
         dumnnodlr=dumnnodlr+1;     dumlnodr(dumnnodlr)=ipnt5
         dumnnodlr=dumnnodlr+1;     dumlnodr(dumnnodlr)=ipnt8

! Left side Sending nodes
        dumnnodls=dumnnodls+1;      dumlnods(dumnnodls)=ipnt2
        dumnnodls=dumnnodls+1;      dumlnods(dumnnodls)=ipnt3
        dumnnodls=dumnnodls+1;      dumlnods(dumnnodls)=ipnt6
        dumnnodls=dumnnodls+1;      dumlnods(dumnnodls)=ipnt7

        xc(dumnelem)=sum(coord(elem(ielem,:),1))/8.0
        yc(dumnelem)=sum(coord(elem(ielem,:),2))/8.0
        zc(dumnelem)=maxval(coord(elem(ielem,:),3))

        label(dumnelem)=myid+left
 
       else if(ipminx<xright.and.ipmaxx>xright) then     

         dumnelem=dumnelem+1;     locelem(dumnelem)=ielem
        
! Right side receviing nodes
        dumnnodrr=dumnnodrr+1;       dumrnodr(dumnnodrr)=ipnt2
        dumnnodrr=dumnnodrr+1;       dumrnodr(dumnnodrr)=ipnt3
        dumnnodrr=dumnnodrr+1;       dumrnodr(dumnnodrr)=ipnt6
        dumnnodrr=dumnnodrr+1;       dumrnodr(dumnnodrr)=ipnt7

! Right side Sending nodes
        dumnnodrs=dumnnodrs+1;       dumrnods(dumnnodrs)=ipnt1
        dumnnodrs=dumnnodrs+1;       dumrnods(dumnnodrs)=ipnt4
        dumnnodrs=dumnnodrs+1;       dumrnods(dumnnodrs)=ipnt5
        dumnnodrs=dumnnodrs+1;       dumrnods(dumnnodrs)=ipnt8

        xc(dumnelem)=sum(coord(elem(ielem,:),1))/8.0
        yc(dumnelem)=sum(coord(elem(ielem,:),2))/8.0
        zc(dumnelem)=maxval(coord(elem(ielem,:),3))

        label(dumnelem)=myid+right
        
        endif 
       enddo !ielem

!-----------------------------------------------------
! End of finding the elements of each processors
!-----------------------------------------------------
!---------------------------------------------------------------
!     Writing label of the mesh
!---------------------------------------------------------------

        if(myid.ne.0) then
! All processors sending the data to processor 0
        call MPI_Send(dumnelem,1,MPI_INTEGER,0,5*(myid-1)+1,comm1d,ierr)
 
        call MPI_Send(xc,dumnelem,MPI_DOUBLE_PRECISION,
     &   0,5*(myid-1)+2,comm1d,ierr)
        call MPI_Send(yc,dumnelem,MPI_DOUBLE_PRECISION,
     &   0,5*(myid-1)+3,comm1d,ierr)
        call MPI_Send(zc,dumnelem,MPI_DOUBLE_PRECISION,
     &   0,5*(myid-1)+4,comm1d,ierr)

        call MPI_Send(label,dumnelem,MPI_INTEGER,
     &   0,5*(myid-1)+5,comm1d,ierr)

        else
! Processor 0 receving the data and writing to a file
        open(25,file="label3d_parallel.dat")

        write(*,*) "writing the data of the processor:",myid

        do ielem = 1,dumnelem
       if(zc(ielem)>1.98e-3) write (25,201)xc(ielem),yc(ielem),zc(ielem)
     &   ,label(ielem)
        end do
        
        deallocate(xc,yc,zc,label)
        
        do iproc=1,nprocs-1 

        call MPI_RECV(dumimin,1,MPI_INTEGER,iproc,5*(iproc-1)+1
     &   ,comm1d,status,ierr)

        allocate(xc(1:dumimin))
        allocate(yc(1:dumimin))
        allocate(zc(1:dumimin))
        allocate(label(1:dumimin))
     
        call MPI_Recv(xc,dumimin,MPI_DOUBLE_PRECISION
     &  ,iproc,5*(iproc-1)+2,comm1d,status,ierr)
        call MPI_Recv(yc,dumimin,MPI_DOUBLE_PRECISION
     &  ,iproc,5*(iproc-1)+3,comm1d,status,ierr)
        call MPI_Recv(zc,dumimin,MPI_DOUBLE_PRECISION
     &  ,iproc,5*(iproc-1)+4,comm1d,status,ierr)

        call MPI_Recv(label,dumimin,MPI_INTEGER
     &  ,iproc,5*(iproc-1)+5,comm1d,status,ierr)
     
        do ielem = 1,dumimin
       if(zc(ielem)>1.98e-3) write (25,201)xc(ielem),yc(ielem),zc(ielem)
     &   ,label(ielem)
        end do
        
        deallocate(xc,yc,zc,label)
        write(*,*) "End of writing the data of the processor:",iproc
        enddo
        close(25)
        endif
!---------------------------------------------------------------
!     End of writing label of the mesh
!---------------------------------------------------------------

! Creating new grid points
        dumnnod=0
        dumimin=minval(elem(locelem(1:dumnelem),:))
        dumimax=maxval(elem(locelem(1:dumnelem),:))
        allocate(mapold2new(1:dumimax))
         call mapping3d(elem,nelem,locelem(1:dumnelem),dumnelem
     &   ,mapold2new,dumimax,dumnnod)

!Creating nodes for submesh 
        allocate(dumcoord(1:dumnnod,1:3))
        dumcoord=0
        genflag=0
        do inod=1,dumimax
        if(mapold2new(inod).ne.0) then
        genflag=genflag+1
        dumcoord(genflag,:)=coord(inod,:)
        endif
        enddo !inod

   
        
        if(dumnnod.ne.genflag) then
            write(*,*) "dumnnod not equal to genflag. Some problem"
            write(*,*) "dumnnod:",dumnnod
            write(*,*) "genflag:",genflag
            write(*,*) "Exiting the program"
            call exit(1)
        endif
                
!Creating elements for submesh
        allocate(dumelem(1:dumnelem,1:8))
        do ielem=1,dumnelem       

        ipnt1=elem(locelem(ielem),1);
        ipnt2=elem(locelem(ielem),2);
        ipnt3=elem(locelem(ielem),3);
        ipnt4=elem(locelem(ielem),4);
        ipnt5=elem(locelem(ielem),5);
        ipnt6=elem(locelem(ielem),6);
        ipnt7=elem(locelem(ielem),7);
        ipnt8=elem(locelem(ielem),8);
        
        dumelem(ielem,1)=mapold2new(ipnt1);
        dumelem(ielem,2)=mapold2new(ipnt2);
        dumelem(ielem,3)=mapold2new(ipnt3);
        dumelem(ielem,4)=mapold2new(ipnt4);
        dumelem(ielem,5)=mapold2new(ipnt5);
        dumelem(ielem,6)=mapold2new(ipnt6);
        dumelem(ielem,7)=mapold2new(ipnt7);
        dumelem(ielem,8)=mapold2new(ipnt8);
        enddo !ielem

!Dividing the postprocessing points into different processors
        allocate(dumppnod(1:pptot))
        dumppnod=0
        dumpptot=0

        do inod=1,pptot

        ipnt1=ppnod(inod)
        xmin=coord(ipnt1,1)
        
        if(xmin>xleft.and.xmin<xright)  then
        dumpptot=dumpptot+1;    dumppnod(dumpptot)=mapold2new(ipnt1)
        endif

        enddo

! Finding nodes in submesh for left side receiving nodes
        if(dumnnodlr==0) then
        dumnnodlr=1;    dumlnodr(1)=dumelem(1,1)
        endif
        
         call intrmvdup(dumlnodr(1:dumnnodlr),dumnnodlr,bndcnt)
        allocate(lnodr(1:bndcnt))
        lnodrcnt=bndcnt
        do inod=1,lnodrcnt
        lnodr(inod)=mapold2new(dumlnodr(inod))
        enddo

! Finding nodes in submesh for left side sending nodes
        if(dumnnodls==0) then
        dumnnodls=1;     dumlnods(1)=dumelem(1,2)
        endif

        call intrmvdup(dumlnods(1:dumnnodls),dumnnodls,bndcnt)
        allocate(lnods(1:bndcnt))
        lnodscnt=bndcnt
        do inod=1,lnodscnt
        lnods(inod)=mapold2new(dumlnods(inod))
        enddo

! Finding nodes in submesh for Right side receiving nodes
        if(dumnnodrr==0) then
        dumnnodrr=1;    dumrnodr(1)=dumelem(1,1)
        endif

        call intrmvdup(dumrnodr(1:dumnnodrr),dumnnodrr,bndcnt)
        allocate(rnodr(1:bndcnt))
        rnodrcnt=bndcnt
        do inod=1,rnodrcnt
        rnodr(inod)=mapold2new(dumrnodr(inod))
        enddo

! Finding nodes in submesh for right side sending nodes
        if(dumnnodrs==0) then
        dumnnodrs=1;     dumrnods(1)=dumelem(1,2)
        endif

        call intrmvdup(dumrnods(1:dumnnodrs),dumnnodrs,bndcnt)
        allocate(rnods(1:bndcnt))
        rnodscnt=bndcnt
        do inod=1,rnodscnt
        rnods(inod)=mapold2new(dumrnods(inod))
        enddo

! dividing the faces into different processors
        dumnface=0
        dumface=0
        do iface=1,nface
        ipnt1=int(face(iface,1))
        ipnt2=int(face(iface,2))
        ipnt3=int(face(iface,3))
        ipnt4=int(face(iface,4))
        
        xmin=min(coord(ipnt1,1),coord(ipnt2,1))
        xmin=min(xmin,coord(ipnt3,1))
        xmin=min(xmin,coord(ipnt4,1))
        
        xmax=max(coord(ipnt1,1),coord(ipnt2,1))
        xmax=max(xmax,coord(ipnt3,1))
        xmax=max(xmax,coord(ipnt4,1))
 
        if((xmin>xleft.and.xmin<xright).or.(xmax>xleft.and.xmax<xright))
     &   then
        dumnface=dumnface+1
!        write(*,*) iprocx,iface,nface,"enterd",dumnface
        dumface(dumnface,:)=face(iface,:)
        dumface(dumnface,1)=dble(mapold2new(ipnt1))
        dumface(dumnface,2)=dble(mapold2new(ipnt2))
        dumface(dumnface,3)=dble(mapold2new(ipnt3))
        dumface(dumnface,4)=dble(mapold2new(ipnt4))
        endif
        enddo !iface

        deallocate(elem,coord,face,ppnod)
        
        allocate(elem(1:dumnelem,1:8))
        allocate(coord(1:dumnnod,1:3))
        allocate(face(1:dumnface,1:9))
        allocate(ppnod(1:dumpptot))
        
        nelem=dumnelem;  nnod=dumnnod   
        nface=dumnface;  pptot=dumpptot 
        
        elem=dumelem(1:nelem,:);    coord=dumcoord(1:nnod,:);
        face=dumface(1:nface,:);    ppnod=dumppnod(1:pptot)
        
        deallocate(dumlnodr,dumlnods,dumrnodr,dumrnods,mapold2new)
        deallocate(dumelem,dumppnod,dumcoord,dumface,locelem)

        etime=MPI_WTIME()
                
        write(*,*) myid,"Nelem     :",nelem
        write(*,*) myid,"nodes     :",nnod
        write(*,*) myid,"nfaces    :",nface
        write(*,*) myid,"pptot     :",pptot
        write(*,*) myid,"lnodrcnt  :",lnodrcnt
        write(*,*) myid,"lnodscnt  :",lnodscnt
        write(*,*) myid,"rnodrcnt  :",rnodrcnt
        write(*,*) myid,"rnodscnt  :",rnodscnt
        write(*,*) myid,"com. time :",etime-stime

!---------------------------------------------------------------
!     Gathering all the geometry information
!---------------------------------------------------------------
        allocate(locnelem(1:nprocs),locnnod(1:nprocs))
        allocate(locnface(1:nprocs),locpptot(1:nprocs))
        
        call MPI_ALLGATHER(nelem,1,MPI_INTEGER,locnelem,1,MPI_INTEGER
     &  ,MPI_COMM_WORLD,ierr)
        call MPI_ALLGATHER(nnod,1,MPI_INTEGER,locnnod,1,MPI_INTEGER
     &  ,MPI_COMM_WORLD,ierr)
        call MPI_ALLGATHER(nface,1,MPI_INTEGER,locnface,1,MPI_INTEGER
     &  ,MPI_COMM_WORLD,ierr)
        call MPI_ALLGATHER(pptot,1,MPI_INTEGER,locpptot,1,MPI_INTEGER
     &  ,MPI_COMM_WORLD,ierr)

         if(myid==0) then
           write(*,*) "Total elements :",sum(locnelem)
           write(*,*) "Total nodes    :",sum(locnnod)
           write(*,*) "Total faces    :",sum(locnface)
           write(*,*) "Total pptot    :",sum(locpptot)
         endif
        
        write(*,*) myid,"End of grid division"

  201  FORMAT(E,E,E,I)         
        end

