        subroutine divisionofgrid3d_2d()
! Subroutine to division the mesh into sub-meshes for parallel code 
        use variables
        
        IMPLICIT NONE
        DOUBLE PRECISION coordmin(spcdim),coordmax(spcdim)
        DOUBLE PRECISION dlenx,dminx,dmaxx,coordlen(spcdim),fact
        DOUBLE PRECISION dleny,dminy,dmaxy,xleft,xright
        INTEGER axis,ispc,ielem,maxb,bndcnt,inod,maxf,maxe
        INTEGER genflag
        INTEGER ipnt1,ipnt2,ipnt3,ipnt4,ipnt5,ipnt6,ipnt7,ipnt8,iface

        INTEGER,dimension (:,:), allocatable :: loclnodr,loclnods
        INTEGER,dimension (:,:), allocatable :: locrnodr,locrnods
        INTEGER,dimension (:),   allocatable :: locnnodls,locnnodlr
        INTEGER,dimension (:),   allocatable :: locnnodrs,locnnodrr

        INTEGER,dimension (:,:), allocatable :: locbnodr,locbnods
        INTEGER,dimension (:,:), allocatable :: loctnodr,loctnods
        INTEGER,dimension (:),   allocatable :: locnnodbs,locnnodbr
        INTEGER,dimension (:),   allocatable :: locnnodts,locnnodtr

        INTEGER,dimension (:,:), allocatable :: locelem

        INTEGER dumnelem,dumnnod,dumpptot,dumnface
        INTEGER,dimension (:), allocatable :: dumppnod
        INTEGER,dimension (:,:), allocatable :: dumelem
        DOUBLE PRECISION, dimension (:,:), allocatable :: dumcoord
        DOUBLE PRECISION, dimension (:,:), allocatable :: dumface
        INTEGER,dimension (:), allocatable :: mapold2new,dummy

        INTEGER iproc,iprocx,iprocy
        INTEGER proccoordall(2,nprocs)     
        DOUBLE PRECISION ix,iy,iz,stime,etime
        DOUBLE PRECISION dumlab,ipmaxx,ipminx,dumd1,dumd2,xmin,xmax
        DOUBLE PRECISION ipmaxy,ipminy,ymin,ymax
        INTEGER thefile,ptype,dumimin,dumimax
        INTEGER(kind=MPI_OFFSET_KIND) disp
        character(len=20) :: labname
        character(len=100) :: line
        
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
!        write(*,*) "procoordtomyid:",proccordtomyid(:,0)
        if(myid==0) then         
        
! First processor dividing the mesh into multiple processors         
        call readgrid3d()

! allocating the memory to variables
        maxe=nelem*2/(nprocx*nprocy)
        maxf=nface*2/(nprocx*nprocy)
        maxb=nnod*2/(nprocx*nprocy)
        
       allocate(loclnodr(1:maxb,0:nprocx-1))
       allocate(loclnods(1:maxb,0:nprocx-1))
       allocate(locrnodr(1:maxb,0:nprocx-1))
       allocate(locrnods(1:maxb,0:nprocx-1))
       
       allocate(locnnodls(0:nprocx-1),locnnodlr(0:nprocx-1))
       allocate(locnnodrs(0:nprocx-1),locnnodrr(0:nprocx-1))

       allocate(locnelem(0:nprocx-1))
       allocate(locelem(1:maxb,0:nprocx-1))
       allocate(dumppnod(1:pptot))
       allocate(dumface(1:maxf,1:9))

       loclnodr=0;   loclnods=0;    locrnodr=0;     locrnods=0; 
       locnnodls=0;  locnnodlr=0;   locnnodrs=0;    locnnodrr=0;
       locnelem=0;   locelem=0 
       
       dlenx=maxval(coord(:,1))-minval(coord(:,1))

       dminx=minval(coord(:,1))-0.01*dlenx;
       dmaxx=maxval(coord(:,1))+0.01*dlenx

       fact=1.0377773
       dlenx=fact*dlenx/dble(nprocx);
       write(*,*) "Number of processors in x-direction:",nprocx
       write(*,*) "min x-coord(actual,used)",minval(coord(:,1)),dminx
       write(*,*) "Length of each domain along x-dir:",dlenx
!-----------------------------------------------------
! Labeling the mesh points with the processor number
!-----------------------------------------------------
       allocate(label(1:nelem))
       label=1000
        
! division of elements
       do ielem=1,nelem
        ipminx=minval(coord(elem(ielem,:),1))
        ipmaxx=maxval(coord(elem(ielem,:),1))

        ipminy=minval(coord(elem(ielem,:),2))
        ipmaxy=maxval(coord(elem(ielem,:),2))

!       write(*,*) xleft,xright
        ipnt1=elem(ielem,1); ipnt2=elem(ielem,2)
        ipnt3=elem(ielem,3); ipnt4=elem(ielem,4)
        ipnt5=elem(ielem,5); ipnt6=elem(ielem,6)
        ipnt7=elem(ielem,7); ipnt8=elem(ielem,8) 

       do iprocy=0,nprocy-1
       do iprocx=0,nprocx-1
       
       iproc=proccordtomyid(iprocx,iprocy)
       
       xleft=dminx+iprocx*dlenx;   xright=dminx+(iprocx+1)*dlenx
       ybottom=dminy+iprocy*dleny; ytop=dminy+(iprocy+1)*dleny
               
       if(ipminx>xleft  .and.ipmaxx<xright.and.
     &    ipminy>ybottom.and.ipmaxy<ytop) then
     
        locnelem(iproc)=locnelem(iproc)+1
        locelem(locnelem(iproc),iproc)=ielem
        label(ielem)=iproc
        
       elseif(ipminx<xleft.and.ipmaxx>xleft) then

        locnelem(iprocx)=locnelem(iprocx)+1
        locelem(locnelem(iprocx),iprocx)=ielem

! Left side receviing nodes
        locnnodlr(iprocx)=locnnodlr(iprocx)+1; 
        loclnodr(locnnodlr(iprocx),iprocx)=ipnt1
        
        locnnodlr(iprocx)=locnnodlr(iprocx)+1;
        loclnodr(locnnodlr(iprocx),iprocx)=ipnt4
        
        locnnodlr(iprocx)=locnnodlr(iprocx)+1;   
        loclnodr(locnnodlr(iprocx),iprocx)=ipnt5
        
        locnnodlr(iprocx)=locnnodlr(iprocx)+1;   
        loclnodr(locnnodlr(iprocx),iprocx)=ipnt8

! Left side Sending nodes
        locnnodls(iprocx)=locnnodls(iprocx)+1;   
        loclnods(locnnodls(iprocx),iprocx)=ipnt2
        
        locnnodls(iprocx)=locnnodls(iprocx)+1;   
        loclnods(locnnodls(iprocx),iprocx)=ipnt3
        
        locnnodls(iprocx)=locnnodls(iprocx)+1;   
        loclnods(locnnodls(iprocx),iprocx)=ipnt6
        
        locnnodls(iprocx)=locnnodls(iprocx)+1;   
        loclnods(locnnodls(iprocx),iprocx)=ipnt7

        label(ielem)=2*iprocx-1

       else if(ipminx<xright.and.ipmaxx>xright) then     

        locnelem(iprocx)=locnelem(iprocx)+1
        locelem(locnelem(iprocx),iprocx)=ielem
        
! Right side receviing nodes
        locnnodrr(iprocx)=locnnodrr(iprocx)+1;
        locrnodr(locnnodrr(iprocx),iprocx)=ipnt2
        
        locnnodrr(iprocx)=locnnodrr(iprocx)+1;   
        locrnodr(locnnodrr(iprocx),iprocx)=ipnt3
        
        locnnodrr(iprocx)=locnnodrr(iprocx)+1;
        locrnodr(locnnodrr(iprocx),iprocx)=ipnt6
        
        locnnodrr(iprocx)=locnnodrr(iprocx)+1;   
        locrnodr(locnnodrr(iprocx),iprocx)=ipnt7

! Right side Sending nodes
        locnnodrs(iprocx)=locnnodrs(iprocx)+1;   
        locrnods(locnnodrs(iprocx),iprocx)=ipnt1
        
        locnnodrs(iprocx)=locnnodrs(iprocx)+1;   
        locrnods(locnnodrs(iprocx),iprocx)=ipnt4
        
        locnnodrs(iprocx)=locnnodrs(iprocx)+1;   
        locrnods(locnnodrs(iprocx),iprocx)=ipnt5
        
        locnnodrs(iprocx)=locnnodrs(iprocx)+1;   
        locrnods(locnnodrs(iprocx),iprocx)=ipnt8

        label(ielem)=2*iprocx+1
        
       endif 
       enddo !iprocx
       enddo
       enddo !ielem

!        write(*,*) locnelem,sum(locnelem)

        open(25,file="labelname3d.dat")
        do ielem=1,nelem
        ix=sum(coord(elem(ielem,:),1))/8.0
        iy=sum(coord(elem(ielem,:),2))/8.0
        iz=maxval(coord(elem(ielem,:),3))
        if(iz>1.99e-3) write(25,*) ix,iy,iz,label(ielem)
        enddo
        close(25)
        deallocate(label)

! End of division of elements
!-----------------------------------------------------
! end of Labeling the mesh points with the processor number
!-----------------------------------------------------

!-----------------------------------------------------
! Creating grid for all the processors and sending it
!-----------------------------------------------------
       do iprocx=nprocx-1,0,-1

       xleft=dminx+iprocx*dlenx; xright=dminx+(iprocx+1)*dlenx

!mapping old nodes to new nodes of submesh
        dumnnod=0
        dumimin=minval(elem(locelem(1:locnelem(iprocx),iprocx),:))
        dumimax=maxval(elem(locelem(1:locnelem(iprocx),iprocx),:))
        allocate(mapold2new(1:dumimax))
        call mapping3d(elem,nelem,locelem(1:locnelem(iprocx),iprocx)
     &   ,locnelem(iprocx),mapold2new,dumimax,dumnnod)

!Creating nodes for submesh 
        allocate(dumcoord(1:dumnnod,1:3))
        dumcoord=0
        genflag=0
        do inod=dumimin,dumimax
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
        allocate(dumelem(1:locnelem(iprocx),1:8))
        do ielem=1,locnelem(iprocx)       

        ipnt1=elem(locelem(ielem,iprocx),1);
        ipnt2=elem(locelem(ielem,iprocx),2);
        ipnt3=elem(locelem(ielem,iprocx),3);
        ipnt4=elem(locelem(ielem,iprocx),4);
        ipnt5=elem(locelem(ielem,iprocx),5);
        ipnt6=elem(locelem(ielem,iprocx),6);
        ipnt7=elem(locelem(ielem,iprocx),7);
        ipnt8=elem(locelem(ielem,iprocx),8);
        
        dumelem(ielem,1)=mapold2new(ipnt1);
        dumelem(ielem,2)=mapold2new(ipnt2);
        dumelem(ielem,3)=mapold2new(ipnt3);
        dumelem(ielem,4)=mapold2new(ipnt4);
        dumelem(ielem,5)=mapold2new(ipnt5);
        dumelem(ielem,6)=mapold2new(ipnt6);
        dumelem(ielem,7)=mapold2new(ipnt7);
        dumelem(ielem,8)=mapold2new(ipnt8);
        enddo !ielem
        dumnelem=locnelem(iprocx)

!Dividing the postprocessing points into different processors
        dumppnod=0
        dumpptot=0
        do inod=1,pptot
        ipnt1=ppnod(inod)
        xmin=coord(ipnt1,1)
        if(xmin>xleft.and.xmin<xright)  then
        dumpptot=dumpptot+1
        dumppnod(dumpptot)=mapold2new(ipnt1)
        endif
        enddo

! Finding nodes in submesh for left side receiving nodes
        if(locnnodlr(iprocx)==0) then
        locnnodlr(iprocx)=1
        loclnodr(1,iprocx)=dumelem(1,1)
        endif

        call intrmvdup(loclnodr(1:locnnodlr(iprocx),iprocx)
     &   ,locnnodlr(iprocx),bndcnt)
        allocate(lnodr(1:bndcnt))
        lnodrcnt=bndcnt
        do inod=1,lnodrcnt
        lnodr(inod)=mapold2new(loclnodr(inod,iprocx))
        enddo

! Finding nodes in submesh for left side sending nodes
        if(locnnodls(iprocx)==0) then
        locnnodls(iprocx)=1
        loclnods(1,iprocx)=dumelem(1,2)
        endif

        call intrmvdup(loclnods(1:locnnodls(iprocx),iprocx)
     &  ,locnnodls(iprocx),bndcnt)
        allocate(lnods(1:bndcnt))
        lnodscnt=bndcnt
        do inod=1,lnodscnt
        lnods(inod)=mapold2new(loclnods(inod,iprocx))
        enddo

! Finding nodes in submesh for right side receiving nodes
        if(locnnodrr(iprocx)==0) then
        locnnodrr(iprocx)=1
        locrnodr(1,iprocx)=dumelem(1,1)
        endif

        call intrmvdup(locrnodr(1:locnnodrr(iprocx),iprocx),
     &   locnnodrr(iprocx),bndcnt)
        allocate(rnodr(1:bndcnt))
        rnodrcnt=bndcnt
        do inod=1,rnodrcnt
        rnodr(inod)=mapold2new(locrnodr(inod,iprocx))
        enddo

! Finding nodes in submesh for right side sending nodes
        if(locnnodrs(iprocx)==0) then
        locnnodrs(iprocx)=1
        locrnods(1,iprocx)=dumelem(1,1)
        endif

        call intrmvdup(locrnods(1:locnnodrs(iprocx),iprocx),
     &   locnnodrs(iprocx),bndcnt)
        allocate(rnods(1:bndcnt))
        rnodscnt=bndcnt
        do inod=1,rnodscnt
        rnods(inod)=mapold2new(locrnods(inod,iprocx))
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
        
! sending t processors other than 0
        if(iprocx.ne.0) then
! Sending element information to other processors than 0       
        call MPI_Send(dumnelem,1,MPI_INTEGER,proccordtomyid(iprocx,0)
     &  ,1,comm1d,ierr)

        call MPI_Send(dumelem(1:dumnelem,1:8),dumnelem*8,MPI_INTEGER
     &  ,proccordtomyid(iprocx,0),2,comm1d,ierr) 

! Sending coordinate information to all processors
        call MPI_Send(dumnnod,1,MPI_INTEGER,proccordtomyid(iprocx,0)
     &  ,3,comm1d,ierr)

        call MPI_Send(dumcoord(1:dumnnod,1:3),dumnnod*3,
     & MPI_DOUBLE_PRECISION,proccordtomyid(iprocx,0),4,comm1d,ierr)

! Sending face information to all processors
        call MPI_Send(dumnface,1,MPI_INTEGER,proccordtomyid(iprocx,0)
     &  ,5,comm1d,ierr)

       call MPI_Send(dumface(1:dumnface,1:9),dumnface*9,
     & MPI_DOUBLE_PRECISION,proccordtomyid(iprocx,0),6,comm1d,ierr)

! Sending post processing point information to all processors
        call MPI_Send(dumpptot,1,MPI_INTEGER,proccordtomyid(iprocx,0),
     &  7,comm1d,ierr)

        call MPI_Send(dumppnod(1:dumpptot),dumpptot,
     & MPI_INTEGER,proccordtomyid(iprocx,0),8,comm1d,ierr)     

! Sending right side exchange nodes information to all processors
        call MPI_Send(rnodrcnt,1,MPI_INTEGER,proccordtomyid(iprocx,0)
     &  ,9,comm1d,ierr)

        call MPI_Send(rnodr(1:rnodrcnt),rnodrcnt,MPI_INTEGER
     &   ,proccordtomyid(iprocx,0),10,comm1d,ierr)

        call MPI_Send(rnodscnt,1,MPI_INTEGER,proccordtomyid(iprocx,0)
     &  ,11,comm1d,ierr)

        call MPI_Send(rnods(1:rnodscnt),rnodscnt,MPI_INTEGER
     &   ,proccordtomyid(iprocx,0),12,comm1d,ierr)

! Sending left side exchange nodes information to all processors
        call MPI_Send(lnodrcnt,1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,0),13,comm1d,ierr)

        call MPI_Send(lnodr(1:lnodrcnt),lnodrcnt,MPI_INTEGER
     &   ,proccordtomyid(iprocx,0),14,comm1d,ierr)

        call MPI_Send(lnodscnt,1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,0),15,comm1d,ierr)

        call MPI_Send(lnods(1:lnodscnt),lnodscnt,MPI_INTEGER
     &   ,proccordtomyid(iprocx,0),16,comm1d,ierr)

        deallocate(lnodr,lnods)    
        deallocate(rnodr,rnods)
! end of processors other than 0
        else 

! Processor 0
        deallocate(elem)
        allocate(elem(1:locnelem(0),1:8))
        elem=dumelem(1:locnelem(0),1:8)
        nelem=locnelem(0)
        
        deallocate(coord)
        allocate(coord(1:dumnnod,1:3))
        coord=dumcoord(1:dumnnod,1:3)
        nnod=dumnnod
        
        deallocate(face)
        allocate(face(1:dumnface,1:9))
        face=dumface
        nface=dumnface

        deallocate(ppnod)
        allocate(ppnod(1:dumpptot))
        ppnod=dumppnod
        pptot=dumpptot

       etime=MPI_WTIME()
    
        write(*,*) "myid      :",myid
        write(*,*) "Nelem     :",nelem
        write(*,*) "nodes     :",nnod
        write(*,*) "nfaces    :",nface
        write(*,*) "pptot     :",pptot
        write(*,*) "lnodrcnt  :",lnodrcnt
        write(*,*) "lnodscnt  :",lnodscnt
        write(*,*) "rnodrcnt  :",rnodrcnt
        write(*,*) "rnodscnt  :",rnodscnt
        write(*,*) "com. time :",etime-stime
        write(*,*) "--------------------------------------"
        endif

        deallocate(mapold2new,dumcoord,dumelem) 

       enddo !iprocx
!-----------------------------------------------------
! Creating grid for all the processors and sending it except 0
!-----------------------------------------------------

       deallocate(loclnodr,loclnods,locrnodr,locrnods)
       deallocate(locnnodls,locnnodlr,locnnodrs,locnnodrr)
       deallocate(locnelem,locelem,dumppnod,dumface)
              
       else 
!-------------------------------------------------------------------
! other than processor 0 receviing the information
!-------------------------------------------------------------------
        call MPI_RECV(nelem,1,MPI_INTEGER,0,1,comm1d,status,ierr)
        
        allocate(elem(1:nelem,8))
        call MPI_RECV(elem(1:nelem,1:8),nelem*8,MPI_INTEGER,0,2,
     &  comm1d,status,ierr) 

! Sending coordinate information to all processors
        call MPI_RECV(nnod,1,MPI_INTEGER,0,3,comm1d,status,ierr)
        
        allocate(coord(1:nnod,1:3))
        call MPI_RECV(coord(1:nnod,1:3),nnod*3,MPI_DOUBLE_PRECISION
     & ,0,4,comm1d,status,ierr)

! Sending face information to all processors
       call MPI_RECV(nface,1,MPI_INTEGER,0,5,comm1d,status,ierr)

       allocate(face(1:nface,1:9))
       call MPI_RECV(face(1:nface,1:9),nface*9,MPI_DOUBLE_PRECISION
     & ,0,6,comm1d,status,ierr)

! Sending post processing point information to all processors
        call MPI_RECV(pptot,1,MPI_INTEGER,0,7,comm1d,status,ierr)
        
        allocate(ppnod(1:pptot))
        call MPI_RECV(ppnod(1:pptot),pptot,MPI_INTEGER,0,8,comm1d,
     &   status,ierr)     

! Sending right side exchange nodes information to all processors
        call MPI_RECV(rnodrcnt,1,MPI_INTEGER,0,9,comm1d,status,ierr)
        
        allocate(rnodr(1:rnodrcnt))
        call MPI_RECV(rnodr(1:rnodrcnt),rnodrcnt,MPI_INTEGER
     &   ,0,10,comm1d,status,ierr)

        call MPI_RECV(rnodscnt,1,MPI_INTEGER,0,11,comm1d,status,ierr)

        allocate(rnods(1:rnodscnt))
        call MPI_RECV(rnods(1:rnodscnt),rnodscnt,MPI_INTEGER
     &   ,0,12,comm1d,status,ierr)

! Sending left side exchange nodes information to all processors
        call MPI_RECV(lnodrcnt,1,MPI_INTEGER,0,13,comm1d,status,ierr)

        allocate(lnodr(1:lnodrcnt))
        call MPI_RECV(lnodr(1:lnodrcnt),lnodrcnt,MPI_INTEGER
     &   ,0,14,comm1d,status,ierr)

        call MPI_RECV(lnodscnt,1,MPI_INTEGER,0,15,comm1d,status,ierr)

        allocate(lnods(1:lnodscnt))
        call MPI_RECV(lnods(1:lnodscnt),lnodscnt,MPI_INTEGER
     &   ,0,16,comm1d,status,ierr)

       etime=MPI_WTIME()
    
        write(*,*) "myid      :",myid
        write(*,*) "Nelem     :",nelem
        write(*,*) "nodes     :",nnod
        write(*,*) "nfaces    :",nface
        write(*,*) "pptot     :",pptot
        write(*,*) "lnodrcnt  :",lnodrcnt
        write(*,*) "lnodscnt  :",lnodscnt
        write(*,*) "rnodrcnt  :",rnodrcnt
        write(*,*) "rnodscnt  :",rnodscnt
        write(*,*) "com. time :",etime-stime
        write(*,*) "-------------------------------------"

       endif 
!-------------------------------------------------------------------
! end of other than processor 0 receviing the information
!-------------------------------------------------------------------
        allocate(locnelem(1:nprocx),locnnod(1:nprocx))
        allocate(locnface(1:nprocx),locpptot(1:nprocx))
        
        call MPI_ALLGATHER(nelem,1,MPI_INTEGER,locnelem,1,MPI_INTEGER
     &  ,MPI_COMM_WORLD,ierr)
        call MPI_ALLGATHER(nnod,1,MPI_INTEGER,locnnod,1,MPI_INTEGER
     &  ,MPI_COMM_WORLD,ierr)
        call MPI_ALLGATHER(nface,1,MPI_INTEGER,locnface,1,MPI_INTEGER
     &  ,MPI_COMM_WORLD,ierr)
        call MPI_ALLGATHER(pptot,1,MPI_INTEGER,locpptot,1,MPI_INTEGER
     &  ,MPI_COMM_WORLD,ierr)

         if(myid==0) then
           write(*,*) "Total elements:",sum(locnelem)
           write(*,*) "Total nodes:",sum(locnnod)
           write(*,*) "Total faces:",sum(locnface)
           write(*,*) "Total pptot:",sum(locpptot)
          
        open(100,file="comp_time.dat",ACCESS="Append")        
        write(100,*) sum(locnnod),"Preprocessing",nprocs,etime-stime
        close(100)
        write(*,*) "End of grid division"
         endif

!  201  FORMAT(E,E,E,I)         
        end

