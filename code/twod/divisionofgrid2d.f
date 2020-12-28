        subroutine divisionofgrid2d()
! Subroutine to division the mesh into sub-meshes for parallel code 
        use variables
        
        IMPLICIT NONE
        DOUBLE PRECISION coordmin(spcdim),coordmax(spcdim)
        DOUBLE PRECISION dlen,dmin,dmax,coordlen(spcdim),fact
        INTEGER axis,ispc,ielem,BUFSIZE,maxb,bndcnt,inod,maxf
        INTEGER flage,flagrs,flagrr,flagls,flaglr,genflag,flagf
        INTEGER ipnt1,ipnt2,ipnt3,ipnt4,iface
        INTEGER,dimension (:), allocatable :: dumlnodr,dumlnods
        INTEGER,dimension (:), allocatable :: dumrnodr,dumrnods,dumppnod
        INTEGER,dimension (:), allocatable :: locelem,mapold2new
        INTEGER,dimension (:,:), allocatable :: dumelem
        DOUBLE PRECISION, dimension (:,:), allocatable :: dumcoord
        DOUBLE PRECISION, dimension (:,:), allocatable :: dumface
        INTEGER dum,dumnnod
        DOUBLE PRECISION ix,iy,dumlab,ipmax,ipmin,dum1,dum2,xmin,xmax
        INTEGER thefile,ptype,dum3
        INTEGER(kind=MPI_OFFSET_KIND) disp
        character(len=20) :: labname
        character(len=100) :: line
        
        maxb=nelem*2/nprocs
        maxf=nface*2/nprocs
        
        allocate(dumlnodr(1:maxb),dumlnods(1:maxb))
        allocate(dumrnodr(1:maxb),dumrnods(1:maxb))
        allocate(locelem(1:maxb),dumface(1:maxf,1:6))
        allocate(dumppnod(1:pptot))
        
        locelem(:)=0;
        dumlnodr(:)=0;   dumlnods(:)=0
        dumrnodr(:)=0;   dumrnods(:)=0
        
c        call MPI_CART_COORDS(comm1d,myid,1,proccoord,ierr)
c!        Write(*,*) "Hello Check this", myid,proccoord
c        if(proccoord(1)==0) left=MPI_PROC_NULL
c        if(proccoord(1)==nprocs-1) right=MPI_PROC_NULL
!-----------------------------------------------------
! Finding the axis of partiotioning the mesh
!-----------------------------------------------------
        do ispc=1,spcdim
        coordmin(ispc)=minval(coord(:,ispc));
        coordmax(ispc)=maxval(coord(:,ispc))
        enddo
        coordlen=coordmax-coordmin
        
        do ispc=1,spcdim
            if(maxval(coordlen)==coordlen(ispc)) then
             axis=ispc;       
            endif
        enddo
        
        dlen=coordlen(axis);
        dmin=coordmin(axis)-0.01*dlen; dmax=coordmax(axis)+0.01*dlen
        
! Length of geometry for each processor. Factor fact can be changed to 
! get proper division of geometry into multiple grids
        fact=1.023
        dlen=fact*dlen/dble(nprocs)
!-----------------------------------------------------
! Labeling the mesh points with the processor number
!-----------------------------------------------------
        allocate(label(1:nelem))
        label(:)=1000

        flage=0
        flagrs=0;    flagrr=0
        flagls=0;    flaglr=0
        
        do ielem=1,nelem
        ipmin=minval(coord(elem(ielem,:),axis))
        ipmax=maxval(coord(elem(ielem,:),axis))
        
        
        if(ipmin>dmin+myid*dlen.and.ipmax<dmin+(myid+1)*dlen) then
        
        flage=flage+1;       locelem(flage)=ielem
        
        label(ielem)=2*(myid)
        
        elseif(ipmin<dmin+myid*dlen.and.ipmax>dmin+myid*dlen) then
        
        flage=flage+1;        locelem(flage)=ielem
        
        ipnt1=elem(ielem,1); ipnt2=elem(ielem,2)
        ipnt3=elem(ielem,3); ipnt4=elem(ielem,4)

        flaglr=flaglr+1; dumlnodr(flaglr)=ipnt1
        flaglr=flaglr+1; dumlnodr(flaglr)=ipnt4

        flagls=flagls+1; dumlnods(flagls)=ipnt2
        flagls=flagls+1; dumlnods(flagls)=ipnt3
        
        label(ielem)=2*(myid)-1
        
        elseif(ipmin<dmin+(myid+1)*dlen.and.
     &    ipmax>dmin+(myid+1)*dlen) then
        
        flage=flage+1;        locelem(flage)=ielem
        
        ipnt1=elem(ielem,1); ipnt2=elem(ielem,2)
        ipnt3=elem(ielem,3); ipnt4=elem(ielem,4)

        flagrr=flagrr+1; dumrnodr(flagrr)=ipnt2
        flagrr=flagrr+1; dumrnodr(flagrr)=ipnt3

        flagrs=flagrs+1; dumrnods(flagrs)=ipnt1
        flagrs=flagrs+1; dumrnods(flagrs)=ipnt4
        
        label(ielem)=2*(myid)+1
        
        endif

        enddo

!mapping old nodes to new nodes of submesh
        
        dum=maxval(elem(locelem(1:flage),:))
        allocate(mapold2new(1:dum))
        call mapping(elem,nelem,locelem,flage,mapold2new,dum,dumnnod)
!        write(*,*) mapold2new
!Creating elements for submesh

        allocate(dumelem(1:flage,1:4),dumcoord(1:dumnnod,1:2))
        do ielem=1,flage        
        ipnt1=elem(locelem(ielem),1);  ipnt2=elem(locelem(ielem),2);
        ipnt3=elem(locelem(ielem),3);  ipnt4=elem(locelem(ielem),4);
        
        
        dumelem(ielem,1)=mapold2new(ipnt1);
        dumelem(ielem,2)=mapold2new(ipnt2);
        dumelem(ielem,3)=mapold2new(ipnt3);
        dumelem(ielem,4)=mapold2new(ipnt4);
        enddo
        
!Creating nodes for submesh 
        genflag=0
        do inod=1,dum
        if(mapold2new(inod).ne.0) then
        genflag=genflag+1
        dumcoord(genflag,:)=coord(inod,:)
        endif
        enddo

! Finding nodes in submesh for right side receiving nodes
        if(flagrr==0) then
        flagrr=1
        dumrnodr(1)=elem(locelem(flage),2)
        endif

        call intrmvdup(dumrnodr(1:flagrr),flagrr,bndcnt)
        
        allocate(rnodr(1:bndcnt))
        rnodrcnt=bndcnt
        do inod=1,rnodrcnt
        rnodr(inod)=mapold2new(dumrnodr(inod))
        enddo

c! Finding nodes in submesh for right side sending nodes
        if(flagrs==0) then
        flagrs=1
        dumrnods(1)=elem(locelem(flage),1)
        endif

        call intrmvdup(dumrnods(1:flagrs),flagrs,bndcnt)
        allocate(rnods(1:bndcnt))
        rnodscnt=bndcnt
        do inod=1,rnodscnt
        rnods(inod)=mapold2new(dumrnods(inod))
        enddo
      
! Finding nodes in submesh for left side receiving nodes
        if(flaglr==0) then
        flaglr=1
        dumlnodr(1)=elem(locelem(1),1)
        endif

        call intrmvdup(dumlnodr(1:flaglr),flaglr,bndcnt)
        allocate(lnodr(1:bndcnt))
        lnodrcnt=bndcnt
        do inod=1,lnodrcnt
        lnodr(inod)=mapold2new(dumlnodr(inod))
        enddo

! Finding nodes in submesh for left side sending nodes
        if(flagls==0) then
        flagls=1
        dumlnods(1)=elem(locelem(1),2)
        endif

        call intrmvdup(dumlnods(1:flagls),flagls,bndcnt)
        allocate(lnods(1:bndcnt))
        lnodscnt=bndcnt
        do inod=1,lnodscnt
        lnods(inod)=mapold2new(dumlnods(inod))
        enddo

! dividing the faces into different processors
        flagf=0
        do iface=1,nface
        ipnt1=int(face(iface,1))
        ipnt2=int(face(iface,2))
        xmin=min(coord(ipnt1,1),coord(ipnt2,1))
        xmax=max(coord(ipnt1,1),coord(ipnt2,1))
        if((xmin>dmin+myid*dlen.and.xmin<dmin+(myid+1)*dlen).or.
     &     (xmax>dmin+myid*dlen.and.xmax<dmin+(myid+1)*dlen)) then
        flagf=flagf+1
        dumface(flagf,:)=face(iface,:)
        dumface(flagf,1)=dble(mapold2new(ipnt1))
        dumface(flagf,2)=dble(mapold2new(ipnt2))
        endif
        enddo
        
        nface=flagf
        deallocate(face)
        allocate(face(1:nface,1:6))
        face(1:nface,1:6)=dumface(1:nface,1:6)
        deallocate(dumface)

!Dividing the postprocessing points into different processors
        flagf=0
        do inod=1,pptot
        ipnt1=ppnod(inod)
        xmin=coord(ipnt1,1)
        if(xmin>dmin+myid*dlen.and.xmin<dmin+(myid+1)*dlen) then
        flagf=flagf+1
        dumppnod(flagf)=mapold2new(ipnt1)
        endif
        enddo
        
        pptot=flagf
        deallocate(ppnod)
        allocate(ppnod(1:pptot))
        ppnod(1:pptot)=dumppnod(1:pptot)
        deallocate(dumppnod)
!-----------------------------------------------------
! Writing the label of each mesh points to a file
!-----------------------------------------------------
        allocate(locnelem(1:nprocs),locnnod(1:nprocs))
        allocate(locnface(1:nprocs),locpptot(1:nprocs))
        
        call MPI_ALLGATHER(flage,1,MPI_INTEGER,locnelem,1,MPI_INTEGER
     &  ,MPI_COMM_WORLD,ierr)
        call MPI_ALLGATHER(dumnnod,1,MPI_INTEGER,locnnod,1,MPI_INTEGER
     &  ,MPI_COMM_WORLD,ierr)
        call MPI_ALLGATHER(nface,1,MPI_INTEGER,locnface,1,MPI_INTEGER
     &  ,MPI_COMM_WORLD,ierr)
        call MPI_ALLGATHER(pptot,1,MPI_INTEGER,locpptot,1,MPI_INTEGER
     &  ,MPI_COMM_WORLD,ierr)
          
c        all MPI_TYPE_CONTIGUOUS(3,MPI_DOUBLE,ptype,ierr)
c        all MPI_TYPE_COMMIT(ptype,ierr)

c        all MPI_File_open(MPI_COMM_WORLD,"prldivmesh",
c     &  MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,thefile,ierr)
        
c        do ielem=1,flage
        
c        ix=sum(coord(elem(locelem(ielem),:),1))/4.0
c        iy=sum(coord(elem(locelem(ielem),:),2))/4.0
c        dumlab=dble(label(locelem(ielem)))
c        write(*,*) myid, ielem,ix,iy,dumlab
        
c        if(myid==0) then
c            disp=(ielem-1)*8*3*0
c        else
c            disp=sum(locnelem(1:myid))*8*3+(ielem-1)*8*3*0
c        endif
        
c!        call MPI_FILE_SET_VIEW(thefile,4,MPI_INTEGER,MPI_INTEGER
c!     &   , 'native', MPI_INFO_NULL, ierr)
c        all MPI_FILE_WRITE(thefile,1,1,MPI_INTEGER,MPI_STATUS_IGNORE
c     &  , ierr)

c        enddo
c        all MPI_FILE_CLOSE(thefile, ierr)

! Writing data to multiple files
        write (labname, "('label',I3.3,'.dat')") myid

        open(25,file=labname)

        do ielem=1,nelem
        ix=sum(coord(elem(ielem,:),1))/4.0
        iy=sum(coord(elem(ielem,:),2))/4.0
        if(label(ielem).ne.1000) then
         write(25,*) ix,iy,label(ielem)
        endif
        enddo
        close(25)

! Processor 0 reading the parallel files and writing to single file
        if(myid==0) then
         write (labname, "('label',I3.3,'.dat')") myid
        open(25,file=labname,status="old",position="append",
     &     action="write")
        
        do inod=1,nprocs-1
        
        write (labname, "('label',I3.3,'.dat')") inod
        open(26,file=labname,status="old",action="read")

        do ielem=1,locnelem(inod+1)
        read(26,*) dum1,dum2,dum3
        write(25,*) dum1,dum2,dum3
        enddo

        close(26)

        enddo
        
        endif
! End of file operation        
        deallocate(elem,coord)
        
        allocate(elem(1:flage,1:4),coord(1:dumnnod,1:2))
        elem=dumelem
        coord=dumcoord
        nelem=flage
        nnod=dumnnod
        
        deallocate(dumlnodr,dumlnods,dumrnodr,dumrnods,locelem)
        deallocate(mapold2new,dumelem,dumcoord)
        
        write(*,*) myid,nelem,nface,nnod,pptot

        if(myid==0) then
         Write(*,*) "Elements in each processor:",locnelem
         write(*,*) "Total elements:",sum(locnelem)
         Write(*,*) "Nodes in each processor:",locnnod
         write(*,*) "Total nodes:",sum(locnnod)
         Write(*,*) "Faces in each processor:",locnface
         write(*,*) "Total faces:",sum(locnface)
         Write(*,*) "PPtot in each processor:",locpptot
         write(*,*) "Total pptot:",sum(locpptot)
         endif
         
        end

