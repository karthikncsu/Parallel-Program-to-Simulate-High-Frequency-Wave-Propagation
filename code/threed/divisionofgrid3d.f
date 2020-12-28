        subroutine divisionofgrid3d()
! Subroutine to division the mesh into sub-meshes for parallel code 
        use variables
        
        IMPLICIT NONE
        DOUBLE PRECISION coordmin(spcdim),coordmax(spcdim)
        DOUBLE PRECISION dlenx,dminx,dmaxx,coordlen(spcdim),fact
        DOUBLE PRECISION dleny,dminy,dmaxy
        INTEGER axis,ispc,ielem,BUFSIZE,maxb,bndcnt,inod,maxf
        INTEGER genflag,flagf
        INTEGER flage,flagls,flaglr,flagrs,flagrr,flagbs,flagbr
        INTEGER flagts,flagtr,flaglbs,flaglbr,flagrbs,flagrbr
        INTEGER flaglts,flagltr,flagrts,flagrtr
        INTEGER ipnt1,ipnt2,ipnt3,ipnt4,ipnt5,ipnt6,ipnt7,ipnt8,iface
        INTEGER,dimension (:), allocatable :: dumlnodr,dumlnods
        INTEGER,dimension (:), allocatable :: dumrnodr,dumrnods
        INTEGER,dimension (:), allocatable :: dumbnodr,dumbnods
        INTEGER,dimension (:), allocatable :: dumtnodr,dumtnods
        INTEGER,dimension (:), allocatable :: dumlbnodr,dumlbnods
        INTEGER,dimension (:), allocatable :: dumrbnodr,dumrbnods
        INTEGER,dimension (:), allocatable :: dumltnodr,dumltnods
        INTEGER,dimension (:), allocatable :: dumrtnodr,dumrtnods
        INTEGER,dimension (:), allocatable :: dumppnod,locelem
        INTEGER,dimension (:), allocatable :: mapold2new,dummy
        INTEGER,dimension (:,:), allocatable :: dumelem
        DOUBLE PRECISION, dimension (:,:), allocatable :: dumcoord
        DOUBLE PRECISION, dimension (:,:), allocatable :: dumface
        INTEGER dum,dumnnod,dumpptot,dumnface,iproc,iprocx,iprocy
        INTEGER proccoordall(2,nprocs)     
        DOUBLE PRECISION ix,iy,iz
        DOUBLE PRECISION dumlab,ipmaxx,ipminx,dum1,dum2,xmin,xmax
        DOUBLE PRECISION ipmaxy,ipminy,ymin,ymax
        INTEGER thefile,ptype,dum3
        INTEGER(kind=MPI_OFFSET_KIND) disp
        character(len=20) :: labname
        character(len=100) :: line
        
        call MPI_CART_COORDS(comm1d,myid,2,proccoord,ierr)
!        Write(*,*) myid,"Processor coordiantes",proccoord
        
        call MPI_ALLGATHER(proccoord,2,MPI_INTEGER,proccoordall,2,
     &   MPI_INTEGER,comm1d,ierr)
     
!         write(*,*) myid,proccoordall
        allocate(proccordtomyid(0:nprocx-1,0:nprocy-1))
        proccordtomyid=0
         
         do iproc=1,nprocs
         ipnt1=proccoordall(1,iproc)
         ipnt2=proccoordall(2,iproc)

         proccordtomyid(ipnt1,ipnt2)=iproc-1

         enddo
!         write(*,*) rtnods
!         write(*,*) proccordtomyid
!         write(*,*) myid,left,right
         
        if(myid==0) then
         call readgrid3d()

! First processor dividing the mesh into multiple processors         
        maxb=nelem*2/(nprocx*nprocy)
        maxf=nface*2/(nprocx*nprocy)
        
        allocate(dumlnodr(1:maxb),dumlnods(1:maxb))
        allocate(dumrnodr(1:maxb),dumrnods(1:maxb))
        allocate(dumbnodr(1:maxb),dumbnods(1:maxb))
        allocate(dumtnodr(1:maxb),dumtnods(1:maxb))
        
        allocate(dumlbnodr(1:maxf),dumlbnods(1:maxf))
        allocate(dumrbnodr(1:maxf),dumrbnods(1:maxf))
        allocate(dumltnodr(1:maxf),dumltnods(1:maxf))
        allocate(dumrtnodr(1:maxf),dumrtnods(1:maxf))

        allocate(locelem(1:maxb),dumface(1:maxf,1:9))
        allocate(dumppnod(1:pptot))

!-----------------------------------------------------
! Finding the axis of partiotioning the mesh
!-----------------------------------------------------
        coordmin=0; coordmax=0; coordlen=0
        
        do ispc=1,spcdim
        coordmin(ispc)=minval(coord(:,ispc))
        coordmax(ispc)=maxval(coord(:,ispc))
        enddo
        coordlen=coordmax-coordmin

        do ispc=1,spcdim
            if(maxval(coordlen)==coordlen(ispc)) then
             axis=ispc;       
            endif
        enddo

        if(axis==3) then
        Write(*,*) "Longest dimension of the geometry is along z-axis"
        Write(*,*) "Redraw the geometry with longest dimension along
     &  x-axis for parallel program to work"
        call exit(5)
        endif
        
        dlenx=coordlen(1);            dleny=coordlen(2) 
        dminx=coordmin(1)-0.01*dlenx; dmaxx=coordmax(1)+0.01*dlenx
        dminy=coordmin(2)-0.01*dleny; dmaxy=coordmax(2)+0.01*dleny
        
! Length of geometry for each processor. Factor fact can be changed to 
! get proper division of geometry into multiple grids
        fact=1.033
        dlenx=fact*dlenx/dble(nprocx);  dleny=fact*dleny/dble(nprocy)
        write(*,*)  dlenx,dleny
!-----------------------------------------------------
! Labeling the mesh points with the processor number
!-----------------------------------------------------
        allocate(label(1:nelem))
        label=1000
        
        do iprocy=nprocy-1,0,-1
        do iprocx=nprocx-1,0,-1
        
        dumlnodr=0;  dumlnods=0;  dumrnodr=0;  dumrnods=0
        dumbnodr=0;  dumbnods=0;  dumtnodr=0;  dumtnods=0
        dumlbnodr=0; dumlbnods=0; dumrbnodr=0; dumrbnods=0
        dumltnodr=0; dumltnods=0; dumrtnodr=0; dumrtnods=0

        locelem=0;  dumface=0;  dumppnod=0; flage=0
        
        flagls=0; flaglr=0; flagrs=0; flagrr=0
        flagbs=0; flagbr=0; flagts=0; flagtr=0
        
        flaglbs=0; flaglbr=0; flagrbs=0; flagrbr=0
        flaglts=0; flagltr=0; flagrts=0; flagrtr=0

        do ielem=1,nelem
        ipminx=minval(coord(elem(ielem,:),1))
        ipmaxx=maxval(coord(elem(ielem,:),1))

        ipminy=minval(coord(elem(ielem,:),2))
        ipmaxy=maxval(coord(elem(ielem,:),2))
        
        ipnt1=elem(ielem,1); ipnt2=elem(ielem,2)
        ipnt3=elem(ielem,3); ipnt4=elem(ielem,4)
        ipnt5=elem(ielem,5); ipnt6=elem(ielem,6)
        ipnt7=elem(ielem,7); ipnt8=elem(ielem,8)
        
        iproc=iprocy*nprocx+iprocx+1

!        write(*,*) iproc,flagrr(iproc),maxb
        
        if(ipminx>dminx+iprocx*dlenx.and.ipmaxx<dminx+(iprocx+1)*dlenx
     &  .and.ipminy>dminy+iprocy*dleny
     &  .and.ipmaxy<dminy+(iprocy+1)*dleny) then
        
        flage=flage+1; locelem(flage)=ielem
        
        label(ielem)=proccordtomyid(iprocx,iprocy)
        
        elseif(ipminx<dminx+iprocx*dlenx.and.ipmaxx>dminx+iprocx*dlenx
     &  .and.ipminy>dminy+iprocy*dleny
     &  .and.ipmaxy<dminy+(iprocy+1)*dleny) then
        
        flage=flage+1;  locelem(flage)=ielem

        flaglr=flaglr+1;    dumlnodr(flaglr)=ipnt1
        flaglr=flaglr+1;    dumlnodr(flaglr)=ipnt4
        flaglr=flaglr+1;    dumlnodr(flaglr)=ipnt5
        flaglr=flaglr+1;    dumlnodr(flaglr)=ipnt8

        flagls=flagls+1;    dumlnods(flagls)=ipnt2
        flagls=flagls+1;    dumlnods(flagls)=ipnt3
        flagls=flagls+1;    dumlnods(flagls)=ipnt6
        flagls=flagls+1;    dumlnods(flagls)=ipnt7
        
        label(ielem)=10
        
        elseif(ipminx<dminx+(iprocx+1)*dlenx.and.ipmaxx>dminx+(iprocx+1)
     &   *dlenx.and.ipminy>dminy+iprocy*dleny
     &        .and.ipmaxy<dminy+(iprocy+1)*dleny) then
        
        flage=flage+1;  locelem(flage)=ielem

        flagrr=flagrr+1;     dumrnodr(flagrr)=ipnt2
        flagrr=flagrr+1;     dumrnodr(flagrr)=ipnt3
        flagrr=flagrr+1;     dumrnodr(flagrr)=ipnt6
        flagrr=flagrr+1;     dumrnodr(flagrr)=ipnt7

        flagrs=flagrs+1;     dumrnods(flagrs)=ipnt1
        flagrs=flagrs+1;     dumrnods(flagrs)=ipnt4
        flagrs=flagrs+1;     dumrnods(flagrs)=ipnt5
        flagrs=flagrs+1;     dumrnods(flagrs)=ipnt8
        
        label(ielem)=10

        elseif(ipminy<dminy+iprocy*dleny.and.ipmaxy>dminy+iprocy*dleny
     &  .and.ipminx>dminx+iprocx*dlenx
     &  .and.ipmaxx<dminx+(iprocx+1)*dlenx) then
        
        flage=flage+1;  locelem(flage)=ielem
        
        flagbr=flagbr+1;      dumbnodr(flagbr)=ipnt1
        flagbr=flagbr+1;      dumbnodr(flagbr)=ipnt2
        flagbr=flagbr+1;      dumbnodr(flagbr)=ipnt5
        flagbr=flagbr+1;      dumbnodr(flagbr)=ipnt6

        flagbs=flagbs+1;      dumbnods(flagbs)=ipnt3
        flagbs=flagbs+1;      dumbnods(flagbs)=ipnt4
        flagbs=flagbs+1;      dumbnods(flagbs)=ipnt7
        flagbs=flagbs+1;      dumbnods(flagbs)=ipnt8
        
        label(ielem)=10
        
        elseif(ipminy<dminy+(iprocy+1)*dleny.and.ipmaxy>dminy+(iprocy+1)
     &   *dleny.and.ipminx>dminx+iprocx*dlenx
     &        .and.ipmaxx<dminx+(iprocx+1)*dlenx) then
        
        flage=flage+1;  locelem(flage)=ielem
        
        flagtr=flagtr+1;      dumtnodr(flagtr)=ipnt3
        flagtr=flagtr+1;      dumtnodr(flagtr)=ipnt4
        flagtr=flagtr+1;      dumtnodr(flagtr)=ipnt7
        flagtr=flagtr+1;      dumtnodr(flagtr)=ipnt8

        flagts=flagts+1;      dumtnods(flagts)=ipnt1
        flagts=flagts+1;      dumtnods(flagts)=ipnt2
        flagts=flagts+1;      dumtnods(flagts)=ipnt5
        flagts=flagts+1;      dumtnods(flagts)=ipnt6
        
        label(ielem)=10
        
        elseif(ipminx<dminx+iprocx*dlenx.and.ipmaxx>dminx+iprocx*dlenx
     &  .and.ipminy<dminy+iprocy*dleny.and.ipmaxy>dminy+iprocy*dleny)
     &  then
        
        flage=flage+1;  locelem(flage)=ielem
        
        flaglbr=flaglbr+1;     dumlbnodr(flaglbr)=ipnt1
        flaglbr=flaglbr+1;     dumlbnodr(flaglbr)=ipnt5
        
        flaglbs=flaglbs+1;     dumlbnods(flaglbs)=ipnt3
        flaglbs=flaglbs+1;     dumlbnods(flaglbs)=ipnt7
        
        label(ielem)=10

        elseif(ipminx<dminx+(iprocx+1)*dlenx.and.ipmaxx>dminx+(iprocx+1)
     &  *dlenx.and.ipminy<dminy+iprocy*dleny
     &        .and.ipmaxy>dminy+iprocy*dleny) then
        
        flage=flage+1;  locelem(flage)=ielem
        
        flagrbr=flagrbr+1;      dumrbnodr(flagrbr)=ipnt2
        flagrbr=flagrbr+1;      dumrbnodr(flagrbr)=ipnt6

        flagrbs=flagrbs+1;      dumrbnods(flagrbs)=ipnt4
        flagrbs=flagrbs+1;      dumrbnods(flagrbs)=ipnt8
        
        label(ielem)=10

        elseif(ipminx<dminx+iprocx*dlenx.and.ipmaxx>dminx+iprocx*dlenx
     &  .and.ipminy<dminy+(iprocy+1)*dleny
     &  .and.ipmaxy>dminy+(iprocy+1)*dleny)
     &  then
        
        flage=flage+1;  locelem(flage)=ielem
        
        flagltr=flagltr+1;       dumltnodr(flagltr)=ipnt4
        flagltr=flagltr+1;       dumltnodr(flagltr)=ipnt8

        flaglts=flaglts+1;       dumltnods(flaglts)=ipnt2
        flaglts=flaglts+1;       dumltnods(flaglts)=ipnt6
        
        label(ielem)=10
        
        elseif(ipminx<dminx+(iprocx+1)*dlenx.and.ipmaxx>dminx+(iprocx+1)
     &  *dlenx.and.ipminy<dminy+(iprocy+1)*dleny
     &        .and.ipmaxy>dminy+(iprocy+1)*dleny) then
        
        flage=flage+1;  locelem(flage)=ielem
        
        flagrtr=flagrtr+1;       dumrtnodr(flagrtr)=ipnt3
        flagrtr=flagrtr+1;       dumrtnodr(flagrtr)=ipnt7

        flagrts=flagrts+1;       dumrtnods(flagrts)=ipnt1
        flagrts=flagrts+1;       dumrtnods(flagrts)=ipnt5
        
        label(ielem)=10
        
        endif
        
        enddo  !ielem

!mapping old nodes to new nodes of submesh
        dum=maxval(elem(locelem(1:flage),:))
        allocate(mapold2new(1:dum))
        call mapping3d(elem,nelem,locelem,flage,mapold2new,dum,dumnnod)
        
! dividing the faces into different processors
        flagf=0
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
        
        ymin=min(coord(ipnt1,2),coord(ipnt2,2))
        ymin=min(ymin,coord(ipnt3,2))
        ymin=min(ymin,coord(ipnt4,2))
        
        ymax=max(coord(ipnt1,2),coord(ipnt2,2))
        ymax=max(ymax,coord(ipnt3,2))
        ymax=max(ymax,coord(ipnt4,2))
        
        if(((xmin>dminx+iprocx*dlenx.and.xmin<dminx+(iprocx+1)*dlenx)
     &  .or.(xmax>dminx+iprocx*dlenx.and.xmax<dminx+(iprocx+1)*dlenx))
     &  .and.((ymin>dminy+iprocy*dleny.and.ymin<dminy+(iprocy+1)*dleny)      
     &  .or.(ymax>dminx+iprocy*dleny.and.ymax<dminy+(iprocy+1)*dleny)))                                )
     &  then
        flagf=flagf+1
        dumface(flagf,:)=face(iface,:)
        dumface(flagf,1)=dble(mapold2new(ipnt1))
        dumface(flagf,2)=dble(mapold2new(ipnt2))
        dumface(flagf,3)=dble(mapold2new(ipnt3))
        dumface(flagf,4)=dble(mapold2new(ipnt4))
        endif
        enddo !iface
        dumnface=flagf
        
!Dividing the postprocessing points into different processors
        flagf=0
        dumppnod=0
        do inod=1,pptot
        ipnt1=ppnod(inod)
        xmin=coord(ipnt1,1)
        ymin=coord(ipnt1,2)
        if(xmin>dminx+iprocx*dlenx.and.xmin<dminx+(iprocx+1)*dlenx
     & .and.ymin>dminy+iprocy*dleny.and.ymin<dminy+(iprocy+1)*dleny)
     &    then
        flagf=flagf+1
        dumppnod(flagf)=mapold2new(ipnt1)
        endif
        enddo
        dumpptot=flagf        

!Creating elements for submesh

        allocate(dumelem(1:flage,1:8),dumcoord(1:dumnnod,1:3))
        do ielem=1,flage        

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
        dumrnodr(1)=elem(locelem(flage/2),2)
        endif

        call intrmvdup(dumrnodr(1:flagrr),flagrr,bndcnt)
        
        allocate(rnodr(1:bndcnt))
        rnodrcnt=bndcnt
        do inod=1,rnodrcnt
        rnodr(inod)=mapold2new(dumrnodr(inod))
        enddo

       
! Finding nodes in submesh for right side sending nodes
        if(flagrs==0) then
        flagrs=1
        dumrnods(1)=elem(locelem(flage/2),1)
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
        dumlnodr(1)=elem(locelem(flage/2),1)
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
        dumlnods(1)=elem(locelem(flage/2),2)
        endif

        call intrmvdup(dumlnods(1:flagls),flagls,bndcnt)
        allocate(lnods(1:bndcnt))
        lnodscnt=bndcnt
        do inod=1,lnodscnt
        lnods(inod)=mapold2new(dumlnods(inod))
        enddo

! Finding nodes in submesh for top side receiving nodes
        if(flagtr==0) then
        flagtr=1
        dumtnodr(1)=elem(locelem(flage/2),2)
        endif

        call intrmvdup(dumtnodr(1:flagtr),flagtr,bndcnt)
        
        allocate(tnodr(1:bndcnt))
        tnodrcnt=bndcnt
        do inod=1,tnodrcnt
        tnodr(inod)=mapold2new(dumtnodr(inod))
        enddo
        
! Finding nodes in submesh for top side sending nodes
        if(flagts==0) then
        flagts=1
        dumtnods(1)=elem(locelem(flage/2),1)
        endif

        call intrmvdup(dumtnods(1:flagts),flagts,bndcnt)
        
        allocate(tnods(1:bndcnt))
        tnodscnt=bndcnt
        do inod=1,tnodscnt
        tnods(inod)=mapold2new(dumtnods(inod))
        enddo
        
! Finding nodes in submesh for bottom side receiving nodes
        if(flagbr==0) then
        flagbr=1
        dumbnodr(1)=elem(locelem(flage/2),1)
        endif

        call intrmvdup(dumbnodr(1:flagbr),flagbr,bndcnt)
        allocate(bnodr(1:bndcnt))
        bnodrcnt=bndcnt
        do inod=1,bnodrcnt
        bnodr(inod)=mapold2new(dumbnodr(inod))
        enddo

! Finding nodes in submesh for bottom side sending nodes
        if(flagbs==0) then
        flagbs=1
        dumbnods(1)=elem(locelem(flage/2),2)
        endif

        call intrmvdup(dumbnods(1:flagbs),flagbs,bndcnt)
        allocate(bnods(1:bndcnt))
        bnodscnt=bndcnt
        do inod=1,bnodscnt
        bnods(inod)=mapold2new(dumbnods(inod))
        enddo

! Finding nodes in submesh for Left Bottom sending nodes
        if(flaglbs==0) then
        flaglbs=1
        dumlbnods(1)=elem(locelem(flage/2),2)
        endif

        call intrmvdup(dumlbnods(1:flaglbs),flaglbs,bndcnt)
        allocate(lbnods(1:bndcnt))
        lbnodscnt=bndcnt
        do inod=1,lbnodscnt
        lbnods(inod)=mapold2new(dumlbnods(inod))
        enddo

! Finding nodes in submesh for Left Bottom receiving nodes
        if(flaglbr==0) then
        flaglbr=1
        dumlbnodr(1)=elem(locelem(flage/2),2)
        endif

        call intrmvdup(dumlbnodr(1:flaglbr),flaglbr,bndcnt)
        allocate(lbnodr(1:bndcnt))
        lbnodrcnt=bndcnt
        do inod=1,lbnodrcnt
        lbnodr(inod)=mapold2new(dumlbnodr(inod))
        enddo

! Finding nodes in submesh for Right Bottom sending nodes
        if(flagrbs==0) then
        flagrbs=1
        dumrbnods(1)=elem(locelem(flage/2),2)
        endif

        call intrmvdup(dumrbnods(1:flagrbs),flagrbs,bndcnt)
        allocate(rbnods(1:bndcnt))
        rbnodscnt=bndcnt
        do inod=1,rbnodscnt
        rbnods(inod)=mapold2new(dumrbnods(inod))
        enddo


! Finding nodes in submesh for Right Bottom receiving nodes
        if(flagrbr==0) then
        flagrbr=1
        dumrbnodr(1)=elem(locelem(flage/2),2)
        endif

        call intrmvdup(dumrbnodr(1:flagrbr),flagrbr,bndcnt)
        allocate(rbnodr(1:bndcnt))
        rbnodrcnt=bndcnt
        do inod=1,rbnodrcnt
        rbnodr(inod)=mapold2new(dumrbnodr(inod))
        enddo
        
! Finding nodes in submesh for Left Top sending nodes
        if(flaglts==0) then
        flaglts=1
        dumltnods(1)=elem(locelem(flage/2),2)
        endif

        call intrmvdup(dumltnods(1:flaglts),flaglts,bndcnt)
        allocate(ltnods(1:bndcnt))
        ltnodscnt=bndcnt
        do inod=1,ltnodscnt
        ltnods(inod)=mapold2new(dumltnods(inod))
        enddo
        
! Finding nodes in submesh for Left Top receiving nodes
        if(flagltr==0) then
        flagltr=1
        dumltnodr(1)=elem(locelem(flage/2),2)
        endif

        call intrmvdup(dumltnodr(1:flagltr),flagltr,bndcnt)
        allocate(ltnodr(1:bndcnt))
        ltnodrcnt=bndcnt
        do inod=1,ltnodrcnt
        ltnodr(inod)=mapold2new(dumltnodr(inod))
        enddo

! Finding nodes in submesh for Right Top sending nodes
        if(flagrts==0) then
        flagrts=1
        dumrtnods(1)=elem(locelem(flage/2),2)
        endif

        call intrmvdup(dumrtnods(1:flagrts),flagrts,bndcnt)
        allocate(rtnods(1:bndcnt))
        rtnodscnt=bndcnt
!        write(*,*) myid,rtnodscnt
        do inod=1,rtnodscnt
        rtnods(inod)=mapold2new(dumrtnods(inod))
        enddo        

! Finding nodes in submesh for Right Top receiving nodes
        if(flagrtr==0) then
        flagrtr=1
        dumrtnodr(1)=elem(locelem(flage/2),2)
        endif

        call intrmvdup(dumrtnodr(1:flagrtr),flagrtr,bndcnt)
        allocate(rtnodr(1:bndcnt))
        rtnodrcnt=bndcnt
        do inod=1,rtnodrcnt
        rtnodr(inod)=mapold2new(dumrtnodr(inod))
        enddo 
        
        if(iproc.ne.1) then
! Sending element information to other processors        
        call MPI_Send(flage,1,MPI_INTEGER,proccordtomyid(iprocx,iprocy)
     &  ,1,comm1d,ierr)

        call MPI_Send(dumelem(1:flage,1:8),flage*8,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),2,comm1d,ierr)
     
! Sending coord information to all processors
        call MPI_Send(dumnnod,1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),3,comm1d,ierr)

        call MPI_Send(dumcoord(1:dumnnod,1:3),dumnnod*3,
     & MPI_DOUBLE_PRECISION,proccordtomyid(iprocx,iprocy),4,comm1d,ierr)

! Sending face information to all processors
        call MPI_Send(dumnface,1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),5,comm1d,ierr)

        call MPI_Send(dumface(1:dumnface,1:9),dumnface*9,
     & MPI_DOUBLE_PRECISION,proccordtomyid(iprocx,iprocy),6,comm1d,ierr)

! Sending post processing point information to all processors
        call MPI_Send(dumpptot,1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),7,comm1d,ierr)

        call MPI_Send(dumppnod(1:dumpptot),dumpptot,
     & MPI_INTEGER,proccordtomyid(iprocx,iprocy),8,comm1d,ierr)

! Sending right side exchange nodes information to all processors
        call MPI_Send(rnodrcnt,1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),9,comm1d,ierr)

        call MPI_Send(rnodr(1:rnodrcnt),rnodrcnt,MPI_INTEGER
     &   ,proccordtomyid(iprocx,iprocy),10,comm1d,ierr)

        call MPI_Send(rnodscnt,1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),11,comm1d,ierr)

        call MPI_Send(rnods(1:rnodscnt),rnodscnt,MPI_INTEGER
     &   ,proccordtomyid(iprocx,iprocy),12,comm1d,ierr)

! Sending left side exchange nodes information to all processors
        call MPI_Send(lnodrcnt,1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),13,comm1d,ierr)

        call MPI_Send(lnodr(1:lnodrcnt),lnodrcnt,MPI_INTEGER
     &   ,proccordtomyid(iprocx,iprocy),14,comm1d,ierr)

        call MPI_Send(lnodscnt,1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),15,comm1d,ierr)

        call MPI_Send(lnods(1:lnodscnt),lnodscnt,MPI_INTEGER
     &   ,proccordtomyid(iprocx,iprocy),16,comm1d,ierr)

! Sending top side exchange nodes information to all processors
        call MPI_Send(tnodrcnt,1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),17,comm1d,ierr)

        call MPI_Send(tnodr(1:tnodrcnt),tnodrcnt,MPI_INTEGER
     &   ,proccordtomyid(iprocx,iprocy),18,comm1d,ierr)

        call MPI_Send(tnodscnt,1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),19,comm1d,ierr)

        call MPI_Send(tnods(1:tnodscnt),tnodscnt,MPI_INTEGER
     &   ,proccordtomyid(iprocx,iprocy),20,comm1d,ierr)

! Sending bottom side exchange nodes information to all processors
        call MPI_Send(bnodrcnt,1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),21,comm1d,ierr)

        call MPI_Send(bnodr(1:bnodrcnt),bnodrcnt,MPI_INTEGER
     &   ,proccordtomyid(iprocx,iprocy),22,comm1d,ierr)

        call MPI_Send(bnodscnt,1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),23,comm1d,ierr)

        call MPI_Send(bnods(1:bnodscnt),bnodscnt,MPI_INTEGER
     &   ,proccordtomyid(iprocx,iprocy),24,comm1d,ierr)


! Sending left bottom exchange nodes information to all processors
        call MPI_Send(lbnodrcnt,1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),25,comm1d,ierr)

        call MPI_Send(lbnodr(1:lbnodrcnt),lbnodrcnt,MPI_INTEGER
     &   ,proccordtomyid(iprocx,iprocy),26,comm1d,ierr)

        call MPI_Send(lbnodscnt,1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),27,comm1d,ierr)

        call MPI_Send(lbnods(1:lbnodscnt),lbnodscnt,MPI_INTEGER
     &   ,proccordtomyid(iprocx,iprocy),28,comm1d,ierr)

! Sending right bottom exchange nodes information to all processors
        call MPI_Send(rbnodrcnt,1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),29,comm1d,ierr)

        call MPI_Send(rbnodr(1:rbnodrcnt),rbnodrcnt,MPI_INTEGER
     &   ,proccordtomyid(iprocx,iprocy),30,comm1d,ierr)

        call MPI_Send(rbnodscnt,1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),31,comm1d,ierr)

        call MPI_Send(rbnods(1:rnodscnt),rbnodscnt,MPI_INTEGER
     &   ,proccordtomyid(iprocx,iprocy),32,comm1d,ierr)

! Sending left top exchange nodes information to all processors
        call MPI_Send(ltnodrcnt,1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),33,comm1d,ierr)

        call MPI_Send(ltnodr(1:ltnodrcnt),ltnodrcnt,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),34,comm1d,ierr)

        call MPI_Send(ltnodscnt,1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),35,comm1d,ierr)

        call MPI_Send(ltnods(1:ltnodscnt),ltnodscnt,MPI_INTEGER
     &   ,proccordtomyid(iprocx,iprocy),36,comm1d,ierr)

! Sending right top exchange nodes information to all processors
        call MPI_Send(rtnodrcnt,1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),37,comm1d,ierr)

        call MPI_Send(rtnodr(1:rtnodrcnt),rtnodrcnt,MPI_INTEGER
     &   ,proccordtomyid(iprocx,iprocy),38,comm1d,ierr)

!        write(*,*) "sending",myid,rtnodscnt
        call MPI_SSend(rtnodscnt,1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),45,comm1d,ierr)
!        write(*,*) "sending",myid,rtnodscnt
        call MPI_SSend(rtnods(1:rtnodscnt),rtnodscnt,MPI_INTEGER
     &   ,proccordtomyid(iprocx,iprocy),46,comm1d,ierr)

        deallocate(rnodr,rnods)
        deallocate(lnodr,lnods)
        deallocate(tnodr,tnods)
        deallocate(bnodr,bnods)

        deallocate(lbnodr,lbnods)
        deallocate(rbnodr,rbnods)
        deallocate(ltnodr,ltnods)
        deallocate(rtnodr,rtnods)

        else
! Elements of processor 0
        nelem=flage
        deallocate(elem); 
        allocate(elem(1:nelem,1:8))
        elem=dumelem(1:nelem,1:8)
        
! Nodes of processor 0
        nnod=dumnnod
        deallocate(coord); 
        allocate(coord(1:nnod,1:3))
        coord=dumcoord(1:dumnnod,1:3)

! Faces of processor 0
        nface=dumnface
        deallocate(face); 
        allocate(face(1:nface,1:9))
        face=dumface(1:dumnface,1:9)

! Post processing of processor 0
        pptot=dumpptot
        deallocate(ppnod); 
        allocate(ppnod(1:pptot))
        ppnod=dumppnod(1:dumpptot)
        
        write(*,*) "myid   :",myid
        write(*,*) "Nelem  :",nelem
        write(*,*) "nodes  :",nnod
        write(*,*) "nfaces :",nface
        write(*,*) "pptot  :",pptot
        write(*,*) "lnodrcnt:",lnodrcnt
        write(*,*) "lnodscnt:",lnodscnt
        write(*,*) "rnodrcnt:",rnodrcnt
        write(*,*) "rnodscnt:",rnodscnt
        write(*,*) "bnodrcnt:",bnodrcnt
        write(*,*) "bnodscnt:",bnodscnt
        write(*,*) "tnodrcnt:",tnodrcnt
        write(*,*) "tnodscnt:",tnodscnt

        write(*,*) "lbnodrcnt:",lbnodrcnt
        write(*,*) "lbnodscnt:",lbnodscnt
        write(*,*) "rbnodrcnt:",rbnodrcnt
        write(*,*) "rbnodscnt:",rbnodscnt
        write(*,*) "ltnodrcnt:",ltnodrcnt
        write(*,*) "ltnodscnt:",ltnodscnt
        write(*,*) "rtnodrcnt:",rtnodrcnt
        write(*,*) "rtnodscnt:",rtnodscnt
        write(*,*) "-------------------------------------"
        
        endif       
                
        deallocate(mapold2new)
        deallocate(dumelem,dumcoord)

        enddo
        enddo
        
        
        open(25,file="labelname3d.dat")
        do ielem=1,nelem
        ix=sum(coord(elem(ielem,:),1))/8.0
        iy=sum(coord(elem(ielem,:),2))/8.0
        iz=maxval(coord(elem(ielem,:),3))
        if(iz>1.97e-3) write(25,*) ix,iy,iz,label(ielem)
        enddo
        close(25)
        
        deallocate(label)
        deallocate(locelem,dumface,dumppnod)
        
        deallocate(dumlnodr,dumlnods,dumrnodr,dumrnods,dumbnodr)
        deallocate(dumbnods,dumtnodr,dumtnods,dumlbnodr,dumlbnods)
        deallocate(dumrbnodr,dumrbnods,dumltnodr,dumltnods)
        deallocate(dumrtnodr,dumrtnods)
! Finding the elements of each processor and sending it to the processor        
        else
! Prococessors receiving elements
        call MPI_RECV(nelem,1,MPI_INTEGER,0,1,comm1d,status,ierr)
        allocate(elem(1:nelem,1:8))
        call MPI_RECV(elem,nelem*8,MPI_INTEGER,0,2,comm1d,status,ierr)

! Prococessors receiving nodes
        call MPI_RECV(nnod,1,MPI_INTEGER,0,3,comm1d,status,ierr)
        allocate(coord(1:nnod,1:3))
        call MPI_RECV(coord,nnod*3,MPI_DOUBLE_PRECISION,0,4,
     &  comm1d,status,ierr)
     
! Prococessors receiving nodes
        call MPI_RECV(nface,1,MPI_INTEGER,0,5,comm1d,status,ierr)
        allocate(face(1:nface,1:9))
        call MPI_RECV(face,nface*9,MPI_DOUBLE_PRECISION,0,6,
     &  comm1d,status,ierr)

! Prococessors receiving postprocessor nodes
        call MPI_RECV(pptot,1,MPI_INTEGER,0,7,comm1d,status,ierr)
        allocate(ppnod(1:pptot))
        call MPI_RECV(ppnod,pptot,MPI_INTEGER,0,8,comm1d,status,ierr)

! processors receving right side exchange nodes
        call MPI_RECV(rnodrcnt,1,MPI_INTEGER,0,9,comm1d,status,ierr)
        allocate(rnodr(1:rnodrcnt))
        call MPI_RECV(rnodr,rnodrcnt,MPI_INTEGER,0,10,comm1d
     &  ,status,ierr)

        call MPI_RECV(rnodscnt,1,MPI_INTEGER,0,11,comm1d,status,ierr)
        allocate(rnods(1:rnodscnt))
        call MPI_RECV(rnods,rnodscnt,MPI_INTEGER,0,12,comm1d
     &  ,status,ierr)

! processors receving left side exchange nodes
        call MPI_RECV(lnodrcnt,1,MPI_INTEGER,0,13,comm1d,status,ierr)
        allocate(lnodr(1:lnodrcnt))
        call MPI_RECV(lnodr,lnodrcnt,MPI_INTEGER,0,14,comm1d
     &  ,status,ierr)

        call MPI_RECV(lnodscnt,1,MPI_INTEGER,0,15,comm1d,status,ierr)
        allocate(lnods(1:lnodscnt))
        call MPI_RECV(lnods,lnodscnt,MPI_INTEGER,0,16,comm1d
     &  ,status,ierr)

! processors receving top side exchange nodes
        call MPI_RECV(tnodrcnt,1,MPI_INTEGER,0,17,comm1d,status,ierr)
        allocate(tnodr(1:tnodrcnt))
        call MPI_RECV(tnodr,tnodrcnt,MPI_INTEGER,0,18,comm1d
     &  ,status,ierr)

        call MPI_RECV(tnodscnt,1,MPI_INTEGER,0,19,comm1d,status,ierr)
        allocate(tnods(1:tnodscnt))
        call MPI_RECV(tnods,tnodscnt,MPI_INTEGER,0,20,comm1d
     &  ,status,ierr)

! processors receving bottom side exchange nodes
        call MPI_RECV(bnodrcnt,1,MPI_INTEGER,0,21,comm1d,status,ierr)
        allocate(bnodr(1:bnodrcnt))
        call MPI_RECV(bnodr,bnodrcnt,MPI_INTEGER,0,22,comm1d
     &  ,status,ierr)

        call MPI_RECV(bnodscnt,1,MPI_INTEGER,0,23,comm1d,status,ierr)
        allocate(bnods(1:bnodscnt))
        call MPI_RECV(bnods,bnodscnt,MPI_INTEGER,0,24,comm1d
     &  ,status,ierr)

! processors receving left bottom side exchange nodes
        call MPI_RECV(lbnodrcnt,1,MPI_INTEGER,0,25,comm1d,status,ierr)
        allocate(lbnodr(1:lbnodrcnt))
        call MPI_RECV(lbnodr,lbnodrcnt,MPI_INTEGER,0,26,comm1d
     &  ,status,ierr)

        call MPI_RECV(lbnodscnt,1,MPI_INTEGER,0,27,comm1d,status,ierr)
        allocate(lbnods(1:lbnodscnt))
        call MPI_RECV(lbnods,lbnodscnt,MPI_INTEGER,0,28,comm1d
     &  ,status,ierr)

! processors receving right bottom side exchange nodes
        call MPI_RECV(rbnodrcnt,1,MPI_INTEGER,0,29,comm1d,status,ierr)
        allocate(rbnodr(1:rbnodrcnt))
        call MPI_RECV(rbnodr,rbnodrcnt,MPI_INTEGER,0,30,comm1d
     &  ,status,ierr)

        call MPI_RECV(rbnodscnt,1,MPI_INTEGER,0,31,comm1d,status,ierr)
        allocate(rbnods(1:rbnodscnt))
        call MPI_RECV(rbnods,rbnodscnt,MPI_INTEGER,0,32,comm1d
     &  ,status,ierr)     

! processors receving left top side exchange nodes
        call MPI_RECV(ltnodrcnt,1,MPI_INTEGER,0,33,comm1d,status,ierr)
        allocate(ltnodr(1:ltnodrcnt))
        call MPI_RECV(ltnodr,ltnodrcnt,MPI_INTEGER,0,34,comm1d
     &  ,status,ierr)

        call MPI_RECV(ltnodscnt,1,MPI_INTEGER,0,35,comm1d,status,ierr)
        allocate(ltnods(1:ltnodscnt))
        call MPI_RECV(ltnods,ltnodscnt,MPI_INTEGER,0,36,comm1d
     &  ,status,ierr)

! processors receving right top side exchange nodes
        call MPI_RECV(rtnodrcnt,1,MPI_INTEGER,0,37,comm1d,status,ierr)
        allocate(rtnodr(1:rtnodrcnt))
        call MPI_RECV(rtnodr,rtnodrcnt,MPI_INTEGER,0,38,comm1d
     &  ,status,ierr)

c~         call MPI_RECV(rtnodscnt,1,MPI_INTEGER,0,45,comm1d,status,ierr)
c~ !        write(*,*) "receiving:",myid,rtnodscnt,size(rtnods)
c~          allocate(rtnods(1:rtnodscnt))
c~         call MPI_RECV(rtnods,rtnodscnt,MPI_INTEGER,0,46,comm1d
c~      &  ,status,ierr)


        call MPI_RECV(dum3,1,MPI_INTEGER,0,45,comm1d,status,ierr)
       allocate(dummy(1:dum3))
        call MPI_RECV(dummy,dum3,MPI_INTEGER,0,46,comm1d
     &  ,status,ierr)

        write(*,*) "myid   :",myid
        write(*,*) "Nelem  :",nelem
        write(*,*) "nodes  :",nnod
        write(*,*) "nfaces :",nface
        write(*,*) "pptot  :",pptot
        
        write(*,*) "lnodrcnt:",lnodrcnt
        write(*,*) "lnodscnt:",lnodscnt
        write(*,*) "rnodrcnt:",rnodrcnt
        write(*,*) "rnodscnt:",rnodscnt
        write(*,*) "bnodrcnt:",bnodrcnt
        write(*,*) "bnodscnt:",bnodscnt
        write(*,*) "tnodrcnt:",tnodrcnt
        write(*,*) "tnodscnt:",tnodscnt

        write(*,*) "lbnodrcnt:",lbnodrcnt
        write(*,*) "lbnodscnt:",lbnodscnt
        write(*,*) "rbnodrcnt:",rbnodrcnt
        write(*,*) "rbnodscnt:",rbnodscnt
        write(*,*) "ltnodrcnt:",ltnodrcnt
        write(*,*) "ltnodscnt:",ltnodscnt
        write(*,*) "rtnodrcnt:",rtnodrcnt
        write(*,*) "rtnodscnt:",rtnodscnt
        write(*,*) "-------------------------------------"

           
        endif

!-----------------------------------------------------
! Writing the label of each mesh points to a file
!-----------------------------------------------------
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

c~         if(myid==0) then
c~          Write(*,*) "Elements in each processor:",locnelem
c~          write(*,*) "Total elements:",sum(locnelem)
c~          Write(*,*) "Nodes in each processor:",locnnod
c~          write(*,*) "Total nodes:",sum(locnnod)
c~          Write(*,*) "Faces in each processor:",locnface
c~          write(*,*) "Total faces:",sum(locnface)
c~          Write(*,*) "PPtot in each processor:",locpptot
c~          write(*,*) "Total pptot:",sum(locpptot)
c~          endif

  201  FORMAT(F7.2,F7.2,F15.10,3I3)         
        end

