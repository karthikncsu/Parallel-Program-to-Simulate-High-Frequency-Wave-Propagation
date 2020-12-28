        module variables
		implicit none
        include 'mpif.h'
       
! GLobal Constants
        REAL, PARAMETER :: PI= 3.1415926,eps=1e-16
 ! Global Variables        
        integer spcdim,nfls,ndmn,ntyp
        integer stdtyp,stdmet,Niter,disctyp,tmeth
        integer exc_typ
        DOUBLE PRECISION t0,tf,dt,tol,tsave
        DOUBLE PRECISION E0,nu,rho0,freq
        character(20) gridname

! Mesh Variables
        DOUBLE PRECISION, dimension (:,:),allocatable :: coord,face
        integer, dimension (:,:),allocatable :: elem
        integer, dimension (:),allocatable :: psp,pspl
        INTEGER nelem,nnod,nface

!Post processing Variables        
        integer, dimension (:),allocatable :: ppnod
        integer pptot
        
! oned problem variables
        integer bcl,bcr
        DOUBLE PRECISION bclv,bcrv,xl,xr,dx

! Physics Variables
        DOUBLE PRECISION, dimension (:),allocatable :: Kstif,Mass,Ma3d

!Parallel progrmaing variables
        INTEGER nprocs,myid,ierr
        INTEGER, dimension (:), allocatable :: label,locnelem,locnnod
        INTEGER, dimension (:), allocatable :: locnface,locpptot
        INTEGER,dimension (:), allocatable :: lnodr,lnods,rnodr,rnods
        INTEGER,dimension (:), allocatable :: bnodr,bnods,tnodr,tnods
        INTEGER,dimension (:),allocatable :: lbnodr,lbnods,rbnodr,rbnods
        INTEGER,dimension (:),allocatable :: ltnodr,ltnods,rtnodr,rtnods
!        INTEGER,dimension (:),allocatable :: dumrtnods
        INTEGER,dimension (:,:),allocatable :: proccordtomyid
        INTEGER lnodrcnt,lnodscnt,rnodrcnt,rnodscnt
        INTEGER bnodrcnt,bnodscnt,tnodrcnt,tnodscnt
        INTEGER lbnodrcnt,lbnodscnt,rbnodrcnt,rbnodscnt
        INTEGER ltnodrcnt,ltnodscnt,rtnodrcnt,rtnodscnt
        INTEGER, dimension (:), allocatable :: dims,isoperiodic
        INTEGER ndims,comm1d,reorder,status
        INTEGER left,right,top,bottom
        INTEGER lefttop,righttop,leftbottom,rightbottom
        INTEGER nprocx,nprocy
        INTEGER proccoord(2)     
        
        end module
	
