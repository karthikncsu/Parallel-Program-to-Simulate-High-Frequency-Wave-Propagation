        subroutine divofdof()        
! Subroutine to division the DOF into different processors
        use variables
        IMPLICIT NONE

        INTEGER,PARAMETER :: maxall=1e4,maxb=1e3
        INTEGER iproc,dum,elemloc,ielem,dumelem(maxall),pnelm
        INTEGER dumrelm
        
c        if(myid==0) then
c        allocate(pnelm(1:nprocs))
c!---------------------------------------------------------------
c! Finding the mesh elements of each processor from label
c!---------------------------------------------------------------
c        do iproc=nprocs,1,-1
        
c        dumelem(:)=0
c        pnelm=0
c        do ielem=1,nelem
c            if(label(ielem)==2*(iproc-1)) then
c            pnelm=pnelm+1
c            dumelem(pnelm)=ielem
c            else if(label(ielem)==2*(iproc-1)+1) then
            
            
c        enddo
        
        
        
c        enddo
        
        
        
        
        
        
        
        
        
c        else

c        endif

c        write(*,*) "Proc",myid,"reached divofdof"
 

        end

