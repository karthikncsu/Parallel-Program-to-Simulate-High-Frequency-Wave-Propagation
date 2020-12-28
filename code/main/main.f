        program main
        use variables
        implicit none
        
        INTEGER iproc,diff


        call MPI_init(ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierr )
        call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierr )
    
        call readinput()
        call checkinput()
        
! Two dimensional problem
       if(spcdim==2) then
        write(*,*) ndmn
         ndims=1
        allocate(dims(1:ndims),isoperiodic(1:ndims))
        
        isoperiodic(1)=0
        reorder=1
        dims(1)=nprocs
        else

! Three dimensional problem

! Domain Decomposition for 3D probelm        
        nprocx=nprocs; nprocy=1;
        diff=abs(nprocx-nprocy)
        do iproc=1,nprocs/2
        if(mod(nprocs,iproc)==0) then
                if(abs(nprocs/iproc-iproc)<diff) then
                    nprocx=iproc; nprocy=nprocs/iproc
                    diff=abs(nprocs/iproc-iproc)
                endif
        endif
        enddo
        
        if(nprocx<nprocy) then
        diff=nprocx; nprocx=nprocy; nprocy=diff
        endif
        
        if(myid==0) then
        write(*,*) "Number of Processors in x-direction:",nprocx
        write(*,*) "Number of Processors in y-direction:",nprocy
        endif
        endif
        
        call phymanager()
       
        call MPI_Finalize ( ierr )


        end

