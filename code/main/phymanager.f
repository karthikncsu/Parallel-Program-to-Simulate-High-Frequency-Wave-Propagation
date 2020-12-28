        subroutine phymanager()
        use variables
        implicit none
        INTEGER dum
        
    
        if (spcdim==1) then
!----------------------------------------------------------
! 1D Space dimension
!----------------------------------------------------------
            call waveeqn1d_fd_expl()
!----------------------------------------------------------
! 2D Space dimension
!----------------------------------------------------------
        else if(spcdim==2) then
            call readgrid2d()
            
            if(nprocs==1) then
                call planestrain_fem_expl()
            else
        
!        write(*,*) '---------------------------------------'

      	call MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,isoperiodic,
     &        reorder,comm1d,ierr)
        call MPI_Comm_rank(comm1d,myid,ierr)
        call MPI_cart_shift(comm1d,0,1,left,right,ierr)
               
            call divisionofgrid2d()

            if(myid==0) then
! BCAST is to make remaining processor wait untill the 1st processor
! reads the grid and labels the grid with the procesor number
                call MPI_BCAST(1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            else
                call MPI_BCAST(dum,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            endif !myid==0

            call planestrain_fem_expl_parallel()
        endif !nprocs==1
!----------------------------------------------------------
! 3D Space dimension
!----------------------------------------------------------
       else if(spcdim==3) then
            if(nprocs==1) then
             call readgrid3d()
             call FEM_expl()
            elseif(nprocs.gt.1) then

! oned domain decomposition            
                if(nprocy==1) then
                ndims=1
                allocate(dims(1:ndims),isoperiodic(1:ndims))
                isoperiodic(1)=0;
                reorder=1;
                dims(1)=nprocx;
      	call MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,isoperiodic,
     &        reorder,comm1d,ierr)
        call MPI_Comm_rank(comm1d,myid,ierr)
        call MPI_cart_shift(comm1d,0,1,left,right,ierr)

!                call divisionofgrid3d_1dx()
                call divisionofgrid3d_1dx_parallel()

        if(myid==0) then
! BCAST is to make all the processor wait ntill 1st processor is done
             call MPI_BCAST(1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        else
             call MPI_BCAST(dum,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        endif !myid==0 

                call FEM_expl_parallel_1dx()
                
                else
! twod domain decomposition
                ndims=2
                allocate(dims(1:ndims),isoperiodic(1:ndims))
                isoperiodic(1)=0; isoperiodic(2)=0;
                reorder=1;
                dims(1)=nprocx; dims(2)=nprocy;
      	call MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,isoperiodic,
     &        reorder,comm1d,ierr)
        call MPI_Comm_rank(comm1d,myid,ierr)
        call MPI_cart_shift(comm1d,0,1,left,right,ierr)

                call divisionofgrid3d_2d_parallel()

        if(myid==0) then
! BCAST is to make remaining processor wait untill the 1st processor
! reads the grid and labels the grid with the procesor number
             call MPI_BCAST(1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        else
             call MPI_BCAST(dum,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        endif !myid==0 

                call FEM_expl_parallel_2d()

               
                endif !nprocy==1

        endif !nprocs.gt.1
        endif !spcdim==3
                        
            
            
        end
