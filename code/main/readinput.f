        subroutine readinput()
        use variables
        implicit none
        logical file_e
        integer dum,n
        character(5) test

        inquire(file='input_file.yaml', exist=file_e)

        if (.not.file_e) then
            write(*,*) "File doesnot exists" 
            call exit(1)
        end if
         
        open(10,file='input_file.yaml')
        
        do while(1>0)
            read(10,*) test
            if(test=="start") then
            goto 99
            end if
        end do
   99   read(10,*) dum,spcdim
        read(10,*) dum,gridname
        read(10,*) dum,xl
        read(10,*) dum,xr
        read(10,*) dum,dx
        
        nnod=int((xr-xl)/dx)+1

        do while(1>0)
            read(10,*) test
            if(test=="start") then
            goto 100
            end if
        end do
  100   read(10,*) dum,dum
        read(10,*) dum,rho0
        read(10,*) dum,E0
        read(10,*) dum,nu
 
         do while(1>0)
            read(10,*) test
            if(test=="start") then
            goto 101
            end if
        end do
  101   read(10,*) dum,dum
        read(10,*) dum,dum
        read(10,*) dum,dum
        
        do while(1>0)
            read(10,*) test
            if(test=="start") then
            goto 102
            end if
        end do
        
  102   read(10,*) dum,bcl
        read(10,*) dum,bcr
        read(10,*) dum,bclv
        read(10,*) dum,bcrv
        read(10,*) dum,exc_typ
        read(10,*) dum,freq

        do while(1>0)
            read(10,*) test
            if(test=="start") then
            goto 103
            end if
        end do
        
  103   read(10,*) dum,dum
        read(10,*) dum,tmeth
        read(10,*) dum,t0
        read(10,*) dum,tf
        read(10,*) dum,dt
        read(10,*) dum,tsave
        read(10,*) dum,Niter
        read(10,*) dum,tol
        
!        write(*,*) myid,"End of reading input file"
        if(tsave<dt) then
        write(*,*) 'Time interval to save need to be more than dt'
        call exit(2)
        endif
        
        end
