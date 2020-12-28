        subroutine waveeqn1d_fd_expl()
        use variables
        implicit NONE
        real *16 unew(nnod),uold(nnod),uold0(nnod)
        real *16 alpha,c0,time,bclv_ti,bcrv_ti,omega
        integer inod,nsol
        
        nsol=0
        c0=E0/rho0
        alpha=c0*dt**2/dx**2
        time=t0
        omega=2.0*PI*freq
        do inod=1,nnod
        uold(inod)=0
        uold0(inod)=0
        unew(inod)=0
        end do
        
        open(25,file='animation_wave.dat')
!        write(25,*) 'zone I=',nnod
!        write(25,*) 'STRANDID=', nsol, 'SOLUTIONTIME=',0.0
        do inod = 1,nnod
        write (25,*) xl+dble(inod-1)*dx,time,unew(inod)
        end do
        write(*,*) bcrv
        write(*,*) bcr
        do while(time<=tf)
               
            write(*,*) 'Time:',time
            time=time+dt
            nsol=nsol+1
            if(time<10.0*PI/omega) then
              bclv_ti=bclv*sin(omega*time)*(sin(omega*time/10.0))**2
            else
              bclv_ti=0
            end if

             if(bcl==1) then
             unew(1)=bclv_ti
             else
             unew(1)=alpha*(2*uold(2)+2*dx/E0*bclv_ti-2*uold(1))
             unew(1)=unew(1)+2*uold(1)-uold0(1)
             end if

             if(bcr==1) then
             unew(nnod)=bcrv
             else
             unew(nnod)=(2*uold(nnod-1)-2*dx/E0*bcrv-2*uold(nnod))
             unew(nnod)=alpha*unew(nnod)
             unew(nnod)=unew(nnod)+2*uold(nnod)-uold0(nnod)
             end if
             
             do inod=2,nnod-1
             unew(inod)=(uold(inod-1)-2*uold(inod)+uold(inod+1))
             unew(inod)=alpha*unew(inod)
             unew(inod)=unew(inod)+2*uold(inod)-uold0(inod)
             end do
             
             uold0=uold
             uold=unew
!            write (25,*) time 
            do inod = 1,nnod
            write (25,*) xl+dble(inod-1)*dx,time,unew(inod)
            end do
       
        end do
        write (*,*) 'CFL: ',sqrt(alpha)
        write (*,*) 'dt: ',dt
        write(*,*) 'omega: ',omega
        write(*,*) 'nnod: ',nnod

        CLOSE(25)
  101   format(E13.6,X,E13.6,X,E13.6)
        
        
        end
