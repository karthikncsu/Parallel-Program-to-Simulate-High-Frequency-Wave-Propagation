        subroutine planestrain_fd_expl()
        use variables
        implicit NONE
        real *16, dimension (:,:),allocatable :: unew,uold,uold0
        real *16 alpha,c0,time,bclv_ti,bcrv_ti,omega,Ac,Bc,Cc,excit
!        real *16 k1(nnod),k2(nnod),k3(nnod),k4(nnod),k5(nnod)
!        real *16 k6(nnod),k7(nnod),k8(nnod),k9(nnod),k10(nnod)
!        real k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,excit
        real K1x,K2x,K3x,K4x,K5x
        real K1y,K2y,K3y,K4y,K5y
        real ku1,ku2,kv1,kv2
        integer inod,nsol,pnt1,pnt2,pnt3,pnt4,pnt5,pnt6,pnt7,pnt8
        
        allocate (unew(0:nnod,1:2),uold(0:nnod,1:2),uold0(0:nnod,1:2))

        Ac=E0/(rho0*(1+nu)*(1-2.0*nu))*(1-nu)
        Bc=E0/(rho0*(1+nu)*(1-2.0*nu))*(0.5-nu)
        Cc=0.5*E0/(rho0*(1+nu)*(1-2.0*nu))
        nsol=0
        
        time=t0
        omega=2.0*PI*freq
        do inod=0,nnod
        uold(inod,1)=0
        uold(inod,2)=0
        uold0(inod,1)=0
        uold0(inod,2)=0
        unew(inod,1)=0
        unew(inod,2)=0
        end do
        
        open(25,file='animation_wave.dat')
!        write(25,*) 'zone I=',nnod
!        write(25,*) 'STRANDID=', nsol, 'SOLUTIONTIME=',0.0
        do inod = 1,nnod
        K1x=coord(inod,1)
        K1y=coord(inod,2)
        write (25,*) K1x,K1y,time,unew(inod,1),unew(inod,2)
        end do
        do while(time<=tf)
               
!            write(*,*) 'Time:',time
            time=time+dt
            nsol=nsol+1
            if(time<10.0*PI/omega) then
              excit=bclv*sin(omega*time)*(sin(omega*time/10.0))**2
            else
              excit=0
            end if
            
             
        do inod=1,nnod

        pnt1=posupo(inod,1)
        pnt2=posupo(inod,2)
        pnt3=posupo(inod,3)
        pnt4=posupo(inod,4)
        pnt5=posupo(inod,5)
        pnt6=posupo(inod,6)
        pnt7=posupo(inod,7)
        pnt8=posupo(inod,8)

        K1x=coord(pnt1,1)-coord(inod,1)
        K2x=coord(inod,1)-coord(pnt5,1)
        K3x=K1x-K2x
        K1y=coord(pnt3,2)-coord(inod,2)
        K2y=coord(inod,2)-coord(pnt7,2)
        K3y=K1y-K2y
            
        K4y=coord(pnt2,2)+coord(pnt3,2)-coord(pnt8,2)-coord(pnt7,2)
        Ku1=uold(pnt2,1)+uold(pnt3,1)-uold(pnt8,1)-uold(pnt7,1)
        Kv1=uold(pnt2,2)+uold(pnt3,2)-uold(pnt8,2)-uold(pnt7,2)
        K5y=coord(pnt3,2)+coord(pnt4,2)-coord(pnt7,2)-coord(pnt6,2)
        Ku2=uold(pnt3,1)+uold(pnt4,1)-uold(pnt7,1)-uold(pnt6,1)
        Kv2=uold(pnt3,2)+uold(pnt4,2)-uold(pnt7,2)-uold(pnt6,2)
        
! Calculting u-displacement values            
        unew(inod,1)=Ac*2.d0*uold(pnt1,1)/(K3x*K1x)
        unew(inod,1)=unew(inod,1)+Ac*2.d0*uold(pnt5,1)/(K3x*K2x)
        unew(inod,1)=unew(inod,1)-Ac*2.d0*uold(inod,1)/(K1x*K2x)
        
        unew(inod,1)=unew(inod,1)+Bc*2.d0*uold(pnt3,1)/(K3y*K1y)
        unew(inod,1)=unew(inod,1)+Bc*2.d0*uold(pnt7,1)/(K3y*K2y)
        unew(inod,1)=unew(inod,1)-Bc*2.d0*uold(inod,1)/(K1y*K2y)
                
        unew(inod,1)=unew(inod,1)+Cc*2.d0/K3x*(Kv1/K4y-Kv2/K5y)
     
        unew(inod,1)=dt*dt*unew(inod,1)+2*uold(inod,1)-uold0(inod,1)
        
! Calculting v-displacement values            
        unew(inod,2)=Bc*2.d0*uold(pnt1,2)/(K3x*K1x)
        unew(inod,2)=unew(inod,2)+Bc*2.d0*uold(pnt5,2)/(K3x*K2x)
        unew(inod,2)=unew(inod,2)-Bc*2.d0*uold(inod,2)/(K1x*K2x)
        
        unew(inod,2)=unew(inod,2)+Ac*2.d0*uold(pnt3,2)/(K3y*K1y)
        unew(inod,2)=unew(inod,2)+Ac*2.d0*uold(pnt7,2)/(K3y*K2y)
        unew(inod,2)=unew(inod,2)-Ac*2.d0*uold(inod,2)/(K1y*K2y)
                
        unew(inod,2)=unew(inod,2)+Cc*2.d0/K3x*(Ku1/K4y-Ku2/K5y)
     
        unew(inod,2)=dt*dt*unew(inod,2)+2*uold(inod,2)-uold0(inod,2)



        end do
             
        uold0=uold
        uold=unew

        do inod = 1,nnod
        K1x=coord(inod,1)
        K1y=coord(inod,2)
        write (25,*) K1x,K1y,time,unew(inod,1),unew(inod,2)
        end do
           
        end do
        write (*,*) 'CFL: ',sqrt(alpha)
        write (*,*) 'dt: ',dt
        write(*,*) 'omega: ',omega
        write(*,*) 'nnod: ',nnod

        CLOSE(25)
  101   format(E13.6,X,E13.6,X,E13.6)
        
        
        end
