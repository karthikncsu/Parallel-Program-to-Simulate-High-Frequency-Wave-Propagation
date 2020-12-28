        subroutine readgrid3d()
        use variables
        implicit none
        logical file_e
        integer dum,n,spcdimc,i
        DOUBLE PRECISION ipnt1,ipnt2,ipnt3,ipnt4
        DOUBLE PRECISION ipnt5,ipnt6,ipnt7,ipnt8
        DOUBLE PRECISION x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,area1,area2
        character(3) test

        inquire(file='mesh3d.dat', exist=file_e)

        if (.not.file_e) then
            write(*,*) "File doesnot exists: mesh3d.dat" 
            call exit(1)
        end if
         
        open(10,file='mesh3d.dat')
        
        read(10,*) test,spcdimc
        read(10,*) test,nnod
        
        if(spcdimc.ne.spcdim) then
        write(*,*) "The Grid file provided is not ",spcdim,"dimesion"
        call exit(2)
        end if
        
!------------------------------------------------------------------
! Reading coordinates of the mesh points from mesh file
!------------------------------------------------------------------
        allocate (coord(1:nnod,1:3))
        do i=1,nnod
        read(10,*) coord(i,1),coord(i,2),coord(i,3)
        enddo

!------------------------------------------------------------------
! Reading faces points from mesh file
!------------------------------------------------------------------
        read(10,*) test,nface
        allocate (face(1:nface,1:9))
        
        do i=1,nface
        read(10,*) ipnt1,ipnt2,ipnt3,ipnt4,ipnt5,ipnt6,ipnt7,ipnt8
        face(i,1)=ipnt1
        face(i,2)=ipnt2
        face(i,3)=ipnt4
        face(i,4)=ipnt3
        face(i,5)=ipnt5
        face(i,6)=ipnt6
        face(i,7)=ipnt7
        face(i,8)=ipnt8

        x1=coord(int(ipnt1),1);       x2=coord(int(ipnt2),1)
        y1=coord(int(ipnt1),2);       y2=coord(int(ipnt2),2)
        z1=coord(int(ipnt1),3);       z2=coord(int(ipnt2),3)

        x3=coord(int(ipnt4),1);       x4=coord(int(ipnt3),1)
        y3=coord(int(ipnt4),2);       y4=coord(int(ipnt3),2)
        z3=coord(int(ipnt4),3);       z4=coord(int(ipnt3),3)
        
        call triarea3d(x1,y1,z1,x2,y2,z2,x3,y3,z3,area1)
        call triarea3d(x1,y1,z1,x4,y4,z4,x3,y3,z3,area2)

        face(i,9)=area1+area2
        enddo
!------------------------------------------------------------------
! Reading elements from mesh file
!------------------------------------------------------------------
        read(10,*) test,nelem
        allocate (elem(1:nelem,1:8))
        
        do i=1,nelem
        read(10,*) ipnt1,ipnt2,ipnt3,ipnt4,ipnt5,ipnt6,ipnt7,ipnt8
        elem(i,1)=int(ipnt1)
        elem(i,2)=int(ipnt2)
        elem(i,3)=int(ipnt4)
        elem(i,4)=int(ipnt3)
        elem(i,5)=int(ipnt5)
        elem(i,6)=int(ipnt6)
        elem(i,7)=int(ipnt8)
        elem(i,8)=int(ipnt7)
        enddo
!------------------------------------------------------------------
! Reading postprocessing from mesh file
!------------------------------------------------------------------
        read(10,*) test,pptot
        allocate(ppnod(1:pptot))
        
        do i=1,pptot
        read(10,*) ppnod(i)
        enddo
        
        close(10)

!        write(*,*) "end of readgrid file"
        write(*,*) myid,"Nelem:",nelem,"Nnode:",nnod,"pptot:",pptot
        write(*,*) "Nface:",nface
          
  102  FORMAT(F7.2,F7.2)
  103  FORMAT(I5,I5,I5,I5,I5)          
        
        end
