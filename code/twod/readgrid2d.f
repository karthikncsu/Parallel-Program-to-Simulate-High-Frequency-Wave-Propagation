        subroutine readgrid2d()
        use variables
        implicit none
        logical file_e
        integer dum,n,spcdimc,i
        DOUBLE PRECISION ipnt1,ipnt2,ipnt3,ipnt4,x1,x2,y1,y2
        DOUBLE PRECISION ipnt5,ipnt6,ipnt7,ipnt8
        character(3) test

        inquire(file='mesh.dat', exist=file_e)

        if (.not.file_e) then
            write(*,*) "File doesnot exists" 
            call exit(1)
        end if
         
        open(10,file='mesh.dat')
        
        read(10,*) test,spcdimc
        read(10,*) test,nnod
        
        if(spcdimc.ne.spcdim) then
        write(*,*) "The Grid file provided is not ",spcdim,"dimesion"
        call exit(2)
        end if
        
!------------------------------------------------------------------
! Reading coordinates of the mesh points from mesh file
!------------------------------------------------------------------
        allocate (coord(1:nnod,1:2))
        do i=1,nnod
        read(10,*) coord(i,1),coord(i,2)
        enddo

!------------------------------------------------------------------
! Reading faces points from mesh file
!------------------------------------------------------------------
        read(10,*) test,nface
        allocate (face(1:nface,1:6))
        
        do i=1,nface
        read(10,*) ipnt1,ipnt2,ipnt3,ipnt4,ipnt5
        face(i,1)=ipnt1
        face(i,2)=ipnt2
        face(i,3)=ipnt3
        face(i,4)=ipnt4
        face(i,5)=ipnt5
        x1=coord(int(ipnt1),1); y1=coord(int(ipnt1),2)
        x2=coord(int(ipnt2),1); y2=coord(int(ipnt2),2)
        face(i,6)=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
        enddo
!------------------------------------------------------------------
! Reading elements from mesh file
!------------------------------------------------------------------
        read(10,*) test,nelem
        allocate (elem(1:nelem,1:4))
        
        do i=1,nelem
        read(10,*) ipnt1,ipnt2,ipnt3,ipnt4
        elem(i,1)=int(ipnt1)
        elem(i,2)=int(ipnt2)
        elem(i,3)=int(ipnt4)
        elem(i,4)=int(ipnt3)
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
