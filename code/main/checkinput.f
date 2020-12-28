        subroutine checkinput()
        use variables
        implicit none

        open(10,file='input_check.yaml')
        
        write(10,100) 'space dim :',spcdim
        write(10,*) 'Grid name :',gridname

        write(10,101) 'Density :',rho0
        write(10,101) 'Youngs Modulus :',E0
        write(10,101) 'Poissons Ratio :',nu
        write(10,100) 'bcl :',bcl
        write(10,100) 'bcr :',bcr
        write(10,101) 'bclv :',bclv
        write(10,101) 'bcrv :',bcrv
        write(10,100) 'Study Method :',tmeth
        write(10,101) 't0 :',t0
        write(10,101) 'tf :',tf
        write(10,101) 'dt :',dt
        write(10,100) 'Niter :',Niter
        write(10,101) 'Tolerance :',tol
        
  100   format(A,I5)
  101   format(A,E8.1)
        end
