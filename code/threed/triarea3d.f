        subroutine triarea3d(x1,y1,z1,x2,y2,z2,x3,y3,z3,area)
        implicit none
        DOUBLE PRECISION x1,y1,z1,x2,y2,z2,x3,y3,z3,area
        Double PRECISION a,b,c,s
        
        a=sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
        b=sqrt((x2-x3)**2+(y2-y3)**2+(z2-z3)**2)
        c=sqrt((x1-x3)**2+(y1-y3)**2+(z1-z3)**2)

        s=0.5*(a+b+c)
        area=sqrt(s*(s-a)*(s-b)*(s-c))
  
        
        end
