c   program to calculate the diff in DG values
c
        dimension nk(100), x1(100), x2(100)
        dimension dx(100)
c 
        ntot=56
c
        open(unit=13,file='output_GB95',status='old')
c
        do i=1,ntot
           read(13,*)nk(i),x1(i)
        end do
c
        open(unit=14,file='output_GA95',status='old')
c
        do i=1,ntot
           read(14,*)nk(i),x2(i)
        end do
c
        open(unit=16,file='output_diff',status='unknown')
c
         do i=1,ntot
            dx(i)=x1(i)-x2(i)
            write(16,*)nk(i),dx(i)
c
         end do
c
         stop
         end
