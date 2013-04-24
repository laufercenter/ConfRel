c      program to calculate free energy from freq
c
       implicit real*8 (a-h,o-z)
c
        dimension nk(100), x1(100), x2(100)
        dimension y1(100), y2(100), y3(100), y4(100)
        dimension z1(100),z2(100),z3(100),z4(100),z5(100)
        dimension z6(100),z7(100),z8(100),z9(100),z10(100)
        dimension g1(100), g2(100), dg(100)
        dimension az(100), bz(100)
c
        ntot=56
c
        open(unit=13,file='confinement_Alpha_GA95',status='old')
c
        do i=1,ntot
           read(13,*)nk(i),y1(i)
        end do
c
        open(unit=14,file='confinement_Beta_GA95',status='old')
c
        do i=1,ntot
           read(14,*)nk(i),y2(i)
        end do
c 

         open(unit=18,file='residue_Alpha_GA95',status='old')
c
        do i=1,ntot
           read(18,10)nk(i),z1(i),z2(i),z3(i),z4(i),z5(i)
 10      format(8x,i2,5f10.3)    
        end do
c
         open(unit=19,file='residue_Beta_GA95',status='old')
c
        do i=1,ntot
           read(19,10)nk(i),z6(i),z7(i),z8(i),z9(i),z10(i)
        end do
c
        open(unit=21,file='enthalpy_spring_Alpha_GA95',status='old')
c
        do i=1,ntot
           read(21,*)nk(i),az(i)
        end do
c
        open(unit=22,file='enthalpy_spring_Beta_GA95',status='old')
c
        do i=1,ntot
           read(22,*)nk(i),bz(i)
        end do
c
         open(unit=16,file='output_GA95',status='unknown')
c
         do i=1,ntot
            g1=z1(i)+z2(i)+z3(i)+z4(i)+(z5(i)*0.0072)+az(i)
            g2=z6(i)+z7(i)+z8(i)+z9(i)+(z10(i)*0.0072)+bz(i)
            dg(i)=y1(i)-y2(i)+g2(i)-g1(i)
            write(16,*)nk(i),dg(i)
         end do
c
         stop
         end
