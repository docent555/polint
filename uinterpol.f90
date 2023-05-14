module uinterpol
   use, intrinsic :: iso_c_binding

   implicit none

   !real(c_double), private :: za(663), rea(663), ima(663)
   real(c_double), private :: zex
   include 'z.inc'
   include 're.inc'
   include 'im.inc'
   integer(c_int), private :: n

contains

   function squval(zz)

      implicit none

      real(c_double), intent(in) :: zz

      complex(c_double_complex) squval
      real(c_double) z, re, im, dz, z1
      integer(c_int) l

      zex = 185.5
      dz = 0.28011
      z = zz/zex*185.5
      l = z/dz

      if (l .eq. 0) then
         l = 2
      elseif (l .ge. 662) then
         l = 662
      else
         if ((z - za(l)) .gt. 0.5*dz) l = l + 1
      end if

      z1 = za(l)
      z = z - 8.5d0

      re = rea(l - 1) + ((-rea(l - 1) + rea(l))*(dz + z - z1))/dz + ((rea(l - 1)/2.0d0 - rea(l) + &
                                                                      rea(l + 1)/2.0d0)*(z - z1)*(dz + z - z1))/dz/dz
      im = ima(l - 1) + ((-ima(l - 1) + ima(l))*(dz + z - z1))/dz + ((ima(l - 1)/2.0d0 - ima(l) + &
                                                                      ima(l + 1)/2.0d0)*(z - z1)*(dz + z - z1))/dz/dz
      !!!NE RABOTAET
      !re = ((rea(l - 1) - 2*rea(l) + rea(l + 1))*z**2)/(2.*dz**2) &
      !     + (z*(dz*(-rea(l - 1) + rea(l + 1)) - 2*(rea(l - 1) - 2*rea(l) + rea(l + 1))*z1))/(2.*dz**2) + &
      !     -(2*dz**2*rea(l) + dz*(rea(l - 1) - rea(l + 1))*z1 + (rea(l - 1) - 2*rea(l) + rea(l + 1))*z1**2)/(2.*dz**2)
      !im = ((ima(l - 1) - 2*ima(l) + ima(l + 1))*z**2)/(2.*dz**2) + &
      !     (z*(dz*(-ima(l - 1) + ima(l + 1)) - 2*(ima(l - 1) - 2*ima(l) + ima(l + 1))*z1))/(2.*dz**2) + &
      !     -(2*dz**2*ima(l) + dz*(ima(l - 1) - ima(l + 1))*z1 + (ima(l - 1) - 2*ima(l) + ima(l + 1))*z1**2)/(2.*dz**2)

      squval = dcmplx(re, im)

   end function squval

   function linuval(zz)

      implicit none

      real(c_double), intent(in) :: zz

      complex(c_double_complex) linuval
      real(c_double) z, re, im, d
      integer(c_int) l

      !zex = 11.258330249197702
      zex = 185.5
      z = zz/zex*185.5 - 8.5
      l = (z + 8.5)/0.280211 + 1
      d = z - za(l)

      !print *, z, l, d

      if (d .gt. 0.0 .and. l /= 663) then
         re = (rea(l)*za(l + 1) - rea(l + 1)*za(l))/(za(l + 1) - za(l)) + &
              (rea(l + 1) - rea(l))/(za(l + 1) - za(l))*z
         im = (ima(l)*za(l + 1) - ima(l + 1)*za(l))/(za(l + 1) - za(l)) + &
              (ima(l + 1) - ima(l))/(za(l + 1) - za(l))*z
      else if (d .lt. 0.0 .and. l /= 1) then
         re = (rea(l - 1)*za(l) - rea(l)*za(l - 1))/(za(l) - za(l - 1)) + &
              (rea(l) - rea(l - 1))/(za(l) - za(l - 1))*z
         im = (ima(l - 1)*za(l) - ima(l)*za(l - 1))/(za(l) - za(l - 1)) + &
              (ima(l) - ima(l - 1))/(za(l) - za(l - 1))*z
      else
         re = rea(l)
         im = ima(l)
      end if

      linuval = dcmplx(re, im)

   end function linuval

   !subroutine init_uinterpol(n)
   !
   !   implicit none
   !
   !   integer(c_int) :: n, i
   !   real(c_double) d1, d2, d3
   !
   !   open (1, file='table.dat', status='old')
   !   do i = 1, 663
   !      read (1, '(6f12.6)') za(i), d1, d2, d3, rea(i), ima(i)
   !      !write (*, '(6f12.6)') za(i), d1, d2, d3, rea(i), ima(i)
   !   end do
   !   close (1)
   !
   !   open (1, file='z.inc')
   !   write (1, '(a)') 'real(c_double), private :: za(663) =(/ &'
   !   do i = 1, 662
   !      write (1, '(f12.6,a)') za(i), ', &'
   !   end do
   !   write (1, '(f12.6,a)') za(663), '/)'
   !   close (1)
   !
   !   open (1, file='re.inc')
   !   write (1, '(a)') 'real(c_double), private :: rea(663) =(/ &'
   !   do i = 1, 662
   !      write (1, '(f12.6,a)') rea(i), ', &'
   !   end do
   !   write (1, '(f12.6,a)') rea(663), '/)'
   !   close (1)
   !
   !   open (1, file='im.inc')
   !   write (1, '(a)') 'real(c_double), private :: ima(663) =(/ &'
   !   do i = 1, 662
   !      write (1, '(f12.6,a)') ima(i), ', &'
   !   end do
   !   write (1, '(f12.6,a)') ima(663), '/)'
   !   close (1)
   !
   !end subroutine init_uinterpol

end module uinterpol
