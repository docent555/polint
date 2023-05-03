module uinterpol
   use, intrinsic :: iso_c_binding

   implicit none

   !real(c_double), private :: za(663), rea(663), ima(663)
   include 'z.inc'
   include 're.inc'
   include 'im.inc'
   integer(c_int), private :: n

contains

   function ucalc(z)

      implicit none

      complex(c_double_complex) ucalc
      real(c_double) z, re, im, d
      integer(c_int) l

      l = (z + 8.5)/0.28021 + 1

      d = z - za(l)

      !print *, z, l, d

      if (d .gt. 0.0) then
         re = (rea(l)*za(l + 1) - rea(l + 1)*za(l))/(za(l + 1) - za(l)) + &
              (rea(l + 1) - rea(l))/(za(l + 1) - za(l))*z
         im = (ima(l)*za(l + 1) - ima(l + 1)*za(l))/(za(l + 1) - za(l)) + &
              (ima(l + 1) - ima(l))/(za(l + 1) - za(l))*z
      else if (d .lt. 0.0) then
         re = (rea(l - 1)*za(l) - rea(l)*za(l - 1))/(za(l) - za(l - 1)) + &
              (rea(l) - rea(l - 1))/(za(l) - za(l - 1))*z
         im = (ima(l - 1)*za(l) - ima(l)*za(l - 1))/(za(l) - za(l - 1)) + &
              (ima(l) - ima(l - 1))/(za(l) - za(l - 1))*z
      else
         re = rea(l)
         im = ima(l)
      end if

      ucalc = dcmplx(re, im)

   end function ucalc

   subroutine init_uinterpol(n)

      implicit none

      integer(c_int) :: n, i
      real(c_double) d1, d2, d3

      open (1, file='table.dat', status='old')
      do i = 1, 663
         read (1, '(6f12.6)') za(i), d1, d2, d3, rea(i), ima(i)
         !write (*, '(6f12.6)') za(i), d1, d2, d3, rea(i), ima(i)
      end do
      close (1)

      open (1, file='z.inc')
      write (1, '(a)') 'real(c_double), private :: za(663) =(/ &'
      do i = 1, 662
         write (1, '(f12.6,a)') za(i), ', &'
      end do
      write (1, '(f12.6,a)') za(663), '/)'
      close (1)

      open (1, file='re.inc')
      write (1, '(a)') 'real(c_double), private :: rea(663) =(/ &'
      do i = 1, 662
         write (1, '(f12.6,a)') rea(i), ', &'
      end do
      write (1, '(f12.6,a)') rea(663), '/)'
      close (1)

      open (1, file='im.inc')
      write (1, '(a)') 'real(c_double), private :: ima(663) =(/ &'
      do i = 1, 662
         write (1, '(f12.6,a)') ima(i), ', &'
      end do
      write (1, '(f12.6,a)') ima(663), '/)'
      close (1)

   end subroutine init_uinterpol

end module uinterpol
