program p
   use uinterpol

   implicit none

   integer(4), parameter :: n = 1000
   real(8) :: dz, z(n)
   complex(8) :: u(n)
   integer(4) i, l

   call init_uinterpol(n)

   dz = 185.5/(n - 1)

   do i = 1, 1000
      z(i) = -8.5 + (i - 1)*dz
      u(i) = ucalc(z(i))
   end do

   open (1, file='interpol.dat')
   do i = 1, 1000
      write (1, '(3f12.6)') z(i), dreal(u(i)), dimag(u(i))
   end do
   close (1)

end program p

