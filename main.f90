program p
   use uinterpol

   implicit none

   integer(4), parameter :: n = 663
   real(4) :: dz, z(n), zex
   complex(8) :: u(n)
   integer(4) i, l

   !call init_uinterpol(n)

   !dz = 185.5/(n - 1)
   !do i = 1, 1000
   !   z(i) = -8.5 + (i - 1)*dz
   !   u(i) = uval(z(i))
   !end do

   !zex = 11.258330249197702
   zex = 185.5
   dz = zex/(n - 1)
   do i = 1, n
      z(i) = (i - 1)*dz
      u(i) = squval(dble(z(i)))
      z(i) = z(i) - 8.5
   end do

   open (1, file='interpol.dat')
   do i = 1, n
      write (1, '(3f12.6)') z(i), dreal(u(i)), dimag(u(i))
   end do
   close (1)

end program p

