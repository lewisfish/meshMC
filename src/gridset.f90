MODULE gridset_mod

implicit none
save

CONTAINS
   subroutine gridset(id)

   use constants, only : nxg, nyg, nzg, xmax, ymax, zmax
   use iarray,    only : xface, yface, zface!, rhokap
   ! use opt_prop,  only : kappa

   implicit none

   integer, intent(IN) :: id
   integer :: i, j, k
   real    :: x, y, z

   if(id == 0)then
      print*, ' '
      print *, 'Setting up density grid....'
   end if
   !**********  Linear Cartesian grid. Set up grid faces ****************
   do i = 1, nxg+1
      xface(i)=(i-1)*2.*xmax/nxg
   end do
   do i = 1, nyg+1
      yface(i)=(i-1)*2.*ymax /nyg
   end do
   do i = 1, nzg+1
      zface(i)=(i-1)*2.*zmax/nzg
   end do
   ! call init_opt4
   !**************  Loop through x, y, and z to set up grid density and refractive index grid.  ****
   do i = 1, nxg
      do j = 1, nyg
         do k = 1, nzg
            x = xface(i)-xmax+xmax/nxg
            y = yface(j)-ymax+ymax/nyg
            z = zface(k)-zmax+zmax/nzg
!***********Call density setup subroutine
            
               ! if(i >= 355 .and. i<= 395)then
               !    if(j >= 355 .and. j<= 395)then
               !       if(k >= 355 .and. k<= 395)then
               !          rhokap(i,j,k)=100000.*kappa
               !       else
               !          rhokap(i,j,k)=kappa
               !       end if
               !    else
               !       rhokap(i,j,k)=kappa
               !    end if
               ! else
               !    rhokap(i,j,k)=kappa
               ! end if
         end do
      end do
   end do

   
   end subroutine gridset
end MODULE gridset_mod
