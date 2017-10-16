MODULE sourceph_mod

    implicit none
    save

    private
    public :: sourceph

    CONTAINS
        subroutine sourceph(xcell, ycell, zcell, iseed)

            use constants,   only : nxg, nyg, nzg, xmax, ymax, zmax
            use photon_vars, only : xp, yp, zp, nxp, nyp, nzp, sint, cost, sinp, cosp

            implicit none

            integer, intent(OUT)   :: xcell, ycell, zcell
            integer, intent(INOUT) :: iseed

            call point(iseed)
            nxp = sint * cosp  
            nyp = sint * sinp
            nzp = cost

            xcell = int(nxg * (xp + xmax) / (2. * xmax)) + 1
            ycell = int(nyg * (yp + ymax) / (2. * ymax)) + 1
            zcell = int(nzg * (zp + zmax) / (2. * zmax)) + 1

        end subroutine sourceph


        subroutine point(iseed)

            use constants,   only : twopi, in
            use photon_vars, only : xp, yp, zp, phi, cosp, sinp, cost, sint

            implicit none

            integer, intent(INOUT) :: iseed

            real :: ran2

            zp = 0.!1.
            xp = 0.
            yp = 0.

            cost = 2. * ran2(iseed) - 1.
            sint = sqrt(1. - cost * cost) 

            phi = twopi * ran2(iseed)
            cosp = cos(phi)
            sinp = sin(phi)
            in = .false.

        end subroutine point


        ! subroutine uniform(iseed)

        !     use photon_vars, only : phi, phase, cost, sint, cosp, sinp, xp, yp, zp
        !     use constants,   only : xmax, ymax, zmax, nzg

        !     implicit none

        !     integer, intent(INOUT) :: iseed
        !     real :: ran2

        !     zp = zmax - (1.e-5*(2.*zmax/nzg))
        !     xp = 2. * xmax * ran2(iseed) - xmax
        !     yp = 2. * ymax * ran2(iseed) - xmax

        !     cost = -1.
        !     sint = 1. - cost**2.

        !     phi = 0.
        !     cosp = cos(phi)
        !     sinp = sin(phi)

        !     phase = 0.
        ! end subroutine uniform


        ! real function rang(avg, sigma, iseed)

        !     implicit none

        !     real,    intent(IN) :: avg, sigma
        !     integer, intent(INOUT) :: iseed
            
        !     real :: s, u, tmp

        !     s = 1.

        !     do while(s >= 1.)
        !         u = ranu(-1., 1., iseed)
        !         s = ranu(-1., 1., iseed)
        !         s = s**2 + u**2
        !     end do

        !     tmp = u*sqrt(-2.*log(s)/s)
        !     rang = avg + sigma*tmp

        ! end function rang


        ! real function ranu(a, b, iseed)

        !     implicit none


        !     real, intent(IN)       :: a, b
        !     integer, intent(INOUT) :: iseed

        !     real :: ran2

        !     ranu = a + ran2(iseed) * (b - a)

        ! end function ranu

end MODULE sourceph_mod
