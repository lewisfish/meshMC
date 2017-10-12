module inttau2

   implicit none
   
   private
   public :: tauint1

CONTAINS

    subroutine tauint1(xcell,ycell,zcell,tflag,iseed,delta)
    !optical depth integration subroutine
    !
    !
        use constants,   only : xmax, ymax, zmax, in
        use photon_vars, only : xp, yp, zp, nxp, nyp, nzp, sint, cost, sinp, cosp, phi
        use iarray,      only : jmeanin, jmeanout, tarray
        use opt_prop,    only : meshkappa, mediumkappa, meshalbedo, mediumalbedo, albedo
        use types,       only : vector, operator(*), operator(+), normal
        use subs,        only : barycentric
        use kdtree_mod, only : kdtree

        implicit none

        real,    intent(IN)    :: delta
        integer, intent(INOUT) :: xcell, ycell, zcell, iseed
        logical, intent(INOUT) :: tflag

        real                   :: tau, taurun, taucell, xcur, ycur, zcur, d, dcell, ran2, dmesh, kappa, n1, n2
        integer                :: celli, cellj, cellk, tri
        logical                :: dir(3), rflag
        type(vector) :: v1, v2, v3, p, bary, intrpNormal, incd

        xcur = xp + xmax
        ycur = yp + ymax
        zcur = zp + zmax

        celli = xcell
        cellj = ycell
        cellk = zcell

        taurun = 0.
        d = 0.
        dir = (/.FALSE., .FALSE., .FALSE./)

        tau = -log(ran2(iseed))
        dmesh = huge(1.)
        call get_mesh_dist_tree(kdtree, xcur, ycur, zcur, dmesh, tri)

        !make sure photon is actually in mesh
        if((in) .and. (dmesh >= huge(1.)) )in = .not.in

        if(in)then
            kappa = meshkappa
            albedo = meshalbedo
        else
            kappa = mediumkappa
            albedo = mediumalbedo
        end if

        do
            dir = (/.FALSE., .FALSE., .FALSE./)
            dcell = wall_dist(celli, cellj, cellk, xcur, ycur, zcur, dir)
            taucell = dcell * kappa
            if(taurun + taucell < tau)then

                if(d + dcell >= dmesh)then
                    dcell = dmesh - d
                    dcell = dcell + delta
                    d = d + dcell
                    taucell = dcell*kappa 
                    taurun = taurun + taucell

                    call update_pos(xcur, ycur, zcur, celli, cellj, cellk, dcell, .TRUE., dir, delta)

                    if(in)then
                        in = .not. in
                        kappa = mediumkappa
                        albedo = mediumalbedo
                    else
                        in = .not. in
                        kappa = meshkappa
                        albedo = meshalbedo
                    end if
                    exit
                else
                    taurun = taurun + taucell
                    d = d + dcell
                    call update_pos(xcur, ycur, zcur, celli, cellj, cellk, dcell, .TRUE., dir, delta)
                end if
            if(celli == -1 .or. cellj == -1 .or. cellk == -1)then
                tflag = .true.
                exit
            end if
                if(in)then
                    jmeanin(celli, cellj, cellk) = jmeanin(celli, cellj, cellk) + dcell
                else
                    jmeanout(celli, cellj, cellk) = jmeanout(celli, cellj, cellk) + dcell
                end if
            else

                dcell = (tau - taurun) / kappa
                if(d + dcell >= dmesh)then
                    dcell = dmesh - d
                    d = d + dcell
                    taucell = dcell*kappa 
                    taurun = taurun + taucell
                    call update_pos(xcur, ycur, zcur, celli, cellj, cellk, dcell, .TRUE., dir, delta)
               
                    if(in)then
                        in = .not. in
                        kappa = mediumkappa
                        albedo = mediumalbedo
                    else
                        in = .not. in
                        kappa = meshkappa
                        albedo = meshalbedo
                    end if

                    exit
                else
                    d = d + dcell
                    call update_pos(xcur, ycur, zcur, celli, cellj, cellk, dcell, .FALSE., dir, delta)
                    exit
                end if
                if(in)then
                    jmeanin(celli, cellj, cellk) = jmeanin(celli, cellj, cellk) + dcell
                else
                    jmeanout(celli, cellj, cellk) = jmeanout(celli, cellj, cellk) + dcell
                end if
            end if
            if(celli == -1 .or. cellj == -1 .or. cellk == -1)then
                tflag = .true.
                exit
            end if
        end do
            if(celli == -1 .or. cellj == -1 .or. cellk == -1)then
                tflag = .true.
            end if
        xp = xcur - xmax
        yp = ycur - ymax
        zp = zcur - zmax
        xcell = celli
        ycell = cellj
        zcell = cellk
    end subroutine tauint1
   

    recursive subroutine get_mesh_dist_tree(element, x, y, z, maxval, tri)

        use kdtree_mod,  only : intersect_bbox, node
        use photon_vars, only : nxp, nyp, nzp
        use constants,   only : xmax, ymax, zmax, in
        use types,       only : vector

        implicit none

        type(node), intent(IN)    :: element
        real,       intent(IN)    :: x, y, z
        real :: distance
        
        integer      :: tri
        type(vector) :: ray, origin
        logical      :: flag
        real :: maxval

        if(element%level == 0)maxval=huge(1.)

        ray    = vector(nxp, nyp, nzp)
        origin = vector(x-xmax, y-ymax, z-zmax)

        flag = intersect_bbox(origin, ray, element%mindim, element%maxdim)
        if(flag)then
            if(associated(element%rn) .or. associated(element%ln))then
                call get_mesh_dist_tree(element%ln, x, y, z, maxval, tri)
                call get_mesh_dist_tree(element%rn, x, y, z, maxval, tri)
            else
                !hit leaf
                call ray_triangle_intersection(element%list, x, y, z, distance, tri)
                if(distance < maxval)then
                    maxval = distance
                    return
                end if
            end if
        else
            distance = huge(1.)
            return
        end if
    end subroutine get_mesh_dist_tree


    subroutine ray_triangle_intersection(list, x, y, z, distance, tri)

        use triangleclass
        use types
        use photon_vars, only : nxp, nyp, nzp
        use iarray,      only : tarray
        use constants,   only : xmax, ymax, zmax, in

        implicit none

        integer, intent(IN)  :: list(:)
        real,    intent(IN)  :: x, y, z
        real,    intent(OUT) :: distance
        integer, intent(OUT) :: tri

        real         :: eps=1d-10, det, invDet, u, v, t, min
        type(vector) :: v0, v1, v2, e1, e2, pvec, s, q, ray, origin
        integer      :: i, counter, loopdx
        logical      :: flag

        counter = 0
        ray    = vector(nxp, nyp, nzp)
        origin = vector(x-xmax, y-ymax, z-zmax)

        flag = .false.
        min  = huge(1.)

        loopdx = 0
        do i = 1, size(list)
            ! v0 = tarray(i)%vert(1) 
            ! v1 = tarray(i)%vert(2) 
            ! v2 = tarray(i)%vert(3) 

            v0 = tarray(list(i))%vert(1) 
            v1 = tarray(list(i))%vert(2) 
            v2 = tarray(list(i))%vert(3) 

            e1 = v1 - v0
            e2 = v2 - v0
            pvec = ray .cross. e2
            det = e1 .dot. pvec 

            if(abs(det) < eps)cycle
            invDet = 1./det
            s = origin - v0
            u = invDet * (s .dot. pvec)
            if(u < 0. .or. u > 1.)cycle
            q = s .cross. e1
            v = invDet * (ray .dot. q)
            if(v < 0.0 .or. (u + v) > 1.)cycle
            t = invDet * (e2 .dot. q)

            if(t > eps)then
                counter = counter + 1
                if(t < min)then
                    loopdx = i
                    min = t
                end if
                flag = .true.
            end if
        end do
        distance = min

        ! if(mod(real(counter), 2.) == 0)then
        !     ! print*,counter,'f'
        !     in = .false.
        ! else
        !     ! print*,counter,'t'
        !     in = .true.
        ! end if

        if(.not. flag)then
            distance = huge(1.)
            ! in = .false.
            tri = 0
        else
            tri = list(loopdx)
        end if
        ! tri = loopdx
    end subroutine ray_triangle_intersection


    real function wall_dist(celli, cellj, cellk, xcur, ycur, zcur, dir)
    !funtion that returns distance to nearest wall and which wall that is (x,y or z)
    !
    !
        use iarray,      only : xface, yface, zface
        use photon_vars, only : nxp, nyp, nzp

        implicit none

        real,    intent(INOUT) :: xcur, ycur, zcur
        logical, intent(INOUT) :: dir(:)
        integer, intent(INOUT) :: celli, cellj, cellk
        real                   :: dx, dy, dz

        if(nxp > 0.)then
            dx = (xface(celli+1) - xcur)/nxp
        elseif(nxp < 0.)then
            dx = (xface(celli) - xcur)/nxp
        elseif(nxp == 0.)then
            dx = 100000.
        end if

        if(nyp > 0.)then
            dy = (yface(cellj+1) - ycur)/nyp
        elseif(nyp < 0.)then
            dy = (yface(cellj) - ycur)/nyp
        elseif(nyp == 0.)then
            dy = 100000.
        end if

        if(nzp > 0.)then
            dz = (zface(cellk+1) - zcur)/nzp
        elseif(nzp < 0.)then
            dz = (zface(cellk) - zcur)/nzp
        elseif(nzp == 0.)then
            dz = 100000.
        end if

        wall_dist = min(dx, dy, dz)
        if(wall_dist < 0.)then
            print*,'dcell < 0.0 warning! ',wall_dist,dx,dy,dz,nxp,nyp,nzp
            error stop 1
        end if
        if(wall_dist == dx)dir=(/.TRUE., .FALSE., .FALSE./)
        if(wall_dist == dy)dir=(/.FALSE., .TRUE., .FALSE./)
        if(wall_dist == dz)dir=(/.FALSE., .FALSE., .TRUE./)
        if(.not.dir(1) .and. .not.dir(2) .and. .not.dir(3))print*,'Error in dir flag'
      
   end function wall_dist


    subroutine update_pos(xcur, ycur, zcur, celli, cellj, cellk, dcell, wall_flag, dir, delta)
    !routine that upates postions of photon and calls fresnel routines if photon leaves current voxel
    !
    !
        use photon_vars, only : nxp, nyp, nzp
        use iarray,      only : xface, yface, zface
        use utils,       only : str, red, bold, colour

        implicit none
      
      real,    intent(INOUT) :: xcur, ycur, zcur
      real,    intent(IN)    :: dcell, delta
      integer, intent(INOUT) :: celli, cellj, cellk
      logical, intent(IN)    :: wall_flag, dir(:)
      character(len=32)      :: tmp  
      
      if(wall_flag)then
      
         if(dir(1))then
            if(nxp > 0.)then
               xcur = xface(celli+1) + delta
            elseif(nxp < 0.)then
               xcur = xface(celli) - delta
            else
               print*,'Error in x dir in update_pos', dir, nxp, nyp, nzp
            end if
            ycur = ycur + nyp*dcell 
            zcur = zcur + nzp*dcell
         elseif(dir(2))then
            xcur = xcur + nxp*dcell
            if(nyp > 0.)then
               ycur = yface(cellj+1) + delta
            elseif(nyp < 0.)then
               ycur = yface(cellj) - delta
            else
               print*,'Error in y dir in update_pos', dir, nxp, nyp, nzp
            end if
            zcur = zcur + nzp*dcell
         elseif(dir(3))then
            xcur = xcur + nxp*dcell
            ycur = ycur + nyp*dcell 
            if(nzp > 0.)then
               zcur = zface(cellk+1) + delta
            elseif(nzp < 0.)then
               zcur = zface(cellk) - delta
            else
               print*,'Error in z dir in update_pos', dir, nxp, nyp, nzp
            end if
         else
            tmp = colour('Error in update_pos... '//str(dir), red, bold)
            error stop 1
         end if
      else
      
         xcur = xcur + nxp*dcell
         ycur = ycur + nyp*dcell 
         zcur = zcur + nzp*dcell
      
      end if


      if(wall_flag)then
         call update_voxels(xcur, ycur, zcur, celli, cellj, cellk)
      end if
      
    end subroutine update_pos


    subroutine update_voxels(xcur, ycur, zcur, celli, cellj, cellk)
    !updates the current voxel based upon position
    !
    !
        use iarray, only : xface, yface, zface
        use subs, only : find
        implicit none

        real,    intent(IN)    :: xcur, ycur, zcur
        integer, intent(INOUT) :: celli, cellj, cellk

        celli = find(xcur, xface) 
        cellj = find(ycur, yface)
        cellk = find(zcur, zface) 

    end subroutine update_voxels


subroutine reflect_refract(I, N, n1, n2, iseed, rflag)

        use types

        implicit none

        type(vector), intent(INOUT) :: I
        type(vector), intent(INOUT) :: N
        real,         intent(IN)    :: n1, n2
        integer,      intent(INOUT) :: iseed
        logical,      intent(OUT)   :: rflag

        real :: ran2

        rflag = .FALSE.

        if(ran2(iseed) <= fresnel(I, N, n1, n2))then
            call reflect(I, N)
            rflag = .true.
        else
            call refract(I, N, n1/n2)
        end if

    end subroutine reflect_refract


    subroutine reflect(I, N)
    !   get vector of reflected photon
    !
    !
        use types

        implicit none

        type(vector), intent(INOUT) :: I
        type(vector), intent(IN)    :: N

        type(vector) :: R

        R = I - 2. * (N .dot. I) * N
        I = R

    end subroutine reflect


    subroutine refract(I, N, eta)
    !   get vector of refracted photon
    !
    !
        use types

        implicit none

        type(vector), intent(INOUT) :: I
        type(vector), intent(IN)    :: N
        real,         intent(IN)    :: eta

        type(vector) :: T, Ntmp

        real :: c1, c2

        Ntmp = N

        c1 = (Ntmp .dot. I)
        if(c1 < 0.)then
            c1 = -c1
        else
            Ntmp = (-1.) * N
        end if
        c2 = sqrt(1. - (eta)**2 * (1.-c1**2))

        T = eta*I + (eta * c1 - c2) * Ntmp 

        I = T

    end subroutine refract


    function fresnel(I, N, n1, n2) result (tir)
    !   calculates the fresnel coefficents
    !
    !
        use types
        use ieee_arithmetic, only : ieee_is_nan

        implicit none

        real, intent(IN)         :: n1, n2
        type(vector), intent(IN) :: I, N

        real             ::  costt, sintt, sint2, cost2, tir, f1, f2

        costt = abs(I .dot. N)

        sintt = sqrt(1. - costt * costt)
        sint2 = n1/n2 * sintt
        if(sint2 > 1.)then
            tir = 1.0
            return
        elseif(costt == 1.)then
            tir = 0.
            return
        else
            sint2 = (n1/n2)*sintt
            cost2 = sqrt(1. - sint2 * sint2)
            f1 = abs((n1*costt - n2*cost2) / (n1*costt + n2*cost2))**2
            f2 = abs((n1*cost2 - n2*costt) / (n1*cost2 + n2*costt))**2

            tir = 0.5 * (f1 + f2)
        if(ieee_is_nan(tir) .or. tir > 1. .or. tir < 0.)print*,'TIR: ', tir, f1, f2, costt,sintt,cost2,sint2
            return
        end if
    end function fresnel

end module inttau2



                    ! interpolatedNormal = bary[0] * p1_normal + bary[2] * p2_normal + bary[1] * p3_normal;

                    ! v1 = tarray(tri)%vert(1)
                    ! v2 = tarray(tri)%vert(2)
                    ! v3 = tarray(tri)%vert(3)

                    ! p = vector(xcur, ycur, zcur)

                    ! bary = barycentric(v1, v2, v3, p)

                    ! intrpNormal = bary%x * tarray(tri)%norms(1) + bary%y * tarray(tri)%norms(2) + bary%z * tarray(tri)%norms(3) 
                    ! incd = vector(nxp, nyp, nzp)

                    ! if(in)then
                    !     intrpNormal = vector(-intrpNormal%x, -intrpNormal%y, -intrpNormal%z)
                    !     n1 = 1.52
                    !     n2 = 1.0
                    ! else
                    !     n1 = 1.0
                    !     n2 = 1.52
                    ! end if
                    ! call reflect_refract(incd, intrpNormal, n1, n2, iseed, rflag)

                    ! ! print*,nxp,nyp,nzp

                    ! incd = normal(incd)
                    ! nxp = incd%x
                    ! nyp = incd%y
                    ! nzp = incd%z

                    ! phi = atan2(nyp, nxp)
                    ! sinp = sin(phi)
                    ! cosp = cos(phi)

                    ! cost = nzp
                    ! sint = sqrt(1.-cost*cost)
                    ! ! print*,nxp,nyp,nzp
                    ! ! stop
                    ! if(.not.rflag)then