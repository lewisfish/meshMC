module subs

implicit none

save

    contains

    type(vector) function barycentric(a, b, c, p)

        use types

        implicit none

        type(vector), intent(IN) :: a, b, c, p

        type(vector) :: v0, v1, v2
        real         :: d00, d01, d11, d20, d21, invDenom, v, u, w

        v0 = b - a
        v1 = c - a
        v2 = p - a
        d00 = v0 .dot. v0
        d01 = v0 .dot. v1
        d11 = v1 .dot. v1
        d20 = v2 .dot. v0
        d21 = v2 .dot. v1
        invDenom = 1. / ((d00 * d11) - (d01 * d01))
        v = (d11 * d20 - d01 * d21) * invDenom
        w = (d00 * d21 - d01 * d20) * invDenom
        u = 1. - v - w

        barycentric = vector(u, v, w)

    end function barycentric


    integer function find(val, a)
    !searchs for bracketing indicies for a value val in an array a
    !
    !
        implicit none

        real, intent(IN) :: val, a(:)
        integer          :: n, lo, mid, hi
        logical          :: ascnd

        n = size(a)
        ascnd = (a(n) >= a(1))
        lo = 0
        hi = n+1
        do
            if (hi-lo <= 1) exit
            mid = (hi+lo)/2
            if (ascnd .eqv. (val >= a(mid))) then
                lo = mid
            else
                hi = mid
            end if
        end do

        if (val == a(1)) then
            find = 1
        else if (val == a(n)) then
            find = n-1
        else if(ascnd.and. (val > a(n) .or. val < a(1))) then
            find = -1
        else if(.not.ascnd.and. (val < a(n) .or. val > a(1))) then
            find = -1
        else
            find = lo
        end if

    end function find


        ! subroutine add_triangle(i,j,k,l)

        !     use iarray, only : trilist

        !     implicit none

        !     integer, intent(IN) :: i, j, k, l
        !     integer, allocatable :: wk(:)

        !     if(.not. allocated(trilist(i,j,k)%triangles))then
        !         allocate(trilist(i,j,k)%triangles(10))
        !         trilist(i,j,k)%triangles(:) = 0
        !     end if

        !     ! print*,trilist(i,j,k)%iptr 
        !     if(trilist(i,j,k)%iptr > size(trilist(i,j,k)%triangles))then
        !         allocate(wk(size(trilist(i,j,k)%triangles)+3))
        !         wk = 0
        !         wk(1:size(trilist(i,j,k)%triangles)) = trilist(i,j,k)%triangles
        !         call move_alloc(wk, trilist(i,j,k)%triangles)
        !     end if
        !     trilist(i,j,k)%triangles(trilist(i,j,k)%iptr) = l 
        !     trilist(i,j,k)%iptr = trilist(i,j,k)%iptr + 1

        !     ! print*,trilist(i,j,k)%triangles(:),trilist(i,j,k)%iptr
        ! end subroutine add_triangle


        subroutine directory
        !  subroutine defines vars to hold paths to various folders   
        !   
        !   
            use constants, only : cwd, homedir, fileplace, resdir

            implicit none

            !get current working directory

            call get_environment_variable('PWD', cwd)

            ! get 'home' dir from cwd
            homedir = trim(cwd(1:len(trim(cwd))-3))
            ! get data dir
            fileplace = trim(homedir)//'data/'
            ! get res dir
            resdir=trim(homedir)//'res/'

        end subroutine directory


        subroutine zarray

            use iarray
            use types, only : assignment(=)

            !sets all arrays to zero
            implicit none

            xface = 0.
            yface = 0.
            zface = 0.
            ! rhokap = 0.
            jmeanout = 0.
            jmeanoutglobal = 0.
            jmeanin = 0.
            jmeaninglobal = 0.
            ! trilist = 0
        end subroutine zarray


        subroutine alloc_array(numproc, id)
        !  subroutine allocates allocatable arrays
        !   
        !   
            use iso_fortran_env, only : int64
            use utils,           only : str, mem_free
            use constants,       only : nxg, nyg, nzg
            use iarray

            implicit none

            integer, intent(IN) :: numproc, id

            integer(int64) :: limit, cnt, i


            limit = mem_free()
            cnt = 0_int64

            allocate(xface(nxg+1))
            inquire(iolength=i)xface(:)
            call chck_mem(cnt, i, limit, 'xface', numproc)

            allocate(yface(nyg+1))
            inquire(iolength=i)yface
            call chck_mem(cnt, i, limit, 'yface', numproc)

            allocate(zface(nzg+1))
            inquire(iolength=i)zface
            call chck_mem(cnt, i, limit, 'zface', numproc)

            ! allocate(trilist(nxg, nyg, nzg))

            ! allocate(rhokap(nxg,nyg,nzg))
            ! inquire(iolength=i)rhokap
            ! call chck_mem(cnt, i, limit, 'rhokap', numproc)

            allocate(jmeanin(nxg,nyg,nzg))
            inquire(iolength=i)jmeanin
            call chck_mem(cnt, i, limit, 'jmeanin', numproc)

            allocate(jmeaninglobal(nxg,nyg,nzg))
            inquire(iolength=i)jmeaninglobal
            call chck_mem(cnt, i, limit, 'jmeaninglobal', numproc)

            allocate(jmeanout(nxg,nyg,nzg))
            inquire(iolength=i)jmeanout
            call chck_mem(cnt, i, limit, 'jmeanout', numproc)

            allocate(jmeanoutglobal(nxg,nyg,nzg))
            inquire(iolength=i)jmeanoutglobal
            call chck_mem(cnt, i, limit, 'jmeanoutglobal', numproc)

            if(id == 0)print'(A,1X,F5.2,A)','allocated:',dble(cnt)/dble(limit)*100.d0,' % of total RAM'
        end subroutine alloc_array


        subroutine chck_mem(cur, new, limit, name, numproc)
        !routine to check if the system has enough RAM available in order to run the simulation
        !cur: current memory assigned, new: new memory to be assigned
        !limit: the limit of RAM available, name: name of array to be assigned, numproc: processor #

            use iso_fortran_env, only : int64
            use utils,           only : str

            implicit none

            integer(int64), intent(IN)    :: new, limit
            integer(int64), intent(INOUT) :: cur 
            integer,        intent(IN)    :: numproc
            character(*),   intent(IN)    :: name

            ! integer :: error

            cur = cur + new * numproc
            if(cur > limit)then
                print*,'Need '//str(cur-limit)//' more memory to run. '//name
                ! call mpi_finalize(error)
                stop
            end if
        end subroutine chck_mem
end MODULE subs