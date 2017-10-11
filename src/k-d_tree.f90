Module kdtree_mod

    use types,         only : vector
    use triangleclass, only : triangle
    
    implicit none

    type :: node
        integer              :: level
        real                 :: mindim(3), maxdim(3)
        logical              :: flag=.false.
        type(node), pointer  :: ln => null(), rn =>null()
        integer, allocatable :: list(:)
    end type node    

    type(node) :: kdtree

    contains

        subroutine get_bbox(triarray, mind, maxd)

            use types,         only : vector
            use triangleclass, only : triangle
            use iarray,        only : tarray

            implicit none

            integer, intent(IN)  :: triarray(:)
            real,    intent(OUT) :: mind(:), maxd(:)

            integer :: i, j
            real    :: x, y, z

            mind = [huge(1.), huge(1.), huge(1.)]
            maxd = [-huge(1.), -huge(1.), -huge(1.)]

            do i = 1, size(triarray)
                do j = 1, 3
                    x = tarray(triarray(i))%vert(j)%x
                    y = tarray(triarray(i))%vert(j)%y
                    z = tarray(triarray(i))%vert(j)%z
                    if(x < mind(1))mind(1) = x
                    if(y < mind(2))mind(2) = y
                    if(z < mind(3))mind(3) = z

                    if(x > maxd(1))maxd(1) = x
                    if(y > maxd(2))maxd(2) = y
                    if(z > maxd(3))maxd(3) = z
                end do
            end do

        end subroutine get_bbox


        recursive function buildTree(trilist, depth) result(tree)

            use triangleclass, only : triangle
            use types,         only : vector

            implicit none

            type(node) :: tree
            integer, intent(IN)  :: trilist(:), depth
            
            real :: mindim(3), maxdim(3), up, low, mid, diff
            integer, allocatable :: left(:), right(:), leftlst(:), rightlst(:)
            integer :: i, axis, counter
            type(node) :: newNode

            newNode%level = depth

            axis = mod(depth, 3) + 1 !so account for fortran arrays stating at 1...

            call get_bbox(trilist, mindim, maxdim)!get bbox for all triangles
            newNode%mindim = mindim
            newNode%maxdim = maxdim
            
            !stop tree at depth of 3 and create leafs with data
            if(depth == 9)then
                newNode%flag = .true.
                allocate(newNode%list(size(trilist)))
                newNode%list = trilist
            else
                !recurse
                !divide down middle of bbox
                diff = (maxdim(axis) - mindim(axis)) / 2.
                up = maxdim(axis)
                low = mindim(axis)
                mid = up - diff

                ! print*,low,mid,up,depth

                !the section is ugly as fuck
                !gets list of triangles that belong either to left or right
                allocate(left(size(trilist)), right(size(trilist)))
                left = 0
                right = 0

                do i = 1, size(trilist)
                    call get_bbox(trilist(i:i), mindim, maxdim)
                    if(mindim(axis) >= mid)then
                        !right
                        right(i) = trilist(i)
                    else
                        if(maxdim(axis) >= mid)then
                            right(i) = trilist(i)
                            left(i) = trilist(i)
                        else
                            left(i) = trilist(i)
                        end if
                    end if
                end do

                counter = 0
                do i = 1, size(left)
                    if(left(i) > 0)counter = counter + 1
                end do

                allocate(leftlst(counter))
                counter = 0
                do i = 1, size(left)
                    if(left(i) > 0)then
                        counter = counter + 1
                        leftlst(counter) = left(i)
                    end if
                end do

                counter = 0
                do i = 1, size(right)
                    if(right(i) > 0)counter = counter + 1
                end do

                allocate(rightlst(counter))
                counter = 0
                do i = 1, size(right)
                    if(right(i) > 0)then
                        counter = counter + 1
                        rightlst(counter) = right(i)
                    end if
                end do

                allocate(newNode%rn, newNode%ln)
                newNode%rn = buildtree(rightlst, depth + 1)
                newNode%ln = buildtree(leftlst, depth + 1)
            end if

            tree = newNode
        end function buildtree


        logical function intersect_bbox(origin, ray, mindim, maxdim)

            use types, only : vector, swap

            implicit none

            type(vector), intent(IN) :: origin, ray
            real, intent(IN) :: mindim(:), maxdim(:)

            real :: tmin, tmax, tymin, tymax, tzmin, tzmax

            if(1./ray%x >= 0)then
                tmin = (mindim(1) - origin%x) / ray%x
                tmax = (maxdim(1) - origin%x) / ray%x
            else
                tmin = (maxdim(2) - origin%y) / ray%y
                tmax = (mindim(2) - origin%y) / ray%y
            end if

            if(tmin > tmax)call swap(tmin, tmax)

            if(1./ray%y >= 0)then
                tymin = (mindim(2) - origin%y) / ray%y
                tymax = (maxdim(2) - origin%y) / ray%y
            else
                tymin = (maxdim(2) - origin%y) / ray%y
                tymax = (mindim(2) - origin%y) / ray%y
            end if
            if(tymin > tymax)call swap(tymin, tymax)

            if((tmin > tymax) .or. (tymin > tmax))then
                intersect_bbox = .FALSE.    
                return
            end if

            if(tymin > tmin)tmin = tymin
            if(tymax < tmax)tmax = tymax

            if(1./ray%z >= 0)then
                tzmin = (mindim(3) - origin%z) / ray%z
                tzmax = (maxdim(3) - origin%z) / ray%z
            else
                tzmin = (maxdim(3) - origin%z) / ray%z
                tzmax = (mindim(3) - origin%z) / ray%z
            end if

            if(tzmin > tzmax)call swap(tzmin, tzmax)

            if((tmin > tzmax) .or. (tzmin > tmax))then
                intersect_bbox = .false.
                return
            end if

            if(tzmin > tmin)tmin = tzmin
            if(tzmax < tmax)tmax = tzmax

            if(tmax > 0)then
                intersect_bbox = .true.
                return
            else
                intersect_bbox = .false.
                return
            end if
        end function intersect_bbox


        recursive subroutine write_bbox(tmp, u)

            implicit none

            integer, intent(INOUT) :: u
            type(node) :: tmp


            write(u,*) tmp%mindim
            write(u,*) tmp%maxdim(1),tmp%mindim(2),tmp%mindim(3)
            write(u,*) tmp%maxdim(1),tmp%maxdim(2),tmp%mindim(3)
            write(u,*) tmp%mindim(1),tmp%maxdim(2),tmp%mindim(3)
            write(u,*) tmp%mindim
            write(u,*) tmp%mindim(1),tmp%mindim(2),tmp%maxdim(3)
            write(u,*) tmp%maxdim(1),tmp%mindim(2),tmp%maxdim(3)
            write(u,*) tmp%maxdim(1),tmp%maxdim(2),tmp%maxdim(3)
            write(u,*) tmp%mindim(1),tmp%maxdim(2),tmp%maxdim(3)
            write(u,*) tmp%mindim(1),tmp%mindim(2),tmp%maxdim(3)
            write(u,*) ''
            write(u,*) ''
            if(.not. tmp%flag)then
                call write_bbox(tmp%rn, u)
                call write_bbox(tmp%ln, u)
            end if

        end subroutine write_bbox


end module kdtree_mod
! program testing

!     use triangleclass, only : triangle
!     use obj_reader,    only : read_obj
!     use types,         only : vector
!     use kdtree_mod,          only : node, buildTree, tarray
!     implicit none
    

!     type(node) :: root
!     integer :: id = 0, i
!     integer, allocatable :: listtri(:)

!     call read_obj('res/gourd.obj', tarray, id, .FALSE.)

!     allocate(listtri(size(tarray)))

!     listtri = 0

!     do i = 1, size(tarray)
!         listtri(i) = i
!     end do


!     root = buildtree(listtri, 0)
!     print*,root%level
!     print*,root%mindim
!     print*,root%maxdim
!     print*,''
!     print*,root%ln%level
!     print*,root%ln%mindim
!     print*,root%ln%maxdim
!     print*,''
!     print*,root%ln%ln%level
!     print*,root%ln%ln%mindim
!     print*,root%ln%ln%maxdim
!     print*,''
!     print*,root%ln%ln%ln%level
!     print*,root%ln%ln%ln%mindim
!     print*,root%ln%ln%ln%maxdim
!     print*,''
!     print*,root%ln%ln%ln%ll%level
!     print*,root%ln%ln%ln%ll%list
!     print*,''
!     print*,root%ln%ln%ln%ll%level
!     print*,root%ln%ln%ln%rl%list
!     ! call build_tree(tarray, root)




! end program testing