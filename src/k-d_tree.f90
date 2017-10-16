Module kdtree_mod

    use types,         only : vector
    use triangleclass, only : triangle
    
    implicit none

    type :: node
        integer              :: level
        real                 :: mindim(3), maxdim(3)
        logical              :: isLeaf=.false.
        type(node), pointer  :: ln => null(), rn =>null()
        integer, allocatable :: list(:)
    end type node    

    !defintion of kdtree for program
    type(node) :: kdtree

    public
    private :: get_bbox

    contains

        subroutine get_bbox(triarray, mind, maxd)
        ! gets bounding box for list of triangles
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
            if(depth == 8)then
                newNode%isLeaf = .true.
                allocate(newNode%list(size(trilist)))
                newNode%list = trilist
            else
                !recurse
                !divide down middle of bbox
                diff = (maxdim(axis) - mindim(axis)) / 2.
                up = maxdim(axis)
                low = mindim(axis)
                mid = up - diff

                ! mid = get_median(trilist, axis)

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
                        ! if(maxdim(axis) >= mid)then
                        !     right(i) = trilist(i)
                        !     left(i) = trilist(i)
                        ! else
                            left(i) = trilist(i)
                        ! end if
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


        function sort(a)
            !  return on call the array that is passed, in order from lowest to highest value
            !  preserves order of array passed in. Slow sort, good for only small arrays
            !  INPUT:
            !           a     real    the array that is to be sorted         
            !
            !  OUTPUT:  a     real    the sorted array
            implicit none

            real, intent(in) :: a(:)
            real             :: sort(size(a))
            integer          :: i, minIndex
            real             :: temp

            sort = a 

            do i = 1, size(sort)
               minIndex = MINLOC(sort(i:), 1) + i - 1
               if(sort(i) .gt. sort(minIndex))then
                  temp = sort(i)
                  sort(i) = sort(minIndex)
                  sort(minIndex) = temp
               end if
            end do

        end function sort


        real function get_median(list, axis)

            use iarray, only : tarray

            implicit none

            integer, intent(IN) :: list(:), axis

            integer :: i,j, counter
            real    :: array(size(list)*3)

            counter = 1

            do i = 1, size(list)
                do j = 1, 3
                    if(axis==1)array(counter) = tarray(list(i))%vert(j)%x
                    if(axis==2)array(counter) = tarray(list(i))%vert(j)%y
                    if(axis==3)array(counter) = tarray(list(i))%vert(j)%z
                    counter = counter + 1
                end do
            end do

            if(mod(size(array), 2) == 0)then
                counter = size(array)/2.
                if(array(counter+1) == array(counter))then
                    get_median = array(counter)
                else
                    get_median = (array(counter+1) - array(counter))/2.
                end if
            else
                get_median =  array(nint(real(size(array))/2.))
            end if
        end function get_median

        logical function intersect_bbox(origin, ray, mindim, maxdim)
        !   calculates if ray intersects with bounding box

            use types, only : vector, swap, operator(/)

            implicit none

            type(vector), intent(IN) :: origin, ray
            real, intent(IN) :: mindim(:), maxdim(:)
            type(vector) :: invRay
            real :: tmin, tmax, tymin, tymax, tzmin, tzmax


            invRay = 1./ray

            if(invRay%x >= 0)then
                tmin = (mindim(1) - origin%x) * invRay%x
                tmax = (maxdim(1) - origin%x) * invRay%x
            else
                tmin = (maxdim(1) - origin%x) * invRay%x
                tmax = (mindim(1) - origin%x) * invRay%x
            end if

            if(tmin > tmax)call swap(tmin, tmax)

            if(invRay%y >= 0)then
                tymin = (mindim(2) - origin%y) * invRay%y
                tymax = (maxdim(2) - origin%y) * invRay%y
            else
                tymin = (maxdim(2) - origin%y) * invRay%y
                tymax = (mindim(2) - origin%y) * invRay%y
            end if
            if(tymin > tymax)call swap(tymin, tymax)

            if((tmin > tymax) .or. (tymin > tmax))then
                intersect_bbox = .FALSE.    
                return
            end if

            if(tymin > tmin)tmin = tymin
            if(tymax < tmax)tmax = tymax

            if(invRay%z >= 0)then
                tzmin = (mindim(3) - origin%z) * invRay%z
                tzmax = (maxdim(3) - origin%z) * invRay%z
            else
                tzmin = (maxdim(3) - origin%z) * invRay%z
                tzmax = (mindim(3) - origin%z) * invRay%z
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
        !   writes out bounding box for each node

            implicit none

            integer, intent(INOUT) :: u
            type(node) :: tmp

            write(u,*) 'v', tmp%mindim
            write(u,*) 'v', tmp%maxdim(1),tmp%mindim(2),tmp%mindim(3)
            write(u,*) 'v', tmp%maxdim(1),tmp%maxdim(2),tmp%mindim(3)
            write(u,*) 'v', tmp%mindim(1),tmp%maxdim(2),tmp%mindim(3)
            write(u,*) 'f ', '1 2 3'
            write(u,*) 'f ', '1 3 4'

            ! write(u,*) tmp%mindim
            ! write(u,*) tmp%mindim(1),tmp%mindim(2),tmp%maxdim(3)
            ! write(u,*) tmp%maxdim(1),tmp%mindim(2),tmp%maxdim(3)
            ! write(u,*) tmp%maxdim(1),tmp%maxdim(2),tmp%maxdim(3)
            ! write(u,*) tmp%mindim(1),tmp%maxdim(2),tmp%maxdim(3)
            ! write(u,*) tmp%mindim(1),tmp%mindim(2),tmp%maxdim(3)
            close(u)
            stop
            ! write(u,*) ''
            ! write(u,*) ''

            if(.not. tmp%isLeaf)then
                call write_bbox(tmp%rn, u)
                call write_bbox(tmp%ln, u)
            end if
        end subroutine write_bbox
end module kdtree_mod