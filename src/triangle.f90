module triangleclass

    use shapeclass

    implicit none

    type, extends(shape) :: triangle
        type(vector) :: vert(3), uvs(3), norms(3)        
    end type triangle

    contains

        subroutine bbox(tri, min, max)

            implicit none

            type(vector), intent(OUT)  :: min, max
            type(triangle), intent(IN) :: tri
            integer :: i

            min = vector(huge(1.), huge(1.), huge(1.))
            max = vector(-huge(1.), -huge(1.), -huge(1.))
            do i = 1, 3
                if(tri%vert(i)%x < min%x)min%x = tri%vert(i)%x
                if(tri%vert(i)%y < min%y)min%y = tri%vert(i)%y
                if(tri%vert(i)%z < min%z)min%z = tri%vert(i)%z

                if(tri%vert(i)%x > max%x)max%x = tri%vert(i)%x
                if(tri%vert(i)%y > max%y)max%y = tri%vert(i)%y
                if(tri%vert(i)%z > max%z)max%z = tri%vert(i)%z
            end do

        end subroutine bbox

end module triangleclass
