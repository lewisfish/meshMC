module triangleclass

    use shapeclass

    implicit none

    type, extends(shape) :: triangle
        type(vector) :: vert(3), uvs(3), norms(3)        
    end type triangle

    ! contains


end module triangleclass
