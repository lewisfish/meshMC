Module shapeclass

    use types

    implicit none

    type, abstract :: shape

    Contains

    ! procedure :: info => info_fn

    end type shape

    ! private :: info_fn

    Contains

    ! subroutine info_fn(this)

    !     class(shape), intent(IN) :: this

    !     print*,'This is a shape. Colour: ', this%colour
    ! end subroutine info_fn

end Module shapeclass