module types

    implicit none

    type :: vector
        real :: x, y, z
    end type vector

    type :: ivec
        integer :: x, y, z
    end type ivec

    type point
        integer :: x, y
    end type point

    type trigrid
        integer, allocatable :: triangles(:)
        integer              :: iptr = 1
    end type trigrid

    interface assignment(=)
        module procedure equal_sub
    end interface
    
    interface operator (.dot.)
        module procedure vecDot
    end interface

    interface operator (.cross.)
        module procedure vecCross
        module procedure ivecCross
    end interface

    interface operator (-)
        module procedure vecSub
        module procedure ivecSub            
    end interface

    interface operator (*)
        module procedure vecMulA
        module procedure ivecMulA
        module procedure vecMulB
    end interface

    interface operator (/)
        module procedure vecdivA
        module procedure ivecdivA
        module procedure vecdivB
        module procedure ivecDivAI
    end interface

    interface operator (+)
        module procedure vecAdd
        module procedure ivecAdd
    end interface

    interface operator(-)
        module procedure pointSub
    end interface

    interface operator(+)
       module procedure pointAdd
    end interface

    interface operator(*)
       module procedure pointmultscal
    end interface

    interface swap
        module procedure swap_I
        module procedure swap_R
    end interface

    private
    public :: vector, ivec, point,trigrid
    public :: operator(.dot.), operator(.cross.), operator(-), operator(*), operator(+), operator(/), assignment(=)
    public :: magnitude, normal, swap

    contains
        subroutine equal_sub(a, b)

            implicit none
            
            type(trigrid), intent(INOUT) :: a(:,:,:)
            integer,       intent(IN)    :: b

            integer :: i, j, k, l

            do i = 1, size(a, 1)
                do j = 1, size(a, 2)
                    do k = 1, size(a, 3)
                        do l = 1, size(a(i,j,k)%triangles)
                            a(i,j,k)%triangles(l) = b
                        end do
                    end do
                end do
            end do
        end subroutine equal_sub


        subroutine swap_I(a, b)

            implicit none

            integer, intent(INOUT) :: a, b
            integer                :: tmp

            tmp = a
            a = b
            b = tmp

        end subroutine swap_I


        subroutine swap_R(a, b)

            implicit none

            real, intent(INOUT) :: a, b
            real                :: tmp

            tmp = a
            a = b
            b = tmp

        end subroutine swap_R


        type(point) function pointSub(a, b)

           implicit none

           type(point), intent(IN) :: a, b

           pointSub = point(a%x - b%x, a%y - b%y)

        end function pointSub


        type(point) function pointAdd(a, b)

            implicit none

            type(point), intent(IN) :: a, b

            pointAdd = point(a%x + b%x, a%y + b%y)

        end function pointAdd


        type(point) function pointmultscal(a, b)

            implicit none

            type(point), intent(IN) :: a
            real,        intent(IN) :: b

            pointmultscal = point(int(a%x * b), int(a%y * b))

        end function pointmultscal


        type(vector) function vecMulA(a, b)

            implicit none

            type(vector), intent(IN) :: a
            real,         intent(IN) :: b

            vecMulA = vector(a%x * b, a%y * b, a%z * b)

        end function vecMulA


        type(ivec) function ivecMulA(a, b)

            implicit none

            type(ivec), intent(IN) :: a
            real,         intent(IN) :: b

            ivecMulA = ivec(int(a%x * b), int(a%y * b), int(a%z * b))

        end function ivecMulA


        type(vector) function vecMulB(a, b)

            implicit none

            type(vector), intent(IN) :: b
            real,         intent(IN) :: a

            vecMulB = vector(a * b%x, a * b%y, a * b%z)

        end function vecMulB


        type(vector) function vecDivA(a, b)

            implicit none

            type(vector), intent(IN) :: a
            real,         intent(IN) :: b

            vecDivA = vector(a%x / b, a%y / b, a%z / b)

        end function vecDivA


        type(ivec) function ivecDivA(a, b)

            implicit none

            type(ivec), intent(IN) :: a
            real,         intent(IN) :: b

            ivecDivA = ivec(int(a%x / b), int(a%y / b), int(a%z / b))

        end function ivecDivA


        type(ivec) function ivecDivAI(a, b)

            implicit none

            type(ivec), intent(IN) :: a
            integer,    intent(IN) :: b

            ivecDivAI = ivec(a%x / b, a%y / b, a%z / b)

        end function ivecDivAI


        type(vector) function vecDivB(a, b)

            implicit none

            type(vector), intent(IN) :: b
            real,         intent(IN) :: a

            vecDivB = vector(a / b%x, a / b%y, a / b%z)

        end function vecDivB


        type(vector) function vecSub(a, b)

            implicit none

            type(vector), intent(IN) :: a, b

            vecSub = vector(a%x - b%x, a%y - b%y , a%z - b%z)

        end function vecSub


            type(ivec) function ivecSub(a, b)

            implicit none

            type(ivec), intent(IN) :: a, b

            ivecSub = ivec(a%x - b%x, a%y - b%y , a%z - b%z)

        end function ivecSub


        type(vector) function vecAdd(a, b)

            implicit none

            type(vector), intent(IN) :: a, b

            vecAdd = vector(a%x + b%x, a%y + b%y , a%z + b%z)

        end function vecAdd


        type(ivec) function ivecAdd(a, b)

            implicit none

            type(ivec), intent(IN) :: a, b

            ivecAdd = ivec(a%x + b%x, a%y + b%y , a%z + b%z)

        end function ivecAdd


        real function vecDot(a, b)

            implicit none

            type(vector), intent(IN) :: a, b
            vecDot = (a%x * b%x) + (a%y * b%y) + (a%z * b%z)

        end function vecDot


        type(vector) function vecCross(a, b)

            implicit none

            type(vector), intent(IN) :: a, b
            real                     :: i, j, k

            i = (a%y * b%z) - (b%y * a%z)
            j = (a%x * b%z) - (b%x * a%z)
            k = (a%x * b%y) - (b%x * a%y)

            j = -j

            vecCross = vector(i, j, k)

        end function vecCross


        type(ivec) function ivecCross(a, b)

            implicit none

            type(ivec), intent(IN) :: a, b
            integer                :: i, j, k

            i = (a%y * b%z) - (b%y * a%z)
            j = -((a%x * b%z) - (b%x * a%z))
            k = (a%x * b%y) - (b%x * a%y)

            j = -j

            ivecCross = ivec(i, j, k)

        end function ivecCross


        real function magnitude(this)

            implicit none

            type(vector) :: this

            magnitude = sqrt(this%x**2 + this%y**2 + this%z**2)

        end function magnitude

        type(vector) function normal(this)

            implicit none

            type(vector), intent(in) :: this
            real :: mag

            mag = magnitude(this)
            if(mag == 0.)stop "magnitude = 0"
            normal = vector(this%x/mag, this%y/mag, this%z/mag)

        end function normal
end module types