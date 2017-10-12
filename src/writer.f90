MODULE writer_mod

implicit none
save

private
public :: writer

interface write_binary
    ! module procedure write_binaryR2
    module procedure write_binaryR3
    module procedure write_binaryR4
    ! module procedure write_binaryI3
end interface

CONTAINS
    subroutine writer(nphotons, numproc, file)

        use iarray,          only : jmeaninGLOBAL, jmeanoutGLOBAL
        use utils,           only : str

        implicit none

        integer,      intent(IN) :: nphotons
        character(*), intent(IN) :: file
        integer :: numproc
        character(len=256) :: filename

        filename = 'jmean/jmeanin-'//trim(file)//'.raw'
        call write_binary(trim(filename), jmeaninGLOBAL(:,:,:))
        call write_bov(trim(filename), trim(file), '-in')
        print*,trim(filename)

        filename = 'jmean/jmeanout-'//trim(file)//'.raw'
        call write_binary(trim(filename), jmeanoutGLOBAL(:,:,:))
        call write_bov(trim(filename), trim(file), '-out')
        print*,trim(filename)
    end subroutine writer

    subroutine write_bov(filename, name, suffix)

        use constants, only : fileplace, xmax, ymax, zmax, nxg, nyg, nzg
        use utils,     only : str

        implicit none

        character(*), intent(IN) :: filename, name, suffix

        integer :: u

        open(newunit=u, file=trim(name)//suffix//'.bov')
        write(u, *)"DATA_FILE: "//trim(filename)
        write(u, *)"DATA_SIZE: "//str(nxg)//' '//str(nyg)//' '//str(nzg)
        write(u, *)"DATA_FORMAT: DOUBLE"
        write(u, *)"VARIABLE: "//trim(name)//trim(suffix)
        write(u, *)"DATA_ENDIAN: LITTLE"
        write(u, *)"CENTERING: ZONAL"
        write(u, *)"BRICK_ORIGIN: "//str(-xmax)//' '//str(-ymax)//' '//str(-zmax)
        write(u, *)"BRICK_SIZE: "//str(2.*xmax)//' '//str(2.*ymax)//' '//str(2.*zmax)
        close(U)


    end subroutine write_bov


    ! subroutine write_binaryR2(filename, array)

    !     use constants, only : fileplace

    !     implicit none

    !     character(len=*), intent(IN) :: filename
    !     real,             intent(IN) :: array(:,:)
        
    !     integer(kind=8) :: u, i

    !     inquire(iolength=i)array
    !     open(newunit=u,file=trim(fileplace)//trim(filename),access='stream',status='REPLACE',form='unformatted')
    !     write(u) array
    !     close(u)

    ! end subroutine write_binaryR2


    subroutine write_binaryR3(filename, array)

        use constants, only : fileplace

        implicit none

        character(len=*), intent(IN) :: filename
        real,             intent(IN) :: array(:,:,:)
        
        integer(kind=8) :: u, i

        inquire(iolength=i)array
        open(newunit=u,file=trim(fileplace)//trim(filename),access='stream',status='REPLACE',form='unformatted')
        write(u) array
        close(u)

    end subroutine write_binaryR3


    subroutine write_binaryR4(filename, array)

        use constants, only : fileplace

        implicit none

        character(len=*), intent(IN) :: filename
        real,             intent(IN) :: array(:,:,:,:)
        
        integer :: u, i

        inquire(iolength=i)array
        open(newunit=u,file=trim(fileplace)//trim(filename),access='direct',status='REPLACE',form='unformatted',&
        recl=i)
        write(u,rec=1) array
        close(u)

    end subroutine write_binaryR4


    ! subroutine write_binaryI3(filename, array)

    !     use constants, only : fileplace

    !     implicit none

    !     character(len=*), intent(IN) :: filename
    !     integer,          intent(IN) :: array(:,:,:)
        
    !     integer :: u, i

    !     inquire(iolength=i)array
    !     open(newunit=u,file=trim(fileplace)//trim(filename),access='direct',status='REPLACE',form='unformatted',&
    !     recl=i)
    !     write(u,rec=1) array
    !     close(u)

    ! end subroutine write_binaryI3

end MODULE writer_mod
