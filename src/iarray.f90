MODULE iarray
!
!  Contains all array var names.
!
    use triangleclass
    use types, only : trigrid
    
    implicit none
    save

    real,           allocatable :: xface(:), yface(:), zface(:)
    type(triangle), allocatable :: tarray(:)
    real,           allocatable :: jmeanin(:,:,:), jmeaninGLOBAL(:,:,:)
    real,           allocatable :: jmeanout(:,:,:), jmeanoutGLOBAL(:,:,:)

    ! type(trigrid),  allocatable :: trilist(:,:,:)
    ! real, allocatable :: rhokap(:,:,:)

end MODULE iarray
