program mcpolar

    use mpi

    !shared data
    use utils, only : green, blue, bold, colour, red, white_b, black, str
    use constants
    use photon_vars
    use iarray
    use opt_prop
    use triangleclass, only : triangle
    
    !subroutines
    use subs
    use gridset_mod
    use sourceph_mod
    use inttau2
    use ch_opt
    use stokes_mod
    use writer_mod
    use obj_reader, only : read_obj
    use kdtree_mod, only : node, buildTree, kdtree, write_bbox

    implicit none

    integer              :: iseed, j, xcell, ycell, zcell, nphotons, u
    integer, allocatable :: listtri(:)
    logical              :: tflag
    double precision     :: nscatt
    real                 :: ran, delta, start, finish, ran2, start2, finish2, tmp
    character(len=512)   :: meshfile, meshname

    integer :: error,  id, numproc

    call MPI_init(error)

    call MPI_Comm_size(MPI_COMM_WORLD, numproc, error)

    call MPI_Comm_rank(MPI_COMM_WORLD, id, error)


    !set directory paths
    call directory()

    !allocate and set arrays to 0
    call alloc_array(numproc, id)
    call zarray()

    !**** Read in parameters from the file input.params
    open(10,file=trim(resdir)//'input.params',status='old')
    read(10,*) nphotons
    read(10,*) xmax
    read(10,*) ymax
    read(10,*) zmax
    read(10,*) beam
    read(10,*) meshname
    close(10)

    meshname = adjustl(meshname)

    meshfile = trim(meshname)
    call read_obj(trim(meshfile), tarray, id, .FALSE.)

    allocate(listtri(size(tarray)))

    listtri = 0

    do j = 1, size(tarray)
        listtri(j) = j
    end do

    kdtree = buildtree(listtri, 0)


    open(newunit=u,file='/home/lewis/phdshizz/meshMC/testkd.dat')

    call write_bbox(kdtree, u)
    close(u)


    stop
    write(u,*) kdtree%mindim
    write(u,*) kdtree%maxdim(1),kdtree%mindim(2),kdtree%mindim(3)
    write(u,*) kdtree%maxdim(1),kdtree%maxdim(2),kdtree%mindim(3)
    write(u,*) kdtree%mindim(1),kdtree%maxdim(2),kdtree%mindim(3)
    write(u,*) kdtree%mindim
    write(u,*) kdtree%mindim(1),kdtree%mindim(2),kdtree%maxdim(3)
    write(u,*) kdtree%maxdim(1),kdtree%mindim(2),kdtree%maxdim(3)
    write(u,*) kdtree%maxdim(1),kdtree%maxdim(2),kdtree%maxdim(3)
    write(u,*) kdtree%mindim(1),kdtree%maxdim(2),kdtree%maxdim(3)
    write(u,*) kdtree%mindim(1),kdtree%mindim(2),kdtree%maxdim(3)

    write(u,*) ''
    write(u,*) ''
    write(u,*) kdtree%ln%mindim
    write(u,*) kdtree%ln%maxdim(1),kdtree%ln%mindim(2),kdtree%ln%mindim(3)
    write(u,*) kdtree%ln%maxdim(1),kdtree%ln%maxdim(2),kdtree%ln%mindim(3)
    write(u,*) kdtree%ln%mindim(1),kdtree%ln%maxdim(2),kdtree%ln%mindim(3)
    write(u,*) kdtree%ln%mindim
    write(u,*) kdtree%ln%mindim(1),kdtree%ln%mindim(2),kdtree%ln%maxdim(3)
    write(u,*) kdtree%ln%maxdim(1),kdtree%ln%mindim(2),kdtree%ln%maxdim(3)
    write(u,*) kdtree%ln%maxdim(1),kdtree%ln%maxdim(2),kdtree%ln%maxdim(3)
    write(u,*) kdtree%ln%mindim(1),kdtree%ln%maxdim(2),kdtree%ln%maxdim(3)
    write(u,*) kdtree%ln%mindim(1),kdtree%ln%mindim(2),kdtree%ln%maxdim(3)

    write(u,*) ''
    write(u,*) ''
    write(u,*) kdtree%ln%ln%mindim
    write(u,*) kdtree%ln%ln%maxdim(1),kdtree%ln%ln%mindim(2),kdtree%ln%ln%mindim(3)
    write(u,*) kdtree%ln%ln%maxdim(1),kdtree%ln%ln%maxdim(2),kdtree%ln%ln%mindim(3)
    write(u,*) kdtree%ln%ln%mindim(1),kdtree%ln%ln%maxdim(2),kdtree%ln%ln%mindim(3)
    write(u,*) kdtree%ln%ln%mindim
    write(u,*) kdtree%ln%ln%mindim(1),kdtree%ln%ln%mindim(2),kdtree%ln%ln%maxdim(3)
    write(u,*) kdtree%ln%ln%maxdim(1),kdtree%ln%ln%mindim(2),kdtree%ln%ln%maxdim(3)
    write(u,*) kdtree%ln%ln%maxdim(1),kdtree%ln%ln%maxdim(2),kdtree%ln%ln%maxdim(3)
    write(u,*) kdtree%ln%ln%mindim(1),kdtree%ln%ln%maxdim(2),kdtree%ln%ln%maxdim(3)
    write(u,*) kdtree%ln%ln%mindim(1),kdtree%ln%ln%mindim(2),kdtree%ln%ln%maxdim(3)

    write(u,*) ''
    write(u,*) ''
    write(u,*) kdtree%ln%rn%mindim
    write(u,*) kdtree%ln%rn%maxdim(1),kdtree%ln%rn%mindim(2),kdtree%ln%rn%mindim(3)
    write(u,*) kdtree%ln%rn%maxdim(1),kdtree%ln%rn%maxdim(2),kdtree%ln%rn%mindim(3)
    write(u,*) kdtree%ln%rn%mindim(1),kdtree%ln%rn%maxdim(2),kdtree%ln%rn%mindim(3)
    write(u,*) kdtree%ln%rn%mindim
    write(u,*) kdtree%ln%rn%mindim(1),kdtree%ln%rn%mindim(2),kdtree%ln%rn%maxdim(3)
    write(u,*) kdtree%ln%rn%maxdim(1),kdtree%ln%rn%mindim(2),kdtree%ln%rn%maxdim(3)
    write(u,*) kdtree%ln%rn%maxdim(1),kdtree%ln%rn%maxdim(2),kdtree%ln%rn%maxdim(3)
    write(u,*) kdtree%ln%rn%mindim(1),kdtree%ln%rn%maxdim(2),kdtree%ln%rn%maxdim(3)
    write(u,*) kdtree%ln%rn%mindim(1),kdtree%ln%rn%mindim(2),kdtree%ln%rn%maxdim(3)





    write(u,*) ''
    write(u,*) ''
    write(u,*) kdtree%rn%mindim
    write(u,*) kdtree%rn%maxdim(1),kdtree%rn%mindim(2),kdtree%rn%mindim(3)
    write(u,*) kdtree%rn%maxdim(1),kdtree%rn%maxdim(2),kdtree%rn%mindim(3)
    write(u,*) kdtree%rn%mindim(1),kdtree%rn%maxdim(2),kdtree%rn%mindim(3)
    write(u,*) kdtree%rn%mindim
    write(u,*) kdtree%rn%mindim(1),kdtree%rn%mindim(2),kdtree%rn%maxdim(3)
    write(u,*) kdtree%rn%maxdim(1),kdtree%rn%mindim(2),kdtree%rn%maxdim(3)
    write(u,*) kdtree%rn%maxdim(1),kdtree%rn%maxdim(2),kdtree%rn%maxdim(3)
    write(u,*) kdtree%rn%mindim(1),kdtree%rn%maxdim(2),kdtree%rn%maxdim(3)
    write(u,*) kdtree%rn%mindim(1),kdtree%rn%mindim(2),kdtree%rn%maxdim(3)


close(u)
stop

    meshname = trim(meshname(:len_trim(meshname)-4))

    call gridset(id)

    iseed = -45345443+id

    iseed = -abs(iseed)  ! Random number seed must be negative for ran2

    call init_opt4

    if(id == 0)then
        print*, ''      
        print*,'# of photons to run',nphotons*int(numproc)
    end if

    !***** Set up density grid *******************************************

    !***** Set small distance for use in optical depth integration routines 
    !***** for roundoff effects when crossing cell walls
    delta = 1.e-8*(2.*zmax/nzg)
    !   deltay = 1.e-5*(2.*ymax/nyg)          !!!!!!!!1! impliment!!!!!!!!!!!!!!!!!!!!!!!!!
    !   deltaz = 1.e-5*(2.*zmax/nzg)
    nscatt=0

    call cpu_time(start)
    call cpu_time(start2)

    !loop over photons 
    call mpi_barrier(MPI_COMM_WORLD, error)
    print*,'Photons now running on core: ',colour(id, green)
    do j=1,nphotons

        call init_opt4

        if(j == 1000 .and. id == 0)then
            call cpu_time(finish2)
            print*,' '
            tmp = (finish2-start2)/1000.*real(nphotons)
            if(tmp >= 60.)then
                tmp = tmp / 60.
                if(tmp > 60)then
                    tmp = tmp / 60.
                    print*,str(tmp),' hrs'
                else
                    print*,str(tmp),' mins'
                end if
            else
                print*,str(tmp),' s'
            end if
            print*,' '
        end if

        tflag=.FALSE.

        !***** Release photon from point source *******************************
        call sourceph(xcell,ycell,zcell,iseed)

        !****** Find scattering location
        call tauint1(xcell,ycell,zcell,tflag,iseed,delta)

        !******** Photon scatters in grid until it exits (tflag=TRUE) 
        do while(tflag.eqv..FALSE.)
            ran = ran2(iseed)
            if(ran < albedo)then!interacts with tissue
                  call stokes(iseed)
                  nscatt = nscatt + 1        
               else
                  tflag=.true.
                  exit
            end if

            !************ Find next scattering location
            call tauint1(xcell,ycell,zcell,tflag,iseed,delta)

        end do
        if(mod(j,1000) == 0)then
            print *, colour(j, blue, bold),' scattered photons completed on core: ',colour(id, str(30+mod(id,7)), bold)
        end if
    end do      ! end loop over nph photons

    call mpi_barrier(mpi_comm_world, error)

    deallocate(xface,yface,zface)

    call MPI_REDUCE(jmeanin,jmeaninGLOBAL,nxg*nyg*nzg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
    call MPI_BARRIER(MPI_COMM_WORLD, error)

    call MPI_REDUCE(jmeanout,jmeanoutGLOBAL,nxg*nyg*nzg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
    call MPI_BARRIER(MPI_COMM_WORLD, error)

    call cpu_time(finish)
    if(finish-start.ge.60.)then
    print*,floor((finish-start)/60.)+mod(finish-start,60.)/100.
    else
    print*, 'time taken ~ ',colour(floor(finish-start/60.),red, bold),'s'
    end if
    if(id == 0)then
    print*,'Average # of scatters per photon:',nscatt/(nphotons*numproc)
    !write out files

    call writer(nphotons, numproc, meshname)
    print*,colour('write done',black,white_b,bold)
    end if
    call mpi_finalize(error)
end program mcpolar