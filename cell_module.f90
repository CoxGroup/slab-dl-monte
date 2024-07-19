! ***************************************************************************
! *   Copyright (C) 2007 - 2015 by J.A.Purton                               *
! *   john.purton[@]stfc.ac.uk                                              *
! *                                                                         *
! *   Contributors:                                                         *
! *   -------------                                                         *
! *   A.V.Brukhno - Free Energy Diff, planar pore & algorithm optimizations *
! *   andrey.brukhno[@]stfc.ac.uk abrukhno[@]gmail.com                      *
! *                                                                         *
! *   T.L.Underwood - Lattice/Phase-Switch MC method & optimizations        *
! *   t.l.Underwood[@]bath.ac.uk                                            *
! *                                                                         *
! *   J.Grant - Orthogonal PBC & Ewald branches, optimization of exclusions *
! *   r.j.grant[@]bath.ac.uk                                                *
! ***************************************************************************

!> @brief
!> - Simulation cell manipulations: configuration, geometry, PBC, some MC moves
!> @usage
!> - @stdusage
!> - @stdspecs
!> @using
!> - `control_type`
!> - `config_type`

!> @modulefor simulation cell: configuration, geometry, etc

module cell_module

    use kinds_f90
    use control_type
    use config_type
    use molecule_type

    use dcd_format_module

    implicit none

        !AB: this is only a temporary solution for 'job'
        !> local copy of global 'job' structure 
    type(control), save :: myjob

        !> array to store configurations read from input
    type(config), allocatable, dimension(:) :: cfgs

        !> number of atoms in ARCHIVE / HISTORY files
    integer, allocatable, dimension(:) :: natoms_hist, matoms_hist

        !> arrays for cavity bias grid
    integer, allocatable, dimension(:) :: ncavity_grid, ncavity_free, cavatt

        !> number of configuration replicas present
    integer, save :: nread, nconfigs

        !> flag for storing only a subset of molecules in trajectoy file(s)
    logical, save :: is_subset = .false.

    integer, save :: nfreed, ntaken

contains


!> @brief
!> - deallocates array(s) for configuration container
!> @using
!> - `kinds_f90`
!> - `slit_module, only : in_slit, deallocate_MFA_arrays`

!> deallocates array(s) for input configurations
subroutine dealloc_configs()

    use slit_module, only : in_slit, deallocate_MFA_arrays

    !implicit none

    integer :: fail(5)

    fail = 0

    if (allocated(cfgs)) deallocate(cfgs, stat = fail(1))
    if (allocated(natoms_hist)) deallocate(natoms_hist, stat = fail(2))
    if (allocated(matoms_hist)) deallocate(matoms_hist, stat = fail(3))

    ! deallocate cavity bias arrays
    if (allocated(ncavity_grid)) deallocate(ncavity_grid, stat = fail(4))
    if (allocated(ncavity_free)) deallocate(ncavity_free, stat = fail(5))

    if ( any(fail > 0) ) then

        call error (300)

    endif

    ! AB: deallocate arrays for slit if needed
    if( in_slit ) then

        call deallocate_MFA_arrays

    end if

end subroutine

!> @brief
!> - allocates array(s) for configurations to be read from input
!> @using
!> - `kinds_f90`

!> allocates array(s) for input configurations
subroutine alloc_configs()

    !use kinds_f90

    !implicit none

    integer :: fail(5), icfg

    fail = 0

    allocate(cfgs(nconfigs), stat = fail(1))
    allocate(natoms_hist(nconfigs), stat = fail(2))
    allocate(matoms_hist(nconfigs), stat = fail(3))

    ! allocate cavity bias bits here as small
    allocate(ncavity_grid(nconfigs), stat = fail(4))
    allocate(ncavity_free(nconfigs), stat = fail(5))

    if ( any(fail > 0) ) then

        call error (300)

    endif

end subroutine

!> @brief
!> - reads in '`nconfig`' configurations (cell, molecules and atoms)
!> @using
!> - `kinds_f90`
!> - `control_type`
!> - `constants_module, only : uout, ucfg`
!> - `config_module, only : inputconfig`
!> - `comms_mpi_module, only : master, is_parallel, idnode, gsync`
!> - `parse_module`

!> reads '`nconfig`' configurations from input CONFIG file
subroutine inputcells(job, nconfig)

    !use kinds_f90
    use control_type
    use constants_module, only : uout, ucfg
    use config_module, only : inputconfig
    use comms_mpi_module, only : master, is_parallel, idnode, gsync, gsum_world
    use parallel_loop_module, only : idgrp, open_nodefiles
    use parse_module

    !AB: planar pore (2D slit) stuff
    use slit_module, only : in_slit, allocate_MFA_arrays

    !implicit none

        !> control structure for the job (various flags and switches, see `control_type.f90`)
    type(control), intent(inout) :: job

        !> number of configurations to read from input
    integer, intent(inout) :: nconfig

    integer :: fail, istat, icfg, cfgfmt, i !, nread
    character :: record*128, word*32, ch_num*10, fnum*3, fname*10

    logical :: safe
    
    nconfigs = nconfig

    if( nconfigs > 999 ) then 

        write(ch_num,'(i10)') nconfigs
        call cry(uout,'', &
                    "ERROR: The number of configs (cells) is too large: "//trim(adjustL(ch_num))//&
                   &" (>999) !!!",999)

    endif

    nread = ucfg
    istat = 0

    fname = 'CONFIG'
    !fname(1:6) = 'CONFIG'

    if (master) then

        if( job%repexch ) then

            write(ch_num,'(i10)') idgrp

            call int_2_char3(idgrp, fnum, safe)

            if( .not.safe ) call cry(uout,'', &
                    "ERROR: The group index (idgrp) is too large: "//trim(adjustL(ch_num))//&
                   &" (>999) !!!",999)

            !fname = 'CONFIG.'//trim(adjustL(ch_num))
            fname = 'CONFIG.'//fnum

            open(nread, file=trim(fname), status = 'old', iostat=istat)

            write(uout,"(/,a,i3,2a,/)")" Group master ",idgrp," reading file '"//trim(fname)//"'" !'CONFIG.'//fnum

            !if( istat /= 0 )  call cry(uout,'', &
            !        "ERROR: Could not read input file : '"//trim(fname)//&
            !       &"' !!!",999)

        else !if( nconfigs == 1 ) then

            open(nread, file=trim(fname), status = 'old', iostat=istat)

            !if( istat /= 0 )  call cry(uout,'', &
            !        "ERROR: Could not read input file : '"//trim(fname)//&
            !       &"' !!!",999)

        endif

    endif

    if( is_parallel ) call gsum_world(istat)

    if( istat /= 0 )  call cry(uout,'', &
                           "ERROR: Could not read input file : '"//trim(fname)//"' !!!",999)

    call alloc_configs

    natoms_hist = 0
    matoms_hist = 0

    is_subset = ( len_trim(adjustL(job%hist_molids)) > 0 )

    do  icfg = 1, nconfigs

        fnum = ""
        safe = .true.

        !AB: the following code is for opening CONFIG files in the case of Gibbs ensemble
        !AB: but currently (as before) the two configurations can be read from a single CONFIG
        goto 101

        if( nconfigs > 1 ) then

            istat = 0
            if( master ) then

                write(fnum,'(i3)') icfg

                fname = 'CONFIG-'//trim(adjustL(fnum))

                open(nread, file=trim(fname), status = 'old', iostat=istat)
                !open(nread, file='CONFIG-'//trim(fname), form='formatted', status = 'old', iostat=istat)

                !if( istat /= 0 )  call cry(uout,'', &
                !    "ERROR: Could not read input file : '"//trim(fname)//"' !!!",999)

            end if

            if( is_parallel ) call gsum_world(istat)

            if( istat /= 0 )  call cry(uout,'', &
                              "ERROR: Could not read input file : '"//trim(fname)//"' !!!",999)

        endif

101     continue

        if( master ) then
            write(uout,"(/,/,1x,50('-'))")
            write(uout,"(a,i10)")" configuration box No. ",icfg !,nread
            write(uout,"(1x,50('-'))")
        endif

        call inputconfig(job, icfg, cfgs(icfg), cfgfmt)

        call invert_latticevectors(icfg)

        ! cfgfmt defines format of the input cell: 
        ! 1/2/3/6  - Cartesian (cubic/orthorhombic/parallelepiped/slit), 
        ! cfgfmt<1 - fractional
        if (abs(cfgfmt) > 3 .and. abs(cfgfmt) /= 6 ) call error(209)

        if (cfgfmt == 0) call fractocart(icfg)

        call pbc_simulation_cell(icfg)

        !AB: initialise the molecule COM:s
        call set_mol_coms(icfg)

        call printbasis(icfg)

        do i = 1, cfgs(icfg)%num_mols

           !AB: if the molecule name is in the list of those to store the coords of...
           if( is_subset .and. &
               index( job%hist_molids, cfgs(icfg)%mols(i)%molname ) < 1 ) cycle

           natoms_hist(icfg) = natoms_hist(icfg)+cfgs(icfg)%mols(i)%natom
           matoms_hist(icfg) = natoms_hist(icfg)+cfgs(icfg)%mols(i)%mxatom

        end do

        cycle
        
        if( nconfigs > 1 .and. master )  close(nread)

    end do

    myjob = job

    !if( job%usecavitybias ) write(uout,*)"inputcells(): cavity distance = ",myjob%cavity_radius

    if( is_parallel ) call gsync

    if( nconfigs == 1 .and. master )  close(nread)

    return

    call error(200)

end subroutine

!> @brief
!> - initialises arrays and parameters for Ewald summation scheme
!> @using
!> - `coul_module`
!> - `control_type`
!> - `latticevectors_module`
!> - `species_module, only : number_of_molecules`

!> prepares Ewald summation scheme
subroutine setup_ewald(ib, job)

    use coul_module
    use control_type
    use latticevectors_module
    use species_module, only : number_of_molecules

    !implicit none

        !> configuration identifier (index in '`cfgs`' array)
    integer, intent(in) :: ib

        !> control structure: various flags and switches, see `control_type.f90`
    type(control), intent(in) :: job

    integer :: i, j

    ! make sure recip lattice vectors are upto date

    call setrcpvec(cfgs(ib)%vec)

    call setewald(job, cfgs(ib)%vec, cfgs(ib)%cl)

    if( cfgs(ib)%newjob ) then 

        call printewald(ib, cfgs(ib)%cl)

        if( job%coultype(ib) == 1) then

          !call estimate_gvector_size(cfgs(ib)%vec, cfgs(ib)%cl, job%dielec, job%calcstress, job%distribgvec)
          call estimate_gvector_size(cfgs(ib)%vec, cfgs(ib)%cl, job%dielec, &
                                    (job%type_vol_move /= NPT_OFF), job%distribgvec)

          call allocate_ewald_arrays(.false., cfgs(ib)%cl, cfgs(ib)%maxno_of_atoms)

        end if

    end if

    if( job%coultype(ib) == 1 ) call findgvector(cfgs(ib)%vec, cfgs(ib)%cl, job%dielec, job%calcstress, job%distribgvec)

    if( cfgs(ib)%newjob ) then 

        call printnumgvec(cfgs(ib)%cl, job%repexch, job%distribgvec)

        call erfcgen(cfgs(ib)%cl, job)

        cfgs(ib)%newjob = .false.

    end if

end subroutine setup_ewald

!AB: <<-- front-end wrappers for input/output of trajectories in DL_MONTE and DL_POLY-2/4 formats

!> reads in configuration from trajectory file
subroutine read_config(ib, typ, is_close, natoms, cfgeng, istep, istat)

    !use control_type
    use constants_module
    use comms_mpi_module, only : master!, idnode
    use parallel_loop_module, only : open_nodefiles

    !implicit none

    !control structure for the job (various flags and switches, see `control_type.f90`)
    !type(control), intent(in) :: job

        !> replica identifier
    integer, intent(in)    :: ib 

        !> trajectory type identifier
    integer, intent(in)    :: typ

        !> flag specifying to close the trajectory file after reading one snapshot (frame)
    logical, intent(in)    :: is_close

        !> number of atoms read in
    integer, intent(inout) :: natoms

        !> MC step to read in
    integer, intent(inout) :: istep 

        !> input error flag
    integer, intent(inout) :: istat

        !> total energy of the configuration
    real(kind=wp), intent(inout) :: cfgeng

    ! flag for the first read, i.e. newly open trajectory file
    logical, save :: is_new = .true.

    ! auxilary variables
    integer :: lb, i, j

    if( master ) then

        if( is_new ) then

            if( typ < DLP2 ) then  ! < 2

                call open_nodefiles('ARCHIVE', utraj, istat)

                !write(uout,*)'read_config(..):: Opened ARCHIVE file for reading...'

            else if( typ < DCD  ) then ! < 6

                call open_nodefiles('HISTORY', utraj, istat)

            else
            !AB: nothing to do here for DCD trajectory
                return
            end if

        else if( is_close ) then 

            is_new = .true.
            close(utraj)
            return

        end if

        if( typ == DLP4 ) then

            call read_config_dlpoly4(utraj, ib, is_subset, is_new, natoms, cfgeng, istep, istat)

            !if(istat/=0) write(uout,*)'Failed reading HISTORY file, istat & is_new: ',istat,is_new
            !flush(uout)

            if( is_new ) is_new = .false.

        else if( typ == DLP2 ) then

            call read_config_dlpoly2(utraj, ib, is_subset, is_new, natoms, cfgeng, istep, istat)

            !if(istat/=0) write(uout,*)'Failed reading HISTORY file, istat & is_new: ',istat,is_new
            !flush(uout)

            if( is_new ) is_new = .false.

        else if( typ == DLP ) then

            call read_config_dlpoly(utraj, ib, is_new, natoms, cfgeng, istep, istat)

            !if(istat/=0) write(uout,*)'Failed reading ARCHIVE file, istat & is_new: ',istat,is_new
            !flush(uout)

            if( is_new ) is_new = .false.

        else if( typ < DLP ) then

            call read_config_dlmonte(utraj, ib, natoms, cfgeng, istep, istat)

            !if(istat/=0) write(uout,*)'Failed reading ARCHIVE file, istat & is_new: ',istat,is_new
            !flush(uout)

            if( is_new ) is_new = .false.

        else 

            is_new = .true.
            close(utraj)
            return

        end if

        if( is_close ) then 
            close(utraj)
            is_new = .true.
        end if

    endif

    cfgeng = cfgeng * myjob%energyunit

    return

    !write(iout,*)'Could not read HISTORY file - failure within a frame)'
    !flush(uout)

end subroutine read_config

!> @brief
!> - writes out configuration data to archive file in either DL_POLY or DL_MONTE format
!> - if '`typ`' == 1 then dl_poly format is assumed (otherwise dl_monte format)
!> @using
!> - `species_module`
!> - `molecule_type`
!> - `constants_module, only : urevc`
!> - `comms_mpi_module, only : idnode`
!> - `atom_module, only : write_atom_pos`

!> writes out configuration to archive file
subroutine dump_config(unit, ib, typ, cfgeng, istep)

    use constants_module
    use comms_mpi_module, only : master, idnode

    !implicit none
 
        !> unit to output configuration data to (i.e. the archive file)
    integer, intent(in)    :: unit

        !> replica identifier
    integer, intent(in)    :: ib 

        !> trajectory type identifier
    integer, intent(in)    :: typ

        !> MC step to read in
    integer, intent(in) :: istep 

        !> total energy of the configuration
    real(kind=wp), intent(in) :: cfgeng

    ! auxilary variables
    integer :: lb, i, j, ierr

    ! flag for the first read, i.e. newly open trajectory file
    logical, save :: is_new = .true.

    !AB: the file(s) are open by not only the root (idnode == 0) but every workgroup master (if job%repexch)
    !AB: in this case it is clearer to simply call this routine by master(s) only

    !if (master) then

        if( typ == VTK ) then 

            call write_config_vtk(unit, ib, istep, ierr)
            
        else if( typ < DLP ) then ! < 0

            call write_config_dlmonte(unit, ib, cfgeng/myjob%energyunit, istep, ierr)

        else if( typ == DLP ) then ! == 1

            call write_config_dlpoly(unit, ib, cfgeng/myjob%energyunit, is_new, istep, ierr)

            if( is_new ) is_new = .false.

        else if( typ < DLP4 ) then ! == 2/3

            call write_config_dlpoly2(unit, ib, cfgeng/myjob%energyunit, is_subset, is_new, istep, ierr)

            if( is_new ) is_new = .false.

        else if( typ < DCD ) then ! == 4/5

            call write_config_dlpoly4(unit, ib, cfgeng/myjob%energyunit, is_subset, is_new, istep, ierr)

            if( is_new ) is_new = .false.

        else ! default

            call write_config_dlmonte(unit, ib, cfgeng/myjob%energyunit, istep, ierr)

        end if

    !endif

end subroutine dump_config

!> @brief
!> - writes out configuration data to REVCON file in either DL_POLY or DL_MONTE format
!> - if '`typ`' == 1 then dl_poly format is assumed (otherwise dl_monte format)
!> @using
!> - `species_module`
!> - `molecule_type`
!> - `constants_module, only : urevc`
!> - `comms_mpi_module, only : idnode`
!> - `atom_module, only : write_atom_pos`

!> writes out configuration to REVCON file
subroutine dump_revcon(typ, istep, cfgeng)

    use constants_module
    use comms_mpi_module, only : master, idnode
    use parallel_loop_module, only : open_nodefiles

    !implicit none

        !> trajectory type identifier
    integer, intent(in)    :: typ

        !> MC step to read in
    integer, intent(in) :: istep 

        !> array of configuration energies
    real(kind=wp), dimension(:), intent(in) :: cfgeng

    ! auxilary variables
    integer :: ib, lb, i, j, ierr

    if( master ) then

    !AB: the file(s) are open by not only the root (idnode == 0) but every workgroup master (if job%repexch)

        !AB: all boxes are written in the same REVCON file, just like input in CONFIG (should be split!!!)
        !open(urevc, file = "REVCON.000")     ! write final config
        call open_nodefiles('REVCON', urevc, ierr)

        do ib = 1, nconfigs

            if( typ == VTK ) then
                
                call write_config_dlmonte(urevc, ib, cfgeng(ib)/myjob%energyunit, istep, ierr)

                close(urevc)

                call open_nodefiles('REVCON-VTK', urevc, ierr)
                
                call write_config_vtk(urevc, ib, istep, ierr)
                
            else if( typ < DLP ) then ! < 0

                call write_config_dlmonte(urevc, ib, cfgeng(ib)/myjob%energyunit, istep, ierr)

            else if( typ == DLP ) then ! == 1

                call write_config_dlmonte(urevc, ib, cfgeng(ib)/myjob%energyunit, istep, ierr)

                close(urevc)

                call open_nodefiles('REVCON-DLP', urevc, ierr)

                call write_config_dlpoly(urevc, ib, cfgeng(ib)/myjob%energyunit, .true., istep, ierr)

            else if( typ < DLP4 ) then ! == 2/3

                call write_config_dlmonte(urevc, ib, cfgeng(ib)/myjob%energyunit, istep, ierr)

                close(urevc)

                call open_nodefiles('REVCON-DLP2', urevc, ierr)

                call write_config_dlpoly2(urevc, ib, cfgeng(ib)/myjob%energyunit, .false., .true., istep, ierr)

            else if( typ < DCD ) then ! == 4/5

                call write_config_dlmonte(urevc, ib, cfgeng(ib)/myjob%energyunit, istep, ierr)

                close(urevc)

                call open_nodefiles('REVCON-DLP4', urevc, ierr)

                call write_config_dlpoly4(urevc, ib, cfgeng(ib)/myjob%energyunit, .false., .true., istep, ierr)

            else ! default

                call write_config_dlmonte(urevc, ib, cfgeng(ib)/myjob%energyunit, istep, ierr)

            end if

        end do

        close(urevc)

    end if

end subroutine dump_revcon

!AB: -->> front-end wrappers for input/output of trajectories in DL_MONTE and DL_POLY-2/4 formats

!AB: <<-- dealing with trajectories in DL_POLY-2/4 formats

!AB: *** WARNING *** 
!AB: ALL TRAJECTORY FORMATS ASSUME ORTHORHOMBIC CELL, imcon = 2 (cubic for DL_MONTE, imcon = 1)
!AB: imcon is not even a parameter here but replaced by the constants (directly set in write's)

!AB: *** WARNING ***
!AB: DL_POLY-2 format is deficient for HARD-SPHERE SYSTEMS due to the restricted precision in coordinates

!> @brief
!> - reads in configuration data in DL_POLY-2 format (HISTORY style)
!> - the input file unit is specified by '`nread`' parameter
!> @using
!> - `kinds_f90`
!> - `species_module, only : set_species_type, element`
!> - `atom_module, only : write_atom_pos`

!> reads in configuration in DL_POLY-2 format
subroutine read_config_dlpoly2(nread, ib, do_subset, is_header, natoms, cfgeng, istep, istat)

    !use kinds_f90
    use species_module, only : set_species_type, element
    use atom_module,    only : write_atom_pos

    !implicit none

        !> file unit for output (the file is opened externally) 
    integer, intent(in)    :: nread

        !> configuration identifier (index in '`cfgs`' array)
    integer, intent(in)    :: ib

        !> flags for writing only a subset of molecules and the file header
    logical, intent(in)    :: do_subset, is_header

        !> total number of atoms present, MC step and output (error) status
    integer, intent(inout) :: natoms, istep, istat

        !> array of configuration energies
    real(kind=wp), intent(inout) :: cfgeng

    integer, save :: imcon  = 2
    !integer, save :: imcon0 = 2   !AB: it must become a parameter (eventually)

    integer, save :: itype  = 0    !AB: this defines the level of HISTORY for DL_POLY (0 for coordiantes only)
    !integer, save :: itype0 = 0    !AB: this defines the level of HISTORY for DL_POLY (0 for coordiantes only)

    ! auxilary variables
    integer   :: nframes,lb, i, j, k, ic
    character :: word*8 !,tag*1

    !AB: would need to update this if molecule is allowed to change in size (natoms)
    !natoms_hist(ib) = natoms_hist(ib)+cfgs(ib)%mols(i)%natom

    if( is_header ) then 

        read(nread,'(a)', err=33, iostat=istat) cfgs(ib)%cfg_title

        cfgs(ib)%cfg_title = trim(adjustL(cfgs(ib)%cfg_title))

        read(nread,*, err=33, iostat=istat) itype,imcon,natoms !,istep

        cfgs(ib)%number_of_atoms = natoms

        istep = 0
        return

    end if

    !AB: ignoring timestep 'dtime'

    read(nread,*, err=66, iostat=istat) word,istep,natoms,itype,imcon !,0.001

    cfgs(ib)%number_of_atoms = natoms

    !if( imcon /= imcon0 ) then ... !AB: check if 'imcon' changed suddenly

    read(nread,*, err=99, iostat=istat) (cfgs(ib)%vec%latvector(1,i), i = 1,3)
    read(nread,*, err=99, iostat=istat) (cfgs(ib)%vec%latvector(2,i), i = 1,3)
    read(nread,*, err=99, iostat=istat) (cfgs(ib)%vec%latvector(3,i), i = 1,3)

    k = 0
    do i = 1, cfgs(ib)%num_mols

        !AB: if the molecule name is in the list of those to store the coords of...
        !if( do_subset .and. index( myjob%hist_molids, cfgs(ib)%mols(i)%molname ) < 1 ) cycle

        do j = 1, cfgs(ib)%mols(i)%natom
           k = k+1

            !AB: it is important to read the correct FIELD first, and allocate the necessary arrays!!!

            lb = cfgs(ib)%mols(i)%atms(j)%atlabel

            read(nread,'(a8,i10,2f12.6)', err=99, iostat=istat) element(lb),k, &
                 cfgs(ib)%mols(i)%atms(j)%mass,cfgs(ib)%mols(i)%atms(j)%charge

            read(nread,*, err=99, iostat=istat) (cfgs(ib)%mols(i)%atms(j)%rpos(ic), ic = 1,3)

        enddo

    enddo

    return

33  istat = 33 !AB: flag the failure in the first record for frame

    return

66  istat = 66 !AB: flag the failure in the first record for frame

    return

99  istat = 99
    !write(iout,*)'Could not read HISTORY file - failure within a frame)'
    !flush(uout)

    return

end subroutine read_config_dlpoly2

!> @brief
!> - writes out configuration data in DL_POLY-2 format (HISTORY style)
!> - the output file unit is specified by '`nwrite`' parameter
!> @using
!> - `kinds_f90`
!> - `constants_module, only : uout, urevc`
!> - `species_module, only : set_species_type, element`
!> - `atom_module, only : write_atom_pos`

!> writes out configuration in DL_POLY-2 format
subroutine write_config_dlpoly2(nwrite, ib, cfgeng, do_subset, is_header, istep, istat)

    !use kinds_f90
    use constants_module, only : DZERO, uout, urevc
    use species_module,   only : set_species_type, element, number_of_molecules
    use atom_module,      only : write_atom_pos
    use molecule_module,  only : reform_molecule_ortho, mol_com, full_mol_com, shift_molecule

    !implicit none

        !> file unit for output (the file must be pre-open externally!) 
    integer, intent(in) :: nwrite

        !> configuration identifier (index in '`cfgs`' array)
    integer, intent(in) :: ib

        !> array of configuration energies
    real(kind=wp), intent(in) :: cfgeng

        !> flags for writing only a subset of molecules and the file header
    logical, intent(in) :: do_subset, is_header

        !> MC step (iteration)
    integer, intent(in) :: istep

        !> output (error) status
    integer, intent(inout) :: istat

    !integer, intent(in) :: imcon   !AB: it must become a parameter eventually
    integer, save :: imcon = 2
    integer, save :: itype = 0      !AB: this defines the level of HISTORY for DL_POLY (0 for coordiantes only)

    ! auxilary variables
    real(kind=wp) :: mcom(3)
    integer :: lb, i, j, k, nvirt,nvmol, natmv
    character :: tag*1, fltfmt*2, intfmt*2, crdfmt*9

    if( myjob%nfltprec == 0 ) then
        write(crdfmt,"('(3g20.10)')")
    else
        write(fltfmt,'(i2)') myjob%nfltprec
        write(intfmt,'(i2)') myjob%nfltprec+8
        write(crdfmt,'(a)')"(3e"//trim(intfmt)//"."//trim(fltfmt)//")"
    end if

    !AB: NOTE - most recent DL_POLY-4 versions read/write only 72 characters!
    if( is_header ) write(nwrite,'(a)', err=66, iostat=istat) cfgs(ib)%cfg_title

    !AB: pretending that timestep = 0.001 (ps) - to allow for DL_POLY to accept it!

    if( nwrite /= urevc .and. do_subset ) then

        !AB: in case it's GCMC reset the number of atoms for output
        natoms_hist(ib) = 0
        do i = 1, cfgs(ib)%num_mols

           !AB: if the molecule name is in the list of those to store the coords of...
           if( index( myjob%hist_molids, cfgs(ib)%mols(i)%molname ) < 1 ) cycle

           natoms_hist(ib) = natoms_hist(ib)+cfgs(ib)%mols(i)%natom
           !matoms_hist(ib) = natoms_hist(ib)+cfgs(ib)%mols(i)%mxatom

        end do

        if( is_header ) write(nwrite,'(2i9,i10)', err=66, iostat=istat) &
                              itype,imcon,natoms_hist(ib)

        write(nwrite,'(a8,4i10,2f15.7,e15.7)', err=99, iostat=istat) &
             "timestep",istep,natoms_hist(ib),itype,imcon,0.001, myjob%systemp, cfgeng 
    else

        natoms_hist(ib) = cfgs(ib)%number_of_atoms
        if( nwrite == urevc .and. myjob%is_gcmc ) natoms_hist(ib) = cfgs(ib)%maxno_of_atoms

        if( is_header ) write(nwrite,'(2i9,i10)', err=66, iostat=istat) &
                              itype,imcon,natoms_hist(ib)

        write(nwrite,'(a8,4i10,2f15.7,e15.7)', err=99, iostat=istat) &
             "timestep",istep,natoms_hist(ib),itype,imcon,0.001, myjob%systemp, cfgeng 
    end if

    ! cell matrix (higher precision format; still read by VMD plugin)
    write(nwrite,crdfmt(1:len(crdfmt)), err=99, iostat=istat) (cfgs(ib)%vec%latvector(1,i), i = 1,3)
    write(nwrite,crdfmt(1:len(crdfmt)), err=99, iostat=istat) (cfgs(ib)%vec%latvector(2,i), i = 1,3)
    write(nwrite,crdfmt(1:len(crdfmt)), err=99, iostat=istat) (cfgs(ib)%vec%latvector(3,i), i = 1,3)
    
    !write(nwrite,'(3e15.7)', err=99, iostat=istat) (cfgs(ib)%vec%latvector(1,i), i = 1,3)
    !write(nwrite,'(3e15.7)', err=99, iostat=istat) (cfgs(ib)%vec%latvector(2,i), i = 1,3)
    !write(nwrite,'(3e15.7)', err=99, iostat=istat) (cfgs(ib)%vec%latvector(3,i), i = 1,3)

    ! cell matrix (DL_POLY-2 format)
    !write(nwrite,'(3g12.6)', err=99, iostat=istat) (cfgs(ib)%vec%latvector(1,i), i = 1,3)
    !write(nwrite,'(3g12.6)', err=99, iostat=istat) (cfgs(ib)%vec%latvector(2,i), i = 1,3)
    !write(nwrite,'(3g12.6)', err=99, iostat=istat) (cfgs(ib)%vec%latvector(3,i), i = 1,3)

    if( cfgs(ib)%num_mols < 1 ) return

    k = 0
    do i = 1, cfgs(ib)%num_mols

        !AB: if the molecule name is in the list of those to store the coords of...
        if( nwrite /= urevc .and. do_subset .and. index( myjob%hist_molids, cfgs(ib)%mols(i)%molname ) < 1 ) cycle

        !AB: unwrap molecule (undo PBC) if it is rigid or has structure due to bonds
        if( cfgs(ib)%mols(i)%rigid_body .or. cfgs(ib)%mols(i)%blist%npairs > 0 ) then

            !AB: before unwrapping molecule store its atom positions and com
            call store_molecule_pos(ib,i)

            if( myjob%useorthogonal ) then

                !if( uniq_mol(ml)%blist%npairs > 0 ) &
                call reform_molecule_ortho(cfgs(ib)%mols(i),cfgs(ib)%vec%latvector)
                call mol_com(cfgs(ib)%mols(i))

            else

                call full_mol_com(cfgs(ib)%mols(i), cfgs(ib)%vec%invlat,cfgs(ib)%vec%latvector)

            end if

            mcom(:) = cfgs(ib)%mols(i)%rcom(:)

            call pbc_cart_vec_calc(ib,mcom)

            mcom(:) = mcom(:) - cfgs(ib)%mols(i)%rcom(:)

            call shift_molecule(cfgs(ib)%mols(i), mcom)

        end if

        do j = 1, cfgs(ib)%mols(i)%natom
           k = k+1

            lb = cfgs(ib)%mols(i)%atms(j)%atlabel
            tag = set_species_type(cfgs(ib)%mols(i)%atms(j)%atype)

            write(nwrite,'(a8,i10,2f12.6)', err=99, iostat=istat) &
                  element(lb),k,cfgs(ib)%mols(i)%atms(j)%mass,cfgs(ib)%mols(i)%atms(j)%charge

            ! coordinates (higher precision format; still read by VMD plugin)
            !write(nwrite,'(3g20.10)', err=99, iostat=istat) &
            !write(nwrite,'(3e15.7)', err=99, iostat=istat) &
            write(nwrite,crdfmt(1:len(crdfmt)), err=99, iostat=istat) &
                  cfgs(ib)%mols(i)%atms(j)%rpos(1:3)

            ! coordinates (DL_POLY-2 format)
            !write(nwrite,'(1p,3e12.4)', err=99, iostat=istat) &

        enddo

        !AB: restore molecule atom positions and com
        if( cfgs(ib)%mols(i)%rigid_body .or. cfgs(ib)%mols(i)%blist%npairs > 0 ) &
            call restore_molecule_pos(ib,i)

    enddo

    if( nwrite == urevc ) then

        if( k /= cfgs(ib)%number_of_atoms ) &
            call cry(uout,'', &
                    "ERROR: failed storing REVCON file (DL_POLY2 style) - "//&
                    "inconsistent number of molecules/atoms in the system !!!",999)

        nvirt = cfgs(ib)%maxno_of_atoms - cfgs(ib)%number_of_atoms

        if( myjob%is_gcmc .and. nvirt > 0 ) then !.and. nwrite == urevc ) then

            if( cfgs(ib)%mtypes(number_of_molecules)%is_field ) then
            !if( cfgs(ib)%mols(cfgs(ib)%num_mols)%is_field ) then

                if( cfgs(ib)%mols(cfgs(ib)%num_mols)%mxatom - cfgs(ib)%mols(cfgs(ib)%num_mols)%natom /= nvirt ) then
                    write(uout,*)
                    write(uout,*)"write_config_dlpoly2(..):: GCMC detected: N_virt = ", nvirt, &
                                 ", last molecule is 'atomic field': N_max - N_atm = ",& 
                             cfgs(ib)%mols(cfgs(ib)%num_mols)%mxatom," - ",cfgs(ib)%mols(cfgs(ib)%num_mols)%natom," = ", &
                             (cfgs(ib)%mols(cfgs(ib)%num_mols)%mxatom - cfgs(ib)%mols(cfgs(ib)%num_mols)%natom)
                    write(uout,*)

                    call cry(uout,'', &
                         "ERROR: cannot store REVCON file (DL_POLY2 style) - "//&
                         "not the last molecule in FIELD being used for GCMC !!!",999)
                endif

                natmv = cfgs(ib)%mols(cfgs(ib)%num_mols)%natom
                lb  = cfgs(ib)%mols(cfgs(ib)%num_mols)%atms(natmv)%atlabel
                tag = set_species_type(cfgs(ib)%mols(cfgs(ib)%num_mols)%atms(natmv)%atype)

                !do i = 1,nvmol

                    !AB: if the molecule name is in the list of those to store the coords of...
                    !if( do_subset .and. index( myjob%hist_molids, cfgs(ib)%mols(i)%molname ) < 1 ) cycle

                    do j = cfgs(ib)%mols(cfgs(ib)%num_mols)%natom+1,cfgs(ib)%mols(cfgs(ib)%num_mols)%mxatom
                       k = k+1

                        !lb = cfgs(ib)%mols(cfgs(ib)%num_mols)%atms(j)%atlabel
                        !tag = set_species_type(cfgs(ib)%mols(cfgs(ib)%num_mols)%atms(j)%atype)

                        write(nwrite,'(a8,i10,2f12.6)', err=99, iostat=istat) &
                              adjustL("V"//trim(element(lb))), k, 0.0_wp, 0.0_wp
                              !"V"//trim(element(lb)),k,cfgs(ib)%mols(cfgs(ib)%num_mols)%atms(natmv)%mass, &
                              !cfgs(ib)%mols(cfgs(ib)%num_mols)%atms(natmv)%charge

                        ! coordinates (higher precision format; still read by VMD plugin)
                        write(nwrite,'(1p,3e15.7)', err=99, iostat=istat) 0.0_wp, 0.0_wp, 0.0_wp !&
                              !cfgs(ib)%mols(cfgs(ib)%num_mols)%atms(j)%rpos(1:3)

                        ! coordinates (DL_POLY-2 format)
                        !write(nwrite,'(1p,3e12.4)', err=99, iostat=istat) &

                    enddo

                !enddo

            else 
            
                if( mod(nvirt,cfgs(ib)%mols(cfgs(ib)%num_mols)%natom) > DZERO ) then
                    write(uout,*)
                    write(uout,*)"write_config_dlpoly2(..):: GCMC detected for last molecule: "//&
                                 "mod[ N_virt / N_atm(last) ] = ",nvirt," / ",& 
                                 cfgs(ib)%mols(cfgs(ib)%num_mols)%natom," = ",&
                                 mod(nvirt,cfgs(ib)%mols(cfgs(ib)%num_mols)%natom)

                                 !real(nvirt) / real(cfgs(ib)%mols(cfgs(ib)%num_mols)%natom)

                    call cry(uout,'', &
                            "ERROR: cannot store REVCON file (DL_POLY2 style) - "//&
                            "not the last molecule in FIELD being used for GCMC !!!",999)
                endif

                nvmol = nvirt / cfgs(ib)%mols(cfgs(ib)%num_mols)%natom

                do i = 1,nvmol

                    !AB: if the molecule name is in the list of those to store the coords of...
                    !if( do_subset .and. index( myjob%hist_molids, cfgs(ib)%mols(i)%molname ) < 1 ) cycle

                    do j = 1, cfgs(ib)%mols(cfgs(ib)%num_mols)%natom
                       k = k+1

                        lb = cfgs(ib)%mols(cfgs(ib)%num_mols)%atms(j)%atlabel
                        tag = set_species_type(cfgs(ib)%mols(cfgs(ib)%num_mols)%atms(j)%atype)

                        write(nwrite,'(a8,i10,2f12.6)', err=99, iostat=istat) &
                              adjustL("V"//trim(element(lb))), k, cfgs(ib)%mols(cfgs(ib)%num_mols)%atms(j)%mass, &
                              cfgs(ib)%mols(cfgs(ib)%num_mols)%atms(j)%charge
                              !element(lb),k,cfgs(ib)%mols(cfgs(ib)%num_mols)%atms(j)%mass, &
                              !cfgs(ib)%mols(cfgs(ib)%num_mols)%atms(j)%charge

                        ! coordinates (higher precision format; still read by VMD plugin)
                        write(nwrite,'(1p,3e15.7)', err=99, iostat=istat) 0.0_wp, 0.0_wp, 0.0_wp !&
                              !cfgs(ib)%mols(cfgs(ib)%num_mols)%atms(j)%rpos(1:3)

                        ! coordinates (DL_POLY-2 format)
                        !write(nwrite,'(1p,3e12.4)', err=99, iostat=istat) &

                    enddo

                enddo

            endif

            if( k /= cfgs(ib)%maxno_of_atoms ) then
                write(uout,*)
                write(uout,*)"write_config_dlpoly2(..):: GCMC detected for last molecule: "//&
                             "N_virt, N_atm(last) =?= ",nvirt,", ", nvmol, " * ", cfgs(ib)%mols(cfgs(ib)%num_mols)%natom

                call cry(uout,'', &
                        "ERROR: failed storing REVCON file (DL_POLY2 style) - "//&
                        "inconsistent number of (virtual) molecules/atoms due to GCMC !!!",999)
            endif

        endif

    endif

    return

66  istat = 66 !AB: flag the failure in the first record for frame

    return

99  istat = 99
    !write(iout,*)'Could not write to HISTORY file - failure within a frame'
    !flush(uout)

end subroutine write_config_dlpoly2

!> @brief
!> - reads in configuration data in DL_POLY-4 format (HISTORY style)
!> - the input file unit is specified by '`nread`' parameter
!> @using
!> - `kinds_f90`
!> - `species_module, only : set_species_type, element`
!> - `atom_module, only : write_atom_pos`

!> reads in configuration in DL_POLY-4 format
subroutine read_config_dlpoly4(nread, ib, do_subset, is_header, natoms, cfgeng, istep, istat)

    !use kinds_f90
    use species_module, only : set_species_type, element
    use atom_module,    only : write_atom_pos

    !implicit none

        !> file unit for output (the file must be pre-open externally!) 
    integer, intent(in)    :: nread

        !> configuration identifier (index in '`cfgs`' array)
    integer, intent(in)    :: ib

        !> flags for writing only a subset of molecules and the file header
    logical, intent(in)    :: do_subset, is_header

        !> total number of atoms present, MC step and output (error) status
    integer, intent(inout) :: natoms, istep, istat

        !> array of configuration energies
    real(kind=wp), intent(inout) :: cfgeng

    !integer, intent(in) :: imcon   !AB: it must become a parameter (eventually)
    integer, save :: imcon = 2
    integer, save :: itype = 0      !AB: this defines the level of HISTORY for DL_POLY (0 for coordiantes only)

    ! auxilary variables
    integer :: lb, i, j, k, ic
    character :: word*8, tag*1

    real(kind=wp) :: dummy3(3), dummy

    !AB: NOTE - most recent DL_POLY-4 versions read/write only 72 characters!

    if( is_header ) then 

        read(nread,'(a)', err=33, iostat=istat) cfgs(ib)%cfg_title

        cfgs(ib)%cfg_title = trim(adjustL(cfgs(ib)%cfg_title))

        read(nread,'(3i10,i21)', err=33, iostat=istat) itype,imcon,natoms,istep!,0

        cfgs(ib)%number_of_atoms = natoms

        istep = 0
        return

    end if

    !AB: ignoring timestep 'dtime' and total time at the end of the header

    read(nread,*, err=66, iostat=istat) word,istep,natoms,itype,imcon !,0.001!,float(istep)*0.001 

    cfgs(ib)%number_of_atoms = natoms

    read(nread,'(3f20.10)', err=99, iostat=istat) (cfgs(ib)%vec%latvector(1,i), i = 1,3)
    read(nread,'(3f20.10)', err=99, iostat=istat) (cfgs(ib)%vec%latvector(2,i), i = 1,3)
    read(nread,'(3f20.10)', err=99, iostat=istat) (cfgs(ib)%vec%latvector(3,i), i = 1,3)

    k = 0
    do i = 1, cfgs(ib)%num_mols

        do j = 1, cfgs(ib)%mols(i)%natom
           k = k+1

            !AB: important to read the correct FIELD first, and allocate the necessary arrays!!!

            lb = cfgs(ib)%mols(i)%atms(j)%atlabel

            !AB: displacements at the end are set to zero (to be added?)
            read(nread,'(a8,i10,3f12.6)', err=99, iostat=istat) element(lb),k, &
                        cfgs(ib)%mols(i)%atms(j)%mass,cfgs(ib)%mols(i)%atms(j)%charge,dummy !,0.0 

            read(nread,'(3g20.10)', err=99, iostat=istat) (cfgs(ib)%mols(i)%atms(j)%rpos(ic), ic = 1,3)

        enddo

    enddo

    return

33  istat = 33 !AB: flag the failure in the first record for frame

    return

66  istat = 66 !AB: flag the failure in the first record for frame

    return

99  istat = 99
    !write(iout,*)'Could not read HISTORY file - failure within a frame)'
    !flush(uout)

    return

end subroutine read_config_dlpoly4

!> @brief
!> - writes out configuration data in DL_POLY-4 format (HISTORY style)
!> - the output file unit is specified by '`nwrite`' parameter
!> @using
!> - `kinds_f90`
!> - `species_module, only : set_species_type, element`
!> - `atom_module, only : write_atom_pos`

!> writes out configuration in DL_POLY-4 format
subroutine write_config_dlpoly4(nwrite, ib, cfgeng, do_subset, is_header, istep, istat)

    !use kinds_f90
    use species_module,  only : set_species_type, element
    use atom_module,     only : write_atom_pos
    use molecule_module, only : reform_molecule_ortho, mol_com, full_mol_com, shift_molecule

    !implicit none

        !> file unit for output (the file must be pre-open externally!) 
    integer, intent(in) :: nwrite

        !> configuration identifier (index in '`cfgs`' array)
    integer, intent(in) :: ib

        !> array of configuration energies
    real(kind=wp), intent(in) :: cfgeng

        !> flags for writing only a subset of molecules and the file header
    logical, intent(in) :: do_subset, is_header

        !> MC step (iteration)
    integer, intent(in) :: istep

        !> output (error) status
    integer, intent(inout) :: istat

    !integer, intent(in) :: imcon   !AB: it must become a parameter (eventually)
    integer, save :: imcon = 2
    integer, save :: itype = 0      !AB: this defines the level of HISTORY for DL_POLY (0 for coordiantes only)

    ! auxilary variables
    real(kind=wp) :: mcom(3)
    integer :: lb, i, j, k
    !character :: tag*1

    !AB: would need to update this if molecule is allowed to change in size (natoms)
    !natoms_hist(ib) = natoms_hist(ib)+cfgs(b)%mols(i)%natom

    !AB: NOTE - most recent DL_POLY-4 versions read/write only 72 characters!
    if( is_header ) write(nwrite,'(a80)', err=66, iostat=istat) cfgs(ib)%cfg_title           !record 1 (a80)

    !AB: pretending that timestep = 0.001 (ps) - to allow for DL_POLY to accept it!

    if( do_subset ) then

        !AB: in case it's GCMC reset the number of atoms for output
        natoms_hist(ib) = 0
        do i = 1, cfgs(ib)%num_mols

           !AB: if the molecule name is in the list of those to store the coords of...
           if( index( myjob%hist_molids, cfgs(ib)%mols(i)%molname ) < 1 ) cycle

           natoms_hist(ib) = natoms_hist(ib)+cfgs(ib)%mols(i)%natom
           !matoms_hist(ib) = natoms_hist(ib)+cfgs(ib)%mols(i)%mxatom

        end do

        !AB: in the header istep must be replaced by nstep - the total number of frames saved!
        if( is_header ) write(nwrite,'(3i10,2i21)', err=66, iostat=istat) &
                              itype,imcon,natoms_hist(ib),istep,0

        write(nwrite,'(a8,2i10,2i2,2f20.6)', err=99, iostat=istat) &
             "timestep",istep,natoms_hist(ib),itype,imcon,0.001,float(istep)*0.001

    else
        if( is_header ) write(nwrite,'(3i10,2i21)', err=66, iostat=istat) &
                              itype,imcon,cfgs(ib)%number_of_atoms,istep,0

        write(nwrite,'(a8,2i10,2i2,2f20.6)', err=99, iostat=istat) &
              "timestep",istep,cfgs(ib)%number_of_atoms,itype,imcon,0.001,float(istep)*0.001 
    end if

    write(nwrite,'(3f20.10,a12)', err=99, iostat=istat) (cfgs(ib)%vec%latvector(1,i), i = 1,3) &
                ,repeat(' ',12)

    write(nwrite,'(3f20.10,a12)', err=99, iostat=istat) (cfgs(ib)%vec%latvector(2,i), i = 1,3) &
                ,repeat(' ',12)

    write(nwrite,'(3f20.10,a12)', err=99, iostat=istat) (cfgs(ib)%vec%latvector(3,i), i = 1,3) &
                ,repeat(' ',12)

    k = 0
    do i = 1, cfgs(ib)%num_mols

        !AB: if the molecule name is in the list of those to store the coords of...
        if( do_subset .and. index( myjob%hist_molids, cfgs(ib)%mols(i)%molname ) < 1 ) cycle

        !AB: unwrap molecule (undo PBC) if it is rigid or has structure due to bonds
        if( cfgs(ib)%mols(i)%rigid_body .or. cfgs(ib)%mols(i)%blist%npairs > 0 ) then

            !AB: before unwrapping molecule store its atom positions and com
            call store_molecule_pos(ib,i)

            if( myjob%useorthogonal ) then

                !if( uniq_mol(ml)%blist%npairs > 0 ) &
                call reform_molecule_ortho(cfgs(ib)%mols(i),cfgs(ib)%vec%latvector)
                call mol_com(cfgs(ib)%mols(i))

            else

                call full_mol_com(cfgs(ib)%mols(i), cfgs(ib)%vec%invlat,cfgs(ib)%vec%latvector)

            end if

            mcom(:) = cfgs(ib)%mols(i)%rcom(:)

            call pbc_cart_vec_calc(ib,mcom)

            mcom(:) = mcom(:) - cfgs(ib)%mols(i)%rcom(:)

            call shift_molecule(cfgs(ib)%mols(i), mcom)

        end if

        do j = 1, cfgs(ib)%mols(i)%natom
           k = k+1

            lb = cfgs(ib)%mols(i)%atms(j)%atlabel

            !AB: displacements at the end are set to zero (to be added?)
            write(nwrite,'(a8,i10,3f12.6,a18)', err=99, iostat=istat) &
                  element(lb),k,cfgs(ib)%mols(i)%atms(j)%mass,cfgs(ib)%mols(i)%atms(j)%charge,0.0 &
                              ,repeat(' ',18)

            write(nwrite,'(3g20.10,a12)', err=99, iostat=istat) &
                  cfgs(ib)%mols(i)%atms(j)%rpos(1:3) &
                        ,repeat(' ',12) !*10.0, i = 1,3) 
    
        enddo
        
        if( cfgs(ib)%mols(i)%rigid_body .or. cfgs(ib)%mols(i)%blist%npairs > 0 ) &
            call restore_molecule_pos(ib,i)

    enddo

    return

66  istat = 66 !AB: flag the failure in the first record for frame

    return

99  istat = 99
    !write(iout,*)'Could not write to HISTORY file - failure within a frame'
    !flush(uout)

end subroutine write_config_dlpoly4

!AB: -->> dealing with trajectories in DL_POLY-2/4 formats

!AB: <<-- dealing with ARCHIVE trajectories in old DL_MONTE/DL_POLY styles

!> @brief
!> - reads in configuration data in ''DL_POLY-style'' format (HISTORY)
!> - the input file unit is specified by '`nread`' parameter
!> @using
!> - `kinds_f90`
!> - `species_module, only : set_species_type, element`
!> - `atom_module, only : write_atom_pos`

!> reads in configuration in ''DL_POLY-style'' format
subroutine read_config_dlpoly(nread, ib, is_header, natoms, cfgeng, istep, istat)

!AB: *** WARNING *** 
!AB: "DL_POLY style" is a hybrid of DL_POLY-2 & DL_MONTE formats

    !use kinds_f90
    use species_module, only : set_species_type, element
    use atom_module,    only : write_atom_pos

    !implicit none

        !> file unit for output (the file must be pre-open externally!) 
    integer, intent(in) :: nread

        !> configuration identifier (index in '`cfgs`' array)
    integer, intent(in) :: ib

        !> flag for writing the file header
    logical, intent(in) :: is_header

        !> array of configuration energies
    real(kind=wp), intent(inout) :: cfgeng

        !> total number of atoms present, MC step and output (error) status
    integer, intent(inout) :: natoms, istep, istat

    ! auxilary variables
    integer, save :: iter  = 0
    integer, save :: itype = 0   !AB: this defines the level of HISTORY for DL_POLY (0 for coordiantes only)
    integer, save :: imcon = 2   !AB: it must become a parameter (eventually)

    integer, save :: matoms = 0
    integer, save :: nmols  = 0

    character :: word*8, tag*1

    real(4)  :: dtime
    integer  :: lb, i, j, k, katm

    if( is_header ) then

        ! title
        read(nread,'(a80)', err=66, iostat=istat) cfgs(ib)%cfg_title

        ! header record
        read(nread,"(2i9,3i10,10i10)", err=66, iostat=istat) itype, imcon, natoms, matoms, nmols, &
                                                             cfgs(ib)%mxmol_type

        !AB: reset atom/molecule numbers to the actually found
        cfgs(ib)%number_of_atoms = natoms
        cfgs(ib)%maxno_of_atoms  = matoms
        cfgs(ib)%num_mols        = nmols

        istep = 0
        return

    end if

    ! timestep record
    read(nread,*, err=99, iostat=istat) word, istep, natoms, itype, imcon, dtime, myjob%systemp, cfgeng 

    ! cell matrix
    read(nread,*, err=99, iostat=istat) (cfgs(ib)%vec%latvector(1,i), i = 1,3)
    read(nread,*, err=99, iostat=istat) (cfgs(ib)%vec%latvector(2,i), i = 1,3)
    read(nread,*, err=99, iostat=istat) (cfgs(ib)%vec%latvector(3,i), i = 1,3)

    k = 0
    do i = 1, cfgs(ib)%num_mols

        do j = 1, cfgs(ib)%mols(i)%natom
           k = k+1

            lb = cfgs(ib)%mols(i)%atms(j)%atlabel
            tag = set_species_type(cfgs(ib)%mols(i)%atms(j)%atype)

            read(nread,'(a8,i10,2f12.6,2x,a1)', err=99, iostat=istat) element(lb), katm, &
                 cfgs(ib)%mols(i)%atms(j)%mass, cfgs(ib)%mols(i)%atms(j)%charge, tag 

            read(nread,*, err=99, iostat=istat) cfgs(ib)%mols(i)%atms(j)%rpos(1:3)
 
        enddo

    enddo

    return

66  istat = 66 !AB: flag the failure in the first record for frame

    return

99  istat = 99
    !write(iout,*)'Could not read ARCHIVE file - failure within a frame'
    !flush(uout)

    return

end subroutine read_config_dlpoly

!> @brief
!> - writes out configuration data in DL_POLY format
!> - the output file unit is specified by '`nwrite`' parameter
!> @using
!> - `kinds_f90`
!> - `species_module, only : set_species_type, element`
!> - `atom_module, only : write_atom_pos`

!> writes out configuration in ''DL_POLY-style'' format
subroutine write_config_dlpoly(nwrite, ib, cfgeng, is_header, istep, istat)!, do_subset)

!AB: *** WARNING *** 
!AB: "DL_POLY style" is a hybrid of DL_POLY-2 & DL_MONTE formats

    !use kinds_f90
    use species_module,  only : set_species_type, element
    use atom_module,     only : write_atom_pos
    use molecule_module, only : reform_molecule_ortho, mol_com, full_mol_com, shift_molecule

    !implicit none

        !> file unit for output (the file must be pre-open externally!) 
    integer, intent(in) :: nwrite

        !> configuration identifier (index in '`cfgs`' array)
    integer, intent(in) :: ib

        !> array of configuration energies
    real(kind=wp), intent(in) :: cfgeng

        !> flag for writing the file header
    logical, intent(in) :: is_header

        !> current MC step  (iteration)
    integer, intent(in) :: istep

        !> output (error) status
    integer, intent(inout) :: istat

    integer, save :: itype = 0   !AB: this defines the level of HISTORY for DL_POLY (0 for coordiantes only)
    integer, save :: imcon = 2   !AB: it must become a parameter (eventually)

    ! auxilary variables
    real(kind=wp) :: mcom(3)
    integer :: lb, i, j, k
    character :: string*32,tag*1
    character :: fltfmt*2, intfmt*2, crdfmt*9

    if( myjob%nfltprec == 0 ) then
        write(crdfmt,"('(3f12.6)')")
    else
        write(fltfmt,'(i2)') myjob%nfltprec
        write(intfmt,'(i2)') myjob%nfltprec+8
        write(crdfmt,'(a)')"(3e"//trim(intfmt)//"."//trim(fltfmt)//")"
    end if

    !AB: NO SUBSET output due to DL_MONTE format incomatibility otherwise

    !AB: DL_POLY-2 format + DL_MONTE extras:
    if( is_header ) then

        ! title
        write(nwrite,'(a80)', err=66, iostat=istat) cfgs(ib)%cfg_title

        ! header record
        write(nwrite,"(2i9,3i10,10i10)", err=66, iostat=istat) &
              itype,imcon,cfgs(ib)%number_of_atoms, cfgs(ib)%maxno_of_atoms, &
                          cfgs(ib)%num_mols, cfgs(ib)%mxmol_type(:)

    end if

    ! timestep record
    write(nwrite,"(a8,4i10,2f15.7,e15.7)", err=66, iostat=istat) &
                 "timestep",istep,cfgs(ib)%number_of_atoms,itype,imcon,0.001, myjob%systemp, cfgeng

    ! cell matrix (higher precision format taking same space)
    write(nwrite,crdfmt(1:len(crdfmt)), err=99, iostat=istat) (cfgs(ib)%vec%latvector(1,i), i = 1,3)
    write(nwrite,crdfmt(1:len(crdfmt)), err=99, iostat=istat) (cfgs(ib)%vec%latvector(2,i), i = 1,3)
    write(nwrite,crdfmt(1:len(crdfmt)), err=99, iostat=istat) (cfgs(ib)%vec%latvector(3,i), i = 1,3)

    !write(nwrite,'(3f12.6)', err=99, iostat=istat) (cfgs(ib)%vec%latvector(1,i), i = 1,3)
    !write(nwrite,'(3f12.6)', err=99, iostat=istat) (cfgs(ib)%vec%latvector(2,i), i = 1,3)
    !write(nwrite,'(3f12.6)', err=99, iostat=istat) (cfgs(ib)%vec%latvector(3,i), i = 1,3)

    !AB: atom positions (DL_POLY-2 format)
    k = 0
    do i = 1, cfgs(ib)%num_mols

        write(string,"(2x,a8,2i8)") &
              cfgs(ib)%mols(i)%molname, cfgs(ib)%mols(i)%natom, cfgs(ib)%mols(i)%mxatom

        !AB: unwrap molecule (undo PBC) if it is rigid or has structure due to bonds
        if( cfgs(ib)%mols(i)%rigid_body .or. cfgs(ib)%mols(i)%blist%npairs > 0 ) then

            !AB: before unwrapping molecule store its atom positions and com
            call store_molecule_pos(ib,i)

            if( myjob%useorthogonal ) then

                !if( uniq_mol(ml)%blist%npairs > 0 ) &
                call reform_molecule_ortho(cfgs(ib)%mols(i),cfgs(ib)%vec%latvector)
                call mol_com(cfgs(ib)%mols(i))

            else

                call full_mol_com(cfgs(ib)%mols(i), cfgs(ib)%vec%invlat,cfgs(ib)%vec%latvector)

            end if

            mcom(:) = cfgs(ib)%mols(i)%rcom(:)

            call pbc_cart_vec_calc(ib,mcom)

            mcom(:) = mcom(:) - cfgs(ib)%mols(i)%rcom(:)

            call shift_molecule(cfgs(ib)%mols(i), mcom)

        end if

        do j = 1, cfgs(ib)%mols(i)%natom
           k = k+1

            lb = cfgs(ib)%mols(i)%atms(j)%atlabel
            tag = set_species_type(cfgs(ib)%mols(i)%atms(j)%atype)

            write(nwrite,'(a8,i10,2f12.6,2x,a)', err=99, iostat=istat) element(lb), k, &
                  cfgs(ib)%mols(i)%atms(j)%mass, cfgs(ib)%mols(i)%atms(j)%charge, tag//trim(string)

            !write(nwrite,'(1p,3f12.6)', err=99, iostat=istat) &
            write(nwrite,crdfmt(1:len(crdfmt)), err=99, iostat=istat) &
                  cfgs(ib)%mols(i)%atms(j)%rpos(1:3)

            string=""
        end do

        if( cfgs(ib)%mols(i)%rigid_body .or. cfgs(ib)%mols(i)%blist%npairs > 0 ) &
            call restore_molecule_pos(ib,i)

    end do

    return

66  istat = 66 !AB: flag the failure in the first record for frame

    return

99  istat = 99
    !write(iout,*)'Could not write to ARCHIVE file - failure within a frame'
    !flush(uout)

end subroutine write_config_dlpoly


!> @brief
!> - reads in configuration data in DL_MONTE format (ARCHIVE)
!> - the input file unit is specified by '`nread`' parameter
!> @using
!> - `kinds_f90`
!> - `parse_module`
!> - `species_module, only : set_species_type, element`
!> - `atom_module, only : write_atom_pos`

!> reads in configuration in DL_MONTE format
subroutine read_config_dlmonte(nread, ib, natoms, cfgeng, istep, istat)

    !use kinds_f90
    use parse_module
    use species_module, only : set_species_type, element
    use atom_module,    only : write_atom_pos

    !implicit none

        !> file unit for output (the file must be pre-open externally!) 
    integer, intent(in) :: nread

        !> configuration identifier (index in '`cfgs`' array)
    integer, intent(in) :: ib

        !> total number of atoms present, MC step and output (error) status
    integer, intent(inout) :: natoms, istep, istat

        !> array of configuration energies
    real(kind=wp), intent(inout) :: cfgeng

    ! auxilary variables
    integer, save :: iter  = 0
    integer, save :: itype = 0   !AB: this defines the level of HISTORY for DL_POLY (0 for coordiantes only)
    integer, save :: imcon = 2   !AB: it must become a parameter (eventually)
    integer, save :: matoms = 0

    integer   :: lb, i, j, k

    character :: record*128, word*17, word2*32, tag*1

    logical   :: safe

    iter  = iter+1
    istep = 0

    call get_line(safe, nread, record)
    if (.not.safe) go to 66

    cfgs(ib)%cfg_title = record(1:80)

    call get_line(safe, nread, record)
    if (.not.safe) go to 66

    call get_word(record, word)
    if( len_trim(word) < 1) goto 66
    itype = nint(word_2_real(word))

    call get_word(record, word)
    if( len_trim(word) < 1) goto 66
    imcon = nint(word_2_real(word))

    !AB: new format convention (if found):
    call get_word(record, word)
    if( len_trim(word) < 1) goto 66

    natoms = nint(word_2_real(word))

    call get_word(record, word)
    if( trim(adjustL(word)) == "MAXATM" .or. trim(adjustL(word)) == "EXTRAS" ) then

    call get_word(record, word)
    if( len_trim(word) < 1) goto 66
        matoms = nint(word_2_real(word))

        call get_word(record, word)
        if( len_trim(word) < 1) goto 66
        istep = nint(word_2_real(word))

        call get_word(record, word)
        if( len_trim(word) < 1) goto 66
        myjob%systemp = word_2_real(word)

        call get_word(record, word)
        if( len_trim(word) < 1) goto 66
        cfgeng = word_2_real(word)

    end if

    read(nread,*, err=99, iostat=istat) (cfgs(ib)%vec%latvector(1,i), i = 1,3)
    read(nread,*, err=99, iostat=istat) (cfgs(ib)%vec%latvector(2,i), i = 1,3)
    read(nread,*, err=99, iostat=istat) (cfgs(ib)%vec%latvector(3,i), i = 1,3)

    read(nread,*, err=99, iostat=istat) word(1:7), cfgs(ib)%num_mols, cfgs(ib)%mxmol_type

    k = 0
    natoms = 0
    do i = 1, cfgs(ib)%num_mols

        read(nread,*, err=99, iostat=istat) word, & !tag, word, &
              cfgs(ib)%mols(i)%molname, cfgs(ib)%mols(i)%natom, cfgs(ib)%mols(i)%mxatom

        natoms = natoms+cfgs(ib)%mols(i)%natom

        do j = 1, cfgs(ib)%mols(i)%natom
           k = k+1

            lb = cfgs(ib)%mols(i)%atms(j)%atlabel

            read(nread, "(1x,a8,1x,a1)", err=99, iostat=istat) element(lb), tag

            read(nread,*, err=99, iostat=istat) cfgs(ib)%mols(i)%atms(j)%rpos(1:3)

        end do

    end do

    !if( istep == 0 ) istep = iter

    return

66  istat = 66 !AB: flag the failure in the first record for frame

    return

99  istat = 99
    !write(iout,*)'Could not read ARCHIVE file - failure within a frame'
    !flush(uout)

    return

end subroutine read_config_dlmonte


!> @brief
!> - writes out configuration data in a Legacy VTK format
!> - the output file unit is specified by '`nwrite`' parameter
subroutine write_config_vtk(nwrite, ib, istep, istat)

        !> file unit for output (the file must be pre-open externally!) 
    integer, intent(in) :: nwrite

        !> configuration identifier (index in '`cfgs`' array)
    integer, intent(in) :: ib

        !> move number
    integer, intent(in) :: istep

        !> output (error) status
    integer, intent(inout) :: istat

    integer :: i, j

    write(nwrite,'(a)') "# vtk DataFile Version 2.0" 
    write(nwrite,'(2a,i10)') trim(cfgs(ib)%cfg_title),";  timestep = ",istep
    write(nwrite,'(a)') "ASCII"
    write(nwrite,'(a)') "DATASET POLYDATA"

    write(nwrite,'(a,i7,a)') "POINTS", cfgs(ib)%number_of_atoms, " FLOAT" 
    do j = 1, cfgs(ib)%num_mols
        do i = 1, cfgs(ib)%mols(j)%natom
            write(nwrite,'(3es15.7)') cfgs(ib)%mols(j)%atms(i)%rpos
        end do
    end do

    write(nwrite,'(a,i7)') "POINT_DATA",cfgs(ib)%number_of_atoms
    write(nwrite,'(a)') "VECTORS spin FLOAT"
    do j = 1, cfgs(ib)%num_mols
        do i = 1, cfgs(ib)%mols(j)%natom
            write(nwrite,'(3es15.7)') cfgs(ib)%mols(j)%atms(i)%spin
        end do
    end do
    
end subroutine write_config_vtk


!> @brief
!> - writes out configuration data in DL_MONTE format
!> - the output file unit is specified by '`nwrite`' parameter
!> @using
!> - `kinds_f90`
!> - `species_module, only : set_species_type, element`
!> - `atom_module, only : write_atom_pos`

!> writes out configuration in DL_MONTE format
subroutine write_config_dlmonte(nwrite, ib, cfgeng, istep, istat)

    !use kinds_f90
    use species_module,  only : set_species_type, element
    use atom_module,     only : write_atom_pos
    use molecule_module, only : reform_molecule_ortho, mol_com, full_mol_com, shift_molecule
    use constants_module, only : ATM_TYPE_SPIN

    !implicit none

        !> file unit for output (the file must be pre-open externally!) 
    integer, intent(in) :: nwrite

        !> configuration identifier (index in '`cfgs`' array)
    integer, intent(in) :: ib, istep

        !> array of configuration energies
    real(kind=wp), intent(in) :: cfgeng

        !> output (error) status
    integer, intent(inout) :: istat

    integer, save :: itype = 0   !AB: this defines the level of HISTORY for DL_POLY (0 for coordiantes only)
    integer, save :: imcon = 2   !AB: it must become a parameter (eventually)

    ! auxilary variables
    real(kind=wp) :: mcom(3)
    integer :: lb, i, j, k
    character :: tag*1, fltfmt*2, intfmt*2, crdfmt*9

    if( myjob%nfltprec == 0 ) then
        write(crdfmt,"('(3f15.7)')")
    else
        write(fltfmt,'(i2)') myjob%nfltprec
        write(intfmt,'(i2)') myjob%nfltprec+8
        write(crdfmt,'(a)')"(3e"//trim(intfmt)//"."//trim(fltfmt)//")"
    end if

    ! title
    write(nwrite,'(a80)', err=66, iostat=istat) cfgs(ib)%cfg_title

    ! header record
    write(nwrite,"(2i3,i10,2x,'EXTRAS',1x,2i10,f15.7,e15.7)", err=66, iostat=istat) &
              itype,imcon,cfgs(ib)%number_of_atoms, cfgs(ib)%maxno_of_atoms, istep, myjob%systemp, cfgeng

    ! cell matrix
    write(nwrite,"("//crdfmt(1:len(crdfmt))//")", err=99, iostat=istat) (cfgs(ib)%vec%latvector(1,i), i = 1,3)
    write(nwrite,"("//crdfmt(1:len(crdfmt))//")", err=99, iostat=istat) (cfgs(ib)%vec%latvector(2,i), i = 1,3)
    write(nwrite,"("//crdfmt(1:len(crdfmt))//")", err=99, iostat=istat) (cfgs(ib)%vec%latvector(3,i), i = 1,3)

    !write(nwrite,'(3f15.7)', err=99, iostat=istat) (cfgs(ib)%vec%latvector(1,i), i = 1,3)
    !write(nwrite,'(3f15.7)', err=99, iostat=istat) (cfgs(ib)%vec%latvector(2,i), i = 1,3)
    !write(nwrite,'(3f15.7)', err=99, iostat=istat) (cfgs(ib)%vec%latvector(3,i), i = 1,3)

    ! all molecules record
    write(nwrite,"('NUMMOL',1x,i10,10i10)", err=99, iostat=istat) &
          cfgs(ib)%num_mols, cfgs(ib)%mxmol_type(:)

    k = 0
    do i = 1, cfgs(ib)%num_mols

        ! individual molecule record
        write(nwrite,"('MOLECULE',1x,a8,2i10)", err=99, iostat=istat) &
              cfgs(ib)%mols(i)%molname, cfgs(ib)%mols(i)%natom, cfgs(ib)%mols(i)%mxatom

        !AB: unwrap molecule (undo PBC) if it is rigid or has structure due to bonds
        if( cfgs(ib)%mols(i)%rigid_body .or. cfgs(ib)%mols(i)%blist%npairs > 0 ) then

            !AB: before unwrapping molecule store its atom positions and com
            call store_molecule_pos(ib,i)

            if( myjob%useorthogonal ) then

                !if( uniq_mol(ml)%blist%npairs > 0 ) &
                call reform_molecule_ortho(cfgs(ib)%mols(i),cfgs(ib)%vec%latvector)
                call mol_com(cfgs(ib)%mols(i))

            else

                call full_mol_com(cfgs(ib)%mols(i), cfgs(ib)%vec%invlat,cfgs(ib)%vec%latvector)

            end if

            mcom(:) = cfgs(ib)%mols(i)%rcom(:)

            call pbc_cart_vec_calc(ib,mcom)

            mcom(:) = mcom(:) - cfgs(ib)%mols(i)%rcom(:)

            call shift_molecule(cfgs(ib)%mols(i), mcom)

        end if

        do j = 1, cfgs(ib)%mols(i)%natom
           k = k+1

            lb = cfgs(ib)%mols(i)%atms(j)%atlabel
            tag = set_species_type(cfgs(ib)%mols(i)%atms(j)%atype)

            write(nwrite, "(1x,a8,1x,a1)", err=99, iostat=istat) element(lb), tag

            write(nwrite,"("//crdfmt(1:len(crdfmt))//",2x,i4)") &
                  cfgs(ib)%mols(i)%atms(j)%rpos(1:3), cfgs(ib)%mols(i)%atms(j)%site

            !write(nwrite, '(3f15.7,2x,i4)') cfgs(ib)%mols(i)%atms(j)%rpos(1:3), cfgs(ib)%mols(i)%atms(j)%site
            !call write_atom_pos(cfgs(ib)%mols(i)%atms(j), nwrite)

            !TU: Write the spins for type 'spin' atoms
            if( cfgs(ib)%mols(i)%atms(j)%atype == ATM_TYPE_SPIN ) then

                write(nwrite,"("//crdfmt(1:len(crdfmt))//")") cfgs(ib)%mols(i)%atms(j)%spin(1:3)
                
            end if
            
        enddo

        if( cfgs(ib)%mols(i)%rigid_body .or. cfgs(ib)%mols(i)%blist%npairs > 0 ) &
            call restore_molecule_pos(ib,i)

    enddo

    return

66  istat = 66 !AB: flag the failure in the first record for frame

    return

99  istat = 99
    !write(iout,*)'Could not write to ARCHIVE file - failure within a frame'
    !flush(uout)

    return

end subroutine write_config_dlmonte

!AB: -->> dealing with ARCHIVE trajectories in old DL_MONTE/DL_POLY styles

!AB: <<-- dealing with trajectories in DCD format

!> @brief
!> - reads in configuration data in DCD format
!> @using
!> - `kinds_f90`
!> - `constants_module`, only : uout
!> - `comms_mpi_module`, only : master
!> - `parse_module`,     only : int_2_char3

!> reads in configuration from DCD file
subroutine read_config_dcd(fname, iodcd, iform, ib, icell, matoms, mframes, iframe, ifstep, iclose)

    use kinds_f90
    use constants_module, only : uout
    use parse_module,     only : int_2_char3
    use comms_mpi_module, only : master !, idnode
    use parallel_loop_module, only : idgrp

    implicit none

    ! name of DCD trajectory file 
    character*(*), intent(in) :: fname

    !> file unit for output (the file must be pre-open externally!) 
    integer, intent(in) :: iodcd

    !> flag for single/double precision (32/64 bit)
    integer, intent(in) :: iform

    !> configuration identifier (index in '`cfgs`' array)
    integer, intent(in) :: ib

    !> flag for storing (or not) the cell info
    integer, intent(inout) :: icell

    !> max number of atoms (expected)
    integer, intent(inout) :: matoms

    !> max number of frames (expected)
    integer, intent(inout) :: mframes

    !> index of current frame
    integer, intent(inout) :: iframe

    !> frame stepping (number of frames to skip)
    integer, intent(inout) :: ifstep

    !> flag for final closing of DCD file
    integer, intent(inout) :: iclose

    ! <- internal variables

    real(4), save, allocatable :: x_dcd4(:), y_dcd4(:), z_dcd4(:)
    real(8), save, allocatable :: x_dcd8(:), y_dcd8(:), z_dcd8(:)

    real(kind=wp), save :: cell_pdb(6)
    real(kind=wp)       :: cell_min, cell_var, rscale
    real(4)             :: dtime

    ! title of DCD trajectory file 
    character*(80) :: dcd_title(2)

    ! name of DCD trajectory file 
    character*(80) :: dcd_fname

    ! four symbols PDB header (normally 'CORD')
    character*(4)  :: dcd_header != 'CORD'
    character*(8)  :: dcd_number
    character*(64) :: dcd_string

    logical, save :: is_new  = .true.
    logical, save :: is_open = .false.
    integer, save :: natoms = 0
    integer, save :: nframe = 0
    integer, save :: jframe = 0
    integer, save :: nwarns = 0
    integer, save :: iopen  = 1

    integer :: fail(3)
    integer :: matomst,natomst

    integer :: ijframe, i, j, k !, lb
    character :: tag*1

    character(3)  :: char_num1, char_num2

    logical :: safe = .true.

    ! -> internal variables

    if( .not.master ) return
    
    write(dcd_number,'(i8)') idgrp

    call int_2_char3(idgrp, char_num2, safe)

    if( .not.safe ) call cry(uout,'', &
                    "ERROR: The replica index is too large: "//trim(adjustL(dcd_number))//" (>999) !!!",999)

    if( myjob%nconfigs > 1 ) then

        write(dcd_number,'(i8)') ib

        if( ib > 999 ) call cry(uout,'', &
                        "ERROR: The cell/box index is too large: "//trim(adjustL(dcd_number))//" (>999) !!!",999)

        write(dcd_fname,'(a)') trim(adjustL( trim(adjustL(fname))//'-'//trim(adjustL(dcd_number))//&
                         &'.'//trim(adjustL(char_num2)) ))

    else

        write(dcd_fname,'(a)') trim(adjustL( trim(adjustL(fname))//'.'//trim(adjustL(char_num2)) ))

    end if

    dcd_title(:) = ""
    dcd_header   = ""

    if( is_new ) then 

        matoms  = 0
        mframes = 0
        iframe  = 0
        ifstep  = 0
        icell   = 0

        is_new  = .false.

        call dcd_open_for_reading(iodcd, dcd_fname, dcd_title, dcd_header, iclose, &
                                  natoms, nframe, iframe, ifstep, icell, dtime)

        jframe  = iframe
        matoms  = natoms
        mframes = nframe

        is_open = (iclose == 0)

        if( iclose /= 0 ) then 
        !AB: DCD file has been closed and only HEADER was supposed to be read

            call cry(uout,'', &
                     "ERROR: the input DCD file '"//trim(dcd_fname)// &
                     "' unexpectedly closed after reading the HEADER !!!",999)

            return

        end if

    fail = 0
    if( iform == 0 ) then

        allocate(x_dcd4(cfgs(ib)%maxno_of_atoms), stat = fail(1))
        allocate(y_dcd4(cfgs(ib)%maxno_of_atoms), stat = fail(2))
        allocate(z_dcd4(cfgs(ib)%maxno_of_atoms), stat = fail(3))

        if( any(fail > 0) ) &
            call cry(uout,'', &
                     "ERROR: failed allocating memory for DCD arrays '"//&
                    &"'- full stop !!!",999)
    else

        allocate(x_dcd8(cfgs(ib)%maxno_of_atoms), stat = fail(1))
        allocate(y_dcd8(cfgs(ib)%maxno_of_atoms), stat = fail(2))
        allocate(z_dcd8(cfgs(ib)%maxno_of_atoms), stat = fail(3))

        if( any(fail > 0) ) &
            call cry(uout,'', &
                     "ERROR: failed allocating memory for DCD arrays '"//&
                    &"'- full stop !!!",999)
    end if

        iopen = 0

        return

    end if

    matoms  = natoms
    mframes = nframe

    if( iclose /= 0 .or. iopen /= 0 ) &
        call cry(uout,'', &
                 "ERROR: the input DCD file '"//trim(dcd_fname)// &
                 "' is not expected to be re-open before reading next FRAME !!!",999)

    cell_pdb = 0.0_wp

    if( iform == 0 ) then

        x_dcd4 = 0.0
        y_dcd4 = 0.0
        z_dcd4 = 0.0

        call dcd_read_next_frame4(iodcd, dcd_fname, iopen, iclose, &
                                 natoms, icell, x_dcd4, y_dcd4, z_dcd4, cell_pdb)
    else

        x_dcd8 = 0.0
        y_dcd8 = 0.0
        z_dcd8 = 0.0

        call dcd_read_next_frame8(iodcd, dcd_fname, iopen, iclose, &
                                 natoms, icell, x_dcd8, y_dcd8, z_dcd8, cell_pdb)
    end if

    !AB: failed to read next frame
    if( iclose > 90 ) goto 101

    iframe = iframe+ifstep

    !AB: if no cell specs present => the cell/box PBC are taken from CONFIG (be wary!)
    !AB: otherwise the cell/box specs are read for each frame from DCD file

    !AB: VMD can't handle non-orthorhombic (non-orthogonal) unit cells! 
    !AB: thus, PDB angles must be set to 90 degress for VMD to show PBC etc
    if( icell == 1 ) then

        cfgs(ib)%vec%latvector(1,1) = cell_pdb(1) ! pbc_cell: A & angle1 (cubic)
        cfgs(ib)%vec%latvector(2,2) = cell_pdb(1) ! pbc_cell: A & angle1 (cubic)
        cfgs(ib)%vec%latvector(3,3) = cell_pdb(1) ! pbc_cell: A & angle1 (cubic)

        !cell_pdb(2) = 90.0
        !cell_pdb(4) = 90.0
        !cell_pdb(5) = 90.0

    else if( icell == 2 ) then

        cfgs(ib)%vec%latvector(1,1) = cell_pdb(1) ! pbc_cell: A & angle1 (orthorhombic)
        cfgs(ib)%vec%latvector(2,2) = cell_pdb(3) ! pbc_cell: B & angle2 (orthorhombic)
        cfgs(ib)%vec%latvector(3,3) = cell_pdb(6) ! pbc_cell: C & angle3 (orthorhombic)

        !cell_pdb(2) = 90.0
        !cell_pdb(4) = 90.0
        !cell_pdb(5) = 90.0

    !else if( icell > 2 ) then
    else if( icell /= 0 ) then ! non-orthorhombic cells, angles /= 90, so let's stick to the cell matrix
        !AB: some sources refer to unit cells in the lower diagonal matrix representation (instead of PDB format: a,b,c + angles)

        cfgs(ib)%vec%latvector(1,1) = cell_pdb(1) ! pbc_cell: lower diagonal
        cfgs(ib)%vec%latvector(1,2) = cell_pdb(2) ! pbc_cell: lower diagonal
        cfgs(ib)%vec%latvector(2,1) = cell_pdb(2) ! pbc_cell: lower diagonal
        cfgs(ib)%vec%latvector(2,2) = cell_pdb(3) ! pbc_cell: lower diagonal
        cfgs(ib)%vec%latvector(1,3) = cell_pdb(4) ! pbc_cell: lower diagonal
        cfgs(ib)%vec%latvector(3,1) = cell_pdb(4) ! pbc_cell: lower diagonal
        cfgs(ib)%vec%latvector(2,3) = cell_pdb(5) ! pbc_cell: lower diagonal
        cfgs(ib)%vec%latvector(3,2) = cell_pdb(5) ! pbc_cell: lower diagonal
        cfgs(ib)%vec%latvector(3,3) = cell_pdb(6) ! pbc_cell: lower diagonal

    !else if( icell > ??? ) then ! crystallographic data (CRYST1 in PDB)

        !AB: TO DO - transform DL_MONTE cell into PDB/PBC representation !!!
        !AB: need to introduce a conversion routine: cell angles -> cell matrix
        !AB: currently only cubic / orthorhombic cells can be dealt with

    end if

    rscale = 0.1 !AB: nm -> Angstrom (coordinate/distance re-scaling needed?)

    k = 0
    if( iform == 0 ) then

        do  i = 1, cfgs(ib)%num_mols
            do  j = 1, cfgs(ib)%mols(i)%natom
                k = k + 1

                if( k > natoms ) then
                    iclose = 9999
                    exit
                end if

                cfgs(ib)%mols(i)%atms(j)%rpos(1) = x_dcd4(k) !*rscale
                cfgs(ib)%mols(i)%atms(j)%rpos(2) = y_dcd4(k) !*rscale
                cfgs(ib)%mols(i)%atms(j)%rpos(3) = z_dcd4(k) !*rscale
            end do
        end do

    else

        do  i = 1, cfgs(ib)%num_mols
            do  j = 1, cfgs(ib)%mols(i)%natom
                k = k + 1

                if( k > natoms ) then
                    iclose = 9999
                    exit
                end if

                cfgs(ib)%mols(i)%atms(j)%rpos(1) = x_dcd8(k) !*rscale
                cfgs(ib)%mols(i)%atms(j)%rpos(2) = y_dcd8(k) !*rscale
                cfgs(ib)%mols(i)%atms(j)%rpos(3) = z_dcd8(k) !*rscale
            end do
        end do

    end if

    if( iclose == 0 ) return

101 fail = 0
    if( iform == 0 ) then
        deallocate(x_dcd4, stat = fail(1))
        deallocate(y_dcd4, stat = fail(2))
        deallocate(z_dcd4, stat = fail(3))

        if( any(fail > 0) ) &
            call cry(uout,'', &
                     "ERROR: failed deallocating memory after DCD arrays '"//&
                    &"'- full stop !!!",999)
    else
        deallocate(x_dcd8, stat = fail(1))
        deallocate(y_dcd8, stat = fail(2))
        deallocate(z_dcd8, stat = fail(3))

        if( any(fail > 0) ) &
            call cry(uout,'', &
                     "ERROR: failed deallocating memory after DCD arrays '"//&
                    &"'- full stop !!!",999)
    endif

end subroutine read_config_dcd

!AB: <- writing trajectory in DCD format

!TU: This does not work for Gibbs because two trajectories, one corresponding to each box,
!TU: and the 'saved' variables in the below function cannot cope with openning two different
!TU: units

!> @brief
!> - writes out configuration data in DCD format
!> @using
!> - `kinds_f90`
!> - `constants_module`, only : uout
!> - `comms_mpi_module`, only : master
!> - `parse_module`,     only : int_2_char3

!> writes out configuration in DCD format
subroutine write_config_dcd(fname, iodcd, iform, ib, icell, matoms, mframes, iframe, ifstep, iclose)

    use kinds_f90
    use constants_module, only : DZERO, uout
    use molecule_module,  only : reform_molecule_ortho, mol_com, full_mol_com, shift_molecule
    use parse_module,     only : int_2_char3
    use comms_mpi_module, only : master !, idnode
    use parallel_loop_module, only : idgrp

    !use species_module, only : set_species_type, element
    !use atom_module, only : write_atom_pos

    implicit none

    ! name of DCD trajectory file 
    character*(*), intent(in) :: fname

    !> file unit for output (the file must be pre-open externally!) 
    integer, intent(in) :: iodcd

    !> flag for single/double precision (32/64 bit)
    integer, intent(in) :: iform

    !> configuration identifier (index in '`cfgs`' array)
    integer, intent(in) :: ib

    !> flag for storing (or not) the cell info
    integer, intent(in) :: icell

    !> max number of atoms (expected)
    integer, intent(in) :: matoms

    !> max number of frames (expected)
    integer, intent(in) :: mframes

    !> index of current frame
    integer, intent(in) :: iframe

    !> frame stepping (number of frames to skip)
    integer, intent(in) :: ifstep

    !> flag for final closing of DCD file
    integer, intent(inout) :: iclose

    ! <- internal variables

    real(4), save, allocatable :: x_dcd4(:), y_dcd4(:), z_dcd4(:)
    real(8), save, allocatable :: x_dcd8(:), y_dcd8(:), z_dcd8(:)

    real(kind=wp), save :: cell_pdb(6)
    real(kind=wp)       :: cell_min, cell_var, rscale, mcom(3)

    ! title of DCD trajectory file 
    character*(80) :: dcd_title(2)

    ! name of DCD trajectory file 
    character*(80) :: dcd_fname

    ! four symbols PDB header (normally 'CORD')
    character*(4)  :: dcd_header = 'CORD'
    character*(8)  :: dcd_number
    character*(64) :: dcd_string

    logical, save :: is_new  = .true.
    logical, save :: is_open = .false.
    integer, save :: natoms = 0
    integer, save :: nframe = 0
    integer, save :: jframe = 0
    integer, save :: nwarn0 = 0
    integer, save :: nwarn1 = 0
    integer, save :: nwarn2 = 0
    integer, save :: nwarn3 = 0
    integer, save :: nwarn4 = 0

    integer :: fail(3)
    integer :: iopen !, iclose
    integer :: matomst,natomst
    integer :: ijframe, i, j, k !, lb

    character(3) :: char_num1, char_num2

    logical :: safe = .true.

    ! -> internal variables

    if( cfgs(ib)%num_mols < 1 ) return

    matomst = cfgs(ib)%maxno_of_atoms
    natomst = cfgs(ib)%number_of_atoms

    if( matomst < 1 .or. natomst < 1 .or. matomst /= matoms ) then

        write(dcd_string,'(3i10)') matoms,matomst,natomst

        call cry(uout,'', &
                 "ERROR: inconsistent number of atoms for DCD file: dcd_atoms, max_atoms, sys_atoms = '"//&
                &trim(dcd_string)//"' - check the setup !!!",999)

    end if

    if( .not.master ) return
    
    write(dcd_number,'(i8)') idgrp

    call int_2_char3(idgrp, char_num2, safe)

    if( .not.safe ) call cry(uout,'', &
                    "ERROR: The replica index is too large: "//trim(adjustL(dcd_number))//" (>999) !!!",999)

    if( myjob%nconfigs > 1 ) then

        write(dcd_number,'(i8)') ib

        if( ib > 999 ) call cry(uout,'', &
                        "ERROR: The cell/box index is too large: "//trim(adjustL(dcd_number))//" (>999) !!!",999)

        write(dcd_fname,'(a)') trim(adjustL( trim(adjustL(fname))//'-'//trim(adjustL(dcd_number))//&
                         &'.'//trim(adjustL(char_num2)) ))

    else

        write(dcd_fname,'(a)') trim(adjustL( trim(adjustL(fname))//'.'//trim(adjustL(char_num2)) ))

    end if

    iopen = 0

    write(dcd_title(1),'(a)') trim(adjustL(cfgs(ib)%cfg_title))

    cell_pdb = 0.0_wp

    if( is_new ) then 

        is_new  = .false.

        call dcd_open_for_writing(iodcd, dcd_fname, dcd_title, dcd_header, iclose, &
                                  matoms, mframes, iframe, ifstep, icell, 1.0)

        !AB: need to introduce a conversion routine: cell matrix -> cell angles
        !AB: currently only cubic / orthorhombic cells can be dealt with

        !AB: VMD can't handle non-orthorhombic (non-orthogonal) unit cells! 
        !AB: thus, PDB angles must be set to 90 degress for VMD to show PBC etc

        if( icell == 1 ) then

            cell_pdb(1) = cfgs(ib)%vec%latvector(1,1) ! pbc_cell: cubic
            cell_pdb(3) = cfgs(ib)%vec%latvector(1,1) ! pbc_cell: cubic
            cell_pdb(6) = cfgs(ib)%vec%latvector(1,1) ! pbc_cell: cubic
            cell_pdb(2) = 90.0
            cell_pdb(4) = 90.0
            cell_pdb(5) = 90.0

        else if( icell == 2 ) then

            cell_pdb(1) = cfgs(ib)%vec%latvector(1,1) ! pbc_cell: orthorhombic
            cell_pdb(3) = cfgs(ib)%vec%latvector(2,2) ! pbc_cell: orthorhombic
            cell_pdb(6) = cfgs(ib)%vec%latvector(3,3) ! pbc_cell: orthorhombic
            cell_pdb(2) = 90.0
            cell_pdb(4) = 90.0
            cell_pdb(5) = 90.0

        !else if( icell > 2 ) then
        else if( icell /= 0 ) then ! non-orthorhombic cells, angles /= 90, so let's stick to the cell matrix
        !AB: some sources refer to unit cells in the lower diagonal matrix representation (instead of PDB format: a,b,c + angles)

            cell_pdb(1) = cfgs(ib)%vec%latvector(1,1) ! pbc_cell: lower diagonal
            cell_pdb(2) = cfgs(ib)%vec%latvector(2,1) ! pbc_cell: lower diagonal
            cell_pdb(3) = cfgs(ib)%vec%latvector(2,2) ! pbc_cell: lower diagonal
            cell_pdb(4) = cfgs(ib)%vec%latvector(3,1) ! pbc_cell: lower diagonal
            cell_pdb(5) = cfgs(ib)%vec%latvector(3,2) ! pbc_cell: lower diagonal
            cell_pdb(6) = cfgs(ib)%vec%latvector(3,3) ! pbc_cell: lower diagonal

        !else if( icell > ??? ) then ! crystallographic data (CRYST1 in PDB)

    end if

    is_open = (iclose == 0)
    if( is_open ) iopen = 1

    natoms = matoms
    jframe = iframe

    fail = 0

    if( iform == 0 ) then

        allocate(x_dcd4(matoms), stat = fail(1))
        allocate(y_dcd4(matoms), stat = fail(2))
        allocate(z_dcd4(matoms), stat = fail(3))

        if( any(fail > 0) ) &
            call cry(uout,'', &
                     "ERROR: failed allocating memory for DCD arrays '"//&
                    &"'- full stop !!!",999)

    else

        allocate(x_dcd8(matoms), stat = fail(1))
        allocate(y_dcd8(matoms), stat = fail(2))
        allocate(z_dcd8(matoms), stat = fail(3))

        if( any(fail > 0) ) &
            call cry(uout,'', &
                     "ERROR: failed allocating memory for DCD arrays '"//&
                    &"'- full stop !!!",999)

    end if

    else

        ijframe = iframe-jframe

        if( ijframe > 0 .and. ijframe /= ifstep ) then

            nwarn0 = nwarn0+1

            if( nwarn0 < 11 ) then 

                call cry(uout,'', &
                     "WARNING: unexpected frames are stored in DCD file '"//trim(dcd_fname)//&
                    &"' - check the frames step and/or sampling frequency !!!",0)

            else if( nwarn0 == 11 ) then 

                call cry(uout,'', &
                     "WARNING: unexpected frames are stored in DCD file '"//trim(dcd_fname)//&
                    &"' - keeping silent about more such events !!!",0)

            end if

        end if

        jframe = iframe

    end if

    nframe = nframe + 1

    if( mframes == nframe ) then

        iclose  = 1
        is_open = .false.

    else if( mframes < nframe .and. is_open ) then

        nwarn1 = nwarn1+1

        if( nwarn1 < 11 ) then 

            call cry(uout,'', &
                 "WARNING: # frames overflow for DCD file '"//trim(dcd_fname)// &
                 "' - closing the file now; no more frames expected !!!",0)

        else if( nwarn1 == 11 ) then 

            call cry(uout,'', &
                 "WARNING: # frames overflow for DCD file '"//trim(dcd_fname)// &
                 "' - keeping silent about more such events !!!",0)
        end if

        close(iodcd)

        iclose  = 1
        is_open = .false.

        return

    end if

    if( icell == 0 ) then

        cell_min = min( cfgs(ib)%vec%latvector(1,1), &
                        cfgs(ib)%vec%latvector(2,2), &
                        cfgs(ib)%vec%latvector(3,3) )

        cell_var = max( abs(cell_pdb(1) - cfgs(ib)%vec%latvector(1,1)), &
                        abs(cell_pdb(3) - cfgs(ib)%vec%latvector(2,2)), &
                        abs(cell_pdb(6) - cfgs(ib)%vec%latvector(3,3)) ) 

        if( cell_var/cell_min  > DZERO ) then

            nwarn2 = nwarn2+1

            if( nwarn2 < 11 ) then 

                call cry(uout,'', &
                     "WARNING: cell variation detected for DCD file '"//trim(dcd_fname)// &
                     "' - must be NPT; consider storing the cell specs !!!",0)

            else if( nwarn2 == 11 ) then 

                call cry(uout,'', &
                     "WARNING: cell variation detected for DCD file '"//trim(dcd_fname)// &
                     "' - keeping silent about more such events !!!",0)
            end if

        end if

    end if

    !AB: TO DO - transform DL_MONTE cell into PDB/PBC representation !!!
    !AB: this only works for cubic & orthorhombic cells/boxes

    if( icell == 1 ) then

        cell_pdb(1) = cfgs(ib)%vec%latvector(1,1) ! pbc_cell: cubic
        cell_pdb(3) = cfgs(ib)%vec%latvector(1,1) ! pbc_cell: cubic
        cell_pdb(6) = cfgs(ib)%vec%latvector(1,1) ! pbc_cell: cubic
        cell_pdb(2) = 90.0
        cell_pdb(4) = 90.0
        cell_pdb(5) = 90.0

    else if( icell == 2 ) then

        cell_pdb(1) = cfgs(ib)%vec%latvector(1,1) ! pbc_cell: orthorhombic
        cell_pdb(3) = cfgs(ib)%vec%latvector(2,2) ! pbc_cell: orthorhombic
        cell_pdb(6) = cfgs(ib)%vec%latvector(3,3) ! pbc_cell: orthorhombic
        cell_pdb(2) = 90.0
        cell_pdb(4) = 90.0
        cell_pdb(5) = 90.0

    !else if( icell > 2 ) then
    else if( icell /= 0 ) then ! non-orthorhombic cells, angles /= 90, so let's stick to the cell matrix
    !AB: some sources refer to unit cells in the lower diagonal matrix representation (instead of PDB format: a,b,c + angles)

        cell_pdb(1) = cfgs(ib)%vec%latvector(1,1) ! pbc_cell: lower diagonal
        cell_pdb(2) = cfgs(ib)%vec%latvector(2,1) ! pbc_cell: lower diagonal
        cell_pdb(3) = cfgs(ib)%vec%latvector(2,2) ! pbc_cell: lower diagonal
        cell_pdb(4) = cfgs(ib)%vec%latvector(3,1) ! pbc_cell: lower diagonal
        cell_pdb(5) = cfgs(ib)%vec%latvector(3,2) ! pbc_cell: lower diagonal
        cell_pdb(6) = cfgs(ib)%vec%latvector(3,3) ! pbc_cell: lower diagonal

    !else if( icell > ??? ) then ! crystallographic data (CRYST1 in PDB)
    end if

    rscale = 10.0 !AB: Angstrom -> nm (coordinate/distance re-scaling needed?)

    fail = 0
    k = 0

    if( iform == 0 ) then

        do  i = 1, cfgs(ib)%num_mols

            !AB: unwrap molecule (undo PBC) if it is rigid or has structure due to bonds
            if( cfgs(ib)%mols(i)%rigid_body .or. cfgs(ib)%mols(i)%blist%npairs > 0 ) then

                !AB: before unwrapping molecule store its atom positions and com
                call store_molecule_pos(ib,i)

                if( myjob%useorthogonal ) then

                    !if( uniq_mol(ml)%blist%npairs > 0 ) &
                    call reform_molecule_ortho(cfgs(ib)%mols(i),cfgs(ib)%vec%latvector)
                    call mol_com(cfgs(ib)%mols(i))

                else

                    call full_mol_com(cfgs(ib)%mols(i), cfgs(ib)%vec%invlat,cfgs(ib)%vec%latvector)

                end if

                mcom(:) = cfgs(ib)%mols(i)%rcom(:)

                call pbc_cart_vec_calc(ib,mcom)

                mcom(:) = mcom(:) - cfgs(ib)%mols(i)%rcom(:)

                call shift_molecule(cfgs(ib)%mols(i), mcom)

            end if

            do  j = 1, cfgs(ib)%mols(i)%natom
                k = k + 1
                x_dcd4(k) = cfgs(ib)%mols(i)%atms(j)%rpos(1) !*rscale
                y_dcd4(k) = cfgs(ib)%mols(i)%atms(j)%rpos(2) !*rscale
                z_dcd4(k) = cfgs(ib)%mols(i)%atms(j)%rpos(3) !*rscale
            end do

            if( cfgs(ib)%mols(i)%rigid_body .or. cfgs(ib)%mols(i)%blist%npairs > 0 ) &
                call restore_molecule_pos(ib,i)

        end do

    else

        do  i = 1, cfgs(ib)%num_mols

            !AB: unwrap molecule (undo PBC) if it is rigid or has structure due to bonds
            if( cfgs(ib)%mols(i)%rigid_body .or. cfgs(ib)%mols(i)%blist%npairs > 0 ) then

                !AB: before unwrapping molecule store its atom positions and com
                call store_molecule_pos(ib,i)

                if( myjob%useorthogonal ) then

                    !if( uniq_mol(ml)%blist%npairs > 0 ) &
                    call reform_molecule_ortho(cfgs(ib)%mols(i),cfgs(ib)%vec%latvector)
                    call mol_com(cfgs(ib)%mols(i))

                else

                    call full_mol_com(cfgs(ib)%mols(i), cfgs(ib)%vec%invlat,cfgs(ib)%vec%latvector)

                end if

                mcom(:) = cfgs(ib)%mols(i)%rcom(:)

                call pbc_cart_vec_calc(ib,mcom)

                mcom(:) = mcom(:) - cfgs(ib)%mols(i)%rcom(:)

                call shift_molecule(cfgs(ib)%mols(i), mcom)

            end if

            do  j = 1, cfgs(ib)%mols(i)%natom
                k = k + 1
                x_dcd8(k) = cfgs(ib)%mols(i)%atms(j)%rpos(1) !*rscale
                y_dcd8(k) = cfgs(ib)%mols(i)%atms(j)%rpos(2) !*rscale
                z_dcd8(k) = cfgs(ib)%mols(i)%atms(j)%rpos(3) !*rscale
            end do

            if( cfgs(ib)%mols(i)%rigid_body .or. cfgs(ib)%mols(i)%blist%npairs > 0 ) &
                call restore_molecule_pos(ib,i)

        end do

    end if

    if( natoms /= k ) then

        nwarn3 = nwarn3+1

        write(dcd_string,'(i10,a,i10)') k,' =?= ',natoms

        if( nwarn3 > 1 .and. nwarn3 < 11 ) then 

            !AB: this cry is not suitable for GCMC - introduce a check for this case!!!
            call cry(uout,'', &
                 "WARNING: # atoms variation detected in DCD file '"//trim(dcd_fname)// &
                 "' : "//trim(adjustL(dcd_string))//" - must be GCMC; padding with (0,0,0) !!!",0)

        else if( nwarn3 == 11 ) then 

            call cry(uout,'', &
                 "WARNING: # atoms variation detected in DCD file '"//trim(dcd_fname)// &
                 "' - keeping silent about more such events !!!",0)
        end if

    end if

    if( natoms < matoms ) then

        nwarn4 = nwarn4+1

        write(dcd_string,'(i10,a,i10)') k,' =?= ',matoms

        if( nwarn4 < 11 ) then 

            call cry(uout,'', &
                 "WARNING: # atoms < maximum # atoms for DCD file '"//trim(dcd_fname)// &
                 "' : "//trim(adjustL(dcd_string))//" - must be GCMC; padding with (0,0,0) !!!",0)

        else if( nwarn4 == 11 ) then 

            call cry(uout,'', &
                 "WARNING: # atoms < maximum # atoms for DCD file '"//trim(dcd_fname)// &
                 "' - keeping silent about more such events !!!",0)
        end if

    else if( matoms < natoms ) then

        write(dcd_string,'(i10,a,i10)') natoms,' > ',matoms

        call cry(uout,'', &
                 "ERROR: # atoms overflow detected while writing DCD file '"//trim(dcd_fname)// &
                 "' : "//trim(adjustL(dcd_string))//" - increase max # atoms in CONFIG file !!!",999)

    end if

    natoms = k

    !AB: all virtual atoms are fixed at the origin (0,0,0)

    k = k + 1

    if( iform == 0 ) then

        do  i = k, matoms
            x_dcd4(i) = 0.0
            y_dcd4(i) = 0.0
            z_dcd4(i) = 0.0
        end do

        call dcd_write_next_frame4(iodcd, dcd_fname, iopen, iclose, &
                                   matoms, mframes, icell, x_dcd4, y_dcd4, z_dcd4, cell_pdb)

    else

        do  i = k, matoms
            x_dcd8(i) = 0.0
            y_dcd8(i) = 0.0
            z_dcd8(i) = 0.0
        end do

        call dcd_write_next_frame8(iodcd, dcd_fname, iopen, iclose, &
                                   matoms, mframes, icell, x_dcd8, y_dcd8, z_dcd8, cell_pdb)

    end if

    if( iclose == 0 ) return

    fail = 0
    if( iform == 0 ) then
        deallocate(x_dcd4, stat = fail(1))
        deallocate(y_dcd4, stat = fail(2))
        deallocate(z_dcd4, stat = fail(3))

        if( any(fail > 0) ) &
            call cry(uout,'', &
                     "ERROR: failed deallocating memory after DCD arrays '"//&
                    &"'- full stop !!!",999)
    else
        deallocate(x_dcd8, stat = fail(1))
        deallocate(y_dcd8, stat = fail(2))
        deallocate(z_dcd8, stat = fail(3))

        if( any(fail > 0) ) &
            call cry(uout,'', &
                     "ERROR: failed deallocating memory after DCD arrays '"//&
                    &"'- full stop !!!",999)
    endif

end subroutine write_config_dcd

!AB: -> writing trajectory in DCD format

!AB: -->> dealing with trajectories in DCD format

!> @brief
!> - writes out configuration sample of up to 10 atom positions
!> @using
!> - `kinds_f90`
!> - `constants_module`, only : uout
!> - `species_module`, only : set_species_type, element
!> - `latticevectors_module`, only : printlatticevectors, setrcpvec, printrcpvec
!> - `atom_module`, only : write_atom_pos
!> - `comms_mpi_module`, only : master

!> prints out basis positions etc
subroutine printbasis(ib)

    use constants_module, only : uout
    use species_module, only : set_species_type, element
    use latticevectors_module, only : printlatticevectors, setrcpvec, printrcpvec
    use atom_module, only : write_atom_pos
    use comms_mpi_module, only : master

    !implicit none

        !> replica configuration & cell identifier
    integer, intent(in) :: ib

    ! auxilary variables
    integer :: lb, i, j, ii
    character :: tag*1

    ! lattice parameters
    call printlatticevectors(cfgs(ib)%vec)

    ! calculate recip lattice vectors for start and print
    call setrcpvec(cfgs(ib)%vec)
    call printrcpvec(cfgs(ib)%vec)

    ii = 0

    if (master) write(uout,"(/,1x,50('-'),/,1x,a,/,1x,50('-'),/)") &
                      "configuration sample (up to 10 atoms):"

    loop: do i = 1, cfgs(ib)%num_mols

        if (master) write(uout, "(1x,'molecule',1x,a8)") cfgs(ib)%mols(i)%molname

        do j = 1, cfgs(ib)%mols(i)%natom

            ii = ii + 1

            if( ii > 10 ) exit loop

            lb = cfgs(ib)%mols(i)%atms(j)%atlabel
            tag = set_species_type(cfgs(ib)%mols(i)%atms(j)%atype)

            if (master) then

                write(uout, "(1x,a8,1x,2a,1x,i8,a,i3)") element(lb)," ",tag, ii,"   type ", lb
                call write_atom_pos(cfgs(ib)%mols(i)%atms(j), uout)

            endif

        enddo

    enddo loop

end subroutine printbasis

!> @brief
!> - check for the total charge in replica 'ib'
!> @using
!> - `kinds_f90`
!> - `constants_module`, only : uout
!> - `comms_mpi_module`, only : is_parallel, gsum

!> check for the total charge in a replica
subroutine cell_charge(ib)

    !use kinds_f90
    use constants_module, only : uout
    use comms_mpi_module, only : is_parallel, gsum

    !implicit none

        !> replica configuration & cell identifier
    integer, intent(in) :: ib

    integer :: i, j

    real(kind = wp) :: tot_charge

    tot_charge = 0

    do i = 1, cfgs(ib)%num_mols

        do j = 1, cfgs(ib)%mols(i)%natom
            tot_charge = tot_charge + cfgs(ib)%mols(i)%atms(j)%charge
        enddo

    enddo

    if( is_parallel ) call gsum(tot_charge)

    if( tot_charge /= 0.0_wp ) then
        call error(201)
    endif

end subroutine cell_charge

!> @brief
!> - resets atom positions to within the primary cell for replica `ib`
!> @using
!> - `kinds_f90`
!> - `constants_module`, only : uout
!> - `latticevectors_type`
!> - `latticevectors_module`
!> - `molecule_module`, only : reform_molecule_frac, full_reform_molecule, carttofrac_mol, fractocart_mol
!> - `slit_module`, only : in_bulk

!> puts all atoms back in the primary simulation cell
subroutine pbc_simulation_cell(ib)

    !use kinds_f90
    use constants_module, only : DZERO, uout
    use latticevectors_type
    use latticevectors_module
    use molecule_module, only : full_reform_molecule 
    use slit_module, only : in_bulk

    !implicit none

        !> replica configuration & cell identifier
    integer, intent(in) :: ib

    integer :: im, ia, ja, nz_out, nz_fail, ntry
    real(kind=wp) :: zia, zij, zmax

    logical :: is_reform = .false.

    character :: num1*8, num2*8, num3*8, text*32

    if( in_bulk ) then

      do im = 1, cfgs(ib)%num_mols

        do ia = 1, cfgs(ib)%mols(im)%natom

            ! correct atom position within box (PBC/MIC)
            call pbc_atom(ib, im, ia)

        end do

      end do

    else !AB: i.e. in_slit (only orthogonal PBC allowed)

      ntry = 0

101   nz_out  = 0
      nz_fail = 0
      zmax = 0.0_wp

      do im = 1, cfgs(ib)%num_mols

        !AB: chek if the molecule is rigid or has a bonded structure/topology to it
        is_reform = cfgs(ib)%mols(im)%rigid_body .and. cfgs(ib)%mols(im)%blist%npairs > 0

        if( is_reform ) call full_reform_molecule(cfgs(ib)%mols(im), &
                             cfgs(ib)%vec%invlat,cfgs(ib)%vec%latvector)

        do ia = 1, cfgs(ib)%mols(im)%natom

            ! correct atom position within box (PBC/MIC)
            call pbc_atom(ib, im, ia)

            !AB: check molecule integrity in Z dimension
            if( is_reform ) then

                zia = cfgs(ib)%mols(im)%atms(ia)%rpos(3)

                if( abs(zia) > cfgs(ib)%vec%latvector(3,3)*0.5_wp ) then 
                    nz_out = nz_out+1
                    zmax = max(zmax,abs(zia))
                end if

                do ja = 1, ia-1

                   zij = abs(zia-cfgs(ib)%mols(im)%atms(ja)%rpos(3))

                   if( zij > cfgs(ib)%vec%latvector(3,3)*0.5_wp ) then

                       nz_fail = nz_fail+1

                       write(num1,'(i8)') ia
                       write(num2,'(i8)') ja
                       write(num3,'(i8)') im
                       write(text,'(e14.7,a,e14.7)') zij,'  > ',cfgs(ib)%vec%latvector(3,3)*0.5_wp

                       call cry(uout,'', &
                           "WARNING: atoms  "//trim(adjustL(num1))//"  &  "//trim(adjustL(num2))//&
                          &"  in molecule  "//trim(adjustL(num3))//&
                          &"  more than half Z-dim apart: "//text//" !!!",0)

                   end if

                end do

            end if

        end do

      end do

      if( nz_out > 0 ) then
          write(num3,'(i8)') nz_out
          write(text,'(e14.7,a,e14.7)') zmax,'  > ',cfgs(ib)%vec%latvector(3,3)*0.5_wp
          call cry(uout,'', &
               "WARNING: found  "//trim(adjustL(num3))//&
              &"  atoms in reformed molecules outside the slit: "//text//" - resetting Z-size !!!",0)

          cfgs(ib)%vec%latvector(3,3) = zmax*2.0_wp + DZERO
          cfgs(ib)%vec%invlat(3,3) = 1.0_wp/cfgs(ib)%vec%latvector(3,3)
          ntry = ntry+1
          if( ntry < 2 ) goto 101

      end if

      if( nz_fail > 0 ) then
          write(num3,'(i8)') nz_fail
          call cry(uout,'', &
               "ERROR: found  "//trim(adjustL(num3))//&
              &"  atom pairs in molecules more than half slit Z-size apart !!!",999)
      end if

    end if

end subroutine pbc_simulation_cell


!> re-calculate COM:s for all rigid molecules in simulation cell
subroutine set_rigid_mol_coms(ib)

    use molecule_module, only : full_mol_com

        !> replica configuration & cell identifier
    integer, intent(in) :: ib

    integer :: im

    real(kind=wp) :: com_var

    do im = 1, cfgs(ib)%num_mols
        if (cfgs(ib)%mols(im)%rigid_body) then
            call set_mol_com(ib, im)
        end if
    enddo

end subroutine set_rigid_mol_coms


!> re-calculate COM:s for all molecules in simulation cell
subroutine set_mol_coms(ib)

    use molecule_module, only : full_mol_com, reform_molecule_frac, mol_com &
                               ,carttofrac_mol,fractocart_mol

        !> replica configuration & cell identifier
    integer, intent(in) :: ib

    integer :: im

    real(kind=wp) :: com_var

    do im = 1, cfgs(ib)%num_mols

        !AB: recalculate molecule mass and COM from scratch (atom positions remain intact!)
        !AB: COM calculation done in one go (loop) over the molecule's atoms
        call full_mol_com(cfgs(ib)%mols(im), cfgs(ib)%vec%invlat, cfgs(ib)%vec%latvector)

        call store_molecule_pos(ib, im)

        cycle

        !AB: the rest is needed only for testing/debugging

        !AB: recalculate molecule COM from scratch (mass is not updated! atom positions affected!)
        !call set_mol_com(ib, im)

        !AB: the original way of recalculating molecule COM
        call carttofrac_mol(cfgs(ib)%mols(im), cfgs(ib)%vec%invlat)
        call reform_molecule_frac(cfgs(ib)%mols(im))
        call fractocart_mol(cfgs(ib)%mols(im), cfgs(ib)%vec%latvector)
        call mol_com(cfgs(ib)%mols(im))

        com_var = max(abs( cfgs(ib)%mols(im)%rcom(1)-cfgs(ib)%mols(im)%store_rcom(1) ), &
                      abs( cfgs(ib)%mols(im)%rcom(2)-cfgs(ib)%mols(im)%store_rcom(2) ), &
                      abs( cfgs(ib)%mols(im)%rcom(3)-cfgs(ib)%mols(im)%store_rcom(3) ))

        if( com_var > 1.d-8 ) then

            write(*,*)
            write(*,*)"COM calculations are inconsistent for molecule ",im," in box ",ib, com_var

        endif

        call restore_molecule_pos(ib, im)

    enddo

end subroutine set_mol_coms


!> fully recalculates the centre of mass of a molecule while unwrapping the atom positions
subroutine set_mol_com(ib, im)

    use molecule_module, only : mol_com, full_reform_molecule 
    !, reform_molecule_frac, carttofrac_mol, fractocart_mol

        !> replica configuration/cell & molecule identifiers
    integer, intent(in) :: ib,im

    integer :: j
    real(kind = wp) :: orig(3)

    if( cfgs(ib)%vec%is_orthogonal ) then

        orig(:) = cfgs(ib)%mols(im)%atms(1)%rpos(:)

        do j = 2, cfgs(ib)%mols(im)%natom

            cfgs(ib)%mols(im)%atms(j)%rpos(1) = cfgs(ib)%mols(im)%atms(j)%rpos(1) - orig(1)
            cfgs(ib)%mols(im)%atms(j)%rpos(2) = cfgs(ib)%mols(im)%atms(j)%rpos(2) - orig(2)
            cfgs(ib)%mols(im)%atms(j)%rpos(3) = cfgs(ib)%mols(im)%atms(j)%rpos(3) - orig(3)

            call pbc_ortho_coor_full(ib,cfgs(ib)%mols(im)%atms(j)%rpos(1), &
                                        cfgs(ib)%mols(im)%atms(j)%rpos(2), &
                                        cfgs(ib)%mols(im)%atms(j)%rpos(3))

            cfgs(ib)%mols(im)%atms(j)%rpos(1) = cfgs(ib)%mols(im)%atms(j)%rpos(1) + orig(1)
            cfgs(ib)%mols(im)%atms(j)%rpos(2) = cfgs(ib)%mols(im)%atms(j)%rpos(2) + orig(2)
            cfgs(ib)%mols(im)%atms(j)%rpos(3) = cfgs(ib)%mols(im)%atms(j)%rpos(3) + orig(3)

        end do

        call mol_com(cfgs(ib)%mols(im))

    else

        !AB: the original way of recalculating molecule COM
        !call carttofrac_mol(cfgs(ib)%mols(im), cfgs(ib)%vec%invlat)
        !call reform_molecule_frac(cfgs(ib)%mols(im))
        !call fractocart_mol(cfgs(ib)%mols(im), cfgs(ib)%vec%latvector)

        !AB: equivalent COM calculation but done in one go (loop) over the molecule's atoms
        call full_reform_molecule(cfgs(ib)%mols(im), cfgs(ib)%vec%invlat, cfgs(ib)%vec%latvector)
        call mol_com(cfgs(ib)%mols(im))

    end if

end subroutine set_mol_com


!> fully recalculates the centre of mass of a molecule while unwrapping the atom positions
subroutine check_mol_com(ib, im)

    use molecule_module, only : full_mol_com

        !> replica configuration/cell & molecule identifiers
    integer, intent(in) :: ib,im

    real(kind=wp), dimension(3) :: rcom

    real(kind=wp) :: com_var

    rcom(:) = cfgs(ib)%mols(im)%rcom(:)

    !AB: recalculate molecule mass and COM from scratch (atom positions remain intact!)
    call full_mol_com(cfgs(ib)%mols(im), cfgs(ib)%vec%invlat, cfgs(ib)%vec%latvector)

    !AB: recalculate molecule COM from scratch (mass is not updated! atom positions affected!)
    !call set_mol_com(ib, im)

    com_var = max(abs( cfgs(ib)%mols(im)%rcom(1)-rcom(1) ), &
                  abs( cfgs(ib)%mols(im)%rcom(2)-rcom(2) ), &
                  abs( cfgs(ib)%mols(im)%rcom(3)-rcom(3) ))

    if( com_var > 1.d-8 ) then

        write(*,*)
        write(*,*)"Running COM calculations are inconsistent for molecule ",im," in box ",ib, com_var

    endif

end subroutine check_mol_com

!> re-calculate COM:s for all molecules in simulation cell
subroutine check_mol_coms(ib)

        !> replica configuration & cell identifier
    integer, intent(in) :: ib

    ! molecule identifier
    integer :: im

    do im = 1, cfgs(ib)%num_mols

       call check_mol_com(ib, im)

    enddo

end subroutine check_mol_coms

!> fully recalculates the centre of mass of a molecule while unwrapping the atom positions
logical function is_mol_com_ok(ib, im)

    use molecule_module, only : full_mol_com

        !> replica configuration/cell & molecule identifiers
    integer, intent(in) :: ib,im

    real(kind=wp), dimension(3) :: rcom

    real(kind=wp) :: com_var
    
    is_mol_com_ok = .true.
    
    rcom(:) = cfgs(ib)%mols(im)%rcom(:)

    !AB: recalculate molecule mass and COM from scratch (atom positions remain intact!)
    call full_mol_com(cfgs(ib)%mols(im), cfgs(ib)%vec%invlat, cfgs(ib)%vec%latvector)

    !AB: recalculate molecule COM from scratch (mass is not updated! atom positions affected!)
    !call set_mol_com(ib, im)

    com_var = max(abs( cfgs(ib)%mols(im)%rcom(1)-rcom(1) ), &
                  abs( cfgs(ib)%mols(im)%rcom(2)-rcom(2) ), &
                  abs( cfgs(ib)%mols(im)%rcom(3)-rcom(3) ))

    if( com_var > 1.d-8 ) then
    
        is_mol_com_ok = .false.

        write(*,*)
        write(*,*)"Running COM calculations are inconsistent for molecule ",im," in box ",ib, com_var

    endif

end function is_mol_com_ok


!> re-calculate COM:s for all molecules in simulation cell
subroutine check_mol_coms_ok(ib,safe)

        !> replica configuration & cell identifier
    integer, intent(in) :: ib
    
    logical, intent(out) :: safe

    ! molecule identifier
    integer :: im

    safe = .true.

    do im = 1, cfgs(ib)%num_mols

       !call check_mol_com(ib, im)
       safe = ( safe .and. is_mol_com_ok(ib, im) )

    enddo

end subroutine check_mol_coms_ok

!> @brief
!> - PBC minimum image convention: puts a \b single \b coordinate back in simulation cell
!> - works on a fractional coordinate (come here after accounting for the cell geometry)
!> - implemented via \b `if(` \b ` rc > 0.5 ) then rc = rc - 1.0 else rc = rc + 1.0 endif`
!> @using
!> - `kinds_f90`

!> PBC/MIC for a single fractional coordinate, via \b `if(` \b ` rc > 0.5 ) then ... else ...`
subroutine pbc_coordinate_call(rc)

    !use kinds_f90

    !implicit none

        !> fractional coordinate: x/y/z
    real (kind = wp), intent(inout) :: rc

    !AB: assuming that a particle did not diffuse away from the origin more than one box length
    !AB: i.e. it relies on putting particles back into the primary box after every single move!
    if( rc > 0.5_wp ) then
        rc = rc - 1.0_wp
    else if( rc <= -0.5_wp ) then
        rc = rc + 1.0_wp
    end if

end subroutine pbc_coordinate_call

!> @brief
!> - PBC minimum image convention: puts a \b single \b atom back in simulation cell
!> - works on fractional corrdinates (come here after accounting for the cell geometry)
!> - implemented via `'call pbc_coordinate_call(..)'` for every coordinate
!> - implements bulk/slit PBC via \b `if( in_bulk ) call pbc_coordinate_call(`\b `z)`
!> @using
!> - `kinds_f90`
!> - `slit_module, only : in_bulk`

!> PBC/MIC for a single atom, via `'call pbc_coordinate_call(..)'` for every coordinate
subroutine pbc_atom_pos_call(rpos)

    !use kinds_f90
    use slit_module, only : in_bulk

    !implicit none

        !> fractional coordinates of atom (x,y,z)
    real (kind = wp), dimension (3), intent(inout) :: rpos

    integer :: k

    call pbc_coordinate_call(rpos(1))
    call pbc_coordinate_call(rpos(2))

    ! don't use PBC in Z for planar pore geometry 
    if( in_bulk ) call pbc_coordinate_call(rpos(3))

end subroutine pbc_atom_pos_call

!> @brief
!> - PBC minimum image convention: puts a \b single \b coordinate back in simulation cell
!> - works on a fractional coordinate (come here after accounting for the cell geometry)
!> - implemented via direct calculation of \b `rc = rc - anint(rc)`\b `*on_off`
!> - \b `on_off` switches on/off the PBC/MIC for a single coordinate (e.g. Z-PBC off in slit)
!> @using
!> - `kinds_f90`

!> PBC/MIC for a single fractional coordinate, via \b `x = x - anint(x)`\b `*on_off`
real(kind=wp) pure function pbc_single_calc(rc, on_off)

    !use kinds_f90

    !implicit none

        !> fractional coordinate: x/y/z
    real (kind = wp), intent(in) :: rc

        !> factor to switch Z-PBC on/off (in bulk/slit)
    real (kind = wp), intent(in) :: on_off

    !AB: one way is via ANINT(X), which takes and returns default Real(*) values - better use DNINT(X)!
    !AB: if X > 0.0, then ANINT(X) := AINT(X+0.5); else if X <= 0.0, then ANINT(X) := AINT(X-0.5). 

    pbc_single_calc = rc - anint(rc)*on_off

end function pbc_single_calc


!> @brief
!> - PBC minimum image convention: puts a \b single \b coordinate back in simulation cell
!> - works on a fractional coordinate (come here after accounting for the cell geometry)
!> - implemented via direct calculation of \b `x = x - anint(x)`\b `*on_off`
!> - \b `on_off` switches no/off the PBC/MIC for a single coordinate
!> @using
!> - `kinds_f90`

!> PBC/MIC for a single fractional coordinate, via \b `x = x - anint(x)`\b `*on_off`
subroutine pbc_coordinate_calc(rc, on_off)

    !use kinds_f90

    !implicit none

        !> fractional coordinate: x/y/z
    real (kind = wp), intent(inout) :: rc

        !> factor to switch PBC on/off (bulk/slit for Z)
    real (kind = wp), intent(in) :: on_off

    !AB: one way is via ANINT(X), which takes and returns default Real(*) values - better use DNINT(X)!
    !AB: if X > 0.0, then ANINT(X) := AINT(X+0.5); else if X <= 0.0, then ANINT(X) := AINT(X-0.5). 

    if( abs(rc) > 0.5_wp ) rc = rc - anint(rc)*on_off

end subroutine pbc_coordinate_calc

!> @brief
!> - PBC minimum image convention 
!> - PBC interparticle distance corrector for orthogonal cell matrix (lattice vectors)
!> - works directly on a \b Cartesian \b coordinates
!> @using
!> - `kinds_f90`
!> - `slit_module, only : in_bulk`

subroutine pbc_ortho_coor_calc(ib, rx, ry, rz)

    !use kinds_f90
    use slit_module, only : in_bulk

        !> Cartesian coordinates: x,y,z
    real (kind = wp), intent(inout) :: rx, ry, rz

        !> replica configuration & cell identifier
    integer, intent(in) :: ib

!AB: This assumes that a particle did not diffuse away from the origin more than one box length
!AB: i.e. it relies on putting particles back into the primary box after every single move!
!JG: Which switching on orthogonal does.

!AB: this does not apply to a CONFIG not satisfying the above condition (whence 'pbc_ortho_coor_full' below)

    if( rx > cfgs(ib)%vec%halfvector(1) ) then
        rx = rx - cfgs(ib)%vec%fullvector(1) 
    !else if( rx <= -cfgs(ib)%vec%halfvector(1) ) then
    else if( rx < -cfgs(ib)%vec%halfvector(1) ) then
        rx = cfgs(ib)%vec%fullvector(1) + rx
    end if

    if( ry > cfgs(ib)%vec%halfvector(2) ) then
        ry = ry - cfgs(ib)%vec%fullvector(2)
    !else if( ry <= -cfgs(ib)%vec%halfvector(2) ) then
    else if( ry < -cfgs(ib)%vec%halfvector(2) ) then
        ry = cfgs(ib)%vec%fullvector(2) + ry
    end if

    if( in_bulk ) then
        if( rz > cfgs(ib)%vec%halfvector(3) ) then
            rz = rz - cfgs(ib)%vec%fullvector(3)
        !else if(rz <= -cfgs(ib)%vec%halfvector(3) ) then
        else if(rz < -cfgs(ib)%vec%halfvector(3) ) then
            rz = cfgs(ib)%vec%fullvector(3) + rz
        end if
    end if

end subroutine pbc_ortho_coor_calc

subroutine pbc_ortho_coor_full(ib, rx, ry, rz)

    !use kinds_f90
    use slit_module, only : im_bulk

        !> replica configuration & cell identifier
    integer, intent(in) :: ib

        !> Cartesian coordinates: x,y,z
    real(kind = wp), intent(inout) :: rx, ry, rz

    real(kind = wp) :: xx, yy, zz

    xx = rx*cfgs(ib)%vec%invlat(1,1)
    yy = ry*cfgs(ib)%vec%invlat(2,2)
    zz = rz*cfgs(ib)%vec%invlat(3,3)

    if( abs(xx) > 0.5_wp ) &
        rx = (xx - anint(xx)) * cfgs(ib)%vec%latvector(1,1)

    if( abs(yy) > 0.5_wp ) &
        ry = (yy - anint(yy)) * cfgs(ib)%vec%latvector(2,2)
        
    if( abs(zz) > 0.5_wp ) &
        rz = (zz - anint(zz)*real(im_bulk,wp)) * cfgs(ib)%vec%latvector(3,3)

end subroutine pbc_ortho_coor_full

!> @brief
!> - PBC minimum image convention: puts a \b single \b atom back in simulation cell
!> - works on fractional corrdinates (come here after accounting for the cell geometry)
!> - implemented via direct calculation of \b `x = x - anint(x)` and so on
!> - implements bulk/slit PBC directly via \b `z = z - anint(z)`\b `*im_bulk`
!> @using
!> - `kinds_f90`
!> - `slit_module, only : im_bulk`

!> PBC/MIC for a single atom, by direct calculation: `'x = x - anint(x)'` etc
subroutine pbc_atom_pos_calc(rx, ry, rz)

    !use kinds_f90
    use slit_module, only : im_bulk

    !implicit none

        !> fractional coordinates of atom (x,y,z)
    real (kind = wp), intent(inout) :: rx, ry, rz

    !AB: with gfortran (Linux/Ubuntu) there is no difference between anint() & dnint()
    !AB: if:s prevent applying costly anint(..) when unnecessary, i.e. make it faster

    if( abs(rx) > 0.5_wp ) &
        rx = rx - anint(rx)

    if( abs(ry) > 0.5_wp ) & 
        ry = ry - anint(ry)

    if( abs(rz) > 0.5_wp ) &
        rz = rz - anint(rz)*real(im_bulk,wp)

end subroutine pbc_atom_pos_calc

!> @brief
!> - PBC minimum image convention: puts a \b Cartesian \b position \b vector back in simulation cell
!> - first converts the coordinates: Cartesian -> fractional, accounting for the cell geometry
!> - applies PBC/MIC and converts the coordinates back to Cartesian form
!> @using
!> - `kinds_f90`

!> puts an atom back in simulation cell (PBC/MIC upon Cartesian -> fractional conversion)
subroutine pbc_cart_vec_calc(ib, rpos)

    !use kinds_f90
    use slit_module, only : im_bulk

    !implicit none

        !> configuration box id/index (replica)
    integer, intent(in) :: ib

        !> Cartesian position vector
    real (kind = wp), intent(inout) :: rpos(3)

    real (kind = wp) :: rx, ry, rz, xx, yy, zz

    ! the Cartesian coordinates (within the cell geomentry)
    rx = rpos(1)
    ry = rpos(2)
    rz = rpos(3)

    !JG: If orthogonal lattice vectors use simple check
    if(cfgs(ib)%vec%is_orthogonal) then

        call pbc_ortho_coor_calc(ib, rx, ry, rz)
        !call pbc_ortho_coor_full(ib, rx, ry, rz)

        rpos(1) = rx
        rpos(2) = ry
        rpos(3) = rz
   
    else

        ! get the fractional coordinates

        xx = cfgs(ib)%vec%invlat(1,1) * rx &
           + cfgs(ib)%vec%invlat(2,1) * ry &
           + cfgs(ib)%vec%invlat(3,1) * rz

        yy = cfgs(ib)%vec%invlat(1,2) * rx &
           + cfgs(ib)%vec%invlat(2,2) * ry &
           + cfgs(ib)%vec%invlat(3,2) * rz

        zz = cfgs(ib)%vec%invlat(1,3) * rx &
           + cfgs(ib)%vec%invlat(2,3) * ry &
           + cfgs(ib)%vec%invlat(3,3) * rz

        ! apply PBC/MIC to fractional coordinates:

        !call pbc_atom_pos_calc(xx, yy, zz)

        if( abs(xx) > 0.5_wp ) &
            xx = xx - anint(xx)

        if( abs(yy) > 0.5_wp ) & 
            yy = yy - anint(yy)

        if( abs(zz) > 0.5_wp ) &
            zz = zz - anint(zz)*real(im_bulk,wp)

        ! get back the Cartesian coordinates (re-apply the cell geomentry)

        rpos(1) = cfgs(ib)%vec%latvector(1,1) * xx &
                + cfgs(ib)%vec%latvector(2,1) * yy &
                + cfgs(ib)%vec%latvector(3,1) * zz

        rpos(2) = cfgs(ib)%vec%latvector(1,2) * xx &
                + cfgs(ib)%vec%latvector(2,2) * yy &
                + cfgs(ib)%vec%latvector(3,2) * zz

        rpos(3) = cfgs(ib)%vec%latvector(1,3) * xx &
                + cfgs(ib)%vec%latvector(2,3) * yy &
                + cfgs(ib)%vec%latvector(3,3) * zz

    end if

end subroutine pbc_cart_vec_calc

!> @brief
!> - PBC minimum image convention: puts \b Cartesian \b position \b coordinates back in simulation cell
!> - first converts the corrdinates: Cartesian -> fractional, accounting for the cell geometry
!> - applies PBC/MIC and converts the coordinates back to Cartesian form
!> @using
!> - `kinds_f90`

!> puts an atom back in simulation cell (PBC/MIC upon Cartesian -> fractional conversion)
subroutine pbc_cart_pos_calc(ib, rx, ry, rz)

    !use kinds_f90
    use slit_module, only : im_bulk

    !implicit none

        !> configuration box id/index (replica)
    integer, intent(in) :: ib

        !> Cartesian coordinates
    real (kind = wp), intent(inout) :: rx, ry, rz

        ! fractional coordinates
    real (kind = wp) :: xx, yy, zz

    !JG: If orthogonal lattice vectors use simple check
    if(cfgs(ib)%vec%is_orthogonal) then

        call pbc_ortho_coor_calc(ib, rx, ry, rz)
        !call pbc_ortho_coor_full(ib, rx, ry, rz)
   
    else

        ! get the fractional coordinates

        xx = cfgs(ib)%vec%invlat(1,1) * rx &
           + cfgs(ib)%vec%invlat(2,1) * ry &
           + cfgs(ib)%vec%invlat(3,1) * rz

        yy = cfgs(ib)%vec%invlat(1,2) * rx &
           + cfgs(ib)%vec%invlat(2,2) * ry &
           + cfgs(ib)%vec%invlat(3,2) * rz

        zz = cfgs(ib)%vec%invlat(1,3) * rx &
           + cfgs(ib)%vec%invlat(2,3) * ry &
           + cfgs(ib)%vec%invlat(3,3) * rz

        ! apply PBC/MIC to fractional coordinates:

        !call pbc_atom_pos_calc(xx, yy, zz)

        if( abs(xx) > 0.5_wp ) &
            xx = xx - anint(xx)

        if( abs(yy) > 0.5_wp ) & 
            yy = yy - anint(yy)

        if( abs(zz) > 0.5_wp ) &
            zz = zz - anint(zz)*real(im_bulk,wp)

        ! get back the Cartesian coordinates (re-apply the cell geomentry)

        rx = cfgs(ib)%vec%latvector(1,1) * xx &
           + cfgs(ib)%vec%latvector(2,1) * yy &
           + cfgs(ib)%vec%latvector(3,1) * zz

        ry = cfgs(ib)%vec%latvector(1,2) * xx &
           + cfgs(ib)%vec%latvector(2,2) * yy &
           + cfgs(ib)%vec%latvector(3,2) * zz

        rz = cfgs(ib)%vec%latvector(1,3) * xx &
           + cfgs(ib)%vec%latvector(2,3) * yy &
           + cfgs(ib)%vec%latvector(3,3) * zz

    end if

end subroutine pbc_cart_pos_calc

!> puts an atom back in simulation cell (PBC/MIC upon Cartesian -> fractional conversion)
subroutine pbc_cart2frac_pos_calc(ib, rx, ry, rz)

    !use kinds_f90
    use slit_module, only : im_bulk

    !implicit none

        !> configuration box id/index (replica)
    integer, intent(in) :: ib

        !> Cartesian coordinates
    real (kind = wp), intent(inout) :: rx, ry, rz

        ! fractional coordinates
    real (kind = wp) :: xx, yy, zz

    ! get the fractional coordinates

    xx = cfgs(ib)%vec%invlat(1,1) * rx &
       + cfgs(ib)%vec%invlat(2,1) * ry &
       + cfgs(ib)%vec%invlat(3,1) * rz

    yy = cfgs(ib)%vec%invlat(1,2) * rx &
       + cfgs(ib)%vec%invlat(2,2) * ry &
       + cfgs(ib)%vec%invlat(3,2) * rz

    zz = cfgs(ib)%vec%invlat(1,3) * rx &
       + cfgs(ib)%vec%invlat(2,3) * ry &
       + cfgs(ib)%vec%invlat(3,3) * rz

    ! apply PBC/MIC to fractional coordinates:

    !call pbc_atom_pos_calc(xx, yy, zz)

    if( abs(xx) > 0.5_wp ) &
        xx = xx - anint(xx)

    if( abs(yy) > 0.5_wp ) & 
        yy = yy - anint(yy)

    if( abs(zz) > 0.5_wp ) &
        zz = zz - anint(zz)*real(im_bulk,wp)

    ! get back the Cartesian coordinates (re-apply the cell geomentry)

    rx = cfgs(ib)%vec%latvector(1,1) * xx &
       + cfgs(ib)%vec%latvector(2,1) * yy &
       + cfgs(ib)%vec%latvector(3,1) * zz

    ry = cfgs(ib)%vec%latvector(1,2) * xx &
       + cfgs(ib)%vec%latvector(2,2) * yy &
       + cfgs(ib)%vec%latvector(3,2) * zz

    rz = cfgs(ib)%vec%latvector(1,3) * xx &
       + cfgs(ib)%vec%latvector(2,3) * yy &
       + cfgs(ib)%vec%latvector(3,3) * zz

end subroutine pbc_cart2frac_pos_calc

!> @brief
!> - PBC minimum image convention: puts a \b single \b atom back in simulation cell
!> - first converts the corrdinates: Cartesian -> fractional, accounting for the cell geometry
!> - applies PBC/MIC and converts the coordinates back to Cartesian form
!> @using
!> - `kinds_f90`

!> puts an atom back in simulation cell (PBC/MIC upon Cartesian -> fractional conversion)
subroutine pbc_atom(ib, im, ia)

    !use kinds_f90

    !implicit none

        !> configuration index (microstate or replica)
    integer, intent(in) :: ib

        !> molecule index
    integer, intent(in) :: im

        !> atom index (within molecule?)
    integer, intent(in) :: ia

    real (kind = wp) :: rx, ry, rz, xx, yy, zz

    ! get the current Cartesian coordinates (within the cell geomentry)

    rx = cfgs(ib)%mols(im)%atms(ia)%rpos(1)
    ry = cfgs(ib)%mols(im)%atms(ia)%rpos(2)
    rz = cfgs(ib)%mols(im)%atms(ia)%rpos(3)

    !JG: If orthogonal lattice vectors use simple check
    !AB: make sure that atoms more than half cell dimention away are brought back into the primary cell
    !AB: - "simple check" does not warranty this!

    if(cfgs(ib)%vec%is_orthogonal) then

        !call pbc_ortho_coor_calc(ib, rx, ry, rz)
        call pbc_ortho_coor_full(ib, rx, ry, rz)
        cfgs(ib)%mols(im)%atms(ia)%rpos(1) = rx
        cfgs(ib)%mols(im)%atms(ia)%rpos(2) = ry
        cfgs(ib)%mols(im)%atms(ia)%rpos(3) = rz
   
    else

        ! get the fractional coordinates

        xx = cfgs(ib)%vec%invlat(1,1) * rx + &
             cfgs(ib)%vec%invlat(2,1) * ry + &
             cfgs(ib)%vec%invlat(3,1) * rz

        yy = cfgs(ib)%vec%invlat(1,2) * rx + &
             cfgs(ib)%vec%invlat(2,2) * ry + &
             cfgs(ib)%vec%invlat(3,2) * rz

        zz = cfgs(ib)%vec%invlat(1,3) * rx + &
             cfgs(ib)%vec%invlat(2,3) * ry + &
             cfgs(ib)%vec%invlat(3,3) * rz

        ! apply PBC/MIC to fractional coordinates:

        call pbc_atom_pos_calc(xx, yy, zz)

        ! get back the Cartesian coordinates (re-apply the cell geomentry)

        cfgs(ib)%mols(im)%atms(ia)%rpos(1) = cfgs(ib)%vec%latvector(1,1) * xx + &
                                             cfgs(ib)%vec%latvector(2,1) * yy + &
                                             cfgs(ib)%vec%latvector(3,1) * zz

        cfgs(ib)%mols(im)%atms(ia)%rpos(2) = cfgs(ib)%vec%latvector(1,2) * xx + &
                                             cfgs(ib)%vec%latvector(2,2) * yy + &
                                             cfgs(ib)%vec%latvector(3,2) * zz

        cfgs(ib)%mols(im)%atms(ia)%rpos(3) = cfgs(ib)%vec%latvector(1,3) * xx + & 
                                             cfgs(ib)%vec%latvector(2,3) * yy + &
                                             cfgs(ib)%vec%latvector(3,3) * zz
    end if

end subroutine pbc_atom

!> @brief
!> - minimum image convention with PBC: puts \b all \b atoms \b in a \b molecule back in simulation cell
!> - first converts the corrdinates: Cartesian -> fractional, accounting for the cell geometry
!> - applies PBC/MIC and converts coordinates back to Cartesian form
!> @using
!> - none

!> puts all atoms within a molecule back in simulation cell (PBC/MIC)
subroutine pbc_molecule(ib, im)

    !implicit none

        !> replica configuration/cell & molecule identifiers
    integer, intent(in) :: ib, im
    integer :: i

    do i = 1, cfgs(ib)%mols(im)%natom

        call pbc_atom(ib, im, i)

    enddo

end subroutine pbc_molecule


!> update molecule COM upon accepting an MC move
subroutine update_mol_com_by(ib, im, dcom)

    !use kinds_f90
    !implicit none

        !> replica configuration/cell & molecule identifiers
    integer, intent(in) :: ib, im

        !> COM displacement vector for molecule `im`
    real(kind=wp), intent(in) :: dcom(3)

    cfgs(ib)%mols(im)%rcom(:) = cfgs(ib)%mols(im)%rcom(:) + dcom(:)

end subroutine update_mol_com_by


!> update molecule COM upon accepting an atom move
subroutine update_mol_com(ib, im, ia)

    !use kinds_f90
    !implicit none

        !> replica configuration/cell, molecule & atom identifiers
    integer, intent(in) :: ib, im, ia

    real(kind=wp) :: wght, xshift, yshift, zshift

    wght = cfgs(ib)%mols(im)%atms(ia)%mass / cfgs(ib)%mols(im)%mass

    xshift = wght*(cfgs(ib)%mols(im)%atms(ia)%rpos(1) &
                 - cfgs(ib)%mols(im)%atms(ia)%store_rpos(1))

    yshift = wght*(cfgs(ib)%mols(im)%atms(ia)%rpos(2) &
                 - cfgs(ib)%mols(im)%atms(ia)%store_rpos(2))

    zshift = wght*(cfgs(ib)%mols(im)%atms(ia)%rpos(3) &
                 - cfgs(ib)%mols(im)%atms(ia)%store_rpos(3))

    call add_mol_com(ib, im, xshift, yshift, zshift)

end subroutine update_mol_com


!> add atom displacement to molecule COM upon moving an atom
subroutine add_mol_com(ib, imol, dxc, dyc, dzc)

    !use kinds_f90
    !implicit none

        !> replica configuration/cell & molecule identifiers
    integer, intent(in) :: ib, imol

        !> COM addition coordinates for molecule `im`
    real(kind = wp), intent(in) :: dxc, dyc, dzc

    cfgs(ib)%mols(imol)%rcom(1) = cfgs(ib)%mols(imol)%rcom(1) + dxc
    cfgs(ib)%mols(imol)%rcom(2) = cfgs(ib)%mols(imol)%rcom(2) + dyc
    cfgs(ib)%mols(imol)%rcom(3) = cfgs(ib)%mols(imol)%rcom(3) + dzc

end subroutine add_mol_com


!> check if atom is still inside the slit
!> i.e. in between the bottom & top walls
logical function atom_inside(zc, zl)

    use kinds_f90
    use slit_module, only : in_slit, slit_zfrac1, slit_zfrac2

    implicit none

        !> atom Cartesian Z coordinate
    real (kind = wp), intent(in) :: zc

        !> cell Cartesian Z dimension: latvector(3,3)
        !> applicable only with Cartesian cell vectors: cfgfmt = 1
    real (kind = wp), intent(in) :: zl

    atom_inside = .true.

    if ( in_slit ) atom_inside = ( zc > zl*slit_zfrac1 .and. zc < zl*slit_zfrac2 )
       !atom_inside = zc <= zl*0.5_wp .and. zc >= -zl*0.5_wp

end function atom_inside


!> check if atom got outside the slit
logical function atom_outside(zc, zl)

    use kinds_f90
    use slit_module, only : in_slit, slit_zfrac1, slit_zfrac2

    implicit none

        !> atom Cartesian Z coordinate
    real (kind = wp), intent(in) :: zc

        !> cell Cartesian Z dimension: latvector(3,3)
        !> applicable only for cubic cells: cfgfmt = 1
    real (kind = wp), intent(in) :: zl

    atom_outside = .false.

    if( in_slit ) atom_outside = ( zc <= zl*slit_zfrac1 .or. zc >= zl*slit_zfrac2 )
       !atom_outside = zc > zl*0.5_wp .and. zc < -zl*0.5_wp

end function atom_outside


!> translate molecule to new COM and check if any of its atoms got outside the slit
logical function position_molecule_out(mol, pos, zdim)

    use molecule_type

    implicit none

        !> new COM vector for molecule `mol`
    real (kind = wp), intent(in) :: pos(3)

        !> Z dimension of cell
    real (kind = wp), intent(in) :: zdim

        !> molecule object
    type(molecule), intent(inout) :: mol

    integer :: j

    position_molecule_out = .false.

    do j = 1, mol%natom

        !remove center of mass
        !mol%atms(j)%rpos(:) = mol%atms(j)%rpos(:) - mol%rcom(:)
        !mol%atms(j)%rpos(:) = mol%atms(j)%rpos(:) + pos(:)

        mol%atms(j)%rpos(:) = mol%atms(j)%rpos(:) - mol%rcom(:) + pos(:)

        position_molecule_out = ( position_molecule_out .or. &
                                  atom_outside(mol%atms(j)%rpos(3), zdim) )

        !if( position_molecule_out ) return

    enddo

    mol%rcom(:) = pos(:)

end function position_molecule_out


!> move atom `atom` in molecule `imol` by vector `rmov` (in replica `ib`)
subroutine move_atom_by(ib, imol, atom, rmov, atom_out)

#ifdef DEBUG
!debugging!
    use slit_module, only : in_bulk
#endif

        !> replica configuration/cell, molecule & atom identifiers
    integer, intent(in) :: ib, imol, atom

        !> displacement vector for atom `atom` in molecule `imol`
    real(kind = wp), intent(in) :: rmov(3)

        !> flag for any atom being outside (if in slit)
    logical, intent(inout) :: atom_out

    cfgs(ib)%mols(imol)%atms(atom)%rpos(:) = cfgs(ib)%mols(imol)%atms(atom)%rpos(:) &
                                           + rmov(:)

    atom_out = atom_outside(cfgs(ib)%mols(imol)%atms(atom)%rpos(3), &
                            cfgs(ib)%vec%latvector(3,3))

    if( cfgs(ib)%vec%is_orthogonal ) &
        call pbc_ortho_coor_calc(ib,  &
             cfgs(ib)%mols(imol)%atms(atom)%rpos(1), &
             cfgs(ib)%mols(imol)%atms(atom)%rpos(2), &
             cfgs(ib)%mols(imol)%atms(atom)%rpos(3))

#ifdef DEBUG
!debugging!
    if(in_bulk .and. atom_out) write(*,*)'cell_module::move_atom(..) - ERROR!!! undue atom_out = ', atom_out
#endif

end subroutine move_atom_by


!> random move for atom `atom` in molecule `imol` by max `distmax` (in replica `ib`)
subroutine move_atom(ib, imol, atom, distmax, atom_out)

    !use kinds_f90
    use random_module, only : duni

#ifdef DEBUG
!debugging!
    use slit_module, only : in_bulk
#endif

        !> replica configuration/cell, molecule & atom identifiers
    integer, intent(in) :: ib, imol, atom

        !> max displacement (scalar) for atom `atom` in molecule `imol`
    real(kind = wp), intent(in) :: distmax

        !> flag for any atom being outside (if in slit)
    logical, intent(inout) :: atom_out

    cfgs(ib)%mols(imol)%atms(atom)%rpos(1) = cfgs(ib)%mols(imol)%atms(atom)%rpos(1) &
                                           + 2.0_wp * (duni() - 0.5_wp) * distmax

    cfgs(ib)%mols(imol)%atms(atom)%rpos(2) = cfgs(ib)%mols(imol)%atms(atom)%rpos(2) &
                                           + 2.0_wp * (duni() - 0.5_wp) * distmax

    cfgs(ib)%mols(imol)%atms(atom)%rpos(3) = cfgs(ib)%mols(imol)%atms(atom)%rpos(3) &
                                           + 2.0_wp * (duni() - 0.5_wp) * distmax

    atom_out = atom_outside(cfgs(ib)%mols(imol)%atms(atom)%rpos(3), &
                            cfgs(ib)%vec%latvector(3,3))

    if( cfgs(ib)%vec%is_orthogonal ) &
        call pbc_ortho_coor_calc(ib,  &
             cfgs(ib)%mols(imol)%atms(atom)%rpos(1), &
             cfgs(ib)%mols(imol)%atms(atom)%rpos(2), &
             cfgs(ib)%mols(imol)%atms(atom)%rpos(3))

#ifdef DEBUG
!debugging!
    if(in_bulk .and. atom_out) write(*,*)'cell_module::move_atom(..) - ERROR!!! undue atom_out = ', atom_out
#endif

end subroutine move_atom


!> translate molecule `im` by a given vector `dcom`
subroutine move_molecule_by(ib, im, dcom, atom_out)

    !use kinds_f90
    use molecule_module
    use random_module, only : duni

#ifdef DEBUG
!debugging!
    use slit_module, only : in_bulk
#endif

    !implicit none

        !> replica configuration/cell & molecule identifiers
    integer, intent(in) :: ib, im

        !> displacement vector for COM of molecule `im`
    real(kind = wp), intent(in) :: dcom(3)

        !> flag for any atom being outside (if in slit)
    logical, intent(inout) :: atom_out

    integer :: j

    ! move atoms
    do j = 1, cfgs(ib)%mols(im)%natom

        ! first, update and check Z coordinate in case atom gets out of the cell/slit
        cfgs(ib)%mols(im)%atms(j)%rpos(3) = cfgs(ib)%mols(im)%atms(j)%rpos(3) + dcom(3)

        atom_out = atom_outside(cfgs(ib)%mols(im)%atms(j)%rpos(3), cfgs(ib)%vec%latvector(3,3))

#ifdef DEBUG
!debugging!
        if(in_bulk .and. atom_out) write(*,*)'cell_module::move_molecule_by(..) ERROR!!! undue atom_out = ', atom_out
#endif

        if( atom_out ) return

        ! if still inside, proceed to update X & Y
        cfgs(ib)%mols(im)%atms(j)%rpos(1) = cfgs(ib)%mols(im)%atms(j)%rpos(1) + dcom(1)
        cfgs(ib)%mols(im)%atms(j)%rpos(2) = cfgs(ib)%mols(im)%atms(j)%rpos(2) + dcom(2)

        if( cfgs(ib)%vec%is_orthogonal ) &
            call pbc_ortho_coor_calc(ib,  &
                 cfgs(ib)%mols(im)%atms(j)%rpos(1), &
                 cfgs(ib)%mols(im)%atms(j)%rpos(2), &
                 cfgs(ib)%mols(im)%atms(j)%rpos(3))

    enddo

    ! move com 
    cfgs(ib)%mols(im)%rcom(:) = cfgs(ib)%mols(im)%rcom(:) + dcom(:)

!AB: putting COM:s back into the primary cell only at molecule moves/rotations 
!AB: is at odds with atom moves where this is not done (or would add extra overheads)
!AB: this also defeats tracking COM's drift upon acceptance (done in all MC moves)
!AB: there is no need to keep COM:s inside the primary cell!

    !if( cfgs(ib)%vec%is_orthogonal ) &
    !    call pbc_ortho_coor_calc(ib, &
    !         cfgs(ib)%mols(im)%rcom(1), & 
    !         cfgs(ib)%mols(im)%rcom(2), &
    !         cfgs(ib)%mols(im)%rcom(3) )

    !AB: not needed here, provided COM:s have been set correctly beforehand
    ! calculate center of mass
    !call mol_com(cfgs(ib)%mols(im))

    ! check within box (done on the system globally, due to "use reset" in CONTROL)
    !call pbc_molecule(ib, im)

end subroutine move_molecule_by


!> random translation of molecule `im` by distance `distmax` at max (in replica `ib`)
subroutine move_molecule(ib, im, distmax, atom_out)

    !use kinds_f90
    use molecule_module
    use random_module, only : duni

#ifdef DEBUG
!debugging!
    use slit_module, only : in_bulk
#endif

    !implicit none

        !> replica configuration/cell & molecule identifiers
    integer, intent(in) :: ib, im

        !> max displacement (scalar) for molecule `im`
    real(kind = wp), intent(in) :: distmax

        !> flag for any atom being outside (if in slit)
    logical, intent(inout) :: atom_out

    integer :: j
    real(kind = wp) :: xshift, yshift, zshift

    ! displacement
    xshift = 2.0_wp * (duni() - 0.5_wp) * distmax
    yshift = 2.0_wp * (duni() - 0.5_wp) * distmax
    zshift = 2.0_wp * (duni() - 0.5_wp) * distmax

    ! store and move atoms
    do j = 1, cfgs(ib)%mols(im)%natom

        ! first, update and check Z coordinate in case atom gets out of the cell/slit
        cfgs(ib)%mols(im)%atms(j)%rpos(3) = cfgs(ib)%mols(im)%atms(j)%rpos(3) + zshift

        atom_out = atom_outside(cfgs(ib)%mols(im)%atms(j)%rpos(3), cfgs(ib)%vec%latvector(3,3))

#ifdef DEBUG
!debugging!
        if(in_bulk .and. atom_out) write(*,*)'cell_module::move_molecule(..) ERROR!!! undue atom_out = ', atom_out
#endif

        if( atom_out ) return

        ! if still inside, proceed to update X & Y
        cfgs(ib)%mols(im)%atms(j)%rpos(1) = cfgs(ib)%mols(im)%atms(j)%rpos(1) + xshift
        cfgs(ib)%mols(im)%atms(j)%rpos(2) = cfgs(ib)%mols(im)%atms(j)%rpos(2) + yshift

        if( cfgs(ib)%vec%is_orthogonal ) &
            call pbc_ortho_coor_calc(ib,  &
                 cfgs(ib)%mols(im)%atms(j)%rpos(1), &
                 cfgs(ib)%mols(im)%atms(j)%rpos(2), &
                 cfgs(ib)%mols(im)%atms(j)%rpos(3))

    end do

    ! move com 
    cfgs(ib)%mols(im)%rcom(1) = cfgs(ib)%mols(im)%rcom(1) + xshift
    cfgs(ib)%mols(im)%rcom(2) = cfgs(ib)%mols(im)%rcom(2) + yshift
    cfgs(ib)%mols(im)%rcom(3) = cfgs(ib)%mols(im)%rcom(3) + zshift

!AB: putting COM:s back into the primary cell only at molecule moves/rotations 
!AB: is at odds with atom moves where this is not done (or would add extra overheads)
!AB: this also defeats tracking COM's drift upon acceptance (done in all MC moves)
!AB: there is no need to keep COM:s inside the primary cell!

!    if( cfgs(ib)%vec%is_orthogonal ) &
!        call pbc_ortho_coor_calc(ib,  &
!             cfgs(ib)%mols(im)%rcom(1), & 
!             cfgs(ib)%mols(im)%rcom(2), &
!             cfgs(ib)%mols(im)%rcom(3) )

    !AB: not needed here, provided COM:s have been set correctly beforehand
    ! calculate center of mass
    !call mol_com(cfgs(ib)%mols(im))

    ! check within box (done on the system globally, due to "use reset" in CONTROL)
    !call pbc_molecule(ib, im)

end subroutine move_molecule


!> Generate a random rotation matrix using the quaternion method
!> The matrix is returned as an array of dimension 9
function generate_rot_matrix_quat(rotmax)

    use random_module, only : duni

        !> max rotation angle (scalar)
    real(kind = wp), intent(in) :: rotmax

    real(kind=wp), dimension(9) :: generate_rot_matrix_quat
    
    real(kind = wp) :: angle, axis(3), qtn(4), raxis

    ! generate a random rotation axis
    !JG: axis generation is fine but needs to be rejected if it lies outside the unit sphere
    ! otherwise the corners of the sphere cause systematic oversampling of diagonal axis

    axis(1) = 2.0_wp * duni() - 1.0_wp
    axis(2) = 2.0_wp * duni() - 1.0_wp
    axis(3) = 2.0_wp * duni() - 1.0_wp

    !JG: first calculate the length of the axis vector
    raxis = dot_product(axis,axis)

    !JG: now since we are considering the unit do not yet need to take the sqaure root
    do while (raxis > 1.0_wp)
        !JG: generate new axis

        axis(1) = 2.0_wp * duni() - 1.0_wp
        axis(2) = 2.0_wp * duni() - 1.0_wp
        axis(3) = 2.0_wp * duni() - 1.0_wp

        !JG: calculate the length of the new axis vector
        raxis = dot_product(axis,axis)
  
    end do

    !JG now normalise axis using square root of already calculated raxis
    axis(:) = axis(:)/sqrt(raxis)

    ! generate a random rotation angle
    angle   = rotmax*(2.0_wp*duni()-1.0_wp)
    qtn(1)   = cos(0.5_wp * angle)
    qtn(2:4) = axis(:) * sin(0.5_wp * angle)

    ! normalise
    !JG: need check but the below should be unecessary since the quaternion is normalised by construction
    qtn(:) = qtn(:)/sqrt(dot_product(qtn,qtn))

    ! construct rotation matrix
    generate_rot_matrix_quat(1) = qtn(1)**2+qtn(2)**2-qtn(3)**2-qtn(4)**2
    generate_rot_matrix_quat(2) = 2.0_wp*(qtn(2)*qtn(3) - qtn(1)*qtn(4))
    generate_rot_matrix_quat(3) = 2.0_wp*(qtn(2)*qtn(4) + qtn(1)*qtn(3))
    generate_rot_matrix_quat(4) = 2.0_wp*(qtn(2)*qtn(3) + qtn(1)*qtn(4))
    generate_rot_matrix_quat(5) = qtn(1)**2-qtn(2)**2+qtn(3)**2-qtn(4)**2
    generate_rot_matrix_quat(6) = 2.0_wp*(qtn(3)*qtn(4) - qtn(1)*qtn(2))
    generate_rot_matrix_quat(7) = 2.0_wp*(qtn(2)*qtn(4) - qtn(1)*qtn(3))
    generate_rot_matrix_quat(8) = 2.0_wp*(qtn(3)*qtn(4) + qtn(1)*qtn(2))
    generate_rot_matrix_quat(9) = qtn(1)**2-qtn(2)**2-qtn(3)**2+qtn(4)**2
  
end function generate_rot_matrix_quat   


!> random rotation of molecule `im` by angle `rotmax` at max (in replica `ib`)
!> (using quaternion approach). Also correspondingly rotates the orientations of 
!> 'spin' type atoms in the molecule
subroutine rotate_molecule_quat(ib, im, rotmax, atom_out, rot)

    !use kinds_f90
    use molecule_module
    use constants_module

#ifdef DEBUG
!debugging!
    use slit_module, only : in_bulk
#endif

    !implicit none

        !> replica configuration/cell & molecule identifiers
    integer, intent(in) :: ib, im

        !> max rotation angle (scalar) for molecule `im`
    real(kind = wp), intent(in) :: rotmax

        !> flag for any atom being outside (if in slit)
    logical, intent(inout) :: atom_out

      !> Rotation matrix applied to the molecule rotated by this procedure
    real(kind = wp), intent(out), optional :: rot(9)

    integer :: j, k

    real(kind = wp) :: angle, axis(3), qtn(4), raxis

    rot = generate_rot_matrix_quat(rotmax) 

    call set_mol_com(ib,im)

    ! subtract centre of mass
    do j = 1, cfgs(ib)%mols(im)%natom

        cfgs(ib)%mols(im)%atms(j)%rpos(:) = cfgs(ib)%mols(im)%atms(j)%rpos(:) - cfgs(ib)%mols(im)%rcom(:)

    enddo

    call molecule_quaternion_rotation(cfgs(ib)%mols(im), rot)

    ! add back on centre of mass
    do j = 1, cfgs(ib)%mols(im)%natom
    
        cfgs(ib)%mols(im)%atms(j)%rpos(:) = cfgs(ib)%mols(im)%atms(j)%rpos(:) + cfgs(ib)%mols(im)%rcom(:)

        atom_out = atom_outside(cfgs(ib)%mols(im)%atms(j)%rpos(3), cfgs(ib)%vec%latvector(3,3))

#ifdef DEBUG
!debugging!
        if(in_bulk .and. atom_out) write(*,*)'cell_module::rotate_molecule_quat(..) ERROR!!! undue atom_out = ', atom_out
#endif

        if( atom_out ) return

    enddo

    if(cfgs(ib)%vec%is_orthogonal) then

        do j = 1, cfgs(ib)%mols(im)%natom
            call pbc_ortho_coor_calc(ib, &
                 cfgs(ib)%mols(im)%atms(j)%rpos(1), &
                 cfgs(ib)%mols(im)%atms(j)%rpos(2), &
                 cfgs(ib)%mols(im)%atms(j)%rpos(3) )
        end do

    end if

end subroutine rotate_molecule_quat


!TU: Added by me...
!> Rotates a specified molecule by a specified rotation matrix.  
!> Also correspondingly rotates the orientations of 'spin' type atoms in the molecule
subroutine rotate_molecule_quat_by(ib, im, atom_out, rot)

    !use kinds_f90
    use molecule_module
    use constants_module
    use random_module, only : duni

#ifdef DEBUG
!debugging!
    use slit_module, only : in_bulk
#endif

    !implicit none

        !> replica configuration/cell & molecule identifiers
    integer, intent(in) :: ib, im

        !> flag for any atom being outside (if in slit)
    logical, intent(inout) :: atom_out

        !> Rotation matrix applied to the molecule rotated by this procedure
    real(kind = wp), intent(in) :: rot(9)

    integer :: j, k

    real(kind = wp) :: angle, axis(3), qtn(4), raxis

    call set_mol_com(ib,im)

    ! subtract centre of mass
    do j = 1, cfgs(ib)%mols(im)%natom

        cfgs(ib)%mols(im)%atms(j)%rpos(:) = cfgs(ib)%mols(im)%atms(j)%rpos(:) - cfgs(ib)%mols(im)%rcom(:)

    enddo

    call molecule_quaternion_rotation(cfgs(ib)%mols(im), rot)

    ! add back on centre of mass
    do j = 1, cfgs(ib)%mols(im)%natom
    
        cfgs(ib)%mols(im)%atms(j)%rpos(:) = cfgs(ib)%mols(im)%atms(j)%rpos(:) + cfgs(ib)%mols(im)%rcom(:)

        atom_out = atom_outside(cfgs(ib)%mols(im)%atms(j)%rpos(3), cfgs(ib)%vec%latvector(3,3))

#ifdef DEBUG
!debugging!
        if(in_bulk .and. atom_out) write(*,*)'cell_module::rotate_molecule_quat(..) ERROR!!! undue atom_out = ', atom_out
#endif

        if( atom_out ) return

    enddo

    if(cfgs(ib)%vec%is_orthogonal) then

        do j = 1, cfgs(ib)%mols(im)%natom
            call pbc_ortho_coor_calc(ib, &
                 cfgs(ib)%mols(im)%atms(j)%rpos(1), &
                 cfgs(ib)%mols(im)%atms(j)%rpos(2), &
                 cfgs(ib)%mols(im)%atms(j)%rpos(3) )
        end do

!AB: putting COM:s back into the primary cell only at molecule moves/rotations 
!AB: is at odds with atom moves where this is not done (or would add extra overheads)
!AB: this also defeats tracking COM's drift upon acceptance (done in all MC moves)
!AB: there is no need to keep COM:s inside the primary cell!

        !call pbc_ortho_coor_calc(ib, &
        !     cfgs(ib)%mols(im)%rcom(1), &
        !     cfgs(ib)%mols(im)%rcom(2), &
        !     cfgs(ib)%mols(im)%rcom(3) )
    end if

end subroutine rotate_molecule_quat_by


!> rotation of molecule `im` by angle `rotmax` at max, around OX/OY/OZ (in replica `ib`)
!> (using Euler's approach).
!> This procedure does not rotate atomic orientations
!> in the molecule (spins), and hence is unsuitable for molecules containing
!> 'spin' type atoms
subroutine rotate_molecule(ib, im, rotmax, atom_out)

    !use kinds_f90
    use molecule_module
    use constants_module
    use random_module, only : duni

#ifdef DEBUG
!debugging!
    use slit_module, only : in_bulk
#endif

    !implicit none

        !> replica configuration/cell & molecule identifiers
    integer, intent(in) :: ib, im

        !> max rotation angle (scalar) for molecule `im`
    real(kind = wp), intent(in) :: rotmax

        !> flag for any atom being outside (if in slit)
    logical, intent(inout) :: atom_out

    integer :: j, k, axis

    real(kind = wp) :: tmpx, tmpy, tmpz, dr, cosdr, sindr

    call set_mol_com(ib,im)

    ! remove center of mass
    do j = 1, cfgs(ib)%mols(im)%natom

        cfgs(ib)%mols(im)%atms(j)%rpos(:) = cfgs(ib)%mols(im)%atms(j)%rpos(:) - cfgs(ib)%mols(im)%rcom(:)

    enddo

    do axis = 1, 3

        ! displacement
        dr = 2.0_wp * (duni() - 0.5_wp) * rotmax

        select case (axis)

            case (1)    !     rotate about x axis

                 cosdr = cos(dr)
                 sindr = sin(dr)

                 do j = 1, cfgs(ib)%mols(im)%natom

                     ! move to origin
                     tmpy = cfgs(ib)%mols(im)%atms(j)%rpos(2)
                     tmpz = cfgs(ib)%mols(im)%atms(j)%rpos(3)

                     ! rotate molecule + shift back from origin
                     cfgs(ib)%mols(im)%atms(j)%rpos(2) = cosdr * tmpy + sindr * tmpz
                     cfgs(ib)%mols(im)%atms(j)%rpos(3) = cosdr * tmpz - sindr * tmpy

                 enddo

            case (2)    !     rotate about y axis

                 cosdr = cos(dr)
                 sindr = sin(dr)

                 do j = 1, cfgs(ib)%mols(im)%natom

                     tmpx = cfgs(ib)%mols(im)%atms(j)%rpos(1)
                     tmpz = cfgs(ib)%mols(im)%atms(j)%rpos(3)

                     ! rotate molecule
                     cfgs(ib)%mols(im)%atms(j)%rpos(1) = cosdr * tmpx - sindr * tmpz
                     cfgs(ib)%mols(im)%atms(j)%rpos(3) = cosdr * tmpz + sindr * tmpx

                 enddo

            case (3)    !    rotate about z

                 cosdr = cos(dr)
                 sindr = sin(dr)

                 do j = 1, cfgs(ib)%mols(im)%natom

                     tmpx = cfgs(ib)%mols(im)%atms(j)%rpos(1)
                     tmpy = cfgs(ib)%mols(im)%atms(j)%rpos(2)

                     ! rotate molecule
                     cfgs(ib)%mols(im)%atms(j)%rpos(1) = cosdr * tmpx + sindr * tmpy
                     cfgs(ib)%mols(im)%atms(j)%rpos(2) = cosdr * tmpy - sindr * tmpx

                 enddo

        end select

    enddo

    ! add center of mass back
    do j = 1, cfgs(ib)%mols(im)%natom

        cfgs(ib)%mols(im)%atms(j)%rpos(:) = cfgs(ib)%mols(im)%atms(j)%rpos(:) + cfgs(ib)%mols(im)%rcom(:)

        atom_out = atom_outside(cfgs(ib)%mols(im)%atms(j)%rpos(3), cfgs(ib)%vec%latvector(3,3))

#ifdef DEBUG
!debugging!
        if(in_bulk .and. atom_out) write(*,*)'cell_module::rotate_molecule(..) ERROR!!! undue atom_out = ', atom_out
#endif

        if( atom_out ) return

    enddo

    if(cfgs(ib)%vec%is_orthogonal) then

        do j = 1, cfgs(ib)%mols(im)%natom
            call pbc_ortho_coor_calc(ib, &
                 cfgs(ib)%mols(im)%atms(j)%rpos(1), &
                 cfgs(ib)%mols(im)%atms(j)%rpos(2), &
                 cfgs(ib)%mols(im)%atms(j)%rpos(3) )
        end do

    end if

    !call mol_com(cfgs(ib)%mols(im))

!   do j = 1, cfgs(ib)%mols(im)%natom

!       ! check within box
!       call pbc_atom(ib, im, j)

!   enddo

end subroutine rotate_molecule


!> back up the molecule (all its atoms) position(s) before an MC move
subroutine store_molecule_pos(ib, im)

    !implicit none

        !> replica configuration/cell & molecule identifiers
    integer, intent(in) :: ib, im

    integer j, i

    ! store com

    cfgs(ib)%mols(im)%store_rcom(:) = cfgs(ib)%mols(im)%rcom(:)

    ! store atoms
    do j = 1, cfgs(ib)%mols(im)%natom

        !store positions first
        cfgs(ib)%mols(im)%atms(j)%store_rpos(:) = cfgs(ib)%mols(im)%atms(j)%rpos(:)

    enddo

end subroutine store_molecule_pos


!>  reset the atom position within a molecule and correct the molecule COM
subroutine restore_mol_atom_pos(ib, im, ia)

        !> replica configuration/cell, molecule & atom identifiers
    integer, intent(in) :: ib, im, ia

    real(kind=wp) :: xshift, yshift, zshift, wght

    xshift = cfgs(ib)%mols(im)%atms(ia)%store_rpos(1) - cfgs(ib)%mols(im)%atms(ia)%rpos(1)
    yshift = cfgs(ib)%mols(im)%atms(ia)%store_rpos(2) - cfgs(ib)%mols(im)%atms(ia)%rpos(2)
    zshift = cfgs(ib)%mols(im)%atms(ia)%store_rpos(3) - cfgs(ib)%mols(im)%atms(ia)%rpos(3)

    wght = cfgs(ib)%mols(im)%atms(ia)%mass / cfgs(ib)%mols(im)%mass

    call add_mol_com(ib,im,xshift*wght,yshift*wght,zshift*wght)

    cfgs(ib)%mols(im)%atms(ia)%rpos(:) = cfgs(ib)%mols(im)%atms(ia)%store_rpos(:)

end subroutine restore_mol_atom_pos


!>  reset the molecule position and COM
subroutine restore_molecule_pos(ib, im)

        !> replica configuration/cell & molecule identifiers
    integer, intent(in) :: ib, im

    integer :: j

    cfgs(ib)%mols(im)%rcom(:) = cfgs(ib)%mols(im)%store_rcom(:)

    do j = 1, cfgs(ib)%mols(im)%natom

       cfgs(ib)%mols(im)%atms(j)%rpos(:) = cfgs(ib)%mols(im)%atms(j)%store_rpos(:)

    enddo

end subroutine


!> stores all the atoms 
subroutine store_all_pos(ib)

    !implicit none

        !> replica configuration/cell identifier
    integer, intent(in) :: ib

    integer :: i, j

    do i = 1, cfgs(ib)%num_mols

        cfgs(ib)%mols(i)%store_rcom(:) = cfgs(ib)%mols(i)%rcom(:)

        do j = 1, cfgs(ib)%mols(i)%natom

            cfgs(ib)%mols(i)%atms(j)%store_rpos(:) = cfgs(ib)%mols(i)%atms(j)%rpos(:)

        enddo

    enddo

end subroutine store_all_pos


!>  resets all the atoms 
subroutine restore_all_pos(ib)

    !implicit none

        !> replica configuration/cell identifier
    integer, intent(in) :: ib
    integer :: i, j

    do i = 1, cfgs(ib)%num_mols

        cfgs(ib)%mols(i)%rcom(:) = cfgs(ib)%mols(i)%store_rcom(:)

        do j = 1, cfgs(ib)%mols(i)%natom

            cfgs(ib)%mols(i)%atms(j)%rpos(:) = cfgs(ib)%mols(i)%atms(j)%store_rpos(:)

        enddo

    enddo

end subroutine restore_all_pos


!> converts fractional to Cartesian coordinates
subroutine fractocart(ib)

    !use kinds_f90
    use latticevectors_module, only : frac_to_cart

    !implicit none

        !> replica configuration/cell identifier
    integer, intent(in) :: ib

    integer :: i, j

    real(kind = wp) :: rx, ry, rz

    call frac_to_cart(cfgs(ib))

end subroutine fractocart


!> inverts the cell matrix (lattice vectors)
subroutine invert_latticevectors(ib)

    !use kinds_f90
    use latticevectors_module

        !> replica configuration/cell identifier
    integer, intent(in) :: ib

    call invertlatvec(cfgs(ib)%vec)

end subroutine invert_latticevectors


!> converts Cartesian to fractional coordinates
subroutine carttofrac(ib)

    !use kinds_f90
    use latticevectors_module, only : cart_to_frac

        !> replica configuration/cell identifier
    integer, intent(in) :: ib

    integer :: i,j

    real(kind = wp) :: rx, ry, rz

    ! should no longer be necessary as a current version of invlat
    ! is held in memory
    ! call invertlatvec(cfgs(ib)%vec)

    call cart_to_frac(cfgs(ib))

    return

    do j = 1, cfgs(ib)%num_mols

        do i = 1, cfgs(ib)%mols(j)%natom

            rx = cfgs(ib)%mols(j)%atms(i)%rpos(1)
            ry = cfgs(ib)%mols(j)%atms(i)%rpos(2)
            rz = cfgs(ib)%mols(j)%atms(i)%rpos(3)

            cfgs(ib)%mols(j)%atms(i)%rpos(1) = cfgs(ib)%vec%invlat(1,1) * rx + cfgs(ib)%vec%invlat(2,1) * ry +    &
                                         cfgs(ib)%vec%invlat(3,1) * rz
            cfgs(ib)%mols(j)%atms(i)%rpos(2) = cfgs(ib)%vec%invlat(1,2) * rx + cfgs(ib)%vec%invlat(2,2) * ry +    &
                                         cfgs(ib)%vec%invlat(3,2) * rz
            cfgs(ib)%mols(j)%atms(i)%rpos(3) = cfgs(ib)%vec%invlat(1,3) * rx + cfgs(ib)%vec%invlat(2,3) * ry +    &
                                         cfgs(ib)%vec%invlat(3,3) * rz

        enddo

    enddo

end subroutine carttofrac


!> finds the number of moving sites for volume move
subroutine setup_volume_moves(job)

    !use kinds_f90
    use control_type
    use constants_module, only : uout
    use species_module, only : number_of_molecules, uniq_mol, eletype, element, number_of_elements
    use comms_mpi_module, only : master
    !use parallel_loop_module, only : master

    !implicit none

    type(control), intent(inout) :: job
    integer :: j, k, ib, movers, num

    do ib = 1, nconfigs

        movers = 0

        if (job%moveatm) then

            do j = 1, number_of_elements

                do k = 1, job%numatmmove

                    if(element(j) == job%atmmover(k) .and. eletype(j) == job%atmmover_type(k)) then

                        call findnumtype(ib, j, num)

                        movers = movers + num

                    endif

                enddo

            enddo

        endif

        if (job%movemol) then

            ! set flag for moveable molecules
            do j = 1, number_of_molecules

                do k = 1, job % nummolmove

                    if (uniq_mol(j)%molname == job % molmover(k)) then

                        num = cfgs(ib)%mtypes(j)%num_mols

                        movers = movers + num

                    endif

                enddo

             enddo

        endif

        cfgs(ib)%num_vol_sites = movers

        if (master) write(uout,"(/,1x,'number of moving particles for box ',i4,i10)")ib, cfgs(ib)%num_vol_sites

    enddo

end subroutine setup_volume_moves


!TU-: This is obsolete and could be deleted during a refactor: expand_cell_ortho
!TU-: could be used to do the same thing as expand_cell

!> isotropic expansion/contraction during V-step in NPT MC for cubic cells
subroutine expand_cell(ib, movetype, maxvolchange, volnew, scaling)

    !use kinds_f90
    use latticevectors_module, only : invertlatvec, expandcell
    !use slit_module,           only : in_slit

    !implicit none

        !> id/index of the configuration box/cell to work on
    integer, intent(in) :: ib

        !> type of V-step to attempt: ln(V), 1/L, L (see `constants_module.f90`)
    integer, intent(in) :: movetype

        !> current max V-step (displacement)
    real (kind = wp), intent(in) :: maxvolchange

        !> resulting attempted volume
    real (kind = wp), intent(out) :: volnew

        !> resulting scaling factor
    real (kind = wp), intent(out) :: scaling

    ! back-up all atom postions in the cell
    call store_all_pos(ib)

    !JG: Currently rigid molecule non-scaling is performed in scale_positions
    ! However when PIC is implemented, as in orthogonal flag
    ! the old lattice vectors should be used to put a molecule in its own frame
    ! not the new ones as occurs in current implementation.
    ! Solution is to set molecule into its own frame, before expanding cell.
    ! Future implementation will avoid this by using relative position of
    ! rigid molecules to generate new configuration but for the time being below is required:

!    if(cfgs(ib)%vec%is_orthogonal) then
    call set_rigid_mol_coms(ib) 
!    end if

    ! attempt the V-step
    !call expandcell(scaling, movetype, maxvolchange, volnew, cfgs(ib)%vec)
    !call expandcell(cfgs(ib)%vec, scaling, volnew, maxvolchange, movetype)
    call expandcell(cfgs(ib)%vec, scaling, volnew, maxvolchange, movetype)

    ! get the attempted inverted cell matrix
    call invertlatvec(cfgs(ib)%vec)

    ! get the new atom positions by rescaling
    call scale_positions(ib, scaling)

end subroutine expand_cell


!TU: I added this. It is almost the same as 'expand_cell'
!> Shape-preserving cell expansion/contraction during V-step in NPT MC for ORTHORHOMBIC cells
subroutine expand_cell_ortho(ib, movetype, maxvolchange, volnew, scaling)

    !use kinds_f90
    use latticevectors_module, only : invertlatvec, expandcell_ortho
    !use slit_module,           only : in_slit

    !implicit none

        !> id/index of the configuration box/cell to work on
    integer, intent(in) :: ib

        !> type of V-step to attempt: ln(V), 1/L, L (see `constants_module.f90`)
    integer, intent(in) :: movetype

        !> current max V-step (displacement)
    real (kind = wp), intent(in) :: maxvolchange

        !> resulting attempted volume
    real (kind = wp), intent(out) :: volnew

        !> resulting scaling factor for each dimension of the cell
    real (kind = wp), intent(out) :: scaling

    ! back-up all atom postions in the cell
    call store_all_pos(ib)

!    if(cfgs(ib)%vec%is_orthogonal) then
    call set_rigid_mol_coms(ib) 
!    end if

    ! attempt the V-step
    call expandcell_ortho(cfgs(ib)%vec, scaling, volnew, maxvolchange, movetype, myjob%is_XY_npt)

    ! get the attempted inverted cell matrix
    call invertlatvec(cfgs(ib)%vec)

    ! get the new atom positions by rescaling
    call scale_positions(ib, scaling)

end subroutine expand_cell_ortho




!> Anisotropic cell expansion/contraction during V-step in NPT MC for ORTHORHOMBIC cells
subroutine expand_cell_orthoani(ib, indx, maxvolchange, volnew, scaling)

    use latticevectors_module, only : invertlatvec, expandcell_orthoani

        !> id/index of the configuration box/cell to work on
    integer, intent(in) :: ib

        !> dimension of the box to expand/contract (1=x,2=y,3=z)
    integer, intent(in) :: indx
    
        !> current max V-step (displacement) for the dimension under consideration
    real (kind = wp), intent(in) :: maxvolchange

        !> resulting attempted volume
    real (kind = wp), intent(out) :: volnew

        !> resulting scaling factor applied to the dimension of the cell
    real (kind = wp), intent(out) :: scaling

    ! back-up all atom postions in the cell
    call store_all_pos(ib)

    call set_rigid_mol_coms(ib) 

    ! attempt the V-step
    call expandcell_orthoani(cfgs(ib)%vec, indx, scaling, volnew, maxvolchange)

    ! get the attempted inverted cell matrix
    call invertlatvec(cfgs(ib)%vec)

    ! get the new atom positions by rescaling
    call scale_positions_orthoani(ib, indx, scaling)

end subroutine expand_cell_orthoani




!TU: I added this...
!> Shape-preserving cell expansion/contraction V-step (NPT for ORTHORHOMBIC cells; possible with slit XY PBC),
!> where the scaling applied to each cell dimension is SPECIFIED (i.e. not random)
subroutine expand_cell_ortho_by(ib, volnew, scaling)

    !use kinds_f90
    use latticevectors_module, only : invertlatvec, expandcell_ortho_by
    !use slit_module,           only : in_slit

    !implicit none

        !> id/index of the configuration box/cell to work on
    integer, intent(in) :: ib

        !> resulting attempted volume
    real (kind = wp), intent(out) :: volnew

        !> scaling factor to be applied to each dimension of the cell
    real (kind = wp), intent(in) :: scaling

    ! back-up all atom postions in the cell
    call store_all_pos(ib)

!    if(cfgs(ib)%vec%is_orthogonal) then
    call set_rigid_mol_coms(ib) 
!    end if

    ! attempt the V-step
    call expandcell_ortho_by(cfgs(ib)%vec, scaling, volnew, myjob%is_XY_npt)

    ! get the attempted inverted cell matrix
    call invertlatvec(cfgs(ib)%vec)

    ! get the new atom positions by rescaling
    call scale_positions(ib, scaling)

end subroutine expand_cell_ortho_by




!> most general cell distortion (NPT ensemble with non-orthogonal cell vectors)
subroutine distort_cell(ib, indx, maxvolchange, volnew)

    !use kinds_f90
    use latticevectors_module, only : invertlatvec, distortcell

    !implicit none

        !> id/index of the configuration box/cell to work on
    integer, intent(in) :: ib, indx

        !> current max V-step (displacement)
    real (kind = wp), intent(in) :: maxvolchange

        !> resulting attempted volume
    real (kind = wp), intent(out) :: volnew

        !> reshaping matrix
    real (kind = wp) :: bulks(9)

    call store_all_pos(ib)

!    if(cfgs(ib)%vec%is_orthogonal) then
    call set_rigid_mol_coms(ib) 
!    end if

    call distortcell(indx, maxvolchange, volnew, cfgs(ib)%vec, bulks)

    call scale_positions_strain(ib, bulks)

    call invertlatvec(cfgs(ib)%vec)

end subroutine distort_cell


!> scales molecule/atom positions by factor `scl`
subroutine scale_positions(ib, scl)

    !use kinds_f90
    use molecule_module, only : shift_molecule 
    !, mol_com, reform_molecule_frac, fractocart_mol, carttofrac_mol
    !use slit_module,     only : in_bulk

    !implicit none

        !> id/index of the configuration box/cell to work on
    integer, intent(in) :: ib

        !> scaling factor
    real(kind = wp), intent(in) :: scl

    integer :: im, i
    logical :: in_bulk

    real(kind = wp) :: shift(3), com_tmp(3)

    in_bulk = ( .not.myjob%is_XY_npt )

    shift(:)   = 0.0_wp
    com_tmp(:) = 0.0_wp

    do  im = 1, cfgs(ib)%num_mols

        ! if this is a rigid molecule translate c.o.m
        if (cfgs(ib)%mols(im)%rigid_body) then

            com_tmp(1) = cfgs(ib)%mols(im)%rcom(1) * scl
            com_tmp(2) = cfgs(ib)%mols(im)%rcom(2) * scl
            com_tmp(3) = cfgs(ib)%mols(im)%rcom(3)

            if( in_bulk ) com_tmp(3) = com_tmp(3) * scl

            shift(:) = com_tmp(:) - cfgs(ib)%mols(im)%rcom(:)

            ! move molecule
            call shift_molecule(cfgs(ib)%mols(im), shift)

        else

            do i = 1, cfgs(ib)%mols(im)%natom

                cfgs(ib)%mols(im)%atms(i)%rpos(1) = cfgs(ib)%mols(im)%atms(i)%rpos(1)*scl
                cfgs(ib)%mols(im)%atms(i)%rpos(2) = cfgs(ib)%mols(im)%atms(i)%rpos(2)*scl
                
                if( in_bulk ) cfgs(ib)%mols(im)%atms(i)%rpos(3) = &
                              cfgs(ib)%mols(im)%atms(i)%rpos(3)*scl

            enddo

            call set_mol_com(ib, im)

        endif

        if(cfgs(ib)%vec%is_orthogonal) then

            do i = 1, cfgs(ib)%mols(im)%natom
                call pbc_ortho_coor_calc(ib, &
                     cfgs(ib)%mols(im)%atms(i)%rpos(1), &
                     cfgs(ib)%mols(im)%atms(i)%rpos(2), &
                     cfgs(ib)%mols(im)%atms(i)%rpos(3) )
            end do

!AB: putting COM:s back into the primary cell only at molecule moves/rotations 
!AB: is at odds with atom moves where this is not done (or would add extra overheads)
!AB: this also defeats tracking COM's drift upon acceptance (done in all MC moves)
!AB: there is no need to keep COM:s inside the primary cell!

!            call pbc_ortho_coor_calc(ib, &
!                 cfgs(ib)%mols(im)%rcom(1), &
!                 cfgs(ib)%mols(im)%rcom(2), &
!                 cfgs(ib)%mols(im)%rcom(3) )
        end if

    enddo

end subroutine scale_positions



!> scales molecule/atom positions along a given dimension by factor `scl`
subroutine scale_positions_orthoani(ib, indx, scl)

    use molecule_module, only : shift_molecule 

        !> id/index of the configuration box/cell to work on
    integer, intent(in) :: ib

        !> dimension of the box to expand/contract (1=x,2=y,3=z)
    integer, intent(in) :: indx
    
        !> scaling factor
    real(kind = wp), intent(in) :: scl

    integer :: im, i

    real(kind = wp) :: shift(3), com_tmp(3)

    shift(:)   = 0.0_wp
    com_tmp(:) = 0.0_wp

    do  im = 1, cfgs(ib)%num_mols

        ! if this is a rigid molecule translate c.o.m
        if (cfgs(ib)%mols(im)%rigid_body) then

            com_tmp = cfgs(ib)%mols(im)%rcom(1)
            com_tmp(indx) = com_tmp(indx) * scl
            

            shift(:) = com_tmp(:) - cfgs(ib)%mols(im)%rcom(:)

            ! move molecule
            call shift_molecule(cfgs(ib)%mols(im), shift)

        else

            do i = 1, cfgs(ib)%mols(im)%natom

                cfgs(ib)%mols(im)%atms(i)%rpos(indx) = cfgs(ib)%mols(im)%atms(i)%rpos(indx)*scl

            enddo

            call set_mol_com(ib, im)

        endif

        if(cfgs(ib)%vec%is_orthogonal) then

            do i = 1, cfgs(ib)%mols(im)%natom
                call pbc_ortho_coor_calc(ib, &
                     cfgs(ib)%mols(im)%atms(i)%rpos(1), &
                     cfgs(ib)%mols(im)%atms(i)%rpos(2), &
                     cfgs(ib)%mols(im)%atms(i)%rpos(3) )
            end do

!AB: putting COM:s back into the primary cell only at molecule moves/rotations 
!AB: is at odds with atom moves where this is not done (or would add extra overheads)
!AB: this also defeats tracking COM's drift upon acceptance (done in all MC moves)
!AB: there is no need to keep COM:s inside the primary cell!

!            call pbc_ortho_coor_calc(ib, &
!                 cfgs(ib)%mols(im)%rcom(1), &
!                 cfgs(ib)%mols(im)%rcom(2), &
!                 cfgs(ib)%mols(im)%rcom(3) )
        end if

    enddo

end subroutine scale_positions_orthoani








subroutine scale_positions_strain(ib, bulks)

    !use kinds_f90
    use molecule_module, only : shift_molecule !, mol_com, reform_molecule_frac, fractocart_mol, carttofrac_mol
    use slit_module,     only : in_bulk

    !implicit none

    integer, intent(in) :: ib
    real(kind = wp), intent(in) :: bulks(6)

    integer :: im, i, j
    real(kind = wp) :: shift(3), com_tmp(3), cm(3), xx(3)

    do im = 1, cfgs(ib)%num_mols

        ! if this a rigid molecule translate c.o.m
        if (cfgs(ib)%mols(im)%rigid_body) then

            cm(:) = cfgs(ib)%mols(im)%rcom(:)

            com_tmp(1) = (1.0_wp + bulks(1)) * cm(1) + 0.5_wp * bulks(6) *    &
                 cm(2) + 0.5_wp * bulks(5)* cm(3)
            com_tmp(2) = (1.0_wp + bulks(2)) * cm(2) + 0.5_wp * bulks(6) *    &
                 cm(1) + 0.5_wp * bulks(4)* cm(3)
            com_tmp(3) = (1.0_wp + bulks(3)) * cm(3) + 0.5_wp * bulks(5) *    &
                 cm(1) + 0.5_wp * bulks(4)* cm(2)

            shift(:) =  com_tmp(:) - cfgs(ib)%mols(im)%rcom(:)

            call shift_molecule(cfgs(ib)%mols(im), shift)

        endif
        
        if(cfgs(ib)%vec%is_orthogonal) then

            do j = 1, cfgs(ib)%mols(im)%natom
                call pbc_ortho_coor_calc(ib, &
                     cfgs(ib)%mols(im)%atms(j)%rpos(1), &
                     cfgs(ib)%mols(im)%atms(j)%rpos(2), &
                     cfgs(ib)%mols(im)%atms(j)%rpos(3) )
            end do

!AB: putting COM:s back into the primary cell only at molecule moves/rotations 
!AB: is at odds with atom moves where this is not done (or would add extra overheads)
!AB: this also defeats tracking COM's drift upon acceptance (done in all MC moves)
!AB: there is no need to keep COM:s inside the primary cell!

!            call pbc_ortho_coor_calc(ib, &
!                 cfgs(ib)%mols(im)%rcom(1), &
!                 cfgs(ib)%mols(im)%rcom(2), &
!                 cfgs(ib)%mols(im)%rcom(3) )
        end if
    enddo

end subroutine scale_positions_strain


!TU-: The below could be improved by exploiting the fact that cfgs(ib)%nums_elemts(typ) gives
!TU-: the number of atoms in box 'ib' belonging to species 'typ'. If this is 0 then the search
!TU-: need not be performed.

!TU-: I need to look at this again

!> selects a movable atom at random
!AB: the success (after a possible short search) is ensured by the initial checks in setupmc(job)
subroutine select_species_atomov(ib, typ, k, j, i, atmovetyp, founda)

    !use kinds_f90
    use random_module, only : duni

    !implicit none

    integer, intent(in)    :: ib
    integer, intent(inout) :: typ
    integer, intent(out)   :: k, j, i
    integer, intent(inout) :: atmovetyp(:)
    logical, intent(out)   :: founda

    integer :: m, mitem, mloop

    typ = 0
    k = 0
    j = 0
    i = 0

    founda = .false.

    !AB: first, select a molecular species with movable atoms
    !AB: then, randomly pick up a molecule of that type and, finally, a movable atom in it

    mitem = size(cfgs(ib)%mtypes)
    mloop = 100*mitem
    do  m = 1, mloop

        k = int(mitem * duni() + 1.0_wp)
        if( cfgs(ib)%mtypes(k)%atoms_move ) exit

    end do

#ifdef DEBUG
!debugging
    if( cfgs(ib)%mtypes(k)%num_mols < 1 ) then
        write(*,*)"select_species_atomov(): movable by atoms molecule type ",k," is empty..."
        return
    endif

!debugging
    if( .not.cfgs(ib)%mtypes(k)%atoms_move ) then 
        write(*,*)"select_species_atomov(): *** SOS *** species ",k," does not contain movable atoms..."
        return
    endif
#endif

    j = int(cfgs(ib)%mtypes(k)%num_mols * duni() + 1.0_wp)
    j = cfgs(ib)%mtypes(k)%mol_id(j)

    mitem = cfgs(ib)%mols(j)%natom

    !TU: I've commented out the following warning. I assume it was intended to be used with the
    !TU: DEBUG flag?
    !if( mitem < 1 ) then 
    !    write(*,*)"select_species_atomov(): movable by atoms molecule ",j," does not contain any atom..."
    !    return
    !endif

    mloop = 100*mitem
    do  m = 1, mloop

        i = int(mitem * duni() + 1.0_wp)

        typ = cfgs(ib)%mols(j)%atms(i)%atlabel

        !AB: make sure the atom is amongst the movable ones
        do k = 1, myjob%numatmmove
            founda = ( typ == atmovetyp(k) )
            if( founda ) return 
            !debugging
            !then
            !    if( i /= k ) &
            !        write(*,*)"select_species_atomov(): moving atom ",i,&
            !                  " found under index ",k," in the atom-move list..."
            !    return
            !endif 
        enddo

    enddo

#ifdef DEBUG
!debugging
    write(*,*)"select_species_atomov(): *** SOS *** timed out while picking up a movable atom (not found)..."
    write(*,*)"m, typ, k, j, i, nmols = ", m, typ, k, j, i, cfgs(ib)%mtypes(typ)%num_mols, mitem
    !write(*,*)
#endif

end subroutine select_species_atomov


!TU-: I need to look at this again - it seems inelegant to call this a lot, if indeed that is the case

!> selects a molecule of a specific type at random
!AB: the success (after a possible short search) is ensured by the initial checks in setupmc(job)
subroutine select_species_molmov(ib, typ, k, j, molmovtyp, foundm)

    !use kinds_f90
    use random_module, only : duni

    !implicit none

    integer, intent(in)    :: ib
    integer, intent(inout) :: typ
    integer, intent(out)   :: k, j
    integer, intent(inout) :: molmovtyp(:)
    logical, intent(out)   :: foundm

    integer :: m, mitem, mloop

    typ = 0
    k = 0
    j = 0

    foundm = .false.

    mitem = size(cfgs(ib)%mtypes)
    mloop = 100*mitem
    do  m = 1, mloop

        typ = int(mitem * duni() + 1.0_wp)
        if( cfgs(ib)%mtypes(typ)%translate ) exit

    end do

    if( cfgs(ib)%mtypes(typ)%num_mols < 1 ) return

#ifdef DEBUG
!debugging
    if( cfgs(ib)%mtypes(typ)%num_mols < 1 ) then
        write(*,*)"select_species_molmov(): movable molecule type ",typ," is empty..."
        return
    endif

!debugging
    if( .not.cfgs(ib)%mtypes(typ)%translate ) then 
        write(*,*)"select_species_molmov() *** SOS *** species ",typ," does not contain movable molecules..."
        return
    endif
#endif

    j = int(cfgs(ib)%mtypes(typ)%num_mols * duni() + 1.0_wp)
    j = cfgs(ib)%mtypes(typ)%mol_id(j)

    !AB: make sure the molecule is amongst the movable ones
    do k = 1, myjob%nummolmove
        foundm = ( typ == molmovtyp(k) )
        if( foundm ) return 
        !debugging
        !then
        !    !if( i /= k ) &
        !        write(*,*)"select_species_molmov(): translated species ",i," found under index ",k," in the mol-move list..."
        !    return
        !endif 
    enddo
    
#ifdef DEBUG
!debugging
    write(*,*)"select_species_molmov() *** SOS *** timed out while picking up a movable molecule (not found)..."
    write(*,*)"m, typ, j, k, nmols = ", m, typ, j, k, cfgs(ib)%mtypes(typ)%num_mols
    !write(*,*)
#endif

end subroutine select_species_molmov


!TU-: I need to look at this again - it seems inelegant to call this a lot, if indeed that is the case

!> selects a molecule of a specific type at random
!AB: the success (after a possible short search) is ensured by the initial checks in setupmc(job)
subroutine select_species_molrot(ib, typ, k, j, molrottyp, foundm)

    !use kinds_f90
    use random_module, only : duni

    !implicit none

    integer, intent(in)    :: ib
    integer, intent(inout) :: typ
    integer, intent(out)   :: k, j
    integer, intent(inout) :: molrottyp(:)
    logical, intent(out)   :: foundm

    integer :: m, mitem, mloop

    typ = 0
    k = 0
    j = 0

    foundm = .false.

    mitem = size(cfgs(ib)%mtypes)
    mloop = 100*mitem
    do  m = 1, mloop

        typ = int(mitem * duni() + 1.0_wp)
        if( cfgs(ib)%mtypes(typ)%rotate ) exit

    end do

    if( cfgs(ib)%mtypes(typ)%num_mols < 1 ) return

#ifdef DEBUG
!debugging
    if( cfgs(ib)%mtypes(typ)%num_mols < 1 ) then
        write(*,*)"select_species_molrot(): rotatable molecule type ",typ," is empty..."
        return
    endif

!debugging
    if( .not.cfgs(ib)%mtypes(typ)%rotate ) then 
        write(*,*)"select_species_molrot(): *** SOS *** species ",typ," does not contain rotatable molecules..."
        return
    endif
#endif

    j = int(cfgs(ib)%mtypes(typ)%num_mols * duni() + 1.0_wp)
    j = cfgs(ib)%mtypes(typ)%mol_id(j)

    !AB: make sure the molecule is amongst the rotatable ones
    do k = 1, myjob%nummolrot
        foundm = ( typ == molrottyp(k) )
        if( foundm ) return 
        !debugging
        !then
        !    !if( i /= k ) &
        !        write(*,*)"select_species_molrot(): rotated species ",i," found under index ",k," in the mol-rot list..."
        !    return
        !endif 
    enddo

#ifdef DEBUG
!debugging
    write(*,*)"select_species_molrot(): *** SOS *** timed out while picking up a rotatable molecule (not found)..."
    write(*,*)"m, typ, j, k, nmols = ", m, typ, j, k, cfgs(ib)%mtypes(typ)%num_mols
    !write(*,*)
#endif

end subroutine select_species_molrot


!TU-: This is a very inelegant way to select an atom of a given type. I recall this is not employed
!TU-: for atom translation moves - the situation is a little better. However if the atom type to be
!TU-: located is very rare in the system then this procedure will take a long time to find such an
!TU-: atom. Also, only the atlabel is checked; we should explicitly forbid different atomic species
!TU-: from having the same label (e.g. forbid having two species 'A core' and 'A spin')
!TU-:
!TU-: Perhaps a better way to do this would be to keep an array holding the location (mol index and atom
!TU-: index within that molecule) of all atoms belonging to each atomic species. Then to pick an atom
!TU-: belonging to species 'A', one simply generates a random integer between 1 and N_A, say 'i', and
!TU-: takes the molecule and atom index from the ith element in the array corresponding to species
!TU-: 'A'.

!TU-: This hangs if there are no atoms of the specified species!

!> selects an atom of a specific type at random
!AB: the success (albeit after potentially long search) is ensured by the initial checks in setupmc(job)
subroutine select_atom(typ, ib, j, k)

    !use kinds_f90
    use random_module, only : duni

    !implicit none

    integer, intent(in) :: ib, typ

    integer, intent(out) :: j, k   

    integer :: i

    j = 0
    k = 0

    do 
        j = int(cfgs(ib)%num_mols * duni() + 1)
        k = int(cfgs(ib)%mols(j)%natom * duni() + 1)
        if (cfgs(ib)%mols(j)%atms(k)%atlabel == typ) return
    enddo

end subroutine select_atom


!TU-: I need to look at this warning again

!TU-: This hangs if there are no atoms of the specified species!

!> selects an atom of a specific type at random from a specific molecule*/
!TU: WARNING: The loop in this procedure will never finish if there is no atom of type 'typ'
!TU: in the molecule! Before calling this procedure one must ensure that this is not the case.
!AB: the success (albeit after potentially long search) is ensured by the initial checks in setupmc(job)
subroutine select_molatom(typ, ib, j, k)

    !use kinds_f90
    use random_module, only : duni

    !implicit none

    integer, intent(in) :: typ, ib, j

    integer, intent(out) :: k

    integer :: i

    k = 0

    do 
        k = int(cfgs(ib)%mols(j)%natom * duni() + 1)
        if (cfgs(ib)%mols(j)%atms(k)%atlabel == typ) return
    enddo

end subroutine select_molatom


!TU-: I need to look at this again - it seems inelegant

!TU: WARNING: The loop in this procedure will never finish if there is no molecule of type 'typ'
!TU: in the configuration! Before calling this procedure one must ensure that this is not the case.
!AB: the success (albeit after potentially long search) is ensured by the initial checks in setupmc(job)
subroutine select_molecule(typ, ib, imol)

    use random_module, only : duni

    !implicit none

    integer, intent(in) :: ib, typ

    integer, intent(out) :: imol

    integer :: i

    imol = 0

    do 
        imol = int(cfgs(ib)%num_mols * duni() + 1)
        if (cfgs(ib)%mols(imol)%mol_label == typ) return
    enddo

end subroutine select_molecule


!TU: This function is now obsolete; the counter array cfgs(ib)%nums_elemts
!TU: now does the same job. For now this function remains, with the aim of
!TU: it being removed in a future refactoring
!> returns the number of atoms of a given integer type in a configuration
subroutine findnumtype(ib, typ, num)

    !use kinds_f90

    !implicit none

    integer, intent(in) :: ib, typ

    integer, intent(out) :: num

    integer :: i, j

    !TU: As a stop-gap we just return the counter array element in this
    !TU: function. Later refactoring will remove this whole procedure, 
    !TU: replacing calls to this function with references to nums_elemts
    num = cfgs(ib)%nums_elemts(typ)

    !TU: Old version of this procedure...
    ! num = 0
    ! 
    ! do i = 1, cfgs(ib)%num_mols
    ! 
    !     do j = 1, cfgs(ib)%mols(i)%natom
    ! 
    !         if (cfgs(ib)%mols(i)%atms(j)%atlabel == typ) num = num + 1
    ! 
    !     enddo
    ! 
    ! enddo

end subroutine



!> Checks that the current counters for the number of atoms of each type in a 
!> configuration are correct, and returns an error if not. This procedure is
!> intended for catching bugs.
subroutine check_species_counters(ib)

    use species_module
    use constants_module
    use comms_mpi_module, only : master

    implicit none

        !> Configuration box number
    integer, intent(in) :: ib

    integer :: i, j, typ

    integer, allocatable :: true_nums_elemts(:)

    if( master) then

        write(uout,*)
        write(uout,*) "Checking counters for atomic types in box ",ib,"..."

    end if

    allocate(true_nums_elemts(number_of_elements))
    true_nums_elemts = 0

    do i = 1, cfgs(ib)%num_mols

        do j = 1, cfgs(ib)%mols(i)%natom

            typ = cfgs(ib)%mols(i)%atms(j)%atlabel
            true_nums_elemts(typ) = true_nums_elemts(typ) + 1

        end do

    enddo

    do typ = 1, number_of_elements

        if( true_nums_elemts(typ) /= cfgs(ib)%nums_elemts(typ) ) then

            if( master ) then

                write(uout,*)
                write(uout,*) "ERROR: Discrepancy in counts for element ",typ," in box ",ib,": true value = ", &
                              true_nums_elemts(typ),", counter = ",cfgs(ib)%nums_elemts(typ)

            end if
            call cry(uout,'', "ERROR: Counter for atomic elements in configuration has diverged " &
                 //"from true value.",999)

        end if

    end do

end subroutine




!TU-: I need to look at this again - it seems inelegant to call this a lot, if indeed that is the case


!TU: Added this to fix atomic GCMC bug.
!> returns the number of atoms of a given integer type in a specific molecule
subroutine findnumtype_in_mol(ib, typ, imol, num)

      !> configuration number
    integer, intent(in) :: ib
    
      !> integer type of atom
    integer, intent(in) :: typ
    
      !> molecule number in configuration
    integer, intent(in) :: imol

      !> number of atoms of the considered type in the molecule
    integer, intent(out) :: num

    integer :: j

    num = 0

    do j = 1, cfgs(ib)%mols(imol)%natom

        if (cfgs(ib)%mols(imol)%atms(j)%atlabel == typ) num = num + 1

    enddo

end subroutine


!*!
! * note that the identities are swapped rather than positions
! * this is so that the interaction lists do not need to be updated.
!*!

!> swap atom positions around 
subroutine swapatompositions(ib, m1, at1, m2, at2)

    !use kinds_f90
    use atom_module, only : swap_atom_types

    !implicit none

    integer, intent(in) :: ib, at1, at2, m1, m2
    integer :: tm, tsite, tatype, tlabel

    real(kind = wp) :: tmass, tcharge

    call swap_atom_types(cfgs(ib)%mols(m1)%atms(at1), cfgs(ib)%mols(m2)%atms(at2))

end subroutine

!> swap molecule positions around - or com's in reality
subroutine swapmolpositions(ib, m1, m2)

    !use kinds_f90
    use molecule_module, only : mol_com, swap_mol_positions

    !implicit none

    integer, intent(in) :: ib, m1, m2

    integer :: i
    real (kind = wp) :: cm(3)

    ! calculate com's of the molecules
    call mol_com(cfgs(ib)%mols(m1))
    call mol_com(cfgs(ib)%mols(m2))

    call swap_mol_positions(cfgs(ib)%mols(m1), cfgs(ib)%mols(m2))

    ! make sure image convention holds
    call pbc_molecule(ib, m1)
    call pbc_molecule(ib, m2)

end subroutine

!> change atom type for semi-grand ensemble
subroutine mutate_atom(ib, im, i, typ)

    use species_module, only : atm_charge, atm_mass, eletype
    use atom_module, only : set_atom_type

    !implicit none

    integer, intent(in) :: ib, i, im, typ

    integer :: oldtyp

    oldtyp = cfgs(ib)%mols(im)%atms(i)%atlabel
    cfgs(ib)%nums_elemts(typ) =  cfgs(ib)%nums_elemts(typ) + 1
    cfgs(ib)%nums_elemts(oldtyp) = cfgs(ib)%nums_elemts(oldtyp) - 1

    call set_atom_type(cfgs(ib)%mols(im)%atms(i), typ, atm_charge(typ), atm_mass(typ), eletype(typ))


end subroutine

!> change the molecule type, also provides random rotation
subroutine mutate_molecule(ib, im, typ, atom_out)

    !use kinds_f90
    use constants_module, only : PI
    use species_module, only : uniq_mol
    use molecule_module, only : copy_molecule, mol_com !, reform_molecule_frac, fractocart_mol

#ifdef DEBUG
!debugging!
    use slit_module, only : in_bulk
#endif

    !implicit none

    integer, intent(in) :: ib, im, typ

    logical, intent(out) :: atom_out

    integer :: i, atlabel
    real(kind = wp) :: cm(3)

    ! store the centre of mass
    call mol_com(cfgs(ib)%mols(im))
    cm(:) = cfgs(ib)%mols(im)%rcom(:)

    !TU: Amend cfgs(ib)%nums_elemts(typ) to reflect removal of 'old' typ of atom
    do i = 1, cfgs(ib)%mols(im)%natom

        atlabel =  cfgs(ib)%mols(im)%atms(i)%atlabel
        cfgs(ib)%nums_elemts(atlabel) = cfgs(ib)%nums_elemts(atlabel) - 1

    end do

    ! copy the molecule across
    call copy_molecule(cfgs(ib)%mols(im), uniq_mol(typ))


    !TU: Amend cfgs(ib)%nums_elemts(typ) to reflect addition of 'new' typ of atom
    do i = 1, cfgs(ib)%mols(im)%natom

        atlabel =  cfgs(ib)%mols(im)%atms(i)%atlabel
        cfgs(ib)%nums_elemts(atlabel) = cfgs(ib)%nums_elemts(atlabel) + 1

    end do

    ! move to the correct c.o.m
    call mol_com(cfgs(ib)%mols(im))

    do i = 1, cfgs(ib)%mols(im)%natom

        cfgs(ib)%mols(im)%atms(i)%rpos(:) = cfgs(ib)%mols(im)%atms(i)%rpos(:) - cfgs(ib)%mols(im)%rcom(:)
        cfgs(ib)%mols(im)%atms(i)%rpos(:) = cfgs(ib)%mols(im)%atms(i)%rpos(:) + cm(:)

    enddo

    call mol_com(cfgs(ib)%mols(im))

    ! apply random rotation
    !TU-: Here we can randomise orientations better by calling this multiple times
    !TU-: - or should we somehow try correlate the orientation of the mutated molecule
    !TU-:   with that of the current molecule. Does applying a random rotation violate
    !TU-:   detailed balance?
    call rotate_molecule(ib, im, PI, atom_out)

#ifdef DEBUG
!debugging!
    if( in_bulk .and. atom_out ) &
        write(*,*)'cell_module::mutate_molecule(..) ERROR!!! - undue atom_out = ', atom_out
#endif

    if( atom_out ) return

    call mol_com(cfgs(ib)%mols(im))


end subroutine


!> Function which converts a random fractional coordinate generated
!> uniformly between -0.5 and 0.5 into one uniformly between -0.5 and
!> 0.5, excluding the range between gcexcludeslab_ubound and gcexclude_ubound.
!> This function is used in insert_atom and insert_molecule if a
!> slab region is excluded from GCMC insertiona.
subroutine frac_coord_gcexclude(fpos)

        !> Fractional coordinate (between -0.5 and 0.5)
    real(kind=wp), intent(inout) :: fpos

    fpos = fpos - (myjob%gcexcludeslab_ubound - myjob%gcexcludeslab_lbound) * (fpos + 0.5_wp)
        
    if( fpos > myjob%gcexcludeslab_lbound ) then
        
        fpos = fpos + (myjob%gcexcludeslab_ubound - myjob%gcexcludeslab_lbound)
        
    end if

end subroutine frac_coord_gcexclude


!> inserts an atom into random position
subroutine insert_atom(ib, typ, imol, ia, lauto, atom_out)

    !use kinds_f90
    use species_module
    use molecule_module, only : add_atom
    use random_module, only : duni
    use slit_module, only : in_bulk, slit_zfrac1, slit_zfrac2 !, slit_gcmc_z
    use atom_module, only : rotate_spin_by
    use constants_module, only : PI, ATM_TYPE_SPIN
    
    !implicit none

    integer, intent(in)  :: ib, typ, imol

    integer, intent(out) :: ia

    logical, intent(in)  :: lauto

    logical, intent(out) :: atom_out

    integer :: ml, natm_cfg

    real(kind = wp) :: pos(3), buff(3), rot(9)

    atom_out = .false.

    natm_cfg = cfgs(ib)%number_of_atoms + 1

    if( natm_cfg > cfgs(ib)%maxno_of_atoms ) call error(203)
    
    ml = cfgs(ib)%mols(imol)%mol_label

    !AB: select a trial position in fractional coords

    buff(1) = duni() - 0.5_wp
    buff(2) = duni() - 0.5_wp
    buff(3) = duni() - 0.5_wp !anna bui 19th July 2024 

    if( myjob%usegcexcludeslab ) then

        call frac_coord_gcexclude(buff(myjob%gcexcludeslab_dim))

    end if

    !AB: convert fractional to Cartesian coords

    !AB: first, update and check Z coordinate in case atom gets out of the cell/slit
    pos(3) = cfgs(ib)%vec%latvector(1,3) * buff(1) &
           + cfgs(ib)%vec%latvector(2,3) * buff(2) &
           + cfgs(ib)%vec%latvector(3,3) * buff(3)

    atom_out = atom_outside(pos(3), cfgs(ib)%vec%latvector(3,3))

#ifdef DEBUG
!debugging!
    if( in_bulk .and. atom_out ) &
        write(*,*)'cell_module::insert_atom(..) ERROR!!! - undue atom_out = ', atom_out
#endif

    !AB: we have not inserted the atom yet!
    !AB: make sure that no atom is removed
    !AB: in the case of atom_out==.true. (see gcmc_moves.f90)

    if( atom_out ) return

    ! if still inside, proceed to update X & Y (here it does work in this order!)
    pos(1) = cfgs(ib)%vec%latvector(1,1) * buff(1) &
           + cfgs(ib)%vec%latvector(2,1) * buff(2) &
           + cfgs(ib)%vec%latvector(3,1) * buff(3)

    pos(2) = cfgs(ib)%vec%latvector(1,2) * buff(1) &
           + cfgs(ib)%vec%latvector(2,2) * buff(2) &
           + cfgs(ib)%vec%latvector(3,2) * buff(3)

    cfgs(ib)%number_of_atoms = natm_cfg
    cfgs(ib)%mtypes(ml)%num_elems(typ) = cfgs(ib)%mtypes(ml)%num_elems(typ)+1

    call add_atom(cfgs(ib)%mols(imol), typ, pos, &
         atm_charge(typ), atm_mass(typ), eletype(typ), 0, lauto, ia)

    cfgs(ib)%mols(imol)%glob_no(ia) = cfgs(ib)%number_of_atoms

    !TU: If the atom has an orientation then randomise it
    if ( cfgs(ib)%mols(imol)%atms(ia)%atype == ATM_TYPE_SPIN ) then

        !TU-: This could be done multiple times to randomise it better
        cfgs(ib)%mols(imol)%atms(ia)%spin = (/ 0.0_wp, 0.0_wp, 1.0_wp /)
        rot = generate_rot_matrix_quat(PI)
        call rotate_spin_by(cfgs(ib)%mols(imol)%atms(ia), rot)        
        
    end if
    
    !TU: 'num_vol_sites' must be correct if Gibbs volume moves are in use, hence the
    !TU: below amendment. Note that I assume the inserted atom is a
    !TU: 'moving particle' for the sake of volume moves - a reasonable assumption
    cfgs(ib)%num_vol_sites = cfgs(ib)%num_vol_sites + 1

    cfgs(ib)%nums_elemts(typ) =  cfgs(ib)%nums_elemts(typ) + 1

end subroutine


!> find a "free cavity" and insert an atom at the cavity position
subroutine insert_atom_atcavity(ib, typ, imol, ia, ig, prob, lauto, atom_out)

    !use kinds_f90
    use constants_module, only : HALF, uout, PI, ATM_TYPE_SPIN
    use species_module
    use molecule_module, only : add_atom
    use random_module, only : duni
    use atom_module, only : rotate_spin_by
    
#ifdef DEBUG
!debugging!
    use slit_module, only : in_bulk
#endif

    !implicit none

    integer, intent(in)    :: ib, typ, imol

    integer, intent(inout) :: ia, ig

    real(kind = wp), intent(inout) :: prob

    logical, intent(in)    :: lauto

    logical, intent(inout) :: atom_out

    integer, save :: natt = 0

    integer, save :: nout = 0

    integer :: ml, natm_cfg, ic, it

    real(kind = wp) :: pos(3), rot(9) 

    logical :: found

    natt = natt + 1

    atom_out = .false.

    ia = 0
    it = 0
    ig = 0

    prob = -1.0_wp

    ! check if there is free space

    if( ncavity_free(ib) < 1 ) then
        nout = nout + 1
        return
    endif

    natm_cfg = cfgs(ib)%number_of_atoms + 1
    !cfgs(ib)%number_of_atoms = cfgs(ib)%number_of_atoms + 1

    if( natm_cfg > cfgs(ib)%maxno_of_atoms ) call error(203)
    !if(cfgs(ib)%number_of_atoms > cfgs(ib)%maxno_of_atoms) call error(203)

    ml = cfgs(ib)%mols(imol)%mol_label

#ifdef DEBUG
!debug<

    found = .false.

    do while( .not.found  )

        !AB: original selection of a random site (required iteration)
        !ig = int(duni() * ncavity_grid(ib) + 1)
        !if( cfgs(ib)%cavity(ig) == 0 ) found = .true.
!debug>
#endif
        !AB: select a free site on the cavity grid (now definite)
        ic = int(duni() * ncavity_free(ib)) + 1
        ig = cfgs(ib)%cavids(ic)

        if( cfgs(ib)%cavity(ig) > 0 ) then

            write(uout,*)"insert_atom_atcavity(): "//&
            "cavity selected for insertion is not free!!! - attempt ",&
            natt," after ",nout," cavity failures"

            prob = -1.0_wp
        endif
#ifdef DEBUG
!debug<
        found = ( cfgs(ib)%cavity(ig) == 0 )

        it = it + 1

        ! only allow 10000 searches
        if(it > 10000) then

            !cfgs(ib)%number_of_atoms = cfgs(ib)%number_of_atoms - 1
            !cfgs(ib)%mtypes(ml)%num_elems(typ) = cfgs(ib)%mtypes(ml)%num_elems(typ)-1

            ig = 0
            prob = -1.0_wp

            write(uout,*)"insert_atom_atcavity(): "//&
            "free cavity for insertion not found after 10000 attempts - attempt ",&
            natt," after ",nout," cavity failures"

            return

        endif

    enddo
!debug>
#endif

    pos(1) = cfgs(ib)%cavgrid(1,ig)
    pos(2) = cfgs(ib)%cavgrid(2,ig)
    pos(3) = cfgs(ib)%cavgrid(3,ig)

    if( myjob % cavity_mesh ) then

        pos(1) = pos(1) + (duni()-HALF)*cfgs(ib)%drcav(1)
        pos(2) = pos(2) + (duni()-HALF)*cfgs(ib)%drcav(2)
        pos(3) = pos(3) + (duni()-HALF)*cfgs(ib)%drcav(3)

    endif

    !AB: first, check Z coordinate in case atom gets out of the cell/slit (here it does work in this order!)
    !atom_out = atom_outside(cfgs(ib)%cavgrid(3,ig), cfgs(ib)%vec%latvector(3,3))
    atom_out = atom_outside(pos(3), cfgs(ib)%vec%latvector(3,3))

#ifdef DEBUG
!debugging!
    if( in_bulk .and. atom_out ) &
        write(uout,*)"cell_module::insert_atom_atcavity(..) "//&
                     "ERROR!!! - undue atom_out = ", atom_out
#endif

    !AB: we have not inserted the atom yet!
    !AB: make sure that no atom is removed, nor cavity grid reset/updated
    !AB: in the case of atom_out=.true. (see gcmc_moves.f90)

    if( atom_out ) return

    ! if still inside, proceed to adding atom (here it does work in this order!)

    cfgs(ib)%number_of_atoms = natm_cfg
    cfgs(ib)%mtypes(ml)%num_elems(typ) = cfgs(ib)%mtypes(ml)%num_elems(typ)+1

    !call add_atom(cfgs(ib)%mols(imol), typ, cfgs(ib)%cavgrid(:,ig), &
    call add_atom(cfgs(ib)%mols(imol), typ, pos, &
         atm_charge(typ), atm_mass(typ), eletype(typ), 0, lauto, ia)

    cfgs(ib)%mols(imol)%glob_no(ia) = cfgs(ib)%number_of_atoms

    !TU: If the atom has an orientation then randomise it
    if ( cfgs(ib)%mols(imol)%atms(ia)%atype == ATM_TYPE_SPIN ) then

        !TU-: This could be done multiple times to randomise it better
        cfgs(ib)%mols(imol)%atms(ia)%spin = (/ 0.0_wp, 0.0_wp, 1.0_wp /)
        rot = generate_rot_matrix_quat(PI)
        call rotate_spin_by(cfgs(ib)%mols(imol)%atms(ia), rot)        
        
    end if
    
    !TU: 'num_vol_sites' must be correct if Gibbs volume moves are in use, hence the
    !TU: below amendment. Note that I assume the inserted atom is a
    !TU: 'moving particle' for the sake of volume moves - a reasonable assumption
    cfgs(ib)%num_vol_sites = cfgs(ib)%num_vol_sites + 1

    ! cavity now occupied

    call calc_pcav_insert(ib, ig, myjob%cavity_radius, pos, prob, .true.)
    !call calc_pcav_insert(ib, imol, ia, ig, myjob%cavity_radius, prob)

!AB: Occupying a single cavity upon insertion is only valid if min{dx,dy,dz} >= D_core (approx),
!AB: otherwise a particle can occupy more than one node on a grid (or more than one cell on a mesh)
!AB: BEWARE, though, letting min{dx,dy,dz} >= D_core will produce undersampling artifacts, 
!AB: unless other (regular) MC moves are also attempted in the same simulation; generally, 
!AB: undersampling of the configurational space is the major danger when applying cavity bias!

    !cfgs(ib)%cavity(ig) = 1
    !cfgs(ib)%cavids(ic) = cfgs(ib)%cavids(ncavity_free(ib))
    !cfgs(ib)%cavids(ncavity_free(ib)) = 0
    !ncavity_free(ib) = ncavity_free(ib) - 1
    !prob = real(ncavity_free(ib),wp) / real(ncavity_grid(ib),wp)


    cfgs(ib)%nums_elemts(typ) =  cfgs(ib)%nums_elemts(typ) + 1

end subroutine


!AB: new variant of molecule insertion routine (see below for the old/original one)
subroutine insert_molecule(ib, ml, imol, atom_out, as_whole)

    !use kinds_f90
    use constants_module, only : uout, PI
    use molecule_type
    use molecule_module, only : copy_molecule, mol_com, reform_molecule_ortho, full_mol_com 
    use species_module
    use slit_module, only : in_bulk, in_slit, slit_zfrac1, slit_zfrac2 !, slit_gcmc_z
    use random_module, only : duni

    !implicit none

    integer, intent(in) :: ib, ml

    integer, intent(out) :: imol

    logical, intent(inout) :: atom_out

    logical, intent(in)    :: as_whole

    integer :: i, counter, atlabel

    real(kind = wp) :: pos(3), buff(3), rot(9)

    atom_out = .false.

    imol = cfgs(ib)%num_mols + 1

    ! check the total no of molecules has not been exceeded
    if( imol > cfgs(ib)%mxmol ) then 

        atom_out = .true.

        call error(204)

    end if

    ! check the total no of atoms has not been exceeded
    if( (cfgs(ib)%number_of_atoms + uniq_mol(ml)%natom) > cfgs(ib)%maxno_of_atoms ) then

        atom_out = .true.

        call error(203)

    endif

    cfgs(ib)%num_mols = imol
    cfgs(ib)%mtypes(ml)%num_mols = cfgs(ib)%mtypes(ml)%num_mols + 1
    cfgs(ib)%mtypes(ml)%mol_id( cfgs(ib)%mtypes(ml)%num_mols ) = cfgs(ib)%num_mols

    buff(:) = 0.0_wp

    !AB: create a molecule by duplicating the template
    call copy_molecule(cfgs(ib)%mols(imol), uniq_mol(ml))

    !TU: Amend cfgs(ib)%num_elemts(typ) to reflect addition of new molecule
    do i = 1, cfgs(ib)%mols(imol)%natom

        atlabel =  cfgs(ib)%mols(imol)%atms(i)%atlabel
        cfgs(ib)%nums_elemts(atlabel) = cfgs(ib)%nums_elemts(atlabel) + 1

    end do


    !AB: check if the molecule has structure/topology to it
    !if( uniq_mol(ml)%natom > 1 .and. &
    !  ( uniq_mol(ml)%rigid_body .or. uniq_mol(ml)%blist%npairs > 0 ) ) then
    if( as_whole ) then

        ! give the atoms a global number
        do  i = 1, cfgs(ib)%mols(imol)%natom

            cfgs(ib)%mols(imol)%glob_no(i) = cfgs(ib)%number_of_atoms + i


        enddo

        !AB: get com-position for a new molecule
        buff(1) = duni() - 0.5_wp
        buff(2) = duni() - 0.5_wp
        buff(3) = duni() - 0.5_wp !Anna Bui 19th July 2024 


        if( myjob%usegcexcludeslab ) then
    
            call frac_coord_gcexclude(buff(myjob%gcexcludeslab_dim))
    
        end if

        if( myjob%useorthogonal ) then

            pos(1) = cfgs(ib)%vec%latvector(1,1) * buff(1)
            pos(2) = cfgs(ib)%vec%latvector(2,2) * buff(2)
            pos(3) = cfgs(ib)%vec%latvector(3,3) * buff(3)

            !if( uniq_mol(ml)%blist%npairs > 0 ) &
            call reform_molecule_ortho(cfgs(ib)%mols(imol),cfgs(ib)%vec%latvector)
            call mol_com(cfgs(ib)%mols(imol))

        else

            pos(1) = cfgs(ib)%vec%latvector(1,1) * buff(1) + & 
                     cfgs(ib)%vec%latvector(2,1) * buff(2) + & 
                     cfgs(ib)%vec%latvector(3,1) * buff(3)

            pos(2) = cfgs(ib)%vec%latvector(1,2) * buff(1) + & 
                     cfgs(ib)%vec%latvector(2,2) * buff(2) + & 
                     cfgs(ib)%vec%latvector(3,2) * buff(3)

            pos(3) = cfgs(ib)%vec%latvector(1,3) * buff(1) + & 
                     cfgs(ib)%vec%latvector(2,3) * buff(2) + & 
                     cfgs(ib)%vec%latvector(3,3) * buff(3)

            call full_mol_com(cfgs(ib)%mols(imol), cfgs(ib)%vec%invlat,cfgs(ib)%vec%latvector)

        end if

        atom_out = position_molecule_out(cfgs(ib)%mols(imol), pos, cfgs(ib)%vec%latvector(3,3))

        ! rotate molecule
        if (myjob%usequaternion) then

            !TU-: Here we can randomise orientations better by calling this multiple times
            call rotate_molecule_quat(ib, imol, PI, atom_out, rot)

        else

            !TU-: Here we can randomise orientations better by calling this multiple times
            call rotate_molecule(ib, imol, PI, atom_out)

        endif

        ! increase the total no of atoms in the simulation by the number of atoms in the new molecule
        cfgs(ib)%number_of_atoms = cfgs(ib)%number_of_atoms + uniq_mol(ml)%natom

    else
    !AB: insert unstructured molecule's atoms at random positions

        do i = 1, cfgs(ib)%mols(imol)%natom
           
           ! give the atoms a global number
           cfgs(ib)%mols(imol)%glob_no(i) = cfgs(ib)%number_of_atoms + i

           buff(1) = duni() - 0.5_wp
           buff(2) = duni() - 0.5_wp
           buff(3) = duni() - 0.5_wp !Anna Bui 19th July 2024 

           
           if( myjob%usegcexcludeslab ) then
        
               call frac_coord_gcexclude(buff(myjob%gcexcludeslab_dim))        

           end if


           if( myjob%useorthogonal ) then

               pos(1) = cfgs(ib)%vec%latvector(1,1) * buff(1)
               pos(2) = cfgs(ib)%vec%latvector(2,2) * buff(2)
               pos(3) = cfgs(ib)%vec%latvector(3,3) * buff(3)

           else

               pos(1) = cfgs(ib)%vec%latvector(1,1) * buff(1) + & 
                        cfgs(ib)%vec%latvector(2,1) * buff(2) + &
                        cfgs(ib)%vec%latvector(3,1) * buff(3)

               pos(2) = cfgs(ib)%vec%latvector(1,2) * buff(1) + & 
                        cfgs(ib)%vec%latvector(2,2) * buff(2) + & 
                        cfgs(ib)%vec%latvector(3,2) * buff(3)

               pos(3) = cfgs(ib)%vec%latvector(1,3) * buff(1) + & 
                        cfgs(ib)%vec%latvector(2,3) * buff(2) + & 
                        cfgs(ib)%vec%latvector(3,3) * buff(3)

           end if

           cfgs(ib)%mols(imol)%atms(i)%rpos(:)       = pos(:)
           cfgs(ib)%mols(imol)%atms(i)%store_rpos(:) = pos(:)

#ifdef DEBUG
!debugging!
!           write(*,*)'cell_module::insert_molecule(..) - attempted insertion of field-molecule ', &
!                     ml," no.",imol,", atom ",i," at z = ",buff(3)," / ",pos(3)
!                     !ml," no.",imol,", atom ",i,", glob ",cfgs(ib)%mols(imol)%glob_no(i)," at z = ",buff(3)
#endif
        enddo

        ! increase the total no of atoms in the simulation by the number of atoms in the new molecule
        cfgs(ib)%number_of_atoms = cfgs(ib)%number_of_atoms + uniq_mol(ml)%natom

        call mol_com(cfgs(ib)%mols(imol))

    end if

    !TU: 'num_vol_sites' must be correct if Gibbs volume moves are in use, hence the
    !TU: below amendment. NOTE THAT I ASSUME THAT THE INSERTED MOLECULE IS A SINGLE 'MOVING
    !TU: PARTICLE' for the sake of volume moves - a reasonable assumption for rigid molecules
    !TU: but not necessarily so for flexible molecules or something else
    cfgs(ib)%num_vol_sites = cfgs(ib)%num_vol_sites + 1

#ifdef DEBUG
!debugging!
    if( in_bulk .and. atom_out ) &
        write(*,*)'cell_module::insert_molecule(..) ERROR!!! - undue atom_out = ', atom_out
#endif


end subroutine


!TU-: The 3rd argument 'typ' is not used and should be removed; 'ml' is the molecular species
!AB: new variant of molecule insertion routine (see below for the old one)
subroutine insert_molecule_atcavity(ib, ml, imol, typ, ig, prob, atom_out, as_whole)

    !use kinds_f90
    use constants_module, only : PI, HALF, DZERO, FZERO, uout
    use molecule_type
    use molecule_module, only : copy_molecule, mol_com, reform_molecule_ortho, full_mol_com!, carttofrac_mol
    use species_module
    use slit_module, only : in_bulk, in_slit, slit_zfrac1, slit_zfrac2 !, slit_gcmc_z
    use random_module, only : duni

    !implicit none

    integer, intent(in) :: ib, ml, typ

    integer, intent(inout) :: imol, ig

    real(kind = wp), intent(inout) :: prob

    logical, intent(inout) :: atom_out

    logical, intent(in)    :: as_whole

    integer :: i, ic, it, counter, atlabel

    real(kind = wp) :: rx, ry, rz, pcav

    real(kind = wp) :: pos(3), buff(3), rot(9)

    logical :: found

    integer, save :: natt = 0

    integer, save :: nout = 0

    natt = natt + 1

    atom_out = .false.

    ig = 0

    prob = -1.0_wp

    ! check if there is free space (AB: what if not? - rejection???)

    if( ncavity_free(ib) < 1 ) then
        nout = nout + 1
        return
    end if

    imol = cfgs(ib)%num_mols + 1

    ! check the total no of molecules has not been exceeded
    if( imol > cfgs(ib)%mxmol ) then 

        atom_out = .true.

        call error(204)

    end if

    ! check the total no of atoms has not been exceeded
    if( (cfgs(ib)%number_of_atoms + uniq_mol(ml)%natom) > cfgs(ib)%maxno_of_atoms ) then

        atom_out = .true.

        call error(203)

    endif

!AB: below the actual insertion takes place

    cfgs(ib)%num_mols = imol
    cfgs(ib)%mtypes(ml)%num_mols = cfgs(ib)%mtypes(ml)%num_mols + 1
    cfgs(ib)%mtypes(ml)%mol_id( cfgs(ib)%mtypes(ml)%num_mols ) = cfgs(ib)%num_mols

    buff(:) = 0.0_wp

    !AB: create a molecule by duplicating the template
    call copy_molecule(cfgs(ib)%mols(imol), uniq_mol(ml))

    !TU: Amend cfgs(ib)%num_elemts(typ) to reflect addition of new molecule
    do i = 1, cfgs(ib)%mols(imol)%natom

        atlabel =  cfgs(ib)%mols(imol)%atms(i)%atlabel
        cfgs(ib)%nums_elemts(atlabel) = cfgs(ib)%nums_elemts(atlabel) + 1

    end do

    prob = 1.0_wp

    !AB: check if the molecule has structure/topology to it
    !if( uniq_mol(ml)%natom > 1 .and. &
    !  ( uniq_mol(ml)%rigid_body .or. uniq_mol(ml)%blist%npairs > 0 ) ) then
    if( as_whole ) then

        ! give the atoms a global number
        do  i = 1, cfgs(ib)%mols(imol)%natom

            cfgs(ib)%mols(imol)%glob_no(i) = cfgs(ib)%number_of_atoms + i

        enddo

        pcav  = 1.0_wp

#ifdef DEBUG
!debug<
        ic = 0
        it = 0

        found = .false.

        do while( .not.found  )

            !AB: original selection of a random site (required iteration)
            !ig = int(duni() * ncavity_grid(ib) + 1)
            !if( cfgs(ib)%cavity(ig) == 0 ) found = .true.
!debug>
#endif
            !AB: select a free site on the cavity grid (now definite)
            ic = int(duni() * ncavity_free(ib)) + 1
            ig = cfgs(ib)%cavids(ic)

            if( cfgs(ib)%cavity(ig) > 0 ) then
                pcav = -1.0_wp
                write(uout,*)"insert_molecule_atcavity():: "//&
                             "cavity selected for molecule insertion by com is not free!!! - attempt ",&
                             natt," after ",nout," cavity failures; pcav = ",pcav
!            else
!                write(uout,*)"insert_molecule_atcavity():: "//&
!                             "selected cavity ",ig," for insertion of molecule ",imol," by com - attempt ",&
!                             natt," after ",nout," cavity failures; pcav = ",pcav
            endif
#ifdef DEBUG
!debug<
            found = ( cfgs(ib)%cavity(ig) == 0 )

            it = it + 1

            ! only allow 10000 searches
            if(it > 10000) then

                !cfgs(ib)%number_of_atoms = cfgs(ib)%number_of_atoms - 1
                !cfgs(ib)%mtypes(ml)%num_elems(typ) = cfgs(ib)%mtypes(ml)%num_elems(typ)-1

                ig = 0
                pcav = -1.0_wp

                write(uout,*)"insert_molecule_atcavity():: "//&
                             "free cavity for insertion not found after 10000 attempts - attempt ",&
                             natt," after ",nout," cavity failures; pcav = ",pcav

                exit
                !return

            endif

        enddo
!debug>
#endif

        if( pcav < 0.0_wp ) then

            !AB: get com-position for a new molecule
            buff(1) = duni() - 0.5_wp
            buff(2) = duni() - 0.5_wp
            buff(3) = duni() - 0.5_wp !Anna Bui 19th July 2024 

            if( myjob%useorthogonal ) then

                pos(1) = cfgs(ib)%vec%latvector(1,1) * buff(1)
                pos(2) = cfgs(ib)%vec%latvector(2,2) * buff(2)
                pos(3) = cfgs(ib)%vec%latvector(3,3) * buff(3)

                call reform_molecule_ortho(cfgs(ib)%mols(imol),cfgs(ib)%vec%latvector)
                call mol_com(cfgs(ib)%mols(imol))

            else

                pos(1) = cfgs(ib)%vec%latvector(1,1) * buff(1) + & 
                         cfgs(ib)%vec%latvector(2,1) * buff(2) + & 
                         cfgs(ib)%vec%latvector(3,1) * buff(3)

                pos(2) = cfgs(ib)%vec%latvector(1,2) * buff(1) + & 
                         cfgs(ib)%vec%latvector(2,2) * buff(2) + & 
                         cfgs(ib)%vec%latvector(3,2) * buff(3)

                pos(3) = cfgs(ib)%vec%latvector(1,3) * buff(1) + & 
                         cfgs(ib)%vec%latvector(2,3) * buff(2) + & 
                         cfgs(ib)%vec%latvector(3,3) * buff(3)

                call full_mol_com(cfgs(ib)%mols(imol), cfgs(ib)%vec%invlat,cfgs(ib)%vec%latvector)

            end if

        else

            pos(1) = cfgs(ib)%cavgrid(1,ig)
            pos(2) = cfgs(ib)%cavgrid(2,ig)
            pos(3) = cfgs(ib)%cavgrid(3,ig)

            if( myjob % cavity_mesh ) then

                pos(1) = pos(1) + (duni()-HALF)*cfgs(ib)%drcav(1)
                pos(2) = pos(2) + (duni()-HALF)*cfgs(ib)%drcav(2)
                pos(3) = pos(3) + (duni()-HALF)*cfgs(ib)%drcav(3)

            endif

            !call calc_pcav_insert(ib, ig, myjob%cavity_radius, pos, pcav, .true.)
            !prob = pcav

            if( myjob%useorthogonal ) then

                !if( uniq_mol(ml)%blist%npairs > 0 ) &
                call reform_molecule_ortho(cfgs(ib)%mols(imol),cfgs(ib)%vec%latvector)
                call mol_com(cfgs(ib)%mols(imol))

            else

                call full_mol_com(cfgs(ib)%mols(imol), cfgs(ib)%vec%invlat,cfgs(ib)%vec%latvector)

            end if

        endif

        atom_out = position_molecule_out(cfgs(ib)%mols(imol), pos, cfgs(ib)%vec%latvector(3,3))

        ! rotate molecule
        if (myjob%usequaternion) then
            !TU-: Here we can randomise orientations better by calling this multiple times
            call rotate_molecule_quat(ib, imol, PI, atom_out, rot)
        else
            !TU-: Here we can randomise orientations better by calling this multiple times
            call rotate_molecule(ib, imol, PI, atom_out)
        endif

!AB: check/update the cavity grid and pcav ...

!#ifdef DEBUG
!debug<
        if( (cfgs(ib)%mols(imol)%rcom(1)-pos(1))**2 > DZERO .or. &
            (cfgs(ib)%mols(imol)%rcom(2)-pos(2))**2 > DZERO .or. &
            (cfgs(ib)%mols(imol)%rcom(3)-pos(3))**2 > DZERO ) then
            write(uout,*)"insert_molecule_atcavity():: "//&
                         "molecule COM if off the selected free cavity position at attempt ", &
                         natt," after ",nout," cavity failures"
!        else
!            write(uout,*)"insert_molecule_atcavity():: molecule COM insertion attempt ", &
!                         natt," after ",nout," cavity failures; pcav = ",pcav
        endif
!debug>
!#endif

        if( pcav > 0.0_wp ) then
            !call calc_pcav_insert(ib, ig, myjob%cavity_radius, &
            !     cfgs(ib)%mols(imol)%rcom, pcav, .true.)

            call calc_pcav_insert_mol(ib, ig, imol, myjob%cavity_radius, pcav, .true.)
            prob = pcav
        end if

        ! increase the total no of atoms in the simulation by the number of atoms in the new molecule
        cfgs(ib)%number_of_atoms = cfgs(ib)%number_of_atoms + uniq_mol(ml)%natom

    else
    !AB: insert unstructured molecule's atoms at random positions

        do i = 1, cfgs(ib)%mols(imol)%natom

            if( cfgs(ib)%mols(imol)%atms(i)%mass < FZERO ) &
                write(*,*)&
                "insert_molecule_atcavity():: ERROR!!! attempt to insert fictitious atom ",i, &
                " in unstructured molecule ",imol," (zero mass within atomic field?)"

            ! give the atoms a global number
            cfgs(ib)%mols(imol)%glob_no(i) = cfgs(ib)%number_of_atoms + i

!AB: cavity bias for atom positions:

            pcav  = 1.0_wp

#ifdef DEBUG
!debug<
            ic = 0
            it = 0

            found = .false.

            do while( .not.found  )

                !AB: original selection of a random site (required iteration)
                !ig = int(duni() * ncavity_grid(ib) + 1)
                !if( cfgs(ib)%cavity(ig) == 0 ) found = .true.
!debug>
#endif
                !AB: select a free site on the cavity grid (now definite)
                ic = int(duni() * ncavity_free(ib)) + 1
                ig = cfgs(ib)%cavids(ic)

                if( cfgs(ib)%cavity(ig) > 0 ) then
                    pcav = -1.0_wp
                    write(uout,*)"insert_molecule_atcavity():: "//&
                                 "cavity selected for atom insertion is not free!!! - attempt ",&
                                 natt," after ",nout," cavity failures; pcav = ",pcav
!                else
!                    write(uout,*)"insert_molecule_atcavity():: "//&
!                                 "selected cavity ",ig," for insertion of atom ",i," in molecule ",imol, &
!                                 " - attempt ",natt," after ",nout," cavity failures; pcav = ",pcav
                endif
#ifdef DEBUG
!debug<
                found = ( cfgs(ib)%cavity(ig) == 0 )

                it = it + 1

                ! only allow 10000 searches
                if(it > 10000) then

                    !cfgs(ib)%number_of_atoms = cfgs(ib)%number_of_atoms - 1
                    !cfgs(ib)%mtypes(ml)%num_elems(typ) = cfgs(ib)%mtypes(ml)%num_elems(typ)-1

                    ig = 0
                    pcav = -1.0_wp

                    write(uout,*)"insert_molecule_atcavity():: "//&
                                 "free cavity for insertion not found after 10000 attempts - attempt ",&
                                 natt," after ",nout," cavity failures; pcav = ",pcav

                    exit
                    !return

                endif

            enddo
!debug>
#endif
            if( pcav < 0.0_wp ) then

               buff(1) = duni() - 0.5_wp
               buff(2) = duni() - 0.5_wp
               buff(3) = duni() - 0.5_wp !Anna Bui 19th July 2024 

               if( myjob%useorthogonal ) then

                   pos(1) = cfgs(ib)%vec%latvector(1,1) * buff(1)
                   pos(2) = cfgs(ib)%vec%latvector(2,2) * buff(2)
                   pos(3) = cfgs(ib)%vec%latvector(3,3) * buff(3)

               else

                   pos(1) = cfgs(ib)%vec%latvector(1,1) * buff(1) + & 
                            cfgs(ib)%vec%latvector(2,1) * buff(2) + &
                            cfgs(ib)%vec%latvector(3,1) * buff(3)

                   pos(2) = cfgs(ib)%vec%latvector(1,2) * buff(1) + & 
                            cfgs(ib)%vec%latvector(2,2) * buff(2) + & 
                            cfgs(ib)%vec%latvector(3,2) * buff(3)

                   pos(3) = cfgs(ib)%vec%latvector(1,3) * buff(1) + & 
                            cfgs(ib)%vec%latvector(2,3) * buff(2) + & 
                            cfgs(ib)%vec%latvector(3,3) * buff(3)

               end if

            else

                pos(1) = cfgs(ib)%cavgrid(1,ig)
                pos(2) = cfgs(ib)%cavgrid(2,ig)
                pos(3) = cfgs(ib)%cavgrid(3,ig)

                if( myjob % cavity_mesh ) then

                    pos(1) = pos(1) + (duni()-HALF)*cfgs(ib)%drcav(1)
                    pos(2) = pos(2) + (duni()-HALF)*cfgs(ib)%drcav(2)
                    pos(3) = pos(3) + (duni()-HALF)*cfgs(ib)%drcav(3)

                endif

!AB: check/update the cavity grid and pcav ...

                if( i == 1 ) then
                    call calc_pcav_insert(ib, ig, myjob%cavity_radius, pos, pcav, .true.)
                else
                    call calc_pcav_insert(ib, ig, myjob%cavity_radius, pos, pcav, .false.)
                endif

                prob = prob*pcav

            endif

            cfgs(ib)%mols(imol)%atms(i)%rpos(:)       = pos(:)
            cfgs(ib)%mols(imol)%atms(i)%store_rpos(:) = pos(:)

#ifdef DEBUG
!debugging!
!            write(*,*)'cell_module::insert_molecule(..) - attempted insertion of field-molecule ', &
!                      ml," no.",imol,", atom ",i," at z = ",buff(3)," / ",pos(3)
!                      !ml," no.",imol,", atom ",i,", glob ",cfgs(ib)%mols(imol)%glob_no(i)," at z = ",buff(3)
#endif
        enddo

        ! increase the total no of atoms in the simulation by the number of atoms in the new molecule
        cfgs(ib)%number_of_atoms = cfgs(ib)%number_of_atoms + uniq_mol(ml)%natom

        call mol_com(cfgs(ib)%mols(imol))

    end if

    !TU: 'num_vol_sites' must be correct if Gibbs volume moves are in use, hence the
    !TU: below amendment. NOTE THAT I ASSUME THAT THE INSERTED MOLECULE IS A SINGLE 'MOVING
    !TU: PARTICLE' for the sake of volume moves - a reasonable assumption for rigid molecules
    !TU: but not necessarily so for flexible molecules or something else
    cfgs(ib)%num_vol_sites = cfgs(ib)%num_vol_sites + 1

#ifdef DEBUG
!debugging!
    if( in_bulk .and. atom_out ) &
        write(*,*)&
        "insert_molecule_atcavity():: ERROR!!! undue atom_out = ", atom_out
#endif

end subroutine




!> removes an atom from the box
subroutine remove_atom(ib, im, i, lusel, lauto)

    !use kinds_f90
    use nbrlist_module
    use atom_module, only : copy_atom, zero_atom

    !implicit none

    integer, intent(in) :: ib, im, i

    logical, intent(in) :: lusel, lauto

    integer :: j, i2, gno, gno2, typ, ml
    
    typ = cfgs(ib)%mols(im)%atms(i)%atlabel

    if(i == cfgs(ib)%mols(im)%natom) then ! this atom is the last in the cell

        ! set global no before it is deleted

        call zero_atom(cfgs(ib)%mols(im)%atms(i))

        ! get rid of this entry in the nbr list
        if(lusel) call clear_atomnbrlist(cfgs(ib)%mols(im)%nlist, i)

    else ! replace with the last atom

        i2 = cfgs(ib)%mols(im)%natom

        ! shift the neighbour list as well
        if(lusel) call swap_atomnbrlist(cfgs(ib)%mols(im)%nlist, i, i2)

        call copy_atom(cfgs(ib)%mols(im)%atms(i), cfgs(ib)%mols(im)%atms(i2))

        ! set the last entry to zero
        j = cfgs(ib)%mols(im)%natom

        if(lusel) call clear_atomnbrlist(cfgs(ib)%mols(im)%nlist, i2)

        call zero_atom(cfgs(ib)%mols(im)%atms(j))

    endif

    cfgs(ib)%mols(im)%natom = cfgs(ib)%mols(im)%natom - 1
    cfgs(ib)%number_of_atoms = cfgs(ib)%number_of_atoms - 1
    
    ml = cfgs(ib)%mols(im)%mol_label
    cfgs(ib)%mtypes(ml)%num_elems(typ) = cfgs(ib)%mtypes(ml)%num_elems(typ) - 1

    !TU: 'num_vol_sites' must be correct if Gibbs volume moves are in use, hence the
    !TU: below amendment. Note that I assume that the removed atom is a 'moving particle'
    !TU: for the sake of volume moves - a reasonable assumption
    cfgs(ib)%num_vol_sites = cfgs(ib)%num_vol_sites - 1

    cfgs(ib)%nums_elemts(typ) =  cfgs(ib)%nums_elemts(typ) - 1


end subroutine


subroutine remove_molecule(ib, imol, fail ) ! lusel, lauto)

    !use kinds_f90
    use molecule_module, only : dealloc_mol_atms_arrays, copy_molecule, zero_molecule
    use constants_module, only : uout

    !implicit none

    integer, intent(in)    :: ib, imol
    integer, intent(inout) :: fail

    integer :: ml1, ml2, im1, im2, im, jm, i, atlabel
    
    logical :: found1, found2

    fail = 0
    
    found1 = .false.
    found2 = .false.

    !TU: Amend cfgs(ib)%nums_elemts(typ) to keep correct count
    do i = 1, cfgs(ib)%mols(imol)%natom

        atlabel =  cfgs(ib)%mols(imol)%atms(i)%atlabel
        cfgs(ib)%nums_elemts(atlabel) = cfgs(ib)%nums_elemts(atlabel) - 1

    end do


    ! if the last in the array then just delete
    if(imol == cfgs(ib)%num_mols) then

        !im1 = imol
        ml1 = cfgs(ib)%mols(imol)%mol_label

        do i = cfgs(ib)%mtypes(ml1)%num_mols,1,-1
           if( imol==cfgs(ib)%mtypes(ml1)%mol_id(i) ) then 
               im = i
               found1 = .true.
               exit
           end if
        end do
        
        if( .not.found1 ) then
            fail = 240
            return
        end if
        
        if( im /= cfgs(ib)%mtypes(ml1)%num_mols ) &
            cfgs(ib)%mtypes(ml1)%mol_id(im) = cfgs(ib)%mtypes(ml1)%mol_id( cfgs(ib)%mtypes(ml1)%num_mols )

        cfgs(ib)%mtypes(ml1)%mol_id( cfgs(ib)%mtypes(ml1)%num_mols ) = 0
        cfgs(ib)%mtypes(ml1)%num_mols = cfgs(ib)%mtypes(ml1)%num_mols - 1
        
        cfgs(ib)%num_mols = cfgs(ib)%num_mols - 1
        cfgs(ib)%number_of_atoms = cfgs(ib)%number_of_atoms - cfgs(ib)%mols(imol)%natom

        ! destroy molecule

        call zero_molecule(cfgs(ib)%mols(imol))

!JG: ?There is no need to spend time allocating and deallocating for single molecule gcmc of specific type.
        call dealloc_mol_atms_arrays(cfgs(ib)%mols(imol))
!        print *,'deallocating mol arrays 1'

    else ! replace mol with the last molecule

        im1 = imol
        im2 = cfgs(ib)%num_mols
        
        ml1 = cfgs(ib)%mols(im1)%mol_label
        ml2 = cfgs(ib)%mols(im2)%mol_label
        
        if( ml1 == ml2 ) then

            do i = cfgs(ib)%mtypes(ml1)%num_mols,1,-1
               if( im1 == cfgs(ib)%mtypes(ml1)%mol_id(i) ) then
                   im = i
                   found1 = .true.
                   if( found2 ) exit
               end if
               if( im2 == cfgs(ib)%mtypes(ml1)%mol_id(i) ) then
                   jm = i
                   found2 = .true.
                   if( found1 ) exit
               end if
            end do
            
            if( .not.found1 .or. .not.found2 ) then
                fail = 241
                return
            end if
            
            !cfgs(ib)%mtypes(ml1)%mol_id(im)=im2
            cfgs(ib)%mtypes(ml2)%mol_id(jm)=im1
            
            if( im /= cfgs(ib)%mtypes(ml1)%num_mols ) &
                cfgs(ib)%mtypes(ml1)%mol_id(im) = cfgs(ib)%mtypes(ml1)%mol_id( cfgs(ib)%mtypes(ml1)%num_mols )
            
            !if( jm /= cfgs(ib)%mtypes(ml1)%num_mols ) &
            !    cfgs(ib)%mtypes(ml1)%mol_id(jm) = cfgs(ib)%mtypes(ml1)%mol_id( cfgs(ib)%mtypes(ml1)%num_mols )
        
        else

            do i = cfgs(ib)%mtypes(ml1)%num_mols,1,-1
               if( im1==cfgs(ib)%mtypes(ml1)%mol_id(i) ) then 
                   im = i
                   found1 = .true.
                   exit
               end if
            end do
            do i = cfgs(ib)%mtypes(ml2)%num_mols,1,-1
               if( im2==cfgs(ib)%mtypes(ml2)%mol_id(i) ) then
                   jm = i
                   found2 = .true.
                   exit
               end if
            end do
            
            if( .not.found1 .or. .not.found2 ) then
                fail = 242
                return
            end if
            
            !cfgs(ib)%mtypes(ml1)%mol_id(im) = im2
            cfgs(ib)%mtypes(ml2)%mol_id(jm) = im1
            
            if( im /= cfgs(ib)%mtypes(ml1)%num_mols ) &
                cfgs(ib)%mtypes(ml1)%mol_id(im) = cfgs(ib)%mtypes(ml1)%mol_id( cfgs(ib)%mtypes(ml1)%num_mols )
        end if

        cfgs(ib)%mtypes(ml1)%mol_id( cfgs(ib)%mtypes(ml1)%num_mols ) = 0
        cfgs(ib)%mtypes(ml1)%num_mols = cfgs(ib)%mtypes(ml1)%num_mols - 1

        !DA: need to adjust these first in case im1 and im2 have different no of atoms
        !DA: if we do copy_molecule first, the code eventually runs into error 203
        cfgs(ib)%num_mols = cfgs(ib)%num_mols - 1
        cfgs(ib)%number_of_atoms = cfgs(ib)%number_of_atoms - cfgs(ib)%mols(imol)%natom

        ! copy data into previous molecule
        call copy_molecule(cfgs(ib)%mols(im1), cfgs(ib)%mols(im2))

        ! destroy molecule

        call zero_molecule(cfgs(ib)%mols(im2))

!JG: ?There is no need to spend time allocating and deallocating for single molecule gcmc of specific type.
        call dealloc_mol_atms_arrays(cfgs(ib)%mols(im2))
!        print *,'deallocating mol arrays 2'

    endif
    
    call update_glob_no(ib)

    !TU: 'num_vol_sites' must be correct if Gibbs volume moves are in use, hence the
    !TU: below amendment. NOTE THAT I ASSUME THAT THE REMOVED MOLECULE IS A SINGLE
    !TU: 'MOVING PARTICLE' FOR THE SAKE OF VOLUME MOVES - a reasonable assumption if
    !TU: the molecule is rigid, but not necessarily so if the molecule is flexible or
    !TU: something else
    cfgs(ib)%num_vol_sites = cfgs(ib)%num_vol_sites - 1

end subroutine

!> updates glob_no array after mutation / removal / insertion
subroutine update_glob_no(ib)

    !implicit none

    integer, intent(in) :: ib

    integer :: im, i, idx

    idx = 0
    do im = 1, cfgs(ib)%num_mols

        do i = 1, cfgs(ib)%mols(im)%natom

            idx = idx + 1
            cfgs(ib)%mols(im)%glob_no(i) = idx

        enddo

    enddo

end subroutine

!> allocate grid for Cavity Bias GCMC
!AB: the routine must be obsolete by now - the allocations are done in build_cavitygrid()
subroutine alloc_cavitygrid0(cfg, dx0, dy0, dz0)

    use kinds_f90
    use config_type

    implicit none

        !> configuration(s) container - array in the case of several replicas being used
    type (config), intent(inout) :: cfg

        !> 3D grid bin sizes (dX, dY, dZ)
    real (kind = wp), intent(in) :: dx0, dy0, dz0

    ! number of 3D grid points
    integer :: nx, ny, nz

    allocate(cfg%nrcav(3))
    allocate(cfg%drcav(3))

    nx = int(cfg%vec%latvector(1,1)/dx0)
    ny = int(cfg%vec%latvector(2,2)/dy0)
    nz = int(cfg%vec%latvector(3,3)/dz0)

!AB: make sure the grid is fully filling the volume

    cfg%nrcav(1) = nx
    cfg%nrcav(2) = ny
    cfg%nrcav(3) = nz

    cfg%drcav(1) = cfg%vec%latvector(1,1)/real(nx,wp)
    cfg%drcav(2) = cfg%vec%latvector(2,2)/real(ny,wp)
    cfg%drcav(3) = cfg%vec%latvector(3,3)/real(nz,wp)

    allocate(cfg%cavgrid(3,nx*ny*nz))
    allocate(cfg%cavity(nx*ny*nz))

end subroutine

!> creates the cavity bias grid
!AB: the routine must be obsolete by now - see below for the new/tested routine
subroutine build_cavitygrid0(ib, dx, dy, dz)

    !use kinds_f90

    !implicit none

    integer, intent(in) :: ib

    real(kind = wp), intent(in) :: dx, dy, dz

    integer :: ni, i, j , k
    integer :: cavx, cavy, cavz

    real(kind = wp) :: gx, gy, gz

    cavx = int(cfgs(ib)%vec%latvector(1,1) / dx)
    cavy = int(cfgs(ib)%vec%latvector(2,2) / dy)
    cavz = int(cfgs(ib)%vec%latvector(3,3) / dz)

    ! zero accumulators
    ni = 0 

    !AB: The grid displacements above (dx,dy,dz) are certainly 
    !    assumed to be given in Cartesian coordinates (Angstrom)
    !    Similarly, the grid coordinates in create_cavitylist(ib, crad)
    !    are assumed to be Cartesian too
    !
    !AB: Why the grid starts at (-0.5, -0.5, -0.5) ???
    !AB: Is it supposed to be defined in fractional coordinates ???
    !
    !AB: If not in fractional coordinates, shouldn't the grid start at 
    !    (-0.5*cfgs(ib)%vec%latvector(1,1),
    !     -0.5*cfgs(ib)%vec%latvector(2,2),
    !     -0.5*cfgs(ib)%vec%latvector(3,3)) ???
    !
    !AB: It seems also that the grid defined below 
    !    could only be correct for orthorhombic cells
    
    gx = -0.5

    do i = 1, cavx

        gy = -0.5

        do j = 1, cavy

            gz = -0.5

            do k = 1, cavz

                ni = ni + 1

                cfgs(ib)%cavgrid(1,ni) = gx
                cfgs(ib)%cavgrid(2,ni) = gy
                cfgs(ib)%cavgrid(3,ni) = gz

                gz = gz + dy !AB: should add dz

            enddo

            gy = gy + dy

        enddo

        gx = gx + dx

    enddo

    if( ni /= cavx*cavy*cavz ) call error(136)

    ncavity_grid(ib) = ni

end subroutine

!> creates the cavity bias grid
!AB: new, fully reworked and tested variant for the cavity grid/mesh
subroutine build_cavitygrid(ib, dx0, dy0, dz0)

    !use kinds_f90
    use constants_module, only : HALF, DZERO, uout
    use arrays_module, only : reallocate_quiet
    use parse_module, only : int_2_word

    !implicit none

    integer, intent(in)         :: ib
    real(kind = wp), intent(in) :: dx0, dy0, dz0

    integer, save :: nc0 = 0
    integer, save :: nmols0  = 0
    integer, save :: natoms0 = 0
    integer, save :: ntypes0 = 0
    integer, save :: melems0 = 0

    integer :: fail(0:3)
    integer :: nc, ni, i, j , k, ign
    integer :: cavx, cavy, cavz

    real(kind = wp) :: dxcav, dycav, dzcav
    real(kind = wp) :: gx, gy, gz

    if( natoms0 == 0 ) natoms0 = cfgs(ib)%number_of_atoms

    cavx = int(cfgs(ib)%vec%latvector(1,1)/dx0)
    cavy = int(cfgs(ib)%vec%latvector(2,2)/dy0)
    cavz = int(cfgs(ib)%vec%latvector(3,3)/dz0)

    if( cavx < 2 .or. cavy < 2 .or. cavz < 2 ) &
        call cry(uout,'', &
             "ERROR: build_cavitygrid():: too few points on cavity grid ("//&
             trim(int_2_word(cavx))//", "//&
             trim(int_2_word(cavy))//", "//&
             trim(int_2_word(cavz))//", "//&
             ") !!!",999)

!AB: make sure the grid is fully filling the volume

    dxcav = cfgs(ib)%vec%latvector(1,1)/real(cavx,wp)
    dycav = cfgs(ib)%vec%latvector(2,2)/real(cavy,wp)
    dzcav = cfgs(ib)%vec%latvector(3,3)/real(cavz,wp)

    nc = cavx*cavy*cavz

    if( nc /= nc0 ) then

        fail = 0

        call reallocate_quiet(cfgs(ib)%nrcav, 3, fail(1))
        call reallocate_quiet(cfgs(ib)%drcav, 3, fail(2))
        call reallocate_quiet(cfgs(ib)%cavgrid, 3, nc, fail(3))

        if( any(fail > 0) ) &
            call cry(uout,'', &
                 "ERROR: build_cavitygrid() failed to (re)allocate arrays for cavity bins !!!",999)

        fail = 0

        call reallocate_quiet(cfgs(ib)%cavity,  nc, fail(1))
        call reallocate_quiet(cfgs(ib)%cavids,  nc, fail(2))
        call reallocate_quiet(cavatt, nc, fail(3))

        if( any(fail > 0) ) &
            call cry(uout,'', &
                 "ERROR: build_cavitygrid() failed to (re)allocate arrays for cavity grid !!!",999)

        nc0 = nc

        fail = 0

!AB: exhaustive allocation for P_tot(max_atoms):
!
!        call reallocate_quiet(cavatt, nc, fail(1))
!        call reallocate_quiet(cfgs(ib)%pcavity_tot, 2, cfgs(ib)%maxno_of_atoms, fail(2))
!
!        nmols0  = size(cfgs(ib)%mtypes(:))
!        ntypes0 = size(cfgs(ib)%mtypes(1)%num_elems)
!        melems0 = 0
!        do j = 1, ntypes0
!            do i = 1, nmols0
!                !melems0 = melems0 + cfgs(ib)%mtypes(i)%max_mols*cfgs(ib)%mtypes(i)%num_elems(j)
!                melems0 = max(melems0,cfgs(ib)%mtypes(i)%max_mols*cfgs(ib)%mtypes(i)%num_elems(j))
!                natoms0 = natoms0 + melems0
!            enddo
!        enddo
!        write(uout,'(/,4(a,i8),a,i3,a,/)')"build_cavitygrid(ib) : max number of elements ",ntypes0," * ",melems0, &
!                                " => ",natoms0," =?= ",cfgs(ib)%maxno_of_atoms," (",ib,")"
!
!AB: exhaustive allocation for P_typ(max_atoms,num_elems):
!
!        ntypes0 = size(cfgs(ib)%mtypes(1)%num_elems)
!        call reallocate_quiet(cfgs(ib)%pcavity_typ, cfgs(ib)%maxno_of_atoms, ntypes0, fail(3))
!
!        if( any(fail > 0) ) &
!            call cry(uout,'', &
!                    "ERROR: build_cavitygrid() failed to (re)allocate arrays for cavity probs !!!",999)
!
!        cfgs(ib)%pcavity_tot = 0.0_wp
!        cfgs(ib)%pcavity_typ = 0.0_wp

    endif

    cfgs(ib)%nrcav(1) = cavx
    cfgs(ib)%nrcav(2) = cavy
    cfgs(ib)%nrcav(3) = cavz

    cfgs(ib)%drcav(1) = dxcav
    cfgs(ib)%drcav(2) = dycav
    cfgs(ib)%drcav(3) = dzcav

    cfgs(ib)%cavgrid(:,:) = 0.0_wp
    cfgs(ib)%cavity(:)    = 0
    cfgs(ib)%cavids(:)    = 0
    cavatt(:) = 0

    !AB: number of grid points
    ni = 0 

    gx = -HALF*(cfgs(ib)%vec%latvector(1,1) - dxcav)
    !gx = -HALF*cfgs(ib)%vec%latvector(1,1)

    do i = 1, cavx

        gy = -HALF*(cfgs(ib)%vec%latvector(2,2) - dycav)
        !gy = -HALF*cfgs(ib)%vec%latvector(2,2)

        do j = 1, cavy

            gz = -HALF*(cfgs(ib)%vec%latvector(3,3) - dzcav)
            !gz = -HALF*cfgs(ib)%vec%latvector(3,3)

            do k = 1, cavz

                ni = ni + 1

                cfgs(ib)%cavgrid(1,ni) = gx
                cfgs(ib)%cavgrid(2,ni) = gy
                cfgs(ib)%cavgrid(3,ni) = gz

                gz = gz + dzcav

            enddo

            gy = gy + dycav

        enddo

        gx = gx + dxcav

    enddo

    if (ni /= nc) call error(136)

    ncavity_grid(ib) = ni

end subroutine build_cavitygrid

!> constructs cavity bias list
!AB: new, fully reworked and tested variant for the cavity look-up list(s)
subroutine create_cavitylist(ib, rc) !, by_molcom)

    !use kinds_f90
    use constants_module, only : MAXINT1, HALF, DZERO, FZERO, uout
    use slit_module, only : in_slit

    !implicit none

    integer, intent(in)         :: ib
    real(kind = wp), intent(in) :: rc
    !logical, intent(in)         :: by_molcom

    logical, save :: new = .true.
    integer, save :: ucavb=137
    integer, save :: ucavf=147

    integer jm, j, i, ni, nfree0

    real(kind = wp) :: xi, yi, zi, xj, yj, zj, rc2 !, rsq
    real(kind = wp) :: xx, yy, zz, xx2, yy2, zz2

    logical, save :: is_gcmol = .true.

    !is_GCMC = ( job%gcmcatom .or. job%gcmcmol .or. job%gibbsatomtran .or. job%gibbsmoltran & 
    !                         .or. job%gibbs_indvol .or. job%gibbsatomexch .or. job%gibbsmolexch )

    if( new ) then
        open(ucavb, file='CAVITYB0')
        open(ucavf, file='CAVITYF0')
        write(ucavb,"(a)")'# Initial cavity grid (busy/free)'
        write(ucavf,"(a)")'# Initial cavity grid (free only)'
        new = .false.
        nfree0 = 0
    else
        open(ucavb, file='CAVITYB')
        open(ucavf, file='CAVITYF')
        write(ucavb,"(a)")'# Latest cavity grid (busy/free)'
        write(ucavf,"(a)")'# Latest cavity grid (free only)'
        nfree0 = ncavity_free(ib)
    endif

    rc2 = rc*rc

    ncavity_free(ib)   = 0
    cfgs(ib)%cavity(:) = 0
    cfgs(ib)%cavids(:) = 0
    cavatt(:) = 0

    ni = 0

    !AB: search for and count occupancies, as well as free cavity sites

    grid: &
    do i = 1, ncavity_grid(ib)

        xi = cfgs(ib)%cavgrid(1,i)
        yi = cfgs(ib)%cavgrid(2,i)
        zi = cfgs(ib)%cavgrid(3,i)

        !AB: make sure no free cavity sites are found outside the 'slab'
        if( atom_outside(zi, cfgs(ib)%vec%latvector(3,3)) ) then
            cfgs(ib)%cavity(i) = MAXINT1/2
            cycle grid
        end if

        mols: &
        do jm = 1, cfgs(ib)%num_mols

!            if( by_molcom ) then
!            else
!            endif

            atoms: &
            do j = 1, cfgs(ib)%mols(jm)%natom

!AB: skip fictitious atoms if any (set mass to zero)
                if( cfgs(ib)%mols(jm)%atms(j)%mass < FZERO ) cycle atoms

                xj = cfgs(ib)%mols(jm)%atms(j)%rpos(1)
                yj = cfgs(ib)%mols(jm)%atms(j)%rpos(2)
                zj = cfgs(ib)%mols(jm)%atms(j)%rpos(3)

                xx = xi - xj
                yy = yi - yj
                zz = zi - zj

                call pbc_cart_pos_calc(ib, xx, yy, zz)

                !if( xx*xx + yy*yy + zz*zz - rc2 > DZERO ) cycle atoms

                xx2 = xx*xx
                if( xx2 - rc2 > DZERO ) cycle atoms
                !if( xx2 > rc2 ) cycle atoms

                yy2 = yy*yy
                if( yy2 - rc2 > DZERO ) cycle atoms
                !if( yy2 > rc2 ) cycle atoms

                zz2 = zz*zz
                if( zz2 - rc2 > DZERO ) cycle atoms
                !if( zz2 > rc2 ) cycle atoms

                if( xx2 + yy2 + zz2 - rc2 > DZERO ) cycle atoms

                !if( xx2 + yy2 + zz2 - rc2 < DZERO ) &
                    cfgs(ib)%cavity(i) = cfgs(ib)%cavity(i) + 1

                    !rsq = xx2 + yy2 + zz2 
                    !if( rc2-rsq > DZERO ) cfgs(ib)%cavity(i) = cfgs(ib)%cavity(i) + 1
                    !if( rc2-rsq < DZERO ) cycle

                    !!if( rsq < rc2 ) then
                    !!    cfgs(ib)%cavity(i) = cfgs(ib)%cavity(i) + 1
                    !!endif

                !endif

             enddo atoms

        enddo mols

        if( cfgs(ib)%cavity(i) == 0 ) then

            ni = ni+1

            cfgs(ib)%cavids(ni) = i

            write(ucavb,"(3(a,f15.7),a,3i8)")'free: ',xi,' ',yi,' ',zi,&
                                             '  # ',cfgs(ib)%cavity(i),i,ni
            write(ucavf,"(3(a,f15.7),a,3i8)")'free: ',xi,' ',yi,' ',zi,&
                                             '  # ',cfgs(ib)%cavity(i),i,ni
        else
            write(ucavb,"(3(a,f15.7),a,3i8)")'busy: ',xi,' ',yi,' ',zi,&
                                             '  # ',cfgs(ib)%cavity(i),i,ncavity_grid(ib)-ni
        endif

    enddo grid

    close(ucavb)
    close(ucavf)
    
    nfreed = 0
    ntaken = 0

    ncavity_free(ib) = ni

    ni = count( cfgs(ib)%cavids(:) > 0 )

    if( ni /= ncavity_free(ib) ) &
        call cry(uout,'', &
                "ERROR: create_cavitylist() - count(s) of the cavity sites inconsistent !!!",999)

    write(uout,'(/,a,i3,3(a,i8),a,3f12.5,a,/)')&
         "free/total cavity sites (",ib,") : ",nfree0," =?= ",ncavity_free(ib)," / ",ncavity_grid(ib),&
         "; (dx,dy,dz) = (",cfgs(ib)%drcav(1),cfgs(ib)%drcav(2),cfgs(ib)%drcav(3),")"

    if( in_slit ) then

        ni = count( cfgs(ib)%cavids(:) == MAXINT1/2 )

        write(uout,'(a,i8,/)')"occup. cavity sites outside slit (ib)  : ",ni

    endif

!    cfgs(ib)%pcavity_tot(1,cfgs(ib)%number_of_atoms) = ncavity_grid(ib)
!    cfgs(ib)%pcavity_tot(2,cfgs(ib)%number_of_atoms) = ncavity_free(ib)

    flush(uout)

end subroutine create_cavitylist

!AB: *** new - highly optimised, debugged and tested - routines for partially updating the cavity list(s) ***

!AB: new, fully reworked and tested variant for updating the cavity list(s) upon atom insertion
subroutine calc_pcav_insert(ib, ig, crad, pos, prob, is_first)

    !use kinds_f90
    use constants_module, only : HALF, DHALF, DZERO, FZERO, uout
    use parse_module, only : int_2_word
    use random_module, only : duni

    !implicit none

    integer, intent(in)            :: ib, ig
    real(kind = wp), intent(in)    :: crad
    real(kind = wp), intent(inout) :: pos(:), prob
    logical, intent(in)            :: is_first

    integer :: igx, igy, igz, ign, ipm
    integer :: ixmin,ixmax,iymin,iymax,izmin,izmax
    integer :: jgx, jgy, jgz, nfree, k, ic

    real(kind = wp) :: rc2, dxcav, dycav, dzcav, xhalf, yhalf, zhalf
    real(kind = wp) :: ix, iy, iz, gx, gy, gz, xx, yy, zz, xx2, yy2, zz2 !, rsq

#ifdef DEBUG
!debug<
    logical :: found

    integer, save :: natt = 0
    integer :: nins !, nfree

    nins  = 0
    nfree = 0
!debug>
#endif

    !prob = -1.0_wp
    if( ncavity_free(ib) < 1 ) return

!AB: for detailed balance the probability is to be estimated before atom insertion
    prob = real(ncavity_free(ib),wp) / real(ncavity_grid(ib),wp)

#ifdef DEBUG
!debug<
    natt  = natt + 1
!debug>
#endif

    rc2 = crad * crad

    dxcav = cfgs(ib)%drcav(1)
    dycav = cfgs(ib)%drcav(2)
    dzcav = cfgs(ib)%drcav(3)

    ix = pos(1)
    iy = pos(2)
    iz = pos(3)

!AB: need to apply PBC to map onto the cavity mesh
    !if( .not.cfgs(ib)%vec%is_orthogonal ) &
        call pbc_cart2frac_pos_calc(ib, ix, iy, iz)

!AB: optimisation - update cavities in box [ix+/-dxcav, iy+/-dycav, iz+/-dzcav]

    xhalf = cfgs(ib)%vec%latvector(1,1)*HALF
    yhalf = cfgs(ib)%vec%latvector(2,2)*HALF
    zhalf = cfgs(ib)%vec%latvector(3,3)*HALF

    igx = int((ix+xhalf)/dxcav)+1
    igy = int((iy+yhalf)/dycav)+1
    igz = int((iz+zhalf)/dzcav)+1

    xhalf = xhalf-dxcav*HALF
    yhalf = yhalf-dycav*HALF
    zhalf = zhalf-dzcav*HALF

    ign = igz + ((igy-1) + (igx-1)*cfgs(ib)%nrcav(2))*cfgs(ib)%nrcav(3)

#ifdef DEBUG
!debug<
    if( ig > 0 ) then 

        if( ign /= ig ) &
            write(*,'(1x,2(a,i5),a,3f15.7,a)') &
            "calc_pcav_insert(1):: recalculated and given indices differ: ",&
            ign," =/= ",ig," for site (",ix,iy,iz,") !!!"

        if( myjob % cavity_mesh ) then

            if( abs( cfgs(ib)%cavgrid(1,ign) - ix ) > dxcav .or. &
                abs( cfgs(ib)%cavgrid(2,ign) - iy ) > dycav .or. &
                abs( cfgs(ib)%cavgrid(3,ign) - iz ) > dzcav ) &
                write(*,'(1x,2(a,3f15.7),a,i6,a)')&
                "calc_pcav_insert(1):: coordinates for new atom "//&
                "given and taken in 1D differ: (",ix,iy,iz,") =?= (",&
                cfgs(ib)%cavgrid(1,ig),cfgs(ib)%cavgrid(2,ig),cfgs(ib)%cavgrid(3,ig),&
                ") for 1D index ",ig," !!!"

        else

            if( abs( cfgs(ib)%cavgrid(1,ign) - ix ) > DZERO .or. &
                abs( cfgs(ib)%cavgrid(2,ign) - iy ) > DZERO .or. &
                abs( cfgs(ib)%cavgrid(3,ign) - iz ) > DZERO) &
                write(*,'(1x,2(a,3f15.7),a,i6,a)')&
                "calc_pcav_insert(1):: coordinates for new atom "//&
                "given and taken in 1D differ: (",ix,iy,iz,") =?= (",&
                cfgs(ib)%cavgrid(1,ig),cfgs(ib)%cavgrid(2,ig),cfgs(ib)%cavgrid(3,ig),&
                ") for 1D index ",ig," !!!"

        endif

!        if( abs( cfgs(ib)%cavgrid(1,ign) - ix ) > dxcav .or. &
!            abs( cfgs(ib)%cavgrid(2,ign) - iy ) > dycav .or. &
!            abs( cfgs(ib)%cavgrid(3,ign) - iz ) > dzcav ) &
!            write(*,'(1x,2(a,3f15.7),a,i6,a)')&
!            "calc_pcav_insert(2):: coordinates for new atom "//&
!            "given and taken in 1D differ: (",ix,iy,iz,") =?= (",&
!            cfgs(ib)%cavgrid(1,ig),cfgs(ib)%cavgrid(2,ig),cfgs(ib)%cavgrid(3,ig),&
!            ") for 1D index ",ig," !!!"

    endif
!debug>
#endif

    ipm   = int(crad/dxcav)+1
    ixmin = igx-ipm
    ixmax = igx+ipm

    ipm   = int(crad/dycav)+1
    iymin = igy-ipm
    iymax = igy+ipm

    ipm   = int(crad/dzcav)+1
    izmin = igz-ipm
    izmax = igz+ipm

    if( is_first ) then
        cavatt(1) = 1
        ntaken = 0

        !if( ig < 0 ) then
        !    ntaken = 1
        !    cavatt(1) = 2
        !    cavatt(2) = -ig
        !endif
    endif

    do jgx = ixmin, ixmax

        igx = jgx
        if( igx < 1 ) then
            igx = igx+cfgs(ib)%nrcav(1)
        elseif( igx > cfgs(ib)%nrcav(1) ) then
            igx = igx-cfgs(ib)%nrcav(1)
        endif

        gx = real(igx-1,wp)*dxcav - xhalf

        xx = (gx - ix)/cfgs(ib)%vec%latvector(1,1)
        xx = (xx - nint(xx))*cfgs(ib)%vec%latvector(1,1)

        xx2 = xx*xx
        if( xx2 - rc2 > DZERO ) cycle
        !if( xx2 > rc2 ) cycle

        do jgy = iymin, iymax

            igy = jgy
            if( igy < 1 ) then
                igy = igy+cfgs(ib)%nrcav(2)
            elseif( igy > cfgs(ib)%nrcav(2) ) then
                igy = igy-cfgs(ib)%nrcav(2)
            endif

            gy = real(igy-1,wp)*dycav - yhalf

            yy = (gy - iy)/cfgs(ib)%vec%latvector(2,2)
            yy = (yy - nint(yy))*cfgs(ib)%vec%latvector(2,2)

            yy2 = yy*yy
            if( yy2 - rc2 > DZERO ) cycle
            !if( yy2 > rc2 ) cycle

            do jgz = izmin, izmax

                igz = jgz
                if( igz < 1 ) then
                    igz = igz+cfgs(ib)%nrcav(3)
                elseif( igz > cfgs(ib)%nrcav(3) ) then
                    igz = igz-cfgs(ib)%nrcav(3)
                endif

                gz = real(igz-1,wp)*dzcav - zhalf

                zz = (gz - iz)/cfgs(ib)%vec%latvector(3,3)
                zz = (zz - nint(zz))*cfgs(ib)%vec%latvector(3,3)

                zz2 = zz*zz
                if( zz2 - rc2 > DZERO ) cycle
                !!if( zz2 > rc2 ) cycle

                if( xx2 + yy2 + zz2 - rc2 > DZERO ) cycle

                ign = igz + ((igy-1) + (igx-1)*cfgs(ib)%nrcav(2))*cfgs(ib)%nrcav(3)

                !rsq = xx2 + yy2 + zz2 
                !if( rc2-rsq < DZERO ) cycle

!!!#ifdef DEBUG
!debug<
                if( abs( cfgs(ib)%cavgrid(1,ign) - gx ) > DZERO .or. &
                    abs( cfgs(ib)%cavgrid(2,ign) - gy ) > DZERO .or. &
                    abs( cfgs(ib)%cavgrid(3,ign) - gz ) > DZERO ) &
                    write(*,'(1x,2(a,3f15.7),a,i6,a)')&
                    "calc_pcav_insert(3):: coordinates for new atom "//&
                    "given and taken in 1D differ: (",gx,gy,gz,") =?= (",&
                    cfgs(ib)%cavgrid(1,ign),cfgs(ib)%cavgrid(2,ign),cfgs(ib)%cavgrid(3,ign),&
                    ") for 1D index ",ign," !!!"
!debug>
!!!#endif

#ifndef DEBUG
!nondebug<
                cfgs(ib)%cavity(ign) = cfgs(ib)%cavity(ign) + 1
                if( cfgs(ib)%cavity(ign) == 1 ) ntaken = ntaken+1

                cavatt(1) = cavatt(1)+1
                cavatt(cavatt(1)) = ign

                cycle
!nondebug>
#endif

#ifdef DEBUG
!debug<
                !if( rc2-rsq > DZERO ) then
                !!if( rsq < rc2 ) then

                    found  = .false.
                    do k = 1,ncavity_free(ib)

                       if( cfgs(ib)%cavids(k) == ign ) then

                           found = .true.
                           ic = k

                           exit
                       end if

                    enddo

                    if( found ) then
                        if( cfgs(ib)%cavity(ign) > 0 ) then
                            call cry(uout,'', &
                                 "ERROR: calc_pcav_insert() - "//&
                                 "cavity count > 0 at new atom insertion "// &
                                 "(yet found amongst free cavities) !!!",999)
                        endif

                        nfree  = nfree+1
                        
                        cfgs(ib)%cavids(ic) = cfgs(ib)%cavids(ncavity_free(ib))
                        cfgs(ib)%cavids(ncavity_free(ib)) = 0
                        ncavity_free(ib) = ncavity_free(ib) - 1
                    end if

                    cfgs(ib)%cavity(ign) = cfgs(ib)%cavity(ign) + 1
                    if( cfgs(ib)%cavity(ign) == 1 ) ntaken = ntaken+1

                    cavatt(1) = cavatt(1)+1
                    cavatt(cavatt(1)) = ign

                !end if !( rsq < rc2 ) then
!debug>
#endif
            enddo

        enddo

    enddo

#ifdef DEBUG
!debug<
    if( ncavity_free(ib) < 1 ) &
        call cry(uout,'', &
             "WARNING: calc_pcav_insert() - "//&
             "free cavity count < 1 upon inserting new atom !!!",0)

    if( nfree < 1 ) &
    !if( ntaken < 1 ) &
        call cry(uout,'', &
             "WARNING: calc_pcav_insert() - "//&
             "no existing free cavity occupied upon inserting new atom !!!",0)
             !"ERROR: calc_pcav_insert() - "//&
             !"no existing free cavity occupied upon inserting new atom !!!",999)


    if( ig > 0 .and. ntaken /= nfree ) &
        call cry(uout,'', &
             "WARNING: calc_pcav_insert() - "//&
             "taken cavity count does not equal "//&
             "the free cavity difference upon inserting new atom !!!",0)
!debug>
#endif

#ifndef DEBUG
!nondebug<

    cfgs(ib)%cavids(:) = 0

    nfree = 0
    do k = 1, ncavity_grid(ib)

        if( cfgs(ib)%cavity(k) > 0 ) cycle

        nfree = nfree+1

        cfgs(ib)%cavids(nfree) = k

    enddo

    if( ig > 0 .and. ncavity_free(ib) <= nfree ) &
        call cry(uout,'', &
             "WARNING: calc_pcav_insert() - "//&
             "free cavity count did not reduce upon inserting new atom !!!",0)

    if( ig > 0 .and. ntaken /= (ncavity_free(ib) - nfree) ) &
        call cry(uout,'', &
             "WARNING: calc_pcav_insert() - "//&
             "taken cavity count is not equal to the free cavity difference upon inserting new atom !!!",0)

    ncavity_free(ib) = nfree

!nondebug>
#endif

end subroutine calc_pcav_insert

!AB: new, fully reworked and tested variant for updating the cavity list(s) upon molecule removal
subroutine calc_pcav_insert_mol(ib, ig0, im, crad, prob, by_com)

    !use kinds_f90
    use constants_module, only : FZERO, uout

    !implicit none

    integer, intent(in)            :: ib, ig0, im
    real(kind = wp), intent(in)    :: crad
    real(kind = wp), intent(inout) :: prob
    logical, intent(in)            :: by_com

    integer :: ig, nc, ia, n, k

    logical :: is_first != .true.

    real(kind = wp) :: pcav

    prob = 1.0_wp

!    ig = ig0
    nc = cfgs(ib)%cavity(ig0)
!    ig =-ig

    is_first = .true.

!    call calc_pcav_insert(ib, ig0, myjob%cavity_radius, &
!         cfgs(ib)%mols(im)%rcom, prob, is_first)
!
!    is_first = .false.

!    ig = 0

        do ia = 1, cfgs(ib)%mols(im)%natom

!AB: skip fictitious atoms if any (set mass to zero)
            if( cfgs(ib)%mols(im)%atms(ia)%mass < FZERO ) cycle

!AB: check/update the cavity grid and pcav ...

            pcav = 1.0_wp

            call calc_pcav_insert(ib, 0, myjob%cavity_radius, &
                 cfgs(ib)%mols(im)%atms(ia)%rpos, pcav, is_first)

            prob = pcav !prob*pcav

            is_first = .false.

!            n = cavatt(1)
!            write(uout,*)"calc_pcav_insert_mol(",im,ia,") - affected cavities: ",cavatt(1),ntaken

!    do k = 2, cavatt(1)
!        if( cfgs(ib)%cavity(cavatt(k)) == 1 ) &
!            write(uout,*)"# ",cavatt(k)," -> ",cfgs(ib)%cavity(cavatt(k))
!    enddo

        enddo

!        n = cavatt(1)
!        write(uout,*)"calc_pcav_insert_mol(",im,") - affected cavities: ",cavatt(1:n)

    !enddo

    if( cfgs(ib)%cavity(ig0) - nc < 1 ) &
        call cry(uout,'', &
         "WARNING: calc_pcav_insert_mol() - "//&
         "occupied cavity count for COM did not increase upon inserting new molecule !!!",0)

end subroutine calc_pcav_insert_mol

!AB: new, fully reworked and tested variant for updating the cavity list(s) upon removal of an inserted atom
subroutine reset_pcav_insert(ib, ig, crad, pos)

    !use kinds_f90
    use constants_module, only : HALF, DZERO, uout
    use parse_module, only : int_2_word

    !implicit none

    integer, intent(in)            :: ib, ig
    real(kind = wp), intent(in)    :: crad
    real(kind = wp), intent(inout) :: pos(:)

    integer :: k, ign, nfree, ncmax, ncmin, n

#ifdef DEBUG
!debug<
    integer :: igx, igy, igz, ipm
    integer :: ixmin,ixmax,iymin,iymax,izmin,izmax
    integer :: jgx, jgy, jgz

    real(kind = wp) :: rc2, dxcav, dycav, dzcav, xhalf, yhalf, zhalf
    real(kind = wp) :: ix, iy, iz, gx, gy, gz, xx, yy, zz, xx2, yy2, zz2 !, rsq

    logical :: found

    integer, save :: natt = 0
    integer :: nrem

    natt = natt + 1
    nrem = 0

!debug>
#endif

    nfreed = 0

    ncmax = 0
    ncmin = ncavity_grid(ib)

#ifndef DEBUG
!nondebug<

!    n = cavatt(1)
!    write(uout,*)"reset_pcav_insert(..) - affected cavities: ",cavatt(1),ntaken

    do k = 2, cavatt(1)
        ign = cavatt(k)

!        if( cfgs(ib)%cavity(cavatt(k)) == 1 ) &
!            write(uout,*)"# ",cavatt(k)," -> ",cfgs(ib)%cavity(cavatt(k))

        ncmax = max(ncmax,cfgs(ib)%cavity(ign))
        ncmin = min(ncmin,cfgs(ib)%cavity(ign))

        cfgs(ib)%cavity(ign) = cfgs(ib)%cavity(ign) - 1

        if( cfgs(ib)%cavity(ign) == 0 ) then

            ncavity_free(ib) = ncavity_free(ib) + 1

            cfgs(ib)%cavids(ncavity_free(ib)) = ign
    
            nfreed = nfreed + 1

        endif

        cavatt(k) = 0

!        if( cfgs(ib)%cavity(cavatt(k)) == 0 ) &
!            write(uout,*)"# ",cavatt(k)," -> ",cfgs(ib)%cavity(cavatt(k))

    enddo

    if( nfreed < 1 ) &
        call cry(uout,'', &
             "WARNING: reset_pcav_insert() - "//&
             "number of freed cavitities < 1 after inserted atom removal "//&
             " (cavities affected: "//trim(int_2_word(cavatt(1)))//&
             " min/max = "//trim(int_2_word(ncmin))//" / "//trim(int_2_word(ncmax))//") !!!",0)
             !") !!!",0)

    if( nfreed /= ntaken ) &
        call cry(uout,'', &
             "WARNING: reset_pcav_insert() - "//&
             "number of freed cavitities "//trim(int_2_word(nfreed))//" =?= "//trim(int_2_word(ntaken))//&
             " taken (after inserted atom removal) !!!",0)
             !"number of freed cavitities is not equal to that taken (after inserted atom removal)"//&

             !"ERROR: reset_pcav_insert() - "//&
             !"cavity count > 1 after new atom removal !!!",999)

    return

!nondebug>
#endif

#ifdef DEBUG
!debug<

    rc2 = crad * crad

    dxcav = cfgs(ib)%drcav(1)
    dycav = cfgs(ib)%drcav(2)
    dzcav = cfgs(ib)%drcav(3)

    ix = pos(1)
    iy = pos(2)
    iz = pos(3)

!AB: need to apply PBC to map onto the cavity mesh
    !if( .not.cfgs(ib)%vec%is_orthogonal ) &
        call pbc_cart2frac_pos_calc(ib, ix, iy, iz)

    if( myjob % cavity_mesh ) then
        if( abs( cfgs(ib)%cavgrid(1,ig) - ix ) > dxcav .or. &
            abs( cfgs(ib)%cavgrid(2,ig) - iy ) > dycav .or. &
            abs( cfgs(ib)%cavgrid(3,ig) - iz ) > dzcav ) &
            write(*,'(1x,2(a,3f15.7),a,i6,a)')&
            "reset_pcav_insert(1):: coordinates for new atom given and taken in 1D differ: (",ix,iy,iz,&
            ") =?= (",cfgs(ib)%cavgrid(1,ig),cfgs(ib)%cavgrid(2,ig),cfgs(ib)%cavgrid(3,ig),&
            ") for 1D index ",ig," !!!"
    else
        if( abs( cfgs(ib)%cavgrid(1,ig) - ix ) > DZERO .or. &
            abs( cfgs(ib)%cavgrid(2,ig) - iy ) > DZERO .or. &
            abs( cfgs(ib)%cavgrid(3,ig) - iz ) > DZERO) &
            write(*,'(1x,2(a,3f15.7),a,i6,a)')&
            "reset_pcav_insert(1):: coordinates for new atom given and taken in 1D differ: (",ix,iy,iz,&
            ") =?= (",cfgs(ib)%cavgrid(1,ig),cfgs(ib)%cavgrid(2,ig),cfgs(ib)%cavgrid(3,ig),&
            ") for 1D index ",ig," !!!"
    endif

!AB: optimisation - update cavities in box [ix+/-cavx, iy+/-cavy, iz+/-cavz]

    xhalf = cfgs(ib)%vec%latvector(1,1)*HALF
    yhalf = cfgs(ib)%vec%latvector(2,2)*HALF
    zhalf = cfgs(ib)%vec%latvector(3,3)*HALF

    igx = int((ix+xhalf)/dxcav)+1
    igy = int((iy+yhalf)/dycav)+1
    igz = int((iz+zhalf)/dzcav)+1

    xhalf = xhalf-dxcav*HALF
    yhalf = yhalf-dycav*HALF
    zhalf = zhalf-dzcav*HALF

    ign = igz + ((igy-1) + (igx-1)*cfgs(ib)%nrcav(2))*cfgs(ib)%nrcav(3)

!#ifdef DEBUG
!debug<
    if( ign /= ig ) &
        write(*,'(1x,2(a,i5),a,3f15.7,a)')&
        "reset_pcav_insert(1):: recalculated and given indices differ: ",&
        ign," =/= ",ig," for site (",ix,iy,iz,") !!!"

    if( abs( cfgs(ib)%cavgrid(1,ign) - ix ) > dxcav .or. &
        abs( cfgs(ib)%cavgrid(2,ign) - iy ) > dycav .or. &
        abs( cfgs(ib)%cavgrid(3,ign) - iz ) > dzcav ) &
        write(*,'(1x,2(a,3f15.7),a,i6,a)')&
        "reset_pcav_insert(2):: coordinates for new atom given and taken in 1D differ: (",ix,iy,iz,&
        ") =?= (",cfgs(ib)%cavgrid(1,ig),cfgs(ib)%cavgrid(2,ig),cfgs(ib)%cavgrid(3,ig),&
        ") for 1D index ",ig," !!!"
!debug>
!#endif

    ipm   = int(crad/dxcav)+1
    ixmin = igx-ipm
    ixmax = igx+ipm

    ipm   = int(crad/dycav)+1
    iymin = igy-ipm
    iymax = igy+ipm

    ipm   = int(crad/dzcav)+1
    izmin = igz-ipm
    izmax = igz+ipm

    do jgx = ixmin, ixmax

        igx = jgx
        if( igx < 1 ) then
            igx = igx+cfgs(ib)%nrcav(1)
        elseif( igx > cfgs(ib)%nrcav(1) ) then
            igx = igx-cfgs(ib)%nrcav(1)
        endif

        gx = real(igx-1,wp)*dxcav - xhalf

        xx = (gx - ix)/cfgs(ib)%vec%latvector(1,1)
        xx = (xx - nint(xx))*cfgs(ib)%vec%latvector(1,1)
                
        xx2 = xx*xx
        if( xx2 - rc2 > DZERO ) cycle
        !if( xx2 > rc2 ) cycle

        do jgy = iymin, iymax

            igy = jgy
            if( igy < 1 ) then
                igy = igy+cfgs(ib)%nrcav(2)
            elseif( igy > cfgs(ib)%nrcav(2) ) then
                igy = igy-cfgs(ib)%nrcav(2)
            endif

            gy = real(igy-1,wp)*dycav - yhalf

            yy = (gy - iy)/cfgs(ib)%vec%latvector(2,2)
            yy = (yy - nint(yy))*cfgs(ib)%vec%latvector(2,2)

            yy2 = yy*yy
            if( yy2 - rc2 > DZERO ) cycle
            !if( yy2 > rc2 ) cycle

            do jgz = izmin, izmax

                igz = jgz
                if( igz < 1 ) then
                    igz = igz+cfgs(ib)%nrcav(3)
                elseif( igz > cfgs(ib)%nrcav(3) ) then
                    igz = igz-cfgs(ib)%nrcav(3)
                endif

                gz = real(igz-1,wp)*dzcav - zhalf

                zz = (gz - iz)/cfgs(ib)%vec%latvector(3,3)
                zz = (zz - nint(zz))*cfgs(ib)%vec%latvector(3,3)

                zz2 = zz*zz
                if( zz2 - rc2 > DZERO ) cycle
                !!if( zz2 > rc2 ) cycle

                !if( xx2 + yy2 + zz2 - rc2 > DZERO ) cycle

                ign = igz + ((igy-1) + (igx-1)*cfgs(ib)%nrcav(2))*cfgs(ib)%nrcav(3)
!#ifdef DEBUG
!debug<
!                if( abs( cfgs(ib)%cavgrid(1,ign) - gx ) > DZERO .or. &
!                    abs( cfgs(ib)%cavgrid(2,ign) - gy ) > DZERO .or. &
!                    abs( cfgs(ib)%cavgrid(3,ign) - gz ) > DZERO) &
!                    write(*,'(1x,2(a,3f15.7),a,i6,a)')&
!                    "reset_pcav_insert(3):: coordinates for new atom "//&
!                    "given and taken in 1D differ: (",gx,gy,gz,") =?= (",&
!                    cfgs(ib)%cavgrid(1,ign),cfgs(ib)%cavgrid(2,ign),cfgs(ib)%cavgrid(3,ign),&
!                    ") for 1D index ",ign," !!!"
!debug>
!#endif

!#ifndef DEBUG
!nondebug<
                if( cfgs(ib)%cavity(ign) > 0 ) then
!nondebug>
!#endif
                    if( xx2 + yy2 + zz2 - rc2 > DZERO ) cycle

                    !rsq = xx2 + yy2 + zz2 

                    !if( rc2-rsq < DZERO ) cycle
                    !!if( rc2-rsq > DZERO ) then
                    !!!if( rsq < rc2 ) then
!#ifdef DEBUG
!debug<
                        nrem = nrem + 1

                        if( cfgs(ib)%cavity(ign) > 0 ) then

                            if( ign /= ig ) then
                                write(uout,*)"reset_pcav_insert(",nrem,&
                                "): NON selected busy cavity freed at atom insertion fail ",&
                                ign," =/=",ig,"(",cfgs(ib)%cavity(ign),&
                                "), free_tot = ",ncavity_free(ib)
                            else
                                write(uout,*)"reset_pcav_insert(",nrem,&
                                "): the selected busy cavity freed at atom insertion fail ",&
                                ign," =?=",ig,"(",cfgs(ib)%cavity(ign),&
                                "), free_tot = ",ncavity_free(ib)
                            endif
!debug>
!#endif
                            cfgs(ib)%cavity(ign) = cfgs(ib)%cavity(ign) - 1

                            if( cfgs(ib)%cavity(ign) == 0 ) then

                                ncavity_free(ib) = ncavity_free(ib) + 1

                                cfgs(ib)%cavids(ncavity_free(ib)) = ign
!#ifdef DEBUG
!debug<
                                nfreed = nfreed + 1

!                            else
!                                call cry(uout,'', &
!                                     "WARNING: reset_pcav_insert() - "//&
!                                     "cavity count > 1 after new atom removal !!!",0)
                                     !"ERROR: reset_pcav_insert() - "//&
                                     !"cavity count > 1 after new atom removal !!!",999)
!debug>
!#endif
                            endif
!#ifdef DEBUG
!debug<
                        else

                            if( ign /= ig ) then
                                write(uout,*)"reset_pcav_insert(",nrem,&
                                "): NON-selected 'busy' cavity is already free at atom insertion fail ",&
                                ign," =/=",ig,"(",cfgs(ib)%cavity(ign),"), free_tot = ",ncavity_free(ib)
                            else 
                                write(uout,*)"reset_pcav_insert(",nrem,&
                                "): the selected 'busy' cavity is already free at atom insertion fail ",&
                                ign," =?=",ig,"(",cfgs(ib)%cavity(ign),"), free_tot = ",ncavity_free(ib)
                            endif

                            found = .false.
                            do k = 1,ncavity_free(ib)

                                 if( cfgs(ib)%cavids(k) == ign ) then
                                     found = .true.
                                     exit
                                 end if

                            enddo

                            if( found ) then
                                call cry(uout,'', &
                                     "WARNING: reset_pcav_insert() - cavity count < 1 "// &
                                     "for free grid point at atom removal (failed insertion) "// &
                                     trim(int_2_word(natt))//" in cavity ("//trim(int_2_word(nrem))//","// &
                                     trim(int_2_word(k))//" -> "//trim(int_2_word(ign))//") !!!",0)
                            else
                                call cry(uout,'', &
                                     "ERROR: reset_pcav_insert() - cavity count < 1 "// &
                                     "for busy grid point at atom removal (failed insertion) "// &
                                     trim(int_2_word(natt))//" in cavity ("//trim(int_2_word(nrem))//","// &
                                     trim(int_2_word(ign))// &
                                     ") (not found amongst free cavities) !!!",999)
                            end if

                        end if
!debug>
!#endif
                    !end if !( rsq < rc2 )
!#ifndef DEBUG
!nondebug<
                endif !( cfgs(ib)%cavity(ign) > 0 )
!nondebug>
!#endif
            enddo

        enddo

    enddo

    if( nfreed < 1 ) &
        call cry(uout,'', &
             "WARNING: reset_pcav_insert() - "//&
             "freed cavity count < 1 after new atom removal !!!",0)
             !"ERROR: reset_pcav_insert() - "//&
             !"cavity count > 1 after new atom removal !!!",999)

    if( nfreed /= ntaken ) &
        call cry(uout,'', &
             "WARNING: reset_pcav_insert() - "//&
             "number of freed cavitities does not equal that taken "//&
             "(after inserted atom removal) !!!",0)

!debug>
#endif

end subroutine reset_pcav_insert

!AB: new, fully reworked and tested variant for updating the cavity list(s) upon atom removal
subroutine calc_pcav_remove(ib, crad, pos, prob, is_first)

    !use kinds_f90
    use constants_module, only : HALF, DZERO, uout
    use parse_module, only : int_2_word

    !implicit none

    integer, intent(in)            :: ib
    !integer, intent(inout)         :: ig
    real(kind = wp), intent(in)    :: crad
    real(kind = wp), intent(inout) :: pos(:), prob
    logical, intent(in)            :: is_first

    integer :: igx, igy, igz, ign, ipm
    integer :: ixmin,ixmax,iymin,iymax,izmin,izmax
    integer :: jgx, jgy, jgz, nfree, k, ncmax, ncmin

    real(kind = wp) :: rc2, dxcav, dycav, dzcav, xhalf, yhalf, zhalf
    real(kind = wp) :: ix, iy, iz, gx, gy, gz, xx, yy, zz, xx2, yy2, zz2 !, rsq

#ifdef DEBUG
!debug<
    logical :: found

    integer, save :: natt = 0
    integer :: nrem

    nrem  = 0
!debug>
#endif

    ncmax =-1
    ncmin = ncavity_grid(ib)

    prob  = 1.0_wp
    if( ncavity_free(ib) < 1 ) return

#ifdef DEBUG
!debug<
    natt = natt + 1
!debug>
#endif

    rc2 = crad * crad

    dxcav = cfgs(ib)%drcav(1)
    dycav = cfgs(ib)%drcav(2)
    dzcav = cfgs(ib)%drcav(3)

    xhalf = cfgs(ib)%vec%latvector(1,1)*HALF
    yhalf = cfgs(ib)%vec%latvector(2,2)*HALF
    zhalf = cfgs(ib)%vec%latvector(3,3)*HALF

    ix = pos(1)
    iy = pos(2)
    iz = pos(3)

!AB: need to apply PBC to map onto the cavity mesh
    !if( .not.cfgs(ib)%vec%is_orthogonal ) &
        call pbc_cart2frac_pos_calc(ib, ix, iy, iz)

!AB: optimisation - update cavities in box [ix+/-cavx, iy+/-cavy, iz+/-cavz]

    igx = int((ix+xhalf)/dxcav)+1
    igy = int((iy+yhalf)/dycav)+1
    igz = int((iz+zhalf)/dzcav)+1

    xhalf = xhalf-dxcav*HALF
    yhalf = yhalf-dycav*HALF
    zhalf = zhalf-dzcav*HALF

    ign = igz + ((igy-1) + (igx-1)*cfgs(ib)%nrcav(2))*cfgs(ib)%nrcav(3)

    ipm   = int(crad/dxcav)+1
    ixmin = igx-ipm
    ixmax = igx+ipm

    ipm   = int(crad/dycav)+1
    iymin = igy-ipm
    iymax = igy+ipm

    ipm   = int(crad/dzcav)+1
    izmin = igz-ipm
    izmax = igz+ipm

    !if( is_first ) cavatt(1) = 1
    if( is_first ) then
        cavatt(1) = 1
        nfreed = 0

        !if( ig < 1 ) then
        !    nfreed = 1
        !    cavatt(1) = 2
        !    cavatt(2) = -ign
        !endif
    endif

    do jgx = ixmin, ixmax

        igx = jgx
        if( igx < 1 ) then
            igx = igx+cfgs(ib)%nrcav(1)
        elseif( igx > cfgs(ib)%nrcav(1) ) then
            igx = igx-cfgs(ib)%nrcav(1)
        endif

        gx = real(igx-1,wp)*dxcav - xhalf

        xx = (gx - ix)/cfgs(ib)%vec%latvector(1,1)
        xx = (xx - nint(xx))*cfgs(ib)%vec%latvector(1,1)
                
        xx2 = xx*xx
        if( xx2 - rc2 > DZERO ) cycle
        !if( xx2 > rc2 ) cycle

        do jgy = iymin, iymax

            igy = jgy
            if( igy < 1 ) then
                igy = igy+cfgs(ib)%nrcav(2)
            elseif( igy > cfgs(ib)%nrcav(2) ) then
                igy = igy-cfgs(ib)%nrcav(2)
            endif

            gy = real(igy-1,wp)*dycav - yhalf

            yy = (gy - iy)/cfgs(ib)%vec%latvector(2,2)
            yy = (yy - nint(yy))*cfgs(ib)%vec%latvector(2,2)

            yy2 = yy*yy
            if( yy2 - rc2 > DZERO ) cycle
            !if( yy2 > rc2 ) cycle

            do jgz = izmin, izmax

                igz = jgz
                if( igz < 1 ) then
                    igz = igz+cfgs(ib)%nrcav(3)
                elseif( igz > cfgs(ib)%nrcav(3) ) then
                    igz = igz-cfgs(ib)%nrcav(3)
                endif

                gz = real(igz-1,wp)*dzcav - zhalf

                zz = (gz - iz)/cfgs(ib)%vec%latvector(3,3)
                zz = (zz - nint(zz))*cfgs(ib)%vec%latvector(3,3)

                zz2 = zz*zz
                if( zz2 - rc2 > DZERO ) cycle
                !!if( zz2 > rc2 ) cycle

                if( xx2 + yy2 + zz2 - rc2 > DZERO ) cycle

                ign = igz + ((igy-1) + (igx-1)*cfgs(ib)%nrcav(2))*cfgs(ib)%nrcav(3)

                ncmax = max(ncmax,cfgs(ib)%cavity(ign))
                ncmin = min(ncmin,cfgs(ib)%cavity(ign))
#ifdef DEBUG
!debug<
                if( abs( cfgs(ib)%cavgrid(1,ign) - gx ) > DZERO .or. &
                    abs( cfgs(ib)%cavgrid(2,ign) - gy ) > DZERO .or. &
                    abs( cfgs(ib)%cavgrid(3,ign) - gz ) > DZERO) &
                    write(*,'(1x,2(a,3f15.7),a,i6,a)')&
                    "calc_pcav_remove(..):: coordinates for new atom "//&
                    "given and taken in 1D differ: (",gx,gy,gz,") =?= (",&
                    cfgs(ib)%cavgrid(1,ign),cfgs(ib)%cavgrid(2,ign),cfgs(ib)%cavgrid(3,ign),&
                    ") for 1D index ",ign," !!!"
!debug>
#endif
                !rsq = xx2 + yy2 + zz2 
#ifndef DEBUG
!nondebug<
                !if( rc2-rsq < DZERO ) cycle

                if( cfgs(ib)%cavity(ign) > 0 ) then
!nondebug>
#endif
!                    if( xx2 + yy2 + zz2 - rc2 > DZERO ) cycle
#ifdef DEBUG
!debug<
                    !if( rc2-rsq < DZERO ) cycle
                    !!if( rc2-rsq > DZERO ) then
                    !!!if( rsq < rc2 ) then
                        nrem = nrem+1

                        write(uout,*)
                        if( cfgs(ib)%cavity(ign) > 0 ) then

                            write(uout,*)"calc_pcav_remove(",nrem,&
                            "): found a busy cavity to free at atom removal ",&
                            ign," (",cfgs(ib)%cavity(ign),"), free_tot = ",ncavity_free(ib)
!debug>
#endif
                            cfgs(ib)%cavity(ign) = cfgs(ib)%cavity(ign) - 1

                            cavatt(1) = cavatt(1)+1
                            cavatt(cavatt(1)) = ign

                            if( cfgs(ib)%cavity(ign) == 0 ) then

                                ncavity_free(ib) = ncavity_free(ib) + 1

                                cfgs(ib)%cavids(ncavity_free(ib)) = ign

                                nfreed = nfreed+1

                            endif
#ifdef DEBUG
!debug<
                        else

                          found = .false.
                          do k = 1,ncavity_free(ib)
                             if( cfgs(ib)%cavids(k) == ign ) then
                                 found = .true.
                                 exit
                             end if
                          enddo

                          if( found ) then

                              write(uout,*)"calc_pcav_remove(",nrem,&
                              "): found an already free cavity for atom removal ",&
                              k," -> ",ign," (",cfgs(ib)%cavity(ign),"), free_tot = ",ncavity_free(ib)

                                call cry(uout,'', &
                                     "ERROR: calc_pcav_remove() - "//&
                                     "cavity count < 1 before atom removal !!!",999)
                          else

                              write(uout,*)"calc_pcav_remove(",nrem,&
                              "): found an already free cavity for atom removal ",&
                              k," -? ",ign," (",cfgs(ib)%cavity(ign),"), free_tot = ",ncavity_free(ib)

                                call cry(uout,'', &
                                     "ERROR: calc_pcav_remove() - "//&
                                     "cavity count < 1 before atom removal "// &
                                     "but not found amongst free cavities !!!",999)
                          end if

                        end if

                   !end if !( rsq < rc2 )
!debug>
#endif

#ifndef DEBUG
!nondebug<
                endif !( cfgs(ib)%cavity(ign) > 0 )
!nondebug>
#endif
            enddo

        enddo

    enddo

#ifdef DEBUG
!debug<
    !if( .not.is_first .and. nfreed < 1 ) &
    if( nfreed < 1 ) &
        call cry(uout,'', &
             "WARNING: calc_pcav_remove() - "//&
             "no cavity freed upon removing an existing atom "//&
             " (cavities affected: "//trim(int_2_word(cavatt(1)))//&
             " min/max = "//trim(int_2_word(ncmin))//" / "//trim(int_2_word(ncmax))//") !!!",0)
!             " (cavities affected: "//trim(int_2_word(cavatt(1)))//") !!!",0)
             !"no cavity freed upon removing an existing atom (off the grid?)...",0)
!debug>
#endif

!AB: for detailed balance the probability is to be estimated after atom removal
    if( ncavity_free(ib) == 0 ) then
        prob = 1.0_wp
    else
        prob = real(ncavity_free(ib),wp) / real(ncavity_grid(ib),wp)
    endif

end subroutine calc_pcav_remove

!AB: new, fully reworked and tested variant for updating the cavity list(s) upon molecule removal
subroutine calc_pcav_remove_mol(ib, im, crad, prob, by_com)

    !use kinds_f90
    use constants_module, only : FZERO, uout

    !implicit none

    integer, intent(in)            :: ib, im
    real(kind = wp), intent(in)    :: crad
    real(kind = wp), intent(inout) :: prob
    logical, intent(in)            :: by_com

    integer :: ig, ia
    logical :: is_first ! = .true.

    real(kind = wp) :: pcav

    prob = 1.0_wp

    !ig = ig0
    !nc = cfgs(ib)%cavity(ig)
    !ig =-ig

    is_first = .true.

    if( by_com ) then

        !call calc_pcav_remove(ib, crad, &
        !     cfgs(ib)%mols(im)%rcom, prob, is_first)

        !is_first = .false.

        do ia = 1, cfgs(ib)%mols(im)%natom

!AB: skip fictitious atoms if any (set mass to zero)
            if( cfgs(ib)%mols(im)%atms(ia)%mass < FZERO ) cycle

            pcav = 1.0_wp

            call calc_pcav_remove(ib, crad, &
                 cfgs(ib)%mols(im)%atms(ia)%rpos, pcav, is_first)

            prob = pcav !prob*pcav

            is_first = .false.

        enddo

    else

        do ia = 1, cfgs(ib)%mols(im)%natom

!AB: skip fictitious atoms if any (set mass to zero)
            if( cfgs(ib)%mols(im)%atms(ia)%mass < FZERO ) cycle

            pcav = 1.0_wp

            call calc_pcav_remove(ib, crad, &
                 cfgs(ib)%mols(im)%atms(ia)%rpos, pcav, is_first)

            prob = prob*pcav

            is_first = .false.

        enddo

    endif

end subroutine calc_pcav_remove_mol

!AB: new, fully reworked and tested variant for updating the cavity list(s) upon atom re-insertion (failed removal)
subroutine reset_pcav_remove(ib, crad, pos)

    !use kinds_f90
    use constants_module, only : HALF, DZERO, uout
    use parse_module, only : int_2_word

    !implicit none

    integer, intent(in) :: ib !, im, ia

    real(kind = wp), intent(in)    :: crad
    real(kind = wp), intent(inout) :: pos(:)

    integer :: k, ign, nfree, ncmax, ncmin

#ifdef DEBUG
!debug<
    real(kind = wp) :: rc2, dxcav, dycav, dzcav, xhalf, yhalf, zhalf
    real(kind = wp) :: ix, iy, iz, gx, gy, gz, xx, yy, zz, xx2, yy2, zz2 !, rsq

    integer :: igx, igy, igz, ipm
    integer :: ixmin,ixmax,iymin,iymax,izmin,izmax
    integer :: jgx, jgy, jgz

    logical :: found

    integer, save :: natt = 0
    integer :: nins

    natt  = natt + 1
    nins  = 0

    if( ncavity_free(ib) < 1 ) & !return
        call cry(uout,'', &
             "ERROR: reset_pcav_remove() - "//&
             "free cavity count < 1 while re-inserting atom back !!!",999)
!debug>
#endif

    ntaken = 0

    ncmax =-1
    ncmin = ncavity_grid(ib)

#ifndef DEBUG
!nondebug<

    do k = 2, cavatt(1)
        ign = cavatt(k)

        ncmax = max(ncmax,cfgs(ib)%cavity(ign))
        ncmin = min(ncmin,cfgs(ib)%cavity(ign))

        if( cfgs(ib)%cavity(ign) == 0 ) then

            cfgs(ib)%cavids(ncavity_free(ib)) = 0

            ncavity_free(ib) = ncavity_free(ib) - 1

            ntaken = ntaken+1

        endif

        cfgs(ib)%cavity(ign) = cfgs(ib)%cavity(ign) + 1

        cavatt(k) = 0

    enddo

!    if( ntaken < 1 ) &
!        call cry(uout,'', &
!             "WARNING: reset_pcav_remove() - "//&
!             "number of re-occupied cavitities < 1 after removed atom insertion "//&
!             " (cavities affected: "//trim(int_2_word(cavatt(1)))//&
!             " min/max = "//trim(int_2_word(ncmin))//" / "//trim(int_2_word(ncmax))//") !!!",0)

    if( ntaken /= nfreed ) &
        call cry(uout,'', &
             "WARNING: reset_pcav_remove() - "//&
             "number of taken cavitities "//trim(int_2_word(ntaken))//" =?= "//trim(int_2_word(nfreed))//&
             " freed (after removed atom insertion) !!!",0)
             !"number of taken cavitities is not equal to that freed (after removed atom insertion)"//&
             !" !!!",0)

    return

!nondebug>
#endif

#ifdef DEBUG
!debug<

    rc2 = crad * crad

    dxcav = cfgs(ib)%drcav(1)
    dycav = cfgs(ib)%drcav(2)
    dzcav = cfgs(ib)%drcav(3)

    xhalf = cfgs(ib)%vec%latvector(1,1)*HALF
    yhalf = cfgs(ib)%vec%latvector(2,2)*HALF
    zhalf = cfgs(ib)%vec%latvector(3,3)*HALF

    ix = pos(1)
    iy = pos(2)
    iz = pos(3)

!AB: need to apply PBC to map onto the cavity mesh
    !if( .not.cfgs(ib)%vec%is_orthogonal ) &
        call pbc_cart2frac_pos_calc(ib, ix, iy, iz)

!AB: optimisation - update cavities in box [ix+/-cavx, iy+/-cavy, iz+/-cavz]

    igx = int((ix+xhalf)/dxcav)+1
    igy = int((iy+yhalf)/dycav)+1
    igz = int((iz+zhalf)/dzcav)+1

    xhalf = xhalf-dxcav*HALF
    yhalf = yhalf-dycav*HALF
    zhalf = zhalf-dzcav*HALF

    ign = igz + ((igy-1) + (igx-1)*cfgs(ib)%nrcav(2))*cfgs(ib)%nrcav(3)

    ipm   = int(crad/dxcav)+1
    ixmin = igx-ipm
    ixmax = igx+ipm

    ipm   = int(crad/dycav)+1
    iymin = igy-ipm
    iymax = igy+ipm

    ipm   = int(crad/dzcav)+1
    izmin = igz-ipm
    izmax = igz+ipm

    do jgx = ixmin, ixmax

        igx = jgx
        if( igx < 1 ) then
            igx = igx+cfgs(ib)%nrcav(1)
        elseif( igx > cfgs(ib)%nrcav(1) ) then
            igx = igx-cfgs(ib)%nrcav(1)
        endif

        gx = real(igx-1,wp)*dxcav - xhalf

        xx = (gx - ix)/cfgs(ib)%vec%latvector(1,1)
        xx = (xx - nint(xx))*cfgs(ib)%vec%latvector(1,1)
                
        xx2 = xx*xx
        if( xx2 - rc2 > DZERO ) cycle
        !if( xx2 > rc2 ) cycle

        do jgy = iymin, iymax

            igy = jgy
            if( igy < 1 ) then
                igy = igy+cfgs(ib)%nrcav(2)
            elseif( igy > cfgs(ib)%nrcav(2) ) then
                igy = igy-cfgs(ib)%nrcav(2)
            endif

            gy = real(igy-1,wp)*dycav - yhalf

            yy = (gy - iy)/cfgs(ib)%vec%latvector(2,2)
            yy = (yy - nint(yy))*cfgs(ib)%vec%latvector(2,2)

            yy2 = yy*yy
            if( yy2 - rc2 > DZERO ) cycle
            !if( yy2 > rc2 ) cycle

            do jgz = izmin, izmax

                igz = jgz
                if( igz < 1 ) then
                    igz = igz+cfgs(ib)%nrcav(3)
                elseif( igz > cfgs(ib)%nrcav(3) ) then
                    igz = igz-cfgs(ib)%nrcav(3)
                endif

                gz = real(igz-1,wp)*dzcav - zhalf

                zz = (gz - iz)/cfgs(ib)%vec%latvector(3,3)
                zz = (zz - nint(zz))*cfgs(ib)%vec%latvector(3,3)

                zz2 = zz*zz
                if( zz2 - rc2 > DZERO ) cycle
                !!if( zz2 > rc2 ) cycle

                if( xx2 + yy2 + zz2 - rc2 > DZERO ) cycle

                ign = igz + ((igy-1) + (igx-1)*cfgs(ib)%nrcav(2))*cfgs(ib)%nrcav(3)
!#ifdef DEBUG
!debug<
                if( abs( cfgs(ib)%cavgrid(1,ign) - gx ) > DZERO .or. &
                    abs( cfgs(ib)%cavgrid(2,ign) - gy ) > DZERO .or. &
                    abs( cfgs(ib)%cavgrid(3,ign) - gz ) > DZERO) &
                    write(*,'(1x,2(a,3f15.7),a,i6,a)')&
                    "reset_pcav_remove(..):: coordinates for new atom given and taken in 1D differ: (",gx,gy,gz,&
                    ") =?= (",cfgs(ib)%cavgrid(1,ign),cfgs(ib)%cavgrid(2,ign),cfgs(ib)%cavgrid(3,ign),&
                    ") for 1D index ",ign," !!!"

                !rsq = xx2 + yy2 + zz2 

                !if( rc2-rsq < DZERO ) cycle
                !!if( rc2-rsq > DZERO ) then
                !!!if( rsq < rc2 ) then
!#ifdef DEBUG
!debug<
                    nins = nins+1
                    write(uout,*)
                    if( ncavity_free(ib) > 0 ) then
!debug>
!#endif
                        found = .false.
                        do k = 1,ncavity_free(ib)

                           if( cfgs(ib)%cavids(k) == ign ) then
                               found = .true.
                               exit
                           end if

                        enddo

                        if( found ) then
!#ifdef DEBUG
!debug<
                            write(uout,*)"reset_pcav_remove(",nins,&
                            "): found a free cavity to reoccupy at atom removal fail ",&
                            k," -> ",ign," (",cfgs(ib)%cavity(ign),"), free_tot = ",ncavity_free(ib)

                            nfree = nfree+1
!debug>
!#endif
                            cfgs(ib)%cavids(k) = cfgs(ib)%cavids(ncavity_free(ib))
                            cfgs(ib)%cavids(ncavity_free(ib)) = 0
                            ncavity_free(ib) = ncavity_free(ib) - 1
!#ifdef DEBUG
!debug<
                        else

                            write(uout,*)"reset_pcav_remove(",nins,&
                            "): found a BUSY cavity to reoccupy at atom removal fail ",&
                            k," -? ",ign," (",cfgs(ib)%cavity(ign),"), free_tot = ",ncavity_free(ib)
!debug>
!#endif
                        end if

                        cfgs(ib)%cavity(ign) = cfgs(ib)%cavity(ign) + 1
!#ifdef DEBUG
!debug<
                    else

                        write(uout,*)"reset_pcav_remove(",nins,&
                        "): did NOT find a free cavity to reoccupy at atom removal fail ",&
                        ign," (",cfgs(ib)%cavity(ign),"), free_tot = ",ncavity_free(ib)

                        call cry(uout,'', &
                             "ERROR: reset_pcav_remove() - "//&
                             "free cavity count < 1 while re-inserting atom back !!!",999)
                    endif
!debug>
!#endif
                !end if !( rsq < rc2 )

            enddo

        enddo

    enddo

!#ifdef DEBUG
!debug<
    if( nfree < 1 ) &
        call cry(uout,'', &
             "ERROR: reset_pcav_remove() - "//&
             "no free cavity found on the list for re-inserting a removed atom !!!",999)
!debug>
#endif

end subroutine reset_pcav_remove

!AB: old (based on the original code) throughout look-up routines for updating the cavity list(s)
!AB: these routines are very inefficient(!) but made correct...

subroutine calc_pcav_insert0(ib, im, ia, crad, ig)

    !use kinds_f90
    use constants_module, only : HALF, DHALF, DZERO, uout
    use parse_module, only : int_2_word

    !implicit none

    integer, intent(in) :: ib, im, ia, ig

    real(kind = wp), intent(in) :: crad

    integer :: i, j, k, ic

    real(kind = wp) :: ix, iy, iz, jx, jy, jz, rad
    real(kind = wp) :: xx, yy, zz, rx, ry, rz, rsq

    logical :: found !,found2

#ifdef DEBUG
!debug<
    integer, save :: natt = 0
    integer :: nins, nfree

    natt  = natt + 1
    nins  = 0
    nfree = 0

!    if (ncavity_free(ib) == 0) then
!        prob = -1.0_wp
!        return
!    else
!        prob = real(ncavity_free(ib),wp) / real(ncavity_grid(ib),wp)
!    endif
!debug>
#endif

#ifndef DEBUG
!nondebug<
    if( ncavity_free(ib) < 1 ) return
!nondebug>
#endif

    rad = crad * crad

    ! search over atom
    ix = cfgs(ib)%mols(im)%atms(ia)%rpos(1)
    iy = cfgs(ib)%mols(im)%atms(ia)%rpos(2)
    iz = cfgs(ib)%mols(im)%atms(ia)%rpos(3)

    do j = 1, ncavity_grid(ib)

        jx = cfgs(ib)%cavgrid(1,j)
        jy = cfgs(ib)%cavgrid(2,j)
        jz = cfgs(ib)%cavgrid(3,j)

        xx = jx - ix
        yy = jy - iy
        zz = jz - iz

        call pbc_cart_pos_calc(ib, xx, yy, zz)

        rsq = xx * xx + yy * yy + zz * zz 

        if( rad-rsq > DZERO ) then
        !if( rsq < rad ) then

#ifdef DEBUG
!debug<
            nins = nins + 1

            !write(uout,*)
            if( ncavity_free(ib) > 0 ) then
!debug>
#endif

                found  = .false.
                !found2 = .false.
                do k = 1,ncavity_free(ib)

                   if( cfgs(ib)%cavids(k) == j ) then
                       !if( found ) then
                       !    found2 = .true.
                       !    ic = k
                       !else

                       found = .true.
                       ic = k

                       !endif

                       exit
                   end if

                enddo

                if( found ) then
#ifdef DEBUG
!debug<
                    if( j /= ig ) then
                        !write(uout,*)
                        write(uout,*)"calc_pcav_insert(",nins,"): NON-selected free cavity taken at atom insertion ",&
                        j," =/=",ig," / ",ic," -> ",cfgs(ib)%cavids(ic)," (",cfgs(ib)%cavity(j), &
                        "), free_tot = ",ncavity_free(ib)
                        !j," =/=",ig," / ",ic," -> ",cfgs(ib)%cavids(ic)," (",cfgs(ib)%cavids(ncavity_free(ib)), &
                    else
                        !write(uout,*)
                        write(uout,*)"calc_pcav_insert(",nins,"): the selected free cavity taken at atom insertion ",&
                        j," =?=",ig," / ",ic," -> ",cfgs(ib)%cavids(ic)," (",cfgs(ib)%cavity(j), &
                        "), free_tot = ",ncavity_free(ib)
                        !j," =?=",ig," / ",ic," -> ",cfgs(ib)%cavids(ic)," (",cfgs(ib)%cavids(ncavity_free(ib)), &
                    endif

                    !if( found2 ) &
                    !    call cry(uout,'', &
                    !            "ERROR: calc_pcav_insert() - free cavity "//trim(int_2_word(j))//&
                    !            " found more than once in cavids(..) !!!",999)

                    if( cfgs(ib)%cavity(j) > 0 ) then
                        call cry(uout,'', &
                                "ERROR: calc_pcav_insert() - cavity count > 0 at new atom insertion "// &
                                "(yet found amongst free cavities) !!!",999)
                    endif

                    nfree = nfree+1
!debug>
#endif

                    cfgs(ib)%cavids(ic) = cfgs(ib)%cavids(ncavity_free(ib))
                    cfgs(ib)%cavids(ncavity_free(ib)) = 0
                    ncavity_free(ib) = ncavity_free(ib) - 1
#ifdef DEBUG
!debug<
                    cfgs(ib)%cavity(j) = cfgs(ib)%cavity(j) + 1
!debug>
#endif

#ifndef DEBUG
!nondebug<
                end if

                cfgs(ib)%cavity(j) = cfgs(ib)%cavity(j) + 1
!nondebug>
#endif

#ifdef DEBUG
!debug<
                else if( cfgs(ib)%cavity(j) > 0 ) then

                    !write(uout,*)
                    write(uout,*)"calc_pcav_insert(",nins,"): NON-selected & NON-free cavity found at atom insertion ",&
                    j," (",cfgs(ib)%cavity(j),"), free_tot = ",ncavity_free(ib)

                    cfgs(ib)%cavity(j) = cfgs(ib)%cavity(j) + 1

                    call cry(uout,'', &
                            "WARNING: calc_pcav_insert() - cavity count > 0 at atom insertion attempt "// &
                            trim(int_2_word(natt))//" in cavity ("//trim(int_2_word(nins))//","//trim(int_2_word(j))// &
                            ") (not found amongst free cavities) !!!",0)
                            !"ERROR: calc_pcav_insert() - cavity count > 0 at atom insertion attempt "// &

                else
                    write(uout,*)"calc_pcav_insert(",nins,"): NON-selected free cavity found at atom insertion ",&
                    j," (",cfgs(ib)%cavity(j),"), free_tot = ",ncavity_free(ib)

                    call cry(uout,'', &
                            "WARNING: calc_pcav_insert() - found an EXTRA free cavity at atom insertion attempt "// &
                            trim(int_2_word(natt))//" in cavity ("//trim(int_2_word(nins))//","//trim(int_2_word(j))// &
                            ") (not found amongst free cavities) !!!",0)
!!debug>
!#endif
                end if
!#ifdef DEBUG
!!debug<
            else

                if( j == ig ) then
                    write(uout,*)"calc_pcav_insert(",nins,"): the selected free cavity for atom insertion is busy ",ig
                endif

                !call create_cavitylist(ib, crad)

                call cry(uout,'', &
                        "ERROR: calc_pcav_insert() - free cavity count < 1 at atom insertion attempt "// &
                        trim(int_2_word(natt))//" in cavity "//trim(int_2_word(nins))//","//trim(int_2_word(j))//" !!!",999)
            endif
!debug>
#endif

        end if

    enddo

#ifdef DEBUG
!debug<
    if( ncavity_free(ib) < 1 ) call cry(uout,'', &
            "WARNING: calc_pcav_insert() - free cavity count < 1 upon inserting new atom !!!",0)

    if( nfree < 1 ) call cry(uout,'', &
            "ERROR: calc_pcav_insert() - no existing free cavity occupied upon inserting new atom !!!",999)
!debug>
#endif

end subroutine

subroutine reset_pcav_insert0(ib, im, i, crad, ig)

    !use kinds_f90
    use constants_module, only : DZERO, uout
    use parse_module, only : int_2_word

    !implicit none

    integer, intent(in) :: ib, im, i, ig

    real(kind = wp), intent(in) :: crad

    integer :: j, k

    real(kind = wp) :: ix, iy, iz, jx, jy, jz, rad
    real(kind = wp) :: xx, yy, zz, rx, ry, rz, rsq

#ifdef DEBUG
!debug<
    logical :: found

    integer, save :: natt = 0
    integer :: nrem, nfree

    natt = natt + 1
    nrem = 0

    !nfree = ncavity_free(ib)
!debug>
#endif

    rad = crad * crad

    ! search over atom
    ix = cfgs(ib)%mols(im)%atms(i)%rpos(1)
    iy = cfgs(ib)%mols(im)%atms(i)%rpos(2)
    iz = cfgs(ib)%mols(im)%atms(i)%rpos(3)

    do j = 1, ncavity_grid(ib)

#ifndef DEBUG
!nondebug<
        if( cfgs(ib)%cavity(j) > 0 ) then
!nondebug>
#endif

            jx = cfgs(ib)%cavgrid(1,j)
            jy = cfgs(ib)%cavgrid(2,j)
            jz = cfgs(ib)%cavgrid(3,j)

            xx = jx - ix
            yy = jy - iy
            zz = jz - iz

            call pbc_cart_pos_calc(ib, xx, yy, zz)

            rsq = xx * xx + yy * yy + zz * zz

            if( rad-rsq > DZERO ) then
            !if( rsq < rad ) then

#ifdef DEBUG
!debug<
                nrem = nrem + 1

                if( cfgs(ib)%cavity(j) > 0 ) then

                    if( j /= ig ) then
                        write(uout,*)"reset_pcav_insert(",nrem,"): NON selected busy cavity freed at atom insertion fail ",&
                            j," =/=",ig,"(",cfgs(ib)%cavity(j),"), free_tot = ",ncavity_free(ib)
                    else
                        write(uout,*)"reset_pcav_insert(",nrem,"): the selected busy cavity freed at atom insertion fail ",&
                        j," =?=",ig,"(",cfgs(ib)%cavity(j),"), free_tot = ",ncavity_free(ib)
                    endif
!debug>
#endif
                    cfgs(ib)%cavity(j) = cfgs(ib)%cavity(j) - 1

                    if( cfgs(ib)%cavity(j) == 0 ) then

                        ncavity_free(ib) = ncavity_free(ib) + 1

                        cfgs(ib)%cavids(ncavity_free(ib)) = j

#ifdef DEBUG
!debug<
                    else
                        call cry(uout,'', &
                             "WARNING: reset_pcav_insert() - cavity count > 1 after new atom removal !!!",0)
                             !"ERROR: reset_pcav_insert() - cavity count > 1 after new atom removal !!!",999)

!debug>
#endif
                    endif
#ifdef DEBUG
!debug<
                else

                    if( j /= ig ) then
                        write(uout,*)"reset_pcav_insert(",nrem,&
                        "): NON-selected 'busy' cavity is already free at atom insertion fail ",&
                        j," =/=",ig,"(",cfgs(ib)%cavity(j),"), free_tot = ",ncavity_free(ib)
                    else 
                        write(uout,*)"reset_pcav_insert(",nrem,&
                        "): the selected 'busy' cavity is already free at atom insertion fail ",&
                        j," =?=",ig,"(",cfgs(ib)%cavity(j),"), free_tot = ",ncavity_free(ib)
                    endif

                    found = .false.
                    do k = 1,ncavity_free(ib)

                         if( cfgs(ib)%cavids(k) == j ) then
                             found = .true.
                             exit
                         end if

                    enddo

                    if( found ) then
                        call cry(uout,'', &
                                "WARNING: reset_pcav_insert() - cavity count < 1 "// &
                                "for free grid point at atom removal (failed insertion) "// &
                                trim(int_2_word(natt))//" in cavity ("//trim(int_2_word(nrem))//","// &
                                trim(int_2_word(k))//" -> "//trim(int_2_word(j))//") !!!",0)
                    else
                        call cry(uout,'', &
                                "ERROR: reset_pcav_insert() - cavity count < 1 "// &
                                "for busy grid point at atom removal (failed insertion) "// &
                                trim(int_2_word(natt))//" in cavity ("//trim(int_2_word(nrem))//","// &
                                trim(int_2_word(j))// &
                                ") (not found amongst free cavities) !!!",999)
                    end if

                end if
!debug>
#endif
            end if

#ifndef DEBUG
!nondebug<
        endif !( cfgs(ib)%cavity(j) > 0 )
!nondebug>
#endif

    enddo

end subroutine


subroutine calc_pcav_remove0(ib, im, i, crad, prob)

    !use kinds_f90
    use constants_module, only : DZERO, uout

    !implicit none

    integer, intent(in) :: ib, im, i

    real(kind = wp), intent(in) :: crad

    real(kind = wp), intent(out) :: prob

    integer :: j, k

    real(kind = wp) :: ix, iy, iz, jx, jy, jz, rad
    real(kind = wp) :: xx, yy, zz, rx, ry, rz, rsq

#ifdef DEBUG
!debug<
    logical :: found

    integer, save :: natt = 0
    integer :: nrem, nfree

    natt = natt + 1
    nrem = 0
    nfree= 0
!debug>
#endif

    rad = crad * crad

    ! search over atom
    ix = cfgs(ib)%mols(im)%atms(i)%rpos(1)
    iy = cfgs(ib)%mols(im)%atms(i)%rpos(2)
    iz = cfgs(ib)%mols(im)%atms(i)%rpos(3)

    do j = 1, ncavity_grid(ib)

#ifndef DEBUG
!nondebug<
        if( cfgs(ib)%cavity(j) > 0 ) then
!nondebug>
#endif

            jx = cfgs(ib)%cavgrid(1,j)
            jy = cfgs(ib)%cavgrid(2,j)
            jz = cfgs(ib)%cavgrid(3,j)

            xx = jx - ix
            yy = jy - iy
            zz = jz - iz

            call pbc_cart_pos_calc(ib, xx, yy, zz)

            rsq = xx * xx + yy * yy + zz * zz

            if( rad-rsq > DZERO ) then
            !if( rsq < rad ) then

#ifdef DEBUG
!debug<
                nrem = nrem+1

                write(uout,*)
                if( cfgs(ib)%cavity(j) > 0 ) then

                    write(uout,*)"calc_pcav_remove(",nrem,&
                    "): found a busy cavity to free at atom removal ",&
                    j," (",cfgs(ib)%cavity(j),"), free_tot = ",ncavity_free(ib)
!debug>
#endif
                    cfgs(ib)%cavity(j) = cfgs(ib)%cavity(j) - 1

                    if( cfgs(ib)%cavity(j) == 0 ) then

                        ncavity_free(ib) = ncavity_free(ib) + 1

                        cfgs(ib)%cavids(ncavity_free(ib)) = j

#ifdef DEBUG
!debug<
                        nfree = nfree+1
#endif

                    endif
#ifdef DEBUG
!debug<
                else

                  found = .false.
                  do k = 1,ncavity_free(ib)
                     if( cfgs(ib)%cavids(k) == j ) then
                         found = .true.
                         exit
                     end if
                  enddo

                  if( found ) then

                      write(uout,*)"calc_pcav_remove(",nrem,&
                      "): found an already free cavity for atom removal ",&
                      k," -> ",j," (",cfgs(ib)%cavity(j),"), free_tot = ",ncavity_free(ib)

                        call cry(uout,'', &
                                "ERROR: calc_pcav_remove() - cavity count < 1 before atom removal !!!",999)
                  else

                      write(uout,*)"calc_pcav_remove(",nrem,&
                      "): found an already free cavity for atom removal ",&
                      k," -? ",j," (",cfgs(ib)%cavity(j),"), free_tot = ",ncavity_free(ib)

                        call cry(uout,'', &
                                "ERROR: calc_pcav_remove() - cavity count < 1 before atom removal"// &
                                " but not found amongst free cavities !!!",999)
                  end if

                end if
!debug>
#endif
           end if

#ifndef DEBUG
!nondebug<
        endif !( cfgs(ib)%cavity(j) > 0 )
!nondebug>
#endif

    enddo

#ifdef DEBUG
!debug<
    if( nfree < 1 ) call cry(uout,'', &
                         "ERROR: calc_pcav_remove() - no new free cavity found upon removing an existing atom !!!",999)
!debug>
#endif

    if (ncavity_free(ib) == 0) then

        prob = 1.0_wp

    else

        prob = real(ncavity_free(ib),wp) / real(ncavity_grid(ib),wp)

    endif

end subroutine


subroutine reset_pcav_remove0(ib, im, i, crad)

    !use kinds_f90
    use constants_module, only : DZERO, uout

    !implicit none

    integer, intent(in) :: ib, im, i

    real(kind = wp), intent(in) :: crad

    integer :: j, k

    real(kind = wp) :: ix, iy, iz, jx, jy, jz, rad
    real(kind = wp) :: xx, yy, zz, rx, ry, rz, rsq

    logical :: found

#ifdef DEBUG
!debug<
    integer, save :: natt = 0
    integer :: nins, nfree

    natt = natt + 1
    nins = 0
    nfree= 0
!debug>
#endif

#ifndef DEBUG
!nondebug<
    if( ncavity_free(ib) < 1 ) return
!nondebug>
#endif

    rad = crad * crad

    ! search over atom
    ix = cfgs(ib)%mols(im)%atms(i)%rpos(1)
    iy = cfgs(ib)%mols(im)%atms(i)%rpos(2)
    iz = cfgs(ib)%mols(im)%atms(i)%rpos(3)

    do j = 1, ncavity_grid(ib)

        jx = cfgs(ib)%cavgrid(1,j)
        jy = cfgs(ib)%cavgrid(2,j)
        jz = cfgs(ib)%cavgrid(3,j)

        xx = jx - ix
        yy = jy - iy
        zz = jz - iz

        call pbc_cart_pos_calc(ib, xx, yy, zz)

        rsq = xx * xx + yy * yy + zz * zz 

        if( rad-rsq > DZERO ) then
        !if( rsq < rad ) then

#ifdef DEBUG
!debug<
            nins = nins+1
            write(uout,*)
            if( ncavity_free(ib) > 0 ) then
!debug>
#endif
                found = .false.
                do k = 1,ncavity_free(ib)

                   if( cfgs(ib)%cavids(k) == j ) then
                       found = .true.
                       exit
                   end if

                enddo

                if( found ) then
#ifdef DEBUG
!debug<
                    write(uout,*)"reset_pcav_remove(",nins,&
                    "): found a free cavity to reoccupy at atom removal fail ",&
                    k," -> ",j," (",cfgs(ib)%cavity(j),"), free_tot = ",ncavity_free(ib)

                    nfree = nfree+1
!debug>
#endif

                    cfgs(ib)%cavids(k) = cfgs(ib)%cavids(ncavity_free(ib))
                    cfgs(ib)%cavids(ncavity_free(ib)) = 0
                    ncavity_free(ib) = ncavity_free(ib) - 1

#ifdef DEBUG
!debug<
                else

                    write(uout,*)"reset_pcav_remove(",nins,&
                    "): found a BUSY cavity to reoccupy at atom removal fail ",&
                    k," -? ",j," (",cfgs(ib)%cavity(j),"), free_tot = ",ncavity_free(ib)
!debug>
#endif
                end if

                cfgs(ib)%cavity(j) = cfgs(ib)%cavity(j) + 1
#ifdef DEBUG
!debug<
            else

                write(uout,*)"reset_pcav_remove(",nins,&
                "): did NOT find a free cavity to reoccupy at atom removal fail ",&
                j," (",cfgs(ib)%cavity(j),"), free_tot = ",ncavity_free(ib)

                call cry(uout,'', &
                        "ERROR: reset_pcav_remove() - free cavity count < 1 while re-inserting atom back !!!",999)
            endif
!debug>
#endif

        end if

    enddo

#ifdef DEBUG
!debug<
    if( nfree < 1 ) &
        call cry(uout,'', &
             "ERROR: reset_pcav_remove() - "//&
             "no free cavity found on the list for re-inserting a removed atom !!!",999)
!debug>
#endif

end subroutine

!AB: *** the old/original/untested and now obsolete variants of cavity list updating routines - start ***

subroutine calc_pcav_molremove(ib, im, crad, prob)

    !use kinds_f90
    use constants_module, only : uout

    !implicit none

    integer, intent(in) :: ib, im

    real(kind = wp), intent(in) :: crad

    real(kind = wp), intent(out) :: prob

    integer :: i, j

    real(kind = wp) :: ix, iy, iz, jx, jy, jz, rad
    real(kind = wp) :: xx, yy, zz, rx, ry, rz, rsq

    rad = crad * crad

    ! search over molecule
    do i = 1, cfgs(ib)%mols(im)%natom

        ix = cfgs(ib)%mols(im)%atms(i)%rpos(1)
        iy = cfgs(ib)%mols(im)%atms(i)%rpos(2)
        iz = cfgs(ib)%mols(im)%atms(i)%rpos(3)

        do j = 1, ncavity_grid(ib)

!nondebug<
!            if( cfgs(ib)%cavity(j) > 0 ) then
!nondebug>
                jx = cfgs(ib)%cavgrid(1,j)
                jy = cfgs(ib)%cavgrid(2,j)
                jz = cfgs(ib)%cavgrid(3,j)

                xx = jx - ix
                yy = jy - iy
                zz = jz - iz

                call pbc_cart_pos_calc(ib, xx, yy, zz)

                rsq = xx * xx + yy * yy + zz * zz 
                
                !if( rsq < rad) cfgs(ib)%cavity(j) = cfgs(ib)%cavity(j) - 1

                if( rsq < rad ) then
!debug<
!AB: NOTE that when inserting a molecule only one cavity it taken (by COM)
!AB: yet, here several cavities can be created (those freed by each atom)
!AB: this discrepancy MUST be eliminated! -- wrt either COM or each atom
!AB: in order to make cavity treatment consistent upon insertion & deletion

                    if( cfgs(ib)%cavity(j) > 0 ) then
!debug>
                        cfgs(ib)%cavity(j) = cfgs(ib)%cavity(j) - 1

                        if( cfgs(ib)%cavity(j) == 0 ) then

                            ncavity_free(ib) = ncavity_free(ib) + 1

                            cfgs(ib)%cavids(ncavity_free(ib)) = j

                        endif
!debug<
                    else

                        call cry(uout,'', &
                                "ERROR: calc_pcav_remove() - cavity count < 1 before atom removal !!!",999)

                    end if
!debug>
                end if

                cycle

!AB: the rest of the loop should be obsolete

                if(cfgs(ib)%vec%is_orthogonal) then

                    call pbc_ortho_coor_calc(ib, xx, yy, zz)
                    rsq = xx * xx + yy * yy + zz * zz 

                else
                    rx = cfgs(ib)%vec%invlat(1,1) * xx + cfgs(ib)%vec%invlat(2,1) * yy + cfgs(ib)%vec%invlat(3,1) * zz
                    ry = cfgs(ib)%vec%invlat(1,2) * xx + cfgs(ib)%vec%invlat(2,2) * yy + cfgs(ib)%vec%invlat(3,2) * zz
                    rz = cfgs(ib)%vec%invlat(1,3) * xx + cfgs(ib)%vec%invlat(2,3) * yy + cfgs(ib)%vec%invlat(3,3) * zz

                ! AB: unified way to apply PBC/MIC
                    xx = rx
                    yy = ry
                    zz = rz
                    call pbc_atom_pos_calc(xx, yy, zz)

                ! AB: original code for PBC/MIC
                !xx = rx - anint(rx)
                !yy = ry - anint(ry)
                !zz = rz - anint(rz)

                rx = cfgs(ib)%vec%latvector(1,1) * xx + cfgs(ib)%vec%latvector(2,1) * yy + cfgs(ib)%vec%latvector(3,1) * zz
                ry = cfgs(ib)%vec%latvector(1,2) * xx + cfgs(ib)%vec%latvector(2,2) * yy + cfgs(ib)%vec%latvector(3,2) * zz
                rz = cfgs(ib)%vec%latvector(1,3) * xx + cfgs(ib)%vec%latvector(2,3) * yy + cfgs(ib)%vec%latvector(3,3) * zz

                    rsq = rx * rx + ry * ry + rz * rz

                end if

                if( rsq < rad) then

                    cfgs(ib)%cavity(j) = cfgs(ib)%cavity(j) - 1

                endif

!nondebug<
!            endif !( cfgs(ib)%cavity(j) > 0 )
!nondebug>

        enddo

    enddo


    ! sum up free spaces on grid
!    ncavity_free(ib)   = 0
!    cfgs(ib)%cavids(:) = 0
!    
!    do j = 1, ncavity_grid(ib)
!
!        if( cfgs(ib)%cavity(j) == 0 ) then
!
!            ncavity_free(ib) = ncavity_free(ib) + 1
!
!            cfgs(ib)%cavids(ncavity_free(ib)) = j
!
!        endif
!
!    enddo

    if (ncavity_free(ib) == 0) then

        prob = 1.0_wp

    else

        prob = real(ncavity_free(ib),wp) / real(ncavity_grid(ib),wp)

    endif

end subroutine


subroutine reset_pcav_molremove(ib, im, crad)

    !use kinds_f90
    use constants_module, only : uout

    !implicit none

    integer, intent(in) :: ib, im

    real(kind = wp), intent(in) :: crad

    integer :: i, j

    real(kind = wp) :: ix, iy, iz, jx, jy, jz, rad
    real(kind = wp) :: xx, yy, zz, rx, ry, rz, rsq

    rad = crad * crad

    ! search over molecule

    do i = 1, cfgs(ib)%mols(im)%natom

        ix = cfgs(ib)%mols(im)%atms(i)%rpos(1)
        iy = cfgs(ib)%mols(im)%atms(i)%rpos(2)
        iz = cfgs(ib)%mols(im)%atms(i)%rpos(3)

        do j = 1, ncavity_grid(ib)

            jx = cfgs(ib)%cavgrid(1,j)
            jy = cfgs(ib)%cavgrid(2,j)
            jz = cfgs(ib)%cavgrid(3,j)

            xx = jx - ix
            yy = jy - iy
            zz = jz - iz

            call pbc_cart_pos_calc(ib, xx, yy, zz)

            rsq = xx * xx + yy * yy + zz * zz 
            
            !if( rsq < rad) cfgs(ib)%cavity(j) = cfgs(ib)%cavity(j) + 1

            if( rsq < rad ) then

                cfgs(ib)%cavity(j) = cfgs(ib)%cavity(j) + 1

                if( ncavity_free(ib) > 0 ) then

                    cfgs(ib)%cavids(ncavity_free(ib)) = 0
                    ncavity_free(ib) = ncavity_free(ib) - 1
!debug<
!AB: NOTE that when inserting a molecule only one cavity it taken (by COM)
!AB: yet, here several cavities can be created (those freed by each atom)
!AB: this discrepancy MUST be eliminated! -- wrt either COM or each atom
!AB: in order to make cavity treatment consistent upon insertion & deletion

                else
                    call cry(uout,'', &
                            "ERROR: reset_pcav_molremove() - free cavity count < 1 while resetting atom back !!!",999)
!debug>
                endif

            end if

            cycle

!AB: the rest of the loop should be obsolete

            if(cfgs(ib)%vec%is_orthogonal) then

                call pbc_ortho_coor_calc(ib, xx, yy, zz)
                rsq = xx * xx + yy * yy + zz * zz 

            else
                rx = cfgs(ib)%vec%invlat(1,1) * xx + cfgs(ib)%vec%invlat(2,1) * yy + cfgs(ib)%vec%invlat(3,1) * zz
                ry = cfgs(ib)%vec%invlat(1,2) * xx + cfgs(ib)%vec%invlat(2,2) * yy + cfgs(ib)%vec%invlat(3,2) * zz
                rz = cfgs(ib)%vec%invlat(1,3) * xx + cfgs(ib)%vec%invlat(2,3) * yy + cfgs(ib)%vec%invlat(3,3) * zz

            ! AB: unified way to apply PBC/MIC
                xx = rx
                yy = ry
                zz = rz
                call pbc_atom_pos_calc(xx, yy, zz)

            ! AB: original code for PBC/MIC
            !xx = rx - anint(rx)
            !yy = ry - anint(ry)
            !zz = rz - anint(rz)

                rx = cfgs(ib)%vec%latvector(1,1) * xx + cfgs(ib)%vec%latvector(2,1) * yy + cfgs(ib)%vec%latvector(3,1) * zz
                ry = cfgs(ib)%vec%latvector(1,2) * xx + cfgs(ib)%vec%latvector(2,2) * yy + cfgs(ib)%vec%latvector(3,2) * zz
                rz = cfgs(ib)%vec%latvector(1,3) * xx + cfgs(ib)%vec%latvector(2,3) * yy + cfgs(ib)%vec%latvector(3,3) * zz

                rsq = rx * rx + ry * ry + rz * rz

            end if

            if( rsq < rad) then

                cfgs(ib)%cavity(j) = cfgs(ib)%cavity(j) + 1

            endif

        enddo

    enddo

end subroutine

!AB: *** the old/original/untested and now obsolete variants of cavity list updating routines - end ***


subroutine check_distance_atoms(ib, distance_atm_max, uselist, shell, ierr, is_reset)

    !use kinds_f90
    use constants_module, only : uout
    use species_module, only : number_of_elements
    use latticevectors_module, only : dcell

    !implicit none

    integer, intent(in) :: ib

    real (kind = wp), intent(inout) :: distance_atm_max(number_of_elements,nconfigs)

    real (kind = wp), intent(in) :: shell

    logical, intent(in) :: uselist, is_reset

    integer, intent(out) :: ierr

    integer, save :: nbox_out = 0

    integer, save :: nshell_out = 0

    character*(32) :: word

    real (kind = wp) :: dist(3), bbb(10)

    integer :: i

    call dcell(cfgs(ib)%vec%latvector, bbb)

    dist(1) = bbb(7) / 2.0_wp
    dist(2) = bbb(8) / 2.0_wp
    dist(3) = bbb(9) / 2.0_wp

    ierr = 0

    if( is_reset ) then 
        nbox_out   = 0
        nshell_out = 0
    end if

    ! first check that the distance is not greater than half the box
    do i = 1, number_of_elements

        if( any( dist(:) < distance_atm_max(i,ib) ) ) then 

            ierr = ierr+1 !call error(347)

            distance_atm_max(i,ib) = minval(dist)-1.e-5_wp

            nbox_out = nbox_out + 1

            if( nbox_out < 11 ) then

                write(word,'(a,i10,a,i5,a)')"(element",i,", cell",ib,")"

                call cry(uout,'(/,1x,a)', &
                     "WARNING: updated max atom step > min(cell dims)/2 - resetting just below it... "&
                     &//trim(word),0)

            else if( nbox_out == 11 ) then

                write(word,'(a,i5,a)')"(more than 10 events, cell",ib,")"

                call cry(uout,'(/,1x,a)', &
                     "WARNING: updated max atom step > min(cell dims)/2 - resetting just below it... "&
                     &//trim(word),0)

            else if( nbox_out == huge(nbox_out) ) then

                 nbox_out = 0

            end if

        end if

    end do

    ! check that it is not greater than verlet shell when using nbr list
    if (uselist) then

        do i = 1, number_of_elements

           if( distance_atm_max(i,ib) > shell ) then 

               ierr = ierr+1 !call error(348)

               distance_atm_max(i,ib) = shell-1.e-5_wp

               nshell_out = nshell_out + 1

               if( nshell_out < 11 ) then

                   write(word,'(a,i10,a,i5,a)')"(element",i,", cell",ib,")"

                   call cry(uout,'(/,1x,a)', &
                        "WARNING: updated max atom step > Verlet shell - resetting just below it... "&
                        &//trim(word),0)

               else if( nshell_out == 11 ) then

                   write(word,'(a,i5,a)')"(more than 10 events, cell",ib,")"

                   call cry(uout,'(/,1x,a)', &
                        "WARNING: updated max atom step > Verlet shell - resetting just below it... "&
                        &//trim(word),0)

               else if( nshell_out == huge(nshell_out) ) then

                 nshell_out = 0

               end if

           end if

        end do

    end if

end subroutine


subroutine check_distance_molecules(ib, distance_mol_max, uselist, shell, ierr, is_reset)

    !use kinds_f90
    use constants_module, only : uout
    use species_module, only : number_of_molecules
    use latticevectors_module, only : dcell

    !implicit none

    integer, intent(in) :: ib

    real (kind = wp), intent(inout) :: distance_mol_max(number_of_molecules)

    real (kind = wp), intent(in) :: shell

    logical, intent(in) :: uselist, is_reset

    integer, intent(out) :: ierr

    integer, save :: nbox_out = 0

    integer, save :: nshell_out = 0

    character*(32) :: word

    integer :: i

    real (kind = wp) :: dist(3), bbb(10)

    call dcell(cfgs(ib)%vec%latvector, bbb)

    dist(1) = bbb(7) / 2.0_wp
    dist(2) = bbb(8) / 2.0_wp
    dist(3) = bbb(9) / 2.0_wp

    ierr = 0

    if( is_reset ) then 
        nbox_out   = 0
        nshell_out = 0
    end if

    ! first check that the distance is not greater than half the box
    do i = 1, number_of_molecules

        if( any( dist(:) < distance_mol_max(i) ) ) then 

            ierr = ierr+1 !call error(347)

            distance_mol_max(i) = minval(dist)-1.e-5_wp

            nbox_out = nbox_out + 1

            if( nbox_out < 11 ) then

                write(word,'(a,i10,a,i5,a)')"(molecule",i,", cell",ib,")"

                call cry(uout,'(/,1x,a)', &
                     "WARNING: updated max molecule step > min(cell dims)/2 - resetting just below it... "&
                     &//trim(word),0)

            else if( nbox_out == 11 ) then

                write(word,'(a,i5,a)')"(more than 10 events, cell",ib,")"

                call cry(uout,'(/,1x,a)', &
                     "WARNING: updated max molecule step > min(cell dims)/2 - resetting just below it... "&
                     &//trim(word),0)

            else if( nbox_out == huge(nbox_out) ) then

                 nbox_out = 0

            end if

        end if

    end do

    ! check that it is not greater than verlet shell when using nbr list
    if (uselist) then

        do i = 1, number_of_molecules

           if( distance_mol_max(i) > shell ) then 

               ierr = ierr+1 !call error(348)

               distance_mol_max(i) = shell-1.e-5_wp

               nshell_out = nshell_out + 1

               if( nshell_out < 11 ) then

                   write(word,'(a,i10,a,i5,a)')"(element",i,", cell",ib,")"

                   call cry(uout,'(/,1x,a)', &
                        "WARNING: updated max molecule step > Verlet shell - resetting just below it... "&
                        &//trim(word),0)

               else if( nshell_out == 11 ) then

                   write(word,'(a,i5,a)')"(more than 10 events, cell",ib,")"

                   call cry(uout,'(/,1x,a)', &
                        "WARNING: updated max molecule step > Verlet shell - resetting just below it... "&
                        &//trim(word),0)

               else if( nshell_out == huge(nshell_out) ) then

                 nshell_out = 0

               end if

           end if

        end do

    end if

end subroutine


subroutine broadcast_atom_pos(ib, im, i)

    !use kinds_f90
    use comms_mpi_module, only : master, is_serial, gsync, msg_bcast
    !use parallel_loop_module, only : master

    !implicit none

    integer, intent(in) :: ib, im, i

    real (kind = wp) :: buf(3)

    if( is_serial ) return

    call gsync

    ! load positions into buffers
    if (master) then

        buf(1) = cfgs(ib)%mols(im)%atms(i)%rpos(1)
        buf(2) = cfgs(ib)%mols(im)%atms(i)%rpos(2)
        buf(3) = cfgs(ib)%mols(im)%atms(i)%rpos(3)

    endif

    ! send to other nodes
    call msg_bcast(buf, 3)


    ! copy back
    cfgs(ib)%mols(im)%atms(i)%rpos(1) = buf(1)
    cfgs(ib)%mols(im)%atms(i)%rpos(2) = buf(2)
    cfgs(ib)%mols(im)%atms(i)%rpos(3) = buf(3)
    
end subroutine

subroutine broadcast_atom_type(ib, im, i)

    !use kinds_f90
    use comms_mpi_module, only : master, is_serial, gsync, msg_bcast
    !use parallel_loop_module, only : master

    !implicit none

    integer, intent(in) :: ib, im, i

    real (kind = wp) :: buf(5)

    return

    if( is_serial ) return

    call gsync

    ! load atom species data to buffer
    if (master) then

        buf(1) = real(cfgs(ib)%mols(im)%atms(i)%atlabel, wp)
        buf(2) = real(cfgs(ib)%mols(im)%atms(i)%atype, wp)
        buf(3) = cfgs(ib)%mols(im)%atms(i)%charge
        buf(4) = cfgs(ib)%mols(im)%atms(i)%mass
        buf(5) = real(cfgs(ib)%mols(im)%atms(i)%site, wp)

    endif

    ! send to other nodes
    call msg_bcast(buf, 5)


    ! copy back
    cfgs(ib)%mols(im)%atms(i)%atlabel = nint(buf(1))
    cfgs(ib)%mols(im)%atms(i)%atype = nint(buf(2))
    cfgs(ib)%mols(im)%atms(i)%charge = buf(3)
    cfgs(ib)%mols(im)%atms(i)%mass = buf(4)
    cfgs(ib)%mols(im)%atms(i)%site = nint(buf(5))

end subroutine

subroutine broadcast_molecule_pos(ib, im)

    !use kinds_f90
    use comms_mpi_module, only : master, is_serial, gsync, msg_bcast
    !use parallel_loop_module, only : master

    !implicit none

    integer, intent(in) :: ib, im

    integer :: i, i2, i3

    real (kind = wp) :: buf(cfgs(ib)%mols(im)%natom*3)

    return

    if( is_serial ) return

    i2 = cfgs(ib)%mols(im)%natom
    i3 = cfgs(ib)%mols(im)%natom * 2

    call gsync

    ! load positions into buffers
    if (master) then

        do i = 1, cfgs(ib)%mols(im)%natom

            buf(i) = cfgs(ib)%mols(im)%atms(i)%rpos(1)
            buf(i+i2) = cfgs(ib)%mols(im)%atms(i)%rpos(2)
            buf(i+i3) = cfgs(ib)%mols(im)%atms(i)%rpos(3)

        enddo

    endif

    ! send to other nodes
    call msg_bcast(buf, cfgs(ib)%mols(im)%natom*3)


    ! copy back
    do i = 1, cfgs(ib)%mols(im)%natom

        cfgs(ib)%mols(im)%atms(i)%rpos(1) = buf(i)
        cfgs(ib)%mols(im)%atms(i)%rpos(2) = buf(i+i2)
        cfgs(ib)%mols(im)%atms(i)%rpos(3) = buf(i+i3)

    enddo

end subroutine

subroutine broadcast_molecule_type(ib, im)

    !use kinds_f90
    use comms_mpi_module, only : master, is_serial, gsync, msg_bcast
    !use parallel_loop_module, only : master

    !implicit none

    integer, intent(in) :: ib, im

    integer :: i, i2, i3, i4, i5

    real (kind = wp) :: buf(cfgs(ib)%mols(im)%natom*3)

    return

    if( is_serial ) return

    i2 = cfgs(ib)%mols(im)%natom
    i3 = cfgs(ib)%mols(im)%natom * 2
    i4 = cfgs(ib)%mols(im)%natom * 3
    i5 = cfgs(ib)%mols(im)%natom * 4

    call gsync

    if (master) then

        do i = 1, cfgs(ib)%mols(im)%natom

            buf(i) = real(cfgs(ib)%mols(im)%atms(i)%atlabel, wp)
            buf(i+i2) = real(cfgs(ib)%mols(im)%atms(i)%atype, wp)
            buf(i+i3) = cfgs(ib)%mols(im)%atms(i)%charge
            buf(i+i4) = cfgs(ib)%mols(im)%atms(i)%mass
            buf(i+i5) = real(cfgs(ib)%mols(im)%atms(i)%site, wp)

        enddo

    endif

    ! send to other nodes
    call msg_bcast(buf, cfgs(ib)%mols(im)%natom*5)

    ! copy back
    do i = 1, cfgs(ib)%mols(im)%natom

        cfgs(ib)%mols(im)%atms(i)%atlabel = nint(buf(i))
        cfgs(ib)%mols(im)%atms(i)%atype = nint(buf(i+i2))
        cfgs(ib)%mols(im)%atms(i)%charge = buf(i+i3)
        cfgs(ib)%mols(im)%atms(i)%mass = buf(i+i4)
        cfgs(ib)%mols(im)%atms(i)%site = nint(buf(i+i5))

    enddo

end subroutine

subroutine broadcast_allatoms(ib)

    !use kinds_f90
    use comms_mpi_module, only : idnode, is_serial, gsync, msg_bcast

    !implicit none

    integer, intent(in) :: ib

    integer :: i, im, i1, i2, i3

    real (kind = wp) :: buf(cfgs(ib)%number_of_atoms*3)

    return

    if( is_serial ) return

    i1 = 1
    i2 = i1 + cfgs(ib)%number_of_atoms
    i3 = i2 + cfgs(ib)%number_of_atoms

    call gsync

    ! load positions into buffers
    if (idnode == 0) then

        do im = 1, cfgs(ib)%num_mols

            do i = 1, cfgs(ib)%mols(im)%natom

                buf(i1) = cfgs(ib)%mols(im)%atms(i)%rpos(1)
                buf(i2) = cfgs(ib)%mols(im)%atms(i)%rpos(2)
                buf(i3) = cfgs(ib)%mols(im)%atms(i)%rpos(3)

                i1 = i1 + 1
                i2 = i2 + 1
                i3 = i3 + 1

            enddo

        enddo

    endif

    ! send to other nodes
    call msg_bcast(buf, cfgs(ib)%number_of_atoms*3)

    i1 = 1
    i2 = i1 + cfgs(ib)%number_of_atoms
    i3 = i2 + cfgs(ib)%number_of_atoms

    ! copy back
    do im = 1, cfgs(ib)%num_mols

        do i = 1, cfgs(ib)%mols(im)%natom

            cfgs(ib)%mols(im)%atms(i)%rpos(1) = buf(i1)
            cfgs(ib)%mols(im)%atms(i)%rpos(2) = buf(i2)
            cfgs(ib)%mols(im)%atms(i)%rpos(3) = buf(i3)

            i1 = i1 + 1
            i2 = i2 + 1
            i3 = i3 + 1

        enddo

    enddo

end subroutine



!TU: Procedures for calculating the separation between particles 


!TU: Added by me...
!> Returns the separation vector between two position vectors in the specified configuration, wrapping
!> the positions into the primary image via the appropriate periodic boundary conditions (including 
!> 'orthogonal' and 'slit' geometries). The separation is from position 1 to 2
function pbc_sepvec(ib, pos1, pos2)

    implicit none

        !> Configuration identifier (index in '`cfgs`' array)
    integer, intent(in) :: ib

        !> Position vector 1
    real(wp), intent(in) :: pos1(3)

        !> Position vector 2
    real(wp), intent(in) :: pos2(3)

    real(wp) :: pbc_sepvec(3)
    
        ! Separation vector in absolute and fractional coordinates
    real(wp) :: r(3), fr(3)

    r = pos2 - pos1

    if(cfgs(ib)%vec%is_orthogonal) then

        ! For the orthogonal case...

        ! Wrap the components of r to the primary image (pbc_ortho_coor_calc takes absolute coordinates as input)
        call pbc_ortho_coor_calc(ib, r(1), r(2), r(3))
       
        pbc_sepvec = r

    else

        ! For the general/non-orthogonal case...

        ! Calculate seperation vector in fractional coordinates
        fr(:) = cfgs(ib)%vec%invlat(1,:) * r(1) + cfgs(ib)%vec%invlat(2,:) * r(2) + cfgs(ib)%vec%invlat(3,:) * r(3)

        ! Wrap the components of fr to the primary image (pbc_atom_pos_calc takes fractional coordinates as input)
        call pbc_atom_pos_calc(fr(1), fr(2), fr(3))

        ! Convert back to absolute coordinates
        pbc_sepvec(:) = cfgs(ib)%vec%latvector(1,:) * fr(1) &
                      + cfgs(ib)%vec%latvector(2,:) * fr(2) &
                      + cfgs(ib)%vec%latvector(3,:) * fr(3)

    end if

end function pbc_sepvec




!TU: Added by me...
!> Returns the squared separation between two position vectors in the specified configuration, wrapping
!> the positions into the primary image via the appropriate periodic boundary conditions (including 
!> 'orthogonal' and 'slit' geometries) 
real(wp) function pbc_sepsq(ib, pos1, pos2)

    implicit none

        !> Configuration identifier (index in '`cfgs`' array)
    integer, intent(in) :: ib

        !> Position vector
    real(wp), intent(in) :: pos1(3)

        !> Position vector
    real(wp), intent(in) :: pos2(3)

    real(wp) :: sepvec(3)

    sepvec = pbc_sepvec(ib, pos1, pos2)
    
    pbc_sepsq = sepvec(1)*sepvec(1) + sepvec(2)*sepvec(2) + sepvec(3)*sepvec(3)

end function pbc_sepsq




!TU: Added by me...
!> Returns the separation vector between two atoms in the specified configuration, wrapping
!> their positions into the primary image via the appropriate periodic boundary conditions (including 
!> 'orthogonal' and 'slit' geometries). The separation is from position 1 to 2
function pbc_sepvec_atoms(ib, m1, a1, m2, a2)

    implicit none

        !> Configuration identifier (index in '`cfgs`' array)
    integer, intent(in) :: ib

        !> Molecule index for atom 1
    integer, intent(in) :: m1

        !> Atom index for atom 1
    integer, intent(in) :: a1

        !> Molecule index for atom 2
    integer, intent(in) :: m2

        !> Atom index for atom 2
    integer, intent(in) :: a2

    real(wp) :: pbc_sepvec_atoms(3)

    pbc_sepvec_atoms = pbc_sepvec(ib, cfgs(ib)%mols(m1)%atms(a1)%rpos, cfgs(ib)%mols(m2)%atms(a2)%rpos)

end function pbc_sepvec_atoms




!TU: Added by me...
!> Returns the squared separation between two atoms in the specified configuration, wrapping
!> their positions into the primary image via the appropriate periodic boundary conditions (including 
!> 'orthogonal' and 'slit' geometries) 
real(wp) function pbc_sepsq_atoms(ib, m1, a1, m2, a2)

    implicit none

        !> Configuration identifier (index in '`cfgs`' array)
    integer, intent(in) :: ib

        !> Molecule index for atom 1
    integer, intent(in) :: m1

        !> Atom index for atom 1
    integer, intent(in) :: a1

        !> Molecule index for atom 2
    integer, intent(in) :: m2

        !> Atom index for atom 2
    integer, intent(in) :: a2

    pbc_sepsq_atoms = pbc_sepsq(ib, cfgs(ib)%mols(m1)%atms(a1)%rpos, cfgs(ib)%mols(m2)%atms(a2)%rpos)

end function pbc_sepsq_atoms




end module
