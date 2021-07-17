! -*-f90-*-
!-----------------------------------------------------------------------------!
!> This module stores arrays used for the evaluation of group thermodynamic properties
  module general_file

    implicit none

    private 
    public :: get_nb_lines, check_file, get_unit

    integer, parameter :: hpc_4 = 4
    integer, parameter :: hpc_8 = 8

    contains 

      !---------------------------------------------------!
      !> This subroutine computes the number of lines of a file
      function get_nb_lines(infile)

        integer(kind=hpc_4), intent(in) :: infile
        integer(kind=hpc_4)             :: get_nb_lines

        get_nb_lines = 0
        do
          read(infile,*,end=30)
          get_nb_lines = get_nb_lines + 1
        end do
   30   rewind(infile)

      end function get_nb_lines

      !---------------------------------------------------!
      !> This subroutine checks if a given file has been opened correctly 
      subroutine check_file(ios, filename)

        integer(kind=hpc_4), intent(in) :: ios
        character*(*), intent(in)       :: filename

        if (ios.ne.0) then 
          write(*,5)'PLATO library -> ERROR!!'
          write(*,5)'File "'//trim(filename)//'" not found...'
          print*
          stop
        endif

5  format(a)

      end subroutine check_file 

      !---------------------------------------------------!
      !> This function returns a 'file unit' integer number by avoiding to select 
      !! a file which is already opened 
      function get_unit() result(iunit)

        integer(kind=hpc_4) :: iunit

        integer(kind=hpc_4) :: i, ios
        logical             :: lopen

        iunit = 0

        do i = 1, 99

           if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

              inquire (unit=i,opened=lopen,iostat=ios)

              if (ios.eq.0) then
                 if (.not.lopen) then
                    iunit = i
                    return
                 endif
              endif

           endif

        enddo

      end function get_unit

  end module general_file
!-----------------------------------------------------------------------------!
!> This program 
  program main

    use general_file,            only: get_unit, get_nb_lines, check_file

    implicit none

    !-----------------------!
    !> Precision parameters
    integer, parameter                             :: hpc_4 = 4
    integer, parameter                             :: hpc_8 = 8

    integer(kind=hpc_4), parameter                 :: s_str = 40
    integer(kind=hpc_4), parameter                 :: l_str = 300

    !-----------------------!
    !> Mixture related data
    integer(kind=hpc_4)                            :: ng                  !< number of groups
    integer(kind=hpc_4)                            :: np                  !< number of points in look-up-tables

    !-----------------------!
    !> Thermodynamic tables 
    ! - internal energy
    ! - specific heat
    ! - partition function
    real(kind=hpc_8), dimension(:,:), allocatable  :: eint_table          !< internal energy table 
    real(kind=hpc_8), dimension(:,:), allocatable  :: cvint_table         !< internal specific heat table
    real(kind=hpc_8), dimension(:,:), allocatable  :: Qint_table          !< internal partition function table

    integer(kind=hpc_4)                             :: astat, file_unit, ios
    integer(kind=hpc_4)                             :: g, p
    integer(kind=hpc_4)                             :: k, kp
    integer(kind=hpc_4)                             :: nd, nexc_mass, nexc_en
    integer(kind=hpc_4)                             :: dummy_int
    character(len=s_str)                            :: solver    
    character(len=s_str)                            :: g_char, ng_char, specie
    character(len=l_str)                            :: path, datafile

    solver = 'CONVERTER'
    specie = 'O2'
    path   = '/home/zanardi/WORKSPACE/CFD/HyperNet/hypernet/hypernet/applications/shock1d/dataGen/database/grouping/'//trim(specie)

    ! Initialize number of species, groups and temperatures
    ng = 3
    nd = ng

    nexc_mass = 0 
    do k = 1,ng - 1
       do kp = k + 1,ng 
          nexc_mass = nexc_mass + 1
       enddo
    enddo

    nexc_en = 0 
    do k = 1,ng
       do kp = 1,ng 
          nexc_en = nexc_en + 1
       enddo
    enddo

    ! Common factor
    write(ng_char,10)ng
    ng_char = trim(adjustl(ng_char))

    !-------------------!
    ! Load group data:
    ! - temperature table
    ! - internal energy, specific heat and partition function tables

    !-------------------!
    ! Temperature table (ascii format)
    datafile = trim(path)//'/thermo/'//trim(ng_char)//'g/Tt_table.dat'

    file_unit = get_unit()
    open(unit=file_unit,file=trim(datafile),status='old',action='read',iostat=ios)

    call check_file(ios, trim(datafile))

    np = get_nb_lines(file_unit)

    close(file_unit)

    ! Allocate memory 
    ! - Thermodynamics
    allocate(eint_table(np,ng),stat=astat)
    if (astat.ne.0) stop "Error while allocating: 'eint_table'"

    allocate(cvint_table(np,ng),stat=astat)
    if (astat.ne.0) stop "Error while allocating: 'cvint_table'"

    allocate(Qint_table(np,ng),stat=astat)
    if (astat.ne.0) stop "Error while allocating: 'Qint_table'"

    !-------------------!
    ! Tables for thermodynamic properties (ascii format)
    print*
    write(*,5)trim(solver)//':: -> loading GROUP thermodynamic data...'

    ! - internal energy 
    datafile = trim(path)//'/thermo/'//trim(ng_char)//'g/'//trim(specie)//'_eint.dat'

    file_unit = get_unit()
    open(unit=file_unit,file=trim(datafile),status='old',action='read',iostat=ios,form='formatted')

    call check_file(ios, trim(datafile)) 

    do g = 1,ng
       read(file_unit,*) dummy_int
       do p = 1,np
          read(file_unit,*) eint_table(p,g)
       enddo
    enddo

    close(file_unit)

    ! - internal specific heat
    datafile = trim(path)//'/thermo/'//trim(ng_char)//'g/'//trim(specie)//'_cvint.dat'

    file_unit = get_unit()
    open(unit=file_unit,file=trim(datafile),status='old',action='read',iostat=ios,form='formatted')

    call check_file(ios, trim(datafile)) 

    do g = 1,ng 
       read(file_unit,*) dummy_int
       do p = 1,np
          read(file_unit,*) cvint_table(p,g)
       enddo
    enddo

    close(file_unit)

    ! - internal partition function
    datafile = trim(path)//'/thermo/'//trim(ng_char)//'g/'//trim(specie)//'_qint.dat'

    file_unit = get_unit()
    open(unit=file_unit,file=trim(datafile),status='old',action='read',iostat=ios,form='formatted')

    call check_file(ios, trim(datafile)) 

    do g = 1,ng 
       read(file_unit,*)dummy_int
       do p = 1,np
          read(file_unit,*) Qint_table(p,g)
       enddo
    enddo

    close(file_unit)

    print*
    write(*,5)'DONE!'

    !-------------------!
    ! Tables for thermodynamic properties (binary format)
    print*
    write(*,5)trim(solver)//':: -> writing GROUP thermodynamic data in BINARY format...'

    do g = 1,ng

       write(g_char,10)g
       g_char = adjustl(g_char)

       ! - internal energy 
       datafile = trim(path)//'/thermo/'//trim(ng_char)//'g/'//trim(specie)//'_'//trim(g_char)//'_eint.dat'

       file_unit = get_unit()
       open(unit=file_unit,file=trim(datafile),status='replace',action='write',iostat=ios,form='unformatted', & 
          & access='direct',recl=np*hpc_8)

       call check_file(ios, trim(datafile)) 

       write(file_unit,rec=1)(eint_table(p,g), p = 1,np)

       close(file_unit)

       ! - internal specific heat
       datafile = trim(path)//'/thermo/'//trim(ng_char)//'g/'//trim(specie)//'_'//trim(g_char)//'_cvint.dat'

       file_unit = get_unit()
       open(unit=file_unit,file=trim(datafile),status='replace',action='write',iostat=ios,form='unformatted', & 
          & access='direct',recl=np*hpc_8)

       call check_file(ios, trim(datafile)) 

       write(file_unit,rec=1)(cvint_table(p,g), p = 1,np)

       close(file_unit)

       ! - internal partition function
       datafile = trim(path)//'/thermo/'//trim(ng_char)//'g/'//trim(specie)//'_'//trim(g_char)//'_qint.dat'

       file_unit = get_unit()
       open(unit=file_unit,file=trim(datafile),status='replace',action='write',iostat=ios,form='unformatted', & 
          & access='direct',recl=np*hpc_8)

       call check_file(ios, trim(datafile)) 

       write(file_unit,rec=1)(qint_table(p,g), p = 1,np)

       close(file_unit)

    enddo

    print*
    write(*,5)'DONE!'

    !-----------------------!
    ! Dellocate memory
    deallocate(eint_table,stat=astat)
    if (astat.ne.0) stop "Error while deallocating: 'eint_table'"

    deallocate(cvint_table,stat=astat)
    if (astat.ne.0) stop "Error while deallocating: 'cvint_table'"

    deallocate(Qint_table,stat=astat)
    if (astat.ne.0) stop "Error while deallocating: 'Qint_table'"

5  format(a)
10 format(i4)

      end program main
!-----------------------------------------------------------------------------!
