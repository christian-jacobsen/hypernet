!------------------------------------------------------------------------------!
!> This program converts the energy levels from [cm^-1] to [eV]
  program convert

    implicit none

    integer, parameter             :: hpc_4 = 4
    integer, parameter             :: hpc_8 = 8
    integer(kind=hpc_4), parameter :: in1  = 10
    integer(kind=hpc_4), parameter :: out1 = 20
    integer(kind=hpc_4)            :: i, gi_int
    real(kind=hpc_8), parameter    :: ukb = 1.380658d-23
    real(kind=hpc_8), parameter    :: ue  = 1.602191d-19
    real(kind=hpc_8), parameter    :: fac = 1.43876866d0*ukb/ue
    real(kind=hpc_4)               :: gi, Ei

    open(unit=in1,file='levels',status='old')
    open(unit=out1,file='levels_new',status='unknown')

    do 
       read(in1,*,end=20)gi_int,Ei
       gi = gi_int*1.d0
       Ei = Ei*fac 
       write(out1,10)gi,Ei
    enddo
 20 rewind(in1)

    close(in1)
    close(out1)

10 format(100e14.6)

  end program convert
!------------------------------------------------------------------------------!
