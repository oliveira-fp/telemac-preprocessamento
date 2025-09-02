program interpolacao_vento
  use netcdf
  implicit none

  ! Parâmetros
  integer, parameter :: dp = selected_real_kind(15)
  integer :: dt
  integer :: ncid, lat_varid, lon_varid, u_varid, v_varid, time_varid,dimid
  integer :: i, j, k, status, nlat, nlon, ntime, nnos
  real(dp), allocatable :: lat(:), lon(:), time(:)
  real(dp), allocatable :: u10m(:,:,:), v10m(:,:,:)
  real(dp), allocatable :: xutm(:), yutm(:), lon_malha(:), lat_malha(:)
  character(len=100) :: ncfile_eta,linha,dimname,uwind,vwind,xlon,ylat,ttime,malha_2d
  real(dp), allocatable :: u_interp(:), v_interp(:)
  namelist /atmosferico/ ncfile_eta,dt,xlon,ylat,ttime,uwind,vwind

  open(unit=10, file="namelist_pre.nml", status="old")
  read(10, nml=atmosferico)
  close(10)
  malha_2d='../resultados/malha_2d_deg.dat'

  ! Arquivos de entrada
  ! === Abrir netCDF ===
  print*,"LENDO ARQUIVO NETCDF: ",ncfile_eta
  call check(nf90_open(ncfile_eta, NF90_NOWRITE, ncid))
  call check(nf90_inq_varid(ncid, ylat, lat_varid))
  call check(nf90_inq_varid(ncid, xlon, lon_varid))
  call check(nf90_inq_varid(ncid, ttime, time_varid))

  call check(nf90_inquire_dimension(ncid,lat_varid,dimname, nlat))
  call check(nf90_inquire_dimension(ncid,lon_varid,dimname, nlon))
  call check(nf90_inquire_dimension(ncid,time_varid,dimname, ntime))

  allocate(lat(nlat), lon(nlon), time(ntime))
  call check(nf90_get_var(ncid, lat_varid, lat))
  call check(nf90_get_var(ncid, lon_varid, lon))
  call check(nf90_get_var(ncid, time_varid, time))
 
  ! === Ler malha.dat ===
  print*,"LENDO A MALHA: ", trim(malha_2d)
  open(unit=10, file=trim(malha_2d), status='old')
  nnos = 0
  do
    read(10, '(A)', iostat=status) linha
    if (status /= 0) exit
    nnos = nnos + 1
  end do
  rewind(10)

  allocate(xutm(nnos), yutm(nnos), lon_malha(nnos), lat_malha(nnos))
  do i = 1, nnos
    read(10, *) j, lon_malha(i), lat_malha(i)
  end do
  close(10)

  ! === Alocar saídas ===
  allocate(u_interp(nnos), v_interp(nnos))

  ! === Abrir saída ===
  PRINT*, "ESCREVENDO AS SAÍDAS"
  open(unit=20, file='../resultados/Malha_Vento_ETA03.dat', status='replace')
  write(20, '(I6, 1X, *(I6, 1X))') nnos, (i, i=1, nnos)

  ! === Loop de tempo ===
  allocate(u10m(nlon, nlat, 1), v10m(nlon, nlat, 1))  ! ler 1 tempo por vez
  do k = 1, ntime
    print *, 'Tempo: ', time(k)*dt, "T= ",k

    call check(nf90_inq_varid(ncid, uwind, u_varid))
    call check(nf90_inq_varid(ncid, vwind, v_varid))

    call check(nf90_get_var(ncid, u_varid, u10m, start=[1,1,k], count=[nlon,nlat,1]))
    call check(nf90_get_var(ncid, v_varid, v10m, start=[1,1,k], count=[nlon,nlat,1]))

    do i = 1, nnos
      call interp_bilinear(lat, lon, u10m(:,:,1), lat_malha(i), lon_malha(i), u_interp(i))
      call interp_bilinear(lat, lon, v10m(:,:,1), lat_malha(i), lon_malha(i), v_interp(i))
    end do

    write(20, '(F10.0, 1X, *(F7.2, 1X))') time(k)*dt, (u_interp(i), i=1, nnos)
    write(20, '(F10.0, 1X, *(F7.2, 1X))') time(k)*dt, (v_interp(i), i=1, nnos)
  end do

  close(20)
  call check(nf90_close(ncid))

  print *, 'Arquivo salvo: Malha_Vento_ETA03.dat'
contains

  subroutine check(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 1
    end if
  end subroutine

  subroutine interp_bilinear(lat, lon, grid, latp, lonp, val)
    real(dp), intent(in) :: lat(:), lon(:), grid(size(lon), size(lat))
    real(dp), intent(in) :: latp, lonp
    real(dp), intent(out) :: val
    integer :: i, j
    real(dp) :: t, s, f00, f01, f10, f11

    do i = 1, size(lat)-1
      if (latp >= lat(i) .and. latp <= lat(i+1)) exit
    end do
    do j = 1, size(lon)-1
      if (lonp >= lon(j) .and. lonp <= lon(j+1)) exit
    end do

    t = (latp - lat(i)) / (lat(i+1) - lat(i))
    s = (lonp - lon(j)) / (lon(j+1) - lon(j))

    f00 = grid(j,i)
    f01 = grid(j,i+1)
    f10 = grid(j+1,i)
    f11 = grid(j+1,i+1)

    val = (1-t)*(1-s)*f00 + t*(1-s)*f01 + (1-t)*s*f10 + t*s*f11
  end subroutine
end program

