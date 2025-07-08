! Programa Fortran para interpolar SSH (nível do mar) em pontos de uma malha
! e extrair os valores apenas nos pontos de contorno de tipo 5 (nível)

program interp_ssh_contorno
  use netcdf
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)
  integer :: ncid, lat_varid, lon_varid, ssh_varid, time_varid,dimid
  integer :: nlat, nlon, nnos
  integer :: i, j, n, m, d, g, L, ntime, nn, mn
  real(dp), allocatable :: ssh(:,:,:), lat(:), lon(:), mesh(:,:), cont(:,:), ssh_interp(:)
  real(dp), allocatable :: mesh_ll(:,:), dif(:), dif1(:), cont_vel(:,:), tempo(:),time(:)
  integer, allocatable :: nodos_cont(:)
  character(len=100) :: ncfile, linha,xlon,ylat,ttime,hsea
  character(len=30) :: dummy,dimname

  ! Ajuste conforme os tamanhos reais esperados
  integer, parameter :: nmax = 100000, tmax = 1000
  namelist /ssh_config/ ncfile,xlon,ylat,ttime,hsea
  open(unit=10, file="namelist_pre.nml", status="old")
  read(10, nml=ssh_config)
  close(10)

  ! === Definição do arquivo NetCDF ===
!  ncfile = '../dados/20070101.ocean_telemac_recorte.nc'

  print*,"LENDO ARQUIVO NETCDF: ",ncfile
  call check(nf90_open(ncfile, NF90_NOWRITE, ncid))
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

  allocate(ssh(nlon, nlat,ntime))
  
  call check(nf90_inq_varid(ncid, hsea, ssh_varid))
  call check(nf90_get_var(ncid, ssh_varid, ssh))

  ! === Leitura da malha ===
  print*, "LENDO A MALHA: ", "malha_deg.dat"

  open(unit=10, file='../resultados/malha_deg.dat', status='old')
  allocate(mesh(nmax, 3))
  i = 0
  do
    read(10, *, iostat=n) mesh(i+1,1:3)
    if (n /= 0) exit
    i = i + 1
  end do
  close(10)
  allocate(mesh_ll(i, 3))

  ! === Leitura do arquivo de contorno ===
  print*, "LENDO ARQUIVO DE CONTORNO: ", "c_contorno_BC.cli"
  open(unit=11, file='../dados/c_contorno_BC.cli', status='old')
  allocate(cont(nmax, 20))
  i = 0
  do
    read(11, *, iostat=n) cont(i+1,1:20)
    if (n /= 0) exit
    i = i + 1
  end do
  close(11)

  allocate(nodos_cont(i))
  g = 0
  do n = 1, i
    if (nint(cont(n,1)) == 5) then
      g = g + 1
      nodos_cont(g) = int(cont(n,12))
    end if
  end do

  ! === Interpolação do ssh para os pontos da malha ===
  allocate(ssh_interp(g))
  allocate(cont_vel(tmax+1, g+1))
  allocate(tempo(tmax))

  print*, "INTERPOLANDO A SSH PARA OS PONTOS DE GRADE"
 
  cont_vel(1,1) = g
  print*,tmax,g,nodos_cont(1)
  do i = 1, g
    cont_vel(1,i+1) = mesh(nodos_cont(i),1)
  end do

  do L = 1, ntime
   PRINT*, "PROCESSANDO O TEMPO: ", (L-1)*3600, 'T= ',L-1
    do i = 1, g
      ! Encontrar a posição mais próxima no grid lat/lon
      allocate(dif(size(lat)))
      allocate(dif1(size(lon)))
      do m = 1, size(lat)
        dif(m) = abs(mesh_ll(nodos_cont(i),2) - lat(m))
      end do
      do m = 1, size(lon)
        dif1(m) = abs(mesh_ll(nodos_cont(i),1) - lon(m))
      end do
      j = minloc(dif,1 )
      n = minloc(dif1, 1)
      ssh_interp(i) = ssh(j,n,L)  ! Nota: inversão de j/n pode ser necessária
      deallocate(dif, dif1)
    end do

    cont_vel(L+1,1) = real((L-1)*3600, dp)  ! tempo em segundos
    do i = 1, g
      cont_vel(L+1,i+1) = ssh_interp(i)
    end do
  end do

  ! === Escrita final ===
  open(unit=12, file='../resultados/nivel_cont.dat', status='replace')
  do i = 1, ntime+1
    write(12,*) (cont_vel(i,j), j=1,g+1)
  end do
  close(12)

  print *, 'Arquivo nivel_cont.dat gerado com sucesso.'

contains
  subroutine check(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 1
    end if
  end subroutine

end program interp_ssh_contorno

