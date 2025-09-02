program interp_ssh_contorno
  use netcdf
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)
  integer :: ncid, lat_varid, lon_varid, ssh_varid, time_varid,dimid
  integer :: nlat, nlon, nnos,ncols,ncontorno,status
  integer :: i, j, n, m, d, g, L, ntime, nn, mn
  real(dp), allocatable :: ssh(:,:,:), lat(:), lon(:), mesh(:,:), cont(:,:), ssh_interp(:)
  real(dp), allocatable :: dif(:), dif1(:), cont_vel(:,:), time(:), &
          rcontorno(:,:)
  integer, allocatable :: nodos_cont(:)
  character(len=100) :: ncfile_mom, linha,xlon,ylat,ttime,hsea, &
          malha_2d,contorno,line,z_col_deg
  character(len=30) :: dummy,dimname

  namelist /ssh_config/ ncfile_mom,xlon,ylat,ttime,hsea
  open(unit=10, file="namelist_pre.nml", status="old")
  read(10, nml=ssh_config)
  close(10)
  z_col_deg='../resultados/zz_coluna_deg.dat'
  malha_2d='../resultados/malha_2d_deg.dat'
  contorno='../resultados/c_contorno.cli'
 
  ! === Definição do arquivo NetCDF ===
  print*,"LENDO ARQUIVO NETCDF: ",ncfile_mom
  call check(nf90_open(ncfile_mom, NF90_NOWRITE, ncid))
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

  call check(nf90_close(ncid))
  do i = 1, nlon
    do j = 1, nlat
      do L = 1, ntime
          if (isnan(ssh(i,j,L)) .or.ssh(i,j,L) .gt. 1000.0 ) ssh(i,j,L) = 1000.0_dp
      end do
    end do
  end do

  ! Primeiro contar o número de colunas
  open(unit=10, file=trim(z_col_deg), status='old', action='read')
  read(10, '(a)', iostat=status) line ! Ler primeira linha para contar colunas
  ncols = 0
  do
    read(line, *, iostat=status) (dummy, i=1,ncols+1)
    if (status /= 0) exit
    ncols = ncols + 1
  end do
  rewind(10)
  print*, "Número de colunas em",trim(z_col_deg),":",ncols
  ! Agora contar número de linhas
  nnos = 0
  do
    read(10, *, iostat=status)
    if (status /= 0) exit
    nnos = nnos + 1
  end do
  rewind(10)
  print*, "Número de linhas em",trim(z_col_deg),":", nnos

  allocate(mesh(nnos, ncols))
  do i = 1, nnos
    read(10, *) mesh(i,:)
  end do
  close(10)

  ! === Leitura do arquivo de contorno ===
  print*, "LENDO ARQUIVO DE CONTORNO: ", trim(contorno)
  open(unit=11, file=trim(contorno), status='old')
  ! Contar número de linhas
  ncontorno = 0
  do
    read(11, *, iostat=status)
    if (status /= 0) exit
    ncontorno = ncontorno + 1
  end do
  rewind(11)
  print*, "Número de linhas em c_contorno_BC.cli: ",ncontorno

  ! Ler dados do contorno
  allocate(cont(ncontorno, 13), rcontorno(ncontorno, 13))
  do i = 1, ncontorno
    read(11, *) rcontorno(i,:)
    cont(i,:) = int(rcontorno(i,:))
  end do
  close(11)
  print*,"Salvando condições de contorno"
  g = 0
  do i = 1, ncontorno
    if ((int(cont(i,2)) .eq. 6) .or. cont(i,2) .eq. 4 .or. cont(i,2) .eq. 5) then
            g = g + 1
    end if
  end do
  allocate(nodos_cont(g))
  g = 0
  do i = 1, ncontorno
    if ((int(cont(i,2)) .eq. 6) .or. cont(i,2) .eq. 4 .or. cont(i,2) .eq. 5) then 
      g = g + 1
      nodos_cont(g) = int(cont(i,12))
    end if
  end do
  allocate(cont_vel(ntime+1, g+1))

  print*, "INTERPOLANDO A SSH PARA OS PONTOS DE GRADE"
 
  cont_vel(1,1) = g
  do i = 1, g
    cont_vel(1,i+1) = real(nodos_cont(i),dp)
  end do
  print*,ntime,g
  do L = 1, ntime
   PRINT*, "PROCESSANDO O TEMPO: ", (L-1)*3600, 'T= ',L-1
    do i = 1, g
      ! Encontrar a posição mais próxima no grid lat/lon
      allocate(dif(size(lat)))
      allocate(dif1(size(lon)))
      do m = 1, size(lat)
        dif(m) = abs(mesh(nodos_cont(i),2) - lat(m))
      end do
      do m = 1, size(lon)
        dif1(m) = abs(mesh(nodos_cont(i),1) - lon(m))
      end do
      j = minloc(dif,1 )
      n = minloc(dif1, 1)
      cont_vel(L+1,i+1)=ssh(n,j,L)
      deallocate(dif, dif1)
    end do

    cont_vel(L+1,1) = real((L-1)*3600, dp)  ! tempo em segundos
  end do

  ! === Escrita final ===
  open(unit=12, file='../resultados/nivel_cont.dat', status='replace')
  do i = 1, ntime+1
    write(12,*) cont_vel(i,:)
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

