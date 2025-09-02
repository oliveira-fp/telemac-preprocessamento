program interp_velocidades
  use netcdf
  implicit none

  ! Parâmetros
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer :: profu = 10  ! Número de níveis da malha
  integer :: ncid, lat_varid, lon_varid, u_varid, v_varid, depth_varid, time_varid
  integer :: nlat, nlon, ndepth, ntime, nnos, ncontorno,nlon2
  integer :: i, j, k, L, m, n, t, u, cont, contt, g, d, status, ncols
  real(dp) :: tempo = 0.0_dp
  real(dp), allocatable :: lat(:), lon(:), depth(:), time(:),lon2(:)
  real(dp), allocatable :: uo(:,:,:,:), vo(:,:,:,:)
  real(dp), allocatable :: mesh(:,:), meshu1(:,:), mesh_utm(:,:)
  real(dp), allocatable :: mesh_u(:,:,:), mesh_v(:,:,:), meshuu(:,:), meshvv(:,:)
  real(dp), allocatable :: cont_vel(:,:)
  real(dp), allocatable :: dif(:), dif1(:), rcontorno(:,:)
  integer, allocatable :: contorno(:,:), nodos_cont(:)
  character(len=256) :: ncfile_mom, dimname, line,xlon,xlon2,ylat,ttime,&
         sigmalev,uwind,vwind,z_col,z_col_deg,cont_file,result_format, &
         integer_format

  real(dp) :: dummy

  namelist /corrente/ ncfile_mom,xlon,xlon2,ylat,ttime,sigmalev,uwind,vwind
  open(unit=10, file="namelist_pre.nml", status="old")
  read(10, nml=corrente)
  close(10)
  z_col_deg='../resultados/zz_coluna_deg.dat'
  z_col='../resultados/zz_coluna.dat'
  cont_file='../resultados/c_contorno.cli'

  ! === Leitura do arquivo NetCDF ===
  print*, "Lendo arquivo NetCDF: ", trim(ncfile_mom)
  
  call check(nf90_open(ncfile_mom, NF90_NOWRITE, ncid))
  
  ! Obter dimensões
  call check(nf90_inq_dimid(ncid, ylat, lat_varid))
  call check(nf90_inquire_dimension(ncid, lat_varid, dimname, nlat))
  call check(nf90_inq_dimid(ncid, xlon, lon_varid))
  call check(nf90_inquire_dimension(ncid, lon_varid, dimname, nlon))
  call check(nf90_inq_dimid(ncid, sigmalev, depth_varid))
  call check(nf90_inquire_dimension(ncid, depth_varid, dimname, ndepth))
  call check(nf90_inq_dimid(ncid, ttime, time_varid))
  call check(nf90_inquire_dimension(ncid, time_varid, dimname, ntime))
  call check(nf90_inq_dimid(ncid, xlon2, lon_varid))
  call check(nf90_inquire_dimension(ncid, lon_varid, dimname, nlon2))

  ! Alocar e ler variáveis
  allocate(lat(nlat), lon(nlon), lon2(nlon2),depth(ndepth), time(ntime))
  allocate(uo(nlon2, nlat, ndepth, ntime), vo(nlon, nlat, ndepth, ntime))

  call check(nf90_inq_varid(ncid, ylat, lat_varid))
  call check(nf90_get_var(ncid, lat_varid, lat))
  call check(nf90_inq_varid(ncid, xlon2, lon_varid))
  call check(nf90_get_var(ncid, lon_varid, lon2))
  call check(nf90_inq_varid(ncid, xlon, lon_varid))
  call check(nf90_get_var(ncid, lon_varid, lon))

  call check(nf90_inq_varid(ncid, sigmalev, depth_varid))
  call check(nf90_get_var(ncid, depth_varid, depth))
  call check(nf90_inq_varid(ncid, ttime, time_varid))
  call check(nf90_get_var(ncid, time_varid, time))
  call check(nf90_inq_varid(ncid, uwind, u_varid))
  call check(nf90_get_var(ncid, u_varid, uo))
  call check(nf90_inq_varid(ncid, vwind, v_varid))
  call check(nf90_get_var(ncid, v_varid, vo))
  
  call check(nf90_close(ncid))

  ! Substituir NaN por 1000
  do i = 1, nlon
    do j = 1, nlat
      do k = 1, ndepth
        do L = 1, ntime
          if (isnan(uo(i,j,k,L)) .or.uo(i,j,k,L) .gt. 1000.0 ) uo(i,j,k,L) = 1000.0_dp
        end do
      end do
    end do
  end do

  do i = 1, nlon
    do j = 1, nlat
      do k = 1, ndepth
        do L = 1, ntime
          if (isnan(vo(i,j,k,L)) .or.vo(i,j,k,L) .gt. 1000.0 ) vo(i,j,k,L) = 1000.0_dp
        end do
      end do
    end do
  end do

  ! === Leitura da malha ===
    print*, "Lendo arquivo:",trim(z_col_deg)
  
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
  print*, "Número de colunas em", trim(z_col_deg), ":",ncols
  ! Agora contar número de linhas
  nnos = 0
  do
    read(10, *, iostat=status)
    if (status /= 0) exit
    nnos = nnos + 1
  end do
  rewind(10)
  print*, "Número de linhas em",trim(z_col_deg),":", nnos
  ! Ler dados da malha
  allocate(mesh(nnos, ncols))

  do i = 1, nnos
    read(10, *) mesh(i,:)
  end do
  close(10)

  print*, "Interpolação horizontal"
  ! === Interpolação horizontal ===
  allocate(meshu1(nnos, 2))
  allocate(dif(nlat), dif1(nlon))
  
  do n = 1, nnos
    ! Encontrar ponto mais próximo em lat
    dif = abs(mesh(n,2) - lat)
    meshu1(n,1) = minloc(dif, dim=1)
    
    ! Encontrar ponto mais próximo em lon
    dif1 = abs(mesh(n,1) - lon2)
    meshu1(n,2) = minloc(dif1, dim=1)
  end do
  deallocate(dif, dif1)
  deallocate(lat, lon, time)

  ! Contar número de colunas
  open(unit=19, file=trim(z_col), status='old', action='read')
  ! Ler dados da malha
  allocate(mesh_utm(nnos, ncols))
  do i = 1, nnos
    read(19, *) mesh_utm(i,:)
  end do
  close(19)
  mesh(:,1:2)=mesh_utm(:,1:2)
  deallocate(mesh_utm)

  ! === Interpolação vertical ===
  print*, "Interpolação vertical"
  allocate(mesh_u(nnos, profu,ntime), mesh_v(nnos, profu,ntime))

  do L = 1,ntime		! Para cada passo de tempo
  do m = profu,1,-1	! Para cada nível da malha TELEMAC
    do t = 1, nnos		! Para cada ponto da malha
    ! Encontrar nível mais próximo no MOM
      allocate(dif(ndepth))
      dif = abs(mesh(t,m) - depth)
      k = minloc(dif, dim=1)
      deallocate(dif)
      
      ! Obter valores interpolados
      if (uo(int(meshu1(t,2)), int(meshu1(t,1)), k, L) == 1000.0_dp) then
        mesh_u(t,m,L) = mesh_u(t-1,m,L)
      else
        mesh_u(t,m,L) = uo(int(meshu1(t,2)), int(meshu1(t,1)), k, L)
      end if
      
      if (vo(int(meshu1(t,2)), int(meshu1(t,1)), k, L) == 1000.0_dp) then
        mesh_v(t,m,L) = mesh_v(t-1,m,L)
      else
        mesh_v(t,m,L) = vo(int(meshu1(t,2)), int(meshu1(t,1)), k, L)
      end if
    end do
  end do
  end do
  deallocate(meshu1)

  ! === Tratamento de valores NaN/1000 ===
  ! Para velocidade u
  do L = 1,ntime
  do i = size(mesh_u,2,1), 1, -1
    do j = 1, size(mesh_u,1,1)
      if (mesh_u(j,i,L) == 1000.0_dp) then
        if (j > 1) then
          mesh_u(j,i,L) = mesh_u(j-1,i,L)
        else
          mesh_u(j,i,L) = 0.0_dp
        end if
      end if
    end do
  end do
  end do
  
  ! Para velocidade v
  do L = 1,ntime
  do i = size(mesh_v,2,1), 1, -1
    do j = 1, size(mesh_v,1,1)
      if (mesh_v(j,i,L) == 1000.0_dp) then
        if (j > 1) then
          mesh_v(j,i,L) = mesh_v(j-1,i,L)
        else
          mesh_v(j,i,L) = 0.0_dp
        end if
      end if
    end do
  end do
  end do

  ! === Condições iniciais ===
  print*, "Salvando condições iniciais"
  allocate(meshuu(nnos*profu, 4), meshvv(nnos*profu, 4))
  
  do u = 1, profu
    do j = 1, nnos
      meshuu((u-1)*nnos + j, 1) = real(j, dp)
      meshuu((u-1)*nnos + j, 2:3) = mesh(j,1:2)
      meshuu((u-1)*nnos + j, 4) = mesh_u(j,u,1)
      
      meshvv((u-1)*nnos + j, 1) = real(j, dp)
      meshvv((u-1)*nnos + j, 2:3) = mesh(j,1:2)
      meshvv((u-1)*nnos + j, 4) = mesh_v(j,u,1)
    end do
  end do

  ! Salvar condições iniciais
  open(unit=20, file='../resultados/Vel_u_ini.dat', status='replace')
  do i = 1, size(meshuu,1)
    write(20, '(F15.1,2F15.6,F15.10)') meshuu(i,:)
  end do
  close(20)
  
  open(unit=21, file='../resultados/Vel_v_ini.dat', status='replace')
  do i = 1, size(meshvv,1)
    write(21, '(F15.1,2F15.6,F15.10)') meshvv(i,:)
  end do
  close(21)
  deallocate(meshuu,meshvv)
  ! === Leitura do contorno ===
  print*, "Lendo arquivo de contorno:",trim(cont_file)
  open(unit=11, file=trim(cont_file), status='old', action='read')

  ! Contar número de linhas
  ncontorno = 0
  do
    read(11, *, iostat=status)
    if (status /= 0) exit
    ncontorno = ncontorno + 1
  end do
  rewind(11)
  print*, "Número de linhas em", trim(cont_file), ":",ncontorno

  ! Ler dados do contorno
  allocate(contorno(ncontorno, 13), rcontorno(ncontorno, 13))
  do i = 1, ncontorno
    read(11, *) rcontorno(i,:)
    contorno(i,:) = int(rcontorno(i,:))
  end do
  close(11)

  ! === Processamento dos contornos ===
  ! Identificar nós de contorno tipo 4
  print*,'Salvando condições de contorno'
  g = 0
  do i = 1, ncontorno
    if (contorno(i,2) == 6 .or. contorno(i,2) == 4 .or. contorno(i,2) == 5) g = g + 1
  end do
  allocate(nodos_cont(g))
  
  g = 0
  do i = 1, ncontorno
    if (contorno(i,2) == 6 .or. contorno(i,2) == 4 .or. contorno(i,2) == 5) then
      g = g + 1
      nodos_cont(g) = contorno(i,12)
    end if
  end do

  ! === Preparar arquivo de saída ===
  allocate(cont_vel(2*ntime+1, g*profu+1))
  
  ! Cabeçalho
  cont_vel(1,1) = real(g, dp)
  do i = 1, g
    cont_vel(1, (i-1)*profu+2 : i*profu+1) = real(nodos_cont(i), dp)
    cont_vel(1, (i-1)*profu+3 : i*profu+1) = 0.0_dp
  end do

  ! Dados temporais
  cont = 1
  contt = 2
  do L = 1, ntime
    print*, "Processando tempo: ", tempo, "s (", L, "/", ntime, ")"
    
    do i = 1, g !nivel
      d = nodos_cont(i)
      cont_vel(2*L-1+cont, (i-1)*profu+2 : i*profu+1) = mesh_u(d,:,L)
      cont_vel(2*L-1+contt, (i-1)*profu+2 : i*profu+1) = mesh_v(d,:,L)
    end do
    
    cont_vel(2*L-1+cont,1) = tempo
    cont_vel(2*L-1+contt,1) = tempo
    tempo = tempo + 3600.0_dp  ! Incremento de 1 hora
  end do

  ! Salvar arquivo final
  open(unit=22, file='../resultados/cont_velo.dat', status='replace')
  cont=1
  contt=2

  write(integer_format, "(I0)") size(cont_vel(1,:))
  result_format = "("//trim(integer_format)//"F16.8)"

  write(22, trim(result_format)) cont_vel(1, : )
  do L = 1,ntime
    write(22, trim(result_format)) cont_vel(2*L-1+cont, : ) 
    write(22, trim(result_format)) cont_vel(2*L-1+contt, : )
  enddo
  close(22)

  print*, "Arquivos gerados com sucesso!"
  print*, "- Vel_u_ini.dat (condições iniciais de velocidade u)"
  print*, "- Vel_v_ini.dat (condições iniciais de velocidade v)"
  print*, "- cont_velo.dat (contornos de velocidade)"

contains
  subroutine check(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
      print*, "Erro NetCDF: ", trim(nf90_strerror(status))
      stop 1
    end if
  end subroutine

  ! Função para verificar NaN (simplificada)
  elemental logical function isnan(x)
    real(dp), intent(in) :: x
    isnan = (x /= x)
  end function isnan
end program interp_velocidades
