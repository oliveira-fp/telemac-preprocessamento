program interp_velocidades
  use netcdf
  implicit none

  ! Parâmetros
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer :: profu = 10  ! Número de níveis da malha
  integer :: ncid, lat_varid, lon_varid, u_varid, v_varid, depth_varid, time_varid
  integer :: nlat, nlon, ndepth, ntime, nnos, ncontorno
  integer :: i, j, k, L, m, n, t, u, cont, contt, g, d, status, ncols
  real(dp) :: tempo = 0.0_dp
  real(dp), allocatable :: lat(:), lon(:), depth(:), time(:)
  real(dp), allocatable :: uo(:,:,:,:), vo(:,:,:,:)
  real(dp), allocatable :: mesh(:,:), mat1(:,:), meshu1(:,:), local(:,:)
  real(dp), allocatable :: mesh_u(:,:), mesh_v(:,:), meshuu(:,:), meshvv(:,:)
  real(dp), allocatable :: contvelu(:,:), contvelv(:,:), cont_vel(:,:)
  real(dp), allocatable :: dif(:), dif1(:), dif2(:), rcontorno(:,:)
  integer, allocatable :: contorno(:,:), nodos_cont(:)
  character(len=256) :: ncfile, dimname, line,xlon,ylat,ttime,sigmalev,uwind,vwind
  logical, allocatable :: mask(:,:,:,:)
  real(dp) :: dummy

  namelist /corrente/ ncfile,xlon,ylat,ttime,sigmalev,uwind,vwind
  open(unit=10, file="namelist_pre.nml", status="old")
  read(10, nml=corrente)
  close(10)

  ! === Leitura do arquivo NetCDF ===
!  ncfile = '../dados/20070101.ocean_telemac_recorte.nc'
  print*, "Lendo arquivo NetCDF: ", trim(ncfile)
  
  call check(nf90_open(ncfile, NF90_NOWRITE, ncid))
  
  ! Obter dimensões
  call check(nf90_inq_dimid(ncid, ylat, lat_varid))
  call check(nf90_inquire_dimension(ncid, lat_varid, dimname, nlat))
  call check(nf90_inq_dimid(ncid, xlon, lon_varid))
  call check(nf90_inquire_dimension(ncid, lon_varid, dimname, nlon))
  call check(nf90_inq_dimid(ncid, sigmalev, depth_varid))
  call check(nf90_inquire_dimension(ncid, depth_varid, dimname, ndepth))
  call check(nf90_inq_dimid(ncid, ttime, time_varid))
  call check(nf90_inquire_dimension(ncid, time_varid, dimname, ntime))

  ! Alocar e ler variáveis
  allocate(lat(nlat), lon(nlon), depth(ndepth), time(ntime))
  allocate(uo(nlon, nlat, ndepth, ntime), vo(nlon, nlat, ndepth, ntime))
  allocate(mask(nlon, nlat, ndepth, ntime))

  call check(nf90_inq_varid(ncid, ylat, lat_varid))
  call check(nf90_get_var(ncid, lat_varid, lat))
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
  mask = isnan(uo)
  do i = 1, nlon
    do j = 1, nlat
      do k = 1, ndepth
        do L = 1, ntime
          if (mask(i,j,k,L)) uo(i,j,k,L) = 1000.0_dp
        end do
      end do
    end do
  end do

  mask = isnan(vo)
  do i = 1, nlon
    do j = 1, nlat
      do k = 1, ndepth
        do L = 1, ntime
          if (mask(i,j,k,L)) vo(i,j,k,L) = 1000.0_dp
        end do
      end do
    end do
  end do

  ! === Leitura da malha ===
  print*, "Lendo arquivo da malha: zz_coluna.dat"
  
  ! Contar número de colunas
  open(unit=10, file='../resultados/zz_coluna_deg.dat', status='old', action='read')
  read(10, '(a)', iostat=status) line
  ncols = 0
  do
    read(line, *, iostat=status) (dummy, i=1,ncols+1)
    if (status /= 0) exit
    ncols = ncols + 1
  end do
  rewind(10)
  
  ! Contar número de linhas
  nnos = 0
  do
    read(10, *, iostat=status)
    if (status /= 0) exit
    nnos = nnos + 1
  end do
  rewind(10)
  
  ! Ler dados da malha
  allocate(mesh(nnos, ncols))
  do i = 1, nnos
    read(10, *) mesh(i,:)
  end do
  close(10)

  allocate(mat1(nnos, 3))
  do n = 1, nnos
    mat1(n,1) = mesh(n,1) ! longitude placeholder
    mat1(n,2) = mesh(n,2) ! latitude placeholder
    mat1(n,3) = real(n, dp)
  end do

  ! === Leitura do contorno ===
  print*, "Lendo arquivo de contorno: c_contorno_BC.cli"
  open(unit=11, file='../dados/c_contorno_BC.cli', status='old', action='read')
  
  ! Contar número de linhas
  ncontorno = 0
  do
    read(11, *, iostat=status)
    if (status /= 0) exit
    ncontorno = ncontorno + 1
  end do
  rewind(11)
  
  ! Ler dados do contorno
  allocate(contorno(ncontorno, 13), rcontorno(ncontorno, 13))
  do i = 1, ncontorno
    read(11, *) rcontorno(i,:)
    contorno(i,:) = int(rcontorno(i,:))
  end do
  close(11)

  ! === Interpolação horizontal ===
  allocate(meshu1(nnos, 4), local(nnos, size(mesh,2)+2))
  allocate(dif(nlat), dif1(nlon))
  
  do n = 1, nnos
    ! Encontrar ponto mais próximo em lat
    dif = abs(mat1(n,2) - lat)
    meshu1(n,3) = minloc(dif, dim=1)
    
    ! Encontrar ponto mais próximo em lon
    dif1 = abs(mat1(n,1) - lon)
    meshu1(n,4) = minloc(dif1, dim=1)
    
    meshu1(n,1:2) = mesh(n,1:2)
  end do
  
  local(:,1:2) = mesh(:,1:2)
  local(:,3:4) = real(meshu1(:,3:4), dp)
  local(:,5:) = mesh(:,3:)

  ! === Interpolação vertical ===
  allocate(mesh_u(nnos, profu+2), mesh_v(nnos, profu+2))
  allocate(dif2(ndepth))
  
  do m = 5, size(local,2)
    do t = 1, nnos
      ! Encontrar nível mais próximo no MOM
      dif2 = abs(local(t,m) - depth)
      k = minloc(dif2, dim=1)
      
      ! Obter valores interpolados
      mesh_u(t,1:2) = mesh(t,1:2)
      mesh_v(t,1:2) = mesh(t,1:2)
      
      if (uo(int(local(t,4)), int(local(t,3)), k, 1) == 1000.0_dp) then
        mesh_u(t,m-2) = 0.0_dp
      else
        mesh_u(t,m-2) = uo(int(local(t,4)), int(local(t,3)), k, 1)
      end if
      
      if (vo(int(local(t,4)), int(local(t,3)), k, 1) == 1000.0_dp) then
        mesh_v(t,m-2) = 0.0_dp
      else
        mesh_v(t,m-2) = vo(int(local(t,4)), int(local(t,3)), k, 1)
      end if
    end do
  end do

  ! === Tratamento de valores NaN/1000 ===
  ! Para velocidade u
  do i = size(mesh_u,2), 3, -1
    do j = 1, size(mesh_u,1)
      if (mesh_u(j,i) == 1000.0_dp) then
        if (j > 1) then
          mesh_u(j,i) = mesh_u(j-1,i)
        else
          mesh_u(j,i) = 0.0_dp
        end if
      end if
    end do
  end do
  
  ! Para velocidade v
  do i = size(mesh_v,2), 3, -1
    do j = 1, size(mesh_v,1)
      if (mesh_v(j,i) == 1000.0_dp) then
        if (j > 1) then
          mesh_v(j,i) = mesh_v(j-1,i)
        else
          mesh_v(j,i) = 0.0_dp
        end if
      end if
    end do
  end do

  ! === Condições iniciais ===
  allocate(meshuu(nnos*profu, 4), meshvv(nnos*profu, 4))
  
  do u = 1, profu
    do j = 1, nnos
      meshuu((u-1)*nnos + j, 1) = real(j, dp)
      meshuu((u-1)*nnos + j, 2:3) = mesh_u(j,1:2)
      meshuu((u-1)*nnos + j, 4) = mesh_u(j,u+2)
      
      meshvv((u-1)*nnos + j, 1) = real(j, dp)
      meshvv((u-1)*nnos + j, 2:3) = mesh_v(j,1:2)
      meshvv((u-1)*nnos + j, 4) = mesh_v(j,u+2)
    end do
  end do

  ! Salvar condições iniciais
  open(unit=20, file='../resultados/Vel_u_ini.dat', status='replace')
  do i = 1, size(meshuu,1)
    write(20, '(4F15.6)') meshuu(i,:)
  end do
  close(20)
  
  open(unit=21, file='../resultados/Vel_v_ini.dat', status='replace')
  do i = 1, size(meshvv,1)
    write(21, '(4F15.6)') meshvv(i,:)
  end do
  close(21)

  ! === Processamento dos contornos ===
  ! Identificar nós de contorno tipo 4
  g = count(contorno(:,2) == 4)
  allocate(nodos_cont(g))
  
  g = 0
  do i = 1, ncontorno
    if (contorno(i,2) == 4) then
      g = g + 1
      nodos_cont(g) = contorno(i,12)
    end if
  end do

  allocate(contvelu(g, profu+3), contvelv(g, profu+3))
  
  do i = 1, g
    d = nodos_cont(i)
    contvelu(i,1:3) = [real(d, dp), mesh(d,1:2)]
    contvelv(i,1:3) = [real(d, dp), mesh(d,1:2)]
    contvelu(i,4:) = mesh_u(d,3:)
    contvelv(i,4:) = mesh_v(d,3:)
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
    
    do i = 1, g
      cont_vel(L+cont, (i-1)*profu+2 : i*profu+1) = contvelu(i,4:)
      cont_vel(L+contt, (i-1)*profu+2 : i*profu+1) = contvelv(i,4:)
    end do
    
    cont_vel(L+cont,1) = tempo
    cont_vel(L+contt,1) = tempo
    
    cont = cont + 1
    contt = contt + 1
    tempo = tempo + 3600.0_dp  ! Incremento de 1 hora
  end do

  ! Salvar arquivo final
  open(unit=22, file='../resultados/cont_velo.dat', status='replace')
  do i = 1, size(cont_vel,1)
    write(22, '(10000F15.6)') cont_vel(i,:)
  end do
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

  elemental logical function isnan(x)
    real(dp), intent(in) :: x
    isnan = (x /= x)
  end function isnan
end program interp_velocidades
