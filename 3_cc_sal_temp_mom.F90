program interp_sal_temp
  use netcdf
  implicit none

  ! Parâmetros
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer :: profu = 10  ! Número de níveis da malha
  integer :: ncid, lat_varid, lon_varid, sal_varid, temp_varid, depth_varid, time_varid
  integer :: nlat, nlon, ndepth, ntime, nnos, ncontorno
  integer :: i, j, k, L, m, n, t, u, cont, contt, g, d, dd, status, ncols
  real(dp) :: tempo = 0.0_dp
  real(dp), allocatable :: lat(:), lon(:), depth(:), time(:)
  real(dp), allocatable :: sal(:,:,:,:), temp(:,:,:,:)
  real(dp), allocatable :: mesh(:,:), mat1(:,:), meshu1(:,:), local(:,:)
  real(dp), allocatable :: meshsal(:,:), meshtem(:,:), meshuu(:,:), meshvv(:,:)
  real(dp), allocatable :: contvelu(:,:), contvelv(:,:), cont_vel(:,:)
  real(dp), allocatable :: dif(:), dif1(:), dif2(:),rcontorno(:,:)
  integer, allocatable :: contorno(:,:), nodos_cont(:)
  character(len=256) :: ncfile, dimname, line,xlon,ylat,ttime,sigmalev,ttemp,salt
  logical, allocatable :: mask(:,:,:,:)
  real(dp) :: dummy

  namelist /temp_salt/ ncfile,xlon,ylat,ttime,sigmalev,ttemp,salt
  open(unit=10, file="namelist_pre.nml", status="old")
  read(10, nml=temp_salt)
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
  allocate(sal(nlon, nlat, ndepth, ntime), temp(nlon, nlat, ndepth, ntime))
  allocate(mask(nlon, nlat, ndepth, ntime))
  
  call check(nf90_inq_varid(ncid, ylat, lat_varid))
  call check(nf90_get_var(ncid, lat_varid, lat))
  call check(nf90_inq_varid(ncid, xlon, lon_varid))
  call check(nf90_get_var(ncid, lon_varid, lon))
  call check(nf90_inq_varid(ncid, sigmalev, depth_varid))
  call check(nf90_get_var(ncid, depth_varid, depth))
  call check(nf90_inq_varid(ncid, ttime, time_varid))
  call check(nf90_get_var(ncid, time_varid, time))
  call check(nf90_inq_varid(ncid, salt, sal_varid))
  call check(nf90_get_var(ncid, sal_varid, sal))
  call check(nf90_inq_varid(ncid, ttemp, temp_varid))
  call check(nf90_get_var(ncid, temp_varid, temp))
  
  call check(nf90_close(ncid))

  ! Substituir NaN por 1000 - método alternativo sem WHERE
  mask = isnan(sal)
  do i = 1, nlon
    do j = 1, nlat
      do k = 1, ndepth
        do L = 1, ntime
          if (mask(i,j,k,L)) sal(i,j,k,L) = 1000.0_dp
        end do
      end do
    end do
  end do

  mask = isnan(temp)
  do i = 1, nlon
    do j = 1, nlat
      do k = 1, ndepth
        do L = 1, ntime
          if (mask(i,j,k,L)) temp(i,j,k,L) = 1000.0_dp
        end do
      end do
    end do
  end do

  ! === Leitura da malha ===
    print*, "Lendo arquivo: zz_coluna.dat"
  
  ! Primeiro contar o número de colunas
  open(unit=10, file='../resultados/zz_coluna_deg.dat', status='old', action='read')
  read(10, '(a)', iostat=status) line ! Ler primeira linha para contar colunas
  ncols = 0
  do
    read(line, *, iostat=status) (dummy, i=1,ncols+1) 
    if (status /= 0) exit
    ncols = ncols + 1
  end do
  rewind(10)
  print*, "Número de colunas em zz_coluna_deg.dat: ",ncols
  ! Agora contar número de linhas
  nnos = 0

  do
    read(10, *, iostat=status)
    if (status /= 0) exit
    nnos = nnos + 1
  end do
  rewind(10)
  print*, "Número de linhas em zz_coluna_deg.dat: ", nnos
  ! Ler dados da malha
  print*, "Carregando malha"
  allocate(mesh(nnos, ncols))
  do i = 1, nnos
    read(10, *) mesh(i,:)
  end do
  close(10)

  ! Converter UTM para lat/lon (simplificado - usar biblioteca apropriada)
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
  print*, "Número de linhas em c_contorno_BC.cli: ",ncontorno
 
  ! Ler dados do contorno
  allocate(contorno(ncontorno, 13),  rcontorno(ncontorno, 13))
  do i = 1, ncontorno
    read(11, *) rcontorno(i,:)
    contorno(i,:)=int(rcontorno(i,:))
!    read(rcontorno(i,:),*) contorno(i,:)
  end do
  close(11)
  print*, "Interpolação horizontal"
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
  local(:,3:4) = real(meshu1(:,3:4), dp)  ! Converter para real
  local(:,5:) = mesh(:,3:)

  ! === Interpolação vertical ===
  print*, "Interpolação vertical"
  allocate(meshsal(nnos, profu+2), meshtem(nnos, profu+2))
  allocate(dif2(ndepth))
  
  do m = 5, size(local,2)  ! Para cada nível da malha TELEMAC
    do t = 1, nnos         ! Para cada ponto da malha
      ! Encontrar nível mais próximo no MOM
      dif2 = abs(local(t,m) - depth)
      k = minloc(dif2, dim=1)
      
      ! Obter valores interpolados
      meshsal(t,1:2) = mesh(t,1:2)
      meshtem(t,1:2) = mesh(t,1:2)
      
      if (sal(int(local(t,4)), int(local(t,3)), k, 1) == 1000.0_dp) then
        meshsal(t,m-2) = 0.0_dp  ! Valor padrão para NaN
      else
        meshsal(t,m-2) = sal(int(local(t,4)), int(local(t,3)), k, 1)
      end if
      
      if (temp(int(local(t,4)), int(local(t,3)), k, 1) == 1000.0_dp) then
        meshtem(t,m-2) = 20.0_dp  ! Valor padrão para NaN
      else
        meshtem(t,m-2) = temp(int(local(t,4)), int(local(t,3)), k, 1)
      end if
    end do
  end do

  ! === Tratamento de valores NaN/1000 ===
  ! Para salinidade
  do i = 3, size(meshsal,2)
    do j = 1, size(meshsal,1)
      if (meshsal(j,i) == 1000.0_dp) meshsal(j,i) = 0.0_dp
    end do
  end do
  
  do i = size(meshsal,2)-1, 3, -1
    do j = 1, size(meshsal,1)
      if (meshsal(j,i) == 1000.0_dp) meshsal(j,i) = meshsal(j,i+1)
    end do
  end do

  ! Para temperatura
  do i = 3, size(meshtem,2)
    do j = 1, size(meshtem,1)
      if (meshtem(j,i) == 1000.0_dp) meshtem(j,i) = 20.0_dp
    end do
  end do
  
  do i = size(meshtem,2)-1, 3, -1
    do j = 1, size(meshtem,1)
      if (meshtem(j,i) == 1000.0_dp) meshtem(j,i) = meshtem(j,i+1)
    end do
  end do

  ! === Condições iniciais ===
  allocate(meshuu(nnos*profu, 4), meshvv(nnos*profu, 4))
  
  do u = 1, profu
    do j = 1, nnos
      meshuu((u-1)*nnos + j, 1) = real(j, dp)
      meshuu((u-1)*nnos + j, 2:3) = meshsal(j,1:2)
      meshuu((u-1)*nnos + j, 4) = meshsal(j,u+2)
      
      meshvv((u-1)*nnos + j, 1) = real(j, dp)
      meshvv((u-1)*nnos + j, 2:3) = meshtem(j,1:2)
      meshvv((u-1)*nnos + j, 4) = meshtem(j,u+2)
    end do
  end do

  ! Salvar condições iniciais
  open(unit=20, file='../resultados/Sal_ini.dat', status='replace')
  do i = 1, size(meshuu,1)
    write(20, '(4F15.6)') meshuu(i,:)
  end do
  close(20)
  
  open(unit=21, file='../resultados/Tem_ini.dat', status='replace')
  do i = 1, size(meshvv,1)
    write(21, '(4F15.6)') meshvv(i,:)
  end do
  close(21)

  ! === Processamento dos contornos ===
  ! Identificar nós de contorno
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

  allocate(contvelu(g, profu+3), contvelv(g, profu+3))
  
  do i = 1, g
    d = nodos_cont(i)
    contvelu(i,1:3) = [real(d, dp), mesh(d,1:2)]
    contvelv(i,1:3) = [real(d, dp), mesh(d,1:2)]
    
    if (any(contorno(:,12) == d .and. contorno(:,2) == 5)) then
      ! Contorno tipo 5 (rios)
      contvelu(i,4:) = 0.0_dp    ! Salinidade zero
      contvelv(i,4:) = 20.0_dp    ! Temperatura 20°C
    else
      ! Outros contornos
      contvelu(i,4:) = meshsal(d,3:)
      contvelv(i,4:) = meshtem(d,3:)
    end if
  end do

  ! === Preparar arquivo de saída ===
  allocate(cont_vel(2*ntime+1, g*profu+1))
  
  ! Cabeçalho
  cont_vel(1,1) = real(g, dp)
  do i = 1, g
    cont_vel(1, (i-1)*profu+2 : i*profu+1) = real(nodos_cont(i), dp)
    cont_vel(1, (i-1)*profu+3 : i*profu+1) = 0.0_dp  ! Zeros adicionais
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
  open(unit=22, file='../resultados/cont_sal_temp.dat', status='replace')
  do i = 1, size(cont_vel,1)
    write(22, '(10000F15.6)') cont_vel(i,:)  ! Ajuste o formato conforme necessário
  end do
  close(22)

  print*, "Arquivos gerados com sucesso!"
  print*, "- Sal_ini.dat (condições iniciais de salinidade)"
  print*, "- Tem_ini.dat (condições iniciais de temperatura)"
  print*, "- cont_sal_temp.dat (contornos de salinidade e temperatura)"

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
end program interp_sal_temp
