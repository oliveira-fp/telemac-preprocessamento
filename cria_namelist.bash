#!/bin/bash

#eta_file='../dados/Eta03_BESM_2007010100_2D.nc'
eta_file=`ls ../dados/*Eta*`
dt_eta=3600
lon_eta='lon'
lat_eta='lat'
time_eta='time'
u10_eta='u10m'
v10_eta='v10m'

#besm_file='../dados/20070101.ocean_telemac_recorte.nc'
besm_file=`ls ../dados/*telemac*`
lon_besm='xh'
lat_besm='yh'
time_besm='time'
sigma_besm='sigma2_l'
ssh_besm='ssh'
temp_besm='thetao'
salt_besm='so'
u_besm='uo'
v_besm='vo'

echo "
&atmosferico
       ncfile = $eta_file,
       dt = $dt_eta,
       xlon = $lon_eta,
       ylat = $lat_eta,
       ttime = $time_eta,
       uwind = $u10_eta,
       vwind = $v10_eta
 /

 &ssh_config
       ncfile = $besm_file,
       xlon = $lon_besm,
       ylat = $lat_besm,
       ttime = $time_besm,
       hsea = $ssh_besm,
 /

 &temp_salt
       ncfile = $besm_file,
       xlon = $lon_besm,
       ylat = $lat_besm,
       ttime = $time_besm,
       sigmalev = $sigma_besm,
       ttemp = $temp_besm,
       salt = $salt_besm,
 /

 &corrente
       ncfile = $besm_file,
       xlon = $lon_besm,
       ylat = $lat_besm,
       ttime = $time_besm,
       sigmalev = $sigma_besm,
       uwind = $u_besm,
       vwind = $v_besm
 /
" > namelist_pre.nml

