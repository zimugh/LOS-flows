pro g2s

datadir='/home/test/ltdata/'
savdir='../'

itstep=0
start=160
endd=161
for ii=start,endd-1,1 do begin
itstep=ii

print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print,'itstep',itstep
print,SYSTIME()
print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

;if (itstep mod 10) eq 0 then begin
;continue
;endif 

getdata_sph,datadir,itstep,b1c,b2c,b3c,v1c,v2c,v3c,dc,tempc,pc,x1c,x2c,x3c,time


dim=size(b1c)
print,dim
;dimx = dim(1)
;dimy = dim(2)
;dimz = dim(3)

;vr=fltarr(dimx,dimy,dimz)
;vth=fltarr(dimx,dimy,dimz)
;vph=fltarr(dimx,dimy,dimz)

;br=fltarr(dimx,dimy,dimz)
;bth=fltarr(dimx,dimy,dimz)
;vph=fltarr(dimx,dimy,dimz)

;dens=dblarr(dimx,dimy,dimz)
;pres=dblarr(dimx,dimy,dimz)
;temp=fltarr(dimx,dimy,dimz)

;r=fltarr(dimx)
;th=fltarr(dimy)
;ph=fltarr(dimz)


;br=b1c
;bth=b2c
;bph=b3c
;dens=dc*double(1.0)
;pres=pc*double(1.0)
;temp=tempc
;r=x1c
;th=x2c
;ph=x3c
;vr=v1c
;vth=v2c
;vph=v3c

;struct0={cube,br:b1c,bth:b2c,bph:b3c, dens:dc*double(1.0), pres:pc*double(1.0), temp:tempc, $
;r:x1c,th:x2c,ph:x3c, vr:v1c,vth:v2c,vph:v3c}
itstep0=160
savdata='../P3DS/mfrdata/P3D/data4sav/0simuresults.sav'
restore,savdata
savrho=fs_rho_t[*,*,*,itstep0]
savtem=fs_tem_t[*,*,*,itstep0]
savvel=fs_v1_t[*,*,*,itstep0]
savpre=fs_p_t[*,*,*,itstep0]

unit_density =   2.341670487800000E-015
unit_pressure =   0.317549224000000     
unit_velocity =    11645084.2956225  
unitte=1e6

savrho=savrho*unit_density/1.0e9
savtem=savtem
savvel=savvel*unit_velocity
savpre=savpre*unit_pressure

b1c=b1c[0:169:2,17:178:2,278:681:2]
b2c=b2c[0:169:2,17:178:2,278:681:2]
b3c=b3c[0:169:2,17:178:2,278:681:2]

x1c=x1c[0:169:2]
x2c=x2c[17:178:2]
x3c=x3c[278:681:2]

dim=size(b1c)
print,dim

dc=TRANSPOSE(savrho, [2, 0, 1])
pc=TRANSPOSE(savpre, [2, 0, 1])
tempc=TRANSPOSE(savtem, [2, 0, 1])
;v1c=fltarr(81, 202, 85)
;v2c=fltarr(81, 202, 85)

v1c=fltarr(85, 81, 202)
v2c=fltarr(85, 81, 202)

v3c=TRANSPOSE(savvel, [2, 0, 1])

dim=size(v3c)
print,dim



cube={br:b1c,bth:b2c,bph:b3c, dens:dc*double(1e24), pres:pc*double(1.0), temp:tempc, $
r:x1c,th:x2c,ph:x3c, vr:v1c/1.0e5,vth:v2c/1.0e5,vph:v3c/1.0e5,global:0,axisym:0,hydro:0,title:'idlsave_'+strtrim(string(itstep0,format='(i4.4)'),1)+'_format',filename:'./'+'for_idlsave_'+strtrim(string(itstep0,format='(i4.4)'),1)+'_format.dat'}


;save, br,bth,bph,vr,vth,vph,dens,temp,pres,r,th,ph,$
save, cube,$
filename=savdir+'nfor_idlsave3_'+strtrim(string(itstep0,format='(i4.4)'),1)+'_format.dat'
print,savdir+'nfor_idlsave3_'+strtrim(string(itstep0,format='(i4.4)'),1)+'_format.dat'
endfor
end


PRO getdata_sph,datadir,itstep,b1c,b2c,b3c,v1c,v2c,v3c,dc,tempc,pc,x1c,x2c,x3c,time
;INPUT:
;datadir: full path of the directory of the data
;itstep: integer time step
;OUTPUt:
;b1c: br array (in G)
;b2c: bth array (in G)
;b3c: bph array (in G)
;dc: density array (in g/cm^3)
;pc: pressue array (in dyn/cm^2)
;tempc: temperature array (in K)
;time: time (in sec)
;v1c: vr (in cm/s)
;v2c: vth (in cm/s)
;v3c: vph (in cm/s)
;x1c: r coordinates (in Rsun)
;x2c: theta coordinates (in rad)
;x3c: phi coordinates (in rad)
;all quantities above are in CGS units

header=datadir

;read grid.dat
in=long(85)
jn=long(85)
kn=long(85)
file=header+'grid.dat'
print,'read ',file
openr,/get_lun,unit,file
readu,unit,in,jn,kn

is=3-1 & js=3-1 & ks=3-1
ie=is+(in-6) & je=js+(jn-6) & ke=ks+(kn-6)
ism1=is-1 & jsm1=js-1 & ksm1=ks-1
ism2=is-2 & jsm2=js-2 & ksm2=ks-2
iep1=ie+1 & jep1=je+1 & kep1=ke+1
iep2=ie+2 & jep2=je+2 & kep2=ke+2
iep3=ie+3 & jep3=je+3 & kep3=ke+3

x1a=dblarr(iep3-ism2+1)
x1b=dblarr(iep2-ism2+1)
x2a=dblarr(jep3-jsm2+1)
x2b=dblarr(jep2-jsm2+1)
x3a=dblarr(kep3-ksm2+1)
x3b=dblarr(kep2-ksm2+1)

dx1a=dblarr(iep2-ism2+1)
dx1b=dblarr(iep3-ism2+1)
dx2a=dblarr(jep2-jsm2+1)
dx2b=dblarr(jep3-jsm2+1)
dx3a=dblarr(kep2-ksm2+1)
dx3b=dblarr(kep3-ksm2+1)

g2a=dblarr(iep3-ism2+1)
g2b=dblarr(iep2-ism2+1)
g31a=dblarr(iep3-ism2+1)
g31b=dblarr(iep2-ism2+1)
g32a=dblarr(jep3-jsm2+1)
g32b=dblarr(jep2-jsm2+1)

dg2bd1=dblarr(iep3-ism2+1)
dg2ad1=dblarr(iep2-ism2+1)
dg31bd1=dblarr(iep3-ism2+1)
dg31ad1=dblarr(iep2-ism2+1)
dg32bd2=dblarr(jep3-ism2+1)
dg32ad2=dblarr(jep2-ism2+1)

dvl1a=dblarr(iep2-ism2+1)
dvl1b=dblarr(iep3-ism2+1)
dvl2a=dblarr(jep2-jsm2+1)
dvl2b=dblarr(jep3-jsm2+1)
dvl3a=dblarr(kep2-ksm2+1)
dvl3b=dblarr(kep3-ksm2+1)

dx1ai=dblarr(iep2-ism2+1)
dx1bi=dblarr(iep3-ism2+1)
dx2ai=dblarr(jep2-jsm2+1)
dx2bi=dblarr(jep3-jsm2+1)
dx3ai=dblarr(kep2-ksm2+1)
dx3bi=dblarr(kep3-ksm2+1)

g2ai=dblarr(iep3-ism2+1)
g2bi=dblarr(iep2-ism2+1)
g31ai=dblarr(iep3-ism2+1)
g31bi=dblarr(iep2-ism2+1)
g32ai=dblarr(jep3-jsm2+1)
g32bi=dblarr(jep2-jsm2+1)

dvl1ai=dblarr(iep2-ism2+1)
dvl1bi=dblarr(iep3-ism2+1)
dvl2ai=dblarr(jep2-jsm2+1)
dvl2bi=dblarr(jep3-jsm2+1)
dvl3ai=dblarr(kep2-ksm2+1)
dvl3bi=dblarr(kep3-ksm2+1)

readu,unit,x1a
readu,unit,x1b
readu,unit,x2a
readu,unit,x2b
readu,unit,x3a
readu,unit,x3b

readu,unit,dx1a
readu,unit,dx1b
readu,unit,dx2a
readu,unit,dx2b
readu,unit,dx3a
readu,unit,dx3b

readu,unit,g2a
readu,unit,g2b
readu,unit,g31a
readu,unit,g31b
readu,unit,g32a
readu,unit,g32b

readu,unit,dg2bd1
readu,unit,dg2ad1
readu,unit,dg31bd1
readu,unit,dg31ad1
readu,unit,dg32bd2
readu,unit,dg32ad2

readu,unit,dvl1a
readu,unit,dvl1b
readu,unit,dvl2a
readu,unit,dvl2b
readu,unit,dvl3a
readu,unit,dvl3b

readu,unit,dx1ai
readu,unit,dx1bi
readu,unit,dx2ai
readu,unit,dx2bi
readu,unit,dx3ai
readu,unit,dx3bi

readu,unit,g2ai
readu,unit,g2bi
readu,unit,g31ai
readu,unit,g31bi
readu,unit,g32ai
readu,unit,g32bi

readu,unit,dvl1ai
readu,unit,dvl1bi
readu,unit,dvl2ai
readu,unit,dvl2bi
readu,unit,dvl3ai
readu,unit,dvl3bi

free_lun,unit

;read physparams.dat
g_const=0.D0
rgas=0.D0
kboltz=0.D0
mproton=0.D0
pi=0.D0
gamma=0.D0
unit_rho=0.D0
unit_len=0.D0
unit_b=0.D0
unit_temp=0.D0
unit_v=0.D0
unit_time=0.D0
msol=0.D0
rsol=0.D0
temp_c=0.D0
rho_c=0.D0
muconst=0.D0

file=header+'physparams.dat'
print,'read ',file
openr,/get_lun,unit,file

readu,unit,g_const,rgas,kboltz,mproton,pi,gamma,$
unit_rho,unit_len,unit_b,unit_temp,unit_v,unit_time,$
msol,rsol,temp_c,rho_c,muconst

free_lun,unit

help,in-5,jn-5,kn-5
help,g_const,rgas,kboltz,mproton,pi,gamma,$
unit_rho,unit_len,unit_b,unit_temp,unit_v,unit_time,$
msol,rsol,temp_c,rho_c,muconst

rstar=kboltz/(muconst*mproton)

tiny=1.E-30

print,'unit_len,unit_rho,unit_v,unit_time,unit_b,muconst,rstar,gamma',$
unit_len,unit_rho,unit_v,unit_time,unit_b,muconst,rstar,gamma

file=header+'b1_'+strtrim(string(itstep,format='(i4.4)'),1)+'.dat'
read_data,data,time,file,in,jn,kn
b1c=0.50*(float(data)+shift(float(data),-1,0,0))
b1c=b1c(is:ie,js:je,ks:ke)

file=header+'v1_'+strtrim(string(itstep,format='(i4.4)'),1)+'.dat'
read_data,data,time,file,in,jn,kn
v1c=0.50*(float(data)+shift(float(data),-1,0,0))
v1c=v1c(is:ie,js:je,ks:ke)

file=header+'b2_'+strtrim(string(itstep,format='(i4.4)'),1)+'.dat'
read_data,data,time,file,in,jn,kn
b2c=0.50*(float(data)+shift(float(data),0,-1,0))
b2c=b2c(is:ie,js:je,ks:ke)

file=header+'v2_'+strtrim(string(itstep,format='(i4.4)'),1)+'.dat'
read_data,data,time,file,in,jn,kn
v2c=0.50*(float(data)+shift(float(data),0,-1,0))
v2c=v2c(is:ie,js:je,ks:ke)

file=header+'b3_'+strtrim(string(itstep,format='(i4.4)'),1)+'.dat'
read_data,data,time,file,in,jn,kn
b3c=0.50*(float(data)+shift(float(data),0,0,-1))
b3c=b3c(is:ie,js:je,ks:ke)

file=header+'v3_'+strtrim(string(itstep,format='(i4.4)'),1)+'.dat'
read_data,data,time,file,in,jn,kn
v3c=0.50*(float(data)+shift(float(data),0,0,-1))
v3c=v3c(is:ie,js:je,ks:ke)

file=header+'d_'+strtrim(string(itstep,format='(i4.4)'),1)+'.dat'
read_data,data,time,file,in,jn,kn
dc=float(data(is:ie,js:je,ks:ke))

file=header+'e_'+strtrim(string(itstep,format='(i4.4)'),1)+'.dat'
read_data,data,time,file,in,jn,kn
ec=float(data(is:ie,js:je,ks:ke))

pc=ec*float(gamma-1.D0)
;tempc=pc/dc*float(unit_v*unit_v/(rstar*unit_temp))
tempc=dc
ss0=size(dc)
for k=0,ss0[3]-1 do begin
    for j=0,ss0[2]-1 do begin
        for i=0,ss0[1]-1 do begin
		if dc[i,j,k] gt 0 then begin
		tempc[i,j,k]=pc[i,j,k]/dc[i,j,k]*float(unit_v*unit_v/(rstar*unit_temp))
		endif else begin
		tempc[i,j,k]=0
		endelse
        endfor
   endfor
endfor

;data to be returned
time=float(time*unit_time)
b1c=b1c*float(unit_b)
b2c=b2c*float(unit_b)
b3c=b3c*float(unit_b)
v1c=v1c*float(unit_v)
v2c=v2c*float(unit_v)
v3c=v3c*float(unit_v)
dc=dc*float(unit_rho)
tempc=tempc*float(unit_temp)
pc=pc*float(unit_rho*unit_v*unit_v)
x1c=float(x1b(is:ie))
x2c=float(x2b(js:je))
x3c=float(x3b(ks:ke))

end

PRO read_data,data,time,file,in,jn,kn,endian=endian

if n_elements(endian) eq 0 then endian='little'
in=long(85)
jn=long(85)
kn=long(85)
time=0.D0

print,'read ',file
if endian eq 'little' then begin
openr,/get_lun,unit,file
endif else begin
openr,/get_lun,unit,file,/swap_if_little_endian
endelse
readu,unit,time
readu,unit,in,jn,kn

is=3-1 & js=3-1 & ks=3-1
ie=is+(in-6) & je=js+(jn-6) & ke=ks+(kn-6)
ism1=is-1 & jsm1=js-1 & ksm1=ks-1
ism2=is-2 & jsm2=js-2 & ksm2=ks-2
iep1=ie+1 & jep1=je+1 & kep1=ke+1
iep2=ie+2 & jep2=je+2 & kep2=ke+2
iep3=ie+3 & jep3=je+3 & kep3=ke+3

data=dblarr(iep2-ism2+1,jep2-jsm2+1,kep2-ksm2+1)
readu,unit,data
free_lun,unit

return
end
