pro P3D_render, MFR_path = MFR_path, out_path = out_path, file_xdr = file_xdr

Tmax = 600
plt_path = file_search(MFR_path + '/MFR*')
Nplt = n_elements(plt_path)
xyz_F = []
rho_F = []
Tem_F = []

print,'============================================='
print,'         Loading simulation data             '

 for i = 0 , Nplt -1 do begin
    print, format='($,"Progress:  %",i3,a)', $
         double((i+1.)*100.)/double(Nplt), $
            string(13b)
    file_xyz = plt_path[i] + '/sav/' + 'xyz' + '.sav'
    file_rho = plt_path[i] + '/sav/' + 'rho' + '.sav'
    file_Te = plt_path[i] + '/sav/' + 'Te' + '.sav'
    restore, file_xyz
    restore, file_rho
    restore, file_Te
    
    print,size(rho)
    
;    if (size(rho))[1] eq 601 then begin
    if (size(rho))[1] eq 801 then begin
      xyz_F = [[xyz_F],[xyz_o]]
      rho_F = [[rho_F],[rho]]
      Tem_F = [[Tem_F],[Te]]
      endif
    endfor
X = reform(xyz_F[0,*])
Y = reform(xyz_F[1,*])
Z = reform(xyz_F[2,*])
print,'  '
print,'Loading completed  √'
help,X,Y,Z,rho_F,Tem_F
print,'============================================='
print,'  '
print,'  '

print,'  '
print,'  '
print,'============================================='
print,'           Grid initialization               '
dx = 4.d-2
dy = dx
dz = 4.d-2
Gx = 7. + findgen(ceil(20./dx)+1)*dx
Gy = 2.5 + findgen(ceil(21.5/dy)+1)*dy
Gz = 0.1 + findgen(ceil(11./dz)+1)*dz;;lt how to change these parameters

print,'The range of X: ' + strmid(strtrim(string(min(Gx)),2),0,5) + $
                           ' --> ' + $
                           strmid(strtrim(string(max(Gx)),2),0,5) + $
                           '    dx = ' + $
                           strmid(strtrim(string(dx/(10.^floor(alog10(dx)))),2),0,3) + $
                           'e' +strtrim(string(floor(alog10(dx))),2)
print,'The range of Y: ' + strmid(strtrim(string(min(Gy)),2),0,5) + $
                           ' --> ' + $
                           strmid(strtrim(string(max(Gy)),2),0,5) + $
                           '    dy = ' + $
                           strmid(strtrim(string(dy/(10.^floor(alog10(dy)))),2),0,3) + $
                           'e' +strtrim(string(floor(alog10(dy))),2)
print,'The range of Z: ' + strmid(strtrim(string(min(Gz)),2),0,5) + $
                           ' --> ' + $
                           strmid(strtrim(string(max(Gz)),2),0,5) + $
                           '    dz = ' + $
                           strmid(strtrim(string(dz/(10.^floor(alog10(dz)))),2),0,3) + $
                           'e' +strtrim(string(floor(alog10(dz))),2)

print,''
print,'  Interpolating scattered data into 3D grid =>'
F_root = GRID3(X, Y, Z, x*0.+1., Gx, Gy, Gz, DTOL=1.d-6, /Grid)
print,'  Intial interpolation completed √'
print,' '
print,'  Using the Nearest Neighbor Mean(NNM) method'
pos=where(F_root eq 1.)
F_root_dims=size(F_root,/dimensions)
coord_root=array_indices(F_root_dims,pos,/dimensions)  
Nroots = N_ELEMENTS(pos)
Fs_rho = dblarr(F_root_dims[0],F_root_dims[1],F_root_dims[2])
Fs_Tem = dblarr(F_root_dims[0],F_root_dims[1],F_root_dims[2])
zgrid_top = ceil((max(coord_root[2,*])+3.))

print,'@Generating the NNM Library =>'
lib_rho = dblarr(Tmax + 1,Nroots)
lib_tem = dblarr(Tmax + 1,Nroots)
for i = 0, Nroots-1 do begin
  print, format='($,"Progress:  %",i3,a)', $
         double((i+1.)*100.)/double(Nroots), $
            string(13b)
  x_root = Gx[coord_root[0,i]]
  y_root = Gy[coord_root[1,i]]
  z_root = Gz[coord_root[2,i]]
  dis_root = sqrt((X - x_root)^2 + (Y - y_root)^2 + (Z - z_root)^2)
  pos_root = where(dis_root eq min(dis_root))
  if n_elements(pos_root) gt 1 then begin
    lib_rho[*,i] = alog10(mean(rho_F[*,pos_root], DIMENSION= 2))
    lib_tem[*,i] = alog10(mean(tem_F[*,pos_root], DIMENSION= 2))
    endif else begin
      lib_rho[*,i] = alog10(rho_F[*,pos_root])
      lib_tem[*,i] = alog10(tem_F[*,pos_root])
      endelse
  endfor

print,' '
print,'The NNM Library completed  √'
print,'Parameter information =>'
help,lib_rho,lib_tem,Fs_rho,Fs_tem
print,'============================================='
print,'  '
print,'  '

print,'  '
print,'  '
print,'============================================='
print,'           Radiation synthesis               '
print,'  In this program, the radiation of prominence '
print,'is synthesized by using pure emission model,   '
print,'and the radiation of filament is synthesized by'
print,'solving the optically thick radiation transfer '
print,'equation.'  
Ibg = 2.d4
segmaH = 5.16d-20   ;cm^2
mH = 1.6726d-24     ;g
zone = 30
print,'  '
print,'Parameter information =>'
print,'Background radiation value: ' + strtrim(string(Ibg),2) + ' DN s^-1'
print,'Hydrogen scattering cross section: ' + $
                           strmid(strtrim(string(segmaH/(10.^floor(alog10(segmaH)))),2),0,4) + $
                           'e' +strtrim(string(floor(alog10(segmaH))),2) + ' cm^2'
print,'Hydrogen atomic mass: ' + $
                           strmid(strtrim(string(mH/(10.^floor(alog10(mH)))),2),0,4) + $
                           'e' +strtrim(string(floor(alog10(mH))),2) + ' g'
print,'Bottom height: ' + strmid(strtrim(string(Gz(zone)*10.),2),0,4) + ' Mm'
I_pro_YZ = dblarr(n_elements(Gy),n_elements(Gz),Tmax + 1)
I_pro_XZ = dblarr(n_elements(Gx),n_elements(Gz),Tmax + 1)
I_fil_XY = dblarr(n_elements(Gx),n_elements(Gy),Tmax + 1)
print,'  '
print,'Radiation synthesis =>'
for itm = 0, Tmax do begin
  print, format='($,"Progress:  %",i3,a)', $
         double((itm+1.)*100.)/double(Tmax), $
            string(13b)
  Fs_rho[pos] = reform(lib_rho[itm,*])
  Fs_tem[pos] = reform(lib_tem[itm,*])
  lnn = double(Fs_rho)
  ltt = double(Fs_Tem)
  channel = [171]
  aiarp = aia_resp(ltt,channel,logn = lnn, $
                       file = file_xdr)
  
  I_pro_YZ[*,*,itm] = total(aiarp.R,1)
  I_pro_XZ[*,*,itm] = total(aiarp.R,2)
  Irad = dblarr(F_root_dims[0],F_root_dims[1],F_root_dims[2])
  Irad[*,*,0:zone-1] = Ibg
  for iz = zone, zgrid_top - 1 do begin
    Irad[*,*,iz] = reform(Irad[*,*,iz-1]) * (1 - $
                   dz * 1.d9 * segmaH * double(10.^reform(Fs_rho[*,*,iz-1]))) + $
                   (dz * 1.d1) * $
                   reform((reform(double(aiarp.R)))[*,*,iz-1])
  endfor
  
  I_fil_XY[*,*,itm] = Irad[*,*,zgrid_top - 1]
  
  
endfor
print,'  '
print,'Radiation synthesis completed!'
print,'  '
print,'Parameter information =>'
help,I_pro_YZ, I_pro_XZ
print,'Variables saving =>'
save, I_pro_YZ, I_pro_XZ, I_fil_XY, file = out_path +'/Iradiation.sav'
print,'All variables have been saved at '
print, out_path +'/Iradiation.sav'
print,'============================================='
print,'  '
print,'  '

end
