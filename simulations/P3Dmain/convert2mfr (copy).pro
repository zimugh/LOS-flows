function convert2mfr, MFR_ID, geo_path = geo_path, plt_path = plt_path, out_path = out_path, $
                      Ngrids = Ngrids, Radsyn = Radsyn, $
                      XYZplane = XYZplane, XV = XV, YV = YV, ZV = ZV


; coordinates of a 3D magnetic field line
geo_fn = geo_path + '/path.dat'
ngeo = numlines(geo_fn)
xyz_out = dblarr(3,ngeo)
openr,lun,geo_fn,/get_lun
readf,lun,format='(3e25.16)',xyz_out
close,lun
free_lun,lun

x_i = reform(xyz_out[0,*],ngeo)
y_i = reform(xyz_out[1,*],ngeo)
z_i = reform(xyz_out[2,*],ngeo)

Ninterpol = ngeo - 1
s_i = dblarr(ngeo)
dpoint = dblarr(3,Ninterpol)
dpoint = xyz_out[*,1:Ninterpol]-xyz_out[*,0:Ninterpol-1]
for i = 1,Ninterpol do begin
  s_i[i] = s_i[i-1] + DOUBLE(sqrt(TOTAL(dpoint[*,i-1]^2, 1)))
  endfor

print,minmax(s_i),s_i[0],s_i[Ninterpol]

; read *.plt
plt_fn = file_search(plt_path + '/*.plt')
nplt = N_elements(plt_fn)

Tgirds = dindgen(nplt)
dt = 1.0d ; savetime
Tgirds = (Tgirds*dt*86)/3600. ; 86 is the alfven time
len = DOUBLE(TOTAL(sqrt(TOTAL(dpoint^2, 1)),1))    ; the length of the loop

rho = dblarr(Ngrids,nplt)
v1 = dblarr(Ngrids,nplt)
p = dblarr(Ngrids,nplt)
Te = dblarr(Ngrids,nplt)

rho_i = dblarr(Ngrids)
v1_i = dblarr(Ngrids)
p_i = dblarr(Ngrids)
Te_i = dblarr(Ngrids)

s_o = (dindgen(Ngrids)/(Ngrids-1))*len
rho_o = dblarr(Ngrids)
v1_o = dblarr(Ngrids)
p_o = dblarr(Ngrids)
Te_o = dblarr(Ngrids)
aia_ra = dblarr(Ngrids,nplt)

for i = 0,nplt-1 do begin
  lines = numlines(plt_fn[i])
  character1=''
  character2=''
  data_plt = dblarr(5,lines-2)
  openr,var_lun,plt_fn[i],/get_lun
  readf,var_lun,character1
  readf,var_lun,character2
  readf,var_lun,data_plt
  close,var_lun
  free_lun,var_lun
  s_b = reform(data_plt[0,*],lines-2)
  rho_b = reform(data_plt[1,*],lines-2)
  v1_b = reform(data_plt[2,*],lines-2)
  p_b = reform(data_plt[3,*],lines-2)
  Te_b = reform(data_plt[4,*],lines-2)
  
  ; Linterp
  ;; 1st step: obtaining the uniform grid physical points
  LINTERP,s_b,rho_b,s_o,rho_o
  LINTERP,s_b,v1_b,s_o,v1_o
  LINTERP,s_b,p_b,s_o,p_o
  LINTERP,s_b,Te_b,s_o,Te_o
  
  ;; 2st step: Matching the coordinates with the physical parameters
  LINTERP,s_i,x_i,s_o,x_o 
  LINTERP,s_i,y_i,s_o,y_o 
  LINTERP,s_i,z_i,s_o,z_o 
  
  rho[*,i] = rho_o
  v1[*,i] = v1_o
  p[*,i] = p_o
  Te[*,i] = Te_o
  
  rho[*,i] = rho[*,i]*1.0d9
  Te[*,i] = Te[*,i]*1.0d6
  
  if keyword_set(Radsyn) then begin
    ; radiation synthesis
    lnn = alog10(double(rho[*,i]))
    ltt = alog10(double(te[*,i]))
    channel = [Radsyn]
    ;aiarp = aia_resp(ltt,channel,logn = lnn, $
    ;                 file = '/ssw/gen/idl_libs/astron/radiation/aia_resp.xdr')
    aiarp = aia_resp(ltt,channel,logn = lnn,file = './radiation/aia_resp.xdr')
    aia_ra[*,i] = aiarp.R
    endif
    
endfor

xyz_o = dblarr(3,Ngrids)
xyz_o[0,*] = x_o
xyz_o[1,*] = y_o
xyz_o[2,*] = z_o
print,min(rho),max(rho)
rho = transpose(rho)
v1 = transpose(v1)
p = transpose(p)
Te = transpose(Te)

if keyword_set(Radsyn) then begin
  aia_ra = transpose(aia_ra)
  ops = aia_ra
  endif else begin
    ops = alog10(rho)
    endelse
    
file_mkdir, out_path + '/sav'
save,filename = out_path + '/sav/'+'xyz'+'.sav', xyz_o
save,filename = out_path + '/sav/'+'Tgirds'+'.sav', Tgirds
save,filename = out_path + '/sav/'+'s_o'+'.sav', s_o
save,filename = out_path + '/sav/'+'rho'+'.sav', rho
save,filename = out_path + '/sav/'+'v1'+'.sav', v1
save,filename = out_path + '/sav/'+'p'+'.sav', p
save,filename = out_path + '/sav/'+'Te'+'.sav', Te


Ngrids_ID = indgen(Ngrids)
for it = 0, nplt-1 do begin
  openw, 1, out_path + '/'+strtrim(string(MFR_ID),2)+'mfr_path'+strtrim(string(it),2) + '.vtk'
  printf, 1, '# vtk DataFile Version 3.0'
  printf, 1, 'helix'
  printf, 1, 'ASCII'
  printf, 1, 'DATASET UNSTRUCTURED_GRID'
  printf, 1, 'POINTS'+' '+strtrim(string(Ngrids),2)+' '+'double'
  printf, 1, format='(3f)', xyz_o
  
  printf, 1, 'CELLS 1' + '  '+strtrim(string(Ngrids+1),2)
  printf, 1, Ngrids, Ngrids_ID
  
  printf, 1, 'CELL_TYPES  1' 
  printf, 1, '4'
  
  printf, 1, 'POINT_DATA'+' '+strtrim(string(Ngrids),2)
  printf, 1, 'SCALARS'+' '+'rho'+' '+'double'+' '+'1'
  printf, 1, 'LOOKUP_TABLE  table1'
  printf, 1, format='(1f)', reform(ops[it,*],Ngrids)
  
  close,1

endfor

return,1
end




