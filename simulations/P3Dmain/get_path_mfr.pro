function get_path_mfr, xyz_in, Ninterpol

jmax=Ninterpol+1
size_xyz=size(xyz_in)
;s_i=dblarr(jmax)

x0=reform(xyz_in[0,*],size_xyz[2])
y0=reform(xyz_in[1,*],size_xyz[2])
z0=reform(xyz_in[2,*],size_xyz[2])
x1=dblarr(jmax)
y1=dblarr(jmax)
z1=dblarr(jmax)

x1=interpol(x0,indgen(size_xyz[2]),indgen(jmax)*size_xyz[2]/(1.*Ninterpol))
y1=interpol(y0,indgen(size_xyz[2]),indgen(jmax)*size_xyz[2]/(1.*Ninterpol))
z1=interpol(z0,indgen(size_xyz[2]),indgen(jmax)*size_xyz[2]/(1.*Ninterpol))

time_hmi = '10-May-2012 00:00:00'
radius_sun = 6.955e10         ; sun radius in cm
ang = pb0r(time_hmi,/arcsec)
arcsec2cm=radius_sun/ang[2]   ; 1 arcsec in cm
dx=4.0d*0.5d*arcsec2cm
dx=dx/1.0d9

;x1=x1*dx
;y1=y1*dx
;z1=z1*dx

;s_i[1:jmax-1]=s_i[1:jmax-1]+ $
;              sqrt((x1[1:jmax-1]-x1[0:jmax-2])^2 + $ 
;                   (y1[1:jmax-1]-y1[0:jmax-2])^2 + $ 
;                   (z1[1:jmax-1]-z1[0:jmax-2])^2) 

xyz_out=dblarr(3,jmax)
xyz_out[0,*]=x1
xyz_out[1,*]=y1
xyz_out[2,*]=z1

return,xyz_out

end
