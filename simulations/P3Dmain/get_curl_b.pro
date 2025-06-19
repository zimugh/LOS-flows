function get_curl_B,Bfile,nx,ny,nz,$
         path=path,alpha_B=alpha_B,Btotal=Btotal,curBtotal=curBtotal,$
         Bx=Bx,By=By,Bz=Bz,curBx=curBx,curBy=curBy,curBz=curBz

;=============================================
; Parameter needed to be set for your own case
;=============================================
time_hmi = '10-May-2012 00:00:00'
ang = pb0r(time_hmi,/arcsec)
print,ang[2]
radius_sun = 69.55          ; sun radius in 10 Mm
arcsec2tenMm=radius_sun/ang[2] ; 1 arcsec in 10 Mm
dx=4.0*0.50d*arcsec2tenMm   ; x spatial resolution of the preprocessed magnetic field
dy=dx                       ; y spatial resolution
dz=dx                       ; z spatial resolution
;==============================================
print,dx,dy,dz

fn_B = Bfile
Bx = dblarr(nx,ny,nz)
By = dblarr(nx,ny,nz)
Bz = dblarr(nx,ny,nz)
openr,lun,fn_B,/get_lun,/F77_UNFORMATTED,ERROR=ERR
readu,lun,Bx
readu,lun,By
readu,lun,Bz
close,lun
free_lun,lun

xx=findgen(nx)*dx
yy=findgen(ny)*dy
zz=findgen(nz)*dz

curl_xyz, Bx, By, Bz, xx, yy, zz, curBx, curBy, curBz

alpha_B=sqrt(((curBx^2)+(curBy^2)+(curBz^2))/((Bx^2)+(By^2)+(Bz^2)))

Btotal=sqrt((Bx^2)+(By^2)+(Bz^2))
curBtotal=sqrt((curBx^2)+(curBy^2)+(curBz^2))

return,1

;save,file=path+'/curl.sav',Bx,By,Bz,curBx,curBy,curBz,alpha_B




end