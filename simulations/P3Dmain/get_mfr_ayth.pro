function get_MFR_ayth,path = path, $
                      Bfile = Bfile, nx = nx, ny = ny, nz = nz, $
                      ayth_p = ayth_p

files=file_search(path+'/MFR_L_*.txt')
z_pos=2.0
get_curl_Bs = get_curl_B(Bfile,nx,ny,nz,$
         path=path,alpha_B=alpha_B,Btotal=Btotal,curBtotal=curBtotal,$
         Bx=Bx,By=By,Bz=Bz,curBx=curBx,curBy=curBy,curBz=curBz)

ayth=dblarr(6,n_elements(files))

 tlb=WIDGET_BASE(xsize=1200,ysize=900)
 DEVICE,get_screen_size=sz
 info=WIDGET_INFO(tlb,/GEOMETRY)
 tlb_XY=[info.SCR_XSIZE,info.SCR_YSIZE]
 offset=[sz-tlb_XY]/2
 WIDGET_CONTROL,tlb,XOFFSET=offset[0],YOFFSET=offset[1]
 WIDGET_CONTROL,tlb,map=0,/real
 prsbar=IDLITWDPROGRESSBAR(Group_LEADER=tlb,title='Progress',Cancel=cancelIn)


for i=0,n_elements(files)-1 do begin
  
  fn=files[i]
  openr,lun,fn,/get_lun
  dos=0
  flag0='Points number:'
  flag='X(pixels) Y(pixels) Z(pixels) '
  while(dos ne 1) do begin
    tmp=''
    readf,lun,tmp
    dos0=strpos(strupcase(tmp),strupcase(flag0))
    if dos0 ge 0 then begin
      dos_tmp=STRSPLIT(tmp, ' ', /EXTRACT)
      num=fix(dos_tmp(n_elements(dos_tmp)-1))
    endif
    dos=strpos(strupcase(tmp),strupcase(flag)) ge 0
  endwhile
;print,tmp

point_data_read=[]
num_PD=0
while(num_PD ne num*3) do begin
  tmp=''
  readf,lun,tmp
  PD_tmp=double(STRSPLIT(tmp, ' ', /EXTRACT))
  point_data_read=[[point_data_read],[PD_tmp]]
  num_PD=n_elements(point_data_read)
  endwhile
;print,tmp
;help,point_data_read
close,lun
free_lun,lun

floc=point_data_read[2,*]-z_pos
ffl=where(floc ge 0)

;;for min
lb=floc[min(ffl)-1];-0.1
rb=floc[min(ffl)];0.05
res=((0-lb)/(rb-lb));0.1/0.15 (x1,y1,1.9) (x2,y2,2.05) (xc,yc,2)
xc=res*(point_data_read[0,min(ffl)]-point_data_read[0,min(ffl)-1])+point_data_read[0,min(ffl)-1]
yc=res*(point_data_read[1,min(ffl)]-point_data_read[1,min(ffl)-1])+point_data_read[1,min(ffl)-1]
zc=z_pos
dx1=xc-floor(xc) & dy1=yc-floor(yc) & dz1=zc-floor(zc) 
    xv=dblarr(3,8)
    corner,floor(xc),floor(yc),floor(zc),curBx,curBy,curBz,xv
    curBx0=xitp(0,dx1,dy1,dz1,xv)
    curBy0=xitp(1,dx1,dy1,dz1,xv)
    curBz0=xitp(2,dx1,dy1,dz1,xv)
    corner,floor(xc),floor(yc),floor(zc),Bx,By,Bz,xv
    Bx0=xitp(0,dx1,dy1,dz1,xv)
    By0=xitp(1,dx1,dy1,dz1,xv)
    Bz0=xitp(2,dx1,dy1,dz1,xv)
ayth[0,i]=sqrt(((curBx0^2)+(curBy0^2)+(curBz0^2))/((Bx0^2)+(By0^2)+(Bz0^2)))
ayth[3,i]=(Bx0^2)+(By0^2)+(Bz0^2)

;;for max
lb=floc[max(ffl)-1]
rb=floc[max(ffl)]
res=((0-lb)/(rb-lb))
xc=res*(point_data_read[0,max(ffl)]-point_data_read[0,max(ffl)-1])+point_data_read[0,max(ffl)-1]
yc=res*(point_data_read[1,max(ffl)]-point_data_read[1,max(ffl)-1])+point_data_read[1,max(ffl)-1]
zc=z_pos
dx1=xc-floor(xc) & dy1=yc-floor(yc) & dz1=zc-floor(zc) 
    xv=dblarr(3,8)
    corner,floor(xc),floor(yc),floor(zc),curBx,curBy,curBz,xv
    curBx0=xitp(0,dx1,dy1,dz1,xv)
    curBy0=xitp(1,dx1,dy1,dz1,xv)
    curBz0=xitp(2,dx1,dy1,dz1,xv)
    corner,floor(xc),floor(yc),floor(zc),Bx,By,Bz,xv
    Bx0=xitp(0,dx1,dy1,dz1,xv)
    By0=xitp(1,dx1,dy1,dz1,xv)
    Bz0=xitp(2,dx1,dy1,dz1,xv)
ayth[1,i]=sqrt(((curBx0^2)+(curBy0^2)+(curBz0^2))/((Bx0^2)+(By0^2)+(Bz0^2)))
ayth[4,i]=(Bx0^2)+(By0^2)+(Bz0^2)

ayth[2,i]=ayth[0,i]/ayth[1,i]
ayth[5,i]=ayth[2,i]*ayth[3,i]/ayth[4,i]
ayth_p=exp(abs(alog(ayth[5,*])))

  IF WIDGET_INFO(prsbar,/valid)THEN BEGIN
    IDLITWDPROGRESSBAR_SETVALUE,prsbar,((i+1.)/n_elements(files))*100.
  ENDIF ELSE BEGIN
    tmp=DIALOG_MESSAGE('Cancel Percent'+STRING(((i+1.)/n_elements(files))*100.)+'%',/info)
    Break
  ENDELSE
  WAIT,0.5

endfor
return,ayth

end
