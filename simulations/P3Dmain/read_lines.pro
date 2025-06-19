pro read_lines,fn,number,nlines,Bfile,outpath = outpath,nx=nx,ny=ny,nz=nz
;;Name:
;;  read_lines 
;;
;;Purpose: 
;;  Get the coordinates and magnetic field of MFRs from the VTK file of Paraview
;;
;;Input Parameters:
;;      fn - the filenames of VTK data
;;      number - the total points number of VTK data
;;      nlines - the total number of magnetic field lines 
;;      Bfile - the filenames of magnetic field data (from such as AMRVAC);;;;.dat file or the file of b1,b2,b3?
;;      outpath - the filepath of output results
;;
;;History:
;;  2020-09-29 Version 2.0 created by Y. W. Ni
;;
;;Contact:
;;  nsfzniyiwei@163.com

num=long(number)
nl=long(nlines)
help,num

out_path = ''
if keyword_set(outpath) then out_path = outpath + '/'

file_mkdir,out_path+'MFR_lines/'
file_mkdir,out_path+'MFR_lines/IDLsave/'
file_mkdir,out_path+'MFR_lines/location_MFRs/'
file_mkdir,out_path+'MFR_lines/Bxyz_MFRs/'

;file='/home/ni/RBSL_EX/QSL/0.04_440_h+/mfr.vtk'
file=fn
openr,lun,file,/get_lun

dos=0
flag='POINTS'
while(dos ne 1) do begin
  tmp=''
  readf,lun,tmp
  dos=strpos(strupcase(tmp),strupcase(flag)) ge 0;;;; the STRPOS function returns the position of the search string from the beginning of the string ( 0 is the first character).
  endwhile
print,tmp

point_data_read=[]

num_PD=0
while(num_PD ne num*3) do begin
  tmp=''
  readf,lun,tmp
;  print,tmp
  PD_tmp=double(STRSPLIT(tmp, ' ', /EXTRACT))
  point_data_read=[point_data_read,PD_tmp]
  num_PD=n_elements(point_data_read)
  endwhile
print,tmp


dos=0
flag='LINES'
while(dos ne 1) do begin
  tmp=''
  readf,lun,tmp
  dos=strpos(strupcase(tmp),strupcase(flag)) ge 0
endwhile

print,tmp

nm='MFRs from Regularized Biot-Safar Laws'

Forwardstr={ID:'1 2 3',IT:'1 2 3'}
Backwardstr={ID:'1 2 3',IT:'1 2 3'}
MFR={Name:nm,ID:0,Forward:Forwardstr,Backward:Backwardstr}
MFRsys=replicate(MFR,nl)

for it=0,2*nl-1 do begin
  str=''
  readf,lun,str
  if it lt nl then begin
    MFRsys[it].ID=it
    MFRsys[it].Forward.ID=str
    endif else begin
      MFRsys[it-nl].Backward.ID=str
      endelse
  endfor

dos=0
flag='IntegrationTime'
while(dos ne 1) do begin
  tmp=''
  readf,lun,tmp
  dos=strpos(strupcase(tmp),strupcase(flag)) ge 0
  endwhile

  IT_data_read=[]
  num_IT=0
  while(num_IT ne num) do begin
    tmp=''
    readf,lun,tmp
    IT_tmp=double(STRSPLIT(tmp, ' ', /EXTRACT))
    IT_data_read=[IT_data_read,IT_tmp]
    num_IT=n_elements(IT_data_read)
    endwhile
  print,tmp
  
  IT_data=reform(IT_data_read,num)


close,lun
free_lun,lun


;print,num
point_data_0=reform(point_data_read,num*3)
point_data=dblarr(num,3)
for ipd=0,num-1 do begin
  point_data[ipd,*]=point_data_read[[3*long(ipd),3*long(ipd)+1,3*long(ipd)+2]]
  endfor


;help,point_data
openw,lun,out_path+'MFR_lines/point_data.txt',/get_lun
printf,lun,Transpose(point_data),format='(3(f9.4))
close,lun
free_lun,lun

for in=0,nl-1 do begin
  MFR_tmp=MFRsys[in]
  FID=long(STRSPLIT(MFR_tmp.Forward.ID, ' ', /EXTRACT))
  BID=long(STRSPLIT(MFR_tmp.Backward.ID, ' ', /EXTRACT))
  F_It=string(IT_data[FID[0]])
  B_It=string(IT_data[BID[0]])
  for iFn=1,n_elements(FID)-1 do begin
    F_IT=strcompress(F_IT+','+string(IT_data[FID[iFn]]),/remove_all)
  endfor
  for iBn=1,n_elements(BID)-1 do begin
    B_IT=strcompress(B_IT+','+string(IT_data[BID[iBn]]),/remove_all)
  endfor
  MFRsys[in].Forward.IT=F_IT
  MFRsys[in].Backward.IT=B_IT
endfor

save,file=out_path+'MFR_lines/IDLsave/MFR_VTK.sav',MFRsys,point_data,IT_data,num,nl





fn_B = Bfile;;;;Bfile has to be prepared?
;nx = 240;180;
;ny = 200;180;
;nz = 200;180;

Bx = dblarr(nx,ny,nz)
By = dblarr(nx,ny,nz)
Bz = dblarr(nx,ny,nz)
openr,lun,fn_B,/get_lun,/F77_UNFORMATTED,ERROR=ERR
readu,lun,Bx
readu,lun,By
readu,lun,Bz
close,lun
free_lun,lun

;Bx=Bx*2.0
;By=By*2.0
;Bz=Bz*2.0;;;;why

loadct,34
TVLCT, red, green, blue, /GET
loadct,0
cb=round(255.*findgen(nl)/(nl-1))
oSym= Obj_New('IDLgrSymbol',data=3)

 tlb=WIDGET_BASE(xsize=1200,ysize=900)
 DEVICE,get_screen_size=sz
 info=WIDGET_INFO(tlb,/GEOMETRY)
 tlb_XY=[info.SCR_XSIZE,info.SCR_YSIZE]
 offset=[sz-tlb_XY]/2
 WIDGET_CONTROL,tlb,XOFFSET=offset[0],YOFFSET=offset[1]
 WIDGET_CONTROL,tlb,map=0,/real
 prsbar=IDLITWDPROGRESSBAR(Group_LEADER=tlb,title='Progress',Cancel=cancelIn)

begintime = systime(1)
for iw=0,nl-1 do begin
  ;iw=0
  MFR_tmp=MFRsys[iw]
  FID=long(STRSPLIT(MFR_tmp.Forward.ID, ' ', /EXTRACT))
  BID=long(STRSPLIT(MFR_tmp.Backward.ID, ' ', /EXTRACT))
  F_It=double(STRSPLIT(MFR_tmp.Forward.IT, ',', /EXTRACT))
  B_It=double(STRSPLIT(MFR_tmp.Backward.IT, ',', /EXTRACT))
  location=[FID[1:n_elements(FID)-1],BID[1:n_elements(BID)-1]]
  It_a=[F_It[1:n_elements(FID)-1],B_It[1:n_elements(BID)-1]]
  
  if iw eq 0 then print,n_elements(FID)+n_elements(BID)-2,FID[0]+BID[0]
  
  index0=sort(It_a)
  location=location[index0]
  It0=It_a[index0]
  
  index1=uniq(It0)
  location=location[index1]
  It1=It0[index1]
  
  cd_MFR=point_data[location,*]
  
  XP=point_data[location,0]
  YP=point_data[location,1]
  ZP=point_data[location,2]
  
  rank_nl=floor(alog10(nl))+3
  if iw ne 0 then begin
    rank_iw=floor(alog10(iw))+1
    endif else begin
      rank_iw=1
      endelse
  rank_zero=rank_nl-rank_iw
  str0=strtrim(string(iw),1)
  str=str0
  
  for izo=0,rank_zero-1 do begin
    str=strcompress('0'+str,/remove_all)
    endfor
  
  openw,lun,out_path+'MFR_lines/location_MFRs/MFR_L_'+str+'.txt',/get_lun
  printf,lun,'MFRs from Regularized Biot-Safar Laws'
  printf,lun,'MFR ID: '+str+' Points number: '+strtrim(string(n_elements(location)),1)
  printf,lun,'Points data coordinate'
  printf,lun,'X(pixels) Y(pixels) Z(pixels) '
  printf,lun,Transpose(point_data[location,*]),format='(3(f10.3))
  close,lun
  free_lun,lun
  
  B_P=dblarr(n_elements(location),3)
  
  for il=0,n_elements(location)-1 do begin
    xc=XP[il] & yc=YP[il] & zc=ZP[il] 
    dx1=xc-floor(xc) & dy1=yc-floor(yc) & dz1=zc-floor(zc) 
    xv=dblarr(3,8)
    corner,floor(xc),floor(yc),floor(zc),Bx,By,Bz,xv
    B_P[il,0]=xitp(0,dx1,dy1,dz1,xv)
    B_P[il,1]=xitp(1,dx1,dy1,dz1,xv)
    B_P[il,2]=xitp(2,dx1,dy1,dz1,xv)
    endfor
    
  openw,lun,out_path+'MFR_lines/Bxyz_MFRs/MFR_B_'+str+'.txt',/get_lun
  printf,lun,'MFRs from Regularized Biot-Safar Laws'
  printf,lun,'MFR ID: '+str+' Points number: '+strtrim(string(n_elements(location)),1)
  printf,lun,'Magnetic Field'
  printf,lun,'    Bx        By      Bz     '
  printf,lun,Transpose(B_P),format='(3(f9.4))
  close,lun
  free_lun,lun
  
  
  Xplot3d,XP,YP,ZP,linestyle=2,thick=2,symbol = oSym,$
    XRANGE=[-8,8],YRANGE=[-10,10],ZRANGE=[0,12],$
    color=[red[cb[iw]],green[cb[iw]],blue[cb[iw]]],$
    /OVERPLOT
  
  endtime = systime(1)
  runtimes = endtime-begintime
  t_hr=floor(runtimes/3600.)
  t_mm=floor((runtimes-t_hr*3600)/60.)
  t_ss=runtimes mod 60  
  bar='Cost Time: '+strcompress(strtrim(string(t_hr),1)+'h'+strtrim(string(t_mm),1)+'m'+strmid(strtrim(string(t_ss),1),0,4)+'s',/remove_all)  
  print,format='($,"'+bar+'",a)', $
          string(13b)
    
  IF WIDGET_INFO(prsbar,/valid)THEN BEGIN
    IDLITWDPROGRESSBAR_SETVALUE,prsbar,((iw+1.)/nl)*100.
  ENDIF ELSE BEGIN
    tmp=DIALOG_MESSAGE('Cancel Percent'+STRING(((iw+1.)/nl)*100.)+'%',/info)
    Break
  ENDELSE
  WAIT,0.5
endfor

 Widget_Control,tlb,/destroy
print,''
print,'Completed!'
end
