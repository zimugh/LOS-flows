pro rad2png
;初始化颜色及色表
TVLCT,0,0,0      
bcolors=!p.background
!p.background=!p.color
!p.color=bcolors
print,!p.background,!p.color

;=============================================
; Parameter needed to be set for your own case
;=============================================
time_hmi = '10-May-2012 00:00:00'
ang = pb0r(time_hmi,/arcsec)
radius_sun = 69.55          ; sun radius in 10 Mm
nx = 240                    ; x-size of the preprocessed magnetic field
ny = 200                    ; y-size
arcsec2tenMm=radius_sun/ang[2] ; 1 arcsec in 10 Mm

fpath = '/mnt/P3D/T500'
fsav = fpath + '/Iradiation.sav'
restore,file = fsav, /ver
aia_lct, rr, gg, bb, wave=171, /load
pos = [.12,.12,.90,.90]
tA = (8.6d1)/3600.
Tpp = strmid(strtrim(string(dindgen(601)*tA),2),0,4)

dx = 4.d-2/arcsec2tenMm
dy = dx/arcsec2tenMm
dz = 4.d-2*10.
Stx = 7./arcsec2tenMm
Sty = 2.5/arcsec2tenMm 
Stz = 0.1*10. 
nt = 601

window,0,xs=1000,ys=1000
for i=0,nt - 1 do begin
    print, format='($,"Progress:  %",i3,a)', $
         double((i+1.)*100.)/double(nt), $
            string(13b)
            
  rank_nt=floor(alog10(nt))+3
  if i ne 0 then begin
    rank_i=floor(alog10(i))+1
    endif else begin
      rank_i=1
      endelse
  rank_zero=rank_nt-rank_i
  str0=strtrim(string(i),2)
  str=str0
  for izo=0,rank_zero-1 do begin
    str=strcompress('0'+str,/remove_all)
    endfor

  set_plot, 'ps'
  device, file = fpath+'/movie/XY/eps/XYpng' + str + '.eps', /inches,/color,bits=8,xsize=12,ysize=12,/TT_FONT
  plot_image, alog10(reform(I_FIL_XY[*,*,i])>1.), min=3., max=4.5, $
              ori = [Stx-dx/2.,Sty-dy/2.], sca = [dx,dy], pos = pos, $
              xtit = 'Solar X (arcsec)', ytit = 'Solar Y (arcsec)', tit = 'Time: ' + Tpp[i] + ' hr'
  device,/close
  set_plot,'x'
  
  set_plot, 'ps'
  device, file = fpath+'/movie/XZ/eps/XZpng' + str + '.eps', /inches,/color,bits=8,xsize=12,ysize=12,/TT_FONT
  plot_image, alog10(reform(I_PRO_XZ[*,*,i])>1.), min=1.5, max=4.5, $
              ori = [Stx-dx/2.,Stz-dz/2.], sca = [dx,dz], pos = pos, $
              xtit = 'Solar X (arcsec)', ytit = 'Height Z (Mm)', tit = 'Time: ' + Tpp[i] + ' hr'
  device,/close
  set_plot,'x'

 
  set_plot, 'ps'
  device, file = fpath+'/movie/YZ/eps/YZpng' + str + '.eps', /inches,/color,bits=8,xsize=12,ysize=12,/TT_FONT  
  plot_image, alog10(reform(I_PRO_YZ[*,*,i])>1.), min=1.5, max=4.5, $
              ori = [Sty-dx/2.,Stz-dz/2.], sca = [dy,dz], pos = pos, $
              xtit = 'Solar Y (arcsec)', ytit = 'Height Z (Mm)', tit = 'Time: ' + Tpp[i] + ' hr'
  device,/close
  set_plot,'x'
  
  endfor


;初始化颜色及色表
TVLCT,0,0,0      
bcolors=!p.background
!p.background=!p.color
!p.color=bcolors
print,!p.background,!p.color
end
