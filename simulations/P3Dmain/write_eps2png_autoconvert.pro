pro write_eps2png_autoconvert
nt = 601
fpath = '/mnt/P3D/T500'
openw,lun, fpath + '/eps2png', /get_lun
for i = 0, nt-1 do begin
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
  printf,lun,'convert '+ fpath +'/movie/XY/eps/XYpng'+str+'.eps '+ fpath +'/movie/XY/png/XYpng'+str+'.png'
  printf,lun,'convert '+ fpath +'/movie/XZ/eps/XZpng'+str+'.eps '+ fpath +'/movie/XZ/png/XZpng'+str+'.png'
  printf,lun,'convert '+ fpath +'/movie/YZ/eps/YZpng'+str+'.eps '+ fpath +'/movie/YZ/png/YZpng'+str+'.png'
  endfor
  close,lun
  free_lun,lun

; ffmpeg -r 60 -f image2 -i XYpng%05d.png XYversion.mp4

end
