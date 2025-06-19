pro make_AMRVACs,outpath = outpath, amrvac_path = amrvac_path, path_MFR = path_MFR,$
                 ASY_Heat=ASY_Heat,AB2=AB2,Bfile = Bfile, nx = nx, ny = ny, nz = nz

fn0=amrvac_path+'/amrvac.par'
fn1=amrvac_path+'/mod_usr.t'

str_amrvac=[]
openr,lun,fn0,/get_lun
dos=0
flag='/IDL_STOP'
while(dos ne 1) do begin
  tmp=''
  readf,lun,tmp
  dos=strpos(strupcase(tmp),strupcase(flag)) ge 0
  dos_base_filename=strpos(strupcase(tmp),strupcase('base_filename')) lt 0
  dos_xprobmax1=strpos(strupcase(tmp),strupcase('xprobmax1        =  34.258009068569379d0')) lt 0
  if (1-dos)*dos_base_filename*dos_xprobmax1 eq 1 then str_amrvac=[[str_amrvac],[tmp]]
  endwhile
close,lun
free_lun,lun

str_mod_usr=[]
openr,lun,fn1,/get_lun
dos=0
dos_rang=1
flag=' /IDL_STOP'
while(dos ne 1) do begin
  tmp=''
  readf,lun,tmp
  dos=strpos(strupcase(tmp),strupcase(flag)) ge 0
  dos_define=strpos(strupcase(tmp),strupcase('double precision :: lQgrid(ixI^S),lambda,lQ0,lH0,tramp')) lt 0
  ;dos_LQ0=strpos(strupcase(tmp),strupcase('lQ0=4.0d-2/heatunit')) lt 0
  dos_LQ0=strpos(strupcase(tmp),strupcase('lQ0=5.0d-3/heatunit')) lt 0
  dos_LQ_start=strpos(strupcase(tmp),strupcase('where(x(ixO^S,1) .lt. lH0)')) ge 0
  dos_LQ_end=strpos(strupcase(tmp),strupcase('end subroutine getlQ')) ge 0
  if dos_LQ_start eq 1 then dos_rang=dos_rang-1
  if dos_LQ_end eq 1 then dos_rang=dos_rang+1
  if (1-dos)*dos_define*dos_LQ0*dos_rang eq 1 then str_mod_usr=[[str_mod_usr],[tmp]]
  endwhile
close,lun
free_lun,lun




file_mkdir,outpath+'/P3D/'
file_MFR=file_search(path_MFR+'/MFR_L_*.txt')
N_MFR=n_elements(file_MFR)

;=============================================
; Parameter needed to be set for your own case
;=============================================
time_hmi = '03-Aug-2023 16:36:00'
ang = pb0r(time_hmi,/arcsec)
radius_sun = 69.55          ; sun radius in 10 Mm
arcsec2tenMm=radius_sun/ang[2] ; 1 arcsec in 10 Mm
dlen=4.0*0.50d*arcsec2tenMm   
;==============================================

for i_MFR=0,N_MFR-1 do begin

  str0=strtrim(string(i_MFR),2)
  str=str0
file_mkdir,outpath+'/P3D/MFR_'+str
file_mkdir,outpath+'/P3D/data/MFR_'+str
fn00=outpath+'/P3D/MFR_'+str+'/amrvac.par'
fn10=outpath+'/P3D/MFR_'+str+'/mod_usr.t'
fn20=outpath+'/P3D/MFR_'+str+'/path.dat'
fn30=outpath+'/P3D/MFR_'+str+'/job'


pos_base_filename=where(strpos(strupcase(str_amrvac),strupcase($
        '&filelist'$
         )) ge 0)
pos_xprobmax1=where(strpos(strupcase(str_amrvac),strupcase($
        'xprobmin1        = 0.d0'$
         )) ge 0)
pos_defined=where(strpos(strupcase(str_mod_usr),strupcase($
        '!> where the lQ defined'$
         )) ge 0)
pos_LQ0=where(strpos(strupcase(str_mod_usr),strupcase($
        '!> where the lQ0 setting'$
         )) ge 0)
pos_LQ=where(strpos(strupcase(str_mod_usr),strupcase($
        '!> where the lQ setting'$
         )) ge 0)





str_base_filename=strarr(1,1)
;str_base_filename[*,0]="        base_filename = '/gpfsdata/home/zhouyuhao/nyw/P3D/P3D_test/P3D/"+$
;                                                "data/MFR_"+str+"/pro'"
;str_base_filename[*,0]="        base_filename = '/share/home/chenpf/nyw/P3D/"+$
;                                                "data/MFR_"+str+"/pro'";;;;need to change path
                                                
;str_base_filename[*,0]="        base_filename = '/home/test/ltjupyter/work/1paperwork/20231108s3d_simulation/P3DS/mfrdata/P3D/"+$
;                                                "data/MFR_"+str+"/pro'";;;;need to change path
                                                
str_base_filename[*,0]="        base_filename = '../data/MFR_"+str+"/pro'";;;;need to change path

;print,file_MFR[i_MFR]
openr,lun,file_MFR[i_MFR],/get_lun
dos_tmp_len=0
dos_tmp_num=0
while(dos_tmp_len ne 1) do begin
  tmp_Len=''
  readf,lun,tmp_Len
  dos_tmp_Len=strpos(strupcase(tmp_Len),strupcase('X(pixels) Y(pixels) Z(pixels) ')) ge 0
  dos_tmp_num=strpos(strupcase(tmp_Len),$
                    strupcase('Points number: ')) ge 0
  if dos_tmp_num ne 0 then begin
    tmp_num_split=STRSPLIT(tmp_Len, ' ', /EXTRACT)
    num_len=LONG(tmp_num_split[n_elements(tmp_num_split)-1])
    endif
  endwhile
point_data=dblarr(3,num_len)
readf,lun,point_data
close,lun
free_lun,lun

Ninterpol=5000;;;;why 5000
xyz_out = get_path_mfr(point_data, Ninterpol)

dpoint=dblarr(3,Ninterpol)
dpoint=xyz_out[*,1:Ninterpol]-xyz_out[*,0:Ninterpol-1]
len=DOUBLE(TOTAL(sqrt(TOTAL(dpoint^2, 1)),1))
Plen=floor(alog10(len))
Nlen=len/10^(Plen)
str_len='        xprobmax1        = '+$
          strcompress(strtrim(string(Nlen),2)+$
          'd'+strtrim(string(Plen),2),/remove_all)
str_len=reform(str_len,1,1)
str_amrvac_in=[[str_amrvac[*,0:pos_base_filename]],$
               [str_base_filename],$
               [str_amrvac[*,pos_base_filename+1:pos_xprobmax1]],$
               [str_len],$
               [str_amrvac[*,pos_xprobmax1+1:n_elements(str_amrvac)-1]]];lt;why +1


Q0m = 1.
Q1m = 1.

if keyword_set(ASY_Heat) then begin
  ayth=get_MFR_ayth(path = path_MFR, $
                    Bfile = Bfile, nx = nx, ny = ny, nz = nz, $
                    ayth_p = ayth_p)

  
if keyword_set(AB2) then begin
  QM = 1828.4027;;;; why this number
  ;QM=mean((ayth[0,*]*ayth[3,*]+ayth[1,*]*ayth[4,*])/2.)
  ;print,QM
  Q0=ayth[0,i_MFR]*ayth[3,i_MFR]/QM> 0.25d-2
  Q1=ayth[1,i_MFR]*ayth[4,i_MFR]/QM> 0.25d-2
endif else begin
  QM=1.
  Q0=1.
  Q1=1.
endelse
      
      
      
  print,'Num: ',i_MFR, '  ', Q0, '  ',Q1 
  if Q0 le Q0m then Q0m = Q0
  if Q1 le Q1m then Q1m = Q1
  str_define=strarr(1,1)
  str_define[*,0]='    double precision :: lQgrid(ixI^S),lambda,lQ0,lQ1,lH0,tramp'
  str_LQ_val=strarr(1,2)
  str_LQ_val[*,0]='    lQ0=('+strtrim(string(Q0),2)+')*4.0d-2/heatunit'
  str_LQ_val[*,1]='    lQ1=('+strtrim(string(Q1),2)+')*4.0d-2/heatunit'
  str_LQ=strarr(1,9)
  str_LQ[*,0]='    where(x(ixO^S,1) .lt. lH0)'
  str_LQ[*,1]='      lQgrid(ixO^S)=lQ0*lQgrid(ixO^S)'
  str_LQ[*,2]='    else where(x(ixO^S,1) .lt. xprobmax1/2.d0)'
  str_LQ[*,3]='      lQgrid(ixO^S)=lQ0*lQgrid(ixO^S)*dexp(-(x(ixO^S,1)-lH0)/lambda)'
  str_LQ[*,4]='    else where(x(ixO^S,1) .lt. xprobmax1-lH0)'
  str_LQ[*,5]='      lQgrid(ixO^S)=lQ1*lQgrid(ixO^S)*dexp(-(xprobmax1-lH0-x(ixO^S,1))/lambda)'
  str_LQ[*,6]='    else where'
  str_LQ[*,7]='      lQgrid(ixO^S)=lQ1*lQgrid(ixO^S)'
  str_LQ[*,8]='    end where'
  endif else begin
        str_define=strarr(1,1)
        str_define[*,0]='    double precision :: lQgrid(ixI^S),lambda,lQ0,lH0,tramp'
        str_LQ_val=strarr(1,1)
        str_LQ_val[*,0]='    lQ0=4.0d-2/heatunit'
        str_LQ=strarr(1,9)
        str_LQ[*,0]='    where(x(ixO^S,1) .lt. lH0)'
        str_LQ[*,1]='      lQgrid(ixO^S)=lQ0*lQgrid(ixO^S)'
        str_LQ[*,2]='    else where(x(ixO^S,1) .lt. xprobmax1/2.d0)'
        str_LQ[*,3]='      lQgrid(ixO^S)=lQ0*lQgrid(ixO^S)*dexp(-(x(ixO^S,1)-lH0)/lambda)'
        str_LQ[*,4]='    else where(x(ixO^S,1) .lt. xprobmax1-lH0)'
        str_LQ[*,5]='      lQgrid(ixO^S)=lQ0*lQgrid(ixO^S)*dexp(-(xprobmax1-lH0-x(ixO^S,1))/lambda)'
        str_LQ[*,6]='    else where'
        str_LQ[*,7]='      lQgrid(ixO^S)=lQ0*lQgrid(ixO^S)'
        str_LQ[*,8]='    end where'
    endelse
    str_mod_usr_in=[[str_mod_usr[*,0:pos_defined]],$
                    [str_define],$
                    [str_mod_usr[*,pos_defined+1:pos_LQ0]],$
                    [str_LQ_val],$
                    [str_mod_usr[*,pos_LQ0+1:pos_LQ]],$
                    [str_LQ],$
                    [str_mod_usr[*,pos_LQ+1:n_elements(str_mod_usr)-1]]]

    str_job = strarr(1,9)
    str_job[*,0] = '#!/bin/bash'
;    str_job[*,1] = '#BSUB -q mpi'
;    str_job[*,2] = '#BSUB -J P3D_500_MFR_'+str
;    str_job[*,3] = '#BSUB -n 24' 
;    str_job[*,4] = '#BSUB -W 7200'
;    str_job[*,5] = '#BSUB -o ./0100.out' 
;    str_job[*,6] = '#BSUB -e ./0100.err'
;    str_job[*,7] = 'cd /home/test/ltjupyter/work/1paperwork/20231108s3d_simulation/P3DS/mfrdata2/P3D/MFR_'+str
    str_job[*,7] = 'cd ./'
    str_job[*,8] = 'mpirun -np 24 ./amrvac'



openw,lun,fn00,/get_lun
printf,lun,str_amrvac_in
close,lun
free_lun,lun

openw,lun,fn10,/get_lun
printf,lun,str_mod_usr_in
close,lun
free_lun,lun

openw,lun,fn20,/get_lun
printf,lun,format='(3e25.16)',xyz_out
close,lun
free_lun,lun

openw,lun,fn30,/get_lun
printf,lun,str_job
close,lun
free_lun,lun

endfor
print,'Min: ',Q0M, Q1M
end
