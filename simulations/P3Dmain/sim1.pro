
;fn='../../202308031640rbsl6/paraview1/rbsl_fl10c.vtk';;;;vtk4.2 maybe paraview5.1——5.5
fn='../datamr/pv54tdmfl200.vtk';;;;vtk4.2 maybe paraview5.1——5.5
;number=752635
;nlines=100

number=120312
nlines=200

Bfile='../datamr/test.dat'
outpath='../P3DS/mfrdata'

nx = 80
ny = 100
nz = 60

read_lines,fn,number,nlines,Bfile,outpath = outpath,nx=nx,ny=ny,nz=nz


amrvac_path='../READ/amrvac_doc'
path_MFR=outpath+'/MFR_lines/location_MFRs'
print,path_MFR

;make_AMRVACs,outpath = outpath, amrvac_path = amrvac_path, path_MFR = path_MFR,$
;                 ASY_Heat=ASY_Heat,AB2=AB2,Bfile = Bfile, nx = nx, ny = ny, nz = nz
                 
make_AMRVACs,outpath = outpath, amrvac_path = amrvac_path, path_MFR = path_MFR

end
                 
