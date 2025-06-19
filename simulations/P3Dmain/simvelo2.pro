

plt_path_in='/home/test/ltjupyter/work/1paperwork/20231108s3d_simulation/P3DS/mfrdata2/P3D/data2para'
out_path_in='/home/test/ltjupyter/work/1paperwork/20231108s3d_simulation/P3DS/mfrdata2/P3D/data2para'

N_MFRs = n_elements(file_search(plt_path_in+'/MFR_1h*'))

Ngrids=10000

for MFR_ID = 0, N_MFRs-1 do begin

	out_path = out_path_in+'/MFR_1h'+strtrim(string(MFR_ID),2)

	plt_path = plt_path_in+'/MFR_1h'+strtrim(string(MFR_ID),2)

	datav1 = file_search(plt_path + '/sav/v1.sav')
	dataxyz = file_search(plt_path + '/sav/xyz.sav')
	restore, datav1
	restore, dataxyz

	help,datav1
	help,dataxyz
	
	print,size(v1)
	print,size(xyz_o)
	Ngrids_ID = indgen(Ngrids)
		for it = 0, (size(v1))[1]-1 do begin
		  openw, 1, out_path + '/'+'v1'+strtrim(string(MFR_ID),2)+'mfr_path'+strtrim(string(it),2) + '.vtk'
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
		  printf, 1, 'SCALARS'+' '+'v1'+' '+'double'+' '+'1'
		  printf, 1, 'LOOKUP_TABLE  table1'
		  printf, 1, format='(1f)', reform(v1[it,*],Ngrids)
		  
		  close,1
		  
		endfor
		
	endfor

end


                 
