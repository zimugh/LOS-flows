pro MFRs2paraview, MFR_path_in, geo_path_in = geo_path_in, plt_path_in = plt_path_in, out_path_in, Ngrids_in = Ngrids_in

 tlb=WIDGET_BASE(xsize=1200,ysize=900)
 DEVICE,get_screen_size=sz
 info=WIDGET_INFO(tlb,/GEOMETRY)
 tlb_XY=[info.SCR_XSIZE,info.SCR_YSIZE]
 offset=[sz-tlb_XY]/2
 WIDGET_CONTROL,tlb,XOFFSET=offset[0],YOFFSET=offset[1]
 WIDGET_CONTROL,tlb,map=0,/real
 prsbar=IDLITWDPROGRESSBAR(Group_LEADER=tlb,title='Progress',Cancel=cancelIn)

N_MFRs = n_elements(file_search(MFR_path_in+'/MFR*'))

for MFR_ID = 0, N_MFRs-1 do begin
     
  geo_path = geo_path_in+'/MFR_'+strtrim(string(MFR_ID),2)
  plt_path = plt_path_in+'/MFR_'+strtrim(string(MFR_ID),2)
  out_path_in2 = out_path_in+'/MFR_'+strtrim(string(MFR_ID),2)
  file_mkdir, out_path_in2
  
  convert_flag = convert2mfr(MFR_ID, geo_path = geo_path, plt_path = plt_path, out_path = out_path_in2, $
                      Ngrids = Ngrids_in)
  
  IF WIDGET_INFO(prsbar,/valid)THEN BEGIN
    IDLITWDPROGRESSBAR_SETVALUE,prsbar,((MFR_ID+1.)/N_MFRs)*100.
  ENDIF ELSE BEGIN
    tmp=DIALOG_MESSAGE('Cancel Percent'+STRING(((MFR_ID+1.)/N_MFRs)*100.)+'%',/info)
    Break
  ENDELSE
  WAIT,0.5
endfor

                    
end
