Home_path = '/data/home/chenpf/IDL/radiation'
.compile -v '/data/home/chenpf/IDL/radiation/aia_resp.pro'
;.compile -v Home_path + '/aia_resp.pro'
file_xdr = Home_path + '/aia_resp.xdr'
restore , file = Home_path + '/grid3_o.sav', /ver
.compile -v '/data/home/chenpf/IDL/P3D_codes/grid3_test.pro'
.compile -v '/data/home/chenpf/IDL/P3D_codes/grid3_emission_test.pro'
.compile -v '/data/home/chenpf/IDL/radiation/gird3_emission_opth_test.pro'


plot_image,img3,xst=4+1,yst=4+1,pos=pos
plot,sx,sy,xr=sx,yr=sy,pos=pos,/nodata,/noerase,color=0,tit=tit3,thick=18,xtit=xtxt,ytit=ytxt


MFR_path = '/data/home/chenpf/nyw/sav'
out_path = '/data/home/chenpf/nyw/P30_render'
file_xdr = '/data/home/chenpf/IDL/radiation/aia_resp.xdr'
.compile -v '/data/home/chenpf/IDL/radiation/aia_resp.pro'
.compile -v '/data/home/chenpf/IDL/radiation/p3d_render.pro'
P3D_render, MFR_path = MFR_path, out_path = out_path, file_xdr = file_xdr

