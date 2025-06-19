function aia_resp, logT, channel, all=all, load=load, file=file, logN=logN
;+
; Calculates response of the AIA channels for given temperature cube.
;
; Uses precalculated response tables interpolated from aia response
; function provided through SSW, aia_get_response(/temp,/dn).
;
; If density values are provided it returns the emissivity at
; each gridpoint.
;
;
; INPUTS:  logT      values of the temperature in logT
;                    can be single number, or 1, 2, or 3D array
;
;          channel   AIA channel(s) for which response is to be calculated
;                    Channel can be single number or 1D array.
;                    It can be either the index of the AIA channel (0...6)
;                    or the wavelength of the channel (94,131,171,...)
;
; KEYWORDS:  /all    if set do all AIA channels
;
;            /load   if set, (reload) the response table
;                    (only needeed if table changed...)
;
;            file    file name of rsponse table (optional)
;                    default is !HOME+'/idl/aia/response/aia_resp.xdr'
;
;            logN    If set returns emissivity = response*density^2
;                    IMPORTANT:  logN has to have the SAME dimension as logT
;                                [N]=cm-3
;                    In this case the units returned are DN s^-1 pix^-1 Mm^-1
;
; RETURNED VALUE:    structure with three entries {r, ch, units}:
;                    .r:  has same dimension as input log, extended by
;                         the dimension for the channels.
;                         e.g. logT=fltarr(100,150,120)   --> returns array
;                              channel=[171,304]          --> [100,150,120,2]
;                    .ch:    array of channels
;                    .units: physical units
;
; Units:  as in the original AIA response files the response is
;         returned in units of DN cm^5 s^-1 pix^-1
;
;         If logN is provided it returns the emissivity [DN s^-1 pix^-1 Mm^-1].
;         In this case one has to integrate through the 3D cube to get
;         the image as AIA would see it (using Mm for the unit along the LOS).
;
;
; USAGE:  resp=aia_resp(logT, [171,304]) -> response in DN cm^5 s^-1 pix^-1
;         resp=aia_resp(logT, [2,5])     -> same (channes indexed)
; 
;         resp=aia_resp(logT, [2,5], logN) -> emissivity in DN s^-1 pix^-1 Mm^-1
;
;         resp = aia_resp(logt, [304,171,193,211], logN=logN) 
;
;
; SIDE EFFECTS: -- Needs aia_resp_prep to calculate response tables.
;               -- Response table stored by default in
;                  !HOME+'/idl/aia/response/aia_resp.xdr'
;               -- uses common block for convenient storage of response table
;
;
; WRITTEN:  v1:  peter@mps.mpg.de  9.Feb.2011
;-


Common aia_resp, ri, loaded

if not keyword_set(file) then file=!HOME+'/idl/aia/response/aia_resp.xdr'

if keyword_set(load) then loaded=0

if not keyword_set(loaded) then $
   begin
   print, '% loading response file'
   restore, file
   end


;
;-----------------------------------------------------------------------------
; get the channel:
; either index (0-6) or the wavelength (94,131,171,...)
; channel can just one numer or an array.
; get the indices of the channels, ich.
;-----------------------------------------------------------------------------
;
if keyword_set(all) then channel=indgen(n_elements(ri.channels))

n_channels=n_elements(channel)
;
if max(channel) gt 7 then $
   begin
   ich=lonarr(n_channels)
   for i=0, n_channels-1 do $
       ich[i]=where(ri.channels eq 'A'+strtrim(fix(channel[i]),2) )
   end $
else $
   ich=fix(channel)
;
;-----------------------------------------------------------------------------
; get indices in the table based on temperature values
;-----------------------------------------------------------------------------
;
idx = round((logT-ri.Trange[0])/(ri.Trange[1]-ri.Trange[0])*(ri.num-1))

;
; should the input temperature be outside the range in the table
; set the index to zero.
; for this the response at the top of the response table (at entry
; width index ri.num instead of the nominammlz highest entry ri.num-1)
; a zero is present.
;
bad = where(idx lt 0 or idx ge ri.num)
if bad[0] ge 0 then idx[bad]=ri.num
;
;-----------------------------------------------------------------------------
; initialze the response function
; Take care of the possible case of the dimension of logT
;-----------------------------------------------------------------------------
;
sz=size(logT)
case sz[0] of
0: resp=fltarr(1,                 n_channels,/nozero)
1: resp=fltarr(sz[1],             n_channels,/nozero)
2: resp=fltarr(sz[1],sz[2],       n_channels,/nozero)
3: resp=fltarr(sz[1],sz[2],sz[3], n_channels,/nozero)
endcase


;
; if logN is set calculate emissivity: response*n^2:
; multiply by N^2, otherwise multiply by 1...
;
if not keyword_set(logN) then logN=0.

for i=0, n_channels-1 do $
    begin
    case sz[0] of
    0: resp          = (ri.resp[*,ich[i]])[idx]*10^(2*logN)  ; ri.logt[idx]
    1: resp[*,i]     = (ri.resp[*,ich[i]])[idx]*10^(2*logN)
    2: resp[*,*,i]   = (ri.resp[*,ich[i]])[idx]*10^(2*logN)
    3: resp[*,*,*,i] = (ri.resp[*,ich[i]])[idx]*10^(2*logN)
    endcase
    end
;
; if calculing emissivity get unit per Mm (instead of per cm)
;
if keyword_set(logN) then $
   begin
   resp=resp*1e8
   units='DN s^-1 pix^-1 Mm^-1'
   end $
else $
   units=ri.units

ch=ri.channels[ich]
lft=ri.lft[ich]
for i=0, n_channels-1 do ch[i]=strmid(ch[i],1,strlen(ch[i])-1)

return, {r:resp, ch:ch, lft:lft, units:units }
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;; dim=1

;; case dim of
;; 0: t= randomu(seed, 1)
;; 1: t= randomu(seed, 10000)
;; 2: t= randomu(seed, 100,100)
;; 3: t= randomu(seed, 256,256,256)
;; endcase

;; ;logN=randomu(seed,100,100,100)*4+7

;; logt=t*5+4

;; channel=[304,171,193,211]

;; date
;; resp = aia_resp(logt, channel, logN=9.0) 
;; date

;; norm=(1e9)^2*1e8  ; assuming density 1e9 cm-3  and  LOS length of 1 Mm

;; usrcol

;; plot, logt,resp[*,1],ps=3,/yl, xr=[4.8,7.5], xs=1

;; restore,!HOME+'/idl/aia/response/aia_resp_prep.xdr'
;; oplot, r.logte, r.a304.tresp*norm, psym=-4, col=11
;; oplot, r.logte, r.a171.tresp*norm, psym=-4, col=12
;; oplot, r.logte, r.a193.tresp*norm, psym=-4, col=13
;; oplot, r.logte, r.a211.tresp*norm, psym=-4, col=14

;; wshow,ic=0

;; end
