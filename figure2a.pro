;combination of :
;rain_pdf.pro
;event_lenght_pdf.pro
;====================================================================================
PRO tendency, filename1, pchf 
COMPILE_OPT IDL2
;====================================================================================
;spCAM5

nlon0=9
nlat0=6
filename1c = '/home/rtm/spCAM5/PCP/convective_rain_in_updrafts_and_downdrafts_TWP_Q_GT1E-4_h2_DEF3_1997.nc'

Id1c  = NCDF_OPEN(filename1c)
NCDF_VARGET, Id1c,   'CONV_RAIN',        pchf

pchf0     = FLTARR(7*4*12240)
pchf0_3h  = FLTARR(7*4*12240/12)
counter0  = 0.  

FOR l0=1,nlon0-2 DO BEGIN
 FOR l1=1,nlat0-2 DO BEGIN
  FOR l2=0,12239 DO BEGIN
   pchf0     [counter0] = pchf[l0,l1,l2]
   counter0=counter0+1.
  ENDFOR
 ENDFOR
ENDFOR

counter0 = 0.
FOR l0=0,n_elements(pchf0)-1, 12 DO BEGIN
   pchf0_3h   [counter0] = mean(pchf0[(l0):(l0+11)])
   counter0=counter0+1.
ENDFOR

;=====================================================================================
;GCM18

filename10a = '/home/rtm/GCM/CanAM4_hires/ensamble/pchf/TWP_sa_jcl_hf_1997_pchf_gs.nc' ; Convective Precipitation  kg m-2 s-1 

Id10a  = NCDF_OPEN(filename10a)

NCDF_VARGET, Id10a,   'PCHF',     pchf18a

nlon0=9
nlat0=4

pchf18 = (pchf18a)*3600. ; change units to kg/m2/h

ntim = 44160
pchf018     = FLTARR(nlon0*nlat0*(ntim))
pchf018_3h  = FLTARR(nlon0*nlat0*(ntim)/18)

counter018  = 0.

FOR l0=0,nlon0-1 DO BEGIN
 FOR l1=0,nlat0-1 DO BEGIN
  FOR l2=0,(ntim)-1 DO BEGIN
   pchf018     [counter018] = pchf18[l0,l1,l2]
   counter018 = counter018+1.
  ENDFOR
 ENDFOR
ENDFOR

counter018 = 0.
FOR l0=0,n_elements(pchf018)-1, 18 DO BEGIN
   pchf018_3h   [counter018] = mean(pchf018[(l0):(l0+17)])
   counter018=counter018+1.
ENDFOR


;===================================================================================
;CAM5 

filename20a = '/home/rtm/convection_paper/data/CAM5/PRECC_CAM5.nc' ; Convective Precipitation  m s-1 
Id20a  = NCDF_OPEN(filename20a)
NCDF_VARGET, Id20a,   'PRECC',     PRECC
NCDF_VARGET, Id20a,   'lon',       lon
NCDF_VARGET, Id20a,   'lat',       lat

PRECC = PRECC*3600*1000. ; change units to mm/h

nlon0=9
nlat0=6

ntim = 3900
cam5 = FLTARR(nlon0*nlat0*(ntim))
cam5[*] = 'NAN'

counter_cam  = 0.

FOR l0=60,68 DO BEGIN
 FOR l1=48,53 DO BEGIN
  FOR l2=0,(ntim)-1 DO BEGIN
   cam5 [counter_cam] = PRECC[l0,l1,l2]
   counter_cam = counter_cam+1.
  ENDFOR
 ENDFOR
ENDFOR

;===================================================================================
;TRMM 3B42 3-hours

filename30a = '/home/rtm/TRMM/TRMM_2deg.nc' ; Convective Precipitation  m s-1 
Id30a  = NCDF_OPEN(filename30a)
NCDF_VARGET, Id30a,   'PREC2',     PREC_TRMM

nlon_trmm = 10
nlat_trmm = 5
ntim_trmm = 3680
trmm = FLTARR(10*5*3680)
trmm[*] = 'NAN'

counter_trmm = 0.

FOR l0=0,9 DO BEGIN
 FOR l1=0,4 DO BEGIN
  FOR l2=0,(ntim_trmm)-1 DO BEGIN
   trmm [counter_trmm] = PREC_TRMM[l0,l1,l2]
   counter_trmm = counter_trmm+1.
  ENDFOR
 ENDFOR
ENDFOR

;===================================================================================
;generate own histogram function

;define min, max, binsize
xmin =  0.
xmax = 10.
bin0 =  0.1
bins = CEIL(xmax/bin0)

;define arrays
hist_arr_spcam_conv = FLTARR(bins)
hist_arr_gcm18_conv = FLTARR(bins)
hist_arr_cam5_conv  = FLTARR(bins)
hist_arr_trmm_conv  = FLTARR(bins)

  FOR bn0=0,bins-1 DO BEGIN
    INDEX1 = WHERE((pchf0_3h[*] GT (bn0*bin0 + xmin+0.01)) AND (pchf0_3h[*] LE ((bn0+1)*bin0 + xmin+0.01)), elem1)
    hist_arr_spcam_conv[bn0] = 100.*elem1/counter0 ; in percent from the total
  ENDFOR

  FOR bn0=0,bins-1 DO BEGIN
    INDEX1 = WHERE((pchf018_3h[*] GT (bn0*bin0 + xmin+0.01)) AND (pchf018_3h[*] LE ((bn0+1)*bin0 + xmin+0.01)), elem1)
    hist_arr_gcm18_conv[bn0] = 100.*elem1/counter018 ; in percent from the total
  ENDFOR

  FOR bn0=0,bins-1 DO BEGIN
    INDEX1 = WHERE((cam5[*] GT (bn0*bin0 + xmin+0.01)) AND (cam5[*] LE ((bn0+1)*bin0 + xmin+0.01)), elem1)
    hist_arr_cam5_conv[bn0] = 100.*elem1/counter_cam ; in percent from the total
  ENDFOR

  FOR bn0=0,bins-1 DO BEGIN
    INDEX1 = WHERE((trmm[*] GT (bn0*bin0 + xmin+0.01)) AND (trmm[*] LE ((bn0+1)*bin0 + xmin+0.01)), elem1)
    hist_arr_trmm_conv[bn0] = 100.*elem1/counter_trmm ; in percent from the total
  ENDFOR


X1 = xmin + bin0 * FINDGEN(bins)

;===================================================================================

rpsopen, 'figure1b_ensemble_3hourly.eps', /encap, xs=8, ys=8, /inches, /color

Device, Decomposed=0
LoadCT, 33, NColors=(100), Bottom=3
!P.MULTI = [0, 1, 1, 0, 0]

D=4.

a1=0.25 ;   x1 left column
a2=0.95 ;   x2 left column

b1=0.15 ;   row 1
b2=0.85 ;   

;================================================================
;panel (b)

plot, X1, hist_arr_spcam_conv, XTHICK=8.0/D, YTHICK=8.0/D,position=[a1, b1, a2, b2], $
    XRANGE = [0,3], YRANGE = [0.1,100.], /YLOG, XTICKINTERVAL=1, XMINOR=2, YTICKINTERVAL=4, YMINOR=1, $
    TITLE='(b)',$
    xtitle='Conv. PCP (mmh!E-1!N)',$
    ytitle='Frequency Density (%)', $
    charsize=10.5/D,charthick=20.5/D,thick=6.
oplot, X1, hist_arr_gcm18_conv, color=20, thick=6.
oplot, X1, hist_arr_cam5_conv, color=55, thick=6.
;oplot, X1, hist_arr_trmm_conv, color=80, thick=6.

rpsclose

END

