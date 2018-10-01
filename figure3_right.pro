;same as CanAM4_rtm_1.pro
;only in this code use hi-res gcm18 runs (see email from Jason on April 2, 2017)
;====================================================================================
PRO tendency, filename1, filename2, ctl_tnt 
COMPILE_OPT IDL2
;====================================================================================
;omega from level
lvl_w=47

filename1a = '/home/rtm/GCM/CanAM4_hires/TWP_sa_rtm_hf2_1997_pchf_gs.001.nc' ; Convective Precipitation	kg m-2 s-1 
filename1c = '/home/rtm/GCM/CanAM4_hires/TWP_sa_rtm_hf2_1997_rCAPE_gs.001.nc' ;rCAPE in J/kg
filename1d = '/home/rtm/GCM/CanAM4_hires/TWP_sa_rtm_hf2_1997_rCAPElse_gs.001.nc' ;rCAPElse 
filename1h = '/home/rtm/GCM/CanAM4_hires/TWP_sa_rtm_hf2_1997_whf_gs.001.nc'
filename1j = '/home/rtm/GCM/CanAM4_hires/TWP_sa_rtm_hf2_1997_RH_gs.001.nc'

Id1a  = NCDF_OPEN(filename1a)
Id1c  = NCDF_OPEN(filename1c)
Id1d  = NCDF_OPEN(filename1d)
Id1h  = NCDF_OPEN(filename1h)
Id1j  = NCDF_OPEN(filename1j)

NCDF_VARGET, Id1a,   'PCHF',     pchf
NCDF_VARGET, Id1c,   'CAPE',    PCAPE
NCDF_VARGET, Id1d,   'CAPElse', PCAPElse
NCDF_VARGET, Id1h,   'WHF',      whf
NCDF_VARGET, Id1j,   'RH',       RH

pchf = pchf*3600. ; convert to kg/m2/h
whf = whf*3600.   ; convert to Pa/h

PCAPElse = (PCAPElse-PCAPE)*4. ; 4* to change units to J/kg/h

;=====================================================================================
; Rain Event Definition
; see: new_event_definition_TWP_STRONG_FORCING_spCAM2000.pro
;=====================================================================================

arr_prc=FLTARR(10000,72)
arr_prc[*,*]='NaN'
arr_CAPE =FLTARR(10000,72)
arr_CAPE [*,*]='NaN'
arr_CIN =FLTARR(10000,72)
arr_CIN [*,*]='NaN'
arr_CAPElse =FLTARR(10000,72)
arr_CAPElse [*,*]='NaN'

arr_whf =FLTARR(10000,72)
arr_whf [*,*]='NaN' 
arr_RH =FLTARR(10000,72,nlev0)
arr_RH [*,*,*]='NaN'


arr_dCAPE = FLTARR(10000,72)
arr_dCAPE[*,*] = 'NaN'

count=0
FOR NLON=0,nlon0-1 DO BEGIN
 FOR NLAT=0,nlat0-1 DO BEGIN
  FOR nt=12,n_elements(pchf[0,0,*])-49 DO BEGIN
  MAX_RAIN = 0.
  MAX_LSE = 0.
; IF - detects time when convective rainfall starts and ends. Also it saves maximum convective rainfall
   IF (pchf[NLON,NLAT,nt-3] LT 0.1) AND (pchf[NLON,NLAT,nt-2] LT 0.1) AND (pchf[NLON,NLAT,nt-1] LT 0.1) AND (pchf[NLON,NLAT,nt] GT 0.1) AND $
      (pchf[NLON,NLAT,nt-6] LT 0.1) AND (pchf[NLON,NLAT,nt-5] LT 0.1) AND (pchf[NLON,NLAT,nt-4] LT 0.1) AND $
      (pchf[NLON,NLAT,nt-9] LT 0.1) AND (pchf[NLON,NLAT,nt-8] LT 0.1) AND (pchf[NLON,NLAT,nt-7] LT 0.1) AND $
      (pchf[NLON,NLAT,nt-12] LT 0.1) AND (pchf[NLON,NLAT,nt-11] LT 0.1) AND (pchf[NLON,NLAT,nt-10] LT 0.1) THEN BEGIN 
    FOR TIME2=nt+1,nt+48 DO BEGIN
     IF (pchf[NLON,NLAT,TIME2] GT MAX_RAIN) THEN BEGIN
      MAX_RAIN = pchf[NLON,NLAT,TIME2]
      MAX_LSE = PCAPElse[NLON,NLAT,TIME2]
     ENDIF
     IF (pchf[NLON,NLAT,TIME2] EQ 0.) THEN GOTO, OUT
    ENDFOR ; TIME2
   ENDIF

   OUT:
;=====================
; IF - defines start and end times of an event that satisfies rainfall and dPCAPElse thresholds
; the thresholds are defined using 3-months (May,June,July)
;=====================
;STRONG FORCING STRONG RAINFALL
;   IF (MAX_RAIN GT 1.0) AND (MAX_LSE LT 50.) THEN BEGIN
    IF (MAX_RAIN GT 1.0) THEN BEGIN
    START = nt 
    FINISH= TIME2
    FOR TIME3=START-12,FINISH DO BEGIN ; from -3 hours to the end of an event time
     arr_prc [count, TIME3-START+12] = pchf[NLON,NLAT,TIME3]
     arr_CAPE[count, TIME3-START+12] = PCAPE[NLON,NLAT,TIME3]
     arr_CIN [count, TIME3-START+12] = PCIN[NLON,NLAT,TIME3]
     arr_CAPElse[count, TIME3-START+12] = PCAPElse[NLON,NLAT,TIME3] ; the forcing at t=0 is for t=0 to t=1 ?
     arr_whf [count, TIME3-START+12] =  whf[NLON,NLAT,lvl_w,TIME3]
     arr_dCAPE[count, TIME3-START+12] = (PCAPE[NLON,NLAT,TIME3] - PCAPE[NLON,NLAT,TIME3-1])*4.
     arr_RH[count, TIME3-START+12,*] = RH[NLON,NLAT,*,TIME3]
    ENDFOR
    MAX_RAIN = 0.
    MAX_LSE = 0.
    count=count+1
   ENDIF
  ENDFOR
 ENDFOR
ENDFOR
print, count
;================================================================
;this part is for -3hours to 1-time step prior to initiation time

mean_arr_prc0=FLTARR(12)
mean_arr_prc0[*]='NaN'
mean_arr_whf0 =FLTARR(12)
mean_arr_whf0 [*]='NaN'
mean_arr_CAPE0 =FLTARR(12)
mean_arr_CAPE0 [*]='NaN'
mean_arr_CIN0 =FLTARR(12)
mean_arr_CIN0 [*]='NaN'
mean_arr_CAPElse0 =FLTARR(12)
mean_arr_CAPElse0 [*]='NaN'

mean_arr_dCAPE0 =FLTARR(12)
mean_arr_dCAPE0 [*]='NaN'
mean_arr_RH0=FLTARR(nlev0,12)
mean_arr_RH0[*,*]='NaN'


nev=20
FOR nt1 = 0,11 DO BEGIN
 INDEX_prc = WHERE(FINITE(arr_prc [*,nt1]))
 IF (n_elements(INDEX_prc) GT nev) THEN BEGIN
  mean_arr_prc0 [nt1] = mean(arr_prc [*,nt1], /NaN)
 ENDIF
 IF (n_elements(INDEX_prc) GT nev) THEN BEGIN
  mean_arr_whf0[nt1] = mean(arr_whf[*,nt1], /NaN)
 ENDIF
 IF (n_elements(INDEX_prc) GT nev) THEN BEGIN
  mean_arr_CAPE0[nt1] = mean(arr_CAPE[*,nt1], /NaN)
  mean_arr_CIN0[nt1] = mean(arr_CIN[*,nt1], /NaN)
 ENDIF
 IF (n_elements(INDEX_prc) GT nev) THEN BEGIN
  mean_arr_CAPElse0[nt1] = mean(arr_CAPElse[*,nt1], /NaN)
 ENDIF
 IF (n_elements(INDEX_prc) GT nev) THEN BEGIN
  mean_arr_dCAPE0[nt1] = mean(arr_dCAPE[*,nt1], /NaN)
 ENDIF
 IF (n_elements(INDEX_prc) GT nev) THEN BEGIN
   FOR NL0=0,nlev0-1 DO BEGIN
    mean_arr_RH0 [NL0,nt1] = mean(arr_RH [*,nt1,NL0], /NaN)
   ENDFOR
 ENDIF
ENDFOR
;=============================================================================

nmbr = 20 ; number of elements in X grid, to which rain event times are interpolated to. 
arr_prc_new = FLTARR(3000, nmbr+1)
arr_prc_new[*,*] = 'NaN'

arr_whf_new = FLTARR(3000, nmbr+1)
arr_whf_new[*,*] = 'NaN'

arr_CAPE_new = FLTARR(3000, nmbr+1)
arr_CAPE_new[*,*] = 'NaN'

arr_CIN_new = FLTARR(3000, nmbr+1)
arr_CIN_new[*,*] = 'NaN'

arr_dCAPE_new = FLTARR(3000, nmbr+1)
arr_dCAPE_new[*,*] = 'NaN'

arr_CAPElse_new = FLTARR(3000, nmbr+1)
arr_CAPElse_new[*,*] = 'NaN'

arr_RH_new = FLTARR(3000, nlev0, nmbr+1)
arr_RH_new[*,*,*] = 'NaN'

FOR newtime1=0,COUNT-1 DO BEGIN
 end_event=0
 FOR newtime2=12,70 DO BEGIN
   IF (arr_prc[newtime1, newtime2] LT 0.01) AND (end_event EQ 0) THEN BEGIN
     end_event=newtime2
     differ = INDGEN(newtime2 - 11 + 1)
     newtime= INDGEN(nmbr+1)*((n_elements(differ)-1)/FLOAT(nmbr))

     ;use INTERPOL
     arr_prc_new     [newtime1, *] = INTERPOL(arr_prc    [newtime1, 11:newtime2], differ, newtime)
     arr_whf_new     [newtime1, *] = INTERPOL(arr_whf    [newtime1, 11:newtime2], differ, newtime)
     arr_CAPE_new    [newtime1, *] = INTERPOL(arr_CAPE   [newtime1, 11:newtime2], differ, newtime)
     arr_CIN_new     [newtime1, *] = INTERPOL(arr_CIN    [newtime1, 11:newtime2], differ, newtime)
     arr_CAPElse_new [newtime1, *] = INTERPOL(arr_CAPElse[newtime1, 11:newtime2], differ, newtime)
     arr_dCAPE_new   [newtime1, *] = INTERPOL(arr_dCAPE  [newtime1, 11:newtime2], differ, newtime)
     FOR NL0=0,nlev0-1 DO BEGIN
       arr_RH_new    [newtime1, NL0, *] = INTERPOL(arr_RH [newtime1, 11:newtime2, NL0], differ, newtime)
     ENDFOR
   ENDIF
 ENDFOR
ENDFOR

;================================================================
mean_arr_prc1=FLTARR(nmbr+1)
mean_arr_prc1[*]='NaN'
mean_arr_whf1 =FLTARR(nmbr+1)
mean_arr_whf1 [*]='NaN'
mean_arr_CAPE1 =FLTARR(nmbr+1)
mean_arr_CAPE1 [*]='NaN'
mean_arr_CIN1 =FLTARR(nmbr+1)
mean_arr_CIN1 [*]='NaN'
mean_arr_dCAPE1 =FLTARR(nmbr+1)
mean_arr_dCAPE1 [*]='NaN'
mean_arr_CAPElse1 =FLTARR(nmbr+1)
mean_arr_CAPElse1 [*]='NaN'

mean_arr_RH1=FLTARR(nlev0,nmbr+1)
mean_arr_RH1[*,*]=1.E30

nev=count/10

FOR nt1=0,nmbr DO BEGIN
 INDEX_prc = WHERE(FINITE(arr_prc_new [*,nt1]))
 IF (n_elements(INDEX_prc) GT nev) THEN BEGIN
  mean_arr_prc1 [nt1] = mean(arr_prc_new [*,nt1], /NaN)
 ENDIF
 IF (n_elements(INDEX_prc) GT nev) THEN BEGIN
  mean_arr_whf1[nt1] = mean(arr_whf_new[*,nt1], /NaN)
 ENDIF
 IF (n_elements(INDEX_prc) GT nev) THEN BEGIN
  mean_arr_CAPE1[nt1] = mean(arr_CAPE_new[*,nt1], /NaN)
  mean_arr_CIN1[nt1] = mean(arr_CIN_new[*,nt1], /NaN)
 ENDIF
 IF (n_elements(INDEX_prc) GT nev) THEN BEGIN
  mean_arr_CAPElse1[nt1] = mean(arr_CAPElse_new[*,nt1], /NaN)
 ENDIF
 IF (n_elements(INDEX_prc) GT nev) THEN BEGIN
  mean_arr_dCAPE1[nt1] = mean(arr_dCAPE_new[*,nt1], /NaN)
 ENDIF
 IF (n_elements(INDEX_prc) GT nev) THEN BEGIN
   FOR NL0=0,nlev0-1 DO BEGIN
    mean_arr_RH1 [NL0,nt1] = mean(arr_RH_new [*,NL0,nt1], /NaN)
   ENDFOR
 ENDIF
ENDFOR
;================================================================
;connect two arrays, array1 prior to t=0, and array 2 after t=0

mean_arr_prc    = FLTARR(12+nmbr+1)
mean_arr_whf    = FLTARR(12+nmbr+1)
mean_arr_CAPE   = FLTARR(12+nmbr+1)
mean_arr_CIN    = FLTARR(12+nmbr+1)
mean_arr_CAPElse= FLTARR(12+nmbr+1)
mean_arr_dCAPE  = FLTARR(12+nmbr+1)
mean_arr_RH     = FLTARR(nlev0,12+nmbr+1)

FOR ntim1=0,11 DO BEGIN
  mean_arr_prc     [ntim1] = mean_arr_prc0     [ntim1]
  mean_arr_whf     [ntim1] = mean_arr_whf0    [ntim1]
  mean_arr_CAPE    [ntim1] = mean_arr_CAPE0    [ntim1]
  mean_arr_CIN     [ntim1] = mean_arr_CIN0    [ntim1]
  mean_arr_CAPElse [ntim1] = mean_arr_CAPElse0 [ntim1]
  mean_arr_dCAPE   [ntim1] = mean_arr_dCAPE0   [ntim1]
  mean_arr_RH      [*,ntim1] = mean_arr_RH0    [*,ntim1]
ENDFOR

FOR ntim2=0,nmbr DO BEGIN
  mean_arr_prc     [ntim2+12] = mean_arr_prc1     [ntim2]
  mean_arr_whf     [ntim2+12] = mean_arr_whf1     [ntim2]
  mean_arr_CAPE    [ntim2+12] = mean_arr_CAPE1    [ntim2]
  mean_arr_CIN     [ntim2+12] = mean_arr_CIN1    [ntim2]
  mean_arr_CAPElse [ntim2+12] = mean_arr_CAPElse1 [ntim2]
  mean_arr_dCAPE   [ntim2+12] = mean_arr_dCAPE1   [ntim2]
  mean_arr_RH      [*,ntim2+12] = mean_arr_RH1    [*,ntim2]
ENDFOR

;====================================================================================
;15=DARK BLUE
;25=MEDIUM BLUE
;35=LIGHT BLUE
;45=MAGENTA
;55=GREEN
;65=YELLOW
;70=ORANGE
;85=LIGHT RED
;95=DARK RED

Device, Decomposed=0
LOADCT, 33, NColors=100, Bottom=3

!P.MULTI = [0, 4, 1, 0, 0]

a1=0.20 ;   x1 left column
a2=0.84 ;   x2 left column

b1=0.77 ;   row 1
b2=0.97 ;   

c1=0.53 ;   row 2
c2=0.73 ; 

d1=0.29
d2=0.49

e1=0.05
e2=0.25

x1=[0,0]
z1=[0,3]
z2=[-1000,5000]

IF (REG EQ 1) THEN BEGIN
rpsopen, 'GCM18_hires_rtm_1_TWP_all_forcing_1mm_time_rCAPE.eps', /encap, xs=10, ys=30, /inches, /color
ENDIF

D=2.2

time = INDGEN(12+nmbr+1)-12.
plot, time , mean_arr_prc, THICK=8., XTHICK=8.0, YTHICK=8.0, position=[a1, b1, a2, b2], $
    XRANGE = [-10,nmbr], YRANGE = [0,1.5], XTICKINTERVAL=nmbr/4., XMINOR=2, YTICKINTERVAL=0.5, YMINOR=1, $
;    xtitle='time',$
    ytitle='Conv PCP (mm h!E-1!N)',$
    title='CanAM4.3', $
    charsize=12.5/D,charthick=20.5/D;,thick=5.0
oplot, x1, z1

AXIS, YAXIS=1, YRANGE = [-150.,0.],YTICKINTERVAL=50,COLOR=20, YMINOR=1, YSTYLE = 1, $
   YTITLE = '!4x!X (Pah!E-1!N)', charsize=12.5/D,charthick=20.5/D
OPLOT, time, mean_arr_whf*(1.5/150.)+1.5, color=20, thick=8.

plot, time , mean_arr_CAPE, THICK=8., XTHICK=8.0, YTHICK=8.0, position=[a1, c1, a2, c2], $
    XRANGE = [-10,nmbr], YRANGE = [300,500], XTICKINTERVAL=nmbr/4., XMINOR=2, YTICKINTERVAL=100, YMINOR=1, $
;    xtitle='time',$
    ytitle='CAPE (Jkg!E-1!N)',$
;    subtitle='',title='Strong Forcing', $
    charsize=12.5/D,charthick=20.5/D;,thick=5.0
oplot, x1, z2

AXIS, YAXIS=1, YRANGE = [-15.,0.],YTICKINTERVAL=5,COLOR=20, YMINOR=1, YSTYLE = 1, $
   YTITLE = 'CIN (Jkg!E-1!N)', charsize=12.5/D,charthick=20.5/D
OPLOT, time, mean_arr_CIN*(200./15.)+500., color=20, thick=8.
;oplot, time, mean_arr_dCAPE, COLOR=90, thick=8.

;oplot, time, (mean_arr_CAPElse), COLOR=20 
;oplot, time, (mean_arr_CAPElsp), COLOR=90

plot, time , mean_arr_CAPElse, THICK=8., XTHICK=8.0, YTHICK=8.0, position=[a1, d1, a2, d2], $
    XRANGE = [-10,nmbr], YRANGE = [0,200], XTICKINTERVAL=nmbr/4., XMINOR=2, YTICKINTERVAL=100, YMINOR=1, $
;    xtitle='time',$
    ytitle='dCAPE!DLSFT!N (Jkg!E-1!Nh!E-1!N)',$
;    subtitle='',title='Strong Forcing', $
    charsize=12.5/D,charthick=20.5/D;,thick=5.0
oplot, x1, z2

GCM_P=[1.09393,1.45874,1.92544,2.53210,3.31254,4.31458,5.57166,7.17193,9.15651,11.1419,15.1906,18.2285,23.2956,28.3602,35.4417,43.5411,52.6482,64.7940, $
       77.9558,93.1324,111.341,132.590,155.858,183.161,213.509,247.882,286.300,328.752,374.231,424.755,479.306,536.885,598.488,655.038,705.521,749.932, $
       788.287,821.604,849.858,874.076,894.261,911.416,925.542,939.668,953.794,967.920,982.048,996.172,1005.25]

;print, min(mean_arr_RH), max(mean_arr_RH)
max_value = 100.
min_value = 60.
levels = 20.  ; Number of contours
LoadCT, 33, NColors=levels, Bottom=3
step = (max_value - min_value) / levels
userLevels = IndGen(levels) * step + min_value
contour,(transpose(mean_arr_RH[*,*])), time,GCM_P, /fill, XRANGE = [-10,nmbr], YRANGE = [1000,500], C_Colors=(Indgen(levels)+3), Background=255, $
LEVELS=userLevels, $
   XSTYLE = 1, YSTYLE = 1, XTHICK=12.5/D, YTHICK=12.5/D, position=[a1, e1, a2, e2], Color=2, $
   CHARSIZE=12.5/D,CHARTHICK=20.5/D,$
;   TITLE = ' RH', $
   XTITLE = 'Time',$
   YTITLE = 'Height (hPa)'
;
Contour, (transpose(mean_arr_RH[*,*])), time, GCM_P, /Overplot, Levels=userLevels, C_THICK=12.5/D, C_CHARSIZE=10.5/D, C_CHARTHICK=12.5/D, /Follow, color=2
;=========================================

rpsclose

END

