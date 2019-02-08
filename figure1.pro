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
filename1a = '/home/rtm/spCAM5/CAPE/reversible_CAPE_h2_1997_buoyancy+0K_and_rho_from_Tv.nc'
filename1c = '/home/rtm/spCAM5/PCP/convective_rain_in_updrafts_and_downdrafts_TWP_Q_GT1E-4_h2_DEF3_1997.nc'

Id1a  = NCDF_OPEN(filename1a)
Id1c  = NCDF_OPEN(filename1c)
NCDF_VARGET, Id1a,   'CAPE_dcapelse',    PCAPElse
NCDF_VARGET, Id1c,   'TOTAL_RAIN',       thf ; total pcp
NCDF_VARGET, Id1c,   'CONV_RAIN',        pchf

print, mean(pchf)/mean(thf)*100.
print, mean(thf)*24.

filename4 = '/home/rtm/spCAM5/GCM_OMEGA/GCM_OMEGA_TWP_start_from_top.nc'
Id4 = NCDF_OPEN(filename4)
NCDF_VARGET, Id4,   'GCM_OMEGA',   OMEGA

pchf0     = FLTARR(7*4*12240)
PCAPElse0 = FLTARR(7*4*12240)
MC        = FLTARR(7*4*12240)

counter0  = 0.  

FOR l0=1,nlon0-2 DO BEGIN
 FOR l1=1,nlat0-2 DO BEGIN
  FOR l2=0,12239 DO BEGIN
   pchf0     [counter0] = pchf[l0,l1,l2]
   ; in spCAM5 PCAPElse is CAPE generation
   PCAPElse0 [counter0] = PCAPElse[l0,l1,l2,0]*6. ; convert to J/kg/h, where 6 stands for 6 10-minute intervals. PCAPElse already is substracted from PCAPE
   MC        [counter0] = OMEGA[l0,l1,1,l2]*3600. ; levels start from the surface Pa/h
   counter0=counter0+1.
  ENDFOR
 ENDFOR
ENDFOR
;=====================================================================================
;GCM18

filename10a = '/home/rtm/GCM/CanAM4_hires/ensamble/pchf/TWP_sa_jcl_hf_1997_pchf_gs.nc' ; Convective Precipitation  kg m-2 s-1 
;filename10b = '/home/rtm/convection_paper/data/CanAM/TWP_sa_rtm_hf2_1997_plhf_gs.001.nc'
filename10h = '/home/rtm/GCM/CanAM4_hires/ensamble/whf/TWP_sa_jcl_hf_1997_whf_gs.nc'
filename10c = '/home/rtm/GCM/CanAM4_hires/ensamble/rCAPE/TWP_sa_jcl_hf_1997_rCAPE_gs.nc'
filename10d = '/home/rtm/GCM/CanAM4_hires/ensamble/rCAPElse/TWP_sa_jcl_hf_1997_rCAPElse_gs.nc'

Id10a  = NCDF_OPEN(filename10a)
;Id10b  = NCDF_OPEN(filename10b)
Id10c  = NCDF_OPEN(filename10c)
Id10d  = NCDF_OPEN(filename10d)
Id10h  = NCDF_OPEN(filename10h)

NCDF_VARGET, Id10a,   'PCHF',     pchf18a
;NCDF_VARGET, Id10b,   'PLHF',     plhf18
NCDF_VARGET, Id10c,   'CAPE',    PCAPE18
NCDF_VARGET, Id10d,   'rCAPElse', PCAPElse18a
NCDF_VARGET, Id10h,   'WHF',      whf18a

;print, mean(pchf18a)/mean(pchf18a+plhf18)*100.
;print, mean(pchf18a+plhf18)*86400.

whf18 = FLTARR(9, 4, 44160)
whf18[*,*,*]= 'NAN'
nlon0=9
nlat0=4

pchf18 = (pchf18a)*3600. ; change units to kg/m2/h
whf18 [*,*,*]= whf18a[*,*,47,*]*3600. ; change units to Pa/h ; level=47 is 1 level above the surface

; in CanAM4 cape generation rate is computed as the difference between PCAPElse and PCAPE
PCAPElse18 = (PCAPElse18a-PCAPE18)*4. ; convert to  J/kg/h, where 4 stands for 4 15-minute time intervals

ntim = 44160
pchf018     = FLTARR(nlon0*nlat0*(ntim))
PCAPElse018 = FLTARR(nlon0*nlat0*(ntim))
whf018        = FLTARR(nlon0*nlat0*(ntim))
counter018  = 0.

FOR l0=0,nlon0-1 DO BEGIN
 FOR l1=0,nlat0-1 DO BEGIN
  FOR l2=0,(ntim)-1 DO BEGIN
   pchf018     [counter018] = pchf18[l0,l1,l2]
   PCAPElse018 [counter018] = PCAPElse18[l0,l1,l2]
   whf018      [counter018] = whf18[l0,l1,l2]
   counter018 = counter018+1.
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
hist_arr_spcam_lse = FLTARR(bins)
hist_arr_gcm18_conv = FLTARR(bins)
hist_arr_gcm18_lse = FLTARR(bins)

hist_arr_spcam_mc  = FLTARR(bins)
hist_arr_gcm18_mc  = FLTARR(bins)


  FOR bn0=0,bins-1 DO BEGIN
    INDEX1 = WHERE((pchf0[*] GT (bn0*bin0 + xmin+0.01)) AND (pchf0[*] LE ((bn0+1)*bin0 + xmin+0.01)), elem1)
    hist_arr_spcam_conv[bn0] = 100.*elem1/counter0 ; in percent from the total
    hist_arr_spcam_lse [bn0] = mean(PCAPElse0[INDEX1], /NAN)
    hist_arr_spcam_mc  [bn0] = mean(MC[INDEX1], /NaN)
  ENDFOR

  FOR bn0=0,bins-1 DO BEGIN
    INDEX1 = WHERE((pchf018[*] GT (bn0*bin0 + xmin)) AND (pchf018[*] LE ((bn0+1)*bin0 + xmin)), elem1)
    hist_arr_gcm18_conv[bn0] = 100.*elem1/counter018 ; in percent from the total
    hist_arr_gcm18_lse [bn0] = mean(PCAPElse018[INDEX1])
    hist_arr_gcm18_mc  [bn0] = mean(whf018[INDEX1])
  ENDFOR

print, total(hist_arr_spcam_conv), total(hist_arr_gcm18_conv)
X1 = xmin + bin0 * FINDGEN(bins)

;===================================================================================
filename5a = '/home/rtm/spCAM5/event_lenght_1mm_1hr.nc'
Id5a  = NCDF_OPEN(filename5a)
NCDF_VARGET, Id5a,   'EL',   EL_spcam

filename5b = '/home/rtm/GCM/CanAM4_hires/event_lenght_GCM18_ensamble_1mm_1hr.nc'
Id5b  = NCDF_OPEN(filename5b)
NCDF_VARGET, Id5b,   'EL',   EL_gcm18

;define min, max, binsize
xmin =  0.
xmax = 12.
bin0 =  0.5
bins = CEIL(xmax/bin0)

;define arrays
hist_arr_spcam = FLTARR(bins)
hist_arr_gcm18 = FLTARR(bins)
TOTAL_spcam = N_ELEMENTS(WHERE((EL_spcam GT 0.) AND (EL_spcam LE 12.)))
TOTAL_gcm18 = N_ELEMENTS(WHERE((EL_gcm18 GT 0.) AND (EL_gcm18 LE 12.)))


  FOR bn0=0,bins-1 DO BEGIN
    INDEX1 = WHERE((EL_spcam[*] GT (bn0*bin0 + xmin)) AND (EL_spcam[*] LE ((bn0+1)*bin0 + xmin)), elem1)
    hist_arr_spcam[bn0] = 100.*elem1/TOTAL_spcam ; in percent from the total
  ENDFOR

  FOR bn0=0,bins-1 DO BEGIN
    INDEX2 = WHERE((EL_gcm18[*] GT (bn0*bin0 + xmin)) AND (EL_gcm18[*] LE ((bn0+1)*bin0 + xmin)), elem2)
    hist_arr_gcm18[bn0] = 100.*elem2/TOTAL_gcm18 ; in percent from the total
  ENDFOR

X2 = xmin + bin0 * FINDGEN(bins) + bin0/2.

tot_spCAM = 0.
tot_CanAM = 0.

FOR III=0,20 DO BEGIN
  tot_spCAM = tot_spCAM + hist_arr_spcam[III]
  tot_CanAM = tot_CanAM + hist_arr_gcm18[III]
  PRINT, III, tot_spCAM, tot_CanAM
ENDFOR
print, 'figure 1(c)', hist_arr_gcm18
;=====================================================================================
; Rain Event Definition
; see: new_event_definition_TWP_STRONG_FORCING_spCAM2000.pro
;=====================================================================================
EVENT_LENGHT = FLTARR(50000)
EVENT_LENGHT[*] = 'NaN'
EVENT_PCP = FLTARR(50000)
EVENT_PCP[*] = 'NaN'

count=0L

FOR nt=12,n_elements(pchf018[*])-49 DO BEGIN
  MAX_RAIN = 0.
; IF - detects time when convective rainfall starts and ends. Also it saves maximum convective rainfall
   IF (pchf018[nt-3] LT 0.1) AND (pchf018[nt-2] LT 0.1) AND (pchf018[nt-1] LT 0.1) AND (pchf018[nt] GT 0.1) THEN BEGIN; AND $
;      (pchf018[nt-6] LT 0.1) AND (pchf018[nt-5] LT 0.1) AND (pchf018[nt-4] LT 0.1) AND $
;      (pchf018[nt-9] LT 0.1) AND (pchf018[nt-8] LT 0.1) AND (pchf018[nt-7] LT 0.1) AND $
;      (pchf018[nt-12] LT 0.1) AND (pchf018[nt-11] LT 0.1) AND (pchf018[nt-10] LT 0.1) THEN BEGIN 
    FOR TIME2=nt+1,nt+48 DO BEGIN
     IF (pchf018[TIME2] GT MAX_RAIN) THEN BEGIN
      MAX_RAIN = pchf018[TIME2]
     ENDIF
     IF (pchf018[TIME2] LT 0.1) THEN GOTO, OUT
    ENDFOR ; TIME2
   ENDIF

   OUT:
;=====================
; IF - defines start and end times of an event that satisfies rainfall and dPCAPElse thresholds
; the thresholds are defined using 3-months (May,June,July)
;=====================
;STRONG FORCING STRONG RAINFALL
    IF (MAX_RAIN GT 1.0) THEN BEGIN
    EVENT_LENGHT[count] = (TIME2-nt)/4.
    EVENT_PCP[count] = MAX_RAIN
    START = nt
    FINISH= TIME2
    MAX_RAIN = 0.
    count=count+1
   ENDIF
ENDFOR
;==============================================================
;Define min and max on X, Y axis
X_MAX =  12.
X_MIN =   0.

;Define number of bins on X and Y axis
N_BINS = 24.

STEP_=0.
STEP_X=(X_MAX-X_MIN)/(N_BINS)
X3 = FLTARR(N_BINS)
PLOT_EVENT = FLTARR(N_BINS)
PLOT_EVENT[*] = 'NAN'

MIN_RAIN=1.
FOR XX=0,N_BINS-1 DO BEGIN
    INDEX = WHERE((EVENT_LENGHT[*] GT (X_MIN+(STEP_X*XX))) AND (EVENT_LENGHT[*] LE (X_MIN+(STEP_X*(XX+1)))) AND (EVENT_PCP[*] GT MIN_RAIN))
    IF (N_ELEMENTS(INDEX) GT 1.) THEN BEGIN
        PLOT_EVENT [XX] = MEAN(EVENT_PCP[INDEX], /NAN)
    ENDIF
ENDFOR

FOR XX=0,N_BINS-1 DO BEGIN
    X3[XX]=(((X_MIN+(STEP_X*XX)) + (X_MIN+(STEP_X*(XX+1))))/2.) ; mean interval value
ENDFOR

;================================================================
;Define Event, and collect data for all events

EVENT_LENGHT_spcam = FLTARR(10000)
EVENT_LENGHT_spcam[*] = 'NaN'
EVENT_PCP_spcam = FLTARR(10000)
EVENT_PCP_spcam[*] = 'NaN'

COUNT_EVENT = 0L

FOR TIME1=18, N_ELEMENTS(pchf0)-73 DO BEGIN
;=====================
        MAX_RAIN = 0.
;=====================
; DEF 2
        IF (pchf0[TIME1+1] GT 0.1) AND $
        (pchf0[TIME1] LT 0.1) AND (pchf0[TIME1-1] LT 0.1) AND (pchf0[TIME1-2] LT 0.1) AND $ ; THEN BEGIN
        (pchf0[TIME1-3] LT 0.1) AND (pchf0[TIME1-4] LT 0.1) AND (pchf0[TIME1-5] LT 0.1) THEN BEGIN
;        (pchf0[TIME1-6] LT 0.1) AND (pchf0[TIME1-7] LT 0.1) AND (pchf0[TIME1-8] LT 0.1) AND $
;        (pchf0[TIME1-9] LT 0.1) AND (pchf0[TIME1-10] LT 0.1) AND (pchf0[TIME1-11] LT 0.1) AND $
;        (pchf0[TIME1-12] LT 0.1) AND (pchf0[TIME1-13] LT 0.1) AND (pchf0[TIME1-14] LT 0.1) AND $
;        (pchf0[TIME1-15] LT 0.1) AND (pchf0[TIME1-16] LT 0.1) AND (pchf0[TIME1-17] LT 0.1) THEN BEGIN

                FOR TIME2=TIME1+1,TIME1+72 DO BEGIN
                        IF (pchf0[TIME2] GT MAX_RAIN) THEN BEGIN
                                MAX_RAIN = pchf0[TIME2]
                        ENDIF
;DEF 2
                        IF (pchf0[TIME2] LT 0.1) THEN GOTO, OUT1
                ENDFOR
        ENDIF
        OUT1:
;=====================
;DEF 2
        IF (MAX_RAIN GT 1.0) THEN BEGIN
;=====================
                EVENT_LENGHT_spcam[COUNT_EVENT] = (TIME2-TIME1)/6. ; convert to hours
                EVENT_PCP_spcam[COUNT_EVENT] = MAX_RAIN
                START = TIME1
                FINISH= TIME2
                MAX_RAIN = 0.
                COUNT_EVENT = COUNT_EVENT + 1
        ENDIF
;=====================
ENDFOR

;================================================================
;Define min and max on X, Y axis
X_MAX =  12.
X_MIN =   0.

;Define number of bins on X and Y axis
N_BINS = 24.

STEP_=0.
STEP_X=(X_MAX-X_MIN)/(N_BINS)
X4 = FLTARR(N_BINS)
PLOT_EVENT_spcam = FLTARR(N_BINS)
PLOT_EVENT_spcam[*] = 'NAN'

MIN_RAIN=1.
FOR XX=0,N_BINS-1 DO BEGIN
    INDEX = WHERE((EVENT_LENGHT_spcam[*] GT (X_MIN+(STEP_X*XX))) AND (EVENT_LENGHT_spcam[*] LE (X_MIN+(STEP_X*(XX+1)))) AND (EVENT_PCP_spcam[*] GT MIN_RAIN))
    IF (N_ELEMENTS(INDEX) GT 1.) THEN BEGIN
        PLOT_EVENT_spcam [XX] = MEAN(EVENT_PCP_spcam[INDEX])
    ENDIF
ENDFOR

FOR XX=0,N_BINS-1 DO BEGIN
    X4[XX]=(((X_MIN+(STEP_X*XX)) + (X_MIN+(STEP_X*(XX+1))))/2.) ; mean interval value
ENDFOR

;===================================================================================


rpsopen, 'figure1_ensemble.eps', /encap, xs=14, ys=12, /inches, /color

Device, Decomposed=0
LoadCT, 33, NColors=(100), Bottom=3
!P.MULTI = [0, 2, 2, 0, 0]

D=5.

a1=0.14 ;   x1 left column
a2=0.44 ;   x2 left column
a3=0.67
a4=0.97

b1=0.62 ;   row 1
b2=0.92 ;   

b3=0.14
b4=0.44

;b5=0.07
;b6=0.31

;b7=0.05
;b8=0.20

;================================================================
;panel (a)
print, hist_arr_spcam_lse

plot, X1, hist_arr_spcam_lse, XTHICK=8.0/D, YTHICK=8.0/D,position=[a1, b1, a2, b2], $
    XRANGE = [0,3], YRANGE = [-200.,400.], XTICKINTERVAL=1, XMINOR=2, YTICKINTERVAL=100, YMINOR=1, $
    TITLE='(a)', $
    xtitle='Conv. PCP (mmh!E-1!N)',$
    ytitle='dCAPE!DLSFT!N (Jkg!E-1!Nh!E-1!N)', $
    charsize=10.5/D,charthick=20.5/D,thick=8.
OPLOT, X1, hist_arr_gcm18_lse, color=20, thick=8. 

AXIS, YAXIS=1, YRANGE = [-50, 100.],  YTICKINTERVAL=25, YMINOR=1, YSTYLE = 1, $
   YTITLE = '-!4x!X (Pah!E-1!N)', $
charsize=10.5/D,charthick=20.5/D
OPLOT, X1, (-1.)*(hist_arr_spcam_mc*4.), linestyle=5, thick=8.
OPLOT, X1, (-1.)*(hist_arr_gcm18_mc*4.), linestyle=5, thick=8., color=20
;================================================================
;panel (b)

plot, X1, hist_arr_spcam_conv, XTHICK=8.0/D, YTHICK=8.0/D,position=[a3, b1, a4, b2], $
    XRANGE = [0,3], YRANGE = [0.1,100.], /YLOG, XTICKINTERVAL=1, XMINOR=2, YTICKINTERVAL=4, YMINOR=1, $
    TITLE='(b)',$
    xtitle='Conv. PCP (mmh!E-1!N)',$
    ytitle='Frequency Density (%)', $
    charsize=10.5/D,charthick=20.5/D,thick=5.
oplot, X1, hist_arr_gcm18_conv, color=20, thick=5.
print, hist_arr_spcam_conv[0:10]
print, hist_arr_gcm18_conv[0:10]

;================================================================
;panel (c)

plot, X2, hist_arr_spcam, XTHICK=8.0/D, YTHICK=8.0/D,position=[a1, b3, a2, b4], $
    XRANGE = [0,5], YRANGE = [0.1,50.], /YLOG, XTICKINTERVAL=1, XMINOR=2, YTICKINTERVAL=1, YMINOR=1, $
    TITLE='(c)',$
    xtitle='Event length (h)',$
    ytitle='Frequency Density (%)', $
    charsize=10.5/D,charthick=20.5/D,thick=8.
oplot, X2, hist_arr_gcm18, color=20, thick=8.

;================================================================
;panel (d)

plot, X4, PLOT_EVENT_spcam, XTHICK=8.0/D, YTHICK=8.0/D,position=[a3, b3, a4, b4], $
    XRANGE = [0,5], YRANGE = [1,3], XTICKINTERVAL=1, XMINOR=2, YTICKINTERVAL=1, YMINOR=1, $
    TITLE='(d)',$
    xtitle='Event length (h)',$
    ytitle='Mean Peak PCP (mmh!E-1!N)', $
    charsize=10.5/D,charthick=20.5/D,thick=8.
oplot, X3, PLOT_EVENT, color=20, thick=8.


rpsclose

END

