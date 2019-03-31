;if you use 3-hr mean spcam5 fields 
;change convert=0 to convert=1

;start with zero rain and identify when rain is initiated and excides 
;threshold value (1.0 mm/h) 

;same as new_event_definition_TWP_STRONG_FORCING_spCAM1997.pro
;only instead hours in the plotting use no units from 0 (initiation time) to 1 (event ends)
;that way all events (short and long) will start at 0 and end at 1
;===============================================================
PRO read_processed_ncdf, filename1, filename3, CAPE_capefld, CAPE_dcapelse, CAPE_dcapelsp, CAPE_dcapels, CAPE_dcapecue, CAPE_dcapecup, CAPE_dcapecu, CONV_RAIN
COMPILE_OPT IDL2
;===============================================================
nlon=9
nlat=6
nlev=30
count1 = 85*144

filename1 = '/home/rtm/spCAM5/CAPE/reversible_CAPE_h2_1997_buoyancy+0K_and_rho_from_Tv.nc'
filename4 = '/home/rtm/spCAM5/GCM_OMEGA/GCM_OMEGA_TWP_start_from_top.nc'
filename3 = '/home/rtm/spCAM5/PCP/convective_rain_in_updrafts_and_downdrafts_TWP_Q_GT1E-4_h2_DEF3_1997.nc'

Id1= NCDF_OPEN(filename1)
NCDF_VARGET, Id1,  'CAPE_capefld',      CAPE_capefld
NCDF_VARGET, Id1,  'CAPE_dcapelse',     CAPE_dcapelse

Id3 = NCDF_OPEN(filename3)
NCDF_VARGET, Id3, 'CONV_RAIN',  CONV_PCP1

Id4 = NCDF_OPEN(filename4)
NCDF_VARGET, Id4,   'GCM_OMEGA',   OMEGA

filename5 = '/home/rtm/spCAM5/RELHUM_h2.1997.nc'
Id5 = NCDF_OPEN(filename5)
NCDF_VARGET, Id5,   'RH',   RH


dCAPElse1 = CAPE_DCAPELSE[*,*,*,0]*6.		
CAPE      = CAPE_capefld[*,*,*,0]
CIN	  = CAPE_capefld[*,*,*,2]

;plev values came from the original spCAM nc files, the levels are reversed
plev = [992., 976., 957., 936., 912., 887., 859., 820., 763., 691., 609., 524., 445., 379., 322., 273., 232., 197., 168., 142., 121., 103., 87., 72., 54., 38., 24., 14., 7., 3.]

MC = OMEGA*3600.

;================================================================
;Define Event, and collect data for all events

CAPE_EVENT = FLTARR(25000,count1/85)
CAPE_EVENT[*,*] = 1.E30
CIN_EVENT = FLTARR(25000,count1/85)
CIN_EVENT[*,*] = 1.E30

dCAPEdt_EVENT = FLTARR(25000,count1/85)
dCAPEdt_EVENT[*,*] = 1.E30
dCAPElse_EVENT = FLTARR(25000,count1/85)
dCAPElse_EVENT[*,*] = 1.E30

CONV_PCP_EVENT = FLTARR(25000,count1/85)
CONV_PCP_EVENT[*,*] = 1.E30

RH_EVENT = FLTARR(25000,nlev,count1/85)
RH_EVENT[*,*,*] = 1.E30

MC_EVENT = FLTARR(25000,nlev,count1/85)
MC_EVENT[*,*,*] = 1.E30

COUNT_EVENT = 0L

FOR LON1=1,nlon-2 DO BEGIN 
  FOR LAT1=1,nlat-2 DO BEGIN
    FOR TIME1=18, count1-144 DO BEGIN
;=====================
	MAX_RAIN = 0.
	LSE_max_RAIN = 0. ; dPCAPElse at time when rain is maximum
	LSE_init_RAIN= 0. ; dPCAPElse at time of initiation
;=====================
; IF - detects when rainfall starts and ends, and what is the maximum rainfall

; DEF 1
;       IF (CONV_PCP1[LON1, LAT1, TIME1+1] GT 0.1) AND $
;       (array_tot_pcp[LON1, LAT1, TIME1] LT 0.1) AND (array_tot_pcp[LON1, LAT1, TIME1-1] LT 0.1) AND (array_tot_pcp[LON1, LAT1, TIME1-2] LT 0.1) THEN BEGIN; AND  $
          ; (CONV_PCP1[LON1, LAT1, TIME1-5] LT 0.1) AND (CONV_PCP1[LON1, LAT1, TIME1-6] LT 0.1) AND (CONV_PCP1[LON1, LAT1, TIME1-7] LT 0.1) AND $
;       (CONV_PCP1[LON1, LAT1, TIME1-8] LT 0.1) AND (CONV_PCP1[LON1, LAT1, TIME1-9] LT 0.1) AND (CONV_PCP1[LON1, LAT1, TIME1-10] LT 0.1) AND  $
;       (CONV_PCP1[LON1, LAT1, TIME1-11] LT 0.1) AND (CONV_PCP1[LON1, LAT1, TIME1-12] LT 0.1) AND (CONV_PCP1[LON1, LAT1, TIME1-13] LT 0.1) AND  $
;       (CONV_PCP1[LON1, LAT1, TIME1-14] LT 0.1) AND (CONV_PCP1[LON1, LAT1, TIME1-15] LT 0.1) AND (CONV_PCP1[LON1, LAT1, TIME1-16] LT 0.1) AND  $
;       (CONV_PCP1[LON1, LAT1, TIME1-17] LT 0.1) THEN BEGIN  

; DEF 2
        IF (CONV_PCP1[LON1, LAT1, TIME1+1] GT 0.1) AND $
        (CONV_PCP1[LON1, LAT1, TIME1] LT 0.1) AND (CONV_PCP1[LON1, LAT1, TIME1-1] LT 0.1) AND (CONV_PCP1[LON1, LAT1, TIME1-2] LT 0.1) AND $ ; THEN BEGIN
        (CONV_PCP1[LON1, LAT1, TIME1-3] LT 0.1) AND (CONV_PCP1[LON1, LAT1, TIME1-4] LT 0.1) AND (CONV_PCP1[LON1, LAT1, TIME1-5] LT 0.1) AND $ ;THEN BEGIN
        (CONV_PCP1[LON1, LAT1, TIME1-6] LT 0.1) AND (CONV_PCP1[LON1, LAT1, TIME1-7] LT 0.1) AND (CONV_PCP1[LON1, LAT1, TIME1-8] LT 0.1) AND $
        (CONV_PCP1[LON1, LAT1, TIME1-9] LT 0.1) AND (CONV_PCP1[LON1, LAT1, TIME1-10] LT 0.1) AND (CONV_PCP1[LON1, LAT1, TIME1-11] LT 0.1) AND $
        (CONV_PCP1[LON1, LAT1, TIME1-12] LT 0.1) AND (CONV_PCP1[LON1, LAT1, TIME1-13] LT 0.1) AND (CONV_PCP1[LON1, LAT1, TIME1-14] LT 0.1) AND $
        (CONV_PCP1[LON1, LAT1, TIME1-15] LT 0.1) AND (CONV_PCP1[LON1, LAT1, TIME1-16] LT 0.1) AND (CONV_PCP1[LON1, LAT1, TIME1-17] LT 0.1) THEN BEGIN

;TM++ NOV 01, 2015
;CAPE and dCAPEs precede rainfall by 1 time step
;dPCAPE at 12:00 means forcing between 12:00-12:10. 
;rainfall at 12:10 means rainfall between 12:00-12:10
;thus compare forcing at 12:00 to rainfall at 12:10

;TM--              ;LSE_init_RAIN = dCAPElse1[LON1, LAT1, TIME1+1]
;TM++
		LSE_init_RAIN = dCAPElse1[LON1, LAT1, TIME1]	
                FOR TIME2=TIME1+1,TIME1+72 DO BEGIN
                        IF (CONV_PCP1[LON1, LAT1, TIME2] GT MAX_RAIN) THEN BEGIN
                                MAX_RAIN = CONV_PCP1[LON1, LAT1, TIME2]
;TM--                                LSE_max_RAIN = dCAPElse1[LON1, LAT1, TIME2]
;TM++
				LSE_max_RAIN = dCAPElse1[LON1, LAT1, TIME2-1]
                        ENDIF
;DEF 1
;                       IF (CONV_PCP1[LON1, LAT1, TIME2] LT 0.1) AND (CONV_PCP1[LON1, LAT1, TIME2+1] LT 0.1) THEN GOTO, OUT

;DEF 2 (for total rain change from EQ 0 , to , LT 0.1
                        IF (CONV_PCP1[LON1, LAT1, TIME2] LT 0.1) THEN GOTO, OUT
                ENDFOR
        ENDIF
        OUT:
;=====================
; IF - defines start and end times of an event that satisfies rainfall and dPCAPElse thresholds
; the thresholds are defined using 3-months (May,June,July)
;=====================
;STRONG FORCING STRONG RAINFALL
;lowered the rain threshold from 2.20 to 2.00 and LSE from 150 (75%) to 50 (25%) as Guang suggested (Sept 24); in order to increase the number of events

;DEF 1
;       IF (MAX_RAIN GT 2.0) AND (LSE_max_RAIN GT 50.0) AND (LSE_init_RAIN GT 50.0) AND (LSE_max_RAIN NE 0.0) AND (LSE_init_RAIN NE 0.0) THEN BEGIN

;DEF 2
        IF (MAX_RAIN GT 1.0) THEN BEGIN
;=====================
;TM++ NOV 01, 2015
;CAPE and dCAPEs precede rainfall by 1 time step
;dPCAPE at 12:00 means forcing between 12:00-12:10. 
;rainfall at 12:10 means rainfall between 12:00-12:10
;thus compare forcing at 12:00 to rainfall at 12:10
		;EVENT_LENGHT[COUNT_EVENT] = (TIME2-TIME1)/6. ; convert to hours
		START = TIME1
		FINISH= TIME2
		FOR TIME3=START-18,FINISH DO BEGIN ; from -3 hours to the end of an event time
;TM--			CAPE_EVENT      [COUNT_EVENT, TIME3-START+18] = CAPE[LON1, LAT1, TIME3]	
;TM++
			CAPE_EVENT      [COUNT_EVENT, TIME3-START+18] = CAPE[LON1, LAT1, TIME3]
			CIN_EVENT       [COUNT_EVENT, TIME3-START+18] = CIN[LON1, LAT1, TIME3]
			;IF (CAPE[LON1, LAT1, TIME3] GT 0.) AND (CAPE[LON1, LAT1, TIME3-1] GT 0.) THEN BEGIN	
			;	dCAPEdt_EVENT   [COUNT_EVENT, TIME3-START+18] = (CAPE[LON1, LAT1, TIME3] - CAPE[LON1, LAT1, TIME3-1])*6.
			;ENDIF
;TM--	
;	                IF (dCAPElse1[LON1, LAT1, TIME3] NE 0.) AND (dCAPElsp[LON1, LAT1, TIME3] NE 0.) AND (dCAPEcue[LON1, LAT1, TIME3] NE 0.) AND (dCAPEcup[LON1, LAT1, TIME3] NE 0.) THEN BEGIN    
;                               dCAPEdt_EVENT[COUNT_EVENT, TIME3-START+18] = (dCAPElse1[LON1, LAT1, TIME3] + dCAPElsp[LON1, LAT1, TIME3] + dCAPEcue[LON1, LAT1, TIME3] + dCAPEcup[LON1, LAT1, TIME3])

;TM++ change TIME3 with TIME3-1, do not for the rainfall. Double checked on Sept 11 2017.
			dCAPElse_EVENT  [COUNT_EVENT, TIME3-START+18] = dCAPElse1[LON1, LAT1, TIME3-1]
                        CONV_PCP_EVENT  [COUNT_EVENT, TIME3-START+18] = CONV_PCP1[LON1, LAT1, TIME3]
                        MC_EVENT        [COUNT_EVENT, *, TIME3-START+18] = MC[LON1, LAT1, *,TIME3] ; MC and LTS are at same time step as PCP		        
                        RH_EVENT        [COUNT_EVENT, *, TIME3-START+18] = RH  [LON1, LAT1, *, TIME3]
		ENDFOR

		MAX_RAIN = 0.	
		COUNT_EVENT = COUNT_EVENT + 1	
	ENDIF
;=====================
    ENDFOR
  ENDFOR
ENDFOR

;================================================================
;this part is for -3hours to 1-time step prior to initiation time

MEAN_CAPE0 = FLTARR(18)
MEAN_CAPE0[*] = 'NaN'
MEAN_CIN0 = FLTARR(18)
MEAN_CIN0[*] = 'NaN'
MEAN_dCAPEdt0 = FLTARR(18)
MEAN_dCAPEdt0[*] = 'NaN'
MEAN_RH0 = FLTARR(nlev,18)
MEAN_RH0[*,*] = 'NaN'
MEAN_MC0 = FLTARR(nlev,18)
MEAN_MC0[*,*] = 'NaN'
MEAN_LSE0 = FLTARR(18)
MEAN_LSE0[*] = 'NaN'
MEAN_PCP0 = FLTARR(18)
MEAN_PCP0[*] = 'NaN'

AV_N = 20 ; min number of events
FOR TIME4 = 0,17 DO BEGIN

        INDEX1000 = WHERE((CAPE_EVENT[*,TIME4] NE 0.) AND (CAPE_EVENT[*,TIME4] LT 10000.) AND (dCAPElse_EVENT[*,TIME4] GT -50000.) AND (dCAPElse_EVENT[*,TIME4] LT 50000.) AND (dCAPElse_EVENT[*,TIME4] NE 0.), icnt1000)

        IF (icnt1000 GT AV_N) THEN BEGIN
                MEAN_CAPE0 [TIME4]  = MEAN(CAPE_EVENT[INDEX1000,TIME4])
		MEAN_CIN0  [TIME4]  = MEAN(CIN_EVENT[INDEX1000,TIME4])
                MEAN_LSE0 [TIME4]  = MEAN(dCAPElse_EVENT[INDEX1000,TIME4])
                MEAN_PCP0 [TIME4]  = MEAN(CONV_PCP_EVENT[INDEX1000,TIME4])
	        FOR NL0=0,nlev-1 DO BEGIN
                  MEAN_MC0 [NL0,TIME4]  = MEAN(MC_EVENT[INDEX1000,NL0,TIME4])
		  MEAN_RH0 [NL0,TIME4]  = MEAN(RH_EVENT[INDEX1000,NL0,TIME4])
		ENDFOR
        ENDIF
ENDFOR
;=============================================================================

;new part Jun 21, 2017

nmbr = 20 ; number of elements in X grid, to which rain event times are interpolated to. used 50 and 100, no difference. However, use 20 (since most events last 3-hrs or 18-time step)
CONV_PCP_EVENT_new = FLTARR(3000, nmbr+1)
CONV_PCP_EVENT_new[*,*] = 1.E30

CAPE_EVENT_new = FLTARR(3000, nmbr+1)
CAPE_EVENT_new[*,*] = 1.E30

CIN_EVENT_new = FLTARR(3000, nmbr+1)
CIN_EVENT_new[*,*] = 1.E30

dCAPElse_EVENT_new = FLTARR(3000, nmbr+1)
dCAPElse_EVENT_new[*,*] = 1.E30

RH_EVENT_new = FLTARR(3000, 30, nmbr+1)
RH_EVENT_new[*,*,*] = 1.E30

MC_EVENT_new = FLTARR(3000, 30, nmbr+1)
MC_EVENT_new[*,*,*] = 'NAN' 

FOR newtime1=0,COUNT_EVENT-1 DO BEGIN
 end_event=0
 FOR newtime2=18,142 DO BEGIN
   IF (CONV_PCP_EVENT[newtime1, newtime2+1] GT 1.E5) AND (end_event EQ 0) THEN BEGIN
     end_event=newtime2
     differ = INDGEN(newtime2 - 18 + 1)
     newtime= INDGEN(nmbr+1)*((n_elements(differ)-1)/FLOAT(nmbr))
     ;use INTERPOL
     CONV_PCP_EVENT_new     [newtime1, *] = INTERPOL(CONV_PCP_EVENT    [newtime1, 18:newtime2], differ, newtime)
     CAPE_EVENT_new         [newtime1, *] = INTERPOL(CAPE_EVENT        [newtime1, 18:newtime2], differ, newtime)
     CIN_EVENT_new          [newtime1, *] = INTERPOL(CIN_EVENT         [newtime1, 18:newtime2], differ, newtime)
     dCAPElse_EVENT_new     [newtime1, *] = INTERPOL(dCAPElse_EVENT    [newtime1, 18:newtime2], differ, newtime)
     FOR LV0=0,nlev-1 DO BEGIN
       MC_EVENT_new          [newtime1, LV0, *] = INTERPOL(MC_EVENT        [newtime1, LV0, 18:newtime2], differ, newtime)
       RH_EVENT_new          [newtime1, LV0, *] = INTERPOL(RH_EVENT        [newtime1, LV0, 18:newtime2], differ, newtime)
     ENDFOR
   ENDIF
 ENDFOR
ENDFOR

;================================================================
;for all events, compute mean

MEAN_CAPE1 = FLTARR(nmbr+1)
MEAN_CAPE1[*] = 'NaN'
MEAN_CIN1 = FLTARR(nmbr+1)
MEAN_CIN1[*] = 'NaN'
MEAN_dCAPEdt1 = FLTARR(nmbr+1)
MEAN_dCAPEdt1[*] = 'NaN'
MEAN_LSE1 = FLTARR(nmbr+1)
MEAN_LSE1[*] = 'NaN'
MEAN_PCP1 = FLTARR(nmbr+1)
MEAN_PCP1[*] = 'NaN'

MEAN_RH1 = FLTARR(nlev, nmbr+1)
MEAN_RH1[*,*] = 1.E30
MEAN_MC1 = FLTARR(nlev, nmbr+1)
MEAN_MC1[*,*] = 'NaN'

NUMBER_EVENTS = FLTARR(nmbr+1)

FOR TIME4 = 0,nmbr DO BEGIN

	INDEX1000 = WHERE((CAPE_EVENT_NEW[*,TIME4] NE 0.) AND (CAPE_EVENT_NEW[*,TIME4] LT 10000.) AND (dCAPElse_EVENT_NEW[*,TIME4] GT -50000.) AND (dCAPElse_EVENT_NEW[*,TIME4] LT 50000.) AND (dCAPElse_EVENT_NEW[*,TIME4] NE 0.), icnt1000) 

	IF (icnt1000 GT AV_N) THEN BEGIN
		MEAN_CAPE1 [TIME4]     = MEAN(CAPE_EVENT_NEW[INDEX1000,TIME4])
                MEAN_CIN1 [TIME4]      = MEAN(CIN_EVENT_NEW[INDEX1000,TIME4])
                MEAN_LSE1 [TIME4]      = MEAN(dCAPElse_EVENT_NEW[INDEX1000,TIME4])
                MEAN_PCP1 [TIME4]      = MEAN(CONV_PCP_EVENT_NEW[INDEX1000,TIME4])
                FOR LEV0=0,nlev-1 DO BEGIN
                  MEAN_MC1 [LEV0, TIME4] = MEAN(MC_EVENT_NEW [INDEX1000, LEV0, TIME4])
		  MEAN_RH1 [LEV0, TIME4] = MEAN(RH_EVENT_NEW [INDEX1000, LEV0, TIME4])
                ENDFOR
	ENDIF
	NUMBER_EVENTS [TIME4] = N_ELEMENTS(INDEX1000)

ENDFOR

;===============================================================
;connect two arrays, array1 prior to t=0, and array 2 after t=0

MEAN_CAPE    = FLTARR(18+nmbr+1)
MEAN_CIN     = FLTARR(18+nmbr+1)
MEAN_LSE     = FLTARR(18+nmbr+1)
MEAN_PCP     = FLTARR(18+nmbr+1)
MEAN_MC      = FLTARR(nlev, 18+nmbr+1)
MEAN_RH       = FLTARR(nlev,18+nmbr+1)
MEAN_mRH	= FLTARR(18+nmbr+1)

FOR ntim1=0,17 DO BEGIN
  MEAN_CAPE    [ntim1] = MEAN_CAPE0    [ntim1]
  MEAN_CIN     [ntim1] = MEAN_CIN0     [ntim1]
  MEAN_LSE     [ntim1] = MEAN_LSE0     [ntim1]
  MEAN_PCP     [ntim1] = MEAN_PCP0     [ntim1]
  MEAN_MC      [*,ntim1] = MEAN_MC0    [*,ntim1]
  MEAN_RH       [*,ntim1] = MEAN_RH0   [*,ntim1]
ENDFOR

FOR ntim2=0,nmbr DO BEGIN
  MEAN_CAPE    [ntim2+18] = MEAN_CAPE1    [ntim2]
  MEAN_CIN     [ntim2+18] = MEAN_CIN1     [ntim2]
  MEAN_LSE     [ntim2+18] = MEAN_LSE1     [ntim2]
  MEAN_PCP     [ntim2+18] = MEAN_PCP1     [ntim2]
  MEAN_MC      [*,ntim2+18] = MEAN_MC1    [*,ntim2]
  MEAN_RH      [*,ntim2+18] = MEAN_RH1    [*,ntim2]
ENDFOR

PCP_from_LSE = FLTARR(18+nmbr+1)
FOR ntim2=0,nmbr DO BEGIN
  PCP_from_LSE[ntim2+18] = (MEAN_LSE1[ntim2] - (MEAN_LSE1[0] + (((MEAN_LSE1[nmbr] - MEAN_LSE1[0])/(nmbr+1.))*ntim2)))/32.8 ; 32.8 for CAPE and 17.6 for WCAPE
ENDFOR

FOR nntim=0,n_elements(MEAN_PCP)-1 DO BEGIN
  MEAN_mRH [nntim] = MEAN(MEAN_RH [27:27, nntim])
ENDFOR
;================================================================
;Plot means

D=2.2

rpsopen, 'all_forcing_initiate_TWP_spCAM1997_time_rCAPE_1mm_rain_1E-4_dCAPEdt.eps', /encap, xs=24, ys=20, /inches, /color

Device, Decomposed=0
LoadCT, 33, NColors=(10), Bottom=3
!P.MULTI = [0, 1, 1, 0, 0]

a1=0.15 ;   x1 left column
a2=0.85 ;   x2 left column

b1=0.15 ;   row 1
b2=0.85 ;   

;e1=0.35
;e2=0.45

;f1=0.20
;f2=0.30

;g1=0.05
;g2=0.15

X = INDGEN(18+nmbr+1)-18.
;Y=REVERSE(MEAN_SIGMA)
X2 = FLTARR(5)
Y2 = [-3000.,-1000.,0.,100.,3000.]

;================================================================
;Convective Rainfall
;Figure (a)
;===============
plot,X, MEAN_PCP, XTHICK=8.0, YTHICK=8.0, position=[a1, b1, a2, b2], $
    XRANGE = [-10,nmbr], YRANGE = [0,1.5], XTICKINTERVAL=5, XMINOR=5, YTICKINTERVAL=0.5, YMINOR=1, $
    TITLE='',$
    xtitle='Time',$
    ytitle='Conv PCP (mm h!E-1!N)', $
    charsize=12.5/D,charthick=20.5/D,thick=8.

OPLOT, X2,Y2

;AXIS, YAXIS=1, YRANGE = [0.,150.],YTICKINTERVAL=50, YMINOR=1, COLOR=11, YSTYLE = 1, $
;   YTITLE = 'dCAPE!DLSFT!N (Jkg!E-1!Nh!E-1!N)', charsize=12.5/D,charthick=20.5/D
;oplot, X, (MEAN_LSE)/(100.), linestyle=0, color=11, thick=8

dCAPE = FLTARR(39)
for l0 = 1,37 DO BEGIN
  dCAPE[l0] = (MEAN_LSE[l0+1]-MEAN_LSE[l0-1])/2.
endfor

;smoothing, too noisy prior to time=0
for l0 = 1,37 DO BEGIN
  dCAPE[l0] = (dCAPE[l0]+dCAPE[l0-1]+dCAPE[l0+1])/3.
endfor

AXIS, YAXIS=1, YRANGE = [-5.,10.],YTICKINTERVAL=5, YMINOR=1, COLOR=11, YSTYLE = 1, $
   YTITLE = 'dCAPE!DLSFT!N/dt ', charsize=12.5/D,charthick=20.5/D
oplot, X, (dCAPE)/(10.)+0.5, linestyle=0, color=11, thick=8
;================================================================
rpsclose

END
