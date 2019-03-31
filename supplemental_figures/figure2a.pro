;2D contour plot, dPCAPElse vs omega, and contours show conv PCP
;====================================================================================
PRO tendency, filename1, filename2, ctl_tnt 
COMPILE_OPT IDL2
;====================================================================================

nlon0=9
nlat0=6

filename1a = '/home/rtm/spCAM5/CAPE/reversible_CAPE_h2_1997_buoyancy+0K_and_rho_from_Tv.nc'
filename1b = '/home/rtm/spCAM5/GCM_OMEGA/GCM_OMEGA_TWP_start_from_top.nc'
filename1c = '/home/rtm/spCAM5/PCP/convective_rain_in_updrafts_and_downdrafts_TWP_Q_GT1E-4_h2_DEF3_1997.nc'

Id1a  = NCDF_OPEN(filename1a)
Id1b  = NCDF_OPEN(filename1b)
Id1c  = NCDF_OPEN(filename1c)

NCDF_VARGET, Id1a,   'CAPE_dcapelse',    PCAPElse
NCDF_VARGET, Id1a,   'CAPE_capefld',     PCAPE0
NCDF_VARGET, Id1b,   'GCM_OMEGA',    	 whf
NCDF_VARGET, Id1c,   'TOTAL_RAIN',        pchf

pchf0     = FLTARR(7*4*12240)
PCAPElse0 = FLTARR(7*4*12240)
PCAPE     = FLTARR(7*4*12240)
CIN       = FLTARR(7*4*12240)
OMEGA0    = FLTARR(7*4*12240)    ; levels start from surface -> top, but MSE levels start from top -> surface

counter0=FLTARR(1)

FOR l0=1,nlon0-2 DO BEGIN
 FOR l1=1,nlat0-2 DO BEGIN
  FOR l2=1,12239 DO BEGIN
   pchf0 [counter0] = pchf[l0,l1,l2]
   PCAPElse0[counter0] = PCAPElse[l0,l1,l2-1,0]*6.
   PCAPE    [counter0] = PCAPE0  [l0,l1,l2,1] + PCAPE0[l0,l1,l2,2]
   CIN      [counter0] = PCAPE0  [l0,l1,l2,2]
   OMEGA0   [counter0] = whf[l0,l1,1,l2]*3600.
   counter0=counter0+1
  ENDFOR
 ENDFOR
ENDFOR

FOR ll0=0,n_elements(PCAPE)-1 DO BEGIN
  IF (PCAPElse0[ll0] EQ 0.) THEN BEGIN
    PCAPElse0[ll0] = 1.E30
  ENDIF
  IF (PCAPE[ll0] EQ 0.) THEN BEGIN
    PCAPE[ll0] = 1.E30
  ENDIF
ENDFOR
;===================================================================================
;Define min and max on X, Y axis

VAR1 = PCAPElse0 
VAR2 = OMEGA0 
INDEX1 = WHERE((VAR1[*] GT -10000.) AND (VAR1[*] LT 10000.) AND (VAR1[*] NE 0.))
INDEX2 = WHERE((VAR2[*] GT -10000.) AND (VAR2[*] LT 10000.) AND (VAR2[*] NE 0.))

X_MAX =  500.
X_MIN = -100.
Y_MAX =  250.
Y_MIN = -250.

;Define number of bins on X and Y axis
N_BINS = 40.

;Array that will be contoured 
RAIN = FLTARR(N_BINS, N_BINS)
RAIN[*,*] = 0.

STEP_X=0.
STEP_X=(X_MAX-X_MIN)/(N_BINS)
X1 = FLTARR(N_BINS)

STEP_Y=0.
STEP_Y=(Y_MAX-Y_MIN)/(N_BINS)
Y1 = FLTARR(N_BINS)

MIN_RAIN=0.1
FOR XX=0,N_BINS-1 DO BEGIN
	FOR YY=0,N_BINS-1 DO BEGIN
		INDEX = WHERE((VAR2[*] GT (Y_MIN+(STEP_Y*YY))) AND (VAR2[*] LE (Y_MIN+(STEP_Y*(YY+1)))) AND (VAR1[*] GT (X_MIN+(STEP_X*XX))) AND (VAR1[*] LE (X_MIN+(STEP_X*(XX+1)))) AND (pchf0[*] GT MIN_RAIN))		
                IF (N_ELEMENTS(INDEX) GT 20.) THEN BEGIN
                        RAIN [YY,XX] = MEAN(pchf0[INDEX])
                ENDIF		
		;RAIN [YY,XX] = MEAN(array_pcp[WHERE((array_cape[*] GT (Y_MIN+(STEP_Y*YY))) AND (array_cape[*] LE (Y_MIN+(STEP_Y*(YY+1)))) AND (array_dcapelsp[*] GT (X_MIN+(STEP_X*XX))) AND (array_dcapelsp[*] LE (X_MIN+(STEP_X*(XX+1)))) AND (array_dcapelse[*] GT LSE_MIN) AND (array_dcapelse[*] LT LSE_MAX) AND (array_cape[*] GT 0.) AND (array_pcp[*] GT MIN_RAIN) AND (array_pcp[*] LT MAX_RAIN) AND (array_dcapelse[*] NE 0.))])
	ENDFOR
print, XX/N_BINS*100., '%'
ENDFOR

FOR XX=0,N_BINS-1 DO BEGIN
	 X1[XX]= ((X_MIN+(STEP_X*XX)) + (X_MIN+(STEP_X*(XX+1))))/2. ; mean interval value
ENDFOR

FOR YY=0,N_BINS-1 DO BEGIN
         Y1[YY]= ((Y_MIN+(STEP_Y*YY)) + (Y_MIN+(STEP_Y*(YY+1))))/2. ; mean interval value
ENDFOR

;===================================================================================
rpsopen, '2D_pdf_reversible_dCAPElse_vs_OMEGA_spCAM5_h2_pcp_e-4_small_binsize_tp', /encap, xs=14, ys=10, /inches, /color

Device, Decomposed=0
!P.MULTI = [0, 1, 1, 0, 0]

D=2.8

a1=0.18 ;   x1 left column
a2=0.78 ;   x2 left column

b1=0.15 ;   row 1
b2=0.90 ;   

;c1=0.65 ;   row 2
;c2=0.75

;d1=0.50 ;   row 3
;d2=0.60

;e1=0.35 ;   row 4
;e2=0.45

;f1=0.20 ;   row 5
;f2=0.30

;g1=0.05
;g2=0.15

;================================================================

levels = 10.  ; Number of contours
LoadCT, 33, NColors=levels, Bottom=3

max_value = 3.
min_value = 0.
step = (max_value - min_value) / levels
userLevels = IndGen(levels) * step + min_value
contour,(transpose(RAIN[*,*])), X1,Y1, /fill, XRANGE = [X_MIN, X_MAX], YRANGE = [Y_MIN, Y_MAX], C_Colors=(Indgen(levels)+3), Background=255, $
LEVELS=userLevels, $
   XSTYLE = 1, YSTYLE = 1, XTHICK=12.5/D, YTHICK=12.5/D, position=[a1, b1, a2, b2], Color=2, $
   CHARSIZE=8.5/D,CHARTHICK=15.5/D,$
   TITLE = 'spCAM5', $
   XTITLE = 'dCAPE!DLSE!N (Jkg!E-1!Nh!E-1!N)',$
;   YTITLE = 'dWCAPE!DLSE!N (Jm!E-3!Nh!E-1!N)'
;   YTITLE = 'CAPE (Jkg!E-1!N)'
   YTITLE = '!4x!X (Pah!E-1!N)'
Contour, (transpose(RAIN[*,*])), x1, y1, /Overplot, Levels=userLevels ;, C_THICK=12.5/D, C_CHARSIZE=10.5/D, C_CHARTHICK=12.5/D, /Follow, color=2

ax1=0.95
ax2=0.98
ay1=0.20
ay2=0.90

  ;COLORBAR,  NColors=levels, REVERSE=REVERSE(Indgen(levels)+3), Bottom=3, Divisions=10, VERTICAL=[0.94, 0.75, 0.98, 0.97], color=2, YTITLE = 'Mc [kg/m2/h]', $
  COLORBAR,  NColors=levels, C_Colors=(Indgen(levels)+3), Bottom=3, Divisions=10, YTITLE = 'Total Precipitation (mmh!E-1!N)', $
  RANGE=[min_value,max_value],  CHARSIZE=8.5/D,CHARTHICK=15.5/D, VERTICAL=[0,0,0,0], $
  TICKNAMES=[''], ax1, ay1, ax2, ay2; LEAVE TICKNAME empty so the ticks will be from min to max values

;===================================================================================

rpsclose

END

