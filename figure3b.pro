;2D contour plot, dWCAPElse vs omega, and contours show conv PCP
;====================================================================================
PRO tendency, filename1, filename2, ctl_tnt 
COMPILE_OPT IDL2
;====================================================================================
filename1a = '/home/rtm/GCM/CanAM4_hires/ensamble/pchf/TWP_sa_jcl_hf_1997_pchf_gs.nc' ; Convective Precipitation  kg m-2 s-1 
filename1h = '/home/rtm/GCM/CanAM4_hires/ensamble/whf/TWP_sa_jcl_hf_1997_whf_gs.nc'
filename1c = '/home/rtm/GCM/CanAM4_hires/ensamble/rCAPE/TWP_sa_jcl_hf_1997_rCAPE_gs.nc'
filename1d = '/home/rtm/GCM/CanAM4_hires/ensamble/rCAPElse/TWP_sa_jcl_hf_1997_rCAPElse_gs.nc'

Id1a  = NCDF_OPEN(filename1a)
Id1c  = NCDF_OPEN(filename1c)
Id1d  = NCDF_OPEN(filename1d)
Id1h  = NCDF_OPEN(filename1h)

NCDF_VARGET, Id1a,   'PCHF',    pchf
NCDF_VARGET, Id1c,   'CAPE',   PCAPE
NCDF_VARGET, Id1d,   'rCAPElse', PCAPElse
NCDF_VARGET, Id1h,   'WHF',     whf

pchf = pchf*3600. ; convert to kg/m2/h
whf = whf*3600.   ; convert to Pa/h

PCAPElse = (PCAPElse-PCAPE)*4. ; convert to  J/kg/h
;PCAPElsp = (PCAPElsp-PCAPE)*4.

nlon = 9 ;
nlat = 4 ;
ntime = 44160 ;

pchf0     = FLTARR(nlon*nlat*ntime)
PCAPE0    = FLTARR(nlon*nlat*ntime)
PCAPElse0 = FLTARR(nlon*nlat*ntime)
OMEGA0    = FLTARR(nlon*nlat*ntime)

counter0=FLTARR(1)

FOR l0=0,nlon-1 DO BEGIN
 FOR l1=0,nlat-1 DO BEGIN
  FOR l2=1,ntime-2 DO BEGIN
   pchf0 [counter0] = pchf[l0,l1,l2]
   PCAPE0[counter0] = PCAPE[l0,l1,l2]
   PCAPElse0[counter0] = PCAPElse[l0,l1,l2]
   OMEGA0   [counter0] = whf[l0,l1,47,l2]
   counter0=counter0+1
  ENDFOR
 ENDFOR
ENDFOR

INDEX = WHERE((pchf0[*] gt .01) and (pchf[*] lt 100.) and (PCAPElse0[*] gt -5000.) and (PCAPElse0[*] lt 5000.))
print, n_elements(INDEX)
print, correlate(pchf0[INDEX], PCAPElse0[INDEX])
INDEX = WHERE((pchf0[*] gt .01) and (pchf[*] lt 100.) and (PCAPE0[*] gt -5000.) and (PCAPE0[*] lt 5000.))
print, correlate(pchf0[INDEX], PCAPE0[INDEX])

;===================================================================================
;Define min and max on X, Y axis
X_MAX =   500.
X_MIN =  -100.
Y_MAX =  600.
Y_MIN = -600.

VAR_X = PCAPElse0
VAR_Y = OMEGA0 

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
		INDEX = WHERE((VAR_Y[*] GT (Y_MIN+(STEP_Y*YY))) AND (VAR_Y[*] LE (Y_MIN+(STEP_Y*(YY+1)))) AND (VAR_X[*] GT (X_MIN+(STEP_X*XX))) AND (VAR_X[*] LE (X_MIN+(STEP_X*(XX+1)))) AND (pchf0[*] GT MIN_RAIN))		
                IF (N_ELEMENTS(INDEX) GT 50.) THEN BEGIN ; use GT 20 for 1 ensemble and GT 50 for 5 enembles
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
rpsopen, 'figure2b_ensemble.eps', /encap, xs=14, ys=10, /inches, /color

Device, Decomposed=0
!P.MULTI = [0, 1, 1, 0, 0]

D=2.8

a1=0.18 ;   x1 left column
a2=0.78 ;   x2 left column

b1=0.15 ;   row 1
b2=0.90 ;   

;================================================================

levels = 12.  ; Number of contours
LoadCT, 33, NColors=levels, Bottom=3

max_value = 2.4
min_value = 0.
step = (max_value - min_value) / levels
userLevels = IndGen(levels) * step + min_value
contour,(transpose(RAIN[*,*])), X1,Y1, /fill, XRANGE = [X_MIN, X_MAX], YRANGE = [Y_MIN, Y_MAX], C_Colors=(Indgen(levels)+3), Background=255, $
LEVELS=userLevels, $
   XSTYLE = 1, YSTYLE = 1, XTHICK=12.5/D, YTHICK=12.5/D, position=[a1, b1, a2, b2], Color=2, $
   CHARSIZE=8.5/D,CHARTHICK=15.5/D,$
   TITLE = 'CanAM4.3', $
   XTITLE = 'dCAPE!DLSE!N (Jkg!E-1!Nh!E-1!N)',$
   YTITLE = '!4x!X (Pah!E-1!N)'
;   YTITLE = 'CAPE (Jkg!E-1!N)'
Contour, (transpose(RAIN[*,*])), x1, y1, /Overplot, Levels=userLevels ;$, C_THICK=12.5/D, C_CHARSIZE=10.5/D, C_CHARTHICK=12.5/D, /Follow, color=2

ax1=0.95
ax2=0.98
ay1=0.20
ay2=0.90

  ;COLORBAR,  NColors=levels, REVERSE=REVERSE(Indgen(levels)+3), Bottom=3, Divisions=10, VERTICAL=[0.94, 0.75, 0.98, 0.97], color=2, YTITLE = 'Mc [kg/m2/h]', $
  COLORBAR,  NColors=levels, C_Colors=(Indgen(levels)+3), Bottom=3, Divisions=12, YTITLE = 'Conv. precipitation (mmh!E-1!N)', $
  RANGE=[min_value,max_value],  CHARSIZE=8.5/D,CHARTHICK=15.5/D, VERTICAL=[0,0,0,0], $
  TICKNAMES=[''], ax1, ay1, ax2, ay2; LEAVE TICKNAME empty so the ticks will be from min to max values

;===================================================================================

rpsclose

END

