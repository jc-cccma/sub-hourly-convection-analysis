PRO buoyan_dilute_rcape, q, t, p, z, zf, capeten
;========================================================
;subroutine buoyan_dilute(q,T,p,z, pf,cape, msg)
;----------------------------------------------------------------------- 
;  
; Purpose: 
; Calculates CAPE the lifting condensation level and the convective top
; where buoyancy is first -ve.
; 
; Method: Calculates the parcel temperature based on a simple constant
; entraining plume model. CAPE is integrated from buoyancy.
; 09/09/04 - Simplest approach using an assumed entrainment rate for 
;            testing (dmpdp). 
; 08/04/05 - Swap to convert dmpdz to dmpdp  
;
; SCAM Logical Switches - DILUTE:RBN - Now Disabled 
; ---------------------
; switch(1) = .T. - Uses the dilute parcel calculation to obtain tendencies.
; switch(2) = .T. - Includes entropy/q changes due to condensate loss and freezing.
; switch(3) = .T. - Adds the PBL Tpert for the parcel temperature at all levels.
; 
; References:
; Raymond and Blythe (1992) JAS 
; 
; Author:
; Richard Neale - September 2004
; 
;-----------------------------------------------------------------------
nlev=n_elements(p)

tpert=0.         ;perturbation temperature by pbl processes

; output arguments
;=================
tp=FLTARR(nlev)       ;parcel temperature
qstp=FLTARR(nlev)     ;saturation mixing ratio of parcel
tl=0.            ;parcel temperature at lcl
lcl=0.        
lel=0.       
lon=0.          ;level of onset of deep convection
mx=0.           ;level of max moist static energy

;--------------------------Local Variables------------------------------

tv1=FLTARR(nlev)       
tpv=FLTARR(nlev)     
buoy=FLTARR(nlev)

a1=0.
a2=0.
estp=0.
pl=0.
plexp=0.
hmax=0.
hmn=0.
y=0.
hmn_g=FLTARR(nlev)
dAtp=0.          ;CAPE change due to Tp change
dAtv=0.	         ;CAPE change due to Tv change
dAlel=0.         ;CAPE change due to L.E.L. change

;logical plge600(ntime)
idoc=0.
knt=0.
lelten=INTARR(5)

e=0.
rhd=0.

;-----------------------------------------------------------------------
;

rd = 287.04
grav = 9.80616
cp = 1004.64

a = 21.656
b = 5418.
c1 = 6.112
c2 = 17.67
c3 = 243.5

rgrav = 1./grav
eps1 = 0.622
qmin = 1.E-20
tfreez = 273.16
rl = 2.5104E6

;=================================
FOR n = 0,4 DO BEGIN
        	lelten[n] = nlev-1	;lelten[n]=30
ENDFOR

	lon = nlev			;level of onset of deep convection
     	knt = 0
      	idoc = 1
      	lel= nlev
      	mx = lon
      	hmax = 0.

	FOR j = 0,nlev-1 DO BEGIN
        	tp[j]=t[j]
         	qstp[j]=q[j]	;Parcel water vapour (sat value above lcl).
         	tv1[j]=t[j]*(1.+1.608*q[j])/(1.+q[j])
         	tpv[j]=tv1[j]
         	buoy[j]=0.0
	ENDFOR
;find maximum moist static energy (hmax) and level where hmax (mx)
;=================================================================
FOR k = nlev-1, 0,-1 DO BEGIN
        	hmn = cp*t[k] + grav*z[k] + rl*q[k]
        	hmn_g[k] = hmn
;k=16 is at 600 hPa        	
		IF (k GE 16) AND (k LE (lon-1)) AND (hmn GT hmax) THEN BEGIN	;lon=level of onset of deep convection
            		hmax = hmn
            		mx = k
		ENDIF
ENDFOR
; LCL dilute calculation - initialize to mx(i)
; Determine lcl in parcel_dilute and get pl,tl after parcel_dilute
; Original code actually sets LCL as level above wher condensate forms.
; Therefore in parcel_dilute lcl(i) will be at first level where qsmix < qtmix.

         lcl= mx
         tl = t[mx]
         pl = p[mx]

;mse_max_level = mx
;====================================================
;print,'lcl', lcl
msg = 0 ; number of missing levels at the top . This is because in CAM 2 levels were missing and msg=2
mx1=nlev-1 ; launch level is at the surface

;when instead (pl, tl, lcl) I am using (pl1, tl1, lcl1) in parcel_dilute, i get same tpv values
;pl1=nlev-1
;tl1=nlev-1
;lcl1=nlev-1
	 parcel_dilute_reversible, msg, mx1, p, t, q,tpert, tp, tpv, qstp, pl, tl, lcl
;        call parcel_dilute(lchnk, ncol, msg, mx, p, t, q, tpert, tp, tpv, qstp, pl, tl, lcl)
;print,'-------------------'
; if lcl is above the nominal level of non-divergence (600 mbs),
; no deep convection is permitted (ensuing calculations
; skipped and cape retains initialized value of zero).
;
      	;FOR i = 0,ntime-1 DO BEGIN
	;**CHECK this...comment this  plge600[i] = pl[i] GE 600.
      	;ENDFOR
	;because I wont use previous IF, the next IF (k LT mx[i]) AND (plge600[i]) THEN BEGIN
	;will now be, IF (k LT mx[i]) AND (pl(i) GE 600.) THEN BEGIN
;
; Main buoyancy calculation.
;
      FOR k = nlev-1, 0,-1 DO BEGIN
    ;IF (k LT mx[i]) AND (plge600[i]) THEN BEGIN
           IF (pl GE 600.) THEN BEGIN 
		  ;Define buoy from launch level to cloud top.
            	tv1[k] = t[k]*(1.+1.608*q[k])/ (1.+q[k])
            	buoy[k] = tpv[k] - tv1[k]; + 0.5  ; +0.5K or not?
	   ENDIF ELSE BEGIN
            qstp[k] = q[k]
            tp[k]   = t[k]
            tpv[k]  = tv1[k]
          ENDELSE
      ENDFOR
;-------------------------------------------------------------------------------
;get LNB

      FOR k =0,nlev-1 DO BEGIN
	   IF (k LT lcl) AND (idoc NE 0) THEN BEGIN
            IF (buoy[k+1] GT 0.) AND (buoy[k] LE 0.) THEN BEGIN
               ;knt[i] = min(5,knt[i] + 1)
		test_min=intarr(2)
		test_min[0]=knt+1
		test_min[1]=5
		knt = min(test_min)
               lelten[(knt-1)] = k
            ENDIF
          ENDIF
       ENDFOR

;get first buoyancy level above the surface
buoy1 = nlev-1
find_b = 1
;nlev-15 is at 3000m
      FOR k =nlev-1, nlev-15,-1 DO BEGIN
           IF (find_b NE 0) AND (idoc NE 0) THEN BEGIN
            IF (buoy[k] LE 0.) AND (buoy[k-1] GT 0.) THEN BEGIN
               buoy1 = k-1
	       find_b = 0
            ENDIF
          ENDIF
      ENDFOR


;if (mx LT 31) then begin
;	print, buoy
;stop	

;endif
;=======================================================
; add air density. See Guang's Email on April 13, 2015.

rho = FLTARR(nlev)
Rd =  287.04                    ; gas constant (J/kg/K)

FOR LEVEL1=0,nlev-1 DO BEGIN
        rho[LEVEL1] = (p[LEVEL1]*100.)/(Rd * tv1[LEVEL1])     ; air density, kg/m3 
ENDFOR
;print, 'rho', rho
;print, 'p', p
;print, 't',t
;=======================================================
;lelten[0] is the level above LNB
      	FOR k = (lelten[0]+1),nlev-1 DO BEGIN
	   IF (idoc NE 0) AND (buoy[k] GT -20.) AND (buoy[k] LT 20.) THEN BEGIN
           	capeten[0] = capeten[0] + grav*(buoy[k]/tv1[k])*(zf[k]-zf[k+1])
	   ENDIF
      	ENDFOR

        FOR k = (lelten[0]+1), buoy1 DO BEGIN
           IF (idoc NE 0) AND (buoy[k] GT 0.) AND (buoy[k] LT 20.) THEN BEGIN
                capeten[1] = capeten[1] + grav*(buoy[k]/tv1[k])*(zf[k]-zf[k+1])
           ENDIF
        ENDFOR

        FOR k = buoy1, nlev-1 DO BEGIN
           IF (idoc NE 0) AND (buoy[k] GT -20.) AND (buoy[k] LT 0.) THEN BEGIN
                capeten[2] = capeten[2] + grav*(buoy[k]/tv1[k])*(zf[k]-zf[k+1])
           ENDIF
        ENDFOR

;if (mx lt 31) then begin
;print, mx
;print, buoy
;print, buoy1
;print, lelten
;print,'=============='
;print, 'capeten'
;print, capeten[2]
;print, '============='
;print, 'rho'
;print, rho
;print, '============'
;print, 'tv1'
;print, tv1
;print, '============'
;print, 'zf'
;print, zf

;stop
;endif


      RETURN ; subroutine also returns values
END
