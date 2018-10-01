PRO ientropy,s,p,qt,T,qsat,Tfg

;========================================================
; p(mb), Tfg/T(K), qt/qv(kg/kg), a(J/kg). 
; Inverts entropy, pressure and total water qt 
; for T and saturated vapor mixing ratio
; 
 qv=0.
 Ts=0.
 dTs=0.
 fs1=0.
 fs2=0.
 esat=0.
 L=0.
 e=0.
 i=0
 cpres=0.
;----------------------
        LOOPMAX = 100                   ;* max number of iteration loops 
        ; Values for entropy
        pref = 1000.0          		; mb ref pressure.
        eref = 6.106           		; sat p at tfreez (mb)
      	c1 = 6.112
      	c2 = 17.67
      	c3 = 243.5
        tfreez=273.16
      	cpwv=2.08 * 1000.  		; specific heat capacity of water vapor at temperature of 373.16K  (J/(kg k))
      	cpliq=4.18 * 1000.
      	rl = 2.5104E6
      	eps1 = 0.622
      	rgas=287.04  			; gas constent of air (J/(kg k)), same as rd mentioned above
      	cpres=1004.64
        rwat=461.5046
;----------------------
        ; Invert the entropy equation -- use Newton's method

        Ts = Tfg           ; Better first guess based on Tprofile from conv.
	FOR i=0, LOOPMAX-1 DO BEGIN
       
           L = rl - (cpliq - cpwv)*(Ts-tfreez)
           esat = c1*exp(c2*(Ts-tfreez)/(c3+Ts-tfreez)) ; Bolton (eq. 10)
           qsat = eps1*esat/(p-esat)
	
	   test_min=FLTARR(2)
	   test_min[0]=qt
	   test_min[1]=qsat
           qv = min(test_min)

           e = qv*p / (eps1 +qv)  ; Bolton (eq. 16)
           fs1=(((cpres + qt*cpliq)*alog((Ts)/tfreez)) - (rgas*alog((p-e)/pref)) + (L*qv/(Ts)) - (qv*rwat*alog(qv/qsat)) - (s))
	   L = (rl) - (cpliq - cpwv)*(Ts-1.-tfreez)
	   esat = c1*exp(c2*(Ts-1.-tfreez)/(c3+Ts-1.-tfreez))
           qsat = eps1*esat/(p-esat)

	   test_min1=FLTARR(2)
           test_min1[0]=qt
           test_min1[1]=qsat
           qv = min(test_min1) 
	   e = qv*p /(eps1 +qv)
           fs2 = ((cpres + qt*cpliq)*alog((Ts-1.)/tfreez)) - (rgas*alog((p-e)/pref)) + (L*qv/(Ts-1.)) - (qv*rwat*alog(qv/qsat)) - (s)
           dTs = fs1/(fs2 - fs1)
	   Ts  = Ts+dTs
	   IF (abs(dTs) LT 0.01) THEN GOTO, OUT

;	   IF (abs(dTs) LT 0.01) exit converge
           IF (i EQ (LOOPMAX-1)) THEN BEGIN
                ;print,'**** IENTROPY: Tmix did not converge ****'
	   stop
           ENDIF
 
	ENDFOR
	OUT: 
	;Replace call to satmixutils.
        esat = c1*exp(c2*(Ts-tfreez)/(c3+Ts-tfreez))
        qsat=eps1*esat/(p-esat)

           test_min2=FLTARR(2)
           test_min2[0]=qt
           test_min2[1]=qsat
           qv = min(test_min2)               ;       /* check for saturation */
T = Ts
RETURN

END
