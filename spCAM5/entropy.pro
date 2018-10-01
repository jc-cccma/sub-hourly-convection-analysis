FUNCTION entropy, TK, p, qtot
;========================================================

; TK(K),p(mb),qtot(kg/kg)
; from Raymond and Blyth 1992
;
	qv=0.
	qsat=0.
	e=0.
	esat=0.
	L=0.

        pref = 1000.0           ; mb
        eref = 6.106            ; sat p at tfreez (mb)
        cpres=1004.64

        tfreez=273.16
        cpwv=2.08 * 1000.  ; specific heat capacity of water vapor at temperature of 373.16K  (J/(kg k))
        cpliq=4.18 * 1000.
        rl = 2.5104E6
        c1 = 6.112
        c2 = 17.67
        c3 = 243.5
        eps1 = 0.622
        rgas=287.04   ; gas constent of air (J/(kg k)), same as rd mentioned above
        rwat=461.5046

;------------------

        L = rl - (cpliq - cpwv)*(TK-tfreez)         ; T IN CENTIGRADE
        ; Replace call to satmixutils.
        esat = c1*exp(c2*(TK-tfreez)/(c3+TK-tfreez))       ;esat(T) in mb
        qsat=eps1*esat/(p-esat)            ;Sat. mixing ratio (in kg/kg).
	min_test=fltarr(2)
	min_test[0]=qtot
	min_test[1]=qsat
	qv = min(min_test)                ; Partition qtot into vapor part only.
	
	e = qv*p / (eps1 +qv)
        entropy = (cpres + qtot*cpliq)*(alog(TK/tfreez))-rgas*(alog((p-e)/pref)) + (L*qv/TK) - (qv*rwat*(alog(qv/qsat)))
return, entropy
END
