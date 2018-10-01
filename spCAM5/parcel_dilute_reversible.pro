PRO parcel_dilute_reversible, limconv, klaunch, p, t, q, tpert, tp, tpv, qstp, pl, tl, lcl
;========================================================
; Define tpv within this routine.
;================================
;Routine  to determine 
;   1. Tp   - Parcel temperature
;   2. qstp - Saturated mixing ratio at the parcel temperature.
nlev=n_elements(p)
;==================================================
; Have to be careful as s is also dry static energy.
; If we are to retain the fact that CAM loops over grid-points in the internal
; loop then we need to dimension sp,atp,mp,xsh2o with ntime.
;==================================================
      tmix=FLTARR(nlev)        ; Tempertaure of the entraining parcel.
      qtmix=FLTARR(nlev)       ; Total water of the entraining parcel.
      qsmix=FLTARR(nlev)       ; Saturated mixing ratio at the tmix.
      smix=FLTARR(nlev)        ; Entropy of the entraining parcel.
      xsh2o=FLTARR(nlev)       ; Precipitate lost from parcel.
      ds_xsh2o=FLTARR(nlev)    ; Entropy change due to loss of condensate.
      ds_freeze=FLTARR(nlev)   ; Entropy change due to freezing of precip.

      mp=0.   ; Parcel mass flux.
      qtp=0.   ; Parcel total water.
      sp=0.    ; Parcel entropy.

      sp0=0.    ; Parcel launch entropy.
      qtp0=0.   ; Parcel launch total water.
      mp0=0.   ; Parcel launch relative mass flux.

      lwmax=0.      ; Maximum condesate that can be held in cloud before rainout.
      dmpdp=0.      ; Parcel fractional mass entrainment rate (/mb).
      ;dmpdpc=0.    ; In cloud parcel mass entrainment rate (/mb).
      dmpdz=0.      ; Parcel fractional mass entrainment rate (/m)
      dpdz=0.
      dzdp=0.	    ; Hydrstatic relation and inverse of.
      senv=0.       ; Environmental entropy at each grid point.
      qtenv=0.      ; Environmental total water "   "   ".
      penv=0.       ; Environmental total pressure "   "   ".
      tenv=0.       ; Environmental total temperature "   "   ".
      new_s=0.	    ; Hold value for entropy after condensation/freezing adjustments.
      new_q=0.	    ; Hold value for total water after condensation/freezing adjustments.
      dp=0.         ;Layer thickness (center to center)
      tfguess=0.    ;First guess for entropy inversion - crucial for efficiency!
      tscool=0.     ;Super cooled temperature offset (in degC) (eg -35).

      	qxsk=0.
	qxskp1=0.        ;LCL excess water (k, k+1)
        dsdp=0.
	dqtdp=0.
	dqxsdp=0.	; LCL s, qt, p gradients (k, k+1)
        slcl=0.
	qtlcl=0.
	qslcl=0. 	; LCL s, qt, qs values.

        ;Number of iterations for condensation/freezing loop.
	ii=0   		; Loop counters.

;======================================================================
;    SUMMARY
;
;  9/9/04 - Assumes parcel is initiated from level of maxh (klaunch)
;           and entrains at each level with a specified entrainment rate.
;
; 15/9/04 - Calculates lcl(i) based on k where qsmix is first < qtmix.          
;
;======================================================================
;
; Set some values that may be changed frequently.
;   
      	cpliq=4.18 * 1000.
        tfreez=273.16
        nit_lheat = 2 		; iterations for ds,dq changes from condensation freezing.
        dmpdz=0.0;-1.e-3;         ; Entrainment rate. (-ve for /m) , if =0, no entrainment
	;dmpdpc = 3.e-2   	; In cloud entrainment rate (/mb).
        lwmax = 1000. ; 1.e-3    	; Need to put formula in for this. if lxmax=1000 keep condensate.
        tscool = -200.0   		; Temp at which water loading freezes in the cloud. if tscool=-200.0 ignore freezing

        FOR k = nlev-1,limconv,-1 DO BEGIN
        	qtmix[k]=0.
        	smix[k]=0.
        ENDFOR

        qtenv = 0.
        senv = 0.
        tenv = 0.
        penv = 0.

        	qtp0 = 0.
        	sp0 = 0.
        	mp0 = 0.

        	qtp = 0.
        	sp = 0.
        	mp = 0.

        new_q = 0.
        new_s = 0.
        latice=334*1000.	;Latice is latent heat of fusion of ice, (J/kg)
      	rgas=287.04   		;gas constent of air (J/(kg k)), same as rd mentioned above
      	grav = 9.80616
;===========================================
; **** Begin loops ****

         tpert=0.


        FOR k = nlev-1,limconv,-1 DO BEGIN
;===========================================
;Initialize parcel values at launch level.
;===========================================
      			IF (k EQ klaunch) THEN BEGIN
         			qtp0 = q[k]				;Parcel launch total water (assuming subsaturated) OK????.
         			sp0  = entropy(t[k],p[k],qtp0)	;Parcel launch entropy.
         			mp0  = 1.					;Parcel launch relative mass (i.e. 1 parcel stays 1 parcel for dmpdp=0, undilute). 

         			smix[k]  = sp0
         			qtmix[k] = qtp0
         			tfguess = t[k]
				TMIX_IN = tmix[k]
				QSMIX_IN= qsmix[k]
				ientropy, smix[k],p[k],qtmix[k],TMIX_IN,QSMIX_IN,tfguess
				tmix[k]=TMIX_IN
				qsmix[k]=QSMIX_IN

      			ENDIF
;===========================================
;Entraining levels

      			IF (k LT klaunch) THEN BEGIN

;Set environmental values for this level.                 

         			dp = (p[k]-p[k+1])
;In -ve mb as p decreasing with height - difference between center of layers.
         			qtenv = 0.5*(q[k]+q[k+1])         ; Total water of environment.
         			tenv  = 0.5*(t[k]+t[k+1])
         			penv  = 0.5*(p[k]+p[k+1])
         			senv  = entropy(tenv,penv,qtenv)   ;Entropy of environment.   
;Determine fractional entrainment rate /pa given value /m.
         			dpdz = -(penv*grav)/(rgas*tenv) ;in mb/m since  p in mb.
         			dzdp = 1./dpdz                  ; in m/mb

;TONI added mixing from Bechtold 2014
;dmpdz=-1.8e-3*(1.3-RH[k])*((rh2mixr(100.,p[k],t[k])/rh2mixr(100.,p[lcl],t[lcl]))^3.)


         			dmpdp = dmpdz*dzdp              ; /mb Fractional entrainment

; Sum entrainment to current level
; entrains q,s out of intervening dp layers, in which linear variation is assumed
; so really it entrains the mean of the 2 stored values.
         			sp  = sp  - dmpdp*dp*senv
         			qtp = qtp - dmpdp*dp*qtenv
         			mp  = mp  - dmpdp*dp

;Entrain s and qt to next level.

         			smix[k]  = (sp0  +  sp) / (mp0 + mp)
         			qtmix[k] = (qtp0 + qtp) / (mp0 + mp)

;Invert entropy from s and q to determine T and saturation-capped q of mixture.
;t(i,k) used as a first guess so that it converges faster.
;        print *,'2 call:', smix(i,k), sp0(i),sp(i), mp0(i), mp(i)
         			tfguess = tmix[k+1]
				TMIX_IN = tmix[k]
                                QSMIX_IN= qsmix[k]
     				ientropy,smix[k],p[k],qtmix[k],TMIX_IN,QSMIX_IN,tfguess
		                tmix[k]=TMIX_IN
                                qsmix[k]=QSMIX_IN
;STOP HERE	print, tmix_in, tfguess	
; Determine if this is lcl of this column if qsmix <= qtmix.
; FIRST LEVEL where this happens on ascending.

         			IF (qsmix[k] LE qtmix[k]) AND (qsmix[k+1] GT qtmix[k+1]) THEN BEGIN
            				lcl = k
            				qxsk   = qtmix[k] - qsmix[k]
            				qxskp1 = qtmix[k+1] - qsmix[k+1]
            				dqxsdp = (qxsk - qxskp1)/dp
            				pl  = p[k+1] - qxskp1/dqxsdp    ;pressure level of actual lcl.
            				dsdp   = (smix[k]  - smix[k+1])/dp
					dqtdp  = (qtmix[k] - qtmix[k+1])/dp
            				slcl   = smix[k+1]  +  dsdp* (pl-p[k+1])
            				qtlcl  = qtmix[k+1] +  dqtdp*(pl-p[k+1])
            				tfguess = tmix[k]
					ientropy,slcl,pl,qtlcl,tl,qslcl,tfguess
;       			        print*,' p',p(i,k+1),pl(i),p(i,lcl(i))
;       			        print*,' t',tmix(i,k+1),tl(i),tmix(i,lcl(i))
;       				print*,' s',smix(i,k+1),slcl,smix(i,lcl(i))
;        				print*,'qt',qtmix(i,k+1),qtlcl,qtmix(i,lcl(i))
;        				print*,'qs',qsmix(i,k+1),qslcl,qsmix(i,lcl(i))

      				ENDIF
;
      			ENDIF ;  k < klaunch
	ENDFOR ; Columns loop
;==========================================================================
;!!!!!!!!!!!!!!!!!!!!!!!!!!END ENTRAINMENT LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!

;; Could stop now and test with this as it will provide some estimate of buoyancy
;; without the effects of freezing/condensation taken into account for tmix.

;; So we now have a profile of entropy and total water of the entraining parcel
;; Varying with height from the launch level klaunch parcel=environment. To the 
;; top allowed level for the existence of convection.

;; Now we have to adjust these values such that the water held in vaopor is < or 
;; = to qsmix. Therefore, we assume that the cloud holds a certain amount of
;; condensate (lwmax) and the rest is rained out (xsh2o). This, obviously 
;; provides latent heating to the mixed parcel and so this has to be added back 
;; to it. But does this also increase qsmix as well? Also freezing processes


       FOR k = nlev-1,limconv,-1 DO BEGIN
        	xsh2o[k] = 0.
        	ds_xsh2o[k] = 0.
        	ds_freeze[k] = 0.
       ENDFOR

;!!!!!!!!!!!!!!!!!!!!!!!!!PRECIPITATION/FREEZING LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!
;!! Iterate solution twice for accuracy

        FOR k = nlev-1,limconv,-1 DO BEGIN
;          FOR k = nlev-1, 5,-1 DO BEGIN

; Initialize variables at k=klaunch
      		  IF (k EQ klaunch) THEN BEGIN

; Set parcel values at launch level assume not liquid water.            
         		tp[k]  = tmix[k]
         		qstp[k]  = q[k]
         		tpv[k]   =  (tp[k] + tpert) * (1.+1.608*qstp[k]) / (1.+qstp[k])
		  ENDIF
		  IF (k LT klaunch) THEN BEGIN

; Initiaite loop if switch(2) = .T. - RBN:DILUTE - TAKEN OUT BUT COULD BE RETURNED LATER.

; Iterate nit_lheat times for s,qt changes.

         	  FOR ii=0,nit_lheat-1 DO BEGIN

;Rain (xsh2o) is excess condensate, bar LWMAX (Accumulated loss from qtmix).
			test_max1=fltarr(2)
        		test_max1[0]=0.
        		test_max1[1]=(qtmix[k] - qsmix[k] - lwmax)
        		xsh2o[k] = max(test_max1)

;Contribution to ds from precip loss of condensate (Accumulated change from smix).(-ve)
                        test_max2=fltarr(2)
                        test_max2[0]=0.
                        test_max2[1]=(xsh2o[k]-xsh2o[k+1])                     
			ds_xsh2o[k] = ds_xsh2o[k+1] - cpliq*alog(tmix[k]/tfreez) * max(test_max2)

; Entropy of freezing: latice times amount of water involved divided by T.


          		IF (tmix[k] LE (tfreez+tscool)) AND (ds_freeze[k+1] EQ 0.) THEN BEGIN  ;One off freezing of condensate. 
				test_max3=fltarr(2)
                        	test_max3[0]=0.
                        	test_max3[1]=(qtmix[k]-qsmix[k]-xsh2o[k])
				ds_freeze[k] = (latice/tmix[k]) * max(test_max3)   ;Gain of LH
          		ENDIF

         		IF (tmix[k] LE (tfreez+tscool)) AND (abs(ds_freeze[k+1]) GT 0.) THEN BEGIN ; Continual freezing of additional condensate.
		          	test_max4=fltarr(2)
                                test_max4[0]=0.
                                test_max4[1]=(qsmix[k+1]-qsmix[k])
				ds_freeze[k] = ds_freeze[k+1]+(latice/tmix[k]) * max(test_max4)
          		ENDIF

; Adjust entropy and accordingly to sum of ds (be careful of signs).

			new_s = smix[k] + ds_xsh2o[k] + ds_freeze[k]

; Adjust liquid water and accordingly to xsh2o.

            		new_q = qtmix[k] - xsh2o[k]

;Invert entropy to get updated Tmix and qsmix of parcel.
            		tfguess = tmix[k]
                        TMIX_IN = tmix[k]
                        QSMIX_IN= qsmix[k]
			ientropy, new_s, p[k], new_q, TMIX_IN, QSMIX_IN, tfguess
                        tmix[k]=TMIX_IN
                        qsmix[k]=QSMIX_IN

		ENDFOR ; Iteration loop for freezing processes.

; tp  - Parcel temp is temp of mixture.
; tpv - Parcel v. temp should be density temp with new_q total water. 

		tp[k] = tmix[k]
; tpv = tprho in the presence of condensate (i.e. when new_q > qsmix)

         	IF(new_q GT qsmix[k]) THEN BEGIN  ; Super-saturated so condensate present - reduces buoyancy.
            		qstp[k] = qsmix[k]
         	ENDIF ELSE BEGIN                    ; Just saturated/sub-saturated - no condensate virtual effects.
            		qstp[k] = new_q
         	ENDELSE
		tpv[k] = (tp[k]+tpert)*(1.+1.608*qstp[k]) / (1.+ new_q)
      		ENDIF  ; k < klaunch
;print,  tp[i,k], qstp[i,k]

   ENDFOR ; Loop for vertical levels.
RETURN

END
