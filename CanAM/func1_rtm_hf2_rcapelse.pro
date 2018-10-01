PRO func1_rtm_hf2_rcapelse, Tfld, qfld, zfld, pfld, divtfld, divqfld, qrlfld, qrsfld, fname

test_size = size(Tfld)
nlon = test_size[1]
nlat = test_size[2] 
nlev = test_size[3]
ntim = test_size[4]

;========================================================
rd   = 287.04
grav = 9.80616
capefld = FLTARR(nlon, nlat, ntim,3)
;==========================================================
FOR i=0,nlon-1 DO BEGIN
  FOR j=0,nlat-1 DO BEGIN
    FOR l=0,ntim-2 DO BEGIN
        T =     FLTARR(nlev)  
        q =     FLTARR(nlev)
        cape =  FLTARR(3)
        p =     FLTARR(nlev)
        zf =    FLTARR(nlev+1)            ;pressure of levels between standard pressure levels
        z =     FLTARR(nlev)
        ;z0 =    0.

; give values for T, Q, P for one time step for all levels
; and calculate pf (pressure for inter levels) and z (height in m)
; P is in mb and Pf in mb
;=========================================================
        T[*]=Tfld[i,j,*,l]
        q[*]=qfld[i,j,*,l]
        p[*]=pfld[i,j,*,l]
        z[*]=zfld[i,j,*,l]

        dt=900. ; time step
	FOR k=0,nlev-3 DO BEGIN 
;        FOR k=nlev-2, nlev-1 DO BEGIN
           T[k]=T[k]+dt*(divtfld[i,j,k,l+1]+qrsfld[i,j,k,l+1]+qrlfld[i,j,k,l+1])
           q[k]=q[k]+dt*(divqfld[i,j,k,l+1])
           IF (q[k] LE 0.) THEN BEGIN
		q[k]=1.e-10
	   ENDIF
        ENDFOR

        zf[0]=z[0] - (z[1]-z[0])*0.5           ;pressure (half-level) 
        FOR k=1,nlev-1 DO BEGIN
          zf[k]=(z[k]+z[k-1])*0.5               ;pressure between standard pressure levels
        ENDFOR
        zf[nlev]=z[nlev-1]+(z[nlev-1] - z[nlev-2])*0.5
        if (zf[nlev] LT 0.) then begin
                zf[nlev] = 0.
        endif

        p=p*0.01				;change from Pa to hPa

        ;Value has to start from top and end at the surface 
;        T = REVERSE(T)
;        q = REVERSE(q)
;        p = REVERSE(p)
;        z = REVERSE(z)
;        pf= REVERSE(pf)

;call buoyant_dilute subroutine and calculate CAPE
;print, 'cape before',cape       
        buoyan_dilute_rcape, q, T, p, z, zf, cape
;print, 'cape after', cape
;==================================================
        capefld[i,j,l,*]=cape[*]
;==================================================
    ENDFOR
  ENDFOR
ENDFOR
;==========================================================================================================
;save into netcdf
file_name = fname +'.nc'     ;Name of file to create

CAPElse = FLTARR (nlon,nlat,ntim,3)
CAPElse[*,*,*,*]  = capefld[*,*,*,*]

id   = NCDF_CREATE(file_name, CLOBBER = clobber)                        ;Create netCDF output file

bid = NCDF_DIMDEF(id, 'Long',      nlon)
cid = NCDF_DIMDEF(id, 'Lat',       nlat)
fid = NCDF_DIMDEF(id, 'Time',      ntim)
gid = NCDF_DIMDEF(id, 'values',     3)

vid = NCDF_VARDEF(id, 'CAPElse',	  [bid,cid,fid,gid],   /FLOAT)

NCDF_ATTPUT, id, 'CAPElse',    'units', 'CAPElse (J/kg)'     ;Write PCP units attribute

NCDF_CONTROL, id, /ENDEF                                    ;Exit define mode

NCDF_VARPUT, id, 'CAPElse',    CAPElse

NCDF_CLOSE, id 

capefld = 0.

RETURN

END