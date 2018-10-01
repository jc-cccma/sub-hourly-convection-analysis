;Jan 10, 2018, code computes CAPE reversible, and dCAPElse 

PRO cape_calc, filename, T_LON_150e_to_170e_LAT_0n_to_10n, Q_LON_150e_to_170e_LAT_0n_to_10n, PRES_LON_150e_to_170e_LAT_0n_to_10n, PHIS_LON_150e_to_170e_LAT_0n_to_10n, DPRES_LON_150e_to_170e_LAT_0n_to_10n, SPDT_LON_150e_to_170e_LAT_0n_to_10n, SPQRS_LON_150e_to_170e_LAT_0n_to_10n, SPQRL_LON_150e_to_170e_LAT_0n_to_10n, SPQTLS_LON_150e_to_170e_LAT_0n_to_10n, SPTLS_LON_150e_to_170e_LAT_0n_to_10n, SPDQ_LON_150e_to_170e_LAT_0n_to_10n

COMPILE_OPT IDL2

;========================================================
nlev=30
ntime=144
nlat=6
nlon=9
count_ntime0 = 0

Tfld     = FLTARR(nlon,nlat,nlev, ntime*85)
qfld     = FLTARR(nlon,nlat,nlev, ntime*85)
pfld     = FLTARR(nlon,nlat,nlev, ntime*85)
phisfld  = FLTARR(nlon,nlat,ntime*85)
dpfld    = FLTARR(nlon,nlat,nlev, ntime*85)
spdtfld  = FLTARR(nlon,nlat,nlev, ntime*85)
qrsfld   = FLTARR(nlon,nlat,nlev, ntime*85)
qrlfld   = FLTARR(nlon,nlat,nlev, ntime*85)
divqfld  = FLTARR(nlon,nlat,nlev, ntime*85)
divtfld  = FLTARR(nlon,nlat,nlev, ntime*85)
spdqfld  = FLTARR(nlon,nlat,nlev, ntime*85)

DAYS = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31']
MONTHS = ['05','06','07']

FOR MONTH=0,2 DO BEGIN

  IF (MONTH EQ 0) THEN BEGIN
    N_DAYS=31
  ENDIF
  IF (MONTH EQ 1) THEN BEGIN
    N_DAYS=30
  ENDIF
  IF (MONTH EQ 2) THEN BEGIN
    N_DAYS=24       ; only 24 days in July for TWP
  ENDIF

  FOR FILE_N=0,(N_DAYS)-1 DO BEGIN

    filename = '/home/rtm/spCAM5/run-spcam_1.9x2.5_convection.cam.h2.0001-' + MONTHS[MONTH] +  '-' + DAYS[FILE_N] +'-00600.nc'
;    print, filename
    ;Open the file for reading

    Id = NCDF_OPEN(filename)

    ;Reading the GFDL Files

    NCDF_VARGET, Id, 'T_LON_150e_to_170e_LAT_0n_to_10n',            Tfld0
    NCDF_VARGET, Id, 'Q_LON_150e_to_170e_LAT_0n_to_10n',            qfld0
    NCDF_VARGET, Id, 'PRES_LON_150e_to_170e_LAT_0n_to_10n',         pfld0
    NCDF_VARGET, Id, 'PHIS_LON_150e_to_170e_LAT_0n_to_10n',         phisfld0
    NCDF_VARGET, Id, 'DPRES_LON_150e_to_170e_LAT_0n_to_10n',        dpfld0
    NCDF_VARGET, Id, 'SPDT_LON_150e_to_170e_LAT_0n_to_10n',         spdtfld0
    NCDF_VARGET, Id, 'SPQRS_LON_150e_to_170e_LAT_0n_to_10n',        qrsfld0
    NCDF_VARGET, Id, 'SPQRL_LON_150e_to_170e_LAT_0n_to_10n',        qrlfld0
    NCDF_VARGET, Id, 'SPQTLS_LON_150e_to_170e_LAT_0n_to_10n',       divqfld0
    NCDF_VARGET, Id, 'SPTLS_LON_150e_to_170e_LAT_0n_to_10n',        divtfld0
    NCDF_VARGET, Id, 'SPDQ_LON_150e_to_170e_LAT_0n_to_10n',         spdqfld0

    FOR ntime0=0,143 DO BEGIN
      Tfld    [*,*,*,count_ntime0] = Tfld0   [*,*,*,ntime0]
      qfld    [*,*,*,count_ntime0] = qfld0   [*,*,*,ntime0]
      pfld    [*,*,*,count_ntime0] = pfld0   [*,*,*,ntime0]
      phisfld [*,*,count_ntime0]   = phisfld0[*,*,ntime0]
      dpfld   [*,*,*,count_ntime0] = dpfld0  [*,*,*,ntime0]
      spdtfld [*,*,*,count_ntime0] = spdtfld0[*,*,*,ntime0]
      qrsfld  [*,*,*,count_ntime0] = qrsfld0 [*,*,*,ntime0]
      qrlfld  [*,*,*,count_ntime0] = qrlfld0 [*,*,*,ntime0]
      divqfld [*,*,*,count_ntime0] = divqfld0[*,*,*,ntime0]
      divtfld [*,*,*,count_ntime0] = divtfld0[*,*,*,ntime0]
      spdqfld [*,*,*,count_ntime0] = spdqfld0[*,*,*,ntime0]
      count_ntime0 = count_ntime0 + 1
    ENDFOR

  ENDFOR
ENDFOR

capefld =       FLTARR(nlon, nlat, ntime*85, 3)
dcapelse =      FLTARR(nlon, nlat, ntime*85, 3)
dcapelsp =      FLTARR(nlon, nlat, ntime*85, 3)

print, count_ntime0

;======================================================================================

msg=2 		 ;ic number of missing moisture levels at the top of model.
limconv=0        ;convection limiting level
;==========================================================

rd = 287.04
grav = 9.80616

FOR i=0,nlon-1 DO BEGIN
  FOR j=0,nlat-1 DO BEGIN
print, j
    FOR l=0,count_ntime0-2 DO BEGIN

       T =     FLTARR(nlev)
       q =     FLTARR(nlev)
       cape =  FLTARR(3)
       p =     FLTARR(nlev)
       zf =    FLTARR(nlev+1)            ;pressure of levels between standard pressure levels
       datp =  0.
       datv =  0.
       dalel = 0.
       divt =  FLTARR(nlev)
       divq =  FLTARR(nlev)
       qrl =   FLTARR(nlev)
       qrs =   FLTARR(nlev)
       cape1 = FLTARR(3)
       cape2 = FLTARR(3)
       tls =   FLTARR(nlev)
       qls =   FLTARR(nlev)
       Tm1 =   FLTARR(nlev)
       qm1 =   FLTARR(nlev)
       z =     FLTARR(nlev)

; give values for T, Q, P for one time step for all levels
; and calculate pf (pressure for inter levels) and z (height in m)
; P is in mb and Pf in mb
;=========================================================
            T[*]=Tfld[i,j,*,l]
            q[*]=qfld[i,j,*,l]
            p[*]=pfld[i,j,*,l]

; calculate height in meters
;===========================
          z[nlev-1]=phisfld[i,j,l]/grav                                ;sfc height

          FOR k=nlev-2,0,-1 DO BEGIN
                z[k] = z[k+1] - (rd/grav*(T[k+1]+T[k])*0.5*(1.+0.608*(q[k+1]+q[k])*0.5)*(alog(p[k])-alog(p[k+1])))
          ENDFOR

            zf[0]=z[0] - (z[1]-z[0])*0.5           ;pressure (half-level) 
            FOR k=1,nlev-1 DO BEGIN
              zf[k]=(z[k]+z[k-1])*0.5               ;pressure between standard pressure levels
            ENDFOR
            zf[nlev]=z[nlev-1]+(z[nlev-1] - z[nlev-2])*0.5
            if (zf[nlev] LT 0.) then begin
                zf[nlev] = 0.
            endif
;print, z
;stop
;print,''
;print, zf
;stop
            p=p*0.01				;change from Pa to hPa
;==================================================     
;print, 'cape before',cape	 
;if (i EQ 4) and (j EQ 4) and (l EQ 284) then begin
	buoyan_dilute_cape, q, T, p, z, zf, cape, msg;, mse_max_level
;endif
;print, 'cape after', cape
;==================================================
        capefld[i,j,l,*]=cape[*]
;==================================================
        dt=600.
        FOR k=0,nlev-3 DO BEGIN
;TM++ sept 22 , instead forcing at time t add forcing at t+1, balance check of temperature
           T[k]=T[k]+dt*(divtfld[i,j,k,l+1]+qrsfld[i,j,k,l+1]+qrlfld[i,j,k,l+1])
           q[k]=q[k]+dt*(divqfld[i,j,k,l+1])
           IF (q[k] LE 0.) THEN BEGIN
		q[k]=1.e-10
	   ENDIF
        ENDFOR
;cape1 is the cape after free tropospheric large-scale forcing is applied
;goto 20
;print, 'cape 1 before', cape1
        buoyan_dilute_cape,q,T,p,z, zf,cape1, msg;, mse_max_level
        dcapelse[i,j,l,*]=cape1[*]-cape[*]
;print, 'cape 1 after', cape1
;==================================================
;reset T and q fields
T=0.
q=0.
       T =     FLTARR(nlev)
       q =     FLTARR(nlev)

                T[*]=Tfld[i,j,*,l]
                q[*]=qfld[i,j,*,l]
;print, 'q',q

            FOR k=nlev-2, nlev-1 DO BEGIN
		T[k]=T[k]+dt*(divtfld[i,j,k,l+1]+qrsfld[i,j,k,l+1]+qrlfld[i,j,k,l+1])
           	q[k]=q[k]+dt*(divqfld[i,j,k,l+1])
            IF (q[k] LE 0.) THEN BEGIN
		 q[k]=1.e-10
	    ENDIF
;print, dt*(divqfld[i,j,k,l])
            ENDFOR
;print, 'q after', q
;stop
;       cape2 is the cape after PBL large-scale forcing is applied

;print, 'cape 2 before', cape2
        buoyan_dilute_cape,q,T,p,z, zf,cape2, msg;, mse_max_level
        dcapelsp[i,j,l,*]=cape2[*]-cape[*]

;==================================================
	ENDFOR
    ENDFOR
ENDFOR

PRINT, MAX(capefld)

;******************************** SAVE RESULTS TO NCDF FORMAT ********************
;Specify Filename

file_name = 'reversible_CAPE_h2_1997_buoyancy+0K_and_rho_from_Tv.nc'     ;Name of file to create
;     Compute dependent variables

CAPE_capefld = FLTARR (nlon,nlat,count_ntime0,3)
CAPE_capefld[*,*,*,*]  = capefld[*,*,*,*]
CAPE_dcapelse = FLTARR (nlon,nlat,count_ntime0,3)
CAPE_dcapelse[*,*,*,*]  = dcapelse[*,*,*,*]
CAPE_dcapelsp = FLTARR (nlon,nlat,count_ntime0,3)
CAPE_dcapelsp[*,*,*,*]  = dcapelsp[*,*,*,*]
;      Create output file and write data

id   = NCDF_CREATE(file_name, CLOBBER = clobber)                        ;Create netCDF output file


bid = NCDF_DIMDEF(id, 'LON',       nlon)
cid = NCDF_DIMDEF(id, 'LAT',       nlat)
fid = NCDF_DIMDEF(id, 'TIME',      count_ntime0)
gid = NCDF_DIMDEF(id, 'VARIABLE',  3)

vid = NCDF_VARDEF(id, 'CAPE_capefld',	  [bid,cid,fid,gid],   /FLOAT)
vid = NCDF_VARDEF(id, 'CAPE_dcapelse',    [bid,cid,fid,gid],   /FLOAT)
vid = NCDF_VARDEF(id, 'CAPE_dcapelsp',    [bid,cid,fid,gid],   /FLOAT)

NCDF_ATTPUT, id, 'CAPE_capefld',    'units', 'CAPE (J/kg)'     ;Write PCP units attribute
NCDF_ATTPUT, id, 'CAPE_dcapelse',   'units', 'CAPE (J/kg)'     ;Write PCP units attribute
NCDF_ATTPUT, id, 'CAPE_dcapelsp',   'units', 'CAPE (J/kg)'     ;Write PCP units attribute

NCDF_CONTROL, id, /ENDEF                                    ;Exit define mode

NCDF_VARPUT, id, 'CAPE_capefld',    CAPE_capefld
NCDF_VARPUT, id, 'CAPE_dcapelse',   CAPE_dcapelse
NCDF_VARPUT, id, 'CAPE_dcapelsp',   CAPE_dcapelsp

NCDF_CLOSE, id                                                          ;Close netCDF output file

;*******************************************************************************************************************
;closing
;    NCDF_CLOSE, Id

PRINT, 'DONE..........'

END
