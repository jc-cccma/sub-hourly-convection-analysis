;spCAM5
;computes OMEGA [Pa/s] from U and V winds
;saves as .nc file

;======================================================================================
PRO omega,filename1, PRES0
COMPILE_OPT IDL2
;======================================================================================

NLON  = 9
NLAT  = 6
NTIME = 144
NLEV  = 30
count_ntime0 = 0

PRES   = FLTARR(NLON, NLAT, NLEV, NTIME*85)
PS     = FLTARR(NLON, NLAT, NTIME*85)
DPRES  = FLTARR(NLON, NLAT, NLEV, NTIME*85)
U      = FLTARR(NLON, NLAT, NLEV, NTIME*85)
V      = FLTARR(NLON, NLAT, NLEV, NTIME*85)

PRES1  = FLTARR(NLON, NLAT, NLEV, NTIME*85)
DPRES1 = FLTARR(NLON, NLAT, NLEV, NTIME*85)
U1     = FLTARR(NLON, NLAT, NLEV, NTIME*85)
V1     = FLTARR(NLON, NLAT, NLEV, NTIME*85)

DAYS = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31']
MONTHS = ['05','06','07']

FOR MONTH=0,2 DO BEGIN
        IF (MONTH EQ 0) THEN BEGIN
                STR=0
                N_DAYS=31
        ENDIF
        IF (MONTH EQ 1) THEN BEGIN
                STR=0
                N_DAYS=30
        ENDIF
        IF (MONTH EQ 2) THEN BEGIN
                STR=0
                N_DAYS=24
        ENDIF

        FOR FILE_N=STR,(N_DAYS)-1 DO BEGIN

        filename1 = '/HOME/rtm/Desktop/spCAM5/run-spcam_1.9x2.5_convection.cam.h2.0001-' + MONTHS[MONTH] +  '-' + DAYS[FILE_N] +'-00600.nc'

        ;Open the file for reading
        Id1 = NCDF_OPEN(filename1)

        NCDF_VARGET, Id1, 'DPRES_LON_150e_to_170e_LAT_0n_to_10n',        DPRES0
        NCDF_VARGET, Id1, 'PRES_LON_150e_to_170e_LAT_0n_to_10n',         PRES0
        NCDF_VARGET, Id1, 'PS_LON_150e_to_170e_LAT_0n_to_10n',           PS0
        NCDF_VARGET, Id1, 'LON_150e_to_170e',                            LON0
        NCDF_VARGET, Id1, 'LAT_0n_to_10n',                               LAT0
        NCDF_VARGET, Id1, 'U_LON_150e_to_170e_LAT_0n_to_10n',            U0
        NCDF_VARGET, Id1, 'V_LON_150e_to_170e_LAT_0n_to_10n',            V0
        NCDF_VARGET, Id1, 'date',            date
        NCDF_VARGET, Id1, 'ndbase',          ndbase
        NCDF_VARGET, Id1, 'nbdate',          nbdate

;print, DPRES0[3,3,29,100], DPRES0[3,3,28,100], DPRES0[3,3,27,100]
;print, (PS0[3,3,100] - PRES0[3,3,29,100])*2.
;print, (PRES0[3,3,29,100] - PRES0[3,3,28,100] - 1504.33/2.)*2. 

    FOR ntime0=0,143 DO BEGIN
      PRES1  [*,*,*,count_ntime0] = PRES0  [*,*,*,ntime0]
      PS     [*,*,count_ntime0]   = PS0    [*,*,ntime0]
      DPRES1 [*,*,*,count_ntime0] = DPRES0 [*,*,*,ntime0]
      U1     [*,*,*,count_ntime0] = U0     [*,*,*,ntime0]
      V1     [*,*,*,count_ntime0] = V0     [*,*,*,ntime0]
      count_ntime0 = count_ntime0 + 1
    ENDFOR

  ENDFOR
ENDFOR

;======================================================================================
; REVERSE LEVELS IN GCM FIELDS 

FOR LON1=0,NLON-1 DO BEGIN
  FOR LAT1=0,NLAT-1 DO BEGIN
    FOR TIME1=0,count_ntime0-1 DO BEGIN

      dP     = FLTARR(30)
      u_wind = FLTARR(30)
      v_wind = FLTARR(30)

      dP[*] = PRES1[LON1,LAT1,*,TIME1]
      PRES[LON1,LAT1,*,TIME1] = REVERSE(dP)

      dP[*] = DPRES1[LON1,LAT1,*,TIME1]
      DPRES[LON1,LAT1,*,TIME1] = REVERSE(dP)

      u_wind[*] = U1[LON1,LAT1,*,TIME1]
      U[LON1,LAT1,*,TIME1] = REVERSE(u_wind) ; compared GCM U and V winds to CRM U winds. GCM_V=MEAN(CRM U) , likely mistake when running model?

      v_wind[*] = V1[LON1,LAT1,*,TIME1]
      V[LON1,LAT1,*,TIME1] = REVERSE(v_wind)

    ENDFOR
  ENDFOR
ENDFOR
;======================================================================================
;CALCULATE OMEGA

OMEGA = FLTARR(NLON, NLAT, NLEV, NTIME*85)
OMEGA[*,*,*,*] = 1.E30
 
FOR LON1=1,NLON-2 DO BEGIN
  FOR LAT1=1,NLAT-2 DO BEGIN
    FOR TIME1=0,count_ntime0-1 DO BEGIN
;      FOR LEVEL1=0,29 DO BEGIN
      FOR LEVEL1=29,0,-1 DO BEGIN
        DISTANCE_X = 0.
        DISTANCE_Y = 0.

        DISTANCE_X=MAP_2POINTS(LON0[lon1],LAT0[lat1],LON0[lon1+1],LAT0[lat1], /METERS)
        DISTANCE_Y=MAP_2POINTS(LON0[lon1],LAT0[lat1],LON0[lon1],LAT0[lat1+1], /METERS)
;print, distance_x, distance_y

        A_COMP=0.
        B_COMP=0.
        A_COMP=((U[LON1+1,LAT1,LEVEL1,TIME1]-U[LON1-1,LAT1,LEVEL1,TIME1])/((2)*DISTANCE_X))
        B_COMP=((V[LON1,LAT1+1,LEVEL1,TIME1]-V[LON1,LAT1-1,LEVEL1,TIME1])/((2)*DISTANCE_Y))

;        IF (LEVEL1 EQ 0) THEN BEGIN
	 IF (LEVEL1 EQ 29) THEN BEGIN
          OMEGA[LON1,LAT1,LEVEL1,TIME1] = (A_COMP+B_COMP)*(DPRES[LON1,LAT1,LEVEL1,TIME1])*(-1.)
        ENDIF

;        IF (LEVEL1 GE 1) THEN BEGIN
	 IF (LEVEL1 NE 29) THEN BEGIN
;          OMEGA[LON1,LAT1,LEVEL1,TIME1] = OMEGA[LON1,LAT1,LEVEL1-1,TIME1] + ((A_COMP+B_COMP)*DPRES[LON1,LAT1,LEVEL1,TIME1])
          OMEGA[LON1,LAT1,LEVEL1,TIME1] = OMEGA[LON1,LAT1,LEVEL1+1,TIME1] + ((A_COMP+B_COMP)*(DPRES[LON1,LAT1,LEVEL1,TIME1]))*(-1.)
        ENDIF
 
     ENDFOR
    ENDFOR
  ENDFOR
ENDFOR

;======================================================================================
;Specify Filename

file_name = 'GCM_OMEGA_TWP_start_from_top.nc'     ;Name of file to create

;     Compute dependent variables

GCM_OMEGA = FLTARR (NLON,NLAT,30,144*85)
GCM_OMEGA[*,*,*,*]  = OMEGA[*,*,*,*]

;      Create output file and write data

id   = NCDF_CREATE(file_name, CLOBBER = clobber)                        ;Create netCDF output file

aid = NCDF_DIMDEF(id, 'LON',       NLON)
bid = NCDF_DIMDEF(id, 'LAT',       NLAT)
cid = NCDF_DIMDEF(id, 'LEVEL',     30)
did = NCDF_DIMDEF(id, 'TIME',      144*85)

vid = NCDF_VARDEF(id, 'GCM_OMEGA',    [aid,bid,cid,did],   /FLOAT)

NCDF_ATTPUT, id, 'GCM_OMEGA', 'units', 'GCM OMEGA (levels reversed and start from bottom -> up) computed from U and V winds starting from top (Pa/s)'     ;Write PCP units attribute

NCDF_CONTROL, id, /ENDEF                                                ;Exit define mode

NCDF_VARPUT, id, 'GCM_OMEGA',  GCM_OMEGA

;======================================================================================

END
