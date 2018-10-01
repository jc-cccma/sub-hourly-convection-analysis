;look for convective rainfall within GCM grid box

;A CRM-column generates convective rainfall if W>1m/s or W<-1m/s (downdrafts)  at least at one level
;and cloud water + cloud ice mixing ratio is > 0.1 g/kg at the same level

PRO convective_PCP, filename, CRM_PREC_LON_150e_to_170e_LAT_0n_to_10n,CRM_W_LON_150e_to_170e_LAT_0n_to_10n,CRM_QC_LON_150e_to_170e_LAT_0n_to_10n,CRM_QI_LON_150e_to_170e_LAT_0n_to_10n

COMPILE_OPT IDL2
;========================================================
nlon=9
nlat=6
ntime=144
nlev=28
count_ntime0 = 0

crm_prec = FLTARR(nlon,nlat,32,ntime*85)
crm_w    = FLTARR(nlon,nlat,32,nlev, ntime*85)
crm_qc   = FLTARR(nlon,nlat,32,nlev, ntime*85)
crm_qi   = FLTARR(nlon,nlat,32,nlev, ntime*85)
CONV_PCP = FLTARR(nlon,nlat,ntime*85)
TOTAL_PCP= FLTARR(nlon,nlat,ntime*85)

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

        filename = '/HOME/rtm/Desktop/spCAM5/run-spcam_1.9x2.5_convection.cam.h2.0001-' + MONTHS[MONTH] +  '-' + DAYS[FILE_N] +'-00600.nc'

        ;Open the file for reading
        Id = NCDF_OPEN(filename)

        ;Reading the GFDL Files
        NCDF_VARGET, Id, 'CRM_PREC_LON_150e_to_170e_LAT_0n_to_10n', crm_prec0
        NCDF_VARGET, Id, 'CRM_W_LON_150e_to_170e_LAT_0n_to_10n',    crm_w0
        NCDF_VARGET, Id, 'CRM_QC_LON_150e_to_170e_LAT_0n_to_10n',   crm_qc0
        NCDF_VARGET, Id, 'CRM_QI_LON_150e_to_170e_LAT_0n_to_10n',   crm_qi0

    FOR ntime0=0,143 DO BEGIN
      crm_prec[*,*,*,count_ntime0]   = crm_prec0[*,*,*,0,ntime0]
      crm_w   [*,*,*,*,count_ntime0] = crm_w0  [*,*,*,0,*,ntime0]
      crm_qc  [*,*,*,*,count_ntime0] = crm_qc0 [*,*,*,0,*,ntime0]
      crm_qi  [*,*,*,*,count_ntime0] = crm_qi0 [*,*,*,0,*,ntime0]
      count_ntime0 = count_ntime0 + 1
    ENDFOR

  ENDFOR
ENDFOR

;========================================================

FOR LON1=0,nlon-1 DO BEGIN
 FOR LAT1=0,nlat-1 DO BEGIN
  FOR TIME1=0,count_ntime0-1 DO BEGIN
    COUNT1=0.
    PCP1=0.
    FOR CRM_POINT=0,31 DO BEGIN
      CHECK = 0	
      FOR LEVEL1=0,27 DO BEGIN
        IF (CHECK EQ 0) AND (((CRM_W[LON1,LAT1,CRM_POINT,LEVEL1,TIME1]) GT 1.) OR ((CRM_W[LON1,LAT1,CRM_POINT,LEVEL1,TIME1]) LT -1.)) THEN BEGIN
          FOR LEVEL1=0,27 DO BEGIN
            IF (CHECK EQ 0) AND ((CRM_QC[LON1,LAT1,CRM_POINT,LEVEL1,TIME1]+CRM_QI[LON1,LAT1,CRM_POINT,LEVEL1,TIME1]) GT 1E-4) THEN BEGIN

					CHECK = 1	
					IF (crm_prec[LON1,LAT1,CRM_POINT,TIME1] GE 1000.) THEN STOP
					IF (MAX(CRM_W[LON1,LAT1,CRM_POINT,*,TIME1]) GT 1000.) THEN STOP
					IF (MIN(CRM_W[LON1,LAT1,CRM_POINT,*,TIME1]) LT -1000.) THEN STOP
					IF (MAX(CRM_QC[LON1,LAT1,CRM_POINT,*,TIME1]) GT 1000.) THEN STOP
					IF (MAX(CRM_QI[LON1,LAT1,CRM_POINT,*,TIME1]) GT 1000.) THEN STOP
					
					COUNT1=COUNT1+1
					PCP1 = PCP1 + crm_prec[LON1,LAT1,CRM_POINT,TIME1]*3600.*1000.
	    ENDIF ;def2
	  ENDFOR ;def2
	ENDIF
      ENDFOR ; LEVELS
    ENDFOR
    IF (COUNT1 GT 0) THEN BEGIN
	CONV_PCP[LON1,LAT1,TIME1] = PCP1/32.
    ENDIF
    TOTAL_PCP[LON1,LAT1,TIME1] = TOTAL(crm_prec[LON1,LAT1,*,TIME1]*3600.*1000.)/32.
  ENDFOR
 ENDFOR
ENDFOR
;========================================================

file_name = 'convective_rain_in_updrafts_and_downdrafts_TWP_Q_GT1E-4_h2_DEF2.1997.nc'     ;Name of file to create

CONV_RAIN = FLTARR (NLON,NLAT,144*85)
CONV_RAIN[*,*,*]  = CONV_PCP[*,*,*]
TOTAL_RAIN = FLTARR (NLON,NLAT,144*85)
TOTAL_RAIN[*,*,*] = TOTAL_PCP[*,*,*]

;      Create output file and write data

id   = NCDF_CREATE(file_name, CLOBBER = clobber)                        ;Create netCDF output file


aid = NCDF_DIMDEF(id, 'LON',       NLON)
bid = NCDF_DIMDEF(id, 'LAT',       NLAT)
did = NCDF_DIMDEF(id, 'TIME',      144*85)

vid = NCDF_VARDEF(id, 'CONV_RAIN',    [aid,bid,did],   /FLOAT)
vid = NCDF_VARDEF(id, 'TOTAL_RAIN',   [aid,bid,did],   /FLOAT)

NCDF_ATTPUT, id, 'CONV_RAIN', 'units', 'Convective rainfall within GCM gridbox (mm/h)'     ;Write PCP units attribute
NCDF_ATTPUT, id, 'TOTAL_RAIN','units', 'TOTAL rainfall within GCM gridbox (mm/h)'     ;Write PCP units attribute

NCDF_CONTROL, id, /ENDEF                                                ;Exit define mode

NCDF_VARPUT, id, 'CONV_RAIN',  CONV_RAIN
NCDF_VARPUT, id, 'TOTAL_RAIN', TOTAL_RAIN

;*******************************************************************************************************************

NCDF_CLOSE, id                                                          ;Close netCDF output file

;*******************************************************************************************************************

END
