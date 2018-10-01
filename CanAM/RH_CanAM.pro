;CanAM
;computes RH from temperature and spec. humidity profiles. 
;saves RH as nc file

PRO RH, filename1, filename2, filename3, filename4, thf, qhf, gzhf, phf
COMPILE_OPT IDL2

;=======================================================================================
rd = 287.04 ; dry air
rv = 461.5  ; water vapor
Lv = 2.5104E6
grav = 9.80665; m/s2
cp   = 1005.7 ; J/kg/K
pi = 3.1415926535897932

;specify month
month = 'm08'

filename1 = '/HOME/rtm/Desktop/CanAM4_hires/data/sa_rtm_hf2_1997_'+ month +'_thf_gs.nc.001'
Id1 = NCDF_OPEN(filename1)
NCDF_VARGET, Id1,  'THF',    thf ; Array[128, 64, 49, 2976] or [lon, lat, level, time]

;use lon and lat to match spCAM5 area, and levels from 16 (5000 Pa) or 20000 m 
thf0 = thf[53:61, 32:35, *, *] - 273.15 ; TWP
thf=0.

filename2 = '/HOME/rtm/Desktop/CanAM4_hires/data/sa_rtm_hf2_1997_'+ month +'_qhf_gs.nc.001'
Id2 = NCDF_OPEN(filename2)
NCDF_VARGET, Id2,  'QHF',    qhf
qhf0=qhf[53:61, 32:35, *, *]
qhf=0.

;phf is pressure in Pa and starts from top of the atmosphere
filename4 = '/HOME/rtm/Desktop/CanAM4_hires/data/sa_rtm_hf2_1997_'+ month +'_phf_gs.nc.001'
Id4 = NCDF_OPEN(filename4)
NCDF_VARGET, Id4,  'PHF',    phf 
phf0=phf[53:61, 32:35, *, *]/100.  ; to mb
phf=0.

test_size = size(thf0)
nlon = test_size[1]
nlat = test_size[2]
nlev = test_size[3]
ntim = test_size[4]

RH0 = FLTARR(nlon,nlat,nlev,ntim)
RH0[*,*,*,*] = 'NaN'

FOR lon0 = 0, nlon-1 DO BEGIN
  FOR lat0 = 0, nlat-1 DO BEGIN
    FOR lev0 = 0, nlev-1 DO BEGIN
      FOR tim0 = 0, ntim-1 DO BEGIN
	RH0[lon0,lat0,lev0,tim0] = 100.*qhf0[lon0,lat0,lev0,tim0]/(qsat(thf0[lon0,lat0,lev0,tim0], phf0[lon0,lat0,lev0,tim0]))
      ENDFOR
    ENDFOR
  ENDFOR
ENDFOR

PRINT, RH0[3,3,48, 100]

fname = '/HOME/rtm/Desktop/CanAM4_hires/data/TWP_sa_rtm_hf2_1997_'+ month +'_RH'
;=======================================================================================
;save into netcdf
file_name = fname +'.nc'     ;Name of file to create

RH = FLTARR (nlon,nlat,nlev,ntim)
RH[*,*,*,*]  = RH0[*,*,*,*]

id   = NCDF_CREATE(file_name, CLOBBER = clobber)                        ;Create netCDF output file

bid = NCDF_DIMDEF(id, 'Long',      nlon)
cid = NCDF_DIMDEF(id, 'Lat',       nlat)
did = NCDF_DIMDEF(id, 'level',     nlev)
fid = NCDF_DIMDEF(id, 'Time',      ntim)

vid = NCDF_VARDEF(id, 'RH',         [bid,cid,did,fid],   /FLOAT)

NCDF_ATTPUT, id, 'RH',    'units', '(%)'     ;Write PCP units attribute

NCDF_CONTROL, id, /ENDEF                                    ;Exit define mode

NCDF_VARPUT, id, 'RH',    RH

NCDF_CLOSE, id



END
