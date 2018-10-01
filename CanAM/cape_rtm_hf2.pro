;computes reversible CAPE 

PRO CAPE, filename1, filename2, filename3, filename4, thf, qhf, gzhf, phf
COMPILE_OPT IDL2

;=======================================================================================
month = ['m06', 'm07' , 'm08']

    FOR mm0=0,2 DO BEGIN

	filename1 = '/HOME/rtm/Desktop/CanAM4_hires/data/sa_rtm_hf2_1997_'+ month[mm0] +'_thf_gs.nc.001'
	Id1 = NCDF_OPEN(filename1)
	NCDF_VARGET, Id1,  'THF',    thf ; Array[128, 64, 49, 2976] or [lon, lat, level, time]

	;use lon and lat to match spCAM5 area, and levels from 16 (5000 Pa) or 20000 m 
	thf0 = thf[53:61, 32:35, 16:48, *] ; TWP
	;thf0 = thf[91:96, 42:48, 16:48, *] ; SGP
	thf=0.

	filename2 = '/HOME/rtm/Desktop/CanAM4_hires/data/sa_rtm_hf2_1997_'+ month[mm0] +'_qhf_gs.nc.001'
	Id2 = NCDF_OPEN(filename2)
	NCDF_VARGET, Id2,  'QHF',    qhf
	qhf0=qhf[53:61, 32:35, 16:48, *]
	;qhf0=qhf[91:96, 42:48, 16:48, *]
	qhf=0.

	filename3 = '/HOME/rtm/Desktop/CanAM4_hires/data/sa_rtm_hf2_1997_'+ month[mm0] +'_gzhf_gs.nc.001'
	Id3 = NCDF_OPEN(filename3)
	NCDF_VARGET, Id3,  'GZHF',   gzhf
	gzhf0=gzhf[53:61, 32:35, 16:48, *]
	;gzhf0=gzhf[91:96, 42:48, 16:48, *]
	gzhf=0.

	;phf is pressure in Pa and starts from top of the atmosphere
	
	filename4 = '/HOME/rtm/Desktop/CanAM4_hires/data/sa_rtm_hf2_1997_'+ month[mm0] +'_phf_gs.nc.001'
	Id4 = NCDF_OPEN(filename4)
	NCDF_VARGET, Id4,  'PHF',    phf 
	phf0=phf[53:61, 32:35, 16:48, *]
	;phf0=phf[91:96, 42:48, 16:48, *]
	phf=0.

        fname = '/HOME/rtm/Desktop/CanAM4_hires/data/TWP_sa_rtm_hf2_1997_'+ month[mm0] +'rCAPE_buoyancy+0K_and_rho_from_Tv_gs.nc.001'
        func1_rtm_hf2_rcape, thf0, qhf0, gzhf0, phf0, fname

   ENDFOR
;=======================================================================================

END
