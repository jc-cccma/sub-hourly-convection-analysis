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
	
	filename5 = '/HOME/rtm/Desktop/CanAM4_hires/data/sa_rtm_hf2_1997_'+ month[mm0] +'_ttd_gs.nc.001'
	Id5 = NCDF_OPEN(filename5)
	NCDF_VARGET, Id5,  'TTD',    ttd
	ttd0=ttd[53:61, 32:35, 16:48, *]
	;ttd0=ttd[91:96, 42:48, 16:48, *]
	ttd=0.

	filename6 = '/HOME/rtm/Desktop/CanAM4_hires/data/sa_rtm_hf2_1997_'+ month[mm0] +'_estd_gs.nc.001'
	Id6 = NCDF_OPEN(filename6)
	NCDF_VARGET, Id6,  'ESTD',   estd
	estd0=estd[53:61, 32:35, 16:48, *]
	;estd0=estd[91:96, 42:48, 16:48, *]
	estd=0.

	filename7 = '/HOME/rtm/Desktop/CanAM4_hires/data/sa_rtm_hf2_1997_'+ month[mm0] +'_hlhf_gs.nc.001'
	Id7 = NCDF_OPEN(filename7)
	NCDF_VARGET, Id7,  'HLHF',    hlhf
	hlhf0=hlhf[53:61, 32:35, 16:48, *]
	;hlhf0=hlhf[91:96, 42:48, 16:48, *]
	hlhf=0.
	
	filename8 = '/HOME/rtm/Desktop/CanAM4_hires/data/sa_rtm_hf2_1997_'+ month[mm0] +'_hshf_gs.nc.001'
	Id8 = NCDF_OPEN(filename8)
	NCDF_VARGET, Id8,  'HSHF',    hshf
	hshf0=hshf[53:61, 32:35, 16:48, *]
	;hshf0=hshf[91:96, 42:48, 16:48, *]
	hshf=0.

        fname = '/HOME/rtm/Desktop/CanAM4_hires/data/TWP_sa_rtm_hf2_1997_'+ month[mm0] +'rCAPElse_buoyancy+0K_and_rho_from_Tv_gs.nc.001'
        func1_rtm_hf2_rcapelse, thf0, qhf0, gzhf0, phf0, ttd0, estd0, hlhf0, hshf0, fname

   ENDFOR

;=======================================================================================

END
