;+
; NAME:
;RPSCLOSE -- close ps, open X.
;     
; PURPOSE:
;       To close the Postscript device and set the graphics output
;       device back to X Windows.
;     
; CALLING SEQUENCE:
;       RPSCLOSE [, x0=x0, y0=y0, x1=x1, y1=y1, $
;              dy0=dy0, dy1=dy1, dx0=dx0, dx1=dx1, /high,
;		/skipsize]
;     
; INPUTS:
;       /high: set for highest quality PDFs
;       /skipsize: use to delete line defining paper size
;  Bounding Box Keywords
;       x0: x position of lower left hand corner
;       y0: y position of lower left hand corner
;       x1: x position of upper right hand corner
;       y1: y position of upper_right hand corner
;       dx0: offset for x position of lower left hand corner
;       dy0: offset for y position of lower left hand corner
;       dx1: offest for x position of upper right hand corner
;       dy1: offest for y position of upper_right hand corner
;     
; OUTPUTS:
;       None.
;
; KEYWORDS:
;       None.
;
; COMMON BLOCKS:
;       None.
;
; SIDE EFFECTS:
;       The device is changed.
;
; RESTRICTIONS:
;       A PostScript file must be open.
;
; RELATED PROCEDURES:
;       RPSOPEN
;
; MODIFICATION HISTORY:
;       Written by Tim Robishaw in ancient times.
;	Edited to include ps_mod to edit the PS file, Rds
;		1-20-11
;		
;-

pro rpsclose, _extra=_extra, thick=thick
common rpsopen_files, filename1

; MAKE SURE WE HAVE POSTSCRIPT DEVICE OPEN...
if (!d.name ne 'PS') then begin
    message, 'DEVICE is not set to PS!', /INFO
    return
endif

; CLOSE THE POSTSCRIPT DEVICE...
device, /CLOSE_FILE

; SET THE GRAPHICS OUTPUT DEVICE TO X WINDOWS...
set_plot, 'X'

if keyword_set(_extra) then ps_mod, filename1, _extra=_extra

delvarx, filename1
delvarx, rpsopen_files
end; psclose

