;
; Batch file for making 2d Ey plots from 2d vpic-3 runs. 
;
; This file is intended to be run from command line:  idl plot2d_ey
; The output filename is ey.ps
; dialog_pickfile is called from read_vpic_data, so one still needs
; to provide some input into the process. 
;
.run read_vpic_data
d = load_field_array(0)
eydata = reform(d.data.ey(*,0,*))
set_plot, 'ps'
device, filename='ey.ps'
contour, eydata, d.xmesh, d.zmesh, /fill, nlevels=20
device, /close
