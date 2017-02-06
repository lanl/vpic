pro movie

; Read in my color table

rgb=fltarr(3,226)
usecolor=fltarr(3,226)
red=fltarr(226)
blue=fltarr(226)
green=fltarr(226)
openr, inunit, '/scratch1/bjax/pproc/color.tbl', /get_lun
readf, inunit, rgb
close, inunit
free_lun, inunit
red=rgb(0,*)
green=rgb(1,*)
blue=rgb(2,*)
tvlct,red,green,blue

; Use X-Window plot

set_plot,'ps'
device, /Color
!x.style=1
!y.style=1
!p.color = 1
!p.background = 0
!p.charsize=2.0
!p.multi=[0,1,3]

;  Create array with some other color tables

; First open file

;????????????????????????????????????????????????
 openr,unit3,'movie_phase_He.bin',/get_lun, /swap_if_big_endian
;openr,unit1,'movie_ey_bsrs.bin',/get_lun, /swap_if_big_endian
 openr,unit1,'movie_ey.bin',/get_lun, /swap_if_big_endian
 openr,unit2,'movie_ex.bin',/get_lun, /swap_if_big_endian

; nframes=6000
 nframes=22
 nvx=200
 nbinx=200
 xmax=6112  ; in Debye Length
 nx=8192
 vthe= 0.037
 dt = 0.0270616
 movie_interval=40

 dvx=20/float(nvx)
 ex_norm = sqrt(2.0*vthe^2)

print,' Frames=',nframes
print,' Vx Bins=',nvx
print,'  x Bins=',nbinx
print,' Box size =',xmax,' Debye'

; Now set Vx grid
 vx=dvx*(findgen(nvx)-nvx/2)+dvx

; Set x grid
xmax = xmax/1000.0
dxx=xmax/float(nbinx)
xx=dxx*(findgen(nbinx)+0.5)
xx2=xmax*findgen(nx)/(nx-1)

; Now attach to file as a direct access binary

struct3 = {fe:fltarr(nvx,nbinx)}
field3 = assoc(unit3,struct3)

struct1 = {fy:fltarr(nx)}
field1 = assoc(unit1,struct1)

struct2 = {ee:fltarr(nx)}
field2 = assoc(unit2,struct2)

; Declare memory
f1v_ele_ave=fltarr(nvx)
f1v0_ele_ave=fltarr(nvx)
f1v_ele=fltarr(nvx,nbinx)
ey=fltarr(nx)
ex=fltarr(nx)

; Save average profile at t=0

struct3=field3[1]
f1v_ele = abs(struct3.fe)
f1v_ele_ave(*)=0.
for i=0,nvx-1 do begin
    for j=0,nbinx-1 do begin
        f1v_ele_ave(i)=f1v_ele_ave(i)+f1v_ele(i,j)
    endfor
endfor
 f1v0_ele_ave=f1v_ele_ave/float(nbinx)

; Now loop through frames
;????????????????????????????????????????????????
 for iframe=1,1 do begin

; Time for this frame
 time=float(iframe)*float(movie_interval)*dt

; Read data

    struct3=field3[iframe]
    f1v_ele = abs(struct3.fe)
    struct1=field1[iframe]
    ey = struct1.fy
    struct2=field2[iframe]
    ex = struct2.ee

; Set xrange and title

 !p.title='t*!7x!X!Dpe!N='+strcompress(string(FORMAT='(f5.0)',time))

; Contour plots

 !p.linestyle=0
 plot,xx2,ey, $
    xtitle='!6x (10!U3!X!N !7k!X!DD!N!6)',$
    ytitle='!6E!Dy!X!N', $
    xrange=[0,xmax],yrange=[-1.2,1.2]

 !p.linestyle=0
 plot,xx2,ex, $
    xtitle='!6x (10!U3!X!N !7k!X!DD!N!6)',$
    ytitle='!6E!Dx!X!N', $
    xrange=[0,xmax],yrange=[-0.03,0.03]

    foo = max(f1v0_ele_ave)
    num=32
    eps=1.0e-6
    facmax=1.19
    facmin=-0.04
    tmp = alog(f1v_ele/foo+eps)
    f1 = min(tmp)
    tmp = tmp - f1
    tmax=max(tmp)*facmax
    tmin=facmin*tmax
    step=(tmax-tmin)/float(num)
    clevels=indgen(num)*step + tmin
    contour,transpose(tmp),xx,vx, $
       ytitle='!6V!Ue!X!N!Dx!X!N/!6V!Ue!X!N!Dth!X!N', $
       xtitle='!6x (10!U3!X!N !7k!X!DD!N!6)',$
       xrange=[0.,xmax], yrange=[-8,8], /fill, levels=clevels

endfor

; Close the binary file

close,unit1
close,unit2
close,unit3

end



