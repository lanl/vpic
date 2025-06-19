pro read_slice,ifield
  common pdata,fulldata,simulation_time,mytitle
  common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
  common choice, rawdata,computedquantity
  common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
  common picdata, field, struct

    fulldata(*,*) = 0.0
    print,'v.time for read =',v.time
    for iave = -v.tave,v.tave do begin
      print,"Reading data --> ",rawdata(ifield)," time slice=",v.time+iave
      field = assoc(ifield,struct)
      struct = field[v.time + iave]
      fulldata(*,*) = fulldata(*,*) + struct.data(*,*,v.ycut)
   endfor
    fulldata = fulldata/(2*v.tave+1)
    fulldata = smooth(fulldata,v.smoothing,/EDGE_TRUNCATE)
end

pro gyrcan
  common pdata,fulldata,simulation_time,mytitle
  common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
  common choice, rawdata,computedquantity
  common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
  common picdata, field, struct

print,'--------------------------------'
print," Time Index=",v.time
print," Gyroviscous cancel. "
mytitle="Gyroviscous cancel."

; Define variables

Pxx = fltarr(nx,nz)
Pyy = fltarr(nx,nz)
Pzz = fltarr(nx,nz)
Ptot = fltarr(nx,nz)
Bx = fltarr(nx,nz)
By = fltarr(nx,nz)
Bz = fltarr(nx,nz)
Uex = fltarr(nx,nz)
Uey = fltarr(nx,nz)
Uez = fltarr(nx,nz)
vpar = fltarr(nx,nz)
dxPtot = fltarr(nx,nz)
dzPtot = fltarr(nx,nz)

dxvpar = fltarr(nx,nz)
dzvpar = fltarr(nx,nz)

; Read in data
read_slice,ipe
Pxx = fulldata
read_slice,ipe+3
Pyy = fulldata
read_slice,ipe+5
Pzz = fulldata

read_slice,ibx
Bx = fulldata
read_slice,iby
By = fulldata
read_slice,ibz
Bz = fulldata

read_slice,iuex
Uex = fulldata
read_slice,iuey
Uey = fulldata
read_slice,iuez
Uez = fulldata

; Calculate total pressure
Ptot = (Pxx + Pyy + Pzz)/3.0

; Calculate ve_parallel
vpar = (Uex*Bx+Uey*By+Uez*Bz)/sqrt(Bx*Bx+By*By+Bz*Bz)

; Calculate derivatives
dx=xmax/nx
dz=zmax/nz

for i=1,nx-2 do begin
   for j=1,nz-2 do begin
      dxPtot(i,j) = (Ptot(i+1,j)-Ptot(i-1,j))/(2.*dx)
      dzPtot(i,j) = (Ptot(i,j+1)-Ptot(i,j-1))/(2.*dz)
      dxvpar(i,j) = (vpar(i+1,j)-vpar(i-1,j))/(2.*dx)
      dzvpar(i,j) = (vpar(i,j+1)-vpar(i,j-1))/(2.*dz)
;         fulldata(i,j) =
;         fulldata(i,j)+(Pyz(i,j+1)-Pyz(i,j-1))/(2.*dz)                                                                                       
   end
end

fulldata = (dxPtot*dzvpar-dzPtot*dxvpar)/By

end

pro gyroterms
  common pdata,fulldata,simulation_time,mytitle
  common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
  common choice, rawdata,computedquantity
  common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
  common picdata, field, struct

print,'--------------------------------'
print," Time Index=",v.time
print," Other gyroviscous terms "
mytitle="Some Gyroviscous terms."

; Define variables                                                                                                                                            

Pxx = fltarr(nx,nz)
Pyy = fltarr(nx,nz)
Pzz = fltarr(nx,nz)
Ptot = fltarr(nx,nz)
Bx = fltarr(nx,nz)
By = fltarr(nx,nz)
Bz = fltarr(nx,nz)
Btot = fltarr(nx,nz)
ByBtot = fltarr(nx,nz)
Uex = fltarr(nx,nz)
Uey = fltarr(nx,nz)
Uez = fltarr(nx,nz)
vpar = fltarr(nx,nz)

dxBtot = fltarr(nx,nz)
dzBtot = fltarr(nx,nz)
dxvpar = fltarr(nx,nz)
dzvpar = fltarr(nx,nz)
dxByBt = fltarr(nx,nz)
dzByBt = fltarr(nx,nz)

dxvex = fltarr(nx,nz)
dzvex = fltarr(nx,nz)
dxvez = fltarr(nx,nz)
dzvez = fltarr(nx,nz)
Ax = fltarr(nx,nz)
Az = fltarr(nx,nz)
dzAx = fltarr(nx,nz)
dxAz = fltarr(nx,nz)

dzBx = fltarr(nx,nz)
dxBz = fltarr(nx,nz)

print,dzBx

ident = fltarr(nx,nz)

gyro1 = fltarr(nx,nz)
gyro2 = fltarr(nx,nz)
gyro3 = fltarr(nx,nz)

; Read in data                                                                                                                                                
read_slice,ipe
Pxx = fulldata
read_slice,ipe+3
Pyy = fulldata
read_slice,ipe+5
Pzz = fulldata

read_slice,ibx
Bx = fulldata
read_slice,iby
By = fulldata
read_slice,ibz
Bz = fulldata

read_slice,iuex
Uex = fulldata
read_slice,iuey
Uey = fulldata
read_slice,iuez
Uez = fulldata

; Calculate total B

Btot = sqrt(Bx*Bx+By*By+Bz*Bz)

; ByovBtot

ByBtot = By/Btot

; Calculate total pressure                                                                                                                                    
Ptot = (Pxx + Pyy + Pzz)/3.0

; Calculate ve_parallel                                                                                                                                       
vpar = (Uex*Bx+Uey*By+Uez*Bz)/sqrt(Bx*Bx+By*By+Bz*Bz)

; Calculate derivatives                                                                                                                                       
dx=xmax/nx
dz=zmax/nz

print,dx,dz

for i=1,nx-2 do begin
   for j=1,nz-2 do begin
      dxBtot(i,j) = (Btot(i+1,j)-Btot(i-1,j))/(2.*dx)
      dzBtot(i,j) = (Btot(i,j+1)-Btot(i,j-1))/(2.*dz)
      dxvpar(i,j) = (vpar(i+1,j)-vpar(i-1,j))/(2.*dx)
      dzvpar(i,j) = (vpar(i,j+1)-vpar(i,j-1))/(2.*dz)
      dxByBt(i,j) = (ByBtot(i+1,j)-ByBtot(i-1,j))/(2.*dx)
      dzByBt(i,j) = (ByBtot(i,j+1)-ByBtot(i,j-1))/(2.*dz)
      
      dxvex(i,j) = (Uex(i+1,j)-Uex(i-1,j))/(2.*dx)
      dzvex(i,j) = (Uex(i,j+1)-Uex(i,j-1))/(2.*dz)
      dxvez(i,j) = (Uez(i+1,j)-Uez(i-1,j))/(2.*dx)
      dzvez(i,j) = (Uez(i,j+1)-Uez(i,j-1))/(2.*dz)

      dxBz(i,j)  = (Bz(i+1,j)-Bz(i-1,j))/(2.*dx)
      dzBx(i,j)  = (Bx(i,j+1)-Bx(i,j-1))/(2.*dz)
   end
end


Ax = (Bx*dxvex + Bz*dzvex)/Btot
Az = (Bx*dxvez + Bz*dzvez)/Btot

for i=1,nx-2 do begin
   for j=1,nz-2 do begin
      dzAx(i,j) = (Ax(i,j+1)-Ax(i,j-1))/(2.*dz)
      dxAz(i,j) = (Az(i+1,j)-Az(i-1,j))/(2.*dx)
   end
end

gyro1 = (Ptot/(Btot*Btot))*(dxBtot*dzvpar - dzBtot*dxvpar)
gyro2 = (Ptot/Btot)*(dzByBt*dxvpar-dxByBt*dzvpar)
gyro3 = (Ptot/Btot)*(dzAx - dxAz)
gyro4 = (Ptot/(Btot*Btot))*(dzBx*dxvex-dxBz*dzvez)

;ident = ident+1.0

print,'About to plot'
;print,dxBz
fulldata = -gyro1

end

pro EHesse
  common pdata,fulldata,simulation_time,mytitle
  common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
  common choice, rawdata,computedquantity
  common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
  common picdata, field, struct

print,'--------------------------------'
print," Time Index=",v.time
print," Hesse electric field. "
mytitle="Hesse field"


; Allocate for arrays

  Uex = fltarr(nx,nz)
  Uey = fltarr(nx,nz)
  Uez = fltarr(nx,nz)
  vpar = fltarr(nx,nz)
  Pxx = fltarr(nx,nz)
  Pyy = fltarr(nx,nz)
  Pzz = fltarr(nx,nz)
  Ptot = fltarr(nx,nz)
  Bx = fltarr(nx,nz)
  By = fltarr(nx,nz)
  Bz = fltarr(nx,nz)
  Btot = fltarr(nx,nz)
  dxUex = fltarr(nx,nz)
  dxxvpar = fltarr(nx,nz)
  dzzvpar = fltarr(nx,nz)

; Load in data

  read_slice,ipe
  Pxx = fulldata
  read_slice,ipe+3
  Pyy = fulldata
  read_slice,ipe+5
  Pzz = fulldata
  
  read_slice,ibx
  Bx = fulldata
  read_slice,iby
  By = fulldata
  read_slice,ibz
  Bz = fulldata
  
  read_slice,iuex
  Uex = fulldata
  read_slice,iuey
  Uey = fulldata
  read_slice,iuez
  Uez = fulldata

; Calculate total B
  Btot = sqrt(Bx*Bx+By*By+Bz*Bz)

; Calculate total pressure                                                                                                                                    
  Ptot = (Pxx + Pyy + Pzz)/3.0

; Calculate ve_parallel                                                                                                                                       
  vpar = (Uex*Bx+Uey*By+Uez*Bz)/sqrt(Bx*Bx+By*By+Bz*Bz)

; Calculate derivatives                                                                                                                                       
  dx=xmax/nx
  dz=zmax/nz

  print,dx,dz

; First derivatives
  for i=1,nx-2 do begin
     for j=1,nz-2 do begin
        dxUex(i,j) = (Uex(i+1,j)-Uex(i-1,j))/(2.*dx)
     end
  end

; Second derivative
  for i=1,nx-2 do begin
     for j=1,nz-2 do begin
        dxxvpar(i,j) = (vpar(i+1,j) - 2.0*vpar(i,j) + vpar(i-1,j))/(dx*dx)
        dzzvpar(i,j) = (vpar(i,j+1) - 2.0*vpar(i,j) + vpar(i,j-1))/(dz*dz)
     end
  end

fulldata = -0.5*(Ptot/(Btot*Btot))*dxUex*(dxxvpar+dzzvpar)


end


pro EUexB_y
  common pdata,fulldata,simulation_time,mytitle
  common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
  common choice, rawdata,computedquantity
  common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
  common picdata, field, struct

print,'--------------------------------'
print," Time Index=",v.time
print," EUexB_y "
mytitle="EUexB_y"

dne = fltarr(nx,nz)
ey  = fltarr(nx,nz)
Uex = fltarr(nx,nz)
Uez = fltarr(nx,nz)
Bx  = fltarr(nx,nz)
Bz  = fltarr(nx,nz)

read_slice,ine
dne = fulldata

read_slice,iey
ey = fulldata

read_slice,iuex
Uex = fulldata
read_slice,iuez
Uez = fulldata

read_slice,ibx
Bx = fulldata
read_slice,ibz
Bz = fulldata

fulldata = -dne*(ey + Uez*Bx-Uex*Bz)

end

pro dtnUey
  common pdata,fulldata,simulation_time,mytitle
  common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
  common choice, rawdata,computedquantity
  common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
  common picdata, field, struct

print,'--------------------------------'
print," Time Index=",v.time
print," dtnUey "
mytitle="dtnUey"

print, 'v.tave=',v.tave
t1 = v.tave

if (t1 eq 0) then t1=1
print,'t1=',t1

snapwci = 0.025 
wpewce = 8.
mime = 25.
dtwpe = snapwci*mime*wpewce

dnet0 = fltarr(nx,nz)
dnet1 = fltarr(nx,nz)

Ueyt0 = fltarr(nx,nz)
Ueyt1 = fltarr(nx,nz)

read_slice,ine
dnet0 = fulldata

read_slice,iuey
Ueyt0 = fulldata

; Go forward in time
v.time = v.time + t1
print,'Set time = ',v.time

read_slice,ine
dnet1 = fulldata

read_slice,iuey
Ueyt1 = fulldata

; Reset time
v.time = v.time - t1
print,'Reset time = ',v.time


fulldata = (dnet1*Ueyt1-dnet0*Ueyt0)/(t1*dtwpe)

end

pro divnUeUey
  common pdata,fulldata,simulation_time,mytitle
  common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
  common choice, rawdata,computedquantity
  common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
  common picdata, field, struct

print,'--------------------------------'
print," Time Index=",v.time
print," divnUeUey "
mytitle="divnUeUey"

; Declare needed variables

dne = fltarr(nx,nz)
Uex = fltarr(nx,nz)
Uey = fltarr(nx,nz)
Uez = fltarr(nx,nz)

dxUey = fltarr(nx,nz)
dzUey = fltarr(nx,nz)

read_slice,ine
dne = fulldata

read_slice,iuex
Uex = fulldata
read_slice,iuey
Uey = fulldata
read_slice,iuez
Uez = fulldata

dx=xmax/nx
dz=zmax/nz
   
for i=1,nx-2 do begin
   for j=1,nz-2 do begin
      dxUey(i,j) = (Uey(i+1,j)-Uey(i-1,j))/(2.*dx)
      dzUey(i,j) = (Uey(i,j+1)-Uey(i,j-1))/(2.*dz)
;         fulldata(i,j) = fulldata(i,j)+(Pyz(i,j+1)-Pyz(i,j-1))/(2.*dz)
   end
end

fulldata = dne*(Uex*dxUey+Uez*dzUey)

end

pro divPey
common pdata,fulldata,simulation_time,mytitle
common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct

print,'--------------------------------'
print," Time Index=",v.time
print," divPey "
mytitle="divPey"

; Declare needed variables

Pxy = fltarr(nx,nz)
Pyz = fltarr(nx,nz)

dxPxy = fltarr(nx,nz)
dzPyz = fltarr(nx,nz)

; Read in needed data

read_slice,ipe+1
Pxy = fulldata

read_slice,ipe+4
Pyz = fulldata

dx=xmax/nx
dz=zmax/nz
   
for i=1,nx-2 do begin
   for j=1,nz-2 do begin
      dxPxy(i,j) = (Pxy(i+1,j)-Pxy(i-1,j))/(2.*dx)
      dzPyz(i,j) = (Pyz(i,j+1)-Pyz(i,j-1))/(2.*dz)
   end
end

fulldata = dxPxy + dzPyz

end

pro EdotJ
common pdata,fulldata,simulation_time,mytitle
common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct

; Declare needed variables

jx = fltarr(nx,nz)
jy = fltarr(nx,nz)
jz = fltarr(nx,nz)

Ex = fltarr(nx,nz)
Ey = fltarr(nx,nz)
Ez = fltarr(nx,nz)


; Read needed data for contour plot

print,'--------------------------------'
print," Time Index=",v.time
print," E*J "
mytitle="E*J"

read_slice,ijx
jx =fulldata
read_slice,ijy
jy =fulldata
read_slice,ijz
jz =fulldata

read_slice,iex
Ex =fulldata
read_slice,iey
Ey =fulldata
read_slice,iez
Ez =fulldata

; Return E*J

fulldata = Ex*jx + Ey*jy + Ez*jz
fulldata = smooth(fulldata,v.smoothing,/EDGE_TRUNCATE)

end


pro beta,isp
common pdata,fulldata,simulation_time,mytitle
common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct

if isp eq 0 then begin
   print, "Computing Electron Beta"
   mytitle="Beta-e"
   ip=ipe
endif
if isp eq 1 then begin
   print, "Computing Ion Beta"
   mytitle="Beta-i"
   ip=ipi
endif
if isp eq 2 then begin
   print, "Computing Total Beta"
   mytitle="Beta-Total"
   ip=ipe
endif

; Declare needed variables

Pxx = fltarr(nx,nz)
Pyy = fltarr(nx,nz)
Pzz = fltarr(nx,nz)
Pixx = fltarr(nx,nz)
Piyy = fltarr(nx,nz)
Pizz = fltarr(nx,nz)

bx = fltarr(nx,nz)
by = fltarr(nx,nz)
bz = fltarr(nx,nz)

bm = fltarr(nx,nz)

; Read needed data for contour plot

print,'--------------------------------'
print,"Time Index=",v.time

; Read MagneticField

read_slice,ibx
bx = fulldata
read_slice,iby
by = fulldata
read_slice,ibz
bz = fulldata

; Now read in pressure the data

read_slice,ip
Pxx = fulldata
read_slice,ip+3
Pyy = fulldata
read_slice,ip+5
Pzz = fulldata

if isp eq 2 then begin
   read_slice,ipi
   Pxx = Pxx + fulldata
   read_slice,ipi+3
   Pyy = Pyy + fulldata
   read_slice,ipi+5
   Pzz = Pzz + fulldata
endif

; Return beta

fulldata = 2.0*(Pxx+Pyy+Pzz)/(bx^2+by^2+bz^2)/3.0

end

pro temperature,isp
  common pdata,fulldata,simulation_time,mytitle
  common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
  common choice, rawdata,computedquantity
  common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
  common picdata, field, struct

  if isp eq 0 then begin
     print, "Computing Electron Temperature"
     mytitle="Te"
     ip=ipe
  endif
  if isp eq 1 then begin
     print, "Computing Ion Temperature"
     mytitle="Ti"
     ip=ipi
  endif
  if isp eq 2 then begin
     print, "Computing Total Temperature"
     mytitle="T"
     ip=ipe
  endif

; Declare needed variables

  Pxx = fltarr(nx,nz)
  Pyy = fltarr(nx,nz)
  Pzz = fltarr(nx,nz)
  Pixx = fltarr(nx,nz)
  Piyy = fltarr(nx,nz)
  Pizz = fltarr(nx,nz)
  dne = fltarr(nx,nz)

; Read needed data for contour plot

  print,'--------------------------------'
  print,"Time Index=",v.time
  simulation_time = v.time

; Density

  read_slice,ine
  dne = fulldata

; Read Pressure

  read_slice,ip
  Pxx = fulldata
  read_slice,ip+3
  Pyy = fulldata
  read_slice,ip+5
  Pzz = fulldata

  if isp eq 2 then begin
     read_slice,ipi
     Pxx = Pxx + fulldata
     read_slice,ipi+3
     Pyy = Pyy + fulldata
     read_slice,ipi+5
     Pzz = Pzz + fulldata
endif

; Return Temperature

  fulldata = (Pxx+Pyy+Pzz)/(3.0*dne)

end

pro anisotropy,isp
common pdata,fulldata,simulation_time,mytitle
common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct

if isp eq 0 then begin
   print, "Computing Electron Anisotropy"
   mytitle="Ane"
   ip=ipe
endif
if isp eq 1 then begin
   print, "Computing Ion Anisotropy"
   mytitle="Ani"
   ip=ipi
endif

; Declare needed variables

Pxx = fltarr(nx,nz)
Pxy = fltarr(nx,nz)
Pxz = fltarr(nx,nz)
Pyy = fltarr(nx,nz)
Pyz = fltarr(nx,nz)
Pzz = fltarr(nx,nz)

ppar = fltarr(nx,nz)
pper = fltarr(nx,nz)

bx = fltarr(nx,nz)
by = fltarr(nx,nz)
bz = fltarr(nx,nz)

bm = fltarr(nx,nz)

; Read needed data for contour plot

print,'--------------------------------'
print,"Time Index=",v.time

read_slice,ibx
bx = fulldata
read_slice,iby
by = fulldata
read_slice,ibz
bz = fulldata

bm = sqrt(bx^2 + by^2 + bz^2)
bx = bx/bm
by = by/bm
bz = bz/bm

; Now read in the pressure tensor data

  read_slice,ip
  Pxx = fulldata
  read_slice,ip+1
  Pxy = fulldata
  read_slice,ip+2
  Pxz = fulldata
  read_slice,ip+3
  Pyy = fulldata
  read_slice,ip+4
  Pyz = fulldata
  read_slice,ip+5
  Pzz = fulldata

; Compute parallel and perpendicular pressure

ppar = Pxx*bx^2 + Pyy*by^2 + Pzz*bz^2 + 2.0*Pxy*bx*by + 2.0*Pxz*bx*bz+ 2.0*Pyz*by*bz
pper = (Pxx + Pyy +Pzz - ppar)/2.0

; Return data

fulldata = pper/ppar

end

pro agyrotropy,isp
common pdata,fulldata,simulation_time,mytitle
common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct

; Choose either electrons or ions

if isp eq 0 then begin
   print, "Computing Electron Agyrotropy"
   mytitle="A0e"
   ip=ipe
endif
if isp eq 1 then begin
   print, "Computing Ion Agyrotropy"
   mytitle="A0i"
   ip=ipi
endif

; Declare needed variables

Pxx = fltarr(nx,nz)
Pxy = fltarr(nx,nz)
Pxz = fltarr(nx,nz)
Pyy = fltarr(nx,nz)
Pyz = fltarr(nx,nz)
Pzz = fltarr(nx,nz)

bx = fltarr(nx,nz)
by = fltarr(nx,nz)
bz = fltarr(nx,nz)

bm = fltarr(nx,nz)

Nxx = fltarr(nx,nz)
Nxy = fltarr(nx,nz)
Nxz = fltarr(nx,nz)
Nyy = fltarr(nx,nz)
Nyz = fltarr(nx,nz)
Nzz = fltarr(nx,nz)
alpha = fltarr(nx,nz)
beta = fltarr(nx,nz)

; Read needed data for contour plot

print,'--------------------------------'
print," Time Index=",v.time
print," Smoothing =",v.smoothing

read_slice,ibx
bx = fulldata
read_slice,iby
by = fulldata
read_slice,ibz
bz = fulldata

bm = sqrt(bx^2 + by^2 + bz^2)
bx = bx/bm
by = by/bm
bz = bz/bm

; Now read in the pressure tensor data

  read_slice,ip
  Pxx = fulldata
  read_slice,ip+1
  Pxy = fulldata
  read_slice,ip+2
  Pxz = fulldata
  read_slice,ip+3
  Pyy = fulldata
  read_slice,ip+4
  Pyz = fulldata
  read_slice,ip+5
  Pzz = fulldata

; Compute parallel and perpendicular pressure

Nxx =  by*by*Pzz - 2.0*by*bz*Pyz + bz*bz*Pyy
Nxy = -by*bx*Pzz + by*bz*Pxz + bz*bx*Pyz - bz*bz*Pxy
Nxz =  by*bx*Pyz - by*by*Pxz - bz*bx*Pyy + bz*by*Pxy
Nyy =  bx*bx*Pzz - 2.0*bx*bz*Pxz + bz*bz*Pxx
Nyz = -bx*bx*Pyz + bx*by*Pxz + bz*bx*Pxy - bz*by*Pxx
Nzz =  bx*bx*Pyy - 2.0*bx*by*Pxy + by*by*Pxx

alpha = Nxx + Nyy + Nzz
beta = -(Nxy^2 + Nxz^2 + Nyz^2 - Nxx*Nyy - Nxx*Nzz - Nyy*Nzz)

; Return agyrotropy data

fulldata = 2.0*sqrt(alpha^2-4.0*beta)/alpha

end


pro firehose,nsp
common pdata,fulldata,simulation_time,mytitle
common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct

if nsp eq 0 then begin
   print, "Computing Electron Firehose"
   mytitle="E-Firehose"
endif
if nsp eq 1 then begin
   print, "Computing Total Firehose"
   mytitle="Firehose"
endif

; Declare needed variables

Pxx = fltarr(nx,nz)
Pxy = fltarr(nx,nz)
Pxz = fltarr(nx,nz)
Pyy = fltarr(nx,nz)
Pyz = fltarr(nx,nz)
Pzz = fltarr(nx,nz)

ppar = fltarr(nx,nz)
pper = fltarr(nx,nz)

bx = fltarr(nx,nz)
by = fltarr(nx,nz)
bz = fltarr(nx,nz)

bm = fltarr(nx,nz)

; Read needed data for contour plot

print,'--------------------------------'
print," Time Index=",v.time
print," Firehose condition"

print,"Reading data --> ",rawdata(ibx)
field = assoc(ibx,struct)
struct = field[v.time]
bx(*,*) = struct.data(*,v.ycut,*)
bx = smooth(bx,v.smoothing,/EDGE_TRUNCATE)
simulation_time = v.time

print,"Reading data --> ",rawdata(iby)
field = assoc(iby,struct)
struct = field[v.time]
by(*,*) = struct.data(*,v.ycut,*)
by = smooth(by,v.smoothing,/EDGE_TRUNCATE)

print,"Reading data --> ",rawdata(ibz)
field = assoc(ibz,struct)
struct = field[v.time]
bz(*,*) = struct.data(*,v.ycut,*)
bz = smooth(bz,v.smoothing,/EDGE_TRUNCATE)

bm = sqrt(bx^2 + by^2 + bz^2)
bx = bx/bm
by = by/bm
bz = bz/bm

; Now read in the data

for isp=0,nsp do begin
   if isp eq 0 then ip=ipe
   if isp eq 1 then ip=ipi
   print,"ip=",ip,ipe,ipi

   print,"Reading data --> ",rawdata(ip)
   field = assoc(ip,struct)
   struct = field[v.time]
   Pxx(*,*) = Pxx(*,*) + struct.data(*,v.ycut,*)

   print,"Reading data --> ",rawdata(ip+1)
   field = assoc(ip+1,struct)
   struct = field[v.time]
   Pxy(*,*) = Pxy(*,*) +  struct.data(*,v.ycut,*)

   print,"Reading data --> ",rawdata(ip+2)
   field = assoc(ip+2,struct)
   struct(*,*) = field[v.time]
   Pxz(*,*) =  Pxz(*,*) + struct.data(*,v.ycut,*)

   print,"Reading data --> ",rawdata(ip+3)
   field = assoc(ip+3,struct)
   struct(*,*) = field[v.time]
   Pyy(*,*) = Pyy(*,*) + struct.data(*,v.ycut,*)

   print,"Reading data --> ",rawdata(ip+4)
   field = assoc(ip+4,struct)
   struct = field[v.time]
   Pyz(*,*) = Pyz(*,*) + struct.data(*,v.ycut,*)

   print,"Reading data --> ",rawdata(ip+5)
   field = assoc(ip+5,struct)
   struct = field[v.time]
   Pzz(*,*) = Pzz(*,*) + struct.data(*,v.ycut,*)

endfor

; Smooth the total pressure

Pxx = smooth(Pxx,v.smoothing,/EDGE_TRUNCATE)
Pxy = smooth(Pxy,v.smoothing,/EDGE_TRUNCATE)
Pxz = smooth(Pxz,v.smoothing,/EDGE_TRUNCATE)
Pyy = smooth(Pyy,v.smoothing,/EDGE_TRUNCATE)
Pyz = smooth(Pyz,v.smoothing,/EDGE_TRUNCATE)
Pzz = smooth(Pzz,v.smoothing,/EDGE_TRUNCATE)

; Compute parallel and perpendicular pressure

ppar = Pxx*bx^2 + Pyy*by^2 + Pzz*bz^2 + 2.0*Pxy*bx*by + 2.0*Pxz*bx*bz+ 2.0*Pyz*by*bz
pper = (Pxx + Pyy + Pzz - ppar)/2.0

; Compute firehose condition - this is the actual phase velocity of the
;                              Alfven wave divided by Va (isotropic)

fulldata = 1+ pper/bm^2-ppar/bm^2

end

pro kelvin_helmholtz
common pdata,fulldata,simulation_time,mytitle
common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct

; Declare needed variables

ux = fltarr(nx,nz)
uy = fltarr(nx,nz)
uz = fltarr(nx,nz)
nee = fltarr(nx,nz)

bx = fltarr(nx,nz)
by = fltarr(nx,nz)
bz = fltarr(nx,nz)

; Read needed data for contour plot

print,'--------------------------------'
print," Time Index=",v.time
print," Kelvin Helmholtz condition"
mytitle="KH Stability"

print,"Reading data --> ",rawdata(ibx)
field = assoc(ibx,struct)
struct = field[v.time]
bx(*,*) = struct.data(*,v.ycut,*)
bx = smooth(bx,v.smoothing,/EDGE_TRUNCATE)
simulation_time = v.time

print,"Reading data --> ",rawdata(iby)
field = assoc(iby,struct)
struct = field[v.time]
by(*,*) = struct.data(*,v.ycut,*)
by = smooth(by,v.smoothing,/EDGE_TRUNCATE)

print,"Reading data --> ",rawdata(ibz)
field = assoc(ibz,struct)
struct = field[v.time]
bz(*,*) = struct.data(*,v.ycut,*)
bz = smooth(bz,v.smoothing,/EDGE_TRUNCATE)

print,"Reading data --> ",rawdata(iuix)
field = assoc(iuix,struct)
struct = field[v.time]
ux(*,*) = struct.data(*,v.ycut,*)
ux = smooth(ux,v.smoothing,/EDGE_TRUNCATE)

print,"Reading data --> ",rawdata(iuiy)
field = assoc(iuiy,struct)
struct = field[v.time]
uy(*,*) = struct.data(*,v.ycut,*)
uy = smooth(uy,v.smoothing,/EDGE_TRUNCATE)

print,"Reading data --> ",rawdata(iuiz)
field = assoc(iuiz,struct)
struct = field[v.time]
uz(*,*) = struct.data(*,v.ycut,*)
uz = smooth(uz,v.smoothing,/EDGE_TRUNCATE)

print,"Reading data --> ",rawdata(ine)
field = assoc(ine,struct)
struct = field[v.time]
nee(*,*) = struct.data(*,v.ycut,*)
nee = smooth(nee,v.smoothing,/EDGE_TRUNCATE)

; Comput KH ratio
;
;  ratio = (k*u)/(2k*Va)
;
;  where u is local fluid velocity
;        Va is local Alfven velocity
;        k is unity vector for KH 
;
;  Need this ratio to change by factor of >1 across a shear layer
;
;  *** Need to specify mass ratio ***

mime=100

;  Assume "k" is in the direction of "u"
;fulldata = (ux^2 + uy^2 + uz^2)*sqrt(nee*mime)/(ux*bx + uy*by + uz*bz)/2

;  Use angle slider to pick direction

pi=3.1415927
kx = cos(v.angle*pi/180)
ky = sin(v.angle*pi/180)
kz = 0.0
print,"k=",kx,ky,kz

; Use this direction to compute KH ratio

fulldata = ((kx*ux + ky*uy + kz*uz)*sqrt(nee*mime))/(2*abs(kx*bx + ky*by + kz*bz))

end

pro total_current
common pdata,fulldata,simulation_time,mytitle
common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct

; Declare needed variables

uix = fltarr(nx,nz)
uiy = fltarr(nx,nz)
uiz = fltarr(nx,nz)

ux = fltarr(nx,nz)
uy = fltarr(nx,nz)
uz = fltarr(nx,nz)
nee = fltarr(nx,nz)

; Read needed data for contour plot

print,'--------------------------------'
print," Time Index=",v.time
print," Total Current "
mytitle="|J|"

print,"Reading data --> ",rawdata(iuix)
field = assoc(iuix,struct)
struct = field[v.time]
ux(*,*) = struct.data(*,v.ycut,*)

print,"Reading data --> ",rawdata(iuiy)
field = assoc(iuiy,struct)
struct = field[v.time]
uy(*,*) = struct.data(*,v.ycut,*)

print,"Reading data --> ",rawdata(iuiz)
field = assoc(iuiz,struct)
struct = field[v.time]
uz(*,*) = struct.data(*,v.ycut,*)

print,"Reading data --> ",rawdata(iuex)
field = assoc(iuex,struct)
struct = field[v.time]
ux(*,*) = ux(*,*) - struct.data(*,v.ycut,*)

print,"Reading data --> ",rawdata(iuey)
field = assoc(iuiy,struct)
struct = field[v.time]
uy(*,*) = uy(*,*) - struct.data(*,v.ycut,*)

print,"Reading data --> ",rawdata(iuez)
field = assoc(iuiz,struct)
struct = field[v.time]
uz(*,*) = uz(*,*) - struct.data(*,v.ycut,*)

print,"Reading data --> ",rawdata(ine)
field = assoc(ine,struct)
struct = field[v.time]
nee(*,*) = struct.data(*,v.ycut,*)

; Total current

fulldata = nee*sqrt(ux^2+uy^2+uz^2)
fulldata = smooth(fulldata,v.smoothing,/EDGE_TRUNCATE)

end

pro shear_angle
common pdata,fulldata,simulation_time,mytitle
common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct

; Declare needed variables

bx = fltarr(nx,nz)
by = fltarr(nx,nz)
bz = fltarr(nx,nz)

; Read needed data for contour plot

print,'--------------------------------'
print," Time Index=",v.time
print," Computing Angle of Magnetic Shear"
mytitle="Magnetic Shear Angle"

print,"Reading data --> ",rawdata(ibx)
field = assoc(ibx,struct)
struct = field[v.time]
bx(*,*) = struct.data(*,v.ycut,*)
bx = smooth(bx,v.smoothing,/EDGE_TRUNCATE)
simulation_time = v.time

print,"Reading data --> ",rawdata(iby)
field = assoc(iby,struct)
struct = field[v.time]
by(*,*) = struct.data(*,v.ycut,*)
by = smooth(by,v.smoothing,/EDGE_TRUNCATE)

print,"Reading data --> ",rawdata(ibz)
field = assoc(ibz,struct)
struct = field[v.time]
bz(*,*) = struct.data(*,v.ycut,*)
bz = smooth(bz,v.smoothing,/EDGE_TRUNCATE)

;  Use angle slider to pick reference direction

pi=3.1415927
kx = cos(v.angle*pi/180)
ky = sin(v.angle*pi/180)
kz = 0.0
print,"k=",kx,ky,kz

; Use this direction to compute magnetic shear angle

fulldata = (180/pi)*acos((kx*bx + ky*by + kz*bz)/sqrt(bx^2+by^2+bz^2))

end


pro overplot_contours,xarr,zarr,lx,lz,imin,imax,jmin,jmax,outform
common choice,rawdata,computedquantity
common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common colortable,rgb,usecolor,red,blue,green,range1,range2,r1,r2,tmax,tmin
common picdata, field, struct
fulldata2 = fltarr(nx,nz)
temp2 = fltarr(lx,lz)	
print,"Overplot with = ",rawdata(v.overplot)
field = assoc(v.overplot,struct)
struct = field[v.time]
fulldata2 = struct.data

;print, "Adding contours",imin,imax,v.ycut,jmin,jmax
;temp2(*,*) = fulldata2(imin:imax,v.ycut,jmin:jmax) 
temp2(*,*) = fulldata2(imin:imax,jmin:jmax)

if (v.smoothing le 8) then temp2 = smooth(temp2,v.smoothing)
if (outform ne 2) then begin
    tmax=max(temp2)
    tmin=min(temp2)
    dr = (tmax-tmin)*0.20
    tmin=tmin-dr/2
    tmax=tmax+dr/2
endif

;temp2 = transpose(temp2)

step=(tmax-tmin)/v.contours
clevels=indgen(v.contours)*step + tmin       
contour,temp2,xarr,zarr,levels=clevels,/overplot,color=0
end

pro doplot,outform
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct
common colortable,rgb,usecolor,red,blue,green,range1,range2,r1,r2,tmax,tmin
common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common imagesize,nxpix,nypix,xoff,yoff,xpic,ypic
common pdata,fulldata,simulation_time,mytitle
common refresh,update

simulation_time = v.time

; Check and see if we need to read in new data

if ( update eq 1) then begin

;  If raw plot is selected - then use this.  Otherwise do the
;                            diagnostic plot

; set how many slices to average over

   if (v.rawplot ge 1) then begin
      read_slice,v.rawplot
;    fulldata(*,*) = 0.0
;    for iave = -v.tave,v.tave do begin
;      print,"Reading data --> ",rawdata(v.rawplot)," time slice=",v.time+iave
;      field = assoc(v.rawplot,struct)
;      struct = field[v.time + iave]
;      fulldata(*,*) = fulldata(*,*) + struct.data(*,v.ycut,*)
;   endfor
;    fulldata = fulldata/(2*v.tave+1)
;    fulldata = smooth(fulldata,v.smoothing,/EDGE_TRUNCATE)
    print,"Maximum=",max(fulldata)
    mytitle = rawdata(v.rawplot)
   endif else begin

; Compute the desired diagnostic

; Electron Anisotropy
      if v.computed eq 0 then anisotropy,0

; Ion Anisotropy
      if v.computed eq 1 then anisotropy,1

; Electron Agyrotropy
      if v.computed eq 2 then agyrotropy,0

; Ion Agyrotropy
      if v.computed eq 3 then agyrotropy,1

; Total Firehose condition
      if v.computed eq 4 then firehose,1

; Electron Firehose condition
      if v.computed eq 5 then firehose,0

; Simple KH condition
      if v.computed eq 6 then kelvin_helmholtz

; Magnetic Rotation Angle
      if v.computed eq 7 then shear_angle

; Magnetic Rotation Angle
      if v.computed eq 8 then total_current

; Beta-e
      if v.computed eq 9 then beta,0

; Beta-i
      if v.computed eq 10 then beta,1

; Beta-Total
      if v.computed eq 11 then beta,2

; EdotJ
      if v.computed eq 12 then EdotJ

; Electron temperature
      if v.computed eq 13 then temperature,0

; Ion temperature
      if v.computed eq 14 then temperature,1

; Total temperature
      if v.computed eq 15 then temperature,2

; y-cpt ideal Ohms
      if v.computed eq 16 then EUexB_y

; dt(nUey)
      if v.computed eq 17 then dtnUey

; div(nUeUey)
      if v.computed eq 18 then divnUeUey

; divPey
      if v.computed eq 19 then divPey

; Gyroviscous cancellation
      if v.computed eq 20 then gyrcan
; Gyroterms
      if v.computed eq 21 then gyroterms

; Hesse electric field
      if v.computed eq 22 then EHesse
   endelse

;  Now shift data the required amount in x direction

fulldata = shift(fulldata,v.shift,0)

endif

print,"Time slice=",v.time,"  Shifted in x by =",v.shift

; Select region of 2d-data to plot - This allows me to create contours just within the specified
; region and ignore the rest of the data

; Declare memory for reduced region and load region into plotting array

imin = fix(v.xmin/xmax*nx)
imax = fix(v.xmax/xmax*nx)-1
lx = imax-imin +1 
jmin = fix(v.zmin/zmax*nz)
jmax = fix(v.zmax/zmax*nz)-1
lz = jmax-jmin +1 

;print,v.zmax,zmax,imin,imax,jmin,jmax,lx,lz

temp = fltarr(lx,lz)
temp(*,*) = fulldata(imin:imax,jmin:jmax)

xarr = ((v.xmax-v.xmin)/lx)*(findgen(lx)+0.5) + v.xmin
zarr = v.zmin + ((v.zmax-v.zmin)/lz)*(findgen(lz)+0.5) -zmax/2

; Output info about data

print,'Maximum Value=',max(abs(temp))

; set output for either X11 or a postscript file

if ( outform eq 0 ) then begin
    set_plot, 'X'
    device,true_color=24,decomposed=0	
    !x.style=1
    !y.style=1
    !P.color =0
    !p.background=1
endif

if (outform eq 1 ) then begin
    set_plot, 'ps'
    !p.color =1
    !p.background = 0
    !p.font = 0	
    !x.style=1
    !y.style=1
    width=6.0
    asp=0.45
    height=width*asp
    if (v.rawplot ge 1) then fname = rawdata(v.rawplot) else fname=computedquantity(v.computed)
    file=strcompress(fname+string(v.time)+'.eps',/remove_all)
    print,"Postscript file --> ",file
    device, /encapsulated, bits_per_pixel=32, filename = file, $  
            /inches, xsize=width, ysize = height, /TIMES, font_size=12, /color, set_font='Times-Roman'
endif

;  Set colors so that red is largest value

if ( outform ne 2 ) then begin
    r1=min(temp);-0.05;min(temp)
    r2=max(temp);5.0;max(temp)
    dr = (r2-r1)*0.40
    r1=r1-dr/2
    r2=r2+dr/2
endif

; Use limited range for KH diagnostic

If (v.computed eq 6 and v.rawplot eq 0) then begin
   r1 = -2.0
   r2 = 2.0
endif

; Also firehose
;If ((v.computed eq 5) or (v.computed eq 4) and v.rawplot eq 0) then begin
;   r1 = -0.5
;   r2 = 1.0
;endif

; ********************************************************************************
; *********************** Contour Plot with Colorbar ****************************
; ********************************************************************************

if ( v.plottype eq 0 ) then begin

; Set position and titles

;    print, size(temp , /dimensions)

    xoff = 0.11
    yoff = 0.10
    xpic = 0.90
    ypic = 0.92
    dx1=0.012
    dx2=0.045
    !x.title=""
    !y.title=""
    !p.position=[xoff,yoff,xpic,ypic] 
    shade_surf, temp, xarr, zarr, ax=90,az=0,shades=bytscl(temp,max=r2,min=r1),zstyle=4,charsize=0.9,pixels=1000
;    xyouts,v.xmin+(v.xmax-v.xmin)/3.1,(v.zmax)*0.51,mytitle+"     t="+strcompress(string(FORMAT='(f8.1)',simulation_time))+"     Index="+strcompress(string(v.time),/remove_all),charsize=1.2

if ( v.overplot gt 0 ) then overplot_contours,xarr,zarr,lx,lz,imin,imax,jmin,jmax,outform

;  Now add the color bar

    !x.title=""
    !y.title=""
    colorbar, Position = [xpic+dx1,yoff,xpic+dx2,ypic], /Vertical, $
              Range=[r1,r2], Format='(f6.2)', /Right, $
              Bottom=5, ncolors=251, Divisions=6, font=9, charsize=0.8

endif

; **************************************************************************
; **********************  Contour Plot with X-average **********************
; **************************************************************************

if ( v.plottype eq 1 ) then begin

    xave = fltarr(lz)
    for j=0,lz-1 do begin
        xave(j) = 0.0
        for i=0,lx-1 do begin
            xave(j) = xave(j) +  temp(i,j)
        endfor
    endfor
    xave = xave/lx

    if ( outform eq 0) then begin
        range1=min(xave)
        range2=max(xave)
        dr = (range2-range1)*0.20
        range1=range1-dr/2
        range2=range2+dr/2
    endif

; Set position and titles

    xoff = 0.09
    yoff = 0.10
    xpic = 0.70
    ypic = 0.94

    !p.position=[xoff,yoff,xpic,ypic]
    !x.title="x "
    !y.title="z"
    !p.position=[xoff,yoff,xpic,ypic]     
    shade_surf, temp, xarr, zarr, ax=90,az=0,shades=bytscl(temp,max=r2,min=r1),zstyle=4,charsize=0.9,pixels=1000
    xyouts,v.xmin+(v.xmax-v.xmin)/3.1,(v.zmax)*0.51,mytitle+"      t="+strcompress(string(FORMAT='(f8.1)',simulation_time))+"     Index="+strcompress(string(v.time),/remove_all),charsize=0.9

if ( v.overplot gt 0 ) then overplot_contours,xarr,zarr,lx,lz,imin,imax,jmin,jmax,outform

    xoff2 = 0.705
    yoff2 = 0.10
    xpic2 = 0.94
    ypic2 = 0.94
    !x.title=" "
    !y.title=' '
    !p.position=[xoff2,yoff2,xpic2,ypic2] 
    plot,xave,zarr,charsize=0.8,yrange=[v.zmin-zmax/2,v.zmax-zmax/2],ytickname=replicate(' ',8),xrange=[range1,range2],/NOERASE

endif

; **************************************************************************
; **********************  Contour Plot with Z-Slice ************************
; **************************************************************************

if ( v.plottype eq 2 ) then begin

; Z-slice of data (fix x and vary z)

    slice = fltarr(lx)
    islice=fix((v.zslice-v.zmin)/(v.zmax-v.zmin)*float(lz)-0.5)
    print,"islice=",islice
    for j=0,lx-1 do begin
        slice(j) = temp(j,islice)
    endfor

    if ( outform eq 0) then begin
        range1=min(slice)
        range2=max(slice)
        dr = (range2-range1)*0.20
        range1=range1-dr/2
        range2=range2+dr/2
    endif

; Set position and titles
    
    xoff = 0.09
    yoff = 0.40
    xpic = 0.94
    ypic = 0.94

    !p.position=[xoff,yoff,xpic,ypic]
    !x.title=" "
    !y.title="z"
    !p.position=[xoff,yoff,xpic,ypic]     
    shade_surf, temp, xarr, zarr, ax=90,az=0,shades=bytscl(temp,max=r2,min=r1),zstyle=4,charsize=0.9,pixels=1000,xtickname=replicate(' ',8)
    xyouts,v.xmin+(v.xmax-v.xmin)/3.1,(v.zmax)*0.51,mytitle+"      t ="+strcompress(string(FORMAT='(f8.1)',simulation_time))+"     Index="+strcompress(string(v.time),/remove_all),charsize=0.9

if ( v.overplot gt 0 ) then overplot_contours,xarr,zarr,lx,lz,imin,imax,jmin,jmax,outform

; Overplot a line to show where the x-slice is located

    xs=fltarr(lx)
    for j=0,lx-1 do begin
        xs(j) = v.zslice-zmax/2
    endfor
    oplot,xarr,xs,color=0

; Add x-slice plot

    xoff2 = 0.09
    yoff2 = 0.10
    xpic2 = 0.94
    ypic2 = 0.395

    !x.title="x"
    !y.title=mytitle
    !p.position=[xoff2,yoff2,xpic2,ypic2] 
    plot,xarr,slice,charsize=0.9,xrange=[v.xmin,v.xmax],yrange=[range1,range2],/NOERASE
;    oplot,xarr,slice2,color=110
;    oplot,xarr,slice3,color=220

endif

; **************************************************************************
; **********************  Contour Plot with X-Slice ************************
; **************************************************************************

if ( v.plottype eq 3 ) then begin

; X-slice of data (fix z and vary x)

    jslice=fix((v.xslice-v.xmin)/(v.xmax-v.xmin)*float(lx)-0.5)
;    print,"jslice=",jslice,v.xslice,v.xmin,v,xmax,lx
    slice = fltarr(lz)
;    slice2 = fltarr(lz)
;    slice3 = fltarr(lz)
    for i=0,lz-1 do begin
        slice(i) = temp(jslice,i)
    endfor

    if ( outform eq 0) then begin
        range1=min(slice)
        range2=max(slice)
        dr = (range2-range1)*0.30
        range1=range1-dr
        range2=range2+dr*1.15
    endif


;    If ((v.computed eq 5) or (v.computed eq 4) and v.rawplot eq 0) then begin
;       range1 = -0.5
;       range2 = 1.0
;   endif

    If (v.computed eq 5 and v.rawplot eq 0) then begin
       range1 = -2.0
       range2 = 2.0
    endif

; Set position and titles

    xoff = 0.09
    yoff = 0.10
    xpic = 0.70
    ypic = 0.94

    !p.position=[xoff,yoff,xpic,ypic]
    !x.title="x"
    !y.title="z"
    !p.position=[xoff,yoff,xpic,ypic]     
    shade_surf, temp, xarr, zarr, ax=90,az=0,shades=bytscl(temp,max=r2,min=r1),zstyle=4,charsize=0.9
    xyouts,v.xmin+(v.xmax-v.xmin)/3.1,(v.zmax)*0.51,mytitle+"      t="+strcompress(string(FORMAT='(f8.1)',simulation_time))+"     Index="+strcompress(string(v.time),/remove_all),charsize=0.9

if ( v.overplot gt 0 ) then overplot_contours,xarr,zarr,lx,lz,imin,imax,jmin,jmax,outform

; Overplot a line to show where the z-slice is located

    xs=fltarr(lz)
    for i=0,lz-1 do begin
        xs(i) = v.xslice
    endfor
    oplot,xs,zarr,color=0

    If (v.computed eq 6 and v.rawplot eq 0) then begin
    range1 = -2.0
    range2 = 2.0
endif

; Add z-slice plot

    xoff2 = 0.705
    yoff2 = 0.10
    xpic2 = 0.94
    ypic2 = 0.94

    !x.title=" "
    !y.title=" "
    !p.position=[xoff2,yoff2,xpic2,ypic2] 
    plot,slice,zarr,charsize=0.9,yrange=[v.zmin-zmax/2,v.zmax-zmax/2],ytickname=replicate(' ',8),xrange=[range1,range2],/NOERASE
;    oplot,slice2,zarr,color=110
;    oplot,slice3,zarr,color=50

 endif

; Set refresh to false (i.e. won't read back in data unless needed)

update = 1

end


pro handle_event, ev
common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq, anames
common imagesize,nxpix,nypix,xoff,yoff,xpic,ypic
;common oldposition,xmax0,zmin0
common oldposition,xmin0,zmax0
common refresh,update
widget_control, ev.id, get_uval = whichbutton

case whichbutton of
    'done' : begin
        close,/All
        free_lun, 1
        widget_control, ev.top, /destroy
    end
    'plot' : begin
        widget_control,ev.top,get_uvalue=v

        v.zmin=0.
        v.zmax=zmax
        v.xmin=0.
        v.xmax=xmax
        widget_control,ev.top,set_uvalue=v    
        doplot,0

    end
    'render' : begin
        widget_control,ev.top,get_uvalue=v
        doplot,1
    end
    'animate' : begin
        widget_control,ev.top,get_uvalue=v
        create_movie
    end
    'region' : begin
        widget
    end

    'mouse' : begin
        widget_control,ev.top,get_uvalue=v
        if (ev.type eq 0) then begin
            xmin0=v.xmin
            zmax0=v.zmax
            v.xmin=(float(ev.x)/float(nxpix)-xoff)*(v.xmax-v.xmin)/(xpic-xoff)+v.xmin
            v.zmax=(float(ev.y)/float(nypix)-yoff)*(v.zmax-v.zmin)/(ypic-yoff)+v.zmin
            if (v.zmax gt zmax) then v.zmax=zmax
            if (v.xmin lt 0.0) then v.xmin=0.
            widget_control,ev.top,set_uvalue=v    
            print,"xmin=",v.xmin," zmax=",v.zmax
        endif
        if (ev.type eq 1) then begin
            v.xmax=(float(ev.x)/float(nxpix)-xoff)*(v.xmax-xmin0)/(xpic-xoff)+xmin0
            v.zmin=(float(ev.y)/float(nypix)-yoff)*(zmax0-v.zmin)/(ypic-yoff)+v.zmin
            if (v.zmin lt 0 ) then v.zmin=0
            if (v.xmax gt xmax) then v.xmax=xmax
            widget_control,ev.top,set_uvalue=v  
            print,"xmax=",v.xmax," zmin=",v.zmin  
            doplot,0
        endif
    end
    
    'computed' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.computed=ev.index
        widget_control,ev.top,set_uvalue=v    
        update = 1
     end

    'angle' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.angle=ev.value
        widget_control,ev.top,set_uvalue=v    
        update = 1
    end

    'tave' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.tave=ev.value+1
        widget_control,ev.top,set_uvalue=v    
        update = 1
    end

    'overplot' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.overplot=ev.index
        widget_control,ev.top,set_uvalue=v    
    end

    'rawplot' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.rawplot=ev.index
        widget_control,ev.top,set_uvalue=v    
        update = 1
    end

    'options' : begin
        widget_control,ev.top,get_uvalue=v                	
        print," Selected Option = ",ev.value       
        if ( ev.value ge 2 ) and ( ev.value le 9 ) then v.smoothing=ev.value-1
        if ( ev.value ge 11 ) and ( ev.value le 30) then v.contours=(ev.value*2)-10
        if ( ev.value eq 31 ) then v.data=0 
        if ( ev.value eq 32 ) then v.data=1
        if ( ev.value ge 35 ) then v.map=35-ev.value
        if ( ev.value eq 36 ) then xloadct,bottom=5
        widget_control,ev.top,set_uvalue=v    
        update = 1
    end

    'ptype' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.plottype=ev.index
        widget_control,ev.top,set_uvalue=v
        doplot,0
    end

    'time' : begin
        widget_control,ev.top,get_uvalue=v                	 	
        v.time=ev.value
        widget_control,ev.top,set_uvalue=v 
        update = 1
        doplot,0
    end

    'xplane' : begin
        widget_control,ev.top,get_uvalue=v                	 	
        v.xcut=ev.value
        widget_control,ev.top,set_uvalue=v 
        update = 1
        doplot,0
    end

    'yplane' : begin
        widget_control,ev.top,get_uvalue=v                	 	
        v.ycut=ev.value
        widget_control,ev.top,set_uvalue=v 
        update = 1
        doplot,0
    end

    'zplane' : begin
        widget_control,ev.top,get_uvalue=v                	 	
        v.zcut=ev.value
        widget_control,ev.top,set_uvalue=v 
        update = 1
        doplot,0
    end

    'shift' : begin
        widget_control,ev.top,get_uvalue=v                	 	
        v.shift=ev.value
        widget_control,ev.top,set_uvalue=v 
        update = 1
        doplot,0
    end

    'zmax' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.zmax=ev.value + zmax/2
        v.zmin=zmax/2  - ev.value
        widget_control,ev.top,set_uvalue=v    
    end

    'xslice' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.xslice=ev.value
        widget_control,ev.top,set_uvalue=v    
    end

    'zslice' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.zslice=ev.value+zmax/2
        widget_control,ev.top,set_uvalue=v    
    end

    'xmax' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.xmax=ev.value
        widget_control,ev.top,set_uvalue=v    
    end

    'xmin' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.xmin=ev.value
        widget_control,ev.top,set_uvalue=v    
    end

    'ymax' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.ymax=ev.value
        widget_control,ev.top,set_uvalue=v    
    end

    'ymin' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.ymin=ev.value
        widget_control,ev.top,set_uvalue=v    
    end


endcase

end

; NAME: COLORBAR
; PURPOSE:
;       The purpose of this routine is to add a color bar to the current
;       graphics window.
; AUTHOR:
;   FANNING SOFTWARE CONSULTING
;   David Fanning, Ph.D.
;   2642 Bradbury Court
;   Fort Collins, CO 80521 USA
;   Phone: 970-221-0438
;   E-mail: davidf@dfanning.com
;   Coyote's Guide to IDL Programming: http://www.dfanning.com/
; CATEGORY:
;       Graphics, Widgets.
; CALLING SEQUENCE:
;       COLORBAR
; INPUTS:
;       None.
; KEYWORD PARAMETERS:
;       BOTTOM:   The lowest color index of the colors to be loaded in
;                 the bar.
;       CHARSIZE: The character size of the color bar annotations. Default is 1.0.
;
;       COLOR:    The color index of the bar outline and characters. Default
;                 is !P.Color..
;       DIVISIONS: The number of divisions to divide the bar into. There will
;                 be (divisions + 1) annotations. The default is 6.
;       FONT:     Sets the font of the annotation. Hershey: -1, Hardware:0, 
;                 True-Type: 1.
;       FORMAT:   The format of the bar annotations. Default is '(I5)'.
;       MAXRANGE: The maximum data value for the bar annotation. Default is
;                 NCOLORS.
;       MINRANGE: The minimum data value for the bar annotation. Default is 0.
;       MINOR:    The number of minor tick divisions. Default is 2.
;       NCOLORS:  This is the number of colors in the color bar.
;
;       POSITION: A four-element array of normalized coordinates in the same
;                 form as the POSITION keyword on a plot. Default is
;                 [0.88, 0.15, 0.95, 0.95] for a vertical bar and
;                 [0.15, 0.88, 0.95, 0.95] for a horizontal bar.
;       RANGE:    A two-element vector of the form [min, max]. Provides an
;                 alternative way of setting the MINRANGE and MAXRANGE keywords.
;       RIGHT:    This puts the labels on the right-hand side of a vertical
;                 color bar. It applies only to vertical color bars.
;       TITLE:    This is title for the color bar. The default is to have
;                 no title.
;       TOP:      This puts the labels on top of the bar rather than under it.
;                 The keyword only applies if a horizontal color bar is rendered.
;       VERTICAL: Setting this keyword give a vertical color bar. The default
;                 is a horizontal color bar.
; COMMON BLOCKS:
;       None.
; SIDE EFFECTS:
;       Color bar is drawn in the current graphics window.
; RESTRICTIONS:
;       The number of colors available on the display device (not the
;       PostScript device) is used unless the NCOLORS keyword is used.
; EXAMPLE:
;       To display a horizontal color bar above a contour plot, type:
;       LOADCT, 5, NCOLORS=100
;       CONTOUR, DIST(31,41), POSITION=[0.15, 0.15, 0.95, 0.75], $
;          C_COLORS=INDGEN(25)*4, NLEVELS=25
;       COLORBAR, NCOLORS=100, POSITION=[0.15, 0.85, 0.95, 0.90]
; MODIFICATION HISTORY:
;       Written by: David Fanning, 10 JUNE 96.
;       10/27/96: Added the ability to send output to PostScript. DWF
;       11/4/96: Substantially rewritten to go to screen or PostScript
;           file without having to know much about the PostScript device
;           or even what the current graphics device is. DWF
;       1/27/97: Added the RIGHT and TOP keywords. Also modified the
;            way the TITLE keyword works. DWF
;       7/15/97: Fixed a problem some machines have with plots that have
;            no valid data range in them. DWF
;       12/5/98: Fixed a problem in how the colorbar image is created that
;            seemed to tickle a bug in some versions of IDL. DWF.
;       1/12/99: Fixed a problem caused by RSI fixing a bug in IDL 5.2. 
;                Sigh... DWF.
;       3/30/99: Modified a few of the defaults. DWF.
;       3/30/99: Used NORMAL rather than DEVICE coords for positioning bar. DWF.
;       3/30/99: Added the RANGE keyword. DWF.
;       3/30/99: Added FONT keyword. DWF
;       5/6/99: Many modifications to defaults. DWF.
;       5/6/99: Removed PSCOLOR keyword. DWF.
;       5/6/99: Improved error handling on position coordinates. DWF.
;       5/6/99. Added MINOR keyword. DWF.
;       5/6/99: Set Device, Decomposed=0 if necessary. DWF.
;       2/9/99: Fixed a problem caused by setting BOTTOM keyword, but not 
;               NCOLORS. DWF.
;       8/17/99. Fixed a problem with ambiguous MIN and MINOR keywords. DWF
;       8/25/99. I think I *finally* got the BOTTOM/NCOLORS thing sorted out. 
;                :-( DWF.
;       10/10/99. Modified the program so that current plot and map 
;                 coordinates are
;            saved and restored after the colorbar is drawn. DWF.
;       3/18/00. Moved a block of code to prevent a problem with color 
;                decomposition. DWF.
;       4/28/00. Made !P.Font default value for FONT keyword. DWF.
;###########################################################################
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright Å© 2000 Fanning Software Consulting.
;
; This software is provided "as-is", without any express or
; implied warranty. In no event will the authors be held liable
; for any damages arising from the use of this software.
; Permission is granted to anyone to use this software for any
; purpose, including commercial applications, and to alter it and
; redistribute it freely, subject to the following restrictions:
; 1. The origin of this software must not be misrepresented; you must
;    not claim you wrote the original software. If you use this software
;    in a product, an acknowledgment in the product documentation
;    would be appreciated, but is not required.
; 2. Altered source versions must be plainly marked as such, and must
;    not be misrepresented as being the original software.
; 3. This notice may not be removed or altered from any source distribution.
; For more information on Open Source Software, visit the Open Source
; web site: http://www.opensource.org.
;###########################################################################
PRO COLORBAR, BOTTOM=bottom, CHARSIZE=charsize, COLOR=color, $
   DIVISIONS=divisions, $
   FORMAT=format, POSITION=position, MAXRANGE=maxrange, MINRANGE=minrange, $
   NCOLORS=ncolors, $
   TITLE=title, VERTICAL=vertical, TOP=top, RIGHT=right, MINOR=minor, $
   RANGE=range, FONT=font, TICKLEN=ticklen, _EXTRA=extra

   ; Return to caller on error.
On_Error, 2

   ; Save the current plot state.

bang_p = !P
bang_x = !X
bang_Y = !Y
bang_Z = !Z
bang_Map = !Map

   ; Is the PostScript device selected?

postScriptDevice = (!D.NAME EQ 'PS' OR !D.NAME EQ 'PRINTER')

   ; Which release of IDL is this?

thisRelease = Float(!Version.Release)

    ; Check and define keywords.

IF N_ELEMENTS(ncolors) EQ 0 THEN BEGIN

   ; Most display devices to not use the 256 colors available to
   ; the PostScript device. This presents a problem when writing
   ; general-purpose programs that can be output to the display or
   ; to the PostScript device. This problem is especially bothersome
   ; if you don't specify the number of colors you are using in the
   ; program. One way to work around this problem is to make the
   ; default number of colors the same for the display device and for
   ; the PostScript device. Then, the colors you see in PostScript are
   ; identical to the colors you see on your display. Here is one way to
   ; do it.

   IF postScriptDevice THEN BEGIN
      oldDevice = !D.NAME

         ; What kind of computer are we using? SET_PLOT to appropriate
         ; display device.

      thisOS = !VERSION.OS_FAMILY
      thisOS = STRMID(thisOS, 0, 3)
      thisOS = STRUPCASE(thisOS)
      CASE thisOS of
         'MAC': SET_PLOT, thisOS
         'WIN': SET_PLOT, thisOS
         ELSE: SET_PLOT, 'X'
      ENDCASE

         ; Here is how many colors we should use.

      ncolors = !D.TABLE_SIZE
      SET_PLOT, oldDevice
    ENDIF ELSE ncolors = !D.TABLE_SIZE
ENDIF
IF N_ELEMENTS(bottom) EQ 0 THEN bottom = 0B
IF N_ELEMENTS(charsize) EQ 0 THEN charsize = 1.0
IF N_ELEMENTS(format) EQ 0 THEN format = '(I5)'
IF N_ELEMENTS(color) EQ 0 THEN color = !P.Color
IF N_ELEMENTS(minrange) EQ 0 THEN minrange = 0
IF N_ELEMENTS(maxrange) EQ 0 THEN maxrange = ncolors
IF N_ELEMENTS(ticklen) EQ 0 THEN ticklen = 0.2
IF N_ELEMENTS(minor) EQ 0 THEN minor = 2
IF N_ELEMENTS(range) NE 0 THEN BEGIN
   minrange = range[0]
   maxrange = range[1]
ENDIF
IF N_ELEMENTS(divisions) EQ 0 THEN divisions = 6
IF N_ELEMENTS(font) EQ 0 THEN font = !P.Font
IF N_ELEMENTS(title) EQ 0 THEN title = ''

IF KEYWORD_SET(vertical) THEN BEGIN
   bar = REPLICATE(1B,20) # BINDGEN(ncolors)
   IF N_ELEMENTS(position) EQ 0 THEN BEGIN
      position = [0.88, 0.1, 0.95, 0.9]
   ENDIF ELSE BEGIN
      IF position[2]-position[0] GT position[3]-position[1] THEN BEGIN
         position = [position[1], position[0], position[3], position[2]]
      ENDIF
      IF position[0] GE position[2] THEN Message, "Position coordinates can't be reconciled."
      IF position[1] GE position[3] THEN Message, "Position coordinates can't be reconciled."
   ENDELSE
ENDIF ELSE BEGIN
   bar = BINDGEN(ncolors) # REPLICATE(1B, 20)
   IF N_ELEMENTS(position) EQ 0 THEN BEGIN
      position = [0.1, 0.88, 0.9, 0.95]
   ENDIF ELSE BEGIN
      IF position[3]-position[1] GT position[2]-position[0] THEN BEGIN
         position = [position[1], position[0], position[3], position[2]]
      ENDIF
      IF position[0] GE position[2] THEN Message, "Position coordinates can't be reconciled."
      IF position[1] GE position[3] THEN Message, "Position coordinates can't be reconciled."
   ENDELSE
ENDELSE

   ; Scale the color bar.

 bar = BYTSCL(bar, TOP=(ncolors-1 < (255-bottom))) + bottom

   ; Get starting locations in NORMAL coordinates.

xstart = position(0)
ystart = position(1)

   ; Get the size of the bar in NORMAL coordinates.

xsize = (position(2) - position(0))
ysize = (position(3) - position(1))

   ; Display the color bar in the window. Sizing is
   ; different for PostScript and regular display.

IF postScriptDevice THEN BEGIN

   TV, bar, xstart, ystart, XSIZE=xsize, YSIZE=ysize, /Normal

ENDIF ELSE BEGIN

   bar = CONGRID(bar, CEIL(xsize*!D.X_VSize), CEIL(ysize*!D.Y_VSize), /INTERP)

        ; Decomposed color off if device supports it.

   CASE  StrUpCase(!D.NAME) OF
        'X': BEGIN
            IF thisRelease GE 5.2 THEN Device, Get_Decomposed=thisDecomposed
            Device, Decomposed=0
            ENDCASE
        'WIN': BEGIN
            IF thisRelease GE 5.2 THEN Device, Get_Decomposed=thisDecomposed
            Device, Decomposed=0
            ENDCASE
        'MAC': BEGIN
            IF thisRelease GE 5.2 THEN Device, Get_Decomposed=thisDecomposed
            Device, Decomposed=0
            ENDCASE
        ELSE:
   ENDCASE

   TV, bar, xstart, ystart, /Normal

      ; Restore Decomposed state if necessary.

   CASE StrUpCase(!D.NAME) OF
      'X': BEGIN
         IF thisRelease GE 5.2 THEN Device, Decomposed=thisDecomposed
         ENDCASE
      'WIN': BEGIN
         IF thisRelease GE 5.2 THEN Device, Decomposed=thisDecomposed
         ENDCASE
      'MAC': BEGIN
         IF thisRelease GE 5.2 THEN Device, Decomposed=thisDecomposed
         ENDCASE
      ELSE:
   ENDCASE

ENDELSE

   ; Annotate the color bar.

IF KEYWORD_SET(vertical) THEN BEGIN

   IF KEYWORD_SET(right) THEN BEGIN

      PLOT, [minrange,maxrange], [minrange,maxrange], /NODATA, XTICKS=1, $
         YTICKS=divisions, XSTYLE=1, YSTYLE=9, $
         POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
         YTICKFORMAT='(A1)', XTICKFORMAT='(A1)', YTICKLEN=ticklen , $
         YRANGE=[minrange, maxrange], FONT=font, _EXTRA=extra, YMINOR=minor

      AXIS, YAXIS=1, YRANGE=[minrange, maxrange], YTICKFORMAT=format, YTICKS=divisions, $
         YTICKLEN=ticklen, YSTYLE=1, COLOR=color, CHARSIZE=charsize, $
         FONT=font, YTITLE=title, _EXTRA=extra, YMINOR=minor

   ENDIF ELSE BEGIN

      PLOT, [minrange,maxrange], [minrange,maxrange], /NODATA, XTICKS=1, $
         YTICKS=divisions, XSTYLE=1, YSTYLE=9, YMINOR=minor, $
         POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
         YTICKFORMAT=format, XTICKFORMAT='(A1)', YTICKLEN=ticklen , $
         YRANGE=[minrange, maxrange], FONT=font, YTITLE=title, _EXTRA=extra

      AXIS, YAXIS=1, YRANGE=[minrange, maxrange], YTICKFORMAT='(A1)', YTICKS=divisions, $
         YTICKLEN=ticklen, YSTYLE=1, COLOR=color, CHARSIZE=charsize, $
         FONT=font, _EXTRA=extra, YMINOR=minor

   ENDELSE

ENDIF ELSE BEGIN

   IF KEYWORD_SET(top) THEN BEGIN

      PLOT, [minrange,maxrange], [minrange,maxrange], /NODATA, XTICKS=divisions, $
         YTICKS=1, XSTYLE=9, YSTYLE=1, $
         POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
         YTICKFORMAT='(A1)', XTICKFORMAT='(A1)', XTICKLEN=ticklen, $
         XRANGE=[minrange, maxrange], FONT=font, _EXTRA=extra, XMINOR=minor

      AXIS, XTICKS=divisions, XSTYLE=1, COLOR=color, CHARSIZE=charsize, $
         XTICKFORMAT=format, XTICKLEN=ticklen, XRANGE=[minrange, maxrange], XAXIS=1, $
         FONT=font, XTITLE=title, _EXTRA=extra, XCHARSIZE=charsize, XMINOR=minor

   ENDIF ELSE BEGIN

      PLOT, [minrange,maxrange], [minrange,maxrange], /NODATA, XTICKS=divisions, $
         YTICKS=1, XSTYLE=1, YSTYLE=1, TITLE=title, $
         POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
         YTICKFORMAT='(A1)', XTICKFORMAT=format, XTICKLEN=ticklen, $
         XRANGE=[minrange, maxrange], FONT=font, XMinor=minor, _EXTRA=extra

    ENDELSE

ENDELSE

   ; Restore the previous plot and map system variables.

!P = bang_p
!X = bang_x
!Y = bang_y
!Z = bang_z
!Map = bang_map

END

pro create_movie
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct
common colortable,rgb,usecolor,red,blue,green,range1,range2,r1,r2,tmax,tmin
common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common imagesize,nxpix,nypix,xoff,yoff,xpic,ypic

; Temporary storage

myimage = bytarr(nxpix,nypix)

; Open MPEG sequence

mpegID = MPEG_OPEN([nxpix,nypix],QUALITY=100)

; Set intital frame and number of padded frames

iframe = 0
ipad = 1
iskip = 1

;  Loop over all time slices

itime=0
for itime = 0,tslices,iskip do begin
;for itime = 0,400,iskip do begin
v.time = itime
doplot,2

; Add image in current window to MPEG movie

for i=1,ipad do begin
    iframe = iframe + 1
    MPEG_PUT, mpegID, WINDOW=!D.Window, FRAME=iframe, /ORDER
endfor

; end of loop

end

; SAVE the MPEG sequence

if (v.rawplot ge 1) then fname = rawdata(v.rawplot) else fname=computedquantity(v.computed)
MPEG_SAVE, mpegID,FILENAME=strcompress(fname+'.mpg',/remove_all)

; Close the MPEG sequence

MPEG_CLOSE, mpegID

;  end of subroutine
end

pro diagohm
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq, anames
common pdata,fulldata,simulation_time,mytitle
common picdata, field, struct
common colortable,rgb,usecolor,red,blue,green,range1,range2,r1,r2,tmax,tmin
common controlplot,v,ipe,ipi,ibx,iby,ibz,ine,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common imagesize,nxpix,nypix,xoff,yoff,xpic,ypic
common refresh,update

; First time starting

update = 0

; First determine the screen size

dimensions = GET_SCREEN_SIZE(RESOLUTION=resolution)
print,"Dimensions of screen=",dimensions

; Pick size of viewer based on fraction of screen

nxpix = dimensions(0)*0.44
nypix = dimensions(1)*0.44

; Hard code in your choice (better for movies)

;nxpix=1000
;nypix=600
;nxpix=900
;nypix=540

; Declare structure for controlling the plots

v={time:0,xmin:0.0,xmax:0.0,ymin:0.0,ymax:0.0,xslice:0.0,zmin:0.0,zmax:0.0,zslice:0.0,smoothing:1,contours:24,plottype:0,overplot:0,rawplot:0,shift:0,data:0,map:0,xcut:1,ycut:0,zcut:1,computed:0,angle:0.0,tave:0}

; Declare strings for menu

ptype = strarr(4)  
porient = strarr(3)
anames = strarr(3,2)

; Read in my color table

nc=256
rgb=fltarr(3,nc)
usecolor=fltarr(3,nc)
red=fltarr(nc)
blue=fltarr(nc)
green=fltarr(nc)
openr, unit, '~/Color_Tables/c5.tbl', /get_lun
readf, unit, rgb
close, unit
free_lun, unit
red=rgb(0,*)
green=rgb(1,*)
blue=rgb(2,*)
tvlct,red,green,blue

; Declare long integers 

nx=0L
ny=0L
nz=0L
numq=0L

; Determine if there are multiple data directories and allow user to select

datafiles = file_search('dat*',count=numq) 
if (numq gt 1) then directory = dialog_pickfile(/directory,filter='data*',TITLE='Choose directory with data') else  directory = 'data'
    
; open binary file for problem description 

if ( file_test(directory+'/info') eq 1 ) then begin
    print," *** Found Data Information File *** "
endif else begin
    print," *** Error - File info is missing ***"
endelse

;  First Try little endian then switch to big endian if this fails

little = 0
on_ioerror, switch_endian
openr,unit,directory+'/info',/f77_unformatted, /get_lun,/swap_if_big_endian
readu,unit,nx,nz,ny
little=1
switch_endian: if (not little) then begin
    print, " ** Little endian failed --> Switch to big endian"
    close,unit
    free_lun,unit
    openr,unit,directory+'/info',/f77_unformatted, /get_lun,/swap_if_little_endian
    readu,unit,nx,nz,ny
endif

; Read the problem desciption

on_ioerror, halt_error1
readu,unit,xmax,zmax,ymax
close,unit
free_lun, unit

;xmax = xmax/10.0
;ymax = ymax/10.0
;zmax = zmax/10.0

; Find the names of data files in the data directory

datafiles = file_search(directory+'/*.gda',count=numq) 

; Now open each file and save the basename to identify later

print," Number of files=",numq
rawdata = strarr(numq+1)  
instring='     '
rawdata(0)='None'
for i=1,numq do begin
    if (not little) then openr,i,datafiles(i-1) else openr,i,datafiles(i-1),/swap_if_big_endian
    rawdata(i) = file_basename(datafiles(i-1),'.gda')
    print,"i=",i," --> ",rawdata(i)
    if ( rawdata(i) eq 'pe-xx' or rawdata(i) eq 'Pe-xx' ) then ipe =i
    if ( rawdata(i) eq 'pi-xx' or rawdata(i) eq 'Pi-xx' ) then ipi =i
    if ( rawdata(i) eq 'bx' or rawdata(i) eq 'Bx' ) then ibx =i
    if ( rawdata(i) eq 'by' or rawdata(i) eq 'By' ) then iby =i
    if ( rawdata(i) eq 'bz' or rawdata(i) eq 'Bz' ) then ibz =i
    if ( rawdata(i) eq 'ex' or rawdata(i) eq 'Ex' ) then iex =i
    if ( rawdata(i) eq 'ey' or rawdata(i) eq 'Ey' ) then iey =i
    if ( rawdata(i) eq 'ez' or rawdata(i) eq 'Ez' ) then iez =i
    if ( rawdata(i) eq 'jx' or rawdata(i) eq 'Jx' ) then ijx =i
    if ( rawdata(i) eq 'jy' or rawdata(i) eq 'Jy' ) then ijy =i
    if ( rawdata(i) eq 'jz' or rawdata(i) eq 'Jz' ) then ijz =i
    if ( rawdata(i) eq 'ne' ) then ine =i
    if ( rawdata(i) eq 'uix' or rawdata(i) eq 'Uix' ) then iuix =i
    if ( rawdata(i) eq 'uiy' or rawdata(i) eq 'Uiy' ) then iuiy =i
    if ( rawdata(i) eq 'uiz' or rawdata(i) eq 'Uiz' ) then iuiz =i
    if ( rawdata(i) eq 'uex' or rawdata(i) eq 'Uex' ) then iuex =i
    if ( rawdata(i) eq 'uey' or rawdata(i) eq 'Uey' ) then iuey =i
    if ( rawdata(i) eq 'uez' or rawdata(i) eq 'Uez' ) then iuez =i
 endfor

; Define different types of reduced diagnostic

computedquantity = strarr(23)  
computedquantity(0) = "Electron Anisotropy"
computedquantity(1) = "Ion Anisotropy"
computedquantity(2) = "Electron Agyrotopy"
computedquantity(3) = "Ion Agyrotopy"
computedquantity(4) = "Firehose Condition"
computedquantity(5) = "Electron Firehose"
computedquantity(6) = "KH Condition"
computedquantity(7) = "Magnetic Shear Angle"
computedquantity(8) = "Total Current"
computedquantity(9) = "Beta-e"
computedquantity(10) = "Beta-i"
computedquantity(11) = "Beta-Total"
computedquantity(12) = "E*J"
computedquantity(13) = "Te"
computedquantity(14) = "Ti"
computedquantity(15) = "Total-T"
computedquantity(16) = "EUexB_y"
computedquantity(17) = "dtnUey"
computedquantity(18) = "divneUeUey"
computedquantity(19) = "divPey"
computedquantity(20) = "Gyroviscous cancellation"
computedquantity(21) = "Gyroterms"
computedquantity(22) = "Hesse field"
; Close the input file

close,unit
free_lun, unit

; Define zcenter 

zcenter=zmax/2.0

; Echo information

print,'nx=',nx,'  ny=',ny,'  nz=',nz
print,'xmax=',xmax,'  ymax=',ymax,'  zmax=',zmax

; Plotting array

fulldata = fltarr(nx,nz)

; Define structure of data files

; Bill's way of saving gda data
;struct = {data:fltarr(nx,nz,ny),time:0.0,it:500000}
;record_length = 4L*(nx*ny*nz+2L)

; Homa's way of saving gda data

struct = {data:fltarr(nx,nz,ny)}
record_length = 4L*(nx*ny*nz)

; Determine number of time slices

information=file_info(datafiles(0))
tslices=information.size/record_length-1
print,"File Size=",information.size
print,"Record Length=",record_length
print,"Time Slices",tslices
if (tslices lt 1) then tslices=1

; Plot type

ptype(0)='Contour'
ptype(1)='Contour+X-Average'
ptype(2)='Contour+Z-Slice'
ptype(3)='Contour+X-Slice'

; plane orientation

porient(0)='XZ'
porient(1)='XY'
porient(2)='ZY'

; axis names
anames(0,0) = 'Z'
anames(0,1) = 'X'

anames(1,0) = 'Y'
anames(1,1) = 'X'

anames(2,0) = 'Z'
anames(2,1) = 'Y'


; Options menu

desc =[ '1\Options' , $
        '1\Smoothing' , $
        '0\1' , $
        '0\2' , $
        '0\3' , $
        '0\4' , $
        '0\5' , $
        '0\6' , $
        '0\7' , $
        '2\8' , $
        '1\Contours' , $
        '0\12' , $
        '0\14' , $
        '0\16' , $
        '0\18' , $
        '0\20' , $
        '0\22' , $
        '0\24' , $
        '0\26' , $
        '0\28' , $
        '0\30' , $
        '0\32' , $
        '0\34' , $
        '0\36' , $
        '0\38' , $
        '0\40' , $
        '0\42' , $
        '0\44' , $
        '0\46' , $
        '2\48' , $
        '1\Transform Data' , $
        '0\Full Data' , $
        '0\Perturbed', $
        '2\',  $
        '1\Color Map', $
        '0\Default' , $
        '2\Load Table' ]

; Setup widgets

base = widget_base(/column)
row1 = widget_base(base,/row,scr_ysize=nypix/7)
row2 = widget_base(base,/row)
row3 = widget_base(base,/row)
button1 = widget_button(row1, value = '  Done  ', uvalue = 'done',/sensitive)
button2 = widget_button(row1, value = '  Plot  ', uvalue = 'plot',/sensitive)
button3 = widget_button(row1, value = 'Render PS', uvalue = 'render')
button4 = widget_button(row1, value = 'Animate', uvalue = 'animate')
list5=widget_list(row1,value=ptype,uvalue='ptype')
list1=widget_droplist(row1,title=' Raw Plot',value=rawdata,uvalue='rawplot')
list2=widget_droplist(row1,title=' Diagnostic Plot',value=computedquantity,uvalue='computed')
;list6=widget_droplist(row1,title=' Orientation',value=porient,uvalue='orientation')
list4=widget_droplist(row1,title=' Overplot',value=rawdata,uvalue='overplot')
opt = cw_pdmenu(row1,desc,/return_index,uvalue='options')
slider10=cw_fslider(row1,title='Reference Angle',value=0.0,format='(f4.0)',uvalue='angle',minimum=0,maximum=180.0,/drag,scroll=1.0)
slider11=cw_fslider(row1,title='Time Average',value=0,format='(f2.0)',uvalue='tave',minimum=0,maximum=9,/drag,scroll=1)
slider2=cw_fslider(row2,title='X-min',value=0.0,format='(f7.1)',uvalue='xmin',minimum=0,maximum=xmax,/drag)
slider3=cw_fslider(row2,title='X-max',value=xmax,format='(f7.1)',uvalue='xmax',minimum=0,maximum=xmax,/drag)
slider2b=cw_fslider(row2,title='Y-min',value=0.0,format='(f7.1)',uvalue='ymin',minimum=0,maximum=ymax,/drag)
slider3b=cw_fslider(row2,title='Y-max',value=ymax,format='(f7.1)',uvalue='ymax',minimum=0,maximum=ymax,/drag)
slider6=cw_fslider(row2,title='Z-max',value=zcenter,format='(f7.1)',uvalue='zmax',minimum=0.0,maximum=zcenter,/drag)
slider4=cw_fslider(row2,title='X-Slice',value=xmax/2,format='(f7.1)',uvalue='xslice',minimum=0,maximum=xmax,/drag,scroll=1.0)
slider5=cw_fslider(row2,title='Z-Slice',value=0,format='(f7.1)',uvalue='zslice',minimum=-zmax/2,maximum=zmax/2,/drag,scroll=1.0)
slider7=widget_slider(row2,title=' X-Shift',scrol=1,value=0,uvalue='shift',minimum=-nx/2,maximum=nx/2)
slider1=widget_slider(row2,title='Time Slice',scrol=1,value=0,uvalue='time',minimum=0,maximum=tslices,scr_xsize=nxpix/3.6)

;slider8  =widget_slider(row3,title='x-plane',scrol=1,value=0,uvalue='xplane',minimum=0,maximum=nx-1,scr_xsize=nxpix/3.6)

; Set this slider to work with 2D or 3D
nymax = ny-1
if (ny eq 1) then nymax = 1 
slider9  =widget_slider(row3,title='y-plane',scrol=1,value=0,uvalue='yplane',minimum=0,maximum=nymax,scr_xsize=nxpix/3.6)

;slider10 =widget_slider(row3,title='z-plane',scrol=1,value=0,uvalue='zplane',minimum=0,maximum=nz-1,scr_xsize=nxpix/3.6)

draw = widget_draw(base,retain=2, xsize = nxpix, ysize = nypix,/button_events,uvalue='mouse')
widget_control, base, /realize
widget_control,list5,set_list_select=0
widget_control,list1,set_droplist=0
widget_control,list2,set_droplist=0
widget_control, base, set_uvalue={time:0,xmin:0.0,xmax:xmax,ymin:0.0,ymax:ymax,xslice:xmax/2,zmin:0.0,zmax:zmax,zslice:zcenter,smoothing:2,contours:24,plottype:0,overplot:0,rawplot:0,shift:0,data:0,map:0,orientation:0,xcut:1,ycut:0,zcut:1,computed:0,angle:0.0,tave:0}
widget_control, draw, get_value = index
xmanager, 'handle', base
halt_error1: print, ' *** Halting Program ***'
end

