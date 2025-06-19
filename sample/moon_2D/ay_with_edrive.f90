!---------------------------------------------------------------------------------------
!  This program computes the flux surfaces - Ay, by reading the the inplane magnetic field
!  components from data/bx.gda and data/bz.gda, taking the curl, and solving the needed
!  Poisson equation
!
!  This version is for a single sheet - periodic in x, with conducting BC in z
!---------------------------------------------------------------------------------------

program translate
  implicit none
  integer input_record,output_record,input_error,output_error,record_length,it1,it2,nx,ny,nz,i,j,it,interval
  real(kind=4)diff,xmax,ymax,zmax,time
  real(kind=4), allocatable, dimension(:,:) :: bx,bz,Ay
  real(kind=8) dx,dz
  real(kind=4) ay0,intervalwci,mime,wpewce,intervalwpe,b0,MA,vflow,theta,E0
  character(22) :: fname
! Range of time slices to convert.  Make big and let error terminate program

  it1=1
  it2=1300
  output_record = 0

! Open data information file to find out the size of the arrays

  open(unit=10,file='data/info',status='old',form='unformatted')

!  Read data from info  file

  read(10)nx,nz,ny
  read(10)xmax,zmax,ymax 

! Cell size

  dx = xmax/real(nx)
  dz = zmax/real(nz)

! Echo this information

  print *,"---------------------------------------------------"
  print *,"xmax=",xmax,"    zmax=",zmax
  print *,"nx=",nx,"   nz=",nz
  print *,"dx=",dx,"   dz=",dz
  print *,"---------------------------------------------------"

! Allocate storage space for fields and moments

  allocate(bx(nx,nz))
  allocate(bz(nx,nz))
  allocate(Ay(nx,nz))

! Set record length for gda files

! Use this for Bill's gda format with extra time and it records
!  inquire(iolength=record_length)bx,time,it

! Use this if you are only saving the matrix in the gda file
  inquire(iolength=record_length)bx

! *** WARNING - also make the read statments below consistent ***

  print *," Setting record length=",record_length

! Open all of the direct access binary files (Bill's preferred data format)

!  open(unit=20,file='data/bx.gda',access='direct',recl=record_length,status='unknown',form='unformatted',action='read')     
!  open(unit=30,file='data/by.gda',access='direct',recl=record_length,status='unknown',form='unformatted',action='read')    

! Loop over time slices

  interval = 819
  ay0 = 0.0
  intervalwci = 0.5
  mime = 100
  wpewce = 2
  intervalwpe = intervalwci*mime*wpewce
  b0 = 1/wpewce
  MA = 6.5
  vflow = MA*b0/sqrt(mime)
  theta = -30.0*3.14159/180.0
  E0 = -vflow*b0*sin(theta)

  print*,"E0=",E0
  do input_record = it1,it2

     write(fname,'(a,i0,a)')'data/bx_',(input_record-1)*interval,'.gda'
     open(unit=20,file=fname,access='direct',recl=record_length,status='unknown',form='unformatted',action='read')
     write(fname,'(a,i0,a)')'data/by_',(input_record-1)*interval,'.gda'
     open(unit=30,file=fname,access='direct',recl=record_length,status='unknown',form='unformatted',action='read')

     print *,"Reading record=",input_record

! Read magnetic field data

     read(20,rec=1)bx
     read(30,rec=1)bz

! Comute instant Ay(x,z) by integration from lower left corner Ay(1,1)=0.0
     ! Apply driving electric field
  Ay(1,1)=ay0
  ay0 = ay0 - E0*intervalwpe
!  Ay(1,1)=ay0
 !calculate the flux function, Ay : B_in= y x grad(Ay)

  do i=2,nx 
  Ay(i,1)=-dx*bz(i-1,1)+Ay(i-1,1)
  enddo

  do i=1,nx 
  do j=2,nz 
  Ay(i,j)=dz*bx(i,j-1)+Ay(i,j-1)
  enddo
  enddo

! Save flux surfaces

       print *,"Saving Ay"
       output_record = output_record +1
       write(fname,'(a,i0,a)')'data/ay_',(input_record-1)*interval,'.gda'
       open(unit=40,file=fname,access='direct',recl=record_length,status='unknown',form='unformatted',action='write')     
       write(40,rec=1)Ay
       close(40)

  enddo

! Close files

  close(20)
  close(30)

end program translate





