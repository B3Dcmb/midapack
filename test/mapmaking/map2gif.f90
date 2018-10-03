MODULE mapOutput

	USE healpix_modules
	USE healpix_types
	USE gifmod

	CONTAINS

	SUBROUTINE map2gif(nside, dpmap, filenamelen, output_file)

		IMPLICIT NONE

		INTEGER(I4B) :: nside
		REAL(SP), DIMENSION(0:12*nside**2-1) :: map
		REAL(DP), DIMENSION(0:12*nside**2-1) :: dpmap

		INTEGER(I4B) :: filenamelen

		CHARACTER(LEN=filenamelen) :: output_file


		map = dpmap !convert double precision to single precision

		print *, ' nside          = ', nside
		print *, ' output_file    = ', output_file

		call inside_map2gif(nside, map, filenamelen, output_file)

		return

	END SUBROUTINE map2gif

	!! ---------------------- MAP2GIF code from HEALpix -------------------------------------

	SUBROUTINE inside_map2gif(nside, map, fnl, output_file)

		IMPLICIT NONE

		INTEGER(I4B) :: nsmax
		INTEGER(I4B) :: npixtot, nside
		INTEGER(I4B) :: ordering
		INTEGER(I4B) :: i

		INTEGER(I4B) :: fnl

		REAL(SP),     DIMENSION(:,:),   ALLOCATABLE :: map_IQU
		REAL(SP),     DIMENSION(:)		    :: map
		REAL(SP),     DIMENSION(:,:),   ALLOCATABLE :: image

		INTEGER(I4B), ALLOCATABLE, DIMENSION(:,:) :: imgint

		INTEGER(I4B), POINTER, DIMENSION(:,:) :: imgbar
		INTEGER(I4B), POINTER, DIMENSION(:,:) :: imgttl

		LOGICAL(LGT), DIMENSION(:,:),   ALLOCATABLE :: mask
		LOGICAL(LGT), DIMENSION(:),     ALLOCATABLE :: observed

		REAL(SP)     :: mindata, maxdata, Tmin, Tmax

		INTEGER(I4B) :: status

		CHARACTER(LEN=7), PARAMETER :: code = 'MAP2GIF'
		character(len=*), parameter :: VERSION = HEALPIX_VERSION

		CHARACTER(LEN=fnl) :: output_file
		CHARACTER(LEN=filenamelen) :: title = ''

		CHARACTER(LEN=3)   :: projection = 'MOL'
		REAL(SP)           :: usermax = 1.1e30
		REAL(SP)           :: usermin = -1.1e30
		REAL(SP)           :: lon0 = 0.0
		REAL(SP)           :: lat0 = 0.0
		REAL(SP)           :: gridresn = 2.0
		INTEGER(I4B)       :: color_table = 5
		INTEGER(I4B)       :: signal = 1         ! I/intensity measurement
		INTEGER(I4B)       :: xsize  = 1280
		LOGICAL(LGT)       :: logflg = .false.
		LOGICAL(LGT)       :: bar = .true.
		LOGICAL(LGT)       :: ashflag = .false.
		REAL(SP)           :: offset = 0.0
		REAL(SP)           :: factor = 1.0

		integer(i4b)       :: n_ext, i_ext, nmaps_sum, icolumn



		!-------------------------------------------------------------------

		npixtot = 12*nside**2;
		bar = .true.
		nsmax = NINT(SQRT( npixtot / 12.) )
		call assert(npixtot == 12*nsmax*nsmax,'wrong number of pixels')



		ordering = 1; !ring



		allocate(observed(0:npixtot-1),stat = status)
		call assert_alloc(status,'map2gif','observed')

		!--- Test for unobserved pixels ---
		observed = map > -1.6375e30
		!!where (observed) map = (map+offset)*factor ! << requires large stack size for large maps
		do i = 0, npixtot-1
		   if (observed(i)) map(i) = (map(i)+offset)*factor
		enddo
		if (ashflag) then
		   do i=0,npixtot-1
		if (observed(i)) then
		  if (map(i)>=0) then
			 map(i) = log(map(i)+sqrt(map(i)**2+1))
		  else
			 map(i) = -log(-map(i)+sqrt(map(i)**2+1))
		  endif
		endif
		   end do
		endif
		IF( logflg )THEN
		   observed = map > 0.
	!       where(observed)
	!         map=log10(map)
	!       elsewhere
	!         map=-1.6375e30
	!       endwhere
		   do i=0, npixtot - 1
		 if (observed(i)) then
			 map(i) = log10(map(i))
		 else
			 map(i) = -1.6375e30
		 endif
		   enddo
		END IF
		mindata = minval(map, mask = observed)
		maxdata = maxval(map, mask = observed)

		if (mindata==maxdata) then
		  print *,"Warning: map contains completely uniform values!"
		  maxdata=mindata+1e-5
		endif

		!--- sets MIN and MAX using linear scaling ---
		if (ABS((maxdata+mindata)/(maxdata-mindata)) < 5.e-2) then
		   !--- center of color scale is 0
		   Tmax = MAXVAL(ABS((/mindata,maxdata/)))
		   Tmin = -Tmax
		else
		   Tmax = maxdata
		   Tmin = mindata
		end if

		!--- allow for user min/max settings ---
		if (usermax <  1.1e30) Tmax = usermax
		if (usermin > -1.1e30) Tmin = usermin

		!--- adjust 'offscale' temperatures ---
		!! WHERE(observed) map = max(min(map,Tmax),Tmin) ! << requires large stack size for large maps
		do i = 0, npixtot - 1
		   if (observed(i)) map(i) = max(min(map(i), Tmax), Tmin)
		enddo
		deallocate(observed)

		!--- create the projection ---
		select case (projection)
			case ('MOL', 'mol')
			   !--- allocate space for image and mask ---
			   ALLOCATE(image(0:xsize-1,0:xsize/2-1),stat = status)
			   call assert_alloc(status,'map2gif','image')
			   ALLOCATE(mask(0:xsize-1,0:xsize/2-1),stat = status)
			   call assert_alloc(status,'map2gif','mask')

			   !--- create the mollweide projection ---
			   call proj_mollw(&
				map,      &
				nsmax,    &
				ordering, &
				xsize,    &
				lon0,     &
				image,    &
				mask      &
				)

			! removed !!! -- case ('GNO', 'gno')

			case default
			   print '(" unknown projection type: ",a3)', projection
			   call fatal_error
		end select

		!--- load the color table ---
		call setcol( color_table )

		!--- scale image ---
		allocate(imgint(size(image,1),size(image,2)))
		call imgscl(&
		image,  &
		imgint, &
		mask    &
		)
		deallocate(image)

		!--- add color bar ---
		if (bar) then
		   call addbar(&
			imgint, &
			imgbar  &
			)
		   deallocate(imgint)
		   allocate(imgint(size(imgbar,1),size(imgbar,2)))
		   imgint = imgbar
		   deallocate(imgbar)
		end if

		!--- add title ---
		if (title /= '') then
		   call addttl(&
			imgint, &
			imgttl, &
			title   &
			)
		   deallocate(imgint)
		   allocate(imgint(size(imgttl,1),size(imgttl,2)))
		   imgint = imgttl
		   deallocate(imgttl)
		end if

		!--- create the GIF image ---
		call gifmap(&
		imgint,           &
		trim(output_file) &
		)

		!--- deallocate memory for arrays ---
		DEALLOCATE(mask,stat = status)

		!--- Exit, stage left ...



	END SUBROUTINE inside_map2gif







	SUBROUTINE PROJ_MOLLW(&
		  map,      &
		  nside,    &
		  ordering, &
		  xsize,    &
		  lon0,     &
		  image,    &
		  mask      &
		  )
		! --------------------------------------------------------------
		! This subroutine takes an input HEALPIX sky map array and
		! generates the (u,v) position on a mollweide projection.
		!
		! INPUTS:
		!       map     = array containing pixel temperatures
		!       nside   = defines number of pixels on the sky
		!       ordering= whether input map is in RING or NESTED format
		!       xsize   = x-dimension of output image array
		!       lon0    = longitude (in degrees) of the center of the plot
		!       image   = output mollweide projection image
		!       mask    = logical array defining usable pixels in output image array
		!
		! PROCEDURES USED:
		!            ang_pix, ang_nest_pix
		!
		!  RELATED LITERATURE:
		!  see the web site http://www.tac.dk/~healpix
		!
		! MODIFICATION HISTORY:
		!    October 1997, Eric Hivon, TAC (original IDL code)
		!       December 1998, A.J. Banday MPA (conversion to F90)
		!
		! --------------------------------------------------------------

		USE healpix_types
		USE pix_tools

		IMPLICIT NONE

		INTEGER(I4B) :: I,J
		INTEGER(I4B) :: id_pix
		INTEGER(I4B) :: ordering
		INTEGER(I4B) :: nside

		INTEGER(I4B) :: xsize, ysize

		REAL(SP), DIMENSION(0:12*nside**2-1) :: map

		REAL(SP),     DIMENSION(0:xsize-1,0:xsize/2-1) :: image
		LOGICAL(LGT), DIMENSION(0:xsize-1,0:xsize/2-1) :: mask

		REAL(SP) :: lon0
		REAL(DP) :: lon0rad
		REAL(DP) :: DtoR = PI/180.0_dp

		INTEGER(I4B) :: xc, dx, yc, dy

		REAL(DP), DIMENSION(:,:), ALLOCATABLE :: u, v
		REAL(DP), DIMENSION(:,:), ALLOCATABLE :: out_lat, out_lon

		INTEGER(I4B) :: status

		! --------------------------------------------------------------

		ysize = xsize/2

		!--- generates the (u,v) position on the mollweide map ---
		image(:,:) = 1.                         ! white background ???

		xc = (xsize-1)/2
		dx = xc
		yc = (ysize-1)/2
		dy = yc

		!--- allocate memory for arrays ---
		ALLOCATE(u(0:xsize-1,0:ysize-1),stat = status)
		call assert_alloc(status,'map2gif','u')
		do i = 0,xsize-1
		   u(i,:) = real(i,dp)
		end do
		ALLOCATE(v(0:xsize-1,0:ysize-1),stat = status)
		call assert_alloc(status,'map2gif','v')
		do i = 0,ysize-1
		   v(:,i) = real(i,dp)
		end do
		ALLOCATE(out_lat(0:xsize-1,0:ysize-1),stat = status)
		call assert_alloc(status,'map2gif','out_lat')
		ALLOCATE(out_lon(0:xsize-1,0:ysize-1),stat = status)
		call assert_alloc(status,'map2gif','out_lon')

		u =  2._dp*(u - real(xc,dp))/(real(dx,dp)/1.02_dp)   ! 1.02 : fudge factor in [-2,2]*1.02
		v =    (v - real(yc,dp))/(real(dy,dp)/1.02_dp)   ! in [-1,1] * 1.02

		!--- for each point on the mollweide map looks for the corresponding
		!    position (lon, lat) on the sphere ---
		mask = (u**2/4._dp + v**2) <= 1._dp
		out_lat = 1.e12_dp ! values for points out of the sphere
		out_lon = 1.e12_dp
		lon0rad = real(lon0,dp) * DtoR
		WHERE( mask )
		   out_lat = PI/2._dp - (ASIN ( 2._dp/PI * (&
			& ASIN(v) + v*SQRT( (1._dp-v)*(1._dp+v) )) ))
		   ! colat in [0,pi]
		   out_lon = -lon0rad - PI/2._dp * u/MAX(SQRT( (1._dp-v)*(1._dp+v) ),1.e-6_dp)
		   ! lon in [-pi,pi], the minus sign for astro convention
		END WHERE

		WHERE(out_lon < 0._dp)
		   out_lon = out_lon  + 2._dp*PI
		   ! lon in RAD in [0,2pi]
		END WHERE

		!--- converts the position on the sphere into pixel number and project the
		!    corresponding data value on the map ---
		DO I = 0,xsize-1
		   DO J = 0,ysize-1
		 IF(ABS(out_lat(I,J)) .LE. 1.e5_dp)THEN ! keep only meaningful points
			 if (ordering .eq. 1) then
			    call ang2pix_ring( nside, out_lat(I,J), out_lon(I,J), id_pix )
			    ! id_pix in [0,12*nside^2-1]
			 else
			    call ang2pix_nest( nside, out_lat(I,J), out_lon(I,J), id_pix )
			    ! id_pix in [0,12*nside^2-1]
			 end if
			 image(I,J) = map(id_pix)

		 END IF
		   END DO
		END DO

		!--- deallocate memory for arrays ---
		deallocate( u )
		deallocate( v )
		deallocate( out_lat )
		deallocate( out_lon )

	END SUBROUTINE PROJ_MOLLW

END MODULE mapOutput
