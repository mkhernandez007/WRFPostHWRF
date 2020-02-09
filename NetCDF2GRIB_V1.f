C========================================================================
C                        NetCDF to Grib 1 writer
C========================================================================
C
C  This program is to write data from the WRF-NMM (or HRS) netcdf output
C  into Grib 1 format.  The output file created by this program will have 
C  the same name as the input netcdf file, but with the extension of .gb.
C  This will take and convert Pressure, u- and v- component wind, 
C  temperature, specific humidity, sea-land mask, accumulated precip
C  amount, albedo, ground roughness, skin temperature, and mean sea level
C  pressure.  With this information this code calculates geopotential 
C  height and potential temperature.
C
C  Updates:
C    10/31/08: Michael Kevin Hernandez wrote the code to read in an ascii
C              file and make a grib file as a result.
C    11/07/08: Erica Collura and Shelly Guare wrote the code to write out
C              the NetCDF data into the ascii file.
C    11/23/08: Michael Kevin Hernandez combined both codes and created 
C              the two subroutine to create geopotential height and 
C              potential temperature.
C    12/01/08: Dr. Eugene Clothiaux created a small interpolation code
C              using closest neighbors approach.
C
C  Note:
C    Without the aid of Erica's and Shelly's NetCDF to ascii, this
C    program will not have been accomplished in a timely matter. With 
C    time more interpolation options will become available. 
C
C  Grib subroutines:
C    The subroutines that were used in this converter came from National
C    Centers for Environmental Precition, which can be found at thier 
C    website:
C               http://www.nco.ncep.noaa.gov/pmb/docs/on388/
C    On their website they provide a file called w3lib-1.6.tar, which is
C    a library containing GRIBI encoder/decoder and utility routines.
C
C  NetCDF subroutines:
C    The subroutines that were used in this converter came from Unidata,
C    which can be found at thier website:
C               http://www.unidata.ucar.edu/software/netcdf/
C    On their website they provide a netcdf.tar.Z, which is a library
C    containing NetCDF encoder/decoder and utility routines.
C
C------------------------------------------------------------------------
C List of argument
C------------------------------------------------------------------------
C
C Z0      = Surface Hieght
C ALBASE  = Base Albedo 
C MSLP    = Mean Sea Level Pressure
C PINT    = Pressure
C T       = Temperature
C U       = U-wind component
C V       = V-wind component
C Q       = Specific Humidity
C SM      = Sea and Land Mask
C TGROUND = Skin Temperature
C ACPREC  = Accumulated Precipitation
C GEOPOT  = PHI = Geopotential Hieght
C THETA   = Potential Temperature
C
C========================================================================
C @ program: NetCDF2GRIB.f
C @ author : Michael Kevin Hernandez, Erica Collura, Shelly Guare, and
C            Dr. Eugene Clothiaux
C @ email  : mkh182@psu.edu, elc5046@psu.edu, sag5085@psu.edu, 
C            cloth@meteo.psu.edu
C @ compile: pgf90 -o NetCDF2GRIB NetCDF2GRIB_V1.f [GRIBI library path] 
C            [NetCDF library path] 
C @ running: ./NetCDF2GRIB netcdfouputfile
C @ output : netcdfoutputfile.gb
C @ version: 1.0.0
C========================================================================
      PROGRAM NetCDF2GRIB

        IMPLICIT NONE
        
C
CC Must include the NetCDF Include file so that the program recognizes 
CC the NetCDF function and subroutine calls.
C

        include 'netcdf.inc'

C
CC NetCDF usage varaibles
C

        INTEGER*4 ncid, ncrcode, ndims, nvars, ngatts, recdim, varid 
        INTEGER*4 start(7), vartype, nvdims, nvatts, vdims(7), ivc, ill
        INTEGER*4 dimsize(7), count(7), status, dimid

        CHARACTER(len=16),ALLOCATABLE, dimension(:) :: varskeep
        CHARACTER(len = 256) dimname(7)

C
CC GRIB usage variables
C

        INTEGER*4, dimension(200) :: KPDS, KGDS
        LOGICAL*1, allocatable :: LB (:)
        REAL*4,    allocatable :: F(:)
        REAL*4,    allocatable :: Ftmp(:)
        INTEGER*4 fid, LUGB, KF, iexists, iret
        

C
CC For reading in variable
C

        INTEGER*4 nargs,IARGC, file1, file2, numbvar
        INTEGER*4 numvars, londim, latdim, modelhieght, time
        INTEGER*4 dim1, dim2, dim3, i, j, k, l, zed, numptsphi
        INTEGER*4 numpts, minpt, maxpt, minpt2, maxpt2, onelevel
        CHARACTER(LEN=256) ncfile, newname
        CHARACTER(LEN=16)  var
        REAL*4, ALLOCATABLE, dimension(:) :: glat, glon
        INTEGER*4, ALLOCATABLE, dimension(:) :: imap
        REAL*4, ALLOCATABLE, dimension(:) :: datum, pres, temp    
        REAL*4    glatmax, glatmin, glonmax, glonmin

C
CC In order to get the date of the file
C

        CHARACTER*4 yy, mo, dd, hh, mm
        INTEGER*4 iyy, imo, idd, ihh, imm

C
CC Variable you want to be gribbed.
C
        numbvar = 11
        ALLOCATE(varskeep(numbvar))

        varskeep(1)  = 'Z0'
        varskeep(2)  = 'ALBASE'
        varskeep(2)  = 'MSLP'
        varskeep(3)  = 'PINT'
        varskeep(4)  = 'T'
        varskeep(5)  = 'U'
        varskeep(6)  = 'V'
        varskeep(8)  = 'Q'
        varskeep(9)  = 'SM'
        varskeep(10) = 'TGROUND'
        varskeep(11) = 'ACPREC'
        
C
CC Retrieving a file
C 
        nargs = IARGC()
        if (nargs < 1) then
          write(*,*) 'You must enter an NetCDF file:'
          stop
        end if

        CALL GETARG(1, ncfile)

C
CC Retrieving the date of the forecast for which it is valid for.
C

        open(unit=130,file='date')
        yy = ncfile(12:15)
        mo = ncfile(17:18)
        dd = ncfile(20:21)
        hh = ncfile(23:24)
        mm = ncfile(26:27)
        write (130,'(5A7)') yy, mo, dd, hh, mm
        write (130,'(A30,A3)') ncfile,".gb"
        close (130)
        open(unit=130,file='date')
        read(130,'(5I7)')  iyy, imo, idd, ihh, imm
        read(130,'(A50)', end = 10) newname
10      write(*,'(A17,I3,A1,I2,A1,I4,A1,I2,A3)') " Gribbing date : ",  
     &       imo,"/",idd,"/",iyy,"/",ihh,":00"
        close(130)

C
CC Opening the NetCDF file.
C

        ncid = NCOPN(ncfile, NCNOWRIT, ncrcode)

C
CC Handle return if there is problems with NCOPN.
C

        if (ncrcode .ne. 0) then
           write(*,*) 'NCOPN return code:', ncrcode
           stop
        end if

C
CC NCINQ returns information about the open NetCDF file.
C

        CALL NCINQ(ncid, ndims, nvars, ngatts, recdim, ncrcode)

C
CC Handle return if problems occur from NCINQ.
C

        if (ncrcode .ne. 0) then
           stop
        else 
           continue
        end if

C
CC NCQINQ returns the name and size of a dimension.
CCC Loop though the seven dimenions of each variable extracting the
CCC dimension ID, name, and size. 
C
       
        do i = 1, ndims, 1
           dimid = i
           CALL NCDINQ(ncid, dimid, dimname(dimid), dimsize(dimid),
     &              ncrcode)
        end do

C
CC NCVINQ returns information about a NetCDF variable.  The information
CC returned is the name, type, number of dimensions, a list of dimension
CC IDs describing the shape of the variable, and the number of variable
CC attributes that have been assigned to the variable.
C

        do i=1, nvars, 1
           varid = i
           CALL NCVINQ (ncid, varid, var, vartype, nvdims, vdims,
     &                 nvatts, ncrcode)
        end do

C
CC Looping to set start equal to one each time program loops through for
CC each number of dimensions.
C

        do i=1, ndims, 1
           start(i) = 1
        end do

C
CC Open Grib1 file and setting LUGB
C

        LUGB = 50
        call baopenw(LUGB,newname,iret)

C
CC Open the interpolation mapping file.  If it does not exist,
CC then create it and save it for future use.
C

        fid = 60

        OPEN(UNIT=fid, FILE='Mapping2CartesianLatLonGrid.data',
     1       STATUS='OLD', ACTION='READ', IOSTAT=iexists)

C
CC Looping to collect information to form the KGDS array.
C

        IF (iexists .NE. 0) THEN

          DO i=1, nvars, 1

            varid = i

            CALL NCVINQ (ncid, varid, var, vartype, nvdims, vdims,
     &                   nvatts, ncrcode)

C
CC Looping through count, it is the number of dimensions of the specified 
CC variable.  
C 

            count = 1
         
            do j=1, nvdims, 1
               count(j) = dimsize(vdims(j))
            end do

C
CC Setting count to their normal dimensions
C

            londim = count(1)
            latdim = count(2)

           numpts = londim * latdim 

            if      (var == 'GLAT') then
               ALLOCATE(glat(numpts))
               CALL NCVGT(ncid, varid, start, count, glat, ncrcode)
            else if (var == 'GLON') then
               ALLOCATE(glon(numpts))
               CALL NCVGT(ncid, varid, start, count, glon, ncrcode)
C
CC             Creating a simple near-neighbor interpolation map.
C
               ALLOCATE(imap(numpts))

               CALL InterpMapCreate(latdim, londim, numpts, glat, glon,
     1              glatmin, glatmax, glonmin, glonmax, imap)

               DEALLOCATE(glat)
               DEALLOCATE(glon)

               GOTO 200

            ENDIF

          END DO

        ELSE

           READ(fid, "(3I10,4F16.4)") latdim, londim, numpts, 
     1                    glatmin, glatmax, glonmin, glonmax
           ALLOCATE(imap(numpts))
           READ(fid, "(I10)") imap
           CLOSE(fid)

        ENDIF

200     write(*,'(A34)') '**********************************'
        write(*,'(A24)') 'Grid has been defined   '
        write(*,'(A24)') '** Commence Gribbing    '
        write(*,'(A34)') '**********************************'

        CALL define_kgds(londim, latdim, modelhieght,
     1                   glonmax, glatmax, glonmin, glatmin, KGDS)

C
CC Looping to collect information about each variable.
C

        do i=1, nvars, 1

          varid = i

          CALL NCVINQ (ncid, varid, var, vartype, nvdims, vdims,
     &               nvatts, ncrcode)
 
          count = 1
         
          do j=1, nvdims, 1
             count(j) = dimsize(vdims(j))
          end do

C
CC Setting count to their normal dimensions
C

          londim = count(1)
          latdim = count(2)
          modelhieght = count(3)
          time =  count(4)
          dim1 =  count(5)
          dim2 =  count(6)
          dim3 =  count(7)

          numpts = londim * latdim * time * modelhieght

C
CC NCVGT reads an array of values from NetCDF variable from the file.
CCC Looping to recognize each variable. The if statement recongnizes 
CCC each of the variabels listed in varkeep(i).
C

          do k = 1, numbvar, 1
             if (var .eq. varskeep(k)) then
                ALLOCATE(datum(numpts))
                CALL NCVGT(ncid, varid, start, count, datum, ncrcode)

                if (var == 'Q') then
                   datum = (datum * 1000)
                   call ToGrib(numpts, modelhieght, datum, iyy, imo, 
     &                         idd, ihh, imm, var, imap, KGDS)
                   DEALLOCATE(datum)
                else if (var == 'PINT') then
                   ALLOCATE(pres(numpts))
                   pres = 0
                   pres = datum
                   call ToGrib(numpts, modelhieght, datum, iyy, imo,
     &                         idd, ihh, imm, var, imap, KGDS)
                   DEALLOCATE(datum)
                else if (var == 'T') then
                   ALLOCATE(temp(numpts))
                   temp = 1
                   temp = datum
                   numptsphi = numpts
                   call ToGrib(numpts, modelhieght, datum, iyy, imo,
     &                         idd, ihh, imm, var, imap, KGDS)
                   DEALLOCATE(datum)
                else
                   call ToGrib(numpts, modelhieght, datum, iyy, imo,
     &                         idd, ihh, imm, var, imap, KGDS)

                   DEALLOCATE(datum)
                end if
             end if
          end do
        end do

C
CC The subroutine to calculate both geopotential and potential 
CC temperature.
C

        call thetas(temp, pres, numptsphi, modelhieght, iyy, imo, idd,
     &              ihh, imm, var, imap, KGDS)

        call geopotential(temp, pres, numptsphi, modelhieght, iyy, imo, 
     &                    idd, ihh, imm, var, imap, KGDS)


        call baclose(LUGB,iret)

        STOP
     
      END

C========================================================================
C Subroutine theta
C========================================================================
C This subroutine calculates potential temperature and sets up the appro-
C priate theta array to be sent into togrib subroutine which gribs the 
C data.
C========================================================================
      subroutine thetas(temp, pres, numpts, modelhieght, iyy, imo, idd,
     &                  ihh, imm, var, imap, KGDS)

        implicit none

        INTEGER*4  numpts, onelevel, modelhieght, g, Rd, minpt, maxpt
        INTEGER*4  iyy, imo, idd, ihh, imm, minpt2, maxpt2, j, alpha, cp
        CHARACTER(len=16) var
        INTEGER*4  imap(numpts/modelhieght)
        INTEGER*4, dimension(200) :: KGDS
        REAL*4, dimension(numpts) :: temp, theta
        REAL*4, dimension(numpts/42*43) :: pres

        g = 9.81
        Rd = 287.15
        cp = 1952

        var = 'THETA'
        modelhieght = 42
        onelevel = numpts/modelhieght
        alpha = 2/(numpts+1)

C
CC After the first level the hydrostatic equation provides the
CC algoritm to calculate this variable.
C

        do j = 1, modelhieght, 1                    ! modelhieght, 1
           minpt = ((j - 1) * onelevel) + 1
           minpt2 = (j*onelevel) + 1
           maxpt = (j * onelevel)
           maxpt2 = (j+1)*onelevel
           theta(minpt:maxpt) = temp(minpt:maxpt)*((((pres(minpt:maxpt)
     1                          +(1-alpha)*pres(minpt:maxpt))/(2-alpha))
     2                          /(10000))**(Rd/cp))
        end do



        call ToGrib(numpts, modelhieght, theta, iyy, imo, idd, ihh, imm,
     &              var, imap, KGDS)

      end subroutine

C========================================================================
C Subroutine geopotential
C========================================================================
C This subroutine calculates geopotential hieght and sets up the appro-
C priate phi array to be sent into togrib subroutine which gribs the
C data.  The data in phi is the sum of the change of phi and the previous
C value of phi.
C========================================================================

      subroutine geopotential(temp, pres, numpts, modelhieght, iyy, imo,
     &                        idd, ihh, imm, var, imap, KGDS)

        implicit none

        INTEGER*4  numpts, onelevel, modelhieght, g, Rd, minpt, maxpt
        INTEGER*4  iyy, imo, idd, ihh, imm, minpt2, maxpt2, j
        INTEGER*4  minpt0, maxpt0
        CHARACTER(len=16) var
        INTEGER*4  imap(numpts/modelhieght)
        INTEGER*4, dimension(200) :: KGDS
        REAL*4, dimension(numpts) :: temp, phi
        REAL*4, dimension(numpts/42*43) :: pres
  
        g = 9.81
        Rd = 287.15

        var = 'GEOPOT'
        modelhieght = 42        
        onelevel = numpts/modelhieght

C
CC Before the first level the hydrostatic equation provides the 
CC algoritm to calculate this variable. To get geopotential (m^2/s^2)
CC into geopotential in meters alone, you must devide by gravity.
C

        phi(1:onelevel) = 0

        do j = 2, modelhieght, 1                    ! modelhieght, 1
           minpt0 = ((j - 2) * onelevel) + 1
           minpt = ((j - 1) * onelevel) + 1
           minpt2 = (j*onelevel) + 1
           maxpt0 =(j-1)*onelevel
           maxpt = (j * onelevel)
           maxpt2 = (j+1)*onelevel
           phi(minpt:maxpt) = ((Rd * temp(minpt:maxpt)*
     1                        (alog(pres(minpt:maxpt)/
     2                             pres(minpt2:maxpt2))))/g)
     3                        + phi(minpt0:maxpt0)
        end do

        call ToGrib(numpts, modelhieght, phi, iyy, imo, idd, ihh, imm,
     &              var, imap, KGDS)

      end subroutine
C========================================================================
C Subroutine ToGrib
C========================================================================
C This subroutine recieves in the data file, creates the KPDS, and splits
C the datum array (the array that holds the data) into multiple layers,
C dictated by the hieght dimension in the netcdf file.  
C========================================================================
      subroutine ToGrib(numpts, modelhieght, datum, iyy, imo, idd, ihh, 
     &                  imm, var, imap, KGDS)

        implicit none

        INTEGER*4  onelevel, numpts, modelhieght, j, zed, KF, minpt
        INTEGER*4  maxpt, LUGB, iyy, imo, idd, ihh, imm, iret
        CHARACTER(len=16) var

        INTEGER*4, dimension(200)   :: KPDS, KGDS
        INTEGER*4  imap(numpts/modelhieght)

        LOGICAL*1, allocatable      :: LB (:)
        REAL*4,    allocatable      :: F(:), Finterp(:)
        REAL*4,    dimension(numpts):: datum

        write(*,*) '** Gribbing variable:  ',var

        LUGB = 50
        onelevel = numpts/modelhieght
        KF = onelevel

        ALLOCATE (F(onelevel))
        ALLOCATE (Finterp(onelevel))
        ALLOCATE (LB(onelevel))

        LB = .false.

        do j = 1, modelhieght, 1                    ! modelhieght, 1
           KF = onelevel
           zed = j
           if (zed == 43) then
              goto 210
           end if
C           write(*,'(2A24)')    ' variable             :', var
C           write(*,'(A24,I10)') ' Current Height       :', zed

           call define_kpds(zed, iyy, imo, idd, ihh, imm, var, KPDS)

           minpt = ((j - 1) * onelevel) + 1
           maxpt = (j * onelevel)
           
           F = datum(minpt:maxpt)

           CALL InterpMapApply(onelevel, imap, F, Finterp)

           CALL putgb(LUGB,KF,KPDS,KGDS,LB,Finterp,iret)
C          CALL irets(iret)
 
       enddo

210    continue
       DEALLOCATE(LB)
       DEALLOCATE(F)
       DEALLOCATE(Finterp)

      end subroutine
C========================================================================
C Subroutine define_kpds
C========================================================================
C This subroutine creates the KPDS array, which describes the variable as
C it is.
C========================================================================
      subroutine define_kpds(zed, iyy, imo, idd, ihh, imm, var, kpds)

        implicit none      

        CHARACTER*16 var
        INTEGER*4, dimension(200) :: KPDS
        INTEGER*4 iyy, imo, idd, ihh, imm, zed, kpds5 

C
CC  The below code is a manual hard coding of the KPDS array needed for 
CC  this program to save the ascii file into a Grib file.
C

        KPDS(1)  = 51   ! Originating Center is Miami
        KPDS(2)  = 112  ! WRF-NMM, generic resolution (HRS is a derivative of
                        ! of WRF-NMM)
        KPDS(3)  = 255  ! non-standard grid - defined in the GDS 
        KPDS(4)  = 00000001   ! Section 2 included, Section 3 omitted
        KPDS(5)  = kpds5(var) ! function that assings units to variables
        KPDS(6)  = 107  ! 107 for sigma level
        KPDS(7)  = zed  ! if KPDS(6) = 107 then zed which is the model 
                        ! number
        KPDS(8)  = iyy  ! Year
        KPDS(9)  = imo  ! Month 
        KPDS(10) = idd  ! Day
        KPDS(11) = ihh  ! Hour 
        KPDS(12) = imm  ! Minute
        KPDS(13) = 10   ! Unit of time is 3 hour
        KPDS(14) = 0    ! Time Range 1
        KPDS(15) = 0    ! Time Range 2
        KPDS(16) = 0    ! Forecast product valid for reference time 
        KPDS(17) = 0    ! Number included in an average (no average)
        KPDS(18) = 1    ! Version of Grib  
        KPDS(19) = 2    ! GRIB parameter table version
        KPDS(20) = 0    ! Number missing from averages or accumulations
        KPDS(21) = 21   ! Current Century of data 
        KPDS(22) = 0    ! Unit scaling power of 10
        KPDS(23) = 0    ! Hurricane Research Division doesn't have a subcenter 

        return

      end subroutine
C========================================================================
C Subroutine irets
C========================================================================
C In togrib subroutine, there is a certain feature that if this file 
C doesn't grib the data, the error messages can be called onto the screen
C so that they can be fixed immediately.
C========================================================================
      subroutine irets(iret)

        implicit none

        integer*4 iret

        if (iret == 0) then
           write(*,*) "COMPLETED MAKING GRIB FIELD WITHOUT ERROR"
        else if (iret == 1) then
           write(*,*) "IPFLAG NOT 0 OR 1"
        else if (iret == 2) then
           write(*,*) "IGFLAG NOT 0 OR 1"
        else if (iret == 3) then
           write(*,*) "ERROR CONVERTING IEEE F.P. NUMBER TO IBM370 F.P."
        else if (iret == 4) then
           write(*,*) "W3FI71 ERROR/IGRID NOT DEFINED"
        else if (iret == 5) then
           write(*,*) "W3FK74 ERROR/GRID REPRESENTATION TYPE NOT VALID"
        else if (iret == 6) then
           write(*,*) "GRID TOO LARGE FOR PACKER DIMENSION ARRAYS"
        else if (iret == 7) then
           write(*,*) "LENGTH OF BIT MAP NOT EQUAL TO SIZE OF FLD/IFLD"
        else if (iret == 8) then
           write(*,*) "W3FI73 ERROR, ALL VALUES IN IBMAP ARE ZERO"
        end if
  
        return

      end subroutine

C========================================================================
C Subroutine define_kgds
C========================================================================
C This subroutine creates the KGDS array, which describes the variable's 
C and its value onto a rotated latitude longitude grid.
C========================================================================
      subroutine define_kgds(londim, latdim, modelhieght, glonmax,
     &                       glatmax, glonmin, glatmin, kgds)

        implicit none

        INTEGER*4, dimension(200) :: KGDS
        INTEGER*4 londim, latdim, modelhieght, total
        REAL*4 glonmax, glatmax, glonmin, glatmin, dlat, dlon, cenlat, 
     &         cenlon
     

C
CC Total number of actual points
C   

        total = londim * latdim

C
CC Calculating the central lat and lon in m^o
C

        cenlon = (-72.751+360) * 1000
        cenlat = 21.0 * 1000

C
CC Resolution in m^o degrees
C

        dlat  = abs((glatmax-glatmin)/latdim)
        dlon  = abs((glonmax-glonmin)/londim)

C
CC  The below code is a manual hard coding of the KGDS array needed for
CC  this program to save the ascii file into a Grib file.
C

        KGDS(1)  = 0        ! regular lat lon
        KGDS(2)  = londim   ! Number of points on Latitude    
        KGDS(3)  = latdim   ! Number of points on Longitude
        KGDS(4)  = glatmin  ! latitude of first grid point
        KGDS(5)  = glonmin  ! longitude of first grid point
        KGDS(6)  = 136      ! Direction increments given, Earth assumed
                            ! spherical, reserved, reserved, u and v are
                            ! resolved relative to the defined direction 
                            ! of increasing x and y coordinates 
                            ! respectively
        KGDS(7)  = glatmax  ! latitude of extreme point
        KGDS(8)  = glonmax  ! longitude of extreme point
        KGDS(9)  = dlon     ! longitudinal direction of increment
        KGDS(10) = dlat     ! latitudinal direction increment
        KGDS(11) = 64       ! points scanned in +i direction, +j
                            ! direction, Adjacent points in i are
                            ! consecutive, reserved
        KGDS(19) = 0        ! number for vertical coordinate parameters
        KGDS(20) = 255      ! no octet number of the list of vertical 
                            ! coordinate parameters
        KGDS(21) = 0        ! no PL
        KGDS(22) = 0        ! number of words in each row

        return

      end subroutine
C========================================================================
C Subroutine kpds5
C========================================================================
C This function will assign the number of the particular  variable that 
C must be sent into the KPDS array inorder for the variable to be gribbed
C correctly.  These values come from a NCEP provided table.
C========================================================================
      integer function kpds5(var)
      
         implicit none 

         character*16 var

         if (var == 'PINT') then
            kpds5 = 1           ! Pressure [Pa]
         else if (var == 'MSLP') then
            kpds5 = 2           ! Pressure reduce to MSL [Pa]
         else if (var == 'GEOPOT') then
            kpds5 = 7           ! Geopotential [m]
         else if (var == 'T') then
            kpds5 = 11          ! Temperature[K]
         else if (var == 'THETA') then
            kpds5 = 13          ! Potential tempertuare [k]
         else if (var == 'U') then
            kpds5 = 33          ! u-wind component [m/s]
         else if (var == 'V') then
            kpds5 = 34          ! v-wind component [m/s]
         else if (var == 'W') then
            kpds5 = 40          ! w-wind component [m/s]
         else if (var == 'Q') then
            kpds5 = 51          ! specific humidity[g/kg] 
         else if (var == 'ACPREC') then
            kpds5 = 61          ! total precip
         else if (var == 'SM') then
            kpds5 = 81          ! Land [=0] Sea [=1] Mask
         else if (var == 'Z0') then
            kpds5 = 83          ! Roughness height
         else if (var == 'ALBASE') then
            kpds5 = 84          ! Base albedo
         else if (var == 'TGROUND') then
            kpds5 = 85          ! Deep ground soil temperature
         else             
            kpds5 = 255         ! Missing value
         end if

         return

      endfunction

C========================================================================
C Subroutine InterpMapCreate
C========================================================================
C This subroutine creates a file with the calculations to be used in 
C interpolating the code onto a regular latitude and longitude grid.
C========================================================================

      SUBROUTINE InterpMapCreate(nlats, nlons, npts, lat, lon, 
     1             latmin, latmax, lonmin, lonmax, imap)

        IMPLICIT NONE

        INTEGER*4 nlats, nlons, npts
        REAL*4 lat(npts), lon(npts)
        REAL*4 latmin, latmax, lonmin, lonmax
        INTEGER*4 imap(npts)

        INTEGER*4 fid, i, j, ilat, ilon
        REAL*4 dlatnew, dlonnew, dist, distmin
        REAL*4, ALLOCATABLE, dimension(:) :: latnew, lonnew

        latmin = -1000.
        DO ilon = 1, nlons, 1
          IF (latmin .LT. lat(ilon)) THEN
            latmin = lat(ilon)
          END IF
        END DO

        latmax = 1000.
         DO ilon = (nlats-1)*nlons + 1, nlats*nlons, 1
          IF (latmax .GT. lat(ilon)) THEN
            latmax = lat(ilon)
          END IF
        END DO

        lonmin = -1000.
        DO ilat = 1, (nlats-1)*nlons + 1, nlons
          IF (lonmin .LT. lon(ilat)) THEN
            lonmin = lon(ilat)
          END IF
        END DO

        lonmax = 1000.
        DO ilat = nlons, nlats*nlons, nlons
          IF (lonmax .GT. lon(ilat)) THEN
            lonmax = lon(ilat)
          END IF
        END DO

        ALLOCATE(latnew(npts))
        ALLOCATE(lonnew(npts))

        dlatnew = (latmax - latmin) / (nlats - 1.)
        dlonnew = (lonmax - lonmin) / (nlons - 1.)

        i = 1
        DO ilat = 1, nlats, 1
          DO ilon = 1, nlons, 1
            latnew(i) = latmin + (ilat-1) * dlatnew
            lonnew(i) = lonmin + (ilon-1) * dlonnew
            i = i + 1
          END DO
        END DO

        DO i = 1, npts, 1
          distmin = 1000.
          DO j = 1, npts, 1
            dist = SQRT((lon(j)-lonnew(i))**2+(lat(j)-latnew(i))**2)
            IF (dist .LT. distmin) THEN
              distmin = dist
              imap(i) = j
            END IF
          END DO
C         WRITE(*,"(I8,A5,I8)") i, ' ==> ', imap(i)
        END DO

        latmin = latmin * 57.2957795 * 1000.  ! 1 radian = 57.2957795^o
        latmax = latmax * 57.2957795 * 1000.  ! Convert to millidegrees
        lonmin = (lonmin * 57.2957795 + 360.) * 1000.
        lonmax = (lonmax * 57.2957795 + 360.) * 1000.

        fid = 60
        OPEN(UNIT = fid, FILE = 'Mapping2CartesianLatLonGrid.data')
        WRITE(fid, "(3I10,4F16.4)") nlats, nlons, npts,
     1                              latmin, latmax, lonmin, lonmax
        WRITE(fid, "(I10)") imap
        CLOSE(fid)

      END SUBROUTINE InterpMapCreate

C========================================================================
C Subroutine InterpMapApply
C========================================================================
C This subroutine applies the file that was created to the data to make 
C the interpolation into the regular latitude and longitude grid work.
C========================================================================
      SUBROUTINE InterpMapApply(npts, imap, F, Finterp)

        IMPLICIT NONE

        INTEGER*4 npts, imap(npts)
        REAL*4 F(npts), Finterp(npts)

        INTEGER*4 i

        DO i = 1, npts, 1
          Finterp(i) = F(imap(i))
        END DO

      END SUBROUTINE InterpMapApply

C========================================================================

