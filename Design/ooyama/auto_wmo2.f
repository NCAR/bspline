C 
C     Program for real-time and post-flight editing of dropsonde data.  
C     Written by James L. Franklin, NOAA/Hurricane Research Division.
C
C     Compiles and loads with: fcl_edits.cmd.
C
C     Requires modules: lib/tempdrop.f, sondelib.f, plotlib.f, skewt.f,
C                       synoptic_map.f, avapslib.f
C     -----------------------------------------------------------------
C
C
C     REVISION HISTORY
C     -----------------------------------------------------------------
C
C     01/09/97    V1.35   - Adds 10166,67 groups to TEMP-DROP message.
C                         - 'Doubtful' flag available for aircraft use.
C
C     12/23/96    V1.34   - Changed AVGDEPTH in DYNTCORR from 5 to 2 mb
C                           to reduce the introduction of super-adiab.
C                           lapse rates just below sharp inversions.
C                         - Adds flag "S" for subjective edit (GND).
C                         - Changed XSD for PR from 3.5 to 4.5 (was
C                           flagging too many points in short drops.
C                         - Bug fix in XXAA that causes bad heights for
C                           mandatory level falling between PRL and 
C                           first sonde level.
C                         - Bug fix in ODEND: failed if no PTH.
C
C     12/19/96    V1.33   - Interpolates VVEL for GPS sondes
C                         - Changes PRGAPMAX to 50 for ground use.
C                         - DSOND replaced with AUTOQC.
C                         - Rings bell if old files present, P bias.
C                         - Flush NEEDOPT for skipped graphics.
C                         - Post-filter check ensures RH>=0.1%
C                         - LVLINT, SNDNG, PRINTAC can write to
C                           log file only, if desired.
C                         - PBIAS now works for ODWs also.
C
C     12/17/96    V1.32   - Revised format for skew_ctrl.dat.
C                         - Displays possible flight GA/PS errors for
C                           upward height integrations.
C                         - Data flagging permitted by time or pressure.
C
C     12/13/96    V1.31   - Displays flags in LVLINT. 
C                         - Resets flags for missing data to '0'. 
C                         - Moves PRINTAC, LVLINT, SNDNG, CFLG, CXFLG,
C                           INDEX to lib/avapslib.f
C                         - Adds NOPLOTS to turn off graphics for NCAR.
C                         - READGSONDE fix for gaps at beginning.
C                         - NEEDOPT(9) turned on for dynamic p fix.
C                         - Bug fix in TEMP-DROP decoder DROP.
C                         - Refinements in synoptic charts.
C                         - Minor bug fix in STATUS.
C                         - Default in PTHADD is 'P'
C                         - SYNOPTIC_MAP now requires LUG passed.
C
C     11/21/96    V1.30   - FASTEX VERSION
C                         - Improves dynamic T correction to reduce 
C                           noise.
C                         - Dynamic correction a default in CWORD.
C                         - READGSONDE should read new AVAPS VER line,
C                           toss CRC bit errors as directed for A/C
C                           use, and read # GPS satellites, RH1 and RH2.
C                         - # satellites inserted in WQ slot.
C                         - Fixes bug in XDATZ if all data bad.
C
C     11/19/96    V1.29   - Adds COMMON PROCESSING to provide record of
C                           all processing steps.
C                         - NEEDOPT data now saved in SAVEDATA.
C                         - Corrects bug which prevented any PTH
C                           filtering if any one parameter failed.
C                         - Displays processing record before terminating
C                           (GND only).
C                         - In SMOOTHSPL, filtered values not used at
C                           first and last points in sounding.
C
C     11/18/96    V1.28   - Adds option to align wind profile w/PTH.
C                         - Fast-fall detect warning.
C                         - Winds flagged for fast falls.
C                         - Can specify all pts w/"A" in STATUS.
C                         - Can flag data from RH sensors 1 or 2.
C                         - Revamp STATUS interface.
C                         - READGSONDE: all vvel tossed < 960928.
C                         - User can flag data doubtful (GND).
C                         - Long interpolations flagged doubtful (GND).
C                         - Cuts RH ambient time if T > -30 C.
C
C     11/04/96    V1.27   - Adds 51515 Additional Data Groups 10190 and
C                           10191: extrapolated heights and sfc press
C                           to parts A and B of WMO message (tempdrop.f)
C                         - Initial RH flagging based on time constant
C                         - 10190 groups are decoded
C                         - Command for missing heights now 'M' (was X)
C                         - Sig levels not plotted on synoptic charts
C
C     11/01/96    V1.26   - Adds dynamic temperature correction option
C                         - Initial T flagging based on time constant
C
C     10/23/96    V1.25   - Automatic encoding on hardwire switch
C                         - Moves pds files to subdirectory AVAPS_fnl
C                         - Raw AVAPS files now in AVAPS_raw (gnd only)
C                         - Some file subroutines moved to avapslib
C
C     10/14/96    V1.24   - Fixes erroneous hardwire wind flag
C                         - HEIGHTS uses default RH in the air also
C                         - Improved skew-t AUTOPAR for PRL<200
C                         - Tightened buddy check on RH for CRC errors
C                         - Est pr calculated in READGSONDE
C                         - Adds VVEL QC on GPS winds
C
C     10/10/96    V1.23   - If no PTH, ODEND looks for last wind
C                         - Check if all data are bad in DSOND
C                         - Trivial change to skewt
C
C     10/7/96     V1.22   - Adds comment line to .pds file
C                         - Refined equilibration flagging
C                         - Can't change flag for missing data to "K"
C                         - Option to skip QC
C                         - Change PRGAPMAX
C                         - Estimated pressure for GPS sondes
C
C     10/04/96    V1.21   - Fix in SYNOPTIC_MAP if no valid data
C                         - Move EDIT data from opt 8 to opt 6
C                         - Can exit program with Q as well as 16
C                         - Restores pause option in SNDNG
C
C     09/27/96    V1.20   - Check for pressure biases
C                         - NEEDOPT(11) false for ground use
C                         - Check if AC record close to launch time
C                         - Check for premature launch
C                         - Doubles time for ambient adj for high drops
C                         - Skips bad lines in READGSONDE
C                           and inserts time gaps (except at start)
C                         - Overhaul of HEIGHTS: gaps are now OK
C                         - Fix bug in TEMPDROP if index of sig lvl = 999.
C                         - Plots VERT VEL on U/V plots
C                         - WQ for GPS sondes set to -999.
C
C     09/23/96    V1.19   - Default EOD is last record w/PTH data
C                         - Can select listing by time
C                         - Does not pause in SHOWLEVEL
C                         - Only write pds file to NRECS, not MXRC
C                           (No V1.18 pds files survive.)
C                         - Move TIMEPLOT subroutines to plotlib
C                         - Nicer PRINTAC
C                         - Change BC in SMOOTHSPL to fit data better
C                           Add warning for excessive smoothing (gnd only)
C                         - Change order of prompts in STATUS
C
C     09/10/96    V1.18   - Refine wind QC
C                         - Rename control files
C                         - Thin wind barbs on tall skew-t's
C                         - Option to read,write "final" .pds file
C                         - Pen fixes for djet time plots
C                         - Add NSITE as .ctl parameter
C
C     09/03/96    V1.17   - Really corrects EOD bug for ODWs
C                         - PRGAPMAX depends on NSITE
C                         - Ht calc uses deflt if RH missing (gnd only)
C                         - Displays full operator comment lines
C
C     08/26/96	  V1.16   - Corrects extra line before 31313 group
C                         - Checks for no obs before plotting syn map
C
C     08/12/96    V1.15 - Updates for ground use:
C                         - Does not pause in SNDNG
C                         - Does not show 999's at end in ODEND
C                         - Can plot U,V components
C                         - User defined pr interval in LVLINT
C                         - Prints tables to editsonde.log
C                         - Plots syn maps even on first sonde
C                         - Checks to see if launch data missing
C                         - Corrects update error in launch RH
C                         - Minor skew-t improvements (WINDBARB)
C                         - Added common block PROGCTL
C                         - Better post-launch auto-flagging
C                         - Adds wind error checking for GPS
C
C     08/09/96    V1.14   - correct bug in skewt.f/AUTOPAR
C     07/17/96    V1.13   - corrects EOD bug for ODWs
C     07/15/96    V1.12   - update to sondelib routine VICSETUP
C     05/28/96    V1.11   - adds 31313 group to message.
C
C     05/17/96	  V1.10   - initial enhanced distribution.
C     ------------------------------------------------------------------
C
C
C
C     DESCRIPTION OF PROGRAM-WIDE VARIABLES
C     ------------------------------------------------------------------
C
C     Program parameters:
C
C     MXRC   	=   Maximum number of data records in raw sonde file
C     NINVAR 	=   Number of sounding variables 
C     NOPT   	=   Number of menu options 
C     NPROC     =   Size of array for processing record
C 
C
C     Program control variables (Common block PROGCTL):
C 
C     NSITE     =   1 for ground use, 2 for aircraft (set by PREFLIGHT)
C     VERSION   =   Version of EDITSONDE
C     AUTOTEMP  =   Normally .FALSE. - set to .T. to automatically
C                   generate TEMP-DROP message.
C
C
C     Logical units (Common block PARAM):
C
C     LUT	=   Terminal
C     LUFI	=   Raw sonde input file (.avp)
C     LUFX      =   Supplemental files (skewt control data)
C     LUFW      =   Binary sonde files (.tmp, .pds)
C     LUFP      =   File for input control parameters (editsonde.ctl)
C     LUFO      =   Output files (tempdrop.dat)
C     LUPR      =   File for printed output (with carriage controls)
C     
C     Sonde characteristics (read in from editsonde.ctl):
C
C     NPLTFORM  =   Aircraft platform (42, 43, or 49)
C     NSNDTYPE  =   Type of sonde (1 = ODW, 2 = GPS AVAPS, 3 = GPS GLASS, 4 = CLASS(Omega or LORAN))
C     FALLRATE  =   Nominal fallrate of sonde (mb/s)
C     DATARATE  =   Desired sampling rate from raw sonde data (Hz)
C
C
C     Launch variables (Common block LAUNCH):
C
C     IYR	Year (yy)
C     IMO	Month (mm)
C     IDY	Day (dd)
C     RLAT	Latitude (deg N)
C     RLON	Longitude (deg W)
C     TIML	Time (hhmmss)
C     PRL	Pressure (mb)
C     TEL	Temperature (C)
C     RHL	Rel. humidity (%)
C     DPL	Dew point (C)
C     HTL	Geopotential height (m)
C     WDL	Wind direction (deg)
C     WSL	Wind speed (m/s)
C     STNS      Omega stations (A4)
C
C
C     Sounding variables (Common block SOUNDING):
C
C     SNDDAT = Array of sounding data values: 
C 
C     SNDDAT(I,1) = Time after launch (s) 
C     SNDDAT(I,2) = Pressure (mb) 
C     SNDDAT(I,3) = Temperature (C) 
C     SNDDAT(I,4) = Relative humidity (%) 
C     SNDDAT(I,5) = Estimated pressure (mb) 
C     SNDDAT(I,6) = Wind direction (deg)
C     SNDDAT(I,7) = Wind speed (m/s)
C     SNDDAT(I,8) = Wind quality (0-9)   	
C     SNDDAT(I,9) = Geopotential height (m) 
C     SNDDAT(I,10) = Vertical sonde vel (m/s) 	(G-sonde only) 
C     SNDDAT(I,11) = Latitude (deg N)    	(G-sonde only) 
C     SNDDAT(I,12) = Longitude (deg W)   	(G-sonde only) 
C 
C     ISNDFLG = Array of data status flags: 
C 
C     ISNDFLG(I,1) = Flag for pressure
C     ISNDFLG(I,2) = Flag for temperature 
C     ISNDFLG(I,3) = Flag for humidity
C     ISNDFLG(I,4) = Flag for wind
C 
C     ISNDFLG:   0 = (K)eep (data ok)
C                1 = (R)eplace with interpolation
C                2 = (M)issing
C                3 = (I)nterpolated value
C                4 = (D)oubtful or questionable accuracy
C                5 = (S)ubjectively determined
C 
C
C
C     Processing variables (COMMON block PROCESSING):
C
C     NPSTP	= Number of elements of PROC currently used
C
C     PROC = Array of processing record
C
C     PROC(1)	= PTH filter (s)
C     PROC(2)   = Wind filter (s)
C     PROC(3)   = P offset (mb)
C     PROC(4)   = T offset (C)
C     PROC(5)   = H offset (%)
C     PROC(6)   = Splash pr (mb)
C     PROC(7)   = Hydrostatic anchor (1 = FL, 2 = SFC)
C     PROC(8)   = Substituted estimated pressures (1 = YES)
C     PROC(9)   = Dynamic T correction (1 = YES)
C
C
C     ------------------------------------------------------------------
C     Program EDITSONDE begins here...
C     ------------------------------------------------------------------
C
C 
      CHARACTER*1 ANS, ANSD
      CHARACTER*4 STNS, VERSION
      CHARACTER*9 RUNSTRING, DIREC
      character*100 FILENAME, FNAME
      CHARACTER*20 INPUT
C
      LOGICAL NEEDOPT,TEMPOK,TEMPFILE,PROCFILE,AUTOTEMP,NOPLOTS
      LOGICAL DISP,PRINT
C 
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
      COMMON /OPTIONS/ NEEDOPT(NOPT)
      COMMON /PROCESSING/ NPSTP,PROC(NPROC)
C 
C
      DATA VERSION/'1.35'/
      DATA NPSTP/9/
      DATA NEEDOPT/.FALSE.,.TRUE.,.TRUE.,.FALSE.,.TRUE., 
     *             .TRUE.,.TRUE.,.FALSE.,.TRUE.,.TRUE.,
     *             .TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.,.FALSE./ 
      DATA BAD/-999./, PBIAS/0.0/, AUTOTEMP/.TRUE./,NOPLOTS/.TRUE./

	integer status

C
C
C 
C     ----------------------------------------------------------------
C     Initialize LU's, open terminal and printer files.
C     Note that library routines use lu's 11-13 and 141-144.
C     ----------------------------------------------------------------
	flag_thr = 2					!thread started

	status = SetThreadPriority(handle_thr, THREAD_PRIORITY_LOWEST)

	if (.not. status) then
		status_msg_thr = "unable to lower priority"
	endif

      LUT  = 0                          
      LUPR = 6                         
      LUFI = 121                       	
      LUFX = 122                       
      LUFO = 123                       
      LUFP = 124                       
      LUFW = 125                       
C      OPEN(LUT,FILE='con',BLOCKSIZE=10000,
      OPEN(LUT,FILE='hrd_v135.out',BLOCKSIZE=10000,
     *	ACTION='READWRITE')        
      OPEN(LUPR,FILE='editsonde.log',STATUS='UNKNOWN', 
     *     IOSTAT=IOS) 
C
C
C     Get control parameters from disk file
C     -------------------------------------
c      OPEN(LUFP,FILE='editsonde.ctl',STATUS='OLD',IOSTAT=IOS,ERR=1001) 
c	goto 1002
c1001	continue
c	status_msg_thr = 'Error opening editsonde.ctl'
c	return
c1002	continue
c      READ(LUFP,*) NSITE
c      READ(LUFP,*) NPLTFORM
c      READ(LUFP,*) NSNDTYPE
c      READ(LUFP,*) FALLRATE
c      READ(LUFP,*) DATARATE
c      CLOSE(LUFP)
	nsite = 2
	npltform = 49
	if (fileType_thr .LT. 3) then
		nsndtype = 2					
		fallrate = 0.975
		datarate = 2.0
	else
		nsndtype = fileType_thr
		fallrate = 0.5
	    if (nsndtype .EQ. 3) then
			datarate = 1.0
		else
			datarate = 0.1
		endif
	endif

C
C
C     Automatic encoding w/ minimal user interaction?
C     If so, go off, never to return...
C     -----------------------------------------------
	RUNSTRING = '         '
C	status_msg_thr = 'Initiating AVAPS_WMO'

      CALL AVAPS_WMO(RUNSTRING)

	handle_thr = 0

	return
      END 
C 
C 
C 
C     ------------------------------------------------------------ 

	 SUBROUTINE AVAPS_WMO(SONDEID)
	
C 
C     Program for automatic processing of AVAPS dropwindsonde 
C     files.  
C     ------------------------------------------------------------
C 
C 
	use thread_common

      PARAMETER (MXRC = 9000, NINVAR = 12, MAXTMPMSG = 500) 

C 
      CHARACTER*4 STNS, VERSION
      CHARACTER*9 SONDEID
C 
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
C
C
C		   
	status_msg_thr = 'reading data'

      CALL READGSONDE(NRECS,SONDEID,flname_thr) 
c
c	snddat now has input sounding
C
	npoints_thr = IXEND
	raw_data_thr = snddat

	status_msg_thr = 'data read'

c
c   if we are not doing the qc, then just put bad vals into the
c   qc data and return

	if (dontQcFlag_thr) then
		DO I = 1, MXRC 
			DO J = 1,NINVAR 
				qc_data_thr(I,J) = -999.
			enddo
		enddo
        status_msg_thr = "normal termination"
	    flag_thr = 3						! succesfull termination
        return
	endif

      IF (RLAT.LT.-90. .OR. RLON.LT.-180.) THEN
         WRITE(LUT,'(////)')
         status_msg_thr = " *** FATAL ERROR: BAD SONDE LAT/LON *** "                          
	     flag_thr = -1							!error in processing
         return
      ENDIF
C
      CALL ODEND(NRECS)                !   Set end of drop
	if (flag_thr .LT. 0) return

 	status_msg_thr = 'auto qc'
      CALL AUTOQC(PBIAS)               !   Error detection
	if (flag_thr .LT. 0) return

	status_msg_thr = 'filling gaps'
      CALL FILLGAPS                    !   Interpolation
	if (flag_thr .LT. 0) return

      if (nsndtype .NE. 3)  then
		status_msg_thr = 'spline smoothing'
		CALL SMOOTHSPL                   !   Filtering
	else
		status_msg_thr = 'lowpass smoothing'
		call smooth
	endif

	qc_data_thr = snddat

	if (flag_thr .LT. 0) return

	status_msg_thr = 'calculating heights'
      CALL HEIGHTS                     !   Geopotential heights
	if (flag_thr .LT. 0) return

	status_msg_thr = 'display sounding'
      CALL LVLINT(10.,.TRUE.,.FALSE.)  !   Display sounding

c	
c	snddat now has QC'd data with calculated heights
c

	qc_data_thr = snddat

c      CALL SNDNG(.FALSE.,.TRUE.)       !   Write to log file
c	status_msg_thr = 'coding TEMPDROP'
c	CALL TEMPDROP(flname_thr)               !   Write tempdrop file
C
      status_msg_thr = "normal termination"
	  flag_thr = 3						! succesfull termination
      return
      END 
C 
C
C     --------------------------------------------------------- 
      SUBROUTINE BIAS(PBIAS)
C 
C     Routine looks for possible biases in p data.
C     --------------------------------------------------------- 
C 
      PARAMETER (MXRC = 9000, NINVAR = 12)
C 
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
C 
      CHARACTER*4 VERSION
      LOGICAL AUTOTEMP
      DATA BAD/-999./
C
      IF (NSNDTYPE.EQ.1) THEN
         T = 10.0
         PDIFFBAR = 9.0
         TOLERANCE = 3.0
         ELSE
         T = 3.0
         PDIFFBAR = 3.6
         TOLERANCE = 1.0
         ENDIF
C
      CALL POLATE(IXEND,SNDDAT(1,1),SNDDAT(1,2),T,P,M,BAD)
      PDIFF = P-PRL-PDIFFBAR
      IF (P.EQ.BAD .OR. PRL.EQ.BAD) RETURN
      IF (NSITE.EQ.2 .AND. ABS(PDIFF).LT.TOLERANCE) RETURN
C
C
      WRITE(LUT,'(" POSSIBLE PRESSURE BIAS OF ",F5.1" MB.")')
     *      PDIFF
      PBIAS = PDIFF
      IF (ABS(PDIFF).GE.TOLERANCE) WRITE(LUT,'(A)') CHAR(7)
C
      RETURN
      END
C
C
C
C     --------------------------------------------------------- 
      SUBROUTINE AUTOQC(PBIAS)
C 
C     Routine performs automatic quality control on sonde data.
C     AUTOQC flags bad points, but does NOT interpolate.
C     --------------------------------------------------------- 
C 
      PARAMETER (MXRC = 9000, NOPT = 16, NINVAR = 12)
      DIMENSION IBAD(MXRC), PRGAPMAX(2), PR_AMBIENT(3,3)
      CHARACTER*1 OPTD(5), ANS
      CHARACTER*4 STNS, VERSION
      CHARACTER*20 INPUT
      LOGICAL LAST, NOGOODDATA, NEEDOPT, AUTOTEMP, NOQC
C 
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
      COMMON /OPTIONS/ NEEDOPT(NOPT)
C 
      DATA BAD/-999./
      DATA PRGAPMAX/50.,50./
      DATA OPTD/'P','T','H','W','Q'/
      DATA PR_AMBIENT/10.,3.,10.,                ! PR ADJ DEPTH: ODW,GPS,GLASS
     *                15.,20.,15.,               ! TE ADJ DEPTH: ODW,GPS,GLASS
     *                20.,50.,20./               ! RH ADJ DEPTH: ODW,GPS,GLASS
C 
C 
C     Select parameter to be QC'd
C     ---------------------------
      IOPT = 0
100   IOPT = IOPT+1
      IF (AUTOTEMP) THEN
         ANS = ' '
         GOTO 110
         ENDIF
105   WRITE(LUT,'(/," Select parameter for quality control:")') 
      WRITE(LUT,'(  " -------------------------------------")') 
      WRITE(LUT,'(  " (P) - Pressure")')
      WRITE(LUT,'(  " (T) - Temperature")') 
      WRITE(LUT,'(  " (H) - Humidity")')
      WRITE(LUT,'(  " (W) - Wind")')
      WRITE(LUT,'(  " (Q) - Quit this option")')
      WRITE(LUT,'(  " ----------------------")') 
      WRITE(LUT,'(/," Enter selection [",A1,"]: ",$)') OPTD(IOPT)
      READ(*,'(A)') ANS
110   IF (ANS.EQ.' ') ANS = OPTD(IOPT)
C 
      CALL INDEX(ANS,IVX,IFX) 
      IF (IVX.EQ.-1) GOTO 105
      IF (IVX.EQ.0) THEN
         RETURN                
         ENDIF
C
C
C     -----------------------------------------
C     Set QC parameters depending on variable.
C     IVX = SNDDAT index number of the variable
C     IFX = ISNDFLG index of the variable
C     -----------------------------------------
      NOQC = .FALSE.
      IF (IVX.EQ.2) GOTO 120
      IF (IVX.EQ.3) GOTO 130
      IF (IVX.EQ.4) GOTO 140
      IF (IVX.EQ.6) GOTO 160
      GOTO 105
C
C     Pressure
C     ------------------------------------
120   ADJ_DEPTH = PR_AMBIENT(NSNDTYPE,IFX)
      if (nsndtype .GE. 3) then
		nhardflaga = 0
	else
		IF (PRL.GT.0. .AND. PRL.LT.350.) ADJ_DEPTH = ADJ_DEPTH*4.0
		NHARDFLAG = NELM(ADJ_DEPTH,FALLRATE,DATARATE)
		NHARDFLAGA = NAMBIENT(NHARDFLAG,PRL,SNDDAT(1,IVX))+1
	endif

      VALMX = 1200.
      VALMN = 10.
      TRDDEF = FALLRATE/DATARATE       
      TRDLL = TRDDEF*0.75              
      TRDUL = TRDDEF*1.25
      DEVB = 1.5
      DEVL = 1.0
      XSD = 4.5 
      CUTOFF = 60.0 
      IF (DATARATE.GE.0.4) CUTOFF = 30.
      IF (nsndtype.EQ.3  ) CUTOFF = 100.
      INTRVL = NINT(1./DATARATE) 
      DEVX = 1.5
      GOTO 200
C
C     Temperature
C     ------------------------------------
130   ADJ_DEPTH = PR_AMBIENT(NSNDTYPE,IFX)
      IF (PRL.GT.0. .AND. PRL.LT.350.) ADJ_DEPTH = ADJ_DEPTH*4.0
      if (nsndtype .GE. 3) then
		nhardflaga = 0
	else
		NHARDFLAG = NELM(ADJ_DEPTH,FALLRATE,DATARATE)
		NHARDFLAGA = NAMBIENT(NHARDFLAG,TEL,SNDDAT(1,IVX))
	endif
      TC = TEMPTC(HTL,BAD)
      IF (TC.NE.BAD) NHARDFLAGA = INT(5.0*TC*DATARATE + 1.)
      VALMX = 50.
      VALMN = -100.
      TRDDEF = (FALLRATE/DATARATE)*0.08            
      TRDLL = 0.
      TRDUL = TRDDEF*2.0
      DEVB = 0.5
      DEVL = 2.0
      XSD = 5.0 
      CUTOFF = 60.0 
      IF (DATARATE.GE.0.4) CUTOFF = 30.
      IF (nsndtype.EQ.3  ) CUTOFF = 100.
      INTRVL = NINT(1./DATARATE)
      DEVX = 1.0
      GOTO 200
C
C     Relative humidity
C     ------------------------------------
140   ADJ_DEPTH = PR_AMBIENT(NSNDTYPE,IFX)
      IF (PRL.GT.0. .AND. PRL.LT.350.) ADJ_DEPTH = ADJ_DEPTH*4.0
      if (nsndtype .GE. 3) then
		nhardflaga = 0
	else
		NHARDFLAG = NELM(ADJ_DEPTH,FALLRATE,DATARATE)
		NHARDFLAGA = NAMBIENTRH(NHARDFLAG,SNDDAT(1,IVX))
	endif

      TC = RHTC(HTL,BAD)
      IF (TC.NE.BAD) THEN
         NHARDFLAGA = INT(4.0*TC*DATARATE + 1.)
         IF (TEL.NE.BAD.AND.TEL.GT.-30.) NHARDFLAGA = NHARDFLAGA/2
      ENDIF

      VALMX = 100.
      VALMN = 0.
      TRDDEF = 0.0                     
      TRDLL = 0.
      TRDUL = 0.
      DEVB = 20.0
      IF (DATARATE.GE.0.4) DEVB = 3.0
      DEVL = 30.0
      IF (DATARATE.GE.0.4) DEVL = 20.0
      XSD = 10.0 
      CUTOFF = 40.0 
      IF (DATARATE.GE.0.4) CUTOFF = 20.
      IF (nsndtype.EQ. 3)  CUTOFF = 100.
      INTRVL = NINT(1./DATARATE)
      DEVX = 20.0
      GOTO 200
C
C     Winds
C     -----
160   NHARDFLAG = 0
      NHARDFLAGA = 0
      VALMX = 360.
      VALMN = 0.
      TRDDEF = 0.
      TRDLL = 0.
      TRDUL = 0.
      DEVB = 1.5
      DEVL = 1.0
      DEVV = 3.0
      DEVH = 2.0
      XSD = 5.0 
      CUTOFF = 60.0 
      IF (DATARATE.GE.0.4) CUTOFF = 30.
      IF (nsndtype .GE. 3) CUTOFF = 100.
      INTRVL = NINT(1./DATARATE) 
      DEVX = 3.0
      GOTO 200
C
C
C
C     ---------------------------
C     Quality control begins here
C     ---------------------------
C
C     Reset all flags, check for values out of limits
C     -----------------------------------------------
200   DO 210 I = 1,IXEND
         ISNDFLG(I,IFX) = 0
         IF (SNDDAT(I,IVX).LT.VALMN .OR. SNDDAT(I,IVX).GT.VALMX) THEN
            SNDDAT(I,IVX) = BAD
            IF (IVX.EQ.6) SNDDAT(I,7) = BAD
            ISNDFLG(I,IFX) = 1
            ENDIF
210      CONTINUE
C
C
C     Assign data to temporary array X
C     --------------------------------
      NOGOODDATA = .TRUE.
      DO 220 I = 1, IXEND 
         IF (IVX .LE. 4) X(I) = SNDDAT(I,IVX) 
         IF (IVX .EQ. 6) THEN
            X(I) = UCMP(SNDDAT(I,6),SNDDAT(I,7)) 
            Y(I) = VCMP(SNDDAT(I,6),SNDDAT(I,7)) 
            ENDIF
         IF (X(I).NE.BAD) NOGOODDATA = .FALSE.
220      CONTINUE
C
C
C     Compute trend in the data, display parameters
C     ---------------------------------------------
      NSKP = NELM(50.,FALLRATE,DATARATE)
      TRD = (X(IXEND-NSKP)-X(NSKP)) / FLOAT(IXEND-NSKP-NSKP)
      TREND = TRD
      IF (TRD.LT.TRDLL .OR. TRD.GT.TRDUL) TREND = TRDDEF
      CALL QCDISP(IFX,XSD,TREND,CUTOFF,DEVX,DEVB)
C
C
C     Check if all the data are bad
C     -----------------------------
      IF (NOGOODDATA) THEN
         WRITE(LUT,'(/,
     *   " *** THERE ARE NO VALID DATA FOR THIS PARAMETER ***")')
         NOQC = .TRUE.
         DO 225 I = 1, IXEND
            ISNDFLG(I,IFX) = 2
225         CONTINUE
         IF (IVX.EQ.2 .AND. .NOT.AUTOTEMP) THEN
            WRITE(LUT,'(
     *      " *** PRESSURES CAN BE ESTIMATED IN OPTION 6 ***")')
            NEEDOPT(6) = .TRUE.
            ENDIF
         ENDIF
      IF (NOQC) GOTO 800
C
C
C     Check for bias against launch value
C     -----------------------------------
      IF (IVX.EQ.2) CALL BIAS(PBIAS)
C
C
C     Run data lock check (must be before detrend)
C     --------------------------------------------
      CALL VALUELOCK(X,IXEND,DEVL,NBAD,IBAD,IFX) 
      IF (NBAD.GT.0) THEN
         WRITE(LUT,'(" Data locked at index =         ",$)')
         WRITE(LUT,'(8I5,/,100(32X,8I5,/))') (IBAD(L),L=1,NBAD)
         ENDIF
      CALL FLAGSET(IFX,NBAD,IBAD)
      IF (IVX.EQ.6) THEN
         CALL VALUELOCK(Y,IXEND,DEVL,NBAD,IBAD,IFX) 
         IF (NBAD.GT.0) THEN
            WRITE(LUT,'(" Data locked at index =         ",$)')
            WRITE(LUT,'(8I5,/,100(32X,8I5,/))') (IBAD(L),L=1,NBAD)
            ENDIF
         CALL FLAGSET(IFX,NBAD,IBAD)
         ENDIF
C
C
C     Set hard-wire flags at top of sounding (PTH only)
C     -------------------------------------------------
      IF (IVX.LE.4) THEN
         WRITE(LUT,'(" Ambient equilibration found at IX = ",
     *   i3)') NHARDFLAGA
         DO 230 I = 1,NHARDFLAGA
            ISNDFLG(I,IFX) = 1
            X(I) = BAD
230         CONTINUE
         ENDIF
C
C
C     Wind quality checks (WQ for ODW, VVEL for GPS)
C     ----------------------------------------------
C
C     ODWS
C     --------------------------------------
      IF (IVX.EQ.6 .AND. NSNDTYPE.EQ.1) THEN
         LWTU = 7
235      WRITE(LUT,'(//," Verify lowest good WQ = ",I1," : ",$)') 
     *   LWTU
         READ(*,'(A)') INPUT 
         IF (INPUT.NE.' ') READ (INPUT,*,ERR=235) LWTU
         DO 240 I = 1,IXEND
            IF (SNDDAT(I,8).LT.LWTU) THEN
               ISNDFLG(I,4) = 1 
               X(I) = BAD
               ENDIF
240         CONTINUE 
         ENDIF
C
C     GPS dropsondes
C     --------------------------------------
      IF (IVX.EQ.6 .AND. NSNDTYPE.EQ.2) THEN
         CALL FFCHECK(NV,NF)
         FF = 0.
         IF (NV.GT.0) FF = FLOAT(NF)/FLOAT(NV)
C
         IF (FF.GT.0.5) THEN
            WRITE(LUT,'(/," Fast fall sonde...all winds tossed.")')
            DO 245 I = 1,IXEND
               ISNDFLG(I,4) = 2
245            CONTINUE
            GOTO 800
            ENDIF                  
C
         CALL VVCHECK(DEVV,DEVH,IXEND,BAD,NBAD,IBAD,AVGW)
         WRITE(LUT,'(" Mean vertical wind component = ",
     *   f7.2," m/s.")') AVGW
         IF (NBAD.GT.0 .AND. IVX.EQ.6) THEN
            WRITE(LUT,'(" VVEL check flags at index =    ",$)')
            WRITE(LUT,'(8I5,/,100(32X,8I5,/))') (IBAD(L),L=1,NBAD)
            ENDIF
         CALL FLAGSET(4,NBAD,IBAD)
         DO 250 I = 1,NBAD
            X(IBAD(I)) = -999.
250         Y(IBAD(I)) = -999.
         ENDIF
C
C
C       Now detrend the data
C       --------------------
        DO 260 I = 1, IXEND 
           IF (X(I).GT.BAD) X(I) = X(I)-TREND*FLOAT(I-1) 
260        CONTINUE
C
C
C       Run buddyckeck
C       --------------
        CALL BUDDYCHK(X,IXEND,DEVB,NBAD,IBAD) 
        IF (NBAD.GT.0) THEN
           WRITE(LUT,'(" Buddy check flags at index =   ",$)')
           WRITE(LUT,'(8I5,/,100(32X,8I5,/))') (IBAD(L),L=1,NBAD)
           ENDIF
        CALL FLAGSET(IFX,NBAD,IBAD)
        IF (IVX.NE.6) GOTO 280
        CALL BUDDYCHK(Y,IXEND,DEVB,NBAD,IBAD) 
        IF (NBAD.GT.0) THEN
           WRITE(LUT,'(" Buddy check flags at index =   ",$)')
           WRITE(LUT,'(8I5,/,100(32X,8I5,/))') (IBAD(L),L=1,NBAD)
           ENDIF
        CALL FLAGSET(IFX,NBAD,IBAD)
C
C
C       Run outlier check
C       -----------------
280     CALL AUTOFLAG(XSD,X,IXEND,NBAD,IBAD)
        IF (NBAD.GT.0) THEN
           WRITE(LUT,'(" Outlier check flags at index = ",$)')
           WRITE(LUT,'(8I5,/,100(32X,8I5,/))') (IBAD(L),L=1,NBAD)
           ENDIF
        CALL FLAGSET(IFX,NBAD,IBAD)
        IF (IVX.NE.6) GOTO 300
        CALL AUTOFLAG(XSD,Y,IXEND,NBAD,IBAD)
        IF (NBAD.GT.0) THEN
           WRITE(LUT,'(" Outlier check flags at index = ",$)')
           WRITE(LUT,'(8I5,/,100(32X,8I5,/))') (IBAD(L),L=1,NBAD)
           ENDIF
        CALL FLAGSET(IFX,NBAD,IBAD)
C
C
C       Run filter check
C       ----------------
300     CALL VSPLFLAG(SNDDAT(1,1),X,IXEND,CUTOFF,DEVX,NBAD,IBAD)
        IF (NBAD.GT.0) THEN
           WRITE(LUT,'(" Filter check flags at index =  ",$)')
           WRITE(LUT,'(8I5,/,100(32X,8I5,/))') (IBAD(L),L=1,NBAD)
           ENDIF
        CALL FLAGSET(IFX,NBAD,IBAD)
        IF (IVX.NE.6) GOTO 350
        CALL VSPLFLAG(SNDDAT(1,1),Y,IXEND,CUTOFF,DEVX,NBAD,IBAD)
        IF (NBAD.GT.0) THEN
           WRITE(LUT,'(" Filter check flags at index =  ",$)')
           WRITE(LUT,'(8I5,/,100(32X,8I5,/))') (IBAD(L),L=1,NBAD)
           ENDIF
        CALL FLAGSET(IFX,NBAD,IBAD)
C 
C
C     -----------------------------------------------------------
C     QC is finished. Check for runs of 'Y' flags longer than 
C     PRGAPMAX and change all the flags to 'X'.  (for T,H, winds)
C     -----------------------------------------------------------
C
350   IF (IVX.EQ.2) GOTO 800
      JYLIM = NELM(PRGAPMAX(NSITE),FALLRATE,DATARATE)
      LAST = .FALSE.
      DO 360 I = 1,IXEND
	 IF (ISNDFLG(I,IFX).EQ.1) THEN
            IF (.NOT. LAST) ISTRT = I
	    LAST = .TRUE.
	    ELSE
	    IF (LAST) THEN
	       ISTOP = I-1
	       IF (ISTOP-ISTRT .GE. JYLIM) THEN
	  	  DO 370 J = ISTRT,ISTOP
370		     ISNDFLG(J,IFX) = 2
   		  ENDIF
	       ENDIF
    	    LAST = .FALSE.
	    ENDIF
360	 CONTINUE
C
C
C     Go back to menu
C     ---------------
800   WRITE(LUT,'(" QC for this parameter is complete.")')
      WRITE(LUT,'(1x,71("-"),//)')
      GOTO 100
      END 
C 
C
C
C     -----------------------------------------------
      FUNCTION NAMBIENT(N,XL,X)
C     Finds index of first increasing value (ambient)
C     -----------------------------------------------
C
      DIMENSION X(N)
C
      NAMBIENT = N
      NLAST = 0
C
      DO 100 I = 1,N-1
         IF (X(I).EQ.-999.) GOTO 100
         IF (X(I).LT.XL) GOTO 100
         XLAST = X(I)
         NLAST = I
         GOTO 150
100      CONTINUE
C
150   IF (NLAST.EQ.0) GOTO 300
C
      DO 200 I = NLAST+1,N
         IF (X(I).EQ.-999.) GOTO 200
         IF (X(I).GT.XLAST) THEN
            NAMBIENT = I-1
            GOTO 300
            ENDIF
         XLAST = X(I)
200      CONTINUE
C
300   RETURN
      END
C
C
C
C     -----------------------------------------------
      FUNCTION NAMBIENTRH(N,X)
C     Finds index of first change in trend (ambient)
C     -----------------------------------------------
C
      DIMENSION X(N)
      LOGICAL INCR
C
      NAMBIENTRH = N
      NLAST = 0
C
      DO 100 I = 3,N-1
         IF (X(I).EQ.-999.) GOTO 100
         XLAST = X(I)
         NLAST = I
         GOTO 150
100      CONTINUE
C
150   IF (NLAST.EQ.0) GOTO 300
C
      INCR = .FALSE.
      DO 200 I = NLAST+1,N
         IF (X(I).EQ.-999. .OR. X(I).EQ.XLAST) GOTO 200
         IF (X(I).GT.XLAST) INCR = .TRUE.
         GOTO 225
200      CONTINUE
      GOTO 300
C
225   DO 250 I = NLAST+1,N
         IF (X(I).EQ.-999.) GOTO 250
         IF ((.NOT.INCR .AND. X(I).GT.XLAST) .OR.
     *        (INCR .AND. X(I).LT.XLAST)) THEN
            NAMBIENTRH = I-1
            GOTO 300
            ENDIF
         XLAST = X(I)
250      CONTINUE
C
300   RETURN
      END
C
C
C
C     -----------------------------------------------------------------
      SUBROUTINE QCDISP(IVAR,XSD,TREND,CUTOFF,DEVX,DEVB)
C     -----------------------------------------------------------------
C
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
C
      LOGICAL AUTOTEMP
      CHARACTER*4 VERSION
      CHARACTER*11 LABEL(4)
      CHARACTER*3 UNITS(4)
      DATA LABEL/'Pressure   ','Temperature','Humidity   ',
     *           'Wind (u,v) '/
      DATA UNITS/'mb','C','%','m/s'/
C
C
100   WRITE(LUT,'(//,1x,71("-"))')
      WRITE(LUT,'(" Quality control for parameter: ",A11)')
     *      LABEL(IVAR)
      WRITE(LUT,'(1x,71("-"))')
      WRITE(LUT,'(" Trend removed = ",F6.2," ",A,"/sample")')
     *      TREND,UNITS(IVAR)
      WRITE(LUT,'(" Buddy check:   Flagged when delta(",A2,") = ",
     *      F6.2," ",A)') LABEL(IVAR),DEVB,UNITS(IVAR)
      WRITE(LUT,'(" Outlier check: Flagged when deviation = ",
     *      F6.2," s.d.")') XSD
      WRITE(LUT,'(" Filter check:  Flagged when deviation = ",F6.2,
     *      " ",A)') DEVX,UNITS(IVAR)
      WRITE(LUT,'("                Filter wavelength      = ",F6.2,
     *      " s")') CUTOFF
C
      WRITE(LUT,'(1x,71("-"))')
C
      RETURN
      END
C
C
C
C     ----------------------------------
      SUBROUTINE FLAGSET(IVAR,NBAD,IBAD)
C     ----------------------------------
C
      PARAMETER (MXRC = 9000, NOPT = 16, NINVAR = 12) 
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      DIMENSION IBAD(NBAD)
C
      IVX = IVAR
      IF (IVX.GE.4) IVX = 4
      IF (NBAD.EQ.0) RETURN 
      DO 110 I = 1,NBAD 
110      ISNDFLG(IBAD(I),IVX) = 1 
      RETURN
      END
C
C
C
C     ----------------------------------------------------------
      SUBROUTINE FLTRFLAG(X,NDIMX,CUTOFF,INTRVL,DEVX,NBAD,IBAD) 
C 
C     Routine filters X and identifies points that deviate from 
C     filtered array. 
C 
C     X       =   Input array to be tested
C     NDIMX   =   Dimension of X
C     CUTOFF  =   Filter cutoff wavelength (s)
C     INTRVL  =   Time interval (s) between elements of X 
C     DEVX    =   Deviation threshhold for flagging 
C     NBAD    =   Number of elements of X identified as bad 
C     IBAD    =   Array containing flagged array element #'s
C     ----------------------------------------------------------
C 
      PARAMETER (MXRC = 9000)
C 
      DIMENSION X(NDIMX),Y(MXRC)
      DIMENSION IBAD(NDIMX) 
      DIMENSION IBGN(100), IEND(100)
C 
C 
      NBAD = 0
      NTERMS = 25 
      WCUT = 2.0*FLOAT(INTRVL)/CUTOFF 
C 
      DO 100 I = 1,NDIMX
         Y(I) = X(I)
100      CONTINUE 
C 
      CALL GAP(X,NDIMX,IBGN,IEND,NGAPS,-999.)     
      DO 150 K = 1,NGAPS+1
         NTG = IEND(K)-IBGN(K)+1
         IF (NTG*INTRVL .LT. INT(CUTOFF)) GOTO 150
         CALL LOPASS(X(IBGN(K)),Y(IBGN(K)),NTG,WCUT,NTERMS) 
150      CONTINUE 
C 
      DO 500 I = 3,NDIMX-2
         IF (X(I).LE.-999.) GOTO 500
         IF (ABS(X(I)-Y(I)).GT.DEVX) THEN 
              NBAD = NBAD+1 
              IBAD(NBAD) = I
              NBADX = NBADX+1 
              ENDIF 
500      CONTINUE 
C 
      RETURN
      END 
C 
C 
C 
C     ----------------------------------------------------------
      SUBROUTINE VSPLFLAG(XT,X,NDIMX,YDCWL,DEVX,NBAD,IBAD) 
C 
C     Routine filters X and identifies points that deviate from 
C     filtered array. 
C 
C     XT      =   Array of times
C     X       =   Data array to be tested
C     NDIMX   =   Dimension of X,XT
C     YDCWL   =   Filter cutoff wavelength (s)
C     DEVX    =   Deviation threshhold for flagging 
C     NBAD    =   Number of elements of X identified as bad 
C     IBAD    =   Array containing flagged array element #'s
C     ----------------------------------------------------------
C 
      PARAMETER (MXRC = 9000)
C 
      DIMENSION XT(NDIMX),X(NDIMX),XW(MXRC)
      DIMENSION IBAD(NDIMX) 
      LOGICAL ECHO
C
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
C
      DATA ECHO/.TRUE./
C 
C 
      NBAD = 0
      KYBC = 2
      YBCWL = YDCWL
      KDAT = 1
      DO 100 I = 1,NDIMX
         XW(I) = 1.
100      CONTINUE
C 
      CALL VICSETUP(XT,NDIMX,YDCWL,NX,YNB,YNT,2.0,IERR,ECHO)
      IF (IERR.NE.0) GOTO 900
      CALL VICSPL(XT,XW,X,NDIMX,MXRC,KDAT,YNB,YNT,NX,YDCWL,
     *            KYBC,KYBC,YBCWL,YBCWL,IERR)
      IF (IERR.NE.0) GOTO 900
      DO 500 I = 3,NDIMX-2
         IF (X(I).LE.-999.) GOTO 500
	 CALL SPOTVAL(XT(I),KDAT,FOUT,FOUTD)
         IF (ABS(X(I)-FOUT).GT.DEVX) THEN 
              NBAD = NBAD+1 
              IBAD(NBAD) = I
              NBADX = NBADX+1 
              ENDIF 
500      CONTINUE 
C 
      RETURN
C
900   WRITE(LUT,'(
     *  " *** WARNING: FILTERING FAILED...CHECK DATA ***")')
      RETURN
C
      END 
C 
C 
C 
C     ----------------------------------------------------------
      SUBROUTINE VVCHECK(DEVV,DEVH,NDIMX,BAD,NBAD,IBAD,AVGW) 
C 
C     Routine checks GPS sonde fall rate against theoretical and
C     hydrostatic fall rates.
C 
C     DEVV    =   Deviation threshhold for GPS difference
C     DEVH    =   Threshold for hydrostatic difference
C     NBAD    =   Number of elements of X identified as bad 
C     IBAD    =   Array containing flagged array element #'s
C     ----------------------------------------------------------
C 
      PARAMETER (MXRC = 9000, NINVAR = 12) 
C 
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
C 
      DIMENSION IBAD(NDIMX) 
      CHARACTER*4 VERSION
      LOGICAL AUTOTEMP,DEBUG
C 
      DATA DEBUG/.FALSE./
C
C 
      NBAD = 0
      NT = 0
      SUM = 0.0
      DO 100 I = 1,IXEND
         EXPFALL = GSNDFALL(SNDDAT(I,5),SNDDAT(I,3),BAD)
         GPSFALL = SNDDAT(I,10)
         IF (GPSFALL.EQ.BAD .OR. EXPFALL.EQ.BAD) GOTO 100
         DIFF = GPSFALL-EXPFALL
         HYDFALL = VVHYD(I)
         DIFFH = BAD
         IF (HYDFALL.NE.BAD) DIFFH = GPSFALL-HYDFALL
         NT = NT+1
         SUM = SUM+DIFF
         IF (DEBUG) WRITE(LUPR,*) I,EXPFALL,GPSFALL,HYDFALL,DIFF,DIFFH
         IF (ABS(DIFF).GT.DEVV .AND. ABS(DIFFH).GT.DEVH) THEN
            NBAD = NBAD+1
            IBAD(NBAD) = I
            ENDIF
100      CONTINUE
C 
C 
      IF (NT.EQ.0) THEN
         AVGW = BAD
         ELSE
         AVGW = SUM/FLOAT(NT)
         ENDIF
C
      RETURN
      END 
C 
C 
C 
C     -----------------------------------------------
      FUNCTION VVHYD(N)
C
C     Calculates hydrostatic vertical velocity
C     -----------------------------------------------
C
      PARAMETER (MXRC = 9000, NINVAR = 12) 
C 
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
C
      CHARACTER*4 VERSION
      LOGICAL AUTOTEMP
      DATA BAD/-999./
C
C
      VVHYD = BAD
      N1 = N-1
      N2 = N+1
      IF (N1.LT.1) N1=1
      IF (N2.GT.IXEND) N2=IXEND
      IF (SNDDAT(N1,9).EQ.BAD .OR. SNDDAT(N2,9).EQ.BAD) GOTO 900
      VVHYD = (SNDDAT(N2,9)-SNDDAT(N1,9))/(SNDDAT(N2,1)-SNDDAT(N1,1))
C
900   RETURN
      END
C
C
C
C     ----------------------------------------------------------
      SUBROUTINE FFCHECK(NV,NF) 
C 
C     Routine checks for unusually fast fall rate.
C 
C     NV = # of valid checks obtained.
C     NF = # of valid checks that were unusually fast.
C     ----------------------------------------------------------
C 
      PARAMETER (MXRC = 9000, NINVAR = 12) 
C 
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
C 
      CHARACTER*4 VERSION
      LOGICAL AUTOTEMP
C
      DATA BAD/-999./
C 
C 
      NV = 0
      NF = 0
      IF (NSNDTYPE .NE. 2) RETURN
C
      DO 100 I = 1,IXEND
         EXPFALL = GSNDFALL(SNDDAT(I,2),SNDDAT(I,3),BAD)
         HYDFALL = VVHYD(I)
         IF (HYDFALL.EQ.BAD .OR. EXPFALL.EQ.BAD) GOTO 100
C
         RATIO = HYDFALL/EXPFALL
         IF (RATIO.LT.0 .OR. RATIO.GT.3.) GOTO 100
         NV = NV+1
         IF (RATIO.GT.1.25) NF = NF+1
100      CONTINUE
C 
C 
      RETURN
      END 
C 
C 
C 
C     -----------------------------------------------
      SUBROUTINE DYNTCORR(IERR)
C
C     Calculates dynamic temperature correction
C     -----------------------------------------------
C
      PARAMETER (MXRC = 9000, NINVAR = 12) 
C 
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
C
      CHARACTER*4 VERSION
      LOGICAL AUTOTEMP,DEBUG
      DATA BAD/-999./
      DATA AVGDEPTH/2.0/
      DATA DEBUG/.FALSE./
C
C
C     Determine where to start the adjustment.
C     We don't want to adjust during initial equilibration
C     because it doesn't work well, even if NLPS = 1.
C     This will also be the earliest lapse rate used.
C     ----------------------------------------------------
      IERR = 0
      NSTART = NELM(30.,FALLRATE,DATARATE)
      NSTART = NAMBIENT(NSTART,TEL,SNDDAT(1,3))
C
C
C     Calculate lapse rate over previous AVGDEPTH mb
C     ----------------------------------------------
      NLPS = NELM(AVGDEPTH,FALLRATE,DATARATE)
      IF (NLPS.LT.1) NLPS = 1
C
C
      TLPLAST = 0.
      DO 100 L = NSTART+1,IXEND
C
C        The default lapse rate DT/Dt is determined from mean 
C        tropical dT/dZ and the theoretical sonde fall rate.
C        ------------------------------------------------------
         DTDZ = TROPLAPSZ(SNDDAT(L,2),BAD)    
         EXPFALL = GSNDFALL(SNDDAT(L,5),SNDDAT(L,3),BAD)
         IF (DTDZ.NE.BAD .AND. EXPFALL.NE.BAD) THEN
            TLAPS = DTDZ * EXPFALL
            ELSE
            TLAPS = 0.0  
            ENDIF
C
C
C        Actually, let's try using the last lapse rate as default
C        --------------------------------------------------------
         TLAPS = TLPLAST
C
C
C        Calculate actual temperature lapse rate (deg/s)
C        If value seems unreasonalble, use default.
C        -----------------------------------------------
         N1 = L-NLPS
         N2 = L
         IF (N1.LT.NSTART) N1=NSTART
         IF (N2.GT.IXEND) N2=IXEND
C
         IF (N1.EQ.N2) GOTO 150
         IF (SNDDAT(N1,3).EQ.BAD .OR. SNDDAT(N2,3).EQ.BAD) GOTO 150
C
         TLP = (SNDDAT(N2,3)-SNDDAT(N1,3))/(SNDDAT(N2,1)-SNDDAT(N1,1))
         IF (TLP.GT.0.25 .OR. TLP.LT.-0.50) GOTO 150
C
C        Value is ok
C        -----------
         TLAPS = TLP
         TLPLAST = TLAPS
C
C
C        Calculate time constant, temperature correction
C        -----------------------------------------------
150      TC = TEMPTC(SNDDAT(L,9),BAD)
         IF (TC.EQ.BAD .OR. TLAPS.EQ.BAD) THEN
            TCORR = BAD
            ELSE
            TCORR = TC * TLAPS
            ENDIF
         IF (DEBUG) WRITE(LUPR,*) L,SNDDAT(L,3),TC,TLAPS,TCORR
         X(L) = TCORR
100      CONTINUE
C
C
C      Go back and fix temperatures.  If correction is bad, use
C      previous correction.
C      --------------------------------------------------------
       TCORR = 0.
       NBAD = 0
       NTOT = 0
       DO 200 L = NSTART+1,IXEND
          IF (SNDDAT(L,3).EQ.BAD) GOTO 200
          IF (X(L).NE.BAD) THEN
             TCORR = X(L)
             ELSE
             NBAD = NBAD+1
             ENDIF
          NTOT = NTOT+1
          SNDDAT(L,3) = SNDDAT(L,3)+TCORR
200       CONTINUE
C
C
      IF (NBAD.EQ.0) 
     *   WRITE(LUT,'(/," Dynamic temperature correction applied.",/)')
      IF (NBAD.GT.0 .AND. NBAD.LT.NTOT) 
     *   WRITE(LUT,'(/," Calculation failed on ",I4," of ",
     *   I4," data points (default correction used).",/)') NBAD,NTOT
      IF (NBAD.EQ.NTOT) THEN
         IERR = 1
         WRITE(LUT,'(/,
     *   " *** WARNING: Dynamic temperature correction failed *** ",/)')
         ENDIF
      RETURN
      END
C
C
C
C     ----------------------------------------------------------
      SUBROUTINE BUDDYCHK(X,NDIMX,DEVX,NBAD,IBAD) 
C 
C     Routine checks array X and identifies points that deviate from 
C     the two neighboring values (in same direction) by at least DEVX.  
C 
C     X       =   Input array to be tested
C     NDIMX   =   Dimension of X
C     DEVX    =   Deviation threshhold for flagging 
C     NBAD    =   Number of elements of X identified as bad 
C     IBAD    =   Array containing flagged array element #'s
C     ----------------------------------------------------------
C 
      PARAMETER (MXRC = 9000)
C 
      DIMENSION X(NDIMX)
      DIMENSION IBAD(NDIMX) 
C 
C 
      NBAD = 0
C 
      DO 500 I = 2,NDIMX-1
         IF (X(I).LE.-999.) GOTO 500
C
         DO 550 L = I-1,I-5,-1
            IF (L.EQ.1) GOTO 555
            IF (X(L).GT.-999.) GOTO 555
550         CONTINUE
C
555      IB = L
         DO 560 L = I+1,I+5
            IF (L.EQ.NDIMX) GOTO 580
            IF (X(L).GT.-999.) GOTO 580
560         CONTINUE
C
580      IA = L
         IF (X(IB).LE.-999. .OR. X(IA).LE.-999.) GOTO 500
C
         DEV1 = X(I)-X(IB)
         DEV2 = X(I)-X(IA)
         IF (ABS(DEV1).LT.DEVX) GOTO 500
         IF (ABS(DEV2).LT.DEVX) GOTO 500
         IF ((DEV1*DEV2) .LT. 0.) GOTO 500
         NBAD = NBAD+1 
         IBAD(NBAD) = I
         NBADX = NBADX+1 
500      CONTINUE 
C 
      IF (NBAD.EQ.0) RETURN
      DO 510 I = 1,NBAD
	 X(IBAD(I)) = -999.
510      CONTINUE
C
      RETURN
      END 
C 
C 
C 
C     ----------------------------------------------------------
      SUBROUTINE VALUELOCK(X,NDIMX,DEVX,NBAD,IBAD,IVAR) 
C 
C     Routine checks array X and identifies regions locked on one
C     value. Flags region if it begins or ends with jump of DEVX.
C     Must work on original data, NOT DETRENDED data.
C 
C     X       =   Input array to be tested
C     NDIMX   =   Dimension of X
C     DEVX    =   Deviation threshhold for flagging 
C     NBAD    =   Number of elements of X identified as bad 
C     IBAD    =   Array containing flagged array element #'s
C     IVAR    =   Variable (1=P, 2=T, 3=H)
C     ----------------------------------------------------------
C 
      PARAMETER (MXRC = 9000)
C 
      DIMENSION X(NDIMX),DEV(NDIMX)
      DIMENSION IBAD(NDIMX) 
C 
C 
      NBAD = 0
C 
C     Find array of increments DEV
C     ----------------------------
      DO 500 I = 2,NDIMX
         DEV(I) = -999.
         IF (X(I).LE.-999.) GOTO 500
C
         DO 550 L = I-1,1,-1
            IF (X(L).LE.-999.) GOTO 550
            DEV(I) = X(I)-X(L)
            GOTO 500
550         CONTINUE
C
500      CONTINUE 
C 
C     Scan increments and find blocks of constant value
C     -------------------------------------------------
      ISKIP = 0
      DO 600 I = 3,NDIMX-1
         IF (I.LE.ISKIP) GOTO 600
         IF (DEV(I).EQ.-999.) GOTO 600
         IF (DEV(I).EQ.0.) THEN
            IA = -1
            IB = -1
            DO 650 L = I-1,1,-1
               IF (DEV(L).EQ.-999.) GOTO 650
               DEVB = DEV(L)
               IB = L
               GOTO 660
650            CONTINUE
C
660         DO 670 L = I+1,NDIMX
               IF (DEV(L).EQ.-999.) GOTO 670
               IF (DEV(L).EQ.0.) GOTO 670
               DEVA = DEV(L)
               IA = L
               GOTO 675
670            CONTINUE
C
C           We have now identified a block of constant values that
C           runs from I = IB to IA-1. Only flag blocks longer than 2.
C           ------------------------------------------------------
675         IF (IB.EQ.-1) IB = 1
            IF (IA.EQ.-1) IA = NDIMX+1
            IF (IA-IB .LT. 3) GOTO 690
C
C           Allows humidity to lock on 100%
C           -------------------------------
            IF (IVAR.EQ.3 .AND. X(IB).EQ.100.) GOTO 690
C
            IF (ABS(DEVB).GT.DEVX .OR. ABS(DEVA).GT.DEVX) THEN
               DO 680 L = IB,IA-1
                  NBAD = NBAD+1 
                  IBAD(NBAD) = L
                  NBADX = NBADX+1 
680               CONTINUE
               ENDIF
C
C           Reset I to the end of this interval
C           -----------------------------------
690         ISKIP = IA-1
            ENDIF
C
C
600      CONTINUE
C         
C
      IF (NBAD.EQ.0) RETURN
      DO 710 I = 1,NBAD
	 X(IBAD(I)) = -999.
710      CONTINUE
C
      RETURN
      END 
C 
C 
C 
C     ----------------------------------------------------------
      SUBROUTINE AUTOFLAG(XSD,X,NDIMX,NBAD,IBAD)
C 
C     Routine identifies points >= than XSD standard deviations 
C     away from the mean.  Iterates until no more bad pts are
C     found.
C 
C     X       =   Input array to be tested
C     NDIMX   =   Dimension of X
C     NBAD    =   Number of elements of X identified as bad 
C     IBAD    =   Array containing flagged array element #'s
C     ----------------------------------------------------------
C 
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
C
      DIMENSION X(NDIMX), IBAD(NDIMX) 
      LOGICAL DEBUG
      DATA DEBUG/.FALSE./
C 
      NBAD = 0
C 
C     Begin first iterration
C     ----------------------
C 
100   NT = 0
      SUMX = 0.0
      SUMX2 = 0.0 
      NBADX = 0 
C 
      DO 110 I = 1,NDIMX
         IF (X(I).GT.-999.) THEN
              NT = NT+1 
              SUMX = SUMX+X(I)
              SUMX2 = SUMX2+X(I)*X(I) 
              ENDIF 
110      CONTINUE               
C 
      IF (NT.LT.2) RETURN 
      RNT = NT
      STDEV = ((RNT*SUMX2-SUMX*SUMX)/(RNT*(RNT-1.0)))**0.5
      AVG = SUMX/RNT
C 
      DO 120 I = 1,NDIMX
         IF (X(I).LE.-999.) GOTO 120
         IF (ABS(X(I)-AVG).GT.STDEV*XSD) THEN 
              IF (DEBUG) 
     *        WRITE(LUPR,*) I,X(I),ABS(X(I)-AVG),STDEV*XSD,STDEV
              NBAD = NBAD+1 
              IBAD(NBAD) = I
              NBADX = NBADX+1 
              X(I) = -999.         ! To be ignored on next pass.
              ENDIF 
120      CONTINUE 
C 
      IF (NBADX.GT.0) GOTO 100     ! Found bad pts, go back again.
C 
      RETURN
      END 
C 
C 
C 
C     ----------------------------------------------------- 
      SUBROUTINE STATUS 
C 
C     Allows user to manually set status of data flags, and 
C     to select the subsequent use of estimated pressure
C     instead of measured pressure. 
C     ----------------------------------------------------- 
C 
      PARAMETER (MXRC = 9000, NOPT = 16, NINVAR = 12) 
C 
      CHARACTER*1 ANS, OPTD(3)
      CHARACTER*4 STNS
      CHARACTER*4 VERSION
      LOGICAL NEEDOPT, WARN, AUTOTEMP, RHTOGGLE, SKIP, SELECTED
C 
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
      COMMON /OPTIONS/ NEEDOPT(NOPT)
C
      DATA OPTD/'T','Q','I'/
C 
C 
      RHTOGGLE = .FALSE.
      WARN = .FALSE.
      NDISP = NELM(25.,FALLRATE,DATARATE)
C
C
C     Select parameter to be flagged
C     ------------------------------
200   WRITE(LUT,'(/," Flag which parameter:")') 
      WRITE(LUT,'(  " ----------------------")') 
      WRITE(LUT,'(  " (P) - Pressure")')
      WRITE(LUT,'(  " (T) - Temperature")') 
      WRITE(LUT,'(  " (H) - Humidity")')
      WRITE(LUT,'(  " (W) - Wind")')
      WRITE(LUT,'(  " (Q) - Quit this option")')
      WRITE(LUT,'(  " ----------------------")') 
      WRITE(LUT,'(/," Enter selection [",A1,"]: ",$)') OPTD(2)
      READ(*,'(A)') ANS
      IF (ANS.EQ.' ') ANS = OPTD(2)
C 
      CALL INDEX(ANS,IXDAT,IXFLG) 
      IF (IXDAT.EQ.-1) GOTO 200    ! Invalid selection
      IF (IXDAT.EQ.0) THEN
         RETURN                ! Quit 
         ENDIF
C
C
C     Specify the interval for flagging
C     ---------------------------------
      NOPTD = 1
      IXSTRT = -1
      IXSTOP = -1
300   IF (IXSTOP.GT.0 .AND. IXSTRT.GT.0) GOTO 350
C
      WRITE(LUT,'(/," Set the interval for flagging:")') 
      WRITE(LUT,'(  " (Current end of drop at index = ",i4.4,")")') 
     *      IXEND
      WRITE(LUT,'(  " ------------------------------------------")') 
      WRITE(LUT,'(  " (A) - Flag ALL the data points")')
      IF (IXDAT.EQ.4) THEN
         WRITE(LUT,'(  " (1) - Flag all RH from SENSOR #1")')
         WRITE(LUT,'(  " (2) - Flag all RH from SENSOR #2")')
         ENDIF
      WRITE(LUT,'(  " (I) - Specify interval by INDEX number")')
      WRITE(LUT,'(  " (P) - Specify interval by PRESSURE")')
      WRITE(LUT,'(  " (T) - Specify interval by TIME")')
      WRITE(LUT,'(  " (L) - LIST the sounding")')
      WRITE(LUT,'(  " (Q) - CANCEL flagging, go to main menu")')
      WRITE(LUT,'(  " ------------------------------------------")') 
      WRITE(LUT,'(/," Enter selection [",A1,"]: ",$)') OPTD(NOPTD)
      READ(*,'(A)') ANS
      IF (ANS.EQ.' ') ANS = OPTD(NOPTD)
C
C
C     Check response
C     --------------
      IF (ANS.EQ.'Q' .OR. ANS.EQ.'q') RETURN
C
      IF (ANS.EQ.'I' .OR. ANS.EQ.'i') THEN
310      WRITE(LUT,'(/" Enter starting, ending index numbers: ",$)')
         READ(*,*,ERR=310) IXSTRT,IXSTOP
         IF (IXSTRT.LE.0 .OR. IXSTOP.GT.IXEND .OR.
     *       IXSTOP.LT.IXSTRT) THEN
            WRITE(LUT,'(" *** INVALID INTERVAL ***")') 
            GOTO 310
            ENDIF
         GOTO 350
         ENDIF
C
      IF (ANS.EQ.'T' .OR. ANS.EQ.'t') THEN
320      WRITE(LUT,'(/" Enter starting, ending times (s): ",$)')
         READ(*,*,ERR=320) TBGN,TEND
         IF (TEND.LT.TBGN) THEN
            WRITE(LUT,'(" *** INVALID INTERVAL ***")') 
            GOTO 320
            ENDIF
         IF (TBGN.LE.0.) TBGN = -2.0
         IXSTRT = IXTIME(TBGN)
         IXSTOP = IXTIME(TEND)
         WRITE(LUT,'(" Data will be flagged from index ",I4.4," to ",
     *               I4.4,".")') IXSTRT,IXSTOP
         GOTO 350
         ENDIF
C
      IF (ANS.EQ.'P' .OR. ANS.EQ.'p') THEN
330      WRITE(LUT,'(/" Enter starting, ending pressures (mb): ",$)')
         READ(*,*,ERR=330) PBGN,PEND
         IF (PEND.LT.PBGN) THEN
            WRITE(LUT,'(" *** INVALID INTERVAL ***")') 
            GOTO 330
            ENDIF
         IXSTRT = IXPRESS(PBGN)
         IXSTOP = IXPRESS(PEND)
         IF (PBGN.LT.PRL) IXSTRT = 1
         WRITE(LUT,'(" Data will be flagged from index ",I4.4," to ",
     *               I4.4,".")') IXSTRT,IXSTOP
         GOTO 350
         ENDIF
C
      IF (IXDAT.EQ.4 .AND. (ANS.EQ.'1' .OR. ANS.EQ.'2')) THEN
         RHTOGGLE = .TRUE.
         READ(ANS,*) IRHSF
         GOTO 350
         ENDIF         
C
      IF (ANS.EQ.'L' .OR. ANS.EQ.'l') THEN
         CALL SHOWSNDNG
         WRITE(LUT,'(" ")')
         NOPTD = 3
         GOTO 300
         ENDIF
C
      IF (ANS.EQ.'A' .OR. ANS.EQ.'a') THEN
         IXSTRT = 1
         IXSTOP = IXEND
         GOTO 350
         ENDIF
C
      GOTO 300
C
C 
C     Specify the status flag value
C     -----------------------------
350   WRITE(LUT,'(/," Set data status flag to:")') 
      WRITE(LUT,'(  " -------------------------------------")')
      WRITE(LUT,'(  " (R)eplace  (with interpolated value")')
      WRITE(LUT,'(  " (M)issing  (do NOT interpolate)")')
      WRITE(LUT,'(  " -------------------------------------")')
      WRITE(LUT,'(  " (K)eep     (data are good)")')
      WRITE(LUT,'(  " (D)oubtful (but keep data anyway)")')
      IF (NSITE.EQ.1) THEN
          WRITE(LUT,'(  " (S)ubjective edit (but keep)")')
          ENDIF
      WRITE(LUT,'(  " -------------------------------------")')
      WRITE(LUT,'(  " (Q)uit     (Cancel data flagging)")')
      WRITE(LUT,'(  " -------------------------------------")')
      WRITE(LUT,'(/," Enter selection: ",$)')
      READ(*,'(A)') ANS 
      CALL CXFLG(ANS,IVFLG) 
      IF (IVFLG.EQ.-2) RETURN
      IF (IVFLG.EQ.-1) GOTO 350    ! Invalid selection
      IF (IVFLG.LT.4) NEEDOPT(10) = .TRUE.

C 
C
C     Reset flags...do not reset bad data to "K"
C     ------------------------------------------
      DO 400 I = 1, IXEND 
         SKIP = .FALSE.
         IF (RHTOGGLE) THEN
            CALL RHTOGL(IRHSF,SNDDAT(I,1),SELECTED)
            SKIP = .NOT. SELECTED
            ELSE
            IF (I.LT.IXSTRT .OR. I.GT.IXSTOP) SKIP = .TRUE.
            ENDIF
         IF (SKIP) GOTO 400
C
         IF (SNDDAT(I,IXDAT).EQ.-999. .AND. 
     *       (IVFLG.EQ.0 .OR. IVFLG.GE.4)) THEN
            WARN = .TRUE.
            GOTO 400
            ENDIF
         ISNDFLG(I,IXFLG) = IVFLG
400      CONTINUE 
C
      IF (WARN) WRITE(LUT,'(/,
     * " Flags for missing data in the interval left unchanged."/)')
      IF (IVFLG.LE.3) NEEDOPT(12) = .TRUE.
      IF (IVFLG.NE.2 .AND. NSITE.EQ.2) NEEDOPT(11) = .TRUE.
C 
      GOTO 200
      END 
C 
C 
C 
C     -----------------------------------------------------
      SUBROUTINE RHTOGL(ISEN,T,SELECTED)
C
C     Determines whether the RH value at time T is from
C     RH sensor ISEN.
C     -----------------------------------------------------
C
      LOGICAL SELECTED
      DATA SLOP/3.0/
C
      SELECTED = .FALSE.
      IF (ISEN.LT.1 .OR. ISEN.GT.2) GOTO 900
      IF (T.EQ.-999.) GOTO 900
C
      REM = AMOD(T,140.)
      IF (ISEN.EQ.1) THEN
         IF (REM.GT.140.-SLOP .OR. REM.LT.70.+SLOP) SELECTED = .TRUE.
         IF (T.LT.0.) SELECTED = .TRUE.
         ELSE
         IF (REM.GT.70.-SLOP .OR. REM.LT.SLOP) SELECTED = .TRUE.
         ENDIF
C
900   RETURN
      END


C     ----------------------------------------------------- 
      SUBROUTINE CWORD
C 
C     Allows user to manually edit data values.
C     Also allows application of relatively rare automatic
C     editing steps.
C     ----------------------------------------------------- 
C 
      PARAMETER (MXRC = 9000, NOPT = 16, NINVAR = 12, NPROC = 50) 
C 
      CHARACTER*20 INPUT
      CHARACTER*1 ANS 
      CHARACTER*4 STNS,STRING
      CHARACTER*1 OPTD(3)
      CHARACTER*4 VERSION
      LOGICAL NEEDOPT, AUTOTEMP
C 
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
      COMMON /OPTIONS/ NEEDOPT(NOPT)
      COMMON /PROCESSING/ NPSTP, PROC(NPROC)
C 
C 
      DATA OPTD /'L','T','Q'/
C 
C 
      BAD = -999.
      NDISP = NELM(25.,FALLRATE,DATARATE)
C
      NOPTD = 0
C
10    NOPTD = NOPTD+1
      IF (PROC(9).EQ.1. .AND. NOPTD.EQ.2) NOPTD = 3
      IF (NOPTD.GT.3) NOPTD = 3
      WRITE(LUT,'(///," Editing menu options:")')
      WRITE(LUT,'(" ------------------------------------")')
      WRITE(LUT,'(" (L) - Edit launch data")')
      WRITE(LUT,'(" (S) - Edit sounding data")')
      WRITE(LUT,'(" (P) - Substitute estimated pressures")')
      WRITE(LUT,'(" (T) - Apply dynamic temp correction")')
      WRITE(LUT,'(" (Q) - Quit this option")')
      WRITE(LUT,'(" ------------------------------------")')
      WRITE(LUT,'(/," Enter selection [",A1,"]: ",$)') OPTD(NOPTD)
C 
      READ(*,'(A)') ANS
      IF (ANS.EQ.' ') ANS = OPTD(NOPTD)
      IF (ANS.EQ.'P' .OR. ANS.EQ.'p') GOTO 100
      IF (ANS.EQ.'L' .OR. ANS.EQ.'l') GOTO 200
      IF (ANS.EQ.'S' .OR. ANS.EQ.'s') GOTO 300
      IF (ANS.EQ.'T' .OR. ANS.EQ.'t') GOTO 400
      IF (ANS.EQ.'Q' .OR. ANS.EQ.'q') RETURN
      NOPTD = 0
      GOTO 10 
C 
C 
C     Substitute estimated pressure for pressure
C     Recalculate ESTPR based on current PRL.
C     ------------------------------------------
C 
100   WRITE(LUT,'(/," *** MEASURED P VALUES WILL BE LOST ***")')
      WRITE(LUT,'(" Are you sure you want to do this (y/n)? ",$)') 
      READ(*,'(A)') ANS 
      IF (ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
C
         IF (NSNDTYPE.EQ.2) THEN
            PRLAST = PRL
            IF (PRLAST.EQ.-999.) THEN
               WRITE(LUT,'(/" *** CANNOT ESTIMATE PRESSURES.",
     *         "  LAUNCH PRESSURE IS BAD. ***")')
               NOPTD = 0
               GOTO 10
               ENDIF
            ENDIF
C
         DO 110 I = 1, IXEND
              IF (NSNDTYPE.EQ.2) THEN
                 ESTPRI = 2.0*(0.3 + (PRLAST*.000410))/DATARATE
                 IF (SNDDAT(I,1).LE.0.) ESTPRI = 0.0
                 IF (SNDDAT(I,1).LE.2.) ESTPRI = ESTPRI/2.0
                 SNDDAT(I,5) = PRLAST+ESTPRI
                 PRLAST = SNDDAT(I,5)
                 ENDIF
C
              SNDDAT(I,2) = SNDDAT(I,5) 
              ISNDFLG(I,1) = 0
110           CONTINUE
         ENDIF
C
      WRITE(LUT,'(/" Last estimated pressure is ",f6.1," mb.")') PRLAST
      NEEDOPT(12) = .TRUE.
      PROC(8) = 1.
      GOTO 10 
C 
C 
C     Edit launch parameters
C     ----------------------
C 
200   WRITE(LUT,'(//," Enter new value, (return) to keep value, or",
     *               " ''M'' to set missing:")')
201   WRITE(LUT,'(/," Verify YEAR (yy)     = ",I6," : ",$)') IYR
      READ(*,'(A)') INPUT 
      IF (INPUT.NE.' ') READ (INPUT,*,ERR=201) IYR
C 
205   WRITE(LUT,'(" Verify MONTH (mm)    = ",I6," : ",$)') IMO
      READ(*,'(A)') INPUT 
      IF (INPUT.NE.' ') READ (INPUT,*,ERR=205) IMO
C 
210   WRITE(LUT,'(" Verify DAY (dd)      = ",I6," : ",$)') IDY
      READ(*,'(A)') INPUT 
      IF (INPUT.NE.' ') READ (INPUT,*,ERR=210) IDY
C 
215   WRITE(LUT,'(" Verify TIME (hhmmss) = ",I6.6," : ",$)') TIML
      READ(*,'(A)') INPUT 
      IF (INPUT.NE.' ') READ (INPUT,*,ERR=215) TIML 
C 
220   LATD = INT(RLAT)
      RLATM = AMOD(RLAT,1.)*60.
      WRITE(LUT,'(" Verify LAT [N] (dddd mm.m) = ",
     *              I4,1x,f4.1," : ",$)') LATD,RLATM
      READ(*,'(A)') INPUT 
      IF (INPUT.NE.' ') READ (INPUT,*,ERR=220) LATD,RLATM
      RLAT = ABS(LATD)+RLATM/60.
      IF (LATD.LT.0) RLAT = -RLAT
C 
225   LOND = INT(RLON)
      RLONM = AMOD(RLON,1.)*60.
      WRITE(LUT,'(" Verify LON [W] (dddd mm.m) = ",
     *              I4,1x,f4.1," : ",$)') LOND,RLONM
      READ(*,'(A)') INPUT 
      IF (INPUT.NE.' ') READ (INPUT,*,ERR=225) LOND,RLONM
      RLON = ABS(LOND)+RLONM/60.
      IF (LOND.LT.0) RLON = -RLON
C 
230   WRITE(LUT,'(" Verify WD (deg)      = ",I6," : ",$)') WDL
      READ(*,'(A)') INPUT 
      IF (INPUT.EQ.'M' .OR. INPUT.EQ.'m') THEN
         WDL = BAD
         WSL = BAD
         GOTO 240
         ENDIF
      IF (INPUT.NE.' ') READ (INPUT,*,ERR=230) WDL
      IF (WDL.EQ.BAD) THEN
         WSL = BAD
         GOTO 240
         ENDIF
C 
235   WSLK = -999.
      IF (WSL.NE.BAD) WSLK = WSL*1.94
      WRITE(LUT,'(" Verify WS (kts)      = ",F6.1," : ",$)') WSLK
      READ(*,'(A)') INPUT 
      IF (INPUT.NE.' ') READ (INPUT,*,ERR=235) WSLK
      WSL = WSLK/1.94
C 
240   WRITE(LUT,'(" Verify T (deg C)     = ",F6.1," : ",$)') TEL
      READ(*,'(A)') INPUT 
      IF (INPUT.EQ.'M' .OR. INPUT.EQ.'m') THEN
         TEL = BAD
         NEEDOPT(12) = .TRUE. 
         GOTO 245
         ENDIF
      IF (INPUT.NE.' ') THEN
         READ (INPUT,*,ERR=240) TEL 
         NEEDOPT(12) = .TRUE. 
         ENDIF
C 
245   WRITE(LUT,'(" Verify Td (deg C)    = ",F6.1," : ",$)') DPL
      READ(*,'(A)') INPUT 
      IF (INPUT.EQ.'M' .OR. INPUT.EQ.'m') THEN
         DPL = BAD
         RHL = BAD
         NEEDOPT(12) = .TRUE. 
         GOTO 250
         ENDIF
      IF (INPUT.NE.' ')  THEN 
         READ (INPUT,*,ERR=245) DPL 
         NEEDOPT(12) = .TRUE. 
         ENDIF
C 
250   WRITE(LUT,'(" Verify PR (mb)       = ",F6.1," : ",$)') PRL
      READ(*,'(A)') INPUT 
      IF (INPUT.EQ.'M' .OR. INPUT.EQ.'m') THEN
         PRL = BAD
         NEEDOPT(12) = .TRUE. 
         GOTO 255
         ENDIF
      IF (INPUT.NE.' ')  THEN 
         READ (INPUT,*,ERR=250) PRL 
         NEEDOPT(12) = .TRUE. 
         ENDIF
C
      RHL =  RELHU(TEL,DPL,PRL) 
C 
255   WRITE(LUT,'(" Verify GA (m)        = ",F6.0," : ",$)') HTL 
      READ(*,'(A)') INPUT 
      IF (INPUT.EQ.'M' .OR. INPUT.EQ.'m') THEN
         HTL = BAD
         NEEDOPT(12) = .TRUE. 
         GOTO 290
         ENDIF
      IF (INPUT.NE.' ')  THEN 
         READ (INPUT,*,ERR=255) HTL 
         NEEDOPT(12) = .TRUE. 
         ENDIF
C 
290   GOTO 10 
C 
C 
C     Edit sounding data
C     ------------------
C 
300   WRITE(LUT,'(/" Current EOD at index = ",I4.4)') IXEND
301   WRITE(LUT,'(/" Enter index of the record to be edited",/,
     *     " (or ''L'' to list sounding by index numbers): ",$)') 
      READ(*,'(a)',ERR=301) STRING
      IF (STRING.EQ.'L' .OR. STRING.EQ.'l') THEN
         CALL SHOWSNDNG
         WRITE(LUT,'(" ")')
         GOTO 300
         ENDIF
C
      READ(STRING,*,ERR=300) IX
      IF (IX.LE.0 .OR. IX.GT.IXEND) THEN
         WRITE(LUT,'(" *** INVALID INDEX NUMBER ***")') 
         GOTO 300
         ENDIF
C 
C 
305   CALL SHOWLEVEL(IX,IX,0) 
      WRITE(LUT,'(  " Enter sounding parameter to be edited:")')
      WRITE(LUT,'(  " --------------------------------------")')
      WRITE(LUT,'(" (1) - Time",/,
     *            " (2) - Pressure",/,
     *            " (3) - Temperature",/, 
     *            " (4) - Relative humidity",/, 
     *            " (5) - Estimated pressure",/,
     *            " (6) - Wind direction",/,
     *            " (7) - Wind speed",/,
     *            " (8) - Wind quality",/,
     *            " (9) - Geopotential height",/, 
     *            " (E) - EDIT another data record",/,
     *            " (Q) - QUIT sounding edits")') 
      WRITE(LUT,'(  " --------------------------------------")')
      WRITE(LUT,'(/," Enter selection: ",$)')
C 
      READ(*,'(A1)',ERR=305) ANS
      IF (ANS.EQ.'Q' .OR. ANS.EQ.'q') GOTO 10
      IF (ANS.EQ.'E' .OR. ANS.EQ.'e') GOTO 300
      READ(ANS,*,ERR=305) IVAR
      IF (IVAR.GT.9) GOTO 305 
C 
      WRITE(LUT,'(/," Current value is ",F10.2)') 
     *      SNDDAT(IX,IVAR) 
      WRITE(LUT,'(" Enter new value: ",$)')
      READ (LUT,*) XVAL 
      IF (XVAL.NE.-999.) THEN 
         IF (IVAR.EQ.2) ISNDFLG(IX,1) = 0 
         IF (IVAR.EQ.3) ISNDFLG(IX,2) = 0 
         IF (IVAR.EQ.4) ISNDFLG(IX,3) = 0 
         IF (IVAR.EQ.5 .OR. IVAR.EQ.6) ISNDFLG(IX,4) = 0
         ENDIF
      IF (XVAL.NE.SNDDAT(IX,IVAR)) THEN 
         IF (NSITE.EQ.2) NEEDOPT(11) = .TRUE. 
         IF (IVAR.GE.2 .AND. IVAR.LE.4) NEEDOPT(12) = .TRUE.
         ENDIF
      IF (NSITE.EQ.1 .AND. XVAL.NE.SNDDAT(IX,IVAR)) THEN 
         IF (IVAR.EQ.2) ISNDFLG(IX,1) = 5 
         IF (IVAR.EQ.3) ISNDFLG(IX,2) = 5 
         IF (IVAR.EQ.4) ISNDFLG(IX,3) = 5 
         IF (IVAR.EQ.5 .OR. IVAR.EQ.6) ISNDFLG(IX,4) = 5
         ENDIF
      SNDDAT(IX,IVAR) = XVAL
      GOTO 305
C 
C
C
C     Apply dynamic temperature correction
C     ------------------------------------
400   CALL DYNTCORR(IERR)
      IF (IERR.EQ.0) PROC(9) = 1.
      NEEDOPT(11) = .TRUE.
      NEEDOPT(12) = .TRUE.
      GOTO 10
C
C
      END 
C 
C 
C 
C     ------------------------------------------------------
      SUBROUTINE FILLGAPS 
C 
C     Interpolates through gaps in data arrays.  Uses launch
C     data as initial point.  Interpolates all sounding 
C     arrays except geopotential (#9).  Before interpolating
C     data whose status is "bad" (1 or 2) is set to -999. 
C     ------------------------------------------------------
C 
      PARAMETER (MXRC = 9000, NINVAR = 12)
C 
      DIMENSION X(MXRC+1), Y(MXRC+1), Z(MXRC+1)
      DIMENSION IFLG(MXRC+1), ACVAL(8)
      CHARACTER*4 STNS
      CHARACTER*4 VERSION
      LOGICAL LAST
C 
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
C
      DATA PRDBTMAX/25./
      DATA BAD/-999./
C 
C 
      ACVAL(1) = 0
      ACVAL(2) = PRL
      ACVAL(3) = TEL
      ACVAL(4) = RHL
      ACVAL(5) = PRL
      ACVAL(6) = UCMP(WDL,WSL) 
      ACVAL(7) = VCMP(WDL,WSL) 
      ACVAL(8) = 9.0
      IF (NSNDTYPE.GT.1) ACVAL(8) = BAD
C 
C     Do PTH first
C     ------------
      DO 100 J = 1,5
         X(1) = ACVAL(J)
         IFLG(1) = 0
         DO 110 I = 1,IXEND 
              X(I+1) = SNDDAT(I,J)
              IF (J.EQ.1 .OR. J.EQ.5) IFLG(I+1) = 0 
              IF (J.EQ.2) IFLG(I+1) = ISNDFLG(I,1)
              IF (J.EQ.3) IFLG(I+1) = ISNDFLG(I,2)
              IF (J.EQ.4) IFLG(I+1) = ISNDFLG(I,3)
C 
C             Check status and set data to missing. 
C             ------------------------------------- 
              IF (IFLG(I+1).EQ.1 .OR. IFLG(I+1).EQ.2) X(I+1) = BAD
C 
110           CONTINUE
C 
         CALL INTRP(X,IFLG,IXEND+1) 
C 
         DO 120 I = 1,IXEND 
              SNDDAT(I,J) = X(I+1)
              IF (J.EQ.2) ISNDFLG(I,1) = IFLG(I+1)  
              IF (J.EQ.3) ISNDFLG(I,2) = IFLG(I+1)
              IF (J.EQ.4) ISNDFLG(I,3) = IFLG(I+1)
120           CONTINUE
100      CONTINUE 
C 
C 
C     Now do winds...convert to u,v components first
C     For ODWs interpolate WQ, for GPS, interpolate VVEL
C     --------------------------------------------------
      X(1) = ACVAL(6) 
      Y(1) = ACVAL(7) 
      Z(1) = ACVAL(8)
      IFLG(1) = 0 
      IF (X(1).EQ.BAD) Z(1) = BAD
C 
      DO 210 I = 1,IXEND
         IFLG(I+1) = ISNDFLG(I,4) 
         X(I+1) = UCMP(SNDDAT(I,6),SNDDAT(I,7))
         Y(I+1) = VCMP(SNDDAT(I,6),SNDDAT(I,7))
         IF (NSNDTYPE.EQ.1) THEN
            Z(I+1) = SNDDAT(I,8)
            ELSE
            Z(I+1) = SNDDAT(I,10)
            ENDIF
         IF (IFLG(I+1).EQ.1 .OR. IFLG(I+1).EQ.2) THEN 
              X(I+1) = BAD
              Y(I+1) = BAD
              Z(I+1) = BAD
              ENDIF
210      CONTINUE 
C 
      CALL INTRP(X,IFLG,IXEND+1)
      CALL INTRP(Y,IFLG,IXEND+1)
      CALL INTRP(Z,IFLG,IXEND+1)
C 
      DO 220 I = 1,IXEND
         SNDDAT(I,6) = WDCOMP(X(I+1),Y(I+1))    
         SNDDAT(I,7) = WSCOMP(X(I+1),Y(I+1))    
         IF (NSNDTYPE.EQ.1) THEN
            SNDDAT(I,8) = Z(I+1)
            ELSE
            SNDDAT(I,10) = Z(I+1)
            ENDIF
         ISNDFLG(I,4) = IFLG(I+1) 
220      CONTINUE 
C 
C
C     Automatically flag runs of PRDBTMAX mb as doubtful
C     --------------------------------------------------
      JYLIM = NELM(PRDBTMAX,FALLRATE,DATARATE)
      DO 300 IV = 1,4
    	 LAST = .FALSE.
	 DO 350 I = 1,IXEND
	    IF (ISNDFLG(I,IV).EQ.3) THEN
		IF (.NOT. LAST) ISTRT = I
		LAST = .TRUE.
		ELSE
		IF (LAST) THEN
		   ISTOP = I-1
		   IF (ISTOP-ISTRT .GE. JYLIM) THEN
			DO 360 J = ISTRT,ISTOP
360			   ISNDFLG(J,IV) = 4
			ENDIF
		   ENDIF
		LAST = .FALSE.
		ENDIF
350	    CONTINUE
300      CONTINUE

C
C
C     Reset flags for all missing data
C     --------------------------------
400   DO 410 I = 1,IXEND
         IF (SNDDAT(I,2).EQ.BAD) ISNDFLG(I,1) = 0
         IF (SNDDAT(I,3).EQ.BAD) ISNDFLG(I,2) = 0
         IF (SNDDAT(I,4).EQ.BAD) ISNDFLG(I,3) = 0
         IF (SNDDAT(I,6).EQ.BAD) ISNDFLG(I,4) = 0
410      CONTINUE
C
C
900   RETURN
      END 
C 
C 
C 
C     --------------------------------------------------------
      SUBROUTINE INTRP(X,IFLG,NT) 
C 
C     Interpolates through gaps in data arrays.  Looks at 
C     data flag first to see if a point is to be interpolated.
C     Based on second part of KCB routine DSPIK.
C     --------------------------------------------------------
C 
      DIMENSION X(NT), IFLG(NT) 
C 
      NB = 0                              ! Index of start of bad data 
C
      DO 500 M = 1,NT                     ! Loop through input data
         IF (IFLG(M).EQ.2) IFLG(M) = 0    ! Reset 'M' flag
         NE = 0                           ! Index of next good point 
C 
C        Check if this point is start of a gap
C        -------------------------------------
         IF (NB.EQ.0) THEN
              IF (X(M).EQ.-999.) NB = M 
              GOTO 500
              ENDIF 
C 
C        Check if this point is end of gap
C        ---------------------------------
         IF (X(M).EQ.-999.) GOTO 500       ! Still in data gap
C 
C        Point is good, gap over.  Do interpolation.
C        -------------------------------------------
         NE = M 
         IF (NB.EQ.1) THEN                 ! Bad from start,
              NB = 0                       ! don't do interpolation 
              GOTO 500
              ENDIF 
         AB = FLOAT(NB-1) 
         AE = FLOAT(NE) 
         AN = AB
         RG = ABS(AE-AB)
         NOW = NB - 1 
C 
50       AN = AN + 1.0
         NOW = NOW + 1
C 
C        Skip unless flag is 'R'.  Need to interpolate if flag is 'I'
C        also, because winds get passed through three times.
C        -------------------------------------------------------------
         IF (IFLG(NOW).EQ.1 .OR. IFLG(NOW).EQ.3) THEN
            R1 = ABS(AN-AB)/RG 
            R2 = 1.0 - R1
            X(NOW) = X(NB-1)*R2 + X(NE)*R1    ! Get interpolated value 
            IFLG(NOW) = 3                     ! Set flag to status 'I' 
            ENDIF
C 
         IF (NOW.LT.NE-1) GOTO 50             ! Go back for next bad pt
         NB = 0 
500      CONTINUE 
C
      RETURN
      END 
C 
C 
C 
C     --------------------------------------------------------- 
      SUBROUTINE SMOOTHSPL
C 
C     Program puts drop PTH or wind data through VICSPL filter.
C     --------------------------------------------------------- 
C 
C
	use thread_common
	 
      PARAMETER (MXRC = 9000, NOPT = 16, NINVAR = 12, NPROC = 50) 
C 
      CHARACTER*20 INPUT
      CHARACTER*1 VAR, OPTD(3)
C
      CHARACTER*4 STNS
      CHARACTER*4 VERSION
      LOGICAL NEEDOPT, AUTOTEMP 
      LOGICAL ECHO, DEBUG
C
      DIMENSION FOUT(2),FOUTD(2),XW(MXRC),DIFL(4)
C 
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
      COMMON /OPTIONS/ NEEDOPT(NOPT)
      COMMON /PROCESSING/ NPSTP,PROC(NPROC)
C 
      DATA ECHO/.TRUE./
      DATA OPTD /'P','W','Q'/
      DATA DEBUG/.FALSE./
      DATA DIFL/0.1,0.1,3.0,1.0/
C 
C 
C 
      DO 5 L = 1,IXEND
5	XW(L) = 1.0

      IF (AUTOTEMP) GOTO 100
      NOPTD = 0
10    NOPTD = NOPTD+1
      IF (NOPTD.GT.3) NOPTD = 1
      WRITE(LUT,'(/," Filter which parameter(s):")')
      WRITE(LUT,'(  " --------------------------")') 
      WRITE(LUT,'(  " (P) - PTH data")')
      WRITE(LUT,'(  " (W) - Wind data")') 
      WRITE(LUT,'(  " (Q) - Quit this option")')
      WRITE(LUT,'(  " --------------------------")') 
      WRITE(LUT,'(/," Enter selection [",A1,"]: ",$)') OPTD(NOPTD)
C 
	READ(*,'(A)') VAR 
      IF (VAR.EQ.' ') VAR = OPTD(NOPTD)
      CALL INDEX(VAR,IXDAT,IXFLG) 
      IF (IXDAT.EQ.0) RETURN
      IF (IXDAT.EQ.2) GOTO 100
      IF (IXDAT.EQ.6) GOTO 200
      GOTO 10 
C 
C 
C     Filter PTH
C     ----------
C 
100   IF (NSNDTYPE.EQ.1) THEN
          YDCWL = 240.0 
      endif

	  if (nsndtype .EQ. 2) then
		YDCWL = 10.0
	  endif

	  if (nsndtype .EQ.3) then
		YDCWL = 100.0
	  ENDIF

	if (nsndtype .EQ.4) then
		YDCWL = 30.0
	  ENDIF

      IF (YDCWL.LT.4.0/DATARATE) YDCWL=4.0/DATARATE
C
      KDAT = 1
      KYBC1 = 2
      KYBC2 = 2
      YBCWL1 = 2.0
      YBCWL2 = 2.0
C
      IF (AUTOTEMP) GOTO 110
105   WRITE(LUT,'(" Verify PTH filter cutoff (s) = ",F4.0," : ",
     *      $)') YDCWL
      READ(*,'(A)') INPUT 
      IF (INPUT.NE.' ') READ (INPUT,*,ERR=100) YDCWL 
110   WRITE (LUT,'(" PTH filter cutoff = ",F4.0," s."/)') YDCWL
C
      CALL VICSETUP(SNDDAT(1,1),IXEND,YDCWL,NX,YNB,YNT,2.0,IERR,ECHO)

	if (AUTOTEMP) then
		if (IERR.NE.0) then
			status_msg_thr = "Problem in PTH filter cuttoff"
			flag_thr = -1
			return
		endif
	endif

      IF (IERR.NE.0) GOTO 105
C 
      DO 175 J = 2,4
	 CALL VICSPL(SNDDAT(1,1),XW,SNDDAT(1,J),IXEND,MXRC,
     *               KDAT,YNB,YNT,NX,YDCWL,KYBC1,KYBC2,
     *               YBCWL1,YBCWL2,IERR)
         IF (IERR.EQ.1) THEN
            WRITE(LUT,'(
     *      " *** WARNING: FILTERING FAILED...CHECK DATA ***")')
            GOTO 175
            ENDIF
c
 	 DO 120 L = 2, IXEND-1
            IF (SNDDAT(L,J).NE.-999.) THEN
	       CALL SPOTVAL(SNDDAT(L,1),KDAT,FOUT,FOUTD)
C
C              Warning for excessive smoothing
C              -------------------------------
               IF (NSITE.EQ.1 .AND. DEBUG) THEN
                  DIFF = SNDDAT(L,J)-FOUT(1)
                  IF (ABS(DIFF).GT.DIFL(J-1)) THEN
                     WRITE(LUT,'(1x,A3,I1,I5,2x,3F10.2)') 'PTH',J-1,L,
     *               SNDDAT(L,J),FOUT(1),DIFF
                     ENDIF
                  ENDIF
C
	       SNDDAT(L,J) = FOUT(1)
	       ENDIF
            IF (J.EQ.4 .AND. SNDDAT(L,J).GT.100.) SNDDAT(L,J) = 100.
            IF (J.EQ.4 .AND. SNDDAT(L,J).LT.0.1
     *          .AND. SNDDAT(L,J).GT.-999.) SNDDAT(L,J) = 0.1
120         CONTINUE
175      CONTINUE 
C 
      PROC(1) = YDCWL
      NEEDOPT(12) = .TRUE.
      IF (AUTOTEMP) GOTO 200
      GOTO 10 
C 
C 
C     Filter winds
C     ------------
C 
200   IF (NSNDTYPE.EQ.1) THEN
          YDCWL = 240.0 
      endif

	if (nsndtype .EQ. 2) then
		YDCWL = 10.0
	endif

	if (nsndtype .EQ.3) then
		YDCWL = 100.0
	ENDIF

	if (nsndtype .EQ.4) then
		YDCWL = 30.0
	ENDIF

      IF (YDCWL.LT.4.0/DATARATE) YDCWL=4.0/DATARATE
C
      KDAT = 2
      KYBC1 = 2
      KYBC2 = 2
      YBCWL1 = 2.0
      YBCWL1 = 2.0
C
      IF (AUTOTEMP) GOTO 210
205   WRITE(LUT,'(" Verify wind filter cutoff (s) = ",f4.0," : ",$)')
     *      YDCWL
      READ(*,'(A)') INPUT 
      IF (INPUT.NE.' ') READ (INPUT,*,ERR=200) YDCWL 
210   WRITE (LUT,'(" Wind filter cutoff = ",F4.0," s."/)') YDCWL 
C 
      DO 215 L = 1,IXEND
         X(L) = UCMP(SNDDAT(L,6),SNDDAT(L,7))
         Y(L) = VCMP(SNDDAT(L,6),SNDDAT(L,7))
215      CONTINUE           
C 
      CALL VICSETUP(SNDDAT(1,1),IXEND,YDCWL,NX,YNB,YNT,2.0,IERR,ECHO)

	if (AUTOTEMP) then
		if (IERR.NE.0) then
			status_msg_thr = "Problem with wind filter cutoff"
			flag_thr = -1
			return
		endif
	endif

      IF (IERR.NE.0) GOTO 205
C 
C     Note---sending array X to VICSPL also sends array Y since KDAT=2
C     ----------------------------------------------------------------
      CALL VICSPL(SNDDAT(1,1),XW,X,IXEND,MXRC,
     *               KDAT,YNB,YNT,NX,YDCWL,KYBC1,KYBC2,
     *               YBCWL1,YBCWL2,IERR)
	if (AUTOTEMP) then
		if (IERR.NE.0) then
			status_msg_thr = "Wind filtering failed"
			flag_thr = -1
			return
		endif
	endif

      IF (IERR.EQ.1) THEN
         WRITE(LUT,'(
     *   " *** WARNING: FILTERING FAILED...CHECK DATA ***")')
         GOTO 10
         ENDIF
C
      DO 220 L = 2, IXEND-1
         IF (SNDDAT(L,6).NE.-999.) THEN
	       CALL SPOTVAL(SNDDAT(L,1),KDAT,FOUT,FOUTD)
C
C              Warning for excessive smoothing
C              -------------------------------
               IF (NSITE.EQ.1 .AND. DEBUG) THEN
                  DIFF1 = X(L)-FOUT(1)
                  DIFF2 = Y(L)-FOUT(2)
                  DIFF = (DIFF1**2.0+DIFF2**2.0)**0.5
                  IF (DIFF.GT.DIFL(4))
     *              WRITE(LUT,'(1x,A4,1X,I5,2x,F10.2)') 'WIND',L,DIFF
                  ENDIF
C
               SNDDAT(L,6) = WDCOMP(FOUT(1),FOUT(2))
               SNDDAT(L,7) = WSCOMP(FOUT(1),FOUT(2))
	       ENDIF
220      CONTINUE
C
      PROC(2) = YDCWL
      IF (AUTOTEMP) RETURN
      GOTO 10 
C
      END 
C 
C 
C 
C     --------------------------------------------------------- 
      SUBROUTINE SMOOTH 
C 
C     Program puts ODW PTH or wind data through low-pass filter.
C     Filter has default cutoff period at 40 seconds for PTH, 
C     240 s for wind. 
C     --------------------------------------------------------- 
C 
C 
      PARAMETER (MXRC = 9000, NOPT = 16, NINVAR = 12) 
C 
      DIMENSION IBGN(100),IEND(100) 
      CHARACTER*20 INPUT
      CHARACTER*1 VAR 
      CHARACTER*4 STNS
      CHARACTER*4 VERSION
      LOGICAL NEEDOPT, AUTOTEMP
C 
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
      COMMON /OPTIONS/ NEEDOPT(NOPT)
C 
C 
      IF (AUTOTEMP) GOTO 100
C
10    WRITE(LUT,'(/," Filter which parameter(s):")')
      WRITE(LUT,'(  " --------------------------")') 
      WRITE(LUT,'(  " (P) - PTH data")')
      WRITE(LUT,'(  " (W) - Wind data")') 
      WRITE(LUT,'(  " (Q) - Quit this option")')
      WRITE(LUT,'(  " --------------------------")') 
      WRITE(LUT,'(/," Enter selection: ",$)')
C 
      if (AUTOTEMP) then
		VAR = 'Q'
	else
		READ(*,'(A)') VAR 
	endif
      CALL INDEX(VAR,IXDAT,IXFLG) 
      IF (IXDAT.EQ.0) RETURN
      IF (IXDAT.EQ.2) GOTO 100
      IF (IXDAT.EQ.6) GOTO 200
      GOTO 10 
C 
C 
C     Filter PTH
C     ----------
C 
100   NTERMS = 25 
      INTERVAL = NINT(1./DATARATE)
      IF (NSNDTYPE.EQ.1) THEN
         CUTOFF = 40.0 
         ELSE
         CUTOFF = 20.0
         ENDIF
	if (nsndtype .ge. 3) cutoff = 20.0

      DO 105 L = 1,MXRC 
105      X(L) = 0.0 
C 
      IF (AUTOTEMP) GOTO 110
      WRITE(LUT,'(" Verify PTH filter cutoff (s) = ",I2," : ",
     *      $)') CUTOFF
      READ(*,'(A)') INPUT 
      IF (INPUT.NE.' ') READ (INPUT,*,ERR=100) CUTOFF 
110   WRITE (LUT,'(" PTH filter cutoff = ",F3.0," s."/)') CUTOFF
      WCUT = 2.0*FLOAT(INTERVAL)/CUTOFF 
C 
      DO 175 J = 2,4
         CALL GAP(SNDDAT(1,J),IXEND,IBGN,IEND,NGAPS,-999.)    
         DO 150 K = 1,NGAPS+1 
              NTG = IEND(K)-IBGN(K)+1 
              IF (NTG*INTERVAL .LT. INT(CUTOFF)) GOTO 150 
              CALL LOPASS(SNDDAT(IBGN(K),J),X,NTG,WCUT,NTERMS)
              DO 120 L = 1,NTG
                IF (J.EQ.4 .AND. X(L).GT.100.) X(L) = 100.
120             SNDDAT(IBGN(K)+L-1,J) = X(L)
150           CONTINUE
175      CONTINUE 
C 
      NEEDOPT(12) = .TRUE.
      IF (AUTOTEMP) GOTO 200
      GOTO 10 
C 
C 
C     Filter winds
C     ------------
C 
200   NTERMS = 25 
      INTERVAL = NINT(1./DATARATE)
      IF (NSNDTYPE.EQ.1) THEN
         CUTOFF = 240. 
         ELSE
         CUTOFF = 40.
         ENDIF

	if (nsndtype .ge. 3) cutoff = 40.0

      DO 205 L = 1,MXRC 
205      Z(L) = 0.0 
C 
      IF (AUTOTEMP) GOTO 210
      WRITE(LUT,'(" Verify wind filter cutoff (s) = ",F4.0," : ",$)')
     *      CUTOFF
      READ(*,'(A)') INPUT 
      IF (INPUT.NE.' ') READ (INPUT,*,ERR=200) CUTOFF 
210   WRITE (LUT,'(" Wind filter cutoff = ",F4.0," s."/)') CUTOFF 
      WCUT = 2.0*FLOAT(INTERVAL)/CUTOFF 
C 
      DO 215 L = 1,IXEND
         X(L) = UCMP(SNDDAT(L,6),SNDDAT(L,7))
         Y(L) = VCMP(SNDDAT(L,6),SNDDAT(L,7))
215      CONTINUE           
C 
      CALL GAP(X,IXEND,IBGN,IEND,NGAPS,-999.)     
      DO 250 K = 1,NGAPS+1
         NTG = IEND(K)-IBGN(K)+1
         IF (NTG*INTERVAL .LT. INT(CUTOFF)) GOTO 250
         CALL LOPASS(X(IBGN(K)),Z,NTG,WCUT,NTERMS)
         DO 220 L = 1,NTG 
220           X(IBGN(K)+L-1) = Z(L) 
250      CONTINUE 
C 
      CALL GAP(Y,IXEND,IBGN,IEND,NGAPS,-999.)     
      DO 255 K = 1,NGAPS+1
         NTG = IEND(K)-IBGN(K)+1
         IF (NTG*INTERVAL .LT. INT(CUTOFF)) GOTO 255
         CALL LOPASS(Y(IBGN(K)),Z,NTG,WCUT,NTERMS)
         DO 225 L = 1,NTG 
225           Y(IBGN(K)+L-1) = Z(L) 
255      CONTINUE 
C 
      DO 260 L = 1,IXEND
         SNDDAT(L,6) = WDCOMP(X(L),Y(L))
         SNDDAT(L,7) = WSCOMP(X(L),Y(L))
260      CONTINUE           
C 
      IF (AUTOTEMP) RETURN
      GOTO 10 
      END 
C 
C 
C 
C     ----------------------------------------------------- 
      SUBROUTINE HEIGHTS
C 
C     Program computes geopotential heights/surface press.
C     Can do computation from the top down, or surface up 
C     ----------------------------------------------------- 
C 
C 
      PARAMETER (MXRC = 9000, NINVAR = 12, NPROC = 50)
C 
      CHARACTER*20 INPUT
      CHARACTER*4 STNS, VERSION
      CHARACTER*1  HITWATER, BEGIN, AGAIN, HITD, BEGIND
      LOGICAL REPHUM, AUTOTEMP
C 
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
      COMMON /PROCESSING/ NPSTP, PROC(NPROC)
C 
      DATA HUMDEF/70./
      DATA BAD/-999./
C
C     Set up initial parameters 
C     ------------------------- 
      SFCPINIT = SNDDAT(IXEND,2)
      REPHUM = .TRUE.
C      IF (NSITE.EQ.2) REPHUM = .FALSE.
C
C
C     Show last value...is this sea level
C     -----------------------------------
      AGAIN = 'N' 
      WRITE(LUT,'(/," The last level in the file is...")')
      CALL SHOWLEVEL(IXEND,IXEND,0) 
      HITD = 'N'
c	if the last pressure was  > 975, or the calculated heights ere within about 20 m of the
c     surface, we assume that we it the surface. We also check the next to last height, since
c     the surface computed height is frequently missing.
      IF ((SNDDAT(IXEND,2)       .GT.975.)  .OR.
	1    (abs(SNDDAT(IXEND,9))  .LT. 20 )  .OR.
     2    (abs(SNDDAT(IXEND-1,9)).LT. 26 ) )     HITD = 'Y'	 
c
c	if glass sonde, then last point in data IS
c     the surface obs
c
      if (nsndtype .ge. 3) HITD = 'Y'
c
c     don't ask the question about sea surface. Instead,
c	  assume that it is if P > 930, else it is not.
c
c      WRITE(LUT,'(" Is this the sea surface [",A1,"]? ",$)') HITD
c      READ(*,'(A)') HITWATER
	  HITWATER = ' '
      IF (HITWATER.EQ.' ') HITWATER = HITD
      IF (HITWATER.EQ.'y') HITWATER = 'Y'
      IF (HITWATER.EQ.'Y' .AND. nsndtype .lt. 3) THEN
         SNDDAT(IXEND,9) = 0.0
         ELSE
         PROC(6) = -999.
         ENDIF
C 
C 
C     Set starting point for hydrostatic calculation
C     ----------------------------------------------
C
      BEGIND = 'F'
      IF (AUTOTEMP) THEN
	     IF (HITWATER .EQ. 'Y') THEN
		    BEGIND = 'S'
		 ELSE
		    BEGIND = 'F'
		 ENDIF
	  ENDIF
C
210   IF (AUTOTEMP) THEN
         BEGIN = BEGIND
         ELSE
         WRITE(LUT,'(/" Specify height calculation option:")') 
         WRITE(LUT,'( " ----------------------------------------")') 
         WRITE(LUT,'( " (F) - Compute heights, from Flight level")') 
         WRITE(LUT,'( " (S) - Compute heights, from Surface")')
         WRITE(LUT,'( " (M) - Set heights missing (",F6.0,")")') BAD
         WRITE(LUT,'( " (Q) - Exit and return to main menu")') 
         WRITE(LUT,'( " ----------------------------------------")') 
         WRITE(LUT,'(/" Enter selection [",A1,"]: ",$)') BEGIND
         READ(*,'(A)') BEGIN 
         IF (BEGIN.EQ.' ') BEGIN = BEGIND
         IF (BEGIN.EQ.'f') BEGIN = 'F' 
         IF (BEGIN.EQ.'s') BEGIN = 'S' 
         IF (BEGIN.EQ.'m') BEGIN = 'M' 
         IF (BEGIN.EQ.'q') BEGIN = 'Q' 
         ENDIF
C
      IF (BEGIN.EQ.'Q') RETURN
C
      IF (BEGIN.NE.'F'.AND.BEGIN.NE.'S'
     *    .AND.BEGIN.NE.'M'.AND.BEGIN.NE.'Q') THEN
	 WRITE(LUT,'(/" *** Invalid selection ***")')
         GOTO 210
      	 ENDIF
C
      IF (BEGIN.EQ.'S' .AND. HITWATER.NE.'Y') THEN
         IF (AUTOTEMP) BEGIND = 'M'
         WRITE(LUT,'(/" *** WARNING ***",/,
     *   " You cannot compute heights from the surface up",/,
     *   " for a sonde that did not reach the surface!"//)') 
         GOTO 210
         ENDIF
C 
      IF (BEGIN.EQ.'M') THEN
         DO 215 L = 1,IXEND 
215           SNDDAT(L,9) = BAD 
         IF (HITWATER.EQ.'Y'.AND. nsndtype .lt. 3) SNDDAT(IXEND,9) = 0.
         RETURN 
         ENDIF
C
C
C     OK, now we actually are going to calculate some heights
C     Set parameters depending on the direction of integration
C     --------------------------------------------------------
      IF (BEGIN.EQ.'F') THEN
         Z1 = HTL 
         T1 = TEL 
         H1 = RHCK(RHL,REPHUM,HUMDEF,BAD)
         P1 = PRL 
         IBGN = 1 
         IEND = IXEND 
         ISTEP = 1
         PROC(7) = 1.
         ELSE 
         Z1 = 0.0 
c
c	for GLASS, we will start with the altitude of the
c     surface record
C
	if (nsndtype .ge. 3) then
		Z1 = snddat(IXEND,9)
	endif

         T1 = SNDDAT(IXEND,3) 
         H1 = RHCK(SNDDAT(IXEND,4),REPHUM,HUMDEF,BAD)
         P1 = SNDDAT(IXEND,2) 
         IBGN = IXEND-1 
         IEND = 1 
         ISTEP = -1 
         PROC(7) = 2.
         ENDIF
C 
C 
C     Reset previous heights to missing
C     ---------------------------------
      DO 240 L = 1,IXEND 
240      SNDDAT(L,9) = BAD 

C	for GLASS sonde, we want to preserve the surface altitude
	if (nsndtype .ge. 3) SNDDAT(IXEND,9) = Z1
C
C
C     Loop through sounding
C     ---------------------
      DO 250 K = IBGN,IEND,ISTEP
C 
         T2 = SNDDAT(K,3)
         H2 = RHCK(SNDDAT(K,4),REPHUM,HUMDEF,BAD)
         P2 = SNDDAT(K,2)
C 
         IF (K.EQ.IEND.AND.BEGIN.EQ.'F'.AND.HITWATER.EQ.'Y') THEN
            Z2 = 0.0
            SFCPR = HYDROP(P1,T1,H1,Z1,T2,H2,Z2,2,BAD)
            ELSE
            Z2 = HYDROZ(P1,T1,H1,Z1,P2,T2,H2,2,BAD) 
            ENDIF 
C 
        SNDDAT(K,9) = Z2
C
C       If calculation failed, just go on to next level
C       -----------------------------------------------
        IF (Z2.EQ.BAD) GOTO 250
C
C       Otherwise, reset point 1 and continue
C       -------------------------------------
        P1 = P2 
        H1 = RHCK(H2,REPHUM,HUMDEF,BAD)
        T1 = T2 
        Z1 = Z2 
250     CONTINUE
C 
C 
C     ----------------------------------------------------------
C     Height calculation terminated...cleanup and adjust profile
C     Adjust pressure profile...only for top-down integrations
C     and only if the sounding hit the surface.
C     ----------------------------------------------------------
C 
300   IF (BEGIN.EQ.'S') THEN
         T1 = T2
         H1 = H2
         P1 = P2
         Z1 = Z2
         T2 = TEL 
         H2 = RHCK(RHL,REPHUM,HUMDEF,BAD)
         P2 = PRL 
         ZAC = HYDROZ(P1,T1,H1,Z1,P2,T2,H2,2,BAD) 
         WRITE(LUT,395) PRL,HTL,ZAC
395      FORMAT(/,' At flight level of ',f5.1,' mb,',
     *            ' aircraft measured GA is ',f5.0,' m.',/,1x,
     *        20('.'),' the hydrostatic estimated GA is ',f5.0,' m.',/)
         IF (ZAC.NE.BAD .AND. HTL.NE.BAD) THEN
            DIFZ = HTL-ZAC
            DIFP = ABS((P2-P1)/(ZAC-Z1))*DIFZ
            WRITE(LUT,396) DIFZ,DIFP
396         FORMAT(' This height error of ',f5.0,' m is equivalent',
     *      ' to a pressure error of ',f6.1, ' mb.',//)
            ENDIF
         ENDIF
C
      IF (HITWATER.EQ.'Y' .AND. nsndtype .lt. 3) SNDDAT(IXEND,9) = 0.0
      IF (BEGIN.EQ.'S' .OR. HITWATER.NE.'Y') GOTO 500 
C 
C
C     Display current and calculated sfc pressures
C     --------------------------------------------
      WRITE(LUT,*)
      IF (SFCPINIT.NE.SNDDAT(IXEND,2)) WRITE(LUT,
     * '(" Initial surface pressure was: ",F6.1," mb")') SFCPINIT
      WRITE(LUT,390) SNDDAT(IXEND,2), SFCPR 
390   FORMAT(' Current surface pressure is:  ',F6.1,' mb',
     *      /,' Hydrostatic estimate is:      ',F6.1,' mb')
      HYDPR = SFCPR
      IF (SFCPR.EQ.BAD) SFCPR = SNDDAT(IXEND,2) 
C
      IF (AUTOTEMP) GOTO 310
301   WRITE(LUT,'(/" Verify surface pressure = ",F6.1,": ",$)') SFCPR
      READ(*,'(A)') INPUT 
      IF (INPUT.NE.' ') READ(INPUT,*,ERR=301) SFCPR 
C
      IF (SFCPR.LT.850.) THEN 
        WRITE(LUT,'(/,1X,F6.1," is not a valid pressure value!")') SFCPR
        WRITE(LUT,'(" Check value and re-enter surface pressure."/)')   
        GOTO 301
        ENDIF 
C 
C     Determine corrections
C     ---------------------
310   PCORR = SFCPR - SNDDAT(IXEND,2) 
      HCORR = SFCPR - HYDPR
      DO 350 K = 1,IXEND
         IF (SNDDAT(K,2).EQ.BAD) GOTO 350
         SNDDAT(K,2) = SNDDAT(K,2) + PCORR*FLOAT(K)/FLOAT(IXEND)
350      CONTINUE 
C 
C
C     Determine what needs to be done next and Return to HEIGHTS menu
C     ---------------------------------------------------------------
C 
500   IF (.NOT. AUTOTEMP) THEN
         WRITE(LUT,'(/," Last 3 levels of data are now...")') 
         CALL SHOWLEVEL(IXEND-2,IXEND,0)
         ENDIF
C
      BEGIND = 'Q'
      IF (BEGIN.EQ.'S' .OR. HITWATER.NE.'Y') GOTO 210
C
      IF (ABS(PCORR).GE.0.1 .AND. HYDPR.EQ.SFCPR) THEN
         BEGIND = 'F'
         WRITE(LUT,'(" *** WARNING ***",//,
     *     " HEIGHTS SHOULD BE RECALCULATED FROM FLIGHT LEVEL",/, 
     *     " DOWN TO REFINE THE HYDROSTATIC SURFACE ESTIMATE.")') 
         ENDIF
C
      IF (ABS(HCORR).GE.0.3 .AND. HYDPR.GT.0.) THEN
         BEGIND = 'S'
         WRITE(LUT, '(" *** WARNING ***",//,
     *     " HEIGHTS MUST BE RECALCULATED FROM THE SURFACE UP",/, 
     *     " TO BE CONSISTENT WITH THE CHOSEN SURFACE PRESSURE.")') 
         ENDIF
C 
      GOTO 210
      END 
C 
C 
C 
C     --------------------------------------
      FUNCTION RHCK(RH,REPHUM,HUMDEF,BAD)
C   
C     Replaces missing RH with default value
C     --------------------------------------
C
      LOGICAL REPHUM
C
      IF (REPHUM .AND. RH.EQ.BAD) THEN
         RHCK = HUMDEF
         ELSE
         RHCK = RH
         ENDIF
      RETURN
      END
C
C
C     ------------------------------------------------- 
      SUBROUTINE ODEND(NRECS)
C 
C     Sets end of drop marker 
C     ------------------------------------------------- 
C 
C 
      PARAMETER (MXRC = 9000, NOPT = 16, NINVAR = 12, NPROC = 50) 
      PARAMETER (NEOD = 100)
C 
      CHARACTER*1 ANS, ANSD
      CHARACTER*4 STNS,STRING
      CHARACTER*4 VERSION
      LOGICAL NEEDOPT, AUTOTEMP
C 
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
      COMMON /OPTIONS/ NEEDOPT(NOPT)
      COMMON /EODDATA/ EOD(NEOD,NINVAR)
      COMMON /PROCESSING/ NPSTP, PROC(NPROC)
C 
C 
      BAD = -999.
      FILERATE = 2.0      
      NDISP = NEOD
      IF (NSNDTYPE.EQ.1) NDISP = 22
      OLDEND = SNDDAT(IXEND,1)
      IENDREC = 0
      IENDWIND = 0

C
C
C     Try and find default EOD: first look for last wind
C     --------------------------------------------------
      DO 105 L = IXEND,1,-1
         IF (SNDDAT(L,6).EQ.BAD) GOTO 105
         IENDREC = L
         IENDWIND = L
         GOTO 110
105      CONTINUE
C
C     Override with last record w/PTH data
C     ------------------------------------
110   DO 115 L = IXEND,1,-1
         IF (SNDDAT(L,2).EQ.BAD .AND.
     *       SNDDAT(L,3).EQ.BAD .AND.
     *       SNDDAT(L,4).EQ.BAD) GOTO 115
         IENDREC = L
         GOTO 120
115      CONTINUE
C
C     There are no valid records in the file
C     --------------------------------------
      IF (IENDREC.EQ.0) THEN
         WRITE(LUT,'(/" *** THIS SONDE HAS NO PTH OR WIND DATA ***"/)')
         GOTO 500
         ENDIF
C
C     
C     Branch for automatic operation
C     ------------------------------
120   IF (AUTOTEMP) THEN
         IXEND = IENDREC
         RETURN
         ENDIF
C
C
C     Display end portion of drop and verify findings
C     -----------------------------------------------
      CALL SHOWLEVEL(IENDREC-20,IENDREC+5,IENDREC)
125   WRITE(LUT,'(" Total number of data records in sounding = ",
     *      I4.4)') NRECS
      WRITE(LUT,'(" End of drop marker is now",
     *            " at index (IX)  = ",I4.4)') IXEND
      WRITE(LUT,'(/" Verify end of drop at IX = ",i4.4,/,
     *     " (or enter ''L'' to list sounding by index numbers): ",
     *     $)') IENDREC
      READ(*,'(a)',ERR=120) STRING
      IF (STRING.EQ.'L' .OR. STRING.EQ.'l') THEN
         CALL SHOWSNDNG
         WRITE(LUT,'(" ")')
         GOTO 125
         ENDIF
C
      IF (STRING.NE.' ') READ(STRING,*,ERR=125) IENDREC
      IF (IENDREC.EQ.0 .OR. IENDREC.GT.NRECS) THEN 
         WRITE(LUT,'(/," INVALID END OF DROP INDEX:",I4,/)') IENDREC
         GO TO 125
	 ENDIF
C
      IF (IENDREC.GT.IXEND) THEN 
         WRITE(LUT,'(/," Are you sure? ",$)')
	 READ(*,'(A)') ANS
	 IF (ANS.EQ.'Y' .OR. ANS.EQ.'y') GOTO 130
	 GOTO 120
         ENDIF
C
130    IXEND = IENDREC
C
C
C     If using full resolution, check to see if we need to slide
C     wind data up or down due to time uncertainty matching to PTH
C     ------------------------------------------------------------
      IWDIFF = IXEND-IENDWIND
      IF (DATARATE.EQ.FILERATE .AND. IWDIFF.NE.0 .AND.
     *   ABS(IWDIFF).LE.3) THEN
         ANSD = 'N'
         IF (IWDIFF.LT.0) ANSD = 'N'
         WRITE(LUT,'(/," Slide wind profile to end at ",I4.4," [",
     *      a1,"]? ",$)') IXEND, ANSD
         READ(*,'(A)') ANS
         IF (ANS.EQ.' ') ANS = ANSD
         IF (ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
C
            IF (IWDIFF.GT.0) THEN
               IBGN1 = IXEND
               IEND1 = 1 + IWDIFF
               INC = -1
               IBGN2 = IEND1-1
               IEND2 = 1
               ELSE
               IBGN1 = 1
               IEND1 = IXEND
               INC = 1
               IBGN2 = IXEND+1
               IEND2 = IXEND-IWDIFF
               ENDIF
C
            DO 140 L = IBGN1,IEND1,INC
               SNDDAT(L,6) = SNDDAT(L-IWDIFF,6)
               SNDDAT(L,7) = SNDDAT(L-IWDIFF,7)
               SNDDAT(L,8) = SNDDAT(L-IWDIFF,8)
               SNDDAT(L,10) = SNDDAT(L-IWDIFF,10)
               SNDDAT(L,11) = SNDDAT(L-IWDIFF,11)
               SNDDAT(L,12) = SNDDAT(L-IWDIFF,12)
               ISNDFLG(L,4) = ISNDFLG(L-IWDIFF,4)
140            CONTINUE
            DO 150 L = IBGN2,IEND2,INC
               SNDDAT(L,6) = BAD
               SNDDAT(L,7) = BAD
               SNDDAT(L,8) = BAD
               SNDDAT(L,10) = BAD
               SNDDAT(L,11) = BAD
               SNDDAT(L,12) = BAD
               ISNDFLG(L,4) = 0
150            CONTINUE
C
            ENDIF
         ENDIF
C 
C
C     Check whether to refine splash time
C     -----------------------------------
C 
200   REFEOD = SNDDAT(IXEND,1)
      DO 210 L = 1,NEOD
	 IF (EOD(L,1).EQ.REFEOD) GOTO 300
210	 CONTINUE
C
C     No match, can't do refinement, just exit
C     ----------------------------------------
C
      IF (DATARATE.NE.FILERATE) 
     *    WRITE(LUT,'(/," Splash time refinement is unavailable.")')
      GOTO 500
C
C
C     We have a match, now check with user to get refined end time
C     Try to find last good record.
C     ------------------------------------------------------------
C
300   IF (NSNDTYPE.EQ.1) THEN
         LDFLT = L
         ELSE
         SKIP = FILERATE/DATARATE
         NSKIP = NINT(FILERATE/DATARATE)
         NBEG = L+1
         NEND = L+NSKIP-1
         IF (NEND.GT.NEOD) NEND = NEOD
         DO 305 K = NBEG,NEND
            IF (EOD(K,2).EQ.BAD) GOTO 310
305         CONTINUE
C
310      LDFLT = K-1
         ENDIF
C
      CALL SHOWEOD(NDISP,LDFLT,LX)
C
C     --------------------------------------------------------
C     Now check refined time to see if IXEND must be changed.
C     Replace data values of IXEND record with refined values.
C     Time value is also changed.
C     --------------------------------------------------------
C
400   DO 410 L = NRECS,1,-1
	 IF (EOD(LX,1).EQ.SNDDAT(IXEND,1)) GOTO 500
	 IF (EOD(LX,1).GT.SNDDAT(L,1)) THEN
	    IXEND = L+1
            IF (IXEND.GT.NRECS) NRECS = IXEND
            DO 420 K = 1,NINVAR
		SNDDAT(IXEND,K) = EOD(LX,K)
420		CONTINUE
	    GOTO 500
	    ENDIF
410	 CONTINUE
C
C
C
C     Exit
C     ----
C
500   WRITE(LUT,'(/," End of drop is now at IX = ",I4.4,
     *          ".  Last 5 levels of data will be...")') IXEND
      CALL SHOWLEVEL(IXEND-4,IXEND,0)
      WRITE(LUT,'(" Hit return for main menu... ",$)')
      READ(*,'(A)') ANS 
      PROC(6) = SNDDAT(IXEND,2)
      IF (SNDDAT(IXEND,1).NE.OLDEND) THEN 
         NEEDOPT(12) = .TRUE. 
         NEEDOPT(13) = .TRUE. 
         NEEDOPT(14) = .TRUE. 
         ENDIF
      RETURN
      END 
C 
C 
C 
C     ------------------------------------
      SUBROUTINE PTHADD(PBIAS) 
C 
C     Corrects PTH offsets
C     ------------------------------------
C 
      PARAMETER (MXRC = 9000, NINVAR = 12, NPROC = 50)
C 
      CHARACTER*1 ANS, OPTD
      CHARACTER*4 STNS
      CHARACTER*20 INPUT
C 
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
      COMMON /PROCESSING/ NPSTP, PROC(NPROC)
C 
C 
      OPTD = 'P'
10    WRITE(LUT,'(/," Enter variable (PTH) to adjust or Q to quit [",
     *      A1,"]: ",$)') OPTD
      READ(*,'(A1)') ANS
      IF (ANS.EQ.' ') ANS = OPTD
      CALL INDEX(ANS,IXDAT,IXFLG) 
      IF (IXDAT.EQ.0) RETURN
      IF (IXDAT.LT.2 .OR. IXDAT.GT.4) GOTO 10 
C 
      CORR = 0.0
      IF (IXDAT.EQ.2) CORR = -PBIAS
C
20    WRITE(LUT,'(" Enter additive correction: [",F5.1,"]: ",$)') CORR
      READ(*,'(A)') INPUT 
      IF (INPUT.NE.' ') READ (INPUT,*,ERR=20) CORR
C 
C 
      DO 50 L = 1,IXEND 
         IF (SNDDAT(L,IXDAT).NE.-999.)  
     *       SNDDAT(L,IXDAT) = SNDDAT(L,IXDAT) + CORR 
50       CONTINUE 
      PROC(IXDAT+1) = PROC(IXDAT+1)+CORR
      RETURN
      END 
C 
C 
C 
C     ----------------------------------------------------- 
      SUBROUTINE SHOWSNDNG
C 
C     Allows user to display part or all of tabulated
C     listing.
C     ----------------------------------------------------- 
C 
      PARAMETER (MXRC = 9000) 
C 
      CHARACTER*1 ANSD, ANS, MODE
C 
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
C 
C
      ANSD = 'N'
      WRITE(LUT,'(/," List entire sounding [",A1,"]? ",$)') ANSD
      READ(*,'(A)') ANS
      IF (ANS.EQ.' ') ANS = ANSD
      IF (ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
         CALL SNDNG(.TRUE.,.FALSE.)
         ELSE
         ANSD = 'T'
         WRITE(LUT,'(" Select interval to list by ",
     *         "(t)ime or (p)ressure [",A1,"]? ",$)') ANSD
         READ(*,'(A)') MODE
         IF (MODE.EQ.' ') MODE = ANSD
C
         IF (MODE.EQ.'P' .OR. MODE.EQ.'p') THEN
110        WRITE(LUT,'(" Enter top of interval to list (mb): ",$)')
           READ(*,*,ERR=110) PTOP              
120        WRITE(LUT,'(" Enter btm of interval to list (mb): ",$)')
           READ(*,*,ERR=120) PBTM
           IXT = IXPRESS(PTOP)
           IXB = IXPRESS(PBTM)
           ELSE
130        WRITE(LUT,'(" Enter top of interval to list (s): ",$)')
           READ(*,*,ERR=130) PTOP              
140        WRITE(LUT,'(" Enter btm of interval to list (s): ",$)')
           READ(*,*,ERR=140) PBTM
           IXT = IXTIME(PTOP)
           IXB = IXTIME(PBTM)
           ENDIF
C
         CALL SHOWLEVEL(IXT,IXB,0)
         ENDIF
C 
      RETURN
      END 
C 
C 
C 
C     ----------------------------------------------------
      FUNCTION IXTIME(TIMX) 
C 
C     Function returns index closest to time TIMX
C     ----------------------------------------------------
C 
      PARAMETER (MXRC = 9000, NINVAR = 12)
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
C 
      IXTIME = 0 
      DIFF = 99999. 
      DO 301 I = 1,IXEND
         IF (ABS(SNDDAT(I,1)-TIMX).LT.DIFF) THEN 
              DIFF = ABS(SNDDAT(I,1)-TIMX) 
              IXTIME = I 
              ENDIF 
301      CONTINUE 
      RETURN
      END 
C 
C 
C 
C     ----------------------------------------------------
      FUNCTION IXPRESS(PRX) 
C 
C     Function returns index closest to pressure PRX
C     ----------------------------------------------------
C 
      PARAMETER (MXRC = 9000, NINVAR = 12)
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
C 
      IXPRESS = 0 
      DIFF = 99999. 
      DO 301 I = 1,IXEND
         IF (ABS(SNDDAT(I,2)-PRX).LT.DIFF) THEN 
              DIFF = ABS(SNDDAT(I,2)-PRX) 
              IXPRESS = I 
              ENDIF 
301      CONTINUE 
      RETURN
      END 
C 
C 
C 
C     ---------------------------------------
      SUBROUTINE SHOWLEVEL(NSTRT,NSTOP,IMARK) 
C 
C     Displays piece of data by time
C     ---------------------------------------
C 
      PARAMETER (MXRC = 9000, NINVAR = 12)
C 
      CHARACTER*1 CFLG
      CHARACTER*1 ANS 
      CHARACTER*1 MARK
      CHARACTER*4 STNS
      LOGICAL PAUSE
C 
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
C
      DATA PAUSE/.FALSE./
C 
C 
      IF (NSTRT.LT.1) NSTRT = 1 
      IF (NSTOP.GT.IXEND) NSTOP = IXEND 
C 
      WRITE(LUT,'(/)')
      WRITE(LUT,320)
320   FORMAT('   IX  TIME      PR       TE      RH        Z',
     *'      WD      WS      WQ     VV',/, 
     *       ' ----------------------------------------------', 
     *'-----------------------------')  
C 
      K = NSTRT-1
      KK = 0
325   K = K+1
C
      MARK = ' '
      IF (K.EQ.IMARK) MARK = '>'
      WRITE(LUT,330) MARK,K,SNDDAT(K,1),SNDDAT(K,2),CFLG(ISNDFLG(K,1)), 
     *               SNDDAT(K,3),CFLG(ISNDFLG(K,2)),SNDDAT(K,4), 
     *               CFLG(ISNDFLG(K,3)),NINT(SNDDAT(K,9)),
     *               NINT(SNDDAT(K,6)),SNDDAT(K,7),
     *               CFLG(ISNDFLG(K,4)),
     *               NINT(SNDDAT(K,8)),SNDDAT(K,10)
330   FORMAT(A1,I4,1X,F6.1,2X,F6.1,1X,A1,F7.1,1X,A1,F6.1,1X,A1,
     *          I7,3X,I4,2X,F6.1,1X,A1,2X,I4,1X,F6.1)
         KK = KK+1
         IF (MOD(KK,30).EQ.0 .AND. K.LT.NSTOP .AND. PAUSE) THEN 
              WRITE(LUT,
     *        '(/," Hit return for more (or (Q) to quit)... ",$)') 
              READ(*,'(A)') ANS 
              IF (ANS.EQ.'Q' .OR. ANS.EQ.'q') RETURN
              WRITE(LUT,320)
              K = K-1 
              ENDIF 
C
      IF (K.LT.NSTOP) GOTO 325
C 
      WRITE(LUT,'(1x,75("-"),/)')
      RETURN
      END 
C 
C 
C 
C     ------------------------------------
      SUBROUTINE SHOWEOD(NDISP,IDFLT,IXX)
C 
C     Displays end of drop refined data.
C     ------------------------------------
C 
      PARAMETER (NEOD = 100, NINVAR = 12)
C 
      CHARACTER*20 INPUT
      CHARACTER*1 LABEL
C
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /EODDATA/ EOD(NEOD,NINVAR)
C 
C
      WRITE(LUT,50)
50    FORMAT(//," High resolution data for splash refinement: ") 
C 
      IF (NDISP.EQ.22) THEN
         NBEG = 1
         NEND = 22
         ELSE
         NBEG = IDFLT-10
         NEND = IDFLT+11
         IF (NBEG.LT.1) NBEG = 1
         IF (NEND.GT.NDISP) NEND = NDISP
         ENDIF
C
300   WRITE(LUT,'(" ")')
      WRITE(LUT,320)
320   FORMAT(' IXX   TIME     PR       TE      RH     WD     WS      Z'
     *     ,/,
     *     ' --------------------------------------------------------') 
C
      DO 325 K = NBEG,NEND
         KIX = K-NBEG+1
         LABEL = ' '
         IF (K.EQ.IDFLT) LABEL = '>'
         WRITE(LUT,330) LABEL,KIX,(EOD(K,L),L=1,4),NINT(EOD(K,6)),
     *                  EOD(K,7),NINT(EOD(K,9))
C         IF (K.EQ.IDFLT) WRITE(LUT,335)
325      CONTINUE
330   FORMAT(1X,A1,I2,4(2X,F6.1),2X,I4,2X,F6.1,2X,I5)
335   FORMAT(1X,56("-"))
      WRITE(LUT,335)
C
C
350   KDFLT = IDFLT-NBEG+1
      WRITE(LUT,'(/," Verify end of drop at IXX = ",I2," : ",$)') KDFLT
      READ(*,'(A)') INPUT 
      IF (INPUT.NE.' ') READ (INPUT,*,ERR=350) KDFLT
C
      IF (KDFLT.LT.1 .OR. KDFLT.GT.NEND-NBEG+1) THEN
         WRITE(LUT,'(/," *** INVALID IXX ***")') 
         GOTO 300
         ENDIF
C
      IXX = KDFLT+NBEG-1
      RETURN
      END 
C 
C 
C 
C     --------------------------------------------------------- 
      FUNCTION NELM(PR,FALLRATE,DATARATE)
C 
C     Returns number of sounding data intervals corresponding
C     to a pressure depth PR (mb).
C     --------------------------------------------------------- 
C 
      PERIOD = 1./DATARATE
      PRINT = PERIOD*FALLRATE
      NELM = NINT(PR/PRINT)
      RETURN
      END 
C
C
C
C     ----------------------------------------------------- 
      SUBROUTINE TEMPDROP(FNAME)
	
C     Writes out WMO TEMP-DROP message: mandatory and sig levels.
C     James L. Franklin, NOAA/Hurricane Research Division
C     --------------------------------------------------------------------
C
C
C     Parameters:
C     --------------------------------------------------------------------
C     MODE_WMO  = 1 for NOAA encoding
C               = 2 for USAF encoding
C     CMTPROMPT = .TRUE. to prompt for comments.
C     LUG       = Graphics device for skew-t (0 to bypass)
C     --------------------------------------------------------------------
C
C
C     Revision history:
C     --------------------------------------------------------------------
C     04/03/98          Structural rewrite to accomodate inclusion of
C                       regional (51515) groups in Part A as well as Part
C                       B.  New subroutines GET_51515, CODE_51515, and 
C                       CODE_61616.
C                       Due to inclusion of 51515 in part A, we now permit
C                       extrapolated surface pressures in part A.
C
C     04/01/98          Pressure of max wind level in part A is now
C                       calculated as EXP(LOG(DBLE(P))), to ensure
C                       consistency with value chosen in part B.
C                       Comments now placed after parts A and B.
C
C     03/26/98          Changes RHMIN from 2.0 to 0.1%.
C
C     03/10/98          Adds sig level within 20 mb of surface, if one
C                       has not been previously chosen.
C
C     02/13/98          Will not encode max wind group if only available
C                       wind above 500 mb is the flight-level wind.
C
C     02/11/98          Aircraft identifiers obtained from acft_ctrl.dat
C                       Geographical ID (NT,PA,PN) from MANOPAA.
C
C     02/09/98          Adds CMTPROMPT parameter to prompt for comments.
C
C     02/06/98          Inversions above the first trop or 300 mb,
C                       whichever is higher, are not coded.
C
C     02/03/98          Includes 61616 mission ID line in Parts A and B.
C                       Includes WMO abbreviated header in Part A.
C                       Includes 62626 comment line.
C                       Adds coding of wind shear group.
C                       Adds tropopause coding.  This required changes to
C                       the way levels are stored prior to coding in Part
C                       A.
C 
C     01/30/98          Sets WS sig level criteria to 5 kt below 850 mb
C                       for NOAA encoding.
C
C     10/02/97          Uses 10 m wind in all surface reports.
C
C     07/29/97          Uses flt level T/H/W if mand level interpolation
C                       fails and flt level is within 0.5 mb of mand lvl
C                       in XXAA.  Upward extrapolations discontinued in
C                       XXBB if less than 0.5 mb.
C
C     07/28/97          Uses 70% for missing RH for height extrapolations.
C                       Adds checks on EXP(bad pressures).
C
C     04/25/97          Identifies superadiabatic layers (code = 16,17)
C                       (NOAA only.)
C
C     03/07/97          Identifies temp gap sig level (code = 15).
C
C                       Does not permit two sig levels at the same encoded
C                       pressure.
C
C     02/04/97          Establishes RHMIN as lower limit for allowable RH.
C                       Currently set for 2%.  LIMIT IS APPLIED ON 
C                       ENCODING ONLY, not on any calculations.
C
C     01/27/97          Check for * characters in message.
C
C     01/15/97          Changes way max wind at top is found.  Corrects PC
C                       roundoff for LVLHI?  Allows short extrapolations
C                       of mandatory heights in NOAA version.
C
C                       Adds MODE parameter to XXAA and XXBB calls to turn
C                       AFRES on and off without recompiling.
C
C     01/09/97          Adds 10166,67 for doubtful data.
C
C     12/24/96          Fix typo that would cause ht error if man level 
C                       fell between launch and first sonde pressure.
C
C     11/04/96          Section 9 (51515 group) added, for 10190 and 
C                       10191 indicators.  Extrapolated hts/prs only
C                       sent in part A for AFRES = .T.
C
C     10/02/96          Use 9999, not 999 for blank sig lvl.
C
C     05/24/96          Section 7 (31313 group) added.
C     ==================================================================== 
C 
C 
      PARAMETER (MXRC = 9000, NINVAR = 12, MAXTMPMSG = 500)
C 
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
C
      CHARACTER*1 ANS
      CHARACTER*65 COMMENT
      CHARACTER*66 LINE
      CHARACTER*4 STNS
      CHARACTER*3 MONTHS(12)
      CHARACTER*4 VERSION
	character*100  F1, F2, F3
	character*100 FNAME
	character*100 FILENAME
      LOGICAL WRCOMMENT, WARN, AUTOTEMP
      DATA MONTHS/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug',
     *            'Sep','Oct','Nov','Dec'/
      DATA WRCOMMENT/.FALSE./

C 
	F1 = ADJUSTL(FNAME)
	F2 = TRIM(F1)
	idot = SCAN(F2,'.')
	if (idot) then
		F3 = F2(1:idot-1)
	else
		F3 = F2
	endif
	iblank = SCAN(F3,' ')
	if (iblank) then
		FILENAME = F3(1:iblank-1)//'.tempdrop'
	else
		FILENAME = F3//'.tempdrop'
	endif
C 
c      OPEN(LUFO,FILE=FILENAME,STATUS='UNKNOWN',IOSTAT=IOS,ERR=900)
c      WRITE(LUFO,'("Sonde # ",I9.9,2X,F4.0," UTC",2X,I2,1X,A,1X,I2)')
c     +    ID,TIML/100,IDY,MONTHS(IMO),IYR
C
      WARN = .FALSE.
      DO 50 L = 1,IXEND
50       IF (SNDDAT(L,2).EQ.-999.) WARN = .TRUE.
      IF (WARN) THEN
         WRITE(LUT,'(//,1x,38("-"))')
         WRITE(LUT,'(1X,"*** WARNING: BAD PRESSURES IN FILE ***")')
         WRITE(LUT,'(1x,38("-"))')
         ENDIF
C
      if (nsndtype .lt. 3) then
		CALL XXAA(2) 
		call code_62626
	else
		call UUAA
	endif

      IF (.NOT. AUTOTEMP) THEN
         WRITE(LUT,'(/," Hit return to continue... ",$)')
         READ(*,'(A)') ANS
         ENDIF
	if (nsndtype .lt. 3) then
		CALL XXBB(2) 
		call code_62626
	else
		call UUBB
	endif

      IF (.NOT. AUTOTEMP) THEN
C         IF (LUG.NE.0) CALL PLOTSKEWT(LUG,1)
         WRITE(LUT,'(/," Hit return to continue... ",$)')
         READ(*,'(A)') ANS
         ENDIF
C
      IF (WRCOMMENT) THEN
        WRITE(LUT,'(/," Append comments to message (y/n)? ",$)')
        READ(*,'(A)') ANS 
        IF (ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
          WRITE(LUT,'(/," Enter comment (64 char max): ",$)')
          READ(*,'(A)') COMMENT
          ELSE
          COMMENT = ''
          ENDIF
        COMMENT(46:46) = ';'
        DO 100 L = 45,1,-1
          IF (COMMENT(L:L).EQ.' ') THEN
              COMMENT(L:L) = COMMENT(L+1:L+1)
              COMMENT(L+1:L+1) = ' '
              ELSE
              GOTO 200
              ENDIF
100       CONTINUE
200   continue  
c	WRITE(LUFO,'(A)') COMMENT
        ENDIF
C
C
c      CLOSE(LUFO) 
C
      WRITE(LUT,'(/,1x,65("-"))')
c      WRITE(LUT,'(" Message written to file ",A,":",/)') FILENAME
c      OPEN(LUFO,FILE=FILENAME,STATUS='UNKNOWN',IOSTAT=IOS,ERR=900)
250   continue
c	READ(LUFO,'(A)',END=910,ERR=900) LINE
c      WRITE(LUT,'(1X,A)') LINE
c      GOTO 250
C
910   WRITE(LUT,'(1x,65("-"))')
c      CLOSE(LUFO)
      RETURN
C 
900   WRITE(LUT,'(" *** ERROR IN SUBROUTINE TEMPDROP ***")')
      RETURN
      END 
C 
C 
C 
C     --------------------------------------------------------- 
      SUBROUTINE SAVEDATA
C 
C     Subroutine saves current sonde data to a holding buffer
C     Data can be retrieved with subroutine RESTORE
C     --------------------------------------------------------- 
C 
      PARAMETER (MXRC = 9000, NINVAR = 12, NOPT = 16, NPROC = 50)
      CHARACTER*4 STNS,STNSH
      LOGICAL NEEDOPT,NEEDOPTH
C 
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /OPTIONS/ NEEDOPT(NOPT)
      COMMON /PROCESSING/ NPSTP,PROC(NPROC)
C
      COMMON /HOLDSND/ IDH,IXENDH,SNDDATH(MXRC,NINVAR),
     *                 ISNDFLGH(MXRC,4)
      COMMON /HOLDLCH/ IYRH,IMOH,IDYH,RLATH,RLONH,TIMLH,PRLH,TELH,
     *                 RHLH,DPLH,HTLH,WDLH,WSLH,STNSH
      COMMON /HOLDOPT/ NEEDOPTH(NOPT)
      COMMON /HOLDPROC/ PROCH(NPROC)
C 
C 
      DO 100 I = 1,MXRC
         DO 110 J = 1,NINVAR
            SNDDATH(I,J) = SNDDAT(I,J)
110         CONTINUE
         DO 120 J = 1,4
            ISNDFLGH(I,J) = ISNDFLG(I,J)
120         CONTINUE
100      CONTINUE
C
      DO 150 I = 1,NOPT
150      NEEDOPTH(I) = NEEDOPT(I)
C
      DO 160 I = 1,NPROC
160      PROCH(I) = PROC(I)
C
      IDH = ID
      IXENDH = IXEND
      IYRH = IYR
      IMOH = IMO
      IDYH = IDY
      RLATH = RLAT
      RLONH = RLON
      TIMLH = TIML
      PRLH = PRL
      TELH = TEL
      RHLH = RHL
      DPLH = DPL
      HTLH = HTL
      WDLH = WDL
      WSLH = WSL
      STNSH = STNS
C
      RETURN
      END 
C
C
C
C     --------------------------------------------------------- 
      SUBROUTINE RESTORE(KOP)
C 
C     Subroutine restores sonde data from holding buffer
C     --------------------------------------------------------- 
C 
      PARAMETER (MXRC = 9000, NINVAR = 12, NOPT = 16, NPROC = 50)
      CHARACTER*4 STNS,STNSH
C 
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /OPTIONS/ NEEDOPT(NOPT)
      COMMON /PROCESSING/ NPSTP,PROC(NPROC)
C
      COMMON /HOLDSND/ IDH,IXENDH,SNDDATH(MXRC,NINVAR),
     *                 ISNDFLGH(MXRC,4)
      COMMON /HOLDLCH/ IYRH,IMOH,IDYH,RLATH,RLONH,TIMLH,PRLH,TELH,
     *                 RHLH,DPLH,HTLH,WDLH,WSLH,STNSH 
      COMMON /HOLDOPT/ NEEDOPTH(NOPT)
      COMMON /HOLDPROC/ PROCH(NPROC)
C
C 
      DO 100 I = 1,MXRC
         DO 110 J = 1,NINVAR
            SNDDAT(I,J) = SNDDATH(I,J)
110         CONTINUE
         DO 120 J = 1,4
            ISNDFLG(I,J) = ISNDFLGH(I,J)
120         CONTINUE
100      CONTINUE
      DO 150 I = 1,NOPT
150      NEEDOPT(I) = NEEDOPTH(I)
C
      DO 160 I = 1,NPROC
160      PROC(I) = PROCH(I)
C
C
      ID = IDH
      IXEND = IXENDH
      IYR = IYRH
      IMO = IMOH
      IDY = IDYH
      RLAT = RLATH
      RLON = RLONH
      TIML = TIMLH
      PRL = PRLH
      TEL = TELH
      RHL = RHLH
      DPL = DPLH
      HTL = HTLH
      WDL = WDLH
      WSL = WSLH
      STNS = STNSH
C
      WRITE(LUT,'(/," Last call to option ",i2.2," has been negated.")')
     *      KOP
      WRITE(LUT,'(1X,"Suggested options may now be inappropriate."/)')
      KOP = 0
      RETURN
      END 
C
C
C
C     -----------------------------------------------
      SUBROUTINE SAVEFILE(NTYPE,NRECS)
	
C
C     Saves working file of sonde data
C     -----------------------------------------------
C 
      PARAMETER (MXRC = 9000, NOPT = 16, NINVAR = 12) 
      PARAMETER (NEOD = 100) 
      PARAMETER (NPROC = 50)
C 
      CHARACTER*4 STNS
      character*100  FILENAME
      CHARACTER*9 DIREC
      CHARACTER*4 FTYPE, VERSION
      character*100  COMMENT
      CHARACTER*32 PROCTIME
      LOGICAL NEEDOPT, AUTOTEMP
C 
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
      COMMON /OPTIONS/ NEEDOPT(NOPT)
      COMMON /EODDATA/ EOD(NEOD,NINVAR)
      COMMON /PROCESSING/ NPSTP,PROC(NPROC)
C
C
      IF (NTYPE.EQ.1) THEN
         FTYPE = '.tmp'
         DIREC = './'
         IF (NSITE.EQ.2) DIREC = 'sondedata'
         ENDIF
      IF (NTYPE.EQ.2) THEN
         FTYPE = '.pds'
         DIREC = 'AVAPS_fnl'
         ENDIF
C
      IF (NSNDTYPE.EQ.1) 
     *   WRITE(FILENAME,'(A9,"/odw",I5.5,A4)') DIREC,ID,FTYPE
      IF (NSNDTYPE.GT.1) 
     *   WRITE(FILENAME,'(A9,"/g",I9.9,A4)') DIREC,ID,FTYPE
      OPEN(LUFW,FILE=FILENAME,FORM='UNFORMATTED',STATUS='UNKNOWN',
     *     ERR=900)
C
C     Write extra lines for processed files
C     -------------------------------------
      NR = MXRC
      IF (FTYPE.EQ.'.pds') THEN
         WRITE(LUFW,ERR=900) VERSION
         WRITE(LUFW,ERR=900) NRECS
         WRITE(LUT,'(" Add comment (A60): ",$)')
         READ(*,'(A)') COMMENT
         WRITE(LUFW,ERR=900) COMMENT
         NR = NRECS
C         CALL CFTIME(PROCTIME)
         WRITE(LUFW,ERR=900) PROCTIME(6:17)
         ENDIF
C
C     Write rest of file
C     ------------------
      WRITE(LUFW,ERR=900) LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *     NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      WRITE(LUFW,ERR=900) ID,IXEND,
     *     ((SNDDAT(I,J),I=1,NR),J=1,NINVAR),
     *     ((ISNDFLG(I,J),I=1,NR),J=1,4)
      WRITE(LUFW,ERR=900) IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,
     *     DPL,HTL,WDL,WSL,STNS 
      WRITE(LUFW,ERR=900) NEEDOPT,EOD,(X(I),I=1,NR),
     *     (Y(I),I=1,NR),(Z(I),I=1,NR)     
      WRITE(LUFW,ERR=900) NPSTP,(PROC(I),I=1,NPROC)
      CLOSE(LUFW)
      RETURN
C
C     Error on file
C     -------------
900   WRITE(LUT,'(/," *** ERROR WRITING FILE ",A,/)') FILENAME 
      RETURN
      END
C
C
C     -----------------------------------------------
      SUBROUTINE RETRIEVE(NTYPE,NRECS)
	
C
C     Retrieves last working file of sonde data
C     -----------------------------------------------
C 
      character*100  FILENAME
C 
      CALL GETFNAME(NTYPE,FILENAME)
      CALL OPENSONDE(FILENAME,NTYPE,NRECS)
      RETURN
      END
C
C
C
C     -------------------------------------------------------------
      SUBROUTINE SETLU(LUG)
C     -------------------------------------------------------------
C     Strictly speaking, this should be set to the graphics lu.
C     For plotting to the terminal however, calls to INITL in
C     hpslib will override LUG and set it to 88.  So leave this as
C     50 and the code will work for both laser printer AND terminal
C     -------------------------------------------------------------
      LUG = 50
      RETURN
      END
C
C
C
c     -------------------------------------------------------------
      SUBROUTINE DROPDECODE(LUT)
c     
c     Decodes accumulated TEMPDROP messages and appends current
c     sonde's message into spline analysis input format.  
c     -------------------------------------------------------------
c
      character*70 line
      character*30 filein,fileout
      logical current, no_prev_messages
      common /output/nrecs
c
c     -----------------------------------------------------
c     open unit 11 for control data
c     -----------------------------------------------------
      open (11,file='synmap.ctl',iostat=istat,status='old',err=99) 
      no_prev_messages = .false.
c
c     read in the beginning year and month 
c     ------------------------------------
      read (11,'(3i2)') iyrs,imns,idys
      read (11,'(a30)') filein
      read (11,'(a30)') fileout
      close (11)
c
c     -----------------------------------------------------
c     open unit 12 for tempdrop input and 13 for SAP output
c     -----------------------------------------------------
      open (12,file=filein,iostat=istat,status='old',err=95)
5     open (13,file=fileout,iostat=istat,status='unknown',
     *      err=99)
c
c
      iwx = 1
      iflag = 1
      nrecs = 0
      current = .false.
      if (no_prev_messages) goto 55
c
c     --------------------------------------------------------
c     read the first line of a message.  This version does not
c     require UZNT header, just the XXAA
c     --------------------------------------------------------
10    read(12,'(a)',end=50)line
      if(index(line,'XXAA').ne.0) 
     *   call drop(iwx,iflag,iyrs,imns,idys,line)
c
c
c     start checking the next message
c     -------------------------------
      goto 10
c
c     ----------------------------------------------------------
c     end of accumulated messages has been reached.  For current
c     sonde we change data type to 50 to distinguish it from
c     others.  Also a flag to decoder that we are at the XXAA 
c     line
c     ----------------------------------------------------------
 50   close(12)
 55   if (current) goto 90
      current = .true.
      open (12,file='tempdrop.dat',iostat=istat,status='old',err=99)
      iwx = 50
      goto 10
c
c
90    close(13)
      return
c
c
c     Can't open file of previous messages
c     ------------------------------------
 95   write(lut,'(/" CANNOT OPEN FILE OF TRANSMITTED MESSAGES...")')
      write(lut,'(" Current sonde assumed to be the first one..."/)')
      no_prev_messages = .true.
      goto 5
c
 99   write(lut,'(" CANNOT OPEN tempdrop.dat FILE ***")')
      close(13)
      return
      end
C
C
C
C     ------------------------------------------------------
C     avapslib.f
C
C     A collection of routines relating to AVAPS dropsondes,
C     including file I/O and sonde characteristics.
C     ------------------------------------------------------
C
C
C
C     ----------------------------------------------------
      SUBROUTINE READGSONDE(NRECS,RUNSTRING,FILENAME)
	
C 
C     Opens file of G-SONDE data and loads up data
C     arrays.  Initializes arrays first.  NRECS is number
C     of sounding records read.
C     ----------------------------------------------------
C 
	use thread_common
	  
	PARAMETER (MXRC = 9000, NOPT = 16, NINVAR = 12)
      PARAMETER (NEOD = 100) 
      PARAMETER (NPROC = 50)
C 
      DIMENSION TIMGAP(MXRC)
      CHARACTER*4 STNS
      CHARACTER*9 RUNSTRING, DIREC
      character*100  FNAME 
	character*100 FILENAME
      CHARACTER*122 LINE
      CHARACTER*3 MTYP
      CHARACTER*4 VERSION
      LOGICAL AUTOTEMP
      LOGICAL SREC,LREC,AREC,ARECW,NEEDOPT,LPREM,TGAP
      LOGICAL CRCTOSS(3)
C
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
      COMMON /EODDATA/ EOD(NEOD,NINVAR)
      COMMON /OPTIONS/ NEEDOPT(NOPT)
      COMMON /PROCESSING/ NPSTP,PROC(NPROC)
C
      DATA BAD/-999./
	DATA EPS/0.01/

	real(4) GlassTime(MXRC)
	real(4) GlassPres(MXRC)
	real(4) GlassTdry(MXRC)
	real(4) GlassRH(MXRC)
	real(4) GlassWspd(MXRC)
	real(4) GlassWdir(MXRC)
	real(4) GlassAlt(MXRC)
	real(4) GlassLat(MXRC)
	real(4) GlassLon(MXRC)
	real(4) GlassDz(MXRC)
	real(4) GlassSats(MXRC)

	integer(4) GlassDataPoints
	
      INTERFACE
		integer FUNCTION getsounding(filename)
			!MS$ATTRIBUTES STDCALL:: getsounding
			!MS$ATTRIBUTES REFERENCE:: filename
			character*100 filename
		END FUNCTION getsounding

		SUBROUTINE delsounding()
			!MS$ATTRIBUTES STDCALL:: delsounding
		END SUBROUTINE delsounding

		integer FUNCTION datapoints()
			!MS$ATTRIBUTES STDCALL:: datapoints
		END FUNCTION datapoints

		SUBROUTINE sndtime(p)
			!MS$ATTRIBUTES STDCALL:: sndtime
			!MS$ATTRIBUTES REFERENCE:: p
			real(4) p(1)
		END SUBROUTINE sndtime

		SUBROUTINE sndpres(p)
			!MS$ATTRIBUTES STDCALL:: sndpres
			!MS$ATTRIBUTES REFERENCE:: p
			real(4) p(1)
		END SUBROUTINE sndpres

		SUBROUTINE sndtdry(p)
			!MS$ATTRIBUTES STDCALL:: sndtdry
			!MS$ATTRIBUTES REFERENCE:: p
			real(4) p(1)
		END SUBROUTINE sndtdry

		SUBROUTINE sndrh(p)
			!MS$ATTRIBUTES STDCALL:: sndrh
			!MS$ATTRIBUTES REFERENCE:: p
			real(4) p(1)
		END SUBROUTINE sndrh

		SUBROUTINE sndwspd(p)
			!MS$ATTRIBUTES STDCALL:: sndwspd
			!MS$ATTRIBUTES REFERENCE:: p
			real(4) p(1)
		END SUBROUTINE sndwspd

		SUBROUTINE sndwdir(p)
			!MS$ATTRIBUTES STDCALL:: sndwdir
			!MS$ATTRIBUTES REFERENCE:: p
			real(4) p(1)
		END SUBROUTINE sndwdir
		
		SUBROUTINE sndlat(p)
			!MS$ATTRIBUTES STDCALL:: sndlat
			!MS$ATTRIBUTES REFERENCE:: p
			real(4) p(1)
		END SUBROUTINE sndlat

		SUBROUTINE sndlon(p)
			!MS$ATTRIBUTES STDCALL:: sndlon
			!MS$ATTRIBUTES REFERENCE:: p
			real(4) p(1)
		END SUBROUTINE sndlon

		SUBROUTINE sndalt(p)
			!MS$ATTRIBUTES STDCALL:: sndalt
			!MS$ATTRIBUTES REFERENCE:: p
			real(4) p(1)
		END SUBROUTINE sndalt

		SUBROUTINE sndsats(p)
			!MS$ATTRIBUTES STDCALL:: sndsats
			!MS$ATTRIBUTES REFERENCE:: p
			real(4) p(1)
		END SUBROUTINE sndsats

		SUBROUTINE snddz(p)
			!MS$ATTRIBUTES STDCALL:: snddz
			!MS$ATTRIBUTES REFERENCE:: p
			real(4) p(1)
		END SUBROUTINE snddz

		INTEGER FUNCTION sndyear()
			!MS$ATTRIBUTES STDCALL:: sndyear
		END FUNCTION sndyear

		INTEGER FUNCTION sndmonth()
			!MS$ATTRIBUTES STDCALL:: sndmonth
		END FUNCTION sndmonth

		INTEGER FUNCTION sndday()
			!MS$ATTRIBUTES STDCALL:: sndday
		END FUNCTION sndday

		INTEGER FUNCTION sndhhmmss()
			!MS$ATTRIBUTES STDCALL:: sndhhmmss
		END FUNCTION sndhhmmss


	END INTERFACE


C 
C 
C     Do we keep or toss telemetry bit errors for PTH?
C     ------------------------------------------------
      IF (NSITE.EQ.1) THEN
         CRCTOSS(1) = .FALSE.  		
         CRCTOSS(2) = .FALSE.		
         CRCTOSS(3) = .TRUE.		
         ELSE
         CRCTOSS(1) = .TRUE.  		
         CRCTOSS(2) = .TRUE.		
         CRCTOSS(3) = .TRUE.		
         ENDIF
C
C
C     Initialize other variables 
C     -------------------------- 
      SREC = .FALSE.
      LREC = .FALSE.
      AREC = .FALSE.
      ARECW = .FALSE.
      LPREM = .FALSE.
C
      ifmt = 1
	if (fileType_thr .EQ. 1) then
		ifmt = 0
	endif
	if (filetype_thr .EQ. 2) then
		ifmt = 1
	endif

      L = 0
      LF = 0
      IXEND = 0 
      IXEOD = 0
      NTG = 0
      CURR = 0.
C
      PRLAST = BAD
      PRL = BAD
      TEL = BAD
      RHL = BAD
      DPL = BAD
      HTL = BAD
      WDL = BAD
      WSL = BAD
      RLAT = BAD
      RLON = BAD
      STNS = 'GPS '
C
      DO 10 I = 1, MXRC 
         DO 15 J = 1,NINVAR 
15            SNDDAT(I,J) = BAD 
         DO 20 J = 1,4
              ISNDFLG(I,J) = 0
20            CONTINUE
10       CONTINUE 
      DO 25 I = 1, NEOD 
         DO 30 J = 1,NINVAR
              EOD(I,J) = BAD
30            CONTINUE
25       CONTINUE 
      DO 35 I = 1, NPSTP
35       PROC(I) = 0.
      DO 40 I = NPSTP+1,NPROC
40       PROC(I) = BAD
C 
C
C     Calculate how often to keep data (every NSKIPth point)
C     ------------------------------------------------------
      FILERATE = 2.0                ! Data are 2/sec
	if (nsndtype .EQ.3) FILERATE = 1.0
	if (nsndtype .EQ.4) FILERATE = 0.1
      SKIP = FILERATE/DATARATE
      NSKIP = NINT(FILERATE/DATARATE)
      TIMRES = 1./FILERATE
      IF (FLOAT(NSKIP).NE.SKIP) THEN
         WRITE(LUT,'(" *** FATAL ERROR: INVALID DATARATE ***")')
c		 stop
		 return
         ENDIF

	if (nsndtype .GE. 3) then

		istatus = getsounding(filename)

		if (istatus) then
			flag_thr = -1
			status_msg_thr = 'illegal GLASS file'
		else
			status_msg_thr = 'GLASS file successfully read'
			GlassDataPoints = datapoints()
			if (GlassDataPoints .GE. MXRC) then
				status_msg_thr = 'GLASS data file to large'
				flag_thr = -1
			else
				call sndtime(GlassTime)
				call sndpres(GlassPres)
				call sndtdry(GlassTdry)
				call sndrh(GlassRH)
				call sndwspd(GlassWspd)
				call sndwdir(GlassWdir)
				call sndalt(GlassAlt)
				call sndlon(GlassLon)
				call sndlat(GlassLat)
				call snddz(GlassDz)
				call sndsats(GlassSats)

				call sndTransfer(GlassTime, GlassDataPoints, 1,  BAD)
				call sndTransfer(GlassPres, GlassDataPoints, 2,  BAD)
				call sndTransfer(GlassTdry, GlassDataPoints, 3,  BAD)
				call sndTransfer(GlassRH,   GlassDataPoints, 4,  BAD)
				call sndTransfer(GlassWdir, GlassDataPoints, 6,  BAD)
				call sndTransfer(GlassWspd, GlassDataPoints, 7,  BAD)
				call sndTransfer(GlassSats, GlassDataPoints, 8,  BAD)
				call sndTransfer(GlassAlt,  GlassDataPoints, 9,  BAD)
				call sndTransfer(GlassDz,   GlassDataPoints, 10, BAD)
				call sndTransfer(GlassLat,  GlassDataPoints, 11, BAD)
				call sndTransfer(GlassLon,  GlassDataPoints, 12, BAD)

				RLAT = GlassLat(1)
				RLON = -GlassLon(1) 
C				NOTICE THAT RLON IS POSITIVE WEST!!!
				
				IYR = sndyear() - 1900
				IMO = sndmonth()  
				IDY = sndday()
				TIML = sndhhmmss()

				IXEND = GlassDataPoints

c			make the launch time consistent
c			with the sampled data (subtract 1 second)
			snddat(IXEND, 1) = snddat(IXEND-1,1) - 1

c			delete the sounding so that the file gets closed.
			call delsounding()

			endif
		endif

		return

	endif
C
C
C     Get sonde id, open input file
C     -----------------------------
c    IF (RUNSTRING.NE.'         ') GOTO 55
c50    WRITE(LUT,'(/," Enter sonde id #: ",$)')
c    READ(*,'(a)',ERR=50) RUNSTRING
c55    READ(RUNSTRING,*,ERR=50) ID
C
      DIREC = '.\'
      IF (NSITE.EQ.2) DIREC = 'sondedata'
c      WRITE(FILENAME,'(A9,"/g",I9.9,".avp")') DIREC,ID
      OPEN(LUFI,FILE=filename,IOSTAT=IERR,STATUS='OLD',ERR=900) 
c	INQUIRE(LUFI,NAME=FILENAME)
C
C 
C     Echo text lines, extract launch record
C     --------------------------------------
      WRITE(LUT,90) ID
90    FORMAT(/,1x,'HEADER AND COMMENTS FOR SONDE ID ',
     *       I9.9,':',/,1x,78("-"))
C
C
C     Read next line in file
C     ----------------------
100   READ(LUFI,'(A)',END=800) LINE
C
      IF (LINE(11:13) .EQ. 'END') GOTO 800
C
      IF (LINE(11:13) .EQ. 'VER') THEN
         IFMT = 1
         ENDIF
C
      IF (LINE(11:13) .EQ. 'STA') THEN
         READ(LINE,195) IDL, IYR, IMO, IDY, TIMS
	   ID = IDL
         IF (IDL.NE.ID) GOTO 910
         GOTO 100
         ENDIF
C
      IF (LINE(11:13) .EQ. 'LAU') THEN
         READ(LINE,195,ERR=100) IDL, IYR, IMO, IDY, TIML
195      FORMAT(14X,I9,1X,3I2,1X,F8.1)
	   ID = IDL
	   IF (IDL.NE.ID) GOTO 910
         LREC = .TRUE.
         IF (HMSTS(TIML).LT.HMSTS(TIMS)) LPREM = .TRUE.
         ENDIF
C
      IF (LINE(1:7) .EQ. 'AVAPS-T') THEN
         IF (LINE(11:30) .EQ. 'COM --------- ------' .OR.
     *       LINE(11:30) .EQ. 'COM    ID     yymmdd' .OR.
     *       LINE(11:30) .EQ. 'COM   Sonde    Date ' .OR.
     *       LINE(11:30) .EQ. 'COM            GMT  ' .OR.
     *       LINE(11:30) .EQ. 'COM                 ') GOTO 100
         WRITE(LUT,'(1x,A)') LINE(11:88)
         DO 105 J=89,122
            IF (LINE(J:J).NE.' ') THEN
               WRITE(LUT,'(40X,A)') LINE(J:122)
               GOTO 100
               ENDIF
105         CONTINUE
         GOTO 100
         ENDIF
C
C
C     Line has data - decode and process
C     ----------------------------------
      IF (IFMT.EQ.0) THEN
         READ(LINE,190,ERR=100) 
     *   ICH,MTYP,ID,IYRX,IMOX,IDYX,TIME,
     *   PR,TE,RH,WD,WS,VV,SLON,SLAT,GA
         NSATS = BAD
         RH1 = BAD
         RH2 = BAD
         ENDIF
C
      IF (IFMT.EQ.1) READ(LINE,191,ERR=100) 
     *   ICH,MTYP,ID,IYRX,IMOX,IDYX,TIME,
     *   PR,TE,RH,WD,WS,VV,SLON,SLAT,GA,
     *   NSATS,RH1,RH2
C
190   FORMAT(7X,I2,1X,A3,1X,I9,1X,3I2,1X,F8.1,1X,
     *       F7.2,1X,5(F6.2,1X),F10.5,1X,F9.5,1X,F8.2)
191   FORMAT(7X,I2,1X,A3,1X,I9,1X,3I2,1X,F8.1,1X,
     *       F7.2,1X,5(F6.2,1X),F10.5,1X,F9.5,1X,F8.2,
     *       1X,I4,1X,F6.2,1X,F6.2)
C
C
C     Check and skip pre-launch records
C     Use pre-launch location as first guess
C     --------------------------------------
      IF (MTYP(1:1).EQ.'P') THEN
         PLAT = SLAT
         PLON = -SLON
         GOTO 100
         ENDIF
C
C     Is this the flight level data?
C     ------------------------------
      IF (MTYP(1:1).EQ.'A') THEN
         AREC = .TRUE.
         TIMA = TIME
         IF (LREC .AND. ABS(HMSTS(TIMA)-HMSTS(TIML)).GT.20.) THEN
            ARECW = .TRUE.
            NEEDOPT(6) = .TRUE.
            ENDIF
         RLAT = SLAT
         RLON = -SLON         
         PRL = PR
         TEL = TE
         RHL = RH
         DPL = DEWPT(TEL,RHL)         
         HTL = GA
         WDL = WD
         WSL = WS
         STNS = 'GPS '
         PRLAST = PRL
         GOTO 100
         ENDIF
C
C
C     Everything that gets to this point is sounding data.
C     Check whether we got launch and aircraft data yet.
C     If not, get these parameters from the first 
C     sounding record or prelaunch.
C     ----------------------------------------------------
      IF (SREC) GOTO 125
C
      IF (.NOT.LREC) THEN
         TIML = TIME      
         IYR = IYRX
         IMO = IMOX
         IDY = IDYX
         ENDIF
C
      IF (.NOT.AREC) THEN
         RLAT = PLAT
         RLON = PLON
         ENDIF
C
C
C     Process sounding data
C     ---------------------
125   SREC = .TRUE.
      TGAP = .FALSE.
      LF = LF+1
      L = L+1
C
C
C     Check if file too long
C     ----------------------
      IF (L.GT.MXRC-1) THEN            
         WRITE(LUT,'(/," *** WARNING: INPUT FILE TOO LONG ***"/)') 
         GOTO 800
         ENDIF
C
C
C     Get time of current file record
C     -------------------------------
      IF (IDYX.NE.IDY) TIME = TIME+240000.
      SNDDAT(L,1) = HMSTS(TIME)-HMSTS(TIML)
      PREV = CURR
      CURR = SNDDAT(L,1)
c	write (*,'(1x, "prev and curr are: ", f5.1,2x,f5.1)')prev,curr
C
C
C     Check and fill time gaps
C     Gap at start.
C     ------------------------
      IF (LF.EQ.1) THEN
         IF (CURR.LE.TIMRES) GOTO 140
         REM = AMOD(CURR,TIMRES)
         PREV = REM-TIMRES
C
130      NTG = NTG+1
         TIMGAP(NTG) = PREV+TIMRES
         SNDDAT(L,1) = PREV+TIMRES
         PREV = SNDDAT(L,1)
         DO 131 J = 2,NINVAR
            SNDDAT(L,J) = BAD
131         CONTINUE
         L = L+1
         IF ((CURR-PREV-TIMRES).GT.EPS) GOTO 130
         SNDDAT(L,1) = CURR
         GOTO 140
         ENDIF
C
C     Gap in middle
C     -------------
C	write(*,'(2x,f5.1,2x,f5.1,2x,f5.1,2x,f5.1)')curr,prev,timres,eps
      IF ((CURR-PREV-TIMRES).GT.EPS) THEN
         NTG = NTG+1
         TIMGAP(NTG) = PREV+TIMRES
         SNDDAT(L,1) = PREV+TIMRES
         CURR = SNDDAT(L,1)
         DO 135 J = 2,NINVAR
            SNDDAT(L,J) = BAD
135         CONTINUE
         TGAP = .TRUE.
         GOTO 145
         ENDIF
C
C
C     Record exists, assign met variables
C     -----------------------------------
140   IF (PR.LT.9998.)  SNDDAT(L,2) = PR
      IF (TE.LT.98.)    SNDDAT(L,3) = TE
      IF (RH.LT.998.)   SNDDAT(L,4) = RH
C
C
C     Check for telemetry bit errors and toss if desired
C     --------------------------------------------------
      IF (MTYP(2:2).EQ.'1') THEN
         DO 141 I = 1,3
            IF (CRCTOSS(I)) SNDDAT(L,I+1) = BAD
141         CONTINUE
         ENDIF
C
C     Calculate estimated pressure
C     ----------------------------
      IF (PRLAST.NE.BAD) THEN
         ESTPRI = 2.0*(0.3 + (PRLAST*.000410))/DATARATE
         IF (SNDDAT(L,1).LE.0.) ESTPRI = 0.0
         IF (SNDDAT(L,1).LE.2.) ESTPRI = ESTPRI/2.0
         SNDDAT(L,5) = PRLAST+ESTPRI
         PRLAST = SNDDAT(L,5)
         ENDIF
C
      IF (WD.LT.998. .AND. WS.LT.998.) THEN
                        SNDDAT(L,6) = WD
                        SNDDAT(L,7) = WS
                        SNDDAT(L,8) = FLOAT(NSATS)
                        ENDIF
      IF (GA.LT.99998.) SNDDAT(L,9) = GA
      IF (VV.LT.98.)    SNDDAT(L,10) = VV
      IF (SLAT.LT.98.)  SNDDAT(L,11) = SLAT
      IF (SLON.LT.998.) SNDDAT(L,12) = -SLON
C
C
C     Fill EOD buffer if DATARATE less than full resolution
C     -----------------------------------------------------
145   IF (DATARATE.GE.2.0) GOTO 200
      DO 150 I = 1,NEOD-1
         DO 155 J = 1,NINVAR
155         EOD(I,J) = EOD(I+1,J)
150      CONTINUE
      DO 160 J = 1,NINVAR
160      EOD(NEOD,J) = SNDDAT(L,J)
C         
C
C     Skip this one?
C     --------------
200   IF (MOD(LF-1,NSKIP).NE.0) THEN
         DO 210 J = 1,NINVAR
            SNDDAT(L,J) = BAD
210         CONTINUE
         L = L-1
         GOTO 300
         ENDIF
C
C
C     Return for next line...
C     But don't read new line till gap over
C     -------------------------------------
300   IF (TGAP) GOTO 125        
      GOTO 100
C
C
C     End of data, normal exit
C     ------------------------
800   IXEND = L
      if (IXEND .gt. MXRC) then
		IXEND = MXRC
	endif

      IF (IXEND.EQ.0) GOTO 920
      
	NRECS = IXEND

C
C     Check flight level values
C     -------------------------
      IF (PRL.GT.9998.) PRL = BAD
      IF (TEL.GT.98.) TEL = BAD
      IF (RHL.GT.998.) THEN
         RHL = BAD
         DPL = BAD
         ENDIF
      IF (HTL.GT.99998.) HTL = BAD
      IF (WDL.GT.998.) WDL = BAD
      IF (WSL.GT.998.) WSL = BAD
      IF (RLAT.GT.98.) RLAT = BAD
      IF (RLON.GT.998.) RLON = BAD
C
      WRITE(LUT,'(1x,78("-"),/)')
C
C
C     Check date, if earlier than 28 Sep 1996, then fall rate
C     was a hydrostatic fall rate...discard.
C     -------------------------------------------------------
      DATE = IDY+100.*IMO+10000.*IYR
      IF (DATE.LT.960928.) THEN
         DO 810 J = 1,IXEND
            SNDDAT(J,10) = BAD
810         CONTINUE
         ENDIF
C
C
C     Display file warning messages
C     -----------------------------
      IF (NTG.GT.0) THEN
         WRITE(LUT,*)
         DO 850 J = 1,NTG
c            WRITE(LUT,'(" *** RECORD BAD OR MISSING AT ", 
c     *      F6.1," S ***")') TIMGAP(J)
            status_msg_thr = " *** RECORD BAD OR MISSING "
	
850         CONTINUE
         ENDIF
      IF (.NOT.AREC) 
	*	status_msg_thr = " *** WARNING: NO AIRCRAFT DATA IN FILE ***"
c     *   WRITE(LUT,'(/" *** WARNING: NO AIRCRAFT DATA IN FILE ***")')
      IF (ARECW)
	*	status_msg_thr = 
	*		" *** WARNING: AIRCRAFT DATA NOT CLOSE TO LAUNCH ***" 
c		WRITE(LUT,
c     *   '(/" *** WARNING: AIRCRAFT DATA NOT CLOSE TO LAUNCH ***")')
      IF (LPREM)
	*	status_msg_thr =
	*		" *** WARNING: LAUNCH TIME EARLIER THAN START TIME ***"
c	WRITE(LUT,
c     *   '(/" *** WARNING: LAUNCH TIME EARLIER THAN START TIME ***")')
      IF (.NOT.LREC)
	*	status_msg_thr = 
	*		" *** WARNING: LAUNCH TIME EARLIER THAN START TIME ***"
c     *   WRITE(LUT,'(/" *** WARNING: NO LAUNCH RECORD IN FILE ***")')
	close(lufi)
      RETURN
C
C
C
C     Error stops
C     -----------
900   continue
c	WRITE(LUT,'(/," *** FATAL ERROR OPENING INPUT FILE ",A,/)') 
c    *      FILENAME 
c      STOP
	status_msg_thr = " *** FATAL ERROR OPENING INPUT FILE "
	flag_thr = -1
	close(lufi)
	return
C
910   continue
c	WRITE(LUT,'(/," *** FATAL ERROR: BAD SONDE ID (",I9,")",A,/)') 
c     *      IDL,FILENAME 
c      STOP
	status_msg_thr = " *** FATAL ERROR: BAD SONDE ID" 
	flag_thr = -1
	close(lufi)
	return
C
920   continue
c	WRITE(LUT,'(/," *** FATAL INPUT ERROR: NO VALID TIMES",/)') 
c      STOP
	status_msg_thr = " *** FATAL INPUT ERROR: NO VALID TIMES"
	flag_thr = -1
	close(lufi)
	return
C
      END 
C
C
C
C     ----------------------------------------------------
      SUBROUTINE READODW(NRECS,RUNSTRING)
	
C 
C     Opens file of raw ODW data and loads up ODW data
C     arrays.  Initializes arrays first.  NRECS is number
C     of 10-s records read.
C     ----------------------------------------------------
C 
      PARAMETER (MXRC = 9000, NINVAR = 12)
      PARAMETER (NEOD = 100, NEODODW = 22) 
      PARAMETER (NPROC = 50)
C 
      CHARACTER*4 STNS
      CHARACTER*50 RUNSTRING, DIREC
      character*100  FILENAME 
      CHARACTER*4 VERSION
      LOGICAL AUTOTEMP
C
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
      COMMON /EODDATA/ EOD(NEOD,NINVAR)
      COMMON /PROCESSING/ NPSTP,PROC(NPROC)
C 
C 
      IXEND = 0 
      IXEOD = 0
C
      IF (RUNSTRING.NE.'         ') GOTO 15
10    WRITE(LUT,'(/,"Enter sonde id #: ",$)')
      READ(*,'(a)',ERR=10) RUNSTRING
15    READ(RUNSTRING,*,ERR=10) ID
C
      DIREC = './'
      IF (NSITE.EQ.2) DIREC = 'sondedata'
      WRITE(FILENAME,'(A9,"/odw",I5.5,".dat")') DIREC,ID
      OPEN(LUFI,FILE=FILENAME,IOSTAT=IERR,STATUS='OLD',ERR=900) 
C 
C     Initialize all the input ODW data array variables 
C     ------------------------------------------------- 
      DO 50 I = 1, MXRC 
         DO 60 J = 1,NINVAR 
60            SNDDAT(I,J) = -999. 
         DO 70 J = 1,4
              ISNDFLG(I,J) = 0
70            CONTINUE
50       CONTINUE 
      DO 80 I = 1, NEOD 
         DO 85 J = 1,NINVAR
              EOD(I,J) = -999.
85            CONTINUE
80       CONTINUE 
      DO 90 I = 1, NPSTP
90       PROC(I) = 0.
      DO 95 I = NPSTP+1,NPROC
95       PROC(I) = -999.
C 
C     ------------------------------------------------------------
C     Program currently set up to read ODW files reformatted
C     on HP workstation (program sonde_input).
C     ------------------------------------------------------------
C 
C     Read launch parameters
C     ----------------------
      READ(LUFI,9000) IDS,IYR,IMO,IDY,TIML,RLAT,RLON,IPA,TEL,
     *                WDL,WSL,STNS
9000  FORMAT(I5,4X,3I2,1X,I6,1X,2F8.3,I5,F6.1,2I5,1X,A)
C 
C     Store permanent launch parameters (N,W are positive)
C     ----------------------------------------------------
      PRL = IPA 
      RHL = -999. 
      DPL = -999. 
      HTL = -999. 
C 
C 
C     Read PTH and wind records
C     -------------------------
100   READ(LUFI,9010) TIME,PRESS,TEMP,RELH,ESTPR,WDIR,WSPD,WQAL
9010  FORMAT(2I5,1X,F5.1,5I5)
      IF (TIME.GT.9990) GOTO 120
C 
C     Assign PTH and wind data to permanent variables
C     -----------------------------------------------
      IX = INT(TIME/10.+0.5)
      IF (IX.LT.1 .OR. IX.GT.MXRC) THEN            !Check if too big    
         WRITE(LUT,'(/,"*** INVALID INPUT TIME: ",F6.0,/)') TIME
      ELSE 
         SNDDAT(IX,1) = TIME
         SNDDAT(IX,2) = PRESS 
         SNDDAT(IX,3) = TEMP
         SNDDAT(IX,4) = RELH
         SNDDAT(IX,5) = ESTPR         ! This position is for est. pr.      
         SNDDAT(IX,6) = WDIR
         SNDDAT(IX,7) = WSPD
         SNDDAT(IX,8) = WQAL
         IF (IX.GT.IXEND) IXEND = IX
         ENDIF
      GOTO 100
C 
C     Read end-of-drop records 
C     ------------------------
120   READ(LUFI,9010) TIME,PRESS,TEMP,RELH
      IF (TIME.GT.9990) GOTO 140
C 
C     Assign end-of-drop data to permanent variables 
C     ----------------------------------------------
      IXEOD = IXEOD + 1
      IF (IXEOD.GT.NEODODW)  THEN            !Check if too big    
         WRITE(LUT,'(/,"*** TOO MANY EODs, TIME: ",F6.0,/)') TIME
      ELSE 
         EOD(IXEOD,1) = TIME
         EOD(IXEOD,2) = PRESS
         EOD(IXEOD,3) = TEMP
         EOD(IXEOD,4) = RELH
         ENDIF
      GOTO 120
C 
140   IF (IXEND.EQ.0) THEN
         WRITE(LUT,'(/," *** FATAL INPUT ERROR: NO VALID TIMES",/)') 
         STOP 
         ENDIF
C
      NRECS = IXEND
      RETURN
C 
900   WRITE(LUT,'(/," *** FATAL ERROR OPENING INPUT FILE ",A,/)') 
     *   FILENAME 
      STOP
      END 
C
C
C
C     -----------------------------------------------
      SUBROUTINE GETFNAME(NTYPE,FILENAME)
	
C
C     Forms file name to open
C     -----------------------------------------------
C 
C 
      PARAMETER (MXRC = 9000, NINVAR = 12) 
C
      CHARACTER*4 FTYPE
      CHARACTER*9 DIREC
      character*100  FILENAME
      CHARACTER*4 VERSION
      LOGICAL AUTOTEMP
C 
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
C
C
      IF (NTYPE.EQ.1) THEN
         FTYPE = '.tmp'
         DIREC = './'
         IF (NSITE.EQ.2) DIREC = 'sondedata'
         ENDIF
      IF (NTYPE.EQ.2) THEN
         FTYPE = '.pds'
         DIREC = 'AVAPS_fnl'
         ENDIF
C
      IF (NSNDTYPE.EQ.1)
     *   WRITE(FILENAME,'(A9,"/odw",I5.5,A4)') DIREC,ID,FTYPE
      IF (NSNDTYPE.GT.1)
     *   WRITE(FILENAME,'(A9,"/g",I9.9,A4)') DIREC,ID,FTYPE
      RETURN
      END
C
C
C
C     -----------------------------------------------
      SUBROUTINE OPENSONDE(FILENAME,NTYPE,NRECS)
	
C
C     Retrieves last working file of sonde data
C     -----------------------------------------------
C 
      PARAMETER (MXRC = 9000, NOPT = 16, NINVAR = 12) 
      PARAMETER (NEOD = 100) 
      PARAMETER (NPROC = 50)
C 
      CHARACTER*4 STNS,VER
      character*100  FILENAME
      character*100  COMMENT
      CHARACTER*12 PROCTIME
      CHARACTER*4 VERSION
      LOGICAL NEEDOPT, PDSFILE
      LOGICAL AUTOTEMP
C 
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
      COMMON /OPTIONS/ NEEDOPT(NOPT)
      COMMON /EODDATA/ EOD(NEOD,NINVAR)
      COMMON /PROCESSING/ NPSTP,PROC(NPROC)
C
C
      OPEN(LUFW,FILE=FILENAME,FORM='UNFORMATTED',STATUS='OLD',
     *     ERR=900)
C
C
C     For processed file, read extra header line with version #
C     ---------------------------------------------------------
      NR = MXRC
      PDSFILE = .FALSE.
      IF (NTYPE.EQ.2) THEN
         PDSFILE = .TRUE.
         READ(LUFW,ERR=900) VER
         READ(LUFW,ERR=900) NRECS
         NR = NRECS
         READ(VER,'(F4.2)') RVER
         IF (RVER.GE.1.22) READ(LUFW,ERR=900) COMMENT         
         IF (COMMENT.EQ.' ') COMMENT = 'None'
         PROCTIME = ' '
         IF (RVER.GE.1.29) READ(LUFW,ERR=900) PROCTIME
         ENDIF
C
C
C     Now read rest of file
C     ---------------------
      READ(LUFW,ERR=900) LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *     NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      READ(LUFW,ERR=900) ID,IXEND,
     *     ((SNDDAT(I,J),I=1,NR),J=1,NINVAR),
     *     ((ISNDFLG(I,J),I=1,NR),J=1,4)
      READ(LUFW,ERR=900) IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,
     *     DPL,HTL,WDL,WSL,STNS 
      READ(LUFW,ERR=900) NEEDOPT,EOD,(X(I),I=1,NR),
     *     (Y(I),I=1,NR),(Z(I),I=1,NR)     
C
C
C     Read processing record, if available.
C     ------------------------------------
      IF (RVER.GE.1.29 .OR. NTYPE.EQ.1) THEN
         READ(LUFW,ERR=900) NPSTP,(PROC(I),I=1,NPROC)
         CALL DISPPROC(PDSFILE,RVER,PROCTIME,COMMENT)
         ENDIF
C
C
      CLOSE(LUFW)
      RETURN
C
C     Error on file
C     -------------
900   WRITE(LUT,'(/," *** ERROR READING FILE ",A,/)') FILENAME 
      RETURN
      END
C
C
C
C     -----------------------------------------------------
      SUBROUTINE DISPPROC(PDSFILE,RVER,PROCTIME,COMMENT)
	
C
C     Displays processing record
C     -----------------------------------------------------
C
      PARAMETER (MXRC = 9000, NINVAR = 12, NPROC = 50)
C
      CHARACTER*12 PROCTIME
      character*100  COMMENT
      CHARACTER*3 HYDANC
      CHARACTER*1 ESTP, DYNT
      LOGICAL PDSFILE
C
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /PROCESSING/ NPSTP, PROC(NPROC)
C
C
      HYDANC = 'MSG'
      IF (PROC(7).EQ.1.) HYDANC = 'FLT'
      IF (PROC(7).EQ.2.) HYDANC = 'SFC'
      ESTP = 'N'
      IF (PROC(8).EQ.1.) ESTP = 'Y'
      DYNT = 'N'
      IF (PROC(9).EQ.1.) DYNT = 'Y'
C
      WRITE(LUT,'(/,1X,75("="))')
      WRITE(LUT,'(" Processing record for sonde ",I9.9,".",/,
     *      " Launch Date: ",3I2.2,"  Time: ",I6.6," UTC.")') 
     *      ID,IYR,IMO,IDY,TIML
      WRITE(LUT,'(1X,75("-"))')
      IF (PDSFILE) WRITE(LUT,'(" Sonde ",I9.9,
     *      " processed by editsonde V",
     *      F4.2," on ",a12,".")') ID,RVER,PROCTIME
C
C
      WRITE(LUT,'(" Bias corrections:  PR = ",F5.1," mb  TE = ",
     *            F5.1," C  RH = ",F5.1," %")') 
     *            PROC(3),PROC(4),PROC(5)
      WRITE(LUT,'(" Filters applied (s):   PTH = ",f5.1,"   WIND = ",
     *            F5.1)') PROC(1),PROC(2)
      WRITE(LUT,'(" Splash PR = ",F6.1,
     *            " mb    Hydrostatic anchor = ",a3)')
     *            PROC(6),HYDANC
      WRITE(LUT,'(" Estimated PR used (",A1,
     *            ")    Dynamic T correction (",A1")")') ESTP, DYNT
C
      IF (PDSFILE) WRITE(LUT,'(" Comments: ",A60)') COMMENT
      WRITE(LUT,'(1X,75("="))')
      WRITE(LUT,*)
C
      IF (.NOT.PDSFILE) RETURN
C
C
C     Echo to log file
C     ----------------
      WRITE(LUPR,'(1X,75("="))')
      WRITE(LUPR,'(" Processing record for sonde ",I9.9,".",/,
     *      " Launch Date: ",3I2.2,"  Time: ",I6.6," UTC.")') 
     *      ID,IYR,IMO,IDY,TIML
      WRITE(LUPR,'(1X,75("-"))')
      IF (PDSFILE) WRITE(LUPR,'(" Sonde ",I9.9,
     *      " processed by editsonde V",
     *      F4.2," on ",a12,".")') ID,RVER,PROCTIME
C
C
      WRITE(LUPR,'(" Bias corrections:  PR = ",F5.1," mb  TE = ",
     *            F5.1," C  RH = ",F5.1," %")') 
     *            PROC(3),PROC(4),PROC(5)
      WRITE(LUPR,'(" Filters applied (s):   PTH = ",f5.1,"   WIND = ",
     *            F5.1)') PROC(1),PROC(2)
      WRITE(LUPR,'(" Splash PR = ",F6.1,
     *            " mb    Hydrostatic anchor = ",a3)')
     *            PROC(6),HYDANC
      WRITE(LUPR,'(" Estimated PR used (",A1,
     *            ")    Dynamic T correction (",A1")")') ESTP, DYNT
C
      IF (PDSFILE) WRITE(LUPR,'(" Comments: ",A60)') COMMENT
      WRITE(LUPR,'(1X,75("="))')
      WRITE(LUPR,*)
C
      RETURN
      END
C
C
C
C     -----------------------------------------------------
      FUNCTION GSNDFALL(PR,TE,BAD)
C
C     Function returns theoretical fall rate of GPS sonde.
C     If TE is missing, mean West Indies hurricane season
C     TE is used.  Pressure must be valid to compute fall.
C     -----------------------------------------------------
C
      PARAMETER (NSTD = 20)
      DIMENSION PSTD(20), TSTD(20)
C
      DATA PSTD/100.,150.,200.,250.,300.,350.,400.,450.,500.,550.,
     1          600.,650.,700.,750.,800.,850.,900.,950.,1000.,1050./
      DATA TSTD/-73.5,-67.6,-55.2,-43.3,-33.2,-24.8,-17.7,-11.9,
     1          -6.9,-2.5,1.4,5.1,8.6,11.8,14.6,17.3,19.8,23.0,
     2          26.0,26.3/
C
C
C     Set constants
C     -------------
      SMASS = 0.390                 !	MASS OF SONDE (KG)
      CD = .63                      !   DRAG COEFF
      AREA = 0.09                   !   AREA OF PARACHUTE (30 CM x 30 CM)
      GRAV = 9.8                    !   GRAVITY
C
C
      IF (PR.LE.BAD) GOTO 900
C
      IF (TE.GT.BAD) THEN
         T = TE
         ELSE
         CALL POLATE(NSTD,PSTD,TSTD,PR,T,M,BAD)
         IF (T.EQ.BAD) GOTO 900
         ENDIF
C
C     Compute density
C     ---------------
      RHO = (PR * 100.)/(287.*(T+273.16))
C
C     Compute fall speed - theoretical
C     --------------------------------
      GSNDFALL = -((2.0*SMASS*GRAV)/(CD*AREA*RHO))**0.5
C
C     Add empirical fudge factor
C     --------------------------
      GSNDFALL = GSNDFALL*1.15                          
      RETURN
C
C
900   GSNDFALL = BAD
      RETURN
      END
C
C
C
C     -----------------------------------------------
      FUNCTION TEMPTC(Z,BAD)
C
C     Temperature time constant of RS-80 sensor.
C     Designed for tropical atmosphere.
C     Equation provided by Hal Cole.  Output in (s).
C     -----------------------------------------------
C
      A = 2.47
      B = 7.81E-5
      C = 3.70E-9
      D = 6.90E-14
C
      IF (Z.EQ.BAD) THEN
         TEMPTC = BAD
         RETURN
         ENDIF
C
      IF (Z.LT.0) THEN
         TEMPTC = A
         RETURN
         ENDIF
C
      TEMPTC = A + B*Z + C*Z**2.0 + D*Z**3.0
C
      RETURN
      END
C
C
C
C     ------------------------------------------------
      FUNCTION RHTC(Z,BAD)
C
C     RH time constant of RS-80 sensor arm (dominant).
C     Designed for tropical atmosphere.
C     Equation provided by Hal Cole.  Output in (s).
C     ------------------------------------------------
C
      A = 8.48
      B = 2.1E-4
      C = 1.91E-8
      D = 1.62E-13
C
      IF (Z.EQ.BAD) THEN
         RHTC = BAD
         RETURN
         ENDIF
C
      IF (Z.LT.0) THEN
         RHTC = A
         RETURN
         ENDIF
C
      RHTC = A + B*Z + C*Z**2.0 + D*Z**3.0
C
      RETURN
      END
C
C
C
C     ------------------------------------- 
      SUBROUTINE PRINTAC(DISP,PRINT)
C 
C     Prints a/c data in list header format 
C     ------------------------------------- 
C 
      PARAMETER (MXRC = 9000, NINVAR = 12)
      CHARACTER*1 ALAT,ALON,OMST1(8)
      CHARACTER*4 STNS
      CHARACTER*2 OMST2(8),OMOUT(4)
      LOGICAL DISP,PRINT
C 
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
C 
      DATA OMST2/'JA','NW','LI','HA','ND','LA','AR','AU'/
      DATA OMST1/'H','A','B','C','D','E','F','G'/
C 
      ALAT = 'N'
      ALON = 'W'
      DATE = 10000.*IYR+100*IMO+1.*IDY
      IF (RLAT.LT.0.) ALAT = 'S'
      IF (RLON.LT.0.) ALON = 'E'
C
      IF (STNS.EQ.'GPS ') THEN
         OMOUT(1) = 'GP'
         OMOUT(2) = '  '
         OMOUT(3) = '  '
         OMOUT(4) = '  '
         ELSE
         DO 100 L=1,4
            OMOUT(L) = '  '
            DO 110 K = 1,8
               IF (STNS(L:L).EQ.OMST1(K)) OMOUT(L) = OMST2(K)
110            CONTINUE
100      CONTINUE
         ENDIF
C 
      IF (DISP) WRITE(LUT,900)
      IF (PRINT) WRITE(LUPR,901)
900   FORMAT(///,1X,75('='))
901   FORMAT(/,1X,75('='))
      IF (DISP) WRITE(LUT,910) DATE,ABS(RLAT),ALAT,TEL,PRL,WDL
      IF (PRINT) WRITE(LUPR,910) DATE,ABS(RLAT),ALAT,TEL,PRL,WDL
910   FORMAT(' Date: ',F6.0,3X,'Lat: ',F6.2,1X,A1,3X,'TA: ',F6.1, 
     *       1X,'C   PS: ',F6.1,1X,'mb   WD: ',f6.0,' deg') 
      IF (DISP)
     * WRITE(LUT,920) NINT(TIML),ABS(RLON),ALON,DPL,NINT(HTL),WSL
      IF (PRINT) 
     * WRITE(LUPR,920) NINT(TIML),ABS(RLON),ALON,DPL,NINT(HTL),WSL
920   FORMAT(' Time: ',I6.6,3X,'Lon: ',F6.2,1X,A1,3X,'TD: ',F6.1, 
     *       1X,'C   GA: ',I6,1X,'m    WS:',F6.1,' m/s')
      IF (DISP) WRITE(LUT,930) ID,RHL,(OMOUT(L),L=1,4)
      IF (PRINT) WRITE(LUPR,930) ID,RHL,(OMOUT(L),L=1,4)
930   FORMAT(' SID:  ',I9.9,16X,'RH: ',F6.1,' %   Navaid:   ',
     *       4(A2,1X))
      IF (DISP) WRITE(LUT,935)
      IF (PRINT) WRITE(LUPR,935)
935   FORMAT(1X,75('='),/)
C 
      RETURN
      END 
C 
C 
C 
C     ------------------------------------
      SUBROUTINE SNDNG(DISP,PRINT)
C 
C     PROGRAM DISPLAYS DATA VS TIME
C     ------------------------------------
C 
      PARAMETER (MXRC = 9000, NINVAR = 12)
C 
      CHARACTER*1 CFLG
      CHARACTER*1 ANS 
      CHARACTER*4 STNS
      LOGICAL PAUSE, PRINT, DISP
C 
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
C
C
C 
      PAUSE = .TRUE.
      CALL PRINTAC(DISP,PRINT)
C 
      IF (DISP) WRITE(LUT,320)
      IF (PRINT) WRITE(LUPR,320)
320   FORMAT('   IX   TIME    PR        TE        RH     GA',
     *'     WD      WS      WQ     VV',/, 
     *       '         (s)   (mb)       (C)       (%)    (m)',
     *'  (deg)   (m/s)          (m/s)',/,
     *       ' ----------------------------------------------', 
     *'-----------------------------')  
C 
      KK = 0
      K = 0
325   K = K+1
C
      IF (DISP)
     * WRITE(LUT,330) K,SNDDAT(K,1),SNDDAT(K,2),CFLG(ISNDFLG(K,1)), 
     *               SNDDAT(K,3),CFLG(ISNDFLG(K,2)),SNDDAT(K,4), 
     *               CFLG(ISNDFLG(K,3)),NINT(SNDDAT(K,9)),
     *               NINT(SNDDAT(K,6)),SNDDAT(K,7),
     *               CFLG(ISNDFLG(K,4)),
     *               NINT(SNDDAT(K,8)),SNDDAT(K,10)
      IF (PRINT)
     * WRITE(LUPR,330) K,SNDDAT(K,1),SNDDAT(K,2),CFLG(ISNDFLG(K,1)), 
     *               SNDDAT(K,3),CFLG(ISNDFLG(K,2)),SNDDAT(K,4), 
     *               CFLG(ISNDFLG(K,3)),NINT(SNDDAT(K,9)),
     *               NINT(SNDDAT(K,6)),SNDDAT(K,7),
     *               CFLG(ISNDFLG(K,4)),
     *               NINT(SNDDAT(K,8)),SNDDAT(K,10)
330      FORMAT(1X,I4,1X,F6.1,2X,F6.1,1X,A1,F8.2,1X,A1,F7.1,1X,A1,
     *          I6,2X,I4,2X,F6.1,1X,A1,2X,I4,1X,F6.1)
         KK = KK+1
         IF (MOD(KK-24,30).EQ.0 .AND. K.LT.IXEND .AND. PAUSE
     *        .AND. DISP) THEN
              WRITE(LUT,
     *        '(/," Enter (CR) for next page, (W) for whole list, ",
     *            "(Q) to quit: ",$)') 
              READ(*,'(A)') ANS 
              IF (ANS.EQ.'Q' .OR. ANS.EQ.'q') RETURN
              IF (ANS.EQ.'W' .OR. ANS.EQ.'w') PAUSE = .FALSE.
              WRITE(LUT,320)
              K = K-1 
              ENDIF 
      IF (K.LT.IXEND) GOTO 325
C 
      IF (DISP) WRITE(LUT,'(1x,75("-"),/)')
      IF (PRINT) WRITE(LUPR,'(1x,75("-"),/)')
C
      RETURN
      END 
C 
C 
C 
C     ------------------------------------------
      SUBROUTINE LVLINT(PINT,DISP,PRINT) 
C 
C     Lists data vs pressure, at PINT intervals.
C     If PINT = 0, then only the mandatory
C     levels are displayed. 
C     ------------------------------------------
C 
      PARAMETER (MXRC = 9000, NINVAR = 12)
      PARAMETER (NMANL = 11)
      LOGICAL MANDATORY, DISP, PRINT, AUTOTEMP
      CHARACTER*1 CFLG
      CHARACTER*3 FLAG
      CHARACTER*4 STNS
      CHARACTER*4 VERSION
      DIMENSION PRMAN(NMANL)
C 
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
C 
      DATA PRMAN/100.,150.,200.,250.,300.,400.,500.,
     *           700.,850.,925.,1000./ 
C
C
C 
      MANDATORY = .FALSE.
      IF (PINT.EQ.0) THEN 
         PINT = 25. 
         MANDATORY = .TRUE. 
         ENDIF
C 
      CALL PRINTAC(DISP,PRINT)
      IF (DISP) WRITE(LUT,140)
      IF (PRINT) WRITE(LUPR,140)
140   FORMAT(
     *   '   PR        TE        TD      RH        GA    ',
     *   '  WD       WS       WQ' 
     *,/,'  (mb)       (C)       (C)     (%)       (m)    ',
     *   '(deg)   (m/s)         '
     *,/,' ',70('-'))
C 
C 
C     Convert WD,WS into U,V for interpolation
C     ----------------------------------------
C 
200   DO 210 L = 1,IXEND
         X(L) = UCMP(SNDDAT(L,6),SNDDAT(L,7))
         Y(L) = VCMP(SNDDAT(L,6),SNDDAT(L,7))
210      CONTINUE 
C 
C 
C     Find initial pressure level...begin looping down
C     ------------------------------------------------
C 
      PNOW = 50. 
      PRBEGIN = 1100.
      PREND = 0.
      DO 250 L = 1,IXEND
         PR = SNDDAT(L,2)
         IF (PR.GT.PREND) PREND = PR
         IF (PR.LT.PRBEGIN .AND. PR.NE.-999.) PRBEGIN = PR
250      CONTINUE
      IF (PRL.LT.PRBEGIN .AND. PRL.GT.0.) PRBEGIN = PRL 
      IF (PRBEGIN.EQ.1100. .OR. PREND.EQ.0.) GOTO 900
C 
220   PNOW = PNOW + PINT
      IF (PNOW .LT. PRBEGIN) GOTO 220 
      IF (PNOW .GE. PREND) GOTO 400 
      IF (MANDATORY) THEN 
         DO 225 I = 1,NMANL
              IF (PNOW.EQ.PRMAN(I)) GOTO 300
225           CONTINUE
         GOTO 220 
         ENDIF
C 
C 
C     Interpolate from data arrays for pressure = PNOW
C     ------------------------------------------------
C 
300   CALL POLATE(IXEND,SNDDAT(1,2),SNDDAT(1,3),PNOW,TNOW,M,-999.)
      FLAG(1:1) = CFLG(ISNDFLG(M,2))
      CALL POLATE(IXEND,SNDDAT(1,2),SNDDAT(1,4),PNOW,HNOW,M,-999.)
      FLAG(2:2) = CFLG(ISNDFLG(M,3))
      CALL POLATE(IXEND,SNDDAT(1,2),SNDDAT(1,9),PNOW,ZNOW,M,-999.)
      CALL POLATE(IXEND,SNDDAT(1,2),X,PNOW,UNOW,M,-999.)
      CALL POLATE(IXEND,SNDDAT(1,2),Y,PNOW,VNOW,M,-999.)
      FLAG(3:3) = CFLG(ISNDFLG(M,4))
      CALL POLATE(IXEND,SNDDAT(1,2),SNDDAT(1,8),PNOW,WQNOW,M,-999.) 
      WSNOW = WSCOMP(UNOW,VNOW) 
      WDNOW = WDCOMP(UNOW,VNOW) 
C 
C     Check if level falls between launch and first ODW pr
C     ----------------------------------------------------
      IF (PNOW.LT.SNDDAT(1,2) .AND. PNOW.GE.PRL .AND. 
     *        ABS(SNDDAT(1,2)-PRL).LE.50.) THEN 
              PRS = SNDDAT(1,2) 
              TNOW = POLATE2(PRL,PRS,TEL,SNDDAT(1,3),PNOW,-999.)
              HNOW = POLATE2(PRL,PRS,RHL,SNDDAT(1,4),PNOW,-999.)
              UNOW = POLATE2(PRL,PRS,UCMP(WDL,WSL),X(1),PNOW,-999.)
              VNOW = POLATE2(PRL,PRS,VCMP(WDL,WSL),Y(1),PNOW,-999.)
              WSNOW = WSCOMP(UNOW,VNOW) 
              WDNOW = WDCOMP(UNOW,VNOW) 
              ZNOW = HYDROZ(PRS,SNDDAT(1,3),SNDDAT(1,4),
     *               SNDDAT(1,9),PNOW,TNOW,HNOW,2,-999.)
              FLAG(1:1) = CFLG(ISNDFLG(1,2))
              FLAG(2:2) = CFLG(ISNDFLG(1,3))
              FLAG(3:3) = CFLG(ISNDFLG(1,4))
              ENDIF 
C 
C 
      DPNOW = DEWPT(TNOW,HNOW)
      IF (DISP) WRITE(LUT,330) PNOW,TNOW,FLAG(1:1),DPNOW,HNOW,FLAG(2:2),
     *               NINT(ZNOW),NINT(WDNOW),WSNOW,FLAG(3:3),NINT(WQNOW)
      IF (PRINT) WRITE(LUPR,330) PNOW,TNOW,FLAG(1:1),DPNOW,
     *           HNOW,FLAG(2:2),NINT(ZNOW),NINT(WDNOW),WSNOW,FLAG(3:3),
     *           NINT(WQNOW)
330   FORMAT(1X,F6.1,2X,F7.1,1X,A1,2(1X,F7.1),1X,A1,3X,I5,3X,I5,3X,F6.1,
     *       1X,A1,3X,I4) 
C
C
      GOTO 220
C 
C 
C     Write out last record 
C     --------------------- 
C 
400   DPNOW = DEWPT(SNDDAT(IXEND,3),SNDDAT(IXEND,4))
      FLAG(1:1) = CFLG(ISNDFLG(IXEND,2))
      FLAG(2:2) = CFLG(ISNDFLG(IXEND,3))
      FLAG(3:3) = CFLG(ISNDFLG(IXEND,4))
      IF (DISP) WRITE(LUT,330) 
     *               SNDDAT(IXEND,2),SNDDAT(IXEND,3),FLAG(1:1),DPNOW,
     *               SNDDAT(IXEND,4),FLAG(2:2),
     *               NINT(SNDDAT(IXEND,9)),NINT(SNDDAT(IXEND,6)),
     *               SNDDAT(IXEND,7),FLAG(3:3),NINT(SNDDAT(IXEND,8))
      IF (DISP) WRITE(LUT,'(1X,70("-"))')
      IF (PRINT) WRITE(LUPR,330) 
     *               SNDDAT(IXEND,2),SNDDAT(IXEND,3),FLAG(1:1),DPNOW,
     *               SNDDAT(IXEND,4),FLAG(2:2),
     *               NINT(SNDDAT(IXEND,9)),NINT(SNDDAT(IXEND,6)),
     *               SNDDAT(IXEND,7),FLAG(3:3),NINT(SNDDAT(IXEND,8))
      IF (PRINT) WRITE(LUPR,'(1X,70("-"))')
C 
C
      RETURN
C
C
900   WRITE(LUT,'(/," *** NO VALID PRESSURES IN SOUNDING ***",/)')
      RETURN
      END 
C
C
C
C     ----------------
      FUNCTION CFLG(I)
C     ----------------
C 
      CHARACTER*1 CFLG
      IF (I.EQ.0) CFLG = ' ' 
      IF (I.EQ.1) CFLG = 'R' 
      IF (I.EQ.2) CFLG = 'M' 
      IF (I.EQ.3) CFLG = 'I' 
      IF (I.EQ.4) CFLG = 'D' 
      IF (I.EQ.5) CFLG = 'S' 
      RETURN
      END 
C 
C 
C 
C     --------------------------------------------------------
      SUBROUTINE CXFLG(ANS,IVFLG) 
C 
C     Returns numerical flag value from character input 
C     --------------------------------------------------------
C 
      CHARACTER*1 ANS 
C 
      IVFLG = -1
C 
      IF (ANS.EQ.'K' .OR. ANS.EQ.'k') IVFLG = 0 
      IF (ANS.EQ.'R' .OR. ANS.EQ.'r') IVFLG = 1 
      IF (ANS.EQ.'M' .OR. ANS.EQ.'m') IVFLG = 2 
      IF (ANS.EQ.'D' .OR. ANS.EQ.'d') IVFLG = 4 
      IF (ANS.EQ.'S' .OR. ANS.EQ.'s') IVFLG = 5 
      IF (ANS.EQ.'Q' .OR. ANS.EQ.'q') IVFLG = -2 
      RETURN
      END 
C
C
C
C     --------------------------------------------------------
      SUBROUTINE INDEX(ANS,IXDAT,IXFLG) 
C 
C     Returns array index positions for SNDDAT and ISNDFLG
C     arrays. 
C     --------------------------------------------------------
C 
      CHARACTER*1 ANS 
C 
      IXDAT = -1
      IXFLG = -1
C 
      IF (ANS.EQ.'Q' .OR. ANS.EQ.'q') THEN
         IXDAT = 0
         IXFLG = 0
         ENDIF
      IF (ANS.EQ.'P' .OR. ANS.EQ.'p') THEN
         IXDAT = 2
         IXFLG = 1
         ENDIF
      IF (ANS.EQ.'T' .OR. ANS.EQ.'t') THEN
         IXDAT = 3
         IXFLG = 2
         ENDIF
      IF (ANS.EQ.'H' .OR. ANS.EQ.'h') THEN
         IXDAT = 4
         IXFLG = 3
         ENDIF
      IF (ANS.EQ.'W' .OR. ANS.EQ.'w') THEN
         IXDAT = 6
         IXFLG = 4
         ENDIF
      IF (ANS.EQ.'U' .OR. ANS.EQ.'u') THEN
         IXDAT = 90
         IXFLG = 4
         ENDIF
      RETURN
      END 

C     --------------------------------------------------------
C     SONDELIB
C
C     Collection of subroutines relating to soundings.
C     This files has three main sections:  
C
C         1) Routines relating to Ooyama's spline filtering.
C         2) Mathematical or meteorological routines
C         3) Routines for message decoding
C
C     --------------------------------------------------------
C
C
C
C
C     --------------------------------------------------------
C     Section 1: Routines relating to Ooyama spline filtering.
C     --------------------------------------------------------
C
C
C     --------------------------------------------------------
      SUBROUTINE VICSPL(XT,XW,XDAT,NXP,NXPDIM,KDAT,YNB,YNT,NX,
     *                  YDCWL,KYBC1,KYBC2,YBCWL1,YBCWL2,IERR)
C
C     Driver subroutine for Ooyama 1-D spline smoothing.
C     VICSPL does the filtering, output obtained through 
C     subroutine SPOTVAL.
C
C     01 Aug 96:  Mean subtracted first to improve accuracy.
C     --------------------------------------------------------
C
C
C     NDIM   = Maximum allowable number of nodes
C     LDXDIM = Maximum allowable number of data points
C     KDTDIM = Maximum number of input data variables
C     KVDIM  = Maximum number of variables to be filtered
C
C     XT     = Array of independent variable (e.g., time)
C     XW     = Array of relative weights (0<=wt<=1). If weight < 0 or > 1,
C              subroutine sets relative weight to be 1.
C     XDAT   = Array of input data to be filtered.  First array to
C              be filtered is XDAT(1...NXP,1), last is XDAT(1...NXP,KDAT)
C     NXP    = Number of data points in arrays XT, XW, and XDAT.
C     NXPDIM = 1st dimension of XDAT array in calling program.
C     KDAT   = Number of variables in XDAT to be filtered (1 through KDAT).
C
C     YNB    = Lower domain boundary (real units, e.g., seconds)
C     YNT    = Upper domain boundary
C     NX     = Number of nodal INTERVALS spanning (YNB,YNT)
C     YDCWL  = Filter cutoff wavelength (real units)
C     KYBC1  = Boundary condition type at YNB boundary
C     KYBC2  = Boundary condition type at YNT boundary
C     YBCWL1 = Scale parameter for YNB boundary condition (real units)
C     YBCWL2 = Scale parameter for YNT boundary condition (real units)
C     IERR   = 0 if smoothing went ok, 1 if error
C
C
      PARAMETER (NDIM=240)
      PARAMETER (LDXDIM=9000)
      PARAMETER (KDTDIM=2, KVDIM=2*KDTDIM)
C
      CHARACTER*121 IDSPC
      CHARACTER*8 IDSPD1,IDSPD2*72
      DIMENSION XT(NXP),XW(NXP),XDAT(NXPDIM,KDAT),WCELL(NDIM)
      LOGICAL DENCORR
C
      COMMON /BIG/ A(KDTDIM,NDIM), FD(KDTDIM,LDXDIM), FKMN(KVDIM,NDIM),
     *             AMM(NDIM,NDIM), DMM(NDIM,NDIM)
      COMMON /LITTLE/ ZS(NDIM,NDIM), ZD(NDIM,7)
      COMMON /SP2DTPT/ LDX, YD(LDXDIM), WD(LDXDIM)
      COMMON /SP2DCHR/ IDSPD1, IDSPD2
C
      DATA KDC/3/
      DATA DENCORR/.TRUE./
C
C
C     Check and load in input data arrays
C     -----------------------------------
      IERR = 0
      CALL XDATZ(KDAT,NXP,NXPDIM,XT,XW,XDAT,YNB,YNT,FD,IERR)
      IF (IERR.EQ.1) RETURN
C
      YN1 = YNB
      YNXX = YNT
      DYN = (YNXX-YN1)/NX
      RATIO = YDCWL/DYN
C
C
C     Correct weights for local density
C     ---------------------------------
      IF (DENCORR) THEN
         DO 110 NCZ = 1,NX
110        WCELL(NCZ)= 0.
         DO 120 L = 1,LDX
           NC = MIN0(INT((YD(L)-YN1)/DYN)+1,NX)
           WCELL(NC) = WCELL(NC)+WD(L)
120        CONTINUE
         DO 140 L = 1,LDX
           NC = MIN0(INT((YD(L)-YN1)/DYN)+1,NX)
           WFAC = 1./WCELL(NC)
           WD(L) = MIN(WD(L)*WFAC,1.0)
140        CONTINUE
         ENDIF
C
C
C     Check if filter is a reasonable multiple of the nodal interval
C     --------------------------------------------------------------
      IF (RATIO.LT.2.0 .OR. RATIO.GT.15.0) THEN
	 WRITE(1,900) YDCWL,DYN,RATIO
900	 FORMAT(/,' ** IMPROPER RATIO OF FILTER TO NODAL SPACING **',/,
     *            ' FILTER CUTOFF  = ',F10.3,/,
     *            ' NODAL INTERVAL = ',F10.3,/,
     *            ' RATIO = ',F10.3, ' (OK IF BETWEEN 2 AND 15)')
	 STOP
	 ENDIF
C
      CALL SETSPID(IDSPC)
      CALL SETSPND(NX,YN1,YNXX)
      CALL SETSPBC(KYBC1,KYBC2,YBCWL1,YBCWL2)
      CALL SETSPDC(KDC,YDCWL)
      CALL SPLOVER
      CALL SPLARGO
C
      KDX = KDAT
      CALL SETKPX(KDX)
      CALL SPLIMBO(1,KDX)
      CALL PRIDSPK(0)
      CALL PRFKMND(A,FKMN,1,KDX)
C
      RETURN
      END
C
C
C
      SUBROUTINE XDATZ(KDAT,NXP,NXPDIM,XT,XW,XDAT,YNB,YNT,FD,IERR)
C
      PARAMETER (NDIM=240)
      PARAMETER (LDXDIM=9000)
      PARAMETER (KDTDIM=2, KVDIM=2*KDTDIM)
C
      COMMON /SP2PRIC/ NX,KDCTYP,KBCTYP(2),YN1,YNXX,YDCWL,
     *                 YBCWL(2),IDSPC
      COMMON /SP2DTPT/ LDX,YD(LDXDIM),WD(LDXDIM)
      COMMON /SP2DCHR/ IDSPD1,IDSPD2
      COMMON /SP2KOMP/ IDVAR(KVDIM),KPX,IDSPK(KVDIM)
      COMMON /MEANS/ FDMEAN(KDTDIM)
C
      CHARACTER*8 IDSPD1,IDSPD2*72,IDSPC*121,IDSPK*40
C
      DIMENSION FD(KDTDIM,LDXDIM)
      DIMENSION XT(NXP),XW(NXP),XDAT(NXPDIM,KDAT)
C
C
      IF (KDAT.GT.KDTDIM) GOTO 930
      IF (NXP.GT.NXPDIM) GOTO 950
C
      DO 10 L = 1,LDXDIM
	 DO 15 KD = 1,KDAT
15	    FD(KD,L) = -999.
	 WD(L) = -999.
	 YD(L) = -999.
10	 CONTINUE
C
      LDX = 0
      DO 100 L = 1, NXP
	 IF (XT(L).EQ.-999.) GOTO 100
         DO 110 KD = 1, KDAT
	    IF (XDAT(L,KD).EQ.-999.) GOTO 100
110	    CONTINUE
	 LDX = LDX + 1
	 IF (LDX.GT.LDXDIM) GOTO 900
	 YD(LDX) = XT(L)
	 WD(LDX) = XW(L)
	 IF (WD(LDX).GT.1. .OR. WD(LDX).LT.0.) WD(LDX) = 1.0
	 DO 120 KD = 1, KDAT
	    FD(KD,LDX) = XDAT(L,KD)
120	    CONTINUE
100      CONTINUE
C
      IF (LDX.LT.2) THEN
         IERR = 1
         RETURN
         ENDIF
C
C
C     Remove mean of data to increase accuracy.
C     Mean replaced in SPOTVAL
C     -----------------------------------------
      DO 220 KD = 1, KDAT
         TOT = 0.
         DO 230 L = 1, LDX
            TOT = TOT + FD(KD,L)
230         CONTINUE
         FDMEAN(KD) = TOT/FLOAT(LDX)
         DO 240 L = 1, LDX
            FD(KD,L) = FD(KD,L) - FDMEAN(KD)
240         CONTINUE
220      CONTINUE
C
      RETURN
C
C
900   WRITE(0,910) LDX,LDXDIM
910   FORMAT(/,' ** LDX EXCEEDS LDXDIM IN XDATZ **',2I5)
      IERR = 1
      RETURN
C
930   WRITE(0,940) KDAT,KDTDIM
940   FORMAT(/,' ** KDAT EXCEEDS KDTDIM IN XDATZ **',2I5)
      IERR = 1
      RETURN
C
950   WRITE(0,960) NXP,NXPDIM
960   FORMAT(/,' ** NXP EXCEEDS NXPDIM IN XDATZ **',2I5)
      IERR = 1
      RETURN
      END
C
C
C
C     -----------------------------------------------------------------
      SUBROUTINE VICSETUP(XT,NXP,YDCWL,NX,YNB,YNT,FMIN,IERR,ECHO)
C
C     Returns "appropriate" values of NX,YNB,YNT for user to
C     sent to VICSPL.  YNB and YNT are the smallest and 
C     largest "times" in the XDATA array, +ADDON%.  NX is calculated
C     so that there are at least FMIN nodes for the desired filter.
C     -----------------------------------------------------------------
C
C     XT     = Array of the independent variable, (e.g., time). 
C     NXP    = Number of data points in XT.
C     YDCWL  = Filter cutoff (real units). 
C     FMIN   = Minimum acceptable YDCWL (delta Y units; usually 2.0)
C     ECHO   = .F. to suppress routine output
C
C     Returned values:
C
C     YNB    = Lower domain boundary - addon% (real units, e.g., seconds)
C     YNT    = Upper domain boundary + addon%
C     NX     = Number of nodal intervals spanning (YNB,YNT)
C     IERR   = Error flag: 0 if OK, 1 if error detected   
C
C
C      PARAMETER (NDIM=240)
      PARAMETER (NDIM=240)
      PARAMETER (LDXDIM=9000)
      PARAMETER (ADDON=0.0)
C
      DIMENSION XT(NXP)
      LOGICAL ECHO
C
C
      IERR = 0
      YNB = 99999999.
      YNT = -99999999.
      DO 10 L = 1,NXP
	IF (XT(L).EQ.-999.) GOTO 10
	IF (XT(L).LT.YNB) YNB = XT(L)
	IF (XT(L).GT.YNT) YNT = XT(L)
10	CONTINUE
C
      ADDONP = ADDON*(YNT-YNB)
      YNB = YNB - ADDONP
      YNT = YNT + ADDONP

      IF (YDCWL.GT.(YNT-YNB)) GOTO 200
C
C
      NX = 10
20    DELTAY = (YNT-YNB)/FLOAT(NX)
      RATIOF = YDCWL/DELTAY
      RATIOD = FLOAT(NXP)/FLOAT(NX+1)
      IF (RATIOD.LT.1.) GOTO 950
      IF (RATIOF.LT.FMIN) THEN
	NX = NX+1
	GOTO 20
	ENDIF
      IF (NX.GE.NDIM) GOTO 950
C
C
C     Found marginally acceptable setup, now optimize
C     -----------------------------------------------
30    NX = NX+1
      DELTAY = (YNT-YNB)/FLOAT(NX)
      RATIOF = YDCWL/DELTAY
      RATIOD = FLOAT(NXP)/FLOAT(NX+1)
      IF (RATIOD.LT.1. .OR. NX.EQ.NDIM .OR. RATIOF.GT.15.) THEN
         NX = NX-1
         GOTO 50
         ENDIF
      IF (RATIOF.LT.4. .OR. RATIOD.GT.2.0) GOTO 30
C      IF (MOD(NXP,NX+1).NE.0) GOTO 30
C
C
C     Found an acceptable NX...calculate final parameters
C     ---------------------------------------------------
50    DELTAY = (YNT-YNB)/FLOAT(NX)
      RATIOF = YDCWL/DELTAY
      RATIOD = FLOAT(NXP)/FLOAT(NX+1)
C
      IF (ECHO) WRITE(0,900) YDCWL,NX+1,DELTAY,YNB,YNT,RATIOF,RATIOD
900   FORMAT(/," For requested filter wavelength = ",f10.3,/,
     *         " VICSETUP found ",i3," nodes.   DELTAY = ",F10.2,/,
     *         " YNB, YNT = ",2F10.2,/,
     *         " Filter/nodal interval ratio = ",f6.2,/,
     *         " # data points / # nodes     = ",f6.2)
      RETURN
C
C
C     TOO MANY NODES
C     --------------------------------------
950   WRITE(0,955) YDCWL,NX+1,YNB,YNT,NDIM,NXP
955   FORMAT(/," *** ERROR IN SUBROUTINE VICSETUP ***",//,
     *         " For filter wavelength = ",f10.3,/,
     *         " VICSETUP has found ",I4," nodes.   YNB, YNT = ",
     *2f10.2,/," # of nodes exceeds maximum allowable value of ",I4,/,
     *         " or exceeds number of input data points NXP = ",i4,//,
     *         " Decrease domain, or increase filter or data density.")
      IERR = 1
      RETURN
C
C
C     Filter is too big or error in domain
C     ------------------------------------
200   WRITE(0,930) YDCWL,YNT-YNB
930   FORMAT(/," *** WARNING IN SUBROUTINE VICSETUP ***",/,
     *         " *** FILTER WAVELENGTH IS TOO LARGE ***",//,
     *         " For filter wavelength YDCWL  = ",f10.2,/,
     *         " VICSETUP has found (YNT-YNB) = ",f10.2,/,
     *         " YDCWL should not exceed domain size.",/)
      IERR = 1
      RETURN
C
      END
C
C
C
      SUBROUTINE SPOTVAL(X,KDAT,FOUT,FOUTD)
C
      PARAMETER (NDIM=240)
      PARAMETER (LDXDIM=9000)
      PARAMETER (KDTDIM=2, KVDIM=2*KDTDIM)
C
      COMMON /BIG/ A(KDTDIM,NDIM), FD(KDTDIM,LDXDIM), FKMN(KVDIM,NDIM),
     *             AMM(NDIM,NDIM), DMM(NDIM,NDIM)
      COMMON /MEANS/ FDMEAN(KDTDIM)
C
      DIMENSION FOUT(KDAT),FOUTD(KDAT)
C
      CALL SPOTINO(X,A,1,KDAT,FOUT)
      CALL SPOTDER(X,A,1,KDAT,FOUTD)
C
C     Add back in means
C     -----------------
      DO 100 KD = 1,KDAT
         FOUT(KD) = FOUT(KD)+FDMEAN(KD)
100      CONTINUE
C
      RETURN
      END
C
C
C
C     ------------------------------------
C     Now we begin with SJL's subroutines.
C     Modified from SPLINE1D.FTN on HP-3.
C     ------------------------------------
C
C
      SUBROUTINE XPRAKMN(A,KP1,KP2) 
C 
      PARAMETER (NDIM=240)
      PARAMETER (KDTDIM=2)
      PARAMETER (KVDIM=2*KDTDIM)
C 
      COMMON /SP2NODE/ NXX,DYN,YN(NDIM) 
C 
      DIMENSION A(KDTDIM,NDIM)
C 
      DO 20 KP=KP1,KP2
D     WRITE(6,21) KP,(YN(N),A(KP,N),N=1,NDIM) 
20    CONTINUE
   21 FORMAT(' .......NODAL AMPLITUDES   KP='I3,/,(10X,F6.1,F13.8,/)) 
      RETURN
      END 
C 
C 
C                 SPLOVER **  11/20/84
C                             SPLOCCO ADDED FOR CORNER CONDITIONS 
C 
C     ++ ASSUMED DIM PARAMETERS ...NDIM=240, LDXDIM=9000, 
C                               ...KDTDIM=2 
C 
C 
C 
      SUBROUTINE SPLOVER
C 
C     ** COMPILATION PARAMETERS **
C 
      PARAMETER (NDIM=240)
      PARAMETER (LDXDIM=9000) 
      PARAMETER (KDTDIM=2)
      PARAMETER (KVDIM=2*KDTDIM)
C 
C 
C     NDIM     ARRAY DIMENSION FOR NODAL INDEX N IN THE Y DIRECTION 
C     LDXDIM   ARRAY DIMENSION FOR DATA POINT INDEX LD
C     KVDIM    NUMBER OF VARIABLES TO BE PARALLEL PROCESSED 
C     KDTDIM   NUMBER OF RAW DATA VARIABLES (T,RH,Z,U,V)
C 
C     ** PRIMARY CONSTANTS FOR SPLINE ANALYSIS PACKAGE
C 
      COMMON /SP2PRIC/ NX,KDCTYP,KBCTYP(2),YN1,YNXX,YDCWL,
     1        YBCWL(2),IDSPC
C 
C     ** SECONDARY CONSTANTS FOR CONVENIENCE, DEFINED BY SPLOVER
C 
      COMMON /SP2NODE/ NXX,DYN,YN(NDIM) 
C 
      COMMON /SP2BCDC/ DY0,DY1,DY2,CBCBF(4),YSDCPI
C 
C       ++ DEFINED BY SPLOCCO,CALLED WITHIN SPLOVER 
      COMMON /SP2CORN/ CORN(4,4)
C 
C     ** DUMMY ARGUMENT 
      CHARACTER*121 IDSPCZ,IDSPC,KENTRY*7 
C 
C     ** LOCAL WORK SPACE 
      DIMENSION YSBCPI(2) 
C 
C     ** NUMERICAL CONSTANTS
      SAVE YSBCPI 
      DATA PI2/6.283185307180/
      DATA KENTRY/'SPLOVER'/
C 
C 
C     ** OVERHEAD CHORES FOR TWO-D SPLINE ANALYSIS AND INVERSE TRANSFORM. 
C 
C     INPUT    ALL THE PRIMARY CONSTANTS IN /SP2PRIC/ 
C              IDSPC .. CODE AND EXPLANATION (80 CHARACTERS)
C 
C     OUTPUT   ALL THE SECONDARY CONSTANTS IN /SP2NODE/ AND /SP2BCDC/ 
C 
C     ENTRIES FOR DEFINING THE PRIMARY CONSTANTS
C              SETSPID(IDSPCZ)
C              SETSPND(NXZ,YN1Z,YNXXZ)
C              SETSPDC(KDCTYPZ,YDCWLZ)
C              SETSPBC(KBCTYP1,KBCTYP2,YBCWL1,YBCWL2) 
C 
C     ENTRIES FOR PRINT OUTPUT
C              PROVERC ... TO PRINT /SP2PRIC/ 
C              PRSBCDC ... TO PRINT /SP2BCDC/ AND XYSBCPI 
C 
C 
C 
C 
C     SECONDARY CONSTANTS 
C 
      NXX=NX+1
      DYN=(YNXX-YN1)/FLOAT(NX)
      YN(1)=YN1 
C 
C     ** CHECK NDIM **
C 
      IF(NXX.GT.NDIM) GO TO 90
C 
      DO 20 N=2,NXX 
   20 YN(N)=YN1+(N-1)*DYN 
C 
C     CONSTANTS FOR BASIS FUNCTION CALCULATION
C 
      DY0=1./DYN
      DY1=3.*DY0
      DY2=6.*DY0**2 
C 
C 
C     ** DC AND BC CONSTANTS ARE NONDIMENSIONALIZED 
C 
C     YDCWL AND YBCWL MUST BE GIVEN IN INTERNAL LENGTH UNITS ... METERS ... 
C 
      Z=1./PI2
      ZC=Z/DYN
      YSDCPI   =ZC*YDCWL
      YSBCPI(1)=ZC*YBCWL(1) 
      YSBCPI(2)=ZC*YBCWL(2) 
C 
C     ** BOUNDARY LINES ... K=1 AT YN1, =2 AT YNXX
C 
C     **BOUNDARY CONDITION TYPES (PERIODIC CONDITION IS NOT INCLUDED, 03/26/79) 
C 
C     KBCTYP= 0 ... YSBF=0
C             1 ... DYSBF=0 
C             2 ... DDYSBF=0
C            10 ... (-/+)YSBCPI*DYSBF+YSBF=0     (- AT YN1, + AT YNXX)
C            20 ... (YSBCPI**2)*DDYSBF+YSBF=0 
C            21 ... (-/+)YSBCPI*DDYSBF+DYSBF=0
C 
C     ** COEFFICIENTS DEFINING YBF AT THE BOUNDARY
C 
      DO 38 K=1,2 
      KK=K+K-1
      IF(KBCTYP(K).EQ. 0) GO TO 31
      IF(KBCTYP(K).EQ. 1) GO TO 32
      IF(KBCTYP(K).EQ. 2) GO TO 33
      IF(KBCTYP(K).EQ.10) GO TO 34
      IF(KBCTYP(K).EQ.20) GO TO 35
      IF(KBCTYP(K).EQ.21) GO TO 36
      GO TO 91
C 
C     TYPE 0
   31 CBCBF(KK)  =-1. 
      CBCBF(KK+1)=-0.25 
        GO TO 38
C     TYPE 1
   32 CBCBF(KK)  = 0. 
      CBCBF(KK+1)= 0.25 
        GO TO 38
C     TYPE 2
   33 CBCBF(KK)  = 0.5
      CBCBF(KK+1)=-0.25 
        GO TO 38
C     TYPE 10 
   34 CBCBF(KK)  =-1./(3.*YSBCPI(K)+1.) 
      CBCBF(KK+1)= 0.25+0.5*CBCBF(KK) 
      GO TO 38
C     TYPE 20 
   35 Z=3.*YSBCPI(K)**2 
      CBCBF(KK)  = 0.5*(Z-1.)/(Z+0.5) 
      CBCBF(KK+1)=-0.25 
        GO TO 38
C     TYPE 21 
   36 CBCBF(KK)  = 0.5*YSBCPI(K)/(0.5+YSBCPI(K))
      CBCBF(KK+1)= 0.25-CBCBF(KK) 
   38 CONTINUE
C 
C 
C     ** SPLOVER DONE **
C 
C 
      ASSIGN 82 TO IDPRINT
      GO TO 80
C 
C 
      ENTRY SETSPID(IDSPCZ) 
      IDSPC=IDSPCZ
      RETURN
C 
      ENTRY SETSPND(NXZ,YN1Z,YNXXZ) 
      NX=NXZ
      YN1 =YN1Z 
      YNXX=YNXXZ
      RETURN
C 
      ENTRY SETSPDC(KDCTYPZ,YDCWLZ) 
      KDCTYP=KDCTYPZ
      YDCWL=YDCWLZ
      RETURN
C 
      ENTRY SETSPBC(KBCTYP1,KBCTYP2,YBCWL1,YBCWL2)
      KBCTYP(1)=KBCTYP1 
      KBCTYP(2)=KBCTYP2 
      YBCWL(1)=YBCWL1 
      YBCWL(2)=YBCWL2 
      RETURN
C 
C 
      ENTRY PROVERC 
C 
      ASSIGN 81 TO IDPRINT
C 
   80 CONTINUE
D     WRITE(6,880) IDSPC,NX,KDCTYP,(KBCTYP(I),I=1,2), 
D    1    YN1,YNXX,DYN,YDCWL/DYN, 
D    2    (YBCWL(I)/DYN,I=1,2)
C 
  880 FORMAT('0.../SP2PRIC/...',/,9X,A121,/,'  NX=',
     1       I3,', (DC)',I2,', (BC)',2I3,',   (Y LIMITS)',2F10.3,/, 
     2       9X,'DELTAY=',F7.3,',  (DCWL [DELTAY UNITS])',F8.1, 
     3       ',  (BCWL [DELTAY UNITS])',2F8.1)
      GO TO IDPRINT,(82,81) 
   81 RETURN
C 
C 
      ENTRY PRSBCDC 
C 
   82 CONTINUE
D     WRITE(6,882) (CBCBF(I),I=1,4),YSDCPI,(YSBCPI(I),I=1,2)
C 
  882 FORMAT('0.../SP2BCDC/...(CBCBF)',4F8.3,',  YSDCPI=',F8.3,/, 
     1                 14X,'(YSBCPI)',2(F8.3,8X)) 
      RETURN
C 
C 
C     ** ERROR EXITS ** 
   90 WRITE(1,890) NXX,NDIM 
  890 FORMAT('0 === SPLOVER FOUND NXX.GT. NDIM) ===',3X,2I6)
      GO TO 99
C 
   91 WRITE(1,891) KBCTYP,YBCWL 
  891 FORMAT('0 === SPLOVER MET UNDEFINED KBCTYP ===',2I12,/,41X,2F12.3)
   99 CALL SPLABORT(KENTRY)
      END 
C 
C
C
      SUBROUTINE SETKPX(KPXZ) 
C 
C          ENTRY PRIDSPC
C          ENTRY PRIDSPD
C          ENTRY PRIDSPK(KP)
C 
      PARAMETER (NDIM=240)
      PARAMETER (LDXDIM=9000) 
      PARAMETER (KDTDIM=2)
      PARAMETER (KVDIM=2*KDTDIM)
C 
C 
      COMMON /SP2PRIC/ NX,KDCTYP,KBCTYP(2),YN1,YNXX,YDCWL,
     1        YBCWL(2),IDSPC
      COMMON/SP2DTPT/LDX,YD(LDXDIM),WD(LDXDIM)
      COMMON/SP2DCHR/IDSPD1,IDSPD2
      COMMON /SP2KOMP/ IDVAR(KVDIM),KPX,IDSPK(KVDIM)
      CHARACTER*8 IDSPD1,IDSPD2*72,IDSPK*40,IDSPC*121,KENTRY*7
C 
C     ** DEFAULT VALUE
      DATA KENTRY/' SETKPX'/
C 
C 
C  ** ENTRY SETKPX(KPXZ)
C 
      KPX=KPXZ
      RETURN
C 
C 
      ENTRY PRIDSPC 
C 
D     WRITE(6,810) IDSPC
  810 FORMAT(9X,'...',A121) 
      RETURN
C 
C 
      ENTRY PRIDSPD 
C 
D     WRITE(6,810) IDSPD
      RETURN
C 
C 
      ENTRY PRIDSPK(KP) 
C 
      KPA=MAX0(KP,1)
      KPB=MIN0(KP,KPX)
C 
      IF(KPA.NE.KPB) GO TO 18 
C 
D     WRITE(6,815) IDSPK(KP)
  815 FORMAT('H',/,30X,A40) 
      RETURN
C 
   18 CONTINUE
      IF(KPX.LT.1.OR.KPX.GT.KVDIM) GO TO 90 
C 
      KPA=1 
      KPB=KPX 
C 
D     WRITE(6,820) KPX
  820 FORMAT('0 === PRIDSPK LISTS ALL ENTRIES IN IDSPK ===   KPX=', I2) 
C 
   22 CONTINUE
D     WRITE(6,822) (K,IDSPK(K),K=KPA,KPB) 
  822 FORMAT(23X,I4,'...',A40)
      RETURN
C 
C 
   90 WRITE(1,890) KPX,KVDIM
  890 FORMAT('0 === PRIDSPK(0) FOUND KPX UNDEFINED OR IN ERROR.',  2I4) 
      CALL SPLABORT(KENTRY)
      END 
C 
C 
C 
      FUNCTION YBF(N,Y) 
C        ENTRY YSBF, DYBF, DYSBF, DDYBF, DDYSBF 
C 
      PARAMETER (NDIM=240)
C 
      COMMON /SP2PRIC/ NX,KDCTYP,KBCTYP(2),YN1,YNXX,YDCWL,
     1        YBCWL(2),IDSPC
      COMMON /SP2NODE/ NXX,DYN,YN(NDIM) 
C 
      COMMON /SP2BCDC/ DY0,DY1,DY2,CBCBF(4),YSDCPI
C 
      CHARACTER*121 IDSPC 
C 
      EQUIVALENCE (NX,LX),(CBCBF(1),CBC11),(CBCBF(2),CBC12),
     1                    (CBCBF(3),CBCX1),(CBCBF(4),CBCX2) 
C 
C     ENTRY YBF 
C     L=L               !WON'T COMPILE WITHOUT THIS LINE--HP(3) VERSION
      ENTRY YSBF(N,Y) 
    5 L=N 
      XS=(Y-YN(L))*DY0
   10 Z=2.-ABS(XS)
      IF(Z)  12, 12, 15
   12 YBF=0.
      RETURN 
   15 F=0.25*Z**3 
      Z=Z-1.
   20 IF(Z)  23, 23, 22 
   22 F=F-Z**3
   23 IF(L-2)  30, 35, 24 
   24 IF(LX-L)  40, 45, 25
   25 YBF=F 
C
C     Need to watch underflow of powers of YBF in calling routines.
C
      RETURN 
C 
C     BOUNDARY ADJUSTMENTS
C 
   30 IF(Z)  25, 25, 32 
   32 F=F+CBC11*(1.-XS)**3
      GO TO 25 
   35 IF(XS)  37, 25, 25
   37 F=F-CBC12*XS**3 
      GO TO  25 
   40 IF(Z)  25, 25, 42 
   42 F=F+CBCX1*(1.+XS)**3
      GO TO 25 
   45 IF(XS)  25, 25, 47
   47 F=F+CBCX2*XS**3 
      GO TO  25 
C 
      ENTRY DYBF(N,Y) 
      D1=DY1
      GO TO 105
      ENTRY DYSBF(N,Y)
      D1=3. 
  105 L=N 
      XS=(Y-YN(L))*DY0
      D1=-SIGN(D1,XS)
  110 Z=2.-ABS(XS)
      IF(Z) 112,112,115
  112 YBF=0.
      RETURN 
  115 F=0.25*Z**2 
      Z=Z-1.
  120 IF(Z) 123,123,122 
  122 F=F-Z**2
  123 IF(L-2) 130,135,124 
  124 IF(LX-L) 140,145,125
  125 YBF=F*D1
      RETURN 
  130 IF(Z) 125,125,132 
  132 F=F+CBC11*SIGN((1.-XS)**2,XS) 
      GO TO 125 
  135 IF(XS) 137,125,125
  137 F=F-CBC12*XS**2 
      GO TO 125 
  140 IF(Z) 125,125,142 
  142 F=F-CBCX1*SIGN((1.+XS)**2,XS) 
      GO TO 125 
  145 IF(XS) 125,125,147
  147 F=F-CBCX2*XS**2 
      GO TO 125 
C 
      ENTRY DDYBF(N,Y)
      D2=DY2
      GO TO 205
      ENTRY DDYSBF(N,Y) 
      D2=6. 
  205 L=N 
      XS=(Y-YN(L))*DY0
  210 Z=2.-ABS(XS)
      IF(Z) 212,212,215
  212 YBF=0.
      RETURN 
  215 F=0.25*Z
      Z=Z-1.
  220 IF(Z) 223,223,222 
  222 F=F-Z 
  223 IF(L-2) 230,235,224 
  224 IF(LX-L) 240,245,225
  225 YBF=F*D2
      RETURN 
  230 IF(Z) 225,225,232 
  232 F=F+CBC11*(1.-XS) 
      GO TO 225 
  235 IF(XS) 237,225,225
  237 F=F-CBC12*XS
      GO TO 225 
  240 IF(Z) 225,225,242 
  242 F=F+CBCX1*(1.+XS) 
      GO TO 225 
  245 IF(XS) 225,225,247
  247 F=F+CBCX2*XS
      GO TO 225 
      END 
C 
C 
      SUBROUTINE XTIMESP(DATE,GMT,IYR,IMO,IDA,HR) 
      IYR=INT(DATE/10000.)
      IMO=INT((DATE-FLOAT(IYR*10000))/100.) 
      IDA=INT(DATE-IYR*10000-IMO*100) 
      HR =FLOAT(NINT(GMT/100.-0.3)) 
      HR =HR+(GMT-HR*100.)/60.
      RETURN
      END 
C 
C 
      SUBROUTINE PRBFTAB(NAME,ZBF)
C 
C 
C     ** ** WARNING ** **  EXTERNAL  ** **
C 
C       THE CALLING ROUTINE MUST DECLARE THE ACTUAL  ZBF  AS  EXTERNAL
C 
C 
C     NAME     CHARACTER NAME OF TH BASIS FUNCTION (8H  ZBF   ) 
C              ZBF IS THE Y COORDINATE
C     ZBF      FUNCTION SUBROUTINE NAME 
C 
C         EXAMPLE OF USAGE ... CALL PRBFTAB(6HDDYSBF,2,DDYSBF)
C 
C 
      PARAMETER (NDIM=240)
C 
C     ** LOCAL PARAMETERS FOR INTERVAL SUBDIVISION
      PARAMETER (ISUBD=10)
      PARAMETER (IXX=4*ISUBD+1) 
C 
      EXTERNAL ZBF
C
      DIMENSION TAB(IXX,5)
      CHARACTER*6 NODE(IXX),KOL(5)*5
C 
      COMMON /SP2NODE/ NXX,DYN,YN(NDIM) 
      DATA NODE/IXX*' '/
      DATA KOL/'  N=1','  N=2','  N=N',' N=NX','N=NXX'/ 
C 
C 
      Z1=YN(1)
      Z4=YN(NXX-4)
      DZ=DYN/ISUBD
      L4=NXX-1
      L5=NXX
C 
      DO 8 I=1,IXX,ISUBD
    8 NODE(I)=' *NODE'
C 
      DO 12 I=1,IXX 
      Z=(I-1)*DZ
      TAB(I,1)=ZBF(1,Z1+Z)
      TAB(I,2)=ZBF(2,Z1+Z)
      TAB(I,3)=ZBF(3,Z1+Z)
      TAB(I,4)=ZBF(L4,Z4+Z) 
   12 TAB(I,5)=ZBF(L5,Z4+Z) 
C 
D     WRITE(6,820) NAME,(KOL(J),J=1,5), 
D    1          (NODE(I),(TAB(I,J),NODE(I),J=1,5),I=1,IXX)
  820 FORMAT('0 ',A8,6X,5(A8,8X),//,(7X,A6,5(F14.9,A2)))
C 
      RETURN
      END 
C 
C
C
      SUBROUTINE XDTPTA(L1) 
C 
      PARAMETER (NDIM=240)
      PARAMETER (LDXDIM=9000) 
C 
      COMMON /SP2NODE/ NXX,DYN,YN(NDIM) 
C 
      COMMON/SP2DTPT/LDX,YD(LDXDIM),WD(LDXDIM)
      COMMON/SP2DCHR/IDSPD1,IDSPD2
      COMMON /DSK/MSKDA1,MSKDA2 
      CHARACTER*8 IDSPD1,IDSPD2*72
C 
      DATA PI/3.14159265359/
C 
C 
C  ** ENTRY XDTPTA(L1)
C 
   10 LDX=L1
      DYL2=(YN(NXX)-YN(1))/(L1-1) 
C 
      DO 14 LB=1,L1 
   14 YD(LB)=YN(1)+(LB-1)*DYL2
C 
C       WD =1.
      DO 16 L=1,LDX 
   16 WD(L)=1.0 
C 
   18 IDSPD1=' XDTPTA ' 
      WRITE(IDSPD2,818) L1
  818 FORMAT('**      LINEAR DATA POINT ARRAY  ',I3,' POINTS WITH', 
     1       '  POINTS ON BOUNDARIES  ')
      RETURN
C 
      END 
C 
C
      SUBROUTINE SPLABORT(KENTRY)
      CHARACTER*7 KENTRY
C 
      WRITE(1,10) KENTRY
   10 FORMAT('0 .... GRIEVOUS ERROR IN ',A7)
      STOP
      END 
C 
C 
C
      SUBROUTINE PRFKMND(A,FKMN,KP1,KP2)
C 
      PARAMETER (NDIM=240)
      PARAMETER (KDTDIM=2)
      PARAMETER (KVDIM=2*KDTDIM)

C 
C     ** ** EXTERNALS MUST BE DECLARED ** **
C 
      EXTERNAL YBF,DYBF,DDYBF 
C 
      COMMON /SP2KOMP/ IDVAR(KVDIM),KPX,IDSPK(KVDIM)
      COMMON /SP2NODE/ NXX,DYN,YN(NDIM) 
C 
      DIMENSION A(KDTDIM,NDIM),FKMN(KVDIM,NDIM) 
C 
      CHARACTER*40 IDSPK,KSCALE*3 
      DATA KSCALE/'   '/
C 
C 
      CALL SPLINTA(A,1,NXX,1) 
      CALL SPLINTD(YBF) 
C 
      DO 12 KP=KP1,KP2
   12 CALL SPLINTO(A,KP,FKMN,KVDIM,NDIM,KP) 
C 
      CALL SPLINTD(DYBF)
      DO 13 KP=KP1,KP2
13    CALL SPLINTO(A,KP,FKMN,KVDIM,NDIM,KP+KDTDIM)

C      DO 20 KP=KP1,KP2 
C   20 WRITE(6,21) KP,IDSPK(KP),(YN(N),FKMN(KP,N),N=1,NDIM) 
C   21 FORMAT('1',15X,'KP=',I3,5X,A40,/,(25X,F6.1,F13.8,/)) 
C 
      RETURN
C 
      END 
C 
C 
C 
      BLOCK DATA
      PARAMETER (KDTDIM=2)
      PARAMETER (KVDIM=2*KDTDIM)

      COMMON /ZETAC/ CPA,CPB,CSTR,CPSTR0,CZ0,CZ1,CZ2,CZ3,CP00 
      COMMON /TITLE/ JOBID,JOBH 
      COMMON /SP2PRIC/ NX,KDCTYP,KBCTYP(2),YN1,YNXX,YDCWL,
     1        YBCWL(2),IDSPC
      COMMON /SP2KOMP/ IDVAR(KVDIM),KPX,IDSPK(KVDIM)
      CHARACTER IDSPK*40,IDSPC*121
      DATA CPA/1012./,CPB/70./,CSTR/5./,CPSTR0/700./,CP00/0.001/
      DATA JOBID/' '/ 
      DATA IDSPC/' '/ 
      DATA KPX/KVDIM/ 
      END 
C 
C                 SPLARGO **  11/20/84
C                             AMM IS SQUARE MATRIX
C 
C     ++ ASSUMED DIM PARAMETERS ... NDIM=240, LDXDIM=9000,
C                               ... KDTDIM=2
C 
C 
      SUBROUTINE SPLARGO
C 
C     ** ONE-D VERSION OF THE EARLIER TWO-D SPLINE ROUTINES **
C          PERIODIC BC IS NOT ALLOWED AT THIS TIME (03/31/79) 
C 
      PARAMETER (NDIM=240)
      PARAMETER (KDTDIM=2)
      PARAMETER (KVDIM=2*KDTDIM)
      PARAMETER (LDXDIM=9000) 
C 
C     ** SECONDARY PARAMETERS FOR COMPILATION 
      PARAMETER (NSDIM=NDIM*NDIM) 
C 
C     ** WORK SPACE IN THE BLANK COMMON **
C 
C        AMM IS SQUARE MATRIX 
C 
      COMMON /BIG/ A(KDTDIM,NDIM),FD(KDTDIM,LDXDIM),FKMN(KVDIM,NDIM), 
     1          AMM(NDIM,NDIM),DMM(NDIM,NDIM) 
C 
C     ** NODAL CONSTANTS
      COMMON /SP2NODE/ NXX,DYN,YN(NDIM) 
C 
C 
C 
C     ** ** SPLARGO BEGINS ** **
C 
C     ZSEC = SECNDS(0.0)
      ZSEC=0
C 
C     ** PREPARE FOR GETAMR.  SPLIMAC ASSUMES /SP2DTPT/ HAS BEEN DEFINED. 
      CALL SPLIMAC
C 
C     ** INITIALIZE THE N-MATRIX ROUTINES 
      CALL PASDEXX(NXX) 
C 
C     ** ZERO MATRICES
      CALL CLRADE 
C 
C     ** MATRIX INVERSION 
      CALL CALDMM 
C 
C     DELTA = SECNDS(ZSEC)
      DELTA=0 
D     WRITE(6,880) DELTA
  880 FORMAT('0 === SPLARGO DONE ===',T91,'ELAPSED TIME=',F10.4)
      RETURN
C 
      END 
C 
C 
C
      SUBROUTINE SPLIMAC
C          ENTRY GETAMR (AMM) 
C 
      PARAMETER (NDIM=240)
      PARAMETER (LDXDIM=9000) 
C 
C     ** DATA POINTS INFORMATION ** 
C 
      COMMON/SP2DTPT/LDX,YD(LDXDIM),WD(LDXDIM)
      COMMON/SP2DCHR/IDSPD1,IDSPD2
C 
C     DUMMY ARGUMENTS ... AMM,IS SQUARE MATRIX  (OCT 79)
      DIMENSION AMM(NDIM,NDIM)
C 
C     LOCAL WORK SPACE
      DIMENSION ROWM(NDIM,4)
C 
C     OTHER CONSTANTS 
      COMMON /SP2PRIC/ NX,KDCTYP,KBCTYP(2),YN1,YNXX,YDCWL,
     1        YBCWL(2),IDSPC
      COMMON /SP2NODE/ NXX,DYN,YN(NDIM) 
C     LOCAL CONSTANTS FOR SPLADDQ 
      COMMON /SP2MACQ/ NXX1,NX1,NX2,NX3 
      CHARACTER*8 IDSPD1,IDSPD2*72,IDSPC*121,KENTRY*7 
C 
      SAVE ROWM,LD,DY0
      DATA KENTRY/'SPLIMAC'/
C 
C     INPUT TO SPLIMAC VIA /SP2DTPT/
C      YD      Y COORDINATES OF DATA POINTS 
C      WD      WEIGHTS ON DATA AT INDIVIDUAL DATA POINTS
C              MUST HAVE BEEN SET BEFORE SPLIMAC IS CALLED
C              ++ TO AVOID POSSIBLE HANG UP IN GETAMR, WD WILL BE 
C                 SET NEGATIVE IF XD OR YD IS OFF LIMITS. 
C              ANY POINTS WITH ZERO OR NEGATIVE WEIGHT WILL NOT BE
C              INCLUDED IN THE ROWM ACCUMULATION BY SPLIMAC/GETAMR
C      LDX     MAX NUMBER OF DATA POINTS.  MUST BE .LE. LDXDIM
C      IDSPD   ID OF THE DATA POINT SET, AND EXPLANATION (80 CHARACTERS)
C 
C     INPUT TO GETAMR 
C 
C     OUT FROM GETAMR 
C      AMM     (N,N) MATRIX AT (M,M) DIAGONAL 
C     +++ +++ ANRTOS0 IS USED TO TRANSFORM REC ARRAYS TO SQUARE ARRAYS +++ +++
C 
C     LOCAL ARRAYS
C      ROWM    WILL BE ACTIVE DURING SPLARGO EXECUTION
C 
C 
C     ** DEFINE /SP2MACQ/ FOR THE USE BY SPLADDQ
      NXX1=NXX+1
      NX1=NX-1
      NX2=NX-2
      NX3=NX-3
C 
C     ** INITIALIZATION FOR GETAMR ** 
C 
C     CHECK DATA POINTS AND WEIGHTS 
C       ONLY DATA INSIDE DOMAIN ARE ACCEPTED
C 
      LDOK=0
      LDOFF=0 
      DO 84 LD=1,LDX
C       IGNORE IF THE WEIGHT IS ALREADY MARKED OFF
      IF(WD(LD).LE.0.) GO TO 84 
C       CHECK FOR OFF-LIMIT POINTS
      IF(YD(LD) .LT. YN1 .OR. YD(LD) .GT. YNXX)  GO TO 82 
C       ACCEPTED POINTS 
      LDOK=LDOK+1 
      GO TO 84
C       OFF THE LIMITS.  WD SET TO NEGATIVE 
   82 WD(LD)=-ABS(WD(LD)) 
      WRITE(1,830) YD(LD) 
  830 FORMAT(' ',10X,'REJECTED DATA POINT...  YD=',F6.2)
      LDOFF=LDOFF+1 
   84 CONTINUE
C 
D     WRITE(6,885) LDX,LDOK,LDOFF 
  885 FORMAT('0 === SPLIMAC DATA POINTS CHECK.     LDX=',I6,
     1    ' (INPUT TOTAL COUNT)',/,36X,'LDOK=',I6,' (ACCEPTED BY '
     2    'SPLIMAC)',/,35X,'LDOFF=',I6,' (REJECTED.  WD REDEFINED)')
C 
C     QUIT, IF LDOK IS ZERO 
      IF(LDOK.LE.0) GO TO 91
C 
C 
C 
      DO 106 NR=1,4 
      DO 106 N =1,NDIM
  106 ROWM(N,NR)=0. 
C 
      DY0=1./DYN
      LD=0
C 
      CALL SPLIQBG
C 
      RETURN
C 
C 
      ENTRY GETAMR(AMM) 
C 
C 
C     CHECK IF UNPROCESSED DATA STILL REMAIN
  110 IF(LD.GE.LDX) GO TO 120 
      LD=LD+1 
      Y=YD(LD)
      W=AMIN1(WD(LD),1.0) 
      IF(W.LE.0.) GO TO 110 
C 
C     ** ROWM ACCUMULATION ** 
      NA1=MAX0(INT((Y-YN1)*DY0),1)
      NAX=MIN0(NA1+3,NXX) 
C 
      DO 115 NA=NA1,NAX 
      WMABNA=YBF(NA,Y)*W
      IF (ABS(WMABNA).LT.1.E-18) WMABNA = 0.
      DO 115 NB=NA,NAX
      NR=NB-NA+1
  115 ROWM(NA,NR)=ROWM(NA,NR)+YBF(NB,Y)*WMABNA
C 
      GO TO 110 
C 
C 
C     NO MORE DATA POINTS LEFT IN THE PRESENT ACCUMULATION RANGE
C 
C     ADD Q 
C 
  120 CALL SPLADDQ(ROWM(1,1)) 
C 
C     TRANSFER RESULTS TO AMM 
C 
      CALL ANRTOS0(ROWM(1,1),AMM(1,1))
C 
      RETURN
C 
C     ** ERROR EXITS ** 
C 
   91 WRITE(1,891)
  891 FORMAT('  === SPLIMAC CALLS EXIT.  NO ACCEPTED DATA POINTS.') 
      CALL SPLABORT(KENTRY)
      END 
C 
C 
      SUBROUTINE SPLIQBG
C 
      PARAMETER (NDIM=240)
C 
      COMMON /SP2PRIC/ NX,KDCTYP,KBCTYP(2),YN1,YNXX,YDCWL,
     1        YBCWL(2),IDSPC
      COMMON /SP2NODE/ NXX,DYN,YN(NDIM) 
      COMMON /SP2BCDC/ DY0,DY1,DY2,CBCBF(4),YSDCPI
C 
      COMMON /SP2QRBG/ QNRBG(5,4,4) 
      CHARACTER*121 IDSPC 
C 
C     ** STRUCTURE OF QNRBG (I,J,K)  (FOR ANY K)
C 
C   I   J= 1         2         3         4
C 
C              QUADRATURE INTEGRAL TERMS               THE ORDER OF DERIVATIVE
C   1     1,1       1,2       1,3       1,4                 IN QUADRATURE 
C   2     2,2       2,3       2,4       2,5 
C   3     3,3       3,4       3,5       3,6              K=1   0TH ORDER
C   4    NX,NX    NX1,NX    NX2,NX    NX3,NX               2   1ST ORDER
C   5   NXX,NXX    NX,NXX   NX1,NXX   NX2,NXX              3   2ND ORDER
C                                                          4   3RD ORDER
C 
C     ORTHO- AND PARA-ELEMENTS ARE PARTIAL QUADRATURES OF BF AND DERIVATIVES
C 
      DIMENSION QO11(4),QO12(4),QO22(4),QP11(4),QP12(4),QP22(4) 
C 
      DATA QO11/0.530357142857143, 0.63750, 2.250, 20.25/,
     1     QO12/0.057589285714286, 0.13125,-1.125, -6.75/,
     2     QO22/0.008928571428571, 0.11250, 0.750,  2.25/,
     3     QP11/0.416517857142857,-0.54375,-1.125,-20.25/,
     4     QP12/0.026785714285714,-0.22500, 0.000,  6.75/,
     5     QP22/0.000446428571429,-0.01875, 0.375, -2.25/ 
C 
C 
C     ** IF PRINT OUTPUT OF /SP2QRBG/, BEFORE BEING MULTIPLIED BY DC FACTORS, 
C        IS DESIRED, CALL SPLIQBP, INSTEAD. 
C 
C     ** ENTRY SPLIQBG, THE NORMAL NO-PRINT ENTRY 
C 
      ASSIGN 40 TO NOPR 
C 
   10 CONTINUE
C 
C     FOR INTERIOR NODES
      DO 12 K=1,4 
      QNRBG(3,1,K)=2.*(QO11(K)+QO22(K)) 
      QNRBG(3,2,K)= 2.*QO12(K)+QP11(K)
      QNRBG(3,3,K)= 2.*QP12(K)
   12 QNRBG(3,4,K)=    QP22(K)
      DO 18 K=1,4 
      DO 16 I=2,4 
   16 QNRBG(I,3,K)=QNRBG(3,3,K) 
      DO 18 I=1,5 
   18 QNRBG(I,4,K)=QNRBG(3,4,K) 
C 
C     FOR BOUNDARY NODES
      ZN11=4.*CBCBF(1)
      ZN12=4.*CBCBF(2)
      ZNX1=4.*CBCBF(3)
      ZNX2=4.*CBCBF(4)
C 
      DO 22 K=1,4 
      QNRBG(1,1,K)=QO11(K)+QO22(K)+ZN11*(2.*QO12(K)+ZN11*QO22(K)) 
      QNRBG(1,2,K)=QP11(K)+QO12(K)+ZN12*(QO12(K)+ZN11*QO22(K))
     1             +ZN11*QP12(K)
      QNRBG(1,3,K)=2.*QP12(K)+ZN11*QP22(K)
      QNRBG(2,1,K)=2.*QO11(K)+QO22(K)+ZN12*(2.*QP12(K)+ZN12*QO22(K))
      QNRBG(2,2,K)=2.*QO12(K)+QP11(K)+ZN12*QP22(K)
      QNRBG(4,1,K)=2.*QO11(K)+QO22(K)+ZNX2*(2.*QP12(K)+ZNX2*QO22(K))
      QNRBG(4,2,K)=2.*QO12(K)+QP11(K)+ZNX2*QP22(K)
      QNRBG(5,1,K)=QO11(K)+QO22(K)+ZNX1*(2.*QO12(K)+ZNX1*QO22(K)) 
      QNRBG(5,2,K)=QP11(K)+QO12(K)+ZNX2*(QO12(K)+ZNX1*QO22(K))
     1             +ZNX1*QP12(K)
      QNRBG(5,3,K)= 2.*QP12(K)+ZNX1*QP22(K) 
   22 CONTINUE
C 
      GO TO NOPR,(40,80)
C 
C     ** PRINT OUTPUT, ONLY IF SPLIQBP WAS CALLED 
   80 CONTINUE
D     WRITE(6,880) KBCTYP(1),KBCTYP(2),ZN11,ZN12,ZNX1,ZNX2
D     WRITE(6,881) (((QNRBG(I,J,K),J=1,4),I=1,5),K=1,4) 
  880 FORMAT('0.../SP2QRBG/... BEFORE DC FACTOR MULTIPLICATION ...',//, 
     1    11X,I3,13X,I3,/,(7X,2F7.3,2X,2F7.3))
  881 FORMAT(//,(5(5X,4F16.10,/),/))
C 
C     ** MULTIPLY DC-FACTORS ** 
C 
   40 DO 42 K=2,4 
      KK=2*(K-1)
      Z2=YSDCPI**KK 
      DO 42 J=1,4 
      DO 42 I=1,5 
      QNRBG(I,J,K)=Z2*QNRBG(I,J,K)
   42 CONTINUE
C 
      RETURN
C 
C 
      ENTRY SPLIQBP 
C 
C     ** THIS ENTRY WILL PRINT /SP2QRBG/, BEFORE DC-FACTOR MULTIPLICATION,
C        BESIDES DOING EVERYTHING SPLIQBG WILL. 
C 
      ASSIGN 80 TO NOPR 
      GO TO 10
C 
      END 
C 
C 
      SUBROUTINE SPLADDQ(PNMR)
C 
C 
      PARAMETER (NDIM=240)
C 
      COMMON /SP2PRIC/ NX,KDCTYP,KBCTYP(2),YN1,YNXX,YDCWL,
     1        YBCWL(2),IDSPC
      COMMON /SP2NODE/ NXX,DYN,YN(NDIM) 
C 
      COMMON /SP2QRBG/ QNRBG(5,4,4) 
      COMMON /SP2MACQ/ NXX1,NX1,NX2,NX3 
      CHARACTER*121 IDSPC,KENTRY*7
C 
      DIMENSION PNMR(NDIM,4)
      DIMENSION QNMG(5,4) 
C 
      DATA KENTRY/'SPLADDQ'/
C 
      DO 8 NR=1,4 
      DO 8 NG=1,5 
    8 QNMG(NG,NR)=0.
C 
C 
C     ** TEST KDCTYP   === TYPE 2 AND 3 ARE PRESENTLY ALLOWED === 
C 
   20 IF(KDCTYP-2) 90,50,30 
C 
C     THIRD ORDER CONSTRAINT ... KDCTYP=3 
   30 CONTINUE
C 
C     SYMMETRIC BODY TERMS
      DO 32 NR=1,4
      DO 32 NG=1,5
   32 QNMG(NG,NR)=QNRBG(NG,NR,4)
C 
C 
C     ** ADD QNMG TO PNMR 
C 
   40 DO 44 NR=1,4
      DO 42 N=1,2 
   42 PNMR(N,NR)=PNMR(N,NR)+QNMG(N,NR)
      DO 44 N=3,NX2 
   44 PNMR(N,NR)=PNMR(N,NR)+QNMG(3,NR)
      DO 46 N=NX1,NXX 
      NRX=NXX1-N
      DO 46 NR=1,NRX
   46 PNMR(N,NR)=PNMR(N,NR)+QNMG(NR+5-NRX,NR) 
C 
      RETURN
C 
C 
C     SECOND ORDER CONSTRAINT ... KDCTYP=2
C 
   50 DO 52 NR=1,4
      DO 52 NG=1,5
   52 QNMG(NG,NR)=QNRBG(NG,NR,3)
   54 CONTINUE
      GO TO 40
C 
C 
   90 WRITE(1,890) KDCTYP 
  890 FORMAT('0 === SPLADDQ DOES NOT ACCEPT   KDCTYP=',I2)
      CALL SPLABORT(KENTRY)
C 
      END 
C 
C
C
      SUBROUTINE ADEXXX 
C 
      PARAMETER (NDIM=240)
      PARAMETER (KDTDIM=2)
      PARAMETER (KVDIM=2*KDTDIM)
      PARAMETER (LDXDIM=9000) 
C 
C     ** SECONDARY PARAMETERS FOR COMPILATION 
      PARAMETER (NSDIM=NDIM*NDIM) 
C 
C     DUMMY ARGUMENTS ARE NOW IN BLANK COMMON FOR VAX-11
      COMMON /BIG/ A(KDTDIM,NDIM),FD(KDTDIM,LDXDIM),FKMN(KVDIM,NDIM), 
     1          AMM(NDIM,NDIM),DMM(NDIM,NDIM) 
C 
      ENTRY CLRADE
C 
      CALL CLRANS(AMM)
      CALL CLRANS(DMM)
      RETURN
C 
C 
C 
      ENTRY CALDMM
C 
      CALL GETAMR(AMM)
C 
      CALL ANEGINV(AMM,DMM) 
C 
      RETURN
C 
      END 
C 
C 
C
      SUBROUTINE PASDEXX(NXX) 
C 
C     ** COLLECTION OF N-MATRIX ROUTINES FOR SPLEFSK
C 
C     ** PASDEXX MUST BE CALLED, ONCE FOR ALL,
C        BEFORE ANY ENTRIES, BELOW, ARE CALLED. 
C 
C 
C  ++++++++ THE LIST OF N-MATRIX ROUTINES ++++++++
C 
C    SUB NAMES   ENTRY NAMES            OPERATIONS
C 
C     PASDEXX    ANRXNSA(AR,BS,CS)      CS=AR*BS+CS 
C                ANRXNS0(AR,BS,CS)      CS=AR*BS
C                ANRTOSA(BR,CS)         CS=   BR+CS 
C                ANRTOS0(BR,CS)         CS=   BR
C                ANRTPSA(BR,CS)         CS=TP(BR)+CS
C                ANRTPS0(BR,CS)         CS=TP(BR) 
C                ANEGINV(BS,CS)         CS=-INV(BS) 
C 
C     PASDESS    ANSXNSA(AS,BS,CS)      CS=AS*BS+CS 
C                ANSXNS0(AS,BS,CS)      CS=AS*BS
C                ANTXNSA(AS,BS,CS)      CS=TP(AS)*BS+CS 
C                ANTXNS0(AS,BS,CS)      CS=TP(AS)*BS
C                ANSTOSA(BS,CS)         CS=   BS+CS 
C                ANSTOS0(BS,CS)         CS=   BS
C                ANSTPS0(AS,CS)         CS=TP(AS) 
C                CLRANS (CS)            CS=0. 
C                 .. VARIANTS IN ADDRESS SHARING .. 
C                ANSAXNS(BS,AS,CS)      CS=BS*AS+CS 
C                ANS0XNS(BS,AS,CS)      CS=BS*AS
C 
C     PASDERR    ANRTOD (BR,CD)         CD=   BR
C                ANRTPD (BR,CD)         CD=TP(BR) 
C                ANRTPR0(BR,CR)         CR=TP(BR) 
C                CLRANR (CR)            CR=0. 
C 
C         AS,BS,CS  NDIM*NDIM SQUARE MATRICES 
C         AR,BR,CR  NDIM*4 RECTANGLE MATRICES 
C               CD  NDIM*7 RECTANGLE MATRICES (DOUBLE R)
C 
C         CS,CR,CD  ARE ALWAYS THE OUTPUT 
C         BS OR BR  MAY SHARE THE SAME ADDRESS WITH CS OR CR, RESPECTIVELY
C         AS OR AR  MUST HAVE ITS OWN ADDRESS DIFFERENT FROM BS,BR,CS,OR CR 
C 
C     ** WORK SPACE IS PROVIDED LOCALLY FOR 
C         ZS(NDIM,NDIM) AND ZD(NDIM,7) IN PASDEXX 
C         ZV(NDIM)                     IN PASDESS 
C 
C     ** DO INDEX CONSTANTS ARE GENERATED IN /PASDECC/ BY PASDEXX 
C 
C 
      PARAMETER (NDIM=240)
C 
C     DUMMY ARGUMENTS 
      DIMENSION BS(NDIM,NDIM),CS(NDIM,NDIM),BSD(NDIM,NDIM)
      DIMENSION AR(NDIM,4)   ,BR(NDIM,4)   ,CSD(NDIM,NDIM)
      DOUBLE PRECISION BSD,CSD,DET
C 
C      LOCAL WORK SPACE 
      COMMON /LITTLE/ ZS(NDIM,NDIM),ZD(NDIM,7)
      CHARACTER*7 KENTRY
C 
      COMMON /PASDECC/ IXX,IX,IX1,IX2 
      DATA KENTRY/'PASDEXX'/
C 
C 
C     ENTRY PASDEXX(NXX)   ... MUST BE CALLED TO DEFINE DO INDEX CONSTANTS ...
C 
      IXX=NXX 
      IX =IXX-1 
      IX1=IX-1
      IX2=IX-2
C 
      IF(NXX.LE.NDIM) RETURN
C 
      WRITE(1,800) NXX,NDIM 
  800 FORMAT('0 === PASDEXX FOUND NXX .GT. NDIM ===',3X,2I6)
      CALL SPLABORT(KENTRY)
C 
C 
      ENTRY ANRXNSA(AR,BS,CS) 
C 
      ASSIGN 16 TO LGO
C 
   10 CALL ANRTOD (AR,ZD) 
      CALL CLRANS (ZS)
C 
      DO 12 I=1,IXX 
      DO 12 J=MAX0(I-3,1),MIN0(I+3,IXX) 
   12 ZS(I,J)=ZD(I,J+4-I) 
C 
      GO TO LGO,(16,18) 
C 
   16 CALL ANSXNSA(ZS,BS,CS)
      RETURN
C 
C 
      ENTRY ANRXNS0(AR,BS,CS) 
C 
      ASSIGN 18 TO LGO
      GO TO 10
C 
   18 CALL ANSXNS0(ZS,BS,CS)
      RETURN
C 
C 
      ENTRY ANRTOSA(BR,CS)
C 
      CALL ANRTOD (BR,ZD) 
C 
   30 DO 32 I=1,IXX 
      DO 32 J=MAX0(I-3,1),MIN0(I+3,IXX) 
   32 CS(I,J)=ZD(I,J+4-I)+CS(I,J) 
      RETURN
C 
C 
      ENTRY ANRTOS0(BR,CS)
C 
      CALL ANRTOD (BR,ZD) 
C 
   35 DO 38 I=1,IXX 
      DO 36 J=1,IXX 
   36 CS(I,J)=0.
      DO 38 J=MAX0(I-3,1),MIN0(I+3,IXX) 
   38 CS(I,J)=ZD(I,J+4-I) 
      RETURN
C 
C 
      ENTRY ANRTPSA(BR,CS)
C 
      ASSIGN 30 TO LGO
C 
   40 CALL ANRTPD (BR,ZD) 
      GO TO LGO,(30,35) 
C 
C 
      ENTRY ANRTPS0(BR,CS)
C 
      ASSIGN 35 TO LGO
      GO TO 40
C 
C 
      ENTRY ANEGINV(BS,CS)
C 
C     COPY BS INTO DOUBLE PRECISION ARRAY 
C 
      DO 50 J=1,IXX 
      DO 50 I=1,IXX 
   50 BSD(I,J)=BS(I,J)
C 
      CALL INVMTX(BSD,NDIM,CSD,NDIM,IXX,DET,ZD,IER) 
C 
      IF(IER.GT.1) GO TO 54 
C 
      DO 52 J=1,IXX 
      DO 52 I=1,IXX 
      IF (ABS(CSD(I,J)).LT.1.E-18) THEN
         CS(I,J) = 0.
         ELSE
         CS(I,J) = CSD(I,J)
         ENDIF
52    CONTINUE
      RETURN
C 
C 
   54 WRITE(1,854) IER,DET
  854 FORMAT('0 === ANEGINV/INVMTX FOUND BS TO BE SINGULAR.   IER=',I2, 
     1    3X, D12.3) 
      KENTRY = 'ANEGINV'
      CALL SPLABORT(KENTRY)
C 
      END 
C 
C
C
      SUBROUTINE PASDESS
C 
      PARAMETER (NDIM=240)
C 
C     DUMMY ARGUMENTS 
      DIMENSION AS(NDIM,NDIM),BS(NDIM,NDIM),CS(NDIM,NDIM) 
C 
C     LOCAL WORK SPACE
      DIMENSION ZV(NDIM)
C 
      COMMON /PASDECC/ IXX,IX,IX1,IX2 
C 
C 
      ENTRY ANSXNSA(AS,BS,CS) 
C 
      DO 16 J=1,IXX 
      DO 12 I=1,IXX 
   12 ZV(I)=CS(I,J) 
      DO 14 K=1,IXX 
      DO 14 I=1,IXX 
   14 ZV(I)=AS(I,K)*BS(K,J)+ZV(I) 
      DO 16 I=1,IXX 
   16 CS(I,J)=ZV(I) 
      RETURN
C 
C 
      ENTRY ANSXNS0(AS,BS,CS) 
C 
      DO 26 J=1,IXX 
      DO 22 I=1,IXX 
   22 ZV(I)=0.
      DO 24 K=1,IXX 
      DO 24 I=1,IXX 
   24 ZV(I)=AS(I,K)*BS(K,J)+ZV(I) 
      DO 26 I=1,IXX 
   26 CS(I,J)=ZV(I) 
      RETURN
C 
C 
      ENTRY ANTXNSA(AS,BS,CS) 
C 
      DO 86 J=1,IXX 
      DO 82 I=1,IXX 
   82 ZV(I)=CS(I,J) 
      DO 84 K=1,IXX 
      DO 84 I=1,IXX 
   84 ZV(I)=AS(K,I)*BS(K,J)+ZV(I) 
      DO 86 I=1,IXX 
   86 CS(I,J)=ZV(I) 
      RETURN
C 
C 
      ENTRY ANTXNS0(AS,BS,CS) 
C 
      DO 96 J=1,IXX 
      DO 92 I=1,IXX 
   92 ZV(I)=0.
      DO 94 K=1,IXX 
      DO 94 I=1,IXX 
   94 ZV(I)=AS(K,I)*BS(K,J)+ZV(I) 
      DO 96 I=1,IXX 
   96 CS(I,J)=ZV(I) 
      RETURN
C 
C 
      ENTRY ANSTOSA(BS,CS)
C 
      DO 32 J=1,IXX 
      DO 32 I=1,IXX 
   32 CS(I,J)=BS(I,J)+CS(I,J) 
      RETURN
C 
C 
      ENTRY ANSTOS0(BS,CS)
C 
      DO 42 J=1,IXX 
      DO 42 I=1,IXX 
   42 CS(I,J)=BS(I,J) 
      RETURN
C 
C 
      ENTRY ANSTPS0(AS,CS)
C 
      DO 45 J=1,IXX 
      DO 45 I=1,IXX 
   45 CS(I,J)=AS(J,I) 
      RETURN
C 
C 
      ENTRY CLRANS (CS) 
C 
      DO 52 J=1,IXX 
      DO 52 I=1,IXX 
   52 CS(I,J)=0.
      RETURN
C 
C 
C     ** VARIATION IN THE ADDRESS SHARING **
C        THE FIRST (LEFT) MATRIX MAY BE OF THE SAME ADDRESS WITH THE RESULT, CS.
C 
C 
      ENTRY ANSAXNS(BS,AS,CS) 
C 
      DO 66 I=1,IXX 
      DO 62 J=1,IXX 
   62 ZV(J)=CS(I,J) 
      DO 64 K=1,IXX 
      DO 64 J=1,IXX 
   64 ZV(J)=BS(I,K)*AS(K,J)+ZV(J) 
      DO 66 J=1,IXX 
   66 CS(I,J)=ZV(J) 
      RETURN
C 
C 
      ENTRY ANS0XNS(BS,AS,CS) 
C 
      DO 76 I=1,IXX 
      DO 72 J=1,IXX 
   72 ZV(J)=0.
      DO 74 K=1,IXX 
      DO 74 J=1,IXX 
   74 ZV(J)=BS(I,K)*AS(K,J)+ZV(J) 
      DO 76 J=1,IXX 
   76 CS(I,J)=ZV(J) 
      RETURN
C 
      END 
C 
      SUBROUTINE PASDERR
C 
      PARAMETER (NDIM=240)
C 
C     DUMMY ARGUMENTS 
      DIMENSION BR(NDIM,4),CR(NDIM,4),CD(NDIM,7)
C 
      COMMON /PASDECC/ IXX,IX,IX1,IX2 
C 
C 
      ENTRY ANRTOD (BR,CD)
C 
      ASSIGN 16 TO LGO
C 
C      EXPAND BR TO DOUBLE RECTANGLE ARRAY, CD(NDIM,7)
C     ** THE UN-USED CORNERS OF CD ARE NOT CLEARED. 
C 
   10 DO 12 J=4,7 
      DO 12 I=1,IXX 
   12 CD(I,J)=BR(I,J-3) 
C 
      DO 14 J=1,3 
      DO 14 I=5-J,IXX 
   14 CD(I,J)=BR(I+J-4,5-J) 
C 
      GO TO LGO,(16,18) 
C 
C     ADD OR SUBTRACT ANTISYMMETRIC TERMS 
   16 CD(1,5)=CD(1,5)+CD(IXX,6) 
      CD(2,3)=CD(2,3)-CD(IXX,6) 
      CD(IXX-1,5)=CD(IXX-1,5)+CD(IXX-1,7) 
      CD(IXX  ,3)=CD(IXX  ,3)-CD(IXX-1,7) 
      RETURN
C 
C 
      ENTRY ANRTPD (BR,CD)
C 
      ASSIGN 18 TO LGO
      GO TO 10
C 
C      REVERSE THE SIGN OF ANTISYMMETRIC TERMS (TRANSPOSE)
   18 CD(IXX  ,6)=-CD(IXX  ,6)
      CD(IXX-1,7)=-CD(IXX-1,7)
      GO TO 16
C 
C 
      ENTRY ANRTPR0(BR,CR)
C 
      DO 22 J=1,4 
      DO 22 I=1,IXX 
   22 CR(I,J)=BR(I,J) 
C 
      CR(IXX  ,3)=-CR(IXX  ,3)
      CR(IXX-1,4)=-CR(IXX-1,4)
      RETURN
C 
C 
      ENTRY CLRANR (CR) 
C 
      DO 32 J=1,4 
      DO 32 I=1,IXX 
   32 CR(I,J)=0.
      RETURN
C 
      END 
C 
C                 SPLIMBO **  11/20/84
C     ++ ASSUMED DIM PARAMETERS ... NDIM=240, LDXDIM=9000,
C                               ... KDTDIM=2
C     ++ FOR FULL PRINTING OF SPLIMBO LOG, SET "LIMBOPR=0"="LIMBOPR=1"
C 
C 
      SUBROUTINE SPLIMBO(KP1,KP2) 
C 
C     ** NODAL SPLINE AMPLITUDES A(K,N) FROM DATA FD(K,L)** 
C 
C 
      PARAMETER (NDIM=240)
      PARAMETER (LDXDIM=9000) 
      PARAMETER (KDTDIM=2)
      PARAMETER (KVDIM=2*KDTDIM)
C 
C     ** SECONDARY PARAMETERS FOR COMPILATION 
C 
      PARAMETER (NSDIM=NDIM*NDIM) 
      PARAMETER (LIMBOPR=0) 
C 
C     ** DUMMY ARGUMENTS
C 
C      DIMENSION FD(KDTDIM,LDXDIM),A(KDTDIM,NDIM) 
C 
C 
C              FD(K,L)   K-TH COMPONENT OF FIELD DATA AT LOCATION L 
C              KP1       FIRST COMPONENT K TO BE PROCESSED
C              KP2       LAST  COMPONENT K TO BE PROCESSED
C 
C     OUTPUT   A(K,N)  NODAL AMPLITUDE OF K-TH COMPONENT AT NODE N
C 
C 
C 
C        $$ SPLOVER WAS CALLED TO DEFINE CONSTANTS IN SEVERAL COMMONS 
C 
C        $$ DATA POINTS WERE DEFINED IN /SP2DTPT/ INCLUDING WEIGHTS 
C 
C     ** EXPLICITLY REQUIRED COMMONS
C 
      COMMON /BIG/ A(KDTDIM,NDIM),FD(KDTDIM,LDXDIM),FKMN(KVDIM,NDIM), 
     1          AMM(NDIM,NDIM),DMM(NDIM,NDIM) 
C 
      COMMON /SP2NODE/ NXX,DYN,YN(NDIM) 
      COMMON/SP2DTPT/LDX,YD(LDXDIM),WD(LDXDIM)
      COMMON/SP2DCHR/IDSPD1,IDSPD2
C 
      CHARACTER*8 IDSPD1,IDSPD2*72,KENTRY*7 
      DIMENSION ZV(NDIM,KDTDIM) 
      DATA KENTRY/'SPLIMBO'/
C 
C 
C     ** LOCAL WORK SPACE 
C 
C 
C 
C 
      KA=MIN0(KP1,KP2)
      KB=MAX0(KP1,KP2)
      IF(KB.GT.KVDIM) GO TO 90
C 
      DO 8 K=KA,KB
      DO 8 N=1,NXX
    8 A(K,N)=0. 
C 
C 
C     ** CONSTANTS FOR CONVENIENCE
      YN1=YN(1) 
      DY0=1./DYN
C 
C 
C     ** CONVERSION OF DISCRETE FIELD DATA TO NODAL FORCING **
C 
      DO18 L=1,LDX
C 
      IF(WD(L).LE.0) GO TO 18 
      NA=MAX0(INT((YD(L)-YN1)*DY0),1) 
      NB=MIN0(NA+3,NXX) 
C 
      DO 15 N=NA,NB 
      WMN=YBF(N,(YD(L)))*AMIN1(WD(L),1.0) 
      DO 15 K=KA,KB 
   15 A(K,N)=A(K,N)+WMN*FD(K,L) 
C 
   18 CONTINUE
C 
C 
C     ** CONVERSION OF NODAL FORCING TO NODAL AMPLITUDES ** 
C       +++ WE RETURN TO THE RIGHT-HANDED FORCING +++ 
C       +++ THE LEFT-HANDED CONVENTION IN PLIB SPLEND IS NOT USED +++ 
C 
C 
      DO 30 K=KA,KB 
      DO 20 I=1,NXX 
   20 ZV(I,K)=0.0 
      DO 30 I=1,NXX 
      DO 30 J=1,NXX 
      ZV(I,K)=ZV(I,K)+DMM(I,J)*A(K,J) 
      IF (ABS(ZV(I,K)).LT.1.E-18) ZV(I,K) = 0.
   30 CONTINUE
      DO 40 K=KA,KB 
      DO 40 N=1,NXX 
   40 A(K,N)=ZV(N,K)
C 
C 
      RETURN
C 
C     ** ERROR EXIT 
C 
   90 WRITE(1,890) KP1,KP2,KVDIM
  890 FORMAT('0 === SPLIMBO FOUND THE COMPONENT RANGE  K=',I2,',',I2, 
     1       '  EXCEEDED  KVDIM=',I2) 
      CALL SPLABORT(KENTRY)
      END 
C 
C 
      SUBROUTINE SPLINT (A,KPA,F,KPFDIM,NFDIM,KPF)
C 
C 
C 
      PARAMETER (NDIM=240)
      PARAMETER (KDTDIM=2)
      PARAMETER (KVDIM=2*KDTDIM)
C
C      EXTERNAL ZYBF     !Added 11/20/95 for HP
C 
C     ** DUMMY ARGUMENTS
C 
C     ** INPUT ... SPLINE AMPLITUDE 
      DIMENSION A(KDTDIM,NDIM)
C 
C     ** OUTPUT ..
      DIMENSION F(KPFDIM,NFDIM) 
C 
C     ** PARAMETER FOR THE LOCAL USE ... MAXIMUM DIVISIONS ALLOWED
      PARAMETER (JXBFDIM=8) 
C 
C     ** LOCAL ARRAYS ... BF TABLES 
      DIMENSION YB1(4,JXBFDIM),YBN(4,JXBFDIM),YBX(4,JXBFDIM),YBXX(4)
C     ** LOCAL WORK SPACE 
      CHARACTER*7 KENTRY
C 
      COMMON /SP2NODE/ NXX,DYN,YN(NDIM) 
      SAVE NX,NX1,NB,JX,NAIL,NBIL1,NBIL2,NFDIM0,JXBF0,YB1,YBN,YBX,YBXX
      SAVE LSFIN,LSFOUT,LGON1,LGONX 
      DATA KENTRY/'SPLINT'/ 
C 
          GO TO 80
C 
      ENTRY SPLINTA(A,N1,N2,N3) 
      KENTRY = 'SPLINTA'
C 
C 
C       LOCALLY DEFINED 
      NX=NXX-1
      NX1=NX-1
      IF(NXX.GT.NDIM) GO TO 900 
C 
C     ** SET CONSTANTS FOR DEFINING THE AREA AND SUBDIVISIONS 
C 
      NA=MIN0(N1,N2)
      NB=MAX0(N1,N2)
      JX=MAX0(N3,1) 
      IF(NA.LT.1)   GO TO 901 
      IF(NB.GT.NXX) GO TO 901 
C 
      IF(JX.GT.JXBFDIM) GO TO 905 
C 
C 
      NAOL=MAX0(NA-1,1) 
      NBOL=MIN0(NB+1,NXX) 
      NAIL=MAX0(NA,2) 
      NBIL1=MIN0(NB,NX1)
      NBIL2=MIN0(NB-1,NX1)
C 
C     ** MINIMUM DIMENSIONS REQUIRED OF THE OUTPUT FIELD
C 
      NFDIM0=JX*(NB-NA)+1 
C 
C     ** SET THE IN AND OUT SWITCHES TO THE NORMAL POSITIONS
C 
      ASSIGN  83 TO LSFIN 
      ASSIGN 320 TO LSFOUT
C 
C     ** ROUTING DECISIONS
C 
      ASSIGN 120 TO LGON1 
      ASSIGN 150 TO LGONX 
C 
      IF(NA.EQ.1)   ASSIGN 110 TO LGON1 
      IF(NB.EQ.NXX) ASSIGN 126 TO LGONX 
      IF(NB.EQ.NX)  ASSIGN 128 TO LGONX 
C 
      RETURN
C 
C 
C     ** INPUT ... THE SPLINE BASIS FUNCTIONS 
C 
C     ZXBF     XBF,DXBF, OR DDXBF 
C     ZYBF     YBF,DYBF, OR DDYBF 
C 
C     ** ** ** WARNING ** ** ** 
C 
C       THERE IS NO DEFAULT SETTING OF ZYBF.
C         SPLINTD MUST BE CALLED, AT LEAST ONCE, TO DEFINE THE BF TABLES
C 
C     ** WARNING ** 
C       DO NOT SHORTCUT SPLINTD WITHIN THIS SUBROUTINE ALTHOUGH IT MAY
C       SEEM TO BE UNNECESSARY WHEN N3 IS UNCHANGED.
C 
C       SINCE ZYBF MAY BE CHANGED, ONLY THE CALLING ROUTINE 
C       CAN DETERMINE THE NEED OF CALLING SPLINTD ANEW. 
C 
C       SIMILARLY, THE TEST OF N3 BY SPLINTO IS NO GUARANTEE FOR
C       CORRECT BF TABLES, IF THERE WERE OMISSION OF CALLING SPLINTD. 
C 
C 
      ENTRY SPLINTD(ZYBF) 
C 
C     ** SAVE THE CURRENT JX(=N3) 
      JXBF0=JX
C 
      DYNJ=DYN/FLOAT(JX)
C 
      DO 64 J=1,JX
      YNJ=(J-1)*DYNJ
C 
      DO 61 ND=0,2
   61 YB1(2+ND,J)=ZYBF(1+ND,YN(1)+YNJ)
C 
      DO 62 ND=-1,2 
   62 YBN(2+ND,J)=ZYBF(3+ND,YN(3)+YNJ)
C 
      DO 63 ND=-1,1 
   63 YBX(2+ND,J)=ZYBF(NX+ND,YN(NX)+YNJ)
   64 CONTINUE
C 
      DO 65 ND=-1,0 
   65 YBXX(2+ND)=ZYBF(NXX+ND,YN(NXX)) 
C 
      RETURN
C 
C 
C 
C 
C     ** INPUT ARGUMENTS
C 
C       ISFIN=-1   SIGN OF INPUT F IS REVERSED (NOT CLEARED)
C            = 0   INPUT F IS CLEARED TO ZERO 
C            = 1   INPUT F IS PASSED AS IS, TO BE ACCUMULATED UPON
C 
C       ISFAC= 1   POSITIVE ACCUMULATION OF THE NEW F ON THE INPUT F
C            =-1   NEGATIVE ACCUMULATION OF THE NEW F ON THE INPUT F
C 
C       NOTE ... LSFIN,LSFOUT DO NOT CORRESPOND TO ISFIN,ISFAC, INDIVIDUALLY. 
C                SINCE LSFOUT FLIPS THE SIGN OF THE FINAL OUTPUT FIELD, 
C                LSFIN MUST BE SET KNOWING WHAT LSFOUT WILL DO AT THE END 
C 
C       ++ SPLINTO, WHEN EXECUTED, ALWAYS RESETS THE INTERNAL SWITCHES
C          TO  ISFIN=0 AND ISFAC=1.  SPLINTS MUST BE CALLED EVERYTIME 
C         WHEN DIFFERENT SETTINGS ARE DESIRED 
C 
C 
      ENTRY SPLINTS(ISFIN,ISFAC)
C 
      IF(ISFAC) 71,72,72
C 
C       FLIP THE SIGN OF THE OUTPUT FIELD 
   71 ASSIGN 310 TO LSFOUT
      IF(ISFIN) 75,74,73
C 
C       PASS THE OUTPUT FIELD AS IS (THE NORMAL OUT SWITCH POSITION)
   72 ASSIGN 320 TO LSFOUT
      IF(ISFIN) 73,74,75
C 
C       FLIP THE SIGN OF THE INPUT FIELD
   73 ASSIGN 81 TO LSFIN
      GO TO 78
C 
C       CLEAR THE INPUT FIELD (THE NORMAL IN SWITCH POSITION) 
   74 ASSIGN 83 TO LSFIN
      GO TO 78
C 
C       PASS THE INPUT FIELD AS IS
   75 ASSIGN 85 TO LSFIN
C 
   78 RETURN
C 
C 
      ENTRY SPLINTO(A,KPA,F,KPFDIM,NFDIM,KPF) 
      KENTRY = 'SPLINTO'
C 
   80     CONTINUE
C 
      KA=KPA
      KF=KPF
C 
C     ** CHECK THE INPUT ARGUMENTS
      IF(KPA.GT.KDTDIM) GO TO 902 
      IF(KPF.GT.KPFDIM) GO TO 910 
      IF(NFDIM.LT.NFDIM0) GO TO 910 
C 
      IF(JX.NE.JXBF0) GO TO 920 
C 
C     ** THE IN SWITCH
C 
      GO TO LSFIN,(81,83,85)
C 
   81 DO 82 N=1,NFDIM 
   82 F(KF,N)=-F(KF,N)
      GO TO 85
   83 DO 84 N=1,NFDIM 
   84 F(KF,N)=0.
C 
C       RESET THE IN SWITCH TO THE NORMAL POSITION ... CLEAR
C 
   85 ASSIGN 83 TO LSFIN
C 
  100 DO 300 J=1,JX 
C 
      IF(J-2) 101,102,103 
C 
C       ON THE N-NODES (J=1)
  101 NDX=1 
      NBIL=NBIL1
      NBJ=NB
      GO TO 103 
C 
C       OFF-NODE FINE MSESH POINTS
  102 NDX=2 
      NBIL=NBIL2
      NBJ=NB-1
C 
  103 CONTINUE
C 
C 
      GO TO LGON1,(110,120) 
C 
C       AT N=1. (FOR NA=1 ONLY) 
  110 DO 112 ND=0,NDX 
  112 F(KF,1)=A(KA,1+ND)*YB1(2+ND,J)+F(KF,1)
C 
C       AT ALL N BETWEEN NAIL AND NBIL, BUT NBIL IS.LE. NX1 
C         NOTE ... THE SECOND DO WILL NOT EXECUTE IF NAIL .GT. NBIL 
C 
  120 DO 122 ND=-1,NDX
      DO 122 N=NAIL,NBIL
  122 F(KF,N)=A(KA,N+ND)*YBN(2+ND,J)+F(KF,N)
C 
C      CHECK THE LEFT LIMIT OF ACCUMULATION 
      GO TO LGONX,(126,128,150) 
  126 IF(J.EQ.1) GO TO 130
      GO TO 140 
  128 IF(J.EQ.1) GO TO 140
      GO TO 150 
C 
C       ONLY FOR J=1 AND NB=NXX 
  130 DO 132 ND=-1,0
  132 F(KF,NXX)=A(KA,NXX+ND)*YBXX(2+ND)+F(KF,NXX) 
C 
C       FOR ALL J IF NB=NXX(ENTERS FROM 132)
C       ONLY FOR J=1 IF NB=NX 
C 
  140 DO 142 ND=-1,1
  142 F(KF,NX)=A(KA,NX+ND)*YBX(2+ND,J)+F(KF,NX) 
C 
  150 CONTINUE
C 
C 
  300 CONTINUE
C 
C 
C     ** THE OUT SWITCH 
C 
      GO TO LSFOUT,(310,320)
C 
  310 DO 312 N=1,NFDIM
  312 F(KF,N)=-F(KF,N)
C 
C       RESET THE OUT SWITCH TO THE NORMAL POSITION (PASS AS IS)
C 
      ASSIGN 320 TO LSFOUT
  320 CONTINUE
C 
      RETURN
C 
C 
C     ** ERROR EXITS
C 
  900 WRITE(1,8900) NXX,NDIM
 8900 FORMAT('0 === SPLINTA FOUND NXX TOO LARGE FOR NDIM' 
     1   ,5X,I5,5X,I5)
      GO TO 990 
C 
  901 WRITE(1,8901) N1,N2,NXX 
 8901 FORMAT('0 === SPLINTA INPUT ARGUMENTS N1,N2 ARE IN ERROR' 
     1    ,5X,I5,3X,I5,5X,'NXX IS',3X,I5) 
      GO TO 990 
C 
  902 WRITE(1,8902) KPA,KDTDIM
 8902 FORMAT('0 === SPLINTO ARGUMENT KPA IS TOO LARGE FOR KDTDIM',5X,2I3
     1) 
      GO TO 990 
C 
  905 WRITE(1,8905) N3,JXBFDIM
 8905 FORMAT('0 === SPLINTA INPUT ARGUMENTS N3 ARE TOO LARGE FOR '
     1    ,'THE LOCAL PARAMETER JXBFDIM',5X,I3,5X,I3) 
      GO TO 990 
C 
  910 WRITE(1,8910) KPFDIM,KPF,NFDIM,NFDIM0 
 8910 FORMAT('0 === SPLINTO ARGUMENTS KPFDIM,NFDIM ARE TOO SMALL '
     1   ,'FOR THE REQUIRED KPF,NFDIM0', 5X,I3,I5,5X,I3,I5) 
      GO TO 990 
C 
  920 WRITE(1,8920) JX,JXBF0
 8920 FORMAT('0 === SPLINTO FOUND N3 HAVE BEEN CHANGED SINCE THE '
     1   ,'LAST CALL TO SPLINTD WAS MADE',5X,I3,5X,I3)
      GO TO 990 
C 
  990 CALL SPLABORT(KENTRY)
      END 
C 
C 
      SUBROUTINE SPOTINO(Y,A,K1,K2,F) 
C 
      PARAMETER (NDIM=240)
      PARAMETER (KDTDIM=2)
      PARAMETER (KVDIM=2*KDTDIM)
C 
C 
C     ** ** SPOT INVERSE TRANSFORMATION OF SPLINE AMPLITUDES TO F AT X,Y
C 
C       X,Y        THE SPOT COORDINATES AT WHICH OUTPUT F IS DESIRED
C       A(K,N)     KTH COMPONET OF NODAL AMPLITUDE AT N 
C       K1,K2      THE  RANGE OF COMPONENTS K OF THE AMPLITUDE (NOT OF OUTPUT F)
C       F(KF)      THE SPOT OUTPUT F.  KF=1 FOR K=K1
C                                      KF=2 FOR K=K1+1,  TILL 
C                                      KF=K2-K1+1 FOR K=K2
C 
C     ** DUMMY ARGUMENTS
      DIMENSION A(KDTDIM,NDIM),F(*) 
C 
      COMMON /SP2NODE/ NXX,DYN,YN(NDIM) 
      EQUIVALENCE (YN(1),Y1)
C 
C     ** LOCAL WORK SPACE 
      DIMENSION BN(4),ZN(4) 
C 
C 
      N0=MAX0(INT((Y-Y1)/DYN),1)-1
      N4=MIN0(NXX-N0,4) 
C 
      DO 13 N=1,N4
   13 BN(N)=YBF(N+N0,Y) 
C 
      KF=0
      DO 26 K=K1,K2 
      KF=KF+1 
C 
      DO 22 N=1,4 
   22 ZN(N)=0.
C 
      DO 24 N=1,N4
   24 ZN(N)=ZN(N)+A(K,N+N0)*BN(N) 
C 
   26 F(KF)=ZN(1)+ZN(2)+ZN(3)+ZN(4) 
C 
      RETURN
      END 
C 
C 
      SUBROUTINE SPOTDER(Y,A,K1,K2,F) 
C 
      PARAMETER (NDIM=240)
      PARAMETER (KDTDIM=2)
      PARAMETER (KVDIM=2*KDTDIM)
C 
C 
C     ** ** SPOT INVERSE TRANSFORMATION OF SPLINE AMPLITUDES TO F AT X,Y
C 
C       X,Y        THE SPOT COORDINATES AT WHICH OUTPUT F IS DESIRED
C       A(K,N)     KTH COMPONET OF NODAL AMPLITUDE AT N 
C       K1,K2      THE  RANGE OF COMPONENTS K OF THE AMPLITUDE (NOT OF OUTPUT F)
C       F(KF)      THE SPOT OUTPUT F.  KF=1 FOR K=K1
C                                      KF=2 FOR K=K1+1,  TILL 
C                                      KF=K2-K1+1 FOR K=K2
C 
C     ** DUMMY ARGUMENTS
      DIMENSION A(KDTDIM,NDIM),F(*) 
C 
      COMMON /SP2NODE/ NXX,DYN,YN(NDIM) 
      EQUIVALENCE (YN(1),Y1)
C 
C     ** LOCAL WORK SPACE 
      DIMENSION BN(4),ZN(4) 
C 
C 
      N0=MAX0(INT((Y-Y1)/DYN),1)-1
      N4=MIN0(NXX-N0,4) 
C 
      DO 13 N=1,N4
   13 BN(N)=DYBF(N+N0,Y)
C 
      KF=0
      DO 26 K=K1,K2 
      KF=KF+1 
C 
      DO 22 N=1,4 
   22 ZN(N)=0.
C 
      DO 24 N=1,N4
   24 ZN(N)=ZN(N)+A(K,N+N0)*BN(N) 
C 
   26 F(KF)=ZN(1)+ZN(2)+ZN(3)+ZN(4) 
C 
      RETURN
      END 
C 
C 
      SUBROUTINE SPOTIKA(Y,A,K1,K2,F) 
C 
      PARAMETER (NDIM=240)
      PARAMETER (KDTDIM=2)
      PARAMETER (KVDIM=2*KDTDIM)
C 
C 
C     ** ** SPOT INVERSE TRANSFORMATION OF SPLINE AMPLITUDES TO F AT X,Y
C 
C     ++ ++ THIS VERSION, SPOTIKA, IS DIFFERENT FROM SPOTINO IN SPLIMBO,
C           IN THAT EXECUTION IS VECTORIZED WITH RESPECT TO KP. 
C 
C       Y          THE SPOT COORDINATES AT WHICH OUTPUT F IS DESIRED
C       A(K,N)     KTH COMPONET OF NODAL AMPLITUDE AT N 
C       K1,K2      THE  RANGE OF COMPONENTS K OF THE AMPLITUDE (NOT OF OUTPUT F)
C       F(KF)      THE SPOT OUTPUT F.  KF=1 FOR K=K1
C                                      KF=2 FOR K=K1+1,  TILL 
C                                      KF=K2-K1+1 FOR K=K2
C 
C     ** DUMMY ARGUMENTS
      DIMENSION A(KDTDIM,NDIM),F(*) 
C 
      COMMON /SP2NODE/ NXX,DYN,YN(NDIM) 
      EQUIVALENCE (YN(1),Y1)
C 
C     ** LOCAL WORK SPACE 
      DIMENSION BN(4),ZN(4,KDTDIM)
C 
C 
      N0=MAX0(INT((Y-Y1)/DYN),1)-1
      N4=MIN0(NXX-N0,4) 
C 
      DO 13 N=1,N4
   13 BN(N)=YBF(N+N0,Y) 
C 
      K0=K1-1 
      KFX=K2-K0 
C 
      DO 22 N=1,4 
      DO 22 KF=1,KFX
   22 ZN(N,KF)=0. 
C 
      DO 24 N=1,N4
      DO 24 KF=1,KFX
   24 ZN(N,KF)=ZN(N,KF)+A(KF+K0,N+N0)*BN(N) 
C 
      DO 26 KF=1,KFX
   26 F(KF)=ZN(1,KF)+ZN(2,KF)+ZN(3,KF)+ZN(4,KF) 
C 
      RETURN
      END 
C 
C 
C
      SUBROUTINE INVMTX (A,NA,V,NV,N,D,IP,IER)
C 
C 
C DIMENSION OF           A(NA,N),V(NV,N),IP(2*N)
C ARGUMENTS 
C 
C LATEST REVISION        AUGUST 1977
C 
C PURPOSE                INVMTX CALCULATES THE INVERSE OF THE MATRIX A, 
C                        OF ORDER N, USING GAUSSIAN ELIMINATION WITH
C                        FULL PIVOTING. 
C 
C ACCESS CARD            *FORTRAN,S=ULIB,N=INVMTX 
C 
C USAGE                  CALL INVMTX (A,NA,V,NV,N,D,IP,IER) 
C 
C ARGUMENTS 
C 
C ON INPUT               A
C                          A TWO-DIMENSIONAL VARIABLE WITH ROW, (FIRST) 
C                          DIMENSION NA AND COLUMN (SECOND) DIMENSION 
C                          .GE. N.  ON INPUT, A CONTAINS THE ELEMENTS OF
C                          THE N BY N MATRIX TO BE INVERTED.
C 
C                        NA 
C                          AN INTEGER INPUT VARIABLE SET EQUAL TO THE 
C                          ROW (FIRST) DIMENSION OF A AS DECLARED IN THE
C                          CALLING PROGRAM.  NA MUST BE .GE. N. 
C 
C                        NV 
C                          AN INTEGER INPUT VARIABLE SET EQUAL TO THE 
C                          ROW (FIRST) DIMENSION OF V AS DECLARED IN THE
C                          CALLING PROGRAM.  NV MUST BE .GE. N. 
C 
C                        N
C                          AN INTEGER INPUT VARIABLE SET EQUAL TO THE 
C                          ORDER OF THE MATRIX TO BE INVERTED.
C 
C                        IP 
C                          AN INTEGER ARRAY USED INTERNALLY FOR WORKING 
C                          STORAGE.  IT MUST HAVE DIMENSION AT LEAST
C                          2*N. 
C 
C ON OUTPUT              V
C                          A TWO-DIMENSIONAL VARIABLE WITH ROW (FIRST)
C                          DIMENSION NV AND COLUMN (SECOND) DIMENSION 
C                          .GE. N.  ON OUTPUT, V CONTAINS THE INVERSE OF
C                          A. 
C 
C                          IF THE MATRIX A NEED NOT BE SAVED, THEN THE
C                          MATRIX A MAY BE USED AS THE INVERSE V IN THE 
C                          CALL; I.E.,
C 
C                            CALL INVMTX(A,NA,A,NA,N,D,IP,IER)
C 
C                          THEN, ON OUTPUT, A WILL CONTAIN THE INVERSE
C                          MATRIX AND THE ORIGINAL A WILL HAVE BEEN 
C                          DESTROYED. 
C 
C                        D
C                          A REAL VARIABLE WHICH, ON OUTPUT, CONTAINS 
C                          THE DETERMINANT OF A.  IF, DURING THE
C                          COMPUTATION OF THE INVERSE, THE DETERMINANT
C                          BECOMES LARGER THAN THE LARGEST FLOATING 
C                          POINT NUMBER THE MACHINE CAN REPRESENT, THE
C                          DETERMINANT CALCULATION IS ABANDONED, A
C                          MESSAGE IS PRINTED, AND A MEANINGLESS NUMBER 
C                          IS RETURNED IN D.
C 
C                        IER
C                          AN INTEGER ERROR FLAG
C                            =  1  IF THE DETERMINANT IS LARGER THAN THE
C                                  LARGEST FLOATING POINT NUMBER THAT 
C                                  THE MACHINE CAN REPRESENT. 
C                            = 33  IF THE MATRIX IS NUMERICALLY 
C                                  SINGULAR.
C                            = 34  IF N .LT. 1  . 
C                            = 35  IF NA .LT. N.
C                            = 36  IF NV .LT. N.
C 
C ENTRY POINTS           INVMTX, IERINV 
C 
C COMMON BLOCKS          NONE 
C 
C I/O                    AN NCAR RESIDENT ROUTINE, ULIBER, IS USED TO 
C                        PRINT ERROR MESSAGES.
C 
C PRECISION              DOUBLE (6/22/84 FOR VAX-11 USAGE)
C 
C REQUIRED ULIB          NONE 
C ROUTINES
C 
C SPECIALISTS            JO WALSH AND BEN DOMENICO, NCAR, 
C                        BOULDER, COLORADO 80307
C 
C LANGUAGE               FORTRAN
C 
C HISTORY                REVISED AND STANDARDIZED FOR NSSL BY JO WALSH, 
C                        APRIL 1974; CERTIFIED BY BEN DOMENICO, 
C                        APRIL 1976.
C 
C ALGORITHM              THE SUBROUTINE INVMTX SOLVES THE MATRIX
C                        EQUATION 
C                            A*V = I
C                        FOR THE MATRIX V (I.E., A**(-1)), WHERE I IS 
C                        THE IDENTITY MATRIX.  THE METHOD USED IS 
C                        GAUSSIAN ELIMINATION WITH FULL PIVOTING.  THE
C                        MATRIX A IS DECOMPOSED INTO THE PRODUCT OF A 
C                        LOWER TRIANGULAR MATRIX AND AN UPPER TRIANGULAR
C                        MATRIX 
C                            A = L*U
C                        THUS 
C                                A*V = I
C                            (L*U)*V = I
C                                U*V = L**(-1). 
C 
C                        THIS UPPER TRIANGULAR MATRIX EQUATION IS SOLVED
C                        FOR THE COLUMNS OF V USING BACK SUBSTITUTION.
C 
C                        IF AT ANY POINT, AFTER PIVOTING IS DONE, A 
C                        PIVOT ELEMENT IS NUMERICALLY ZERO RELATIVE TO
C                        THE LARGEST ELEMENT OF THE MATRIX, THE MATRIX A
C                        IS DECLARED SINGULAR AND THE SUBROUTINE
C                        TERMINATES.
C 
C SPACE REQUIRED         535 (OCTAL) = 349 (DECIMAL)
C 
C ACCURACY               FOR ANY INPUT MATRIX, A, 
C 
C                          MAX1((NORM(A*V(I)-E(I))/(NORM(A)*NORM(V(I))) 
C                              .LE. TOL*N*EPS 
C 
C                        WHERE V(I) IS THE ITH COLUMN OF THE MATRIX V,
C                        E(I) IS THE ITH COLUMN OF THE IDENTITY MATRIX, 
C                        EPS IS THE LARGEST FLOATING POINT NUMBER SUCH
C                        THAT 1.+EPS = 1., AND ALL THE NORMS ARE
C                        L-INFINITY NORMS.  TOL IS A TOLERANCE FACTOR 
C                        SET TO 100 FOR TESTING PURPOSES. 
C 
C                        FOR A SPECIFIC INPUT MATRIX A, AN UPPER BOUND
C                        ON THE RELATIVE ERROR IN THE COMPUTED INVERSE V
C                        CAN BE CALCULATED IN TERMS OF THE RESIDUALS
C                        MATRIX R AS FOLLOWS: 
C 
C                            ((NORM(E))/(NORM(A**(-1)))) .LE. NORM(R) 
C 
C                        WHERE V = A**(-1)+E (THE EXACT INVERSE,
C                        A**(-1), PLUS AN ERROR MATRIX, E).  THUS, THE
C                        RELATIVE ERROR IN V, AS MEASURED BY NORM 
C                        (E)/NORM (A**(-1)), IS BOUNDED BY THE NORM OF
C                        THE RESIDUAL MATRIX WHICH CAN BE COMPUTED AS 
C                        R = A*V**(-1). 
C 
C TIMING                 THE TIME REQUIRED BY INVMTX IS PROPORTIONAL TO 
C                        THE QUANTITY N**3 .   ON THE NCAR
C                        CONTROL DATA 7600 WITH N = 100, INVMTX TAKES 
C                        ABOUT 1.6 TO 1.7 SECONDS DEPENDING ON THE TIME 
C                        REQUIRED TO FIND THE PIVOT ELEMENTS. 
C 
C PORTABILITY            THE NCAR ROUTINE, ULIBER, IS CALLED.  THERE IS 
C                        A FORTRAN VERSION OF ULIBER ON ULIB.  THERE IS 
C                        A MACHINE-DEPENDENT CONSTANT IN A DATA 
C                        STATEMENT IN INVMTX.  THIS CONSTANT SHOULD BE
C                        SET TO THE APPROPRIATE VALUE FOR THE MACHINE 
C                        BEING USED, OR THE CODE WHICH USES THE CONSTANT
C                        CAN BE DELETED AS INDICATED IN THE COMMENT 
C                        CARDS. 
C 
C REQUIRED RESIDENT      ULIBER 
C ROUTINES
C 
C NOTE                   IN ALMOST ALL APPLICATIONS (E.G., SOLVING
C                        LINEAR SYSTEMS OR DETERMINING PARTICULAR 
C                        ELEMENTS OF THE INVERSE), THE NSSL PACKAGE,
C                        LINEQSV, CAN PROVIDE THE SAME CAPABILITIES AS
C                        INVMTX APPROXIMATELY THREE TIMES AS FAST.
C                        HOWEVER, INVMTX PERFORMS FULL PIVOTING AND 
C                        REQUIRES LESS STORAGE SINCE IT CAN INVERT THE
C                        MATRIX IN PLACE. 
C 
C 
C 
C 
C 
      INTEGER         NA         ,NV         ,N          ,         IER
      DIMENSION     A(NA,N)    ,V(NV,N), IP(*)
      DOUBLE PRECISION A,V,D,HOLD,PVT,VH,VMAX
C 
C FLTMAX IS SET TO THE LARGEST FLOATING POINT NUMBER REPRESENTABLE IN 
C THE MACHINE 
C 
      CHARACTER *34 ERMESS
      DATA IEXMAX/38/ 
C 
C CALL INVMTA TO CHECK VALIDITY OF ARGUMENT RANGES
C 
      IER = IERINV(N,NA,NV) 
      IF (IER .NE. 0) RETURN
C 
C STORE A IN V
C 
      DO 102 J=1,N
         IP(J) = 0
         DO 101 I=1,N 
            V(I,J) = A(I,J) 
  101    CONTINUE 
  102 CONTINUE
      D = 1.
      IEX = 0 
      DO 110 M=1,N
         VMAX = 0.
         DO 104 J=1,N 
            IF (IP(J) .NE. 0) GO TO 104 
C 
C FIND MAXIMUM PIVOT ELEMENT
C 
            DO 103 I=1,N
               IF (IP(I) .NE. 0) GO TO 103
               VH = ABS(V(I,J)) 
               IF (VMAX .GE. VH) GO TO 103
               VMAX = VH
               K = I
               L = J
  103       CONTINUE
  104    CONTINUE 
         IP(L) = K
         NPM = N+M
         IP(NPM) = L
         D = D*V(K,L) 
  105    IF (ABS(D) .LE. 1.0) GO TO 106 
         D = D*0.1
         IEX = IEX+1
         GO TO 105
  106    CONTINUE 
         PVT = V(K,L) 
C 
C CHECK FOR NUMERCALLY SINGULAR MATRIX
C 
         IF (M .EQ. 1) PVTMX = ABS(PVT) 
         IF (ABS(PVT/FLOAT(M))+PVTMX .EQ. PVTMX) GO TO 113
C 
C INTERCHANGE ROWS, PLACING PIVOT ELEMENT ON THE DIAGONAL.
C THEN PROCEDE WITH PIVOTING. 
C 
         V(K,L) = 1.
         DO 107 J=1,N 
            HOLD = V(K,J) 
            V(K,J) = V(L,J) 
            V(L,J) = HOLD/PVT 
  107    CONTINUE 
         DO 109 I=1,N 
            IF (I .EQ. L) GO TO 109 
            HOLD = V(I,L) 
            IF (ABS(HOLD).LT.1.E-18) HOLD = 0.
            V(I,L) = 0. 
            DO 108 J=1,N
               V(I,J) = V(I,J)-V(L,J)*HOLD
  108       CONTINUE
  109    CONTINUE 
  110 CONTINUE
C 
C PERMUTE FINAL INVERSE MATRIX
C 
      M = N+N+1 
      DO 112 J=1,N
         M = M-1
         L = IP(M)
         K = IP(L)
         IF (K .EQ. L) GO TO 112
         D = -D 
         DO 111 I=1,N 
            HOLD = V(I,L) 
            V(I,L) = V(I,K) 
            V(I,K) = HOLD 
  111    CONTINUE 
  112 CONTINUE
C 
C CHECK FOR OVERFLOW IN DETERMINANT CALCULATION.
C 
      IF (IEX .GT. IEXMAX) GO TO 114
      D = D*10.**IEX
      RETURN
  113 IER = 33
      ERMESS = 'MATRIX SINGULAR IN INVMTX'
      GO TO 115 
  114 IER = 1 
      D = FLOAT(IEX)
      ERMESS = 'DETERMINANT TOO LARGE IN INVMTX'
      RETURN
  115 WRITE(1,116) IER,ERMESS 
  116 FORMAT('0  ....  ERROR NO. ',I4,' IN INVMTX',4X,A34)
      RETURN
      END 
      FUNCTION IERINV (N,NA,NV) 
C 
C THIS FUNCTION TESTS THE INPUT ARGUMENTS N, NA, AND NV TO INSURE 
C THAT THEY ARE IN THE PROPER RANGES.  THE INPUT ARGUMENTS HAVE THE 
C SAME MEANING AS IN INVMTX.
C 
C IF ONE OR MORE OF THE FOLLOWING CONDITIONS OCCURS, A FATAL ERROR
C MESSAGE IS WRITTEN TO UNIT 6 AND THE FUNCTION RETURNS WITH IERINV 
C EQUAL TO THE VALUE INDICATED. 
C 
C        FOR N .LT. 1, IERINV = 34. 
C        FOR NA .LT. N, IERINV = 35.
C        FOR NV .LT. N, IERINV = 36.
C 
C IF NO ERROR CONDITION IS DETECTED, THE FUNCTION RETURNS WITH
C IERINV = 0. 
C 
      CHARACTER*24 ERMESS 
C 
      IERINV = 0
      IF (N .GE. 1) GO TO 101 
      IERINV = 34 
      ERMESS = ' N .LT. 1 IN INVMTX ' 
      GO TO 103 
  101 IF (NA .GE. N) GO TO 102
      IERINV = 35 
      ERMESS ='NA .LT. N IN INVMTX '
      GO TO 103 
  102 IF (NV .GE. N) RETURN 
      IERINV = 36 
      ERMESS = ' NV .LT. N IN INVMTX '
  103 WRITE(1,104) IERINV,ERMESS
  104 FORMAT('0 ....  ERROR IN IERINV NO. ',I3,4X,A24)
      RETURN
      END 
C
C
C
C     --------------------------------------------------------
C     Section 2: Meteorological or mathematical routines.
C     --------------------------------------------------------
C
C
C     ------------------------------------------- 
      SUBROUTINE POLATE(N,X,Y,XIN,YOUT,M,BAD) 
C     ------------------------------------------- 
C 
      DIMENSION X(N),Y(N) 
      LOGICAL INC 
C 
      IF (N.LE.0) GOTO 500
C 
C     Determine if X increases or decreases 
C     ------------------------------------- 
C 
      XMAX = -9999999.
      XMIN = 9999999.
C
      DO 100 L = 1,N
         IF (X(L).NE.BAD .AND. X(L).GT.XMAX) THEN
            XMAX = X(L)
            LMAX = L
            ENDIF
         IF (X(L).NE.BAD .AND. X(L).LT.XMIN) THEN
            XMIN = X(L)
            LMIN = L
            ENDIF
100      CONTINUE
C
      IF (XMIN.GT.XMAX) GOTO 500
      IF (XMIN.EQ.XMAX) GOTO 500
      IF (XIN.LT.XMIN.OR.XIN.GT.XMAX) GOTO 500
C
      INC = .FALSE. 
      IF (LMAX.GT.LMIN) INC = .TRUE.
C 
C 
C     Check for exact match 
C     --------------------- 
C 
      DO 15 L=1,N 
         IF (XIN.EQ.X(L)) THEN
            YOUT=Y(L) 
            M = L 
            RETURN
            ENDIF 
15       CONTINUE
C
      IF (N.EQ.1) GOTO 500
C 
C     Interpolate
C     -----------
      DO 10 L=1,N-1 
         IF (XIN.GT.X(L) .AND. XIN.LT.X(L+1)) GOTO 50
         IF (XIN.LT.X(L) .AND. XIN.GT.X(L+1)) GOTO 50
10       CONTINUE          
C 
50    M=L+1 
      DUM=((X(M)-XIN)*(Y(M)-Y(L)))/(X(M)-X(L))
      YOUT=Y(M)-DUM 
C 
      IF(Y(M).EQ.BAD .OR. Y(L).EQ.BAD) YOUT = BAD
      IF (.NOT. INC) M = L
      RETURN
C 
C 
500   YOUT = BAD
      M = -1
      RETURN
      END 
C 
C 
C 
C     ----------------------------------------------- 
      FUNCTION POLATE2(X1,X2,Y1,Y2,XX,BAD)
C 
C     Interpolates between 2 points.
C     ----------------------------------------------- 
C 
      POLATE2 = BAD+1.
      IF (X1.EQ.BAD .OR. X2.EQ.BAD .OR. Y1.EQ.BAD .OR.
     *        Y2.EQ.BAD) POLATE2 = BAD
      IF (XX.GT.X1 .AND. XX.GT.X2) POLATE2 = BAD
      IF (XX.LT.X1 .AND. XX.LT.X2) POLATE2 = BAD
      IF (POLATE2.EQ.BAD) RETURN
C 
      IF (XX.EQ.X1) THEN
              POLATE2 = Y1
              RETURN
              ENDIF 
      IF (XX.EQ.X2) THEN
              POLATE2 = Y2
              RETURN
              ENDIF 
C 
      POLATE2 = Y1 + ((XX-X1)/(X2-X1))*(Y2-Y1)
      RETURN
      END 
C 
C 
C 
C     ------------------------------------
      FUNCTION UCMP(WDD,WSS) 
C 
C     Returns u-wind from direction, speed
C     ------------------------------------
C 
      WD = WDD
      WS = WSS
      IF (WD.LT.0. .OR. WS.EQ.-999.) THEN 
         UCMP = -999.
         ELSE 
         WD = 270. - WD 
         IF (WD .LT. 0.) WD = WD+360. 
         UCMP = WS*COS(WD/57.2958) 
         ENDIF
      RETURN
      END 
C 
C 
C 
C     ------------------------------------
      FUNCTION VCMP(WDD,WSS) 
C 
C     Returns v-wind from direction, speed
C     ------------------------------------
C 
      WD = WDD
      WS = WSS
      IF (WD.LT.0. .OR. WS.EQ.-999.) THEN 
         VCMP = -999.
         ELSE 
         WD = 270. - WD 
         IF (WD .LT. 0.) WD = WD+360. 
         VCMP = WS*SIN(WD/57.2958) 
         ENDIF
      RETURN
      END 
C 
C 
C 
C     ------------------------------------
      FUNCTION WDCOMP(UU,VV)
C 
C     Returns wind direction from u,v 
C     ------------------------------------
C 
      U = UU
      V = VV
      IF (U.EQ.-999. .OR. V.EQ.-999.) THEN
         WDCOMP = -999. 
         ELSE 
         IF (V.NE.1.E-3) V = V-1.E-3
         WDCOMP = 57.2958*ATAN2(-U,-V)  
         IF (WDCOMP.LT.0.) WDCOMP = WDCOMP+360. 
         ENDIF
      RETURN
      END 
C 
C 
C 
C     ------------------------------------
      FUNCTION WSCOMP(U,V)  
C 
C     Returns wind speed from u,v 
C     ------------------------------------
C 
      IF (U.EQ.-999. .OR. V.EQ.-999.) THEN
         WSCOMP = -999. 
         ELSE 
         WSCOMP = SQRT(U*U + V*V) 
         ENDIF
      RETURN
      END 
C 
C 
C 
C     ----------------------------------
      FUNCTION DEWPT(TE,RH) 
C     ----------------------------------
C 
      IF (TE.LT.-998. .OR. RH.GT.100. .OR. RH.LE.0.) THEN 
              DEWPT=-999. 
              RETURN
              ENDIF 
      TM=0.   
      TK=TE+273.16
      TN=TK 
      ES=10.**(9.4051-2353./TK) 
      EO=0.01*RH*ES 
442   EP=5417.98*ES/(TN*TN) 
      TN=TN+(EO-ES)/EP
      TS=ABS(TM-TN) 
      TM=TN 
      IF (TS.LT.0.005) GO TO 443
      ES=10.**(9.4051-2353./TN) 
      GO TO 442 
443   DD=TK-TN
      DEWPT=TK-273.16-DD
      IF (DEWPT.GT.TE) DEWPT=TE 
      END 
C 
C 
C 
C     ------------------------------------
      FUNCTION RELHU(T,TD,P)
C 
C     Returns RH (%) from T,TD (C),P (MB) 
C     ------------------------------------
C 
      IF (T.LE.-999. .OR. TD.LE.-999. .OR. P.LE.-999.) GOTO 100 
      IF (TD.GT.T) GOTO 100 
      RELHU = 100.0 * FMXR(TD,P)/FMXR(T,P)
      RETURN
C 
100   RELHU = -999. 
      RETURN
      END 
C 
C 
C 
C     ------------------------------------
      FUNCTION FMXR(T,P) 
C 
C     Returns mixing ratio (g/kg).  
C     If T=temp FMXR is saturated m.r. 
C     If T=Td then FMXR is ambient m.r.
C     ------------------------------------
C 
      TK = T + 273.16 
      E = 10.**(22.5518-(2937.4/TK)-4.9283*ALOG10(TK))*10. 
      FMXR = 621.98 * E/(P-E)
      RETURN
      END 
C 
C 
C 
C     ------------------------------------- 
      FUNCTION WLAPR(T,P)
C     ------------------------------------- 
C 
C     DOUBLE PRECISION F0,F1,F2 
C 
      H = (2500.-2.274*T)*1000. 
      TK = T+273.16 
      ES = 10.**(22.5518-(2937.4/TK)-4.9283*ALOG10(TK))*10.
      WS = 0.62198*ES/(P-ES)
      F0 = H/TK 
      F1 = 286.998+WS*F0
      F2 = 1004.*286.998+0.62198*WS*F0*F0 
      WLAPR = F1*286.998*TK/(P*F2) 
      RETURN
      END 
C 
C 
C 
C     --------------------------------------------------- 
      FUNCTION HYDROZ(P1,T1,H1,Z1,P2,T2,H2,IOPT,BAD)
C 
C     Geopotential height algorithm.  Returns height of 
C     level 2 given PTH of both levels and height of
C     level 1.
C 
C     P1      Pressure of level 1 (mb)
C     T1      Temperature of level 1 (c)
C     H1      Humidity of level 1 (DP if IOPT=1, RH if OPT=2) 
C     Z1      Geopotential height of level 1
C     P2      Pressure of level 2 (mb)
C     T2      Temperature of level 2 (c)
C     H2      Humidity of level 2 (DP if IOPT=1, RH if OPT=2) 
C     --------------------------------------------------- 
C 
C     Check for missing data
C     ----------------------
C 
      X = BAD
      IF (P1.LE.X .OR. T1.LE.X .OR. H1.LE.X .OR. Z1.LE.X
     *   .OR. P2.LE.X .OR. T2.LE.X .OR. H2.LE.X) THEN 
              HYDROZ = X
              RETURN
              ENDIF 
C 
C 
      IF (IOPT.EQ.2) THEN       !INPUT IS RH
              W1=QSATW(T1+273.16,P1)*H1*.01 
              W2=QSATW(T2+273.16,P2)*H2*.01 
              ELSE              !INPUT IS DEW PT
              W1=QSATW(H1+273.16,P1)
              W2=QSATW(H2+273.16,P2)
              ENDIF 
C 
C 
C     Virtual Temperature 
C     ------------------- 
C 
      TK1=(T1+273.16)*(1.0+0.609*W1)
      TK2=(T2+273.16)*(1.0+0.609*W2)
      TB = (TK1*ALOG(ABS(P1))+TK2*ALOG(ABS(P2))) / !FROM 19.43 NOTES
     *        (ALOG(ABS(P1))+ALOG(ABS(P2))) 
C 
C 
C     Compute height of level 2 
C     ------------------------- 
C 
      HYDROZ = Z1-(TB/0.03414)*ALOG(ABS(P2)/P1) 
      RETURN
      END 
C 
C 
C 
C     --------------------------------------------------- 
      FUNCTION HYDROP(P1,T1,H1,Z1,T2,H2,Z2,IOPT,BAD)
C 
C     Geopotential height algorithm.  Returns pressure of 
C     level 2 given THZ of both levels and pressure of
C     level 1.
C 
C     P1      Pressure of level 1 (mb)
C     T1      Temperature of level 1 (c)
C     H1      Humidity of level 1 (DP if IOPT=1, RH if OPT=2) 
C     Z1      Geopotential height of level 1
C     T2      Temperature of level 2 (c)
C     H2      Humidity of level 2 (DP if IOPT=1, RH if OPT=2) 
C     Z2      Geopotential height of level 2
C     --------------------------------------------------- 
C 
C     Check for missing data
C     ----------------------
C 
      X = BAD
      IF (P1.LE.X .OR. T1.LE.X .OR. H1.LE.X .OR. Z1.LE.X
     *   .OR. T2.LE.X .OR. H2.LE.X .OR. Z2.LE.X) THEN 
              HYDROP = X
              RETURN
              ENDIF 
C 
C 
C     Need first crude estimate of P2 to calculate  
C     mixing ratio and mean virtual temperature 
C     --------------------------------------------
C 
      TK1 = T1+273.16 
      TK2 = T2+273.16 
      TB = 2.0*TK1*TK2/(TK1+TK2)
      P2 = P1*EXP(0.03414*(Z1-Z2)/TB) 
C 
      IF (IOPT.EQ.2) THEN       !INPUT IS RH
              W1=QSATW(T1+273.16,P1)*H1*.01 
              W2=QSATW(T2+273.16,P2)*H2*.01 
              ELSE              !INPUT IS DEW PT
              W1=QSATW(H1+273.16,P1)
              W2=QSATW(H2+273.16,P2)
              ENDIF 
C 
C 
C     Virtual Temperature 
C     ------------------- 
C 
      TK1=(T1+273.16)*(1.0+0.609*W1)
      TK2=(T2+273.16)*(1.0+0.609*W2)
      TB = (TK1*ALOG(ABS(P1))+TK2*ALOG(ABS(P2))) / !FROM 19.43 NOTES
     *        (ALOG(ABS(P1))+ALOG(ABS(P2))) 
C 
C 
C     Compute P2
C     ----------
C 
      HYDROP = P1*EXP(0.03414*(Z1-Z2)/TB) 
      RETURN
      END 
C 
C 
C 
C     ----------------------------------
      FUNCTION QSATW(TA,P)
C     ----------------------------------
C 
      DATA PS/1013.246/,TS/373.16/
      E1=11.344*(1.0-TA/TS) 
      E2=-3.49149*(TS/TA-1.0) 
      F1=-7.90298*(TS/TA-1.0) 
      F2=5.02808*ALOG10(TS/TA)
      F3=-1.3816*(10.0**E1-1.0)*1.E-7 
      F4=8.1328*(10.0**E2-1.0)*1.E-3
      F5=ALOG10(PS) 
      F=F1+F2+F3+F4+F5
      ES=10.0**F
      QSATW=.62197*ES/(P-ES)
      RETURN
      END 
C 
C 
C 
C     ----------------------------------
      FUNCTION QSATI(TA,P)
C     ----------------------------------
C 
      DATA PO/6.1071/,TO/273.16/
      F1=-9.09718*(TO/TA-1.0) 
      F2=-3.56654*ALOG10(TO/TA) 
      F3=0.876793*(1.0-TA/TO) 
      F4=ALOG10(PO) 
      ES=10.0**(F1+F2+F3+F4)
      QSATI=0.62197*ES/(P-ES) 
      RETURN
      END 
C 
C 
C 
C     ----------------------------------------------------------
      FUNCTION TROPLAPSP(P,BAD)
C 
C     This routines returns the mean tropical (West Indies) 
C     lapse rate for July-October according to Jordan, 1958 
C     (Journal of Meteorology, Vol 15, p 94). 
C 
C     P= PRESSURE (MB)
C     TROPLAPSP = LAPSE RATE AT P  (dT/dP  deg/mb)
C     ----------------------------------------------------------
C 
C 
      PARAMETER (NT=19) 
      DIMENSION PR(NT), RT(NT)
      DATA PR/100.,150.,200.,250.,300.,350.,400.,450.,500.,550.,
     *        600.,650.,700.,750.,800.,850.,900.,950.,1000./
      DATA RT/-0.053,.222,.243,.220,.185,.155,.129,.108,.094, 
     *        .083,.076,.072,.067,.060,.055,.052,.057,.062,.051/
C 
C 
      CALL POLATE(NT,PR,RT,P,RATE,IFLAG,BAD)
C
      IF (P.EQ.BAD) THEN
         TROPLAPSP = BAD
         RETURN
         ENDIF
C
      IF(IFLAG.EQ.-1) THEN
        RATE=RT(1)
        IF(P.GT.PR(NT)) RATE=RT(NT) 
        ENDIF 
C
      TROPLAPSP = RATE 
C
      RETURN
      END 
C
C
C
C     ----------------------------------------------------------
      FUNCTION TROPLAPSZ(P,BAD)
C 
C     This routines returns the mean tropical (West Indies) 
C     lapse rate for July-October according to Jordan, 1958 
C     (Journal of Meteorology, Vol 15, p 94). 
C 
C     P= PRESSURE (MB)
C     TROPLAPSZ = LAPSE RATE AT P  (dT/dZ  deg/m)
C     ----------------------------------------------------------
C 
C 
      PARAMETER (NT=19) 
      DIMENSION PR(NT), RT(NT)
      DATA PR/100.,150.,200.,250.,300.,350.,400.,450.,500.,550.,
     *        600.,650.,700.,750.,800.,850.,900.,950.,1000./
      DATA RT/+.00091,-.00529,-.00750,-.00811,-.00786,-.00743,-.00687,
     *        -.00633,-.00601,-.00574,-.00565,-.00571,-.00566,-.00536,
     *        -.00518,-.00515,-.00591,-.00672,-.00566/
C 
C 
      CALL POLATE(NT,PR,RT,P,RATE,IFLAG,BAD)
C
      IF (P.EQ.BAD) THEN
         TROPLAPSZ = BAD
         RETURN
         ENDIF
C
      IF(IFLAG.EQ.-1) THEN
        RATE=RT(1)
        IF(P.GT.PR(NT)) RATE=RT(NT) 
        ENDIF 
C
      TROPLAPSZ = RATE 
C
      RETURN
      END 
C
C
C
C     ----------------------
      FUNCTION HMSTS(X) 
C     ----------------------
C 
      KHR=(X+0.1)*0.0001
      KMIN=(X-KHR*10000 + 1)*.01 
      SEC=(X-KHR*10000-KMIN*100)
      HMSTS=KHR*3600.+KMIN*60.+SEC 
C	write (*,'(2x,f10.1,2x,I3,2x,f10.1, 2x, i3,2x,f4.0,2x,f10.1)')
C     *  x, khr, (x-khr*10000), kmin, sec, hmsts
      RETURN
      END 
C 
C 
C 
C     ----------------- 
      FUNCTION STHMS(X) 
C     ----------------- 
C 
      IHR = INT(X/3600.)
      IMN = INT((X-(IHR*3600.))/60.)
      SC = X-(IHR*3600.)-(IMN*60.)
      STHMS = IHR*10000.+IMN*100.+SC 
      RETURN
      END 
C
C
C
C
C     ----------------------------------------------------------------- 
      SUBROUTINE STMREL(DAYMN,DAYOB,GMTOB,XD,YD,NFD,DAYSTM,GMTSTM,
     *                  XSTM,YSTM,DAY0,GMT0)
C     ----------------------------------------------------------------- 
C 
C 
      DIMENSION DAYSTM(NFD),GMTSTM(NFD),XSTM(NFD),YSTM(NFD) 
      DIMENSION DTFIX(50) 
      LOGICAL MATCH 
C 
      MATCH=.FALSE. 
      NFIX = NFD
C 
C 
C     Find storm position at reference time 
C     ------------------------------------- 
C 
      DO 110 L=1,NFIX 
              IF (DAY0.EQ.DAYSTM(L).AND.GMT0.EQ.GMTSTM(L)) THEN 
                 XSTM0=XSTM(L)
                 YSTM0=YSTM(L)
                 MATCH=.TRUE. 
                 GOTO 120 
                 ENDIF
110           CONTINUE
120   IF (.NOT.MATCH) THEN
              WRITE(1,'("*** BAD REFERENCE TIME IN STMREL ***")') 
              STOP
              ENDIF 
C 
C 
200   CALL XTRTIME(DAYMN,000.0,IYMN,IMMN,IDMN,HRMN) 
      DO 210 N=1,NFIX 
        CALL XTRTIME(DAYSTM(N),GMTSTM(N),IYR,IMO,IDA,HR)
        CALL DIFTIME(IYMN,IMMN,IDMN,HRMN,IYR,IMO,IDA,HR,DTFIX(N))   
210     CONTINUE
C 
C 
      CALL XTRTIME(DAYOB,GMTOB,IYR,IMO,IDA,HR)
      CALL DIFTIME(IYMN,IMMN,IDMN,HRMN,IYR,IMO,IDA,HR,DTOB) 
      NFIX1=NFIX-1
C 
      DO 300 N=1,NFIX1
        NF=N
        IF (DTOB.GE.DTFIX(N) .AND. DTOB.LE.DTFIX(N+1)) GOTO 310 
300     CONTINUE
C 
      WRITE(1,301) DAYOB,GMTOB
301   FORMAT(' **** CANNOT INTERPOLATE FIXES FOR ',I6,1X,I4,' ****')
      STOP
C 
C 
310   DT=DTFIX(NF+1)-DTFIX(NF)
      IF (ABS(DT) .LT. 1.E-4) GOTO 910
      RAT = (DTOB-DTFIX(NF))/(DTFIX(NF+1)-DTFIX(NF))
      XSTOB = XSTM(NF)+(XSTM(NF+1)-XSTM(NF))*RAT
      YSTOB = YSTM(NF)+(YSTM(NF+1)-YSTM(NF))*RAT
      XD = XD - XSTOB + XSTM0 
      YD = YD - YSTOB + YSTM0       
      RETURN
C 
C 
C     Errors
C     ------
C 
910   WRITE(1,911) NF,DTFIX(NF),DTFIX(NF+1) 
911   FORMAT("**** FIX TIMES ARE EQUAL IN STMREL***",I5,2F15.5) 
      STOP
      END 
C 
C 
C 
C     --------------------------------------------
      SUBROUTINE XTRTIME(DATE,GMT,IYR,IMO,IDA,HR) 
C     --------------------------------------------
C 
      IYR=INT(DATE/10000.)
      IMO=INT((DATE-FLOAT(IYR)*10000.)/100.)
      IDA=INT(DATE-FLOAT(IYR)*10000.-IMO*100.)
      HR=FLOAT(NINT(GMT/100.-0.3))
      HR=HR+(GMT-HR*100.)/60. 
      RETURN
      END 
C 
C 
C 
C     --------------------------------------------------------
      SUBROUTINE DIFTIME(IYZ,IMZ,IDZ,HRZ,IYX,IMX,IDX,HRX,DHR) 
C     --------------------------------------------------------
C 
      DIMENSION MODA(12)
      DATA MODA/31,28,31,30,31,30,31,31,30,31,30,31/
C 
      DHR=0.0 
      IDZZ=IDZ
      IMZZ=IMZ
      HRZZ=HRZ
C 
      IF(IYZ.EQ.IYX) GOTO 100 
      DHR=24.-HRZZ+FLOAT((MODA(IMZZ)-IDZZ)*24+MAX0(IYX-(IYZ+1),0)*8760) 
      IMZZ=MOD(IMZZ+1,12) 
      IDZZ=1
      HRZZ=0.0
C 
      IF (IMZZ.EQ.1) GOTO 100 
      DO 10 I=IMZZ,12 
10    DHR=DHR+FLOAT(MODA(I)*24) 
      IMZZ=1
C 
100   IF (IMZZ.EQ.IMX) GOTO 200 
      DHR=DHR+FLOAT((MODA(IMZZ)-IDZZ)*24)+24.-HRZZ
      IMZZ1=IMZZ+1
C 
      IF (IMZZ1 .GE. IMX) GOTO 120
      IMX1=IMX-1
      DO 110 I=IMZZ1,IMX1 
110   DHR=DHR+FLOAT(MODA(I)*24) 
C 
120   IMZZ=IMZZ1
      IDZZ=1
      HRZZ=0.0
C 
200   DHR=DHR+FLOAT((IDX-IDZZ-1)*24)+24.-HRZZ+HRX 
      RETURN
      END 
C 
C 
C 
C     --------------------------------------------------
      SUBROUTINE DAYCMPR(FLID,DATE,TODAY,TOMORROW)
C     --------------------------------------------------
C 
C 
      DIMENSION IMONTH(12)
      LOGICAL TODAY,TOMORROW
      DATA IMONTH/31,28,31,30,31,30,31,31,30,31,30,31/
C     
      TODAY = .FALSE. 
      TOMORROW = .FALSE.
C 
C     SEPERATE DATES INTO YEAR,MONTH,DAY
C     ----------------------------------
C 
      IYR = INT(FLID/10000.)
      IMN = INT(FLID/100.) - (IYR*100.) 
      IDY = FLID - (IYR*10000.) - (IMN*100.)
      JYR = INT(DATE/10000.)
      JMN = INT(DATE/100.) - (JYR*100.) 
      JDY = DATE - (JYR*10000.) - (JMN*100.)
      IF (MOD(IYR,4).EQ.0) IMONTH(2) = 29 
C 
C     FIND IF SAME DAY OR NEXT
C     ------------------------
C 
      IF (FLID.EQ.DATE) TODAY = .TRUE.
      IF (FLID+1.EQ.DATE) TOMORROW = .TRUE. 
C 
C     CHECK FOR MONTHLY OR YEARLY FLIP OVERS
C     --------------------------------------
C 
      IF (JDY.EQ.1.AND.IDY.EQ.IMONTH(IMN)) THEN 
        IF (JMN.EQ.IMN+1) TOMORROW = .TRUE. 
        IF (JMN.EQ.1.AND.IMN.EQ.12) THEN
          IF (JYR.EQ.IYR+1) TOMORROW = .TRUE. 
          ENDIF 
        ENDIF 
      RETURN
      END 
C
C
C
C     --------------------------------------------------------
      SUBROUTINE LOPASS(GIN,GOUT,NMAX,FRAC,NTRM)
C     Complete symmetric filter with multiplicative adjustment
C     --------------------------------------------------------
C 
C 
      PARAMETER(MTERMS=100) 
      DIMENSION GIN(NMAX),GOUT(NMAX),WT(MTERMS) 
      DATA PI/3.1415926535/ 
      IPRT=1
      NTRM1=NTRM+1
      NMAX1=NMAX-NTRM 
      OMCUT=PI*FRAC 
      WTO=FRAC
C 
C     Endpoints 
C     --------- 
C 
      GOUT(1)=GIN(1)
      GOUT(NMAX)=GIN(NMAX)
C 
C 
      DO 40 L=2,NMAX-1
              SUM=0.0 
              MTRM=NTRM 
              IF(L.LE.NTRM) MTRM=L-1
              IF(L.GE.NMAX+1-MTRM) MTRM=NMAX-L    !FIX 2/17/88
C 
C     Calculate Weights 
C     ----------------- 
C 
      WTO = FRAC
      WFSUM = 0.0 
      DO 20 N=1,MTRM
              FAC=FLOAT(N)*2.*PI/(2.*MTRM+1)
              WT(N)=SIN(FLOAT(N)*OMCUT)/(PI*FLOAT(N))*SIN(FAC)/FAC
              WFSUM=WFSUM+2.0*WT(N) 
20            CONTINUE
C 
      WFSUM=WFSUM+WTO 
      WTO=WTO/WFSUM 
C 
      DO 25 N=1,MTRM
25            WT(N)=WT(N)/WFSUM 
              DO 30 N=1,MTRM
30                  SUM=SUM+(GIN(L-N)+GIN(L+N))*WT(N) 
              GOUT(L)=SUM+GIN(L)*WTO    
40            CONTINUE
      RETURN
      END 
C 
C 
C 
C     --------------------------------------------
      SUBROUTINE GAP(X,NT,IBGN,IEND,NGAPS,BAD)
C     --------------------------------------------
C 
      PARAMETER (MXGPS = 100) 
C 
      DIMENSION X(NT),IBGN(MXGPS),IEND(MXGPS) 
      LOGICAL LAST
C
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
C 
      IGAP=0
      LAST=.FALSE.
      NGAPS=0 
C 
      DO 5 L=1,MXGPS
         IBGN(L)=0
         IEND(L)=0
5        CONTINUE 
C 
      DO 10 L=1,NT
      IF (X(L).GT.BAD) THEN
              IF (LAST) GO TO 10
              IGAP=IGAP+1 
              IBGN(IGAP)=L
              LAST=.TRUE. 
              GO TO 10
      ELSE
              IF(.NOT.LAST) GO TO 10
              IEND(IGAP)=L-1
              LAST=.FALSE.
              GO TO 10
      ENDIF 
10    CONTINUE
C 
      IF(X(NT).GT.BAD) IEND(IGAP)=NT 
      NGAPS=IGAP-1
      IF (NGAPS.LT.0) NGAPS=0 
C
      IF (NGAPS.GT.MXGPS) THEN
         WRITE(LUT,'("* FATAL ERROR: TOO MANY GAPS IN SOUNDING *")')
         STOP
         ENDIF
C
      RETURN
      END 
C
C
C
C     --------------------------------------------------------
C     Section 3: Routines relating to message decoding.
C     --------------------------------------------------------
C
C
c     ----------------------------------------------
      subroutine drop(iwx,iflag,iyrs,imns,idys,line)
c
c     decodes tempdrop message into spline format
c     ----------------------------------------------
c
      character*2 header(12)
      character*30 blank
      character*70 line
      dimension prs(12)
      logical plev(12),knots,skpwind
c
      data header /'99','00','92','85','70','50','40','30',
     *             '25','20','15','10'/
      data prs /1070.,1000.,925.,850.,700.,500.,400.,300.,
     *          250.,200.,150.,100./
      data plev /12*.false./
      data blank /'                             '/
c
      lvl=0
      sfcp = 9999.
c
c     ------------------------------------------------
c     if iflag=1 then we are already at the XXAA line
c     ------------------------------------------------
c
c     read the line with the mission number
c     -------------------------------------
 10   if (iflag.ne.1) then
         read(12,'(a)',end=99)
         if(line(1:30).eq.blank)goto 10
         endif
c
c     read the first data line
c     ------------------------
 40   if (iflag.ne.1) then
         read(12,'(a)',end=99)line
         if(line(1:30).eq.blank)goto 40
         endif
c
c
c     read the value of the day
c     check if the winds are in knots or m/s
c     --------------------------------------
      read (line(7:8),'(i2)') iday
      if (iday.gt.50) then
         iday=iday-50
         knots = .true.
         else
         knots = .false.
         endif
c
c     check for month, year flips
c     ---------------------------
      yy = iyrs
      mm = imns
      if (idys.gt.27.and.iday.eq.1) then
         mm = imns+1
         if (mm.eq.13) then
            mm = 1
            yy = iyrs + 1.
            if (yy.eq.100.) yy = 00.
            endif
         endif
c
      yymmdd = yy * 10000. + float(mm) * 100. + float(iday)
c
c     read the value of the hour
c     --------------------------
      read (line(9:10),'(i2)') ihour
      igmt = ihour * 100.
c
c     set the value of the highest mandatory wind level reporting
c     -----------------------------------------------------------
      if(line(11:11).eq.'/')then
        maxlvl=1
      elseif(line(11:11).eq.'0')then
        maxlvl=2
      elseif(line(11:11).eq.'9')then
        maxlvl=3
      elseif(line(11:11).eq.'8')then
        maxlvl=4
      elseif(line(11:11).eq.'7')then
        maxlvl=5
      elseif(line(11:11).eq.'5')then
        maxlvl=6
      elseif(line(11:11).eq.'4')then
        maxlvl=7
      elseif(line(11:11).eq.'3')then
        maxlvl=8
      elseif(line(11:11).eq.'2')then
        maxlvl=10
      elseif(line(11:11).eq.'1')then
        maxlvl=12
      endif
c
c     read the latitude
c     -----------------
      read (line(15:17),'(f3.1)') alat
c
c     read the longitude
c     ------------------
      read (line(20:23),'(f4.1)') alon
c
c     go to column 31 to read the surface group
c     -----------------------------------------
      itag=31
c
c
c     Go on to next level
c     --------------------------------------------
200   do 205 l = 1,12
         plev(l)=.false.
205      continue
      press = -99.
      temp = -99.
      rh = -99.
      geopot = -99.
      wdir = -99.
      wspd = -99.
      skpwind = .false.
c
c     count the number of the mandatory level
c
      lvl=lvl+1
c
c     check to see if 925 level is missing
c
      if (lvl.eq.3 .and. line(itag:itag+1).eq.'85') lvl=lvl+1
c
c     check the data for each mandatory level
c
      if(line(itag:itag+1).eq.header(lvl))then
        plev(lvl)=.true.
        press=prs(lvl)
        call geo (line,itag,plev,geopot,sfcp)
        itag=itag+6
        call tagtst(itag,line)
        pressx = press
        if (press.eq.1070. .and. sfcp.le.1070.) pressx = sfcp
        call temdew (line,itag,pressx,temp,rh)
        if(lvl.le.maxlvl)then
          itag=itag+6
          call tagtst(itag,line)

c         check if sfc wind group is missing
c         ----------------------------------
          if (lvl.eq.1 .and. line(itag:itag+1) .eq. '00') then
             skpwind = .true.
             else
             call wind (line,itag,wdir,wspd,*99)
             call dstouv (wdir,wspd,alat,alon,knots)
             skpwind = .false.
             endif
c
          endif
        if (temp .ne. -99. .or. rh .ne. -99. .or. geopot .ne. -99.
     1        .or. wdir .ne. -99. .or. wspd .ne. -99.) call out
     2        (iwx,yymmdd,igmt,alat,alon,press,temp,rh,geopot,wdir,
     3         wspd,1)
        if (.not. skpwind) then
           itag=itag+6
           call tagtst(itag,line)
           endif
        go to 200
      elseif (line(itag:itag+1) .eq. '88') then
        itag=itag+6
        call tagtst(itag,line)
        if(line(itag+2:itag+4).ne.'999')then
          read (line(itag+2:itag+4),'(f3.0)') press
          itag=itag+6
          call tagtst(itag,line)
          read (line(itag:itag+1),'(f2.0)') wdir
          read (line(itag+2:itag+4),'(f3.0)') wspd
          if(wspd.ge.500.0)then
            wspd=wspd-500
            wdir=wdir+0.5
          endif
          call dstouv (wdir,wspd,alat,alon,knots)
          if (temp .ne. -99. .or. rh .ne. -99. .or. geopot .ne. -99.
     1          .or. wdir .ne. -99. .or. wspd .ne. -99.) call out (iwx,
     2          yymmdd,igmt,alat,alon,press,temp,rh,geopot,wdir,wspd,6)
        endif
      endif
c
c     look for significant level data
c
 60   read(12,'(a)',end=99)line
c
c     check if the line has data
c
      if(line(1:30).eq.blank)goto 60
c
c     check significant level data
c
 75   if(line(1:4).eq.'XXBB')then
        itag=31
        ihead=-11
 70     ihead=ihead+11
        if(ihead.gt.99)ihead=11
        read(line(itag:itag+1),'(i2)')jhead
        if(jhead.eq.ihead)then
          if(line(itag+2:itag+4).eq.'///') then
             itag=itag+6
             call tagtst(itag,line)
             itag=itag+6          
             goto 71
             endif
          read(line(itag+2:itag+4),'(i3)')iprs
          if(iprs.lt.100)iprs=iprs+1000.
          press=iprs
          itag=itag+6
          call tagtst(itag,line)
          if(line(itag:itag+4).eq.'21212')goto 75
          call temdew (line,itag,press,temp,rh)
          if (temp .ne. -99. .or. rh .ne. -99. .or. geopot .ne. -99.
     1          .or. wdir .ne. -99. .or. wspd .ne. -99.) call out (iwx,
     2          yymmdd,igmt,alat,alon,press,temp,rh,geopot,wdir,wspd,2)
          itag=itag+6
 71       call tagtst(itag,line)
          if(line(itag:itag+4).eq.'21212')goto 75
          goto 70
        endif
c
      elseif(line(itag:itag+4).eq.'21212')then
        itag=19
        ihead=0
 30     ihead=ihead+11
        if(ihead.gt.99)ihead=11
        read(line(itag:itag+1),'(i2)')jhead
        if(jhead.eq.ihead)then
          read(line(itag+2:itag+4),'(i3)')iprs
          if(iprs.lt.100)iprs=iprs+1000.
          press=iprs
          itag=itag+6
          call tagtst(itag,line)
          read (line(itag:itag+1),'(f2.0)') wdir
          read (line(itag+2:itag+4),'(f3.0)') wspd
          if(wspd.ge.500.0)then
            wspd=wspd-500
            wdir=wdir+0.5
          endif
          call dstouv (wdir,wspd,alat,alon,knots)
          temp=-99.
          rh=-99.
          geopot=-99.
          if (temp .ne. -99. .or. rh .ne. -99. .or. geopot .ne. -99.
     1          .or. wdir .ne. -99. .or. wspd .ne. -99.) call out (iwx,
     2          yymmdd,igmt,alat,alon,press,temp,rh,geopot,wdir,wspd,2)
          itag=itag+6
          call tagtst(itag,line)
          goto 30
        endif
        goto 75
c
      elseif(line(itag:itag+4).eq.'31313')then
        itag = itag+19
        call tagtst(itag,line)
        goto 75
c
      elseif(line(itag:itag+4).eq.'51515')then
        itag = itag+6
        call tagtst(itag,line)
500     if (line(itag:itag+4).eq.'10190') then
          itag = itag+6
          call tagtst(itag,line)
          do 505 l = 1,12
            plev(l)=.false.
505         continue
          press = -99.
          temp = -99.
          rh = -99.
          geopot = -99.
          wdir = -99.
          wspd = -99.
c
          do 510 l = 1,12
            if(line(itag:itag+1).eq.header(l))then
               plev(l)=.true.
               press=prs(l)
               call geo (line,itag,plev,geopot,sfcp)
               if (geopot .ne. -99.) call out(iwx,yymmdd,igmt,
     *         alat,alon,press,temp,rh,geopot,wdir,wspd,3)
               endif
510         continue
          itag = itag+6
          call tagtst(itag,line)
          goto 500
          endif
c
      elseif(line(1:4).eq.'NNNN')then
        return
      endif
c
 99   return
      end
c
c
c
      subroutine uvcomp (dir,spd)
c 
c     this subroutine changes dir to u, and spd to v, where dir is
c     given in meteorological degrees.  The original values of dir
c     and spd are destroyed.
c
      degrad = atan(1.0) / 45.
      dirdg = 270.0 - dir
      if (dirdg .lt. 0.0) dirdg = dirdg + 360.
      dirrd = dirdg * degrad
      dir = spd * cos(dirrd)
      spd = spd * sin(dirrd)
      return
      end
c
c
c
      subroutine temdew (line,lptr,press,temp,rh)
      character*70 line
c
c     extract the temperature
c
      if (line(lptr:lptr+2) .ne. '///') then
        read (line(lptr:lptr+2),'(f3.1)') atemp
        read (line(lptr+2:lptr+2),'(i1)') ifrac
        if (mod(ifrac,2) .eq. 0) then
          temp = atemp
        else
          temp = -atemp
        endif
      endif
c
c     extract the dewpoint depression
c
      if (line(lptr+3:lptr+4) .ne. '//') then
        read (line(lptr+3:lptr+4),'(i2)') idd
        if (idd .gt. 50) then
          dd = float (idd - 50)
        else
          dd = float (idd) / 10.
        endif
        dewpt = temp - dd
        call relhum (press,temp,dewpt,rh)
      endif
      return
      end
c
c
c
      subroutine relhum (press,temp,dewpt,rh)
      parameter (tkelvn = 273.16)
      parameter (em = 9.4051)
      parameter (e = 2353.)
c
c     compute the relative humidity using the vapor pressure vprs
c     and the saturation vapor pressure svprs
c
      vprs = 10**(em - e / (dewpt + tkelvn))
      svprs = 10**(em - e / (temp + tkelvn))
      fmixr = vprs / (press - vprs)
      smixr = svprs / (press - svprs)
      rh = 100. * fmixr / smixr
      if(rh.gt.100.)rh=100.
      return
      end
c
c
c
      subroutine geo (line,lptr,plev,geopot,sfcp)
      character*70 line
      logical plev
      dimension plev(12)
c
c     extract the geopential height (modifications by JLF 11/92)
c
      if (line(lptr+2:lptr+4) .ne. '///') then
        read (line(lptr+2:lptr+4),'(f3.0)') geopot 
c
        if (plev(1)) then                          ! Surface
          if (geopot .lt. 100.) geopot = geopot + 1000.
          sfcp = geopot
          endif
c
        if (plev(2)) then                          ! 1000 mb
          if (geopot .ge. 500.) geopot = -(geopot-500.)
          endif
c
        if (plev(3)) then
          if (sfcp.le.925..and.geopot.ge.500.) geopot=-(geopot-500.)
          endif
c
        if (plev(4)) then
          geopot = geopot+1000.
          if (sfcp.le.950..and.geopot.gt.1500.) geopot=geopot-1000.
          endif
c
        if (plev(5)) then                          ! 700 mb
          add = 2000.
          if (geopot .lt. 500.) add = 3000.
          if (sfcp.lt.960.) add = 2000.    
          geopot = geopot + add
          endif
c
        if (plev(6) .or. plev(7)) then             ! 500, 400 mb
          geopot = geopot * 10.
          endif
c
        if (plev(8) .or. plev(9) .or. plev(10)     ! >= 300 mb
     *      .or. plev(11) .or. plev(12)) then
          geopot = geopot * 10.
          if (geopot.lt.8500.) geopot = geopot + 10000.
          endif
c
      endif
      return
      end
c
c
c
      subroutine wind (line,lptr,wdir,wspd,*)
      character*70 line
c
c     extract the wind direction and speed
c
      if (line(lptr:lptr+4) .ne. '/////') then
        read (line(lptr:lptr+1),'(f2.0)') wdir
        read (line(lptr+2:lptr+4),'(f3.0)') wspd
      else
        wdir = -99.
        wspd = -99.
      endif
      return
      end
c
c
c
      subroutine dstouv (wdir,wspd,alat,alon,knots)
      logical knots
      real alat,alon,wdir,wspd
c
c     convert wind direction and speed to u, v (in m/s)
c
      if (wdir .ne. -99.) then
        wdir = wdir * 10.
        if(wspd.ge.500.0)then
          wspd=wspd-500.
          wdir=wdir+5
        endif
        if (knots) wspd = 0.514444 * wspd
        call uvcomp (wdir,wspd)
      endif
      return
      end
c
c
c
      subroutine tagtst(itag,line)
      character*70 line
c
c     check if the end of the line has been reached and the next line should 
c     be read
c  
c      if(itag.gt.47)then
        if(itag.lt.66)then
          do i=itag,itag+5
            if(line(i:i).ne.' ')return
          end do
        endif
        read(12,'(70a)',end=99)line
        itag=1
c     endif
 99   return
      end
c
c
c
      subroutine recco(iyrs,imns,idys,line)
      character*1 quad
      character*30 blank
      character*70 line
      logical knots
      data blank /'                             '/
c
c
c     Read the day
c     ------------
      read(line(13:14),'(i2)') iday
      if(line(13:13).eq.' ') read (line(14:15),'(i2)') iday
      if (iday.gt.50) then
         iday=iday-50
         knots = .true.
         else
         knots = .false.
         endif
c
c     check for month, year flips
c     ---------------------------
      yy = iyrs
      mm = imns
      if (idys.gt.27.and.iday.eq.1) then
         mm = imns+1
         if (mm.eq.13) then
            mm = 1
            yy = iyrs + 1.
            if (yy.eq.100.) yy = 00.
            endif
         endif
c
      yymmdd = yy * 10000. + float(mm) * 100. + float(iday)
c
c
c     read the next line
c
 20   read(12,'(a)',end=99)line
c
c     check if the line has information
c
      if(line(1:30).eq.blank)goto 20
c
c     define the data type
c
      if(index(line,'AF').ne.0)iwx=6
      if(index(line,'NOAA').ne.0)iwx=3
c
c     read the data line
c
 10   read(12,'(a)',end=99)line
c
c     if line is NNNN return
c
      if(index(line,'NNNN').ne.0)return
c
c     check if the line has data
c
      if(line(1:30).eq.blank)goto 10
c
c     recco's begin with 97779 or 95559
c
      if (line(1:5) .eq. '97779'.or. line(1:5).eq.'95559') then
        read (line(7:8),'(i2)') ihour
        read (line(9:10),'(i2)') min
        quad = line(14:14)
        read (line(15:17),'(f3.1)') alat
        read (line(19:21),'(f3.1)') alon
        if (quad .eq. '1' .and. line(19:19) .ne. '9')alon = alon + 100.
        read (line(31:32),'(f2.0)') wdir
        read (line(33:35),'(f3.0)') wspd
        if(wspd.ge.500.0)then
          wspd=wspd-500.
          wdir=wdir+0.5
        endif
        call rtmdew (line,temp,dewpt,*99)
        call rpress (line,press,geopot)
        call relhum (press,temp,dewpt,rh)
        call rdstuv (wdir,wspd,alat,alon)
        igmt = float(ihour) * 100. + min
        call out (iwx,yymmdd,igmt,alat,alon,press,
     1            temp,rh,geopot,wdir,wspd,4)
        go to 10
      endif
 99   return
      end
c
c
c
      subroutine rtmdew (line,temp,dwpt,*)
      character*70 line,line2
c
c     extract the temperature from the RECCO
c
      read (line(37:38),'(f2.0)',err=99) atemp
      if (atemp .lt. 50.) then
        temp = atemp
      else
        atemp = atemp - 50.
        temp = -atemp
        endif
c
c     if the dewpoint is missing, it may be in plain text on line 2
c
      if (line(39:40) .eq. '//') then
        read (12,'(a)',end=99) line2
        if (line2(1:3) .eq. 'DEW') then
          read (line2(11:15),'(f5.1)') dewpt
        else
          dewpt = 0.0
          endif
        dwpt = dewpt
        go to 20
        endif
c
c     otherwise, get the dewpoint from the main line
      read (line(39:40),'(f2.0)',err=99) dewpt
      if (dewpt .lt. 50.) then
        dwpt = dewpt
      else
        dewpt = dewpt - 50.
        dwpt = -dewpt
        endif
20    continue
      return
99    return 1
      end
c
c
c
      subroutine rpress (line,press,geopot)
      character*70 line
      integer prsind
      dimension sprs(7)
      data sprs /200.,850.,700.,500.,400.,300.,250./
c
c     extract the pressure and geopotential from the RECCO message
c
      read (line(44:44),'(i1)') prsind
      if (prsind .eq. 0) then
        read (line(45:47),'(f3.0)') geopot
        if (geopot .lt. 800.) geopot = geopot + 1000.
        press = 1070.
      elseif (prsind .eq. 9) then
        geopot = -99.
        read (line(25:27),'(f3.0)') tralt
        pralt = tralt * 10.
        press = 1013.25 * (1. - (pralt / 44331.)) ** (1. / 0.190263)
      elseif (prsind .ge. 1 .and. prsind .le. 7) then
        press = sprs(prsind)
        read (line(45:47),'(f3.0)') geopot
        if (prsind .gt. 3 .and. prsind .lt. 7) then
          geopot = geopot * 10.
        elseif (prsind .eq. 2) then
          if (geopot .lt. 800.) geopot = geopot + 1000.
        elseif (prsind .eq. 1 .or. prsind .eq. 7) then
          geopot = geopot * 10.
          if (geopot .lt. 800.) geopot = geopot + 1000.
        elseif (prsind.eq.3)then
          geopot=geopot+3000.
        endif
      else if (prsind .eq. 8) then
        read (line(25:27),'(f3.0)') tralt   ! true alt in decameters
        read (line(45:47),'(f3.0)') dvalue  ! d-value in decameters
        if (dvalue .gt. 500.) dvalue = -(dvalue - 500.)
        pralt = tralt * 10. - dvalue * 10.
        press = 1013.25 * (1. - (pralt / 44331.))
     1    ** (1. / 0.190263)
        geopot = pralt
        if (geopot .lt. 0.) geopot = geopot + 500.
      else
        press = 0.
        geopot = 0.
        endif
      return
      end
c
c
c
      subroutine rdstuv (wdir,wspd,alat,alon)
      real alat,alon,wdir,wspd
c
c     convert wind direction and speed to u, v for RECCOs
c
      wdir = wdir * 10.
      wspd = 0.514444 * wspd
      call uvcomp (wdir,wspd)
      return
      end
c
c
c
      subroutine vortex(iyrs,imns,idys,line)
      character*4 itime1,itime2
      character*30 blank
      character*70 line
      logical knots
      dimension a(10,8)
      data blank /'                             '/
c
c
      i=0
      iwx=7
c
c
c     Read the day
c     ------------
      read(line(13:14),'(i2)') iday
      if(line(13:13).eq.' ') read (line(14:15),'(i2)') iday
      if (iday.gt.50) then
         iday=iday-50
         knots = .true.
         else
         knots = .false.
         endif
c
c     check for month, year flips
c     ---------------------------
      yy = iyrs
      mm = imns
      if (idys.gt.27.and.iday.eq.1) then
         mm = imns+1
         if (mm.eq.13) then
            mm = 1
            yy = iyrs + 1.
            if (yy.eq.100.) yy = 00.
            endif
         endif
c
      yymmdd = yy * 10000. + float(mm) * 100. + float(iday)
c
c
c
c     read the line with the mission number
c
 10   read(12,'(a)',end=99)line
c
c     check if the line has information
c
      if(line(1:30).eq.blank)goto 10
c
c     read the line 'SUPPLEMENTARY VORTEX DATA MESSAGE'
c
 20   read(12,'(a)',end=99)line
c
c     check if the line has the information
c
      if(line(1:30).eq.blank)goto 20
c
c     read the data line
c
3330  read(12,'(a)',end=99)line
c
c     check if line is data, remarks, or has no information at all
c
      if(index(line,'MF').ne.0.or.index(line,'REM').ne.0.or.
     1   line(1:30).eq.blank)goto 3330
c
c     if the line with the observation times has been read, write data
c
      if(index(line,'OB').ne.0)goto 90
c
c     count the number of the data point
c
      i=i+1
c
c     check for end of message
c
      if(index(line,'NNNN').ne.0)return
c
c     read the latitude
c
      read(line(3:5),'(i3)')ilat
c
c     read the longitude
c
      read(line(8:11),'(i4)')ilon
c
c     save the value of the latitude
c
      alat=ilat/10.
c
c     save the value of the longitude
c
      alon=ilon/10.
c
c     read the pressure level
c
      read(line(14:14),'(i1)',err=190)ihgt
c
c     save the value of the pressure
c
      if(ihgt.eq.5)press=400.
      if(ihgt.eq.4)press=500.
      if(ihgt.eq.3)press=700.
      if(ihgt.eq.2)press=850.
      if(ihgt.eq.1)press=1000.
      if(ihgt.eq.0)press=1070.
c
c     read the value of the geopotential
c
 190  read(line(15:17),'(i3)',err=191)ihgt
c
c     save the value of the geopotential
c
      geopot=ihgt
c
c     adjust the value of the geopotential
c
      if(press.eq.1070..and.geopot.lt.100)geopot=geopot+1000.
      if(press.eq.850.0)geopot=geopot+1000.
      if(press.eq.700.0)geopot=geopot+3000.
      if(press.eq.500.0)geopot=geopot+5000.
      if(press.eq.400.0)geopot=geopot+7000.
c
c     read the value of the temperature
c
 191  read(line(20:21),'(f2.0)',err=192)temp
c
c     correct the value of the temperature
c
      if(temp.ge.50.0)temp=50.0-temp
c
c     read the value of the dewpoint
c
 192  read(line(22:23),'(f2.0)',err=193)dewpt
c
c     correct the value of the dewpoint
c
      if(dewpt.ge.50.0)dewpt=50.0-dewpt
c
c     read the value of the wind direction
c
 193  read(line(25:26),'(i2)',err=194)iwdir
c
c     save the value of the wind direction
c
      wdir=iwdir
c
c     read the value of the wind speed
c
      read(line(27:29),'(i3)',err=194)iwspd
c
c     save the value of the wind speed
c
      wspd=iwspd
c
c     calculate the relative humidity
c
      call relhum (press,temp,dewpt,rh)
c
c     calculate the u and v components of the wind
c
      call rdstuv (wdir,wspd,alat,alon)
c
c     save values until time has been calculated
c
      a(i,1)=alat
      a(i,2)=alon
      a(i,3)=press
      a(i,4)=temp
      a(i,5)=rh
      a(i,6)=geopot
      a(i,7)=wdir
      a(i,8)=wspd
c
c     reinitialize values of the variables
c
 194  press=-99.
      temp=-99.
      rh=-99.
      geopot=-99.
      wdir=-99.
      wspd=-99.
c
c     continue loop
c
      goto 3330
c
c     read time of the first observation
c
 90   if(index(line,'SFC').ne.0)goto 3330
      idx1=index(line,'AT')
      if(idxchk.ne.2)then
        read(line(idx1+3:idx1+6),'(a)')itime1
        read(itime1(1:2),'(i2)')ihour1
        read(itime1(3:4),'(i2)')imin1
        jtime1=ihour1*60+imin1
      else
        read(line(idx1+3:idx1+6),'(a)')itime2
        read(itime2(1:2),'(i2)')ihour2
        read(itime2(3:4),'(i2)')imin2
        jtime2=ihour2*60+imin2
      endif
c
c     read time of the second observation
c
      if(idxchk.ne.2)then
        idx2=index(line(idx1+1:70),'AT')+idx1
        if(idx2.eq.0)then
          read(12,'(a)',end=99)line
          idxchk=2
          goto 90
        endif
        read(line(idx2+3:idx2+6),'(a)')itime2
        read(itime2(1:2),'(i2)')ihour2
        read(itime2(3:4),'(i2)')imin2
        jtime2=ihour2*60+imin2
      endif
      idxchk=1
c
c     calculate the times and write out the data
c
      do j=1,i
        itime=j*(jtime2-jtime1)/i+jtime1
        itime=itime/60*100+mod(itime,60)
        alat=a(j,1)
        alon=a(j,2)
        press=a(j,3)
        temp=a(j,4)
        rh=a(j,5)
        geopot=a(j,6)
        wdir=a(j,7)
        wspd=a(j,8)
        call out (iwx,yymmdd,itime,alat,alon,press,temp,rh,geopot,
     1            wdir,wspd,5)
      end do
      i=0
      goto 3330
99    continue
      stop
      end
c
c
c
      subroutine out (iwx,yymmdd,igmt,alat,alon,press,
     1                temp,rh,geopot,wdir,wspd,msgtype)
      character*4 tail
      real alat,alon,wdir,wspd
c
      common /output/nrecs

      tail='0000'
      if (msgtype.eq.1) tail = 'MANL'     ! DROP/Mandatory
      if (msgtype.eq.2) tail = 'SIGL'     ! DROP/Significant
      if (msgtype.eq.3) tail = 'ADDL'     ! DROP/Additional (51515)
      if (msgtype.eq.4) tail = 'RECO'     ! RECCO
      if (msgtype.eq.5) tail = 'SUPV'     ! SUPPL VTX
      if (msgtype.eq.6) tail = 'MWND'     ! DROP/Max wind
      nrecs = nrecs+1
c
c     print the data and write it to the output file
c
c      write (1,511) iwx,yymmdd,igmt,alat,alon,press,temp,rh,geopot,
c     1              wdir,wspd,tail
      write (13,510) iwx,yymmdd,igmt,alat,alon,press,temp,rh,geopot,
     1               wdir,wspd,tail
      return
510   format (i2,1x,f7.0,1x,i4,1x,2(f7.3,1x),3(f6.1,1x),
     1  f7.1,2(f6.1,1x),a4)
511   format (1x,i2,1x,f7.0,1x,i4,1x,2(f7.3,1x),3(f6.1,1x),
     1  f7.1,2(f6.1,1x),a4)
      end

C
C
C
C     -----------------------------------------------------
      SUBROUTINE XXAA(MODE_WMO)

	use thread_common

C 
C     Program writes XXAA portion of dropwindsonde code.
C
C     MODE_WMO = 1 for NOAA encoding
C              = 2 for USAF encoding
C     ----------------------------------------------------- 
C 
      PARAMETER (MXRC = 9000, ninvar = 12)
      PARAMETER (LNWMO = 66)
      PARAMETER (NDTMAX = 100)
      PARAMETER (NMANL = 12)
      PARAMETER (MXT = 10)
C 
      DIMENSION MONTH(12)
      DIMENSION PRMANL(NMANL)
      DIMENSION SDATA(NMANL,6)
      DIMENSION IXTROP(MXT)
C
      CHARACTER*1 BLANK
      CHARACTER*2 WMOB
      CHARACTER*4 STNS, ICAO
      CHARACTER*5 HEAD(5),LEVL(3,NMANL+MXT+1),NODAT,MRSBLK
      CHARACTER*20 MISSIONID
      CHARACTER*80 REMARK
      CHARACTER*(LNWMO) CODE(10) 
C
      LOGICAL WINDMAX, WINDTOP, TERMEXT, AFRES, ASTERISK, GOODWINDS
      LOGICAL AG_51515, AG_10190, AG_10191, AG_10166, AG_10167
C 
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
      COMMON /TEMPCTL/ MISSIONID, ICAO, RHMIN, REMARK
      COMMON /REGDATA/ AG_51515, AG_10190(2), AG_10191, AG_10166, 
     *                 AG_10167, PR_10190(2), GA_10190(2), PR_10191,
     *                 PZBTM, PZTOP, PDBTM(NDTMAX), PDTOP(NDTMAX),
     *                 NDT
C 
      DATA PRMANL/1070.,1000.,925.,850.,700.,500.,400.,
     *            300.,250.,200.,150.,100./
      DATA MONTH/31,28,31,30,31,30,31,31,30,31,30,31/ 
      DATA BLANK/' '/ 
      DATA NODAT/'*****'/ 
      DATA BAD/-999./ 
	DATA RHMIN/0.1/
C
C
C     Initialize variables
C     --------------------
      AFRES = .FALSE.
      IF (MODE_WMO.EQ.2) AFRES = .TRUE.
      ASTERISK = .FALSE.
C 
C 
C     Get day and hour of drop
C     ------------------------
      IDAY = IDY
      IF (MOD(IYR,4).EQ.0) MONTH(2) = 29
      IHOUR = NINT(HMSTS(TIML)/3600.) 
      IF (IHOUR.EQ.24) THEN 
         IHOUR = 0
         IDAY = IDAY+1
         IF (IDAY.GT.MONTH(IMO)) IDAY = 1 
         ENDIF
C 
C 
C     Convert WD,WS into U,V for interpolation
C     ----------------------------------------
100   DO 110 L = 1,IXEND
         X(L) = UCMP(SNDDAT(L,6),SNDDAT(L,7))
         Y(L) = VCMP(SNDDAT(L,6),SNDDAT(L,7))
110      CONTINUE 
C 
C 
C     Collect data at the mandatory levels.
C     First pick off surface data (or last level w/height)
C     ----------------------------------------------------
      IF (SNDDAT(IXEND,9).EQ.0.) THEN 
         TERMEXT = .FALSE.
         SDATA(1,1) = SNDDAT(IXEND,2) 
         SDATA(1,2) = SNDDAT(IXEND,9) 
         SDATA(1,3) = SNDDAT(IXEND,3) 
         RH = RHCHK(SNDDAT(IXEND,4),RHMIN)
         SDATA(1,4) = DEWPT(SNDDAT(IXEND,3),RH)
         SDATA(1,5) = SNDDAT(IXEND,6) 
         SDATA(1,6) = SNDDAT(IXEND,7) 
         CALL WINDZINT(10.,WD10,WS10,WQ10,M10)
         IF (WD10.GE.0.) THEN
            SDATA(1,5) = WD10
            SDATA(1,6) = WS10
            ENDIF
         ELSE 
         DO 150 I = 1,6 
              SDATA(1,I) = BAD
150           CONTINUE
         DO 160 I = IXEND,1,-1
            IF (SNDDAT(I,9).NE.BAD) THEN
               TERMEXT = .TRUE.
               TERMPR = SNDDAT(I,2) 
               TERMGA = SNDDAT(I,9) 
               TERMTE = SNDDAT(I,3) 
               TERMRH = SNDDAT(I,4)
               GOTO 175
               ENDIF
160         CONTINUE
         ENDIF
C
C 
C     Then interpolate for other mandatory levels 
C     ------------------------------------------- 
175   DO 200 I = 2,NMANL
         SDATA(I,1) = PRMANL(I)
         PNOW = SDATA(I,1)
         CALL POLATE(IXEND,SNDDAT(1,2),SNDDAT(1,9),PNOW,ZNOW,M,BAD) 
         CALL POLATE(IXEND,SNDDAT(1,2),SNDDAT(1,3),PNOW,TNOW,M,BAD) 
         CALL POLATE(IXEND,SNDDAT(1,2),SNDDAT(1,4),PNOW,HNOW,M,BAD) 
         CALL POLATE(IXEND,SNDDAT(1,2),X,PNOW,UNOW,M,BAD) 
         CALL POLATE(IXEND,SNDDAT(1,2),Y,PNOW,VNOW,M,BAD) 
         SDATA(I,2) = ZNOW
         SDATA(I,3) = TNOW
         SDATA(I,4) = DEWPT(TNOW,RHCHK(HNOW,RHMIN))
         SDATA(I,5) = WDCOMP(UNOW,VNOW) 
         SDATA(I,6) = WSCOMP(UNOW,VNOW) 
C 
C
C        Check if mand level falls between launch and first ODW pr
C        ---------------------------------------------------------
         IF (PNOW.LT.SNDDAT(1,2) .AND. PNOW.GE.FLOAT(NINT(PRL)) .AND.
     *   ABS(SNDDAT(1,2)-PRL).LE.50.) THEN 
              PRS = SNDDAT(1,2) 
              TNOW = POLATE2(PRL,PRS,TEL,SNDDAT(1,3),PNOW,BAD)
              HNOW = POLATE2(PRL,PRS,RHL,SNDDAT(1,4),PNOW,BAD)
              UNOW = POLATE2(PRL,PRS,UCMP(WDL,WSL),X(1),PNOW,BAD)
              VNOW = POLATE2(PRL,PRS,VCMP(WDL,WSL),Y(1),PNOW,BAD)
              SDATA(I,3) = TNOW 
              SDATA(I,4) = DEWPT(TNOW,RHCHK(HNOW,RHMIN)) 
              SDATA(I,5) = WDCOMP(UNOW,VNOW)
              SDATA(I,6) = WSCOMP(UNOW,VNOW)
C
C             If interpolation fails and launch is at a mand level
C             ----------------------------------------------------
              IF (NINT(PNOW).EQ.NINT(PRL)) THEN
                 IF (SDATA(I,3).EQ.BAD) SDATA(I,3) = TEL
                 IF (SDATA(I,4).EQ.BAD) THEN
                    HNOW = RHL
                    SDATA(I,4) = DEWPT(SDATA(I,3),RHCHK(HNOW,RHMIN))
                    ENDIF
                 IF (SDATA(I,5).EQ.BAD) SDATA(I,5) = WDL
                 IF (SDATA(I,6).EQ.BAD) SDATA(I,6) = WSL
                 ENDIF
C
C             Try and calculate height
C             ------------------------
              RH1 = RHCK(SNDDAT(1,4),.TRUE.,70.,BAD)
              RH2 = RHCK(HNOW,.TRUE.,70.,BAD)
              SDATA(I,2) = HYDROZ(PRS,SNDDAT(1,3),RH1,SNDDAT(1,9),
     *                     PNOW,SDATA(I,3),RH2,2,BAD)
              ENDIF
C 
200      CONTINUE 
C
C
C     Check if sfc pressure is between 950 and 1000 mb.
C     If so, estimate ht of 1000 mb surface
C     ----------------------------------------------------
      IF (SDATA(1,1).LT.1000.0 .AND. SDATA(1,1).GE.950.0) THEN
         SDATA(2,2) = HYDROZ(SDATA(1,1),SDATA(1,3),SDATA(1,4),
     *                SDATA(1,2),1000.,SDATA(1,3),SDATA(1,4),1,BAD)
         ENDIF
C
C 
C     For early terminations, extrap ht of nearest std sfc if 
C     within 25 mb.  
C     -------------------------------------------------------
      IF (.NOT. TERMEXT) GOTO 220
      DO 210 I = 2,NMANL
         PRDIFF = PRMANL(I)-TERMPR
         IF (PRDIFF.GT.0. .AND. PRDIFF.LE.25.) THEN
            CALL EXTRAP_TEMP(TERMPR,PRMANL(I),TX,BAD)
            RHX = RHCK(TERMRH,.TRUE.,70.,BAD)
            SDATA(I,2) = HYDROZ(TERMPR,TERMTE,RHX,TERMGA,
     *                   SDATA(I,1),TX,RHX,2,BAD)
            GOTO 220
            ENDIF
210      CONTINUE
C
C
C     Get information for regional groups, including sfc extrap.
C     ----------------------------------------------------------
220   CALL GET_51515(MODE_WMO)
      IF (AG_10191) THEN
         SDATA(1,1) = PR_10191
         SDATA(1,2) = 0.
         ENDIF
C
C
C     Find SDATA index of highest level with a wind: LVLHI
C     ----------------------------------------------------
230   LVLHI = -1
      DO 235 I = NMANL,1,-1 
         IF (SDATA(I,5).NE.BAD) THEN
              LVLHI = I 
              GOTO 240
              ENDIF 
235      CONTINUE 
C 
C 
C     Find SDATA index of highest level with THZ: LASTLVL 
C     ----------------------------------------------------
240   LASTLVL = -1
      DO 245 I = NMANL,1,-1 
         IF (SDATA(I,3).NE.BAD .OR. SDATA(I,4).NE.BAD 
     *       .OR. SDATA(I,2).NE.BAD) THEN 
              LASTLVL = I 
              GOTO 250
              ENDIF 
245      CONTINUE 
C 
C 
C     Search for level of maximum wind.  Peculiar expression
C     for the pressure is to ensure consistency between
C     this level and the level chosen in part B, which works
C     in log P coordinates and double precision
C     ------------------------------------------------------
250   WINDMAX = .FALSE. 
      WINDTOP = .TRUE.
      GOODWINDS = .FALSE.
      WSMAX = 0.0 
      DO 260 I = 1,IXEND
         IF (SNDDAT(I,2).GE.500. .OR. SNDDAT(I,2).EQ.BAD) GOTO 260
         IF (SNDDAT(I,6).EQ.BAD) GOTO 260
         GOODWINDS = .TRUE.
         IF (SNDDAT(I,7).GT.WSMAX .AND. SNDDAT(I,7).GE.31.) THEN
              WSMAX = SNDDAT(I,7) 
              WDMAX = SNDDAT(I,6) 
              PMAX  = EXP(LOG(DBLE(SNDDAT(I,2))))
              ZMAX  = SNDDAT(I,9)
              WINDMAX = .TRUE.
              WINDTOP = .FALSE. 
              ENDIF 
260   CONTINUE
C
      IF (PRL.NE.BAD .AND. PRL.LT.500. .AND. WSL.GE.31. .AND.
     *    WSL.GT.WSMAX .AND. GOODWINDS) THEN
         WSMAX = WSL
         WDMAX = WDL
         PMAX = EXP(LOG(DBLE(PRL)))
         ZMAX = HTL
         WINDMAX = .TRUE.
         WINDTOP = .TRUE.
         ENDIF
C
C
C     Calculate wind shear 3000 ft above and below max
C     ------------------------------------------------
      SHRA = BAD
      SHRB = BAD
      IF (.NOT.WINDMAX .OR. ZMAX.EQ.BAD) GOTO 270
      CALL WINDZINT(ZMAX+915.,WDA,WSA,WQA,MDUM)
      CALL WINDZINT(ZMAX-915.,WDB,WSB,WQB,MDUM)
      IF (WDA.NE.BAD .AND. WDMAX.NE.BAD) THEN
         USHRA = UCMP(WDA,WSA) - UCMP(WDMAX,WSMAX)
         VSHRA = VCMP(WDA,WSA) - VCMP(WDMAX,WSMAX)
         SHRA = ((USHRA**2.0 + VSHRA**2.0)**0.5)*1.94
         ENDIF
      IF (WDB.NE.BAD .AND. WDMAX.NE.BAD) THEN
         USHRB = UCMP(WDB,WSB) - UCMP(WDMAX,WSMAX)
         VSHRB = VCMP(WDB,WSB) - VCMP(WDMAX,WSMAX)
         SHRB = ((USHRB**2.0 + VSHRB**2.0)**0.5)*1.94
         ENDIF
C
C
C     Find tropopause
C     ---------------
270   CALL FINDTROP(NTROP,IXTROP)
C 
C 
C     Begin to encode message 
C     ----------------------- 
1     FORMAT(I1.1)
2     FORMAT(I2.2)
3     FORMAT(I3.3)
4     FORMAT(I4.4)
5     FORMAT(I5.5)
6     FORMAT(I6.6)
C 
C
C     Start with message header info
C     ------------------------------
300   HEAD(1)(1:5) = 'XXAA '
      WRITE(HEAD(2)(1:2),2) IDAY+50     
      WRITE(HEAD(2)(3:4),2) IHOUR 
C
      IF (LVLHI.EQ.-1) THEN
          HEAD(2)(5:5) = '/' 
          ELSE
          LVL = (SDATA(LVLHI,1)+1.)/100.
          IF (LVLHI.EQ.2) LVL = 0 
          WRITE(HEAD(2)(5:5),1) LVL 
          ENDIF
C
      LAT = NINT(ABS(RLAT*10))
      LON = NINT(ABS(RLON*10))
      IF (RLAT.GE.0. .AND. RLON.GE.0.) IQ = 7 
      IF (RLAT.LT.0. .AND. RLON.GE.0.) IQ = 5 
      IF (RLAT.LT.0. .AND. RLON.LT.0.) IQ = 3 
      IF (RLAT.GE.0. .AND. RLON.LT.0.) IQ = 1 
      HEAD(3)(1:2) = '99' 
      WRITE(HEAD(3)(3:5),3) LAT 
      WRITE(HEAD(4)(1:1),1) IQ
      WRITE(HEAD(4)(2:5),4) LON 
      LATN = LAT/10
      LONN = LON/10
      CALL MARSDEN(LATN,LONN,IQ,MRSBLK)
      HEAD(5)(1:5) = MRSBLK
C 
C 
C     Now start to code up the level data 
C     ----------------------------------- 
      DO 350 I = 1,NMANL
C 
C        Code pressure/height block 
C        -------------------------- 
         IF (I.EQ.1) THEN          ! Surface
              IP = 99 
              IHT = MOD(NINT(SDATA(1,1)),1000)
              WRITE(LEVL(1,I)(1:2),2) IP
              WRITE(LEVL(1,I)(3:5),3) IHT 
              IF (SDATA(1,1).EQ.BAD .OR. SDATA(1,2).NE.0.)
     *            LEVL(1,I)(3:5) = '///'
              ELSE                 ! All other levels 
              IP = SDATA(I,1)/10. 
              IF (IP.EQ.100) IP = 0 
C
C             Code up negative heights
C             ------------------------
              IF (SDATA(I,2).LT.0. .AND. SDATA(I,2).NE.BAD)
     *            SDATA(I,2) = 500. + ABS(SDATA(I,2))
C
              IHT = NINT(SDATA(I,2))
              IF (SDATA(I,1).LE.500.) IHT = NINT(SDATA(I,2)/10) 
              IHT = MOD(IHT,1000) 
              WRITE(LEVL(1,I)(1:2),2) IP
              IF (SDATA(I,1).EQ.BAD) LEVL(1,I)(1:2) = '//'
              WRITE(LEVL(1,I)(3:5),3) IHT 
              IF (SDATA(I,2).EQ.BAD) LEVL(1,I)(3:5) = '///' 
              ENDIF 
C 
C        Code temperature/dewpoint block
C        -------------------------------
         ITP = 2 * NINT(5.0*SDATA(I,3)) 
         IF (ITP.GE.0) GOTO 351 
         IF (SDATA(I,3)*10.0 .GE. FLOAT(ITP)) THEN
              ITP = ITP+1 
              ELSE
              ITP = ITP-1 
              ENDIF 
C 
351      DD = SDATA(I,3) - SDATA(I,4) 
         IF (DD.LE.5.) THEN 
              IDD = NINT(DD*10.)
              ELSE
              IDD = 50 + NINT(DD) 
              ENDIF 
         IF (IDD.GT.50 .AND. IDD.LE.55) IDD = 50
C
C        We will allow DPD > 80
C        ----------------------
         IF (IDD.GT.99) IDD = 99
C        IF (IDD.GT.80) IDD = 80
C 
         IF (SDATA(I,3).EQ.BAD .OR. SDATA(I,4).EQ.BAD) IDD = -999 
         WRITE(LEVL(2,I)(1:3),3) ABS(ITP) 
         IF (SDATA(I,3).EQ.BAD) LEVL(2,I)(1:3) = '///'
         WRITE(LEVL(2,I)(4:5),2) IDD
         IF (SDATA(I,4).EQ.BAD) LEVL(2,I)(4:5) = '//' 
C 
C        Code wind block
C        ---------------
         CALL WNDCODE(SDATA(I,5),SDATA(I,6),IWD,IWS)
         WRITE(LEVL(3,I)(1:2),2) IWD
         IF (SDATA(I,5).EQ.BAD) LEVL(3,I)(1:2) = '//' 
         WRITE(LEVL(3,I)(3:5),3) IWS
         IF (SDATA(I,6).EQ.BAD) LEVL(3,I)(3:5) = '///'
C
C        --------------------------------------------------------
C        Following WMO rules - sfc wind group is always included.
C        Otherwise, no wind groups are included after the last 
C        good wind.
C        --------------------------------------------------------
         IF (I.GT.LVLHI .AND. I.GT.1) LEVL(3,I) = NODAT 
C
C        Special rules for 200 and 100 mb.  Wind groups will be 
C        included for these levels if 250 or 150 mb, respectively,
C        has a wind group included.
C        ----------------------------------------------------------
         IF ((I.EQ.10 .AND. LVLHI.EQ.9) .OR. 
     *       (I.EQ.12 .AND. LVLHI.EQ.11))  LEVL(3,I) = '/////'
C
350      CONTINUE 
C
C
C     Code up tropopause
C     ------------------
      NLX = NMANL
      IF (NTROP.EQ.0) THEN
         NLX = NLX+1
         LEVL(1,NLX) = '88999' 
         LEVL(2,NLX) = NODAT
         LEVL(3,NLX) = NODAT
         ELSE
         DO 360 I = 1,NTROP
            NLX = NLX+1
            PTROP = SNDDAT(IXTROP(I),2)
            TTROP = SNDDAT(IXTROP(I),3)
            HTROP = SNDDAT(IXTROP(I),4)
            WDTROP = SNDDAT(IXTROP(I),6)
            WSTROP = SNDDAT(IXTROP(I),7)
            IP = MOD(NINT(PTROP),1000)
            CALL WNDCODE(WDTROP,WSTROP,IWD,IWS)
            LEVL(1,NLX)(1:2) = '88'
            WRITE(LEVL(1,NLX)(3:5),3) IP 
            IF (PTROP.EQ.BAD) LEVL(1,NLX)(3:5) = '///'
C
            ITP = 2 * NINT(5.0*TTROP) 
            IF (ITP.GE.0) GOTO 365
            IF (TTROP*10.0 .GE. FLOAT(ITP)) THEN
               ITP = ITP+1 
               ELSE
               ITP = ITP-1 
               ENDIF 
C 
365         RH = RHCHK(HTROP,RHMIN)
            DPTROP = DEWPT(TTROP,RH)
            DD = TTROP - DPTROP
            IF (DD.LE.5.) THEN 
               IDD = NINT(DD*10.)
               ELSE
               IDD = 50 + NINT(DD) 
               ENDIF 
            IF (IDD.GT.50 .AND. IDD.LE.55) IDD = 50
            IF (IDD.GT.99) IDD = 99
            IF (TTROP.EQ.BAD .OR. HTROP.EQ.BAD) IDD = -999 
            WRITE(LEVL(2,NLX)(1:3),3) ABS(ITP) 
            IF (TTROP.EQ.BAD) LEVL(2,NLX)(1:3) = '///'
            WRITE(LEVL(2,NLX)(4:5),2) IDD
            IF (HTROP.EQ.BAD) LEVL(2,NLX)(4:5) = '//' 
C
            WRITE(LEVL(3,NLX)(1:2),2) IWD
            WRITE(LEVL(3,NLX)(3:5),3) IWS
            IF (WDTROP.EQ.BAD) LEVL(3,NLX)(1:5) = '/////'
360         CONTINUE
         ENDIF
C 
C 
C     Code up max wind groups 
C     ----------------------- 
      NLX = NLX+1
      IF (WINDMAX) THEN 
         IP = MOD(NINT(PMAX),1000)
         CALL WNDCODE(WDMAX,WSMAX,IWD,IWS)
         IF (WINDTOP) THEN
              LEVL(1,NLX)(1:2) = '66' 
              ELSE
              LEVL(1,NLX)(1:2) = '77' 
              ENDIF 
         WRITE(LEVL(1,NLX)(3:5),3) IP 
         IF (PMAX.EQ.BAD) LEVL(1,NLX)(3:5) = '///'
         WRITE(LEVL(2,NLX)(1:2),2) IWD
         WRITE(LEVL(2,NLX)(3:5),3) IWS
         IF (SHRA.GT.99.) SHRA = 99.
         IF (SHRB.GT.99.) SHRB = 99.
         LEVL(3,NLX) = '4////'
         IF (SHRB.GE.0.) WRITE(LEVL(3,NLX)(2:3),2) NINT(SHRB)
         IF (SHRA.GE.0.) WRITE(LEVL(3,NLX)(4:5),2) NINT(SHRA)
         ELSE 
         LEVL(1,NLX) = '77999'
         LEVL(2,NLX) = NODAT
         LEVL(3,NLX) = NODAT
         ENDIF
C 
C 
C     ----------------------------------------------
C     Place individual blocks into final code array.
C     First do WMO abbreviated header, then XXAA.
C     ----------------------------------------------
500   ICNT = 1
      LINE = 1
      CODE(LINE)(1:6) = 'UZNT13'
      CALL MANOPAA(RLAT,RLON,WMOB)
      CODE(LINE)(ICNT+2:ICNT+3) = WMOB
c      IF (NPLTFORM.EQ.49) CODE(LINE)(ICNT+5:ICNT+5) = '4'
      CODE(LINE)(7:7) = BLANK
	do i = lnwmo-7, 1
		if (abrheader_thr(i:i) .ne. blank .and. 
	*		abrheader_thr(i:i) .ne. char(0)) goto 501
	end do
501	ml = i

      CODE(LINE)(8:ml+8) = abrheader_thr(1:ml)
	do i = ml+9, lnwmo
		code(line)(i:i) = blank
	end do
      LINE = LINE+1
      ICNT = 1
C      
      DO 530 I = 1,5
         CODE(LINE)(ICNT:ICNT+4) = HEAD(I)
         CODE(LINE)(ICNT+5:ICNT+5) = BLANK
         CALL INCR(ICNT,LINE,LNWMO) 
530      CONTINUE 
C 
C 
C     Now do level data 
C     ----------------- 
      ISTOP = LVLHI 
      IF (LASTLVL.GT.ISTOP) ISTOP = LASTLVL 
C 
      DO 545 IK = 1, NLX
         IF (IK.GT.ISTOP .AND. IK.LE.NMANL) GOTO 545
         DO 540  IL = 1,3 
              IF (LEVL(IL,IK).EQ.NODAT) GO TO 540 
              CODE(LINE)(ICNT:ICNT+4) = LEVL(IL,IK) 
              CODE(LINE)(ICNT+5:ICNT+5) = BLANK 
              CALL INCR(ICNT,LINE,LNWMO)
540           CONTINUE
545      CONTINUE 
C
C 
C     Add in regional data groups (51515).  Start on new line.
C     --------------------------------------------------------
      CALL CODE_51515(CODE,LINE,ICNT)
C
C
C     Add in mission 61616 Mission ID line
C     ------------------------------------
      CALL CODE_61616(CODE,LINE,ICNT)
C
C
C     All through...Blank out remaining code groups 
C     ---------------------------------------------
      IF (ICNT.EQ.1) THEN 
         LINE = LINE - 1
         ELSE 
         DO 590 IK = ICNT,LNWMO 
590           CODE(LINE)(IK:IK) = BLANK 
         ENDIF
C
C
C     Write out data to output file and check for bad characters
C     ----------------------------------------------------------
c600   WRITE(LUT,'(/," Mandatory levels: ",NN)') 
c      PINT = 0. 
c      CALL LVLINT(PINT,.TRUE.,.FALSE.) 
c      WRITE(LUT,*)
C
c      IF (NTROP.GT.0) THEN
c         DO 620 I = 1,NTROP
c            WRITE(LUT,'(" Found tropopause at ",f6.1," mb.")') 
c     *      SNDDAT(IXTROP(I),2)
c620         CONTINUE
c         WRITE(LUT,*)
c         ENDIF
C
c      DO 665 IJ = 1,LINE
c         DO 670 IK = 1,LNWMO
c            IF (CODE(IJ)(IK:IK).EQ.'*') ASTERISK = .TRUE.
c670         CONTINUE
c         WRITE(LUFO,'(A)',ERR=900) CODE(IJ) 
c665      CONTINUE 
C 
c      IF (ASTERISK) WRITE(LUT,'(/,A,
c     *" *** XXAA MESSAGE FLAWED - DO NOT TRANSMIT ***"/)') CHAR(7)
c      RETURN
C 
C 
C     Errors
C     ------
c900   WRITE(LUT,'(/," *** ERROR WRITING XXAA MESSAGE FILE ***",/)')
c      RETURN
C
c905   WRITE(LUT,'(/," *** INVALID NPLTFORM IN XXAA ***",/)')
c      CLOSE(LUFX)
c      RETURN
C 
c910   WRITE(LUT,'(/," *** ERROR ON AIRCRAFT CONTROL FILE ***",/)')
c      RETURN
       DO IJ = 1,LINE
		 nMsgLines_thr = nMsgLines_thr + 1
		 temp_msg_thr(nMsgLines_thr) = code(IJ)

         WRITE(LUT,'(1X,A)') CODE(IJ) 
       END DO
       END 
C 
C 
C 
C     ----------------------------------------------------- 
      SUBROUTINE XXBB(MODE_WMO)

	use thread_common

C 
C     Program writes XXBB portion of dropwindsonde code.
C     XXBB contains significant levels for T/H and wind.
C
C     MODE_WMO = 1 for NOAA encoding
C              = 2 for USAF encoding
C     ----------------------------------------------------- 
C 
      PARAMETER (MXRC = 9000, ninvar = 12)
      PARAMETER (LNWMO = 66)
      PARAMETER (NMANL = 12)
      PARAMETER (NDTMAX = 100)
      PARAMETER (NPROC = 50)
      PARAMETER (MXT = 10)
C 
      DIMENSION PRMANL(NMANL)
      DIMENSION MONTH(12) 
      DIMENSION IXTROP(MXT)
      DIMENSION IBGN(100),IEND(100) 
      DIMENSION ISIGT(MXRC+1), ISIGW(MXRC+1)
      REAL*8 PRDBL
C 
      LOGICAL SOK, INVERSION, INVPEND, AFRES, ASTERISK
      LOGICAL SUPADB
      LOGICAL AG_51515, AG_10190, AG_10191, AG_10166, AG_10167
C 
      CHARACTER*1 BLANK
      CHARACTER*4 STNS, ICAO
      CHARACTER*5 HEAD(5),LEVL(2,2*(MXRC+1)),NODAT,MRSBLK,SYSTEM(3)
      CHARACTER*20 MISSIONID
      CHARACTER*80 REMARK
      CHARACTER*(LNWMO) CODE((MXRC*2)/(LNWMO/12))   
C 
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /PROCESSING/ NPSTP, PROC(NPROC)
C
      COMMON /SIGLEVEL/ SGMET(MXRC+1,7), NT, PRDBL(MXRC+1)
      COMMON /SIGTYPES/ NSIGL,ISIGL(MXRC+1),ISIGTYPE(MXRC+1)
      COMMON /TEMPCTL/ MISSIONID, ICAO, RHMIN, REMARK
      COMMON /REGDATA/ AG_51515, AG_10190(2), AG_10191, AG_10166, 
     *                 AG_10167, PR_10190(2), GA_10190(2), PR_10191,
     *                 PZBTM, PZTOP, PDBTM(NDTMAX), PDTOP(NDTMAX),
     *                 NDT
C 
      DATA PRMANL/1070.,1000.,925.,850.,700.,500.,400.,
     *            300.,250.,200.,150.,100./
      DATA MONTH/31,28,31,30,31,30,31,31,30,31,30,31/ 
      DATA BLANK/' '/ 
      DATA NODAT/'*****'/ 
      DATA BAD/-999./ 
C
C
C     Do we follow standard NOAA or AFRES practices?
C     ----------------------------------------------
      AFRES = .FALSE.
      IF (MODE_WMO.EQ.2) AFRES = .TRUE.
      ASTERISK = .FALSE.
C 
C 
C     Set up local met arrays 
C     ----------------------- 
      CALL LOADSGMET
C
C
C     Get regional coding data (including extrapolations)
C     ---------------------------------------------------
      CALL GET_51515(MODE_WMO)
C
C
C     Determine the primary significant levels for T/H
C     ------------------------------------------------
C 
90    NSIGT = 0          ! # of significant T/H levels
      NSIGW = 0          ! # of significant wind levels 
      NSIGL = 0          ! # of real levels to plot
C 
C 
C     First primary level is the surface. If no sfc, send 9999.
C     --------------------------------------------------------
      IF (SNDDAT(IXEND,9).EQ.0.) THEN 
         ISX = NT 
         IF (SGMET(NT,1).EQ.BAD) ISX = 9999
         CALL INSERT(NSIGT,ISIGT,ISX,SOK,10) 
         ELSE 
         CALL INSERT(NSIGT,ISIGT,9999,SOK,10) 
         ENDIF
C 
C 
C     Tropopause level(s)
C     -------------------
      CALL FINDTROP(NTROP,IXTROP)
      IF (NTROP.GT.0) THEN
         DO 100 L = 1,NTROP
            ISX = IXTROP(L)+1
            CALL INSERT(NSIGT,ISIGT,ISX,SOK,18)
100         CONTINUE
         ENDIF
C
C
C     Look for all gaps in temperature
C     --------------------------------
      CALL GAP(SGMET(1,2),NT,IBGN,IEND,NGAPS,BAD) 
      DO 110 L = 1,NGAPS+1
         IB = IBGN(L) 
         IE = IEND(L) 
         CALL INSERT(NSIGT,ISIGT,IB,SOK,10)
         CALL INSERT(NSIGT,ISIGT,IE,SOK,10)
         IF (L.EQ.1) GOTO 110 
         CALL INSERT(NSIGT,ISIGT,IB-1,SOK,15)
110      CONTINUE 
C
C 
C     Identify bases and tops of inversions at least 20 mb thick
C     or of magnitude > 2.5 deg or 20% RH.
C     NC20 is the number of levels corresponding to 20 mb depth.
C     ----------------------------------------------------------
      NC20 = NELM(20.,FALLRATE,DATARATE)
      INVERSION = .FALSE.
      INVPEND = .FALSE.
      DO 120 J = 2,NT
         TECUR = SGMET(J,2)
         TELST = SGMET(J-1,2)
         IF (TECUR.EQ.BAD .OR. TELST.EQ.BAD) GOTO 120
C
C        Mark the top of inversion
C        -------------------------
         IF (TECUR.LE.TELST .AND. .NOT.INVERSION) THEN
            INVERSION = .TRUE.
C
C           If shallow inversion is still pending, don't reset top
C           ------------------------------------------------------
            IF (INVPEND) GOTO 120
            ITOP = J-1
            TTOP = TELST
            RHTOP = SGMET(J-1,3)
            ENDIF
C
C        Mark the bottom of inversion.  Only include it if the base
C        is below 300 mb or the first tropopause, whichever is higher
C        ------------------------------------------------------------
         IF ((TECUR.GT.TELST .OR. J.EQ.NT) .AND. INVERSION) THEN
            INVERSION = .FALSE.
            IBTM = J-1
            IF (TECUR.LE.TELST) IBTM = J
            TBTM = SGMET(IBTM,2)
            RHBTM = SGMET(IBTM,3)
            TDIFF = TTOP - TBTM
            RHDIFF = RHBTM - RHTOP
            IF (RHBTM.EQ.BAD .OR. RHTOP.EQ.BAD) RHDIFF = BAD
            IF (IBTM-ITOP.GE.NC20 .OR. TDIFF.GE.2.5 .OR.
     *          RHDIFF.GT.20) THEN
               CALL CHECKINV(IBTM-1,NTROP,IXTROP(1),INVOK)
               IF (INVOK.EQ.1) THEN
                  CALL INSERT(NSIGT,ISIGT,ITOP,SOK,13) 
                  CALL INSERT(NSIGT,ISIGT,IBTM,SOK,14) 
                  ENDIF
               INVPEND = .FALSE.
               ELSE
               INVPEND = .TRUE.
               ENDIF
            ENDIF
C
         IF (INVPEND .AND. TECUR.GT.TTOP) INVPEND = .FALSE.
120      CONTINUE 
C 
C
C     Identify bases and tops of superadiabatic layers at least
C     20 mb thick.  (NOAA ONLY)
C     ----------------------------------------------------------
      IF (MODE_WMO.EQ.2) GOTO 126
      SUPADB = .FALSE.
      TOL = 1.00
      DO 125 J = 2,NT
         TECUR = SGMET(J,2)
         TELST = SGMET(J-1,2)
         PRCUR = BAD
         IF (PRDBL(J).NE.BAD) PRCUR = EXP(PRDBL(J))
         PRLST = BAD
         IF (PRDBL(J-1).NE.BAD) PRLST = EXP(PRDBL(J-1))
         IF (TECUR.EQ.BAD .OR. TELST.EQ.BAD) GOTO 125
         IF (PRCUR.EQ.BAD .OR. PRLST.EQ.BAD) GOTO 125
         DRYLAPSE = BAD
         IF (PRCUR.NE.0.) DRYLAPSE = 0.2857143*(TECUR+273.16)/PRCUR
         DTDP = (TECUR-TELST)/(PRCUR-PRLST)
C
C        Mark the top of layer
C        ---------------------
         IF (DTDP.GT.DRYLAPSE*TOL .AND. .NOT.SUPADB) THEN
            SUPADB = .TRUE.
            ITOP = J-1
            TTOP = TELST
            RHTOP = SGMET(J-1,3)
            ENDIF
C
C        Mark the bottom of layer
C        ------------------------
         IF ((DTDP.LE.DRYLAPSE*TOL .OR. J.EQ.NT) .AND. SUPADB) THEN
            SUPADB = .FALSE.
            IBTM = J-1
            IF (J.EQ.NT) IBTM = J
            TBTM = SGMET(IBTM,2)
            RHBTM = SGMET(IBTM,3)
            IF (IBTM-ITOP.GE.NC20) THEN
               CALL INSERT(NSIGT,ISIGT,ITOP,SOK,16) 
               CALL INSERT(NSIGT,ISIGT,IBTM,SOK,17) 
               ENDIF
            ENDIF
C
125      CONTINUE 
126   CONTINUE
C 
C
C     Find first and last good humidity levels
C     ----------------------------------------
      CALL GAP(SGMET(1,3),NT,IBGN,IEND,NGAPS,BAD) 
      IB = IBGN(1)
      IE = IEND(NGAPS+1)
      CALL INSERT(NSIGT,ISIGT,IB,SOK,20) 
      CALL INSERT(NSIGT,ISIGT,IE,SOK,20) 
C
C
C     Identify bases of all cloud decks...defined as a minimum
C     of 15 mb  w/ RH>=95 above a minimum of 15 mb with RH<95.
C     NC15 is the number of levels corresponding to 15 mb depth.
C     ----------------------------------------------------------
      NC15 = NELM(15.,FALLRATE,DATARATE)
      DO 130 L = NC15,NT-NC15 
         DO 132 J = L-(NC15-1),L+NC15 
              IF (SGMET(J,3).EQ.BAD) GOTO 130 
C
C             RH is converted to be wrt/ice for T<0 for this purpose.
C             -------------------------------------------------------
              RH = SGMET(J,3)
              TE = SGMET(J,2)
              PR = BAD
              IF (PRDBL(J).NE.BAD) PR = EXP(PRDBL(J))
              IF (TE.LT.0. .AND. TE.NE.BAD) THEN
                 Q = 0.01*RH*QSATW(TE+273.16,PR)
                 RH = Q/QSATI(TE+273.16,PR)*100.
                 ENDIF
C
              IF (J.LE.L .AND. RH.LT.95.) GOTO 130
              IF (J.GT.L .AND. RH.GE.95.) GOTO 130
132           CONTINUE
         CALL INSERT(NSIGT,ISIGT,L,SOK,23) 
130      CONTINUE 
C 
C     Find highest and lowest T,H 
C     --------------------------- 
      CALL MINMAX(SGMET(1,2),NT,XMIN,XMAX,IN,IX)
      CALL INSERT(NSIGT,ISIGT,IX,SOK,11) 
      CALL INSERT(NSIGT,ISIGT,IN,SOK,11) 
C 
      CALL MINMAX(SGMET(1,3),NT,XMIN,XMAX,IN,IX)
      CALL INSERT(NSIGT,ISIGT,IX,SOK,21) 
      CALL INSERT(NSIGT,ISIGT,IN,SOK,21) 
C 
C 
C     Now find all secondary significant levels for temperature 
C     --------------------------------------------------------- 
      CALL SECSIG(NSIGT,ISIGT,2,MODE_WMO)
C 
C     Find secondary levels for humidity
C     ----------------------------------
      CALL SECSIG(NSIGT,ISIGT,3,MODE_WMO)
C
C
C     Do we have a level within 20 mb of surface?
C     -------------------------------------------
      IF (SNDDAT(IXEND,9).NE.0.) GOTO 150
      IF (PRDBL(ISIGT(1)).EQ.BAD .OR. PRDBL(ISIGT(2)).EQ.BAD) GOTO 150
      PDIFF = EXP(PRDBL(ISIGT(1))) - EXP(PRDBL(ISIGT(2)))
      IF (PDIFF.LE.20.) GOTO 150
C
C     No we don't.  Select one.
C     -------------------------
      CALL SIG20(NSIGT,ISIGT)
150   CONTINUE
C 
C 
C     Find primary significant levels for wind
C     ----------------------------------------
C 
C     First primary level is the surface. If no sfc, send 9999.
C     --------------------------------------------------------
      IF (SNDDAT(IXEND,9).EQ.0.) THEN 
         ISX = NT 
         IF (SGMET(NT,1).EQ.BAD) ISX = 9999
         CALL INSERT(NSIGW,ISIGW,ISX,SOK,30) 
         ELSE 
         CALL INSERT(NSIGW,ISIGW,9999,SOK,30) 
         ENDIF
C 
C     Now find first and last good wind levels, also max
C     --------------------------------------------------
      CALL GAP(SGMET(1,5),NT,IBGN,IEND,NGAPS,BAD) 
      IB = IBGN(1)
      IE = IEND(NGAPS+1)
      CALL INSERT(NSIGW,ISIGW,IB,SOK,30) 
      CALL INSERT(NSIGW,ISIGW,IE,SOK,30) 
      CALL MINMAX(SGMET(1,5),NT,XMIN,XMAX,IN,IX)
      CALL INSERT(NSIGW,ISIGW,IX,SOK,31) 
C 
C     Add in secondary wind levels...First check on speed.
C     ----------------------------------------------------
      CALL SECSIG(NSIGW,ISIGW,5,MODE_WMO)
C 
C     Now add in those levels due to direction
C     ----------------------------------------
      CALL SECSIG(NSIGW,ISIGW,4,MODE_WMO)
C
C 
C     ------------------- 
C     Now code up message 
C     ------------------- 
1     FORMAT(I1.1)
2     FORMAT(I2.2)
3     FORMAT(I3.3)
4     FORMAT(I4.4)
5     FORMAT(I5.5)
6     FORMAT(I6.6)
C
C 
C     Get day and hour of drop
C     ------------------------
300   IDAY = IDY
      IF (MOD(IYR,4).EQ.0) MONTH(2) = 29
      IHOUR = NINT(HMSTS(TIML)/3600.) 
      IF (IHOUR.EQ.24) THEN 
         IHOUR = 0
         IDAY = IDAY+1
         IF (IDAY.GT.MONTH(IMO)) IDAY = 1 
         ENDIF
C 
C 
C     Start with message header info, sounding system info
C     ----------------------------------------------------
      HEAD(1)(1:5) = 'XXBB '
      WRITE(HEAD(2)(1:2),2) IDAY+50     
      WRITE(HEAD(2)(3:4),2) IHOUR 
      HEAD(2)(5:5) = '/'
      LAT = NINT(ABS(RLAT*10))
      LON = NINT(ABS(RLON*10))
      IF (RLAT.GE.0. .AND. RLON.GE.0.) IQ = 7 
      IF (RLAT.LT.0. .AND. RLON.GE.0.) IQ = 5 
      IF (RLAT.LT.0. .AND. RLON.LT.0.) IQ = 3 
      IF (RLAT.GE.0. .AND. RLON.LT.0.) IQ = 1 
      HEAD(3)(1:2) = '99' 
      WRITE(HEAD(3)(3:5),3) LAT 
      WRITE(HEAD(4)(1:1),1) IQ
      WRITE(HEAD(4)(2:5),4) LON 
      LATN = LAT/10
      LONN = LON/10
      CALL MARSDEN(LATN,LONN,IQ,MRSBLK)
      HEAD(5)(1:5) = MRSBLK
C
      SYSTEM(1)(1:5) = '31313'
      SYSTEM(2)(1:1) = '0'
      SYSTEM(2)(2:3) = '96'
      IF (NSNDTYPE.EQ.1) SYSTEM(2)(4:5) = '05'
      IF (NSNDTYPE.EQ.2) SYSTEM(2)(4:5) = '08'
      SYSTEM(3)(1:1) = '8'
      WRITE(SYSTEM(3)(2:5),4) INT(TIML/100.)
C 
C 
C     Now start to code up the T/H sig levels 
C     --------------------------------------- 
      IF (NSIGT.LE.1) THEN
         NSIGT = 0
         GOTO 355 
         ENDIF              
C 
      IP = -11
      DO 350 I = 1,NSIGT
C 
C        Grab values of P,T,TD
C        ---------------------
         J = ISIGT(I) 
         IF (J.EQ.9999) THEN 
              PR = BAD
              TE = BAD
              TD = BAD
              DD = BAD
              ELSE
              PR = BAD
              IF (PRDBL(J).NE.BAD) PR = EXP(PRDBL(J))  ! EXP(SGMET(J,1))
              TE = SGMET(J,2) 
              TD = DEWPT(TE,RHCHK(SGMET(J,3),RHMIN)) 
C 
C             Is this level identifying a data gap? 
C             ------------------------------------- 
              IF (TE.EQ.BAD .AND. TD.EQ.BAD .AND. I.NE.1) PR = BAD
C 
              DD = TE - TD
              IF (TE.EQ.BAD .OR. TD.EQ.BAD) DD = BAD
              ENDIF 

C        Code indicator/pressure
C        -----------------------
         IP = IP + 11 
         IF (IP.GT.99) IP = 11  
         WRITE(LEVL(1,I)(1:2),2) IP 
         IF (IP.EQ.0 .AND. AG_10191) PR = PR_10191
         IF (PR.EQ.BAD) THEN
              LEVL(1,I)(3:5) = '///'
              ELSE
              IPR = MOD(NINT(PR),1000)
              WRITE(LEVL(1,I)(3:5),3) IPR 
              ENDIF 
C 
C        Code temperature/dewpoint block
C        -------------------------------
         ITP = 2 * NINT(5.0*TE) 
         IF (ITP.GE.0) GOTO 351 
         IF (TE*10.0 .GE. FLOAT(ITP)) THEN
              ITP = ITP+1 
              ELSE
              ITP = ITP-1 
              ENDIF 
351      IF (TE.EQ.BAD) THEN
              LEVL(2,I)(1:3) = '///'
              ELSE
              WRITE(LEVL(2,I)(1:3),3) ABS(ITP)
              ENDIF 
C 
         IF (DD.EQ.BAD) THEN
              LEVL(2,I)(4:5) = '//' 
              ELSE
              IF (DD.LE.5.) THEN
                 IDD = NINT(DD*10.) 
                 ELSE 
                 IDD = 50 + NINT(DD)
                 ENDIF
              IF (IDD.GT.50 .AND. IDD.LE.55) IDD = 50 
              IF (IDD.GT.99) IDD = 99 
C             IF (IDD.GT.80) IDD = 80 
              WRITE(LEVL(2,I)(4:5),2) IDD 
              ENDIF 
C 
350      CONTINUE 
C 
C 
C     Do we have any significant wind levels? 
C     --------------------------------------- 
355   IF (NSIGW.LE.1) THEN
         NSIGW = 0
         K = NSIGT
         GOTO 500 
         ENDIF
C 
C 
C     Separate T/H levels from wind levels
C     ------------------------------------
      LEVL(1,NSIGT+1) = '21212' 
      LEVL(2,NSIGT+1) = NODAT 
C 
C 
C     Now start to code up the wind sig levels
C     ----------------------------------------
      IP = -11
      DO 360 I = 1,NSIGW
C 
C        Grab values of P,WD,WS 
C        ---------------------- 
         K = NSIGT+1+I                 ! Index of LEVL
         J = ISIGW(I)                  ! Index of met arrays
         IF (J.EQ.9999) THEN 
              PR = BAD
              WD = BAD
              WS = BAD
              ELSE
              PR = BAD
              IF (PRDBL(J).NE.BAD) PR = EXP(PRDBL(J))  ! EXP(SGMET(J,1))
              WD = SGMET(J,4) 
              WS = SGMET(J,5) 
              ENDIF 
C
C        Code indicator/pressure.
C        ------------------------
         IP = IP + 11 
         IF (IP.GT.99) IP = 11  
         WRITE(LEVL(1,K)(1:2),2) IP 
         IF (IP.EQ.0 .AND. AG_10191) PR = PR_10191
         IF (PR.EQ.BAD) THEN
              LEVL(1,K)(3:5) = '///'
              ELSE
              IPR = MOD(NINT(PR),1000)
              WRITE(LEVL(1,K)(3:5),3) IPR 
              ENDIF 
C 
C        Code wind block
C        ---------------
         CALL WNDCODE(WD,WS,IWD,IWS)
         WRITE(LEVL(2,K)(1:2),2) IWD
         IF (WD.EQ.BAD) LEVL(2,K)(1:2) = '//' 
         WRITE(LEVL(2,K)(3:5),3) IWS
         IF (WS.EQ.BAD) LEVL(2,K)(3:5) = '///'
C 
360      CONTINUE 
C 
C 
C     ----------------------------------------------
C     Place individual blocks into final code array.
C     First do header.
C     ----------------------------------------------
500   NLEVELS = K 
      ICNT = 1
      LINE = 1
C 
      DO 530 I = 1,5
         CODE(LINE)(ICNT:ICNT+4) = HEAD(I)
         CODE(LINE)(ICNT+5:ICNT+5) = BLANK
         CALL INCR(ICNT,LINE,LNWMO) 
530      CONTINUE 
C
      IF (NLEVELS.EQ.0) GOTO 600
C 
C 
C     Now do level data 
C     ----------------- 
      DO 545 IK = 1, NLEVELS
         DO 540  IL = 1,2 
C 
C             If start of winds, skip to next line
C             ------------------------------------
              IF (IL.EQ.1 .AND. LEVL(IL,IK).EQ.'21212'
     *            .AND. ICNT.NE.1) THEN
                 DO 546 I = ICNT,LNWMO
546                 CODE(LINE)(I:I) = BLANK 
                 LINE = LINE+1
                 ICNT = 1 
                 ENDIF
C 
              IF (LEVL(IL,IK).EQ.NODAT) GO TO 540 
              CODE(LINE)(ICNT:ICNT+4) = LEVL(IL,IK) 
              CODE(LINE)(ICNT+5:ICNT+5) = BLANK 
              CALL INCR(ICNT,LINE,LNWMO)
540           CONTINUE
545      CONTINUE 
C 
C
C     Add in system status groups (31313).  Start on new line.
C     --------------------------------------------------------
600   IF (ICNT.EQ.1) THEN 
         LINE = LINE - 1
         ELSE 
         DO 610 I = ICNT,LNWMO
            CODE(LINE)(I:I) = BLANK 
610         CONTINUE
         ENDIF
C
      LINE = LINE+1
      ICNT = 1 
      DO 620 J = 1,3
         CODE(LINE)(ICNT:ICNT+4) = SYSTEM(J)(1:5)
         CODE(LINE)(ICNT+5:ICNT+5) = BLANK
         CALL INCR(ICNT,LINE,LNWMO)
620      CONTINUE   
C
C 
C     Add in additional data groups (51515).  Start on new line.
C     ----------------------------------------------------------
      CALL CODE_51515(CODE,LINE,ICNT)
C
C
C     Add in mission 61616 Mission ID line
C     ------------------------------------
      CALL CODE_61616(CODE,LINE,ICNT)
C
C
C     All through...Blank out remaining code groups 
C     ---------------------------------------------
700   IF (ICNT.EQ.1) THEN 
         LINE = LINE - 1
         ELSE 
         DO 750 IK = ICNT,LNWMO 
750           CODE(LINE)(IK:IK) = BLANK 
         ENDIF
C 
C 
C     Write out data to output file and check for bad characters
C     ----------------------------------------------------------
c      CALL PRINTSIG
c      WRITE(LUT,*)
C 
c      DO 765 IJ = 1,LINE
c         DO 770 IK = 1,LNWMO
c            IF (CODE(IJ)(IK:IK).EQ.'*') ASTERISK = .TRUE.
c770         CONTINUE
c         WRITE(LUFO,'(A)',ERR=900) CODE(IJ) 
c765      CONTINUE 
C 
c      IF (ASTERISK) WRITE(LUT,'(/,A,
c     *" *** XXBB MESSAGE FLAWED - DO NOT TRANSMIT ***"/)') CHAR(7)
c      RETURN
C 
C 
C     Errors
C
       DO IJ = 1,LINE
		 nMsgLines_thr = nMsgLines_thr + 1
		 temp_msg_thr(nMsgLines_thr) = code(IJ)

         WRITE(LUT,'(1X,A)') CODE(IJ) 
       END DO
c     ------
c900   WRITE(LUT,'(/," *** ERROR WRITING XXBB MESSAGE FILE ***",/)')
c      RETURN
C
c905   WRITE(LUT,'(/," *** INVALID NPLTFORM IN XXBB ***",/)')
c      CLOSE(LUFX)
c      RETURN
C
c910   WRITE(LUT,'(/," *** ERROR ON AIRCRAFT CONTROL FILE ***",/)')
c      RETURN
C 
      END 
C 
C 
C     =====================================================
C     ROUTINES FOR ENCODING WMO TEMP-DROP MESSAGE.
C     ENCODES PARTS A AND B OF MESSAGE.
C
C     Revisions:
C     05/24/96	Part B, section 7 (31313 group) added.
C     10/02/96  Part B: use 9999, not 999 for blank sig lvl
C     11/04/96  Part B: section 9 (51515 group) added, for
C                       10190 and 10191 indicators.  
C                       Extrapolated hts/prs only sent in
C                       part A for AFRES = .T.
C     12/24/96  Part A: Fix typo that would cause ht error
C                       if man level fell between launch 
C                       and first sonde pressure.
C     01/09/97  Part B: Adds 10166,67 for doubtful data.
C     ===================================================== 
C
C     -----------------------------------------------------
      SUBROUTINE UUAA 
C 
C     Program writes UUAA portion of dropwindsonde code.
C     ----------------------------------------------------- 
C 
C 
      use thread_common

      PARAMETER (MXRC = 9000, NINVAR = 12, MAXTMPMSG = 500)
      PARAMETER (LNWMO = 66)
      PARAMETER (NMANL = 12)
C 
      DIMENSION MONTH(12)
C 
      DIMENSION PRMANL(NMANL)
      DIMENSION SDATA(NMANL,6)
      CHARACTER*1 BLANK 
      CHARACTER*4 STNS
      CHARACTER*5 HEAD(6),LEVL(3,NMANL+1),NODAT,WINDX,MRSBLK
      CHARACTER*(LNWMO) CODE(6) 
      LOGICAL WINDMAX, WINDTOP, TERMEXT
C 
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
C 
      DATA PRMANL/1070.,1000.,925.,850.,700.,500.,400.,
     *            300.,250.,200.,150.,100./
      DATA MONTH/31,28,31,30,31,30,31,31,30,31,30,31/ 
      DATA BLANK/' '/ 
      DATA NODAT/'*****'/ 
      DATA BAD/-999./ 


C 
C 
C     Get day and hour of drop
C     ------------------------
      IDAY = IDY
      IF (MOD(IYR,4).EQ.0) MONTH(2) = 29
      IHOUR = NINT(HMSTS(TIML)/3600.) 
      IF (IHOUR.EQ.24) THEN 
         IHOUR = 0
         IDAY = IDAY+1
         IF (IDAY.GT.MONTH(IMO)) IDAY = 1 
         ENDIF
C 
C 
C     Convert WD,WS into U,V for interpolation
C     ----------------------------------------
100   DO 110 L = 1,IXEND
         X(L) = UCMP(SNDDAT(L,6),SNDDAT(L,7))
         Y(L) = VCMP(SNDDAT(L,6),SNDDAT(L,7))
110      CONTINUE 
C 
C 
C     ------------------------------------
C     Collect data at the mandatory levels
C     ------------------------------------
C 
C     First pick off surface data (or last level w/height)
C     ----------------------------------------------------
      IF (SNDDAT(IXEND,9).EQ.0. .OR. nsndtype .GE. 3) THEN 
         TERMEXT = .FALSE.
         SDATA(1,1) = SNDDAT(IXEND,2) 
         SDATA(1,2) = SNDDAT(IXEND,9) 
         SDATA(1,3) = SNDDAT(IXEND,3) 
         SDATA(1,4) = DEWPT(SNDDAT(IXEND,3),SNDDAT(IXEND,4))
         SDATA(1,5) = SNDDAT(IXEND,6) 
         SDATA(1,6) = SNDDAT(IXEND,7) 
         ELSE 
         DO 150 I = 1,6 
              SDATA(1,I) = BAD
150           CONTINUE
         DO 160 I = IXEND,1,-1
            IF (SNDDAT(I,9).NE.BAD) THEN
               TERMEXT = .TRUE.
               TERMPR = SNDDAT(I,2) 
               TERMGA = SNDDAT(I,9) 
               TERMTE = SNDDAT(I,3) 
               TERMRH = SNDDAT(I,4) 
               GOTO 175
               ENDIF
160         CONTINUE
         ENDIF
C
C 
C     Then interpolate for other mandatory levels 
C     ------------------------------------------- 
175   DO 200 I = 2,NMANL
         SDATA(I,1) = PRMANL(I)
         PNOW = SDATA(I,1)
         CALL POLATE(IXEND,SNDDAT(1,2),SNDDAT(1,9),PNOW,ZNOW,M,BAD) 
         CALL POLATE(IXEND,SNDDAT(1,2),SNDDAT(1,3),PNOW,TNOW,M,BAD) 
         CALL POLATE(IXEND,SNDDAT(1,2),SNDDAT(1,4),PNOW,HNOW,M,BAD) 
         CALL POLATE(IXEND,SNDDAT(1,2),X,PNOW,UNOW,M,BAD) 
         CALL POLATE(IXEND,SNDDAT(1,2),Y,PNOW,VNOW,M,BAD) 
         SDATA(I,2) = ZNOW
         SDATA(I,3) = TNOW
         SDATA(I,4) = DEWPT(TNOW,HNOW)
         SDATA(I,5) = WDCOMP(UNOW,VNOW) 
         SDATA(I,6) = WSCOMP(UNOW,VNOW) 
C 
C
C        Check if mand level falls between launch and first ODW pr
C        ---------------------------------------------------------
         IF (PNOW.LT.SNDDAT(1,2) .AND. PNOW.GE.PRL .AND.
     *        ABS(SNDDAT(1,2)-PRL).LE.50.) THEN 
              PRS = SNDDAT(1,2) 
              TNOW = POLATE2(PRL,PRS,TEL,SNDDAT(1,3),PNOW,BAD)
              HNOW = POLATE2(PRL,PRS,RHL,SNDDAT(1,4),PNOW,BAD)
              UNOW = POLATE2(PRL,PRS,UCMP(WDL,WSL),X(1),PNOW,BAD)
              VNOW = POLATE2(PRL,PRS,VCMP(WDL,WSL),Y(1),PNOW,BAD)
              SDATA(I,2) = HYDROZ(PRS,SNDDAT(1,3),SNDDAT(1,4),
     *                     SNDDAT(1,9),PNOW,TNOW,HNOW,2,BAD)
              SDATA(I,3) = TNOW 
              SDATA(I,4) = DEWPT(TNOW,HNOW) 
              SDATA(I,5) = WDCOMP(UNOW,VNOW)
              SDATA(I,6) = WSCOMP(UNOW,VNOW)
              ENDIF 
C 
200      CONTINUE 
C
C
C     Check if sfc pressure is between 950 and 1000 mb.
C     If so, estimate ht of 1000 mb surface
C     ----------------------------------------------------
      IF (SDATA(1,1).LT.1000.0 .AND. SDATA(1,1).GE.950.0) THEN
	 SDATA(2,2) = HYDROZ(SDATA(1,1),SDATA(1,3),SDATA(1,4),
     *                SDATA(1,2),1000.,SDATA(1,3),SDATA(1,4),1,BAD)
         ENDIF
C
C 
C     The following is done in PART A ONLY for AFRES messages:
C     For early terminations, extrap ht of nearest std sfc if 
C     within 25 mb.  Also extrap sfcp if term below 850 mb.
C     -------------------------------------------------------
C
C
C     Find SDATA index of highest level with a wind: LVLHI
C     ----------------------------------------------------
230   LVLHI = -1
      DO 235 I = NMANL,1,-1 
         IF (SDATA(I,5).NE.BAD) THEN
              LVLHI = I 
              GOTO 240
              ENDIF 
235      CONTINUE 
C 
C 
C     Find SDATA index of highest level with THZ: LASTLVL 
C     ----------------------------------------------------
240   LASTLVL = -1
      DO 245 I = NMANL,1,-1 
         IF (SDATA(I,3).NE.BAD .OR. SDATA(I,4).NE.BAD 
     *       .OR. SDATA(I,2).NE.BAD) THEN 
              LASTLVL = I 
              GOTO 250
              ENDIF 
245      CONTINUE 
C 
C 
C     Search for level of maximum wind
C     --------------------------------
250   WINDMAX = .FALSE. 
      WINDTOP = .TRUE.
      WSMAX = 0.0 
      DO 260 I = 1,IXEND
         IF (SNDDAT(I,2).GE.500.) GOTO 260
         IF (SNDDAT(I,7).GT.WSMAX .AND. SNDDAT(I,7).GE.31.) THEN
              WSMAX = SNDDAT(I,7) 
              WDMAX = SNDDAT(I,6) 
              PMAX  = SNDDAT(I,2) 
              WINDMAX = .TRUE.
              WINDTOP = .FALSE. 
              IF (I.EQ.1) WINDTOP = .TRUE.
              ENDIF 
260   CONTINUE
C 
C 
C     Begin to encode message 
C     ----------------------- 
1     FORMAT(I1.1)
2     FORMAT(I2.2)
3     FORMAT(I3.3)
4     FORMAT(I4.4)
5     FORMAT(I5.5)
6     FORMAT(I6.6)
C 
C     Start with message header info
C     ------------------------------
300   HEAD(1)(1:5) = 'UUAA '
	  HEAD(2)(1:5) = 'CGDX '
      WRITE(HEAD(3)(1:2),2) IDAY+50     
      WRITE(HEAD(3)(3:4),2) IHOUR 
C
      IF (LVLHI.EQ.-1) THEN
          HEAD(3)(5:5) = '/' 
          ELSE
          LVL =(SDATA(LVLHI,1)+1.0)/100. 
          IF (LVLHI.EQ.2) LVL = 0 
          WRITE(HEAD(3)(5:5),1) LVL 
          ENDIF
C
      LAT = NINT(ABS(RLAT*10))
      LON = NINT(ABS(RLON*10))
      IF (RLAT.GE.0. .AND. RLON.GE.0.) IQ = 7 
      IF (RLAT.LT.0. .AND. RLON.GE.0.) IQ = 5 
      IF (RLAT.LT.0. .AND. RLON.LT.0.) IQ = 3 
      IF (RLAT.GE.0. .AND. RLON.LT.0.) IQ = 1 
      HEAD(4)(1:2) = '99' 
      WRITE(HEAD(4)(3:5),3) LAT 
      WRITE(HEAD(5)(1:1),1) IQ
      WRITE(HEAD(5)(2:5),4) LON 
      LATN = LAT/10
      LONN = LON/10
      CALL MARSDEN(LATN,LONN,IQ,MRSBLK)
      HEAD(6)(1:5) = MRSBLK
C 
C 
C     Now start to code up the level data 
C     ----------------------------------- 
      DO 350 I = 1,NMANL
C 
C        Code pressure/height block 
C        -------------------------- 
         IF (I.EQ.1) THEN          ! Surface
              IP = 99 
	        IHT = MOD(NINT(SDATA(1,1)),1000)
              WRITE(LEVL(1,I)(1:2),2) IP
              WRITE(LEVL(1,I)(3:5),3) IHT 
              IF (SDATA(1,1).EQ.BAD)
     *            LEVL(1,I)(3:5) = '///'
         ELSE                 ! All other levels 
              IP = SDATA(I,1)/10. 
              IF (IP.EQ.100) IP = 0 
C
C             Code up negative heights
C             ------------------------
              IF (SDATA(I,2).LT.0. .AND. SDATA(I,2).NE.BAD)
     *            SDATA(I,2) = 500. + ABS(SDATA(I,2))
C
              IHT = NINT(SDATA(I,2))
              IF (SDATA(I,1).LE.500.) IHT = NINT(SDATA(I,2)/10) 
              IHT = MOD(IHT,1000) 
              WRITE(LEVL(1,I)(1:2),2) IP
              IF (SDATA(I,1).EQ.BAD) LEVL(1,I)(1:2) = '//'
              WRITE(LEVL(1,I)(3:5),3) IHT 
              IF (SDATA(I,2).EQ.BAD) LEVL(1,I)(3:5) = '///' 
         ENDIF 
C 
C        Code temperature/dewpoint block
C        -------------------------------
         ITP = 2 * NINT(5.0*SDATA(I,3)) 
         IF (ITP.GE.0) GOTO 351 
         IF (SDATA(I,3)*10.0 .GE. FLOAT(ITP)) THEN
              ITP = ITP+1 
              ELSE
              ITP = ITP-1 
              ENDIF 
C 
351      DD = SDATA(I,3) - SDATA(I,4) 
         IF (DD.LE.5.) THEN 
              IDD = NINT(DD*10.)
              ELSE
              IDD = 50 + NINT(DD) 
              ENDIF 
         IF (IDD.GT.50 .AND. IDD.LE.55) IDD = 50
C
C        We will allow DPD > 80
C        ----------------------
         IF (IDD.GT.99) IDD = 99
C        IF (IDD.GT.80) IDD = 80
C 
         IF (SDATA(I,3).EQ.BAD .OR. SDATA(I,4).EQ.BAD) IDD = -999 
         WRITE(LEVL(2,I)(1:3),3) ABS(ITP) 
         IF (SDATA(I,3).EQ.BAD) LEVL(2,I)(1:3) = '///'
         WRITE(LEVL(2,I)(4:5),2) IDD
         IF (SDATA(I,4).EQ.BAD) LEVL(2,I)(4:5) = '//' 
C 
C        Code wind block
C        ---------------
         CALL WNDCODE(SDATA(I,5),SDATA(I,6),IWD,IWS)
         WRITE(LEVL(3,I)(1:2),2) IWD
         IF (SDATA(I,5).EQ.BAD) LEVL(3,I)(1:2) = '//' 
         WRITE(LEVL(3,I)(3:5),3) IWS
         IF (SDATA(I,6).EQ.BAD) LEVL(3,I)(3:5) = '///'
C
C        --------------------------------------------------------
C        Following WMO rules - sfc wind group is always included.
C        Otherwise, no wind groups are included after the last 
C        good wind.
C        --------------------------------------------------------
         IF (I.GT.LVLHI .AND. I.GT.1) LEVL(3,I) = NODAT 
C
C        The following commented line would instead omit the sfc
C        wind if no other level had a wind.  This is how the AFRES
C        does it, but is not kosher.
C        ---------------------------------------------------------
C         IF (I.GT.LVLHI) LEVL(3,I) = NODAT
C 
C        Special rules for 200 and 100 mb.  Wind groups will be 
C        included for these levels if 250 or 150 mb, respectively,
C        has a wind group included.
C        ----------------------------------------------------------
	 IF ((I.EQ.10 .AND. LVLHI.EQ.9) .OR. 
     *       (I.EQ.12 .AND. LVLHI.EQ.11))  LEVL(3,I) = '/////'
C
350      CONTINUE 
C 
C 
C     Code up max wind groups 
C     ----------------------- 
C 
      LEVL(1,NMANL+1) = '88999' 
      IF (WINDMAX) THEN 
C        IP = NINT(AMOD(PMAX,1000.))
	 IP = MOD(NINT(PMAX),1000)
         CALL WNDCODE(WDMAX,WSMAX,IWD,IWS)
         IF (WINDTOP) THEN
              LEVL(2,NMANL+1)(1:2) = '66' 
              ELSE
              LEVL(2,NMANL+1)(1:2) = '77' 
              ENDIF 
         WRITE(LEVL(2,NMANL+1)(3:5),3) IP 
         WRITE(LEVL(3,NMANL+1)(1:2),2) IWD
         WRITE(LEVL(3,NMANL+1)(3:5),3) IWS
         WINDX = '4////'
         ELSE 
         LEVL(2,NMANL+1) = '77999'
         LEVL(3,NMANL+1) = NODAT
         WINDX = NODAT
         ENDIF
C 
C 
C     ----------------------------------------------
C     Place individual blocks into final code array.
C     First do header.
C     ----------------------------------------------
C 
      ICNT = 1
      LINE = 1
      DO 530 I = 1,6
         CODE(LINE)(ICNT:ICNT+4) = HEAD(I)
         CODE(LINE)(ICNT+5:ICNT+5) = BLANK
         CALL INCR(ICNT,LINE,LNWMO) 
530      CONTINUE 
C 
C 
C     Now do level data 
C     ----------------- 
C 
      ISTOP = LVLHI 
      IF (LASTLVL.GT.ISTOP) ISTOP = LASTLVL 
C 
      DO 545 IK = 1, NMANL+1
         IF (IK.GT.ISTOP .AND. IK.NE.(NMANL+1)) GOTO 545
         DO 540  IL = 1,3 
              IF (LEVL(IL,IK).EQ.NODAT) GO TO 540 
              CODE(LINE)(ICNT:ICNT+4) = LEVL(IL,IK) 
              CODE(LINE)(ICNT+5:ICNT+5) = BLANK 
              CALL INCR(ICNT,LINE,LNWMO)
540           CONTINUE
545      CONTINUE 
C 
C     Tack on last wind group 
C     ----------------------- 
C 
      IF (WINDX.NE.NODAT) THEN
         CODE(LINE)(ICNT:ICNT+4) = WINDX
         CODE(LINE)(ICNT+5:ICNT+5) = BLANK
         CALL INCR(ICNT,LINE,LNWMO) 
         ENDIF
C
C	add trailing =
C
C	Al Mungeon authorizes turning off the trailing = sign
C   on UUAA part
C	CODE(LAST_LINE)(LAST_ICNT+5:LAST_ICNT+5) = '='
C 
C 
C     All through...
C     Blank out remaining code groups 
C     ------------------------------- 
C 
      IF (ICNT.EQ.1) THEN 
         LINE = LINE - 1
         ELSE 
         DO 550 IK = ICNT,LNWMO 
550           CODE(LINE)(IK:IK) = BLANK 
         ENDIF
C 
C 
C     Write out data to output file 
C     ----------------------------- 
C 
      WRITE(LUT,'(//," Mandatory levels: ",$)') 
      PINT = 0. 
      CALL LVLINT(PINT,.TRUE.,.FALSE.) 
      WRITE(LUT,*)
      DO 560 IJ = 1,LINE
		 nMsgLines_thr = nMsgLines_thr + 1
		 temp_msg_thr(nMsgLines_thr) = code(IJ)

         WRITE(LUT,'(1X,A)') CODE(IJ) 
560      CONTINUE 
C 
c      DO 565 IJ = 1,LINE
c         WRITE(LUFO,'(A)',ERR=900) CODE(IJ) 
565      CONTINUE 
C 
      RETURN
C 
C 
900   WRITE(LUT,'(/," *** ERROR WRITING UUAA MESSAGE ***",/)')
      RETURN
      END 
C 
C 
C 
C     ----------------------------------------------------- 
      SUBROUTINE UUBB
C 
C     Program writes UUBB portion of dropwindsonde code.
C     UUBB contains significant levels for T/H and wind.
C     ----------------------------------------------------- 
C 
	  use thread_common

      PARAMETER (MXRC = 9000, NINVAR = 12, MAXTMPMSG = 500)
      PARAMETER (LNWMO = 66)
      PARAMETER (NMANL = 12)
      PARAMETER (NDTMAX = 100)
      PARAMETER (NPROC = 50)
C 


C      DIMENSION PRMANL(NMANL)
      DIMENSION MONTH(12) 
      DIMENSION IBGN(100),IEND(100) 
      DIMENSION ISIGT(MXRC+1), ISIGW(MXRC+1)
C      DIMENSION PR_10190(2), GA_10190(2)
C      DIMENSION PDTOP(NDTMAX),PDBTM(NDTMAX)
      REAL*8 PRDBL
C 
C      LOGICAL SOK, INVERSION, INVPEND, AFRES, TERMEXT
      LOGICAL SOK, INVERSION, INVPEND
C      LOGICAL AG_51515, AG_10190(2), AG_10191, AG_10166, AG_10167
C 
      CHARACTER*1 BLANK 
      CHARACTER*4 STNS
      CHARACTER*5 HEAD(6),LEVL(2,2*(MXRC+1)),NODAT,MRSBLK,SYSTEM(3)
      CHARACTER*(LNWMO) CODE((MXRC*2)/(LNWMO/12))   
C 
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /SIGLEVEL/ SGMET(MXRC+1,7), NT, PRDBL(MXRC+1)
      COMMON /SIGTYPES/ NSIGL,ISIGL(MXRC+1),ISIGTYPE(MXRC+1)
      COMMON /PROCESSING/ NPSTP, PROC(NPROC)
C 
C      DATA PRMANL/1070.,1000.,925.,850.,700.,500.,400.,
C     *            300.,250.,200.,150.,100./
      DATA MONTH/31,28,31,30,31,30,31,31,30,31,30,31/ 
      DATA BLANK/' '/ 
      DATA NODAT/'*****'/ 
      DATA BAD/-999./ 
C      DATA AFRES/.TRUE./

C 
C 
C     Set up local met arrays 
C     ----------------------- 
      CALL LOADSGMET
C
C
C     Determine first and last height data for possible extrapolation
C     ---------------------------------------------------------------
C      TERMEXT = .FALSE.
C      AG_10190(1) = .FALSE.
C      AG_10190(2) = .FALSE.
C      AG_10191 = .FALSE.
C
C      IF (SNDDAT(IXEND,9).EQ.0 .OR. nsndtype .GE. 3) GOTO 55
C      DO 50 I = IXEND,1,-1
C         IF (SNDDAT(I,9).NE.BAD) THEN
C            TERMEXT = .TRUE.
C           TERMPR = SNDDAT(I,2) 
C           TERMGA = SNDDAT(I,9) 
C            TERMTE = SNDDAT(I,3) 
C            TERMRH = SNDDAT(I,4) 
C            GOTO 55
C            ENDIF
C50       CONTINUE
C
55    DO 60 I = 1,IXEND
         IF (SNDDAT(I,9).NE.BAD .AND. SNDDAT(I,9).NE.0.) THEN
            TOPPR = SNDDAT(I,2) 
            TOPGA = SNDDAT(I,9) 
            TOPTE = SNDDAT(I,3) 
            TOPRH = SNDDAT(I,4) 
            GOTO 65
            ENDIF
60       CONTINUE
C
C
C     Do we extrapolate upwards
C     -------------------------
65		continue
C65    DO 70 I = 2,NMANL
C         PRDIFF = PRL-PRMANL(I)
C         IF (PRDIFF.GT.0. .AND. PRDIFF.LE.25.) THEN
C            CALL EXTRAP_TEMP(TOPPR,PRMANL(I),TX,BAD)
C            GA = HYDROZ(TOPPR,TOPTE,TOPRH,TOPGA,
C     *                  PRMANL(I),TX,TOPRH,2,BAD)
C            IF (GA.NE.BAD) THEN
C               AG_10190(1) = .TRUE.
C               PR_10190(1) = PRMANL(I)
C               GA_10190(1) = GA
C               ENDIF
C            GOTO 75
C            ENDIF
C70       CONTINUE
C
C
C     Do we extrapolate downwards
C     ---------------------------
C75    IF (.NOT. TERMEXT) GOTO 90
C      DO 80 I = 2,NMANL
C         PRDIFF = PRMANL(I)-TERMPR
C         IF (PRDIFF.GT.0. .AND. PRDIFF.LE.25.) THEN
C            CALL EXTRAP_TEMP(TERMPR,PRMANL(I),TX,BAD)
C            GA = HYDROZ(TERMPR,TERMTE,TERMRH,TERMGA,
C     *                  PRMANL(I),TX,TERMRH,2,BAD)
C            IF (GA.NE.BAD) THEN
C               AG_10190(2) = .TRUE.
C               PR_10190(2) = PRMANL(I)
C               GA_10190(2) = GA
C              ENDIF
C            GOTO 85
C            ENDIF
C80       CONTINUE
C
C
C85    IF (TERMPR.LT.850.) GOTO 90
C      ESTSFCP = HYDROP(TERMPR,TERMTE,TERMRH,TERMGA,
C     *                 TERMTE,TERMRH,0.,2,BAD)
C      CALL EXTRAP_TEMP(TERMPR,ESTSFCP,TX,BAD)
C      ESTSFCP = HYDROP(TERMPR,TERMTE,TERMRH,TERMGA,
C     *                 TX,TERMRH,0.,2,BAD)
C      IF (ESTSFCP.NE.BAD) THEN
C         AG_10191 = .TRUE.
C         PR_10191 = ESTSFCP
C         ENDIF
C
C
C 
C     ------------------------------------------------
C     Determine the primary significant levels for T/H
C     ------------------------------------------------
C 
90    NSIGT = 0          ! # of significant T/H levels
      NSIGW = 0          ! # of significant wind levels 
      NSIGL = 0          ! # of real levels to plot
C 
C 
C     First primary level is the surface. If no sfc, send 9999.
C     --------------------------------------------------------
C      IF (SNDDAT(IXEND,9).EQ.0. .OR. nsndtype .GE. 3) THEN 
         ISX = NT 
         IF (SGMET(NT,1).EQ.BAD) ISX = 9999
         CALL INSERT(NSIGT,ISIGT,ISX,SOK,10) 
C         ELSE 
C         CALL INSERT(NSIGT,ISIGT,9999,SOK,10) 
C         ENDIF
C 
C 
C     Look for all gaps in temperature
C     --------------------------------
      CALL GAP(SGMET(1,2),NT,IBGN,IEND,NGAPS,BAD) 
      DO 100 L = 1,NGAPS+1
         IB = IBGN(L) 
         IE = IEND(L) 
         CALL INSERT(NSIGT,ISIGT,IB,SOK,10)
         CALL INSERT(NSIGT,ISIGT,IE,SOK,10)
         IF (L.EQ.1) GOTO 100 
         CALL INSERT(NSIGT,ISIGT,IB-1,SOK,10)
100      CONTINUE 
C 
C     Identify bases and tops of inversions at least 20 mb thick
C     or of magnitude > 2.5 deg or 20% RH.
C     NC20 is the number of levels corresponding to 20 mb depth.
C     ----------------------------------------------------------
      NC20 = NELM(20.,FALLRATE,DATARATE)
      INVERSION = .FALSE.
      INVPEND = .FALSE.
      DO 120 J = 2,NT
         TECUR = SGMET(J,2)
         TELST = SGMET(J-1,2)
         IF (TECUR.EQ.BAD .OR. TELST.EQ.BAD) GOTO 120
C
C        Mark the top of inversion
C        -------------------------
         IF (TECUR.LE.TELST .AND. .NOT.INVERSION) THEN
            INVERSION = .TRUE.
C
C           If shallow inversion is still pending, don't reset top
C           ------------------------------------------------------
            IF (INVPEND) GOTO 120
            ITOP = J-1
            TTOP = TELST
            RHTOP = SGMET(J-1,3)
            ENDIF
C
C        Mark the bottom of inversion
C        ----------------------------
         IF ((TECUR.GT.TELST .OR. J.EQ.NT) .AND. INVERSION) THEN
            INVERSION = .FALSE.
            IBTM = J-1
            IF (TECUR.LE.TELST) IBTM = J
            TBTM = SGMET(IBTM,2)
            RHBTM = SGMET(IBTM,3)
            TDIFF = TTOP - TBTM
            RHDIFF = RHBTM - RHTOP
            IF (RHBTM.EQ.BAD .OR. RHTOP.EQ.BAD) RHDIFF = BAD
            IF (IBTM-ITOP.GE.NC20 .OR. TDIFF.GE.2.5 .OR.
     *          RHDIFF.GT.20) THEN
               CALL INSERT(NSIGT,ISIGT,ITOP,SOK,13) 
               CALL INSERT(NSIGT,ISIGT,IBTM,SOK,14) 
               INVPEND = .FALSE.
               ELSE
               INVPEND = .TRUE.
               ENDIF
            ENDIF
C
         IF (INVPEND .AND. TECUR.GT.TTOP) INVPEND = .FALSE.
120      CONTINUE 
C 
C
C     Find first and last good humidity levels
C     ----------------------------------------
      CALL GAP(SGMET(1,3),NT,IBGN,IEND,NGAPS,BAD) 
      IB = IBGN(1)
      IE = IEND(NGAPS+1)
      CALL INSERT(NSIGT,ISIGT,IB,SOK,20) 
      CALL INSERT(NSIGT,ISIGT,IE,SOK,20) 
C
C
C     Identify bases of all cloud decks...defined as a minimum
C     of 15 mb  w/ RH>=95 above a minimum of 15 mb with RH<95.
C     NC15 is the number of levels corresponding to 15 mb depth.
C     ----------------------------------------------------------
      NC15 = NELM(15.,FALLRATE,DATARATE)
      DO 130 L = NC15,NT-NC15 
         DO 132 J = L-(NC15-1),L+NC15 
              IF (SGMET(J,3).EQ.BAD) GOTO 130 
C
C             RH is converted to be wrt/ice for T<0 for this purpose.
C             -------------------------------------------------------
              RH = SGMET(J,3)
              TE = SGMET(J,2)
              PR = EXP(PRDBL(J))
              IF (TE.LT.0. .AND. TE.NE.BAD) THEN
                 Q = 0.01*RH*QSATW(TE+273.16,PR)
                 RH = Q/QSATI(TE+273.16,PR)*100.
                 ENDIF
C
              IF (J.LE.L .AND. RH.LT.95.) GOTO 130
              IF (J.GT.L .AND. RH.GE.95.) GOTO 130
132           CONTINUE
         CALL INSERT(NSIGT,ISIGT,L,SOK,23) 
130      CONTINUE 
C 
C     Find highest and lowest T,H 
C     --------------------------- 
      CALL MINMAX(SGMET(1,2),NT,XMIN,XMAX,IN,IX)
      CALL INSERT(NSIGT,ISIGT,IX,SOK,11) 
      CALL INSERT(NSIGT,ISIGT,IN,SOK,11) 
C 
      CALL MINMAX(SGMET(1,3),NT,XMIN,XMAX,IN,IX)
      CALL INSERT(NSIGT,ISIGT,IX,SOK,21) 
      CALL INSERT(NSIGT,ISIGT,IN,SOK,21) 
C 
C 
C     Now find all secondary significant levels for temperature 
C     --------------------------------------------------------- 
      CALL SECSIG(NSIGT,ISIGT,2,1)
C 
C     Find secondary levels for humidity
C     ----------------------------------
      CALL SECSIG(NSIGT,ISIGT,3,1)
C 
C 
C     Find primary significant levels for wind
C     ----------------------------------------
C 
C     First primary level is the surface. If no sfc, send 9999.
C     --------------------------------------------------------
C      IF (SNDDAT(IXEND,9).EQ.0. .OR. nsndtype .GE. 3) THEN 
         ISX = NT 
         IF (SGMET(NT,1).EQ.BAD) ISX = 9999
         CALL INSERT(NSIGW,ISIGW,ISX,SOK,30) 
C         ELSE 
C         CALL INSERT(NSIGW,ISIGW,9999,SOK,30) 
C         ENDIF
C 
C     Now find first and last good wind levels, also max
C     --------------------------------------------------
      CALL GAP(SGMET(1,5),NT,IBGN,IEND,NGAPS,BAD) 
      IB = IBGN(1)
      IE = IEND(NGAPS+1)
      CALL INSERT(NSIGW,ISIGW,IB,SOK,30) 
      CALL INSERT(NSIGW,ISIGW,IE,SOK,30) 
      CALL MINMAX(SGMET(1,5),NT,XMIN,XMAX,IN,IX)
      CALL INSERT(NSIGW,ISIGW,IX,SOK,31) 
C 
C     Add in secondary wind levels...First check on speed.
C     ----------------------------------------------------
      CALL SECSIG(NSIGW,ISIGW,5,1)
C 
C     Now add in those levels due to direction
C     ----------------------------------------
      CALL SECSIG(NSIGW,ISIGW,4,1)
C 
C
C     Check for doubtful temperatures
C     -------------------------------
C      ISTART = 1
C      NDT = 0
C      AG_10167 = .FALSE.
C200   DO 210 J=ISTART,IXEND
C         IF (ISNDFLG(J,2).EQ.4) THEN
C            AG_10167 = .TRUE.
C            NDT = NDT+1
C            IF (NDT.GT.NDTMAX) THEN
C               WRITE(LUT,'(" *** TOO MANY DOUBTFUL LAYERS ***")')
C               GOTO 230
C               ENDIF
C            PDTOP(NDT) = SNDDAT(J,2)
C            DO 220 L = J,IXEND
C              IF (ISNDFLG(L,2).NE.4) THEN
C                  PDBTM(NDT) = SNDDAT(L-1,2)
C                  ISTART = L+1
C                  IF (ISTART.GT.IXEND) GOTO 230
C                 GOTO 200
C                  ENDIF
C220            CONTINUE
C               PDBTM(NDT) = SNDDAT(IXEND,2)
C               GOTO 230
C            ENDIF
C210      CONTINUE
C
C
C     Check for doubtful geopotentials (AFRES ONLY).
C     Data marked doubtful if T doubtful or missing,
C     also if sfc was hydrostatic anchor.
C     ----------------------------------------------
C230   AG_10166 = .FALSE.
C      IF (.NOT. AFRES) GOTO 290
C      PZTOP = 1100.
C      PZBTM = 0.
C      DO 235 J = 1,IXEND
C         IF (SNDDAT(J,9) .NE. BAD) THEN
C            IF (SNDDAT(J,2).LT.PZTOP) PZTOP = SNDDAT(J,2)
C            IF (SNDDAT(J,2).GT.PZBTM) PZBTM = SNDDAT(J,2)
C            ENDIF
C235      CONTINUE
C      IF (PZTOP.EQ.1100. .OR. PZBTM.EQ.0.) GOTO 290
C
C      IF (PROC(7).EQ.2.) THEN
C         AG_10166 = .TRUE.
C         GOTO 290
C         ENDIF         
C      IF (AG_10167 .AND. PROC(7).EQ.1) THEN
C         AG_10166 = .TRUE.
C         IF (PZTOP.LT.PDTOP(1)) PZTOP = PDTOP(1)
C         GOTO 290
C         ENDIF
C      IF (PROC(7).EQ.1) THEN
C         DO 240 J = 1,IXEND
C            IF (SNDDAT(J,3).EQ.BAD) THEN
C               AG_10166 = .TRUE.
C               IF (PZTOP.LT.SNDDAT(J,2)) PZTOP = SNDDAT(J,2)
C               GOTO 290
C               ENDIF
C240         CONTINUE
C         GOTO 290
C         ENDIF
C
C290   CONTINUE
C
C 
C 
C     ------------------- 
C     Now code up message 
C     ------------------- 
C 
1     FORMAT(I1.1)
2     FORMAT(I2.2)
3     FORMAT(I3.3)
4     FORMAT(I4.4)
5     FORMAT(I5.5)
6     FORMAT(I6.6)
C 
C     Get day and hour of drop
C     ------------------------
C 
300   IDAY = IDY
      IF (MOD(IYR,4).EQ.0) MONTH(2) = 29
      IHOUR = NINT(HMSTS(TIML)/3600.) 
      IF (IHOUR.EQ.24) THEN 
         IHOUR = 0
         IDAY = IDAY+1
         IF (IDAY.GT.MONTH(IMO)) IDAY = 1 
         ENDIF
C 
C 
C     Start with message header info, sounding system info
C     ----------------------------------------------------
C 
      HEAD(1)(1:5) = 'UUBB '
	  HEAD(2)(1:5) = 'CGDX '
      WRITE(HEAD(3)(1:2),2) IDAY+50     
      WRITE(HEAD(3)(3:4),2) IHOUR 
      HEAD(3)(5:5) = '/'
      LAT = NINT(ABS(RLAT*10))
      LON = NINT(ABS(RLON*10))
      IF (RLAT.GE.0. .AND. RLON.GE.0.) IQ = 7 
      IF (RLAT.LT.0. .AND. RLON.GE.0.) IQ = 5 
      IF (RLAT.LT.0. .AND. RLON.LT.0.) IQ = 3 
      IF (RLAT.GE.0. .AND. RLON.LT.0.) IQ = 1 
      HEAD(4)(1:2) = '99' 
      WRITE(HEAD(4)(3:5),3) LAT 
      WRITE(HEAD(5)(1:1),1) IQ
      WRITE(HEAD(5)(2:5),4) LON 
      LATN = LAT/10
      LONN = LON/10
      CALL MARSDEN(LATN,LONN,IQ,MRSBLK)
      HEAD(6)(1:5) = MRSBLK
C
      SYSTEM(1)(1:5) = '31313'
      SYSTEM(2)(1:1) = '0'
      SYSTEM(2)(2:3) = '90'
	  SYSTEM(2)(4:5) = '08'
      SYSTEM(3)(1:1) = '8'
      WRITE(SYSTEM(3)(2:5),4) INT(TIML/100.)
C 
C 
C     Now start to code up the T/H sig levels 
C     --------------------------------------- 
C 
      IF (NSIGT.LE.1) THEN
         NSIGT = 0
         GOTO 355 
         ENDIF              
C 
      IP = -11
      DO 350 I = 1,NSIGT
C 
C        Grab values of P,T,TD
C        ---------------------
         J = ISIGT(I) 
         IF (J.EQ.9999) THEN 
              PR = BAD
              TE = BAD
              TD = BAD
              DD = BAD
              ELSE
              PR = EXP(PRDBL(J))  ! EXP(SGMET(J,1))
              TE = SGMET(J,2) 
              TD = DEWPT(TE,SGMET(J,3)) 
C 
C             Is this level identifying a data gap? 
C             ------------------------------------- 
              IF (TE.EQ.BAD .AND. TD.EQ.BAD .AND. J.NE.1) PR = BAD
C 
              DD = TE - TD
              IF (TE.EQ.BAD .OR. TD.EQ.BAD) DD = BAD
              ENDIF 

C
C	UUBB contains only levels below 100 MB
C
		if (PR .LT. 100.0) then
			levl(1,i) = NODAT
			levl(2,i) = NODAT
		else

C        Code indicator/pressure
C        -----------------------
         IP = IP + 11 
         IF (IP.GT.99) IP = 11  
         WRITE(LEVL(1,I)(1:2),2) IP 
C         IF (IP.EQ.0 .AND. AG_10191) PR = PR_10191
         IF (PR.EQ.BAD) THEN
              LEVL(1,I)(3:5) = '///'
              ELSE
 	      IPR = MOD(NINT(PR),1000)
              WRITE(LEVL(1,I)(3:5),3) IPR 
              ENDIF 
C 
C        Code temperature/dewpoint block
C        -------------------------------
         ITP = 2 * NINT(5.0*TE) 
         IF (ITP.GE.0) GOTO 351 
         IF (TE*10.0 .GE. FLOAT(ITP)) THEN
              ITP = ITP+1 
              ELSE
              ITP = ITP-1 
              ENDIF 
351      IF (TE.EQ.BAD) THEN
              LEVL(2,I)(1:3) = '///'
              ELSE
              WRITE(LEVL(2,I)(1:3),3) ABS(ITP)
              ENDIF 
C 
         IF (DD.EQ.BAD) THEN
              LEVL(2,I)(4:5) = '//' 
              ELSE
              IF (DD.LE.5.) THEN
                 IDD = NINT(DD*10.) 
                 ELSE 
                 IDD = 50 + NINT(DD)
                 ENDIF
              IF (IDD.GT.50 .AND. IDD.LE.55) IDD = 50 
              IF (IDD.GT.99) IDD = 99 
C             IF (IDD.GT.80) IDD = 80 
              WRITE(LEVL(2,I)(4:5),2) IDD 
              ENDIF 
C 
		endif

350      CONTINUE 
C 
C 
C     Do we have any significant wind levels? 
C     --------------------------------------- 
355   IF (NSIGW.LE.1) THEN
         NSIGW = 0
         K = NSIGT
         GOTO 500 
         ENDIF
C 
C 
C     Separate T/H levels from wind levels
C     ------------------------------------
      LEVL(1,NSIGT+1) = '21212' 
      LEVL(2,NSIGT+1) = NODAT 
C 
C 
C     Now start to code up the wind sig levels
C     ----------------------------------------
C 
      IP = -11
      DO 360 I = 1,NSIGW
C 
C        Grab values of P,WD,WS 
C        ---------------------- 
         K = NSIGT+1+I                 ! Index of LEVL
         J = ISIGW(I)                  ! Index of met arrays
         IF (J.EQ.9999) THEN 
              PR = BAD
              WD = BAD
              WS = BAD
              ELSE
              PR = EXP(PRDBL(J))  ! EXP(SGMET(J,1))
              WD = SGMET(J,4) 
              WS = SGMET(J,5) 
              ENDIF 

C
C	UUBB contains levels only below 100 mb
C

		if (PR .LT. 100.0) then
			levl(1,k) = NODAT
			levl(2,k) = NODAT
		else


C        Code indicator/pressure
C        -----------------------
         IP = IP + 11 
         IF (IP.GT.99) IP = 11  
         WRITE(LEVL(1,K)(1:2),2) IP 
C         IF (IP.EQ.0 .AND. AG_10191) PR = PR_10191
C         IF (PR.EQ.BAD .OR. WD.EQ.BAD .OR. WS.EQ.BAD) THEN
         IF (PR.EQ.BAD) THEN
              LEVL(1,K)(3:5) = '///'
              ELSE
	      IPR = MOD(NINT(PR),1000)
              WRITE(LEVL(1,K)(3:5),3) IPR 
              ENDIF 
C 
C        Code wind block
C        ---------------
         CALL WNDCODE(WD,WS,IWD,IWS)
         WRITE(LEVL(2,K)(1:2),2) IWD
         IF (WD.EQ.BAD) LEVL(2,K)(1:2) = '//' 
         WRITE(LEVL(2,K)(3:5),3) IWS
         IF (WS.EQ.BAD) LEVL(2,K)(3:5) = '///'
C
 
		endif

360      CONTINUE 
C 
C 
C     ----------------------------------------------
C     Place individual blocks into final code array.
C     First do header.
C     ----------------------------------------------
C 
500   NLEVELS = K 
      ICNT = 1
      LINE = 1
C 
      DO 530 I = 1,6
         CODE(LINE)(ICNT:ICNT+4) = HEAD(I)
         CODE(LINE)(ICNT+5:ICNT+5) = BLANK
         CALL INCR(ICNT,LINE,LNWMO) 
530      CONTINUE 
C
      IF (NLEVELS.EQ.0) GOTO 600
C 
C 
C     Now do level data 
C     ----------------- 
      DO 545 IK = 1, NLEVELS
         DO 540  IL = 1,2 
C 
C             If start of winds, skip to next line
C             ------------------------------------
              IF (IL.EQ.1 .AND. LEVL(IL,IK).EQ.'21212'
     *            .AND. ICNT.NE.1) THEN
                 DO 546 I = ICNT,LNWMO
546                 CODE(LINE)(I:I) = BLANK 
                 LINE = LINE+1
                 ICNT = 1 
                 ENDIF
C 
              IF (LEVL(IL,IK).EQ.NODAT) GO TO 540 
              CODE(LINE)(ICNT:ICNT+4) = LEVL(IL,IK) 
              CODE(LINE)(ICNT+5:ICNT+5) = BLANK 
              CALL INCR(ICNT,LINE,LNWMO)
540           CONTINUE
545      CONTINUE 
C 
C
C     Add in system status groups (31313).  Start on new line.
C     --------------------------------------------------------
600   IF (ICNT.EQ.1) THEN 
         LINE = LINE - 1
         ELSE 
         DO 610 I = ICNT,LNWMO
            CODE(LINE)(I:I) = BLANK 
610         CONTINUE
         ENDIF
C
      LINE = LINE+1
      ICNT = 1 
      DO 620 J = 1,3
         CODE(LINE)(ICNT:ICNT+4) = SYSTEM(J)(1:5)
         CODE(LINE)(ICNT+5:ICNT+5) = BLANK
         CALL INCR(ICNT,LINE,LNWMO)
620      CONTINUE   
C
C 
C     Add in additional data groups (51515).  Start on new line.
C     ----------------------------------------------------------
C      AG_51515 = .FALSE.
C      IF (AG_10191    .OR. 
C     *    AG_10190(1) .OR. 
C     *    AG_10190(2) .OR.
C     *    AG_10166    .OR.
C     *    AG_10167)    AG_51515 = .TRUE.
C     IF (.NOT. AG_51515) GOTO 700
C
C      IF (ICNT.EQ.1) THEN 
C         LINE = LINE - 1
C         ELSE 
C         DO 630 I = ICNT,LNWMO
C            CODE(LINE)(I:I) = BLANK 
C630         CONTINUE
C         ENDIF
C
C      LINE = LINE+1
C      ICNT = 1 
C      CODE(LINE)(ICNT:ICNT+4) = '51515'
C      CODE(LINE)(ICNT+5:ICNT+5) = BLANK
C      CALL INCR(ICNT,LINE,LNWMO,LAST_LINE,LAST_ICNT)
C
C      IF (AG_10166) THEN
C         CODE(LINE)(ICNT:ICNT+4) = '10166'
C         CODE(LINE)(ICNT+5:ICNT+5) = BLANK
C         CALL INCR(ICNT,LINE,LNWMO,LAST_LINE,LAST_ICNT)
C         IP1 = NINT(PZBTM/10.)
C         IF (IP1.GE.100) IP1 = IP1-100
C         IP2 = NINT(PZTOP/10.)
C         IF (IP2.GE.100) IP2 = IP2-100
C         CODE(LINE)(ICNT:ICNT) = '0'
C         WRITE(CODE(LINE)(ICNT+1:ICNT+2),2) IP1
C         IF (PZBTM.EQ.BAD) CODE(LINE)(ICNT+1:ICNT+2) = '//'
C         WRITE(CODE(LINE)(ICNT+3:ICNT+4),2) IP2
C         IF (PZTOP.EQ.BAD) CODE(LINE)(ICNT+3:ICNT+4) = '//'
C         CODE(LINE)(ICNT+5:ICNT+5) = BLANK
C         CALL INCR(ICNT,LINE,LNWMO,LAST_LINE,LAST_ICNT)
C         ENDIF
C
C      IF (AG_10167) THEN
C         CODE(LINE)(ICNT:ICNT+4) = '10167'
C         CODE(LINE)(ICNT+5:ICNT+5) = BLANK
C         CALL INCR(ICNT,LINE,LNWMO,LAST_LINE,LAST_ICNT)
C         DO 640 I = 1,NDT
C            IP1 = NINT(PDBTM(I)/10.)
C            IF (IP1.GE.100) IP1 = IP1-100
C            IP2 = NINT(PDTOP(I)/10.)
C            IF (IP2.GE.100) IP2 = IP2-100
C            CODE(LINE)(ICNT:ICNT) = '0'
C            WRITE(CODE(LINE)(ICNT+1:ICNT+2),2) IP1
C            IF (PDBTM(I).EQ.BAD) CODE(LINE)(ICNT+1:ICNT+2) = '//'
C            WRITE(CODE(LINE)(ICNT+3:ICNT+4),2) IP2
C            IF (PDTOP(I).EQ.BAD) CODE(LINE)(ICNT+3:ICNT+4) = '//'
C            CODE(LINE)(ICNT+5:ICNT+5) = BLANK
C            CALL INCR(ICNT,LINE,LNWMO,LAST_LINE,LAST_ICNT)
C640         CONTINUE
C         ENDIF
C
C      DO 670 I = 1,2
C         IF (AG_10190(I)) THEN
C            CODE(LINE)(ICNT:ICNT+4) = '10190'
C            CODE(LINE)(ICNT+5:ICNT+5) = BLANK
C            CALL INCR(ICNT,LINE,LNWMO,LAST_LINE,LAST_ICNT)
C            IP = PR_10190(I)/10. 
C            IF (IP.EQ.100) IP = 0 
C            IF (GA_10190(I).LT.0.) GA_10190(I) = 500.+ABS(GA_10190(I))
C            IHT = NINT(GA_10190(I))
C            IF (PR_10190(I).LE.500.) IHT = NINT(GA_10190(I)/10) 
C            IHT = MOD(IHT,1000) 
C            WRITE(CODE(LINE)(ICNT:ICNT+1),2) IP
C            IF (PR_10190(I).EQ.BAD) CODE(LINE)(ICNT:ICNT+1) = '//'
C            WRITE(CODE(LINE)(ICNT+2:ICNT+4),3) IHT 
C            IF (GA_10190(I).EQ.BAD) CODE(LINE)(ICNT+2:ICNT+4) = '///' 
C            CODE(LINE)(ICNT+5:ICNT+5) = BLANK
C            CALL INCR(ICNT,LINE,LNWMO,LAST_LINE,LAST_ICNT)
C            ENDIF
C670      CONTINUE
C
C      IF (AG_10191) THEN
C         CODE(LINE)(ICNT:ICNT+4) = '10191'
C         CODE(LINE)(ICNT+5:ICNT+5) = BLANK
C         CALL INCR(ICNT,LINE,LNWMO,LAST_LINE,LAST_ICNT)
C         ENDIF
C
C700	CONTINUE
C	add trailing =
c	CODE(LAST_LINE)(LAST_ICNT+5:LAST_ICNT+5) = '='
C
C
C     -------------------------------
C     All through...
C     Blank out remaining code groups 
C     ------------------------------- 
C 
	   IF (ICNT.EQ.1) THEN 
         LINE = LINE - 1
         ELSE 
         DO 750 IK = ICNT,LNWMO 
750           CODE(LINE)(IK:IK) = BLANK 
         ENDIF
C 
C 
C     Write out data to output file 
C     ----------------------------- 
C 
      CALL PRINTSIG
C
      WRITE(LUT,*)
      DO 760 IJ = 1,LINE
		 nMsgLines_thr = nMsgLines_thr + 1
		 temp_msg_thr(nMsgLines_thr) = CODE(IJ)
         WRITE(LUT,'(1X,A)') CODE(IJ) 
760      CONTINUE 
C 
c      DO 765 IJ = 1,LINE
c         WRITE(LUFO,'(A)',ERR=900) CODE(IJ) 
765      CONTINUE 
C 
      RETURN
C 
C 
900   WRITE(LUT,'(/," *** ERROR WRITING UUBB MESSAGE ***",/)')
      RETURN
      END 
C 
C 
C 
C     ----------------------------------------------------- 
      SUBROUTINE LOADSGMET
C     ----------------------------------------------------- 
C 
      PARAMETER (MXRC = 9000, ninvar = 12)
C
      REAL*8 PRDBL
C
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /SIGLEVEL/ SGMET(MXRC+1,7), NT, PRDBL(MXRC+1)
C
      DATA BAD/-999./ 
C 
      IF (PRL.GT.0.) THEN 
         SGMET(1,1) = ALOG(PRL) 
         PRDBL(1) = LOG(DBLE(PRL))
         ELSE 
         SGMET(1,1) = BAD 
         PRDBL(1) = BAD
         ENDIF
      SGMET(1,2) = TEL
      SGMET(1,3) = RHL
      SGMET(1,4) = WDL
      SGMET(1,5) = WSL
      SGMET(1,6) = UCMP(WDL,WSL) 
      SGMET(1,7) = VCMP(WDL,WSL) 
C
C
C     Throw out RH if T is bad
C     ------------------------
      DO 50 L = 1, IXEND
         IF (SNDDAT(L,2).GT.0.) THEN
              SGMET(L+1,1) = LOG(SNDDAT(L,2))
              PRDBL(L+1) = LOG(DBLE(SNDDAT(L,2)))
              ELSE
              SGMET(L+1,1) = BAD
              PRDBL(L+1) = BAD
              ENDIF 
         SGMET(L+1,2) = SNDDAT(L,3) 
         SGMET(L+1,3) = SNDDAT(L,4)
         IF (SGMET(L+1,2).EQ.BAD) SGMET(L+1,3) = BAD
         SGMET(L+1,4) = SNDDAT(L,6) 
         SGMET(L+1,5) = SNDDAT(L,7) 
         SGMET(L+1,6) = UCMP(SGMET(L+1,4),SGMET(L+1,5))
         SGMET(L+1,7) = VCMP(SGMET(L+1,4),SGMET(L+1,5))
50       CONTINUE 
C
C
C     Substitute 10 m wind for surface value, and toss any other
C     winds below 10 m elevation.
C     ----------------------------------------------------------
      NT = IXEND+1
      IF (SNDDAT(IXEND,9).NE.0.) GOTO 100
      CALL WINDZINT(10.,WD10,WS10,WQ10,M10)
      IF (WD10.GE.0) THEN
         SGMET(NT,4) = WD10
         SGMET(NT,5) = WS10
         SGMET(NT,6) = UCMP(SGMET(NT,4),SGMET(NT,5))
         SGMET(NT,7) = VCMP(SGMET(NT,4),SGMET(NT,5))
         DO 60 L = IXEND-1,IXEND-10,-1
            IF (SNDDAT(L,9).NE.BAD .AND. SNDDAT(L,9).LT.10.) THEN
               SGMET(L+1,4) = BAD
               SGMET(L+1,5) = BAD
               SGMET(L+1,6) = BAD
               SGMET(L+1,7) = BAD
               ENDIF
60          CONTINUE
         ENDIF
C
C
100   RETURN
      END
C
C
C
C     ----------------------------------------------------- 
      SUBROUTINE CHECKINV(IXINV,NTROP,IXTRP,INVOK)
C     ----------------------------------------------------- 
C 
      PARAMETER (MXRC = 9000, ninvar = 12)
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      DATA BAD/-999./ 
C 
      INVOK = 0
      IF (NTROP.EQ.0) THEN
         PRTROP = 1000.
         ELSE
         PRTROP = SNDDAT(IXTRP,2)
         IF (PRTROP.EQ.BAD) PRTROP = 1000.
         ENDIF
      PLIM = AMIN1(PRTROP,300.)
C
      IF (SNDDAT(IXINV,2).GT.PLIM) INVOK = 1
      RETURN
      END
C
C
C
C     ----------------------------------------
      SUBROUTINE MARSDEN(LAT,LON,IQ,MRSBLK)
C
C     CALCULATES MARSDEN SQ BLOCK
C     ----------------------------------------
C
      CHARACTER*5 MRSBLK
C
      MSQ = -1
C     Determine latitude, longitude bands
C
      LATBAND = INT(LAT/10)+1
      LONBAND = INT(LON/10)+1
C
C     Special fix for dateline, pole
      IF (LONBAND.EQ.19) LONBAND=18
      IF (LATBAND.EQ.10) LATBAND=9
C
C     NW Quadrant
      IF (IQ.EQ.7) MSQ = 36*(LATBAND-1) + LONBAND
C
C     NE Quadrant
      IF (IQ.EQ.1) MSQ = 36*LATBAND - (LONBAND-1)
C
C     SW Quadrant
      IF (IQ.EQ.5) MSQ = 299 + 36*(LATBAND-1) + LONBAND 
C
C     SE Quadrant
      IF (IQ.EQ.3) MSQ = 299 + (36*LATBAND) - (LONBAND-1)
C
C     Special fix for N of 80N
      IF (LAT.GE.80 .AND. IQ.EQ.7) MSQ = 900 + LONBAND
      IF (LAT.GE.80 .AND. IQ.EQ.1) MSQ = 937 - LONBAND
C
C     Now get units digit of lat and lon

      LATBAND = INT(LAT/10)+1
      LONBAND = INT(LON/10)+1
      LATI = INT(LAT)
      LONI = INT(LON)
      LATU = LATI - 10*(LATBAND-1)
      LONU = LONI - 10*(LONBAND-1)
C
C     Now code up block          
C
1     FORMAT(I1.1)
3     FORMAT(I3.3)
C
      IF (MSQ.EQ.-1) THEN
          MRSBLK = '/////'
          ELSE
          WRITE(MRSBLK(1:3),3) MSQ
          WRITE(MRSBLK(4:4),1) LATU
          WRITE(MRSBLK(5:5),1) LONU
          ENDIF
      RETURN
      END
C
C
C
C     ----------------------------------------
      SUBROUTINE INCR(ICNT,LINE,LNWMO)
C 
C     INCREMENT LINE AND COUNT FOR CODE BUFFER
C     ----------------------------------------
C 
      ICNT = ICNT + 6 
      IF(ICNT.GE.LNWMO)THEN 
        ICNT = 1
        LINE = LINE + 1 
        ENDIF 
      RETURN
      END 
C 
C 						  C 
C 
C     --------------------------------------------------
      SUBROUTINE WNDCODE(WDIR,WSPD,IWD,IWS) 
C 
C     CODE UP WIND DIRECTION AND SPEED INTO INTEGER FORM
C     CHECK FOR BAD DATA
C     --------------------------------------------------
C 
      BAD = -999. 
      IF(WSPD.EQ.BAD.OR.WDIR.EQ.BAD) THEN 
        IWS = 0 
        IWD = 00
        RETURN
        ENDIF 
C 
C 
C     FIND WIND DIRECTION TO NEAREST 10 DEGREES 
C     ----------------------------------------- 
      IWD = INT(WDIR/10.) 
C 
C     FIND REMAINDER OF WIND DIRECTION
C     --------------------------------
      REM = WDIR - (IWD*10.)
      IF(REM.LT.2.5) REM = 0. 
      IF(REM.GE.2.5.AND.REM.LE.7.5) REM = 5.
      IF(REM.GT.7.5)THEN
        REM = 0.
        IWD = IWD + 1 
        ENDIF 
C 
C     FIND WIND SPEED IN KNOTS PLUS DIRECTION REMAINDER 
C     ------------------------------------------------- 
      IWS = NINT(WSPD*1.94) 
      IF (IWS.EQ.0) THEN
        IWD = 0 
        RETURN
        ENDIF 
      IWS = IWS + (REM*100.)
C 
C     CORRECT FOR NORTH WIND
C     ----------------------
      IF (IWD.EQ.0 .AND. IWS.LT.500) IWD = 36 
C 
C     CORRECT FOR CALM WIND 
C     --------------------- 
      IF(IWS .EQ. 0) IWD = 00 
C 
      RETURN
      END
C 
C 
C 
C     ----------------------------------------------------- 
      SUBROUTINE SECSIG(NSIG,ISIG,IVAR,MODE_WMO) 
C 
C     Subroutine identifies secondary significant levels. 
C     These are levels that deviate from staight lines by 
C     some specified amount (DEVSIG).  See FMH #4, pg B2-2 for
C     requirements for T and H.  Routine is for any scalar. 
C 
C     NSIG  = current number of sig levels
C     ISIG  = array of sig level data 
C     IVAR  = variable being checked. Corresponds to
C             elements in SGMET data arrays:
C 
C             2 = TEMP
C             3 = RH
C             4 = WD
C             5 = WS
C
C     MODE_WMO  = 1 for NOAA, 2 for AF conventions.
C     ----------------------------------------------------- 
C 
C 
      PARAMETER (MXRC = 9000)
      DIMENSION ISIG(*) 
      DIMENSION DEVSIG(5) 
      REAL*8 PRDBL
      LOGICAL SOK,VECTOR,DEBUG
C 
      COMMON /SIGLEVEL/ SGMET(MXRC+1,7), NT, PRDBL(MXRC+1)
C 
      DATA DEVSIG/-999.,1.0,10.,10.,5.15/ 
      DATA BAD/-999./ 
      DATA DEBUG/.FALSE./
C 
C 
C     Can't have any secondary sig levels if only one primary 
C     ------------------------------------------------------- 
      IF (NSIG.LE.1) RETURN     
C 
C     Are we dealing with a vector computation (WD)?
C     ----------------------------------------------
      VECTOR = .FALSE.
      IF (IVAR.EQ.4) VECTOR = .TRUE.
C 
C     Loop through all primary sig level intervals
C     --------------------------------------------
C 
      L = 0 
100   L = L + 1 
      IF (L.EQ.NSIG) RETURN                ! Done 
C 
      IF (ISIG(L).EQ.9999) GOTO 100         ! Sfc is missing 
C 
      KBTM = ISIG(L)                       ! Bottom primary index 
      KTOP = ISIG(L+1)                     ! Top primary index
      PBTM = SGMET(KBTM,1)                 ! Top of primary interval
      PTOP = SGMET(KTOP,1)                 ! Bottom of primary interval 
      PTOPT = PTOP                         ! Temporary top
      KTOPT = KTOP                         ! Temporary top index
      XTOPT = SGMET(KTOP,IVAR)             ! X val at top 
      XBTM = SGMET(KBTM,IVAR)              ! X val at btm 
C 
      IF (VECTOR) THEN
         XTOPT = UCMP(SGMET(KTOP,4),1.0) 
         YTOPT = VCMP(SGMET(KTOP,4),1.0) 
         XBTM  = UCMP(SGMET(KBTM,4),1.0) 
         YBTM  = VCMP(SGMET(KBTM,4),1.0) 
         ENDIF
C 
C 
C     If data are missing at either end, there can be no
C     secondary sig levels in this primary interval.
C     --------------------------------------------------
      IF (XTOPT.EQ.BAD .OR. XBTM.EQ.BAD) GOTO 100 
      IF (VECTOR .AND. (YTOPT.EQ.BAD .OR. YBTM.EQ.BAD)) GOTO 100
C 
C     Search all met data that fall between PBTM and PTOPT, and 
C     find the point which has the maximum deviation
C     --------------------------------------------------------- 
110   DEVMAX = 0. 
      IF (DEBUG) WRITE(1,*) IVAR,PBTM,XBTM,PTOPT,XTOPT
      DO 120 K = 1,NT       
         IF (SGMET(K,1).GE.PBTM) GOTO 120 
         IF (SGMET(K,1).LE.PTOPT) GOTO 120
         IF (SGMET(K,IVAR).EQ.BAD) GOTO 120 
C 
C        Special check for direction when speed < 10 kts
C        -----------------------------------------------
         IF (VECTOR .AND. SGMET(K,5).LT.5.15) GOTO 120
C 
         XLINE = POLATE2(PBTM,PTOPT,XBTM,XTOPT,SGMET(K,1),BAD)
         DEV = ABS(XLINE-SGMET(K,IVAR)) 
C 
         IF (VECTOR) THEN 
              YLINE = POLATE2(PBTM,PTOPT,YBTM,YTOPT,SGMET(K,1),BAD) 
              WDLINE = WDCOMP(XLINE,YLINE)
              WDPOINT = SGMET(K,IVAR) 
              DEV = ABS(WDLINE-WDPOINT) 
              IF (DEV.GT.180.) DEV = 360.-DEV 
              ENDIF 
C 
         IF (DEV.GT.DEVMAX) THEN
              DEVMAX = DEV
              KX = K
              ENDIF 
120      CONTINUE 
C 
C     Does the max deviation exceed the sig level criterion?
C     If so, this point is our new temporary top, and we go back
C     to search below this new top. 
C     ----------------------------------------------------------
      DEVLIM = DEVSIG(IVAR)
C
C     Special loosening of WS threshhold below 850 mb to
C     get more levels (NOAA only).
C     ------------------------------------------------
      IF (IVAR.EQ.5 .AND. MODE_WMO.EQ.1) THEN
         PRCUR = BAD
         IF (PRDBL(KX).NE.BAD) PRCUR = EXP(PRDBL(KX))
         IF (PRCUR.GT.850.) DEVLIM = DEVLIM/2.
         ENDIF
C
      IF (DEBUG) WRITE(1,*) DEVMAX,DEVLIM
      IF (DEVMAX.GE.DEVLIM) THEN
         KTOPT = KX 
         PTOPT = SGMET(KTOPT,1)            ! New temporary top
         XTOPT = SGMET(KTOPT,IVAR)         ! Value at new ttop. 
         IF (VECTOR) THEN 
              XTOPT = UCMP(SGMET(KTOPT,4),1.0) 
              YTOPT = VCMP(SGMET(KTOPT,4),1.0) 
              ENDIF 
         GOTO 110 
         ENDIF
C 
C     We did not find a sig level within the current interval.
C     If the temporary top = the upper primary sig level, 
C     then there are no secondary sig levels in this primary
C     interval: go back to main loop and look at next primary 
C     interval. 
C     --------------------------------------------------------- 
      IF (PTOP.EQ.PTOPT) GOTO 100 
C 
C     The temporary top was not equal to the upper primary
C     sig level, so we make the temporary top a new sig level,
C     and the bottom of the next primary level
C     ------------------------------------------------------- 
      IF (IVAR.EQ.2) ISTYP = 12
      IF (IVAR.EQ.3) ISTYP = 22
      IF (IVAR.EQ.4) ISTYP = 33
      IF (IVAR.EQ.5) ISTYP = 32
      CALL INSERT(NSIG,ISIG,KTOPT,SOK,ISTYP)
      GOTO 100
C 
      END 
C 
C 
C 
C     ----------------------------------------------------- 
      SUBROUTINE SIG20(NSIG,ISIG) 
C 
C     Subroutine identifies significant levels within 20 mb
C     of the surface, by taking the point of maximum 
C     deviation from linearity for temperature.
C
C     Stripped down version of SECSIG.
C 
C     NSIG  = current number of sig levels
C     ISIG  = array of sig level data 
C     ----------------------------------------------------- 
C 
C 
      PARAMETER (MXRC = 9000)
      DIMENSION ISIG(*) 
      REAL*8 PRDBL
      LOGICAL SOK
C 
      COMMON /SIGLEVEL/ SGMET(MXRC+1,7), NT, PRDBL(MXRC+1)
C 
      DATA BAD/-999./ 
C 
C 
C     Can't have any secondary sig levels if only one primary 
C     ------------------------------------------------------- 
      IF (NSIG.LE.1) RETURN     
C 
C
C     Set up search domain (within 20 mb of surface).
C     -----------------------------------------------
      IF (ISIG(1).EQ.9999) RETURN          ! Sfc is missing 
      KBTM = ISIG(1)                       ! Bottom index
      IF (PRDBL(KBTM).EQ.BAD) THEN
         RETURN
         ELSE
         PRSFC = EXP(PRDBL(KBTM))
         ENDIF
C
      DO 50 L = KBTM-1,1,-1
         IF (PRDBL(L).EQ.BAD) RETURN
         PRCUR = EXP(PRDBL(L))
         IF (PRSFC-PRCUR .GT. 20.) GOTO 55
50       CONTINUE
C
55    IF (PRSFC-PRCUR.LT.20.) RETURN
      KTOP = L+1                           ! Top primary index
      PBTM = SGMET(KBTM,1)                 ! Bottom of primary interval
      PTOP = SGMET(KTOP,1)                 ! Top of primary interval 
      PTOPT = PTOP                         ! Temporary top
      KTOPT = KTOP                         ! Temporary top index
      XTOPT = SGMET(KTOP,2)                ! T val at top 
      XBTM = SGMET(KBTM,2)                 ! T val at btm 
C 
C 
C     If data are missing at either end, there can be no
C     secondary sig levels in this primary interval.
C     --------------------------------------------------
      IF (XTOPT.EQ.BAD .OR. XBTM.EQ.BAD) RETURN
C 
C
C     Search all met data that fall between PBTM and PTOPT, and 
C     find the point which has the maximum deviation
C     --------------------------------------------------------- 
      DEVMAX = 0. 
      DO 120 K = 1,NT       
         IF (SGMET(K,1).GE.PBTM) GOTO 120 
         IF (SGMET(K,1).LE.PTOPT) GOTO 120
         IF (SGMET(K,2).EQ.BAD) GOTO 120 
         XLINE = POLATE2(PBTM,PTOPT,XBTM,XTOPT,SGMET(K,1),BAD)
         DEV = ABS(XLINE-SGMET(K,2)) 
         IF (DEV.GT.DEVMAX) THEN
              DEVMAX = DEV
              KX = K
              ENDIF 
120      CONTINUE 
C
C 
C     Accept the point of max deviation
C     ---------------------------------
      ISTYP = 19
      CALL INSERT(NSIG,ISIG,KX,SOK,ISTYP)
      RETURN
      END 
C 
C 
C 
C     --------------------------------------------------- 
      SUBROUTINE INSERT(NSIG,ISIG,ISX,SOK,ISTYP)
C 
C     Subroutine to insert sig level at index = ISX into
C     previously defined levels in ISIG.  Keeps lowest
C     level (surface) as the first element in ISIG.  Total
C     number of levels prior to insertion is received as
C     N, if level is accepted, subroutine increments N, 
C     and sets SOK to .TRUE.
C     --------------------------------------------------- 
C 
      PARAMETER (MXRC = 9000)
      DIMENSION ISIG(NSIG+1)
      REAL*8 PRDBL
      LOGICAL SOK 
C 
      COMMON /SIGLEVEL/ SGMET(MXRC+1,7), NT, PRDBL(MXRC+1)
      COMMON /SIGTYPES/ NSIGL,ISIGL(MXRC+1),ISIGTYPE(MXRC+1)
C
      DATA BAD/-999./ 
C 
C 
      SOK = .FALSE. 
C 
C     If this is the first (sfc) level, skip normal checks
C     ----------------------------------------------------
      IF (NSIG.EQ.0) THEN 
         L = 1
         GOTO 220 
         ENDIF
C 
C     Do not accept invalid indices 
C     ----------------------------- 
      IF (ISX.LE.0. .OR. ISX.GT.NT) RETURN
C 
C     Do not accept bad pressures
C     ---------------------------
      PRNOW = SGMET(ISX,1)
      IF (PRNOW.EQ.BAD) RETURN
C 
C 
C     Check to see where the new level fits into array
C     Check on pressure, not on index number..........
C     ------------------------------------------------
      DO 100 L = 1,NSIG 
C 
C        If we already have this level, return
C        -------------------------------------
         IF (ISX.EQ.ISIG(L)) RETURN 
C 
         IF (ISIG(L).EQ.9999) THEN 
              PR = 99999. 
              ELSE
              PR = SGMET(ISIG(L),1) 
              ENDIF 
         LX = L
         IF (PRNOW.GT.PR) GOTO 200
C 
100      CONTINUE 
C 
C 
C     Do not allow a level if there is already one at
C     NINT(PR), so we don't have two levels with the
C     same CODED pressure.  Exception is for temp gaps
C     (type 15), which have P coded as /// so these
C     are OK.  (Added 3/7/97 JLF.)
C     ------------------------------------------------
200   IF (ISIG(LX).NE.9999 .AND. ISTYP.NE.15) THEN
         IPNOW = NINT(EXP(PRDBL(ISX)))
         DO 205 K = 1,NSIG
            IPOLD = NINT(EXP(PRDBL(ISIG(K))))
            IF (IPNOW.EQ.IPOLD) RETURN
205         CONTINUE
         ENDIF
C
C
C     New level fits in at index L.  Bump remaining levels down 
C     and insert the new level. 
C     ---------------------------------------------------------         
      DO 210 K = NSIG,L,-1
         ISIG(K+1) = ISIG(K)
210      CONTINUE 
C 
C
220   NSIG = NSIG + 1 
      SOK = .TRUE.
      ISIG(L) = ISX 
C
C
C     Keep track of sig level types to plot
C     ------------------------------------- 
      IF (ISX.LE.0 .OR. ISX.GT.NT) RETURN
      NSIGL = NSIGL+1
      ISIGL(NSIGL) = ISX
      ISIGTYPE(NSIGL) = ISTYP
C
      RETURN
      END 
	
     
C 
C 
C 
C     ------------------------------------------------------- 
      SUBROUTINE MINMAX(X,NT,XMIN,XMAX,IN,IX) 
C 
C     Finds minimum and maximum values of array X.  Returns 
C     min and max values, and element #s for these values.
C     If all data are bad, -1 is returned for element #.
C     ------------------------------------------------------- 
C     
      DIMENSION X(NT) 
      DATA BAD/-999./ 
C 
      IX = -1 
      IN = -1 
      XMIN = 999999.
      XMAX = -999999. 
C 
      DO 100 L = 1,NT 
         IF (X(L).EQ.BAD) GOTO 100
         IF (X(L).GT.XMAX) THEN 
              XMAX = X(L) 
              IX = L
              ENDIF 
         IF (X(L).LT.XMIN) THEN 
              XMIN = X(L) 
              IN = L
              ENDIF 
100      CONTINUE 
C 
      IF (IX.EQ.-1) XMAX = BAD
      IF (IN.EQ.-1) XMIN = BAD
C 
      RETURN
      END 
C 
C 
C 
C     --------------------------------------------------- 
      SUBROUTINE PRINTSIG
C 
C     Displays sig level data 
C     --------------------------------------------------- 
C 
      PARAMETER (MXRC = 9000)
C 
      CHARACTER*3 TYPE
      REAL*8 PRDBL
C 
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SIGLEVEL/ SGMET(MXRC+1,7), NT, PRDBL(MXRC+1)
      COMMON /SIGTYPES/ NSIGL,ISIGL(MXRC+1),ISIGTYPE(MXRC+1)
C
C 
      WRITE(LUT,'(//," Significant levels: ")') 
      WRITE(LUT,'(/,
     * "   #   TYPE      PR       TE       RH       WD       WS ")') 
      WRITE(LUT,'(  
     * " -------------------------------------------------------")')
C 
      DO 100 L = 1, NSIGL
         I = ISIGL(L) 
         IF (I.EQ.9999) THEN 
              PR = -999. 
              TE = -999.
              RH = -999.
              WD = -999.
              WS = -999.
              TYPE = '   '
              ELSE
              PR = EXP(SGMET(I,1)) 
              TE = SGMET(I,2) 
              RH = SGMET(I,3) 
              WD = SGMET(I,4) 
              WS = SGMET(I,5) 
              TYPE = '   '
              IF (ISIGTYPE(L).EQ.10) TYPE = 'T:E'
              IF (ISIGTYPE(L).EQ.11) TYPE = 'T:X'
              IF (ISIGTYPE(L).EQ.12) TYPE = 'T:S'
              IF (ISIGTYPE(L).EQ.13) TYPE = 'T:I'
              IF (ISIGTYPE(L).EQ.14) TYPE = 'T:I'
              IF (ISIGTYPE(L).EQ.15) TYPE = 'T:G'
              IF (ISIGTYPE(L).EQ.16) TYPE = 'T:L'
              IF (ISIGTYPE(L).EQ.17) TYPE = 'T:L'
              IF (ISIGTYPE(L).EQ.18) TYPE = 'T:T'
              IF (ISIGTYPE(L).EQ.19) TYPE = 'T:2'
              IF (ISIGTYPE(L).EQ.20) TYPE = 'H:E'
              IF (ISIGTYPE(L).EQ.21) TYPE = 'H:X'
              IF (ISIGTYPE(L).EQ.22) TYPE = 'H:S'
              IF (ISIGTYPE(L).EQ.23) TYPE = 'H:C'
              IF (ISIGTYPE(L).EQ.30) TYPE = 'W:E'
              IF (ISIGTYPE(L).EQ.31) TYPE = 'WSX'
              IF (ISIGTYPE(L).EQ.32) TYPE = 'WSS'
              IF (ISIGTYPE(L).EQ.33) TYPE = 'WDS'
              ENDIF 
C 
         WRITE(LUT,990) L,TYPE,PR,TE,RH,WD,WS   
100      CONTINUE 
C 
990   FORMAT(1X,I3,4X,A3,3F9.1,F9.0,F9.1) 
C 
      WRITE(LUT,'(
     *  " -------------------------------------------------------")')
      RETURN
      END 
C
C
C
C     ---------------------------------------------------------
      SUBROUTINE EXTRAP_TEMP(P1,PX,TX,BAD)
C 
C     Extrapolate from pressure P1 to get the temperature (TX)
C     at pressure PX.  Routine follows WMO procedure of going
C     back P1-PX mb to do the extrapolation.
C     ---------------------------------------------------------
C 
C 
      PARAMETER (MXRC = 9000, ninvar = 12)
C 
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
C 
C
C     Check incoming data
C     -------------------
      IF (P1.EQ.BAD .OR. PX.EQ.BAD) GOTO 900
      DO 100 I = 1,IXEND
         IF (SNDDAT(I,2).EQ.P1) THEN
            IP1 = I
            GOTO 150
            ENDIF
100      CONTINUE
      GOTO 900
C
C
C     Find out how many array elements to go back
C     -------------------------------------------
150   PDELTA = ABS(PX - P1)
      NDELTA = NELM(PDELTA,FALLRATE,DATARATE)
      IF (NDELTA.LT.1) NDELTA = 1
      IF (PX.LT.P1) NDELTA = -NDELTA
      IPE = IP1 - NDELTA
      IF (IPE.LT.1) IPE = 1
      IF (IPE.GT.IXEND) IPE = IXEND
C
C
C     Calculate lapse rate (deg/mb)
C     -----------------------------
      IF (SNDDAT(IPE,2).EQ.BAD .OR. SNDDAT(IPE,3).EQ.BAD) GOTO 900
      IF (SNDDAT(IP1,2).EQ.BAD .OR. SNDDAT(IP1,3).EQ.BAD) GOTO 900
      TDELTA = SNDDAT(IP1,3)-SNDDAT(IPE,3)
      PDELTA = SNDDAT(IP1,2)-SNDDAT(IPE,2)
      IF (PDELTA.EQ.0.) GOTO 900
      TLAPS = TDELTA/PDELTA
C
C
C     Extrapolate temperature
C     -----------------------
      TX = SNDDAT(IP1,3) + (PX-P1)*TLAPS
      RETURN
C
C
900   TX = BAD
      RETURN
      END	
	
C
C
C
C     -------------------------------------------------
      FUNCTION RHCHK(RH,RHMIN)
C
C     Routine sets lower limit to allowable RH
C     -------------------------------------------------
C
      RHCHK = RH
      IF (RH .LT. 0.) THEN
         RHCHK = -999.
         RETURN
         ENDIF
      IF (RH .GT. 100.) THEN
         RHCHK = 100.
         RETURN
         ENDIF
      IF (RH.LT.RHMIN) RHCHK = RHMIN
C
      RETURN
      END
C
C
C
C     -------------------------------------------------
      SUBROUTINE GET_62626(LNWMO)
C
C     Appends comment line (62626) to WMO message
C     -------------------------------------------------
C
      PARAMETER (MXRC = 9000, ninvar = 12)
C 
      DIMENSION IBGN(100), IEND(100)
      CHARACTER*1 ANS, OPTD, ALAT, ALON, OPTWD
      CHARACTER*4 VERSION, STNS, ICAO
      CHARACTER*20 MISSIONID
      CHARACTER*80 REMARK
      LOGICAL AUTOTEMP, LOCALACD
C
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP, LOCALACD
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /TEMPCTL/ MISSIONID, ICAO, RHMIN, REMARK
C
      DATA BAD/-999./
C
C
C     Set parameters
C     --------------
      MAXLENGTH = LNWMO-6
      OPTD = 'Q'
      REMARK = ' '
      LC = 1
C
C
C     Get selection
C     -------------
100   WRITE(LUT,'(//," Comment menu options:")')
      WRITE(LUT,'(" ------------------------------------")')
      WRITE(LUT,'(" (W) - Present weather")')
      WRITE(LUT,'(" (S) - Splash location")')
      WRITE(LUT,'(" (H) - Height of lowest wind")')
      WRITE(LUT,'(" (T) - Sea surface temperature [sst]")')
      WRITE(LUT,'(" (R) - Retransmitted message [ob #]")')
      WRITE(LUT,'(" (L) - Last report")')
      WRITE(LUT,'(" (Q) - Quit comments")')
      WRITE(LUT,'(" ------------------------------------")')
      WRITE(LUT,'(/," Enter selection [",A1,"]: ",NN)') OPTD
      READ(LUT,'(A)') ANS
      IF (ANS.EQ.' ') ANS = OPTD
      IF (ANS.EQ.'W' .OR. ANS.EQ.'w') GOTO 200
      IF (ANS.EQ.'S' .OR. ANS.EQ.'s') GOTO 300
      IF (ANS.EQ.'T' .OR. ANS.EQ.'t') GOTO 400
      IF (ANS.EQ.'L' .OR. ANS.EQ.'l') GOTO 500
      IF (ANS.EQ.'R' .OR. ANS.EQ.'r') GOTO 600
      IF (ANS.EQ.'H' .OR. ANS.EQ.'h') GOTO 700
      IF (ANS.EQ.'Q' .OR. ANS.EQ.'q') GOTO 9000
      GOTO 100
C
C
C     Specify present weather
C     -----------------------
200   CONTINUE
210   OPTWD = 'Q'
      WRITE(LUT,'(//," Present weather options:")')
      WRITE(LUT,'(" ------------------------------------")')
      WRITE(LUT,'(" (E) - Eye")')
      WRITE(LUT,'(" (W) - Eyewall [azimuth]")')
      WRITE(LUT,'(" (R) - Rainband")')
      WRITE(LUT,'(" (Q) - Quit (none apply)")')
      WRITE(LUT,'(" ------------------------------------")')
      WRITE(LUT,'(/," Enter selection [",A1,"]: ",NN)') OPTWD
      READ(LUT,'(A)') ANS
      IF (ANS.EQ.' ') ANS = OPTWD
      IF (ANS.EQ.'W' .OR. ANS.EQ.'w') GOTO 220
      IF (ANS.EQ.'E' .OR. ANS.EQ.'e') GOTO 230
      IF (ANS.EQ.'R' .OR. ANS.EQ.'r') GOTO 240
      GOTO 210
C
220   WRITE(LUT,'(/" Enter approximate azimuth of drop (deg): ",NN)')
      READ(LUT,*) AZIM
      IAZIM = NINT(AZIM)
      IF (IAZIM.GT.359 .OR. IAZIM.LT.0) GOTO 220
      LOPT = 11
      IF (LC+LOPT-1 .GT. MAXLENGTH) GOTO 9910
      WRITE(REMARK(LC:LC+LOPT-1),'(A8,I3.3)') 'EYEWALL ',IAZIM
      GOTO 7000
C
230   LOPT = 3
      IF (LC+LOPT-1 .GT. MAXLENGTH) GOTO 9910
      WRITE(REMARK(LC:LC+LOPT-1),'(A3)') 'EYE'
      GOTO 7000
C
240   LOPT = 8
      IF (LC+LOPT-1 .GT. MAXLENGTH) GOTO 9910
      WRITE(REMARK(LC:LC+LOPT-1),'(A8)') 'RAINBAND'
      GOTO 7000
C
C
C     Specify splash location
C     -----------------------
300   LOPT = 15
      IF (LC+LOPT-1 .GT. MAXLENGTH) GOTO 9910
C
      IF (SNDDAT(IXEND,11).EQ.BAD) THEN
         WRITE(LUT,'(/" Splash location not available.")')
         GOTO 100
         ENDIF
C
      ILAT = ABS(NINT(SNDDAT(IXEND,11)*100.))
      ILON = ABS(NINT(SNDDAT(IXEND,12)*100.))
      ALAT = 'N'
      IF (RLAT.LT.0.) ALAT = 'S'
      ALON = 'W'
      IF (RLON.LT.0.) ALON = 'W'
      WRITE(REMARK(LC:LC+LOPT-1),'(A4,I4.4,A1,I5.5,A1)') 'SPL ',
     *      ILAT,ALAT,ILON,ALON
      GOTO 7000
C
C
C     Specify SST
C     -----------
400   LOPT = 7
      IF (LC+LOPT-1 .GT. MAXLENGTH) GOTO 9910
C
      WRITE(LUT,'(/" Enter SST (C): ",NN)')
      READ(LUT,*) SST
      ISST = NINT(SST*10.)
      WRITE(REMARK(LC:LC+LOPT-1),'(A4,I3.3)') 'SST ',ISST
      GOTO 7000
C
C
C     Last report
C     -----------
500   LOPT = 19
      IF (LC+LOPT-1 .GT. MAXLENGTH) GOTO 9910
      WRITE(REMARK(LC:LC+LOPT-1),'(A15,A4)') 'LAST REPORT TO ',ICAO
      GOTO 7000
C
C
C     Retransmission of previous ob
C     -----------------------------
600   LOPT = 14
      IF (LC+LOPT-1 .GT. MAXLENGTH) GOTO 9910
C
      WRITE(LUT,'(/" Enter ob # being retransmitted: ",NN)')
      READ(LUT,*) NOB
      WRITE(REMARK(LC:LC+LOPT-1),'(A12,I2.2)') 'REXMT OF OB ',NOB
      GOTO 7000
C
C
C     Height of lowest transmitted wind
C     ---------------------------------
700   LOPT = 11
      IF (LC+LOPT-1 .GT. MAXLENGTH) GOTO 9910
      CALL WINDZINT(10.,WD10,WS10,WQ10,M10)
      IF (WD10.GE.0.) THEN
         GALAST = 10.
         ELSE
         CALL GAP(SNDDAT(1,7),IXEND,IBGN,IEND,NGAPS,BAD) 
         IE = IEND(NGAPS+1)
         GALAST = SNDDAT(IE,9)
         ENDIF
      IF (NINT(GALAST).GE.1000) GOTO 9920
      WRITE(REMARK(LC:LC+LOPT-1),'(A8,I3.3)') 'LST WND ',NINT(GALAST)
      GOTO 7000
C
C
C     Echo comment to screen, add blank at end
C     ----------------------------------------
7000  WRITE(LUT,'(/" Appending comment: ",A)') REMARK(LC:LC+LOPT)
      LC = LC + LOPT
      REMARK(LC:LC) = ' '
      LC = LC+1
      GOTO 100
C
C
C     All done, return
C     ----------------
9000  RETURN
C
C
C     Not enough room
C     ---------------
9910  WRITE(LUT,'(/" *** NOT ENOUGH ROOM FOR THIS COMMENT ***")')
      GOTO 100
C
C
C     Inappropriate comment
C     ---------------------
9920  WRITE(LUT,'(/" *** COMMENT INAPPROPRIATE OR UNNECESSARY ***")')
      GOTO 100
      END
C
C
C
C     -------------------------------------------------
c      SUBROUTINE CODE_62626
C
C     Appends comment line (62626) to WMO message.
C     -------------------------------------------------
C
c      PARAMETER (MXRC = 2400, ninvar = 12)
C 
c      CHARACTER*4 ICAO
c      CHARACTER*5 A62626
c      CHARACTER*20 MISSIONID
c      CHARACTER*80 REMARK
c      LOGICAL HAVE_COMMENT
C
c      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
c     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
c      COMMON /TEMPCTL/ MISSIONID, ICAO, RHMIN, REMARK
C
c      DATA A62626/'62626'/
C
C
c      HAVE_COMMENT = .FALSE.
c      DO 110 I = 1,80
c         IF (REMARK(I:I).EQ.'*') GOTO 900
c         IF (REMARK(I:I).NE.' ') HAVE_COMMENT = .TRUE.
c110      CONTINUE
C
c      IF (HAVE_COMMENT) WRITE(LUFO,'(A5,1X,A)') A62626,REMARK
c      RETURN
C
C
C     Bad character(s) present
C     ------------------------
c900   WRITE(LUT,'(/,A,
c     *  " *** COMMENT LINE FLAWED - NOT ACCEPTED ***"/)') CHAR(7)
c      RETURN
C
c      END

      SUBROUTINE CODE_62626()

	use thread_common

      PARAMETER (LNWMO = 66)
C 
C     Routine encodes mission ID line.
C     ----------------------------------------------------- 
C 
C 
      CHARACTER*1 BLANK
      DATA BLANK/' '/ 

      DO 155 I = lnwmo-6,1,-1
155      IF (group62626_thr(I:I).NE.BLANK 
	*	.and. group62626_thr(I:I).NE.char(0)) GOTO 160
C
160   ML = I

	nMsgLines_thr = nMsgLines_thr + 1
	temp_msg_thr(nMsgLines_thr)(1:5) = '62626'
	temp_msg_thr(nMsgLines_thr)(6:6) = BLANK
	temp_msg_thr(nMsgLines_thr)(7:7+ML-1) = group62626_thr(1:ml)
	do i = 7+ml, 100
		temp_msg_thr(nMsgLines_thr)(i:i) = blank
	end do


C     All done, return
C     ----------------
c900   RETURN
C
C
C     Error returns
C     -------------
c905   WRITE(LUT,'(/," *** INVALID NPLTFORM IN CODE_61616 ***",/)')
c      CLOSE(LUFX)
c      RETURN
C 
c910   WRITE(LUT,'(/," *** ERROR ON AIRCRAFT CONTROL FILE ***",/)')
      RETURN
      END
C
C
C
C     ------------------------------------------------------
      SUBROUTINE FINDTROP(NTROP,IXTROP)
C
C     Finds tropopause in sounding data.  
C
C     NTROP     Number of tropopauses found.
C     IXTROP    Index number of each one.
C     ------------------------------------------------------
C
      PARAMETER (MXRC = 9000, ninvar = 12)
      PARAMETER (MXT = 10)
C 
      DIMENSION IXTROP(MXT)
      CHARACTER*4 VERSION, STNS
      LOGICAL AUTOTEMP, LOCALACD, ABOVE500
C
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP, LOCALACD
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
C
      DATA BAD/-999./
C
C
C     Set parameters
C     --------------
      NTROP = 0
      DELTAP = 2.5
      NINC = NELM(DELTAP,FALLRATE,DATARATE)
      NINC = 1
      IBGN = IXEND
C
C
C     Valid sounding?
C     ---------------
      IF (PRL.EQ.BAD .OR. PRL.GT.200.0) RETURN
C
C
C     Begin search from 500 mb upward
C     -------------------------------
100   ABOVE500 = .TRUE.
      DO 110 I = IBGN,1,-1
         IF (SNDDAT(I,2).EQ.BAD) GOTO 110
         IF (SNDDAT(I,2).GT.500.) GOTO 110
         IF (SNDDAT(I,2).LT.30.) GOTO 110
         NA = I - NINC
         IF (NA.LT.1) NA = 1
         NB = I + NINC
         IF (NB.GT.IXEND) NB = IXEND
C
C        Calculate lapse rate (deg/m)
C        ----------------------------
         IF (SNDDAT(I,9).EQ.BAD .OR. SNDDAT(I,3).EQ.BAD) GOTO 110
         IF (SNDDAT(NA,9).EQ.BAD .OR. SNDDAT(NA,3).EQ.BAD) GOTO 110
         IF (SNDDAT(NB,9).EQ.BAD .OR. SNDDAT(NB,3).EQ.BAD) GOTO 110
         TDELTA = SNDDAT(NA,3)-SNDDAT(I,3)
         ZDELTA = (SNDDAT(NA,9)-SNDDAT(I,9))/1000.
         IF (ZDELTA.EQ.0.) GOTO 110
         TLAPSA = -TDELTA/ZDELTA
         TDELTA = SNDDAT(I,3)-SNDDAT(NB,3)
         ZDELTA = (SNDDAT(I,9)-SNDDAT(NB,9))/1000.
         IF (ZDELTA.EQ.0.) GOTO 110
         TLAPSB = -TDELTA/ZDELTA
C
C        Found lapse rate < 2.0 deg/km.  Check lapse rate 2 km above
C        -----------------------------------------------------------
         IF (TLAPSA.LT.2.0 .AND. TLAPSB.GE.2.0) THEN
            ZPTROP = SNDDAT(I,9)
            TPTROP = SNDDAT(I,3)
            IF (ZPTROP.EQ.BAD .OR. TPTROP.EQ.BAD) GOTO 110
            Z2KM = ZPTROP+2000.
            DO 120 J = I-1,1,-1
               Z = SNDDAT(J,9)
               T = SNDDAT(J,3)
               IF (Z.GE.Z2KM) GOTO 130
               IF (Z.EQ.BAD .OR. T.EQ.BAD) GOTO 110
               IF (Z-ZPTROP.EQ.0.) GOTO 120
               TLAPS = -(T-TPTROP)/((Z-ZPTROP)/1000.)
               IF (TLAPS.GT.2.) GOTO 110
120            CONTINUE
            GOTO 110
C
C           Survived check of lapse rate 2 km above.  If multiple trop,
C           need to make sure that 1 km below exceeds 2.0
C           -----------------------------------------------------------
130         IF (NTROP.EQ.0) GOTO 1000
            ZPTROP = SNDDAT(I,9)
            TPTROP = SNDDAT(I,3)
            IF (ZPTROP.EQ.BAD .OR. TPTROP.EQ.BAD) GOTO 110
            Z1KM = ZPTROP-1000.
            DO 140 J = I+1,IXEND
               Z = SNDDAT(J,9)
               T = SNDDAT(J,3)
               IF (Z.GE.Z1KM) GOTO 140
               IF (Z.EQ.BAD .OR. T.EQ.BAD) GOTO 110
               TLAPS = -(T-TPTROP)/((Z-ZPTROP)/1000.)
               IF (TLAPS.GT.2.) THEN
                  GOTO 1000
                  ELSE
                  GOTO 110
                  ENDIF
140            CONTINUE
C
            ENDIF
C
110      CONTINUE
C
C
C     No lapse rates less than 2 deg/km found above 500 mb
C     Search below 500 mb with different depth restriction.
C     This search is from 500 mb downward.
C     -----------------------------------------------------
200   IBGN = 1
205   ABOVE500 = .FALSE.
      DO 210 I = IBGN,IXEND
         IF (SNDDAT(I,2).EQ.BAD) GOTO 210
         IF (SNDDAT(I,2).LT.500.) GOTO 210
         NA = I - NINC
         IF (NA.LT.1) NA = 1
         NB = I + NINC
         IF (NB.GT.IXEND) NB = IXEND
C
C        Calculate lapse rate (deg/m)
C        ----------------------------
         IF (SNDDAT(I,9).EQ.BAD .OR. SNDDAT(I,3).EQ.BAD) GOTO 210
         IF (SNDDAT(NA,9).EQ.BAD .OR. SNDDAT(NA,3).EQ.BAD) GOTO 210
         IF (SNDDAT(NB,9).EQ.BAD .OR. SNDDAT(NB,3).EQ.BAD) GOTO 210
         TDELTA = SNDDAT(NA,3)-SNDDAT(I,3)
         ZDELTA = (SNDDAT(NA,9)-SNDDAT(I,9))/1000.
         IF (ZDELTA.EQ.0.) GOTO 210
         TLAPSA = -TDELTA/ZDELTA
         TDELTA = SNDDAT(I,3)-SNDDAT(NB,3)
         ZDELTA = (SNDDAT(I,9)-SNDDAT(NB,9))/1000.
         IF (ZDELTA.EQ.0.) GOTO 210
         TLAPSB = -TDELTA/ZDELTA
C
C        Found lapse rate < 2.0 deg/km.  Check lapse rate 1 km above
C        -----------------------------------------------------------
         IF (TLAPSA.LT.2.0 .AND. TLAPSB.GE.2.0) THEN
            ZPTROP = SNDDAT(I,9)
            TPTROP = SNDDAT(I,3)
            IF (ZPTROP.EQ.BAD .OR. TPTROP.EQ.BAD) GOTO 210
            Z1KM = ZPTROP+1000.
            DO 220 J = I-1,1,-1
               Z = SNDDAT(J,9)
               T = SNDDAT(J,3)
               IF (Z.GE.Z1KM) GOTO 250
               IF (Z.EQ.BAD .OR. T.EQ.BAD) GOTO 210
               IF (Z-ZPTROP.EQ.0.) GOTO 220
               TLAPS = -(T-TPTROP)/((Z-ZPTROP)/1000.)
               IF (TLAPS.GT.3.) GOTO 210
220            CONTINUE
            ENDIF
C
210      CONTINUE
      GOTO 1500
C
C
C     Tentatively found a low level tropopause
C     Make sure no 1 km layer has lapse rate exceeding 3 deg/km
C     ---------------------------------------------------------
250   DO 260 J = I-1,1,-1
         Z = SNDDAT(J,9)
         T = SNDDAT(J,3)
         IF (Z.EQ.BAD .OR. T.EQ.BAD) GOTO 1500
         Z1KM = Z+1000.
         DO 270 K = J,1,-1
            ZT = SNDDAT(K,9)
            TT = SNDDAT(K,3)
            IF (ZT.EQ.BAD .OR. TT.EQ.BAD) GOTO 1500
            IF (ZT.LT.Z1KM .AND. K.NE.1) GOTO 270
            TLAPS = -(TT-T)/((ZT-Z)/1000.)
            IF (TLAPS.GT.3.) GOTO 1500
            GOTO 260
270         CONTINUE
260      CONTINUE
C
C     Level survived, count it
C     ------------------------
      GOTO 1000
C
C
C     Found one!
C     ----------
1000  NTROP = NTROP+1
      IXTROP(NTROP) = I
      IF (NTROP.EQ.MXT) RETURN
C
C     Need to decide where to begin search for next one
C     For upper trops, must move at least 1 km upwards.
C     -------------------------------------------------
      IF (ABOVE500) THEN
         ZTROP = SNDDAT(I,9)
         DO 1050 L = I-1,1,-1
            IF (SNDDAT(L,9).EQ.BAD) GOTO 1050
            IF (SNDDAT(L,9).GT.ZTROP+1000.) THEN
               IBGN = L
               GOTO 100
               ENDIF
1050        CONTINUE
         RETURN
         ELSE
         IBGN = I+1
         GOTO 205
         ENDIF
C
C
C     None found
C     ----------
1500  RETURN
      END
C
C
C
C     ------------------------------------------------
      SUBROUTINE MANOPAA(RLAT,RLON,AA)
C
C     Returns RECCO MANOP geographical identifier.
C     RLAT = DEG N
C     RLON = DEG W
C     ------------------------------------------------
C
      CHARACTER*2 AA
C
      AA = 'NT'
C
      IF (RLON.LT.-40.) THEN
         AA = 'PA'
         RETURN
         ENDIF
C
      IF ((RLAT.LE.8.75 .AND. RLON.GT.70.0) .OR.
     *    (RLAT.LE.10.0 .AND. RLON.GT.83.5) .OR.
     *    (RLAT.LE.14.0 .AND. RLON.GT.84.5) .OR.
     *    (RLAT.LE.17.5 .AND. RLON.GT.90.0) .OR.
     *    (RLON.GT.100.)) THEN
         AA = 'PN'
         RETURN
         ENDIF
C
      RETURN
      END
C
C
C
C     ----------------------------------------------------- 
      SUBROUTINE GET_51515(MODE_WMO)
C 
C     Routine determines regional practice groups (51515)
C
C     MODE_WMO = 1 for NOAA encoding
C              = 2 for USAF encoding
C     ----------------------------------------------------- 
C 
      PARAMETER (MXRC = 9000, ninvar = 12)
      PARAMETER (NMANL = 12)
      PARAMETER (NDTMAX = 100)
      PARAMETER (NPROC = 50)
C 
      DIMENSION PRMANL(NMANL)
      LOGICAL AFRES, TERMEXT
      LOGICAL AG_51515, AG_10190, AG_10191, AG_10166, AG_10167
      CHARACTER*4 STNS, ICAO
      CHARACTER*20 MISSIONID
      CHARACTER*80 REMARK
C 
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /PROCESSING/ NPSTP, PROC(NPROC)
      COMMON /TEMPCTL/ MISSIONID, ICAO, RHMIN, REMARK
      COMMON /REGDATA/ AG_51515, AG_10190(2), AG_10191, AG_10166, 
     *                 AG_10167, PR_10190(2), GA_10190(2), PR_10191,
     *                 PZBTM, PZTOP, PDBTM(NDTMAX), PDTOP(NDTMAX),
     *                 NDT
C 
      DATA PRMANL/1070.,1000.,925.,850.,700.,500.,400.,
     *            300.,250.,200.,150.,100./
      DATA BAD/-999./ 
C
C
C     Do we follow standard NOAA or AFRES practices?
C     ----------------------------------------------
      AFRES = .FALSE.
      IF (MODE_WMO.EQ.2) AFRES = .TRUE.
C
C
C     Determine first and last height data for possible extrapolation
C     ---------------------------------------------------------------
      TERMEXT = .FALSE.
      AG_10190(1) = .FALSE.
      AG_10190(2) = .FALSE.
      AG_10191 = .FALSE.
C
      IF (SNDDAT(IXEND,9).EQ.0) GOTO 55
      DO 50 I = IXEND,1,-1
         IF (SNDDAT(I,9).NE.BAD) THEN
            TERMEXT = .TRUE.
            TERMPR = SNDDAT(I,2) 
            TERMGA = SNDDAT(I,9) 
            TERMTE = SNDDAT(I,3) 
            TERMRH = SNDDAT(I,4)
            GOTO 55
            ENDIF
50       CONTINUE
C
55    DO 60 I = 1,IXEND
         IF (SNDDAT(I,9).NE.BAD .AND. SNDDAT(I,9).NE.0.) THEN
            TOPPR = SNDDAT(I,2) 
            TOPGA = SNDDAT(I,9) 
            TOPTE = SNDDAT(I,3) 
            TOPRH = SNDDAT(I,4)
            GOTO 65
            ENDIF
60       CONTINUE
C
C
C     Do we extrapolate upwards
C     -------------------------
65    DO 70 I = 2,NMANL
         PRDIFF = PRL-PRMANL(I)
         IF (PRDIFF.GE.0.5 .AND. PRDIFF.LE.25.) THEN
            CALL EXTRAP_TEMP(TOPPR,PRMANL(I),TX,BAD)
            RHX = RHCK(TOPRH,.TRUE.,70.,BAD)
            GA = HYDROZ(TOPPR,TOPTE,RHX,TOPGA,
     *                  PRMANL(I),TX,RHX,2,BAD)
            IF (GA.NE.BAD) THEN
               AG_10190(1) = .TRUE.
               PR_10190(1) = PRMANL(I)
               GA_10190(1) = GA
               ENDIF
            GOTO 75
            ENDIF
70       CONTINUE
C
C
C     Do we extrapolate downwards
C     ---------------------------
75    IF (.NOT. TERMEXT) GOTO 100
      DO 80 I = 2,NMANL
         PRDIFF = PRMANL(I)-TERMPR
         IF (PRDIFF.GT.0. .AND. PRDIFF.LE.25.) THEN
            CALL EXTRAP_TEMP(TERMPR,PRMANL(I),TX,BAD)
            RHX = RHCK(TERMRH,.TRUE.,70.,BAD)
            GA = HYDROZ(TERMPR,TERMTE,RHX,TERMGA,
     *                  PRMANL(I),TX,RHX,2,BAD)
            IF (GA.NE.BAD) THEN
               AG_10190(2) = .TRUE.
               PR_10190(2) = PRMANL(I)
               GA_10190(2) = GA
               ENDIF
            GOTO 85
            ENDIF
80       CONTINUE
C
C
85    IF (TERMPR.LT.850.) GOTO 100
      RHX = RHCK(TERMRH,.TRUE.,70.,BAD)
      ESTSFCP = HYDROP(TERMPR,TERMTE,RHX,TERMGA,
     *                 TERMTE,RHX,0.,2,BAD)
      CALL EXTRAP_TEMP(TERMPR,ESTSFCP,TX,BAD)
      ESTSFCP = HYDROP(TERMPR,TERMTE,RHX,TERMGA,
     *                 TX,RHX,0.,2,BAD)
      IF (ESTSFCP.NE.BAD) THEN
         AG_10191 = .TRUE.
         PR_10191 = ESTSFCP
         ENDIF
C
C
C     Check for doubtful temperatures
C     -------------------------------
100   ISTART = 1
      NDT = 0
      AG_10167 = .FALSE.
200   DO 210 J=ISTART,IXEND
         IF (ISNDFLG(J,2).EQ.4) THEN
            AG_10167 = .TRUE.
            NDT = NDT+1
            IF (NDT.GT.NDTMAX) THEN
               WRITE(LUT,'(" *** TOO MANY DOUBTFUL LAYERS ***")')
               GOTO 230
               ENDIF
            PDTOP(NDT) = SNDDAT(J,2)
            DO 220 L = J,IXEND
               IF (ISNDFLG(L,2).NE.4) THEN
                  PDBTM(NDT) = SNDDAT(L-1,2)
                  ISTART = L+1
                  IF (ISTART.GT.IXEND) GOTO 230
                  GOTO 200
                  ENDIF
220            CONTINUE
               PDBTM(NDT) = SNDDAT(IXEND,2)
               GOTO 230
            ENDIF
210      CONTINUE
C
C
C     Check for doubtful geopotentials (AFRES ONLY).
C     Data marked doubtful if T doubtful or missing,
C     also if sfc was hydrostatic anchor.
C     ----------------------------------------------
230   AG_10166 = .FALSE.
      IF (.NOT. AFRES) GOTO 290
      PZTOP = 1100.
      PZBTM = 0.
      DO 235 J = 1,IXEND
         IF (SNDDAT(J,9) .NE. BAD) THEN
            IF (SNDDAT(J,2).LT.PZTOP) PZTOP = SNDDAT(J,2)
            IF (SNDDAT(J,2).GT.PZBTM) PZBTM = SNDDAT(J,2)
            ENDIF
235      CONTINUE
      IF (PZTOP.EQ.1100. .OR. PZBTM.EQ.0.) GOTO 290
C
      IF (PROC(7).EQ.2.) THEN
         AG_10166 = .TRUE.
         GOTO 290
         ENDIF         
      IF (AG_10167 .AND. PROC(7).EQ.1) THEN
         AG_10166 = .TRUE.
         IF (PZTOP.LT.PDTOP(1)) PZTOP = PDTOP(1)
         GOTO 290
         ENDIF
      IF (PROC(7).EQ.1) THEN
         DO 240 J = 1,IXEND
            IF (SNDDAT(J,3).EQ.BAD) THEN
               AG_10166 = .TRUE.
               IF (PZTOP.LT.SNDDAT(J,2)) PZTOP = SNDDAT(J,2)
               GOTO 290
               ENDIF
240         CONTINUE
         GOTO 290
         ENDIF
C
290   CONTINUE
C
C
C     All done, return
C     ----------------
500   RETURN
      END
C
C
C
C     ----------------------------------------------------- 
      SUBROUTINE CODE_51515(CODE,LINE,ICNT)
C 
C     Routine encodes regional practice groups (51515)
C     ----------------------------------------------------- 
C 
      PARAMETER (MXRC = 9000)
      PARAMETER (LNWMO = 66)
      PARAMETER (NDTMAX = 100)
C 
      CHARACTER*1 BLANK
      CHARACTER*20 MISSIONID
      CHARACTER*80 REMARK
      CHARACTER*(LNWMO) CODE((MXRC*2)/(LNWMO/12))   
C
      LOGICAL AG_51515, AG_10190, AG_10191, AG_10166, AG_10167
C 
      COMMON /TEMPCTL/ MISSIONID, ICAO, RHMIN, REMARK
      COMMON /REGDATA/ AG_51515, AG_10190(2), AG_10191, AG_10166, 
     *                 AG_10167, PR_10190(2), GA_10190(2), PR_10191,
     *                 PZBTM, PZTOP, PDBTM(NDTMAX), PDTOP(NDTMAX),
     *                 NDT
C 
      DATA BLANK/' '/ 
      DATA BAD/-999./ 
C
C
C     Required formats
C     ----------------
2     FORMAT(I2.2)
3     FORMAT(I3.3)
C
C
C     Do we have any valid groups?
C     ----------------------------
      AG_51515 = .FALSE.
      IF (AG_10191    .OR. 
     *    AG_10190(1) .OR. 
     *    AG_10190(2) .OR.
     *    AG_10166    .OR.
     *    AG_10167)    AG_51515 = .TRUE.
      IF (.NOT. AG_51515) RETURN
C
      IF (ICNT.EQ.1) THEN 
         LINE = LINE - 1
         ELSE 
         DO 630 I = ICNT,LNWMO
            CODE(LINE)(I:I) = BLANK 
630         CONTINUE
         ENDIF
C
      LINE = LINE+1
      ICNT = 1 
      CODE(LINE)(ICNT:ICNT+4) = '51515'
      CODE(LINE)(ICNT+5:ICNT+5) = BLANK
      CALL INCR(ICNT,LINE,LNWMO)
C
      IF (AG_10166) THEN
         CODE(LINE)(ICNT:ICNT+4) = '10166'
         CODE(LINE)(ICNT+5:ICNT+5) = BLANK
         CALL INCR(ICNT,LINE,LNWMO)
         IP1 = NINT(PZBTM/10.)
         IF (IP1.GE.100) IP1 = IP1-100
         IP2 = NINT(PZTOP/10.)
         IF (IP2.GE.100) IP2 = IP2-100
         CODE(LINE)(ICNT:ICNT) = '0'
         WRITE(CODE(LINE)(ICNT+1:ICNT+2),2) IP1
         IF (PZBTM.EQ.BAD) CODE(LINE)(ICNT+1:ICNT+2) = '//'
         WRITE(CODE(LINE)(ICNT+3:ICNT+4),2) IP2
         IF (PZTOP.EQ.BAD) CODE(LINE)(ICNT+3:ICNT+4) = '//'
         CODE(LINE)(ICNT+5:ICNT+5) = BLANK
         CALL INCR(ICNT,LINE,LNWMO)
         ENDIF
C
      IF (AG_10167) THEN
         CODE(LINE)(ICNT:ICNT+4) = '10167'
         CODE(LINE)(ICNT+5:ICNT+5) = BLANK
         CALL INCR(ICNT,LINE,LNWMO)
         DO 640 I = 1,NDT
            IP1 = NINT(PDBTM(I)/10.)
            IF (IP1.GE.100) IP1 = IP1-100
            IP2 = NINT(PDTOP(I)/10.)
            IF (IP2.GE.100) IP2 = IP2-100
            CODE(LINE)(ICNT:ICNT) = '0'
            WRITE(CODE(LINE)(ICNT+1:ICNT+2),2) IP1
            IF (PDBTM(I).EQ.BAD) CODE(LINE)(ICNT+1:ICNT+2) = '//'
            WRITE(CODE(LINE)(ICNT+3:ICNT+4),2) IP2
            IF (PDTOP(I).EQ.BAD) CODE(LINE)(ICNT+3:ICNT+4) = '//'
            CODE(LINE)(ICNT+5:ICNT+5) = BLANK
            CALL INCR(ICNT,LINE,LNWMO)
640         CONTINUE
         ENDIF
C
      DO 670 I = 1,2
         IF (AG_10190(I)) THEN
            CODE(LINE)(ICNT:ICNT+4) = '10190'
            CODE(LINE)(ICNT+5:ICNT+5) = BLANK
            CALL INCR(ICNT,LINE,LNWMO)
            IP = PR_10190(I)/10. 
            IF (IP.EQ.100) IP = 0 
            IF (GA_10190(I).LT.0.) GA_10190(I) = 500.+ABS(GA_10190(I))
            IHT = NINT(GA_10190(I))
            IF (PR_10190(I).LE.500.) IHT = NINT(GA_10190(I)/10) 
            IHT = MOD(IHT,1000) 
            WRITE(CODE(LINE)(ICNT:ICNT+1),2) IP
            IF (PR_10190(I).EQ.BAD) CODE(LINE)(ICNT:ICNT+1) = '//'
            WRITE(CODE(LINE)(ICNT+2:ICNT+4),3) IHT 
            IF (GA_10190(I).EQ.BAD) CODE(LINE)(ICNT+2:ICNT+4) = '///' 
            CODE(LINE)(ICNT+5:ICNT+5) = BLANK
            CALL INCR(ICNT,LINE,LNWMO)
            ENDIF
670      CONTINUE
C
      IF (AG_10191) THEN
         CODE(LINE)(ICNT:ICNT+4) = '10191'
         CODE(LINE)(ICNT+5:ICNT+5) = BLANK
         CALL INCR(ICNT,LINE,LNWMO)
         ENDIF
C
C
C     All done, return
C     ----------------
900   RETURN
      END
C
C
C
C     ----------------------------------------------------- 
      SUBROUTINE CODE_61616(CODE,LINE,ICNT)

	use thread_common

C 
C     Routine encodes mission ID line.
C     ----------------------------------------------------- 
C 
      PARAMETER (MXRC = 9000)
      PARAMETER (LNWMO = 66)
C 
      CHARACTER*1 BLANK
      CHARACTER*4 ICAO
      CHARACTER*5 WMOID
      CHARACTER*20 MISSIONID
      CHARACTER*80 REMARK
      CHARACTER*(LNWMO) CODE((MXRC*2)/(LNWMO/12))   
C 
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /TEMPCTL/ MISSIONID, ICAO, RHMIN, REMARK
 
      DATA BLANK/' '/ 
C
C
C     Begin on new line
C     -----------------
100   IF (ICNT.NE.1) THEN
         DO 150 I = ICNT,LNWMO
150         CODE(LINE)(I:I) = BLANK 
            LINE = LINE+1
            ICNT = 1 
         ENDIF
C
      DO 155 I = 100,1,-1
155      IF (group61616_thr(I:I).NE.BLANK
	*		.and. group61616_thr(I:I).NE.char(0)) GOTO 160
C
160   ML = I
c      OPEN(LUFX,FILE='acft_ctrl.dat',STATUS='OLD',ERR=910)
c165   READ(LUFX,'(I3,1X,A5)',ERR=910,END=905)
c     *     NAC,WMOID
c      IF (NPLTFORM.EQ.NAC) THEN
c         CODE(LINE)(7:11) = WMOID
c         GOTO 170
c         ENDIF
c      GOTO 165
C
c170   CLOSE(LUFX)
      CODE(LINE)(1:5) = '61616'
      CODE(LINE)(6:6) = BLANK
c      CODE(LINE)(12:12) = BLANK
      CODE(LINE)(7:7+ML-1) = group61616_thr(1:ml)
      ICNT = 7+ML
	 do i = 7+ml, lnwmo
		code(line)(i:i) = blank
	 end do

	icnt = lnwmo
c      CODE(LINE)(ICNT:ICNT) = BLANK
c      ICNT = ICNT+1
c      CODE(LINE)(ICNT:ICNT+4) = 'OB 00'
c      ICNT = ICNT+5
c      CODE(LINE)(ICNT:ICNT) = BLANK
c      ICNT = ICNT+1
c      CODE(LINE)(ICNT:ICNT+3) = ICAO
c      ICNT = ICNT+4
C
C
C     All done, return
C     ----------------
900   RETURN
C
C
C     Error returns
C     -------------
905   WRITE(LUT,'(/," *** INVALID NPLTFORM IN CODE_61616 ***",/)')
      CLOSE(LUFX)
      RETURN
C 
910   WRITE(LUT,'(/," *** ERROR ON AIRCRAFT CONTROL FILE ***",/)')
      RETURN
      END

C
C
C
C     ------------------------------------------
      SUBROUTINE WINDZINT(ZIN,WD,WS,WQ,M) 
C 
C     Extracts wind at Z m.
C     ------------------------------------------
C 
      PARAMETER (MXRC = 9000, ninvar = 12)
C
      LOGICAL AUTOTEMP, LOCALACD
      CHARACTER*4 STNS
      CHARACTER*4 VERSION
C 
      COMMON /PROGCTL/ NSITE, VERSION, AUTOTEMP, LOCALACD
      COMMON /PARAM/ LUT,LUFI,LUFX,LUFW,LUFP,LUFO,LUPR,
     *               NPLTFORM,NSNDTYPE,FALLRATE,DATARATE
      COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)
      COMMON /LAUNCH/ IYR,IMO,IDY,RLAT,RLON,TIML,PRL,TEL,RHL,DPL,HTL, 
     *                WDL,WSL,STNS 
      COMMON /WORK/ X(MXRC),Y(MXRC),Z(MXRC) 
C
      DATA BAD/-999./
C 
C
C     Default values are missing
C     --------------------------
      WS = BAD
      WD = BAD
      WQ = BAD
C
C
C     Convert WD,WS into U,V for interpolation
C     ----------------------------------------
100   DO 110 L = 1,IXEND
         X(L) = UCMP(SNDDAT(L,6),SNDDAT(L,7))
         Y(L) = VCMP(SNDDAT(L,6),SNDDAT(L,7))
110      CONTINUE 
C 
C 
C     Interpolate from data arrays
C     ----------------------------
      CALL POLATE(IXEND,SNDDAT(1,9),X,ZIN,UNOW,M,BAD)
      CALL POLATE(IXEND,SNDDAT(1,9),Y,ZIN,VNOW,M,BAD)
      CALL POLATE(IXEND,SNDDAT(1,9),SNDDAT(1,8),ZIN,WQNOW,M,BAD)
      IF (UNOW.EQ.BAD .OR. VNOW.EQ.BAD) RETURN
      WS = WSCOMP(UNOW,VNOW) 
      WD = WDCOMP(UNOW,VNOW) 
      WQ = WQNOW
      RETURN
      END 
C
C

     
     
      subroutine sndTransfer(data, n, index, badVal)

	PARAMETER (MXRC = 9000, NINVAR = 12)

	COMMON /SOUNDING/ ID,IXEND,SNDDAT(MXRC,NINVAR),ISNDFLG(MXRC,4)

		integer n
		integer index
		real(4) data(n)

		integer i

		do i = 1, n
			pt = data(i)
			if (pt .LT. -998.0) pt = badVal
			SNDDAT(n-i+1,INDEX) = pt
		end do
	
      end subroutine 

