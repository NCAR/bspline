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
