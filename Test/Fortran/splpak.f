C PACKAGE SPLPAK         DOCUMENTATION FOR USER ENTRIES FOLLOWS
C                        THE GENERAL PACKAGE INFORMATION.
C
C LATEST REVISION        MARCH  1985
C
C PURPOSE                THIS PACKAGE CONTAINS ROUTINES FOR FITTING
C                        (LEAST SQUARES) A MULTIDIMENSIONAL CUBIC SPLINE
C                        TO ARBITRARILY LOCATED DATA.  IT ALSO CONTAINS
C                        ROUTINES FOR EVALUATING THIS SPLINE (OR ITS
C                        PARTIAL DERIVATIVES) AT ANY POINT.
C
C                        COEFFICIENT CALCULATION IS PERFORMED IN
C                        SUBROUTINES SPLCC OR SPLCW AND EVALUATION IS
C                        PERFORMED BY FUNCTIONS SPLFE OR SPLDE.
C
C USAGE                  PACKAGE SPLPAK CONTAINS FOUR USER ENTRIES --
C                        SPLCC, SPLCW, SPLFE, AND SPLDE.
C
C                        THE USER FIRST CALLS SPLCC BY
C
C                          CALL SPLCC (NDIM,XDATA,L1XDAT,YDATA,NDATA,
C                                      XMIN,XMAX,NODES,XTRAP,COEF,NCF,
C                                      WORK,NWRK,IERROR)
C
C                        OR SPLCW BY
C
C                          CALL SPLCW (NDIM,XDATA,L1XDATA,YDATA,WDATA,
C                                      NDATA,XMIN,XMAX,NODES,XTRAP,
C                                      COEF,NCF,WORK,NWRK,IERROR)
C
C                        THE PARAMETER NDATA IN THE CALL TO SPLCW
C                        ENABLES THE USER TO WEIGHT SONE OF THE DATA
C                        POINTS MORE HEAVILY THAN OTHERS.  BOTH
C                        ROUTINES RETURN A SET OF COEFFICIENTS IN THE
C                        ARRAY COEF.  THESE COEFFICIENTS ARE
C                        SUBSEQUENTLY USED IN THE COMPUTATION OF
C                        FUNCTION VALUES AND PARTIAL DERIVATIVES.
C                        THE USER THEN CALLS SPLFE OR SPLDE ANY
C                        NUMBER OF TIMES IN ANY ORDER PROVIDED THAT
C                        THE VALUES OF THE INPUTS, NDIM,COEF,XMIN,
C                        XMAX, AND NODES, ARE PRESERVED BETWEEN CALLS.
C                        THE ROUTINES ARE CALLED IN THE FOLLOWING WAY.
C
C                          F = SPLFE (NDIM,X,COEF,XMIN,XMAX,NODES,
C                                     IERROR)
C
C                        OR
C
C                          F = SPLDE (NDIM,X,NDERIV,COEF,XMIN,XMAX,
C                                     NODES,IERROR)
C
C                        THE ROUTINE, SPLFE,  RETURNS AN INTERPOLATED
C                        VALUE AT THE POINT DEFINED BY THE ARRAY X.
C                        XPLDE AFFORDS THE USER THE ADDITIONAL
C                        CAPABILITY OF CALCULATING AN INTERPOLATED
C                        VALUE FOR ONE OF SEVERAL PARTIAL DERIVATIVES
C                        SPECIFIED BY THE ARRAY NDERIV.
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE, EXCEPT FOR ERROR MESSAGES PRINTED BY
C                        CALLS TO ULIBER, IF AN ERROR IS DETECTED.
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       SUPRLS, ULIBER
C FILES
C
C LANGUAGE               FORTRAN
C
C HISTORY                DEVELOPED IN 1972-73 BY DAVE FULKER OF NCAR'S
C                        SCIENTIFIC COMPUTING DIVISION.
C
C PORTABILITY            FORTRAN 66
C***********************************************************************
C
C SUBROUTINE SPLCC (NDIM,XDATA,L1XDAT,YDATA,NDATA,XMIN,XMAX,NODES,
C                   XTRAP,COEF,NCF,WORK,NWRK,IERROR)
C
C DIMENSION OF           XDATA(L1XDAT,NDATA),YDATA(NDATA),XMIN(NDIM),
C ARGUMENTS              XMAX(NDIM),NODES(NDIM),COEF(NCF),WORK(NWRK)
C
C PURPOSE                N-DIMENSIONAL CUBIC SPLINE COEFFICIENT
C                        CALCULATION BY LEAST SQUARES.
C
C USAGE                  THE USAGE AND ARGUMENTS OF THIS ROUTINE ARE
C                        IDENTICAL TO THOSE FOR SPLCW EXCEPT FOR THE
C                        OMISSION OF THE ARRAY OF WEIGHTS, WDATA.  SEE
C                        ENTRY SPLCW DESCRIPTION IMMEDIATELY BELOW FOR A
C                        COMPLETE DESCRIPTION.
C
C                        CALL SPLCC (NDIM,XDATA,L1XDAT,YDATA,NDATA,XMIN,
C                                    XMAX,NODES,XTRAP,COEF,NCF,WORK,
C                                    NWRK,IERROR)
C***********************************************************************
C
C SUBROUTINE SPLCW (NDIM,XDATA,L1XDAT,YDATA,WDATA,NDATA,XMIN,XMAX,
C                   NODES,XTRAP,COEF,NCF,WORK,NWRK,IERROR)
C
C
C DIMENSION OF           XDATA(L1XDAT,NDATA),YDATA(NDATA),WDATA(NDATA),
C ARGUMENTS              XMIN(NDIM),XMAX(NDIM),NODES(NDIM),COEF(NCF),
C                        WORK(NWRK)
C
C PURPOSE                N-DIMENSIONAL CUBIC SPLINE COEFFICIENT
C                        CALCULATION BY WEIGHTED LEAST SQUARES ON
C                        ARBITRARILY LOCATED DATA.
C
C                        A GRID OF EVENLY SPACED NODES IN NDIM SPACE IS
C                        DEFINED BY THE ARGUMENTS XMIN, XMAX AND NODES.
C                        A LINEAR BASIS FOR THE CLASS OF NATURAL SPLINES
C                        ON THESE NODES IS FORMED, AND A SET OF
C                        CORRESPONDING COEFFICIENTS IS COMPUTED IN THE
C                        ARRAY COEF.  THESE COEFFICIENTS ARE CHOSEN TO
C                        MINIMIZE THE WEIGHTED SUM OF SQUARED ERRORS
C                        BETWEEN THE SPLINE AND THE ARBITRARILY LOCATED
C                        DATA VALUES DESCRIBED BY THE ARGUMENTS XDATA,
C                        YDATA AND NDATA.  THE SMOOTHNESS OF THE SPLINE
C                        IN DATA SPARSE AREAS IS CONTROLLED BY THE
C                        ARGUMENT XTRAP.
C
C NOTE                   IN ORDER TO UNDERSTAND THE ARGUMENTS OF THIS
C                        ROUTINE, ONE SHOULD REALIZE THAT THE NODE GRID
C                        NEED NOT BEAR ANY PARTICULAR RELATION TO THE
C                        DATA POINTS.  IN THE THEORY OF EXACT FIT
C                        INTERPOLATORY SPLINES THE NODES WOULD IN FACT
C                        BE DATA LOCATIONS, BUT IN THIS CASE THEY SERVE
C                        ONLY TO DEFINE THE CLASS OF SPLINES FROM WHICH
C                        THE APPROXIMATING FUNCTION IS CHOSEN.  THIS
C                        NODE GRID IS A RECTANGULAR ARRANGEMENT OF
C                        POINTS IN NDIM SPACE, WITH THE RESTRICTION THAT
C                        ALONG ANY COORDINATE DIRECTION THE NODES ARE
C                        EQUALLY SPACED.  THE CLASS OF NATURAL SPLINES
C                        ON THIS GRID OF NODES (NDIM-CUBIC SPLINES WHOSE
C                        2ND DERIVATIVES NORMAL TO THE BOUNDARIES ARE 0)
C                        HAS AS MANY DEGREES OF FREEDOM AS THE GRID HAS
C                        NODES.  THUS THE SMOOTHNESS OR FLEXIBILITY OF
C                        THE SPLINES IS DETERMINED BY THE CHOICE OF THE
C                        NODE GRID.
C
C USAGE                  CALL SPLCW (NDIM,XDATA,L1XDAT,YDATA,WDATA,
C                                    NDATA,XMIN,XMAX,NODES,XTRAP,COEF,
C                                    NCF,WORK,NWRK,IERROR)
C                        THE SPLINE (OR ITS DERIVATIVES) MAY THEN BE
C                        EVALUATED BY USING FUNCTION SPLFE (OR SPLDE).
C
C ARGUMENTS
C
C ON INPUT               NDIM
C                          THE DIMENSIONALITY OF THE PROBLEM.  THE
C                          SPLINE IS A FUNCTION OF NDIM VARIABLES OR
C                          COORDINATES AND THUS A POINT IN THE
C                          INDEPENDENT VARIABLE SPACE IS AN NDIM VECTOR.
C                          NDIM MUST BE IN THE RANGE 1 .LE. NDIM .LE. 4.
C
C                        XDATA
C                          A COLLECTION OF LOCATIONS FOR THE DATA
C                          VALUES, I.E., POINTS FROM THE INDEPENDENT
C                          VARIABLE SPACE.  THIS COLLECTION IS A
C                          2-DIMENSIONAL ARRAY WHOSE 1ST DIMENSION
C                          INDEXES THE NDIM COORDINATES OF A GIVEN POINT
C                          AND WHOSE 2ND DIMENSION LABELS THE DATA
C                          POINT.  FOR EXAMPLE, THE DATA POINT WITH
C                          LABEL IDATA IS LOCATED AT THE POINT
C                          (XDATA(1,IDATA),...,XDATA(NDIM,IDATA)) WHERE
C                          THE ELEMENTS OF THIS VECTOR ARE THE VALUES OF
C                          THE NDIM COORDINATES.  THE LOCATION, NUMBER
C                          AND ORDERING OF THE DATA POINTS IS ARBITRARY.
C                               THE DIMENSION IS ASSUMED TO BE
C                               XDATA(L1XDAT,NDATA).
C
C                        L1XDAT
C                          THE LENGTH OF THE 1ST DIMENSION OF XDATA IN
C                          THE CALLING PROGRAM.  L1XDAT MUST BE .GE.
C                          NDIM.
C                               NOTE:  FOR 1 DIMENSIONAL PROBLEMS L1XDAT
C                                      IS USUALLY 1.
C
C                        YDATA
C                          A COLLECTION OF DATA VALUES CORRESPONDING TO
C                          THE POINTS IN XDATA.  YDATA(IDATA) IS THE
C                          DATA VALUE ASSOCIATED WITH THE POINT
C                          (XDATA(1,IDATA),...,XDATA(NDIM,IDATA)) IN THE
C                          INDEPENDENT VARIABLE SPACE.  THE SPLINE WHOSE
C                          COEFFICIENTS ARE COMPUTED BY THIS ROUTINE
C                          APPROXIMATES THESE DATA VALUES IN THE LEAST
C                          SQUARES SENSE.
C                               THE DIMENSION IS ASSUMED TO BE
C                               YDATA(NDATA).
C
C                        WDATA
C                          A COLLECTION OF WEIGHTS.  WDATA(IDATA) IS A
C                          WEIGHT ASSOCIATED WITH THE DATA POINT
C                          LABELLED IDATA.  IT SHOULD BE NON-NEGATIVE,
C                          BUT MAY BE OF ANY MAGNITUDE.  THE WEIGHTS
C                          HAVE THE EFFECT OF FORCING GREATER OR LESSER
C                          ACCURACY AT A GIVEN POINT AS FOLLOWS--THIS
C                          ROUTINE CHOOSES COEFFICIENTS TO MINIMIZE THE
C                          SUM OVER ALL DATA POINTS OF THE QUANTITY
C                          (WDATA(IDATA)*(YDATA(IDATA) - SPLINE VALUE AT
C                          XDATA(IDATA)))**2.  THUS, IF THE RELIABILITY
C                          OF A DATA POINT IS KNOWN TO BE LOW, THE
C                          CORRESPONDING WEIGHT MAY BE MADE SMALL
C                          (RELATIVE TO THE OTHER WEIGHTS) SO THAT THE
C                          SUM OVER ALL DATA POINTS IS AFFECTED LESS BY
C                          DISCREPENCIES AT THE UNRELIABLE POINT.  DATA
C                          POINTS WITH ZERO WEIGHT ARE COMPLETELY
C                          IGNORED.
C                               NOTE:  IF WDATA(1) IS .LT. 0, THE OTHER
C                                      ELEMENTS OF WDATA ARE NOT
C                                      REFERENCED, AND ALL WEIGHTS ARE
C                                      ASSUMED TO BE UNITY.
C                          THE DIMENSION IS ASSUMED TO BE WDATA(NDATA)
C                          UNLESS WDATA(1) .LT. 0 IN WHICH CASE THE
C                          DIMENSION IS ASSUMED TO BE 1.
C
C                        NDATA
C                          THE NUMBER OF DATA POINTS MENTIONED IN THE
C                          ABOVE ARGUMENTS.
C
C                        XMIN
C                          A VECTOR DESCRIBING THE LOWER EXTREME CORNER
C                          OF THE NODE GRID.  A SET OF EVENLY SPACED
C                          NODES IS FORMED ALONG EACH COORDINATE AXIS
C                          AND XMIN(IDIM) IS THE LOCATION OF THE FIRST
C                          NODE ALONG THE IDIM AXIS.
C                               THE DIMENSION IS ASSUMED TO BE
C                               XMIN(NDIM).
C
C                        XMAX
C                          A VECTOR DESCRIBING THE UPPER EXTREME CORNER
C                          OF THE NODE GRID.  A SET OF EVENLY SPACED
C                          NODES IS FORMED ALONG EACH COORDINATE AXIS
C                          AND XMAX(IDIM) IS THE LOCATION OF THE LAST
C                          NODE ALONG THE IDIM AXIS.
C                               THE DIMENSION IS ASSUMED TO BE
C                               XMAX(NDIM).
C
C                        NODES
C                          A VECTOR OF INTEGERS DESCRIBING THE NUMBER OF
C                          NODES ALONG EACH AXIS.  NODES(IDIM) IS THE
C                          NUMBER OF NODES (COUNTING ENDPOINTS) ALONG
C                          THE IDIM AXIS AND DETERMINES THE FLEXIBILITY
C                          OF THE SPLINE IN THAT COORDINATE DIRECTION.
C                          NODES(IDIM) MUST BE .GE. 4 BUT MAY BE AS
C                          LARGE AS THE ARRAYS COEF AND WORK ALLOW.
C                               THE DIMENSION IS ASSUMED TO BE
C                               NODES(NDIM).
C
C                          NOTE:  THE NODE GRID IS COMPLETELY DEFINED BY
C                                 THE ARGUMENTS XMIN, XMAX AND NODES.
C                                 THE SPACING OF THIS GRID IN THE IDIM
C                                 COORDINATE DIRECTION IS
C                                 DX(IDIM) = (XMAX(IDIM)-XMIN(IDIM)) /
C                                             (NODES(IDIM)-1).
C                                 A NODE IN THIS GRID MAY BE INDEXED BY
C                                 AN NDIM VECTOR OF INTEGERS
C                                 (IN(1),...,IN(NDIM)) WHERE
C                                 1 .LE. IN(IDIM) .LE. NODES(IDIM).
C                                 THE LOCATION OF SUCH A NODE MAY BE
C                                 REPRESENTED BY AN NDIM VECTOR
C                                 (X(1),...,X(NDIM)) WHERE
C                                 X(IDIM) = XMIN(IDIM) + (IN(IDIM)-1) *
C                                                         DX(IDIM).
C
C                        XTRAP
C                          A PARAMETER TO CONTROL EXTRAPOLATION TO DATA
C                          SPARSE AREAS.  THE REGION DESCRIBED BY XMIN
C                          AND XMAX IS DIVIDED INTO RECTANGLES, THE
C                          NUMBER OF WHICH IS DETERMINED BY NODES, AND
C                          ANY RECTANGLE CONTAINING A DISPROPORTIONATELY
C                          SMALL NUMBER OF DATA POINTS IS CONSIDERED TO
C                          BE DATA SPARSE (RECTANGLE IS USED HERE TO
C                          MEAN NDIM-DIMENSIONAL RECTANGLE).  IF XTRAP
C                          IS NONZERO THE LEAST SQUARES PROBLEM IS
C                          ARGUMENTED WITH DERIVATIVE CONSTRAINTS IN THE
C                          DATA SPARSE AREAS TO PREVENT THE MATRIX FROM
C                          BECOMING POORLY CONDITIONED.  XTRAP SERVES AS
C                          A WEIGHT FOR THESE CONSTRAINTS, AND THUS MAY
C                          BE USED TO CONTROL SMOOTHNESS IN DATA SPARSE
C                          AREAS.  THE EXPERIENCE OF THE AUTHOR
C                          INDICATES THAT UNITY IS A GOOD FIRST GUESS
C                          FOR THIS PARAMETER, AND THAT ITS SIZE IS NOT
C                          VERY CRITICAL.
C
C                               NOTE:  IF XTRAP IS 0 SUBSTANTIAL
C                                      PORTIONS OF THE ROUTINE WILL BE
C                                      SKIPPED, BUT A SINGULAR MATRIX
C                                      CAN RESULT IF LARGE PORTIONS OF
C                                      THE REGION ARE WITHOUT DATA.
C
C                        NCF
C                          THE LENGTH OF THE ARRAY COEF IN THE CALLING
C                          PROGRAM.  IF NCF IS .LT.
C                          NODES(1)*...*NODES(NDIM) A FATAL ERROR IS
C                          DIAGNOSED.
C
C                        WORK
C                          A WORKSPACE ARRAY FOR SOLVING THE LEAST
C                          SQUARES MATRIX GENERATED BY THIS ROUTINE.
C                          ITS REQUIRED SIZE IS A FUNCTION OF THE TOTAL
C                          NUMBER OF NODES IN THE NODE GRID.  THIS
C                          TOTAL, NCOL = NODES(1)*...*NODES(NDIM), IS
C                          ALSO THE NUMBER OF COLUMNS IN THE LEAST
C                          SQUARES MATRIX.  THE LENGTH OF THE ARRAY WORK
C                          MUST EQUAL OR EXCEED NCOL*(NCOL+8)/2
C                          (OR NCOL*(NCOL+6)/2 IF XTRAP IS 0).
C                               THE DIMENSION IS ASSUMED TO BE
C                               WORK(NWRK).
C
C                        NWRK
C                          THE LENGTH OF THE ARRAY WORK IN THE CALLING
C                          PROGRAM.  IF
C                          NCOL = NODES(1)*...*NODES(NDIM) IS THE TOTAL
C                          NUMBER OF NODES, THEN A FATAL ERROR IS
C                          DIAGNOSED IF NWRK IS LESS THAN
C                          NCOL*(NCOL+8)/2 (OR NCOL*(NCOL+6)/2 IF
C                          XTRAP IS 0).
C
C ON OUTPUT              COEF
C                          THE ARRAY OF COEFFICIENTS COMPUTED BY THIS
C                          ROUTINE.  EACH COEFFICIENT CORRESPONDS TO A
C                          PARTICULAR BASIS FUNCTION WHICH IN TURN
C                          CORRESPONDS TO A NODE IN THE NODE GRID.  THIS
C                          CORRESPONDENCE BETWEEN THE NODE GRID AND THE
C                          ARRAY COEF IS AS IF COEF WERE AN
C                          NDIM-DIMENSIONAL FORTRAN ARRAY WITH
C                          DIMENSIONS NODES(1),...,NODES(NDIM), I.E., TO
C                          STORE THE ARRAY LINEARLY, THE LEFTMOST
C                          INDICES ARE INCREMENTED MOST FREQUENTLY.
C                          HENCE THE LENGTH OF THE COEF ARRAY MUST EQUAL
C                          OR EXCEED THE TOTAL NUMBER OF NODES, WHICH IS
C                          NODES(1)*...*NODES(NDIM).  THE COMPUTED ARRAY
C                          COEF MAY BE USED WITH FUNCTION SPLFE
C                          (OR SPLDE) TO EVALUATE THE SPLINE (OR ITS
C                          DERIVATIVES) AT AN ARBITRARY POINT IN NDIM
C                          SPACE.
C                               THE DIMENSION IS ASSUMED TO BE
C                               COEF(NCF).
C
C                        WORK
C                          THE WORKSPACE CONTAINING INTERMEDIATE
C                          CALCULATIONS.  IT NEED NOT BE SAVED.
C
C                        IERROR
C                          AN ERROR FLAG WITH THE FOLLOWING MEANINGS--
C                              0  NO ERROR.
C                            101  NDIM IS .LT. 1 OR IS .GT. 4.
C                            102  NODES(IDIM) IS .LT. 4 FOR SOME IDIM.
C                            103  XMIN(IDIM) = XMAX(IDIM) FOR SOME IDIM.
C                            104  NCF (SIZE OF COEF) IS
C                                 .LT. NODES(1)*...*NODES(NDIM).
C                            105  NDATA IS .LT. 1.
C                            106  NWRK (SIZE OF WORK) IS TOO
C                                 SMALL--SEE COMMENTS.
C                            107  SUPRLS FAILURE (USUALLY INSUFFICIENT
C                                 DATA) -- ORDINARILY OCCURS ONLY IF
C                                 XTRAP IS ZERO OR WDATA CONTAINS ALL
C                                 ZEROS.
C
C ALGORITHM              AN OVERDETERMINED SYSTEM OF LINEAR EQUATIONS
C                        IS FORMED--ONE EQUATION FOR EACH DATA POINT
C                        PLUS EQUATIONS FOR DERIVATIVE CONSTRAINTS.
C                        THIS SYSTEM IS SOLVED USING PACKAGE SUPRLS.
C
C ACCURACY               IF THERE IS EXACTLY ONE DATA POINT IN THE
C                        NEAR VICINITY OF EACH NODE AND NO EXTRA DATA,
C                        THE RESULTING SPLINE WILL AGREE WITH THE
C                        DATA VALUES TO 12 OR 13 DIGITS.  HOWEVER, IF
C                        THE PROBLEM IS OVERDETERMINED OR THE SPARSE
C                        DATA OPTION IS UTILIZED, THE ACCURACY IS HARD
C                        TO PREDICT.  BASICALLY, SMOOTH FUNCTIONS
C                        REQUIRE FEWER NODES THAN ROUGH ONES FOR THE
C                        SAME ACCURACY.
C
C TIMING                 THE EXECUTION TIME IS ROUGHLY PROPORTIONAL
C                        TO NDATA*NCOF**2 WHERE NCOF = NODES(1)*...*
C                        NODES(NDIM).
C***********************************************************************
      SUBROUTINE SPLCC (NDIM,XDATA,L1XDAT,YDATA,NDATA,XMIN,XMAX,NODES,
     1                  XTRAP,COEF,NCF,WORK,NWRK,IERROR)
      DIMENSION       XDATA(L1XDAT,NDATA)    ,YDATA(NDATA)           ,
     1                XMIN(NDIM) ,XMAX(NDIM) ,NODES(NDIM),COEF(NCF)  ,
     2                WORK(NWRK)
      DIMENSION       W(1)
      SAVE
C
      W(1) = -1.
      CALL SPLCW (NDIM,XDATA,L1XDAT,YDATA,W,NDATA,XMIN,XMAX,NODES,
     1            XTRAP,COEF,NCF,WORK,NWRK,IERROR)
      RETURN
      END
C PACKAGE SPLPAK         DOCUMENTATION FOR USER ENTRIES FOLLOWS
C                        THE GENERAL PACKAGE INFORMATION.
C
C LATEST REVISION        MARCH  1985
C
C PURPOSE                THIS PACKAGE CONTAINS ROUTINES FOR FITTING
C                        (LEAST SQUARES) A MULTIDIMENSIONAL CUBIC SPLINE
C                        TO ARBITRARILY LOCATED DATA.  IT ALSO CONTAINS
C                        ROUTINES FOR EVALUATING THIS SPLINE (OR ITS
C                        PARTIAL DERIVATIVES) AT ANY POINT.
C
C                        COEFFICIENT CALCULATION IS PERFORMED IN
C                        SUBROUTINES SPLCC OR SPLCW AND EVALUATION IS
C                        PERFORMED BY FUNCTIONS SPLFE OR SPLDE.
C
C USAGE                  PACKAGE SPLPAK CONTAINS FOUR USER ENTRIES --
C                        SPLCC, SPLCW, SPLFE, AND SPLDE.
C
C                        THE USER FIRST CALLS SPLCC BY
C
C                          CALL SPLCC (NDIM,XDATA,L1XDAT,YDATA,NDATA,
C                                      XMIN,XMAX,NODES,XTRAP,COEF,NCF,
C                                      WORK,NWRK,IERROR)
C
C                        OR SPLCW BY
C
C                          CALL SPLCW (NDIM,XDATA,L1XDATA,YDATA,WDATA,
C                                      NDATA,XMIN,XMAX,NODES,XTRAP,
C                                      COEF,NCF,WORK,NWRK,IERROR)
C
C                        THE PARAMETER NDATA IN THE CALL TO SPLCW
C                        ENABLES THE USER TO WEIGHT SONE OF THE DATA
C                        POINTS MORE HEAVILY THAN OTHERS.  BOTH
C                        ROUTINES RETURN A SET OF COEFFICIENTS IN THE
C                        ARRAY COEF.  THESE COEFFICIENTS ARE
C                        SUBSEQUENTLY USED IN THE COMPUTATION OF
C                        FUNCTION VALUES AND PARTIAL DERIVATIVES.
C                        THE USER THEN CALLS SPLFE OR SPLDE ANY
C                        NUMBER OF TIMES IN ANY ORDER PROVIDED THAT
C                        THE VALUES OF THE INPUTS, NDIM,COEF,XMIN,
C                        XMAX, AND NODES, ARE PRESERVED BETWEEN CALLS.
C                        THE ROUTINES ARE CALLED IN THE FOLLOWING WAY.
C
C                          F = SPLFE (NDIM,X,COEF,XMIN,XMAX,NODES,
C                                     IERROR)
C
C                        OR
C
C                          F = SPLDE (NDIM,X,NDERIV,COEF,XMIN,XMAX,
C                                     NODES,IERROR)
C
C                        THE ROUTINE, SPLFE,  RETURNS AN INTERPOLATED
C                        VALUE AT THE POINT DEFINED BY THE ARRAY X.
C                        XPLDE AFFORDS THE USER THE ADDITIONAL
C                        CAPABILITY OF CALCULATING AN INTERPOLATED
C                        VALUE FOR ONE OF SEVERAL PARTIAL DERIVATIVES
C                        SPECIFIED BY THE ARRAY NDERIV.
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE, EXCEPT FOR ERROR MESSAGES PRINTED BY
C                        CALLS TO ULIBER, IF AN ERROR IS DETECTED.
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       SUPRLS, ULIBER
C FILES
C
C LANGUAGE               FORTRAN
C
C HISTORY                DEVELOPED IN 1972-73 BY DAVE FULKER OF NCAR'S
C                        SCIENTIFIC COMPUTING DIVISION.
C
C PORTABILITY            FORTRAN 66
C***********************************************************************
C
C SUBROUTINE SPLCW (NDIM,XDATA,L1XDAT,YDATA,WDATA,NDATA,XMIN,XMAX,
C                   NODES,XTRAP,COEF,NCF,WORK,NWRK,IERROR)
C
C
C DIMENSION OF           XDATA(L1XDAT,NDATA),YDATA(NDATA),WDATA(NDATA),
C ARGUMENTS              XMIN(NDIM),XMAX(NDIM),NODES(NDIM),COEF(NCF),
C                        WORK(NWRK)
C
C PURPOSE                N-DIMENSIONAL CUBIC SPLINE COEFFICIENT
C                        CALCULATION BY WEIGHTED LEAST SQUARES ON
C                        ARBITRARILY LOCATED DATA.
C
C                        A GRID OF EVENLY SPACED NODES IN NDIM SPACE IS
C                        DEFINED BY THE ARGUMENTS XMIN, XMAX AND NODES.
C                        A LINEAR BASIS FOR THE CLASS OF NATURAL SPLINES
C                        ON THESE NODES IS FORMED, AND A SET OF
C                        CORRESPONDING COEFFICIENTS IS COMPUTED IN THE
C                        ARRAY COEF.  THESE COEFFICIENTS ARE CHOSEN TO
C                        MINIMIZE THE WEIGHTED SUM OF SQUARED ERRORS
C                        BETWEEN THE SPLINE AND THE ARBITRARILY LOCATED
C                        DATA VALUES DESCRIBED BY THE ARGUMENTS XDATA,
C                        YDATA AND NDATA.  THE SMOOTHNESS OF THE SPLINE
C                        IN DATA SPARSE AREAS IS CONTROLLED BY THE
C                        ARGUMENT XTRAP.
C
C NOTE                   IN ORDER TO UNDERSTAND THE ARGUMENTS OF THIS
C                        ROUTINE, ONE SHOULD REALIZE THAT THE NODE GRID
C                        NEED NOT BEAR ANY PARTICULAR RELATION TO THE
C                        DATA POINTS.  IN THE THEORY OF EXACT FIT
C                        INTERPOLATORY SPLINES THE NODES WOULD IN FACT
C                        BE DATA LOCATIONS, BUT IN THIS CASE THEY SERVE
C                        ONLY TO DEFINE THE CLASS OF SPLINES FROM WHICH
C                        THE APPROXIMATING FUNCTION IS CHOSEN.  THIS
C                        NODE GRID IS A RECTANGULAR ARRANGEMENT OF
C                        POINTS IN NDIM SPACE, WITH THE RESTRICTION THAT
C                        ALONG ANY COORDINATE DIRECTION THE NODES ARE
C                        EQUALLY SPACED.  THE CLASS OF NATURAL SPLINES
C                        ON THIS GRID OF NODES (NDIM-CUBIC SPLINES WHOSE
C                        2ND DERIVATIVES NORMAL TO THE BOUNDARIES ARE 0)
C                        HAS AS MANY DEGREES OF FREEDOM AS THE GRID HAS
C                        NODES.  THUS THE SMOOTHNESS OR FLEXIBILITY OF
C                        THE SPLINES IS DETERMINED BY THE CHOICE OF THE
C                        NODE GRID.
C
C USAGE                  CALL SPLCW (NDIM,XDATA,L1XDAT,YDATA,WDATA,
C                                    NDATA,XMIN,XMAX,NODES,XTRAP,COEF,
C                                    NCF,WORK,NWRK,IERROR)
C                        THE SPLINE (OR ITS DERIVATIVES) MAY THEN BE
C                        EVALUATED BY USING FUNCTION SPLFE (OR SPLDE).
C
C ARGUMENTS
C
C ON INPUT               NDIM
C                          THE DIMENSIONALITY OF THE PROBLEM.  THE
C                          SPLINE IS A FUNCTION OF NDIM VARIABLES OR
C                          COORDINATES AND THUS A POINT IN THE
C                          INDEPENDENT VARIABLE SPACE IS AN NDIM VECTOR.
C                          NDIM MUST BE IN THE RANGE 1 .LE. NDIM .LE. 4.
C
C                        XDATA
C                          A COLLECTION OF LOCATIONS FOR THE DATA
C                          VALUES, I.E., POINTS FROM THE INDEPENDENT
C                          VARIABLE SPACE.  THIS COLLECTION IS A
C                          2-DIMENSIONAL ARRAY WHOSE 1ST DIMENSION
C                          INDEXES THE NDIM COORDINATES OF A GIVEN POINT
C                          AND WHOSE 2ND DIMENSION LABELS THE DATA
C                          POINT.  FOR EXAMPLE, THE DATA POINT WITH
C                          LABEL IDATA IS LOCATED AT THE POINT
C                          (XDATA(1,IDATA),...,XDATA(NDIM,IDATA)) WHERE
C                          THE ELEMENTS OF THIS VECTOR ARE THE VALUES OF
C                          THE NDIM COORDINATES.  THE LOCATION, NUMBER
C                          AND ORDERING OF THE DATA POINTS IS ARBITRARY.
C                               THE DIMENSION IS ASSUMED TO BE
C                               XDATA(L1XDAT,NDATA).
C
C                        L1XDAT
C                          THE LENGTH OF THE 1ST DIMENSION OF XDATA IN
C                          THE CALLING PROGRAM.  L1XDAT MUST BE .GE.
C                          NDIM.
C                               NOTE:  FOR 1 DIMENSIONAL PROBLEMS L1XDAT
C                                      IS USUALLY 1.
C
C                        YDATA
C                          A COLLECTION OF DATA VALUES CORRESPONDING TO
C                          THE POINTS IN XDATA.  YDATA(IDATA) IS THE
C                          DATA VALUE ASSOCIATED WITH THE POINT
C                          (XDATA(1,IDATA),...,XDATA(NDIM,IDATA)) IN THE
C                          INDEPENDENT VARIABLE SPACE.  THE SPLINE WHOSE
C                          COEFFICIENTS ARE COMPUTED BY THIS ROUTINE
C                          APPROXIMATES THESE DATA VALUES IN THE LEAST
C                          SQUARES SENSE.
C                               THE DIMENSION IS ASSUMED TO BE
C                               YDATA(NDATA).
C
C                        WDATA
C                          A COLLECTION OF WEIGHTS.  WDATA(IDATA) IS A
C                          WEIGHT ASSOCIATED WITH THE DATA POINT
C                          LABELLED IDATA.  IT SHOULD BE NON-NEGATIVE,
C                          BUT MAY BE OF ANY MAGNITUDE.  THE WEIGHTS
C                          HAVE THE EFFECT OF FORCING GREATER OR LESSER
C                          ACCURACY AT A GIVEN POINT AS FOLLOWS--THIS
C                          ROUTINE CHOOSES COEFFICIENTS TO MINIMIZE THE
C                          SUM OVER ALL DATA POINTS OF THE QUANTITY
C                          (WDATA(IDATA)*(YDATA(IDATA) - SPLINE VALUE AT
C                          XDATA(IDATA)))**2.  THUS, IF THE RELIABILITY
C                          OF A DATA POINT IS KNOWN TO BE LOW, THE
C                          CORRESPONDING WEIGHT MAY BE MADE SMALL
C                          (RELATIVE TO THE OTHER WEIGHTS) SO THAT THE
C                          SUM OVER ALL DATA POINTS IS AFFECTED LESS BY
C                          DISCREPENCIES AT THE UNRELIABLE POINT.  DATA
C                          POINTS WITH ZERO WEIGHT ARE COMPLETELY
C                          IGNORED.
C                               NOTE:  IF WDATA(1) IS .LT. 0, THE OTHER
C                                      ELEMENTS OF WDATA ARE NOT
C                                      REFERENCED, AND ALL WEIGHTS ARE
C                                      ASSUMED TO BE UNITY.
C                          THE DIMENSION IS ASSUMED TO BE WDATA(NDATA)
C                          UNLESS WDATA(1) .LT. 0 IN WHICH CASE THE
C                          DIMENSION IS ASSUMED TO BE 1.
C
C                        NDATA
C                          THE NUMBER OF DATA POINTS MENTIONED IN THE
C                          ABOVE ARGUMENTS.
C
C                        XMIN
C                          A VECTOR DESCRIBING THE LOWER EXTREME CORNER
C                          OF THE NODE GRID.  A SET OF EVENLY SPACED
C                          NODES IS FORMED ALONG EACH COORDINATE AXIS
C                          AND XMIN(IDIM) IS THE LOCATION OF THE FIRST
C                          NODE ALONG THE IDIM AXIS.
C                               THE DIMENSION IS ASSUMED TO BE
C                               XMIN(NDIM).
C
C                        XMAX
C                          A VECTOR DESCRIBING THE UPPER EXTREME CORNER
C                          OF THE NODE GRID.  A SET OF EVENLY SPACED
C                          NODES IS FORMED ALONG EACH COORDINATE AXIS
C                          AND XMAX(IDIM) IS THE LOCATION OF THE LAST
C                          NODE ALONG THE IDIM AXIS.
C                               THE DIMENSION IS ASSUMED TO BE
C                               XMAX(NDIM).
C
C                        NODES
C                          A VECTOR OF INTEGERS DESCRIBING THE NUMBER OF
C                          NODES ALONG EACH AXIS.  NODES(IDIM) IS THE
C                          NUMBER OF NODES (COUNTING ENDPOINTS) ALONG
C                          THE IDIM AXIS AND DETERMINES THE FLEXIBILITY
C                          OF THE SPLINE IN THAT COORDINATE DIRECTION.
C                          NODES(IDIM) MUST BE .GE. 4 BUT MAY BE AS
C                          LARGE AS THE ARRAYS COEF AND WORK ALLOW.
C                               THE DIMENSION IS ASSUMED TO BE
C                               NODES(NDIM).
C
C                          NOTE:  THE NODE GRID IS COMPLETELY DEFINED BY
C                                 THE ARGUMENTS XMIN, XMAX AND NODES.
C                                 THE SPACING OF THIS GRID IN THE IDIM
C                                 COORDINATE DIRECTION IS
C                                 DX(IDIM) = (XMAX(IDIM)-XMIN(IDIM)) /
C                                             (NODES(IDIM)-1).
C                                 A NODE IN THIS GRID MAY BE INDEXED BY
C                                 AN NDIM VECTOR OF INTEGERS
C                                 (IN(1),...,IN(NDIM)) WHERE
C                                 1 .LE. IN(IDIM) .LE. NODES(IDIM).
C                                 THE LOCATION OF SUCH A NODE MAY BE
C                                 REPRESENTED BY AN NDIM VECTOR
C                                 (X(1),...,X(NDIM)) WHERE
C                                 X(IDIM) = XMIN(IDIM) + (IN(IDIM)-1) *
C                                                         DX(IDIM).
C
C                        XTRAP
C                          A PARAMETER TO CONTROL EXTRAPOLATION TO DATA
C                          SPARSE AREAS.  THE REGION DESCRIBED BY XMIN
C                          AND XMAX IS DIVIDED INTO RECTANGLES, THE
C                          NUMBER OF WHICH IS DETERMINED BY NODES, AND
C                          ANY RECTANGLE CONTAINING A DISPROPORTIONATELY
C                          SMALL NUMBER OF DATA POINTS IS CONSIDERED TO
C                          BE DATA SPARSE (RECTANGLE IS USED HERE TO
C                          MEAN NDIM-DIMENSIONAL RECTANGLE).  IF XTRAP
C                          IS NONZERO THE LEAST SQUARES PROBLEM IS
C                          ARGUMENTED WITH DERIVATIVE CONSTRAINTS IN THE
C                          DATA SPARSE AREAS TO PREVENT THE MATRIX FROM
C                          BECOMING POORLY CONDITIONED.  XTRAP SERVES AS
C                          A WEIGHT FOR THESE CONSTRAINTS, AND THUS MAY
C                          BE USED TO CONTROL SMOOTHNESS IN DATA SPARSE
C                          AREAS.  THE EXPERIENCE OF THE AUTHOR
C                          INDICATES THAT UNITY IS A GOOD FIRST GUESS
C                          FOR THIS PARAMETER, AND THAT ITS SIZE IS NOT
C                          VERY CRITICAL.
C
C                               NOTE:  IF XTRAP IS 0 SUBSTANTIAL
C                                      PORTIONS OF THE ROUTINE WILL BE
C                                      SKIPPED, BUT A SINGULAR MATRIX
C                                      CAN RESULT IF LARGE PORTIONS OF
C                                      THE REGION ARE WITHOUT DATA.
C
C                        NCF
C                          THE LENGTH OF THE ARRAY COEF IN THE CALLING
C                          PROGRAM.  IF NCF IS .LT.
C                          NODES(1)*...*NODES(NDIM) A FATAL ERROR IS
C                          DIAGNOSED.
C
C                        WORK
C                          A WORKSPACE ARRAY FOR SOLVING THE LEAST
C                          SQUARES MATRIX GENERATED BY THIS ROUTINE.
C                          ITS REQUIRED SIZE IS A FUNCTION OF THE TOTAL
C                          NUMBER OF NODES IN THE NODE GRID.  THIS
C                          TOTAL, NCOL = NODES(1)*...*NODES(NDIM), IS
C                          ALSO THE NUMBER OF COLUMNS IN THE LEAST
C                          SQUARES MATRIX.  THE LENGTH OF THE ARRAY WORK
C                          MUST EQUAL OR EXCEED NCOL*(NCOL+8)/2
C                          (OR NCOL*(NCOL+6)/2 IF XTRAP IS 0).
C                               THE DIMENSION IS ASSUMED TO BE
C                               WORK(NWRK).
C
C                        NWRK
C                          THE LENGTH OF THE ARRAY WORK IN THE CALLING
C                          PROGRAM.  IF
C                          NCOL = NODES(1)*...*NODES(NDIM) IS THE TOTAL
C                          NUMBER OF NODES, THEN A FATAL ERROR IS
C                          DIAGNOSED IF NWRK IS LESS THAN
C                          NCOL*(NCOL+8)/2 (OR NCOL*(NCOL+6)/2 IF
C                          XTRAP IS 0).
C
C ON OUTPUT              COEF
C                          THE ARRAY OF COEFFICIENTS COMPUTED BY THIS
C                          ROUTINE.  EACH COEFFICIENT CORRESPONDS TO A
C                          PARTICULAR BASIS FUNCTION WHICH IN TURN
C                          CORRESPONDS TO A NODE IN THE NODE GRID.  THIS
C                          CORRESPONDENCE BETWEEN THE NODE GRID AND THE
C                          ARRAY COEF IS AS IF COEF WERE AN
C                          NDIM-DIMENSIONAL FORTRAN ARRAY WITH
C                          DIMENSIONS NODES(1),...,NODES(NDIM), I.E., TO
C                          STORE THE ARRAY LINEARLY, THE LEFTMOST
C                          INDICES ARE INCREMENTED MOST FREQUENTLY.
C                          HENCE THE LENGTH OF THE COEF ARRAY MUST EQUAL
C                          OR EXCEED THE TOTAL NUMBER OF NODES, WHICH IS
C                          NODES(1)*...*NODES(NDIM).  THE COMPUTED ARRAY
C                          COEF MAY BE USED WITH FUNCTION SPLFE
C                          (OR SPLDE) TO EVALUATE THE SPLINE (OR ITS
C                          DERIVATIVES) AT AN ARBITRARY POINT IN NDIM
C                          SPACE.
C                               THE DIMENSION IS ASSUMED TO BE
C                               COEF(NCF).
C
C                        WORK
C                          THE WORKSPACE CONTAINING INTERMEDIATE
C                          CALCULATIONS.  IT NEED NOT BE SAVED.
C
C                        IERROR
C                          AN ERROR FLAG WITH THE FOLLOWING MEANINGS--
C                              0  NO ERROR.
C                            101  NDIM IS .LT. 1 OR IS .GT. 4.
C                            102  NODES(IDIM) IS .LT. 4 FOR SOME IDIM.
C                            103  XMIN(IDIM) = XMAX(IDIM) FOR SOME IDIM.
C                            104  NCF (SIZE OF COEF) IS
C                                 .LT. NODES(1)*...*NODES(NDIM).
C                            105  NDATA IS .LT. 1.
C                            106  NWRK (SIZE OF WORK) IS TOO
C                                 SMALL--SEE COMMENTS.
C                            107  SUPRLS FAILURE (USUALLY INSUFFICIENT
C                                 DATA) -- ORDINARILY OCCURS ONLY IF
C                                 XTRAP IS ZERO OR WDATA CONTAINS ALL
C                                 ZEROS.
C
C ALGORITHM              AN OVERDETERMINED SYSTEM OF LINEAR EQUATIONS
C                        IS FORMED--ONE EQUATION FOR EACH DATA POINT
C                        PLUS EQUATIONS FOR DERIVATIVE CONSTRAINTS.
C                        THIS SYSTEM IS SOLVED USING PACKAGE SUPRLS.
C
C ACCURACY               IF THERE IS EXACTLY ONE DATA POINT IN THE
C                        NEAR VICINITY OF EACH NODE AND NO EXTRA DATA,
C                        THE RESULTING SPLINE WILL AGREE WITH THE
C                        DATA VALUES TO 12 OR 13 DIGITS.  HOWEVER, IF
C                        THE PROBLEM IS OVERDETERMINED OR THE SPARSE
C                        DATA OPTION IS UTILIZED, THE ACCURACY IS HARD
C                        TO PREDICT.  BASICALLY, SMOOTH FUNCTIONS
C                        REQUIRE FEWER NODES THAN ROUGH ONES FOR THE
C                        SAME ACCURACY.
C
C TIMING                 THE EXECUTION TIME IS ROUGHLY PROPORTIONAL
C                        TO NDATA*NCOF**2 WHERE NCOF = NODES(1)*...*
C                        NODES(NDIM).
C***********************************************************************
      SUBROUTINE SPLCW (NDIM,XDATA,L1XDAT,YDATA,WDATA,NDATA,XMIN,XMAX,
     1                  NODES,XTRAP,COEF,NCF,WORK,NWRK,IERROR)
      DIMENSION       XDATA(L1XDAT,NDATA)    ,YDATA(NDATA)           ,
     1                WDATA(NDATA)           ,XMIN(NDIM) ,XMAX(NDIM) ,
     2                NODES(NDIM),COEF(NCF)  ,WORK(NWRK)
      DIMENSION       X(4)       ,NDERIV(4)  ,IN(4)      ,INMX(4)
      COMMON /SPLCOM/ MDIM       ,DX(4)      ,DXIN(4)    ,IB(4)      ,
     1                IBMN(4)    ,IBMX(4)
      SAVE
C
C THE RESTRICTION FOR NDIM TO BE LE 4 CAN BE ELIMINATED BY INCREASING
C THE ABOVE DIMENSIONS, BUT THE REQUIRED LENGTH OF WORK BECOMES
C PROHIBITIVE ON CURRENT COMPUTERS.
C
C SPCRIT IS USED TO DETERMINE DATA SPARSENESS AS FOLLOWS -
C THE WEIGHTS ASSIGNED TO ALL DATA POINTS ARE TOTALED INTO THE
C VARIABLE TOTLWT. (IF NO WEIGHTS ARE ENTERED IT IS SET TO
C NDATA.)  EACH NODE OF THE NODE NETWORK IS ASSIGNED A
C RECTANGLE (IN WHICH IT IS CONTAINED) AND THE WEIGHTS OF ALL
C DATA POINTS WHICH FALL IN THAT RECTANGLE ARE TOTALED.  IF THAT
C TOTAL IS LESS THAN SPCRIT*EXPECT (EXPECT IS DEFINED BELOW)
C THEN THE NODE IS ASCERTAINED TO BE IN A DATA SPARSE LOCATION.
C EXPECT IS THAT FRACTION OF TOTLWT THAT WOULD BE EXPECTED BY
C COMPARING THE AREA OF THE RECTANGLE WITH THE TOTAL AREA UNDER
C CONSIDERATION.
C
      DATA SPCRIT/.75/
C
      IERROR = 0
      MDIM = NDIM
      IF (MDIM.LT.1 .OR. MDIM.GT.4) GO TO 127
      NCOL = 1
      DO 101 IDIM=1,MDIM
         NOD = NODES(IDIM)
         IF (NOD .LT. 4) GO TO 128
C
C NUMBER OF COLUMNS IN LEAST SQUARES MATRIX = NUMBER OF COEFFICIENTS =
C PRODUCT OF NODES OVER ALL DIMENSIONS
C
         NCOL = NCOL*NOD
         XRNG = XMAX(IDIM)-XMIN(IDIM)
         IF (XRNG .EQ. 0.) GO TO 129
C
C DX(IDIM) IS THE NODE SPACING ALONG THE IDIM COORDINATE
C
         DX(IDIM) = XRNG/FLOAT(NOD-1)
         DXIN(IDIM) = 1./DX(IDIM)
         NDERIV(IDIM) = 0
  101 CONTINUE
      IF (NCOL .GT. NCF) GO TO 130
      NWRK1 = 1
      MDATA = NDATA
      IF (MDATA .LT. 1) GO TO 131
C
C SWGHT IS A LOCAL VARIABLE = XTRAP, AND CAN BE CONSIDERED A SMOOTHING
C WEIGHT FOR DATA SPARSE AREAS.  IF SWGHT .EQ. 0 NO SMOOTHING
C COMPUTATIONS ARE PERFORMED
C
      SWGHT = XTRAP
C
C SET ASIDE WORKSPACE FOR COUNTING DATA POINTS
C
      IF (SWGHT .NE. 0.) NWRK1 = NCOL+1
C
C NWLFT IS THE LENGTH OF THE REMAINING WORKSPACE
C
      NWLFT = NWRK-NWRK1+1
      IF (NWLFT .LT. 1) GO TO 132
      IROW = 0
C
C ROWWT IS USED TO WEIGHT ROWS OF THE LEAST SQUARES MATRIX
C
      ROWWT = 1.
C
C LOOP THRU ALL DATA POINTS, COMPUTING A ROW FOR EACH
C
      DO 108 IDATA=1,MDATA
C
C WDATA(1) LT 0 MEANS WEIGHTS HAVE NOT BEEN ENTERED.  IN THAT CASE
C ROWWT IS LEFT EQUAL TO  1. FOR ALL POINTS.  OTHERWISE ROWWT IS
C EQUAL TO WDATA(IDATA).
C EVERY ELEMENT OF THE ROW, AS WELL AS THE CORRESPONDING RIGHT HAND
C SIDE, IS MULTIPLIED BY ROWWT
C
         IF (WDATA(1) .LT. 0.) GO TO 102
         ROWWT = WDATA(IDATA)
C
C DATA POINTS WITH 0 WEIGHT ARE IGNORED
C
         IF (ROWWT .EQ. 0.) GO TO 108
  102    IROW = IROW+1
C
C ONE ROW OF THE LEAST SQUARES MATRIX CORRESPONDS TO EACH DATA POINT
C THE RIGHT HAND FOR THAT ROW WILL CORRESPOND TO THE FUNCTION VALUE
C YDATA AT THAT POINT
C
         RHS = ROWWT*YDATA(IDATA)
         DO 103 IDIM=1,MDIM
            X(IDIM) = XDATA(IDIM,IDATA)
  103    CONTINUE
C
C COEF ARRAY SERVES AS ROW OF LEAST SQUARES MATRIX.  ITS VALUE IS ZERO
C EXCEPT FOR COLUMNS CORRESPONDING TO FUNCTIONS WHICH ARE
C NONZERO AT X.
C
         DO 104 ICOL=1,NCOL
            COEF(ICOL) = 0.
  104    CONTINUE
         DO 105 IDIM=1,MDIM
            NOD = NODES(IDIM)
C
C COMPUTE INDICES OF BASIS FUNCTIONS WHICH ARE NONZERO AT X
C
            IT = DXIN(IDIM)*(X(IDIM)-XMIN(IDIM))
C
C IBMN MUST BE IN RANGE 0 TO NODES-2
C
            IBMN(IDIM) = MIN0(MAX0(IT-1,0),NOD-2)
C
C INITIALIZE BASIS INDICES
C
            IB(IDIM) = IBMN(IDIM)
C
C IBMX MUST BE IN RANGE 1 TO NODES-1
C
            IBMX(IDIM) = MAX0(MIN0(IT+2,NOD-1),1)
  105    CONTINUE
C
C BEGINING OF BASIS INDEX LOOP- TRAVERSE ALL INDICES CORRESPONDING
C TO BASIS FUNCTIONS WHICH ARE NONZERO AT X
C THE INDICES ARE IN IB AND PASS THRU COMMON TO BASCMP
C
  106    CALL BASCMP (X,NDERIV,XMIN,NODES,ICOL,BASM)
C
C BASCMP COMPUTES ICOL AND BASM WHERE BASM IS THE VALUE (AT X) OF THE
C N-DIMENSIONAL BASIS FUNCTION CORRESPONDING TO COLUMN ICOL
C
         COEF(ICOL) = ROWWT*BASM
C
C INCREMENT BASIS INDICES
C
         DO 107 IDIM=1,MDIM
            IB(IDIM) = IB(IDIM)+1
            IF (IB(IDIM) .LE. IBMX(IDIM)) GO TO 106
            IB(IDIM) = IBMN(IDIM)
  107    CONTINUE
C
C END OF BASIS INDEX LOOP
C
C SEND ROW OF LEAST SQUARES MATRIX TO REDUCTION ROUTINE
C
         CALL SUPRLS (IROW,COEF,NCOL,RHS,WORK(NWRK1),NWLFT,COEF,RESERR,
     1                LSERR)
         IF (LSERR .NE. 0) GO TO 133
  108 CONTINUE
C
C ROW COMPUTATIONS FOR DATA POINTS COMPLETE
C
C IF SWGHT EQ 0 THE LEAST SQUARES MATRIX IS COMPLETE AND NO SMOOTHING
C ROWS ARE COMPUTED
C
      IF (SWGHT .EQ. 0.) GO TO 126
C
C INITIALIZE SMOOTHING COMPUTATIONS FOR DATA SPARSE AREAS
C DERIVATIVE CONSTRAINTS WILL ALWAYS HAVE ZERO RIGHT HAND SIDE
C
      RHS = 0.
      NRECT = 1
      DO 109 IDIM=1,MDIM
C
C INITIALIZE NODE INDICES
C
         IN(IDIM) = 0
         INMX(IDIM) = NODES(IDIM)-1
C
C COMPUTE NUMBER OF RECTANGLES FORMED BY THE NODE NETWORK
C
         NRECT = NRECT*INMX(IDIM)
  109 CONTINUE
C
C EVERY NODE IS ASSIGNED AN ELEMENT OF THE WORKSPACE (SET ASIDE
C PREVIOUSLY) IN WHICH DATA POINTS ARE COUNTED
C
      DO 110 IIN=1,NCOL
         WORK(IIN) = 0.
  110 CONTINUE
C
C ASSIGN EACH DATA POINT TO A NODE, TOTAL THE ASSIGNMENTS FOR EACH
C NODE, AND SAVE IN THE WORKSPACE
C
      TOTLWT = 0.
      DO 112 IDATA=1,MDATA
C
C BUMP IS THE WEIGHT ASSOCIATED WITH THE DATA POINT
C
         BUMP = 1.
         IF (WDATA(1) .GE. 0.) BUMP = WDATA(IDATA)
         IF (BUMP .EQ. 0.) GO TO 112
C
C FIND THE NEAREST NODE
C
         IIN = 0
         DO 111 IDIMC=1,MDIM
            IDIM = MDIM+1-IDIMC
            INIDIM = INT(DXIN(IDIM)*(XDATA(IDIM,IDATA)-XMIN(IDIM))+.5)
C
C POINTS NOT IN RANGE (+ OR - 1/2 NODE SPACING) ARE NOT COUNTED
C
            IF (INIDIM.LT.0 .OR. INIDIM.GT.INMX(IDIM)) GO TO 112
C
C COMPUTE LINEAR ADDRESS OF NODE IN WORKSPACE BY HORNERS METHOD
C
            IIN = (INMX(IDIM)+1)*IIN+INIDIM
  111    CONTINUE
C
C BUMP COUNTER FOR THAT NODE
C
         WORK(IIN+1) = WORK(IIN+1)+BUMP
         TOTLWT = TOTLWT+BUMP
  112 CONTINUE
C
C COMPUTE EXPECTED WEIGHT PER RECTANGLE
C
      WTPRRC = TOTLWT/FLOAT(NRECT)
C
C IN CONTAINS INDICES OF THE NODE (PREVIOUSLY INITIALIZED)
C IIN WILL BE LINEAR ADDRESS OF NODE IN WORKSPACE
C
      IIN = 0
C
C LOOP THRU ALL NODES, COMPUTING DERIVATIVE CONSTRAINT ROWS FOR THOSE
C IN DATA SPARSE LOCATIONS
C
C BEGINING OF NODE INDEX LOOP- TRAVERSE ALL NODE INDICES
C THE INDICES ARE IN IN
C
  113 IIN = IIN+1
      EXPECT = WTPRRC
      DO 114 IDIM=1,MDIM
C
C RECTANGLES AT EDGE OF NETWORK ARE SMALLER AND HENCE LESS WEIGHT
C SHOULD BE EXPECTED
C
         IF (IN(IDIM).EQ.0 .OR. IN(IDIM).EQ.INMX(IDIM))
     1       EXPECT = .5*EXPECT
  114 CONTINUE
C
C THE EXPECTED WEIGHT MINUS THE ACTUAL WEIGHT SERVES TO DEFINE DATA
C SPARSENESS AND IS ALSO USED TO WEIGHT THE DERIVATIVE
C CONSTRAINT ROWS
C NO CONSTRAINT IF NOT DATA SPARCE
C
      IF (WORK(IIN) .GE. SPCRIT*EXPECT) GO TO 124
      DCWGHT = EXPECT-WORK(IIN)
      DO 115 IDIM=1,MDIM
         INIDIM = IN(IDIM)
C
C COMPUTE LOCATION OF NODE
C
         X(IDIM) = XMIN(IDIM)+FLOAT(INIDIM)*DX(IDIM)
C
C COMPUTE INDICES OF BASIS FUNCTIONS WHICH ARE NONZERO AT THE NODE
C
         IBMN(IDIM) = INIDIM-1
         IBMX(IDIM) = INIDIM+1
C
C DISTINGUISH BOUNDARIES
C
         IF (INIDIM .EQ. 0) IBMN(IDIM) = 0
         IF (INIDIM .EQ. INMX(IDIM)) IBMX(IDIM) = INMX(IDIM)
C
C INITIALIZE BASIS INDICES
C
         IB(IDIM) = IBMN(IDIM)
  115 CONTINUE
C
C MULTIPLY BY EXTRAPOLATION PARAMETER (ACTS AS SMOOTHING WEIGHT)
C
      DCWGHT = SWGHT*DCWGHT
C
C COEF ARRAY SERVES AS ROW OF LEAST SQUARES MATRIX.  ITS VALUE IS ZERO
C EXCEPT FOR COLUMNS CORRESPONDING TO FUNCTIONS WHICH ARE
C NONZERO AT THE NODE
C
      DO 116 ICOL=1,NCOL
         COEF(ICOL) = 0.
  116 CONTINUE
C
C 2ND DERIVATIVE FOR FUNCTION OF MDIM VARIABLES MAY BE THOUGHT OF AS A
C SYMMETRIC MDIM BY MDIM MATRIX OF 2ND ORDER PARTIAL DERIVATIVES
C TRAVERSE THE UPPER TRIANGLE OF THIS MATRIX AND, FOR EACH ELEMENT,
C COMPUTE A ROW OF THE LEAST SQUARES MATRIX
C
      DO 123 IDM=1,MDIM
         DO 122 JDM=IDM,MDIM
            DO 117 IDIM=1,MDIM
               NDERIV(IDIM) = 0
  117       CONTINUE
C
C OFF DIAGONAL ELEMENTS APPEAR TWICE BY SYMMETRY SO THE CORRESPONDING
C ROW IS WEIGHTED BY A FACTOR OF 2
C
            ROWWT = 2.*DCWGHT
            IF (JDM .NE. IDM) GO TO 118
C
C DIAGONAL
C
            ROWWT = DCWGHT
            NDERIV(JDM) = 2
            IF (IN(IDM).NE.0 .AND. IN(IDM).NE.INMX(IDM)) GO TO 119
C
C NODE IS AT BOUNDARY
C NORMAL 2ND DERIVATIVE CONSTRAINT AT BOUNDARY IS NOT APPROPRIATE FOR
C NATURAL SPLINES (2ND DERIVATIVE 0 BY DEFINITION).  SUBSTITUTE
C 1ST DERIVATIVE CONSTRAINT.
C
  118       NDERIV(IDM) = 1
            NDERIV(JDM) = 1
  119       IROW = IROW+1
C
C BEGINING OF BASIS INDEX LOOP- TRAVERSE ALL INDICES CORRESPONDING
C TO BASIS FUNCTIONS WHICH ARE NONZERO AT X
C THE INDICES ARE IN IB AND PASS THRU COMMON TO BASCMP
C
  120       CALL BASCMP (X,NDERIV,XMIN,NODES,ICOL,BASM)
C
C BASCMP COMPUTES ICOL AND BASM WHERE BASM IS THE VALUE (AT X) OF THE
C N-DIMENSIONAL BASIS FUNCTION CORRESPONDING TO COLUMN ICOL
C
            COEF(ICOL) = ROWWT*BASM
C
C INCREMENT BASIS INDICES
C
            DO 121 IDIM=1,MDIM
               IB(IDIM) = IB(IDIM)+1
               IF (IB(IDIM) .LE. IBMX(IDIM)) GO TO 120
               IB(IDIM) = IBMN(IDIM)
  121       CONTINUE
C
C END OF BASIS INDEX LOOP
C
C SEND ROW OF LEAST SQUARES MATRIX TO REDUCTION ROUTINE
C
         CALL SUPRLS (IROW,COEF,NCOL,RHS,WORK(NWRK1),NWLFT,COEF,RESERR,
     1                LSERR)
         IF (LSERR .NE. 0) GO TO 133
  122    CONTINUE
  123 CONTINUE
C
C INCREMENT NODE INDICES
C
  124 DO 125 IDIM=1,MDIM
         IN(IDIM) = IN(IDIM)+1
         IF (IN(IDIM) .LE. INMX(IDIM)) GO TO 113
         IN(IDIM) = 0
  125 CONTINUE
C
C END OF NODE INDEX LOOP
C
C CALL FOR LEAST SQUARES SOLUTION IN COEF ARRAY
C
  126 IROW = 0
         CALL SUPRLS (IROW,COEF,NCOL,RHS,WORK(NWRK1),NWLFT,COEF,RESERR,
     1                LSERR)
         IF (LSERR .NE. 0) GO TO 133
      RETURN
C
C ERROR SECTION
C
  127 CONTINUE
      IERROR = 101
      CALL ULIBER (IERROR,60H SPLCC OR SPLCW - NDIM IS LT 1 OR IS GT 4
     1                 ,60)
      GO TO 134
  128 CONTINUE
      IERROR = 102
      CALL ULIBER (IERROR,60H SPLCC OR SPLCW - NODES(IDIM) IS LT 4 FOR S
     1OME IDIM         ,60)
      GO TO 134
  129 CONTINUE
      IERROR = 103
      CALL ULIBER (IERROR,60H SPLCC OR SPLCW - XMIN(IDIM) EQ XMAX(IDIM)
     1FOR SOME IDIM    ,60)
      GO TO 134
  130 CONTINUE
      IERROR = 104
      CALL ULIBER (IERROR,60H SPLCC OR SPLCW - NCF (SIZE OF COEF) IS TOO
     1 SMALL           ,60)
      GO TO 134
  131 CONTINUE
      IERROR = 105
      CALL ULIBER (IERROR,60H SPLCC OR SPLCW - NDATA IS LT 1
     1                 ,60)
      GO TO 134
  132 CONTINUE
      IERROR = 106
      CALL ULIBER (IERROR,60H SPLCC OR SPLCW - NWRK (SIZE OF WORK) IS TO
     1O SMALL          ,60)
      GO TO 134
  133 CONTINUE
      IERROR = 107
      CALL ULIBER (IERROR,60H SPLCC OR SPLCW - SUPRLS FAILURE (USUALLY I
     1NSUFFICIENT DATA),60)
  134 RETURN
      END
C PACKAGE SPLPAK         DOCUMENTATION FOR USER ENTRIES FOLLOWS
C                        THE GENERAL PACKAGE INFORMATION.
C
C LATEST REVISION        MARCH  1985
C
C PURPOSE                THIS PACKAGE CONTAINS ROUTINES FOR FITTING
C                        (LEAST SQUARES) A MULTIDIMENSIONAL CUBIC SPLINE
C                        TO ARBITRARILY LOCATED DATA.  IT ALSO CONTAINS
C                        ROUTINES FOR EVALUATING THIS SPLINE (OR ITS
C                        PARTIAL DERIVATIVES) AT ANY POINT.
C
C                        COEFFICIENT CALCULATION IS PERFORMED IN
C                        SUBROUTINES SPLCC OR SPLCW AND EVALUATION IS
C                        PERFORMED BY FUNCTIONS SPLFE OR SPLDE.
C
C USAGE                  PACKAGE SPLPAK CONTAINS FOUR USER ENTRIES --
C                        SPLCC, SPLCW, SPLFE, AND SPLDE.
C
C                        THE USER FIRST CALLS SPLCC BY
C
C                          CALL SPLCC (NDIM,XDATA,L1XDAT,YDATA,NDATA,
C                                      XMIN,XMAX,NODES,XTRAP,COEF,NCF,
C                                      WORK,NWRK,IERROR)
C
C                        OR SPLCW BY
C
C                          CALL SPLCW (NDIM,XDATA,L1XDATA,YDATA,WDATA,
C                                      NDATA,XMIN,XMAX,NODES,XTRAP,
C                                      COEF,NCF,WORK,NWRK,IERROR)
C
C                        THE PARAMETER NDATA IN THE CALL TO SPLCW
C                        ENABLES THE USER TO WEIGHT SONE OF THE DATA
C                        POINTS MORE HEAVILY THAN OTHERS.  BOTH
C                        ROUTINES RETURN A SET OF COEFFICIENTS IN THE
C                        ARRAY COEF.  THESE COEFFICIENTS ARE
C                        SUBSEQUENTLY USED IN THE COMPUTATION OF
C                        FUNCTION VALUES AND PARTIAL DERIVATIVES.
C                        THE USER THEN CALLS SPLFE OR SPLDE ANY
C                        NUMBER OF TIMES IN ANY ORDER PROVIDED THAT
C                        THE VALUES OF THE INPUTS, NDIM,COEF,XMIN,
C                        XMAX, AND NODES, ARE PRESERVED BETWEEN CALLS.
C                        THE ROUTINES ARE CALLED IN THE FOLLOWING WAY.
C
C                          F = SPLFE (NDIM,X,COEF,XMIN,XMAX,NODES,
C                                     IERROR)
C
C                        OR
C
C                          F = SPLDE (NDIM,X,NDERIV,COEF,XMIN,XMAX,
C                                     NODES,IERROR)
C
C                        THE ROUTINE, SPLFE,  RETURNS AN INTERPOLATED
C                        VALUE AT THE POINT DEFINED BY THE ARRAY X.
C                        XPLDE AFFORDS THE USER THE ADDITIONAL
C                        CAPABILITY OF CALCULATING AN INTERPOLATED
C                        VALUE FOR ONE OF SEVERAL PARTIAL DERIVATIVES
C                        SPECIFIED BY THE ARRAY NDERIV.
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE, EXCEPT FOR ERROR MESSAGES PRINTED BY
C                        CALLS TO ULIBER, IF AN ERROR IS DETECTED.
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       SUPRLS, ULIBER
C FILES
C
C LANGUAGE               FORTRAN
C
C HISTORY                DEVELOPED IN 1972-73 BY DAVE FULKER OF NCAR'S
C                        SCIENTIFIC COMPUTING DIVISION.
C
C PORTABILITY            FORTRAN 66
C***********************************************************************
C FUNCTION SPLDE (NDIM,X,NDERIV,COEF,XMIN,XMAX,NODES,IERROR)
C
C DIMENSION OF           X(NDIM),NDERIV(NDIM),
C ARGUMENTS              COEF(NODES(1)*...*NODES(NDIM)),XMIN(NDIM),
C                        XMAX(NDIM),NODES(NDIM)
C
C PURPOSE                N-DIMENSIONAL CUBIC SPLINE DERIVATIVE
C                        EVALUATION.
C
C                        A GRID OF EVENLY SPACED NODES IN NDIM SPACE IS
C                        DEFINED BY THE ARGUMENTS XMIN, XMAX AND NODES.
C                        A LINEAR BASIS FOR THE CLASS OF NATURAL SPLINES
C                        ON THESE NODES IS FORMED, AND TO EACH BASIS
C                        FUNCTION CORRESPONDS A COEFFICIENT IN THE ARRAY
C                        COEF (COMPUTED IN SPLCC OR SPLCW).  USING
C                        NDERIV TO INDICATE THE APPROPRIATE PARTIAL
C                        DERIVATIVES, EACH BASIS FUNCTION IS EVALUATED
C                        AT THE POINT X IN NDIM SPACE.  THESE VALUES ARE
C                        THEN MULTIPLIED BY THE CORRESPONDING
C                        COEFFICIENT AND SUMMED TO FORM THE FUNCTION
C                        RESULT.
C
C USAGE                  COEF, XMIN, XMAX AND NODES MUST BE EXACTLY
C                        RETAINED FROM THE CALL TO SPLCC (OR SPLCW).
C
C                        THEN EVALUATE PARTIAL DERIVATIVES OF THE SPLINE
C                        SET BY
C
C                        Y = SPLDE (NDIM,X,NDERIV,COEF,XMIN,XMAX,NODES,
C                                   IERROR).
C
C ARGUMENTS
C
C ON INPUT               NDIM
C                          THE DIMENSIONALITY OF THE PROBLEM.  THE
C                          SPLINE IS A FUNCTION OF NDIM VARIABLES OR
C                          COORDINATES AND THUS A POINT IN THE
C                          INDEPENDENT VARIABLE SPACE IS AN NDIM VECTOR.
C                          NDIM MUST BE IN THE RANGE 1 .LE. NDIM .LE. 4.
C
C                        X
C                          AN NDIM VECTOR DESCRIBING THE POINT IN THE
C                          INDEPENDENT VARIABLE SPACE AT WHICH THE
C                          SPLINE IS TO BE EVALUATED.
C                               THE DIMENSION IS ASSUMED TO BE X(NDIM).
C
C                        NDERIV
C                          AN NDIM VECTOR OF INTEGERS SPECIFYING THE
C                          PARTIAL DERIVATIVE TO BE EVALUATED.  THE
C                          ORDER OF THE DERIVATIVE ALONG THE IDIM AXIS
C                          IS NDERIV(IDIM).  THESE INTEGERS MUST BE IN
C                          THE RANGE 0 .LE. NDERIV(IDIM) .LE. 2.
C
C                        COEF
C                          THE ARRAY OF COEFFICIENTS WHICH DETERMINE THE
C                          SPLINE.  EACH COEFFICIENT CORRESPONDS TO A
C                          PARTICULAR BASIS FUNCTION WHICH IN TURN
C                          CORRESPONDS TO A NODE IN THE NODE GRID.  THIS
C                          CORRESPONDENCE BETWEEN THE NODE GRID AND THE
C                          ARRAY COEF IS AS IF COEF WERE AN
C                          NDIM-DIMENSIONAL FORTRAN ARRAY WITH
C                          DIMENSIONS NODES(1),...,NODES(NDIM), I.E., TO
C                          STORE THE ARRAY LINEARLY, THE LEFTMOST
C                          INDICES ARE INCREMENTED MOST FREQUENTLY.
C                          COEF MAY BE COMPUTED BY USING ROUTINES SPLCC
C                          OR SPLCW.
C                               THE DIMENSION IS ASSUMED TO BE
C                                 COEF(NODES(1)*...*NODES(NDIM)).
C
C                        XMIN
C                          A VECTOR DESCRIBING THE LOWER EXTREME CORNER
C                          OF THE NODE GRID.  A SET OF EVENLY SPACED
C                          NODES IS FORMED ALONG EACH COORDINATE AXIS
C                          AND XMIN(IDIM) IS THE LOCATION OF THE FIRST
C                          NODE ALONG THE IDIM AXIS.
C                               THE DIMENSION IS ASSUMED TO BE
C                               XMIN(NDIM).
C
C                        XMAX
C                          A VECTOR DESCRIBING THE UPPER EXTREME CORNER
C                          OF THE NODE GRID.  A SET OF EVENLY SPACED
C                          NODES IS FORMED ALONG EACH COORDINATE AXIS
C                          AND XMAX(IDIM) IS THE LOCATION OF THE LAST
C                          NODE ALONG THE IDIM AXIS.
C                               THE DIMENSION IS ASSUMED TO BE
C                               XMAX(NDIM).
C
C                        NODES
C                          A VECTOR OF INTEGERS DESCRIBING THE NUMBER OF
C                          NODES ALONG EACH AXIS.  NODES(IDIM) IS THE
C                          THE NUMBER OF NODES (COUNTING ENDPOINTS)
C                          ALONG THE IDIM AXIS AND DETERMINES THE
C                          FLEXIBILITY OF THE SPLINE IN THAT COORDINATE
C                          DIRECTION.  NODES(IDIM) MUST BE .GE. 4 BUT
C                          MAY BE AS LARGE AS THE ARRAYS COEF AND WORK
C                          ALLOW.
C                               THE DIMENSION IS ASSUMED TO BE
C                               NODES(NDIM).
C                          NOTE:  THE NODE GRID IS COMPLETELY DEFINED BY
C                                 THE ARGUMENTS XMIN, XMAX AND NODES.
C                                 THE SPACING OF THIS GRID IN THE IDIM
C                                 COORDINATE DIRECTION IS
C                                 DX(IDIM) = (XMAX(IDIM)-XMIN(IDIM)) /
C                                             (NODES(IDIM)-1).
C                                 A NODE IN THIS GRID MAY BE INDEXED BY
C                                 AN NDIM VECTOR OF INTEGERS
C                                 (IN(1),...,IN(NDIM)) WHERE
C                                 1 .LE. IN(IDIM) .LE. NODES(IDIM).  THE
C                                 LOCATION OF SUCH A NODE MAY BE
C                                 REPRESENTED BY AN NDIM VECTOR
C                                 (X(1),...,X(NDIM))  WHERE
C                                 X(IDIM) = XMIN(IDIM)+(IN(IDIM)-1) *
C                                           DX(IDIM).
C
C ON OUTPUT              SPLDE
C                          THE FUNCTION VALUE RETURNED IS THE PARTIAL
C                          DERIVATIVE (INDICATED BY NDERIV) OF THE
C                          SPLINE EVALUATED AT X.
C
C                        IERROR
C                          AN ERROR FLAG WITH THE FOLLOWING MEANINGS--
C                              0  NO ERROR.
C                            101  NDIM IS .LT. 1 OR IS .GT. 4.
C                            102  NODES(IDIM) IS .LT. 4 FOR SOME IDIM.
C                            103  XMIN(IDIM) = XMAX(IDIM) FOR SOME IDIM.
C                            104  NDERIV(IDIM) IS .LT. 0 OR IS .GT. 2
C                                 FOR SOME IDIM.
C
C ALGORITHM              THE MULTI-DIMENSIONAL BASIS FUNCTIONS ARE
C                        FORMED AS TENSOR PRODUCTS OF ONE-DIMENSIONAL
C                        BASIS FUNCTIONS.  THESE ARE EVALUATED AT X,
C                        MULTIPLIED BY THE APPROPRIATE COEFFICIENT,
C                        AND SUMMED TO FORM THE RESULT.
C
C ACCURACY               ESSENTIALLY MACHINE PRECISION
C
C TIMING                 ROUGHLY PROPORTIONAL TO (NDIM-1).
C***********************************************************************
      FUNCTION SPLDE (NDIM,X,NDERIV,COEF,XMIN,XMAX,NODES,IERROR)
      DIMENSION       X(NDIM)    ,NDERIV(NDIM)           ,COEF(1)    ,
     1                XMIN(NDIM) ,XMAX(NDIM) ,NODES(NDIM)
      COMMON /SPLCOM/ MDIM       ,DX(4)      ,DXIN(4)    ,IB(4)      ,
     1                IBMN(4)    ,IBMX(4)
      SAVE
C
C
C THE RESTRICTION FOR NDIM TO BE LE 4 CAN BE ELIMINATED BY INCREASING
C THE ABOVE DIMENSIONS.
C
      IERROR = 0
      MDIM = NDIM
      IF (MDIM.LT.1 .OR. MDIM.GT.4) GO TO 105
      IIBMX = 1
      DO 101 IDIM=1,MDIM
         NOD = NODES(IDIM)
         IF (NOD .LT. 4) GO TO 106
         XRNG = XMAX(IDIM)-XMIN(IDIM)
         IF (XRNG .EQ. 0.) GO TO 107
         IF (NDERIV(IDIM).LT.0 .OR. NDERIV(IDIM).GT.2) GO TO 108
C
C DX(IDIM) IS THE NODE SPACING ALONG THE IDIM COORDINATE
C
         DX(IDIM) = XRNG/FLOAT(NOD-1)
         DXIN(IDIM) = 1./DX(IDIM)
C
C COMPUTE INDICES OF BASIS FUNCTIONS WHICH ARE NONZERO AT X
C
         IT = DXIN(IDIM)*(X(IDIM)-XMIN(IDIM))
C
C IBMN MUST BE IN RANGE 0 TO NODES-2
C
         IBMN(IDIM) = MIN0(MAX0(IT-1,0),NOD-2)
C
C IBMX MUST BE IN RANGE 1 TO NODES-1
C
         IBMX(IDIM) = MAX0(MIN0(IT+2,NOD-1),1)
         IIBMX = IIBMX*(IBMX(IDIM)-IBMN(IDIM)+1)
         IB(IDIM) = IBMN(IDIM)
  101 CONTINUE
      SUM = 0.
      IIB = 0
C
C BEGINING OF BASIS INDEX LOOP- TRAVERSE ALL INDICES CORRESPONDING
C TO BASIS FUNCTIONS WHICH ARE NONZERO AT X
C
  102 IIB = IIB+1
C
C THE INDICES ARE IN IB AND PASS THRU COMMON TO BASCMP
C
      CALL BASCMP (X,NDERIV,XMIN,NODES,ICOF,BASM)
C
C BASCMP COMPUTES ICOF AND BASM WHERE BASM IS THE VALUE (AT X) OF THE
C N-DIMENSIONAL BASIS FUNCTION CORRESPONDING TO COEF(ICOF)
C
      SUM = SUM+COEF(ICOF)*BASM
      IF (IIB .GE. IIBMX) GO TO 104
C
C INCREMENT BASIS INDICES
C
      DO 103 IDIM=1,MDIM
         IB(IDIM) = IB(IDIM)+1
         IF (IB(IDIM) .LE. IBMX(IDIM)) GO TO 102
         IB(IDIM) = IBMN(IDIM)
  103 CONTINUE
C
C END OF BASIS INDEX LOOP
C
  104 SPLDE = SUM
      RETURN
  105 CONTINUE
      IERROR = 101
      CALL ULIBER (IERROR,60H SPLFE OR SPLDE - NDIM IS LT 1 OR IS GT 4
     1                 ,60)
      GO TO 109
  106 CONTINUE
      IERROR = 102
      CALL ULIBER (IERROR,60H SPLFE OR SPLDE - NODES(IDIM) IS LT 4 FOR S
     1OME IDIM         ,60)
      GO TO 109
  107 CONTINUE
      IERROR = 103
      CALL ULIBER (IERROR,60H SPLFE OR SPLDE - XMIN(IDIM) EQ XMAX(IDIM)
     1FOR SOME IDIM    ,60)
      GO TO 109
  108 CONTINUE
      IERROR = 104
      CALL ULIBER (IERROR,60H SPLDE - NDERIV(IDIM) IS LT 0 OR IS GT 2 FO
     1R SOME IDIM      ,60)
  109 STOP
      END
C PACKAGE SPLPAK         DOCUMENTATION FOR USER ENTRIES FOLLOWS
C                        THE GENERAL PACKAGE INFORMATION.
C
C LATEST REVISION        MARCH  1985
C
C PURPOSE                THIS PACKAGE CONTAINS ROUTINES FOR FITTING
C                        (LEAST SQUARES) A MULTIDIMENSIONAL CUBIC SPLINE
C                        TO ARBITRARILY LOCATED DATA.  IT ALSO CONTAINS
C                        ROUTINES FOR EVALUATING THIS SPLINE (OR ITS
C                        PARTIAL DERIVATIVES) AT ANY POINT.
C
C                        COEFFICIENT CALCULATION IS PERFORMED IN
C                        SUBROUTINES SPLCC OR SPLCW AND EVALUATION IS
C                        PERFORMED BY FUNCTIONS SPLFE OR SPLDE.
C
C USAGE                  PACKAGE SPLPAK CONTAINS FOUR USER ENTRIES --
C                        SPLCC, SPLCW, SPLFE, AND SPLDE.
C
C                        THE USER FIRST CALLS SPLCC BY
C
C                          CALL SPLCC (NDIM,XDATA,L1XDAT,YDATA,NDATA,
C                                      XMIN,XMAX,NODES,XTRAP,COEF,NCF,
C                                      WORK,NWRK,IERROR)
C
C                        OR SPLCW BY
C
C                          CALL SPLCW (NDIM,XDATA,L1XDATA,YDATA,WDATA,
C                                      NDATA,XMIN,XMAX,NODES,XTRAP,
C                                      COEF,NCF,WORK,NWRK,IERROR)
C
C                        THE PARAMETER NDATA IN THE CALL TO SPLCW
C                        ENABLES THE USER TO WEIGHT SONE OF THE DATA
C                        POINTS MORE HEAVILY THAN OTHERS.  BOTH
C                        ROUTINES RETURN A SET OF COEFFICIENTS IN THE
C                        ARRAY COEF.  THESE COEFFICIENTS ARE
C                        SUBSEQUENTLY USED IN THE COMPUTATION OF
C                        FUNCTION VALUES AND PARTIAL DERIVATIVES.
C                        THE USER THEN CALLS SPLFE OR SPLDE ANY
C                        NUMBER OF TIMES IN ANY ORDER PROVIDED THAT
C                        THE VALUES OF THE INPUTS, NDIM,COEF,XMIN,
C                        XMAX, AND NODES, ARE PRESERVED BETWEEN CALLS.
C                        THE ROUTINES ARE CALLED IN THE FOLLOWING WAY.
C
C                          F = SPLFE (NDIM,X,COEF,XMIN,XMAX,NODES,
C                                     IERROR)
C
C                        OR
C
C                          F = SPLDE (NDIM,X,NDERIV,COEF,XMIN,XMAX,
C                                     NODES,IERROR)
C
C                        THE ROUTINE, SPLFE,  RETURNS AN INTERPOLATED
C                        VALUE AT THE POINT DEFINED BY THE ARRAY X.
C                        XPLDE AFFORDS THE USER THE ADDITIONAL
C                        CAPABILITY OF CALCULATING AN INTERPOLATED
C                        VALUE FOR ONE OF SEVERAL PARTIAL DERIVATIVES
C                        SPECIFIED BY THE ARRAY NDERIV.
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE, EXCEPT FOR ERROR MESSAGES PRINTED BY
C                        CALLS TO ULIBER, IF AN ERROR IS DETECTED.
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       SUPRLS, ULIBER
C FILES
C
C LANGUAGE               FORTRAN
C
C HISTORY                DEVELOPED IN 1972-73 BY DAVE FULKER OF NCAR'S
C                        SCIENTIFIC COMPUTING DIVISION.
C
C PORTABILITY            FORTRAN 66
C***********************************************************************
C FUNCTION SPLFE (NDIM,X,COEF,XMIN,XMAX,NODES,IERROR)
C
C DIMENSION OF           X(NDIM),NDERIV(NDIM),
C ARGUMENTS              COEF(NODES(1)*...*NODES(NDIM)),XMIN(NDIM),
C                        XMAX(NDIM),NODES(NDIM)
C
C NOTE                   COEF, XMIN, XMAX AND NODES MUST BE EXACTLY
C                        RETAINED FROM THE CALL TO SPLCC (OR SPLCW).
C
C PURPOSE                N-DIMENSIONAL CUBIC SPLINE FUNCTION EVALUATION.
C
C USAGE                  EXCEPT FOR LACK OF DERIVATIVE CAPABILITY, THIS
C                        FUNCTION IS IDENTICAL TO FUNCTION SPLDE IN
C                        USAGE.  THE ARGUMENT LIST IS ALSO IDENTICAL
C                        EXCEPT FOR THE OMISSION OF NDERIV.  SEE
C                        FUNCTION SPLDE DESCRIPTION IMMEDIATELY BELOW
C                        FOR A COMPLETE DESCRIPTION.
C
C                        Y = SPLFE (NDIM,X,COEF,XMIN,XMAX,NODES,IERROR)
C***********************************************************************
C
C FUNCTION SPLDE (NDIM,X,NDERIV,COEF,XMIN,XMAX,NODES,IERROR)
C
C DIMENSION OF           X(NDIM),NDERIV(NDIM),
C ARGUMENTS              COEF(NODES(1)*...*NODES(NDIM)),XMIN(NDIM),
C                        XMAX(NDIM),NODES(NDIM)
C
C PURPOSE                N-DIMENSIONAL CUBIC SPLINE DERIVATIVE
C                        EVALUATION.
C
C                        A GRID OF EVENLY SPACED NODES IN NDIM SPACE IS
C                        DEFINED BY THE ARGUMENTS XMIN, XMAX AND NODES.
C                        A LINEAR BASIS FOR THE CLASS OF NATURAL SPLINES
C                        ON THESE NODES IS FORMED, AND TO EACH BASIS
C                        FUNCTION CORRESPONDS A COEFFICIENT IN THE ARRAY
C                        COEF (COMPUTED IN SPLCC OR SPLCW).  USING
C                        NDERIV TO INDICATE THE APPROPRIATE PARTIAL
C                        DERIVATIVES, EACH BASIS FUNCTION IS EVALUATED
C                        AT THE POINT X IN NDIM SPACE.  THESE VALUES ARE
C                        THEN MULTIPLIED BY THE CORRESPONDING
C                        COEFFICIENT AND SUMMED TO FORM THE FUNCTION
C                        RESULT.
C
C USAGE                  COEF, XMIN, XMAX AND NODES MUST BE EXACTLY
C                        RETAINED FROM THE CALL TO SPLCC (OR SPLCW).
C
C                        THEN EVALUATE PARTIAL DERIVATIVES OF THE SPLINE
C                        SET BY
C
C                        Y = SPLDE (NDIM,X,NDERIV,COEF,XMIN,XMAX,NODES,
C                                   IERROR).
C
C ARGUMENTS
C
C ON INPUT               NDIM
C                          THE DIMENSIONALITY OF THE PROBLEM.  THE
C                          SPLINE IS A FUNCTION OF NDIM VARIABLES OR
C                          COORDINATES AND THUS A POINT IN THE
C                          INDEPENDENT VARIABLE SPACE IS AN NDIM VECTOR.
C                          NDIM MUST BE IN THE RANGE 1 .LE. NDIM .LE. 4.
C
C                        X
C                          AN NDIM VECTOR DESCRIBING THE POINT IN THE
C                          INDEPENDENT VARIABLE SPACE AT WHICH THE
C                          SPLINE IS TO BE EVALUATED.
C                               THE DIMENSION IS ASSUMED TO BE X(NDIM).
C
C                        NDERIV
C                          AN NDIM VECTOR OF INTEGERS SPECIFYING THE
C                          PARTIAL DERIVATIVE TO BE EVALUATED.  THE
C                          ORDER OF THE DERIVATIVE ALONG THE IDIM AXIS
C                          IS NDERIV(IDIM).  THESE INTEGERS MUST BE IN
C                          THE RANGE 0 .LE. NDERIV(IDIM) .LE. 2.
C
C                        COEF
C                          THE ARRAY OF COEFFICIENTS WHICH DETERMINE THE
C                          SPLINE.  EACH COEFFICIENT CORRESPONDS TO A
C                          PARTICULAR BASIS FUNCTION WHICH IN TURN
C                          CORRESPONDS TO A NODE IN THE NODE GRID.  THIS
C                          CORRESPONDENCE BETWEEN THE NODE GRID AND THE
C                          ARRAY COEF IS AS IF COEF WERE AN
C                          NDIM-DIMENSIONAL FORTRAN ARRAY WITH
C                          DIMENSIONS NODES(1),...,NODES(NDIM), I.E., TO
C                          STORE THE ARRAY LINEARLY, THE LEFTMOST
C                          INDICES ARE INCREMENTED MOST FREQUENTLY.
C                          COEF MAY BE COMPUTED BY USING ROUTINES SPLCC
C                          OR SPLCW.
C                               THE DIMENSION IS ASSUMED TO BE
C                                 COEF(NODES(1)*...*NODES(NDIM)).
C
C                        XMIN
C                          A VECTOR DESCRIBING THE LOWER EXTREME CORNER
C                          OF THE NODE GRID.  A SET OF EVENLY SPACED
C                          NODES IS FORMED ALONG EACH COORDINATE AXIS
C                          AND XMIN(IDIM) IS THE LOCATION OF THE FIRST
C                          NODE ALONG THE IDIM AXIS.
C                               THE DIMENSION IS ASSUMED TO BE
C                               XMIN(NDIM).
C
C                        XMAX
C                          A VECTOR DESCRIBING THE UPPER EXTREME CORNER
C                          OF THE NODE GRID.  A SET OF EVENLY SPACED
C                          NODES IS FORMED ALONG EACH COORDINATE AXIS
C                          AND XMAX(IDIM) IS THE LOCATION OF THE LAST
C                          NODE ALONG THE IDIM AXIS.
C                               THE DIMENSION IS ASSUMED TO BE
C                               XMAX(NDIM).
C
C                        NODES
C                          A VECTOR OF INTEGERS DESCRIBING THE NUMBER OF
C                          NODES ALONG EACH AXIS.  NODES(IDIM) IS THE
C                          THE NUMBER OF NODES (COUNTING ENDPOINTS)
C                          ALONG THE IDIM AXIS AND DETERMINES THE
C                          FLEXIBILITY OF THE SPLINE IN THAT COORDINATE
C                          DIRECTION.  NODES(IDIM) MUST BE .GE. 4 BUT
C                          MAY BE AS LARGE AS THE ARRAYS COEF AND WORK
C                          ALLOW.
C                               THE DIMENSION IS ASSUMED TO BE
C                               NODES(NDIM).
C                          NOTE:  THE NODE GRID IS COMPLETELY DEFINED BY
C                                 THE ARGUMENTS XMIN, XMAX AND NODES.
C                                 THE SPACING OF THIS GRID IN THE IDIM
C                                 COORDINATE DIRECTION IS
C                                 DX(IDIM) = (XMAX(IDIM)-XMIN(IDIM)) /
C                                             (NODES(IDIM)-1).
C                                 A NODE IN THIS GRID MAY BE INDEXED BY
C                                 AN NDIM VECTOR OF INTEGERS
C                                 (IN(1),...,IN(NDIM)) WHERE
C                                 1 .LE. IN(IDIM) .LE. NODES(IDIM).  THE
C                                 LOCATION OF SUCH A NODE MAY BE
C                                 REPRESENTED BY AN NDIM VECTOR
C                                 (X(1),...,X(NDIM))  WHERE
C                                 X(IDIM) = XMIN(IDIM)+(IN(IDIM)-1) *
C                                           DX(IDIM).
C
C ON OUTPUT              SPLDE
C                          THE FUNCTION VALUE RETURNED IS THE PARTIAL
C                          DERIVATIVE (INDICATED BY NDERIV) OF THE
C                          SPLINE EVALUATED AT X.
C
C                        IERROR
C                          AN ERROR FLAG WITH THE FOLLOWING MEANINGS--
C                              0  NO ERROR.
C                            101  NDIM IS .LT. 1 OR IS .GT. 4.
C                            102  NODES(IDIM) IS .LT. 4 FOR SOME IDIM.
C                            103  XMIN(IDIM) = XMAX(IDIM) FOR SOME IDIM.
C                            104  NDERIV(IDIM) IS .LT. 0 OR IS .GT. 2
C                                 FOR SOME IDIM.
C
C ALGORITHM              THE MULTI-DIMENSIONAL BASIS FUNCTIONS ARE
C                        FORMED AS TENSOR PRODUCTS OF ONE-DIMENSIONAL
C                        BASIS FUNCTIONS.  THESE ARE EVALUATED AT X,
C                        MULTIPLIED BY THE APPROPRIATE COEFFICIENT,
C                        AND SUMMED TO FORM THE RESULT.
C
C ACCURACY               ESSENTIALLY MACHINE PRECISION
C
C TIMING                 ROUGHLY PROPORTIONAL TO (NDIM-1).
C***********************************************************************
      FUNCTION SPLFE (NDIM,X,COEF,XMIN,XMAX,NODES,IERROR)
      DIMENSION       X(NDIM)    ,COEF(1)    ,XMIN(NDIM) ,XMAX(NDIM) ,
     1                NODES(NDIM)
      DIMENSION       NDERIV(4)
      SAVE
C
      DATA NDERIV(1),NDERIV(2),NDERIV(3),NDERIV(4)/0,0,0,0/
C
C THE RESTRICTION FOR NDIM TO BE LE 4 CAN BE ELIMINATED BY INCREASING
C THE ABOVE DIMENSION AND THOSE IN SPLDE.
C
      SPLFE = SPLDE(NDIM,X,NDERIV,COEF,XMIN,XMAX,NODES,IERROR)
      RETURN
      END
C     SUBROUTINE SUPRLS (I,ROWI,N,BI,A,NN,SOLN,ERR,IER)
C
C
C DIMENSION OF           ROWI(N), A(NN), SOLN(N)
C ARGUMENTS
C
C LATEST REVISION        MARCH 1985
C
C PURPOSE                TO DETERMINE THE LEAST SQUARES SOLUTION OF A
C                        LARGE OVERDETERMINED LINEAR SYSTEM.  GIVEN THE
C                        M BY N MATRIX R (M .GE. N) AND THE M-VECTOR B,
C                        THIS ROUTINE CALCULATES THE N-VECTOR X SUCH
C                        THAT THE EUCLIDEAN NORM OF THE RESIDUE (R*X-B)
C                        IS MINIMIZED.  THE SUBROUTINE ACCEPTS ROWS OF
C                        THE MATRIX ONE BY ONE SO THAT THE ENTIRE MATRIX
C                        NEED NOT BE STORED AT ONE TIME.  THIS ALLOWS
C                        LARGE PROBLEMS TO BE SOLVED WITHOUT PERIPHERAL
C                        STORAGE.  THE LENGTH OF THE ROWS IS LIMITED BY
C                        THE AMOUNT OF SCRATCH STORAGE WHICH CAN BE SET
C                        ASIDE FOR USE BY THE ROUTINE.  THERE IS NO
C                        RESTRICTION ON THE NUMBER OF ROWS.
C
C USAGE                  THIS PACKAGE HAS ONLY ONE USER ENTRY,
C                        SUBROUTINE SUPRLS. SUPRLS IS CALLED ONCE FOR
C                        EACH ROW OF THE MATRIX.  A FINAL CALL RETURNS
C                        THE SOLUTION VECTOR AND THE EUCLIDEAN NORM
C                        OF THE RESIDUAL. THIS FOLLOWING SEQUENCE WOULD
C                        PROCESS THE M BY N MATRIX R AND THE RIGHT HAND
C                        SIDE M-VECTOR B
C                               DO  1  I = 1,M
C                               DO  2  J = 1,N
C                        HERE SET ROWI(J) TO THE (I,J) ELEMENT OF R
C                            2  CONTINUE
C                        HERE SET BI TO THE ITH COMPONENT OF B.
C                               CALL SUPRLS(I,ROWI,N,BI,A,NN,SOLN,ERR,
C                                           IER)
C                            1  CONTINUE
C                               CALL SUPRLS (0,ROWI,N,BI,A,NN,SOLN,ERR,
C                                           IER)
C
C ARGUMENTS
C
C ON INPUT               I
C                          THE INDEX OF THE ROW BEING ENTERED.  (I IS 1
C                          FOR THE FIRST CALL, INCREASES BY ONE FOR EACH
C                          CALL, AND IS M WHEN THE FINAL ROW IS
C                          ENTERED).  AFTER THE FINAL ROW HAS BEEN
C                          ENTERED, SUPRLS IS CALLED WITH I = 0 TO
C                          COMPLETE THE REDUCTION AND SOLUTION.
C
C                        ROWI
C                          A VECTOR WHICH ON THE ITH CALL CONTAINS THE N
C                          COMPONENTS OF THE ITH ROW OF THE MATRIX.  THE
C                          DIMENSION OF ROWI IN CALLING PROGRAM MUST BE
C                          AT LEAST N.
C
C                        N
C                          THE LENGTH OF THE ROWS OF THE MATRIX (I.E.,
C                          THE NUMBER OF COLUMNS).  N .LE. M, WHERE M IS
C                          THE NUMBER OF ROWS.
C
C                        BI
C                          ON THE ITH CALL, BI CONTAINS THE ITH ELEMENT
C                          OF THE RIGHT HAND SIDE VECTOR B.
C
C                        A
C                          A WORKING ARRAY WHICH MUST NOT BE CHANGED
C                          BETWEEN THE SUCCESSIVE CALLS TO SUPRLS.
C                          DIMENSION OF A IN CALLING PROGRAM IS NN.
C
C                        NN
C                          LENGTH OF SCRATCH ARRAY A.  NN MUST BE AT
C                          LEAST N*(N+5)/2+1.  FOR SPEED, NN SHOULD BE
C                          AS LARGE AS POSSIBLE UP TO A MAXIMUM OF
C                          (N+1)*M.
C
C ON OUTPUT              SOLN
C                          THE N-COMPONENTS OF THE SOLUTION VECTOR ARE
C                          RETURNED IN THIS ARRAY AFTER THE FINAL CALL
C                          TO SUPRLS.  DIMENSION OF SOLN IN CALLING
C                          PROGRAM MUST BE AT LEAST N.
C
C                        ERR
C                          THE EUCLIDEAN NORM OF THE RESIDUAL IS
C                          RETURNED IN ERR AFTER THE FINAL CALL TO
C                          SUPRLS.
C
C                        IER
C                          ERROR PARAMETER.
C                          FATAL ERRORS.
C                          = 32  INSUFFICIENT SCRATCH STORAGE PROVIDED,
C                                MUST HAVE NN .GE. N*(N+5)/2+1.
C                          = 33  ARRAY HAS TOO FEW ROWS.  MUST HAVE
C                                M .GE. N.
C                          = 34  SYSTEM IS SINGULAR.
C                          = 35  VALUES OF I NOT IN SEQUENCE.
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    ERROR MESSAGES ARE WRITTEN TO UNIT 6 BY
C                        ROUTINE ULIBER.
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       ULIBER
C FILES
C
C LANGUAGE               FORTRAN
C
C HISTORY                WRITTEN IN MAY 1972 BY A.K. CLINE OF NCAR'S
C                        SCIENTIFIC COMPUTING DIVISION.
C
C PORTABILITY            FORTRAN 66
C
C ALGORITHM              GIVEN THE  M BY N  MATRIX  R (M.GE.N) AND THE
C                        M-VECTOR B, WE WISH TO FIND AN N-VECTOR X
C                        SUCH THAT
C
C                            E = L2-NORM OF  R*X-B
C
C                        IS MINIMIZED.  SINCE THE EUCLIDEAN NORM IS
C                        INVARIANT UNDER ORTHOGONAL TRANSFORMATION,
C                        R AND B MAY BE PREMULTIPLIED BY ANY
C                        ORTHOGONAL MATRIX WITHOUT CHANGING THE
C                        NORM OF THE RESIDUAL (R*X-B).  R IS REDUCED
C                        TO UPPER TRIANGULAR FORM BY PREMULTIPLYING
C                        R AND B BY A SEQUENCE OF HOUSEHOLDER AND
C                        ROTATION MATRICES.  WHEN THE REDUCTION IS
C                        COMPLETE, THE NORM OF THE RESIDUAL TAKES THE
C                        FORM
C
C                          E =  L2 NORM(T*X-B(N))+L2 NORM(B(M-N))
C
C                        WHERE T IS AN  N BY N  UPPER TRIANGULAR
C                        MATRIX, B(N) IS A VECTOR OF THE FIRST N
C                        COMPONENTS OF B, B(M-N) IS A VECTOR OF
C                        THE REMAINING (M-N) COMPONENTS OF B.  E IS
C                        MINIMIZED BY TAKING X TO BE THE SOLUTION
C                        OF THE SYSTEM  T*X=B(N).  THIS TRIANGULAR
C                        SYSTEM IS THEREFORE SOLVED TO GIVE THE
C                        REQUIRED LEAST SQUARES SOLUTION.  THE NORM
C                        OF THE RESIDUAL IS THEN THE L2-NORM OF B(M-N).
C
C                        AT EACH PHASE OF THE REDUCTION, AS MANY ROWS
C                        AS SPACE PERMITS ARE ENTERED INTO THE SCRATCH
C                        AREA.  HOUSEHOLDER TRANSFORMATIONS ARE THEN
C                        USED TO ZERO OUT SUBDIAGONAL ELEMENTS.  SPACE
C                        IS SAVED BY ELIMINATING STORAGE FOR THE
C                        ZERO SUBDIAGONAL TERMS.  IF THERE IS ROOM
C                        FOR ONLY ONE NEW ROW, ROTATION RATHER THAN
C                        HOUSEHOLDER MATRICES ARE USED FOR GREATER
C                        SPEED.  WHEN ALL  M  ROWS HAVE BEEN ENTERED,
C                        REDUCTION IS COMPLETED AND THE TRIANGULAR
C                        SYSTEM SOLVED.
C
C                        FOR GREATER DETAIL SEE HANSON, R.J., AND
C                        LAWSON, C.L., 1969--EXTENSIONS AND
C                        APPLICATIONS OF THE HOUSEHOLDER ALGORITHM
C                        FOR SOLVING LINEAR LEAST SQUARES PROBLEMS.
C                        MATH. OF COMP. VOL.23, PP. 787-812.
C
C ACCURACY               THIS WILL DEPEND UPON THE SIZE AND CONDITION
C                        OF THE MATRIX.  NEAR MACHINE ACCURACY MAY BE
C                        EXPECTED FOR WELL CONDITIONED SYSTEMS OF
C                        MODERATE SIZE.  IF ILL CONDITIONING IS
C                        SUSPECT, A VERSION USING PIVOTING MAY BE
C                        NECESSARY.
C
C TIMING                 THIS DEPENDS NOT ONLY UPON THE DIMENSIONS OF
C                        THE MATRIX, BUT ALSO UPON THE AMOUNT OF
C                        SCRATCH STORAGE AVAILABLE.
C***********************************************************************
      SUBROUTINE SUPRLS (I,ROWI,N,BI,A,NN,SOLN,ERR,IER)
      DIMENSION       ROWI(N)    ,A(NN)      ,SOLN(N)
      SAVE
C
      IER = 0
      IF (I .GT. 1) GO TO 101
C
C ROUTINE ENTERED WITH I.LE.0 MEANS COMPLETE THE REDUCTION AND STORE
C THE SOLUTION IN SOLN
C
      IF (I .LE. 0) GO TO 125
C
C SET UP QUANTITIES ON FIRST CALL
C
      IOLD = 0
      NP1 = N+1
C
C COMPUTE HOW MANY ROWS CAN BE INPUT NOW
C
      L = NN/NP1
      ILAST = 0
      IL1 = 0
      K = 0
      K1 = 0
      ERRSUM = 0.
      NREQ = ((N+5)*N+2)/2
C
C ERROR EXIT IF INSUFFICIENT SCRATCH STORAGE PROVIDED.
C
      IF (NN .GE. NREQ) GO TO 101
      IER = 32
      CALL ULIBER (IER,86H SUPRLS-INSUFFICIENT SCRATCH STORAGE PROVIDED.
     1 AT LEAST ((N+5)*N+2)/2 LOCATIONS NEEDED,86)
      RETURN
C
C STORE THE ROW IN THE SCRATCH STORAGE
C
  101 CONTINUE
C
C ERROR EXIT IF (I-IOLD).NE.1
C
      IF ((I-IOLD) .EQ. 1) GO TO 102
      IER = 35
      CALL ULIBER (IER,35H SUPRLS-VALUES OF I NOT IN SEQUENCE,35)
      RETURN
  102 CONTINUE
      IOLD = I
      DO 103 J=1,N
         ILJ = ILAST+J
         A(ILJ) = ROWI(J)
  103 CONTINUE
      ILNP = ILAST+NP1
      A(ILNP) = BI
      ILAST = ILAST+NP1
      ISAV = I
      IF (I .LT. L) RETURN
  104 CONTINUE
      IF (K .EQ. 0) GO TO 115
      K1 = MIN0(K,N)
      IDIAG = -NP1
      IF (L-K .EQ. 1) GO TO 110
C
C APPLY HOUSEHOLDER TRANSFORMATIONS TO ZERO OUT NEW ROWS
C
      DO 109 J=1,K1
         IDIAG = IDIAG+(NP1-J+2)
         I1 = IL1+J
         I2 = I1+NP1*(L-K-1)
         S = A(IDIAG)*A(IDIAG)
         DO 105 II=I1,I2,NP1
            S = S+A(II)*A(II)
  105    CONTINUE
         IF (S .EQ. 0.) GO TO 109
         TEMP = A(IDIAG)
         A(IDIAG) = SQRT(S)
         IF (TEMP .GT. 0.) A(IDIAG) = -A(IDIAG)
         TEMP = TEMP-A(IDIAG)
         TEMP1 = 1./(TEMP*A(IDIAG))
         JP1 = J+1
         DO 108 J1=JP1,NP1
            JDEL = J1-J
            IDJ = IDIAG+JDEL
            S = TEMP*A(IDJ)
            DO 106 II=I1,I2,NP1
               IIJD = II+JDEL
               S = S+A(II)*A(IIJD)
  106       CONTINUE
            S = S*TEMP1
            A(IDJ) = A(IDJ)+S*TEMP
            DO 107 II=I1,I2,NP1
               IIJD = II+JDEL
               A(IIJD) = A(IIJD)+S*A(II)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
      GO TO 113
C
C APPLY ROTATIONS TO ZERO OUT THE SINGLE NEW ROW
C
  110 DO 112 J=1,K1
         IDIAG = IDIAG+(NP1-J+2)
         I1 = IL1+J
         S = SQRT(A(IDIAG)*A(IDIAG)+A(I1)*A(I1))
         IF (S .EQ. 0.) GO TO 112
         TEMP = A(IDIAG)
         A(IDIAG) = S
         S = 1./S
         CN = TEMP*S
         SN = A(I1)*S
         JP1 = J+1
         DO 111 J1=JP1,NP1
            JDEL = J1-J
            IDJ = IDIAG+JDEL
            TEMP = A(IDJ)
            I1JD = I1+JDEL
            A(IDJ) = CN*TEMP+SN*A(I1JD)
            A(I1JD) = -SN*TEMP+CN*A(I1JD)
  111    CONTINUE
  112 CONTINUE
  113 IF (K .LT. N) GO TO 115
      LMKM1 = L-K
C
C ACCUMULATE RESIDUAL SUM OF SQUARES
C
      DO 114 II=1,LMKM1
         ILNP = IL1+II*NP1
         ERRSUM = ERRSUM+A(ILNP)*A(ILNP)
  114 CONTINUE
      IF (I .LE. 0) GO TO 127
      K = L
      ILAST = IL1
C
C DETERMINE HOW MANY NEW ROWS MAY BE INPUT ON NEXT ITERATION
C
      L = K+(NN-ILAST)/NP1
      RETURN
  115 K11 = K1+1
      K1 = MIN0(L,N)
      IF (L-K .EQ. 1) GO TO 122
      K1M1 = K1-1
      IF (L .GT. N) K1M1 = N
      I1 = IL1+K11-NP1-1
C
C PERFORM HOUSEHOLDER TRANSFORMATIONS TO REDUCE ROWS TO UPPER TRIANGULAR
C FORM
C
      DO 120 J=K11,K1M1
         I1 = I1+(NP1+1)
         I2 = I1+(L-J)*NP1
         S = 0.
         DO 116 II=I1,I2,NP1
            S = S+A(II)*A(II)
  116    CONTINUE
         IF (S .EQ. 0.) GO TO 120
         TEMP = A(I1)
         A(I1) = SQRT(S)
         IF (TEMP .GT. 0.) A(I1) = -A(I1)
         TEMP = TEMP-A(I1)
         TEMP1 = 1./(TEMP*A(I1))
         JP1 = J+1
         I11 = I1+NP1
         DO 119 J1=JP1,NP1
            JDEL = J1-J
            I1JD = I1+JDEL
            S = TEMP*A(I1JD)
            DO 117 II=I11,I2,NP1
               IIJD = II+JDEL
               S = S+A(II)*A(IIJD)
  117       CONTINUE
            S = S*TEMP1
            I1JD = I1+JDEL
            A(I1JD) = A(I1JD)+S*TEMP
            DO 118 II=I11,I2,NP1
               IIJD = II+JDEL
               A(IIJD) = A(IIJD)+S*A(II)
  118       CONTINUE
  119    CONTINUE
  120 CONTINUE
      IF (L .LE. N) GO TO 122
      NP1MK = NP1-K
      LMK = L-K
C
C ACCUMULATE RESIDUAL SUM OF SQUARES
C
      DO 121 II=NP1MK,LMK
         ILNP = IL1+II*NP1
         ERRSUM = ERRSUM+A(ILNP)*A(ILNP)
  121 CONTINUE
  122 IMOV = 0
      I1 = IL1+K11-NP1-1
C
C SQUEEZE THE UNNECESSARY ELEMENTS OUT OF SCRATCH STORAGE TO ALLOW SPACE
C FOR MORE ROWS
C
      DO 124 II=K11,K1
         IMOV = IMOV+(II-1)
         I1 = I1+NP1+1
         I2 = I1+NP1-II
         DO 123 III=I1,I2
            IIIM = III-IMOV
            A(IIIM) = A(III)
  123    CONTINUE
  124 CONTINUE
      ILAST = I2-IMOV
      IL1 = ILAST
      IF (I .LE. 0) GO TO 127
      K = L
C
C DETERMINE HOW MANY NEW ROWS MAY BE INPUT ON NEXT ITERATION
C
      L = K+(NN-ILAST)/NP1
      RETURN
C
C COMPLETE REDUCTION AND STORE SOLUTION IN SOLN
C
  125 L = ISAV
C
C ERROR EXIT IF L LESS THAN N
C
      IF (L .GE. N) GO TO 126
      IER = 33
      CALL ULIBER (IER,31H SUPRLS-ARRAY HAS TOO FEW ROWS.,31)
  126 CONTINUE
C
C K NE ISAV MEANS FURTHER REDUCTION NEEDED
C
      IF (K .NE. ISAV) GO TO 104
  127 ILAST = (NP1*(NP1+1))/2-1
      IF (A(ILAST-1) .EQ. 0.) GO TO 130
C
C SOLVE TRIANGULAR SYSTEM INTO ROWI
C
      SOLN(N) = A(ILAST)/A(ILAST-1)
      DO 129 II=2,N
         IIM1 = II-1
         ILAST = ILAST-II
         S = A(ILAST)
         DO 128 K=1,IIM1
            ILK = ILAST-K
            NPK = NP1-K
            S = S-A(ILK)*SOLN(NPK)
  128    CONTINUE
         ILII = ILAST-II
         IF (A(ILII) .EQ. 0.) GO TO 130
         NPII = NP1-II
         SOLN(NPII) = S/A(ILII)
  129 CONTINUE
C
C STORE RESIDUAL NORM
C
      ERR = SQRT(ERRSUM)
      RETURN
C
C ERROR RETURN IF SYSTEM IS SINGULAR
C
  130 CONTINUE
      IER = 34
      CALL ULIBER (IER,27H SUPRLS-SYSTEM IS SINGULAR.,27)
      RETURN
C
C REVISION HISTORY---
C
C AUGUST 1977      MODIFIED DOCUMENTATION TO MATCH THAT IN
C                  VOLUME 1 OF THE NSSL MANUALS.  MODIFIED CODE
C                  TO ENHANCE PORTABILITY.
C
C JANUARY 1978     DELETED REFERENCES TO THE  *COSY  CARDS
C                  AND MOVED THE REVISION HISTORIES TO APPEAR BEFORE
C                  THE FINAL END CARD
C-----------------------------------------------------------------------
      END
C SUBROUTINE ULIBER (IERR,MESS,LMESS)                                   
C                                                                       
C PURPOSE                TO PRINT AN ERROR NUMBER AND AN ERROR MESSAGE  
C                        OR JUST AN ERROR MESSAGE.                      
C                                                                       
C USAGE                  CALL ULIBER (IERR,MESS,LMESS)                  
C                                                                       
C ARGUMENTS                                                             
C ON INPUT               IERR                                           
C                          THE ERROR NUMBER (PRINTED ONLY IF NON-ZERO). 
C                                                                       
C                        MESS                                           
C                          MESSAGE TO BE PRINTED.                       
C                                                                       
C                        LMESS                                          
C                          NUMBER OF CHARACTERS IN MESS (.LE. 130).     
C                                                                       
C ARGUMENTS                                                             
C ON OUTPUT              NONE                                           
C                                                                       
C I/O                    THE MESSAGE IS WRITEN TO UNIT 101.             
C ******************************************************************    
      SUBROUTINE ULIBER (IERR,MESS,LMESS)                               
      REAL MESS(1)                                                      
      SAVE
C                                                                       
      IF (IERR.NE.0) WRITE (101,1001) IERR                              
      NWORDS=(LMESS+7)/8                                                
      WRITE (101,1002) (MESS(I),I=1,NWORDS)                             
      RETURN                                                            
C                                                                       
 1001 FORMAT (6H0IERR=,I5)                                              
 1002 FORMAT (16A8,A2)                                                  
C                                                                       
      END                                                               

