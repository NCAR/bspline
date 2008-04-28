C     --------------------------------------------------------
      SUBROUTINE FLOPASS(GIN,GOUT,NMAX,FRAC,NTRM)
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

#ifdef notdef
      DIMENSION GIN(1000),GOUT(1000)
      nmax = 0

 601  format(R)
      do 60 i=1,1000
         read (*, 601) gin(i)
         nmax = nmax + 1
 60   continue


 701  format(R)
      do 70 i=1,nmax
         write (*, 601) gin(i)
 70   continue

      END
#endif
