      PROGRAM SIS
      IMPLICIT NONE
      INTEGER*4 SECTION,N,E,SEED,NPAS,I,J,K,L,NINF,EACT,RNDI,Z
      PARAMETER(SECTION=2,N=1039,E=1305,SEED=987654321,NPAS=10000)
      INTEGER*4 NI(N),EDGELIST(2*E,2),POINTERS(N,10),KI(N)
      INTEGER*4 INFLIST(N),ACTKI(N),COUNTER,COUNTER2,KAV
      REAL*8 T,LAMBDA,GENRAND_REAL2,INFECTED0,X,Y,RO1(N),RO2(N),AVERAGE
      REAL*8 DT,TMAX,LMBDLIST(18),SUM1,SUM2,MEANK,FULLYMIXEDRO
      PARAMETER(INFECTED0=0.1D0,KAV=5000,MEANK=2.512030798845043)
      OPEN(100,FILE="xc_connected_edgelist.dat")
      OPEN(200,FILE="xc_infected-time_10.dat")
      OPEN(300,FILE="xc_ro-time_10.dat")
      OPEN(400,FILE="xc_infectedeq-lambda.dat")
      CALL INIT_GENRAND(SEED)

      DO J=1,2*E ! READING THE EDGELIST
        READ(100,*) EDGELIST(J,1),EDGELIST(J,2)
      ENDDO
      CLOSE(100)
      DO J=1,2*E ! CONTRUCTING POINTERS AND DEGREES
        KI(EDGELIST(J,1))=KI(EDGELIST(J,1))+1 ! DEGREE OF EACH NODE
        POINTERS(EDGELIST(J,1),KI(EDGELIST(J,1)))=EDGELIST(J,2) ! NODES TO WHICH EACH NODE IS CONNECTED
      ENDDO
C           #INFECTED AS A FUNCTION OF TIME FOR A GIVEN LAMBDA
      IF (SECTION.EQ.1) THEN ! 
        NINF=0
        DO I=1,N ! INITIAL CONDITIONS AND COMPUTATION OF NUMBER OF INFECTED NODES
          X=GENRAND_REAL2()
          IF (X.LE.INFECTED0) THEN
            NI(I)=1
            NINF=NINF+1 ! NUMBER OF INFECTED NODES
            INFLIST(NINF)=I ! LIST OF INFECTED NODES
          ELSE
            NI(I)=0
          ENDIF
          RO1(I)=INFECTED0 ! AVERAGE INFECTED NODES OVER REALIZATIONS
        ENDDO
        EACT=0
        DO I=1,N ! COMPUTATION OF ACTIVE LINKS
          ACTKI(I)=0
          IF (NI(I).EQ.1) THEN
            DO J=1,KI(I)
              IF (NI(POINTERS(I,J)).EQ.0) THEN
                ACTKI(I)=ACTKI(I)+1 ! ACTIVE LINKS STARTING IN EACH NODE 
                !(OR NUMBER OF SUSCEPTIBLE NEIGHBOURS IF I IS INFECTED, 0 OTHERWISE)
                EACT=EACT+1 ! NUMBER ACTIVE LINKS
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        LAMBDA=10.D0
        T=0.D0
        WRITE(200,*) T,DBLE(NINF)/DBLE(N)
        DO K=1,NPAS ! LOOP OVER TIMESTEPS
          Y=(NINF+LAMBDA*EACT)*GENRAND_REAL2()
          IF (Y.LT.DBLE(NINF)) THEN ! RECOVERY
            Z=INT(NINF*GENRAND_REAL2())+1
            RNDI=INFLIST(Z) ! NODE THAT GETS RECOVERED
            NI(RNDI)=0
            DO I=Z,NINF
              INFLIST(I)=INFLIST(I+1) ! IT IS NOT OPTIMIZED, BUT THERE IS NO NEED FOR THAT N
            ENDDO
            NINF=NINF-1
            EACT=EACT+KI(RNDI)-2*ACTKI(RNDI)
            ! #LINKS THAT CHANGE 1->1 TO 0->1 AND START BEING ACTIVE: KI(RNDI)-ACTKI(RNDI)
            ! #LINKS THAT CHANGE 1->1 TO 1->0 AND CEASES TO BE ACTIVE: ACTKI(RNDI)
            ACTKI(RNDI)=0
            DO J=1,KI(RNDI)
              IF (NI(POINTERS(RNDI,J)).EQ.1) THEN
                ACTKI(POINTERS(RNDI,J))=ACTKI(POINTERS(RNDI,J))+1
              ENDIF
            ENDDO
          ELSE ! INFECTION
            Z=INT(EACT*GENRAND_REAL2())+1
            COUNTER=0
            DO I=1,NINF ! LOOP OVER THE NODES THAT CAN INFECT
              COUNTER=COUNTER+ACTKI(INFLIST(I))
              IF (Z.LE.COUNTER) THEN
                COUNTER2=0
                DO J=1,KI(INFLIST(I)) ! LOOP OVER THE NODES THAT CAN GET INFECTED, GIVEN THE ONE THAT INFECTS
                  IF (NI(POINTERS(INFLIST(I),J)).EQ.0) THEN
                    IF (COUNTER2.EQ.(COUNTER-Z)) THEN
                      RNDI=POINTERS(INFLIST(I),J) ! NODE THAT GETS INFECTED
                      GO TO 10
                    ELSE
                      COUNTER2=COUNTER2+1
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
 10         CONTINUE
            NI(RNDI)=1
            NINF=NINF+1
            INFLIST(NINF)=RNDI
            DO J=1,KI(RNDI)
              IF (NI(POINTERS(RNDI,J)).EQ.0) THEN
                ACTKI(RNDI)=ACTKI(RNDI)+1
              ELSE
                ACTKI(POINTERS(RNDI,J))=ACTKI(POINTERS(RNDI,J))-1
              ENDIF
            ENDDO
            EACT=EACT+2*ACTKI(RNDI)-KI(RNDI)
            ! #LINKS THAT CHANGE 0->1 TO 1->1 AND CEASES TO BE ACTIVE: KI(RNDI)-ACTKI(RNDI)
            ! #LINKS THAT CHANGE 1->0 TO 1->1 AND START BEING ACTIVE: ACTKI(RNDI)
          ENDIF
          IF ((NINF.EQ.0).AND.(EACT.EQ.0)) THEN ! EXIT THE LOOP IF THERE ARE NO INFECTED NODES
            GO TO 20
          ENDIF
          X=GENRAND_REAL2()
          T=T-DLOG(X)/(NINF+LAMBDA*EACT) ! TIME UPDATE FOLLOWING GILLESPIE ALGORITHM
          WRITE(200,*) T,DBLE(NINF)/DBLE(N)
        ENDDO
 20     CONTINUE
        TMAX=T
        T=0.D0
        DT=1.D-2/(0.5D0*(N+LAMBDA*E))
        WRITE(300,*) T,AVERAGE(RO1,N)
        DO K=1,INT(TMAX/DT) ! NUMERIC INTEGRATION OF THE MEAN FIELD APPROX.
          T=T+DT
          CALL RUNGEKUTTA4_SIS(T,DT,RO1,N,RO2,KI,POINTERS,LAMBDA)
          DO I=1,N
            RO1(I)=RO2(I)
          ENDDO
          IF (MOD(K,100).EQ.0) THEN
            WRITE(300,*) T,AVERAGE(RO1,N)
          ENDIF
        ENDDO
C               STATIONARY #INFECTED AS A FUNCTION OF LAMBDA
      ELSE IF (SECTION.EQ.2) THEN
        LMBDLIST=(/0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.5,2.0,3.0,4.0,
     &  5.0,6.0,8.0,10.0,12.0/) ! LIST OF LAMBDAS TO TRY
        DO L=1,SIZE(LMBDLIST) ! LOOP OVER LAMBDAS
          LAMBDA=LMBDLIST(L)
          ! INITIAL CONDITIONS, SAME AS ABOVE BUT FOR EACH LAMBDA
          NINF=0
          DO I=1,N
            X=GENRAND_REAL2()
            IF (X.LE.INFECTED0) THEN
              NI(I)=1
              NINF=NINF+1
              INFLIST(NINF)=I
            ELSE
              NI(I)=0
            ENDIF
          ENDDO
          EACT=0
          DO I=1,N
            ACTKI(I)=0
            IF (NI(I).EQ.1) THEN
              DO J=1,KI(I)
                IF (NI(POINTERS(I,J)).EQ.0) THEN
                  ACTKI(I)=ACTKI(I)+1
                  EACT=EACT+1
                ENDIF
              ENDDO
            ENDIF
          ENDDO
          T=0.D0
          SUM1=0.D0
          SUM2=0.D0
          DO K=1,NPAS ! LOOP OVER TIMESTEPS, SIMILAR CODE AS ABOVE
            Y=(NINF+LAMBDA*EACT)*GENRAND_REAL2()
            IF (Y.LT.DBLE(NINF)) THEN ! RECOVERY
              Z=INT(NINF*GENRAND_REAL2())+1
              RNDI=INFLIST(Z)
              NI(RNDI)=0
              DO I=Z,NINF
                INFLIST(I)=INFLIST(I+1)
              ENDDO
              NINF=NINF-1
              EACT=EACT+KI(RNDI)-2*ACTKI(RNDI)
              ACTKI(RNDI)=0
              DO J=1,KI(RNDI)
                IF (NI(POINTERS(RNDI,J)).EQ.1) THEN
                  ACTKI(POINTERS(RNDI,J))=ACTKI(POINTERS(RNDI,J))+1
                ENDIF
              ENDDO
            ELSE ! INFECTION
              Z=INT(EACT*GENRAND_REAL2())+1
              COUNTER=0
              DO I=1,NINF
                COUNTER=COUNTER+ACTKI(INFLIST(I))
                IF (Z.LE.COUNTER) THEN
                  COUNTER2=0
                  DO J=1,KI(INFLIST(I))
                    IF (NI(POINTERS(INFLIST(I),J)).EQ.0) THEN
                      IF (COUNTER2.EQ.(COUNTER-Z)) THEN
                        RNDI=POINTERS(INFLIST(I),J)
                        GO TO 30
                      ELSE
                        COUNTER2=COUNTER2+1
                      ENDIF
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
 30           CONTINUE
              NI(RNDI)=1
              NINF=NINF+1
              INFLIST(NINF)=RNDI
              DO J=1,KI(RNDI)
                IF (NI(POINTERS(RNDI,J)).EQ.0) THEN
                  ACTKI(RNDI)=ACTKI(RNDI)+1
                ELSE
                  ACTKI(POINTERS(RNDI,J))=ACTKI(POINTERS(RNDI,J))-1
                ENDIF
              ENDDO
              EACT=EACT+2*ACTKI(RNDI)-KI(RNDI)
            ENDIF
            IF ((NINF.EQ.0).AND.(EACT.EQ.0)) THEN
              GO TO 40
            ENDIF
            X=GENRAND_REAL2()
            T=T-DLOG(X)/(NINF+LAMBDA*EACT)
            IF (K.GE.(NPAS+1-KAV)) THEN ! LOOP OVER (MORE OR LESS) EQUILIBRIUM STATES
              SUM1=SUM1+DBLE(NINF)/N/KAV
              SUM2=SUM2+DBLE(NINF*NINF)/N/N/KAV ! FOR COMPUTING THE STANDARD DEVIATION
            ENDIF
          ENDDO
 40       CONTINUE
          IF (LAMBDA.GT.1.D0/MEANK) THEN ! VALUE OF THE EQUILIBRIM FRACTION OF INFECTED IN A FULLY MIXED NETWORK
            FULLYMIXEDRO=1.D0-1.D0/MEANK/LAMBDA
          ELSE
            FULLYMIXEDRO=0.D0
          ENDIF
          WRITE(400,*) LAMBDA,SUM1,DSQRT(SUM2-SUM1*SUM1),FULLYMIXEDRO
        ENDDO
      ENDIF
      END PROGRAM

C     FUNCTIONS AND SUBROUTINES
      
      REAL*8 FUNCTION AVERAGE(LIST,LEN)
      ! TAKES THE AVERAGE OF A LIST
      IMPLICIT NONE
      INTEGER*4 LEN,I
      REAL*8 LIST(LEN),SUM
      SUM=0.D0
      DO I=1,LEN
        SUM=SUM+LIST(I)
      ENDDO
      AVERAGE=SUM/DBLE(LEN)
      RETURN
      END FUNCTION

      SUBROUTINE RUNGEKUTTA4_SIS(T,DT,YYIN,NEQUS,YYOUT,KI,PNTRS,LMBD)
      ! MAKES ONE STEP OF THE RUNGE-KUTTA 4 METHOD
      IMPLICIT NONE
      REAL*8 YYIN(NEQUS),YYOUT(NEQUS),T,DT,LMBD
      REAL*8 K0(NEQUS),K1(NEQUS),K2(NEQUS),K3(NEQUS)
      INTEGER*4 NEQUS,I,KI(NEQUS),PNTRS(NEQUS,10)
      CALL DERIV_SIS(T,YYIN,K0,NEQUS,KI,PNTRS,LMBD)
      CALL DERIV_SIS(T+0.5D0*DT,YYIN+0.5D0*DT*K0,K1,NEQUS,KI,PNTRS,LMBD)
      CALL DERIV_SIS(T+0.5D0*DT,YYIN+0.5D0*DT*K1,K2,NEQUS,KI,PNTRS,LMBD)
      CALL DERIV_SIS(T,YYIN+DT*K2,K3,NEQUS,KI,PNTRS,LMBD)
      DO I=1,NEQUS
        YYOUT(I)=YYIN(I)+(DT/6.D0)*(K0(I)+2.D0*K1(I)+2.D0*K2(I)+K3(I))
      ENDDO
      END SUBROUTINE

      SUBROUTINE DERIV_SIS(T,YIN,DYOUT,NEQU,KI,POINTERS,LMBD)
      ! DERIVATIVES OF RO_I
      IMPLICIT NONE
      REAL*8 T,YIN(NEQU),DYOUT(NEQU),LMBD
      INTEGER*4 NEQU,I,J,KI(NEQU),POINTERS(NEQU,10)
      DO I=1,NEQU
        DYOUT(I)=0.D0
        DO J=1,KI(I)
          DYOUT(I)=DYOUT(I)+YIN(POINTERS(I,J))
        ENDDO
        DYOUT(I)=(1.D0-YIN(I))*LMBD*DYOUT(I)-YIN(I)
      ENDDO
      RETURN
      END SUBROUTINE

c     initialize mt(0:N-1) with a seed
      subroutine init_genrand(s)
      integer s
      integer N
      integer DONE
      integer ALLBIT_MASK
      parameter (N=624)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask1/ ALLBIT_MASK
c
      call mt_initln
      mt(0)=iand(s,ALLBIT_MASK)
      do 100 mti=1,N-1
        mt(mti)=1812433253*
     &          ieor(mt(mti-1),ishft(mt(mti-1),-30))+mti
        mt(mti)=iand(mt(mti),ALLBIT_MASK)
  100 continue
      initialized=DONE
c
      return
      end

c     generates a random number on [0,0xffffffff]-interval
      function genrand_int32()
      integer genrand_int32
      integer N,M
      integer DONE
      integer UPPER_MASK,LOWER_MASK,MATRIX_A
      integer T1_MASK,T2_MASK
      parameter (N=624)
      parameter (M=397)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      integer y,kk
      integer mag01(0:1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
c
      if(initialized.ne.DONE)then
        call init_genrand(21641)
      endif
c
      if(mti.ge.N)then
        do 100 kk=0,N-M-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
  100   continue
        do 200 kk=N-M,N-1-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
  200   continue
        y=ior(iand(mt(N-1),UPPER_MASK),iand(mt(0),LOWER_MASK))
        mt(kk)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti=0
      endif
c
      y=mt(mti)
      mti=mti+1
c
      y=ieor(y,ishft(y,-11))
      y=ieor(y,iand(ishft(y,7),T1_MASK))
      y=ieor(y,iand(ishft(y,15),T2_MASK))
      y=ieor(y,ishft(y,-18))
c
      genrand_int32=y
      return
      end

c     generates a random number on [0,1)-real-interval
      function genrand_real2()
      double precision genrand_real2,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real2=r/4294967296.d0
      return
      end

c     initialize large number (over 32-bit constant number)
      subroutine mt_initln
      integer ALLBIT_MASK
      integer TOPBIT_MASK
      integer UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      integer mag01(0:1)
      common /mt_mask1/ ALLBIT_MASK
      common /mt_mask2/ TOPBIT_MASK
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
CC    TOPBIT_MASK = Z'80000000'
CC    ALLBIT_MASK = Z'ffffffff'
CC    UPPER_MASK  = Z'80000000'
CC    LOWER_MASK  = Z'7fffffff'
CC    MATRIX_A    = Z'9908b0df'
CC    T1_MASK     = Z'9d2c5680'
CC    T2_MASK     = Z'efc60000'
      TOPBIT_MASK=1073741824
      TOPBIT_MASK=ishft(TOPBIT_MASK,1)
      ALLBIT_MASK=2147483647
      ALLBIT_MASK=ior(ALLBIT_MASK,TOPBIT_MASK)
      UPPER_MASK=TOPBIT_MASK
      LOWER_MASK=2147483647
      MATRIX_A=419999967
      MATRIX_A=ior(MATRIX_A,TOPBIT_MASK)
      T1_MASK=489444992
      T1_MASK=ior(T1_MASK,TOPBIT_MASK)
      T2_MASK=1875247104
      T2_MASK=ior(T2_MASK,TOPBIT_MASK)
      mag01(0)=0
      mag01(1)=MATRIX_A
      return
      end