PROGRAM Main
    !! Hi
    USE, INTRINSIC:: ISO_FORTRAN_ENV

    USE GetMeanValues
    USE GetMean_And_Pet

    IMPLICIT NONE
    INTEGER, PARAMETER:: rDef = REAL64

    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE ::  Q_Mean              ,&
                                                        Q_Exact             ,&
                                                        Q_Mean_Store        ,&
                                                        Q_Mean_And_Pet_Store, &
                                                        Q_Mean_And_Pet      ,&
                                                        Q_Peturbation

    REAL(KIND = rDef):: p_Exit  ,&
                        gam

    REAL(KIND = rDef), DIMENSION(:), ALLOCATABLE::  x           ,&
                                                    Jac_Curv

    INTEGER::   nFlow   ,&
                i       ,&
                j       ,&
                iMin    ,&
                iMax    ,&
                bMin    ,&
                bMax    ,&
                nL      ,&
                nD      ,&
                DS

    LOGICAL:: debug

    debug = .TRUE.
    debug = .FALSE.

    nFlow = 1

    CALL MeanValues(Q_Mean  = Q_Mean    ,&
                    Q_Exact = Q_Exact   ,&
                    nFlow   = nFlow     ,&
                    p_Exit  = p_Exit    ,&
                    nL      = nL        ,&
                    nD      = nD        ,&
                    DS      = DS)
    STOP

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    OPEN(14, FILE = '/Users/zaidhs/Documents/PhD/&
                &Third_CAA_C1_P1/FinalSol/Initialize_CAA.dat'&
                , FORM = 'FORMATTED', STATUS = 'OLD')
        READ(14, *)    p_Exit   ,&
                       iMax     ,&
                       iMin     ,&
                       bMin     ,&
                       bMax     ,&
                       gam      ,&
                       nD       ,&
                       nL       ,&
                       DS
    CLOSE(14)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ALLOCATE(   Q_Mean_Store(iMax, bMax)            ,&
                Q_Mean_And_Pet_Store(iMax, bMax)    ,&
                Q_Mean_And_Pet(iMax, bMax)          ,&
                Q_Peturbation(iMax, bMax)           ,&
                x(iMax)                             ,&
                Jac_Curv(iMax))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ALLOCATE(   Q_Mean(iMax, bMax)     ,&
                Q_Exact(iMax, bMax))

    OPEN(15, FILE = '/Users/zaidhs/Documents/PhD/&
                &Third_CAA_C1_P1/FinalSol/Mean_Flow_CAA.dat'&
                , FORM = 'FORMATTED', STATUS = 'OLD')
        DO i = iMin, iMax
            READ(15, *) Q_Mean(i, 1)               ,&
                        Q_Mean(i, 2)               ,&
                        Q_Mean(i, 3)               ,&
                        Q_Exact(i, 1)              ,&
                        Q_Exact(i, 2)              ,&
                        Q_Exact(i, 3)              ,&
                        x(i)                       ,&
                        Jac_Curv(i)
        END DO
    CLOSE(15)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    nL = nL
    nD = nD
    DS = DS

    DO i = iMin, iMax
        DO j = bMin, bMax
            Q_Mean_Store(i, j) = Q_Mean(i, j)
        END DO
    END DO

    nFlow = 2
    CALL Mean_And_Pet( Q_Mean_And_Pet   = Q_Mean_And_Pet    ,&
                       Q_Mean           = Q_Mean            ,&
                       Q_Peturbation    = Q_Peturbation     ,&
                       Q_Exact          = Q_Exact           ,&
                       nFlow            = nFlow             ,&
                       p_Exit           = p_Exit            ,&
                       nL               = nL                ,&
                       nD               = nD                ,&
                       DS               = DS)


END PROGRAM Main
