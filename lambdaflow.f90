MODULE LAM_MOD
CONTAINS
SUBROUTINE DLAM_CAL(LAM,dtLAM,expuISr,expdISr,expuISs,expdISs  &
        ,RqSu,dRqSu,RqRu,dRqRu,RqSd,dRqSd,RqRd,dRqRd           &
        ,tanhur,tanhdr,tanhus,tanhds,qxy2,qzt2,ss,sr,m2)
  USE CONST_SAVE, ONLY : PI
  USE TMU_MOD

  IMPLICIT NONE

  REAL*8,INTENT(IN) :: m2
  REAL*8,INTENT(IN) :: LAM(:)
  REAL*8,INTENT(IN) :: expuISr,expdISr,expuISs,expdISs
  REAL*8,INTENT(IN) :: RqSu, dRqSu, RqRu, dRqRu
  REAL*8,INTENT(IN) :: RqSd, dRqSd, RqRd, dRqRd
  REAL*8,INTENT(IN) :: qxy2,qzt2,ss,sr
  REAL*8,INTENT(IN) :: tanhur,tanhdr,tanhus,tanhds
  REAL*8,INTENT(OUT) :: dtLAM(:)

  dtLAM(1)=                                                               &
   (2*dRqrd*expdISr*quB*qxy2**2*(1 + Rqrd)**2*tanhdr**3*                        &
       (expuISs*(1 + Rqsu)*(-1 + tanhus**2)*                                    &
          (16*(3*lam(1) - 2*lam(3))*(15*lam(1) + 2*lam(3)) +                    &
            45*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                         &
            96*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) - 64*lam(8)**2)       &
+ 4*expdISs*(1 + Rqsd)*(-1 + tanhds**2)*                                        &
          (2736*lam(1)**2 + 1944*lam(2)**2 + 960*lam(1)*lam(3) +                &
            88*lam(3)**2 + 9*(lam(4)**2 + 4*lam(5)**2 + lam(5)*lam(6) +         &
               lam(6)**2 + (8*lam(5) + lam(6))*lam(7) + 4*lam(7)**2 +           &
               lam(4)*(lam(5) + 2*lam(6) + lam(7))) -                           &
            54*lam(2)*(3*(lam(4) + 3*lam(5) + lam(6) + 3*lam(7)) -              &
               16*lam(8)) - 48*(lam(4) + 2*lam(5) + lam(6) + 2*lam(7))*         &
             lam(8) + 88*lam(8)**2)) +                                          &
      2*dRqru*expuISr*qdB*qxy2**2*(1 + Rqru)**2*tanhur**3*                      &
       (expdISs*(1 + Rqsd)*(-1 + tanhds**2)*                                    &
          (16*(3*lam(1) - 2*lam(3))*(15*lam(1) + 2*lam(3)) +                    &
            45*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                         &
            96*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) - 64*lam(8)**2)       &
+ 4*expuISs*(1 + Rqsu)*(-1 + tanhus**2)*                                        &
          (2736*lam(1)**2 + 1944*lam(2)**2 + 960*lam(1)*lam(3) +                &
            88*lam(3)**2 + 9*(lam(4)**2 + 4*lam(5)**2 + lam(5)*lam(6) +         &
               lam(6)**2 + (8*lam(5) + lam(6))*lam(7) + 4*lam(7)**2 +           &
               lam(4)*(lam(5) + 2*lam(6) + lam(7))) -                           &
            54*lam(2)*(3*(lam(4) + 3*lam(5) + lam(6) + 3*lam(7)) -              &
               16*lam(8)) - 48*(lam(4) + 2*lam(5) + lam(6) + 2*lam(7))*         &
             lam(8) + 88*lam(8)**2)) -                                          &
      expdISr*qxy2*tanhdr**2*(dRqrd*quB*                                        &
          (-(expuISs*(1 + Rqsu)*                                                &
               (-2*qzt2*(1 + Rqrd)**2*tanhus +                                  &
                 qdB*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*(-1 + tanhus**2))*          &
               (16*(3*lam(1) - 2*lam(3))*(15*lam(1) + 2*lam(3)) +               &
                 45*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                    &
                 96*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -                &
                 64*lam(8)**2)) -                                               &
            4*expdISs*(1 + Rqsd)*                                               &
             (-2*qzt2*(1 + Rqrd)**2*tanhds +                                    &
               qdB*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*(-1 + tanhds**2))*            &
             (2736*lam(1)**2 + 1944*lam(2)**2 + 960*lam(1)*lam(3) +             &
               88*lam(3)**2 +                                                   &
               9*(lam(4)**2 + 4*lam(5)**2 + lam(5)*lam(6) + lam(6)**2 +         &
                  (8*lam(5) + lam(6))*lam(7) + 4*lam(7)**2 +                    &
                  lam(4)*(lam(5) + 2*lam(6) + lam(7))) -                        &
               54*lam(2)*(3*(lam(4) + 3*lam(5) + lam(6) + 3*lam(7)) -           &
                  16*lam(8)) -                                                  &
               48*(lam(4) + 2*lam(5) + lam(6) + 2*lam(7))*lam(8) +              &
               88*lam(8)**2)) +                                                 &
         (1 + Rqrd)*(dRqsu*expuISs*qdB*(-1 + tanhus)*(1 + tanhus)*              &
             (quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +                              &
               2*qxy2*(1 + Rqsu)**2*tanhus)*                                    &
             (64*(45*lam(1)**2 - 24*lam(1)*lam(3) + 14*lam(3)**2) +             &
               9*(5*lam(4)**2 + 29*lam(5)**2 + 5*lam(6)**2 +                    &
                  2*lam(4)*(7*lam(5) - 5*lam(6) - 7*lam(7)) +                   &
                  14*lam(6)*lam(7) + 29*lam(7)**2 -                             &
                  2*lam(5)*(7*lam(6) + 29*lam(7))) -                            &
               288*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +             &
               768*lam(9)**2 +                                                  &
               192*(lam(4) + 2*lam(5) - lam(6) - 2*lam(7) - 4*lam(9))*          &
                lam(10) + 224*lam(10)**2) +                                     &
            4*dRqsd*expdISs*quB*                                                &
             (qdB*(-1 + 2*qzt2*(1 + Rqsd)**2*ss) +                              &
               2*qxy2*(1 + Rqsd)**2*tanhds)*(-1 + tanhds**2)*                   &
             (64*(3*lam(1) + lam(3))**2 +                                       &
               (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))**2)))      &
- expuISr*qxy2*tanhur**2*(dRqru*qdB*                                            &
          (-(expdISs*(1 + Rqsd)*                                                &
               (-2*qzt2*(1 + Rqru)**2*tanhds +                                  &
                 quB*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*(-1 + tanhds**2))*          &
               (16*(3*lam(1) - 2*lam(3))*(15*lam(1) + 2*lam(3)) +               &
                 45*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                    &
                 96*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -                &
                 64*lam(8)**2)) -                                               &
            4*expuISs*(1 + Rqsu)*                                               &
             (-2*qzt2*(1 + Rqru)**2*tanhus +                                    &
               quB*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*(-1 + tanhus**2))*            &
             (2736*lam(1)**2 + 1944*lam(2)**2 + 960*lam(1)*lam(3) +             &
               88*lam(3)**2 +                                                   &
               9*(lam(4)**2 + 4*lam(5)**2 + lam(5)*lam(6) + lam(6)**2 +         &
                  (8*lam(5) + lam(6))*lam(7) + 4*lam(7)**2 +                    &
                  lam(4)*(lam(5) + 2*lam(6) + lam(7))) -                        &
               54*lam(2)*(3*(lam(4) + 3*lam(5) + lam(6) + 3*lam(7)) -           &
                  16*lam(8)) -                                                  &
               48*(lam(4) + 2*lam(5) + lam(6) + 2*lam(7))*lam(8) +              &
               88*lam(8)**2)) +                                                 &
         (1 + Rqru)*(dRqsd*expdISs*quB*                                         &
             (qdB*(-1 + 2*qzt2*(1 + Rqsd)**2*ss) +                              &
               2*qxy2*(1 + Rqsd)**2*tanhds)*(-1 + tanhds**2)*                   &
             (64*(45*lam(1)**2 - 24*lam(1)*lam(3) + 14*lam(3)**2) +             &
               9*(5*lam(4)**2 + 29*lam(5)**2 + 5*lam(6)**2 +                    &
                  2*lam(4)*(7*lam(5) - 5*lam(6) - 7*lam(7)) +                   &
                  14*lam(6)*lam(7) + 29*lam(7)**2 -                             &
                  2*lam(5)*(7*lam(6) + 29*lam(7))) -                            &
               288*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +             &
               768*lam(9)**2 +                                                  &
               192*(lam(4) + 2*lam(5) - lam(6) - 2*lam(7) - 4*lam(9))*          &
                lam(10) + 224*lam(10)**2) +                                     &
            4*dRqsu*expuISs*qdB*(-1 + tanhus)*(1 + tanhus)*                     &
             (quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +                              &
               2*qxy2*(1 + Rqsu)**2*tanhus)*                                    &
             (64*(3*lam(1) + lam(3))**2 +                                       &
               (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))**2)))      &
- expdISs*expuISr*quB*(-(dRqru*qdB*                                             &
            (qxy2*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                    &
               (16*(3*lam(1) - 2*lam(3))*(15*lam(1) + 2*lam(3)) +               &
                 45*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                    &
                 96*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -                &
                 64*lam(8)**2) +                                                &
              qzt2*(720*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*              &
                  lam(1)**2 -                                                   &
                 64*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                  &
                  lam(3)**2 +                                                   &
                 128*m2*(1 + Rqru)*sr*lam(3)*                                   &
                  (3*(lam(4) - lam(5) + lam(6) - lam(7)) + 4*lam(8)) +          &
                 (1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                     &
                  (3*(lam(4) - lam(5) + lam(6) - lam(7)) - 8*lam(8))*           &
                  (15*(lam(4) - lam(5) + lam(6) - lam(7)) + 8*lam(8)) +         &
                 96*lam(1)*(-4*(1 + Rqsd)*                                      &
                     (-1 + 2*qzt2*(1 + Rqru)**2*sr)*lam(3) +                    &
                    m2*(1 + Rqru)*sr*                                           &
                     (-15*(lam(4) - lam(5) + lam(6) - lam(7)) +                 &
                       16*lam(8)))))) -                                         &
         2*dRqsd*qxy2**2*(1 + Rqru)*(1 + Rqsd)**2*tanhds**3*                    &
          (64*(45*lam(1)**2 - 24*lam(1)*lam(3) + 14*lam(3)**2) +                &
            9*(5*lam(4)**2 + 29*lam(5)**2 + 5*lam(6)**2 +                       &
               2*lam(4)*(7*lam(5) - 5*lam(6) - 7*lam(7)) +                      &
               14*lam(6)*lam(7) + 29*lam(7)**2 -                                &
               2*lam(5)*(7*lam(6) + 29*lam(7))) -                               &
            288*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +                &
            768*lam(9)**2 + 192*                                                &
             (lam(4) + 2*lam(5) - lam(6) - 2*lam(7) - 4*lam(9))*lam(10)         &
+ 224*lam(10)**2) + 2*dRqsd*qxy2*(1 + Rqsd)*tanhds*                             &
          (16*m2*(9*lam(1)*(5*lam(4) + 7*lam(5) - 5*lam(6) -                    &
                  7*lam(7) - 16*lam(9)) +                                       &
               12*lam(3)*(-lam(4) - 5*lam(5) + lam(6) + 5*lam(7) +              &
                  8*lam(9)) + 8*(12*lam(1) - 5*lam(3))*lam(10)) +               &
            qxy2*(1 + Rqru)*(1 + Rqsd)*                                         &
             (64*(45*lam(1)**2 - 24*lam(1)*lam(3) + 14*lam(3)**2) +             &
               9*(5*lam(4)**2 + 29*lam(5)**2 + 5*lam(6)**2 +                    &
                  2*lam(4)*(7*lam(5) - 5*lam(6) - 7*lam(7)) +                   &
                  14*lam(6)*lam(7) + 29*lam(7)**2 -                             &
                  2*lam(5)*(7*lam(6) + 29*lam(7))) -                            &
               288*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +             &
               768*lam(9)**2 +                                                  &
               192*(lam(4) + 2*lam(5) - lam(6) - 2*lam(7) - 4*lam(9))*          &
                lam(10) + 224*lam(10)**2) +                                     &
            qzt2*(1 + Rqru)*(1 + Rqsd)*                                         &
             (64*(45*lam(1)**2 - 24*lam(1)*lam(3) + 14*lam(3)**2) +             &
               9*(5*lam(4)**2 + 29*lam(5)**2 + 5*lam(6)**2 +                    &
                  2*lam(4)*(7*lam(5) - 5*lam(6) - 7*lam(7)) +                   &
                  14*lam(6)*lam(7) + 29*lam(7)**2 -                             &
                  2*lam(5)*(7*lam(6) + 29*lam(7))) -                            &
               288*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +             &
               768*lam(9)**2 +                                                  &
               192*(lam(4) + 2*lam(5) - lam(6) - 2*lam(7) - 4*lam(9))*          &
                lam(10) + 224*lam(10)**2)) +                                    &
         qdB*qxy2*tanhds**2*(dRqru*(1 + Rqsd)*                                  &
             (-1 + 2*qzt2*(1 + Rqru)**2*sr)*                                    &
             (16*(3*lam(1) - 2*lam(3))*(15*lam(1) + 2*lam(3)) +                 &
               45*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                      &
               96*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -                  &
               64*lam(8)**2) -                                                  &
            dRqsd*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                    &
             (64*(45*lam(1)**2 - 24*lam(1)*lam(3) + 14*lam(3)**2) +             &
               9*(5*lam(4)**2 + 29*lam(5)**2 + 5*lam(6)**2 +                    &
                  2*lam(4)*(7*lam(5) - 5*lam(6) - 7*lam(7)) +                   &
                  14*lam(6)*lam(7) + 29*lam(7)**2 -                             &
                  2*lam(5)*(7*lam(6) + 29*lam(7))) -                            &
               288*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +             &
               768*lam(9)**2 +                                                  &
               192*(lam(4) + 2*lam(5) - lam(6) - 2*lam(7) - 4*lam(9))*          &
                lam(10) + 224*lam(10)**2)) +                                    &
         dRqsd*qdB*(qxy2*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*             &
             (64*(45*lam(1)**2 - 24*lam(1)*lam(3) + 14*lam(3)**2) +             &
               9*(5*lam(4)**2 + 29*lam(5)**2 + 5*lam(6)**2 +                    &
                  2*lam(4)*(7*lam(5) - 5*lam(6) - 7*lam(7)) +                   &
                  14*lam(6)*lam(7) + 29*lam(7)**2 -                             &
                  2*lam(5)*(7*lam(6) + 29*lam(7))) -                            &
               288*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +             &
               768*lam(9)**2 +                                                  &
               192*(lam(4) + 2*lam(5) - lam(6) - 2*lam(7) - 4*lam(9))*          &
                lam(10) + 224*lam(10)**2) +                                     &
            qzt2*(2880*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*               &
                lam(1)**2 + 896*(1 + Rqru)*                                     &
                (-1 + 2*qzt2*(1 + Rqsd)**2*ss)*lam(3)**2 -                      &
               128*m2*(1 + Rqsd)*ss*lam(3)*                                     &
                (3*lam(4) + 15*lam(5) -                                         &
                  3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10)) +              &
               (1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                       &
                (9*(5*lam(4)**2 + 29*lam(5)**2 + 5*lam(6)**2 +                  &
                     2*lam(4)*(7*lam(5) - 5*lam(6) - 7*lam(7)) +                &
                     14*lam(6)*lam(7) + 29*lam(7)**2 -                          &
                     2*lam(5)*(7*lam(6) + 29*lam(7))) -                         &
                  288*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +          &
                  768*lam(9)**2 +                                               &
                  192*(lam(4) + 2*lam(5) - lam(6) - 2*lam(7) -                  &
                     4*lam(9))*lam(10) + 224*lam(10)**2) +                      &
               96*lam(1)*(-16*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*        &
                   lam(3) + m2*(1 + Rqsd)*ss*                                   &
                   (15*lam(4) + 21*lam(5) -                                     &
                     3*(5*lam(6) + 7*lam(7) + 16*lam(9)) + 32*lam(10))))))      &
- expdISr*tanhdr*(-2*dRqrd*quB*qxy2*(1 + Rqrd)*                                 &
          (expuISs*(-48*m2*(15*lam(1) - 4*lam(3))*                              &
                (lam(4) - lam(5) + lam(6) - lam(7)) +                           &
               256*m2*(3*lam(1) + lam(3))*lam(8) +                              &
               qxy2*(1 + Rqrd)*(1 + Rqsu)*                                      &
                (16*(3*lam(1) - 2*lam(3))*(15*lam(1) + 2*lam(3)) +              &
                  45*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                   &
                  96*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -               &
                  64*lam(8)**2) +                                               &
               qzt2*(1 + Rqrd)*(1 + Rqsu)*                                      &
                (16*(3*lam(1) - 2*lam(3))*(15*lam(1) + 2*lam(3)) +              &
                  45*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                   &
                  96*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -               &
                  64*lam(8)**2)) +                                              &
            4*expdISs*(qxy2*(1 + Rqrd)*(1 + Rqsd)*                              &
                (2736*lam(1)**2 + 1944*lam(2)**2 + 960*lam(1)*lam(3) +          &
                  88*lam(3)**2 +                                                &
                  9*(lam(4)**2 + 4*lam(5)**2 + lam(5)*lam(6) +                  &
                     lam(6)**2 + (8*lam(5) + lam(6))*lam(7) +                   &
                     4*lam(7)**2 + lam(4)*(lam(5) + 2*lam(6) + lam(7)))         &
- 54*lam(2)*(3*(lam(4) + 3*lam(5) + lam(6) + 3*lam(7)) - 16*lam(8)) -           &
                  48*(lam(4) + 2*lam(5) + lam(6) + 2*lam(7))*lam(8) +           &
                  88*lam(8)**2) +                                               &
               qzt2*(1 + Rqrd)*(1 + Rqsd)*                                      &
                (2736*lam(1)**2 + 1944*lam(2)**2 + 960*lam(1)*lam(3) +          &
                  88*lam(3)**2 +                                                &
                  9*(lam(4)**2 + 4*lam(5)**2 + lam(5)*lam(6) +                  &
                     lam(6)**2 + (8*lam(5) + lam(6))*lam(7) +                   &
                     4*lam(7)**2 + lam(4)*(lam(5) + 2*lam(6) + lam(7)))         &
- 54*lam(2)*(3*(lam(4) + 3*lam(5) + lam(6) + 3*lam(7)) - 16*lam(8)) -           &
                  48*(lam(4) + 2*lam(5) + lam(6) + 2*lam(7))*lam(8) +           &
                  88*lam(8)**2) +                                               &
               4*m2*(4*lam(3)*                                                  &
                   (108*lam(2) -                                                &
                     3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +                &
                     26*lam(8)) +                                               &
                  3*lam(1)*(756*lam(2) - 39*lam(4) - 87*lam(5) -                &
                     39*lam(6) - 87*lam(7) + 176*lam(8))))) -                   &
         2*expuISs*qxy2*(1 + Rqsu)*tanhus**2*                                   &
          (-(dRqrd*quB*qxy2*(1 + Rqrd)**2*                                      &
               (16*(3*lam(1) - 2*lam(3))*(15*lam(1) + 2*lam(3)) +               &
                 45*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                    &
                 96*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -                &
                 64*lam(8)**2)) +                                               &
            dRqsu*qdB*(16*m2*                                                   &
                (9*lam(1)*(5*lam(4) + 7*lam(5) - 5*lam(6) - 7*lam(7) -          &
                     16*lam(9)) +                                               &
                  12*lam(3)*(-lam(4) - 5*lam(5) + lam(6) + 5*lam(7) +           &
                     8*lam(9)) + 8*(12*lam(1) - 5*lam(3))*lam(10)) +            &
               qzt2*(1 + Rqrd)*(1 + Rqsu)*                                      &
                (64*(45*lam(1)**2 - 24*lam(1)*lam(3) + 14*lam(3)**2) +          &
                  9*(5*lam(4)**2 + 29*lam(5)**2 + 5*lam(6)**2 +                 &
                     2*lam(4)*(7*lam(5) - 5*lam(6) - 7*lam(7)) +                &
                     14*lam(6)*lam(7) + 29*lam(7)**2 -                          &
                     2*lam(5)*(7*lam(6) + 29*lam(7))) -                         &
                  288*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +          &
                  768*lam(9)**2 +                                               &
                  192*(lam(4) + 2*lam(5) - lam(6) - 2*lam(7) -                  &
                     4*lam(9))*lam(10) + 224*lam(10)**2))) +                    &
         4*expdISs*qdB*quB*qzt2*tanhds*                                         &
          (dRqrd*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                     &
             (2736*lam(1)**2 + 1944*lam(2)**2 + 960*lam(1)*lam(3) +             &
               88*lam(3)**2 +                                                   &
               9*(lam(4)**2 + 4*lam(5)**2 + lam(5)*lam(6) +                     &
                  lam(6)**2 + (8*lam(5) + lam(6))*lam(7) +                      &
                  4*lam(7)**2 + lam(4)*(lam(5) + 2*lam(6) + lam(7))) -          &
               54*lam(2)*(3*(lam(4) + 3*lam(5) + lam(6) + 3*lam(7)) -           &
                  16*lam(8)) -                                                  &
               48*(lam(4) + 2*lam(5) + lam(6) + 2*lam(7))*lam(8) +              &
               88*lam(8)**2) -                                                  &
            dRqsd*(576*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*               &
                lam(1)**2 + 64*(1 + Rqrd)*                                      &
                (-1 + 2*qzt2*(1 + Rqsd)**2*ss)*lam(3)**2 +                      &
               32*m2*(1 + Rqsd)*ss*lam(3)*                                      &
                (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)) +           &
               (1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                       &
                (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10))**2 +        &
               96*lam(1)*(4*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*          &
                   lam(3) + m2*(1 + Rqsd)*ss*                                   &
                   (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)))))       &
- 8*expdISs*quB*qxy2*(1 + Rqsd)*tanhds**2*                                      &
          (-(dRqrd*qxy2*(1 + Rqrd)**2*                                          &
               (2736*lam(1)**2 + 1944*lam(2)**2 + 960*lam(1)*lam(3) +           &
                 88*lam(3)**2 +                                                 &
                 9*(lam(4)**2 + 4*lam(5)**2 + lam(5)*lam(6) +                   &
                    lam(6)**2 + (8*lam(5) + lam(6))*lam(7) +                    &
                    4*lam(7)**2 + lam(4)*(lam(5) + 2*lam(6) + lam(7)))          &
- 54*lam(2)*(3*(lam(4) + 3*lam(5) + lam(6) + 3*lam(7)) - 16*lam(8)) -           &
                 48*(lam(4) + 2*lam(5) + lam(6) + 2*lam(7))*lam(8) +            &
                 88*lam(8)**2)) +                                               &
            dRqsd*(16*m2*(3*lam(1) + lam(3))*                                   &
                (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)) +           &
               qzt2*(1 + Rqrd)*(1 + Rqsd)*                                      &
                (64*(3*lam(1) + lam(3))**2 +                                    &
                  (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10))**2)       &
)) + expuISs*qdB*quB*qzt2*tanhus*                                               &
          (dRqrd*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                     &
             (16*(3*lam(1) - 2*lam(3))*(15*lam(1) + 2*lam(3)) +                 &
               45*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                      &
               96*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -                  &
               64*lam(8)**2) -                                                  &
            dRqsu*(2880*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*              &
                lam(1)**2 + 896*(1 + Rqrd)*                                     &
                (-1 + 2*qzt2*(1 + Rqsu)**2*ss)*lam(3)**2 -                      &
               128*m2*(1 + Rqsu)*ss*lam(3)*                                     &
                (3*lam(4) + 15*lam(5) -                                         &
                  3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10)) +              &
               (1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                       &
                (9*(5*lam(4)**2 + 29*lam(5)**2 + 5*lam(6)**2 +                  &
                     2*lam(4)*(7*lam(5) - 5*lam(6) - 7*lam(7)) +                &
                     14*lam(6)*lam(7) + 29*lam(7)**2 -                          &
                     2*lam(5)*(7*lam(6) + 29*lam(7))) -                         &
                  288*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +          &
                  768*lam(9)**2 +                                               &
                  192*(lam(4) + 2*lam(5) - lam(6) - 2*lam(7) -                  &
                     4*lam(9))*lam(10) + 224*lam(10)**2) +                      &
               96*lam(1)*(-16*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*        &
                   lam(3) + m2*(1 + Rqsu)*ss*                                   &
                   (15*lam(4) + 21*lam(5) -                                     &
                     3*(5*lam(6) + 7*lam(7) + 16*lam(9)) + 32*lam(10))))))      &
- expdISr*expuISs*qdB*(-(dRqrd*quB*                                             &
            (qxy2*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                    &
               (16*(3*lam(1) - 2*lam(3))*(15*lam(1) + 2*lam(3)) +               &
                 45*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                    &
                 96*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -                &
                 64*lam(8)**2) +                                                &
              qzt2*(720*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*              &
                  lam(1)**2 -                                                   &
                 64*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                  &
                  lam(3)**2 +                                                   &
                 128*m2*(1 + Rqrd)*sr*lam(3)*                                   &
                  (3*(lam(4) - lam(5) + lam(6) - lam(7)) + 4*lam(8)) +          &
                 (1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                     &
                  (3*(lam(4) - lam(5) + lam(6) - lam(7)) - 8*lam(8))*           &
                  (15*(lam(4) - lam(5) + lam(6) - lam(7)) + 8*lam(8)) +         &
                 96*lam(1)*(-4*(1 + Rqsu)*                                      &
                     (-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(3) +                    &
                    m2*(1 + Rqrd)*sr*                                           &
                     (-15*(lam(4) - lam(5) + lam(6) - lam(7)) +                 &
                       16*lam(8)))))) -                                         &
         2*dRqsu*qxy2**2*(1 + Rqrd)*(1 + Rqsu)**2*tanhus**3*                    &
          (64*(45*lam(1)**2 - 24*lam(1)*lam(3) + 14*lam(3)**2) +                &
            9*(5*lam(4)**2 + 29*lam(5)**2 + 5*lam(6)**2 +                       &
               2*lam(4)*(7*lam(5) - 5*lam(6) - 7*lam(7)) +                      &
               14*lam(6)*lam(7) + 29*lam(7)**2 -                                &
               2*lam(5)*(7*lam(6) + 29*lam(7))) -                               &
            288*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +                &
            768*lam(9)**2 + 192*                                                &
             (lam(4) + 2*lam(5) - lam(6) - 2*lam(7) - 4*lam(9))*lam(10)         &
+ 224*lam(10)**2) + 2*dRqsu*qxy2*(1 + Rqsu)*tanhus*                             &
          (16*m2*(9*lam(1)*(5*lam(4) + 7*lam(5) - 5*lam(6) -                    &
                  7*lam(7) - 16*lam(9)) +                                       &
               12*lam(3)*(-lam(4) - 5*lam(5) + lam(6) + 5*lam(7) +              &
                  8*lam(9)) + 8*(12*lam(1) - 5*lam(3))*lam(10)) +               &
            qxy2*(1 + Rqrd)*(1 + Rqsu)*                                         &
             (64*(45*lam(1)**2 - 24*lam(1)*lam(3) + 14*lam(3)**2) +             &
               9*(5*lam(4)**2 + 29*lam(5)**2 + 5*lam(6)**2 +                    &
                  2*lam(4)*(7*lam(5) - 5*lam(6) - 7*lam(7)) +                   &
                  14*lam(6)*lam(7) + 29*lam(7)**2 -                             &
                  2*lam(5)*(7*lam(6) + 29*lam(7))) -                            &
               288*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +             &
               768*lam(9)**2 +                                                  &
               192*(lam(4) + 2*lam(5) - lam(6) - 2*lam(7) - 4*lam(9))*          &
                lam(10) + 224*lam(10)**2) +                                     &
            qzt2*(1 + Rqrd)*(1 + Rqsu)*                                         &
             (64*(45*lam(1)**2 - 24*lam(1)*lam(3) + 14*lam(3)**2) +             &
               9*(5*lam(4)**2 + 29*lam(5)**2 + 5*lam(6)**2 +                    &
                  2*lam(4)*(7*lam(5) - 5*lam(6) - 7*lam(7)) +                   &
                  14*lam(6)*lam(7) + 29*lam(7)**2 -                             &
                  2*lam(5)*(7*lam(6) + 29*lam(7))) -                            &
               288*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +             &
               768*lam(9)**2 +                                                  &
               192*(lam(4) + 2*lam(5) - lam(6) - 2*lam(7) - 4*lam(9))*          &
                lam(10) + 224*lam(10)**2)) +                                    &
         quB*qxy2*tanhus**2*(dRqrd*(1 + Rqsu)*                                  &
             (-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                                    &
             (16*(3*lam(1) - 2*lam(3))*(15*lam(1) + 2*lam(3)) +                 &
               45*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                      &
               96*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -                  &
               64*lam(8)**2) -                                                  &
            dRqsu*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                    &
             (64*(45*lam(1)**2 - 24*lam(1)*lam(3) + 14*lam(3)**2) +             &
               9*(5*lam(4)**2 + 29*lam(5)**2 + 5*lam(6)**2 +                    &
                  2*lam(4)*(7*lam(5) - 5*lam(6) - 7*lam(7)) +                   &
                  14*lam(6)*lam(7) + 29*lam(7)**2 -                             &
                  2*lam(5)*(7*lam(6) + 29*lam(7))) -                            &
               288*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +             &
               768*lam(9)**2 +                                                  &
               192*(lam(4) + 2*lam(5) - lam(6) - 2*lam(7) - 4*lam(9))*          &
                lam(10) + 224*lam(10)**2)) +                                    &
         dRqsu*quB*(qxy2*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*             &
             (64*(45*lam(1)**2 - 24*lam(1)*lam(3) + 14*lam(3)**2) +             &
               9*(5*lam(4)**2 + 29*lam(5)**2 + 5*lam(6)**2 +                    &
                  2*lam(4)*(7*lam(5) - 5*lam(6) - 7*lam(7)) +                   &
                  14*lam(6)*lam(7) + 29*lam(7)**2 -                             &
                  2*lam(5)*(7*lam(6) + 29*lam(7))) -                            &
               288*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +             &
               768*lam(9)**2 +                                                  &
               192*(lam(4) + 2*lam(5) - lam(6) - 2*lam(7) - 4*lam(9))*          &
                lam(10) + 224*lam(10)**2) +                                     &
            qzt2*(2880*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*               &
                lam(1)**2 + 896*(1 + Rqrd)*                                     &
                (-1 + 2*qzt2*(1 + Rqsu)**2*ss)*lam(3)**2 -                      &
               128*m2*(1 + Rqsu)*ss*lam(3)*                                     &
                (3*lam(4) + 15*lam(5) -                                         &
                  3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10)) +              &
               (1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                       &
                (9*(5*lam(4)**2 + 29*lam(5)**2 + 5*lam(6)**2 +                  &
                     2*lam(4)*(7*lam(5) - 5*lam(6) - 7*lam(7)) +                &
                     14*lam(6)*lam(7) + 29*lam(7)**2 -                          &
                     2*lam(5)*(7*lam(6) + 29*lam(7))) -                         &
                  288*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +          &
                  768*lam(9)**2 +                                               &
                  192*(lam(4) + 2*lam(5) - lam(6) - 2*lam(7) -                  &
                     4*lam(9))*lam(10) + 224*lam(10)**2) +                      &
               96*lam(1)*(-16*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*        &
                   lam(3) + m2*(1 + Rqsu)*ss*                                   &
                   (15*lam(4) + 21*lam(5) -                                     &
                     3*(5*lam(6) + 7*lam(7) + 16*lam(9)) + 32*lam(10))))))      &
+ 4*expdISr*expdISs*quB*(2*dRqsd*qxy2**2*(1 + Rqrd)*(1 + Rqsd)**2*              &
          tanhds**3*(64*(3*lam(1) + lam(3))**2 +                                &
            (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))**2) -         &
         qdB*qxy2*tanhds**2*(dRqrd*(1 + Rqsd)*                                  &
             (-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                                    &
             (2736*lam(1)**2 + 1944*lam(2)**2 + 960*lam(1)*lam(3) +             &
               88*lam(3)**2 +                                                   &
               9*(lam(4)**2 + 4*lam(5)**2 + lam(5)*lam(6) +                     &
                  lam(6)**2 + (8*lam(5) + lam(6))*lam(7) +                      &
                  4*lam(7)**2 + lam(4)*(lam(5) + 2*lam(6) + lam(7))) -          &
               54*lam(2)*(3*(lam(4) + 3*lam(5) + lam(6) + 3*lam(7)) -           &
                  16*lam(8)) -                                                  &
               48*(lam(4) + 2*lam(5) + lam(6) + 2*lam(7))*lam(8) +              &
               88*lam(8)**2) -                                                  &
            dRqsd*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                    &
             (64*(3*lam(1) + lam(3))**2 +                                       &
               (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))**2))       &
- 2*dRqsd*qxy2*(1 + Rqsd)*tanhds*                                               &
          (16*m2*(3*lam(1) + lam(3))*                                           &
             (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)) +              &
            qxy2*(1 + Rqrd)*(1 + Rqsd)*                                         &
             (64*(3*lam(1) + lam(3))**2 +                                       &
               (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))**2)        &
+ qzt2*(1 + Rqrd)*(1 + Rqsd)*(64*(3*lam(1) + lam(3))**2 +                       &
               (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10))**2)) +       &
         qdB*(dRqrd*(qxy2*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*            &
                (2736*lam(1)**2 + 1944*lam(2)**2 + 960*lam(1)*lam(3) +          &
                  88*lam(3)**2 +                                                &
                  9*(lam(4)**2 + 4*lam(5)**2 + lam(5)*lam(6) +                  &
                     lam(6)**2 + (8*lam(5) + lam(6))*lam(7) +                   &
                     4*lam(7)**2 + lam(4)*(lam(5) + 2*lam(6) + lam(7)))         &
- 54*lam(2)*(3*(lam(4) + 3*lam(5) + lam(6) + 3*lam(7)) - 16*lam(8)) -           &
                  48*(lam(4) + 2*lam(5) + lam(6) + 2*lam(7))*lam(8) +           &
                  88*lam(8)**2) +                                               &
               qzt2*(2736*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*            &
                   lam(1)**2 +                                                  &
                  1944*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*               &
                   lam(2)**2 - 88*lam(3)**2 - 88*Rqsd*lam(3)**2 +               &
                  176*qzt2*sr*lam(3)**2 + 352*qzt2*Rqrd*sr*lam(3)**2 +          &
                  176*qzt2*Rqrd**2*sr*lam(3)**2 +                               &
                  176*qzt2*Rqsd*sr*lam(3)**2 +                                  &
                  352*qzt2*Rqrd*Rqsd*sr*lam(3)**2 +                             &
                  176*qzt2*Rqrd**2*Rqsd*sr*lam(3)**2 -                          &
                  96*m2*sr*lam(3)*lam(4) - 96*m2*Rqrd*sr*lam(3)*lam(4) -        &
                  9*lam(4)**2 - 9*Rqsd*lam(4)**2 +                              &
                  18*qzt2*sr*lam(4)**2 + 36*qzt2*Rqrd*sr*lam(4)**2 +            &
                  18*qzt2*Rqrd**2*sr*lam(4)**2 +                                &
                  18*qzt2*Rqsd*sr*lam(4)**2 +                                   &
                  36*qzt2*Rqrd*Rqsd*sr*lam(4)**2 +                              &
                  18*qzt2*Rqrd**2*Rqsd*sr*lam(4)**2 -                           &
                  480*m2*sr*lam(3)*lam(5) -                                     &
                  480*m2*Rqrd*sr*lam(3)*lam(5) - 9*lam(4)*lam(5) -              &
                  9*Rqsd*lam(4)*lam(5) + 18*qzt2*sr*lam(4)*lam(5) +             &
                  36*qzt2*Rqrd*sr*lam(4)*lam(5) +                               &
                  18*qzt2*Rqrd**2*sr*lam(4)*lam(5) +                            &
                  18*qzt2*Rqsd*sr*lam(4)*lam(5) +                               &
                  36*qzt2*Rqrd*Rqsd*sr*lam(4)*lam(5) +                          &
                  18*qzt2*Rqrd**2*Rqsd*sr*lam(4)*lam(5) - 36*lam(5)**2 -        &
                  36*Rqsd*lam(5)**2 + 72*qzt2*sr*lam(5)**2 +                    &
                  144*qzt2*Rqrd*sr*lam(5)**2 +                                  &
                  72*qzt2*Rqrd**2*sr*lam(5)**2 +                                &
                  72*qzt2*Rqsd*sr*lam(5)**2 +                                   &
                  144*qzt2*Rqrd*Rqsd*sr*lam(5)**2 +                             &
                  72*qzt2*Rqrd**2*Rqsd*sr*lam(5)**2 -                           &
                  96*m2*sr*lam(3)*lam(6) - 96*m2*Rqrd*sr*lam(3)*lam(6) -        &
                  18*lam(4)*lam(6) - 18*Rqsd*lam(4)*lam(6) +                    &
                  36*qzt2*sr*lam(4)*lam(6) +                                    &
                  72*qzt2*Rqrd*sr*lam(4)*lam(6) +                               &
                  36*qzt2*Rqrd**2*sr*lam(4)*lam(6) +                            &
                  36*qzt2*Rqsd*sr*lam(4)*lam(6) +                               &
                  72*qzt2*Rqrd*Rqsd*sr*lam(4)*lam(6) +                          &
                  36*qzt2*Rqrd**2*Rqsd*sr*lam(4)*lam(6) -                       &
                  9*lam(5)*lam(6) - 9*Rqsd*lam(5)*lam(6) +                      &
                  18*qzt2*sr*lam(5)*lam(6) +                                    &
                  36*qzt2*Rqrd*sr*lam(5)*lam(6) +                               &
                  18*qzt2*Rqrd**2*sr*lam(5)*lam(6) +                            &
                  18*qzt2*Rqsd*sr*lam(5)*lam(6) +                               &
                  36*qzt2*Rqrd*Rqsd*sr*lam(5)*lam(6) +                          &
                  18*qzt2*Rqrd**2*Rqsd*sr*lam(5)*lam(6) - 9*lam(6)**2 -         &
                  9*Rqsd*lam(6)**2 + 18*qzt2*sr*lam(6)**2 +                     &
                  36*qzt2*Rqrd*sr*lam(6)**2 +                                   &
                  18*qzt2*Rqrd**2*sr*lam(6)**2 +                                &
                  18*qzt2*Rqsd*sr*lam(6)**2 +                                   &
                  36*qzt2*Rqrd*Rqsd*sr*lam(6)**2 +                              &
                  18*qzt2*Rqrd**2*Rqsd*sr*lam(6)**2 -                           &
                  480*m2*sr*lam(3)*lam(7) -                                     &
                  480*m2*Rqrd*sr*lam(3)*lam(7) - 9*lam(4)*lam(7) -              &
                  9*Rqsd*lam(4)*lam(7) + 18*qzt2*sr*lam(4)*lam(7) +             &
                  36*qzt2*Rqrd*sr*lam(4)*lam(7) +                               &
                  18*qzt2*Rqrd**2*sr*lam(4)*lam(7) +                            &
                  18*qzt2*Rqsd*sr*lam(4)*lam(7) +                               &
                  36*qzt2*Rqrd*Rqsd*sr*lam(4)*lam(7) +                          &
                  18*qzt2*Rqrd**2*Rqsd*sr*lam(4)*lam(7) -                       &
                  72*lam(5)*lam(7) - 72*Rqsd*lam(5)*lam(7) +                    &
                  144*qzt2*sr*lam(5)*lam(7) +                                   &
                  288*qzt2*Rqrd*sr*lam(5)*lam(7) +                              &
                  144*qzt2*Rqrd**2*sr*lam(5)*lam(7) +                           &
                  144*qzt2*Rqsd*sr*lam(5)*lam(7) +                              &
                  288*qzt2*Rqrd*Rqsd*sr*lam(5)*lam(7) +                         &
                  144*qzt2*Rqrd**2*Rqsd*sr*lam(5)*lam(7) -                      &
                  9*lam(6)*lam(7) - 9*Rqsd*lam(6)*lam(7) +                      &
                  18*qzt2*sr*lam(6)*lam(7) +                                    &
                  36*qzt2*Rqrd*sr*lam(6)*lam(7) +                               &
                  18*qzt2*Rqrd**2*sr*lam(6)*lam(7) +                            &
                  18*qzt2*Rqsd*sr*lam(6)*lam(7) +                               &
                  36*qzt2*Rqrd*Rqsd*sr*lam(6)*lam(7) +                          &
                  18*qzt2*Rqrd**2*Rqsd*sr*lam(6)*lam(7) - 36*lam(7)**2 -        &
                  36*Rqsd*lam(7)**2 + 72*qzt2*sr*lam(7)**2 +                    &
                  144*qzt2*Rqrd*sr*lam(7)**2 +                                  &
                  72*qzt2*Rqrd**2*sr*lam(7)**2 +                                &
                  72*qzt2*Rqsd*sr*lam(7)**2 +                                   &
                  144*qzt2*Rqrd*Rqsd*sr*lam(7)**2 +                             &
                  72*qzt2*Rqrd**2*Rqsd*sr*lam(7)**2 +                           &
                  54*lam(2)*(64*m2*(1 + Rqrd)*sr*lam(3) -                       &
                     (1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                 &
                      (3*(lam(4) + 3*lam(5) + lam(6) + 3*lam(7)) -              &
                        16*lam(8))) +                                           &
                  16*(52*m2*(1 + Rqrd)*sr*lam(3) -                              &
                     3*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*               &
                      (lam(4) + 2*lam(5) + lam(6) + 2*lam(7)))*lam(8) +         &
                  88*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                 &
                   lam(8)**2 +                                                  &
                  24*lam(1)*(40*(1 + Rqsd)*                                     &
                      (-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(3) +                   &
                     m2*(1 + Rqrd)*sr*                                          &
                      (756*lam(2) - 39*lam(4) - 87*lam(5) - 39*lam(6) -         &
                        87*lam(7) + 176*lam(8))))) -                            &
            dRqsd*(qxy2*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*              &
                (64*(3*lam(1) + lam(3))**2 +                                    &
                  (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))**       &
                   2) + qzt2*(576*(1 + Rqrd)*                                   &
                   (-1 + 2*qzt2*(1 + Rqsd)**2*ss)*lam(1)**2 +                   &
                  64*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                 &
                   lam(3)**2 +                                                  &
                  32*m2*(1 + Rqsd)*ss*lam(3)*                                   &
                   (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)) +        &
                  (1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                    &
                   (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10))**2       &
+ 96*lam(1)*(4*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*lam(3) +               &
                     m2*(1 + Rqsd)*ss*                                          &
                      (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)))      &
)))) + 4*expuISr*expuISs*qdB*(2*dRqsu*qxy2**2*(1 + Rqru)*(1 + Rqsu)**2*         &
          tanhus**3*(64*(3*lam(1) + lam(3))**2 +                                &
            (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))**2) -         &
         quB*qxy2*tanhus**2*(dRqru*(1 + Rqsu)*                                  &
             (-1 + 2*qzt2*(1 + Rqru)**2*sr)*                                    &
             (2736*lam(1)**2 + 1944*lam(2)**2 + 960*lam(1)*lam(3) +             &
               88*lam(3)**2 +                                                   &
               9*(lam(4)**2 + 4*lam(5)**2 + lam(5)*lam(6) +                     &
                  lam(6)**2 + (8*lam(5) + lam(6))*lam(7) +                      &
                  4*lam(7)**2 + lam(4)*(lam(5) + 2*lam(6) + lam(7))) -          &
               54*lam(2)*(3*(lam(4) + 3*lam(5) + lam(6) + 3*lam(7)) -           &
                  16*lam(8)) -                                                  &
               48*(lam(4) + 2*lam(5) + lam(6) + 2*lam(7))*lam(8) +              &
               88*lam(8)**2) -                                                  &
            dRqsu*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                    &
             (64*(3*lam(1) + lam(3))**2 +                                       &
               (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))**2))       &
- 2*dRqsu*qxy2*(1 + Rqsu)*tanhus*                                               &
          (16*m2*(3*lam(1) + lam(3))*                                           &
             (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)) +              &
            qxy2*(1 + Rqru)*(1 + Rqsu)*                                         &
             (64*(3*lam(1) + lam(3))**2 +                                       &
               (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))**2)        &
+ qzt2*(1 + Rqru)*(1 + Rqsu)*(64*(3*lam(1) + lam(3))**2 +                       &
               (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10))**2)) +       &
         quB*(dRqru*(qxy2*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*            &
                (2736*lam(1)**2 + 1944*lam(2)**2 + 960*lam(1)*lam(3) +          &
                  88*lam(3)**2 +                                                &
                  9*(lam(4)**2 + 4*lam(5)**2 + lam(5)*lam(6) +                  &
                     lam(6)**2 + (8*lam(5) + lam(6))*lam(7) +                   &
                     4*lam(7)**2 + lam(4)*(lam(5) + 2*lam(6) + lam(7)))         &
- 54*lam(2)*(3*(lam(4) + 3*lam(5) + lam(6) + 3*lam(7)) - 16*lam(8)) -           &
                  48*(lam(4) + 2*lam(5) + lam(6) + 2*lam(7))*lam(8) +           &
                  88*lam(8)**2) +                                               &
               qzt2*(2736*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*            &
                   lam(1)**2 +                                                  &
                  1944*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*               &
                   lam(2)**2 - 88*lam(3)**2 - 88*Rqsu*lam(3)**2 +               &
                  176*qzt2*sr*lam(3)**2 + 352*qzt2*Rqru*sr*lam(3)**2 +          &
                  176*qzt2*Rqru**2*sr*lam(3)**2 +                               &
                  176*qzt2*Rqsu*sr*lam(3)**2 +                                  &
                  352*qzt2*Rqru*Rqsu*sr*lam(3)**2 +                             &
                  176*qzt2*Rqru**2*Rqsu*sr*lam(3)**2 -                          &
                  96*m2*sr*lam(3)*lam(4) - 96*m2*Rqru*sr*lam(3)*lam(4) -        &
                  9*lam(4)**2 - 9*Rqsu*lam(4)**2 +                              &
                  18*qzt2*sr*lam(4)**2 + 36*qzt2*Rqru*sr*lam(4)**2 +            &
                  18*qzt2*Rqru**2*sr*lam(4)**2 +                                &
                  18*qzt2*Rqsu*sr*lam(4)**2 +                                   &
                  36*qzt2*Rqru*Rqsu*sr*lam(4)**2 +                              &
                  18*qzt2*Rqru**2*Rqsu*sr*lam(4)**2 -                           &
                  480*m2*sr*lam(3)*lam(5) -                                     &
                  480*m2*Rqru*sr*lam(3)*lam(5) - 9*lam(4)*lam(5) -              &
                  9*Rqsu*lam(4)*lam(5) + 18*qzt2*sr*lam(4)*lam(5) +             &
                  36*qzt2*Rqru*sr*lam(4)*lam(5) +                               &
                  18*qzt2*Rqru**2*sr*lam(4)*lam(5) +                            &
                  18*qzt2*Rqsu*sr*lam(4)*lam(5) +                               &
                  36*qzt2*Rqru*Rqsu*sr*lam(4)*lam(5) +                          &
                  18*qzt2*Rqru**2*Rqsu*sr*lam(4)*lam(5) - 36*lam(5)**2 -        &
                  36*Rqsu*lam(5)**2 + 72*qzt2*sr*lam(5)**2 +                    &
                  144*qzt2*Rqru*sr*lam(5)**2 +                                  &
                  72*qzt2*Rqru**2*sr*lam(5)**2 +                                &
                  72*qzt2*Rqsu*sr*lam(5)**2 +                                   &
                  144*qzt2*Rqru*Rqsu*sr*lam(5)**2 +                             &
                  72*qzt2*Rqru**2*Rqsu*sr*lam(5)**2 -                           &
                  96*m2*sr*lam(3)*lam(6) - 96*m2*Rqru*sr*lam(3)*lam(6) -        &
                  18*lam(4)*lam(6) - 18*Rqsu*lam(4)*lam(6) +                    &
                  36*qzt2*sr*lam(4)*lam(6) +                                    &
                  72*qzt2*Rqru*sr*lam(4)*lam(6) +                               &
                  36*qzt2*Rqru**2*sr*lam(4)*lam(6) +                            &
                  36*qzt2*Rqsu*sr*lam(4)*lam(6) +                               &
                  72*qzt2*Rqru*Rqsu*sr*lam(4)*lam(6) +                          &
                  36*qzt2*Rqru**2*Rqsu*sr*lam(4)*lam(6) -                       &
                  9*lam(5)*lam(6) - 9*Rqsu*lam(5)*lam(6) +                      &
                  18*qzt2*sr*lam(5)*lam(6) +                                    &
                  36*qzt2*Rqru*sr*lam(5)*lam(6) +                               &
                  18*qzt2*Rqru**2*sr*lam(5)*lam(6) +                            &
                  18*qzt2*Rqsu*sr*lam(5)*lam(6) +                               &
                  36*qzt2*Rqru*Rqsu*sr*lam(5)*lam(6) +                          &
                  18*qzt2*Rqru**2*Rqsu*sr*lam(5)*lam(6) - 9*lam(6)**2 -         &
                  9*Rqsu*lam(6)**2 + 18*qzt2*sr*lam(6)**2 +                     &
                  36*qzt2*Rqru*sr*lam(6)**2 +                                   &
                  18*qzt2*Rqru**2*sr*lam(6)**2 +                                &
                  18*qzt2*Rqsu*sr*lam(6)**2 +                                   &
                  36*qzt2*Rqru*Rqsu*sr*lam(6)**2 +                              &
                  18*qzt2*Rqru**2*Rqsu*sr*lam(6)**2 -                           &
                  480*m2*sr*lam(3)*lam(7) -                                     &
                  480*m2*Rqru*sr*lam(3)*lam(7) - 9*lam(4)*lam(7) -              &
                  9*Rqsu*lam(4)*lam(7) + 18*qzt2*sr*lam(4)*lam(7) +             &
                  36*qzt2*Rqru*sr*lam(4)*lam(7) +                               &
                  18*qzt2*Rqru**2*sr*lam(4)*lam(7) +                            &
                  18*qzt2*Rqsu*sr*lam(4)*lam(7) +                               &
                  36*qzt2*Rqru*Rqsu*sr*lam(4)*lam(7) +                          &
                  18*qzt2*Rqru**2*Rqsu*sr*lam(4)*lam(7) -                       &
                  72*lam(5)*lam(7) - 72*Rqsu*lam(5)*lam(7) +                    &
                  144*qzt2*sr*lam(5)*lam(7) +                                   &
                  288*qzt2*Rqru*sr*lam(5)*lam(7) +                              &
                  144*qzt2*Rqru**2*sr*lam(5)*lam(7) +                           &
                  144*qzt2*Rqsu*sr*lam(5)*lam(7) +                              &
                  288*qzt2*Rqru*Rqsu*sr*lam(5)*lam(7) +                         &
                  144*qzt2*Rqru**2*Rqsu*sr*lam(5)*lam(7) -                      &
                  9*lam(6)*lam(7) - 9*Rqsu*lam(6)*lam(7) +                      &
                  18*qzt2*sr*lam(6)*lam(7) +                                    &
                  36*qzt2*Rqru*sr*lam(6)*lam(7) +                               &
                  18*qzt2*Rqru**2*sr*lam(6)*lam(7) +                            &
                  18*qzt2*Rqsu*sr*lam(6)*lam(7) +                               &
                  36*qzt2*Rqru*Rqsu*sr*lam(6)*lam(7) +                          &
                  18*qzt2*Rqru**2*Rqsu*sr*lam(6)*lam(7) - 36*lam(7)**2 -        &
                  36*Rqsu*lam(7)**2 + 72*qzt2*sr*lam(7)**2 +                    &
                  144*qzt2*Rqru*sr*lam(7)**2 +                                  &
                  72*qzt2*Rqru**2*sr*lam(7)**2 +                                &
                  72*qzt2*Rqsu*sr*lam(7)**2 +                                   &
                  144*qzt2*Rqru*Rqsu*sr*lam(7)**2 +                             &
                  72*qzt2*Rqru**2*Rqsu*sr*lam(7)**2 +                           &
                  54*lam(2)*(64*m2*(1 + Rqru)*sr*lam(3) -                       &
                     (1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                 &
                      (3*(lam(4) + 3*lam(5) + lam(6) + 3*lam(7)) -              &
                        16*lam(8))) +                                           &
                  16*(52*m2*(1 + Rqru)*sr*lam(3) -                              &
                     3*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*               &
                      (lam(4) + 2*lam(5) + lam(6) + 2*lam(7)))*lam(8) +         &
                  88*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                 &
                   lam(8)**2 +                                                  &
                  24*lam(1)*(40*(1 + Rqsu)*                                     &
                      (-1 + 2*qzt2*(1 + Rqru)**2*sr)*lam(3) +                   &
                     m2*(1 + Rqru)*sr*                                          &
                      (756*lam(2) - 39*lam(4) - 87*lam(5) - 39*lam(6) -         &
                        87*lam(7) + 176*lam(8))))) -                            &
            dRqsu*(qxy2*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*              &
                (64*(3*lam(1) + lam(3))**2 +                                    &
                  (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))**       &
                   2) + qzt2*(576*(1 + Rqru)*                                   &
                   (-1 + 2*qzt2*(1 + Rqsu)**2*ss)*lam(1)**2 +                   &
                  64*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                 &
                   lam(3)**2 +                                                  &
                  32*m2*(1 + Rqsu)*ss*lam(3)*                                   &
                   (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)) +        &
                  (1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                    &
                   (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10))**2       &
+ 96*lam(1)*(4*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*lam(3) +               &
                     m2*(1 + Rqsu)*ss*                                          &
                      (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)))      &
)))) - expuISr*tanhur*(-2*expdISs*qxy2*(1 + Rqsd)*tanhds**2*                    &
          (-(dRqru*qdB*qxy2*(1 + Rqru)**2*                                      &
               (16*(3*lam(1) - 2*lam(3))*(15*lam(1) + 2*lam(3)) +               &
                 45*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                    &
                 96*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -                &
                 64*lam(8)**2)) +                                               &
            dRqsd*quB*(16*m2*(9*lam(1)*                                         &
                   (5*lam(4) + 7*lam(5) - 5*lam(6) - 7*lam(7) -                 &
                     16*lam(9)) +                                               &
                  12*lam(3)*(-lam(4) - 5*lam(5) + lam(6) + 5*lam(7) +           &
                     8*lam(9)) + 8*(12*lam(1) - 5*lam(3))*lam(10)) +            &
               qzt2*(1 + Rqru)*(1 + Rqsd)*                                      &
                (64*(45*lam(1)**2 - 24*lam(1)*lam(3) + 14*lam(3)**2) +          &
                  9*(5*lam(4)**2 + 29*lam(5)**2 + 5*lam(6)**2 +                 &
                     2*lam(4)*(7*lam(5) - 5*lam(6) - 7*lam(7)) +                &
                     14*lam(6)*lam(7) + 29*lam(7)**2 -                          &
                     2*lam(5)*(7*lam(6) + 29*lam(7))) -                         &
                  288*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +          &
                  768*lam(9)**2 +                                               &
                  192*(lam(4) + 2*lam(5) - lam(6) - 2*lam(7) - 4*lam(9))*       &
                   lam(10) + 224*lam(10)**2))) +                                &
         expdISs*qdB*quB*qzt2*tanhds*                                           &
          (dRqru*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                     &
             (16*(3*lam(1) - 2*lam(3))*(15*lam(1) + 2*lam(3)) +                 &
               45*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                      &
               96*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -                  &
               64*lam(8)**2) -                                                  &
            dRqsd*(2880*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*              &
                lam(1)**2 + 896*(1 + Rqru)*                                     &
                (-1 + 2*qzt2*(1 + Rqsd)**2*ss)*lam(3)**2 -                      &
               128*m2*(1 + Rqsd)*ss*lam(3)*                                     &
                (3*lam(4) + 15*lam(5) -                                         &
                  3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10)) +              &
               (1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                       &
                (9*(5*lam(4)**2 + 29*lam(5)**2 + 5*lam(6)**2 +                  &
                     2*lam(4)*(7*lam(5) - 5*lam(6) - 7*lam(7)) +                &
                     14*lam(6)*lam(7) + 29*lam(7)**2 -                          &
                     2*lam(5)*(7*lam(6) + 29*lam(7))) -                         &
                  288*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +          &
                  768*lam(9)**2 +                                               &
                  192*(lam(4) + 2*lam(5) - lam(6) - 2*lam(7) -                  &
                     4*lam(9))*lam(10) + 224*lam(10)**2) +                      &
               96*lam(1)*(-16*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*        &
                   lam(3) + m2*(1 + Rqsd)*ss*                                   &
                   (15*lam(4) + 21*lam(5) -                                     &
                     3*(5*lam(6) + 7*lam(7) + 16*lam(9)) + 32*lam(10)))))       &
- 2*qdB*(dRqru*qxy2*(1 + Rqru)*                                                 &
             (expdISs*(-48*m2*(15*lam(1) - 4*lam(3))*                           &
                   (lam(4) - lam(5) + lam(6) - lam(7)) +                        &
                  256*m2*(3*lam(1) + lam(3))*lam(8) +                           &
                  qxy2*(1 + Rqru)*(1 + Rqsd)*                                   &
                   (16*(3*lam(1) - 2*lam(3))*(15*lam(1) + 2*lam(3)) +           &
                     45*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                &
                     96*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -            &
                     64*lam(8)**2) +                                            &
                  qzt2*(1 + Rqru)*(1 + Rqsd)*                                   &
                   (16*(3*lam(1) - 2*lam(3))*(15*lam(1) + 2*lam(3)) +           &
                     45*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                &
                     96*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -            &
                     64*lam(8)**2)) +                                           &
               4*expuISs*(qxy2*(1 + Rqru)*(1 + Rqsu)*                           &
                   (2736*lam(1)**2 + 1944*lam(2)**2 +                           &
                     960*lam(1)*lam(3) + 88*lam(3)**2 +                         &
                     9*(lam(4)**2 + 4*lam(5)**2 + lam(5)*lam(6) +               &
                        lam(6)**2 + (8*lam(5) + lam(6))*lam(7) +                &
                        4*lam(7)**2 +                                           &
                        lam(4)*(lam(5) + 2*lam(6) + lam(7))) -                  &
                     54*lam(2)*                                                 &
                      (3*(lam(4) + 3*lam(5) + lam(6) + 3*lam(7)) -              &
                        16*lam(8)) -                                            &
                     48*(lam(4) + 2*lam(5) + lam(6) + 2*lam(7))*lam(8) +        &
                     88*lam(8)**2) +                                            &
                  qzt2*(1 + Rqru)*(1 + Rqsu)*                                   &
                   (2736*lam(1)**2 + 1944*lam(2)**2 +                           &
                     960*lam(1)*lam(3) + 88*lam(3)**2 +                         &
                     9*(lam(4)**2 + 4*lam(5)**2 + lam(5)*lam(6) +               &
                        lam(6)**2 + (8*lam(5) + lam(6))*lam(7) +                &
                        4*lam(7)**2 +                                           &
                        lam(4)*(lam(5) + 2*lam(6) + lam(7))) -                  &
                     54*lam(2)*                                                 &
                      (3*(lam(4) + 3*lam(5) + lam(6) + 3*lam(7)) -              &
                        16*lam(8)) -                                            &
                     48*(lam(4) + 2*lam(5) + lam(6) + 2*lam(7))*lam(8) +        &
                     88*lam(8)**2) +                                            &
                  4*m2*(4*lam(3)*                                               &
                      (108*lam(2) -                                             &
                        3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +             &
                        26*lam(8)) +                                            &
                     3*lam(1)*                                                  &
                      (756*lam(2) - 39*lam(4) - 87*lam(5) - 39*lam(6) -         &
                        87*lam(7) + 176*lam(8))))) -                            &
            2*expuISs*quB*qzt2*tanhus*                                          &
             (dRqru*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                  &
                (2736*lam(1)**2 + 1944*lam(2)**2 + 960*lam(1)*lam(3) +          &
                  88*lam(3)**2 +                                                &
                  9*(lam(4)**2 + 4*lam(5)**2 + lam(5)*lam(6) +                  &
                     lam(6)**2 + (8*lam(5) + lam(6))*lam(7) +                   &
                     4*lam(7)**2 + lam(4)*(lam(5) + 2*lam(6) + lam(7)))         &
- 54*lam(2)*(3*(lam(4) + 3*lam(5) + lam(6) + 3*lam(7)) - 16*lam(8)) -           &
                  48*(lam(4) + 2*lam(5) + lam(6) + 2*lam(7))*lam(8) +           &
                  88*lam(8)**2) -                                               &
               dRqsu*(576*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*            &
                   lam(1)**2 +                                                  &
                  64*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                 &
                   lam(3)**2 +                                                  &
                  32*m2*(1 + Rqsu)*ss*lam(3)*                                   &
                   (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)) +        &
                  (1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                    &
                   (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10))**2       &
+ 96*lam(1)*(4*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*lam(3) +               &
                     m2*(1 + Rqsu)*ss*                                          &
                      (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)))      &
)) + 4*expuISs*qxy2*(1 + Rqsu)*tanhus**2*                                       &
             (-(dRqru*qxy2*(1 + Rqru)**2*                                       &
                  (2736*lam(1)**2 + 1944*lam(2)**2 + 960*lam(1)*lam(3) +        &
                    88*lam(3)**2 +                                              &
                    9*(lam(4)**2 + 4*lam(5)**2 + lam(5)*lam(6) +                &
                       lam(6)**2 + (8*lam(5) + lam(6))*lam(7) +                 &
                       4*lam(7)**2 + lam(4)*(lam(5) + 2*lam(6) + lam(7)))       &
- 54*lam(2)*(3*(lam(4) + 3*lam(5) + lam(6) + 3*lam(7)) - 16*lam(8)) -           &
                    48*(lam(4) + 2*lam(5) + lam(6) + 2*lam(7))*lam(8) +         &
                    88*lam(8)**2)) +                                            &
               dRqsu*(16*m2*(3*lam(1) + lam(3))*                                &
                   (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)) +        &
                  qzt2*(1 + Rqru)*(1 + Rqsu)*                                   &
                   (64*(3*lam(1) + lam(3))**2 +                                 &
                     (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10))**2     &
))))))/1296.
  dtLAM(2)=                                                       &
   (2*dRqrd*expdISr*quB*qxy2**2*(1 + Rqrd)**2*tanhdr**3*                        &
       (expuISs*(1 + Rqsu)*(-1 + tanhus**2)*                                    &
          (576*lam(2)**2 + 576*lam(2)*lam(8) +                                  &
            144*lam(8)*(-lam(1) + lam(8)) -                                     &
            12*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +                  &
               14*lam(8)) + (3*lam(4) + 3*lam(5) - 3*lam(6) -                   &
               3*lam(7) + 10*lam(9) - 12*lam(10))**2) +                         &
         2*expdISs*(1 + Rqsd)*(-1 + tanhds**2)*                                 &
          (576*lam(2)**2 - 18*lam(3)*lam(4) + 9*lam(4)**2 -                     &
            18*lam(3)*lam(5) - 36*lam(4)*lam(5) + 63*lam(5)**2 -                &
            18*lam(3)*lam(6) - 18*lam(4)*lam(6) + 36*lam(5)*lam(6) +            &
            9*lam(6)**2 - 18*lam(3)*lam(7) + 36*lam(4)*lam(7) -                 &
            126*lam(5)*lam(7) - 36*lam(6)*lam(7) + 63*lam(7)**2 +               &
            288*lam(2)*(2*lam(3) - lam(8)) + 96*lam(3)*lam(8) +                 &
            144*lam(8)**2 + 108*lam(1)*                                         &
             (28*lam(2) - lam(4) - 3*lam(5) - lam(6) - 3*lam(7) +               &
               6*lam(8)) - 30*lam(4)*lam(9) + 150*lam(5)*lam(9) +               &
            30*lam(6)*lam(9) - 150*lam(7)*lam(9) + 100*lam(9)**2 +              &
            18*(-5*lam(4) + 11*lam(5) + 5*lam(6) - 11*lam(7) +                  &
               10*lam(9))*lam(10) + 228*lam(10)**2)) +                          &
      2*dRqru*expuISr*qdB*qxy2**2*(1 + Rqru)**2*tanhur**3*                      &
       (expdISs*(1 + Rqsd)*(-1 + tanhds**2)*                                    &
          (576*lam(2)**2 + 576*lam(2)*lam(8) +                                  &
            144*lam(8)*(-lam(1) + lam(8)) -                                     &
            12*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +                  &
               14*lam(8)) + (3*lam(4) + 3*lam(5) - 3*lam(6) -                   &
               3*lam(7) + 10*lam(9) - 12*lam(10))**2) +                         &
         2*expuISs*(1 + Rqsu)*(-1 + tanhus**2)*                                 &
          (576*lam(2)**2 - 18*lam(3)*lam(4) + 9*lam(4)**2 -                     &
            18*lam(3)*lam(5) - 36*lam(4)*lam(5) + 63*lam(5)**2 -                &
            18*lam(3)*lam(6) - 18*lam(4)*lam(6) + 36*lam(5)*lam(6) +            &
            9*lam(6)**2 - 18*lam(3)*lam(7) + 36*lam(4)*lam(7) -                 &
            126*lam(5)*lam(7) - 36*lam(6)*lam(7) + 63*lam(7)**2 +               &
            288*lam(2)*(2*lam(3) - lam(8)) + 96*lam(3)*lam(8) +                 &
            144*lam(8)**2 + 108*lam(1)*                                         &
             (28*lam(2) - lam(4) - 3*lam(5) - lam(6) - 3*lam(7) +               &
               6*lam(8)) - 30*lam(4)*lam(9) + 150*lam(5)*lam(9) +               &
            30*lam(6)*lam(9) - 150*lam(7)*lam(9) + 100*lam(9)**2 +              &
            18*(-5*lam(4) + 11*lam(5) + 5*lam(6) - 11*lam(7) +                  &
               10*lam(9))*lam(10) + 228*lam(10)**2)) +                          &
      2*expdISr*expdISs*quB*(18*dRqsd*qxy2*(1 + Rqsd)*tanhds*                   &
          (-16*m2*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) -                  &
            qzt2*(1 + Rqrd)*(1 + Rqsd)*                                         &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            qxy2*(1 + Rqrd)*(1 + Rqsd)*(-1 + tanhds**2)*                        &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2)) -         &
         9*dRqsd*qdB*(-(qxy2*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*         &
               (-1 + tanhds**2)*                                                &
               (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2))         &
+ qzt2*(16*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*lam(2)**2 +                &
               32*m2*(1 + Rqsd)*ss*lam(2)*                                      &
                (lam(4) + lam(5) + lam(6) + lam(7)) +                           &
               (1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                       &
                (lam(4) + lam(5) + lam(6) + lam(7))**2)) +                      &
         dRqrd*qdB*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                   &
          (qxy2 + qzt2 - qxy2*tanhds**2)*                                       &
          (576*lam(2)**2 - 18*lam(3)*lam(4) + 9*lam(4)**2 -                     &
            18*lam(3)*lam(5) - 36*lam(4)*lam(5) + 63*lam(5)**2 -                &
            18*lam(3)*lam(6) - 18*lam(4)*lam(6) + 36*lam(5)*lam(6) +            &
            9*lam(6)**2 - 18*lam(3)*lam(7) + 36*lam(4)*lam(7) -                 &
            126*lam(5)*lam(7) - 36*lam(6)*lam(7) + 63*lam(7)**2 +               &
            288*lam(2)*(2*lam(3) - lam(8)) + 96*lam(3)*lam(8) +                 &
            144*lam(8)**2 + 108*lam(1)*                                         &
             (28*lam(2) - lam(4) - 3*lam(5) - lam(6) - 3*lam(7) +               &
               6*lam(8)) - 30*lam(4)*lam(9) + 150*lam(5)*lam(9) +               &
            30*lam(6)*lam(9) - 150*lam(7)*lam(9) + 100*lam(9)**2 +              &
            18*(-5*lam(4) + 11*lam(5) + 5*lam(6) - 11*lam(7) +                  &
               10*lam(9))*lam(10) + 228*lam(10)**2) +                           &
         48*dRqrd*m2*qdB*qzt2*(1 + Rqrd)*sr*                                    &
          (108*lam(2)**2 + (6*lam(1) + lam(3))*(24*lam(1) + 5*lam(3)) +         &
            lam(8)*(-15*lam(5) - 3*lam(6) + 7*lam(8) - 10*lam(9) -              &
               9*lam(10)) + lam(2)*                                             &
             (-15*lam(4) - 15*lam(5) - 3*lam(6) - 39*lam(7) + 48*lam(8) +       &
               10*lam(9) + 30*lam(10)))) +                                      &
      2*expuISr*expuISs*qdB*(18*dRqsu*qxy2*(1 + Rqsu)*tanhus*                   &
          (-16*m2*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) -                  &
            qzt2*(1 + Rqru)*(1 + Rqsu)*                                         &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            qxy2*(1 + Rqru)*(1 + Rqsu)*(-1 + tanhus**2)*                        &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2)) -         &
         9*dRqsu*quB*(-(qxy2*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*         &
               (-1 + tanhus**2)*                                                &
               (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2))         &
+ qzt2*(16*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*lam(2)**2 +                &
               32*m2*(1 + Rqsu)*ss*lam(2)*                                      &
                (lam(4) + lam(5) + lam(6) + lam(7)) +                           &
               (1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                       &
                (lam(4) + lam(5) + lam(6) + lam(7))**2)) +                      &
         dRqru*quB*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                   &
          (qxy2 + qzt2 - qxy2*tanhus**2)*                                       &
          (576*lam(2)**2 - 18*lam(3)*lam(4) + 9*lam(4)**2 -                     &
            18*lam(3)*lam(5) - 36*lam(4)*lam(5) + 63*lam(5)**2 -                &
            18*lam(3)*lam(6) - 18*lam(4)*lam(6) + 36*lam(5)*lam(6) +            &
            9*lam(6)**2 - 18*lam(3)*lam(7) + 36*lam(4)*lam(7) -                 &
            126*lam(5)*lam(7) - 36*lam(6)*lam(7) + 63*lam(7)**2 +               &
            288*lam(2)*(2*lam(3) - lam(8)) + 96*lam(3)*lam(8) +                 &
            144*lam(8)**2 + 108*lam(1)*                                         &
             (28*lam(2) - lam(4) - 3*lam(5) - lam(6) - 3*lam(7) +               &
               6*lam(8)) - 30*lam(4)*lam(9) + 150*lam(5)*lam(9) +               &
            30*lam(6)*lam(9) - 150*lam(7)*lam(9) + 100*lam(9)**2 +              &
            18*(-5*lam(4) + 11*lam(5) + 5*lam(6) - 11*lam(7) +                  &
               10*lam(9))*lam(10) + 228*lam(10)**2) +                           &
         48*dRqru*m2*quB*qzt2*(1 + Rqru)*sr*                                    &
          (108*lam(2)**2 + (6*lam(1) + lam(3))*(24*lam(1) + 5*lam(3)) +         &
            lam(8)*(-15*lam(5) - 3*lam(6) + 7*lam(8) - 10*lam(9) -              &
               9*lam(10)) + lam(2)*                                             &
             (-15*lam(4) - 15*lam(5) - 3*lam(6) - 39*lam(7) + 48*lam(8) +       &
               10*lam(9) + 30*lam(10)))) +                                      &
      expdISr*tanhdr*(-2*expuISs*qxy2*(1 + Rqrd)*(1 + Rqsu)*tanhus**2*          &
          (9*dRqsu*qdB*qzt2*(1 + Rqsu)*                                         &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            dRqrd*quB*qxy2*(1 + Rqrd)*                                          &
             (576*lam(2)**2 + 576*lam(2)*lam(8) +                               &
               144*lam(8)*(-lam(1) + lam(8)) -                                  &
               12*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +               &
                  14*lam(8)) +                                                  &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -         &
                  12*lam(10))**2)) -                                            &
         4*expdISs*quB*qxy2*(1 + Rqrd)*(1 + Rqsd)*tanhds**2*                    &
          (9*dRqsd*qzt2*(1 + Rqsd)*                                             &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            dRqrd*qxy2*(1 + Rqrd)*                                              &
             (576*lam(2)**2 - 18*lam(3)*lam(4) + 9*lam(4)**2 -                  &
               18*lam(3)*lam(5) - 36*lam(4)*lam(5) + 63*lam(5)**2 -             &
               18*lam(3)*lam(6) - 18*lam(4)*lam(6) + 36*lam(5)*lam(6) +         &
               9*lam(6)**2 - 18*lam(3)*lam(7) + 36*lam(4)*lam(7) -              &
               126*lam(5)*lam(7) - 36*lam(6)*lam(7) + 63*lam(7)**2 +            &
               288*lam(2)*(2*lam(3) - lam(8)) + 96*lam(3)*lam(8) +              &
               144*lam(8)**2 +                                                  &
               108*lam(1)*(28*lam(2) - lam(4) - 3*lam(5) - lam(6) -             &
                  3*lam(7) + 6*lam(8)) - 30*lam(4)*lam(9) +                     &
               150*lam(5)*lam(9) + 30*lam(6)*lam(9) -                           &
               150*lam(7)*lam(9) + 100*lam(9)**2 +                              &
               18*(-5*lam(4) + 11*lam(5) + 5*lam(6) - 11*lam(7) +               &
                  10*lam(9))*lam(10) + 228*lam(10)**2)) -                       &
         2*expdISs*qdB*quB*qzt2*tanhds*                                         &
          (9*dRqsd*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                   &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            dRqrd*(-576*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*              &
                lam(2)**2 + 576*(1 + Rqsd)*                                     &
                (-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(2)*lam(3) +                  &
               288*m2*(1 + Rqrd)*sr*lam(2)*lam(4) + 18*lam(3)*lam(4) +          &
               18*Rqsd*lam(3)*lam(4) - 36*qzt2*sr*lam(3)*lam(4) -               &
               72*qzt2*Rqrd*sr*lam(3)*lam(4) -                                  &
               36*qzt2*Rqrd**2*sr*lam(3)*lam(4) -                               &
               36*qzt2*Rqsd*sr*lam(3)*lam(4) -                                  &
               72*qzt2*Rqrd*Rqsd*sr*lam(3)*lam(4) -                             &
               36*qzt2*Rqrd**2*Rqsd*sr*lam(3)*lam(4) + 9*lam(4)**2 +            &
               9*Rqsd*lam(4)**2 - 18*qzt2*sr*lam(4)**2 -                        &
               36*qzt2*Rqrd*sr*lam(4)**2 -                                      &
               18*qzt2*Rqrd**2*sr*lam(4)**2 -                                   &
               18*qzt2*Rqsd*sr*lam(4)**2 -                                      &
               36*qzt2*Rqrd*Rqsd*sr*lam(4)**2 -                                 &
               18*qzt2*Rqrd**2*Rqsd*sr*lam(4)**2 -                              &
               576*m2*(1 + Rqrd)*sr*lam(2)*lam(5) + 18*lam(3)*lam(5) +          &
               18*Rqsd*lam(3)*lam(5) - 36*qzt2*sr*lam(3)*lam(5) -               &
               72*qzt2*Rqrd*sr*lam(3)*lam(5) -                                  &
               36*qzt2*Rqrd**2*sr*lam(3)*lam(5) -                               &
               36*qzt2*Rqsd*sr*lam(3)*lam(5) -                                  &
               72*qzt2*Rqrd*Rqsd*sr*lam(3)*lam(5) -                             &
               36*qzt2*Rqrd**2*Rqsd*sr*lam(3)*lam(5) -                          &
               36*lam(4)*lam(5) - 36*Rqsd*lam(4)*lam(5) +                       &
               72*qzt2*sr*lam(4)*lam(5) +                                       &
               144*qzt2*Rqrd*sr*lam(4)*lam(5) +                                 &
               72*qzt2*Rqrd**2*sr*lam(4)*lam(5) +                               &
               72*qzt2*Rqsd*sr*lam(4)*lam(5) +                                  &
               144*qzt2*Rqrd*Rqsd*sr*lam(4)*lam(5) +                            &
               72*qzt2*Rqrd**2*Rqsd*sr*lam(4)*lam(5) + 63*lam(5)**2 +           &
               63*Rqsd*lam(5)**2 - 126*qzt2*sr*lam(5)**2 -                      &
               252*qzt2*Rqrd*sr*lam(5)**2 -                                     &
               126*qzt2*Rqrd**2*sr*lam(5)**2 -                                  &
               126*qzt2*Rqsd*sr*lam(5)**2 -                                     &
               252*qzt2*Rqrd*Rqsd*sr*lam(5)**2 -                                &
               126*qzt2*Rqrd**2*Rqsd*sr*lam(5)**2 -                             &
               288*m2*(1 + Rqrd)*sr*lam(2)*lam(6) + 18*lam(3)*lam(6) +          &
               18*Rqsd*lam(3)*lam(6) - 36*qzt2*sr*lam(3)*lam(6) -               &
               72*qzt2*Rqrd*sr*lam(3)*lam(6) -                                  &
               36*qzt2*Rqrd**2*sr*lam(3)*lam(6) -                               &
               36*qzt2*Rqsd*sr*lam(3)*lam(6) -                                  &
               72*qzt2*Rqrd*Rqsd*sr*lam(3)*lam(6) -                             &
               36*qzt2*Rqrd**2*Rqsd*sr*lam(3)*lam(6) -                          &
               18*lam(4)*lam(6) - 18*Rqsd*lam(4)*lam(6) +                       &
               36*qzt2*sr*lam(4)*lam(6) +                                       &
               72*qzt2*Rqrd*sr*lam(4)*lam(6) +                                  &
               36*qzt2*Rqrd**2*sr*lam(4)*lam(6) +                               &
               36*qzt2*Rqsd*sr*lam(4)*lam(6) +                                  &
               72*qzt2*Rqrd*Rqsd*sr*lam(4)*lam(6) +                             &
               36*qzt2*Rqrd**2*Rqsd*sr*lam(4)*lam(6) +                          &
               36*lam(5)*lam(6) + 36*Rqsd*lam(5)*lam(6) -                       &
               72*qzt2*sr*lam(5)*lam(6) -                                       &
               144*qzt2*Rqrd*sr*lam(5)*lam(6) -                                 &
               72*qzt2*Rqrd**2*sr*lam(5)*lam(6) -                               &
               72*qzt2*Rqsd*sr*lam(5)*lam(6) -                                  &
               144*qzt2*Rqrd*Rqsd*sr*lam(5)*lam(6) -                            &
               72*qzt2*Rqrd**2*Rqsd*sr*lam(5)*lam(6) + 9*lam(6)**2 +            &
               9*Rqsd*lam(6)**2 - 18*qzt2*sr*lam(6)**2 -                        &
               36*qzt2*Rqrd*sr*lam(6)**2 -                                      &
               18*qzt2*Rqrd**2*sr*lam(6)**2 -                                   &
               18*qzt2*Rqsd*sr*lam(6)**2 -                                      &
               36*qzt2*Rqrd*Rqsd*sr*lam(6)**2 -                                 &
               18*qzt2*Rqrd**2*Rqsd*sr*lam(6)**2 +                              &
               576*m2*(1 + Rqrd)*sr*lam(2)*lam(7) + 18*lam(3)*lam(7) +          &
               18*Rqsd*lam(3)*lam(7) - 36*qzt2*sr*lam(3)*lam(7) -               &
               72*qzt2*Rqrd*sr*lam(3)*lam(7) -                                  &
               36*qzt2*Rqrd**2*sr*lam(3)*lam(7) -                               &
               36*qzt2*Rqsd*sr*lam(3)*lam(7) -                                  &
               72*qzt2*Rqrd*Rqsd*sr*lam(3)*lam(7) -                             &
               36*qzt2*Rqrd**2*Rqsd*sr*lam(3)*lam(7) +                          &
               36*lam(4)*lam(7) + 36*Rqsd*lam(4)*lam(7) -                       &
               72*qzt2*sr*lam(4)*lam(7) -                                       &
               144*qzt2*Rqrd*sr*lam(4)*lam(7) -                                 &
               72*qzt2*Rqrd**2*sr*lam(4)*lam(7) -                               &
               72*qzt2*Rqsd*sr*lam(4)*lam(7) -                                  &
               144*qzt2*Rqrd*Rqsd*sr*lam(4)*lam(7) -                            &
               72*qzt2*Rqrd**2*Rqsd*sr*lam(4)*lam(7) -                          &
               126*lam(5)*lam(7) - 126*Rqsd*lam(5)*lam(7) +                     &
               252*qzt2*sr*lam(5)*lam(7) +                                      &
               504*qzt2*Rqrd*sr*lam(5)*lam(7) +                                 &
               252*qzt2*Rqrd**2*sr*lam(5)*lam(7) +                              &
               252*qzt2*Rqsd*sr*lam(5)*lam(7) +                                 &
               504*qzt2*Rqrd*Rqsd*sr*lam(5)*lam(7) +                            &
               252*qzt2*Rqrd**2*Rqsd*sr*lam(5)*lam(7) -                         &
               36*lam(6)*lam(7) - 36*Rqsd*lam(6)*lam(7) +                       &
               72*qzt2*sr*lam(6)*lam(7) +                                       &
               144*qzt2*Rqrd*sr*lam(6)*lam(7) +                                 &
               72*qzt2*Rqrd**2*sr*lam(6)*lam(7) +                               &
               72*qzt2*Rqsd*sr*lam(6)*lam(7) +                                  &
               144*qzt2*Rqrd*Rqsd*sr*lam(6)*lam(7) +                            &
               72*qzt2*Rqrd**2*Rqsd*sr*lam(6)*lam(7) + 63*lam(7)**2 +           &
               63*Rqsd*lam(7)**2 - 126*qzt2*sr*lam(7)**2 -                      &
               252*qzt2*Rqrd*sr*lam(7)**2 -                                     &
               126*qzt2*Rqrd**2*sr*lam(7)**2 -                                  &
               126*qzt2*Rqsd*sr*lam(7)**2 -                                     &
               252*qzt2*Rqrd*Rqsd*sr*lam(7)**2 -                                &
               126*qzt2*Rqrd**2*Rqsd*sr*lam(7)**2 +                             &
               288*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(2)*            &
                lam(8) - 96*lam(3)*lam(8) - 96*Rqsd*lam(3)*lam(8) +             &
               192*qzt2*sr*lam(3)*lam(8) +                                      &
               384*qzt2*Rqrd*sr*lam(3)*lam(8) +                                 &
               192*qzt2*Rqrd**2*sr*lam(3)*lam(8) +                              &
               192*qzt2*Rqsd*sr*lam(3)*lam(8) +                                 &
               384*qzt2*Rqrd*Rqsd*sr*lam(3)*lam(8) +                            &
               192*qzt2*Rqrd**2*Rqsd*sr*lam(3)*lam(8) -                         &
               72*m2*sr*lam(4)*lam(8) - 72*m2*Rqrd*sr*lam(4)*lam(8) +           &
               360*m2*sr*lam(5)*lam(8) + 360*m2*Rqrd*sr*lam(5)*lam(8) +         &
               72*m2*sr*lam(6)*lam(8) + 72*m2*Rqrd*sr*lam(6)*lam(8) -           &
               360*m2*sr*lam(7)*lam(8) - 360*m2*Rqrd*sr*lam(7)*lam(8) +         &
               144*lam(8)**2 + 144*Rqsd*lam(8)**2 -                             &
               288*qzt2*sr*lam(8)**2 - 576*qzt2*Rqrd*sr*lam(8)**2 -             &
               288*qzt2*Rqrd**2*sr*lam(8)**2 -                                  &
               288*qzt2*Rqsd*sr*lam(8)**2 -                                     &
               576*qzt2*Rqrd*Rqsd*sr*lam(8)**2 -                                &
               288*qzt2*Rqrd**2*Rqsd*sr*lam(8)**2 +                             &
               108*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(1)*            &
                (28*lam(2) - lam(4) - 3*lam(5) - lam(6) - 3*lam(7) +            &
                  6*lam(8)) - 480*m2*(1 + Rqrd)*sr*lam(2)*lam(9) -              &
               30*lam(4)*lam(9) - 30*Rqsd*lam(4)*lam(9) +                       &
               60*qzt2*sr*lam(4)*lam(9) +                                       &
               120*qzt2*Rqrd*sr*lam(4)*lam(9) +                                 &
               60*qzt2*Rqrd**2*sr*lam(4)*lam(9) +                               &
               60*qzt2*Rqsd*sr*lam(4)*lam(9) +                                  &
               120*qzt2*Rqrd*Rqsd*sr*lam(4)*lam(9) +                            &
               60*qzt2*Rqrd**2*Rqsd*sr*lam(4)*lam(9) +                          &
               150*lam(5)*lam(9) + 150*Rqsd*lam(5)*lam(9) -                     &
               300*qzt2*sr*lam(5)*lam(9) -                                      &
               600*qzt2*Rqrd*sr*lam(5)*lam(9) -                                 &
               300*qzt2*Rqrd**2*sr*lam(5)*lam(9) -                              &
               300*qzt2*Rqsd*sr*lam(5)*lam(9) -                                 &
               600*qzt2*Rqrd*Rqsd*sr*lam(5)*lam(9) -                            &
               300*qzt2*Rqrd**2*Rqsd*sr*lam(5)*lam(9) +                         &
               30*lam(6)*lam(9) + 30*Rqsd*lam(6)*lam(9) -                       &
               60*qzt2*sr*lam(6)*lam(9) -                                       &
               120*qzt2*Rqrd*sr*lam(6)*lam(9) -                                 &
               60*qzt2*Rqrd**2*sr*lam(6)*lam(9) -                               &
               60*qzt2*Rqsd*sr*lam(6)*lam(9) -                                  &
               120*qzt2*Rqrd*Rqsd*sr*lam(6)*lam(9) -                            &
               60*qzt2*Rqrd**2*Rqsd*sr*lam(6)*lam(9) -                          &
               150*lam(7)*lam(9) - 150*Rqsd*lam(7)*lam(9) +                     &
               300*qzt2*sr*lam(7)*lam(9) +                                      &
               600*qzt2*Rqrd*sr*lam(7)*lam(9) +                                 &
               300*qzt2*Rqrd**2*sr*lam(7)*lam(9) +                              &
               300*qzt2*Rqsd*sr*lam(7)*lam(9) +                                 &
               600*qzt2*Rqrd*Rqsd*sr*lam(7)*lam(9) +                            &
               300*qzt2*Rqrd**2*Rqsd*sr*lam(7)*lam(9) +                         &
               480*m2*sr*lam(8)*lam(9) + 480*m2*Rqrd*sr*lam(8)*lam(9) +         &
               100*lam(9)**2 + 100*Rqsd*lam(9)**2 -                             &
               200*qzt2*sr*lam(9)**2 - 400*qzt2*Rqrd*sr*lam(9)**2 -             &
               200*qzt2*Rqrd**2*sr*lam(9)**2 -                                  &
               200*qzt2*Rqsd*sr*lam(9)**2 -                                     &
               400*qzt2*Rqrd*Rqsd*sr*lam(9)**2 -                                &
               200*qzt2*Rqrd**2*Rqsd*sr*lam(9)**2 -                             &
               18*(8*m2*(1 + Rqrd)*sr*(10*lam(2) - 3*lam(8)) -                  &
                  (1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                    &
                   (5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -               &
                     10*lam(9)))*lam(10) -                                      &
               228*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(10)**2))       &
- expuISs*qdB*quB*qzt2*tanhus*                                                  &
          (9*dRqsu*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                   &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) -          &
            dRqrd*(576*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*               &
                lam(2)**2 - 288*m2*(1 + Rqrd)*sr*lam(2)*lam(4) -                &
               9*lam(4)**2 - 9*Rqsu*lam(4)**2 + 18*qzt2*sr*lam(4)**2 +          &
               36*qzt2*Rqrd*sr*lam(4)**2 +                                      &
               18*qzt2*Rqrd**2*sr*lam(4)**2 +                                   &
               18*qzt2*Rqsu*sr*lam(4)**2 +                                      &
               36*qzt2*Rqrd*Rqsu*sr*lam(4)**2 +                                 &
               18*qzt2*Rqrd**2*Rqsu*sr*lam(4)**2 -                              &
               288*m2*(1 + Rqrd)*sr*lam(2)*lam(5) - 18*lam(4)*lam(5) -          &
               18*Rqsu*lam(4)*lam(5) + 36*qzt2*sr*lam(4)*lam(5) +               &
               72*qzt2*Rqrd*sr*lam(4)*lam(5) +                                  &
               36*qzt2*Rqrd**2*sr*lam(4)*lam(5) +                               &
               36*qzt2*Rqsu*sr*lam(4)*lam(5) +                                  &
               72*qzt2*Rqrd*Rqsu*sr*lam(4)*lam(5) +                             &
               36*qzt2*Rqrd**2*Rqsu*sr*lam(4)*lam(5) - 9*lam(5)**2 -            &
               9*Rqsu*lam(5)**2 + 18*qzt2*sr*lam(5)**2 +                        &
               36*qzt2*Rqrd*sr*lam(5)**2 +                                      &
               18*qzt2*Rqrd**2*sr*lam(5)**2 +                                   &
               18*qzt2*Rqsu*sr*lam(5)**2 +                                      &
               36*qzt2*Rqrd*Rqsu*sr*lam(5)**2 +                                 &
               18*qzt2*Rqrd**2*Rqsu*sr*lam(5)**2 +                              &
               288*m2*(1 + Rqrd)*sr*lam(2)*lam(6) + 18*lam(4)*lam(6) +          &
               18*Rqsu*lam(4)*lam(6) - 36*qzt2*sr*lam(4)*lam(6) -               &
               72*qzt2*Rqrd*sr*lam(4)*lam(6) -                                  &
               36*qzt2*Rqrd**2*sr*lam(4)*lam(6) -                               &
               36*qzt2*Rqsu*sr*lam(4)*lam(6) -                                  &
               72*qzt2*Rqrd*Rqsu*sr*lam(4)*lam(6) -                             &
               36*qzt2*Rqrd**2*Rqsu*sr*lam(4)*lam(6) +                          &
               18*lam(5)*lam(6) + 18*Rqsu*lam(5)*lam(6) -                       &
               36*qzt2*sr*lam(5)*lam(6) -                                       &
               72*qzt2*Rqrd*sr*lam(5)*lam(6) -                                  &
               36*qzt2*Rqrd**2*sr*lam(5)*lam(6) -                               &
               36*qzt2*Rqsu*sr*lam(5)*lam(6) -                                  &
               72*qzt2*Rqrd*Rqsu*sr*lam(5)*lam(6) -                             &
               36*qzt2*Rqrd**2*Rqsu*sr*lam(5)*lam(6) - 9*lam(6)**2 -            &
               9*Rqsu*lam(6)**2 + 18*qzt2*sr*lam(6)**2 +                        &
               36*qzt2*Rqrd*sr*lam(6)**2 +                                      &
               18*qzt2*Rqrd**2*sr*lam(6)**2 +                                   &
               18*qzt2*Rqsu*sr*lam(6)**2 +                                      &
               36*qzt2*Rqrd*Rqsu*sr*lam(6)**2 +                                 &
               18*qzt2*Rqrd**2*Rqsu*sr*lam(6)**2 +                              &
               288*m2*(1 + Rqrd)*sr*lam(2)*lam(7) + 18*lam(4)*lam(7) +          &
               18*Rqsu*lam(4)*lam(7) - 36*qzt2*sr*lam(4)*lam(7) -               &
               72*qzt2*Rqrd*sr*lam(4)*lam(7) -                                  &
               36*qzt2*Rqrd**2*sr*lam(4)*lam(7) -                               &
               36*qzt2*Rqsu*sr*lam(4)*lam(7) -                                  &
               72*qzt2*Rqrd*Rqsu*sr*lam(4)*lam(7) -                             &
               36*qzt2*Rqrd**2*Rqsu*sr*lam(4)*lam(7) +                          &
               18*lam(5)*lam(7) + 18*Rqsu*lam(5)*lam(7) -                       &
               36*qzt2*sr*lam(5)*lam(7) -                                       &
               72*qzt2*Rqrd*sr*lam(5)*lam(7) -                                  &
               36*qzt2*Rqrd**2*sr*lam(5)*lam(7) -                               &
               36*qzt2*Rqsu*sr*lam(5)*lam(7) -                                  &
               72*qzt2*Rqrd*Rqsu*sr*lam(5)*lam(7) -                             &
               36*qzt2*Rqrd**2*Rqsu*sr*lam(5)*lam(7) -                          &
               18*lam(6)*lam(7) - 18*Rqsu*lam(6)*lam(7) +                       &
               36*qzt2*sr*lam(6)*lam(7) +                                       &
               72*qzt2*Rqrd*sr*lam(6)*lam(7) +                                  &
               36*qzt2*Rqrd**2*sr*lam(6)*lam(7) +                               &
               36*qzt2*Rqsu*sr*lam(6)*lam(7) +                                  &
               72*qzt2*Rqrd*Rqsu*sr*lam(6)*lam(7) +                             &
               36*qzt2*Rqrd**2*Rqsu*sr*lam(6)*lam(7) - 9*lam(7)**2 -            &
               9*Rqsu*lam(7)**2 + 18*qzt2*sr*lam(7)**2 +                        &
               36*qzt2*Rqrd*sr*lam(7)**2 +                                      &
               18*qzt2*Rqrd**2*sr*lam(7)**2 +                                   &
               18*qzt2*Rqsu*sr*lam(7)**2 +                                      &
               36*qzt2*Rqrd*Rqsu*sr*lam(7)**2 +                                 &
               18*qzt2*Rqrd**2*Rqsu*sr*lam(7)**2 - 144*lam(1)*lam(8) -          &
               144*Rqsu*lam(1)*lam(8) + 288*qzt2*sr*lam(1)*lam(8) +             &
               576*qzt2*Rqrd*sr*lam(1)*lam(8) +                                 &
               288*qzt2*Rqrd**2*sr*lam(1)*lam(8) +                              &
               288*qzt2*Rqsu*sr*lam(1)*lam(8) +                                 &
               576*qzt2*Rqrd*Rqsu*sr*lam(1)*lam(8) +                            &
               288*qzt2*Rqrd**2*Rqsu*sr*lam(1)*lam(8) +                         &
               576*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(2)*            &
                lam(8) - 144*m2*sr*lam(4)*lam(8) -                              &
               144*m2*Rqrd*sr*lam(4)*lam(8) - 144*m2*sr*lam(5)*lam(8) -         &
               144*m2*Rqrd*sr*lam(5)*lam(8) + 144*m2*sr*lam(6)*lam(8) +         &
               144*m2*Rqrd*sr*lam(6)*lam(8) + 144*m2*sr*lam(7)*lam(8) +         &
               144*m2*Rqrd*sr*lam(7)*lam(8) - 144*lam(8)**2 -                   &
               144*Rqsu*lam(8)**2 + 288*qzt2*sr*lam(8)**2 +                     &
               576*qzt2*Rqrd*sr*lam(8)**2 +                                     &
               288*qzt2*Rqrd**2*sr*lam(8)**2 +                                  &
               288*qzt2*Rqsu*sr*lam(8)**2 +                                     &
               576*qzt2*Rqrd*Rqsu*sr*lam(8)**2 +                                &
               288*qzt2*Rqrd**2*Rqsu*sr*lam(8)**2 +                             &
               12*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(3)*             &
                (3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) + 14*lam(8)) -         &
               960*m2*(1 + Rqrd)*sr*lam(2)*lam(9) - 60*lam(4)*lam(9) -          &
               60*Rqsu*lam(4)*lam(9) + 120*qzt2*sr*lam(4)*lam(9) +              &
               240*qzt2*Rqrd*sr*lam(4)*lam(9) +                                 &
               120*qzt2*Rqrd**2*sr*lam(4)*lam(9) +                              &
               120*qzt2*Rqsu*sr*lam(4)*lam(9) +                                 &
               240*qzt2*Rqrd*Rqsu*sr*lam(4)*lam(9) +                            &
               120*qzt2*Rqrd**2*Rqsu*sr*lam(4)*lam(9) -                         &
               60*lam(5)*lam(9) - 60*Rqsu*lam(5)*lam(9) +                       &
               120*qzt2*sr*lam(5)*lam(9) +                                      &
               240*qzt2*Rqrd*sr*lam(5)*lam(9) +                                 &
               120*qzt2*Rqrd**2*sr*lam(5)*lam(9) +                              &
               120*qzt2*Rqsu*sr*lam(5)*lam(9) +                                 &
               240*qzt2*Rqrd*Rqsu*sr*lam(5)*lam(9) +                            &
               120*qzt2*Rqrd**2*Rqsu*sr*lam(5)*lam(9) +                         &
               60*lam(6)*lam(9) + 60*Rqsu*lam(6)*lam(9) -                       &
               120*qzt2*sr*lam(6)*lam(9) -                                      &
               240*qzt2*Rqrd*sr*lam(6)*lam(9) -                                 &
               120*qzt2*Rqrd**2*sr*lam(6)*lam(9) -                              &
               120*qzt2*Rqsu*sr*lam(6)*lam(9) -                                 &
               240*qzt2*Rqrd*Rqsu*sr*lam(6)*lam(9) -                            &
               120*qzt2*Rqrd**2*Rqsu*sr*lam(6)*lam(9) +                         &
               60*lam(7)*lam(9) + 60*Rqsu*lam(7)*lam(9) -                       &
               120*qzt2*sr*lam(7)*lam(9) -                                      &
               240*qzt2*Rqrd*sr*lam(7)*lam(9) -                                 &
               120*qzt2*Rqrd**2*sr*lam(7)*lam(9) -                              &
               120*qzt2*Rqsu*sr*lam(7)*lam(9) -                                 &
               240*qzt2*Rqrd*Rqsu*sr*lam(7)*lam(9) -                            &
               120*qzt2*Rqrd**2*Rqsu*sr*lam(7)*lam(9) -                         &
               480*m2*sr*lam(8)*lam(9) - 480*m2*Rqrd*sr*lam(8)*lam(9) -         &
               100*lam(9)**2 - 100*Rqsu*lam(9)**2 +                             &
               200*qzt2*sr*lam(9)**2 + 400*qzt2*Rqrd*sr*lam(9)**2 +             &
               200*qzt2*Rqrd**2*sr*lam(9)**2 +                                  &
               200*qzt2*Rqsu*sr*lam(9)**2 +                                     &
               400*qzt2*Rqrd*Rqsu*sr*lam(9)**2 +                                &
               200*qzt2*Rqrd**2*Rqsu*sr*lam(9)**2 +                             &
               24*(24*m2*(1 + Rqrd)*sr*(2*lam(2) + lam(8)) -                    &
                  (1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                    &
                   (3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                 &
                     10*lam(9)))*lam(10) +                                      &
               144*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(10)**2))       &
+ 2*dRqrd*quB*qxy2*(1 + Rqrd)*                                                  &
          (expuISs*qxy2*(1 + Rqrd)*(1 + Rqsu)*                                  &
             (576*lam(2)**2 + 576*lam(2)*lam(8) +                               &
               144*lam(8)*(-lam(1) + lam(8)) -                                  &
               12*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +               &
                  14*lam(8)) +                                                  &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -         &
                  12*lam(10))**2) +                                             &
            expuISs*qzt2*(1 + Rqrd)*(1 + Rqsu)*                                 &
             (576*lam(2)**2 + 576*lam(2)*lam(8) +                               &
               144*lam(8)*(-lam(1) + lam(8)) -                                  &
               12*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +               &
                  14*lam(8)) +                                                  &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -         &
                  12*lam(10))**2) +                                             &
            2*expdISs*(qxy2 + qzt2)*(1 + Rqrd)*(1 + Rqsd)*                      &
             (576*lam(2)**2 - 18*lam(3)*lam(4) + 9*lam(4)**2 -                  &
               18*lam(3)*lam(5) - 36*lam(4)*lam(5) + 63*lam(5)**2 -             &
               18*lam(3)*lam(6) - 18*lam(4)*lam(6) + 36*lam(5)*lam(6) +         &
               9*lam(6)**2 - 18*lam(3)*lam(7) + 36*lam(4)*lam(7) -              &
               126*lam(5)*lam(7) - 36*lam(6)*lam(7) + 63*lam(7)**2 +            &
               288*lam(2)*(2*lam(3) - lam(8)) + 96*lam(3)*lam(8) +              &
               144*lam(8)**2 +                                                  &
               108*lam(1)*(28*lam(2) - lam(4) - 3*lam(5) - lam(6) -             &
                  3*lam(7) + 6*lam(8)) - 30*lam(4)*lam(9) +                     &
               150*lam(5)*lam(9) + 30*lam(6)*lam(9) -                           &
               150*lam(7)*lam(9) + 100*lam(9)**2 +                              &
               18*(-5*lam(4) + 11*lam(5) + 5*lam(6) - 11*lam(7) +               &
                  10*lam(9))*lam(10) + 228*lam(10)**2) +                        &
            24*expuISs*m2*(12*lam(1)*lam(3) + 7*lam(3)**2 +                     &
               lam(8)*(-6*lam(5) + 6*lam(6) + 7*lam(8) - 10*lam(9) +            &
                  12*lam(10)) +                                                 &
               lam(2)*(-6*lam(4) - 6*lam(5) + 6*lam(6) + 6*lam(7) -             &
                  20*lam(9) + 24*lam(10))) +                                    &
            48*expdISs*m2*(108*lam(2)**2 +                                      &
               (6*lam(1) + lam(3))*(24*lam(1) + 5*lam(3)) +                     &
               lam(8)*(-15*lam(5) - 3*lam(6) + 7*lam(8) - 10*lam(9) -           &
                  9*lam(10)) +                                                  &
               lam(2)*(-15*lam(4) - 15*lam(5) - 3*lam(6) - 39*lam(7) +          &
                  48*lam(8) + 10*lam(9) + 30*lam(10))))) +                      &
      expdISr*qxy2*tanhdr**2*(-36*dRqsd*expdISs*quB*qxy2*(1 + Rqrd)*            &
          (1 + Rqsd)**2*tanhds**3*                                              &
          (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) -             &
         18*dRqsu*expuISs*qdB*qxy2*(1 + Rqrd)*(1 + Rqsu)**2*tanhus**3*          &
          (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +             &
         2*expuISs*(1 + Rqrd)*tanhus*                                           &
          (9*dRqsu*qdB*qxy2*(1 + Rqsu)**2*                                      &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            dRqrd*quB*(qzt2*(1 + Rqrd)*(1 + Rqsu)*                              &
                (576*lam(2)**2 + 576*lam(2)*lam(8) +                            &
                  144*lam(8)*(lam(1) + lam(8)) +                                &
                  12*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +            &
                     14*lam(8)) +                                               &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) -                              &
               24*m2*(2*lam(2) + lam(8))*                                       &
                (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -        &
                  12*lam(10)))) +                                               &
         expuISs*qdB*quB*tanhus**2*                                             &
          (-9*dRqsu*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                  &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            dRqrd*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                    &
             (576*lam(2)**2 + 576*lam(2)*lam(8) +                               &
               144*lam(8)*(-lam(1) + lam(8)) -                                  &
               12*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +               &
                  14*lam(8)) +                                                  &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -         &
                  12*lam(10))**2)) +                                            &
         2*expdISs*qdB*quB*tanhds**2*                                           &
          (-9*dRqsd*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                  &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            dRqrd*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                    &
             (576*lam(2)**2 - 18*lam(3)*lam(4) + 9*lam(4)**2 -                  &
               18*lam(3)*lam(5) - 36*lam(4)*lam(5) + 63*lam(5)**2 -             &
               18*lam(3)*lam(6) - 18*lam(4)*lam(6) + 36*lam(5)*lam(6) +         &
               9*lam(6)**2 - 18*lam(3)*lam(7) + 36*lam(4)*lam(7) -              &
               126*lam(5)*lam(7) - 36*lam(6)*lam(7) + 63*lam(7)**2 +            &
               288*lam(2)*(2*lam(3) - lam(8)) + 96*lam(3)*lam(8) +              &
               144*lam(8)**2 +                                                  &
               108*lam(1)*(28*lam(2) - lam(4) - 3*lam(5) - lam(6) -             &
                  3*lam(7) + 6*lam(8)) - 30*lam(4)*lam(9) +                     &
               150*lam(5)*lam(9) + 30*lam(6)*lam(9) -                           &
               150*lam(7)*lam(9) + 100*lam(9)**2 +                              &
               18*(-5*lam(4) + 11*lam(5) + 5*lam(6) - 11*lam(7) +               &
                  10*lam(9))*lam(10) + 228*lam(10)**2)) +                       &
         qdB*quB*(9*(1 + Rqrd)*                                                 &
             (2*dRqsd*expdISs*(-1 + 2*qzt2*(1 + Rqsd)**2*ss) +                  &
               dRqsu*expuISs*(-1 + 2*qzt2*(1 + Rqsu)**2*ss))*                   &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) -          &
            dRqrd*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                               &
             (expuISs*(1 + Rqsu)*                                               &
                (576*lam(2)**2 + 576*lam(2)*lam(8) +                            &
                  144*lam(8)*(-lam(1) + lam(8)) -                               &
                  12*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +            &
                     14*lam(8)) +                                               &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) +                              &
               2*expdISs*(1 + Rqsd)*                                            &
                (576*lam(2)**2 - 18*lam(3)*lam(4) + 9*lam(4)**2 -               &
                  18*lam(3)*lam(5) - 36*lam(4)*lam(5) + 63*lam(5)**2 -          &
                  18*lam(3)*lam(6) - 18*lam(4)*lam(6) +                         &
                  36*lam(5)*lam(6) + 9*lam(6)**2 - 18*lam(3)*lam(7) +           &
                  36*lam(4)*lam(7) - 126*lam(5)*lam(7) -                        &
                  36*lam(6)*lam(7) + 63*lam(7)**2 +                             &
                  288*lam(2)*(2*lam(3) - lam(8)) + 96*lam(3)*lam(8) +           &
                  144*lam(8)**2 +                                               &
                  108*lam(1)*                                                   &
                   (28*lam(2) - lam(4) - 3*lam(5) - lam(6) - 3*lam(7) +         &
                     6*lam(8)) - 30*lam(4)*lam(9) + 150*lam(5)*lam(9) +         &
                  30*lam(6)*lam(9) - 150*lam(7)*lam(9) + 100*lam(9)**2 +        &
                  18*(-5*lam(4) + 11*lam(5) + 5*lam(6) - 11*lam(7) +            &
                     10*lam(9))*lam(10) + 228*lam(10)**2))) -                   &
         4*expdISs*quB*(1 + Rqrd)*tanhds*                                       &
          (-9*dRqsd*qxy2*(1 + Rqsd)**2*                                         &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            dRqrd*(qzt2*(1 + Rqrd)*(1 + Rqsd)*                                  &
                (-576*lam(2)**2 - 18*lam(3)*lam(4) - 9*lam(4)**2 -              &
                  18*lam(3)*lam(5) + 36*lam(4)*lam(5) - 63*lam(5)**2 -          &
                  18*lam(3)*lam(6) + 18*lam(4)*lam(6) -                         &
                  36*lam(5)*lam(6) - 9*lam(6)**2 - 18*lam(3)*lam(7) -           &
                  36*lam(4)*lam(7) + 126*lam(5)*lam(7) +                        &
                  36*lam(6)*lam(7) - 63*lam(7)**2 + 96*lam(3)*lam(8) -          &
                  144*lam(8)**2 + 288*lam(2)*(2*lam(3) + lam(8)) +              &
                  108*lam(1)*                                                   &
                   (28*lam(2) - lam(4) - 3*lam(5) - lam(6) - 3*lam(7) +         &
                     6*lam(8)) + 30*lam(4)*lam(9) - 150*lam(5)*lam(9) -         &
                  30*lam(6)*lam(9) + 150*lam(7)*lam(9) - 100*lam(9)**2 +        &
                  18*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -             &
                     10*lam(9))*lam(10) - 228*lam(10)**2) +                     &
               12*m2*(lam(8)*(-3*lam(4) + 15*lam(5) + 3*lam(6) -                &
                     15*lam(7) + 20*lam(9) + 18*lam(10)) +                      &
                  4*lam(2)*(3*lam(4) - 6*lam(5) - 3*lam(6) + 6*lam(7) -         &
                     5*(lam(9) + 3*lam(10))))))) +                              &
      expdISr*expuISs*qdB*(18*dRqsu*qxy2**2*(1 + Rqrd)*(1 + Rqsu)**2*           &
          tanhus**3*(16*lam(2)**2 +                                             &
            (lam(4) + lam(5) + lam(6) + lam(7))**2) -                           &
         18*dRqsu*qxy2*(1 + Rqsu)*tanhus*                                       &
          (16*m2*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                   &
            qxy2*(1 + Rqrd)*(1 + Rqsu)*                                         &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            qzt2*(1 + Rqrd)*(1 + Rqsu)*                                         &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2)) +         &
         quB*qxy2*tanhus**2*(9*dRqsu*(1 + Rqrd)*                                &
             (-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                                    &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) -          &
            dRqrd*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                    &
             (576*lam(2)**2 + 576*lam(2)*lam(8) +                               &
               144*lam(8)*(-lam(1) + lam(8)) -                                  &
               12*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +               &
                  14*lam(8)) +                                                  &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -         &
                  12*lam(10))**2)) +                                            &
         quB*(-9*dRqsu*(qxy2*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*         &
                (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2)         &
+ qzt2*(16*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*lam(2)**2 +                &
                  32*m2*(1 + Rqsu)*ss*lam(2)*                                   &
                   (lam(4) + lam(5) + lam(6) + lam(7)) +                        &
                  (1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                    &
                   (lam(4) + lam(5) + lam(6) + lam(7))**2)) +                   &
            dRqrd*(qxy2*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*              &
                (576*lam(2)**2 + 576*lam(2)*lam(8) +                            &
                  144*lam(8)*(-lam(1) + lam(8)) -                               &
                  12*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +            &
                     14*lam(8)) +                                               &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) +                              &
               qzt2*(576*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*             &
                   lam(2)**2 +                                                  &
                  96*lam(2)*(6*(1 + Rqsu)*                                      &
                      (-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(8) -                   &
                     m2*(1 + Rqrd)*sr*                                          &
                      (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +              &
                        10*lam(9) - 12*lam(10))) -                              &
                  (1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                    &
                   (144*(lam(1) - lam(8))*lam(8) +                              &
                     12*lam(3)*                                                 &
                      (3*(lam(4) - lam(5) + lam(6) - lam(7)) +                  &
                        14*lam(8)) -                                            &
                     (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +               &
                        10*lam(9) - 12*lam(10))**2) +                           &
                  48*m2*(1 + Rqrd)*sr*                                          &
                   (12*lam(1)*lam(3) + 7*lam(3)**2 +                            &
                     lam(8)*(-6*lam(5) + 6*lam(6) + 7*lam(8) -                  &
                        10*lam(9) + 12*lam(10))))))) +                          &
      expdISs*expuISr*quB*(18*dRqsd*qxy2**2*(1 + Rqru)*(1 + Rqsd)**2*           &
          tanhds**3*(16*lam(2)**2 +                                             &
            (lam(4) + lam(5) + lam(6) + lam(7))**2) -                           &
         18*dRqsd*qxy2*(1 + Rqsd)*tanhds*                                       &
          (16*m2*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                   &
            qxy2*(1 + Rqru)*(1 + Rqsd)*                                         &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            qzt2*(1 + Rqru)*(1 + Rqsd)*                                         &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2)) +         &
         qdB*qxy2*tanhds**2*(9*dRqsd*(1 + Rqru)*                                &
             (-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                                    &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) -          &
            dRqru*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                    &
             (576*lam(2)**2 + 576*lam(2)*lam(8) +                               &
               144*lam(8)*(-lam(1) + lam(8)) -                                  &
               12*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +               &
                  14*lam(8)) +                                                  &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -         &
                  12*lam(10))**2)) +                                            &
         qdB*(-9*dRqsd*(qxy2*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*         &
                (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2)         &
+ qzt2*(16*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*lam(2)**2 +                &
                  32*m2*(1 + Rqsd)*ss*lam(2)*                                   &
                   (lam(4) + lam(5) + lam(6) + lam(7)) +                        &
                  (1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                    &
                   (lam(4) + lam(5) + lam(6) + lam(7))**2)) +                   &
            dRqru*(qxy2*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*              &
                (576*lam(2)**2 + 576*lam(2)*lam(8) +                            &
                  144*lam(8)*(-lam(1) + lam(8)) -                               &
                  12*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +            &
                     14*lam(8)) +                                               &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) +                              &
               qzt2*(576*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*             &
                   lam(2)**2 +                                                  &
                  96*lam(2)*(6*(1 + Rqsd)*                                      &
                      (-1 + 2*qzt2*(1 + Rqru)**2*sr)*lam(8) -                   &
                     m2*(1 + Rqru)*sr*                                          &
                      (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +              &
                        10*lam(9) - 12*lam(10))) -                              &
                  (1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                    &
                   (144*(lam(1) - lam(8))*lam(8) +                              &
                     12*lam(3)*                                                 &
                      (3*(lam(4) - lam(5) + lam(6) - lam(7)) +                  &
                        14*lam(8)) -                                            &
                     (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +               &
                        10*lam(9) - 12*lam(10))**2) +                           &
                  48*m2*(1 + Rqru)*sr*                                          &
                   (12*lam(1)*lam(3) + 7*lam(3)**2 +                            &
                     lam(8)*(-6*lam(5) + 6*lam(6) + 7*lam(8) -                  &
                        10*lam(9) + 12*lam(10))))))) +                          &
      expuISr*tanhur*(-2*expdISs*qxy2*(1 + Rqru)*(1 + Rqsd)*tanhds**2*          &
          (9*dRqsd*quB*qzt2*(1 + Rqsd)*                                         &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            dRqru*qdB*qxy2*(1 + Rqru)*                                          &
             (576*lam(2)**2 + 576*lam(2)*lam(8) +                               &
               144*lam(8)*(-lam(1) + lam(8)) -                                  &
               12*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +               &
                  14*lam(8)) +                                                  &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -         &
                  12*lam(10))**2)) -                                            &
         expdISs*qdB*quB*qzt2*tanhds*                                           &
          (9*dRqsd*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                   &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) -          &
            dRqru*(576*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*               &
                lam(2)**2 - 288*m2*(1 + Rqru)*sr*lam(2)*lam(4) -                &
               9*lam(4)**2 - 9*Rqsd*lam(4)**2 + 18*qzt2*sr*lam(4)**2 +          &
               36*qzt2*Rqru*sr*lam(4)**2 +                                      &
               18*qzt2*Rqru**2*sr*lam(4)**2 +                                   &
               18*qzt2*Rqsd*sr*lam(4)**2 +                                      &
               36*qzt2*Rqru*Rqsd*sr*lam(4)**2 +                                 &
               18*qzt2*Rqru**2*Rqsd*sr*lam(4)**2 -                              &
               288*m2*(1 + Rqru)*sr*lam(2)*lam(5) - 18*lam(4)*lam(5) -          &
               18*Rqsd*lam(4)*lam(5) + 36*qzt2*sr*lam(4)*lam(5) +               &
               72*qzt2*Rqru*sr*lam(4)*lam(5) +                                  &
               36*qzt2*Rqru**2*sr*lam(4)*lam(5) +                               &
               36*qzt2*Rqsd*sr*lam(4)*lam(5) +                                  &
               72*qzt2*Rqru*Rqsd*sr*lam(4)*lam(5) +                             &
               36*qzt2*Rqru**2*Rqsd*sr*lam(4)*lam(5) - 9*lam(5)**2 -            &
               9*Rqsd*lam(5)**2 + 18*qzt2*sr*lam(5)**2 +                        &
               36*qzt2*Rqru*sr*lam(5)**2 +                                      &
               18*qzt2*Rqru**2*sr*lam(5)**2 +                                   &
               18*qzt2*Rqsd*sr*lam(5)**2 +                                      &
               36*qzt2*Rqru*Rqsd*sr*lam(5)**2 +                                 &
               18*qzt2*Rqru**2*Rqsd*sr*lam(5)**2 +                              &
               288*m2*(1 + Rqru)*sr*lam(2)*lam(6) + 18*lam(4)*lam(6) +          &
               18*Rqsd*lam(4)*lam(6) - 36*qzt2*sr*lam(4)*lam(6) -               &
               72*qzt2*Rqru*sr*lam(4)*lam(6) -                                  &
               36*qzt2*Rqru**2*sr*lam(4)*lam(6) -                               &
               36*qzt2*Rqsd*sr*lam(4)*lam(6) -                                  &
               72*qzt2*Rqru*Rqsd*sr*lam(4)*lam(6) -                             &
               36*qzt2*Rqru**2*Rqsd*sr*lam(4)*lam(6) +                          &
               18*lam(5)*lam(6) + 18*Rqsd*lam(5)*lam(6) -                       &
               36*qzt2*sr*lam(5)*lam(6) -                                       &
               72*qzt2*Rqru*sr*lam(5)*lam(6) -                                  &
               36*qzt2*Rqru**2*sr*lam(5)*lam(6) -                               &
               36*qzt2*Rqsd*sr*lam(5)*lam(6) -                                  &
               72*qzt2*Rqru*Rqsd*sr*lam(5)*lam(6) -                             &
               36*qzt2*Rqru**2*Rqsd*sr*lam(5)*lam(6) - 9*lam(6)**2 -            &
               9*Rqsd*lam(6)**2 + 18*qzt2*sr*lam(6)**2 +                        &
               36*qzt2*Rqru*sr*lam(6)**2 +                                      &
               18*qzt2*Rqru**2*sr*lam(6)**2 +                                   &
               18*qzt2*Rqsd*sr*lam(6)**2 +                                      &
               36*qzt2*Rqru*Rqsd*sr*lam(6)**2 +                                 &
               18*qzt2*Rqru**2*Rqsd*sr*lam(6)**2 +                              &
               288*m2*(1 + Rqru)*sr*lam(2)*lam(7) + 18*lam(4)*lam(7) +          &
               18*Rqsd*lam(4)*lam(7) - 36*qzt2*sr*lam(4)*lam(7) -               &
               72*qzt2*Rqru*sr*lam(4)*lam(7) -                                  &
               36*qzt2*Rqru**2*sr*lam(4)*lam(7) -                               &
               36*qzt2*Rqsd*sr*lam(4)*lam(7) -                                  &
               72*qzt2*Rqru*Rqsd*sr*lam(4)*lam(7) -                             &
               36*qzt2*Rqru**2*Rqsd*sr*lam(4)*lam(7) +                          &
               18*lam(5)*lam(7) + 18*Rqsd*lam(5)*lam(7) -                       &
               36*qzt2*sr*lam(5)*lam(7) -                                       &
               72*qzt2*Rqru*sr*lam(5)*lam(7) -                                  &
               36*qzt2*Rqru**2*sr*lam(5)*lam(7) -                               &
               36*qzt2*Rqsd*sr*lam(5)*lam(7) -                                  &
               72*qzt2*Rqru*Rqsd*sr*lam(5)*lam(7) -                             &
               36*qzt2*Rqru**2*Rqsd*sr*lam(5)*lam(7) -                          &
               18*lam(6)*lam(7) - 18*Rqsd*lam(6)*lam(7) +                       &
               36*qzt2*sr*lam(6)*lam(7) +                                       &
               72*qzt2*Rqru*sr*lam(6)*lam(7) +                                  &
               36*qzt2*Rqru**2*sr*lam(6)*lam(7) +                               &
               36*qzt2*Rqsd*sr*lam(6)*lam(7) +                                  &
               72*qzt2*Rqru*Rqsd*sr*lam(6)*lam(7) +                             &
               36*qzt2*Rqru**2*Rqsd*sr*lam(6)*lam(7) - 9*lam(7)**2 -            &
               9*Rqsd*lam(7)**2 + 18*qzt2*sr*lam(7)**2 +                        &
               36*qzt2*Rqru*sr*lam(7)**2 +                                      &
               18*qzt2*Rqru**2*sr*lam(7)**2 +                                   &
               18*qzt2*Rqsd*sr*lam(7)**2 +                                      &
               36*qzt2*Rqru*Rqsd*sr*lam(7)**2 +                                 &
               18*qzt2*Rqru**2*Rqsd*sr*lam(7)**2 - 144*lam(1)*lam(8) -          &
               144*Rqsd*lam(1)*lam(8) + 288*qzt2*sr*lam(1)*lam(8) +             &
               576*qzt2*Rqru*sr*lam(1)*lam(8) +                                 &
               288*qzt2*Rqru**2*sr*lam(1)*lam(8) +                              &
               288*qzt2*Rqsd*sr*lam(1)*lam(8) +                                 &
               576*qzt2*Rqru*Rqsd*sr*lam(1)*lam(8) +                            &
               288*qzt2*Rqru**2*Rqsd*sr*lam(1)*lam(8) +                         &
               576*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*lam(2)*            &
                lam(8) - 144*m2*sr*lam(4)*lam(8) -                              &
               144*m2*Rqru*sr*lam(4)*lam(8) - 144*m2*sr*lam(5)*lam(8) -         &
               144*m2*Rqru*sr*lam(5)*lam(8) + 144*m2*sr*lam(6)*lam(8) +         &
               144*m2*Rqru*sr*lam(6)*lam(8) + 144*m2*sr*lam(7)*lam(8) +         &
               144*m2*Rqru*sr*lam(7)*lam(8) - 144*lam(8)**2 -                   &
               144*Rqsd*lam(8)**2 + 288*qzt2*sr*lam(8)**2 +                     &
               576*qzt2*Rqru*sr*lam(8)**2 +                                     &
               288*qzt2*Rqru**2*sr*lam(8)**2 +                                  &
               288*qzt2*Rqsd*sr*lam(8)**2 +                                     &
               576*qzt2*Rqru*Rqsd*sr*lam(8)**2 +                                &
               288*qzt2*Rqru**2*Rqsd*sr*lam(8)**2 +                             &
               12*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*lam(3)*             &
                (3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) + 14*lam(8)) -         &
               960*m2*(1 + Rqru)*sr*lam(2)*lam(9) - 60*lam(4)*lam(9) -          &
               60*Rqsd*lam(4)*lam(9) + 120*qzt2*sr*lam(4)*lam(9) +              &
               240*qzt2*Rqru*sr*lam(4)*lam(9) +                                 &
               120*qzt2*Rqru**2*sr*lam(4)*lam(9) +                              &
               120*qzt2*Rqsd*sr*lam(4)*lam(9) +                                 &
               240*qzt2*Rqru*Rqsd*sr*lam(4)*lam(9) +                            &
               120*qzt2*Rqru**2*Rqsd*sr*lam(4)*lam(9) -                         &
               60*lam(5)*lam(9) - 60*Rqsd*lam(5)*lam(9) +                       &
               120*qzt2*sr*lam(5)*lam(9) +                                      &
               240*qzt2*Rqru*sr*lam(5)*lam(9) +                                 &
               120*qzt2*Rqru**2*sr*lam(5)*lam(9) +                              &
               120*qzt2*Rqsd*sr*lam(5)*lam(9) +                                 &
               240*qzt2*Rqru*Rqsd*sr*lam(5)*lam(9) +                            &
               120*qzt2*Rqru**2*Rqsd*sr*lam(5)*lam(9) +                         &
               60*lam(6)*lam(9) + 60*Rqsd*lam(6)*lam(9) -                       &
               120*qzt2*sr*lam(6)*lam(9) -                                      &
               240*qzt2*Rqru*sr*lam(6)*lam(9) -                                 &
               120*qzt2*Rqru**2*sr*lam(6)*lam(9) -                              &
               120*qzt2*Rqsd*sr*lam(6)*lam(9) -                                 &
               240*qzt2*Rqru*Rqsd*sr*lam(6)*lam(9) -                            &
               120*qzt2*Rqru**2*Rqsd*sr*lam(6)*lam(9) +                         &
               60*lam(7)*lam(9) + 60*Rqsd*lam(7)*lam(9) -                       &
               120*qzt2*sr*lam(7)*lam(9) -                                      &
               240*qzt2*Rqru*sr*lam(7)*lam(9) -                                 &
               120*qzt2*Rqru**2*sr*lam(7)*lam(9) -                              &
               120*qzt2*Rqsd*sr*lam(7)*lam(9) -                                 &
               240*qzt2*Rqru*Rqsd*sr*lam(7)*lam(9) -                            &
               120*qzt2*Rqru**2*Rqsd*sr*lam(7)*lam(9) -                         &
               480*m2*sr*lam(8)*lam(9) - 480*m2*Rqru*sr*lam(8)*lam(9) -         &
               100*lam(9)**2 - 100*Rqsd*lam(9)**2 +                             &
               200*qzt2*sr*lam(9)**2 + 400*qzt2*Rqru*sr*lam(9)**2 +             &
               200*qzt2*Rqru**2*sr*lam(9)**2 +                                  &
               200*qzt2*Rqsd*sr*lam(9)**2 +                                     &
               400*qzt2*Rqru*Rqsd*sr*lam(9)**2 +                                &
               200*qzt2*Rqru**2*Rqsd*sr*lam(9)**2 +                             &
               24*(24*m2*(1 + Rqru)*sr*(2*lam(2) + lam(8)) -                    &
                  (1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                    &
                   (3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                 &
                     10*lam(9)))*lam(10) +                                      &
               144*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*lam(10)**2))       &
+ 2*qdB*(-9*dRqsu*expuISs*qzt2*(1 + Rqru)*tanhus*                               &
             (quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +                              &
               2*qxy2*(1 + Rqsu)**2*tanhus)*                                    &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            dRqru*(expdISs*qxy2*(1 + Rqru)*                                     &
                (qxy2*(1 + Rqru)*(1 + Rqsd)*                                    &
                   (576*lam(2)**2 + 576*lam(2)*lam(8) +                         &
                     144*lam(8)*(-lam(1) + lam(8)) -                            &
                     12*lam(3)*                                                 &
                      (3*(lam(4) - lam(5) + lam(6) - lam(7)) +                  &
                        14*lam(8)) +                                            &
                     (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +               &
                        10*lam(9) - 12*lam(10))**2) +                           &
                  qzt2*(1 + Rqru)*(1 + Rqsd)*                                   &
                   (576*lam(2)**2 + 576*lam(2)*lam(8) +                         &
                     144*lam(8)*(-lam(1) + lam(8)) -                            &
                     12*lam(3)*                                                 &
                      (3*(lam(4) - lam(5) + lam(6) - lam(7)) +                  &
                        14*lam(8)) +                                            &
                     (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +               &
                        10*lam(9) - 12*lam(10))**2) +                           &
                  24*m2*(12*lam(1)*lam(3) + 7*lam(3)**2 +                       &
                     lam(8)*(-6*lam(5) + 6*lam(6) + 7*lam(8) -                  &
                        10*lam(9) + 12*lam(10)) +                               &
                     lam(2)*(-6*lam(4) - 6*lam(5) + 6*lam(6) +                  &
                        6*lam(7) - 20*lam(9) + 24*lam(10)))) +                  &
               expuISs*(-((1 + Rqsu)*                                           &
                     (quB*qzt2*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*tanhus*           &
                        (-576*lam(2)**2 - 18*lam(3)*lam(4) -                    &
                         9*lam(4)**2 - 18*lam(3)*lam(5) +                       &
                         36*lam(4)*lam(5) - 63*lam(5)**2 -                      &
                         18*lam(3)*lam(6) + 18*lam(4)*lam(6) -                  &
                         36*lam(5)*lam(6) - 9*lam(6)**2 -                       &
                         18*lam(3)*lam(7) - 36*lam(4)*lam(7) +                  &
                         126*lam(5)*lam(7) + 36*lam(6)*lam(7) -                 &
                         63*lam(7)**2 + 96*lam(3)*lam(8) -                      &
                         144*lam(8)**2 +                                        &
                         288*lam(2)*(2*lam(3) + lam(8)) +                       &
                         108*lam(1)*                                            &
                        (28*lam(2) - lam(4) - 3*lam(5) - lam(6) -               &
                        3*lam(7) + 6*lam(8)) + 30*lam(4)*lam(9) -               &
                         150*lam(5)*lam(9) - 30*lam(6)*lam(9) +                 &
                         150*lam(7)*lam(9) - 100*lam(9)**2 +                    &
                         18*(5*lam(4) - 11*lam(5) - 5*lam(6) +                  &
                        11*lam(7) - 10*lam(9))*lam(10) - 228*lam(10)**2)        &
- 2*qxy2*qzt2*(1 + Rqru)**2*(576*lam(2)**2 - 18*lam(3)*lam(4) +                 &
                         9*lam(4)**2 - 18*lam(3)*lam(5) -                       &
                         36*lam(4)*lam(5) + 63*lam(5)**2 -                      &
                         18*lam(3)*lam(6) - 18*lam(4)*lam(6) +                  &
                         36*lam(5)*lam(6) + 9*lam(6)**2 -                       &
                         18*lam(3)*lam(7) + 36*lam(4)*lam(7) -                  &
                         126*lam(5)*lam(7) - 36*lam(6)*lam(7) +                 &
                         63*lam(7)**2 +                                         &
                         288*lam(2)*(2*lam(3) - lam(8)) +                       &
                         96*lam(3)*lam(8) + 144*lam(8)**2 +                     &
                         108*lam(1)*                                            &
                        (28*lam(2) - lam(4) - 3*lam(5) - lam(6) -               &
                        3*lam(7) + 6*lam(8)) - 30*lam(4)*lam(9) +               &
                         150*lam(5)*lam(9) + 30*lam(6)*lam(9) -                 &
                         150*lam(7)*lam(9) + 100*lam(9)**2 +                    &
                         18*(-5*lam(4) + 11*lam(5) + 5*lam(6) -                 &
                        11*lam(7) + 10*lam(9))*lam(10) + 228*lam(10)**2)        &
+ 2*qxy2**2*(1 + Rqru)**2*(-1 + tanhus**2)*                                     &
                        (576*lam(2)**2 - 18*lam(3)*lam(4) +                     &
                         9*lam(4)**2 - 18*lam(3)*lam(5) -                       &
                         36*lam(4)*lam(5) + 63*lam(5)**2 -                      &
                         18*lam(3)*lam(6) - 18*lam(4)*lam(6) +                  &
                         36*lam(5)*lam(6) + 9*lam(6)**2 -                       &
                         18*lam(3)*lam(7) + 36*lam(4)*lam(7) -                  &
                         126*lam(5)*lam(7) - 36*lam(6)*lam(7) +                 &
                         63*lam(7)**2 + 288*lam(2)*(2*lam(3) - lam(8)) +        &
                         96*lam(3)*lam(8) + 144*lam(8)**2 +                     &
                         108*lam(1)*                                            &
                         (28*lam(2) - lam(4) - 3*lam(5) - lam(6) -              &
                         3*lam(7) + 6*lam(8)) - 30*lam(4)*lam(9) +              &
                         150*lam(5)*lam(9) + 30*lam(6)*lam(9) -                 &
                         150*lam(7)*lam(9) + 100*lam(9)**2 +                    &
                         18*(-5*lam(4) + 11*lam(5) + 5*lam(6) -                 &
                        11*lam(7) + 10*lam(9))*lam(10) + 228*lam(10)**2))       &
) + 24*m2*(1 + Rqru)*(quB*qzt2*sr*tanhus*                                       &
                      (lam(8)*                                                  &
                         (3*lam(4) - 15*lam(5) - 3*lam(6) + 15*lam(7) -         &
                         20*lam(9) - 18*lam(10)) +                              &
                        4*lam(2)*                                               &
                         (-3*lam(4) + 6*lam(5) + 3*lam(6) - 6*lam(7) +          &
                         5*lam(9) + 15*lam(10))) +                              &
                     2*qxy2*(108*lam(2)**2 +                                    &
                        (6*lam(1) + lam(3))*(24*lam(1) + 5*lam(3)) +            &
                        lam(8)*                                                 &
                         (-15*lam(5) - 3*lam(6) + 7*lam(8) - 10*lam(9) -        &
                         9*lam(10)) +                                           &
                        lam(2)*                                                 &
                         (-15*lam(4) - 15*lam(5) - 3*lam(6) - 39*lam(7) +       &
                         48*lam(8) + 10*lam(9) + 30*lam(10)))))))) +            &
      expuISr*qxy2*tanhur**2*(-18*dRqsd*expdISs*quB*qxy2*(1 + Rqru)*            &
          (1 + Rqsd)**2*tanhds**3*                                              &
          (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +             &
         2*expdISs*(1 + Rqru)*tanhds*                                           &
          (9*dRqsd*quB*qxy2*(1 + Rqsd)**2*                                      &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            dRqru*qdB*(qzt2*(1 + Rqru)*(1 + Rqsd)*                              &
                (576*lam(2)**2 + 576*lam(2)*lam(8) +                            &
                  144*lam(8)*(lam(1) + lam(8)) +                                &
                  12*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +            &
                     14*lam(8)) +                                               &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) -                              &
               24*m2*(2*lam(2) + lam(8))*                                       &
                (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -        &
                  12*lam(10)))) +                                               &
         expdISs*qdB*quB*tanhds**2*                                             &
          (-9*dRqsd*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                  &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            dRqru*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                    &
             (576*lam(2)**2 + 576*lam(2)*lam(8) +                               &
               144*lam(8)*(-lam(1) + lam(8)) -                                  &
               12*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +               &
                  14*lam(8)) +                                                  &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -         &
                  12*lam(10))**2)) +                                            &
         qdB*(9*(1 + Rqru)*(dRqsd*expdISs*quB*                                  &
                (-1 + 2*qzt2*(1 + Rqsd)**2*ss) -                                &
               2*dRqsu*expuISs*                                                 &
                (quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +                           &
                  2*qxy2*(1 + Rqsu)**2*tanhus)*(-1 + tanhus**2))*               &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            dRqru*(-(expdISs*quB*(1 + Rqsd)*                                    &
                  (-1 + 2*qzt2*(1 + Rqru)**2*sr)*                               &
                  (576*lam(2)**2 + 576*lam(2)*lam(8) +                          &
                    144*lam(8)*(-lam(1) + lam(8)) -                             &
                    12*lam(3)*                                                  &
                     (3*(lam(4) - lam(5) + lam(6) - lam(7)) + 14*lam(8))        &
+ (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) - 12*lam(10))**2)      &
) + 2*expuISs*(quB*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                   &
                   (-1 + tanhus**2)*                                            &
                   (576*lam(2)**2 - 18*lam(3)*lam(4) + 9*lam(4)**2 -            &
                     18*lam(3)*lam(5) - 36*lam(4)*lam(5) + 63*lam(5)**2 -       &
                     18*lam(3)*lam(6) - 18*lam(4)*lam(6) +                      &
                     36*lam(5)*lam(6) + 9*lam(6)**2 - 18*lam(3)*lam(7) +        &
                     36*lam(4)*lam(7) - 126*lam(5)*lam(7) -                     &
                     36*lam(6)*lam(7) + 63*lam(7)**2 +                          &
                     288*lam(2)*(2*lam(3) - lam(8)) + 96*lam(3)*lam(8) +        &
                     144*lam(8)**2 +                                            &
                     108*lam(1)*                                                &
                      (28*lam(2) - lam(4) - 3*lam(5) - lam(6) -                 &
                        3*lam(7) + 6*lam(8)) - 30*lam(4)*lam(9) +               &
                     150*lam(5)*lam(9) + 30*lam(6)*lam(9) -                     &
                     150*lam(7)*lam(9) + 100*lam(9)**2 +                        &
                     18*(-5*lam(4) + 11*lam(5) + 5*lam(6) - 11*lam(7) +         &
                        10*lam(9))*lam(10) + 228*lam(10)**2) -                  &
                  2*(1 + Rqru)*tanhus*                                          &
                   (qzt2*(1 + Rqru)*(1 + Rqsu)*                                 &
                      (-576*lam(2)**2 - 18*lam(3)*lam(4) - 9*lam(4)**2 -        &
                        18*lam(3)*lam(5) + 36*lam(4)*lam(5) -                   &
                        63*lam(5)**2 - 18*lam(3)*lam(6) +                       &
                        18*lam(4)*lam(6) - 36*lam(5)*lam(6) -                   &
                        9*lam(6)**2 - 18*lam(3)*lam(7) -                        &
                        36*lam(4)*lam(7) + 126*lam(5)*lam(7) +                  &
                        36*lam(6)*lam(7) - 63*lam(7)**2 +                       &
                        96*lam(3)*lam(8) - 144*lam(8)**2 +                      &
                        288*lam(2)*(2*lam(3) + lam(8)) +                        &
                        108*lam(1)*                                             &
                         (28*lam(2) - lam(4) - 3*lam(5) - lam(6) -              &
                         3*lam(7) + 6*lam(8)) + 30*lam(4)*lam(9) -              &
                        150*lam(5)*lam(9) - 30*lam(6)*lam(9) +                  &
                        150*lam(7)*lam(9) - 100*lam(9)**2 +                     &
                        18*(5*lam(4) - 11*lam(5) - 5*lam(6) +                   &
                         11*lam(7) - 10*lam(9))*lam(10) - 228*lam(10)**2)       &
+ 12*m2*(lam(8)*(-3*lam(4) + 15*lam(5) + 3*lam(6) - 15*lam(7) +                 &
                         20*lam(9) + 18*lam(10)) +                              &
                        4*lam(2)*                                               &
                         (3*lam(4) - 6*lam(5) - 3*lam(6) + 6*lam(7) -           &
                          5*(lam(9) + 3*lam(10))))))))))/432.
  dtLAM(3)=(2*dRqrd*expdISr*quB*qxy2**2*(1 + Rqrd)**2*tanhdr**3*             &
       (-(expuISs*(1 + Rqsu)*(-1 + tanhus**2)*                                  &
            (144*lam(1)**2 - 336*lam(1)*lam(3) - 164*lam(3)**2 +                &
              9*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                        &
              84*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -                   &
              164*lam(8)**2)) +                                                 &
         expdISs*(1 + Rqsd)*(-1 + tanhds**2)*                                   &
          (144*lam(1)**2 + 528*lam(1)*lam(3) + 232*lam(3)**2 +                  &
            9*(lam(4)**2 + 13*lam(5)**2 + 10*lam(5)*lam(6) + lam(6)**2 +        &
               26*lam(5)*lam(7) + 10*lam(6)*lam(7) + 13*lam(7)**2 +             &
               2*lam(4)*(5*lam(5) + lam(6) + 5*lam(7))) +                       &
            12*(5*lam(4) + lam(5) + 5*lam(6) + lam(7))*lam(8) +                 &
            88*lam(8)**2)) - 2*dRqru*expuISr*qdB*qxy2**2*(1 + Rqru)**2*         &
       tanhur**3*(expdISs*(1 + Rqsd)*(-1 + tanhds**2)*                          &
          (144*lam(1)**2 - 336*lam(1)*lam(3) - 164*lam(3)**2 +                  &
            9*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                          &
            84*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) - 164*lam(8)**2       &
) - expuISs*(1 + Rqsu)*(-1 + tanhus**2)*                                        &
          (144*lam(1)**2 + 528*lam(1)*lam(3) + 232*lam(3)**2 +                  &
            9*(lam(4)**2 + 13*lam(5)**2 + 10*lam(5)*lam(6) + lam(6)**2 +        &
               26*lam(5)*lam(7) + 10*lam(6)*lam(7) + 13*lam(7)**2 +             &
               2*lam(4)*(5*lam(5) + lam(6) + 5*lam(7))) +                       &
            12*(5*lam(4) + lam(5) + 5*lam(6) + lam(7))*lam(8) +                 &
            88*lam(8)**2)) - expdISr*qxy2*tanhdr**2*                            &
       (dRqrd*quB*(expuISs*(1 + Rqsu)*                                          &
             (-2*qzt2*(1 + Rqrd)**2*tanhus +                                    &
               qdB*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*(-1 + tanhus**2))*            &
             (144*lam(1)**2 - 336*lam(1)*lam(3) - 164*lam(3)**2 +               &
               9*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                       &
               84*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -                  &
               164*lam(8)**2) -                                                 &
            expdISs*(1 + Rqsd)*                                                 &
             (-2*qzt2*(1 + Rqrd)**2*tanhds +                                    &
               qdB*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*(-1 + tanhds**2))*            &
             (144*lam(1)**2 + 528*lam(1)*lam(3) + 232*lam(3)**2 +               &
               9*(lam(4)**2 + 13*lam(5)**2 + 10*lam(5)*lam(6) +                 &
                  lam(6)**2 + 26*lam(5)*lam(7) + 10*lam(6)*lam(7) +             &
                  13*lam(7)**2 + 2*lam(4)*(5*lam(5) + lam(6) + 5*lam(7))        &
) + 12*(5*lam(4) + lam(5) + 5*lam(6) + lam(7))*lam(8) + 88*lam(8)**2)) +        &
         (1 + Rqrd)*(-(dRqsu*expuISs*qdB*(-1 + tanhus)*(1 + tanhus)*            &
               (quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +                            &
                 2*qxy2*(1 + Rqsu)**2*tanhus)*                                  &
               (32*(18*lam(1)**2 - 42*lam(1)*lam(3) + 11*lam(3)**2) +           &
                 9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                      &
                    2*lam(4)*(5*lam(5) - lam(6) - 5*lam(7)) +                   &
                    10*lam(6)*lam(7) + 13*lam(7)**2 -                           &
                    2*lam(5)*(5*lam(6) + 13*lam(7))) -                          &
                 144*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +           &
                 384*lam(9)**2 +                                                &
                 12*(5*lam(4) + 19*lam(5) - 5*lam(6) - 19*lam(7) -              &
                    32*lam(9))*lam(10) + 88*lam(10)**2)) +                      &
            dRqsd*expdISs*quB*                                                  &
             (qdB*(-1 + 2*qzt2*(1 + Rqsd)**2*ss) +                              &
               2*qxy2*(1 + Rqsd)**2*tanhds)*(-1 + tanhds**2)*                   &
             (64*(3*lam(1) + lam(3))**2 +                                       &
               (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))**2)))      &
- expuISr*qxy2*tanhur**2*(dRqru*qdB*                                            &
          (expdISs*(1 + Rqsd)*                                                  &
             (-2*qzt2*(1 + Rqru)**2*tanhds +                                    &
               quB*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*(-1 + tanhds**2))*            &
             (144*lam(1)**2 - 336*lam(1)*lam(3) - 164*lam(3)**2 +               &
               9*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                       &
               84*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -                  &
               164*lam(8)**2) -                                                 &
            expuISs*(1 + Rqsu)*                                                 &
             (-2*qzt2*(1 + Rqru)**2*tanhus +                                    &
               quB*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*(-1 + tanhus**2))*            &
             (144*lam(1)**2 + 528*lam(1)*lam(3) + 232*lam(3)**2 +               &
               9*(lam(4)**2 + 13*lam(5)**2 + 10*lam(5)*lam(6) +                 &
                  lam(6)**2 + 26*lam(5)*lam(7) + 10*lam(6)*lam(7) +             &
                  13*lam(7)**2 + 2*lam(4)*(5*lam(5) + lam(6) + 5*lam(7))        &
) + 12*(5*lam(4) + lam(5) + 5*lam(6) + lam(7))*lam(8) + 88*lam(8)**2)) -        &
         (1 + Rqru)*(dRqsd*expdISs*quB*                                         &
             (qdB*(-1 + 2*qzt2*(1 + Rqsd)**2*ss) +                              &
               2*qxy2*(1 + Rqsd)**2*tanhds)*(-1 + tanhds**2)*                   &
             (32*(18*lam(1)**2 - 42*lam(1)*lam(3) + 11*lam(3)**2) +             &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(5*lam(5) - lam(6) - 5*lam(7)) +                     &
                  10*lam(6)*lam(7) + 13*lam(7)**2 -                             &
                  2*lam(5)*(5*lam(6) + 13*lam(7))) -                            &
               144*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +             &
               384*lam(9)**2 +                                                  &
               12*(5*lam(4) + 19*lam(5) - 5*lam(6) - 19*lam(7) -                &
                  32*lam(9))*lam(10) + 88*lam(10)**2) -                         &
            dRqsu*expuISs*qdB*(-1 + tanhus)*(1 + tanhus)*                       &
             (quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +                              &
               2*qxy2*(1 + Rqsu)**2*tanhus)*                                    &
             (64*(3*lam(1) + lam(3))**2 +                                       &
               (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))**2)))      &
- expdISr*tanhdr*(-2*dRqrd*quB*qxy2*(1 + Rqrd)*                                 &
          (-(expuISs*(8*m2*(3*(-6*lam(1) + 7*lam(3))*                           &
                     (lam(4) - lam(5) + lam(6) - lam(7)) +                      &
                    2*(42*lam(1) + 41*lam(3))*lam(8)) +                         &
                 qxy2*(1 + Rqrd)*(1 + Rqsu)*                                    &
                  (144*lam(1)**2 - 336*lam(1)*lam(3) - 164*lam(3)**2 +          &
                    9*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                  &
                    84*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -             &
                    164*lam(8)**2) +                                            &
                 qzt2*(1 + Rqrd)*(1 + Rqsu)*                                    &
                  (144*lam(1)**2 - 336*lam(1)*lam(3) - 164*lam(3)**2 +          &
                    9*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                  &
                    84*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -             &
                    164*lam(8)**2))) +                                          &
            expdISs*(qxy2*(1 + Rqrd)*(1 + Rqsd)*                                &
                (144*lam(1)**2 + 528*lam(1)*lam(3) + 232*lam(3)**2 +            &
                  9*(lam(4)**2 + 13*lam(5)**2 + 10*lam(5)*lam(6) +              &
                     lam(6)**2 + 26*lam(5)*lam(7) + 10*lam(6)*lam(7) +          &
                     13*lam(7)**2 +                                             &
                     2*lam(4)*(5*lam(5) + lam(6) + 5*lam(7))) +                 &
                  12*(5*lam(4) + lam(5) + 5*lam(6) + lam(7))*lam(8) +           &
                  88*lam(8)**2) +                                               &
               qzt2*(1 + Rqrd)*(1 + Rqsd)*                                      &
                (8*(18*lam(1)**2 + 66*lam(1)*lam(3) + 29*lam(3)**2) +           &
                  9*(lam(4)**2 + 13*lam(5)**2 + 10*lam(5)*lam(6) +              &
                     lam(6)**2 + 26*lam(5)*lam(7) + 10*lam(6)*lam(7) +          &
                     13*lam(7)**2 +                                             &
                     2*lam(4)*(5*lam(5) + lam(6) + 5*lam(7))) +                 &
                  12*(5*lam(4) + lam(5) + 5*lam(6) + lam(7))*lam(8) +           &
                  88*lam(8)**2) -                                               &
               8*m2*(6*lam(1)*                                                  &
                   (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +                 &
                     10*lam(8)) +                                               &
                  lam(3)*(33*lam(4) + 57*lam(5) + 33*lam(6) +                   &
                     57*lam(7) + 56*lam(8))))) -                                &
         2*expdISs*quB*qxy2*(1 + Rqsd)*tanhds**2*                               &
          (-(dRqrd*qxy2*(1 + Rqrd)**2*                                          &
               (144*lam(1)**2 + 528*lam(1)*lam(3) + 232*lam(3)**2 +             &
                 9*(lam(4)**2 + 13*lam(5)**2 + 10*lam(5)*lam(6) +               &
                    lam(6)**2 + 26*lam(5)*lam(7) + 10*lam(6)*lam(7) +           &
                    13*lam(7)**2 +                                              &
                    2*lam(4)*(5*lam(5) + lam(6) + 5*lam(7))) +                  &
                 12*(5*lam(4) + lam(5) + 5*lam(6) + lam(7))*lam(8) +            &
                 88*lam(8)**2)) +                                               &
            16*dRqsd*m2*(3*lam(1) + lam(3))*                                    &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +            &
            dRqsd*qzt2*(1 + Rqrd)*(1 + Rqsd)*                                   &
             (64*(3*lam(1) + lam(3))**2 +                                       &
               (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))**2))       &
+ expdISs*qdB*quB*qzt2*tanhds*                                                  &
          (dRqrd*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                     &
             (144*lam(1)**2 + 528*lam(1)*lam(3) + 232*lam(3)**2 +               &
               9*(lam(4)**2 + 13*lam(5)**2 + 10*lam(5)*lam(6) +                 &
                  lam(6)**2 + 26*lam(5)*lam(7) + 10*lam(6)*lam(7) +             &
                  13*lam(7)**2 +                                                &
                  2*lam(4)*(5*lam(5) + lam(6) + 5*lam(7))) +                    &
               12*(5*lam(4) + lam(5) + 5*lam(6) + lam(7))*lam(8) +              &
               88*lam(8)**2) -                                                  &
            dRqsd*(576*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*               &
                lam(1)**2 + 64*(1 + Rqrd)*                                      &
                (-1 + 2*qzt2*(1 + Rqsd)**2*ss)*lam(3)**2 +                      &
               32*m2*(1 + Rqsd)*ss*lam(3)*                                      &
                (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)) +           &
               (1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                       &
                (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10))**2 +        &
               96*lam(1)*(4*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*          &
                   lam(3) + m2*(1 + Rqsd)*ss*                                   &
                   (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)))))       &
- expuISs*qdB*quB*qzt2*tanhus*                                                  &
          (dRqrd*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                     &
             (144*lam(1)**2 - 336*lam(1)*lam(3) - 164*lam(3)**2 +               &
               9*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                       &
               84*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -                  &
               164*lam(8)**2) -                                                 &
            dRqsu*(576*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*               &
                lam(1)**2 + 352*(1 + Rqrd)*                                     &
                (-1 + 2*qzt2*(1 + Rqsu)**2*ss)*lam(3)**2 -                      &
               16*m2*(1 + Rqsu)*ss*lam(3)*                                      &
                (21*lam(4) + 51*lam(5) - 21*lam(6) - 51*lam(7) -                &
                  96*lam(9) + 52*lam(10)) +                                     &
               (1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                       &
                (9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                      &
                     2*lam(4)*(5*lam(5) - lam(6) - 5*lam(7)) +                  &
                     10*lam(6)*lam(7) + 13*lam(7)**2 -                          &
                     2*lam(5)*(5*lam(6) + 13*lam(7))) -                         &
                  144*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +          &
                  384*lam(9)**2 +                                               &
                  12*(5*lam(4) + 19*lam(5) - 5*lam(6) - 19*lam(7) -             &
                     32*lam(9))*lam(10) + 88*lam(10)**2) +                      &
               96*lam(1)*(-14*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*        &
                   lam(3) + m2*(1 + Rqsu)*ss*                                   &
                   (3*lam(4) + 15*lam(5) -                                      &
                     3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10))))) +        &
         2*expuISs*qxy2*(1 + Rqsu)*tanhus**2*                                   &
          (-(dRqrd*quB*qxy2*(1 + Rqrd)**2*                                      &
               (144*lam(1)**2 - 336*lam(1)*lam(3) - 164*lam(3)**2 +             &
                 9*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                     &
                 84*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -                &
                 164*lam(8)**2)) +                                              &
            dRqsu*qdB*(qzt2*(1 + Rqrd)*(1 + Rqsu)*                              &
                (576*lam(1)**2 - 1344*lam(1)*lam(3) + 352*lam(3)**2 +           &
                  9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                     &
                     2*lam(4)*(5*lam(5) - lam(6) - 5*lam(7)) +                  &
                     10*lam(6)*lam(7) + 13*lam(7)**2 -                          &
                     2*lam(5)*(5*lam(6) + 13*lam(7))) -                         &
                  144*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +          &
                  384*lam(9)**2 +                                               &
                  12*(5*lam(4) + 19*lam(5) - 5*lam(6) - 19*lam(7) -             &
                     32*lam(9))*lam(10) + 88*lam(10)**2) +                      &
               8*m2*(lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) +               &
                     51*lam(7) + 96*lam(9) - 52*lam(10)) +                      &
                  6*lam(1)*(3*lam(4) + 15*lam(5) -                              &
                     3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10)))))) +       &
      expdISr*expdISs*quB*(2*dRqsd*qxy2**2*(1 + Rqrd)*(1 + Rqsd)**2*            &
          tanhds**3*(64*(3*lam(1) + lam(3))**2 +                                &
            (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))**2) -         &
         qdB*qxy2*tanhds**2*(dRqrd*(1 + Rqsd)*                                  &
             (-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                                    &
             (144*lam(1)**2 + 528*lam(1)*lam(3) + 232*lam(3)**2 +               &
               9*(lam(4)**2 + 13*lam(5)**2 + 10*lam(5)*lam(6) +                 &
                  lam(6)**2 + 26*lam(5)*lam(7) + 10*lam(6)*lam(7) +             &
                  13*lam(7)**2 +                                                &
                  2*lam(4)*(5*lam(5) + lam(6) + 5*lam(7))) +                    &
               12*(5*lam(4) + lam(5) + 5*lam(6) + lam(7))*lam(8) +              &
               88*lam(8)**2) -                                                  &
            dRqsd*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                    &
             (64*(3*lam(1) + lam(3))**2 +                                       &
               (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))**2))       &
- 2*dRqsd*qxy2*(1 + Rqsd)*tanhds*                                               &
          (16*m2*(3*lam(1) + lam(3))*                                           &
             (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)) +              &
            qxy2*(1 + Rqrd)*(1 + Rqsd)*                                         &
             (64*(3*lam(1) + lam(3))**2 +                                       &
               (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))**2)        &
+ qzt2*(1 + Rqrd)*(1 + Rqsd)*(64*(3*lam(1) + lam(3))**2 +                       &
               (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10))**2)) +       &
         qdB*(dRqrd*(qxy2*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*            &
                (144*lam(1)**2 + 528*lam(1)*lam(3) + 232*lam(3)**2 +            &
                  9*(lam(4)**2 + 13*lam(5)**2 + 10*lam(5)*lam(6) +              &
                     lam(6)**2 + 26*lam(5)*lam(7) + 10*lam(6)*lam(7) +          &
                     13*lam(7)**2 +                                             &
                     2*lam(4)*(5*lam(5) + lam(6) + 5*lam(7))) +                 &
                  12*(5*lam(4) + lam(5) + 5*lam(6) + lam(7))*lam(8) +           &
                  88*lam(8)**2) +                                               &
               qzt2*(144*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*             &
                   lam(1)**2 +                                                  &
                  232*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                &
                   lam(3)**2 -                                                  &
                  16*m2*(1 + Rqrd)*sr*lam(3)*                                   &
                   (33*lam(4) + 57*lam(5) + 33*lam(6) + 57*lam(7) +             &
                     56*lam(8)) +                                               &
                  (1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                    &
                   (9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                   &
                        10*lam(6)*lam(7) + 13*lam(7)**2 +                       &
                        2*lam(4)*(5*lam(5) + lam(6) + 5*lam(7)) +               &
                        2*lam(5)*(5*lam(6) + 13*lam(7))) +                      &
                     12*(5*lam(4) + lam(5) + 5*lam(6) + lam(7))*                &
                      lam(8) + 88*lam(8)**2) +                                  &
                  48*lam(1)*(11*(1 + Rqsd)*                                     &
                      (-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(3) -                   &
                     2*m2*(1 + Rqrd)*sr*                                        &
                      (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +              &
                        10*lam(8))))) -                                         &
            dRqsd*(qxy2*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*              &
                (64*(3*lam(1) + lam(3))**2 +                                    &
                  (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))**       &
                   2) + qzt2*(576*(1 + Rqrd)*                                   &
                   (-1 + 2*qzt2*(1 + Rqsd)**2*ss)*lam(1)**2 +                   &
                  64*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                 &
                   lam(3)**2 +                                                  &
                  32*m2*(1 + Rqsd)*ss*lam(3)*                                   &
                   (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)) +        &
                  (1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                    &
                   (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10))**2       &
+ 96*lam(1)*(4*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*lam(3) +               &
                     m2*(1 + Rqsd)*ss*                                          &
                      (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)))      &
)))) + expuISr*tanhur*(expdISs*qdB*quB*qzt2*tanhds*                             &
          (dRqru*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                     &
             (144*lam(1)**2 - 336*lam(1)*lam(3) - 164*lam(3)**2 +               &
               9*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                       &
               84*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -                  &
               164*lam(8)**2) -                                                 &
            dRqsd*(576*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*               &
                lam(1)**2 + 352*(1 + Rqru)*                                     &
                (-1 + 2*qzt2*(1 + Rqsd)**2*ss)*lam(3)**2 -                      &
               16*m2*(1 + Rqsd)*ss*lam(3)*                                      &
                (21*lam(4) + 51*lam(5) - 21*lam(6) - 51*lam(7) -                &
                  96*lam(9) + 52*lam(10)) +                                     &
               (1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                       &
                (9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                      &
                     2*lam(4)*(5*lam(5) - lam(6) - 5*lam(7)) +                  &
                     10*lam(6)*lam(7) + 13*lam(7)**2 -                          &
                     2*lam(5)*(5*lam(6) + 13*lam(7))) -                         &
                  144*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +          &
                  384*lam(9)**2 +                                               &
                  12*(5*lam(4) + 19*lam(5) - 5*lam(6) - 19*lam(7) -             &
                     32*lam(9))*lam(10) + 88*lam(10)**2) +                      &
               96*lam(1)*(-14*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*        &
                   lam(3) + m2*(1 + Rqsd)*ss*                                   &
                   (3*lam(4) + 15*lam(5) -                                      &
                     3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10))))) -        &
         2*expdISs*qxy2*(1 + Rqsd)*tanhds**2*                                   &
          (-(dRqru*qdB*qxy2*(1 + Rqru)**2*                                      &
               (144*lam(1)**2 - 336*lam(1)*lam(3) - 164*lam(3)**2 +             &
                 9*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                     &
                 84*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -                &
                 164*lam(8)**2)) +                                              &
            dRqsd*quB*(qzt2*(1 + Rqru)*(1 + Rqsd)*                              &
                (576*lam(1)**2 - 1344*lam(1)*lam(3) + 352*lam(3)**2 +           &
                  9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                     &
                     2*lam(4)*(5*lam(5) - lam(6) - 5*lam(7)) +                  &
                     10*lam(6)*lam(7) + 13*lam(7)**2 -                          &
                     2*lam(5)*(5*lam(6) + 13*lam(7))) -                         &
                  144*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +          &
                  384*lam(9)**2 +                                               &
                  12*(5*lam(4) + 19*lam(5) - 5*lam(6) - 19*lam(7) -             &
                     32*lam(9))*lam(10) + 88*lam(10)**2) +                      &
               8*m2*(lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) +               &
                     51*lam(7) + 96*lam(9) - 52*lam(10)) +                      &
                  6*lam(1)*(3*lam(4) + 15*lam(5) -                              &
                     3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10))))) +        &
         qdB*(-2*dRqru*qxy2*(1 + Rqru)*                                         &
             (expdISs*(8*m2*(3*(-6*lam(1) + 7*lam(3))*                          &
                      (lam(4) - lam(5) + lam(6) - lam(7)) +                     &
                     2*(42*lam(1) + 41*lam(3))*lam(8)) +                        &
                  qxy2*(1 + Rqru)*(1 + Rqsd)*                                   &
                   (144*lam(1)**2 - 336*lam(1)*lam(3) -                         &
                     164*lam(3)**2 +                                            &
                     9*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                 &
                     84*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -            &
                     164*lam(8)**2) +                                           &
                  qzt2*(1 + Rqru)*(1 + Rqsd)*                                   &
                   (144*lam(1)**2 - 336*lam(1)*lam(3) - 164*lam(3)**2 +         &
                     9*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                 &
                     84*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -            &
                     164*lam(8)**2)) -                                          &
               expuISs*(qxy2*(1 + Rqru)*(1 + Rqsu)*                             &
                   (144*lam(1)**2 + 528*lam(1)*lam(3) + 232*lam(3)**2 +         &
                     9*(lam(4)**2 + 13*lam(5)**2 + 10*lam(5)*lam(6) +           &
                        lam(6)**2 + 26*lam(5)*lam(7) +                          &
                        10*lam(6)*lam(7) + 13*lam(7)**2 +                       &
                        2*lam(4)*(5*lam(5) + lam(6) + 5*lam(7))) +              &
                     12*(5*lam(4) + lam(5) + 5*lam(6) + lam(7))*                &
                      lam(8) + 88*lam(8)**2) +                                  &
                  qzt2*(1 + Rqru)*(1 + Rqsu)*                                   &
                   (8*(18*lam(1)**2 + 66*lam(1)*lam(3) +                        &
                       29*lam(3)**2) +                                          &
                     9*(lam(4)**2 + 13*lam(5)**2 + 10*lam(5)*lam(6) +           &
                        lam(6)**2 + 26*lam(5)*lam(7) +                          &
                        10*lam(6)*lam(7) + 13*lam(7)**2 +                       &
                        2*lam(4)*(5*lam(5) + lam(6) + 5*lam(7))) +              &
                     12*(5*lam(4) + lam(5) + 5*lam(6) + lam(7))*                &
                      lam(8) + 88*lam(8)**2) -                                  &
                  8*m2*(6*lam(1)*                                               &
                      (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +              &
                        10*lam(8)) +                                            &
                     lam(3)*(33*lam(4) + 57*lam(5) + 33*lam(6) +                &
                        57*lam(7) + 56*lam(8))))) -                             &
            2*expuISs*qxy2*(1 + Rqsu)*tanhus**2*                                &
             (dRqru*qxy2*(1 + Rqru)**2*                                         &
                (144*lam(1)**2 + 528*lam(1)*lam(3) + 232*lam(3)**2 +            &
                  9*(lam(4)**2 + 13*lam(5)**2 + 10*lam(5)*lam(6) +              &
                     lam(6)**2 + 26*lam(5)*lam(7) + 10*lam(6)*lam(7) +          &
                     13*lam(7)**2 +                                             &
                     2*lam(4)*(5*lam(5) + lam(6) + 5*lam(7))) +                 &
                  12*(5*lam(4) + lam(5) + 5*lam(6) + lam(7))*lam(8) +           &
                  88*lam(8)**2) -                                               &
               dRqsu*(16*m2*(3*lam(1) + lam(3))*                                &
                   (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))        &
+ qzt2*(1 + Rqru)*(1 + Rqsu)*(64*(3*lam(1) + lam(3))**2 +                       &
                     (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) +                 &
                        4*lam(10))**2))) -                                      &
            expuISs*quB*qzt2*tanhus*                                            &
             (dRqru*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                  &
                (144*lam(1)**2 + 528*lam(1)*lam(3) + 232*lam(3)**2 +            &
                  9*(lam(4)**2 + 13*lam(5)**2 + 10*lam(5)*lam(6) +              &
                     lam(6)**2 + 26*lam(5)*lam(7) + 10*lam(6)*lam(7) +          &
                     13*lam(7)**2 +                                             &
                     2*lam(4)*(5*lam(5) + lam(6) + 5*lam(7))) +                 &
                  12*(5*lam(4) + lam(5) + 5*lam(6) + lam(7))*lam(8) +           &
                  88*lam(8)**2) -                                               &
               dRqsu*(576*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*            &
                   lam(1)**2 +                                                  &
                  64*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                 &
                   lam(3)**2 +                                                  &
                  32*m2*(1 + Rqsu)*ss*lam(3)*                                   &
                   (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)) +        &
                  (1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                    &
                   (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10))**2       &
+ 96*lam(1)*(4*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*lam(3) +               &
                     m2*(1 + Rqsu)*ss*                                          &
                      (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)))      &
)))) + expuISr*expuISs*qdB*(2*dRqsu*qxy2**2*(1 + Rqru)*(1 + Rqsu)**2*           &
          tanhus**3*(64*(3*lam(1) + lam(3))**2 +                                &
            (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))**2) -         &
         quB*qxy2*tanhus**2*(dRqru*(1 + Rqsu)*                                  &
             (-1 + 2*qzt2*(1 + Rqru)**2*sr)*                                    &
             (144*lam(1)**2 + 528*lam(1)*lam(3) + 232*lam(3)**2 +               &
               9*(lam(4)**2 + 13*lam(5)**2 + 10*lam(5)*lam(6) +                 &
                  lam(6)**2 + 26*lam(5)*lam(7) + 10*lam(6)*lam(7) +             &
                  13*lam(7)**2 +                                                &
                  2*lam(4)*(5*lam(5) + lam(6) + 5*lam(7))) +                    &
               12*(5*lam(4) + lam(5) + 5*lam(6) + lam(7))*lam(8) +              &
               88*lam(8)**2) -                                                  &
            dRqsu*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                    &
             (64*(3*lam(1) + lam(3))**2 +                                       &
               (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))**2))       &
- 2*dRqsu*qxy2*(1 + Rqsu)*tanhus*                                               &
          (16*m2*(3*lam(1) + lam(3))*                                           &
             (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)) +              &
            qxy2*(1 + Rqru)*(1 + Rqsu)*                                         &
             (64*(3*lam(1) + lam(3))**2 +                                       &
               (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))**2)        &
+ qzt2*(1 + Rqru)*(1 + Rqsu)*(64*(3*lam(1) + lam(3))**2 +                       &
               (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10))**2)) +       &
         quB*(dRqru*(qxy2*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*            &
                (144*lam(1)**2 + 528*lam(1)*lam(3) + 232*lam(3)**2 +            &
                  9*(lam(4)**2 + 13*lam(5)**2 + 10*lam(5)*lam(6) +              &
                     lam(6)**2 + 26*lam(5)*lam(7) + 10*lam(6)*lam(7) +          &
                     13*lam(7)**2 +                                             &
                     2*lam(4)*(5*lam(5) + lam(6) + 5*lam(7))) +                 &
                  12*(5*lam(4) + lam(5) + 5*lam(6) + lam(7))*lam(8) +           &
                  88*lam(8)**2) +                                               &
               qzt2*(144*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*             &
                   lam(1)**2 +                                                  &
                  232*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                &
                   lam(3)**2 -                                                  &
                  16*m2*(1 + Rqru)*sr*lam(3)*                                   &
                   (33*lam(4) + 57*lam(5) + 33*lam(6) + 57*lam(7) +             &
                     56*lam(8)) +                                               &
                  (1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                    &
                   (9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                   &
                        10*lam(6)*lam(7) + 13*lam(7)**2 +                       &
                        2*lam(4)*(5*lam(5) + lam(6) + 5*lam(7)) +               &
                        2*lam(5)*(5*lam(6) + 13*lam(7))) +                      &
                     12*(5*lam(4) + lam(5) + 5*lam(6) + lam(7))*                &
                      lam(8) + 88*lam(8)**2) +                                  &
                  48*lam(1)*(11*(1 + Rqsu)*                                     &
                      (-1 + 2*qzt2*(1 + Rqru)**2*sr)*lam(3) -                   &
                     2*m2*(1 + Rqru)*sr*                                        &
                      (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +              &
                        10*lam(8))))) -                                         &
            dRqsu*(qxy2*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*              &
                (64*(3*lam(1) + lam(3))**2 +                                    &
                  (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))**       &
                   2) + qzt2*(576*(1 + Rqru)*                                   &
                   (-1 + 2*qzt2*(1 + Rqsu)**2*ss)*lam(1)**2 +                   &
                  64*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                 &
                   lam(3)**2 +                                                  &
                  32*m2*(1 + Rqsu)*ss*lam(3)*                                   &
                   (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)) +        &
                  (1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                    &
                   (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10))**2       &
+ 96*lam(1)*(4*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*lam(3) +               &
                     m2*(1 + Rqsu)*ss*                                          &
                      (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)))      &
)))) - expdISs*expuISr*quB*(2*dRqsd*qxy2**2*(1 + Rqru)*(1 + Rqsd)**2*           &
          tanhds**3*(576*lam(1)**2 - 1344*lam(1)*lam(3) +                       &
            352*lam(3)**2 + 9*                                                  &
             (lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                            &
               2*lam(4)*(5*lam(5) - lam(6) - 5*lam(7)) +                        &
               10*lam(6)*lam(7) + 13*lam(7)**2 -                                &
               2*lam(5)*(5*lam(6) + 13*lam(7))) -                               &
            144*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +                &
            384*lam(9)**2 + 12*                                                 &
             (5*lam(4) + 19*lam(5) - 5*lam(6) - 19*lam(7) - 32*lam(9))*         &
             lam(10) + 88*lam(10)**2) +                                         &
         qdB*qxy2*tanhds**2*(-(dRqru*(1 + Rqsd)*                                &
               (-1 + 2*qzt2*(1 + Rqru)**2*sr)*                                  &
               (144*lam(1)**2 - 336*lam(1)*lam(3) - 164*lam(3)**2 +             &
                 9*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                     &
                 84*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -                &
                 164*lam(8)**2)) +                                              &
            dRqsd*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                    &
             (576*lam(1)**2 - 1344*lam(1)*lam(3) + 352*lam(3)**2 +              &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(5*lam(5) - lam(6) - 5*lam(7)) +                     &
                  10*lam(6)*lam(7) + 13*lam(7)**2 -                             &
                  2*lam(5)*(5*lam(6) + 13*lam(7))) -                            &
               144*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +             &
               384*lam(9)**2 +                                                  &
               12*(5*lam(4) + 19*lam(5) - 5*lam(6) - 19*lam(7) -                &
                  32*lam(9))*lam(10) + 88*lam(10)**2)) -                        &
         2*dRqsd*qxy2*(1 + Rqsd)*tanhds*                                        &
          (qxy2*(1 + Rqru)*(1 + Rqsd)*                                          &
             (576*lam(1)**2 - 1344*lam(1)*lam(3) + 352*lam(3)**2 +              &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(5*lam(5) - lam(6) - 5*lam(7)) +                     &
                  10*lam(6)*lam(7) + 13*lam(7)**2 -                             &
                  2*lam(5)*(5*lam(6) + 13*lam(7))) -                            &
               144*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +             &
               384*lam(9)**2 +                                                  &
               12*(5*lam(4) + 19*lam(5) - 5*lam(6) - 19*lam(7) -                &
                  32*lam(9))*lam(10) + 88*lam(10)**2) +                         &
            qzt2*(1 + Rqru)*(1 + Rqsd)*                                         &
             (576*lam(1)**2 - 1344*lam(1)*lam(3) + 352*lam(3)**2 +              &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(5*lam(5) - lam(6) - 5*lam(7)) +                     &
                  10*lam(6)*lam(7) + 13*lam(7)**2 -                             &
                  2*lam(5)*(5*lam(6) + 13*lam(7))) -                            &
               144*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +             &
               384*lam(9)**2 +                                                  &
               12*(5*lam(4) + 19*lam(5) - 5*lam(6) - 19*lam(7) -                &
                  32*lam(9))*lam(10) + 88*lam(10)**2) +                         &
            8*m2*(lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) +                  &
                  51*lam(7) + 96*lam(9) - 52*lam(10)) +                         &
               6*lam(1)*(3*lam(4) + 15*lam(5) -                                 &
                  3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10)))) +            &
         qdB*(dRqru*(qxy2*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*            &
                (144*lam(1)**2 - 336*lam(1)*lam(3) - 164*lam(3)**2 +            &
                  9*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                    &
                  84*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -               &
                  164*lam(8)**2) +                                              &
               qzt2*(144*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*             &
                   lam(1)**2 -                                                  &
                  164*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                &
                   lam(3)**2 +                                                  &
                  16*m2*(1 + Rqru)*sr*lam(3)*                                   &
                   (21*(lam(4) - lam(5) + lam(6) - lam(7)) + 82*lam(8))         &
+ (1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                                    &
                   (9*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                  &
                     84*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -            &
                     164*lam(8)**2) +                                           &
                  48*lam(1)*(-7*(1 + Rqsd)*                                     &
                      (-1 + 2*qzt2*(1 + Rqru)**2*sr)*lam(3) +                   &
                     2*m2*(1 + Rqru)*sr*                                        &
                      (-3*lam(4) + 3*(lam(5) - lam(6) + lam(7)) +               &
                        14*lam(8))))) -                                         &
            dRqsd*(qxy2*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*              &
                (576*lam(1)**2 - 1344*lam(1)*lam(3) + 352*lam(3)**2 +           &
                  9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                     &
                     2*lam(4)*(5*lam(5) - lam(6) - 5*lam(7)) +                  &
                     10*lam(6)*lam(7) + 13*lam(7)**2 -                          &
                     2*lam(5)*(5*lam(6) + 13*lam(7))) -                         &
                  144*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +          &
                  384*lam(9)**2 +                                               &
                  12*(5*lam(4) + 19*lam(5) - 5*lam(6) - 19*lam(7) -             &
                     32*lam(9))*lam(10) + 88*lam(10)**2) +                      &
               qzt2*(576*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*             &
                   lam(1)**2 +                                                  &
                  352*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                &
                   lam(3)**2 -                                                  &
                  16*m2*(1 + Rqsd)*ss*lam(3)*                                   &
                   (21*lam(4) + 51*lam(5) - 21*lam(6) - 51*lam(7) -             &
                     96*lam(9) + 52*lam(10)) +                                  &
                  (1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                    &
                   (9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                   &
                        2*lam(4)*(5*lam(5) - lam(6) - 5*lam(7)) +               &
                        10*lam(6)*lam(7) + 13*lam(7)**2 -                       &
                        2*lam(5)*(5*lam(6) + 13*lam(7))) -                      &
                     144*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*               &
                      lam(9) + 384*lam(9)**2 +                                  &
                     12*(5*lam(4) + 19*lam(5) - 5*lam(6) - 19*lam(7) -          &
                        32*lam(9))*lam(10) + 88*lam(10)**2) +                   &
                  96*lam(1)*(-14*(1 + Rqru)*                                    &
                      (-1 + 2*qzt2*(1 + Rqsd)**2*ss)*lam(3) +                   &
                     m2*(1 + Rqsd)*ss*                                          &
                      (3*lam(4) + 15*lam(5) -                                   &
                        3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10))))))      &
) - expdISr*expuISs*qdB*(2*dRqsu*qxy2**2*(1 + Rqrd)*(1 + Rqsu)**2*              &
          tanhus**3*(576*lam(1)**2 - 1344*lam(1)*lam(3) + 352*lam(3)**2 +       &
            9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                           &
               2*lam(4)*(5*lam(5) - lam(6) - 5*lam(7)) +                        &
               10*lam(6)*lam(7) + 13*lam(7)**2 -                                &
               2*lam(5)*(5*lam(6) + 13*lam(7))) -                               &
            144*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +                &
            384*lam(9)**2 + 12*                                                 &
             (5*lam(4) + 19*lam(5) - 5*lam(6) - 19*lam(7) - 32*lam(9))*         &
             lam(10) + 88*lam(10)**2) +                                         &
         quB*qxy2*tanhus**2*(-(dRqrd*(1 + Rqsu)*                                &
               (-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                                  &
               (144*lam(1)**2 - 336*lam(1)*lam(3) - 164*lam(3)**2 +             &
                 9*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                     &
                 84*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -                &
                 164*lam(8)**2)) +                                              &
            dRqsu*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                    &
             (576*lam(1)**2 - 1344*lam(1)*lam(3) + 352*lam(3)**2 +              &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(5*lam(5) - lam(6) - 5*lam(7)) +                     &
                  10*lam(6)*lam(7) + 13*lam(7)**2 -                             &
                  2*lam(5)*(5*lam(6) + 13*lam(7))) -                            &
               144*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +             &
               384*lam(9)**2 +                                                  &
               12*(5*lam(4) + 19*lam(5) - 5*lam(6) - 19*lam(7) -                &
                  32*lam(9))*lam(10) + 88*lam(10)**2)) -                        &
         2*dRqsu*qxy2*(1 + Rqsu)*tanhus*                                        &
          (qxy2*(1 + Rqrd)*(1 + Rqsu)*                                          &
             (576*lam(1)**2 - 1344*lam(1)*lam(3) + 352*lam(3)**2 +              &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(5*lam(5) - lam(6) - 5*lam(7)) +                     &
                  10*lam(6)*lam(7) + 13*lam(7)**2 -                             &
                  2*lam(5)*(5*lam(6) + 13*lam(7))) -                            &
               144*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +             &
               384*lam(9)**2 +                                                  &
               12*(5*lam(4) + 19*lam(5) - 5*lam(6) - 19*lam(7) -                &
                  32*lam(9))*lam(10) + 88*lam(10)**2) +                         &
            qzt2*(1 + Rqrd)*(1 + Rqsu)*                                         &
             (576*lam(1)**2 - 1344*lam(1)*lam(3) + 352*lam(3)**2 +              &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(5*lam(5) - lam(6) - 5*lam(7)) +                     &
                  10*lam(6)*lam(7) + 13*lam(7)**2 -                             &
                  2*lam(5)*(5*lam(6) + 13*lam(7))) -                            &
               144*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +             &
               384*lam(9)**2 +                                                  &
               12*(5*lam(4) + 19*lam(5) - 5*lam(6) - 19*lam(7) -                &
                  32*lam(9))*lam(10) + 88*lam(10)**2) +                         &
            8*m2*(lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) +                  &
                  51*lam(7) + 96*lam(9) - 52*lam(10)) +                         &
               6*lam(1)*(3*lam(4) + 15*lam(5) -                                 &
                  3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10)))) +            &
         quB*(dRqrd*(qxy2*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*            &
                (144*lam(1)**2 - 336*lam(1)*lam(3) - 164*lam(3)**2 +            &
                  9*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                    &
                  84*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -               &
                  164*lam(8)**2) +                                              &
               qzt2*(144*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*             &
                   lam(1)**2 -                                                  &
                  164*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                &
                   lam(3)**2 +                                                  &
                  16*m2*(1 + Rqrd)*sr*lam(3)*                                   &
                   (21*(lam(4) - lam(5) + lam(6) - lam(7)) + 82*lam(8)) +       &
                  (1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                    &
                   (9*(lam(4) - lam(5) + lam(6) - lam(7))**2 -                  &
                     84*(lam(4) - lam(5) + lam(6) - lam(7))*lam(8) -            &
                     164*lam(8)**2) +                                           &
                  48*lam(1)*(-7*(1 + Rqsu)*                                     &
                      (-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(3) +                   &
                     2*m2*(1 + Rqrd)*sr*                                        &
                      (-3*lam(4) + 3*(lam(5) - lam(6) + lam(7)) +               &
                        14*lam(8))))) -                                         &
            dRqsu*(qxy2*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*              &
                (576*lam(1)**2 - 1344*lam(1)*lam(3) + 352*lam(3)**2 +           &
                  9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                     &
                     2*lam(4)*(5*lam(5) - lam(6) - 5*lam(7)) +                  &
                     10*lam(6)*lam(7) + 13*lam(7)**2 -                          &
                     2*lam(5)*(5*lam(6) + 13*lam(7))) -                         &
                  144*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +          &
                  384*lam(9)**2 +                                               &
                  12*(5*lam(4) + 19*lam(5) - 5*lam(6) - 19*lam(7) -             &
                     32*lam(9))*lam(10) + 88*lam(10)**2) +                      &
               qzt2*(576*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*             &
                   lam(1)**2 +                                                  &
                  352*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                &
                   lam(3)**2 -                                                  &
                  16*m2*(1 + Rqsu)*ss*lam(3)*                                   &
                   (21*lam(4) + 51*lam(5) - 21*lam(6) - 51*lam(7) -             &
                     96*lam(9) + 52*lam(10)) +                                  &
                  (1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                    &
                   (9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                   &
                        2*lam(4)*(5*lam(5) - lam(6) - 5*lam(7)) +               &
                        10*lam(6)*lam(7) + 13*lam(7)**2 -                       &
                        2*lam(5)*(5*lam(6) + 13*lam(7))) -                      &
                     144*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +       &
                     384*lam(9)**2 +                                            &
                     12*(5*lam(4) + 19*lam(5) - 5*lam(6) - 19*lam(7) -          &
                        32*lam(9))*lam(10) + 88*lam(10)**2) +                   &
                  96*lam(1)*(-14*(1 + Rqrd)*                                    &
                      (-1 + 2*qzt2*(1 + Rqsu)**2*ss)*lam(3) +                   &
                     m2*(1 + Rqsu)*ss*                                          &
                      (3*lam(4) + 15*lam(5) -                                   &
                        3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10)))))))     &
)/216.
  dtLAM(4)=                                                             &
   (-(expdISs*expuISr*quB*(-(dRqru*qdB*(1 + Rqsd)*                              &
              (-1 + 2*qzt2*(1 + Rqru)**2*sr)*                                   &
              (qxy2 + qzt2 - qxy2*tanhds**2)*                                   &
              (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +                 &
                72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -            &
                   2*lam(8)) +                                                  &
                24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                           &
                   lam(8)*(17*lam(3) + 6*lam(8))) +                             &
                (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                    &
                   10*lam(9) - 12*lam(10))**2)) +                               &
           24*dRqru*m2*qdB*qzt2*(1 + Rqru)*sr*                                  &
            (216*lam(1)**2 - 24*lam(1)*lam(3) + 34*lam(3)**2 +                  &
              27*lam(5)**2 + 27*lam(6)**2 + 27*(lam(4) - lam(7))**2 -           &
              12*lam(6)*lam(8) + 34*lam(8)**2 +                                 &
              6*lam(5)*(-9*lam(6) + 2*lam(8)) + 20*lam(8)*lam(9) +              &
              48*lam(9)**2 + 4*lam(2)*                                          &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -         &
                 12*lam(10)) - 24*lam(8)*lam(10)) +                             &
           3*dRqsd*(1 + Rqru)*                                                  &
            (qdB*(-1 + 2*qzt2*(1 + Rqsd)**2*ss) +                               &
              2*qxy2*(1 + Rqsd)**2*tanhds)*                                     &
            (qxy2 + qzt2 - qxy2*tanhds**2)*                                     &
            (48*lam(2)**2 + 216*lam(2)*                                         &
               (lam(4) - lam(5) + lam(6) - lam(7)) +                            &
              3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                  &
                 lam(6)**2 + 2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +         &
                 26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +          &
              96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -               &
              64*lam(8)**2 - 128*lam(3)*                                        &
               (3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) +                     &
              48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -            &
                 16*lam(9) + 16*lam(10))) +                                     &
           72*dRqsd*m2*(1 + Rqsd)*(qdB*qzt2*ss + qxy2*tanhds)*                  &
            (144*lam(1)**2 + 72*lam(2)**2 +                                     &
              9*(lam(4)**2 + 5*lam(5)**2 - 2*lam(5)*lam(6) + lam(6)**2 -        &
                 2*lam(4)*lam(7) + 5*lam(7)**2) + 32*lam(8)**2 +                &
              4*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7) +                 &
                 16*lam(8)) +                                                   &
              16*(2*(lam(3)**2 + lam(9)**2) - 2*lam(9)*lam(10) +                &
                 lam(10)**2)))) -                                               &
      expdISr*expuISs*qdB*(-(dRqrd*quB*(1 + Rqsu)*                              &
            (-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                                     &
            (qxy2 + qzt2 - qxy2*tanhus**2)*                                     &
            (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +                   &
              72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -              &
                 2*lam(8)) + 24*                                                &
               (24*lam(2)**2 + 24*lam(2)*lam(8) +                               &
                 lam(8)*(17*lam(3) + 6*lam(8))) +                               &
              (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -          &
                 12*lam(10))**2)) +                                             &
         24*dRqrd*m2*quB*qzt2*(1 + Rqrd)*sr*                                    &
          (216*lam(1)**2 - 24*lam(1)*lam(3) + 34*lam(3)**2 +                    &
            27*lam(5)**2 + 27*lam(6)**2 + 27*(lam(4) - lam(7))**2 -             &
            12*lam(6)*lam(8) + 34*lam(8)**2 +                                   &
            6*lam(5)*(-9*lam(6) + 2*lam(8)) + 20*lam(8)*lam(9) +                &
            48*lam(9)**2 + 4*lam(2)*                                            &
             (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -           &
               12*lam(10)) - 24*lam(8)*lam(10)) +                               &
         3*dRqsu*(1 + Rqrd)*(quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +               &
            2*qxy2*(1 + Rqsu)**2*tanhus)*(qxy2 + qzt2 - qxy2*tanhus**2)*        &
          (48*lam(2)**2 + 216*lam(2)*                                           &
             (lam(4) - lam(5) + lam(6) - lam(7)) +                              &
            3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                    &
               lam(6)**2 + 2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +           &
               26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +            &
            96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -                 &
            64*lam(8)**2 - 128*lam(3)*                                          &
             (3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) +                       &
            48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -              &
               16*lam(9) + 16*lam(10))) +                                       &
         72*dRqsu*m2*(1 + Rqsu)*(quB*qzt2*ss + qxy2*tanhus)*                    &
          (144*lam(1)**2 + 72*lam(2)**2 +                                       &
            9*(lam(4)**2 + 5*lam(5)**2 - 2*lam(5)*lam(6) + lam(6)**2 -          &
               2*lam(4)*lam(7) + 5*lam(7)**2) + 32*lam(8)**2 +                  &
            4*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7) +                   &
               16*lam(8)) + 16*                                                 &
             (2*(lam(3)**2 + lam(9)**2) - 2*lam(9)*lam(10) + lam(10)**2)))      &
- 2*dRqrd*expdISr*quB*qxy2**2*(1 + Rqrd)**2*tanhdr**3*                          &
       (-(expuISs*(1 + Rqsu)*(-1 + tanhus**2)*                                  &
            (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +                   &
              72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -              &
                 2*lam(8)) + 24*                                                &
               (24*lam(2)**2 + 24*lam(2)*lam(8) +                               &
                 lam(8)*(17*lam(3) + 6*lam(8))) +                               &
              (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -          &
                 12*lam(10))**2)) +                                             &
         expdISs*(1 + Rqsd)*(-1 + tanhds**2)*                                   &
          (576*lam(2)**2 + 14265*lam(4)**2 + 3798*lam(4)*lam(5) +               &
            117*lam(5)**2 + 1278*lam(4)*lam(6) + 90*lam(5)*lam(6) +             &
            9*lam(6)**2 - 3798*lam(4)*lam(7) - 234*lam(5)*lam(7) -              &
            90*lam(6)*lam(7) + 117*lam(7)**2 +                                  &
            72*lam(1)*(-9*(lam(4) + lam(5) + lam(6) + lam(7)) -                 &
               2*lam(8)) - 288*lam(8)**2 -                                      &
            12*lam(3)*(21*lam(4) + 69*lam(5) + 21*lam(6) + 69*lam(7) +          &
               32*lam(8)) + 6684*lam(4)*lam(9) + 564*lam(5)*lam(9) +            &
            228*lam(6)*lam(9) - 564*lam(7)*lam(9) + 568*lam(9)**2 -             &
            12*(1089*lam(4) + 177*lam(5) + 63*lam(6) - 177*lam(7) +             &
               286*lam(9))*lam(10) + 2616*lam(10)**2 -                          &
            144*lam(2)*(69*lam(4) + 9*lam(5) + 3*lam(6) - 9*lam(7) -            &
               4*(lam(8) - 4*lam(9) + 8*lam(10))))) +                           &
      2*dRqru*expuISr*qdB*qxy2**2*(1 + Rqru)**2*tanhur**3*                      &
       (expdISs*(1 + Rqsd)*(-1 + tanhds**2)*                                    &
          (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +                     &
            72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -                &
               2*lam(8)) + 24*                                                  &
             (24*lam(2)**2 + 24*lam(2)*lam(8) +                                 &
               lam(8)*(17*lam(3) + 6*lam(8))) +                                 &
            (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -            &
               12*lam(10))**2) -                                                &
         expuISs*(1 + Rqsu)*(-1 + tanhus**2)*                                   &
          (576*lam(2)**2 + 14265*lam(4)**2 + 3798*lam(4)*lam(5) +               &
            117*lam(5)**2 + 1278*lam(4)*lam(6) + 90*lam(5)*lam(6) +             &
            9*lam(6)**2 - 3798*lam(4)*lam(7) - 234*lam(5)*lam(7) -              &
            90*lam(6)*lam(7) + 117*lam(7)**2 +                                  &
            72*lam(1)*(-9*(lam(4) + lam(5) + lam(6) + lam(7)) -                 &
               2*lam(8)) - 288*lam(8)**2 -                                      &
            12*lam(3)*(21*lam(4) + 69*lam(5) + 21*lam(6) + 69*lam(7) +          &
               32*lam(8)) + 6684*lam(4)*lam(9) + 564*lam(5)*lam(9) +            &
            228*lam(6)*lam(9) - 564*lam(7)*lam(9) + 568*lam(9)**2 -             &
            12*(1089*lam(4) + 177*lam(5) + 63*lam(6) - 177*lam(7) +             &
               286*lam(9))*lam(10) + 2616*lam(10)**2 -                          &
            144*lam(2)*(69*lam(4) + 9*lam(5) + 3*lam(6) - 9*lam(7) -            &
               4*(lam(8) - 4*lam(9) + 8*lam(10))))) +                           &
      expdISr*expdISs*quB*(3*dRqsd*                                             &
          ((1 + Rqrd)*(qdB*(-1 + 2*qzt2*(1 + Rqsd)**2*ss) +                     &
               2*qxy2*(1 + Rqsd)**2*tanhds)*                                    &
             (qxy2 + qzt2 - qxy2*tanhds**2)*                                    &
             (3*(16*lam(2)**2 -                                                 &
                  192*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) -              &
                  64*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) -               &
                  72*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +               &
                  (lam(4) + lam(5) + lam(6) + lam(7))**2) -                     &
               256*(3*lam(1) + lam(3))*lam(10)) -                               &
            8*m2*(1 + Rqsd)*(qdB*qzt2*ss + qxy2*tanhds)*                        &
             (64*(3*lam(1) + lam(3))**2 +                                       &
               3*(72*lam(2)**2 -                                                &
                  4*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                &
                  9*((lam(5) + lam(6))**2 + (lam(4) + lam(7))**2) +             &
                  16*lam(10)**2))) +                                            &
         dRqrd*qdB*(-24*m2*qzt2*(1 + Rqrd)*sr*                                  &
             (-144*lam(2)**2 +                                                  &
               4*(54*lam(1)**2 + 42*lam(1)*lam(3) + 23*lam(3)**2) -             &
               1161*lam(4)**2 + 135*lam(5)**2 + 54*lam(5)*lam(6) +              &
               27*lam(6)**2 + 135*lam(7)**2 + 60*lam(5)*lam(8) +                &
               12*lam(6)*lam(8) + 20*lam(8)**2 + 40*lam(8)*lam(9) +             &
               48*lam(9)**2 +                                                   &
               4*lam(2)*(213*lam(4) + 15*lam(5) + 3*lam(6) -                    &
                  15*lam(7) + 38*lam(9) - 126*lam(10)) -                        &
               18*lam(4)*(18*lam(5) + 6*lam(6) - 21*lam(7) +                    &
                  32*lam(9) - 64*lam(10)) +                                     &
               12*(3*lam(8) - 8*lam(9))*lam(10) + 96*lam(10)**2) -              &
            (1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                          &
             (qxy2 + qzt2 - qxy2*tanhds**2)*                                    &
             (576*lam(2)**2 + 14265*lam(4)**2 + 3798*lam(4)*lam(5) +            &
               117*lam(5)**2 + 1278*lam(4)*lam(6) + 90*lam(5)*lam(6) +          &
               9*lam(6)**2 - 3798*lam(4)*lam(7) - 234*lam(5)*lam(7) -           &
               90*lam(6)*lam(7) + 117*lam(7)**2 +                               &
               72*lam(1)*(-9*(lam(4) + lam(5) + lam(6) + lam(7)) -              &
                  2*lam(8)) - 288*lam(8)**2 -                                   &
               12*lam(3)*(21*lam(4) + 69*lam(5) + 21*lam(6) +                   &
                  69*lam(7) + 32*lam(8)) + 6684*lam(4)*lam(9) +                 &
               564*lam(5)*lam(9) + 228*lam(6)*lam(9) -                          &
               564*lam(7)*lam(9) + 568*lam(9)**2 -                              &
               12*(1089*lam(4) + 177*lam(5) + 63*lam(6) - 177*lam(7) +          &
                  286*lam(9))*lam(10) + 2616*lam(10)**2 -                       &
               144*lam(2)*(69*lam(4) + 9*lam(5) + 3*lam(6) - 9*lam(7) -         &
                  4*(lam(8) - 4*lam(9) + 8*lam(10)))))) +                       &
      expuISr*expuISs*qdB*(3*dRqsu*                                             &
          ((1 + Rqru)*(quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +                     &
               2*qxy2*(1 + Rqsu)**2*tanhus)*                                    &
             (qxy2 + qzt2 - qxy2*tanhus**2)*                                    &
             (3*(16*lam(2)**2 -                                                 &
                  192*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) -              &
                  64*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) -               &
                  72*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +               &
                  (lam(4) + lam(5) + lam(6) + lam(7))**2) -                     &
               256*(3*lam(1) + lam(3))*lam(10)) -                               &
            8*m2*(1 + Rqsu)*(quB*qzt2*ss + qxy2*tanhus)*                        &
             (64*(3*lam(1) + lam(3))**2 +                                       &
               3*(72*lam(2)**2 -                                                &
                  4*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                &
                  9*((lam(5) + lam(6))**2 + (lam(4) + lam(7))**2) +             &
                  16*lam(10)**2))) +                                            &
         dRqru*quB*(-24*m2*qzt2*(1 + Rqru)*sr*                                  &
             (-144*lam(2)**2 +                                                  &
               4*(54*lam(1)**2 + 42*lam(1)*lam(3) + 23*lam(3)**2) -             &
               1161*lam(4)**2 + 135*lam(5)**2 + 54*lam(5)*lam(6) +              &
               27*lam(6)**2 + 135*lam(7)**2 + 60*lam(5)*lam(8) +                &
               12*lam(6)*lam(8) + 20*lam(8)**2 + 40*lam(8)*lam(9) +             &
               48*lam(9)**2 +                                                   &
               4*lam(2)*(213*lam(4) + 15*lam(5) + 3*lam(6) -                    &
                  15*lam(7) + 38*lam(9) - 126*lam(10)) -                        &
               18*lam(4)*(18*lam(5) + 6*lam(6) - 21*lam(7) +                    &
                  32*lam(9) - 64*lam(10)) +                                     &
               12*(3*lam(8) - 8*lam(9))*lam(10) + 96*lam(10)**2) -              &
            (1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                          &
             (qxy2 + qzt2 - qxy2*tanhus**2)*                                    &
             (576*lam(2)**2 + 14265*lam(4)**2 + 3798*lam(4)*lam(5) +            &
               117*lam(5)**2 + 1278*lam(4)*lam(6) + 90*lam(5)*lam(6) +          &
               9*lam(6)**2 - 3798*lam(4)*lam(7) - 234*lam(5)*lam(7) -           &
               90*lam(6)*lam(7) + 117*lam(7)**2 +                               &
               72*lam(1)*(-9*(lam(4) + lam(5) + lam(6) + lam(7)) -              &
                  2*lam(8)) - 288*lam(8)**2 -                                   &
               12*lam(3)*(21*lam(4) + 69*lam(5) + 21*lam(6) +                   &
                  69*lam(7) + 32*lam(8)) + 6684*lam(4)*lam(9) +                 &
               564*lam(5)*lam(9) + 228*lam(6)*lam(9) -                          &
               564*lam(7)*lam(9) + 568*lam(9)**2 -                              &
               12*(1089*lam(4) + 177*lam(5) + 63*lam(6) - 177*lam(7) +          &
                  286*lam(9))*lam(10) + 2616*lam(10)**2 -                       &
               144*lam(2)*(69*lam(4) + 9*lam(5) + 3*lam(6) - 9*lam(7) -         &
                  4*(lam(8) - 4*lam(9) + 8*lam(10)))))) -                       &
      expdISr*qxy2*tanhdr**2*(-6*dRqsd*expdISs*quB*qxy2*(1 + Rqrd)*             &
          (1 + Rqsd)**2*tanhds**3*                                              &
          (48*lam(2)**2 - 216*lam(2)*                                           &
             (lam(4) + lam(5) + lam(6) + lam(7)) +                              &
            3*(lam(4) + lam(5) + lam(6) + lam(7))**2 -                          &
            192*lam(1)*(3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) +               &
               4*lam(10)) - 64*lam(3)*                                          &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))) +           &
         6*dRqsu*expuISs*qdB*qxy2*(1 + Rqrd)*(1 + Rqsu)**2*tanhus**3*           &
          (48*lam(2)**2 + 216*lam(2)*                                           &
             (lam(4) - lam(5) + lam(6) - lam(7)) +                              &
            3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                    &
               lam(6)**2 + 2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +           &
               26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +            &
            96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -                 &
            64*lam(8)**2 - 128*lam(3)*                                          &
             (3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) +                       &
            48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -              &
               16*lam(9) + 16*lam(10))) -                                       &
         expuISs*qdB*quB*tanhus**2*                                             &
          (dRqrd*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                     &
             (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +                  &
               72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -             &
                  2*lam(8)) +                                                   &
               24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                            &
                  lam(8)*(17*lam(3) + 6*lam(8))) +                              &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                     &
                  10*lam(9) - 12*lam(10))**2) -                                 &
            3*dRqsu*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                  &
             (48*lam(2)**2 + 216*lam(2)*                                        &
                (lam(4) - lam(5) + lam(6) - lam(7)) +                           &
               3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                 &
                  lam(6)**2 +                                                   &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                    &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -              &
               64*lam(8)**2 -                                                   &
               128*lam(3)*(3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) +          &
               48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -           &
                  16*lam(9) + 16*lam(10)))) -                                   &
         2*expuISs*(1 + Rqrd)*tanhus*                                           &
          (3*dRqsu*qdB*qxy2*(1 + Rqsu)**2*                                      &
             (48*lam(2)**2 +                                                    &
               216*lam(2)*(lam(4) - lam(5) + lam(6) - lam(7)) +                 &
               3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                 &
                  lam(6)**2 +                                                   &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                    &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -              &
               64*lam(8)**2 -                                                   &
               128*lam(3)*(3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) +          &
               48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -           &
                  16*lam(9) + 16*lam(10))) +                                    &
            dRqrd*quB*(qzt2*(1 + Rqrd)*(1 + Rqsu)*                              &
                (36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -                &
                  72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -          &
                     2*lam(8)) +                                                &
                  24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                         &
                     lam(8)*(-17*lam(3) + 6*lam(8))) +                          &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) -                              &
               6*m2*(-9*lam(4)**2 - 9*lam(5)**2 + 18*lam(5)*lam(6) -            &
                  9*lam(6)**2 + 18*lam(5)*lam(7) - 18*lam(6)*lam(7) -           &
                  9*lam(7)**2 + 12*lam(5)*lam(8) - 12*lam(6)*lam(8) -           &
                  12*lam(7)*lam(8) +                                            &
                  6*lam(4)*(-3*lam(5) + 3*(lam(6) + lam(7)) +                   &
                     2*lam(8)) + 40*lam(8)*lam(9) - 32*lam(9)**2 +              &
                  8*lam(2)*(3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +         &
                     10*lam(9) - 12*lam(10)) - 48*lam(8)*lam(10)))) +           &
         expdISs*qdB*quB*tanhds**2*                                             &
          (-3*dRqsd*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                  &
             (48*lam(2)**2 -                                                    &
               216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                 &
               3*(lam(4) + lam(5) + lam(6) + lam(7))**2 -                       &
               192*lam(1)*(3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) +            &
                  4*lam(10)) -                                                  &
               64*lam(3)*(3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) +             &
                  4*lam(10))) +                                                 &
            dRqrd*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                    &
             (576*lam(2)**2 + 14265*lam(4)**2 + 3798*lam(4)*lam(5) +            &
               117*lam(5)**2 + 1278*lam(4)*lam(6) + 90*lam(5)*lam(6) +          &
               9*lam(6)**2 - 3798*lam(4)*lam(7) - 234*lam(5)*lam(7) -           &
               90*lam(6)*lam(7) + 117*lam(7)**2 +                               &
               72*lam(1)*(-9*(lam(4) + lam(5) + lam(6) + lam(7)) -              &
                  2*lam(8)) - 288*lam(8)**2 -                                   &
               12*lam(3)*(21*lam(4) + 69*lam(5) + 21*lam(6) +                   &
                  69*lam(7) + 32*lam(8)) + 6684*lam(4)*lam(9) +                 &
               564*lam(5)*lam(9) + 228*lam(6)*lam(9) -                          &
               564*lam(7)*lam(9) + 568*lam(9)**2 -                              &
               12*(1089*lam(4) + 177*lam(5) + 63*lam(6) - 177*lam(7) +          &
                  286*lam(9))*lam(10) + 2616*lam(10)**2 -                       &
               144*lam(2)*(69*lam(4) + 9*lam(5) + 3*lam(6) - 9*lam(7) -         &
                  4*(lam(8) - 4*lam(9) + 8*lam(10))))) +                        &
         qdB*quB*(3*(1 + Rqrd)*                                                 &
             (dRqsd*expdISs*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                     &
                (48*lam(2)**2 -                                                 &
                  216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +              &
                  3*(lam(4) + lam(5) + lam(6) + lam(7))**2 -                    &
                  192*lam(1)*                                                   &
                   (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) +                   &
                     4*lam(10)) -                                               &
                  64*lam(3)*(3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) +          &
                     4*lam(10))) -                                              &
               dRqsu*expuISs*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                    &
                (48*lam(2)**2 +                                                 &
                  216*lam(2)*(lam(4) - lam(5) + lam(6) - lam(7)) +              &
                  3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +              &
                     lam(6)**2 +                                                &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                 &
                     26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2         &
) + 96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) - 64*lam(8)**2 -          &
                  128*lam(3)*                                                   &
                   (3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) +                 &
                  48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) -                   &
                     3*lam(7) - 16*lam(9) + 16*lam(10)))) -                     &
            dRqrd*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                               &
             (-(expuISs*(1 + Rqsu)*                                             &
                  (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +             &
                    72*lam(1)*                                                  &
                     (9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -                 &
                       2*lam(8)) +                                              &
                    24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                       &
                       lam(8)*(17*lam(3) + 6*lam(8))) +                         &
                    (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                &
                       10*lam(9) - 12*lam(10))**2)) +                           &
               expdISs*(1 + Rqsd)*                                              &
                (576*lam(2)**2 + 14265*lam(4)**2 + 3798*lam(4)*lam(5) +         &
                  117*lam(5)**2 + 1278*lam(4)*lam(6) +                          &
                  90*lam(5)*lam(6) + 9*lam(6)**2 - 3798*lam(4)*lam(7) -         &
                  234*lam(5)*lam(7) - 90*lam(6)*lam(7) + 117*lam(7)**2 +        &
                  72*lam(1)*(-9*(lam(4) + lam(5) + lam(6) + lam(7)) -           &
                     2*lam(8)) - 288*lam(8)**2 -                                &
                  12*lam(3)*(21*lam(4) + 69*lam(5) + 21*lam(6) +                &
                     69*lam(7) + 32*lam(8)) + 6684*lam(4)*lam(9) +              &
                  564*lam(5)*lam(9) + 228*lam(6)*lam(9) -                       &
                  564*lam(7)*lam(9) + 568*lam(9)**2 -                           &
                  12*(1089*lam(4) + 177*lam(5) + 63*lam(6) -                    &
                     177*lam(7) + 286*lam(9))*lam(10) +                         &
                  2616*lam(10)**2 -                                             &
                  144*lam(2)*(69*lam(4) + 9*lam(5) + 3*lam(6) -                 &
                     9*lam(7) - 4*(lam(8) - 4*lam(9) + 8*lam(10)))))) +         &
         2*expdISs*quB*(1 + Rqrd)*tanhds*                                       &
          (3*dRqsd*qxy2*(1 + Rqsd)**2*                                          &
             (48*lam(2)**2 - 216*lam(2)*                                        &
                (lam(4) + lam(5) + lam(6) + lam(7)) +                           &
               3*(lam(4) + lam(5) + lam(6) + lam(7))**2 -                       &
               192*lam(1)*(3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) +            &
                  4*lam(10)) -                                                  &
               64*lam(3)*(3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) +             &
                  4*lam(10))) +                                                 &
            dRqrd*(-2*m2*(864*lam(2)**2 + 7173*lam(4)**2 +                      &
                  297*lam(5)**2 + 162*lam(5)*lam(6) + 45*lam(6)**2 -            &
                  594*lam(5)*lam(7) - 162*lam(6)*lam(7) +                       &
                  297*lam(7)**2 - 180*lam(5)*lam(8) - 36*lam(6)*lam(8) +        &
                  180*lam(7)*lam(8) + 576*lam(5)*lam(9) +                       &
                  192*lam(6)*lam(9) - 576*lam(7)*lam(9) -                       &
                  240*lam(8)*lam(9) + 608*lam(9)**2 +                           &
                  6*lam(4)*(297*lam(5) + 93*lam(6) - 297*lam(7) +               &
                     6*lam(8) + 544*(lam(9) - 2*lam(10))) -                     &
                  24*lam(2)*(213*lam(4) + 15*lam(5) + 3*lam(6) -                &
                     15*lam(7) + 38*lam(9) - 126*lam(10)) -                     &
                  8*(144*lam(5) + 48*lam(6) - 144*lam(7) + 27*lam(8) +          &
                     280*lam(9))*lam(10) + 2240*lam(10)**2) +                   &
               qzt2*(1 + Rqrd)*(1 + Rqsd)*                                      &
                (576*lam(2)**2 + 14265*lam(4)**2 + 3798*lam(4)*lam(5) +         &
                  117*lam(5)**2 + 1278*lam(4)*lam(6) + 90*lam(5)*lam(6) +       &
                  9*lam(6)**2 - 3798*lam(4)*lam(7) - 234*lam(5)*lam(7) -        &
                  90*lam(6)*lam(7) + 117*lam(7)**2 - 288*lam(8)**2 +            &
                  72*lam(1)*(9*(lam(4) + lam(5) + lam(6) + lam(7)) +            &
                     2*lam(8)) +                                                &
                  12*lam(3)*(21*lam(4) + 69*lam(5) + 21*lam(6) +                &
                     69*lam(7) + 32*lam(8)) + 6684*lam(4)*lam(9) +              &
                  564*lam(5)*lam(9) + 228*lam(6)*lam(9) -                       &
                  564*lam(7)*lam(9) + 568*lam(9)**2 -                           &
                  12*(1089*lam(4) + 177*lam(5) + 63*lam(6) -                    &
                     177*lam(7) + 286*lam(9))*lam(10) + 2616*lam(10)**2 -       &
                  144*lam(2)*(69*lam(4) + 9*lam(5) + 3*lam(6) -                 &
                     9*lam(7) - 4*(lam(8) - 4*lam(9) + 8*lam(10))))))) +        &
      expdISr*tanhdr*(2*expdISs*quB*qxy2*(1 + Rqsd)*tanhds**2*                  &
          (dRqsd*(6*m2*(128*(3*lam(1) + lam(3))**2 +                            &
                  15*(lam(4) - lam(5) - lam(6) + lam(7))**2) +                  &
               384*m2*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +             &
               64*m2*lam(10)**2 +                                               &
               3*qzt2*(1 + Rqrd)*(1 + Rqsd)*                                    &
                (48*lam(2)**2 -                                                 &
                  216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +              &
                  3*(192*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) +           &
                     64*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +            &
                     (lam(4) + lam(5) + lam(6) + lam(7))**2) +                  &
                  256*(3*lam(1) + lam(3))*lam(10))) +                           &
            dRqrd*qxy2*(1 + Rqrd)**2*                                           &
             (576*lam(2)**2 + 14265*lam(4)**2 + 3798*lam(4)*lam(5) +            &
               117*lam(5)**2 + 1278*lam(4)*lam(6) + 90*lam(5)*lam(6) +          &
               9*lam(6)**2 - 3798*lam(4)*lam(7) - 234*lam(5)*lam(7) -           &
               90*lam(6)*lam(7) + 117*lam(7)**2 +                               &
               72*lam(1)*(-9*(lam(4) + lam(5) + lam(6) + lam(7)) -              &
                  2*lam(8)) - 288*lam(8)**2 -                                   &
               12*lam(3)*(21*lam(4) + 69*lam(5) + 21*lam(6) +                   &
                  69*lam(7) + 32*lam(8)) + 6684*lam(4)*lam(9) +                 &
               564*lam(5)*lam(9) + 228*lam(6)*lam(9) -                          &
               564*lam(7)*lam(9) + 568*lam(9)**2 -                              &
               12*(1089*lam(4) + 177*lam(5) + 63*lam(6) - 177*lam(7) +          &
                  286*lam(9))*lam(10) + 2616*lam(10)**2 -                       &
               144*lam(2)*(69*lam(4) + 9*lam(5) + 3*lam(6) - 9*lam(7) -         &
                  4*(lam(8) - 4*lam(9) + 8*lam(10))))) +                        &
         2*dRqrd*quB*qxy2*(1 + Rqrd)*                                           &
          (-12*expuISs*m2*(216*lam(1)**2 - 24*lam(1)*lam(3) +                   &
               34*lam(3)**2 + 27*lam(5)**2 + 27*lam(6)**2 +                     &
               27*(lam(4) - lam(7))**2 - 12*lam(6)*lam(8) +                     &
               34*lam(8)**2 + 6*lam(5)*(-9*lam(6) + 2*lam(8)) +                 &
               20*lam(8)*lam(9) + 48*lam(9)**2 +                                &
               4*lam(2)*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +            &
                  10*lam(9))) +                                                 &
            expuISs*(qxy2 + qzt2)*(1 + Rqrd)*(1 + Rqsu)*                        &
             (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +                  &
               72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -             &
                  2*lam(8)) +                                                   &
               24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                            &
                  lam(8)*(17*lam(3) + 6*lam(8))) +                              &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                     &
                  10*lam(9) - 12*lam(10))**2) +                                 &
            288*expuISs*m2*(2*lam(2) + lam(8))*lam(10) -                        &
            12*expdISs*m2*(-144*lam(2)**2 +                                     &
               4*(54*lam(1)**2 + 42*lam(1)*lam(3) + 23*lam(3)**2) -             &
               1161*lam(4)**2 + 135*lam(5)**2 + 54*lam(5)*lam(6) +              &
               27*lam(6)**2 + 135*lam(7)**2 + 60*lam(5)*lam(8) +                &
               12*lam(6)*lam(8) + 20*lam(8)**2 + 40*lam(8)*lam(9) +             &
               48*lam(9)**2 +                                                   &
               4*lam(2)*(213*lam(4) + 15*lam(5) + 3*lam(6) -                    &
                  15*lam(7) + 38*lam(9) - 126*lam(10)) -                        &
               18*lam(4)*(18*lam(5) + 6*lam(6) - 21*lam(7) +                    &
                  32*lam(9) - 64*lam(10)) +                                     &
               12*(3*lam(8) - 8*lam(9))*lam(10) + 96*lam(10)**2) -              &
            expdISs*(qxy2 + qzt2)*(1 + Rqrd)*(1 + Rqsd)*                        &
             (576*lam(2)**2 + 14265*lam(4)**2 + 3798*lam(4)*lam(5) +            &
               117*lam(5)**2 + 1278*lam(4)*lam(6) + 90*lam(5)*lam(6) +          &
               9*lam(6)**2 - 3798*lam(4)*lam(7) - 234*lam(5)*lam(7) -           &
               90*lam(6)*lam(7) + 117*lam(7)**2 +                               &
               72*lam(1)*(-9*(lam(4) + lam(5) + lam(6) + lam(7)) -              &
                  2*lam(8)) - 288*lam(8)**2 -                                   &
               12*lam(3)*(21*lam(4) + 69*lam(5) + 21*lam(6) +                   &
                  69*lam(7) + 32*lam(8)) + 6684*lam(4)*lam(9) +                 &
               564*lam(5)*lam(9) + 228*lam(6)*lam(9) -                          &
               564*lam(7)*lam(9) + 568*lam(9)**2 -                              &
               12*(1089*lam(4) + 177*lam(5) + 63*lam(6) - 177*lam(7) +          &
                  286*lam(9))*lam(10) + 2616*lam(10)**2 -                       &
               144*lam(2)*(69*lam(4) + 9*lam(5) + 3*lam(6) - 9*lam(7) -         &
                  4*(lam(8) - 4*lam(9) + 8*lam(10))))) -                        &
         2*expuISs*qxy2*(1 + Rqsu)*tanhus**2*                                   &
          (dRqrd*quB*qxy2*(1 + Rqrd)**2*                                        &
             (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +                  &
               72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -             &
                  2*lam(8)) +                                                   &
               24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                            &
                  lam(8)*(17*lam(3) + 6*lam(8))) +                              &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                     &
                  10*lam(9) - 12*lam(10))**2) +                                 &
            dRqsu*qdB*(3*qzt2*(1 + Rqrd)*(1 + Rqsu)*                            &
                (48*lam(2)**2 +                                                 &
                  216*lam(2)*(lam(4) - lam(5) + lam(6) - lam(7)) +              &
                  3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +              &
                     lam(6)**2 +                                                &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                 &
                     26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2         &
) + 96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) - 64*lam(8)**2 +          &
                  128*lam(3)*                                                   &
                   (3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) -                 &
                  48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) -                   &
                     3*lam(7) - 16*lam(9) + 16*lam(10))) +                      &
               2*m2*(-2592*lam(1)**2 - 576*lam(3)**2 - 27*lam(4)**2 -           &
                  63*lam(5)**2 - 27*lam(6)**2 + 18*lam(6)*lam(7) -              &
                  63*lam(7)**2 - 192*lam(6)*lam(9) - 576*lam(7)*lam(9) -        &
                  320*lam(9)**2 +                                               &
                  6*lam(4)*(3*lam(5) + 9*lam(6) - 3*lam(7) +                    &
                     32*lam(9) - 32*lam(10)) +                                  &
                  64*(3*(lam(6) + lam(7)) + 5*lam(9))*lam(10) -                 &
                  96*lam(10)**2 -                                               &
                  6*lam(5)*(3*lam(6) - 21*lam(7) +                              &
                     32*(-3*lam(9) + lam(10)))))) +                             &
         expuISs*qdB*quB*qzt2*tanhus*                                           &
          (-576*dRqrd*lam(2)**2 +                                               &
            dRqrd*(648*lam(1)*lam(4) - 288*m2*sr*lam(2)*lam(4) -                &
               288*m2*Rqrd*sr*lam(2)*lam(4) - 36*lam(3)*lam(4) -                &
               9*lam(4)**2 + 108*m2*sr*lam(4)**2 +                              &
               108*m2*Rqrd*sr*lam(4)**2 - 288*m2*sr*lam(2)*lam(5) -             &
               288*m2*Rqrd*sr*lam(2)*lam(5) + 36*lam(3)*lam(5) -                &
               18*lam(4)*lam(5) + 216*m2*sr*lam(4)*lam(5) +                     &
               216*m2*Rqrd*sr*lam(4)*lam(5) - 9*lam(5)**2 +                     &
               108*m2*sr*lam(5)**2 + 108*m2*Rqrd*sr*lam(5)**2 +                 &
               288*m2*sr*lam(2)*lam(6) + 288*m2*Rqrd*sr*lam(2)*lam(6) -         &
               36*lam(3)*lam(6) + 18*lam(4)*lam(6) -                            &
               216*m2*sr*lam(4)*lam(6) - 216*m2*Rqrd*sr*lam(4)*lam(6) +         &
               18*lam(5)*lam(6) - 216*m2*sr*lam(5)*lam(6) -                     &
               216*m2*Rqrd*sr*lam(5)*lam(6) - 9*lam(6)**2 +                     &
               108*m2*sr*lam(6)**2 + 108*m2*Rqrd*sr*lam(6)**2 +                 &
               288*m2*sr*lam(2)*lam(7) + 288*m2*Rqrd*sr*lam(2)*lam(7) +         &
               36*lam(3)*lam(7) + 18*lam(4)*lam(7) -                            &
               216*m2*sr*lam(4)*lam(7) - 216*m2*Rqrd*sr*lam(4)*lam(7) +         &
               18*lam(5)*lam(7) - 216*m2*sr*lam(5)*lam(7) -                     &
               216*m2*Rqrd*sr*lam(5)*lam(7) - 18*lam(6)*lam(7) +                &
               216*m2*sr*lam(6)*lam(7) + 216*m2*Rqrd*sr*lam(6)*lam(7) -         &
               9*lam(7)**2 + 108*m2*sr*lam(7)**2 +                              &
               108*m2*Rqrd*sr*lam(7)**2 - 576*lam(2)*lam(8) +                   &
               408*lam(3)*lam(8) - 144*m2*sr*lam(4)*lam(8) -                    &
               144*m2*Rqrd*sr*lam(4)*lam(8) - 144*m2*sr*lam(5)*lam(8) -         &
               144*m2*Rqrd*sr*lam(5)*lam(8) + 144*m2*sr*lam(6)*lam(8) +         &
               144*m2*Rqrd*sr*lam(6)*lam(8) + 144*m2*sr*lam(7)*lam(8) +         &
               144*m2*Rqrd*sr*lam(7)*lam(8) - 144*lam(8)**2 -                   &
               72*lam(1)*(9*(lam(5) - lam(6) + lam(7)) + 2*lam(8)) -            &
               960*m2*sr*lam(2)*lam(9) - 960*m2*Rqrd*sr*lam(2)*lam(9) -         &
               60*lam(4)*lam(9) - 60*lam(5)*lam(9) + 60*lam(6)*lam(9) +         &
               60*lam(7)*lam(9) - 480*m2*sr*lam(8)*lam(9) -                     &
               480*m2*Rqrd*sr*lam(8)*lam(9) - 100*lam(9)**2 +                   &
               384*m2*sr*lam(9)**2 + 384*m2*Rqrd*sr*lam(9)**2 +                 &
               2*qzt2*(1 + Rqrd)**2*sr*                                         &
                (36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -                &
                  72*lam(1)*(9*lam(4) -                                         &
                     9*(lam(5) - lam(6) + lam(7)) - 2*lam(8)) +                 &
                  24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                         &
                     lam(8)*(-17*lam(3) + 6*lam(8))) +                          &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) +                              &
               Rqsu*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                             &
                (36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -                &
                  72*lam(1)*(9*lam(4) -                                         &
                     9*(lam(5) - lam(6) + lam(7)) - 2*lam(8)) +                 &
                  24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                         &
                     lam(8)*(-17*lam(3) + 6*lam(8))) +                          &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) +                              &
               24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                  &
                  24*m2*(1 + Rqrd)*sr*(2*lam(2) + lam(8)) + 10*lam(9))*         &
                lam(10) - 144*lam(10)**2) +                                     &
            dRqsu*(-3*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                &
                (48*lam(2)**2 +                                                 &
                  216*lam(2)*(lam(4) - lam(5) + lam(6) - lam(7)) +              &
                  3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +              &
                     lam(6)**2 +                                                &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                 &
                     26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2         &
) + 96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) - 64*lam(8)**2 +          &
                  128*lam(3)*                                                   &
                   (3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) -                 &
                  48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) -                   &
                     3*lam(7) - 16*lam(9) + 16*lam(10))) +                      &
               4*m2*(1 + Rqsu)*ss*                                              &
                (2592*lam(1)**2 + 576*lam(3)**2 + 27*lam(4)**2 +                &
                  63*lam(5)**2 + 27*lam(6)**2 - 18*lam(6)*lam(7) +              &
                  63*lam(7)**2 + 192*lam(6)*lam(9) + 576*lam(7)*lam(9) +        &
                  320*lam(9)**2 -                                               &
                  6*lam(4)*(3*lam(5) + 9*lam(6) - 3*lam(7) +                    &
                     32*lam(9) - 32*lam(10)) -                                  &
                  64*(3*(lam(6) + lam(7)) + 5*lam(9))*lam(10) +                 &
                  96*lam(10)**2 +                                               &
                  6*lam(5)*(3*lam(6) - 21*lam(7) +                              &
                     32*(-3*lam(9) + lam(10)))))) +                             &
         expdISs*qdB*quB*qzt2*tanhds*                                           &
          (576*dRqrd*lam(2)**2 +                                                &
            dRqsd*(3*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                 &
                (3*(16*lam(2)**2 +                                              &
                     192*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) +           &
                     64*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) -            &
                     72*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +            &
                     (lam(4) + lam(5) + lam(6) + lam(7))**2) +                  &
                  256*(3*lam(1) + lam(3))*lam(10)) +                            &
               4*m2*(1 + Rqsd)*ss*                                              &
                (384*(3*lam(1) + lam(3))**2 +                                   &
                  45*(lam(4) - lam(5) - lam(6) + lam(7))**2 +                   &
                  192*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +             &
                  32*lam(10)**2)) -                                             &
            dRqrd*(1152*qzt2*sr*lam(2)**2 + 2304*qzt2*Rqrd*sr*lam(2)**2 +       &
               1152*qzt2*Rqrd**2*sr*lam(2)**2 - 648*lam(1)*lam(4) +             &
               1296*qzt2*sr*lam(1)*lam(4) +                                     &
               2592*qzt2*Rqrd*sr*lam(1)*lam(4) +                                &
               1296*qzt2*Rqrd**2*sr*lam(1)*lam(4) + 9936*lam(2)*lam(4) -        &
               19872*qzt2*sr*lam(2)*lam(4) -                                    &
               39744*qzt2*Rqrd*sr*lam(2)*lam(4) -                               &
               19872*qzt2*Rqrd**2*sr*lam(2)*lam(4) - 252*lam(3)*lam(4) +        &
               504*qzt2*sr*lam(3)*lam(4) +                                      &
               1008*qzt2*Rqrd*sr*lam(3)*lam(4) +                                &
               504*qzt2*Rqrd**2*sr*lam(3)*lam(4) - 14265*lam(4)**2 +            &
               28530*qzt2*sr*lam(4)**2 + 57060*qzt2*Rqrd*sr*lam(4)**2 +         &
               28530*qzt2*Rqrd**2*sr*lam(4)**2 - 648*lam(1)*lam(5) +            &
               1296*qzt2*sr*lam(1)*lam(5) +                                     &
               2592*qzt2*Rqrd*sr*lam(1)*lam(5) +                                &
               1296*qzt2*Rqrd**2*sr*lam(1)*lam(5) + 1296*lam(2)*lam(5) -        &
               2592*qzt2*sr*lam(2)*lam(5) -                                     &
               5184*qzt2*Rqrd*sr*lam(2)*lam(5) -                                &
               2592*qzt2*Rqrd**2*sr*lam(2)*lam(5) - 828*lam(3)*lam(5) +         &
               1656*qzt2*sr*lam(3)*lam(5) +                                     &
               3312*qzt2*Rqrd*sr*lam(3)*lam(5) +                                &
               1656*qzt2*Rqrd**2*sr*lam(3)*lam(5) - 3798*lam(4)*lam(5) +        &
               7596*qzt2*sr*lam(4)*lam(5) +                                     &
               15192*qzt2*Rqrd*sr*lam(4)*lam(5) +                               &
               7596*qzt2*Rqrd**2*sr*lam(4)*lam(5) - 117*lam(5)**2 +             &
               234*qzt2*sr*lam(5)**2 + 468*qzt2*Rqrd*sr*lam(5)**2 +             &
               234*qzt2*Rqrd**2*sr*lam(5)**2 - 648*lam(1)*lam(6) +              &
               1296*qzt2*sr*lam(1)*lam(6) +                                     &
               2592*qzt2*Rqrd*sr*lam(1)*lam(6) +                                &
               1296*qzt2*Rqrd**2*sr*lam(1)*lam(6) + 432*lam(2)*lam(6) -         &
               864*qzt2*sr*lam(2)*lam(6) -                                      &
               1728*qzt2*Rqrd*sr*lam(2)*lam(6) -                                &
               864*qzt2*Rqrd**2*sr*lam(2)*lam(6) - 252*lam(3)*lam(6) +          &
               504*qzt2*sr*lam(3)*lam(6) +                                      &
               1008*qzt2*Rqrd*sr*lam(3)*lam(6) +                                &
               504*qzt2*Rqrd**2*sr*lam(3)*lam(6) - 1278*lam(4)*lam(6) +         &
               2556*qzt2*sr*lam(4)*lam(6) +                                     &
               5112*qzt2*Rqrd*sr*lam(4)*lam(6) +                                &
               2556*qzt2*Rqrd**2*sr*lam(4)*lam(6) - 90*lam(5)*lam(6) +          &
               180*qzt2*sr*lam(5)*lam(6) +                                      &
               360*qzt2*Rqrd*sr*lam(5)*lam(6) +                                 &
               180*qzt2*Rqrd**2*sr*lam(5)*lam(6) - 9*lam(6)**2 +                &
               18*qzt2*sr*lam(6)**2 + 36*qzt2*Rqrd*sr*lam(6)**2 +               &
               18*qzt2*Rqrd**2*sr*lam(6)**2 - 648*lam(1)*lam(7) +               &
               1296*qzt2*sr*lam(1)*lam(7) +                                     &
               2592*qzt2*Rqrd*sr*lam(1)*lam(7) +                                &
               1296*qzt2*Rqrd**2*sr*lam(1)*lam(7) - 1296*lam(2)*lam(7) +        &
               2592*qzt2*sr*lam(2)*lam(7) +                                     &
               5184*qzt2*Rqrd*sr*lam(2)*lam(7) +                                &
               2592*qzt2*Rqrd**2*sr*lam(2)*lam(7) - 828*lam(3)*lam(7) +         &
               1656*qzt2*sr*lam(3)*lam(7) +                                     &
               3312*qzt2*Rqrd*sr*lam(3)*lam(7) +                                &
               1656*qzt2*Rqrd**2*sr*lam(3)*lam(7) + 3798*lam(4)*lam(7) -        &
               7596*qzt2*sr*lam(4)*lam(7) -                                     &
               15192*qzt2*Rqrd*sr*lam(4)*lam(7) -                               &
               7596*qzt2*Rqrd**2*sr*lam(4)*lam(7) + 234*lam(5)*lam(7) -         &
               468*qzt2*sr*lam(5)*lam(7) -                                      &
               936*qzt2*Rqrd*sr*lam(5)*lam(7) -                                 &
               468*qzt2*Rqrd**2*sr*lam(5)*lam(7) + 90*lam(6)*lam(7) -           &
               180*qzt2*sr*lam(6)*lam(7) -                                      &
               360*qzt2*Rqrd*sr*lam(6)*lam(7) -                                 &
               180*qzt2*Rqrd**2*sr*lam(6)*lam(7) - 117*lam(7)**2 +              &
               234*qzt2*sr*lam(7)**2 + 468*qzt2*Rqrd*sr*lam(7)**2 +             &
               234*qzt2*Rqrd**2*sr*lam(7)**2 - 144*lam(1)*lam(8) +              &
               288*qzt2*sr*lam(1)*lam(8) +                                      &
               576*qzt2*Rqrd*sr*lam(1)*lam(8) +                                 &
               288*qzt2*Rqrd**2*sr*lam(1)*lam(8) - 576*lam(2)*lam(8) +          &
               1152*qzt2*sr*lam(2)*lam(8) +                                     &
               2304*qzt2*Rqrd*sr*lam(2)*lam(8) +                                &
               1152*qzt2*Rqrd**2*sr*lam(2)*lam(8) - 384*lam(3)*lam(8) +         &
               768*qzt2*sr*lam(3)*lam(8) +                                      &
               1536*qzt2*Rqrd*sr*lam(3)*lam(8) +                                &
               768*qzt2*Rqrd**2*sr*lam(3)*lam(8) + 288*lam(8)**2 -              &
               576*qzt2*sr*lam(8)**2 - 1152*qzt2*Rqrd*sr*lam(8)**2 -            &
               576*qzt2*Rqrd**2*sr*lam(8)**2 + 2304*lam(2)*lam(9) -             &
               4608*qzt2*sr*lam(2)*lam(9) -                                     &
               9216*qzt2*Rqrd*sr*lam(2)*lam(9) -                                &
               4608*qzt2*Rqrd**2*sr*lam(2)*lam(9) - 6684*lam(4)*lam(9) +        &
               13368*qzt2*sr*lam(4)*lam(9) +                                    &
               26736*qzt2*Rqrd*sr*lam(4)*lam(9) +                               &
               13368*qzt2*Rqrd**2*sr*lam(4)*lam(9) - 564*lam(5)*lam(9) +        &
               1128*qzt2*sr*lam(5)*lam(9) +                                     &
               2256*qzt2*Rqrd*sr*lam(5)*lam(9) +                                &
               1128*qzt2*Rqrd**2*sr*lam(5)*lam(9) - 228*lam(6)*lam(9) +         &
               456*qzt2*sr*lam(6)*lam(9) +                                      &
               912*qzt2*Rqrd*sr*lam(6)*lam(9) +                                 &
               456*qzt2*Rqrd**2*sr*lam(6)*lam(9) + 564*lam(7)*lam(9) -          &
               1128*qzt2*sr*lam(7)*lam(9) -                                     &
               2256*qzt2*Rqrd*sr*lam(7)*lam(9) -                                &
               1128*qzt2*Rqrd**2*sr*lam(7)*lam(9) - 568*lam(9)**2 +             &
               1136*qzt2*sr*lam(9)**2 + 2272*qzt2*Rqrd*sr*lam(9)**2 +           &
               1136*qzt2*Rqrd**2*sr*lam(9)**2 +                                 &
               12*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                               &
                (384*lam(2) -                                                   &
                  3*(363*lam(4) + 59*lam(5) + 21*lam(6) - 59*lam(7)) -          &
                  286*lam(9))*lam(10) +                                         &
               2616*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(10)**2 -                 &
               4*m2*(1 + Rqrd)*sr*                                              &
                (864*lam(2)**2 + 7173*lam(4)**2 + 297*lam(5)**2 +               &
                  162*lam(5)*lam(6) + 45*lam(6)**2 - 594*lam(5)*lam(7) -        &
                  162*lam(6)*lam(7) + 297*lam(7)**2 -                           &
                  180*lam(5)*lam(8) - 36*lam(6)*lam(8) +                        &
                  180*lam(7)*lam(8) + 576*lam(5)*lam(9) +                       &
                  192*lam(6)*lam(9) - 576*lam(7)*lam(9) -                       &
                  240*lam(8)*lam(9) + 608*lam(9)**2 +                           &
                  6*lam(4)*(297*lam(5) + 93*lam(6) - 297*lam(7) +               &
                     6*lam(8) + 544*(lam(9) - 2*lam(10))) -                     &
                  24*lam(2)*(213*lam(4) + 15*lam(5) + 3*lam(6) -                &
                     15*lam(7) + 38*lam(9) - 126*lam(10)) -                     &
                  8*(144*lam(5) + 48*lam(6) - 144*lam(7) + 27*lam(8) +          &
                     280*lam(9))*lam(10) + 2240*lam(10)**2) +                   &
               Rqsd*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                             &
                (576*lam(2)**2 + 14265*lam(4)**2 + 3798*lam(4)*lam(5) +         &
                  117*lam(5)**2 + 1278*lam(4)*lam(6) + 90*lam(5)*lam(6) +       &
                  9*lam(6)**2 - 3798*lam(4)*lam(7) - 234*lam(5)*lam(7) -        &
                  90*lam(6)*lam(7) + 117*lam(7)**2 - 288*lam(8)**2 +            &
                  72*lam(1)*(9*(lam(4) + lam(5) + lam(6) + lam(7)) +            &
                     2*lam(8)) +                                                &
                  12*lam(3)*(21*lam(4) + 69*lam(5) + 21*lam(6) +                &
                     69*lam(7) + 32*lam(8)) + 6684*lam(4)*lam(9) +              &
                  564*lam(5)*lam(9) + 228*lam(6)*lam(9) -                       &
                  564*lam(7)*lam(9) + 568*lam(9)**2 -                           &
                  12*(1089*lam(4) + 177*lam(5) + 63*lam(6) -                    &
                     177*lam(7) + 286*lam(9))*lam(10) + 2616*lam(10)**2 -       &
                  144*lam(2)*(69*lam(4) + 9*lam(5) + 3*lam(6) -                 &
                     9*lam(7) - 4*(lam(8) - 4*lam(9) + 8*lam(10))))))) +        &
      expuISr*tanhur*(-2*expdISs*qxy2*(1 + Rqsd)*tanhds**2*                     &
          (dRqru*qdB*qxy2*(1 + Rqru)**2*                                        &
             (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +                  &
               72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -             &
                  2*lam(8)) +                                                   &
               24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                            &
                  lam(8)*(17*lam(3) + 6*lam(8))) +                              &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                     &
                  10*lam(9) - 12*lam(10))**2) +                                 &
            dRqsd*quB*(3*qzt2*(1 + Rqru)*(1 + Rqsd)*                            &
                (48*lam(2)**2 +                                                 &
                  216*lam(2)*(lam(4) - lam(5) + lam(6) - lam(7)) +              &
                  3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +              &
                     lam(6)**2 +                                                &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                 &
                     26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2         &
) + 96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) - 64*lam(8)**2 +          &
                  128*lam(3)*                                                   &
                   (3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) -                 &
                  48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) -                   &
                     3*lam(7) - 16*lam(9) + 16*lam(10))) +                      &
               2*m2*(-2592*lam(1)**2 - 576*lam(3)**2 - 27*lam(4)**2 -           &
                  63*lam(5)**2 - 27*lam(6)**2 + 18*lam(6)*lam(7) -              &
                  63*lam(7)**2 - 192*lam(6)*lam(9) - 576*lam(7)*lam(9) -        &
                  320*lam(9)**2 +                                               &
                  6*lam(4)*(3*lam(5) + 9*lam(6) - 3*lam(7) +                    &
                     32*lam(9) - 32*lam(10)) +                                  &
                  64*(3*(lam(6) + lam(7)) + 5*lam(9))*lam(10) -                 &
                  96*lam(10)**2 -                                               &
                  6*lam(5)*(3*lam(6) - 21*lam(7) +                              &
                     32*(-3*lam(9) + lam(10)))))) +                             &
         expdISs*qdB*quB*qzt2*tanhds*                                           &
          (-576*dRqru*lam(2)**2 +                                               &
            dRqru*(648*lam(1)*lam(4) - 288*m2*sr*lam(2)*lam(4) -                &
               288*m2*Rqru*sr*lam(2)*lam(4) - 36*lam(3)*lam(4) -                &
               9*lam(4)**2 + 108*m2*sr*lam(4)**2 +                              &
               108*m2*Rqru*sr*lam(4)**2 - 288*m2*sr*lam(2)*lam(5) -             &
               288*m2*Rqru*sr*lam(2)*lam(5) + 36*lam(3)*lam(5) -                &
               18*lam(4)*lam(5) + 216*m2*sr*lam(4)*lam(5) +                     &
               216*m2*Rqru*sr*lam(4)*lam(5) - 9*lam(5)**2 +                     &
               108*m2*sr*lam(5)**2 + 108*m2*Rqru*sr*lam(5)**2 +                 &
               288*m2*sr*lam(2)*lam(6) + 288*m2*Rqru*sr*lam(2)*lam(6) -         &
               36*lam(3)*lam(6) + 18*lam(4)*lam(6) -                            &
               216*m2*sr*lam(4)*lam(6) - 216*m2*Rqru*sr*lam(4)*lam(6) +         &
               18*lam(5)*lam(6) - 216*m2*sr*lam(5)*lam(6) -                     &
               216*m2*Rqru*sr*lam(5)*lam(6) - 9*lam(6)**2 +                     &
               108*m2*sr*lam(6)**2 + 108*m2*Rqru*sr*lam(6)**2 +                 &
               288*m2*sr*lam(2)*lam(7) + 288*m2*Rqru*sr*lam(2)*lam(7) +         &
               36*lam(3)*lam(7) + 18*lam(4)*lam(7) -                            &
               216*m2*sr*lam(4)*lam(7) - 216*m2*Rqru*sr*lam(4)*lam(7) +         &
               18*lam(5)*lam(7) - 216*m2*sr*lam(5)*lam(7) -                     &
               216*m2*Rqru*sr*lam(5)*lam(7) - 18*lam(6)*lam(7) +                &
               216*m2*sr*lam(6)*lam(7) + 216*m2*Rqru*sr*lam(6)*lam(7) -         &
               9*lam(7)**2 + 108*m2*sr*lam(7)**2 +                              &
               108*m2*Rqru*sr*lam(7)**2 - 576*lam(2)*lam(8) +                   &
               408*lam(3)*lam(8) - 144*m2*sr*lam(4)*lam(8) -                    &
               144*m2*Rqru*sr*lam(4)*lam(8) - 144*m2*sr*lam(5)*lam(8) -         &
               144*m2*Rqru*sr*lam(5)*lam(8) + 144*m2*sr*lam(6)*lam(8) +         &
               144*m2*Rqru*sr*lam(6)*lam(8) + 144*m2*sr*lam(7)*lam(8) +         &
               144*m2*Rqru*sr*lam(7)*lam(8) - 144*lam(8)**2 -                   &
               72*lam(1)*(9*(lam(5) - lam(6) + lam(7)) + 2*lam(8)) -            &
               960*m2*sr*lam(2)*lam(9) - 960*m2*Rqru*sr*lam(2)*lam(9) -         &
               60*lam(4)*lam(9) - 60*lam(5)*lam(9) + 60*lam(6)*lam(9) +         &
               60*lam(7)*lam(9) - 480*m2*sr*lam(8)*lam(9) -                     &
               480*m2*Rqru*sr*lam(8)*lam(9) - 100*lam(9)**2 +                   &
               384*m2*sr*lam(9)**2 + 384*m2*Rqru*sr*lam(9)**2 +                 &
               2*qzt2*(1 + Rqru)**2*sr*                                         &
                (36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -                &
                  72*lam(1)*(9*lam(4) -                                         &
                     9*(lam(5) - lam(6) + lam(7)) - 2*lam(8)) +                 &
                  24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                         &
                     lam(8)*(-17*lam(3) + 6*lam(8))) +                          &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) +                              &
               Rqsd*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                             &
                (36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -                &
                  72*lam(1)*(9*lam(4) -                                         &
                     9*(lam(5) - lam(6) + lam(7)) - 2*lam(8)) +                 &
                  24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                         &
                     lam(8)*(-17*lam(3) + 6*lam(8))) +                          &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) +                              &
               24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                  &
                  24*m2*(1 + Rqru)*sr*(2*lam(2) + lam(8)) + 10*lam(9))*         &
                lam(10) - 144*lam(10)**2) +                                     &
            dRqsd*(-3*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                &
                (48*lam(2)**2 +                                                 &
                  216*lam(2)*(lam(4) - lam(5) + lam(6) - lam(7)) +              &
                  3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +              &
                     lam(6)**2 +                                                &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                 &
                     26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2         &
) + 96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) - 64*lam(8)**2 +          &
                  128*lam(3)*                                                   &
                   (3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) -                 &
                  48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) -                   &
                     3*lam(7) - 16*lam(9) + 16*lam(10))) +                      &
               4*m2*(1 + Rqsd)*ss*                                              &
                (2592*lam(1)**2 + 576*lam(3)**2 + 27*lam(4)**2 +                &
                  63*lam(5)**2 + 27*lam(6)**2 - 18*lam(6)*lam(7) +              &
                  63*lam(7)**2 + 192*lam(6)*lam(9) + 576*lam(7)*lam(9) +        &
                  320*lam(9)**2 -                                               &
                  6*lam(4)*(3*lam(5) + 9*lam(6) - 3*lam(7) +                    &
                     32*lam(9) - 32*lam(10)) -                                  &
                  64*(3*(lam(6) + lam(7)) + 5*lam(9))*lam(10) +                 &
                  96*lam(10)**2 +                                               &
                  6*lam(5)*(3*lam(6) - 21*lam(7) +                              &
                     32*(-3*lam(9) + lam(10)))))) +                             &
         qdB*(dRqsu*expuISs*tanhus*                                             &
             (3*qzt2*(1 + Rqru)*                                                &
                (quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +                           &
                  2*qxy2*(1 + Rqsu)**2*tanhus)*                                 &
                (3*(16*lam(2)**2 +                                              &
                     192*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) +           &
                     64*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) -            &
                     72*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +            &
                     (lam(4) + lam(5) + lam(6) + lam(7))**2) +                  &
                  256*(3*lam(1) + lam(3))*lam(10)) +                            &
               4*m2*(1 + Rqsu)*(quB*qzt2*ss + qxy2*tanhus)*                     &
                (384*(3*lam(1) + lam(3))**2 +                                   &
                  45*(lam(4) - lam(5) - lam(6) + lam(7))**2 +                   &
                  192*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +             &
                  32*lam(10)**2)) +                                             &
            dRqru*(2*expdISs*qxy2*(1 + Rqru)*                                   &
                (-12*m2*(216*lam(1)**2 - 24*lam(1)*lam(3) +                     &
                     34*lam(3)**2 + 27*lam(5)**2 + 27*lam(6)**2 +               &
                     27*(lam(4) - lam(7))**2 - 12*lam(6)*lam(8) +               &
                     34*lam(8)**2 + 6*lam(5)*(-9*lam(6) + 2*lam(8)) +           &
                     20*lam(8)*lam(9) + 48*lam(9)**2 +                          &
                     4*lam(2)*                                                  &
                      (3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +              &
                        10*lam(9))) +                                           &
                  (qxy2 + qzt2)*(1 + Rqru)*(1 + Rqsd)*                          &
                   (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +            &
                     72*lam(1)*                                                 &
                      (9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -                &
                        2*lam(8)) +                                             &
                     24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                      &
                        lam(8)*(17*lam(3) + 6*lam(8))) +                        &
                     (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +               &
                        10*lam(9) - 12*lam(10))**2) +                           &
                  288*m2*(2*lam(2) + lam(8))*lam(10)) -                         &
               4*expuISs*m2*(1 + Rqru)*                                         &
                (6*qxy2*(-144*lam(2)**2 +                                       &
                     4*(54*lam(1)**2 + 42*lam(1)*lam(3) +                       &
                        23*lam(3)**2) - 1161*lam(4)**2 +                        &
                     135*lam(5)**2 + 54*lam(5)*lam(6) + 27*lam(6)**2 +          &
                     135*lam(7)**2 + 60*lam(5)*lam(8) +                         &
                     12*lam(6)*lam(8) + 20*lam(8)**2 +                          &
                     40*lam(8)*lam(9) + 48*lam(9)**2 +                          &
                     4*lam(2)*                                                  &
                      (213*lam(4) + 15*lam(5) + 3*lam(6) - 15*lam(7) +          &
                        38*lam(9) - 126*lam(10)) -                              &
                     18*lam(4)*                                                 &
                      (18*lam(5) + 6*lam(6) - 21*lam(7) + 32*lam(9) -           &
                        64*lam(10)) +                                           &
                     12*(3*lam(8) - 8*lam(9))*lam(10) + 96*lam(10)**2) -        &
                  quB*qzt2*sr*tanhus*                                           &
                   (864*lam(2)**2 + 7173*lam(4)**2 + 297*lam(5)**2 +            &
                     162*lam(5)*lam(6) + 45*lam(6)**2 -                         &
                     594*lam(5)*lam(7) - 162*lam(6)*lam(7) +                    &
                     297*lam(7)**2 - 180*lam(5)*lam(8) -                        &
                     36*lam(6)*lam(8) + 180*lam(7)*lam(8) +                     &
                     576*lam(5)*lam(9) + 192*lam(6)*lam(9) -                    &
                     576*lam(7)*lam(9) - 240*lam(8)*lam(9) +                    &
                     608*lam(9)**2 +                                            &
                     6*lam(4)*                                                  &
                      (297*lam(5) + 93*lam(6) - 297*lam(7) + 6*lam(8) +         &
                        544*(lam(9) - 2*lam(10))) -                             &
                     24*lam(2)*                                                 &
                      (213*lam(4) + 15*lam(5) + 3*lam(6) - 15*lam(7) +          &
                        38*lam(9) - 126*lam(10)) -                              &
                     8*(144*lam(5) + 48*lam(6) - 144*lam(7) +                   &
                        27*lam(8) + 280*lam(9))*lam(10) + 2240*lam(10)**2       &
)) + expuISs*(1 + Rqsu)*(-2*qxy2*qzt2*(1 + Rqru)**2*                            &
                   (576*lam(2)**2 + 14265*lam(4)**2 +                           &
                     3798*lam(4)*lam(5) + 117*lam(5)**2 +                       &
                     1278*lam(4)*lam(6) + 90*lam(5)*lam(6) +                    &
                     9*lam(6)**2 - 3798*lam(4)*lam(7) -                         &
                     234*lam(5)*lam(7) - 90*lam(6)*lam(7) +                     &
                     117*lam(7)**2 +                                            &
                     72*lam(1)*                                                 &
                      (-9*(lam(4) + lam(5) + lam(6) + lam(7)) -                 &
                        2*lam(8)) - 288*lam(8)**2 -                             &
                     12*lam(3)*                                                 &
                      (21*lam(4) + 69*lam(5) + 21*lam(6) + 69*lam(7) +          &
                        32*lam(8)) + 6684*lam(4)*lam(9) +                       &
                     564*lam(5)*lam(9) + 228*lam(6)*lam(9) -                    &
                     564*lam(7)*lam(9) + 568*lam(9)**2 -                        &
                     12*(1089*lam(4) + 177*lam(5) + 63*lam(6) -                 &
                        177*lam(7) + 286*lam(9))*lam(10) +                      &
                     2616*lam(10)**2 -                                          &
                     144*lam(2)*                                                &
                      (69*lam(4) + 9*lam(5) + 3*lam(6) - 9*lam(7) -             &
                        4*(lam(8) - 4*lam(9) + 8*lam(10)))) +                   &
                  2*qxy2**2*(1 + Rqru)**2*(-1 + tanhus**2)*                     &
                   (576*lam(2)**2 + 14265*lam(4)**2 +                           &
                     3798*lam(4)*lam(5) + 117*lam(5)**2 +                       &
                     1278*lam(4)*lam(6) + 90*lam(5)*lam(6) +                    &
                     9*lam(6)**2 - 3798*lam(4)*lam(7) -                         &
                     234*lam(5)*lam(7) - 90*lam(6)*lam(7) +                     &
                     117*lam(7)**2 +                                            &
                     72*lam(1)*                                                 &
                      (-9*(lam(4) + lam(5) + lam(6) + lam(7)) -                 &
                        2*lam(8)) - 288*lam(8)**2 -                             &
                     12*lam(3)*                                                 &
                      (21*lam(4) + 69*lam(5) + 21*lam(6) + 69*lam(7) +          &
                        32*lam(8)) + 6684*lam(4)*lam(9) +                       &
                     564*lam(5)*lam(9) + 228*lam(6)*lam(9) -                    &
                     564*lam(7)*lam(9) + 568*lam(9)**2 -                        &
                     12*(1089*lam(4) + 177*lam(5) + 63*lam(6) -                 &
                        177*lam(7) + 286*lam(9))*lam(10) +                      &
                     2616*lam(10)**2 -                                          &
                     144*lam(2)*                                                &
                      (69*lam(4) + 9*lam(5) + 3*lam(6) - 9*lam(7) -             &
                        4*(lam(8) - 4*lam(9) + 8*lam(10)))) -                   &
                  quB*qzt2*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*tanhus*               &
                   (576*lam(2)**2 + 14265*lam(4)**2 +                           &
                     3798*lam(4)*lam(5) + 117*lam(5)**2 +                       &
                     1278*lam(4)*lam(6) + 90*lam(5)*lam(6) +                    &
                     9*lam(6)**2 - 3798*lam(4)*lam(7) -                         &
                     234*lam(5)*lam(7) - 90*lam(6)*lam(7) +                     &
                     117*lam(7)**2 - 288*lam(8)**2 +                            &
                     72*lam(1)*                                                 &
                      (9*(lam(4) + lam(5) + lam(6) + lam(7)) + 2*lam(8))        &
+ 12*lam(3)*(21*lam(4) + 69*lam(5) + 21*lam(6) + 69*lam(7) + 32*lam(8)) +       &
                     6684*lam(4)*lam(9) + 564*lam(5)*lam(9) +                   &
                     228*lam(6)*lam(9) - 564*lam(7)*lam(9) +                    &
                     568*lam(9)**2 -                                            &
                     12*(1089*lam(4) + 177*lam(5) + 63*lam(6) -                 &
                        177*lam(7) + 286*lam(9))*lam(10) +                      &
                     2616*lam(10)**2 -                                          &
                     144*lam(2)*                                                &
                      (69*lam(4) + 9*lam(5) + 3*lam(6) - 9*lam(7) -             &
                        4*(lam(8) - 4*lam(9) + 8*lam(10)))))))) -               &
      expuISr*qxy2*tanhur**2*(6*dRqsd*expdISs*quB*qxy2*(1 + Rqru)*              &
          (1 + Rqsd)**2*tanhds**3*                                              &
          (48*lam(2)**2 + 216*lam(2)*                                           &
             (lam(4) - lam(5) + lam(6) - lam(7)) +                              &
            3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) + lam(6)**2 +        &
               2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                       &
               26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +            &
            96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -                 &
            64*lam(8)**2 - 128*lam(3)*                                          &
             (3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) +                       &
            48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -              &
               16*lam(9) + 16*lam(10))) -                                       &
         expdISs*qdB*quB*tanhds**2*                                             &
          (dRqru*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                     &
             (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +                  &
               72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -             &
                  2*lam(8)) +                                                   &
               24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                            &
                  lam(8)*(17*lam(3) + 6*lam(8))) +                              &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -         &
                  12*lam(10))**2) -                                             &
            3*dRqsd*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                  &
             (48*lam(2)**2 + 216*lam(2)*                                        &
                (lam(4) - lam(5) + lam(6) - lam(7)) +                           &
               3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                 &
                  lam(6)**2 + 2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +        &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -              &
               64*lam(8)**2 -                                                   &
               128*lam(3)*(3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) +          &
               48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -           &
                  16*lam(9) + 16*lam(10)))) -                                   &
         2*expdISs*(1 + Rqru)*tanhds*                                           &
          (3*dRqsd*quB*qxy2*(1 + Rqsd)**2*                                      &
             (48*lam(2)**2 + 216*lam(2)*                                        &
                (lam(4) - lam(5) + lam(6) - lam(7)) +                           &
               3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                 &
                  lam(6)**2 +                                                   &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                    &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -              &
               64*lam(8)**2 -                                                   &
               128*lam(3)*(3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) +          &
               48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -           &
                  16*lam(9) + 16*lam(10))) +                                    &
            dRqru*qdB*(qzt2*(1 + Rqru)*(1 + Rqsd)*                              &
                (36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -                &
                  72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -          &
                     2*lam(8)) +                                                &
                  24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                         &
                     lam(8)*(-17*lam(3) + 6*lam(8))) +                          &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) -                              &
               6*m2*(-9*lam(4)**2 - 9*lam(5)**2 + 18*lam(5)*lam(6) -            &
                  9*lam(6)**2 + 18*lam(5)*lam(7) - 18*lam(6)*lam(7) -           &
                  9*lam(7)**2 + 12*lam(5)*lam(8) - 12*lam(6)*lam(8) -           &
                  12*lam(7)*lam(8) +                                            &
                  6*lam(4)*(-3*lam(5) + 3*(lam(6) + lam(7)) + 2*lam(8)) +       &
                  40*lam(8)*lam(9) - 32*lam(9)**2 +                             &
                  8*lam(2)*(3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +         &
                     10*lam(9) - 12*lam(10)) - 48*lam(8)*lam(10)))) +           &
         qdB*(-3*(1 + Rqru)*(dRqsu*expuISs*(-1 + tanhus)*(1 + tanhus)*          &
                (quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +                           &
                  2*qxy2*(1 + Rqsu)**2*tanhus)*                                 &
                (3*(16*lam(2)**2 -                                              &
                     192*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) -           &
                     64*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) -            &
                     72*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +            &
                     (lam(4) + lam(5) + lam(6) + lam(7))**2) -                  &
                  256*(3*lam(1) + lam(3))*lam(10)) +                            &
               dRqsd*expdISs*quB*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                &
                (48*lam(2)**2 +                                                 &
                  216*lam(2)*(lam(4) - lam(5) + lam(6) - lam(7)) +              &
                  3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +              &
                     lam(6)**2 +                                                &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                 &
                     26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2)        &
+ 96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) - 64*lam(8)**2 -            &
                  128*lam(3)*(3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) +       &
                  48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -        &
                     16*lam(9) + 16*lam(10)))) +                                &
            dRqru*(expdISs*quB*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*       &
                (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +               &
                  72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -          &
                     2*lam(8)) +                                                &
                  24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                         &
                     lam(8)*(17*lam(3) + 6*lam(8))) +                           &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) +                              &
               expuISs*(quB*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*          &
                   (-1 + tanhus**2)*                                            &
                   (576*lam(2)**2 + 14265*lam(4)**2 +                           &
                     3798*lam(4)*lam(5) + 117*lam(5)**2 +                       &
                     1278*lam(4)*lam(6) + 90*lam(5)*lam(6) +                    &
                     9*lam(6)**2 - 3798*lam(4)*lam(7) -                         &
                     234*lam(5)*lam(7) - 90*lam(6)*lam(7) +                     &
                     117*lam(7)**2 +                                            &
                     72*lam(1)*                                                 &
                      (-9*(lam(4) + lam(5) + lam(6) + lam(7)) - 2*lam(8))       &
- 288*lam(8)**2 - 12*lam(3)*(21*lam(4) + 69*lam(5) + 21*lam(6) +                &
                        69*lam(7) + 32*lam(8)) + 6684*lam(4)*lam(9) +           &
                     564*lam(5)*lam(9) + 228*lam(6)*lam(9) -                    &
                     564*lam(7)*lam(9) + 568*lam(9)**2 -                        &
                     12*(1089*lam(4) + 177*lam(5) + 63*lam(6) -                 &
                        177*lam(7) + 286*lam(9))*lam(10) +                      &
                     2616*lam(10)**2 -                                          &
                     144*lam(2)*                                                &
                      (69*lam(4) + 9*lam(5) + 3*lam(6) - 9*lam(7) -             &
                        4*(lam(8) - 4*lam(9) + 8*lam(10)))) +                   &
                  2*(1 + Rqru)*tanhus*                                          &
                   (2*m2*(-864*lam(2)**2 - 7173*lam(4)**2 -                     &
                        297*lam(5)**2 - 162*lam(5)*lam(6) -                     &
                        45*lam(6)**2 + 594*lam(5)*lam(7) +                      &
                        162*lam(6)*lam(7) - 297*lam(7)**2 +                     &
                        180*lam(5)*lam(8) + 36*lam(6)*lam(8) -                  &
                        180*lam(7)*lam(8) - 576*lam(5)*lam(9) -                 &
                        192*lam(6)*lam(9) + 576*lam(7)*lam(9) +                 &
                        240*lam(8)*lam(9) - 608*lam(9)**2 -                     &
                        6*lam(4)*                                               &
                         (297*lam(5) + 93*lam(6) - 297*lam(7) +                 &
                         6*lam(8) + 544*(lam(9) - 2*lam(10))) +                 &
                        24*lam(2)*                                              &
                         (213*lam(4) + 15*lam(5) + 3*lam(6) -                   &
                         15*lam(7) + 38*lam(9) - 126*lam(10)) +                 &
                        8*(144*lam(5) + 48*lam(6) - 144*lam(7) +                &
                         27*lam(8) + 280*lam(9))*lam(10) - 2240*lam(10)**2      &
) + qzt2*(1 + Rqru)*(1 + Rqsu)*                                                 &
                      (576*lam(2)**2 + 14265*lam(4)**2 +                        &
                        3798*lam(4)*lam(5) + 117*lam(5)**2 +                    &
                        1278*lam(4)*lam(6) + 90*lam(5)*lam(6) +                 &
                        9*lam(6)**2 - 3798*lam(4)*lam(7) -                      &
                        234*lam(5)*lam(7) - 90*lam(6)*lam(7) +                  &
                        117*lam(7)**2 - 288*lam(8)**2 +                         &
                        72*lam(1)*                                              &
                         (9*(lam(4) + lam(5) + lam(6) + lam(7)) +               &
                         2*lam(8)) +                                            &
                        12*lam(3)*                                              &
                         (21*lam(4) + 69*lam(5) + 21*lam(6) + 69*lam(7) +       &
                         32*lam(8)) + 6684*lam(4)*lam(9) +                      &
                        564*lam(5)*lam(9) + 228*lam(6)*lam(9) -                 &
                        564*lam(7)*lam(9) + 568*lam(9)**2 -                     &
                        12*(1089*lam(4) + 177*lam(5) + 63*lam(6) -              &
                         177*lam(7) + 286*lam(9))*lam(10) +                     &
                        2616*lam(10)**2 -                                       &
                        144*lam(2)*                                             &
                         (69*lam(4) + 9*lam(5) + 3*lam(6) - 9*lam(7) -          &
                          4*(lam(8) - 4*lam(9) + 8*lam(10))))))))))/1296.
  dtLAM(5)=(6*dRqrd*expdISr*quB*qxy2**2*(1 + Rqrd)**2*tanhdr**3*             &
       (-(expuISs*(1 + Rqsu)*(-1 + tanhus**2)*                                  &
            (576*lam(2)**2 + 9*lam(4)**2 + 882*lam(4)*lam(5) +                  &
              11241*lam(5)**2 - 18*lam(4)*lam(6) - 882*lam(5)*lam(6) +          &
              9*lam(6)**2 - 18*lam(4)*lam(7) - 882*lam(5)*lam(7) +              &
              18*lam(6)*lam(7) + 9*lam(7)**2 - 144*lam(8)**2 +                  &
              72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +              &
                 2*lam(8)) + 36*lam(3)*                                         &
               (lam(4) - lam(5) + lam(6) - lam(7) + 10*lam(8)) +                &
              132*lam(4)*lam(9) + 4740*lam(5)*lam(9) -                          &
              132*lam(6)*lam(9) - 132*lam(7)*lam(9) + 412*lam(9)**2 -           &
              96*lam(2)*(3*lam(4) + 75*lam(5) -                                 &
                 3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +                  &
              24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                   &
                 10*lam(9))*lam(10) - 144*lam(10)**2)) +                        &
         expdISs*(1 + Rqsd)*(-1 + tanhds**2)*                                   &
          (576*lam(2)**2 + 9*lam(4)**2 - 522*lam(4)*lam(5) -                    &
            5499*lam(5)**2 - 18*lam(4)*lam(6) + 522*lam(5)*lam(6) +             &
            9*lam(6)**2 + 90*lam(4)*lam(7) + 198*lam(5)*lam(7) -                &
            90*lam(6)*lam(7) + 117*lam(7)**2 +                                  &
            36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +           &
            288*lam(8)**2 + 72*lam(1)*                                          &
             (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) + 2*lam(8)) -           &
            156*lam(4)*lam(9) - 2100*lam(5)*lam(9) + 156*lam(6)*lam(9) -        &
            204*lam(7)*lam(9) - 56*lam(9)**2 +                                  &
            48*lam(2)*(3*lam(4) + 75*lam(5) -                                   &
               3*(lam(6) + lam(7) + 4*lam(8)) + 16*lam(9)) -                    &
            36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) - 10*lam(9))*       &
             lam(10) + 456*lam(10)**2)) -                                       &
      6*dRqru*expuISr*qdB*qxy2**2*(1 + Rqru)**2*tanhur**3*                      &
       (expdISs*(1 + Rqsd)*(-1 + tanhds**2)*                                    &
          (576*lam(2)**2 + 9*lam(4)**2 + 882*lam(4)*lam(5) +                    &
            11241*lam(5)**2 - 18*lam(4)*lam(6) - 882*lam(5)*lam(6) +            &
            9*lam(6)**2 - 18*lam(4)*lam(7) - 882*lam(5)*lam(7) +                &
            18*lam(6)*lam(7) + 9*lam(7)**2 - 144*lam(8)**2 +                    &
            72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +                &
               2*lam(8)) + 36*lam(3)*                                           &
             (lam(4) - lam(5) + lam(6) - lam(7) + 10*lam(8)) +                  &
            132*lam(4)*lam(9) + 4740*lam(5)*lam(9) - 132*lam(6)*lam(9) -        &
            132*lam(7)*lam(9) + 412*lam(9)**2 -                                 &
            96*lam(2)*(3*lam(4) + 75*lam(5) -                                   &
               3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +                    &
            24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) + 10*lam(9))*         &
             lam(10) - 144*lam(10)**2) -                                        &
         expuISs*(1 + Rqsu)*(-1 + tanhus**2)*                                   &
          (576*lam(2)**2 + 9*lam(4)**2 - 522*lam(4)*lam(5) -                    &
            5499*lam(5)**2 - 18*lam(4)*lam(6) + 522*lam(5)*lam(6) +             &
            9*lam(6)**2 + 90*lam(4)*lam(7) + 198*lam(5)*lam(7) -                &
            90*lam(6)*lam(7) + 117*lam(7)**2 +                                  &
            36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +           &
            288*lam(8)**2 + 72*lam(1)*                                          &
             (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) + 2*lam(8)) -           &
            156*lam(4)*lam(9) - 2100*lam(5)*lam(9) + 156*lam(6)*lam(9) -        &
            204*lam(7)*lam(9) - 56*lam(9)**2 +                                  &
            48*lam(2)*(3*lam(4) + 75*lam(5) -                                   &
               3*(lam(6) + lam(7) + 4*lam(8)) + 16*lam(9)) -                    &
            36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) - 10*lam(9))*       &
             lam(10) + 456*lam(10)**2)) -                                       &
      6*dRqsu*expuISs*qdB*qxy2**2*(1 + Rqsu)**2*tanhus**3*                      &
       (-(expuISr*(1 + Rqru)*(144*lam(2)**2 -                                   &
              576*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) -                  &
              192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +                  &
              216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                  &
              9*(lam(4) + lam(5) + lam(6) + lam(7))**2 -                        &
              256*(3*lam(1) + lam(3))*lam(10))) +                               &
         expdISr*(1 + Rqrd)*(144*lam(2)**2 +                                    &
            216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +                &
            9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                           &
               2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                       &
               10*lam(6)*lam(7) + 13*lam(7)**2 +                                &
               lam(5)*(-10*lam(6) + 26*lam(7))) +                               &
            288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -                &
            192*lam(8)**2 + 128*lam(3)*                                         &
             (3*lam(4) + 6*lam(5) - 3*(lam(6) + 2*lam(7) + 4*lam(9)) +          &
               7*lam(10)) - 48*lam(1)*                                          &
             (3*lam(4) + 33*lam(5) - 3*(lam(6) + 11*lam(7) + 16*lam(9)) +       &
               16*lam(10)))) +                                                  &
      6*dRqsd*expdISs*quB*qxy2**2*(1 + Rqsd)**2*tanhds**3*                      &
       (expdISr*(1 + Rqrd)*(144*lam(2)**2 -                                     &
            576*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) -                    &
            192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +                    &
            216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                    &
            9*(lam(4) + lam(5) + lam(6) + lam(7))**2 -                          &
            256*(3*lam(1) + lam(3))*lam(10)) -                                  &
         expuISr*(1 + Rqru)*(144*lam(2)**2 +                                    &
            216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +                &
            9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                           &
               2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                       &
               10*lam(6)*lam(7) + 13*lam(7)**2 +                                &
               lam(5)*(-10*lam(6) + 26*lam(7))) +                               &
            288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -                &
            192*lam(8)**2 + 128*lam(3)*                                         &
             (3*lam(4) + 6*lam(5) - 3*(lam(6) + 2*lam(7) + 4*lam(9)) +          &
               7*lam(10)) - 48*lam(1)*                                          &
             (3*lam(4) + 33*lam(5) - 3*(lam(6) + 11*lam(7) + 16*lam(9)) +       &
               16*lam(10)))) -                                                  &
      3*dRqrd*expdISr*qdB*quB*                                                  &
       (expuISs*(qxy2 + qzt2)*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*        &
          (576*lam(2)**2 + 9*lam(4)**2 + 882*lam(4)*lam(5) +                    &
            11241*lam(5)**2 - 18*lam(4)*lam(6) - 882*lam(5)*lam(6) +            &
            9*lam(6)**2 - 18*lam(4)*lam(7) - 882*lam(5)*lam(7) +                &
            18*lam(6)*lam(7) + 9*lam(7)**2 - 144*lam(8)**2 +                    &
            72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +                &
               2*lam(8)) + 36*lam(3)*                                           &
             (lam(4) - lam(5) + lam(6) - lam(7) + 10*lam(8)) +                  &
            132*lam(4)*lam(9) + 4740*lam(5)*lam(9) - 132*lam(6)*lam(9) -        &
            132*lam(7)*lam(9) + 412*lam(9)**2 -                                 &
            96*lam(2)*(3*lam(4) + 75*lam(5) -                                   &
               3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +                    &
            24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) + 10*lam(9))*         &
             lam(10) - 144*lam(10)**2) +                                        &
         48*expdISs*m2*qzt2*(1 + Rqrd)*sr*                                      &
          (24*lam(2)**2 + 6*(6*lam(1)**2 + 14*lam(1)*lam(3) +                   &
               5*lam(3)**2) + 234*lam(5)**2 - 18*lam(5)*lam(6) +                &
            9*lam(4)*(7*lam(5) + lam(6)) + 99*lam(5)*lam(7) +                   &
            45*lam(6)*lam(7) + 30*lam(5)*lam(8) + 6*lam(6)*lam(8) -             &
            6*lam(8)**2 + 96*lam(5)*lam(9) + 20*lam(8)*lam(9) -                 &
            8*lam(9)**2 + 2*lam(2)*                                             &
             (3*lam(4) - 87*lam(5) - 3*lam(6) + 15*lam(7) - 26*lam(9) -         &
               30*lam(10)) + 6*(3*lam(8) + 8*lam(9))*lam(10) -                  &
            48*lam(10)**2) - expdISs*(qxy2 + qzt2)*(1 + Rqsd)*                  &
          (-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                                       &
          (576*lam(2)**2 + 9*lam(4)**2 - 522*lam(4)*lam(5) -                    &
            5499*lam(5)**2 - 18*lam(4)*lam(6) + 522*lam(5)*lam(6) +             &
            9*lam(6)**2 + 90*lam(4)*lam(7) + 198*lam(5)*lam(7) -                &
            90*lam(6)*lam(7) + 117*lam(7)**2 +                                  &
            36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +           &
            288*lam(8)**2 + 72*lam(1)*                                          &
             (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) + 2*lam(8)) -           &
            156*lam(4)*lam(9) - 2100*lam(5)*lam(9) + 156*lam(6)*lam(9) -        &
            204*lam(7)*lam(9) - 56*lam(9)**2 +                                  &
            48*lam(2)*(3*lam(4) + 75*lam(5) -                                   &
               3*(lam(6) + lam(7) + 4*lam(8)) + 16*lam(9)) -                    &
            36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -                   &
               10*lam(9))*lam(10) + 456*lam(10)**2) -                           &
         48*expuISs*m2*qzt2*(1 + Rqrd)*sr*                                      &
          (-48*lam(2)**2 + 3*(12*lam(1)**2 + 4*lam(1)*lam(3) +                  &
               5*lam(3)**2) - 468*lam(5)**2 + 36*lam(5)*lam(6) +                &
            9*lam(4)*(-5*lam(5) + lam(6)) + 45*lam(5)*lam(7) -                  &
            9*lam(6)*lam(7) - 6*lam(5)*lam(8) + 6*lam(6)*lam(8) +               &
            15*lam(8)**2 - 192*lam(5)*lam(9) - 10*lam(8)*lam(9) -               &
            8*lam(9)**2 + 12*lam(8)*lam(10) +                                   &
            lam(2)*(6*lam(4) + 294*lam(5) - 6*lam(6) - 6*lam(7) +               &
               44*lam(9) + 24*lam(10)))) +                                      &
      3*dRqru*expuISr*qdB*quB*                                                  &
       (-(expdISs*(qxy2 + qzt2)*(1 + Rqsd)*                                     &
            (-1 + 2*qzt2*(1 + Rqru)**2*sr)*                                     &
            (576*lam(2)**2 + 9*lam(4)**2 + 882*lam(4)*lam(5) +                  &
              11241*lam(5)**2 - 18*lam(4)*lam(6) - 882*lam(5)*lam(6) +          &
              9*lam(6)**2 - 18*lam(4)*lam(7) - 882*lam(5)*lam(7) +              &
              18*lam(6)*lam(7) + 9*lam(7)**2 - 144*lam(8)**2 +                  &
              72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +              &
                 2*lam(8)) + 36*lam(3)*                                         &
               (lam(4) - lam(5) + lam(6) - lam(7) + 10*lam(8)) +                &
              132*lam(4)*lam(9) + 4740*lam(5)*lam(9) -                          &
              132*lam(6)*lam(9) - 132*lam(7)*lam(9) + 412*lam(9)**2 -           &
              96*lam(2)*(3*lam(4) + 75*lam(5) -                                 &
                 3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +                  &
              24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                   &
                 10*lam(9))*lam(10) - 144*lam(10)**2)) -                        &
         48*expuISs*m2*qzt2*(1 + Rqru)*sr*                                      &
          (24*lam(2)**2 + 6*(6*lam(1)**2 + 14*lam(1)*lam(3) +                   &
               5*lam(3)**2) + 234*lam(5)**2 - 18*lam(5)*lam(6) +                &
            9*lam(4)*(7*lam(5) + lam(6)) + 99*lam(5)*lam(7) +                   &
            45*lam(6)*lam(7) + 30*lam(5)*lam(8) + 6*lam(6)*lam(8) -             &
            6*lam(8)**2 + 96*lam(5)*lam(9) + 20*lam(8)*lam(9) -                 &
            8*lam(9)**2 + 2*lam(2)*                                             &
             (3*lam(4) - 87*lam(5) - 3*lam(6) + 15*lam(7) - 26*lam(9) -         &
               30*lam(10)) + 6*(3*lam(8) + 8*lam(9))*lam(10) -                  &
            48*lam(10)**2) + expuISs*(qxy2 + qzt2)*(1 + Rqsu)*                  &
          (-1 + 2*qzt2*(1 + Rqru)**2*sr)*                                       &
          (576*lam(2)**2 + 9*lam(4)**2 - 522*lam(4)*lam(5) -                    &
            5499*lam(5)**2 - 18*lam(4)*lam(6) + 522*lam(5)*lam(6) +             &
            9*lam(6)**2 + 90*lam(4)*lam(7) + 198*lam(5)*lam(7) -                &
            90*lam(6)*lam(7) + 117*lam(7)**2 +                                  &
            36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +           &
            288*lam(8)**2 + 72*lam(1)*                                          &
             (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) + 2*lam(8)) -           &
            156*lam(4)*lam(9) - 2100*lam(5)*lam(9) + 156*lam(6)*lam(9) -        &
            204*lam(7)*lam(9) - 56*lam(9)**2 +                                  &
            48*lam(2)*(3*lam(4) + 75*lam(5) -                                   &
               3*(lam(6) + lam(7) + 4*lam(8)) + 16*lam(9)) -                    &
            36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -                   &
               10*lam(9))*lam(10) + 456*lam(10)**2) +                           &
         48*expdISs*m2*qzt2*(1 + Rqru)*sr*                                      &
          (-48*lam(2)**2 + 3*(12*lam(1)**2 + 4*lam(1)*lam(3) +                  &
               5*lam(3)**2) - 468*lam(5)**2 + 36*lam(5)*lam(6) +                &
            9*lam(4)*(-5*lam(5) + lam(6)) + 45*lam(5)*lam(7) -                  &
            9*lam(6)*lam(7) - 6*lam(5)*lam(8) + 6*lam(6)*lam(8) +               &
            15*lam(8)**2 - 192*lam(5)*lam(9) - 10*lam(8)*lam(9) -               &
            8*lam(9)**2 + 12*lam(8)*lam(10) +                                   &
            lam(2)*(6*lam(4) + 294*lam(5) - 6*lam(6) - 6*lam(7) +               &
               44*lam(9) + 24*lam(10)))) +                                      &
      3*dRqsd*expdISr*expdISs*qdB*quB*                                          &
       (-((qxy2 + qzt2)*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*              &
            (144*lam(2)**2 - 576*lam(1)*                                        &
               (lam(4) - lam(5) - lam(6) + lam(7)) -                            &
              192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +                  &
              216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                  &
              9*(lam(4) + lam(5) + lam(6) + lam(7))**2 -                        &
              256*(3*lam(1) + lam(3))*lam(10))) +                               &
         16*m2*qzt2*(1 + Rqsd)*ss*                                              &
          (32*(3*lam(1) + lam(3))**2 -                                          &
            3*(36*lam(2)**2 + 9*(lam(5) + lam(6))*(lam(4) + lam(7)) +           &
               6*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) - 8*lam(10)**2      &
))) + 3*dRqsu*expuISr*expuISs*qdB*quB*                                          &
       (-((qxy2 + qzt2)*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*              &
            (144*lam(2)**2 - 576*lam(1)*                                        &
               (lam(4) - lam(5) - lam(6) + lam(7)) -                            &
              192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +                  &
              216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                  &
              9*(lam(4) + lam(5) + lam(6) + lam(7))**2 -                        &
              256*(3*lam(1) + lam(3))*lam(10))) +                               &
         16*m2*qzt2*(1 + Rqsu)*ss*                                              &
          (32*(3*lam(1) + lam(3))**2 -                                          &
            3*(36*lam(2)**2 + 9*(lam(5) + lam(6))*(lam(4) + lam(7)) +           &
               6*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) - 8*lam(10)**2      &
))) - 6*dRqsu*expuISs*qdB*qxy2*(1 + Rqsu)*tanhus*                               &
       (expuISr*(qxy2 + qzt2)*(1 + Rqru)*(1 + Rqsu)*                            &
          (144*lam(2)**2 - 576*lam(1)*                                          &
             (lam(4) - lam(5) - lam(6) + lam(7)) -                              &
            192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +                    &
            216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                    &
            9*(lam(4) + lam(5) + lam(6) + lam(7))**2 -                          &
            256*(3*lam(1) + lam(3))*lam(10)) -                                  &
         expdISr*(qxy2 + qzt2)*(1 + Rqrd)*(1 + Rqsu)*                           &
          (144*lam(2)**2 + 216*lam(2)*                                          &
             (lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +                          &
            9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                           &
               2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                       &
               10*lam(6)*lam(7) + 13*lam(7)**2 +                                &
               lam(5)*(-10*lam(6) + 26*lam(7))) +                               &
            288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -                &
            192*lam(8)**2 + 128*lam(3)*                                         &
             (3*lam(4) + 6*lam(5) - 3*(lam(6) + 2*lam(7) + 4*lam(9)) +          &
               7*lam(10)) - 48*lam(1)*                                          &
             (3*lam(4) + 33*lam(5) -                                            &
               3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))) -              &
         8*expuISr*m2*(32*(3*lam(1) + lam(3))**2 -                              &
            3*(36*lam(2)**2 + 9*(lam(5) + lam(6))*(lam(4) + lam(7)) +           &
               6*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) -                   &
               8*lam(10)**2)) +                                                 &
         8*expdISr*m2*(72*lam(1)**2 - 384*lam(1)*lam(3) + 80*lam(3)**2 -        &
            9*(12*lam(2)**2 + 3*lam(4)*(-5*lam(5) + lam(6)) +                   &
               39*lam(5)*lam(7) - 15*lam(6)*lam(7) + 16*lam(8)**2 +             &
               2*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7) +                &
                  16*lam(8))) +                                                 &
            24*(6*lam(9)**2 - 6*lam(9)*lam(10) + lam(10)**2))) +                &
      6*dRqsd*expdISs*quB*qxy2*(1 + Rqsd)*tanhds*                               &
       (-(expdISr*(qxy2 + qzt2)*(1 + Rqrd)*(1 + Rqsd)*                          &
            (144*lam(2)**2 - 576*lam(1)*                                        &
               (lam(4) - lam(5) - lam(6) + lam(7)) -                            &
              192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +                  &
              216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                  &
              9*(lam(4) + lam(5) + lam(6) + lam(7))**2 -                        &
              256*(3*lam(1) + lam(3))*lam(10))) +                               &
         expuISr*(qxy2 + qzt2)*(1 + Rqru)*(1 + Rqsd)*                           &
          (144*lam(2)**2 + 216*lam(2)*                                          &
             (lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +                          &
            9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                           &
               2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                       &
               10*lam(6)*lam(7) + 13*lam(7)**2 +                                &
               lam(5)*(-10*lam(6) + 26*lam(7))) +                               &
            288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -                &
            192*lam(8)**2 + 128*lam(3)*                                         &
             (3*lam(4) + 6*lam(5) - 3*(lam(6) + 2*lam(7) + 4*lam(9)) +          &
               7*lam(10)) - 48*lam(1)*                                          &
             (3*lam(4) + 33*lam(5) -                                            &
               3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))) +              &
         8*expdISr*m2*(32*(3*lam(1) + lam(3))**2 -                              &
            3*(36*lam(2)**2 + 9*(lam(5) + lam(6))*(lam(4) + lam(7)) +           &
               6*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) -                   &
               8*lam(10)**2)) -                                                 &
         8*expuISr*m2*(72*lam(1)**2 - 384*lam(1)*lam(3) + 80*lam(3)**2 -        &
            9*(12*lam(2)**2 + 3*lam(4)*(-5*lam(5) + lam(6)) +                   &
               39*lam(5)*lam(7) - 15*lam(6)*lam(7) + 16*lam(8)**2 +             &
               2*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7) +                &
                  16*lam(8))) +                                                 &
            24*(6*lam(9)**2 - 6*lam(9)*lam(10) + lam(10)**2))) -                &
      3*dRqsd*expdISs*expuISr*qdB*quB*                                          &
       (-((qxy2 + qzt2)*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*              &
            (144*lam(2)**2 + 216*lam(2)*                                        &
               (lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +                        &
              9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                         &
                 2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                     &
                 10*lam(6)*lam(7) + 13*lam(7)**2 +                              &
                 lam(5)*(-10*lam(6) + 26*lam(7))) +                             &
              288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -              &
              192*lam(8)**2 +                                                   &
              128*lam(3)*(3*lam(4) + 6*lam(5) -                                 &
                 3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) -                &
              48*lam(1)*(3*lam(4) + 33*lam(5) -                                 &
                 3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10)))) +           &
         16*m2*qzt2*(1 + Rqsd)*ss*                                              &
          (72*lam(1)**2 - 384*lam(1)*lam(3) + 80*lam(3)**2 -                    &
            9*(12*lam(2)**2 + 3*lam(4)*(-5*lam(5) + lam(6)) +                   &
               39*lam(5)*lam(7) - 15*lam(6)*lam(7) + 16*lam(8)**2 +             &
               2*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7) +                &
                  16*lam(8))) +                                                 &
            24*(6*lam(9)**2 - 6*lam(9)*lam(10) + lam(10)**2))) -                &
      3*dRqsu*expdISr*expuISs*qdB*quB*                                          &
       (-((qxy2 + qzt2)*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*              &
            (144*lam(2)**2 + 216*lam(2)*                                        &
               (lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +                        &
              9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                         &
                 2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                     &
                 10*lam(6)*lam(7) + 13*lam(7)**2 +                              &
                 lam(5)*(-10*lam(6) + 26*lam(7))) +                             &
              288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -              &
              192*lam(8)**2 +                                                   &
              128*lam(3)*(3*lam(4) + 6*lam(5) -                                 &
                 3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) -                &
              48*lam(1)*(3*lam(4) + 33*lam(5) -                                 &
                 3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10)))) +           &
         16*m2*qzt2*(1 + Rqsu)*ss*                                              &
          (72*lam(1)**2 - 384*lam(1)*lam(3) + 80*lam(3)**2 -                    &
            9*(12*lam(2)**2 + 3*lam(4)*(-5*lam(5) + lam(6)) +                   &
               39*lam(5)*lam(7) - 15*lam(6)*lam(7) + 16*lam(8)**2 +             &
               2*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7) +                &
                  16*lam(8))) +                                                 &
            24*(6*lam(9)**2 - 6*lam(9)*lam(10) + lam(10)**2))) +                &
      3*expuISs*qdB*quB*qxy2*tanhus**2*                                         &
       (dRqrd*expdISr*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                &
          (576*lam(2)**2 + 9*lam(4)**2 + 882*lam(4)*lam(5) +                    &
            11241*lam(5)**2 - 18*lam(4)*lam(6) - 882*lam(5)*lam(6) +            &
            9*lam(6)**2 - 18*lam(4)*lam(7) - 882*lam(5)*lam(7) +                &
            18*lam(6)*lam(7) + 9*lam(7)**2 - 144*lam(8)**2 +                    &
            72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +                &
               2*lam(8)) + 36*lam(3)*                                           &
             (lam(4) - lam(5) + lam(6) - lam(7) + 10*lam(8)) +                  &
            132*lam(4)*lam(9) + 4740*lam(5)*lam(9) - 132*lam(6)*lam(9) -        &
            132*lam(7)*lam(9) + 412*lam(9)**2 -                                 &
            96*lam(2)*(3*lam(4) + 75*lam(5) -                                   &
               3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +                    &
            24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) + 10*lam(9))*         &
             lam(10) - 144*lam(10)**2) -                                        &
         dRqru*expuISr*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*               &
          (576*lam(2)**2 + 9*lam(4)**2 - 522*lam(4)*lam(5) -                    &
            5499*lam(5)**2 - 18*lam(4)*lam(6) + 522*lam(5)*lam(6) +             &
            9*lam(6)**2 + 90*lam(4)*lam(7) + 198*lam(5)*lam(7) -                &
            90*lam(6)*lam(7) + 117*lam(7)**2 +                                  &
            36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +           &
            288*lam(8)**2 + 72*lam(1)*                                          &
             (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) + 2*lam(8)) -           &
            156*lam(4)*lam(9) - 2100*lam(5)*lam(9) + 156*lam(6)*lam(9) -        &
            204*lam(7)*lam(9) - 56*lam(9)**2 +                                  &
            48*lam(2)*(3*lam(4) + 75*lam(5) -                                   &
               3*(lam(6) + lam(7) + 4*lam(8)) + 16*lam(9)) -                    &
            36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -                   &
               10*lam(9))*lam(10) + 456*lam(10)**2) -                           &
         dRqsu*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                                  &
          (-(expuISr*(1 + Rqru)*                                                &
               (144*lam(2)**2 -                                                 &
                 576*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) -               &
                 192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +               &
                 216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +               &
                 9*(lam(4) + lam(5) + lam(6) + lam(7))**2 -                     &
                 256*(3*lam(1) + lam(3))*lam(10))) +                            &
            expdISr*(1 + Rqrd)*                                                 &
             (144*lam(2)**2 +                                                   &
               216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +             &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                    &
                  10*lam(6)*lam(7) + 13*lam(7)**2 +                             &
                  lam(5)*(-10*lam(6) + 26*lam(7))) +                            &
               288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -             &
               192*lam(8)**2 +                                                  &
               128*lam(3)*(3*lam(4) + 6*lam(5) -                                &
                  3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) -               &
               48*lam(1)*(3*lam(4) + 33*lam(5) -                                &
                  3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))))) -         &
      3*expdISs*qdB*quB*qxy2*tanhds**2*                                         &
       (-(dRqru*expuISr*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*              &
            (576*lam(2)**2 + 9*lam(4)**2 + 882*lam(4)*lam(5) +                  &
              11241*lam(5)**2 - 18*lam(4)*lam(6) - 882*lam(5)*lam(6) +          &
              9*lam(6)**2 - 18*lam(4)*lam(7) - 882*lam(5)*lam(7) +              &
              18*lam(6)*lam(7) + 9*lam(7)**2 - 144*lam(8)**2 +                  &
              72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +              &
                 2*lam(8)) + 36*lam(3)*                                         &
               (lam(4) - lam(5) + lam(6) - lam(7) + 10*lam(8)) +                &
              132*lam(4)*lam(9) + 4740*lam(5)*lam(9) -                          &
              132*lam(6)*lam(9) - 132*lam(7)*lam(9) + 412*lam(9)**2 -           &
              96*lam(2)*(3*lam(4) + 75*lam(5) -                                 &
                 3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +                  &
              24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                   &
                 10*lam(9))*lam(10) - 144*lam(10)**2)) +                        &
         dRqrd*expdISr*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*               &
          (576*lam(2)**2 + 9*lam(4)**2 - 522*lam(4)*lam(5) -                    &
            5499*lam(5)**2 - 18*lam(4)*lam(6) + 522*lam(5)*lam(6) +             &
            9*lam(6)**2 + 90*lam(4)*lam(7) + 198*lam(5)*lam(7) -                &
            90*lam(6)*lam(7) + 117*lam(7)**2 +                                  &
            36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +           &
            288*lam(8)**2 + 72*lam(1)*                                          &
             (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) + 2*lam(8)) -           &
            156*lam(4)*lam(9) - 2100*lam(5)*lam(9) + 156*lam(6)*lam(9) -        &
            204*lam(7)*lam(9) - 56*lam(9)**2 +                                  &
            48*lam(2)*(3*lam(4) + 75*lam(5) -                                   &
               3*(lam(6) + lam(7) + 4*lam(8)) + 16*lam(9)) -                    &
            36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -                   &
               10*lam(9))*lam(10) + 456*lam(10)**2) -                           &
         dRqsd*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                                  &
          (expdISr*(1 + Rqrd)*                                                  &
             (144*lam(2)**2 -                                                   &
               576*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) -                 &
               192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +                 &
               216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                 &
               9*(lam(4) + lam(5) + lam(6) + lam(7))**2 -                       &
               256*(3*lam(1) + lam(3))*lam(10)) -                               &
            expuISr*(1 + Rqru)*                                                 &
             (144*lam(2)**2 +                                                   &
               216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +             &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                    &
                  10*lam(6)*lam(7) + 13*lam(7)**2 +                             &
                  lam(5)*(-10*lam(6) + 26*lam(7))) +                            &
               288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -             &
               192*lam(8)**2 +                                                  &
               128*lam(3)*(3*lam(4) + 6*lam(5) -                                &
                  3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) -               &
               48*lam(1)*(3*lam(4) + 33*lam(5) -                                &
                  3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))))) -         &
      expdISr*qxy2*tanhdr**2*(6*dRqsd*expdISs*quB*qxy2*(1 + Rqrd)*              &
          (1 + Rqsd)**2*tanhds**3*                                              &
          (144*lam(2)**2 - 576*lam(1)*                                          &
             (lam(4) - lam(5) - lam(6) + lam(7)) -                              &
            192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +                    &
            216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                    &
            9*(lam(4) + lam(5) + lam(6) + lam(7))**2 -                          &
            256*(3*lam(1) + lam(3))*lam(10)) -                                  &
         6*dRqsu*expuISs*qdB*qxy2*(1 + Rqrd)*(1 + Rqsu)**2*tanhus**3*           &
          (144*lam(2)**2 + 216*lam(2)*                                          &
             (lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +                          &
            9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                           &
               2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                       &
               10*lam(6)*lam(7) + 13*lam(7)**2 +                                &
               lam(5)*(-10*lam(6) + 26*lam(7))) +                               &
            288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -                &
            192*lam(8)**2 + 128*lam(3)*                                         &
             (3*lam(4) + 6*lam(5) - 3*(lam(6) + 2*lam(7) + 4*lam(9)) +          &
               7*lam(10)) - 48*lam(1)*                                          &
             (3*lam(4) + 33*lam(5) -                                            &
               3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))) -              &
         3*expdISs*qdB*quB*tanhds**2*                                           &
          (-(dRqsd*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                   &
               (144*lam(2)**2 -                                                 &
                 576*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) -               &
                 192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +               &
                 216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +               &
                 9*(lam(4) + lam(5) + lam(6) + lam(7))**2 -                     &
                 256*(3*lam(1) + lam(3))*lam(10))) +                            &
            dRqrd*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                    &
             (576*lam(2)**2 + 9*lam(4)**2 - 522*lam(4)*lam(5) -                 &
               5499*lam(5)**2 - 18*lam(4)*lam(6) + 522*lam(5)*lam(6) +          &
               9*lam(6)**2 + 90*lam(4)*lam(7) + 198*lam(5)*lam(7) -             &
               90*lam(6)*lam(7) + 117*lam(7)**2 +                               &
               36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +        &
               288*lam(8)**2 +                                                  &
               72*lam(1)*(3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +           &
                  2*lam(8)) - 156*lam(4)*lam(9) - 2100*lam(5)*lam(9) +          &
               156*lam(6)*lam(9) - 204*lam(7)*lam(9) - 56*lam(9)**2 +           &
               48*lam(2)*(3*lam(4) + 75*lam(5) -                                &
                  3*(lam(6) + lam(7) + 4*lam(8)) + 16*lam(9)) -                 &
               36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -                &
                  10*lam(9))*lam(10) + 456*lam(10)**2)) +                       &
         3*expuISs*qdB*quB*tanhus**2*                                           &
          (dRqrd*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                     &
             (576*lam(2)**2 + 9*lam(4)**2 + 882*lam(4)*lam(5) +                 &
               11241*lam(5)**2 - 18*lam(4)*lam(6) - 882*lam(5)*lam(6) +         &
               9*lam(6)**2 - 18*lam(4)*lam(7) - 882*lam(5)*lam(7) +             &
               18*lam(6)*lam(7) + 9*lam(7)**2 - 144*lam(8)**2 +                 &
               72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +             &
                  2*lam(8)) +                                                   &
               36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                   &
                  10*lam(8)) + 132*lam(4)*lam(9) + 4740*lam(5)*lam(9) -         &
               132*lam(6)*lam(9) - 132*lam(7)*lam(9) + 412*lam(9)**2 -          &
               96*lam(2)*(3*lam(4) + 75*lam(5) -                                &
                  3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +                 &
               24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                  &
                  10*lam(9))*lam(10) - 144*lam(10)**2) -                        &
            dRqsu*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                    &
             (144*lam(2)**2 +                                                   &
               216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +             &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                    &
                  10*lam(6)*lam(7) + 13*lam(7)**2 +                             &
                  lam(5)*(-10*lam(6) + 26*lam(7))) +                            &
               288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -             &
               192*lam(8)**2 +                                                  &
               128*lam(3)*(3*lam(4) + 6*lam(5) -                                &
                  3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) -               &
               48*lam(1)*(3*lam(4) + 33*lam(5) -                                &
                  3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10)))) +          &
         2*expdISs*quB*(1 + Rqrd)*tanhds*                                       &
          (-3*dRqsd*qxy2*(1 + Rqsd)**2*                                         &
             (144*lam(2)**2 -                                                   &
               576*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) -                 &
               192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +                 &
               216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                 &
               9*(lam(4) + lam(5) + lam(6) + lam(7))**2 -                       &
               256*(3*lam(1) + lam(3))*lam(10)) +                               &
            dRqrd*(-3*qzt2*(1 + Rqrd)*(1 + Rqsd)*                               &
                (576*lam(2)**2 + 9*lam(4)**2 - 522*lam(4)*lam(5) -              &
                  5499*lam(5)**2 - 18*lam(4)*lam(6) +                           &
                  522*lam(5)*lam(6) + 9*lam(6)**2 + 90*lam(4)*lam(7) +          &
                  198*lam(5)*lam(7) - 90*lam(6)*lam(7) +                        &
                  117*lam(7)**2 -                                               &
                  36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) +                  &
                     23*lam(7)) + 288*lam(8)**2 -                               &
                  72*lam(1)*(3*                                                 &
                      (lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +                 &
                     2*lam(8)) - 156*lam(4)*lam(9) -                            &
                  2100*lam(5)*lam(9) + 156*lam(6)*lam(9) -                      &
                  204*lam(7)*lam(9) - 56*lam(9)**2 +                            &
                  48*lam(2)*(3*lam(4) + 75*lam(5) -                             &
                     3*(lam(6) + lam(7) + 4*lam(8)) + 16*lam(9)) -              &
                  36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -             &
                     10*lam(9))*lam(10) + 456*lam(10)**2) +                     &
               2*m2*(864*lam(2)**2 + 45*lam(4)**2 + 8793*lam(5)**2 -            &
                  414*lam(5)*lam(6) + 45*lam(6)**2 -                            &
                  1386*lam(5)*lam(7) - 234*lam(6)*lam(7) +                      &
                  369*lam(7)**2 + 540*lam(5)*lam(8) +                           &
                  108*lam(6)*lam(8) - 540*lam(7)*lam(8) +                       &
                  3648*lam(5)*lam(9) - 192*lam(6)*lam(9) -                      &
                  192*lam(7)*lam(9) + 720*lam(8)*lam(9) +                       &
                  608*lam(9)**2 +                                               &
                  6*lam(4)*(69*lam(5) - 15*lam(6) + 39*lam(7) -                 &
                     18*lam(8) + 32*lam(9)) +                                   &
                  72*lam(2)*(3*lam(4) - 87*lam(5) - 3*lam(6) +                  &
                     15*lam(7) - 26*lam(9) - 30*lam(10)) +                      &
                  72*(9*lam(8) - 8*lam(9))*lam(10) + 576*lam(10)**2))) +        &
         2*expuISs*(1 + Rqrd)*tanhus*                                           &
          (3*dRqsu*qdB*qxy2*(1 + Rqsu)**2*                                      &
             (144*lam(2)**2 +                                                   &
               216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +             &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                    &
                  10*lam(6)*lam(7) + 13*lam(7)**2 +                             &
                  lam(5)*(-10*lam(6) + 26*lam(7))) +                            &
               288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -             &
               192*lam(8)**2 +                                                  &
               128*lam(3)*(3*lam(4) + 6*lam(5) -                                &
                  3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) -               &
               48*lam(1)*(3*lam(4) + 33*lam(5) -                                &
                  3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))) +           &
            dRqrd*quB*(3*qzt2*(1 + Rqrd)*(1 + Rqsu)*                            &
                (576*lam(2)**2 + 9*lam(4)**2 + 882*lam(4)*lam(5) +              &
                  11241*lam(5)**2 - 18*lam(4)*lam(6) -                          &
                  882*lam(5)*lam(6) + 9*lam(6)**2 - 18*lam(4)*lam(7) -          &
                  882*lam(5)*lam(7) + 18*lam(6)*lam(7) + 9*lam(7)**2 -          &
                  144*lam(8)**2 -                                               &
                  72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +          &
                     2*lam(8)) -                                                &
                  36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                &
                     10*lam(8)) + 132*lam(4)*lam(9) +                           &
                  4740*lam(5)*lam(9) - 132*lam(6)*lam(9) -                      &
                  132*lam(7)*lam(9) + 412*lam(9)**2 -                           &
                  96*lam(2)*(3*lam(4) + 75*lam(5) -                             &
                     3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +              &
                  24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +               &
                     10*lam(9))*lam(10) - 144*lam(10)**2) +                     &
               2*m2*(1728*lam(2)**2 + 9*lam(4)**2 + 16857*lam(5)**2 -           &
                  1314*lam(5)*lam(6) + 9*lam(6)**2 -                            &
                  1314*lam(5)*lam(7) + 18*lam(6)*lam(7) + 9*lam(7)**2 +         &
                  108*lam(5)*lam(8) - 108*lam(6)*lam(8) -                       &
                  108*lam(7)*lam(8) + 7296*lam(5)*lam(9) -                      &
                  384*lam(6)*lam(9) - 384*lam(7)*lam(9) +                       &
                  360*lam(8)*lam(9) + 928*lam(9)**2 +                           &
                  6*lam(4)*(219*lam(5) -                                        &
                     3*(lam(6) + lam(7) - 6*lam(8)) + 64*lam(9)) -              &
                  432*lam(8)*lam(10) -                                          &
                  72*lam(2)*(3*lam(4) + 147*lam(5) - 3*lam(6) -                 &
                     3*lam(7) + 22*lam(9) + 12*lam(10))))) +                    &
         3*qdB*quB*(dRqrd*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                       &
             (-(expuISs*(1 + Rqsu)*                                             &
                  (576*lam(2)**2 + 9*lam(4)**2 + 882*lam(4)*lam(5) +            &
                    11241*lam(5)**2 - 18*lam(4)*lam(6) -                        &
                    882*lam(5)*lam(6) + 9*lam(6)**2 -                           &
                    18*lam(4)*lam(7) - 882*lam(5)*lam(7) +                      &
                    18*lam(6)*lam(7) + 9*lam(7)**2 - 144*lam(8)**2 +            &
                    72*lam(1)*                                                  &
                     (3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +                 &
                       2*lam(8)) +                                              &
                    36*lam(3)*                                                  &
                     (lam(4) - lam(5) + lam(6) - lam(7) + 10*lam(8)) +          &
                    132*lam(4)*lam(9) + 4740*lam(5)*lam(9) -                    &
                    132*lam(6)*lam(9) - 132*lam(7)*lam(9) +                     &
                    412*lam(9)**2 -                                             &
                    96*lam(2)*                                                  &
                     (3*lam(4) + 75*lam(5) -                                    &
                       3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +            &
                    24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +             &
                       10*lam(9))*lam(10) - 144*lam(10)**2)) +                  &
               expdISs*(1 + Rqsd)*                                              &
                (576*lam(2)**2 + 9*lam(4)**2 - 522*lam(4)*lam(5) -              &
                  5499*lam(5)**2 - 18*lam(4)*lam(6) +                           &
                  522*lam(5)*lam(6) + 9*lam(6)**2 + 90*lam(4)*lam(7) +          &
                  198*lam(5)*lam(7) - 90*lam(6)*lam(7) + 117*lam(7)**2 +        &
                  36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) +                  &
                     23*lam(7)) + 288*lam(8)**2 +                               &
                  72*lam(1)*(3*                                                 &
                      (lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) + 2*lam(8)        &
) - 156*lam(4)*lam(9) - 2100*lam(5)*lam(9) + 156*lam(6)*lam(9) -                &
                  204*lam(7)*lam(9) - 56*lam(9)**2 +                            &
                  48*lam(2)*(3*lam(4) + 75*lam(5) -                             &
                     3*(lam(6) + lam(7) + 4*lam(8)) + 16*lam(9)) -              &
                  36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -             &
                     10*lam(9))*lam(10) + 456*lam(10)**2)) -                    &
            (1 + Rqrd)*(dRqsd*expdISs*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*           &
                (144*lam(2)**2 -                                                &
                  576*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) -              &
                  192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +              &
                  216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +              &
                  9*(lam(4) + lam(5) + lam(6) + lam(7))**2 -                    &
                  256*(3*lam(1) + lam(3))*lam(10)) -                            &
               dRqsu*expuISs*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                    &
                (144*lam(2)**2 +                                                &
                  216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +          &
                  9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                     &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                 &
                     10*lam(6)*lam(7) + 13*lam(7)**2 +                          &
                     lam(5)*(-10*lam(6) + 26*lam(7))) +                         &
                  288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -          &
                  192*lam(8)**2 +                                               &
                  128*lam(3)*(3*lam(4) + 6*lam(5) -                             &
                     3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) -            &
                  48*lam(1)*(3*lam(4) + 33*lam(5) -                             &
                     3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))))))       &
- expdISr*tanhdr*(-6*dRqrd*quB*qxy2*(1 + Rqrd)*                                 &
          (-(expuISs*(qxy2 + qzt2)*(1 + Rqrd)*(1 + Rqsu)*                       &
               (576*lam(2)**2 + 9*lam(4)**2 + 882*lam(4)*lam(5) +               &
                 11241*lam(5)**2 - 18*lam(4)*lam(6) -                           &
                 882*lam(5)*lam(6) + 9*lam(6)**2 - 18*lam(4)*lam(7) -           &
                 882*lam(5)*lam(7) + 18*lam(6)*lam(7) + 9*lam(7)**2 -           &
                 144*lam(8)**2 +                                                &
                 72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +           &
                    2*lam(8)) +                                                 &
                 36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                 &
                    10*lam(8)) + 132*lam(4)*lam(9) +                            &
                 4740*lam(5)*lam(9) - 132*lam(6)*lam(9) -                       &
                 132*lam(7)*lam(9) + 412*lam(9)**2 -                            &
                 96*lam(2)*(3*lam(4) + 75*lam(5) -                              &
                    3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +               &
                 24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                &
                    10*lam(9))*lam(10) - 144*lam(10)**2)) -                     &
            24*expdISs*m2*(24*lam(2)**2 +                                       &
               6*(6*lam(1)**2 + 14*lam(1)*lam(3) + 5*lam(3)**2) +               &
               234*lam(5)**2 - 18*lam(5)*lam(6) +                               &
               9*lam(4)*(7*lam(5) + lam(6)) + 99*lam(5)*lam(7) +                &
               45*lam(6)*lam(7) + 30*lam(5)*lam(8) + 6*lam(6)*lam(8) -          &
               6*lam(8)**2 + 96*lam(5)*lam(9) + 20*lam(8)*lam(9) -              &
               8*lam(9)**2 +                                                    &
               2*lam(2)*(3*lam(4) - 87*lam(5) - 3*lam(6) + 15*lam(7) -          &
                  26*lam(9) - 30*lam(10)) +                                     &
               6*(3*lam(8) + 8*lam(9))*lam(10) - 48*lam(10)**2) +               &
            expdISs*(qxy2 + qzt2)*(1 + Rqrd)*(1 + Rqsd)*                        &
             (576*lam(2)**2 + 9*lam(4)**2 - 522*lam(4)*lam(5) -                 &
               5499*lam(5)**2 - 18*lam(4)*lam(6) + 522*lam(5)*lam(6) +          &
               9*lam(6)**2 + 90*lam(4)*lam(7) + 198*lam(5)*lam(7) -             &
               90*lam(6)*lam(7) + 117*lam(7)**2 +                               &
               36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) +                     &
                  23*lam(7)) + 288*lam(8)**2 +                                  &
               72*lam(1)*(3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +           &
                  2*lam(8)) - 156*lam(4)*lam(9) - 2100*lam(5)*lam(9) +          &
               156*lam(6)*lam(9) - 204*lam(7)*lam(9) - 56*lam(9)**2 +           &
               48*lam(2)*(3*lam(4) + 75*lam(5) -                                &
                  3*(lam(6) + lam(7) + 4*lam(8)) + 16*lam(9)) -                 &
               36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -                &
                  10*lam(9))*lam(10) + 456*lam(10)**2) +                        &
            24*expuISs*m2*(-48*lam(2)**2 +                                      &
               3*(12*lam(1)**2 + 4*lam(1)*lam(3) + 5*lam(3)**2) -               &
               468*lam(5)**2 + 36*lam(5)*lam(6) +                               &
               9*lam(4)*(-5*lam(5) + lam(6)) + 45*lam(5)*lam(7) -               &
               9*lam(6)*lam(7) - 6*lam(5)*lam(8) + 6*lam(6)*lam(8) +            &
               15*lam(8)**2 - 192*lam(5)*lam(9) - 10*lam(8)*lam(9) -            &
               8*lam(9)**2 + 12*lam(8)*lam(10) +                                &
               lam(2)*(6*lam(4) + 294*lam(5) - 6*lam(6) - 6*lam(7) +            &
                  44*lam(9) + 24*lam(10)))) +                                   &
         2*expdISs*quB*qxy2*(1 + Rqsd)*tanhds**2*                               &
          (3*dRqrd*qxy2*(1 + Rqrd)**2*                                          &
             (576*lam(2)**2 + 9*lam(4)**2 - 522*lam(4)*lam(5) -                 &
               5499*lam(5)**2 - 18*lam(4)*lam(6) + 522*lam(5)*lam(6) +          &
               9*lam(6)**2 + 90*lam(4)*lam(7) + 198*lam(5)*lam(7) -             &
               90*lam(6)*lam(7) + 117*lam(7)**2 +                               &
               36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) +                     &
                  23*lam(7)) + 288*lam(8)**2 +                                  &
               72*lam(1)*(3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +           &
                  2*lam(8)) - 156*lam(4)*lam(9) - 2100*lam(5)*lam(9) +          &
               156*lam(6)*lam(9) - 204*lam(7)*lam(9) - 56*lam(9)**2 +           &
               48*lam(2)*(3*lam(4) + 75*lam(5) -                                &
                  3*(lam(6) + lam(7) + 4*lam(8)) + 16*lam(9)) -                 &
               36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -                &
                  10*lam(9))*lam(10) + 456*lam(10)**2) +                        &
            dRqsd*(6*m2*(128*(3*lam(1) + lam(3))**2 +                           &
                  15*(lam(4) - lam(5) - lam(6) + lam(7))**2) +                  &
               384*m2*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +             &
               64*m2*lam(10)**2 +                                               &
               3*qzt2*(1 + Rqrd)*(1 + Rqsd)*                                    &
                (144*lam(2)**2 +                                                &
                  576*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) +              &
                  192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +              &
                  216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +              &
                  9*(lam(4) + lam(5) + lam(6) + lam(7))**2 +                    &
                  256*(3*lam(1) + lam(3))*lam(10)))) -                          &
         expuISs*qdB*quB*qzt2*tanhus*                                           &
          (1728*dRqrd*lam(2)**2 -                                               &
            3*dRqrd*Rqsu*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                        &
             (576*lam(2)**2 + 9*lam(4)**2 + 882*lam(4)*lam(5) +                 &
               11241*lam(5)**2 - 18*lam(4)*lam(6) - 882*lam(5)*lam(6) +         &
               9*lam(6)**2 - 18*lam(4)*lam(7) - 882*lam(5)*lam(7) +             &
               18*lam(6)*lam(7) + 9*lam(7)**2 - 144*lam(8)**2 -                 &
               72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +             &
                  2*lam(8)) -                                                   &
               36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                   &
                  10*lam(8)) + 132*lam(4)*lam(9) + 4740*lam(5)*lam(9) -         &
               132*lam(6)*lam(9) - 132*lam(7)*lam(9) + 412*lam(9)**2 -          &
               96*lam(2)*(3*lam(4) + 75*lam(5) -                                &
                  3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +                 &
               24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                  &
                  10*lam(9))*lam(10) - 144*lam(10)**2) +                        &
            4*dRqsu*m2*(1 + Rqsu)*ss*                                           &
             (96*(9*lam(1)**2 - 48*lam(1)*lam(3) + 10*lam(3)**2) -              &
               9*(lam(4)**2 - 11*lam(5)**2 + 14*lam(5)*lam(6) +                 &
                  lam(6)**2 - 2*lam(4)*(7*lam(5) + lam(6))) -                   &
               18*(7*lam(4) + 11*lam(5) - 7*lam(6))*lam(7) +                    &
               99*lam(7)**2 -                                                   &
               576*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +             &
               960*lam(9)**2 +                                                  &
               192*(lam(4) + 5*lam(5) - lam(6) - 5*(lam(7) + lam(9)))*          &
                lam(10) + 224*lam(10)**2) -                                     &
            4*dRqrd*m2*(1 + Rqrd)*sr*                                           &
             (1728*lam(2)**2 + 9*lam(4)**2 + 16857*lam(5)**2 -                  &
               1314*lam(5)*lam(6) + 9*lam(6)**2 - 1314*lam(5)*lam(7) +          &
               18*lam(6)*lam(7) + 9*lam(7)**2 + 108*lam(5)*lam(8) -             &
               108*lam(6)*lam(8) - 108*lam(7)*lam(8) +                          &
               7296*lam(5)*lam(9) - 384*lam(6)*lam(9) -                         &
               384*lam(7)*lam(9) + 360*lam(8)*lam(9) + 928*lam(9)**2 +          &
               6*lam(4)*(219*lam(5) - 3*(lam(6) + lam(7) - 6*lam(8)) +          &
                  64*lam(9)) - 432*lam(8)*lam(10) -                             &
               72*lam(2)*(3*lam(4) + 147*lam(5) - 3*lam(6) - 3*lam(7) +         &
                  22*lam(9) + 12*lam(10))) +                                    &
            3*dRqsu*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                  &
             (144*lam(2)**2 +                                                   &
               216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +             &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                    &
                  10*lam(6)*lam(7) + 13*lam(7)**2 +                             &
                  lam(5)*(-10*lam(6) + 26*lam(7))) +                            &
               288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -             &
               192*lam(8)**2 -                                                  &
               128*lam(3)*(3*lam(4) + 6*lam(5) -                                &
                  3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) +               &
               48*lam(1)*(3*lam(4) + 33*lam(5) -                                &
                  3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))) +           &
            3*dRqrd*(-216*lam(1)*lam(4) - 36*lam(3)*lam(4) +                    &
               9*lam(4)**2 + 36*lam(3)*lam(5) + 882*lam(4)*lam(5) +             &
               11241*lam(5)**2 - 36*lam(3)*lam(6) - 18*lam(4)*lam(6) -          &
               882*lam(5)*lam(6) + 9*lam(6)**2 + 36*lam(3)*lam(7) -             &
               18*lam(4)*lam(7) - 882*lam(5)*lam(7) + 18*lam(6)*lam(7) +        &
               9*lam(7)**2 + 72*lam(1)*                                         &
                (3*(lam(5) - lam(6) + lam(7)) - 2*lam(8)) -                     &
               360*lam(3)*lam(8) - 144*lam(8)**2 + 132*lam(4)*lam(9) +          &
               4740*lam(5)*lam(9) - 132*lam(6)*lam(9) -                         &
               132*lam(7)*lam(9) + 412*lam(9)**2 -                              &
               96*lam(2)*(3*lam(4) + 75*lam(5) -                                &
                  3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +                 &
               24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                  &
                  10*lam(9))*lam(10) - 144*lam(10)**2 -                         &
               2*qzt2*(1 + Rqrd)**2*sr*                                         &
                (576*lam(2)**2 + 9*lam(4)**2 + 882*lam(4)*lam(5) +              &
                  11241*lam(5)**2 - 18*lam(4)*lam(6) -                          &
                  882*lam(5)*lam(6) + 9*lam(6)**2 - 18*lam(4)*lam(7) -          &
                  882*lam(5)*lam(7) + 18*lam(6)*lam(7) + 9*lam(7)**2 -          &
                  144*lam(8)**2 -                                               &
                  72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +          &
                     2*lam(8)) -                                                &
                  36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                &
                     10*lam(8)) + 132*lam(4)*lam(9) +                           &
                  4740*lam(5)*lam(9) - 132*lam(6)*lam(9) -                      &
                  132*lam(7)*lam(9) + 412*lam(9)**2 -                           &
                  96*lam(2)*(3*lam(4) + 75*lam(5) -                             &
                     3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +              &
                  24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +               &
                     10*lam(9))*lam(10) - 144*lam(10)**2))) +                   &
         expdISs*qdB*quB*qzt2*tanhds*                                           &
          (1728*dRqrd*lam(2)**2 +                                               &
            3*dRqsd*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                  &
             (3*(48*lam(2)**2 +                                                 &
                  192*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) +              &
                  64*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +               &
                  72*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +               &
                  3*(lam(4) + lam(5) + lam(6) + lam(7))**2) +                   &
               256*(3*lam(1) + lam(3))*lam(10)) +                               &
            4*dRqsd*m2*(1 + Rqsd)*ss*                                           &
             (384*(3*lam(1) + lam(3))**2 +                                      &
               45*(lam(4) - lam(5) - lam(6) + lam(7))**2 +                      &
               192*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +                &
               32*lam(10)**2) -                                                 &
            3*dRqrd*Rqsd*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                        &
             (576*lam(2)**2 + 9*lam(4)**2 - 522*lam(4)*lam(5) -                 &
               5499*lam(5)**2 - 18*lam(4)*lam(6) + 522*lam(5)*lam(6) +          &
               9*lam(6)**2 + 90*lam(4)*lam(7) + 198*lam(5)*lam(7) -             &
               90*lam(6)*lam(7) + 117*lam(7)**2 -                               &
               36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) +                     &
                  23*lam(7)) + 288*lam(8)**2 -                                  &
               72*lam(1)*(3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +           &
                  2*lam(8)) - 156*lam(4)*lam(9) - 2100*lam(5)*lam(9) +          &
               156*lam(6)*lam(9) - 204*lam(7)*lam(9) - 56*lam(9)**2 +           &
               48*lam(2)*(3*lam(4) + 75*lam(5) -                                &
                  3*(lam(6) + lam(7) + 4*lam(8)) + 16*lam(9)) -                 &
               36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -                &
                  10*lam(9))*lam(10) + 456*lam(10)**2) +                        &
            4*dRqrd*m2*(1 + Rqrd)*sr*                                           &
             (864*lam(2)**2 + 45*lam(4)**2 + 8793*lam(5)**2 -                   &
               414*lam(5)*lam(6) + 45*lam(6)**2 - 1386*lam(5)*lam(7) -          &
               234*lam(6)*lam(7) + 369*lam(7)**2 + 540*lam(5)*lam(8) +          &
               108*lam(6)*lam(8) - 540*lam(7)*lam(8) +                          &
               3648*lam(5)*lam(9) - 192*lam(6)*lam(9) -                         &
               192*lam(7)*lam(9) + 720*lam(8)*lam(9) + 608*lam(9)**2 +          &
               6*lam(4)*(69*lam(5) - 15*lam(6) + 39*lam(7) -                    &
                  18*lam(8) + 32*lam(9)) +                                      &
               72*lam(2)*(3*lam(4) - 87*lam(5) - 3*lam(6) +                     &
                  15*lam(7) - 26*lam(9) - 30*lam(10)) +                         &
               72*(9*lam(8) - 8*lam(9))*lam(10) + 576*lam(10)**2) +             &
            3*dRqrd*(-216*lam(1)*lam(4) - 252*lam(3)*lam(4) +                   &
               9*lam(4)**2 - 828*lam(3)*lam(5) - 522*lam(4)*lam(5) -            &
               5499*lam(5)**2 - 252*lam(3)*lam(6) - 18*lam(4)*lam(6) +          &
               522*lam(5)*lam(6) + 9*lam(6)**2 - 828*lam(3)*lam(7) +            &
               90*lam(4)*lam(7) + 198*lam(5)*lam(7) - 90*lam(6)*lam(7) +        &
               117*lam(7)**2 + 288*lam(8)**2 -                                  &
               72*lam(1)*(3*(5*lam(5) + lam(6) + 5*lam(7)) + 2*lam(8)) -        &
               156*lam(4)*lam(9) - 2100*lam(5)*lam(9) +                         &
               156*lam(6)*lam(9) - 204*lam(7)*lam(9) - 56*lam(9)**2 +           &
               48*lam(2)*(3*lam(4) + 75*lam(5) -                                &
                  3*(lam(6) + lam(7) + 4*lam(8)) + 16*lam(9)) -                 &
               36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -                &
                  10*lam(9))*lam(10) + 456*lam(10)**2 -                         &
               2*qzt2*(1 + Rqrd)**2*sr*                                         &
                (576*lam(2)**2 + 9*lam(4)**2 - 522*lam(4)*lam(5) -              &
                  5499*lam(5)**2 - 18*lam(4)*lam(6) +                           &
                  522*lam(5)*lam(6) + 9*lam(6)**2 + 90*lam(4)*lam(7) +          &
                  198*lam(5)*lam(7) - 90*lam(6)*lam(7) + 117*lam(7)**2 -        &
                  36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) +                  &
                     23*lam(7)) + 288*lam(8)**2 -                               &
                  72*lam(1)*(3*                                                 &
                      (lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) + 2*lam(8)        &
) - 156*lam(4)*lam(9) - 2100*lam(5)*lam(9) + 156*lam(6)*lam(9) -                &
                  204*lam(7)*lam(9) - 56*lam(9)**2 +                            &
                  48*lam(2)*(3*lam(4) + 75*lam(5) -                             &
                     3*(lam(6) + lam(7) + 4*lam(8)) + 16*lam(9)) -              &
                  36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -             &
                     10*lam(9))*lam(10) + 456*lam(10)**2))) -                   &
         2*expuISs*qxy2*(1 + Rqsu)*tanhus**2*                                   &
          (3*dRqrd*quB*qxy2*(1 + Rqrd)**2*                                      &
             (576*lam(2)**2 + 9*lam(4)**2 + 882*lam(4)*lam(5) +                 &
               11241*lam(5)**2 - 18*lam(4)*lam(6) - 882*lam(5)*lam(6) +         &
               9*lam(6)**2 - 18*lam(4)*lam(7) - 882*lam(5)*lam(7) +             &
               18*lam(6)*lam(7) + 9*lam(7)**2 - 144*lam(8)**2 +                 &
               72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +             &
                  2*lam(8)) +                                                   &
               36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                   &
                  10*lam(8)) + 132*lam(4)*lam(9) + 4740*lam(5)*lam(9) -         &
               132*lam(6)*lam(9) - 132*lam(7)*lam(9) + 412*lam(9)**2 -          &
               96*lam(2)*(3*lam(4) + 75*lam(5) -                                &
                  3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +                 &
               24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                  &
                  10*lam(9))*lam(10) - 144*lam(10)**2) +                        &
            dRqsu*qdB*(6*m2*(32*                                                &
                   (9*lam(1)**2 - 48*lam(1)*lam(3) + 10*lam(3)**2) -            &
                  3*lam(4)**2 + 33*lam(5)**2 - 3*lam(6)**2 +                    &
                  6*lam(4)*(7*lam(5) + lam(6) - 7*lam(7)) +                     &
                  42*lam(6)*lam(7) + 33*lam(7)**2 -                             &
                  6*lam(5)*(7*lam(6) + 11*lam(7)) -                             &
                  192*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +          &
                  320*lam(9)**2) +                                              &
               384*m2*(lam(4) + 5*lam(5) - lam(6) - 5*(lam(7) + lam(9)))*       &
                lam(10) + 448*m2*lam(10)**2 +                                   &
               3*qzt2*(1 + Rqrd)*(1 + Rqsu)*                                    &
                (144*lam(2)**2 +                                                &
                  216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +          &
                  9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                     &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                 &
                     10*lam(6)*lam(7) + 13*lam(7)**2 +                          &
                     lam(5)*(-10*lam(6) + 26*lam(7))) +                         &
                  288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -          &
                  192*lam(8)**2 -                                               &
                  128*lam(3)*(3*lam(4) + 6*lam(5) -                             &
                     3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) +            &
                  48*lam(1)*(3*lam(4) + 33*lam(5) -                             &
                     3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))))))       &
- expuISr*qxy2*tanhur**2*(-6*dRqsd*expdISs*quB*qxy2*(1 + Rqru)*                 &
          (1 + Rqsd)**2*tanhds**3*                                              &
          (144*lam(2)**2 + 216*lam(2)*                                          &
             (lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +                          &
            9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                           &
               2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                       &
               10*lam(6)*lam(7) + 13*lam(7)**2 +                                &
               lam(5)*(-10*lam(6) + 26*lam(7))) +                               &
            288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -                &
            192*lam(8)**2 + 128*lam(3)*                                         &
             (3*lam(4) + 6*lam(5) - 3*(lam(6) + 2*lam(7) + 4*lam(9)) +          &
               7*lam(10)) - 48*lam(1)*                                          &
             (3*lam(4) + 33*lam(5) -                                            &
               3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))) +              &
         3*expdISs*qdB*quB*tanhds**2*                                           &
          (dRqru*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                     &
             (576*lam(2)**2 + 9*lam(4)**2 + 882*lam(4)*lam(5) +                 &
               11241*lam(5)**2 - 18*lam(4)*lam(6) - 882*lam(5)*lam(6) +         &
               9*lam(6)**2 - 18*lam(4)*lam(7) - 882*lam(5)*lam(7) +             &
               18*lam(6)*lam(7) + 9*lam(7)**2 - 144*lam(8)**2 +                 &
               72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +             &
                  2*lam(8)) +                                                   &
               36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                   &
                  10*lam(8)) + 132*lam(4)*lam(9) + 4740*lam(5)*lam(9) -         &
               132*lam(6)*lam(9) - 132*lam(7)*lam(9) + 412*lam(9)**2 -          &
               96*lam(2)*(3*lam(4) + 75*lam(5) -                                &
                  3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +                 &
               24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                  &
                  10*lam(9))*lam(10) - 144*lam(10)**2) -                        &
            dRqsd*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                    &
             (144*lam(2)**2 +                                                   &
               216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +             &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                    &
                  10*lam(6)*lam(7) + 13*lam(7)**2 +                             &
                  lam(5)*(-10*lam(6) + 26*lam(7))) +                            &
               288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -             &
               192*lam(8)**2 +                                                  &
               128*lam(3)*(3*lam(4) + 6*lam(5) -                                &
                  3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) -               &
               48*lam(1)*(3*lam(4) + 33*lam(5) -                                &
                  3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10)))) +          &
         2*expdISs*(1 + Rqru)*tanhds*                                           &
          (3*dRqsd*quB*qxy2*(1 + Rqsd)**2*                                      &
             (144*lam(2)**2 +                                                   &
               216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +             &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                    &
                  10*lam(6)*lam(7) + 13*lam(7)**2 +                             &
                  lam(5)*(-10*lam(6) + 26*lam(7))) +                            &
               288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -             &
               192*lam(8)**2 +                                                  &
               128*lam(3)*(3*lam(4) + 6*lam(5) -                                &
                  3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) -               &
               48*lam(1)*(3*lam(4) + 33*lam(5) -                                &
                  3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))) +           &
            dRqru*qdB*(3*qzt2*(1 + Rqru)*(1 + Rqsd)*                            &
                (576*lam(2)**2 + 9*lam(4)**2 + 882*lam(4)*lam(5) +              &
                  11241*lam(5)**2 - 18*lam(4)*lam(6) -                          &
                  882*lam(5)*lam(6) + 9*lam(6)**2 - 18*lam(4)*lam(7) -          &
                  882*lam(5)*lam(7) + 18*lam(6)*lam(7) + 9*lam(7)**2 -          &
                  144*lam(8)**2 -                                               &
                  72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +          &
                     2*lam(8)) -                                                &
                  36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                &
                     10*lam(8)) + 132*lam(4)*lam(9) +                           &
                  4740*lam(5)*lam(9) - 132*lam(6)*lam(9) -                      &
                  132*lam(7)*lam(9) + 412*lam(9)**2 -                           &
                  96*lam(2)*(3*lam(4) + 75*lam(5) -                             &
                     3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +              &
                  24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +               &
                     10*lam(9))*lam(10) - 144*lam(10)**2) +                     &
               2*m2*(1728*lam(2)**2 + 9*lam(4)**2 + 16857*lam(5)**2 -           &
                  1314*lam(5)*lam(6) + 9*lam(6)**2 -                            &
                  1314*lam(5)*lam(7) + 18*lam(6)*lam(7) + 9*lam(7)**2 +         &
                  108*lam(5)*lam(8) - 108*lam(6)*lam(8) -                       &
                  108*lam(7)*lam(8) + 7296*lam(5)*lam(9) -                      &
                  384*lam(6)*lam(9) - 384*lam(7)*lam(9) +                       &
                  360*lam(8)*lam(9) + 928*lam(9)**2 +                           &
                  6*lam(4)*(219*lam(5) -                                        &
                     3*(lam(6) + lam(7) - 6*lam(8)) + 64*lam(9)) -              &
                  432*lam(8)*lam(10) -                                          &
                  72*lam(2)*(3*lam(4) + 147*lam(5) - 3*lam(6) -                 &
                     3*lam(7) + 22*lam(9) + 12*lam(10))))) +                    &
         qdB*(3*(1 + Rqru)*(dRqsu*expuISs*(-1 + tanhus)*(1 + tanhus)*           &
                (quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +                           &
                  2*qxy2*(1 + Rqsu)**2*tanhus)*                                 &
                (3*(48*lam(2)**2 -                                              &
                     192*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) -           &
                     64*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +            &
                     72*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +            &
                     3*(lam(4) + lam(5) + lam(6) + lam(7))**2) -                &
                  256*(3*lam(1) + lam(3))*lam(10)) +                            &
               dRqsd*expdISs*quB*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                &
                (144*lam(2)**2 +                                                &
                  216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +          &
                  9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                     &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                 &
                     10*lam(6)*lam(7) + 13*lam(7)**2 +                          &
                     lam(5)*(-10*lam(6) + 26*lam(7))) +                         &
                  288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -          &
                  192*lam(8)**2 +                                               &
                  128*lam(3)*                                                   &
                   (3*lam(4) + 6*lam(5) -                                       &
                     3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) -            &
                  48*lam(1)*(3*lam(4) + 33*lam(5) -                             &
                     3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10)))) +       &
            dRqru*(-3*expdISs*quB*(1 + Rqsd)*                                   &
                (-1 + 2*qzt2*(1 + Rqru)**2*sr)*                                 &
                (576*lam(2)**2 + 9*lam(4)**2 + 882*lam(4)*lam(5) +              &
                  11241*lam(5)**2 - 18*lam(4)*lam(6) -                          &
                  882*lam(5)*lam(6) + 9*lam(6)**2 - 18*lam(4)*lam(7) -          &
                  882*lam(5)*lam(7) + 18*lam(6)*lam(7) + 9*lam(7)**2 -          &
                  144*lam(8)**2 +                                               &
                  72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +          &
                     2*lam(8)) +                                                &
                  36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                &
                     10*lam(8)) + 132*lam(4)*lam(9) +                           &
                  4740*lam(5)*lam(9) - 132*lam(6)*lam(9) -                      &
                  132*lam(7)*lam(9) + 412*lam(9)**2 -                           &
                  96*lam(2)*(3*lam(4) + 75*lam(5) -                             &
                     3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +              &
                  24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +               &
                     10*lam(9))*lam(10) - 144*lam(10)**2) +                     &
               expuISs*(-3*quB*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*       &
                   (-1 + tanhus**2)*                                            &
                   (576*lam(2)**2 + 9*lam(4)**2 - 522*lam(4)*lam(5) -           &
                     5499*lam(5)**2 - 18*lam(4)*lam(6) +                        &
                     522*lam(5)*lam(6) + 9*lam(6)**2 +                          &
                     90*lam(4)*lam(7) + 198*lam(5)*lam(7) -                     &
                     90*lam(6)*lam(7) + 117*lam(7)**2 +                         &
                     36*lam(3)*                                                 &
                      (7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +           &
                     288*lam(8)**2 +                                            &
                     72*lam(1)*                                                 &
                      (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +              &
                        2*lam(8)) - 156*lam(4)*lam(9) -                         &
                     2100*lam(5)*lam(9) + 156*lam(6)*lam(9) -                   &
                     204*lam(7)*lam(9) - 56*lam(9)**2 +                         &
                     48*lam(2)*                                                 &
                      (3*lam(4) + 75*lam(5) -                                   &
                        3*(lam(6) + lam(7) + 4*lam(8)) + 16*lam(9)) -           &
                     36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -          &
                        10*lam(9))*lam(10) + 456*lam(10)**2) +                  &
                  2*(1 + Rqru)*tanhus*                                          &
                   (-3*qzt2*(1 + Rqru)*(1 + Rqsu)*                              &
                      (576*lam(2)**2 + 9*lam(4)**2 - 522*lam(4)*lam(5) -        &
                        5499*lam(5)**2 - 18*lam(4)*lam(6) +                     &
                        522*lam(5)*lam(6) + 9*lam(6)**2 +                       &
                        90*lam(4)*lam(7) + 198*lam(5)*lam(7) -                  &
                        90*lam(6)*lam(7) + 117*lam(7)**2 -                      &
                        36*lam(3)*                                              &
                         (7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +        &
                        288*lam(8)**2 -                                         &
                        72*lam(1)*                                              &
                         (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +           &
                         2*lam(8)) - 156*lam(4)*lam(9) -                        &
                        2100*lam(5)*lam(9) + 156*lam(6)*lam(9) -                &
                        204*lam(7)*lam(9) - 56*lam(9)**2 +                      &
                        48*lam(2)*                                              &
                         (3*lam(4) + 75*lam(5) -                                &
                         3*(lam(6) + lam(7) + 4*lam(8)) + 16*lam(9)) -          &
                        36*(5*lam(4) - 11*lam(5) - 5*lam(6) +                   &
                        11*lam(7) - 10*lam(9))*lam(10) + 456*lam(10)**2)        &
+ 2*m2*(864*lam(2)**2 + 45*lam(4)**2 + 8793*lam(5)**2 -                         &
                        414*lam(5)*lam(6) + 45*lam(6)**2 -                      &
                        1386*lam(5)*lam(7) - 234*lam(6)*lam(7) +                &
                        369*lam(7)**2 + 540*lam(5)*lam(8) +                     &
                        108*lam(6)*lam(8) - 540*lam(7)*lam(8) +                 &
                        3648*lam(5)*lam(9) - 192*lam(6)*lam(9) -                &
                        192*lam(7)*lam(9) + 720*lam(8)*lam(9) +                 &
                        608*lam(9)**2 +                                         &
                        6*lam(4)*                                               &
                         (69*lam(5) - 15*lam(6) + 39*lam(7) -                   &
                         18*lam(8) + 32*lam(9)) +                               &
                        72*lam(2)*                                              &
                         (3*lam(4) - 87*lam(5) - 3*lam(6) + 15*lam(7) -         &
                         26*lam(9) - 30*lam(10)) +                              &
                        72*(9*lam(8) - 8*lam(9))*lam(10) + 576*lam(10)**2)      &
))))) - expuISr*tanhur*(-(expdISs*qdB*quB*qzt2*tanhds*                          &
            (1728*dRqru*lam(2)**2 -                                             &
              3*dRqru*Rqsd*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                      &
               (576*lam(2)**2 + 9*lam(4)**2 + 882*lam(4)*lam(5) +               &
                 11241*lam(5)**2 - 18*lam(4)*lam(6) -                           &
                 882*lam(5)*lam(6) + 9*lam(6)**2 - 18*lam(4)*lam(7) -           &
                 882*lam(5)*lam(7) + 18*lam(6)*lam(7) + 9*lam(7)**2 -           &
                 144*lam(8)**2 -                                                &
                 72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +           &
                    2*lam(8)) -                                                 &
                 36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                 &
                    10*lam(8)) + 132*lam(4)*lam(9) +                            &
                 4740*lam(5)*lam(9) - 132*lam(6)*lam(9) -                       &
                 132*lam(7)*lam(9) + 412*lam(9)**2 -                            &
                 96*lam(2)*(3*lam(4) + 75*lam(5) -                              &
                    3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +               &
                 24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                &
                    10*lam(9))*lam(10) - 144*lam(10)**2) +                      &
              4*dRqsd*m2*(1 + Rqsd)*ss*                                         &
               (96*(9*lam(1)**2 - 48*lam(1)*lam(3) + 10*lam(3)**2) -            &
                 9*(lam(4)**2 - 11*lam(5)**2 + 14*lam(5)*lam(6) +               &
                    lam(6)**2 - 2*lam(4)*(7*lam(5) + lam(6))) -                 &
                 18*(7*lam(4) + 11*lam(5) - 7*lam(6))*lam(7) +                  &
                 99*lam(7)**2 -                                                 &
                 576*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +           &
                 960*lam(9)**2 +                                                &
                 192*(lam(4) + 5*lam(5) - lam(6) - 5*(lam(7) + lam(9)))*        &
                  lam(10) + 224*lam(10)**2) -                                   &
              4*dRqru*m2*(1 + Rqru)*sr*                                         &
               (1728*lam(2)**2 + 9*lam(4)**2 + 16857*lam(5)**2 -                &
                 1314*lam(5)*lam(6) + 9*lam(6)**2 - 1314*lam(5)*lam(7) +        &
                 18*lam(6)*lam(7) + 9*lam(7)**2 + 108*lam(5)*lam(8) -           &
                 108*lam(6)*lam(8) - 108*lam(7)*lam(8) +                        &
                 7296*lam(5)*lam(9) - 384*lam(6)*lam(9) -                       &
                 384*lam(7)*lam(9) + 360*lam(8)*lam(9) + 928*lam(9)**2 +        &
                 6*lam(4)*(219*lam(5) -                                         &
                    3*(lam(6) + lam(7) - 6*lam(8)) + 64*lam(9)) -               &
                 432*lam(8)*lam(10) -                                           &
                 72*lam(2)*(3*lam(4) + 147*lam(5) - 3*lam(6) -                  &
                    3*lam(7) + 22*lam(9) + 12*lam(10))) +                       &
              3*dRqsd*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                &
               (144*lam(2)**2 +                                                 &
                 216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +           &
                 9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                      &
                    2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                  &
                    10*lam(6)*lam(7) + 13*lam(7)**2 +                           &
                    lam(5)*(-10*lam(6) + 26*lam(7))) +                          &
                 288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -           &
                 192*lam(8)**2 -                                                &
                 128*lam(3)*(3*lam(4) + 6*lam(5) -                              &
                    3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) +             &
                 48*lam(1)*(3*lam(4) + 33*lam(5) -                              &
                    3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))) +         &
              3*dRqru*(-216*lam(1)*lam(4) - 36*lam(3)*lam(4) +                  &
                 9*lam(4)**2 + 36*lam(3)*lam(5) + 882*lam(4)*lam(5) +           &
                 11241*lam(5)**2 - 36*lam(3)*lam(6) - 18*lam(4)*lam(6) -        &
                 882*lam(5)*lam(6) + 9*lam(6)**2 + 36*lam(3)*lam(7) -           &
                 18*lam(4)*lam(7) - 882*lam(5)*lam(7) +                         &
                 18*lam(6)*lam(7) + 9*lam(7)**2 +                               &
                 72*lam(1)*(3*(lam(5) - lam(6) + lam(7)) - 2*lam(8)) -          &
                 360*lam(3)*lam(8) - 144*lam(8)**2 + 132*lam(4)*lam(9) +        &
                 4740*lam(5)*lam(9) - 132*lam(6)*lam(9) -                       &
                 132*lam(7)*lam(9) + 412*lam(9)**2 -                            &
                 96*lam(2)*(3*lam(4) + 75*lam(5) -                              &
                    3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +               &
                 24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                &
                    10*lam(9))*lam(10) - 144*lam(10)**2 -                       &
                 2*qzt2*(1 + Rqru)**2*sr*                                       &
                  (576*lam(2)**2 + 9*lam(4)**2 + 882*lam(4)*lam(5) +            &
                    11241*lam(5)**2 - 18*lam(4)*lam(6) -                        &
                    882*lam(5)*lam(6) + 9*lam(6)**2 - 18*lam(4)*lam(7) -        &
                    882*lam(5)*lam(7) + 18*lam(6)*lam(7) + 9*lam(7)**2 -        &
                    144*lam(8)**2 -                                             &
                    72*lam(1)*                                                  &
                     (3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) + 2*lam(8))       &
- 36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) + 10*lam(8)) +                   &
                    132*lam(4)*lam(9) + 4740*lam(5)*lam(9) -                    &
                    132*lam(6)*lam(9) - 132*lam(7)*lam(9) +                     &
                    412*lam(9)**2 -                                             &
                    96*lam(2)*                                                  &
                     (3*lam(4) + 75*lam(5) -                                    &
                       3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +            &
                    24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +             &
                       10*lam(9))*lam(10) - 144*lam(10)**2)))) -                &
         2*expdISs*qxy2*(1 + Rqsd)*tanhds**2*                                   &
          (3*dRqru*qdB*qxy2*(1 + Rqru)**2*                                      &
             (576*lam(2)**2 + 9*lam(4)**2 + 882*lam(4)*lam(5) +                 &
               11241*lam(5)**2 - 18*lam(4)*lam(6) - 882*lam(5)*lam(6) +         &
               9*lam(6)**2 - 18*lam(4)*lam(7) - 882*lam(5)*lam(7) +             &
               18*lam(6)*lam(7) + 9*lam(7)**2 - 144*lam(8)**2 +                 &
               72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +             &
                  2*lam(8)) +                                                   &
               36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                   &
                  10*lam(8)) + 132*lam(4)*lam(9) + 4740*lam(5)*lam(9) -         &
               132*lam(6)*lam(9) - 132*lam(7)*lam(9) + 412*lam(9)**2 -          &
               96*lam(2)*(3*lam(4) + 75*lam(5) -                                &
                  3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +                 &
               24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                  &
                  10*lam(9))*lam(10) - 144*lam(10)**2) +                        &
            dRqsd*quB*(6*m2*(32*                                                &
                   (9*lam(1)**2 - 48*lam(1)*lam(3) + 10*lam(3)**2) -            &
                  3*lam(4)**2 + 33*lam(5)**2 - 3*lam(6)**2 +                    &
                  6*lam(4)*(7*lam(5) + lam(6) - 7*lam(7)) +                     &
                  42*lam(6)*lam(7) + 33*lam(7)**2 -                             &
                  6*lam(5)*(7*lam(6) + 11*lam(7)) -                             &
                  192*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +          &
                  320*lam(9)**2) +                                              &
               384*m2*(lam(4) + 5*lam(5) - lam(6) - 5*(lam(7) + lam(9)))*       &
                lam(10) + 448*m2*lam(10)**2 +                                   &
               3*qzt2*(1 + Rqru)*(1 + Rqsd)*                                    &
                (144*lam(2)**2 +                                                &
                  216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +          &
                  9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                     &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                 &
                     10*lam(6)*lam(7) + 13*lam(7)**2 +                          &
                     lam(5)*(-10*lam(6) + 26*lam(7))) +                         &
                  288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -          &
                  192*lam(8)**2 -                                               &
                  128*lam(3)*(3*lam(4) + 6*lam(5) -                             &
                     3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) +            &
                  48*lam(1)*(3*lam(4) + 33*lam(5) -                             &
                     3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))))) +      &
         qdB*(dRqsu*expuISs*tanhus*                                             &
             (3*qzt2*(1 + Rqru)*                                                &
                (quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +                           &
                  2*qxy2*(1 + Rqsu)**2*tanhus)*                                 &
                (3*(48*lam(2)**2 +                                              &
                     192*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) +           &
                     64*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +            &
                     72*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +            &
                     3*(lam(4) + lam(5) + lam(6) + lam(7))**2) +                &
                  256*(3*lam(1) + lam(3))*lam(10)) +                            &
               4*m2*(1 + Rqsu)*(quB*qzt2*ss + qxy2*tanhus)*                     &
                (384*(3*lam(1) + lam(3))**2 +                                   &
                  45*(lam(4) - lam(5) - lam(6) + lam(7))**2 +                   &
                  192*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +             &
                  32*lam(10)**2)) +                                             &
            dRqru*(6*expdISs*qxy2*(1 + Rqru)*                                   &
                ((qxy2 + qzt2)*(1 + Rqru)*(1 + Rqsd)*                           &
                   (576*lam(2)**2 + 9*lam(4)**2 + 882*lam(4)*lam(5) +           &
                     11241*lam(5)**2 - 18*lam(4)*lam(6) -                       &
                     882*lam(5)*lam(6) + 9*lam(6)**2 -                          &
                     18*lam(4)*lam(7) - 882*lam(5)*lam(7) +                     &
                     18*lam(6)*lam(7) + 9*lam(7)**2 - 144*lam(8)**2 +           &
                     72*lam(1)*                                                 &
                      (3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +                &
                        2*lam(8)) +                                             &
                     36*lam(3)*                                                 &
                      (lam(4) - lam(5) + lam(6) - lam(7) + 10*lam(8)) +         &
                     132*lam(4)*lam(9) + 4740*lam(5)*lam(9) -                   &
                     132*lam(6)*lam(9) - 132*lam(7)*lam(9) +                    &
                     412*lam(9)**2 -                                            &
                     96*lam(2)*                                                 &
                      (3*lam(4) + 75*lam(5) -                                   &
                        3*(lam(6) + lam(7) - 2*lam(8)) + 16*lam(9)) +           &
                     24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +            &
                        10*lam(9))*lam(10) - 144*lam(10)**2) -                  &
                  24*m2*(-48*lam(2)**2 +                                        &
                     3*(12*lam(1)**2 + 4*lam(1)*lam(3) + 5*lam(3)**2) -         &
                     468*lam(5)**2 + 36*lam(5)*lam(6) +                         &
                     9*lam(4)*(-5*lam(5) + lam(6)) + 45*lam(5)*lam(7) -         &
                     9*lam(6)*lam(7) - 6*lam(5)*lam(8) +                        &
                     6*lam(6)*lam(8) + 15*lam(8)**2 - 192*lam(5)*lam(9) -       &
                     10*lam(8)*lam(9) - 8*lam(9)**2 + 12*lam(8)*lam(10) +       &
                     lam(2)*(6*lam(4) + 294*lam(5) - 6*lam(6) -                 &
                        6*lam(7) + 44*lam(9) + 24*lam(10)))) +                  &
               expuISs*(3*(1 + Rqsu)*                                           &
                   (-(quB*qzt2*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*tanhus*           &
                        (576*lam(2)**2 + 9*lam(4)**2 -                          &
                         522*lam(4)*lam(5) - 5499*lam(5)**2 -                   &
                         18*lam(4)*lam(6) + 522*lam(5)*lam(6) +                 &
                         9*lam(6)**2 + 90*lam(4)*lam(7) +                       &
                         198*lam(5)*lam(7) - 90*lam(6)*lam(7) +                 &
                         117*lam(7)**2 -                                        &
                         36*lam(3)*                                             &
                         (7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +        &
                         288*lam(8)**2 -                                        &
                         72*lam(1)*                                             &
                         (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +           &
                         2*lam(8)) - 156*lam(4)*lam(9) -                        &
                         2100*lam(5)*lam(9) + 156*lam(6)*lam(9) -               &
                         204*lam(7)*lam(9) - 56*lam(9)**2 +                     &
                         48*lam(2)*                                             &
                         (3*lam(4) + 75*lam(5) -                                &
                         3*(lam(6) + lam(7) + 4*lam(8)) + 16*lam(9)) -          &
                         36*(5*lam(4) - 11*lam(5) - 5*lam(6) +                  &
                        11*lam(7) - 10*lam(9))*lam(10) + 456*lam(10)**2))       &
- 2*qxy2*qzt2*(1 + Rqru)**2*(576*lam(2)**2 + 9*lam(4)**2 -                      &
                        522*lam(4)*lam(5) - 5499*lam(5)**2 -                    &
                        18*lam(4)*lam(6) + 522*lam(5)*lam(6) +                  &
                        9*lam(6)**2 + 90*lam(4)*lam(7) +                        &
                        198*lam(5)*lam(7) - 90*lam(6)*lam(7) +                  &
                        117*lam(7)**2 +                                         &
                        36*lam(3)*                                              &
                         (7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +        &
                        288*lam(8)**2 +                                         &
                        72*lam(1)*                                              &
                         (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +           &
                         2*lam(8)) - 156*lam(4)*lam(9) -                        &
                        2100*lam(5)*lam(9) + 156*lam(6)*lam(9) -                &
                        204*lam(7)*lam(9) - 56*lam(9)**2 +                      &
                        48*lam(2)*                                              &
                         (3*lam(4) + 75*lam(5) -                                &
                         3*(lam(6) + lam(7) + 4*lam(8)) + 16*lam(9)) -          &
                        36*(5*lam(4) - 11*lam(5) - 5*lam(6) +                   &
                        11*lam(7) - 10*lam(9))*lam(10) + 456*lam(10)**2)        &
+ 2*qxy2**2*(1 + Rqru)**2*(-1 + tanhus**2)*                                     &
                      (576*lam(2)**2 + 9*lam(4)**2 - 522*lam(4)*lam(5) -        &
                        5499*lam(5)**2 - 18*lam(4)*lam(6) +                     &
                        522*lam(5)*lam(6) + 9*lam(6)**2 +                       &
                        90*lam(4)*lam(7) + 198*lam(5)*lam(7) -                  &
                        90*lam(6)*lam(7) + 117*lam(7)**2 +                      &
                        36*lam(3)*                                              &
                         (7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +        &
                        288*lam(8)**2 +                                         &
                        72*lam(1)*                                              &
                         (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +           &
                         2*lam(8)) - 156*lam(4)*lam(9) -                        &
                        2100*lam(5)*lam(9) + 156*lam(6)*lam(9) -                &
                        204*lam(7)*lam(9) - 56*lam(9)**2 +                      &
                        48*lam(2)*                                              &
                         (3*lam(4) + 75*lam(5) -                                &
                         3*(lam(6) + lam(7) + 4*lam(8)) + 16*lam(9)) -          &
                        36*(5*lam(4) - 11*lam(5) - 5*lam(6) +                   &
                         11*lam(7) - 10*lam(9))*lam(10) + 456*lam(10)**2))      &
+ 4*m2*(1 + Rqru)*(36*qxy2*(24*lam(2)**2 +                                      &
                        6*(6*lam(1)**2 + 14*lam(1)*lam(3) +                     &
                        5*lam(3)**2) + 234*lam(5)**2 - 18*lam(5)*lam(6) +       &
                        9*lam(4)*(7*lam(5) + lam(6)) + 99*lam(5)*lam(7) +       &
                        45*lam(6)*lam(7) + 30*lam(5)*lam(8) +                   &
                        6*lam(6)*lam(8) - 6*lam(8)**2 +                         &
                        96*lam(5)*lam(9) + 20*lam(8)*lam(9) -                   &
                        8*lam(9)**2 +                                           &
                        2*lam(2)*                                               &
                         (3*lam(4) - 87*lam(5) - 3*lam(6) + 15*lam(7) -         &
                         26*lam(9) - 30*lam(10)) +                              &
                        6*(3*lam(8) + 8*lam(9))*lam(10) - 48*lam(10)**2) +      &
                     quB*qzt2*sr*tanhus*                                        &
                      (864*lam(2)**2 + 45*lam(4)**2 + 8793*lam(5)**2 -          &
                        414*lam(5)*lam(6) + 45*lam(6)**2 -                      &
                        1386*lam(5)*lam(7) - 234*lam(6)*lam(7) +                &
                        369*lam(7)**2 + 540*lam(5)*lam(8) +                     &
                        108*lam(6)*lam(8) - 540*lam(7)*lam(8) +                 &
                        3648*lam(5)*lam(9) - 192*lam(6)*lam(9) -                &
                        192*lam(7)*lam(9) + 720*lam(8)*lam(9) +                 &
                        608*lam(9)**2 +                                         &
                        6*lam(4)*                                               &
                         (69*lam(5) - 15*lam(6) + 39*lam(7) - 18*lam(8) +       &
                         32*lam(9)) +                                           &
                        72*lam(2)*                                              &
                         (3*lam(4) - 87*lam(5) - 3*lam(6) + 15*lam(7) -         &
                         26*lam(9) - 30*lam(10)) +                              &
                        72*(9*lam(8) - 8*lam(9))*lam(10) + 576*lam(10)**2))     &
)))))/3888.
  dtLAM(6)=                                                        &
   (-2*dRqrd*expdISr*quB*qxy2**2*(1 + Rqrd)**2*tanhdr**3*                       &
       (-(expuISs*(1 + Rqsu)*(-1 + tanhus**2)*                                  &
            (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +                   &
              72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -              &
                 2*lam(8)) + 24*                                                &
               (24*lam(2)**2 + 24*lam(2)*lam(8) +                               &
                 lam(8)*(17*lam(3) + 6*lam(8))) +                               &
              (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -          &
                 12*lam(10))**2)) +                                             &
         expdISs*(1 + Rqsd)*(-1 + tanhds**2)*                                   &
          (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +                     &
            117*lam(5)**2 + 1278*lam(4)*lam(6) - 3798*lam(5)*lam(6) +           &
            14265*lam(6)**2 + 90*lam(4)*lam(7) - 234*lam(5)*lam(7) +            &
            3798*lam(6)*lam(7) + 117*lam(7)**2 +                                &
            72*lam(1)*(-9*(lam(4) + lam(5) + lam(6) + lam(7)) -                 &
               2*lam(8)) - 288*lam(8)**2 -                                      &
            12*lam(3)*(21*lam(4) + 69*lam(5) + 21*lam(6) + 69*lam(7) +          &
               32*lam(8)) - 228*lam(4)*lam(9) + 564*lam(5)*lam(9) -             &
            6684*lam(6)*lam(9) - 564*lam(7)*lam(9) + 568*lam(9)**2 -            &
            144*lam(2)*(3*lam(4) - 9*lam(5) + 69*lam(6) + 9*lam(7) -            &
               4*(lam(8) + 4*lam(9) - 8*lam(10))) +                             &
            12*(63*lam(4) - 177*lam(5) + 1089*lam(6) + 177*lam(7) -             &
               286*lam(9))*lam(10) + 2616*lam(10)**2)) +                        &
      2*dRqru*expuISr*qdB*qxy2**2*(1 + Rqru)**2*tanhur**3*                      &
       (expdISs*(1 + Rqsd)*(-1 + tanhds**2)*                                    &
          (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +                     &
            72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -                &
               2*lam(8)) + 24*                                                  &
             (24*lam(2)**2 + 24*lam(2)*lam(8) +                                 &
               lam(8)*(17*lam(3) + 6*lam(8))) +                                 &
            (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -            &
               12*lam(10))**2) -                                                &
         expuISs*(1 + Rqsu)*(-1 + tanhus**2)*                                   &
          (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +                     &
            117*lam(5)**2 + 1278*lam(4)*lam(6) - 3798*lam(5)*lam(6) +           &
            14265*lam(6)**2 + 90*lam(4)*lam(7) - 234*lam(5)*lam(7) +            &
            3798*lam(6)*lam(7) + 117*lam(7)**2 +                                &
            72*lam(1)*(-9*(lam(4) + lam(5) + lam(6) + lam(7)) -                 &
               2*lam(8)) - 288*lam(8)**2 -                                      &
            12*lam(3)*(21*lam(4) + 69*lam(5) + 21*lam(6) + 69*lam(7) +          &
               32*lam(8)) - 228*lam(4)*lam(9) + 564*lam(5)*lam(9) -             &
            6684*lam(6)*lam(9) - 564*lam(7)*lam(9) + 568*lam(9)**2 -            &
            144*lam(2)*(3*lam(4) - 9*lam(5) + 69*lam(6) + 9*lam(7) -            &
               4*(lam(8) + 4*lam(9) - 8*lam(10))) +                             &
            12*(63*lam(4) - 177*lam(5) + 1089*lam(6) + 177*lam(7) -             &
               286*lam(9))*lam(10) + 2616*lam(10)**2)) +                        &
      expdISr*expdISs*quB*(3*dRqsd*(1 + Rqrd)*                                  &
          (qdB*(-1 + 2*qzt2*(1 + Rqsd)**2*ss) +                                 &
            2*qxy2*(1 + Rqsd)**2*tanhds)*(qxy2 + qzt2 - qxy2*tanhds**2)*        &
          (3*(16*lam(2)**2 +                                                    &
               192*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) +                 &
               64*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) -                  &
               72*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                  &
               (lam(4) + lam(5) + lam(6) + lam(7))**2) +                        &
            256*(3*lam(1) + lam(3))*lam(10)) +                                  &
         48*dRqsd*m2*(1 + Rqsd)*(qdB*qzt2*ss + qxy2*tanhds)*                    &
          (-108*lam(2)**2 + 32*(3*lam(1) + lam(3))**2 -                         &
            27*(lam(5) + lam(6))*(lam(4) + lam(7)) +                            &
            6*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) + 24*lam(10)**2)       &
- dRqrd*qdB*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                          &
          (qxy2 + qzt2 - qxy2*tanhds**2)*                                       &
          (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +                     &
            117*lam(5)**2 + 1278*lam(4)*lam(6) - 3798*lam(5)*lam(6) +           &
            14265*lam(6)**2 + 90*lam(4)*lam(7) - 234*lam(5)*lam(7) +            &
            3798*lam(6)*lam(7) + 117*lam(7)**2 +                                &
            72*lam(1)*(-9*(lam(4) + lam(5) + lam(6) + lam(7)) -                 &
               2*lam(8)) - 288*lam(8)**2 -                                      &
            12*lam(3)*(21*lam(4) + 69*lam(5) + 21*lam(6) + 69*lam(7) +          &
               32*lam(8)) - 228*lam(4)*lam(9) + 564*lam(5)*lam(9) -             &
            6684*lam(6)*lam(9) - 564*lam(7)*lam(9) + 568*lam(9)**2 -            &
            144*lam(2)*(3*lam(4) - 9*lam(5) + 69*lam(6) + 9*lam(7) -            &
               4*(lam(8) + 4*lam(9) - 8*lam(10))) +                             &
            12*(63*lam(4) - 177*lam(5) + 1089*lam(6) + 177*lam(7) -             &
               286*lam(9))*lam(10) + 2616*lam(10)**2) -                         &
         48*dRqrd*m2*qdB*qzt2*(1 + Rqrd)*sr*                                    &
          (108*lam(1)**2 + 72*lam(2)**2 + 84*lam(1)*lam(3) +                    &
            46*lam(3)**2 - 162*lam(5)*lam(6) + 594*lam(6)**2 +                  &
            27*lam(4)*(lam(5) + 3*lam(6)) + 135*lam(5)*lam(7) +                 &
            189*lam(6)*lam(7) + 30*lam(5)*lam(8) + 6*lam(6)*lam(8) +            &
            10*lam(8)**2 - 288*lam(6)*lam(9) + 20*lam(8)*lam(9) -               &
            24*lam(9)**2 + 6*(96*lam(6) + 3*lam(8) + 8*lam(9))*lam(10) -        &
            48*lam(10)**2 - 2*lam(2)*                                           &
             (3*lam(4) - 15*lam(5) + 213*lam(6) + 15*lam(7) - 38*lam(9) +       &
               126*lam(10)))) +                                                 &
      expuISr*expuISs*qdB*(3*dRqsu*(1 + Rqru)*                                  &
          (quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +                                 &
            2*qxy2*(1 + Rqsu)**2*tanhus)*(qxy2 + qzt2 - qxy2*tanhus**2)*        &
          (3*(16*lam(2)**2 +                                                    &
               192*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) +                 &
               64*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) -                  &
               72*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                  &
               (lam(4) + lam(5) + lam(6) + lam(7))**2) +                        &
            256*(3*lam(1) + lam(3))*lam(10)) +                                  &
         48*dRqsu*m2*(1 + Rqsu)*(quB*qzt2*ss + qxy2*tanhus)*                    &
          (-108*lam(2)**2 + 32*(3*lam(1) + lam(3))**2 -                         &
            27*(lam(5) + lam(6))*(lam(4) + lam(7)) +                            &
            6*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) + 24*lam(10)**2)       &
- dRqru*quB*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                          &
          (qxy2 + qzt2 - qxy2*tanhus**2)*                                       &
          (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +                     &
            117*lam(5)**2 + 1278*lam(4)*lam(6) - 3798*lam(5)*lam(6) +           &
            14265*lam(6)**2 + 90*lam(4)*lam(7) - 234*lam(5)*lam(7) +            &
            3798*lam(6)*lam(7) + 117*lam(7)**2 +                                &
            72*lam(1)*(-9*(lam(4) + lam(5) + lam(6) + lam(7)) -                 &
               2*lam(8)) - 288*lam(8)**2 -                                      &
            12*lam(3)*(21*lam(4) + 69*lam(5) + 21*lam(6) + 69*lam(7) +          &
               32*lam(8)) - 228*lam(4)*lam(9) + 564*lam(5)*lam(9) -             &
            6684*lam(6)*lam(9) - 564*lam(7)*lam(9) + 568*lam(9)**2 -            &
            144*lam(2)*(3*lam(4) - 9*lam(5) + 69*lam(6) + 9*lam(7) -            &
               4*(lam(8) + 4*lam(9) - 8*lam(10))) +                             &
            12*(63*lam(4) - 177*lam(5) + 1089*lam(6) + 177*lam(7) -             &
               286*lam(9))*lam(10) + 2616*lam(10)**2) -                         &
         48*dRqru*m2*quB*qzt2*(1 + Rqru)*sr*                                    &
          (108*lam(1)**2 + 72*lam(2)**2 + 84*lam(1)*lam(3) +                    &
            46*lam(3)**2 - 162*lam(5)*lam(6) + 594*lam(6)**2 +                  &
            27*lam(4)*(lam(5) + 3*lam(6)) + 135*lam(5)*lam(7) +                 &
            189*lam(6)*lam(7) + 30*lam(5)*lam(8) + 6*lam(6)*lam(8) +            &
            10*lam(8)**2 - 288*lam(6)*lam(9) + 20*lam(8)*lam(9) -               &
            24*lam(9)**2 + 6*(96*lam(6) + 3*lam(8) + 8*lam(9))*lam(10) -        &
            48*lam(10)**2 - 2*lam(2)*                                           &
             (3*lam(4) - 15*lam(5) + 213*lam(6) + 15*lam(7) - 38*lam(9) +       &
               126*lam(10)))) +                                                 &
      expdISs*expuISr*quB*(dRqru*qdB*(1 + Rqsd)*                                &
          (-1 + 2*qzt2*(1 + Rqru)**2*sr)*(qxy2 + qzt2 - qxy2*tanhds**2)*        &
          (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +                     &
            72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -                &
               2*lam(8)) + 24*                                                  &
             (24*lam(2)**2 + 24*lam(2)*lam(8) +                                 &
               lam(8)*(17*lam(3) + 6*lam(8))) +                                 &
            (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -            &
               12*lam(10))**2) -                                                &
         48*dRqru*m2*qdB*qzt2*(1 + Rqru)*sr*                                    &
          (108*lam(1)**2 - 12*lam(1)*lam(3) + 17*lam(3)**2 +                    &
            27*lam(4)*(-lam(5) + lam(6)) + 27*lam(5)*lam(7) -                   &
            27*lam(6)*lam(7) + 6*lam(5)*lam(8) - 6*lam(6)*lam(8) +              &
            17*lam(8)**2 + 10*lam(8)*lam(9) - 24*lam(9)**2 +                    &
            lam(2)*(6*lam(4) + 6*lam(5) - 6*lam(6) - 6*lam(7) +                 &
               20*lam(9) - 24*lam(10)) - 12*lam(8)*lam(10)) -                   &
         3*dRqsd*(1 + Rqru)*(qdB*(-1 + 2*qzt2*(1 + Rqsd)**2*ss) +               &
            2*qxy2*(1 + Rqsd)**2*tanhds)*(qxy2 + qzt2 - qxy2*tanhds**2)*        &
          (48*lam(2)**2 + 216*lam(2)*                                           &
             (lam(4) - lam(5) + lam(6) - lam(7)) +                              &
            3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                    &
               lam(6)**2 + 2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +           &
               26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +            &
            96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -                 &
            64*lam(8)**2 + 128*lam(3)*                                          &
             (3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) -                       &
            48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -              &
               16*lam(9) + 16*lam(10))) +                                       &
         144*dRqsd*m2*(1 + Rqsd)*(qdB*qzt2*ss + qxy2*tanhds)*                   &
          (72*lam(1)**2 - 36*lam(2)**2 + 9*lam(4)*(lam(5) - lam(6)) +           &
            9*(-5*lam(5) + lam(6))*lam(7) - 16*lam(8)**2 -                      &
            2*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7) +                   &
               16*lam(8)) + 8*                                                  &
             (2*(lam(3)**2 + lam(9)**2) - 2*lam(9)*lam(10) + lam(10)**2)))      &
+ expdISr*expuISs*qdB*(dRqrd*quB*(1 + Rqsu)*                                    &
          (-1 + 2*qzt2*(1 + Rqrd)**2*sr)*(qxy2 + qzt2 - qxy2*tanhus**2)*        &
          (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +                     &
            72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -                &
               2*lam(8)) + 24*                                                  &
             (24*lam(2)**2 + 24*lam(2)*lam(8) +                                 &
               lam(8)*(17*lam(3) + 6*lam(8))) +                                 &
            (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -            &
               12*lam(10))**2) -                                                &
         48*dRqrd*m2*quB*qzt2*(1 + Rqrd)*sr*                                    &
          (108*lam(1)**2 - 12*lam(1)*lam(3) + 17*lam(3)**2 +                    &
            27*lam(4)*(-lam(5) + lam(6)) + 27*lam(5)*lam(7) -                   &
            27*lam(6)*lam(7) + 6*lam(5)*lam(8) - 6*lam(6)*lam(8) +              &
            17*lam(8)**2 + 10*lam(8)*lam(9) - 24*lam(9)**2 +                    &
            lam(2)*(6*lam(4) + 6*lam(5) - 6*lam(6) - 6*lam(7) +                 &
               20*lam(9) - 24*lam(10)) - 12*lam(8)*lam(10)) -                   &
         3*dRqsu*(1 + Rqrd)*(quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +               &
            2*qxy2*(1 + Rqsu)**2*tanhus)*(qxy2 + qzt2 - qxy2*tanhus**2)*        &
          (48*lam(2)**2 + 216*lam(2)*                                           &
             (lam(4) - lam(5) + lam(6) - lam(7)) +                              &
            3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                    &
               lam(6)**2 + 2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +           &
               26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +            &
            96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -                 &
            64*lam(8)**2 + 128*lam(3)*                                          &
             (3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) -                       &
            48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -              &
               16*lam(9) + 16*lam(10))) +                                       &
         144*dRqsu*m2*(1 + Rqsu)*(quB*qzt2*ss + qxy2*tanhus)*                   &
          (72*lam(1)**2 - 36*lam(2)**2 + 9*lam(4)*(lam(5) - lam(6)) +           &
            9*(-5*lam(5) + lam(6))*lam(7) - 16*lam(8)**2 -                      &
            2*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7) +                   &
               16*lam(8)) + 8*                                                  &
             (2*(lam(3)**2 + lam(9)**2) - 2*lam(9)*lam(10) + lam(10)**2)))      &
- expdISr*qxy2*tanhdr**2*(-6*dRqsd*expdISs*quB*qxy2*(1 + Rqrd)*                 &
          (1 + Rqsd)**2*tanhds**3*                                              &
          (48*lam(2)**2 - 216*lam(2)*                                           &
             (lam(4) + lam(5) + lam(6) + lam(7)) +                              &
            3*(192*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) +                 &
               64*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +                  &
               (lam(4) + lam(5) + lam(6) + lam(7))**2) +                        &
            256*(3*lam(1) + lam(3))*lam(10)) +                                  &
         6*dRqsu*expuISs*qdB*qxy2*(1 + Rqrd)*(1 + Rqsu)**2*tanhus**3*           &
          (48*lam(2)**2 + 216*lam(2)*                                           &
             (lam(4) - lam(5) + lam(6) - lam(7)) +                              &
            3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                    &
               lam(6)**2 + 2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +           &
               26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +            &
            96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -                 &
            64*lam(8)**2 + 128*lam(3)*                                          &
             (3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) -                       &
            48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -              &
               16*lam(9) + 16*lam(10))) +                                       &
         expdISs*qdB*quB*tanhds**2*                                             &
          (-3*dRqsd*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                  &
             (48*lam(2)**2 -                                                    &
               216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                 &
               3*(192*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) +              &
                  64*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +               &
                  (lam(4) + lam(5) + lam(6) + lam(7))**2) +                     &
               256*(3*lam(1) + lam(3))*lam(10)) +                               &
            dRqrd*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                    &
             (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +                  &
               117*lam(5)**2 + 1278*lam(4)*lam(6) - 3798*lam(5)*lam(6) +        &
               14265*lam(6)**2 + 90*lam(4)*lam(7) - 234*lam(5)*lam(7) +         &
               3798*lam(6)*lam(7) + 117*lam(7)**2 +                             &
               72*lam(1)*(-9*(lam(4) + lam(5) + lam(6) + lam(7)) -              &
                  2*lam(8)) - 288*lam(8)**2 -                                   &
               12*lam(3)*(21*lam(4) + 69*lam(5) + 21*lam(6) +                   &
                  69*lam(7) + 32*lam(8)) - 228*lam(4)*lam(9) +                  &
               564*lam(5)*lam(9) - 6684*lam(6)*lam(9) -                         &
               564*lam(7)*lam(9) + 568*lam(9)**2 -                              &
               144*lam(2)*(3*lam(4) - 9*lam(5) + 69*lam(6) + 9*lam(7) -         &
                  4*(lam(8) + 4*lam(9) - 8*lam(10))) +                          &
               12*(63*lam(4) - 177*lam(5) + 1089*lam(6) + 177*lam(7) -          &
                  286*lam(9))*lam(10) + 2616*lam(10)**2)) -                     &
         expuISs*qdB*quB*tanhus**2*                                             &
          (dRqrd*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                     &
             (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +                  &
               72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -             &
                  2*lam(8)) +                                                   &
               24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                            &
                  lam(8)*(17*lam(3) + 6*lam(8))) +                              &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                     &
                  10*lam(9) - 12*lam(10))**2) -                                 &
            3*dRqsu*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                  &
             (48*lam(2)**2 + 216*lam(2)*                                        &
                (lam(4) - lam(5) + lam(6) - lam(7)) +                           &
               3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                 &
                  lam(6)**2 +                                                   &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                    &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -              &
               64*lam(8)**2 +                                                   &
               128*lam(3)*(3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) -          &
               48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -           &
                  16*lam(9) + 16*lam(10)))) -                                   &
         2*expuISs*(1 + Rqrd)*tanhus*                                           &
          (3*dRqsu*qdB*qxy2*(1 + Rqsu)**2*                                      &
             (48*lam(2)**2 +                                                    &
               216*lam(2)*(lam(4) - lam(5) + lam(6) - lam(7)) +                 &
               3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                 &
                  lam(6)**2 +                                                   &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                    &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -              &
               64*lam(8)**2 +                                                   &
               128*lam(3)*(3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) -          &
               48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -           &
                  16*lam(9) + 16*lam(10))) +                                    &
            dRqrd*quB*(qzt2*(1 + Rqrd)*(1 + Rqsu)*                              &
                (36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -                &
                  72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -          &
                     2*lam(8)) +                                                &
                  24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                         &
                     lam(8)*(-17*lam(3) + 6*lam(8))) +                          &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) -                              &
               6*m2*(9*lam(4)**2 + 9*lam(5)**2 - 18*lam(5)*lam(6) +             &
                  9*lam(6)**2 - 18*lam(5)*lam(7) + 18*lam(6)*lam(7) +           &
                  9*lam(7)**2 + 12*lam(5)*lam(8) - 12*lam(6)*lam(8) -           &
                  12*lam(7)*lam(8) +                                            &
                  6*lam(4)*(3*lam(5) - 3*(lam(6) + lam(7)) + 2*lam(8)) +        &
                  40*lam(8)*lam(9) + 32*lam(9)**2 +                             &
                  8*lam(2)*(3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +         &
                     10*lam(9) - 12*lam(10)) - 48*lam(8)*lam(10)))) +           &
         qdB*quB*(-(dRqrd*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                       &
               (-(expuISs*(1 + Rqsu)*                                           &
                    (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +           &
                      72*lam(1)*                                                &
                       (9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -               &
                        2*lam(8)) +                                             &
                      24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                     &
                        lam(8)*(17*lam(3) + 6*lam(8))) +                        &
                      (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +              &
                        10*lam(9) - 12*lam(10))**2)) +                          &
                 expdISs*(1 + Rqsd)*                                            &
                  (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +             &
                    117*lam(5)**2 + 1278*lam(4)*lam(6) -                        &
                    3798*lam(5)*lam(6) + 14265*lam(6)**2 +                      &
                    90*lam(4)*lam(7) - 234*lam(5)*lam(7) +                      &
                    3798*lam(6)*lam(7) + 117*lam(7)**2 +                        &
                    72*lam(1)*                                                  &
                     (-9*(lam(4) + lam(5) + lam(6) + lam(7)) -                  &
                       2*lam(8)) - 288*lam(8)**2 -                              &
                    12*lam(3)*                                                  &
                     (21*lam(4) + 69*lam(5) + 21*lam(6) + 69*lam(7) +           &
                       32*lam(8)) - 228*lam(4)*lam(9) +                         &
                    564*lam(5)*lam(9) - 6684*lam(6)*lam(9) -                    &
                    564*lam(7)*lam(9) + 568*lam(9)**2 -                         &
                    144*lam(2)*                                                 &
                     (3*lam(4) - 9*lam(5) + 69*lam(6) + 9*lam(7) -              &
                       4*(lam(8) + 4*lam(9) - 8*lam(10))) +                     &
                    12*(63*lam(4) - 177*lam(5) + 1089*lam(6) +                  &
                       177*lam(7) - 286*lam(9))*lam(10) +                       &
                    2616*lam(10)**2))) +                                        &
            3*(1 + Rqrd)*(dRqsd*expdISs*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*         &
                (48*lam(2)**2 -                                                 &
                  216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +              &
                  3*(192*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) +           &
                     64*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +            &
                     (lam(4) + lam(5) + lam(6) + lam(7))**2) +                  &
                  256*(3*lam(1) + lam(3))*lam(10)) -                            &
               dRqsu*expuISs*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                    &
                (48*lam(2)**2 +                                                 &
                  216*lam(2)*(lam(4) - lam(5) + lam(6) - lam(7)) +              &
                  3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +              &
                     lam(6)**2 +                                                &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                 &
                     26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2)        &
+ 96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) - 64*lam(8)**2 +            &
                  128*lam(3)*                                                   &
                   (3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) -                 &
                  48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -        &
                     16*lam(9) + 16*lam(10))))) +                               &
         2*expdISs*quB*(1 + Rqrd)*tanhds*                                       &
          (3*dRqsd*qxy2*(1 + Rqsd)**2*                                          &
             (48*lam(2)**2 - 216*lam(2)*                                        &
                (lam(4) + lam(5) + lam(6) + lam(7)) +                           &
               3*(192*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) +              &
                  64*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +               &
                  (lam(4) + lam(5) + lam(6) + lam(7))**2) +                     &
               256*(3*lam(1) + lam(3))*lam(10)) +                               &
            dRqrd*(qzt2*(1 + Rqrd)*(1 + Rqsd)*                                  &
                (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +               &
                  117*lam(5)**2 + 1278*lam(4)*lam(6) -                          &
                  3798*lam(5)*lam(6) + 14265*lam(6)**2 +                        &
                  90*lam(4)*lam(7) - 234*lam(5)*lam(7) +                        &
                  3798*lam(6)*lam(7) + 117*lam(7)**2 - 288*lam(8)**2 +          &
                  72*lam(1)*(9*(lam(4) + lam(5) + lam(6) + lam(7)) +            &
                     2*lam(8)) +                                                &
                  12*lam(3)*(21*lam(4) + 69*lam(5) + 21*lam(6) +                &
                     69*lam(7) + 32*lam(8)) - 228*lam(4)*lam(9) +               &
                  564*lam(5)*lam(9) - 6684*lam(6)*lam(9) -                      &
                  564*lam(7)*lam(9) + 568*lam(9)**2 -                           &
                  144*lam(2)*                                                   &
                   (3*lam(4) - 9*lam(5) + 69*lam(6) + 9*lam(7) -                &
                     4*(lam(8) + 4*lam(9) - 8*lam(10))) +                       &
                  12*(63*lam(4) - 177*lam(5) + 1089*lam(6) +                    &
                     177*lam(7) - 286*lam(9))*lam(10) + 2616*lam(10)**2)        &
+ 2*m2*(864*lam(2)**2 + 45*lam(4)**2 + 297*lam(5)**2 -                          &
                  1782*lam(5)*lam(6) + 7173*lam(6)**2 -                         &
                  594*lam(5)*lam(7) + 1782*lam(6)*lam(7) +                      &
                  297*lam(7)**2 + 180*lam(5)*lam(8) + 36*lam(6)*lam(8) -        &
                  180*lam(7)*lam(8) + 576*lam(5)*lam(9) -                       &
                  3264*lam(6)*lam(9) - 576*lam(7)*lam(9) +                      &
                  240*lam(8)*lam(9) + 608*lam(9)**2 -                           &
                  6*lam(4)*(27*lam(5) - 93*lam(6) - 27*lam(7) +                 &
                     6*lam(8) + 32*lam(9) - 64*lam(10)) -                       &
                  8*(144*lam(5) -                                               &
                     3*(272*lam(6) + 48*lam(7) + 9*lam(8)) + 280*lam(9))*       &
                   lam(10) + 2240*lam(10)**2 -                                  &
                  24*lam(2)*(3*lam(4) - 15*lam(5) + 213*lam(6) +                &
                     15*lam(7) - 38*lam(9) + 126*lam(10)))))) -                 &
      expdISr*tanhdr*(-2*dRqrd*quB*qxy2*(1 + Rqrd)*                             &
          (expuISs*(qxy2 + qzt2)*(1 + Rqrd)*(1 + Rqsu)*                         &
             (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +                  &
               72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -             &
                  2*lam(8)) +                                                   &
               24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                            &
                  lam(8)*(17*lam(3) + 6*lam(8))) +                              &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                     &
                  10*lam(9) - 12*lam(10))**2) -                                 &
            expdISs*(qxy2 + qzt2)*(1 + Rqrd)*(1 + Rqsd)*                        &
             (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +                  &
               117*lam(5)**2 + 1278*lam(4)*lam(6) -                             &
               3798*lam(5)*lam(6) + 14265*lam(6)**2 +                           &
               90*lam(4)*lam(7) - 234*lam(5)*lam(7) +                           &
               3798*lam(6)*lam(7) + 117*lam(7)**2 +                             &
               72*lam(1)*(-9*(lam(4) + lam(5) + lam(6) + lam(7)) -              &
                  2*lam(8)) - 288*lam(8)**2 -                                   &
               12*lam(3)*(21*lam(4) + 69*lam(5) + 21*lam(6) +                   &
                  69*lam(7) + 32*lam(8)) - 228*lam(4)*lam(9) +                  &
               564*lam(5)*lam(9) - 6684*lam(6)*lam(9) -                         &
               564*lam(7)*lam(9) + 568*lam(9)**2 -                              &
               144*lam(2)*(3*lam(4) - 9*lam(5) + 69*lam(6) +                    &
                  9*lam(7) - 4*(lam(8) + 4*lam(9) - 8*lam(10))) +               &
               12*(63*lam(4) - 177*lam(5) + 1089*lam(6) +                       &
                  177*lam(7) - 286*lam(9))*lam(10) + 2616*lam(10)**2) +         &
            24*expuISs*m2*(-108*lam(1)**2 + 12*lam(1)*lam(3) -                  &
               17*lam(3)**2 + 27*lam(4)*(lam(5) - lam(6)) -                     &
               27*lam(5)*lam(7) + 27*lam(6)*lam(7) - 6*lam(5)*lam(8) +          &
               6*lam(6)*lam(8) - 17*lam(8)**2 - 10*lam(8)*lam(9) +              &
               24*lam(9)**2 + 12*lam(8)*lam(10) +                               &
               lam(2)*(-6*lam(4) - 6*lam(5) + 6*lam(6) + 6*lam(7) -             &
                  20*lam(9) + 24*lam(10))) -                                    &
            24*expdISs*m2*(108*lam(1)**2 + 72*lam(2)**2 +                       &
               84*lam(1)*lam(3) + 46*lam(3)**2 - 162*lam(5)*lam(6) +            &
               594*lam(6)**2 + 27*lam(4)*(lam(5) + 3*lam(6)) +                  &
               135*lam(5)*lam(7) + 189*lam(6)*lam(7) +                          &
               30*lam(5)*lam(8) + 6*lam(6)*lam(8) + 10*lam(8)**2 -              &
               288*lam(6)*lam(9) + 20*lam(8)*lam(9) - 24*lam(9)**2 +            &
               6*(96*lam(6) + 3*lam(8) + 8*lam(9))*lam(10) -                    &
               48*lam(10)**2 -                                                  &
               2*lam(2)*(3*lam(4) - 15*lam(5) + 213*lam(6) + 15*lam(7) -        &
                  38*lam(9) + 126*lam(10)))) +                                  &
         2*expdISs*quB*qxy2*(1 + Rqsd)*tanhds**2*                               &
          (-(dRqrd*qxy2*(1 + Rqrd)**2*                                          &
               (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +                &
                 117*lam(5)**2 + 1278*lam(4)*lam(6) -                           &
                 3798*lam(5)*lam(6) + 14265*lam(6)**2 +                         &
                 90*lam(4)*lam(7) - 234*lam(5)*lam(7) +                         &
                 3798*lam(6)*lam(7) + 117*lam(7)**2 +                           &
                 72*lam(1)*(-9*(lam(4) + lam(5) + lam(6) + lam(7)) -            &
                    2*lam(8)) - 288*lam(8)**2 -                                 &
                 12*lam(3)*(21*lam(4) + 69*lam(5) + 21*lam(6) +                 &
                    69*lam(7) + 32*lam(8)) - 228*lam(4)*lam(9) +                &
                 564*lam(5)*lam(9) - 6684*lam(6)*lam(9) -                       &
                 564*lam(7)*lam(9) + 568*lam(9)**2 -                            &
                 144*lam(2)*(3*lam(4) - 9*lam(5) + 69*lam(6) +                  &
                    9*lam(7) - 4*(lam(8) + 4*lam(9) - 8*lam(10))) +             &
                 12*(63*lam(4) - 177*lam(5) + 1089*lam(6) +                     &
                    177*lam(7) - 286*lam(9))*lam(10) + 2616*lam(10)**2))        &
+ dRqsd*(6*m2*(128*(3*lam(1) + lam(3))**2 +                                     &
                  15*(lam(4) - lam(5) - lam(6) + lam(7))**2) +                  &
               384*m2*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +             &
               64*m2*lam(10)**2 -                                               &
               3*qzt2*(1 + Rqrd)*(1 + Rqsd)*                                    &
                (48*lam(2)**2 -                                                 &
                  216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +              &
                  3*(lam(4) + lam(5) + lam(6) + lam(7))**2 -                    &
                  192*lam(1)*                                                   &
                   (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))        &
- 64*lam(3)*(3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))))) +          &
         expdISs*qdB*quB*qzt2*tanhds*                                           &
          (-576*dRqrd*lam(2)**2 +                                               &
            dRqsd*(-3*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                &
                (3*(16*lam(2)**2 -                                              &
                     192*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) -           &
                     64*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) -            &
                     72*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +            &
                     (lam(4) + lam(5) + lam(6) + lam(7))**2) -                  &
                  256*(3*lam(1) + lam(3))*lam(10)) +                            &
               4*m2*(1 + Rqsd)*ss*                                              &
                (384*(3*lam(1) + lam(3))**2 +                                   &
                  45*(lam(4) - lam(5) - lam(6) + lam(7))**2 +                   &
                  192*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +             &
                  32*lam(10)**2)) +                                             &
            dRqrd*(1152*qzt2*sr*lam(2)**2 +                                     &
               2304*qzt2*Rqrd*sr*lam(2)**2 +                                    &
               1152*qzt2*Rqrd**2*sr*lam(2)**2 - 648*lam(1)*lam(4) +             &
               1296*qzt2*sr*lam(1)*lam(4) +                                     &
               2592*qzt2*Rqrd*sr*lam(1)*lam(4) +                                &
               1296*qzt2*Rqrd**2*sr*lam(1)*lam(4) + 432*lam(2)*lam(4) -         &
               864*qzt2*sr*lam(2)*lam(4) -                                      &
               1728*qzt2*Rqrd*sr*lam(2)*lam(4) -                                &
               864*qzt2*Rqrd**2*sr*lam(2)*lam(4) - 252*lam(3)*lam(4) +          &
               504*qzt2*sr*lam(3)*lam(4) +                                      &
               1008*qzt2*Rqrd*sr*lam(3)*lam(4) +                                &
               504*qzt2*Rqrd**2*sr*lam(3)*lam(4) - 9*lam(4)**2 +                &
               18*qzt2*sr*lam(4)**2 + 36*qzt2*Rqrd*sr*lam(4)**2 +               &
               18*qzt2*Rqrd**2*sr*lam(4)**2 - 648*lam(1)*lam(5) +               &
               1296*qzt2*sr*lam(1)*lam(5) +                                     &
               2592*qzt2*Rqrd*sr*lam(1)*lam(5) +                                &
               1296*qzt2*Rqrd**2*sr*lam(1)*lam(5) - 1296*lam(2)*lam(5) +        &
               2592*qzt2*sr*lam(2)*lam(5) +                                     &
               5184*qzt2*Rqrd*sr*lam(2)*lam(5) +                                &
               2592*qzt2*Rqrd**2*sr*lam(2)*lam(5) - 828*lam(3)*lam(5) +         &
               1656*qzt2*sr*lam(3)*lam(5) +                                     &
               3312*qzt2*Rqrd*sr*lam(3)*lam(5) +                                &
               1656*qzt2*Rqrd**2*sr*lam(3)*lam(5) + 90*lam(4)*lam(5) -          &
               180*qzt2*sr*lam(4)*lam(5) -                                      &
               360*qzt2*Rqrd*sr*lam(4)*lam(5) -                                 &
               180*qzt2*Rqrd**2*sr*lam(4)*lam(5) - 117*lam(5)**2 +              &
               234*qzt2*sr*lam(5)**2 + 468*qzt2*Rqrd*sr*lam(5)**2 +             &
               234*qzt2*Rqrd**2*sr*lam(5)**2 - 648*lam(1)*lam(6) +              &
               1296*qzt2*sr*lam(1)*lam(6) +                                     &
               2592*qzt2*Rqrd*sr*lam(1)*lam(6) +                                &
               1296*qzt2*Rqrd**2*sr*lam(1)*lam(6) + 9936*lam(2)*lam(6) -        &
               19872*qzt2*sr*lam(2)*lam(6) -                                    &
               39744*qzt2*Rqrd*sr*lam(2)*lam(6) -                               &
               19872*qzt2*Rqrd**2*sr*lam(2)*lam(6) - 252*lam(3)*lam(6) +        &
               504*qzt2*sr*lam(3)*lam(6) +                                      &
               1008*qzt2*Rqrd*sr*lam(3)*lam(6) +                                &
               504*qzt2*Rqrd**2*sr*lam(3)*lam(6) - 1278*lam(4)*lam(6) +         &
               2556*qzt2*sr*lam(4)*lam(6) +                                     &
               5112*qzt2*Rqrd*sr*lam(4)*lam(6) +                                &
               2556*qzt2*Rqrd**2*sr*lam(4)*lam(6) + 3798*lam(5)*lam(6) -        &
               7596*qzt2*sr*lam(5)*lam(6) -                                     &
               15192*qzt2*Rqrd*sr*lam(5)*lam(6) -                               &
               7596*qzt2*Rqrd**2*sr*lam(5)*lam(6) - 14265*lam(6)**2 +           &
               28530*qzt2*sr*lam(6)**2 + 57060*qzt2*Rqrd*sr*lam(6)**2 +         &
               28530*qzt2*Rqrd**2*sr*lam(6)**2 - 648*lam(1)*lam(7) +            &
               1296*qzt2*sr*lam(1)*lam(7) +                                     &
               2592*qzt2*Rqrd*sr*lam(1)*lam(7) +                                &
               1296*qzt2*Rqrd**2*sr*lam(1)*lam(7) + 1296*lam(2)*lam(7) -        &
               2592*qzt2*sr*lam(2)*lam(7) -                                     &
               5184*qzt2*Rqrd*sr*lam(2)*lam(7) -                                &
               2592*qzt2*Rqrd**2*sr*lam(2)*lam(7) - 828*lam(3)*lam(7) +         &
               1656*qzt2*sr*lam(3)*lam(7) +                                     &
               3312*qzt2*Rqrd*sr*lam(3)*lam(7) +                                &
               1656*qzt2*Rqrd**2*sr*lam(3)*lam(7) - 90*lam(4)*lam(7) +          &
               180*qzt2*sr*lam(4)*lam(7) +                                      &
               360*qzt2*Rqrd*sr*lam(4)*lam(7) +                                 &
               180*qzt2*Rqrd**2*sr*lam(4)*lam(7) + 234*lam(5)*lam(7) -          &
               468*qzt2*sr*lam(5)*lam(7) -                                      &
               936*qzt2*Rqrd*sr*lam(5)*lam(7) -                                 &
               468*qzt2*Rqrd**2*sr*lam(5)*lam(7) - 3798*lam(6)*lam(7) +         &
               7596*qzt2*sr*lam(6)*lam(7) +                                     &
               15192*qzt2*Rqrd*sr*lam(6)*lam(7) +                               &
               7596*qzt2*Rqrd**2*sr*lam(6)*lam(7) - 117*lam(7)**2 +             &
               234*qzt2*sr*lam(7)**2 + 468*qzt2*Rqrd*sr*lam(7)**2 +             &
               234*qzt2*Rqrd**2*sr*lam(7)**2 - 144*lam(1)*lam(8) +              &
               288*qzt2*sr*lam(1)*lam(8) +                                      &
               576*qzt2*Rqrd*sr*lam(1)*lam(8) +                                 &
               288*qzt2*Rqrd**2*sr*lam(1)*lam(8) - 576*lam(2)*lam(8) +          &
               1152*qzt2*sr*lam(2)*lam(8) +                                     &
               2304*qzt2*Rqrd*sr*lam(2)*lam(8) +                                &
               1152*qzt2*Rqrd**2*sr*lam(2)*lam(8) - 384*lam(3)*lam(8) +         &
               768*qzt2*sr*lam(3)*lam(8) +                                      &
               1536*qzt2*Rqrd*sr*lam(3)*lam(8) +                                &
               768*qzt2*Rqrd**2*sr*lam(3)*lam(8) + 288*lam(8)**2 -              &
               576*qzt2*sr*lam(8)**2 - 1152*qzt2*Rqrd*sr*lam(8)**2 -            &
               576*qzt2*Rqrd**2*sr*lam(8)**2 - 2304*lam(2)*lam(9) +             &
               4608*qzt2*sr*lam(2)*lam(9) +                                     &
               9216*qzt2*Rqrd*sr*lam(2)*lam(9) +                                &
               4608*qzt2*Rqrd**2*sr*lam(2)*lam(9) + 228*lam(4)*lam(9) -         &
               456*qzt2*sr*lam(4)*lam(9) -                                      &
               912*qzt2*Rqrd*sr*lam(4)*lam(9) -                                 &
               456*qzt2*Rqrd**2*sr*lam(4)*lam(9) - 564*lam(5)*lam(9) +          &
               1128*qzt2*sr*lam(5)*lam(9) +                                     &
               2256*qzt2*Rqrd*sr*lam(5)*lam(9) +                                &
               1128*qzt2*Rqrd**2*sr*lam(5)*lam(9) + 6684*lam(6)*lam(9) -        &
               13368*qzt2*sr*lam(6)*lam(9) -                                    &
               26736*qzt2*Rqrd*sr*lam(6)*lam(9) -                               &
               13368*qzt2*Rqrd**2*sr*lam(6)*lam(9) + 564*lam(7)*lam(9) -        &
               1128*qzt2*sr*lam(7)*lam(9) -                                     &
               2256*qzt2*Rqrd*sr*lam(7)*lam(9) -                                &
               1128*qzt2*Rqrd**2*sr*lam(7)*lam(9) - 568*lam(9)**2 +             &
               1136*qzt2*sr*lam(9)**2 + 2272*qzt2*Rqrd*sr*lam(9)**2 +           &
               1136*qzt2*Rqrd**2*sr*lam(9)**2 -                                 &
               12*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                               &
                (384*lam(2) -                                                   &
                  3*(21*lam(4) - 59*lam(5) + 363*lam(6) + 59*lam(7)) +          &
                  286*lam(9))*lam(10) +                                         &
               2616*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(10)**2 +                 &
               Rqsd*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                             &
                (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +               &
                  117*lam(5)**2 + 1278*lam(4)*lam(6) -                          &
                  3798*lam(5)*lam(6) + 14265*lam(6)**2 +                        &
                  90*lam(4)*lam(7) - 234*lam(5)*lam(7) +                        &
                  3798*lam(6)*lam(7) + 117*lam(7)**2 - 288*lam(8)**2 +          &
                  72*lam(1)*(9*(lam(4) + lam(5) + lam(6) + lam(7)) +            &
                     2*lam(8)) +                                                &
                  12*lam(3)*(21*lam(4) + 69*lam(5) + 21*lam(6) +                &
                     69*lam(7) + 32*lam(8)) - 228*lam(4)*lam(9) +               &
                  564*lam(5)*lam(9) - 6684*lam(6)*lam(9) -                      &
                  564*lam(7)*lam(9) + 568*lam(9)**2 -                           &
                  144*lam(2)*                                                   &
                   (3*lam(4) - 9*lam(5) + 69*lam(6) + 9*lam(7) -                &
                     4*(lam(8) + 4*lam(9) - 8*lam(10))) +                       &
                  12*(63*lam(4) - 177*lam(5) + 1089*lam(6) +                    &
                     177*lam(7) - 286*lam(9))*lam(10) + 2616*lam(10)**2)        &
+ 4*m2*(1 + Rqrd)*sr*(864*lam(2)**2 + 45*lam(4)**2 + 297*lam(5)**2 -            &
                  1782*lam(5)*lam(6) + 7173*lam(6)**2 -                         &
                  594*lam(5)*lam(7) + 1782*lam(6)*lam(7) +                      &
                  297*lam(7)**2 + 180*lam(5)*lam(8) + 36*lam(6)*lam(8) -        &
                  180*lam(7)*lam(8) + 576*lam(5)*lam(9) -                       &
                  3264*lam(6)*lam(9) - 576*lam(7)*lam(9) +                      &
                  240*lam(8)*lam(9) + 608*lam(9)**2 -                           &
                  6*lam(4)*(27*lam(5) - 93*lam(6) - 27*lam(7) +                 &
                     6*lam(8) + 32*lam(9) - 64*lam(10)) -                       &
                  8*(144*lam(5) -                                               &
                     3*(272*lam(6) + 48*lam(7) + 9*lam(8)) + 280*lam(9)         &
)*lam(10) + 2240*lam(10)**2 -                                                   &
                  24*lam(2)*(3*lam(4) - 15*lam(5) + 213*lam(6) +                &
                     15*lam(7) - 38*lam(9) + 126*lam(10))))) +                  &
         2*expuISs*qxy2*(1 + Rqsu)*tanhus**2*                                   &
          (dRqrd*quB*qxy2*(1 + Rqrd)**2*                                        &
             (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +                  &
               72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -             &
                  2*lam(8)) +                                                   &
               24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                            &
                  lam(8)*(17*lam(3) + 6*lam(8))) +                              &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                     &
                  10*lam(9) - 12*lam(10))**2) +                                 &
            dRqsu*qdB*(3*qzt2*(1 + Rqrd)*(1 + Rqsu)*                            &
                (48*lam(2)**2 +                                                 &
                  216*lam(2)*(lam(4) - lam(5) + lam(6) - lam(7)) +              &
                  3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +              &
                     lam(6)**2 +                                                &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                 &
                     26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2         &
) + 96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) - 64*lam(8)**2 -          &
                  128*lam(3)*                                                   &
                   (3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) +                 &
                  48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) -                   &
                     3*lam(7) - 16*lam(9) + 16*lam(10))) +                      &
               2*m2*(2592*lam(1)**2 + 576*lam(3)**2 + 27*lam(4)**2 +            &
                  63*lam(5)**2 + 27*lam(6)**2 - 18*lam(6)*lam(7) +              &
                  63*lam(7)**2 + 192*lam(6)*lam(9) + 576*lam(7)*lam(9) +        &
                  320*lam(9)**2 -                                               &
                  6*lam(4)*(3*lam(5) + 9*lam(6) - 3*lam(7) +                    &
                     32*lam(9) - 32*lam(10)) -                                  &
                  64*(3*(lam(6) + lam(7)) + 5*lam(9))*lam(10) +                 &
                  96*lam(10)**2 +                                               &
                  6*lam(5)*(3*lam(6) - 21*lam(7) +                              &
                     32*(-3*lam(9) + lam(10)))))) +                             &
         expuISs*qdB*quB*qzt2*tanhus*                                           &
          (576*dRqrd*lam(2)**2 -                                                &
            dRqrd*(648*lam(1)*lam(4) - 288*m2*sr*lam(2)*lam(4) -                &
               288*m2*Rqrd*sr*lam(2)*lam(4) - 36*lam(3)*lam(4) -                &
               9*lam(4)**2 - 108*m2*sr*lam(4)**2 -                              &
               108*m2*Rqrd*sr*lam(4)**2 - 288*m2*sr*lam(2)*lam(5) -             &
               288*m2*Rqrd*sr*lam(2)*lam(5) + 36*lam(3)*lam(5) -                &
               18*lam(4)*lam(5) - 216*m2*sr*lam(4)*lam(5) -                     &
               216*m2*Rqrd*sr*lam(4)*lam(5) - 9*lam(5)**2 -                     &
               108*m2*sr*lam(5)**2 - 108*m2*Rqrd*sr*lam(5)**2 +                 &
               288*m2*sr*lam(2)*lam(6) + 288*m2*Rqrd*sr*lam(2)*lam(6) -         &
               36*lam(3)*lam(6) + 18*lam(4)*lam(6) +                            &
               216*m2*sr*lam(4)*lam(6) + 216*m2*Rqrd*sr*lam(4)*lam(6) +         &
               18*lam(5)*lam(6) + 216*m2*sr*lam(5)*lam(6) +                     &
               216*m2*Rqrd*sr*lam(5)*lam(6) - 9*lam(6)**2 -                     &
               108*m2*sr*lam(6)**2 - 108*m2*Rqrd*sr*lam(6)**2 +                 &
               288*m2*sr*lam(2)*lam(7) + 288*m2*Rqrd*sr*lam(2)*lam(7) +         &
               36*lam(3)*lam(7) + 18*lam(4)*lam(7) +                            &
               216*m2*sr*lam(4)*lam(7) + 216*m2*Rqrd*sr*lam(4)*lam(7) +         &
               18*lam(5)*lam(7) + 216*m2*sr*lam(5)*lam(7) +                     &
               216*m2*Rqrd*sr*lam(5)*lam(7) - 18*lam(6)*lam(7) -                &
               216*m2*sr*lam(6)*lam(7) - 216*m2*Rqrd*sr*lam(6)*lam(7) -         &
               9*lam(7)**2 - 108*m2*sr*lam(7)**2 -                              &
               108*m2*Rqrd*sr*lam(7)**2 - 576*lam(2)*lam(8) +                   &
               408*lam(3)*lam(8) - 144*m2*sr*lam(4)*lam(8) -                    &
               144*m2*Rqrd*sr*lam(4)*lam(8) - 144*m2*sr*lam(5)*lam(8) -         &
               144*m2*Rqrd*sr*lam(5)*lam(8) + 144*m2*sr*lam(6)*lam(8) +         &
               144*m2*Rqrd*sr*lam(6)*lam(8) + 144*m2*sr*lam(7)*lam(8) +         &
               144*m2*Rqrd*sr*lam(7)*lam(8) - 144*lam(8)**2 -                   &
               72*lam(1)*(9*(lam(5) - lam(6) + lam(7)) + 2*lam(8)) -            &
               960*m2*sr*lam(2)*lam(9) - 960*m2*Rqrd*sr*lam(2)*lam(9) -         &
               60*lam(4)*lam(9) - 60*lam(5)*lam(9) + 60*lam(6)*lam(9) +         &
               60*lam(7)*lam(9) - 480*m2*sr*lam(8)*lam(9) -                     &
               480*m2*Rqrd*sr*lam(8)*lam(9) - 100*lam(9)**2 -                   &
               384*m2*sr*lam(9)**2 - 384*m2*Rqrd*sr*lam(9)**2 +                 &
               2*qzt2*(1 + Rqrd)**2*sr*                                         &
                (36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -                &
                  72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -          &
                     2*lam(8)) +                                                &
                  24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                         &
                     lam(8)*(-17*lam(3) + 6*lam(8))) +                          &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) +                              &
               Rqsu*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                             &
                (36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -                &
                  72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -          &
                     2*lam(8)) +                                                &
                  24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                         &
                     lam(8)*(-17*lam(3) + 6*lam(8))) +                          &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) +                              &
               24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                  &
                  24*m2*(1 + Rqrd)*sr*(2*lam(2) + lam(8)) + 10*lam(9))*         &
                lam(10) - 144*lam(10)**2) +                                     &
            dRqsu*(3*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                 &
                (48*lam(2)**2 +                                                 &
                  216*lam(2)*(lam(4) - lam(5) + lam(6) - lam(7)) +              &
                  3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +              &
                     lam(6)**2 +                                                &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                 &
                     26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2)        &
+ 96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) - 64*lam(8)**2 -            &
                  128*lam(3)*                                                   &
                   (3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) +                 &
                  48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -        &
                     16*lam(9) + 16*lam(10))) +                                 &
               4*m2*(1 + Rqsu)*ss*                                              &
                (2592*lam(1)**2 + 576*lam(3)**2 + 27*lam(4)**2 +                &
                  63*lam(5)**2 + 27*lam(6)**2 - 18*lam(6)*lam(7) +              &
                  63*lam(7)**2 + 192*lam(6)*lam(9) + 576*lam(7)*lam(9) +        &
                  320*lam(9)**2 -                                               &
                  6*lam(4)*(3*lam(5) + 9*lam(6) - 3*lam(7) + 32*lam(9) -        &
                     32*lam(10)) -                                              &
                  64*(3*(lam(6) + lam(7)) + 5*lam(9))*lam(10) +                 &
                  96*lam(10)**2 +                                               &
                  6*lam(5)*(3*lam(6) - 21*lam(7) +                              &
                     32*(-3*lam(9) + lam(10))))))) -                            &
      expuISr*qxy2*tanhur**2*(6*dRqsd*expdISs*quB*qxy2*(1 + Rqru)*              &
          (1 + Rqsd)**2*tanhds**3*                                              &
          (48*lam(2)**2 + 216*lam(2)*                                           &
             (lam(4) - lam(5) + lam(6) - lam(7)) +                              &
            3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                    &
               lam(6)**2 + 2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +           &
               26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +            &
            96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -                 &
            64*lam(8)**2 + 128*lam(3)*                                          &
             (3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) -                       &
            48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -              &
               16*lam(9) + 16*lam(10))) -                                       &
         expdISs*qdB*quB*tanhds**2*                                             &
          (dRqru*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                     &
             (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +                  &
               72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -             &
                  2*lam(8)) +                                                   &
               24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                            &
                  lam(8)*(17*lam(3) + 6*lam(8))) +                              &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                     &
                  10*lam(9) - 12*lam(10))**2) -                                 &
            3*dRqsd*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                  &
             (48*lam(2)**2 + 216*lam(2)*                                        &
                (lam(4) - lam(5) + lam(6) - lam(7)) +                           &
               3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                 &
                  lam(6)**2 +                                                   &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                    &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -              &
               64*lam(8)**2 +                                                   &
               128*lam(3)*(3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) -          &
               48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -           &
                  16*lam(9) + 16*lam(10)))) -                                   &
         2*expdISs*(1 + Rqru)*tanhds*                                           &
          (3*dRqsd*quB*qxy2*(1 + Rqsd)**2*                                      &
             (48*lam(2)**2 +                                                    &
               216*lam(2)*(lam(4) - lam(5) + lam(6) - lam(7)) +                 &
               3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                 &
                  lam(6)**2 +                                                   &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                    &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -              &
               64*lam(8)**2 +                                                   &
               128*lam(3)*(3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) -          &
               48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -           &
                  16*lam(9) + 16*lam(10))) +                                    &
            dRqru*qdB*(qzt2*(1 + Rqru)*(1 + Rqsd)*                              &
                (36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -                &
                  72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -          &
                     2*lam(8)) +                                                &
                  24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                         &
                     lam(8)*(-17*lam(3) + 6*lam(8))) +                          &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) -                              &
               6*m2*(9*lam(4)**2 + 9*lam(5)**2 - 18*lam(5)*lam(6) +             &
                  9*lam(6)**2 - 18*lam(5)*lam(7) + 18*lam(6)*lam(7) +           &
                  9*lam(7)**2 + 12*lam(5)*lam(8) - 12*lam(6)*lam(8) -           &
                  12*lam(7)*lam(8) +                                            &
                  6*lam(4)*(3*lam(5) - 3*(lam(6) + lam(7)) + 2*lam(8)) +        &
                  40*lam(8)*lam(9) + 32*lam(9)**2 +                             &
                  8*lam(2)*(3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +         &
                     10*lam(9) - 12*lam(10)) - 48*lam(8)*lam(10)))) +           &
         qdB*(-3*(1 + Rqru)*(dRqsu*expuISs*(-1 + tanhus)*(1 + tanhus)*          &
                (quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +                           &
                  2*qxy2*(1 + Rqsu)**2*tanhus)*                                 &
                (3*(16*lam(2)**2 +                                              &
                     192*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) +           &
                     64*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) -            &
                     72*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +            &
                     (lam(4) + lam(5) + lam(6) + lam(7))**2) +                  &
                  256*(3*lam(1) + lam(3))*lam(10)) +                            &
               dRqsd*expdISs*quB*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                &
                (48*lam(2)**2 +                                                 &
                  216*lam(2)*(lam(4) - lam(5) + lam(6) - lam(7)) +              &
                  3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +              &
                     lam(6)**2 +                                                &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                 &
                     26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2)        &
+ 96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) - 64*lam(8)**2 +            &
                  128*lam(3)*                                                   &
                   (3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) -                 &
                  48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -        &
                     16*lam(9) + 16*lam(10)))) +                                &
            dRqru*(expdISs*quB*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*       &
                (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +               &
                  72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -          &
                     2*lam(8)) +                                                &
                  24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                         &
                     lam(8)*(17*lam(3) + 6*lam(8))) +                           &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) +                              &
               expuISs*(quB*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*          &
                   (-1 + tanhus**2)*                                            &
                   (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +            &
                     117*lam(5)**2 + 1278*lam(4)*lam(6) -                       &
                     3798*lam(5)*lam(6) + 14265*lam(6)**2 +                     &
                     90*lam(4)*lam(7) - 234*lam(5)*lam(7) +                     &
                     3798*lam(6)*lam(7) + 117*lam(7)**2 +                       &
                     72*lam(1)*                                                 &
                      (-9*(lam(4) + lam(5) + lam(6) + lam(7)) -                 &
                        2*lam(8)) - 288*lam(8)**2 -                             &
                     12*lam(3)*                                                 &
                      (21*lam(4) + 69*lam(5) + 21*lam(6) + 69*lam(7) +          &
                        32*lam(8)) - 228*lam(4)*lam(9) +                        &
                     564*lam(5)*lam(9) - 6684*lam(6)*lam(9) -                   &
                     564*lam(7)*lam(9) + 568*lam(9)**2 -                        &
                     144*lam(2)*                                                &
                      (3*lam(4) - 9*lam(5) + 69*lam(6) + 9*lam(7) -             &
                        4*(lam(8) + 4*lam(9) - 8*lam(10))) +                    &
                     12*(63*lam(4) - 177*lam(5) + 1089*lam(6) +                 &
                        177*lam(7) - 286*lam(9))*lam(10) +                      &
                     2616*lam(10)**2) +                                         &
                  2*(1 + Rqru)*tanhus*                                          &
                   (qzt2*(1 + Rqru)*(1 + Rqsu)*                                 &
                      (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +         &
                        117*lam(5)**2 + 1278*lam(4)*lam(6) -                    &
                        3798*lam(5)*lam(6) + 14265*lam(6)**2 +                  &
                        90*lam(4)*lam(7) - 234*lam(5)*lam(7) +                  &
                        3798*lam(6)*lam(7) + 117*lam(7)**2 -                    &
                        288*lam(8)**2 +                                         &
                        72*lam(1)*                                              &
                         (9*(lam(4) + lam(5) + lam(6) + lam(7)) +               &
                         2*lam(8)) +                                            &
                        12*lam(3)*                                              &
                         (21*lam(4) + 69*lam(5) + 21*lam(6) +                   &
                         69*lam(7) + 32*lam(8)) - 228*lam(4)*lam(9) +           &
                        564*lam(5)*lam(9) - 6684*lam(6)*lam(9) -                &
                        564*lam(7)*lam(9) + 568*lam(9)**2 -                     &
                        144*lam(2)*                                             &
                         (3*lam(4) - 9*lam(5) + 69*lam(6) + 9*lam(7) -          &
                         4*(lam(8) + 4*lam(9) - 8*lam(10))) +                   &
                        12*(63*lam(4) - 177*lam(5) + 1089*lam(6) +              &
                        177*lam(7) - 286*lam(9))*lam(10) +                      &
                        2616*lam(10)**2) +                                      &
                     2*m2*(864*lam(2)**2 + 45*lam(4)**2 + 297*lam(5)**2 -       &
                        1782*lam(5)*lam(6) + 7173*lam(6)**2 -                   &
                        594*lam(5)*lam(7) + 1782*lam(6)*lam(7) +                &
                        297*lam(7)**2 + 180*lam(5)*lam(8) +                     &
                        36*lam(6)*lam(8) - 180*lam(7)*lam(8) +                  &
                        576*lam(5)*lam(9) - 3264*lam(6)*lam(9) -                &
                        576*lam(7)*lam(9) + 240*lam(8)*lam(9) +                 &
                        608*lam(9)**2 -                                         &
                        6*lam(4)*                                               &
                         (27*lam(5) - 93*lam(6) - 27*lam(7) + 6*lam(8) +        &
                         32*lam(9) - 64*lam(10)) -                              &
                        8*(144*lam(5) -                                         &
                         3*(272*lam(6) + 48*lam(7) + 9*lam(8)) +                &
                         280*lam(9))*lam(10) + 2240*lam(10)**2 -                &
                        24*lam(2)*                                              &
                         (3*lam(4) - 15*lam(5) + 213*lam(6) + 15*lam(7) -       &
                         38*lam(9) + 126*lam(10)))))))) -                       &
      expuISr*tanhur*(2*expdISs*qxy2*(1 + Rqsd)*tanhds**2*                      &
          (dRqru*qdB*qxy2*(1 + Rqru)**2*                                        &
             (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +                  &
               72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -             &
                  2*lam(8)) +                                                   &
               24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                            &
                  lam(8)*(17*lam(3) + 6*lam(8))) +                              &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -         &
                  12*lam(10))**2) +                                             &
            dRqsd*quB*(3*qzt2*(1 + Rqru)*(1 + Rqsd)*                            &
                (48*lam(2)**2 +                                                 &
                  216*lam(2)*(lam(4) - lam(5) + lam(6) - lam(7)) +              &
                  3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +              &
                     lam(6)**2 +                                                &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                 &
                     26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2)        &
+ 96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) - 64*lam(8)**2 -            &
                  128*lam(3)*                                                   &
                   (3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) +                 &
                  48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -        &
                     16*lam(9) + 16*lam(10))) +                                 &
               2*m2*(2592*lam(1)**2 + 576*lam(3)**2 + 27*lam(4)**2 +            &
                  63*lam(5)**2 + 27*lam(6)**2 - 18*lam(6)*lam(7) +              &
                  63*lam(7)**2 + 192*lam(6)*lam(9) + 576*lam(7)*lam(9) +        &
                  320*lam(9)**2 -                                               &
                  6*lam(4)*(3*lam(5) + 9*lam(6) - 3*lam(7) + 32*lam(9) -        &
                     32*lam(10)) -                                              &
                  64*(3*(lam(6) + lam(7)) + 5*lam(9))*lam(10) +                 &
                  96*lam(10)**2 +                                               &
                  6*lam(5)*(3*lam(6) - 21*lam(7) +                              &
                     32*(-3*lam(9) + lam(10)))))) +                             &
         expdISs*qdB*quB*qzt2*tanhds*                                           &
          (576*dRqru*lam(2)**2 -                                                &
            dRqru*(648*lam(1)*lam(4) - 288*m2*sr*lam(2)*lam(4) -                &
               288*m2*Rqru*sr*lam(2)*lam(4) - 36*lam(3)*lam(4) -                &
               9*lam(4)**2 - 108*m2*sr*lam(4)**2 -                              &
               108*m2*Rqru*sr*lam(4)**2 - 288*m2*sr*lam(2)*lam(5) -             &
               288*m2*Rqru*sr*lam(2)*lam(5) + 36*lam(3)*lam(5) -                &
               18*lam(4)*lam(5) - 216*m2*sr*lam(4)*lam(5) -                     &
               216*m2*Rqru*sr*lam(4)*lam(5) - 9*lam(5)**2 -                     &
               108*m2*sr*lam(5)**2 - 108*m2*Rqru*sr*lam(5)**2 +                 &
               288*m2*sr*lam(2)*lam(6) + 288*m2*Rqru*sr*lam(2)*lam(6) -         &
               36*lam(3)*lam(6) + 18*lam(4)*lam(6) +                            &
               216*m2*sr*lam(4)*lam(6) + 216*m2*Rqru*sr*lam(4)*lam(6) +         &
               18*lam(5)*lam(6) + 216*m2*sr*lam(5)*lam(6) +                     &
               216*m2*Rqru*sr*lam(5)*lam(6) - 9*lam(6)**2 -                     &
               108*m2*sr*lam(6)**2 - 108*m2*Rqru*sr*lam(6)**2 +                 &
               288*m2*sr*lam(2)*lam(7) + 288*m2*Rqru*sr*lam(2)*lam(7) +         &
               36*lam(3)*lam(7) + 18*lam(4)*lam(7) +                            &
               216*m2*sr*lam(4)*lam(7) + 216*m2*Rqru*sr*lam(4)*lam(7) +         &
               18*lam(5)*lam(7) + 216*m2*sr*lam(5)*lam(7) +                     &
               216*m2*Rqru*sr*lam(5)*lam(7) - 18*lam(6)*lam(7) -                &
               216*m2*sr*lam(6)*lam(7) - 216*m2*Rqru*sr*lam(6)*lam(7) -         &
               9*lam(7)**2 - 108*m2*sr*lam(7)**2 -                              &
               108*m2*Rqru*sr*lam(7)**2 - 576*lam(2)*lam(8) +                   &
               408*lam(3)*lam(8) - 144*m2*sr*lam(4)*lam(8) -                    &
               144*m2*Rqru*sr*lam(4)*lam(8) - 144*m2*sr*lam(5)*lam(8) -         &
               144*m2*Rqru*sr*lam(5)*lam(8) + 144*m2*sr*lam(6)*lam(8) +         &
               144*m2*Rqru*sr*lam(6)*lam(8) + 144*m2*sr*lam(7)*lam(8) +         &
               144*m2*Rqru*sr*lam(7)*lam(8) - 144*lam(8)**2 -                   &
               72*lam(1)*(9*(lam(5) - lam(6) + lam(7)) + 2*lam(8)) -            &
               960*m2*sr*lam(2)*lam(9) - 960*m2*Rqru*sr*lam(2)*lam(9) -         &
               60*lam(4)*lam(9) - 60*lam(5)*lam(9) + 60*lam(6)*lam(9) +         &
               60*lam(7)*lam(9) - 480*m2*sr*lam(8)*lam(9) -                     &
               480*m2*Rqru*sr*lam(8)*lam(9) - 100*lam(9)**2 -                   &
               384*m2*sr*lam(9)**2 - 384*m2*Rqru*sr*lam(9)**2 +                 &
               2*qzt2*(1 + Rqru)**2*sr*                                         &
                (36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -                &
                  72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -          &
                     2*lam(8)) +                                                &
                  24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                         &
                     lam(8)*(-17*lam(3) + 6*lam(8))) +                          &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) +                              &
               Rqsd*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                             &
                (36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -                &
                  72*lam(1)*(9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -          &
                     2*lam(8)) +                                                &
                  24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                         &
                     lam(8)*(-17*lam(3) + 6*lam(8))) +                          &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) +                              &
               24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                  &
                  24*m2*(1 + Rqru)*sr*(2*lam(2) + lam(8)) + 10*lam(9))*         &
                lam(10) - 144*lam(10)**2) +                                     &
            dRqsd*(3*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                 &
                (48*lam(2)**2 +                                                 &
                  216*lam(2)*(lam(4) - lam(5) + lam(6) - lam(7)) +              &
                  3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +              &
                     lam(6)**2 +                                                &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                 &
                     26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2)        &
+ 96*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) - 64*lam(8)**2 -            &
                  128*lam(3)*                                                   &
                   (3*lam(5) - 3*lam(7) - 4*lam(9) + lam(10)) +                 &
                  48*lam(1)*(9*lam(4) + 3*lam(5) - 9*lam(6) - 3*lam(7) -        &
                     16*lam(9) + 16*lam(10))) +                                 &
               4*m2*(1 + Rqsd)*ss*                                              &
                (2592*lam(1)**2 + 576*lam(3)**2 + 27*lam(4)**2 +                &
                  63*lam(5)**2 + 27*lam(6)**2 - 18*lam(6)*lam(7) +              &
                  63*lam(7)**2 + 192*lam(6)*lam(9) + 576*lam(7)*lam(9) +        &
                  320*lam(9)**2 -                                               &
                  6*lam(4)*(3*lam(5) + 9*lam(6) - 3*lam(7) + 32*lam(9) -        &
                     32*lam(10)) -                                              &
                  64*(3*(lam(6) + lam(7)) + 5*lam(9))*lam(10) +                 &
                  96*lam(10)**2 +                                               &
                  6*lam(5)*(3*lam(6) - 21*lam(7) +                              &
                     32*(-3*lam(9) + lam(10)))))) +                             &
         qdB*(dRqsu*expuISs*tanhus*                                             &
             (-3*qzt2*(1 + Rqru)*                                               &
                (quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +                           &
                  2*qxy2*(1 + Rqsu)**2*tanhus)*                                 &
                (3*(16*lam(2)**2 -                                              &
                     192*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) -           &
                     64*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) -            &
                     72*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +            &
                     (lam(4) + lam(5) + lam(6) + lam(7))**2) -                  &
                  256*(3*lam(1) + lam(3))*lam(10)) +                            &
               4*m2*(1 + Rqsu)*(quB*qzt2*ss + qxy2*tanhus)*                     &
                (384*(3*lam(1) + lam(3))**2 +                                   &
                  45*(lam(4) - lam(5) - lam(6) + lam(7))**2 +                   &
                  192*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +             &
                  32*lam(10)**2)) +                                             &
            dRqru*(-2*expdISs*qxy2*(1 + Rqru)*                                  &
                ((qxy2 + qzt2)*(1 + Rqru)*(1 + Rqsd)*                           &
                   (-36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) +            &
                     72*lam(1)*                                                 &
                      (9*lam(4) - 9*(lam(5) - lam(6) + lam(7)) -                &
                        2*lam(8)) +                                             &
                     24*(24*lam(2)**2 + 24*lam(2)*lam(8) +                      &
                        lam(8)*(17*lam(3) + 6*lam(8))) +                        &
                     (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +               &
                        10*lam(9) - 12*lam(10))**2) +                           &
                  24*m2*(-108*lam(1)**2 + 12*lam(1)*lam(3) -                    &
                     17*lam(3)**2 + 27*lam(4)*(lam(5) - lam(6)) -               &
                     27*lam(5)*lam(7) + 27*lam(6)*lam(7) -                      &
                     6*lam(5)*lam(8) + 6*lam(6)*lam(8) - 17*lam(8)**2 -         &
                     10*lam(8)*lam(9) + 24*lam(9)**2 +                          &
                     12*lam(8)*lam(10) +                                        &
                     lam(2)*(-6*lam(4) - 6*lam(5) + 6*lam(6) + 6*lam(7) -       &
                        20*lam(9) + 24*lam(10)))) +                             &
               expuISs*(-((1 + Rqsu)*                                           &
                     (-2*qxy2*qzt2*(1 + Rqru)**2*                               &
                        (576*lam(2)**2 + 9*lam(4)**2 -                          &
                         90*lam(4)*lam(5) + 117*lam(5)**2 +                     &
                         1278*lam(4)*lam(6) - 3798*lam(5)*lam(6) +              &
                         14265*lam(6)**2 + 90*lam(4)*lam(7) -                   &
                         234*lam(5)*lam(7) + 3798*lam(6)*lam(7) +               &
                         117*lam(7)**2 +                                        &
                         72*lam(1)*                                             &
                         (-9*(lam(4) + lam(5) + lam(6) + lam(7)) -              &
                         2*lam(8)) - 288*lam(8)**2 -                            &
                         12*lam(3)*                                             &
                         (21*lam(4) + 69*lam(5) + 21*lam(6) +                   &
                         69*lam(7) + 32*lam(8)) - 228*lam(4)*lam(9) +           &
                         564*lam(5)*lam(9) - 6684*lam(6)*lam(9) -               &
                         564*lam(7)*lam(9) + 568*lam(9)**2 -                    &
                         144*lam(2)*                                            &
                         (3*lam(4) - 9*lam(5) + 69*lam(6) + 9*lam(7) -          &
                         4*(lam(8) + 4*lam(9) - 8*lam(10))) +                   &
                         12*(63*lam(4) - 177*lam(5) + 1089*lam(6) +             &
                        177*lam(7) - 286*lam(9))*lam(10) +                      &
                         2616*lam(10)**2) +                                     &
                       2*qxy2**2*(1 + Rqru)**2*(-1 + tanhus**2)*                &
                        (576*lam(2)**2 + 9*lam(4)**2 -                          &
                         90*lam(4)*lam(5) + 117*lam(5)**2 +                     &
                         1278*lam(4)*lam(6) - 3798*lam(5)*lam(6) +              &
                         14265*lam(6)**2 + 90*lam(4)*lam(7) -                   &
                         234*lam(5)*lam(7) + 3798*lam(6)*lam(7) +               &
                         117*lam(7)**2 +                                        &
                         72*lam(1)*                                             &
                         (-9*(lam(4) + lam(5) + lam(6) + lam(7)) -              &
                         2*lam(8)) - 288*lam(8)**2 -                            &
                         12*lam(3)*                                             &
                         (21*lam(4) + 69*lam(5) + 21*lam(6) +                   &
                         69*lam(7) + 32*lam(8)) - 228*lam(4)*lam(9) +           &
                         564*lam(5)*lam(9) - 6684*lam(6)*lam(9) -               &
                         564*lam(7)*lam(9) + 568*lam(9)**2 -                    &
                         144*lam(2)*                                            &
                         (3*lam(4) - 9*lam(5) + 69*lam(6) + 9*lam(7) -          &
                         4*(lam(8) + 4*lam(9) - 8*lam(10))) +                   &
                         12*(63*lam(4) - 177*lam(5) + 1089*lam(6) +             &
                        177*lam(7) - 286*lam(9))*lam(10) +                      &
                         2616*lam(10)**2) -                                     &
                       quB*qzt2*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*tanhus*          &
                        (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +       &
                         117*lam(5)**2 + 1278*lam(4)*lam(6) -                   &
                         3798*lam(5)*lam(6) + 14265*lam(6)**2 +                 &
                         90*lam(4)*lam(7) - 234*lam(5)*lam(7) +                 &
                         3798*lam(6)*lam(7) + 117*lam(7)**2 -                   &
                         288*lam(8)**2 +                                        &
                         72*lam(1)*                                             &
                         (9*(lam(4) + lam(5) + lam(6) + lam(7)) +               &
                         2*lam(8)) +                                            &
                         12*lam(3)*                                             &
                         (21*lam(4) + 69*lam(5) + 21*lam(6) +                   &
                         69*lam(7) + 32*lam(8)) - 228*lam(4)*lam(9) +           &
                         564*lam(5)*lam(9) - 6684*lam(6)*lam(9) -               &
                         564*lam(7)*lam(9) + 568*lam(9)**2 -                    &
                         144*lam(2)*                                            &
                         (3*lam(4) - 9*lam(5) + 69*lam(6) + 9*lam(7) -          &
                         4*(lam(8) + 4*lam(9) - 8*lam(10))) +                   &
                         12*(63*lam(4) - 177*lam(5) + 1089*lam(6) +             &
                         177*lam(7) - 286*lam(9))*lam(10) +                     &
                         2616*lam(10)**2))) +                                   &
                  4*m2*(1 + Rqru)*                                              &
                   (quB*qzt2*sr*tanhus*                                         &
                      (864*lam(2)**2 + 45*lam(4)**2 + 297*lam(5)**2 -           &
                        1782*lam(5)*lam(6) + 7173*lam(6)**2 -                   &
                        594*lam(5)*lam(7) + 1782*lam(6)*lam(7) +                &
                        297*lam(7)**2 + 180*lam(5)*lam(8) +                     &
                        36*lam(6)*lam(8) - 180*lam(7)*lam(8) +                  &
                        576*lam(5)*lam(9) - 3264*lam(6)*lam(9) -                &
                        576*lam(7)*lam(9) + 240*lam(8)*lam(9) +                 &
                        608*lam(9)**2 -                                         &
                        6*lam(4)*                                               &
                         (27*lam(5) - 93*lam(6) - 27*lam(7) + 6*lam(8) +        &
                         32*lam(9) - 64*lam(10)) -                              &
                        8*(144*lam(5) -                                         &
                         3*(272*lam(6) + 48*lam(7) + 9*lam(8)) +                &
                         280*lam(9))*lam(10) + 2240*lam(10)**2 -                &
                        24*lam(2)*                                              &
                         (3*lam(4) - 15*lam(5) + 213*lam(6) + 15*lam(7) -       &
                         38*lam(9) + 126*lam(10))) +                            &
                     12*qxy2*(108*lam(1)**2 + 72*lam(2)**2 +                    &
                        84*lam(1)*lam(3) + 46*lam(3)**2 -                       &
                        162*lam(5)*lam(6) + 594*lam(6)**2 +                     &
                        27*lam(4)*(lam(5) + 3*lam(6)) +                         &
                        135*lam(5)*lam(7) + 189*lam(6)*lam(7) +                 &
                        30*lam(5)*lam(8) + 6*lam(6)*lam(8) +                    &
                        10*lam(8)**2 - 288*lam(6)*lam(9) +                      &
                        20*lam(8)*lam(9) - 24*lam(9)**2 +                       &
                        6*(96*lam(6) + 3*lam(8) + 8*lam(9))*lam(10) -           &
                        48*lam(10)**2 -                                         &
                        2*lam(2)*                                               &
                         (3*lam(4) - 15*lam(5) + 213*lam(6) + 15*lam(7) -       &
                          38*lam(9) + 126*lam(10)))))))))/1296.
  dtLAM(7)=(-3*dRqsd*expdISr*expdISs*qdB*quB*                                &
       ((qxy2 + qzt2)*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                &
          (144*lam(2)**2 + 576*lam(1)*                                          &
             (lam(4) - lam(5) - lam(6) + lam(7)) +                              &
            192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +                    &
            216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                    &
            9*(lam(4) + lam(5) + lam(6) + lam(7))**2 +                          &
            256*(3*lam(1) + lam(3))*lam(10)) +                                  &
         8*m2*qzt2*(1 + Rqsd)*ss*                                               &
          (216*lam(2)**2 + 64*(3*lam(1) + lam(3))**2 +                          &
            36*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                     &
            27*((lam(5) + lam(6))**2 + (lam(4) + lam(7))**2) +                  &
            48*lam(10)**2)) - 3*dRqsu*expuISr*expuISs*qdB*quB*                  &
       ((qxy2 + qzt2)*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                &
          (144*lam(2)**2 + 576*lam(1)*                                          &
             (lam(4) - lam(5) - lam(6) + lam(7)) +                              &
            192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +                    &
            216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                    &
            9*(lam(4) + lam(5) + lam(6) + lam(7))**2 +                          &
            256*(3*lam(1) + lam(3))*lam(10)) +                                  &
         8*m2*qzt2*(1 + Rqsu)*ss*                                               &
          (216*lam(2)**2 + 64*(3*lam(1) + lam(3))**2 +                          &
            36*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                     &
            27*((lam(5) + lam(6))**2 + (lam(4) + lam(7))**2) +                  &
            48*lam(10)**2)) + 6*dRqrd*expdISr*quB*qxy2**2*(1 + Rqrd)**2*        &
       tanhdr**3*(-(expuISs*(1 + Rqsu)*(-1 + tanhus**2)*                        &
            (576*lam(2)**2 + 9*lam(4)**2 + 18*lam(4)*lam(5) +                   &
              9*lam(5)**2 - 18*lam(4)*lam(6) - 18*lam(5)*lam(6) +               &
              9*lam(6)**2 - 882*lam(4)*lam(7) - 882*lam(5)*lam(7) +             &
              882*lam(6)*lam(7) + 11241*lam(7)**2 - 144*lam(8)**2 +             &
              72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +              &
                 2*lam(8)) + 36*lam(3)*                                         &
               (lam(4) - lam(5) + lam(6) - lam(7) + 10*lam(8)) +                &
              132*lam(4)*lam(9) + 132*lam(5)*lam(9) -                           &
              132*lam(6)*lam(9) - 4740*lam(7)*lam(9) + 412*lam(9)**2 +          &
              96*lam(2)*(3*lam(4) + 3*lam(5) -                                  &
                 3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) +               &
              24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                   &
                 10*lam(9))*lam(10) - 144*lam(10)**2)) +                        &
         expdISs*(1 + Rqsd)*(-1 + tanhds**2)*                                   &
          (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +                     &
            117*lam(5)**2 - 18*lam(4)*lam(6) + 90*lam(5)*lam(6) +               &
            9*lam(6)**2 + 522*lam(4)*lam(7) + 198*lam(5)*lam(7) -               &
            522*lam(6)*lam(7) - 5499*lam(7)**2 +                                &
            36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +           &
            288*lam(8)**2 + 72*lam(1)*                                          &
             (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) + 2*lam(8)) -           &
            156*lam(4)*lam(9) + 204*lam(5)*lam(9) + 156*lam(6)*lam(9) +         &
            2100*lam(7)*lam(9) - 56*lam(9)**2 -                                 &
            48*lam(2)*(3*(lam(4) + lam(5) - lam(6) - 25*lam(7) +                &
                  4*lam(8)) + 16*lam(9)) -                                      &
            36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) - 10*lam(9))*       &
             lam(10) + 456*lam(10)**2)) -                                       &
      6*dRqru*expuISr*qdB*qxy2**2*(1 + Rqru)**2*tanhur**3*                      &
       (expdISs*(1 + Rqsd)*(-1 + tanhds**2)*                                    &
          (576*lam(2)**2 + 9*lam(4)**2 + 18*lam(4)*lam(5) +                     &
            9*lam(5)**2 - 18*lam(4)*lam(6) - 18*lam(5)*lam(6) +                 &
            9*lam(6)**2 - 882*lam(4)*lam(7) - 882*lam(5)*lam(7) +               &
            882*lam(6)*lam(7) + 11241*lam(7)**2 - 144*lam(8)**2 +               &
            72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +                &
               2*lam(8)) + 36*lam(3)*                                           &
             (lam(4) - lam(5) + lam(6) - lam(7) + 10*lam(8)) +                  &
            132*lam(4)*lam(9) + 132*lam(5)*lam(9) - 132*lam(6)*lam(9) -         &
            4740*lam(7)*lam(9) + 412*lam(9)**2 +                                &
            96*lam(2)*(3*lam(4) + 3*lam(5) -                                    &
               3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) +                 &
            24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) + 10*lam(9))*         &
             lam(10) - 144*lam(10)**2) -                                        &
         expuISs*(1 + Rqsu)*(-1 + tanhus**2)*                                   &
          (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +                     &
            117*lam(5)**2 - 18*lam(4)*lam(6) + 90*lam(5)*lam(6) +               &
            9*lam(6)**2 + 522*lam(4)*lam(7) + 198*lam(5)*lam(7) -               &
            522*lam(6)*lam(7) - 5499*lam(7)**2 +                                &
            36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +           &
            288*lam(8)**2 + 72*lam(1)*                                          &
             (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) + 2*lam(8)) -           &
            156*lam(4)*lam(9) + 204*lam(5)*lam(9) + 156*lam(6)*lam(9) +         &
            2100*lam(7)*lam(9) - 56*lam(9)**2 -                                 &
            48*lam(2)*(3*(lam(4) + lam(5) - lam(6) - 25*lam(7) +                &
                  4*lam(8)) + 16*lam(9)) -                                      &
            36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) - 10*lam(9))*       &
             lam(10) + 456*lam(10)**2)) -                                       &
      3*dRqrd*expdISr*qdB*quB*                                                  &
       (expuISs*(qxy2 + qzt2)*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*        &
          (576*lam(2)**2 + 9*lam(4)**2 + 18*lam(4)*lam(5) +                     &
            9*lam(5)**2 - 18*lam(4)*lam(6) - 18*lam(5)*lam(6) +                 &
            9*lam(6)**2 - 882*lam(4)*lam(7) - 882*lam(5)*lam(7) +               &
            882*lam(6)*lam(7) + 11241*lam(7)**2 - 144*lam(8)**2 +               &
            72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +                &
               2*lam(8)) + 36*lam(3)*                                           &
             (lam(4) - lam(5) + lam(6) - lam(7) + 10*lam(8)) +                  &
            132*lam(4)*lam(9) + 132*lam(5)*lam(9) - 132*lam(6)*lam(9) -         &
            4740*lam(7)*lam(9) + 412*lam(9)**2 +                                &
            96*lam(2)*(3*lam(4) + 3*lam(5) -                                    &
               3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) +                 &
            24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) + 10*lam(9))*         &
             lam(10) - 144*lam(10)**2) +                                        &
         24*expdISs*m2*qzt2*(1 + Rqrd)*sr*                                      &
          (72*lam(1)**2 - 48*lam(2)**2 + 168*lam(1)*lam(3) +                    &
            60*lam(3)**2 + 9*lam(4)**2 + 117*lam(5)**2 +                        &
            90*lam(5)*lam(6) + 9*lam(6)**2 + 126*lam(4)*lam(7) +                &
            36*lam(5)*lam(7) - 36*lam(6)*lam(7) - 351*lam(7)**2 +               &
            60*lam(5)*lam(8) + 12*lam(6)*lam(8) - 12*lam(8)**2 +                &
            192*lam(7)*lam(9) + 40*lam(8)*lam(9) + 16*lam(9)**2 +               &
            4*lam(2)*(3*lam(4) - 15*lam(5) - 3*lam(6) + 87*lam(7) -             &
               26*lam(9) - 30*lam(10)) +                                        &
            12*(3*lam(8) - 8*lam(9))*lam(10) + 96*lam(10)**2) -                 &
         expdISs*(qxy2 + qzt2)*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*       &
          (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +                     &
            117*lam(5)**2 - 18*lam(4)*lam(6) + 90*lam(5)*lam(6) +               &
            9*lam(6)**2 + 522*lam(4)*lam(7) + 198*lam(5)*lam(7) -               &
            522*lam(6)*lam(7) - 5499*lam(7)**2 +                                &
            36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +           &
            288*lam(8)**2 + 72*lam(1)*                                          &
             (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) + 2*lam(8)) -           &
            156*lam(4)*lam(9) + 204*lam(5)*lam(9) + 156*lam(6)*lam(9) +         &
            2100*lam(7)*lam(9) - 56*lam(9)**2 -                                 &
            48*lam(2)*(3*(lam(4) + lam(5) - lam(6) - 25*lam(7) +                &
                  4*lam(8)) + 16*lam(9)) -                                      &
            36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -                   &
               10*lam(9))*lam(10) + 456*lam(10)**2) -                           &
         24*expuISs*m2*qzt2*(1 + Rqrd)*sr*                                      &
          (72*lam(1)**2 + 96*lam(2)**2 + 24*lam(1)*lam(3) +                     &
            30*lam(3)**2 + 9*lam(4)**2 + 9*(lam(5) - lam(6))**2 -               &
            90*lam(4)*lam(7) - 72*lam(5)*lam(7) + 72*lam(6)*lam(7) +            &
            945*lam(7)**2 - 12*lam(5)*lam(8) + 12*lam(6)*lam(8) +               &
            30*lam(8)**2 - 384*lam(7)*lam(9) - 20*lam(8)*lam(9) +               &
            16*lam(9)**2 + 24*lam(8)*lam(10) +                                  &
            4*lam(2)*(3*lam(4) + 3*lam(5) - 3*lam(6) - 147*lam(7) +             &
               22*lam(9) + 12*lam(10)))) +                                      &
      3*dRqru*expuISr*qdB*quB*                                                  &
       (-(expdISs*(qxy2 + qzt2)*(1 + Rqsd)*                                     &
            (-1 + 2*qzt2*(1 + Rqru)**2*sr)*                                     &
            (576*lam(2)**2 + 9*lam(4)**2 + 18*lam(4)*lam(5) +                   &
              9*lam(5)**2 - 18*lam(4)*lam(6) - 18*lam(5)*lam(6) +               &
              9*lam(6)**2 - 882*lam(4)*lam(7) - 882*lam(5)*lam(7) +             &
              882*lam(6)*lam(7) + 11241*lam(7)**2 - 144*lam(8)**2 +             &
              72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +              &
                 2*lam(8)) + 36*lam(3)*                                         &
               (lam(4) - lam(5) + lam(6) - lam(7) + 10*lam(8)) +                &
              132*lam(4)*lam(9) + 132*lam(5)*lam(9) -                           &
              132*lam(6)*lam(9) - 4740*lam(7)*lam(9) + 412*lam(9)**2 +          &
              96*lam(2)*(3*lam(4) + 3*lam(5) -                                  &
                 3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) +               &
              24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                   &
                 10*lam(9))*lam(10) - 144*lam(10)**2)) -                        &
         24*expuISs*m2*qzt2*(1 + Rqru)*sr*                                      &
          (72*lam(1)**2 - 48*lam(2)**2 + 168*lam(1)*lam(3) +                    &
            60*lam(3)**2 + 9*lam(4)**2 + 117*lam(5)**2 +                        &
            90*lam(5)*lam(6) + 9*lam(6)**2 + 126*lam(4)*lam(7) +                &
            36*lam(5)*lam(7) - 36*lam(6)*lam(7) - 351*lam(7)**2 +               &
            60*lam(5)*lam(8) + 12*lam(6)*lam(8) - 12*lam(8)**2 +                &
            192*lam(7)*lam(9) + 40*lam(8)*lam(9) + 16*lam(9)**2 +               &
            4*lam(2)*(3*lam(4) - 15*lam(5) - 3*lam(6) + 87*lam(7) -             &
               26*lam(9) - 30*lam(10)) +                                        &
            12*(3*lam(8) - 8*lam(9))*lam(10) + 96*lam(10)**2) +                 &
         expuISs*(qxy2 + qzt2)*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*       &
          (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +                     &
            117*lam(5)**2 - 18*lam(4)*lam(6) + 90*lam(5)*lam(6) +               &
            9*lam(6)**2 + 522*lam(4)*lam(7) + 198*lam(5)*lam(7) -               &
            522*lam(6)*lam(7) - 5499*lam(7)**2 +                                &
            36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +           &
            288*lam(8)**2 + 72*lam(1)*                                          &
             (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) + 2*lam(8)) -           &
            156*lam(4)*lam(9) + 204*lam(5)*lam(9) + 156*lam(6)*lam(9) +         &
            2100*lam(7)*lam(9) - 56*lam(9)**2 -                                 &
            48*lam(2)*(3*(lam(4) + lam(5) - lam(6) - 25*lam(7) +                &
                  4*lam(8)) + 16*lam(9)) -                                      &
            36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -                   &
               10*lam(9))*lam(10) + 456*lam(10)**2) +                           &
         24*expdISs*m2*qzt2*(1 + Rqru)*sr*                                      &
          (72*lam(1)**2 + 96*lam(2)**2 + 24*lam(1)*lam(3) +                     &
            30*lam(3)**2 + 9*lam(4)**2 + 9*(lam(5) - lam(6))**2 -               &
            90*lam(4)*lam(7) - 72*lam(5)*lam(7) + 72*lam(6)*lam(7) +            &
            945*lam(7)**2 - 12*lam(5)*lam(8) + 12*lam(6)*lam(8) +               &
            30*lam(8)**2 - 384*lam(7)*lam(9) - 20*lam(8)*lam(9) +               &
            16*lam(9)**2 + 24*lam(8)*lam(10) +                                  &
            4*lam(2)*(3*lam(4) + 3*lam(5) - 3*lam(6) - 147*lam(7) +             &
               22*lam(9) + 12*lam(10)))) -                                      &
      6*dRqsu*expuISs*qdB*qxy2**2*(1 + Rqsu)**2*tanhus**3*                      &
       (-(expuISr*(1 + Rqru)*(144*lam(2)**2 +                                   &
              576*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) +                  &
              192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +                  &
              216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                  &
              9*(lam(4) + lam(5) + lam(6) + lam(7))**2 +                        &
              256*(3*lam(1) + lam(3))*lam(10))) +                               &
         expdISr*(1 + Rqrd)*(144*lam(2)**2 +                                    &
            216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +                &
            9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                           &
               2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                       &
               10*lam(6)*lam(7) + 13*lam(7)**2 +                                &
               lam(5)*(-10*lam(6) + 26*lam(7))) +                               &
            288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -                &
            192*lam(8)**2 - 128*lam(3)*                                         &
             (3*lam(4) + 6*lam(5) - 3*(lam(6) + 2*lam(7) + 4*lam(9)) +          &
               7*lam(10)) + 48*lam(1)*                                          &
             (3*lam(4) + 33*lam(5) - 3*(lam(6) + 11*lam(7) + 16*lam(9)) +       &
               16*lam(10)))) +                                                  &
      6*dRqsd*expdISs*quB*qxy2**2*(1 + Rqsd)**2*tanhds**3*                      &
       (expdISr*(1 + Rqrd)*(144*lam(2)**2 +                                     &
            576*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) +                    &
            192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +                    &
            216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                    &
            9*(lam(4) + lam(5) + lam(6) + lam(7))**2 +                          &
            256*(3*lam(1) + lam(3))*lam(10)) -                                  &
         expuISr*(1 + Rqru)*(144*lam(2)**2 +                                    &
            216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +                &
            9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                           &
               2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                       &
               10*lam(6)*lam(7) + 13*lam(7)**2 +                                &
               lam(5)*(-10*lam(6) + 26*lam(7))) +                               &
            288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -                &
            192*lam(8)**2 - 128*lam(3)*                                         &
             (3*lam(4) + 6*lam(5) - 3*(lam(6) + 2*lam(7) + 4*lam(9)) +          &
               7*lam(10)) + 48*lam(1)*                                          &
             (3*lam(4) + 33*lam(5) - 3*(lam(6) + 11*lam(7) + 16*lam(9)) +       &
               16*lam(10)))) +                                                  &
      6*dRqsu*expuISs*qdB*qxy2*(1 + Rqsu)*tanhus*                               &
       (-(expuISr*(qxy2 + qzt2)*(1 + Rqru)*(1 + Rqsu)*                          &
            (144*lam(2)**2 + 576*lam(1)*                                        &
               (lam(4) - lam(5) - lam(6) + lam(7)) +                            &
              192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +                  &
              216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                  &
              9*(lam(4) + lam(5) + lam(6) + lam(7))**2 +                        &
              256*(3*lam(1) + lam(3))*lam(10))) -                               &
         4*expuISr*m2*(216*lam(2)**2 + 64*(3*lam(1) + lam(3))**2 +              &
            36*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                     &
            27*((lam(5) + lam(6))**2 + (lam(4) + lam(7))**2) +                  &
            48*lam(10)**2) + expdISr*(qxy2 + qzt2)*(1 + Rqrd)*(1 + Rqsu)*       &
          (144*lam(2)**2 + 216*lam(2)*                                          &
             (lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +                          &
            9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                           &
               2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                       &
               10*lam(6)*lam(7) + 13*lam(7)**2 +                                &
               lam(5)*(-10*lam(6) + 26*lam(7))) +                               &
            288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -                &
            192*lam(8)**2 - 128*lam(3)*                                         &
             (3*lam(4) + 6*lam(5) - 3*(lam(6) + 2*lam(7) + 4*lam(9)) +          &
               7*lam(10)) + 48*lam(1)*                                          &
             (3*lam(4) + 33*lam(5) -                                            &
               3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))) +              &
         4*expdISr*m2*(144*lam(1)**2 - 768*lam(1)*lam(3) +                      &
            160*lam(3)**2 + 9*                                                  &
             (24*lam(2)**2 + 3*                                                 &
                (lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                  &
                  lam(6)**2 - 10*lam(4)*lam(7) + 13*lam(7)**2) +                &
               32*lam(8)**2 +                                                   &
               4*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7) +                &
                  16*lam(8))) +                                                 &
            48*(6*lam(9)**2 - 6*lam(9)*lam(10) + lam(10)**2))) -                &
      6*dRqsd*expdISs*quB*qxy2*(1 + Rqsd)*tanhds*                               &
       (expdISr*(qxy2 + qzt2)*(1 + Rqrd)*(1 + Rqsd)*                            &
          (144*lam(2)**2 + 576*lam(1)*                                          &
             (lam(4) - lam(5) - lam(6) + lam(7)) +                              &
            192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +                    &
            216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                    &
            9*(lam(4) + lam(5) + lam(6) + lam(7))**2 +                          &
            256*(3*lam(1) + lam(3))*lam(10)) +                                  &
         4*expdISr*m2*(216*lam(2)**2 + 64*(3*lam(1) + lam(3))**2 +              &
            36*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                     &
            27*((lam(5) + lam(6))**2 + (lam(4) + lam(7))**2) +                  &
            48*lam(10)**2) - expuISr*(qxy2 + qzt2)*(1 + Rqru)*(1 + Rqsd)*       &
          (144*lam(2)**2 + 216*lam(2)*                                          &
             (lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +                          &
            9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                           &
               2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                       &
               10*lam(6)*lam(7) + 13*lam(7)**2 +                                &
               lam(5)*(-10*lam(6) + 26*lam(7))) +                               &
            288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -                &
            192*lam(8)**2 - 128*lam(3)*                                         &
             (3*lam(4) + 6*lam(5) - 3*(lam(6) + 2*lam(7) + 4*lam(9)) +          &
               7*lam(10)) + 48*lam(1)*                                          &
             (3*lam(4) + 33*lam(5) -                                            &
               3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))) -              &
         4*expuISr*m2*(144*lam(1)**2 - 768*lam(1)*lam(3) +                      &
            160*lam(3)**2 + 9*                                                  &
             (24*lam(2)**2 + 3*                                                 &
                (lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                  &
                  lam(6)**2 - 10*lam(4)*lam(7) + 13*lam(7)**2) +                &
               32*lam(8)**2 +                                                   &
               4*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7) +                &
                  16*lam(8))) +                                                 &
            48*(6*lam(9)**2 - 6*lam(9)*lam(10) + lam(10)**2))) +                &
      3*dRqsd*expdISs*expuISr*qdB*quB*                                          &
       ((qxy2 + qzt2)*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                &
          (144*lam(2)**2 + 216*lam(2)*                                          &
             (lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +                          &
            9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                           &
               2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                       &
               10*lam(6)*lam(7) + 13*lam(7)**2 +                                &
               lam(5)*(-10*lam(6) + 26*lam(7))) +                               &
            288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -                &
            192*lam(8)**2 - 128*lam(3)*                                         &
             (3*lam(4) + 6*lam(5) - 3*(lam(6) + 2*lam(7) + 4*lam(9)) +          &
               7*lam(10)) + 48*lam(1)*                                          &
             (3*lam(4) + 33*lam(5) -                                            &
               3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))) +              &
         8*m2*qzt2*(1 + Rqsd)*ss*                                               &
          (144*lam(1)**2 - 768*lam(1)*lam(3) + 160*lam(3)**2 +                  &
            9*(24*lam(2)**2 +                                                   &
               3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                 &
                  lam(6)**2 - 10*lam(4)*lam(7) + 13*lam(7)**2) +                &
               32*lam(8)**2 +                                                   &
               4*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7) +                &
                  16*lam(8))) +                                                 &
            48*(6*lam(9)**2 - 6*lam(9)*lam(10) + lam(10)**2))) +                &
      3*dRqsu*expdISr*expuISs*qdB*quB*                                          &
       ((qxy2 + qzt2)*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                &
          (144*lam(2)**2 + 216*lam(2)*                                          &
             (lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +                          &
            9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                           &
               2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                       &
               10*lam(6)*lam(7) + 13*lam(7)**2 +                                &
               lam(5)*(-10*lam(6) + 26*lam(7))) +                               &
            288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -                &
            192*lam(8)**2 - 128*lam(3)*                                         &
             (3*lam(4) + 6*lam(5) - 3*(lam(6) + 2*lam(7) + 4*lam(9)) +          &
               7*lam(10)) + 48*lam(1)*                                          &
             (3*lam(4) + 33*lam(5) -                                            &
               3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))) +              &
         8*m2*qzt2*(1 + Rqsu)*ss*                                               &
          (144*lam(1)**2 - 768*lam(1)*lam(3) + 160*lam(3)**2 +                  &
            9*(24*lam(2)**2 +                                                   &
               3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                 &
                  lam(6)**2 - 10*lam(4)*lam(7) + 13*lam(7)**2) +                &
               32*lam(8)**2 +                                                   &
               4*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7) +                &
                  16*lam(8))) +                                                 &
            48*(6*lam(9)**2 - 6*lam(9)*lam(10) + lam(10)**2))) +                &
      3*expuISs*qdB*quB*qxy2*tanhus**2*                                         &
       (dRqrd*expdISr*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                &
          (576*lam(2)**2 + 9*lam(4)**2 + 18*lam(4)*lam(5) +                     &
            9*lam(5)**2 - 18*lam(4)*lam(6) - 18*lam(5)*lam(6) +                 &
            9*lam(6)**2 - 882*lam(4)*lam(7) - 882*lam(5)*lam(7) +               &
            882*lam(6)*lam(7) + 11241*lam(7)**2 - 144*lam(8)**2 +               &
            72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +                &
               2*lam(8)) + 36*lam(3)*                                           &
             (lam(4) - lam(5) + lam(6) - lam(7) + 10*lam(8)) +                  &
            132*lam(4)*lam(9) + 132*lam(5)*lam(9) - 132*lam(6)*lam(9) -         &
            4740*lam(7)*lam(9) + 412*lam(9)**2 +                                &
            96*lam(2)*(3*lam(4) + 3*lam(5) -                                    &
               3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) +                 &
            24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) + 10*lam(9))*         &
             lam(10) - 144*lam(10)**2) -                                        &
         dRqru*expuISr*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*               &
          (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +                     &
            117*lam(5)**2 - 18*lam(4)*lam(6) + 90*lam(5)*lam(6) +               &
            9*lam(6)**2 + 522*lam(4)*lam(7) + 198*lam(5)*lam(7) -               &
            522*lam(6)*lam(7) - 5499*lam(7)**2 +                                &
            36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +           &
            288*lam(8)**2 + 72*lam(1)*                                          &
             (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) + 2*lam(8)) -           &
            156*lam(4)*lam(9) + 204*lam(5)*lam(9) + 156*lam(6)*lam(9) +         &
            2100*lam(7)*lam(9) - 56*lam(9)**2 -                                 &
            48*lam(2)*(3*(lam(4) + lam(5) - lam(6) - 25*lam(7) +                &
                  4*lam(8)) + 16*lam(9)) -                                      &
            36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -                   &
               10*lam(9))*lam(10) + 456*lam(10)**2) -                           &
         dRqsu*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                                  &
          (-(expuISr*(1 + Rqru)*                                                &
               (144*lam(2)**2 +                                                 &
                 576*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) +               &
                 192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +               &
                 216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +               &
                 9*(lam(4) + lam(5) + lam(6) + lam(7))**2 +                     &
                 256*(3*lam(1) + lam(3))*lam(10))) +                            &
            expdISr*(1 + Rqrd)*                                                 &
             (144*lam(2)**2 +                                                   &
               216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +             &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                    &
                  10*lam(6)*lam(7) + 13*lam(7)**2 +                             &
                  lam(5)*(-10*lam(6) + 26*lam(7))) +                            &
               288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -             &
               192*lam(8)**2 -                                                  &
               128*lam(3)*(3*lam(4) + 6*lam(5) -                                &
                  3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) +               &
               48*lam(1)*(3*lam(4) + 33*lam(5) -                                &
                  3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))))) -         &
      3*expdISs*qdB*quB*qxy2*tanhds**2*                                         &
       (-(dRqru*expuISr*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*              &
            (576*lam(2)**2 + 9*lam(4)**2 + 18*lam(4)*lam(5) +                   &
              9*lam(5)**2 - 18*lam(4)*lam(6) - 18*lam(5)*lam(6) +               &
              9*lam(6)**2 - 882*lam(4)*lam(7) - 882*lam(5)*lam(7) +             &
              882*lam(6)*lam(7) + 11241*lam(7)**2 - 144*lam(8)**2 +             &
              72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +              &
                 2*lam(8)) + 36*lam(3)*                                         &
               (lam(4) - lam(5) + lam(6) - lam(7) + 10*lam(8)) +                &
              132*lam(4)*lam(9) + 132*lam(5)*lam(9) -                           &
              132*lam(6)*lam(9) - 4740*lam(7)*lam(9) + 412*lam(9)**2 +          &
              96*lam(2)*(3*lam(4) + 3*lam(5) -                                  &
                 3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) +               &
              24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                   &
                 10*lam(9))*lam(10) - 144*lam(10)**2)) +                        &
         dRqrd*expdISr*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*               &
          (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +                     &
            117*lam(5)**2 - 18*lam(4)*lam(6) + 90*lam(5)*lam(6) +               &
            9*lam(6)**2 + 522*lam(4)*lam(7) + 198*lam(5)*lam(7) -               &
            522*lam(6)*lam(7) - 5499*lam(7)**2 +                                &
            36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +           &
            288*lam(8)**2 + 72*lam(1)*                                          &
             (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) + 2*lam(8)) -           &
            156*lam(4)*lam(9) + 204*lam(5)*lam(9) + 156*lam(6)*lam(9) +         &
            2100*lam(7)*lam(9) - 56*lam(9)**2 -                                 &
            48*lam(2)*(3*(lam(4) + lam(5) - lam(6) - 25*lam(7) +                &
                  4*lam(8)) + 16*lam(9)) -                                      &
            36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -                   &
               10*lam(9))*lam(10) + 456*lam(10)**2) -                           &
         dRqsd*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                                  &
          (expdISr*(1 + Rqrd)*                                                  &
             (144*lam(2)**2 +                                                   &
               576*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) +                 &
               192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +                 &
               216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                 &
               9*(lam(4) + lam(5) + lam(6) + lam(7))**2 +                       &
               256*(3*lam(1) + lam(3))*lam(10)) -                               &
            expuISr*(1 + Rqru)*                                                 &
             (144*lam(2)**2 +                                                   &
               216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +             &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                    &
                  10*lam(6)*lam(7) + 13*lam(7)**2 +                             &
                  lam(5)*(-10*lam(6) + 26*lam(7))) +                            &
               288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -             &
               192*lam(8)**2 -                                                  &
               128*lam(3)*(3*lam(4) + 6*lam(5) -                                &
                  3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) +               &
               48*lam(1)*(3*lam(4) + 33*lam(5) -                                &
                  3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))))) -         &
      expdISr*tanhdr*(-6*dRqrd*quB*qxy2*(1 + Rqrd)*                             &
          (-(expuISs*(qxy2 + qzt2)*(1 + Rqrd)*(1 + Rqsu)*                       &
               (576*lam(2)**2 + 9*lam(4)**2 + 18*lam(4)*lam(5) +                &
                 9*lam(5)**2 - 18*lam(4)*lam(6) - 18*lam(5)*lam(6) +            &
                 9*lam(6)**2 - 882*lam(4)*lam(7) - 882*lam(5)*lam(7) +          &
                 882*lam(6)*lam(7) + 11241*lam(7)**2 - 144*lam(8)**2 +          &
                 72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +           &
                    2*lam(8)) +                                                 &
                 36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                 &
                    10*lam(8)) + 132*lam(4)*lam(9) +                            &
                 132*lam(5)*lam(9) - 132*lam(6)*lam(9) -                        &
                 4740*lam(7)*lam(9) + 412*lam(9)**2 +                           &
                 96*lam(2)*(3*lam(4) + 3*lam(5) -                               &
                    3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) +            &
                 24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                &
                    10*lam(9))*lam(10) - 144*lam(10)**2)) -                     &
            12*expdISs*m2*(72*lam(1)**2 - 48*lam(2)**2 +                        &
               168*lam(1)*lam(3) + 60*lam(3)**2 + 9*lam(4)**2 +                 &
               117*lam(5)**2 + 90*lam(5)*lam(6) + 9*lam(6)**2 +                 &
               126*lam(4)*lam(7) + 36*lam(5)*lam(7) -                           &
               36*lam(6)*lam(7) - 351*lam(7)**2 + 60*lam(5)*lam(8) +            &
               12*lam(6)*lam(8) - 12*lam(8)**2 + 192*lam(7)*lam(9) +            &
               40*lam(8)*lam(9) + 16*lam(9)**2 +                                &
               4*lam(2)*(3*lam(4) - 15*lam(5) - 3*lam(6) + 87*lam(7) -          &
                  26*lam(9) - 30*lam(10)) +                                     &
               12*(3*lam(8) - 8*lam(9))*lam(10) + 96*lam(10)**2) +              &
            expdISs*(qxy2 + qzt2)*(1 + Rqrd)*(1 + Rqsd)*                        &
             (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +                  &
               117*lam(5)**2 - 18*lam(4)*lam(6) + 90*lam(5)*lam(6) +            &
               9*lam(6)**2 + 522*lam(4)*lam(7) + 198*lam(5)*lam(7) -            &
               522*lam(6)*lam(7) - 5499*lam(7)**2 +                             &
               36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) +                     &
                  23*lam(7)) + 288*lam(8)**2 +                                  &
               72*lam(1)*(3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +           &
                  2*lam(8)) - 156*lam(4)*lam(9) + 204*lam(5)*lam(9) +           &
               156*lam(6)*lam(9) + 2100*lam(7)*lam(9) - 56*lam(9)**2 -          &
               48*lam(2)*(3*(lam(4) + lam(5) - lam(6) - 25*lam(7) +             &
                     4*lam(8)) + 16*lam(9)) -                                   &
               36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -                &
                  10*lam(9))*lam(10) + 456*lam(10)**2) +                        &
            12*expuISs*m2*(72*lam(1)**2 + 96*lam(2)**2 +                        &
               24*lam(1)*lam(3) + 30*lam(3)**2 + 9*lam(4)**2 +                  &
               9*(lam(5) - lam(6))**2 - 90*lam(4)*lam(7) -                      &
               72*lam(5)*lam(7) + 72*lam(6)*lam(7) + 945*lam(7)**2 -            &
               12*lam(5)*lam(8) + 12*lam(6)*lam(8) + 30*lam(8)**2 -             &
               384*lam(7)*lam(9) - 20*lam(8)*lam(9) + 16*lam(9)**2 +            &
               24*lam(8)*lam(10) +                                              &
               4*lam(2)*(3*lam(4) + 3*lam(5) - 3*lam(6) - 147*lam(7) +          &
                  22*lam(9) + 12*lam(10)))) -                                   &
         2*expdISs*quB*qxy2*(1 + Rqsd)*tanhds**2*                               &
          (-3*dRqrd*qxy2*(1 + Rqrd)**2*                                         &
             (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +                  &
               117*lam(5)**2 - 18*lam(4)*lam(6) + 90*lam(5)*lam(6) +            &
               9*lam(6)**2 + 522*lam(4)*lam(7) + 198*lam(5)*lam(7) -            &
               522*lam(6)*lam(7) - 5499*lam(7)**2 +                             &
               36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) +                     &
                  23*lam(7)) + 288*lam(8)**2 +                                  &
               72*lam(1)*(3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +           &
                  2*lam(8)) - 156*lam(4)*lam(9) + 204*lam(5)*lam(9) +           &
               156*lam(6)*lam(9) + 2100*lam(7)*lam(9) - 56*lam(9)**2 -          &
               48*lam(2)*(3*(lam(4) + lam(5) - lam(6) - 25*lam(7) +             &
                     4*lam(8)) + 16*lam(9)) -                                   &
               36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -                &
                  10*lam(9))*lam(10) + 456*lam(10)**2) +                        &
            dRqsd*(6*m2*(128*(3*lam(1) + lam(3))**2 +                           &
                  15*(lam(4) - lam(5) - lam(6) + lam(7))**2) +                  &
               384*m2*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +             &
               64*m2*lam(10)**2 -                                               &
               3*qzt2*(1 + Rqrd)*(1 + Rqsd)*                                    &
                (144*lam(2)**2 -                                                &
                  576*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) -              &
                  192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +              &
                  216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +              &
                  9*(lam(4) + lam(5) + lam(6) + lam(7))**2 -                    &
                  256*(3*lam(1) + lam(3))*lam(10)))) +                          &
         expuISs*qdB*quB*qzt2*tanhus*                                           &
          (-1728*dRqrd*lam(2)**2 +                                              &
            3*dRqrd*Rqsu*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                        &
             (576*lam(2)**2 + 9*lam(4)**2 + 18*lam(4)*lam(5) +                  &
               9*lam(5)**2 - 18*lam(4)*lam(6) - 18*lam(5)*lam(6) +              &
               9*lam(6)**2 - 882*lam(4)*lam(7) - 882*lam(5)*lam(7) +            &
               882*lam(6)*lam(7) + 11241*lam(7)**2 - 144*lam(8)**2 -            &
               72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +             &
                  2*lam(8)) -                                                   &
               36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                   &
                  10*lam(8)) + 132*lam(4)*lam(9) + 132*lam(5)*lam(9) -          &
               132*lam(6)*lam(9) - 4740*lam(7)*lam(9) + 412*lam(9)**2 +         &
               96*lam(2)*(3*lam(4) + 3*lam(5) -                                 &
                  3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) +              &
               24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                  &
                  10*lam(9))*lam(10) - 144*lam(10)**2) +                        &
            4*dRqsu*m2*(1 + Rqsu)*ss*                                           &
             (96*(9*lam(1)**2 - 48*lam(1)*lam(3) + 10*lam(3)**2) -              &
               9*(lam(4)**2 - 11*lam(5)**2 + 14*lam(5)*lam(6) +                 &
                  lam(6)**2 - 2*lam(4)*(7*lam(5) + lam(6))) -                   &
               18*(7*lam(4) + 11*lam(5) - 7*lam(6))*lam(7) +                    &
               99*lam(7)**2 -                                                   &
               576*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +             &
               960*lam(9)**2 +                                                  &
               192*(lam(4) + 5*lam(5) - lam(6) - 5*(lam(7) + lam(9)))*          &
                lam(10) + 224*lam(10)**2) -                                     &
            4*dRqrd*m2*(1 + Rqrd)*sr*                                           &
             (1728*lam(2)**2 + 9*lam(4)**2 + 9*lam(5)**2 -                      &
               18*lam(5)*lam(6) + 9*lam(6)**2 - 1314*lam(5)*lam(7) +            &
               1314*lam(6)*lam(7) + 16857*lam(7)**2 -                           &
               108*lam(5)*lam(8) + 108*lam(6)*lam(8) +                          &
               108*lam(7)*lam(8) + 384*lam(5)*lam(9) -                          &
               384*lam(6)*lam(9) - 7296*lam(7)*lam(9) -                         &
               360*lam(8)*lam(9) + 928*lam(9)**2 +                              &
               6*lam(4)*(3*lam(5) -                                             &
                  3*(lam(6) + 73*lam(7) + 6*lam(8)) + 64*lam(9)) +              &
               432*lam(8)*lam(10) +                                             &
               72*lam(2)*(3*lam(4) + 3*lam(5) - 3*lam(6) - 147*lam(7) +         &
                  22*lam(9) + 12*lam(10))) -                                    &
            3*dRqsu*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                  &
             (144*lam(2)**2 +                                                   &
               216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +             &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                    &
                  10*lam(6)*lam(7) + 13*lam(7)**2 +                             &
                  lam(5)*(-10*lam(6) + 26*lam(7))) +                            &
               288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -             &
               192*lam(8)**2 +                                                  &
               128*lam(3)*(3*lam(4) + 6*lam(5) -                                &
                  3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) -               &
               48*lam(1)*(3*lam(4) + 33*lam(5) -                                &
                  3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))) +           &
            3*dRqrd*(216*lam(1)*lam(4) + 36*lam(3)*lam(4) -                     &
               9*lam(4)**2 - 36*lam(3)*lam(5) - 18*lam(4)*lam(5) -              &
               9*lam(5)**2 + 36*lam(3)*lam(6) + 18*lam(4)*lam(6) +              &
               18*lam(5)*lam(6) - 9*lam(6)**2 - 36*lam(3)*lam(7) +              &
               882*lam(4)*lam(7) + 882*lam(5)*lam(7) -                          &
               882*lam(6)*lam(7) - 11241*lam(7)**2 + 360*lam(3)*lam(8) +        &
               144*lam(8)**2 +                                                  &
               72*lam(1)*(-3*(lam(5) - lam(6) + lam(7)) + 2*lam(8)) -           &
               132*lam(4)*lam(9) - 132*lam(5)*lam(9) +                          &
               132*lam(6)*lam(9) + 4740*lam(7)*lam(9) - 412*lam(9)**2 -         &
               96*lam(2)*(3*lam(4) + 3*lam(5) -                                 &
                  3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) -              &
               24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                  &
                  10*lam(9))*lam(10) + 144*lam(10)**2 +                         &
               2*qzt2*(1 + Rqrd)**2*sr*                                         &
                (576*lam(2)**2 + 9*lam(4)**2 + 18*lam(4)*lam(5) +               &
                  9*lam(5)**2 - 18*lam(4)*lam(6) - 18*lam(5)*lam(6) +           &
                  9*lam(6)**2 - 882*lam(4)*lam(7) - 882*lam(5)*lam(7) +         &
                  882*lam(6)*lam(7) + 11241*lam(7)**2 - 144*lam(8)**2 -         &
                  72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +          &
                     2*lam(8)) -                                                &
                  36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                &
                     10*lam(8)) + 132*lam(4)*lam(9) +                           &
                  132*lam(5)*lam(9) - 132*lam(6)*lam(9) -                       &
                  4740*lam(7)*lam(9) + 412*lam(9)**2 +                          &
                  96*lam(2)*(3*lam(4) + 3*lam(5) -                              &
                     3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) +           &
                  24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +               &
                     10*lam(9))*lam(10) - 144*lam(10)**2))) -                   &
         expdISs*qdB*quB*qzt2*tanhds*                                           &
          (-1728*dRqrd*lam(2)**2 -                                              &
            3*dRqsd*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                  &
             (3*(48*lam(2)**2 -                                                 &
                  192*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) -              &
                  64*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +               &
                  72*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +               &
                  3*(lam(4) + lam(5) + lam(6) + lam(7))**2) -                   &
               256*(3*lam(1) + lam(3))*lam(10)) +                               &
            4*dRqsd*m2*(1 + Rqsd)*ss*                                           &
             (384*(3*lam(1) + lam(3))**2 +                                      &
               45*(lam(4) - lam(5) - lam(6) + lam(7))**2 +                      &
               192*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +                &
               32*lam(10)**2) +                                                 &
            3*dRqrd*Rqsd*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                        &
             (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +                  &
               117*lam(5)**2 - 18*lam(4)*lam(6) + 90*lam(5)*lam(6) +            &
               9*lam(6)**2 + 522*lam(4)*lam(7) + 198*lam(5)*lam(7) -            &
               522*lam(6)*lam(7) - 5499*lam(7)**2 -                             &
               36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) +                     &
                  23*lam(7)) + 288*lam(8)**2 -                                  &
               72*lam(1)*(3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +           &
                  2*lam(8)) - 156*lam(4)*lam(9) + 204*lam(5)*lam(9) +           &
               156*lam(6)*lam(9) + 2100*lam(7)*lam(9) - 56*lam(9)**2 -          &
               48*lam(2)*(3*(lam(4) + lam(5) - lam(6) - 25*lam(7) +             &
                     4*lam(8)) + 16*lam(9)) -                                   &
               36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -                &
                  10*lam(9))*lam(10) + 456*lam(10)**2) +                        &
            4*dRqrd*m2*(1 + Rqrd)*sr*                                           &
             (864*lam(2)**2 + 45*lam(4)**2 + 369*lam(5)**2 +                    &
               234*lam(5)*lam(6) + 45*lam(6)**2 - 1386*lam(5)*lam(7) +          &
               414*lam(6)*lam(7) + 8793*lam(7)**2 - 540*lam(5)*lam(8) -         &
               108*lam(6)*lam(8) + 540*lam(7)*lam(8) -                          &
               6*lam(4)*(39*lam(5) + 15*lam(6) + 69*lam(7) -                    &
                  18*lam(8) - 32*lam(9)) + 192*lam(5)*lam(9) -                  &
               192*lam(6)*lam(9) - 3648*lam(7)*lam(9) -                         &
               720*lam(8)*lam(9) + 608*lam(9)**2 -                              &
               72*lam(2)*(3*lam(4) - 15*lam(5) - 3*lam(6) +                     &
                  87*lam(7) - 26*lam(9) - 30*lam(10)) -                         &
               72*(9*lam(8) + 8*lam(9))*lam(10) + 576*lam(10)**2) +             &
            3*dRqrd*(216*lam(1)*lam(4) + 252*lam(3)*lam(4) -                    &
               9*lam(4)**2 + 828*lam(3)*lam(5) + 90*lam(4)*lam(5) -             &
               117*lam(5)**2 + 252*lam(3)*lam(6) + 18*lam(4)*lam(6) -           &
               90*lam(5)*lam(6) - 9*lam(6)**2 + 828*lam(3)*lam(7) -             &
               522*lam(4)*lam(7) - 198*lam(5)*lam(7) +                          &
               522*lam(6)*lam(7) + 5499*lam(7)**2 +                             &
               72*(-4*lam(8)**2 +                                               &
                  lam(1)*(3*(5*lam(5) + lam(6) + 5*lam(7)) + 2*lam(8)))         &
+ 156*lam(4)*lam(9) - 204*lam(5)*lam(9) - 156*lam(6)*lam(9) -                   &
               2100*lam(7)*lam(9) + 56*lam(9)**2 +                              &
               48*lam(2)*(3*(lam(4) + lam(5) - lam(6) - 25*lam(7) +             &
                     4*lam(8)) + 16*lam(9)) +                                   &
               36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -                &
                  10*lam(9))*lam(10) - 456*lam(10)**2 +                         &
               2*qzt2*(1 + Rqrd)**2*sr*                                         &
                (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +               &
                  117*lam(5)**2 - 18*lam(4)*lam(6) + 90*lam(5)*lam(6) +         &
                  9*lam(6)**2 + 522*lam(4)*lam(7) + 198*lam(5)*lam(7) -         &
                  522*lam(6)*lam(7) - 5499*lam(7)**2 -                          &
                  36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) +                  &
                     23*lam(7)) + 288*lam(8)**2 -                               &
                  72*lam(1)*(3*                                                 &
                      (lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) + 2*lam(8)        &
) - 156*lam(4)*lam(9) + 204*lam(5)*lam(9) + 156*lam(6)*lam(9) +                 &
                  2100*lam(7)*lam(9) - 56*lam(9)**2 -                           &
                  48*lam(2)*(3*                                                 &
                      (lam(4) + lam(5) - lam(6) - 25*lam(7) + 4*lam(8))         &
+ 16*lam(9)) - 36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -                &
                     10*lam(9))*lam(10) + 456*lam(10)**2))) -                   &
         2*expuISs*qxy2*(1 + Rqsu)*tanhus**2*                                   &
          (3*dRqrd*quB*qxy2*(1 + Rqrd)**2*                                      &
             (576*lam(2)**2 + 9*lam(4)**2 + 18*lam(4)*lam(5) +                  &
               9*lam(5)**2 - 18*lam(4)*lam(6) - 18*lam(5)*lam(6) +              &
               9*lam(6)**2 - 882*lam(4)*lam(7) - 882*lam(5)*lam(7) +            &
               882*lam(6)*lam(7) + 11241*lam(7)**2 - 144*lam(8)**2 +            &
               72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +             &
                  2*lam(8)) +                                                   &
               36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                   &
                  10*lam(8)) + 132*lam(4)*lam(9) + 132*lam(5)*lam(9) -          &
               132*lam(6)*lam(9) - 4740*lam(7)*lam(9) + 412*lam(9)**2 +         &
               96*lam(2)*(3*lam(4) + 3*lam(5) -                                 &
                  3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) +              &
               24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                  &
                  10*lam(9))*lam(10) - 144*lam(10)**2) +                        &
            dRqsu*qdB*(6*m2*(-32*                                               &
                   (9*lam(1)**2 - 48*lam(1)*lam(3) + 10*lam(3)**2) +            &
                  3*(lam(4)**2 - 11*lam(5)**2 + lam(6)**2 -                     &
                     2*lam(4)*(7*lam(5) + lam(6) - 7*lam(7)) -                  &
                     14*lam(6)*lam(7) - 11*lam(7)**2 +                          &
                     2*lam(5)*(7*lam(6) + 11*lam(7))) +                         &
                  192*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) -          &
                  320*lam(9)**2) -                                              &
               384*m2*(lam(4) + 5*lam(5) - lam(6) - 5*(lam(7) + lam(9)))*       &
                lam(10) - 448*m2*lam(10)**2 +                                   &
               3*qzt2*(1 + Rqrd)*(1 + Rqsu)*                                    &
                (144*lam(2)**2 +                                                &
                  216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +          &
                  9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                     &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                 &
                     10*lam(6)*lam(7) + 13*lam(7)**2 +                          &
                     lam(5)*(-10*lam(6) + 26*lam(7))) +                         &
                  288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -          &
                  192*lam(8)**2 +                                               &
                  128*lam(3)*(3*lam(4) + 6*lam(5) -                             &
                     3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) -            &
                  48*lam(1)*(3*lam(4) + 33*lam(5) -                             &
                     3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))))))       &
- expdISr*qxy2*tanhdr**2*(6*dRqsd*expdISs*quB*qxy2*(1 + Rqrd)*                  &
          (1 + Rqsd)**2*tanhds**3*                                              &
          (144*lam(2)**2 + 576*lam(1)*                                          &
             (lam(4) - lam(5) - lam(6) + lam(7)) +                              &
            192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +                    &
            216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                    &
            9*(lam(4) + lam(5) + lam(6) + lam(7))**2 +                          &
            256*(3*lam(1) + lam(3))*lam(10)) -                                  &
         6*dRqsu*expuISs*qdB*qxy2*(1 + Rqrd)*(1 + Rqsu)**2*tanhus**3*           &
          (144*lam(2)**2 + 216*lam(2)*                                          &
             (lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +                          &
            9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                           &
               2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                       &
               10*lam(6)*lam(7) + 13*lam(7)**2 +                                &
               lam(5)*(-10*lam(6) + 26*lam(7))) +                               &
            288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -                &
            192*lam(8)**2 - 128*lam(3)*                                         &
             (3*lam(4) + 6*lam(5) - 3*(lam(6) + 2*lam(7) + 4*lam(9)) +          &
               7*lam(10)) + 48*lam(1)*                                          &
             (3*lam(4) + 33*lam(5) -                                            &
               3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))) -              &
         3*expdISs*qdB*quB*tanhds**2*                                           &
          (-(dRqsd*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                   &
               (144*lam(2)**2 +                                                 &
                 576*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) +               &
                 192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +               &
                 216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +               &
                 9*(lam(4) + lam(5) + lam(6) + lam(7))**2 +                     &
                 256*(3*lam(1) + lam(3))*lam(10))) +                            &
            dRqrd*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                    &
             (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +                  &
               117*lam(5)**2 - 18*lam(4)*lam(6) + 90*lam(5)*lam(6) +            &
               9*lam(6)**2 + 522*lam(4)*lam(7) + 198*lam(5)*lam(7) -            &
               522*lam(6)*lam(7) - 5499*lam(7)**2 +                             &
               36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +        &
               288*lam(8)**2 +                                                  &
               72*lam(1)*(3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +           &
                  2*lam(8)) - 156*lam(4)*lam(9) + 204*lam(5)*lam(9) +           &
               156*lam(6)*lam(9) + 2100*lam(7)*lam(9) - 56*lam(9)**2 -          &
               48*lam(2)*(3*(lam(4) + lam(5) - lam(6) - 25*lam(7) +             &
                     4*lam(8)) + 16*lam(9)) -                                   &
               36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -                &
                  10*lam(9))*lam(10) + 456*lam(10)**2)) +                       &
         3*expuISs*qdB*quB*tanhus**2*                                           &
          (dRqrd*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                     &
             (576*lam(2)**2 + 9*lam(4)**2 + 18*lam(4)*lam(5) +                  &
               9*lam(5)**2 - 18*lam(4)*lam(6) - 18*lam(5)*lam(6) +              &
               9*lam(6)**2 - 882*lam(4)*lam(7) - 882*lam(5)*lam(7) +            &
               882*lam(6)*lam(7) + 11241*lam(7)**2 - 144*lam(8)**2 +            &
               72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +             &
                  2*lam(8)) +                                                   &
               36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                   &
                  10*lam(8)) + 132*lam(4)*lam(9) + 132*lam(5)*lam(9) -          &
               132*lam(6)*lam(9) - 4740*lam(7)*lam(9) + 412*lam(9)**2 +         &
               96*lam(2)*(3*lam(4) + 3*lam(5) -                                 &
                  3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) +              &
               24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                  &
                  10*lam(9))*lam(10) - 144*lam(10)**2) -                        &
            dRqsu*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                    &
             (144*lam(2)**2 +                                                   &
               216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +             &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                    &
                  10*lam(6)*lam(7) + 13*lam(7)**2 +                             &
                  lam(5)*(-10*lam(6) + 26*lam(7))) +                            &
               288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -             &
               192*lam(8)**2 -                                                  &
               128*lam(3)*(3*lam(4) + 6*lam(5) -                                &
                  3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) +               &
               48*lam(1)*(3*lam(4) + 33*lam(5) -                                &
                  3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10)))) +          &
         2*expuISs*(1 + Rqrd)*tanhus*                                           &
          (3*dRqsu*qdB*qxy2*(1 + Rqsu)**2*                                      &
             (144*lam(2)**2 +                                                   &
               216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +             &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                    &
                  10*lam(6)*lam(7) + 13*lam(7)**2 +                             &
                  lam(5)*(-10*lam(6) + 26*lam(7))) +                            &
               288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -             &
               192*lam(8)**2 -                                                  &
               128*lam(3)*(3*lam(4) + 6*lam(5) -                                &
                  3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) +               &
               48*lam(1)*(3*lam(4) + 33*lam(5) -                                &
                  3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))) +           &
            dRqrd*quB*(2*m2*(9*                                                 &
                   (-192*lam(2)**2 - (lam(4) + lam(5) - lam(6))**2 -            &
                     24*lam(2)*                                                 &
                      (lam(4) + lam(5) - lam(6) - 49*lam(7)) +                  &
                     146*(lam(4) + lam(5) - lam(6))*lam(7) -                    &
                     1873*lam(7)**2 +                                           &
                     12*(lam(4) + lam(5) - lam(6) - lam(7))*lam(8)) -           &
                  24*(66*lam(2) +                                               &
                     16*(lam(4) + lam(5) - lam(6) - 19*lam(7)) -                &
                     15*lam(8))*lam(9) - 928*lam(9)**2 -                        &
                  432*(2*lam(2) + lam(8))*lam(10)) +                            &
               3*qzt2*(1 + Rqrd)*(1 + Rqsu)*                                    &
                (576*lam(2)**2 + 9*lam(4)**2 + 18*lam(4)*lam(5) +               &
                  9*lam(5)**2 - 18*lam(4)*lam(6) - 18*lam(5)*lam(6) +           &
                  9*lam(6)**2 - 882*lam(4)*lam(7) - 882*lam(5)*lam(7) +         &
                  882*lam(6)*lam(7) + 11241*lam(7)**2 - 144*lam(8)**2 -         &
                  72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +          &
                     2*lam(8)) -                                                &
                  36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                &
                     10*lam(8)) + 132*lam(4)*lam(9) +                           &
                  132*lam(5)*lam(9) - 132*lam(6)*lam(9) -                       &
                  4740*lam(7)*lam(9) + 412*lam(9)**2 +                          &
                  96*lam(2)*(3*lam(4) + 3*lam(5) -                              &
                     3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) +           &
                  24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +               &
                     10*lam(9))*lam(10) - 144*lam(10)**2))) -                   &
         2*expdISs*quB*(1 + Rqrd)*tanhds*                                       &
          (3*dRqsd*qxy2*(1 + Rqsd)**2*                                          &
             (144*lam(2)**2 +                                                   &
               576*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) +                 &
               192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +                 &
               216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                 &
               9*(lam(4) + lam(5) + lam(6) + lam(7))**2 +                       &
               256*(3*lam(1) + lam(3))*lam(10)) +                               &
            dRqrd*(3*qzt2*(1 + Rqrd)*(1 + Rqsd)*                                &
                (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +               &
                  117*lam(5)**2 - 18*lam(4)*lam(6) + 90*lam(5)*lam(6) +         &
                  9*lam(6)**2 + 522*lam(4)*lam(7) + 198*lam(5)*lam(7) -         &
                  522*lam(6)*lam(7) - 5499*lam(7)**2 -                          &
                  36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) +                  &
                     23*lam(7)) + 288*lam(8)**2 -                               &
                  72*lam(1)*(3*                                                 &
                      (lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +                 &
                     2*lam(8)) - 156*lam(4)*lam(9) +                            &
                  204*lam(5)*lam(9) + 156*lam(6)*lam(9) +                       &
                  2100*lam(7)*lam(9) - 56*lam(9)**2 -                           &
                  48*lam(2)*(3*                                                 &
                      (lam(4) + lam(5) - lam(6) - 25*lam(7) +                   &
                        4*lam(8)) + 16*lam(9)) -                                &
                  36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -             &
                     10*lam(9))*lam(10) + 456*lam(10)**2) +                     &
               2*m2*(864*lam(2)**2 + 45*lam(4)**2 + 369*lam(5)**2 +             &
                  234*lam(5)*lam(6) + 45*lam(6)**2 -                            &
                  1386*lam(5)*lam(7) + 414*lam(6)*lam(7) +                      &
                  8793*lam(7)**2 - 540*lam(5)*lam(8) -                          &
                  108*lam(6)*lam(8) + 540*lam(7)*lam(8) -                       &
                  6*lam(4)*(39*lam(5) + 15*lam(6) + 69*lam(7) -                 &
                     18*lam(8) - 32*lam(9)) + 192*lam(5)*lam(9) -               &
                  192*lam(6)*lam(9) - 3648*lam(7)*lam(9) -                      &
                  720*lam(8)*lam(9) + 608*lam(9)**2 -                           &
                  72*lam(2)*(3*lam(4) - 15*lam(5) - 3*lam(6) +                  &
                     87*lam(7) - 26*lam(9) - 30*lam(10)) -                      &
                  72*(9*lam(8) + 8*lam(9))*lam(10) + 576*lam(10)**2))) +        &
         3*qdB*quB*(dRqrd*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                       &
             (-(expuISs*(1 + Rqsu)*                                             &
                  (576*lam(2)**2 + 9*lam(4)**2 + 18*lam(4)*lam(5) +             &
                    9*lam(5)**2 - 18*lam(4)*lam(6) - 18*lam(5)*lam(6) +         &
                    9*lam(6)**2 - 882*lam(4)*lam(7) -                           &
                    882*lam(5)*lam(7) + 882*lam(6)*lam(7) +                     &
                    11241*lam(7)**2 - 144*lam(8)**2 +                           &
                    72*lam(1)*                                                  &
                     (3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +                 &
                       2*lam(8)) +                                              &
                    36*lam(3)*                                                  &
                     (lam(4) - lam(5) + lam(6) - lam(7) + 10*lam(8)) +          &
                    132*lam(4)*lam(9) + 132*lam(5)*lam(9) -                     &
                    132*lam(6)*lam(9) - 4740*lam(7)*lam(9) +                    &
                    412*lam(9)**2 +                                             &
                    96*lam(2)*                                                  &
                     (3*lam(4) + 3*lam(5) -                                     &
                       3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) +         &
                    24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +             &
                       10*lam(9))*lam(10) - 144*lam(10)**2)) +                  &
               expdISs*(1 + Rqsd)*                                              &
                (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +               &
                  117*lam(5)**2 - 18*lam(4)*lam(6) + 90*lam(5)*lam(6) +         &
                  9*lam(6)**2 + 522*lam(4)*lam(7) + 198*lam(5)*lam(7) -         &
                  522*lam(6)*lam(7) - 5499*lam(7)**2 +                          &
                  36*lam(3)*(7*lam(4) + 23*lam(5) + 7*lam(6) +                  &
                     23*lam(7)) + 288*lam(8)**2 +                               &
                  72*lam(1)*(3*                                                 &
                      (lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) + 2*lam(8)        &
) - 156*lam(4)*lam(9) + 204*lam(5)*lam(9) + 156*lam(6)*lam(9) +                 &
                  2100*lam(7)*lam(9) - 56*lam(9)**2 -                           &
                  48*lam(2)*(3*                                                 &
                      (lam(4) + lam(5) - lam(6) - 25*lam(7) + 4*lam(8))         &
+ 16*lam(9)) - 36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -                &
                     10*lam(9))*lam(10) + 456*lam(10)**2)) -                    &
            (1 + Rqrd)*(dRqsd*expdISs*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*           &
                (144*lam(2)**2 +                                                &
                  576*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) +              &
                  192*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +              &
                  216*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +              &
                  9*(lam(4) + lam(5) + lam(6) + lam(7))**2 +                    &
                  256*(3*lam(1) + lam(3))*lam(10)) -                            &
               dRqsu*expuISs*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                    &
                (144*lam(2)**2 +                                                &
                  216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +          &
                  9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                     &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                 &
                     10*lam(6)*lam(7) + 13*lam(7)**2 +                          &
                     lam(5)*(-10*lam(6) + 26*lam(7))) +                         &
                  288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -          &
                  192*lam(8)**2 -                                               &
                  128*lam(3)*(3*lam(4) + 6*lam(5) -                             &
                     3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) +            &
                  48*lam(1)*(3*lam(4) + 33*lam(5) -                             &
                     3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))))))       &
- expuISr*tanhur*(expdISs*qdB*quB*qzt2*tanhds*                                  &
          (-1728*dRqru*lam(2)**2 +                                              &
            3*dRqru*Rqsd*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                        &
             (576*lam(2)**2 + 9*lam(4)**2 + 18*lam(4)*lam(5) +                  &
               9*lam(5)**2 - 18*lam(4)*lam(6) - 18*lam(5)*lam(6) +              &
               9*lam(6)**2 - 882*lam(4)*lam(7) - 882*lam(5)*lam(7) +            &
               882*lam(6)*lam(7) + 11241*lam(7)**2 - 144*lam(8)**2 -            &
               72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +             &
                  2*lam(8)) -                                                   &
               36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                   &
                  10*lam(8)) + 132*lam(4)*lam(9) + 132*lam(5)*lam(9) -          &
               132*lam(6)*lam(9) - 4740*lam(7)*lam(9) + 412*lam(9)**2 +         &
               96*lam(2)*(3*lam(4) + 3*lam(5) -                                 &
                  3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) +              &
               24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                  &
                  10*lam(9))*lam(10) - 144*lam(10)**2) +                        &
            4*dRqsd*m2*(1 + Rqsd)*ss*                                           &
             (96*(9*lam(1)**2 - 48*lam(1)*lam(3) + 10*lam(3)**2) -              &
               9*(lam(4)**2 - 11*lam(5)**2 + 14*lam(5)*lam(6) +                 &
                  lam(6)**2 - 2*lam(4)*(7*lam(5) + lam(6))) -                   &
               18*(7*lam(4) + 11*lam(5) - 7*lam(6))*lam(7) +                    &
               99*lam(7)**2 -                                                   &
               576*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) +             &
               960*lam(9)**2 +                                                  &
               192*(lam(4) + 5*lam(5) - lam(6) - 5*(lam(7) + lam(9)))*          &
                lam(10) + 224*lam(10)**2) -                                     &
            4*dRqru*m2*(1 + Rqru)*sr*                                           &
             (1728*lam(2)**2 + 9*lam(4)**2 + 9*lam(5)**2 -                      &
               18*lam(5)*lam(6) + 9*lam(6)**2 - 1314*lam(5)*lam(7) +            &
               1314*lam(6)*lam(7) + 16857*lam(7)**2 -                           &
               108*lam(5)*lam(8) + 108*lam(6)*lam(8) +                          &
               108*lam(7)*lam(8) + 384*lam(5)*lam(9) -                          &
               384*lam(6)*lam(9) - 7296*lam(7)*lam(9) -                         &
               360*lam(8)*lam(9) + 928*lam(9)**2 +                              &
               6*lam(4)*(3*lam(5) -                                             &
                  3*(lam(6) + 73*lam(7) + 6*lam(8)) + 64*lam(9)) +              &
               432*lam(8)*lam(10) +                                             &
               72*lam(2)*(3*lam(4) + 3*lam(5) - 3*lam(6) - 147*lam(7) +         &
                  22*lam(9) + 12*lam(10))) -                                    &
            3*dRqsd*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                  &
             (144*lam(2)**2 +                                                   &
               216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +             &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                    &
                  10*lam(6)*lam(7) + 13*lam(7)**2 +                             &
                  lam(5)*(-10*lam(6) + 26*lam(7))) +                            &
               288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -             &
               192*lam(8)**2 +                                                  &
               128*lam(3)*(3*lam(4) + 6*lam(5) -                                &
                  3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) -               &
               48*lam(1)*(3*lam(4) + 33*lam(5) -                                &
                  3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))) +           &
            3*dRqru*(216*lam(1)*lam(4) + 36*lam(3)*lam(4) -                     &
               9*lam(4)**2 - 36*lam(3)*lam(5) - 18*lam(4)*lam(5) -              &
               9*lam(5)**2 + 36*lam(3)*lam(6) + 18*lam(4)*lam(6) +              &
               18*lam(5)*lam(6) - 9*lam(6)**2 - 36*lam(3)*lam(7) +              &
               882*lam(4)*lam(7) + 882*lam(5)*lam(7) -                          &
               882*lam(6)*lam(7) - 11241*lam(7)**2 + 360*lam(3)*lam(8) +        &
               144*lam(8)**2 +                                                  &
               72*lam(1)*(-3*(lam(5) - lam(6) + lam(7)) + 2*lam(8)) -           &
               132*lam(4)*lam(9) - 132*lam(5)*lam(9) +                          &
               132*lam(6)*lam(9) + 4740*lam(7)*lam(9) - 412*lam(9)**2 -         &
               96*lam(2)*(3*lam(4) + 3*lam(5) -                                 &
                  3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) -              &
               24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                  &
                  10*lam(9))*lam(10) + 144*lam(10)**2 +                         &
               2*qzt2*(1 + Rqru)**2*sr*                                         &
                (576*lam(2)**2 + 9*lam(4)**2 + 18*lam(4)*lam(5) +               &
                  9*lam(5)**2 - 18*lam(4)*lam(6) - 18*lam(5)*lam(6) +           &
                  9*lam(6)**2 - 882*lam(4)*lam(7) - 882*lam(5)*lam(7) +         &
                  882*lam(6)*lam(7) + 11241*lam(7)**2 - 144*lam(8)**2 -         &
                  72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +          &
                     2*lam(8)) -                                                &
                  36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                &
                     10*lam(8)) + 132*lam(4)*lam(9) +                           &
                  132*lam(5)*lam(9) - 132*lam(6)*lam(9) -                       &
                  4740*lam(7)*lam(9) + 412*lam(9)**2 +                          &
                  96*lam(2)*(3*lam(4) + 3*lam(5) -                              &
                     3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) +           &
                  24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +               &
                     10*lam(9))*lam(10) - 144*lam(10)**2))) -                   &
         2*expdISs*qxy2*(1 + Rqsd)*tanhds**2*                                   &
          (3*dRqru*qdB*qxy2*(1 + Rqru)**2*                                      &
             (576*lam(2)**2 + 9*lam(4)**2 + 18*lam(4)*lam(5) +                  &
               9*lam(5)**2 - 18*lam(4)*lam(6) - 18*lam(5)*lam(6) +              &
               9*lam(6)**2 - 882*lam(4)*lam(7) - 882*lam(5)*lam(7) +            &
               882*lam(6)*lam(7) + 11241*lam(7)**2 - 144*lam(8)**2 +            &
               72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +             &
                  2*lam(8)) +                                                   &
               36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                   &
                  10*lam(8)) + 132*lam(4)*lam(9) + 132*lam(5)*lam(9) -          &
               132*lam(6)*lam(9) - 4740*lam(7)*lam(9) + 412*lam(9)**2 +         &
               96*lam(2)*(3*lam(4) + 3*lam(5) -                                 &
                  3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) +              &
               24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                  &
                  10*lam(9))*lam(10) - 144*lam(10)**2) +                        &
            dRqsd*quB*(6*m2*(-32*                                               &
                   (9*lam(1)**2 - 48*lam(1)*lam(3) + 10*lam(3)**2) +            &
                  3*(lam(4)**2 - 11*lam(5)**2 + lam(6)**2 -                     &
                     2*lam(4)*(7*lam(5) + lam(6) - 7*lam(7)) -                  &
                     14*lam(6)*lam(7) - 11*lam(7)**2 +                          &
                     2*lam(5)*(7*lam(6) + 11*lam(7))) +                         &
                  192*(lam(4) + 3*lam(5) - lam(6) - 3*lam(7))*lam(9) -          &
                  320*lam(9)**2) -                                              &
               384*m2*(lam(4) + 5*lam(5) - lam(6) -                             &
                  5*(lam(7) + lam(9)))*lam(10) - 448*m2*lam(10)**2 +            &
               3*qzt2*(1 + Rqru)*(1 + Rqsd)*                                    &
                (144*lam(2)**2 +                                                &
                  216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +          &
                  9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                     &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                 &
                     10*lam(6)*lam(7) + 13*lam(7)**2 +                          &
                     lam(5)*(-10*lam(6) + 26*lam(7))) +                         &
                  288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -          &
                  192*lam(8)**2 +                                               &
                  128*lam(3)*                                                   &
                   (3*lam(4) + 6*lam(5) -                                       &
                     3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) -            &
                  48*lam(1)*(3*lam(4) + 33*lam(5) -                             &
                     3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10)))))        &
+ qdB*(dRqsu*expuISs*tanhus*(3*qzt2*(1 + Rqru)*                                 &
                (quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +                           &
                  2*qxy2*(1 + Rqsu)**2*tanhus)*                                 &
                (3*(48*lam(2)**2 -                                              &
                     192*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) -           &
                     64*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +            &
                     72*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +            &
                     3*(lam(4) + lam(5) + lam(6) + lam(7))**2) -                &
                  256*(3*lam(1) + lam(3))*lam(10)) -                            &
               4*m2*(1 + Rqsu)*(quB*qzt2*ss + qxy2*tanhus)*                     &
                (384*(3*lam(1) + lam(3))**2 +                                   &
                  45*(lam(4) - lam(5) - lam(6) + lam(7))**2 +                   &
                  192*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +             &
                  32*lam(10)**2)) +                                             &
            6*dRqru*expdISs*qxy2*(1 + Rqru)*                                    &
             ((qxy2 + qzt2)*(1 + Rqru)*(1 + Rqsd)*                              &
                (576*lam(2)**2 + 9*lam(4)**2 + 18*lam(4)*lam(5) +               &
                  9*lam(5)**2 - 18*lam(4)*lam(6) - 18*lam(5)*lam(6) +           &
                  9*lam(6)**2 - 882*lam(4)*lam(7) - 882*lam(5)*lam(7) +         &
                  882*lam(6)*lam(7) + 11241*lam(7)**2 - 144*lam(8)**2 +         &
                  72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +          &
                     2*lam(8)) +                                                &
                  36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                &
                     10*lam(8)) + 132*lam(4)*lam(9) +                           &
                  132*lam(5)*lam(9) - 132*lam(6)*lam(9) -                       &
                  4740*lam(7)*lam(9) + 412*lam(9)**2 +                          &
                  96*lam(2)*(3*lam(4) + 3*lam(5) -                              &
                     3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) +           &
                  24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +               &
                     10*lam(9))*lam(10) - 144*lam(10)**2) -                     &
               12*m2*(72*lam(1)**2 + 96*lam(2)**2 + 24*lam(1)*lam(3) +          &
                  30*lam(3)**2 + 9*lam(4)**2 + 9*(lam(5) - lam(6))**2 -         &
                  90*lam(4)*lam(7) - 72*lam(5)*lam(7) +                         &
                  72*lam(6)*lam(7) + 945*lam(7)**2 - 12*lam(5)*lam(8) +         &
                  12*lam(6)*lam(8) + 30*lam(8)**2 - 384*lam(7)*lam(9) -         &
                  20*lam(8)*lam(9) + 16*lam(9)**2 + 24*lam(8)*lam(10) +         &
                  4*lam(2)*(3*lam(4) + 3*lam(5) - 3*lam(6) -                    &
                     147*lam(7) + 22*lam(9) + 12*lam(10)))) +                   &
            dRqru*expuISs*(3*(1 + Rqsu)*                                        &
                (-(quB*qzt2*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*tanhus*              &
                     (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +          &
                       117*lam(5)**2 - 18*lam(4)*lam(6) +                       &
                       90*lam(5)*lam(6) + 9*lam(6)**2 +                         &
                       522*lam(4)*lam(7) + 198*lam(5)*lam(7) -                  &
                       522*lam(6)*lam(7) - 5499*lam(7)**2 -                     &
                       36*lam(3)*                                               &
                        (7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +         &
                       288*lam(8)**2 -                                          &
                       72*lam(1)*                                               &
                        (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +            &
                        2*lam(8)) - 156*lam(4)*lam(9) +                         &
                       204*lam(5)*lam(9) + 156*lam(6)*lam(9) +                  &
                       2100*lam(7)*lam(9) - 56*lam(9)**2 -                      &
                       48*lam(2)*                                               &
                        (3*(lam(4) + lam(5) - lam(6) - 25*lam(7) +              &
                        4*lam(8)) + 16*lam(9)) -                                &
                       36*(5*lam(4) - 11*lam(5) - 5*lam(6) +                    &
                        11*lam(7) - 10*lam(9))*lam(10) + 456*lam(10)**2)        &
) - 2*qxy2*qzt2*(1 + Rqru)**2*                                                  &
                   (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +            &
                     117*lam(5)**2 - 18*lam(4)*lam(6) +                         &
                     90*lam(5)*lam(6) + 9*lam(6)**2 +                           &
                     522*lam(4)*lam(7) + 198*lam(5)*lam(7) -                    &
                     522*lam(6)*lam(7) - 5499*lam(7)**2 +                       &
                     36*lam(3)*                                                 &
                      (7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +           &
                     288*lam(8)**2 +                                            &
                     72*lam(1)*                                                 &
                      (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +              &
                        2*lam(8)) - 156*lam(4)*lam(9) +                         &
                     204*lam(5)*lam(9) + 156*lam(6)*lam(9) +                    &
                     2100*lam(7)*lam(9) - 56*lam(9)**2 -                        &
                     48*lam(2)*                                                 &
                      (3*(lam(4) + lam(5) - lam(6) - 25*lam(7) +                &
                        4*lam(8)) + 16*lam(9)) -                                &
                     36*(5*lam(4) - 11*lam(5) - 5*lam(6) +                      &
                        11*lam(7) - 10*lam(9))*lam(10) + 456*lam(10)**2)        &
+ 2*qxy2**2*(1 + Rqru)**2*(-1 + tanhus**2)*                                     &
                   (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +            &
                     117*lam(5)**2 - 18*lam(4)*lam(6) +                         &
                     90*lam(5)*lam(6) + 9*lam(6)**2 +                           &
                     522*lam(4)*lam(7) + 198*lam(5)*lam(7) -                    &
                     522*lam(6)*lam(7) - 5499*lam(7)**2 +                       &
                     36*lam(3)*                                                 &
                      (7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +           &
                     288*lam(8)**2 +                                            &
                     72*lam(1)*                                                 &
                      (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +              &
                        2*lam(8)) - 156*lam(4)*lam(9) +                         &
                     204*lam(5)*lam(9) + 156*lam(6)*lam(9) +                    &
                     2100*lam(7)*lam(9) - 56*lam(9)**2 -                        &
                     48*lam(2)*                                                 &
                      (3*(lam(4) + lam(5) - lam(6) - 25*lam(7) +                &
                        4*lam(8)) + 16*lam(9)) -                                &
                     36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -          &
                        10*lam(9))*lam(10) + 456*lam(10)**2)) +                 &
               4*m2*(1 + Rqru)*                                                 &
                (18*qxy2*(72*lam(1)**2 - 48*lam(2)**2 +                         &
                     168*lam(1)*lam(3) + 60*lam(3)**2 + 9*lam(4)**2 +           &
                     117*lam(5)**2 + 90*lam(5)*lam(6) + 9*lam(6)**2 +           &
                     126*lam(4)*lam(7) + 36*lam(5)*lam(7) -                     &
                     36*lam(6)*lam(7) - 351*lam(7)**2 +                         &
                     60*lam(5)*lam(8) + 12*lam(6)*lam(8) -                      &
                     12*lam(8)**2 + 192*lam(7)*lam(9) +                         &
                     40*lam(8)*lam(9) + 16*lam(9)**2 +                          &
                     4*lam(2)*                                                  &
                      (3*lam(4) - 15*lam(5) - 3*lam(6) + 87*lam(7) -            &
                        26*lam(9) - 30*lam(10)) +                               &
                     12*(3*lam(8) - 8*lam(9))*lam(10) + 96*lam(10)**2) -        &
                  quB*qzt2*sr*tanhus*                                           &
                   (864*lam(2)**2 + 45*lam(4)**2 + 369*lam(5)**2 +              &
                     234*lam(5)*lam(6) + 45*lam(6)**2 -                         &
                     1386*lam(5)*lam(7) + 414*lam(6)*lam(7) +                   &
                     8793*lam(7)**2 - 540*lam(5)*lam(8) -                       &
                     108*lam(6)*lam(8) + 540*lam(7)*lam(8) -                    &
                     6*lam(4)*                                                  &
                      (39*lam(5) + 15*lam(6) + 69*lam(7) - 18*lam(8) -          &
                        32*lam(9)) + 192*lam(5)*lam(9) -                        &
                     192*lam(6)*lam(9) - 3648*lam(7)*lam(9) -                   &
                     720*lam(8)*lam(9) + 608*lam(9)**2 -                        &
                     72*lam(2)*                                                 &
                      (3*lam(4) - 15*lam(5) - 3*lam(6) + 87*lam(7) -            &
                        26*lam(9) - 30*lam(10)) -                               &
                     72*(9*lam(8) + 8*lam(9))*lam(10) + 576*lam(10)**2))))      &
) - expuISr*qxy2*tanhur**2*(-6*dRqsd*expdISs*quB*qxy2*(1 + Rqru)*               &
          (1 + Rqsd)**2*tanhds**3*                                              &
          (144*lam(2)**2 + 216*lam(2)*                                          &
             (lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +                          &
            9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                           &
               2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                       &
               10*lam(6)*lam(7) + 13*lam(7)**2 +                                &
               lam(5)*(-10*lam(6) + 26*lam(7))) +                               &
            288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -                &
            192*lam(8)**2 - 128*lam(3)*                                         &
             (3*lam(4) + 6*lam(5) - 3*(lam(6) + 2*lam(7) + 4*lam(9)) +          &
               7*lam(10)) + 48*lam(1)*                                          &
             (3*lam(4) + 33*lam(5) - 3*(lam(6) + 11*lam(7) + 16*lam(9)) +       &
               16*lam(10))) + 3*expdISs*qdB*quB*tanhds**2*                      &
          (dRqru*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                     &
             (576*lam(2)**2 + 9*lam(4)**2 + 18*lam(4)*lam(5) +                  &
               9*lam(5)**2 - 18*lam(4)*lam(6) - 18*lam(5)*lam(6) +              &
               9*lam(6)**2 - 882*lam(4)*lam(7) - 882*lam(5)*lam(7) +            &
               882*lam(6)*lam(7) + 11241*lam(7)**2 - 144*lam(8)**2 +            &
               72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +             &
                  2*lam(8)) +                                                   &
               36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                   &
                  10*lam(8)) + 132*lam(4)*lam(9) + 132*lam(5)*lam(9) -          &
               132*lam(6)*lam(9) - 4740*lam(7)*lam(9) + 412*lam(9)**2 +         &
               96*lam(2)*(3*lam(4) + 3*lam(5) -                                 &
                  3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) +              &
               24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                  &
                  10*lam(9))*lam(10) - 144*lam(10)**2) -                        &
            dRqsd*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                    &
             (144*lam(2)**2 +                                                   &
               216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +             &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                    &
                  10*lam(6)*lam(7) + 13*lam(7)**2 +                             &
                  lam(5)*(-10*lam(6) + 26*lam(7))) +                            &
               288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -             &
               192*lam(8)**2 -                                                  &
               128*lam(3)*(3*lam(4) + 6*lam(5) -                                &
                  3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) +               &
               48*lam(1)*(3*lam(4) + 33*lam(5) -                                &
                  3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10)))) +          &
         2*expdISs*(1 + Rqru)*tanhds*                                           &
          (3*dRqsd*quB*qxy2*(1 + Rqsd)**2*                                      &
             (144*lam(2)**2 +                                                   &
               216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +             &
               9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                        &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                    &
                  10*lam(6)*lam(7) + 13*lam(7)**2 +                             &
                  lam(5)*(-10*lam(6) + 26*lam(7))) +                            &
               288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -             &
               192*lam(8)**2 -                                                  &
               128*lam(3)*(3*lam(4) + 6*lam(5) -                                &
                  3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) +               &
               48*lam(1)*(3*lam(4) + 33*lam(5) -                                &
                  3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10))) +           &
            dRqru*qdB*(2*m2*(9*                                                 &
                   (-192*lam(2)**2 - (lam(4) + lam(5) - lam(6))**2 -            &
                     24*lam(2)*(lam(4) + lam(5) - lam(6) - 49*lam(7)) +         &
                     146*(lam(4) + lam(5) - lam(6))*lam(7) -                    &
                     1873*lam(7)**2 +                                           &
                     12*(lam(4) + lam(5) - lam(6) - lam(7))*lam(8)) -           &
                  24*(66*lam(2) +                                               &
                     16*(lam(4) + lam(5) - lam(6) - 19*lam(7)) -                &
                     15*lam(8))*lam(9) - 928*lam(9)**2 -                        &
                  432*(2*lam(2) + lam(8))*lam(10)) +                            &
               3*qzt2*(1 + Rqru)*(1 + Rqsd)*                                    &
                (576*lam(2)**2 + 9*lam(4)**2 + 18*lam(4)*lam(5) +               &
                  9*lam(5)**2 - 18*lam(4)*lam(6) - 18*lam(5)*lam(6) +           &
                  9*lam(6)**2 - 882*lam(4)*lam(7) - 882*lam(5)*lam(7) +         &
                  882*lam(6)*lam(7) + 11241*lam(7)**2 - 144*lam(8)**2 -         &
                  72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +          &
                     2*lam(8)) -                                                &
                  36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                &
                     10*lam(8)) + 132*lam(4)*lam(9) + 132*lam(5)*lam(9) -       &
                  132*lam(6)*lam(9) - 4740*lam(7)*lam(9) +                      &
                  412*lam(9)**2 +                                               &
                  96*lam(2)*(3*lam(4) + 3*lam(5) -                              &
                     3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) +           &
                  24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +               &
                     10*lam(9))*lam(10) - 144*lam(10)**2))) +                   &
         qdB*(3*(1 + Rqru)*(dRqsu*expuISs*(-1 + tanhus)*(1 + tanhus)*           &
                (quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +                           &
                  2*qxy2*(1 + Rqsu)**2*tanhus)*                                 &
                (3*(48*lam(2)**2 +                                              &
                     192*lam(1)*(lam(4) - lam(5) - lam(6) + lam(7)) +           &
                     64*lam(3)*(lam(4) - lam(5) - lam(6) + lam(7)) +            &
                     72*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +            &
                     3*(lam(4) + lam(5) + lam(6) + lam(7))**2) +                &
                  256*(3*lam(1) + lam(3))*lam(10)) +                            &
               dRqsd*expdISs*quB*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                &
                (144*lam(2)**2 +                                                &
                  216*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7)) +          &
                  9*(lam(4)**2 + 13*lam(5)**2 + lam(6)**2 +                     &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) -                 &
                     10*lam(6)*lam(7) + 13*lam(7)**2 +                          &
                     lam(5)*(-10*lam(6) + 26*lam(7))) +                         &
                  288*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) -          &
                  192*lam(8)**2 -                                               &
                  128*lam(3)*(3*lam(4) + 6*lam(5) -                             &
                     3*(lam(6) + 2*lam(7) + 4*lam(9)) + 7*lam(10)) +            &
                  48*lam(1)*(3*lam(4) + 33*lam(5) -                             &
                     3*(lam(6) + 11*lam(7) + 16*lam(9)) + 16*lam(10)))) -       &
            dRqru*(3*expdISs*quB*(1 + Rqsd)*                                    &
                (-1 + 2*qzt2*(1 + Rqru)**2*sr)*                                 &
                (576*lam(2)**2 + 9*lam(4)**2 + 18*lam(4)*lam(5) +               &
                  9*lam(5)**2 - 18*lam(4)*lam(6) - 18*lam(5)*lam(6) +           &
                  9*lam(6)**2 - 882*lam(4)*lam(7) - 882*lam(5)*lam(7) +         &
                  882*lam(6)*lam(7) + 11241*lam(7)**2 - 144*lam(8)**2 +         &
                  72*lam(1)*(3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +          &
                     2*lam(8)) +                                                &
                  36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7) +                &
                     10*lam(8)) + 132*lam(4)*lam(9) + 132*lam(5)*lam(9) -       &
                  132*lam(6)*lam(9) - 4740*lam(7)*lam(9) +                      &
                  412*lam(9)**2 +                                               &
                  96*lam(2)*(3*lam(4) + 3*lam(5) -                              &
                     3*(lam(6) + 25*lam(7) + 2*lam(8)) + 16*lam(9)) +           &
                  24*(3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +               &
                     10*lam(9))*lam(10) - 144*lam(10)**2) +                     &
               expuISs*(3*quB*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*        &
                   (-1 + tanhus**2)*                                            &
                   (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +            &
                     117*lam(5)**2 - 18*lam(4)*lam(6) +                         &
                     90*lam(5)*lam(6) + 9*lam(6)**2 + 522*lam(4)*lam(7) +       &
                     198*lam(5)*lam(7) - 522*lam(6)*lam(7) -                    &
                     5499*lam(7)**2 +                                           &
                     36*lam(3)*                                                 &
                      (7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +           &
                     288*lam(8)**2 +                                            &
                     72*lam(1)*                                                 &
                      (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +              &
                        2*lam(8)) - 156*lam(4)*lam(9) +                         &
                     204*lam(5)*lam(9) + 156*lam(6)*lam(9) +                    &
                     2100*lam(7)*lam(9) - 56*lam(9)**2 -                        &
                     48*lam(2)*                                                 &
                      (3*(lam(4) + lam(5) - lam(6) - 25*lam(7) +                &
                         4*lam(8)) + 16*lam(9)) -                               &
                     36*(5*lam(4) - 11*lam(5) - 5*lam(6) + 11*lam(7) -          &
                        10*lam(9))*lam(10) + 456*lam(10)**2) +                  &
                  2*(1 + Rqru)*tanhus*                                          &
                   (3*qzt2*(1 + Rqru)*(1 + Rqsu)*                               &
                      (576*lam(2)**2 + 9*lam(4)**2 - 90*lam(4)*lam(5) +         &
                        117*lam(5)**2 - 18*lam(4)*lam(6) +                      &
                        90*lam(5)*lam(6) + 9*lam(6)**2 +                        &
                        522*lam(4)*lam(7) + 198*lam(5)*lam(7) -                 &
                        522*lam(6)*lam(7) - 5499*lam(7)**2 -                    &
                        36*lam(3)*                                              &
                         (7*lam(4) + 23*lam(5) + 7*lam(6) + 23*lam(7)) +        &
                        288*lam(8)**2 -                                         &
                        72*lam(1)*                                              &
                         (3*(lam(4) + 5*lam(5) + lam(6) + 5*lam(7)) +           &
                         2*lam(8)) - 156*lam(4)*lam(9) +                        &
                        204*lam(5)*lam(9) + 156*lam(6)*lam(9) +                 &
                        2100*lam(7)*lam(9) - 56*lam(9)**2 -                     &
                        48*lam(2)*                                              &
                         (3*(lam(4) + lam(5) - lam(6) - 25*lam(7) +             &
                         4*lam(8)) + 16*lam(9)) -                               &
                        36*(5*lam(4) - 11*lam(5) - 5*lam(6) +                   &
                         11*lam(7) - 10*lam(9))*lam(10) + 456*lam(10)**2)       &
+ 2*m2*(864*lam(2)**2 + 45*lam(4)**2 + 369*lam(5)**2 + 234*lam(5)*lam(6) +      &
                        45*lam(6)**2 - 1386*lam(5)*lam(7) +                     &
                        414*lam(6)*lam(7) + 8793*lam(7)**2 -                    &
                        540*lam(5)*lam(8) - 108*lam(6)*lam(8) +                 &
                        540*lam(7)*lam(8) -                                     &
                        6*lam(4)*                                               &
                         (39*lam(5) + 15*lam(6) + 69*lam(7) - 18*lam(8) -       &
                         32*lam(9)) + 192*lam(5)*lam(9) -                       &
                        192*lam(6)*lam(9) - 3648*lam(7)*lam(9) -                &
                        720*lam(8)*lam(9) + 608*lam(9)**2 -                     &
                        72*lam(2)*                                              &
                         (3*lam(4) - 15*lam(5) - 3*lam(6) + 87*lam(7) -         &
                         26*lam(9) - 30*lam(10)) -                              &
                        72*(9*lam(8) + 8*lam(9))*lam(10) + 576*lam(10)**2))     &
)))))/3888.
  dtLAM(8)=                                                        &
   (expdISr*expuISs*qdB*(6*dRqsu*qxy2**2*(1 + Rqrd)*(1 + Rqsu)**2*              &
          tanhus**3*(48*lam(2)**2 +                                             &
            3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                    &
               lam(6)**2 + 2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +           &
               26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +            &
            144*lam(2)*lam(8) -                                                 &
            12*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) +                 &
            80*lam(8)**2) - 6*dRqsu*qxy2*(1 + Rqsu)*tanhus*                     &
          (24*m2*(2*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7) -             &
                  2*lam(8)) +                                                   &
               (3*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7)) - 4*lam(8))*          &
                lam(8)) + qxy2*(1 + Rqrd)*(1 + Rqsu)*                           &
             (3*(16*lam(2)**2 + lam(4)**2 + 13*lam(5)**2 -                      &
                  10*lam(5)*lam(6) + lam(6)**2 +                                &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                    &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               12*(12*lam(2) - lam(4) + 3*lam(5) - lam(6) + 3*lam(7))*          &
                lam(8) + 80*lam(8)**2) +                                        &
            qzt2*(1 + Rqrd)*(1 + Rqsu)*                                         &
             (3*(16*lam(2)**2 + lam(4)**2 + 13*lam(5)**2 -                      &
                  10*lam(5)*lam(6) + lam(6)**2 +                                &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                    &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               12*(12*lam(2) - lam(4) + 3*lam(5) - lam(6) + 3*lam(7))*          &
                lam(8) + 80*lam(8)**2)) +                                       &
         quB*(-3*dRqsu*(qxy2*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*         &
                (48*lam(2)**2 +                                                 &
                  3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +              &
                     lam(6)**2 +                                                &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                 &
                     26*lam(5)*lam(7) - 10*lam(6)*lam(7) +                      &
                     13*lam(7)**2) + 144*lam(2)*lam(8) -                        &
                  12*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) +           &
                  80*lam(8)**2) +                                               &
               qzt2*(48*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*              &
                   lam(2)**2 + 96*m2*(1 + Rqsu)*ss*lam(2)*lam(4) +              &
                  3*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                  &
                   lam(4)**2 - 480*m2*(1 + Rqsu)*ss*lam(2)*lam(5) -             &
                  39*lam(5)**2 - 39*Rqrd*lam(5)**2 +                            &
                  78*qzt2*ss*lam(5)**2 + 78*qzt2*Rqrd*ss*lam(5)**2 +            &
                  156*qzt2*Rqsu*ss*lam(5)**2 +                                  &
                  156*qzt2*Rqrd*Rqsu*ss*lam(5)**2 +                             &
                  78*qzt2*Rqsu**2*ss*lam(5)**2 +                                &
                  78*qzt2*Rqrd*Rqsu**2*ss*lam(5)**2 +                           &
                  96*m2*(1 + Rqsu)*ss*lam(2)*lam(6) +                           &
                  30*lam(5)*lam(6) + 30*Rqrd*lam(5)*lam(6) -                    &
                  60*qzt2*ss*lam(5)*lam(6) -                                    &
                  60*qzt2*Rqrd*ss*lam(5)*lam(6) -                               &
                  120*qzt2*Rqsu*ss*lam(5)*lam(6) -                              &
                  120*qzt2*Rqrd*Rqsu*ss*lam(5)*lam(6) -                         &
                  60*qzt2*Rqsu**2*ss*lam(5)*lam(6) -                            &
                  60*qzt2*Rqrd*Rqsu**2*ss*lam(5)*lam(6) - 3*lam(6)**2 -         &
                  3*Rqrd*lam(6)**2 + 6*qzt2*ss*lam(6)**2 +                      &
                  6*qzt2*Rqrd*ss*lam(6)**2 +                                    &
                  12*qzt2*Rqsu*ss*lam(6)**2 +                                   &
                  12*qzt2*Rqrd*Rqsu*ss*lam(6)**2 +                              &
                  6*qzt2*Rqsu**2*ss*lam(6)**2 +                                 &
                  6*qzt2*Rqrd*Rqsu**2*ss*lam(6)**2 -                            &
                  480*m2*(1 + Rqsu)*ss*lam(2)*lam(7) -                          &
                  78*lam(5)*lam(7) - 78*Rqrd*lam(5)*lam(7) +                    &
                  156*qzt2*ss*lam(5)*lam(7) +                                   &
                  156*qzt2*Rqrd*ss*lam(5)*lam(7) +                              &
                  312*qzt2*Rqsu*ss*lam(5)*lam(7) +                              &
                  312*qzt2*Rqrd*Rqsu*ss*lam(5)*lam(7) +                         &
                  156*qzt2*Rqsu**2*ss*lam(5)*lam(7) +                           &
                  156*qzt2*Rqrd*Rqsu**2*ss*lam(5)*lam(7) +                      &
                  30*lam(6)*lam(7) + 30*Rqrd*lam(6)*lam(7) -                    &
                  60*qzt2*ss*lam(6)*lam(7) -                                    &
                  60*qzt2*Rqrd*ss*lam(6)*lam(7) -                               &
                  120*qzt2*Rqsu*ss*lam(6)*lam(7) -                              &
                  120*qzt2*Rqrd*Rqsu*ss*lam(6)*lam(7) -                         &
                  60*qzt2*Rqsu**2*ss*lam(6)*lam(7) -                            &
                  60*qzt2*Rqrd*Rqsu**2*ss*lam(6)*lam(7) -                       &
                  39*lam(7)**2 - 39*Rqrd*lam(7)**2 +                            &
                  78*qzt2*ss*lam(7)**2 + 78*qzt2*Rqrd*ss*lam(7)**2 +            &
                  156*qzt2*Rqsu*ss*lam(7)**2 +                                  &
                  156*qzt2*Rqrd*Rqsu*ss*lam(7)**2 +                             &
                  78*qzt2*Rqsu**2*ss*lam(7)**2 +                                &
                  78*qzt2*Rqrd*Rqsu**2*ss*lam(7)**2 -                           &
                  6*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*lam(4)*           &
                   (5*lam(5) - lam(6) + 5*lam(7)) +                             &
                  (48*(-3*(1 + Rqrd) +                                          &
                       2*(1 + Rqsu)*                                            &
                       (-2*m2 + 3*qzt2*(1 + Rqrd)*(1 + Rqsu))*ss)*              &
                      lam(2) -                                                  &
                     12*(-1 - Rqrd +                                            &
                        2*(1 + Rqsu)*                                           &
                        (-6*m2 + qzt2*(1 + Rqrd)*(1 + Rqsu))*ss)*               &
                      (lam(4) - 3*lam(5) + lam(6) - 3*lam(7)))*lam(8) +         &
                  16*(-5*(1 + Rqrd) +                                           &
                     2*(1 + Rqsu)*                                              &
                      (-6*m2 + 5*qzt2*(1 + Rqrd)*(1 + Rqsu))*ss)*               &
                   lam(8)**2)) +                                                &
            dRqrd*(qzt2*(576*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*         &
                   lam(2)**2 +                                                  &
                  96*lam(2)*(6*(1 + Rqsu)*                                      &
                      (-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(8) -                   &
                     m2*(1 + Rqrd)*sr*                                          &
                      (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +              &
                        10*lam(9) - 12*lam(10))) +                              &
                  (1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                    &
                   (144*lam(8)*(2*lam(1) + lam(8)) +                            &
                     24*lam(3)*                                                 &
                      (3*(lam(4) - lam(5) + lam(6) - lam(7)) +                  &
                        14*lam(8)) +                                            &
                     (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +               &
                        10*lam(9) - 12*lam(10))**2) -                           &
                  48*m2*(1 + Rqrd)*sr*                                          &
                   (24*lam(1)*lam(3) + 14*lam(3)**2 +                           &
                     lam(8)*(9*lam(4) - 3*lam(5) + 3*lam(6) -                   &
                        9*lam(7) + 14*lam(8) + 10*lam(9) - 12*lam(10))))        &
+ qxy2*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                               &
                (576*lam(2)**2 + 576*lam(2)*lam(8) +                            &
                  144*lam(8)*(2*lam(1) + lam(8)) +                              &
                  24*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +            &
                     14*lam(8)) +                                               &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2))) -                            &
         quB*qxy2*tanhus**2*(-3*dRqsu*(1 + Rqrd)*                               &
             (-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                                    &
             (48*lam(2)**2 + 3*                                                 &
                (lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                  &
                  lam(6)**2 +                                                   &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                    &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               144*lam(2)*lam(8) -                                              &
               12*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) +              &
               80*lam(8)**2) +                                                  &
            dRqrd*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                    &
             (576*lam(2)**2 + 576*lam(2)*lam(8) +                               &
               144*lam(8)*(2*lam(1) + lam(8)) +                                 &
               24*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +               &
                  14*lam(8)) +                                                  &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -         &
                  12*lam(10))**2))) +                                           &
      expdISs*expuISr*quB*(6*dRqsd*qxy2**2*(1 + Rqru)*(1 + Rqsd)**2*            &
          tanhds**3*(48*lam(2)**2 +                                             &
            3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                    &
               lam(6)**2 + 2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +           &
               26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +            &
            144*lam(2)*lam(8) -                                                 &
            12*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) +                 &
            80*lam(8)**2) - 6*dRqsd*qxy2*(1 + Rqsd)*tanhds*                     &
          (24*m2*(2*lam(2)*(lam(4) - 5*lam(5) + lam(6) - 5*lam(7) -             &
                  2*lam(8)) +                                                   &
               (3*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7)) - 4*lam(8))*          &
                lam(8)) + qxy2*(1 + Rqru)*(1 + Rqsd)*                           &
             (3*(16*lam(2)**2 + lam(4)**2 + 13*lam(5)**2 -                      &
                  10*lam(5)*lam(6) + lam(6)**2 +                                &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                    &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               12*(12*lam(2) - lam(4) + 3*lam(5) - lam(6) + 3*lam(7))*          &
                lam(8) + 80*lam(8)**2) +                                        &
            qzt2*(1 + Rqru)*(1 + Rqsd)*                                         &
             (3*(16*lam(2)**2 + lam(4)**2 + 13*lam(5)**2 -                      &
                  10*lam(5)*lam(6) + lam(6)**2 +                                &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                    &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               12*(12*lam(2) - lam(4) + 3*lam(5) - lam(6) + 3*lam(7))*          &
                lam(8) + 80*lam(8)**2)) +                                       &
         qdB*(-3*dRqsd*(qxy2*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*         &
                (48*lam(2)**2 +                                                 &
                  3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +              &
                     lam(6)**2 +                                                &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                 &
                     26*lam(5)*lam(7) - 10*lam(6)*lam(7) +                      &
                     13*lam(7)**2) + 144*lam(2)*lam(8) -                        &
                  12*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) +           &
                  80*lam(8)**2) +                                               &
               qzt2*(48*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*              &
                   lam(2)**2 + 96*m2*(1 + Rqsd)*ss*lam(2)*lam(4) +              &
                  3*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                  &
                   lam(4)**2 - 480*m2*(1 + Rqsd)*ss*lam(2)*lam(5) -             &
                  39*lam(5)**2 - 39*Rqru*lam(5)**2 +                            &
                  78*qzt2*ss*lam(5)**2 + 78*qzt2*Rqru*ss*lam(5)**2 +            &
                  156*qzt2*Rqsd*ss*lam(5)**2 +                                  &
                  156*qzt2*Rqru*Rqsd*ss*lam(5)**2 +                             &
                  78*qzt2*Rqsd**2*ss*lam(5)**2 +                                &
                  78*qzt2*Rqru*Rqsd**2*ss*lam(5)**2 +                           &
                  96*m2*(1 + Rqsd)*ss*lam(2)*lam(6) +                           &
                  30*lam(5)*lam(6) + 30*Rqru*lam(5)*lam(6) -                    &
                  60*qzt2*ss*lam(5)*lam(6) -                                    &
                  60*qzt2*Rqru*ss*lam(5)*lam(6) -                               &
                  120*qzt2*Rqsd*ss*lam(5)*lam(6) -                              &
                  120*qzt2*Rqru*Rqsd*ss*lam(5)*lam(6) -                         &
                  60*qzt2*Rqsd**2*ss*lam(5)*lam(6) -                            &
                  60*qzt2*Rqru*Rqsd**2*ss*lam(5)*lam(6) - 3*lam(6)**2 -         &
                  3*Rqru*lam(6)**2 + 6*qzt2*ss*lam(6)**2 +                      &
                  6*qzt2*Rqru*ss*lam(6)**2 +                                    &
                  12*qzt2*Rqsd*ss*lam(6)**2 +                                   &
                  12*qzt2*Rqru*Rqsd*ss*lam(6)**2 +                              &
                  6*qzt2*Rqsd**2*ss*lam(6)**2 +                                 &
                  6*qzt2*Rqru*Rqsd**2*ss*lam(6)**2 -                            &
                  480*m2*(1 + Rqsd)*ss*lam(2)*lam(7) -                          &
                  78*lam(5)*lam(7) - 78*Rqru*lam(5)*lam(7) +                    &
                  156*qzt2*ss*lam(5)*lam(7) +                                   &
                  156*qzt2*Rqru*ss*lam(5)*lam(7) +                              &
                  312*qzt2*Rqsd*ss*lam(5)*lam(7) +                              &
                  312*qzt2*Rqru*Rqsd*ss*lam(5)*lam(7) +                         &
                  156*qzt2*Rqsd**2*ss*lam(5)*lam(7) +                           &
                  156*qzt2*Rqru*Rqsd**2*ss*lam(5)*lam(7) +                      &
                  30*lam(6)*lam(7) + 30*Rqru*lam(6)*lam(7) -                    &
                  60*qzt2*ss*lam(6)*lam(7) -                                    &
                  60*qzt2*Rqru*ss*lam(6)*lam(7) -                               &
                  120*qzt2*Rqsd*ss*lam(6)*lam(7) -                              &
                  120*qzt2*Rqru*Rqsd*ss*lam(6)*lam(7) -                         &
                  60*qzt2*Rqsd**2*ss*lam(6)*lam(7) -                            &
                  60*qzt2*Rqru*Rqsd**2*ss*lam(6)*lam(7) -                       &
                  39*lam(7)**2 - 39*Rqru*lam(7)**2 +                            &
                  78*qzt2*ss*lam(7)**2 + 78*qzt2*Rqru*ss*lam(7)**2 +            &
                  156*qzt2*Rqsd*ss*lam(7)**2 +                                  &
                  156*qzt2*Rqru*Rqsd*ss*lam(7)**2 +                             &
                  78*qzt2*Rqsd**2*ss*lam(7)**2 +                                &
                  78*qzt2*Rqru*Rqsd**2*ss*lam(7)**2 -                           &
                  6*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*lam(4)*           &
                   (5*lam(5) - lam(6) + 5*lam(7)) +                             &
                  (48*(-3*(1 + Rqru) +                                          &
                       2*(1 + Rqsd)*                                            &
                       (-2*m2 + 3*qzt2*(1 + Rqru)*(1 + Rqsd))*ss)*              &
                      lam(2) -                                                  &
                     12*(-1 - Rqru +                                            &
                        2*(1 + Rqsd)*                                           &
                        (-6*m2 + qzt2*(1 + Rqru)*(1 + Rqsd))*ss)*               &
                      (lam(4) - 3*lam(5) + lam(6) - 3*lam(7)))*lam(8) +         &
                  16*(-5*(1 + Rqru) +                                           &
                     2*(1 + Rqsd)*                                              &
                      (-6*m2 + 5*qzt2*(1 + Rqru)*(1 + Rqsd))*ss)*               &
                   lam(8)**2)) +                                                &
            dRqru*(qzt2*(576*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*         &
                   lam(2)**2 +                                                  &
                  96*lam(2)*(6*(1 + Rqsd)*                                      &
                      (-1 + 2*qzt2*(1 + Rqru)**2*sr)*lam(8) -                   &
                     m2*(1 + Rqru)*sr*                                          &
                      (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +              &
                        10*lam(9) - 12*lam(10))) +                              &
                  (1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                    &
                   (144*lam(8)*(2*lam(1) + lam(8)) +                            &
                     24*lam(3)*                                                 &
                      (3*(lam(4) - lam(5) + lam(6) - lam(7)) +                  &
                        14*lam(8)) +                                            &
                     (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +               &
                        10*lam(9) - 12*lam(10))**2) -                           &
                  48*m2*(1 + Rqru)*sr*                                          &
                   (24*lam(1)*lam(3) + 14*lam(3)**2 +                           &
                     lam(8)*(9*lam(4) - 3*lam(5) + 3*lam(6) -                   &
                        9*lam(7) + 14*lam(8) + 10*lam(9) - 12*lam(10))))        &
+ qxy2*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                               &
                (576*lam(2)**2 + 576*lam(2)*lam(8) +                            &
                  144*lam(8)*(2*lam(1) + lam(8)) +                              &
                  24*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +            &
                     14*lam(8)) +                                               &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2))) -                            &
         qdB*qxy2*tanhds**2*(-3*dRqsd*(1 + Rqru)*                               &
             (-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                                    &
             (48*lam(2)**2 + 3*                                                 &
                (lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                  &
                  lam(6)**2 +                                                   &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                    &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               144*lam(2)*lam(8) -                                              &
               12*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) +              &
               80*lam(8)**2) +                                                  &
            dRqru*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                    &
             (576*lam(2)**2 + 576*lam(2)*lam(8) +                               &
               144*lam(8)*(2*lam(1) + lam(8)) +                                 &
               24*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +               &
                  14*lam(8)) +                                                  &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -         &
                  12*lam(10))**2))) -                                           &
      2*dRqrd*expdISr*quB*qxy2**2*(1 + Rqrd)**2*tanhdr**3*                      &
       (-(expuISs*(1 + Rqsu)*(-1 + tanhus**2)*                                  &
            (576*lam(2)**2 + 576*lam(2)*lam(8) +                                &
              144*lam(8)*(2*lam(1) + lam(8)) +                                  &
              24*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +                &
                 14*lam(8)) +                                                   &
              (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -          &
                 12*lam(10))**2)) +                                             &
         expdISs*(1 + Rqsd)*(-1 + tanhds**2)*                                   &
          (9*lam(4)**2 + 117*lam(5)**2 + 90*lam(5)*lam(6) + 9*lam(6)**2 -       &
            36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -                     &
            234*lam(5)*lam(7) - 90*lam(6)*lam(7) + 117*lam(7)**2 -              &
            144*lam(1)*lam(8) +                                                 &
            24*(24*lam(2)**2 - 48*lam(2)*lam(8) +                               &
               lam(8)*(-7*lam(3) + 6*lam(8))) + 240*lam(5)*lam(9) +             &
            120*lam(6)*lam(9) - 240*lam(7)*lam(9) + 100*lam(9)**2 +             &
            12*(39*lam(5) + 9*lam(6) - 39*lam(7) + 50*lam(9))*lam(10) +         &
            312*lam(10)**2 - 6*lam(4)*                                          &
             (15*lam(5) + 3*lam(6) - 15*lam(7) + 20*lam(9) + 18*lam(10))))      &
+ 2*dRqru*expuISr*qdB*qxy2**2*(1 + Rqru)**2*tanhur**3*                          &
       (expdISs*(1 + Rqsd)*(-1 + tanhds**2)*                                    &
          (576*lam(2)**2 + 576*lam(2)*lam(8) +                                  &
            144*lam(8)*(2*lam(1) + lam(8)) +                                    &
            24*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +                  &
               14*lam(8)) + (3*lam(4) + 3*lam(5) - 3*lam(6) -                   &
               3*lam(7) + 10*lam(9) - 12*lam(10))**2) -                         &
         expuISs*(1 + Rqsu)*(-1 + tanhus**2)*                                   &
          (9*lam(4)**2 + 117*lam(5)**2 + 90*lam(5)*lam(6) + 9*lam(6)**2 -       &
            36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -                     &
            234*lam(5)*lam(7) - 90*lam(6)*lam(7) + 117*lam(7)**2 -              &
            144*lam(1)*lam(8) +                                                 &
            24*(24*lam(2)**2 - 48*lam(2)*lam(8) +                               &
               lam(8)*(-7*lam(3) + 6*lam(8))) + 240*lam(5)*lam(9) +             &
            120*lam(6)*lam(9) - 240*lam(7)*lam(9) + 100*lam(9)**2 +             &
            12*(39*lam(5) + 9*lam(6) - 39*lam(7) + 50*lam(9))*lam(10) +         &
            312*lam(10)**2 - 6*lam(4)*                                          &
             (15*lam(5) + 3*lam(6) - 15*lam(7) + 20*lam(9) + 18*lam(10))))      &
- expdISr*qxy2*tanhdr**2*(-18*dRqsd*expdISs*quB*qxy2*(1 + Rqrd)*                &
          (1 + Rqsd)**2*tanhds**3*                                              &
          (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +             &
         6*dRqsu*expuISs*qdB*qxy2*(1 + Rqrd)*(1 + Rqsu)**2*tanhus**3*           &
          (48*lam(2)**2 + 3*(lam(4)**2 + 13*lam(5)**2 -                         &
               10*lam(5)*lam(6) + lam(6)**2 +                                   &
               2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                       &
               26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +            &
            144*lam(2)*lam(8) -                                                 &
            12*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) +                 &
            80*lam(8)**2) - 2*expuISs*(1 + Rqrd)*tanhus*                        &
          (3*dRqsu*qdB*qxy2*(1 + Rqsu)**2*                                      &
             (48*lam(2)**2 +                                                    &
               3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                 &
                  lam(6)**2 +                                                   &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                    &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               144*lam(2)*lam(8) -                                              &
               12*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) +              &
               80*lam(8)**2) +                                                  &
            dRqrd*quB*(qzt2*(1 + Rqrd)*(1 + Rqsu)*                              &
                (576*lam(2)**2 + 576*lam(2)*lam(8) +                            &
                  144*lam(8)*(-2*lam(1) + lam(8)) -                             &
                  24*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +            &
                     14*lam(8)) +                                               &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) -                              &
               24*m2*(2*lam(2) + lam(8))*                                       &
                (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -        &
                  12*lam(10)))) -                                               &
         expuISs*qdB*quB*tanhus**2*                                             &
          (-3*dRqsu*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                  &
             (48*lam(2)**2 +                                                    &
               3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                 &
                  lam(6)**2 +                                                   &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                    &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               144*lam(2)*lam(8) -                                              &
               12*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) +              &
               80*lam(8)**2) +                                                  &
            dRqrd*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                    &
             (576*lam(2)**2 + 576*lam(2)*lam(8) +                               &
               144*lam(8)*(2*lam(1) + lam(8)) +                                 &
               24*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +               &
                  14*lam(8)) +                                                  &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -         &
                  12*lam(10))**2)) +                                            &
         expdISs*qdB*quB*tanhds**2*                                             &
          (-9*dRqsd*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                  &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            dRqrd*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                    &
             (9*lam(4)**2 + 117*lam(5)**2 + 90*lam(5)*lam(6) +                  &
               9*lam(6)**2 - 36*lam(3)*                                         &
                (lam(4) - lam(5) + lam(6) - lam(7)) -                           &
               234*lam(5)*lam(7) - 90*lam(6)*lam(7) + 117*lam(7)**2 -           &
               144*lam(1)*lam(8) +                                              &
               24*(24*lam(2)**2 - 48*lam(2)*lam(8) +                            &
                  lam(8)*(-7*lam(3) + 6*lam(8))) + 240*lam(5)*lam(9) +          &
               120*lam(6)*lam(9) - 240*lam(7)*lam(9) + 100*lam(9)**2 +          &
               12*(39*lam(5) + 9*lam(6) - 39*lam(7) + 50*lam(9))*               &
                lam(10) + 312*lam(10)**2 -                                      &
               6*lam(4)*(15*lam(5) + 3*lam(6) - 15*lam(7) + 20*lam(9) +         &
                  18*lam(10)))) +                                               &
         qdB*quB*(3*(1 + Rqrd)*                                                 &
             (3*dRqsd*expdISs*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                   &
                (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2)         &
- dRqsu*expuISs*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                                 &
                (3*(16*lam(2)**2 + lam(4)**2 + 13*lam(5)**2 -                   &
                     10*lam(5)*lam(6) + lam(6)**2 +                             &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                 &
                     26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2         &
) + 12*(12*lam(2) - lam(4) + 3*lam(5) - lam(6) + 3*lam(7))*lam(8) +             &
                  80*lam(8)**2)) -                                              &
            dRqrd*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                               &
             (-(expuISs*(1 + Rqsu)*                                             &
                  (576*lam(2)**2 + 576*lam(2)*lam(8) +                          &
                    144*lam(8)*(2*lam(1) + lam(8)) +                            &
                    24*lam(3)*                                                  &
                     (3*(lam(4) - lam(5) + lam(6) - lam(7)) +                   &
                       14*lam(8)) +                                             &
                    (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                &
                       10*lam(9) - 12*lam(10))**2)) +                           &
               expdISs*(1 + Rqsd)*                                              &
                (9*lam(4)**2 + 117*lam(5)**2 + 90*lam(5)*lam(6) +               &
                  9*lam(6)**2 -                                                 &
                  36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -               &
                  234*lam(5)*lam(7) - 90*lam(6)*lam(7) + 117*lam(7)**2 -        &
                  144*lam(1)*lam(8) +                                           &
                  24*(24*lam(2)**2 - 48*lam(2)*lam(8) +                         &
                     lam(8)*(-7*lam(3) + 6*lam(8))) +                           &
                  240*lam(5)*lam(9) + 120*lam(6)*lam(9) -                       &
                  240*lam(7)*lam(9) + 100*lam(9)**2 +                           &
                  12*(39*lam(5) + 9*lam(6) - 39*lam(7) + 50*lam(9))*            &
                   lam(10) + 312*lam(10)**2 -                                   &
                  6*lam(4)*(15*lam(5) + 3*lam(6) - 15*lam(7) +                  &
                     20*lam(9) + 18*lam(10))))) +                               &
         2*expdISs*quB*(1 + Rqrd)*tanhds*                                       &
          (9*dRqsd*qxy2*(1 + Rqsd)**2*                                          &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            dRqrd*(-48*m2*(lam(2)*                                              &
                   (3*lam(4) - 15*lam(5) - 3*lam(6) + 15*lam(7) -               &
                     20*lam(9) - 18*lam(10)) +                                  &
                  lam(8)*(-3*lam(4) + 6*lam(5) + 3*lam(6) - 6*lam(7) +          &
                     5*lam(9) + 15*lam(10))) +                                  &
               qzt2*(1 + Rqrd)*(1 + Rqsd)*                                      &
                (9*lam(4)**2 + 117*lam(5)**2 + 90*lam(5)*lam(6) +               &
                  9*lam(6)**2 +                                                 &
                  36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -               &
                  234*lam(5)*lam(7) - 90*lam(6)*lam(7) + 117*lam(7)**2 +        &
                  144*lam(1)*lam(8) +                                           &
                  24*(24*lam(2)**2 - 48*lam(2)*lam(8) +                         &
                     lam(8)*(7*lam(3) + 6*lam(8))) + 240*lam(5)*lam(9) +        &
                  120*lam(6)*lam(9) - 240*lam(7)*lam(9) + 100*lam(9)**2 +       &
                  12*(39*lam(5) + 9*lam(6) - 39*lam(7) + 50*lam(9))*            &
                   lam(10) + 312*lam(10)**2 -                                   &
                  6*lam(4)*(15*lam(5) + 3*lam(6) - 15*lam(7) +                  &
                     20*lam(9) + 18*lam(10)))))) +                              &
      expdISr*tanhdr*(-2*expuISs*qxy2*(1 + Rqrd)*(1 + Rqsu)*tanhus**2*          &
          (3*dRqsu*qdB*qzt2*(1 + Rqsu)*                                         &
             (48*lam(2)**2 +                                                    &
               3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                 &
                  lam(6)**2 +                                                   &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                    &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               144*lam(2)*lam(8) -                                              &
               12*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) +              &
               80*lam(8)**2) +                                                  &
            dRqrd*quB*qxy2*(1 + Rqrd)*                                          &
             (576*lam(2)**2 + 576*lam(2)*lam(8) +                               &
               144*lam(8)*(2*lam(1) + lam(8)) +                                 &
               24*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +               &
                  14*lam(8)) +                                                  &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -         &
                  12*lam(10))**2)) +                                            &
         expdISs*qdB*quB*qzt2*tanhds*                                           &
          (9*dRqsd*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                   &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) -          &
            dRqrd*(576*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*               &
                lam(2)**2 - 288*m2*(1 + Rqrd)*sr*lam(2)*lam(4) -                &
               9*lam(4)**2 - 9*Rqsd*lam(4)**2 + 18*qzt2*sr*lam(4)**2 +          &
               36*qzt2*Rqrd*sr*lam(4)**2 +                                      &
               18*qzt2*Rqrd**2*sr*lam(4)**2 +                                   &
               18*qzt2*Rqsd*sr*lam(4)**2 +                                      &
               36*qzt2*Rqrd*Rqsd*sr*lam(4)**2 +                                 &
               18*qzt2*Rqrd**2*Rqsd*sr*lam(4)**2 +                              &
               1440*m2*(1 + Rqrd)*sr*lam(2)*lam(5) + 90*lam(4)*lam(5) +         &
               90*Rqsd*lam(4)*lam(5) - 180*qzt2*sr*lam(4)*lam(5) -              &
               360*qzt2*Rqrd*sr*lam(4)*lam(5) -                                 &
               180*qzt2*Rqrd**2*sr*lam(4)*lam(5) -                              &
               180*qzt2*Rqsd*sr*lam(4)*lam(5) -                                 &
               360*qzt2*Rqrd*Rqsd*sr*lam(4)*lam(5) -                            &
               180*qzt2*Rqrd**2*Rqsd*sr*lam(4)*lam(5) - 117*lam(5)**2 -         &
               117*Rqsd*lam(5)**2 + 234*qzt2*sr*lam(5)**2 +                     &
               468*qzt2*Rqrd*sr*lam(5)**2 +                                     &
               234*qzt2*Rqrd**2*sr*lam(5)**2 +                                  &
               234*qzt2*Rqsd*sr*lam(5)**2 +                                     &
               468*qzt2*Rqrd*Rqsd*sr*lam(5)**2 +                                &
               234*qzt2*Rqrd**2*Rqsd*sr*lam(5)**2 +                             &
               288*m2*(1 + Rqrd)*sr*lam(2)*lam(6) + 18*lam(4)*lam(6) +          &
               18*Rqsd*lam(4)*lam(6) - 36*qzt2*sr*lam(4)*lam(6) -               &
               72*qzt2*Rqrd*sr*lam(4)*lam(6) -                                  &
               36*qzt2*Rqrd**2*sr*lam(4)*lam(6) -                               &
               36*qzt2*Rqsd*sr*lam(4)*lam(6) -                                  &
               72*qzt2*Rqrd*Rqsd*sr*lam(4)*lam(6) -                             &
               36*qzt2*Rqrd**2*Rqsd*sr*lam(4)*lam(6) -                          &
               90*lam(5)*lam(6) - 90*Rqsd*lam(5)*lam(6) +                       &
               180*qzt2*sr*lam(5)*lam(6) +                                      &
               360*qzt2*Rqrd*sr*lam(5)*lam(6) +                                 &
               180*qzt2*Rqrd**2*sr*lam(5)*lam(6) +                              &
               180*qzt2*Rqsd*sr*lam(5)*lam(6) +                                 &
               360*qzt2*Rqrd*Rqsd*sr*lam(5)*lam(6) +                            &
               180*qzt2*Rqrd**2*Rqsd*sr*lam(5)*lam(6) - 9*lam(6)**2 -           &
               9*Rqsd*lam(6)**2 + 18*qzt2*sr*lam(6)**2 +                        &
               36*qzt2*Rqrd*sr*lam(6)**2 +                                      &
               18*qzt2*Rqrd**2*sr*lam(6)**2 +                                   &
               18*qzt2*Rqsd*sr*lam(6)**2 +                                      &
               36*qzt2*Rqrd*Rqsd*sr*lam(6)**2 +                                 &
               18*qzt2*Rqrd**2*Rqsd*sr*lam(6)**2 -                              &
               1440*m2*(1 + Rqrd)*sr*lam(2)*lam(7) - 90*lam(4)*lam(7) -         &
               90*Rqsd*lam(4)*lam(7) + 180*qzt2*sr*lam(4)*lam(7) +              &
               360*qzt2*Rqrd*sr*lam(4)*lam(7) +                                 &
               180*qzt2*Rqrd**2*sr*lam(4)*lam(7) +                              &
               180*qzt2*Rqsd*sr*lam(4)*lam(7) +                                 &
               360*qzt2*Rqrd*Rqsd*sr*lam(4)*lam(7) +                            &
               180*qzt2*Rqrd**2*Rqsd*sr*lam(4)*lam(7) +                         &
               234*lam(5)*lam(7) + 234*Rqsd*lam(5)*lam(7) -                     &
               468*qzt2*sr*lam(5)*lam(7) -                                      &
               936*qzt2*Rqrd*sr*lam(5)*lam(7) -                                 &
               468*qzt2*Rqrd**2*sr*lam(5)*lam(7) -                              &
               468*qzt2*Rqsd*sr*lam(5)*lam(7) -                                 &
               936*qzt2*Rqrd*Rqsd*sr*lam(5)*lam(7) -                            &
               468*qzt2*Rqrd**2*Rqsd*sr*lam(5)*lam(7) +                         &
               90*lam(6)*lam(7) + 90*Rqsd*lam(6)*lam(7) -                       &
               180*qzt2*sr*lam(6)*lam(7) -                                      &
               360*qzt2*Rqrd*sr*lam(6)*lam(7) -                                 &
               180*qzt2*Rqrd**2*sr*lam(6)*lam(7) -                              &
               180*qzt2*Rqsd*sr*lam(6)*lam(7) -                                 &
               360*qzt2*Rqrd*Rqsd*sr*lam(6)*lam(7) -                            &
               180*qzt2*Rqrd**2*Rqsd*sr*lam(6)*lam(7) - 117*lam(7)**2 -         &
               117*Rqsd*lam(7)**2 + 234*qzt2*sr*lam(7)**2 +                     &
               468*qzt2*Rqrd*sr*lam(7)**2 +                                     &
               234*qzt2*Rqrd**2*sr*lam(7)**2 +                                  &
               234*qzt2*Rqsd*sr*lam(7)**2 +                                     &
               468*qzt2*Rqrd*Rqsd*sr*lam(7)**2 +                                &
               234*qzt2*Rqrd**2*Rqsd*sr*lam(7)**2 - 144*lam(1)*lam(8) -         &
               144*Rqsd*lam(1)*lam(8) + 288*qzt2*sr*lam(1)*lam(8) +             &
               576*qzt2*Rqrd*sr*lam(1)*lam(8) +                                 &
               288*qzt2*Rqrd**2*sr*lam(1)*lam(8) +                              &
               288*qzt2*Rqsd*sr*lam(1)*lam(8) +                                 &
               576*qzt2*Rqrd*Rqsd*sr*lam(1)*lam(8) +                            &
               288*qzt2*Rqrd**2*Rqsd*sr*lam(1)*lam(8) -                         &
               1152*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(2)*           &
                lam(8) + 288*m2*sr*lam(4)*lam(8) +                              &
               288*m2*Rqrd*sr*lam(4)*lam(8) - 576*m2*sr*lam(5)*lam(8) -         &
               576*m2*Rqrd*sr*lam(5)*lam(8) - 288*m2*sr*lam(6)*lam(8) -         &
               288*m2*Rqrd*sr*lam(6)*lam(8) + 576*m2*sr*lam(7)*lam(8) +         &
               576*m2*Rqrd*sr*lam(7)*lam(8) - 144*lam(8)**2 -                   &
               144*Rqsd*lam(8)**2 + 288*qzt2*sr*lam(8)**2 +                     &
               576*qzt2*Rqrd*sr*lam(8)**2 +                                     &
               288*qzt2*Rqrd**2*sr*lam(8)**2 +                                  &
               288*qzt2*Rqsd*sr*lam(8)**2 +                                     &
               576*qzt2*Rqrd*Rqsd*sr*lam(8)**2 +                                &
               288*qzt2*Rqrd**2*Rqsd*sr*lam(8)**2 +                             &
               12*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(3)*             &
                (3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) + 14*lam(8)) +         &
               1920*m2*(1 + Rqrd)*sr*lam(2)*lam(9) + 120*lam(4)*lam(9) +        &
               120*Rqsd*lam(4)*lam(9) - 240*qzt2*sr*lam(4)*lam(9) -             &
               480*qzt2*Rqrd*sr*lam(4)*lam(9) -                                 &
               240*qzt2*Rqrd**2*sr*lam(4)*lam(9) -                              &
               240*qzt2*Rqsd*sr*lam(4)*lam(9) -                                 &
               480*qzt2*Rqrd*Rqsd*sr*lam(4)*lam(9) -                            &
               240*qzt2*Rqrd**2*Rqsd*sr*lam(4)*lam(9) -                         &
               240*lam(5)*lam(9) - 240*Rqsd*lam(5)*lam(9) +                     &
               480*qzt2*sr*lam(5)*lam(9) +                                      &
               960*qzt2*Rqrd*sr*lam(5)*lam(9) +                                 &
               480*qzt2*Rqrd**2*sr*lam(5)*lam(9) +                              &
               480*qzt2*Rqsd*sr*lam(5)*lam(9) +                                 &
               960*qzt2*Rqrd*Rqsd*sr*lam(5)*lam(9) +                            &
               480*qzt2*Rqrd**2*Rqsd*sr*lam(5)*lam(9) -                         &
               120*lam(6)*lam(9) - 120*Rqsd*lam(6)*lam(9) +                     &
               240*qzt2*sr*lam(6)*lam(9) +                                      &
               480*qzt2*Rqrd*sr*lam(6)*lam(9) +                                 &
               240*qzt2*Rqrd**2*sr*lam(6)*lam(9) +                              &
               240*qzt2*Rqsd*sr*lam(6)*lam(9) +                                 &
               480*qzt2*Rqrd*Rqsd*sr*lam(6)*lam(9) +                            &
               240*qzt2*Rqrd**2*Rqsd*sr*lam(6)*lam(9) +                         &
               240*lam(7)*lam(9) + 240*Rqsd*lam(7)*lam(9) -                     &
               480*qzt2*sr*lam(7)*lam(9) -                                      &
               960*qzt2*Rqrd*sr*lam(7)*lam(9) -                                 &
               480*qzt2*Rqrd**2*sr*lam(7)*lam(9) -                              &
               480*qzt2*Rqsd*sr*lam(7)*lam(9) -                                 &
               960*qzt2*Rqrd*Rqsd*sr*lam(7)*lam(9) -                            &
               480*qzt2*Rqrd**2*Rqsd*sr*lam(7)*lam(9) -                         &
               480*m2*sr*lam(8)*lam(9) - 480*m2*Rqrd*sr*lam(8)*lam(9) -         &
               100*lam(9)**2 - 100*Rqsd*lam(9)**2 +                             &
               200*qzt2*sr*lam(9)**2 + 400*qzt2*Rqrd*sr*lam(9)**2 +             &
               200*qzt2*Rqrd**2*sr*lam(9)**2 +                                  &
               200*qzt2*Rqsd*sr*lam(9)**2 +                                     &
               400*qzt2*Rqrd*Rqsd*sr*lam(9)**2 +                                &
               200*qzt2*Rqrd**2*Rqsd*sr*lam(9)**2 +                             &
               12*(24*m2*(1 + Rqrd)*sr*(6*lam(2) - 5*lam(8)) -                  &
                  (1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                    &
                   (9*lam(4) - 39*lam(5) - 9*lam(6) + 39*lam(7) -               &
                     50*lam(9)))*lam(10) +                                      &
               312*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(10)**2))       &
+ expuISs*qdB*quB*qzt2*tanhus*                                                  &
          (-3*dRqsu*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                  &
             (48*lam(2)**2 +                                                    &
               3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                 &
                  lam(6)**2 +                                                   &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                    &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               144*lam(2)*lam(8) -                                              &
               12*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) +              &
               80*lam(8)**2) +                                                  &
            dRqrd*(576*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*               &
                lam(2)**2 - 288*m2*(1 + Rqrd)*sr*lam(2)*lam(4) -                &
               9*lam(4)**2 - 9*Rqsu*lam(4)**2 + 18*qzt2*sr*lam(4)**2 +          &
               36*qzt2*Rqrd*sr*lam(4)**2 +                                      &
               18*qzt2*Rqrd**2*sr*lam(4)**2 +                                   &
               18*qzt2*Rqsu*sr*lam(4)**2 +                                      &
               36*qzt2*Rqrd*Rqsu*sr*lam(4)**2 +                                 &
               18*qzt2*Rqrd**2*Rqsu*sr*lam(4)**2 -                              &
               288*m2*(1 + Rqrd)*sr*lam(2)*lam(5) - 18*lam(4)*lam(5) -          &
               18*Rqsu*lam(4)*lam(5) + 36*qzt2*sr*lam(4)*lam(5) +               &
               72*qzt2*Rqrd*sr*lam(4)*lam(5) +                                  &
               36*qzt2*Rqrd**2*sr*lam(4)*lam(5) +                               &
               36*qzt2*Rqsu*sr*lam(4)*lam(5) +                                  &
               72*qzt2*Rqrd*Rqsu*sr*lam(4)*lam(5) +                             &
               36*qzt2*Rqrd**2*Rqsu*sr*lam(4)*lam(5) - 9*lam(5)**2 -            &
               9*Rqsu*lam(5)**2 + 18*qzt2*sr*lam(5)**2 +                        &
               36*qzt2*Rqrd*sr*lam(5)**2 +                                      &
               18*qzt2*Rqrd**2*sr*lam(5)**2 +                                   &
               18*qzt2*Rqsu*sr*lam(5)**2 +                                      &
               36*qzt2*Rqrd*Rqsu*sr*lam(5)**2 +                                 &
               18*qzt2*Rqrd**2*Rqsu*sr*lam(5)**2 +                              &
               288*m2*(1 + Rqrd)*sr*lam(2)*lam(6) + 18*lam(4)*lam(6) +          &
               18*Rqsu*lam(4)*lam(6) - 36*qzt2*sr*lam(4)*lam(6) -               &
               72*qzt2*Rqrd*sr*lam(4)*lam(6) -                                  &
               36*qzt2*Rqrd**2*sr*lam(4)*lam(6) -                               &
               36*qzt2*Rqsu*sr*lam(4)*lam(6) -                                  &
               72*qzt2*Rqrd*Rqsu*sr*lam(4)*lam(6) -                             &
               36*qzt2*Rqrd**2*Rqsu*sr*lam(4)*lam(6) +                          &
               18*lam(5)*lam(6) + 18*Rqsu*lam(5)*lam(6) -                       &
               36*qzt2*sr*lam(5)*lam(6) -                                       &
               72*qzt2*Rqrd*sr*lam(5)*lam(6) -                                  &
               36*qzt2*Rqrd**2*sr*lam(5)*lam(6) -                               &
               36*qzt2*Rqsu*sr*lam(5)*lam(6) -                                  &
               72*qzt2*Rqrd*Rqsu*sr*lam(5)*lam(6) -                             &
               36*qzt2*Rqrd**2*Rqsu*sr*lam(5)*lam(6) - 9*lam(6)**2 -            &
               9*Rqsu*lam(6)**2 + 18*qzt2*sr*lam(6)**2 +                        &
               36*qzt2*Rqrd*sr*lam(6)**2 +                                      &
               18*qzt2*Rqrd**2*sr*lam(6)**2 +                                   &
               18*qzt2*Rqsu*sr*lam(6)**2 +                                      &
               36*qzt2*Rqrd*Rqsu*sr*lam(6)**2 +                                 &
               18*qzt2*Rqrd**2*Rqsu*sr*lam(6)**2 +                              &
               288*m2*(1 + Rqrd)*sr*lam(2)*lam(7) + 18*lam(4)*lam(7) +          &
               18*Rqsu*lam(4)*lam(7) - 36*qzt2*sr*lam(4)*lam(7) -               &
               72*qzt2*Rqrd*sr*lam(4)*lam(7) -                                  &
               36*qzt2*Rqrd**2*sr*lam(4)*lam(7) -                               &
               36*qzt2*Rqsu*sr*lam(4)*lam(7) -                                  &
               72*qzt2*Rqrd*Rqsu*sr*lam(4)*lam(7) -                             &
               36*qzt2*Rqrd**2*Rqsu*sr*lam(4)*lam(7) +                          &
               18*lam(5)*lam(7) + 18*Rqsu*lam(5)*lam(7) -                       &
               36*qzt2*sr*lam(5)*lam(7) -                                       &
               72*qzt2*Rqrd*sr*lam(5)*lam(7) -                                  &
               36*qzt2*Rqrd**2*sr*lam(5)*lam(7) -                               &
               36*qzt2*Rqsu*sr*lam(5)*lam(7) -                                  &
               72*qzt2*Rqrd*Rqsu*sr*lam(5)*lam(7) -                             &
               36*qzt2*Rqrd**2*Rqsu*sr*lam(5)*lam(7) -                          &
               18*lam(6)*lam(7) - 18*Rqsu*lam(6)*lam(7) +                       &
               36*qzt2*sr*lam(6)*lam(7) +                                       &
               72*qzt2*Rqrd*sr*lam(6)*lam(7) +                                  &
               36*qzt2*Rqrd**2*sr*lam(6)*lam(7) +                               &
               36*qzt2*Rqsu*sr*lam(6)*lam(7) +                                  &
               72*qzt2*Rqrd*Rqsu*sr*lam(6)*lam(7) +                             &
               36*qzt2*Rqrd**2*Rqsu*sr*lam(6)*lam(7) - 9*lam(7)**2 -            &
               9*Rqsu*lam(7)**2 + 18*qzt2*sr*lam(7)**2 +                        &
               36*qzt2*Rqrd*sr*lam(7)**2 +                                      &
               18*qzt2*Rqrd**2*sr*lam(7)**2 +                                   &
               18*qzt2*Rqsu*sr*lam(7)**2 +                                      &
               36*qzt2*Rqrd*Rqsu*sr*lam(7)**2 +                                 &
               18*qzt2*Rqrd**2*Rqsu*sr*lam(7)**2 + 288*lam(1)*lam(8) +          &
               288*Rqsu*lam(1)*lam(8) - 576*qzt2*sr*lam(1)*lam(8) -             &
               1152*qzt2*Rqrd*sr*lam(1)*lam(8) -                                &
               576*qzt2*Rqrd**2*sr*lam(1)*lam(8) -                              &
               576*qzt2*Rqsu*sr*lam(1)*lam(8) -                                 &
               1152*qzt2*Rqrd*Rqsu*sr*lam(1)*lam(8) -                           &
               576*qzt2*Rqrd**2*Rqsu*sr*lam(1)*lam(8) +                         &
               576*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(2)*            &
                lam(8) - 144*m2*sr*lam(4)*lam(8) -                              &
               144*m2*Rqrd*sr*lam(4)*lam(8) - 144*m2*sr*lam(5)*lam(8) -         &
               144*m2*Rqrd*sr*lam(5)*lam(8) + 144*m2*sr*lam(6)*lam(8) +         &
               144*m2*Rqrd*sr*lam(6)*lam(8) + 144*m2*sr*lam(7)*lam(8) +         &
               144*m2*Rqrd*sr*lam(7)*lam(8) - 144*lam(8)**2 -                   &
               144*Rqsu*lam(8)**2 + 288*qzt2*sr*lam(8)**2 +                     &
               576*qzt2*Rqrd*sr*lam(8)**2 +                                     &
               288*qzt2*Rqrd**2*sr*lam(8)**2 +                                  &
               288*qzt2*Rqsu*sr*lam(8)**2 +                                     &
               576*qzt2*Rqrd*Rqsu*sr*lam(8)**2 +                                &
               288*qzt2*Rqrd**2*Rqsu*sr*lam(8)**2 -                             &
               24*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(3)*             &
                (3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) + 14*lam(8)) -         &
               960*m2*(1 + Rqrd)*sr*lam(2)*lam(9) - 60*lam(4)*lam(9) -          &
               60*Rqsu*lam(4)*lam(9) + 120*qzt2*sr*lam(4)*lam(9) +              &
               240*qzt2*Rqrd*sr*lam(4)*lam(9) +                                 &
               120*qzt2*Rqrd**2*sr*lam(4)*lam(9) +                              &
               120*qzt2*Rqsu*sr*lam(4)*lam(9) +                                 &
               240*qzt2*Rqrd*Rqsu*sr*lam(4)*lam(9) +                            &
               120*qzt2*Rqrd**2*Rqsu*sr*lam(4)*lam(9) -                         &
               60*lam(5)*lam(9) - 60*Rqsu*lam(5)*lam(9) +                       &
               120*qzt2*sr*lam(5)*lam(9) +                                      &
               240*qzt2*Rqrd*sr*lam(5)*lam(9) +                                 &
               120*qzt2*Rqrd**2*sr*lam(5)*lam(9) +                              &
               120*qzt2*Rqsu*sr*lam(5)*lam(9) +                                 &
               240*qzt2*Rqrd*Rqsu*sr*lam(5)*lam(9) +                            &
               120*qzt2*Rqrd**2*Rqsu*sr*lam(5)*lam(9) +                         &
               60*lam(6)*lam(9) + 60*Rqsu*lam(6)*lam(9) -                       &
               120*qzt2*sr*lam(6)*lam(9) -                                      &
               240*qzt2*Rqrd*sr*lam(6)*lam(9) -                                 &
               120*qzt2*Rqrd**2*sr*lam(6)*lam(9) -                              &
               120*qzt2*Rqsu*sr*lam(6)*lam(9) -                                 &
               240*qzt2*Rqrd*Rqsu*sr*lam(6)*lam(9) -                            &
               120*qzt2*Rqrd**2*Rqsu*sr*lam(6)*lam(9) +                         &
               60*lam(7)*lam(9) + 60*Rqsu*lam(7)*lam(9) -                       &
               120*qzt2*sr*lam(7)*lam(9) -                                      &
               240*qzt2*Rqrd*sr*lam(7)*lam(9) -                                 &
               120*qzt2*Rqrd**2*sr*lam(7)*lam(9) -                              &
               120*qzt2*Rqsu*sr*lam(7)*lam(9) -                                 &
               240*qzt2*Rqrd*Rqsu*sr*lam(7)*lam(9) -                            &
               120*qzt2*Rqrd**2*Rqsu*sr*lam(7)*lam(9) -                         &
               480*m2*sr*lam(8)*lam(9) - 480*m2*Rqrd*sr*lam(8)*lam(9) -         &
               100*lam(9)**2 - 100*Rqsu*lam(9)**2 +                             &
               200*qzt2*sr*lam(9)**2 + 400*qzt2*Rqrd*sr*lam(9)**2 +             &
               200*qzt2*Rqrd**2*sr*lam(9)**2 +                                  &
               200*qzt2*Rqsu*sr*lam(9)**2 +                                     &
               400*qzt2*Rqrd*Rqsu*sr*lam(9)**2 +                                &
               200*qzt2*Rqrd**2*Rqsu*sr*lam(9)**2 +                             &
               24*(24*m2*(1 + Rqrd)*sr*(2*lam(2) + lam(8)) -                    &
                  (1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                    &
                   (3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                 &
                     10*lam(9)))*lam(10) +                                      &
               144*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(10)**2))       &
+ 2*expdISs*quB*qxy2*(1 + Rqrd)*(1 + Rqsd)*tanhds**2*                           &
          (9*dRqsd*qzt2*(1 + Rqsd)*                                             &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            dRqrd*qxy2*(1 + Rqrd)*                                              &
             (9*lam(4)**2 + 117*lam(5)**2 + 90*lam(5)*lam(6) +                  &
               9*lam(6)**2 - 36*lam(3)*                                         &
                (lam(4) - lam(5) + lam(6) - lam(7)) -                           &
               234*lam(5)*lam(7) - 90*lam(6)*lam(7) + 117*lam(7)**2 -           &
               144*lam(1)*lam(8) +                                              &
               24*(24*lam(2)**2 - 48*lam(2)*lam(8) +                            &
                  lam(8)*(-7*lam(3) + 6*lam(8))) + 240*lam(5)*lam(9) +          &
               120*lam(6)*lam(9) - 240*lam(7)*lam(9) + 100*lam(9)**2 +          &
               12*(39*lam(5) + 9*lam(6) - 39*lam(7) + 50*lam(9))*               &
                lam(10) + 312*lam(10)**2 -                                      &
               6*lam(4)*(15*lam(5) + 3*lam(6) - 15*lam(7) + 20*lam(9) +         &
                  18*lam(10)))) -                                               &
         2*dRqrd*quB*qxy2*(1 + Rqrd)*                                           &
          (-(expuISs*(-24*m2*                                                   &
                  (24*lam(1)*lam(3) + 14*lam(3)**2 +                            &
                    lam(8)*(9*lam(4) - 3*lam(5) + 3*lam(6) -                    &
                       9*lam(7) + 14*lam(8) + 10*lam(9)) +                      &
                    lam(2)*(6*lam(4) + 6*lam(5) - 6*(lam(6) + lam(7)) +         &
                       20*lam(9))) +                                            &
                 qxy2*(1 + Rqrd)*(1 + Rqsu)*                                    &
                  (576*lam(2)**2 + 576*lam(2)*lam(8) +                          &
                    144*lam(8)*(2*lam(1) + lam(8)) +                            &
                    24*lam(3)*                                                  &
                     (3*(lam(4) - lam(5) + lam(6) - lam(7)) +                   &
                       14*lam(8)) +                                             &
                    (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                &
                       10*lam(9) - 12*lam(10))**2) +                            &
                 qzt2*(1 + Rqrd)*(1 + Rqsu)*                                    &
                  (576*lam(2)**2 + 576*lam(2)*lam(8) +                          &
                    144*lam(8)*(2*lam(1) + lam(8)) +                            &
                    24*lam(3)*                                                  &
                     (3*(lam(4) - lam(5) + lam(6) - lam(7)) +                   &
                       14*lam(8)) +                                             &
                    (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                &
                       10*lam(9) - 12*lam(10))**2) +                            &
                 288*m2*(2*lam(2) + lam(8))*lam(10))) +                         &
            expdISs*(qxy2*(1 + Rqrd)*(1 + Rqsd)*                                &
                (9*lam(4)**2 + 117*lam(5)**2 + 90*lam(5)*lam(6) +               &
                  9*lam(6)**2 -                                                 &
                  36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -               &
                  234*lam(5)*lam(7) - 90*lam(6)*lam(7) + 117*lam(7)**2 -        &
                  144*lam(1)*lam(8) +                                           &
                  24*(24*lam(2)**2 - 48*lam(2)*lam(8) +                         &
                     lam(8)*(-7*lam(3) + 6*lam(8))) +                           &
                  240*lam(5)*lam(9) + 120*lam(6)*lam(9) -                       &
                  240*lam(7)*lam(9) + 100*lam(9)**2 +                           &
                  12*(39*lam(5) + 9*lam(6) - 39*lam(7) + 50*lam(9))*            &
                   lam(10) + 312*lam(10)**2 -                                   &
                  6*lam(4)*(15*lam(5) + 3*lam(6) - 15*lam(7) +                  &
                     20*lam(9) + 18*lam(10))) +                                 &
               qzt2*(1 + Rqrd)*(1 + Rqsd)*                                      &
                (9*lam(4)**2 + 117*lam(5)**2 + 90*lam(5)*lam(6) +               &
                  9*lam(6)**2 -                                                 &
                  36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -               &
                  234*lam(5)*lam(7) - 90*lam(6)*lam(7) + 117*lam(7)**2 -        &
                  144*lam(1)*lam(8) +                                           &
                  24*(24*lam(2)**2 - 48*lam(2)*lam(8) +                         &
                     lam(8)*(-7*lam(3) + 6*lam(8))) +                           &
                  240*lam(5)*lam(9) + 120*lam(6)*lam(9) -                       &
                  240*lam(7)*lam(9) + 100*lam(9)**2 +                           &
                  12*(39*lam(5) + 9*lam(6) - 39*lam(7) + 50*lam(9))*            &
                   lam(10) + 312*lam(10)**2 -                                   &
                  6*lam(4)*(15*lam(5) + 3*lam(6) - 15*lam(7) +                  &
                     20*lam(9) + 18*lam(10))) +                                 &
               24*m2*(12*lam(1)*lam(3) + 7*lam(3)**2 +                          &
                  lam(2)*(-6*lam(4) + 30*lam(5) + 6*lam(6) - 30*lam(7) +        &
                     40*lam(9) + 36*lam(10)) +                                  &
                  lam(8)*(9*lam(4) - 15*lam(5) - 3*lam(6) + 9*lam(7) +          &
                     7*lam(8) - 10*(lam(9) + 3*lam(10))))))) +                  &
      expdISr*expdISs*quB*(-18*dRqsd*qxy2**2*(1 + Rqrd)*(1 + Rqsd)**2*          &
          tanhds**3*(16*lam(2)**2 +                                             &
            (lam(4) + lam(5) + lam(6) + lam(7))**2) +                           &
         18*dRqsd*qxy2*(1 + Rqsd)*tanhds*                                       &
          (16*m2*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                   &
            qxy2*(1 + Rqrd)*(1 + Rqsd)*                                         &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            qzt2*(1 + Rqrd)*(1 + Rqsd)*                                         &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2)) +         &
         qdB*qxy2*tanhds**2*(-9*dRqsd*(1 + Rqrd)*                               &
             (-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                                    &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            dRqrd*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                    &
             (9*lam(4)**2 + 117*lam(5)**2 + 90*lam(5)*lam(6) +                  &
               9*lam(6)**2 - 36*lam(3)*                                         &
                (lam(4) - lam(5) + lam(6) - lam(7)) -                           &
               234*lam(5)*lam(7) - 90*lam(6)*lam(7) + 117*lam(7)**2 -           &
               144*lam(1)*lam(8) +                                              &
               24*(24*lam(2)**2 - 48*lam(2)*lam(8) +                            &
                  lam(8)*(-7*lam(3) + 6*lam(8))) + 240*lam(5)*lam(9) +          &
               120*lam(6)*lam(9) - 240*lam(7)*lam(9) + 100*lam(9)**2 +          &
               12*(39*lam(5) + 9*lam(6) - 39*lam(7) + 50*lam(9))*               &
                lam(10) + 312*lam(10)**2 -                                      &
               6*lam(4)*(15*lam(5) + 3*lam(6) - 15*lam(7) + 20*lam(9) +         &
                  18*lam(10)))) +                                               &
         qdB*(9*dRqsd*(qxy2*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*          &
                (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2)         &
+ qzt2*(16*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*lam(2)**2 +                &
                  32*m2*(1 + Rqsd)*ss*lam(2)*                                   &
                   (lam(4) + lam(5) + lam(6) + lam(7)) +                        &
                  (1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                    &
                   (lam(4) + lam(5) + lam(6) + lam(7))**2)) -                   &
            dRqrd*(qxy2*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*              &
                (9*lam(4)**2 + 117*lam(5)**2 + 90*lam(5)*lam(6) +               &
                  9*lam(6)**2 -                                                 &
                  36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -               &
                  234*lam(5)*lam(7) - 90*lam(6)*lam(7) + 117*lam(7)**2 -        &
                  144*lam(1)*lam(8) +                                           &
                  24*(24*lam(2)**2 - 48*lam(2)*lam(8) +                         &
                     lam(8)*(-7*lam(3) + 6*lam(8))) +                           &
                  240*lam(5)*lam(9) + 120*lam(6)*lam(9) -                       &
                  240*lam(7)*lam(9) + 100*lam(9)**2 +                           &
                  12*(39*lam(5) + 9*lam(6) - 39*lam(7) + 50*lam(9))*            &
                   lam(10) + 312*lam(10)**2 -                                   &
                  6*lam(4)*(15*lam(5) + 3*lam(6) - 15*lam(7) +                  &
                     20*lam(9) + 18*lam(10))) +                                 &
               qzt2*(576*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*             &
                   lam(2)**2 -                                                  &
                  96*lam(2)*(12*(1 + Rqsd)*                                     &
                      (-1 + 2*qzt2*(1 + Rqrd)**2*sr)*lam(8) +                   &
                     m2*(1 + Rqrd)*sr*                                          &
                      (3*lam(4) - 15*lam(5) - 3*lam(6) + 15*lam(7) -            &
                        20*lam(9) - 18*lam(10))) -                              &
                  (1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                    &
                   (-9*lam(4)**2 - 117*lam(5)**2 - 90*lam(5)*lam(6) -           &
                     9*lam(6)**2 + 234*lam(5)*lam(7) +                          &
                     90*lam(6)*lam(7) - 117*lam(7)**2 +                         &
                     144*(lam(1) - lam(8))*lam(8) +                             &
                     12*lam(3)*                                                 &
                      (3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +                &
                        14*lam(8)) - 240*lam(5)*lam(9) -                        &
                     120*lam(6)*lam(9) + 240*lam(7)*lam(9) -                    &
                     100*lam(9)**2 -                                            &
                     12*(39*lam(5) + 9*lam(6) - 39*lam(7) + 50*lam(9))*         &
                      lam(10) - 312*lam(10)**2 +                                &
                     6*lam(4)*                                                  &
                      (15*lam(5) + 3*lam(6) - 15*lam(7) + 20*lam(9) +           &
                        18*lam(10))) +                                          &
                  48*m2*(1 + Rqrd)*sr*                                          &
                   (12*lam(1)*lam(3) + 7*lam(3)**2 +                            &
                     lam(8)*(9*lam(4) - 15*lam(5) - 3*lam(6) + 9*lam(7) +       &
                        7*lam(8) - 10*(lam(9) + 3*lam(10)))))))) +              &
      expuISr*expuISs*qdB*(-18*dRqsu*qxy2**2*(1 + Rqru)*(1 + Rqsu)**2*          &
          tanhus**3*(16*lam(2)**2 +                                             &
            (lam(4) + lam(5) + lam(6) + lam(7))**2) +                           &
         18*dRqsu*qxy2*(1 + Rqsu)*tanhus*                                       &
          (16*m2*lam(2)*(lam(4) + lam(5) + lam(6) + lam(7)) +                   &
            qxy2*(1 + Rqru)*(1 + Rqsu)*                                         &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            qzt2*(1 + Rqru)*(1 + Rqsu)*                                         &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2)) +         &
         quB*qxy2*tanhus**2*(-9*dRqsu*(1 + Rqru)*                               &
             (-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                                    &
             (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +          &
            dRqru*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                    &
             (9*lam(4)**2 + 117*lam(5)**2 + 90*lam(5)*lam(6) +                  &
               9*lam(6)**2 - 36*lam(3)*                                         &
                (lam(4) - lam(5) + lam(6) - lam(7)) -                           &
               234*lam(5)*lam(7) - 90*lam(6)*lam(7) + 117*lam(7)**2 -           &
               144*lam(1)*lam(8) +                                              &
               24*(24*lam(2)**2 - 48*lam(2)*lam(8) +                            &
                  lam(8)*(-7*lam(3) + 6*lam(8))) + 240*lam(5)*lam(9) +          &
               120*lam(6)*lam(9) - 240*lam(7)*lam(9) + 100*lam(9)**2 +          &
               12*(39*lam(5) + 9*lam(6) - 39*lam(7) + 50*lam(9))*               &
                lam(10) + 312*lam(10)**2 -                                      &
               6*lam(4)*(15*lam(5) + 3*lam(6) - 15*lam(7) + 20*lam(9) +         &
                  18*lam(10)))) +                                               &
         quB*(9*dRqsu*(qxy2*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*          &
                (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2)         &
+ qzt2*(16*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*lam(2)**2 +                &
                  32*m2*(1 + Rqsu)*ss*lam(2)*                                   &
                   (lam(4) + lam(5) + lam(6) + lam(7)) +                        &
                  (1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                    &
                   (lam(4) + lam(5) + lam(6) + lam(7))**2)) -                   &
            dRqru*(qxy2*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*              &
                (9*lam(4)**2 + 117*lam(5)**2 + 90*lam(5)*lam(6) +               &
                  9*lam(6)**2 -                                                 &
                  36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -               &
                  234*lam(5)*lam(7) - 90*lam(6)*lam(7) + 117*lam(7)**2 -        &
                  144*lam(1)*lam(8) +                                           &
                  24*(24*lam(2)**2 - 48*lam(2)*lam(8) +                         &
                     lam(8)*(-7*lam(3) + 6*lam(8))) +                           &
                  240*lam(5)*lam(9) + 120*lam(6)*lam(9) -                       &
                  240*lam(7)*lam(9) + 100*lam(9)**2 +                           &
                  12*(39*lam(5) + 9*lam(6) - 39*lam(7) + 50*lam(9))*            &
                   lam(10) + 312*lam(10)**2 -                                   &
                  6*lam(4)*(15*lam(5) + 3*lam(6) - 15*lam(7) +                  &
                     20*lam(9) + 18*lam(10))) +                                 &
               qzt2*(576*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*             &
                   lam(2)**2 -                                                  &
                  96*lam(2)*(12*(1 + Rqsu)*                                     &
                      (-1 + 2*qzt2*(1 + Rqru)**2*sr)*lam(8) +                   &
                     m2*(1 + Rqru)*sr*                                          &
                      (3*lam(4) - 15*lam(5) - 3*lam(6) + 15*lam(7) -            &
                        20*lam(9) - 18*lam(10))) -                              &
                  (1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                    &
                   (-9*lam(4)**2 - 117*lam(5)**2 - 90*lam(5)*lam(6) -           &
                     9*lam(6)**2 + 234*lam(5)*lam(7) +                          &
                     90*lam(6)*lam(7) - 117*lam(7)**2 +                         &
                     144*(lam(1) - lam(8))*lam(8) +                             &
                     12*lam(3)*                                                 &
                      (3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) +                &
                        14*lam(8)) - 240*lam(5)*lam(9) -                        &
                     120*lam(6)*lam(9) + 240*lam(7)*lam(9) -                    &
                     100*lam(9)**2 -                                            &
                     12*(39*lam(5) + 9*lam(6) - 39*lam(7) + 50*lam(9))*         &
                      lam(10) - 312*lam(10)**2 +                                &
                     6*lam(4)*                                                  &
                      (15*lam(5) + 3*lam(6) - 15*lam(7) + 20*lam(9) +           &
                        18*lam(10))) +                                          &
                  48*m2*(1 + Rqru)*sr*                                          &
                   (12*lam(1)*lam(3) + 7*lam(3)**2 +                            &
                     lam(8)*(9*lam(4) - 15*lam(5) - 3*lam(6) + 9*lam(7) +       &
                        7*lam(8) - 10*(lam(9) + 3*lam(10)))))))) +              &
      expuISr*tanhur*(-2*expdISs*qxy2*(1 + Rqru)*(1 + Rqsd)*tanhds**2*          &
          (3*dRqsd*quB*qzt2*(1 + Rqsd)*                                         &
             (48*lam(2)**2 +                                                    &
               3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                 &
                  lam(6)**2 +                                                   &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                    &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               144*lam(2)*lam(8) -                                              &
               12*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) +              &
               80*lam(8)**2) +                                                  &
            dRqru*qdB*qxy2*(1 + Rqru)*                                          &
             (576*lam(2)**2 + 576*lam(2)*lam(8) +                               &
               144*lam(8)*(2*lam(1) + lam(8)) +                                 &
               24*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +               &
                  14*lam(8)) +                                                  &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -         &
                  12*lam(10))**2)) +                                            &
         expdISs*qdB*quB*qzt2*tanhds*                                           &
          (-3*dRqsd*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                  &
             (48*lam(2)**2 +                                                    &
               3*(lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                 &
                  lam(6)**2 +                                                   &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                    &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               144*lam(2)*lam(8) -                                              &
               12*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) +              &
               80*lam(8)**2) +                                                  &
            dRqru*(576*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*               &
                lam(2)**2 - 288*m2*(1 + Rqru)*sr*lam(2)*lam(4) -                &
               9*lam(4)**2 - 9*Rqsd*lam(4)**2 + 18*qzt2*sr*lam(4)**2 +          &
               36*qzt2*Rqru*sr*lam(4)**2 +                                      &
               18*qzt2*Rqru**2*sr*lam(4)**2 +                                   &
               18*qzt2*Rqsd*sr*lam(4)**2 +                                      &
               36*qzt2*Rqru*Rqsd*sr*lam(4)**2 +                                 &
               18*qzt2*Rqru**2*Rqsd*sr*lam(4)**2 -                              &
               288*m2*(1 + Rqru)*sr*lam(2)*lam(5) - 18*lam(4)*lam(5) -          &
               18*Rqsd*lam(4)*lam(5) + 36*qzt2*sr*lam(4)*lam(5) +               &
               72*qzt2*Rqru*sr*lam(4)*lam(5) +                                  &
               36*qzt2*Rqru**2*sr*lam(4)*lam(5) +                               &
               36*qzt2*Rqsd*sr*lam(4)*lam(5) +                                  &
               72*qzt2*Rqru*Rqsd*sr*lam(4)*lam(5) +                             &
               36*qzt2*Rqru**2*Rqsd*sr*lam(4)*lam(5) - 9*lam(5)**2 -            &
               9*Rqsd*lam(5)**2 + 18*qzt2*sr*lam(5)**2 +                        &
               36*qzt2*Rqru*sr*lam(5)**2 +                                      &
               18*qzt2*Rqru**2*sr*lam(5)**2 +                                   &
               18*qzt2*Rqsd*sr*lam(5)**2 +                                      &
               36*qzt2*Rqru*Rqsd*sr*lam(5)**2 +                                 &
               18*qzt2*Rqru**2*Rqsd*sr*lam(5)**2 +                              &
               288*m2*(1 + Rqru)*sr*lam(2)*lam(6) + 18*lam(4)*lam(6) +          &
               18*Rqsd*lam(4)*lam(6) - 36*qzt2*sr*lam(4)*lam(6) -               &
               72*qzt2*Rqru*sr*lam(4)*lam(6) -                                  &
               36*qzt2*Rqru**2*sr*lam(4)*lam(6) -                               &
               36*qzt2*Rqsd*sr*lam(4)*lam(6) -                                  &
               72*qzt2*Rqru*Rqsd*sr*lam(4)*lam(6) -                             &
               36*qzt2*Rqru**2*Rqsd*sr*lam(4)*lam(6) +                          &
               18*lam(5)*lam(6) + 18*Rqsd*lam(5)*lam(6) -                       &
               36*qzt2*sr*lam(5)*lam(6) -                                       &
               72*qzt2*Rqru*sr*lam(5)*lam(6) -                                  &
               36*qzt2*Rqru**2*sr*lam(5)*lam(6) -                               &
               36*qzt2*Rqsd*sr*lam(5)*lam(6) -                                  &
               72*qzt2*Rqru*Rqsd*sr*lam(5)*lam(6) -                             &
               36*qzt2*Rqru**2*Rqsd*sr*lam(5)*lam(6) - 9*lam(6)**2 -            &
               9*Rqsd*lam(6)**2 + 18*qzt2*sr*lam(6)**2 +                        &
               36*qzt2*Rqru*sr*lam(6)**2 +                                      &
               18*qzt2*Rqru**2*sr*lam(6)**2 +                                   &
               18*qzt2*Rqsd*sr*lam(6)**2 +                                      &
               36*qzt2*Rqru*Rqsd*sr*lam(6)**2 +                                 &
               18*qzt2*Rqru**2*Rqsd*sr*lam(6)**2 +                              &
               288*m2*(1 + Rqru)*sr*lam(2)*lam(7) + 18*lam(4)*lam(7) +          &
               18*Rqsd*lam(4)*lam(7) - 36*qzt2*sr*lam(4)*lam(7) -               &
               72*qzt2*Rqru*sr*lam(4)*lam(7) -                                  &
               36*qzt2*Rqru**2*sr*lam(4)*lam(7) -                               &
               36*qzt2*Rqsd*sr*lam(4)*lam(7) -                                  &
               72*qzt2*Rqru*Rqsd*sr*lam(4)*lam(7) -                             &
               36*qzt2*Rqru**2*Rqsd*sr*lam(4)*lam(7) +                          &
               18*lam(5)*lam(7) + 18*Rqsd*lam(5)*lam(7) -                       &
               36*qzt2*sr*lam(5)*lam(7) -                                       &
               72*qzt2*Rqru*sr*lam(5)*lam(7) -                                  &
               36*qzt2*Rqru**2*sr*lam(5)*lam(7) -                               &
               36*qzt2*Rqsd*sr*lam(5)*lam(7) -                                  &
               72*qzt2*Rqru*Rqsd*sr*lam(5)*lam(7) -                             &
               36*qzt2*Rqru**2*Rqsd*sr*lam(5)*lam(7) -                          &
               18*lam(6)*lam(7) - 18*Rqsd*lam(6)*lam(7) +                       &
               36*qzt2*sr*lam(6)*lam(7) +                                       &
               72*qzt2*Rqru*sr*lam(6)*lam(7) +                                  &
               36*qzt2*Rqru**2*sr*lam(6)*lam(7) +                               &
               36*qzt2*Rqsd*sr*lam(6)*lam(7) +                                  &
               72*qzt2*Rqru*Rqsd*sr*lam(6)*lam(7) +                             &
               36*qzt2*Rqru**2*Rqsd*sr*lam(6)*lam(7) - 9*lam(7)**2 -            &
               9*Rqsd*lam(7)**2 + 18*qzt2*sr*lam(7)**2 +                        &
               36*qzt2*Rqru*sr*lam(7)**2 +                                      &
               18*qzt2*Rqru**2*sr*lam(7)**2 +                                   &
               18*qzt2*Rqsd*sr*lam(7)**2 +                                      &
               36*qzt2*Rqru*Rqsd*sr*lam(7)**2 +                                 &
               18*qzt2*Rqru**2*Rqsd*sr*lam(7)**2 + 288*lam(1)*lam(8) +          &
               288*Rqsd*lam(1)*lam(8) - 576*qzt2*sr*lam(1)*lam(8) -             &
               1152*qzt2*Rqru*sr*lam(1)*lam(8) -                                &
               576*qzt2*Rqru**2*sr*lam(1)*lam(8) -                              &
               576*qzt2*Rqsd*sr*lam(1)*lam(8) -                                 &
               1152*qzt2*Rqru*Rqsd*sr*lam(1)*lam(8) -                           &
               576*qzt2*Rqru**2*Rqsd*sr*lam(1)*lam(8) +                         &
               576*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*lam(2)*            &
                lam(8) - 144*m2*sr*lam(4)*lam(8) -                              &
               144*m2*Rqru*sr*lam(4)*lam(8) - 144*m2*sr*lam(5)*lam(8) -         &
               144*m2*Rqru*sr*lam(5)*lam(8) + 144*m2*sr*lam(6)*lam(8) +         &
               144*m2*Rqru*sr*lam(6)*lam(8) + 144*m2*sr*lam(7)*lam(8) +         &
               144*m2*Rqru*sr*lam(7)*lam(8) - 144*lam(8)**2 -                   &
               144*Rqsd*lam(8)**2 + 288*qzt2*sr*lam(8)**2 +                     &
               576*qzt2*Rqru*sr*lam(8)**2 +                                     &
               288*qzt2*Rqru**2*sr*lam(8)**2 +                                  &
               288*qzt2*Rqsd*sr*lam(8)**2 +                                     &
               576*qzt2*Rqru*Rqsd*sr*lam(8)**2 +                                &
               288*qzt2*Rqru**2*Rqsd*sr*lam(8)**2 -                             &
               24*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*lam(3)*             &
                (3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) + 14*lam(8)) -         &
               960*m2*(1 + Rqru)*sr*lam(2)*lam(9) - 60*lam(4)*lam(9) -          &
               60*Rqsd*lam(4)*lam(9) + 120*qzt2*sr*lam(4)*lam(9) +              &
               240*qzt2*Rqru*sr*lam(4)*lam(9) +                                 &
               120*qzt2*Rqru**2*sr*lam(4)*lam(9) +                              &
               120*qzt2*Rqsd*sr*lam(4)*lam(9) +                                 &
               240*qzt2*Rqru*Rqsd*sr*lam(4)*lam(9) +                            &
               120*qzt2*Rqru**2*Rqsd*sr*lam(4)*lam(9) -                         &
               60*lam(5)*lam(9) - 60*Rqsd*lam(5)*lam(9) +                       &
               120*qzt2*sr*lam(5)*lam(9) +                                      &
               240*qzt2*Rqru*sr*lam(5)*lam(9) +                                 &
               120*qzt2*Rqru**2*sr*lam(5)*lam(9) +                              &
               120*qzt2*Rqsd*sr*lam(5)*lam(9) +                                 &
               240*qzt2*Rqru*Rqsd*sr*lam(5)*lam(9) +                            &
               120*qzt2*Rqru**2*Rqsd*sr*lam(5)*lam(9) +                         &
               60*lam(6)*lam(9) + 60*Rqsd*lam(6)*lam(9) -                       &
               120*qzt2*sr*lam(6)*lam(9) -                                      &
               240*qzt2*Rqru*sr*lam(6)*lam(9) -                                 &
               120*qzt2*Rqru**2*sr*lam(6)*lam(9) -                              &
               120*qzt2*Rqsd*sr*lam(6)*lam(9) -                                 &
               240*qzt2*Rqru*Rqsd*sr*lam(6)*lam(9) -                            &
               120*qzt2*Rqru**2*Rqsd*sr*lam(6)*lam(9) +                         &
               60*lam(7)*lam(9) + 60*Rqsd*lam(7)*lam(9) -                       &
               120*qzt2*sr*lam(7)*lam(9) -                                      &
               240*qzt2*Rqru*sr*lam(7)*lam(9) -                                 &
               120*qzt2*Rqru**2*sr*lam(7)*lam(9) -                              &
               120*qzt2*Rqsd*sr*lam(7)*lam(9) -                                 &
               240*qzt2*Rqru*Rqsd*sr*lam(7)*lam(9) -                            &
               120*qzt2*Rqru**2*Rqsd*sr*lam(7)*lam(9) -                         &
               480*m2*sr*lam(8)*lam(9) - 480*m2*Rqru*sr*lam(8)*lam(9) -         &
               100*lam(9)**2 - 100*Rqsd*lam(9)**2 +                             &
               200*qzt2*sr*lam(9)**2 + 400*qzt2*Rqru*sr*lam(9)**2 +             &
               200*qzt2*Rqru**2*sr*lam(9)**2 +                                  &
               200*qzt2*Rqsd*sr*lam(9)**2 +                                     &
               400*qzt2*Rqru*Rqsd*sr*lam(9)**2 +                                &
               200*qzt2*Rqru**2*Rqsd*sr*lam(9)**2 +                             &
               24*(24*m2*(1 + Rqru)*sr*(2*lam(2) + lam(8)) -                    &
                  (1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                    &
                   (3*lam(4) + 3*lam(5) - 3*(lam(6) + lam(7)) +                 &
                     10*lam(9)))*lam(10) +                                      &
               144*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*lam(10)**2))       &
+ qdB*(expuISs*quB*qzt2*tanhus*                                                 &
             (9*dRqsu*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                &
                (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2)         &
- dRqru*(576*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*lam(2)**2 -              &
                  288*m2*(1 + Rqru)*sr*lam(2)*lam(4) - 9*lam(4)**2 -            &
                  9*Rqsu*lam(4)**2 + 18*qzt2*sr*lam(4)**2 +                     &
                  36*qzt2*Rqru*sr*lam(4)**2 +                                   &
                  18*qzt2*Rqru**2*sr*lam(4)**2 +                                &
                  18*qzt2*Rqsu*sr*lam(4)**2 +                                   &
                  36*qzt2*Rqru*Rqsu*sr*lam(4)**2 +                              &
                  18*qzt2*Rqru**2*Rqsu*sr*lam(4)**2 +                           &
                  1440*m2*(1 + Rqru)*sr*lam(2)*lam(5) +                         &
                  90*lam(4)*lam(5) + 90*Rqsu*lam(4)*lam(5) -                    &
                  180*qzt2*sr*lam(4)*lam(5) -                                   &
                  360*qzt2*Rqru*sr*lam(4)*lam(5) -                              &
                  180*qzt2*Rqru**2*sr*lam(4)*lam(5) -                           &
                  180*qzt2*Rqsu*sr*lam(4)*lam(5) -                              &
                  360*qzt2*Rqru*Rqsu*sr*lam(4)*lam(5) -                         &
                  180*qzt2*Rqru**2*Rqsu*sr*lam(4)*lam(5) -                      &
                  117*lam(5)**2 - 117*Rqsu*lam(5)**2 +                          &
                  234*qzt2*sr*lam(5)**2 + 468*qzt2*Rqru*sr*lam(5)**2 +          &
                  234*qzt2*Rqru**2*sr*lam(5)**2 +                               &
                  234*qzt2*Rqsu*sr*lam(5)**2 +                                  &
                  468*qzt2*Rqru*Rqsu*sr*lam(5)**2 +                             &
                  234*qzt2*Rqru**2*Rqsu*sr*lam(5)**2 +                          &
                  288*m2*(1 + Rqru)*sr*lam(2)*lam(6) +                          &
                  18*lam(4)*lam(6) + 18*Rqsu*lam(4)*lam(6) -                    &
                  36*qzt2*sr*lam(4)*lam(6) -                                    &
                  72*qzt2*Rqru*sr*lam(4)*lam(6) -                               &
                  36*qzt2*Rqru**2*sr*lam(4)*lam(6) -                            &
                  36*qzt2*Rqsu*sr*lam(4)*lam(6) -                               &
                  72*qzt2*Rqru*Rqsu*sr*lam(4)*lam(6) -                          &
                  36*qzt2*Rqru**2*Rqsu*sr*lam(4)*lam(6) -                       &
                  90*lam(5)*lam(6) - 90*Rqsu*lam(5)*lam(6) +                    &
                  180*qzt2*sr*lam(5)*lam(6) +                                   &
                  360*qzt2*Rqru*sr*lam(5)*lam(6) +                              &
                  180*qzt2*Rqru**2*sr*lam(5)*lam(6) +                           &
                  180*qzt2*Rqsu*sr*lam(5)*lam(6) +                              &
                  360*qzt2*Rqru*Rqsu*sr*lam(5)*lam(6) +                         &
                  180*qzt2*Rqru**2*Rqsu*sr*lam(5)*lam(6) - 9*lam(6)**2 -        &
                  9*Rqsu*lam(6)**2 + 18*qzt2*sr*lam(6)**2 +                     &
                  36*qzt2*Rqru*sr*lam(6)**2 +                                   &
                  18*qzt2*Rqru**2*sr*lam(6)**2 +                                &
                  18*qzt2*Rqsu*sr*lam(6)**2 +                                   &
                  36*qzt2*Rqru*Rqsu*sr*lam(6)**2 +                              &
                  18*qzt2*Rqru**2*Rqsu*sr*lam(6)**2 -                           &
                  1440*m2*(1 + Rqru)*sr*lam(2)*lam(7) -                         &
                  90*lam(4)*lam(7) - 90*Rqsu*lam(4)*lam(7) +                    &
                  180*qzt2*sr*lam(4)*lam(7) +                                   &
                  360*qzt2*Rqru*sr*lam(4)*lam(7) +                              &
                  180*qzt2*Rqru**2*sr*lam(4)*lam(7) +                           &
                  180*qzt2*Rqsu*sr*lam(4)*lam(7) +                              &
                  360*qzt2*Rqru*Rqsu*sr*lam(4)*lam(7) +                         &
                  180*qzt2*Rqru**2*Rqsu*sr*lam(4)*lam(7) +                      &
                  234*lam(5)*lam(7) + 234*Rqsu*lam(5)*lam(7) -                  &
                  468*qzt2*sr*lam(5)*lam(7) -                                   &
                  936*qzt2*Rqru*sr*lam(5)*lam(7) -                              &
                  468*qzt2*Rqru**2*sr*lam(5)*lam(7) -                           &
                  468*qzt2*Rqsu*sr*lam(5)*lam(7) -                              &
                  936*qzt2*Rqru*Rqsu*sr*lam(5)*lam(7) -                         &
                  468*qzt2*Rqru**2*Rqsu*sr*lam(5)*lam(7) +                      &
                  90*lam(6)*lam(7) + 90*Rqsu*lam(6)*lam(7) -                    &
                  180*qzt2*sr*lam(6)*lam(7) -                                   &
                  360*qzt2*Rqru*sr*lam(6)*lam(7) -                              &
                  180*qzt2*Rqru**2*sr*lam(6)*lam(7) -                           &
                  180*qzt2*Rqsu*sr*lam(6)*lam(7) -                              &
                  360*qzt2*Rqru*Rqsu*sr*lam(6)*lam(7) -                         &
                  180*qzt2*Rqru**2*Rqsu*sr*lam(6)*lam(7) -                      &
                  117*lam(7)**2 - 117*Rqsu*lam(7)**2 +                          &
                  234*qzt2*sr*lam(7)**2 + 468*qzt2*Rqru*sr*lam(7)**2 +          &
                  234*qzt2*Rqru**2*sr*lam(7)**2 +                               &
                  234*qzt2*Rqsu*sr*lam(7)**2 +                                  &
                  468*qzt2*Rqru*Rqsu*sr*lam(7)**2 +                             &
                  234*qzt2*Rqru**2*Rqsu*sr*lam(7)**2 -                          &
                  144*lam(1)*lam(8) - 144*Rqsu*lam(1)*lam(8) +                  &
                  288*qzt2*sr*lam(1)*lam(8) +                                   &
                  576*qzt2*Rqru*sr*lam(1)*lam(8) +                              &
                  288*qzt2*Rqru**2*sr*lam(1)*lam(8) +                           &
                  288*qzt2*Rqsu*sr*lam(1)*lam(8) +                              &
                  576*qzt2*Rqru*Rqsu*sr*lam(1)*lam(8) +                         &
                  288*qzt2*Rqru**2*Rqsu*sr*lam(1)*lam(8) -                      &
                  1152*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*lam(2)*        &
                   lam(8) + 288*m2*sr*lam(4)*lam(8) +                           &
                  288*m2*Rqru*sr*lam(4)*lam(8) -                                &
                  576*m2*sr*lam(5)*lam(8) -                                     &
                  576*m2*Rqru*sr*lam(5)*lam(8) -                                &
                  288*m2*sr*lam(6)*lam(8) -                                     &
                  288*m2*Rqru*sr*lam(6)*lam(8) +                                &
                  576*m2*sr*lam(7)*lam(8) +                                     &
                  576*m2*Rqru*sr*lam(7)*lam(8) - 144*lam(8)**2 -                &
                  144*Rqsu*lam(8)**2 + 288*qzt2*sr*lam(8)**2 +                  &
                  576*qzt2*Rqru*sr*lam(8)**2 +                                  &
                  288*qzt2*Rqru**2*sr*lam(8)**2 +                               &
                  288*qzt2*Rqsu*sr*lam(8)**2 +                                  &
                  576*qzt2*Rqru*Rqsu*sr*lam(8)**2 +                             &
                  288*qzt2*Rqru**2*Rqsu*sr*lam(8)**2 +                          &
                  12*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*lam(3)*          &
                   (3*lam(4) - 3*(lam(5) - lam(6) + lam(7)) + 14*lam(8))        &
+ 1920*m2*(1 + Rqru)*sr*lam(2)*lam(9) + 120*lam(4)*lam(9) +                     &
                  120*Rqsu*lam(4)*lam(9) - 240*qzt2*sr*lam(4)*lam(9) -          &
                  480*qzt2*Rqru*sr*lam(4)*lam(9) -                              &
                  240*qzt2*Rqru**2*sr*lam(4)*lam(9) -                           &
                  240*qzt2*Rqsu*sr*lam(4)*lam(9) -                              &
                  480*qzt2*Rqru*Rqsu*sr*lam(4)*lam(9) -                         &
                  240*qzt2*Rqru**2*Rqsu*sr*lam(4)*lam(9) -                      &
                  240*lam(5)*lam(9) - 240*Rqsu*lam(5)*lam(9) +                  &
                  480*qzt2*sr*lam(5)*lam(9) +                                   &
                  960*qzt2*Rqru*sr*lam(5)*lam(9) +                              &
                  480*qzt2*Rqru**2*sr*lam(5)*lam(9) +                           &
                  480*qzt2*Rqsu*sr*lam(5)*lam(9) +                              &
                  960*qzt2*Rqru*Rqsu*sr*lam(5)*lam(9) +                         &
                  480*qzt2*Rqru**2*Rqsu*sr*lam(5)*lam(9) -                      &
                  120*lam(6)*lam(9) - 120*Rqsu*lam(6)*lam(9) +                  &
                  240*qzt2*sr*lam(6)*lam(9) +                                   &
                  480*qzt2*Rqru*sr*lam(6)*lam(9) +                              &
                  240*qzt2*Rqru**2*sr*lam(6)*lam(9) +                           &
                  240*qzt2*Rqsu*sr*lam(6)*lam(9) +                              &
                  480*qzt2*Rqru*Rqsu*sr*lam(6)*lam(9) +                         &
                  240*qzt2*Rqru**2*Rqsu*sr*lam(6)*lam(9) +                      &
                  240*lam(7)*lam(9) + 240*Rqsu*lam(7)*lam(9) -                  &
                  480*qzt2*sr*lam(7)*lam(9) -                                   &
                  960*qzt2*Rqru*sr*lam(7)*lam(9) -                              &
                  480*qzt2*Rqru**2*sr*lam(7)*lam(9) -                           &
                  480*qzt2*Rqsu*sr*lam(7)*lam(9) -                              &
                  960*qzt2*Rqru*Rqsu*sr*lam(7)*lam(9) -                         &
                  480*qzt2*Rqru**2*Rqsu*sr*lam(7)*lam(9) -                      &
                  480*m2*sr*lam(8)*lam(9) -                                     &
                  480*m2*Rqru*sr*lam(8)*lam(9) - 100*lam(9)**2 -                &
                  100*Rqsu*lam(9)**2 + 200*qzt2*sr*lam(9)**2 +                  &
                  400*qzt2*Rqru*sr*lam(9)**2 +                                  &
                  200*qzt2*Rqru**2*sr*lam(9)**2 +                               &
                  200*qzt2*Rqsu*sr*lam(9)**2 +                                  &
                  400*qzt2*Rqru*Rqsu*sr*lam(9)**2 +                             &
                  200*qzt2*Rqru**2*Rqsu*sr*lam(9)**2 +                          &
                  12*(24*m2*(1 + Rqru)*sr*(6*lam(2) - 5*lam(8)) -               &
                     (1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                 &
                      (9*lam(4) - 39*lam(5) - 9*lam(6) + 39*lam(7) -            &
                        50*lam(9)))*lam(10) +                                   &
                  312*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                &
                   lam(10)**2)) +                                               &
            2*expuISs*qxy2*(1 + Rqru)*(1 + Rqsu)*tanhus**2*                     &
             (9*dRqsu*qzt2*(1 + Rqsu)*                                          &
                (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2)         &
+ dRqru*qxy2*(1 + Rqru)*(9*lam(4)**2 + 117*lam(5)**2 +                          &
                  90*lam(5)*lam(6) + 9*lam(6)**2 -                              &
                  36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -               &
                  234*lam(5)*lam(7) - 90*lam(6)*lam(7) + 117*lam(7)**2 -        &
                  144*lam(1)*lam(8) +                                           &
                  24*(24*lam(2)**2 - 48*lam(2)*lam(8) +                         &
                     lam(8)*(-7*lam(3) + 6*lam(8))) +                           &
                  240*lam(5)*lam(9) + 120*lam(6)*lam(9) -                       &
                  240*lam(7)*lam(9) + 100*lam(9)**2 +                           &
                  12*(39*lam(5) + 9*lam(6) - 39*lam(7) + 50*lam(9))*            &
                   lam(10) + 312*lam(10)**2 -                                   &
                  6*lam(4)*(15*lam(5) + 3*lam(6) - 15*lam(7) +                  &
                     20*lam(9) + 18*lam(10)))) +                                &
            2*dRqru*qxy2*(1 + Rqru)*                                            &
             (expdISs*(-24*m2*                                                  &
                   (24*lam(1)*lam(3) + 14*lam(3)**2 +                           &
                     lam(8)*(9*lam(4) - 3*lam(5) + 3*lam(6) -                   &
                        9*lam(7) + 14*lam(8) + 10*lam(9)) +                     &
                     lam(2)*(6*lam(4) + 6*lam(5) -                              &
                        6*(lam(6) + lam(7)) + 20*lam(9))) +                     &
                  qxy2*(1 + Rqru)*(1 + Rqsd)*                                   &
                   (576*lam(2)**2 + 576*lam(2)*lam(8) +                         &
                     144*lam(8)*(2*lam(1) + lam(8)) +                           &
                     24*lam(3)*                                                 &
                      (3*(lam(4) - lam(5) + lam(6) - lam(7)) +                  &
                        14*lam(8)) +                                            &
                     (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +               &
                        10*lam(9) - 12*lam(10))**2) +                           &
                  qzt2*(1 + Rqru)*(1 + Rqsd)*                                   &
                   (576*lam(2)**2 + 576*lam(2)*lam(8) +                         &
                     144*lam(8)*(2*lam(1) + lam(8)) +                           &
                     24*lam(3)*                                                 &
                      (3*(lam(4) - lam(5) + lam(6) - lam(7)) +                  &
                        14*lam(8)) +                                            &
                     (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +               &
                        10*lam(9) - 12*lam(10))**2) +                           &
                  288*m2*(2*lam(2) + lam(8))*lam(10)) -                         &
               expuISs*(qxy2*(1 + Rqru)*(1 + Rqsu)*                             &
                   (9*lam(4)**2 + 117*lam(5)**2 + 90*lam(5)*lam(6) +            &
                     9*lam(6)**2 -                                              &
                     36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -            &
                     234*lam(5)*lam(7) - 90*lam(6)*lam(7) +                     &
                     117*lam(7)**2 - 144*lam(1)*lam(8) +                        &
                     24*(24*lam(2)**2 - 48*lam(2)*lam(8) +                      &
                        lam(8)*(-7*lam(3) + 6*lam(8))) +                        &
                     240*lam(5)*lam(9) + 120*lam(6)*lam(9) -                    &
                     240*lam(7)*lam(9) + 100*lam(9)**2 +                        &
                     12*(39*lam(5) + 9*lam(6) - 39*lam(7) + 50*lam(9))*         &
                      lam(10) + 312*lam(10)**2 -                                &
                     6*lam(4)*                                                  &
                      (15*lam(5) + 3*lam(6) - 15*lam(7) + 20*lam(9) +           &
                        18*lam(10))) +                                          &
                  qzt2*(1 + Rqru)*(1 + Rqsu)*                                   &
                   (9*lam(4)**2 + 117*lam(5)**2 + 90*lam(5)*lam(6) +            &
                     9*lam(6)**2 -                                              &
                     36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -            &
                     234*lam(5)*lam(7) - 90*lam(6)*lam(7) +                     &
                     117*lam(7)**2 - 144*lam(1)*lam(8) +                        &
                     24*(24*lam(2)**2 - 48*lam(2)*lam(8) +                      &
                        lam(8)*(-7*lam(3) + 6*lam(8))) +                        &
                     240*lam(5)*lam(9) + 120*lam(6)*lam(9) -                    &
                     240*lam(7)*lam(9) + 100*lam(9)**2 +                        &
                     12*(39*lam(5) + 9*lam(6) - 39*lam(7) + 50*lam(9))*         &
                      lam(10) + 312*lam(10)**2 -                                &
                     6*lam(4)*                                                  &
                      (15*lam(5) + 3*lam(6) - 15*lam(7) + 20*lam(9) +           &
                        18*lam(10))) +                                          &
                  24*m2*(12*lam(1)*lam(3) + 7*lam(3)**2 +                       &
                     lam(2)*(-6*lam(4) + 30*lam(5) + 6*lam(6) -                 &
                        30*lam(7) + 40*lam(9) + 36*lam(10)) +                   &
                     lam(8)*(9*lam(4) - 15*lam(5) - 3*lam(6) + 9*lam(7) +       &
                        7*lam(8) - 10*(lam(9) + 3*lam(10)))))))) -              &
      expuISr*qxy2*tanhur**2*(6*dRqsd*expdISs*quB*qxy2*(1 + Rqru)*              &
          (1 + Rqsd)**2*tanhds**3*                                              &
          (48*lam(2)**2 + 3*(lam(4)**2 + 13*lam(5)**2 -                         &
               10*lam(5)*lam(6) + lam(6)**2 +                                   &
               2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                       &
               26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +            &
            144*lam(2)*lam(8) -                                                 &
            12*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) +                 &
            80*lam(8)**2) - 2*expdISs*(1 + Rqru)*tanhds*                        &
          (3*dRqsd*quB*qxy2*(1 + Rqsd)**2*                                      &
             (48*lam(2)**2 + 3*                                                 &
                (lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                  &
                  lam(6)**2 +                                                   &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                    &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               144*lam(2)*lam(8) -                                              &
               12*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) +              &
               80*lam(8)**2) +                                                  &
            dRqru*qdB*(qzt2*(1 + Rqru)*(1 + Rqsd)*                              &
                (576*lam(2)**2 + 576*lam(2)*lam(8) +                            &
                  144*lam(8)*(-2*lam(1) + lam(8)) -                             &
                  24*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +            &
                     14*lam(8)) +                                               &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) -                              &
               24*m2*(2*lam(2) + lam(8))*                                       &
                (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -        &
                  12*lam(10)))) -                                               &
         expdISs*qdB*quB*tanhds**2*                                             &
          (-3*dRqsd*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                  &
             (48*lam(2)**2 + 3*                                                 &
                (lam(4)**2 + 13*lam(5)**2 - 10*lam(5)*lam(6) +                  &
                  lam(6)**2 +                                                   &
                  2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                    &
                  26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2) +         &
               144*lam(2)*lam(8) -                                              &
               12*(lam(4) - 3*lam(5) + lam(6) - 3*lam(7))*lam(8) +              &
               80*lam(8)**2) +                                                  &
            dRqru*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                    &
             (576*lam(2)**2 + 576*lam(2)*lam(8) +                               &
               144*lam(8)*(2*lam(1) + lam(8)) +                                 &
               24*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +               &
                  14*lam(8)) +                                                  &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -         &
                  12*lam(10))**2)) +                                            &
         qdB*(-3*(1 + Rqru)*(3*dRqsu*expuISs*(-1 + tanhus)*(1 + tanhus)*        &
                (quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +                           &
                  2*qxy2*(1 + Rqsu)**2*tanhus)*                                 &
                (16*lam(2)**2 + (lam(4) + lam(5) + lam(6) + lam(7))**2) +       &
               dRqsd*expdISs*quB*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                &
                (3*(16*lam(2)**2 + lam(4)**2 + 13*lam(5)**2 -                   &
                     10*lam(5)*lam(6) + lam(6)**2 +                             &
                     2*lam(4)*(-5*lam(5) + lam(6) - 5*lam(7)) +                 &
                     26*lam(5)*lam(7) - 10*lam(6)*lam(7) + 13*lam(7)**2)        &
+ 12*(12*lam(2) - lam(4) + 3*lam(5) - lam(6) + 3*lam(7))*lam(8) +               &
                  80*lam(8)**2)) +                                              &
            dRqru*(expdISs*quB*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*       &
                (576*lam(2)**2 + 576*lam(2)*lam(8) +                            &
                  144*lam(8)*(2*lam(1) + lam(8)) +                              &
                  24*lam(3)*(3*(lam(4) - lam(5) + lam(6) - lam(7)) +            &
                     14*lam(8)) +                                               &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                     10*lam(9) - 12*lam(10))**2) +                              &
               expuISs*(quB*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*          &
                   (-1 + tanhus**2)*                                            &
                   (9*lam(4)**2 + 117*lam(5)**2 + 90*lam(5)*lam(6) +            &
                     9*lam(6)**2 -                                              &
                     36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -            &
                     234*lam(5)*lam(7) - 90*lam(6)*lam(7) +                     &
                     117*lam(7)**2 - 144*lam(1)*lam(8) +                        &
                     24*(24*lam(2)**2 - 48*lam(2)*lam(8) +                      &
                        lam(8)*(-7*lam(3) + 6*lam(8))) +                        &
                     240*lam(5)*lam(9) + 120*lam(6)*lam(9) -                    &
                     240*lam(7)*lam(9) + 100*lam(9)**2 +                        &
                     12*(39*lam(5) + 9*lam(6) - 39*lam(7) + 50*lam(9))*         &
                      lam(10) + 312*lam(10)**2 -                                &
                     6*lam(4)*                                                  &
                      (15*lam(5) + 3*lam(6) - 15*lam(7) + 20*lam(9) +           &
                        18*lam(10))) +                                          &
                  2*(1 + Rqru)*tanhus*                                          &
                   (qzt2*(1 + Rqru)*(1 + Rqsu)*                                 &
                      (9*lam(4)**2 + 117*lam(5)**2 + 90*lam(5)*lam(6) +         &
                        9*lam(6)**2 +                                           &
                        36*lam(3)*(lam(4) - lam(5) + lam(6) - lam(7)) -         &
                        234*lam(5)*lam(7) - 90*lam(6)*lam(7) +                  &
                        117*lam(7)**2 + 144*lam(1)*lam(8) +                     &
                        24*(24*lam(2)**2 - 48*lam(2)*lam(8) +                   &
                         lam(8)*(7*lam(3) + 6*lam(8))) +                        &
                        240*lam(5)*lam(9) + 120*lam(6)*lam(9) -                 &
                        240*lam(7)*lam(9) + 100*lam(9)**2 +                     &
                        12*(39*lam(5) + 9*lam(6) - 39*lam(7) +                  &
                        50*lam(9))*lam(10) + 312*lam(10)**2 -                   &
                        6*lam(4)*                                               &
                         (15*lam(5) + 3*lam(6) - 15*lam(7) + 20*lam(9) +        &
                         18*lam(10))) +                                         &
                     48*m2*(lam(2)*                                             &
                         (-3*lam(4) + 15*lam(5) + 3*lam(6) - 15*lam(7) +        &
                         20*lam(9) + 18*lam(10)) +                              &
                        lam(8)*                                                 &
                         (3*lam(4) - 6*lam(5) - 3*lam(6) + 6*lam(7) -           &
                          5*(lam(9) + 3*lam(10))))))))))/216.
  dtLAM(9)=(6*dRqsd*expdISs*quB*qxy2**2*(1 + Rqsd)**2*tanhds**3*             &
       (2*expdISr*(1 + Rqrd)*(3*lam(1) + lam(3))*                               &
          (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) -               &
         expuISr*(1 + Rqru)*(lam(3)*                                            &
             (-21*lam(4) - 51*lam(5) + 21*lam(6) + 51*lam(7) +                  &
               96*lam(9) - 52*lam(10)) +                                        &
            6*lam(1)*(3*lam(4) + 15*lam(5) -                                    &
               3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10)))) -               &
      18*dRqrd*expdISr*quB*qxy2**2*(1 + Rqrd)**2*tanhdr**3*                     &
       (-(expuISs*(1 + Rqsu)*(-1 + tanhus**2)*(2*lam(2) + lam(8))*              &
            (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -            &
              12*lam(10))) + 2*expdISs*(1 + Rqsd)*(-1 + tanhds**2)*             &
          (lam(2)*(3*lam(4) - 15*lam(5) - 3*lam(6) + 15*lam(7) -                &
               20*lam(9) - 18*lam(10)) +                                        &
            lam(8)*(-3*lam(4) + 6*lam(5) + 3*lam(6) - 6*lam(7) +                &
               5*lam(9) + 15*lam(10)))) +                                       &
      18*dRqru*expuISr*qdB*qxy2**2*(1 + Rqru)**2*tanhur**3*                     &
       (expdISs*(1 + Rqsd)*(-1 + tanhds**2)*(2*lam(2) + lam(8))*                &
          (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -              &
            12*lam(10)) - 2*expuISs*(1 + Rqsu)*(-1 + tanhus**2)*                &
          (lam(2)*(3*lam(4) - 15*lam(5) - 3*lam(6) + 15*lam(7) -                &
               20*lam(9) - 18*lam(10)) +                                        &
            lam(8)*(-3*lam(4) + 6*lam(5) + 3*lam(6) - 6*lam(7) +                &
               5*lam(9) + 15*lam(10)))) -                                       &
      6*dRqsd*expdISs*quB*qxy2*(1 + Rqsd)*tanhds*                               &
       (expdISr*(2*(qxy2 + qzt2)*(1 + Rqrd)*(1 + Rqsd)*                         &
             (3*lam(1) + lam(3))*                                               &
             (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)) +              &
            m2*(8*(3*lam(1) + lam(3))**2 +                                      &
               9*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) -                  &
               6*lam(10)**2)) -                                                 &
         expuISr*(m2*(72*lam(1)**2 - 168*lam(1)*lam(3) + 44*lam(3)**2 -         &
               9*(3*lam(4) + 9*lam(5) - 3*lam(6) - 9*lam(7) - 4*lam(9))*        &
                lam(9) + 9*(lam(4) + 5*lam(5) - lam(6) - 5*lam(7) -             &
                  4*lam(9))*lam(10) + 12*lam(10)**2) +                          &
            (qxy2 + qzt2)*(1 + Rqru)*(1 + Rqsd)*                                &
             (lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) + 51*lam(7) +          &
                  96*lam(9) - 52*lam(10)) +                                     &
               6*lam(1)*(3*lam(4) + 15*lam(5) -                                 &
                  3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10))))) +           &
      3*expdISs*qdB*quB*qxy2*tanhds**2*                                         &
       (dRqsd*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                                   &
          (2*expdISr*(1 + Rqrd)*(3*lam(1) + lam(3))*                            &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) -            &
            expuISr*(1 + Rqru)*                                                 &
             (lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) + 51*lam(7) +          &
                  96*lam(9) - 52*lam(10)) +                                     &
               6*lam(1)*(3*lam(4) + 15*lam(5) -                                 &
                  3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10)))) +            &
         3*(1 + Rqsd)*(-(dRqru*expuISr*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*          &
               (2*lam(2) + lam(8))*                                             &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -         &
                 12*lam(10))) +                                                 &
            2*dRqrd*expdISr*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                     &
             (lam(2)*(3*lam(4) - 15*lam(5) - 3*lam(6) + 15*lam(7) -             &
                  20*lam(9) - 18*lam(10)) +                                     &
               lam(8)*(-3*lam(4) + 6*lam(5) + 3*lam(6) - 6*lam(7) +             &
                  5*lam(9) + 15*lam(10))))) +                                   &
      expdISr*tanhdr*(-(expuISs*qdB*quB*qzt2*tanhus*                            &
            (-9*dRqrd*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                &
               (2*lam(2) + lam(8))*                                             &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -         &
                 12*lam(10)) +                                                  &
              dRqsu*m2*(1 + Rqsu)*ss*                                           &
               (432*lam(1)**2 - 1008*lam(1)*lam(3) + 264*lam(3)**2 +            &
                 9*lam(4)**2 +                                                  &
                 9*(13*lam(5)**2 + lam(6)**2 + 10*lam(6)*lam(7) +               &
                    13*lam(7)**2 - 2*lam(5)*(5*lam(6) + 13*lam(7))) +           &
                 90*(-3*lam(5) + lam(6) + 3*lam(7))*lam(9) +                    &
                 312*lam(9)**2 +                                                &
                 18*lam(4)*(5*lam(5) - lam(6) - 5*(lam(7) + lam(9))) +          &
                 6*(7*lam(4) + 23*lam(5) - 7*lam(6) - 23*lam(7) -               &
                    52*lam(9))*lam(10) + 64*lam(10)**2) +                       &
              dRqrd*m2*(1 + Rqrd)*sr*                                           &
               (9*lam(4)**2 + 9*lam(5)**2 + 9*lam(6)**2 +                       &
                 18*lam(6)*lam(7) + 9*lam(7)**2 +                               &
                 108*(2*lam(2) + lam(8))**2 - 6*lam(6)*lam(9) -                 &
                 6*lam(7)*lam(9) - 8*lam(9)**2 +                                &
                 6*lam(4)*(3*lam(5) - 3*lam(6) - 3*lam(7) + lam(9) -            &
                    9*lam(10)) +                                                &
                 18*(3*(lam(6) + lam(7)) - 10*lam(9))*lam(10) +                 &
                 108*lam(10)**2 -                                               &
                 6*lam(5)*(3*lam(6) + 3*lam(7) - lam(9) + 9*lam(10))) +         &
              3*dRqsu*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                &
               (lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) +                    &
                    51*lam(7) + 96*lam(9) - 52*lam(10)) +                       &
                 6*lam(1)*(3*lam(4) + 15*lam(5) -                               &
                    3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10))))) +         &
         expdISs*qdB*quB*qzt2*tanhds*                                           &
          (6*dRqsd*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                   &
             (3*lam(1) + lam(3))*                                               &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +            &
            dRqsd*m2*(1 + Rqsd)*ss*                                             &
             (48*(3*lam(1) + lam(3))**2 +                                       &
               9*(lam(4) - lam(5) - lam(6) + lam(7))**2 +                       &
               6*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +                  &
               28*lam(10)**2) -                                                 &
            18*dRqrd*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                 &
             (lam(2)*(3*lam(4) - 15*lam(5) - 3*lam(6) + 15*lam(7) -             &
                  20*lam(9) - 18*lam(10)) +                                     &
               lam(8)*(-3*lam(4) + 6*lam(5) + 3*lam(6) - 6*lam(7) +             &
                  5*lam(9) + 15*lam(10))) +                                     &
            dRqrd*m2*(1 + Rqrd)*sr*                                             &
             (9*lam(4)**2 + 117*lam(5)**2 + 9*lam(6)**2 -                       &
               90*lam(6)*lam(7) + 117*lam(7)**2 +                               &
               108*(4*lam(2)**2 - 8*lam(2)*lam(8) + lam(8)**2) +                &
               66*lam(6)*lam(9) - 186*lam(7)*lam(9) + 118*lam(9)**2 +           &
               6*(15*lam(6) - 63*lam(7) + 74*lam(9))*lam(10) +                  &
               240*lam(10)**2 -                                                 &
               6*lam(4)*(15*lam(5) + 3*lam(6) - 15*lam(7) + 11*lam(9) +         &
                  15*lam(10)) +                                                 &
               6*lam(5)*(15*lam(6) - 39*lam(7) + 31*lam(9) + 63*lam(10)))       &
) - 18*dRqrd*quB*qxy2*(1 + Rqrd)*                                               &
          (-(expdISs*m2*(6*(4*lam(2)**2 - 8*lam(2)*lam(8) +                     &
                    lam(8)**2) -                                                &
                 3*lam(9)*(3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) +            &
                    lam(9)) +                                                   &
                 (-3*lam(4) + 3*(5*lam(5) + lam(6) - 5*lam(7)) +                &
                    26*lam(9))*lam(10) + 12*lam(10)**2)) +                      &
            2*expdISs*(qxy2 + qzt2)*(1 + Rqrd)*(1 + Rqsd)*                      &
             (lam(2)*(3*lam(4) - 15*lam(5) - 3*lam(6) + 15*lam(7) -             &
                  20*lam(9) - 18*lam(10)) +                                     &
               lam(8)*(-3*lam(4) + 6*lam(5) + 3*lam(6) - 6*lam(7) +             &
                  5*lam(9) + 15*lam(10))) +                                     &
            expuISs*(-((qxy2 + qzt2)*(1 + Rqrd)*(1 + Rqsu)*                     &
                  (2*lam(2) + lam(8))*                                          &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                    10*lam(9) - 12*lam(10))) +                                  &
               m2*(6*(2*lam(2) + lam(8))**2 +                                   &
                  9*lam(9)*(lam(4) + lam(5) - lam(6) - lam(7) +                 &
                     2*lam(9)) +                                                &
                  (-3*lam(4) - 3*lam(5) + 3*(lam(6) + lam(7)) -                 &
                     10*lam(9))*lam(10) + 6*lam(10)**2))) +                     &
         expdISs*quB*qxy2*(1 + Rqsd)*tanhds**2*                                 &
          (36*dRqrd*qxy2*(1 + Rqrd)**2*                                         &
             (lam(2)*(3*lam(4) - 15*lam(5) - 3*lam(6) + 15*lam(7) -             &
                  20*lam(9) - 18*lam(10)) +                                     &
               lam(8)*(-3*lam(4) + 6*lam(5) + 3*lam(6) - 6*lam(7) +             &
                  5*lam(9) + 15*lam(10))) +                                     &
            dRqsd*(12*qzt2*(1 + Rqrd)*(1 + Rqsd)*(3*lam(1) + lam(3))*           &
                (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)) +           &
               m2*(48*(3*lam(1) + lam(3))**2 +                                  &
                  9*(lam(4) - lam(5) - lam(6) + lam(7))**2 +                    &
                  6*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +               &
                  28*lam(10)**2))) -                                            &
         expuISs*qxy2*(1 + Rqsu)*tanhus**2*                                     &
          (18*dRqrd*quB*qxy2*(1 + Rqrd)**2*(2*lam(2) + lam(8))*                 &
             (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -           &
               12*lam(10)) + dRqsu*qdB*                                         &
             (m2*(24*(18*lam(1)**2 - 42*lam(1)*lam(3) + 11*lam(3)**2) +         &
                  9*lam(4)**2 +                                                 &
                  9*(13*lam(5)**2 + lam(6)**2 + 10*lam(6)*lam(7) +              &
                     13*lam(7)**2 - 2*lam(5)*(5*lam(6) + 13*lam(7))) +          &
                  90*(-3*lam(5) + lam(6) + 3*lam(7))*lam(9) +                   &
                  312*lam(9)**2 +                                               &
                  18*lam(4)*(5*lam(5) - lam(6) - 5*(lam(7) + lam(9))) +         &
                  6*(7*lam(4) + 23*lam(5) - 7*lam(6) - 23*lam(7) -              &
                     52*lam(9))*lam(10) + 64*lam(10)**2) +                      &
               6*qzt2*(1 + Rqrd)*(1 + Rqsu)*                                    &
                (lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) +                   &
                     51*lam(7) + 96*lam(9) - 52*lam(10)) +                      &
                  6*lam(1)*(3*lam(4) + 15*lam(5) -                              &
                     3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10)))))) +       &
      expdISr*qxy2*tanhdr**2*(-12*dRqsd*expdISs*quB*qxy2*(1 + Rqrd)*            &
          (1 + Rqsd)**2*tanhds**3*(3*lam(1) + lam(3))*                          &
          (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +               &
         6*dRqsu*expuISs*qdB*qxy2*(1 + Rqrd)*(1 + Rqsu)**2*tanhus**3*           &
          (lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) + 51*lam(7) +             &
               96*lam(9) - 52*lam(10)) +                                        &
            6*lam(1)*(3*lam(4) + 15*lam(5) -                                    &
               3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10))) +                &
         3*expuISs*qdB*quB*tanhus**2*                                           &
          (3*dRqrd*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                   &
             (2*lam(2) + lam(8))*                                               &
             (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -           &
               12*lam(10)) + dRqsu*(1 + Rqrd)*                                  &
             (-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                                    &
             (lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) + 51*lam(7) +          &
                  96*lam(9) - 52*lam(10)) +                                     &
               6*lam(1)*(3*lam(4) + 15*lam(5) -                                 &
                  3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10)))) -            &
         6*expdISs*qdB*quB*tanhds**2*                                           &
          (dRqsd*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                     &
             (3*lam(1) + lam(3))*                                               &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +            &
            3*dRqrd*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                  &
             (lam(2)*(3*lam(4) - 15*lam(5) - 3*lam(6) + 15*lam(7) -             &
                  20*lam(9) - 18*lam(10)) +                                     &
               lam(8)*(-3*lam(4) + 6*lam(5) + 3*lam(6) - 6*lam(7) +             &
                  5*lam(9) + 15*lam(10)))) +                                    &
         expdISs*quB*(1 + Rqrd)*tanhds*                                         &
          (12*dRqsd*qxy2*(1 + Rqsd)**2*(3*lam(1) + lam(3))*                     &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) -            &
            36*dRqrd*qzt2*(1 + Rqrd)*(1 + Rqsd)*                                &
             (lam(2)*(3*lam(4) - 15*lam(5) - 3*lam(6) + 15*lam(7) -             &
                  20*lam(9) - 18*lam(10)) +                                     &
               lam(8)*(-3*lam(4) + 6*lam(5) + 3*lam(6) - 6*lam(7) +             &
                  5*lam(9) + 15*lam(10))) +                                     &
            dRqrd*m2*(9*lam(4)**2 + 117*lam(5)**2 + 9*lam(6)**2 -               &
               90*lam(6)*lam(7) + 117*lam(7)**2 +                               &
               108*(4*lam(2)**2 - 8*lam(2)*lam(8) + lam(8)**2) +                &
               66*lam(6)*lam(9) - 186*lam(7)*lam(9) + 118*lam(9)**2 +           &
               6*(15*lam(6) - 63*lam(7) + 74*lam(9))*lam(10) +                  &
               240*lam(10)**2 -                                                 &
               6*lam(4)*(15*lam(5) + 3*lam(6) - 15*lam(7) + 11*lam(9) +         &
                  15*lam(10)) +                                                 &
               6*lam(5)*(15*lam(6) - 39*lam(7) + 31*lam(9) + 63*lam(10)))       &
) - expuISs*(1 + Rqrd)*tanhus*                                                  &
          (6*dRqsu*qdB*qxy2*(1 + Rqsu)**2*                                      &
             (lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) + 51*lam(7) +          &
                  96*lam(9) - 52*lam(10)) +                                     &
               6*lam(1)*(3*lam(4) + 15*lam(5) -                                 &
                  3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10))) +             &
            dRqrd*quB*(-18*qzt2*(1 + Rqrd)*(1 + Rqsu)*                          &
                (2*lam(2) + lam(8))*                                            &
                (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                    &
                  10*lam(9) - 12*lam(10)) +                                     &
               m2*(9*lam(4)**2 + 9*lam(5)**2 + 9*lam(6)**2 +                    &
                  18*lam(6)*lam(7) + 9*lam(7)**2 +                              &
                  108*(2*lam(2) + lam(8))**2 - 6*lam(6)*lam(9) -                &
                  6*lam(7)*lam(9) - 8*lam(9)**2 +                               &
                  6*lam(4)*(3*lam(5) - 3*lam(6) - 3*lam(7) + lam(9) -           &
                     9*lam(10)) +                                               &
                  18*(3*(lam(6) + lam(7)) - 10*lam(9))*lam(10) +                &
                  108*lam(10)**2 -                                              &
                  6*lam(5)*(3*lam(6) + 3*lam(7) - lam(9) + 9*lam(10)))))        &
+ 3*qdB*quB*(2*dRqsd*expdISs*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*         &
             (3*lam(1) + lam(3))*                                               &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) -            &
            dRqsu*expuISs*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*            &
             (lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) + 51*lam(7) +          &
                  96*lam(9) - 52*lam(10)) +                                     &
               6*lam(1)*(3*lam(4) + 15*lam(5) -                                 &
                  3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10))) +             &
            3*dRqrd*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                             &
             (-(expuISs*(1 + Rqsu)*(2*lam(2) + lam(8))*                         &
                  (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                  &
                    10*lam(9) - 12*lam(10))) +                                  &
               2*expdISs*(1 + Rqsd)*                                            &
                (lam(2)*(3*lam(4) - 15*lam(5) - 3*lam(6) + 15*lam(7) -          &
                     20*lam(9) - 18*lam(10)) +                                  &
                  lam(8)*(-3*lam(4) + 6*lam(5) + 3*lam(6) - 6*lam(7) +          &
                     5*lam(9) + 15*lam(10)))))) +                               &
      expuISr*qxy2*tanhur**2*(6*dRqsd*expdISs*quB*qxy2*(1 + Rqru)*              &
          (1 + Rqsd)**2*tanhds**3*                                              &
          (lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) + 51*lam(7) +             &
               96*lam(9) - 52*lam(10)) +                                        &
            6*lam(1)*(3*lam(4) + 15*lam(5) -                                    &
               3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10))) +                &
         3*expdISs*qdB*quB*tanhds**2*                                           &
          (3*dRqru*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                   &
             (2*lam(2) + lam(8))*                                               &
             (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -           &
               12*lam(10)) + dRqsd*(1 + Rqru)*                                  &
             (-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                                    &
             (lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) + 51*lam(7) +          &
                  96*lam(9) - 52*lam(10)) +                                     &
               6*lam(1)*(3*lam(4) + 15*lam(5) -                                 &
                  3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10)))) -            &
         expdISs*(1 + Rqru)*tanhds*                                             &
          (6*dRqsd*quB*qxy2*(1 + Rqsd)**2*                                      &
             (lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) + 51*lam(7) +          &
                  96*lam(9) - 52*lam(10)) +                                     &
               6*lam(1)*(3*lam(4) + 15*lam(5) -                                 &
                  3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10))) +             &
            dRqru*qdB*(-18*qzt2*(1 + Rqru)*(1 + Rqsd)*                          &
                (2*lam(2) + lam(8))*                                            &
                (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                    &
                  10*lam(9) - 12*lam(10)) +                                     &
               m2*(9*lam(4)**2 + 9*lam(5)**2 + 9*lam(6)**2 +                    &
                  18*lam(6)*lam(7) + 9*lam(7)**2 +                              &
                  108*(2*lam(2) + lam(8))**2 - 6*lam(6)*lam(9) -                &
                  6*lam(7)*lam(9) - 8*lam(9)**2 +                               &
                  6*lam(4)*(3*lam(5) - 3*lam(6) - 3*lam(7) + lam(9) -           &
                     9*lam(10)) +                                               &
                  18*(3*(lam(6) + lam(7)) - 10*lam(9))*lam(10) +                &
                  108*lam(10)**2 -                                              &
                  6*lam(5)*(3*lam(6) + 3*lam(7) - lam(9) + 9*lam(10)))))        &
+ qdB*(-3*(1 + Rqru)*(2*dRqsu*expuISs*(-1 + tanhus)*(1 + tanhus)*               &
                (quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +                           &
                  2*qxy2*(1 + Rqsu)**2*tanhus)*(3*lam(1) + lam(3))*             &
                (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +         &
               dRqsd*expdISs*quB*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                &
                (lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) +                   &
                     51*lam(7) + 96*lam(9) - 52*lam(10)) +                      &
                  6*lam(1)*(3*lam(4) + 15*lam(5) -                              &
                     3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10)))) +         &
            dRqru*(-9*expdISs*quB*(1 + Rqsd)*                                   &
                (-1 + 2*qzt2*(1 + Rqru)**2*sr)*(2*lam(2) + lam(8))*             &
                (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -        &
                  12*lam(10)) -                                                 &
               18*expuISs*(1 + Rqsu)*                                           &
                (2*qzt2*(1 + Rqru)**2*tanhus +                                  &
                  quB*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*(-1 + tanhus**2))*         &
                (lam(2)*(3*lam(4) - 15*lam(5) - 3*lam(6) + 15*lam(7) -          &
                     20*lam(9) - 18*lam(10)) +                                  &
                  lam(8)*(-3*lam(4) + 6*lam(5) + 3*lam(6) - 6*lam(7) +          &
                     5*lam(9) + 15*lam(10))) +                                  &
               expuISs*m2*(1 + Rqru)*tanhus*                                    &
                (9*lam(4)**2 + 117*lam(5)**2 + 9*lam(6)**2 -                    &
                  90*lam(6)*lam(7) + 117*lam(7)**2 +                            &
                  108*(4*lam(2)**2 - 8*lam(2)*lam(8) + lam(8)**2) +             &
                  66*lam(6)*lam(9) - 186*lam(7)*lam(9) + 118*lam(9)**2 +        &
                  6*(15*lam(6) - 63*lam(7) + 74*lam(9))*lam(10) +               &
                  240*lam(10)**2 -                                              &
                  6*lam(4)*(15*lam(5) + 3*lam(6) - 15*lam(7) +                  &
                     11*lam(9) + 15*lam(10)) +                                  &
                  6*lam(5)*(15*lam(6) - 39*lam(7) + 31*lam(9) +                 &
                     63*lam(10)))))) -                                          &
      3*qdB*(2*dRqsu*expuISs*qxy2**2*(1 + Rqsu)**2*tanhus**3*                   &
          (-2*expuISr*(1 + Rqru)*(3*lam(1) + lam(3))*                           &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +            &
            expdISr*(1 + Rqrd)*                                                 &
             (lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) + 51*lam(7) +          &
                  96*lam(9) - 52*lam(10)) +                                     &
               6*lam(1)*(3*lam(4) + 15*lam(5) -                                 &
                  3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10)))) -            &
         2*dRqsu*expuISs*qxy2*(1 + Rqsu)*tanhus*                                &
          (-2*expuISr*(qxy2 + qzt2)*(1 + Rqru)*(1 + Rqsu)*                      &
             (3*lam(1) + lam(3))*                                               &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +            &
            expuISr*m2*(-8*(3*lam(1) + lam(3))**2 +                             &
               9*(-lam(4) + lam(5) + lam(6) - lam(7))*lam(10) +                 &
               6*lam(10)**2) +                                                  &
            expdISr*m2*(72*lam(1)**2 - 168*lam(1)*lam(3) +                      &
               44*lam(3)**2 -                                                   &
               9*(3*lam(4) + 9*lam(5) - 3*lam(6) - 9*lam(7) -                   &
                  4*lam(9))*lam(9) +                                            &
               9*(lam(4) + 5*lam(5) - lam(6) - 5*lam(7) - 4*lam(9))*            &
                lam(10) + 12*lam(10)**2) +                                      &
            expdISr*(qxy2 + qzt2)*(1 + Rqrd)*(1 + Rqsu)*                        &
             (lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) + 51*lam(7) +          &
                  96*lam(9) - 52*lam(10)) +                                     &
               6*lam(1)*(3*lam(4) + 15*lam(5) -                                 &
                  3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10)))) +            &
         expuISs*quB*qxy2*tanhus**2*                                            &
          (dRqsu*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*                                &
             (-2*expuISr*(1 + Rqru)*(3*lam(1) + lam(3))*                        &
                (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +         &
               expdISr*(1 + Rqrd)*                                              &
                (lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) +                   &
                     51*lam(7) + 96*lam(9) - 52*lam(10)) +                      &
                  6*lam(1)*(3*lam(4) + 15*lam(5) -                              &
                     3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10)))) +         &
            3*(1 + Rqsu)*(dRqrd*expdISr*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*         &
                (2*lam(2) + lam(8))*                                            &
                (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                    &
                  10*lam(9) - 12*lam(10)) -                                     &
               2*dRqru*expuISr*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                  &
                (lam(2)*(3*lam(4) - 15*lam(5) - 3*lam(6) + 15*lam(7) -          &
                     20*lam(9) - 18*lam(10)) +                                  &
                  lam(8)*(-3*lam(4) + 6*lam(5) + 3*lam(6) - 6*lam(7) +          &
                     5*lam(9) + 15*lam(10))))) +                                &
         quB*(-3*dRqrd*expdISr*                                                 &
             (expuISs*(qxy2 + qzt2)*(1 + Rqsu)*                                 &
                (-1 + 2*qzt2*(1 + Rqrd)**2*sr)*(2*lam(2) + lam(8))*             &
                (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +                    &
                  10*lam(9) - 12*lam(10)) -                                     &
               2*expuISs*m2*qzt2*(1 + Rqrd)*sr*                                 &
                (6*(2*lam(2) + lam(8))**2 +                                     &
                  9*lam(9)*(lam(4) + lam(5) - lam(6) - lam(7) +                 &
                     2*lam(9)) +                                                &
                  (-3*lam(4) - 3*lam(5) + 3*(lam(6) + lam(7)) -                 &
                     10*lam(9))*lam(10) + 6*lam(10)**2) +                       &
               2*expdISs*m2*qzt2*(1 + Rqrd)*sr*                                 &
                (6*(4*lam(2)**2 - 8*lam(2)*lam(8) + lam(8)**2) -                &
                  3*lam(9)*(3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) +           &
                     lam(9)) +                                                  &
                  (-3*lam(4) + 3*(5*lam(5) + lam(6) - 5*lam(7)) +               &
                     26*lam(9))*lam(10) + 12*lam(10)**2) -                      &
               2*expdISs*(qxy2 + qzt2)*(1 + Rqsd)*                              &
                (-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                                 &
                (lam(2)*(3*lam(4) - 15*lam(5) - 3*lam(6) + 15*lam(7) -          &
                     20*lam(9) - 18*lam(10)) +                                  &
                  lam(8)*(-3*lam(4) + 6*lam(5) + 3*lam(6) - 6*lam(7) +          &
                     5*lam(9) + 15*lam(10)))) +                                 &
            3*dRqru*expuISr*(-2*expuISs*m2*qzt2*(1 + Rqru)*sr*                  &
                (6*(4*lam(2)**2 - 8*lam(2)*lam(8) + lam(8)**2) -                &
                  3*lam(9)*(3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) +           &
                     lam(9)) +                                                  &
                  (-3*lam(4) + 3*(5*lam(5) + lam(6) - 5*lam(7)) +               &
                     26*lam(9))*lam(10) + 12*lam(10)**2) +                      &
               2*expuISs*(qxy2 + qzt2)*(1 + Rqsu)*                              &
                (-1 + 2*qzt2*(1 + Rqru)**2*sr)*                                 &
                (lam(2)*(3*lam(4) - 15*lam(5) - 3*lam(6) + 15*lam(7) -          &
                     20*lam(9) - 18*lam(10)) +                                  &
                  lam(8)*(-3*lam(4) + 6*lam(5) + 3*lam(6) - 6*lam(7) +          &
                     5*lam(9) + 15*lam(10))) +                                  &
               expdISs*(-((qxy2 + qzt2)*(1 + Rqsd)*                             &
                     (-1 + 2*qzt2*(1 + Rqru)**2*sr)*(2*lam(2) + lam(8))*        &
                     (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +               &
                       10*lam(9) - 12*lam(10))) +                               &
                  2*m2*qzt2*(1 + Rqru)*sr*                                      &
                   (6*(2*lam(2) + lam(8))**2 +                                  &
                     9*lam(9)*                                                  &
                      (lam(4) + lam(5) - lam(6) - lam(7) + 2*lam(9)) +          &
                     (-3*lam(4) - 3*lam(5) + 3*(lam(6) + lam(7)) -              &
                        10*lam(9))*lam(10) + 6*lam(10)**2))) +                  &
            dRqsd*expdISs*(2*expdISr*                                           &
                ((qxy2 + qzt2)*(1 + Rqrd)*                                      &
                   (-1 + 2*qzt2*(1 + Rqsd)**2*ss)*(3*lam(1) + lam(3))*          &
                   (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10))          &
+ m2*qzt2*(1 + Rqsd)*ss*(8*(3*lam(1) + lam(3))**2 +                             &
                     9*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) -            &
                     6*lam(10)**2)) -                                           &
               expuISr*(2*m2*qzt2*(1 + Rqsd)*ss*                                &
                   (72*lam(1)**2 - 168*lam(1)*lam(3) + 44*lam(3)**2 -           &
                     9*(3*lam(4) + 9*lam(5) - 3*lam(6) - 9*lam(7) -             &
                        4*lam(9))*lam(9) +                                      &
                     9*(lam(4) + 5*lam(5) - lam(6) - 5*lam(7) -                 &
                        4*lam(9))*lam(10) + 12*lam(10)**2) +                    &
                  (qxy2 + qzt2)*(1 + Rqru)*                                     &
                   (-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                              &
                   (lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) +                &
                        51*lam(7) + 96*lam(9) - 52*lam(10)) +                   &
                     6*lam(1)*                                                  &
                      (3*lam(4) + 15*lam(5) -                                   &
                        3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10)))))       &
+ dRqsu*expuISs*(2*expuISr*((qxy2 + qzt2)*(1 + Rqru)*                           &
                   (-1 + 2*qzt2*(1 + Rqsu)**2*ss)*(3*lam(1) + lam(3))*          &
                   (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)) +        &
                  m2*qzt2*(1 + Rqsu)*ss*                                        &
                   (8*(3*lam(1) + lam(3))**2 +                                  &
                     9*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) -            &
                     6*lam(10)**2)) -                                           &
               expdISr*(2*m2*qzt2*(1 + Rqsu)*ss*                                &
                   (72*lam(1)**2 - 168*lam(1)*lam(3) + 44*lam(3)**2 -           &
                     9*(3*lam(4) + 9*lam(5) - 3*lam(6) - 9*lam(7) -             &
                        4*lam(9))*lam(9) +                                      &
                     9*(lam(4) + 5*lam(5) - lam(6) - 5*lam(7) -                 &
                        4*lam(9))*lam(10) + 12*lam(10)**2) +                    &
                  (qxy2 + qzt2)*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*      &
                   (lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) +                &
                        51*lam(7) + 96*lam(9) - 52*lam(10)) +                   &
                     6*lam(1)*                                                  &
                      (3*lam(4) + 15*lam(5) -                                   &
                        3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10))))))      &
) + expuISr*tanhur*(-(expdISs*qdB*quB*qzt2*tanhds*                              &
            (-9*dRqru*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                &
               (2*lam(2) + lam(8))*                                             &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -         &
                 12*lam(10)) +                                                  &
              dRqsd*m2*(1 + Rqsd)*ss*                                           &
               (432*lam(1)**2 - 1008*lam(1)*lam(3) + 264*lam(3)**2 +            &
                 9*lam(4)**2 +                                                  &
                 9*(13*lam(5)**2 + lam(6)**2 + 10*lam(6)*lam(7) +               &
                    13*lam(7)**2 - 2*lam(5)*(5*lam(6) + 13*lam(7))) +           &
                 90*(-3*lam(5) + lam(6) + 3*lam(7))*lam(9) +                    &
                 312*lam(9)**2 +                                                &
                 18*lam(4)*(5*lam(5) - lam(6) - 5*(lam(7) + lam(9))) +          &
                 6*(7*lam(4) + 23*lam(5) - 7*lam(6) - 23*lam(7) -               &
                    52*lam(9))*lam(10) + 64*lam(10)**2) +                       &
              dRqru*m2*(1 + Rqru)*sr*                                           &
               (9*lam(4)**2 + 9*lam(5)**2 + 9*lam(6)**2 +                       &
                 18*lam(6)*lam(7) + 9*lam(7)**2 +                               &
                 108*(2*lam(2) + lam(8))**2 - 6*lam(6)*lam(9) -                 &
                 6*lam(7)*lam(9) - 8*lam(9)**2 +                                &
                 6*lam(4)*(3*lam(5) - 3*lam(6) - 3*lam(7) + lam(9) -            &
                    9*lam(10)) +                                                &
                 18*(3*(lam(6) + lam(7)) - 10*lam(9))*lam(10) +                 &
                 108*lam(10)**2 -                                               &
                 6*lam(5)*(3*lam(6) + 3*lam(7) - lam(9) + 9*lam(10))) +         &
              3*dRqsd*(1 + Rqru)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                &
               (lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) + 51*lam(7) +        &
                    96*lam(9) - 52*lam(10)) +                                   &
                 6*lam(1)*(3*lam(4) + 15*lam(5) -                               &
                    3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10))))) -         &
         expdISs*qxy2*(1 + Rqsd)*tanhds**2*                                     &
          (18*dRqru*qdB*qxy2*(1 + Rqru)**2*(2*lam(2) + lam(8))*                 &
             (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -           &
               12*lam(10)) + dRqsd*quB*                                         &
             (m2*(24*(18*lam(1)**2 - 42*lam(1)*lam(3) + 11*lam(3)**2) +         &
                  9*lam(4)**2 +                                                 &
                  9*(13*lam(5)**2 + lam(6)**2 + 10*lam(6)*lam(7) +              &
                     13*lam(7)**2 - 2*lam(5)*(5*lam(6) + 13*lam(7))) +          &
                  90*(-3*lam(5) + lam(6) + 3*lam(7))*lam(9) +                   &
                  312*lam(9)**2 +                                               &
                  18*lam(4)*(5*lam(5) - lam(6) - 5*(lam(7) + lam(9))) +         &
                  6*(7*lam(4) + 23*lam(5) - 7*lam(6) - 23*lam(7) -              &
                     52*lam(9))*lam(10) + 64*lam(10)**2) +                      &
               6*qzt2*(1 + Rqru)*(1 + Rqsd)*                                    &
                (lam(3)*(-21*lam(4) - 51*lam(5) + 21*lam(6) +                   &
                     51*lam(7) + 96*lam(9) - 52*lam(10)) +                      &
                  6*lam(1)*(3*lam(4) + 15*lam(5) -                              &
                     3*(lam(6) + 5*lam(7) + 8*lam(9)) + 10*lam(10))))) +        &
         qdB*(dRqsu*expuISs*tanhus*                                             &
             (6*qzt2*(1 + Rqru)*                                                &
                (quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +                           &
                  2*qxy2*(1 + Rqsu)**2*tanhus)*(3*lam(1) + lam(3))*             &
                (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +         &
               m2*(1 + Rqsu)*(quB*qzt2*ss + qxy2*tanhus)*                       &
                (48*(3*lam(1) + lam(3))**2 +                                    &
                  9*(lam(4) - lam(5) - lam(6) + lam(7))**2 +                    &
                  6*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +               &
                  28*lam(10)**2)) +                                             &
            dRqru*(-18*expdISs*qxy2*(1 + Rqru)*                                 &
                (-((qxy2 + qzt2)*(1 + Rqru)*(1 + Rqsd)*                         &
                     (2*lam(2) + lam(8))*                                       &
                     (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +               &
                       10*lam(9) - 12*lam(10))) +                               &
                  m2*(6*(2*lam(2) + lam(8))**2 +                                &
                     9*lam(9)*                                                  &
                      (lam(4) + lam(5) - lam(6) - lam(7) + 2*lam(9)) +          &
                     (-3*lam(4) - 3*lam(5) + 3*(lam(6) + lam(7)) -              &
                        10*lam(9))*lam(10) + 6*lam(10)**2)) +                   &
               expuISs*(18*(1 + Rqsu)*                                          &
                   (-2*qxy2*(qxy2 + qzt2)*(1 + Rqru)**2 -                       &
                     quB*qzt2*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*tanhus +           &
                     2*qxy2**2*(1 + Rqru)**2*tanhus**2)*                        &
                   (lam(2)*(3*lam(4) - 15*lam(5) - 3*lam(6) +                   &
                        15*lam(7) - 20*lam(9) - 18*lam(10)) +                   &
                     lam(8)*(-3*lam(4) + 6*lam(5) + 3*lam(6) - 6*lam(7) +       &
                        5*lam(9) + 15*lam(10))) +                               &
                  m2*(1 + Rqru)*                                                &
                   (18*qxy2*(6*                                                 &
                         (4*lam(2)**2 - 8*lam(2)*lam(8) + lam(8)**2) -          &
                        3*lam(9)*                                               &
                         (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) +             &
                         lam(9)) +                                              &
                        (-3*lam(4) + 3*(5*lam(5) + lam(6) - 5*lam(7)) +         &
                         26*lam(9))*lam(10) + 12*lam(10)**2) +                  &
                     quB*qzt2*sr*tanhus*                                        &
                      (9*lam(4)**2 + 117*lam(5)**2 + 9*lam(6)**2 -              &
                        90*lam(6)*lam(7) + 117*lam(7)**2 +                      &
                        108*(4*lam(2)**2 - 8*lam(2)*lam(8) + lam(8)**2) +       &
                        66*lam(6)*lam(9) - 186*lam(7)*lam(9) +                  &
                        118*lam(9)**2 +                                         &
                        6*(15*lam(6) - 63*lam(7) + 74*lam(9))*lam(10) +         &
                        240*lam(10)**2 -                                        &
                        6*lam(4)*                                               &
                         (15*lam(5) + 3*lam(6) - 15*lam(7) + 11*lam(9) +        &
                         15*lam(10)) +                                          &
                        6*lam(5)*                                               &
                         (15*lam(6) - 39*lam(7) + 31*lam(9) + 63*lam(10))))     &
)))))/81.
  dtLAM(10)=                                                         &
   (12*dRqsd*expdISs*quB*qxy2**2*                                               &
       (2*expdISr*(1 + Rqrd) + expuISr*(1 + Rqru))*(1 + Rqsd)**2*tanhds**3*     &
       (3*lam(1) + lam(3))*(3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) +           &
         4*lam(10)) - 6*dRqsd*expdISs*quB*qxy2*(1 + Rqsd)*tanhds*               &
       (4*expdISr*(qxy2 + qzt2)*(1 + Rqrd)*(1 + Rqsd)*(3*lam(1) + lam(3))*      &
          (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +               &
         2*expuISr*(qxy2 + qzt2)*(1 + Rqru)*(1 + Rqsd)*                         &
          (3*lam(1) + lam(3))*                                                  &
          (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +               &
         2*expdISr*m2*(8*(3*lam(1) + lam(3))**2 +                               &
            9*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) - 6*lam(10)**2)       &
+ expuISr*m2*(8*(3*lam(1) + lam(3))**2 +                                        &
            9*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) - 6*lam(10)**2))      &
- 18*dRqrd*expdISr*quB*qxy2**2*(1 + Rqrd)**2*tanhdr**3*                         &
       (expuISs*(1 + Rqsu)*(-1 + tanhus**2)*(2*lam(2) + lam(8))*                &
          (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -              &
            12*lam(10)) + expdISs*(1 + Rqsd)*(-1 + tanhds**2)*                  &
          (lam(8)*(-3*lam(4) + 15*lam(5) + 3*lam(6) - 15*lam(7) +               &
               20*lam(9) + 18*lam(10)) +                                        &
            4*lam(2)*(3*lam(4) - 6*lam(5) - 3*lam(6) + 6*lam(7) -               &
               5*(lam(9) + 3*lam(10))))) -                                      &
      18*dRqru*expuISr*qdB*qxy2**2*(1 + Rqru)**2*tanhur**3*                     &
       (expdISs*(1 + Rqsd)*(-1 + tanhds**2)*(2*lam(2) + lam(8))*                &
          (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -              &
            12*lam(10)) + expuISs*(1 + Rqsu)*(-1 + tanhus**2)*                  &
          (lam(8)*(-3*lam(4) + 15*lam(5) + 3*lam(6) - 15*lam(7) +               &
               20*lam(9) + 18*lam(10)) +                                        &
            4*lam(2)*(3*lam(4) - 6*lam(5) - 3*lam(6) + 6*lam(7) -               &
               5*(lam(9) + 3*lam(10))))) +                                      &
      3*expdISs*qdB*quB*qxy2*tanhds**2*                                         &
       (2*dRqsd*(2*expdISr*(1 + Rqrd) + expuISr*(1 + Rqru))*                    &
          (-1 + 2*qzt2*(1 + Rqsd)**2*ss)*(3*lam(1) + lam(3))*                   &
          (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +               &
         3*(1 + Rqsd)*(dRqru*expuISr*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*            &
             (2*lam(2) + lam(8))*                                               &
             (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -           &
               12*lam(10)) + dRqrd*expdISr*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*      &
             (lam(8)*(-3*lam(4) + 15*lam(5) + 3*lam(6) - 15*lam(7) +            &
                  20*lam(9) + 18*lam(10)) +                                     &
               4*lam(2)*(3*lam(4) - 6*lam(5) - 3*lam(6) + 6*lam(7) -            &
                  5*(lam(9) + 3*lam(10)))))) +                                  &
      expdISr*tanhdr*(expuISs*qdB*quB*qzt2*tanhus*                              &
          (-9*dRqrd*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                  &
             (2*lam(2) + lam(8))*                                               &
             (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -           &
               12*lam(10)) + 6*dRqsu*(1 + Rqrd)*                                &
             (-1 + 2*qzt2*(1 + Rqsu)**2*ss)*(3*lam(1) + lam(3))*                &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +            &
            dRqsu*m2*(1 + Rqsu)*ss*                                             &
             (48*(3*lam(1) + lam(3))**2 +                                       &
               9*(lam(4) - lam(5) - lam(6) + lam(7))**2 +                       &
               6*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +                  &
               28*lam(10)**2) +                                                 &
            dRqrd*m2*(1 + Rqrd)*sr*                                             &
             (9*lam(4)**2 + 9*lam(5)**2 + 9*lam(6)**2 +                         &
               18*lam(6)*lam(7) + 9*lam(7)**2 +                                 &
               108*(2*lam(2) + lam(8))**2 - 60*lam(6)*lam(9) -                  &
               60*lam(7)*lam(9) + 118*lam(9)**2 +                               &
               6*lam(4)*(3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -           &
                  9*lam(10)) +                                                  &
               18*(3*(lam(6) + lam(7)) - 10*lam(9))*lam(10) +                   &
               108*lam(10)**2 -                                                 &
               6*lam(5)*(3*lam(6) + 3*lam(7) - 10*lam(9) + 9*lam(10)))) +       &
         expuISs*qxy2*(1 + Rqsu)*tanhus**2*                                     &
          (18*dRqrd*quB*qxy2*(1 + Rqrd)**2*(2*lam(2) + lam(8))*                 &
             (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -           &
               12*lam(10)) + dRqsu*qdB*                                         &
             (12*qzt2*(1 + Rqrd)*(1 + Rqsu)*(3*lam(1) + lam(3))*                &
                (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)) +           &
               m2*(48*(3*lam(1) + lam(3))**2 +                                  &
                  9*(lam(4) - lam(5) - lam(6) + lam(7))**2 +                    &
                  6*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +               &
                  28*lam(10)**2))) +                                            &
         2*expdISs*quB*qxy2*(1 + Rqsd)*tanhds**2*                               &
          (dRqsd*(12*qzt2*(1 + Rqrd)*(1 + Rqsd)*(3*lam(1) + lam(3))*            &
                (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)) +           &
               m2*(48*(3*lam(1) + lam(3))**2 +                                  &
                  9*(lam(4) - lam(5) - lam(6) + lam(7))**2 +                    &
                  6*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +               &
                  28*lam(10)**2)) +                                             &
            9*dRqrd*qxy2*(1 + Rqrd)**2*                                         &
             (lam(8)*(-3*lam(4) + 15*lam(5) + 3*lam(6) - 15*lam(7) +            &
                  20*lam(9) + 18*lam(10)) +                                     &
               4*lam(2)*(3*lam(4) - 6*lam(5) - 3*lam(6) + 6*lam(7) -            &
                  5*(lam(9) + 3*lam(10))))) +                                   &
         18*dRqrd*quB*qxy2*(1 + Rqrd)*                                          &
          (-(expuISs*(qxy2 + qzt2)*(1 + Rqrd)*(1 + Rqsu)*                       &
               (2*lam(2) + lam(8))*                                             &
               (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -         &
                 12*lam(10))) +                                                 &
            expuISs*m2*(6*(2*lam(2) + lam(8))**2 - 3*lam(9)**2 -                &
               10*lam(9)*lam(10) +                                              &
               3*lam(10)*(-lam(4) - lam(5) + lam(6) + lam(7) + 2*lam(10))       &
) + expdISs*m2*(12*(4*lam(2)**2 - 2*lam(2)*lam(8) + lam(8)**2) -                &
               6*lam(9)**2 + 58*lam(9)*lam(10) -                                &
               3*lam(10)*(5*lam(4) - 13*lam(5) - 5*lam(6) + 13*lam(7) +         &
                  8*lam(10))) -                                                 &
            expdISs*(qxy2 + qzt2)*(1 + Rqrd)*(1 + Rqsd)*                        &
             (lam(8)*(-3*lam(4) + 15*lam(5) + 3*lam(6) - 15*lam(7) +            &
                  20*lam(9) + 18*lam(10)) +                                     &
               4*lam(2)*(3*lam(4) - 6*lam(5) - 3*lam(6) + 6*lam(7) -            &
                  5*(lam(9) + 3*lam(10))))) +                                   &
         expdISs*qdB*quB*qzt2*tanhds*                                           &
          (12*dRqsd*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                  &
             (3*lam(1) + lam(3))*                                               &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +            &
            2*dRqsd*m2*(1 + Rqsd)*ss*                                           &
             (48*(3*lam(1) + lam(3))**2 +                                       &
               9*(lam(4) - lam(5) - lam(6) + lam(7))**2 +                       &
               6*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +                  &
               28*lam(10)**2) +                                                 &
            2*dRqrd*m2*(1 + Rqrd)*sr*                                           &
             (9*lam(4)**2 + 63*lam(5)**2 + 9*lam(6)**2 -                        &
               36*lam(6)*lam(7) + 63*lam(7)**2 +                                &
               108*(4*lam(2)**2 - 2*lam(2)*lam(8) + lam(8)**2) +                &
               30*lam(6)*lam(9) - 150*lam(7)*lam(9) + 118*lam(9)**2 +           &
               3*(15*lam(6) - 27*lam(7) + 2*lam(9))*lam(10) +                   &
               300*lam(10)**2 -                                                 &
               3*lam(4)*(12*lam(5) + 6*lam(6) - 12*lam(7) + 10*lam(9) +         &
                  15*lam(10)) +                                                 &
               3*lam(5)*(12*lam(6) - 42*lam(7) + 50*lam(9) + 27*lam(10)))       &
- 9*dRqrd*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                            &
             (lam(8)*(-3*lam(4) + 15*lam(5) + 3*lam(6) - 15*lam(7) +            &
                  20*lam(9) + 18*lam(10)) +                                     &
               4*lam(2)*(3*lam(4) - 6*lam(5) - 3*lam(6) + 6*lam(7) -            &
                  5*(lam(9) + 3*lam(10)))))) +                                  &
      expdISr*qxy2*tanhdr**2*(-24*dRqsd*expdISs*quB*qxy2*(1 + Rqrd)*            &
          (1 + Rqsd)**2*tanhds**3*(3*lam(1) + lam(3))*                          &
          (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) -               &
         12*dRqsu*expuISs*qdB*qxy2*(1 + Rqrd)*(1 + Rqsu)**2*tanhus**3*          &
          (3*lam(1) + lam(3))*                                                  &
          (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) -               &
         3*expuISs*qdB*quB*tanhus**2*                                           &
          (3*dRqrd*(1 + Rqsu)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                   &
             (2*lam(2) + lam(8))*                                               &
             (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -           &
               12*lam(10)) + 2*dRqsu*(1 + Rqrd)*                                &
             (-1 + 2*qzt2*(1 + Rqsu)**2*ss)*(3*lam(1) + lam(3))*                &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))) +           &
         2*expdISs*quB*(1 + Rqrd)*tanhds*                                       &
          (12*dRqsd*qxy2*(1 + Rqsd)**2*(3*lam(1) + lam(3))*                     &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +            &
            dRqrd*m2*(9*lam(4)**2 + 63*lam(5)**2 + 9*lam(6)**2 -                &
               36*lam(6)*lam(7) + 63*lam(7)**2 +                                &
               108*(4*lam(2)**2 - 2*lam(2)*lam(8) + lam(8)**2) +                &
               30*lam(6)*lam(9) - 150*lam(7)*lam(9) + 118*lam(9)**2 +           &
               3*(15*lam(6) - 27*lam(7) + 2*lam(9))*lam(10) +                   &
               300*lam(10)**2 -                                                 &
               3*lam(4)*(12*lam(5) + 6*lam(6) - 12*lam(7) + 10*lam(9) +         &
                  15*lam(10)) +                                                 &
               3*lam(5)*(12*lam(6) - 42*lam(7) + 50*lam(9) + 27*lam(10)))       &
- 9*dRqrd*qzt2*(1 + Rqrd)*(1 + Rqsd)*                                           &
             (lam(8)*(-3*lam(4) + 15*lam(5) + 3*lam(6) - 15*lam(7) +            &
                  20*lam(9) + 18*lam(10)) +                                     &
               4*lam(2)*(3*lam(4) - 6*lam(5) - 3*lam(6) + 6*lam(7) -            &
                  5*(lam(9) + 3*lam(10))))) -                                   &
         3*expdISs*qdB*quB*tanhds**2*                                           &
          (4*dRqsd*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*                   &
             (3*lam(1) + lam(3))*                                               &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +            &
            3*dRqrd*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                  &
             (lam(8)*(-3*lam(4) + 15*lam(5) + 3*lam(6) - 15*lam(7) +            &
                  20*lam(9) + 18*lam(10)) +                                     &
               4*lam(2)*(3*lam(4) - 6*lam(5) - 3*lam(6) + 6*lam(7) -            &
                  5*(lam(9) + 3*lam(10))))) +                                   &
         expuISs*(1 + Rqrd)*tanhus*                                             &
          (12*dRqsu*qdB*qxy2*(1 + Rqsu)**2*(3*lam(1) + lam(3))*                 &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +            &
            dRqrd*quB*(-18*qzt2*(1 + Rqrd)*(1 + Rqsu)*                          &
                (2*lam(2) + lam(8))*                                            &
                (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -        &
                  12*lam(10)) +                                                 &
               m2*(9*lam(4)**2 + 9*lam(5)**2 + 9*lam(6)**2 +                    &
                  18*lam(6)*lam(7) + 9*lam(7)**2 +                              &
                  108*(2*lam(2) + lam(8))**2 - 60*lam(6)*lam(9) -               &
                  60*lam(7)*lam(9) + 118*lam(9)**2 +                            &
                  6*lam(4)*(3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -        &
                     9*lam(10)) +                                               &
                  18*(3*(lam(6) + lam(7)) - 10*lam(9))*lam(10) +                &
                  108*lam(10)**2 -                                              &
                  6*lam(5)*(3*lam(6) + 3*lam(7) - 10*lam(9) + 9*lam(10))))      &
) + 3*qdB*quB*(4*dRqsd*expdISs*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsd)**2*ss)*       &
             (3*lam(1) + lam(3))*                                               &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +            &
            2*dRqsu*expuISs*(1 + Rqrd)*(-1 + 2*qzt2*(1 + Rqsu)**2*ss)*          &
             (3*lam(1) + lam(3))*                                               &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +            &
            3*dRqrd*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                             &
             (expuISs*(1 + Rqsu)*(2*lam(2) + lam(8))*                           &
                (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -        &
                  12*lam(10)) +                                                 &
               expdISs*(1 + Rqsd)*                                              &
                (lam(8)*(-3*lam(4) + 15*lam(5) + 3*lam(6) - 15*lam(7) +         &
                     20*lam(9) + 18*lam(10)) +                                  &
                  4*lam(2)*(3*lam(4) - 6*lam(5) - 3*lam(6) + 6*lam(7) -         &
                     5*(lam(9) + 3*lam(10))))))) +                              &
      expuISr*qxy2*tanhur**2*(-12*dRqsd*expdISs*quB*qxy2*(1 + Rqru)*            &
          (1 + Rqsd)**2*tanhds**3*(3*lam(1) + lam(3))*                          &
          (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) -               &
         3*expdISs*qdB*quB*tanhds**2*                                           &
          (3*dRqru*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                   &
             (2*lam(2) + lam(8))*                                               &
             (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -           &
               12*lam(10)) + 2*dRqsd*(1 + Rqru)*                                &
             (-1 + 2*qzt2*(1 + Rqsd)**2*ss)*(3*lam(1) + lam(3))*                &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10))) +           &
         expdISs*(1 + Rqru)*tanhds*                                             &
          (12*dRqsd*quB*qxy2*(1 + Rqsd)**2*(3*lam(1) + lam(3))*                 &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +            &
            dRqru*qdB*(-18*qzt2*(1 + Rqru)*(1 + Rqsd)*                          &
                (2*lam(2) + lam(8))*                                            &
                (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -        &
                  12*lam(10)) +                                                 &
               m2*(9*lam(4)**2 + 9*lam(5)**2 + 9*lam(6)**2 +                    &
                  18*lam(6)*lam(7) + 9*lam(7)**2 +                              &
                  108*(2*lam(2) + lam(8))**2 - 60*lam(6)*lam(9) -               &
                  60*lam(7)*lam(9) + 118*lam(9)**2 +                            &
                  6*lam(4)*(3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -        &
                     9*lam(10)) +                                               &
                  18*(3*(lam(6) + lam(7)) - 10*lam(9))*lam(10) +                &
                  108*lam(10)**2 -                                              &
                  6*lam(5)*(3*lam(6) + 3*lam(7) - 10*lam(9) + 9*lam(10))))      &
) + qdB*(6*(1 + Rqru)*(dRqsd*expdISs*quB*                                       &
                (-1 + 2*qzt2*(1 + Rqsd)**2*ss) -                                &
               2*dRqsu*expuISs*                                                 &
                (quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +                           &
                  2*qxy2*(1 + Rqsu)**2*tanhus)*(-1 + tanhus**2))*               &
             (3*lam(1) + lam(3))*                                               &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +            &
            dRqru*(9*expdISs*quB*(1 + Rqsd)*                                    &
                (-1 + 2*qzt2*(1 + Rqru)**2*sr)*(2*lam(2) + lam(8))*             &
                (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -        &
                  12*lam(10)) +                                                 &
               2*expuISs*m2*(1 + Rqru)*tanhus*                                  &
                (9*lam(4)**2 + 63*lam(5)**2 + 9*lam(6)**2 -                     &
                  36*lam(6)*lam(7) + 63*lam(7)**2 +                             &
                  108*(4*lam(2)**2 - 2*lam(2)*lam(8) + lam(8)**2) +             &
                  30*lam(6)*lam(9) - 150*lam(7)*lam(9) + 118*lam(9)**2 +        &
                  3*(15*lam(6) - 27*lam(7) + 2*lam(9))*lam(10) +                &
                  300*lam(10)**2 -                                              &
                  3*lam(4)*(12*lam(5) + 6*lam(6) - 12*lam(7) +                  &
                     10*lam(9) + 15*lam(10)) +                                  &
                  3*lam(5)*(12*lam(6) - 42*lam(7) + 50*lam(9) +                 &
                     27*lam(10))) -                                             &
               9*expuISs*(1 + Rqsu)*                                            &
                (2*qzt2*(1 + Rqru)**2*tanhus +                                  &
                  quB*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*(-1 + tanhus**2))*         &
                (lam(8)*(-3*lam(4) + 15*lam(5) + 3*lam(6) - 15*lam(7) +         &
                     20*lam(9) + 18*lam(10)) +                                  &
                  4*lam(2)*(3*lam(4) - 6*lam(5) - 3*lam(6) + 6*lam(7) -         &
                     5*(lam(9) + 3*lam(10))))))) -                              &
      3*qdB*(-4*dRqsu*expuISs*qxy2**2*                                          &
          (expdISr*(1 + Rqrd) + 2*expuISr*(1 + Rqru))*(1 + Rqsu)**2*            &
          tanhus**3*(3*lam(1) + lam(3))*                                        &
          (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +               &
         2*dRqsu*expuISs*qxy2*(1 + Rqsu)*tanhus*                                &
          (2*expdISr*(qxy2 + qzt2)*(1 + Rqrd)*(1 + Rqsu)*                       &
             (3*lam(1) + lam(3))*                                               &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +            &
            4*expuISr*(qxy2 + qzt2)*(1 + Rqru)*(1 + Rqsu)*                      &
             (3*lam(1) + lam(3))*                                               &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +            &
            expdISr*m2*(8*(3*lam(1) + lam(3))**2 +                              &
               9*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) -                  &
               6*lam(10)**2) +                                                  &
            2*expuISr*m2*(8*(3*lam(1) + lam(3))**2 +                            &
               9*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) -                  &
               6*lam(10)**2)) -                                                 &
         expuISs*quB*qxy2*tanhus**2*                                            &
          (2*dRqsu*(expdISr*(1 + Rqrd) + 2*expuISr*(1 + Rqru))*                 &
             (-1 + 2*qzt2*(1 + Rqsu)**2*ss)*(3*lam(1) + lam(3))*                &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +            &
            3*(1 + Rqsu)*(dRqrd*expdISr*(-1 + 2*qzt2*(1 + Rqrd)**2*sr)*         &
                (2*lam(2) + lam(8))*                                            &
                (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -        &
                  12*lam(10)) +                                                 &
               dRqru*expuISr*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                    &
                (lam(8)*(-3*lam(4) + 15*lam(5) + 3*lam(6) - 15*lam(7) +         &
                     20*lam(9) + 18*lam(10)) +                                  &
                  4*lam(2)*(3*lam(4) - 6*lam(5) - 3*lam(6) + 6*lam(7) -         &
                     5*(lam(9) + 3*lam(10)))))) +                               &
         quB*(2*dRqsd*expdISs*                                                  &
             (2*expdISr*(qxy2 + qzt2)*(1 + Rqrd)*                               &
                (-1 + 2*qzt2*(1 + Rqsd)**2*ss)*(3*lam(1) + lam(3))*             &
                (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +         &
               expuISr*(qxy2 + qzt2)*(1 + Rqru)*                                &
                (-1 + 2*qzt2*(1 + Rqsd)**2*ss)*(3*lam(1) + lam(3))*             &
                (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +         &
               2*expdISr*m2*qzt2*(1 + Rqsd)*ss*                                 &
                (8*(3*lam(1) + lam(3))**2 +                                     &
                  9*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) -               &
                  6*lam(10)**2) +                                               &
               expuISr*m2*qzt2*(1 + Rqsd)*ss*                                   &
                (8*(3*lam(1) + lam(3))**2 +                                     &
                  9*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) -               &
                  6*lam(10)**2)) +                                              &
            2*dRqsu*expuISs*(expdISr*(qxy2 + qzt2)*(1 + Rqrd)*                  &
                (-1 + 2*qzt2*(1 + Rqsu)**2*ss)*(3*lam(1) + lam(3))*             &
                (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +         &
               2*expuISr*(qxy2 + qzt2)*(1 + Rqru)*                              &
                (-1 + 2*qzt2*(1 + Rqsu)**2*ss)*(3*lam(1) + lam(3))*             &
                (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +         &
               expdISr*m2*qzt2*(1 + Rqsu)*ss*                                   &
                (8*(3*lam(1) + lam(3))**2 +                                     &
                  9*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) -               &
                  6*lam(10)**2) +                                               &
               2*expuISr*m2*qzt2*(1 + Rqsu)*ss*                                 &
                (8*(3*lam(1) + lam(3))**2 +                                     &
                  9*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) -               &
                  6*lam(10)**2)) -                                              &
            3*(dRqrd*expdISr*(-(expuISs*(qxy2 + qzt2)*(1 + Rqsu)*               &
                     (-1 + 2*qzt2*(1 + Rqrd)**2*sr)*(2*lam(2) + lam(8))*        &
                     (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +               &
                       10*lam(9) - 12*lam(10))) +                               &
                  2*expuISs*m2*qzt2*(1 + Rqrd)*sr*                              &
                   (6*(2*lam(2) + lam(8))**2 - 3*lam(9)**2 -                    &
                     10*lam(9)*lam(10) +                                        &
                     3*lam(10)*                                                 &
                      (-lam(4) - lam(5) + lam(6) + lam(7) + 2*lam(10))) +       &
                  2*expdISs*m2*qzt2*(1 + Rqrd)*sr*                              &
                   (12*(4*lam(2)**2 - 2*lam(2)*lam(8) + lam(8)**2) -            &
                     6*lam(9)**2 + 58*lam(9)*lam(10) -                          &
                     3*lam(10)*                                                 &
                      (5*lam(4) - 13*lam(5) - 5*lam(6) + 13*lam(7) +            &
                        8*lam(10))) -                                           &
                  expdISs*(qxy2 + qzt2)*(1 + Rqsd)*                             &
                   (-1 + 2*qzt2*(1 + Rqrd)**2*sr)*                              &
                   (lam(8)*(-3*lam(4) + 15*lam(5) + 3*lam(6) -                  &
                        15*lam(7) + 20*lam(9) + 18*lam(10)) +                   &
                     4*lam(2)*                                                  &
                      (3*lam(4) - 6*lam(5) - 3*lam(6) + 6*lam(7) -              &
                        5*(lam(9) + 3*lam(10))))) +                             &
               dRqru*expuISr*(-(expdISs*(qxy2 + qzt2)*(1 + Rqsd)*               &
                     (-1 + 2*qzt2*(1 + Rqru)**2*sr)*(2*lam(2) + lam(8))*        &
                     (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +               &
                       10*lam(9) - 12*lam(10))) +                               &
                  2*expdISs*m2*qzt2*(1 + Rqru)*sr*                              &
                   (6*(2*lam(2) + lam(8))**2 - 3*lam(9)**2 -                    &
                     10*lam(9)*lam(10) +                                        &
                     3*lam(10)*                                                 &
                      (-lam(4) - lam(5) + lam(6) + lam(7) + 2*lam(10))) +       &
                  2*expuISs*m2*qzt2*(1 + Rqru)*sr*                              &
                   (12*(4*lam(2)**2 - 2*lam(2)*lam(8) + lam(8)**2) -            &
                     6*lam(9)**2 + 58*lam(9)*lam(10) -                          &
                     3*lam(10)*                                                 &
                      (5*lam(4) - 13*lam(5) - 5*lam(6) + 13*lam(7) +            &
                        8*lam(10))) -                                           &
                  expuISs*(qxy2 + qzt2)*(1 + Rqsu)*                             &
                   (-1 + 2*qzt2*(1 + Rqru)**2*sr)*                              &
                   (lam(8)*(-3*lam(4) + 15*lam(5) + 3*lam(6) -                  &
                        15*lam(7) + 20*lam(9) + 18*lam(10)) +                   &
                     4*lam(2)*(3*lam(4) - 6*lam(5) - 3*lam(6) + 6*lam(7) -      &
                        5*(lam(9) + 3*lam(10)))))))) +                          &
      expuISr*tanhur*(expdISs*qdB*quB*qzt2*tanhds*                              &
          (-9*dRqru*(1 + Rqsd)*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*                  &
             (2*lam(2) + lam(8))*                                               &
             (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -           &
               12*lam(10)) + 6*dRqsd*(1 + Rqru)*                                &
             (-1 + 2*qzt2*(1 + Rqsd)**2*ss)*(3*lam(1) + lam(3))*                &
             (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +            &
            dRqsd*m2*(1 + Rqsd)*ss*                                             &
             (48*(3*lam(1) + lam(3))**2 +                                       &
               9*(lam(4) - lam(5) - lam(6) + lam(7))**2 +                       &
               6*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +                  &
               28*lam(10)**2) +                                                 &
            dRqru*m2*(1 + Rqru)*sr*                                             &
             (9*lam(4)**2 + 9*lam(5)**2 + 9*lam(6)**2 + 18*lam(6)*lam(7) +      &
               9*lam(7)**2 + 108*(2*lam(2) + lam(8))**2 -                       &
               60*lam(6)*lam(9) - 60*lam(7)*lam(9) + 118*lam(9)**2 +            &
               6*lam(4)*(3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -           &
                  9*lam(10)) +                                                  &
               18*(3*(lam(6) + lam(7)) - 10*lam(9))*lam(10) +                   &
               108*lam(10)**2 -                                                 &
               6*lam(5)*(3*lam(6) + 3*lam(7) - 10*lam(9) + 9*lam(10)))) +       &
         expdISs*qxy2*(1 + Rqsd)*tanhds**2*                                     &
          (18*dRqru*qdB*qxy2*(1 + Rqru)**2*(2*lam(2) + lam(8))*                 &
             (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) + 10*lam(9) -           &
               12*lam(10)) + dRqsd*quB*                                         &
             (12*qzt2*(1 + Rqru)*(1 + Rqsd)*(3*lam(1) + lam(3))*                &
                (3*(lam(4) - lam(5) - lam(6) + lam(7)) + 4*lam(10)) +           &
               m2*(48*(3*lam(1) + lam(3))**2 +                                  &
                  9*(lam(4) - lam(5) - lam(6) + lam(7))**2 +                    &
                  6*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +               &
                  28*lam(10)**2))) +                                            &
         qdB*(2*dRqsu*expuISs*tanhus*                                           &
             (6*qzt2*(1 + Rqru)*                                                &
                (quB*(-1 + 2*qzt2*(1 + Rqsu)**2*ss) +                           &
                  2*qxy2*(1 + Rqsu)**2*tanhus)*(3*lam(1) + lam(3))*             &
                (3*lam(4) - 3*(lam(5) + lam(6) - lam(7)) + 4*lam(10)) +         &
               m2*(1 + Rqsu)*(quB*qzt2*ss + qxy2*tanhus)*                       &
                (48*(3*lam(1) + lam(3))**2 +                                    &
                  9*(lam(4) - lam(5) - lam(6) + lam(7))**2 +                    &
                  6*(lam(4) - lam(5) - lam(6) + lam(7))*lam(10) +               &
                  28*lam(10)**2)) +                                             &
            dRqru*(18*expdISs*qxy2*(1 + Rqru)*                                  &
                (-((qxy2 + qzt2)*(1 + Rqru)*(1 + Rqsd)*                         &
                     (2*lam(2) + lam(8))*                                       &
                     (3*lam(4) + 3*lam(5) - 3*lam(6) - 3*lam(7) +               &
                       10*lam(9) - 12*lam(10))) +                               &
                  m2*(6*(2*lam(2) + lam(8))**2 - 3*lam(9)**2 -                  &
                     10*lam(9)*lam(10) +                                        &
                     3*lam(10)*                                                 &
                      (-lam(4) - lam(5) + lam(6) + lam(7) + 2*lam(10)))) +      &
               expuISs*(9*(1 + Rqsu)*                                           &
                   (-2*qxy2*(qxy2 + qzt2)*(1 + Rqru)**2 -                       &
                     quB*qzt2*(-1 + 2*qzt2*(1 + Rqru)**2*sr)*tanhus +           &
                     2*qxy2**2*(1 + Rqru)**2*tanhus**2)*                        &
                   (lam(8)*(-3*lam(4) + 15*lam(5) + 3*lam(6) -                  &
                        15*lam(7) + 20*lam(9) + 18*lam(10)) +                   &
                     4*lam(2)*(3*lam(4) - 6*lam(5) - 3*lam(6) + 6*lam(7) -      &
                        5*(lam(9) + 3*lam(10)))) +                              &
                  2*m2*(1 + Rqru)*                                              &
                   (9*qxy2*(12*                                                 &
                         (4*lam(2)**2 - 2*lam(2)*lam(8) + lam(8)**2) -          &
                        6*lam(9)**2 + 58*lam(9)*lam(10) -                       &
                        3*lam(10)*                                              &
                         (5*lam(4) - 13*lam(5) - 5*lam(6) + 13*lam(7) +         &
                          8*lam(10))) +                                         &
                     quB*qzt2*sr*tanhus*                                        &
                      (9*lam(4)**2 + 63*lam(5)**2 + 9*lam(6)**2 -               &
                        36*lam(6)*lam(7) + 63*lam(7)**2 +                       &
                        108*(4*lam(2)**2 - 2*lam(2)*lam(8) + lam(8)**2) +       &
                        30*lam(6)*lam(9) - 150*lam(7)*lam(9) +                  &
                        118*lam(9)**2 +                                         &
                        3*(15*lam(6) - 27*lam(7) + 2*lam(9))*lam(10) +          &
                        300*lam(10)**2 -                                        &
                        3*lam(4)*                                               &
                         (12*lam(5) + 6*lam(6) - 12*lam(7) + 10*lam(9) +        &
                          15*lam(10)) +                                         &
                        3*lam(5)*                                               &
          (12*lam(6) - 42*lam(7) + 50*lam(9) + 27*lam(10)))))))))/81.

     dtLAM=dtLAM/(quB*qdB)
  END SUBROUTINE DLAM_CAL
END MODULE LAM_MOD
