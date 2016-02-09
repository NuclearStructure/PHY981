!     Minnesota NN potential in k-space 
!     Units are such that v00*sqrt(wt(i)*wt(j))*p*ps (ditto for v02 etc.)
!     is in MeV, where p,ps are momentum points in 1/fm and wt(i),wt(j) are
!     the corresponding weights.  

      MODULE Minnforce 
        IMPLICIT NONE
        PRIVATE
        PUBLIC:: Minnp1p2 
        
        contains
      

        SUBROUTINE Minnp1p2(ll,S,T,coup,p,ps,v00,v02,v20,v22)
        implicit none 
        real(8):: p,ps, mu(3),VC(3),WC(3),VS(3),VCtild(3), fkkw
        real(8):: pi
        real(8) :: v00, v02, v20, v22
        integer:: ll, S, T,i
        logical :: coup
        real(8):: V_i(3),mu_i(3),W_i(3), M_i(3), B_i(3),H_i(3)
        real(8) :: pref_i,plusST, minusST, aplus, aminus
        data V_i/200.,-178.,-91.85/
        data mu_i/1.487,0.639,0.465/
        data W_i/.5,.25,.25/
        data M_i/.5,.25,.25/
        data B_i/0.,.25,-.25/
        data H_i/0.,.25,-.25/
        intent(in) :: coup, ll, p, ps, S,T
        intent(out) :: v00, v02, v20,v22


        v02 = 0.0
        v20 = 0.0
        v00 = 0.0
        v22 = 0.0

        pi = acos(-1.0)
        IF((-1)**(ll+S+T) .GT. 0)THEN
             WRITE(*,*)'INVALID LST', ll,S,T
             CALL ABORT
        ENDIF

        DO i = 1, 3
          VC(i) = V_i(i)*(W_i(i)+ .5d0*(B_i(i)-H_i(i)))
          VS(i) = V_i(i)*B_i(i)/2.d0
          WC(i) = -V_i(i)*H_i(i)/2.d0
          VCtild(i) = V_i(i)*M_i(i)
        ENDDO 

        IF(ll.gt.8)return

    
      IF(.not. coup)THEN
           v00=0.d0
           DO i = 1, 3
                  pref_i = 2./pi*dexp(-(p**2+ps**2)/4/mu_i(i))*(pi/mu_i(i))**1.5/(8.*pi)
                  plusST = VC(i) + WC(i)*(2.*T*(T+1.)-3.) + VS(i)*(2.*S*(S+1.)-3.)
                  minusST = VCtild(i)
                  aplus = p*ps/2./mu_i(i)
                  aminus = -p*ps/2./mu_i(i) 
                  v00 = v00 + pref_i*(plusST*xlegendint(ll,aplus) + minusST*xlegendint(ll,aminus))
!                  IF(ll.gt.5) v00 = fkkw*(pref_a*xlegendint(ll,alpha) + pref_b*xlegendint(ll,beta))
!                  IF(ll.le.5) v00 = fkkw*(pref_a*flalpha(ll, alpha) 
!     $                            + pref_b*flalpha(ll, beta))
          ENDDO
          RETURN
        ELSE
           v00=0.d0;v22=0.d0
           DO i = 1, 3
                  pref_i = 2./pi*dexp(-(p**2+ps**2)/4/mu_i(i))*(pi/mu_i(i))**1.5/(8.*pi)
                  plusST = VC(i) + WC(i)*(2.*T*(T+1.)-3.) + VS(i)*(2.*S*(S+1.)-3.)
                  minusST = VCtild(i)
                  aplus = p*ps/2./mu_i(i)
                  aminus = -p*ps/2./mu_i(i) 
                  v00 = v00 + pref_i*(plusST*xlegendint(ll,aplus) + minusST*xlegendint(ll,aminus))
                  v22 = v22 + pref_i*(plusST*xlegendint(ll+2,aplus) + minusST*xlegendint(ll+2,aminus))
                  
           ENDDO
           RETURN
        ENDIF
        end subroutine
	


        real(8) function xlegendint(ll, alpha)
          implicit none
          integer :: ll, i, nx
          real(8) :: alpha, xm(200), xwt(200), vec(15),pl(200)

           nx = 30
           call gauleg(-1.d0, 1.d0, xm(1:nx), xwt(1:nx),nx)

           do i = 1, nx
              call fleg(xm(i), vec,10)
              pl(i) = vec(ll+1)
           enddo


           xlegendint = 0.

           do i = 1, nx

            xlegendint = xlegendint + exp(xm(i)*alpha)*pl(i)*xwt(i)

           enddo
        end function

      SUBROUTINE fleg(x,pl,nl)
      implicit none
      INTEGER nl
      REAL(8) x,pl(nl)
      INTEGER j
      REAL(8) d,f1,f2,twox
      pl(1)=1.
      pl(2)=x
      if(nl.gt.2) then
        twox=2.*x
        f2=x
        d=1.
        do 11 j=3,nl
          f1=d
          f2=f2+twox
          d=d+1.
          pl(j)=(f2*pl(j-1)-f1*pl(j-2))/d
11      continue
      endif
      return
      END subroutine

      SUBROUTINE gauleg(x1,x2,x,w,n)
      INTEGER n
      REAL*8 x1,x2,x(n),w(n)
      REAL*8 EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      REAL*8 p1,p2,p3,pp,xl,xm,z,z1
      
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END SUBROUTINE GAULEG

      END MODULE Minnforce 
