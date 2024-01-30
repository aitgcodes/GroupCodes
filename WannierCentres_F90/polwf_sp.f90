MODULE runparam

 IMPLICIT NONE
 INTEGER :: tbeg, nt, nsp, nst, nskp, corrs, ncmax, crdstyle, inpunits, epq
 LOGICAL :: rmode, wfadj, tspol
 INTEGER, ALLOCATABLE :: na(:)
 CHARACTER(LEN=40) :: crdfile, wffile
 CHARACTER(LEN=2), ALLOCATABLE :: sym(:)
 DOUBLE PRECISION,ALLOCATABLE :: zv(:)
 DOUBLE PRECISION,PARAMETER :: rau = 0.529177d0, eau = 1.6022d0, debye = eau*rau*10.d0/3.33564d0
 DOUBLE PRECISION :: rloc

END MODULE runparam

MODULE box

 INTEGER :: natm
 DOUBLE PRECISION, ALLOCATABLE :: s(:,:,:,:), h(:,:,:), r(:,:,:,:)
 DOUBLE PRECISION, ALLOCATABLE :: wf(:,:,:), sprd(:,:)
 DOUBLE PRECISION, ALLOCATABLE ::  vol(:)
 INTEGER, ALLOCATABLE :: nep(:,:), wftg(:) 
 DOUBLE PRECISION :: href(3,3)

END MODULE box

PROGRAM polarizn

USE box
USE runparam

IMPLICIT NONE
INTEGER :: i, j, it, nax, ntot, tot_ep
DOUBLE PRECISION,ALLOCATABLE :: pol(:,:)
DOUBLE PRECISION :: dmion(3), dmel(3), z, dmconv, avpol(3), avdmion(3), avdmel(3)
DOUBLE PRECISION :: pmag, polq(3), f, f0, avvol, ufact, volm

OPEN(UNIT=1,FILE='polwf.in')

 READ(1,*) crdfile
 READ(1,*) rmode, crdstyle, inpunits
 READ(1,*) nst, tspol
 READ(1,*) wffile
 READ(1,*) tbeg, nt, nskp, nsp
 ALLOCATE(na(nsp),zv(nsp),sym(nsp))

 READ(1,*) na(:)
 READ(1,*) sym(:)
 READ(1,*) zv(:)
 READ(1,*) rloc
 READ(1,*) corrs, ncmax
 READ(1,*) wfadj
 DO i = 1,3
  READ(1,*) (href(i,j),j=1,3)
  WRITE(6,*) (href(i,j),j=1,3)
 END DO

CLOSE(UNIT=1)

natm = 0
nax = 0
ntot = nt/nskp
WRITE(6,*) 'ntot', ntot
IF (inpunits == 0) THEN
 ufact = 1.0             !******* Atomic units on input
ELSE 
 ufact = rau             !******* Angstroms on input
END IF

epq = 2
IF (tspol) epq = 1       !****** Set electron count quantum to 1 if
                         !****** spin polarized wannier centres

DO i = 1,nsp
 natm = natm + na(i)
 nax = MAX(nax,na(i))
END DO

ALLOCATE(h(nt,3,3), s(nt,nsp,nax,3), r(nt,nsp,nax,3), wf(nt,nst,3), sprd(nt,nst), vol(nt), STAT=i)
IF (i /= 0) THEN
 WRITE(6,*) 'ALLOCATION FAULT IN MAIN (1)!'
 STOP
END IF

CALL readpos !***** Atoms periodized into -L/2 to L/2 box
CALL readwf  !***** WFC read in as periodic into -L/2 to L/2 box

ALLOCATE(pol(nt,3),STAT=i)
IF (i /= 0) THEN
 WRITE(6,*) 'ALLOCATION FAULT IN MAIN (2)!'
  STOP
END IF

!**** COUNTING ELECRON PAIRS AND MOVING WFC IF NEEDED *******************

ALLOCATE(nep(nsp,nax), wftg(nst))
nep = 0
wftg = 0
IF (wfadj) THEN
 CALL epcount(1,0,0)

 DO i = 1,nst
!  WRITE(6,*) "State of WF ", i, "is :", wftg(i)
  IF (wftg(i) == 0) THEN
!   WRITE(6,*) i
   CALL epcount(1,1,i)
   IF (wftg(i) == 0) WRITE(6,*) "State ",i," not assigned!!!"
  END IF
 END DO
END IF

tot_ep = 0
WRITE(6,*) 'AVERAGE ELECTRON COUNT FOR EACH SPECIES :'
!OPEN(UNIT=88,FILE='at_cnt.dat')
DO i = 1,nsp
  WRITE(6,*) 'nsp = ',i,' nel = ',SUM(nep(i,:))/FLOAT(na(i)),' Zion =',zv(i)-SUM(nep(i,:))/FLOAT(na(i))
!  WRITE(88,"(2(I5))") (nep(i,j),j=1,na(i))
  tot_ep = tot_ep + SUM(nep(i,:))
END DO

IF (tot_ep /= nst) WRITE(6,*) "WARNING!! ",nst-tot_ep," states not assigned!" 

dmconv = eau*10.0d0/(rau**2)    !*** Conversion of Polzn. from a.u. to C/m2
WRITE(6,*) 'POLARIZATION QUANTUM (A.U. UNITS):'
DO i = 1,3
 polq(:) = h(1,i,:)/vol(1) 
 WRITE(6,*) polq(:)
END DO
WRITE(6,*) 'POLARIZATION QUANTUM (SI UNITS):'
DO i = 1,3
 polq(:) = h(1,i,:)*dmconv/vol(1) 
 WRITE(6,*) polq(:)
END DO

OPEN(UNIT=11,FILE='wfc.xyz')
DO it = 1,nt
 WRITE(11,*) nst
 WRITE(11,*) 
 DO i = 1,nst
  WRITE(11,"(A2,3(F12.8,2X))") 'H',wf(it,i,:)
 END DO 
END DO
CLOSE(UNIT=11)

!**** CALCULATION OF THE POLARIZATION ***********************************
WRITE(0,*) 'Calculating polarization ...'

avpol = 0.0d0
avdmion = 0.0d0
avdmel = 0.0d0
avvol = 0.0d0

OPEN(UNIT=3,FILE='polwf.out')
DO it = 1,nt

 dmion = 0.0d0
 dmel = 0.0d0
 
 DO i = 1,nsp
  z = zv(i)
  DO j = 1,na(i)
   dmion(:) = dmion(:) + z*r(it,i,j,:) 
  END DO
 END DO 

! dmion(:) = dmion(:)                !*** IONIC POLARIZATION
 avdmion = avdmion + dmion

 DO i = 1,nst
  dmel(:) = dmel(:) - FLOAT(epq)*wf(it,i,:)
 END DO

! dmel = dmel                        !*** ELECTRONIC POLARIZATION
 avdmel = avdmel + dmel

 pol(it,:) = dmion(:) + dmel(:)
 avpol = avpol + pol(it,:)
 pmag = DSQRT(DOT_PRODUCT(pol(it,:),pol(it,:)))
 IF (DABS(pmag) > 1.0d-8) pol(it,:) = pol(it,:)/pmag  

 WRITE(3,"(I8,1X,10(F10.6,2X))") it, dmion(1:3), dmel(1:3), pol(it,1:3)/vol(it), pmag/vol(it)

 avvol = avvol + vol(it)

END DO
CLOSE(UNIT=3)
 
avdmion = avdmion/FLOAT(nt)/ufact
avdmel = avdmel/FLOAT(nt)/ufact
avpol = avpol/FLOAT(nt)/ufact
avvol = avvol/FLOAT(nt)
IF (inpunits == 0)  avvol = avvol*rau**3

!OPEN(UNIT=4,FILE='polcorr.out')
!CALL opcorrv(nt,corrs,pol,4,ncmax)
!CLOSE(UNIT=4)

WRITE(6,*) 'AVERAGE DIPOLE MOMENT (a.u.):'
WRITE(6,*) 'IONIC -'
WRITE(6,"(3(F18.6,2X))") avdmion(:)
WRITE(6,*) 'ELECTRONIC -'
WRITE(6,"(3(F18.6,2X))") avdmel(:)
!pmag = 2.0d0*dmconv/vol(1)
!WRITE(6,*) avdmel(:)/pmag
WRITE(6,*) 'TOTAL -'
WRITE(6,"(3(F18.6,2X))") avpol(:)

!avpol = avpol*debye                       !**** CONVERT TO Debye UNITS

WRITE(6,*) 'AVERAGE DIPOLE MOMENT IN Debye :'
WRITE(6,"(3(F10.6,2X))") avpol(:)*debye
WRITE(6,*) 'MAGNITUDE -'
WRITE(6,"(F10.6)") DSQRT(DOT_PRODUCT(avpol,avpol))*debye

avpol = avpol*dmconv                       !**** CONVERT TO SI UNITS
volm = vol(1)/ufact**3
WRITE(6,*) 'AVERAGE POLARIZATION IN C/m^2 :'
WRITE(6,"(3(F10.6,2X))") avpol(:)/volm
WRITE(6,*) 'MAGNITUDE -'
WRITE(6,"(F10.6)") DSQRT(DOT_PRODUCT(avpol,avpol))/volm
WRITE(6,*) 'AVERAGE VOLUME (\AA^3) :' 
WRITE(6,"(F18.7)") avvol !, vol(1)

END PROGRAM polarizn

SUBROUTINE readpos

USE box
USE runparam

IMPLICIT NONE
INTEGER :: nline, it, ex, is, ia, i, j, ctr, dum
DOUBLE PRECISION :: htmp(3,3), detr, x(3), iden(3,3)
CHARACTER(LEN=2) :: jnk

iden = 0.0d0
DO i =1,3
iden(i,i) = 1.0d0
END DO

WRITE(0,*) 'Reading coordinates from ', TRIM(ADJUSTL(crdfile))
OPEN(UNIT=1,FILE=TRIM(ADJUSTL(crdfile)))
OPEN(UNIT=12,FILE='ions.xyz')

WRITE(12,*) natm
WRITE(12,*) 
it = 0

IF (rmode) THEN

 WRITE(0,*) 'CPR Output format !'
 ex = 0
 IF (MOD(natm,3) /= 0) ex = 1
 nline = 1 + ex + natm/3

 DO ctr = 1,tbeg-1
  DO i = 1,nline
   READ(1,*)
  END DO
 END DO
 
 DO ctr = 1,nt,nskp
  it = it + 1
  READ(1,100) ((h(it,i,j),i=1,3),j=1,3)
  READ(1,100) (((s(it,is,ia,j),j=1,3),ia=1,na(is)),is=1,nsp)

  DO is = 1,nsp
   DO ia = 1,na(is)
!    CALL pbc(s(it,is,ia,:),iden)
    IF (crdstyle == 0) THEN
     r(it,is,ia,:) = MATMUL(h(it,:,:),s(it,is,ia,:))
    ELSE
     r(it,is,ia,:) = s(it,is,ia,:)
    END IF
!    CALL pbc(r(it,is,ia,:),h(it,:,:))
    WRITE(12,*) sym(is), r(it,is,ia,:)
 
   END DO
  END DO

  htmp = h(it,:,:)
  vol(it) = detr(htmp,3)
 END DO

ELSE

 WRITE(0,*) 'Movie format !' !**** This means a sequence of XYZ coordinates
 nline = 2 + natm

 DO is=1,3
  DO ia=1,3
   h(:,is,ia) = href(is,ia)
  END DO
 END DO

 DO ctr = 1,tbeg-1
  DO i = 1,nline
   READ(1,*)
  END DO
 END DO

 DO ctr = 1,nt,nskp
  it = it + 1
  READ(1,*)
  READ(1,*)
  DO is = 1,nsp
   DO ia = 1,na(is)
    READ(1,*) jnk,s(it,is,ia,1:3)
!    CALL pbc(s(it,is,ia,:),iden)
    IF (crdstyle == 0) THEN
     r(it,is,ia,:) = MATMUL(h(it,:,:),s(it,is,ia,:))
    ELSE  
     r(it,is,ia,:) = s(it,is,ia,:)
    END IF
!    CALL pbc(r(it,is,ia,:),h(it,:,:))
    WRITE(12,*) sym(is), r(it,is,ia,:)
   END DO
  END DO
  htmp = h(it,:,:)
  vol(it) = detr(htmp,3)
 END DO

END IF
CLOSE(UNIT=12)
CLOSE(UNIT=1)

100 FORMAT(9(F10.5))

END SUBROUTINE readpos

SUBROUTINE readwf

USE box
USE runparam

IMPLICIT NONE
INTEGER :: i, j, it, nline, ctr
DOUBLE PRECISION :: a(3)

nline = nst

WRITE(0,*) 'Reading Wannier Function centres ',TRIM(ADJUSTL(wffile))

OPEN(UNIT=2,FILE=TRIM(ADJUSTL(wffile)))

DO ctr = 1,tbeg-1
 DO i = 1,nline
  READ(2,*)
 END DO
END DO

it = 0

DO ctr = 1,nt,nskp

 it = it + 1
! READ(2,*)
! READ(2,*) (sprd(it,i),i=1,nst)
! READ(2,*)
 READ(2,*) ((wf(it,i,j),j=1,3),i=1,nst)         !****Centres in input units
 
! DO i = 1,nst
!  CALL pbc(wf(it,i,:),h(it,:,:))
! END DO

END DO
CLOSE(UNIT=2)

END SUBROUTINE readwf

RECURSIVE FUNCTION detr(a,n) RESULT(res)

INTEGER n,i,j,k
!COMPLEX(KIND=8) :: a(n,n),b(n-1,n-1),res
REAL(KIND=8) :: a(n,n),b(n-1,n-1),res

IF (n == 2) THEN

 res = a(1,1)*a(2,2) - a(1,2)*a(2,1)

ELSE

! res = DCMPLX(0.0d0)
 res = 0.0d0
 b(1:n-1,1:n-1) = a(2:n,2:n)
 res = res + a(1,1)*detr(b,n-1)

 DO j = 2,n-1

  b(1:n-1,1:j-1) = a(2:n,1:j-1)
  b(1:n-1,j:n-1) = a(2:n,j+1:n)
  res = res + ((-1)**(j-1))*a(1,j)*detr(b,n-1)

 END DO

 b(1:n-1,1:n-1) = a(2:n,1:n-1)
 res = res + ((-1)**(n-1))*a(n,n)*detr(b,n-1)

END IF

END FUNCTION detr

SUBROUTINE pbc(v,hmat)

!******** Periodizes into -L/2 to L/2 box

IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: hmat(3,3)
DOUBLE PRECISION, INTENT(INOUT) :: v(3)
DOUBLE PRECISION :: b(3), proj, hinv(3,3),det
INTEGER :: i

CALL inv3(hmat,hinv,det)

b = MATMUL(hinv,v)
DO i=1,3
 proj=DSQRT(DOT_PRODUCT(hmat(:,i),hmat(:,i)))
! WRITE(6,*) b(i), proj, NINT(b(i))
 IF ((b(i) <= -0.5) .OR. (b(i) > 0.5)) THEN
  b(i) = b(i) - NINT(b(i))
 END IF
END DO

v = MATMUL(hmat,b)
!WRITE(6,*) "Changed v :", v(:)

END SUBROUTINE pbc

SUBROUTINE epcount(xt,sw,xst)

USE box
USE runparam

IMPLICIT NONE
INTEGER,INTENT(IN) :: xt, sw, xst
INTEGER :: is, ia, ist, i, j
DOUBLE PRECISION :: rsep, rion(3), rel(3)

DO is = 1,nsp
 DO ia = 1,na(is)

  rion(:) = r(xt,is,ia,:) 

  IF (sw == 0) THEN

   DO ist = 1,nst
    IF (wftg(ist) /= 0) CYCLE
    rel(:) = rion(:) - wf(xt,ist,:)
    rsep = DSQRT(DOT_PRODUCT(rel,rel))
    IF (rsep <= rloc) THEN
     nep(is,ia) = nep(is,ia) + epq
     wftg(ist) = 1
     IF (rsep > 1.0d0) WRITE(6,*) ist, 'IS ', rsep,'AWAY FROM (1)', is, ia
    END IF
!    WRITE(6,*) "WF ",ist," :",rsep
   END DO
  
  ELSE
  
    rel(:) = rion(:) - wf(xt,xst,:)
    CALL pbc(rel,h(xt,:,:))
    rsep = DSQRT(DOT_PRODUCT(rel,rel))
!    WRITE(6,*) "WF ",xst," :",rsep, rloc

    IF (rsep <= rloc) THEN
     nep(is,ia) = nep(is,ia) + epq
     wftg(xst) = 1
     wf(xt,xst,:) = rion(:) - rel(:)
     IF (rsep > 2.0d0) WRITE(6,*) xst, 'IS ', rsep,'AWAY FROM (2)', is, ia
     RETURN
    END IF

!    WRITE(6,*) "WF ",xst," :",wf(xt,xst,:)
!    WRITE(6,*)

  END IF

 END DO
END DO

END SUBROUTINE epcount

SUBROUTINE opcorrv(npts,ns,op,un,nmax)
 
IMPLICIT NONE
INTEGER,INTENT(IN) :: npts,ns,un,nmax
DOUBLE PRECISION, INTENT(IN) :: op(npts,3)
INTEGER :: i,mi,m,mmax,j,nin,ib,n
DOUBLE PRECISION :: opc(nmax),t,nu,nusp,rp,ip
 
opc = 0.0d0
 
WRITE(0,*) 'calculated velocities!'

mmax = (npts/ns-1)*ns + 1
ib=0
 
DO i=1,nmax
 
 nin=0
 ib=ib+1

 DO m=1,mmax,ns
 
  mi=m+i-1
  nin=nin+1
 
  opc(i)=opc(i) + DOT_PRODUCT(op(mi,:),op(m,:))
 
 END DO
 
 opc(i)=opc(i)/float(nin)
 WRITE(un,*) i,opc(i)/opc(1)

 IF (ib == ns) THEN
  mmax=mmax-ns
  ib=0
 END IF
 
END DO

opc=opc/opc(1)

END SUBROUTINE opcorrv
