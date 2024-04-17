	subroutine sigmodel_calc(e1pass,e2pass,thpass,zpass,apass,mpass,sig_dis_pass,sig_qe_pass,sig,xflag,factpass)
C       +______________________________________________________________________________
c	
C       Calculate cross section using Peter's F1F209.f routine
c       
c       ARGUMENTS:
c       
c       E1:		-	Incident energy in GeV.
c       E2:		- Scattered energy in GeV.
c       TH:		- Scattering angle in Degrees.
c       A:		- 'A' of nucleus.
c       Z:		- Number of protons in nucleus.
c       M_TGT:	- Mass of target nucleus in GeV/c2.
c       M_REC:	- Mass of recoiling nucleon in GeV/c2.
c       E_SEP:	- Separation energy for target nucleus in GeV/c2.
c       SIG  :	- Calculated cross section in nb/(MeV-ster).
C       ______________________________________________________________________________

        implicit none
	include 'math_physics.inc'
	include 'include.inc'

C       Declare arguments.

	real*8		e1,e2,th,a,z,m_tgt
	real*4          e1pass,e2pass,thpass,mpass
	real*4          sig,factpass,sig_dis_pass,sig_qe_pass
	integer         zpass,apass
	logical		modelhack
	real*8          ff1,deltae,e1cc,e2cc,aux(7) ! added by smoran
	real*8          dum0,dum1,dum2,dum3,dum4,dum5
C       Declare locals.

	real*8	 	sig_qe,sig_dis,y,normfac,fact
        real*8		thr,cs,sn,tn,elastic_peak
	real*8          Q2,nu,WSQ, x
	real*8          F1,F2,W1,W2,sigmott,r
	real*8          inelastic_it
	integer         xflag !flag for which xsec to calculate 1=both 2=QE only 3=DIS only
	logical         first


	save

        real*8 emc_func_xem
	external emc_func_xem

	real*8 emc_func_slac
	external emc_func_slac

	data first/.true./

	e1=dble(e1pass)
	e2=dble(e2pass)
	th=dble(thpass)
	a=dble(apass)
	z=dble(zpass)
	m_tgt=dble(mpass)


        deltae=0.0
        ff1 = (e1+deltae)/e1
        e1cc = e1 + deltae
        e2cc = e2 + deltae


	sig =0.0
	sig_qe=0.0
	sig_dis=0.0
	fact = 0.0		!added by smoran 
C       If at or beyond elastic peak, return ZERO!

	thr = th*d_r
	cs = cos(thr/2.)
	sn = sin(thr/2.)
	tn = tan(thr/2.)
	elastic_peak = e1/(1.+2.*e1*sn**2/m_tgt)
      	if (e2.ge.elastic_peak) then
       	   sig = 0.0
       	   return
       	endif
c	write(6,*) 'got to 1'

	Q2 = 4.*e1*e2*sn**2
	nu=e1-e2
	WSQ = -Q2 + m_p**2 + 2.0*m_p*nu 
        x = Q2/2/m_p/nu



	F1=0
	F2=0
	r=0
	if((xflag.eq.1).or.(xflag.eq.3)) then
c----------------------------------------------------------------
c       
c       do inelastic stuff
	   call F1F2IN09(Z, A, Q2, WSQ, F1, F2, r)
C       Convert F1,F2 to W1,W2
	   W1 = F1/m_p
	   W2 = F2/nu
C       Mott cross section
	   sigmott=(19732.0/(2.0*137.0388*e1*sn**2))**2*cs**2/1.d6
	   sig_dis = 1d3*sigmott*(W2+2.0*W1*tn**2)
CDG apply "iteration" correction (used for XEM analysis)
CDG DO not use this for more "generic" stuff.
CDG	   sig_dis = sig_dis*inelastic_it(x,A)
	   write(6,*) 'IN part: ',  sig_dis, sigmott
	endif


	if((xflag.eq.1).or.(xflag.eq.2)) then
	   if(a.gt.1.5) then
	      call target_info(0.0,dum0,dum1,dum2,z,a,dum3,dum4,aux,dum5)
	      write(6,*)  'before ', e1cc,e2cc,th,a,z,m_tgt,aux,sig_qe,y,fact
	      call f_to_sig(e1cc,e2cc,th,a,z,m_tgt,aux,sig_qe,y,fact)
	      write(6,*)  'after ', e1cc,e2cc,th,a,z,m_tgt,aux,sig_qe,y,fact
	   else
	      sig_qe=0.0
CC	      PRINT *, "Defaulting to zero!: ", sig_qe
	   endif
CC	   sig_qe = sig_qe*ff1**2  ! added by smoran
	endif

	sig = sig_qe + sig_dis !sig is already real*4
	write(6,*) 'sig, sig_qe, sig_dis :', sig , sig_qe, sig_dis
CC       DG this is no longer needed.
Cc	sig = sig*1000.0 !to be consistent with Peter's units

Cc	write(6,*) 'leaving sigmodel_calc_2013'

Cc	sig=sig_qe

	sig_qe_pass = sig_qe ! pass back as real*4
	sig_dis_pass = sig_dis

	return
	end


	real*8 function inelastic_it(x,A)
C DJG: Correction to inelastic cross section to F1F209 from XEM
C DJG: 40 degree data. Just a simple one-pass iteration for use
C DJG: to check our model dependence.
	real*8 A,x
	real*8 p1,p2,p3,p4
	real*8 x1,x2,xit


	inelastic_it = 1.0 ! set to 1 by default

	if(A.gt.1.5.and.A.lt.2.5) then
	   x1=0.3172d0
	   x2=0.7275d0
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=1.2394d0 - 2.1805d0*xit + 5.6853d0*xit**2
     >     -4.3908d0*xit**3
	endif

	if(A.gt.2.5 .and. A.lt.3.5) then
	   x1=0.3172
	   x2=1.055
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=2.9235d0 - 16.075d0*xit + 46.426d0*xit**2
     >     -56.779d0*xit**3 + 25.007d0*xit**4
	endif

	if(A.gt.3.5 .and. A.lt.4.5) then
	   x1=0.3172
	   x2=1.0927
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=1.505d0 - 4.8103d0*xit + 15.221d0*xit**2
     >     -19.713d0*xit**3 + 8.9693d0*xit**4
	endif

	if(A.gt.8.5 .and. A.lt.9.5) then
	   x1=0.6478
	   x2=1.0929
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=0.89139d0 + 0.44009d0*xit + -0.44163d0*xit**2
	endif

	if(A.gt.11.5 .and. A.lt.13.5) then
	   x1=0.3172
	   x2=1.005
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=0.88964d0 + 0.39884d0*xit - 0.36051d0*xit**2
	endif

	if(A.gt.61.0 .and. A.lt.66.0) then
	   x1=0.3172
	   x2=1.005
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=0.77646d0 + 0.90481d0*xit - 0.83815d0*xit**2
	endif


	if(A.gt.196.0 .and. A.lt.198.0) then
	   x1=0.3172
	   x2=1.005
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=0.70439d0 + 1.0510d0*xit - 0.91679d0*xit**2
	endif

	return
	end

	   
	      
	   

	
