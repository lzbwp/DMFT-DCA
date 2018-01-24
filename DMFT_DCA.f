cc**************************************************************************** 
        program fmain
	use Global
c	use omp_lib
	implicit none
 
c**************************************************************************** 
 
c       This program simulates the 3-dimensional disorded Anderson-hubbard 
c	model in the Cluster Typical Medium Theory Approximation
c
c                        H=H1+H2
c                             +
c                H1= -t SUM ((C (i,s)C(j,s)) + h.c. )
c                       <ij>sproperty of
c                               +  
c                    -t' SUM ((C (i,s)C(j,s)) + h.c. )
c                       [ij]s                           i,j n.n.n.
c
c                H2=   SUM (ed + V ) * n(i,s)
c                       i         i    
c
c	V is a random binary variable with <V>=0.  H is particle-hole 
c	symmetric whenever ed=0 and t'=0.  For a two-dimensional system, 
c	this means that 
c
c               +             +    i
c       C  --> C exp(iQ.i) = C (-1)
c        i      i             i
c                                                +              +
c       C  = SUM exp(ik.i) C  --> SUM exp(ik.i) C exp(-iQ.i) = C
c        k    i             i                    i              Q-k
c
c       or if ed=0, G (tau) = -G   (-tau) = G   (beta-tau)
c                    k          Q-k          Q-k
c       and
c                   G (iwn) = - G   (-iwn)  real parts opposite sign
c                    k           Q-k        imaginary parts equal sign
c
c	in real-frequencies, this means that 
c
c                    R          A            *
c                   G (w) = - G   (-w)  =  -G   (-w)
c                    k         Q-k           Q-k
c

c**************************************************************************** 
c	Original serial code by Mark Jarrell 1999-2000
c
c	The modification to implement Typical Medium Theory within a cluster (CTMT) is
c	by Chinedu Ekuma, Zi Yang Meng, and Hanna Terletska 2012
c
c       Checks:
c       ( ) cat *.f > onefile ; mv onefile onefile.f
c           ftnchek -usage=303 -portability=all onefile.f
c
c
c**************************************************************************** 

        integer n,m,i,j,ic,jc,info
        real(kind) ::  time1,time2,mm,r1,r2,r3,r4,r5,r6,e
        character*72 linemc
cH:
        complex*16 ::  Gsc_temp_odd1(2,2),Det
        complex*16 :: A1,B1,C1,D1,C2
      
c**************************************************************************** 
c        call profinit()
c        call profstart('main')

	
	time1= secnds(0.0)
	                   
	
                                                                                                                                                                         
	
	
cc	wtime= OMP_get_wtime()
c       open the output file
        open(unit=40,file='sigma_aa.dat',status='unknown')
        open(unit=41,file='sigma_bb.dat',status='unknown')
        open(unit=39,file='sigma_ab.dat',status='unknown')
c        
c        open(unit=603,file='Determinant.dat',status='unknown')
c        open(unit=41,file='dos.dat',status='unknown')
        open(unit=42,file='run.dat',status='unknown')
c        open(unit=404,file='Gloc_inv_11_22.dat',status='unknown')
c        open(unit=44,file='chieta.dat',status='unknown')
c        open(unit=47,file='gamma_wt.dat',status='unknown')
c        open(unit=4005,file='pdos.dat',status='unknown')
c        open(unit=48,file='sigma_only.dat',status='unknown')
c       open(unit=200,file='dos-vs-iter.dat',access='append')
c        open(unit=300,file='Gamma_old_aa_bb.dat',status='unknown') 
        open(unit=136,file='Gamma_aa.dat',status='unknown')       
        open(unit=137,file='Gamma_bb.dat',status='unknown') 
        open(unit=138,file='Gamma_ab.dat',status='unknown') 
   
        open(unit=131,file='GK_aa_bb.dat',status='unknown') 
        open(unit=132,file='Gc_aa_bb.dat',status='unknown') 
c       open(unit=129,file='odd_aa_bb.dat',status='unknown') 
c        open(unit=124,file='error_vs_iter.dat',status='unknown')
c	open(unit=201,file='gamma-vs-iter.dat',access='append')
c	open(unit=555,file='norm-gsc-vs-iter',access='append')
c        open(unit=5555,file='norm-gaver-vs-iter',access='append')

c	open(unit=201,file='gamma-vs-iter.dat',access='append')
c       Readin the parameters
        call readin
        eta=ii*eta1
c        aa=1.5d0
c        bb=0.5d0
c        ab=1.5
c        ab=(aa+bb)/2.
c        write(*,*) 'taa,tbb,tab===',aa,bb,ab
c        write(*,*) 'V1,V2,V12===',V1a,V2a,V12a
        ba=ab
        E0_aa=0.
        E0_bb=0.
        E0_ab=0.
c	Allocate arrays
        
	call allocate_arrays
	write(*,*) 'allocate arrays-2'
	
c       we must now make the lookup tables
        call tables
 	write(6,*) '***************************************************************'					      
	write(6,*) "table completed successfully."
	isf=0
	test=0.d0	!convergence criteria
	testl=0.d0
        histo_meas=zeror
	do 999 iter=1,niter					
	  if(isi.ge.1.and.iter.gt.1) meas=int(meas*ifac)		
c	Added to enable the sampling of Large meas for histogram
	  if(histo_enable) then
	  if(isi.ge.1.and.iter.gt.1.and.iter.ne.niter) meas=int(meas*ifac)		
	  if(isi.ge.1.and.iter.eq.niter) meas=int(meas*10)	
	end if ! End histogram 	
	  write(6,*) 'iteration',iter,' runs=',run,' meas=',meas
	  write(42,*) 'iteration',iter,' runs=',run,' meas=',meas
	
	if(iter.eq.1) then ! set initial sigma

       if (.not. odd .and. .not. SC) then
	sigma(:,:)=cmplx(0.d0,-0.02d0)
	GcKf = zeroc
	GcKfsc = zeroc
	do ic=1,Nc
	  do n =-nwn,nwn
            do i=0,e_div
              GcKf(ic,n)=GcKf(ic,n)+ 
     &       PDOS(i,ic)*estep/(wn(n)-(estart+float(i)*estep)- ed - 
     & 	     sigma(ic,n))
             end do
	     GcKf(ic,n) = GcKf(ic,n)*real(Nc)
	  end do
	end do
c	  GcKfsc=GcKf
	do n=-nwn,nwn
	  GcKfsc(:,n)=1.d0/(sigma(:,n)+1.d0/GcKf(:,n))
	end do    

	  Gamma_old = zeroc
	  Gamma_new = zeroc
	  Gammarun_old = zeroc
	  Gammarun_new = zeroc
	do n=-nwn,nwn
	  Gamma_old(:,n) = wn(n)-Epsbar(:)-1./GcKfsc(:,n)
	  Gammarun_old(:,n) = wn(n)-Epsbar(:)-1./GcKfsc(:,n)
c	  Gamma_new(:,n) = wn(n)-sigma(:,n)-Epsbar(:)-1.d0/GcKf(:,n)
	end do

        end if

c*********************************************************************
c                         ODD block
c*********************************************************************
         if (odd) then

c 1. initialize Gamma_old
      
        Gamma_old_aa(:,:)=cmplx(0.d0,-0.1d0)
	Gamma_old_bb(:,:)=cmplx(0.d0,-0.100d0)
        Gamma_old_ab(:,:)=cmplx(0.d0,-0.100d0)
        sigma(:,:)=cmplx(0.d0,-0.02d0)


        Gamma_new_aa(:,:)=zeroc
        Gamma_new_bb(:,:)=zeroc      
        Gamma_new_ab(:,:)=zeroc
        Gamma_new_ba(:,:)=zeroc
c another way:

          D1=zeroc
          GcKf_aa = zeroc
          GcKf_bb = zeroc
          GcKf_ab = zeroc
          GcKf_ba = zeroc
        do ic=1,Nc
         do n=-nwn,nwn

           do i=1,e_div
             e=estart+float(i)*estep
            A1=wn(n)+0.01*ii-aa*e-sigma(ic,n)
            B1=wn(n)+0.01*ii-bb*e-sigma(ic,n)
            C1=ab*e+sigma(ic,n)
            D1=A1*B1-C1**2
     
             GcKf_aa(ic,n)=GcKf_aa(ic,n)+PDOS(i,ic)*estep*B1/D1
             GcKf_bb(ic,n)=GcKf_bb(ic,n)+PDOS(i,ic)*estep*A1/D1
             GcKf_ab(ic,n)=GcKf_ab(ic,n)+PDOS(i,ic)*estep*C1/D1
             GcKf_ba(ic,n)=GcKf_ab(ic,n)
             
           end do !eps-loop
c normalize
          GcKf_aa(ic,n) = GcKf_aa(ic,n)*real(Nc)
          GcKf_bb(ic,n) = GcKf_bb(ic,n)*real(Nc)
          GcKf_ab(ic,n) = GcKf_ab(ic,n)*real(Nc)
          GcKf_ba(ic,n) = GcKf_ba(ic,n)*real(Nc)
          end do
          end do

c invert it to construct Delta=w-sigma-G^-1

         D1=zeroc
         Gloc_inv=zeroc
         do ic=1,Nc
         do n=-nwn,nwn
c        D1(ic,n)=GcKf_aa(ic,n)*GcKf_bb(ic,n)-GcKf_ab(ic,n)**2
c         D1(ic,n)=GcKf_aa(ic,n)*GcKf_bb(ic,n)-GcKf_ab(ic,n)*GcKf_ba(ic,n)
   
         Gloc_inv(1,1,ic,n)=GcKf_bb(ic,n)/(GcKf_aa(ic,n)*GcKf_bb(ic,n)-GcKf_ab(ic,n)*GcKf_ba(ic,n))
         Gloc_inv(2,2,ic,n)=GcKf_aa(ic,n)/(GcKf_aa(ic,n)*GcKf_bb(ic,n)-GcKf_ab(ic,n)*GcKf_ba(ic,n))
         Gloc_inv(1,2,ic,n)=-GcKf_ab(ic,n)/(GcKf_aa(ic,n)*GcKf_bb(ic,n)-GcKf_ab(ic,n)*GcKf_ba(ic,n))
         Gloc_inv(2,1,ic,n)=-GcKf_ba(ic,n)/(GcKf_aa(ic,n)*GcKf_bb(ic,n)-GcKf_ab(ic,n)*GcKf_ba(ic,n))

	end do
        end do

   
c  calculate Gamma_old
        do ic=1,Nc
          do n=-nwn,nwn
            Gamma_old_aa(ic,n)=wn(n)-sigma(ic,n)-Gloc_inv(1,1,ic,n)
            Gamma_old_bb(ic,n)=wn(n)-sigma(ic,n)-Gloc_inv(2,2,ic,n)
            Gamma_old_ab(ic,n)=-sigma(ic,n)-Gloc_inv(1,2,ic,n)
            Gamma_old_ba(ic,n)= Gamma_old_ab(ic,n)
c            Gamma_old_ba(ic,n)=wn(n)-1.*ba*sigma(ic,n)-Gloc_inv(2,1,ic,n)

         end do
        end do
        
        
c        Gamma_old_aa(:,:)=cmplx(0.d0,-0.5d0)
c	Gamma_old_bb(:,:)=cmplx(0.d0,-0.2500d0)
c        Gamma_old_ab(:,:)=cmplx(0.d0,-0.100d0)

! works for V0.3, 0.8
c        Gamma_old_aa(:,:)=cmplx(0.d0,-0.01d0)
c	Gamma_old_bb(:,:)=cmplx(0.d0,-0.02500d0)
c        Gamma_old_ab(:,:)=cmplx(0.d0,-0.0100d0)

! works for V0.3, 0.8
cc        Gamma_old_aa(:,:)=cmplx(0.d0,-0.25d0)
cc	Gamma_old_bb(:,:)=cmplx(0.d0,-0.100d0)
cc        Gamma_old_ab(:,:)=cmplx(0.d0,-0.0100d0)

c        Gamma_old_aa(:,:)=cmplx(0.d0,-0.25d0)
c	Gamma_old_bb(:,:)=cmplx(0.d0,-0.100d0)
c        Gamma_old_ab(:,:)=cmplx(0.d0,-0.500d0)

c also works
c        Gamma_old_aa(:,:)=cmplx(0.d0,-0.1d0)
c	Gamma_old_bb(:,:)=cmplx(0.d0,-0.100d0)
c       Gamma_old_ab(:,:)=cmplx(0.d0,-0.0100d0)
c tested below

       Gamma_old_aa(:,:)=cmplx(0.d0,-0.5d0)
	Gamma_old_bb(:,:)=cmplx(0.d0,-0.500d0)
        Gamma_old_ab(:,:)=cmplx(0.d0,-0.0100d0)


c         Gamma_old_aa(ic,n)=wn(n)+ii*eta-V0-0.5/GcKf_aa(ic,n)
c	 Gamma_old_bb(ic,n)=wn(n)+ii*eta+V0-0.5/GcKf_bb(ic,n)
c         Gamma_old_ab(ic,n)=cmplx(wn(n),-0.01d0)



c 2.calculate Gscript_aa,bb,ab components
c***************************************************************************************
c	form G_script to invert
c***************************************************************************************
        if (Nc.eq.1) then
        Epsbar(1)=0.
        end if


         GcKfsc_aa=zeroc
         GcKfsc_bb=zeroc
         GcKfsc_ab=zeroc
         GcKfsc_ba=zeroc


         IF (new1) then
         do ic=1,Nc
         do n=-nwn,nwn
           GcKfsc_aa(ic,n)= 1.d0/(wn(n)+ii*eta*0.-Epsbar(ic)*aa-Gamma_old_aa(ic,n))
           GcKfsc_bb(ic,n)= 1.d0/(wn(n)+ii*eta*0.-Epsbar(ic)*bb-Gamma_old_bb(ic,n))
           GcKfsc_ab(ic,n)= 1.d0/(wn(n)+ii*eta*0.-Epsbar(ic)*ab-Gamma_old_ab(ic,n))
           GcKfsc_ba(ic,n)=1.d0/(wn(n)+ii*eta*0.-Epsbar(ic)*ba-Gamma_old_ba(ic,n))

c           GcKfsc_aa(ic,n)= 1.d0/(sigma(ic,n)+1./GcKf_aa(ic,n))
c           GcKfsc_bb(ic,n)= 1.d0/(sigma(ic,n)+1./GcKf_bb(ic,n))
c           GcKfsc_ab(ic,n)= 1.d0/(sigma(ic,n)+1./GcKf_ab(ic,n))
c          GcKfsc_ba(ic,n)= GcKfsc_ab(ic,n)         
	 end do
         end do
        end if !new1

      
         IF (new2) then
         do ic=1,Nc
         do n=-nwn,nwn
           GcKfsc_aa(ic,n)= (wn(n)+ii*eta*0.-Epsbar(ic)*aa-Gamma_old_aa(ic,n))
           GcKfsc_bb(ic,n)= (wn(n)+ii*eta*0.-Epsbar(ic)*bb-Gamma_old_bb(ic,n))
           GcKfsc_ab(ic,n)= (wn(n)+ii*eta*0.-Epsbar(ic)*ab-Gamma_old_ab(ic,n))
           GcKfsc_ba(ic,n)= (wn(n)+ii*eta*0.-Epsbar(ic)*ba-Gamma_old_ba(ic,n))     
	 end do
         end do
        END IF !new2
   

c   go to geod.f to make it in real space
c *************************************************************
        end if !if odd
        
c ----  MB
        if (MB) then


        Gamma_old_aa(:,:)=cmplx(0.d0,-0.1d0)
	Gamma_old_bb(:,:)=cmplx(0.d0,-0.100d0)
        Gamma_old_ab(:,:)=cmplx(0.d0,-0.100d0)
        sigma(:,:)=cmplx(0.d0,-0.02d0)

c ----- end MB        
       end if 
        
        if (SC) then
        ed=-U/2.
         ni(:)=1.
c         sigma_aa(:,:)=cmplx(0.d0,-0.04d0)
c         sigma_bb(:,:)=cmplx(0.d0,-0.04d0)
	

	if(iso.eq.1) then ! readin sigma from a file	
	
	   do ic=1,Nc 					      
	      do n=-nwn,nwn				      
	  	read(40,*) mm,r1,r2
	  	read(41,*) mm,r3,r4
	  	read(39,*) mm,r5,r6

	  	sigma_aa(ic,n)=cmplx(r1,r2)
		sigma_bb(ic,n)=cmplx(r3,r4)
		sigma_ab(ic,n)=cmplx(r5,r6)
              end do						      
              end do
              
              write(6,*) '**********************************************************'						      
	      write(6,*) 'code initialized with last sigma_aa_bb_ab'
 	write(6,*) '**********************************************************'	      
	  close(39)
	  close(40)
	  close(41)
          open(unit=40,file='sigma_aa.dat',status='unknown')
          open(unit=41,file='sigma_bb.dat',status='unknown')
          open(unit=39,file='sigma_ab.dat',status='unknown')
       end if	
	if(iso.eq.0) then ! readin sigma from a guess	
	

          sigma_aa(:,:)=cmplx(0,0.)
          sigma_bb(:,:)=cmplx(0,0.)
	  sigma_ab(:,:)=cmplx(0.001,0.0)
	  
 	write(6,*) '***************************************************************'					      
	      write(6,*) 'code initialized with a guess of sigma'	     
 	write(6,*) '***************************************************************' 
	write(42,*) 'code initialized with a guess of the sigma'	 
	end if
	 

c  c.5.2---Construct inverted lattice GF GcKf(k,K,w)GcKf=[Gc^(-1)+\Delta-Eps(k)+Eps_bar(K)]^(-1) which is then coarse-grained
	  if(iso.le.1) then 
          D1=zeroc
          GcKf_aa = zeroc
          GcKf_bb = zeroc
          GcKf_ab = zeroc
          GcKf_ba = zeroc
c        if (.not. pd_int) then
        do ic=1,Nc
         do n=-nwn,nwn
           do i=0,e_div
             e=estart+float(i)*estep
            A1=wn(n)+eta+ed-e-sigma_aa(ic,n)
            B1=wn(n)+eta-ed+e-sigma_bb(ic,n)
            C1=-sigma_ab(ic,n)
            C2=-sigma_ab(ic,n)
           D1=A1*B1-(C1)**2

             GcKf_aa(ic,n)=GcKf_aa(ic,n)+PDOS(i,ic)*estep*B1/D1
c             GcKf_bb(ic,n)=GcKf_bb(ic,n)+PDOS(i,ic)*estep*(A1-Epsbar(ic) )/D1
             GcKf_ab(ic,n)=GcKf_ab(ic,n)-PDOS(i,ic)*estep*(C1 )/D1

             
           end do !eps-loop
c normalize
          GcKf_aa(ic,n) = GcKf_aa(ic,n)*real(Nc)
          GcKf_ab(ic,n) = GcKf_ab(ic,n)*real(Nc)
          GcKf_ba(ic,n) = GcKf_ab(ic,n)*real(Nc)


         end do !wn-loop
        end do !ic-loop
        
        do ic=1,Nc
          do n=-nwn,nwn
           GcKf_bb(ic,n) = -CONJG(GcKf_aa(ic,-n))*real(Nc)
             end do
        end do  
               
                
c  invert local lattice GF which is needed to update Delta: Delta_new=(w-Sigma-Gc^-1)

        D1=zeroc
        Gloc_inv=zeroc
        do ic=1,Nc
        do n=-nwn,nwn
         D1=GcKf_aa(ic,n)*GcKf_bb(ic,n)-GcKf_ab(ic,n)**2
   
         Gloc_inv(1,1,ic,n)=GcKf_bb(ic,n)/D1
         Gloc_inv(2,2,ic,n)=GcKf_aa(ic,n)/D1
         Gloc_inv(1,2,ic,n)=-GcKf_ab(ic,n)/D1
         Gloc_inv(2,1,ic,n)=-GcKf_ba(ic,n)/D1

	end do
        end do
        
c                  do n=-nwn,nwn
          
c          write(404,*) wn(n),-aimag(Gloc_inv(1,1,ic,n)),-real(Gloc_inv(1,1,ic,n)),
c     &    	 -aimag(Gloc_inv(2,2,ic,n)),-real(Gloc_inv(2,2,ic,n))
c     &		 -aimag(Gloc_inv(1,2,ic,n)),-real(Gloc_inv(1,2,ic,n))
c     &		 -aimag(Gloc_inv(2,1,ic,n)),-real(Gloc_inv(2,1,ic,n))
c          end do
          
c   invert back by hand
         D1=zeroc
        do ic=1,Nc
        do n=-nwn,nwn
         D1=Gloc_inv(1,1,ic,n)*Gloc_inv(2,2,ic,n)-Gloc_inv(1,2,ic,n)**2
   
         GcKf_aa(ic,n)=Gloc_inv(2,2,ic,n)/D1
         GcKf_bb(ic,n)=Gloc_inv(1,1,ic,n)/D1
         GcKf_ab(ic,n)=-Gloc_inv(1,2,ic,n)/D1
         GcKf_ba(ic,n)=-Gloc_inv(2,1,ic,n)/D1
c         write(603,*) wn(n),-real(D1)/pi

	end do
        end do      
       
        
        Gamma_old_aa(:,:)=zeroc
	Gamma_old_bb(:,:)=zeroc
        Gamma_old_ab(:,:)=zeroc
        
        do ic=1,Nc
          do n=-nwn,nwn
           Gamma_old_aa(ic,n)=wn(n)+ed+eta-sigma_aa(ic,n)-Gloc_inv(1,1,ic,n)
c           Gamma_old_bb(ic,n)=wn(n)-ed+eta-sigma_bb(ic,n)-Gloc_inv(2,2,ic,n)
           Gamma_old_ab(ic,n)=-Gloc_inv(1,2,ic,n)-sigma_ab(ic,n)
          end do
        end do  
       
        do ic=1,Nc
          do n=-nwn,nwn
          Gamma_old_bb(ic,n)=-CONJG(Gamma_old_aa(ic,-n))
         end do
        end do         

       end if
  
  	if(iso.eq.2) then ! readin gamma from a file	
	
	   do ic=1,Nc 					      
	      do n=-nwn,nwn				      
	  	read(136,*) mm,r1,r2
	  	read(137,*) mm,r3,r4
	  	read(138,*) mm,r5,r6

	  	Gamma_old_aa(ic,n)=cmplx(r1,r2)
		Gamma_old_bb(ic,n)=cmplx(r3,r4)
		Gamma_old_ab(ic,n)=cmplx(r5,r6)
              end do						      
              end do
              
              write(6,*) '**********************************************************'						      
	      write(6,*) 'code initialized with last Gamma_aa_bb_ab'
 	write(6,*) '**********************************************************'	      
	  close(136)
	  close(137)
	  close(138)
        open(unit=136,file='Gamma_aa.dat',status='unknown')       
        open(unit=137,file='Gamma_bb.dat',status='unknown') 
        open(unit=138,file='Gamma_ab.dat',status='unknown') 
        end if
  
        end if ! sc
c          stop
c	    Write out a data header for the binned data
	    call output(0)
	  end if
 1000	  format(a72)
	  

c         Make Gcsfsc, the greens function for all V(l)=0 and initialize all the accumulators
          call geodgen

c         Initialize the greens functions and the fields.
          call ginit

 2000	  continue
 
c         The following is the main part of the main program.  It directs
c         measurement, and updating the lattice and the d greens functions
c         between the measurements.
          do 30 nrun=1,run
            do 40 nmeas=1,meas
              call update
              call measure 

 40         continue
c           Now load the data into bins and accumulators.
	   
            call sumup
 30       continue
	  
          call output(1)

	if(iter.eq.niter) then
	
	!	Try to calculate spectra
c	if(interpolation) then
c	call spectra
c	end if
	
 	write(6,*) "Writing Data to Output"
	end if
          time2= secnds(time1)
 	  write(6,*) 'done iteration at time ',time2
 	  write(6,*) '***********************************************************'		      
 	  write(42,*) 'done iteration at time ',time2
 	  write(42,*) '**********************************************************'		      

	if(iter.eq.niter-1) then
 	write(6,*) "Last Iteration. Major output will be written into dos.dat"
 	write(6,*) "The TotalDOS and PDOS will be written into p_rho.dat"
 	write(6,*) "To understand other outputs, see README, output section"
 	  write(6,*) '***********************************************************'		      
 	write(42,*) "Last Iteration. Major output will be written into dos.dat"
 	write(42,*) "The TotalDOS and PDOS will be written into p_rho.dat"
 	write(42,*) "To understand the many outputs, see README output section"
 	  write(42,*) '**********************************************************'		      
	end if

 999	continue	! done iteration loop

	call deallocate_arrays
c        time2= secnds(time1)
cc	wtime= OMP_get_wtime()-wtime
c	write(6,*) 'cpu ,time=  ',time2
c        write(42,*) 'cpu ,time=  ',time2

c	write(6,*) 'cpu ,time=  ',wtime
c	write(42,*) 'cpu ,time=  ',wtime

c        call profend('main')
c        call profstat()

        stop
        end
