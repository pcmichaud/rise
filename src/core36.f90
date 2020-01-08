! Fonseca, Michaud, Galama and Kapteyn
! This code solves the DP model, performs simulations and estimation
!-----------------------------------------------------------------------------------*/

module dphealth
	use newuoa_module
	implicit none
	include 'mpif.h'
  	include 'params36.h'

	contains

		subroutine prepare_estimation
			integer i, n, a
            logical fexist


			if (ismaster) write(10,*) '* loading parameters of the estimation problem'
			! loading parameters from etimation.info
			open(1,file='../params/common/estimation.info')
				read(1,*) buffer, npar
				read(1,*) buffer, alg
				read(1,*) buffer, nsim
				read(1,*) buffer, iweight
				read(1,*) buffer, iload
			close(1)
			if (ismaster) write(10,*) '- estimation using ', alg
			if (ismaster) write(10,*) '- number of simulated persons :', nsim
                        if (ismaster) then
                           if (iweight) then
                              write(10,*) '- moments weighted '
                           else
                              write(10,*) '- moments unweighted'
                           end if
                        end if

                        if (ismaster) then
                           INQUIRE(FILE='../output/common/hdraws.dat', EXIST=fexist)
                           if (fexist .eqv. .false.) then
                           call produce_draws
                           end if
                        end if

                        call MPI_barrier(MPI_COMM_WORLD,ier)
			! shocks (earnings and health)
			allocate(hdraws(nsim,nages))
			allocate(edraws(nsim,nages))
			allocate(idraws(nsim))
			open(1,file='../output/common/hdraws.dat')
			do i = 1, nsim, 1
				do a = 1, nages, 1
					read(1,*) hdraws(i,a)
				end do
			end do
			close(1)

			open(1,file='../output/common/edraws.dat')
			do i = 1, nsim, 1
				do a = 1, nages, 1
                    read(1,*) edraws(i,a)
				end do
			end do
			close(1)

			open(1,file='../output/common/idraws.dat')
			do i = 1, nsim, 1
                read(1,*) idraws(i)
			end do
			close(1)

			if (ismaster) write(10,*) '- draws for earnings and health shocks loaded'
			! load moments from data
			call loadmoments
		end subroutine prepare_estimation

		subroutine produce_draws
			integer i, n, a
			integer, allocatable :: seed(:)
			double precision draw
   			open(1,file='../params/common/estimation.info')
				read(1,*) buffer, npar
				read(1,*) buffer, alg
				read(1,*) buffer, nsim
			close(1)
			call random_seed(size = n)
			allocate(seed(n))
			seed(:) = 12345
			call random_seed(put = seed)
			open(111,file='../output/common/hdraws.dat')
			open(222,file='../output/common/edraws.dat')
			open(333,file='../output/common/idraws.dat')
			do i = 1, nsim, 1
				do a = 1, nages, 1
					call random_number(draw)
					write(111,*) draw

					call random_number(draw)
					write(222,*) draw
				end do
				call random_number(draw)
				write(333,*) draw
			end do
			close(333)
			close(222)
			close(111)

		end subroutine produce_draws

		subroutine prepare_simulation
			integer i, n, a
			integer, allocatable :: seed(:)
			double precision draw
			logical fexist

			if (ismaster) write(10,*) '* loading parameters of the estimation problem'
			! loading parameters from etimation.info
			open(1,file='../params/common/estimation.info')
				read(1,*) buffer, npar
				read(1,*) buffer, alg
				read(1,*) buffer, nsim
			close(1)
			if (ismaster) write(10,*) '- number of simulated persons :', nsim
			! shocks (earnings and health)
                        if (ismaster) then
                           INQUIRE(FILE='../output/common/hdraws.dat', EXIST=fexist)
                           if (fexist .eqv. .false.) then
                           call produce_draws
                           end if
                        end if

                        call MPI_barrier(MPI_COMM_WORLD,ier)

			! shocks (earnings and health)
			allocate(hdraws(nsim,nages))
			allocate(edraws(nsim,nages))
			allocate(idraws(nsim))
			open(1,file='../output/common/hdraws.dat')
			do i = 1, nsim, 1
				do a = 1, nages, 1
					read(1,*) hdraws(i,a)
				end do
			end do
			close(1)

			open(1,file='../output/common/edraws.dat')
			do i = 1, nsim, 1
				do a = 1, nages, 1
                    read(1,*) edraws(i,a)
				end do
			end do
			close(1)

			open(1,file='../output/common/idraws.dat')
			do i = 1, nsim, 1
                read(1,*) idraws(i)
			end do
			close(1)

			if (ismaster) write(10,*) '- draws for earnings and health shocks loaded'
		end subroutine prepare_simulation

		subroutine loadscenario
			character*20 buffer
			open(1,file='../params/scenarios/scenario_'//adjustl(trim(scenario))//'.info')
				read(1,*) buffer, scn_sheet%progress, scn_sheet%iprogress
				read(1,*) buffer, scn_sheet%other, scn_sheet%iother
				read(1,*) buffer, scn_sheet%wages, scn_sheet%iwages
				read(1,*) buffer, scn_sheet%taxes_yr, scn_sheet%itaxes
				read(1,*) buffer, scn_sheet%cashmin_yr, scn_sheet%icashmin
				read(1,*) buffer, scn_sheet%ssben_yr, scn_sheet%issben
				read(1,*) buffer, scn_sheet%coverage, scn_sheet%icoverage
				read(1,*) buffer, scn_sheet%copay_mc, scn_sheet%icopay_mc
				read(1,*) buffer, scn_sheet%copay_ei, scn_sheet%icopay_ei
				read(1,*) buffer, scn_sheet%prices, scn_sheet%iprices
				read(1,*) buffer, scn_sheet%dwage, scn_sheet%idwage
				read(1,*) buffer, scn_sheet%dcopay, scn_sheet%idcopay
			close(1)
		end subroutine

		subroutine auxiliary
			integer i, z, h, a, j, e, hh, ibase, growth
			double precision parh_age(3), parh_lag(3,3), parh_const(3)
			double precision pard_age, pard_lag(3), pard_const, sage, probh(nh)
			double precision span, ne1, rho, sige, sig, sigv, sigu, gap, ageff
			character*10 fmt_age
			fmt = '(6F9.3)'
			fmt_age = '(I4,6F9.3)'
			if (ismaster) write(10,*) '* loading auxiliary parameters and process'
			! productivity parameters for health process
			if (ismaster) write(10,*) '- loading theta (productivity)'
			allocate(theta(2,nh-1))
			open(1,file='../params/input/production-theta.csv')
				read(1,*) theta(1,:)
				if (ismaster) write(10,fmt) theta(1,:)
				read(1,*) theta(2,:)
				if (ismaster) write(10,fmt) theta(2,:)
			close(1)

			! other parameters of health process
			allocate(parh(nparh,nh-1))
			if (ismaster) write(10,*) '- loading other parameters of health process'
			open(1,file='../params/input/production-health.csv')
				do i = 1, nparh, 1
					read(1,*) parh(i,:)
					if (ismaster) write(10,fmt) parh(i,:)
				end do
			close(1)
			parh_age = parh(1,:)
			parh_lag = parh(2:4,:)
			parh_risk = parh(5:6,:)
			parh_const = parh(7,:)

			! parameters of mortality process
			allocate(pard(npard))
			if (ismaster) write(10,*) '- loading parameters of mortality process'
			open(1,file='../params/input/production-death.csv')
				do i = 1, npard, 1
					read(1,*) pard(i)
					if (ismaster) write(10,fmt) pard(i)
				end do
			close(1)
			pard_age = pard(1)
			pard_lag = pard(2:4)
			pard_const = pard(5)
			! read obesity and smoking data
			allocate(obese(nz,nages))
			allocate(smoke(nz,nages))
			if (ismaster) write(10,*) '- loading smoking (1) and obesity (2) data (1940 cohort)'
			open(1, file='../params/input/risk.csv')
				do a = 1, nages, 1
					read(1,*) smoke(zbase,a), obese(zbase,a)
					if (scn_sheet%iother) then
						smoke(zbase,a) = smoke(zbase,1)
						obese(zbase,a) = obese(zbase,1)
					end if								
  					if (ismaster) write(10,'(I4,2F9.3)') ages(a), smoke(zbase,a), obese(zbase,a)
				end do
			close(1)

            if (nz .gt. 1) then
                do z = 1, nz, 1
                    if (z .ne. zbase) then
                        do i = 1, nages, 1
                            ibase = i + gridz(z) - gridz(zbase)
                            if (ibase .lt. 1) then
                                ibase = 1
                            end if
                            if (ibase .gt. nages) then
                                ibase = nages
                            end if
                            smoke(z,i) = smoke(zbase,ibase)
							obese(z,i) = obese(zbase,ibase)
							if (ismaster) write(10,'(2I4,2F9.3)') gridz(z),ages(i), smoke(z,i), obese(z,i)
						end do
                    end if 
                end do
            end if


			! read adjustment for mortality
			allocate(adjust(nages))
			if (ismaster) write(10,*) '- loading adjustment parameters for mortality'
			open(1, file='../params/input/adjustdead.csv')
				read(1,*) par_adjust
				if (ismaster) write(10,fmt) par_adjust
			close(1)
			if (ismaster) write(10,*) 'yields adjustment by age'
			do a = 1, nages, 1
				adjust(a) = dexp(par_adjust(3))
				if (ages(a).lt.60) then
					adjust(a) = adjust(a)*dexp(par_adjust(1)*dble(ages(a)))
				else
					adjust(a) = adjust(a)*dexp(par_adjust(1)*dble(60) + par_adjust(2)*dble(ages(a)-60))
				end if
				if (ismaster) write(10,fmt_age) ages(a), adjust(a)
			end do
			! create index for health processes
			allocate(lambda(nz,nages,nh-1,nh-1))
			do z= 1, nz, 1
				do a = 1, nages, 1
					sage = dble(ages(a))
					do h = 2, nh, 1
						do hh = 2, nh, 1
								lambda(z,a,h-1,hh-1) = dexp(parh_age(hh-1)*sage  &
									+ parh_lag(h-1,hh-1) + parh_risk(1,hh-1)*obese(z,a) &
									+ parh_risk(2,hh-1)*smoke(z,a) + parh_const(hh-1))
						end do
					end do
				end do
			end do

			if (ismaster) write(10,*) '- survival probabilities'
			! survival probabilities
			allocate(probs(nh-1,nages))
			do z = 1, nz, 1
				do h = 2, nh, 1
					if (ismaster) write(10,'(A,I4)') 'base health = ', h
					do a = 1, nages, 1
						sage = dble(ages(a))
						probs(h-1,a) = dexp(-dexp(pard_age*sage &
								+ pard_lag(h-1) + pard_const))
						probs(h-1,a) = 1.0d0 - adjust(a)*(1.0d0 - probs(h-1,a))
						if (probs(h-1,a).lt.0.0d0) then
							probs(h-1,a) = 0.0d0
						end if
					end do
				end do
			end do

			! deterministic earnings
			if (ismaster) write(10,*) '- loading earnings parameters'
			open(2,file='../params/input/earnings.csv')
			do i = 1, 5, 1
				read(2,*) parearn(i)
				if (ismaster) write(10,fmt) parearn(i)
			end do
			close(2)
			allocate(earn(nz, nages))
			if (ismaster) write(10,*) 'yields age profile'
			do z = 1, nz, 1
				do a = 1, nages, 1
					if (scn_sheet%iwages) then
						growth = gridz(z) - gridz(zbase) + a - 1
						growth = - growth
						if (growth .gt. 0) then
							growth = 0
						end if
					else
						growth =  gridz(z) - gridz(zbase)
						if (growth .lt. 0) then
							growth = 0
						end if
					end if
					earn(z,a) = dexp(dble(growth)*scn_sheet%wages)*dexp(parearn(1)*dble(ages(a)) &
							+ parearn(2)*(dble(ages(a))**2) + parearn(3))
					if (scn_sheet%idwage) then
						earn(z,a) = earn(z,a)*(1.0d0 + scn_sheet%dwage)
					end if	
					if (ismaster) write(10,'(I4,I4,F12.3)') gridz(z),ages(a),earn(z,a)
				end do
			end do

			! deterministic other income parameters
			if (ismaster) write(10,*) '- other income profile parameters'
			open(3,file='../params/input/otherinc.csv')
			allocate(otherinc(nz,nages))
			do i = 1, 4, 1
				read(3,*) parother(i)
				if (ismaster) write(10,fmt) parother(i)
			end do
			if (ismaster) write(10,*) 'yields age profile (at zero earning of person)'
			do a= 1, nages, 1
				do z = 1, nz, 1
					if (scn_sheet%iwages) then
						growth = gridz(z) - gridz(zbase) + a - 1
						growth = - growth
						if (growth .gt. 0) then
							growth = 0
						end if
					else
						growth =  gridz(z) - gridz(zbase)
						if (growth .lt. 0) then
							growth = 0
						end if
					end if
					ageff = dble(ages(a))
					if (ageff .gt. 90.0d0) then
						ageff = 90.0d0
					end if
					otherinc(z,a) = parother(2)*ageff + parother(3)*ageff**2 + parother(4)
					otherinc(z,a) = otherinc(z,a)*dexp(scn_sheet%wages*dble(growth))
					if (otherinc(z,a) .lt. 0.0d0) then
						otherinc(z,a) = 0.0D0
					end if
					if (ismaster) write(10,'(I4,2F12.3)') ages(a), otherinc(z,a), otherinc(z,a) + parother(1)*earn(z,a)
				end do
			end do
			close(3)

			! earnings shocks: approximation of Tauchen (1986)
			allocate(probe(ne,ne))
			allocate(pte(ne))
			allocate(cumprobe(ne,ne))
			span = 3.0d0
			if (ne.eq.1) then
				pte(1) = 0.0d0
				probe(1,1) = 1.0d0
				cumprobe(1,1) = 1.0d0
				if (ismaster) write(10,*) '* no uncertainty in earnings'
			else
				if (ismaster) write(10,*) '- computing transition matrices for earnings'
				ne1 = dble(ne) - 1.0d0
				! compute unconditional std
				rho = parearn(4)
				sige = dsqrt(parearn(5))
				! always set sigv to zero
				sigv = 0.0d0
				sigu = dsqrt(sige**2 + sigv**2)
				sig =  dsqrt(sigv**2 + (sige**2)/(1.0d0-rho**2))
				if (ismaster) write(10,'(A,F9.5)') 'unconditional sd of log earnings', sig
				incsig = sig**2
				! grid points
				pte(1) = -span*sig
				pte(ne) = span*sig
				gap = (2.0d0*span*sig)/ne1
				! Fill in between
				do i = 2, ne-1, 1
					pte(i) = pte(1) + dble(i-1)*gap
				end do
				! compute probabilities
				do i = 1, ne, 1
					do j = 1, ne, 1
						if (j.eq.1) then
							probe(i,1) = probn((pte(j)+gap/2.0d0 - rho*pte(i))/sigu)
							cumprobe(i,1) = probe(i,1)

						else if (j.gt.1 .and. j.lt.ne) then
							probe(i,j) = probn((pte(j)+gap/2.0d0 - rho*pte(i))/sigu) &
								- probn((pte(j)-gap/2.0d0 - rho*pte(i))/sigu)
							cumprobe(i,j) = probe(i,j) + cumprobe(i,j-1)
						else
							probe(i,ne) = 1.0d0 - probn((pte(j)-gap/2.0d0 - rho*pte(i))/sigu)
							cumprobe(i,ne) = 1.0d0
						end if
					end do
				end do
			end if
			! scale points (to apply multiplicatively)
			if (ismaster) write(10,*) 'yields distribution of earnings (i, value, transition prob)'
			do i = 1, ne, 1
				pte(i) = dexp(pte(i) - 0.5d0 * incsig**2)
				if (ismaster) write(10,'(I4,12F9.4)') i, pte(i), probe(i,:)
			end do

			allocate(progress(nz,nages))
			do z = 1, nz, 1
				do a = 1, nages, 1
					progress(z,a) = 2005 - (gridz(z) + agemin + a - 1)
					if (progress(z,a) .gt. 40) then
						progress(z,a) = 40
					end if 
					if (scn_sheet%iprogress) then
						progress(z,a) = 40
					end if
					if (ismaster) write(10,'(2I4,I4)') gridz(z), ages(a), progress(z,a)
				end do
				
			end do

		end subroutine auxiliary


		subroutine statespace
			integer i,t,a,h,e,z,f,d,m,w,s,hh, j, ibase

			double precision incopay(4)
			double precision uumin, zip
			! load scenario
			call loadscenario

			if (ismaster) write(10,*) '* creating state-space'
			! load parameters for statespace
			if (ismaster) write(10,*) '- loading parameters for state-space'
			open(1,file='../params/common/solve.info')
				read(1,*) buffer, nstatevars
				read(1,*) buffer, ndecisionvars
				read(1,*) buffer, nages
				read(1,*) buffer, agemin
				read(1,*) buffer, nh
				read(1,*) buffer, ne
				read(1,*) buffer, nz
                read(1,*) buffer, zmin
                read(1,*) buffer, zmax 
				read(1,*) buffer, nf
				read(1,*) buffer, nd
				read(1,*) buffer, nw
				read(1,*) buffer, wmin
				read(1,*) buffer, wmax
				read(1,*) buffer, name
				read(1,*) buffer, amemin
				read(1,*) buffer, amemax
				read(1,*) buffer, nparh
				read(1,*) buffer, npard
				read(1,*) buffer, umin
				read(1,*) buffer, medmax
				read(1,*) buffer, curv
				read(1,*) buffer, isave
				read(1,*) buffer, nsave
				allocate(saveages(nsave))
				read(1,*) buffer, saveages(:)
				read(1,*) buffer, ngridc
				read(1,*) buffer, ngridm
			close(1)
			if (ismaster) then
				write(10,'(A30,I4)') '# of state vars', nstatevars
				write(10,'(A30,I4)') '# of decision vars', ndecisionvars
				write(10,'(A30,I4)') '# of ages', nages
				write(10,'(A30,I4)') 'starting age', agemin
				write(10,'(A30,I4)') '# health states', nh
				write(10,'(A30,I4)') '# earnings shocks', ne
				write(10,'(A30,I4)') '# cohorts', nz
                write(10,'(A30,I4)') 'first cohort', zmin
                write(10,'(A30,I4)') 'last cohort', zmax
				write(10,'(A30,I4)') '# insurance status', nf
				write(10,'(A30,I4)') '# claming states', nd
				write(10,'(A30,I4)') '# wealth points', nw
				write(10,'(A30,F12.1)') 'minimum wealth', wmin
				write(10,'(A30,F12.1)') 'maximum wealth', wmax
				write(10,'(A30,I4)') '# ame points', name
				write(10,'(A30,F12.1)') 'minimum ame', amemin
				write(10,'(A30,F12.1)') 'maximum ame', amemax
				write(10,'(A30,I4)') '# params health', nparh
				write(10,'(A30,I4)') '# params death', npard
				write(10,'(A30,F12.1)') 'max health spend', medmax
				write(10,'(A30,I4)') '# grid points cons', ngridc
				write(10,'(A30,I4)') '# grid points medexp', ngridm

				if (isave) then
					write(10,*) 'saving decision rules in output/...'
				end if
			end if
			! load institutions
			if (ismaster) write(10,*) '- loading institutional parameters (before applying scenario)'
			open(1,file='../params/common/institutions.info')
				read(1,*) buffer, rrate
				read(1,*) buffer, npiabend
				read(1,*) buffer, npiarate
				read(1,*) buffer, nearntax
				read(1,*) buffer, arf
				read(1,*) buffer, drc
				read(1,*) buffer, nra
				read(1,*) buffer, era
				read(1,*) buffer, mcare
				allocate(piabend(npiabend))
				do i = 1, npiabend, 1
					read(1,*) buffer, piabend(i)
				end do
				allocate(piarate(npiarate))
				do i = 1, npiarate, 1
					read(1,*) buffer, piarate(i)
				end do
				allocate(earntaxrate(nearntax))
				do i = 1, nearntax, 1
					read(1,*) buffer, earntaxrate(i)
				end do
				read(1,*) buffer, earntaxgain
				read(1,*) buffer, rmin
				read(1,*) buffer, rmax
				read(1,*) buffer, hours_worked
			close(1)

			! grid for ages and birth year
			allocate(ages(nages))
			do i = 1, nages, 1
				ages(i) = agemin - 1 + i
			end do
            ! birth year
			allocate(gridz(nz))
            if (nz .gt. 1) then
                gapz = int((zmax - zmin)/dble(nz-1))
			    do i = 1, nz, 1
				    gridz(i) = zmin + gapz*(i-1)
			    end do
                zbase = int((1940 - zmin)/gapz) + 1 
            else 
                gapz = 1
                gridz(1) = zmin
                zbase = 1
            end if

			! load social security and medicare taxes
			open(1,file='../params/input/sstax.csv')
			allocate(sstax(nz,nages),ssmax(nz,nages))
            ! read in first for the 1940 cohort
                do i = 1, nages, 1
					read(1,*) sstax(zbase,i),ssmax(zbase,i)
					if (scn_sheet%itaxes) then
						sstax(zbase,i) = sstax(zbase,1)
						ssmax(zbase,i) = ssmax(zbase,1)
					end if
                end do
			close(1)

            if (nz .gt. 1) then
                do z = 1, nz, 1
                    if (z .ne. zbase) then
                        do i = 1, nages, 1
                            ibase = i + gridz(z) - gridz(zbase)
                            if (ibase .lt. 1) then
                                ibase = 1
                            end if
                            if (ibase .gt. nages) then
                                ibase = nages
                            end if
                            sstax(z,i) = sstax(zbase,ibase)
							ssmax(z,i) = ssmax(zbase,ibase)
							if (ismaster) write(10,'(2I4,2F9.3)') gridz(z),ages(i),sstax(z,i), ssmax(z,i)
                        end do
                    end if 
                end do
            end if



			open(1,file='../params/input/mctax.csv')
			allocate(mctax(nz,nages),mcmax(nz,nages))
			do i = 1, nages, 1
				read(1,*) mctax(zbase,i),mcmax(zbase,i)
				if (scn_sheet%itaxes) then
					mctax(zbase,i) = mctax(zbase,1)
					mcmax(zbase,i) = mcmax(zbase,1)
				end if
			end do
			close(1)
            if (nz .gt. 1) then
                do z = 1, nz, 1
                    if (z .ne. zbase) then
                        do i = 1, nages, 1
                            ibase = i + gridz(z) - gridz(zbase)
                            if (ibase .lt. 1) then
                                ibase = 1
                            end if
                            if (ibase .gt. nages) then
                                ibase = nages
                            end if
                            mctax(z,i) = mctax(zbase,ibase)
							mcmax(z,i) = mcmax(zbase,ibase)
							if (ismaster) write(10,'(2I4,2F9.3)') gridz(z),ages(i),mctax(z,i), mcmax(z,i)
                        end do
                    end if 
                end do
            end if

			! generosity of social security
			allocate(ssgen(nz))
			open(1,file='../params/input/pia.csv')
			do i = 1, nages, 1
				read(1,*) zip
				do z = 1, nz, 1
					if (gridz(z) + nra .eq. 1965 + i -1) then
						ssgen(z) = zip
					end if
				end do
			end do
			do z = 1, nz, 1
				if (gridz(z) + nra .lt. 1965) then
					ssgen(z) = 0.657d0
				end if
				if (gridz(z) + nra .gt. 1965 + nages -1) then
					ssgen(z) = 1.0d0
				end if
				if (scn_sheet%issben) then
					ssgen(z) = ssgen(zbase)
				end if
				if (ismaster) write(10,'(I4,F9.4)') gridz(z), ssgen(z)
			end do
			close(1)

			if (ismaster) then
				write(10,'(A30,F9.3)') 'rate of return (before tax)', rrate
				write(10,'(A30,F9.3)') 'ss arf (annual)', arf
				write(10,'(A30,F9.3)') 'ss drc (annual)', drc
				write(10,'(A30,I4)') 'ss era (age)', era
				write(10,'(A30,I4)') 'ss nra (age)', nra
				write(10,'(A30,I4)') 'eligibility medicare (age)', mcare
				write(10,'(A30,2F9.1)') 'ss pia bend points', piabend(:)
				write(10,'(A30,3F9.3)') 'ss pia rates', piarate(:)
				write(10,'(A30,F12.1)') 'ss earnings cap', drc
				write(10,'(A30,2F9.3)') 'ss earnings test rates', earntaxrate(:)
				write(10,'(A30,F12.1)') 'ss earnings test disreg.', earntaxgain
				write(10,'(A30,I4)') 'min age can work', rmin
				write(10,'(A30,I4)') 'max age can work', rmax
				write(10,'(A30,F12.1)') 'hours full-time', hours_worked				
			end if

			! grid for wealth
			if (ismaster) write(10,*) '- creating grid for wealth'
			allocate(gridw(nw))
			allocate(c_gridw(nw))
			c_wmax = dsqrt(wmax)
			c_wmin = dsqrt(wmin)
			gapw = (c_wmax - c_wmin)/dble(nw-1)
			do i = 1, nw, 1
				c_gridw(i) = c_wmin + dble(i-1)*gapw
				gridw(i) = c_gridw(i)**2
				if (ismaster) write(10,'(I4,2F12.1)') i, gridw(i), c_gridw(i)
			end do

			! grid for ame
			allocate(gridame(name))
			if (ismaster) write(10,*) '- creating grid for ame'
			gapame = (amemax - amemin)/dble(name-1)
			do i = 1, name, 1
				gridame(i) = amemin + dble(i-1)*gapame
				gridame(i) = gridame(i)
				if (ismaster) write(10,'(I4,F12.1)') i, gridame(i)
			end do
			ngrid = name*nw
			if (ismaster) write(10,'(A,I4)') '# points ame x wealth', ngrid


			! create the state space matrix (nyears,nstates,nstatevars)
			nspace = nh*ne*nz*nf*nd*nw*name
			if (ismaster) write(10,'(A,I12)') '# potential grid points per age', nspace
			allocate(states(nages,nspace,nstatevars))
			allocate(nstates(nages))
			nmaxstates = 1
			if (ismaster) write(10,*) '# grid points (actual) by age'
			do a = 1, nages, 1
				! initialize counter for state number
				s = 1
				do h = 1, nh, 1
					if (h.eq.1) then
						cycle
					end if
					do e = 1, ne, 1
						if (e.gt.1 .and. ages(a).ge.rmax) then
							cycle
						end if
						! cohort
						do z = 1, nz, 1
							do f = 1, nf, 1
								! medicare once reached 65
								if (ages(a).ge.mcare .and. f.lt.3) then
									cycle
								end if
								do d = 1, nd, 1
									! can only reach ages 25 to 62 with d==1
									if (ages(a).le.era .and. d.eq.2) then
										cycle
									end if
									! can only reach 70 with SS claimed at 70
									if (ages(a).ge.rmax .and. d.eq.1) then
										cycle
									end if
									do m = 1, name, 1
										do w = 1, nw, 1
											states(a,s,:) = (/m,w,h,e,f,d,z/)
											! next state
											s = s+1
										end do
									end do
								end do
							end do
						end do
					end do
				end do
				nstates(a) = s-1
				if (ismaster) write(10,'(I4,I12)') ages(a), nstates(a)
				if (nstates(a).gt.nmaxstates) then
					nmaxstates = nstates(a)
				end if
			end do
			if (ismaster) write(10,'(A,I12)') '# grid points (actual)', sum(nstates(:))

			! other parameters
			allocate(cashmin(nz,nages))
			if (ismaster) write(10,*) '- minimum cash-on-hand by age for 1940 cohort (Moffitt, 2002)'
			open(1,file='../params/input/cmin.csv')
			do a = 1, nages, 1
				read(1,*) cashmin(zbase,a)
				if (scn_sheet%icashmin) then
					cashmin(zbase,a) = cashmin(zbase,1)
				end if
				if (ismaster) write(10,'(I4,F12.1)') 1965 + a - 1, cashmin(zbase,a)
			end do
			close(1)

            if (nz .gt. 1) then
                do z = 1, nz, 1
                    if (z .ne. zbase) then
                        do i = 1, nages, 1
                            ibase = i + gridz(z) - gridz(zbase)
                            if (ibase .lt. 1) then
                                ibase = 1
                            end if
                            if (ibase .gt. nages) then
                                ibase = nages
                            end if
							cashmin(z,i) = cashmin(zbase,ibase)
							if (ismaster) write(10,'(2I4,F9.3)') gridz(z),ages(i),cashmin(z,i)
                        end do
                    end if 
                end do
            end if

			! health insurance
			allocate(copay(nz,nages,nf))
			if (ismaster) write(10,*) '- co-pays by insurance status'
			open(1, file='../params/input/coinsurance.csv')
			do i = 1, 4, 1
				read(1,*) incopay(i)
			end do
			close(1)
			allocate(insgen(nz,nages))
			open(1, file='../params/input/oop.csv')
			do i = 1, nages, 1
				read(1,*) insgen(zbase,i)
				if (scn_sheet%icopay_ei) then
					insgen(zbase,i) = insgen(zbase,1)
				end if 
			end do
			close(1)

            if (nz .gt. 1) then
                do z = 1, nz, 1
                    if (z .ne. zbase) then
                        do i = 1, nages, 1
                            ibase = i + gridz(z) - gridz(zbase)
                            if (ibase .lt. 1) then
                                ibase = 1
                            end if
                            if (ibase .gt. nages) then
                                ibase = nages
                            end if
							insgen(z,i) = insgen(zbase,ibase)
							if (ismaster) write(10,'(2I4,F9.3)') gridz(z),ages(i),insgen(z,i)
                        end do
                    end if 
                end do
            end if
	
			if (scn_sheet%idcopay) then
				incopay(1) = scn_sheet%dcopay
				incopay(2) = scn_sheet%dcopay
				incopay(4) = scn_sheet%dcopay
			end if	
			do i = 1, 4, 1
				if (ismaster) write(10,'(I4,F9.3)') i, incopay(i)
			end do
			mcdcopay = incopay(2)
			if (ismaster) write(10,*) 'yields by age'
			do z = 1, nz, 1
				do a = 1, nages, 1
					do f = 1, nf, 1
						if (ages(a).lt.mcare) then
							if (f.eq.1) then
								copay(z,a,f) = incopay(3)
							else if (f.eq.2) then
								copay(z,a,f) = incopay(1)*insgen(z,a)
							else
								copay(z,a,f) = incopay(1)*insgen(z,a)
							end if
						else
							copay(z,a,f) = incopay(4)*insgen(z,a)
						end if
					end do
					if (ismaster) write(10,'(I4,I4,3F9.3)') gridz(z),ages(a), copay(z,a,:)
				end do
			end do

			! AIME computation (replacement rates come from simulations)
			allocate(rep(nages))
			if (ismaster) write(10,*) '- probability earnings do not update AME (French, 2005)'
			open(1,file='../params/input/rep.csv')
			do a = 1, nages, 1
				if (ages(a).le.60) then
					rep(a) = 0.0d0
				else if (ages(a).gt.60 .and. ages(a).lt.71) then
					read(1,*) rep(a)
				else
					rep(a) = 1.0d0
				end if
				if (ismaster) write(10,'(I4,F9.3)') ages(a), rep(a)
			end do
			close(1)

			! tax parameters Gouveia and Strauss
			allocate(taxpar(nz,nages,3))
			if (ismaster) write(10,*) '- tax function parameters (Gouveia and Strauss)'
			open(1,file='../params/input/gouveia.csv')
			do a = 1, nages, 1
				read(1,*) taxpar(zbase,a,:)
				if (scn_sheet%itaxes) then
					taxpar(zbase,a,:) = taxpar(zbase,1,:)
				end if
				if (ismaster) write(10,'(I4,4F9.3)') 1965 + a - 1, taxpar(zbase,a,:)
			end do
			close(1)

            if (nz .gt. 1) then
                do z = 1, nz, 1
                    if (z .ne. zbase) then
                        do i = 1, nages, 1
                            ibase = i + gridz(z) - gridz(zbase)
                            if (ibase .lt. 1) then
                                ibase = 1
                            end if
                            if (ibase .gt. nages) then
                                ibase = nages
                            end if
							taxpar(z,i,:) = taxpar(zbase,ibase,:)
							if (ismaster) write(10,'(I4,I4,3F9.3)') gridz(z),ages(i), taxpar(z,i,:)
                        end do
                    end if 
                end do
            end if

			

		end subroutine statespace

		subroutine initialize(initpar)
			integer i
			type (params), intent(inout) :: initpar(npar)
			nfreepar = 0
			open(1,file='../params/common/init_params_'//adjustl(trim(scenario))//'.info')
			if (ismaster) write(10,*) '* initial parameters '
			do i = 1, npar, 1
 				read(1,*) initpar(i)%label, initpar(i)%value, initpar(i)%lb, initpar(i)%ub, initpar(i)%free
 				if (ismaster) write(10,'(A20,F9.3,F9.3,F9.3,L4)') initpar(i)%label, initpar(i)%value, &
 					initpar(i)%lb, initpar(i)%ub, initpar(i)%free
 				if (initpar(i)%free) then
 					initpar(i)%step = 	0.75d0*(initpar(i)%ub - initpar(i)%lb)
 					nfreepar = nfreepar + 1
 					initpar(i)%ipos = nfreepar
 				else
 					initpar(i)%step = 0.0d0
 				end if
			end do
			close(1)

			if (ismaster) write(10,'(A20,I4)') '# of free parameters ', nfreepar

			allocate(g_initpar(npar))
			g_initpar = initpar
		end subroutine initialize

		subroutine initialize_simulation(initpar)
			integer i
			type (params), intent(inout) :: initpar(npar)
			open(1,file='../params/common/esti_params.info')
			if (ismaster) write(10,*) '* using parameters '
			do i = 1, npar, 1
 				read(1,*) initpar(i)%label, initpar(i)%value, initpar(i)%serr
 				if (ismaster) write(10,'(A20,F9.3,F9.3)') initpar(i)%label, initpar(i)%value, &
 					initpar(i)%serr

 			end do
			close(1)
		end subroutine initialize_simulation

		subroutine loadmoments
			integer i, a, h, j, n, errorflag, amin, amax

			open(1,file='../params/common/moments.info')
			read(1,*) ngroup
			! getting information on moments
			nmoments = 0
			nfreemoments = 0
			allocate(moments(ngroup))
			do i = 1, ngroup, 1
				read(1,*) moments(i)%label, moments(i)%nages, &
					moments(i)%nh, moments(i)%free
				allocate(moments(i)%ages(moments(i)%nages))
				read(1,*) amin,amax
				do n = 1, moments(i)%nages, 1
					moments(i)%ages(n) = amin + n - 1
				end do	
				if (moments(i)%nh) then
					allocate(moments(i)%sim(moments(i)%nages,3))
					allocate(moments(i)%nsim(moments(i)%nages,3))					
					moments(i)%nmom = moments(i)%nages * 3
				else
					allocate(moments(i)%sim(moments(i)%nages,1))
					allocate(moments(i)%nsim(moments(i)%nages,1))
					moments(i)%nmom = moments(i)%nages
				end if
				nmoments = nmoments + moments(i)%nmom
				if (moments(i)%free) then
					nfreemoments = nfreemoments + moments(i)%nmom
				end if
			end do
			close(1)

			! loading number of observations in both datasets
			open(1,file='../params/input/nobs.dat')
				read(1,*) nind
				read(1,*) ndata
			close(1)

			! load adjusted data
			allocate(data(ndata,nfreemoments))
			open(1,file='../params/input/data_adjusted.csv')
				do i = 1, ndata, 1
					read(1,*) data(i,:)
				end do
			close(1)
			allocate(select(ndata,nfreemoments))
			open(1,file='../params/input/select_adjusted.csv')
				do i = 1, ndata, 1
					read(1,*) select(i,:)
				end do
			close(1)

			! compute optimal weighting matrix
			if (iweight) then
				allocate(datacov(nfreemoments,nfreemoments))
				allocate(optw(nfreemoments,nfreemoments))
				allocate(datamean(nfreemoments))
				allocate(datans(nfreemoments))
				datacov(:,:) = 0.0d0
				datamean(:) = 0.0d0
				datans(:) = 0.0d0
				! get ns with non zero data and means
				do i = 1, nfreemoments
					do j = 1, ndata, 1
						if (select(j,i) .eq. 1) then
							datans(i) = datans(i) + 1.0d0
						end if	
					end do	
					datamean(i) = sum(data(:,i))/datans(i)
				end do
				!Covariance matrix
				do i = 1, nfreemoments
					do j = 1, nfreemoments, 1
						do n = 1, ndata, 1
							! whether observation contributes to both moments
							if (select(n,i) .eq. 1 .and. select(n,j) .eq. 1) then
								datacov(i,j) = datacov(i,j) + (data(n,i) - datamean(i))*(data(n,j) - datamean(j))/dble(nind)
							end if
						end do
					end do
				end do
				! call inverse of covariance matrix of data to get optimal weight matrix
				call invert(datacov, optw, nfreemoments, errorflag)

				! print to screen
			else 
				allocate(optw(nfreemoments,nfreemoments))
				do i = 1, nfreemoments
					do j = 1, nfreemoments, 1
						if (i .eq. j) then
							optw(i,j) = 1.0d0
						else 
							optw(i,j) = 0.0d0
						end if
					end do
				end do
			end if
			! ratio obs for the criterion function
			ratio_obs = dble(nind)/dble(nsim)
			if (ismaster) write(10,'(A, I4)') 'number of obs :', nind
			if (ismaster) write(10,'(A, F9.3)') 'ratio of obs to simulated persons', ratio_obs
			if (ismaster) write(*,'(A, I4)') 'error flag on optimal weight matrix', errorflag

		end subroutine loadmoments

		subroutine estimate
			type (params), allocatable :: initpar(:)

			call initmpi
			if (ismaster) then
				open(10,file='../output/estimation/control_estimate_'//adjustl(trim(scenario))//'.log')
			end if

			! ran by both master and slaves
			call statespace
			call auxiliary
			! initialize for simulations
			call prepare_estimation
			! initialize parameters
			allocate(initpar(npar))
			call initialize(initpar)

			call mpi_barrier(MPI_COMM_WORLD, ier)

			if (ismaster) then
			
				call doestimation(initpar)
				call stop_slaves
			else
				! go to slave job
				call slave_program
			end if

			if (ismaster) then
				close(10)
			end if


			call stopmpi

		end subroutine estimate

		subroutine generate
			type (params), allocatable :: par(:)
			type (rules), allocatable :: optrules(:, :)
			type (sim), allocatable :: pop(:)

			call initmpi
			if (ismaster) then
				open(10,file='../output/simulation/control_simulation_'//adjustl(trim(scenario))//'.log')
			end if
			! ran by both master and slaves
			call statespace
			call auxiliary
			! initialize for simulations
			call prepare_simulation
			! load estimated parameters
			allocate(par(npar))
			call initialize_simulation(par)

			call mpi_barrier(MPI_COMM_WORLD, ier)

			if (ismaster) then
				allocate(optrules(nages, nmaxstates))
				!  solve for decision rules
				call recursion(par, optrules)
				if (isave) then
					call saverules(optrules)
				end if
				! simulate data
				allocate(pop(nsim*nages))
				call simulate(optrules,pop, par)
				call stop_slaves
			else
				! go to slave job
				call slave_program
			end if

			if (ismaster) then
				close(10)
			end if

			call stopmpi

		end subroutine generate

		subroutine target(pop, sim_targ)
				type (sim) pop(nsim*nages)
				double precision sim_targ(2)
		end subroutine

        subroutine extractfreepar(par, freepar)
            type (params) par(npar)
            double precision freepar(nfreepar)
            integer j, k
            k = 1
            do j = 1, npar
                if (par(j)%free) then
                    freepar(k) = par(j)%value
                    k = k + 1
                end if
             end do
         end subroutine extractfreepar

         subroutine setfreepar(par, freepar)
            type (params) par(npar)
            double precision freepar(nfreepar)
            integer j, k
            k = 1
            do j = 1, npar
                if (par(j)%free) then
                     par(j)%value = freepar(k)
                    k = k + 1
                end if
            end do
        end subroutine setfreepar

		subroutine doestimation(initpar)
			type (params) initpar(npar)

			double precision initfreepar(nfreepar),freepar(nfreepar), se(npar), par(npar)
			double precision c(nfreepar), vm(nfreepar), t, eps, rt, fopt, stepfree(nfreepar)
			double precision lbfree(nfreepar), ubfree(nfreepar)
			integer neps, ns, nt, nfcnev, ier, iseed1, iseed2, i, maxevl, iprint,  &
               nacc, nobds, ii, j, jj, n
			logical max
			double precision cum, pvalue, df, iter
			integer ( kind = 4 ) icount
			integer ( kind = 4 ) ifault
			integer ( kind = 4 ) kcount
			integer ( kind = 4 ) konvge
			integer ( kind = 4 ) numres
			real ( kind = 8 ) reqmin
	  		double precision, allocatable :: grid(:,:), temp(:)
	  		double precision gap, minfunc, parmin(npar), ipar(npar), funcvalue
			  integer i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, k
			  integer npt
			  double precision rhobeg, rhoend
		  
			! extract free parameters
			call extractfreepar(initpar, initfreepar)
			g_initpar = initpar
			!  choose which algorithm to use				
            if (alg .eq. 'nm') then
    			write(10,*) '* starting estimation using Nelder-Mead'
    			! nelder-mead
    			reqmin = 1.0D-06
    			konvge = 8
  				kcount = 10000
  				freepar = initfreepar
  				ii = 1
				do i = 1, npar, 1
					if (initpar(i)%free) then
						vm(ii) = initpar(i)%step
						ii = ii + 1
					end if
				end do
    			call nelmin_params ( criterion, nfreepar, initfreepar, freepar, fopt, reqmin, vm, &
    			konvge, kcount, icount, numres, ifault )
    			call setfreepar(initpar,freepar)

			! compute standard errors
			write(10,*) 'termination code for optimization (0 ok) : ', ier
                        write(10,*) ' computing standard errors ...'
				call docovariance(initpar)
			else if (alg .eq. 'po') then
				write(10,*) '* starting estimation using Powell '
    			! nelder-mead
				rhobeg = 0.1d0
    			rhoend = 1.0D-06
  				kcount = 10000
  				freepar = initfreepar
				iprint = 3
				npt = 2*nfreepar + 1
				! set in article of powell
				allocate(temp(10*(npt + nfreepar)))
				call minimize_with_newuoa(calfun, freepar, rhobeg, rhoend, &
							  	npt, iprint, kcount)
				fopt = criterion(freepar)
				call setfreepar(initpar,freepar)
				write(10,*) ' computing standard errors ...'
				call docovariance(initpar)				
            else if (alg .eq. 'sd') then
				write(10,*) ' computing standard errors ...'
				fopt = criterion(initfreepar)
				write(*,*) 'criterion at optimum is : ', fopt
                call docovariance(initpar)
            end if


			! write output
			write(*,*)
			write(10,*) ''
			write(10,*) ''
			write(10,*) '------------------------------------------'
			write(10,'(A20,A9,A9)') 'param','estimate','std.err'
			open(1,file='../output/estimation/estimated-params_'//adjustl(trim(scenario))//'.info')
			do i = 1, npar, 1
				write(1,*) initpar(i)%label, initpar(i)%value,initpar(i)%serr
				if (initpar(i)%free) then
					write(10,'(A20,F12.4,F12.4)') initpar(i)%label, initpar(i)%value,initpar(i)%serr
				else
					write(10,'(A20,F12.4,A9)') initpar(i)%label, initpar(i)%value, '(fixed)'
				end if
			end do
			close(1)
			df = nfreemoments - nfreepar
			write(10,'(A20,F9.3)') 'OID test stat  ', fopt
			write(10,'(A20,I8)') 'degrees of freedom ', int(df)
			call cumchi(fopt,dble(df),cum,pvalue)
			write(10,'(A20,F9.3)') 'pvalue ', pvalue
			write(10,*) '-------------------------------------------'
			write(10,*) ''
			write(10,*) ''


		end subroutine doestimation

		subroutine calfun(theta, h)
			implicit none
			double precision, intent(in)  :: theta(:)
			double precision, intent(out) :: h
			h = criterion(theta)
		end subroutine calfun

		double precision function match(freepar)
			double precision freepar(nfreepar), fval, distance
			type (params) par(npar)
			type (rules) optrules(nages, nmaxstates)
			type (sim) pop(nsim*nages)
			integer i
			double precision sim_e25, sim_medexp
			! save current structural parameters to file
			par(:) = g_initpar(:)
			call setfreepar(par,freepar)
			!  solve for decision rules
			call recursion(par, optrules)
			if (isave) then
				call saverules(optrules)
			end if
			! simulate data
			call simulate(optrules, pop, par)

			! compute moments from simulated dat
			call computetarget(pop(1:nobs), sim_e25, sim_medexp)
			! compare to moments from data
			distance = 0.0d0
			match = dabs(sim_e25 - targ_e25) + dabs(sim_medexp - targ_medexp)
			write(10,*) '* match  : ', match
			write(10,*) 'medexp (sim, target)', sim_medexp, targ_medexp
			write(10,*) 'e25 (sim, target)', sim_e25, targ_e25
			write(10,*) 'at freepars: ', freepar
			write(*,*) '* match  : ', match
			write(*,*) 'medexp (sim, target)', sim_medexp, targ_medexp
			write(*,*) 'e25 (sim, target)', sim_e25, targ_e25
			write(*,*) 'at freepars: ', freepar

		end function
		
		double precision function criterion(freepar)
			double precision freepar(nfreepar), fval, distance(nfreemoments), dd(nfreemoments,1), cc(1,1)
			type (params) par(npar)
			type (rules) optrules(nages, nmaxstates)
			type (sim) pop(nsim*nages)
			integer i
			! save current structural parameters to file
			par(:) = g_initpar(:)
			call setfreepar(par,freepar)
            write(*,*) freepar
			!  solve for decision rules
			call recursion(par, optrules)
			if (isave) then
				call saverules(optrules)
			end if
			! simulate data
			call simulate(optrules, pop, par)

			! compute moments from simulated dat
			call computesimmoments(pop(1:nobs))

			! compare to moments from data
			distance(:) = 0.0d0
			call comparemoments(distance)
			dd(:,1) = distance(:)
			cc =  matmul(matmul(transpose(dd),optw),dd)
			criterion = cc(1,1) * dble(nind) /(1.0d0 + dble(nind)/dble(nsim))
			write(10,*) 'criterion  : ', criterion
			write(10,*) 'at freepars: ', freepar
		end function

		subroutine recursion(par, optrules)
			type (params) par(npar)
			type (rules) optrules(nages, nmaxstates)
			integer a, s, j, k
			integer sent, received, jobsize, index(5)
			integer dum, size
			integer sender, tag, i, nactive, printage(9)
			double precision nextvalue(nmaxstates), thispar(npar), input(npar + nmaxstates)
			double precision, allocatable :: result(:,:)
			printage = (/25,35,45,55,65,75,85,95,105/)
			do j = 1, nmaxstates
				nextvalue(j) = 0.0d0
			end do
			do i = 1, npar, 1
				thispar(i) = par(i)%value
			end do
			dum = 1
			size = npar + nmaxstates
			write(*,*) 'starting recursion ...'
			do a=nages, 1, -1
					nactive = numworkers
					jobsize = floor(dble(nstates(a))/dble(nactive))
					do while (jobsize.lt.nw)
						nactive = nactive - 1
						if (nactive.eq.0) then
							nactive = 1
							exit
						end if
						jobsize = floor(dble(nstates(a))/dble(nactive))
					end do
					if (any(printage.eq.ages(a))) then
						write(*,*) 'at age ', ages(a)
					end if
					! broadcast the value of next period and parameters
					input(1:npar) = thispar
					input(npar + 1: nmaxstates + npar) = nextvalue
					do j = 1, numworkers, 1
						call mpi_send(dum, 1, MPI_INTEGER, j, 1, MPI_COMM_WORLD, ier)
					end do
					call mpi_bcast(input,size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

					! send jobs
					sent = 0
					index(1) = 1
					do j = 1, nactive, 1
						if (j.eq.nactive) then
							index(2) = nstates(a)
							index(3) = nstates(a) - index(1) + 1
						else
							index(2) = index(1) + jobsize - 1
							index(3) = jobsize
						end if
						index(4) = a
						index(5) = 1
						sent = sent + index(3)
						call mpi_send(index, 5, MPI_INTEGER, j, 2, MPI_COMM_WORLD,ier)
						!write(*,*) 'master sent ', index(3), 'to ', j
						call flush(0)
						index(1) = index(2)+1
					end do
					if (sent.lt.nstates(a)) then
						write(*,*) 'WARNING: Did not send all points'
					end if

					! receive jobs and assign to arrays
					received = 0

					do while (received .lt. sent)
						! get first the indices from worker
						call mpi_recv(index, 5, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG,MPI_COMM_WORLD,status,ier)
						! get who sent it
						sender = status(MPI_SOURCE)
						! then get data from that worker
						allocate(result(index(3),10))
						call mpi_recv(result,10*index(3),MPI_DOUBLE_PRECISION,sender, MPI_ANY_TAG, MPI_COMM_WORLD, status, ier)
						received = received + index(3)

						! Assign his data to array of results
						i = 1
						do k = index(1), index(2), 1
							optrules(a,k)%value  = result(i,1)
							optrules(a,k)%cons   = result(i,2)
							optrules(a,k)%medexp = result(i,3)
							optrules(a,k)%claim  = result(i,4)
							optrules(a,k)%work   = result(i,5)
							optrules(a,k)%options(1)   = result(i,6)
							optrules(a,k)%options(2)   = result(i,7)
							optrules(a,k)%options(3)   = result(i,8)
							optrules(a,k)%options(4)   = result(i,9)
							optrules(a,k)%cash = result(i,10)
							i = i + 1
						end do

						deallocate(result)
					end do
					! write value function to nextvalue for next year
					do j = 1, nmaxstates, 1
						nextvalue(j) = optrules(a,j)%value
					end do
			end do


		end subroutine recursion

		subroutine stop_slaves
			integer j, dum
			dum = 1
			! send kill message
			do j = 1, numworkers, 1
				call mpi_send(dum, 1, MPI_INTEGER, j, 9, MPI_COMM_WORLD,ier)
			end do
		end subroutine stop_slaves

		subroutine slave_program
			logical work, flag
			integer index(5), todo, dum, size
			double precision, allocatable :: result(:,:)
			double precision par(npar), input(npar + nmaxstates), nextvalue(nmaxstates)
			type (preferences) prefs
			work = .true.
			size = npar + nmaxstates
			do while (work)
				call mpi_probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ier)
				if (status(MPI_TAG) .eq. 1) then
					call mpi_recv(dum, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, status, ier)
					! retrieve information on parameters
					call mpi_bcast(input,size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
					par = input(1:npar)
					call assign_prefs(par, prefs)
					nextvalue = input(npar + 1: size)
				else if (status(MPI_TAG) .eq. 2) then
					! block receive for the job
					call mpi_recv(index,5,MPI_INTEGER,0,2,MPI_COMM_WORLD,status,ier)
					!perform job
					allocate(result(index(3),10))
					call get_opt(index(4) ,index(1), index(2), index(3), result, nextvalue, prefs)
					call mpi_send(index,5, MPI_INTEGER,0,rank,MPI_COMM_WORLD,ier)
					call mpi_send(result,10*index(3),MPI_DOUBLE_PRECISION,0,rank,MPI_COMM_WORLD,ier)
					deallocate(result)
				else if (status(MPI_TAG) .eq. 9) then
					work = .false.
				end if
			end do
		end subroutine slave_program

		subroutine assign_prefs(par, prefs)
			double precision par(npar)
			type (preferences) prefs
			double precision uumin, good, leisure
			integer z, a, h, hh
			double precision sage
			prefs%sigma = par(1)
			prefs%psi = par(2)
			prefs%phi(1) = par(3)
			prefs%phi(2) = par(4)
			prefs%beta = par(5)
			prefs%alpha(1) = 0.0d0
			prefs%alpha(2) = par(6)
			prefs%alpha(3) = prefs%alpha(2) + dexp(par(7))
			prefs%alpha(4) = prefs%alpha(3) + dexp(par(8))

			prefs%kapa(1) = par(9)
			prefs%kapa(2) = par(10)
			prefs%kapa(3) = 0.0
			prefs%leisure_endow = par(11)

 			! adjust hours penalty to hours
 			prefs%phi(1) = (prefs%leisure_endow - hours_worked)*prefs%phi(1)
 			prefs%phi(2) = (prefs%leisure_endow - hours_worked)*prefs%phi(2)

			do h = 2, nh, 1
				leisure = prefs%leisure_endow
				if (h.eq.3) then
					leisure = leisure - prefs%phi(1)
				else if (h.eq.2) then
					leisure = leisure - prefs%phi(2)
				end if
				good = (cashmin(1,1)**prefs%psi)*(leisure **(1.0d0-prefs%psi))
				prefs%alpha(h) = prefs%alpha(h)*(good**(1.0d0 - prefs%sigma))
			end do

			! assign age parameters
			prefs%age(:) = par(12:14)
			! assign productivity parameters
			prefs%theta(:) = par(15:17)
			prefs%theta2(:) = par(18:20)
			! assign lags parameters
			prefs%gamma(:,1) = par(21:23)
			prefs%gamma(:,2) = par(24:26)
			prefs%gamma(:,3) = par(27:29)

			! assign base
			prefs%base(:) = par(30:32)

			! assign progress
			  prefs%prog_med = par(33)
			  prefs%prog_good = par(34)
			  prefs%prog_vgood = par(35)
      		  prefs%prog_surv = par(36)

      		! recompute lambdas
			do z= 1, nz, 1
				do a = 1, nages, 1
					sage = dble(ages(a))
					do h = 2, nh, 1
						do hh = 2, nh, 1
								lambda(z,a,h-1,hh-1) = dexp(prefs%age(hh-1)*sage  + parh_risk(1,hh-1)*obese(z,a) + &
									parh_risk(2,hh-1)*smoke(z,a) + prefs%gamma(h-1,hh-1) + prefs%base(hh-1) )
						end do
					end do
				end do
			end do
		end subroutine assign_prefs

		subroutine get_opt(a, sbegin, send, ns, rules, nextvalue, prefs)
			integer sbegin, send, ns, s, k, j
			type (preferences) prefs
			!double precision mvalue(nh,ne,nz,nf,nd,ngrid)
			double precision mvalue(ngrid,nh,ne,nf,nd,nz)
			double precision nextvalue(nmaxstates), rules(ns,10)
			double precision inc, nextame,netinc,cash, tr, cstar, mstar, vstar, ben, oopexp
			double precision cguess, mguess, vguess, cashstar
			integer a,h,e,z,f,d,m,w,r,ddd
			integer solution, rstar,dstar, opt
			integer aa, hh, ee, zz,ff, dd, mm, ww, rr, fff

			type (spendargs) myargs

			if (a.ne.nages) then
				! unpack value function
				do j = 1, nstates(a+1), 1
					mm = states(a+1,j,1)
					ww = states(a+1,j,2)
					hh = states(a+1,j,3)
					ee = states(a+1,j,4)
					ff = states(a+1,j,5)
					dd = states(a+1,j,6)
					zz = states(a+1,j,7)
					mvalue((mm-1)*nw + ww,hh,ee,ff,dd,zz) = nextvalue(j)
				end do
			end if


			do s = 1, ns, 1
				! Find state
				k = sbegin + s - 1
				m = states(a,k,1)
				w = states(a,k,2)
				h = states(a,k,3)
				e = states(a,k,4)
				f = states(a,k,5)
				d = states(a,k,6)
				z = states(a,k,7)

				! define type of solution needed
				! no retirement possible
				if (ages(a).lt.rmin) then
					solution = 1
				end if
				! all retired
				if (ages(a).ge.rmax) then
					solution = 2
				end if
				! only retirement possible
				if (ages(a).ge.rmin .and. ages(a).lt.era) then
					solution = 3
				end if
				! SS possible, not medicare eligible
				if (ages(a).ge.era .and. ages(a).lt.mcare) then
					solution = 4
				end if
				! SS possible, medicare eligible
				if (ages(a).ge.mcare .and. ages(a).lt.rmax) then
					solution = 5
				end if

				! arguments constant (state)
				myargs%a = a
				myargs%h = h
				myargs%e = e
				myargs%z = z
				myargs%f = f
				myargs%d = d
				! assign value for each of work x claim value to -infinity
				rules(s,6:10) = -1.0d10

				select case (solution)

					case (1) ! Worker not eligible for pension or SS benefits

						! insurance status will remain same
						ff = f
						! claiming status is same
						dd = d
						! individual does not retire
						rr = 1
						! find earnings
						inc = earn(z,a)*pte(e)
						! update his aime
						nextame = newame(gridame(m),inc,a,z)
						! calculate net income
						inc = inc + spearn(a,z,inc)
						netinc = inctax(inc,rr,a,z)
						! compute cash on hand
						cash = gridw(w) + netinc
						! check if eligible for transfers
						tr = transfer(cash, a, z)
						cash = cash + tr
						! prepare arguments for solving spending allocation
						myargs%cash = cash
						myargs%tr = tr
						myargs%nextame = nextame
						myargs%ff = ff
						myargs%dd = dd
						myargs%rr = rr
						! obtain solution
						call optspend(myargs, mvalue, cstar, mstar, vstar, prefs)
						rules(s,1) = vstar
						rules(s,2) = cstar
						rules(s,3) = mstar
						rules(s,4) = dble(dd)
						rules(s,5) = dble(rr)
						rules(s,6) = vstar
						rules(s,10) = cash

					case (2) ! retired folks
						! do not change insurance status
						ff = 3
						! do not change SS status
						dd = 2
						! is retired
						rr = 2
						! earnings shock defaults to 1
						ee = 1
						! find retirement income (SS)
						inc = 12.0d0*pia(gridame(m))
						inc = inc +  spearn(a,z,inc)
						nextame = gridame(m)
						! calculate net income
						netinc = inctax(inc,rr,a,z)
						! compute cash on hand
						cash = gridw(w) + netinc
						! check if eligible for transfers
						tr = transfer(cash, a, z)
						cash = cash + tr
						if (ages(a).lt.ages(nages)) then
							! prepare arguments for solving spending allocation
							myargs%cash = cash
							myargs%tr = tr
							myargs%nextame = nextame
							myargs%ff = ff
							myargs%dd = dd
							myargs%rr = rr
							! obtain solution
							call optspend(myargs, mvalue, cstar, mstar, vstar, prefs)
							rules(s,1) = vstar
							rules(s,2) = cstar
							rules(s,3) = mstar
							rules(s,4) = dble(dd)
							rules(s,5) = dble(rr)
							rules(s,9) = vstar
							rules(s,10) = cash
						else
							! terminal value at last age
							rules(s,1) = bequest(gridw(w), a, h, prefs)
							rules(s,2)  = 0.0d0
							rules(s,3)   = 0.0d0
							rules(s,4) = dble(dd)
							rules(s,5)= dble(rr)
							rules(s,9) = rules(s,1)
							rules(s,10) = 0.0d0
						end if

					case (3) ! Age 50, need check if decides not to work
						!***
						! if does not work
						! need check whether looses coverage
						if (f.eq.2) then
							ff = 1
						else
							ff = f
						end if
						!! SS status
							dd = d
						!! is retired
							rr = 2
						! spouse income
						inc = spearn(a,z,0.0d0)
						! there is no update of AIME, does not work
						nextame = gridame(m)
						! calculate net income
						netinc = inctax(inc,rr,a,z)
						! compute cash on hand
						cash = gridw(w) + netinc
						! check if eligible for transfers
						tr = transfer(cash, a, z)
						cash = cash + tr
						! prepare arguments for solving spending allocation
						myargs%cash = cash
						cashstar = cash
						myargs%tr = tr
						myargs%nextame = nextame
						myargs%ff = ff
						myargs%dd = dd
						myargs%rr = rr
						! obtain solution
						call optspend(myargs, mvalue, cstar, mstar, vstar, prefs)

						! default to not working
						rstar = rr
						dstar = dd
						opt = (rr-1)*2 + dd
						rules(s,5+opt) = vstar
						! if works

						! keeps his insurance coverage if already has one, gets coverage if uninsured 
						! and takes job
						ff = f
						!! Social Security
						dd = d
						!! is not retired
						rr = 1
						! find earnings
						inc = earn(z,a)*pte(e)
						! update his aime
						nextame = newame(gridame(m),inc,a,z)
						! calculate net income
						inc = inc + spearn(a,z,inc)
						netinc = inctax(inc,rr,a,z)
						! compute cash on hand
						cash = gridw(w) + netinc
						! check if eligible for transfers
						tr = transfer(cash, a, z)
						cash = cash + tr
						! prepare arguments for solving spending allocation
						myargs%cash = cash
						myargs%tr = tr
						myargs%nextame = nextame
						myargs%ff = ff
						myargs%dd = dd
						myargs%rr = rr
						! obtain solution
						call optspend(myargs, mvalue, cguess, mguess, vguess, prefs)
						opt = (rr-1)*2 + dd
						rules(s,5+opt) = vguess
						if (vguess .gt. vstar) then
							vstar = vguess
							cstar = cguess
							mstar = mguess
							rstar = rr
							cashstar = cash
						end if
						! take maximum and assign solution
						rules(s,1) = vstar
						rules(s,2) = cstar
						rules(s,3) = mstar
						rules(s,4) = dble(dd)
						rules(s,5) = dble(rstar)
						rules(s,10) = cashstar

					case (4) ! Age 62, SS possible, not Medicare eligible
						cstar = 0.0d0
						mstar = 0.0d0
						vstar = umin

						! if does not work
						rr = 2

						! need check whether looses coverage
						if (f.eq.2) then
							ff = 1
						else
							ff = f
						end if
						if (ages(a).eq.(mcare-1)) then
							fff = 3
						else
							fff = ff
						end if

						!if does not claims SS
						if (d.eq.1) then
							dd = 1
							! compute income (nothing otherwise)
							inc = spearn(a,z,0.0d0)
							! update aime (no)
							nextame = gridame(m)
							! calculate net income
							netinc = inctax(inc,rr,a,z)
							! compute cash on hand
							cash = gridw(w) + netinc
							! check if eligible for transfers
							tr = transfer(cash, a, z)
							cash = cash + tr
							! prepare arguments for solving spending allocation
							myargs%cash = cash

							myargs%tr = tr
							myargs%nextame = nextame
							myargs%ff = ff
							myargs%dd = dd
							myargs%rr = rr
							myargs%fff = fff
							! obtain solution
							cstar = 0.6d0*myargs%cash
							call optspend(myargs, mvalue, cstar, mstar, vstar, prefs)
							opt = (rr-1)*2 + dd
							rules(s,5+opt) = vstar
							cashstar = cash
							rstar = rr
							dstar = dd
						end if
						! if claims SS
						dd = 2
						if (d.eq.1) then
							! compute income (ssben)
							inc = 12.0d0*ssben(gridame(m),a,z)
							! adjust aime so that pia(aime) = ssben
							nextame = invpia(inc/12.0d0)
						else
							inc = 12.0d0*pia(gridame(m))
							nextame = gridame(m)
						end if
						! calculate net income
						inc = inc + spearn(a,z,inc)
						netinc = inctax(inc,rr,a,z)
						! compute cash on hand
						cash = gridw(w) + netinc
						! check if eligible for transfers
						tr = transfer(cash, a, z)
						cash = cash + tr
						! prepare arguments for solving spending allocation
						myargs%cash = cash
						myargs%tr = tr
						myargs%nextame = nextame
						myargs%ff = ff
						myargs%fff = fff
						myargs%dd = dd
						myargs%rr = rr
						! obtain solution
						if (d .eq. 1) then
							cguess = cstar
						else
							cguess = 0.8d0*myargs%cash
						end if		
						call optspend(myargs, mvalue, cguess, mguess, vguess, prefs)
						opt = (rr-1)*2 + dd
						rules(s,5+opt) = vguess

						if (vguess .gt. vstar) then
							vstar = vguess
							cstar = cguess
							mstar = mguess
							dstar = dd
							rstar = rr
							cashstar = cash
						end if
						! if does works
						rr = 1
						if (f .eq. 1) then
							ff = 2
						else 
							ff = f
						end if		
						if (ages(a).eq.(mcare-1)) then
							fff = 3
						else
							fff = ff
						end if

						! if does not claim SS
						if (d.eq.1) then
							dd = 1
							! compute income ( earnings)
							inc = earn(z,a)*pte(e)
							! update aime
							nextame = newame(gridame(m),inc,a,z)
							! calculate net income
							inc = inc + spearn(a,z,inc)
							netinc = inctax(inc,rr,a,z)
							! compute cash on hand
							cash = gridw(w) + netinc
							! check if eligible for transfers
							tr = transfer(cash, a, z)
							cash = cash + tr
							! prepare arguments for solving spending allocation
							myargs%cash = cash
							myargs%tr = tr
							myargs%nextame = nextame
							myargs%ff = ff
							myargs%dd = dd
							myargs%rr = rr
							myargs%fff = fff
							! obtain solution
							cguess = cstar
							call optspend(myargs, mvalue, cguess, mguess, vguess, prefs)
							opt = (rr-1)*2 + dd
							rules(s,5+opt) = vguess

							if (vguess .gt. vstar) then
								vstar = vguess
								cstar = cguess
								mstar = mguess
								dstar = dd
								rstar = rr
								cashstar = cash
							end if
						end if

						! if claims SS
						dd = 2
						! compute income (ssben and earnings)
						inc = earn(z,a)*pte(e)
						! benefits
						if (d.eq.1) then
							! compute income (ssben)
							ben = 12.0d0*ssben(gridame(m),a,z)
							! adjust aime so that pia(aime) = ssben
							nextame = invpia(ben/12.0d0)
						else
							ben = 12.0d0*pia(gridame(m))
							nextame = gridame(m)
						end if
						! Apply earnings test (if applicable)
						inc = inc + earntest(ben,inc,a)
						inc = inc + spearn(a,z,inc)
						! calculate net income
						netinc = inctax(inc,rr,a,z)
						! compute cash on hand
						cash = gridw(w) + netinc
						! check if eligible for transfers
						tr = transfer(cash, a, z)
						cash = cash + tr
						! prepare arguments for solving spending allocation
						myargs%cash = cash
						myargs%tr = tr
						myargs%nextame = nextame
						myargs%ff = ff
						myargs%dd = dd
						myargs%rr = rr
						myargs%fff = fff
						! obtain solution
						cguess = cstar
						call optspend(myargs, mvalue, cguess, mguess, vguess, prefs)
						opt = (rr-1)*2 + dd
						rules(s,5+opt) = vguess

						if (vguess .gt. vstar) then
							vstar = vguess
							cstar = cguess
							mstar = mguess
							dstar = dd
							rstar = rr
							cashstar = cash
						end if

						! take maximum
						rules(s,1) = vstar
						rules(s,2) = cstar
						rules(s,3) = mstar
						rules(s,4) = dble(dstar)
						rules(s,5) = dble(rstar)
						rules(s,10) = cashstar

					case (5) ! Age 65, SS claiming and Medicare eligible

						ff = nf

							! if does not work
							rr = 2
								vstar = -10.d6
								!if does not claims SS
								if (d.eq.1 .and. ages(a).lt.(rmax-1)) then
									dd = 1
									! compute income (nothing otherwise)
									inc = spearn(a,z,0.0d0)
									! update aime (no)
									nextame = gridame(m)
									! calculate net income
									netinc = inctax(inc,rr,a,z)
									! compute cash on hand
									cash = gridw(w) + netinc
									! check if eligible for transfers
									tr = transfer(cash, a, z)
									cash = cash + tr

									! prepare arguments for solving spending allocation
									myargs%cash = cash
									myargs%tr = tr
									myargs%nextame = nextame
									myargs%ff = ff
									myargs%dd = dd
									myargs%rr = rr
									! obtain solution
									cguess = 0.6d0*myargs%cash
									call optspend(myargs, mvalue, cguess, mguess, vguess, prefs)
									opt = (rr-1)*2 + dd
									rules(s,5+opt) = vguess

									if (vguess .gt. vstar) then
										vstar = vguess
										cstar = cguess
										mstar = mguess
										dstar = dd
										rstar = rr
										cashstar = cash
									end if

								end if

								! if claims SS
								dd = 2

									! benefits
									if (d.eq.1) then
										! compute income (ssben)
										ben = 12.0d0*ssben(gridame(m),a,z)
										! adjust aime so that pia(aime) = ssben
										nextame = invpia(ben/12.0d0)
									else
										ben = 12.0d0*pia(gridame(m))
										nextame = gridame(m)
									end if
									inc = ben + spearn(a,z,ben)
									! calculate net income
									netinc = inctax(inc,rr,a,z)
									! compute cash on hand
									cash = gridw(w) + netinc

									! check if eligible for transfers
									tr = transfer(cash, a, z)
									cash = cash + tr
									! prepare arguments for solving spending allocation
									myargs%cash = cash
									myargs%tr = tr
									myargs%nextame = nextame
									myargs%ff = ff
									myargs%dd = dd
									myargs%rr = rr
									! obtain solution
									cguess = cstar
									call optspend(myargs, mvalue, cguess, mguess, vguess, prefs)
									opt = (rr-1)*2 + dd
									rules(s,5+opt) = vguess
									if (vguess .gt. vstar) then
										vstar = vguess
										cstar = cguess
										mstar = mguess
										dstar = dd
										rstar = rr
										cashstar = cash
									end if



							! if does works
							rr = 1
								! if does not claim SS
								if (d.eq.1 .and. ages(a).lt.(rmax-1)) then
									dd = 1
									! earnings only
									inc =  earn(z,a)*pte(e)
									nextame = newame(gridame(m),inc,a,z)
									! calculate net income
									inc = inc + spearn(a,z,inc)
									netinc = inctax(inc,rr,a,z)
									! compute cash on hand
									cash = gridw(w) + netinc
									! check if eligible for transfers
									tr = transfer(cash, a, z)
									cash = cash + tr
									! prepare arguments for solving spending allocation
									myargs%cash = cash
									myargs%tr = tr
									myargs%nextame = nextame
									myargs%ff = ff
									myargs%dd = dd
									myargs%rr = rr
									! obtain solution
									cguess = cstar
									call optspend(myargs, mvalue, cguess, mguess, vguess, prefs)
									opt = (rr-1)*2 + dd
									rules(s,5+opt) = vguess
									if (vguess .gt. vstar) then
										vstar = vguess
										cstar = cguess
										mstar = mguess
										dstar = dd
										rstar = rr
										cashstar = cash
									end if
								end if

								! if claims SS
								dd = 2
									! compute income (ssben and earnings)
									inc = earn(z,a)*pte(e)
									! benefits
									if (d.eq.1) then
										! compute income (ssben)
										ben = 12.0d0*ssben(gridame(m),a,z)
										! adjust aime so that pia(aime) = ssben
										nextame = invpia(ben/12.0d0)
									else
										ben = 12.0d0*pia(gridame(m))
										nextame = gridame(m)
									end if
									inc = inc + earntest(ben,inc,a)
									inc = inc + spearn(a,z,inc)
									! calculate net income
									netinc = inctax(inc,rr,a,z)
									! compute cash on hand
									cash = gridw(w) + netinc
									! check if eligible for transfers
									tr = transfer(cash, a, z)
									cash = cash + tr

									! prepare arguments for solving spending allocation
									myargs%cash = cash
									myargs%tr = tr
									myargs%nextame = nextame
									myargs%ff = ff
									myargs%dd = dd
									myargs%rr = rr
									! obtain solution
									cguess = cstar
									call optspend(myargs, mvalue, cguess, mguess, vguess, prefs)
									opt = (rr-1)*2 + dd
									rules(s,5+opt) = vguess
									if (vguess .gt. vstar) then
										vstar = vguess
										cstar = cguess
										mstar = mguess
										dstar = dd
										rstar = rr
										cashstar = cash
									end if

						! take maximum
						rules(s,1) = vstar
						rules(s,2) = cstar
						rules(s,3) = mstar
						rules(s,4) = dble(dstar)
						rules(s,5) = dble(rstar)
						rules(s,10) = cashstar

	
				end select

			end do

		end subroutine get_opt

		subroutine optspend(myargs, mvalue, cstar, mstar, vstar, prefs)
			! arguments of subroutine
			type (spendargs) myargs
			type (preferences) prefs
			double precision mvalue(ngrid,nh,ne,nf,nd,nz)
			double precision cstar, mstar, vstar
			! initial values
		  	double precision ax, bx, cx, x0, x1, x2, x3, f0, f1, f2, f3, q1, q2
		  	double precision r, c, tol
	      	parameter (tol = 10.0d0, r= 0.61803399d0, c= 1.0d0 - r)
			double precision par(2), value
			double precision mmax
			! find bounds on consumption and medexp
			myargs%cmax = myargs%cash
			myargs%cmin = 0.0d0
			myargs%mmax = invoop(myargs%cash - myargs%cmin,myargs%a,myargs%ff, myargs%z, myargs%tr)
			if (myargs%mmax .gt. medmax) then
				myargs%mmax = medmax
			end if
			myargs%mmin = 0.0D0
			ax = myargs%cmin
			cx = myargs%cmax
			if (cstar .gt. ax .and. cstar .lt. cx) then
				bx = cstar
			else	
				bx = (cx + ax) * 0.6d0
			end if
			x0 = ax
			x3 = cx
			if (abs(cx-bx).gt.abs(bx-ax)) then
			    x1 = bx
			    x2 = bx + c*(cx-bx)
			else
			    x2 = bx
			    x1 = bx - c*(bx-ax)
			end if
			call solve_m(x1, f1, myargs, mvalue, prefs, q1)
			call solve_m(x2, f2, myargs, mvalue, prefs, q2)
			do while (abs(x3-x0).gt.10.0d0)
			    if (f2.lt.f1) then
			          x0 = x1
			          x1 = x2
			          x2 = r*x1 + c*x3
			          f0 = f1
			          f1 = f2
			          call solve_m(x2, f2, myargs, mvalue, prefs, q2)
			    else
			          x3 = x2
			          x2 = x1
			          x1 = r*x2 + c*x0
			          f3 = f2
			          f2 = f1
			          call solve_m(x1, f1, myargs, mvalue, prefs, q1)
			    end if
			end do
			if (f1.lt.f2) then
			    cstar = x1
			    mstar = q1
			    vstar = -f1
			else
			    cstar = x2
			    mstar = q2
			    vstar = -f2
			end if
		end subroutine optspend

		subroutine solve_m(cons, value, myargs, mvalue, prefs, medexp)
			type (spendargs) myargs
			type (preferences) prefs
			!double precision mvalue(nh,ne,nz,nf,nd,ngrid)
			double precision mvalue(ngrid,nh,ne,nf,nd,nz)
			double precision medexp, cons, value, spend(2)
			! initial values
		  	double precision ax, bx, cx, x0, x1, x2, x3, f0, f1, f2, f3
		  	double precision r, c, tol
	      	parameter (tol = 10.0d0, r= 0.61803399d0, c= 1.0d0 - r)
			ax = 0.0d0
			cx = invoop(myargs%cash - cons,myargs%a,myargs%ff, myargs%z, myargs%tr)
			if (cx .gt. myargs%mmax) then
				cx = myargs%mmax
			end if
			bx = (cx + ax) * 0.33d0
			x0 = ax
			x3 = cx
			if (abs(cx-bx).gt.abs(bx-ax)) then
			    x1 = bx
			    x2 = bx + c*(cx-bx)
			else
			    x2 = bx
			    x1 = bx - c*(bx-ax)
			end if
			spend(1) = cons
			spend(2) = x1
			f1 = negvalue(spend, myargs, mvalue, prefs)
			spend(2) = x2
			f2 = negvalue(spend, myargs, mvalue, prefs)
			do while (abs(x3-x0).gt.10.0d0)
			    if (f2.lt.f1) then
			          x0 = x1
			          x1 = x2
			          x2 = r*x1 + c*x3
			          f0 = f1
			          f1 = f2
			          spend(2) = x2
					  f2 = negvalue(spend, myargs, mvalue, prefs)
			    else
			          x3 = x2
			          x2 = x1
			          x1 = r*x2 + c*x0
			          f3 = f2
			          f2 = f1
			          spend(2) = x1
					  f1 = negvalue(spend, myargs, mvalue, prefs)
			    end if
			end do

			if (f1.lt.f2) then
			    medexp = x1
			    value = f1
			else
			    medexp = x2
			    value = f2
			end if
		end subroutine solve_m


		double precision function negvalue(spend, myargs, mvalue, prefs)
			double precision ctry, mtry, spend(2)
			type (spendargs) myargs
			type (preferences) prefs
			double precision mvalue(ngrid,nh,ne,nf,nd,nz)
			integer fff
			double precision nextw, wu, probh(nh)
			double precision v(2,2), ev, pv, ameu, mmax
			integer w0, w1, ame0, ame1, hh, ee
			integer q00, q01, q10, q11
			! assign current spending variables
			ctry = spend(1)
			mtry = spend(2)
			! transfer arguments
			if (ages(myargs%a).eq.(mcare - 1)) then
				fff = myargs%fff
			else
				fff = myargs%ff
			end if
			! compute end of period assets
			nextw = (1.0d0+rrate)*(myargs%cash - ctry - oop(mtry,myargs%a,myargs%ff, myargs%z, myargs%tr))
			if (nextw .lt. wmin) then
				nextw = wmin
			end if
			if (nextw .gt. wmax) then
				nextw = wmax
			end if
			call scale(dsqrt(nextw),c_wmin,c_wmax,nw,c_gridw,w0,w1,wu)
			! scale for interpolation ame
			call scale(myargs%nextame,amemin,amemax,name,gridame,ame0,ame1,ameu)
			! get probability of health shock, depends on spending
			call getprobh(myargs%z,myargs%a,myargs%h,mtry,probh,prefs)

			ev = 0.0d0
			q00 = (ame0-1)*nw + w0
			q01 = (ame0-1)*nw + w1
			q10 = (ame1-1)*nw + w0
			q11 = (ame1-1)*nw + w1
			if 	(ages(myargs%a).lt.(rmax-1)) then
				do ee = 1, ne
					do hh = 2, nh, 1
						v(1,1) = mvalue(q00,hh,ee,fff,myargs%dd,myargs%z)
						v(1,2) = mvalue(q01,hh,ee,fff,myargs%dd,myargs%z)
						v(2,1) = mvalue(q10,hh,ee,fff,myargs%dd,myargs%z)
						v(2,2) = mvalue(q11,hh,ee,fff,myargs%dd,myargs%z)
						pv = v(1,1) + ameu*(-v(1,1)+v(2,1)) + wu*(-v(1,1)+v(1,2)) &
							+ ameu*wu*(v(1,1) - v(2,1) - v(1,2) + v(2,2))
						ev = ev + probe(myargs%e,ee)*probh(hh)*pv
					end do
				end do

			else
				ee = 1
				do hh = 2, nh, 1
					v(1,1) = mvalue(q00,hh,ee,fff,myargs%dd,myargs%z)
					v(1,2) = mvalue(q01,hh,ee,fff,myargs%dd,myargs%z)
					v(2,1) = mvalue(q10,hh,ee,fff,myargs%dd,myargs%z)
					v(2,2) = mvalue(q11,hh,ee,fff,myargs%dd,myargs%z)
					pv = v(1,1) + ameu*(-v(1,1)+v(2,1)) + wu*(-v(1,1)+v(1,2)) &
							+ ameu*wu*(v(1,1) - v(2,1) - v(1,2) + v(2,2))
					ev = ev + probh(hh)*pv
				end do
			end if

			! bequest
			ev = ev + probh(1)*bequest(nextw, myargs%a,myargs%h, prefs)

			! get utility
			negvalue = utility(ctry,myargs%rr,myargs%h,myargs%a,myargs%z,mtry, prefs) + prefs%beta*ev
			! negative since we minimize
			negvalue = - negvalue
		end function

		subroutine saverules(optrules)
			type (rules) optrules(nages, nmaxstates)
			double precision cons(nw), medexp(nw)
			character*3 agee
			character*100 fmt
			integer j, i, a
			fmt = '(I2,A,I2,A,I2,A,I2,A,I2,A,I2,A,I2,A,F12.3,A,F12.3,A,F12.3,A,I2,A,I2,A,F16.9,A,F12.3)'
			if (scenario .eq. 'reference') then
				! saving sizes
				do i = 1, nsave, 1
					do a = 1, nages, 1
						if (ages(a) .eq. saveages(i)) then
							write(agee,'(I2)') ages(a)
							open(1,file='../output/rules/rules_age'//trim(agee)//'_'//adjustl(trim(scenario))//'.csv')
								do j = 1, nstates(a), 1
									write(1,fmt) states(a,j,3),',',states(a,j,4),',',states(a,j,7),',', &
										states(a,j,5),',',states(a,j,6),',',states(a,j,1),',', &
										states(a,j,2),',',gridw(states(a,j,2)),',', optrules(a,j)%cons,',',optrules(a,j)%medexp,',',&
										int(optrules(a,j)%work),',',int(optrules(a,j)%claim),',',optrules(a,j)%value, ',', optrules(a,j)%cash
								end do
							close(1)
						end if
					end do
				end do
			end if
		end subroutine saverules

		subroutine simulate(optrules, pop, par)
			type (rules) optrules(nages,nmaxstates), optr(2,2)
			type (sim) pop(nsim*nages)
			type (preferences) prefs
			type (params) par(npar)
			double precision pars(npar)
			integer a,i, j , n, z
			character*120 fmt
			integer fsim, esim, hsim, dead, dsim, zsim, rsim, lastdsim, lasthsim
			double precision vsim, wsim, csim, oopsim, amesim, msim,incsim, rrsim, ddsim
			double precision logearn, ben, diff, draw,  vprobe(ne),cumprob, pi
			integer w0,w1,m0,m1,s00,s01,s10,s11, solution
			double precision mu,wu, cash, tr, probh(nh)
			character*1 scn1
			character*1 scn2
			character*1 scn3
			character*1 scn4
			integer hh, ee, zz, ff, dd, mm, ww
			double precision options(4)
			integer choice, loc
			! open file to store simulated data
			pars(:) = par(:)%value
			call assign_prefs(pars,prefs)

			! initial population
			open(1,file='../params/input/initsample.csv')
			nobs = 0
			do i = 1, nsim, 1

				! buffer (logriearn rskill insured ass health )
				read(1,*) logearn, zsim, fsim, wsim, hsim

				! assign cohort
				if (nz .eq. 1) then
					zsim = 1
				else
					pi = 1.0d0/dble(nz)
					do z = 1, nz, 1
						if (idraws(i) .le. pi) then
							zsim = z
							exit
						else
							pi = pi + 1.0d0/dble(nz) 
						end if
					end do 
				end if


				if (wsim .lt. wmin) then
					wsim = wmin
				end if
				if (wsim .gt. wmax) then
					wsim = wmax
				end if
				! find out where fits on error term grid
				diff = dexp(logearn - dlog(earn(zsim,1)))
				call point(diff,pte, ne, esim)

				! remaining states fixed.
				dsim = 1
				amesim = 0.0d0

				! fill in optimal solution and update state-space shocks
				do a = 1, nages, 1

				    if (ages(a).eq.mcare) then
				    	fsim = 3
				    end if
				    ! automatic transitions
				    if (ages(a).ge.rmax) then
				    	esim = 1
				    end if

				    ! find out where we are in state space for m and w
				    call scale(amesim,amemin,amemax,name,gridame,m0,m1,mu)
					call scale(dsqrt(wsim),c_wmin,c_wmax,nw,c_gridw,w0,w1,wu)
					loc = getindex(a,hsim,esim,zsim,fsim,dsim,1,1)
				    optr(1,1) = optrules(a,loc + (m0-1)*nw + w0)
				    optr(1,2) = optrules(a,loc + (m0-1)*nw + w1)
				    optr(2,1) = optrules(a,loc + (m1-1)*nw + w0)
				    optr(2,2) = optrules(a,loc + (m1-1)*nw + w1)

				    ! impute values


				    ! consumption
					call blend_102(mu, wu, optr(:,:)%cons, csim)
				    ! medical exp
					call blend_102(mu, wu, optr(:,:)%medexp, msim)
					if (msim .gt. medmax) then
						msim = medmax
					end if
					if (msim .lt. 0.0d0) then
						msim = 0.0D0
					end if
					! Solution for discrete decisions
					if (ages(a).lt.rmin) then
						solution = 1
					end if
					! all retired
					if (ages(a).ge.rmax) then
						solution = 2
					end if
					! only retirement possible
					if (ages(a).ge.rmin .and. ages(a).lt.era) then
						solution = 3
					end if
					! SS possible, not medicare eligible
					if (ages(a).ge.era .and. ages(a).lt.mcare) then
						solution = 4
					end if
					! SS possible, medicare eligible
					if (ages(a).ge.mcare .and. ages(a).lt.rmax) then
						solution = 5
					end if
					! choices
					! 1 = work, no claim
					! 2 = work, claim
					! 3 = no work, no claim
					! 4 = no work, claim

					select case (solution)
					! several cases : a)
						case (1) ! Worker not eligible for pension or SS benefits

							rsim = 1
							dsim = 1
							choice = (rsim - 1)*nd + dsim
							call blend_102(mu, wu, optr(:,:)%options(choice), vsim)
							! insurance status will remain same
							fsim = fsim

							! find earnings
							incsim = earn(zsim,a)*pte(esim)
							! update his aime
							amesim = newame(amesim,incsim,a,zsim)
							if (amesim.gt.amemax) then
								amesim = amemax
							end if

							! calculate net income
							incsim = incsim + spearn(a,zsim,incsim)
							incsim = inctax(incsim,rsim,a,zsim)

							! compute cash on hand
							cash = wsim + incsim

							! check if eligible for transfers
							tr = transfer(cash,a, zsim)
							cash = cash + tr
							! out-of-pocket
							oopsim = oop(msim,a,fsim, zsim, tr)


							! update wealth
							wsim = (1.0d0+rrate)*(cash - csim - oopsim)
							if (wsim.gt.wmax) then
								wsim = wmax
							end if
							if (wsim.lt.wmin) then
								wsim = wmin
							end if

							! next period's earnings shock
							if (ages(a) .lt. rmax) then
							vprobe = cumprobe(esim,:)
								do j = 1, ne, 1
									if (edraws(i,a).lt.vprobe(j)) then
										esim = j
										exit
									end if
								end do
							end if

						case (2) ! retired folks

							rsim = 2
							dsim = 2
							choice = (rsim - 1)*nd + dsim
							call blend_102(mu, wu, optr(:,:)%options(choice), vsim)
							! do not change insurance status
							fsim = fsim

							! find retirement income (SS)
							incsim = 12.0d0*pia(amesim)
							incsim = incsim + spearn(a,zsim, incsim)

							! calculate net income
							incsim = inctax(incsim, rsim,a,zsim)

							! compute cash on hand
							cash = wsim + incsim

							! check if eligible for transfers
							tr = transfer(cash,a, zsim)
							cash = cash + tr

							! out-of-pocket
							oopsim = oop(msim,a,fsim, zsim, tr)

							! update wealth
							wsim = (1.0d0+rrate)*(cash - csim - oopsim)
							if (wsim.gt.wmax) then
								wsim = wmax
							end if
							if (wsim.lt.wmin) then
								wsim = wmin
							end if

							esim = 1

						case (3) ! Age 55, need check if decides not to work

							call blend_102(mu, wu, optr(:,:)%options(1), options(1))
							call blend_102(mu, wu, optr(:,:)%options(3), options(3))
						    ! work
						    if (options(3) .gt. options(1)) then
						    	rsim = 2
						    	vsim = options(3)
						    else
						    	rsim = 1
						    	vsim = options(1)
						    end if
						    ! claim
						    dsim = 1

							if (rsim.eq.2) then
								! need check whether looses coverage
								if (fsim.eq.2) then
									fsim = 1
								else
									fsim = fsim
								end if

								incsim = spearn(a,zsim,0.0d0)

							else
								! insurance status will remain same if already has it
								! will get tied coverage if works and did not have coverage
								if (fsim .eq. 1) then
									fsim = 2
								else 
									fsim = fsim
								end if

								! find earnings
								incsim = earn(zsim,a)*pte(esim)

								! update his aime
								amesim = newame(amesim,incsim,a,zsim)
								if (amesim.gt.amemax) then
									amesim = amemax
								end if

								! calculate gross income
								incsim = incsim + spearn(a,zsim,incsim)

								! next period's earnings shock
								vprobe = cumprobe(esim,:)
								do j = 1, ne, 1
									if (edraws(i,a).lt.vprobe(j)) then
										esim = j
										exit
									end if
								end do

							end if

							! calculate net income
							incsim = inctax(incsim, rsim,a,zsim)

							! compute cash on hand
							cash = wsim + incsim

							! check if eligible for transfers
							tr = transfer(cash, a, zsim)
							cash = cash + tr

							! out-of-pocket
							oopsim = oop(msim,a,fsim, zsim, tr)

							! update wealth
							wsim = (1.0d0+rrate)*(cash - csim - oopsim)
							if (wsim.gt.wmax) then
								wsim = wmax
							end if
							if (wsim.lt.wmin) then
								wsim = wmin
							end if


						case (4) ! Age 62, SS possible, not Medicare eligible

						    lastdsim = dsim
					    	if (lastdsim.eq.1) then
								! claiming
								call blend_102(mu, wu, optr(:,:)%options(1), options(1))
								call blend_102(mu, wu, optr(:,:)%options(2), options(2))
								call blend_102(mu, wu, optr(:,:)%options(3), options(3))
								call blend_102(mu, wu, optr(:,:)%options(4), options(4))
								choice = maxloc(options(:),dim=1)
								if (choice .eq. 1) then
									dsim = 1
									rsim = 1
									vsim = options(choice)
								else if (choice .eq. 2) then
									dsim = 2
									rsim = 1
									vsim = options(choice)
								else if (choice .eq. 3) then
									dsim = 1
									rsim = 2
									vsim = options(choice)
								else if (choice .eq. 4) then
									dsim = 2
									rsim = 2
									vsim = options(choice)
								end if
							else
								dsim = 2
								call blend_102(mu, wu, optr(:,:)%options(2), options(2))
								call blend_102(mu, wu, optr(:,:)%options(4), options(4))
								if (options(4) .gt. options(2)) then
									rsim = 2
									vsim = options(4)
								else
									rsim = 1
									vsim = options(2)
								end if
							end if

							if (dsim.eq.1) then
								if (rsim.eq.1) then

									! compute income ( earnings)
									incsim = earn(zsim,a)*pte(esim)
									! update aime
									amesim = newame(amesim,incsim,a,zsim)
									! calculate net income
									incsim = incsim + spearn(a,zsim,incsim)

									! next period's earnings shock
									vprobe = cumprobe(esim,:)
									do j = 1, ne, 1
										if (edraws(i,a).lt.vprobe(j)) then
											esim = j
											exit
										end if
									end do


								else
									! check if losses coverage
									if (fsim.eq.2) then
										fsim = 1
									else
										fsim = fsim
									end if

									! compute income (nothing otherwise)
									incsim = spearn(a,zsim,0.0d0)

								end if
							else
								if (rsim.eq.1) then
									! compute income (ssben and earnings)
									incsim = earn(zsim,a)*pte(esim)

									! benefits
									if (lastdsim.eq.1) then
										! compute income (ssben)
										ben = 12.0d0*ssben(amesim,a,zsim)
										! adjust aime so that pia(aime) = ssben
										amesim = invpia(ben/12.0d0)
									else
										ben = 12.0d0*pia(amesim)
									end if
									! Apply earnings test (if applicable)
									incsim = incsim + earntest(ben,incsim,a)
									incsim = incsim + spearn(a,zsim,incsim)

									! next period's earnings shock
									vprobe = cumprobe(esim,:)
									do j = 1, ne, 1
										if (edraws(i,a).lt.vprobe(j)) then
											esim = j
											exit
										end if
									end do

								else
									! check if looses coverage
									if (fsim.eq.2) then
										fsim = 1
									else
										fsim = fsim
									end if

									if (lastdsim.eq.1) then
										! compute income (ssben)
										incsim = 12.0d0*ssben(amesim,a,zsim)
										! adjust aime so that pia(aime) = ssben
										amesim = invpia(incsim/12.0d0)
									else
										incsim = 12.0d0*pia(amesim)
									end if

									! calculate net income
									incsim = incsim + spearn(a,zsim,incsim)

								end if

							end if


							! calculate net income
							incsim = inctax(incsim, rsim,a,zsim)

							! compute cash on hand
							cash = wsim + incsim

							! check if eligible for transfers
							tr = transfer(cash, a, zsim)
							cash = cash + tr

							! out-of-pocket
							oopsim = oop(msim,a,fsim, zsim, tr)

							! update wealth
							wsim = (1.0d0+rrate)*(cash - csim - oopsim)
							if (wsim.gt.wmax) then
								wsim = wmax
							end if
							if (wsim.lt.wmin) then
								wsim = wmin
							end if


						case (5) ! Age 65, SS claiming and Medicare eligible

					    	if (lastdsim.eq.1) then
								! claiming
								call blend_102(mu, wu, optr(:,:)%options(1), options(1))
								call blend_102(mu, wu, optr(:,:)%options(2), options(2))
								call blend_102(mu, wu, optr(:,:)%options(3), options(3))
								call blend_102(mu, wu, optr(:,:)%options(4), options(4))
								choice = maxloc(options(:),dim=1)
								if (choice .eq. 1) then
									dsim = 1
									rsim = 1
									vsim = options(choice)
								else if (choice .eq. 2) then
									dsim = 2
									rsim = 1
									vsim = options(choice)
								else if (choice .eq. 3) then
									dsim = 1
									rsim = 2
									vsim = options(choice)
								else if (choice .eq. 4) then
									dsim = 2
									rsim = 2
									vsim = options(choice)
								end if
							else
								dsim = 2
								call blend_102(mu, wu, optr(:,:)%options(2), options(2))
								call blend_102(mu, wu, optr(:,:)%options(4), options(4))
								if (options(4) .gt. options(2)) then
									rsim = 2
									vsim = options(4)
								else
									rsim = 1
									vsim = options(2)
								end if
							end if

							if (dsim.eq.1) then

								if (rsim.eq.1) then
									! earnings only
									incsim =  earn(zsim,a)*pte(esim)
									amesim = newame(amesim,incsim,a,zsim)
									! calculate net income
									incsim = incsim + spearn(a,zsim,incsim)
									! next period's earnings shock
									vprobe = cumprobe(esim,:)
									do j = 1, ne, 1
										if (edraws(i,a).lt.vprobe(j)) then
											esim = j
											exit
										end if
									end do

								else

									! compute income (nothing otherwise)
									incsim = spearn(a,zsim,0.0d0)

								end if
							else

								if (rsim.eq.1) then

									! compute income (ssben and earnings)
									incsim = earn(zsim,a)*pte(esim)

									! benefits
									if (lastdsim.eq.1) then
										! compute income (ssben)
										ben = 12.0d0*ssben(amesim,a,zsim)
										! adjust aime so that pia(aime) = ssben
										amesim = invpia(ben/12.0d0)
									else
										ben = 12.0d0*pia(amesim)
									end if
									! Apply earnings test (if applicable)
									incsim = incsim + earntest(ben,incsim,a)
									incsim = incsim + spearn(a,zsim, incsim)

									! next period's earnings shock
									vprobe = cumprobe(esim,:)
									do j = 1, ne, 1
										!write(*,*) 'shock prob ', s, j, draw, vprob(j)
										if (edraws(i,a).lt.vprobe(j)) then
											esim = j
											exit
										end if
									end do

								else

									! benefits
									if (lastdsim.eq.1) then
										! compute income (ssben)
										ben = 12.0d0*ssben(amesim,a,zsim)
										! adjust aime so that pia(aime) = ssben
										amesim = invpia(ben/12.0d0)
									else
										ben = 12.0d0*pia(amesim)
									end if
									incsim = ben + spearn(a,zsim,ben)

								end if

							end if

							! calculate net income
							incsim = inctax(incsim, rsim,a,zsim)

							! compute cash on hand
							cash = wsim + incsim

							! check if eligible for transfers
							tr = transfer(cash, a, zsim)
							cash = cash + tr

							! out-of-pocket
							oopsim = oop(msim,a,fsim, zsim, tr)

							! update wealth
							wsim = (1.0d0+rrate)*(cash - csim - oopsim)
							if (wsim.gt.wmax) then
								wsim = wmax
							end if
							if (wsim.lt.wmin) then
								wsim = wmin
							end if

					end select

					! Updating Health
					call getprobh(zsim,a,hsim,msim,probh,prefs)
					lasthsim = hsim
					cumprob = 0.0d0
					do j = 1, nh, 1
						cumprob = cumprob + probh(j)
						if (hdraws(i,a).lt.cumprob) then
							hsim = j
							exit
						end if
					end do

					if (hsim.eq.1) then
						dead = 2
					else
						dead = 1
					end if

					! save data
					nobs = nobs + 1
					pop(nobs)%id = i
					pop(nobs)%age = ages(a)
					pop(nobs)%z = gridz(zsim)
					pop(nobs)%wealth = wsim
					pop(nobs)%work = rsim
					pop(nobs)%claim = dsim
					pop(nobs)%income = incsim
					pop(nobs)%insurance = fsim
					pop(nobs)%ame = amesim
					pop(nobs)%consumption = csim
					pop(nobs)%dead = dead
					pop(nobs)%health = hsim
					pop(nobs)%base_health = lasthsim
					pop(nobs)%oop = oopsim
					pop(nobs)%medexp = msim
					pop(nobs)%value = vsim
					pop(nobs)%tr = tr

					if (dead.eq.2) then
						exit
					end if

				end do

			end do

			close(1)

			open(1,file='../output/simulation/simulated_cohort_'//adjustl(trim(scenario))//'.csv')
			fmt = '(I8,A,I3,A,I4,A,F12.3,A,I1,A,I1,A,F12.3,A,I1,A,F12.3,A,F12.3,A,I1,A,I1,A,I1,A,F12.3,A,F12.3,A,F16.9,A,F12.3)'
			do i = 1, nobs, 1
				write(1,fmt)  pop(i)%id, ',', &
							pop(i)%age, ',', &
							pop(i)%z, ',', &
							pop(i)%wealth, ',', &
							pop(i)%work, ',', &
							pop(i)%claim, ',', &
							pop(i)%income, ',', &
							pop(i)%insurance, ',', &
							pop(i)%ame, ',', &
							pop(i)%consumption, ',', &
							pop(i)%dead, ',', &
							pop(i)%health, ',', &
							pop(i)%base_health, ',', &
							pop(i)%oop, ',', &
							pop(i)%medexp, ',', &
							pop(i)%value, ',', &
							pop(i)%tr
			end do
			close(1)

		end subroutine simulate

		subroutine simulate_fast(optrules, par)
			type (rules) optrules(nages,nmaxstates), optr(2,2)
			type (preferences) prefs
			type (params) par(npar)
			double precision pars(npar)
			integer a,i, j , n, z
			character*120 fmt
			integer fsim, esim, hsim, dead, dsim, zsim, rsim, lastdsim, lasthsim
			double precision vsim, wsim, csim, oopsim, amesim, msim,incsim, rrsim, ddsim
			double precision logearn, ben, diff, draw,  vprobe(ne),cumprob, pi
			integer w0,w1,m0,m1,s00,s01,s10,s11, solution
			double precision mu,wu, cash, tr, probh(nh)
			character*1 scn1
			character*1 scn2
			character*1 scn3
			character*1 scn4
			integer hh, ee, zz, ff, dd, mm, ww
			double precision options(4)
			integer choice, loc
			! open file to store simulated data
			pars(:) = par(:)%value
			call assign_prefs(pars,prefs)

			open(2,file='../output/simulation/simulated_pop_'//adjustl(trim(scenario))//'.csv')
			fmt = '(I8,A,I3,A,I4,A,F12.3,A,I1,A,I1,A,F12.3,A,I1,A,F12.3,A,F12.3,A,I1,A,I1,A,I1,A,F12.3,A,F12.3,A,F16.9,A,F12.3)'

			! initial population
			open(1,file='../params/input/initsample.csv')
			nobs = 0
			do i = 1, nsim, 1

				! buffer (logriearn rskill insured ass health )
				read(1,*) logearn, zsim, fsim, wsim, hsim

				! assign cohort
				if (nz .eq. 1) then
					zsim = 1
				else
					pi = 1.0d0/dble(nz)
					do z = 1, nz, 1
						if (idraws(i) .le. pi) then
							zsim = z
							exit
						else
							pi = pi + 1.0d0/dble(nz) 
						end if
					end do 
				end if


				if (wsim .lt. wmin) then
					wsim = wmin
				end if
				if (wsim .gt. wmax) then
					wsim = wmax
				end if
				! find out where fits on error term grid
				diff = dexp(logearn - dlog(earn(zsim,1)))
				call point(diff,pte, ne, esim)

				! remaining states fixed.
				dsim = 1
				amesim = 0.0d0

				! fill in optimal solution and update state-space shocks
				do a = 1, nages, 1

				    if (ages(a).eq.mcare) then
				    	fsim = 3
				    end if
				    ! automatic transitions
				    if (ages(a).ge.rmax) then
				    	esim = 1
				    end if

				    ! find out where we are in state space for m and w
				    call scale(amesim,amemin,amemax,name,gridame,m0,m1,mu)
					call scale(dsqrt(wsim),c_wmin,c_wmax,nw,c_gridw,w0,w1,wu)
					loc = getindex(a,hsim,esim,zsim,fsim,dsim,1,1)
				    optr(1,1) = optrules(a,loc + (m0-1)*nw + w0)
				    optr(1,2) = optrules(a,loc + (m0-1)*nw + w1)
				    optr(2,1) = optrules(a,loc + (m1-1)*nw + w0)
				    optr(2,2) = optrules(a,loc + (m1-1)*nw + w1)

				    ! impute values
				    ! consumption
					call blend_102(mu, wu, optr(:,:)%cons, csim)
				    ! medical exp
					call blend_102(mu, wu, optr(:,:)%medexp, msim)
					if (msim .gt. medmax) then
						msim = medmax
					end if
					if (msim .lt. 0.0d0) then
						msim = 0.0D0
					end if
					! Solution for discrete decisions
					if (ages(a).lt.rmin) then
						solution = 1
					end if
					! all retired
					if (ages(a).ge.rmax) then
						solution = 2
					end if
					! only retirement possible
					if (ages(a).ge.rmin .and. ages(a).lt.era) then
						solution = 3
					end if
					! SS possible, not medicare eligible
					if (ages(a).ge.era .and. ages(a).lt.mcare) then
						solution = 4
					end if
					! SS possible, medicare eligible
					if (ages(a).ge.mcare .and. ages(a).lt.rmax) then
						solution = 5
					end if
					! choices
					! 1 = work, no claim
					! 2 = work, claim
					! 3 = no work, no claim
					! 4 = no work, claim

					select case (solution)
					! several cases : a)
						case (1) ! Worker not eligible for pension or SS benefits

							rsim = 1
							dsim = 1
							choice = (rsim - 1)*nd + dsim
							call blend_102(mu, wu, optr(:,:)%options(choice), vsim)
							! insurance status will remain same
							fsim = fsim

							! find earnings
							incsim = earn(zsim,a)*pte(esim)
							! update his aime
							amesim = newame(amesim,incsim,a,zsim)
							if (amesim.gt.amemax) then
								amesim = amemax
							end if

							! calculate net income
							incsim = incsim + spearn(a,zsim,incsim)
							incsim = inctax(incsim,rsim,a,zsim)

							! compute cash on hand
							cash = wsim + incsim

							! check if eligible for transfers
							tr = transfer(cash,a, zsim)
							cash = cash + tr
							! out-of-pocket
							oopsim = oop(msim,a,fsim, zsim, tr)


							! update wealth
							wsim = (1.0d0+rrate)*(cash - csim - oopsim)
							if (wsim.gt.wmax) then
								wsim = wmax
							end if
							if (wsim.lt.wmin) then
								wsim = wmin
							end if

							! next period's earnings shock
							if (ages(a) .lt. rmax) then
							vprobe = cumprobe(esim,:)
								do j = 1, ne, 1
									if (edraws(i,a).lt.vprobe(j)) then
										esim = j
										exit
									end if
								end do
							end if

						case (2) ! retired folks

							rsim = 2
							dsim = 2
							choice = (rsim - 1)*nd + dsim
							call blend_102(mu, wu, optr(:,:)%options(choice), vsim)
							! do not change insurance status
							fsim = fsim

							! find retirement income (SS)
							incsim = 12.0d0*pia(amesim)
							incsim = incsim + spearn(a,zsim, incsim)

							! calculate net income
							incsim = inctax(incsim, rsim,a,zsim)

							! compute cash on hand
							cash = wsim + incsim

							! check if eligible for transfers
							tr = transfer(cash,a, zsim)
							cash = cash + tr

							! out-of-pocket
							oopsim = oop(msim,a,fsim, zsim, tr)

							! update wealth
							wsim = (1.0d0+rrate)*(cash - csim - oopsim)
							if (wsim.gt.wmax) then
								wsim = wmax
							end if
							if (wsim.lt.wmin) then
								wsim = wmin
							end if

							esim = 1

						case (3) ! Age 55, need check if decides not to work

							call blend_102(mu, wu, optr(:,:)%options(1), options(1))
							call blend_102(mu, wu, optr(:,:)%options(3), options(3))
						    ! work
						    if (options(3) .gt. options(1)) then
						    	rsim = 2
						    	vsim = options(3)
						    else
						    	rsim = 1
						    	vsim = options(1)
						    end if
						    ! claim
						    dsim = 1

							if (rsim.eq.2) then
								! need check whether looses coverage
								if (fsim.eq.2) then
									fsim = 1
								else
									fsim = fsim
								end if

								incsim = spearn(a,zsim,0.0d0)

							else
								! insurance status will remain same if already has it
								! will get tied coverage if works and did not have coverage
								if (fsim .eq. 1) then
									fsim = 2
								else 
									fsim = fsim
								end if

								! find earnings
								incsim = earn(zsim,a)*pte(esim)

								! update his aime
								amesim = newame(amesim,incsim,a,zsim)
								if (amesim.gt.amemax) then
									amesim = amemax
								end if

								! calculate gross income
								incsim = incsim + spearn(a,zsim,incsim)

								! next period's earnings shock
								vprobe = cumprobe(esim,:)
								do j = 1, ne, 1
									if (edraws(i,a).lt.vprobe(j)) then
										esim = j
										exit
									end if
								end do

							end if

							! calculate net income
							incsim = inctax(incsim, rsim,a,zsim)

							! compute cash on hand
							cash = wsim + incsim

							! check if eligible for transfers
							tr = transfer(cash, a, zsim)
							cash = cash + tr

							! out-of-pocket
							oopsim = oop(msim,a,fsim, zsim, tr)

							! update wealth
							wsim = (1.0d0+rrate)*(cash - csim - oopsim)
							if (wsim.gt.wmax) then
								wsim = wmax
							end if
							if (wsim.lt.wmin) then
								wsim = wmin
							end if


						case (4) ! Age 62, SS possible, not Medicare eligible

						    lastdsim = dsim
					    	if (lastdsim.eq.1) then
								! claiming
								call blend_102(mu, wu, optr(:,:)%options(1), options(1))
								call blend_102(mu, wu, optr(:,:)%options(2), options(2))
								call blend_102(mu, wu, optr(:,:)%options(3), options(3))
								call blend_102(mu, wu, optr(:,:)%options(4), options(4))
								choice = maxloc(options(:),dim=1)
								if (choice .eq. 1) then
									dsim = 1
									rsim = 1
									vsim = options(choice)
								else if (choice .eq. 2) then
									dsim = 2
									rsim = 1
									vsim = options(choice)
								else if (choice .eq. 3) then
									dsim = 1
									rsim = 2
									vsim = options(choice)
								else if (choice .eq. 4) then
									dsim = 2
									rsim = 2
									vsim = options(choice)
								end if
							else
								dsim = 2
								call blend_102(mu, wu, optr(:,:)%options(2), options(2))
								call blend_102(mu, wu, optr(:,:)%options(4), options(4))
								if (options(4) .gt. options(2)) then
									rsim = 2
									vsim = options(4)
								else
									rsim = 1
									vsim = options(2)
								end if
							end if

							if (dsim.eq.1) then
								if (rsim.eq.1) then

									! compute income ( earnings)
									incsim = earn(zsim,a)*pte(esim)
									! update aime
									amesim = newame(amesim,incsim,a,zsim)
									! calculate net income
									incsim = incsim + spearn(a,zsim,incsim)

									! next period's earnings shock
									vprobe = cumprobe(esim,:)
									do j = 1, ne, 1
										if (edraws(i,a).lt.vprobe(j)) then
											esim = j
											exit
										end if
									end do


								else
									! check if losses coverage
									if (fsim.eq.2) then
										fsim = 1
									else
										fsim = fsim
									end if

									! compute income (nothing otherwise)
									incsim = spearn(a,zsim,0.0d0)

								end if
							else
								if (rsim.eq.1) then
									! compute income (ssben and earnings)
									incsim = earn(zsim,a)*pte(esim)

									! benefits
									if (lastdsim.eq.1) then
										! compute income (ssben)
										ben = 12.0d0*ssben(amesim,a,zsim)
										! adjust aime so that pia(aime) = ssben
										amesim = invpia(ben/12.0d0)
									else
										ben = 12.0d0*pia(amesim)
									end if
									! Apply earnings test (if applicable)
									incsim = incsim + earntest(ben,incsim,a)
									incsim = incsim + spearn(a,zsim,incsim)

									! next period's earnings shock
									vprobe = cumprobe(esim,:)
									do j = 1, ne, 1
										if (edraws(i,a).lt.vprobe(j)) then
											esim = j
											exit
										end if
									end do

								else
									! check if looses coverage
									if (fsim.eq.2) then
										fsim = 1
									else
										fsim = fsim
									end if

									if (lastdsim.eq.1) then
										! compute income (ssben)
										incsim = 12.0d0*ssben(amesim,a,zsim)
										! adjust aime so that pia(aime) = ssben
										amesim = invpia(incsim/12.0d0)
									else
										incsim = 12.0d0*pia(amesim)
									end if

									! calculate net income
									incsim = incsim + spearn(a,zsim,incsim)

								end if

							end if


							! calculate net income
							incsim = inctax(incsim, rsim,a,zsim)

							! compute cash on hand
							cash = wsim + incsim

							! check if eligible for transfers
							tr = transfer(cash, a, zsim)
							cash = cash + tr

							! out-of-pocket
							oopsim = oop(msim,a,fsim, zsim, tr)

							! update wealth
							wsim = (1.0d0+rrate)*(cash - csim - oopsim)
							if (wsim.gt.wmax) then
								wsim = wmax
							end if
							if (wsim.lt.wmin) then
								wsim = wmin
							end if


						case (5) ! Age 65, SS claiming and Medicare eligible

					    	if (lastdsim.eq.1) then
								! claiming
								call blend_102(mu, wu, optr(:,:)%options(1), options(1))
								call blend_102(mu, wu, optr(:,:)%options(2), options(2))
								call blend_102(mu, wu, optr(:,:)%options(3), options(3))
								call blend_102(mu, wu, optr(:,:)%options(4), options(4))
								choice = maxloc(options(:),dim=1)
								if (choice .eq. 1) then
									dsim = 1
									rsim = 1
									vsim = options(choice)
								else if (choice .eq. 2) then
									dsim = 2
									rsim = 1
									vsim = options(choice)
								else if (choice .eq. 3) then
									dsim = 1
									rsim = 2
									vsim = options(choice)
								else if (choice .eq. 4) then
									dsim = 2
									rsim = 2
									vsim = options(choice)
								end if
							else
								dsim = 2
								call blend_102(mu, wu, optr(:,:)%options(2), options(2))
								call blend_102(mu, wu, optr(:,:)%options(4), options(4))
								if (options(4) .gt. options(2)) then
									rsim = 2
									vsim = options(4)
								else
									rsim = 1
									vsim = options(2)
								end if
							end if

							if (dsim.eq.1) then

								if (rsim.eq.1) then
									! earnings only
									incsim =  earn(zsim,a)*pte(esim)
									amesim = newame(amesim,incsim,a,zsim)
									! calculate net income
									incsim = incsim + spearn(a,zsim,incsim)
									! next period's earnings shock
									vprobe = cumprobe(esim,:)
									do j = 1, ne, 1
										if (edraws(i,a).lt.vprobe(j)) then
											esim = j
											exit
										end if
									end do

								else

									! compute income (nothing otherwise)
									incsim = spearn(a,zsim,0.0d0)

								end if
							else

								if (rsim.eq.1) then

									! compute income (ssben and earnings)
									incsim = earn(zsim,a)*pte(esim)

									! benefits
									if (lastdsim.eq.1) then
										! compute income (ssben)
										ben = 12.0d0*ssben(amesim,a,zsim)
										! adjust aime so that pia(aime) = ssben
										amesim = invpia(ben/12.0d0)
									else
										ben = 12.0d0*pia(amesim)
									end if
									! Apply earnings test (if applicable)
									incsim = incsim + earntest(ben,incsim,a)
									incsim = incsim + spearn(a,zsim, incsim)

									! next period's earnings shock
									vprobe = cumprobe(esim,:)
									do j = 1, ne, 1
										if (edraws(i,a).lt.vprobe(j)) then
											esim = j
											exit
										end if
									end do

								else

									! benefits
									if (lastdsim.eq.1) then
										! compute income (ssben)
										ben = 12.0d0*ssben(amesim,a,zsim)
										! adjust aime so that pia(aime) = ssben
										amesim = invpia(ben/12.0d0)
									else
										ben = 12.0d0*pia(amesim)
									end if
									incsim = ben + spearn(a,zsim,ben)

								end if

							end if

							! calculate net income
							incsim = inctax(incsim, rsim,a,zsim)

							! compute cash on hand
							cash = wsim + incsim

							! check if eligible for transfers
							tr = transfer(cash, a, zsim)
							cash = cash + tr

							! out-of-pocket
							oopsim = oop(msim,a,fsim, zsim, tr)

							! update wealth
							wsim = (1.0d0+rrate)*(cash - csim - oopsim)
							if (wsim.gt.wmax) then
								wsim = wmax
							end if
							if (wsim.lt.wmin) then
								wsim = wmin
							end if

					end select

					! Updating Health
					call getprobh(zsim,a,hsim,msim,probh,prefs)
					lasthsim = hsim
					cumprob = 0.0d0
					do j = 1, nh, 1
						cumprob = cumprob + probh(j)
						if (hdraws(i,a).lt.cumprob) then
							hsim = j
							exit
						end if
					end do

					if (hsim.eq.1) then
						dead = 2
					else
						dead = 1
					end if

					! save data
					write(2,fmt)  i, ',', &
					ages(a), ',', &
					gridz(zsim), ',', &
					wsim, ',', &
							rsim, ',', &
							dsim, ',', &
							incsim, ',', &
							fsim, ',', &
							amesim, ',', &
							csim, ',', &
							dead, ',', &
							hsim, ',', &
							lasthsim, ',', &
							oopsim, ',', &
							msim, ',', &
							vsim, ',', &
							tr
					if (dead.eq.2) then
						exit
					end if

				end do

			end do

			close(1)


		end subroutine simulate_fast

		subroutine initmpi
			integer ier
			call mpi_init(ier)
			call mpi_comm_size(MPI_COMM_WORLD,numprocs,ier)
			call mpi_comm_rank(MPI_COMM_WORLD,rank,ier)
			numworkers = numprocs - 1

			if (rank.eq.0) then
				ismaster = .true.
				write(10,'(A,I4)') 'running MPI for number of cores', numprocs
			else
				ismaster = .false.
			end if

			starttime = MPI_WTIME()
		end subroutine initmpi

		! Stopping mpi so exit graciously
		subroutine stopmpi
			! CODE TO BE TIMED
			endtime = MPI_WTIME()
			totaltime = endtime - starttime
			if (ismaster) then
				write(*,*) 'total time taken for code =',totaltime
			end if
			call mpi_finalize(ier)
		end subroutine stopmpi

		! need fix this one and the  function
		subroutine computesimmoments(pop)
			implicit none
			type (sim) pop(nobs)
			integer a, i, g, asim, zsim, j, dead, hsim, rsim, fsim, h, f, ii, nexthsim, poor, vgood
			double precision wsim, msim

			! re-initialize simulated moments in moments
			do g = 1, ngroup, 1
				moments(g)%sim(:,:) = 0.0d0
				moments(g)%nsim(:,:) = 0
			end do

				do i = 1, nobs, 1

						! find age and education

						asim = pop(i)%age
						zsim = pop(i)%z
						wsim = pop(i)%wealth
						msim = pop(i)%medexp

						!5. working
						if (pop(i)%work.eq.2) then
							rsim = 0
						else
							rsim = 1
						end if
						!8. insurance (whether insurance or not)
						if (pop(i)%insurance.gt.1) then
							fsim = 2
						else
							fsim = 1
						end if
						!11. mortality
						if (pop(i)%dead.eq.2) then
							dead = 1
						else
							dead = 0
						end if

						if (pop(i)%base_health.eq.2) then
							poor = 1
						else
							poor = 0
						end if

						if (pop(i)%base_health.eq.4) then
							vgood = 1
						else
							vgood = 0
						end if
						!12. health (indicator 3 cat)
						if (pop(i)%base_health.eq.4) then
							hsim = 3
						else if (pop(i)%base_health.eq.3) then
							hsim = 2
						else if (pop(i)%base_health.eq.2) then
							hsim = 1
						end if
						if (pop(i)%health.eq.4) then
							nexthsim = 3
						else if (pop(i)%health.eq.3) then
							nexthsim = 2
						else if (pop(i)%health.eq.2) then
							nexthsim = 1
						end if
						do g = 1, ngroup, 1
							if (moments(g)%label .eq. 'poor') then
								do j = 1, moments(g)%nages, 1
									if (asim .eq. moments(g)%ages(j)) then
										moments(g)%nsim(j,1) = moments(g)%nsim(j,1) + 1
										if (poor .eq. 1) then
											moments(g)%sim(j,1) = moments(g)%sim(j,1) + 1.0d0
										end if
									end if
								end do
							else if (moments(g)%label .eq. 'vgood') then
								do j = 1, moments(g)%nages, 1
									if (asim .eq. moments(g)%ages(j)) then
										moments(g)%nsim(j,1) = moments(g)%nsim(j,1) + 1
										if (vgood .eq. 1) then
											moments(g)%sim(j,1) = moments(g)%sim(j,1) + 1.0d0
										end if
									end if
								end do
							else if (moments(g)%label .eq. 'mortality') then
								do j = 1, moments(g)%nages, 1
									if (asim .eq. moments(g)%ages(j)) then
										moments(g)%nsim(j,1) = moments(g)%nsim(j,1) + 1
										if (dead .eq. 1) then
											moments(g)%sim(j,1) = moments(g)%sim(j,1) + 1.0d0
										end if
									end if
								end do
							else if (moments(g)%label .eq. 'wealth') then
								do j = 1, moments(g)%nages, 1
									if (asim .eq. moments(g)%ages(j)) then
										moments(g)%nsim(j,1) = moments(g)%nsim(j,1) + 1
										moments(g)%sim(j,1) = moments(g)%sim(j,1) + wsim
									end if
								end do
							else if (moments(g)%label .eq. 'medexp') then
								do j = 1, moments(g)%nages, 1
									if (asim .eq. moments(g)%ages(j)) then
										moments(g)%nsim(j,1) = moments(g)%nsim(j,1) + 1
										moments(g)%sim(j,1) = moments(g)%sim(j,1) + msim
									end if
								end do
							else if (moments(g)%label .eq. 'work') then
								do j = 1, moments(g)%nages, 1
									if (asim .eq. moments(g)%ages(j)) then
										moments(g)%nsim(j,hsim) = moments(g)%nsim(j,hsim) + 1
										if (rsim .eq. 1) then
											moments(g)%sim(j,hsim) = moments(g)%sim(j,hsim) + 1.0D0
										end if
									end if
								end do
							end if
						end do
				end do

				! compute means
				do g = 1, ngroup, 1
					do j = 1, moments(g)%nages, 1
						if (moments(g)%nh) then
							do h = 1, 3, 1
								if (moments(g)%nsim(j,h) .eq. 0) then
									moments(g)%sim(j,h) = 0.0d0
								else
									moments(g)%sim(j,h) = moments(g)%sim(j,h)/dble(moments(g)%nsim(j,h))
								end if
							end do
						else
							h = 1
							if (moments(g)%nsim(j,h) .eq. 0) then
								moments(g)%sim(j,h) = 0.0d0
							else
								moments(g)%sim(j,h) = moments(g)%sim(j,h)/dble(moments(g)%nsim(j,h))
							end if
						end if
					end do
				end do

		end subroutine computesimmoments

		! need fix this one and the  function
		subroutine computetarget(pop, sim_e25, sim_medexp)
			implicit none
			type (sim) pop(nobs)
			integer i, asim, a
			double precision life_sim(nages,2)
			double precision sim_e25, sim_medexp

			! target for medical spending
			sim_medexp = sum(pop(:)%medexp)/dble(nobs)
			life_sim(:,:) = 0.0d0
			! get life table
			do i = 1, nobs, 1
				asim = pop(i)%age - agemin + 1
				life_sim(asim,1) = life_sim(asim,1) + 1.0d0
			end do
			! compute survival rates
			do a = 1, nages, 1
				life_sim(a,2) = life_sim(a,1)/life_sim(1,1)
			end do
			sim_e25 = 25.0d0 + sum(life_sim(:,2))

		end subroutine computetarget



		subroutine comparemoments(distance)
			double precision distance(nfreemoments)
			double precision simmom(nfreemoments)
			integer j,h, g, m, i

			! collect simulated moments in a vector
			m = 1
			do g = 1, ngroup, 1
				if (moments(g)%free) then
       				if (moments(g)%nh) then
						do h = 1, nh-1, 1
							do j = 1, moments(g)%nages, 1
								simmom(m) = moments(g)%sim(j,h)
								m = m + 1
							end do
						end do
					else
						h = 1
						do j = 1, moments(g)%nages, 1
							simmom(m) = moments(g)%sim(j,h)
							m = m + 1
						end do	
					end if
				end if
			end do

			! compute distance vector
                        distance(:) = 0.0d0
			do m = 1, nfreemoments
				do i = 1, ndata, 1
					if (select(i,m) .eq. 1) then
						distance(m) = distance(m) + (data(i,m) - simmom(m))
					end if
				end do
				distance(m) = distance(m) / dble(nind)
			end do
		end subroutine comparemoments

		subroutine dmoments(freepar, gm)
			double precision freepar(nfreepar), gm(nfreemoments), simmom(nfreemoments)
			type (params) par(npar)
			type (rules) optrules(nages, nmaxstates)
			type (sim) pop(nsim*nages)
			integer i, n, g, h, j
			! save current structural parameters to file
			par(:) = g_initpar(:)
			call setfreepar(par,freepar)

			!  solve for decision rules
			call recursion(par, optrules)
			if (isave) then
				call saverules(optrules)
			end if
			! simulate data
			call simulate(optrules, pop, par)

			! compute moments from simulated dat
			call computesimmoments(pop(1:nobs))

			! gather moments from simulations in a vector
			simmom(:) = 0.0d0
			n = 1
			do g = 1, ngroup, 1
				if (moments(g)%free) then
					do j = 1, moments(g)%nages, 1
						if (moments(g)%nh) then
							do h = 1, nh-1, 1
								simmom(n) = moments(g)%sim(j,h)
								n = n + 1
							end do
						else
							h = 1
							simmom(n) = moments(g)%sim(j,h)
							n = n + 1
						end if
					end do
				end if
			end do

			! compute distance vector
			gm(:) = 0.0D0 
			do n = 1, nfreemoments
				do i = 1, ndata, 1
					if (select(i,n) .eq. 1) then
						gm(n) = gm(n) + (data(i,n) - simmom(n))
					end if
				end do
				gm(n) = gm(n) / dble(nind)
			end do
	
		end subroutine dmoments

 		subroutine docovariance(par)
 			type (params) par(npar)
 			double precision freepar(nfreepar), eps(nfreepar), freepar_up(nfreepar), freepar_down(nfreepar)
 			double precision dmom(nfreemoments,nfreepar), weight(nfreemoments,nfreemoments), covar(nfreepar,nfreepar)
 			double precision gm_up(nfreemoments), gm_down(nfreemoments), simmom(nfreemoments)
 			double precision invcovar(nfreepar,nfreepar)
 			integer i, ii, j, jj, errorflag, n, h, g, m

			! extract freepar
			call extractfreepar(par, freepar)

 			! getting new optimal weighting matrix at simulated values
 			simmom(:) = 0.0d0
			n = 1
			do g = 1, ngroup, 1
				if (moments(g)%free) then
					do j = 1, moments(g)%nages, 1
						if (moments(g)%nh) then
							do h = 1, nh-1, 1
								simmom(n) = moments(g)%sim(j,h)
								n = n + 1
							end do
						else
							h = 1
							simmom(n) = moments(g)%sim(j,h)
							n = n + 1
						end if
					end do
				end if
			end do
			!Recompute covariance matrix of data
			do i = 1, nfreemoments
				do j = 1, nfreemoments, 1
					do n = 1, ndata, 1
						! whether observation contributes to both moments
						if (select(n,i) .eq. 1 .and. select(n,j) .eq. 1 ) then
							datacov(i,j) = datacov(i,j) + (data(n,i) - simmom(i))*(data(n,j) - simmom(j))/dble(nind)
						end if
					end do
				end do
			end do
			! call inverse of covariance matrix of data to get optimal weight matrix
			call invert(datacov, optw,nfreemoments,errorflag)
			open(1,file='../output/estimation/weight_'//adjustl(trim(scenario))//'.csv')
			do i = 1, nfreemoments, 1
				write(1,*) optw(i,:)
			end do 
			close(1)
			! steps for derivative (using same step as French, 2005)
 			! dh = (.1)*maxc((ax0~(1e-17)*ones(rows(x0),1))').*dax0;/* the .05 part is my modification -- ideally,
			! it should be 1e-8 instead of .05*/ dax0 = x0/dabs(x0)
			write(10,*) ' step size for standard errors :'
 			do i = 1, nfreepar
				 eps(i) = 0.1d0 * dabs(freepar(i))
				 if (eps(i) .gt. 5.0) then
					eps(i) = 5.0d0
				 end if
				 if (eps(i) .lt. 1.0d-3) then
					eps(i) = 1.0d-3
				 end if
 				write(10,*) i, eps(i)
			 end do
			if (iload) then
				write(10,*) ' gradient of moment vector (D matrix): '
				open(1,file='../output/estimation/gradient_'//adjustl(trim(scenario))//'.csv')
				do i = 1,nfreemoments,1
					read(1,*) dmom(i,:)
					write(*,*) dmom(i,:)
				end do
				close(1)
			else
				! for each parameter, compute derivative of moment vector
				dmom(:,:) = 0.0d0
				do i = 1, nfreepar
						! start from estimated vector
						freepar_up = freepar
						freepar_down = freepar
						! add eps to par(i)
						freepar_up(i) = freepar_up(i) + eps(i)
						freepar_down(i) = freepar_down(i) - eps(i)
						gm_up(:) = 0.0d0
						gm_down(:) = 0.0d0
						! compute moments
						call dmoments(freepar_up, gm_up)
						call dmoments(freepar_down, gm_down)
						! compute derivative
						do j = 1, nfreemoments, 1
							dmom(j,i) = (gm_up(j)-gm_down(j))/(2.0d0*eps(i))
						end do
				end do
				write(10,*) ' gradient of moment vector (D matrix): '
				open(1,file='../output/estimation/gradient_'//adjustl(trim(scenario))//'.csv')
				do j = 1, nfreemoments, 1
					write(10,*) dmom(j,:)
					write(1,*) dmom(j,:)
				end do
				close(1)
			end if
 			! compute covariance matrix
			do i = 1,nfreepar, 1
				do j = 1, nfreepar,1
					covar(i,j) = 0.0d0
					do m = 1, nfreemoments, 1
						covar(i,j) = covar(i,j) + dmom(m,j)*dmom(m,i)*optw(m,m)
					end do
				end do
			end do
			write(10,*) ' DWD matrix : '
			do j = 1, nfreepar, 1
				write(10,*) covar(j,:)
			end do
			call invert(covar, invcovar,nfreepar,errorflag)
 			write(10,*) 'check for error flag on  inversion : ', errorflag
 			if (errorflag .eq. -1) then
 				do i = 1, nfreepar, 1
 					write(10,*) covar(i,:)
 				end do
 			end if
             do i = 1, nfreepar, 1
 				do j = 1, nfreepar, 1
 					covar(i,j) = (1.0d0+dble(nind)/dble(nsim))*invcovar(i,j)/dble(nind)
 				end do
 			end do

 			! get standard errors
 			ii = 1
 			do i = 1, npar, 1
 				if (par(i)%free) then
 					if (covar(ii,ii) .gt. 0.0d0) then
 						par(i)%serr = dsqrt(covar(ii,ii))
 						par(i)%pvalue = 1.0d0 - probn(dabs(par(i)%value/par(i)%serr))
 					else
						par(i)%serr = 0.0d0
 						par(i)%pvalue = 0.0d0
 					end if
 					ii = ii + 1
 				else
 					par(i)%serr = 0.0D0
 					par(i)%pvalue = 0.0D0
 				end if
 			end do

			! saving covariance matrix
			open(1,file='../output/estimation/covar_'//adjustl(trim(scenario))//'.csv')
			do i = 1, npar, 1
				write(1,*) covar(i,:)
			end do
			close(1)


 		end subroutine docovariance

		 ! Returns the inverse of a matrix calculated by finding the LU
		! decomposition.  Depends on LAPACK.
		!  subroutine invert(A,Ainv) 
		! 	double precision, dimension(:,:), intent(in) :: A
		! 	double precision, dimension(size(A,1),size(A,2)) :: Ainv
		
		! 	double precision, dimension(size(A,1)) :: work  ! work array for LAPACK
		! 	integer, dimension(size(A,1)) :: ipiv   ! pivot indices
		! 	integer :: n, info
		
		! 	! External procedures defined in LAPACK
		! 	external DGETRF
		! 	external DGETRI
		
		! 	! Store A in Ainv to prevent it from being overwritten by LAPACK
		! 	Ainv = A
		! 	n = size(A,1)
		
		! 	! DGETRF computes an LU factorization of a general M-by-N matrix A
		! 	! using partial pivoting with row interchanges.
		! 	call DGETRF(n, n, Ainv, n, ipiv, info)
		
		! 	if (info /= 0) then
		! 	stop 'Matrix is numerically singular!'
		! 	end if
		
		! 	! DGETRI computes the inverse of a matrix using the LU factorization
		! 	! computed by DGETRF.
		! 	call DGETRI(n, Ainv, n, ipiv, work, n, info)
		
		! 	if (info /= 0) then
		! 	stop 'Matrix inversion failed!'
		! 	end if
		! end subroutine invert

! 			!Subroutine to find the inverse of a square matrix
! 			!Author : Louisda16th a.k.a Ashwith J. Rego
! 			!Reference : Algorithm has been well explained in:
! 			!http://math.uww.edu/~mcfarlat/inverse.htm
! 			!http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
			subroutine invert(matrix, inverse, n, errorflag)
				implicit none
				!declarations
				integer, intent(in) :: n
				integer, intent(out) :: errorflag  !return error status. -1 for error, 0 for normal
				double precision, intent(in), dimension(n,n) :: matrix  !input matrix
				double precision, intent(out), dimension(n,n) :: inverse !inverted matrix

				logical :: flag = .true.
				integer :: i, j, k, ll
				double precision :: m
				double precision, dimension(n,2*n) :: augmatrix !augmented matrix

				!augment input matrix with an identity matrix
				do i = 1, n
						do j = 1, 2*n
								if (j <= n ) then
										augmatrix(i,j) = matrix(i,j)
								else if ((i+n) == j) then
										augmatrix(i,j) = 1.0d0
								else
										augmatrix(i,j) = 0.0d0
								endif
						end do
				end do

				!reduce augmented matrix to upper traingular form
				do k =1, n-1
						if (augmatrix(k,k) == 0.0d0) then
								flag = .false.
								do i = k+1, n
										if (augmatrix(i,k) /= 0.0d0) then
												do j = 1,2*n
														augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
												end do
												flag = .true.
												exit
										endif
										if (flag .eqv. .false.) then
												print*, "matrix is non - invertible"
												inverse = 0.0d0
												errorflag = -1
												return
										endif
								end do
						endif
						do j = k+1, n
								m = augmatrix(j,k)/augmatrix(k,k)
								do i = k, 2*n
										augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
								end do
						end do
				end do

				!test for invertibility
				do i = 1, n
						if (augmatrix(i,i) == 0.0d0) then
								print*, "matrix is non - invertible"
								inverse = 0.0d0
								errorflag = -1
								return
						endif
				end do

				!make diagonal elements as 1
				do i = 1 , n
						m = augmatrix(i,i)
						do j = i , (2 * n)
								   augmatrix(i,j) = (augmatrix(i,j) / m)
						end do
				end do

				!reduced right side half of augmented matrix to identity matrix
				do k = n-1, 1, -1
						do i =1, k
						m = augmatrix(i,k+1)
								do j = k, (2*n)
										augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
								end do
						end do
				end do

				!store answer
				do i =1, n
						do j = 1, n
								inverse(i,j) = augmatrix(i,j+n)
						end do
				end do
				errorflag = 0
			end subroutine invert

		double precision function probn(value)
			double precision value,  q, bound
			integer ifail, status
			call cdfnor(1,probn,q, value, 0.0d0, 1.0d0, status, bound)
		end function

		double precision function quann(prob)
			double precision prob, bound
			integer ifail, status
			ifail = 1
			call cdfnor(2,prob,1.0d0-prob, quann, 0.0d0, 1.0d0, status, bound)
		end function

        integer function year(a,z)
			integer a, z
            year = gridz(z) + agemin + a - 1
        end function

		! Utility function: French (2005) with intercept
		double precision function utility(c,r,h,a,z, m, prefs)
			double precision c,leisure,g, m, x
			integer r,h, a, z
			type (preferences) prefs
			leisure = prefs%leisure_endow
			! loss of leisure due to work
			if (r.eq.1) then
				leisure = leisure - hours_worked
			end if
			! disutility of bad health
			if (h.eq.3) then
				leisure = leisure - prefs%phi(1)
			end if
			if (h.eq.2) then
				leisure = leisure - prefs%phi(2)
			end if
			! aggregate composite good
			x = (c **prefs%psi) * (leisure**(1.0d0-prefs%psi))
			utility =  (x**(1.0d0-prefs%sigma) - prefs%alpha(h))/(1.0d0-prefs%sigma)
			utility = utility*1d9
		end function


		! ! A la French
		double precision function bequest(ass, a, h,  prefs)
			double precision ass, cb, phib, kapa0
			integer a, h
			double precision leisure, compose
			type (preferences) prefs
			!kapa0 = prefs%kapa(1)*dexp(prefs%kapa(3)*dble(a-40))
			cb = prefs%kapa(2)
			!compose = ((cb + ass)**prefs%psi) * (leisure**(1.0d0 - prefs%psi))
			compose = ((cb + ass)**prefs%psi)
			bequest =  prefs%kapa(1)  * (compose**(1.0d0 - prefs%sigma) - prefs%kapa(3)) &
				/ (1.0d0 - prefs%sigma)
			bequest = bequest*1d9
		end function


		!* transfer function (standard Hubbard, Skinner and Zeldes)
		double precision function transfer(x, a, z)
			double precision x
			integer a, z
				if (x.ge.cashmin(z,a)) then
					transfer = 0.0d0
				else
					transfer = cashmin(z,a) - x
				end if
		end function

		subroutine getprobh(z,a,h,med,probh,prefs)
			integer a, h, j, z,yrs
			double precision med, adjmed, probh(nh), exb(nh-1), inv,cexb, theta_a(2), imp_m, imp_h, imp_g(3)
			type (preferences) prefs
			! probability of death
			probh(1) = 1.0d0 - probs(h-1,a)
			if (probh(1).gt.1.0d0) then
				probh(1) = 1.0d0
			end if
			adjmed = med
			! test
			!adjmed = medmax
			if (adjmed .gt. medmax) then
				adjmed = medmax
			end if
            if (adjmed .lt. 0.0d0) then
                adjmed = 0.0d0
            end if
			imp_m = dexp(-prefs%prog_med * dble(progress(z,a)))

			! probability of health state conditional on survival
			cexb = 0.0d0
			do j = 1, nh-1, 1
				theta_a(1) = prefs%theta(j)*imp_m
				theta_a(2) = prefs%theta2(j)*imp_m
				exb(j) = lambda(z,a,h-1,j)*dexp(theta_a(1)*dlog(1.0d0+ adjmed) &
						+ theta_a(2)*(dlog(1.0d0 + adjmed)**2))
				cexb = cexb + exb(j)
			end do
			do j = 2, nh, 1
				probh(j) = (1.0d0-probh(1))*exb(j-1)/cexb
			end do

		end subroutine getprobh

		double precision function spearn(a, z, income)
			double precision income
			integer a, z
			spearn = otherinc(z,a) + parother(1)*income
			if (spearn.lt.0.0d0) then
				spearn = 0.0d0
			end if

		end function

		! out-of-pocket formula
		double precision function oop(m, a, f, z, tr)
			double precision m, tr
			integer a, f, z
			if (f.eq.1) then
				! check whether receiving transfers
				if (tr.gt.0.0d0) then
					! if so medicaid pays
					oop = mcdcopay*insgen(z,a)*m
				else
					oop = copay(z,a,f)*m
				end if
			else
				oop = copay(z,a,f)*m
			end if

		end function

		double precision function invoop(oop, a, f, z, tr)
			double precision oop, tr
			integer a, f, z
			if (f.eq.1) then
				! check whether categorically needy
				if (tr.gt.0.0d0) then
					! if so medicaid pays
					invoop = oop/(mcdcopay*insgen(z,a))
				else
					invoop = oop/copay(z,a,f)
				end if
			else
				invoop = oop/copay(z,a,f)
			end if
		end function

		! Formula to update AIME (comes from French and Jones, 2011)
		double precision function newame(oldame,earn,a, z)
			double precision oldame,earn,monthly,change,g
			integer a,i, z
			! cap earnings and transform to monthly
			if (earn.gt.ssmax(z,a)) then
				monthly = ssmax(z,a)/12.0d0
			else
				monthly = earn/12.0d0
			end if
			! update aime using French and Jones Formula
			! whether real earnings growth (prior to age 60 SS)
			if (ages(a).le.60) then
				g = 1.0d0
			else
				g = 1.0d0
			end if
			! compute change
			change = (monthly - rep(a)*g*oldame)/35.0d0
			! check if change positive, if so update, if not keep old
			if (change.gt.0.0d0) then
				newame = g*oldame + change
			else
				newame = g*oldame
			end if
		end function

		! PIA formula (depends on calendar year when reach age 62)
		double precision function pia(ame)
			double precision ame
			if (ame.le.piabend(1)) then
				pia = piarate(1)*ame
			else if (ame.gt.piabend(1) .and. ame.le.piabend(2)) then
				pia = piarate(1)*piabend(1) + piarate(2)*(ame - piabend(1))
			else
				pia = piarate(1)*piabend(1) + piarate(2)*(piabend(2)-piabend(1)) + piarate(3)*(ame - piabend(2))
			end if

		end function

		! Social Security benefit formula, used when ind. just claimed to bet benefit
		double precision function ssben(ame,a,z)
			integer a, z
			double precision ame,newpia
			newpia = ssgen(z)*pia(ame)
			if (ages(a).lt.nra) then
				ssben = newpia*((1.0d0-arf)**(nra-ages(a)))
			else
				ssben = newpia*(1.0d0+drc)**(ages(a) - nra)
			end if
		end function

		double precision function earntest(ben,earn,a)
			double precision ben, earn, tax
			integer a
			! earnings test
			tax = 0.0d0
			if (earn.gt.earntaxgain) then
				tax = earn - earntaxgain
				if (ages(a).lt.nra) then
					tax = tax*earntaxrate(1)
				else if (ages(a).ge.nra) then
					tax = tax*earntaxrate(2)
				end if
			end if
			if (tax.lt.ben) then
				earntest = ben - tax
			else
				earntest = 0.0d0
			end if
		end function

		! Invert PIA back to AIME (used by equivame, pia depends on year turns 62)
		double precision function invpia(oldpia)
			double precision oldpia,bend(2),rate(3)
			rate = piarate(:)
			bend(1) = pia(piabend(1))
			bend(2) = pia(piabend(2))
			if (oldpia.le.bend(1)) then
				invpia = oldpia/rate(1)
			else if (oldpia.gt.bend(1) .and. oldpia.le.bend(2)) then
				invpia = bend(1)/rate(1) + (oldpia - bend(1))/rate(2)
			else
				invpia = bend(1)/rate(1) + (bend(2) - bend(1))/rate(2) + (oldpia - bend(2))/rate(3)
			end if
		end function


		! Tax function (from Gouveia and Strauss, social security taxes and medicare taxes)
		double precision function inctax(inc, ret,a,z)
			double precision inc, inck, atr,sstaxpaid,mctaxpaid
			integer ret, a, z
			inck = 0.001d0*inc
			! Federal income tax
			if (inck.gt.0.0d0) then
				atr = taxpar(z,a,1) - taxpar(z,a,1)*((taxpar(z,a,3)*(inck**taxpar(z,a,2)) + 1.0d0)**(-1.0d0/taxpar(z,a,2)))
				inctax = inc*(1.0d0-atr)
				! Social security and medicare tax
				sstaxpaid = 0.0
				mctaxpaid = 0.0
				if (ret.eq.1 ) then
					if (inc .lt. ssmax(z,a)) then
						sstaxpaid = sstax(z,a)*inc 
					else 
						sstaxpaid = sstax(z,a)*ssmax(z,a)
					end if	
					if (inc .lt. mcmax(z,a)) then
						mctaxpaid = mctax(z,a)*inc 
					else 
						mctaxpaid = mctax(z,a)*mcmax(z,a)
					end if	
				end if
				inctax = inctax - sstaxpaid - mctaxpaid
			else
				inctax = 0.0d0
			end if

			if (inctax.lt.0.0d0) then
				inctax = 0.0d0
			end if

		end function

	subroutine nelmin_params ( fn, n, start, xmin, ynewlo, reqmin, step, konvge, kcount, &
	  icount, numres, ifault )
	!
	!    Output, integer ( kind = 4 ) IFAULT, error indicator.
	!    0, no errors detected.
	!    1, REQMIN, N, or KONVGE has an illegal value.
	!    2, iteration terminated because KCOUNT was exceeded without convergence.
	!
	  implicit none

	  integer ( kind = 4 ) n

	  real ( kind = 8 ), parameter :: ccoeff = 0.5D+00
	  real ( kind = 8 ) del
	  real ( kind = 8 ), parameter :: ecoeff = 2.0D+00
	  real ( kind = 8 ), parameter :: eps = 0.001D+00
	  real ( kind = 8 ), external :: fn
	  integer ( kind = 4 ) i
	  integer ( kind = 4 ) icount
	  integer ( kind = 4 ) ifault
	  integer ( kind = 4 ) ihi
	  integer ( kind = 4 ) ilo
	  integer ( kind = 4 ) j
	  integer ( kind = 4 ) jcount
	  integer ( kind = 4 ) kcount
	  integer ( kind = 4 ) konvge
	  integer ( kind = 4 ) l
	  integer ( kind = 4 ) numres
	  real ( kind = 8 ) p(n,n+1)
	  real ( kind = 8 ) p2star(n)
	  real ( kind = 8 ) pbar(n)
	  real ( kind = 8 ) pstar(n)
	  real ( kind = 8 ), parameter :: rcoeff = 1.0D+00
	  real ( kind = 8 ) reqmin
	  real ( kind = 8 ) rq
	  real ( kind = 8 ) start(n)
	  real ( kind = 8 ) step(n)
	  real ( kind = 8 ) x
	  real ( kind = 8 ) xmin(n)
	  real ( kind = 8 ) y(n+1)
	  real ( kind = 8 ) y2star
	  real ( kind = 8 ) ylo
	  real ( kind = 8 ) ynewlo
	  real ( kind = 8 ) ystar
	  real ( kind = 8 ) z
	!
	!  Check the input parameters.
	!
	  if ( reqmin <= 0.0D+00 ) then
	    ifault = 1
	    return
	  end if

	  if ( n < 1 ) then
	    ifault = 1
	    return
	  end if

	  if ( konvge < 1 ) then
	    ifault = 1
	    return
	  end if
	!
	!  Initialization.
	!
	  icount = 0
	  numres = 0
	  jcount = konvge
	  del = 1.0D+00
	  rq = reqmin * real ( n, kind = 8 )
	!
	!  Initial or restarted loop.
	!
	  do

	    p(1:n,n+1) = start(1:n)
	    y(n+1) = fn ( start )
	    icount = icount + 1
	!
	!  Define the initial simplex.
	!
	    do j = 1, n
	      x = start(j)
	      start(j) = start(j) + step(j) * del
	      p(1:n,j) = start(1:n)
	      y(j) = fn ( start )
	      icount = icount + 1
	      start(j) = x
	    end do
	!
	!  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
	!  the vertex of the simplex to be replaced.
	!
	    ilo = minloc ( y(1:n+1), 1 )
	    ylo = y(ilo)
	!
	!  Inner loop.
	!
	    do while ( icount < kcount )
	!
	!  YNEWLO is, of course, the HIGHEST value???
	!
	      ihi = maxloc ( y(1:n+1), 1 )
	      ynewlo = y(ihi)
	!
	!  Calculate PBAR, the centroid of the simplex vertices
	!  excepting the vertex with Y value YNEWLO.
	!
	      do i = 1, n
	        pbar(i) = ( sum ( p(i,1:n+1) ) - p(i,ihi) ) / real ( n, kind = 8 )
	      end do
	!
	!  Reflection through the centroid.
	!
	      pstar(1:n) = pbar(1:n) + rcoeff * ( pbar(1:n) - p(1:n,ihi) )
	      ystar = fn ( pstar )
	      icount = icount + 1
	!
	!  Successful reflection, so extension.
	!
	      if ( ystar < ylo ) then

	        p2star(1:n) = pbar(1:n) + ecoeff * ( pstar(1:n) - pbar(1:n) )
	        y2star = fn ( p2star )
	        icount = icount + 1
	!
	!  Retain extension or contraction.
	!
	        if ( ystar < y2star ) then
	          p(1:n,ihi) = pstar(1:n)
	          y(ihi) = ystar
	        else
	          p(1:n,ihi) = p2star(1:n)
	          y(ihi) = y2star
	        end if
	!
	!  No extension.
	!
	      else

	        l = 0
	        do i = 1, n + 1
	          if ( ystar < y(i) ) then
	            l = l + 1
	          end if
	        end do

	        if ( 1 < l ) then

	          p(1:n,ihi) = pstar(1:n)
	          y(ihi) = ystar
	!
	!  Contraction on the Y(IHI) side of the centroid.
	!
	        else if ( l == 0 ) then

	          p2star(1:n) = pbar(1:n) + ccoeff * ( p(1:n,ihi) - pbar(1:n) )
	          y2star = fn ( p2star )
	          icount = icount + 1
	!
	!  Contract the whole simplex.
	!
	          if ( y(ihi) < y2star ) then

	            do j = 1, n + 1
	              p(1:n,j) = ( p(1:n,j) + p(1:n,ilo) ) * 0.5D+00
	              xmin(1:n) = p(1:n,j)
	              y(j) = fn ( xmin )
	              icount = icount + 1
	            end do

	            ilo = minloc ( y(1:n+1), 1 )
	            ylo = y(ilo)

	            cycle
	!
	!  Retain contraction.
	!
	          else
	            p(1:n,ihi) = p2star(1:n)
	            y(ihi) = y2star
	          end if
	!
	!  Contraction on the reflection side of the centroid.
	!
	        else if ( l == 1 ) then

	          p2star(1:n) = pbar(1:n) + ccoeff * ( pstar(1:n) - pbar(1:n) )
	          y2star = fn ( p2star )
	          icount = icount + 1
	!
	!  Retain reflection?
	!
	          if ( y2star <= ystar ) then
	            p(1:n,ihi) = p2star(1:n)
	            y(ihi) = y2star
	          else
	            p(1:n,ihi) = pstar(1:n)
	            y(ihi) = ystar
	          end if

	        end if

	      end if
	!
	!  Check if YLO improved.
	!
	      if ( y(ihi) < ylo ) then
	        ylo = y(ihi)
	        ilo = ihi
	      end if

	      jcount = jcount - 1

	      if ( 0 < jcount ) then
	        cycle
	      end if
	!
	!  Check to see if minimum reached.
	!
	      if ( icount <= kcount ) then

	        jcount = konvge

	        x = sum ( y(1:n+1) ) / real ( n + 1, kind = 8 )
	        z = sum ( ( y(1:n+1) - x )**2 )

	        if ( z <= rq ) then
	          exit
	        end if

	      end if

	    end do
	!
	!  Factorial tests to check that YNEWLO is a local minimum.
	!
	    xmin(1:n) = p(1:n,ilo)
	    ynewlo = y(ilo)

	    if ( kcount < icount ) then
	      ifault = 2
	      exit
	    end if

	    ifault = 0

	    do i = 1, n
	      del = step(i) * eps
	      xmin(i) = xmin(i) + del
	      z = fn ( xmin )
	      icount = icount + 1
	      if ( z < ynewlo ) then
	        ifault = 2
	        exit
	      end if
	      xmin(i) = xmin(i) - del - del
	      z = fn ( xmin )
	      icount = icount + 1
	      if ( z < ynewlo ) then
	        ifault = 2
	        exit
	      end if
	      xmin(i) = xmin(i) + del
	    end do

	    if ( ifault == 0 ) then
	      exit
	    end if
	!
	!  Restart the procedure.
	!
	    start(1:n) = xmin(1:n)
	    del = eps
	    numres = numres + 1

	  end do

	  return
	end	subroutine nelmin_params

   subroutine scale(q,qmin,qmax,nq,gridq,q0,q1,qu)
   		integer nq,q0,q1
   		double precision q, qmin, qmax, gridq(nq),gapq, qu
   		gapq = (qmax - qmin)/dble(nq-1)
   		q0 = floor((q - qmin)/gapq) + 1
		if (q0.lt.1) then
			q = qmin
			q0 = 1
			q1 = 2
		else if (q0.ge.nq) then
			q1 = nq
			q0 = nq - 1
		else
			q1 = q0 + 1
		end if
		qu = (q - gridq(q0))/(gridq(q1)-gridq(q0))

   end subroutine scale


	! routine to linearly intrapolate in two dimensions
	subroutine blend_102 ( r, s, xm, x )
	  double precision r, s, x, xm(2,2)
	  x = xm(1,1) + r*(-xm(1,1)+xm(2,1)) + s*(-xm(1,1)+xm(1,2)) + r*s*(xm(1,1) - xm(2,1) - xm(1,2) + xm(2,2))
	end subroutine

	integer function getindex(a, h,e, z, f, d, ame, w)
		integer h, e, z, f, d, ame, w, a, i
		integer s(nstatevars)
		s = (/ame,w,h,e,f,d,z/)
		do i = 1, nstates(a),1
			if (all(states(a,i,:) .eq. s)) then
				getindex = i
				exit
			end if
		end do
	end function getindex

	elemental double precision function tf(level)
		double precision, intent(in):: level
		tf = level ** curv
	end function tf

	elemental double precision function itf(trans)
		double precision, intent(in):: trans
		itf = trans ** (1.0d0/curv)
	end function itf

	! routine to find the closest point index on a grid
	subroutine closest(val, minval, step, n, low, up)
		double precision val, step,minval
		integer n, low, up
		low = int(floor((val - minval)/step))+1
		if (low.lt.1) then
			low = 1
			up = 2
		else if (low.eq.n) then
			up = low
			low = low - 1
		else
			up = low + 1
		end if
	end subroutine

	! routine to find the closest point on a grid
	subroutine point(val, grid, n, p)
		integer n,i,imin,p
		double precision val, grid(n),diff,dmin
		dmin = 1.0d20; imin = 1
		do i = 1,n,1
			diff = dabs(val - grid(i))
			if (diff.lt.dmin) then
				dmin = diff
				imin = i
			end if
		end do
		p = imin
	end	subroutine

end module dphealth
