program gen
	use dphealth
	! reading arguments from command line and assign to scenario name
    narg = command_argument_count()
    call get_command_argument(1,scenario)

 	! run simulation for this scenario
 	call generate

end program gen