program main
	use dphealth
	! reading arguments from command line and assign to scenario name
    narg = command_argument_count()
    call get_command_argument(1,scenario)

 	call produce_draws

end program main