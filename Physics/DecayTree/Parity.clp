(defrule check-parity
	(declare (salience 99))
	?decay <- (Decay 
					(quantum_number_name "parity") (mother ?parity_mother) 
					(daughters ?parity_daughter1 ?parity_daughter2 $?others)
					(required_variable_names "angular-momentum" $?other_rvns) 
					(violating_quantum_number_list $?violating_quantum_number_list)
			  )
	(test (not (member$ "parity" ?violating_quantum_number_list)))
	=>
	;get the required information
	(bind ?angular_momentum (get-spin-qn-with-unique-id (get-required-variable "angular-momentum" ?decay)))
	(bind ?L (/ (fact-slot-value ?angular_momentum numerator) (fact-slot-value ?angular_momentum denominator)))
	;(printout t ?L " " ?parity_mother " " ?parity_daughter1 " "  ?parity_daughter2 crlf)
	(if (<> ?parity_mother (* (* ?parity_daughter1 ?parity_daughter2) (** -1 ?L)))
	then
		(if (is-qn-conserved "parity")
		then
			(retract ?decay)
	  	else
	  		(modify ?decay (violating_quantum_number_list ?violating_quantum_number_list "parity"))
	  	)
	  	;(printout t "decay violates parity!" crlf)
	)
)

(defrule check-cparity
	(declare (salience 99))
	?charge_decay <- (Decay 
						(quantum_number_name "charge") (mother ?charge_mother)
						(daughters ?charge_daughter1 ?charge_daughter2 $?others)
					 )
	?decay <- (Decay 
				(quantum_number_name "cparity") (mother ?cparity_mother)
				(daughters ?cparity_daughter1 ?cparity_daughter2 $?others)
				(required_decays ?charge_decay $?required_decays)
				(violating_quantum_number_list $?violating_quantum_number_list)
			  )
	(test (not (member$ "cparity" ?violating_quantum_number_list)))
	=>
	(if (or (or (<> ?charge_mother 0) (<> ?charge_daughter1 0)) (<> ?charge_daughter2 0))
	then
		(retract ?decay)
	else
		(if (<> ?cparity_mother (* ?cparity_daughter1 ?cparity_daughter2))
		then		 
	 	 	(if (is-qn-conserved "cparity")
	  		then
				(retract ?decay)
				;(printout t "decay violates cparity!" crlf)
	  		else
	  			(modify ?decay 
	  				(violating_quantum_number_list ?violating_quantum_number_list "cparity")
	  			)
	  		)
			;(printout t "decay violates cparity!" crlf)
		)
	)
)