;;;*************************
;;;* DECAY VIOLATION RULES *
;;;*************************

(defrule check-isospin-z
	(declare (salience 99))
	(SpinQuantumNumber (unique_id ?mother_id) (z_component_numerator ?z_num_mother) (denominator ?denom_mother))
	(SpinQuantumNumber (unique_id ?daughter1_id) (z_component_numerator ?z_num_daughter1) (denominator ?denom_daughter1))
	(SpinQuantumNumber (unique_id ?daughter2_id) (z_component_numerator ?z_num_daughter2) (denominator ?denom_daughter2))
	?decay <- (Decay 
				(quantum_number_name "isospin") (mother ?mother_id)
	 			(daughters ?daughter1_id ?daughter2_id $?others)
	 			(violating_quantum_number_list $?violating_quantum_number_list)
	 		  )
	=>
	(if (<> (/ ?z_num_mother ?denom_mother) 
			(+ (/ ?z_num_daughter1 ?denom_daughter1) (/ ?z_num_daughter2 ?denom_daughter2))
		)
	then
		(if (is-qn-conserved "isospin")
			(retract ?decay)
	  	else
	  		(modify ?decay (violating_quantum_number_list ?violating_quantum_number_list "isospin-z"))
	  	)
		;(printout t "decay violates isospin z!" crlf)
	)
)

(defrule check-isospin
	(declare (salience 99))
	(SpinQuantumNumber (unique_id ?mother_id) (numerator ?num_mother) (denominator ?denom_mother))
	(SpinQuantumNumber (unique_id ?daughter1_id) (numerator ?num_daughter1) (denominator ?denom_daughter1))
	(SpinQuantumNumber (unique_id ?daughter2_id) (numerator ?num_daughter2) (denominator ?denom_daughter2))
	?decay <- (Decay 
				(quantum_number_name "isospin") (mother ?mother_id)
	 			(daughters ?daughter1_id ?daughter2_id $?others)
	 			(violating_quantum_number_list $?violating_quantum_number_list)
	 		  )
	=>
	(if (or
			(> (/ ?num_mother ?denom_mother) 
				(+ (/ ?num_daughter1 ?denom_daughter1) (/ ?num_daughter2 ?denom_daughter2))
			)
			(< (/ ?num_mother ?denom_mother) 
				(abs (- (/ ?num_daughter1 ?denom_daughter1) (/ ?num_daughter2 ?denom_daughter2)))
			)
		)
	then
		(if (is-qn-conserved "isospin")
			(retract ?decay)
	  	else
	  		(modify ?decay (violating_quantum_number_list ?violating_quantum_number_list "isospin"))
	  	)
		;(printout t "decay violates isospin!" crlf)
	)
)