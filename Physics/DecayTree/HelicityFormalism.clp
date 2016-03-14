;;;*************************
;;;* DECAY VIOLATION RULES *
;;;*************************

(defrule check-spin
	(declare (salience 99))
	(SpinQuantumNumber (unique_id ?mother_id) (numerator ?num_mother) (denominator ?denom_mother))
	(SpinQuantumNumber (unique_id ?daughter1_id) (numerator ?num_daughter1) (denominator ?denom_daughter1))
	(SpinQuantumNumber (unique_id ?daughter2_id) (numerator ?num_daughter2) (denominator ?denom_daughter2))
	?decay <- (Decay 
				(quantum_number_name "spin") 
				(mother ?mother_id) 
				(daughters ?daughter1_id ?daughter2_id $?others)
				(required_variable_names "angular-momentum" $?other_rvns) 
				(violating_quantum_number_list $?violating_quantum_number_list)
			  )	
	=>
	;get the required information
	(bind ?angular_momentum (get-spin-qn-with-unique-id (get-required-variable "angular-momentum" ?decay)))
	(bind ?L (/ (fact-slot-value ?angular_momentum numerator) (fact-slot-value ?angular_momentum denominator)))
	
	; in the helicity formalism we throw out all z components of L except 0
	(if (<> 0 (fact-slot-value ?angular_momentum z_component_numerator))
	then
		(retract ?decay)
	else
	
	(bind ?violated TRUE)
	(bind ?comb1 (create$))
	(loop-for-count 
		(?i 
			(abs (- (/ ?num_daughter1 ?denom_daughter1) (/ ?num_daughter2 ?denom_daughter2)))
			(+ (/ ?num_daughter1 ?denom_daughter1) (/ ?num_daughter2 ?denom_daughter2))
		)
	do
		(bind ?comb1 (insert$ ?comb1 1 ?i))
	)
	;(printout t ?L " " ?comb1 crlf)
	;(printout t ?L " " (/ ?num_mother ?denom_mother) " " (/ ?num_daughter1 ?denom_daughter1) " "  (/ ?num_daughter2 ?denom_daughter2) crlf)
	(foreach ?val ?comb1
	do
		(if (and
				(<= (/ ?num_mother ?denom_mother) 
					(+ ?val ?L)
				)
				(>= (/ ?num_mother ?denom_mother) 
					(abs (- ?val ?L))
				)
			)
		then
			(bind ?violated FALSE)
			(break)
		)
	)
	
	(if ?violated
	then
		(if (is-qn-conserved "spin")
		then
			;(printout t "decay violates angular momentum conservation!" crlf)
			(retract ?decay)
	  	else
	  		(modify ?decay (violating_quantum_number_list ?violating_quantum_number_list "spin"))
	  	)
	)
	)
)

(defrule check-helicity
	(declare (salience 99))
	(SpinQuantumNumber (unique_id ?mother_id) (numerator ?num_mother) (denominator ?denom_mother))
	(SpinQuantumNumber (unique_id ?daughter1_id) (z_component_numerator ?z_num_daughter1) (denominator ?denom_daughter1))
	(SpinQuantumNumber (unique_id ?daughter2_id) (z_component_numerator ?z_num_daughter2) (denominator ?denom_daughter2))
	?decay <- (Decay (quantum_number_name "spin") (mother ?mother_id) (daughters ?daughter1_id ?daughter2_id $?others))
	=>
	(if (< (/ ?num_mother ?denom_mother) 
			(abs (- (/ ?z_num_daughter1 ?denom_daughter1) (/ ?z_num_daughter2 ?denom_daughter2)))
		)
	then
	  ;(printout t "decay violates angular momentum conservation!" crlf)
	  (retract ?decay)
	)
)