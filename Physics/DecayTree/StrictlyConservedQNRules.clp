;;;***********************************************************
;;;* Stricly Conserved Quantum Number Conservation Law Rules *
;;;***********************************************************

(defrule check-charge
	(declare (salience 99))
	?decay <- (Decay (quantum_number_name "charge") (mother ?charge_mother) (daughters ?charge_daughter1 ?charge_daughter2 $?others))
	=>
	(if (<> ?charge_mother (+ ?charge_daughter1 ?charge_daughter2))
	then
	  ;(printout t "decay charge conservation!" crlf)
	  (retract ?decay)
	)
)

(defrule check-identical-particle-symmetry
	(declare (salience 99))
	?decay <- (Decay (quantum_number_name "spinwave") (mother ?mother_index) 
				(daughters ?d1_index ?d2_index $?others))
	?decay_tree <- (DecayTree (decays ?decay $?decays) (all_occuring_waves $?all_occuring_waves))
	=>
	;(printout t "checking bose fermi statistics" crlf)
	; check for identical particles in decay
	;(if (not ?decay_tree)
	;then
		;(printout t "wtf steve" crlf)
		
		;	(bind ?mother_wave (get-wave-from-unique-index ?mother_index ?decay_tree))
		;	(bind ?d1_wave (get-wave-from-unique-index ?d1_index ?decay_tree))
		;	(bind ?d2_wave (get-wave-from-unique-index ?d2_index ?decay_tree))
			
		;	(printout t (fact-slot-value ?mother_wave quantum_number_names) " " (fact-slot-value ?mother_wave quantum_number_values) crlf)
		;	(printout t (fact-slot-value ?d1_wave quantum_number_names) " " (fact-slot-value ?d1_wave quantum_number_values) crlf)
		;	(printout t (fact-slot-value ?d2_wave quantum_number_names) " " (fact-slot-value ?d2_wave quantum_number_values) crlf)
		
		(if (check-spinwave-equality ?d1_index ?d2_index ?all_occuring_waves)
		then
			;(bind ?mother_wave (get-wave-from-unique-index ?mother_index ?decay_tree))
			;(bind ?d1_wave (get-wave-from-unique-index ?d1_index ?decay_tree))
			;(bind ?d2_wave (get-wave-from-unique-index ?d2_index ?decay_tree))
			;(printout t "FOUND IDENTIDAL PARTICLES -------------------------------------------------------------" crlf)
		
			;(printout t (fact-slot-value ?mother_wave quantum_number_names) " " (fact-slot-value ?mother_wave quantum_number_values) crlf)
			;(printout t (fact-slot-value ?d1_wave quantum_number_names) " " (fact-slot-value ?d1_wave quantum_number_values) crlf)
			;(printout t (fact-slot-value ?d2_wave quantum_number_names) " " (fact-slot-value ?d2_wave quantum_number_values) crlf)
		
			;(printout t (fact-slot-value ?decay required_variable_names) crlf)
			;(bind ?angular_momentum (get-spin-qn-with-unique-id (get-required-variable "angular-momentum" ?decay)))
			;(bind ?L (/ (fact-slot-value ?angular_momentum numerator) (fact-slot-value ?angular_momentum denominator)))
			;(printout t ?L crlf)
		
			; if they are bosons then check for symmetric wavefunction -> mother parity = +1 
			(if (is-boson ?d1_index ?decay_tree)
			then
				(if (= -1 (get-qn-value ?mother_index ?decay_tree "parity"))
				then 
	  				(printout t "decay violates bose statistics!" crlf)
	  				(retract ?decay_tree)
	  				(retract ?decay)
	  			)
	  		else	  	
	  			; if they are fermions then check for antisymmetric wavefunction -> mother parity = -1 
				(if (is-fermion ?d1_index ?decay_tree)
				then
					(if (= +1 (get-qn-value ?mother_index ?decay_tree "parity"))
					then 
	 					;(printout t "decay violates fermi statistics!" crlf)
	  					(retract ?decay_tree)
	 					(retract ?decay)
	  				)
	  			)
	  		)
		)
	;)
	;(printout t "exiting that check" crlf)
)