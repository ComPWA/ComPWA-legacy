;;;*************************
;;;* DECAY VIOLATION RULES *
;;;*************************

(defrule remove-bad-L
	(declare (salience 99))
	?decay <- (Decay 
				(quantum_number_name ?qn_name) 
				(mother ?mother_id) 
				(daughters ?daughter1_id ?daughter2_id $?others)
				(required_variable_names "angular-momentum" $?other_rvns) 
				(violating_quantum_number_list $?violating_quantum_number_list)
			  )
	(test (not (= 0 (str-compare ?qn_name "spinwave"))))
	=>
	(bind ?angular_momentum (get-spin-qn-with-unique-id (get-required-variable "angular-momentum" ?decay)))
	; in the helicity formalism we throw out all z components of L except 0
	(if (<> 0 (fact-slot-value ?angular_momentum z_component_numerator))
	then
		;(printout t "removing " ?decay crlf)
		(retract ?decay)
    )
)

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
	(test (not (member$ "spin" ?violating_quantum_number_list)))
	=>
	;get the required information
	(bind ?angular_momentum (get-spin-qn-with-unique-id (get-required-variable "angular-momentum" ?decay)))
	(bind ?L (/ (fact-slot-value ?angular_momentum numerator) (fact-slot-value ?angular_momentum denominator)))
	
	
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

(defrule check-helicity
	(declare (salience 99))
	(SpinQuantumNumber (unique_id ?mother_id) (numerator ?num_mother) (denominator ?denom_mother))
	(SpinQuantumNumber (unique_id ?daughter1_id) (z_component_numerator ?z_num_daughter1) (denominator ?denom_daughter1))
	(SpinQuantumNumber (unique_id ?daughter2_id) (z_component_numerator ?z_num_daughter2) (denominator ?denom_daughter2))
	?decay <- (Decay (quantum_number_name "spin") (mother ?mother_id) (daughters ?daughter1_id ?daughter2_id $?others))
	=>
	(if (< (/ ?num_mother ?denom_mother) (abs (+ (/ ?z_num_daughter1 ?denom_daughter1) (/ ?z_num_daughter2 ?denom_daughter2))))
	then
	  ;(printout t "decay violates angular momentum conservation!" crlf)
	  (retract ?decay)
	)
)


;(deffunction find-related-decay (?decay)
;	;invert the helicities of the decay
;	
;	
;	(bind ?decay_relation_list (find-all-facts ((?f NameList)) (= 0 (str-compare ?f:name "DecayRelationList"))))
;	
;	
;	;go through the list and try to find a match
;	
;)

; in the helicity formalism we can relate amplitudes from parity conservation
;(defrule parity-relate-helicity-amplitudes
;	(declare (salience 99))
;	?decay <- (Decay (quantum_number_name "spinwave") (mother ?mother_index) 
;				(daughters ?d1_index ?d2_index $?others))
;	?decay_tree <- (DecayTree (decays ?decay $?decays) (all_occuring_waves $?all_occuring_waves))
;	
;	=>
;	
;	;lets get the quantum numbers
;	(bind ?parity_mother (get-qn-value ?mother_index ?decay_tree "parity"))
;	(bind ?parity_d1 (get-qn-value ?d1_index ?decay_tree "parity"))
;	(bind ?parity_d2 (get-qn-value ?d2_index ?decay_tree "parity"))
;	(bind ?spin_mother 
;		(/ 
;			(fact-slot-value (get-spin-qn-with-unique-id (get-qn-value ?mother_index ?decay_tree "spin")) numerator)
;			(fact-slot-value (get-spin-qn-with-unique-id (get-qn-value ?mother_index ?decay_tree "spin")) denominator)
;		)
;	)
;	(bind ?spin_d1 
;		(/ 
;			(fact-slot-value (get-spin-qn-with-unique-id (get-qn-value ?d1_index ?decay_tree "spin")) numerator)
;			(fact-slot-value (get-spin-qn-with-unique-id (get-qn-value ?d1_index ?decay_tree "spin")) denominator)
;		)
;	)
;	(bind ?spin_d2 
;		(/ 
;			(fact-slot-value (get-spin-qn-with-unique-id (get-qn-value ?d2_index ?decay_tree "spin")) numerator)
;			(fact-slot-value (get-spin-qn-with-unique-id (get-qn-value ?d2_index ?decay_tree "spin")) denominator)
;		)
;	)
;	(bind ?helicity_d1 
;		(/ 
;			(fact-slot-value (get-spin-qn-with-unique-id (get-qn-value ?d1_index ?decay_tree "spin")) z_component_numerator)
;			(fact-slot-value (get-spin-qn-with-unique-id (get-qn-value ?d1_index ?decay_tree "spin")) denominator)
;		)
;	)
;	(bind ?helicity_d2 
;		(/ 
;			(fact-slot-value (get-spin-qn-with-unique-id (get-qn-value ?d2_index ?decay_tree "spin")) z_component_numerator)
;			(fact-slot-value (get-spin-qn-with-unique-id (get-qn-value ?d2_index ?decay_tree "spin")) denominator)
;		)
;	)
;
;	;now try to find the other decay in the relation storage
;	(bind ?make related_decay (find-related-decay ?decay))
;	(if (not FALSE)
;	then
;	 	;if we found one make the relation
;		;compute relation to other decay
;	else
;		;otherwise just add this decay
;		(bind ?decay_relation_list (find-all-facts ((?f NameList)) (= 0 (str-compare ?f:name "DecayRelationList"))))
;	)
;)