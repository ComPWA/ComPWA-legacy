;;;*************************
;;;* DECAY VIOLATION RULES *
;;;*************************


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

(defrule check-isospin-z
	(declare (salience 99))
	(SpinQuantumNumber (unique_id ?mother_id) (z_component_numerator ?z_num_mother) (denominator ?denom_mother))
	(SpinQuantumNumber (unique_id ?daughter1_id) (z_component_numerator ?z_num_daughter1) (denominator ?denom_daughter1))
	(SpinQuantumNumber (unique_id ?daughter2_id) (z_component_numerator ?z_num_daughter2) (denominator ?denom_daughter2))
	?decay <- (Decay (quantum_number_name "isospin") (mother ?mother_id) (daughters ?daughter1_id ?daughter2_id $?others))
	=>
	(if (<> (/ ?z_num_mother ?denom_mother) 
			(+ (/ ?z_num_daughter1 ?denom_daughter1) (/ ?z_num_daughter2 ?denom_daughter2))
		)
	then
	  ;(printout t "decay violates isospin z conservation!" crlf)
	  (retract ?decay)
	)
)

(defrule check-helicity
	(declare (salience 99))
	(SpinQuantumNumber (unique_id ?mother_id) (numerator ?num_mother) (denominator ?denom_mother))
	(SpinQuantumNumber (unique_id ?daughter1_id) (z_component_numerator ?z_num_daughter1) (denominator ?denom_daughter1))
	(SpinQuantumNumber (unique_id ?daughter2_id) (z_component_numerator ?z_num_daughter2) (denominator ?denom_daughter2))
	?decay <- (Decay (quantum_number_name ?name) (mother ?mother_id) (daughters ?daughter1_id ?daughter2_id $?others))
	(test (or (= (str-compare ?name "spin") 0) (= (str-compare ?name "isospin") 0)))
	=>
	(if (< (/ ?num_mother ?denom_mother) 
			(abs (- (/ ?z_num_daughter1 ?denom_daughter1) (/ ?z_num_daughter2 ?denom_daughter2)))
		)
	then
	  ;(printout t "decay violates angular momentum conservation!" crlf)
	  (retract ?decay)
	)
)

(defrule check-parity
	(declare (salience 99))
	?decay <- (Decay (quantum_number_name "parity") (mother ?parity_mother) (daughters ?parity_daughter1 ?parity_daughter2 $?others))
	=>
	(if (<> ?parity_mother (* ?parity_daughter1 ?parity_daughter2))
	then
	  ;(printout t "decay violates parity!" crlf)
	  (retract ?decay)
	)
)

(defrule check-cparity
	(declare (salience 99))
	?decay <- (Decay (quantum_number_name "cparity") (mother ?cparity_mother) (daughters ?cparity_daughter1 ?cparity_daughter2 $?others))
	=>
	(if (<> ?cparity_mother (* ?cparity_daughter1 ?cparity_daughter2))
	then
	  ;(printout t "decay violates cparity!" crlf)
	  (retract ?decay)
	)
)