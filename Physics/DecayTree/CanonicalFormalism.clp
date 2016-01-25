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

;(defrule check-spin
;	(declare (salience 99))
;	(SpinQuantumNumber (unique_id ?mother_id) (numerator ?num_mother) (denominator ?denom_mother))
;	(SpinQuantumNumber (unique_id ?daughter1_id) (numerator ?num_daughter1) (denominator ?denom_daughter1))
;	(SpinQuantumNumber (unique_id ?daughter2_id) (numerator ?num_daughter2) (denominator ?denom_daughter2))
;	?decay <- (Decay (quantum_number_name ?name) (mother ?mother_id) (daughters ?daughter1_id ?daughter2_id $?others))
;	(test (or (= (str-compare ?name "spin") 0) (= (str-compare ?name "isospin") 0)))
;	=>
;	(if (or (> (/ ?num_mother ?denom_mother) 
;				(+ (/ ?num_daughter1 ?denom_daughter1) (/ ?num_daughter2 ?denom_daughter2))
;			)
;			(< (/ ?num_mother ?denom_mother) 
;				(abs (- (/ ?num_daughter1 ?denom_daughter1) (/ ?num_daughter2 ?denom_daughter2)))
; 			)
;		)
;	then
;	  ;(printout t "decay violates spin conservation!" crlf)
;	  (retract ?decay)
;	)
;)

(defrule check-spin-z
	(declare (salience 99))
	(SpinQuantumNumber (unique_id ?mother_id) (z_component_numerator ?z_num_mother) (denominator ?denom_mother))
	(SpinQuantumNumber (unique_id ?daughter1_id) (z_component_numerator ?z_num_daughter1) (denominator ?denom_daughter1))
	(SpinQuantumNumber (unique_id ?daughter2_id) (z_component_numerator ?z_num_daughter2) (denominator ?denom_daughter2))
	?decay <- (Decay (quantum_number_name ?name) (mother ?mother_id) (daughters ?daughter1_id ?daughter2_id $?others))
	(test (or (= (str-compare ?name "spin") 0) (= (str-compare ?name "isospin") 0)))
	=>
	(if (<> (/ ?z_num_mother ?denom_mother) 
			(+ (/ ?z_num_daughter1 ?denom_daughter1) (/ ?z_num_daughter2 ?denom_daughter2))
		)
	then
	  ;(printout t "decay violates isospin z conservation!" crlf)
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





;(defrule check-isospin
;	(declare (salience 99))
;	(SpinWave (unique_id ?mother_id) (isospin_num ?isospin_num_mother) (isospin_denom ?isospin_denom_mother))
;	(SpinWave (unique_id ?daughter1_id) (isospin_num ?isospin_num_daughter1) (isospin_denom ?isospin_denom_daughter1))
;	(SpinWave (unique_id ?daughter2_id) (isospin_num ?isospin_num_daughter2) (isospin_denom ?isospin_denom_daughter2))
;	?decay <- (TwoBodyDecay (mother ?mother_id) (daughter1 ?daughter1_id) (daughter2 ?daughter2_id))
;	?violated_rules <- (ViolatingRulesForDecay (list_of_violated_rules $?violated_rule_ids))
;	(test (not (member$ 1 ?violated_rule_ids)))
;	=>
;	(if (or (> (/ ?isospin_num_mother ?isospin_denom_mother) 
;				(+ (/ ?isospin_num_daughter1 ?isospin_denom_daughter1) (/ ?isospin_num_daughter2 ?isospin_denom_daughter2))
;			)
;			(< (/ ?isospin_num_mother ?isospin_denom_mother) 
;				(abs (- (/ ?isospin_num_daughter1 ?isospin_denom_daughter1) (/ ?isospin_num_daughter2 ?isospin_denom_daughter2)))
;			)
;		)
;	then
;	  ;(printout t "decay violates isospin conservation!" crlf)
;	  ;(retract ?decay)
;	  (modify ?violated_rules (list_of_violated_rules 1 ?violated_rule_ids))
;	)
;)

;(defrule check-isospin-z
;	(declare (salience 99))
;	(SpinWave (unique_id ?mother_id) (isospin_z_num ?isospin_z_num_mother) (isospin_denom ?isospin_denom_mother))
;	(SpinWave (unique_id ?daughter1_id) (isospin_z_num ?isospin_z_num_daughter1) (isospin_denom ?isospin_denom_daughter1))
;	(SpinWave (unique_id ?daughter2_id) (isospin_z_num ?isospin_z_num_daughter2) (isospin_denom ?isospin_denom_daughter2))
;	?decay <- (TwoBodyDecay (mother ?mother_id) (daughter1 ?daughter1_id) (daughter2 ?daughter2_id))
;	?violated_rules <- (ViolatingRulesForDecay (list_of_violated_rules $?violated_rule_ids))
;	(test (not (member$ 2 ?violated_rule_ids)))
;	=>
;	(if (<> (/ ?isospin_z_num_mother ?isospin_denom_mother) 
;			(+ (/ ?isospin_z_num_daughter1 ?isospin_denom_daughter1) (/ ?isospin_z_num_daughter2 ?isospin_denom_daughter2))
;		)
;	then
;	  ;(printout t "decay violates isospin z conservation!" crlf)
;	  ;(retract ?decay)
;	  (modify ?violated_rules (list_of_violated_rules 2 ?violated_rule_ids))
;	)
;)


;(defrule check-charge
;	(declare (salience 99))
;	(SpinWave (unique_id ?mother_id) (charge ?charge_mother))
;	(SpinWave (unique_id ?daughter1_id) (charge ?charge_daughter1))
;	(SpinWave (unique_id ?daughter2_id) (charge ?charge_daughter2))
;	?decay <- (TwoBodyDecay (mother ?mother_id) (daughter1 ?daughter1_id) (daughter2 ?daughter2_id))
;	?violated_rules <- (ViolatingRulesForDecay (list_of_violated_rules $?violated_rule_ids))
;	(test (not (member$ 3 ?violated_rule_ids)))
;	=>
;	(if (<> ?charge_mother (+ ?charge_daughter1 ?charge_daughter2))
;	then
;	  ;(printout t "decay charge conservation!" crlf)
;	  ;(retract ?decay)
;	  (modify ?violated_rules (list_of_violated_rules 3 ?violated_rule_ids))
;	)
;)

;(defrule check-parity
;	(declare (salience 99))
;	(SpinWave (unique_id ?mother_id) (parity ?parity_mother))
;	(SpinWave (unique_id ?daughter1_id) (parity ?parity_daughter1))
;	(SpinWave (unique_id ?daughter2_id) (parity ?parity_daughter2))
;	?decay <- (TwoBodyDecay (mother ?mother_id) (daughter1 ?daughter1_id) (daughter2 ?daughter2_id))
;	?violated_rules <- (ViolatingRulesForDecay (list_of_violated_rules $?violated_rule_ids))
;	(test (not (member$ 4 ?violated_rule_ids)))
;	=>
;	(if (<> ?parity_mother (* ?parity_daughter1 ?parity_daughter2))
;	then
;	  ;(printout t "decay violates parity!" crlf)
;	  ;(retract ?decay)
;	  (modify ?violated_rules (list_of_violated_rules 4 ?violated_rule_ids))
;	)
;)

;(defrule check-cparity
;	(declare (salience 99))
;	(SpinWave (unique_id ?mother_id) (cparity ?cparity_mother) (charge ?charge_mother))
;	(SpinWave (unique_id ?daughter1_id) (cparity ?cparity_daughter1) (charge ?charge_daughter1))
;	(SpinWave (unique_id ?daughter2_id) (cparity ?cparity_daughter2) (charge ?charge_daughter2))
;	?decay <- (TwoBodyDecay (mother ?mother_id) (daughter1 ?daughter1_id) (daughter2 ?daughter2_id))
;	?violated_rules <- (ViolatingRulesForDecay (list_of_violated_rules $?violated_rule_ids))
;	(test (not (member$ 5 ?violated_rule_ids)))
;	=>
;	(if (and (= 0 ?charge_mother ?charge_daughter1 ?charge_daughter2) (<> ?cparity_mother (* ?cparity_daughter1 ?cparity_daughter2)))
;	then
;	  ;(printout t "decay violates c-parity!" crlf)
;	  ;(retract ?decay))
;	  (modify ?violated_rules (list_of_violated_rules 5 ?violated_rule_ids))
;	)
;)

;(defrule check-helicity
;	?mymother <- (HelicityWave (spin ?spin_mother) (helicity ?helicity_mother))
;	?mydaughter1 <- (HelicityWave (helicity ?helicity_daughter1))
;	?mydaughter2 <- (HelicityWave (helicity ?helicity_daughter2))
;	?decay <- (TwoBodyDecay (mother ?mymother) (daughter1 ?mydaughter1) (daughter2 ?mydaughter2))
;	=>
;	(if (or (<> ?helicity_mother (- ?helicity_daughter1 ?helicity_daughter2)) (> (abs (- ?helicity_daughter1 ?helicity_daughter2)) ?spin_mother))
;	then
;	  (printout t "decay violates helicity conservation!" crfl)
;	  (retract ?decay)
;	)
;)

;(defrule check-mass
;	?mymother <- (Particle (mass ?mass_mother))
;	?mydaughter1 <- (Particle (mass ?mass_daughter1))
;	?mydaughter2 <- (Particle (mass ?mass_daughter2))
;	?decay <- (TwoBodyDecay (mother ?mymother) (daughter1 ?mydaughter1) (daughter2 ?mydaughter2))
;	=>
;	(if (<= ?mass_mother (+ ?mass_daughter1 ?mass_daughter2))
;	then
;	  (printout t "decay violates mass/energy conservation!" crfl)
;	  (retract ?decay))
;)
	  
