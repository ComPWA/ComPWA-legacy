;;;*************
;;;* TEMPLATES *
;;;*************

; define particle template
;(deftemplate Particle
;	(slot unique_id (type INTEGER))
;	(slot name (type STRING))
;	(slot pid (type INTEGER))
;	(slot mass (type FLOAT))
;	(slot charge (type INTEGER))
;	(slot isospin (type INTEGER))
;	(slot isospin_z (type INTEGER))
;	(slot spin (type INTEGER))
;	(slot parity (type INTEGER))
;	(slot cparity (type INTEGER))
;)

; define helicity wave final or initial state template
;(deftemplate SpinWaveMultiplet
;    (slot unique_id (type INTEGER))
;	(slot charge (type INTEGER))
;	(slot isospin_num (type INTEGER))
;	(slot isospin_denom (type INTEGER))
;	(slot isospin_z_num (type INTEGER))
;	(slot spin_num (type INTEGER))
;	(slot spin_denom (type INTEGER))
;	; we have multislot of spin z to allow for specific components
;	; in the initial or final state
;	(multislot spin_z_num)
;	(slot parity (type INTEGER))
;	(slot cparity (type INTEGER))
;)
	
; define spin wave	
(deftemplate SpinWave
    (slot unique_id (type INTEGER))
	(slot charge (type INTEGER))
	(slot isospin_num (type INTEGER))
	(slot isospin_denom (type INTEGER))
	(slot isospin_z_num (type INTEGER))
	(slot spin_num (type INTEGER))
	(slot spin_denom (type INTEGER))
	(slot spin_z_num (type INTEGER))
	(slot parity (type INTEGER))
	(slot cparity (type INTEGER))
)
	
; define two body decay template
(deftemplate TwoBodyDecay
	(slot mother)
	(slot daughter1)
	(slot daughter2)
)

; define two body decay sequence or tree
;(deftemplate TwoBodyDecayTree
;    (multislot final_state_multiplets)
;	(multislot available_particles)
;	(multislot two_body_decays)
;)

; allowed intermediate state spins
;(deftemplate AllowedQN
;	(multislot spin_nums)
;	(slot spin_denom (type INTEGER))
;	(multislot isospin_nums)
;	(slot isospin_denom (type INTEGER))
;	(multislot charge)
;	(multislot parity)
;	(multislot cparity)
;)
	
; seed decay trees template
;(deftemplate SeedTwoBodyDecayTrees
;	(multislot two_body_decay_trees)
;)
	
; result template
(deftemplate ViolatingRulesForDecay
	(multislot list_of_violated_rules)
)
	
;;;*******************
;;;* USER CONDITIONS *
;;;*******************

; allowed intermediate state spins
;(deffacts user-conditions
;	(AllowedQN
;		(spin_nums 0 1 2) (spin_denom 1 ) (isospin_nums 0 1) (isospin_denom 1)
;		(charge 0) (parity -1 1) (cparity -1 1)
;	)
;)

;;;*****************
;;;* INITIAL STATE *
;;;*****************

;(defglobal
;   ?*total_unique_id_counter* = 1
;)

;(deffacts initial-state
;	(SpinWaveMultiplet (unique_id 0) (spin_num 1) (spin_denom 1) (spin_z_num -1 1)
;		(isospin_num 0) (isospin_denom 1) (isospin_z_num 0)
;	)
;)

;(deffacts final-state-list
;	(SpinWaveMultiplet (unique_id 1) (spin_num 1) (spin_denom 1) (spin_z_num -1 1)
;		(isospin_num 0) (isospin_denom 1) (isospin_z_num 0)
;	)
;	(SpinWaveMultiplet (unique_id 2) (spin_num 0) (spin_denom 1) (spin_z_num 0)
;		(isospin_num 1) (isospin_denom 1) (isospin_z_num 0)
;	)
;	(SpinWaveMultiplet (unique_id 3) (spin_num 0) (spin_denom 1) (spin_z_num 0)
;		(isospin_num 1) (isospin_denom 1) (isospin_z_num 0)
;	)
;	(SpinWaveMultiplet (unique_id 4) (spin_num 0) (spin_denom 1) (spin_z_num 0)
;		(isospin_num 1) (isospin_denom 1) (isospin_z_num 0)
;	)
;)

;(deffacts initial-state
;  (Particle (unique_id 10) (name "jpsi") (pid 110) (mass 3.097) (parity -1)))

;(deffacts final-state-list
;  (Particle (unique_id 1) (name "gamma") (pid 22) (mass 0.0) (parity -1))
;  (Particle (unique_id 2) (name "pi0") (pid 110) (mass 0.135) (parity -1))
;  (Particle (unique_id 3) (name "pi0") (pid 110) (mass 0.135) (parity -1)))
  
;(deffacts intermediate-state-list
;  (Particle (unique_id 4) (name "f0") (pid 22) (mass 0.99) (parity 1))
;  (Particle (unique_id 5) (name "f2") (pid 110) (mass 1.275) (parity 1))
;  (Particle (unique_id 6) (name "omega") (pid 110) (mass 0.783) (parity -1)))

;(deffacts seed-decay-tree
;  (TwoBodyDecayTree))

;;;***********
;;;* RESULTS *
;;;***********

(deffacts violating-rules-list-for-decay
	(ViolatingRulesForDecay))
	

;;;*************
;;;* FUNCTIONS *
;;;*************

;(deffunction are-spins-equal (?first ?second)
;	(and
;		(= (fact-slot-value ?first spin_num) (fact-slot-value ?second spin_num))
;		(= (fact-slot-value ?first spin_denom) (fact-slot-value ?second spin_denom))
;		(= (fact-slot-value ?first spin_z_num) (fact-slot-value ?second spin_z_num))
;	)
;)
		
	
;;;*******************************
;;;* DECAY TREE GENERATION RULES *
;;;*******************************

; generate all possible allowed two body decays
;(defrule create-reverse-decay
;	(Particle (unique_id ?mymother))
;	(Particle (unique_id ?mydaughter1))
;	(Particle (unique_id ?mydaughter2))
;	=>
;	(if (and (<> ?mymother ?mydaughter1) (<> ?mymother ?mydaughter2) (<> ?mydaughter1 ?mydaughter2))
;	then
;	  (assert(TwoBodyDecay (mother ?mymother) (daughter1 ?mydaughter1) (daughter2 ?mydaughter2)))))
	
; generate all possible two body decay trees
;(defrule grow-decay-tree
;	?decay <- (TwoBodyDecay (mother ?mymother) (daughter1 ?mydaughter1) (daughter2 ?mydaughter2))
;	?decay_tree <- (TwoBodyDecayTree (available_particles $?remaining_particles ?mydaughter1 ?mydaughter2) (two_body_decays $?mytwo_body_decays))
;	=>
;	(duplicate ?decay_tree (two_body_decays ?mytwo_body_decays ?decay)
;						   (available_particles ?remaining_particles ?mymother)))

; filter out decay trees with correct initial state
;(defrule filter-decay-trees
;	?decay_tree <- (TwoBodyDecayTree (available_particles $?remaining_particles ?topnode))
;	(test (and (= ?*initial_state_unique_id* ?topnode) (= (length ?remaining_particles) 0)))
;   ?valid_two_body_decay_trees <- (ValidTwoBodyDecayTrees (two_body_decay_trees $?existing_two_body_decay_trees))
;	(test (not (member$ ?decay_tree ?existing_two_body_decay_trees)))
;	=>
;	  (modify ?valid_two_body_decay_trees (two_body_decay_trees ?decay_tree ?existing_two_body_decay_trees)))
	
	
	

	
; create all spin waves
;(defrule create-initial-spin-waves
;	(AllowedQN
;		(spin_nums $?spin_nums) (spin_denom ?spin_denom) 
;		(isospin_nums $?isospin_nums) (isospin_denom ?isospin_denom)
;		(charge $?charges)
;		(parity $?parities)
;		(cparity $?cparities)
;	)
;	=>
;	(foreach ?charge ?charges
;	(foreach ?parity ?parities
;	(foreach ?cparity ?cparities
;	
;	(foreach ?isospin_num ?isospin_nums
;		(bind ?isospin_z_num (* -1 ?isospin_num))
;		(while (<= ?isospin_z_num ?isospin_num)
;			(foreach ?spin_num ?spin_nums
;		   		(bind ?spin_z_num (* -1 ?spin_num))
;				(while (<= ?spin_z_num ?spin_num) 
;					(assert
;						(SpinWave (unique_id ?*total_unique_id_counter*) 
;						(spin_num ?spin_num) (spin_denom ?spin_denom) (spin_z_num ?spin_z_num)
;						(isospin_num ?isospin_num) (isospin_denom ?isospin_denom) (isospin_z_num ?isospin_z_num)
;						(charge ?charge) (parity ?parity) (cparity ?cparity)
;						)
;					)
;					(bind ?*total_unique_id_counter* (+ ?*total_unique_id_counter* 1))
;					(bind ?spin_z_num (+ ?spin_z_num ?spin_denom))
;				)
;			)
;			(bind ?isospin_z_num (+ ?isospin_z_num ?isospin_denom))
;		)
;	)
;	
;	)
;	)
;	)
;)

; generate seed decay trees
;(defrule create-seed-decay-trees
;	?condition <- (SpinWaveMultiplet (unique_id ?unique_id) 
;		(spin_num ?spin_num) (spin_denom ?spin_denom) (spin_z_num $?other_spin_z ?spin_z_num)
;		(isospin_num ?isospin_num) (isospin_denom ?isospin_denom) (isospin_z_num ?isospin_z_num)
;		)
;	; multiplet unique id 0 is the top node, which we do not need
;	(test (<> ?unique_id 0))
;	(SpinWave (unique_id ?id) 
;		(spin_num ?spin_num) (spin_denom ?spin_denom) (spin_z_num ?spin_z_num)
;		(isospin_num ?isospin_num) (isospin_denom ?isospin_denom) (isospin_z_num ?isospin_z_num)
;	)
;	?decay_tree <- (TwoBodyDecayTree (final_state_multiplets $?final_state_multiplets) (available_particles $?available_particles))
;	(test (not (member$ ?unique_id ?final_state_multiplets)))
;	=>
;	(duplicate ?decay_tree (available_particles ?available_particles ?id)
;		(final_state_multiplets ?final_state_multiplets ?unique_id)
;	)
;	(modify ?condition (spin_z_num $?other_spin_z))
;)

	
; generate all possible allowed two body decays
;(defrule create-two-body-decay
;	(SpinWave (unique_id ?mymother))
;	(SpinWave (unique_id ?mydaughter1))
;	(SpinWave (unique_id ?mydaughter2))
;	=>
;	(assert(TwoBodyDecay (mother ?mymother) (daughter1 ?mydaughter1) (daughter2 ?mydaughter2)))
;)
	
; generate all possible two body decay trees
;(defrule grow-decay-tree
;	?decay <- (TwoBodyDecay (mother ?mymother) (daughter1 ?mydaughter1) (daughter2 ?mydaughter2))
;	?decay_tree <- (TwoBodyDecayTree (available_particles $?remaining_particles ?mydaughter1 ?mydaughter2) (two_body_decays $?mytwo_body_decays))
;	=>
;	(duplicate ?decay_tree (two_body_decays ?mytwo_body_decays ?decay)
;						   (available_particles ?remaining_particles ?mymother))
;)

; filter out decay trees with correct initial state
;(defrule filter-decay-trees
;	(SpinWaveMultiplet (unique_id 0) 
;		(spin_num ?spin_num) (spin_denom ?spin_denom) (spin_z_num $?other_spin_z ?spin_z_num)
;		(isospin_num ?isospin_num) (isospin_denom ?isospin_denom) (isospin_z_num ?isospin_z_num)
;	)
;	(SpinWave (unique_id ?id) 
;		(spin_num ?spin_num) (spin_denom ?spin_denom) (spin_z_num ?spin_z_num)
;		(isospin_num ?isospin_num) (isospin_denom ?isospin_denom) (isospin_z_num ?isospin_z_num)
;	)
;	?decay_tree <- (TwoBodyDecayTree (available_particles $?remaining_particles ?id))
;	(test (= (length ?remaining_particles) 0))
;    ?valid_two_body_decay_trees <- (ValidTwoBodyDecayTrees (two_body_decay_trees $?existing_two_body_decay_trees))
;	(test (not (member$ ?decay_tree ?existing_two_body_decay_trees)))
;	=>
;		(modify ?valid_two_body_decay_trees (two_body_decay_trees ?decay_tree ?existing_two_body_decay_trees))
;)


;(defrule printresult
;	(ValidTwoBodyDecayTrees (two_body_decay_trees $?existing_two_body_decay_trees ?tbdt))
;	=>
;	(foreach ?tbdt ?existing_two_body_decay_trees
;		(printout t ?tbdt crlf)
;	)
;)

;;;*************************
;;;* DECAY VIOLATION RULES *
;;;*************************

(defrule check-isospin
	(declare (salience 99))
	(SpinWave (unique_id ?mother_id) (isospin_num ?isospin_num_mother) (isospin_denom ?isospin_denom_mother))
	(SpinWave (unique_id ?daughter1_id) (isospin_num ?isospin_num_daughter1) (isospin_denom ?isospin_denom_daughter1))
	(SpinWave (unique_id ?daughter2_id) (isospin_num ?isospin_num_daughter2) (isospin_denom ?isospin_denom_daughter2))
	?decay <- (TwoBodyDecay (mother ?mother_id) (daughter1 ?daughter1_id) (daughter2 ?daughter2_id))
	?violated_rules <- (ViolatingRulesForDecay (list_of_violated_rules $?violated_rule_ids))
	(test (not (member$ 1 ?violated_rule_ids)))
	=>
	(if (or (> (/ ?isospin_num_mother ?isospin_denom_mother) 
				(+ (/ ?isospin_num_daughter1 ?isospin_denom_daughter1) (/ ?isospin_num_daughter2 ?isospin_denom_daughter2))
			)
			(< (/ ?isospin_num_mother ?isospin_denom_mother) 
				(abs (- (/ ?isospin_num_daughter1 ?isospin_denom_daughter1) (/ ?isospin_num_daughter2 ?isospin_denom_daughter2)))
			)
		)
	then
	  ;(printout t "decay violates isospin conservation!" crlf)
	  ;(retract ?decay)
	  (modify ?violated_rules (list_of_violated_rules 1 ?violated_rule_ids))
	)
)

(defrule check-isospin-z
	(declare (salience 99))
	(SpinWave (unique_id ?mother_id) (isospin_z_num ?isospin_z_num_mother) (isospin_denom ?isospin_denom_mother))
	(SpinWave (unique_id ?daughter1_id) (isospin_z_num ?isospin_z_num_daughter1) (isospin_denom ?isospin_denom_daughter1))
	(SpinWave (unique_id ?daughter2_id) (isospin_z_num ?isospin_z_num_daughter2) (isospin_denom ?isospin_denom_daughter2))
	?decay <- (TwoBodyDecay (mother ?mother_id) (daughter1 ?daughter1_id) (daughter2 ?daughter2_id))
	?violated_rules <- (ViolatingRulesForDecay (list_of_violated_rules $?violated_rule_ids))
	(test (not (member$ 2 ?violated_rule_ids)))
	=>
	(if (<> (/ ?isospin_z_num_mother ?isospin_denom_mother) 
			(+ (/ ?isospin_z_num_daughter1 ?isospin_denom_daughter1) (/ ?isospin_z_num_daughter2 ?isospin_denom_daughter2))
		)
	then
	  ;(printout t "decay violates isospin z conservation!" crlf)
	  ;(retract ?decay)
	  (modify ?violated_rules (list_of_violated_rules 2 ?violated_rule_ids))
	)
)


(defrule check-charge
	(declare (salience 99))
	(SpinWave (unique_id ?mother_id) (charge ?charge_mother))
	(SpinWave (unique_id ?daughter1_id) (charge ?charge_daughter1))
	(SpinWave (unique_id ?daughter2_id) (charge ?charge_daughter2))
	?decay <- (TwoBodyDecay (mother ?mother_id) (daughter1 ?daughter1_id) (daughter2 ?daughter2_id))
	?violated_rules <- (ViolatingRulesForDecay (list_of_violated_rules $?violated_rule_ids))
	(test (not (member$ 3 ?violated_rule_ids)))
	=>
	(if (<> ?charge_mother (+ ?charge_daughter1 ?charge_daughter2))
	then
	  ;(printout t "decay charge conservation!" crlf)
	  ;(retract ?decay)
	  (modify ?violated_rules (list_of_violated_rules 3 ?violated_rule_ids))
	)
)

(defrule check-parity
	(declare (salience 99))
	(SpinWave (unique_id ?mother_id) (parity ?parity_mother))
	(SpinWave (unique_id ?daughter1_id) (parity ?parity_daughter1))
	(SpinWave (unique_id ?daughter2_id) (parity ?parity_daughter2))
	?decay <- (TwoBodyDecay (mother ?mother_id) (daughter1 ?daughter1_id) (daughter2 ?daughter2_id))
	?violated_rules <- (ViolatingRulesForDecay (list_of_violated_rules $?violated_rule_ids))
	(test (not (member$ 4 ?violated_rule_ids)))
	=>
	(if (<> ?parity_mother (* ?parity_daughter1 ?parity_daughter2))
	then
	  ;(printout t "decay violates parity!" crlf)
	  ;(retract ?decay)
	  (modify ?violated_rules (list_of_violated_rules 4 ?violated_rule_ids))
	)
)

(defrule check-cparity
	(declare (salience 99))
	(SpinWave (unique_id ?mother_id) (cparity ?cparity_mother) (charge ?charge_mother))
	(SpinWave (unique_id ?daughter1_id) (cparity ?cparity_daughter1) (charge ?charge_daughter1))
	(SpinWave (unique_id ?daughter2_id) (cparity ?cparity_daughter2) (charge ?charge_daughter2))
	?decay <- (TwoBodyDecay (mother ?mother_id) (daughter1 ?daughter1_id) (daughter2 ?daughter2_id))
	?violated_rules <- (ViolatingRulesForDecay (list_of_violated_rules $?violated_rule_ids))
	(test (not (member$ 5 ?violated_rule_ids)))
	=>
	(if (and (= 0 ?charge_mother ?charge_daughter1 ?charge_daughter2) (<> ?cparity_mother (* ?cparity_daughter1 ?cparity_daughter2)))
	then
	  ;(printout t "decay violates c-parity!" crlf)
	  ;(retract ?decay))
	  (modify ?violated_rules (list_of_violated_rules 5 ?violated_rule_ids))
	)
)

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
	  
