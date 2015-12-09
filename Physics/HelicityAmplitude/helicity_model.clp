;;;*************
;;;* TEMPLATES *
;;;*************

; define particle template
(deftemplate Particle
	(slot unique_id (type INTEGER))
	(slot name (type STRING))
	(slot pid (type INTEGER))
	(slot mass (type FLOAT))
	(slot charge (type INTEGER))
	(slot isospin (type INTEGER))
	(slot isospin_z (type INTEGER))
	(slot spin (type INTEGER))
	(slot parity (type INTEGER))
	(slot cparity (type INTEGER)))
	
; define two body decay template
(deftemplate TwoBodyDecay
	(slot mother)
	(slot daughter1)
	(slot daughter2))

; define two body decay sequence or tree
(deftemplate TwoBodyDecayTree
	(multislot available_particles)
	(multislot two_body_decays))
	
; result template
(deftemplate ValidTwoBodyDecayTrees
	(multislot two_body_decay_trees))
	
;;;*************
;;;* FUNCTIONS *
;;;*************

; print decay tree function

;;;*****************
;;;* INITIAL STATE *
;;;*****************

(defglobal
   ?*initial_state_unique_id* = 10
)

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
;  (TwoBodyDecayTree 
;  (available_particles 1 2 3)))

;;;***********
;;;* RESULTS *
;;;***********

(deffacts valid-decay-trees
	(ValidTwoBodyDecayTrees))

;;;*******************************
;;;* DECAY TREE GENERATION RULES *
;;;*******************************

; generate all possible allowed two body decays
(defrule create-reverse-decay
	(Particle (unique_id ?mymother))
	(Particle (unique_id ?mydaughter1))
	(Particle (unique_id ?mydaughter2))
	=>
	(if (and (<> ?mymother ?mydaughter1) (<> ?mymother ?mydaughter2) (<> ?mydaughter1 ?mydaughter2))
	then
	  (assert(TwoBodyDecay (mother ?mymother) (daughter1 ?mydaughter1) (daughter2 ?mydaughter2)))))
	
; generate all possible two body decay trees
(defrule grow-decay-tree
	?decay <- (TwoBodyDecay (mother ?mymother) (daughter1 ?mydaughter1) (daughter2 ?mydaughter2))
	?decay_tree <- (TwoBodyDecayTree (available_particles $?remaining_particles ?mydaughter1 ?mydaughter2) (two_body_decays $?mytwo_body_decays))
	=>
	(duplicate ?decay_tree (two_body_decays ?mytwo_body_decays ?decay)
						   (available_particles ?remaining_particles ?mymother)))


; filter out decay trees with correct initial state
(defrule filter-decay-trees
	?decay_tree <- (TwoBodyDecayTree (available_particles $?remaining_particles ?topnode))
	(test (and (= ?*initial_state_unique_id* ?topnode) (= (length ?remaining_particles) 0)))
    ?valid_two_body_decay_trees <- (ValidTwoBodyDecayTrees (two_body_decay_trees $?existing_two_body_decay_trees))
	(test (not (member$ ?decay_tree ?existing_two_body_decay_trees)))
	=>
	  (modify ?valid_two_body_decay_trees (two_body_decay_trees ?decay_tree ?existing_two_body_decay_trees))))

;;;*************************
;;;* DECAY VIOLATION RULES *
;;;*************************

(defrule check-mass
	?mymother <- (Particle (mass ?mass_mother))
	?mydaughter1 <- (Particle (mass ?mass_daughter1))
	?mydaughter2 <- (Particle (mass ?mass_daughter2))
	?decay <- (TwoBodyDecay (mother ?mymother) (daughter1 ?mydaughter1) (daughter2 ?mydaughter2))
	=>
	(if (<= ?mass_mother (+ ?mass_daughter1 ?mass_daughter2))
	then
	  (printout t "decay violates mass/energy conservation!" crfl)
	  (retract ?decay)))
	  
(defrule check-cparity
	?mymother <- (Particle (cparity ?cparity_mother) (charge ?charge_mother))
	?mydaughter1 <- (Particle (cparity ?cparity_daughter1) (charge ?charge_daughter1))
	?mydaughter2 <- (Particle (cparity ?cparity_daughter2) (charge ?charge_daughter2))
	?decay <- (TwoBodyDecay (mother ?mymother) (daughter1 ?mydaughter1) (daughter2 ?mydaughter2))
	=>
	(if (and (= 0 ?charge_mother ?charge_daughter1 ?charge_daughter2) (<> ?cparity_mother (* ?cparity_daughter1 ?cparity_daughter2)))
	then
	  (printout t "decay violates c-parity!" crfl)
	  (retract ?decay)))
	