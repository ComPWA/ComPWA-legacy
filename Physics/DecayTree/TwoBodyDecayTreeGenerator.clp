;;;*******************************
;;;* DECAY TREE GENERATION RULES *
;;;*******************************


; generate all possible allowed decays
(defrule create-decays
	(declare (salience 10))
	(AllowedQuantumNumbers (name ?quantum_number_name) (values $?values))
	?quantum_number_name_list <- (QuantumNumberNameList (names $?current_names))
	=>
	(foreach ?mother_value ?values
		(foreach ?daughter1_value ?values
			(foreach ?daughter2_value ?values	
				(assert(Decay (quantum_number_name ?quantum_number_name) (mother ?mother_value) (daughters ?daughter1_value ?daughter2_value)))
			)
		)
	)
	(if (not (member$ ?quantum_number_name ?current_names))
	then
		(modify ?quantum_number_name_list (names ?current_names ?quantum_number_name))
	)
)

; generate seed decay trees
(defrule create-seed-decay-trees
	(InitialAndFinalState 
		(initial_state ?initial_state) 
		(final_state $?final_state)
	)
	=>
	(bind ?index_pair_list (create$))
	(bind ?availables_waves (create$))
	(bind ?counter 1)
	(foreach ?fs ?final_state
		(bind ?index_pair_list (insert$ ?index_pair_list 1 (assert-index-pair ?counter (member$ ?fs ?final_state))))
		(bind ?availables_waves (insert$ ?availables_waves 1 ?counter))
		(bind ?counter (+ ?counter 1))
	)
	(assert 
		(DecayTree 
			(initial_state_wave ?initial_state)
			(all_occuring_waves ?final_state)
			(available_waves ?availables_waves)
			(unique_index_wave_mapping ?index_pair_list)
		)
	)
)

	
; generate all possible two body decay trees
(defrule create-all-valid-two-body-decay-trees
	?decay_tree <- (DecayTree (available_waves $?available_waves) (decays $?decays)
					(initial_state_wave ?initial_state_wave) (all_occuring_waves $?all_occuring_waves)
					(unique_index_wave_mapping $?unique_index_wave_mapping))
	=>
	;(printout t "asdfasfd " ?all_occuring_waves crlf)
	(bind ?dummy_wave (assert (SpinWave)))
	(if (> (length ?available_waves) 2)
	then
	  (foreach ?wave_unique_id ?available_waves
	  	(bind ?wave (get-wave-from-unique-index ?wave_unique_id ?decay_tree))
	  	
	    (bind ?reduced_available_waves (subseq$ ?available_waves 
	    	(+ 1 (member$ ?wave_unique_id ?available_waves)) (length ?available_waves))
	    )
	  	(bind ?unique_available_waves (get-unique-available-waves ?reduced_available_waves ?decay_tree))
		
		(foreach ?wave2_unique_id ?unique_available_waves
			(bind ?wave2 (get-wave-from-unique-index ?wave2_unique_id ?decay_tree))
			
			(bind ?spin_waves (create$ ?dummy_wave))
			(bind ?temp_spin_waves (create$))
			(bind ?qn_name_list (get-list-of-qn-names ?wave ?wave2))
			
			;(printout t ?qn_name_list crlf)
			(foreach ?qn_name ?qn_name_list
				;(printout t ?qn_name crlf)
				;(printout t ?spin_waves crlf)
				;(printout t "wthasdfas " ?temp_spin_waves crlf)
				(foreach ?spin_wave ?spin_waves
					(bind ?single_qn_decays (find-all-facts ((?d Decay)) (= 0 (str-compare ?d:quantum_number_name ?qn_name))))
					;(printout t "we have "(length ?single_qn_decays)" decays" crlf)
					(foreach ?single_qn_decay ?single_qn_decays
						(if (is-decay-valid ?single_qn_decay ?wave ?wave2)
						then
							;(printout t (fact-slot-value ?spin_wave quantum_number_names) crlf)
							(bind ?spin_wave_fact TRUE)
							(if (not (fact-slot-value ?spin_wave quantum_number_names))
							then
							    (bind ?spin_wave_fact
									(duplicate ?spin_wave 
										(quantum_number_names ?qn_name) 
										(quantum_number_values (fact-slot-value ?single_qn_decay mother))
									)
								)
							else
								;(printout t "updating wave " (fact-slot-value ?spin_wave quantum_number_names) ?qn_name crlf)
								(bind ?spin_wave_fact
									(duplicate ?spin_wave 
										(quantum_number_names (fact-slot-value ?spin_wave quantum_number_names) ?qn_name) 
										(quantum_number_values (fact-slot-value ?spin_wave quantum_number_values) 
											(fact-slot-value ?single_qn_decay mother)
										)
									)
								)
							)
							(if (not ?spin_wave_fact)
							then
								(bind ?spin_wave_fact 
									(find-spinwave-fact (insert$ (fact-slot-value ?spin_wave quantum_number_names) 1 ?qn_name)
										(insert$ (fact-slot-value ?spin_wave quantum_number_values) 1 (fact-slot-value ?single_qn_decay mother))
									)
								)
								;(printout t (fact-slot-value ?spin_wave quantum_number_names) " " (fact-slot-value ?spin_wave quantum_number_values) crlf)
								;(printout t "this wave already exists in form of " ?spin_wave_fact crlf)
							)
							
							(bind ?temp_spin_waves 
								(insert$ ?temp_spin_waves 1 ?spin_wave_fact)
							)
						)
					)
					(retract ?spin_wave)
				)
				;(facts)
				;(printout t "wth " ?temp_spin_waves crlf)
				(bind ?spin_waves ?temp_spin_waves)
				(bind ?temp_spin_waves (create$))
			)
			
			(foreach ?spin_wave ?spin_waves
				(bind ?new_all_occuring_waves ?all_occuring_waves)
				(bind ?mother_index (member$ ?spin_wave ?all_occuring_waves))
				(if (not ?mother_index) then
				    (bind ?new_all_occuring_waves (insert$ ?all_occuring_waves (+ 1 (length ?all_occuring_waves)) ?spin_wave))
					(bind ?mother_index (length ?new_all_occuring_waves))
				)
				(bind ?mother_wave_unique_id (+ 1 (length ?unique_index_wave_mapping)))
				
				(bind ?decay (assert-index-decay ?mother_wave_unique_id ?wave_unique_id ?wave2_unique_id))

				;(printout t ?available_waves ?wave_unique_id ?wave2_unique_id crlf)
				(bind ?other-remaining-particles (delete$ (create$ ?available_waves) 
					(member$ ?wave_unique_id ?available_waves) (member$ ?wave_unique_id ?available_waves))
				)
				(bind ?other-remaining-particles (delete$ (create$ ?other-remaining-particles) 
					(member$ ?wave2_unique_id ?other-remaining-particles) (member$ ?wave2_unique_id ?other-remaining-particles))
				)
				(duplicate ?decay_tree 
					(decays ?decays ?decay)
					(available_waves ?other-remaining-particles ?mother_wave_unique_id)
					(all_occuring_waves ?new_all_occuring_waves)
					(unique_index_wave_mapping ?unique_index_wave_mapping 
						(assert-index-pair ?mother_wave_unique_id ?mother_index)
					)
				)
			)
		)
	  )
	  (retract ?decay_tree)
	else 
		(if (= 2 (length ?available_waves))
		then
			;(printout t "we are about to finish the decaytree here" crlf)
		
			(bind ?wave_unique_id (nth$ 1 ?available_waves))
			(bind ?wave2_unique_id (nth$ 2 ?available_waves))
			
			(bind ?dwave1 (get-wave-from-unique-index ?wave_unique_id ?decay_tree))
			(bind ?dwave2 (get-wave-from-unique-index ?wave2_unique_id ?decay_tree))
		
			(bind ?qn_name_list (get-list-of-qn-names ?dwave1 ?dwave2))
		
			(bind ?valid_tree TRUE)	
			; go through qn name list
			(foreach ?qn_name ?qn_name_list
				(bind ?single_qn_decays (find-all-facts ((?d Decay)) 
					(= 0 (str-compare ?d:quantum_number_name ?qn_name))
				))
				(bind ?found_something FALSE)
				(foreach ?single_qn_decay ?single_qn_decays
					;(printout t (fact-slot-value ?single_qn_decay mother) " "  ?initial_state_wave crlf)
					(if (= (fact-slot-value ?single_qn_decay mother) 
						(nth$ (member$ ?qn_name (fact-slot-value ?initial_state_wave quantum_number_names))
							(fact-slot-value ?initial_state_wave quantum_number_values)
						)
						)
					then
						(if (is-decay-valid ?single_qn_decay ?dwave1 ?dwave2)
						then
							(bind ?found_something TRUE)
							(break)
						)
					)
				)
				(if (not ?found_something) 
				then 
					(bind ?valid_tree FALSE)
					(break)
				)
			)
			
			(if ?valid_tree	then
				(bind ?new_all_occuring_waves ?all_occuring_waves)
				(bind ?mother_index (member$ ?initial_state_wave ?all_occuring_waves))
				(if (not ?mother_index) then
				    (bind ?new_all_occuring_waves (insert$ ?all_occuring_waves (+ 1 (length ?all_occuring_waves)) ?initial_state_wave))
					(bind ?mother_index (length ?new_all_occuring_waves))
				)
				(bind ?mother_wave_unique_id (+ 1 (length ?unique_index_wave_mapping)))
				
				(bind ?decay (assert-index-decay ?mother_wave_unique_id ?wave_unique_id ?wave2_unique_id))

				;(printout t ?all_occuring_waves ?new_all_occuring_waves crlf)

				(modify ?decay_tree 
					(decays ?decays ?decay)
					(available_waves ?mother_wave_unique_id)
					(all_occuring_waves ?new_all_occuring_waves)
					(unique_index_wave_mapping ?unique_index_wave_mapping 
						(assert-index-pair ?mother_wave_unique_id ?mother_index)
					)
				)
			else
				(retract ?decay_tree)
			)
		)
    )
	(retract ?dummy_wave)
)
