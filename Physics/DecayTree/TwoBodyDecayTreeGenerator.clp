;;;*******************************
;;;* DECAY TREE GENERATION RULES *
;;;*******************************


; generate all possible allowed decays
(deffunction create-decays ()
	(printout t "creating decays..." crlf)
	
	(bind ?aqn_facts (find-all-facts ((?f AllowedQuantumNumbers)) TRUE))
	(bind ?list_shrunk TRUE)
	
	(while ?list_shrunk
	do
		(bind ?list_shrunk FALSE)
		(bind ?temp_remaining_aqn_facts ?aqn_facts)
		(foreach ?aqn_fact ?aqn_facts
		do
			(bind ?quantum_number_name (fact-slot-value ?aqn_fact name))
			(bind ?values (fact-slot-value ?aqn_fact values))
			(bind ?required_variable_names (fact-slot-value ?aqn_fact required_variable_names))
			(if (required-variables-exist ?required_variable_names)
			then
				(foreach ?mother_value ?values
					(foreach ?daughter1_value ?values
						(foreach ?daughter2_value ?values
							(bind ?decay 
								(assert
									(Decay 
										(quantum_number_name ?quantum_number_name) 
										(mother ?mother_value)
										(daughters ?daughter1_value ?daughter2_value)
									)
								)
							)
							(duplicate-decays-on-requirements ?decay ?required_variable_names)
						)
					)
				)
				(insert-into-qn-hierarchy ?quantum_number_name)
				(bind ?temp_remaining_aqn_facts (delete-member$ ?temp_remaining_aqn_facts ?aqn_fact))
				(bind ?list_shrunk TRUE)
			)
		)
		(bind ?aqn_facts ?temp_remaining_aqn_facts)
	)
)

; generate seed decay trees
(defrule create-seed-decay-trees
	(declare (salience 5))
	(InitialAndFinalState 
		(initial_state ?initial_state) 
		(final_state $?final_state)
	)
	=>
	(if (= (length (fact-slot-value 
					(nth$ 1 (find-fact ((?f NameList)) (= 0 (str-compare ?f:name "AllQuantumNumbers"))))
					names
				   )
			)
			0
		)
	then
		(create-decays)
	)
	
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
	(declare (salience 5))
	?decay_tree <- (DecayTree (available_waves $?available_waves) (decays $?decays)
					(initial_state_wave ?initial_state_wave) (all_occuring_waves $?all_occuring_waves)
					(unique_index_wave_mapping $?unique_index_wave_mapping))
	=>
	;(printout t "number of decay trees: " (length (find-all-facts ((?d DecayTree)) TRUE)) crlf)
	;(printout t "available waves in this tree: " (length ?available_waves) crlf)
	(bind ?dummy_wave (assert (SpinWave)))
	(bind ?dummy_list (create-list (create$)))
	(if (> (length ?available_waves) 2)
	then
	  (bind ?unique_available_waves (get-unique-available-waves ?available_waves ?decay_tree))
	  (bind ?non_unique_left_overs (get-unique-available-waves-remainder ?available_waves ?decay_tree))
	  
	  (foreach ?wave_unique_id ?unique_available_waves
	  	(bind ?wave (get-wave-from-unique-index ?wave_unique_id ?decay_tree))
	  	
	    (bind ?reduced_available_waves (insert$ (delete-member$ ?unique_available_waves ?wave_unique_id) 1 ?non_unique_left_overs))
	  	(bind ?unique_available_waves (get-unique-available-waves ?reduced_available_waves ?decay_tree))
		
		(foreach ?wave2_unique_id ?unique_available_waves
			(bind ?wave2 (get-wave-from-unique-index ?wave2_unique_id ?decay_tree))
			
			(bind ?spin_waves (create$ ?dummy_wave))
			(bind ?temp_spin_waves (create$))
			
			(bind ?single_decay_lists (create$ ?dummy_list))
			(bind ?decay_violating_lists (create$ ?dummy_list))
			(bind ?temp_single_decay_lists (create$ ?dummy_list))
			(bind ?temp_decay_violating_lists (create$ ?dummy_list))
			
			(bind ?qn_name_list (get-list-of-qn-names ?wave ?wave2))
			
			(foreach ?qn_name ?qn_name_list
				(bind ?extended_tree FALSE)
				(bind ?spin_wave_index 0)
				;(printout t "qn name " ?qn_name crlf)
				;(printout t "spin waves " ?spin_waves crlf)
				(foreach ?spin_wave ?spin_waves
					(bind ?spin_wave_index (+ 1 ?spin_wave_index))
					(bind ?single_qn_decays (find-all-facts ((?d Decay)) (= 0 (str-compare ?d:quantum_number_name ?qn_name))))
					;(printout t "qn names for spin wave " (fact-slot-value ?spin_wave quantum_number_names) crlf)
					(foreach ?single_qn_decay ?single_qn_decays
						;(printout t ?qn_name " of " ?qn_name_list crlf)
						;(printout t (fact-slot-value ?single_qn_decay quantum_number_name) crlf)
						;(printout t (fact-slot-value ?single_qn_decay mother) " " (fact-slot-value ?single_qn_decay daughters) crlf)
						(if (is-decay-valid ?single_qn_decay ?wave ?wave2)
						then
							(bind ?requirements_satisfied TRUE)
							(bind ?spin_wave_fact TRUE)
							(if (not (fact-slot-value ?spin_wave quantum_number_names))
							then
							    (bind ?spin_wave_fact
									(duplicate ?spin_wave 
										(quantum_number_names ?qn_name) 
										(quantum_number_values (fact-slot-value ?single_qn_decay mother))
									)
								)
								
								(if (not ?spin_wave_fact)
								then
									(bind ?spin_wave_fact 
										(find-spinwave-fact (create$ ?qn_name)
											(create$ (fact-slot-value ?single_qn_decay mother))
										)
									)
								)
							else			
								;check the requirements of the new decay
								(bind ?requirements_satisfied (check-decay-requirements 
																(nth$ ?spin_wave_index ?single_decay_lists) ?single_qn_decay
															  )
								)
								(if ?requirements_satisfied
								then
									(bind ?spin_wave_fact
										(duplicate ?spin_wave 
											(quantum_number_names (fact-slot-value ?spin_wave quantum_number_names) ?qn_name) 
											(quantum_number_values (fact-slot-value ?spin_wave quantum_number_values) 
												(fact-slot-value ?single_qn_decay mother)
											)
										)
									)
									;(printout t "this is good!" crlf)
									(if (not ?spin_wave_fact)
									then
										;(printout t (fact-slot-value ?spin_wave quantum_number_names) " " (fact-slot-value ?spin_wave quantum_number_values) crlf)
										(bind ?spin_wave_fact 
											(find-spinwave-fact 
												(insert$ (fact-slot-value ?spin_wave quantum_number_names) 1 ?qn_name)
												(insert$ (fact-slot-value ?spin_wave quantum_number_values) 1 (fact-slot-value ?single_qn_decay mother))
											)
										)
									)
								)
							)
							;(printout t ?qn_name " " ?requirements_satisfied crlf)
							(if ?requirements_satisfied
							then
								(bind ?temp_spin_waves 
									(insert$ ?temp_spin_waves 1 ?spin_wave_fact)
								)
								
								(bind ?temp_decay_violating_list (nth$ ?spin_wave_index ?decay_violating_lists))
							
								(bind ?temp_violated_qns_for_decay (fact-slot-value ?single_qn_decay violating_quantum_number_list))					    
								(bind ?list_values (fact-slot-value ?temp_decay_violating_list values))
								(foreach ?violated_qn ?temp_violated_qns_for_decay
									(bind ?list_values (insert$ ?list_values 1 ?violated_qn))
								)
								
								(bind ?temp_decay_violating_lists 
									(insert$ 
										?temp_decay_violating_lists 
										1 
										(create-list ?list_values)
									)
								)
								

								(bind ?temp_single_decay_list (nth$ ?spin_wave_index ?single_decay_lists))

								(bind ?temp_single_decay_lists 
									(insert$ 
										?temp_single_decay_lists 
										1
										(create-list
											(insert$ (fact-slot-value ?temp_single_decay_list values) 1 ?single_qn_decay)
										)
									)
								)
								(bind ?extended_tree TRUE)
							
							else
								;(printout t "did not meet requirements" crlf)
							)
						)
					)
				)
				;(foreach ?spin_wave ?spin_waves
				;do
				;	(retract ?spin_wave)
				;)
				
				;(if ?extended_tree then
				  (bind ?spin_waves ?temp_spin_waves)
				  (bind ?single_decay_lists ?temp_single_decay_lists)
				  (bind ?decay_violating_lists ?temp_decay_violating_lists)
				  (bind ?temp_spin_waves (create$))
				  (bind ?temp_single_decay_lists (create$ ?dummy_list))
				  (bind ?temp_decay_violating_lists (create$ ?dummy_list))
				;)
			)
			
			;(printout t "we have that many solutions: " (length ?spin_waves) crlf)
			(bind ?spin_wave_index 0)
			(foreach ?spin_wave ?spin_waves
				;(printout t (fact-slot-value ?spin_wave quantum_number_names) crlf)
				
				;spawn new qn here (like c-parity)
				
				(bind ?spin_wave_index (+ 1 ?spin_wave_index))
				(bind ?new_all_occuring_waves ?all_occuring_waves)
				(bind ?mother_index (member$ ?spin_wave ?all_occuring_waves))
				(if (not ?mother_index) then
				    (bind ?new_all_occuring_waves (insert$ ?all_occuring_waves (+ 1 (length ?all_occuring_waves)) ?spin_wave))
					(bind ?mother_index (length ?new_all_occuring_waves))
				)
				(bind ?mother_wave_unique_id (+ 1 (length ?unique_index_wave_mapping)))
				
				;(printout t "adding decay to tree " (nth$ ?spin_wave_index ?decay_violating_lists) crlf)
				(bind ?decay 
					(assert-index-decay 
						?mother_wave_unique_id 
						?wave_unique_id 
						?wave2_unique_id 
						(nth$ ?spin_wave_index ?decay_violating_lists)
						(nth$ ?spin_wave_index ?single_decay_lists)
					)
				)

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
			;(printout t "number of decay trees now...: " (length (find-all-facts ((?d DecayTree)) TRUE)) crlf)
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
			;(printout t ?qn_name_list crlf)
			(bind ?single_decay_list ?dummy_list)
			(bind ?decay_violating_list ?dummy_list)
		
			(bind ?valid_tree TRUE)
			; go through qn name list
			(foreach ?qn_name ?qn_name_list
				(bind ?single_qn_decays (find-all-facts ((?d Decay))
					(= 0 (str-compare ?d:quantum_number_name ?qn_name))
				))
				(printout t ?qn_name crlf)
				(bind ?found_something FALSE)
				(foreach ?single_qn_decay ?single_qn_decays
				    (printout t (fact-slot-value ?single_qn_decay mother) " " (nth$ (member$ ?qn_name (fact-slot-value ?initial_state_wave quantum_number_names))
							(fact-slot-value ?initial_state_wave quantum_number_values)
						) crlf)
						
					(if (= (fact-slot-value ?single_qn_decay mother)
						(nth$ (member$ ?qn_name (fact-slot-value ?initial_state_wave quantum_number_names))
							(fact-slot-value ?initial_state_wave quantum_number_values)
						)
						)
					then
						(if (is-decay-valid ?single_qn_decay ?dwave1 ?dwave2)
						then
							;check the requirements of the new decay
							(printout t "checking if this decay is good " ?single_qn_decay crlf)
							(if (check-decay-requirements ?single_decay_list ?single_qn_decay)
							then
								(bind ?single_decay_list (create-list (insert$ (fact-slot-value ?single_decay_list values) 1 ?single_qn_decay)))
								(bind ?found_something TRUE)
							
							
								(bind ?temp_violated_qns_for_decay (fact-slot-value ?single_qn_decay violating_quantum_number_list))
								(foreach ?violated_qn ?temp_violated_qns_for_decay
									(bind ?decay_violating_list
										(create-list (insert$ (fact-slot-value ?decay_violating_list values) 1 ?violated_qn))
									)
								)
								
								(printout t "this should never happend then..." crlf)
								
								(break)
							)
						)
					)
				)
				(if (not ?found_something)
				then 
				    (printout t "found nothing for " ?qn_name)
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
				
				(bind ?decay 
					(assert-index-decay 
						?mother_wave_unique_id 
						?wave_unique_id 
						?wave2_unique_id 
						?decay_violating_list 
						?single_decay_list
					)
				)

				;(printout t "wuuttt" crlf)
				
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
    (retract ?dummy_list)
	(retract ?dummy_wave)
)
