;;;*************
;;;* TEMPLATES *
;;;*************

(deftemplate IndexPair
	(slot unique_id (type INTEGER))
	(slot list_index (type INTEGER))
)

(deftemplate QuantumNumberNameList
	(multislot names)
)

(deftemplate SpinQuantumNumber
	(slot unique_id (type INTEGER))
	(slot numerator (type INTEGER))
	(slot denominator (type INTEGER))
	(slot z_component_numerator (type INTEGER))
)

; define spin wave	
(deftemplate SpinWave
    (multislot quantum_number_names)
	(multislot quantum_number_values)
)

; allowed quantum number list
(deftemplate AllowedQuantumNumbers
	(slot name (type STRING))
	(multislot values)
)

; initial and final state of a single IF instance
(deftemplate InitialAndFinalState
	(slot initial_state)
	(multislot final_state)
)
	
; define decay template
(deftemplate Decay
	(slot quantum_number_name (type STRING))
	(slot mother)
	(multislot daughters)
	(slot angular_momentum)
)

; define two body decay sequence or tree
(deftemplate DecayTree
	(slot initial_state_wave)
	(multislot all_occuring_waves)
	(multislot available_waves)
	(multislot unique_index_wave_mapping)
	(multislot decays)
)

;;;*************
;;;* FUNCTIONS *
;;;*************

(deffunction get-unique-available-waves
	(?reduced_available_waves ?decay_tree)
	
	(bind ?unique_index_wave_mapping (fact-slot-value ?decay_tree unique_index_wave_mapping))
	(bind ?all_occuring_waves (fact-slot-value ?decay_tree all_occuring_waves))
	
	(bind ?unique_available_waves (create$))
	(bind ?unique_available_waves_indices (create$))
	(foreach ?wave_id ?reduced_available_waves
		(foreach ?index_pair ?unique_index_wave_mapping
			(if (= ?wave_id (fact-slot-value ?index_pair unique_id))
			then
				(bind ?wave (nth$ (fact-slot-value ?index_pair list_index) ?all_occuring_waves))
				(if (not (member$ ?wave ?unique_available_waves))
				then
					(bind ?unique_available_waves (insert$ ?unique_available_waves 1 ?wave))
					(bind ?unique_available_waves_indices (insert$ ?unique_available_waves_indices 1 ?wave_id))
				)
				(break)
			)
		)
	)
	(subseq$ ?unique_available_waves_indices 1 (length ?unique_available_waves_indices))
)

(deffunction assert-index-pair
	(?unique_decay_index ?spin_wave_list_index)
	(bind ?returnvalue (create$))
	(bind ?result 
		(find-fact ((?ip IndexPair)) 
			(and (= ?ip:unique_id ?unique_decay_index) (= ?ip:list_index ?spin_wave_list_index))
		)
	)
	(if (= 0 (length ?result)) then
		(bind ?returnvalue 
			(insert$ (create$) 1 
				(assert (IndexPair (unique_id ?unique_decay_index) (list_index ?spin_wave_list_index)))
			)
		)
	else
		(bind ?returnvalue 
			(insert$ (create$) 1 ?result)
		)
	)
	(nth$ 1 ?returnvalue)
)

(deffunction get-wave-from-unique-index
	(?unique_decay_index ?decay_tree)
	
	(bind ?unique_index_wave_mapping (fact-slot-value ?decay_tree unique_index_wave_mapping))
	(bind ?all_occuring_waves (fact-slot-value ?decay_tree all_occuring_waves))
	(bind ?result (create$))

	(foreach ?index_pair ?unique_index_wave_mapping
		(if (= (fact-slot-value ?index_pair unique_id) ?unique_decay_index) then
			(bind ?result 
				(insert$ (create$) 1 
					(nth$ (fact-slot-value ?index_pair list_index) ?all_occuring_waves)
				)
			)
			(break)
		)
	)
	
	(nth$ 1 ?result)
)
	
(deffunction find-decay-fact
	(?mother_wave ?daughter_wave1 ?daughter_wave2)
	(bind ?return_results (create$))
	(bind ?results 
		(find-all-facts ((?sw Decay)) (= 0 (str-compare ?sw:quantum_number_name "spinwave")))
	)
	(foreach ?result ?results
		;(printout t ?result " " (fact-slot-value ?result mother) " =? " ?mother_wave crlf)
		(if (and 
				(= (fact-slot-value ?result mother) ?mother_wave) 
				(and (member$ ?daughter_wave1 (fact-slot-value ?result daughters))
					(member$ ?daughter_wave2 (fact-slot-value ?result daughters))
				)
			)
		then
			(bind ?return_results (insert$ ?return_results 1 ?result))
			(break)
		)
	)
	(nth$ 1 ?return_results)
)
	
(deffunction assert-index-decay
	(?mother_wave_id ?wave_unique_id ?wave2_unique_id)
	
	(bind ?decay 
		(assert 
			(Decay (quantum_number_name "spinwave") 
				(mother ?mother_wave_id) (daughters ?wave_unique_id ?wave2_unique_id)
			)
		)
	)
	(if (not ?decay)
	then
		(bind ?decay (find-decay-fact ?mother_wave_id ?wave_unique_id ?wave2_unique_id))
	)
	(nth$ 1 (insert$ (create$) 1 ?decay))
)

(deffunction find-spinwave-fact-list
	(?quantum_number_names ?quantum_number_values)
		
	(bind ?found_facts (find-all-facts ((?sw SpinWave)) 
		(and (subsetp ?sw:quantum_number_names ?quantum_number_names) 
			(subsetp ?quantum_number_names ?sw:quantum_number_names)
		)
	))
	
	(bind ?result (create$))
	
	(foreach ?sw_fact ?found_facts
		(bind ?fact_qn_names (fact-slot-value ?sw_fact quantum_number_names))
		(bind ?fact_qn_values (fact-slot-value ?sw_fact quantum_number_values))
		(bind ?correct_fact TRUE)
		(foreach ?quantum_number_name ?quantum_number_names
			(bind ?fact_qn_index (member$ ?quantum_number_name ?fact_qn_names))
			(bind ?index (member$ ?quantum_number_name ?quantum_number_names))
			(if (<> (nth$ ?fact_qn_index ?fact_qn_values) (nth$ ?index ?quantum_number_values))
			then
			  (bind ?correct_fact FALSE)
			  (break)
			)
		)
		(if ?correct_fact
		then
			(bind ?result (insert$ ?result 1 ?sw_fact))
		    (break)
		)
	)
	
	(create$ ?result)
)

(deffunction find-spinwave-fact
	(?quantum_number_names ?quantum_number_values)
	
	(nth$ 1 (find-spinwave-fact-list ?quantum_number_names ?quantum_number_values))
)

(deffunction get-list-of-qn-names
	(?spin_wave1 ?spin_wave2)
	(bind ?names (create$))
	(if (and (fact-slot-value ?spin_wave1 quantum_number_names) (fact-slot-value ?spin_wave2 quantum_number_names))
	then
		(foreach ?name (fact-slot-value ?spin_wave1 quantum_number_names)
			;(printout t (fact-slot-value ?spin_wave1 quantum_number_names) " " (fact-slot-value ?spin_wave2 quantum_number_names) " " ?name crlf)
			(if (member$ ?name (fact-slot-value ?spin_wave2 quantum_number_names))
			then
				(bind ?names (insert$ ?names 1 ?name))
			)
		)
	)
	(subseq$ ?names 1 (length ?names))
) 

(deffunction is-decay-valid 
	(?single_qn_decay ?spin_wave1 ?spin_wave2)
	(and 
		(= (nth$ 1 (fact-slot-value ?single_qn_decay daughters))
			(nth$ 
				(member$ (fact-slot-value ?single_qn_decay quantum_number_name) 
					(fact-slot-value ?spin_wave1 quantum_number_names)
				)
				(fact-slot-value ?spin_wave1 quantum_number_values)
			)
		)
		(= (nth$ 2 (fact-slot-value ?single_qn_decay daughters))
			(nth$ 
				(member$ (fact-slot-value ?single_qn_decay quantum_number_name) 
					(fact-slot-value ?spin_wave2 quantum_number_names)
				)
				(fact-slot-value ?spin_wave2 quantum_number_values)
			)
		)
	)
)

;;;*******************
;;;* USER CONDITIONS *
;;;*******************

(deffacts initial-facts-setup
	(QuantumNumberNameList)
)  
