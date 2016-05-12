;;;*************
;;;* TEMPLATES *
;;;*************

(deftemplate IndexPair
	(slot unique_id (type INTEGER))
	(slot list_index (type INTEGER))
)

(deftemplate List
	(multislot values)
)

(deftemplate NameList
	(slot name (type STRING))
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

(deftemplate AllowedDecayQuantumNumbers
	(slot name (type STRING))
	(multislot values)
	(multislot required_variable_names)
)

; allowed quantum number list
(deftemplate AllowedQuantumNumbers
	(slot name (type STRING))
	(multislot values)
	(multislot required_variable_names)
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
	(multislot required_variable_names)
	(multislot required_variables)
	(multislot required_decays)
	(multislot violating_quantum_number_list)
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

(deffunction int> (?a ?b)
	(> ?a ?b)
)

(deffunction check-multifields-content-identical (?list1 ?list2)
	(and 
		(subsetp ?list1 ?list2) 
		(subsetp ?list2 ?list1)
	)
)

(deffunction check-multifields-identical (?list1 ?list2)
    (if (<> (length ?list1) (length ?list2)) then (return FALSE))
    
	(loop-for-count (?index 1 (length ?list1))
	do
	  (bind ?val1 (nth$ ?index ?list1))
	  (bind ?val2 (nth$ ?index ?list2))
	  (if (stringp ?val1)
	  then
	  	(if (not (stringp ?val2)) then (return FALSE))
	  	(if (<> 0 (str-compare ?val1 ?val2)) then (return FALSE))
	  else
	    (if (stringp ?val2) then (return FALSE))
	    (if (<> ?val1 ?val2) then (return FALSE))
	  )
	)
	(return TRUE)
)

(deffunction get-spin-qn-with-unique-id (?uid)
	(nth$ 1 (find-fact ((?f SpinQuantumNumber)) (= ?f:unique_id ?uid)))
)

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

(deffunction get-unique-available-waves-remainder
	(?reduced_available_waves ?decay_tree)
	
	(bind ?unique_index_wave_mapping (fact-slot-value ?decay_tree unique_index_wave_mapping))
	(bind ?all_occuring_waves (fact-slot-value ?decay_tree all_occuring_waves))
	
	(bind ?unique_available_waves (create$))
	(bind ?unique_available_waves_indices_remainder (create$))
	(foreach ?wave_id ?reduced_available_waves
		(foreach ?index_pair ?unique_index_wave_mapping
			(if (= ?wave_id (fact-slot-value ?index_pair unique_id))
			then
				(bind ?wave (nth$ (fact-slot-value ?index_pair list_index) ?all_occuring_waves))
				(if (not (member$ ?wave ?unique_available_waves))
				then
					(bind ?unique_available_waves (insert$ ?unique_available_waves 1 ?wave))
				else
					(bind ?unique_available_waves_indices_remainder (insert$ ?unique_available_waves_indices_remainder 1 ?wave_id))
				)
				(break)
			)
		)
	)
	(return ?unique_available_waves_indices_remainder)
)

(deffunction create-list (?values)
	(bind ?return_value (create$ (assert (List (values ?values)))))
	(if (not (nth$ 1 ?return_value))
	then 
		(bind ?return_value (find-fact ((?l List)) (check-multifields-content-identical ?l:values ?values)))
	)
	(nth$ 1 ?return_value)
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
	(bind ?result TRUE)

	(foreach ?index_pair ?unique_index_wave_mapping
		(if (= (fact-slot-value ?index_pair unique_id) ?unique_decay_index) then
			(bind ?result (nth$ (fact-slot-value ?index_pair list_index) ?all_occuring_waves))
			(break)
		)
	)
	
	(return ?result)
)
	
(deffunction find-decay-fact
	(?mother_wave ?daughter_wave1 ?daughter_wave2 ?violating_decay_list ?req_var_name ?req_var_values)
	(bind ?return_result TRUE)
	(bind ?results 
		(find-all-facts ((?sw Decay)) (= 0 (str-compare ?sw:quantum_number_name "spinwave")))
	)
	(foreach ?result ?results
		;(printout t ?result " " (fact-slot-value ?result violating_quantum_number_list) " =? " (fact-slot-value ?violating_decay_list values) crlf)
		(if (and 
				(and 
					(= (fact-slot-value ?result mother) ?mother_wave) 
					(and (member$ ?daughter_wave1 (fact-slot-value ?result daughters))
						(member$ ?daughter_wave2 (fact-slot-value ?result daughters))
					)
				)
				(check-multifields-identical (fact-slot-value ?result violating_quantum_number_list) (fact-slot-value ?violating_decay_list values))
			)
		then
			(if 
				(and 
					(check-multifields-identical (fact-slot-value ?result required_variable_names) ?req_var_name) 
					(check-multifields-identical (fact-slot-value ?result required_variables) ?req_var_values)
				)
			then
				(bind ?return_result ?result)
				(break)
			)
		)
	)
	(return ?return_result)
)
	
(deffunction assert-index-decay
	(?mother_wave_id ?wave_unique_id ?wave2_unique_id ?violating_decay_list ?single_decay_list)
	
	;first create list of required variables
	(bind ?req_var_names (create$))
	(bind ?req_var_values (create$))
	(foreach ?single_decay (fact-slot-value ?single_decay_list values)
		(bind ?single_decay_req_var_names (fact-slot-value ?single_decay required_variable_names))
		(bind ?single_decay_req_var_values (fact-slot-value ?single_decay required_variables))
		(loop-for-count (?index 1 (length ?single_decay_req_var_names))
		do
			(bind ?sd_req_var_name (nth$ ?index ?single_decay_req_var_names))
			(bind ?found_index (member$ sd_req_var_name ?req_var_names))
			(if ?found_index
			then
				;we could perform some check that the variable values match but that should never happen
				;per definition...
			else
				(bind ?req_var_names (insert$ ?req_var_names 1 ?sd_req_var_name))
				(bind ?req_var_values (insert$ ?req_var_names 1 (nth$ ?index ?single_decay_req_var_values)))
			)
		)
	)
	
	(bind ?decay 
		(assert 
			(Decay (quantum_number_name "spinwave") 
				(mother ?mother_wave_id) (daughters ?wave_unique_id ?wave2_unique_id)
				(violating_quantum_number_list (fact-slot-value ?violating_decay_list values))
				(required_variable_names ?req_var_names)
				(required_variables ?req_var_values)
			)
		)
	)
	(if (not ?decay)
	then
		(bind ?decay 
			(find-decay-fact 
				?mother_wave_id 
				?wave_unique_id 
				?wave2_unique_id 
				?violating_decay_list
				?req_var_names
				?req_var_values
			)
		)
	)
	(return ?decay)
)

(deffunction compare-spinwave-values
	(?quantum_number_names1 ?quantum_number_names2 ?quantum_number_values1 ?quantum_number_values2)
	
	(bind ?correct_fact TRUE)
	(foreach ?quantum_number_name ?quantum_number_names1
		(bind ?index1 (member$ ?quantum_number_name ?quantum_number_names1))
		(bind ?index2 (member$ ?quantum_number_name ?quantum_number_names2))
		(if (<> (nth$ ?index1 ?quantum_number_values1) (nth$ ?index2 ?quantum_number_values2))
		then
		  (bind ?correct_fact FALSE)
		  (break)
		)
	)
	(return ?correct_fact)
)

(deffunction find-spinwave-fact-list
	(?quantum_number_names ?quantum_number_values)
		
	(bind ?found_facts (find-all-facts ((?sw SpinWave)) 
		(check-multifields-content-identical ?sw:quantum_number_names ?quantum_number_names)
	))
	
	(bind ?result (create$))
	
	(foreach ?sw_fact ?found_facts
		(bind ?fact_qn_names (fact-slot-value ?sw_fact quantum_number_names))
		(bind ?fact_qn_values (fact-slot-value ?sw_fact quantum_number_values))
		
		(if (compare-spinwave-values ?quantum_number_names ?fact_qn_names ?quantum_number_values ?fact_qn_values)
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

(deffunction check-spinwaves-equal (?spinwave1 ?spinwave2)	
	(bind ?return_value TRUE)
	(compare-spinwave-values 
	    (fact-slot-value ?spinwave1 quantum_number_names)
	    (fact-slot-value ?spinwave2 quantum_number_names)
		(fact-slot-value ?spinwave1 quantum_number_values)
		(fact-slot-value ?spinwave2 quantum_number_values)
	)
)

(deffunction check-spinwave-equality (?spinwave_index1 ?spinwave_index2 ?all_occuring_waves)
	(bind ?spinwave1 (nth$ ?spinwave_index1 ?all_occuring_waves))
	(bind ?spinwave2 (nth$ ?spinwave_index2 ?all_occuring_waves))
	
	(check-spinwaves-equal ?spinwave1 ?spinwave2)
)

(deffunction get-list-of-qn-names
	(?spin_wave1 ?spin_wave2)
	
	(bind ?all_qn_names (fact-slot-value (nth$ 1 (find-fact ((?f NameList)) (= 0 (str-compare ?f:name "AllQuantumNumbers")))) names))
	(bind ?names (create$))
	(if (and (fact-slot-value ?spin_wave1 quantum_number_names) (fact-slot-value ?spin_wave2 quantum_number_names))
	then
		(foreach ?name (fact-slot-value ?spin_wave1 quantum_number_names)
			;(printout t (fact-slot-value ?spin_wave1 quantum_number_names) " " (fact-slot-value ?spin_wave2 quantum_number_names) " " ?name crlf)
			(if (member$ ?name (fact-slot-value ?spin_wave2 quantum_number_names))
			then
				(if (member$ ?name ?all_qn_names) then
					(bind ?names (insert$ ?names (+ (length ?names) 1) ?name))
				)
			)
		)
	)
	(subseq$ ?names 1 (length ?names))
) 

(deffunction is-decay-valid 
	(?single_qn_decay ?spin_wave1 ?spin_wave2)
	;(printout t (nth$ 1 (fact-slot-value ?single_qn_decay daughters))  " ==? " (nth$ 
	;			(member$ (fact-slot-value ?single_qn_decay quantum_number_name) 
	;				(fact-slot-value ?spin_wave1 quantum_number_names)
	;			)
	;			(fact-slot-value ?spin_wave1 quantum_number_values)
	;		)
	;		crlf)
	;(printout t (nth$ 2 (fact-slot-value ?single_qn_decay daughters)) " ==? " (nth$ 
	;			(member$ (fact-slot-value ?single_qn_decay quantum_number_name) 
	;				(fact-slot-value ?spin_wave2 quantum_number_names)
	;			)
	;			(fact-slot-value ?spin_wave2 quantum_number_values)
	;		)
	;		crlf)
			
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

(deffunction get-required-variable (?variable_name ?decay)
	(bind ?index (member$ ?variable_name (fact-slot-value ?decay required_variable_names)))
	(nth$ ?index (fact-slot-value ?decay required_variables))
)

(deffunction check-decay-requirements (?single_qn_decay_list ?single_qn_decay)
	; this function checks if all the required information is there
	
	;(printout t (fact-slot-value ?single_qn_decay quantum_number_name) crlf)
	;(printout t (fact-slot-value ?single_qn_decay mother) (fact-slot-value ?single_qn_decay daughters) crlf)
	
	(bind ?required_decays (fact-slot-value ?single_qn_decay required_decays))
	(foreach ?required_decay ?required_decays
		(bind ?is_decay_in_list FALSE)
		(foreach ?single_available_decay (fact-slot-value ?single_qn_decay_list values)
		do
			(if (eq ?single_available_decay ?required_decay)
			then
				(bind ?is_decay_in_list TRUE)
				(break)
			)
		)
		(if (not ?is_decay_in_list)
		then
			(return FALSE)
		)
	)
	
	;now check for conflicts of the required variables with other decays
	(bind ?required_var_names (fact-slot-value ?single_qn_decay required_variable_names))
	(bind ?required_var_values (fact-slot-value ?single_qn_decay required_variables))
	(loop-for-count (?i 1 (length ?required_var_names))
	do
		;(printout t (nth$ ?i ?required_var_names) crlf)
		(foreach ?single_available_decay (fact-slot-value ?single_qn_decay_list values)
			;find required variable name
			(bind ?found_index
				(member$
					(nth$ ?i ?required_var_names) 
					(fact-slot-value ?single_available_decay required_variable_names)
				)
			)
			
			;(bind ?angular_momentum (get-spin-qn-with-unique-id (get-required-variable "angular-momentum" ?single_qn_decay)))
			;(bind ?L (/ (fact-slot-value ?angular_momentum numerator) (fact-slot-value ?angular_momentum denominator)))
			;(printout t ?L " " (fact-slot-value ?angular_momentum z_component_numerator) crlf)
			;(printout t ?found_index " in " (fact-slot-value ?single_available_decay quantum_number_name) crlf)
			(if ?found_index
			then
			
				;(bind ?angular_momentum (get-spin-qn-with-unique-id (get-required-variable "angular-momentum" ?single_available_decay)))
				;(bind ?L (/ (fact-slot-value ?angular_momentum numerator) (fact-slot-value ?angular_momentum denominator)))
				;(printout t "other " ?L crlf)
			
				;(printout t (get-required-variable (nth$ ?i ?required_var_names) ?single_available_decay) " ==? " (nth$ ?i ?required_var_values) crlf)
				(if (<> (get-required-variable (nth$ ?i ?required_var_names) ?single_available_decay) (nth$ ?i ?required_var_values))
				then
					(return FALSE)
				)
			)
		)
	)
	
	(return TRUE)
)

(deffunction get-qn-value (?index ?decay_tree ?qn_name)
	(bind ?spinwave_fact (get-wave-from-unique-index ?index ?decay_tree))
	
	(bind ?val_index (member$ ?qn_name (fact-slot-value ?spinwave_fact quantum_number_names)))
	(nth$ ?val_index (fact-slot-value ?spinwave_fact quantum_number_values))
)

(deffunction is-boson (?index ?decay_tree)
	(bind ?spin_unique_id (get-qn-value ?index ?decay_tree "spin"))
	(bind ?spin_qn (nth$ 1 (find-fact ((?sqn SpinQuantumNumber)) (= ?sqn:unique_id ?spin_unique_id))))
	
	(= 0 (mod (fact-slot-value ?spin_qn numerator) (fact-slot-value ?spin_qn denominator)))
)

(deffunction is-fermion (?index ?decay_tree)
	(bind ?spin_unique_id (get-qn-value ?index ?decay_tree "spin"))
	(bind ?spin_qn (nth$ 1 (find-fact ((?sqn SpinQuantumNumber)) (= ?sqn:unique_id ?spin_unique_id))))
	
	(<> 0 (mod (fact-slot-value ?spin_qn numerator) (fact-slot-value ?spin_qn denominator)))
)

(deffunction is-qn-conserved (?qn_name)
    (bind ?conserved_qn_list 
    	(fact-slot-value 
    		(nth$ 1 (find-fact ((?f NameList)) (= 0 (str-compare ?f:name "ConservedQuantumNumbers"))))
    		names
    	)
    )
	(member$ ?qn_name ?conserved_qn_list)
)


(deffunction get-insert-index-for-correct-requirements (?qn_name)
	(bind ?quantum_number_name_list 
    		(nth$ 1 (find-fact ((?f NameList)) (= 0 (str-compare ?f:name "AllQuantumNumbers"))))
    )
    (bind ?current_names (fact-slot-value ?quantum_number_name_list names))
	
	(bind ?insert_index (+ (length ?current_names) 1))
	(loop-for-count (?i 1 (length ?current_names))
	do
		(bind ?allowd_qn_fact 
			(nth$ 1 (find-fact ((?g AllowedQuantumNumbers)) (= 0 (str-compare ?g:name (nth$ ?i ?current_names)))))
		)
		(bind ?required_variable_names (fact-slot-value ?allowd_qn_fact required_variable_names))

		(if (> (length ?required_variable_names) 0)
		then
			(if (member$ ?qn_name ?required_variable_names)
			then
				(bind ?insert_index ?i)
				(break)
			)
		)
	)
	(return ?insert_index)
)

(deffunction insert-into-qn-hierarchy (?qn_name)
	(bind ?quantum_number_name_list 
    		(nth$ 1 (find-fact ((?f NameList)) (= 0 (str-compare ?f:name "AllQuantumNumbers"))))
    )
    (bind ?current_names (fact-slot-value ?quantum_number_name_list names))
    
	(if (not (member$ ?qn_name ?current_names))
	then
		(modify ?quantum_number_name_list 
			(names (insert$ 
					?current_names 
					(get-insert-index-for-correct-requirements ?qn_name) 
					?qn_name
				   )
			)
		)
	)
)



(deffunction get-required-variable-values (?required_variable)
	(bind ?value_list (create$))
	
	(bind ?fact 
			(find-fact ((?f AllowedDecayQuantumNumbers)) (= 0 (str-compare ?f:name ?required_variable)))
	)

	(if (> (length ?fact) 0)
	then
		(bind ?value_list (fact-slot-value (nth$ 1 ?fact) values))
	)
	(create$ ?value_list)
)

(deffunction get-required-decays (?required_variable)
	(find-all-facts ((?f Decay)) (= 0 (str-compare ?f:quantum_number_name ?required_variable)))
)

(deffunction required-variables-exist (?required_variable_names)
	(bind ?return_value (create$ TRUE))
	(foreach ?required_variable ?required_variable_names
	do
		(bind ?variable_values (get-required-variable-values ?required_variable))
		(if (= (length ?variable_values) 0)
		then
			(if (= (length (get-required-decays ?required_variable)) 0)
			then
				(bind ?return_value (create$ FALSE))
				(break)
			)
		)
	)
	(nth$ 1 ?return_value)
)

(deffunction generate-decay-requirements-dummys (?required_variable_names)
	(bind ?initial_decay_dummy (assert (Decay (quantum_number_name "dummy-decay"))))
	(foreach ?required_variable_name ?required_variable_names
		(bind ?decay_dummies (find-all-facts ((?f Decay)) (= 0 (str-compare ?f:quantum_number_name "dummy-decay"))))
		(bind ?variables (get-required-variable-values ?required_variable_name))
		(foreach ?variable ?variables
			(foreach ?decay_dummy ?decay_dummies
				(duplicate ?decay_dummy 
					(required_variable_names (insert$ (fact-slot-value ?decay_dummy required_variable_names) 1 ?required_variable_name)) 
					(required_variables (insert$ (fact-slot-value ?decay_dummy required_variables) 1 ?variable))
				)
			)
		)
		
		(if (= (length ?variables) 0)
		then
			(bind ?decays (get-required-decays ?required_variable_name))
			(foreach ?decay ?decays
				(foreach ?decay_dummy ?decay_dummies
					(duplicate ?decay_dummy (required_decays (insert$ (fact-slot-value ?decay_dummy required_decays) 1 ?decay)))
				)
			)
		)
		
		;delete old decay dummies
		(foreach ?decay_dummy ?decay_dummies
			(retract ?decay_dummy)
		)
	)
	(retract ?initial_decay_dummy)
)

(deffunction duplicate-decays-on-requirements (?decay ?required_variable_names)
	;(printout t "duplicating..." crlf)
	(generate-decay-requirements-dummys ?required_variable_names)
	(bind ?decay_dummies (find-all-facts ((?f Decay)) (= 0 (str-compare ?f:quantum_number_name "dummy-decay"))))
	(foreach ?decay_dummy ?decay_dummies
	do
		(duplicate ?decay
			(required_variable_names (fact-slot-value ?decay_dummy required_variable_names))
			(required_variables (fact-slot-value ?decay_dummy required_variables))
			(required_decays (fact-slot-value ?decay_dummy required_decays))
		)
	)
	(if (> (length ?decay_dummies) 0)
	then
		(retract ?decay)
		(foreach ?decay_dummy ?decay_dummies
		do
			(retract ?decay_dummy)
		)
	)
)

;;;*******************
;;;* USER CONDITIONS *
;;;*******************

(deffacts initial-facts-setup
	(NameList (name "AllQuantumNumbers"))
	;(NameList (name "ConservedQuantumNumbers"))
)

