label2tissue = {}

label2tissue[0] = 'background'
label2tissue[1] = 'body'
label2tissue[2] = 'skin'
label2tissue[3] = 'thyroid'
label2tissue[9] = 'lesn'

#HEART
label2tissue[12] = 'hrt_myoLV'
label2tissue[13] = 'hrt_myoRV'
label2tissue[14] = 'hrt_myoLA'
label2tissue[15] = 'hrt_myoRA'
label2tissue[16] = 'hrt_bldplLV'
label2tissue[17] = 'hrt_bldplRV'
label2tissue[18] = 'hrt_bldplLA'
label2tissue[19] = 'hrt_bldplRA'

label2tissue[22] = 'lung'
label2tissue[23] = 'airway'
label2tissue[24] = 'bladder'
label2tissue[25] = 'vas_def'
label2tissue[26] = 'testicular'

label2tissue[29] = 'liver'
label2tissue[30] = 'gall_bladder'
label2tissue[31] = 'st_wall'  #stomac
label2tissue[32] = 'st_cnts'
label2tissue[33] = 'pancreas'
label2tissue[34] = 'kidney'

label2tissue[37] = 'spleen'
label2tissue[38] = 'sm_intest'
label2tissue[39] = 'large_intest'

label2tissue[41] = 'li_air'
label2tissue[42] = 'si_air'

# BRAIN
label2tissue[43] = 'brain'
label2tissue[44] = 'cerebral_cortex'
label2tissue[45] = 'cerebellum'
label2tissue[46] = 'corpus_callosum'
label2tissue[47] = 'brainstem'

label2tissue[49] = 'striatum'
label2tissue[50] = 'thal'
label2tissue[51] = 'hippo'
label2tissue[52] = 'hypothalamus'
label2tissue[53] = 'amygdala'
label2tissue[54] = 'lateral_septal_nuclei'
label2tissue[55] = 'anterior_commissure'
label2tissue[56] = 'anterior_pretectal_nucleus'
label2tissue[57] = 'periaqueductal_gray'
label2tissue[58] = 'aqueduct'
label2tissue[59] = 'cerebral_peduncle'
label2tissue[60] = 'cochlear_nuclei'
label2tissue[61] = 'deep_mesencephalic_nuclei'
label2tissue[62] = 'fimbria'
label2tissue[63] = 'fornix'
label2tissue[64] = 'globus_pallidus'
label2tissue[65] = 'inferior_colliculus'
label2tissue[66] = 'internal_capsule'
label2tissue[67] = 'interpeduncular_nucleus'
label2tissue[68] = 'lateral_dorsal_nucleus_of_thalamus_actucleus_of_thalamus'
label2tissue[69] = 'lateral_geniculate'
label2tissue[70] = 'lateral_lemniscus'
label2tissue[71] = 'medial_geniculate'
label2tissue[72] = 'nucleus_accumbens'
label2tissue[73] = 'olfactory_areas'
label2tissue[74] = 'optic_tract'
label2tissue[75] = 'pontine_gray'
label2tissue[76] = 'spinal_trigeminal_tract'
label2tissue[77] = 'substantia_nigra'
label2tissue[78] = 'superior_colliculus'
label2tissue[79] = 'pineal_gland'
label2tissue[80] = 'ventral_thalamic_nuclei'
label2tissue[81] = 'ventricular_system'

# BONES
label2tissue[82] = 'rib'
label2tissue[83] = 'skull'

label2tissue[85] = 'humerus'
label2tissue[86] = 'radius'
label2tissue[87] = 'ulna'
label2tissue[88] = 'femur'
label2tissue[89] = 'fibula'
label2tissue[90] = 'tibia'
label2tissue[91] = 'patella'
label2tissue[92] = 'bone'
label2tissue[93] = 'marrow'
label2tissue[94] = 'spine'

label2tissue[96] = 'vein'
label2tissue[97] = 'artery'
label2tissue[98] = 'tumor'

tissue2label = {v: k for k,v in label2tissue.items() }

merge_labels = {}
merge_labels["intestin"] = ([38, 39, 41, 42], 38)
merge_labels["brain"]    = ([v for v in  range(43,48)] + [v for v in range(49, 82)], 43)
merge_labels["bones"]    = ([82,83] + [v for v in range(85,95)], 92)