function ESN_global_variables(flag_pair_list)
%% Global variables
global tag_name_list data_type_eye_list data_type_BEHAVE_list data_type_neuro_list data_type_EPHYS_list event_type_list ...
    waveform_inds_span length_trace inds_span ...
    ang_step ang_edges ang_values amp_edges vel_edges ...
    range_cell_with_4dir_behave expand_index ...
    corr_data_list
tag_name_list = { ...
    'prim_success', ... % tag 1
    'prim_attempt', ... % tag 2
    'prim_fail', ... % tag 3
    'corr_success', ... % tag 4
    'corr_fail', ... % tag 5
    'back_center_success', ... % tag 6
    'back_center_prim', ... % tag 7
    'back_center_irrelev', ... % tag 8
    'target_irrelev', ... % tag 9
    'other_irrelev', ... % tag 10
    };
data_type_eye_list    = {'eye_vx', 'eye_vy'};
data_type_BEHAVE_list = {'eye_r_vx_filt', 'eye_r_vy_filt'};
data_type_neuro_list  = {'neuro_SS', 'neuro_CS'};
data_type_EPHYS_list  = {'EPHYS_SS_train_1K', 'EPHYS_CS_train_1K'};
event_type_list       = {'visual', 'onset', 'vmax', 'offset', 'auditory'};
waveform_inds_span = ((-60+1) : 1 : (120));
length_trace = 500;
inds_span    = ((-(length_trace/2)+1) : 1 : (length_trace/2))';
ang_step     = 45;
ang_edges    = (0 - (ang_step/2)) : ang_step : (360 + (ang_step/2));
ang_values   = (0) : ang_step : (360 - ang_step);
amp_edges    = [0 1.5 2.5 3.5 4.5 5.5 7.5 100];
vel_edges    = [0 150 250 350 450 550 650 750 10000];
expand_index = 1;
corr_data_list = {'SS1', 'SS2', 'CS1', 'CS2'};

%% if there is no flag_pair_list then return
if ~exist('flag_pair_list','var')
    return;
end
if ~flag_pair_list
    range_cell_with_4dir_behave = [1 51]; % This is to correct the data for the first round of recordings from Mirza, pCell_list_Mirza_pre201906
    % plese set the range_cell_with_4dir_behave to [-1 -1] if you do not have 4dir sessions
else
    range_cell_with_4dir_behave = [1 32]; % This is to correct the data for the first round of recordings from Mirza, pair_list_full_Mirza_pre201906
    % plese set the range_cell_with_4dir_behave to [-1 -1] if you do not have 4dir sessions
end
end
