%% function buildLickData
function lick_data_dir = PGH_buildLickData(LICKS_ALL_DATA,tag_id,alignment,params,funcs)

length_trace = params.lick.length_trace;
ang_edges    = params.lick.ang_edges;
ang_values   = params.lick.ang_values;


% Build lick_data
idx_licks    = ismember(LICKS_ALL_DATA.tag, tag_id);
idx_licks    = idx_licks & LICKS_ALL_DATA.validity;

lick_data.time    = LICKS_ALL_DATA.(['time_',alignment])(:,idx_licks);
lick_data.SS      = LICKS_ALL_DATA.(['neuro_SS_',alignment])(:,idx_licks);
lick_data.CS      = LICKS_ALL_DATA.(['neuro_CS_',alignment])(:,idx_licks);
lick_data.tongue_dm  = LICKS_ALL_DATA.(['tongue_dm_',alignment])(:,idx_licks);
lick_data.tongue_vm  = LICKS_ALL_DATA.(['tongue_vm_',alignment])(:,idx_licks);
lick_data.tongue_ang = LICKS_ALL_DATA.(['tongue_ang_',alignment])(:,idx_licks);

lick_data.tongue_dm_max =  LICKS_ALL_DATA.tongue_dm_max(    :,idx_licks);
lick_data.tongue_vm_max =  LICKS_ALL_DATA.tongue_vm_max(    :,idx_licks);
lick_data.tongue_vm_min =  LICKS_ALL_DATA.tongue_vm_min(    :,idx_licks);
lick_data.tongue_ang_max =  LICKS_ALL_DATA.tongue_ang_max(    :,idx_licks);

lick_data.time_onset   = LICKS_ALL_DATA.time_onset (     :,idx_licks);
lick_data.time_vmax   = LICKS_ALL_DATA.time_vmax (     :,idx_licks);
lick_data.time_dmax  = LICKS_ALL_DATA.time_dmax   (     :,idx_licks);
lick_data.time_vmin   = LICKS_ALL_DATA.time_vmin (     :,idx_licks);
lick_data.time_offset  = LICKS_ALL_DATA.time_offset(     :,idx_licks);

lick_data.time_diff_onset_offset = lick_data.time_offset - lick_data.time_onset;
lick_data.time_diff_onset_dmax= lick_data.time_dmax - lick_data.time_onset;
lick_data.time_diff_offset_dmax= lick_data.time_offset - lick_data.time_onset;

% Build lick_data_dir
lick_data.tongue_ang_bin = 6 - discretize(lick_data.tongue_ang_max, ang_edges);

% 5: -90deg % 4: -45deg % 3: 0deg % 2: 45deg % 1: 90deg

if length(ang_values) ~= 5
    error('plot_single_session_modulation: length ang_values is not 5. Please modify the code.')
end
if length_trace ~= 2000
    error('lick_modulation_index: length_trace is not 2000. Please modify the code.')
end
field_names_lick_data = fieldnames(lick_data);
lick_data_dir = struct;
for counter_dir = 1 : length(ang_values)
    for counter_field = 1 : length(field_names_lick_data)
        field_name = field_names_lick_data{counter_field};
        idx_ang = (lick_data.tongue_ang_bin == counter_dir);
        lick_data_dir(counter_dir).(field_name) = lick_data.(field_name)(:,idx_ang);
    end
end
lick_data_dir = rmfield(lick_data_dir,'tongue_ang_bin');
end
