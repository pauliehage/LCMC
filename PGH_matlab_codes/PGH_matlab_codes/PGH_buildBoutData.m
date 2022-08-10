%% function buildBoutData
function bout_data_dir = PGH_buildBoutData(LICKS_ALL_DATA, alignment, params, funcs)

tag_id = 1:3;

if strcmp(alignment,'offset')
    tag_id = tag_id + 3;
end

% Build lick_data
idx_bout    = ismember(LICKS_ALL_DATA.tag_bout, tag_id);
% idx_bout    = idx_bout & LICKS_ALL_DATA.validity;

bout_data.time       = LICKS_ALL_DATA.(['time_',alignment])(:,idx_bout);
bout_data.SS         = LICKS_ALL_DATA.(['neuro_SS_',alignment])(:,idx_bout);
bout_data.CS         = LICKS_ALL_DATA.(['neuro_CS_',alignment])(:,idx_bout);
bout_data.tongue_dm  = LICKS_ALL_DATA.(['tongue_dm_',alignment])(:,idx_bout);
bout_data.tongue_vm  = LICKS_ALL_DATA.(['tongue_vm_',alignment])(:,idx_bout);
bout_data.tongue_ang = LICKS_ALL_DATA.(['tongue_ang_',alignment])(:,idx_bout);

bout_data.tongue_dm_max  =  LICKS_ALL_DATA.tongue_dm_max(    :,idx_bout);
bout_data.tongue_vm_max  =  LICKS_ALL_DATA.tongue_vm_max(    :,idx_bout);
bout_data.tongue_vm_min  =  LICKS_ALL_DATA.tongue_vm_min(    :,idx_bout);
bout_data.tongue_ang_max =  LICKS_ALL_DATA.tongue_ang_max(    :,idx_bout);

bout_data.time_onset   = LICKS_ALL_DATA.time_onset (     :,idx_bout);
bout_data.time_vmax    = LICKS_ALL_DATA.time_vmax (     :,idx_bout);
bout_data.time_dmax    = LICKS_ALL_DATA.time_dmax   (     :,idx_bout);
bout_data.time_vmin    = LICKS_ALL_DATA.time_vmin (     :,idx_bout);
bout_data.time_offset  = LICKS_ALL_DATA.time_offset(     :,idx_bout);

bout_data.time_diff_onset_offset = bout_data.time_offset - bout_data.time_onset;
bout_data.time_diff_onset_dmax   = bout_data.time_dmax - bout_data.time_onset;
bout_data.time_diff_offset_dmax  = bout_data.time_offset - bout_data.time_onset;

idx_bout_ = LICKS_ALL_DATA.tag_bout(LICKS_ALL_DATA.tag_bout ~= 0);
idx_bout_ = idx_bout_(ismember(idx_bout_, tag_id)) ;

field_names_bout_data = fieldnames(bout_data);
bout_data_dir = struct;
for counter_dir = 1 : length(tag_id)
    for counter_field = 1 : length(field_names_bout_data)
        field_name = field_names_bout_data{counter_field};
        idx_tag = (idx_bout_ == counter_dir);
        bout_data_dir(counter_dir).(field_name) = bout_data.(field_name)(:,idx_tag);
    end
end

end

