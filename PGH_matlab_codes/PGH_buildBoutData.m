%% function buildBoutData
function bout_data_dir = PGH_buildBoutData(LICKS_ALL_DATA, alignment, params, funcs)

if strcmp(alignment,'offset')
    tag_id = 2;
else
    tag_id = 1;
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

% determine if left or right bout
ind_bout_str = find(LICKS_ALL_DATA.tag_bout == 1);
ind_bout_end = find(LICKS_ALL_DATA.tag_bout == 2);
bout_dir = ones(1,numel(ind_bout_str)); % 1 groom 2 left, 3 right
for counter_bout = 1 : sum(idx_bout)
    inds = ind_bout_str(counter_bout) : ind_bout_end(counter_bout);
    if sum(LICKS_ALL_DATA.tongue_ang_max(inds) > 0) > sum(LICKS_ALL_DATA.tongue_ang_max(inds) < 0)
        bout_dir(1,counter_bout) = 3;
    elseif sum(LICKS_ALL_DATA.tongue_ang_max(inds) > 0) < sum(LICKS_ALL_DATA.tongue_ang_max(inds) < 0)
        bout_dir(1,counter_bout) = 2;
    else
        bout_dir(1,counter_bout) = 1;
    end
end

field_names_bout_data = fieldnames(bout_data);
bout_data_dir = struct;
for counter_dir = 1 : 3
    for counter_field = 1 : length(field_names_bout_data)
        field_name = field_names_bout_data{counter_field};
        bout_data_dir(counter_dir).(field_name) = bout_data.(field_name)(:,bout_dir==counter_dir);
    end
end

end

