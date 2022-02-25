%% function buildSacData
function sac_data_dir = ESN_buildSacData(SACS_ALL_DATA,tag_id,alignment,params,funcs)

ang_edges    = params.sac.ang_edges;
ang_values   = params.sac.ang_values;
length_trace = params.sac.length_trace;

% Build sac_data
idx_sacs    = ismember(SACS_ALL_DATA.tag, tag_id);
idx_sacs    = idx_sacs & SACS_ALL_DATA.validity;

switch(alignment)
    case('visual')
        sac_data.time    = SACS_ALL_DATA.time_visual(    :,idx_sacs);
        sac_data.SS      = SACS_ALL_DATA.neuro_SS_visual(:,idx_sacs);
        sac_data.CS      = SACS_ALL_DATA.neuro_CS_visual(:,idx_sacs);
        sac_data.eye_vx  = SACS_ALL_DATA.eye_vx_visual(:,idx_sacs);
        sac_data.eye_vy  = SACS_ALL_DATA.eye_vy_visual(:,idx_sacs);
    case('onset')
        sac_data.time    = SACS_ALL_DATA.time_onset(     :,idx_sacs);
        sac_data.SS      = SACS_ALL_DATA.neuro_SS_onset( :,idx_sacs);
        sac_data.CS      = SACS_ALL_DATA.neuro_CS_onset( :,idx_sacs);
        sac_data.eye_vx  = SACS_ALL_DATA.eye_vx_onset( :,idx_sacs);
        sac_data.eye_vy  = SACS_ALL_DATA.eye_vy_onset( :,idx_sacs);
    case('offset')
        sac_data.time    = SACS_ALL_DATA.time_offset(     :,idx_sacs);
        sac_data.SS      = SACS_ALL_DATA.neuro_SS_offset( :,idx_sacs);
        sac_data.CS      = SACS_ALL_DATA.neuro_CS_offset( :,idx_sacs);
        sac_data.eye_vx  = SACS_ALL_DATA.eye_vx_offset( :,idx_sacs);
        sac_data.eye_vy  = SACS_ALL_DATA.eye_vy_offset( :,idx_sacs);
    case('vmax')
        sac_data.time    = SACS_ALL_DATA.time_vmax(     :,idx_sacs);
        sac_data.SS      = SACS_ALL_DATA.neuro_SS_vmax( :,idx_sacs);
        sac_data.CS      = SACS_ALL_DATA.neuro_CS_vmax( :,idx_sacs);
        sac_data.eye_vx  = SACS_ALL_DATA.eye_vx_vmax( :,idx_sacs);
        sac_data.eye_vy  = SACS_ALL_DATA.eye_vy_vmax( :,idx_sacs);
    case('auditory')
        sac_data.time    = SACS_ALL_DATA.time_auditory(     :,idx_sacs);
        sac_data.SS      = SACS_ALL_DATA.neuro_SS_auditory( :,idx_sacs);
        sac_data.CS      = SACS_ALL_DATA.neuro_CS_auditory( :,idx_sacs);
        sac_data.eye_vx  = SACS_ALL_DATA.eye_vx_auditory( :,idx_sacs);
        sac_data.eye_vy  = SACS_ALL_DATA.eye_vy_auditory( :,idx_sacs);
end

sac_data.eye_vm = sqrt(sac_data.eye_vx.^2 + sac_data.eye_vy.^2);

sac_data.eye_r_px_onset   = SACS_ALL_DATA.eye_r_px_onset(:,idx_sacs);
sac_data.eye_r_py_onset   = SACS_ALL_DATA.eye_r_py_onset(:,idx_sacs);

sac_data.visual_px_offset = SACS_ALL_DATA.visual_px_offset(:,idx_sacs);
sac_data.visual_py_offset = SACS_ALL_DATA.visual_py_offset(:,idx_sacs);

sac_data.time_onset   = SACS_ALL_DATA.time_onset(     :,idx_sacs);
sac_data.time_offset  = SACS_ALL_DATA.time_offset(     :,idx_sacs);
sac_data.time_visual  = SACS_ALL_DATA.time_visual(     :,idx_sacs);

sac_data.time_diff_visual_onset = sac_data.time_onset - sac_data.time_visual;
sac_data.time_diff_onset_offset = sac_data.time_offset - sac_data.time_onset;

% Build sac_data_dir
sac_data.delta_x = sac_data.visual_px_offset - sac_data.eye_r_px_onset;
sac_data.delta_y = sac_data.visual_py_offset - sac_data.eye_r_py_onset;
sac_data.visual_ang = wrapTo360(atan2d(sac_data.delta_y, sac_data.delta_x));
sac_data.visual_ang_bin = discretize(sac_data.visual_ang, ang_edges);
last_bin_id = length(ang_edges) - 1;
sac_data.visual_ang_bin(sac_data.visual_ang_bin == last_bin_id) = 1; % wrap the circle around
% 1: 0deg % 2: 45deg % 3: 90deg % 4: 135deg % 5: 180deg % 6: 225deg % 7: 270deg % 8: 315deg

if length(ang_values) ~= 8
    error('plot_single_session_modulation: length ang_values is not 8. Please modify the code.')
end
if length_trace ~= 500
    error('sac_modulation_index: length_trace is not 500. Please modify the code.')
end

field_names_sac_data = fieldnames(sac_data);
sac_data_dir = struct;
for counter_dir = 1 : 8
    for counter_field = 1 : length(field_names_sac_data)
        field_name = field_names_sac_data{counter_field};
        idx_ang = (sac_data.visual_ang_bin == counter_dir);
        sac_data_dir(counter_dir).(field_name) = sac_data.(field_name)(:,idx_ang);
    end
end
end
