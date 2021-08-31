function PH_plot_tongue_modulation(num_data_set, params)
if nargin < 1
    num_data_set = 1;
    params.manual = 1;
    params.auto = 1;
elseif nargin == 1
    params.manual = 1;
    params.auto = 1;
elseif nargin > 1
    params.auto = 1;
    params.manual = 0;
    
end
%% Build EPHYS_, VID_for each single dataset
clearvars EPHYS VID BEHAVE
for counter_dataset = 1 : 1 : num_data_set
    [EPHYS_(counter_dataset), VID_(counter_dataset)] = build_EPHYS_VID_single_dataset(num_data_set,params);
end

EPHYS.CH_sorted_file_name = EPHYS_(1).CH_sorted_file_name;
EPHYS.CH_sorted_file_path = EPHYS_(1).CH_sorted_file_path;

[~, file_name, ~]  = fileparts(EPHYS.CH_sorted_file_name);
if num_data_set > 1
    file_name = [file_name '_combine_' num2str(num_data_set)];
end
Neural_Properties_data.file_name = file_name;
Neural_Properties_data.file_path = EPHYS.CH_sorted_file_path;

%% Plot-1 SS & CS train lick_onset class
clearvars raster_data_
if counter_dataset > 1
    session_type = VID_(1).DLC.FILE.session_type;
else
    session_type = VID_.DLC.FILE.session_type;
end
for counter_dataset = 1 : 1 : num_data_set
    EPHYS_inds_event            = EPHYS_(counter_dataset).CH_EVE.EPHYS_LED_aligned_inds_lick_onset_100;
    VID_inds_event            = EPHYS_(counter_dataset).CH_EVE.VID_LED_aligned_inds_lick_onset_100;
    raster_data_(counter_dataset).data = single_dataset_raster(EPHYS_(counter_dataset), EPHYS_inds_event,...
        VID_(counter_dataset), VID_inds_event, session_type );
end

raster_data_lick_onset = concatenate_dataset(raster_data_, 'data', @vertcat);
plot_data_lick_onset.fig_num_               = 1;
plot_data_lick_onset.xlabel_text_raster_    = {'Lick onset (ms)'};
plot_data_lick_onset.xlabel_text_raster_bout_   = {'Lick onset (s)'};
plot_data_lick_onset.xlabel_text_CS_probab_ = {'CS prob. [-50 50]ms'};
plot_data_lick_onset.inds_span = EPHYS_(1).CH_EVE.inds_span_lick_onset;
fig_handle_(plot_data_lick_onset.fig_num_)  = plot_rasters_data_class(raster_data_lick_onset, plot_data_lick_onset, session_type);

sgtitle(fig_handle_(plot_data_lick_onset.fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');

%% Plot-2 SS & CS train lick_onset kin
clearvars raster_data_
if counter_dataset > 1
    session_type = VID_(1).DLC.FILE.session_type;
else
    session_type = VID_.DLC.FILE.session_type;
end
for counter_dataset = 1 : 1 : num_data_set
    EPHYS_inds_event            = EPHYS_(counter_dataset).CH_EVE.EPHYS_LED_aligned_inds_lick_onset_100;
    VID_inds_event            = EPHYS_(counter_dataset).CH_EVE.VID_LED_aligned_inds_lick_onset_100;
    raster_data_(counter_dataset).data = single_dataset_raster(EPHYS_(counter_dataset), EPHYS_inds_event,...
        VID_(counter_dataset), VID_inds_event, session_type );
end

raster_data_lick_onset = concatenate_dataset(raster_data_, 'data', @vertcat);
plot_data_lick_onset.fig_num_               = 2;
plot_data_lick_onset.xlabel_text_raster_    = {'Lick onset (ms)'};
plot_data_lick_onset.xlabel_text_raster_bout_   = {'Lick onset (s)'};
plot_data_lick_onset.xlabel_text_CS_probab_ = {'CS prob. [-50 50]ms'};
plot_data_lick_onset.inds_span = EPHYS_(1).CH_EVE.inds_span_lick_onset;
fig_handle_(plot_data_lick_onset.fig_num_)  = plot_rasters_data_kin(raster_data_lick_onset, plot_data_lick_onset, session_type);

sgtitle(fig_handle_(plot_data_lick_onset.fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');

%% Plot-3 SS & CS train lick_vmax
clearvars raster_data_
if counter_dataset > 1
    session_type = VID_(1).DLC.FILE.session_type;
else
    session_type = VID_.DLC.FILE.session_type;
end
for counter_dataset = 1 : 1 : num_data_set
    EPHYS_inds_event            = EPHYS_(counter_dataset).CH_EVE.EPHYS_LED_aligned_inds_lick_vmax_100;
    VID_inds_event            = EPHYS_(counter_dataset).CH_EVE.VID_LED_aligned_inds_lick_vmax_100;
    raster_data_(counter_dataset).data = single_dataset_raster(EPHYS_(counter_dataset), EPHYS_inds_event,...
        VID_(counter_dataset), VID_inds_event, session_type );
end

raster_data_lick_vmax= concatenate_dataset(raster_data_, 'data', @vertcat);

plot_data_lick_vmax.fig_num_               = 3;
plot_data_lick_vmax.xlabel_text_raster_    = {'Lick vmax (ms)'};
plot_data_lick_vmax.xlabel_text_raster_bout_   = {'Lick vmax (s)'};
plot_data_lick_vmax.xlabel_text_CS_probab_ = {'CS prob. [-50 50]ms'};
plot_data_lick_vmax.inds_span = EPHYS_(1).CH_EVE.inds_span_lick_vmax;
fig_handle_(plot_data_lick_vmax.fig_num_)  = plot_rasters_data_class(raster_data_lick_vmax, plot_data_lick_vmax, session_type);

sgtitle(fig_handle_(plot_data_lick_vmax.fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');

%% Plot-4 SS & CS train lick_vmax kin
clearvars raster_data_
if counter_dataset > 1
    session_type = VID_(1).DLC.FILE.session_type;
else
    session_type = VID_.DLC.FILE.session_type;
end
for counter_dataset = 1 : 1 : num_data_set
    EPHYS_inds_event            = EPHYS_(counter_dataset).CH_EVE.EPHYS_LED_aligned_inds_lick_vmax_100;
    VID_inds_event            = EPHYS_(counter_dataset).CH_EVE.VID_LED_aligned_inds_lick_vmax_100;
    raster_data_(counter_dataset).data = single_dataset_raster(EPHYS_(counter_dataset), EPHYS_inds_event,...
        VID_(counter_dataset), VID_inds_event, session_type );
end

raster_data_lick_vmax = concatenate_dataset(raster_data_, 'data', @vertcat);
plot_data_lick_vmax.fig_num_               = 4;
plot_data_lick_vmax.xlabel_text_raster_    = {'Lick vmax (ms)'};
plot_data_lick_vmax.xlabel_text_raster_bout_   = {'Lick vmax (s)'};
plot_data_lick_vmax.xlabel_text_CS_probab_ = {'CS prob. [-50 50]ms'};
plot_data_lick_vmax.inds_span = EPHYS_(1).CH_EVE.inds_span_lick_vmax;
fig_handle_(plot_data_lick_vmax.fig_num_)  = plot_rasters_data_kin(raster_data_lick_vmax, plot_data_lick_vmax, session_type);

sgtitle(fig_handle_(plot_data_lick_vmax.fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');

%% Plot-5 SS & CS train lick_dmax
clearvars raster_data_
if counter_dataset > 1
    session_type = VID_(1).DLC.FILE.session_type;
else
    session_type = VID_.DLC.FILE.session_type;
end
for counter_dataset = 1 : 1 : num_data_set
    EPHYS_inds_event            = EPHYS_(counter_dataset).CH_EVE.EPHYS_LED_aligned_inds_lick_dmax_100;
    VID_inds_event            = EPHYS_(counter_dataset).CH_EVE.VID_LED_aligned_inds_lick_dmax_100;
    raster_data_(counter_dataset).data = single_dataset_raster(EPHYS_(counter_dataset), EPHYS_inds_event,...
        VID_(counter_dataset), VID_inds_event, session_type);
end

raster_data_lick_dmax = concatenate_dataset(raster_data_, 'data', @vertcat);

plot_data_lick_dmax.fig_num_               = 5;
plot_data_lick_dmax.xlabel_text_raster_    = {'Lick dmax (ms)'};
plot_data_lick_dmax.xlabel_text_raster_bout_    = {'Lick dmax (s)'};

plot_data_lick_dmax.xlabel_text_CS_probab_ = {'CS prob. [-50 50]ms'};
plot_data_lick_dmax.inds_span = EPHYS_(1).CH_EVE.inds_span_lick_dmax;
fig_handle_(plot_data_lick_dmax.fig_num_)  = plot_rasters_data_class(raster_data_lick_dmax, plot_data_lick_dmax, session_type);

sgtitle(fig_handle_(plot_data_lick_dmax.fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');

%% Plot-6 SS & CS train lick_dmax kin
clearvars raster_data_
if counter_dataset > 1
    session_type = VID_(1).DLC.FILE.session_type;
else
    session_type = VID_.DLC.FILE.session_type;
end
for counter_dataset = 1 : 1 : num_data_set
    EPHYS_inds_event            = EPHYS_(counter_dataset).CH_EVE.EPHYS_LED_aligned_inds_lick_dmax_100;
    VID_inds_event            = EPHYS_(counter_dataset).CH_EVE.VID_LED_aligned_inds_lick_dmax_100;
    raster_data_(counter_dataset).data = single_dataset_raster(EPHYS_(counter_dataset), EPHYS_inds_event,...
        VID_(counter_dataset), VID_inds_event, session_type );
end

raster_data_lick_dmax = concatenate_dataset(raster_data_, 'data', @vertcat);
plot_data_lick_dmax.fig_num_               = 6;
plot_data_lick_dmax.xlabel_text_raster_    = {'Lick dmax (ms)'};
plot_data_lick_dmax.xlabel_text_raster_bout_   = {'Lick dmax (s)'};
plot_data_lick_dmax.xlabel_text_CS_probab_ = {'CS prob. [-50 50]ms'};
plot_data_lick_dmax.inds_span = EPHYS_(1).CH_EVE.inds_span_lick_dmax;
fig_handle_(plot_data_lick_dmax.fig_num_)  = plot_rasters_data_kin(raster_data_lick_dmax, plot_data_lick_dmax, session_type);

sgtitle(fig_handle_(plot_data_lick_dmax.fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');

%% Plot-7 SS & CS train lick_vmin
clearvars raster_data_
if counter_dataset > 1
    session_type = VID_(1).DLC.FILE.session_type;
else
    session_type = VID_.DLC.FILE.session_type;
end
for counter_dataset = 1 : 1 : num_data_set
    EPHYS_inds_event            = EPHYS_(counter_dataset).CH_EVE.EPHYS_LED_aligned_inds_lick_vmin_100;
    VID_inds_event            = EPHYS_(counter_dataset).CH_EVE.VID_LED_aligned_inds_lick_vmin_100;
    raster_data_(counter_dataset).data = single_dataset_raster(EPHYS_(counter_dataset), EPHYS_inds_event,...
        VID_(counter_dataset), VID_inds_event, session_type );
end

raster_data_lick_vmin= concatenate_dataset(raster_data_, 'data', @vertcat);

plot_data_lick_vmin.fig_num_               = 7;
plot_data_lick_vmin.xlabel_text_raster_    = {'Lick vmin (ms)'};
plot_data_lick_vmin.xlabel_text_raster_bout_   = {'Lick vmin (s)'};
plot_data_lick_vmin.xlabel_text_CS_probab_ = {'CS prob. [-50 50]ms'};
plot_data_lick_vmin.inds_span = EPHYS_(1).CH_EVE.inds_span_lick_vmin;
fig_handle_(plot_data_lick_vmin.fig_num_)  = plot_rasters_data_class(raster_data_lick_vmin, plot_data_lick_vmin, session_type);

sgtitle(fig_handle_(plot_data_lick_vmin.fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');

%% Plot-8 SS & CS train lick_vmin kin
clearvars raster_data_
if counter_dataset > 1
    session_type = VID_(1).DLC.FILE.session_type;
else
    session_type = VID_.DLC.FILE.session_type;
end
for counter_dataset = 1 : 1 : num_data_set
    EPHYS_inds_event            = EPHYS_(counter_dataset).CH_EVE.EPHYS_LED_aligned_inds_lick_vmin_100;
    VID_inds_event            = EPHYS_(counter_dataset).CH_EVE.VID_LED_aligned_inds_lick_vmin_100;
    raster_data_(counter_dataset).data = single_dataset_raster(EPHYS_(counter_dataset), EPHYS_inds_event,...
        VID_(counter_dataset), VID_inds_event, session_type );
end

raster_data_lick_vmin = concatenate_dataset(raster_data_, 'data', @vertcat);
plot_data_lick_vmin.fig_num_               = 8;
plot_data_lick_vmin.xlabel_text_raster_    = {'Lick vmin (ms)'};
plot_data_lick_vmin.xlabel_text_raster_bout_   = {'Lick vmin (s)'};
plot_data_lick_vmin.xlabel_text_CS_probab_ = {'CS prob. [-50 50]ms'};
plot_data_lick_vmin.inds_span = EPHYS_(1).CH_EVE.inds_span_lick_vmin;
fig_handle_(plot_data_lick_vmin.fig_num_)  = plot_rasters_data_kin(raster_data_lick_vmin, plot_data_lick_vmin, session_type);

sgtitle(fig_handle_(plot_data_lick_vmin.fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');

%% Plot-9 SS & CS train lick_offset
clearvars raster_data_
if counter_dataset > 1
    session_type = VID_(1).DLC.FILE.session_type;
else
    session_type = VID_.DLC.FILE.session_type;
end
for counter_dataset = 1 : 1 : num_data_set
    EPHYS_inds_event            = EPHYS_(counter_dataset).CH_EVE.EPHYS_LED_aligned_inds_lick_offset_100;
    VID_inds_event            = EPHYS_(counter_dataset).CH_EVE.VID_LED_aligned_inds_lick_offset_100;
    raster_data_(counter_dataset).data = single_dataset_raster(EPHYS_(counter_dataset), EPHYS_inds_event,...
        VID_(counter_dataset), VID_inds_event, session_type );
end

raster_data_lick_offset = concatenate_dataset(raster_data_, 'data', @vertcat);

plot_data_lick_offset.fig_num_               = 9;
plot_data_lick_offset.xlabel_text_raster_    = {'Lick offset (ms)'};
plot_data_lick_offset.xlabel_text_raster_bout_    = {'Lick offset (s)'};
plot_data_lick_offset.xlabel_text_CS_probab_ = {'CS prob. [-50 50]ms'};
plot_data_lick_offset.inds_span = EPHYS_(1).CH_EVE.inds_span_lick_offset;
fig_handle_(plot_data_lick_offset.fig_num_)  = plot_rasters_data_class(raster_data_lick_offset, plot_data_lick_offset, session_type);

sgtitle(fig_handle_(plot_data_lick_offset.fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');

%% Plot-10 SS & CS train lick_offset kin
clearvars raster_data_
if counter_dataset > 1
    session_type = VID_(1).DLC.FILE.session_type;
else
    session_type = VID_.DLC.FILE.session_type;
end
for counter_dataset = 1 : 1 : num_data_set
    EPHYS_inds_event            = EPHYS_(counter_dataset).CH_EVE.EPHYS_LED_aligned_inds_lick_offset_100;
    VID_inds_event            = EPHYS_(counter_dataset).CH_EVE.VID_LED_aligned_inds_lick_offset_100;
    raster_data_(counter_dataset).data = single_dataset_raster(EPHYS_(counter_dataset), EPHYS_inds_event,...
        VID_(counter_dataset), VID_inds_event, session_type );
end

raster_data_lick_offset = concatenate_dataset(raster_data_, 'data', @vertcat);
plot_data_lick_offset.fig_num_               = 10;
plot_data_lick_offset.xlabel_text_raster_    = {'Lick offset (ms)'};
plot_data_lick_offset.xlabel_text_raster_bout_   = {'Lick offset (s)'};
plot_data_lick_offset.xlabel_text_CS_probab_ = {'CS prob. [-50 50]ms'};
plot_data_lick_offset.inds_span = EPHYS_(1).CH_EVE.inds_span_lick_offset;
fig_handle_(plot_data_lick_offset.fig_num_)  = plot_rasters_data_kin(raster_data_lick_offset, plot_data_lick_offset, session_type);

sgtitle(fig_handle_(plot_data_lick_offset.fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');

%% Plot-11 Neural Properties
CH_sorted_ = concatenate_dataset(EPHYS_, 'CH_sorted', @horzcat);
EPHYS.CH_sorted.SS_data   = concatenate_dataset(CH_sorted_.SS_data, [], @vertcat);
EPHYS.CH_sorted.CS_data   = concatenate_dataset(CH_sorted_.CS_data, [], @vertcat);
EPHYS.CH_sorted.Corr_data = concatenate_dataset(CH_sorted_.Corr_data, [], @vertcat);

if isfield(EPHYS.CH_sorted.SS_data, 'SS_waveform_hipass')
    EPHYS.CH_sorted.SS_data.SS_waveform = EPHYS.CH_sorted.SS_data.SS_waveform_hipass;
    EPHYS.CH_sorted.CS_data.CS_waveform = EPHYS.CH_sorted.CS_data.CS_waveform_hipass;
end

Neural_Properties_data.CH_sorted.SS_data   = EPHYS.CH_sorted.SS_data;
Neural_Properties_data.CH_sorted.CS_data   = EPHYS.CH_sorted.CS_data;
Neural_Properties_data.CH_sorted.Corr_data = EPHYS.CH_sorted.Corr_data;

% figure
Line_Color = lines(7);
plot_data.fig_num_ = 11;
fig_handle_(plot_data.fig_num_) = figure(plot_data.fig_num_);
clf(fig_handle_(plot_data.fig_num_))

% subplot(2,2,1) Waveform
plot_handle_(1) = subplot(2,2,1);
hold on

plot((1:180)/30-2, mean(Neural_Properties_data.CH_sorted.SS_data.SS_waveform)+std(Neural_Properties_data.CH_sorted.SS_data.SS_waveform), ...
    'LineWidth', 1, 'Color', Line_Color(1,:))
plot((1:180)/30-2, mean(Neural_Properties_data.CH_sorted.SS_data.SS_waveform)-std(Neural_Properties_data.CH_sorted.SS_data.SS_waveform), ...
    'LineWidth', 1, 'Color', Line_Color(1,:))
plot((1:180)/30-2, mean(Neural_Properties_data.CH_sorted.SS_data.SS_waveform), ...
    'LineWidth', 2, 'Color', Line_Color(1,:))
plot((1:180)/30-2, mean(Neural_Properties_data.CH_sorted.CS_data.CS_waveform)+std(Neural_Properties_data.CH_sorted.CS_data.CS_waveform), ...
    'LineWidth', 1, 'Color', Line_Color(7,:))
plot((1:180)/30-2, mean(Neural_Properties_data.CH_sorted.CS_data.CS_waveform)-std(Neural_Properties_data.CH_sorted.CS_data.CS_waveform), ...
    'LineWidth', 1, 'Color', Line_Color(7,:))
plot((1:180)/30-2, mean(Neural_Properties_data.CH_sorted.CS_data.CS_waveform), ...
    'LineWidth', 2, 'Color', Line_Color(7,:))

xlabel('Time (ms)')
ylabel('Voltage (uv)')
title('CS & SS Waveform')

% subplot(2,2,2) Probablities
plot_handle_(2) = subplot(2,2,2);
hold on

SSxSS_AUTO = Neural_Properties_data.CH_sorted.Corr_data.SS_SSxSS_AUTO;
if size(SSxSS_AUTO, 1) > 1
    inds_span = mean(Neural_Properties_data.CH_sorted.Corr_data.SS_inds_span);
else
    inds_span =      Neural_Properties_data.CH_sorted.Corr_data.SS_inds_span;
end
bin_size_time = mean(Neural_Properties_data.CH_sorted.Corr_data.SS_bin_size_time);
if (~isempty(SSxSS_AUTO))
    x_axis_ = (inds_span * bin_size_time * 1000)';
    if size(SSxSS_AUTO, 1) > 1
        prob_value_ = mean(SSxSS_AUTO);
    else
        prob_value_ =      SSxSS_AUTO;
    end
    prob_value_(round(length(inds_span)/2)) = 0;
    y_axis_mean_ = prob_value_;
    y_axis_stdv_ = sqrt(size(SSxSS_AUTO, 1) .* ( prob_value_ .* (1 - prob_value_) ) ) ./ size(SSxSS_AUTO, 1);
    xlim([x_axis_(1) x_axis_(end)]);
else
    x_axis_      = nan(size(SSxSS_AUTO));
    y_axis_mean_ = nan(size(SSxSS_AUTO));
    y_axis_stdv_ = nan(size(SSxSS_AUTO));
end
plot(x_axis_, y_axis_mean_+y_axis_stdv_, 'LineWidth', 1, 'Color', Line_Color(1,:))
plot(x_axis_, y_axis_mean_-y_axis_stdv_, 'LineWidth', 1, 'Color', Line_Color(1,:))
plot(x_axis_, y_axis_mean_,              'LineWidth', 2, 'Color', Line_Color(1,:))

CSxSS_AUTO    = Neural_Properties_data.CH_sorted.Corr_data.CS_CSxSS_AUTO;
if size(CSxSS_AUTO, 1) > 1
    inds_span = mean(Neural_Properties_data.CH_sorted.Corr_data.CS_inds_span);
else
    inds_span      = Neural_Properties_data.CH_sorted.Corr_data.CS_inds_span;
end
bin_size_time = mean(Neural_Properties_data.CH_sorted.Corr_data.CS_bin_size_time);
if (~isempty(CSxSS_AUTO))
    x_axis_ = (inds_span * bin_size_time * 1000)';
    if size(CSxSS_AUTO, 1) > 1
        prob_value_ = mean(CSxSS_AUTO);
    else
        prob_value_ =      CSxSS_AUTO;
    end
    y_axis_mean_ = prob_value_;
    y_axis_stdv_ = sqrt(size(CSxSS_AUTO, 1) .* ( prob_value_ .* (1 - prob_value_) ) ) ./ size(CSxSS_AUTO, 1);
    xlim([x_axis_(1) x_axis_(end)]);
else
    x_axis_      = nan(size(CSxSS_AUTO));
    y_axis_mean_ = nan(size(CSxSS_AUTO));
    y_axis_stdv_ = nan(size(CSxSS_AUTO));
end
plot(x_axis_, y_axis_mean_+y_axis_stdv_, 'LineWidth', 1, 'Color', Line_Color(7,:))
plot(x_axis_, y_axis_mean_-y_axis_stdv_, 'LineWidth', 1, 'Color', Line_Color(7,:))
plot(x_axis_, y_axis_mean_,              'LineWidth', 2, 'Color', Line_Color(7,:))

ylim([0 inf])
xlabel('Time (ms)')
ylabel('Probability')
title('X Probability')

% subplot(2,2,3) ISI SS
plot_handle_(3) = subplot(2,2,3);
hold on

edges_SS = (0 : 0.002 : 0.050) *1000;
ISI_SS = diff(Neural_Properties_data.CH_sorted.SS_data.SS_time) * 1000;
histogram(ISI_SS,edges_SS, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', Line_Color(1,:));
set(plot_handle_(3), 'XTick', [0 0.025 0.050]*1000)
xlabel('Time (ms)')
ylabel('Probability')
title('SS Inter-Spike-Interval')

% subplot(2,2,4) ISI CS
plot_handle_(4) = subplot(2,2,4);
hold on
edges_CS = 0 : 0.200 : 5.000;
ISI_CS = diff(Neural_Properties_data.CH_sorted.CS_data.CS_time);
histogram(ISI_CS,edges_CS, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', Line_Color(7,:));
set(plot_handle_(4), 'XTick', [0 2.5 5.0])
xlabel('Time (s)')
ylabel('Probability')
title('CS Inter-Spike-Interval')

ESN_Beautify_Plot
sgtitle(fig_handle_(plot_data.fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');

%% Plot-12 Kinematic Properties
clearvars kinematic_data_

if counter_dataset > 1
    session_type = VID_(1).DLC.FILE.session_type;
else
    session_type = VID_.DLC.FILE.session_type;
end

for counter_dataset = 1 : 1 : num_data_set
    kinematic_data_(counter_dataset).data = single_dataset_kinematic(EPHYS_(counter_dataset), VID_(counter_dataset), session_type );
end

kinematic_data_ = concatenate_dataset(kinematic_data_, [], @vertcat);
kinematic_data = concatenate_dataset(kinematic_data_.data, [],  @vertcat);

% figure
plot_data.fig_num_ = 12;
fig_handle_(plot_data.fig_num_) = figure(plot_data.fig_num_);
fig_handle_(plot_data.fig_num_).WindowState = 'maximized';
clf(fig_handle_(plot_data.fig_num_))

plot_handle_(1) = subplot(3,4,1);
hold on
plot(kinematic_data.VID_d_max_100_grooming, '.g', 'MarkerSize',5)
% plot(fitlm(1:length(kinematic_data.VID_d_max_100_grooming),(kinematic_data.VID_d_max_100_grooming), 'linear'))
plot(kinematic_data.VID_d_max_100_r, '.r', 'MarkerSize', 5)
if (session_type  == 2)
    plot(kinematic_data.VID_d_max_100_l, '.b', 'MarkerSize', 5)
end
title('Displacement')
xlabel('Lick #')
ylabel('mm')
ylim([0 25])

plot_handle_(2) = subplot(3,4,2);
hold on
plot(kinematic_data.VID_v_max_100_grooming, '.g', 'MarkerSize', 5)
plot(kinematic_data.VID_v_max_100_r, '.r', 'MarkerSize', 5)
if (session_type  == 2)
    plot(kinematic_data.VID_v_max_100_l, '.b', 'MarkerSize', 5)
end
plot(kinematic_data.VID_v_min_100_grooming, '.g', 'MarkerSize', 5)
plot(kinematic_data.VID_v_min_100_r, '.r', 'MarkerSize', 5)
if (session_type  == 2)
    plot(kinematic_data.VID_v_min_100_l, '.b', 'MarkerSize', 5)
end
title('Velocity')
ylabel('mm/s')
xlabel('Lick #')
ylim([-600 600])

plot_handle_(3) = subplot(3,4,3);
hold on
plot(abs(kinematic_data.VID_angle_max_100_grooming), '.g', 'MarkerSize', 5)
plot(kinematic_data.VID_angle_max_100_r, '.r', 'MarkerSize', 5)
if (session_type  == 2)
    plot(abs(kinematic_data.VID_angle_max_100_l), '.b', 'MarkerSize', 5)
end
title('Angle')
ylabel('deg')
xlabel('Lick #')
ylim([-20 120])

plot_handle_(4) = subplot(3,4,4);
hold on
plot(kinematic_data.VID_lick_duration_grooming*1000, '.g', 'MarkerSize', 5)
plot(kinematic_data.VID_lick_duration_r*1000, '.r', 'MarkerSize', 5)
if (session_type  == 2)
    plot(kinematic_data.VID_lick_duration_l*1000, '.b', 'MarkerSize', 5)
end
title('Duration')
ylabel('ms')
xlabel('Lick #')
ylim([0 500])

plot_handle_(5) = subplot(3,4,5);
hold on
plot(kinematic_data.VID_d_max_100_grooming,kinematic_data.VID_v_max_100_grooming, '.g', 'MarkerSize', 5)
plot(kinematic_data.VID_d_max_100_r,kinematic_data.VID_v_max_100_r, '.r', 'MarkerSize', 5)
if (session_type  == 2)
    plot(kinematic_data.VID_d_max_100_l,kinematic_data.VID_v_max_100_l, '.b', 'MarkerSize', 5)
end
plot(kinematic_data.VID_d_max_100_grooming,kinematic_data.VID_v_min_100_grooming, '.g', 'MarkerSize', 5)
plot(kinematic_data.VID_d_max_100_r,kinematic_data.VID_v_min_100_r, '.r', 'MarkerSize', 5)
if (session_type  == 2)
    plot(kinematic_data.VID_d_max_100_l,kinematic_data.VID_v_min_100_l, '.b', 'MarkerSize', 5)
end
title('Displacement vs Velocity')
ylabel('mm/s')
xlabel('mm')
ylim([-600 600])
xlim([0 25])

plot_handle_(6) = subplot(3,4,6);
hold on
plot(kinematic_data.VID_d_max_100_grooming,abs(kinematic_data.VID_angle_max_100_grooming), '.g', 'MarkerSize', 5)
plot(kinematic_data.VID_d_max_100_r,kinematic_data.VID_angle_max_100_r, '.r', 'MarkerSize', 5)
if (session_type  == 2)
    plot(kinematic_data.VID_d_max_100_l,abs(kinematic_data.VID_angle_max_100_l), '.b', 'MarkerSize', 5)
end
title('Displacement vs Angle')
ylabel('deg')
xlabel('mm')
ylim([0 100])
xlim([0 25])

plot_handle_(7) = subplot(3,4,7);
hold on
plot(abs(kinematic_data.VID_angle_max_100_grooming),kinematic_data.VID_v_max_100_grooming, '.g', 'MarkerSize', 5)
plot(kinematic_data.VID_angle_max_100_r,kinematic_data.VID_v_max_100_r, '.r', 'MarkerSize', 5)
if (session_type  == 2)
    plot(abs(kinematic_data.VID_angle_max_100_l),kinematic_data.VID_v_max_100_l, '.b', 'MarkerSize', 5)
end
plot(abs(kinematic_data.VID_angle_max_100_grooming),kinematic_data.VID_v_min_100_grooming, '.g', 'MarkerSize', 5)
plot(kinematic_data.VID_angle_max_100_r,kinematic_data.VID_v_min_100_r, '.r', 'MarkerSize', 5)
if (session_type  == 2)
    plot(abs(kinematic_data.VID_angle_max_100_l),kinematic_data.VID_v_min_100_l, '.b', 'MarkerSize', 5)
end
title('Angle vs Velocity')
ylabel('mm/s')
xlabel('deg')
ylim([-600 600])
xlim([0 100])

plot_handle_(8) = subplot(3,4,8);
hold on
plot(kinematic_data.VID_d_max_100_grooming,kinematic_data.VID_lick_duration_grooming*1000, '.g', 'MarkerSize', 5)
plot(kinematic_data.VID_d_max_100_r,kinematic_data.VID_lick_duration_r*1000, '.r', 'MarkerSize', 5)
if (session_type  == 2)
    plot(kinematic_data.VID_d_max_100_l,kinematic_data.VID_lick_duration_l*1000, '.b', 'MarkerSize', 5)
end
title('Displacement vs Duration')
ylabel('ms')
xlabel('mm')
ylim([0 500])
xlim([0 25])

plot_handle_(9) = subplot(3,4,9);
hold on
plot(kinematic_data.VID_num_lick_bout(logical(kinematic_data.is_bout_r)), '.r')
plot(kinematic_data.VID_num_lick_bout(logical(kinematic_data.is_bout_l)), '.b')
title('Bout # of licks')
xlabel('Bout #')
ylabel('Licks')

plot_handle_(10) = subplot(3,4,10);
hold on
plot(kinematic_data.VID_bout_duration(logical(kinematic_data.is_bout_r)), '.r' )
plot(kinematic_data.VID_bout_duration(logical(kinematic_data.is_bout_l)), '.b' )
title('Bout duration')
xlabel('Bout #')
ylabel('Time (s)')

plot_handle_(11) = subplot(3,4,11);
hold on
edges_ILI_all = (0 : 20 : 500);
% ILI_all = kinematic_data.VID_ILI_all*1000;
ILI_grooming = kinematic_data.VID_ILI_grooming*1000;
ILI_r = kinematic_data.VID_ILI_r*1000;
ILI_l = kinematic_data.VID_ILI_l*1000;
% histogram(ILI_all,edges_ILI_all, 'Displaystyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', 'k');
histogram(ILI_grooming,edges_ILI_all, 'Displaystyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', 'g');
histogram(ILI_r,edges_ILI_all, 'Displaystyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', 'r');
histogram(ILI_l,edges_ILI_all, 'Displaystyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', 'b');
set(plot_handle_(11), 'XTick', [0 250 500])
xlabel('Time (ms)')
ylabel('Probability')
title('Inter lick interval')
if isempty(ILI_grooming) == 0
xline(nanmedian(ILI_grooming), 'g', 'LineWidth', 2);
end
if isempty(ILI_r) == 0
xline(nanmedian(ILI_r), 'r', 'LineWidth', 2); 
end
if isempty(ILI_l) == 0
xline(nanmedian(ILI_l), 'b', 'LineWidth', 2);
end

plot_handle_(12) = subplot(3,4,12);
hold on
edges_ILR_all = (0 : 0.25 : 5);
% ILR_all = kinematic_data.VID_ILR_all;
ILR_grooming = kinematic_data.VID_ILR_grooming;
ILR_r = kinematic_data.VID_ILR_r;
ILR_l = kinematic_data.VID_ILR_l;
% histogram(ILR_all,edges_ILR_all, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', 'k');
histogram(ILR_grooming,edges_ILR_all, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', 'g');
histogram(ILR_r,edges_ILR_all, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', 'r');
histogram(ILR_l,edges_ILR_all, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', 'b');
set(plot_handle_(12), 'XTick', [0 2.5 5])
xlabel('Frequency (Hz)')
ylabel('Probability')
title('Instantaneous lick rate')
if isempty(ILR_grooming) == 0
xline(nanmedian(ILR_grooming), 'g', 'LineWidth', 2);
end
if isempty(ILR_r) == 0
xline(nanmedian(ILR_r), 'r', 'LineWidth', 2); 
end
if isempty(ILR_l) == 0
xline(nanmedian(ILR_l), 'b', 'LineWidth', 2);
end

sgtitle(fig_handle_(plot_data.fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');
ESN_Beautify_Plot
%% Plot-13 Signals
% clearvars signal_data_
%
% if counter_dataset > 1
%     session_type = VID_(1).DLC.FILE.session_type;
% else
%     session_type = VID_.DLC.FILE.session_type;
% end
%
% for counter_dataset = 1 : 1 : num_data_set
%     signal_data_(counter_dataset).data = single_dataset_signal(EPHYS_(counter_dataset), VID_(counter_dataset), session_type );
% end
%
% signal_data_ = concatenate_dataset(signal_data_, [], @vertcat);
% signal_data = concatenate_dataset(signal_data_.data, [],  @vertcat);
%
% if length(signal_data.EPHYS_xcorr_time_100) > length(signal_data.VID_d_tip_100)
%     signal_data.EPHYS_xcorr_time_100 = signal_data.EPHYS_xcorr_time_100(1:length(signal_data.VID_d_tip_100));
% elseif length(signal_data.EPHYS_xcorr_time_100) < length(signal_data.VID_d_tip_100)
%     signal_data.VID_d_tip_100 =  signal_data.VID_d_tip_100(1:length(signal_data.EPHYS_xcorr_time_100));
% end
%
% figure
% plot_data.fig_num_ = 13;
% fig_handle_(plot_data.fig_num_) = figure(plot_data.fig_num_);
% clf(fig_handle_(plot_data.fig_num_))
%
% hold on
% plot(signal_data.EPHYS_time_30K,signal_data.EPHYS_signal, 'k')
% plot(signal_data.EPHYS_CS_time, 1000, '*r')
% % plot(EPHYS_SS_time, 2000, '.k')
% ylabel('Voltage uV')
% ylim([-5000 2000])
% yticks([-1000 0 1000])
% yyaxis right;
% plot(signal_data.EPHYS_xcorr_time_100,signal_data.VID_d_tip_100,'b', 'LineWidth', 2)
% ylabel('Tongue displacement (mm)')
% xlabel('Time (s)')
% yticks([0 10 20])
% ylim([0 35])
% set(gca, 'YColor', 'b')
%
% sgtitle(fig_handle_(plot_data.fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');

%% Save Fig
if params.manual
    response_save_fig = questdlg('Do you want to save the figures?',...
        'Question Dialog','Yes','No','Yes');
else
    response_save_fig = 'No';
    close('all')
end
if contains(response_save_fig, 'Yes')
    fprintf(['Saving plots', ' ...'])
    save_file_path = uigetdir(EPHYS.CH_sorted_file_path, 'Select where to save the figures.');
    if ~isequal(save_file_path,0)
        saveas(fig_handle_(1),[save_file_path filesep file_name '_modulation_lick_onset'], 'pdf');
        saveas(fig_handle_(2),[save_file_path filesep file_name '_modulation_lick_onset_kin'], 'pdf');
        saveas(fig_handle_(3),[save_file_path filesep file_name '_modulation_lick_vmax'], 'pdf');
        saveas(fig_handle_(4),[save_file_path filesep file_name '_modulation_lick_vmax_kin'], 'pdf');
        saveas(fig_handle_(5),[save_file_path filesep file_name '_modulation_lick_dmax'], 'pdf');
        saveas(fig_handle_(6),[save_file_path filesep file_name '_modulation_lick_dmax_kin'], 'pdf');
        saveas(fig_handle_(7),[save_file_path filesep file_name '_modulation_lick_vmin'], 'pdf');
        saveas(fig_handle_(8),[save_file_path filesep file_name '_modulation_lick_vmin_kin'], 'pdf');        
        saveas(fig_handle_(9),[save_file_path filesep file_name '_modulation_lick_offset'], 'pdf');
        saveas(fig_handle_(10),[save_file_path filesep file_name '_modulation_lick_offset_kin'], 'pdf');
        %saveas(fig_handle_(11),[save_file_path filesep file_name '_neural_properties'], 'pdf');
        saveas(fig_handle_(12),[save_file_path filesep file_name '_lick_properties'], 'pdf');
        %close(fig_handle_(1))

    end
    fprintf(' --> Completed. \n')
end % if contains(response_save_fig, 'Yes')

%% Report Properties
for counter_dataset = 1 : 1 : num_data_set
    [~, file_name, ~]  = fileparts(EPHYS_(counter_dataset).CH_sorted_file_name);
    duration = (EPHYS_(counter_dataset).CH_EVE.EPHYS_time_30K(end)-EPHYS_(counter_dataset).CH_EVE.EPHYS_time_30K(1));
    numCS = length(EPHYS_(counter_dataset).CH_sorted.CS_data.CS_ind);
    freqCS = numCS/duration;
    numSS = length(EPHYS_(counter_dataset).CH_sorted.SS_data.SS_ind);
    freqSS = numSS/duration;
    numLick = VID_(counter_dataset).DLC.IND.num_lick;
    fprintf(['*******************************************' '\n'])
    fprintf([file_name '\n'])
    fprintf([       'dur'   '\t'        'numCS'   '\t'        'freqCS'   '\t'        'numSS'   '\t'        'freqSS'   '\t'        'numLick'   '\n'])
    fprintf([num2str(duration/60,'%.1f') '\t'  num2str(numCS,'%.0f') '\t' num2str(freqCS,'%.2f') '\t' num2str(numSS,'%.0f') '\t' num2str(freqSS,'%.2f') '\t' num2str(numLick,'%.0f') '\n'])
end

%% Save plot data
if params.manual
    response_save_data = questdlg('Do you want to save the plot_data?',...
        'Question Dialog','Yes','No','Yes');
else
    response_save_data = 'Yes';
end
if contains(response_save_data, 'Yes')
    file_name = [Neural_Properties_data.file_name '_plot_data_lick.mat'];
    file_path = Neural_Properties_data.file_path;
    if ~params.auto
        [save_file_name,save_file_path] = uiputfile([file_path filesep file_name], 'Select where to save the plot data.');
    else
        save_file_name = file_name;
        save_file_path = [file_path filesep '..' filesep 'analyzed_figs'];
    end
    fprintf(['Saving ' save_file_name ' ... ']);
    if ~isequal(save_file_name,0)
        save([save_file_path filesep save_file_name], 'Neural_Properties_data', ...
            'raster_data_lick_onset',    'plot_data_lick_onset',...
            'raster_data_lick_vmax',  'plot_data_lick_vmax', ...
            'raster_data_lick_dmax',  'plot_data_lick_dmax', ...
            'raster_data_lick_vmin', 'plot_data_lick_vmin' , ...
            'raster_data_lick_offset',  'plot_data_lick_offset', ...
            'kinematic_data', ...
            '-v7.3');
    end
    fprintf(' --> Completed. \n');
    
    
end % if contains(response_save_data, 'Yes')

end

%% function ESN_raster_plot_axes
function [x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, x_axis_values, line_half_len)
if nargin < 2
    x_axis_values = 1 : 1 : size(train_data_logic, 2);
    line_half_len = 0.5;
end
if nargin < 3
    line_half_len = 0.5;
end
train_data_logic = train_data_logic > 0.1;
train_data_row_number = nan(size(train_data_logic));
for counter_row = 1 : size(train_data_logic, 1)
    train_data_row_number(counter_row, train_data_logic(counter_row,:)) = counter_row;
end
train_data_col_number = repmat(x_axis_values(:)', size(train_data_logic,1), 1);
x_axis = [train_data_col_number(:)'; train_data_col_number(:)'; nan(length(train_data_col_number(:)), 1)'];
y_axis = [(train_data_row_number(:)-line_half_len)'; (train_data_row_number(:)+line_half_len)'; nan(length(train_data_row_number(:)), 1)'];
x_axis = x_axis(:);
y_axis = y_axis(:);
end

%% function ESN_smooth
function smooth_data_ = ESN_smooth(data_, size)
% method = 'moving';  % Moving average. A lowpass filter with filter coefficients equal to the reciprocal of the span.
% method = 'lowess';  % Local regression using weighted linear least squares and a 1st degree polynomial model.
% method = 'loess';   % Local regression using weighted linear least squares and a 2nd degree polynomial model.
% method = 'sgolay';  % Savitzky-Golay filter. A generalized moving average with filter coefficients determined by an unweighted linear least-squares regression and a polynomial model of specified degree (default is 2). The method can accept nonuniform predictor data.
% method = 'rlowess'; % A robust version of 'lowess' that assigns lower weight to outliers in the regression. The method assigns zero weight to data outside six mean absolute deviations.
% method = 'rloess';  % A robust version of 'loess' that assigns lower weight to outliers in the regression. The method assigns zero weight to data outside six mean absolute deviations.
% smooth_data_ = smooth(data_, method);
smooth_data_ = smooth(data_, size, 'sgolay', 2);
end

%% function build_EPHYS_VID_single_dataset
function [EPHYS, VID, BEHAVE] = build_EPHYS_VID_single_dataset(num_data_set, params)

%% load EPHYS sorted DATA
if ~params.auto && num_data_set == 1
    file_path = [pwd filesep];
    [file_name, file_path] = uigetfile([file_path '.psort'], 'Select psort file');
elseif params.auto && num_data_set > 1
    [params.file_name, params.file_path] = uigetfile('.psort', 'Select psort file');
    file_name = params.file_name;
    file_path = params.file_path;
elseif params.auto && ~params.manual && num_data_set == 1
    file_name = [params.file_name '.psort'];
    file_path = params.file_path;
elseif params.auto && params.manual && num_data_set == 1
    [params.file_name, params.file_path] = uigetfile('.psort', 'Select psort file');
    file_name = params.file_name;
    file_path = params.file_path;
    
end
fprintf(['Loading ', file_name, ' ... ']);
% EPHYS.CH_sorted = load([file_path filesep file_name], 'CS_data', 'SS_data');
DATA_PSORT = Psort_read_psort([file_path file_name]);
EPHYS.CH_sorted_file_name = file_name;
EPHYS.CH_sorted_file_path = file_path;
fprintf(' --> Completed. \n')

%% load PH_aligned DATA
if ~params.auto
    %     file_name = EPHYS.CH_sorted_file_name(1:13);
    [file_name,file_path] = uigetfile([file_path filesep '_EVE1_PH_aligned.mat'], 'Select PH_aligned DATA file');
else
    file_name = [EPHYS.CH_sorted_file_name(1:13) '_EVE1_PH_aligned.mat'];
    file_path = params.file_path;
end
fprintf(['Loading ', file_name, ' ... ']);
EPHYS.CH_EVE = load([file_path file_name]);

% EPHYS.CH_EVE.EPHYS_time_vid  = EPHYS.CH_EVE.EPHYS_time_vid(:);
% EPHYS.CH_EVE.VID_time_vid = EPHYS.CH_EVE.VID_time_vid(:);

EPHYS.CH_EVE.EPHYS_time_100  = EPHYS.CH_EVE.EPHYS_time_100(:);
EPHYS.CH_EVE.VID_time_100 = EPHYS.CH_EVE.VID_time_100(:);
fprintf(' --> Completed. \n')

%% load ESN_aligned DATA
if ~params.auto
    %     file_name = EPHYS.CH_sorted_file_name(1:13);
    [file_name,file_path] = uigetfile([file_path filesep '_EVE1_ESN_aligned.mat'], 'Select ESN_aligned DATA file');
else
    file_name = [EPHYS.CH_sorted_file_name(1:13) '_EVE1_ESN_aligned.mat'];
    file_path = params.file_path;
end
fprintf(['Loading ', file_name, ' ... ']);
EPHYS.CH_EVE_ = load([file_path file_name]);
if isfield(EPHYS.CH_EVE, 'EPHYS_time_15K')
    EPHYS.CH_EVE_.EPHYS_time_30K = EPHYS.CH_EVE_.EPHYS_time_15K(:);
else
    EPHYS.CH_EVE_.EPHYS_time_30K = EPHYS.CH_EVE_.EPHYS_time_30K(:);
end
EPHYS.CH_EVE.BEHAVE_time_100 = EPHYS.CH_EVE_.BEHAVE_time_100(:);
EPHYS.CH_EVE.BEHAVE_time_1K = EPHYS.CH_EVE_.BEHAVE_time_1K(:);

fprintf(' --> Completed. \n')

%% load DLC DATA
if ~params.auto
    %     file_name = EPHYS.CH_sorted_file_name(1:13);
    [file_name,file_path] = uigetfile([file_path filesep '_DLC.mat'], 'Select _DLC file');
else
    file_name = [EPHYS.CH_sorted_file_name(1:13) '_DLC.mat'];
    file_path = params.file_path;
end
fprintf(['Loading ', file_name, ' ... ']);
VID = load([file_path file_name]);
fprintf(' --> Completed. \n')

%% load BEHAVE DATA
if ~params.auto
    %     file_name = EPHYS.CH_sorted_file_name(1:13);
    [file_name,file_path] = uigetfile([file_path filesep '_ANALYZED.mat'], 'Select _ANALYZED file');
else
    file_name = [EPHYS.CH_sorted_file_name(1:13) '_ANALYZED.mat'];
    file_path = params.file_path;
end
fprintf(['Loading ', file_name, ' ... ']);
BEHAVE = load([file_path file_name]);
fprintf(' --> Completed. \n')

%% build EPHYS.CH_sorted from DATA_PSORT
ch_data = double(DATA_PSORT.topLevel_data.ch_data);
ch_time = double(DATA_PSORT.topLevel_data.ch_time);
SS_index = find(logical(double(DATA_PSORT.topLevel_data.ss_index)));
CS_index = find(logical(double(DATA_PSORT.topLevel_data.cs_index)));
SS_time = ch_time(SS_index);
CS_time = ch_time(CS_index);

waveform_inds_span = ((-60+1) : 1 : (120));
SS_inds = repmat(waveform_inds_span(:)', length(SS_index), 1) + repmat(SS_index(:), 1, length(waveform_inds_span));
SS_inds(SS_inds < 1) = 1;
SS_inds(SS_inds > length(ch_data)) = length(ch_data);
CS_inds = repmat(waveform_inds_span(:)', length(CS_index), 1) + repmat(CS_index(:), 1, length(waveform_inds_span));
CS_inds(CS_inds < 1) = 1;
CS_inds(CS_inds > length(ch_data)) = length(ch_data);
SS_waveform = ch_data(SS_inds);
CS_waveform = ch_data(CS_inds);

if length(SS_index) == 1
    SS_waveform = SS_waveform(:)';
end

if length(CS_index) == 1
    CS_waveform = CS_waveform(:)';
end

EPHYS.CH_sorted.SS_data.SS_ind = SS_index;
EPHYS.CH_sorted.CS_data.CS_ind = CS_index;
EPHYS.CH_sorted.SS_data.SS_time = SS_time;
EPHYS.CH_sorted.CS_data.CS_time = CS_time;
EPHYS.CH_sorted.SS_data.SS_waveform = SS_waveform;
EPHYS.CH_sorted.CS_data.CS_waveform = CS_waveform;
EPHYS.signal.ch_data = ch_data;
EPHYS.signal.ch_time = ch_time;

%% build SSxSS AUTO PROBABILITY
clearvars -except EPHYS BEHAVE VID
fprintf(['Building SSxSS_AUTO & CSxSS_AUTO PROBABILITY ' ' ...'])
SS_time   = EPHYS.CH_sorted.SS_data.SS_time;
CS_time   = EPHYS.CH_sorted.CS_data.CS_time;
Corr_data = ESN_correlogram(SS_time, CS_time);
EPHYS.CH_sorted.Corr_data.CS_inds_span     = Corr_data.CS_inds_span;
EPHYS.CH_sorted.Corr_data.CS_bin_size_time = Corr_data.CS_bin_size_time;
EPHYS.CH_sorted.Corr_data.SS_inds_span     = Corr_data.SS_inds_span;
EPHYS.CH_sorted.Corr_data.SS_bin_size_time = Corr_data.SS_bin_size_time;
EPHYS.CH_sorted.Corr_data.SS_SSxSS_AUTO    = Corr_data.SS_SSxSS_AUTO;
EPHYS.CH_sorted.Corr_data.CS_CSxSS_AUTO    = Corr_data.CS_CSxSS_AUTO;

fprintf(' --> Completed. \n')

%% SS & CS train_aligned
clearvars -except EPHYS BEHAVE VID
fprintf(['Building CS & SS train_aligned', ' ... ']);
EPHYS_time_1K     = EPHYS.CH_EVE.EPHYS_time_1K;
length_time_ = length(EPHYS_time_1K);
CS_time = EPHYS.CH_sorted.CS_data.CS_time;
if isempty(CS_time)
    CS_time = EPHYS_time_1K(1);
end
CS_time(end+1) = max([EPHYS_time_1K(end), CS_time(end)])+1;
SS_time = EPHYS.CH_sorted.SS_data.SS_time;
if isempty(SS_time)
    SS_time = EPHYS_time_1K(1);
end
SS_time(end+1) = max([EPHYS_time_1K(end), SS_time(end)])+1;
EPHYS_CS_train_1K = false(size(EPHYS_time_1K));
EPHYS_SS_train_1K = false(size(EPHYS_time_1K));
counter_CS = find(CS_time >= EPHYS_time_1K(1), 1, 'first');
counter_SS = find(SS_time >= EPHYS_time_1K(1), 1, 'first');
for counter_time_point = 1 : length_time_
    time_ponit_ = EPHYS_time_1K(counter_time_point);
    if time_ponit_>=CS_time(counter_CS)
        EPHYS_CS_train_1K(counter_time_point) = true;
        counter_CS = counter_CS + 1;
    end
    if time_ponit_>=SS_time(counter_SS)
        EPHYS_SS_train_1K(counter_time_point) = true;
        counter_SS = counter_SS + 1;
    end
end
EPHYS.CH_EVE.EPHYS_CS_train_1K = EPHYS_CS_train_1K;
EPHYS.CH_EVE.EPHYS_SS_train_1K = EPHYS_SS_train_1K;

EPHYS_time_100  = EPHYS.CH_EVE.EPHYS_time_100  ;

EPHYS_CS_train_100 = nan(length(EPHYS_time_100), 1);
EPHYS_SS_train_100 = nan(length(EPHYS_time_100), 1);

inds_str = 1:10:(length(EPHYS_CS_train_1K) - rem(length(EPHYS_CS_train_1K),10) - 10);
inds_end = inds_str + 9;

for counter_time_100 = 1 : length(inds_str)
    inds_ = inds_str(counter_time_100) : inds_end(counter_time_100);
    EPHYS_CS_train_100(counter_time_100,1) = nansum(EPHYS_CS_train_1K(inds_));
    EPHYS_SS_train_100(counter_time_100,1) = nansum(EPHYS_SS_train_1K(inds_));
    %disp(num2str(counter_time_100))
end

EPHYS.CH_EVE.EPHYS_CS_train_100 = EPHYS_CS_train_100;
EPHYS.CH_EVE.EPHYS_SS_train_100 = EPHYS_SS_train_100;


fprintf(' --> Completed. \n')

%% Build kinematics
clearvars -except EPHYS BEHAVE VID
fprintf(['Building d_tip, v_tip, & angle_tip', ' ... ']);

VID_time_vid = VID.DLC.TIME.time_vid';
VID_time_100 = EPHYS.CH_EVE.VID_time_100;


d_tip_vid = VID.DLC.KINEMATIC.d_tip;
v_tip_vid = ([0 (diff(d_tip_vid)./ VID_time_vid(1))'])';
angle_tip_vid = VID.DLC.KINEMATIC.angle_midtip;

d_tip_100 = interp1(VID_time_vid,  d_tip_vid, VID_time_100);
v_tip_100 = interp1(VID_time_vid,  v_tip_vid, VID_time_100);
angle_tip_100 = interp1(VID_time_vid,  angle_tip_vid, VID_time_100);


EPHYS.CH_EVE.VID_d_tip_100 = d_tip_100;
EPHYS.CH_EVE.VID_v_tip_100 = v_tip_100;
EPHYS.CH_EVE.VID_angle_tip_100 = angle_tip_100;


fprintf(' --> Completed. \n')

%% Build BEHAVE_eye_r_vm_filt
% clearvars -except EPHYS BEHAVE VID
% fprintf(['Building BEHAVE_eye_r_vm_filt', ' ... ']);
% eye_r_vm_filt = cell2mat(BEHAVE.TRIALS_DATA.eye_r_vm_filt(:));
% time_1K_cell2mat = cell2mat(BEHAVE.TRIALS_DATA.time_1K(:));
% BEHAVE_time_1K     = EPHYS.CH_EVE.BEHAVE_time_1K;
% BEHAVE_time_100     = EPHYS.CH_EVE.BEHAVE_time_100;
% length_time_ = length(BEHAVE_time_1K);
% BEHAVE_eye_r_vm_filt_1K = nan(size(BEHAVE_time_1K));
% time_1K_cell2mat(end+1)    = max([BEHAVE_time_1K(end), time_1K_cell2mat(end)])+1;
% counter_time_1K_cell2mat = find(time_1K_cell2mat >= BEHAVE_time_1K(1), 1, 'first');
% for counter_time_point = 1 : length_time_
%     time_ponit_ = BEHAVE_time_1K(counter_time_point);
%     if time_ponit_>=time_1K_cell2mat(counter_time_1K_cell2mat)
%         BEHAVE_eye_r_vm_filt_1K(counter_time_point) = eye_r_vm_filt(counter_time_1K_cell2mat);
%         counter_time_1K_cell2mat = counter_time_1K_cell2mat + 1;
%     end
% end
% EPHYS.CH_EVE.BEHAVE_eye_r_vm_filt_1K = BEHAVE_eye_r_vm_filt_1K;
% BEHAVE_eye_r_vm_filt_100 = (interp1(BEHAVE_time_1K, BEHAVE_eye_r_vm_filt_1K, BEHAVE_time_100));
%
% % if length(BEHAVE_eye_r_vm_filt_100) < length(EPHYS.CH_EVE.VID_time_100)
% %     BEHAVE_eye_r_vm_filt_100 = [BEHAVE_eye_r_vm_filt_100; nan((length(EPHYS.CH_EVE.VID_time_100) - length(BEHAVE_eye_r_vm_filt_100)),1)];
% % end
%
% EPHYS.CH_EVE.BEHAVE_eye_r_vm_filt_100 = BEHAVE_eye_r_vm_filt_100;
% fprintf(' --> Completed. \n')

%% Build Raster Data (inds_span)
clearvars -except EPHYS VID BEHAVE
fprintf(['Building Raster Plot Data', ' ... ']);

LB =30 ; UB = 30;

VID_LED_xcorr_time_100  = EPHYS.CH_EVE.align_LED.VID_LED_xcorr_time_100  ;
length_time_ = length(VID_LED_xcorr_time_100);
VID_LED_xcorr_time_lick_onset =  VID.DLC.TIME.time_lick_onset;
VID_LED_xcorr_time_lick_vmax = VID.DLC.TIME.time_v_lick_max_abs;
VID_LED_xcorr_time_lick_dmax =  VID.DLC.TIME.time_d_lick_max_abs;
VID_LED_xcorr_time_lick_vmin =  VID.DLC.TIME.time_v_lick_min_abs;
VID_LED_xcorr_time_lick_offset =  VID.DLC.TIME.time_lick_offset;

num_licks = nansum(VID_LED_xcorr_time_lick_offset < VID_LED_xcorr_time_100(end));

VID_LED_xcorr_time_lick_onset = VID_LED_xcorr_time_lick_onset(1:num_licks);
VID_LED_xcorr_time_lick_vmax = VID_LED_xcorr_time_lick_vmax(1:num_licks);
VID_LED_xcorr_time_lick_dmax = VID_LED_xcorr_time_lick_dmax(1:num_licks);
VID_LED_xcorr_time_lick_vmin = VID_LED_xcorr_time_lick_vmin(1:num_licks);
VID_LED_xcorr_time_lick_offset = VID_LED_xcorr_time_lick_offset(1:num_licks);

VID_LED_xcorr_time_lick_onset = VID_LED_xcorr_time_lick_onset(1:num_licks);
VID_LED_xcorr_time_lick_vmax = VID_LED_xcorr_time_lick_vmax(1:num_licks);
VID_LED_xcorr_time_lick_dmax = VID_LED_xcorr_time_lick_dmax(1:num_licks);
VID_LED_xcorr_time_lick_vmin = VID_LED_xcorr_time_lick_vmin(1:num_licks);
VID_LED_xcorr_time_lick_offset = VID_LED_xcorr_time_lick_offset(1:num_licks);

VID_LED_xcorr_time_lick_onset(end+1)    = max([VID_LED_xcorr_time_100(end), VID_LED_xcorr_time_lick_onset(end)])+1;
VID_LED_xcorr_time_lick_vmax(end+1)    = max([VID_LED_xcorr_time_100(end), VID_LED_xcorr_time_lick_vmax(end)])+1;
VID_LED_xcorr_time_lick_dmax(end+1)    = max([VID_LED_xcorr_time_100(end), VID_LED_xcorr_time_lick_dmax(end)])+1;
VID_LED_xcorr_time_lick_vmin(end+1)    = max([VID_LED_xcorr_time_100(end), VID_LED_xcorr_time_lick_vmin(end)])+1;
VID_LED_xcorr_time_lick_offset(end+1)    = max([VID_LED_xcorr_time_100(end), VID_LED_xcorr_time_lick_offset(end)])+1;

VID_LED_xcorr_ind_lick_onset = nan(num_licks, 1);
VID_LED_xcorr_ind_lick_vmax   = nan(num_licks, 1);
VID_LED_xcorr_ind_lick_dmax = nan(num_licks, 1);
VID_LED_xcorr_ind_lick_vmin = nan(num_licks, 1);
VID_LED_xcorr_ind_lick_offset = nan(num_licks, 1);

counter_lick_onset   = find(VID_LED_xcorr_time_lick_onset    >= VID_LED_xcorr_time_100(1), 1, 'first');
counter_lick_vmax  = find(VID_LED_xcorr_time_lick_vmax  >= VID_LED_xcorr_time_100(1), 1, 'first');
counter_lick_dmax   = find(VID_LED_xcorr_time_lick_dmax   >= VID_LED_xcorr_time_100(1), 1, 'first');
counter_lick_vmin = find(VID_LED_xcorr_time_lick_vmin >= VID_LED_xcorr_time_100(1), 1, 'first');
counter_lick_offset  = find(VID_LED_xcorr_time_lick_offset  >= VID_LED_xcorr_time_100(1), 1, 'first');

for counter_time_point = 1:length_time_
    time_point_ = VID_LED_xcorr_time_100(counter_time_point);
    if time_point_ >= VID_LED_xcorr_time_lick_onset(counter_lick_onset)
        VID_LED_xcorr_ind_lick_onset(counter_lick_onset) = counter_time_point;
        counter_lick_onset = counter_lick_onset + 1;
    end
    if time_point_ >= VID_LED_xcorr_time_lick_vmax(counter_lick_vmax)
        VID_LED_xcorr_ind_lick_vmax(counter_lick_vmax) = counter_time_point;
        counter_lick_vmax = counter_lick_vmax + 1;
    end
    if time_point_ >= VID_LED_xcorr_time_lick_dmax(counter_lick_dmax)
        VID_LED_xcorr_ind_lick_dmax(counter_lick_dmax) = counter_time_point;
        counter_lick_dmax = counter_lick_dmax + 1;
    end
    if time_point_ >= VID_LED_xcorr_time_lick_vmin(counter_lick_vmin)
        VID_LED_xcorr_ind_lick_vmin(counter_lick_vmin) = counter_time_point;
        counter_lick_vmin = counter_lick_vmin + 1;
    end
    if time_point_ >= VID_LED_xcorr_time_lick_offset(counter_lick_offset)
        VID_LED_xcorr_ind_lick_offset(counter_lick_offset) = counter_time_point;
        counter_lick_offset = counter_lick_offset + 1;
    end
    
    
end

% convert xcorr to aligned for EPHYS. We find the events on VID and then
% should convert it to EPHYS for SS & CS (EPHYS related events). We should not convert VID for
% VID related events.

EPHYS_LED_aligned_ind_lick_onset_100    = EPHYS.CH_EVE.align_LED.EPHYS_LED_aligned_ind_100(VID_LED_xcorr_ind_lick_onset);
EPHYS_LED_aligned_ind_lick_vmax_100    = EPHYS.CH_EVE.align_LED.EPHYS_LED_aligned_ind_100(VID_LED_xcorr_ind_lick_vmax);
EPHYS_LED_aligned_ind_lick_dmax_100    = EPHYS.CH_EVE.align_LED.EPHYS_LED_aligned_ind_100(VID_LED_xcorr_ind_lick_dmax);
EPHYS_LED_aligned_ind_lick_vmin_100    = EPHYS.CH_EVE.align_LED.EPHYS_LED_aligned_ind_100(VID_LED_xcorr_ind_lick_vmin);
EPHYS_LED_aligned_ind_lick_offset_100    = EPHYS.CH_EVE.align_LED.EPHYS_LED_aligned_ind_100(VID_LED_xcorr_ind_lick_offset);

inds_span_lick_onset    = ((-LB+1) : 1 : (UB))';
inds_span_lick_vmax  = ((-LB+1) : 1 : (UB))';
inds_span_lick_dmax   = ((-LB+1) : 1 : (UB))';
inds_span_lick_vmin = ((-LB+1) : 1 : (UB))';
inds_span_lick_offset  = ((-LB+1) : 1 : (UB))';

% Build EPHYS_LED_aligned_inds
EPHYS_time_100     = EPHYS.CH_EVE.EPHYS_time_100;
length_time_      = length(EPHYS_time_100);

EPHYS_LED_aligned_inds_lick_onset_100 = repmat( EPHYS_LED_aligned_ind_lick_onset_100(:), 1, length(inds_span_lick_onset)) + repmat(inds_span_lick_onset(:)', length(VID_LED_xcorr_ind_lick_onset), 1);
EPHYS_LED_aligned_inds_lick_onset_100( EPHYS_LED_aligned_inds_lick_onset_100 < 1 ) = 1;
EPHYS_LED_aligned_inds_lick_onset_100( EPHYS_LED_aligned_inds_lick_onset_100 > length_time_ ) = length_time_;
EPHYS_LED_aligned_inds_lick_vmax_100 = repmat( EPHYS_LED_aligned_ind_lick_vmax_100(:), 1, length(inds_span_lick_vmax)) + repmat(inds_span_lick_vmax(:)', length(VID_LED_xcorr_ind_lick_vmax), 1);
EPHYS_LED_aligned_inds_lick_vmax_100( EPHYS_LED_aligned_inds_lick_vmax_100 < 1 ) = 1;
EPHYS_LED_aligned_inds_lick_vmax_100( EPHYS_LED_aligned_inds_lick_vmax_100 > length_time_ ) = length_time_;
EPHYS_LED_aligned_inds_lick_dmax_100 = repmat( EPHYS_LED_aligned_ind_lick_dmax_100(:), 1, length(inds_span_lick_dmax)) + repmat(inds_span_lick_dmax(:)', length(VID_LED_xcorr_ind_lick_dmax), 1);
EPHYS_LED_aligned_inds_lick_dmax_100( EPHYS_LED_aligned_inds_lick_dmax_100 < 1 ) = 1;
EPHYS_LED_aligned_inds_lick_dmax_100( EPHYS_LED_aligned_inds_lick_dmax_100 > length_time_ ) = length_time_;
EPHYS_LED_aligned_inds_lick_vmin_100 = repmat( EPHYS_LED_aligned_ind_lick_vmin_100(:), 1, length(inds_span_lick_vmin)) + repmat(inds_span_lick_vmin(:)', length(VID_LED_xcorr_ind_lick_vmin), 1);
EPHYS_LED_aligned_inds_lick_vmin_100( EPHYS_LED_aligned_inds_lick_vmin_100 < 1 ) = 1;
EPHYS_LED_aligned_inds_lick_vmin_100( EPHYS_LED_aligned_inds_lick_vmin_100 > length_time_ ) = length_time_;
EPHYS_LED_aligned_inds_lick_offset_100 = repmat( EPHYS_LED_aligned_ind_lick_offset_100(:), 1, length(inds_span_lick_offset)) + repmat(inds_span_lick_offset(:)', length(VID_LED_xcorr_ind_lick_offset), 1);
EPHYS_LED_aligned_inds_lick_offset_100( EPHYS_LED_aligned_inds_lick_offset_100 < 1 ) = 1;
EPHYS_LED_aligned_inds_lick_offset_100( EPHYS_LED_aligned_inds_lick_offset_100 > length_time_ ) = length_time_;

% We should not convert VID for VID related events.
VID_LED_aligned_ind_lick_onset_100    = EPHYS.CH_EVE.align_LED.VID_LED_xcorr_ind_100(VID_LED_xcorr_ind_lick_onset);
VID_LED_aligned_ind_lick_vmax_100    = EPHYS.CH_EVE.align_LED.VID_LED_xcorr_ind_100(VID_LED_xcorr_ind_lick_vmax);
VID_LED_aligned_ind_lick_dmax_100    = EPHYS.CH_EVE.align_LED.VID_LED_xcorr_ind_100(VID_LED_xcorr_ind_lick_dmax);
VID_LED_aligned_ind_lick_vmin_100    = EPHYS.CH_EVE.align_LED.VID_LED_xcorr_ind_100(VID_LED_xcorr_ind_lick_vmin);
VID_LED_aligned_ind_lick_offset_100    = EPHYS.CH_EVE.align_LED.VID_LED_xcorr_ind_100(VID_LED_xcorr_ind_lick_offset);

% Build VID_LED_aligned_inds
VID_time_100     = EPHYS.CH_EVE.VID_time_100;
length_time_      = length(VID_time_100);

VID_LED_aligned_inds_lick_onset_100 = repmat( VID_LED_aligned_ind_lick_onset_100(:), 1, length(inds_span_lick_onset)) + repmat(inds_span_lick_onset(:)', length(VID_LED_xcorr_ind_lick_onset), 1);
VID_LED_aligned_inds_lick_onset_100( VID_LED_aligned_inds_lick_onset_100 < 1 ) = 1;
VID_LED_aligned_inds_lick_onset_100( VID_LED_aligned_inds_lick_onset_100 > length_time_ ) = length_time_;
VID_LED_aligned_inds_lick_vmax_100 = repmat( VID_LED_aligned_ind_lick_vmax_100(:), 1, length(inds_span_lick_vmax)) + repmat(inds_span_lick_vmax(:)', length(VID_LED_xcorr_ind_lick_vmax), 1);
VID_LED_aligned_inds_lick_vmax_100( VID_LED_aligned_inds_lick_vmax_100 < 1 ) = 1;
VID_LED_aligned_inds_lick_vmax_100( VID_LED_aligned_inds_lick_vmax_100 > length_time_ ) = length_time_;
VID_LED_aligned_inds_lick_dmax_100 = repmat( VID_LED_aligned_ind_lick_dmax_100(:), 1, length(inds_span_lick_dmax)) + repmat(inds_span_lick_dmax(:)', length(VID_LED_xcorr_ind_lick_dmax), 1);
VID_LED_aligned_inds_lick_dmax_100( VID_LED_aligned_inds_lick_dmax_100 < 1 ) = 1;
VID_LED_aligned_inds_lick_dmax_100( VID_LED_aligned_inds_lick_dmax_100 > length_time_ ) = length_time_;
VID_LED_aligned_inds_lick_vmin_100 = repmat( VID_LED_aligned_ind_lick_vmin_100(:), 1, length(inds_span_lick_vmin)) + repmat(inds_span_lick_vmin(:)', length(VID_LED_xcorr_ind_lick_vmin), 1);
VID_LED_aligned_inds_lick_vmin_100( VID_LED_aligned_inds_lick_vmin_100 < 1 ) = 1;
VID_LED_aligned_inds_lick_vmin_100( VID_LED_aligned_inds_lick_vmin_100 > length_time_ ) = length_time_;
VID_LED_aligned_inds_lick_offset_100 = repmat( VID_LED_aligned_ind_lick_offset_100(:), 1, length(inds_span_lick_offset)) + repmat(inds_span_lick_offset(:)', length(VID_LED_xcorr_ind_lick_offset), 1);
VID_LED_aligned_inds_lick_offset_100( VID_LED_aligned_inds_lick_offset_100 < 1 ) = 1;
VID_LED_aligned_inds_lick_offset_100( VID_LED_aligned_inds_lick_offset_100 > length_time_ ) = length_time_;

EPHYS.CH_EVE.EPHYS_LED_aligned_ind_lick_onset_100=  EPHYS_LED_aligned_ind_lick_onset_100;
EPHYS.CH_EVE.EPHYS_LED_aligned_ind_lick_vmax_100 = EPHYS_LED_aligned_ind_lick_vmax_100;
EPHYS.CH_EVE.EPHYS_LED_aligned_ind_lick_dmax_100 = EPHYS_LED_aligned_ind_lick_dmax_100;
EPHYS.CH_EVE.EPHYS_LED_aligned_ind_lick_vmin_100 = EPHYS_LED_aligned_ind_lick_vmin_100;
EPHYS.CH_EVE.EPHYS_LED_aligned_ind_lick_offset_100 = EPHYS_LED_aligned_ind_lick_offset_100;
EPHYS.CH_EVE.EPHYS_LED_aligned_inds_lick_onset_100=  EPHYS_LED_aligned_inds_lick_onset_100;
EPHYS.CH_EVE.EPHYS_LED_aligned_inds_lick_vmax_100 = EPHYS_LED_aligned_inds_lick_vmax_100;
EPHYS.CH_EVE.EPHYS_LED_aligned_inds_lick_dmax_100 = EPHYS_LED_aligned_inds_lick_dmax_100;
EPHYS.CH_EVE.EPHYS_LED_aligned_inds_lick_vmin_100 = EPHYS_LED_aligned_inds_lick_vmin_100;
EPHYS.CH_EVE.EPHYS_LED_aligned_inds_lick_offset_100 = EPHYS_LED_aligned_inds_lick_offset_100;
EPHYS.CH_EVE.VID_LED_aligned_ind_lick_onset_100=  VID_LED_aligned_ind_lick_onset_100;
EPHYS.CH_EVE.VID_LED_aligned_ind_lick_vmax_100 = VID_LED_aligned_ind_lick_vmax_100;
EPHYS.CH_EVE.VID_LED_aligned_ind_lick_dmax_100 = VID_LED_aligned_ind_lick_dmax_100;
EPHYS.CH_EVE.VID_LED_aligned_ind_lick_vmin_100 = VID_LED_aligned_ind_lick_vmin_100;
EPHYS.CH_EVE.VID_LED_aligned_ind_lick_offset_100 = VID_LED_aligned_ind_lick_offset_100;
EPHYS.CH_EVE.VID_LED_aligned_inds_lick_onset_100=  VID_LED_aligned_inds_lick_onset_100;
EPHYS.CH_EVE.VID_LED_aligned_inds_lick_vmax_100 = VID_LED_aligned_inds_lick_vmax_100;
EPHYS.CH_EVE.VID_LED_aligned_inds_lick_dmax_100 = VID_LED_aligned_inds_lick_dmax_100;
EPHYS.CH_EVE.VID_LED_aligned_inds_lick_vmin_100 = VID_LED_aligned_inds_lick_vmin_100;
EPHYS.CH_EVE.VID_LED_aligned_inds_lick_offset_100 = VID_LED_aligned_inds_lick_offset_100;

EPHYS.CH_EVE.inds_span_lick_onset = inds_span_lick_onset(:)';
EPHYS.CH_EVE.inds_span_lick_vmax = inds_span_lick_vmax(:)';
EPHYS.CH_EVE.inds_span_lick_dmax = inds_span_lick_dmax(:)';
EPHYS.CH_EVE.inds_span_lick_vmin = inds_span_lick_vmin(:)';
EPHYS.CH_EVE.inds_span_lick_offset = inds_span_lick_offset(:)';

% EYE
num_trials = length(BEHAVE.TRIALS_DATA.time_end);
BEHAVE_EB_xcorr_time_100     = EPHYS.CH_EVE_.align_states.BEHAVE_EB_xcorr_time_100;
length_time_ = length(BEHAVE_EB_xcorr_time_100);
BEHAVE_EB_xcorr_time_cue_present    = nan(num_trials, 1);
BEHAVE_EB_xcorr_time_primSac_onset  = nan(num_trials, 1);
BEHAVE_EB_xcorr_time_primSac_vmax  = nan(num_trials, 1);
BEHAVE_EB_xcorr_time_primSac_offset = nan(num_trials, 1);
BEHAVE_EB_xcorr_time_corrSac_onset  = nan(num_trials, 1);
BEHAVE_EB_xcorr_time_corrSac_vmax  = nan(num_trials, 1);
BEHAVE_EB_xcorr_time_corrSac_offset  = nan(num_trials, 1);
for counter_trial = 1 : 1 : num_trials
    BEHAVE_EB_xcorr_time_cue_present(counter_trial) = BEHAVE.TRIALS_DATA.time_state_cue_present{1,counter_trial}(end);
    
    ind_primSac_onset_  = (BEHAVE.SACS_PRIM_DATA.ind_start(1,counter_trial)  == BEHAVE.SACS_PRIM_DATA.inds(:,counter_trial));
    if nansum(ind_primSac_onset_) ~= 1
        ind_primSac_onset_ = 1;
    end
    BEHAVE_EB_xcorr_time_primSac_onset(counter_trial)  = BEHAVE.SACS_PRIM_DATA.time(ind_primSac_onset_, counter_trial);
    
    ind_primSac_vmax_  = (BEHAVE.SACS_PRIM_DATA.ind_vmax(1,counter_trial)  == BEHAVE.SACS_PRIM_DATA.inds(:,counter_trial));
    if nansum(ind_primSac_vmax_) ~= 1
        ind_primSac_vmax_ = 60;
    end
    BEHAVE_EB_xcorr_time_primSac_vmax(counter_trial)  = BEHAVE.SACS_PRIM_DATA.time(ind_primSac_vmax_, counter_trial);
    
    ind_primSac_offset_ = (BEHAVE.SACS_PRIM_DATA.ind_finish(1,counter_trial) == BEHAVE.SACS_PRIM_DATA.inds(:,counter_trial));
    if nansum(ind_primSac_offset_) ~= 1
        ind_primSac_offset_ = 150;
    end
    BEHAVE_EB_xcorr_time_primSac_offset(counter_trial) = BEHAVE.SACS_PRIM_DATA.time(ind_primSac_offset_, counter_trial);
    
    ind_corrSac_onset_  = (BEHAVE.SACS_CORR_DATA.ind_start(1,counter_trial)  == BEHAVE.SACS_CORR_DATA.inds(:,counter_trial));
    if nansum(ind_corrSac_onset_) ~= 1
        ind_corrSac_onset_ = 1;
    end
    BEHAVE_EB_xcorr_time_corrSac_onset(counter_trial)  = BEHAVE.SACS_CORR_DATA.time(ind_corrSac_onset_, counter_trial);
    
    ind_corrSac_vmax_  = (BEHAVE.SACS_CORR_DATA.ind_vmax(1,counter_trial)  == BEHAVE.SACS_CORR_DATA.inds(:,counter_trial));
    if nansum(ind_corrSac_vmax_) ~= 1
        ind_corrSac_vmax_ = 60;
    end
    BEHAVE_EB_xcorr_time_corrSac_vmax(counter_trial)  = BEHAVE.SACS_CORR_DATA.time(ind_corrSac_vmax_, counter_trial);
    
    ind_corrSac_offset_ = (BEHAVE.SACS_CORR_DATA.ind_finish(1,counter_trial) == BEHAVE.SACS_CORR_DATA.inds(:,counter_trial));
    if nansum(ind_corrSac_offset_) ~= 1
        ind_corrSac_offset_ = 150;
    end
    BEHAVE_EB_xcorr_time_corrSac_offset(counter_trial) = BEHAVE.SACS_CORR_DATA.time(ind_corrSac_offset_, counter_trial);
end

% num_trials = nansum(BEHAVE_EB_xcorr_time_corrSac_offset < BEHAVE_EB_xcorr_time_100(end));
last_trial_num = find(BEHAVE_EB_xcorr_time_corrSac_offset < BEHAVE_EB_xcorr_time_100(end), 1, 'last');
first_trial_num = find(BEHAVE_EB_xcorr_time_cue_present > BEHAVE_EB_xcorr_time_100(1), 1, 'first');
range_trials = first_trial_num : last_trial_num;
num_trials = length(range_trials);

BEHAVE_EB_xcorr_time_cue_present    = BEHAVE_EB_xcorr_time_cue_present(   range_trials);
BEHAVE_EB_xcorr_time_primSac_onset  = BEHAVE_EB_xcorr_time_primSac_onset( range_trials);
BEHAVE_EB_xcorr_time_primSac_vmax   = BEHAVE_EB_xcorr_time_primSac_vmax(  range_trials);
BEHAVE_EB_xcorr_time_primSac_offset = BEHAVE_EB_xcorr_time_primSac_offset(range_trials);
BEHAVE_EB_xcorr_time_corrSac_onset  = BEHAVE_EB_xcorr_time_corrSac_onset( range_trials);
BEHAVE_EB_xcorr_time_corrSac_vmax   = BEHAVE_EB_xcorr_time_corrSac_vmax(  range_trials);
BEHAVE_EB_xcorr_time_corrSac_offset = BEHAVE_EB_xcorr_time_corrSac_offset(range_trials);

BEHAVE_EB_xcorr_time_cue_present(end+1)    = max([BEHAVE_EB_xcorr_time_100(end), BEHAVE_EB_xcorr_time_cue_present(end)])+1;
BEHAVE_EB_xcorr_time_primSac_onset(end+1)  = max([BEHAVE_EB_xcorr_time_100(end), BEHAVE_EB_xcorr_time_primSac_onset(end)])+1;
BEHAVE_EB_xcorr_time_primSac_vmax(end+1)   = max([BEHAVE_EB_xcorr_time_100(end), BEHAVE_EB_xcorr_time_primSac_vmax(end)])+1;
BEHAVE_EB_xcorr_time_primSac_offset(end+1) = max([BEHAVE_EB_xcorr_time_100(end), BEHAVE_EB_xcorr_time_primSac_offset(end)])+1;
BEHAVE_EB_xcorr_time_corrSac_onset(end+1)  = max([BEHAVE_EB_xcorr_time_100(end), BEHAVE_EB_xcorr_time_corrSac_onset(end)])+1;
BEHAVE_EB_xcorr_time_corrSac_vmax(end+1)   = max([BEHAVE_EB_xcorr_time_100(end), BEHAVE_EB_xcorr_time_corrSac_vmax(end)])+1;
BEHAVE_EB_xcorr_time_corrSac_offset(end+1) = max([BEHAVE_EB_xcorr_time_100(end), BEHAVE_EB_xcorr_time_corrSac_offset(end)])+1;
BEHAVE_EB_xcorr_ind_cue_present    = nan(num_trials, 1);
BEHAVE_EB_xcorr_ind_primSac_onset  = nan(num_trials, 1);
BEHAVE_EB_xcorr_ind_primSac_vmax   = nan(num_trials, 1);
BEHAVE_EB_xcorr_ind_primSac_offset = nan(num_trials, 1);
BEHAVE_EB_xcorr_ind_corrSac_onset  = nan(num_trials, 1);
BEHAVE_EB_xcorr_ind_corrSac_vmax   = nan(num_trials, 1);
BEHAVE_EB_xcorr_ind_corrSac_offset = nan(num_trials, 1);
counter_cue_present    = find(BEHAVE_EB_xcorr_time_cue_present    >= BEHAVE_EB_xcorr_time_100(1), 1, 'first');
counter_primSac_onset  = find(BEHAVE_EB_xcorr_time_primSac_onset  >= BEHAVE_EB_xcorr_time_100(1), 1, 'first');
counter_primSac_vmax   = find(BEHAVE_EB_xcorr_time_primSac_vmax   >= BEHAVE_EB_xcorr_time_100(1), 1, 'first');
counter_primSac_offset = find(BEHAVE_EB_xcorr_time_primSac_offset >= BEHAVE_EB_xcorr_time_100(1), 1, 'first');
counter_corrSac_onset  = find(BEHAVE_EB_xcorr_time_corrSac_onset  >= BEHAVE_EB_xcorr_time_100(1), 1, 'first');
counter_corrSac_vmax   = find(BEHAVE_EB_xcorr_time_corrSac_vmax   >= BEHAVE_EB_xcorr_time_100(1), 1, 'first');
counter_corrSac_offset = find(BEHAVE_EB_xcorr_time_corrSac_offset >= BEHAVE_EB_xcorr_time_100(1), 1, 'first');
counter_trial_cue_present    = 1;
counter_trial_primSac_onset  = 1;
counter_trial_primSac_vmax   = 1;
counter_trial_primSac_offset = 1;
counter_trial_corrSac_onset  = 1;
counter_trial_corrSac_vmax   = 1;
counter_trial_corrSac_offset = 1;

% BEHAVE_EB_xcorr_ind_cue_present(isnan(BEHAVE_EB_xcorr_ind_cue_present)) = [];
% BEHAVE_EB_xcorr_ind_primSac_onset(isnan(BEHAVE_EB_xcorr_ind_primSac_onset)) = [];
% BEHAVE_EB_xcorr_ind_primSac_vmax(isnan(BEHAVE_EB_xcorr_ind_primSac_vmax)) = [];
% BEHAVE_EB_xcorr_ind_primSac_offset(isnan(BEHAVE_EB_xcorr_ind_primSac_offset)) = [];
% BEHAVE_EB_xcorr_ind_corrSac_onset(isnan(BEHAVE_EB_xcorr_ind_corrSac_onset)) = [];
% BEHAVE_EB_xcorr_ind_corrSac_vmax(isnan(BEHAVE_EB_xcorr_ind_corrSac_vmax)) = [];
% BEHAVE_EB_xcorr_ind_corrSac_offset(isnan(BEHAVE_EB_xcorr_ind_corrSac_offset)) = [];

for counter_time_point = 1 : length_time_
    time_ponit_ = BEHAVE_EB_xcorr_time_100(counter_time_point);
    if time_ponit_>=BEHAVE_EB_xcorr_time_cue_present(counter_cue_present)
        BEHAVE_EB_xcorr_ind_cue_present(counter_trial_cue_present) = counter_time_point;
        counter_cue_present = counter_cue_present + 1;
        counter_trial_cue_present = counter_trial_cue_present + 1;
    end
    if time_ponit_>=BEHAVE_EB_xcorr_time_primSac_onset(counter_primSac_onset)
        BEHAVE_EB_xcorr_ind_primSac_onset(counter_trial_primSac_onset) = counter_time_point;
        counter_primSac_onset = counter_primSac_onset + 1;
        counter_trial_primSac_onset = counter_trial_primSac_onset + 1;
    end
    if time_ponit_>=BEHAVE_EB_xcorr_time_primSac_vmax(counter_primSac_vmax)
        BEHAVE_EB_xcorr_ind_primSac_vmax(counter_trial_primSac_vmax) = counter_time_point;
        counter_primSac_vmax = counter_primSac_vmax + 1;
        counter_trial_primSac_vmax = counter_trial_primSac_vmax + 1;
    end
    if time_ponit_>=BEHAVE_EB_xcorr_time_primSac_offset(counter_primSac_offset)
        BEHAVE_EB_xcorr_ind_primSac_offset(counter_trial_primSac_offset) = counter_time_point;
        counter_primSac_offset = counter_primSac_offset + 1;
        counter_trial_primSac_offset = counter_trial_primSac_offset + 1;
    end
    if time_ponit_>=BEHAVE_EB_xcorr_time_corrSac_onset(counter_corrSac_onset)
        BEHAVE_EB_xcorr_ind_corrSac_onset(counter_trial_corrSac_onset) = counter_time_point;
        counter_corrSac_onset = counter_corrSac_onset + 1;
        counter_trial_corrSac_onset = counter_trial_corrSac_onset + 1;
    end
    if time_ponit_>=BEHAVE_EB_xcorr_time_corrSac_vmax(counter_corrSac_vmax)
        BEHAVE_EB_xcorr_ind_corrSac_vmax(counter_trial_corrSac_vmax) = counter_time_point;
        counter_corrSac_vmax = counter_corrSac_vmax + 1;
        counter_trial_corrSac_vmax = counter_trial_corrSac_vmax + 1;
    end
    if time_ponit_>=BEHAVE_EB_xcorr_time_corrSac_offset(counter_corrSac_offset)
        BEHAVE_EB_xcorr_ind_corrSac_offset(counter_trial_corrSac_offset) = counter_time_point;
        counter_corrSac_offset = counter_corrSac_offset + 1;
        counter_trial_corrSac_offset = counter_trial_corrSac_offset + 1;
    end
end
% convert xcorr to aligned for EPHYS. We find the events on BEHAVE and then
% should convert it to EPHYS for SS & CS (EPHYS related events). We should not convert BEHAVE for
% BEHAVE related events.
EPHYS_EB_aligned_ind_cue_present_100    = EPHYS.CH_EVE_.align_states.EPHYS_EB_aligned_ind_100(BEHAVE_EB_xcorr_ind_cue_present);
EPHYS_EB_aligned_ind_primSac_onset_100  = EPHYS.CH_EVE_.align_states.EPHYS_EB_aligned_ind_100(BEHAVE_EB_xcorr_ind_primSac_onset);
EPHYS_EB_aligned_ind_primSac_vmax_100  = EPHYS.CH_EVE_.align_states.EPHYS_EB_aligned_ind_100(BEHAVE_EB_xcorr_ind_primSac_vmax);
EPHYS_EB_aligned_ind_primSac_offset_100 = EPHYS.CH_EVE_.align_states.EPHYS_EB_aligned_ind_100(BEHAVE_EB_xcorr_ind_primSac_offset);
EPHYS_EB_aligned_ind_corrSac_onset_100  = EPHYS.CH_EVE_.align_states.EPHYS_EB_aligned_ind_100(BEHAVE_EB_xcorr_ind_corrSac_onset);
EPHYS_EB_aligned_ind_corrSac_vmax_100  = EPHYS.CH_EVE_.align_states.EPHYS_EB_aligned_ind_100(BEHAVE_EB_xcorr_ind_corrSac_vmax);
EPHYS_EB_aligned_ind_corrSac_offset_100 = EPHYS.CH_EVE_.align_states.EPHYS_EB_aligned_ind_100(BEHAVE_EB_xcorr_ind_corrSac_offset);
inds_span_cue_present    = ((-LB+1) : 1 : (UB))';
inds_span_primSac_onset  = ((-LB+1) : 1 : (UB))';
inds_span_primSac_vmax   =((-LB+1) : 1 : (UB))';
inds_span_primSac_offset = ((-LB+1) : 1 : (UB))';
inds_span_corrSac_onset  = ((-LB+1) : 1 : (UB))';
inds_span_corrSac_vmax   = ((-LB+1) : 1 : (UB))';
inds_span_corrSac_offset = ((-LB+1) : 1 : (UB))';
% Build EPHYS_EB_aligned_inds
EPHYS_time_100     = EPHYS.CH_EVE.EPHYS_time_100;
length_time_      = length(EPHYS_time_100);

EPHYS_EB_aligned_inds_cue_present_100 = repmat( EPHYS_EB_aligned_ind_cue_present_100(:), 1, length(inds_span_cue_present)) + repmat(inds_span_cue_present(:)', length(BEHAVE_EB_xcorr_ind_cue_present), 1);
EPHYS_EB_aligned_inds_cue_present_100( EPHYS_EB_aligned_inds_cue_present_100 < 1 ) = 1;
EPHYS_EB_aligned_inds_cue_present_100( EPHYS_EB_aligned_inds_cue_present_100 > length_time_ ) = length_time_;

EPHYS_EB_aligned_inds_primSac_onset_100 = repmat( EPHYS_EB_aligned_ind_primSac_onset_100(:), 1, length(inds_span_primSac_onset)) + repmat(inds_span_primSac_onset(:)', length(BEHAVE_EB_xcorr_ind_primSac_onset), 1);
EPHYS_EB_aligned_inds_primSac_onset_100( EPHYS_EB_aligned_inds_primSac_onset_100 < 1 ) = 1;
EPHYS_EB_aligned_inds_primSac_onset_100( EPHYS_EB_aligned_inds_primSac_onset_100 > length_time_ ) = length_time_;

EPHYS_EB_aligned_inds_primSac_vmax_100 = repmat( EPHYS_EB_aligned_ind_primSac_vmax_100(:), 1, length(inds_span_primSac_vmax)) + repmat(inds_span_primSac_vmax(:)', length(BEHAVE_EB_xcorr_ind_primSac_vmax), 1);
EPHYS_EB_aligned_inds_primSac_vmax_100( EPHYS_EB_aligned_inds_primSac_vmax_100 < 1 ) = 1;
EPHYS_EB_aligned_inds_primSac_vmax_100( EPHYS_EB_aligned_inds_primSac_vmax_100 > length_time_ ) = length_time_;

EPHYS_EB_aligned_inds_primSac_offset_100 = repmat( EPHYS_EB_aligned_ind_primSac_offset_100(:), 1, length(inds_span_primSac_offset)) + repmat(inds_span_primSac_offset(:)', length(BEHAVE_EB_xcorr_ind_primSac_offset), 1);
EPHYS_EB_aligned_inds_primSac_offset_100( EPHYS_EB_aligned_inds_primSac_offset_100 < 1 ) = 1;
EPHYS_EB_aligned_inds_primSac_offset_100( EPHYS_EB_aligned_inds_primSac_offset_100 > length_time_ ) = length_time_;

EPHYS_EB_aligned_inds_corrSac_onset_100 = repmat( EPHYS_EB_aligned_ind_corrSac_onset_100(:), 1, length(inds_span_corrSac_onset)) + repmat(inds_span_corrSac_onset(:)', length(BEHAVE_EB_xcorr_ind_corrSac_onset), 1);
EPHYS_EB_aligned_inds_corrSac_onset_100( EPHYS_EB_aligned_inds_corrSac_onset_100 < 1 ) = 1;
EPHYS_EB_aligned_inds_corrSac_onset_100( EPHYS_EB_aligned_inds_corrSac_onset_100 > length_time_ ) = length_time_;

EPHYS_EB_aligned_inds_corrSac_vmax_100 = repmat( EPHYS_EB_aligned_ind_corrSac_vmax_100(:), 1, length(inds_span_corrSac_vmax)) + repmat(inds_span_corrSac_vmax(:)', length(BEHAVE_EB_xcorr_ind_corrSac_vmax), 1);
EPHYS_EB_aligned_inds_corrSac_vmax_100( EPHYS_EB_aligned_inds_corrSac_vmax_100 < 1 ) = 1;
EPHYS_EB_aligned_inds_corrSac_vmax_100( EPHYS_EB_aligned_inds_corrSac_vmax_100 > length_time_ ) = length_time_;

EPHYS_EB_aligned_inds_corrSac_offset_100 = repmat( EPHYS_EB_aligned_ind_corrSac_offset_100(:), 1, length(inds_span_corrSac_offset)) + repmat(inds_span_corrSac_offset(:)', length(BEHAVE_EB_xcorr_ind_corrSac_offset), 1);
EPHYS_EB_aligned_inds_corrSac_offset_100( EPHYS_EB_aligned_inds_corrSac_offset_100 < 1 ) = 1;
EPHYS_EB_aligned_inds_corrSac_offset_100( EPHYS_EB_aligned_inds_corrSac_offset_100 > length_time_ ) = length_time_;

% We should not convert BEHAVE for BEHAVE related events.
BEHAVE_EB_aligned_ind_cue_present_100    = EPHYS.CH_EVE_.align_states.BEHAVE_EB_xcorr_ind_100(BEHAVE_EB_xcorr_ind_cue_present);
BEHAVE_EB_aligned_ind_primSac_onset_100  = EPHYS.CH_EVE_.align_states.BEHAVE_EB_xcorr_ind_100(BEHAVE_EB_xcorr_ind_primSac_onset);
BEHAVE_EB_aligned_ind_primSac_vmax_100  = EPHYS.CH_EVE_.align_states.BEHAVE_EB_xcorr_ind_100(BEHAVE_EB_xcorr_ind_primSac_vmax);
BEHAVE_EB_aligned_ind_primSac_offset_100 = EPHYS.CH_EVE_.align_states.BEHAVE_EB_xcorr_ind_100(BEHAVE_EB_xcorr_ind_primSac_offset);
BEHAVE_EB_aligned_ind_corrSac_onset_100  = EPHYS.CH_EVE_.align_states.BEHAVE_EB_xcorr_ind_100(BEHAVE_EB_xcorr_ind_corrSac_onset);
BEHAVE_EB_aligned_ind_corrSac_vmax_100  = EPHYS.CH_EVE_.align_states.BEHAVE_EB_xcorr_ind_100(BEHAVE_EB_xcorr_ind_corrSac_vmax);
BEHAVE_EB_aligned_ind_corrSac_offset_100  = EPHYS.CH_EVE_.align_states.BEHAVE_EB_xcorr_ind_100(BEHAVE_EB_xcorr_ind_corrSac_offset);
% Build BEHAVE_EB_aligned_inds
BEHAVE_time_100     = EPHYS.CH_EVE.BEHAVE_time_100;
length_time_      = length(BEHAVE_time_100);

BEHAVE_EB_aligned_inds_cue_present_100 = repmat( BEHAVE_EB_aligned_ind_cue_present_100(:), 1, length(inds_span_cue_present)) + repmat(inds_span_cue_present(:)', length(BEHAVE_EB_xcorr_ind_cue_present), 1);
BEHAVE_EB_aligned_inds_cue_present_100( BEHAVE_EB_aligned_inds_cue_present_100 < 1 ) = 1;
BEHAVE_EB_aligned_inds_cue_present_100( BEHAVE_EB_aligned_inds_cue_present_100 > length_time_ ) = length_time_;

BEHAVE_EB_aligned_inds_primSac_onset_100 = repmat( BEHAVE_EB_aligned_ind_primSac_onset_100(:), 1, length(inds_span_primSac_onset)) + repmat(inds_span_primSac_onset(:)', length(BEHAVE_EB_xcorr_ind_primSac_onset), 1);
BEHAVE_EB_aligned_inds_primSac_onset_100( BEHAVE_EB_aligned_inds_primSac_onset_100 < 1 ) = 1;
BEHAVE_EB_aligned_inds_primSac_onset_100( BEHAVE_EB_aligned_inds_primSac_onset_100 > length_time_ ) = length_time_;

BEHAVE_EB_aligned_inds_primSac_vmax_100 = repmat( BEHAVE_EB_aligned_ind_primSac_vmax_100(:), 1, length(inds_span_primSac_vmax)) + repmat(inds_span_primSac_vmax(:)', length(BEHAVE_EB_xcorr_ind_primSac_vmax), 1);
BEHAVE_EB_aligned_inds_primSac_vmax_100( BEHAVE_EB_aligned_inds_primSac_vmax_100 < 1 ) = 1;
BEHAVE_EB_aligned_inds_primSac_vmax_100( BEHAVE_EB_aligned_inds_primSac_vmax_100 > length_time_ ) = length_time_;

BEHAVE_EB_aligned_inds_primSac_offset_100 = repmat( BEHAVE_EB_aligned_ind_primSac_offset_100(:), 1, length(inds_span_primSac_offset)) + repmat(inds_span_primSac_offset(:)', length(BEHAVE_EB_xcorr_ind_primSac_offset), 1);
BEHAVE_EB_aligned_inds_primSac_offset_100( BEHAVE_EB_aligned_inds_primSac_offset_100 < 1 ) = 1;
BEHAVE_EB_aligned_inds_primSac_offset_100( BEHAVE_EB_aligned_inds_primSac_offset_100 > length_time_ ) = length_time_;

BEHAVE_EB_aligned_inds_corrSac_onset_100 = repmat( BEHAVE_EB_aligned_ind_corrSac_onset_100(:), 1, length(inds_span_corrSac_onset)) + repmat(inds_span_corrSac_onset(:)', length(BEHAVE_EB_xcorr_ind_corrSac_onset), 1);
BEHAVE_EB_aligned_inds_corrSac_onset_100( BEHAVE_EB_aligned_inds_corrSac_onset_100 < 1 ) = 1;
BEHAVE_EB_aligned_inds_corrSac_onset_100( BEHAVE_EB_aligned_inds_corrSac_onset_100 > length_time_ ) = length_time_;

BEHAVE_EB_aligned_inds_corrSac_vmax_100 = repmat( BEHAVE_EB_aligned_ind_corrSac_vmax_100(:), 1, length(inds_span_corrSac_vmax)) + repmat(inds_span_corrSac_vmax(:)', length(BEHAVE_EB_xcorr_ind_corrSac_vmax), 1);
BEHAVE_EB_aligned_inds_corrSac_vmax_100( BEHAVE_EB_aligned_inds_corrSac_vmax_100 < 1 ) = 1;
BEHAVE_EB_aligned_inds_corrSac_vmax_100( BEHAVE_EB_aligned_inds_corrSac_vmax_100 > length_time_ ) = length_time_;

BEHAVE_EB_aligned_inds_corrSac_offset_100 = repmat( BEHAVE_EB_aligned_ind_corrSac_offset_100(:), 1, length(inds_span_corrSac_offset)) + repmat(inds_span_corrSac_offset(:)', length(BEHAVE_EB_xcorr_ind_corrSac_offset), 1);
BEHAVE_EB_aligned_inds_corrSac_offset_100( BEHAVE_EB_aligned_inds_corrSac_offset_100 < 1 ) = 1;
BEHAVE_EB_aligned_inds_corrSac_offset_100( BEHAVE_EB_aligned_inds_corrSac_offset_100 > length_time_ ) = length_time_;

EPHYS.CH_EVE.EPHYS_EB_aligned_ind_cue_present_100     = EPHYS_EB_aligned_ind_cue_present_100;
EPHYS.CH_EVE.EPHYS_EB_aligned_ind_primSac_onset_100   = EPHYS_EB_aligned_ind_primSac_onset_100;
EPHYS.CH_EVE.EPHYS_EB_aligned_ind_primSac_vmax_100    = EPHYS_EB_aligned_ind_primSac_vmax_100;
EPHYS.CH_EVE.EPHYS_EB_aligned_ind_primSac_offset_100  = EPHYS_EB_aligned_ind_primSac_offset_100;
EPHYS.CH_EVE.EPHYS_EB_aligned_ind_corrSac_onset_100   = EPHYS_EB_aligned_ind_corrSac_onset_100;
EPHYS.CH_EVE.EPHYS_EB_aligned_ind_corrSac_vmax_100    = EPHYS_EB_aligned_ind_corrSac_vmax_100;
EPHYS.CH_EVE.EPHYS_EB_aligned_ind_corrSac_offset_100  = EPHYS_EB_aligned_ind_corrSac_offset_100;
EPHYS.CH_EVE.EPHYS_EB_aligned_inds_cue_present_100    = EPHYS_EB_aligned_inds_cue_present_100;
EPHYS.CH_EVE.EPHYS_EB_aligned_inds_primSac_onset_100  = EPHYS_EB_aligned_inds_primSac_onset_100;
EPHYS.CH_EVE.EPHYS_EB_aligned_inds_primSac_vmax_100   = EPHYS_EB_aligned_inds_primSac_vmax_100;
EPHYS.CH_EVE.EPHYS_EB_aligned_inds_primSac_offset_100 = EPHYS_EB_aligned_inds_primSac_offset_100;
EPHYS.CH_EVE.EPHYS_EB_aligned_inds_corrSac_onset_100  = EPHYS_EB_aligned_inds_corrSac_onset_100;
EPHYS.CH_EVE.EPHYS_EB_aligned_inds_corrSac_vmax_100   = EPHYS_EB_aligned_inds_corrSac_vmax_100;
EPHYS.CH_EVE.EPHYS_EB_aligned_inds_corrSac_offset_100 = EPHYS_EB_aligned_inds_corrSac_offset_100;
EPHYS.CH_EVE.BEHAVE_EB_aligned_ind_cue_present_100     = BEHAVE_EB_aligned_ind_cue_present_100;
EPHYS.CH_EVE.BEHAVE_EB_aligned_ind_primSac_onset_100   = BEHAVE_EB_aligned_ind_primSac_onset_100;
EPHYS.CH_EVE.BEHAVE_EB_aligned_ind_primSac_vmax_100    = BEHAVE_EB_aligned_ind_primSac_vmax_100;
EPHYS.CH_EVE.BEHAVE_EB_aligned_ind_primSac_offset_100  = BEHAVE_EB_aligned_ind_primSac_offset_100;
EPHYS.CH_EVE.BEHAVE_EB_aligned_ind_corrSac_onset_100   = BEHAVE_EB_aligned_ind_corrSac_onset_100;
EPHYS.CH_EVE.BEHAVE_EB_aligned_ind_corrSac_vmax_100    = BEHAVE_EB_aligned_ind_corrSac_vmax_100;
EPHYS.CH_EVE.BEHAVE_EB_aligned_ind_corrSac_offset_100  = BEHAVE_EB_aligned_ind_corrSac_offset_100;
EPHYS.CH_EVE.BEHAVE_EB_aligned_inds_cue_present_100    = BEHAVE_EB_aligned_inds_cue_present_100;
EPHYS.CH_EVE.BEHAVE_EB_aligned_inds_primSac_onset_100  = BEHAVE_EB_aligned_inds_primSac_onset_100;
EPHYS.CH_EVE.BEHAVE_EB_aligned_inds_primSac_vmax_100   = BEHAVE_EB_aligned_inds_primSac_vmax_100;
EPHYS.CH_EVE.BEHAVE_EB_aligned_inds_primSac_offset_100 = BEHAVE_EB_aligned_inds_primSac_offset_100;
EPHYS.CH_EVE.BEHAVE_EB_aligned_inds_corrSac_onset_100  = BEHAVE_EB_aligned_inds_corrSac_onset_100;
EPHYS.CH_EVE.BEHAVE_EB_aligned_inds_corrSac_vmax_100   = BEHAVE_EB_aligned_inds_corrSac_vmax_100;
EPHYS.CH_EVE.BEHAVE_EB_aligned_inds_corrSac_offset_100 = BEHAVE_EB_aligned_inds_corrSac_offset_100;
EPHYS.CH_EVE.inds_span_cue_present      = inds_span_cue_present(:)';
EPHYS.CH_EVE.inds_span_primSac_onset    = inds_span_primSac_onset(:)';
EPHYS.CH_EVE.inds_span_primSac_vmax     = inds_span_primSac_vmax(:)';
EPHYS.CH_EVE.inds_span_primSac_offset   = inds_span_primSac_offset(:)';
EPHYS.CH_EVE.inds_span_corrSac_onset    = inds_span_corrSac_onset(:)';
EPHYS.CH_EVE.inds_span_corrSac_vmax     = inds_span_corrSac_vmax(:)';
EPHYS.CH_EVE.inds_span_corrSac_offset   = inds_span_corrSac_offset(:)';

fprintf(' --> Completed. \n')

end

%% function concatenate_dataset
function upper_field_struct = concatenate_dataset(dataset_, upper_field_name, horz_OR_vert)
if ~isempty(upper_field_name)
    dataset = struct(upper_field_name,struct());
    field_names_ = fieldnames(dataset_(1).(upper_field_name));
    for counter_fields = 1 : 1 : length(field_names_)
        for counter_dataset = 1 : 1 : length(dataset_)
            variable_LICKS_ALL_ = dataset_(counter_dataset).(upper_field_name).(field_names_{counter_fields});
            % the field does not exist in LICKS_ALL
            if ~isfield(dataset.(upper_field_name), field_names_{counter_fields})
                dataset.(upper_field_name).(field_names_{counter_fields}) = [];
            end
            if(size(variable_LICKS_ALL_, 2) == 1) && (size(variable_LICKS_ALL_,1) == 60 || 100 || 200)
                variable_LICKS_ALL_ = variable_LICKS_ALL_';
            end
            variable_LICKS_DATA_ = dataset.(upper_field_name).(field_names_{counter_fields});
            variable_LICKS_DATA_ = horz_OR_vert(variable_LICKS_DATA_, variable_LICKS_ALL_);
            dataset.(upper_field_name).(field_names_{counter_fields}) = variable_LICKS_DATA_;
            %                    [num2str(counter_dataset) ', ' num2str(counter_fields)]
        end
    end
    upper_field_struct = dataset.(upper_field_name);
elseif isempty(upper_field_name)
    dataset = struct();
    field_names_ = fieldnames(dataset_);
    for counter_fields = 1 : 1 : length(field_names_)
        for counter_dataset = 1 : 1 : length(dataset_)
            variable_LICKS_ALL_ = dataset_(counter_dataset).(field_names_{counter_fields});
            % the field does not exist in TRIALS_DATA
            if ~isfield(dataset, field_names_{counter_fields})
                dataset.(field_names_{counter_fields}) = [];
            end
            variable_LICKS_DATA_ = dataset.(field_names_{counter_fields});
            variable_LICKS_DATA_ = horz_OR_vert(variable_LICKS_DATA_, variable_LICKS_ALL_);
            dataset.(field_names_{counter_fields}) = variable_LICKS_DATA_;
        end
    end
    upper_field_struct = dataset;
end
end

%% function plot_rasters_data_class
function fig_handle_ = plot_rasters_data_class(raster_data, plot_data, session_type)
%%
fig_num_               = plot_data.fig_num_;
xlabel_text_raster_    = plot_data.xlabel_text_raster_;
xlabel_text_raster_bout_   = plot_data.xlabel_text_raster_bout_;

xlabel_text_CS_probab_ = plot_data.xlabel_text_CS_probab_;
inds_span              = plot_data.inds_span * 10;

%bout: [-1 1]s
inds_span_LB           =(inds_span(1)-700) :10: inds_span(1)-10;
inds_span_UB           = (inds_span(end)+10) :10: (inds_span(end)+700);
inds_span_bout      = [inds_span_LB inds_span inds_span_UB]/1000;

% %corr: [0 1]s
% inds_span_UB           = (inds_span(end)+10) :10: (inds_span(end)+700);
% inds_span_corr = [inds_span(31:60) inds_span_UB]/1000;

%corr: [-1 1]s
inds_span_LB           =(inds_span(1)-700) :10: inds_span(1)-10;
inds_span_UB           = (inds_span(end)+10) :10: (inds_span(end)+700);
inds_span_corr      = [inds_span_LB inds_span inds_span_UB]/1000;

% All licks
train_data_logic_SS_all_corr = raster_data.train_data_logic_SS_all_corr;
train_data_logic_CS_all_corr = raster_data.train_data_logic_CS_all_corr;
train_data_logic_lick_onset_all_corr = raster_data.train_data_logic_lick_onset_all_corr;
train_data_logic_lick_vmax_all_corr = raster_data.train_data_logic_lick_vmax_all_corr;
train_data_logic_lick_dmax_all_corr = raster_data.train_data_logic_lick_dmax_all_corr;
train_data_logic_lick_vmin_all_corr = raster_data.train_data_logic_lick_vmin_all_corr;
train_data_logic_lick_offset_all_corr = raster_data.train_data_logic_lick_offset_all_corr;
train_data_logic_SS_all = raster_data.train_data_logic_SS_all ;
train_data_logic_CS_all  = raster_data.train_data_logic_CS_all ;
train_data_logic_lick_onset_all     = raster_data.train_data_logic_lick_onset_all ;
train_data_logic_lick_dmax_all    = raster_data.train_data_logic_lick_dmax_all ;
train_data_logic_lick_offset_all     = raster_data.train_data_logic_lick_offset_all ;

% Grooming licks
train_data_logic_SS_grooming_corr = raster_data.train_data_logic_SS_grooming_corr;
train_data_logic_CS_grooming_corr = raster_data.train_data_logic_CS_grooming_corr;
train_data_logic_lick_onset_grooming_corr = raster_data.train_data_logic_lick_onset_grooming_corr;
train_data_logic_lick_vmax_grooming_corr = raster_data.train_data_logic_lick_vmax_grooming_corr;
train_data_logic_lick_dmax_grooming_corr = raster_data.train_data_logic_lick_dmax_grooming_corr;
train_data_logic_lick_vmin_grooming_corr = raster_data.train_data_logic_lick_vmin_grooming_corr;
train_data_logic_lick_offset_grooming_corr = raster_data.train_data_logic_lick_offset_grooming_corr;
train_data_logic_SS_grooming = raster_data.train_data_logic_SS_grooming ;
% train_data_logic_SS_grooming_510= raster_data.train_data_logic_SS_grooming_510 ;
% train_data_logic_SS_grooming_10= raster_data.train_data_logic_SS_grooming_10 ;
train_data_logic_CS_grooming  = raster_data.train_data_logic_CS_grooming ;
train_data_logic_lick_onset_grooming     = raster_data.train_data_logic_lick_onset_grooming ;
train_data_logic_lick_vmax_grooming   = raster_data.train_data_logic_lick_vmax_grooming ;
train_data_logic_lick_dmax_grooming    = raster_data.train_data_logic_lick_dmax_grooming ;
train_data_logic_lick_vmin_grooming   = raster_data.train_data_logic_lick_vmin_grooming ;
train_data_logic_lick_offset_grooming     = raster_data.train_data_logic_lick_offset_grooming ;

% R licks
train_data_logic_SS_r_corr = raster_data.train_data_logic_SS_r_corr;
train_data_logic_CS_r_corr = raster_data.train_data_logic_CS_r_corr;
train_data_logic_lick_onset_r_corr = raster_data.train_data_logic_lick_onset_r_corr;
train_data_logic_lick_vmax_r_corr = raster_data.train_data_logic_lick_vmax_r_corr;
train_data_logic_lick_dmax_r_corr = raster_data.train_data_logic_lick_dmax_r_corr;
train_data_logic_lick_vmin_r_corr = raster_data.train_data_logic_lick_vmin_r_corr;
train_data_logic_lick_offset_r_corr = raster_data.train_data_logic_lick_offset_r_corr;
train_data_logic_SS_r = raster_data.train_data_logic_SS_r ;
% train_data_logic_SS_r_15 = raster_data.train_data_logic_SS_r_15 ;
% train_data_logic_SS_r_1520 = raster_data.train_data_logic_SS_r_1520 ;
train_data_logic_CS_r  = raster_data.train_data_logic_CS_r ;
train_data_logic_CS_r_1  = raster_data.train_data_logic_CS_r_1 ;
train_data_logic_CS_r_0  = raster_data.train_data_logic_CS_r_0 ;
train_data_logic_lick_onset_r     = raster_data.train_data_logic_lick_onset_r ;
train_data_logic_lick_vmax_r   = raster_data.train_data_logic_lick_vmax_r ;
train_data_logic_lick_dmax_r    = raster_data.train_data_logic_lick_dmax_r ;
train_data_logic_lick_vmin_r   = raster_data.train_data_logic_lick_vmin_r ;
train_data_logic_lick_offset_r     = raster_data.train_data_logic_lick_offset_r ;

% L licks
if (session_type  == 2)
    train_data_logic_SS_l_corr = raster_data.train_data_logic_SS_l_corr;
    train_data_logic_CS_l_corr = raster_data.train_data_logic_CS_l_corr;
    train_data_logic_lick_onset_l_corr = raster_data.train_data_logic_lick_onset_l_corr;
    train_data_logic_lick_vmax_l_corr = raster_data.train_data_logic_lick_vmax_l_corr;
    train_data_logic_lick_dmax_l_corr = raster_data.train_data_logic_lick_dmax_l_corr;
    train_data_logic_lick_vmin_l_corr = raster_data.train_data_logic_lick_vmin_l_corr;
    train_data_logic_lick_offset_l_corr = raster_data.train_data_logic_lick_offset_l_corr;
    train_data_logic_SS_l = raster_data.train_data_logic_SS_l ;
    %     train_data_logic_SS_l_15 = raster_data.train_data_logic_SS_l_15 ;
    %     train_data_logic_SS_l_1520 = raster_data.train_data_logic_SS_l_1520 ;
    train_data_logic_CS_l  = raster_data.train_data_logic_CS_l ;
    train_data_logic_CS_l_1  = raster_data.train_data_logic_CS_l_1 ;
    train_data_logic_CS_l_0  = raster_data.train_data_logic_CS_l_0 ;
    train_data_logic_lick_onset_l     = raster_data.train_data_logic_lick_onset_l ;
    train_data_logic_lick_vmax_l   = raster_data.train_data_logic_lick_vmax_l ;
    train_data_logic_lick_dmax_l    = raster_data.train_data_logic_lick_dmax_l ;
    train_data_logic_lick_vmin_l   = raster_data.train_data_logic_lick_vmin_l ;
    train_data_logic_lick_offset_l     = raster_data.train_data_logic_lick_offset_l ;
end


% Str bout licks
train_data_logic_SS_str_bout = raster_data.train_data_logic_SS_str_bout ;
train_data_logic_CS_str_bout  = raster_data.train_data_logic_CS_str_bout ;
train_data_logic_lick_onset_str_bout     = raster_data.train_data_logic_lick_onset_str_bout ;
train_data_logic_lick_vmax_str_bout   = raster_data.train_data_logic_lick_vmax_str_bout ;
train_data_logic_lick_dmax_str_bout    = raster_data.train_data_logic_lick_dmax_str_bout ;
train_data_logic_lick_vmin_str_bout   = raster_data.train_data_logic_lick_vmin_str_bout ;
train_data_logic_lick_offset_str_bout     = raster_data.train_data_logic_lick_offset_str_bout ;
% % Left bout
% if (session_type  == 2)
%     train_data_logic_SS_str_bout_l= raster_data.train_data_logic_SS_str_bout_l ;
%     train_data_logic_CS_str_bout_l  = raster_data.train_data_logic_CS_str_bout_l ;
%     train_data_logic_lick_onset_str_bout_l     = raster_data.train_data_logic_lick_onset_str_bout_l ;
% end
% % Right bout
% train_data_logic_SS_str_bout_r = raster_data.train_data_logic_SS_str_bout_r ;
% train_data_logic_CS_str_bout_r  = raster_data.train_data_logic_CS_str_bout_r ;
% train_data_logic_lick_onset_str_bout_r     = raster_data.train_data_logic_lick_onset_str_bout_r ;

% End bout licks
train_data_logic_SS_end_bout = raster_data.train_data_logic_SS_end_bout ;
train_data_logic_CS_end_bout  = raster_data.train_data_logic_CS_end_bout ;
train_data_logic_lick_onset_end_bout     = raster_data.train_data_logic_lick_onset_end_bout ;
train_data_logic_lick_vmax_end_bout   = raster_data.train_data_logic_lick_vmax_end_bout ;
train_data_logic_lick_dmax_end_bout    = raster_data.train_data_logic_lick_dmax_end_bout ;
train_data_logic_lick_vmin_end_bout   = raster_data.train_data_logic_lick_vmin_end_bout ;
train_data_logic_lick_offset_end_bout     = raster_data.train_data_logic_lick_offset_end_bout ;
% % Left bout
% if (session_type  == 2)
%     train_data_logic_SS_end_bout_l= raster_data.train_data_logic_SS_end_bout_l ;
%     train_data_logic_CS_end_bout_l  = raster_data.train_data_logic_CS_end_bout_l ;
%     train_data_logic_lick_onset_end_bout_l     = raster_data.train_data_logic_lick_onset_end_bout_l ;
% end
% % Right bout
% train_data_logic_SS_end_bout_r = raster_data.train_data_logic_SS_end_bout_r ;
% train_data_logic_CS_end_bout_r  = raster_data.train_data_logic_CS_end_bout_r ;
% train_data_logic_lick_onset_end_bout_r     = raster_data.train_data_logic_lick_onset_end_bout_r ;

% Saccades
% train_data_logic_primSac_onset   = raster_data.train_data_logic_primSac_onset ;
% train_data_logic_corrSac_onset   = raster_data.train_data_logic_corrSac_onset ;
train_data_logic_primSac_onset_grooming   = raster_data.train_data_logic_primSac_onset_grooming;
train_data_logic_corrSac_onset_grooming   = raster_data.train_data_logic_corrSac_onset_grooming ;
train_data_logic_primSac_onset_r   = raster_data.train_data_logic_primSac_onset_r;
train_data_logic_corrSac_onset_r   = raster_data.train_data_logic_corrSac_onset_r ;
if (session_type  == 2)
    train_data_logic_primSac_onset_l   = raster_data.train_data_logic_primSac_onset_l;
    train_data_logic_corrSac_onset_l   = raster_data.train_data_logic_corrSac_onset_l ;
end
train_data_logic_primSac_onset_str_bout   = raster_data.train_data_logic_primSac_onset_str_bout;
train_data_logic_corrSac_onset_str_bout   = raster_data.train_data_logic_corrSac_onset_str_bout ;
train_data_logic_primSac_onset_end_bout   = raster_data.train_data_logic_primSac_onset_end_bout;
train_data_logic_corrSac_onset_end_bout   = raster_data.train_data_logic_corrSac_onset_end_bout ;

% SS baseline based on bout [-1 0.5]s
SS_mean = nanmean(train_data_logic_SS_str_bout) * 100;
SS_baseline_mean = nanmean(SS_mean(1:50));

% % CS Probab
if contains(xlabel_text_CS_probab_, '+200')
    range_inds_probability = 31:50;
elseif contains(xlabel_text_CS_probab_, '-200')
    range_inds_probability = 11:30;
elseif contains(xlabel_text_CS_probab_, '100')
    range_inds_probability = 21:40;
elseif contains(xlabel_text_CS_probab_, '50')
    range_inds_probability = 26:35;
end

prob_grooming = nansum( nansum(train_data_logic_CS_grooming(:,range_inds_probability) ,2) > 0 ) / size(train_data_logic_CS_grooming, 1);
prob_r = nansum( nansum(train_data_logic_CS_r(:,range_inds_probability) ,2) > 0 ) / size(train_data_logic_CS_r, 1);
prob_all = nansum( nansum(train_data_logic_CS_all(:,range_inds_probability) ,2) > 0 ) / size(train_data_logic_CS_all, 1);
prob_mean  = nanmean(train_data_logic_CS_all);

% prob_r_1 = nansum( nansum(train_data_logic_CS_r_1(:,range_inds_probability) ,2) > 0 ) / size(train_data_logic_CS_r_1, 1);
% prob_r_0 = nansum( nansum(train_data_logic_CS_r_0(:,range_inds_probability) ,2) > 0 ) / size(train_data_logic_CS_r_0, 1);
% prob_r_1 = 0;
% prob_r_0 = 0;
% prob_amplitude_r = [prob_r_1 prob_r_0];

if (session_type  == 2)
    prob_l = nansum( nansum(train_data_logic_CS_l(:,range_inds_probability) ,2) > 0 ) / size(train_data_logic_CS_l, 1);
    %     prob_l_1 = nansum( nansum(train_data_logic_CS_l_1(:,range_inds_probability) ,2) > 0 ) / size(train_data_logic_CS_l_1, 1);
    %     prob_l_0 = nansum( nansum(train_data_logic_CS_l_0(:,range_inds_probability) ,2) > 0 ) / size(train_data_logic_CS_l_0, 1);
    %     prob_l_1 = 0;
    %     prob_l_0 = 0;
    prob_amplitude = [prob_l prob_grooming prob_r prob_all prob_mean];
    %     prob_amplitude_l = [prob_l_1 prob_l_0];
else
    prob_amplitude = [prob_grooming prob_r prob_all prob_mean];
    
end


% % plot xlim and ylim
range_SS_Firing = [0 300];

Line_Color = lines(7);
color_SS = Line_Color(1,:);
color_CS = Line_Color(7,:);
color_lick_onset = Line_Color(3,:);
color_lick_vmax = Line_Color(4,:);
color_lick_dmax = Line_Color(5,:);
color_lick_vmin = Line_Color(2,:);
color_lick_offset = Line_Color(6,:);

fig_handle_ = figure(fig_num_);
fig_handle_.WindowState = 'maximized';
clf(fig_handle_)

%% Grooming licks
train_data_logic_SS_corr_ = train_data_logic_SS_grooming_corr;
train_data_logic_CS_corr_ = train_data_logic_CS_grooming_corr;
train_data_logic_lick_onset_corr_ = train_data_logic_lick_onset_grooming_corr ;
train_data_logic_lick_vmax_corr_ = train_data_logic_lick_vmax_grooming_corr ;
train_data_logic_lick_dmax_corr_ = train_data_logic_lick_dmax_grooming_corr ;
train_data_logic_lick_vmin_corr_ = train_data_logic_lick_vmin_grooming_corr ;
train_data_logic_lick_offset_corr_ = train_data_logic_lick_offset_grooming_corr ;
train_data_logic_SS_ = train_data_logic_SS_grooming;
% train_data_logic_SS_510 = train_data_logic_SS_grooming_510;
% train_data_logic_SS_10 = train_data_logic_SS_grooming_10;
train_data_logic_CS_ = train_data_logic_CS_grooming;
train_data_logic_lick_onset_    = train_data_logic_lick_onset_grooming;
train_data_logic_lick_vmax_   = train_data_logic_lick_vmax_grooming;
train_data_logic_lick_dmax_ = train_data_logic_lick_dmax_grooming;
train_data_logic_lick_vmin_  = train_data_logic_lick_vmin_grooming;
train_data_logic_lick_offset_ = train_data_logic_lick_offset_grooming;

% SS
subplot(4,9,14)
yyaxis left;
hold on
[x_axis_SS_, y_axis_SS_] = ESN_raster_plot_axes(train_data_logic_SS_, inds_span, 0.5);
plot(x_axis_SS_(:), y_axis_SS_(:), 'LineWidth', 1.5, 'Color', color_SS)
if contains(xlabel_text_raster_, 'onset')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_onset_, inds_span, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_onset)
elseif contains(xlabel_text_raster_, 'vmax')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_vmax_, inds_span, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_vmax)
elseif contains(xlabel_text_raster_, 'dmax')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_dmax_, inds_span, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_dmax)
elseif contains(xlabel_text_raster_, 'vmin')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_vmin_, inds_span, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_vmin)
elseif contains(xlabel_text_raster_, 'offset')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_offset_, inds_span, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_offset)
end
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic_CS_,1)+3)])
ylabel('Licks')
title('SS | G licks')
set(gca, 'YColor', color_lick_onset)
yyaxis right;
yline(SS_baseline_mean,'m', 'LineWidth', 2);
firing_SS_ = nanmean(train_data_logic_SS_) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_,1));
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);

% subplot(4,9,23)
% firing_SS_510 =  nanmean(train_data_logic_SS_510) * 100;
% firing_SS_10 =  nanmean(train_data_logic_SS_10) * 100;
% hold on
% plot(inds_span, ESN_smooth(firing_SS_510,5), '--k', 'LineWidth', 2)
% plot(inds_span, ESN_smooth(firing_SS_10,5), 'k', 'LineWidth', 2)
% ylabel('SS Firing (spk/s)')
% ylim(range_SS_Firing)
% xlim([min(inds_span)-1 max(inds_span)+1])
% xlabel(xlabel_text_raster_);
% title('SS | G licks | -- < -')
% xline(0,'LineWidth', 2, 'Color', color_lick_onset);

% CS
subplot(4,9,17)
hold on
yyaxis left;
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_CS_, inds_span, 3);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 3, 'Color', color_CS)
if contains(xlabel_text_raster_, 'onset')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_onset_, inds_span, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_onset)
elseif contains(xlabel_text_raster_, 'vmax')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_vmax_, inds_span, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_vmax)
elseif contains(xlabel_text_raster_, 'dmax')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_dmax_, inds_span, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_dmax)
elseif contains(xlabel_text_raster_, 'vmin')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_vmin_, inds_span, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_vmin)
elseif contains(xlabel_text_raster_, 'offset')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_offset_, inds_span, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_offset)
end
ylabel('Licks')
ylim([(1-3) (size(train_data_logic_CS_grooming,1)+3)])
set(gca, 'YColor', color_lick_onset)
yyaxis right;
firing_CS_= nanmean(train_data_logic_CS_) * 100;
firing_CS_sd = nanstd(train_data_logic_CS_) * 100;
firing_CS_se = firing_CS_sd/sqrt(size(train_data_logic_CS_,1));
plot(inds_span, ESN_smooth(firing_CS_,10), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_CS_ + firing_CS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_CS_ - firing_CS_se,5), '-k', 'LineWidth', 1);
ylabel('CS Firing (spk/s)')
xlabel(xlabel_text_raster_);
title('CS | G licks')
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([0 3])
xlabel(xlabel_text_raster_);
set(gca, 'YColor', 'k')

% auto corr
if contains(xlabel_text_raster_, 'onset')
    lick_auto_corr = nansum(train_data_logic_lick_onset_corr_)/nansum(nansum(train_data_logic_lick_onset_corr_));
elseif contains(xlabel_text_raster_, 'vmax')
    lick_auto_corr = nansum(train_data_logic_lick_vmax_corr_)/nansum(nansum(train_data_logic_lick_vmax_corr_));
elseif contains(xlabel_text_raster_, 'dmax')
    lick_auto_corr = nansum(train_data_logic_lick_dmax_corr_)/nansum(nansum(train_data_logic_lick_dmax_corr_));
elseif contains(xlabel_text_raster_, 'dvmin')
    lick_auto_corr = nansum(train_data_logic_lick_vmin_corr_)/nansum(nansum(train_data_logic_lick_vmin_corr_));
elseif contains(xlabel_text_raster_, 'offset')
    lick_auto_corr = nansum(train_data_logic_lick_offset_corr_)/nansum(nansum(train_data_logic_lick_offset_corr_));
end
lick_auto_corr(100) = 0;
if (length(lick_auto_corr) < length(inds_span_corr))
    lick_auto_corr = zeros(1,length(inds_span_corr));
end

% cross corr
lick_cross_corr_SS = nansum(train_data_logic_SS_corr_)/nansum(nansum(train_data_logic_SS_corr_));
lick_cross_corr_CS = nansum(train_data_logic_CS_corr_)/nansum(nansum(train_data_logic_CS_corr_));

subplot(4,9,29)
% plot(inds_span_corr, ESN_smooth(lick_auto_corr,5), 'k','LineWidth', 2)
area(inds_span_corr,ESN_smooth(lick_auto_corr,5), 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 1.5)
title('Prob lick | G licks')
xlabel(xlabel_text_raster_bout_)
xlim([-1 1])
ylabel('Probability')
if (isnan(max(lick_auto_corr)) || max(lick_auto_corr) == 0)
    ylim([0 0.1])
else
    ylim([0 (max(lick_auto_corr) + max(lick_auto_corr)/10)])
end
set(gca, 'YColor', 'k')
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'vmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmax);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'vmin')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmin);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end

subplot(4,9,32)
% plot(inds_span_corr, ESN_smooth(lick_cross_corr_SS,5), 'k','LineWidth', 2)
area(inds_span_corr,ESN_smooth(lick_cross_corr_SS,5), 'FaceColor', 'k', 'LineWidth', 1.5)
title('Prob SS | G licks')
xlabel(xlabel_text_raster_bout_)
xlim([-1 1])
ylabel('Probability')
if (isnan(max(lick_cross_corr_SS)))
    ylim([0 0.1])
else
    ylim([0 (max(lick_cross_corr_SS) + max(lick_cross_corr_SS)/10)])
end
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'vmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmax);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'vmin')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmin);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end

subplot(4,9,35)
% plot(inds_span_corr, ESN_smooth(lick_cross_corr_CS,5), 'k','LineWidth', 2)
area(inds_span_corr,ESN_smooth(lick_cross_corr_CS,10), 'FaceColor', 'k', 'LineWidth', 1)
title('Prob CS | G licks')
xlabel(xlabel_text_raster_bout_)
xlim([-1 1])
ylabel('Probability')
if (isnan(max(lick_cross_corr_CS)))
    ylim([0 0.1])
else
    ylim([0 (max(lick_cross_corr_CS) + max(lick_cross_corr_CS)/10)])
end
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'vmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmax);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'vmin')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmin);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
%% Grooming Sacs
train_data_logic_primSac_onset_ = train_data_logic_primSac_onset_grooming;
train_data_logic_corrSac_onset_ = train_data_logic_corrSac_onset_grooming;

train_data_logic_Sac_onset_ = train_data_logic_primSac_onset_ + train_data_logic_corrSac_onset_;
mean_sacc_ = nanmean(train_data_logic_Sac_onset_);


subplot(4,9,11)
hold on
% plot(inds_span, ESN_smooth(mean_sacc_,5), 'LineWidth', 2, 'Color', 'k')
area(inds_span,ESN_smooth(mean_sacc_,5), 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 1.5)
ylabel('Prob sacc onset')
xlabel(xlabel_text_raster_);
set(gca, 'YColor', 'k')
title('Prob sacc | G licks')
xlim([min(inds_span)-1 max(inds_span)+1])
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'vmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmax);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'vmin')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmin);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
ylim([0 0.05])

%% R licks
train_data_logic_SS_corr_ = train_data_logic_SS_r_corr;
train_data_logic_CS_corr_ = train_data_logic_CS_r_corr;
train_data_logic_lick_onset_corr_ = train_data_logic_lick_onset_r_corr ;
train_data_logic_lick_vmax_corr_ = train_data_logic_lick_vmax_r_corr ;
train_data_logic_lick_dmax_corr_ = train_data_logic_lick_dmax_r_corr ;
train_data_logic_lick_vmin_corr_ = train_data_logic_lick_vmin_r_corr ;
train_data_logic_lick_offset_corr_ = train_data_logic_lick_offset_r_corr ;
train_data_logic_SS_ = train_data_logic_SS_r;
% train_data_logic_SS_15 = train_data_logic_SS_r_15;
% train_data_logic_SS_1520 = train_data_logic_SS_r_1520;
train_data_logic_CS_ = train_data_logic_CS_r;
train_data_logic_lick_onset_    = train_data_logic_lick_onset_r;
train_data_logic_lick_vmax_   = train_data_logic_lick_vmax_r;
train_data_logic_lick_dmax_ = train_data_logic_lick_dmax_r;
train_data_logic_lick_vmin_  = train_data_logic_lick_vmin_r;
train_data_logic_lick_offset_ = train_data_logic_lick_offset_r;

% SS
subplot(4,9,15)
yyaxis left
hold on
[x_axis_SS_, y_axis_SS_] = ESN_raster_plot_axes(train_data_logic_SS_, inds_span, 0.5);
plot(x_axis_SS_(:), y_axis_SS_(:), 'LineWidth', 1.5, 'Color', color_SS)
if contains(xlabel_text_raster_, 'onset')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_onset_, inds_span, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_onset)
elseif contains(xlabel_text_raster_, 'vmax')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_vmax_, inds_span, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_vmax)
elseif contains(xlabel_text_raster_, 'dmax')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_dmax_, inds_span, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_dmax)
elseif contains(xlabel_text_raster_, 'vmin')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_vmin_, inds_span, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_vmin)
elseif contains(xlabel_text_raster_, 'offset')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_offset_, inds_span, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_offset)
end
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic_CS_,1)+3)])
ylabel('Licks')
title('SS | R licks')
set(gca, 'YColor', color_lick_onset)
yyaxis right;
yline(SS_baseline_mean,'m', 'LineWidth', 2);
firing_SS_ = nanmean(train_data_logic_SS_) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_,1));
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
set(gca, 'YColor', 'k')

% subplot(4,9,24)
% firing_SS_15 =  nanmean(train_data_logic_SS_15) * 100;
% firing_SS_1520 =  nanmean(train_data_logic_SS_1520) * 100;
% hold on
% plot(inds_span, ESN_smooth(firing_SS_15,5),  '--k',  'LineWidth', 2 )
% plot(inds_span, ESN_smooth(firing_SS_1520,5), 'k', 'LineWidth', 2)
% ylabel('SS Firing (spk/s)')
% ylim(range_SS_Firing)
% xlim([min(inds_span)-1 max(inds_span)+1])
% xlabel(xlabel_text_raster_);
% title('SS | R licks | -- < -')
% xline(0,'LineWidth', 2, 'Color', color_lick_onset);

% CS
subplot(4,9,18)
hold on
yyaxis left;
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_CS_, inds_span, 3);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 3, 'Color', color_CS)
if contains(xlabel_text_raster_, 'onset')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_onset_, inds_span, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_onset)
elseif contains(xlabel_text_raster_, 'vmax')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_vmax_, inds_span, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_vmax)
elseif contains(xlabel_text_raster_, 'dmax')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_dmax_, inds_span, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_dmax)
elseif contains(xlabel_text_raster_, 'vmin')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_vmin_, inds_span, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_vmin)
elseif contains(xlabel_text_raster_, 'offset')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_offset_, inds_span, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_offset)
end
ylabel('Licks')
ylim([(1-3) (size(train_data_logic_CS_r,1)+3)])
title('CS | R licks')
set(gca, 'YColor', color_lick_onset)
yyaxis right;
firing_CS_ = nanmean(train_data_logic_CS_) * 100;
firing_CS_sd = nanstd(train_data_logic_CS_) * 100;
firing_CS_se = firing_CS_sd/sqrt(size(train_data_logic_CS_,1));
plot(inds_span, ESN_smooth(firing_CS_,10), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_CS_ + firing_CS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_CS_ - firing_CS_se,5), '-k', 'LineWidth', 1);
ylabel('CS Firing (spk/s)')
xlabel(xlabel_text_raster_);
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([0 3])
xlabel(xlabel_text_raster_);
set(gca, 'YColor', 'k')

% auto corr
if contains(xlabel_text_raster_, 'onset')
    lick_auto_corr = nansum(train_data_logic_lick_onset_corr_)/nansum(nansum(train_data_logic_lick_onset_corr_));
elseif contains(xlabel_text_raster_, 'vmax')
    lick_auto_corr = nansum(train_data_logic_lick_vmax_corr_)/nansum(nansum(train_data_logic_lick_vmax_corr_));
elseif contains(xlabel_text_raster_, 'dmax')
    lick_auto_corr = nansum(train_data_logic_lick_dmax_corr_)/nansum(nansum(train_data_logic_lick_dmax_corr_));
elseif contains(xlabel_text_raster_, 'vmin')
    lick_auto_corr = nansum(train_data_logic_lick_vmin_corr_)/nansum(nansum(train_data_logic_lick_vmin_corr_));
elseif contains(xlabel_text_raster_, 'offset')
    lick_auto_corr = nansum(train_data_logic_lick_offset_corr_)/nansum(nansum(train_data_logic_lick_offset_corr_));
end
lick_auto_corr(100) = 0;
if (length(lick_auto_corr) < length(inds_span_corr))
    lick_auto_corr = zeros(1,length(inds_span_corr));
end

% cross corr
lick_cross_corr_SS = nansum(train_data_logic_SS_corr_)/nansum(nansum(train_data_logic_SS_corr_));
lick_cross_corr_CS = nansum(train_data_logic_CS_corr_)/nansum(nansum(train_data_logic_CS_corr_));

subplot(4,9,30)
% plot(inds_span_corr, ESN_smooth(lick_auto_corr,5), 'k','LineWidth', 1)
area(inds_span_corr,ESN_smooth(lick_auto_corr,5), 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 1.5)
title('Prob lick | R licks')
xlabel(xlabel_text_raster_bout_);
xlim([-1 1])
ylabel('Probability')
if (isnan(max(lick_auto_corr)) || max(lick_auto_corr) == 0)
    ylim([0 0.1])
else
    ylim([0 (max(lick_auto_corr) + max(lick_auto_corr)/10)])
end
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'vmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmax);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'vmin')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmin);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end

subplot(4,9,33)
% plot(inds_span_corr, ESN_smooth(lick_cross_corr_SS,5), 'k','LineWidth', 1)
area(inds_span_corr,ESN_smooth(lick_cross_corr_SS,5), 'FaceColor', 'k', 'LineWidth', 1.5)
title('Prob SS | R licks')
xlabel(xlabel_text_raster_bout_);
xlim([-1 1])
ylabel('Probability')
if (isnan(max(lick_cross_corr_SS)))
    ylim([0 0.1])
else
    ylim([0 (max(lick_cross_corr_SS) + max(lick_cross_corr_SS)/10)])
end

if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'vmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmax);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'vmin')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmin);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end

subplot(4,9,36)
% plot(inds_span_corr, ESN_smooth(lick_cross_corr_CS,5), 'k','LineWidth', 1)
area(inds_span_corr,ESN_smooth(lick_cross_corr_CS,10), 'FaceColor', 'k', 'LineWidth', 1)
title('Prob CS | R licks')
xlabel(xlabel_text_raster_bout_);
xlim([-1 1])
ylabel('Probability')
if (isnan(max(lick_cross_corr_CS)))
    ylim([0 0.1])
else
    ylim([0 (max(lick_cross_corr_CS) + max(lick_cross_corr_CS)/10)])
end
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'vmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmax);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'vmin')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmin);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
%% R Sacs
train_data_logic_primSac_onset_ = train_data_logic_primSac_onset_r;
train_data_logic_corrSac_onset_ = train_data_logic_corrSac_onset_r;

train_data_logic_Sac_onset_ = train_data_logic_primSac_onset_ + train_data_logic_corrSac_onset_;

mean_sacc_ = nanmean(train_data_logic_Sac_onset_);
subplot(4,9,12)
hold on
% plot(inds_span, ESN_smooth(mean_sacc_, 5), 'LineWidth', 2, 'Color', 'k')
area(inds_span,ESN_smooth(mean_sacc_,5), 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 1.5)
ylabel('Prob sacc onset')
xlabel(xlabel_text_raster_);
set(gca, 'YColor', 'k')
title('Prob sacc | R licks')
xlim([min(inds_span)-1 max(inds_span)+1])
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'vmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmax);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'vmin')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmin);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
ylim([0 0.05])

%% L licks
if (session_type  == 2)
    train_data_logic_SS_corr_ = train_data_logic_SS_l_corr;
    train_data_logic_CS_corr_ = train_data_logic_CS_l_corr;
    train_data_logic_lick_onset_corr_ = train_data_logic_lick_onset_l_corr ;
    train_data_logic_lick_vmax_corr_ = train_data_logic_lick_vmax_l_corr ;
    train_data_logic_lick_dmax_corr_ = train_data_logic_lick_dmax_l_corr ;
    train_data_logic_lick_vmin_corr_ = train_data_logic_lick_vmin_l_corr ;
    train_data_logic_lick_offset_corr_ = train_data_logic_lick_offset_l_corr ;
    train_data_logic_SS_ = train_data_logic_SS_l;
    %     train_data_logic_SS_15 = train_data_logic_SS_l_15;
    %     train_data_logic_SS_1520 = train_data_logic_SS_l_1520;
    train_data_logic_CS_ = train_data_logic_CS_l;
    train_data_logic_lick_onset_    = train_data_logic_lick_onset_l;
    train_data_logic_lick_vmax_   = train_data_logic_lick_vmax_l;
    train_data_logic_lick_dmax_ = train_data_logic_lick_dmax_l;
    train_data_logic_lick_vmin_  = train_data_logic_lick_vmin_l;
    train_data_logic_lick_offset_ = train_data_logic_lick_offset_l;
    
    % SS
    subplot(4,9,13)
    yyaxis left;
    hold on
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_SS_, inds_span, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_SS)
    if contains(xlabel_text_raster_, 'onset')
        [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_onset_, inds_span, 0.5);
        plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_onset)
    elseif contains(xlabel_text_raster_, 'vmax')
        [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_vmax_, inds_span, 0.5);
        plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_vmax)
    elseif contains(xlabel_text_raster_, 'dmax')
        [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_dmax_, inds_span, 0.5);
        plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_dmax)
    elseif contains(xlabel_text_raster_, 'vmin')
        [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_vmin_, inds_span, 0.5);
        plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_vmin)
    elseif contains(xlabel_text_raster_, 'offset')
        [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_offset_, inds_span, 0.5);
        plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_offset)
    end
    xlim([min(inds_span)-1 max(inds_span)+1])
    ylim([(1-3) (size(train_data_logic_CS_,1)+3)])
    ylabel('Licks')
    title('SS | L licks')
    set(gca, 'YColor',color_lick_onset)
    yyaxis right;
    yline(SS_baseline_mean,'m', 'LineWidth', 2);
    firing_SS_ = nanmean(train_data_logic_SS_) * 100;
    firing_SS_sd = nanstd(train_data_logic_SS_) * 100;
    firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_,1));
    plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
    plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
    plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
    ylabel('SS Firing (spk/s)')
    ylim(range_SS_Firing)
    xlim([min(inds_span)-1 max(inds_span)+1])
    xlabel(xlabel_text_raster_);
    set(gca, 'YColor', 'k')
    
    %     subplot(4,9,22)
    %     firing_SS_15 =  nanmean(train_data_logic_SS_15) * 100;
    %     firing_SS_1520 =  nanmean(train_data_logic_SS_1520) * 100;
    %     %     firing_SS_20 =  nanmean(train_data_logic_SS_20) * 100;
    %     hold on
    %     plot(inds_span, ESN_smooth(firing_SS_15,5), '--k', 'LineWidth', 2 )
    %     plot(inds_span, ESN_smooth(firing_SS_1520,5), 'k', 'LineWidth', 2)
    %     ylabel('SS Firing (spk/s)')
    %     ylim(range_SS_Firing)
    %     xlim([min(inds_span)-1 max(inds_span)+1])
    %     xlabel(xlabel_text_raster_);
    %     title('SS | L licks | -- < -')
    %     xline(0,'LineWidth', 2, 'Color', color_lick_onset);
    %
    
    % CS
    subplot(4,9,16)
    hold on
    yyaxis left;
    [x_axis_SS_, y_axis_SS_] = ESN_raster_plot_axes(train_data_logic_CS_, inds_span, 3);
    plot(x_axis_SS_(:), y_axis_SS_(:), 'LineWidth', 3, 'Color', color_CS)
    if contains(xlabel_text_raster_, 'onset')
        [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_onset_, inds_span, 0.5);
        plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_onset)
    elseif contains(xlabel_text_raster_, 'vmax')
        [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_vmax_, inds_span, 0.5);
        plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_vmax)
    elseif contains(xlabel_text_raster_, 'dmax')
        [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_dmax_, inds_span, 0.5);
        plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_dmax)
    elseif contains(xlabel_text_raster_, 'vmin')
        [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_vmin_, inds_span, 0.5);
        plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_vmin)
    elseif contains(xlabel_text_raster_, 'offset')
        [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_offset_, inds_span, 0.5);
        plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_offset)
    end
    ylabel('Licks')
    ylim([(1-3) (size(train_data_logic_CS_,1)+3)])
    title('CS | L licks')
    set(gca, 'YColor', color_lick_onset)
    yyaxis right;
    firing_CS_ = nanmean(train_data_logic_CS_) * 100;
    firing_CS_sd = nanstd(train_data_logic_CS_) * 100;
    firing_CS_se = firing_CS_sd/sqrt(size(train_data_logic_CS_,1));
    plot(inds_span, ESN_smooth(firing_CS_,10), '-k', 'LineWidth', 2);
    % plot(inds_span, ESN_smooth(firing_CS_ + firing_CS_se,5), '-k', 'LineWidth', 1);
    % plot(inds_span, ESN_smooth(firing_CS_ - firing_CS_se,5), '-k', 'LineWidth', 1);
    ylabel('CS Firing (spk/s)')
    xlabel(xlabel_text_raster_);
    xlim([min(inds_span)-1 max(inds_span)+1])
    ylim([0 3])
    xlabel(xlabel_text_raster_);
    set(gca, 'YColor', 'k')
    
    % auto corr
    if contains(xlabel_text_raster_, 'onset')
        lick_auto_corr = nansum(train_data_logic_lick_onset_corr_)/nansum(nansum(train_data_logic_lick_onset_corr_));
    elseif contains(xlabel_text_raster_, 'vmax')
        lick_auto_corr = nansum(train_data_logic_lick_vmax_corr_)/nansum(nansum(train_data_logic_lick_vmax_corr_));
    elseif contains(xlabel_text_raster_, 'dmax')
        lick_auto_corr = nansum(train_data_logic_lick_dmax_corr_)/nansum(nansum(train_data_logic_lick_dmax_corr_));
    elseif contains(xlabel_text_raster_, 'vmin')
        lick_auto_corr = nansum(train_data_logic_lick_vmin_corr_)/nansum(nansum(train_data_logic_lick_vmin_corr_));
    elseif contains(xlabel_text_raster_, 'offset')
        lick_auto_corr = nansum(train_data_logic_lick_offset_corr_)/nansum(nansum(train_data_logic_lick_offset_corr_));
    end
    lick_auto_corr(100) = 0;
    if (length(lick_auto_corr) < length(inds_span_corr))
        lick_auto_corr = zeros(1,length(inds_span_corr));
    end
    
    
    % cross corr
    lick_cross_corr_SS = nansum(train_data_logic_SS_corr_)/nansum(nansum(train_data_logic_SS_corr_));
    lick_cross_corr_CS = nansum(train_data_logic_CS_corr_)/nansum(nansum(train_data_logic_CS_corr_));
    
    subplot(4,9,28)
    % plot(inds_span_corr, ESN_smooth(lick_auto_corr,5), 'k','LineWidth', 1)
    area(inds_span_corr,ESN_smooth(lick_auto_corr,5), 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 1.5)
    title('Prob Lick | L licks')
    xlabel(xlabel_text_raster_bout_);
    xlim([-1 1])
    ylabel('Probability')
    if (isnan(max(lick_auto_corr)) || max(lick_auto_corr) == 0)
        ylim([0 0.1])
    else
        ylim([0 (max(lick_auto_corr) + max(lick_auto_corr)/10)])
    end
    
    if contains(xlabel_text_raster_, 'onset')
        xline(0,'LineWidth', 2, 'Color', color_lick_onset);
    elseif contains(xlabel_text_raster_, 'vmax')
        xline(0,'LineWidth', 2, 'Color', color_lick_vmax);
    elseif contains(xlabel_text_raster_, 'dmax')
        xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
    elseif contains(xlabel_text_raster_, 'vmin')
        xline(0,'LineWidth', 2, 'Color', color_lick_vmin);
    elseif contains(xlabel_text_raster_, 'offset')
        xline(0,'LineWidth', 2, 'Color', color_lick_offset);
    end
    
    subplot(4,9,31)
    % plot(inds_span_corr, ESN_smooth(lick_cross_corr_SS,5), 'k','LineWidth', 1)
    area(inds_span_corr,ESN_smooth(lick_cross_corr_SS,5), 'FaceColor', 'k', 'LineWidth', 1.5)
    title('Prob SS | L licks')
    xlabel(xlabel_text_raster_bout_);
    xlim([-1 1])
    ylabel('Probability')
    if (isnan(max(lick_cross_corr_SS)))
        ylim([0 0.1])
    else
        ylim([0 (max(lick_cross_corr_SS) + max(lick_cross_corr_SS)/10)])
    end
    
    if contains(xlabel_text_raster_, 'onset')
        xline(0,'LineWidth', 2, 'Color', color_lick_onset);
    elseif contains(xlabel_text_raster_, 'vmax')
        xline(0,'LineWidth', 2, 'Color', color_lick_vmax);
    elseif contains(xlabel_text_raster_, 'dmax')
        xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
    elseif contains(xlabel_text_raster_, 'vmin')
        xline(0,'LineWidth', 2, 'Color', color_lick_vmin);
    elseif contains(xlabel_text_raster_, 'offset')
        xline(0,'LineWidth', 2, 'Color', color_lick_offset);
    end
    
    subplot(4,9,34)
    % plot(inds_span_corr, ESN_smooth(lick_cross_corr_CS,5), 'k','LineWidth', 1)
    area(inds_span_corr,ESN_smooth(lick_cross_corr_CS,10), 'FaceColor', 'k', 'LineWidth', 1)
    title('Prob CS | L licks')
    xlabel(xlabel_text_raster_bout_);
    xlim([-1 1])
    ylabel('Probability')
    if (isnan(max(lick_cross_corr_CS)))
        ylim([0 0.1])
    else
        ylim([0 (max(lick_cross_corr_CS) + max(lick_cross_corr_CS)/10)])
    end
end
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'vmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmax);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'vmin')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmin);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
%% L Sacs
if (session_type  == 2)
    train_data_logic_primSac_onset_ = train_data_logic_primSac_onset_l;
    train_data_logic_corrSac_onset_ = train_data_logic_corrSac_onset_l;
    
    train_data_logic_Sac_onset_ = train_data_logic_primSac_onset_ + train_data_logic_corrSac_onset_;
    
    mean_sacc_ = nanmean(train_data_logic_Sac_onset_);
    subplot(4,9,10)
    hold on
    %     plot(inds_span, ESN_smooth(mean_sacc_,5), 'LineWidth', 2, 'Color', 'k')
    area(inds_span,ESN_smooth(mean_sacc_,5), 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 1.5)
    ylabel('Prob sacc onset')
    xlabel(xlabel_text_raster_);
    set(gca, 'YColor', 'k')
    title('Prob sacc | L licks')
    xlim([min(inds_span)-1 max(inds_span)+1])
    if contains(xlabel_text_raster_, 'onset')
        xline(0,'LineWidth', 2, 'Color', color_lick_onset);
    elseif contains(xlabel_text_raster_, 'vmax')
        xline(0,'LineWidth', 2, 'Color', color_lick_vmax);
    elseif contains(xlabel_text_raster_, 'dmax')
        xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
    elseif contains(xlabel_text_raster_, 'vmin')
        xline(0,'LineWidth', 2, 'Color', color_lick_vmin);
    elseif contains(xlabel_text_raster_, 'offset')
        xline(0,'LineWidth', 2, 'Color', color_lick_offset);
    end
    ylim([0 0.05])
    
end

%% Str bout licks
train_data_logic_SS_ = train_data_logic_SS_str_bout;
train_data_logic_CS_ = train_data_logic_CS_str_bout;
train_data_logic_lick_onset_    = train_data_logic_lick_onset_str_bout;
train_data_logic_lick_vmax_   = train_data_logic_lick_vmax_str_bout;
train_data_logic_lick_dmax_ = train_data_logic_lick_dmax_str_bout;
train_data_logic_lick_vmin_  = train_data_logic_lick_vmin_str_bout;
train_data_logic_lick_offset_ = train_data_logic_lick_offset_str_bout;

% SS
subplot(4,9,4)
hold on
yyaxis left;
[x_axis_SS_, y_axis_SS_] = ESN_raster_plot_axes(train_data_logic_SS_, inds_span_bout, 0.5);
plot(x_axis_SS_(:), y_axis_SS_(:), 'LineWidth', 0.5, 'Color', color_SS)
if contains(xlabel_text_raster_bout_, 'onset')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_onset_, inds_span_bout, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_onset)
elseif contains(xlabel_text_raster_bout_, 'vmax')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_vmax_, inds_span_bout, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_vmax)
elseif contains(xlabel_text_raster_bout_, 'dmax')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_dmax_, inds_span_bout, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_dmax)
elseif contains(xlabel_text_raster_bout_, 'vmin')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_vmin_, inds_span_bout, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_vmin)
elseif contains(xlabel_text_raster_bout_, 'offset')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_offset_, inds_span_bout, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_offset)
end
xlim([min(inds_span_bout)-0.01 max(inds_span_bout)])
ylim([(1-3) (size(train_data_logic_CS_,1)+3)])
ylabel('Licks')
title('SS | Start bout licks')
xlabel(xlabel_text_raster_bout_);
set(gca, 'YColor', color_lick_onset)
yyaxis right;
yline(SS_baseline_mean,'m', 'LineWidth', 2);
firing_SS_ = nanmean(train_data_logic_SS_) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_,1));
plot(inds_span_bout, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
plot(inds_span_bout, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
plot(inds_span_bout, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
xlim([min(inds_span_bout)-0.01 max(inds_span_bout)])
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')

% CS
subplot(4,9,7)
hold on
yyaxis left;
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_CS_, inds_span_bout, 1);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_CS)
if contains(xlabel_text_raster_bout_, 'onset')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_onset_, inds_span_bout, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_onset)
elseif contains(xlabel_text_raster_bout_, 'vmax')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_vmax_, inds_span_bout, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_vmax)
elseif contains(xlabel_text_raster_bout_, 'dmax')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_dmax_, inds_span_bout, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_dmax)
elseif contains(xlabel_text_raster_bout_, 'vmin')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_vmin_, inds_span_bout, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_vmin)
elseif contains(xlabel_text_raster_bout_, 'offset')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_offset_, inds_span_bout, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_offset)
end
xlim([min(inds_span_bout)-0.01 max(inds_span_bout)])
ylim([(1-3) (size(train_data_logic_CS_,1)+3)])
ylabel('Licks')
title(' CS | Str bout licks')
xlabel(xlabel_text_raster_bout_);
set(gca, 'YColor', color_lick_onset)
yyaxis right;
firing_CS_ = nanmean(train_data_logic_CS_) * 100;
firing_CS_sd = nanstd(train_data_logic_CS_) * 100;
firing_CS_se = firing_CS_sd/sqrt(size(train_data_logic_CS_,1));
plot(inds_span_bout, ESN_smooth(firing_CS_,50), '-k', 'LineWidth', 2);
% plot(inds_span_bout, ESN_smooth(firing_CS_ + firing_CS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span_bout, ESN_smooth(firing_CS_ - firing_CS_se,5), '-k', 'LineWidth', 1);
ylabel('CS Firing (spk/s)')
ylim([0 3])
xlim([min(inds_span_bout)-0.01 max(inds_span_bout)])
set(gca, 'YColor', 'k')
%% End bout licks
train_data_logic_SS_ = train_data_logic_SS_end_bout;
train_data_logic_CS_ = train_data_logic_CS_end_bout;
train_data_logic_lick_onset_    = train_data_logic_lick_onset_end_bout;
train_data_logic_lick_vmax_   = train_data_logic_lick_vmax_end_bout;
train_data_logic_lick_dmax_ = train_data_logic_lick_dmax_end_bout;
train_data_logic_lick_vmin_  = train_data_logic_lick_vmin_end_bout;
train_data_logic_lick_offset_ = train_data_logic_lick_offset_end_bout;

% SS
subplot(4,9,6)
yyaxis left;
hold on
[x_axis_SS_, y_axis_SS_] = ESN_raster_plot_axes(train_data_logic_SS_, inds_span_bout, 0.5);
plot(x_axis_SS_(:), y_axis_SS_(:), 'LineWidth', 0.5, 'Color', color_SS)
if contains(xlabel_text_raster_bout_, 'onset')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_onset_, inds_span_bout, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_onset)
elseif contains(xlabel_text_raster_bout_, 'vmax')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_vmax_, inds_span_bout, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_vmax)
elseif contains(xlabel_text_raster_bout_, 'dmax')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_dmax_, inds_span_bout, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_dmax)
elseif contains(xlabel_text_raster_bout_, 'vmin')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_vmin_, inds_span_bout, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_vmin)
elseif contains(xlabel_text_raster_bout_, 'offset')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_offset_, inds_span_bout, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_offset)
end
xlim([min(inds_span_bout)-0.01 max(inds_span_bout)])
ylim([(1-3) (size(train_data_logic_CS_,1)+3)])
ylabel('Licks')
title('SS | End bout licks')
set(gca, 'YColor', color_lick_onset)
yyaxis right;
yline(SS_baseline_mean,'m', 'LineWidth', 2);
firing_SS_ = nanmean(train_data_logic_SS_) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_,1));
plot(inds_span_bout, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
plot(inds_span_bout, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
plot(inds_span_bout, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span_bout)-0.01 max(inds_span_bout)])
xlabel(xlabel_text_raster_bout_);

% CS
subplot(4,9,9)
yyaxis left;
hold on
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_CS_, inds_span_bout, 1);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_CS)
if contains(xlabel_text_raster_bout_, 'onset')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_onset_, inds_span_bout, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_onset)
elseif contains(xlabel_text_raster_bout_, 'vmax')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_vmax_, inds_span_bout, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_vmax)
elseif contains(xlabel_text_raster_bout_, 'dmax')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_dmax_, inds_span_bout, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_dmax)
elseif contains(xlabel_text_raster_bout_, 'vmin')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_vmin_, inds_span_bout, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_vmin)
elseif contains(xlabel_text_raster_bout_, 'offset')
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_lick_offset_, inds_span_bout, 0.5);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 1.5, 'Color', color_lick_offset)
end
xlim([min(inds_span_bout)-0.01 max(inds_span_bout)])
ylim([(1-3) (size(train_data_logic_CS_,1)+3)])
ylabel('Licks')
title('CS | End bout licks')
xlabel(xlabel_text_raster_bout_);
set(gca, 'YColor', color_lick_onset)
yyaxis right;
firing_CS_ = nanmean(train_data_logic_CS_) * 100;
firing_CS_sd = nanstd(train_data_logic_CS_) * 100;
firing_CS_se = firing_CS_sd/sqrt(size(train_data_logic_CS_,1));
plot(inds_span_bout, ESN_smooth(firing_CS_,50), '-k', 'LineWidth', 2);
% plot(inds_span_bout, ESN_smooth(firing_CS_ + firing_CS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span_bout, ESN_smooth(firing_CS_ - firing_CS_se,5), '-k', 'LineWidth', 1);
ylabel('CS Firing (spk/s)')
xlim([min(inds_span_bout)-0.01 max(inds_span_bout)])
ylim([0 3])
set(gca, 'YColor', 'k')
%% In/Out bout comparison SS & CS
train_data_logic_SS_str_bout_L = train_data_logic_SS_str_bout(:,1:100);
mean_train_data_logic_SS_str_bout_L = nanmean(nanmean(train_data_logic_SS_str_bout_L));
train_data_logic_SS_str_bout_U = train_data_logic_SS_str_bout(:,101:200);
mean_train_data_logic_SS_str_bout_U = nanmean( nanmean(train_data_logic_SS_str_bout_U));

mean_train_data_logic_SS_str_bout_diff = mean_train_data_logic_SS_str_bout_U - mean_train_data_logic_SS_str_bout_L;

train_data_logic_CS_str_bout_L = train_data_logic_CS_str_bout(:,1:100);
mean_train_data_logic_CS_str_bout_L = nanmean(nanmean(train_data_logic_CS_str_bout_L));
train_data_logic_CS_str_bout_U = train_data_logic_CS_str_bout(:,101:200);
mean_train_data_logic_CS_str_bout_U = nanmean( nanmean(train_data_logic_CS_str_bout_U));

mean_train_data_logic_CS_str_bout_diff = mean_train_data_logic_CS_str_bout_U - mean_train_data_logic_CS_str_bout_L;

train_data_logic_SS_end_bout_L = train_data_logic_SS_end_bout(:,1:100);
mean_train_data_logic_SS_end_bout_L = nanmean(nanmean(train_data_logic_SS_end_bout_L));
train_data_logic_SS_end_bout_U = train_data_logic_SS_end_bout(:,101:200);
mean_train_data_logic_SS_end_bout_U = nanmean( nanmean(train_data_logic_SS_end_bout_U));

mean_train_data_logic_SS_end_bout_diff = mean_train_data_logic_SS_end_bout_U - mean_train_data_logic_SS_end_bout_L;

train_data_logic_CS_end_bout_L = train_data_logic_CS_end_bout(:,1:100);
mean_train_data_logic_CS_end_bout_L = nanmean(nanmean(train_data_logic_CS_end_bout_L));
train_data_logic_CS_end_bout_U = train_data_logic_CS_end_bout(:,101:200);
mean_train_data_logic_CS_end_bout_U = nanmean( nanmean(train_data_logic_CS_end_bout_U));

mean_train_data_logic_CS_end_bout_diff = mean_train_data_logic_CS_end_bout_U - mean_train_data_logic_CS_end_bout_L;

% mean_diff = [mean_train_data_logic_SS_str_bout_diff, mean_train_data_logic_SS_end_bout_diff, ...
%     mean_train_data_logic_CS_str_bout_diff , mean_train_data_logic_CS_end_bout_diff] * 100;

mean_diff_SS = [mean_train_data_logic_SS_str_bout_diff, mean_train_data_logic_SS_end_bout_diff] * 100;
mean_diff_CS = [mean_train_data_logic_CS_str_bout_diff, mean_train_data_logic_CS_end_bout_diff] * 100;
x_axis = categorical({'Str', 'End'});
x_axis = reordercats(x_axis, {'Str', 'End'});

subplot(4,9,5)
bar(x_axis, mean_diff_SS, 'k')
ylabel('Change in SS Firing (spk/s)')
ylim([-100 100])
title('SS | Bout difference')
set(gca, 'YColor', 'k')

subplot(4,9,8)
bar(x_axis, mean_diff_CS, 'k')
ylabel('Change in CS Firing (spk/s)')
ylim([-0.4 0.4])
title('CS | Bout difference')
set(gca, 'YColor', 'k')
%% In/Out bout comparison saccs
train_data_logic_Sac_onset_str_bout = train_data_logic_primSac_onset_str_bout + train_data_logic_corrSac_onset_str_bout;
train_data_logic_Sac_onset_end_bout = train_data_logic_primSac_onset_end_bout + train_data_logic_corrSac_onset_end_bout;

mean_sacc_str_bout = nanmean(train_data_logic_Sac_onset_str_bout);
mean_sacc_end_bout = nanmean(train_data_logic_Sac_onset_end_bout);

subplot(4,9,19)
% plot(inds_span_bout, ESN_smooth(mean_sacc_str_bout,30), 'LineWidth', 2, 'Color', 'k')
area(inds_span_bout,ESN_smooth(mean_sacc_str_bout,10), 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 1.5)
ylabel('Prob sacc onset')
xlabel(xlabel_text_raster_bout_);
set(gca, 'YColor', 'k')
title('Prob sacc | Str bout licks')
xlim([min(inds_span_bout)-0.01 max(inds_span_bout)])
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'vmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmax);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'vmin')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmin);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
ylim([0 0.05])


subplot(4,9,21)
hold on
% plot(inds_span_bout, ESN_smooth(mean_sacc_end_bout,30), 'LineWidth', 2, 'Color', 'k')
area(inds_span_bout,ESN_smooth(mean_sacc_end_bout,10), 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 1.5)
ylabel('Prob sacc onset')
xlabel(xlabel_text_raster_bout_);
set(gca, 'YColor', 'k')
title('Prob sacc | End bout licks')
xlim([min(inds_span_bout)-0.01 max(inds_span_bout)])
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'vmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmax);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'vmin')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmin);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
ylim([0 0.05])

train_data_logic_Sac_onset_str_bout_L = train_data_logic_Sac_onset_str_bout(:,1:100);
mean_train_data_logic_Sac_onset_str_bout_L = nanmean(nanmean(train_data_logic_Sac_onset_str_bout_L));
train_data_logic_Sac_onset_str_bout_U = train_data_logic_Sac_onset_str_bout(:,101:200);
mean_train_data_logic_Sac_onset_str_bout_U = nanmean( nanmean(train_data_logic_Sac_onset_str_bout_U));

train_data_logic_Sac_onset_end_bout_L = train_data_logic_Sac_onset_end_bout(:,1:100);
mean_train_data_logic_Sac_onset_end_bout_L = nanmean(nanmean(train_data_logic_Sac_onset_end_bout_L));
train_data_logic_Sac_onset_end_bout_U = train_data_logic_Sac_onset_end_bout(:,101:200);
mean_train_data_logic_Sac_onset_end_bout_U = nanmean( mean(train_data_logic_Sac_onset_end_bout_U));

mean_train_data_logic_Sac_onset_str_bout_diff = mean_train_data_logic_Sac_onset_str_bout_U - mean_train_data_logic_Sac_onset_str_bout_L;
mean_train_data_logic_Sac_onset_end_bout_diff = mean_train_data_logic_Sac_onset_end_bout_U - mean_train_data_logic_Sac_onset_end_bout_L;

mean_sac_diff = [mean_train_data_logic_Sac_onset_str_bout_diff mean_train_data_logic_Sac_onset_end_bout_diff];
x_axis = categorical({'Str', 'End'});
x_axis = reordercats(x_axis, {'Str', 'End'});


subplot(4,9,20)
bar(x_axis, mean_sac_diff, 'k')
ylabel('Change in prob sac onset')
ylim([-0.025 0.025])
title('Prob sacc | Bout diff')
% set(gca, 'YColor', color_SS)

%% d_tip
VID_d_tip_100_grooming = raster_data.VID_d_tip_100_grooming;
VID_d_tip_100_r = raster_data.VID_d_tip_100_r;

VID_d_tip_100_grooming_mean = nanmean(VID_d_tip_100_grooming);
VID_d_tip_100_grooming_sem = nanstd(VID_d_tip_100_grooming)/sqrt(size(VID_d_tip_100_grooming,1));
VID_d_tip_100_r_mean = nanmean(VID_d_tip_100_r);
VID_d_tip_100_r_sem = nanstd(VID_d_tip_100_r)/sqrt(size(VID_d_tip_100_r,1));

if (session_type  == 2)
    VID_d_tip_100_l = raster_data.VID_d_tip_100_l;
    VID_d_tip_100_l_mean = nanmean(VID_d_tip_100_l);
    VID_d_tip_100_l_sem = nanstd(VID_d_tip_100_l)/sqrt(size(VID_d_tip_100_l,1));
    
end

subplot(4,9,1)
hold on
plot(inds_span,VID_d_tip_100_grooming_mean,'g' , 'LineWidth', 2)
plot(inds_span,VID_d_tip_100_grooming_mean+VID_d_tip_100_grooming_sem, '--g' , 'LineWidth', 1)
plot(inds_span,VID_d_tip_100_grooming_mean-VID_d_tip_100_grooming_sem, '--g' , 'LineWidth', 1)

plot(inds_span,VID_d_tip_100_r_mean,'r' , 'LineWidth', 2)
plot(inds_span,VID_d_tip_100_r_mean+VID_d_tip_100_r_sem, '--r', 'LineWidth', 1)
plot(inds_span,VID_d_tip_100_r_mean-VID_d_tip_100_r_sem, '--r', 'LineWidth', 1)

if (session_type  == 2)
    plot(inds_span,VID_d_tip_100_l_mean ,'b' , 'LineWidth', 2)
    plot(inds_span,VID_d_tip_100_l_mean+VID_d_tip_100_l_sem, '--b' , 'LineWidth', 1)
    plot(inds_span,VID_d_tip_100_l_mean-VID_d_tip_100_l_sem, '--b','LineWidth', 1)
end
ylabel('mm')
ylim([0 20])
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'vmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmax);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'vmin')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmin);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
title('Lick displacement')
%% v_tip
VID_v_tip_100_grooming = raster_data.VID_v_tip_100_grooming;
VID_v_tip_100_r = raster_data.VID_v_tip_100_r;

VID_v_tip_100_grooming_mean = nanmean(VID_v_tip_100_grooming);
VID_v_tip_100_grooming_sem = nanstd(VID_v_tip_100_grooming)/sqrt(size(VID_v_tip_100_grooming,1));
VID_v_tip_100_r_mean = nanmean(VID_v_tip_100_r);
VID_v_tip_100_r_sem = nanstd(VID_v_tip_100_r)/sqrt(size(VID_v_tip_100_r,1));

if (session_type  == 2)
    VID_v_tip_100_l = raster_data.VID_v_tip_100_l;
    VID_v_tip_100_l_mean = nanmean(VID_v_tip_100_l);
    VID_v_tip_100_l_sem = nanstd(VID_v_tip_100_l)/sqrt(size(VID_v_tip_100_l,1));
end

subplot(4,9,2)
hold on
plot(inds_span,VID_v_tip_100_grooming_mean,'g' , 'LineWidth', 2)
plot(inds_span,VID_v_tip_100_grooming_mean+VID_v_tip_100_grooming_sem, '--g' , 'LineWidth', 1)
plot(inds_span,VID_v_tip_100_grooming_mean-VID_v_tip_100_grooming_sem, '--g' , 'LineWidth', 1)

plot(inds_span,VID_v_tip_100_r_mean,'r' , 'LineWidth', 2)
plot(inds_span,VID_v_tip_100_r_mean+VID_v_tip_100_r_sem, '--r', 'LineWidth', 1)
plot(inds_span,VID_v_tip_100_r_mean-VID_v_tip_100_r_sem, '--r', 'LineWidth', 1)

if (session_type  == 2)
    plot(inds_span,VID_v_tip_100_l_mean ,'b' , 'LineWidth', 2)
    plot(inds_span,VID_v_tip_100_l_mean+VID_v_tip_100_l_sem, '--b' , 'LineWidth', 1)
    plot(inds_span,VID_v_tip_100_l_mean-VID_v_tip_100_l_sem, '--b','LineWidth', 1)
end
ylabel('mm/s')
ylim([-400 400])
xlim([min(inds_span)-1 max(inds_span)+1])
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'vmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmax);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'vmin')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmin);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
xlabel(xlabel_text_raster_);
title('Lick velocity')
%% angle_tip
VID_angle_tip_100_grooming = raster_data.VID_angle_tip_100_grooming;
VID_angle_tip_100_r = raster_data.VID_angle_tip_100_r;

VID_angle_tip_100_grooming_mean = nanmean(VID_angle_tip_100_grooming);
VID_angle_tip_100_grooming_sem = nanstd(VID_angle_tip_100_grooming)/sqrt(size(VID_angle_tip_100_grooming,1));
VID_angle_tip_100_r_mean = nanmean(VID_angle_tip_100_r);
VID_angle_tip_100_r_sem = nanstd(VID_angle_tip_100_r)/sqrt(size(VID_angle_tip_100_r,1));

if (session_type  == 2)
    VID_angle_tip_100_l = raster_data.VID_angle_tip_100_l;
    VID_angle_tip_100_l_mean = nanmean(VID_angle_tip_100_l);
    VID_angle_tip_100_l_sem = nanstd(VID_angle_tip_100_l)/sqrt(size(VID_angle_tip_100_l,1));
end

subplot(4,9,3)
hold on
plot(inds_span,VID_angle_tip_100_grooming_mean,'g' , 'LineWidth', 2)
plot(inds_span,VID_angle_tip_100_grooming_mean+VID_angle_tip_100_grooming_sem, '--g' , 'LineWidth', 1)
plot(inds_span,VID_angle_tip_100_grooming_mean-VID_angle_tip_100_grooming_sem, '--g' , 'LineWidth', 1)

plot(inds_span,VID_angle_tip_100_r_mean,'r' , 'LineWidth', 2)
plot(inds_span,VID_angle_tip_100_r_mean+VID_angle_tip_100_r_sem, '--r', 'LineWidth', 1)
plot(inds_span,VID_angle_tip_100_r_mean-VID_angle_tip_100_r_sem, '--r', 'LineWidth', 1)

if (session_type  == 2)
    plot(inds_span,VID_angle_tip_100_l_mean ,'b' , 'LineWidth', 2)
    plot(inds_span,VID_angle_tip_100_l_mean+VID_angle_tip_100_l_sem, '--b' , 'LineWidth', 1)
    plot(inds_span,VID_angle_tip_100_l_mean-VID_angle_tip_100_l_sem, '--b','LineWidth', 1)
end
ylabel('deg')
ylim([0 90])
xlim([min(inds_span)-1 max(inds_span)+1])
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'vmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmax);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'vmin')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmin);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
xlabel(xlabel_text_raster_);
title('Lick angle')

%% CS probability bar plots

% subplot(4,9,[25 26 27])
% hold on;
% if (session_type  == 2)
%     x_axis = categorical({'L', 'G', 'R', 'All'});
%     x_axis = reordercats(x_axis,{'L', 'G', 'R', 'All'});
% else
%     x_axis = categorical({'G', 'R', 'All', 'Mean'});
%     x_axis = reordercats(x_axis,{'G', 'R', 'All'});
% end
%
% if (session_type  == 2)
%     bar(x_axis(1), prob_amplitude(1), 'k')
%     bar(x_axis(2), prob_amplitude(2), 'k')
%     bar(x_axis(3), prob_amplitude(3), 'k')
%     bar(x_axis(4), prob_amplitude(4), 'k')
% else
%     bar(x_axis(1), prob_amplitude(1), 'k')
%     bar(x_axis(2), prob_amplitude(2), 'k')
%     bar(x_axis(3), prob_amplitude(3), 'k')
%
% end
% kinematic_title = char(xlabel_text_raster_);
% title(['Prob CS [-50 50]ms | ' kinematic_title(1:end-5)])
% ylabel('Probability')
% if (max(prob_amplitude) > 0)
%     ylim([0 max(prob_amplitude) + max(prob_amplitude)/10])
% else
%     ylim([0 0.1])
%
% end

% if (session_type  == 2)
%     subplot(4,9,25)
%     x_axis = categorical({'L_1', 'L_0',});
%     x_axis = reordercats(x_axis,{'L_1', 'L_0'});
%     bar(x_axis, prob_amplitude_l, 'k' )
%     title('Left: L_1:rew | L_0:no rew')
%     ylabel(xlabel_text_CS_probab_)
%     ylim([0 0.4])
% end


% subplot(4,9,27)
% x_axis = categorical({'R_1', 'R_0',});
% x_axis = reordercats(x_axis,{'R_1', 'R_0'});
% bar(x_axis, prob_amplitude_r, 'k' )
% ylabel(xlabel_text_CS_probab_)
% title('Right: R_1:rew | R_0:no rew')
% ylim([0 0.4])

%% Frequency and phase analysis
Fs = 100;

%%% Grooming %%%
train_data_logic_SS_ = train_data_logic_SS_grooming_corr;
train_data_logic_CS_ = train_data_logic_CS_grooming_corr;
d_tip_ = raster_data.VID_d_tip_100_grooming_corr;

mean_d_tip_ = nanmean(d_tip_);
se_d_tip_ = nanstd(d_tip_)/sqrt(length(d_tip_));
mean_SS_ = ESN_smooth(nanmean(train_data_logic_SS_)*100,5);
se_SS_ = ESN_smooth(nanstd(train_data_logic_SS_)/sqrt(length(train_data_logic_SS_))*100,5);
mean_CS_ = ESN_smooth(nanmean(train_data_logic_CS_)*100,10);
se_CS_ = ESN_smooth(nanstd(train_data_logic_CS_)/sqrt(length(train_data_logic_CS_))*100,10);

phase_diff_ = wrapTo180(rad2deg(phdiffmeasure(mean_d_tip_,mean_SS_)));


subplot(4,9,23)
N = length(mean_d_tip_);
freq = 0:Fs/length(mean_SS_):Fs/2;
mean_d_tip_subtracted = mean_d_tip_- mean(mean_d_tip_);
xdft_d_tip_ = fft(mean_d_tip_subtracted);
xdft_d_tip_ = xdft_d_tip_(1:N/2+1);
psdx_d_tip_ = ESN_smooth((1/(Fs*N)) * abs(xdft_d_tip_).^2,3);
psdx_d_tip_(2:end-1) = 2*psdx_d_tip_(2:end-1);
[max_d_tip_freq,ind_max_d_tip_freq] = max(10*log10(psdx_d_tip_(2:end)));
freq_d_tip = freq(ind_max_d_tip_freq+1);
mean_SS_subtracted = mean_SS_ - mean(mean_SS_);
xdft_SS_ = fft(mean_SS_subtracted);
xdft_SS_ = xdft_SS_(1:N/2+1);
psdx_SS_ = ESN_smooth((1/(Fs*N)) * abs(xdft_SS_).^2,3);
psdx_SS_(2:end-1) = 2*psdx_SS_(2:end-1);
[max_SS_freq,ind_max_SS_freq] = max(10*log10(psdx_SS_(2:end)));
freq_SS = freq(ind_max_SS_freq+1);
hold on
plot(freq,10*log10(psdx_d_tip_), 'Color', [0.7 0.7 0.7], 'LineWidth', 2)
plot(freq_d_tip, max_d_tip_freq,'o', 'Color', [0.7 0.7 0.7], 'LineWidth', 2)
plot(freq,10*log10(psdx_SS_), 'Color', [0 0 0], 'LineWidth', 2)
plot(freq_SS, max_SS_freq, 'o', 'Color', [0 0 0], 'LineWidth', 2)
xlabel('Freq (Hz)')
xlim([0 10])
ylim([-40 40])
ylabel('Pow/Freq (dB/Hz)')
if  abs(freq_SS - freq_d_tip) <= 0.5
    title(['R(SS:' num2str(freq_SS) ', Lick: ' num2str(freq_d_tip) ') | G'])
else
    title(['NR(SS:' num2str(freq_SS) ', Lick: ' num2str(freq_d_tip) ') | G'])
end

subplot(4,9,26)
hold on
yyaxis left
plot(inds_span_corr, mean_d_tip_, '-' ,'Color', [0.7 0.7 0.7] , 'LineWidth', 2)
plot(inds_span_corr, mean_d_tip_ + se_d_tip_, '--', 'Color', [0.7 0.7 0.7] , 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_ - se_d_tip_, '--' , 'Color', [0.7 0.7 0.7] , 'LineWidth', 1)
ylabel('Disp. (mm)')
ylim([0 20])
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
plot(inds_span_corr, mean_SS_, '-' , 'Color', [0 0 0] , 'LineWidth', 2)
plot(inds_span_corr, mean_SS_ + se_SS_ , '-' , 'Color', [0 0 0] , 'LineWidth', 1)
plot(inds_span_corr, mean_SS_ - se_SS_ , '-' , 'Color', [0 0 0] , 'LineWidth', 1)
ylabel('Disp. (mm)')
set(gca, 'YColor', 'k')
ylabel('SS Firing(spks/s)')
xlabel(xlabel_text_raster_bout_);
ylim([0 300])
if  abs(freq_SS - freq_d_tip) > 0.5
    title('N/A | G' )
elseif abs(freq_SS - freq_d_tip) <= 0.5 && phase_diff_ >= -30 && phase_diff_ <= 30
    title(['in-phase(' num2str(round(phase_diff_)) 'deg) | G'])
elseif abs(freq_SS - freq_d_tip) <= 0.5 &&(phase_diff_ > 150 || phase_diff_ < -150)
    title(['anti-phase(' num2str(round(phase_diff_)) 'deg) | G'])
elseif abs(freq_SS - freq_d_tip) <= 0.5 && (phase_diff_ > 30 && phase_diff_ <= 150)
    title(['lead(' num2str(round(phase_diff_)) 'deg) | G'])
elseif abs(freq_SS - freq_d_tip) <= 0.5 && (phase_diff_ <-30 && phase_diff_ >= -150)
    title(['lag(' num2str(round(phase_diff_)) 'deg) | G'])
end
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'vmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmax);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'vmin')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmin);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end


%%% R %%%
train_data_logic_SS_ = train_data_logic_SS_r_corr;
train_data_logic_CS_ = train_data_logic_CS_r_corr;
d_tip_ = raster_data.VID_d_tip_100_r_corr;

mean_d_tip_ = nanmean(d_tip_);
se_d_tip_ = nanstd(d_tip_)/sqrt(length(d_tip_));
mean_SS_ = ESN_smooth(nanmean(train_data_logic_SS_)*100,5);
se_SS_ = ESN_smooth(nanstd(train_data_logic_SS_)/sqrt(length(train_data_logic_SS_))*100,5);
mean_CS_ = ESN_smooth(nanmean(train_data_logic_CS_)*100,10);
se_CS_ = ESN_smooth(nanstd(train_data_logic_CS_)/sqrt(length(train_data_logic_CS_))*100,10);

phase_diff_ = wrapTo180(rad2deg(phdiffmeasure(mean_d_tip_,mean_SS_)));

subplot(4,9,24)
N = length(mean_d_tip_);
freq = 0:Fs/length(mean_SS_):Fs/2;
mean_d_tip_subtracted = mean_d_tip_- mean(mean_d_tip_);
xdft_d_tip_ = fft(mean_d_tip_subtracted);
xdft_d_tip_ = xdft_d_tip_(1:N/2+1);
psdx_d_tip_ = ESN_smooth((1/(Fs*N)) * abs(xdft_d_tip_).^2,3);
psdx_d_tip_(2:end-1) = 2*psdx_d_tip_(2:end-1);
[max_d_tip_freq,ind_max_d_tip_freq] = max(10*log10(psdx_d_tip_(2:end)));
freq_d_tip = freq(ind_max_d_tip_freq+1);
mean_SS_subtracted = mean_SS_ - mean(mean_SS_);
xdft_SS_ = fft(mean_SS_subtracted);
xdft_SS_ = xdft_SS_(1:N/2+1);
psdx_SS_ = ESN_smooth((1/(Fs*N)) * abs(xdft_SS_).^2,3);
psdx_SS_(2:end-1) = 2*psdx_SS_(2:end-1);
[max_SS_freq,ind_max_SS_freq] = max(10*log10(psdx_SS_(2:end)));
freq_SS = freq(ind_max_SS_freq+1);
hold on
plot(freq,10*log10(psdx_d_tip_), 'Color', [0.7 0.7 0.7], 'LineWidth', 2)
plot(freq_d_tip, max_d_tip_freq,'o', 'Color', [0.7 0.7 0.7], 'LineWidth', 2)
plot(freq,10*log10(psdx_SS_), 'Color', [0 0 0], 'LineWidth', 2)
plot(freq_SS, max_SS_freq, 'o', 'Color', [0 0 0], 'LineWidth', 2)
xlabel('Freq (Hz)')
xlim([0 10])
ylim([-40 40])
ylabel('Pow/Freq (dB/Hz)')
if  abs(freq_SS - freq_d_tip) <= 0.5
    title(['R(SS:' num2str(freq_SS) ', Lick: ' num2str(freq_d_tip) ') | R'])
else
    title(['NR(SS:' num2str(freq_SS) ', Lick: ' num2str(freq_d_tip) ') | R'])
end

subplot(4,9,27)
hold on
yyaxis left
plot(inds_span_corr, mean_d_tip_, '-' ,'Color', [0.7 0.7 0.7] , 'LineWidth', 2)
plot(inds_span_corr, mean_d_tip_ + se_d_tip_, '--', 'Color', [0.7 0.7 0.7] , 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_ - se_d_tip_, '--' , 'Color', [0.7 0.7 0.7] , 'LineWidth', 1)
ylabel('Disp. (mm)')
ylim([0 20])
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
plot(inds_span_corr, mean_SS_, '-' , 'Color', [0 0 0] , 'LineWidth', 2)
plot(inds_span_corr, mean_SS_ + se_SS_ , '-' , 'Color', [0 0 0] , 'LineWidth', 1)
plot(inds_span_corr, mean_SS_ - se_SS_ , '-' , 'Color', [0 0 0] , 'LineWidth', 1)
ylabel('Disp. (mm)')
set(gca, 'YColor', 'k')
ylabel('SS Firing(spks/s)')
xlabel(xlabel_text_raster_bout_);
ylim([0 300])
if  abs(freq_SS - freq_d_tip) > 0.5
    title('N/A | R' )
elseif abs(freq_SS - freq_d_tip) <= 0.5 && phase_diff_ >= -30 && phase_diff_ <= 30
    title(['in-phase(' num2str(round(phase_diff_)) 'deg) | R'])
elseif abs(freq_SS - freq_d_tip) <= 0.5 &&(phase_diff_ > 150 || phase_diff_ < -150)
    title(['anti-phase(' num2str(round(phase_diff_)) 'deg) | R'])
elseif abs(freq_SS - freq_d_tip) <= 0.5 && (phase_diff_ > 30 && phase_diff_ <= 150)
    title(['lead(' num2str(round(phase_diff_)) 'deg) | R'])
elseif abs(freq_SS - freq_d_tip) <= 0.5 && (phase_diff_ <-30 && phase_diff_ >= -150)
    title(['lag(' num2str(round(phase_diff_)) 'deg) | R'])
end
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'vmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmax);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'vmin')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmin);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end


%%% L %%%
train_data_logic_SS_ = train_data_logic_SS_l_corr;
train_data_logic_CS_ = train_data_logic_CS_l_corr;
d_tip_ = raster_data.VID_d_tip_100_l_corr;

mean_d_tip_ = nanmean(d_tip_);
se_d_tip_ = nanstd(d_tip_)/sqrt(length(d_tip_));
mean_SS_ = ESN_smooth(nanmean(train_data_logic_SS_)*100,5);
se_SS_ = ESN_smooth(nanstd(train_data_logic_SS_)/sqrt(length(train_data_logic_SS_))*100,5);
mean_CS_ = ESN_smooth(nanmean(train_data_logic_CS_)*100,10);
se_CS_ = ESN_smooth(nanstd(train_data_logic_CS_)/sqrt(length(train_data_logic_CS_))*100,10);

phase_diff_ = wrapTo180(rad2deg(phdiffmeasure(mean_d_tip_,mean_SS_)));

subplot(4,9,22)
N = length(mean_d_tip_);
freq = 0:Fs/length(mean_SS_):Fs/2;
mean_d_tip_subtracted = mean_d_tip_- mean(mean_d_tip_);
xdft_d_tip_ = fft(mean_d_tip_subtracted);
xdft_d_tip_ = xdft_d_tip_(1:N/2+1);
psdx_d_tip_ = ESN_smooth((1/(Fs*N)) * abs(xdft_d_tip_).^2,3);
psdx_d_tip_(2:end-1) = 2*psdx_d_tip_(2:end-1);
[max_d_tip_freq,ind_max_d_tip_freq] = max(10*log10(psdx_d_tip_(2:end)));
freq_d_tip = freq(ind_max_d_tip_freq+1);
mean_SS_subtracted = mean_SS_ - mean(mean_SS_);
xdft_SS_ = fft(mean_SS_subtracted);
xdft_SS_ = xdft_SS_(1:N/2+1);
psdx_SS_ = ESN_smooth((1/(Fs*N)) * abs(xdft_SS_).^2,3);
psdx_SS_(2:end-1) = 2*psdx_SS_(2:end-1);
[max_SS_freq,ind_max_SS_freq] = max(10*log10(psdx_SS_(2:end)));
freq_SS = freq(ind_max_SS_freq+1);
hold on
plot(freq,10*log10(psdx_d_tip_), 'Color', [0.7 0.7 0.7], 'LineWidth', 2)
plot(freq_d_tip, max_d_tip_freq,'o', 'Color', [0.7 0.7 0.7], 'LineWidth', 2)
plot(freq,10*log10(psdx_SS_), 'Color', [0 0 0], 'LineWidth', 2)
plot(freq_SS, max_SS_freq, 'o', 'Color', [0 0 0], 'LineWidth', 2)
xlabel('Freq (Hz)')
xlim([0 10])
ylim([-40 40])
ylabel('Pow/Freq (dB/Hz)')
if  abs(freq_SS - freq_d_tip) <= 0.5
    title(['R(SS:' num2str(freq_SS) ', Lick: ' num2str(freq_d_tip) ') | L'])
else
    title(['NR(SS:' num2str(freq_SS) ', Lick: ' num2str(freq_d_tip) ') | L'])
end

subplot(4,9,25)
hold on
yyaxis left
plot(inds_span_corr, mean_d_tip_, '-' ,'Color', [0.7 0.7 0.7] , 'LineWidth', 2)
plot(inds_span_corr, mean_d_tip_ + se_d_tip_, '--', 'Color', [0.7 0.7 0.7] , 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_ - se_d_tip_, '--' , 'Color', [0.7 0.7 0.7] , 'LineWidth', 1)
ylabel('Disp. (mm)')
ylim([0 20])
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
plot(inds_span_corr, mean_SS_, '-' , 'Color', [0 0 0] , 'LineWidth', 2)
plot(inds_span_corr, mean_SS_ + se_SS_ , '-' , 'Color', [0 0 0] , 'LineWidth', 1)
plot(inds_span_corr, mean_SS_ - se_SS_ , '-' , 'Color', [0 0 0] , 'LineWidth', 1)
ylabel('Disp. (mm)')
set(gca, 'YColor', 'k')
ylabel('SS Firing(spks/s)')
xlabel(xlabel_text_raster_bout_);
ylim([0 300])
if  abs(freq_SS - freq_d_tip) > 0.5
    title('N/A | L' )
elseif ~isempty(freq_SS) && abs(freq_SS - freq_d_tip) <= 0.5 && phase_diff_ >= -30 && phase_diff_ <= 30
    title(['in-phase(' num2str(round(phase_diff_)) 'deg) | L'])
elseif ~isempty(freq_SS) && abs(freq_SS - freq_d_tip) <= 0.5 &&(phase_diff_ > 150 || phase_diff_ < -150)
    title(['anti-phase(' num2str(round(phase_diff_)) 'deg) | L'])
elseif ~isempty(freq_SS) && abs(freq_SS - freq_d_tip) <= 0.5 && (phase_diff_ > 30 && phase_diff_ <= 150)
    title(['lead(' num2str(round(phase_diff_)) 'deg) | L'])
elseif ~isempty(freq_SS) && abs(freq_SS - freq_d_tip) <= 0.5 && (phase_diff_ <-30 && phase_diff_ >= -150)
    title(['lag(' num2str(round(phase_diff_)) 'deg) | L'])
end
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'vmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmax);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'vmin')
    xline(0,'LineWidth', 2, 'Color', color_lick_vmin);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end

% % Response thresholds
% %all
% SS_corr_mean_all = ESN_smooth(nanmean(train_data_logic_SS_all_corr)*100 - SS_baseline_mean,5);
% SS_corr_sd_all = nanstd(train_data_logic_SS_all_corr)*100;
% SS_corr_mean_threshold = nanmean(SS_corr_mean_all);
% SS_corr_sd_threshold = nanstd(SS_corr_mean_all);
%
% % Response vs threshold for strong rhythmic modulation
% SS_corr_all_max_pass_bool = SS_corr_mean_all > (SS_corr_mean_threshold + SS_corr_sd_threshold);
% SS_corr_all_min_pass_bool = SS_corr_mean_all < (SS_corr_mean_threshold - SS_corr_sd_threshold);
%
% % String response
% if sum(SS_corr_all_max_pass_bool) > 0 || sum(SS_corr_all_min_pass_bool) > 0
%     response = 'Rhythmic';
% else
%     response = 'Non Rhythmic';
% end
%
% subplot(4,9,[22 23 24])
% hold on;
% plot(inds_span_corr, SS_corr_mean_all, 'k', 'LineWIdth', 2);
% plot(inds_span_corr(SS_corr_all_max_pass_bool), SS_corr_mean_all(SS_corr_all_max_pass_bool), '*r')
% plot(inds_span_corr(SS_corr_all_min_pass_bool), SS_corr_mean_all(SS_corr_all_min_pass_bool), '*r')
% % plot(inds_span_corr, SS_corr_mean_grooming, 'g', 'LineWIdth', 2);
% % plot(inds_span_corr, SS_corr_mean_r, 'r', 'LineWIdth', 2);
% % plot(inds_span_corr, SS_corr_mean_l, 'b', 'LineWIdth', 2);
% yline(SS_corr_mean_threshold,'--m', 'LineWidth', 1.5);
% yline(SS_corr_mean_threshold + SS_corr_sd_threshold,'m', 'LineWidth', 1.5);
% yline(SS_corr_mean_threshold - SS_corr_sd_threshold,'m', 'LineWidth', 1.5);
% ylabel('Change in SS firing from baseline (spk/s)')
% ylim([(min(SS_corr_mean_all(SS_corr_all_min_pass_bool)) - 20) (max(SS_corr_mean_all(SS_corr_all_max_pass_bool)) + 20)])
% xlim([-1 1])
% if contains(xlabel_text_raster_, 'onset')
%     xline(0,'LineWidth', 2, 'Color', color_lick_onset);
% elseif contains(xlabel_text_raster_, 'dmax')
%     xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
% elseif contains(xlabel_text_raster_, 'offset')
%     xline(0,'LineWidth', 2, 'Color', color_lick_offset);
% end
% xlabel(xlabel_text_raster_bout_);
% title(['Rhymicity test: ' response ' | Change in SS | All licks'])
% set(gca, 'YColor', 'k')

%ESN_Beautify_Plot
end

%% function plot_rasters_data_kin
function fig_handle_ = plot_rasters_data_kin(raster_data, plot_data, session_type)
%%
fig_num_               = plot_data.fig_num_;
xlabel_text_raster_    = plot_data.xlabel_text_raster_;
xlabel_text_raster_bout_   = plot_data.xlabel_text_raster_bout_;

xlabel_text_CS_probab_ = plot_data.xlabel_text_CS_probab_;
inds_span              = plot_data.inds_span * 10;

% % plot xlim and ylim
range_SS_Firing = [0 300];

Line_Color = lines(7);
color_SS = Line_Color(1,:);
color_CS = Line_Color(7,:);
color_lick_onset = Line_Color(3,:);
color_lick_vmax = Line_Color(4,:);
color_lick_dmax = Line_Color(5,:);
color_lick_vmin = Line_Color(2,:);
color_lick_offset = Line_Color(6,:);

% 80toINF licks
train_data_logic_SS_80toINF_0to5 = raster_data.train_data_logic_SS_80toINF_0to5;
train_data_logic_SS_80toINF_5to10 = raster_data.train_data_logic_SS_80toINF_5to10;
train_data_logic_SS_80toINF_10to15 = raster_data.train_data_logic_SS_80toINF_10to15;
train_data_logic_SS_80toINF_15toINF= raster_data.train_data_logic_SS_80toINF_15toINF;
train_data_logic_CS_80toINF_0to5 = raster_data.train_data_logic_CS_80toINF_0to5;
train_data_logic_CS_80toINF_5to10 = raster_data.train_data_logic_CS_80toINF_5to10;
train_data_logic_CS_80toINF_10to15 = raster_data.train_data_logic_CS_80toINF_10to15;
train_data_logic_CS_80toINF_15toINF = raster_data.train_data_logic_CS_80toINF_15toINF;
VID_d_tip_100_80toINF_0to5 = raster_data.VID_d_tip_100_80toINF_0to5;
VID_d_tip_100_80toINF_5to10 = raster_data.VID_d_tip_100_80toINF_5to10;
VID_d_tip_100_80toINF_10to15 = raster_data.VID_d_tip_100_80toINF_10to15;
VID_d_tip_100_80toINF_15toINF = raster_data.VID_d_tip_100_80toINF_15toINF;
VID_v_tip_100_80toINF_0to5 = raster_data.VID_v_tip_100_80toINF_0to5;
VID_v_tip_100_80toINF_5to10 = raster_data.VID_v_tip_100_80toINF_5to10;
VID_v_tip_100_80toINF_10to15 = raster_data.VID_v_tip_100_80toINF_10to15;
VID_v_tip_100_80toINF_15toINF = raster_data.VID_v_tip_100_80toINF_15toINF;
VID_angle_tip_100_80toINF_0to5 = raster_data.VID_angle_tip_100_80toINF_0to5;
VID_angle_tip_100_80toINF_5to10 = raster_data.VID_angle_tip_100_80toINF_5to10;
VID_angle_tip_100_80toINF_10to15 = raster_data.VID_angle_tip_100_80toINF_10to15;
VID_angle_tip_100_80toINF_15toINF = raster_data.VID_angle_tip_100_80toINF_15toINF;

% 60to80 licks
train_data_logic_SS_60to80_0to5 = raster_data.train_data_logic_SS_60to80_0to5;
train_data_logic_SS_60to80_5to10 = raster_data.train_data_logic_SS_60to80_5to10;
train_data_logic_SS_60to80_10to15 = raster_data.train_data_logic_SS_60to80_10to15;
train_data_logic_SS_60to80_15toINF= raster_data.train_data_logic_SS_60to80_15toINF;
train_data_logic_CS_60to80_0to5 = raster_data.train_data_logic_CS_60to80_0to5;
train_data_logic_CS_60to80_5to10 = raster_data.train_data_logic_CS_60to80_5to10;
train_data_logic_CS_60to80_10to15 = raster_data.train_data_logic_CS_60to80_10to15;
train_data_logic_CS_60to80_15toINF = raster_data.train_data_logic_CS_60to80_15toINF;
VID_d_tip_100_60to80_0to5 = raster_data.VID_d_tip_100_60to80_0to5;
VID_d_tip_100_60to80_5to10 = raster_data.VID_d_tip_100_60to80_5to10;
VID_d_tip_100_60to80_10to15 = raster_data.VID_d_tip_100_60to80_10to15;
VID_d_tip_100_60to80_15toINF = raster_data.VID_d_tip_100_60to80_15toINF;
VID_v_tip_100_60to80_0to5 = raster_data.VID_v_tip_100_60to80_0to5;
VID_v_tip_100_60to80_5to10 = raster_data.VID_v_tip_100_60to80_5to10;
VID_v_tip_100_60to80_10to15 = raster_data.VID_v_tip_100_60to80_10to15;
VID_v_tip_100_60to80_15toINF = raster_data.VID_v_tip_100_60to80_15toINF;
VID_angle_tip_100_60to80_0to5 = raster_data.VID_angle_tip_100_60to80_0to5;
VID_angle_tip_100_60to80_5to10 = raster_data.VID_angle_tip_100_60to80_5to10;
VID_angle_tip_100_60to80_10to15 = raster_data.VID_angle_tip_100_60to80_10to15;
VID_angle_tip_100_60to80_15toINF = raster_data.VID_angle_tip_100_60to80_15toINF;

% 20to60 licks
train_data_logic_SS_20to60_0to5 = raster_data.train_data_logic_SS_20to60_0to5;
train_data_logic_SS_20to60_5to10 = raster_data.train_data_logic_SS_20to60_5to10;
train_data_logic_SS_20to60_10to15 = raster_data.train_data_logic_SS_20to60_10to15;
train_data_logic_SS_20to60_15toINF= raster_data.train_data_logic_SS_20to60_15toINF;
train_data_logic_CS_20to60_0to5 = raster_data.train_data_logic_CS_20to60_0to5;
train_data_logic_CS_20to60_5to10 = raster_data.train_data_logic_CS_20to60_5to10;
train_data_logic_CS_20to60_10to15 = raster_data.train_data_logic_CS_20to60_10to15;
train_data_logic_CS_20to60_15toINF = raster_data.train_data_logic_CS_20to60_15toINF;
VID_d_tip_100_20to60_0to5 = raster_data.VID_d_tip_100_20to60_0to5;
VID_d_tip_100_20to60_5to10 = raster_data.VID_d_tip_100_20to60_5to10;
VID_d_tip_100_20to60_10to15 = raster_data.VID_d_tip_100_20to60_10to15;
VID_d_tip_100_20to60_15toINF = raster_data.VID_d_tip_100_20to60_15toINF;
VID_v_tip_100_20to60_0to5 = raster_data.VID_v_tip_100_20to60_0to5;
VID_v_tip_100_20to60_5to10 = raster_data.VID_v_tip_100_20to60_5to10;
VID_v_tip_100_20to60_10to15 = raster_data.VID_v_tip_100_20to60_10to15;
VID_v_tip_100_20to60_15toINF = raster_data.VID_v_tip_100_20to60_15toINF;
VID_angle_tip_100_20to60_0to5 = raster_data.VID_angle_tip_100_20to60_0to5;
VID_angle_tip_100_20to60_5to10 = raster_data.VID_angle_tip_100_20to60_5to10;
VID_angle_tip_100_20to60_10to15 = raster_data.VID_angle_tip_100_20to60_10to15;
VID_angle_tip_100_20to60_15toINF = raster_data.VID_angle_tip_100_20to60_15toINF;


% neg20to20 licks
train_data_logic_SS_neg20to20_0to5 = raster_data.train_data_logic_SS_neg20to20_0to5;
train_data_logic_SS_neg20to20_5to10 = raster_data.train_data_logic_SS_neg20to20_5to10;
train_data_logic_SS_neg20to20_10to15 = raster_data.train_data_logic_SS_neg20to20_10to15;
train_data_logic_SS_neg20to20_15toINF= raster_data.train_data_logic_SS_neg20to20_15toINF;
train_data_logic_CS_neg20to20_0to5 = raster_data.train_data_logic_CS_neg20to20_0to5;
train_data_logic_CS_neg20to20_5to10 = raster_data.train_data_logic_CS_neg20to20_5to10;
train_data_logic_CS_neg20to20_10to15 = raster_data.train_data_logic_CS_neg20to20_10to15;
train_data_logic_CS_neg20to20_15toINF = raster_data.train_data_logic_CS_neg20to20_15toINF;
VID_d_tip_100_neg20to20_0to5 = raster_data.VID_d_tip_100_neg20to20_0to5;
VID_d_tip_100_neg20to20_5to10 = raster_data.VID_d_tip_100_neg20to20_5to10;
VID_d_tip_100_neg20to20_10to15 = raster_data.VID_d_tip_100_neg20to20_10to15;
VID_d_tip_100_neg20to20_15toINF = raster_data.VID_d_tip_100_neg20to20_15toINF;
VID_v_tip_100_neg20to20_0to5 = raster_data.VID_v_tip_100_neg20to20_0to5;
VID_v_tip_100_neg20to20_5to10 = raster_data.VID_v_tip_100_neg20to20_5to10;
VID_v_tip_100_neg20to20_10to15 = raster_data.VID_v_tip_100_neg20to20_10to15;
VID_v_tip_100_neg20to20_15toINF = raster_data.VID_v_tip_100_neg20to20_15toINF;
VID_angle_tip_100_neg20to20_0to5 = raster_data.VID_angle_tip_100_neg20to20_0to5;
VID_angle_tip_100_neg20to20_5to10 = raster_data.VID_angle_tip_100_neg20to20_5to10;
VID_angle_tip_100_neg20to20_10to15 = raster_data.VID_angle_tip_100_neg20to20_10to15;
VID_angle_tip_100_neg20to20_15toINF = raster_data.VID_angle_tip_100_neg20to20_15toINF;

% neg20toneg60 licks
train_data_logic_SS_neg20toneg60_0to5 = raster_data.train_data_logic_SS_neg20toneg60_0to5;
train_data_logic_SS_neg20toneg60_5to10 = raster_data.train_data_logic_SS_neg20toneg60_5to10;
train_data_logic_SS_neg20toneg60_10to15 = raster_data.train_data_logic_SS_neg20toneg60_10to15;
train_data_logic_SS_neg20toneg60_15toINF= raster_data.train_data_logic_SS_neg20toneg60_15toINF;
train_data_logic_CS_neg20toneg60_0to5 = raster_data.train_data_logic_CS_neg20toneg60_0to5;
train_data_logic_CS_neg20toneg60_5to10 = raster_data.train_data_logic_CS_neg20toneg60_5to10;
train_data_logic_CS_neg20toneg60_10to15 = raster_data.train_data_logic_CS_neg20toneg60_10to15;
train_data_logic_CS_neg20toneg60_15toINF = raster_data.train_data_logic_CS_neg20toneg60_15toINF;
VID_d_tip_100_neg20toneg60_0to5 = raster_data.VID_d_tip_100_neg20toneg60_0to5;
VID_d_tip_100_neg20toneg60_5to10 = raster_data.VID_d_tip_100_neg20toneg60_5to10;
VID_d_tip_100_neg20toneg60_10to15 = raster_data.VID_d_tip_100_neg20toneg60_10to15;
VID_d_tip_100_neg20toneg60_15toINF = raster_data.VID_d_tip_100_neg20toneg60_15toINF;
VID_v_tip_100_neg20toneg60_0to5 = raster_data.VID_v_tip_100_neg20toneg60_0to5;
VID_v_tip_100_neg20toneg60_5to10 = raster_data.VID_v_tip_100_neg20toneg60_5to10;
VID_v_tip_100_neg20toneg60_10to15 = raster_data.VID_v_tip_100_neg20toneg60_10to15;
VID_v_tip_100_neg20toneg60_15toINF = raster_data.VID_v_tip_100_neg20toneg60_15toINF;
VID_angle_tip_100_neg20toneg60_0to5 = raster_data.VID_angle_tip_100_neg20toneg60_0to5;
VID_angle_tip_100_neg20toneg60_5to10 = raster_data.VID_angle_tip_100_neg20toneg60_5to10;
VID_angle_tip_100_neg20toneg60_10to15 = raster_data.VID_angle_tip_100_neg20toneg60_10to15;
VID_angle_tip_100_neg20toneg60_15toINF = raster_data.VID_angle_tip_100_neg20toneg60_15toINF;

% neg60toneg80 licks
train_data_logic_SS_neg60toneg80_0to5 = raster_data.train_data_logic_SS_neg60toneg80_0to5;
train_data_logic_SS_neg60toneg80_5to10 = raster_data.train_data_logic_SS_neg60toneg80_5to10;
train_data_logic_SS_neg60toneg80_10to15 = raster_data.train_data_logic_SS_neg60toneg80_10to15;
train_data_logic_SS_neg60toneg80_15toINF= raster_data.train_data_logic_SS_neg60toneg80_15toINF;
train_data_logic_CS_neg60toneg80_0to5 = raster_data.train_data_logic_CS_neg60toneg80_0to5;
train_data_logic_CS_neg60toneg80_5to10 = raster_data.train_data_logic_CS_neg60toneg80_5to10;
train_data_logic_CS_neg60toneg80_10to15 = raster_data.train_data_logic_CS_neg60toneg80_10to15;
train_data_logic_CS_neg60toneg80_15toINF = raster_data.train_data_logic_CS_neg60toneg80_15toINF;
VID_d_tip_100_neg60toneg80_0to5 = raster_data.VID_d_tip_100_neg60toneg80_0to5;
VID_d_tip_100_neg60toneg80_5to10 = raster_data.VID_d_tip_100_neg60toneg80_5to10;
VID_d_tip_100_neg60toneg80_10to15 = raster_data.VID_d_tip_100_neg60toneg80_10to15;
VID_d_tip_100_neg60toneg80_15toINF = raster_data.VID_d_tip_100_neg60toneg80_15toINF;
VID_v_tip_100_neg60toneg80_0to5 = raster_data.VID_v_tip_100_neg60toneg80_0to5;
VID_v_tip_100_neg60toneg80_5to10 = raster_data.VID_v_tip_100_neg60toneg80_5to10;
VID_v_tip_100_neg60toneg80_10to15 = raster_data.VID_v_tip_100_neg60toneg80_10to15;
VID_v_tip_100_neg60toneg80_15toINF = raster_data.VID_v_tip_100_neg60toneg80_15toINF;
VID_angle_tip_100_neg60toneg80_0to5 = raster_data.VID_angle_tip_100_neg60toneg80_0to5;
VID_angle_tip_100_neg60toneg80_5to10 = raster_data.VID_angle_tip_100_neg60toneg80_5to10;
VID_angle_tip_100_neg60toneg80_10to15 = raster_data.VID_angle_tip_100_neg60toneg80_10to15;
VID_angle_tip_100_neg60toneg80_15toINF = raster_data.VID_angle_tip_100_neg60toneg80_15toINF;

% neg80tonegINF licks
train_data_logic_SS_neg80tonegINF_0to5 = raster_data.train_data_logic_SS_neg80tonegINF_0to5;
train_data_logic_SS_neg80tonegINF_5to10 = raster_data.train_data_logic_SS_neg80tonegINF_5to10;
train_data_logic_SS_neg80tonegINF_10to15 = raster_data.train_data_logic_SS_neg80tonegINF_10to15;
train_data_logic_SS_neg80tonegINF_15toINF= raster_data.train_data_logic_SS_neg80tonegINF_15toINF;
train_data_logic_CS_neg80tonegINF_0to5 = raster_data.train_data_logic_CS_neg80tonegINF_0to5;
train_data_logic_CS_neg80tonegINF_5to10 = raster_data.train_data_logic_CS_neg80tonegINF_5to10;
train_data_logic_CS_neg80tonegINF_10to15 = raster_data.train_data_logic_CS_neg80tonegINF_10to15;
train_data_logic_CS_neg80tonegINF_15toINF = raster_data.train_data_logic_CS_neg80tonegINF_15toINF;
VID_d_tip_100_neg80tonegINF_0to5 = raster_data.VID_d_tip_100_neg80tonegINF_0to5;
VID_d_tip_100_neg80tonegINF_5to10 = raster_data.VID_d_tip_100_neg80tonegINF_5to10;
VID_d_tip_100_neg80tonegINF_10to15 = raster_data.VID_d_tip_100_neg80tonegINF_10to15;
VID_d_tip_100_neg80tonegINF_15toINF = raster_data.VID_d_tip_100_neg80tonegINF_15toINF;
VID_v_tip_100_neg80tonegINF_0to5 = raster_data.VID_v_tip_100_neg80tonegINF_0to5;
VID_v_tip_100_neg80tonegINF_5to10 = raster_data.VID_v_tip_100_neg80tonegINF_5to10;
VID_v_tip_100_neg80tonegINF_10to15 = raster_data.VID_v_tip_100_neg80tonegINF_10to15;
VID_v_tip_100_neg80tonegINF_15toINF = raster_data.VID_v_tip_100_neg80tonegINF_15toINF;
VID_angle_tip_100_neg80tonegINF_0to5 = raster_data.VID_angle_tip_100_neg80tonegINF_0to5;
VID_angle_tip_100_neg80tonegINF_5to10 = raster_data.VID_angle_tip_100_neg80tonegINF_5to10;
VID_angle_tip_100_neg80tonegINF_10to15 = raster_data.VID_angle_tip_100_neg80tonegINF_10to15;
VID_angle_tip_100_neg80tonegINF_15toINF = raster_data.VID_angle_tip_100_neg80tonegINF_15toINF;

fig_handle_ = figure(fig_num_);
fig_handle_.WindowState = 'maximized';
clf(fig_handle_)

%% 80toINF licks
train_data_logic_SS_0to5 = train_data_logic_SS_80toINF_0to5;
train_data_logic_SS_5to10 = train_data_logic_SS_80toINF_5to10;
train_data_logic_SS_10to15 = train_data_logic_SS_80toINF_10to15;
train_data_logic_SS_15toINF = train_data_logic_SS_80toINF_15toINF;
train_data_logic_CS_0to5 = train_data_logic_CS_80toINF_0to5;
train_data_logic_CS_5to10 = train_data_logic_CS_80toINF_5to10;
train_data_logic_CS_10to15 = train_data_logic_CS_80toINF_10to15;
train_data_logic_CS_15toINF = train_data_logic_CS_80toINF_15toINF;
VID_d_tip_100_0to5 = VID_d_tip_100_80toINF_0to5;
VID_d_tip_100_5to10 = VID_d_tip_100_80toINF_5to10;
VID_d_tip_100_10to15 = VID_d_tip_100_80toINF_10to15;
VID_d_tip_100_15toINF = VID_d_tip_100_80toINF_15toINF;
VID_v_tip_100_0to5 = VID_v_tip_100_80toINF_0to5;
VID_v_tip_100_5to10 = VID_v_tip_100_80toINF_5to10;
VID_v_tip_100_10to15 = VID_v_tip_100_80toINF_10to15;
VID_v_tip_100_15toINF = VID_v_tip_100_80toINF_15toINF;
VID_angle_tip_100_0to5 = VID_angle_tip_100_80toINF_0to5;
VID_angle_tip_100_5to10 = VID_angle_tip_100_80toINF_5to10;
VID_angle_tip_100_10to15 = VID_angle_tip_100_80toINF_10to15;
VID_angle_tip_100_15toINF = VID_angle_tip_100_80toINF_15toINF;

subplot(4,7,7)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_0to5) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_0to5) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_0to5,1));
if ~isempty(train_data_logic_SS_0to5) && size(train_data_logic_SS_0to5,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['80toINF deg | 0to5 mm | n = ' num2str(size(train_data_logic_SS_0to5,1)) ])

subplot(4,7,14)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_5to10) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_5to10) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_5to10,1));
if ~isempty(train_data_logic_SS_5to10) && size(train_data_logic_SS_5to10,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['80toINF deg | 5to10 mm | n = ' num2str(size(train_data_logic_SS_5to10,1)) ])

subplot(4,7,21)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_10to15) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_10to15) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_10to15,1));
if ~isempty(train_data_logic_SS_10to15) && size(train_data_logic_SS_10to15,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['80toINF deg | 10to15 mm | n = ' num2str(size(train_data_logic_SS_10to15,1)) ])

subplot(4,7,28)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_15toINF) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_15toINF) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_15toINF,1));
if ~isempty(train_data_logic_SS_15toINF) && size(train_data_logic_SS_15toINF,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['80toINF deg | 15toINF mm | n = ' num2str(size(train_data_logic_SS_15toINF,1)) ])

%% 60to80 licks
train_data_logic_SS_0to5 = train_data_logic_SS_60to80_0to5;
train_data_logic_SS_5to10 = train_data_logic_SS_60to80_5to10;
train_data_logic_SS_10to15 = train_data_logic_SS_60to80_10to15;
train_data_logic_SS_15toINF = train_data_logic_SS_60to80_15toINF;
train_data_logic_CS_0to5 = train_data_logic_CS_60to80_0to5;
train_data_logic_CS_5to10 = train_data_logic_CS_60to80_5to10;
train_data_logic_CS_10to15 = train_data_logic_CS_60to80_10to15;
train_data_logic_CS_15toINF = train_data_logic_CS_60to80_15toINF;
VID_d_tip_100_0to5 = VID_d_tip_100_60to80_0to5;
VID_d_tip_100_5to10 = VID_d_tip_100_60to80_5to10;
VID_d_tip_100_10to15 = VID_d_tip_100_60to80_10to15;
VID_d_tip_100_15toINF = VID_d_tip_100_60to80_15toINF;
VID_v_tip_100_0to5 = VID_v_tip_100_60to80_0to5;
VID_v_tip_100_5to10 = VID_v_tip_100_60to80_5to10;
VID_v_tip_100_10to15 = VID_v_tip_100_60to80_10to15;
VID_v_tip_100_15toINF = VID_v_tip_100_60to80_15toINF;
VID_angle_tip_100_0to5 = VID_angle_tip_100_60to80_0to5;
VID_angle_tip_100_5to10 = VID_angle_tip_100_60to80_5to10;
VID_angle_tip_100_10to15 = VID_angle_tip_100_60to80_10to15;
VID_angle_tip_100_15toINF = VID_angle_tip_100_60to80_15toINF;

subplot(4,7,6)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_0to5) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_0to5) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_0to5,1));
if ~isempty(train_data_logic_SS_0to5) && size(train_data_logic_SS_0to5,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['60to80 deg | 0to5 mm | n = ' num2str(size(train_data_logic_SS_0to5,1)) ])

subplot(4,7,13)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_5to10) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_5to10) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_5to10,1));
if ~isempty(train_data_logic_SS_5to10) && size(train_data_logic_SS_5to10,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['60to80 deg | 5to10 mm | n = ' num2str(size(train_data_logic_SS_5to10,1)) ])

subplot(4,7,20)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_10to15) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_10to15) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_10to15,1));
if ~isempty(train_data_logic_SS_10to15) && size(train_data_logic_SS_10to15,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['60to80 deg | 10to15 mm | n = ' num2str(size(train_data_logic_SS_10to15,1)) ])

subplot(4,7,27)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_15toINF) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_15toINF) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_15toINF,1));
if ~isempty(train_data_logic_SS_15toINF) && size(train_data_logic_SS_15toINF,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_ , 5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['60to80 deg | 15toINF mm | n = ' num2str(size(train_data_logic_SS_15toINF,1)) ])

%% 20to60 licks
train_data_logic_SS_0to5 = train_data_logic_SS_20to60_0to5;
train_data_logic_SS_5to10 = train_data_logic_SS_20to60_5to10;
train_data_logic_SS_10to15 = train_data_logic_SS_20to60_10to15;
train_data_logic_SS_15toINF = train_data_logic_SS_20to60_15toINF;
train_data_logic_CS_0to5 = train_data_logic_CS_20to60_0to5;
train_data_logic_CS_5to10 = train_data_logic_CS_20to60_5to10;
train_data_logic_CS_10to15 = train_data_logic_CS_20to60_10to15;
train_data_logic_CS_15toINF = train_data_logic_CS_20to60_15toINF;
VID_d_tip_100_0to5 = VID_d_tip_100_20to60_0to5;
VID_d_tip_100_5to10 = VID_d_tip_100_20to60_5to10;
VID_d_tip_100_10to15 = VID_d_tip_100_20to60_10to15;
VID_d_tip_100_15toINF = VID_d_tip_100_20to60_15toINF;
VID_v_tip_100_0to5 = VID_v_tip_100_20to60_0to5;
VID_v_tip_100_5to10 = VID_v_tip_100_20to60_5to10;
VID_v_tip_100_10to15 = VID_v_tip_100_20to60_10to15;
VID_v_tip_100_15toINF = VID_v_tip_100_20to60_15toINF;
VID_angle_tip_100_0to5 = VID_angle_tip_100_20to60_0to5;
VID_angle_tip_100_5to10 = VID_angle_tip_100_20to60_5to10;
VID_angle_tip_100_10to15 = VID_angle_tip_100_20to60_10to15;
VID_angle_tip_100_15toINF = VID_angle_tip_100_20to60_15toINF;

subplot(4,7,5)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_0to5) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_0to5) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_0to5,1));
if ~isempty(train_data_logic_SS_0to5) && size(train_data_logic_SS_0to5,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
%plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['20to60 deg | 0to5 mm | n = ' num2str(size(train_data_logic_SS_0to5,1)) ])

subplot(4,7,12)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_5to10) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_5to10) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_5to10,1));
if ~isempty(train_data_logic_SS_5to10) && size(train_data_logic_SS_5to10,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['20to60 deg | 5to10 mm | n = ' num2str(size(train_data_logic_SS_5to10,1)) ])

subplot(4,7,19)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_10to15) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_10to15) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_10to15,1));
if ~isempty(train_data_logic_SS_10to15) && size(train_data_logic_SS_10to15,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['20to60 deg | 10to15 mm | n = ' num2str(size(train_data_logic_SS_10to15,1)) ])

subplot(4,7,26)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_15toINF) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_15toINF) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_15toINF,1));
if ~isempty(train_data_logic_SS_15toINF) && size(train_data_logic_SS_15toINF,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['20to60 deg | 15toINF mm | n = ' num2str(size(train_data_logic_SS_15toINF,1)) ])

%% neg20to20 licks
train_data_logic_SS_0to5 = train_data_logic_SS_neg20to20_0to5;
train_data_logic_SS_5to10 = train_data_logic_SS_neg20to20_5to10;
train_data_logic_SS_10to15 = train_data_logic_SS_neg20to20_10to15;
train_data_logic_SS_15toINF = train_data_logic_SS_neg20to20_15toINF;
train_data_logic_CS_0to5 = train_data_logic_CS_neg20to20_0to5;
train_data_logic_CS_5to10 = train_data_logic_CS_neg20to20_5to10;
train_data_logic_CS_10to15 = train_data_logic_CS_neg20to20_10to15;
train_data_logic_CS_15toINF = train_data_logic_CS_neg20to20_15toINF;
VID_d_tip_100_0to5 = VID_d_tip_100_neg20to20_0to5;
VID_d_tip_100_5to10 = VID_d_tip_100_neg20to20_5to10;
VID_d_tip_100_10to15 = VID_d_tip_100_neg20to20_10to15;
VID_d_tip_100_15toINF = VID_d_tip_100_neg20to20_15toINF;
VID_v_tip_100_0to5 = VID_v_tip_100_neg20to20_0to5;
VID_v_tip_100_5to10 = VID_v_tip_100_neg20to20_5to10;
VID_v_tip_100_10to15 = VID_v_tip_100_neg20to20_10to15;
VID_v_tip_100_15toINF = VID_v_tip_100_neg20to20_15toINF;
VID_angle_tip_100_0to5 = VID_angle_tip_100_neg20to20_0to5;
VID_angle_tip_100_5to10 = VID_angle_tip_100_neg20to20_5to10;
VID_angle_tip_100_10to15 = VID_angle_tip_100_neg20to20_10to15;
VID_angle_tip_100_15toINF = VID_angle_tip_100_neg20to20_15toINF;

subplot(4,7,4)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_0to5) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_0to5) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_0to5,1));
if ~isempty(train_data_logic_SS_0to5) && size(train_data_logic_SS_0to5,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['-20to20 deg | 0to5 mm | n = ' num2str(size(train_data_logic_SS_0to5,1)) ])

subplot(4,7,11)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_5to10) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_5to10) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_5to10,1));
if ~isempty(train_data_logic_SS_5to10) && size(train_data_logic_SS_5to10,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['-20to20 deg | 5to10 mm | n = ' num2str(size(train_data_logic_SS_5to10,1)) ])

subplot(4,7,18)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_10to15) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_10to15) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_10to15,1));
if ~isempty(train_data_logic_SS_10to15) && size(train_data_logic_SS_10to15,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['-20to20 deg | 10to15 mm | n = ' num2str(size(train_data_logic_SS_10to15,1)) ])

subplot(4,7,25)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_15toINF) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_15toINF) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_15toINF,1));
if ~isempty(train_data_logic_SS_15toINF) && size(train_data_logic_SS_15toINF,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['-20to20 deg | 15toINF mm | n = ' num2str(size(train_data_logic_SS_15toINF,1)) ])

%% neg20toneg60 licks
train_data_logic_SS_0to5 = train_data_logic_SS_neg20toneg60_0to5;
train_data_logic_SS_5to10 = train_data_logic_SS_neg20toneg60_5to10;
train_data_logic_SS_10to15 = train_data_logic_SS_neg20toneg60_10to15;
train_data_logic_SS_15toINF = train_data_logic_SS_neg20toneg60_15toINF;
train_data_logic_CS_0to5 = train_data_logic_CS_neg20toneg60_0to5;
train_data_logic_CS_5to10 = train_data_logic_CS_neg20toneg60_5to10;
train_data_logic_CS_10to15 = train_data_logic_CS_neg20toneg60_10to15;
train_data_logic_CS_15toINF = train_data_logic_CS_neg20toneg60_15toINF;
VID_d_tip_100_0to5 = VID_d_tip_100_neg20toneg60_0to5;
VID_d_tip_100_5to10 = VID_d_tip_100_neg20toneg60_5to10;
VID_d_tip_100_10to15 = VID_d_tip_100_neg20toneg60_10to15;
VID_d_tip_100_15toINF = VID_d_tip_100_neg20toneg60_15toINF;
VID_v_tip_100_0to5 = VID_v_tip_100_neg20toneg60_0to5;
VID_v_tip_100_5to10 = VID_v_tip_100_neg20toneg60_5to10;
VID_v_tip_100_10to15 = VID_v_tip_100_neg20toneg60_10to15;
VID_v_tip_100_15toINF = VID_v_tip_100_neg20toneg60_15toINF;
VID_angle_tip_100_0to5 = VID_angle_tip_100_neg20toneg60_0to5;
VID_angle_tip_100_5to10 = VID_angle_tip_100_neg20toneg60_5to10;
VID_angle_tip_100_10to15 = VID_angle_tip_100_neg20toneg60_10to15;
VID_angle_tip_100_15toINF = VID_angle_tip_100_neg20toneg60_15toINF;

subplot(4,7,3)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_0to5) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_0to5) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_0to5,1));
if ~isempty(train_data_logic_SS_0to5) && size(train_data_logic_SS_0to5,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['-20to-60 deg | 0to5 mm | n = ' num2str(size(train_data_logic_SS_0to5,1)) ])

subplot(4,7,10)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_5to10) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_5to10) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_5to10,1));
if ~isempty(train_data_logic_SS_5to10) && size(train_data_logic_SS_5to10,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['-20to-60 deg | 5to10 mm | n = ' num2str(size(train_data_logic_SS_5to10,1)) ])

subplot(4,7,17)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_10to15) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_10to15) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_10to15,1));
if ~isempty(train_data_logic_SS_10to15) && size(train_data_logic_SS_10to15,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['-20to-60 deg | 10to15 mm | n = ' num2str(size(train_data_logic_SS_10to15,1)) ])

subplot(4,7,24)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_15toINF) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_15toINF) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_15toINF,1));
if ~isempty(train_data_logic_SS_15toINF) && size(train_data_logic_SS_15toINF,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['-20to-60 deg | 15toINF mm | n = ' num2str(size(train_data_logic_SS_15toINF,1)) ])

%% neg60toneg80 licks
train_data_logic_SS_0to5 = train_data_logic_SS_neg60toneg80_0to5;
train_data_logic_SS_5to10 = train_data_logic_SS_neg60toneg80_5to10;
train_data_logic_SS_10to15 = train_data_logic_SS_neg60toneg80_10to15;
train_data_logic_SS_15toINF = train_data_logic_SS_neg60toneg80_15toINF;
train_data_logic_CS_0to5 = train_data_logic_CS_neg60toneg80_0to5;
train_data_logic_CS_5to10 = train_data_logic_CS_neg60toneg80_5to10;
train_data_logic_CS_10to15 = train_data_logic_CS_neg60toneg80_10to15;
train_data_logic_CS_15toINF = train_data_logic_CS_neg60toneg80_15toINF;
VID_d_tip_100_0to5 = VID_d_tip_100_neg60toneg80_0to5;
VID_d_tip_100_5to10 = VID_d_tip_100_neg60toneg80_5to10;
VID_d_tip_100_10to15 = VID_d_tip_100_neg60toneg80_10to15;
VID_d_tip_100_15toINF = VID_d_tip_100_neg60toneg80_15toINF;
VID_v_tip_100_0to5 = VID_v_tip_100_neg60toneg80_0to5;
VID_v_tip_100_5to10 = VID_v_tip_100_neg60toneg80_5to10;
VID_v_tip_100_10to15 = VID_v_tip_100_neg60toneg80_10to15;
VID_v_tip_100_15toINF = VID_v_tip_100_neg60toneg80_15toINF;
VID_angle_tip_100_0to5 = VID_angle_tip_100_neg60toneg80_0to5;
VID_angle_tip_100_5to10 = VID_angle_tip_100_neg60toneg80_5to10;
VID_angle_tip_100_10to15 = VID_angle_tip_100_neg60toneg80_10to15;
VID_angle_tip_100_15toINF = VID_angle_tip_100_neg60toneg80_15toINF;

subplot(4,7,2)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_0to5) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_0to5) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_0to5,1));
if ~isempty(train_data_logic_SS_0to5) && size(train_data_logic_SS_0to5,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['-60to-80 deg | 0to5 mm | n = ' num2str(size(train_data_logic_SS_0to5,1)) ])

subplot(4,7,9)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_5to10) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_5to10) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_5to10,1));
if ~isempty(train_data_logic_SS_5to10) && size(train_data_logic_SS_5to10,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['-60to-80 deg | 5to10 mm | n = ' num2str(size(train_data_logic_SS_5to10,1)) ])

subplot(4,7,16)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_10to15) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_10to15) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_10to15,1));
if ~isempty(train_data_logic_SS_10to15) && size(train_data_logic_SS_10to15,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['-60to-80 deg | 10to15 mm | n = ' num2str(size(train_data_logic_SS_10to15,1)) ])

subplot(4,7,23)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_15toINF) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_15toINF) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_15toINF,1));
if ~isempty(train_data_logic_SS_15toINF) && size(train_data_logic_SS_15toINF,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['-60to-80 deg | 15toINF mm | n = ' num2str(size(train_data_logic_SS_15toINF,1)) ])

%% neg80tonegINF licks
train_data_logic_SS_0to5 = train_data_logic_SS_neg80tonegINF_0to5;
train_data_logic_SS_5to10 = train_data_logic_SS_neg80tonegINF_5to10;
train_data_logic_SS_10to15 = train_data_logic_SS_neg80tonegINF_10to15;
train_data_logic_SS_15toINF = train_data_logic_SS_neg80tonegINF_15toINF;
train_data_logic_CS_0to5 = train_data_logic_CS_neg80tonegINF_0to5;
train_data_logic_CS_5to10 = train_data_logic_CS_neg80tonegINF_5to10;
train_data_logic_CS_10to15 = train_data_logic_CS_neg80tonegINF_10to15;
train_data_logic_CS_15toINF = train_data_logic_CS_neg80tonegINF_15toINF;
VID_d_tip_100_0to5 = VID_d_tip_100_neg80tonegINF_0to5;
VID_d_tip_100_5to10 = VID_d_tip_100_neg80tonegINF_5to10;
VID_d_tip_100_10to15 = VID_d_tip_100_neg80tonegINF_10to15;
VID_d_tip_100_15toINF = VID_d_tip_100_neg80tonegINF_15toINF;
VID_v_tip_100_0to5 = VID_v_tip_100_neg80tonegINF_0to5;
VID_v_tip_100_5to10 = VID_v_tip_100_neg80tonegINF_5to10;
VID_v_tip_100_10to15 = VID_v_tip_100_neg80tonegINF_10to15;
VID_v_tip_100_15toINF = VID_v_tip_100_neg80tonegINF_15toINF;
VID_angle_tip_100_0to5 = VID_angle_tip_100_neg80tonegINF_0to5;
VID_angle_tip_100_5to10 = VID_angle_tip_100_neg80tonegINF_5to10;
VID_angle_tip_100_10to15 = VID_angle_tip_100_neg80tonegINF_10to15;
VID_angle_tip_100_15toINF = VID_angle_tip_100_neg80tonegINF_15toINF;

subplot(4,7,1)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_0to5) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_0to5) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_0to5,1));
if ~isempty(train_data_logic_SS_0to5) && size(train_data_logic_SS_0to5,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['-80to-INF deg | 0to5 mm | n = ' num2str(size(train_data_logic_SS_0to5,1)) ])

subplot(4,7,8)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_5to10) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_5to10) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_5to10,1));
if ~isempty(train_data_logic_SS_5to10) && size(train_data_logic_SS_5to10,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['-80to-INF deg | 5to10 mm | n = ' num2str(size(train_data_logic_SS_5to10,1)) ])

subplot(4,7,15)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_10to15) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_10to15) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_10to15,1));
if ~isempty(train_data_logic_SS_10to15) && size(train_data_logic_SS_10to15,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['-80to-INF deg | 10to15 mm | n = ' num2str(size(train_data_logic_SS_10to15,1)) ])

subplot(4,7,22)
% SS
hold on
firing_SS_ = nanmean(train_data_logic_SS_15toINF) * 100;
firing_SS_sd = nanstd(train_data_logic_SS_15toINF) * 100;
firing_SS_se = firing_SS_sd/sqrt(size(train_data_logic_SS_15toINF,1));
if ~isempty(train_data_logic_SS_15toINF) && size(train_data_logic_SS_15toINF,1) > 1 
shade(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), 'k', inds_span,  ESN_smooth(firing_SS_ - firing_SS_se,5), 'k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(inds_span, ESN_smooth(firing_SS_,5), '-k', 'LineWidth', 2);
% plot(inds_span, ESN_smooth(firing_SS_ + firing_SS_se,5), '-k', 'LineWidth', 1);
% plot(inds_span, ESN_smooth(firing_SS_ - firing_SS_se,5), '-k', 'LineWidth', 1);
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', 'k')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
xline(0, 'k', 'LineWidth', 1.5);
title(['-80to-INF deg | 15toINF mm | n = ' num2str(size(train_data_logic_SS_15toINF,1)) ])

ESN_Beautify_Plot


end

%% function single_dataset_raster
function raster_data = single_dataset_raster(EPHYS, EPHYS_inds_event, VID, VID_inds_event, session_type)
% function raster_data = single_dataset_raster(EPHYS, VID, EPHYS_inds_event, VID_inds_event, onset_OR_vmax_OR_dmax_OR_vmin_OR_offset)

% % time
EPHYS_time_100     = EPHYS.CH_EVE.EPHYS_time_100;
EPHYS_time_30K     = EPHYS.CH_EVE.EPHYS_time_30K;


% % CS and SS data
SS_train_aligned     = EPHYS.CH_EVE.EPHYS_SS_train_100;
CS_train_aligned     = EPHYS.CH_EVE.EPHYS_CS_train_100;

% % Kinematics
VID_d_tip_100 = EPHYS.CH_EVE.VID_d_tip_100;
VID_v_tip_100 = EPHYS.CH_EVE.VID_v_tip_100;
VID_angle_tip_100 = EPHYS.CH_EVE.VID_angle_tip_100;

% % Eye events
EPHYS_primSac_onset_train_aligned  = false(size(SS_train_aligned));
EPHYS_corrSac_onset_train_aligned  = false(size(SS_train_aligned));
EPHYS_primSac_onset_train_aligned( EPHYS.CH_EVE.EPHYS_EB_aligned_ind_primSac_onset_100)  = true;
EPHYS_corrSac_onset_train_aligned( EPHYS.CH_EVE.EPHYS_EB_aligned_ind_corrSac_onset_100)  = true;

% % Lick events
EPHYS_lick_onset_train_aligned    = false(size(SS_train_aligned));
EPHYS_lick_vmax_train_aligned  = false(size(SS_train_aligned));
EPHYS_lick_dmax_train_aligned   = false(size(SS_train_aligned));
EPHYS_lick_vmin_train_aligned = false(size(SS_train_aligned));
EPHYS_lick_offset_train_aligned  = false(size(SS_train_aligned));
EPHYS_lick_onset_train_aligned( EPHYS.CH_EVE.EPHYS_LED_aligned_ind_lick_onset_100)      = true;
EPHYS_lick_vmax_train_aligned( EPHYS.CH_EVE.EPHYS_LED_aligned_ind_lick_vmax_100)        = true;
EPHYS_lick_dmax_train_aligned(  EPHYS.CH_EVE.EPHYS_LED_aligned_ind_lick_dmax_100)       = true;
EPHYS_lick_vmin_train_aligned(EPHYS.CH_EVE.EPHYS_LED_aligned_ind_lick_vmin_100)         = true;
EPHYS_lick_offset_train_aligned( EPHYS.CH_EVE.EPHYS_LED_aligned_ind_lick_offset_100)    = true;

num_licks = size(EPHYS_inds_event,1);

%inds bout
inds_str_bout = ismember(VID.DLC.IND.ind_lick_onset,VID.DLC.IND.ind_lick_onset_str_bout);
inds_str_bout = inds_str_bout(1:num_licks-1);
inds_end_bout = ismember(VID.DLC.IND.ind_lick_onset,VID.DLC.IND.ind_lick_onset_end_bout);
inds_end_bout = inds_end_bout(1:num_licks-1);

% correct for bout str & end offset
if (find(inds_str_bout == 1, 1, 'last')) > (find(inds_end_bout == 1, 1, 'last'))
    inds_str_bout(find(inds_str_bout == 1, 1, 'last')) = 0;
end

%inds class of bout
if (session_type  == 2 && sum(contains(fields(VID.DLC.CLASS), 'is_bout_l') > 0))
    is_bout_l = logical(VID.DLC.CLASS.is_bout_l(1:sum(inds_str_bout)));
end
if (sum(contains(fields(VID.DLC.CLASS), 'is_bout_r') > 0))
    is_bout_r = logical(VID.DLC.CLASS.is_bout_r(1:sum(inds_str_bout)));
end

% inds class of lick
is_grooming = logical(VID.DLC.CLASS.is_grooming_lick);
is_grooming = is_grooming(1:num_licks);
is_r = logical(VID.DLC.CLASS.is_r_reward_inner_tube_lick + VID.DLC.CLASS.is_r_reward_outer_tube_lick + VID.DLC.CLASS.is_r_noreward_inner_tube_lick + VID.DLC.CLASS.is_r_noreward_outer_tube_lick);
is_r = is_r(1:num_licks);
is_r_1 = logical(VID.DLC.CLASS.is_r_reward_inner_tube_lick + VID.DLC.CLASS.is_r_reward_outer_tube_lick);
is_r_1 = is_r_1(1:num_licks);
is_r_0 =  logical(VID.DLC.CLASS.is_r_noreward_inner_tube_lick + VID.DLC.CLASS.is_r_noreward_outer_tube_lick);
is_r_0 = is_r_0(1:num_licks);
if (session_type  == 2)
    is_l =  logical(VID.DLC.CLASS.is_l_reward_inner_tube_lick + VID.DLC.CLASS.is_l_reward_outer_tube_lick + VID.DLC.CLASS.is_l_noreward_inner_tube_lick + VID.DLC.CLASS.is_l_noreward_outer_tube_lick);
    is_l = is_l(1:num_licks);
    inds_l_1 = logical( VID.DLC.CLASS.is_l_reward_inner_tube_lick + VID.DLC.CLASS.is_l_reward_outer_tube_lick);
    inds_l_1 = inds_l_1(1:num_licks);
    inds_l_0 = logical(VID.DLC.CLASS.is_l_noreward_inner_tube_lick + VID.DLC.CLASS.is_l_noreward_outer_tube_lick);
    inds_l_0 = inds_l_0(1:num_licks);
end

% inds kinematics - angles
is_80toINF = VID.DLC.KINEMATIC.angle_lick_max > 80;
is_80toINF = is_80toINF(1:num_licks);
is_60to80 = VID.DLC.KINEMATIC.angle_lick_max > 60 & VID.DLC.KINEMATIC.angle_lick_max <= 80;
is_60to80 = is_60to80(1:num_licks);
is_20to60 = VID.DLC.KINEMATIC.angle_lick_max > 20 & VID.DLC.KINEMATIC.angle_lick_max <= 60;
is_20to60 = is_20to60(1:num_licks);
is_neg20to20 = VID.DLC.KINEMATIC.angle_lick_max >= -20 & VID.DLC.KINEMATIC.angle_lick_max <= 20;
is_neg20to20 = is_neg20to20(1:num_licks);
is_neg20toneg60 = VID.DLC.KINEMATIC.angle_lick_max >= -60 & VID.DLC.KINEMATIC.angle_lick_max < -20;
is_neg20toneg60 = is_neg20toneg60(1:num_licks);
is_neg60toneg80 = VID.DLC.KINEMATIC.angle_lick_max >= -80 & VID.DLC.KINEMATIC.angle_lick_max < -60;
is_neg60toneg80 = is_neg60toneg80(1:num_licks);
is_neg80tonegINF = VID.DLC.KINEMATIC.angle_lick_max < -80;
is_neg80tonegINF = is_neg80tonegINF(1:num_licks);

% inds kinematics - displacement
is_0to5 = VID.DLC.KINEMATIC.d_lick_max > 0 & VID.DLC.KINEMATIC.d_lick_max <=5;
is_0to5 = is_0to5(1:num_licks);
is_5to10 = VID.DLC.KINEMATIC.d_lick_max > 5 & VID.DLC.KINEMATIC.d_lick_max <=10;
is_5to10 = is_5to10(1:num_licks);
is_10to15 = VID.DLC.KINEMATIC.d_lick_max > 10 & VID.DLC.KINEMATIC.d_lick_max <=15;
is_10to15 = is_10to15(1:num_licks);
is_15toINF = VID.DLC.KINEMATIC.d_lick_max > 15;
is_15toINF = is_15toINF(1:num_licks);



% Expand EPHYS_inds_event to include corr: [0 1]s
% for counter_lick = 1 : length(EPHYS_inds_event)
%     EPHYS_inds_event_corr_UB(counter_lick, :) = EPHYS_inds_event(counter_lick,end)+1 : EPHYS_inds_event(counter_lick,end)+ 70;
% end
% EPHYS_inds_event_corr = [EPHYS_inds_event(:,31:60) EPHYS_inds_event_corr_UB];
% EPHYS_inds_event_corr(EPHYS_inds_event_corr > length(SS_train_aligned)) = 1;

% Expand EPHYS_inds_event to include corr: [-1 1]s
for counter_lick = 1 : length(EPHYS_inds_event)
    EPHYS_inds_event_corr_LB(counter_lick, :) = EPHYS_inds_event(counter_lick,1)- 70 : EPHYS_inds_event(counter_lick,1)-1;
    EPHYS_inds_event_corr_UB(counter_lick, :) = EPHYS_inds_event(counter_lick,end)+1 : EPHYS_inds_event(counter_lick,end)+ 70;
end
EPHYS_inds_event_corr = [EPHYS_inds_event_corr_LB EPHYS_inds_event EPHYS_inds_event_corr_UB];
EPHYS_inds_event_corr(EPHYS_inds_event_corr > length(SS_train_aligned)) = 1;

% Expand VID_inds_event to include corr: [-1 1]s
for counter_lick = 1 : length(VID_inds_event)
    VID_inds_event_corr_LB(counter_lick, :) = VID_inds_event(counter_lick,1)- 70 : VID_inds_event(counter_lick,1)-1;
    VID_inds_event_corr_UB(counter_lick, :) = VID_inds_event(counter_lick,end)+1 : VID_inds_event(counter_lick,end)+ 70;
end
VID_inds_event_corr = [VID_inds_event_corr_LB VID_inds_event VID_inds_event_corr_UB];
VID_inds_event_corr(VID_inds_event_corr > length(VID_d_tip_100)) = 1;
VID_inds_event_corr(VID_inds_event_corr > length(SS_train_aligned)) = 1;
VID_inds_event_corr(VID_inds_event_corr<1) = 1;

% 80toINF licks
train_data_logic_SS_80toINF_0to5 = SS_train_aligned(EPHYS_inds_event(is_80toINF & is_0to5,:));
train_data_logic_SS_80toINF_5to10 = SS_train_aligned(EPHYS_inds_event(is_80toINF & is_5to10,:));
train_data_logic_SS_80toINF_10to15 = SS_train_aligned(EPHYS_inds_event(is_80toINF & is_10to15,:));
train_data_logic_SS_80toINF_15toINF = SS_train_aligned(EPHYS_inds_event(is_80toINF & is_15toINF,:));
train_data_logic_CS_80toINF_0to5 = CS_train_aligned(EPHYS_inds_event(is_80toINF & is_0to5,:));
train_data_logic_CS_80toINF_5to10 = CS_train_aligned(EPHYS_inds_event(is_80toINF & is_5to10,:));
train_data_logic_CS_80toINF_10to15 = CS_train_aligned(EPHYS_inds_event(is_80toINF & is_10to15,:));
train_data_logic_CS_80toINF_15toINF = CS_train_aligned(EPHYS_inds_event(is_80toINF & is_15toINF,:));
VID_d_tip_100_80toINF_0to5 = VID_d_tip_100(VID_inds_event(is_80toINF & is_0to5,:));
VID_d_tip_100_80toINF_5to10 = VID_d_tip_100(VID_inds_event(is_80toINF & is_5to10,:));
VID_d_tip_100_80toINF_10to15 = VID_d_tip_100(VID_inds_event(is_80toINF & is_10to15,:));
VID_d_tip_100_80toINF_15toINF = VID_d_tip_100(VID_inds_event(is_80toINF & is_15toINF,:));
VID_v_tip_100_80toINF_0to5 = VID_v_tip_100(VID_inds_event(is_80toINF & is_0to5,:));
VID_v_tip_100_80toINF_5to10 = VID_v_tip_100(VID_inds_event(is_80toINF & is_5to10,:));
VID_v_tip_100_80toINF_10to15 = VID_v_tip_100(VID_inds_event(is_80toINF & is_10to15,:));
VID_v_tip_100_80toINF_15toINF = VID_v_tip_100(VID_inds_event(is_80toINF & is_15toINF,:));
VID_angle_tip_100_80toINF_0to5 = VID_angle_tip_100(VID_inds_event(is_80toINF & is_0to5,:));
VID_angle_tip_100_80toINF_5to10 = VID_angle_tip_100(VID_inds_event(is_80toINF & is_5to10,:));
VID_angle_tip_100_80toINF_10to15 = VID_angle_tip_100(VID_inds_event(is_80toINF & is_10to15,:));
VID_angle_tip_100_80toINF_15toINF = VID_angle_tip_100(VID_inds_event(is_80toINF & is_15toINF,:));

% 60to80 licks
train_data_logic_SS_60to80_0to5 = SS_train_aligned(EPHYS_inds_event(is_60to80 & is_0to5,:));
train_data_logic_SS_60to80_5to10 = SS_train_aligned(EPHYS_inds_event(is_60to80 & is_5to10,:));
train_data_logic_SS_60to80_10to15 = SS_train_aligned(EPHYS_inds_event(is_60to80 & is_10to15,:));
train_data_logic_SS_60to80_15toINF = SS_train_aligned(EPHYS_inds_event(is_60to80 & is_15toINF,:));
train_data_logic_CS_60to80_0to5 = CS_train_aligned(EPHYS_inds_event(is_60to80 & is_0to5,:));
train_data_logic_CS_60to80_5to10 = CS_train_aligned(EPHYS_inds_event(is_60to80 & is_5to10,:));
train_data_logic_CS_60to80_10to15 = CS_train_aligned(EPHYS_inds_event(is_60to80 & is_10to15,:));
train_data_logic_CS_60to80_15toINF = CS_train_aligned(EPHYS_inds_event(is_60to80 & is_15toINF,:));
VID_d_tip_100_60to80_0to5 = VID_d_tip_100(VID_inds_event(is_60to80 & is_0to5,:));
VID_d_tip_100_60to80_5to10 = VID_d_tip_100(VID_inds_event(is_60to80 & is_5to10,:));
VID_d_tip_100_60to80_10to15 = VID_d_tip_100(VID_inds_event(is_60to80 & is_10to15,:));
VID_d_tip_100_60to80_15toINF = VID_d_tip_100(VID_inds_event(is_60to80 & is_15toINF,:));
VID_v_tip_100_60to80_0to5 = VID_v_tip_100(VID_inds_event(is_60to80 & is_0to5,:));
VID_v_tip_100_60to80_5to10 = VID_v_tip_100(VID_inds_event(is_60to80 & is_5to10,:));
VID_v_tip_100_60to80_10to15 = VID_v_tip_100(VID_inds_event(is_60to80 & is_10to15,:));
VID_v_tip_100_60to80_15toINF = VID_v_tip_100(VID_inds_event(is_60to80 & is_15toINF,:));
VID_angle_tip_100_60to80_0to5 = VID_angle_tip_100(VID_inds_event(is_60to80 & is_0to5,:));
VID_angle_tip_100_60to80_5to10 = VID_angle_tip_100(VID_inds_event(is_60to80 & is_5to10,:));
VID_angle_tip_100_60to80_10to15 = VID_angle_tip_100(VID_inds_event(is_60to80 & is_10to15,:));
VID_angle_tip_100_60to80_15toINF = VID_angle_tip_100(VID_inds_event(is_60to80 & is_15toINF,:));

% 20to60 licks
train_data_logic_SS_20to60_0to5 = SS_train_aligned(EPHYS_inds_event(is_20to60 & is_0to5,:));
train_data_logic_SS_20to60_5to10 = SS_train_aligned(EPHYS_inds_event(is_20to60 & is_5to10,:));
train_data_logic_SS_20to60_10to15 = SS_train_aligned(EPHYS_inds_event(is_20to60 & is_10to15,:));
train_data_logic_SS_20to60_15toINF = SS_train_aligned(EPHYS_inds_event(is_20to60 & is_15toINF,:));
train_data_logic_CS_20to60_0to5 = CS_train_aligned(EPHYS_inds_event(is_20to60 & is_0to5,:));
train_data_logic_CS_20to60_5to10 = CS_train_aligned(EPHYS_inds_event(is_20to60 & is_5to10,:));
train_data_logic_CS_20to60_10to15 = CS_train_aligned(EPHYS_inds_event(is_20to60 & is_10to15,:));
train_data_logic_CS_20to60_15toINF = CS_train_aligned(EPHYS_inds_event(is_20to60 & is_15toINF,:));
VID_d_tip_100_20to60_0to5 = VID_d_tip_100(VID_inds_event(is_20to60 & is_0to5,:));
VID_d_tip_100_20to60_5to10 = VID_d_tip_100(VID_inds_event(is_20to60 & is_5to10,:));
VID_d_tip_100_20to60_10to15 = VID_d_tip_100(VID_inds_event(is_20to60 & is_10to15,:));
VID_d_tip_100_20to60_15toINF = VID_d_tip_100(VID_inds_event(is_20to60 & is_15toINF,:));
VID_v_tip_100_20to60_0to5 = VID_v_tip_100(VID_inds_event(is_20to60 & is_0to5,:));
VID_v_tip_100_20to60_5to10 = VID_v_tip_100(VID_inds_event(is_20to60 & is_5to10,:));
VID_v_tip_100_20to60_10to15 = VID_v_tip_100(VID_inds_event(is_20to60 & is_10to15,:));
VID_v_tip_100_20to60_15toINF = VID_v_tip_100(VID_inds_event(is_20to60 & is_15toINF,:));
VID_angle_tip_100_20to60_0to5 = VID_angle_tip_100(VID_inds_event(is_20to60 & is_0to5,:));
VID_angle_tip_100_20to60_5to10 = VID_angle_tip_100(VID_inds_event(is_20to60 & is_5to10,:));
VID_angle_tip_100_20to60_10to15 = VID_angle_tip_100(VID_inds_event(is_20to60 & is_10to15,:));
VID_angle_tip_100_20to60_15toINF = VID_angle_tip_100(VID_inds_event(is_20to60 & is_15toINF,:));

% neg20to20 licks
train_data_logic_SS_neg20to20_0to5 = SS_train_aligned(EPHYS_inds_event(is_neg20to20 & is_0to5,:));
train_data_logic_SS_neg20to20_5to10 = SS_train_aligned(EPHYS_inds_event(is_neg20to20 & is_5to10,:));
train_data_logic_SS_neg20to20_10to15 = SS_train_aligned(EPHYS_inds_event(is_neg20to20 & is_10to15,:));
train_data_logic_SS_neg20to20_15toINF = SS_train_aligned(EPHYS_inds_event(is_neg20to20 & is_15toINF,:));
train_data_logic_CS_neg20to20_0to5 = CS_train_aligned(EPHYS_inds_event(is_neg20to20 & is_0to5,:));
train_data_logic_CS_neg20to20_5to10 = CS_train_aligned(EPHYS_inds_event(is_neg20to20 & is_5to10,:));
train_data_logic_CS_neg20to20_10to15 = CS_train_aligned(EPHYS_inds_event(is_neg20to20 & is_10to15,:));
train_data_logic_CS_neg20to20_15toINF = CS_train_aligned(EPHYS_inds_event(is_neg20to20 & is_15toINF,:));
VID_d_tip_100_neg20to20_0to5 = VID_d_tip_100(VID_inds_event(is_neg20to20 & is_0to5,:));
VID_d_tip_100_neg20to20_5to10 = VID_d_tip_100(VID_inds_event(is_neg20to20 & is_5to10,:));
VID_d_tip_100_neg20to20_10to15 = VID_d_tip_100(VID_inds_event(is_neg20to20 & is_10to15,:));
VID_d_tip_100_neg20to20_15toINF = VID_d_tip_100(VID_inds_event(is_neg20to20 & is_15toINF,:));
VID_v_tip_100_neg20to20_0to5 = VID_v_tip_100(VID_inds_event(is_neg20to20 & is_0to5,:));
VID_v_tip_100_neg20to20_5to10 = VID_v_tip_100(VID_inds_event(is_neg20to20 & is_5to10,:));
VID_v_tip_100_neg20to20_10to15 = VID_v_tip_100(VID_inds_event(is_neg20to20 & is_10to15,:));
VID_v_tip_100_neg20to20_15toINF = VID_v_tip_100(VID_inds_event(is_neg20to20 & is_15toINF,:));
VID_angle_tip_100_neg20to20_0to5 = VID_angle_tip_100(VID_inds_event(is_neg20to20 & is_0to5,:));
VID_angle_tip_100_neg20to20_5to10 = VID_angle_tip_100(VID_inds_event(is_neg20to20 & is_5to10,:));
VID_angle_tip_100_neg20to20_10to15 = VID_angle_tip_100(VID_inds_event(is_neg20to20 & is_10to15,:));
VID_angle_tip_100_neg20to20_15toINF = VID_angle_tip_100(VID_inds_event(is_neg20to20 & is_15toINF,:));

% neg20toneg60 licks
train_data_logic_SS_neg20toneg60_0to5 = SS_train_aligned(EPHYS_inds_event(is_neg20toneg60 & is_0to5,:));
train_data_logic_SS_neg20toneg60_5to10 = SS_train_aligned(EPHYS_inds_event(is_neg20toneg60 & is_5to10,:));
train_data_logic_SS_neg20toneg60_10to15 = SS_train_aligned(EPHYS_inds_event(is_neg20toneg60 & is_10to15,:));
train_data_logic_SS_neg20toneg60_15toINF = SS_train_aligned(EPHYS_inds_event(is_neg20toneg60 & is_15toINF,:));
train_data_logic_CS_neg20toneg60_0to5 = CS_train_aligned(EPHYS_inds_event(is_neg20toneg60 & is_0to5,:));
train_data_logic_CS_neg20toneg60_5to10 = CS_train_aligned(EPHYS_inds_event(is_neg20toneg60 & is_5to10,:));
train_data_logic_CS_neg20toneg60_10to15 = CS_train_aligned(EPHYS_inds_event(is_neg20toneg60 & is_10to15,:));
train_data_logic_CS_neg20toneg60_15toINF = CS_train_aligned(EPHYS_inds_event(is_neg20toneg60 & is_15toINF,:));
VID_d_tip_100_neg20toneg60_0to5 = VID_d_tip_100(VID_inds_event(is_neg20toneg60 & is_0to5,:));
VID_d_tip_100_neg20toneg60_5to10 = VID_d_tip_100(VID_inds_event(is_neg20toneg60 & is_5to10,:));
VID_d_tip_100_neg20toneg60_10to15 = VID_d_tip_100(VID_inds_event(is_neg20toneg60 & is_10to15,:));
VID_d_tip_100_neg20toneg60_15toINF = VID_d_tip_100(VID_inds_event(is_neg20toneg60 & is_15toINF,:));
VID_v_tip_100_neg20toneg60_0to5 = VID_v_tip_100(VID_inds_event(is_neg20toneg60 & is_0to5,:));
VID_v_tip_100_neg20toneg60_5to10 = VID_v_tip_100(VID_inds_event(is_neg20toneg60 & is_5to10,:));
VID_v_tip_100_neg20toneg60_10to15 = VID_v_tip_100(VID_inds_event(is_neg20toneg60 & is_10to15,:));
VID_v_tip_100_neg20toneg60_15toINF = VID_v_tip_100(VID_inds_event(is_neg20toneg60 & is_15toINF,:));
VID_angle_tip_100_neg20toneg60_0to5 = VID_angle_tip_100(VID_inds_event(is_neg20toneg60 & is_0to5,:));
VID_angle_tip_100_neg20toneg60_5to10 = VID_angle_tip_100(VID_inds_event(is_neg20toneg60 & is_5to10,:));
VID_angle_tip_100_neg20toneg60_10to15 = VID_angle_tip_100(VID_inds_event(is_neg20toneg60 & is_10to15,:));
VID_angle_tip_100_neg20toneg60_15toINF = VID_angle_tip_100(VID_inds_event(is_neg20toneg60 & is_15toINF,:));

% neg60toneg80 licks
train_data_logic_SS_neg60toneg80_0to5 = SS_train_aligned(EPHYS_inds_event(is_neg60toneg80 & is_0to5,:));
train_data_logic_SS_neg60toneg80_5to10 = SS_train_aligned(EPHYS_inds_event(is_neg60toneg80 & is_5to10,:));
train_data_logic_SS_neg60toneg80_10to15 = SS_train_aligned(EPHYS_inds_event(is_neg60toneg80 & is_10to15,:));
train_data_logic_SS_neg60toneg80_15toINF = SS_train_aligned(EPHYS_inds_event(is_neg60toneg80 & is_15toINF,:));
train_data_logic_CS_neg60toneg80_0to5 = CS_train_aligned(EPHYS_inds_event(is_neg60toneg80 & is_0to5,:));
train_data_logic_CS_neg60toneg80_5to10 = CS_train_aligned(EPHYS_inds_event(is_neg60toneg80 & is_5to10,:));
train_data_logic_CS_neg60toneg80_10to15 = CS_train_aligned(EPHYS_inds_event(is_neg60toneg80 & is_10to15,:));
train_data_logic_CS_neg60toneg80_15toINF = CS_train_aligned(EPHYS_inds_event(is_neg60toneg80 & is_15toINF,:));
VID_d_tip_100_neg60toneg80_0to5 = VID_d_tip_100(VID_inds_event(is_neg60toneg80 & is_0to5,:));
VID_d_tip_100_neg60toneg80_5to10 = VID_d_tip_100(VID_inds_event(is_neg60toneg80 & is_5to10,:));
VID_d_tip_100_neg60toneg80_10to15 = VID_d_tip_100(VID_inds_event(is_neg60toneg80 & is_10to15,:));
VID_d_tip_100_neg60toneg80_15toINF = VID_d_tip_100(VID_inds_event(is_neg60toneg80 & is_15toINF,:));
VID_v_tip_100_neg60toneg80_0to5 = VID_v_tip_100(VID_inds_event(is_neg60toneg80 & is_0to5,:));
VID_v_tip_100_neg60toneg80_5to10 = VID_v_tip_100(VID_inds_event(is_neg60toneg80 & is_5to10,:));
VID_v_tip_100_neg60toneg80_10to15 = VID_v_tip_100(VID_inds_event(is_neg60toneg80 & is_10to15,:));
VID_v_tip_100_neg60toneg80_15toINF = VID_v_tip_100(VID_inds_event(is_neg60toneg80 & is_15toINF,:));
VID_angle_tip_100_neg60toneg80_0to5 = VID_angle_tip_100(VID_inds_event(is_neg60toneg80 & is_0to5,:));
VID_angle_tip_100_neg60toneg80_5to10 = VID_angle_tip_100(VID_inds_event(is_neg60toneg80 & is_5to10,:));
VID_angle_tip_100_neg60toneg80_10to15 = VID_angle_tip_100(VID_inds_event(is_neg60toneg80 & is_10to15,:));
VID_angle_tip_100_neg60toneg80_15toINF = VID_angle_tip_100(VID_inds_event(is_neg20toneg60 & is_15toINF,:));

% neg80tonegINF licks
train_data_logic_SS_neg80tonegINF_0to5 = SS_train_aligned(EPHYS_inds_event(is_neg80tonegINF & is_0to5,:));
train_data_logic_SS_neg80tonegINF_5to10 = SS_train_aligned(EPHYS_inds_event(is_neg80tonegINF & is_5to10,:));
train_data_logic_SS_neg80tonegINF_10to15 = SS_train_aligned(EPHYS_inds_event(is_neg80tonegINF & is_10to15,:));
train_data_logic_SS_neg80tonegINF_15toINF = SS_train_aligned(EPHYS_inds_event(is_neg80tonegINF & is_15toINF,:));
train_data_logic_CS_neg80tonegINF_0to5 = CS_train_aligned(EPHYS_inds_event(is_neg80tonegINF & is_0to5,:));
train_data_logic_CS_neg80tonegINF_5to10 = CS_train_aligned(EPHYS_inds_event(is_neg80tonegINF & is_5to10,:));
train_data_logic_CS_neg80tonegINF_10to15 = CS_train_aligned(EPHYS_inds_event(is_neg80tonegINF & is_10to15,:));
train_data_logic_CS_neg80tonegINF_15toINF = CS_train_aligned(EPHYS_inds_event(is_neg80tonegINF & is_15toINF,:));
VID_d_tip_100_neg80tonegINF_0to5 = VID_d_tip_100(VID_inds_event(is_neg80tonegINF & is_0to5,:));
VID_d_tip_100_neg80tonegINF_5to10 = VID_d_tip_100(VID_inds_event(is_neg80tonegINF & is_5to10,:));
VID_d_tip_100_neg80tonegINF_10to15 = VID_d_tip_100(VID_inds_event(is_neg80tonegINF & is_10to15,:));
VID_d_tip_100_neg80tonegINF_15toINF = VID_d_tip_100(VID_inds_event(is_neg80tonegINF & is_15toINF,:));
VID_v_tip_100_neg80tonegINF_0to5 = VID_v_tip_100(VID_inds_event(is_neg80tonegINF & is_0to5,:));
VID_v_tip_100_neg80tonegINF_5to10 = VID_v_tip_100(VID_inds_event(is_neg80tonegINF & is_5to10,:));
VID_v_tip_100_neg80tonegINF_10to15 = VID_v_tip_100(VID_inds_event(is_neg80tonegINF & is_10to15,:));
VID_v_tip_100_neg80tonegINF_15toINF = VID_v_tip_100(VID_inds_event(is_neg80tonegINF & is_15toINF,:));
VID_angle_tip_100_neg80tonegINF_0to5 = VID_angle_tip_100(VID_inds_event(is_neg80tonegINF & is_0to5,:));
VID_angle_tip_100_neg80tonegINF_5to10 = VID_angle_tip_100(VID_inds_event(is_neg80tonegINF & is_5to10,:));
VID_angle_tip_100_neg80tonegINF_10to15 = VID_angle_tip_100(VID_inds_event(is_neg80tonegINF & is_10to15,:));
VID_angle_tip_100_neg80tonegINF_15toINF = VID_angle_tip_100(VID_inds_event(is_neg80tonegINF & is_15toINF,:));

% All licks
train_data_logic_SS_all_corr = SS_train_aligned(EPHYS_inds_event_corr);
train_data_logic_CS_all_corr = CS_train_aligned(EPHYS_inds_event_corr);
VID_d_tip_100_all_corr = VID_d_tip_100(VID_inds_event_corr);
VID_v_tip_100_all_corr = VID_v_tip_100(VID_inds_event_corr);
VID_angle_tip_100_all_corr = VID_angle_tip_100(VID_inds_event_corr);
train_data_logic_lick_onset_all_corr = EPHYS_lick_onset_train_aligned(EPHYS_inds_event_corr);
train_data_logic_lick_vmax_all_corr = EPHYS_lick_vmax_train_aligned(EPHYS_inds_event_corr);
train_data_logic_lick_dmax_all_corr = EPHYS_lick_dmax_train_aligned(EPHYS_inds_event_corr);
train_data_logic_lick_vmin_all_corr = EPHYS_lick_vmin_train_aligned(EPHYS_inds_event_corr);
train_data_logic_lick_offset_all_corr = EPHYS_lick_offset_train_aligned(EPHYS_inds_event_corr);
train_data_logic_SS_all = SS_train_aligned(EPHYS_inds_event);
train_data_logic_CS_all = CS_train_aligned(EPHYS_inds_event);
train_data_logic_lick_onset_all = EPHYS_lick_onset_train_aligned(EPHYS_inds_event);
train_data_logic_lick_vmax_all = EPHYS_lick_vmax_train_aligned(EPHYS_inds_event);
train_data_logic_lick_dmax_all = EPHYS_lick_dmax_train_aligned(EPHYS_inds_event);
train_data_logic_lick_vmin_all = EPHYS_lick_vmin_train_aligned(EPHYS_inds_event);
train_data_logic_lick_offset_all = EPHYS_lick_offset_train_aligned(EPHYS_inds_event);
VID_d_tip_100_all = VID_d_tip_100(VID_inds_event);
VID_v_tip_100_all = VID_v_tip_100(VID_inds_event);
VID_angle_tip_100_all = VID_angle_tip_100(VID_inds_event);
train_data_logic_primSac_onset_all = EPHYS_primSac_onset_train_aligned(EPHYS_inds_event);
train_data_logic_corrSac_onset_all = EPHYS_corrSac_onset_train_aligned(EPHYS_inds_event);

% Grooming licks
train_data_logic_SS_grooming_corr = SS_train_aligned(EPHYS_inds_event_corr(is_grooming,:));
train_data_logic_CS_grooming_corr = CS_train_aligned(EPHYS_inds_event_corr(is_grooming,:));
VID_d_tip_100_grooming_corr = VID_d_tip_100(VID_inds_event_corr(is_grooming,:));
VID_v_tip_100_grooming_corr = VID_v_tip_100(VID_inds_event_corr(is_grooming,:));
VID_angle_tip_100_grooming_corr = VID_angle_tip_100(VID_inds_event_corr(is_grooming,:));
train_data_logic_lick_onset_grooming_corr = EPHYS_lick_onset_train_aligned(EPHYS_inds_event_corr(is_grooming,:));
train_data_logic_lick_vmax_grooming_corr = EPHYS_lick_vmax_train_aligned(EPHYS_inds_event_corr(is_grooming,:));
train_data_logic_lick_dmax_grooming_corr = EPHYS_lick_dmax_train_aligned(EPHYS_inds_event_corr(is_grooming,:));
train_data_logic_lick_vmin_grooming_corr = EPHYS_lick_vmin_train_aligned(EPHYS_inds_event_corr(is_grooming,:));
train_data_logic_lick_offset_grooming_corr = EPHYS_lick_offset_train_aligned(EPHYS_inds_event_corr(is_grooming,:));
train_data_logic_SS_grooming = SS_train_aligned(EPHYS_inds_event(is_grooming,:));
train_data_logic_CS_grooming = CS_train_aligned(EPHYS_inds_event(is_grooming,:));
train_data_logic_lick_onset_grooming = EPHYS_lick_onset_train_aligned(EPHYS_inds_event(is_grooming,:));
train_data_logic_lick_vmax_grooming = EPHYS_lick_vmax_train_aligned(EPHYS_inds_event(is_grooming,:));
train_data_logic_lick_dmax_grooming = EPHYS_lick_dmax_train_aligned(EPHYS_inds_event(is_grooming,:));
train_data_logic_lick_vmin_grooming = EPHYS_lick_vmin_train_aligned(EPHYS_inds_event(is_grooming,:));
train_data_logic_lick_offset_grooming = EPHYS_lick_offset_train_aligned(EPHYS_inds_event(is_grooming,:));
VID_d_tip_100_grooming = VID_d_tip_100(VID_inds_event(is_grooming,:));
VID_v_tip_100_grooming = VID_v_tip_100(VID_inds_event(is_grooming,:));
VID_angle_tip_100_grooming = VID_angle_tip_100(VID_inds_event(is_grooming,:));
train_data_logic_primSac_onset_grooming = EPHYS_primSac_onset_train_aligned(EPHYS_inds_event(is_grooming,:));
train_data_logic_corrSac_onset_grooming = EPHYS_corrSac_onset_train_aligned(EPHYS_inds_event(is_grooming,:));

% R licks
train_data_logic_SS_r_corr = SS_train_aligned(EPHYS_inds_event_corr(is_r,:));
train_data_logic_CS_r_corr = CS_train_aligned(EPHYS_inds_event_corr(is_r,:));
VID_d_tip_100_r_corr = VID_d_tip_100(VID_inds_event_corr(is_r,:));
VID_v_tip_100_r_corr = VID_v_tip_100(VID_inds_event_corr(is_r,:));
VID_angle_tip_100_r_corr = VID_angle_tip_100(VID_inds_event_corr(is_r,:));
train_data_logic_lick_onset_r_corr = EPHYS_lick_onset_train_aligned(EPHYS_inds_event_corr(is_r,:));
train_data_logic_lick_vmax_r_corr = EPHYS_lick_vmax_train_aligned(EPHYS_inds_event_corr(is_r,:));
train_data_logic_lick_dmax_r_corr = EPHYS_lick_dmax_train_aligned(EPHYS_inds_event_corr(is_r,:));
train_data_logic_lick_vmin_r_corr = EPHYS_lick_vmin_train_aligned(EPHYS_inds_event_corr(is_r,:));
train_data_logic_lick_offset_r_corr = EPHYS_lick_offset_train_aligned(EPHYS_inds_event_corr(is_r,:));
train_data_logic_SS_r = SS_train_aligned(EPHYS_inds_event(is_r,:));
train_data_logic_CS_r = CS_train_aligned(EPHYS_inds_event(is_r,:));
train_data_logic_CS_r_1 = CS_train_aligned(EPHYS_inds_event(is_r_1,:));
train_data_logic_CS_r_0 = CS_train_aligned(EPHYS_inds_event(is_r_0,:));
train_data_logic_lick_onset_r = EPHYS_lick_onset_train_aligned(EPHYS_inds_event(is_r,:));
train_data_logic_lick_vmax_r = EPHYS_lick_vmax_train_aligned(EPHYS_inds_event(is_r,:));
train_data_logic_lick_dmax_r = EPHYS_lick_dmax_train_aligned(EPHYS_inds_event(is_r,:));
train_data_logic_lick_vmin_r = EPHYS_lick_vmin_train_aligned(EPHYS_inds_event(is_r,:));
train_data_logic_lick_offset_r = EPHYS_lick_offset_train_aligned(EPHYS_inds_event(is_r,:));
VID_d_tip_100_r = VID_d_tip_100(VID_inds_event(is_r,:));
VID_v_tip_100_r = VID_v_tip_100(VID_inds_event(is_r,:));
VID_angle_tip_100_r = VID_angle_tip_100(VID_inds_event(is_r,:));
train_data_logic_primSac_onset_r = EPHYS_primSac_onset_train_aligned(EPHYS_inds_event(is_r,:));
train_data_logic_corrSac_onset_r = EPHYS_corrSac_onset_train_aligned(EPHYS_inds_event(is_r,:));

% L licks
if (VID.DLC.FILE.session_type  == 2)
    train_data_logic_SS_l_corr = SS_train_aligned(EPHYS_inds_event_corr(is_l,:));
    train_data_logic_CS_l_corr = CS_train_aligned(EPHYS_inds_event_corr(is_l,:));
    VID_d_tip_100_l_corr = VID_d_tip_100(VID_inds_event_corr(is_l,:));
    VID_v_tip_100_l_corr = VID_v_tip_100(VID_inds_event_corr(is_l,:));
    VID_angle_tip_100_l_corr = VID_angle_tip_100(VID_inds_event_corr(is_l,:))*-1;
    train_data_logic_lick_onset_l_corr = EPHYS_lick_onset_train_aligned(EPHYS_inds_event_corr(is_l,:));
    train_data_logic_lick_vmax_l_corr = EPHYS_lick_vmax_train_aligned(EPHYS_inds_event_corr(is_l,:));
    train_data_logic_lick_dmax_l_corr = EPHYS_lick_dmax_train_aligned(EPHYS_inds_event_corr(is_l,:));
    train_data_logic_lick_vmin_l_corr = EPHYS_lick_vmin_train_aligned(EPHYS_inds_event_corr(is_l,:));
    train_data_logic_lick_offset_l_corr = EPHYS_lick_offset_train_aligned(EPHYS_inds_event_corr(is_l,:));
    train_data_logic_SS_l = SS_train_aligned(EPHYS_inds_event(is_l,:));
    train_data_logic_CS_l = CS_train_aligned(EPHYS_inds_event(is_l,:));
    train_data_logic_CS_l_1 = CS_train_aligned(EPHYS_inds_event(inds_l_1,:));
    train_data_logic_CS_l_0 = CS_train_aligned(EPHYS_inds_event(inds_l_0,:));
    if  nansum(inds_l_0) ==1
        train_data_logic_CS_l_0 = [];
    end
    train_data_logic_lick_onset_l = EPHYS_lick_onset_train_aligned(EPHYS_inds_event(is_l,:));
    train_data_logic_lick_vmax_l = EPHYS_lick_vmax_train_aligned(EPHYS_inds_event(is_l,:));
    train_data_logic_lick_dmax_l = EPHYS_lick_dmax_train_aligned(EPHYS_inds_event(is_l,:));
    train_data_logic_lick_vmin_l = EPHYS_lick_vmin_train_aligned(EPHYS_inds_event(is_l,:));
    train_data_logic_lick_offset_l = EPHYS_lick_offset_train_aligned(EPHYS_inds_event(is_l,:));
    VID_d_tip_100_l = VID_d_tip_100(VID_inds_event(is_l,:));
    VID_v_tip_100_l = VID_v_tip_100(VID_inds_event(is_l,:));
    VID_angle_tip_100_l = VID_angle_tip_100(VID_inds_event(is_l,:))*-1;
    train_data_logic_primSac_onset_l = EPHYS_primSac_onset_train_aligned(EPHYS_inds_event(is_l,:));
    train_data_logic_corrSac_onset_l = EPHYS_corrSac_onset_train_aligned(EPHYS_inds_event(is_l,:));
end

% Expand EPHYS_inds_event to include bout
EPHYS_inds_event_bout_LB = nan(length(EPHYS_inds_event), 70);
EPHYS_inds_event_bout_UB = nan(length(EPHYS_inds_event), 70);
for counter_lick = 1 : length(EPHYS_inds_event)
    EPHYS_inds_event_bout_LB(counter_lick, :) = EPHYS_inds_event(counter_lick,1)- 70 : EPHYS_inds_event(counter_lick,1)-1;
    EPHYS_inds_event_bout_UB(counter_lick, :) = EPHYS_inds_event(counter_lick,end)+1 : EPHYS_inds_event(counter_lick,end)+ 70;
end
EPHYS_inds_event_bout = [EPHYS_inds_event_bout_LB EPHYS_inds_event EPHYS_inds_event_bout_UB];
EPHYS_inds_event_bout(EPHYS_inds_event_bout > length(SS_train_aligned)) = 1;

% Expand VID_inds_event to include bout
VID_inds_event_bout_LB = nan(length(VID_inds_event), 70);
VID_inds_event_bout_UB = nan(length(VID_inds_event), 70);
for counter_lick = 1 : length(VID_inds_event)
    VID_inds_event_bout_LB(counter_lick, :) = VID_inds_event(counter_lick,1)- 70 : VID_inds_event(counter_lick,1)-1;
    VID_inds_event_bout_UB(counter_lick, :) = VID_inds_event(counter_lick,end)+1 : VID_inds_event(counter_lick,end)+ 70;
end
VID_inds_event_bout = [VID_inds_event_bout_LB VID_inds_event VID_inds_event_bout_UB];
VID_inds_event_bout (VID_inds_event_bout > length(SS_train_aligned)) = 1;
VID_inds_event_bout(VID_inds_event_bout<1) = 1;

% Str bout licks
train_data_logic_SS_str_bout = SS_train_aligned(EPHYS_inds_event_bout(inds_str_bout,:));
train_data_logic_CS_str_bout = CS_train_aligned(EPHYS_inds_event_bout(inds_str_bout,:));
train_data_logic_lick_onset_str_bout = EPHYS_lick_onset_train_aligned(EPHYS_inds_event_bout(inds_str_bout,:));
train_data_logic_lick_vmax_str_bout = EPHYS_lick_vmax_train_aligned(EPHYS_inds_event_bout(inds_str_bout,:));
train_data_logic_lick_dmax_str_bout = EPHYS_lick_dmax_train_aligned(EPHYS_inds_event_bout(inds_str_bout,:));
train_data_logic_lick_vmin_str_bout = EPHYS_lick_vmin_train_aligned(EPHYS_inds_event_bout(inds_str_bout,:));
train_data_logic_lick_offset_str_bout = EPHYS_lick_offset_train_aligned(EPHYS_inds_event_bout(inds_str_bout,:));
train_data_logic_primSac_onset_str_bout = EPHYS_primSac_onset_train_aligned(EPHYS_inds_event_bout(inds_str_bout,:));
train_data_logic_corrSac_onset_str_bout = EPHYS_corrSac_onset_train_aligned(EPHYS_inds_event_bout(inds_str_bout,:));
VID_d_tip_100_str_bout = VID_d_tip_100(VID_inds_event_bout(inds_str_bout,:));
VID_v_tip_100_str_bout = VID_v_tip_100(VID_inds_event_bout(inds_str_bout,:));
VID_angle_tip_100_str_bout = VID_angle_tip_100(VID_inds_event_bout(inds_str_bout,:));
% Left bout
if (VID.DLC.FILE.session_type  == 2)
    train_data_logic_SS_str_bout_l = train_data_logic_SS_str_bout(is_bout_l,:);
    train_data_logic_CS_str_bout_l = train_data_logic_CS_str_bout(is_bout_l,:);
    train_data_logic_lick_onset_str_bout_l = train_data_logic_lick_onset_str_bout(is_bout_l,:);
    train_data_logic_lick_vmax_str_bout_l = train_data_logic_lick_vmax_str_bout(is_bout_l,:);
    train_data_logic_lick_dmax_str_bout_l = train_data_logic_lick_dmax_str_bout(is_bout_l,:);
    train_data_logic_lick_vmin_str_bout_l = train_data_logic_lick_vmin_str_bout(is_bout_l,:);
    train_data_logic_lick_offset_str_bout_l = train_data_logic_lick_offset_str_bout(is_bout_l,:);
    train_data_logic_primSac_onset_str_bout_l = train_data_logic_primSac_onset_str_bout(is_bout_l,:);
    train_data_logic_corrSac_onset_str_bout_l = train_data_logic_corrSac_onset_str_bout(is_bout_l,:);
    VID_d_tip_100_str_bout_l = VID_d_tip_100_str_bout(is_bout_l,:);
    VID_v_tip_100_str_bout_l = VID_v_tip_100_str_bout(is_bout_l,:);
    VID_angle_tip_100_str_bout_l = VID_angle_tip_100_str_bout(is_bout_l,:);
    
end
% Right bout
train_data_logic_SS_str_bout_r = train_data_logic_SS_str_bout(is_bout_r,:);
train_data_logic_CS_str_bout_r = train_data_logic_CS_str_bout(is_bout_r,:);
train_data_logic_lick_onset_str_bout_r = train_data_logic_lick_onset_str_bout(is_bout_r,:);
train_data_logic_lick_vmax_str_bout_r = train_data_logic_lick_vmax_str_bout(is_bout_r,:);
train_data_logic_lick_dmax_str_bout_r = train_data_logic_lick_dmax_str_bout(is_bout_r,:);
train_data_logic_lick_vmin_str_bout_r = train_data_logic_lick_vmin_str_bout(is_bout_r,:);
train_data_logic_lick_offset_str_bout_r = train_data_logic_lick_offset_str_bout(is_bout_r,:);
train_data_logic_primSac_onset_str_bout_r = train_data_logic_primSac_onset_str_bout(is_bout_r,:);
train_data_logic_corrSac_onset_str_bout_r = train_data_logic_corrSac_onset_str_bout(is_bout_r,:);
VID_d_tip_100_str_bout_r = VID_d_tip_100_str_bout(is_bout_r,:);
VID_v_tip_100_str_bout_r = VID_v_tip_100_str_bout(is_bout_r,:);
VID_angle_tip_100_str_bout_r = VID_angle_tip_100_str_bout(is_bout_r,:);

% End bout licks
train_data_logic_SS_end_bout = SS_train_aligned(EPHYS_inds_event_bout(inds_end_bout,:));
train_data_logic_CS_end_bout = CS_train_aligned(EPHYS_inds_event_bout(inds_end_bout,:));
train_data_logic_lick_onset_end_bout = EPHYS_lick_onset_train_aligned(EPHYS_inds_event_bout(inds_end_bout,:));
train_data_logic_lick_vmax_end_bout = EPHYS_lick_vmax_train_aligned(EPHYS_inds_event_bout(inds_end_bout,:));
train_data_logic_lick_dmax_end_bout = EPHYS_lick_dmax_train_aligned(EPHYS_inds_event_bout(inds_end_bout,:));
train_data_logic_lick_vmin_end_bout = EPHYS_lick_vmin_train_aligned(EPHYS_inds_event_bout(inds_end_bout,:));
train_data_logic_lick_offset_end_bout = EPHYS_lick_offset_train_aligned(EPHYS_inds_event_bout(inds_end_bout,:));
train_data_logic_primSac_onset_end_bout = EPHYS_primSac_onset_train_aligned(EPHYS_inds_event_bout(inds_end_bout,:));
train_data_logic_corrSac_onset_end_bout = EPHYS_corrSac_onset_train_aligned(EPHYS_inds_event_bout(inds_end_bout,:));
VID_d_tip_100_end_bout = VID_d_tip_100(VID_inds_event_bout(inds_end_bout,:));
VID_v_tip_100_end_bout = VID_v_tip_100(VID_inds_event_bout(inds_end_bout,:));
VID_angle_tip_100_end_bout = VID_angle_tip_100(VID_inds_event_bout(inds_end_bout,:));
% Left bout
if (VID.DLC.FILE.session_type  == 2)
    train_data_logic_SS_end_bout_l = train_data_logic_SS_end_bout(is_bout_l,:);
    train_data_logic_CS_end_bout_l = train_data_logic_CS_end_bout(is_bout_l,:);
    train_data_logic_lick_onset_end_bout_l = train_data_logic_lick_onset_end_bout(is_bout_l,:);
    train_data_logic_lick_vmax_end_bout_l = train_data_logic_lick_vmax_end_bout(is_bout_l,:);
    train_data_logic_lick_dmax_end_bout_l = train_data_logic_lick_dmax_end_bout(is_bout_l,:);
    train_data_logic_lick_vmin_end_bout_l = train_data_logic_lick_vmin_end_bout(is_bout_l,:);
    train_data_logic_lick_offset_end_bout_l = train_data_logic_lick_offset_end_bout(is_bout_l,:);
    train_data_logic_primSac_onset_end_bout_l = train_data_logic_primSac_onset_end_bout(is_bout_l,:);
    train_data_logic_corrSac_onset_end_bout_l = train_data_logic_corrSac_onset_end_bout(is_bout_l,:);
    VID_d_tip_100_end_bout_l = VID_d_tip_100_end_bout(is_bout_l,:);
    VID_v_tip_100_end_bout_l = VID_v_tip_100_end_bout(is_bout_l,:);
    VID_angle_tip_100_end_bout_l = VID_angle_tip_100_end_bout(is_bout_l,:);
end
% Right bout
train_data_logic_SS_end_bout_r = train_data_logic_SS_end_bout(is_bout_r,:);
train_data_logic_CS_end_bout_r = train_data_logic_CS_end_bout(is_bout_r,:);
train_data_logic_lick_onset_end_bout_r = train_data_logic_lick_onset_end_bout(is_bout_r,:);
train_data_logic_lick_vmax_end_bout_r = train_data_logic_lick_vmax_end_bout(is_bout_r,:);
train_data_logic_lick_dmax_end_bout_r = train_data_logic_lick_dmax_end_bout(is_bout_r,:);
train_data_logic_lick_vmin_end_bout_r = train_data_logic_lick_vmin_end_bout(is_bout_r,:);
train_data_logic_lick_offset_end_bout_r = train_data_logic_lick_offset_end_bout(is_bout_r,:);
train_data_logic_primSac_onset_end_bout_r = train_data_logic_primSac_onset_end_bout(is_bout_r,:);
train_data_logic_corrSac_onset_end_bout_r = train_data_logic_corrSac_onset_end_bout(is_bout_r,:);
VID_d_tip_100_end_bout_r = VID_d_tip_100_end_bout(is_bout_r,:);
VID_v_tip_100_end_bout_r = VID_v_tip_100_end_bout(is_bout_r,:);
VID_angle_tip_100_end_bout_r = VID_angle_tip_100_end_bout(is_bout_r,:);

% % time
% raster_data.EPHYS_time_100 = EPHYS_time_100;
% raster_data.EPHYS_time_30K = EPHYS_time_30K;
%
% % SS and CS entire boolean
% raster_data.SS_train_aligned = SS_train_aligned;
% raster_data.CS_train_aligned = CS_train_aligned;

% 80toINF licks
raster_data.train_data_logic_SS_80toINF_0to5 = train_data_logic_SS_80toINF_0to5;
raster_data.train_data_logic_SS_80toINF_5to10 = train_data_logic_SS_80toINF_5to10;
raster_data.train_data_logic_SS_80toINF_10to15 = train_data_logic_SS_80toINF_10to15;
raster_data.train_data_logic_SS_80toINF_15toINF = train_data_logic_SS_80toINF_15toINF;
raster_data.train_data_logic_CS_80toINF_0to5 = train_data_logic_CS_80toINF_0to5;
raster_data.train_data_logic_CS_80toINF_5to10 = train_data_logic_CS_80toINF_5to10;
raster_data.train_data_logic_CS_80toINF_10to15 = train_data_logic_CS_80toINF_10to15;
raster_data.train_data_logic_CS_80toINF_15toINF = train_data_logic_CS_80toINF_15toINF;
raster_data.VID_d_tip_100_80toINF_0to5 = VID_d_tip_100_80toINF_0to5;
raster_data.VID_d_tip_100_80toINF_5to10 = VID_d_tip_100_80toINF_5to10;
raster_data.VID_d_tip_100_80toINF_10to15 = VID_d_tip_100_80toINF_10to15;
raster_data.VID_d_tip_100_80toINF_15toINF = VID_d_tip_100_80toINF_15toINF;
raster_data.VID_v_tip_100_80toINF_0to5 = VID_v_tip_100_80toINF_0to5;
raster_data.VID_v_tip_100_80toINF_5to10 = VID_v_tip_100_80toINF_5to10;
raster_data.VID_v_tip_100_80toINF_10to15 = VID_v_tip_100_80toINF_10to15;
raster_data.VID_v_tip_100_80toINF_15toINF = VID_v_tip_100_80toINF_15toINF;
raster_data.VID_angle_tip_100_80toINF_0to5 = VID_angle_tip_100_80toINF_0to5;
raster_data.VID_angle_tip_100_80toINF_5to10 = VID_angle_tip_100_80toINF_5to10;
raster_data.VID_angle_tip_100_80toINF_10to15 = VID_angle_tip_100_80toINF_10to15;
raster_data.VID_angle_tip_100_80toINF_15toINF = VID_angle_tip_100_80toINF_15toINF;

% 60to80 licks
raster_data.train_data_logic_SS_60to80_0to5 = train_data_logic_SS_60to80_0to5;
raster_data.train_data_logic_SS_60to80_5to10 = train_data_logic_SS_60to80_5to10;
raster_data.train_data_logic_SS_60to80_10to15= train_data_logic_SS_60to80_10to15;
raster_data.train_data_logic_SS_60to80_15toINF = train_data_logic_SS_60to80_15toINF;
raster_data.train_data_logic_CS_60to80_0to5 = train_data_logic_CS_60to80_0to5;
raster_data.train_data_logic_CS_60to80_5to10 = train_data_logic_CS_60to80_5to10;
raster_data.train_data_logic_CS_60to80_10to15 = train_data_logic_CS_60to80_10to15;
raster_data.train_data_logic_CS_60to80_15toINF = train_data_logic_CS_60to80_15toINF;
raster_data.VID_d_tip_100_60to80_0to5 = VID_d_tip_100_60to80_0to5;
raster_data.VID_d_tip_100_60to80_5to10 = VID_d_tip_100_60to80_5to10;
raster_data.VID_d_tip_100_60to80_10to15 = VID_d_tip_100_60to80_10to15;
raster_data.VID_d_tip_100_60to80_15toINF = VID_d_tip_100_60to80_15toINF;
raster_data.VID_v_tip_100_60to80_0to5 = VID_v_tip_100_60to80_0to5;
raster_data.VID_v_tip_100_60to80_5to10 = VID_v_tip_100_60to80_5to10;
raster_data.VID_v_tip_100_60to80_10to15 = VID_v_tip_100_60to80_10to15;
raster_data.VID_v_tip_100_60to80_15toINF = VID_v_tip_100_60to80_15toINF;
raster_data.VID_angle_tip_100_60to80_0to5 = VID_angle_tip_100_60to80_0to5;
raster_data.VID_angle_tip_100_60to80_5to10 = VID_angle_tip_100_60to80_5to10;
raster_data.VID_angle_tip_100_60to80_10to15 = VID_angle_tip_100_60to80_10to15;
raster_data.VID_angle_tip_100_60to80_15toINF = VID_angle_tip_100_60to80_15toINF;

% 20to60 licks
raster_data.train_data_logic_SS_20to60_0to5 = train_data_logic_SS_20to60_0to5;
raster_data.train_data_logic_SS_20to60_5to10 = train_data_logic_SS_20to60_5to10;
raster_data.train_data_logic_SS_20to60_10to15 = train_data_logic_SS_20to60_10to15;
raster_data.train_data_logic_SS_20to60_15toINF = train_data_logic_SS_20to60_15toINF;
raster_data.train_data_logic_CS_20to60_0to5 = train_data_logic_CS_20to60_0to5;
raster_data.train_data_logic_CS_20to60_5to10 = train_data_logic_CS_20to60_5to10;
raster_data.train_data_logic_CS_20to60_10to15 = train_data_logic_CS_20to60_10to15;
raster_data.train_data_logic_CS_20to60_15toINF = train_data_logic_CS_20to60_15toINF;
raster_data.VID_d_tip_100_20to60_0to5 = VID_d_tip_100_20to60_0to5;
raster_data.VID_d_tip_100_20to60_5to10 = VID_d_tip_100_20to60_5to10;
raster_data.VID_d_tip_100_20to60_10to15 = VID_d_tip_100_20to60_10to15;
raster_data.VID_d_tip_100_20to60_15toINF = VID_d_tip_100_20to60_15toINF;
raster_data.VID_v_tip_100_20to60_0to5 = VID_v_tip_100_20to60_0to5;
raster_data.VID_v_tip_100_20to60_5to10 = VID_v_tip_100_20to60_5to10;
raster_data.VID_v_tip_100_20to60_10to15 = VID_v_tip_100_20to60_10to15;
raster_data.VID_v_tip_100_20to60_15toINF = VID_v_tip_100_20to60_15toINF;
raster_data.VID_angle_tip_100_20to60_0to5 = VID_angle_tip_100_20to60_0to5;
raster_data.VID_angle_tip_100_20to60_5to10 = VID_angle_tip_100_20to60_5to10;
raster_data.VID_angle_tip_100_20to60_10to15 = VID_angle_tip_100_20to60_10to15;
raster_data.VID_angle_tip_100_20to60_15toINF = VID_angle_tip_100_20to60_15toINF;

% neg20to20 licks
raster_data.train_data_logic_SS_neg20to20_0to5 = train_data_logic_SS_neg20to20_0to5;
raster_data.train_data_logic_SS_neg20to20_5to10 = train_data_logic_SS_neg20to20_5to10;
raster_data.train_data_logic_SS_neg20to20_10to15 = train_data_logic_SS_neg20to20_10to15;
raster_data.train_data_logic_SS_neg20to20_15toINF = train_data_logic_SS_neg20to20_15toINF;
raster_data.train_data_logic_CS_neg20to20_0to5 = train_data_logic_CS_neg20to20_0to5;
raster_data.train_data_logic_CS_neg20to20_5to10 = train_data_logic_CS_neg20to20_5to10;
raster_data.train_data_logic_CS_neg20to20_10to15 = train_data_logic_CS_neg20to20_10to15;
raster_data.train_data_logic_CS_neg20to20_15toINF = train_data_logic_CS_neg20to20_15toINF;
raster_data.VID_d_tip_100_neg20to20_0to5 = VID_d_tip_100_neg20to20_0to5;
raster_data.VID_d_tip_100_neg20to20_5to10 = VID_d_tip_100_neg20to20_5to10;
raster_data.VID_d_tip_100_neg20to20_10to15 = VID_d_tip_100_neg20to20_10to15;
raster_data.VID_d_tip_100_neg20to20_15toINF = VID_d_tip_100_neg20to20_15toINF;
raster_data.VID_v_tip_100_neg20to20_0to5 = VID_v_tip_100_neg20to20_0to5;
raster_data.VID_v_tip_100_neg20to20_5to10 = VID_v_tip_100_neg20to20_5to10;
raster_data.VID_v_tip_100_neg20to20_10to15 = VID_v_tip_100_neg20to20_10to15;
raster_data.VID_v_tip_100_neg20to20_15toINF = VID_v_tip_100_neg20to20_15toINF;
raster_data.VID_angle_tip_100_neg20to20_0to5 = VID_angle_tip_100_neg20to20_0to5;
raster_data.VID_angle_tip_100_neg20to20_5to10 = VID_angle_tip_100_neg20to20_5to10;
raster_data.VID_angle_tip_100_neg20to20_10to15 = VID_angle_tip_100_neg20to20_10to15;
raster_data.VID_angle_tip_100_neg20to20_15toINF = VID_angle_tip_100_neg20to20_15toINF;

% neg20toneg60 licks
raster_data.train_data_logic_SS_neg20toneg60_0to5 = train_data_logic_SS_neg20toneg60_0to5;
raster_data.train_data_logic_SS_neg20toneg60_5to10 = train_data_logic_SS_neg20toneg60_5to10;
raster_data.train_data_logic_SS_neg20toneg60_10to15 = train_data_logic_SS_neg20toneg60_10to15;
raster_data.train_data_logic_SS_neg20toneg60_15toINF = train_data_logic_SS_neg20toneg60_15toINF;
raster_data.train_data_logic_CS_neg20toneg60_0to5 = train_data_logic_CS_neg20toneg60_0to5;
raster_data.train_data_logic_CS_neg20toneg60_5to10 = train_data_logic_CS_neg20toneg60_5to10;
raster_data.train_data_logic_CS_neg20toneg60_10to15 = train_data_logic_CS_neg20toneg60_10to15;
raster_data.train_data_logic_CS_neg20toneg60_15toINF = train_data_logic_CS_neg20toneg60_15toINF;
raster_data.VID_d_tip_100_neg20toneg60_0to5 = VID_d_tip_100_neg20toneg60_0to5;
raster_data.VID_d_tip_100_neg20toneg60_5to10 = VID_d_tip_100_neg20toneg60_5to10;
raster_data.VID_d_tip_100_neg20toneg60_10to15 = VID_d_tip_100_neg20toneg60_10to15;
raster_data.VID_d_tip_100_neg20toneg60_15toINF = VID_d_tip_100_neg20toneg60_15toINF;
raster_data.VID_v_tip_100_neg20toneg60_0to5 = VID_v_tip_100_neg20toneg60_0to5;
raster_data.VID_v_tip_100_neg20toneg60_5to10 = VID_v_tip_100_neg20toneg60_5to10;
raster_data.VID_v_tip_100_neg20toneg60_10to15 = VID_v_tip_100_neg20toneg60_10to15;
raster_data.VID_v_tip_100_neg20toneg60_15toINF = VID_v_tip_100_neg20toneg60_15toINF;
raster_data.VID_angle_tip_100_neg20toneg60_0to5 = VID_angle_tip_100_neg20toneg60_0to5;
raster_data.VID_angle_tip_100_neg20toneg60_5to10 = VID_angle_tip_100_neg20toneg60_5to10;
raster_data.VID_angle_tip_100_neg20toneg60_10to15 = VID_angle_tip_100_neg20toneg60_10to15;
raster_data.VID_angle_tip_100_neg20toneg60_15toINF = VID_angle_tip_100_neg20toneg60_15toINF;

% neg60toneg80 licks
raster_data.train_data_logic_SS_neg60toneg80_0to5 = train_data_logic_SS_neg60toneg80_0to5;
raster_data.train_data_logic_SS_neg60toneg80_5to10 = train_data_logic_SS_neg60toneg80_5to10;
raster_data.train_data_logic_SS_neg60toneg80_10to15 = train_data_logic_SS_neg60toneg80_10to15;
raster_data.train_data_logic_SS_neg60toneg80_15toINF = train_data_logic_SS_neg60toneg80_15toINF;
raster_data.train_data_logic_CS_neg60toneg80_0to5 = train_data_logic_CS_neg60toneg80_0to5;
raster_data.train_data_logic_CS_neg60toneg80_5to10 = train_data_logic_CS_neg60toneg80_5to10;
raster_data.train_data_logic_CS_neg60toneg80_10to15 = train_data_logic_CS_neg60toneg80_10to15;
raster_data.train_data_logic_CS_neg60toneg80_15toINF = train_data_logic_CS_neg60toneg80_15toINF;
raster_data.VID_d_tip_100_neg60toneg80_0to5 = VID_d_tip_100_neg60toneg80_0to5;
raster_data.VID_d_tip_100_neg60toneg80_5to10 = VID_d_tip_100_neg60toneg80_5to10;
raster_data.VID_d_tip_100_neg60toneg80_10to15 = VID_d_tip_100_neg60toneg80_10to15;
raster_data.VID_d_tip_100_neg60toneg80_15toINF = VID_d_tip_100_neg60toneg80_15toINF;
raster_data.VID_v_tip_100_neg60toneg80_0to5 = VID_v_tip_100_neg60toneg80_0to5;
raster_data.VID_v_tip_100_neg60toneg80_5to10 = VID_v_tip_100_neg60toneg80_5to10;
raster_data.VID_v_tip_100_neg60toneg80_10to15 = VID_v_tip_100_neg60toneg80_10to15;
raster_data.VID_v_tip_100_neg60toneg80_15toINF = VID_v_tip_100_neg60toneg80_15toINF;
raster_data.VID_angle_tip_100_neg60toneg80_0to5 = VID_angle_tip_100_neg60toneg80_0to5;
raster_data.VID_angle_tip_100_neg60toneg80_5to10 = VID_angle_tip_100_neg60toneg80_5to10;
raster_data.VID_angle_tip_100_neg60toneg80_10to15 = VID_angle_tip_100_neg60toneg80_10to15;
raster_data.VID_angle_tip_100_neg60toneg80_15toINF = VID_angle_tip_100_neg60toneg80_15toINF;

% neg80tonegINF licks
raster_data.train_data_logic_SS_neg80tonegINF_0to5 = train_data_logic_SS_neg80tonegINF_0to5;
raster_data.train_data_logic_SS_neg80tonegINF_5to10 = train_data_logic_SS_neg80tonegINF_5to10;
raster_data.train_data_logic_SS_neg80tonegINF_10to15 = train_data_logic_SS_neg80tonegINF_10to15;
raster_data.train_data_logic_SS_neg80tonegINF_15toINF = train_data_logic_SS_neg80tonegINF_15toINF;
raster_data.train_data_logic_CS_neg80tonegINF_0to5 = train_data_logic_CS_neg80tonegINF_0to5;
raster_data.train_data_logic_CS_neg80tonegINF_5to10 = train_data_logic_CS_neg80tonegINF_5to10;
raster_data.train_data_logic_CS_neg80tonegINF_10to15 = train_data_logic_CS_neg80tonegINF_10to15;
raster_data.train_data_logic_CS_neg80tonegINF_15toINF = train_data_logic_CS_neg80tonegINF_15toINF;
raster_data.VID_d_tip_100_neg80tonegINF_0to5 = VID_d_tip_100_neg80tonegINF_0to5;
raster_data.VID_d_tip_100_neg80tonegINF_5to10 = VID_d_tip_100_neg80tonegINF_5to10;
raster_data.VID_d_tip_100_neg80tonegINF_10to15 = VID_d_tip_100_neg80tonegINF_10to15;
raster_data.VID_d_tip_100_neg80tonegINF_15toINF = VID_d_tip_100_neg80tonegINF_15toINF;
raster_data.VID_v_tip_100_neg80tonegINF_0to5 = VID_v_tip_100_neg80tonegINF_0to5;
raster_data.VID_v_tip_100_neg80tonegINF_5to10 = VID_v_tip_100_neg80tonegINF_5to10;
raster_data.VID_v_tip_100_neg80tonegINF_10to15 = VID_v_tip_100_neg80tonegINF_10to15;
raster_data.VID_v_tip_100_neg80tonegINF_15toINF = VID_v_tip_100_neg80tonegINF_15toINF;
raster_data.VID_angle_tip_100_neg80tonegINF_0to5 = VID_angle_tip_100_neg80tonegINF_0to5;
raster_data.VID_angle_tip_100_neg80tonegINF_5to10 = VID_angle_tip_100_neg80tonegINF_5to10;
raster_data.VID_angle_tip_100_neg80tonegINF_10to15 = VID_angle_tip_100_neg80tonegINF_10to15;
raster_data.VID_angle_tip_100_neg80tonegINF_15toINF = VID_angle_tip_100_neg80tonegINF_15toINF;

% All licks
raster_data.train_data_logic_SS_all_corr = train_data_logic_SS_all_corr ;
raster_data.train_data_logic_CS_all_corr = train_data_logic_CS_all_corr ;
raster_data.VID_d_tip_100_all_corr  = VID_d_tip_100_all_corr;
raster_data.VID_v_tip_100_all_corr  = VID_v_tip_100_all_corr;
raster_data.VID_angle_tip_100_all_corr  = VID_angle_tip_100_all_corr;
raster_data.train_data_logic_lick_onset_all_corr = train_data_logic_lick_onset_all_corr ;
raster_data.train_data_logic_lick_vmax_all_corr = train_data_logic_lick_vmax_all_corr ;
raster_data.train_data_logic_lick_dmax_all_corr = train_data_logic_lick_dmax_all_corr ;
raster_data.train_data_logic_lick_vmin_all_corr = train_data_logic_lick_vmin_all_corr ;
raster_data.train_data_logic_lick_offset_all_corr = train_data_logic_lick_offset_all_corr ;
raster_data.train_data_logic_SS_all = train_data_logic_SS_all;
raster_data.train_data_logic_CS_all = train_data_logic_CS_all;
raster_data.train_data_logic_lick_onset_all    = train_data_logic_lick_onset_all;
raster_data.train_data_logic_lick_vmax_all  = train_data_logic_lick_vmax_all;
raster_data.train_data_logic_lick_dmax_all   = train_data_logic_lick_dmax_all;
raster_data.train_data_logic_lick_vmin_all = train_data_logic_lick_vmin_all;
raster_data.train_data_logic_lick_offset_all  = train_data_logic_lick_offset_all;
raster_data.VID_d_tip_100_all  = VID_d_tip_100_all;
raster_data.VID_v_tip_100_all  = VID_v_tip_100_all;
raster_data.VID_angle_tip_100_all  = VID_angle_tip_100_all;
raster_data.train_data_logic_primSac_onset_all = train_data_logic_primSac_onset_all;
raster_data.train_data_logic_corrSac_onset_all = train_data_logic_corrSac_onset_all;

% Grooming licks
raster_data.train_data_logic_SS_grooming_corr = train_data_logic_SS_grooming_corr ;
raster_data.train_data_logic_CS_grooming_corr = train_data_logic_CS_grooming_corr ;
raster_data.VID_d_tip_100_grooming_corr  = VID_d_tip_100_grooming_corr;
raster_data.VID_v_tip_100_grooming_corr  = VID_v_tip_100_grooming_corr;
raster_data.VID_angle_tip_100_grooming_corr  = VID_angle_tip_100_grooming_corr;
raster_data.train_data_logic_lick_onset_grooming_corr = train_data_logic_lick_onset_grooming_corr ;
raster_data.train_data_logic_lick_vmax_grooming_corr = train_data_logic_lick_vmax_grooming_corr ;
raster_data.train_data_logic_lick_dmax_grooming_corr = train_data_logic_lick_dmax_grooming_corr ;
raster_data.train_data_logic_lick_vmin_grooming_corr = train_data_logic_lick_vmin_grooming_corr ;
raster_data.train_data_logic_lick_offset_grooming_corr = train_data_logic_lick_offset_grooming_corr ;
raster_data.train_data_logic_SS_grooming = train_data_logic_SS_grooming;
% raster_data.train_data_logic_SS_grooming_510 = train_data_logic_SS_grooming_510;
% raster_data.train_data_logic_SS_grooming_10 = train_data_logic_SS_grooming_10;
raster_data.train_data_logic_CS_grooming = train_data_logic_CS_grooming;
raster_data.train_data_logic_lick_onset_grooming    = train_data_logic_lick_onset_grooming;
raster_data.train_data_logic_lick_vmax_grooming  = train_data_logic_lick_vmax_grooming;
raster_data.train_data_logic_lick_dmax_grooming   = train_data_logic_lick_dmax_grooming;
raster_data.train_data_logic_lick_vmin_grooming = train_data_logic_lick_vmin_grooming;
raster_data.train_data_logic_lick_offset_grooming  = train_data_logic_lick_offset_grooming;
raster_data.VID_d_tip_100_grooming  = VID_d_tip_100_grooming;
raster_data.VID_v_tip_100_grooming  = VID_v_tip_100_grooming;
raster_data.VID_angle_tip_100_grooming  = VID_angle_tip_100_grooming;
%     raster_data.BEHAVE_eye_r_vm_grooming = BEHAVE_eye_r_vm_grooming;
raster_data.train_data_logic_primSac_onset_grooming = train_data_logic_primSac_onset_grooming;
raster_data.train_data_logic_corrSac_onset_grooming = train_data_logic_corrSac_onset_grooming;

% R reward licks
raster_data.train_data_logic_SS_r_corr = train_data_logic_SS_r_corr ;
raster_data.train_data_logic_CS_r_corr = train_data_logic_CS_r_corr ;
raster_data.VID_d_tip_100_r_corr  = VID_d_tip_100_r_corr;
raster_data.VID_v_tip_100_r_corr  = VID_v_tip_100_r_corr;
raster_data.VID_angle_tip_100_r_corr  = VID_angle_tip_100_r_corr;
raster_data.train_data_logic_lick_onset_r_corr = train_data_logic_lick_onset_r_corr ;
raster_data.train_data_logic_lick_vmax_r_corr = train_data_logic_lick_vmax_r_corr ;
raster_data.train_data_logic_lick_dmax_r_corr = train_data_logic_lick_dmax_r_corr ;
raster_data.train_data_logic_lick_vmin_r_corr = train_data_logic_lick_vmin_r_corr ;
raster_data.train_data_logic_lick_offset_r_corr = train_data_logic_lick_offset_r_corr ;
raster_data.train_data_logic_SS_r = train_data_logic_SS_r;
% raster_data.train_data_logic_SS_r_15 = train_data_logic_SS_r_15;
% raster_data.train_data_logic_SS_r_1520 = train_data_logic_SS_r_1520;
raster_data.train_data_logic_CS_r = train_data_logic_CS_r;
raster_data.train_data_logic_CS_r_1 = train_data_logic_CS_r_1;
raster_data.train_data_logic_CS_r_0 = train_data_logic_CS_r_0;
raster_data.train_data_logic_lick_onset_r    = train_data_logic_lick_onset_r;
raster_data.train_data_logic_lick_vmax_r  = train_data_logic_lick_vmax_r;
raster_data.train_data_logic_lick_dmax_r   = train_data_logic_lick_dmax_r;
raster_data.train_data_logic_lick_vmin_r = train_data_logic_lick_vmin_r;
raster_data.train_data_logic_lick_offset_r  = train_data_logic_lick_offset_r;
raster_data.VID_d_tip_100_r  = VID_d_tip_100_r;
raster_data.VID_v_tip_100_r  = VID_v_tip_100_r;
raster_data.VID_angle_tip_100_r  = VID_angle_tip_100_r;
%     raster_data.BEHAVE_eye_r_vm_r = BEHAVE_eye_r_vm_r;
raster_data.train_data_logic_primSac_onset_r = train_data_logic_primSac_onset_r;
raster_data.train_data_logic_corrSac_onset_r = train_data_logic_corrSac_onset_r;

% L reward licks
if (VID.DLC.FILE.session_type  == 2)
    raster_data.train_data_logic_SS_l_corr = train_data_logic_SS_l_corr ;
    raster_data.train_data_logic_CS_l_corr = train_data_logic_CS_l_corr ;
    raster_data.VID_d_tip_100_l_corr  = VID_d_tip_100_l_corr;
    raster_data.VID_v_tip_100_l_corr  = VID_v_tip_100_l_corr;
    raster_data.VID_angle_tip_100_l_corr  = VID_angle_tip_100_l_corr;
    raster_data.train_data_logic_lick_onset_l_corr = train_data_logic_lick_onset_l_corr ;
    raster_data.train_data_logic_lick_vmax_l_corr = train_data_logic_lick_vmax_l_corr ;
    raster_data.train_data_logic_lick_dmax_l_corr = train_data_logic_lick_dmax_l_corr ;
    raster_data.train_data_logic_lick_vmin_l_corr = train_data_logic_lick_vmin_l_corr ;
    raster_data.train_data_logic_lick_offset_l_corr = train_data_logic_lick_offset_l_corr ;
    raster_data.train_data_logic_SS_l = train_data_logic_SS_l;
    %     raster_data.train_data_logic_SS_l_15 = train_data_logic_SS_l_15;
    %     raster_data.train_data_logic_SS_l_1520 = train_data_logic_SS_l_1520;
    raster_data.train_data_logic_CS_l = train_data_logic_CS_l;
    raster_data.train_data_logic_CS_l_1 = train_data_logic_CS_l_1;
    raster_data.train_data_logic_CS_l_0 = train_data_logic_CS_l_0;
    raster_data.train_data_logic_lick_onset_l    = train_data_logic_lick_onset_l;
    raster_data.train_data_logic_lick_vmax_l  = train_data_logic_lick_vmax_l;
    raster_data.train_data_logic_lick_dmax_l   = train_data_logic_lick_dmax_l;
    raster_data.train_data_logic_lick_vmin_l = train_data_logic_lick_vmin_l;
    raster_data.train_data_logic_lick_offset_l  = train_data_logic_lick_offset_l;
    raster_data.VID_d_tip_100_l  = VID_d_tip_100_l;
    raster_data.VID_v_tip_100_l  = VID_v_tip_100_l;
    raster_data.VID_angle_tip_100_l  = VID_angle_tip_100_l;
    %         raster_data.BEHAVE_eye_r_vm_l = BEHAVE_eye_r_vm_l;
    raster_data.train_data_logic_primSac_onset_l = train_data_logic_primSac_onset_l;
    raster_data.train_data_logic_corrSac_onset_l = train_data_logic_corrSac_onset_l;
end

% Str bout licks
raster_data.train_data_logic_SS_str_bout = train_data_logic_SS_str_bout;
raster_data.train_data_logic_CS_str_bout = train_data_logic_CS_str_bout;
raster_data.train_data_logic_lick_onset_str_bout    = train_data_logic_lick_onset_str_bout;
raster_data.train_data_logic_lick_vmax_str_bout  = train_data_logic_lick_vmax_str_bout;
raster_data.train_data_logic_lick_dmax_str_bout   = train_data_logic_lick_dmax_str_bout;
raster_data.train_data_logic_lick_vmin_str_bout = train_data_logic_lick_vmin_str_bout;
raster_data.train_data_logic_lick_offset_str_bout  = train_data_logic_lick_offset_str_bout;
raster_data.train_data_logic_primSac_onset_str_bout = train_data_logic_primSac_onset_str_bout;
raster_data.train_data_logic_corrSac_onset_str_bout = train_data_logic_corrSac_onset_str_bout;
raster_data.VID_d_tip_100_str_bout  = VID_d_tip_100_str_bout;
raster_data.VID_v_tip_100_str_bout  = VID_v_tip_100_str_bout;
raster_data.VID_angle_tip_100_str_bout  = VID_angle_tip_100_str_bout;
% Left bout
if (VID.DLC.FILE.session_type  == 2)
    raster_data.train_data_logic_SS_str_bout_l = train_data_logic_SS_str_bout_l;
    raster_data.train_data_logic_CS_str_bout_l = train_data_logic_CS_str_bout_l;
    raster_data.train_data_logic_lick_onset_str_bout_l    = train_data_logic_lick_onset_str_bout_l;
    raster_data.train_data_logic_lick_vmax_str_bout_l  = train_data_logic_lick_vmax_str_bout_l;
    raster_data.train_data_logic_lick_dmax_str_bout_l   = train_data_logic_lick_dmax_str_bout_l;
    raster_data.train_data_logic_lick_vmin_str_bout_l = train_data_logic_lick_vmin_str_bout_l;
    raster_data.train_data_logic_lick_offset_str_bout_l = train_data_logic_lick_offset_str_bout_l;
    raster_data.train_data_logic_primSac_onset_str_bout_l = train_data_logic_primSac_onset_str_bout_l;
    raster_data.train_data_logic_corrSac_onset_str_bout_l = train_data_logic_corrSac_onset_str_bout_l;
    raster_data.VID_d_tip_100_str_bout_l  = VID_d_tip_100_str_bout_l;
    raster_data.VID_v_tip_100_str_bout_l  = VID_v_tip_100_str_bout_l;
    raster_data.VID_angle_tip_100_str_bout_l  = VID_angle_tip_100_str_bout_l;
end
% Right bout
raster_data.train_data_logic_SS_str_bout_r = train_data_logic_SS_str_bout_r;
raster_data.train_data_logic_CS_str_bout_r = train_data_logic_CS_str_bout_r;
raster_data.train_data_logic_lick_onset_str_bout_r    = train_data_logic_lick_onset_str_bout_r;
raster_data.train_data_logic_lick_vmax_str_bout_r = train_data_logic_lick_vmax_str_bout_r;
raster_data.train_data_logic_lick_dmax_str_bout_r   = train_data_logic_lick_dmax_str_bout_r;
raster_data.train_data_logic_lick_vmin_str_bout_r = train_data_logic_lick_vmin_str_bout_r;
raster_data.train_data_logic_lick_offset_str_bout_r = train_data_logic_lick_offset_str_bout_r;
raster_data.train_data_logic_primSac_onset_str_bout_r = train_data_logic_primSac_onset_str_bout_r;
raster_data.train_data_logic_corrSac_onset_str_bout_r = train_data_logic_corrSac_onset_str_bout_r;
raster_data.VID_d_tip_100_str_bout_r  = VID_d_tip_100_str_bout_r;
raster_data.VID_v_tip_100_str_bout_r  = VID_v_tip_100_str_bout_r;
raster_data.VID_angle_tip_100_str_bout_r  = VID_angle_tip_100_str_bout_r;

% End bout licks
raster_data.train_data_logic_SS_end_bout = train_data_logic_SS_end_bout;
raster_data.train_data_logic_CS_end_bout = train_data_logic_CS_end_bout;
raster_data.train_data_logic_lick_onset_end_bout    = train_data_logic_lick_onset_end_bout;
raster_data.train_data_logic_lick_vmax_end_bout  = train_data_logic_lick_vmax_end_bout;
raster_data.train_data_logic_lick_dmax_end_bout   = train_data_logic_lick_dmax_end_bout;
raster_data.train_data_logic_lick_vmin_end_bout = train_data_logic_lick_vmin_end_bout;
raster_data.train_data_logic_lick_offset_end_bout  = train_data_logic_lick_offset_end_bout;
raster_data.train_data_logic_primSac_onset_end_bout = train_data_logic_primSac_onset_end_bout;
raster_data.train_data_logic_corrSac_onset_end_bout = train_data_logic_corrSac_onset_end_bout;
raster_data.VID_d_tip_100_end_bout  = VID_d_tip_100_end_bout;
raster_data.VID_v_tip_100_end_bout  = VID_v_tip_100_end_bout;
raster_data.VID_angle_tip_100_end_bout  = VID_angle_tip_100_end_bout;
% Left bout
if (VID.DLC.FILE.session_type  == 2)
    raster_data.train_data_logic_SS_end_bout_l = train_data_logic_SS_end_bout_l;
    raster_data.train_data_logic_CS_end_bout_l = train_data_logic_CS_end_bout_l;
    raster_data.train_data_logic_lick_onset_end_bout_l = train_data_logic_lick_onset_end_bout_l;
    raster_data.train_data_logic_lick_vmax_end_bout_l = train_data_logic_lick_vmax_end_bout_l;
    raster_data.train_data_logic_lick_dmax_end_bout_l   = train_data_logic_lick_dmax_end_bout_l;
    raster_data.train_data_logic_lick_vmin_end_bout_l = train_data_logic_lick_vmin_end_bout_l;
    raster_data.train_data_logic_lick_offset_end_bout_l  = train_data_logic_lick_offset_end_bout_l;
    raster_data.train_data_logic_primSac_onset_end_bout_l = train_data_logic_primSac_onset_end_bout_l;
    raster_data.train_data_logic_corrSac_onset_end_bout_l = train_data_logic_corrSac_onset_end_bout_l;
    raster_data.VID_d_tip_100_end_bout_l  = VID_d_tip_100_end_bout_l;
    raster_data.VID_v_tip_100_end_bout_l  = VID_v_tip_100_end_bout_l;
    raster_data.VID_angle_tip_100_end_bout_l  = VID_angle_tip_100_end_bout_l;
end
% Right bout
raster_data.train_data_logic_SS_end_bout_r = train_data_logic_SS_end_bout_r;
raster_data.train_data_logic_CS_end_bout_r = train_data_logic_CS_end_bout_r;
raster_data.train_data_logic_lick_onset_end_bout_r = train_data_logic_lick_onset_end_bout_r;
raster_data.train_data_logic_lick_vmax_end_bout_r = train_data_logic_lick_vmax_end_bout_r;
raster_data.train_data_logic_lick_dmax_end_bout_r   = train_data_logic_lick_dmax_end_bout_r;
raster_data.train_data_logic_lick_vmin_end_bout_r = train_data_logic_lick_vmin_end_bout_r;
raster_data.train_data_logic_lick_offset_end_bout_r  = train_data_logic_lick_offset_end_bout_r;
raster_data.train_data_logic_primSac_onset_end_bout_r = train_data_logic_primSac_onset_end_bout_r;
raster_data.train_data_logic_corrSac_onset_end_bout_r = train_data_logic_corrSac_onset_end_bout_r;
raster_data.VID_d_tip_100_end_bout_r  = VID_d_tip_100_end_bout_r;
raster_data.VID_v_tip_100_end_bout_r  = VID_v_tip_100_end_bout_r;
raster_data.VID_angle_tip_100_end_bout_r  = VID_angle_tip_100_end_bout_r;

end

%% function single_dataset_kinematic
function kinematic_data = single_dataset_kinematic(EPHYS, VID, session_type)

% % Kinematics
VID_d_tip_100 = EPHYS.CH_EVE.VID_d_tip_100;
VID_v_tip_100 = EPHYS.CH_EVE.VID_v_tip_100;
VID_angle_tip_100 = EPHYS.CH_EVE.VID_angle_tip_100;
VID_time_100 = EPHYS.CH_EVE.VID_time_100;
VID_lick_duration = VID.DLC.TIME.time_lick_duration;
VID_time_lick_onset = VID.DLC.TIME.time_lick_onset;

ind_lick_onset_100  = EPHYS.CH_EVE.VID_LED_aligned_ind_lick_onset_100 ;
ind_lick_dmax_100  = EPHYS.CH_EVE.VID_LED_aligned_ind_lick_dmax_100 ;
ind_lick_vmax_100  = EPHYS.CH_EVE.VID_LED_aligned_ind_lick_vmax_100 ;
ind_lick_vmin_100  = EPHYS.CH_EVE.VID_LED_aligned_ind_lick_vmin_100 ;
ind_lick_anglemax_100  = EPHYS.CH_EVE.VID_LED_aligned_ind_lick_dmax_100;
ind_lick_offset_100  = EPHYS.CH_EVE.VID_LED_aligned_ind_lick_offset_100 ;

num_licks = length(ind_lick_dmax_100);

% inds class of lick
is_grooming = logical(VID.DLC.CLASS.is_grooming_lick);
is_grooming = is_grooming(1:num_licks);
is_r = logical(VID.DLC.CLASS.is_r_reward_inner_tube_lick + VID.DLC.CLASS.is_r_reward_outer_tube_lick + VID.DLC.CLASS.is_r_noreward_inner_tube_lick + VID.DLC.CLASS.is_r_noreward_outer_tube_lick);
is_r = is_r(1:num_licks);
is_r_1 = logical(VID.DLC.CLASS.is_r_reward_inner_tube_lick + VID.DLC.CLASS.is_r_reward_outer_tube_lick);
is_r_1 = is_r_1(1:num_licks);
is_r_0 =  logical(VID.DLC.CLASS.is_r_noreward_inner_tube_lick + VID.DLC.CLASS.is_r_noreward_outer_tube_lick);
is_r_0 = is_r_0(1:num_licks);
if (session_type  == 2)
    is_l =  logical(VID.DLC.CLASS.is_l_reward_inner_tube_lick + VID.DLC.CLASS.is_l_reward_outer_tube_lick + VID.DLC.CLASS.is_l_noreward_inner_tube_lick + VID.DLC.CLASS.is_l_noreward_outer_tube_lick);
    is_l = is_l(1:num_licks);
    inds_l_1 = logical( VID.DLC.CLASS.is_l_reward_inner_tube_lick + VID.DLC.CLASS.is_l_reward_outer_tube_lick);
    inds_l_1 = inds_l_1(1:num_licks);
    inds_l_0 = logical(VID.DLC.CLASS.is_l_noreward_inner_tube_lick + VID.DLC.CLASS.is_l_noreward_outer_tube_lick);
    inds_l_0 = inds_l_0(1:num_licks);
end
ind_str_bout = ismember(VID.DLC.IND.ind_lick_onset,VID.DLC.IND.ind_lick_onset_str_bout);
ind_str_bout = ind_str_bout(1:num_licks);
ind_end_bout = ismember(VID.DLC.IND.ind_lick_onset,VID.DLC.IND.ind_lick_onset_end_bout);
ind_end_bout = ind_end_bout(1:num_licks);

%inds class of bout
if (session_type  == 2)
    is_bout_l = logical(VID.DLC.CLASS.is_bout_l(1:sum(ind_str_bout)));
end
is_bout_r = logical(VID.DLC.CLASS.is_bout_r(1:sum(ind_str_bout)));

% Bout
VID_num_lick_bout_all = VID.DLC.IND.num_lick_bout;
VID_num_lick_bout_r = VID_num_lick_bout_all(is_bout_r);
VID_num_lick_bout_l = VID_num_lick_bout_all(is_bout_l);

VID_bout_duration_all = VID.DLC.TIME.time_bout_duration;
VID_bout_duration_r = VID_bout_duration_all(is_bout_r);
VID_bout_duration_l = VID_bout_duration_all(is_bout_l);

% All licks
VID_d_max_100_all = VID_d_tip_100(ind_lick_dmax_100);
VID_v_max_100_all = VID_v_tip_100(ind_lick_vmax_100);
VID_v_min_100_all = VID_v_tip_100(ind_lick_vmin_100);
VID_angle_max_100_all = VID_angle_tip_100(ind_lick_anglemax_100);

VID_duration_d_max_100_all = VID_time_100(ind_lick_dmax_100) - VID_time_100(ind_lick_onset_100);
VID_duration_v_max_100_all = VID_time_100(ind_lick_vmax_100)  - VID_time_100(ind_lick_onset_100);
VID_duration_v_min_100_all = VID_time_100(ind_lick_vmin_100)  - VID_time_100(ind_lick_onset_100);

VID_lick_duration_all = VID_lick_duration;

% Grooming licks
VID_d_max_100_grooming = VID_d_tip_100(ind_lick_dmax_100(is_grooming,:));
VID_v_max_100_grooming = VID_v_tip_100(ind_lick_vmax_100(is_grooming,:));
VID_v_min_100_grooming = VID_v_tip_100(ind_lick_vmin_100(is_grooming,:));
VID_angle_max_100_grooming = VID_angle_tip_100(ind_lick_anglemax_100(is_grooming,:));

VID_duration_d_max_100_grooming = VID_time_100(ind_lick_dmax_100(is_grooming,:)) - VID_time_100(ind_lick_onset_100(is_grooming,:));
VID_duration_v_max_100_grooming = VID_time_100(ind_lick_vmax_100(is_grooming,:)) - VID_time_100(ind_lick_onset_100(is_grooming,:));
VID_duration_v_min_100_grooming = VID_time_100(ind_lick_vmin_100(is_grooming,:)) - VID_time_100(ind_lick_onset_100(is_grooming,:));

VID_lick_duration_grooming = VID_lick_duration(is_grooming);

% R licks
VID_d_max_100_r = VID_d_tip_100(ind_lick_dmax_100(is_r,:));
VID_v_max_100_r = VID_v_tip_100(ind_lick_vmax_100(is_r,:));
VID_v_min_100_r = VID_v_tip_100(ind_lick_vmin_100(is_r,:));
VID_angle_max_100_r = VID_angle_tip_100(ind_lick_anglemax_100(is_r,:));

VID_duration_d_max_100_r = VID_time_100(ind_lick_dmax_100(is_r,:))- VID_time_100(ind_lick_onset_100(is_r,:));
VID_duration_v_max_100_r = VID_time_100(ind_lick_vmax_100(is_r,:))- VID_time_100(ind_lick_onset_100(is_r,:));
VID_duration_v_min_100_r = VID_time_100(ind_lick_vmin_100(is_r,:))- VID_time_100(ind_lick_onset_100(is_r,:));

VID_lick_duration_r = VID_lick_duration(is_r);

% L licks
if (VID.DLC.FILE.session_type  == 2)
    VID_d_max_100_l = VID_d_tip_100(ind_lick_dmax_100(is_l,:));
    VID_v_max_100_l = VID_v_tip_100(ind_lick_vmax_100(is_l,:));
    VID_v_min_100_l = VID_v_tip_100(ind_lick_vmin_100(is_l,:));
    VID_angle_max_100_l = VID_angle_tip_100(ind_lick_anglemax_100(is_l,:));
    
    VID_duration_d_max_100_l = VID_time_100(ind_lick_dmax_100(is_l,:))- VID_time_100(ind_lick_onset_100(is_l,:));
    VID_duration_v_max_100_l = VID_time_100(ind_lick_vmax_100(is_l,:))- VID_time_100(ind_lick_onset_100(is_l,:));
    VID_duration_v_min_100_l = VID_time_100(ind_lick_vmin_100(is_l,:))- VID_time_100(ind_lick_onset_100(is_l,:));
    
    VID_lick_duration_l = VID_lick_duration(is_l);
    
end

% ILI & IFR
% all
VID_ILI_all = [diff(VID_time_lick_onset(1:num_licks)); nan];
% VID_ILI(VID_ILI > 500) = nan;
VID_ILR_all = 1./VID_ILI_all;

% grooming
VID_ILI_grooming = VID_ILI_all(is_grooming);
VID_ILR_grooming = VID_ILR_all(is_grooming);
% r
VID_ILI_r = VID_ILI_all(is_r);
VID_ILR_r = VID_ILR_all(is_r);
% l
VID_ILI_l = VID_ILI_all(is_l);
VID_ILR_l = VID_ILR_all(is_l);


kinematic_data.is_bout_l = is_bout_l;
kinematic_data.is_bout_r = is_bout_r;

kinematic_data.VID_ILI_all = VID_ILI_all;
kinematic_data.VID_ILI_grooming = VID_ILI_grooming;
kinematic_data.VID_ILI_r = VID_ILI_r;
kinematic_data.VID_ILI_l = VID_ILI_l;

kinematic_data.VID_ILR_all = VID_ILR_all;
kinematic_data.VID_ILR_grooming = VID_ILR_grooming;
kinematic_data.VID_ILR_r = VID_ILR_r;
kinematic_data.VID_ILR_l = VID_ILR_l;

kinematic_data.VID_num_lick_bout = VID_num_lick_bout_all;
kinematic_data.VID_bout_duration = VID_bout_duration_all;

kinematic_data.VID_time_100 = VID_time_100;
kinematic_data.VID_d_max_100_all = VID_d_max_100_all;
kinematic_data.VID_v_max_100_all = VID_v_max_100_all;
kinematic_data.VID_v_min_100_all = VID_v_min_100_all;
kinematic_data.VID_angle_max_100_all = VID_angle_max_100_all;
kinematic_data.VID_duration_d_max_100_all = VID_duration_d_max_100_all;
kinematic_data.VID_duration_v_max_100_all = VID_duration_v_max_100_all;
kinematic_data.VID_duration_v_min_100_all = VID_duration_v_min_100_all;
kinematic_data.VID_lick_duration_all = VID_lick_duration_all;

kinematic_data.VID_d_max_100_grooming = VID_d_max_100_grooming;
kinematic_data.VID_v_max_100_grooming = VID_v_max_100_grooming;
kinematic_data.VID_v_min_100_grooming = VID_v_min_100_grooming;
kinematic_data.VID_angle_max_100_grooming = VID_angle_max_100_grooming;
kinematic_data.VID_duration_d_max_100_grooming = VID_duration_d_max_100_grooming;
kinematic_data.VID_duration_v_max_100_grooming = VID_duration_v_max_100_grooming;
kinematic_data.VID_duration_v_min_100_grooming = VID_duration_v_min_100_grooming;
kinematic_data.VID_lick_duration_grooming = VID_lick_duration_grooming;

kinematic_data.VID_d_max_100_r = VID_d_max_100_r;
kinematic_data.VID_v_max_100_r = VID_v_max_100_r;
kinematic_data.VID_v_min_100_r = VID_v_min_100_r;
kinematic_data.VID_angle_max_100_r = VID_angle_max_100_r;
kinematic_data.VID_duration_d_max_100_r = VID_duration_d_max_100_r;
kinematic_data.VID_duration_v_max_100_r = VID_duration_v_max_100_r;
kinematic_data.VID_duration_v_min_100_r = VID_duration_v_min_100_r;
kinematic_data.VID_lick_duration_r = VID_lick_duration_r;
kinematic_data.VID_num_lick_bout_r = VID_num_lick_bout_r;
kinematic_data.VID_bout_duration_r = VID_bout_duration_r;


if (VID.DLC.FILE.session_type  == 2)
    kinematic_data.VID_d_max_100_l = VID_d_max_100_l;
    kinematic_data.VID_v_max_100_l = VID_v_max_100_l;
    kinematic_data.VID_v_min_100_l = VID_v_min_100_l;
    kinematic_data.VID_angle_max_100_l = VID_angle_max_100_l;
    
    kinematic_data.VID_duration_d_max_100_l = VID_duration_d_max_100_l;
    kinematic_data.VID_duration_v_max_100_l = VID_duration_v_max_100_l;
    kinematic_data.VID_duration_v_min_100_l = VID_duration_v_min_100_l;
    kinematic_data.VID_lick_duration_l = VID_lick_duration_l;
    
    kinematic_data.VID_num_lick_bout_l = VID_num_lick_bout_l;
    kinematic_data.VID_bout_duration_l = VID_bout_duration_l;
end
end

%% function single_dataset_signal
function signal_data = single_dataset_signal(EPHYS, BEHAVE, VID, session_type)

EPHYS_signal = EPHYS.signal.ch_data;
EPHYS_SS_ind = EPHYS.CH_sorted.SS_data.SS_ind;
EPHYS_CS_ind = EPHYS.CH_sorted.CS_data.CS_ind;
EPHYS_SS_time = EPHYS.CH_sorted.SS_data.SS_time;
EPHYS_CS_time = EPHYS.CH_sorted.CS_data.CS_time;
EPHYS_xcorr_time_100 = EPHYS.CH_EVE.align_LED.EPHYS_LED_xcorr_time_100;
EPHYS_time_30K   = EPHYS.CH_EVE_.EPHYS_time_30K;

% VID_xcorr_time_100 = EPHYS.CH_EVE.align_LED.VID_LED_xcorr_time_100;
VID_d_tip_100 = EPHYS.CH_EVE.VID_d_tip_100  ;

signal_data.EPHYS_signal = EPHYS_signal;
signal_data.EPHYS_SS_ind = EPHYS_SS_ind;
signal_data.EPHYS_CS_ind = EPHYS_CS_ind;
signal_data.EPHYS_SS_time = EPHYS_SS_time;
signal_data.EPHYS_CS_time = EPHYS_CS_time;
signal_data.EPHYS_xcorr_time_100 = EPHYS_xcorr_time_100;
signal_data.EPHYS_time_30K = EPHYS_time_30K;
signal_data.VID_d_tip_100 = VID_d_tip_100;

end

%% function ESN_correlogram
function Corr_data = ESN_correlogram(SS_time, CS_time)
bin_size_time = 1e-3; % seconds
span_window_size = (1 / bin_size_time) * (100 / 1000);
span_window_size_half = round(span_window_size / 2);
inds_span = ((-span_window_size_half+1) : 1 : (span_window_size_half))';

if (~isempty(CS_time)) && (~isempty(SS_time))
    ch_time_min = min([SS_time(1) CS_time(1)]);
    ch_time_min = max([(ch_time_min-2.0) 0]);
    ch_time_max = max([SS_time(end) CS_time(end)]) + 2.0;
    
    CH__.SS_data.SS_time =  SS_time - ch_time_min;
    CH__.CS_data.CS_time =  CS_time - ch_time_min;
    CH__.SS_data.SS_time = ESN_Round(CH__.SS_data.SS_time, bin_size_time);
    CH__.CS_data.CS_time = ESN_Round(CH__.CS_data.CS_time, bin_size_time);
    CH__.SS_data.SS_ind  = round(CH__.SS_data.SS_time .* (1/bin_size_time));
    CH__.SS_data.SS_ind( CH__.SS_data.SS_ind < 1 ) = 1;
    CH__.CS_data.CS_ind  = round(CH__.CS_data.CS_time .* (1/bin_size_time));
    CH__.CS_data.CS_ind( CH__.CS_data.CS_ind < 1 ) = 1;
    CH__.CH_data.ch_time_reconstruct = ch_time_min : bin_size_time : ch_time_max;
elseif (~isempty(CS_time))
    ch_time_min = min(  CS_time(1) );
    ch_time_min = max([(ch_time_min-2.0) 0]);
    ch_time_max = max(  CS_time(end) ) + 2.0;
    
    CH__.CS_data.CS_time =  CS_time - ch_time_min;
    CH__.CS_data.CS_time = ESN_Round(CH__.CS_data.CS_time, bin_size_time);
    CH__.CS_data.CS_ind  = round(CH__.CS_data.CS_time .* (1/bin_size_time));
    CH__.CS_data.CS_ind( CH__.CS_data.CS_ind < 1 ) = 1;
    CH__.CH_data.ch_time_reconstruct = ch_time_min : bin_size_time : ch_time_max;
elseif (~isempty(SS_time))
    ch_time_min = min( SS_time(1)  );
    ch_time_min = max([(ch_time_min-2.0) 0]);
    ch_time_max = max( SS_time(end)  ) + 2.0;
    
    CH__.SS_data.SS_time =  SS_time - ch_time_min;
    CH__.SS_data.SS_time = ESN_Round(CH__.SS_data.SS_time, bin_size_time);
    CH__.SS_data.SS_ind  = round(CH__.SS_data.SS_time .* (1/bin_size_time));
    CH__.SS_data.SS_ind( CH__.SS_data.SS_ind < 1 ) = 1;
    CH__.CH_data.ch_time_reconstruct = ch_time_min : bin_size_time : ch_time_max;
end

% SSxSS_AUTO
if (~isempty(SS_time))
    CH__.SS_data.SS_inds_reconstruct = repmat( CH__.SS_data.SS_ind(:), 1, length(inds_span)) + repmat(inds_span(:)', length(CH__.SS_data.SS_ind), 1);
    CH__.SS_data.SS_inds_reconstruct( CH__.SS_data.SS_inds_reconstruct < 1 ) = 1;
    CH__.SS_data.SS_inds_reconstruct( CH__.SS_data.SS_inds_reconstruct > length( CH__.CH_data.ch_time_reconstruct ) ) = length( CH__.CH_data.ch_time_reconstruct );
    
    CH__.SS_data.SS_event_trace = false( size(CH__.CH_data.ch_time_reconstruct) );
    CH__.SS_data.SS_event_trace( CH__.SS_data.SS_ind ) = true ;
    CH__.SS_data.SS_event_trace( 1   ) = false;
    CH__.SS_data.SS_event_trace( end ) = false;
    
    CH__.SS_data.SS_event_reconstruct = CH__.SS_data.SS_event_trace( CH__.SS_data.SS_inds_reconstruct );
    % SSxSS correlogram
    SSxSS_AUTO       = CH__.SS_data.SS_event_reconstruct;
    ss_inds_span     = repmat(inds_span(:)',     size(SS_time(:),1), 1);
    ss_bin_size_time = repmat(bin_size_time(:)', size(SS_time(:),1), 1);
else
    SSxSS_AUTO       = false(0, length(inds_span(:)'));
    ss_inds_span     = nan(0, length(inds_span(:)'));
    ss_bin_size_time = nan(0, 1);
end

% CSxSS_WITHIN
if (~isempty(CS_time)) && (~isempty(SS_time))
    CH__.CS_data.CS_inds_reconstruct = repmat( CH__.CS_data.CS_ind(:), 1, length(inds_span)) + repmat(inds_span(:)', length(CH__.CS_data.CS_ind), 1);
    CH__.CS_data.CS_inds_reconstruct( CH__.CS_data.CS_inds_reconstruct < 1 ) = 1;
    CH__.CS_data.CS_inds_reconstruct( CH__.CS_data.CS_inds_reconstruct > length( CH__.CH_data.ch_time_reconstruct ) ) = length( CH__.CH_data.ch_time_reconstruct );
    
    CH__.SS_data.SS_event_trace = false( size(CH__.CH_data.ch_time_reconstruct) );
    CH__.SS_data.SS_event_trace( CH__.SS_data.SS_ind ) = true ;
    CH__.SS_data.SS_event_trace( 1   ) = false;
    CH__.SS_data.SS_event_trace( end ) = false;
    
    CH__.SS_data.SS_event_reconstruct = CH__.SS_data.SS_event_trace( CH__.CS_data.CS_inds_reconstruct );
    % CSxSS correlogram
    CSxSS_AUTO       = CH__.SS_data.SS_event_reconstruct;
    cs_inds_span     = repmat(inds_span(:)',     size(CS_time(:),1), 1);
    cs_bin_size_time = repmat(bin_size_time(:)', size(CS_time(:),1), 1);
else
    CSxSS_AUTO       = false(0, length(inds_span(:)'));
    cs_inds_span     = nan(0, length(inds_span(:)'));
    cs_bin_size_time = nan(0, 1);
end

Corr_data = struct;
Corr_data.CS_inds_span     = cs_inds_span;
Corr_data.CS_bin_size_time = cs_bin_size_time;
Corr_data.SS_inds_span     = ss_inds_span;
Corr_data.SS_bin_size_time = ss_bin_size_time;
Corr_data.SS_SSxSS_AUTO    = SSxSS_AUTO;
Corr_data.CS_CSxSS_AUTO    = CSxSS_AUTO;
end






%% ---------- OLD CODE ---------- %%
%% Plot-9 Kinematic & SS/CS lick_onset -- FIX TO MAKE ADAPTABLE TO ANY KINEMATIC EVENT
% clearvars raster_data_
% if counter_dataset > 1
%     session_type = VID_(1).DLC.FILE.session_type;
% else
%     session_type = VID_.DLC.FILE.session_type;
% end
% for counter_dataset = 1 : 1 : num_data_set
%     EPHYS_inds_event            = EPHYS_(counter_dataset).CH_EVE.EPHYS_LED_aligned_inds_lick_onset_100;
%     VID_inds_event            = EPHYS_(counter_dataset).CH_EVE.VID_LED_aligned_inds_lick_onset_100;
%     raster_data_(counter_dataset).data = single_dataset_raster(EPHYS_(counter_dataset), EPHYS_inds_event,...
%         VID_(counter_dataset), VID_inds_event, session_type );
% end
%
% raster_data_lick_onset = concatenate_dataset(raster_data_, 'data', @vertcat);
% plot_data_lick_onset.fig_num_               = 9;
% plot_data_lick_onset.xlabel_text_raster_    = {'Lick onset (ms)'};
% plot_data_lick_onset.inds_span = EPHYS_(1).CH_EVE.inds_span_lick_onset;
% fig_handle_(plot_data_lick_onset.fig_num_)  = plot_kinematic_sscs_data(raster_data_lick_onset, plot_data_lick_onset, session_type);
%
% sgtitle(fig_handle_(plot_data_lick_onset.fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');
%% Plot-11 Frequency & SS/CS lick_onset corrs -- FIX TO MAKE ADAPTABLE TO ANY KINEMATIC EVENT
% clearvars raster_data_
% if counter_dataset > 1
%     session_type = VID_(1).DLC.FILE.session_type;
% else
%     session_type = VID_.DLC.FILE.session_type;
% end
% for counter_dataset = 1 : 1 : num_data_set
%     EPHYS_inds_event            = EPHYS_(counter_dataset).CH_EVE.EPHYS_LED_aligned_inds_lick_onset_100;
%     VID_inds_event            = EPHYS_(counter_dataset).CH_EVE.VID_LED_aligned_inds_lick_onset_100;
%     raster_data_(counter_dataset).data = single_dataset_raster(EPHYS_(counter_dataset), EPHYS_inds_event,...
%         VID_(counter_dataset), VID_inds_event, session_type );
% end
%
% raster_data_lick_onset = concatenate_dataset(raster_data_, 'data', @vertcat);
% plot_data_lick_onset.fig_num_               = 11;
% plot_data_lick_onset.xlabel_text_raster_    = {'Lick onset (ms)'};
% plot_data_lick_onset.inds_span = EPHYS_(1).CH_EVE.inds_span_lick_onset;
% fig_handle_(plot_data_lick_onset.fig_num_)  = plot_frequency_data_corr(raster_data_lick_onset, plot_data_lick_onset, session_type);
%
% sgtitle(fig_handle_(plot_data_lick_onset.fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');
%% Plot-13 Frequency & SS/CS lick_onset tip -- FIX TO MAKE ADAPTABLE TO ANY KINEMATIC EVENT
% clearvars raster_data_
% if counter_dataset > 1
%     session_type = VID_(1).DLC.FILE.session_type;
% else
%     session_type = VID_.DLC.FILE.session_type;
% end
% for counter_dataset = 1 : 1 : num_data_set
%     EPHYS_inds_event            = EPHYS_(counter_dataset).CH_EVE.EPHYS_LED_aligned_inds_lick_onset_100;
%     VID_inds_event            = EPHYS_(counter_dataset).CH_EVE.VID_LED_aligned_inds_lick_onset_100;
%     raster_data_(counter_dataset).data = single_dataset_raster(EPHYS_(counter_dataset), EPHYS_inds_event,...
%         VID_(counter_dataset), VID_inds_event, session_type );
% end
%
% raster_data_lick_onset = concatenate_dataset(raster_data_, 'data', @vertcat);
% plot_data_lick_onset.fig_num_               = 13;
% plot_data_lick_onset.xlabel_text_raster_    = {'Lick onset (ms)'};
% plot_data_lick_onset.inds_span = EPHYS_(1).CH_EVE.inds_span_lick_onset;
% fig_handle_(plot_data_lick_onset.fig_num_)  = plot_frequency_data_tip(raster_data_lick_onset, plot_data_lick_onset, session_type);
%
% sgtitle(fig_handle_(plot_data_lick_onset.fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');
%% function plot_kinematic_sscs_data
function fig_handle_ = plot_kinematic_sscs_data(raster_data, plot_data, session_type)
%%
fig_num_               = plot_data.fig_num_;
xlabel_text_raster_    = plot_data.xlabel_text_raster_;
inds_span              = plot_data.inds_span * 10;

% Grooming licks
train_data_logic_SS_grooming = raster_data.train_data_logic_SS_grooming ;
train_data_logic_CS_grooming = raster_data.train_data_logic_CS_grooming ;

% R licks
train_data_logic_SS_r = raster_data.train_data_logic_SS_r ;
train_data_logic_CS_r = raster_data.train_data_logic_CS_r ;

% L licks
if (session_type  == 2)
    train_data_logic_SS_l = raster_data.train_data_logic_SS_l ;
    train_data_logic_CS_l = raster_data.train_data_logic_CS_l ;
    
end

% d_tip
VID_d_tip_100_grooming = raster_data.VID_d_tip_100_grooming;
VID_d_tip_100_grooming_mean = nanmean(VID_d_tip_100_grooming);
VID_d_tip_100_grooming_sem = nanstd(VID_d_tip_100_grooming)/sqrt(size(VID_d_tip_100_grooming,1));

VID_d_tip_100_r = raster_data.VID_d_tip_100_r;
VID_d_tip_100_r_mean = nanmean(VID_d_tip_100_r);
VID_d_tip_100_r_sem = nanstd(VID_d_tip_100_r)/sqrt(size(VID_d_tip_100_r,1));

if (session_type  == 2)
    VID_d_tip_100_l = raster_data.VID_d_tip_100_l;
    VID_d_tip_100_l_mean = nanmean(VID_d_tip_100_l);
    VID_d_tip_100_l_sem = nanstd(VID_d_tip_100_l)/sqrt(size(VID_d_tip_100_l,1));
    
end

% v_tip
VID_v_tip_100_grooming = raster_data.VID_v_tip_100_grooming;
VID_v_tip_100_r = raster_data.VID_v_tip_100_r;
VID_v_tip_100_grooming_mean = nanmean(VID_v_tip_100_grooming);
VID_v_tip_100_grooming_sem = nanstd(VID_v_tip_100_grooming)/sqrt(size(VID_v_tip_100_grooming,1));
VID_v_tip_100_r_mean = nanmean(VID_v_tip_100_r);
VID_v_tip_100_r_sem = nanstd(VID_v_tip_100_r)/sqrt(size(VID_v_tip_100_r,1));
if (session_type  == 2)
    VID_v_tip_100_l = raster_data.VID_v_tip_100_l;
    VID_v_tip_100_l_mean = nanmean(VID_v_tip_100_l);
    VID_v_tip_100_l_sem = nanstd(VID_v_tip_100_l)/sqrt(size(VID_v_tip_100_l,1));
end

VID_angle_tip_100_grooming = abs(raster_data.VID_angle_tip_100_grooming);
VID_angle_tip_100_r = abs(raster_data.VID_angle_tip_100_r);
VID_angle_tip_100_grooming_mean = nanmean(VID_angle_tip_100_grooming);
VID_angle_tip_100_grooming_sem = nanstd(VID_angle_tip_100_grooming)/sqrt(size(VID_angle_tip_100_grooming,1));
VID_angle_tip_100_r_mean = nanmean(VID_angle_tip_100_r);
VID_angle_tip_100_r_sem = nanstd(VID_angle_tip_100_r)/sqrt(size(VID_angle_tip_100_r,1));
if (session_type  == 2)
    VID_angle_tip_100_l = abs(raster_data.VID_angle_tip_100_l);
    VID_angle_tip_100_l_mean = nanmean(VID_angle_tip_100_l);
    VID_angle_tip_100_l_sem = nanstd(VID_angle_tip_100_l)/sqrt(size(VID_angle_tip_100_l,1));
end

fig_handle_ = figure(fig_num_);
fig_handle_.WindowState = 'maximized';
clf(fig_handle_)

Line_Color = lines(7);
color_SS = Line_Color(1,:);
color_CS = Line_Color(7,:);
color_lick_onset = Line_Color(3,:);
color_lick_vmax = Line_Color(4,:);
color_lick_dmax = Line_Color(5,:);
color_lick_vmin = Line_Color(2,:);
color_lick_offset = Line_Color(4,:);

%% Grooming licks
train_data_logic_SS_ = train_data_logic_SS_grooming;
train_data_logic_CS_ = train_data_logic_CS_grooming;


subplot(3,6,2)
yyaxis left;
hold on
plot(inds_span,VID_d_tip_100_grooming_mean,'k' , 'LineWidth', 2)
plot(inds_span,VID_d_tip_100_grooming_mean+VID_d_tip_100_grooming_sem, '--k' , 'LineWidth', 1)
plot(inds_span,VID_d_tip_100_grooming_mean-VID_d_tip_100_grooming_sem, '--k' , 'LineWidth', 1)
ylabel('Lick displacement (mm)')
ylim([0 25])
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
set(gca, 'YColor', 'k')
yyaxis right;
firing_SS_ = nanmean(train_data_logic_SS_) * 100;
plot(inds_span, ESN_smooth(firing_SS_,5), 'LineWidth', 2, 'Color', 'b')
ylabel('SS Firing (spk/s)')
ylim([0 300])
set(gca, 'YColor', 'b')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
title('SS | G Licks Displacement')

subplot(3,6,5)
yyaxis left;
hold on
plot(inds_span,VID_d_tip_100_grooming_mean,'k' , 'LineWidth', 2)
plot(inds_span,VID_d_tip_100_grooming_mean+VID_d_tip_100_grooming_sem, '--k' , 'LineWidth', 1)
plot(inds_span,VID_d_tip_100_grooming_mean-VID_d_tip_100_grooming_sem, '--k' , 'LineWidth', 1)
ylabel('Lick displacement (mm)')
ylim([0 25])
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
set(gca, 'YColor', 'k')
yyaxis right;
firing_CS_ = nanmean(train_data_logic_CS_) * 100;
plot(inds_span, ESN_smooth(firing_CS_,5), 'LineWidth', 2, 'Color', 'r')
ylabel('CS Firing (spk/s)')
ylim([0 3])
set(gca, 'YColor', 'r')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
title('CS | G Licks Displacement')

subplot(3,6,8)
yyaxis left;
hold on
plot(inds_span,VID_v_tip_100_grooming_mean,'k' , 'LineWidth', 2)
plot(inds_span,VID_v_tip_100_grooming_mean+VID_v_tip_100_grooming_sem, '--k' , 'LineWidth', 1)
plot(inds_span,VID_v_tip_100_grooming_mean-VID_v_tip_100_grooming_sem, '--k' , 'LineWidth', 1)
ylabel('Lick velocity (mm/s)')
ylim([-400 400])
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
set(gca, 'YColor', 'k')
yyaxis right;
firing_SS_ = nanmean(train_data_logic_SS_) * 100;
plot(inds_span, ESN_smooth(firing_SS_,5), 'LineWidth', 2, 'Color', 'b')
ylabel('SS Firing (spk/s)')
ylim([0 300])
set(gca, 'YColor', 'b')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
title('SS | G Licks Velocity')

subplot(3,6,11)
yyaxis left;
hold on
plot(inds_span,VID_v_tip_100_grooming_mean,'k' , 'LineWidth', 2)
plot(inds_span,VID_v_tip_100_grooming_mean+VID_v_tip_100_grooming_sem, '--k' , 'LineWidth', 1)
plot(inds_span,VID_v_tip_100_grooming_mean-VID_v_tip_100_grooming_sem, '--k' , 'LineWidth', 1)
ylabel('Lick velocity (mm/s)')
ylim([-400 400])
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
set(gca, 'YColor', 'k')
yyaxis right;
firing_CS_ = nanmean(train_data_logic_CS_) * 100;
plot(inds_span, ESN_smooth(firing_CS_,5), 'LineWidth', 2, 'Color', 'r')
ylabel('CS Firing (spk/s)')
ylim([0 3])
set(gca, 'YColor', 'r')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
title('CS | G Licks Velocity')

subplot(3,6,14)
yyaxis left;
hold on
plot(inds_span,VID_angle_tip_100_grooming_mean,'k' , 'LineWidth', 2)
plot(inds_span,VID_angle_tip_100_grooming_mean+VID_angle_tip_100_grooming_sem, '--k' , 'LineWidth', 1)
plot(inds_span,VID_angle_tip_100_grooming_mean-VID_angle_tip_100_grooming_sem, '--k' , 'LineWidth', 1)
ylabel('Lick angle (deg)')
ylim([0 100])
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
set(gca, 'YColor', 'k')
yyaxis right;
firing_SS_ = nanmean(train_data_logic_SS_) * 100;
plot(inds_span, ESN_smooth(firing_SS_,5), 'LineWidth', 2, 'Color', 'b')
ylabel('SS Firing (spk/s)')
ylim([0 300])
set(gca, 'YColor', 'b')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
title('SS | G Licks Angle')

subplot(3,6,17)
yyaxis left;
hold on
plot(inds_span,VID_angle_tip_100_grooming_mean,'k' , 'LineWidth', 2)
plot(inds_span,VID_angle_tip_100_grooming_mean+VID_angle_tip_100_grooming_sem, '--k' , 'LineWidth', 1)
plot(inds_span,VID_angle_tip_100_grooming_mean-VID_angle_tip_100_grooming_sem, '--k' , 'LineWidth', 1)
ylabel('Lick angle (deg)')
ylim([0 100])
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
set(gca, 'YColor', 'k')
yyaxis right;
firing_CS_ = nanmean(train_data_logic_CS_) * 100;
plot(inds_span, ESN_smooth(firing_CS_,5), 'LineWidth', 2, 'Color', 'r')
ylabel('CS Firing (spk/s)')
ylim([0 3])
set(gca, 'YColor', 'r')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
title('CS | G Licks Angle')

%% R licks
train_data_logic_SS_ = train_data_logic_SS_r;
train_data_logic_CS_ = train_data_logic_CS_r;

subplot(3,6,3)
yyaxis left;
hold on
plot(inds_span,VID_d_tip_100_r_mean,'k' , 'LineWidth', 2)
plot(inds_span,VID_d_tip_100_r_mean+VID_d_tip_100_r_sem, '--k' , 'LineWidth', 1)
plot(inds_span,VID_d_tip_100_r_mean-VID_d_tip_100_r_sem, '--k' , 'LineWidth', 1)
ylabel('Lick displacement (mm)')
ylim([0 25])
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
set(gca, 'YColor', 'k')
yyaxis right;
firing_SS_ = nanmean(train_data_logic_SS_) * 100;
plot(inds_span, ESN_smooth(firing_SS_,5), 'LineWidth', 2, 'Color', 'b')
ylabel('SS Firing (spk/s)')
ylim([0 300])
set(gca, 'YColor', 'b')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
title('SS | R Licks Displacement')

subplot(3,6,6)
yyaxis left;
hold on
plot(inds_span,VID_d_tip_100_r_mean,'k' , 'LineWidth', 2)
plot(inds_span,VID_d_tip_100_r_mean+VID_d_tip_100_r_sem, '--k' , 'LineWidth', 1)
plot(inds_span,VID_d_tip_100_r_mean-VID_d_tip_100_r_sem, '--k' , 'LineWidth', 1)
ylabel('Lick displacement (mm)')
ylim([0 25])
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
set(gca, 'YColor', 'k')
yyaxis right;
firing_CS_ = nanmean(train_data_logic_CS_) * 100;
plot(inds_span, ESN_smooth(firing_CS_,5), 'LineWidth', 2, 'Color', 'r')
ylabel('CS Firing (spk/s)')
ylim([0 3])
set(gca, 'YColor', 'r')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
title('CS | R Licks Displacement')

subplot(3,6,9)
yyaxis left;
hold on
plot(inds_span,VID_v_tip_100_r_mean,'k' , 'LineWidth', 2)
plot(inds_span,VID_v_tip_100_r_mean+VID_v_tip_100_r_sem, '--k' , 'LineWidth', 1)
plot(inds_span,VID_v_tip_100_r_mean-VID_v_tip_100_r_sem, '--k' , 'LineWidth', 1)
ylabel('Lick velocity (mm/s)')
ylim([-400 400])
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
set(gca, 'YColor', 'k')
yyaxis right;
firing_SS_ = nanmean(train_data_logic_SS_) * 100;
plot(inds_span, ESN_smooth(firing_SS_,5), 'LineWidth', 2, 'Color', 'b')
ylabel('SS Firing (spk/s)')
ylim([0 300])
set(gca, 'YColor', 'b')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
title('SS | R Licks Velocity')

subplot(3,6,12)
yyaxis left;
hold on
plot(inds_span,VID_v_tip_100_r_mean,'k' , 'LineWidth', 2)
plot(inds_span,VID_v_tip_100_r_mean+VID_v_tip_100_r_sem, '--k' , 'LineWidth', 1)
plot(inds_span,VID_v_tip_100_r_mean-VID_v_tip_100_r_sem, '--k' , 'LineWidth', 1)
ylabel('Lick velocity (mm/s)')
ylim([-400 400])
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
set(gca, 'YColor', 'k')
yyaxis right;
firing_CS_ = nanmean(train_data_logic_CS_) * 100;
plot(inds_span, ESN_smooth(firing_CS_,5), 'LineWidth', 2, 'Color', 'r')
ylabel('CS Firing (spk/s)')
ylim([0 3])
set(gca, 'YColor', 'r')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
title('CS | R Licks Velocity')

subplot(3,6,15)
yyaxis left;
hold on
plot(inds_span,VID_angle_tip_100_r_mean,'k' , 'LineWidth', 2)
plot(inds_span,VID_angle_tip_100_r_mean+VID_angle_tip_100_r_sem, '--k' , 'LineWidth', 1)
plot(inds_span,VID_angle_tip_100_r_mean-VID_angle_tip_100_r_sem, '--k' , 'LineWidth', 1)
ylabel('Lick angle (deg)')
ylim([0 100])
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
set(gca, 'YColor', 'k')
yyaxis right;
firing_SS_ = nanmean(train_data_logic_SS_) * 100;
plot(inds_span, ESN_smooth(firing_SS_,5), 'LineWidth', 2, 'Color', 'b')
ylabel('SS Firing (spk/s)')
ylim([0 300])
set(gca, 'YColor', 'b')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
title('SS | R Licks Angle')


subplot(3,6,18)
yyaxis left;
hold on
plot(inds_span,VID_angle_tip_100_r_mean,'k' , 'LineWidth', 2)
plot(inds_span,VID_angle_tip_100_r_mean+VID_angle_tip_100_r_sem, '--k' , 'LineWidth', 1)
plot(inds_span,VID_angle_tip_100_r_mean-VID_angle_tip_100_r_sem, '--k' , 'LineWidth', 1)
ylabel('Lick angle (deg)')
ylim([0 100])
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
set(gca, 'YColor', 'k')
yyaxis right;
firing_CS_ = nanmean(train_data_logic_CS_) * 100;
plot(inds_span, ESN_smooth(firing_CS_,5), 'LineWidth', 2, 'Color', 'r')
ylabel('CS Firing (spk/s)')
ylim([0 3])
set(gca, 'YColor', 'r')
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel(xlabel_text_raster_);
title('CS | R Licks Angle')

%% L licks
if (session_type  == 2)
    train_data_logic_SS_ = train_data_logic_SS_l;
    train_data_logic_CS_ = train_data_logic_CS_l;
    
    subplot(3,6,1)
    yyaxis left;
    hold on
    plot(inds_span,VID_d_tip_100_l_mean,'k' , 'LineWidth', 2)
    plot(inds_span,VID_d_tip_100_l_mean+VID_d_tip_100_r_sem, '--k' , 'LineWidth', 1)
    plot(inds_span,VID_d_tip_100_l_mean-VID_d_tip_100_r_sem, '--k' , 'LineWidth', 1)
    ylabel('Lick displacement (mm)')
    ylim([0 25])
    if contains(xlabel_text_raster_, 'onset')
        xline(0,'LineWidth', 2, 'Color', color_lick_onset);
    elseif contains(xlabel_text_raster_, 'dmax')
        xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
    elseif contains(xlabel_text_raster_, 'offset')
        xline(0,'LineWidth', 2, 'Color', color_lick_offset);
    end
    xlim([min(inds_span)-1 max(inds_span)+1])
    xlabel(xlabel_text_raster_);
    set(gca, 'YColor', 'k')
    yyaxis right;
    firing_SS_ = nanmean(train_data_logic_SS_) * 100;
    plot(inds_span, ESN_smooth(firing_SS_,5), 'LineWidth', 2, 'Color', 'b')
    ylabel('SS Firing (spk/s)')
    ylim([0 300])
    set(gca, 'YColor', 'b')
    xlim([min(inds_span)-1 max(inds_span)+1])
    xlabel(xlabel_text_raster_);
    title('SS | L Licks Displacement')
    
    subplot(3,6,4)
    yyaxis left;
    hold on
    plot(inds_span,VID_d_tip_100_l_mean,'k' , 'LineWidth', 2)
    plot(inds_span,VID_d_tip_100_l_mean+VID_d_tip_100_r_sem, '--k' , 'LineWidth', 1)
    plot(inds_span,VID_d_tip_100_l_mean-VID_d_tip_100_r_sem, '--k' , 'LineWidth', 1)
    ylabel('Lick displacement (mm)')
    ylim([0 25])
    if contains(xlabel_text_raster_, 'onset')
        xline(0,'LineWidth', 2, 'Color', color_lick_onset);
    elseif contains(xlabel_text_raster_, 'dmax')
        xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
    elseif contains(xlabel_text_raster_, 'offset')
        xline(0,'LineWidth', 2, 'Color', color_lick_offset);
    end
    xlim([min(inds_span)-1 max(inds_span)+1])
    xlabel(xlabel_text_raster_);
    set(gca, 'YColor', 'k')
    yyaxis right;
    firing_CS_ = nanmean(train_data_logic_CS_) * 100;
    plot(inds_span, ESN_smooth(firing_CS_,5), 'LineWidth', 2, 'Color', 'r')
    ylabel('CS Firing (spk/s)')
    ylim([0 3])
    set(gca, 'YColor', 'r')
    xlim([min(inds_span)-1 max(inds_span)+1])
    xlabel(xlabel_text_raster_);
    title('CS | L Licks Displacement')
    
    subplot(3,6,7)
    yyaxis left;
    hold on
    plot(inds_span,VID_v_tip_100_l_mean,'k' , 'LineWidth', 2)
    plot(inds_span,VID_v_tip_100_l_mean+VID_v_tip_100_r_sem, '--k' , 'LineWidth', 1)
    plot(inds_span,VID_v_tip_100_l_mean-VID_v_tip_100_r_sem, '--k' , 'LineWidth', 1)
    ylabel('Lick velocity (mm/s)')
    ylim([-400 400])
    if contains(xlabel_text_raster_, 'onset')
        xline(0,'LineWidth', 2, 'Color', color_lick_onset);
    elseif contains(xlabel_text_raster_, 'dmax')
        xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
    elseif contains(xlabel_text_raster_, 'offset')
        xline(0,'LineWidth', 2, 'Color', color_lick_offset);
    end
    xlim([min(inds_span)-1 max(inds_span)+1])
    xlabel(xlabel_text_raster_);
    set(gca, 'YColor', 'k')
    yyaxis right;
    firing_SS_ = nanmean(train_data_logic_SS_) * 100;
    plot(inds_span, ESN_smooth(firing_SS_,5), 'LineWidth', 2, 'Color', 'b')
    ylabel('SS Firing (spk/s)')
    ylim([0 300])
    set(gca, 'YColor', 'b')
    xlim([min(inds_span)-1 max(inds_span)+1])
    xlabel(xlabel_text_raster_);
    title('SS | L Licks Velocity')
    
    subplot(3,6,10)
    yyaxis left;
    hold on
    plot(inds_span,VID_v_tip_100_l_mean,'k' , 'LineWidth', 2)
    plot(inds_span,VID_v_tip_100_l_mean+VID_v_tip_100_r_sem, '--k' , 'LineWidth', 1)
    plot(inds_span,VID_v_tip_100_l_mean-VID_v_tip_100_r_sem, '--k' , 'LineWidth', 1)
    ylabel('Lick velocity (mm/s)')
    ylim([-400 400])
    if contains(xlabel_text_raster_, 'onset')
        xline(0,'LineWidth', 2, 'Color', color_lick_onset);
    elseif contains(xlabel_text_raster_, 'dmax')
        xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
    elseif contains(xlabel_text_raster_, 'offset')
        xline(0,'LineWidth', 2, 'Color', color_lick_offset);
    end
    xlim([min(inds_span)-1 max(inds_span)+1])
    xlabel(xlabel_text_raster_);
    set(gca, 'YColor', 'k')
    yyaxis right;
    firing_CS_ = nanmean(train_data_logic_CS_) * 100;
    plot(inds_span, ESN_smooth(firing_CS_,5), 'LineWidth', 2, 'Color', 'r')
    ylabel('CS Firing (spk/s)')
    ylim([0 3])
    set(gca, 'YColor', 'r')
    xlim([min(inds_span)-1 max(inds_span)+1])
    xlabel(xlabel_text_raster_);
    title('CS | L Licks Velocity')
    
    subplot(3,6,13)
    yyaxis left;
    hold on
    plot(inds_span,VID_angle_tip_100_l_mean,'k' , 'LineWidth', 2)
    plot(inds_span,VID_angle_tip_100_l_mean+VID_angle_tip_100_r_sem, '--k' , 'LineWidth', 1)
    plot(inds_span,VID_angle_tip_100_l_mean-VID_angle_tip_100_r_sem, '--k' , 'LineWidth', 1)
    ylabel('Lick angle (deg)')
    ylim([0 100])
    if contains(xlabel_text_raster_, 'onset')
        xline(0,'LineWidth', 2, 'Color', color_lick_onset);
    elseif contains(xlabel_text_raster_, 'dmax')
        xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
    elseif contains(xlabel_text_raster_, 'offset')
        xline(0,'LineWidth', 2, 'Color', color_lick_offset);
    end
    xlim([min(inds_span)-1 max(inds_span)+1])
    xlabel(xlabel_text_raster_);
    set(gca, 'YColor', 'k')
    yyaxis right;
    firing_SS_ = nanmean(train_data_logic_SS_) * 100;
    plot(inds_span, ESN_smooth(firing_SS_,5), 'LineWidth', 2, 'Color', 'b')
    ylabel('SS Firing (spk/s)')
    ylim([0 300])
    set(gca, 'YColor', 'b')
    xlim([min(inds_span)-1 max(inds_span)+1])
    xlabel(xlabel_text_raster_);
    title('SS | L Licks Angle')
    
    subplot(3,6,16)
    yyaxis left;
    hold on
    plot(inds_span,VID_angle_tip_100_l_mean,'k' , 'LineWidth', 2)
    plot(inds_span,VID_angle_tip_100_l_mean+VID_angle_tip_100_r_sem, '--k' , 'LineWidth', 1)
    plot(inds_span,VID_angle_tip_100_l_mean-VID_angle_tip_100_r_sem, '--k' , 'LineWidth', 1)
    ylabel('Lick angle (deg)')
    ylim([0 100])
    if contains(xlabel_text_raster_, 'onset')
        xline(0,'LineWidth', 2, 'Color', color_lick_onset);
    elseif contains(xlabel_text_raster_, 'dmax')
        xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
    elseif contains(xlabel_text_raster_, 'offset')
        xline(0,'LineWidth', 2, 'Color', color_lick_offset);
    end
    xlim([min(inds_span)-1 max(inds_span)+1])
    xlabel(xlabel_text_raster_);
    set(gca, 'YColor', 'k')
    yyaxis right;
    firing_CS_ = nanmean(train_data_logic_CS_) * 100;
    plot(inds_span, ESN_smooth(firing_CS_,5), 'LineWidth', 2, 'Color', 'r')
    ylabel('CS Firing (spk/s)')
    ylim([0 3])
    set(gca, 'YColor', 'r')
    xlim([min(inds_span)-1 max(inds_span)+1])
    xlabel(xlabel_text_raster_);
    title('CS | L Licks Angle')
    
    
end


end
%% function plot_frequency_data_corr
function fig_handle_ = plot_frequency_data_corr(raster_data, plot_data, session_type)
%%
fig_num_               = plot_data.fig_num_;
xlabel_text_raster_    = plot_data.xlabel_text_raster_;
xlabel_text_raster_bout_    = plot_data.xlabel_text_raster_bout_;

inds_span              = plot_data.inds_span * 10;

Fs = 100;

%bout
inds_span_LB           =(inds_span(1)-700) :10: inds_span(1)-10;
inds_span_UB           = (inds_span(end)+10) :10: (inds_span(end)+700);
inds_span_bout      = [inds_span_LB inds_span inds_span_UB]/1000;

% %corr: [0 1]s
% inds_span_UB           = (inds_span(end)+10) :10: (inds_span(end)+700);
% inds_span_corr = [inds_span(31:60) inds_span_UB]/1000;

%corr: [-1 1]s
inds_span_LB           =(inds_span(1)-700) :10: inds_span(1)-10;
inds_span_UB           = (inds_span(end)+10) :10: (inds_span(end)+700);
inds_span_corr      = [inds_span_LB inds_span inds_span_UB]/1000;

% All licks
train_data_logic_SS_all_corr = raster_data.train_data_logic_SS_all_corr;
train_data_logic_CS_all_corr = raster_data.train_data_logic_CS_all_corr;
train_data_logic_lick_onset_all_corr = raster_data.train_data_logic_lick_onset_all_corr;
train_data_logic_lick_dmax_all_corr = raster_data.train_data_logic_lick_dmax_all_corr;
train_data_logic_lick_offset_all_corr = raster_data.train_data_logic_lick_offset_all_corr;
train_data_logic_SS_all = raster_data.train_data_logic_SS_all ;
train_data_logic_CS_all  = raster_data.train_data_logic_CS_all ;
train_data_logic_lick_onset_all     = raster_data.train_data_logic_lick_onset_all ;
train_data_logic_lick_dmax_all    = raster_data.train_data_logic_lick_dmax_all ;
train_data_logic_lick_offset_all     = raster_data.train_data_logic_lick_offset_all ;
VID_d_tip_100_all = raster_data.VID_d_tip_100_all ;
VID_v_tip_100_all = raster_data.VID_v_tip_100_all ;
VID_angle_tip_100_all = raster_data.VID_angle_tip_100_all;


% Grooming licks
train_data_logic_SS_grooming_corr = raster_data.train_data_logic_SS_grooming_corr;
train_data_logic_CS_grooming_corr = raster_data.train_data_logic_CS_grooming_corr;
train_data_logic_lick_onset_grooming_corr = raster_data.train_data_logic_lick_onset_grooming_corr;
train_data_logic_lick_dmax_grooming_corr = raster_data.train_data_logic_lick_dmax_grooming_corr;
train_data_logic_lick_offset_grooming_corr = raster_data.train_data_logic_lick_offset_grooming_corr;

% R licks
train_data_logic_SS_r_corr = raster_data.train_data_logic_SS_r_corr;
train_data_logic_CS_r_corr = raster_data.train_data_logic_CS_r_corr;
train_data_logic_lick_onset_r_corr = raster_data.train_data_logic_lick_onset_r_corr;
train_data_logic_lick_dmax_r_corr = raster_data.train_data_logic_lick_dmax_r_corr;
train_data_logic_lick_offset_r_corr = raster_data.train_data_logic_lick_offset_r_corr;

% L licks
if (session_type  == 2)
    train_data_logic_SS_l_corr = raster_data.train_data_logic_SS_l_corr;
    train_data_logic_CS_l_corr = raster_data.train_data_logic_CS_l_corr;
    train_data_logic_lick_onset_l_corr = raster_data.train_data_logic_lick_onset_l_corr;
    train_data_logic_lick_dmax_l_corr = raster_data.train_data_logic_lick_dmax_l_corr;
    train_data_logic_lick_offset_l_corr = raster_data.train_data_logic_lick_offset_l_corr;
end

Line_Color = lines(7);
color_SS = Line_Color(1,:);
color_CS = Line_Color(7,:);
color_lick_onset = Line_Color(3,:);
color_lick_vmax = Line_Color(4,:);
color_lick_dmax = Line_Color(5,:);
color_lick_vmin = Line_Color(2,:);
color_lick_offset = Line_Color(4,:);


fig_handle_ = figure(fig_num_);
fig_handle_.WindowState = 'maximized';
clf(fig_handle_)

%% Grooming licks
train_data_logic_SS_corr_ = train_data_logic_SS_grooming_corr;
train_data_logic_CS_corr_ = train_data_logic_CS_grooming_corr;
train_data_logic_lick_onset_corr_ = train_data_logic_lick_onset_grooming_corr ;
train_data_logic_lick_dmax_corr_ = train_data_logic_lick_dmax_grooming_corr ;
train_data_logic_lick_offset_corr_ = train_data_logic_lick_offset_grooming_corr ;

% auto corr
if contains(xlabel_text_raster_, 'onset')
    lick_auto_corr = nansum(train_data_logic_lick_onset_corr_)/nansum(nansum(train_data_logic_lick_onset_corr_));
elseif contains(xlabel_text_raster_, 'dmax')
    lick_auto_corr = nansum(train_data_logic_lick_dmax_corr_)/nansum(nansum(train_data_logic_lick_dmax_corr_));
elseif contains(xlabel_text_raster_, 'offset')
    lick_auto_corr = nansum(train_data_logic_lick_offset_corr_)/nansum(nansum(train_data_logic_lick_offset_corr_));
end
lick_auto_corr(100) = 0;
if (length(lick_auto_corr) < length(inds_span_corr))
    lick_auto_corr = zeros(1,length(inds_span_corr));
end


% cross corr
lick_cross_corr_SS = nansum(train_data_logic_SS_corr_)/nansum(nansum(train_data_logic_SS_corr_));
lick_cross_corr_CS = nansum(train_data_logic_CS_corr_)/nansum(nansum(train_data_logic_CS_corr_));

% fs = 100;
% y = fft(lick_auto_corr);
% n = length(lick_auto_corr);
% f = (0:n-1)*(fs/n);
% power = abs(y).^2/n;
%
% plot(f,power);
% xlabel('Freq')
% ylabel('Power')
% xlim([0 10])

subplot(6,5,4)
% plot(inds_span_corr, ESN_smooth(lick_auto_corr,5), 'k','LineWidth', 2)
area(inds_span_corr,lick_auto_corr, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 1.5)
title('Prob lick | G licks')
xlabel(xlabel_text_raster_bout_);
xlim([-1 1])
ylabel('Probability')
if (isnan(max(lick_auto_corr)))
    ylim([0 0.1])
else
    ylim([0 (max(lick_auto_corr) + max(lick_auto_corr)/10)])
end
set(gca, 'YColor', 'k')
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end

subplot(6,5,9)
% periodogram(lick_auto_corr,rectwin(length(lick_auto_corr)), length(lick_auto_corr), Fs)
N = length(lick_auto_corr);
xdft = fft(lick_auto_corr);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(lick_auto_corr):Fs/2;
plot(freq,10*log10(psdx), 'Color', [0.7 0.7 0.7], 'LineWidth', 2)
title('Periodogram | G Licks')
xlabel('Freq (Hz)')
xlim([0 10])
ylabel('Pow/Freq (dB/Hz)')


subplot(6,5,14)
% plot(inds_span_corr, ESN_smooth(lick_cross_corr_SS,5), 'k','LineWidth', 2)
area(inds_span_corr,lick_cross_corr_SS, 'FaceColor', 'k', 'LineWidth', 1.5)
title('Prob SS | G licks')
xlabel(xlabel_text_raster_bout_);
xlim([-1 1])
ylabel('Probability')
if (isnan(max(lick_cross_corr_SS)))
    ylim([0 0.1])
else
    ylim([0 (max(lick_cross_corr_SS) + max(lick_cross_corr_SS)/10)])
end
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end

subplot(6,5,19)
% periodogram(lick_cross_corr_SS,rectwin(length(lick_cross_corr_SS)), length(lick_cross_corr_SS), Fs)
N = length(lick_cross_corr_SS);
xdft = fft(lick_cross_corr_SS);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(lick_cross_corr_SS):Fs/2;
plot(freq,10*log10(psdx), 'k', 'LineWidth', 2)
title('Periodogram | G Licks - SS')
xlabel('Freq (Hz)')
xlim([0 10])
ylabel('Pow/Freq (dB/Hz)')

subplot(6,5,24)
% plot(inds_span_corr, ESN_smooth(lick_cross_corr_CS,5), 'k','LineWidth', 2)
area(inds_span_corr,ESN_smooth(lick_cross_corr_CS,5), 'FaceColor', 'k', 'LineWidth', 1.5)
title('Prob CS | G licks')
xlabel(xlabel_text_raster_bout_);
xlim([-1 1])
ylabel('Probability')
if (isnan(max(lick_cross_corr_CS)))
    ylim([0 0.1])
else
    ylim([0 (max(lick_cross_corr_CS) + max(lick_cross_corr_CS)/10)])
end
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end

subplot(6,5,29)
% periodogram(lick_cross_corr_CS,rectwin(length(lick_cross_corr_CS)), length(lick_cross_corr_CS), Fs)
N = length(lick_cross_corr_CS);
xdft = fft(lick_cross_corr_CS);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(lick_cross_corr_CS):Fs/2;
plot(freq,10*log10(psdx),'k', 'LineWidth', 2)
title('Periodogram | G Licks - CS')
xlabel('Freq (Hz)')
xlim([0 10])
ylabel('Pow/Freq (dB/Hz)')

%% R Licks
train_data_logic_SS_corr_ = train_data_logic_SS_r_corr;
train_data_logic_CS_corr_ = train_data_logic_CS_r_corr;
train_data_logic_lick_onset_corr_ = train_data_logic_lick_onset_r_corr ;
train_data_logic_lick_dmax_corr_ = train_data_logic_lick_dmax_r_corr ;
train_data_logic_lick_offset_corr_ = train_data_logic_lick_offset_r_corr ;
% auto corr
if contains(xlabel_text_raster_, 'onset')
    lick_auto_corr = nansum(train_data_logic_lick_onset_corr_)/nansum(nansum(train_data_logic_lick_onset_corr_));
elseif contains(xlabel_text_raster_, 'dmax')
    lick_auto_corr = nansum(train_data_logic_lick_dmax_corr_)/nansum(nansum(train_data_logic_lick_dmax_corr_));
elseif contains(xlabel_text_raster_, 'offset')
    lick_auto_corr = nansum(train_data_logic_lick_offset_corr_)/nansum(nansum(train_data_logic_lick_offset_corr_));
end
lick_auto_corr(100) = 0;
if (length(lick_auto_corr) < length(inds_span_corr))
    lick_auto_corr = zeros(1,length(inds_span_corr));
end

% cross corr
lick_cross_corr_SS = nansum(train_data_logic_SS_corr_)/nansum(nansum(train_data_logic_SS_corr_));
lick_cross_corr_CS = nansum(train_data_logic_CS_corr_)/nansum(nansum(train_data_logic_CS_corr_));

subplot(6,5,5)
% plot(inds_span_corr, ESN_smooth(lick_auto_corr,5), 'k','LineWidth', 1)
area(inds_span_corr,lick_auto_corr, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 1.5)
title('Prob lick | R licks')
xlabel(xlabel_text_raster_bout_);
xlim([-1 1])
ylabel('Probability')
ylim([0 (max(lick_auto_corr) + max(lick_auto_corr)/10)])
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end

subplot(6,5,10)
% periodogram(lick_auto_corr,rectwin(length(lick_auto_corr)), length(lick_auto_corr), Fs)
N = length(lick_auto_corr);
xdft = fft(lick_auto_corr);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(lick_auto_corr):Fs/2;
plot(freq,10*log10(psdx), 'Color', [0.7 0.7 0.7], 'LineWidth', 2)
title('Periodogram | R Licks')
xlabel('Freq (Hz)')
xlim([0 10])
ylabel('Pow/Freq (dB/Hz)')

subplot(6,5,15)
% plot(inds_span_corr, ESN_smooth(lick_cross_corr_SS,5), 'k','LineWidth', 1)
area(inds_span_corr,lick_cross_corr_SS, 'FaceColor', 'k', 'LineWidth', 1.5)
title('Prob SS | R licks')
xlabel(xlabel_text_raster_bout_);
xlim([-1 1])
ylabel('Probability')
if (isnan(max(lick_cross_corr_SS)))
    ylim([0 0.1])
else
    ylim([0 (max(lick_cross_corr_SS) + max(lick_cross_corr_SS)/10)])
end
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end

subplot(6,5,20)
% periodogram(lick_cross_corr_SS,rectwin(length(lick_cross_corr_SS)), length(lick_cross_corr_SS), Fs)
N = length(lick_cross_corr_SS);
xdft = fft(lick_cross_corr_SS);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(lick_cross_corr_SS):Fs/2;
plot(freq,10*log10(psdx), 'k', 'LineWidth', 2)
title('Periodogram | R Licks - SS')
xlabel('Freq (Hz)')
xlim([0 10])
ylabel('Pow/Freq (dB/Hz)')

subplot(6,5,25)
% plot(inds_span_corr, ESN_smooth(lick_cross_corr_CS,5), 'k','LineWidth', 1)
area(inds_span_corr,ESN_smooth(lick_cross_corr_CS,5), 'FaceColor','k', 'LineWidth', 1.5)
title('Prob CS | R licks')
xlabel(xlabel_text_raster_bout_);
xlim([-1 1])
ylabel('Probability')
if (isnan(max(lick_cross_corr_CS)))
    ylim([0 0.1])
else
    ylim([0 (max(lick_cross_corr_CS) + max(lick_cross_corr_CS)/10)])
end
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end

subplot(6,5,30)
% periodogram(lick_cross_corr_CS,rectwin(length(lick_cross_corr_CS)), length(lick_cross_corr_CS), Fs)
N = length(lick_cross_corr_CS);
xdft = fft(lick_cross_corr_CS);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(lick_cross_corr_CS):Fs/2;
plot(freq,10*log10(psdx),'k', 'LineWidth', 2)
title('Periodogram | R Licks - CS')
xlabel('Freq (Hz)')
xlim([0 10])
ylabel('Pow/Freq (dB/Hz)')

%% L Licks
if (session_type  == 2)
    train_data_logic_SS_corr_ = train_data_logic_SS_l_corr;
    train_data_logic_CS_corr_ = train_data_logic_CS_l_corr;
    train_data_logic_lick_onset_corr_ = train_data_logic_lick_onset_l_corr ;
    train_data_logic_lick_dmax_corr_ = train_data_logic_lick_dmax_l_corr ;
    train_data_logic_lick_offset_corr_ = train_data_logic_lick_offset_l_corr ;
    % auto corr
    if contains(xlabel_text_raster_, 'onset')
        lick_auto_corr = nansum(train_data_logic_lick_onset_corr_)/nansum(nansum(train_data_logic_lick_onset_corr_));
    elseif contains(xlabel_text_raster_, 'dmax')
        lick_auto_corr = nansum(train_data_logic_lick_dmax_corr_)/nansum(nansum(train_data_logic_lick_dmax_corr_));
    elseif contains(xlabel_text_raster_, 'offset')
        lick_auto_corr = nansum(train_data_logic_lick_offset_corr_)/nansum(nansum(train_data_logic_lick_offset_corr_));
    end
    lick_auto_corr(100) = 0;
    if (length(lick_auto_corr) < length(inds_span_corr))
        lick_auto_corr = zeros(1,length(inds_span_corr));
    end
    
    % cross corr
    lick_cross_corr_SS = nansum(train_data_logic_SS_corr_)/nansum(nansum(train_data_logic_SS_corr_));
    lick_cross_corr_CS = nansum(train_data_logic_CS_corr_)/nansum(nansum(train_data_logic_CS_corr_));
    
    subplot(6,5,3)
    % plot(inds_span_corr, ESN_smooth(lick_auto_corr,5), 'k','LineWidth', 1)
    area(inds_span_corr,lick_auto_corr, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 1.5)
    title('Prob Lick | L licks')
    xlabel(xlabel_text_raster_bout_);
    xlim([-1 1])
    ylabel('Probability')
    if (isnan(max(lick_auto_corr)) || max(lick_auto_corr) == 0)
        ylim([0 0.1])
    else
        ylim([0 (max(lick_auto_corr) + max(lick_auto_corr)/10)])
    end
    if contains(xlabel_text_raster_, 'onset')
        xline(0,'LineWidth', 2, 'Color', color_lick_onset);
    elseif contains(xlabel_text_raster_, 'dmax')
        xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
    elseif contains(xlabel_text_raster_, 'offset')
        xline(0,'LineWidth', 2, 'Color', color_lick_offset);
    end
    
    subplot(6,5,8)
    % periodogram(lick_auto_corr,rectwin(length(lick_auto_corr)), length(lick_auto_corr), Fs)
    N = length(lick_auto_corr);
    xdft = fft(lick_auto_corr);
    xdft = xdft(1:N/2+1);
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:Fs/length(lick_auto_corr):Fs/2;
    plot(freq,10*log10(psdx), 'Color', [0.7 0.7 0.7], 'LineWidth', 2)
    title('Periodogram | L Licks')
    xlabel('Freq (Hz)')
    xlim([0 10])
    ylabel('Pow/Freq (dB/Hz)')
    
    
    subplot(6,5,13)
    % plot(inds_span_corr, ESN_smooth(lick_cross_corr_SS,5), 'k','LineWidth', 1)
    area(inds_span_corr,lick_cross_corr_SS, 'FaceColor','k', 'LineWidth', 1.5)
    title('Prob SS | L licks')
    xlabel(xlabel_text_raster_bout_);
    xlim([-1 1])
    ylabel('Probability')
    if (isnan(max(lick_cross_corr_SS)))
        ylim([0 0.1])
    else
        ylim([0 (max(lick_cross_corr_SS) + max(lick_cross_corr_SS)/10)])
    end
    if contains(xlabel_text_raster_, 'onset')
        xline(0,'LineWidth', 2, 'Color', color_lick_onset);
    elseif contains(xlabel_text_raster_, 'dmax')
        xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
    elseif contains(xlabel_text_raster_, 'offset')
        xline(0,'LineWidth', 2, 'Color', color_lick_offset);
    end
    
    subplot(6,5,18)
    % periodogram(lick_cross_corr_SS,rectwin(length(lick_cross_corr_SS)), length(lick_cross_corr_SS), Fs)
    N = length(lick_cross_corr_SS);
    xdft = fft(lick_cross_corr_SS);
    xdft = xdft(1:N/2+1);
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:Fs/length(lick_cross_corr_SS):Fs/2;
    plot(freq,10*log10(psdx),'k', 'LineWidth', 2)
    title('Periodogram | L Licks - SS')
    xlabel('Freq (Hz)')
    xlim([0 10])
    ylabel('Pow/Freq (dB/Hz)')
    
    subplot(6,5,23)
    % plot(inds_span_corr, ESN_smooth(lick_cross_corr_CS,5), 'k','LineWidth', 1)
    area(inds_span_corr,ESN_smooth(lick_cross_corr_CS,5), 'FaceColor', 'k', 'LineWidth', 1.5)
    title('Prob CS | L licks')
    xlabel(xlabel_text_raster_bout_);
    xlim([-1 1])
    ylabel('Probability')
    if (isnan(max(lick_cross_corr_CS)))
        ylim([0 0.1])
    else
        ylim([0 (max(lick_cross_corr_CS) + max(lick_cross_corr_CS)/10)])
    end
    if contains(xlabel_text_raster_, 'onset')
        xline(0,'LineWidth', 2, 'Color', color_lick_onset);
    elseif contains(xlabel_text_raster_, 'dmax')
        xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
    elseif contains(xlabel_text_raster_, 'offset')
        xline(0,'LineWidth', 2, 'Color', color_lick_offset);
    end
    
    subplot(6,5,28)
    % periodogram(lick_cross_corr_CS,rectwin(length(lick_cross_corr_CS)), length(lick_cross_corr_CS), Fs)
    N = length(lick_cross_corr_CS);
    xdft = fft(lick_cross_corr_CS);
    xdft = xdft(1:N/2+1);
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:Fs/length(lick_cross_corr_CS):Fs/2;
    plot(freq,10*log10(psdx), 'k' , 'LineWidth', 2)
    title('Periodogram | L Licks - CS')
    xlabel('Freq (Hz)')
    xlim([0 10])
    ylabel('Pow/Freq (dB/Hz)')
    
    
end

%% All licks phase relationships
train_data_logic_SS_corr_ = train_data_logic_SS_all_corr;
train_data_logic_CS_corr_ = train_data_logic_CS_all_corr;
train_data_logic_lick_onset_corr_ = train_data_logic_lick_onset_all_corr ;
train_data_logic_lick_dmax_corr_ = train_data_logic_lick_dmax_all_corr ;
train_data_logic_lick_offset_corr_ = train_data_logic_lick_offset_all_corr ;

% auto corr
if contains(xlabel_text_raster_, 'onset')
    lick_auto_corr = nansum(train_data_logic_lick_onset_corr_)/nansum(nansum(train_data_logic_lick_onset_corr_));
elseif contains(xlabel_text_raster_, 'dmax')
    lick_auto_corr = nansum(train_data_logic_lick_dmax_corr_)/nansum(nansum(train_data_logic_lick_dmax_corr_));
elseif contains(xlabel_text_raster_, 'offset')
    lick_auto_corr = nansum(train_data_logic_lick_offset_corr_)/nansum(nansum(train_data_logic_lick_offset_corr_));
end
lick_auto_corr(100) = 0;
if (length(lick_auto_corr) < length(inds_span_corr))
    lick_auto_corr = zeros(1,length(inds_span_corr));
end

% cross corr
lick_cross_corr_SS = nansum(train_data_logic_SS_corr_)/nansum(nansum(train_data_logic_SS_corr_));
lick_cross_corr_CS = nansum(train_data_logic_CS_corr_)/nansum(nansum(train_data_logic_CS_corr_));

subplot(6,5,[1 2 6 7])
hold on
yyaxis left;
area(inds_span_corr,ESN_smooth(lick_cross_corr_SS,5), 'FaceColor', 'k', 'LineWidth', 1.5)
ylabel('Prob. SS')
if (isnan(max(lick_cross_corr_SS)))
    ylim([0 0.1])
else
    ylim([0 (max(lick_cross_corr_SS) + max(lick_cross_corr_SS)/10)])
end
set(gca, 'YColor', 'k')
yyaxis right;
area(inds_span_corr,ESN_smooth(lick_auto_corr,5), 'FaceColor', [0.7 0.7 0.7], 'LineWidth', 1.5)
% plot(inds_span_corr, ESN_smooth(lick_auto_corr,5), 'Color' , [0.7 0.7 0.7], 'LineWIdth', 2)
ylabel('Prob. lick onset')
xlabel(xlabel_text_raster_')
title('All licks')
set(gca, 'YColor', [0.7 0.7 0.7])
xline(0,'--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2 );
if (isnan(max(lick_auto_corr)))
    ylim([0 0.1])
else
    ylim([0 (max(lick_auto_corr) + max(lick_auto_corr)/10)])
end

subplot(6,5,[11 12 16 17])
N_lick_auto_corr = length(lick_auto_corr);
xdft_lick_auto_corr = fft(lick_auto_corr);
xdft_lick_auto_corr = xdft_lick_auto_corr(1:N_lick_auto_corr/2+1);
psdx_lick_auto_corr = (1/(Fs*N_lick_auto_corr)) * abs(xdft_lick_auto_corr).^2;
psdx_lick_auto_corr(2:end-1) = 2*psdx_lick_auto_corr(2:end-1);
N_lick_cross_corr_SS = length(lick_cross_corr_SS);
xdft_lick_cross_corr_SS = fft(lick_cross_corr_SS);
xdft_lick_cross_corr_SS = xdft_lick_cross_corr_SS(1:N_lick_cross_corr_SS/2+1);
psdx_lick_cross_corr_SS = (1/(Fs*N_lick_cross_corr_SS)) * abs(xdft_lick_cross_corr_SS).^2;
psdx_lick_cross_corr_SS(2:end-1) = 2*psdx_lick_cross_corr_SS(2:end-1);
freq_lick_auto_corr = 0:Fs/length(lick_auto_corr):Fs/2;
freq_lick_cross_corr_SS = 0:Fs/length(lick_cross_corr_SS):Fs/2;
hold on;
plot(freq_lick_auto_corr,10*log10(psdx_lick_auto_corr), 'Color', [0.7 0.7 0.7], 'LineWidth', 2)
plot(freq_lick_cross_corr_SS,10*log10(psdx_lick_cross_corr_SS), 'k', 'LineWidth', 2)
xlabel('Freq (Hz)')
xlim([0 10])
ylabel('Pow/Freq (dB/Hz)')

subplot(6,5,[21 22 26 27])
x_axis = categorical({'Phase diff (deg)'});
x_axis = reordercats(x_axis, {'Phase diff (deg)'});
phase_diff = rad2deg(phdiffmeasure(lick_auto_corr,lick_cross_corr_SS));
yyaxis left;
bar(x_axis(1),phase_diff, 'k' )
ylabel('Phase diff (deg)')
if(phase_diff > 0)
    ylim([0 (phase_diff+40)])
else
    ylim([(phase_diff-40) 0])
    
end
set(gca, 'YColor', 'k')
yyaxis right;
set(gca, 'YColor', 'k')



end
%% function plot_frequency_data_tip
function fig_handle_ = plot_frequency_data_tip(raster_data, plot_data, session_type)
%%
fig_num_               = plot_data.fig_num_;
xlabel_text_raster_    = plot_data.xlabel_text_raster_;
xlabel_text_raster_bout_    = plot_data.xlabel_text_raster_bout_;

inds_span              = plot_data.inds_span * 10;

Fs = 100;

% All licks
train_data_logic_SS_all = raster_data.train_data_logic_SS_all ;
train_data_logic_CS_all  = raster_data.train_data_logic_CS_all ;
VID_d_tip_100_all = raster_data.VID_d_tip_100_all ;
VID_v_tip_100_all = raster_data.VID_v_tip_100_all ;
VID_angle_tip_100_all = raster_data.VID_angle_tip_100_all;


% Grooming licks
train_data_logic_SS_grooming = raster_data.train_data_logic_SS_grooming ;
train_data_logic_CS_grooming  = raster_data.train_data_logic_CS_grooming ;
VID_d_tip_100_grooming = raster_data.VID_d_tip_100_grooming ;


% R licks
train_data_logic_SS_r = raster_data.train_data_logic_SS_r ;
train_data_logic_CS_r  = raster_data.train_data_logic_CS_r ;
VID_d_tip_100_r = raster_data.VID_d_tip_100_r ;

% L licks
if (session_type  == 2)
    train_data_logic_SS_l = raster_data.train_data_logic_SS_l ;
    train_data_logic_CS_l  = raster_data.train_data_logic_CS_l ;
    VID_d_tip_100_l = raster_data.VID_d_tip_100_l ;
end

Line_Color = lines(7);
color_SS = Line_Color(1,:);
color_CS = Line_Color(7,:);
color_lick_onset = Line_Color(3,:);
color_lick_vmax = Line_Color(4,:);
color_lick_dmax = Line_Color(5,:);
color_lick_vmin = Line_Color(2,:);
color_lick_offset = Line_Color(4,:);


fig_handle_ = figure(fig_num_);
fig_handle_.WindowState = 'maximized';
clf(fig_handle_)

%% Grooming licks
train_data_logic_SS_ = train_data_logic_SS_grooming;
train_data_logic_CS_ = train_data_logic_CS_grooming;
d_tip_ = VID_d_tip_100_grooming;

mean_d_tip_ = nanmean(d_tip_);
se_d_tip_ = nanstd(d_tip_)/sqrt(length(d_tip_));
mean_SS_ = ESN_smooth(nanmean(train_data_logic_SS_)*100,5);
se_SS_ = ESN_smooth(nanstd(train_data_logic_SS_)/sqrt(length(train_data_logic_SS_))*100,5);
mean_CS_ = ESN_smooth(nanmean(train_data_logic_CS_)*100,10);
se_CS_ = ESN_smooth(nanstd(train_data_logic_CS_)/sqrt(length(train_data_logic_CS_))*100,10);


subplot(6,5,4)
hold on
plot(inds_span, mean_d_tip_, 'Color', [0.7 0.7 0.7],'LineWidth', 2)
plot(inds_span, mean_d_tip_ + se_d_tip_, 'Color', [0.7 0.7 0.7],'LineWidth',1)
plot(inds_span, mean_d_tip_ - se_d_tip_, 'Color', [0.7 0.7 0.7],'LineWidth', 1)
title('Displacement | G licks')
xlabel(xlabel_text_raster_);
xlim([-300 300])
ylabel('Displacement (mm)')
if (isnan(max(mean_d_tip_)))
    ylim([0 10])
else
    ylim([0 (max(mean_d_tip_) + 1)])
end
set(gca, 'YColor', 'k')
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end

subplot(6,5,9)
N = length(mean_d_tip_);
xdft = fft(mean_d_tip_);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(mean_d_tip_):Fs/2;
plot(freq,10*log10(psdx), 'Color', [0.7 0.7 0.7], 'LineWidth', 2)
title('Periodogram | G Licks')
xlabel('Freq (Hz)')
xlim([0 10])
ylabel('Pow/Freq (dB/Hz)')


subplot(6,5,14)
hold on
plot(inds_span, mean_SS_, 'Color', 'k','LineWidth', 2)
plot(inds_span, mean_SS_ + se_SS_, 'Color', 'k' ,'LineWidth',1)
plot(inds_span, mean_SS_ - se_SS_, 'Color', 'k' ,'LineWidth', 1)
title('SS firing | G licks')
xlabel(xlabel_text_raster_);
xlim([-300 300])
ylabel('SS firing (spks/s)')
if (isnan(max(mean_SS_)))
    ylim([0 300])
else
    ylim([0 (max(mean_SS_) + max(mean_SS_)/2)])
end
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end

subplot(6,5,19)
N = length(mean_SS_);
xdft = fft(mean_SS_);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(mean_SS_):Fs/2;
plot(freq,10*log10(psdx), 'k', 'LineWidth', 2)
title('Periodogram | G Licks - SS')
xlabel('Freq (Hz)')
xlim([0 10])
ylabel('Pow/Freq (dB/Hz)')

subplot(6,5,24)
hold on
plot(inds_span, mean_CS_, 'Color', 'k','LineWidth', 2)
plot(inds_span, mean_CS_ + se_CS_, 'Color', 'k' ,'LineWidth',1)
plot(inds_span, mean_CS_ - se_CS_, 'Color', 'k' ,'LineWidth', 1)
title('CS firing | G licks')
xlabel(xlabel_text_raster_);
xlim([-300 300])
ylabel('CS firing (spks/s)')
if (isnan(max(mean_SS_)))
    ylim([0 3])
else
    ylim([0 (max(mean_CS_) + max(mean_CS_)/2)])
end

if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end

subplot(6,5,29)
N = length(mean_CS_);
xdft = fft(mean_CS_);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(mean_CS_):Fs/2;
plot(freq,10*log10(psdx), 'k', 'LineWidth', 2)
title('Periodogram | G Licks - CS')
xlabel('Freq (Hz)')
xlim([0 10])
ylabel('Pow/Freq (dB/Hz)')

%% R Licks
train_data_logic_SS_ = train_data_logic_SS_r;
train_data_logic_CS_ = train_data_logic_CS_r;
d_tip_ = VID_d_tip_100_r;

mean_d_tip_ = nanmean(d_tip_);
se_d_tip_ = nanstd(d_tip_)/sqrt(length(d_tip_));
mean_SS_ = ESN_smooth(nanmean(train_data_logic_SS_)*100,5);
se_SS_ = ESN_smooth(nanstd(train_data_logic_SS_)/sqrt(length(train_data_logic_SS_))*100,5);
mean_CS_ = ESN_smooth(nanmean(train_data_logic_CS_)*100,10);
se_CS_ = ESN_smooth(nanstd(train_data_logic_CS_)/sqrt(length(train_data_logic_CS_))*100,10);


subplot(6,5,5)
hold on
plot(inds_span, mean_d_tip_, 'Color', [0.7 0.7 0.7],'LineWidth', 2)
plot(inds_span, mean_d_tip_ + se_d_tip_, 'Color', [0.7 0.7 0.7],'LineWidth',1)
plot(inds_span, mean_d_tip_ - se_d_tip_, 'Color', [0.7 0.7 0.7],'LineWidth', 1)
title('Displacement | R licks')
xlabel(xlabel_text_raster_);
xlim([-300 300])
ylabel('Displacement (mm)')
if (isnan(max(mean_d_tip_)))
    ylim([0 10])
else
    ylim([0 (max(mean_d_tip_) + 1)])
end
set(gca, 'YColor', 'k')
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end

subplot(6,5,10)
N = length(mean_d_tip_);
xdft = fft(mean_d_tip_);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(mean_d_tip_):Fs/2;
plot(freq,10*log10(psdx), 'Color', [0.7 0.7 0.7], 'LineWidth', 2)
title('Periodogram | R Licks')
xlabel('Freq (Hz)')
xlim([0 10])
ylabel('Pow/Freq (dB/Hz)')


subplot(6,5,15)
hold on
plot(inds_span, mean_SS_, 'Color', 'k','LineWidth', 2)
plot(inds_span, mean_SS_ + se_SS_, 'Color', 'k' ,'LineWidth',1)
plot(inds_span, mean_SS_ - se_SS_, 'Color', 'k' ,'LineWidth', 1)
title('SS firing | R licks')
xlabel(xlabel_text_raster_);
xlim([-300 300])
ylabel('SS firing (spks/s)')
if (isnan(max(mean_SS_)))
    ylim([0 300])
else
    ylim([0 (max(mean_SS_) + max(mean_SS_)/2)])
end
if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end

subplot(6,5,20)
N = length(mean_SS_);
xdft = fft(mean_SS_);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(mean_SS_):Fs/2;
plot(freq,10*log10(psdx), 'k', 'LineWidth', 2)
title('Periodogram | R Licks - SS')
xlabel('Freq (Hz)')
xlim([0 10])
ylabel('Pow/Freq (dB/Hz)')

subplot(6,5,25)
hold on
plot(inds_span, mean_CS_, 'Color', 'k','LineWidth', 2)
plot(inds_span, mean_CS_ + se_CS_, 'Color', 'k' ,'LineWidth',1)
plot(inds_span, mean_CS_ - se_CS_, 'Color', 'k' ,'LineWidth', 1)
title('CS firing | R licks')
xlabel(xlabel_text_raster_);
xlim([-300 300])
ylabel('CS firing (spks/s)')
if (isnan(max(mean_SS_)))
    ylim([0 3])
else
    ylim([0 (max(mean_CS_) + max(mean_CS_)/2)])
end

if contains(xlabel_text_raster_, 'onset')
    xline(0,'LineWidth', 2, 'Color', color_lick_onset);
elseif contains(xlabel_text_raster_, 'dmax')
    xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
elseif contains(xlabel_text_raster_, 'offset')
    xline(0,'LineWidth', 2, 'Color', color_lick_offset);
end

subplot(6,5,30)
N = length(mean_CS_);
xdft = fft(mean_CS_);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(mean_CS_):Fs/2;
plot(freq,10*log10(psdx), 'k', 'LineWidth', 2)
title('Periodogram | R Licks - CS')
xlabel('Freq (Hz)')
xlim([0 10])
ylabel('Pow/Freq (dB/Hz)')

%% L Licks
if (session_type  == 2)
    train_data_logic_SS_ = train_data_logic_SS_l;
    train_data_logic_CS_ = train_data_logic_CS_l;
    d_tip_ = VID_d_tip_100_l;
    
    mean_d_tip_ = nanmean(d_tip_);
    se_d_tip_ = nanstd(d_tip_)/sqrt(length(d_tip_));
    mean_SS_ = ESN_smooth(nanmean(train_data_logic_SS_)*100,5);
    se_SS_ = ESN_smooth(nanstd(train_data_logic_SS_)/sqrt(length(train_data_logic_SS_))*100,5);
    mean_CS_ = ESN_smooth(nanmean(train_data_logic_CS_)*100,10);
    se_CS_ = ESN_smooth(nanstd(train_data_logic_CS_)/sqrt(length(train_data_logic_CS_))*100,10);
    
    
    subplot(6,5,3)
    hold on
    plot(inds_span, mean_d_tip_, 'Color', [0.7 0.7 0.7],'LineWidth', 2)
    plot(inds_span, mean_d_tip_ + se_d_tip_, 'Color', [0.7 0.7 0.7],'LineWidth',1)
    plot(inds_span, mean_d_tip_ - se_d_tip_, 'Color', [0.7 0.7 0.7],'LineWidth', 1)
    title('Displacement | L licks')
    xlabel(xlabel_text_raster_);
    xlim([-300 300])
    ylabel('Displacement (mm)')
    if (isnan(max(mean_d_tip_)))
        ylim([0 10])
    else
        ylim([0 (max(mean_d_tip_) + 1)])
    end
    set(gca, 'YColor', 'k')
    if contains(xlabel_text_raster_, 'onset')
        xline(0,'LineWidth', 2, 'Color', color_lick_onset);
    elseif contains(xlabel_text_raster_, 'dmax')
        xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
    elseif contains(xlabel_text_raster_, 'offset')
        xline(0,'LineWidth', 2, 'Color', color_lick_offset);
    end
    
    subplot(6,5,8)
    N = length(mean_d_tip_);
    xdft = fft(mean_d_tip_);
    xdft = xdft(1:N/2+1);
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:Fs/length(mean_d_tip_):Fs/2;
    plot(freq,10*log10(psdx), 'Color', [0.7 0.7 0.7], 'LineWidth', 2)
    title('Periodogram | L Licks')
    xlabel('Freq (Hz)')
    xlim([0 10])
    ylabel('Pow/Freq (dB/Hz)')
    
    
    subplot(6,5,13)
    hold on
    plot(inds_span, mean_SS_, 'Color', 'k','LineWidth', 2)
    plot(inds_span, mean_SS_ + se_SS_, 'Color', 'k' ,'LineWidth',1)
    plot(inds_span, mean_SS_ - se_SS_, 'Color', 'k' ,'LineWidth', 1)
    title('SS firing | L licks')
    xlabel(xlabel_text_raster_);
    xlim([-300 300])
    ylabel('SS firing (spks/s)')
    if (isnan(max(mean_SS_)))
        ylim([0 300])
    else
        ylim([0 (max(mean_SS_) + max(mean_SS_)/2)])
    end
    if contains(xlabel_text_raster_, 'onset')
        xline(0,'LineWidth', 2, 'Color', color_lick_onset);
    elseif contains(xlabel_text_raster_, 'dmax')
        xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
    elseif contains(xlabel_text_raster_, 'offset')
        xline(0,'LineWidth', 2, 'Color', color_lick_offset);
    end
    
    subplot(6,5,18)
    N = length(mean_SS_);
    xdft = fft(mean_SS_);
    xdft = xdft(1:N/2+1);
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:Fs/length(mean_SS_):Fs/2;
    plot(freq,10*log10(psdx), 'k', 'LineWidth', 2)
    title('Periodogram | L Licks - SS')
    xlabel('Freq (Hz)')
    xlim([0 10])
    ylabel('Pow/Freq (dB/Hz)')
    
    subplot(6,5,23)
    hold on
    plot(inds_span, mean_CS_, 'Color', 'k','LineWidth', 2)
    plot(inds_span, mean_CS_ + se_CS_, 'Color', 'k' ,'LineWidth',1)
    plot(inds_span, mean_CS_ - se_CS_, 'Color', 'k' ,'LineWidth', 1)
    title('CS firing | L licks')
    xlabel(xlabel_text_raster_);
    xlim([-300 300])
    ylabel('CS firing (spks/s)')
    if (isnan(max(mean_CS_)) || max(mean_CS_) == 0 )
        ylim([0 3])
    else
        ylim([0 (max(mean_CS_) + max(mean_CS_)/2)])
    end
    
    if contains(xlabel_text_raster_, 'onset')
        xline(0,'LineWidth', 2, 'Color', color_lick_onset);
    elseif contains(xlabel_text_raster_, 'dmax')
        xline(0,'LineWidth', 2, 'Color', color_lick_dmax);
    elseif contains(xlabel_text_raster_, 'offset')
        xline(0,'LineWidth', 2, 'Color', color_lick_offset);
    end
    
    subplot(6,5,28)
    N = length(mean_CS_);
    xdft = fft(mean_CS_);
    xdft = xdft(1:N/2+1);
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:Fs/length(mean_CS_):Fs/2;
    plot(freq,10*log10(psdx), 'k', 'LineWidth', 2)
    title('Periodogram | L Licks - CS')
    xlabel('Freq (Hz)')
    xlim([0 10])
    ylabel('Pow/Freq (dB/Hz)')
    
end

%% All licks phase relationships
train_data_logic_SS_ = train_data_logic_SS_all;
train_data_logic_CS_ = train_data_logic_CS_all;
d_tip_ = VID_d_tip_100_all;

mean_d_tip_ = nanmean(d_tip_);
se_d_tip_ = nanstd(d_tip_)/sqrt(length(d_tip_));
mean_SS_ = ESN_smooth(nanmean(train_data_logic_SS_)*100,5);
se_SS_ = ESN_smooth(nanstd(train_data_logic_SS_)/sqrt(length(train_data_logic_SS_))*100,5);
mean_CS_ = ESN_smooth(nanmean(train_data_logic_CS_)*100,10);
se_CS_ = ESN_smooth(nanstd(train_data_logic_CS_)/sqrt(length(train_data_logic_CS_))*100,10);

phase_diff = rad2deg(phdiffmeasure(mean_d_tip_,mean_SS_));
[max_mean_d_tip,ind_max_mean_d_tip] = max(mean_d_tip_);
[max_mean_SS, ind_max_mean_SS] = max(mean_SS_);
time_diff_peak = (ind_max_mean_SS - ind_max_mean_d_tip)*10;

subplot(6,5,[1 2 6 7])
hold on
yyaxis left;
plot(inds_span, mean_SS_, 'Color', 'k','LineWidth', 2)
plot(inds_span, mean_SS_ + se_SS_, 'Color', 'k' ,'LineWidth',1)
plot(inds_span, mean_SS_ - se_SS_, 'Color', 'k' ,'LineWidth', 1)
plot(inds_span(ind_max_mean_SS), max_mean_SS, 'ok', 'LineWidth', 5)
xlim([-300 300])
ylabel('SS firing (spks/s)')
set(gca, 'YColor', 'k')
if (isnan(max(mean_SS_)))
    ylim([0 300])
else
    ylim([0 (max(mean_SS_) + max(mean_SS_)/2)])
end
yyaxis right;
plot(inds_span, mean_d_tip_, 'Color', [0.7 0.7 0.7],'LineWidth', 2)
plot(inds_span, mean_d_tip_ + se_d_tip_, 'Color', [0.7 0.7 0.7],'LineWidth',1)
plot(inds_span, mean_d_tip_ - se_d_tip_, 'Color', [0.7 0.7 0.7],'LineWidth', 1)
plot(inds_span(ind_max_mean_d_tip), max_mean_d_tip, 'o', 'Color' , [0.7 0.7 0.7], 'LineWidth', 5)
title('All licks')
if (isnan(max(mean_d_tip_)))
    ylim([0 10])
else
    ylim([0 (max(mean_d_tip_) + 1)])
end
xline(0,'--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2 );
set(gca, 'YColor', [0.7 0.7 0.7])
ylabel('Displacement (mm)')


subplot(6,5,[11 12 16 17])
N_mean_d_tip = length(mean_d_tip_);
xdft_mean_d_tip = fft(mean_d_tip_);
xdft_mean_d_tip = xdft_mean_d_tip(1:N_mean_d_tip/2+1);
psdx_mean_d_tip = (1/(Fs*N_mean_d_tip)) * abs(xdft_mean_d_tip).^2;
psdx_mean_d_tip(2:end-1) = 2*psdx_mean_d_tip(2:end-1);
N_mean_SS = length(mean_SS_);
xdft_mean_SS = fft(mean_SS_);
xdft_mean_SS = xdft_mean_SS(1:N_mean_SS/2+1);
psdx_mean_SS = (1/(Fs*N_mean_SS)) * abs(xdft_mean_SS).^2;
psdx_mean_SS(2:end-1) = 2*psdx_mean_SS(2:end-1);
freq_mean_d_tip = 0:Fs/length(mean_d_tip_):Fs/2;
freq_mean_SS = 0:Fs/length(mean_SS_):Fs/2;
hold on;
plot(freq_mean_d_tip,10*log10(psdx_mean_d_tip), 'Color', [0.7 0.7 0.7], 'LineWidth', 2)
plot(freq_mean_SS,10*log10(psdx_mean_SS), 'k', 'LineWidth', 2)
xlabel('Freq (Hz)')
xlim([0 10])
ylabel('Pow/Freq (dB/Hz)')

subplot(6,5,[21 22 26 27])
x_axis = categorical({'Phase diff (deg)', 'Time diff peak (ms)'});
x_axis = reordercats(x_axis, {'Phase diff (deg)', 'Time diff peak (ms)'});
yyaxis left;
bar(x_axis(1),phase_diff, 'k' )
ylabel('Phase diff (deg)')
if(phase_diff > 0)
    ylim([0 (phase_diff+40)])
else
    ylim([(phase_diff-40) 0])
end
set(gca, 'YColor', 'k')

yyaxis right;
set(gca, 'YColor', 'k')
bar(x_axis(2),time_diff_peak, 'k' )
ylabel('Time diff (ms)')
if(time_diff_peak > 0)
    ylim([0 (time_diff_peak+40)])
else
    ylim([(time_diff_peak-40) 0])
end
set(gca, 'YColor', 'k')


end
%% ------------------------------ %%



