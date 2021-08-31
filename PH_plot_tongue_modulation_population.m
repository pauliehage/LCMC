clc; clear all; close all;
%% LOAD DATA
[file_name, file_path] = uigetfile([pwd], 'Select ALL_PCELL_COMPRESSED_DATA file');
load([file_path file_name ] )

%% FIND SIMULTANEOUS SESSIONS

%% SPECIFY CELL LOCATION

%% SPECIFY SUBJECT

%% LICK PROPERTIES
clearvars -except ALL_PCELL_COMPRESSED_DATA POPULATION

alignment = 'dmax';
% alignment = 'onset';

for counter_pCell = 1: length(ALL_PCELL_COMPRESSED_DATA)
    %     if (sum(contains(fields(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset), 'onset_l'), 1) > 0 )
    id(counter_pCell, 1) = string(ALL_PCELL_COMPRESSED_DATA(counter_pCell).id{1, 1});
    if contains(alignment,'onset') == 1
        % COUNTS
        num_g(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_grooming_corr_numLick;
        num_r(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_r_corr_numLick;
        num_l(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_l_corr_numLick;
        num_all(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_grooming_corr_numLick + ...
            ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_r_corr_numLick + ...
            ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_l_corr_numLick;
        
        num_bout_r(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_str_bout_r_numLick;
        num_bout_l(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_str_bout_l_numLick;
        num_bout_all(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_str_bout_r_numLick + ...
            ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_str_bout_l_numLick;
        
    elseif contains(alignment,'dmax') == 1
        num_g(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_grooming_corr_numLick;
        num_r(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_r_corr_numLick;
        num_l(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_l_corr_numLick;
        num_all(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_grooming_corr_numLick + ...
            ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_r_corr_numLick + ...
            ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_l_corr_numLick;
        
        num_bout_r(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_str_bout_r_numLick;
        num_bout_l(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_str_bout_l_numLick;
        num_bout_all(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_str_bout_r_numLick + ...
            ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_str_bout_l_numLick;
    end
    
    % ILI & ILR
    %filter
    ILI_all_ = ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_ILI_all ;
    ILI_all_(ILI_all_ > 0.5) = nan;
    ILI_g_ = ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_ILI_grooming ;
    ILI_g_(ILI_g_ > 0.5) = nan;
    ILI_r_ = ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_ILI_r ;
    ILI_r_(ILI_r_ > 0.5) = nan;
    ILI_l_ = ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_ILI_l ;
    ILI_l_(ILI_l_ > 0.5) = nan;
    ILR_all_ = ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_ILR_all ;
    ILR_all_(ILR_all_ < 2 ) = nan;
    ILR_g_ = ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_ILR_grooming ;
    ILR_g_(ILR_g_ < 2) = nan;
    ILR_r_ = ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_ILR_r ;
    ILR_r_(ILR_r_ < 2) = nan;
    ILR_l_ = ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_ILR_l ;
    ILR_l_(ILR_l_ < 2) = nan;
    
    ILI_all(counter_pCell,1) = nanmean(ILI_all_);
    ILI_g(counter_pCell,1) = nanmean(ILI_g_);
    ILI_r(counter_pCell,1) = nanmean(ILI_r_);
    ILI_l(counter_pCell,1) = nanmean(ILI_l_);
    ILR_all(counter_pCell,1) = nanmean(ILR_all_);
    ILR_g(counter_pCell,1) = nanmean(ILR_g_);
    ILR_r(counter_pCell,1) = nanmean(ILR_r_);
    ILR_l(counter_pCell,1) = nanmean(ILR_l_);
    
    % bout duration and # of licks
    % filter
    num_lick_bout_all_ = ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_num_lick_bout;
    num_lick_bout_all_(num_lick_bout_all_<3) = nan;
    num_lick_bout_r_ = ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_num_lick_bout_r;
    num_lick_bout_r_(num_lick_bout_r_<3) = nan;
    num_lick_bout_l_ = ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_num_lick_bout_l;
    num_lick_bout_l_(num_lick_bout_l_<3) = nan;
    num_lick_bout_all(counter_pCell,1) = nanmean(num_lick_bout_all_);
    num_lick_bout_r(counter_pCell,1) = nanmean(num_lick_bout_r_);
    num_lick_bout_l(counter_pCell,1) = nanmean(num_lick_bout_l_);
    
    duration_bout_all_ = ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_bout_duration;
    duration_bout_all_(duration_bout_all_<0.8) = nan;
    duration_bout_r_ = ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_bout_duration_r;
    duration_bout_r_(duration_bout_r_<0.8) = nan;
    duration_bout_l_ = ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_bout_duration_l;
    duration_bout_l_(duration_bout_l_<0.8) = nan;
    duration_bout_all(counter_pCell,1) = nanmean(duration_bout_all_);
    duration_bout_r(counter_pCell,1) = nanmean(duration_bout_r_);
    duration_bout_l(counter_pCell,1) = nanmean(duration_bout_l_);
    
    % lick durations
    % filter
    duration_all_ = ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_lick_duration_all;
    duration_all_(duration_all_ > 0.6) = nan;
    duration_g_ = ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_lick_duration_grooming;
    duration_g_(duration_g_ > 0.6) = nan;
    duration_r_ = ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_lick_duration_r;
    duration_r_(duration_r_ > 0.6) = nan;
    duration_l_ = ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_lick_duration_l;
    duration_l_(duration_l_ > 0.6) = nan;
    
    duration_all(counter_pCell,1) = nanmean(duration_all_);
    duration_g(counter_pCell,1) = nanmean(duration_g_);
    duration_r(counter_pCell,1) = nanmean(duration_r_);
    duration_l(counter_pCell,1) = nanmean(duration_l_);
    
    duration_d_max_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_duration_d_max_100_grooming);
    duration_d_max_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_duration_d_max_100_r);
    duration_d_max_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_duration_d_max_100_l);
    
    duration_v_max_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_duration_v_max_100_grooming);
    duration_v_max_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_duration_v_max_100_r);
    duration_v_max_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_duration_v_max_100_l);
    
    duration_v_min_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_duration_v_min_100_grooming);
    duration_v_min_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_duration_v_min_100_r);
    duration_v_min_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_duration_v_min_100_l);
end

sum_g = nansum(num_g);
sum_r = nansum(num_r);
sum_l = nansum(num_l);
sum_bout_r = nansum(num_bout_r);
sum_bout_l = nansum(num_bout_l);

%%% PLOTS %%%
figure
subplot(4,7,1)
hold on;
edges_all = (0 : 500 : max(num_all)+500);
histogram((num_all), edges_all, 'FaceColor', 'k');
xline(nanmean(num_all), 'k', 'LineWidth', 2);
xlabel('Num. of licks')
ylabel('Count')
title(['All licks | Total: ' num2str(sum_g+sum_r+sum_l)])

subplot(4,7,3)
hold on;
% edges_grl = (0 : 500 : max([max(num_g) max(num_r) max(num_l)]) + 500);
edges_all = (0 : 500 : max(num_all)+500);
histogram((num_g), edges_all, 'FaceColor', 'g');
xline(nanmean(num_g), 'g', 'LineWidth', 2);
xlabel('Num. of licks')
ylabel('Count')
title(['G licks | Total: '  num2str(sum_g)])

subplot(4,7,4)
hold on;
edges_all = (0 : 500 : max(num_all)+500);
histogram((num_r), edges_all, 'FaceColor', 'r');
xline(nanmean(num_r), 'r', 'LineWidth', 2);
xlabel('Num. of licks')
ylabel('Count')
title(['R licks | Total: ' num2str(sum_r)])

subplot(4,7,2)
hold on;
edges_all = (0 : 500 : max(num_all)+500);
histogram((num_l), edges_all, 'FaceColor', 'b');
xline(nanmean(num_l), 'b', 'LineWidth', 2);
xlabel('Num. of licks')
ylabel('Count')
title(['L licks | Total: ' num2str(sum_l)])

subplot(4,7,5)
hold on;
edges_all = (0 : 25 : max(num_bout_all)+50);
histogram((num_bout_all), edges_all, 'FaceColor', 'k');
xline(nanmean(num_bout_all), 'k', 'LineWidth', 2);
xlabel('Num. of bouts')
ylabel('Count')
title(['All bouts | Total: ' num2str(sum_bout_l+sum_bout_r) ])

subplot(4,7,6)
hold on;
edges_all = (0 : 25 : max(num_bout_all)+50);
% edges_bout_grl = (0 : 50 : max([max(num_bout_r) max(num_bout_l)]) + 50);
histogram((num_bout_l), edges_all, 'FaceColor', 'b');
xline(nanmean(num_bout_l), 'b', 'LineWidth', 2);
xlabel('Num. of bouts')
ylabel('Count')
title(['L bouts | Total: ' num2str(sum_bout_l)])

subplot(4,7,7)
hold on;
edges_all = (0 : 25 : max(num_bout_all)+50);
histogram((num_bout_r), edges_all, 'FaceColor', 'r');
xline(nanmean(num_bout_r), 'r', 'LineWidth', 2);
xlabel('Num. of bouts')
ylabel('Count')
title(['R bouts | Total: ' num2str(sum_bout_r)])

subplot(4,7,8)
hold on;
edges_duration_all = (0.2 : 0.01 : 0.4);
histogram((duration_all), edges_duration_all, 'FaceColor', 'k');
xline(nanmean(duration_all), 'k', 'LineWidth', 2);
xlabel('Duration (s)')
ylabel('Count')

subplot(4,7,9)
hold on;
edges_duration_l = (0.2 : 0.01 : 0.4);
histogram((duration_l), edges_duration_l, 'FaceColor', 'b');
xline(nanmean(duration_l), 'b', 'LineWidth', 2);
xlabel('Duration (s)')
ylabel('Count')

subplot(4,7,10)
hold on;
edges_duration_g = (0.2 : 0.01 : 0.4);
histogram((duration_g), edges_duration_g, 'FaceColor', 'g');
xline(nanmean(duration_g), 'g', 'LineWidth', 2);
xlabel('Duration (s)')
ylabel('Count')

subplot(4,7,11)
hold on;
edges_duration_r = (0.2 : 0.01 : 0.4);
histogram((duration_r), edges_duration_r, 'FaceColor', 'r');
xline(nanmean(duration_r), 'r', 'LineWidth', 2);
xlabel('Duration (s)')
ylabel('Count')

subplot(4,7,15)
hold on;
edges_ILI_all = (0.2 : 0.01 : 0.4);
histogram((ILI_all), edges_ILI_all, 'FaceColor', 'k');
xline(nanmean(ILI_all), 'k', 'LineWidth', 2);
xlabel('ILI (s)')
ylabel('Count')

subplot(4,7,16)
hold on;
edges_ILI_l = (0.2 : 0.01 : 0.4);
histogram((ILI_l), edges_ILI_l, 'FaceColor', 'b');
xline(nanmean(ILI_l), 'b', 'LineWidth', 2);
xlabel('ILI (s)')
ylabel('Count')

subplot(4,7,17)
hold on;
edges_ILI_g = (0.2 : 0.01 : 0.4);
histogram((ILI_g), edges_ILI_g, 'FaceColor', 'g');
xline(nanmean(ILI_g), 'g', 'LineWidth', 2);
xlabel('ILI (s)')
ylabel('Count')

subplot(4,7,18)
hold on;
edges_ILI_r = (0.2 : 0.01 : 0.4);
histogram((ILI_r), edges_ILI_r, 'FaceColor', 'r');
xline(nanmean(ILI_r), 'r', 'LineWidth', 2);
xlabel('ILI (s)')
ylabel('Count')
title('R licks')

subplot(4,7,22)
hold on;
edges_ILR_all = (2.5 : 0.1 : 4);
histogram((ILR_all), edges_ILR_all, 'FaceColor', 'k');
xline(nanmean(ILR_all), 'k', 'LineWidth', 2);
xlabel('ILR (Hz)')
ylabel('Count')

subplot(4,7,23)
hold on;
edges_ILR_l = (2.5 : 0.1 : 4);
histogram((ILR_l), edges_ILR_l, 'FaceColor', 'b');
xline(nanmean(ILR_l), 'b', 'LineWidth', 2);
xlabel('ILR (Hz)')
ylabel('Count')

subplot(4,7,24)
hold on;
edges_ILR_g = (2.5 : 0.1 : 4);
histogram((ILR_g), edges_ILR_g, 'FaceColor', 'g');
xline(nanmean(ILR_g), 'g', 'LineWidth', 2);
xlabel('ILR (Hz)')
ylabel('Count')

subplot(4,7,25)
hold on;
edges_ILR_r = (2.5 : 0.1 : 4);
histogram((ILR_r), edges_ILR_r, 'FaceColor', 'r');
xline(nanmean(ILR_r), 'r', 'LineWidth', 2);
xlabel('ILR (Hz)')
ylabel('Count')

subplot(4,7,12)
hold on;
edges_num_lick_bout_all = (0 : 5 : max(num_lick_bout_all) + 10);
histogram((num_lick_bout_all), edges_num_lick_bout_all, 'FaceColor', 'k');
xline(nanmean(num_lick_bout_all), 'k', 'LineWidth', 2);
xlabel('Num. of licks in a bout')
ylabel('Count')

subplot(4,7,13)
hold on;
edges_num_lick_bout_l = (0 : 5 : max(num_lick_bout_l) + 10);
histogram((num_lick_bout_l), edges_num_lick_bout_l, 'FaceColor',' b');
xline(nanmean(num_lick_bout_l), 'b', 'LineWidth', 2);
xlabel('Num. of licks in a bout')
ylabel('Count')

subplot(4,7,14)
hold on;
edges_num_lick_bout_r = (0 : 5 : max(num_lick_bout_r) + 10);
histogram((num_lick_bout_r), edges_num_lick_bout_r, 'FaceColor', 'r');
xline(nanmean(num_lick_bout_r), 'r', 'LineWidth', 2);
xlabel('Num. of licks in a bout')
ylabel('Count')

subplot(4,7,19)
hold on;
edges_duration_bout_all = (0 : 1 : max(duration_bout_all) + 1);
histogram((duration_bout_all), edges_duration_bout_all, 'FaceColor', 'k');
xline(nanmean(duration_bout_all), 'k', 'LineWidth', 2);
xlabel('Duration of a bout (s)')
ylabel('Count')

subplot(4,7,20)
hold on;
edges_duration_bout_l = (0 : 1 : max(duration_bout_l) + 1);
histogram((duration_bout_l), edges_duration_bout_l, 'FaceColor', 'b');
xline(nanmean(duration_bout_l), 'b', 'LineWidth', 2);
xlabel('Duration of a bout (s)')
ylabel('Count')

subplot(4,7,21)
hold on;
edges_duration_bout_r = (0 : 1 : max(duration_bout_r) + 1);
histogram((duration_bout_r), edges_duration_bout_r, 'FaceColor', 'r');
xline(nanmean(duration_bout_r), 'r', 'LineWidth', 2);
xlabel('Duration of a bout (s)')
ylabel('Count')

%% KINEMATIC PROPERTIES & TRACES
clearvars -except ALL_PCELL_COMPRESSED_DATA POPULATION

num_pCell = 1:62;

alignment = 'dmax';
% alignment = 'onset';

for counter_pCell = num_pCell
    id(counter_pCell, 1) = string(ALL_PCELL_COMPRESSED_DATA(counter_pCell).id{1, 1});
    if contains(alignment,'onset') == 1
        %%% CALCULATE MEANS OF KINEMATIC TRACES %%%
        %%% d_tip %%%
        d_tip_str_bout_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_str_bout);
        d_tip_str_bout_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_str_bout_r);
        d_tip_str_bout_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_str_bout_l);
        d_tip_end_bout_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_end_bout);
        d_tip_end_bout_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_end_bout_r);
        d_tip_end_bout_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_end_bout_l);
        d_tip_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_all);
        d_tip_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_grooming);
        d_tip_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_r);
        d_tip_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_l);
        d_max_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_d_max_100_all);
        d_max_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_d_max_100_grooming);
        d_max_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_d_max_100_r);
        d_max_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_d_max_100_l);
        
        %%% v_tip %%%
        v_tip_str_bout_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_v_tip_100_str_bout);
        v_tip_str_bout_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_v_tip_100_str_bout_r);
        v_tip_str_bout_l(counter_pCell,:) =  nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_v_tip_100_str_bout_l);
        v_tip_end_bout_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_v_tip_100_end_bout);
        v_tip_end_bout_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_v_tip_100_end_bout_r);
        v_tip_end_bout_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_v_tip_100_end_bout_l);
        v_tip_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_v_tip_100_all);
        v_tip_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_v_tip_100_grooming);
        v_tip_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_v_tip_100_r);
        v_tip_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_v_tip_100_l);
        v_max_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_v_max_100_all);
        v_max_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_v_max_100_grooming);
        v_max_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_v_max_100_r);
        v_max_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_v_max_100_l);
        v_min_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_v_min_100_all);
        v_min_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_v_min_100_grooming);
        v_min_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_v_min_100_r);
        v_min_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_v_min_100_l);
        
        %%% ang_tip %%%
        ang_tip_str_bout_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_angle_tip_100_str_bout);
        ang_tip_str_bout_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_angle_tip_100_str_bout_r);
        ang_tip_str_bout_l(counter_pCell,:) = nanmean(abs(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_angle_tip_100_str_bout_l));
        ang_tip_end_bout_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_angle_tip_100_end_bout);
        ang_tip_end_bout_r(counter_pCell,:) =nanmean( ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_angle_tip_100_end_bout_r);
        ang_tip_end_bout_l(counter_pCell,:) = nanmean(abs(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_angle_tip_100_end_bout_l));
        ang_tip_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_angle_tip_100_all);
        ang_tip_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_angle_tip_100_grooming);
        ang_tip_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_angle_tip_100_r);
        ang_tip_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_angle_tip_100_l);
        ang_tip_max_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_angle_max_100_all);
        ang_tip_max_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_angle_max_100_grooming);
        ang_tip_max_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_angle_max_100_r);
        ang_tip_max_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_angle_max_100_l);
        
        %%% Counts %%%
        num_g(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_grooming_corr_numLick;
        num_r(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_r_corr_numLick;
        num_l(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_l_corr_numLick;
        num_all(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_grooming_corr_numLick + ...
            ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_r_corr_numLick + ...
            ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_l_corr_numLick;
        
        num_bout_r(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_str_bout_r_numLick;
        num_bout_l(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_str_bout_l_numLick;
        num_bout_all(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_str_bout_r_numLick + ...
            ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_str_bout_l_numLick;
    elseif contains(alignment,'dmax') == 1
        
        %%% CALCULATE MEANS OF KINEMATIC TRACES %%%
        %%% d_tip %%%
        d_tip_str_bout_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_str_bout);
        d_tip_str_bout_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_str_bout_r);
        d_tip_str_bout_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_str_bout_l);
        d_tip_end_bout_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_end_bout);
        d_tip_end_bout_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_end_bout_r);
        d_tip_end_bout_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_end_bout_l);
        d_tip_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_all);
        d_tip_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_grooming);
        d_tip_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_r);
        d_tip_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_l);
        d_max_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_d_max_100_all);
        d_max_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_d_max_100_grooming);
        d_max_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_d_max_100_r);
        d_max_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_d_max_100_l);
        
        %%% v_tip %%%
        v_tip_str_bout_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_v_tip_100_str_bout);
        v_tip_str_bout_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_v_tip_100_str_bout_r);
        v_tip_str_bout_l(counter_pCell,:) =  nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_v_tip_100_str_bout_l);
        v_tip_end_bout_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_v_tip_100_end_bout);
        v_tip_end_bout_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_v_tip_100_end_bout_r);
        v_tip_end_bout_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_v_tip_100_end_bout_l);
        v_tip_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_v_tip_100_all);
        v_tip_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_v_tip_100_grooming);
        v_tip_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_v_tip_100_r);
        v_tip_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_v_tip_100_l);
        v_max_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_v_max_100_all);
        v_max_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_v_max_100_grooming);
        v_max_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_v_max_100_r);
        v_max_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_v_max_100_l);
        v_min_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_v_min_100_all);
        v_min_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_v_min_100_grooming);
        v_min_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_v_min_100_r);
        v_min_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_v_min_100_l);
        
        %%% ang_tip %%%
        ang_tip_str_bout_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_angle_tip_100_str_bout);
        ang_tip_str_bout_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_angle_tip_100_str_bout_r);
        ang_tip_str_bout_l(counter_pCell,:) = nanmean(abs(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_angle_tip_100_str_bout_l));
        ang_tip_end_bout_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_angle_tip_100_end_bout);
        ang_tip_end_bout_r(counter_pCell,:) =nanmean( ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_angle_tip_100_end_bout_r);
        ang_tip_end_bout_l(counter_pCell,:) = nanmean(abs(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_angle_tip_100_end_bout_l));
        ang_tip_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_angle_tip_100_all);
        ang_tip_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_angle_tip_100_grooming);
        ang_tip_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_angle_tip_100_r);
        ang_tip_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_angle_tip_100_l);
        ang_tip_max_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_angle_max_100_all);
        ang_tip_max_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_angle_max_100_grooming);
        ang_tip_max_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_angle_max_100_r);
        ang_tip_max_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).kinematic_data.VID_angle_max_100_l);
        
        %%% Counts %%%
        num_g(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_grooming_corr_numLick;
        num_r(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_r_corr_numLick;
        num_l(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_l_corr_numLick;
        num_all(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_grooming_corr_numLick + ...
            ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_r_corr_numLick + ...
            ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_l_corr_numLick;
        
        num_bout_r(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_str_bout_r_numLick;
        num_bout_l(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_str_bout_l_numLick;
        num_bout_all(counter_pCell,1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_str_bout_r_numLick + ...
            ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_str_bout_l_numLick;
    end
end

%%% inds_span %%%
inds_span = -290:10:300;

%%% inds_span_bout %%%
inds_span_bout = -0.99:0.01:1;

%%% Sums %%%
sum_g = nansum(num_g);
sum_r = nansum(num_r);
sum_l = nansum(num_l);
sum_bout_r = nansum(num_bout_r);
sum_bout_l = nansum(num_bout_l);

%%% CALCULATE MEAN & SE TRACES %%%
%%% d_tip %%%
% str
mean_d_tip_str_bout_all = nanmean(d_tip_str_bout_all);
se_d_tip_str_bout_all = std(d_tip_str_bout_all)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
mean_d_tip_str_bout_r = nanmean(d_tip_str_bout_r);
se_d_tip_str_bout_r = std(d_tip_str_bout_r)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
mean_d_tip_str_bout_l = nanmean(d_tip_str_bout_l);
se_d_tip_str_bout_l = std(d_tip_str_bout_l)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
% end
mean_d_tip_end_bout_all = nanmean(d_tip_end_bout_all);
se_d_tip_end_bout_all = std(d_tip_end_bout_all)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
mean_d_tip_end_bout_r = nanmean(d_tip_end_bout_r);
se_d_tip_end_bout_r = std(d_tip_end_bout_r)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
mean_d_tip_end_bout_l = nanmean(d_tip_end_bout_l);
se_d_tip_end_bout_l = std(d_tip_end_bout_l)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
% lick
mean_d_tip_all = nanmean(d_tip_all);
se_d_tip_all = std(d_tip_all)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
mean_d_tip_g = nanmean(d_tip_g);
se_d_tip_g = std(d_tip_g)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
mean_d_tip_r = nanmean(d_tip_r);
se_d_tip_r = std(d_tip_r)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
mean_d_tip_l = nanmean(d_tip_l);
se_d_tip_l = std(d_tip_l)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));

%%% v_tip %%%
% str
mean_v_tip_str_bout_all = nanmean(v_tip_str_bout_all);
se_v_tip_str_bout_all = std(v_tip_str_bout_all)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
mean_v_tip_str_bout_r = nanmean(v_tip_str_bout_r);
se_v_tip_str_bout_r = std(v_tip_str_bout_r)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
mean_v_tip_str_bout_l = nanmean(v_tip_str_bout_l);
se_v_tip_str_bout_l = std(v_tip_str_bout_l)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
% end
mean_v_tip_end_bout_all = nanmean(v_tip_end_bout_all);
se_v_tip_end_bout_all = std(v_tip_end_bout_all)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
mean_v_tip_end_bout_r = nanmean(v_tip_end_bout_r);
se_v_tip_end_bout_r = std(v_tip_end_bout_r)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
mean_v_tip_end_bout_l = nanmean(v_tip_end_bout_l);
se_v_tip_end_bout_l = std(v_tip_end_bout_r)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
% lick
mean_v_tip_all = nanmean(v_tip_all);
se_v_tip_all = std(v_tip_all)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
mean_v_tip_g = nanmean(v_tip_g);
se_v_tip_g = std(v_tip_g)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
mean_v_tip_r = nanmean(v_tip_r);
se_v_tip_r = std(v_tip_r)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
mean_v_tip_l = nanmean(v_tip_l);
se_v_tip_l = std(v_tip_l)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));

%%% ang_tip %%%
mean_ang_tip_str_bout_all = nanmean(ang_tip_str_bout_all);
se_ang_tip_str_bout_all = std(ang_tip_str_bout_all)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
mean_ang_tip_str_bout_r = nanmean(ang_tip_str_bout_r);
se_ang_tip_str_bout_r = std(ang_tip_str_bout_r)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
mean_ang_tip_str_bout_l = nanmean(ang_tip_str_bout_l);
se_ang_tip_str_bout_l = std(ang_tip_str_bout_l)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
% end
mean_ang_tip_end_bout_all = nanmean(ang_tip_end_bout_all);
se_ang_tip_end_bout_all = std(ang_tip_end_bout_all)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
mean_ang_tip_end_bout_r = nanmean(ang_tip_end_bout_r);
se_ang_tip_end_bout_r = std(ang_tip_end_bout_r)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
mean_ang_tip_end_bout_l = nanmean(ang_tip_end_bout_l);
se_ang_tip_end_bout_l = std(ang_tip_end_bout_r)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
% lick
mean_ang_tip_all = nanmean(ang_tip_all);
se_ang_tip_all = std(ang_tip_all)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
mean_ang_tip_g = nanmean(ang_tip_g);
se_ang_tip_g = std(ang_tip_g)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
mean_ang_tip_r = nanmean(ang_tip_r);
se_ang_tip_r = std(ang_tip_r)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));
mean_ang_tip_l = nanmean(ang_tip_l);
se_ang_tip_l = std(ang_tip_l)/sqrt(length(ALL_PCELL_COMPRESSED_DATA));

%%% KINEMATIC PROPERTIES %%%
figure
subplot(4,4,1)
hold on;
edges_all = (0 : 0.5 : max(d_max_all)+ 5);
histogram((d_max_all), edges_all, 'FaceColor', 'k');
xline(nanmean(d_max_all), 'k', 'LineWidth', 2);
xlabel('Num. of licks')
ylabel('Count')
title(['Max displacement | All licks'])

subplot(4,4,2)
hold on;
edges_all = (0 : 0.5 : max(d_max_all)+ 5);
histogram((d_max_l), edges_all, 'FaceColor', 'b');
xline(nanmean(d_max_l), 'b', 'LineWidth', 2);
xlabel('Num. of licks')
ylabel('Count')
title(['Max displacement | L licks'])

subplot(4,4,3)
hold on;
edges_all = (0 : 0.5 : max(d_max_all)+ 5);
histogram((d_max_g), edges_all, 'FaceColor', 'g');
xline(nanmean(d_max_g), 'g', 'LineWidth', 2);
xlabel('Num. of licks')
ylabel('Count')
title(['Max displacement | G licks'])

subplot(4,4,4)
hold on;
edges_all = (0 : 0.5 : max(d_max_all)+ 5);
histogram((d_max_r), edges_all, 'FaceColor', 'r');
xline(nanmean(d_max_r), 'r', 'LineWidth', 2);
xlabel('Num. of licks')
ylabel('Count')
title(['Max displacement | R licks'])

subplot(4,4,5)
hold on;
edges_all = (0 : 10 : max(abs(v_min_all))+ 100);
histogram((v_max_all), edges_all, 'FaceColor', 'k');
xline(nanmean(v_max_all), 'k', 'LineWidth', 2);
xlabel('Num. of licks')
xlabel('Speed (mm/s)')
title(['Max protrusion speed | All licks'])

subplot(4,4,6)
hold on;
edges_all = (0 : 10 : max(abs(v_min_all))+ 100);
histogram((v_max_l), edges_all, 'FaceColor', 'b');
xline(nanmean(v_max_l), 'b', 'LineWidth', 2);
xlabel('Speed (mm/s)')
ylabel('Count')
title(['Max protrusion speed | L licks'])

subplot(4,4,7)
hold on;
edges_all = (0 : 10 : max(abs(v_min_all))+ 100);
histogram((v_max_g), edges_all, 'FaceColor', 'g');
xline(nanmean(v_max_g), 'g', 'LineWidth', 2);
xlabel('Speed (mm/s)')
ylabel('Count')
title(['Max protrusion speed | G licks'])

subplot(4,4,8)
hold on;
edges_all = (0 : 10 : max(abs(v_min_all))+ 100);
histogram((v_max_r), edges_all, 'FaceColor', 'r');
xline(nanmean(v_max_r), 'r', 'LineWidth', 2);
xlabel('Speed (mm/s)')
ylabel('Count')
title(['Max protrusion speed | R licks'])

subplot(4,4,9)
hold on;
edges_all = (0 : 10 : max(abs(v_min_all))+ 100);
histogram((abs(v_min_all)), edges_all, 'FaceColor', 'k');
xline(nanmean(abs(v_min_all)), 'k', 'LineWidth', 2);
ylabel('Count')
xlabel('Speed (mm/s)')
title(['Max retraction speed | All licks'])

subplot(4,4,10)
hold on;
edges_all = (0 : 10 : max(abs(v_min_all))+ 100);
histogram((abs(v_min_l)), edges_all, 'FaceColor', 'b');
xline(nanmean(abs(v_min_l)), 'b', 'LineWidth', 2);
ylabel('Count')
xlabel('Speed (mm/s)')
title(['Max retraction speed | L licks'])

subplot(4,4,11)
hold on;
edges_all = (0 : 10 : max(abs(v_min_all))+ 100);
histogram((abs(v_min_g)), edges_all, 'FaceColor', 'g');
xline(nanmean(abs(v_min_g)), 'g', 'LineWidth', 2);
ylabel('Count')
xlabel('Speed (mm/s)')
title(['Max retraction speed | G licks'])

subplot(4,4,12)
hold on;
edges_all = (0 : 10 : max(abs(v_min_all))+ 100);
histogram((abs(v_min_r)), edges_all, 'FaceColor', 'r');
xline(nanmean(abs(v_min_r)), 'r', 'LineWidth', 2);
ylabel('Count')
xlabel('Speed (mm/s)')
title(['Max retraction speed | R licks'])

subplot(4,4,13)
hold on;
edges_all = (0 : 5 : 100);
histogram((abs(ang_tip_max_all)), edges_all, 'FaceColor', 'k');
xline(nanmean(abs(ang_tip_max_all)), 'k', 'LineWidth', 2);
ylabel('Count')
xlabel('Angle (deg)')
title(['Max angle | All licks'])

subplot(4,4,14)
hold on;
edges_all = (0 : 5 : 100);
histogram((abs(ang_tip_max_l)), edges_all, 'FaceColor', 'b');
xline(nanmean(abs(ang_tip_max_l)), 'b', 'LineWidth', 2);
ylabel('Count')
xlabel('Angle (deg)')
title(['Max angle | L licks'])

subplot(4,4,15)
hold on;
edges_all = (0 : 5 : 100);
histogram((abs(ang_tip_max_g)), edges_all, 'FaceColor', 'g');
xline(nanmean(abs(ang_tip_max_g)), 'g', 'LineWidth', 2);
ylabel('Count')
xlabel('Angle (deg)')
title(['Max angle | G licks'])

subplot(4,4,16)
hold on;
edges_all = (0 : 5 : 100);
histogram((abs(ang_tip_max_r)), edges_all, 'FaceColor', 'r');
xline(nanmean(abs(ang_tip_max_r)), 'r', 'LineWidth', 2);
ylabel('Count')
xlabel('Angle (deg)')
title(['Max angle | R licks'])

%%% KINEMATIC TRACES %%%
figure
subplot(2,3,1)
hold on
plot(inds_span_bout, mean_d_tip_str_bout_r(:), 'r', 'LineWidth', 2)
plot(inds_span_bout, mean_d_tip_str_bout_r + se_d_tip_str_bout_r , 'r', 'LineWidth', 1)
plot(inds_span_bout, mean_d_tip_str_bout_r - se_d_tip_str_bout_r , 'r', 'LineWidth', 1)
plot(inds_span_bout, mean_d_tip_str_bout_l, 'b', 'LineWidth', 2)
plot(inds_span_bout, mean_d_tip_str_bout_l + se_d_tip_str_bout_l , 'b', 'LineWidth', 1)
plot(inds_span_bout, mean_d_tip_str_bout_l - se_d_tip_str_bout_l , 'b', 'LineWidth', 1)
xlabel('Bout onset (s)')
ylabel('Displacement (mm)')
xlim([-0.5 1])
xline(0, 'k', 'LineWidth', 2);
subplot(2,3,2)
hold on
plot(inds_span_bout, mean_v_tip_str_bout_r, 'r', 'LineWidth', 2)
plot(inds_span_bout, mean_v_tip_str_bout_r + se_v_tip_str_bout_r , 'r', 'LineWidth', 1)
plot(inds_span_bout, mean_v_tip_str_bout_r - se_v_tip_str_bout_r , 'r', 'LineWidth', 1)
plot(inds_span_bout, mean_v_tip_str_bout_l, 'b', 'LineWidth', 2)
plot(inds_span_bout, mean_v_tip_str_bout_l + se_v_tip_str_bout_l , 'b', 'LineWidth', 1)
plot(inds_span_bout, mean_v_tip_str_bout_l - se_v_tip_str_bout_l , 'b', 'LineWidth', 1)
xlabel('Bout onset (s)')
ylabel('Velocity (mm/s)')
xlim([-0.5 1])
xline(0, 'k', 'LineWidth', 2);
title(['Number of licking bouts: ' num2str(sum_bout_l) '(L) & ' num2str(sum_bout_r) '(R) | Total: ' num2str(sum_bout_l+sum_bout_r)])
subplot(2,3,3)
hold on
plot(inds_span_bout, mean_ang_tip_str_bout_r, 'r', 'LineWidth', 2)
plot(inds_span_bout, mean_ang_tip_str_bout_r + se_ang_tip_str_bout_r , 'r', 'LineWidth', 1)
plot(inds_span_bout, mean_ang_tip_str_bout_r - se_ang_tip_str_bout_r , 'r', 'LineWidth', 1)
plot(inds_span_bout, mean_ang_tip_str_bout_l, 'b', 'LineWidth', 2)
plot(inds_span_bout, mean_ang_tip_str_bout_l + se_ang_tip_str_bout_l , 'b', 'LineWidth', 1)
plot(inds_span_bout, mean_ang_tip_str_bout_l - se_ang_tip_str_bout_l , 'b', 'LineWidth', 1)
xlabel('Bout onset (s)')
ylabel('Angle (deg)')
xlim([-0.5 1])
xline(0, 'k', 'LineWidth', 2);

subplot(2,3,4)
hold on
plot(inds_span, mean_d_tip_g, 'g', 'LineWidth', 2)
plot(inds_span, mean_d_tip_g + se_d_tip_g , 'g', 'LineWidth', 1)
plot(inds_span, mean_d_tip_g - se_d_tip_g , 'g', 'LineWidth', 1)
plot(inds_span, mean_d_tip_r, 'r', 'LineWidth', 2)
plot(inds_span, mean_d_tip_r + se_d_tip_r , 'r', 'LineWidth', 1)
plot(inds_span, mean_d_tip_r - se_d_tip_r , 'r', 'LineWidth', 1)
plot(inds_span, mean_d_tip_l, 'b', 'LineWidth', 2)
plot(inds_span, mean_d_tip_l + se_d_tip_l , 'b', 'LineWidth', 1)
plot(inds_span, mean_d_tip_l - se_d_tip_l , 'b', 'LineWidth', 1)
xlabel(['Lick ' alignment  ' (ms)'])
ylabel('Displacement (mm)')
xlim([-400 400])
xline(0, 'k', 'LineWidth', 2);

subplot(2,3,5)
hold on
plot(inds_span, mean_v_tip_g, 'g', 'LineWidth', 2)
plot(inds_span, mean_v_tip_g + se_v_tip_g , 'g', 'LineWidth', 1)
plot(inds_span, mean_v_tip_g - se_v_tip_g , 'g', 'LineWidth', 1)
plot(inds_span, mean_v_tip_r, 'r', 'LineWidth', 2)
plot(inds_span, mean_v_tip_r + se_v_tip_r , 'r', 'LineWidth', 1)
plot(inds_span, mean_v_tip_r - se_v_tip_r , 'r', 'LineWidth', 1)
plot(inds_span, mean_v_tip_l, 'b', 'LineWidth', 2)
plot(inds_span, mean_v_tip_l + se_v_tip_l , 'b', 'LineWidth', 1)
plot(inds_span, mean_v_tip_l - se_v_tip_l , 'b', 'LineWidth', 1)
xlabel(['Lick ' alignment  ' (ms)'])
ylabel('Velocity (mm/s)')
xlim([-400 400])
xline(0, 'k', 'LineWidth', 2);
title(['Number of licks: ' num2str(sum_g) '(G) & ' num2str(sum_l) '(L) & ' num2str(sum_r) '(R) | Total: ' num2str(sum_g+sum_r+sum_l)])

subplot(2,3,6)
hold on
plot(inds_span, mean_ang_tip_g, 'g', 'LineWidth', 2)
plot(inds_span, mean_ang_tip_g + se_ang_tip_g , 'g', 'LineWidth', 1)
plot(inds_span, mean_ang_tip_g - se_ang_tip_g , 'g', 'LineWidth', 1)
plot(inds_span, mean_ang_tip_r, 'r', 'LineWidth', 2)
plot(inds_span, mean_ang_tip_r + se_ang_tip_r , 'r', 'LineWidth', 1)
plot(inds_span, mean_ang_tip_r - se_ang_tip_r , 'r', 'LineWidth', 1)
plot(inds_span, mean_ang_tip_l, 'b', 'LineWidth', 2)
plot(inds_span, mean_ang_tip_l + se_ang_tip_l , 'b', 'LineWidth', 1)
plot(inds_span, mean_ang_tip_l - se_ang_tip_l , 'b', 'LineWidth', 1)
xlabel(['Lick ' alignment  ' (ms)'])
ylabel('Angle (deg)')
xlim([-400 400])
xline(0, 'k', 'LineWidth', 2);

%% BOUT SS RESPONSE
clearvars -except ALL_PCELL_COMPRESSED_DATA POPULATION

number_pCell =  length(ALL_PCELL_COMPRESSED_DATA);

%%% d_tip_str_bout %%%
for counter_pCell = 1 : number_pCell
    id(counter_pCell, 1) = string(ALL_PCELL_COMPRESSED_DATA(counter_pCell).id{1, 1});
    
    %%% d_tip_str_bout %%%
    d_tip_str_bout_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_str_bout);
    d_tip_str_bout_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_str_bout_r);
    d_tip_str_bout_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_str_bout_l);
    
    %%% SS_str_bout %%%
    SS_str_bout_all(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_str_bout * 100), 15, 'sgolay', 2);
    SS_str_bout_baseline_all(counter_pCell, 1) = nanmean(SS_str_bout_all(counter_pCell, 1:50));
    SS_str_bout_post_all(counter_pCell, 1) = nanmean(SS_str_bout_all(counter_pCell, 101:200));
    change_SS_str_bout_all(counter_pCell,:) = SS_str_bout_all(counter_pCell,:) - SS_str_bout_baseline_all(counter_pCell,1);
    mean_change_SS_str_bout_all(counter_pCell, 1) = nanmean(change_SS_str_bout_all(counter_pCell,101:200));
    se_change_SS_str_bout_all(counter_pCell, 1) = nanstd(change_SS_str_bout_all(counter_pCell,101:200))/sqrt(number_pCell);
    
    SS_str_bout_r(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_str_bout_r * 100), 15, 'sgolay', 2);
    SS_str_bout_baseline_r(counter_pCell, 1) = nanmean(SS_str_bout_r(counter_pCell, 1:50));
    SS_str_bout_post_r(counter_pCell, 1) = nanmean(SS_str_bout_r(counter_pCell, 101:200));
    change_SS_str_bout_r(counter_pCell,:) = SS_str_bout_r(counter_pCell,:) - SS_str_bout_baseline_r(counter_pCell,1);
    mean_change_SS_str_bout_r(counter_pCell, 1) = nanmean(change_SS_str_bout_r(counter_pCell,101:200));
    se_change_SS_str_bout_r(counter_pCell, 1) = nanstd(change_SS_str_bout_r(counter_pCell,101:200))/sqrt(number_pCell);
    
    SS_str_bout_l(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_str_bout_l * 100), 15, 'sgolay', 2);
    SS_str_bout_baseline_l(counter_pCell, 1) = nanmean(SS_str_bout_l(counter_pCell, 1:50));
    SS_str_bout_post_l(counter_pCell, 1) = nanmean(SS_str_bout_l(counter_pCell, 101:200));
    change_SS_str_bout_l(counter_pCell,:) = SS_str_bout_l(counter_pCell,:) - SS_str_bout_baseline_l(counter_pCell,1);
    mean_change_SS_str_bout_l(counter_pCell, 1) = nanmean(change_SS_str_bout_l(counter_pCell,101:200));
    se_change_SS_str_bout_l(counter_pCell, 1) = nanstd(change_SS_str_bout_l(counter_pCell,101:200))/sqrt(number_pCell);
    
end

%%% Mean of d_tip_str_bout %%%
mean_d_tip_str_bout_all = nanmean(d_tip_str_bout_all);
se_d_tip_str_bout_all = nanstd(d_tip_str_bout_all)/sqrt(number_pCell);
mean_d_tip_str_bout_r = nanmean(d_tip_str_bout_r);
se_d_tip_str_bout_r = nanstd(d_tip_str_bout_r)/sqrt(number_pCell);
mean_d_tip_str_bout_l = nanmean(d_tip_str_bout_l);
se_d_tip_str_bout_l = nanstd(d_tip_str_bout_l)/sqrt(number_pCell);

%%% Mean of SS_str_bout %%%
mean_SS_str_bout_all = nanmean(SS_str_bout_all);
se_SS_str_bout_all= nanstd(SS_str_bout_all)/sqrt(number_pCell);
change_mean_SS_str_bout_all = nanmean(change_SS_str_bout_all);
change_se_SS_str_bout_all= nanstd(change_SS_str_bout_all)/sqrt(number_pCell);
mean_SS_str_bout_r = nanmean(SS_str_bout_r);
se_SS_str_bout_r = nanstd(SS_str_bout_r)/sqrt(number_pCell);
change_mean_SS_str_bout_r = nanmean(change_SS_str_bout_r);
change_se_SS_str_bout_r= nanstd(change_SS_str_bout_r)/sqrt(number_pCell);
mean_SS_str_bout_l = nanmean(SS_str_bout_l);
se_SS_str_bout_l = nanstd(SS_str_bout_l)/sqrt(number_pCell);
change_mean_SS_str_bout_l = nanmean(change_SS_str_bout_l);
change_se_SS_str_bout_l = nanstd(change_SS_str_bout_l)/sqrt(number_pCell);

%%% Calculate diff in change based on direction
diff_mean_change_SS_str_bout_rl = mean_change_SS_str_bout_r - mean_change_SS_str_bout_l;

%%% bout labels %%%
label_large_increase_all       = mean_change_SS_str_bout_all >= 30;
label_medium_increase_all        = mean_change_SS_str_bout_all >= 15 & mean_change_SS_str_bout_all < 30 ;
label_small_increase_all      = mean_change_SS_str_bout_all >= 5 & mean_change_SS_str_bout_all < 15;
label_none_all                   = mean_change_SS_str_bout_all > -5 & mean_change_SS_str_bout_all < 5;
label_small_decrease_all         = mean_change_SS_str_bout_all > -15 & mean_change_SS_str_bout_all <= -5;
label_medium_decrease_all        = mean_change_SS_str_bout_all > -30 & mean_change_SS_str_bout_all <= -15;
label_large_decrease_all         = mean_change_SS_str_bout_all <= -30;

%%% bout labels R - L %%%
label_large_r       = diff_mean_change_SS_str_bout_rl >= 30;
label_medium_r        = diff_mean_change_SS_str_bout_rl >= 15 & diff_mean_change_SS_str_bout_rl < 30 ;
label_small_r      = diff_mean_change_SS_str_bout_rl >= 5 & diff_mean_change_SS_str_bout_rl < 15;
label_none_dir                   = diff_mean_change_SS_str_bout_rl > -5 & diff_mean_change_SS_str_bout_rl < 5;
label_small_l         = diff_mean_change_SS_str_bout_rl > -15 & diff_mean_change_SS_str_bout_rl <= -5;
label_medium_l        = diff_mean_change_SS_str_bout_rl > -30 & diff_mean_change_SS_str_bout_rl <= -15;
label_large_l         = diff_mean_change_SS_str_bout_rl <= -30;


%%% Mean d_tip_str_bout based on labels all
mean_d_tip_str_bout_large_increase_all = nanmean(d_tip_str_bout_all(label_large_increase_all,:));
se_d_tip_str_bout_large_increase_all = nanstd(d_tip_str_bout_all(label_large_increase_all,:))/sqrt(sum(label_large_increase_all));
mean_d_tip_str_bout_medium_increase_all = nanmean(d_tip_str_bout_all(label_medium_increase_all,:));
se_d_tip_str_bout_medium_increase_all = nanstd(d_tip_str_bout_all(label_medium_increase_all,:))/sqrt(sum(label_medium_increase_all));
mean_d_tip_str_bout_small_increase_all = nanmean(d_tip_str_bout_all(label_small_increase_all,:));
se_d_tip_str_bout_small_increase_all = nanstd(d_tip_str_bout_all(label_small_increase_all,:))/sqrt(sum(label_small_increase_all));
mean_d_tip_str_bout_none_all = nanmean(d_tip_str_bout_all(label_none_all,:));
se_d_tip_str_bout_none_all = nanstd(d_tip_str_bout_all(label_none_all,:))/sqrt(sum(label_none_all));
mean_d_tip_str_bout_large_decrease_all = nanmean(d_tip_str_bout_all(label_large_decrease_all,:));
se_d_tip_str_bout_large_decrease_all = nanstd(d_tip_str_bout_all(label_large_decrease_all,:))/sqrt(sum(label_large_decrease_all));
mean_d_tip_str_bout_medium_decrease_all = nanmean(d_tip_str_bout_all(label_medium_decrease_all,:));
se_d_tip_str_bout_medium_decrease_all = nanstd(d_tip_str_bout_all(label_medium_decrease_all,:))/sqrt(sum(label_medium_decrease_all));
mean_d_tip_str_bout_small_decrease_all = nanmean(d_tip_str_bout_all(label_small_decrease_all,:));
se_d_tip_str_bout_small_decrease_all = nanstd(d_tip_str_bout_all(label_small_decrease_all,:))/sqrt(sum(label_small_decrease_all));

%%% Mean d_tip_str_bout based on labels r
mean_d_tip_str_bout_large_r_r = d_tip_str_bout_r(label_large_r,:);
se_d_tip_str_bout_large_r_r = nanstd(d_tip_str_bout_r(label_large_r,:))/sqrt(sum(label_large_r));
mean_d_tip_str_bout_large_r_l = (d_tip_str_bout_l(label_large_r,:));
se_d_tip_str_bout_large_r_l = nanstd(d_tip_str_bout_l(label_large_r,:))/sqrt(sum(label_large_r));
mean_d_tip_str_bout_medium_r_r = nanmean(d_tip_str_bout_r(label_medium_r,:));
se_d_tip_str_bout_medium_r_r = nanstd(d_tip_str_bout_r(label_medium_r,:))/sqrt(sum(label_medium_r));
mean_d_tip_str_bout_medium_r_l = nanmean(d_tip_str_bout_l(label_medium_r,:));
se_d_tip_str_bout_medium_r_l = nanstd(d_tip_str_bout_l(label_medium_r,:))/sqrt(sum(label_medium_r));
mean_d_tip_str_bout_small_r_r = nanmean(d_tip_str_bout_r(label_small_r,:));
se_d_tip_str_bout_small_r_r = nanstd(d_tip_str_bout_r(label_small_r,:))/sqrt(sum(label_small_r));
mean_d_tip_str_bout_small_r_l = nanmean(d_tip_str_bout_l(label_small_r,:));
se_d_tip_str_bout_small_r_l = nanstd(d_tip_str_bout_l(label_small_r,:))/sqrt(sum(label_small_r));

%%% Mean d_tip_str_bout based on labels l
mean_d_tip_str_bout_large_l_r = nanmean(d_tip_str_bout_r(label_large_l,:));
se_d_tip_str_bout_large_l_r = nanstd(d_tip_str_bout_r(label_large_l,:))/sqrt(sum(label_large_l));
mean_d_tip_str_bout_large_l_l = nanmean(d_tip_str_bout_l(label_large_l,:));
se_d_tip_str_bout_large_l_l = nanstd(d_tip_str_bout_l(label_large_l,:))/sqrt(sum(label_large_l));
mean_d_tip_str_bout_medium_l_r = nanmean(d_tip_str_bout_r(label_medium_l,:));
se_d_tip_str_bout_medium_l_r = nanstd(d_tip_str_bout_r(label_medium_l,:))/sqrt(sum(label_medium_l));
mean_d_tip_str_bout_medium_l_l = nanmean(d_tip_str_bout_l(label_medium_l,:));
se_d_tip_str_bout_medium_l_l = nanstd(d_tip_str_bout_l(label_medium_l,:))/sqrt(sum(label_medium_l));
mean_d_tip_str_bout_small_l_r = nanmean(d_tip_str_bout_r(label_small_l,:));
se_d_tip_str_bout_small_l_r = nanstd(d_tip_str_bout_r(label_small_l,:))/sqrt(sum(label_small_l));
mean_d_tip_str_bout_small_l_l = nanmean(d_tip_str_bout_l(label_small_l,:));
se_d_tip_str_bout_small_l_l = nanstd(d_tip_str_bout_l(label_small_l,:))/sqrt(sum(label_small_l));

%%% Mean d_tip_str_bout based on labels none
mean_d_tip_str_bout_none_r = nanmean(d_tip_str_bout_r(label_none_dir,:));
se_d_tip_str_bout_none_r = nanstd(d_tip_str_bout_r(label_none_dir,:))/sqrt(sum(label_none_dir));
mean_d_tip_str_bout_none_l = nanmean(d_tip_str_bout_l(label_none_dir,:));
se_d_tip_str_bout_none_l = nanstd(d_tip_str_bout_l(label_none_dir,:))/sqrt(sum(label_none_dir));

%%% Mean change_SS_str_bout based on labels all
mean_change_SS_str_bout_large_increase_all = nanmean(change_SS_str_bout_all(label_large_increase_all,:));
se_change_SS_str_bout_large_increase_all = nanstd(change_SS_str_bout_all(label_large_increase_all,:))/sqrt(sum(label_large_increase_all));
mean_change_SS_str_bout_medium_increase_all = nanmean(change_SS_str_bout_all(label_medium_increase_all,:));
se_change_SS_str_bout_medium_increase_all = nanstd(change_SS_str_bout_all(label_medium_increase_all,:))/sqrt(sum(label_medium_increase_all));
mean_change_SS_str_bout_small_increase_all = nanmean(change_SS_str_bout_all(label_small_increase_all,:));
se_change_SS_str_bout_small_increase_all = nanstd(change_SS_str_bout_all(label_small_increase_all,:))/sqrt(sum(label_small_increase_all));
mean_change_SS_str_bout_none_all = nanmean(change_SS_str_bout_all(label_none_all,:));
se_change_SS_str_bout_none_all = nanstd(change_SS_str_bout_all(label_none_all,:))/sqrt(sum(label_none_all));
mean_change_SS_str_bout_large_decrease_all = nanmean(change_SS_str_bout_all(label_large_decrease_all,:));
se_change_SS_str_bout_large_decrease_all = nanstd(change_SS_str_bout_all(label_large_decrease_all,:))/sqrt(sum(label_large_decrease_all));
mean_change_SS_str_bout_medium_decrease_all = nanmean(change_SS_str_bout_all(label_medium_decrease_all,:));
se_change_SS_str_bout_medium_decrease_all = nanstd(change_SS_str_bout_all(label_medium_decrease_all,:))/sqrt(sum(label_medium_decrease_all));
mean_change_SS_str_bout_small_decrease_all = nanmean(change_SS_str_bout_all(label_small_decrease_all,:));
se_change_SS_str_bout_small_decrease_all = nanstd(change_SS_str_bout_all(label_small_decrease_all,:))/sqrt(sum(label_small_decrease_all));

%%% Mean d_tip_str_bout based on labels r
mean_change_SS_str_bout_large_r_r = (change_SS_str_bout_r(label_large_r,:));
se_change_SS_str_bout_large_r_r = nanstd(change_SS_str_bout_r(label_large_r,:))/sqrt(sum(label_large_r));
mean_change_SS_str_bout_large_r_l = (change_SS_str_bout_l(label_large_r,:));
se_change_SS_str_bout_large_r_l = nanstd(change_SS_str_bout_l(label_large_r,:))/sqrt(sum(label_large_r));
mean_change_SS_str_bout_medium_r_r = nanmean(change_SS_str_bout_r(label_medium_r,:));
se_change_SS_str_bout_medium_r_r = nanstd(change_SS_str_bout_r(label_medium_r,:))/sqrt(sum(label_medium_r));
mean_change_SS_str_bout_medium_r_l = nanmean(change_SS_str_bout_l(label_medium_r,:));
se_change_SS_str_bout_medium_r_l = nanstd(change_SS_str_bout_l(label_medium_r,:))/sqrt(sum(label_medium_r));
mean_change_SS_str_bout_small_r_r = nanmean(change_SS_str_bout_r(label_small_r,:));
se_change_SS_str_bout_small_r_r = nanstd(change_SS_str_bout_r(label_small_r,:))/sqrt(sum(label_small_r));
mean_change_SS_str_bout_small_r_l = nanmean(change_SS_str_bout_l(label_small_r,:));
se_change_SS_str_bout_small_r_l = nanstd(change_SS_str_bout_l(label_small_r,:))/sqrt(sum(label_small_r));

%%% Mean change_SS_str_bout based on labels l
mean_change_SS_str_bout_large_l_r = nanmean(change_SS_str_bout_r(label_large_l,:));
se_change_SS_str_bout_large_l_r = nanstd(change_SS_str_bout_r(label_large_l,:))/sqrt(sum(label_large_l));
mean_change_SS_str_bout_large_l_l = nanmean(change_SS_str_bout_l(label_large_l,:));
se_change_SS_str_bout_large_l_l = nanstd(change_SS_str_bout_l(label_large_l,:))/sqrt(sum(label_large_l));
mean_change_SS_str_bout_medium_l_r = nanmean(change_SS_str_bout_r(label_medium_l,:));
se_change_SS_str_bout_medium_l_r = nanstd(change_SS_str_bout_r(label_medium_l,:))/sqrt(sum(label_medium_l));
mean_change_SS_str_bout_medium_l_l = nanmean(change_SS_str_bout_l(label_medium_l,:));
se_change_SS_str_bout_medium_l_l = nanstd(change_SS_str_bout_l(label_medium_l,:))/sqrt(sum(label_medium_l));
mean_change_SS_str_bout_small_l_r = nanmean(change_SS_str_bout_r(label_small_l,:));
se_change_SS_str_bout_small_l_r = nanstd(change_SS_str_bout_r(label_small_l,:))/sqrt(sum(label_small_l));
mean_change_SS_str_bout_small_l_l = nanmean(change_SS_str_bout_l(label_small_l,:));
se_change_SS_str_bout_small_l_l = nanstd(change_SS_str_bout_l(label_small_l,:))/sqrt(sum(label_small_l));

%%% Mean change_SS_str_bout based on labels none
mean_change_SS_str_bout_none_r = nanmean(change_SS_str_bout_r(label_none_dir,:));
se_change_SS_str_bout_none_r = nanstd(change_SS_str_bout_r(label_none_dir,:))/sqrt(sum(label_none_dir));
mean_change_SS_str_bout_none_l = nanmean(change_SS_str_bout_l(label_none_dir,:));
se_change_SS_str_bout_none_l = nanstd(change_SS_str_bout_l(label_none_dir,:))/sqrt(sum(label_none_dir));

%%% inds_span_bout %%%
inds_span_bout = -0.99 : 0.01 : 1;

%%% xaxis for bar plot %%%
xaxis_pCell = 1:62;

%%% Save to POPULATION structure %%%
POPULATION.bout.label_large_increase_all = label_large_increase_all;
POPULATION.bout.label_medium_increase_all = label_medium_increase_all;
POPULATION.bout.label_small_increase_all = label_small_increase_all;
POPULATION.bout.label_none_all = label_none_all;
POPULATION.bout.label_small_decrease_all = label_small_decrease_all;
POPULATION.bout.label_medium_decrease_all = label_medium_decrease_all;
POPULATION.bout.label_large_decrease_all = label_large_decrease_all;
POPULATION.bout.label_large_r = label_large_r;
POPULATION.bout.label_medium_r = label_medium_r;
POPULATION.bout.label_small_r = label_small_r;
POPULATION.bout.label_none_dir = label_none_dir;
POPULATION.bout.label_small_l = label_small_l;
POPULATION.bout.label_medium_l = label_medium_l;
POPULATION.bout.label_large_l = label_large_l;
%% Raw
num_pCell = [6];
%     num_pCell = [6 2 1 31 35];
color = [1 0 0];
figure
sgtitle(['pCells: ' num2str(number_pCell)])
subplot(4,3,1)
hold on
for counter_pCell = num_pCell
    plot(inds_span_bout, d_tip_str_bout_all(counter_pCell, :),'Color',[0.7 0.7 0.7] , 'LineWidth' , 2 )
end
%     plot(inds_span_bout, mean_d_tip_str_bout_all , 'Color', [0.7 0.7 0.7], 'LineWidth' , 2 )
%     plot(inds_span_bout, mean_d_tip_str_bout_all + se_d_tip_str_bout_all , 'Color', [0.7 0.7 0.7], 'LineWidth' , 1 )
%     plot(inds_span_bout, mean_d_tip_str_bout_all - se_d_tip_str_bout_all , 'Color', [0.7 0.7 0.7], 'LineWidth' ,1 )
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-1 1])
ylim([-1 18])
xlabel('Bout start (s)')
ylabel('Displacement (mm)')
title(['Displacement | All | pCell: ' num2str(counter_pCell)])

subplot(4,3,2)
hold on
for counter_pCell = num_pCell
    plot(inds_span_bout, d_tip_str_bout_l(counter_pCell, :),'Color',[0.7 0.7 0.7] , 'LineWidth' , 2 )
end
%     plot(inds_span_bout, mean_d_tip_str_bout_l , 'Color', [0.7 0.7 0.7], 'LineWidth' , 2 )
%     plot(inds_span_bout, mean_d_tip_str_bout_l + se_d_tip_str_bout_l , 'Color', [0.7 0.7 0.7], 'LineWidth' , 1 )
%     plot(inds_span_bout, mean_d_tip_str_bout_l - se_d_tip_str_bout_l , 'Color', [0.7 0.7 0.7], 'LineWidth' ,1 )
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-1 1])
ylim([-1 18])
xlabel('Bout start (s)')
ylabel('Displacement (mm)')
title(['Displacement | L | pCell: ' num2str(counter_pCell)])

subplot(4,3,3)
hold on
for counter_pCell = num_pCell
    plot(inds_span_bout, d_tip_str_bout_r(counter_pCell, :),'Color',[0.7 0.7 0.7] , 'LineWidth' , 2 )
end
%     plot(inds_span_bout, mean_d_tip_str_bout_r , 'Color',[0.7 0.7 0.7], 'LineWidth' , 2 )
%     plot(inds_span_bout, mean_d_tip_str_bout_r + se_d_tip_str_bout_r , 'Color', [0.7 0.7 0.7], 'LineWidth' , 1 )
%     plot(inds_span_bout, mean_d_tip_str_bout_r - se_d_tip_str_bout_r , 'Color', [0.7 0.7 0.7], 'LineWidth' ,1 )
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-1 1])
ylim([-1 18])
xlabel('Bout start (s)')
ylabel('Displacement (mm)')
title(['Displacement | R | pCell: ' num2str(counter_pCell)])

subplot(4,3,4)
hold on
for counter_pCell = num_pCell
    plot(inds_span_bout, smooth(change_SS_str_bout_all(counter_pCell,:),10,'sgolay', 2), 'Color', [1 0 0] , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-1 1])
ylim([-30 200])
xlabel('Bout start (s)')
ylabel('Change in SS firing (spks/s)')
title(['SS firing | All | pCell: ' num2str(counter_pCell)])

subplot(4,3,5)
hold on
for counter_pCell = num_pCell
    plot(inds_span_bout, smooth(change_SS_str_bout_l(counter_pCell,:),10,'sgolay', 2),'Color', color , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-1 1])
ylim([-30 200])
xlabel('Bout start (s)')
ylabel('Change in SS firing (spks/s)')
title(['SS firing | L | pCell: ' num2str(counter_pCell)])

subplot(4,3,6)
hold on
for counter_pCell = num_pCell
    plot(inds_span_bout, smooth(change_SS_str_bout_r(counter_pCell,:),10,'sgolay', 2),'Color', color , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-1 1])
ylim([-30 200])
xlabel('Bout start (s)')
ylabel('Change in SS firing (spks/s)')
title(['SS firing | R | pCell: ' num2str(counter_pCell)])

subplot(4,3,7)
hold on
for counter_pCell = 1 : number_pCell
    plot(SS_str_bout_baseline_all(counter_pCell), SS_str_bout_post_all(counter_pCell), 'ok')
end
plot(SS_str_bout_baseline_all(num_pCell), SS_str_bout_post_all(num_pCell), 'or')
vert = refline(1,0);
vert.Color = 'k';
xlabel('Mean SS firing before bout (spks/s)')
ylabel('Mean SS firing during bout (spks/s)')
xlim([0 250])
ylim([0 250])
title('Bout effect on SS firing | All licks')

subplot(4,3,8)
hold on
for counter_pCell = 1 : number_pCell
    plot(SS_str_bout_baseline_l(counter_pCell), SS_str_bout_post_l(counter_pCell), 'ok')
end
plot(SS_str_bout_baseline_l(num_pCell), SS_str_bout_post_l(num_pCell), 'or')
vert = refline(1,0);
vert.Color = 'k';
xlabel('Mean SS firing before bout (spks/s)')
ylabel('Mean SS firing during bout (spks/s)')
xlim([0 250])
ylim([0 250])
title('Bout effect on SS firing | L licks')

subplot(4,3,9)
hold on
for counter_pCell = 1 : number_pCell
    plot(SS_str_bout_baseline_r(counter_pCell), SS_str_bout_post_r(counter_pCell), 'ok')
end
plot(SS_str_bout_baseline_r(num_pCell), SS_str_bout_post_r(num_pCell), 'or')
vert = refline(1,0);
vert.Color = 'k';
xlabel('Mean SS firing before bout (spks/s)')
ylabel('Mean SS firing during bout (spks/s)')
xlim([0 250])
ylim([0 250])
title('Bout effect on SS firing | R licks')

subplot(4,3,10)
hold on;
bar(mean_change_SS_str_bout_all, 'k')
bar(xaxis_pCell(num_pCell), mean_change_SS_str_bout_all(num_pCell), 'FaceColor', color)
% plot(mean_change_SS_str_bout_all + se_change_SS_str_bout_all, '*r')
% plot(mean_change_SS_str_bout_all - se_change_SS_str_bout_all, '*r')
ylim([-20 115])
xlabel('pCell')
ylabel('Change in SS firing (spks/s)')

subplot(4,3,11)
hold on;
bar(mean_change_SS_str_bout_l, 'k')
bar(xaxis_pCell(num_pCell), mean_change_SS_str_bout_l(num_pCell), 'FaceColor', color)
% plot(mean_change_SS_str_bout_l + se_change_SS_str_bout_l, '*r')
% plot(mean_change_SS_str_bout_l - se_change_SS_str_bout_l, '*r')
ylim([-20 115])
xlabel('pCell')
ylabel('Change in SS firing (spks/s)')

subplot(4,3,12)
hold on;
bar(mean_change_SS_str_bout_r, 'k')
bar(xaxis_pCell(num_pCell), mean_change_SS_str_bout_r(num_pCell), 'FaceColor', color)
% plot(mean_change_SS_str_bout_r + se_change_SS_str_bout_r, '*r')
% plot(mean_change_SS_str_bout_r - se_change_SS_str_bout_r, '*r')
ylim([-20 115])
xlabel('pCell')
ylabel('Change in SS firing (spks/s)')
%% Labeled properties
%%% PROPERTIES %%%
figure
sgtitle(['pCells: ' num2str(number_pCell)])

subplot(2,5,[1 2])
hold on;
bar(xaxis_pCell(label_large_increase_all ), mean_change_SS_str_bout_all(label_large_increase_all ), 'r', 'Barwidth', 1)
bar(xaxis_pCell(label_medium_increase_all ), mean_change_SS_str_bout_all(label_medium_increase_all ), 'b', 'Barwidth', 1)
bar(xaxis_pCell(label_small_increase_all ), mean_change_SS_str_bout_all(label_small_increase_all ), 'g', 'Barwidth', 1)
bar(xaxis_pCell(label_none_all), mean_change_SS_str_bout_all(label_none_all ), 'k', 'Barwidth', 1)
bar(xaxis_pCell(label_small_decrease_all ), mean_change_SS_str_bout_all(label_small_decrease_all ), 'm', 'Barwidth', 0.20)
yline(30,'--r', 'LineWidth', 2);
yline(15,'--b', 'LineWidth', 2);
yline(5,'--g', 'LineWidth', 2);
yline(-5,'--m', 'LineWidth', 2);
ylim([-20 115])
xlabel('pCell')
ylabel('Change in SS firing (spks/s)')
title('Effect of bout on SS firing | All bouts')

subplot(2,5,6)
hold on;
plot(SS_str_bout_baseline_all(label_large_increase_all), SS_str_bout_post_all(label_large_increase_all), 'or')
plot(SS_str_bout_baseline_all(label_medium_increase_all), SS_str_bout_post_all(label_medium_increase_all), 'ob')
plot(SS_str_bout_baseline_all(label_small_increase_all), SS_str_bout_post_all(label_small_increase_all), 'og')
plot(SS_str_bout_baseline_all(label_none_all), SS_str_bout_post_all(label_none_all), 'ok')
plot(SS_str_bout_baseline_all(label_small_decrease_all), SS_str_bout_post_all(label_small_decrease_all), 'om')
vert = refline(1,0);
vert.Color = 'k';
xlabel('Mean SS firing before bout (spks/s)')
ylabel('Mean SS firing during bout (spks/s)')
xlim([0 250])
ylim([0 250])

pie_1 = subplot(2,5,7);
pie_values = [sum(label_large_increase_all) sum(label_medium_increase_all) sum(label_small_increase_all) sum(label_none_all) sum(label_small_decrease_all)];
pie_labels = {num2str(sum(label_large_increase_all)) ...
    num2str(sum(label_medium_increase_all)) ...
    num2str(sum(label_small_increase_all)) ...
    num2str(sum(label_none_all)) ...
    num2str(sum(label_small_decrease_all))};
explode_1 = [1 1 1 1 1];
pie(pie_values, explode_1, pie_labels)
colormap(pie_1, [1 0 0; 0 0 1; 0 1 0; 0 0 0; 1 0 1]);
%     title('Modulation classification breakdown')

subplot(2,5,[3 4 8 9])
hold on;
barh(xaxis_pCell(label_large_r), diff_mean_change_SS_str_bout_rl(label_large_r), 'r', 'BaseValue', 0, 'Barwidth', 1)
barh(xaxis_pCell(label_medium_r), diff_mean_change_SS_str_bout_rl(label_medium_r), 'b', 'BaseValue', 0, 'Barwidth', 0.85)
barh(xaxis_pCell(label_small_r), diff_mean_change_SS_str_bout_rl(label_small_r), 'g', 'BaseValue', 0, 'Barwidth', 0.5)
barh(xaxis_pCell(label_none_dir), diff_mean_change_SS_str_bout_rl(label_none_dir), 'k', 'BaseValue', 0, 'Barwidth', 0.9)
barh(xaxis_pCell(label_small_l), diff_mean_change_SS_str_bout_rl(label_small_l), 'g', 'BaseValue', 0, 'Barwidth', 1)
barh(xaxis_pCell(label_medium_l), diff_mean_change_SS_str_bout_rl(label_medium_l), 'b', 'BaseValue', 0, 'Barwidth', 0.15)
barh(xaxis_pCell(label_large_l), diff_mean_change_SS_str_bout_rl(label_large_l), 'r', 'BaseValue', 0, 'Barwidth', 0.03)
xline(30,'--r', 'LineWidth', 2);
xline(-30,'--r', 'LineWidth', 2);
xline(15,'--b', 'LineWidth', 2);
xline(-15,'--b', 'LineWidth', 2);
xline(5,'--g', 'LineWidth', 2);
xline(-5,'--g', 'LineWidth', 2);
xlabel('Difference in change of SS firing (spks/s)')
ylabel('pCell')
xlim([-70 70])
ylim([-3 65])
title('Effect of bout direction on SS firing | R - L bouts')

subplot(2,5,5)
hold on;
plot(mean_change_SS_str_bout_r(label_large_r), mean_change_SS_str_bout_l(label_large_r), 'or')
plot(mean_change_SS_str_bout_r(label_medium_r), mean_change_SS_str_bout_l(label_medium_r), 'ob')
plot(mean_change_SS_str_bout_r(label_small_r), mean_change_SS_str_bout_l(label_small_r), 'og')
plot(mean_change_SS_str_bout_r(label_none_dir), mean_change_SS_str_bout_l(label_none_dir), 'ok')
plot(mean_change_SS_str_bout_r(label_small_l), mean_change_SS_str_bout_l(label_small_l), 'og')
plot(mean_change_SS_str_bout_r(label_medium_l), mean_change_SS_str_bout_l(label_medium_l), 'ob')
plot(mean_change_SS_str_bout_r(label_large_l), mean_change_SS_str_bout_l(label_large_l), 'or')
vert = refline(1,0);
vert.Color = 'k';
xlabel('Mean change in SS firing R bout (spks/s)')
ylabel('Mean change in SS firing L bout (spks/s)')
xlim([-30 125])
ylim([-30 125])

pie_2 = subplot(2,5,10);
pie_values = [sum(label_large_r) sum(label_medium_r) sum(label_small_r) sum(label_none_dir) ...
    sum(label_small_l) sum(label_medium_l) sum(label_large_l)];
pie_labels = {['R: ' num2str(sum(label_large_r))] ...
    ['R: ' num2str(sum(label_medium_r))] ...
    ['R: '  num2str(sum(label_small_r))] ...
    ['None: ' num2str(sum(label_none_dir))] ...
    ['L: ' num2str(sum(label_small_l))] ...
    ['L: ' num2str(sum(label_medium_l))] ...
    ['L: ' num2str(sum(label_large_l))]};
explode_2 = [1 1 1 1 1 1 1];
pie(pie_values, explode_2, pie_labels)
colormap(pie_2, [1 0 0; 0 0 1; 0 1 0; 0 0 0; 0 1 0; 0 0 1; 1 0 0]);
%     title('Modulation classification breakdown - directional')
%% Labeled traces
%%% TRACES %%%
figure
sgtitle(['Bout | pCells: ' num2str(number_pCell)])

subplot(2,5,1)
hold on
plot(inds_span_bout, mean_d_tip_str_bout_large_increase_all(:), 'r', 'LineWidth', 2)
%     plot(inds_span_bout, mean_d_tip_str_bout_large_increase_all + se_d_tip_str_bout_large_increase_all, 'r', 'LineWidth', 1)
%     plot(inds_span_bout, mean_d_tip_str_bout_large_increase_all - se_d_tip_str_bout_large_increase_all, 'r', 'LineWidth', 1)
plot(inds_span_bout, mean_d_tip_str_bout_medium_increase_all(:), 'b', 'LineWidth', 2)
%     plot(inds_span_bout, mean_d_tip_str_bout_medium_increase_all + se_d_tip_str_bout_medium_increase_all, 'b', 'LineWidth', 1)
%     plot(inds_span_bout, mean_d_tip_str_bout_medium_increase_all - se_d_tip_str_bout_medium_increase_all, 'b', 'LineWidth', 1)
plot(inds_span_bout, mean_d_tip_str_bout_small_increase_all(:), 'g', 'LineWidth', 2)
%     plot(inds_span_bout, mean_d_tip_str_bout_small_increase_all + se_d_tip_str_bout_small_increase_all, 'g', 'LineWidth', 1)
%     plot(inds_span_bout, mean_d_tip_str_bout_small_increase_all - se_d_tip_str_bout_small_increase_all, 'g', 'LineWidth', 1)
plot(inds_span_bout, mean_d_tip_str_bout_none_all(:), 'k', 'LineWidth', 2)
%     plot(inds_span_bout, mean_d_tip_str_bout_none_all + se_d_tip_str_bout_none_all, 'k', 'LineWidth', 1)
%     plot(inds_span_bout, mean_d_tip_str_bout_none_all - se_d_tip_str_bout_none_all, 'k', 'LineWidth', 1)
plot(inds_span_bout, mean_d_tip_str_bout_small_decrease_all(:), 'm', 'LineWidth', 2)
%     plot(inds_span_bout, mean_d_tip_str_bout_small_decrease_all + se_d_tip_str_bout_small_decrease_all, 'm', 'LineWidth', 1)
%     plot(inds_span_bout, mean_d_tip_str_bout_small_decrease_all - se_d_tip_str_bout_small_decrease_all, 'm', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-1 1])
ylabel('Displacement (mm)')
xlabel('Bout start (s)')
title(['Bout trace | All bouts'])

subplot(2,5,2)
hold on
plot(inds_span_bout, mean_d_tip_str_bout_large_l_l(:), 'r', 'LineWidth', 2)
plot(inds_span_bout, mean_d_tip_str_bout_medium_l_l(:), 'b', 'LineWidth', 2)
plot(inds_span_bout, mean_d_tip_str_bout_small_l_l(:), 'g', 'LineWidth', 2)
plot(inds_span_bout, mean_d_tip_str_bout_none_l(:), 'k', 'LineWidth', 2)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-1 1])
ylabel('Displacement (mm)')
xlabel('Bout start (s)')
title(['Bout trace | L-tuned, L bouts'])

subplot(2,5,3)
hold on
plot(inds_span_bout, mean_d_tip_str_bout_large_l_r(:), 'r', 'LineWidth', 2)
plot(inds_span_bout, mean_d_tip_str_bout_medium_l_r(:), 'b', 'LineWidth', 2)
plot(inds_span_bout, mean_d_tip_str_bout_small_l_r(:), 'g', 'LineWidth', 2)
plot(inds_span_bout, mean_d_tip_str_bout_none_r(:), 'k', 'LineWidth', 2)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-1 1])
ylabel('Displacement (mm)')
xlabel('Bout start (s)')
title(['Bout trace | L-tuned, R bouts'])

subplot(2,5,4)
hold on
plot(inds_span_bout, mean_d_tip_str_bout_large_r_l(:), 'r', 'LineWidth', 2)
plot(inds_span_bout, mean_d_tip_str_bout_medium_r_l(:), 'b', 'LineWidth', 2)
plot(inds_span_bout, mean_d_tip_str_bout_small_r_l(:), 'g', 'LineWidth', 2)
plot(inds_span_bout, mean_d_tip_str_bout_none_r(:), 'k', 'LineWidth', 2)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-1 1])
ylabel('Displacement (mm)')
xlabel('Bout start (s)')
title(['Bout trace | R-tuned, L bouts'])

subplot(2,5,5)
hold on
plot(inds_span_bout, mean_d_tip_str_bout_large_r_r(:), 'r', 'LineWidth', 2)
plot(inds_span_bout, mean_d_tip_str_bout_medium_r_r(:), 'b', 'LineWidth', 2)
plot(inds_span_bout, mean_d_tip_str_bout_small_r_r(:), 'g', 'LineWidth', 2)
plot(inds_span_bout, mean_d_tip_str_bout_none_r(:), 'k', 'LineWidth', 2)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-1 1])
ylabel('Displacement (mm)')
xlabel('Bout start (s)')
title(['Bout trace | R-tuned, R bouts'])

subplot(2,5,6)
hold on
plot(inds_span_bout, mean_change_SS_str_bout_large_increase_all(:), 'r', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_large_increase_all + se_change_SS_str_bout_large_increase_all, 'r', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_large_increase_all - se_change_SS_str_bout_large_increase_all, 'r', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_medium_increase_all(:), 'b', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_medium_increase_all + se_change_SS_str_bout_medium_increase_all, 'b', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_medium_increase_all - se_change_SS_str_bout_medium_increase_all, 'b', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_small_increase_all(:), 'g', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_small_increase_all + se_change_SS_str_bout_small_increase_all, 'g', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_small_increase_all - se_change_SS_str_bout_small_increase_all, 'g', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_none_all(:), 'k', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_none_all + se_change_SS_str_bout_none_all, 'k', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_none_all - se_change_SS_str_bout_none_all, 'k', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_small_decrease_all(:), 'm', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_small_decrease_all + se_change_SS_str_bout_small_decrease_all, 'm', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_small_decrease_all - se_change_SS_str_bout_small_decrease_all, 'm', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(30,'--r', 'LineWidth', 1);
yline(-30,'--r', 'LineWidth', 1);
yline(15,'--b', 'LineWidth', 1);
yline(-15,'--b', 'LineWidth', 1);
yline(5,'--g', 'LineWidth', 1);
yline(-5,'--m', 'LineWidth', 1);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
xlim([-1 1])
ylim([-30 150])
xlabel('Bout start (s)')
ylabel('Change in SS firing (spks/s)')
title(['Effect of bout on SS firing | All bouts'])

subplot(2,5,7)
hold on
plot(inds_span_bout, mean_change_SS_str_bout_large_l_l(:), 'r', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_large_l_l(:) + se_change_SS_str_bout_large_l_l(:), 'r', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_large_l_l(:) - se_change_SS_str_bout_large_l_l(:), 'r', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_medium_l_l(:), 'b', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_medium_l_l(:) + se_change_SS_str_bout_medium_l_l(:), 'b', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_medium_l_l(:) - se_change_SS_str_bout_medium_l_l(:), 'b', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_small_l_l(:), 'g', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_small_l_l(:) + se_change_SS_str_bout_small_l_l(:), 'g', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_small_l_l(:) - se_change_SS_str_bout_small_l_l(:), 'g', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_none_l(:), 'k', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_none_l(:) + se_change_SS_str_bout_none_l(:) , 'k', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_none_l(:) - se_change_SS_str_bout_none_l(:) , 'k', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
xlim([-1 1])
ylim([-30 150])
xlabel('Bout start (s)')
ylabel('Change in SS firing (spks/s)')
title(['Effect of bout on SS firing | L-tuned L bouts'])

subplot(2,5,8)
hold on
plot(inds_span_bout, mean_change_SS_str_bout_large_l_r(:), 'r', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_large_l_r(:) + se_change_SS_str_bout_large_l_r(:), 'r', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_large_l_r(:) - se_change_SS_str_bout_large_l_r(:), 'r', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_medium_l_r(:), 'b', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_medium_l_r(:) + se_change_SS_str_bout_medium_l_r(:), 'b', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_medium_l_r(:) - se_change_SS_str_bout_medium_l_r(:), 'b', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_small_l_r(:), 'g', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_small_l_r(:) + se_change_SS_str_bout_small_l_r(:), 'g', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_small_l_r(:) - se_change_SS_str_bout_small_l_r(:), 'g', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_none_r(:), 'k', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_none_r(:) + se_change_SS_str_bout_none_r(:) , 'k', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_none_r(:) - se_change_SS_str_bout_none_r(:) , 'k', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
xlim([-1 1])
ylim([-30 150])
xlabel('Bout start (s)')
ylabel('Change in SS firing (spks/s)')
title(['Effect of bout on SS firing | L-tuned R bouts'])

subplot(2,5,9)
hold on
plot(inds_span_bout, mean_change_SS_str_bout_large_r_l(:), 'r', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_large_r_l(:) + se_change_SS_str_bout_large_r_l(:), 'r', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_large_r_l(:) - se_change_SS_str_bout_large_r_l(:), 'r', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_medium_r_l(:), 'b', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_medium_r_l(:) + se_change_SS_str_bout_medium_r_l(:), 'b', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_medium_r_l(:) - se_change_SS_str_bout_medium_r_l(:), 'b', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_small_r_l(:), 'g', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_small_r_l(:) + se_change_SS_str_bout_small_r_l(:), 'g', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_small_r_l(:) - se_change_SS_str_bout_small_r_l(:), 'g', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_none_l(:), 'k', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_none_l(:) + se_change_SS_str_bout_none_l(:) , 'k', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_none_l(:) - se_change_SS_str_bout_none_l(:) , 'k', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
xlim([-1 1])
ylim([-30 150])
xlabel('Bout start (s)')
ylabel('Change in SS firing (spks/s)')
title(['Effect of bout on SS firing | R-tuned L bouts'])

subplot(2,5,10)
hold on
plot(inds_span_bout, mean_change_SS_str_bout_large_r_r(:), 'r', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_large_r_r(:) + se_change_SS_str_bout_large_r_r(:), 'r', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_large_r_r(:) - se_change_SS_str_bout_large_r_r(:), 'r', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_medium_r_r(:), 'b', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_medium_r_r(:) + se_change_SS_str_bout_medium_r_r(:), 'b', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_medium_r_r(:) - se_change_SS_str_bout_medium_r_r(:), 'b', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_small_r_r(:), 'g', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_small_r_r(:) + se_change_SS_str_bout_small_r_r(:), 'g', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_small_r_r(:) - se_change_SS_str_bout_small_r_r(:), 'g', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_none_r(:), 'k', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_none_r(:) + se_change_SS_str_bout_none_r(:) , 'k', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_none_r(:) - se_change_SS_str_bout_none_r(:) , 'k', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
xlim([-1 1])
ylim([-30 150])
xlabel('Bout start (s)')
ylabel('Change in SS firing (spks/s)')
title(['Effect of bout on SS firing | R-tuned R bouts'])

%% LICK SS RESPONSE
clearvars -except ALL_PCELL_COMPRESSED_DATA POPULATION

alignment = 'dmax';
% alignment = 'onset';

number_pCell =  length(ALL_PCELL_COMPRESSED_DATA);

for counter_pCell = 1 : number_pCell
    id(counter_pCell, 1) = string(ALL_PCELL_COMPRESSED_DATA(counter_pCell).id{1, 1});
    %%% Initialize based on alignment event %%%
    if contains(alignment,'onset') == 1
        %%% d_tip %%%
        d_tip_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_all);
        d_tip_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_grooming);
        d_tip_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_r);
        d_tip_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_l);
        
        %%% SS_str_bout %%%
        SS_str_bout_all(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_str_bout * 100), 15, 'sgolay', 2);
        SS_str_bout_baseline_all(counter_pCell, :) = nanmean(SS_str_bout_all(counter_pCell, 1:50));
        SS_str_bout_r(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_str_bout_r * 100), 15, 'sgolay', 2);
        SS_str_bout_baseline_r(counter_pCell, :) = nanmean(SS_str_bout_r(counter_pCell, 1:50));
        SS_str_bout_l(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_str_bout_l * 100), 15, 'sgolay', 2);
        SS_str_bout_baseline_l(counter_pCell, :) = nanmean(SS_str_bout_l(counter_pCell, 1:50));
        
        %%% SS_lick %%%
        SS_lick_all(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_all * 100), 15, 'sgolay', 2);
        change_SS_lick_all(counter_pCell,:) = SS_lick_all(counter_pCell,:) - SS_str_bout_baseline_all(counter_pCell, 1);
        mean_window_SS_lick_all(counter_pCell,1) = nanmean(SS_lick_all(counter_pCell, :));
        se_window_SS_lick_all(counter_pCell,1) = nanstd(SS_lick_all(counter_pCell, :))/sqrt(number_pCell);
        mean_change_window_SS_lick_all(counter_pCell,1) = nanmean(change_SS_lick_all(counter_pCell, :));
        se_change_window_SS_lick_all(counter_pCell,1) = nanstd(change_SS_lick_all(counter_pCell, :))/sqrt(number_pCell);
        
        SS_lick_r(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_r * 100), 15, 'sgolay', 2);
        change_SS_lick_r(counter_pCell,:) = SS_lick_r(counter_pCell,:) - SS_str_bout_baseline_r(counter_pCell, 1);
        mean_window_SS_lick_r(counter_pCell,1) = nanmean(SS_lick_r(counter_pCell, :));
        se_window_SS_lick_r(counter_pCell,1) = nanstd(SS_lick_r(counter_pCell, :))/sqrt(number_pCell);
        mean_change_window_SS_lick_r(counter_pCell,1) = nanmean(change_SS_lick_r(counter_pCell, :));
        se_change_window_SS_lick_r(counter_pCell,1) = nanstd(change_SS_lick_r(counter_pCell, :))/sqrt(number_pCell);
        
        SS_lick_g(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_grooming * 100), 15, 'sgolay', 2);
        change_SS_lick_g(counter_pCell,:) = SS_lick_g(counter_pCell,:) - SS_str_bout_baseline_all(counter_pCell, 1);
        mean_window_SS_lick_g(counter_pCell,1) = nanmean(SS_lick_g(counter_pCell, :));
        se_window_SS_lick_g(counter_pCell,1) = nanstd(SS_lick_g(counter_pCell, :))/sqrt(number_pCell);
        mean_change_window_SS_lick_g(counter_pCell,1) = nanmean(change_SS_lick_g(counter_pCell, :));
        se_change_window_SS_lick_g(counter_pCell,1) = nanstd(change_SS_lick_g(counter_pCell, :))/sqrt(number_pCell);
        
        SS_lick_l(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_l * 100), 15, 'sgolay', 2);
        change_SS_lick_l(counter_pCell,:) = SS_lick_l(counter_pCell,:) - SS_str_bout_baseline_l(counter_pCell, 1);
        mean_window_SS_lick_l(counter_pCell,1) = nanmean(SS_lick_l(counter_pCell, :));
        se_window_SS_lick_l(counter_pCell,1) = nanstd(SS_lick_l(counter_pCell, :))/sqrt(number_pCell);
        mean_change_window_SS_lick_l(counter_pCell,1) = nanmean(change_SS_lick_l(counter_pCell, :));
        se_change_window_SS_lick_l(counter_pCell,1) = nanstd(change_SS_lick_l(counter_pCell, :))/sqrt(number_pCell);
        
    elseif contains(alignment,'dmax') == 1
        %%% d_tip %%%
        d_tip_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_all);
        d_tip_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_grooming);
        d_tip_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_r);
        d_tip_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_l);
        
        %%% SS_str_bout %%%
        SS_str_bout_all(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_SS_str_bout * 100), 15, 'sgolay', 2);
        SS_str_bout_baseline_all(counter_pCell, :) = nanmean(SS_str_bout_all(counter_pCell, 1:50));
        
        SS_str_bout_r(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_SS_str_bout_r * 100), 15, 'sgolay', 2);
        SS_str_bout_baseline_r(counter_pCell, :) = nanmean(SS_str_bout_r(counter_pCell, 1:50));
        
        SS_str_bout_l(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_SS_str_bout_l * 100), 15, 'sgolay', 2);
        SS_str_bout_baseline_l(counter_pCell, :) = nanmean(SS_str_bout_l(counter_pCell, 1:50));
        
        %%% SS_lick %%%
        SS_lick_all(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_SS_all * 100), 15, 'sgolay', 2);
        change_SS_lick_all(counter_pCell,:) = SS_lick_all(counter_pCell,:) - SS_str_bout_baseline_all(counter_pCell, 1);
        mean_window_SS_lick_all(counter_pCell,1) = nanmean(SS_lick_all(counter_pCell, :));
        se_window_SS_lick_all(counter_pCell,1) = nanstd(SS_lick_all(counter_pCell, :))/sqrt(number_pCell);
        mean_change_window_SS_lick_all(counter_pCell,1) = nanmean(change_SS_lick_all(counter_pCell, :));
        se_change_window_SS_lick_all(counter_pCell,1) = nanstd(change_SS_lick_all(counter_pCell, :))/sqrt(number_pCell);
        
        SS_lick_r(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_SS_r * 100), 15, 'sgolay', 2);
        change_SS_lick_r(counter_pCell,:) = SS_lick_r(counter_pCell,:) - SS_str_bout_baseline_r(counter_pCell, 1);
        mean_window_SS_lick_r(counter_pCell,1) = nanmean(SS_lick_r(counter_pCell, :));
        se_window_SS_lick_r(counter_pCell,1) = nanstd(SS_lick_r(counter_pCell, :))/sqrt(number_pCell);
        mean_change_window_SS_lick_r(counter_pCell,1) = nanmean(change_SS_lick_r(counter_pCell, :));
        se_change_window_SS_lick_r(counter_pCell,1) = nanstd(change_SS_lick_r(counter_pCell, :))/sqrt(number_pCell);
        
        SS_lick_g(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_SS_grooming * 100), 15, 'sgolay', 2);
        change_SS_lick_g(counter_pCell,:) = SS_lick_g(counter_pCell,:) - SS_str_bout_baseline_all(counter_pCell, 1);
        mean_window_SS_lick_g(counter_pCell,1) = nanmean(SS_lick_g(counter_pCell, :));
        se_window_SS_lick_g(counter_pCell,1) = nanstd(SS_lick_g(counter_pCell, :))/sqrt(number_pCell);
        mean_change_window_SS_lick_g(counter_pCell,1) = nanmean(change_SS_lick_g(counter_pCell, :));
        se_change_window_SS_lick_g(counter_pCell,1) = nanstd(change_SS_lick_g(counter_pCell, :))/sqrt(number_pCell);
        
        SS_lick_l(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_SS_l * 100), 15, 'sgolay', 2);
        change_SS_lick_l(counter_pCell,:) = SS_lick_l(counter_pCell,:) - SS_str_bout_baseline_l(counter_pCell, 1);
        mean_window_SS_lick_l(counter_pCell,1) = nanmean(SS_lick_l(counter_pCell, :));
        se_window_SS_lick_l(counter_pCell,1) = nanstd(SS_lick_l(counter_pCell, :))/sqrt(number_pCell);
        mean_change_window_SS_lick_l(counter_pCell,1) = nanmean(change_SS_lick_l(counter_pCell, :));
        se_change_window_SS_lick_l(counter_pCell,1) = nanstd(change_SS_lick_l(counter_pCell, :))/sqrt(number_pCell);
    end
end

%%% Mean of d_tip_lick %%%
mean_d_tip_all = nanmean(d_tip_all);
se_d_tip_all = nanstd(d_tip_all)/sqrt(number_pCell);
mean_d_tip_l = nanmean(d_tip_l);
se_d_tip_l = nanstd(d_tip_all)/sqrt(number_pCell);
mean_d_tip_g = nanmean(d_tip_g);
se_d_tip_g = nanstd(d_tip_g)/sqrt(number_pCell);
mean_d_tip_r = nanmean(d_tip_r);
se_d_tip_r = nanstd(d_tip_r)/sqrt(number_pCell);

%%% Mean of SS_lick %%%
mean_SS_lick_all = nanmean(SS_lick_all);
se_SS_lick_all= nanstd(SS_lick_all)/sqrt(number_pCell);
mean_SS_lick_l = nanmean(SS_lick_l);
se_SS_lick_l= nanstd(SS_lick_l)/sqrt(number_pCell);
mean_SS_lick_g = nanmean(SS_lick_g);
se_SS_lick_g = nanstd(SS_lick_g)/sqrt(number_pCell);
mean_SS_lick_r = nanmean(SS_lick_r);
se_SS_lick_r = nanstd(SS_lick_r)/sqrt(number_pCell);

%%% Calculate diff in change based on direction
diff_mean_change_window_SS_lick_rl = mean_change_window_SS_lick_r - mean_change_window_SS_lick_l;
diff_mean_change_window_SS_lick_rg = mean_change_window_SS_lick_r - mean_change_window_SS_lick_g;
diff_mean_change_window_SS_lick_lg = mean_change_window_SS_lick_l - mean_change_window_SS_lick_g;

%%% lick labels %%%
label_large_increase_all       = mean_change_window_SS_lick_all >= 30;
label_medium_increase_all        = mean_change_window_SS_lick_all >= 15 & mean_change_window_SS_lick_all < 30 ;
label_small_increase_all      = mean_change_window_SS_lick_all >= 5 & mean_change_window_SS_lick_all < 15;
label_none_all                   = mean_change_window_SS_lick_all > -5 & mean_change_window_SS_lick_all < 5;
label_small_decrease_all         = mean_change_window_SS_lick_all > -15 & mean_change_window_SS_lick_all <= -5;
label_medium_decrease_all        = mean_change_window_SS_lick_all > -30 & mean_change_window_SS_lick_all <= -15;
label_large_decrease_all         = mean_change_window_SS_lick_all <= -30;

%%% lick labels R - L %%%
label_large_r       = diff_mean_change_window_SS_lick_rl >= 30;
label_medium_r        = diff_mean_change_window_SS_lick_rl >= 15 & diff_mean_change_window_SS_lick_rl < 30 ;
label_small_r      = diff_mean_change_window_SS_lick_rl >= 5 & diff_mean_change_window_SS_lick_rl < 15;
label_none_dir                   = diff_mean_change_window_SS_lick_rl > -5 & diff_mean_change_window_SS_lick_rl < 5;
label_small_l         = diff_mean_change_window_SS_lick_rl > -15 & diff_mean_change_window_SS_lick_rl <= -5;
label_medium_l        = diff_mean_change_window_SS_lick_rl > -30 & diff_mean_change_window_SS_lick_rl <= -15;
label_large_l         = diff_mean_change_window_SS_lick_rl <= -30;

%%% Mean d_tip based on labels
mean_d_tip_large_increase_all = nanmean(d_tip_all(label_large_increase_all,:));
se_d_tip_large_increase_all = nanstd(d_tip_all(label_large_increase_all,:))/sqrt(sum(label_large_increase_all));
mean_d_tip_medium_increase_all = nanmean(d_tip_all(label_medium_increase_all,:));
se_d_tip_medium_increase_all = nanstd(d_tip_all(label_medium_increase_all,:))/sqrt(sum(label_medium_increase_all));
mean_d_tip_small_increase_all = nanmean(d_tip_all(label_small_increase_all,:));
se_d_tip_small_increase_all = nanstd(d_tip_all(label_small_increase_all,:))/sqrt(sum(label_small_increase_all));
mean_d_tip_none_all = nanmean(d_tip_all(label_none_all,:));
se_d_tip_none_all = nanstd(d_tip_all(label_none_all,:))/sqrt(sum(label_none_all));
mean_d_tip_large_decrease_all = nanmean(d_tip_all(label_large_decrease_all,:));
se_d_tip_large_decrease_all = nanstd(d_tip_all(label_large_decrease_all,:))/sqrt(sum(label_large_decrease_all));
mean_d_tip_medium_decrease_all = nanmean(d_tip_all(label_medium_decrease_all,:));
se_d_tip_medium_decrease_all = nanstd(d_tip_all(label_medium_decrease_all,:))/sqrt(sum(label_medium_decrease_all));
mean_d_tip_small_decrease_all = nanmean(d_tip_all(label_small_decrease_all,:));
se_d_tip_small_decrease_all = nanstd(d_tip_all(label_small_decrease_all,:))/sqrt(sum(label_small_decrease_all));

%%% Mean change_SS_lick based on labels
mean_change_SS_lick_large_increase_all = nanmean(change_SS_lick_all(label_large_increase_all,:));
se_change_SS_lick_large_increase_all = nanstd(change_SS_lick_all(label_large_increase_all,:))/sqrt(sum(label_large_increase_all));
mean_change_SS_lick_medium_increase_all = nanmean(change_SS_lick_all(label_medium_increase_all,:));
se_change_SS_lick_medium_increase_all = nanstd(change_SS_lick_all(label_medium_increase_all,:))/sqrt(sum(label_medium_increase_all));
mean_change_SS_lick_small_increase_all = nanmean(change_SS_lick_all(label_small_increase_all,:));
se_change_SS_lick_small_increase_all = nanstd(change_SS_lick_all(label_small_increase_all,:))/sqrt(sum(label_small_increase_all));
mean_change_SS_lick_none_all = nanmean(change_SS_lick_all(label_none_all,:));
se_change_SS_lick_none_all = nanstd(change_SS_lick_all(label_none_all,:))/sqrt(sum(label_none_all));
mean_change_SS_lick_large_decrease_all = nanmean(change_SS_lick_all(label_large_decrease_all,:));
se_change_SS_lick_large_decrease_all = nanstd(change_SS_lick_all(label_large_decrease_all,:))/sqrt(sum(label_large_decrease_all));
mean_change_SS_lick_medium_decrease_all = nanmean(change_SS_lick_all(label_medium_decrease_all,:));
se_change_SS_lick_medium_decrease_all = nanstd(change_SS_lick_all(label_medium_decrease_all,:))/sqrt(sum(label_medium_decrease_all));
mean_change_SS_lick_small_decrease_all = nanmean(change_SS_lick_all(label_small_decrease_all,:));
se_change_SS_lick_small_decrease_all = nanstd(change_SS_lick_all(label_small_decrease_all,:))/sqrt(sum(label_small_decrease_all));

%%% Mean d_tip_bout based on labels r
mean_d_tip_large_r_r = nanmean(d_tip_r(label_large_r,:));
se_d_tip_large_r_r = nanstd(d_tip_r(label_large_r,:))/sqrt(sum(label_large_r));
mean_d_tip_large_r_l = nanmean((d_tip_l(label_large_r,:)));
se_d_tip_large_r_l = nanstd(d_tip_l(label_large_r,:))/sqrt(sum(label_large_r));
mean_d_tip_medium_r_r = nanmean(d_tip_r(label_medium_r,:));
se_d_tip_medium_r_r = nanstd(d_tip_r(label_medium_r,:))/sqrt(sum(label_medium_r));
mean_d_tip_medium_r_l = nanmean(d_tip_l(label_medium_r,:));
se_d_tip_medium_r_l = nanstd(d_tip_l(label_medium_r,:))/sqrt(sum(label_medium_r));
mean_d_tip_small_r_r = nanmean(d_tip_r(label_small_r,:));
se_d_tip_small_r_r = nanstd(d_tip_r(label_small_r,:))/sqrt(sum(label_small_r));
mean_d_tip_small_r_l = nanmean(d_tip_l(label_small_r,:));
se_d_tip_small_r_l = nanstd(d_tip_l(label_small_r,:))/sqrt(sum(label_small_r));

%%% Mean d_tip_str_bout based on labels l
mean_d_tip_large_l_r = nanmean(d_tip_r(label_large_l,:));
se_d_tip_large_l_r = nanstd(d_tip_r(label_large_l,:))/sqrt(sum(label_large_l));
mean_d_tip_large_l_l = nanmean(d_tip_l(label_large_l,:));
se_d_tip_large_l_l = nanstd(d_tip_l(label_large_l,:))/sqrt(sum(label_large_l));
mean_d_tip_medium_l_r = nanmean(d_tip_r(label_medium_l,:));
se_d_tip_medium_l_r = nanstd(d_tip_r(label_medium_l,:))/sqrt(sum(label_medium_l));
mean_d_tip_medium_l_l = nanmean(d_tip_l(label_medium_l,:));
se_d_tip_medium_l_l = nanstd(d_tip_l(label_medium_l,:))/sqrt(sum(label_medium_l));
mean_d_tip_small_l_r = nanmean(d_tip_r(label_small_l,:));
se_d_tip_small_l_r = nanstd(d_tip_r(label_small_l,:))/sqrt(sum(label_small_l));
mean_d_tip_small_l_l = nanmean(d_tip_l(label_small_l,:));
se_d_tip_small_l_l = nanstd(d_tip_l(label_small_l,:))/sqrt(sum(label_small_l));

%%% Mean d_tip_str_bout based on labels none
mean_d_tip_none_r = nanmean(d_tip_r(label_none_dir,:));
se_d_tip_none_r = nanstd(d_tip_r(label_none_dir,:))/sqrt(sum(label_none_dir));
mean_d_tip_none_l = nanmean(d_tip_l(label_none_dir,:));
se_d_tip_none_l = nanstd(d_tip_l(label_none_dir,:))/sqrt(sum(label_none_dir));
%%% Mean change_SS_lick based on labels all
mean_change_SS_lick_large_increase_all = nanmean(change_SS_lick_all(label_large_increase_all,:));
se_change_SS_lick_large_increase_all = nanstd(change_SS_lick_all(label_large_increase_all,:))/sqrt(sum(label_large_increase_all));
mean_change_SS_lick_medium_increase_all = nanmean(change_SS_lick_all(label_medium_increase_all,:));
se_change_SS_lick_medium_increase_all = nanstd(change_SS_lick_all(label_medium_increase_all,:))/sqrt(sum(label_medium_increase_all));
mean_change_SS_lick_small_increase_all = nanmean(change_SS_lick_all(label_small_increase_all,:));
se_change_SS_lick_small_increase_all = nanstd(change_SS_lick_all(label_small_increase_all,:))/sqrt(sum(label_small_increase_all));
mean_change_SS_lick_none_all = nanmean(change_SS_lick_all(label_none_all,:));
se_change_SS_lick_none_all = nanstd(change_SS_lick_all(label_none_all,:))/sqrt(sum(label_none_all));
mean_change_SS_lick_large_decrease_all = nanmean(change_SS_lick_all(label_large_decrease_all,:));
se_change_SS_lick_large_decrease_all = nanstd(change_SS_lick_all(label_large_decrease_all,:))/sqrt(sum(label_large_decrease_all));
mean_change_SS_lick_medium_decrease_all = nanmean(change_SS_lick_all(label_medium_decrease_all,:));
se_change_SS_lick_medium_decrease_all = nanstd(change_SS_lick_all(label_medium_decrease_all,:))/sqrt(sum(label_medium_decrease_all));
mean_change_SS_lick_small_decrease_all = nanmean(change_SS_lick_all(label_small_decrease_all,:));
se_change_SS_lick_small_decrease_all = nanstd(change_SS_lick_all(label_small_decrease_all,:))/sqrt(sum(label_small_decrease_all));

%%% Mean d_tip_lick based on labels r
mean_change_SS_lick_large_r_r = nanmean((change_SS_lick_r(label_large_r,:)));
se_change_SS_lick_large_r_r = nanstd(change_SS_lick_r(label_large_r,:))/sqrt(sum(label_large_r));
mean_change_SS_lick_large_r_l = nanmean((change_SS_lick_l(label_large_r,:)));
se_change_SS_lick_large_r_l = nanstd(change_SS_lick_l(label_large_r,:))/sqrt(sum(label_large_r));
mean_change_SS_lick_medium_r_r = nanmean(change_SS_lick_r(label_medium_r,:));
se_change_SS_lick_medium_r_r = nanstd(change_SS_lick_r(label_medium_r,:))/sqrt(sum(label_medium_r));
mean_change_SS_lick_medium_r_l = nanmean(change_SS_lick_l(label_medium_r,:));
se_change_SS_lick_medium_r_l = nanstd(change_SS_lick_l(label_medium_r,:))/sqrt(sum(label_medium_r));
mean_change_SS_lick_small_r_r = nanmean(change_SS_lick_r(label_small_r,:));
se_change_SS_lick_small_r_r = nanstd(change_SS_lick_r(label_small_r,:))/sqrt(sum(label_small_r));
mean_change_SS_lick_small_r_l = nanmean(change_SS_lick_l(label_small_r,:));
se_change_SS_lick_small_r_l = nanstd(change_SS_lick_l(label_small_r,:))/sqrt(sum(label_small_r));

%%% Mean change_SS_lick based on labels l
mean_change_SS_lick_large_l_r = nanmean(change_SS_lick_r(label_large_l,:));
se_change_SS_lick_large_l_r = nanstd(change_SS_lick_r(label_large_l,:))/sqrt(sum(label_large_l));
mean_change_SS_lick_large_l_l = nanmean(change_SS_lick_l(label_large_l,:));
se_change_SS_lick_large_l_l = nanstd(change_SS_lick_l(label_large_l,:))/sqrt(sum(label_large_l));
mean_change_SS_lick_medium_l_r = nanmean(change_SS_lick_r(label_medium_l,:));
se_change_SS_lick_medium_l_r = nanstd(change_SS_lick_r(label_medium_l,:))/sqrt(sum(label_medium_l));
mean_change_SS_lick_medium_l_l = nanmean(change_SS_lick_l(label_medium_l,:));
se_change_SS_lick_medium_l_l = nanstd(change_SS_lick_l(label_medium_l,:))/sqrt(sum(label_medium_l));
mean_change_SS_lick_small_l_r = nanmean(change_SS_lick_r(label_small_l,:));
se_change_SS_lick_small_l_r = nanstd(change_SS_lick_r(label_small_l,:))/sqrt(sum(label_small_l));
mean_change_SS_lick_small_l_l = nanmean(change_SS_lick_l(label_small_l,:));
se_change_SS_lick_small_l_l = nanstd(change_SS_lick_l(label_small_l,:))/sqrt(sum(label_small_l));

%%% Mean change_SS_lick based on labels none
mean_change_SS_lick_none_r = nanmean(change_SS_lick_r(label_none_dir,:));
se_change_SS_lick_none_r = nanstd(change_SS_lick_r(label_none_dir,:))/sqrt(sum(label_none_dir));
mean_change_SS_lick_none_l = nanmean(change_SS_lick_l(label_none_dir,:));
se_change_SS_lick_none_l = nanstd(change_SS_lick_l(label_none_dir,:))/sqrt(sum(label_none_dir));

%%% inds_span_bout %%%
inds_span = -290 : 10 : 300;

%%% xaxis for bar plot %%%
xaxis_pCell = 1:62;

%%% SAVE TO POPULATION STRUCTURE %%%
if contains(alignment,'onset') == 1
    POPULATION.onset.label_large_increase_all = label_large_increase_all;
    POPULATION.onset.label_medium_increase_all = label_medium_increase_all;
    POPULATION.onset.label_small_increase_all = label_small_increase_all;
    POPULATION.onset.label_none_all = label_none_all;
    POPULATION.onset.label_small_decrease_all = label_small_decrease_all;
    POPULATION.onset.label_medium_decrease_all = label_medium_decrease_all;
    POPULATION.onset.label_large_decrease_all = label_large_decrease_all;
    POPULATION.onset.label_large_r = label_large_r;
    POPULATION.onset.label_medium_r = label_medium_r;
    POPULATION.onset.label_small_r = label_small_r;
    POPULATION.onset.label_none_dir = label_none_dir;
    POPULATION.onset.label_small_l = label_small_l;
    POPULATION.onset.label_medium_l = label_medium_l;
    POPULATION.onset.label_large_l = label_large_l;
elseif contains(alignment,'dmax') == 1
    POPULATION.dmax.label_large_increase_all = label_large_increase_all;
    POPULATION.dmax.label_medium_increase_all = label_medium_increase_all;
    POPULATION.dmax.label_small_increase_all = label_small_increase_all;
    POPULATION.dmax.label_none_all = label_none_all;
    POPULATION.dmax.label_small_decrease_all = label_small_decrease_all;
    POPULATION.dmax.label_medium_decrease_all = label_medium_decrease_all;
    POPULATION.dmax.label_large_decrease_all = label_large_decrease_all;
    POPULATION.dmax.label_large_r = label_large_r;
    POPULATION.dmax.label_medium_r = label_medium_r;
    POPULATION.dmax.label_small_r = label_small_r;
    POPULATION.dmax.label_none_dir = label_none_dir;
    POPULATION.dmax.label_small_l = label_small_l;
    POPULATION.dmax.label_medium_l = label_medium_l;
    POPULATION.dmax.label_large_l = label_large_l;
end
%% Raw
num_pCell = 7;

figure
sgtitle(['pCells: ' num2str(number_pCell)])
subplot(4,4,1)
hold on
for counter_pCell = num_pCell
    plot(inds_span, d_tip_all(counter_pCell,:),'Color', [0.7 0.7 0.7] , 'LineWidth' , 2 )
end
%     plot(inds_span, mean_d_tip_all , 'Color', [0.7 0.7 0.7], 'LineWidth' , 2 )
%     plot(inds_span, mean_d_tip_all + se_d_tip_all , 'Color', [0.7 0.7 0.7], 'LineWidth' , 1 )
%     plot(inds_span, mean_d_tip_all - se_d_tip_all , 'Color', [0.7 0.7 0.7], 'LineWidth' ,1 )
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-400 400])
ylim([0 18])
xlabel(['Lick ' alignment  ' (ms)'])
ylabel('Displacement (mm)')
title(['Displacement | All | pCell: ' num2str(counter_pCell)])

subplot(4,4,2)
hold on
for counter_pCell = num_pCell
    plot(inds_span, d_tip_l(counter_pCell,:),'Color', [0.7 0.7 0.7] , 'LineWidth' , 2 )
end
%     plot(inds_span, mean_d_tip_l , 'Color', [0.7 0.7 0.7], 'LineWidth' , 2 )
%     plot(inds_span, mean_d_tip_l + se_d_tip_l , 'Color', [0.7 0.7 0.7], 'LineWidth' , 1 )
%     plot(inds_span, mean_d_tip_l - se_d_tip_l , 'Color', [0.7 0.7 0.7], 'LineWidth' ,1 )
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-400 400])
ylim([0 18])
xlabel(['Lick ' alignment  ' (ms)'])
ylabel('Displacement (mm)')
title(['Displacement | L | pCell: ' num2str(counter_pCell)])

subplot(4,4,3)
hold on
for counter_pCell = num_pCell
    plot(inds_span, d_tip_g(counter_pCell,:),'Color', [0.7 0.7 0.7] , 'LineWidth' , 2 )
end
%     plot(inds_span, mean_d_tip_g , 'Color',[0.7 0.7 0.7], 'LineWidth' , 2 )
%     plot(inds_span, mean_d_tip_g + se_d_tip_g , 'Color', [0.7 0.7 0.7], 'LineWidth' , 1 )
%     plot(inds_span, mean_d_tip_g - se_d_tip_g , 'Color', [0.7 0.7 0.7], 'LineWidth' ,1 )
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-400 400])
ylim([0 18])
xlabel(['Lick ' alignment  ' (ms)'])
ylabel('Displacement (mm)')
title(['Displacement | G  | pCell: ' num2str(counter_pCell)])

subplot(4,4,4)
hold on
for counter_pCell = num_pCell
    plot(inds_span, d_tip_r(counter_pCell,:),'Color', [0.7 0.7 0.7] , 'LineWidth' , 2 )
end
%     plot(inds_span, mean_d_tip_r , 'Color',[0.7 0.7 0.7], 'LineWidth' , 2 )
%     plot(inds_span, mean_d_tip_r + se_d_tip_r , 'Color', [0.7 0.7 0.7], 'LineWidth' , 1 )
%     plot(inds_span, mean_d_tip_r - se_d_tip_r , 'Color', [0.7 0.7 0.7], 'LineWidth' ,1 )
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-400 400])
ylim([0 18])
xlabel(['Lick ' alignment  ' (ms)'])
ylabel('Displacement (mm)')
title(['Displacement |  R | pCell: ' num2str(counter_pCell)])

subplot(4,4,5)
hold on
for counter_pCell = num_pCell
    plot(inds_span, smooth(change_SS_lick_all(counter_pCell,:),5,'sgolay', 2),'r' , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-400 400])
ylim([-50 200])
xlabel(['Lick ' alignment  ' (ms)'])
ylabel('Change in SS firing (spks/s)')
title(['SS firing | All | pCell: ' num2str(counter_pCell)])

subplot(4,4,6)
hold on
for counter_pCell = num_pCell
    plot(inds_span, change_SS_lick_l(counter_pCell,:),'r' , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-400 400])
ylim([-50 200])
xlabel(['Lick ' alignment  ' (ms)'])
ylabel('Change in SS firing (spks/s)')
title(['SS firing | L | pCell: ' num2str(counter_pCell)])

subplot(4,4,7)
hold on
for counter_pCell = num_pCell
    plot(inds_span, change_SS_lick_g(counter_pCell,:),'r' , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-400 400])
ylim([-50 200])
xlabel(['Lick ' alignment  ' (ms)'])
ylabel('Change in SS firing (spks/s)')
title(['SS firing | G | pCell: ' num2str(counter_pCell)])

subplot(4,4,8)
hold on
for counter_pCell = num_pCell
    plot(inds_span, change_SS_lick_r(counter_pCell,:),'r' , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-400 400])
ylim([-50 200])
xlabel(['Lick ' alignment  ' (ms)'])
ylabel('Change in SS firing (spks/s)')
title(['SS firing | R | pCell: ' num2str(counter_pCell)])

subplot(4,4,9)
hold on
for counter_pCell = 1 : number_pCell
    plot(SS_str_bout_baseline_all(counter_pCell), mean_window_SS_lick_all(counter_pCell), 'ok')
end
plot(SS_str_bout_baseline_all(num_pCell), mean_window_SS_lick_all(num_pCell), 'or')
vert = refline(1,0);
vert.Color = 'k';
xlabel('Mean SS firing before bout (spks/s)')
ylabel('Mean SS firing during lick (spks/s)')
xlim([0 250])
ylim([0 250])
title('Lick effect on SS firing | All licks')

subplot(4,4,10)
hold on
for counter_pCell = 1 : number_pCell
    plot(SS_str_bout_baseline_l(counter_pCell), mean_window_SS_lick_l(counter_pCell), 'ok')
end
plot(SS_str_bout_baseline_l(num_pCell), mean_window_SS_lick_l(num_pCell), 'or')
vert = refline(1,0);
vert.Color = 'k';
xlabel('Mean SS firing before bout (spks/s)')
ylabel('Mean SS firing during lick (spks/s)')
xlim([0 250])
ylim([0 250])
title('Lick effect on SS firing | L licks')

subplot(4,4,11)
hold on
for counter_pCell = 1 : number_pCell
    plot(SS_str_bout_baseline_all(counter_pCell), mean_window_SS_lick_g(counter_pCell), 'ok')
end
plot(SS_str_bout_baseline_all(num_pCell), mean_window_SS_lick_g(num_pCell), 'or')
vert = refline(1,0);
vert.Color = 'k';
xlabel('Mean SS firing before bout (spks/s)')
ylabel('Mean SS firing during lick (spks/s)')
xlim([0 250])
ylim([0 250])
title('Lick effect on SS firing | G licks')

subplot(4,4,12)
hold on
for counter_pCell = 1 : number_pCell
    plot(SS_str_bout_baseline_r(counter_pCell), mean_window_SS_lick_r(counter_pCell), 'ok')
end
plot(SS_str_bout_baseline_r(num_pCell), mean_window_SS_lick_r(num_pCell), 'or')
vert = refline(1,0);
vert.Color = 'k';
xlabel('Mean SS firing before bout (spks/s)')
ylabel('Mean SS firing during lick (spks/s)')
xlim([0 250])
ylim([0 250])
title('Lick effect on SS firing | R licks')

subplot(4,4,13)
hold on;
bar(mean_change_window_SS_lick_all, 'k')
bar(xaxis_pCell(num_pCell), mean_change_window_SS_lick_all(num_pCell), 'r')

% plot(mean_change_window_SS_lick_all + se_change_window_SS_lick_all, '*r')
% plot(mean_change_window_SS_lick_all - se_change_window_SS_lick_all, '*r')
ylim([-20 115])
xlabel('pCell')
ylabel('Change in SS firing (spks/s)')

subplot(4,4,14)
hold on;
bar(mean_change_window_SS_lick_l, 'k')
bar(xaxis_pCell(num_pCell), mean_change_window_SS_lick_l(num_pCell), 'r')
% plot(mean_change_window_SS_lick_all + se_change_window_SS_lick_all, '*r')
% plot(mean_change_window_SS_lick_all - se_change_window_SS_lick_all, '*r')
ylim([-20 115])
xlabel('pCell')
ylabel('Change in SS firing (spks/s)')

subplot(4,4,15)
hold on;
bar(mean_change_window_SS_lick_g, 'k')
bar(xaxis_pCell(num_pCell), mean_change_window_SS_lick_g(num_pCell), 'r')
% plot(mean_change_window_SS_lick_all + se_change_window_SS_lick_all, '*r')
% plot(mean_change_window_SS_lick_all - se_change_window_SS_lick_all, '*r')
ylim([-20 115])
xlabel('pCell')
ylabel('Change in SS firing (spks/s)')

subplot(4,4,16)
hold on;
bar(mean_change_window_SS_lick_r, 'k')
bar(xaxis_pCell(num_pCell), mean_change_window_SS_lick_r(num_pCell), 'r')
% plot(mean_change_window_SS_lick_all + se_change_window_SS_lick_all, '*r')
% plot(mean_change_window_SS_lick_all - se_change_window_SS_lick_all, '*r')
ylim([-20 115])
xlabel('pCell')
ylabel('Change in SS firing (spks/s)')
%% Labeled properties
xaxis_pCell = 1:62;

figure
subplot(2,5,[1 2])
hold on;
bar(xaxis_pCell(label_large_increase_all ), mean_change_window_SS_lick_all(label_large_increase_all ), 'r', 'Barwidth', 1)
bar(xaxis_pCell(label_medium_increase_all ), mean_change_window_SS_lick_all(label_medium_increase_all ), 'b', 'Barwidth', 1)
bar(xaxis_pCell(label_small_increase_all ), mean_change_window_SS_lick_all(label_small_increase_all ), 'g', 'Barwidth', 1)
bar(xaxis_pCell(label_none_all), mean_change_window_SS_lick_all(label_none_all ), 'k', 'Barwidth', 1)
bar(xaxis_pCell(label_small_decrease_all ), mean_change_window_SS_lick_all(label_small_decrease_all ), 'm', 'Barwidth', 0.20)
yline(30,'--r', 'LineWidth', 2);
yline(15,'--b', 'LineWidth', 2);
yline(5,'--g', 'LineWidth', 2);
yline(-5,'--m', 'LineWidth', 2);
ylim([-20 115])
xlabel('pCell')
ylabel('Change in SS firing (spks/s)')
title('Effect of lick on SS firing | All licks')

subplot(2,5,6)
hold on;
plot(SS_str_bout_baseline_all(label_large_increase_all), mean_window_SS_lick_all(label_large_increase_all), 'or')
plot(SS_str_bout_baseline_all(label_medium_increase_all), mean_window_SS_lick_all(label_medium_increase_all), 'ob')
plot(SS_str_bout_baseline_all(label_small_increase_all), mean_window_SS_lick_all(label_small_increase_all), 'og')
plot(SS_str_bout_baseline_all(label_none_all), mean_window_SS_lick_all(label_none_all), 'ok')
plot(SS_str_bout_baseline_all(label_small_decrease_all), mean_window_SS_lick_all(label_small_decrease_all), 'om')
vert = refline(1,0);
vert.Color = 'k';
xlabel('Mean SS firing before bout (spks/s)')
ylabel('Mean SS firing during lick window (spks/s)')
xlim([0 250])
ylim([0 250])

pie_1 = subplot(2,5,7);
pie_values = [sum(label_large_increase_all) sum(label_medium_increase_all) sum(label_small_increase_all) sum(label_none_all) sum(label_small_decrease_all)];
pie_labels = {num2str(sum(label_large_increase_all)) ...
    num2str(sum(label_medium_increase_all)) ...
    num2str(sum(label_small_increase_all)) ...
    num2str(sum(label_none_all)) ...
    num2str(sum(label_small_decrease_all))};
explode_1 = [1 1 1 1 1];
pie(pie_values, explode_1, pie_labels)
colormap(pie_1, [1 0 0; 0 0 1; 0 1 0; 0 0 0; 1 0 1]);

subplot(2,5,[3 4 8 9])
hold on;
barh(xaxis_pCell(label_large_r), diff_mean_change_window_SS_lick_rl(label_large_r), 'r', 'BaseValue', 0, 'Barwidth', 0.06)
barh(xaxis_pCell(label_medium_r), diff_mean_change_window_SS_lick_rl(label_medium_r), 'b', 'BaseValue', 0, 'Barwidth', 1)
barh(xaxis_pCell(label_small_r), diff_mean_change_window_SS_lick_rl(label_small_r), 'g', 'BaseValue', 0, 'Barwidth', 1)
barh(xaxis_pCell(label_none_dir), diff_mean_change_window_SS_lick_rl(label_none_dir), 'k', 'BaseValue', 0, 'Barwidth', 1)
barh(xaxis_pCell(label_small_l), diff_mean_change_window_SS_lick_rl(label_small_l), 'g', 'BaseValue', 0, 'Barwidth', 1)
barh(xaxis_pCell(label_medium_l), diff_mean_change_window_SS_lick_rl(label_medium_l), 'b', 'BaseValue', 0, 'Barwidth', 1)
barh(xaxis_pCell(label_large_l), diff_mean_change_window_SS_lick_rl(label_large_l), 'r', 'BaseValue', 0, 'Barwidth', 0.2)
xline(30,'--r', 'LineWidth', 2);
xline(-30,'--r', 'LineWidth', 2);
xline(15,'--b', 'LineWidth', 2);
xline(-15,'--b', 'LineWidth', 2);
xline(5,'--g', 'LineWidth', 2);
xline(-5,'--g', 'LineWidth', 2);
xlabel('Difference in change of SS firing (spks/s)')
ylabel('pCell')
xlim([-70 70])
ylim([-3 65])
title('Effect of lick direction on SS firing | R - L licks')

subplot(2,5,5)
hold on;
plot(mean_change_window_SS_lick_r(label_large_r), mean_change_window_SS_lick_l(label_large_r), 'or')
plot(mean_change_window_SS_lick_r(label_medium_r), mean_change_window_SS_lick_l(label_medium_r), 'ob')
plot(mean_change_window_SS_lick_r(label_small_r), mean_change_window_SS_lick_l(label_small_r), 'og')
plot(mean_change_window_SS_lick_r(label_none_dir), mean_change_window_SS_lick_l(label_none_dir), 'ok')
plot(mean_change_window_SS_lick_r(label_small_l), mean_change_window_SS_lick_l(label_small_l), 'og')
plot(mean_change_window_SS_lick_r(label_medium_l), mean_change_window_SS_lick_l(label_medium_l), 'ob')
plot(mean_change_window_SS_lick_r(label_large_l), mean_change_window_SS_lick_l(label_large_l), 'or')
vert = refline(1,0);
vert.Color = 'k';
xlabel('Mean change in SS firing R lick (spks/s)')
ylabel('Mean change in SS firing L lick (spks/s)')
xlim([-30 125])
ylim([-30 125])

pie_2 = subplot(2,5,10);
pie_values = [sum(label_large_r) sum(label_medium_r) sum(label_small_r) sum(label_none_dir) ...
    sum(label_small_l) sum(label_medium_l) sum(label_large_l)];
pie_labels = {['R: ' num2str(sum(label_large_r))] ...
    ['R: ' num2str(sum(label_medium_r))] ...
    ['R: '  num2str(sum(label_small_r))] ...
    ['None: ' num2str(sum(label_none_dir))] ...
    ['L: ' num2str(sum(label_small_l))] ...
    ['L: ' num2str(sum(label_medium_l))] ...
    ['L: ' num2str(sum(label_large_l))]};
explode_2 = [1 1 1 1 1 1 1];
pie(pie_values, explode_2, pie_labels)
colormap(pie_2, [1 0 0; 0 0 1; 0 1 0; 0 0 0; 0 1 0; 0 0 1; 1 0 0]);
%% Labeled traces
%%% TRACES %%%
figure
sgtitle(['Lick | pCells: ' num2str(number_pCell)])
subplot(2,5,1)
hold on
plot(inds_span, mean_d_tip_large_increase_all(:), 'r', 'LineWidth', 2)
%     plot(inds_span, mean_d_tip_large_increase_all + se_d_tip_large_increase_all, 'r', 'LineWidth', 1)
%     plot(inds_span, mean_d_tip_large_increase_all - se_d_tip_large_increase_all, 'r', 'LineWidth', 1)
plot(inds_span, mean_d_tip_medium_increase_all(:), 'b', 'LineWidth', 2)
%     plot(inds_span, mean_d_tip_medium_increase_all + se_d_tip_medium_increase_all, 'b', 'LineWidth', 1)
%     plot(inds_span, mean_d_tip_medium_increase_all - se_d_tip_medium_increase_all, 'b', 'LineWidth', 1)
plot(inds_span, mean_d_tip_small_increase_all(:), 'g', 'LineWidth', 2)
%     plot(inds_span, mean_d_tip_small_increase_all + se_d_tip_small_increase_all, 'g', 'LineWidth', 1)
%     plot(inds_span, mean_d_tip_small_increase_all - se_d_tip_small_increase_all, 'g', 'LineWidth', 1)
plot(inds_span, mean_d_tip_none_all(:), 'k', 'LineWidth', 2)
%     plot(inds_span, mean_d_tip_none_all + se_d_tip_none_all, 'k', 'LineWidth', 1)
%     plot(inds_span, mean_d_tip_none_all - se_d_tip_none_all, 'k', 'LineWidth', 1)
plot(inds_span, mean_d_tip_small_decrease_all(:), 'm', 'LineWidth', 2)
%     plot(inds_span, mean_d_tip_small_decrease_all + se_d_tip_small_decrease_all, 'm', 'LineWidth', 1)
%     plot(inds_span, mean_d_tip_small_decrease_all - se_d_tip_small_decrease_all, 'm', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-400 400])
ylabel('Displacement (mm)')
xlabel(['Lick ' alignment  ' (ms)'])

subplot(2,5,2)
hold on
plot(inds_span, mean_d_tip_large_l_l(:), 'r', 'LineWidth', 2)
plot(inds_span, mean_d_tip_medium_l_l(:), 'b', 'LineWidth', 2)
plot(inds_span, mean_d_tip_small_l_l(:), 'g', 'LineWidth', 2)
plot(inds_span, mean_d_tip_none_l(:), 'k', 'LineWidth', 2)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-400 400])
ylabel('Displacement (mm)')
xlabel(['Lick ' alignment  ' (ms)'])
title(['Lick trace | L-tuned, L bouts'])

subplot(2,5,3)
hold on
plot(inds_span, mean_d_tip_large_l_r(:), 'r', 'LineWidth', 2)
plot(inds_span, mean_d_tip_medium_l_r(:), 'b', 'LineWidth', 2)
plot(inds_span, mean_d_tip_small_l_r(:), 'g', 'LineWidth', 2)
plot(inds_span, mean_d_tip_none_r(:), 'k', 'LineWidth', 2)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-400 400])
ylabel('Displacement (mm)')
xlabel(['Lick ' alignment  ' (ms)'])
title(['Lick trace | L-tuned, R bouts'])

subplot(2,5,4)
hold on
plot(inds_span, mean_d_tip_large_r_l(:), 'r', 'LineWidth', 2)
plot(inds_span, mean_d_tip_medium_r_l(:), 'b', 'LineWidth', 2)
plot(inds_span, mean_d_tip_small_r_l(:), 'g', 'LineWidth', 2)
plot(inds_span, mean_d_tip_none_r(:), 'k', 'LineWidth', 2)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-400 400])
ylabel('Displacement (mm)')
xlabel(['Lick ' alignment  ' (ms)'])
title(['Lick trace | R-tuned, L bouts'])

subplot(2,5,5)
hold on
plot(inds_span, mean_d_tip_large_r_r(:), 'r', 'LineWidth', 2)
plot(inds_span, mean_d_tip_medium_r_r(:), 'b', 'LineWidth', 2)
plot(inds_span, mean_d_tip_small_r_r(:), 'g', 'LineWidth', 2)
plot(inds_span, mean_d_tip_none_r(:), 'k', 'LineWidth', 2)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-400 400])
ylabel('Displacement (mm)')
xlabel(['Lick ' alignment  ' (ms)'])
title(['Lick trace | R-tuned, R bouts'])

subplot(2,5,6)
hold on
plot(inds_span, mean_change_SS_lick_large_increase_all(:), 'r', 'LineWidth', 2)
plot(inds_span, mean_change_SS_lick_large_increase_all + se_change_SS_lick_large_increase_all, 'r', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_large_increase_all - se_change_SS_lick_large_increase_all, 'r', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_medium_increase_all(:), 'b', 'LineWidth', 2)
plot(inds_span, mean_change_SS_lick_medium_increase_all + se_change_SS_lick_medium_increase_all, 'b', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_medium_increase_all - se_change_SS_lick_medium_increase_all, 'b', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_small_increase_all(:), 'g', 'LineWidth', 2)
plot(inds_span, mean_change_SS_lick_small_increase_all + se_change_SS_lick_small_increase_all, 'g', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_small_increase_all - se_change_SS_lick_small_increase_all, 'g', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_none_all(:), 'k', 'LineWidth', 2)
plot(inds_span, mean_change_SS_lick_none_all + se_change_SS_lick_none_all, 'k', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_none_all - se_change_SS_lick_none_all, 'k', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_small_decrease_all(:), 'm', 'LineWidth', 2)
plot(inds_span, mean_change_SS_lick_small_decrease_all + se_change_SS_lick_small_decrease_all, 'm', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_small_decrease_all - se_change_SS_lick_small_decrease_all, 'm', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
xlim([-400 400])
ylim([-30 120])
xlabel(['Lick ' alignment  ' (ms)'])
ylabel('Change in SS firing (spks/s)')
title(['Effect of lick on SS firing | All'])

subplot(2,5,7)
hold on
plot(inds_span, mean_change_SS_lick_large_l_l(:), 'r', 'LineWidth', 2)
plot(inds_span, mean_change_SS_lick_large_l_l(:) + se_change_SS_lick_large_l_l(:), 'r', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_large_l_l(:) - se_change_SS_lick_large_l_l(:), 'r', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_medium_l_l(:), 'b', 'LineWidth', 2)
plot(inds_span, mean_change_SS_lick_medium_l_l(:) + se_change_SS_lick_medium_l_l(:), 'b', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_medium_l_l(:) - se_change_SS_lick_medium_l_l(:), 'b', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_small_l_l(:), 'g', 'LineWidth', 2)
plot(inds_span, mean_change_SS_lick_small_l_l(:) + se_change_SS_lick_small_l_l(:), 'g', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_small_l_l(:) - se_change_SS_lick_small_l_l(:), 'g', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_none_l(:), 'k', 'LineWidth', 2)
plot(inds_span, mean_change_SS_lick_none_l(:) + se_change_SS_lick_none_l(:) , 'k', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_none_l(:) - se_change_SS_lick_none_l(:) , 'k', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
xlim([-400 400])
ylim([-30 150])
xlabel(['Lick ' alignment  ' (ms)'])
ylabel('Change in SS firing (spks/s)')
title(['Effect of lick on SS firing | L-tuned L bouts'])

subplot(2,5,8)
hold on
plot(inds_span, mean_change_SS_lick_large_l_r(:), 'r', 'LineWidth', 2)
plot(inds_span, mean_change_SS_lick_large_l_r(:) + se_change_SS_lick_large_l_r(:), 'r', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_large_l_r(:) - se_change_SS_lick_large_l_r(:), 'r', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_medium_l_r(:), 'b', 'LineWidth', 2)
plot(inds_span, mean_change_SS_lick_medium_l_r(:) + se_change_SS_lick_medium_l_r(:), 'b', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_medium_l_r(:) - se_change_SS_lick_medium_l_r(:), 'b', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_small_l_r(:), 'g', 'LineWidth', 2)
plot(inds_span, mean_change_SS_lick_small_l_r(:) + se_change_SS_lick_small_l_r(:), 'g', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_small_l_r(:) - se_change_SS_lick_small_l_r(:), 'g', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_none_r(:), 'k', 'LineWidth', 2)
plot(inds_span, mean_change_SS_lick_none_r(:) + se_change_SS_lick_none_r(:) , 'k', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_none_r(:) - se_change_SS_lick_none_r(:) , 'k', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
xlim([-400 400])
ylim([-30 150])
xlabel(['Lick ' alignment  ' (ms)'])
ylabel('Change in SS firing (spks/s)')
title(['Effect of lick on SS firing | L-tuned R bouts'])

subplot(2,5,9)
hold on
plot(inds_span, mean_change_SS_lick_large_r_l(:), 'r', 'LineWidth', 2)
plot(inds_span, mean_change_SS_lick_large_r_l(:) + se_change_SS_lick_large_r_l(:), 'r', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_large_r_l(:) - se_change_SS_lick_large_r_l(:), 'r', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_medium_r_l(:), 'b', 'LineWidth', 2)
plot(inds_span, mean_change_SS_lick_medium_r_l(:) + se_change_SS_lick_medium_r_l(:), 'b', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_medium_r_l(:) - se_change_SS_lick_medium_r_l(:), 'b', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_small_r_l(:), 'g', 'LineWidth', 2)
plot(inds_span, mean_change_SS_lick_small_r_l(:) + se_change_SS_lick_small_r_l(:), 'g', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_small_r_l(:) - se_change_SS_lick_small_r_l(:), 'g', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_none_l(:), 'k', 'LineWidth', 2)
plot(inds_span, mean_change_SS_lick_none_l(:) + se_change_SS_lick_none_l(:) , 'k', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_none_l(:) - se_change_SS_lick_none_l(:) , 'k', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
xlim([-400 400])
ylim([-30 150])
xlabel(['Lick ' alignment  ' (ms)'])
ylabel('Change in SS firing (spks/s)')
title(['Effect of lick on SS firing | R-tuned L bouts'])

subplot(2,5,10)
hold on
plot(inds_span, mean_change_SS_lick_large_r_r(:), 'r', 'LineWidth', 2)
plot(inds_span, mean_change_SS_lick_large_r_r(:) + se_change_SS_lick_large_r_r(:), 'r', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_large_r_r(:) - se_change_SS_lick_large_r_r(:), 'r', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_medium_r_r(:), 'b', 'LineWidth', 2)
plot(inds_span, mean_change_SS_lick_medium_r_r(:) + se_change_SS_lick_medium_r_r(:), 'b', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_medium_r_r(:) - se_change_SS_lick_medium_r_r(:), 'b', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_small_r_r(:), 'g', 'LineWidth', 2)
plot(inds_span, mean_change_SS_lick_small_r_r(:) + se_change_SS_lick_small_r_r(:), 'g', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_small_r_r(:) - se_change_SS_lick_small_r_r(:), 'g', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_none_r(:), 'k', 'LineWidth', 2)
plot(inds_span, mean_change_SS_lick_none_r(:) + se_change_SS_lick_none_r(:) , 'k', 'LineWidth', 1)
plot(inds_span, mean_change_SS_lick_none_r(:) - se_change_SS_lick_none_r(:) , 'k', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
xlim([-400 400])
ylim([-30 150])
xlabel(['Lick ' alignment  ' (ms)'])
ylabel('Change in SS firing (spks/s)')
title(['Effect of lick on SS firing | R-tuned R bouts'])

%% LICK CS RESPONSE
clearvars -except ALL_PCELL_COMPRESSED_DATA POPULATION

number_pCell =  length(ALL_PCELL_COMPRESSED_DATA);

%%% range_inds_prob [-50 50] %%%
range_inds_probability = 26:35;

for counter_pCell = 1 : number_pCell
    id(counter_pCell, 1) = string(ALL_PCELL_COMPRESSED_DATA(counter_pCell).id{1, 1});
   
    %%% Initialize based on ONSET event %%%
        %%% d_tip %%%
        d_tip_all_onset(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_all);
        d_tip_g_onset(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_grooming);
        d_tip_r_onset(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_r);
        d_tip_l_onset(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_l);
                
        %%% CS_lick FR %%%
        CS_lick_all_onset(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_all * 100), 5, 'sgolay', 2);      
        CS_lick_r_onset(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_r * 100), 5, 'sgolay', 2);
        CS_lick_g_onset(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_grooming * 100), 5, 'sgolay', 2);
        CS_lick_l_onset(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_l * 100), 5, 'sgolay', 2);

        %%% CS_lick prob %%%
        CS_lick_prob_all_onset(counter_pCell,:) = nansum(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_all) /  nansum(nansum(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_all));      
        CS_lick_prob_r_onset(counter_pCell,:) = nansum(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_r) /  nansum(nansum(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_r));      
        CS_lick_prob_g_onset(counter_pCell,:) = nansum(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_grooming) /  nansum(nansum(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_grooming));      
        CS_lick_prob_l_onset(counter_pCell,:) = nansum(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_l) /  nansum(nansum(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_CS_l));      
  
     %%% Initialize based on DMAX event %%%
        %%% d_tip %%%
        d_tip_all_dmax(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_all);
        d_tip_g_dmax(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_grooming);
        d_tip_r_dmax(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_r);
        d_tip_l_dmax(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_l);
        
        %%% CS_lick FR %%%
        CS_lick_all_dmax(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_all * 100), 5, 'sgolay', 2);      
        CS_lick_r_dmax(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_r * 100), 5, 'sgolay', 2);
        CS_lick_g_dmax(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_grooming * 100), 5, 'sgolay', 2);
        CS_lick_l_dmax(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_l * 100), 5, 'sgolay', 2);

        %%% CS_lick prob %%%
        CS_lick_prob_all_dmax(counter_pCell,:) = nansum(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_all) / nansum(nansum(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_all));      
        CS_lick_prob_r_dmax(counter_pCell,:) = nansum(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_r) / nansum(nansum(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_r));      
        CS_lick_prob_g_dmax(counter_pCell,:) = nansum(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_grooming) / nansum(nansum(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_grooming));      
        CS_lick_prob_l_dmax(counter_pCell,:) = nansum(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_l) / nansum(nansum(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_CS_l));      
  
end


%%% inds_span %%%
inds_span = -290 : 10 : 300;

%%% xaxis for bar plot %%%
xaxis_pCell = 1:62;
%% Raw
num_pCell = 34;

figure
sgtitle(['pCells: ' num2str(number_pCell)])
subplot(4,4,1)
hold on
for counter_pCell = num_pCell
    plot(inds_span, d_tip_all_onset(counter_pCell,:),'Color', [0.7 0.7 0.7] , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-300 300])
ylim([0 18])
xlabel(['Lick onset (ms)'])
ylabel('Displacement (mm)')
title(['Displacement | All | pCell: ' num2str(counter_pCell)])

subplot(4,4,2)
hold on
for counter_pCell = num_pCell
    plot(inds_span, d_tip_l_onset(counter_pCell,:),'Color', [0.7 0.7 0.7] , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-300 300])
ylim([0 18])
xlabel(['Lick onset (ms)'])
ylabel('Displacement (mm)')
title(['Displacement | L | pCell: ' num2str(counter_pCell)])

subplot(4,4,3)
hold on
for counter_pCell = num_pCell
    plot(inds_span, d_tip_g_onset(counter_pCell,:),'Color', [0.7 0.7 0.7] , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-300 300])
ylim([0 18])
xlabel(['Lick onset (ms)'])
ylabel('Displacement (mm)')
title(['Displacement | G  | pCell: ' num2str(counter_pCell)])

subplot(4,4,4)
hold on
for counter_pCell = num_pCell
    plot(inds_span, d_tip_r_onset(counter_pCell,:),'Color', [0.7 0.7 0.7] , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-300 300])
ylim([0 18])
xlabel(['Lick onset (ms)'])
ylabel('Displacement (mm)')
title(['Displacement |  R | pCell: ' num2str(counter_pCell)])

subplot(4,4,5)
hold on
for counter_pCell = num_pCell
    plot(inds_span, CS_lick_prob_all_onset(counter_pCell,:),'r' , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-300 300])
ylim([0 0.06])
xlabel(['Lick onset (ms)'])
ylabel('CS probability')
title(['CS firing | All | pCell: ' num2str(counter_pCell)])

subplot(4,4,6)
hold on
for counter_pCell = num_pCell
    plot(inds_span, CS_lick_prob_l_onset(counter_pCell,:),'r' , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-300 300])
ylim([0 0.06])
xlabel(['Lick onset (ms)'])
ylabel('CS probability')
title(['SS firing | L | pCell: ' num2str(counter_pCell)])

subplot(4,4,7)
hold on
for counter_pCell = num_pCell
    plot(inds_span, CS_lick_prob_g_onset(counter_pCell,:),'r' , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-300 300])
ylim([0 0.06])
xlabel(['Lick onset (ms)'])
ylabel('CS probability')
title(['SS firing | G | pCell: ' num2str(counter_pCell)])

subplot(4,4,8)
hold on
for counter_pCell = num_pCell
    plot(inds_span, CS_lick_prob_r_onset(counter_pCell,:),'r' , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-300 300])
ylim([0 0.06])
xlabel(['Lick onset (ms)'])
ylabel('CS probability')
title(['SS firing | R | pCell: ' num2str(counter_pCell)])

subplot(4,4,9)
hold on
for counter_pCell = num_pCell
    plot(inds_span, d_tip_all_dmax(counter_pCell,:),'Color', [0.7 0.7 0.7] , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-300 300])
ylim([0 18])
xlabel(['Lick dmax (ms)'])
ylabel('Displacement (mm)')
title(['Displacement | All | pCell: ' num2str(counter_pCell)])

subplot(4,4,10)
hold on
for counter_pCell = num_pCell
    plot(inds_span, d_tip_l_dmax(counter_pCell,:),'Color', [0.7 0.7 0.7] , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-300 300])
ylim([0 18])
xlabel(['Lick dmax (ms)'])
ylabel('Displacement (mm)')
title(['Displacement | L | pCell: ' num2str(counter_pCell)])

subplot(4,4,11)
hold on
for counter_pCell = num_pCell
    plot(inds_span, d_tip_g_dmax(counter_pCell,:),'Color', [0.7 0.7 0.7] , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-300 300])
ylim([0 18])
xlabel(['Lick dmax (ms)'])
ylabel('Displacement (mm)')
title(['Displacement | G  | pCell: ' num2str(counter_pCell)])

subplot(4,4,12)
hold on
for counter_pCell = num_pCell
    plot(inds_span, d_tip_r_dmax(counter_pCell,:),'Color', [0.7 0.7 0.7] , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-300 300])
ylim([0 18])
xlabel(['Lick dmax (ms)'])
ylabel('Displacement (mm)')
title(['Displacement |  R | pCell: ' num2str(counter_pCell)])

subplot(4,4,13)
hold on
for counter_pCell = num_pCell
    plot(inds_span, CS_lick_prob_all_dmax(counter_pCell,:),'r' , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-300 300])
ylim([0 0.06])
xlabel(['Lick dmax (ms)'])
ylabel('CS probability')
title(['CS firing | All | pCell: ' num2str(counter_pCell)])

subplot(4,4,14)
hold on
for counter_pCell = num_pCell
    plot(inds_span, CS_lick_prob_l_dmax(counter_pCell,:),'r' , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-300 300])
ylim([0 0.06])
xlabel(['Lick dmax (ms)'])
ylabel('CS probability')
title(['SS firing | L | pCell: ' num2str(counter_pCell)])

subplot(4,4,15)
hold on
for counter_pCell = num_pCell
    plot(inds_span, CS_lick_prob_g_dmax(counter_pCell,:),'r' , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-300 300])
ylim([0 0.06])
xlabel(['Lick dmax (ms)'])
ylabel('CS probability')
title(['SS firing | G | pCell: ' num2str(counter_pCell)])

subplot(4,4,16)
hold on
for counter_pCell = num_pCell
    plot(inds_span, CS_lick_prob_r_dmax(counter_pCell,:),'r' , 'LineWidth' , 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([-300 300])
ylim([0 0.06])
xlabel(['Lick dmax (ms)'])
ylabel('CS probability')
title(['SS firing | R | pCell: ' num2str(counter_pCell)])


%% LICK SS PHASE & FREQUENCY
clearvars -except ALL_PCELL_COMPRESSED_DATA POPULATION

number_pCell =  length(ALL_PCELL_COMPRESSED_DATA);

alignment = 'dmax';
% alignment = 'onset';

for counter_pCell = 1 : number_pCell
    if contains(alignment,'onset') == 1
        id(counter_pCell, 1) = string(ALL_PCELL_COMPRESSED_DATA(counter_pCell).id{1, 1});
        
        %%% d_tip %%%
        d_tip_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_all_corr);
        d_tip_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_grooming_corr);
        d_tip_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_r_corr);
        d_tip_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_l_corr);
        
        %%% SS_str_bout %%%
        SS_str_bout_all(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_str_bout * 100), 15, 'sgolay', 2);
        SS_str_bout_baseline_all(counter_pCell, :) = nanmean(SS_str_bout_all(counter_pCell, 1:50));
        SS_str_bout_r(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_str_bout_r * 100), 15, 'sgolay', 2);
        SS_str_bout_baseline_r(counter_pCell, :) = nanmean(SS_str_bout_r(counter_pCell, 1:50));
        SS_str_bout_l(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_str_bout_l * 100), 15, 'sgolay', 2);
        SS_str_bout_baseline_l(counter_pCell, :) = nanmean(SS_str_bout_l(counter_pCell, 1:50));
        
        %%% SS_lick %%%
        SS_lick_all(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_all_corr * 100), 15, 'sgolay', 2);
        change_SS_lick_all(counter_pCell,:) = SS_lick_all(counter_pCell,:) - SS_str_bout_baseline_all(counter_pCell, 1);
        SS_lick_r(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_r_corr * 100), 15, 'sgolay', 2);
        change_SS_lick_r(counter_pCell,:) = SS_lick_r(counter_pCell,:) - SS_str_bout_baseline_r(counter_pCell, 1);
        SS_lick_g(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_grooming_corr * 100), 15, 'sgolay', 2);
        change_SS_lick_g(counter_pCell,:) = SS_lick_g(counter_pCell,:) - SS_str_bout_baseline_all(counter_pCell, 1);
        SS_lick_l(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_l_corr * 100), 15, 'sgolay', 2);
        change_SS_lick_l(counter_pCell,:) = SS_lick_l(counter_pCell,:) - SS_str_bout_baseline_l(counter_pCell, 1);
        
    elseif contains(alignment,'dmax') == 1
        %%% d_tip %%%
        d_tip_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_all_corr);
        d_tip_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_grooming_corr);
        d_tip_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_r_corr);
        d_tip_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_l_corr);
        
        %%% SS_str_bout %%%
        SS_str_bout_all(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_SS_str_bout * 100), 15, 'sgolay', 2);
        SS_str_bout_baseline_all(counter_pCell, :) = nanmean(SS_str_bout_all(counter_pCell, 1:50));
        SS_str_bout_r(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_SS_str_bout_r * 100), 15, 'sgolay', 2);
        SS_str_bout_baseline_r(counter_pCell, :) = nanmean(SS_str_bout_r(counter_pCell, 1:50));
        SS_str_bout_l(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_SS_str_bout_l * 100), 15, 'sgolay', 2);
        SS_str_bout_baseline_l(counter_pCell, :) = nanmean(SS_str_bout_l(counter_pCell, 1:50));
        
        %%% SS_lick %%%
        SS_lick_all(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_SS_all_corr * 100), 15, 'sgolay', 2);
        change_SS_lick_all(counter_pCell,:) = SS_lick_all(counter_pCell,:) - SS_str_bout_baseline_all(counter_pCell, 1);
        SS_lick_r(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_SS_r_corr * 100), 15, 'sgolay', 2);
        change_SS_lick_r(counter_pCell,:) = SS_lick_r(counter_pCell,:) - SS_str_bout_baseline_r(counter_pCell, 1);
        SS_lick_g(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_SS_grooming_corr * 100), 15, 'sgolay', 2);
        change_SS_lick_g(counter_pCell,:) = SS_lick_g(counter_pCell,:) - SS_str_bout_baseline_all(counter_pCell, 1);
        SS_lick_l(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_SS_l_corr * 100), 15, 'sgolay', 2);
        change_SS_lick_l(counter_pCell,:) = SS_lick_l(counter_pCell,:) - SS_str_bout_baseline_l(counter_pCell, 1);
        
    end
end

%%% Mean of d_tip_lick %%%
mean_d_tip_all = nanmean(d_tip_all);
se_d_tip_all = nanstd(d_tip_all)/sqrt(number_pCell);
mean_d_tip_l = nanmean(d_tip_l);
se_d_tip_l = nanstd(d_tip_l)/sqrt(number_pCell);
mean_d_tip_g = nanmean(d_tip_g);
se_d_tip_g = nanstd(d_tip_g)/sqrt(number_pCell);
mean_d_tip_r = nanmean(d_tip_r);
se_d_tip_r = nanstd(d_tip_r)/sqrt(number_pCell);

%%% Power spectrum of d_tip_lick %%%
% cutoff_freq = 0.5;
% [b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
% EPHYS_EB_xcorr_state_combined_100_filt  = filtfilt(b_butter,a_butter,EPHYS_EB_xcorr_state_combined_100);
Fs = 100;
N = length(mean_d_tip_all);
freq = 0:Fs/length(mean_d_tip_all):Fs/10;

for counter_pCell = 1 : number_pCell
    xdft_d_tip_all = fft(d_tip_all(counter_pCell, :));
    xdft_d_tip_all = xdft_d_tip_all(1:N/10+1);
    psdx_d_tip_all = (1/(Fs*N)) * abs(xdft_d_tip_all).^2;
    psdx_d_tip_all(2:end-1) = 2*psdx_d_tip_all(2:end-1);
    power_d_tip_all(counter_pCell,:) = 10*log10(psdx_d_tip_all);
    [peak_power_d_tip_all(counter_pCell,1), ind_peak_power_d_tip_all(counter_pCell,1)] = max(power_d_tip_all(counter_pCell,2:end));
    
    xdft_d_tip_r = fft(d_tip_r(counter_pCell,:));
    xdft_d_tip_r = xdft_d_tip_r(1:N/10+1);
    psdx_d_tip_r = (1/(Fs*N)) * abs(xdft_d_tip_r).^2;
    psdx_d_tip_r(2:end-1) = 2*psdx_d_tip_r(2:end-1);
    power_d_tip_r(counter_pCell,:) = 10*log10(psdx_d_tip_r);
    [peak_power_d_tip_r(counter_pCell,1), ind_peak_power_d_tip_r(counter_pCell,1)] = max(power_d_tip_r(counter_pCell,2:end));
    
    xdft_d_tip_g = fft(d_tip_g(counter_pCell,:));
    xdft_d_tip_g = xdft_d_tip_g(1:N/10+1);
    psdx_d_tip_g = (1/(Fs*N)) * abs(xdft_d_tip_g).^2;
    psdx_d_tip_g(2:end-1) = 2*psdx_d_tip_g(2:end-1);
    power_d_tip_g(counter_pCell,:) = 10*log10(psdx_d_tip_g);
    [peak_power_d_tip_g(counter_pCell,1), ind_peak_power_d_tip_g(counter_pCell,1)] = max(power_d_tip_g(counter_pCell,2:end));
    
    xdft_d_tip_l = fft(d_tip_l(counter_pCell,:));
    xdft_d_tip_l = xdft_d_tip_l(1:N/10+1);
    psdx_d_tip_l = (1/(Fs*N)) * abs(xdft_d_tip_l).^2;
    psdx_d_tip_l(2:end-1) = 2*psdx_d_tip_l(2:end-1);
    power_d_tip_l(counter_pCell,:) = 10*log10(psdx_d_tip_l);
    [peak_power_d_tip_l(counter_pCell,1), ind_peak_power_d_tip_l(counter_pCell,1)] = max(power_d_tip_l(counter_pCell,2:end));
end

%%% Mean of power_d_tip %%%
mean_power_d_tip_all = nanmean(power_d_tip_all);
se_power_d_tip_all = nanstd(power_d_tip_all)/sqrt(number_pCell);
mean_power_d_tip_r = nanmean(power_d_tip_r);
se_power_d_tip_r = nanstd(power_d_tip_r)/sqrt(number_pCell);
mean_power_d_tip_g = nanmean(power_d_tip_g);
se_power_d_tip_g = nanstd(power_d_tip_g)/sqrt(number_pCell);
mean_power_d_tip_l = nanmean(power_d_tip_l);
se_power_d_tip_l = nanstd(power_d_tip_l)/sqrt(number_pCell);

%%% Mean of change_SS_lick %%%
mean_change_SS_lick_all = nanmean(change_SS_lick_all);
se_change_SS_lick_all= nanstd(change_SS_lick_all)/sqrt(number_pCell);
mean_change_SS_lick_l = nanmean(change_SS_lick_l);
se_change_SS_lick_l= nanstd(change_SS_lick_l)/sqrt(number_pCell);
mean_change_SS_lick_g = nanmean(change_SS_lick_g);
se_change_SS_lick_g = nanstd(change_SS_lick_g)/sqrt(number_pCell);
mean_change_SS_lick_r = nanmean(change_SS_lick_r);
se_change_SS_lick_r = nanstd(change_SS_lick_r)/sqrt(number_pCell);

%%% Power spectrum of change_SS_lick %%%
Fs = 100;
N = length(mean_change_SS_lick_all);
freq = 0:Fs/length(mean_change_SS_lick_all):Fs/10;

for counter_pCell = 1 : number_pCell
    xdft_SS_all = fft(change_SS_lick_all(counter_pCell,:));
    xdft_SS_all = xdft_SS_all(1:N/10+1);
    psdx_SS_all = (1/(Fs*N)) * abs(xdft_SS_all).^2;
    psdx_SS_all(2:end-1) = 2*psdx_SS_all(2:end-1);
    power_SS_all(counter_pCell,:) = 10*log10(psdx_SS_all);
    [peak_power_SS_all(counter_pCell,1), ind_peak_power_SS_all(counter_pCell,1)] = max(power_SS_all(counter_pCell,2:end));
    
    xdft_SS_r = fft(change_SS_lick_r(counter_pCell,:));
    xdft_SS_r = xdft_SS_r(1:N/10+1);
    psdx_SS_r = (1/(Fs*N)) * abs(xdft_SS_r).^2;
    psdx_SS_r(2:end-1) = 2*psdx_SS_r(2:end-1);
    power_SS_r(counter_pCell,:) = 10*log10(psdx_SS_r);
    [peak_power_SS_r(counter_pCell,1), ind_peak_power_SS_r(counter_pCell,1)] = max(power_SS_r(counter_pCell,2:end));
    
    xdft_SS_g = fft(change_SS_lick_g(counter_pCell,:));
    xdft_SS_g = xdft_SS_g(1:N/10+1);
    psdx_SS_g = (1/(Fs*N)) * abs(xdft_SS_g).^2;
    psdx_SS_g(2:end-1) = 2*psdx_SS_g(2:end-1);
    power_SS_g(counter_pCell,:) = 10*log10(psdx_SS_g);
    [peak_power_SS_g(counter_pCell,1), ind_peak_power_SS_g(counter_pCell,1)] = max(power_SS_g(counter_pCell,2:end));
    
    xdft_SS_l = fft(change_SS_lick_l(counter_pCell,:));
    xdft_SS_l = xdft_SS_l(1:N/10+1);
    psdx_SS_l = (1/(Fs*N)) * abs(xdft_SS_l).^2;
    psdx_SS_l(2:end-1) = 2*psdx_SS_l(2:end-1);
    power_SS_l(counter_pCell,:) = 10*log10(psdx_SS_l);
    [peak_power_SS_l(counter_pCell,1), ind_peak_power_SS_l(counter_pCell,1)] = max(power_SS_l(counter_pCell,2:end));
end

%%% Mean of power_d_tip %%%
mean_power_SS_all = nanmean(power_SS_all);
se_power_SS_all = nanstd(power_SS_all)/sqrt(number_pCell);
mean_power_SS_r = nanmean(power_SS_r);
se_power_SS_r = nanstd(power_SS_r)/sqrt(number_pCell);
mean_power_SS_g = nanmean(power_SS_g);
se_power_SS_g = nanstd(power_SS_g)/sqrt(number_pCell);
mean_power_SS_l = nanmean(power_SS_l);
se_power_SS_l = nanstd(power_SS_l)/sqrt(number_pCell);

%%% Phase diff of d_tip_lick X change_SS_lick %%%
for counter_pCell = 1 : number_pCell
    phase_diff_all(counter_pCell,1) = wrapTo180(rad2deg(phdiffmeasure(d_tip_all(counter_pCell,:),change_SS_lick_all(counter_pCell,:))));
    phase_diff_r(counter_pCell,1) = wrapTo180(rad2deg(phdiffmeasure(d_tip_r(counter_pCell,:),change_SS_lick_r(counter_pCell,:))));
    phase_diff_g(counter_pCell,1) = wrapTo180(rad2deg(phdiffmeasure(d_tip_g(counter_pCell,:),change_SS_lick_g(counter_pCell,:))));
    phase_diff_l(counter_pCell,1) = wrapTo180(rad2deg(phdiffmeasure(d_tip_l(counter_pCell,:),change_SS_lick_l(counter_pCell,:))));
end

%%% Time diff of d_tip X change_SS_lick using xcorr %%%
for counter_pCell = 1 : number_pCell
    [xcorr_value,xcorr_lag] = xcorr(d_tip_all(counter_pCell,:), change_SS_lick_all(counter_pCell,:));
    [~,ind_max_xcross] = max(abs(xcorr_value));
    time_diff_all(counter_pCell,1) = xcorr_lag(ind_max_xcross) * 10;
    
    [xcorr_value,xcorr_lag] = xcorr(d_tip_l(counter_pCell,:), change_SS_lick_l(counter_pCell,:)); % cross-correlate signals with each other
    [~,ind_max_xcross] = max(abs(xcorr_value));
    time_diff_l(counter_pCell,1) = xcorr_lag(ind_max_xcross) * 10;
    
    [xcorr_value,xcorr_lag] = xcorr( d_tip_g(counter_pCell,:), change_SS_lick_g(counter_pCell,:)); % cross-correlate signals with each other
    [~,ind_max_xcross] = max(abs(xcorr_value));
    time_diff_g(counter_pCell,1) = xcorr_lag(ind_max_xcross) * 10;
    
    [xcorr_value,xcorr_lag] = xcorr(d_tip_r(counter_pCell,:), change_SS_lick_r(counter_pCell,:)); % cross-correlate signals with each other
    [~,ind_max_xcross] = max(abs(xcorr_value));
    time_diff_r(counter_pCell,1) = xcorr_lag(ind_max_xcross) *10;
end

%%% peak_power_SS_all - peak_power_d_tip_all %%%
diff_peak_power_all = (((freq(ind_peak_power_SS_all) + 1) - (freq(ind_peak_power_d_tip_all) + 1)))' ;
diff_peak_power_l = (((freq(ind_peak_power_SS_l) + 1) - (freq(ind_peak_power_d_tip_l) + 1)))' ;
diff_peak_power_g = (((freq(ind_peak_power_SS_g) + 1) - (freq(ind_peak_power_d_tip_g) + 1)))';
diff_peak_power_r = (((freq(ind_peak_power_SS_r) + 1) - (freq(ind_peak_power_d_tip_r) + 1)))' ;

%%% rhymicity labels %%%
label_rhythmic_all = (diff_peak_power_all >= -0.5 & diff_peak_power_all <= 0.5) ;
label_weakrhythmic_all = ~label_rhythmic_all;
label_rhythmic_r = diff_peak_power_r >= -0.5 & diff_peak_power_r <= 0.5 ;
label_weakrhythmic_r = ~label_rhythmic_r;
label_rhythmic_g = diff_peak_power_g >= -0.5 & diff_peak_power_g <= 0.5;
label_weakrhythmic_g = ~label_rhythmic_g;
label_rhythmic_l = diff_peak_power_l >= -0.5 & diff_peak_power_l <= 0.5 ;
label_weakrhythmic_l = ~label_rhythmic_l;

%%% phase labels &&&
label_phase_all = phase_diff_all >= -30 & phase_diff_all <= 30;
label_antiphase_all = phase_diff_all > 150 | phase_diff_all < -150;
label_lead_all = phase_diff_all > 30 & phase_diff_all <= 150;
label_lead30to90_all = phase_diff_all > 30 & phase_diff_all <= 90;
label_lead90to150_all = phase_diff_all > 90 & phase_diff_all <= 150;
label_lag_all = phase_diff_all <-30  & phase_diff_all >= -150;
label_lag30to90_all = phase_diff_all <-30 & phase_diff_all >= -90;
label_lag90to150_all = phase_diff_all <-90 & phase_diff_all >= -150;

label_phase_l = phase_diff_l >= -30 & phase_diff_l <= 30;
label_antiphase_l = phase_diff_l > 150 | phase_diff_l < -150;
label_lead_l = phase_diff_l > 30 & phase_diff_l <= 150;
label_lead30to90_l = phase_diff_l > 30 & phase_diff_l <= 90;
label_lead90to150_l = phase_diff_l > 90 & phase_diff_l <= 150;
label_lag_l = phase_diff_l <-30 & phase_diff_l >= -150;
label_lag30to90_l = phase_diff_l <-30 & phase_diff_l >= -90;
label_lag90to150_l = phase_diff_l <-90 & phase_diff_l >= -150;

label_phase_g = phase_diff_g >= -30 & phase_diff_g <= 30;
label_antiphase_g = phase_diff_g > 150 | phase_diff_g < -150;
label_lead_g = phase_diff_g > 30 & phase_diff_g <= 150;
label_lead30to90_g = phase_diff_g > 30 & phase_diff_g <= 90;
label_lead90to150_g = phase_diff_g > 90 & phase_diff_g <= 150;
label_lag_g = phase_diff_g <-30 & phase_diff_g >= -150;
label_lag30to90_g = phase_diff_g <-30 & phase_diff_g >= -90;
label_lag90to150_g = phase_diff_g <-90 & phase_diff_g >= -150;

label_phase_r = phase_diff_r >= -30 & phase_diff_r <= 30;
label_antiphase_r = phase_diff_r > 150 | phase_diff_r < -150;
label_lead_r = phase_diff_r > 30 & phase_diff_r <= 150;
label_lead30to90_r = phase_diff_r > 30 & phase_diff_r <= 90;
label_lead90to150_r = phase_diff_r > 90 & phase_diff_r <= 150;
label_lag_r = phase_diff_r <-30 & phase_diff_r >= -150;
label_lag30to90_r = phase_diff_r <-30 & phase_diff_r >= -90;
label_lag90to150_r = phase_diff_r <-90 & phase_diff_r >= -150;

%%% Mean d_tip based on labels
mean_d_tip_rhythmic_all = nanmean(d_tip_all(label_rhythmic_all,:));
se_d_tip_rhythmic_all = nanstd(d_tip_all(label_rhythmic_all,:))/sqrt(sum(label_rhythmic_all));
mean_d_tip_weakrhythmic_all = nanmean(d_tip_all(label_weakrhythmic_all,:));
se_d_tip_weakrhythmic_all = nanstd(d_tip_all(label_weakrhythmic_all,:))/sqrt(sum(label_weakrhythmic_all));
mean_d_tip_phase_all = nanmean(d_tip_all(label_phase_all,:));
se_d_tip_phase_all = nanstd(d_tip_all(label_phase_all,:))/sqrt(sum(label_phase_all));
mean_d_tip_antiphase_all = nanmean(d_tip_all(label_antiphase_all,:));
se_d_tip_antiphase_all = nanstd(d_tip_all(label_antiphase_all,:))/sqrt(sum(label_antiphase_all));
mean_d_tip_lead_all = nanmean(d_tip_all(label_lead_all,:));
se_d_tip_lead_all = nanstd(d_tip_all(label_lead_all,:))/sqrt(sum(label_lead_all));
mean_d_tip_lag_all = nanmean(d_tip_all(label_lag_all,:));
se_d_tip_lag_all = nanstd(d_tip_all(label_lag_all,:))/sqrt(sum(label_lag_all));
mean_d_tip_rhythmic_l = nanmean(d_tip_l(label_rhythmic_l,:));
se_d_tip_rhythmic_l = nanstd(d_tip_l(label_rhythmic_l,:))/sqrt(sum(label_rhythmic_l));
mean_d_tip_weakrhythmic_l = nanmean(d_tip_l(label_weakrhythmic_l,:));
se_d_tip_weakrhythmic_l = nanstd(d_tip_l(label_weakrhythmic_l,:))/sqrt(sum(label_weakrhythmic_l));
mean_d_tip_phase_l = nanmean(d_tip_l(label_phase_l,:));
se_d_tip_phase_l = nanstd(d_tip_l(label_phase_l,:))/sqrt(sum(label_phase_l));
mean_d_tip_antiphase_l = nanmean(d_tip_l(label_antiphase_l,:));
se_d_tip_antiphase_l = nanstd(d_tip_l(label_antiphase_l,:))/sqrt(sum(label_antiphase_l));
mean_d_tip_lead_l = nanmean(d_tip_l(label_lead_l,:));
se_d_tip_lead_l = nanstd(d_tip_l(label_lead_l,:))/sqrt(sum(label_lead_l));
mean_d_tip_lag_l = nanmean(d_tip_l(label_lag_l,:));
se_d_tip_lag_l = nanstd(d_tip_l(label_lag_l,:))/sqrt(sum(label_lag_l));
mean_d_tip_rhythmic_g = nanmean(d_tip_g(label_rhythmic_g,:));
se_d_tip_rhythmic_g = nanstd(d_tip_g(label_rhythmic_g,:))/sqrt(sum(label_rhythmic_g));
mean_d_tip_weakrhythmic_g = nanmean(d_tip_g(label_weakrhythmic_g,:));
se_d_tip_weakrhythmic_g = nanstd(d_tip_g(label_weakrhythmic_g,:))/sqrt(sum(label_weakrhythmic_g));
mean_d_tip_phase_g = nanmean(d_tip_g(label_phase_g,:));
se_d_tip_phase_g = nanstd(d_tip_g(label_phase_g,:))/sqrt(sum(label_phase_g));
mean_d_tip_antiphase_g = nanmean(d_tip_g(label_antiphase_g,:));
se_d_tip_antiphase_g = nanstd(d_tip_g(label_antiphase_g,:))/sqrt(sum(label_antiphase_g));
mean_d_tip_lead_g = nanmean(d_tip_g(label_lead_g,:));
se_d_tip_lead_g = nanstd(d_tip_g(label_lead_g,:))/sqrt(sum(label_lead_g));
mean_d_tip_lag_g = nanmean(d_tip_g(label_lag_g,:));
se_d_tip_lag_g = nanstd(d_tip_g(label_lag_g,:))/sqrt(sum(label_lag_g));
mean_d_tip_rhythmic_r = nanmean(d_tip_r(label_rhythmic_r,:));
se_d_tip_rhythmic_r = nanstd(d_tip_r(label_rhythmic_r,:))/sqrt(sum(label_rhythmic_r));
mean_d_tip_weakrhythmic_r = nanmean(d_tip_r(label_weakrhythmic_r,:));
se_d_tip_weakrhythmic_r = nanstd(d_tip_r(label_weakrhythmic_r,:))/sqrt(sum(label_weakrhythmic_r));
mean_d_tip_phase_r = nanmean(d_tip_r(label_phase_r,:));
se_d_tip_phase_r = nanstd(d_tip_r(label_phase_r,:))/sqrt(sum(label_phase_r));
mean_d_tip_antiphase_r = nanmean(d_tip_r(label_antiphase_r,:));
se_d_tip_antiphase_r = nanstd(d_tip_r(label_antiphase_r,:))/sqrt(sum(label_antiphase_r));
mean_d_tip_lead_r = nanmean(d_tip_r(label_lead_r,:));
se_d_tip_lead_r = nanstd(d_tip_r(label_lead_r,:))/sqrt(sum(label_lead_r));
mean_d_tip_lag_r = nanmean(d_tip_r(label_lag_r,:));
se_d_tip_lag_r = nanstd(d_tip_r(label_lag_r,:))/sqrt(sum(label_lag_r));

%%% Mean change_SS_lick based on labels
mean_change_SS_lick_rhythmic_all = nanmean(change_SS_lick_all(label_rhythmic_all,:));
se_change_SS_lick_rhythmic_all = nanstd(change_SS_lick_all(label_rhythmic_all,:))/sqrt(sum(label_rhythmic_all));
mean_change_SS_lick_weakrhythmic_all = nanmean(change_SS_lick_all(label_weakrhythmic_all,:));
se_change_SS_lick_weakrhythmic_all = nanstd(change_SS_lick_all(label_weakrhythmic_all,:))/sqrt(sum(label_weakrhythmic_all));
mean_change_SS_lick_phase_all = nanmean(change_SS_lick_all(label_phase_all,:));
se_change_SS_lick_phase_all = nanstd(change_SS_lick_all(label_phase_all,:))/sqrt(sum(label_phase_all));
mean_change_SS_lick_antiphase_all = nanmean(change_SS_lick_all(label_antiphase_all,:));
se_change_SS_lick_antiphase_all = nanstd(change_SS_lick_all(label_antiphase_all,:))/sqrt(sum(label_antiphase_all));
mean_change_SS_lick_lead_all = nanmean(change_SS_lick_all(label_lead_all,:));
se_change_SS_lick_lead_all = nanstd(change_SS_lick_all(label_lead_all,:))/sqrt(sum(label_lead_all));
mean_change_SS_lick_lag_all = nanmean(change_SS_lick_all(label_lag_all,:));
se_change_SS_lick_lag_all = nanstd(change_SS_lick_all(label_lag_all,:))/sqrt(sum(label_lag_all));
mean_change_SS_lick_rhythmic_l = nanmean(change_SS_lick_l(label_rhythmic_l,:));
se_change_SS_lick_rhythmic_l = nanstd(change_SS_lick_l(label_rhythmic_l,:))/sqrt(sum(label_rhythmic_l));
mean_change_SS_lick_weakrhythmic_l = nanmean(change_SS_lick_l(label_weakrhythmic_l,:));
se_change_SS_lick_weakrhythmic_l = nanstd(change_SS_lick_l(label_weakrhythmic_l,:))/sqrt(sum(label_weakrhythmic_l));
mean_change_SS_lick_phase_l = nanmean(change_SS_lick_l(label_phase_l,:));
se_change_SS_lick_phase_l = nanstd(change_SS_lick_l(label_phase_l,:))/sqrt(sum(label_phase_l));
mean_change_SS_lick_antiphase_l = nanmean(change_SS_lick_l(label_antiphase_l,:));
se_change_SS_lick_antiphase_l = nanstd(change_SS_lick_l(label_antiphase_l,:))/sqrt(sum(label_antiphase_l));
mean_change_SS_lick_lead_l = nanmean(change_SS_lick_l(label_lead_l,:));
se_change_SS_lick_lead_l = nanstd(change_SS_lick_l(label_lead_l,:))/sqrt(sum(label_lead_l));
mean_change_SS_lick_lag_l = nanmean(change_SS_lick_l(label_lag_l,:));
se_change_SS_lick_lag_l = nanstd(change_SS_lick_l(label_lag_l,:))/sqrt(sum(label_lag_l));
mean_change_SS_lick_rhythmic_g = nanmean(change_SS_lick_g(label_rhythmic_g,:));
se_change_SS_lick_rhythmic_g = nanstd(change_SS_lick_g(label_rhythmic_g,:))/sqrt(sum(label_rhythmic_g));
mean_change_SS_lick_weakrhythmic_g = nanmean(change_SS_lick_g(label_weakrhythmic_g,:));
se_change_SS_lick_weakrhythmic_g = nanstd(change_SS_lick_g(label_weakrhythmic_g,:))/sqrt(sum(label_weakrhythmic_g));
mean_change_SS_lick_phase_g = nanmean(change_SS_lick_g(label_phase_g,:));
se_change_SS_lick_phase_g = nanstd(change_SS_lick_g(label_phase_g,:))/sqrt(sum(label_phase_g));
mean_change_SS_lick_antiphase_g = nanmean(change_SS_lick_g(label_antiphase_g,:));
se_change_SS_lick_antiphase_g = nanstd(change_SS_lick_g(label_antiphase_g,:))/sqrt(sum(label_antiphase_g));
mean_change_SS_lick_lead_g = nanmean(change_SS_lick_g(label_lead_g,:));
se_change_SS_lick_lead_g = nanstd(change_SS_lick_g(label_lead_g,:))/sqrt(sum(label_lead_g));
mean_change_SS_lick_lag_g = nanmean(change_SS_lick_g(label_lag_g,:));
se_change_SS_lick_lag_g = nanstd(change_SS_lick_g(label_lag_g,:))/sqrt(sum(label_lag_g));
mean_change_SS_lick_rhythmic_r = nanmean(change_SS_lick_r(label_rhythmic_r,:));
se_change_SS_lick_rhythmic_r = nanstd(change_SS_lick_r(label_rhythmic_r,:))/sqrt(sum(label_rhythmic_r));
mean_change_SS_lick_weakrhythmic_r = nanmean(change_SS_lick_r(label_weakrhythmic_r,:));
se_change_SS_lick_weakrhythmic_r = nanstd(change_SS_lick_r(label_weakrhythmic_r,:))/sqrt(sum(label_weakrhythmic_r));
mean_change_SS_lick_phase_r = nanmean(change_SS_lick_r(label_phase_r,:));
se_change_SS_lick_phase_r = nanstd(change_SS_lick_r(label_phase_r,:))/sqrt(sum(label_phase_r));
mean_change_SS_lick_antiphase_r = nanmean(change_SS_lick_r(label_antiphase_r,:));
se_change_SS_lick_antiphase_r = nanstd(change_SS_lick_r(label_antiphase_r,:))/sqrt(sum(label_antiphase_r));
mean_change_SS_lick_lead_r = nanmean(change_SS_lick_r(label_lead_r,:));
se_change_SS_lick_lead_r = nanstd(change_SS_lick_r(label_lead_r,:))/sqrt(sum(label_lead_r));
mean_change_SS_lick_lag_r = nanmean(change_SS_lick_r(label_lag_r,:));
se_change_SS_lick_lag_r = nanstd(change_SS_lick_r(label_lag_r,:))/sqrt(sum(label_lag_r));

%%% inds_span_corr %%%
inds_span_corr = -0.99 : 0.01 : 1;

%%% xaxis %%%
xaxis_pCell = 1:62;

%%% SAVE TO POPULATION STRUCTURE %%%
POPULATION.dmax.label_rhythmic_all = label_rhythmic_all;
POPULATION.dmax.label_weakrhythmic_all = label_weakrhythmic_all;
POPULATION.dmax.label_rhythmic_l = label_rhythmic_l;
POPULATION.dmax.label_weakrhythmic_l = label_weakrhythmic_l;
POPULATION.dmax.label_rhythmic_g = label_rhythmic_g;
POPULATION.dmax.label_weakrhythmic_g = label_weakrhythmic_g;
POPULATION.dmax.label_rhythmic_r = label_rhythmic_r;
POPULATION.dmax.label_weakrhythmic_r = label_weakrhythmic_r;

POPULATION.dmax.label_lead_all = label_lead_all;
POPULATION.dmax.label_lead30to90_all = label_lead30to90_all;
POPULATION.dmax.label_lead90to150_all = label_lead90to150_all;
POPULATION.dmax.label_lag_all = label_lag_all;
POPULATION.dmax.label_lag30to90_all =  label_lag30to90_all;
POPULATION.dmax.label_lag90to150_all = label_lag90to150_all;
POPULATION.dmax.label_phase_all = label_phase_all;
POPULATION.dmax.label_antiphase_all = label_antiphase_all;

POPULATION.dmax.label_lead_l = label_lead_l;
POPULATION.dmax.label_lead30to90_l = label_lead30to90_l;
POPULATION.dmax.label_lead90to150_l = label_lead90to150_l;
POPULATION.dmax.label_lag_l = label_lag_l;
POPULATION.dmax.label_lag30to90_l =  label_lag30to90_l;
POPULATION.dmax.label_lag90to150_l = label_lag90to150_l;
POPULATION.dmax.label_phase_l = label_phase_l;
POPULATION.dmax.label_antiphase_l = label_antiphase_l;

POPULATION.dmax.label_lead_g = label_lead_g;
POPULATION.dmax.label_lead30to90_g = label_lead30to90_g;
POPULATION.dmax.label_lead90to150_g = label_lead90to150_g;
POPULATION.dmax.label_lag_g = label_lag_g;
POPULATION.dmax.label_lag30to90_g =  label_lag30to90_g;
POPULATION.dmax.label_lag90to150_g = label_lag90to150_g;
POPULATION.dmax.label_phase_g = label_phase_g;
POPULATION.dmax.label_antiphase_g = label_antiphase_g;

POPULATION.dmax.label_lead_r = label_lead_r;
POPULATION.dmax.label_lead30to90_r = label_lead30to90_r;
POPULATION.dmax.label_lead90to150_r = label_lead90to150_r;
POPULATION.dmax.label_lag_r = label_lag_r;
POPULATION.dmax.label_lag30to90_r =  label_lag30to90_r;
POPULATION.dmax.label_lag90to150_r = label_lag90to150_r;
POPULATION.dmax.label_phase_r = label_phase_r;
POPULATION.dmax.label_antiphase_r = label_antiphase_r;
%% RAW
num_pCell = 10;

figure
subplot(4,4,1)
hold on
yyaxis left
for counter_pCell = num_pCell
    plot(inds_span_corr, d_tip_all(counter_pCell, :), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
ylim([0 18])
ylabel('Displacement (mm)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right
for counter_pCell = num_pCell
    plot(inds_span_corr, change_SS_lick_all(counter_pCell, :), '-', 'Color', 'r', 'LineWidth', 2 )
end
ylim([(min(change_SS_lick_all(counter_pCell, :)) - 30) (max(change_SS_lick_all(counter_pCell, :)) + 30)])
ylabel('Change in SS firing (spks/s)')
set(gca, 'YColor', 'k')
xlabel('Lick dmax (s)')
title(['Displacement vs SS firing | All | pCell: ' num2str(counter_pCell)])

subplot(4,4,5)
hold on
for counter_pCell = num_pCell
    plot(freq, power_d_tip_all(counter_pCell, :), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 2 )
    plot(freq, power_SS_all(counter_pCell, :), '-', 'Color', [1 0 0], 'LineWidth', 2 )
end
xlabel('Freq (Hz)')
ylabel('Pow/Freq (dB/Hz)')
title(['Power spectrum | All | pCell: ' num2str(counter_pCell)])

subplot(4,4,9)
hold on
bar(freq(ind_peak_power_d_tip_all + 1), 'FaceColor', [0.7 0.7 0.7], 'Barwidth', 0.8)
bar(freq(ind_peak_power_SS_all + 1), 'k', 'Barwidth', 0.5)
bar(xaxis_pCell(num_pCell),freq(ind_peak_power_SS_all(num_pCell) + 1), 'r', 'Barwidth', 0.5)
ylim([0 7])
xlabel('pCell')
ylabel('Freq of peak power')

subplot(4,4,13)
hold on
bar(phase_diff_all, 'FaceColor','k', 'Barwidth', 0.8)
bar(xaxis_pCell(num_pCell), phase_diff_all(num_pCell), 'FaceColor','r', 'Barwidth', 0.8)
ylim([-190 190])
xlabel('pCell')
ylabel('Phase difference (deg)')

%     subplot(5,4,17)
%     bar(time_diff_all, 'FaceColor','k', 'Barwidth', 0.5)
% %     ylim([-360 360])
%     xlabel('pCell')
%     ylabel('Time difference (ms)')

subplot(4,4,2)
hold on
yyaxis left
for counter_pCell = num_pCell
    plot(inds_span_corr, d_tip_l(counter_pCell, :), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
ylim([0  18])
ylabel('Displacement (mm)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right
for counter_pCell = num_pCell
    plot(inds_span_corr, change_SS_lick_l(counter_pCell, :), '-', 'Color', 'r', 'LineWidth', 2 )
end
ylim([(min(change_SS_lick_all(counter_pCell, :)) - 30) (max(change_SS_lick_all(counter_pCell, :)) + 30)])
ylabel('Change in SS firing (spks/s)')
set(gca, 'YColor', 'k')
xlabel('Lick dmax (ms)')
title(['Displacement vs SS firing | L | pCell: ' num2str(counter_pCell)])

subplot(4,4,6)
hold on
for counter_pCell = num_pCell
    plot(freq, power_d_tip_l(counter_pCell, :), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 2 )
    plot(freq, power_SS_l(counter_pCell, :), '-', 'Color', [1 0 0], 'LineWidth', 2 )
end

xlabel('Freq (Hz)')
ylabel('Pow/Freq (dB/Hz)')
title(['Power spectrum | L | pCell: ' num2str(counter_pCell)])

subplot(4,4,10)
hold on
bar(freq(ind_peak_power_d_tip_l + 1), 'FaceColor', [0.7 0.7 0.7], 'Barwidth', 0.8)
bar(freq(ind_peak_power_SS_l + 1), 'k', 'Barwidth', 0.5)
bar(xaxis_pCell(num_pCell), freq(ind_peak_power_SS_l(num_pCell) + 1), 'r', 'Barwidth', 0.5)
ylim([0 7])
xlabel('pCell')
ylabel('Freq of peak power')

subplot(4,4,14)
hold on
bar(phase_diff_l, 'FaceColor','k', 'Barwidth', 0.8)
bar(xaxis_pCell(num_pCell), phase_diff_l(num_pCell), 'FaceColor','r', 'Barwidth', 0.8)
ylim([-190 190])
xlabel('pCell')
ylabel('Phase difference (deg)')

%     subplot(5,4,18)
%     bar(time_diff_l, 'FaceColor','k', 'Barwidth', 0.5)
% %     ylim([-360 360])
%     xlabel('pCell')
%     ylabel('Time difference (ms)')

subplot(4,4,3)
hold on
yyaxis left
for counter_pCell = num_pCell
    plot(inds_span_corr, d_tip_g(counter_pCell, :), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
ylim([0 18])
ylabel('Displacement (mm)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right
for counter_pCell = num_pCell
    plot(inds_span_corr, change_SS_lick_g(counter_pCell, :), '-', 'Color', 'r', 'LineWidth', 2 )
end
ylim([(min(change_SS_lick_all(counter_pCell, :)) - 30) (max(change_SS_lick_all(counter_pCell, :)) + 30)])
ylabel('Change in SS firing (spks/s)')
set(gca, 'YColor', 'k')
xlabel('Lick dmax (s)')
title(['Displacement vs SS firing | G | pCell: ' num2str(counter_pCell)])

subplot(4,4,7)
hold on
for counter_pCell = num_pCell
    plot(freq, power_d_tip_g(counter_pCell, :), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 2 )
    plot(freq, power_SS_g(counter_pCell, :), '-', 'Color', [1 0 0], 'LineWidth', 2 )
end
xlabel('Freq (Hz)')
ylabel('Pow/Freq (dB/Hz)')
title(['Power spectrum | G | pCell: ' num2str(counter_pCell)])

subplot(4,4,11)
hold on
bar(freq(ind_peak_power_d_tip_g + 1), 'FaceColor', [0.7 0.7 0.7], 'Barwidth', 0.8)
bar(freq(ind_peak_power_SS_g + 1), 'k', 'Barwidth', 0.5)
bar(xaxis_pCell(num_pCell), freq(ind_peak_power_SS_g(num_pCell) + 1), 'r', 'Barwidth', 0.5)
ylim([0 7])
xlabel('pCell')
ylabel('Freq of peak power')

subplot(4,4,15)
hold on
bar(phase_diff_g, 'FaceColor','k', 'Barwidth', 0.8)
bar(xaxis_pCell(num_pCell), phase_diff_g(num_pCell), 'FaceColor','r', 'Barwidth', 0.8)
ylim([-190 190])
xlabel('pCell')
ylabel('Phase difference (deg)')

%     subplot(5,4,19)
%     bar(time_diff_g, 'FaceColor','k', 'Barwidth', 0.5)
% %     ylim([-360 360])
%     xlabel('pCell')
%     ylabel('Time difference (ms)')

subplot(4,4,4)
hold on
yyaxis left
for counter_pCell = num_pCell
    plot(inds_span_corr, d_tip_r(counter_pCell, :), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 2 )
end
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
ylim([0 18])
ylabel('Displacement (mm)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right
for counter_pCell = num_pCell
    plot(inds_span_corr, change_SS_lick_r(counter_pCell, :), '-', 'Color', 'r', 'LineWidth', 2 )
end
ylim([(min(change_SS_lick_all(counter_pCell, :)) - 30) (max(change_SS_lick_all(counter_pCell, :)) + 30)])
ylabel('Change in SS firing (spks/s)')
set(gca, 'YColor', 'k')
xlabel('Lick dmax (s)')
title(['Displacement vs SS firing | R | pCell: ' num2str(counter_pCell)])

subplot(4,4,8)
hold on
for counter_pCell = num_pCell
    plot(freq, power_d_tip_r(counter_pCell, :), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 2 )
    plot(freq, power_SS_r(counter_pCell, :), '-', 'Color', [1 0 0], 'LineWidth', 2 )
end
xlabel('Freq (Hz)')
ylabel('Pow/Freq (dB/Hz)')
title(['Power spectrum | R | pCell: ' num2str(counter_pCell)])

subplot(4,4,12)
hold on
bar(freq(ind_peak_power_d_tip_r + 1), 'FaceColor', [0.7 0.7 0.7], 'Barwidth', 0.8)
bar(freq(ind_peak_power_SS_r + 1), 'k', 'Barwidth', 0.5)
bar(xaxis_pCell(num_pCell), freq(ind_peak_power_SS_r(num_pCell) + 1), 'r', 'Barwidth', 0.5)
ylim([0 7])
xlabel('pCell')
ylabel('Freq of peak power')

subplot(4,4,16)
hold on
bar(phase_diff_r, 'FaceColor','k', 'Barwidth', 0.8)
bar(xaxis_pCell(num_pCell),phase_diff_r(num_pCell), 'FaceColor','r', 'Barwidth', 0.8)
ylim([-190 190])
xlabel('pCell')
ylabel('Phase difference (deg)')

%     subplot(5,4,20)
%     bar(time_diff_g, 'FaceColor','k', 'Barwidth', 0.5)
% %     ylim([-360 360])
%     xlabel('pCell')
%     ylabel('Time difference (ms)')
%% LABELED
xaxis_pCell = 1:62;

figure
subplot(4, 8, [1 2])
hold on
bar(xaxis_pCell(label_rhythmic_all), freq(ind_peak_power_SS_all(label_rhythmic_all)) + 1, 'FaceColor', 'r','Barwidth', 1)
bar(xaxis_pCell(label_weakrhythmic_all), freq(ind_peak_power_SS_all(label_weakrhythmic_all)) + 1, 'FaceColor', 'k','Barwidth', 1)
bar(freq(ind_peak_power_d_tip_all + 1), 'FaceColor', [0.7 0.7 0.7], 'Barwidth', 0.5)
ylabel('Frequency of peak power (Hz)')
title('Comparison of frequency: SS firing - lick | All licks')
ylim([0 7])

pie_1 = subplot(4,8,9);
pie_values = [sum(label_rhythmic_all) sum(label_weakrhythmic_all) ];
pie_labels = {num2str(sum(label_rhythmic_all)) ...
    num2str(sum(label_weakrhythmic_all))};
explode_1 = [1 1];
pie(pie_values, explode_1, pie_labels)
colormap(pie_1, [1 0 0; 0 0 0]);

subplot(4,8,10)
hold on;
plot(freq(ind_peak_power_d_tip_all(label_rhythmic_all)) + 1, freq(ind_peak_power_SS_all(label_rhythmic_all)) + 1, 'or');
plot(freq(ind_peak_power_d_tip_all(label_weakrhythmic_all)) + 1, freq(ind_peak_power_SS_all(label_weakrhythmic_all)) + 1, 'ok');
vert = refline(1,0);
vert.Color = 'k';
xlabel('Lick frequency (Hz)')
ylabel('SS firing frequency (Hz)')
xlim([0 7])
ylim([0 7])

subplot(4, 8, [3 4])
hold on
bar(xaxis_pCell(label_rhythmic_l), freq(ind_peak_power_SS_l(label_rhythmic_l)) + 1, 'FaceColor', 'r','Barwidth', 1)
bar(xaxis_pCell(label_weakrhythmic_l), freq(ind_peak_power_SS_l(label_weakrhythmic_l)) + 1, 'FaceColor', 'k','Barwidth', 1)
bar(freq(ind_peak_power_d_tip_l + 1), 'FaceColor', [0.7 0.7 0.7], 'Barwidth', 0.5)
ylabel('Frequency of peak power (Hz)')
title('Comparison of frequency: SS firing - lick | L licks')
ylim([0 7])

pie_2 = subplot(4,8,11);
pie_values = [sum(label_rhythmic_l) sum(label_weakrhythmic_l) ];
pie_labels = {num2str(sum(label_rhythmic_l)) ...
    num2str(sum(label_weakrhythmic_l))};
explode_2 = [1 1];
pie(pie_values, explode_1, pie_labels)
colormap(pie_2, [1 0 0; 0 0 0]);

subplot(4,8,12)
hold on;
plot(freq(ind_peak_power_d_tip_l(label_rhythmic_l)) + 1, freq(ind_peak_power_SS_l(label_rhythmic_l)) + 1, 'or');
plot(freq(ind_peak_power_d_tip_l(label_weakrhythmic_l)) + 1, freq(ind_peak_power_SS_l(label_weakrhythmic_l)) + 1, 'ok');
vert = refline(1,0);
vert.Color = 'k';
xlabel('Lick frequency (Hz)')
ylabel('SS firing frequency (Hz)')
xlim([0 7])
ylim([0 7])

subplot(4, 8, [5 6])
hold on
bar(xaxis_pCell(label_rhythmic_g), freq(ind_peak_power_SS_g(label_rhythmic_g)) + 1, 'FaceColor', 'r','Barwidth', 1)
bar(xaxis_pCell(label_weakrhythmic_g), freq(ind_peak_power_SS_g(label_weakrhythmic_g)) + 1, 'FaceColor', 'k','Barwidth', 1)
bar(freq(ind_peak_power_d_tip_g + 1), 'FaceColor', [0.7 0.7 0.7], 'Barwidth', 0.5)
ylabel('Frequency of peak power (Hz)')
title('Comparison of frequency: SS firing - lick | G licks')
ylim([0 7])

pie_3 = subplot(4,8,13);
pie_values = [sum(label_rhythmic_g) sum(label_weakrhythmic_g) ];
pie_labels = {num2str(sum(label_rhythmic_g)) ...
    num2str(sum(label_weakrhythmic_g))};
explode_3 = [1 1];
pie(pie_values, explode_3, pie_labels)
colormap(pie_3, [1 0 0; 0 0 0]);

subplot(4,8,14)
hold on;
plot(freq(ind_peak_power_d_tip_g(label_rhythmic_g)) + 1, freq(ind_peak_power_SS_g(label_rhythmic_g)) + 1, 'or');
plot(freq(ind_peak_power_d_tip_g(label_weakrhythmic_g)) + 1, freq(ind_peak_power_SS_g(label_weakrhythmic_g)) + 1, 'ok');
vert = refline(1,0);
vert.Color = 'k';
xlabel('Lick frequency (Hz)')
ylabel('SS firing frequency (Hz)')
xlim([0 7])
ylim([0 7])

subplot(4, 8, [7 8])
hold on
bar(xaxis_pCell(label_rhythmic_r), freq(ind_peak_power_SS_r(label_rhythmic_r)) + 1, 'FaceColor', 'r','Barwidth', 1)
bar(xaxis_pCell(label_weakrhythmic_r), freq(ind_peak_power_SS_r(label_weakrhythmic_r)) + 1, 'FaceColor', 'k','Barwidth', 1)
bar(freq(ind_peak_power_d_tip_r + 1), 'FaceColor', [0.7 0.7 0.7], 'Barwidth', 0.5)
ylabel('Frequency of peak power (Hz)')
title('Comparison of frequency: SS firing - lick | R licks')
ylim([0 7])

pie_4 = subplot(4,8,15);
pie_values = [sum(label_rhythmic_r) sum(label_weakrhythmic_r) ];
pie_labels = {num2str(sum(label_rhythmic_r)) ...
    num2str(sum(label_weakrhythmic_r))};
explode_4 = [1 1];
pie(pie_values, explode_4, pie_labels)
colormap(pie_4, [1 0 0; 0 0 0]);

subplot(4,8,16)
hold on;
plot(freq(ind_peak_power_d_tip_r(label_rhythmic_r)) + 1, freq(ind_peak_power_SS_r(label_rhythmic_r)) + 1, 'or');
plot(freq(ind_peak_power_d_tip_r(label_weakrhythmic_r)) + 1, freq(ind_peak_power_SS_r(label_weakrhythmic_r)) + 1, 'ok');
vert = refline(1,0);
vert.Color = 'k';
xlabel('Lick frequency (Hz)')
ylabel('SS firing frequency (Hz)')
xlim([0 7])
ylim([0 7])

subplot(4, 8, [17 18])
hold on
bar(xaxis_pCell(label_lead_all), phase_diff_all(label_lead_all), 'FaceColor', 'r','Barwidth', 1)
bar(xaxis_pCell(label_lag_all), phase_diff_all(label_lag_all), 'FaceColor', 'b','Barwidth', 1)
bar(xaxis_pCell(label_phase_all), phase_diff_all(label_phase_all), 'FaceColor', 'k','Barwidth', 1)
bar(xaxis_pCell(label_antiphase_all), phase_diff_all(label_antiphase_all), 'FaceColor', 'm','Barwidth', 0.5)
ylabel('Phase difference (deg)')
title('Phase difference: SS firing - lick | All licks')
ylim([-190 190])

pie_5 = subplot(4,8,[25 26]);
pie_values = [sum(label_lead_all) sum(label_lag_all) sum(label_phase_all) sum(label_antiphase_all)];
pie_labels = {num2str(sum(label_lead_all)) ...
    num2str(sum(label_lag_all)) ...
    num2str(sum(label_phase_all)) ...
    num2str(sum(label_antiphase_all)) };
explode_5 = [1 1 1 1];
pie(pie_values, explode_5, pie_labels)
colormap(pie_5, [1 0 0; 0 0 1; 0 0 0; 1 0 1]);

subplot(4, 8, [19 20])
hold on
bar(xaxis_pCell(label_lead_l), phase_diff_l(label_lead_l), 'FaceColor', 'r','Barwidth', 1)
bar(xaxis_pCell(label_lag_l), phase_diff_l(label_lag_l), 'FaceColor', 'b','Barwidth', 1)
bar(xaxis_pCell(label_phase_l), phase_diff_l(label_phase_l), 'FaceColor', 'k','Barwidth', 1)
bar(xaxis_pCell(label_antiphase_l), phase_diff_l(label_antiphase_l), 'FaceColor', 'm','Barwidth', 1)
ylabel('Phase difference (deg)')
title('Phase difference: SS firing - lick | L licks')
ylim([-190 190])

pie_6 = subplot(4,8,[27 28]);
pie_values = [sum(label_lead_l) sum(label_lag_l) sum(label_phase_l) sum(label_antiphase_l) ];
pie_labels = {num2str(sum(label_lead_l)) ...
    num2str(sum(label_lag_l)) ...
    num2str(sum(label_phase_l)) ...
    num2str(sum(label_antiphase_l))};
explode_6 = [1 1 1 1];
pie(pie_values, explode_6, pie_labels)
colormap(pie_6, [1 0 0; 0 0 1; 0 0 0; 1 0 1]);

subplot(4, 8, [21 22])
hold on
bar(xaxis_pCell(label_lead_g), phase_diff_g(label_lead_g), 'FaceColor', 'r','Barwidth', 1)
bar(xaxis_pCell(label_lag_g), phase_diff_g(label_lag_g), 'FaceColor', 'b','Barwidth', 1)
bar(xaxis_pCell(label_phase_g), phase_diff_g(label_phase_g), 'FaceColor', 'k','Barwidth', 0.5)
bar(xaxis_pCell(label_antiphase_g), phase_diff_g(label_antiphase_g), 'FaceColor', 'm','Barwidth', 1)
ylabel('Phase difference (deg)')
title('Phase difference: SS firing - lick | G licks')
ylim([-190 190])

pie_7 = subplot(4,8,[29 30]);
pie_values = [sum(label_lead_g) sum(label_lag_g) sum(label_phase_g) sum(label_antiphase_g)];
pie_labels = {num2str(sum(label_lead_g)) ...
    num2str(sum(label_lag_g)) ...
    num2str(sum(label_phase_g)) ...
    num2str(sum(label_antiphase_g))};
explode_7 = [1 1 1 1];
pie(pie_values, explode_7, pie_labels)
colormap(pie_7, [1 0 0; 0 0 1; 0 0 0; 1 0 1]);

subplot(4, 8, [23 24])
hold on
bar(xaxis_pCell(label_lead_r), phase_diff_r(label_lead_r), 'FaceColor', 'r','Barwidth', 1)
bar(xaxis_pCell(label_lag_r), phase_diff_r(label_lag_r), 'FaceColor', 'b','Barwidth', 1)
bar(xaxis_pCell(label_phase_r), phase_diff_r(label_phase_r), 'FaceColor', 'k','Barwidth', 0.5)
bar(xaxis_pCell(label_antiphase_r), phase_diff_r(label_antiphase_r), 'FaceColor', 'm','Barwidth', 0.5)
ylabel('Phase difference (deg)')
title('Phase difference: SS firing - lick | R licks')
ylim([-190 190])

pie_8 = subplot(4,8,[31 32]);
pie_values = [sum(label_lead_r) sum(label_lag_r) sum(label_phase_r)  sum(label_antiphase_r) ];
pie_labels = {num2str(sum(label_lead_r)) ...
    num2str(sum(label_lag_r)) ...
    num2str(sum(label_phase_r)) ...
    num2str(sum(label_antiphase_r))};
explode_8 = [1 1 1 1 ];
pie(pie_values, explode_8, pie_labels)
colormap(pie_8, [1 0 0; 0 0 1; 0 0 0; 1 0 1]);
%% LABELED TRACES
%%% TRACES %%%
figure
sgtitle(['Lick: freq & phase | pCells: ' num2str(number_pCell)])

subplot(4,4,1)
hold on
plot(inds_span_corr, mean_d_tip_rhythmic_all(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_all + se_d_tip_rhythmic_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_all - se_d_tip_rhythmic_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_weakrhythmic_all(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_weakrhythmic_all + se_d_tip_weakrhythmic_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_weakrhythmic_all - se_d_tip_weakrhythmic_all, 'r', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-1 1])
ylabel('Displacement (mm)')
xlabel('Lick dmax (ms)')
title('Lick trace | All licks')

subplot(4,4,2)
hold on
plot(inds_span_corr, mean_d_tip_rhythmic_l(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_l + se_d_tip_rhythmic_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_l - se_d_tip_rhythmic_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_weakrhythmic_l(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_weakrhythmic_l + se_d_tip_weakrhythmic_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_weakrhythmic_l - se_d_tip_weakrhythmic_l, 'r', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-1 1])
ylabel('Displacement (mm)')
xlabel('Lick dmax (ms)')
title('Lick trace | L licks')

subplot(4,4,3)
hold on
plot(inds_span_corr, mean_d_tip_rhythmic_g(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_g + se_d_tip_rhythmic_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_g - se_d_tip_rhythmic_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_weakrhythmic_g(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_weakrhythmic_g + se_d_tip_weakrhythmic_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_weakrhythmic_g - se_d_tip_weakrhythmic_g, 'r', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-1 1])
ylabel('Displacement (mm)')
xlabel('Lick dmax (ms)')
title('Lick trace | G licks')

subplot(4,4,4)
hold on
plot(inds_span_corr, mean_d_tip_rhythmic_r(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_r + se_d_tip_rhythmic_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_r - se_d_tip_rhythmic_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_weakrhythmic_r(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_weakrhythmic_r + se_d_tip_weakrhythmic_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_weakrhythmic_r - se_d_tip_weakrhythmic_r, 'r', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-1 1])
ylabel('Displacement (mm)')
xlabel('Lick dmax (ms)')
title('Lick trace | R licks')

subplot(4,4,5)
hold on
plot(inds_span_corr, mean_change_SS_lick_rhythmic_all(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_rhythmic_all + se_change_SS_lick_rhythmic_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_rhythmic_all - se_change_SS_lick_rhythmic_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_change_SS_lick_weakrhythmic_all(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_weakrhythmic_all + se_change_SS_lick_weakrhythmic_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_weakrhythmic_all - se_change_SS_lick_weakrhythmic_all, 'r', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-30 120])
xlim([-1 1])
ylabel('Change in SS firing (spks/s)')
xlabel('Lick dmax (ms)')
title('SS firing | All licks')

subplot(4,4,6)
hold on
plot(inds_span_corr, mean_change_SS_lick_rhythmic_l(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_rhythmic_l + se_change_SS_lick_rhythmic_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_rhythmic_l - se_change_SS_lick_rhythmic_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_change_SS_lick_weakrhythmic_l(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_weakrhythmic_l + se_change_SS_lick_weakrhythmic_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_weakrhythmic_l - se_change_SS_lick_weakrhythmic_l, 'r', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-30 120])
xlim([-1 1])
ylabel('Change in SS firing (spks/s)')
xlabel('Lick dmax (ms)')
title('SS firing | L licks')

subplot(4,4,7)
hold on
plot(inds_span_corr, mean_change_SS_lick_rhythmic_g(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_rhythmic_g + se_change_SS_lick_rhythmic_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_rhythmic_g - se_change_SS_lick_rhythmic_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_change_SS_lick_weakrhythmic_g(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_weakrhythmic_g + se_change_SS_lick_weakrhythmic_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_weakrhythmic_g - se_change_SS_lick_weakrhythmic_g, 'r', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-30 120])
xlim([-1 1])
ylabel('Change in SS firing (spks/s)')
xlabel('Lick dmax (ms)')
title('SS firing | G licks')

subplot(4,4,8)
hold on
plot(inds_span_corr, mean_change_SS_lick_rhythmic_r(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_rhythmic_r + se_change_SS_lick_rhythmic_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_rhythmic_r - se_change_SS_lick_rhythmic_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_change_SS_lick_weakrhythmic_r(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_weakrhythmic_r + se_change_SS_lick_weakrhythmic_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_weakrhythmic_r - se_change_SS_lick_weakrhythmic_r, 'r', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-30 120])
xlim([-1 1])
ylabel('Change in SS firing (spks/s)')
xlabel('Lick dmax (ms)')
title('SS firing | R licks')

subplot(4,4,9)
hold on
plot(inds_span_corr, mean_d_tip_lead_all(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_lead_all + se_d_tip_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_lead_all - se_d_tip_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_lag_all(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_lag_all + se_d_tip_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_lag_all - se_d_tip_lag_all, 'b', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_phase_all(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_phase_all + se_d_tip_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_phase_all - se_d_tip_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_antiphase_all(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_antiphase_all + se_d_tip_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_antiphase_all - se_d_tip_antiphase_all, 'm', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-1 1])
ylabel('Displacement (mm)')
xlabel('Lick dmax (ms)')
title('Lick trace | All licks')

subplot(4,4,10)
hold on
plot(inds_span_corr, mean_d_tip_lead_l(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_lead_l + se_d_tip_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_lead_l - se_d_tip_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_lag_l(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_lag_l + se_d_tip_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_lag_l - se_d_tip_lag_l, 'b', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_phase_l(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_phase_l + se_d_tip_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_phase_l - se_d_tip_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_antiphase_l(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_antiphase_l + se_d_tip_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_antiphase_l - se_d_tip_antiphase_l, 'm', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-1 1])
ylabel('Displacement (mm)')
xlabel('Lick dmax (ms)')
title('Lick trace | L licks')

subplot(4,4,11)
hold on
plot(inds_span_corr, mean_d_tip_lead_g(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_lead_g + se_d_tip_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_lead_g - se_d_tip_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_lag_g(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_lag_g + se_d_tip_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_lag_g - se_d_tip_lag_g, 'b', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_phase_g(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_phase_g + se_d_tip_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_phase_g - se_d_tip_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_antiphase_g(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_antiphase_g + se_d_tip_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_antiphase_g - se_d_tip_antiphase_g, 'm', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-1 1])
ylabel('Displacement (mm)')
xlabel('Lick dmax (ms)')
title('Lick trace | G licks')

subplot(4,4,12)
hold on
plot(inds_span_corr, mean_d_tip_lead_r(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_lead_r + se_d_tip_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_lead_r - se_d_tip_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_lag_r(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_lag_r + se_d_tip_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_lag_r - se_d_tip_lag_r, 'b', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_phase_r(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_phase_r + se_d_tip_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_phase_r - se_d_tip_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_antiphase_r(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_antiphase_r + se_d_tip_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_antiphase_r - se_d_tip_antiphase_r, 'm', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-1 1])
ylabel('Displacement (mm)')
xlabel('Lick dmax (ms)')
title('Lick trace | R licks')

subplot(4,4,13)
hold on
plot(inds_span_corr, mean_change_SS_lick_lead_all(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_lead_all + se_change_SS_lick_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_lead_all - se_change_SS_lick_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_change_SS_lick_lag_all(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_lag_all + se_change_SS_lick_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_lag_all - se_change_SS_lick_lag_all, 'b', 'LineWidth', 1)
plot(inds_span_corr, mean_change_SS_lick_phase_all(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_phase_all + se_change_SS_lick_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_phase_all - se_change_SS_lick_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_change_SS_lick_antiphase_all(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_antiphase_all + se_change_SS_lick_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_antiphase_all - se_change_SS_lick_antiphase_all, 'm', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-10 50])
xlim([-1 1])
ylabel('Change in SS firing (spks/s)')
xlabel('Lick dmax (ms)')
title('SS firing | All licks')

subplot(4,4,14)
hold on
plot(inds_span_corr, mean_change_SS_lick_lead_l(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_lead_l + se_change_SS_lick_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_lead_l - se_change_SS_lick_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_change_SS_lick_lag_l(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_lag_l + se_change_SS_lick_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_lag_l - se_change_SS_lick_lag_l, 'b', 'LineWidth', 1)
plot(inds_span_corr, mean_change_SS_lick_phase_l(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_phase_l + se_change_SS_lick_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_phase_l - se_change_SS_lick_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_change_SS_lick_antiphase_l(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_antiphase_l + se_change_SS_lick_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_antiphase_l - se_change_SS_lick_antiphase_l, 'm', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-10 50])
xlim([-1 1])
ylabel('Change in SS firing (spks/s)')
xlabel('Lick dmax (ms)')
title('SS firing | L licks')

subplot(4,4,15)
hold on
plot(inds_span_corr, mean_change_SS_lick_lead_g(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_lead_g + se_change_SS_lick_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_lead_g - se_change_SS_lick_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_change_SS_lick_lag_g(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_lag_g + se_change_SS_lick_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_lag_g - se_change_SS_lick_lag_g, 'b', 'LineWidth', 1)
plot(inds_span_corr, mean_change_SS_lick_phase_g(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_phase_g + se_change_SS_lick_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_phase_g - se_change_SS_lick_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_change_SS_lick_antiphase_g(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_antiphase_g + se_change_SS_lick_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_antiphase_g - se_change_SS_lick_antiphase_g, 'm', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-10 50])
xlim([-1 1])
ylabel('Change in SS firing (spks/s)')
xlabel('Lick dmax (ms)')
title('SS firing | G licks')

subplot(4,4,16)
hold on
plot(inds_span_corr, mean_change_SS_lick_lead_r(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_lead_r + se_change_SS_lick_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_lead_r - se_change_SS_lick_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_change_SS_lick_lag_r(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_lag_r + se_change_SS_lick_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_lag_r - se_change_SS_lick_lag_r, 'b', 'LineWidth', 1)
plot(inds_span_corr, mean_change_SS_lick_phase_r(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_phase_r + se_change_SS_lick_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_phase_r - se_change_SS_lick_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_change_SS_lick_antiphase_r(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_change_SS_lick_antiphase_r + se_change_SS_lick_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_change_SS_lick_antiphase_r - se_change_SS_lick_antiphase_r, 'm', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-10 50])
xlim([-1 1])
ylabel('Change in SS firing (spks/s)')
xlabel('Lick dmax (ms)')
title('SS firing | R licks')
%% Labeled power
%%% Power traces %%%

figure
subplot(1,4,1)
hold on;
plot(freq,mean_power_d_tip_all, 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
plot(freq,mean_power_d_tip_all + se_power_d_tip_all , 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
plot(freq,mean_power_d_tip_all - se_power_d_tip_all , 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
plot(freq,nanmean(power_SS_all(label_rhythmic_all,:)), 'Color', [1 0 0], 'LineWidth', 2);
plot(freq,nanmean(power_SS_all(label_rhythmic_all,:)) + nanstd(power_SS_all(label_rhythmic_all,:))/sqrt(number_pCell), 'Color', [1 0 0], 'LineWidth', 1);
plot(freq,nanmean(power_SS_all(label_rhythmic_all,:)) - nanstd(power_SS_all(label_rhythmic_all,:))/sqrt(number_pCell), 'Color', [1 0 0], 'LineWidth', 1);
plot(freq,nanmean(power_SS_all(label_weakrhythmic_all,:)), 'Color', [0 0 0], 'LineWidth', 2);
plot(freq,nanmean(power_SS_all(label_weakrhythmic_all,:)) + nanstd(power_SS_all(label_weakrhythmic_all,:))/sqrt(number_pCell), 'Color', [0 0 0], 'LineWidth', 1);
plot(freq,nanmean(power_SS_all(label_weakrhythmic_all,:)) - nanstd(power_SS_all(label_weakrhythmic_all,:))/sqrt(number_pCell), 'Color', [0 0 0], 'LineWidth', 1);
%     plot(freq,power_d_tip_all(label_rhythmic_all,:), 'Color', [0.7 0.7 0.7 0.1]);
%     plot(freq,power_SS_all(label_weakrhythmic_all,:), 'Color', [0 0 0 0.1]);
%     plot(freq,power_SS_all(label_rhythmic_all,:), 'Color', [1 0 0 0.1]);
xlabel('Freq (Hz)')
ylabel('Pow/Freq (dB/Hz)')
title('Power spectrum | All licks')

subplot(1,4,2)
hold on;
plot(freq,mean_power_d_tip_l, 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
plot(freq,mean_power_d_tip_l + se_power_d_tip_l , 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
plot(freq,mean_power_d_tip_l - se_power_d_tip_l , 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
plot(freq,nanmean(power_SS_l(label_rhythmic_l,:)), 'Color', [1 0 0], 'LineWidth', 2);
plot(freq,nanmean(power_SS_l(label_rhythmic_l,:)) + nanstd(power_SS_all(label_rhythmic_all,:))/sqrt(number_pCell), 'Color', [1 0 0], 'LineWidth', 1);
plot(freq,nanmean(power_SS_l(label_rhythmic_l,:)) - nanstd(power_SS_all(label_rhythmic_all,:))/sqrt(number_pCell), 'Color', [1 0 0], 'LineWidth', 1);
plot(freq,nanmean(power_SS_l(label_weakrhythmic_l,:)), 'Color', [0 0 0], 'LineWidth', 2);
plot(freq,nanmean(power_SS_l(label_weakrhythmic_l,:)) + nanstd(power_SS_l(label_weakrhythmic_l,:))/sqrt(number_pCell), 'Color', [0 0 0], 'LineWidth', 1);
plot(freq,nanmean(power_SS_l(label_weakrhythmic_l,:)) - nanstd(power_SS_l(label_weakrhythmic_l,:))/sqrt(number_pCell), 'Color', [0 0 0], 'LineWidth', 1);
%     plot(freq,power_d_tip_l(label_rhythmic_l,:), 'Color', [0.7 0.7 0.7 0.1]);
%     plot(freq,power_SS_l(label_weakrhythmic_l,:), 'Color', [0 0 0 0.1]);
%     plot(freq,power_SS_l(label_rhythmic_l,:), 'Color', [1 0 0 0.1]);
xlabel('Freq (Hz)')
ylabel('Pow/Freq (dB/Hz)')
title('Power spectrum | L licks')


subplot(1,4,3)
hold on;
plot(freq,mean_power_d_tip_g, 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
plot(freq,mean_power_d_tip_g + se_power_d_tip_g , 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
plot(freq,mean_power_d_tip_g - se_power_d_tip_g , 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
plot(freq,nanmean(power_SS_g(label_rhythmic_g,:)), 'Color', [1 0 0], 'LineWidth', 2);
plot(freq,nanmean(power_SS_g(label_rhythmic_g,:)) + nanstd(power_SS_all(label_rhythmic_all,:))/sqrt(number_pCell), 'Color', [1 0 0], 'LineWidth', 1);
plot(freq,nanmean(power_SS_g(label_rhythmic_g,:)) - nanstd(power_SS_all(label_rhythmic_all,:))/sqrt(number_pCell), 'Color', [1 0 0], 'LineWidth', 1);
plot(freq,nanmean(power_SS_g(label_weakrhythmic_g,:)), 'Color', [0 0 0], 'LineWidth', 2);
plot(freq,nanmean(power_SS_g(label_weakrhythmic_g,:)) + nanstd(power_SS_g(label_weakrhythmic_g,:))/sqrt(number_pCell), 'Color', [0 0 0], 'LineWidth', 1);
plot(freq,nanmean(power_SS_g(label_weakrhythmic_g,:)) - nanstd(power_SS_g(label_weakrhythmic_g,:))/sqrt(number_pCell), 'Color', [0 0 0], 'LineWidth', 1);
%     plot(freq,power_d_tip_g(label_rhythmic_g,:), 'Color', [0.7 0.7 0.7 0.1]);
%     plot(freq,power_SS_g(label_weakrhythmic_g,:), 'Color', [0 0 0 0.1]);
%     plot(freq,power_SS_g(label_rhythmic_g,:), 'Color', [1 0 0 0.1]);
xlabel('Freq (Hz)')
ylabel('Pow/Freq (dB/Hz)')
title('Power spectrum | G licks')


subplot(1,4,4)
hold on;
plot(freq,mean_power_d_tip_r, 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
plot(freq,mean_power_d_tip_r + se_power_d_tip_r , 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
plot(freq,mean_power_d_tip_r - se_power_d_tip_r , 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
plot(freq,nanmean(power_SS_r(label_rhythmic_r,:)), 'Color', [1 0 0], 'LineWidth', 2);
plot(freq,nanmean(power_SS_r(label_rhythmic_r,:)) + nanstd(power_SS_all(label_rhythmic_all,:))/sqrt(number_pCell), 'Color', [1 0 0], 'LineWidth', 1);
plot(freq,nanmean(power_SS_r(label_rhythmic_r,:)) - nanstd(power_SS_all(label_rhythmic_all,:))/sqrt(number_pCell), 'Color', [1 0 0], 'LineWidth', 1);
plot(freq,nanmean(power_SS_r(label_weakrhythmic_r,:)), 'Color', [0 0 0], 'LineWidth', 2);
plot(freq,nanmean(power_SS_r(label_weakrhythmic_r,:)) + nanstd(power_SS_r(label_weakrhythmic_r,:))/sqrt(number_pCell), 'Color', [0 0 0], 'LineWidth', 1);
plot(freq,nanmean(power_SS_r(label_weakrhythmic_r,:)) - nanstd(power_SS_r(label_weakrhythmic_r,:))/sqrt(number_pCell), 'Color', [0 0 0], 'LineWidth', 1);
%     plot(freq,power_d_tip_r(label_rhythmic_r,:), 'Color', [0.7 0.7 0.7 0.1]);
%     plot(freq,power_SS_r(label_rhythmic_r,:), 'Color', [1 0 0 0.1]);
%     plot(freq,power_SS_r(label_weakrhythmic_r,:), 'Color', [0 0 0 0.1]);
xlabel('Freq (Hz)')
ylabel('Pow/Freq (dB/Hz)')
title('Power spectrum | R licks')
%% Rhythmicity label histogram
figure
edges_all = 0 : 0.5 : 10;
subplot(2,4,1)
hold on
histogram(freq(ind_peak_power_SS_all + 1), edges_all, 'FaceColor', 'k')
histogram(freq(ind_peak_power_d_tip_all + 1), edges_all, 'FaceColor', [0.7 0.7 0.7])
xline(nanmean(freq(ind_peak_power_SS_all + 1)), 'k', 'LineWidth', 2);
xline(nanmean(freq(ind_peak_power_d_tip_all + 1)), 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([0 6])
ylim([0 50])
ylabel('Count')
xlabel('Frequency (Hz)')
title('All licks')

subplot(2,4,2)
hold on
histogram(freq(ind_peak_power_SS_l + 1), edges_all, 'FaceColor', 'b')
histogram(freq(ind_peak_power_d_tip_l + 1), edges_all, 'FaceColor', [0.7 0.7 0.7])
xline(nanmean(freq(ind_peak_power_SS_l + 1)), 'b', 'LineWidth', 2);
xline(nanmean(freq(ind_peak_power_d_tip_l + 1)), 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([0 6])
ylim([0 50])
ylabel('Count')
xlabel('Frequency (Hz)')
title('L licks')

subplot(2,4,3)
hold on
histogram(freq(ind_peak_power_SS_g + 1), edges_all, 'FaceColor', 'g')
histogram(freq(ind_peak_power_d_tip_g + 1), edges_all, 'FaceColor', [0.7 0.7 0.7])
xline(nanmean(freq(ind_peak_power_SS_g + 1)), 'g', 'LineWidth', 2);
xline(nanmean(freq(ind_peak_power_d_tip_g + 1)), 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([0 6])
ylim([0 50])
ylabel('Count')
xlabel('Frequency (Hz)')
title('G licks')

subplot(2,4,4)
hold on
histogram(freq(ind_peak_power_SS_r + 1), edges_all, 'FaceColor', 'r')
histogram(freq(ind_peak_power_d_tip_r + 1), edges_all, 'FaceColor', [0.7 0.7 0.7])
xline(nanmean(freq(ind_peak_power_SS_r + 1)), 'r', 'LineWidth', 2);
xline(nanmean(freq(ind_peak_power_d_tip_r + 1)), 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlim([0 6])
ylim([0 50])
ylabel('Count')
xlabel('Frequency (Hz)')
title('R licks')

edges_all = -6 : 0.5 : 6;

subplot(2,4,5)
hold on
histogram(freq(ind_peak_power_SS_all(label_rhythmic_all) + 1) - freq(ind_peak_power_d_tip_all(label_rhythmic_all) + 1), edges_all, 'FaceColor', [1 0 0])
histogram(freq(ind_peak_power_SS_all(label_weakrhythmic_all) + 1) - freq(ind_peak_power_d_tip_all(label_weakrhythmic_all) + 1), edges_all, 'FaceColor', [0 0 0])
xline(nanmean(freq(ind_peak_power_SS_all(label_rhythmic_all) + 1) - freq(ind_peak_power_d_tip_all(label_rhythmic_all) + 1)), 'r', 'LineWidth', 2);
xline(nanmean(freq(ind_peak_power_SS_all(label_weakrhythmic_all) + 1) - freq(ind_peak_power_d_tip_all(label_weakrhythmic_all) + 1)), 'k', 'LineWidth', 2);
xlim([-6 6])
ylim([0 50])
ylabel('Count')
xlabel('Diff in frequency (Hz)')
title('freq(SS) - freq(lick) | All licks ')

subplot(2,4,6)
hold on
histogram(freq(ind_peak_power_SS_l(label_rhythmic_l) + 1) - freq(ind_peak_power_d_tip_l(label_rhythmic_l) + 1), edges_all, 'FaceColor', [1 0 0])
histogram(freq(ind_peak_power_SS_l(label_weakrhythmic_l) + 1) - freq(ind_peak_power_d_tip_l(label_weakrhythmic_l) + 1), edges_all, 'FaceColor', [0 0 0])
xline(nanmean(freq(ind_peak_power_SS_l(label_rhythmic_l) + 1) - freq(ind_peak_power_d_tip_l(label_rhythmic_l) + 1)), 'r', 'LineWidth', 2);
xline(nanmean(freq(ind_peak_power_SS_l(label_weakrhythmic_l) + 1) - freq(ind_peak_power_d_tip_l(label_weakrhythmic_l) + 1)), 'k', 'LineWidth', 2);
xlim([-6 6])
ylim([0 50])
ylabel('Count')
xlabel('Diff in frequency (Hz)')
title('freq(SS) - freq(lick) | L licks ')

subplot(2,4,7)
hold on
histogram(freq(ind_peak_power_SS_g(label_rhythmic_g) + 1) - freq(ind_peak_power_d_tip_g(label_rhythmic_g) + 1), edges_all, 'FaceColor', [1 0 0])
histogram(freq(ind_peak_power_SS_g(label_weakrhythmic_g) + 1) - freq(ind_peak_power_d_tip_g(label_weakrhythmic_g) + 1), edges_all, 'FaceColor', [0 0 0])
xline(nanmean(freq(ind_peak_power_SS_g(label_rhythmic_g) + 1) - freq(ind_peak_power_d_tip_g(label_rhythmic_g) + 1)), 'r', 'LineWidth', 2);
xline(nanmean(freq(ind_peak_power_SS_g(label_weakrhythmic_g) + 1) - freq(ind_peak_power_d_tip_g(label_weakrhythmic_g) + 1)), 'k', 'LineWidth', 2);
xlim([-6 6])
ylim([0 50])
ylabel('Count')
xlabel('Diff in frequency (Hz)')
title('freq(SS) - freq(lick) | G licks ')

subplot(2,4,8)
hold on
histogram(freq(ind_peak_power_SS_r(label_rhythmic_r) + 1) - freq(ind_peak_power_d_tip_r(label_rhythmic_r) + 1), edges_all, 'FaceColor', [1 0 0])
histogram(freq(ind_peak_power_SS_r(label_weakrhythmic_r) + 1) - freq(ind_peak_power_d_tip_r(label_weakrhythmic_r) + 1), edges_all, 'FaceColor', [0 0 0])
xline(nanmean(freq(ind_peak_power_SS_r(label_rhythmic_r) + 1) - freq(ind_peak_power_d_tip_r(label_rhythmic_r) + 1)), 'r', 'LineWidth', 2);
xline(nanmean(freq(ind_peak_power_SS_r(label_weakrhythmic_r) + 1) - freq(ind_peak_power_d_tip_r(label_weakrhythmic_r) + 1)), 'k', 'LineWidth', 2);
xlim([-6 6])
ylim([0 50])
ylabel('Count')
xlabel('Diff in frequency (Hz)')
title('freq(SS) - freq(lick) | R licks ')
%% Phase label histogram
figure
edges_all = (0 : 15 : 180);

subplot(1,4,1)
hold on
histogram(abs(phase_diff_all(label_phase_all)), edges_all, 'FaceColor', 'k')
histogram(abs(phase_diff_all(label_antiphase_all)), edges_all, 'FaceColor', 'm')
histogram(abs(phase_diff_all(label_lead_all)), edges_all, 'FaceColor', 'r')
histogram(abs(phase_diff_all(label_lag_all)), edges_all, 'FaceColor', 'b')
xline(nanmean(abs(phase_diff_all(label_phase_all))), 'k', 'LineWidth', 2);
xline(nanmean(abs(phase_diff_all(label_antiphase_all))), 'm', 'LineWidth', 2);
xline(nanmean(abs(phase_diff_all(label_lead_all))), 'r', 'LineWidth', 2);
xline(nanmean(abs(phase_diff_all(label_lag_all))), 'b', 'LineWidth', 2);
xlim([-10 190])
ylabel('Count')
xlabel('Phase diff (deg)')
title('All licks')

subplot(1,4,2)
hold on
histogram(abs(phase_diff_l(label_phase_l)), edges_all, 'FaceColor', 'k')
histogram(abs(phase_diff_l(label_antiphase_l)), edges_all, 'FaceColor', 'm')
histogram(abs(phase_diff_l(label_lead_l)), edges_all, 'FaceColor', 'r')
histogram(abs(phase_diff_l(label_lag_l)), edges_all, 'FaceColor', 'b')
xline(nanmean(abs(phase_diff_l(label_phase_l))), 'k', 'LineWidth', 2);
xline(nanmean(abs(phase_diff_l(label_antiphase_l))), 'm', 'LineWidth', 2);
xline(nanmean(abs(phase_diff_l(label_lead_l))), 'r', 'LineWidth', 2);
xline(nanmean(abs(phase_diff_l(label_lag_l))), 'b', 'LineWidth', 2);
xlim([-10 190])
ylabel('Count')
xlabel('Phase diff (deg)')
title('L licks')

subplot(1,4,3)
hold on
histogram(abs(phase_diff_g(label_phase_g)), edges_all, 'FaceColor', 'k')
histogram(abs(phase_diff_g(label_antiphase_g)), edges_all, 'FaceColor', 'm')
histogram(abs(phase_diff_g(label_lead_g)), edges_all, 'FaceColor', 'r')
histogram(abs(phase_diff_g(label_lag_g)), edges_all, 'FaceColor', 'b')
xline(nanmean(abs(phase_diff_g(label_phase_g))), 'k', 'LineWidth', 2);
xline(nanmean(abs(phase_diff_g(label_antiphase_g))), 'm', 'LineWidth', 2);
xline(nanmean(abs(phase_diff_g(label_lead_g))), 'r', 'LineWidth', 2);
xline(nanmean(abs(phase_diff_g(label_lag_g))), 'b', 'LineWidth', 2);
xlim([-10 190])
ylabel('Count')
xlabel('Phase diff (deg)')
title('G licks')

subplot(1,4,4)
hold on
histogram(abs(phase_diff_r(label_phase_r)), edges_all, 'FaceColor', 'k')
histogram(abs(phase_diff_r(label_antiphase_r)), edges_all, 'FaceColor', 'm')
histogram(abs(phase_diff_r(label_lead_r)), edges_all, 'FaceColor', 'r')
histogram(abs(phase_diff_r(label_lag_r)), edges_all, 'FaceColor', 'b')
xline(nanmean(abs(phase_diff_r(label_phase_r))), 'k', 'LineWidth', 2);
xline(nanmean(abs(phase_diff_r(label_antiphase_r))), 'm', 'LineWidth', 2);
xline(nanmean(abs(phase_diff_r(label_lead_r))), 'r', 'LineWidth', 2);
xline(nanmean(abs(phase_diff_r(label_lag_r))), 'b', 'LineWidth', 2);
xlim([-10 190])
ylabel('Count')
xlabel('Phase diff (deg)')
title('R licks')

%%  ANALYSIS PLOTS BASED ON LABELS - ONSET
clearvars -except ALL_PCELL_COMPRESSED_DATA POPULATION

number_pCell =  length(ALL_PCELL_COMPRESSED_DATA);

for counter_pCell = 1 : number_pCell
    id(counter_pCell, 1) = string(ALL_PCELL_COMPRESSED_DATA(counter_pCell).id{1, 1});
    
    %%% d_tip %%%
    d_tip_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_all_corr);
    d_tip_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_grooming_corr);
    d_tip_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_r_corr);
    d_tip_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_l_corr);
    
    %%% v_tip %%%
    v_tip_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_v_tip_100_all_corr);
    v_tip_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_v_tip_100_grooming_corr);
    v_tip_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_v_tip_100_r_corr);
    v_tip_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_v_tip_100_l_corr);
    
    %%% SS_str_bout %%%
    SS_str_bout_all(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_str_bout * 100), 15, 'sgolay', 2);
    SS_str_bout_baseline_all(counter_pCell, :) = nanmean(SS_str_bout_all(counter_pCell, 1:50));
    SS_str_bout_r(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_str_bout_r * 100), 15, 'sgolay', 2);
    SS_str_bout_baseline_r(counter_pCell, :) = nanmean(SS_str_bout_r(counter_pCell, 1:50));
    SS_str_bout_l(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_str_bout_l * 100), 15, 'sgolay', 2);
    SS_str_bout_baseline_l(counter_pCell, :) = nanmean(SS_str_bout_l(counter_pCell, 1:50));
    
    %%% SS_lick %%%
    SS_lick_all(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_all_corr * 100), 15, 'sgolay', 2);
    change_SS_lick_all(counter_pCell,:) = SS_lick_all(counter_pCell,:) - SS_str_bout_baseline_all(counter_pCell, 1);
    mean_window_SS_lick_all(counter_pCell,1) = nanmean(SS_lick_all(counter_pCell, :));
    mean_change_window_SS_lick_all(counter_pCell,1) = nanmean(change_SS_lick_all(counter_pCell, :));
    SS_lick_all(counter_pCell,:) = SS_lick_all(counter_pCell,:)/mean_window_SS_lick_all(counter_pCell,1);
    
    SS_lick_r(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_r_corr * 100), 15, 'sgolay', 2);
    change_SS_lick_r(counter_pCell,:) = SS_lick_r(counter_pCell,:) - SS_str_bout_baseline_r(counter_pCell, 1);
    mean_window_SS_lick_r(counter_pCell,1) = nanmean(SS_lick_r(counter_pCell, :));
    mean_change_window_SS_lick_r(counter_pCell,1) = nanmean(change_SS_lick_r(counter_pCell, :));
    SS_lick_r(counter_pCell,:) = SS_lick_r(counter_pCell,:)/mean_window_SS_lick_r(counter_pCell,1);
    
    SS_lick_g(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_grooming_corr * 100), 15, 'sgolay', 2);
    change_SS_lick_g(counter_pCell,:) = SS_lick_g(counter_pCell,:) - SS_str_bout_baseline_all(counter_pCell, 1);
    mean_window_SS_lick_g(counter_pCell,1) = nanmean(SS_lick_g(counter_pCell, :));
    mean_change_window_SS_lick_g(counter_pCell,1) = nanmean(change_SS_lick_g(counter_pCell, :));
    SS_lick_g(counter_pCell,:) = SS_lick_g(counter_pCell,:)/mean_window_SS_lick_g(counter_pCell,1);
    
    SS_lick_l(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_l_corr * 100), 15, 'sgolay', 2);
    change_SS_lick_l(counter_pCell,:) = SS_lick_l(counter_pCell,:) - SS_str_bout_baseline_l(counter_pCell, 1);
    mean_window_SS_lick_l(counter_pCell,1) = nanmean(SS_lick_l(counter_pCell, :));
    mean_change_window_SS_lick_l(counter_pCell,1) = nanmean(change_SS_lick_l(counter_pCell, :));
    SS_lick_l(counter_pCell,:) = SS_lick_l(counter_pCell,:)/mean_window_SS_lick_l(counter_pCell,1);
    
end

%%% LOAD LABELS FROM POPULATION STRUCTURE %%%
label_large_increase_all = POPULATION.dmax.label_large_increase_all;
label_medium_increase_all  = POPULATION.dmax.label_medium_increase_all;
label_small_increase_all = POPULATION.dmax.label_small_increase_all;
label_phase_all = POPULATION.dmax.label_phase_all;
label_small_decrease_all = POPULATION.dmax.label_small_decrease_all;
label_medium_decrease_all = POPULATION.dmax.label_medium_decrease_all;
label_large_decrease_all = POPULATION.dmax.label_large_decrease_all;
label_large_r = POPULATION.dmax.label_large_r;
label_medium_r = POPULATION.dmax.label_medium_r;
label_small_r = POPULATION.dmax.label_small_r;
label_none_dir = POPULATION.dmax.label_none_dir;
label_small_l = POPULATION.dmax.label_small_l;
label_medium_l = POPULATION.dmax.label_medium_l;
label_large_l = POPULATION.dmax.label_large_l;

label_rhythmic_all = POPULATION.dmax.label_rhythmic_all;
label_weakrhythmic_all = POPULATION.dmax.label_weakrhythmic_all;
label_rhythmic_l = POPULATION.dmax.label_rhythmic_l;
label_weakrhythmic_l = POPULATION.dmax.label_weakrhythmic_l;
label_rhythmic_g = POPULATION.dmax.label_rhythmic_g;
label_weakrhythmic_g = POPULATION.dmax.label_weakrhythmic_g;
label_rhythmic_r = POPULATION.dmax.label_rhythmic_r;
label_weakrhythmic_r = POPULATION.dmax.label_weakrhythmic_r;

label_lead_all = POPULATION.dmax.label_lead_all;
label_lead30to90_all  = POPULATION.dmax.label_lead30to90_all;
label_lead90to150_all  = POPULATION.dmax.label_lead90to150_all;
label_lag_all = POPULATION.dmax.label_lag_all;
label_lag30to90_all = POPULATION.dmax.label_lag30to90_all;
label_lag90to150_all = POPULATION.dmax.label_lag90to150_all;
label_phase_all = POPULATION.dmax.label_phase_all;
label_antiphase_all = POPULATION.dmax.label_antiphase_all;

label_lead_l = POPULATION.dmax.label_lead_l;
label_lead30to90_l  = POPULATION.dmax.label_lead30to90_l;
label_lead90to150_l  = POPULATION.dmax.label_lead90to150_l;
label_lag_l = POPULATION.dmax.label_lag_l;
label_lag30to90_l = POPULATION.dmax.label_lag30to90_l;
label_lag90to150_l = POPULATION.dmax.label_lag90to150_l;
label_phase_l = POPULATION.dmax.label_phase_l;
label_antiphase_l = POPULATION.dmax.label_antiphase_l;

label_lead_g = POPULATION.dmax.label_lead_g;
label_lead30to90_g  = POPULATION.dmax.label_lead30to90_g;
label_lead90to150_g  = POPULATION.dmax.label_lead90to150_g;
label_lag_g = POPULATION.dmax.label_lag_g;
label_lag30to90_g = POPULATION.dmax.label_lag30to90_g;
label_lag90to150_g = POPULATION.dmax.label_lag90to150_g;
label_phase_g = POPULATION.dmax.label_phase_g;
label_antiphase_g = POPULATION.dmax.label_antiphase_g;

label_lead_r = POPULATION.dmax.label_lead_r;
label_lead30to90_r  = POPULATION.dmax.label_lead30to90_r;
label_lead90to150_r  = POPULATION.dmax.label_lead90to150_r;
label_lag_r = POPULATION.dmax.label_lag_r;
label_lag30to90_r = POPULATION.dmax.label_lag30to90_r;
label_lag90to150_r = POPULATION.dmax.label_lag90to150_r;
label_phase_r = POPULATION.dmax.label_phase_r;
label_antiphase_r = POPULATION.dmax.label_antiphase_r;

%%% Mean d_tip based on labels
mean_d_tip_all = nanmean(d_tip_all);
mean_d_tip_l = nanmean(d_tip_l);
mean_d_tip_g = nanmean(d_tip_g);
mean_d_tip_r = nanmean(d_tip_r);

mean_d_tip_rhythmic_phase_all = nanmean(d_tip_all(label_rhythmic_all & label_phase_all,:));
se_d_tip_rhythmic_phase_all = nanstd(d_tip_all(label_rhythmic_all & label_phase_all,:))/sqrt(sum(label_rhythmic_all & label_phase_all));

mean_d_tip_rhythmic_antiphase_all = nanmean(d_tip_all(label_rhythmic_all & label_antiphase_all,:));
se_d_tip_rhythmic_antiphase_all = nanstd(d_tip_all(label_rhythmic_all & label_antiphase_all,:))/sqrt(sum(label_rhythmic_all & label_antiphase_all));

mean_d_tip_rhythmic_lead_all = nanmean(d_tip_all(label_rhythmic_all & label_lead_all,:));
se_d_tip_rhythmic_lead_all = nanstd(d_tip_all(label_rhythmic_all & label_lead_all,:))/sqrt(sum(label_rhythmic_all & label_lead_all));
mean_d_tip_rhythmic_lead30to90_all = nanmean(d_tip_all(label_rhythmic_all & label_lead30to90_all,:));
se_d_tip_rhythmic_lead30to90_all = nanstd(d_tip_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_d_tip_rhythmic_lead90to150_all = nanmean(d_tip_all(label_rhythmic_all & label_lead90to150_all,:));
se_d_tip_rhythmic_lead90to150_all = nanstd(d_tip_all(label_rhythmic_all & label_lead90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lead90to150_all));

mean_d_tip_rhythmic_lag_all = nanmean(d_tip_all(label_rhythmic_all & label_lag_all,:));
se_d_tip_rhythmic_lag_all = nanstd(d_tip_all(label_rhythmic_all & label_lag_all,:))/sqrt(sum(label_rhythmic_all & label_lag_all));
mean_d_tip_rhythmic_lead30to90_all = nanmean(d_tip_all(label_rhythmic_all & label_lead30to90_all,:));
se_d_tip_rhythmic_lead30to90_all = nanstd(d_tip_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_d_tip_rhythmic_lag90to150_all = nanmean(d_tip_all(label_rhythmic_all & label_lag90to150_all,:));
se_d_tip_rhythmic_lag90to150_all = nanstd(d_tip_all(label_rhythmic_all & label_lag90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lag90to150_all));

mean_d_tip_rhythmic_phase_l = nanmean(d_tip_l(label_rhythmic_l & label_phase_l,:));
se_d_tip_rhythmic_phase_l = nanstd(d_tip_l(label_rhythmic_l & label_phase_l,:))/sqrt(sum(label_rhythmic_l & label_phase_l));

mean_d_tip_rhythmic_antiphase_l = nanmean(d_tip_l(label_rhythmic_l & label_antiphase_l,:));
se_d_tip_rhythmic_antiphase_l = nanstd(d_tip_l(label_rhythmic_l & label_antiphase_l,:))/sqrt(sum(label_rhythmic_l & label_antiphase_l));

mean_d_tip_rhythmic_lead_l = nanmean(d_tip_l(label_rhythmic_l & label_lead_l,:));
se_d_tip_rhythmic_lead_l = nanstd(d_tip_l(label_rhythmic_l & label_lead_l,:))/sqrt(sum(label_rhythmic_l & label_lead_l));
mean_d_tip_rhythmic_lead30to90_l = nanmean(d_tip_l(label_rhythmic_l & label_lead30to90_l,:));
se_d_tip_rhythmic_lead30to90_l = nanstd(d_tip_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_d_tip_rhythmic_lead90to150_l = nanmean(d_tip_l(label_rhythmic_l & label_lead90to150_l,:));
se_d_tip_rhythmic_lead90to150_l = nanstd(d_tip_l(label_rhythmic_l & label_lead90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lead90to150_l));

mean_d_tip_rhythmic_lag_l = nanmean(d_tip_l(label_rhythmic_l & label_lag_l,:));
se_d_tip_rhythmic_lag_l = nanstd(d_tip_l(label_rhythmic_l & label_lag_l,:))/sqrt(sum(label_rhythmic_l & label_lag_l));
mean_d_tip_rhythmic_lead30to90_l = nanmean(d_tip_l(label_rhythmic_l & label_lead30to90_l,:));
se_d_tip_rhythmic_lead30to90_l = nanstd(d_tip_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_d_tip_rhythmic_lag90to150_l = nanmean(d_tip_l(label_rhythmic_l & label_lag90to150_l,:));
se_d_tip_rhythmic_lag90to150_l = nanstd(d_tip_l(label_rhythmic_l & label_lag90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lag90to150_l));

mean_d_tip_rhythmic_phase_g = nanmean(d_tip_g(label_rhythmic_g & label_phase_g,:));
se_d_tip_rhythmic_phase_g = nanstd(d_tip_g(label_rhythmic_g & label_phase_g,:))/sqrt(sum(label_rhythmic_g & label_phase_g));

mean_d_tip_rhythmic_antiphase_g = nanmean(d_tip_g(label_rhythmic_g & label_antiphase_g,:));
se_d_tip_rhythmic_antiphase_g = nanstd(d_tip_g(label_rhythmic_g & label_antiphase_g,:))/sqrt(sum(label_rhythmic_g & label_antiphase_g));

mean_d_tip_rhythmic_lead_g = nanmean(d_tip_g(label_rhythmic_g & label_lead_g,:));
se_d_tip_rhythmic_lead_g = nanstd(d_tip_g(label_rhythmic_g & label_lead_g,:))/sqrt(sum(label_rhythmic_g & label_lead_g));
mean_d_tip_rhythmic_lead30to90_g = nanmean(d_tip_g(label_rhythmic_g & label_lead30to90_g,:));
se_d_tip_rhythmic_lead30to90_g = nanstd(d_tip_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_d_tip_rhythmic_lead90to150_g = nanmean(d_tip_g(label_rhythmic_g & label_lead90to150_g,:));
se_d_tip_rhythmic_lead90to150_g = nanstd(d_tip_g(label_rhythmic_g & label_lead90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lead90to150_g));

mean_d_tip_rhythmic_lag_g = nanmean(d_tip_g(label_rhythmic_g & label_lag_g,:));
se_d_tip_rhythmic_lag_g = nanstd(d_tip_g(label_rhythmic_g & label_lag_g,:))/sqrt(sum(label_rhythmic_g & label_lag_g));
mean_d_tip_rhythmic_lead30to90_g = nanmean(d_tip_g(label_rhythmic_g & label_lead30to90_g,:));
se_d_tip_rhythmic_lead30to90_g = nanstd(d_tip_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_d_tip_rhythmic_lag90to150_g = nanmean(d_tip_g(label_rhythmic_g & label_lag90to150_g,:));
se_d_tip_rhythmic_lag90to150_g = nanstd(d_tip_g(label_rhythmic_g & label_lag90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lag90to150_g));

mean_d_tip_rhythmic_phase_r = nanmean(d_tip_r(label_rhythmic_r & label_phase_r,:));
se_d_tip_rhythmic_phase_r = nanstd(d_tip_r(label_rhythmic_r & label_phase_r,:))/sqrt(sum(label_rhythmic_r & label_phase_r));

mean_d_tip_rhythmic_antiphase_r = nanmean(d_tip_g(label_rhythmic_r & label_antiphase_r,:));
se_d_tip_rhythmic_antiphase_r = nanstd(d_tip_g(label_rhythmic_r & label_antiphase_r,:))/sqrt(sum(label_rhythmic_r & label_antiphase_r));

mean_d_tip_rhythmic_lead_r = nanmean(d_tip_r(label_rhythmic_r & label_lead_r,:));
se_d_tip_rhythmic_lead_r = nanstd(d_tip_r(label_rhythmic_r & label_lead_r,:))/sqrt(sum(label_rhythmic_r & label_lead_r));
mean_d_tip_rhythmic_lead30to90_r = nanmean(d_tip_r(label_rhythmic_r & label_lead30to90_r,:));
se_d_tip_rhythmic_lead30to90_r = nanstd(d_tip_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_d_tip_rhythmic_lead90to150_r = nanmean(d_tip_r(label_rhythmic_r & label_lead90to150_r,:));
se_d_tip_rhythmic_lead90to150_r = nanstd(d_tip_r(label_rhythmic_r & label_lead90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lead90to150_r));

mean_d_tip_rhythmic_lag_r = nanmean(d_tip_r(label_rhythmic_r & label_lag_r,:));
se_d_tip_rhythmic_lag_r = nanstd(d_tip_r(label_rhythmic_r & label_lag_r,:))/sqrt(sum(label_rhythmic_r & label_lag_r));
mean_d_tip_rhythmic_lead30to90_r = nanmean(d_tip_r(label_rhythmic_r & label_lead30to90_r,:));
se_d_tip_rhythmic_lead30to90_r = nanstd(d_tip_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_d_tip_rhythmic_lag90to150_r = nanmean(d_tip_r(label_rhythmic_r & label_lag90to150_r,:));
se_d_tip_rhythmic_lag90to150_r = nanstd(d_tip_r(label_rhythmic_r & label_lag90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lag90to150_r));

%%% Mean v_tip based on labels
mean_v_tip_all = nanmean(v_tip_all);
mean_v_tip_l = nanmean(v_tip_l);
mean_v_tip_g = nanmean(v_tip_g);
mean_v_tip_r = nanmean(v_tip_r);

mean_v_tip_rhythmic_phase_all = nanmean(v_tip_all(label_rhythmic_all & label_phase_all,:));
se_v_tip_rhythmic_phase_all = nanstd(v_tip_all(label_rhythmic_all & label_phase_all,:))/sqrt(sum(label_rhythmic_all & label_phase_all));

mean_v_tip_rhythmic_antiphase_all = nanmean(v_tip_all(label_rhythmic_all & label_antiphase_all,:));
se_v_tip_rhythmic_antiphase_all = nanstd(v_tip_all(label_rhythmic_all & label_antiphase_all,:))/sqrt(sum(label_rhythmic_all & label_antiphase_all));

mean_v_tip_rhythmic_lead_all = nanmean(v_tip_all(label_rhythmic_all & label_lead_all,:));
se_v_tip_rhythmic_lead_all = nanstd(v_tip_all(label_rhythmic_all & label_lead_all,:))/sqrt(sum(label_rhythmic_all & label_lead_all));
mean_v_tip_rhythmic_lead30to90_all = nanmean(v_tip_all(label_rhythmic_all & label_lead30to90_all,:));
se_v_tip_rhythmic_lead30to90_all = nanstd(v_tip_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_v_tip_rhythmic_lead90to150_all = nanmean(v_tip_all(label_rhythmic_all & label_lead90to150_all,:));
se_v_tip_rhythmic_lead90to150_all = nanstd(v_tip_all(label_rhythmic_all & label_lead90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lead90to150_all));

mean_v_tip_rhythmic_lag_all = nanmean(v_tip_all(label_rhythmic_all & label_lag_all,:));
se_v_tip_rhythmic_lag_all = nanstd(v_tip_all(label_rhythmic_all & label_lag_all,:))/sqrt(sum(label_rhythmic_all & label_lag_all));
mean_v_tip_rhythmic_lead30to90_all = nanmean(v_tip_all(label_rhythmic_all & label_lead30to90_all,:));
se_v_tip_rhythmic_lead30to90_all = nanstd(v_tip_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_v_tip_rhythmic_lag90to150_all = nanmean(v_tip_all(label_rhythmic_all & label_lag90to150_all,:));
se_v_tip_rhythmic_lag90to150_all = nanstd(v_tip_all(label_rhythmic_all & label_lag90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lag90to150_all));

mean_v_tip_rhythmic_phase_l = nanmean(v_tip_l(label_rhythmic_l & label_phase_l,:));
se_v_tip_rhythmic_phasell = nanstd(v_tip_l(label_rhythmic_l & label_phase_l,:))/sqrt(sum(label_rhythmic_l & label_phase_l));

mean_v_tip_rhythmic_antiphase_l = nanmean(v_tip_l(label_rhythmic_l & label_antiphase_l,:));
se_v_tip_rhythmic_antiphase_l = nanstd(v_tip_l(label_rhythmic_l & label_antiphase_l,:))/sqrt(sum(label_rhythmic_l & label_antiphase_l));

mean_v_tip_rhythmic_lead_l = nanmean(v_tip_l(label_rhythmic_l & label_lead_l,:));
se_v_tip_rhythmic_lead_l = nanstd(v_tip_l(label_rhythmic_l & label_lead_l,:))/sqrt(sum(label_rhythmic_l & label_lead_l));
mean_v_tip_rhythmic_lead30to90_l = nanmean(v_tip_l(label_rhythmic_l & label_lead30to90_l,:));
se_v_tip_rhythmic_lead30to90_l = nanstd(v_tip_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_v_tip_rhythmic_lead90to150_l = nanmean(v_tip_l(label_rhythmic_l & label_lead90to150_l,:));
se_v_tip_rhythmic_lead90to150_l = nanstd(v_tip_l(label_rhythmic_l & label_lead90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lead90to150_l));

mean_v_tip_rhythmic_lag_l = nanmean(v_tip_l(label_rhythmic_l & label_lag_l,:));
se_v_tip_rhythmic_lag_l = nanstd(v_tip_l(label_rhythmic_l & label_lag_l,:))/sqrt(sum(label_rhythmic_l & label_lag_l));
mean_v_tip_rhythmic_lead30to90_l = nanmean(v_tip_l(label_rhythmic_l & label_lead30to90_l,:));
se_v_tip_rhythmic_lead30to90_l = nanstd(v_tip_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_v_tip_rhythmic_lag90to150_l = nanmean(v_tip_l(label_rhythmic_l & label_lag90to150_l,:));
se_v_tip_rhythmic_lag90to150_l = nanstd(v_tip_l(label_rhythmic_l & label_lag90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lag90to150_l));

mean_v_tip_rhythmic_phase_g = nanmean(v_tip_g(label_rhythmic_g & label_phase_g,:));
se_v_tip_rhythmic_phase_g = nanstd(v_tip_g(label_rhythmic_g & label_phase_g,:))/sqrt(sum(label_rhythmic_g & label_phase_g));

mean_v_tip_rhythmic_antiphase_g = nanmean(v_tip_g(label_rhythmic_g & label_antiphase_g,:));
se_v_tip_rhythmic_antiphase_g = nanstd(v_tip_g(label_rhythmic_g & label_antiphase_g,:))/sqrt(sum(label_rhythmic_g & label_antiphase_g));

mean_v_tip_rhythmic_lead_g = nanmean(v_tip_g(label_rhythmic_g & label_lead_g,:));
se_v_tip_rhythmic_lead_g = nanstd(v_tip_g(label_rhythmic_g & label_lead_g,:))/sqrt(sum(label_rhythmic_g & label_lead_g));
mean_v_tip_rhythmic_lead30to90_g = nanmean(v_tip_g(label_rhythmic_g & label_lead30to90_g,:));
se_v_tip_rhythmic_lead30to90_g = nanstd(v_tip_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_v_tip_rhythmic_lead90to150_g = nanmean(v_tip_g(label_rhythmic_g & label_lead90to150_g,:));
se_v_tip_rhythmic_lead90to150_g = nanstd(v_tip_g(label_rhythmic_g & label_lead90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lead90to150_g));

mean_v_tip_rhythmic_lag_g = nanmean(v_tip_g(label_rhythmic_g & label_lag_g,:));
se_v_tip_rhythmic_lag_g = nanstd(v_tip_g(label_rhythmic_g & label_lag_g,:))/sqrt(sum(label_rhythmic_g & label_lag_g));
mean_v_tip_rhythmic_lead30to90_g = nanmean(v_tip_g(label_rhythmic_g & label_lead30to90_g,:));
se_v_tip_rhythmic_lead30to90_g = nanstd(v_tip_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_v_tip_rhythmic_lag90to150_g = nanmean(v_tip_g(label_rhythmic_g & label_lag90to150_g,:));
se_v_tip_rhythmic_lag90to150_g = nanstd(v_tip_g(label_rhythmic_g & label_lag90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lag90to150_g));

mean_v_tip_rhythmic_phase_r = nanmean(v_tip_r(label_rhythmic_r & label_phase_r,:));
se_v_tip_rhythmic_phase_r = nanstd(v_tip_r(label_rhythmic_r & label_phase_r,:))/sqrt(sum(label_rhythmic_r & label_phase_r));

mean_v_tip_rhythmic_antiphase_r = nanmean(v_tip_r(label_rhythmic_r & label_antiphase_r,:));
se_v_tip_rhythmic_antiphaser = nanstd(v_tip_r(label_rhythmic_r & label_antiphase_r,:))/sqrt(sum(label_rhythmic_r & label_antiphase_r));

mean_v_tip_rhythmic_lead_r = nanmean(v_tip_r(label_rhythmic_r & label_lead_r,:));
se_v_tip_rhythmic_lead_r = nanstd(v_tip_r(label_rhythmic_r & label_lead_r,:))/sqrt(sum(label_rhythmic_r & label_lead_r));
mean_v_tip_rhythmic_lead30to90_r = nanmean(v_tip_r(label_rhythmic_r & label_lead30to90_r,:));
se_v_tip_rhythmic_lead30to90_r = nanstd(v_tip_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_v_tip_rhythmic_lead90to150_r = nanmean(v_tip_r(label_rhythmic_r & label_lead90to150_r,:));
se_v_tip_rhythmic_lead90to150_r = nanstd(v_tip_r(label_rhythmic_r & label_lead90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lead90to150_r));

mean_v_tip_rhythmic_lag_r = nanmean(v_tip_r(label_rhythmic_r & label_lag_r,:));
se_v_tip_rhythmic_lag_r = nanstd(v_tip_r(label_rhythmic_r & label_lag_r,:))/sqrt(sum(label_rhythmic_r & label_lag_r));
mean_v_tip_rhythmic_lead30to90_r = nanmean(v_tip_r(label_rhythmic_r & label_lead30to90_r,:));
se_v_tip_rhythmic_lead30to90_r = nanstd(v_tip_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_v_tip_rhythmic_lag90to150_r = nanmean(v_tip_r(label_rhythmic_r & label_lag90to150_r,:));
se_v_tip_rhythmic_lag90to150_r = nanstd(v_tip_r(label_rhythmic_r & label_lag90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lag90to150_r));

%%% Mean change_SS_lick based on labels
mean_SS_lick_rhythmic_phase_all = nanmean(SS_lick_all(label_rhythmic_all & label_phase_all,:));
se_SS_lick_rhythmic_phase_all = nanstd(SS_lick_all(label_rhythmic_all & label_phase_all,:))/sqrt(sum(label_rhythmic_all & label_phase_all));

mean_SS_lick_rhythmic_antiphase_all = nanmean(SS_lick_all(label_rhythmic_all & label_antiphase_all,:));
se_SS_lick_rhythmic_antiphase_all = nanstd(SS_lick_all(label_rhythmic_all & label_antiphase_all,:))/sqrt(sum(label_rhythmic_all & label_antiphase_all));

mean_SS_lick_rhythmic_lead_all = nanmean(SS_lick_all(label_rhythmic_all & label_lead_all,:));
se_SS_lick_rhythmic_lead_all = nanstd(SS_lick_all(label_rhythmic_all & label_lead_all,:))/sqrt(sum(label_rhythmic_all & label_lead_all));
mean_SS_lick_rhythmic_lead30to90_all = nanmean(SS_lick_all(label_rhythmic_all & label_lead30to90_all,:));
se_SS_lick_rhythmic_lead30to90_all = nanstd(SS_lick_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_SS_lick_rhythmic_lead90to150_all = nanmean(SS_lick_all(label_rhythmic_all & label_lead90to150_all,:));
se_SS_lick_rhythmic_lead90to150_all = nanstd(SS_lick_all(label_rhythmic_all & label_lead90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lead90to150_all));

mean_SS_lick_rhythmic_lag_all = nanmean(SS_lick_all(label_rhythmic_all & label_lag_all,:));
se_SS_lick_rhythmic_lag_all = nanstd(SS_lick_all(label_rhythmic_all & label_lag_all,:))/sqrt(sum(label_rhythmic_all & label_lag_all));
mean_SS_lick_rhythmic_lead30to90_all = nanmean(SS_lick_all(label_rhythmic_all & label_lead30to90_all,:));
se_SS_lick_rhythmic_lead30to90_all = nanstd(SS_lick_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_SS_lick_rhythmic_lag90to150_all = nanmean(SS_lick_all(label_rhythmic_all & label_lag90to150_all,:));
se_SS_lick_rhythmic_lag90to150_all = nanstd(SS_lick_all(label_rhythmic_all & label_lag90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lag90to150_all));

mean_SS_lick_rhythmic_phase_l = nanmean(SS_lick_l(label_rhythmic_l & label_phase_l,:));
se_SS_lick_rhythmic_phase_l = nanstd(SS_lick_l(label_rhythmic_l & label_phase_l,:))/sqrt(sum(label_rhythmic_l & label_phase_l));

mean_SS_lick_rhythmic_antiphase_l = nanmean(SS_lick_l(label_rhythmic_l & label_antiphase_l,:));
se_SS_lick_rhythmic_antiphase_l = nanstd(SS_lick_l(label_rhythmic_l & label_antiphase_l,:))/sqrt(sum(label_rhythmic_l & label_antiphase_l));

mean_SS_lick_rhythmic_lead_l = nanmean(SS_lick_l(label_rhythmic_l & label_lead_l,:));
se_SS_lick_rhythmic_lead_l = nanstd(SS_lick_l(label_rhythmic_l & label_lead_l,:))/sqrt(sum(label_rhythmic_l & label_lead_l));
mean_SS_lick_rhythmic_lead30to90_l = nanmean(SS_lick_l(label_rhythmic_l & label_lead30to90_l,:));
se_SS_lick_rhythmic_lead30to90_l = nanstd(SS_lick_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_SS_lick_rhythmic_lead90to150_l = nanmean(SS_lick_l(label_rhythmic_l & label_lead90to150_l,:));
se_SS_lick_rhythmic_lead90to150_l = nanstd(SS_lick_l(label_rhythmic_l & label_lead90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lead90to150_l));

mean_SS_lick_rhythmic_lag_l = nanmean(SS_lick_l(label_rhythmic_l & label_lag_l,:));
se_SS_lick_rhythmic_lag_l = nanstd(SS_lick_l(label_rhythmic_l & label_lag_l,:))/sqrt(sum(label_rhythmic_l & label_lag_l));
mean_SS_lick_rhythmic_lead30to90_l = nanmean(SS_lick_l(label_rhythmic_l & label_lead30to90_l,:));
se_SS_lick_rhythmic_lead30to90_l = nanstd(SS_lick_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_SS_lick_rhythmic_lag90to150_l = nanmean(SS_lick_l(label_rhythmic_l & label_lag90to150_l,:));
se_SS_lick_rhythmic_lag90to150_l = nanstd(SS_lick_l(label_rhythmic_l & label_lag90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lag90to150_l));

mean_SS_lick_rhythmic_phase_g = nanmean(SS_lick_g(label_rhythmic_g & label_phase_g,:));
se_SS_lick_rhythmic_phase_g = nanstd(SS_lick_g(label_rhythmic_g & label_phase_g,:))/sqrt(sum(label_rhythmic_g & label_phase_g));

mean_SS_lick_rhythmic_antiphase_g = nanmean(SS_lick_g(label_rhythmic_g & label_antiphase_g,:));
se_SS_lick_rhythmic_antiphase_g = nanstd(SS_lick_g(label_rhythmic_g & label_antiphase_g,:))/sqrt(sum(label_rhythmic_g & label_antiphase_g));

mean_SS_lick_rhythmic_lead_g = nanmean(SS_lick_g(label_rhythmic_g & label_lead_g,:));
se_SS_lick_rhythmic_lead_g = nanstd(SS_lick_g(label_rhythmic_g & label_lead_g,:))/sqrt(sum(label_rhythmic_g & label_lead_g));
mean_SS_lick_rhythmic_lead30to90_g = nanmean(SS_lick_g(label_rhythmic_g & label_lead30to90_g,:));
se_SS_lick_rhythmic_lead30to90_g = nanstd(SS_lick_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_SS_lick_rhythmic_lead90to150_g = nanmean(SS_lick_g(label_rhythmic_g & label_lead90to150_g,:));
se_SS_lick_rhythmic_lead90to150_g = nanstd(SS_lick_g(label_rhythmic_g & label_lead90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lead90to150_g));

mean_SS_lick_rhythmic_lag_g = nanmean(SS_lick_g(label_rhythmic_g & label_lag_g,:));
se_SS_lick_rhythmic_lag_g = nanstd(SS_lick_g(label_rhythmic_g & label_lag_g,:))/sqrt(sum(label_rhythmic_g & label_lag_g));
mean_SS_lick_rhythmic_lead30to90_g = nanmean(SS_lick_g(label_rhythmic_g & label_lead30to90_g,:));
se_SS_lick_rhythmic_lead30to90_g = nanstd(SS_lick_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_SS_lick_rhythmic_lag90to150_g = nanmean(SS_lick_g(label_rhythmic_g & label_lag90to150_g,:));
se_SS_lick_rhythmic_lag90to150_g = nanstd(SS_lick_g(label_rhythmic_g & label_lag90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lag90to150_g));

mean_SS_lick_rhythmic_phase_r = nanmean(SS_lick_r(label_rhythmic_r & label_phase_r,:));
se_SS_lick_rhythmic_phase_r = nanstd(SS_lick_r(label_rhythmic_r & label_phase_r,:))/sqrt(sum(label_rhythmic_r & label_phase_r));

mean_SS_lick_rhythmic_antiphase_r = nanmean(SS_lick_r(label_rhythmic_r & label_antiphase_r,:));
se_SS_lick_rhythmic_antiphase_r = nanstd(SS_lick_r(label_rhythmic_r & label_antiphase_r,:))/sqrt(sum(label_rhythmic_r & label_antiphase_r));

mean_SS_lick_rhythmic_lead_r = nanmean(SS_lick_r(label_rhythmic_r & label_lead_r,:));
se_SS_lick_rhythmic_lead_r = nanstd(SS_lick_r(label_rhythmic_r & label_lead_r,:))/sqrt(sum(label_rhythmic_r & label_lead_r));
mean_SS_lick_rhythmic_lead30to90_r = nanmean(SS_lick_r(label_rhythmic_r & label_lead30to90_r,:));
se_SS_lick_rhythmic_lead30to90_r = nanstd(SS_lick_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_SS_lick_rhythmic_lead90to150_r = nanmean(SS_lick_r(label_rhythmic_r & label_lead90to150_r,:));
se_SS_lick_rhythmic_lead90to150_r = nanstd(SS_lick_r(label_rhythmic_r & label_lead90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lead90to150_r));

mean_SS_lick_rhythmic_lag_r = nanmean(SS_lick_r(label_rhythmic_r & label_lag_r,:));
se_SS_lick_rhythmic_lag_r = nanstd(SS_lick_r(label_rhythmic_r & label_lag_r,:))/sqrt(sum(label_rhythmic_r & label_lag_r));
mean_SS_lick_rhythmic_lead30to90_r = nanmean(SS_lick_r(label_rhythmic_r & label_lead30to90_r,:));
se_SS_lick_rhythmic_lead30to90_r = nanstd(SS_lick_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_SS_lick_rhythmic_lag90to150_r = nanmean(SS_lick_r(label_rhythmic_r & label_lag90to150_r,:));
se_SS_lick_rhythmic_lag90to150_r = nanstd(SS_lick_r(label_rhythmic_r & label_lag90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lag90to150_r));

%%% inds_span_corr %%%
inds_span_corr = -0.99 : 0.01 : 1;
%% PLOT LABELED TRACES
figure
subplot(3, 4, 1)
hold on
plot(inds_span_corr, mean_d_tip_rhythmic_phase_all(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_all + se_d_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_all - se_d_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_all(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_all + se_d_tip_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_all - se_d_tip_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lead_all(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_all + se_d_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_all - se_d_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lag_all(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_all + se_d_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_all - se_d_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-0.01 0.35])
ylabel('Displacement (mm)')
xlabel('Lick onset (ms)')
title('Lick trace | All licks')

subplot(3, 4, 2)
hold on
plot(inds_span_corr, mean_d_tip_rhythmic_phase_l(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_l + se_d_tip_rhythmic_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_l - se_d_tip_rhythmic_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_l(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_l + se_d_tip_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_l - se_d_tip_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lead_l(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_l + se_d_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_l - se_d_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lag_l(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_l + se_d_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_l - se_d_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-0.01 0.35])
ylabel('Displacement (mm)')
xlabel('Lick onset (ms)')
title('Lick trace | L licks')

subplot(3, 4, 3)
hold on
plot(inds_span_corr, mean_d_tip_rhythmic_phase_g(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_g + se_d_tip_rhythmic_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_g - se_d_tip_rhythmic_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_g(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_g + se_d_tip_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_g - se_d_tip_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lead_g(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_g + se_d_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_g - se_d_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lag_g(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_g + se_d_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_g - se_d_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-0.01 0.35])
ylabel('Displacement (mm)')
xlabel('Lick onset (ms)')
title('Lick trace | G licks')

subplot(3, 4, 4)
hold on
plot(inds_span_corr, mean_d_tip_rhythmic_phase_r(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_r + se_d_tip_rhythmic_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_r - se_d_tip_rhythmic_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_r(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_r + se_d_tip_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_r - se_d_tip_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lead_r(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_r + se_d_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_r - se_d_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lag_r(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_r + se_d_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_r - se_d_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-0.01 0.35])
ylabel('Displacement (mm)')
xlabel('Lick onset (ms)')
title('Lick trace | R licks')


subplot(3, 4, 5)
hold on
plot(inds_span_corr, mean_v_tip_rhythmic_phase_all(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_all + se_v_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_all - se_v_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_all(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_all + se_v_tip_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_all - se_v_tip_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lead_all(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_all +se_v_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_all - se_v_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lag_all(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_all + se_v_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_all - se_v_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
xlim([-0.01 0.35])
ylabel('Velocity (mm/s)')
xlabel('Lick onset (ms)')
title('Lick trace | All licks')

subplot(3, 4, 6)
hold on
plot(inds_span_corr, mean_v_tip_rhythmic_phase_l(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_l + se_v_tip_rhythmic_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_l - se_v_tip_rhythmic_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_l(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_l + se_v_tip_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_l - se_v_tip_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lead_l(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_l + se_v_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_l - se_v_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lag_l(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_l + se_v_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_l - se_v_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
xlim([-0.01 0.35])
ylabel('Velocity (mm/s)')
xlabel('Lick onset (ms)')
title('Lick trace | L licks')

subplot(3, 4, 7)
hold on
plot(inds_span_corr, mean_v_tip_rhythmic_phase_g(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_g + se_v_tip_rhythmic_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_g - se_v_tip_rhythmic_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_g(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_g + se_v_tip_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_g - se_v_tip_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lead_g(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_g + se_v_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_g - se_v_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lag_g(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_g + se_v_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_g - se_v_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
xlim([-0.01 0.35])
ylabel('Velocity (mm/s)')
xlabel('Lick onset (ms)')
title('Lick trace | G licks')

subplot(3, 4,8)
hold on
plot(inds_span_corr, mean_v_tip_rhythmic_phase_r(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_r + se_v_tip_rhythmic_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_r - se_v_tip_rhythmic_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_r(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_r + se_v_tip_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_r - se_v_tip_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lead_r(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_r + se_v_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_r - se_v_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lag_r(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_r + se_v_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_r - se_v_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
xlim([-0.01 0.35])
ylabel('Velocity (mm/s)')
xlabel('Lick onset (ms)')
title('Lick trace | R licks')

subplot(3, 4, 9)
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all + se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all - se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_all + se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_all - se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all + se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all - se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all + se_SS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all -se_SS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.01 0.35])
%     ylabel('Change in SS firing (spks/s)')
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick onset (ms)')
title('SS firing | All licks')

subplot(3, 4, 10)
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l + se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l - se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_l + se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_l - se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l + se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l - se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l + se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l - se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.01 0.35])
%     ylabel('Change in SS firing (spks/s)')
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick onset (ms)')
title('SS firing | L licks')

subplot(3, 4, 11)
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g + se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g - se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_g(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_g + se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_g - se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g + se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g - se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g + se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g - se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.01 0.35])
%     ylabel('Change in SS firing (spks/s)')
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick onset (ms)')
title('SS firing | G licks')

subplot(3, 4, 12)
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r + se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r - se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_r + se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_r - se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r + se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r - se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r + se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r - se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.01 0.35])
%     ylabel('Change in SS firing (spks/s)')
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick onset (ms)')
title('SS firing | R licks')
%% PLOT LABELED TRACES - Kinematics
figure
subplot(2, 4, 1)
hold on
yyaxis left
plot(inds_span_corr, mean_d_tip_all, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_all + se_d_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_all - se_d_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_all(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_all + se_d_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_all - se_d_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_all(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_all + se_d_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_all - se_d_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
ylabel('Displacement (mm)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all + se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all - se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all + se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all - se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all + se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all - se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all + se_SS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all - se_SS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.02 0.3])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick onset (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | All licks')

subplot(2, 4, 2)
hold on
yyaxis left
plot(inds_span_corr, mean_d_tip_l, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_l + se_d_tip_rhythmic_none_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_l - se_d_tip_rhythmic_none_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_l(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_l + se_d_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_l - se_d_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_l(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_l + se_d_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_l - se_d_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
ylabel('Displacement (mm)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l + se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l - se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l + se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l - se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l + se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l - se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l + se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l - se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.02 0.3])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick onset (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | L licks')


subplot(2, 4, 3)
hold on
yyaxis left
plot(inds_span_corr, mean_d_tip_g, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_g + se_d_tip_rhythmic_none_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_g - se_d_tip_rhythmic_none_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_g(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_g + se_d_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_g - se_d_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_g(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_g + se_d_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_g - se_d_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
ylabel('Displacement (mm)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g + se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g - se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_g(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_g + se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_g - se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g + se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g - se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g + se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g - se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.02 0.3])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick onset (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | G licks')


subplot(2, 4, 4)
hold on
yyaxis left
plot(inds_span_corr, mean_d_tip_r, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_r + se_d_tip_rhythmic_none_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_r - se_d_tip_rhythmic_none_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_r(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_r + se_d_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_r - se_d_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_r(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_r + se_d_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_r - se_d_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
ylabel('Displacement (mm)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r + se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r - se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r + se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r - se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r + se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r - se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r + se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r - se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.02 0.3])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick onset (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | R licks')

subplot(2, 4, 5)
hold on
yyaxis left
plot(inds_span_corr, mean_v_tip_all, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_all + se_v_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_all - se_v_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_all(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_all + se_v_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_all - se_v_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_all(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_all + se_v_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_all - se_v_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
ylabel('Velocity (mm/s)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all + se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all - se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all + se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all - se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all + se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all - se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all + seSS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all - se_SS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.02 0.3])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick onset (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | All licks')

subplot(2, 4, 6)
hold on
yyaxis left
plot(inds_span_corr, mean_v_tip_l, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_l + se_v_tip_rhythmic_none_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_l - se_v_tip_rhythmic_none_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_l(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_l + se_v_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_l - se_v_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_l(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_l + se_v_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_l - se_v_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
ylabel('Velocity (mm/s)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l + se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l - se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l + se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l - se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l + se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l - se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l + se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l - se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.02 0.3])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick onset (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | L licks')

subplot(2, 4, 7)
hold on
yyaxis left
plot(inds_span_corr, mean_v_tip_g, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_g + se_v_tip_rhythmic_none_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_g - se_v_tip_rhythmic_none_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_g(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_g + se_v_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_g - se_v_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_g(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_g + se_v_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_g - se_v_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
ylabel('Velocity (mm/s)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g + se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g - se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l + se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l - se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g + se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g - se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g + se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g - se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.02 0.3])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick onset (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | G licks')

subplot(2, 4, 8)
hold on
yyaxis left
plot(inds_span_corr, mean_v_tip_r, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_r + se_v_tip_rhythmic_none_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_r - se_v_tip_rhythmic_none_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_r(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_r + se_v_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_r - se_v_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_r(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_r + se_v_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_r - se_v_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
ylabel('Velocity (mm/s)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r + se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r - se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r + se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r - se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r + se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r - se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r + se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r - se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.02 0.3])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick onset (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | R licks')

%%  ANALYSIS PLOTS BASED ON LABELS - DMAX
clearvars -except ALL_PCELL_COMPRESSED_DATA POPULATION

number_pCell =  length(ALL_PCELL_COMPRESSED_DATA);

for counter_pCell = 1 : number_pCell
    id(counter_pCell, 1) = string(ALL_PCELL_COMPRESSED_DATA(counter_pCell).id{1, 1});
    
    %%% d_tip %%%
    d_tip_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_all_corr);
    d_tip_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_grooming_corr);
    d_tip_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_r_corr);
    d_tip_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_d_tip_100_l_corr);
    
    %%% v_tip %%%
    v_tip_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_v_tip_100_all_corr);
    v_tip_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_v_tip_100_grooming_corr);
    v_tip_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_v_tip_100_r_corr);
    v_tip_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.VID_v_tip_100_l_corr);
    
    %%% SS_str_bout %%%
    SS_str_bout_all(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_SS_str_bout * 100), 15, 'sgolay', 2);
    SS_str_bout_baseline_all(counter_pCell, :) = nanmean(SS_str_bout_all(counter_pCell, 1:50));
    SS_str_bout_r(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_SS_str_bout_r * 100), 15, 'sgolay', 2);
    SS_str_bout_baseline_r(counter_pCell, :) = nanmean(SS_str_bout_r(counter_pCell, 1:50));
    SS_str_bout_l(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_SS_str_bout_l * 100), 15, 'sgolay', 2);
    SS_str_bout_baseline_l(counter_pCell, :) = nanmean(SS_str_bout_l(counter_pCell, 1:50));
    
    %%% SS_lick %%%
    SS_lick_all(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_SS_all_corr * 100), 15, 'sgolay', 2);
    change_SS_lick_all(counter_pCell,:) = SS_lick_all(counter_pCell,:) - SS_str_bout_baseline_all(counter_pCell, 1);
    mean_window_SS_lick_all(counter_pCell,1) = nanmean(SS_lick_all(counter_pCell, :));
    mean_change_window_SS_lick_all(counter_pCell,1) = nanmean(change_SS_lick_all(counter_pCell, :));
    SS_lick_all(counter_pCell,:) = SS_lick_all(counter_pCell,:)/mean_window_SS_lick_all(counter_pCell,1);
    
    SS_lick_r(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_SS_r_corr * 100), 15, 'sgolay', 2);
    change_SS_lick_r(counter_pCell,:) = SS_lick_r(counter_pCell,:) - SS_str_bout_baseline_r(counter_pCell, 1);
    mean_window_SS_lick_r(counter_pCell,1) = nanmean(SS_lick_r(counter_pCell, :));
    mean_change_window_SS_lick_r(counter_pCell,1) = nanmean(change_SS_lick_r(counter_pCell, :));
    SS_lick_r(counter_pCell,:) = SS_lick_r(counter_pCell,:)/mean_window_SS_lick_r(counter_pCell,1);
    
    SS_lick_g(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_SS_grooming_corr * 100), 15, 'sgolay', 2);
    change_SS_lick_g(counter_pCell,:) = SS_lick_g(counter_pCell,:) - SS_str_bout_baseline_all(counter_pCell, 1);
    mean_window_SS_lick_g(counter_pCell,1) = nanmean(SS_lick_g(counter_pCell, :));
    mean_change_window_SS_lick_g(counter_pCell,1) = nanmean(change_SS_lick_g(counter_pCell, :));
    SS_lick_g(counter_pCell,:) = SS_lick_g(counter_pCell,:)/mean_window_SS_lick_g(counter_pCell,1);
    
    SS_lick_l(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_dmax.train_data_logic_SS_l_corr * 100), 15, 'sgolay', 2);
    change_SS_lick_l(counter_pCell,:) = SS_lick_l(counter_pCell,:) - SS_str_bout_baseline_l(counter_pCell, 1);
    mean_window_SS_lick_l(counter_pCell,1) = nanmean(SS_lick_l(counter_pCell, :));
    mean_change_window_SS_lick_l(counter_pCell,1) = nanmean(change_SS_lick_l(counter_pCell, :));
    SS_lick_l(counter_pCell,:) = SS_lick_l(counter_pCell,:)/mean_window_SS_lick_l(counter_pCell,1);
    
end

%%% LOAD LABELS FROM POPULATION STRUCTURE %%%
label_large_increase_all = POPULATION.dmax.label_large_increase_all;
label_medium_increase_all  = POPULATION.dmax.label_medium_increase_all;
label_small_increase_all = POPULATION.dmax.label_small_increase_all;
label_phase_all = POPULATION.dmax.label_phase_all;
label_small_decrease_all = POPULATION.dmax.label_small_decrease_all;
label_medium_decrease_all = POPULATION.dmax.label_medium_decrease_all;
label_large_decrease_all = POPULATION.dmax.label_large_decrease_all;
label_large_r = POPULATION.dmax.label_large_r;
label_medium_r = POPULATION.dmax.label_medium_r;
label_small_r = POPULATION.dmax.label_small_r;
label_none_dir = POPULATION.dmax.label_none_dir;
label_small_l = POPULATION.dmax.label_small_l;
label_medium_l = POPULATION.dmax.label_medium_l;
label_large_l = POPULATION.dmax.label_large_l;

label_rhythmic_all = POPULATION.dmax.label_rhythmic_all;
label_weakrhythmic_all = POPULATION.dmax.label_weakrhythmic_all;
label_rhythmic_l = POPULATION.dmax.label_rhythmic_l;
label_weakrhythmic_l = POPULATION.dmax.label_weakrhythmic_l;
label_rhythmic_g = POPULATION.dmax.label_rhythmic_g;
label_weakrhythmic_g = POPULATION.dmax.label_weakrhythmic_g;
label_rhythmic_r = POPULATION.dmax.label_rhythmic_r;
label_weakrhythmic_r = POPULATION.dmax.label_weakrhythmic_r;

label_lead_all = POPULATION.dmax.label_lead_all;
label_lead30to90_all  = POPULATION.dmax.label_lead30to90_all;
label_lead90to150_all  = POPULATION.dmax.label_lead90to150_all;
label_lag_all = POPULATION.dmax.label_lag_all;
label_lag30to90_all = POPULATION.dmax.label_lag30to90_all;
label_lag90to150_all = POPULATION.dmax.label_lag90to150_all;
label_phase_all = POPULATION.dmax.label_phase_all;
label_antiphase_all = POPULATION.dmax.label_antiphase_all;

label_lead_l = POPULATION.dmax.label_lead_l;
label_lead30to90_l  = POPULATION.dmax.label_lead30to90_l;
label_lead90to150_l  = POPULATION.dmax.label_lead90to150_l;
label_lag_l = POPULATION.dmax.label_lag_l;
label_lag30to90_l = POPULATION.dmax.label_lag30to90_l;
label_lag90to150_l = POPULATION.dmax.label_lag90to150_l;
label_phase_l = POPULATION.dmax.label_phase_l;
label_antiphase_l = POPULATION.dmax.label_antiphase_l;

label_lead_g = POPULATION.dmax.label_lead_g;
label_lead30to90_g  = POPULATION.dmax.label_lead30to90_g;
label_lead90to150_g  = POPULATION.dmax.label_lead90to150_g;
label_lag_g = POPULATION.dmax.label_lag_g;
label_lag30to90_g = POPULATION.dmax.label_lag30to90_g;
label_lag90to150_g = POPULATION.dmax.label_lag90to150_g;
label_phase_g = POPULATION.dmax.label_phase_g;
label_antiphase_g = POPULATION.dmax.label_antiphase_g;

label_lead_r = POPULATION.dmax.label_lead_r;
label_lead30to90_r  = POPULATION.dmax.label_lead30to90_r;
label_lead90to150_r  = POPULATION.dmax.label_lead90to150_r;
label_lag_r = POPULATION.dmax.label_lag_r;
label_lag30to90_r = POPULATION.dmax.label_lag30to90_r;
label_lag90to150_r = POPULATION.dmax.label_lag90to150_r;
label_phase_r = POPULATION.dmax.label_phase_r;
label_antiphase_r = POPULATION.dmax.label_antiphase_r;

%%% Mean d_tip based on labels
mean_d_tip_all = nanmean(d_tip_all);
mean_d_tip_l = nanmean(d_tip_l);
mean_d_tip_g = nanmean(d_tip_g);
mean_d_tip_r = nanmean(d_tip_r);

mean_d_tip_rhythmic_phase_all = nanmean(d_tip_all(label_rhythmic_all & label_phase_all,:));
se_d_tip_rhythmic_phase_all = nanstd(d_tip_all(label_rhythmic_all & label_phase_all,:))/sqrt(sum(label_rhythmic_all & label_phase_all));

mean_d_tip_rhythmic_antiphase_all = nanmean(d_tip_all(label_rhythmic_all & label_antiphase_all,:));
se_d_tip_rhythmic_antiphase_all = nanstd(d_tip_all(label_rhythmic_all & label_antiphase_all,:))/sqrt(sum(label_rhythmic_all & label_antiphase_all));

mean_d_tip_rhythmic_lead_all = nanmean(d_tip_all(label_rhythmic_all & label_lead_all,:));
se_d_tip_rhythmic_lead_all = nanstd(d_tip_all(label_rhythmic_all & label_lead_all,:))/sqrt(sum(label_rhythmic_all & label_lead_all));
mean_d_tip_rhythmic_lead30to90_all = nanmean(d_tip_all(label_rhythmic_all & label_lead30to90_all,:));
se_d_tip_rhythmic_lead30to90_all = nanstd(d_tip_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_d_tip_rhythmic_lead90to150_all = nanmean(d_tip_all(label_rhythmic_all & label_lead90to150_all,:));
se_d_tip_rhythmic_lead90to150_all = nanstd(d_tip_all(label_rhythmic_all & label_lead90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lead90to150_all));

mean_d_tip_rhythmic_lag_all = nanmean(d_tip_all(label_rhythmic_all & label_lag_all,:));
se_d_tip_rhythmic_lag_all = nanstd(d_tip_all(label_rhythmic_all & label_lag_all,:))/sqrt(sum(label_rhythmic_all & label_lag_all));
mean_d_tip_rhythmic_lead30to90_all = nanmean(d_tip_all(label_rhythmic_all & label_lead30to90_all,:));
se_d_tip_rhythmic_lead30to90_all = nanstd(d_tip_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_d_tip_rhythmic_lag90to150_all = nanmean(d_tip_all(label_rhythmic_all & label_lag90to150_all,:));
se_d_tip_rhythmic_lag90to150_all = nanstd(d_tip_all(label_rhythmic_all & label_lag90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lag90to150_all));

mean_d_tip_rhythmic_phase_l = nanmean(d_tip_l(label_rhythmic_l & label_phase_l,:));
se_d_tip_rhythmic_phase_l = nanstd(d_tip_l(label_rhythmic_l & label_phase_l,:))/sqrt(sum(label_rhythmic_l & label_phase_l));

mean_d_tip_rhythmic_antiphase_l = nanmean(d_tip_l(label_rhythmic_l & label_antiphase_l,:));
se_d_tip_rhythmic_antiphase_l = nanstd(d_tip_l(label_rhythmic_l & label_antiphase_l,:))/sqrt(sum(label_rhythmic_l & label_antiphase_l));

mean_d_tip_rhythmic_lead_l = nanmean(d_tip_l(label_rhythmic_l & label_lead_l,:));
se_d_tip_rhythmic_lead_l = nanstd(d_tip_l(label_rhythmic_l & label_lead_l,:))/sqrt(sum(label_rhythmic_l & label_lead_l));
mean_d_tip_rhythmic_lead30to90_l = nanmean(d_tip_l(label_rhythmic_l & label_lead30to90_l,:));
se_d_tip_rhythmic_lead30to90_l = nanstd(d_tip_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_d_tip_rhythmic_lead90to150_l = nanmean(d_tip_l(label_rhythmic_l & label_lead90to150_l,:));
se_d_tip_rhythmic_lead90to150_l = nanstd(d_tip_l(label_rhythmic_l & label_lead90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lead90to150_l));

mean_d_tip_rhythmic_lag_l = nanmean(d_tip_l(label_rhythmic_l & label_lag_l,:));
se_d_tip_rhythmic_lag_l = nanstd(d_tip_l(label_rhythmic_l & label_lag_l,:))/sqrt(sum(label_rhythmic_l & label_lag_l));
mean_d_tip_rhythmic_lead30to90_l = nanmean(d_tip_l(label_rhythmic_l & label_lead30to90_l,:));
se_d_tip_rhythmic_lead30to90_l = nanstd(d_tip_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_d_tip_rhythmic_lag90to150_l = nanmean(d_tip_l(label_rhythmic_l & label_lag90to150_l,:));
se_d_tip_rhythmic_lag90to150_l = nanstd(d_tip_l(label_rhythmic_l & label_lag90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lag90to150_l));

mean_d_tip_rhythmic_phase_g = nanmean(d_tip_g(label_rhythmic_g & label_phase_g,:));
se_d_tip_rhythmic_phase_g = nanstd(d_tip_g(label_rhythmic_g & label_phase_g,:))/sqrt(sum(label_rhythmic_g & label_phase_g));

mean_d_tip_rhythmic_antiphase_g = nanmean(d_tip_g(label_rhythmic_g & label_antiphase_g,:));
se_d_tip_rhythmic_antiphase_g = nanstd(d_tip_g(label_rhythmic_g & label_antiphase_g,:))/sqrt(sum(label_rhythmic_g & label_antiphase_g));

mean_d_tip_rhythmic_lead_g = nanmean(d_tip_g(label_rhythmic_g & label_lead_g,:));
se_d_tip_rhythmic_lead_g = nanstd(d_tip_g(label_rhythmic_g & label_lead_g,:))/sqrt(sum(label_rhythmic_g & label_lead_g));
mean_d_tip_rhythmic_lead30to90_g = nanmean(d_tip_g(label_rhythmic_g & label_lead30to90_g,:));
se_d_tip_rhythmic_lead30to90_g = nanstd(d_tip_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_d_tip_rhythmic_lead90to150_g = nanmean(d_tip_g(label_rhythmic_g & label_lead90to150_g,:));
se_d_tip_rhythmic_lead90to150_g = nanstd(d_tip_g(label_rhythmic_g & label_lead90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lead90to150_g));

mean_d_tip_rhythmic_lag_g = nanmean(d_tip_g(label_rhythmic_g & label_lag_g,:));
se_d_tip_rhythmic_lag_g = nanstd(d_tip_g(label_rhythmic_g & label_lag_g,:))/sqrt(sum(label_rhythmic_g & label_lag_g));
mean_d_tip_rhythmic_lead30to90_g = nanmean(d_tip_g(label_rhythmic_g & label_lead30to90_g,:));
se_d_tip_rhythmic_lead30to90_g = nanstd(d_tip_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_d_tip_rhythmic_lag90to150_g = nanmean(d_tip_g(label_rhythmic_g & label_lag90to150_g,:));
se_d_tip_rhythmic_lag90to150_g = nanstd(d_tip_g(label_rhythmic_g & label_lag90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lag90to150_g));

mean_d_tip_rhythmic_phase_r = nanmean(d_tip_r(label_rhythmic_r & label_phase_r,:));
se_d_tip_rhythmic_phase_r = nanstd(d_tip_r(label_rhythmic_r & label_phase_r,:))/sqrt(sum(label_rhythmic_r & label_phase_r));

mean_d_tip_rhythmic_antiphase_r = nanmean(d_tip_g(label_rhythmic_r & label_antiphase_r,:));
se_d_tip_rhythmic_antiphase_r = nanstd(d_tip_g(label_rhythmic_r & label_antiphase_r,:))/sqrt(sum(label_rhythmic_r & label_antiphase_r));

mean_d_tip_rhythmic_lead_r = nanmean(d_tip_r(label_rhythmic_r & label_lead_r,:));
se_d_tip_rhythmic_lead_r = nanstd(d_tip_r(label_rhythmic_r & label_lead_r,:))/sqrt(sum(label_rhythmic_r & label_lead_r));
mean_d_tip_rhythmic_lead30to90_r = nanmean(d_tip_r(label_rhythmic_r & label_lead30to90_r,:));
se_d_tip_rhythmic_lead30to90_r = nanstd(d_tip_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_d_tip_rhythmic_lead90to150_r = nanmean(d_tip_r(label_rhythmic_r & label_lead90to150_r,:));
se_d_tip_rhythmic_lead90to150_r = nanstd(d_tip_r(label_rhythmic_r & label_lead90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lead90to150_r));

mean_d_tip_rhythmic_lag_r = nanmean(d_tip_r(label_rhythmic_r & label_lag_r,:));
se_d_tip_rhythmic_lag_r = nanstd(d_tip_r(label_rhythmic_r & label_lag_r,:))/sqrt(sum(label_rhythmic_r & label_lag_r));
mean_d_tip_rhythmic_lead30to90_r = nanmean(d_tip_r(label_rhythmic_r & label_lead30to90_r,:));
se_d_tip_rhythmic_lead30to90_r = nanstd(d_tip_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_d_tip_rhythmic_lag90to150_r = nanmean(d_tip_r(label_rhythmic_r & label_lag90to150_r,:));
se_d_tip_rhythmic_lag90to150_r = nanstd(d_tip_r(label_rhythmic_r & label_lag90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lag90to150_r));

%%% Mean v_tip based on labels
mean_v_tip_all = nanmean(v_tip_all);
mean_v_tip_l = nanmean(v_tip_l);
mean_v_tip_g = nanmean(v_tip_g);
mean_v_tip_r = nanmean(v_tip_r);

mean_v_tip_rhythmic_phase_all = nanmean(v_tip_all(label_rhythmic_all & label_phase_all,:));
se_v_tip_rhythmic_phase_all = nanstd(v_tip_all(label_rhythmic_all & label_phase_all,:))/sqrt(sum(label_rhythmic_all & label_phase_all));

mean_v_tip_rhythmic_antiphase_all = nanmean(v_tip_all(label_rhythmic_all & label_antiphase_all,:));
se_v_tip_rhythmic_antiphase_all = nanstd(v_tip_all(label_rhythmic_all & label_antiphase_all,:))/sqrt(sum(label_rhythmic_all & label_antiphase_all));

mean_v_tip_rhythmic_lead_all = nanmean(v_tip_all(label_rhythmic_all & label_lead_all,:));
se_v_tip_rhythmic_lead_all = nanstd(v_tip_all(label_rhythmic_all & label_lead_all,:))/sqrt(sum(label_rhythmic_all & label_lead_all));
mean_v_tip_rhythmic_lead30to90_all = nanmean(v_tip_all(label_rhythmic_all & label_lead30to90_all,:));
se_v_tip_rhythmic_lead30to90_all = nanstd(v_tip_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_v_tip_rhythmic_lead90to150_all = nanmean(v_tip_all(label_rhythmic_all & label_lead90to150_all,:));
se_v_tip_rhythmic_lead90to150_all = nanstd(v_tip_all(label_rhythmic_all & label_lead90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lead90to150_all));

mean_v_tip_rhythmic_lag_all = nanmean(v_tip_all(label_rhythmic_all & label_lag_all,:));
se_v_tip_rhythmic_lag_all = nanstd(v_tip_all(label_rhythmic_all & label_lag_all,:))/sqrt(sum(label_rhythmic_all & label_lag_all));
mean_v_tip_rhythmic_lead30to90_all = nanmean(v_tip_all(label_rhythmic_all & label_lead30to90_all,:));
se_v_tip_rhythmic_lead30to90_all = nanstd(v_tip_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_v_tip_rhythmic_lag90to150_all = nanmean(v_tip_all(label_rhythmic_all & label_lag90to150_all,:));
se_v_tip_rhythmic_lag90to150_all = nanstd(v_tip_all(label_rhythmic_all & label_lag90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lag90to150_all));

mean_v_tip_rhythmic_phase_l = nanmean(v_tip_l(label_rhythmic_l & label_phase_l,:));
se_v_tip_rhythmic_phasell = nanstd(v_tip_l(label_rhythmic_l & label_phase_l,:))/sqrt(sum(label_rhythmic_l & label_phase_l));

mean_v_tip_rhythmic_antiphase_l = nanmean(v_tip_l(label_rhythmic_l & label_antiphase_l,:));
se_v_tip_rhythmic_antiphase_l = nanstd(v_tip_l(label_rhythmic_l & label_antiphase_l,:))/sqrt(sum(label_rhythmic_l & label_antiphase_l));

mean_v_tip_rhythmic_lead_l = nanmean(v_tip_l(label_rhythmic_l & label_lead_l,:));
se_v_tip_rhythmic_lead_l = nanstd(v_tip_l(label_rhythmic_l & label_lead_l,:))/sqrt(sum(label_rhythmic_l & label_lead_l));
mean_v_tip_rhythmic_lead30to90_l = nanmean(v_tip_l(label_rhythmic_l & label_lead30to90_l,:));
se_v_tip_rhythmic_lead30to90_l = nanstd(v_tip_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_v_tip_rhythmic_lead90to150_l = nanmean(v_tip_l(label_rhythmic_l & label_lead90to150_l,:));
se_v_tip_rhythmic_lead90to150_l = nanstd(v_tip_l(label_rhythmic_l & label_lead90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lead90to150_l));

mean_v_tip_rhythmic_lag_l = nanmean(v_tip_l(label_rhythmic_l & label_lag_l,:));
se_v_tip_rhythmic_lag_l = nanstd(v_tip_l(label_rhythmic_l & label_lag_l,:))/sqrt(sum(label_rhythmic_l & label_lag_l));
mean_v_tip_rhythmic_lead30to90_l = nanmean(v_tip_l(label_rhythmic_l & label_lead30to90_l,:));
se_v_tip_rhythmic_lead30to90_l = nanstd(v_tip_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_v_tip_rhythmic_lag90to150_l = nanmean(v_tip_l(label_rhythmic_l & label_lag90to150_l,:));
se_v_tip_rhythmic_lag90to150_l = nanstd(v_tip_l(label_rhythmic_l & label_lag90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lag90to150_l));

mean_v_tip_rhythmic_phase_g = nanmean(v_tip_g(label_rhythmic_g & label_phase_g,:));
se_v_tip_rhythmic_phase_g = nanstd(v_tip_g(label_rhythmic_g & label_phase_g,:))/sqrt(sum(label_rhythmic_g & label_phase_g));

mean_v_tip_rhythmic_antiphase_g = nanmean(v_tip_g(label_rhythmic_g & label_antiphase_g,:));
se_v_tip_rhythmic_antiphase_g = nanstd(v_tip_g(label_rhythmic_g & label_antiphase_g,:))/sqrt(sum(label_rhythmic_g & label_antiphase_g));

mean_v_tip_rhythmic_lead_g = nanmean(v_tip_g(label_rhythmic_g & label_lead_g,:));
se_v_tip_rhythmic_lead_g = nanstd(v_tip_g(label_rhythmic_g & label_lead_g,:))/sqrt(sum(label_rhythmic_g & label_lead_g));
mean_v_tip_rhythmic_lead30to90_g = nanmean(v_tip_g(label_rhythmic_g & label_lead30to90_g,:));
se_v_tip_rhythmic_lead30to90_g = nanstd(v_tip_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_v_tip_rhythmic_lead90to150_g = nanmean(v_tip_g(label_rhythmic_g & label_lead90to150_g,:));
se_v_tip_rhythmic_lead90to150_g = nanstd(v_tip_g(label_rhythmic_g & label_lead90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lead90to150_g));

mean_v_tip_rhythmic_lag_g = nanmean(v_tip_g(label_rhythmic_g & label_lag_g,:));
se_v_tip_rhythmic_lag_g = nanstd(v_tip_g(label_rhythmic_g & label_lag_g,:))/sqrt(sum(label_rhythmic_g & label_lag_g));
mean_v_tip_rhythmic_lead30to90_g = nanmean(v_tip_g(label_rhythmic_g & label_lead30to90_g,:));
se_v_tip_rhythmic_lead30to90_g = nanstd(v_tip_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_v_tip_rhythmic_lag90to150_g = nanmean(v_tip_g(label_rhythmic_g & label_lag90to150_g,:));
se_v_tip_rhythmic_lag90to150_g = nanstd(v_tip_g(label_rhythmic_g & label_lag90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lag90to150_g));

mean_v_tip_rhythmic_phase_r = nanmean(v_tip_r(label_rhythmic_r & label_phase_r,:));
se_v_tip_rhythmic_phase_r = nanstd(v_tip_r(label_rhythmic_r & label_phase_r,:))/sqrt(sum(label_rhythmic_r & label_phase_r));

mean_v_tip_rhythmic_antiphase_r = nanmean(v_tip_r(label_rhythmic_r & label_antiphase_r,:));
se_v_tip_rhythmic_antiphaser = nanstd(v_tip_r(label_rhythmic_r & label_antiphase_r,:))/sqrt(sum(label_rhythmic_r & label_antiphase_r));

mean_v_tip_rhythmic_lead_r = nanmean(v_tip_r(label_rhythmic_r & label_lead_r,:));
se_v_tip_rhythmic_lead_r = nanstd(v_tip_r(label_rhythmic_r & label_lead_r,:))/sqrt(sum(label_rhythmic_r & label_lead_r));
mean_v_tip_rhythmic_lead30to90_r = nanmean(v_tip_r(label_rhythmic_r & label_lead30to90_r,:));
se_v_tip_rhythmic_lead30to90_r = nanstd(v_tip_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_v_tip_rhythmic_lead90to150_r = nanmean(v_tip_r(label_rhythmic_r & label_lead90to150_r,:));
se_v_tip_rhythmic_lead90to150_r = nanstd(v_tip_r(label_rhythmic_r & label_lead90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lead90to150_r));

mean_v_tip_rhythmic_lag_r = nanmean(v_tip_r(label_rhythmic_r & label_lag_r,:));
se_v_tip_rhythmic_lag_r = nanstd(v_tip_r(label_rhythmic_r & label_lag_r,:))/sqrt(sum(label_rhythmic_r & label_lag_r));
mean_v_tip_rhythmic_lead30to90_r = nanmean(v_tip_r(label_rhythmic_r & label_lead30to90_r,:));
se_v_tip_rhythmic_lead30to90_r = nanstd(v_tip_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_v_tip_rhythmic_lag90to150_r = nanmean(v_tip_r(label_rhythmic_r & label_lag90to150_r,:));
se_v_tip_rhythmic_lag90to150_r = nanstd(v_tip_r(label_rhythmic_r & label_lag90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lag90to150_r));

%%% Mean change_SS_lick based on labels
mean_SS_lick_rhythmic_phase_all = nanmean(SS_lick_all(label_rhythmic_all & label_phase_all,:));
se_SS_lick_rhythmic_phase_all = nanstd(SS_lick_all(label_rhythmic_all & label_phase_all,:))/sqrt(sum(label_rhythmic_all & label_phase_all));

mean_SS_lick_rhythmic_antiphase_all = nanmean(SS_lick_all(label_rhythmic_all & label_antiphase_all,:));
se_SS_lick_rhythmic_antiphase_all = nanstd(SS_lick_all(label_rhythmic_all & label_antiphase_all,:))/sqrt(sum(label_rhythmic_all & label_antiphase_all));

mean_SS_lick_rhythmic_lead_all = nanmean(SS_lick_all(label_rhythmic_all & label_lead_all,:));
se_SS_lick_rhythmic_lead_all = nanstd(SS_lick_all(label_rhythmic_all & label_lead_all,:))/sqrt(sum(label_rhythmic_all & label_lead_all));
mean_SS_lick_rhythmic_lead30to90_all = nanmean(SS_lick_all(label_rhythmic_all & label_lead30to90_all,:));
se_SS_lick_rhythmic_lead30to90_all = nanstd(SS_lick_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_SS_lick_rhythmic_lead90to150_all = nanmean(SS_lick_all(label_rhythmic_all & label_lead90to150_all,:));
se_SS_lick_rhythmic_lead90to150_all = nanstd(SS_lick_all(label_rhythmic_all & label_lead90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lead90to150_all));

mean_SS_lick_rhythmic_lag_all = nanmean(SS_lick_all(label_rhythmic_all & label_lag_all,:));
se_SS_lick_rhythmic_lag_all = nanstd(SS_lick_all(label_rhythmic_all & label_lag_all,:))/sqrt(sum(label_rhythmic_all & label_lag_all));
mean_SS_lick_rhythmic_lead30to90_all = nanmean(SS_lick_all(label_rhythmic_all & label_lead30to90_all,:));
se_SS_lick_rhythmic_lead30to90_all = nanstd(SS_lick_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_SS_lick_rhythmic_lag90to150_all = nanmean(SS_lick_all(label_rhythmic_all & label_lag90to150_all,:));
se_SS_lick_rhythmic_lag90to150_all = nanstd(SS_lick_all(label_rhythmic_all & label_lag90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lag90to150_all));

mean_SS_lick_rhythmic_phase_l = nanmean(SS_lick_l(label_rhythmic_l & label_phase_l,:));
se_SS_lick_rhythmic_phase_l = nanstd(SS_lick_l(label_rhythmic_l & label_phase_l,:))/sqrt(sum(label_rhythmic_l & label_phase_l));

mean_SS_lick_rhythmic_antiphase_l = nanmean(SS_lick_l(label_rhythmic_l & label_antiphase_l,:));
se_SS_lick_rhythmic_antiphase_l = nanstd(SS_lick_l(label_rhythmic_l & label_antiphase_l,:))/sqrt(sum(label_rhythmic_l & label_antiphase_l));

mean_SS_lick_rhythmic_lead_l = nanmean(SS_lick_l(label_rhythmic_l & label_lead_l,:));
se_SS_lick_rhythmic_lead_l = nanstd(SS_lick_l(label_rhythmic_l & label_lead_l,:))/sqrt(sum(label_rhythmic_l & label_lead_l));
mean_SS_lick_rhythmic_lead30to90_l = nanmean(SS_lick_l(label_rhythmic_l & label_lead30to90_l,:));
se_SS_lick_rhythmic_lead30to90_l = nanstd(SS_lick_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_SS_lick_rhythmic_lead90to150_l = nanmean(SS_lick_l(label_rhythmic_l & label_lead90to150_l,:));
se_SS_lick_rhythmic_lead90to150_l = nanstd(SS_lick_l(label_rhythmic_l & label_lead90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lead90to150_l));

mean_SS_lick_rhythmic_lag_l = nanmean(SS_lick_l(label_rhythmic_l & label_lag_l,:));
se_SS_lick_rhythmic_lag_l = nanstd(SS_lick_l(label_rhythmic_l & label_lag_l,:))/sqrt(sum(label_rhythmic_l & label_lag_l));
mean_SS_lick_rhythmic_lead30to90_l = nanmean(SS_lick_l(label_rhythmic_l & label_lead30to90_l,:));
se_SS_lick_rhythmic_lead30to90_l = nanstd(SS_lick_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_SS_lick_rhythmic_lag90to150_l = nanmean(SS_lick_l(label_rhythmic_l & label_lag90to150_l,:));
se_SS_lick_rhythmic_lag90to150_l = nanstd(SS_lick_l(label_rhythmic_l & label_lag90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lag90to150_l));

mean_SS_lick_rhythmic_phase_g = nanmean(SS_lick_g(label_rhythmic_g & label_phase_g,:));
se_SS_lick_rhythmic_phase_g = nanstd(SS_lick_g(label_rhythmic_g & label_phase_g,:))/sqrt(sum(label_rhythmic_g & label_phase_g));

mean_SS_lick_rhythmic_antiphase_g = nanmean(SS_lick_g(label_rhythmic_g & label_antiphase_g,:));
se_SS_lick_rhythmic_antiphase_g = nanstd(SS_lick_g(label_rhythmic_g & label_antiphase_g,:))/sqrt(sum(label_rhythmic_g & label_antiphase_g));

mean_SS_lick_rhythmic_lead_g = nanmean(SS_lick_g(label_rhythmic_g & label_lead_g,:));
se_SS_lick_rhythmic_lead_g = nanstd(SS_lick_g(label_rhythmic_g & label_lead_g,:))/sqrt(sum(label_rhythmic_g & label_lead_g));
mean_SS_lick_rhythmic_lead30to90_g = nanmean(SS_lick_g(label_rhythmic_g & label_lead30to90_g,:));
se_SS_lick_rhythmic_lead30to90_g = nanstd(SS_lick_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_SS_lick_rhythmic_lead90to150_g = nanmean(SS_lick_g(label_rhythmic_g & label_lead90to150_g,:));
se_SS_lick_rhythmic_lead90to150_g = nanstd(SS_lick_g(label_rhythmic_g & label_lead90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lead90to150_g));

mean_SS_lick_rhythmic_lag_g = nanmean(SS_lick_g(label_rhythmic_g & label_lag_g,:));
se_SS_lick_rhythmic_lag_g = nanstd(SS_lick_g(label_rhythmic_g & label_lag_g,:))/sqrt(sum(label_rhythmic_g & label_lag_g));
mean_SS_lick_rhythmic_lead30to90_g = nanmean(SS_lick_g(label_rhythmic_g & label_lead30to90_g,:));
se_SS_lick_rhythmic_lead30to90_g = nanstd(SS_lick_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_SS_lick_rhythmic_lag90to150_g = nanmean(SS_lick_g(label_rhythmic_g & label_lag90to150_g,:));
se_SS_lick_rhythmic_lag90to150_g = nanstd(SS_lick_g(label_rhythmic_g & label_lag90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lag90to150_g));

mean_SS_lick_rhythmic_phase_r = nanmean(SS_lick_r(label_rhythmic_r & label_phase_r,:));
se_SS_lick_rhythmic_phase_r = nanstd(SS_lick_r(label_rhythmic_r & label_phase_r,:))/sqrt(sum(label_rhythmic_r & label_phase_r));

mean_SS_lick_rhythmic_antiphase_r = nanmean(SS_lick_r(label_rhythmic_r & label_antiphase_r,:));
se_SS_lick_rhythmic_antiphase_r = nanstd(SS_lick_r(label_rhythmic_r & label_antiphase_r,:))/sqrt(sum(label_rhythmic_r & label_antiphase_r));

mean_SS_lick_rhythmic_lead_r = nanmean(SS_lick_r(label_rhythmic_r & label_lead_r,:));
se_SS_lick_rhythmic_lead_r = nanstd(SS_lick_r(label_rhythmic_r & label_lead_r,:))/sqrt(sum(label_rhythmic_r & label_lead_r));
mean_SS_lick_rhythmic_lead30to90_r = nanmean(SS_lick_r(label_rhythmic_r & label_lead30to90_r,:));
se_SS_lick_rhythmic_lead30to90_r = nanstd(SS_lick_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_SS_lick_rhythmic_lead90to150_r = nanmean(SS_lick_r(label_rhythmic_r & label_lead90to150_r,:));
se_SS_lick_rhythmic_lead90to150_r = nanstd(SS_lick_r(label_rhythmic_r & label_lead90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lead90to150_r));

mean_SS_lick_rhythmic_lag_r = nanmean(SS_lick_r(label_rhythmic_r & label_lag_r,:));
se_SS_lick_rhythmic_lag_r = nanstd(SS_lick_r(label_rhythmic_r & label_lag_r,:))/sqrt(sum(label_rhythmic_r & label_lag_r));
mean_SS_lick_rhythmic_lead30to90_r = nanmean(SS_lick_r(label_rhythmic_r & label_lead30to90_r,:));
se_SS_lick_rhythmic_lead30to90_r = nanstd(SS_lick_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_SS_lick_rhythmic_lag90to150_r = nanmean(SS_lick_r(label_rhythmic_r & label_lag90to150_r,:));
se_SS_lick_rhythmic_lag90to150_r = nanstd(SS_lick_r(label_rhythmic_r & label_lag90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lag90to150_r));

%%% inds_span_corr %%%
inds_span_corr = -0.99 : 0.01 : 1;
%% PLOT LABELED TRACES
figure
subplot(3, 4, 1)
hold on
plot(inds_span_corr, mean_d_tip_rhythmic_phase_all(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_all + se_d_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_all - se_d_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_all(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_all + se_d_tip_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_all - se_d_tip_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lead_all(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_all + se_d_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_all - se_d_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lag_all(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_all + se_d_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_all - se_d_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-0.17 0.17])
ylabel('Displacement (mm)')
xlabel('Lick dmax (ms)')
title('Lick trace | All licks')

subplot(3, 4, 2)
hold on
plot(inds_span_corr, mean_d_tip_rhythmic_phase_l(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_l + se_d_tip_rhythmic_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_l - se_d_tip_rhythmic_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_l(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_l + se_d_tip_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_l - se_d_tip_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lead_l(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_l + se_d_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_l - se_d_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lag_l(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_l + se_d_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_l - se_d_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-0.17 0.17])
ylabel('Displacement (mm)')
xlabel('Lick dmax (ms)')
title('Lick trace | L licks')

subplot(3, 4, 3)
hold on
plot(inds_span_corr, mean_d_tip_rhythmic_phase_g(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_g + se_d_tip_rhythmic_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_g - se_d_tip_rhythmic_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_g(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_g + se_d_tip_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_g - se_d_tip_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lead_g(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_g + se_d_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_g - se_d_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lag_g(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_g + se_d_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_g - se_d_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-0.17 0.17])
ylabel('Displacement (mm)')
xlabel('Lick dmax (ms)')
title('Lick trace | G licks')

subplot(3, 4, 4)
hold on
plot(inds_span_corr, mean_d_tip_rhythmic_phase_r(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_r + se_d_tip_rhythmic_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_r - se_d_tip_rhythmic_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_r(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_r + se_d_tip_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_r - se_d_tip_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lead_r(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_r + se_d_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_r - se_d_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lag_r(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_r + se_d_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_r - se_d_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-0.17 0.17])
ylabel('Displacement (mm)')
xlabel('Lick dmax (ms)')
title('Lick trace | R licks')


subplot(3, 4, 5)
hold on
plot(inds_span_corr, mean_v_tip_rhythmic_phase_all(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_all + se_v_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_all - se_v_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_all(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_all + se_v_tip_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_all - se_v_tip_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lead_all(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_all +se_v_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_all - se_v_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lag_all(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_all + se_v_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_all - se_v_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
xlim([-0.17 0.17])
ylabel('Velocity (mm/s)')
xlabel('Lick dmax (ms)')
title('Lick trace | All licks')

subplot(3, 4, 6)
hold on
plot(inds_span_corr, mean_v_tip_rhythmic_phase_l(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_l + se_v_tip_rhythmic_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_l - se_v_tip_rhythmic_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_l(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_l + se_v_tip_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_l - se_v_tip_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lead_l(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_l + se_v_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_l - se_v_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lag_l(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_l + se_v_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_l - se_v_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
xlim([-0.17 0.17])
ylabel('Velocity (mm/s)')
xlabel('Lick dmax (ms)')
title('Lick trace | L licks')

subplot(3, 4, 7)
hold on
plot(inds_span_corr, mean_v_tip_rhythmic_phase_g(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_g + se_v_tip_rhythmic_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_g - se_v_tip_rhythmic_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_g(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_g + se_v_tip_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_g - se_v_tip_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lead_g(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_g + se_v_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_g - se_v_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lag_g(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_g + se_v_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_g - se_v_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
xlim([-0.17 0.17])
ylabel('Velocity (mm/s)')
xlabel('Lick dmax (ms)')
title('Lick trace | G licks')

subplot(3, 4,8)
hold on
plot(inds_span_corr, mean_v_tip_rhythmic_phase_r(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_r + se_v_tip_rhythmic_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_r - se_v_tip_rhythmic_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_r(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_r + se_v_tip_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_r - se_v_tip_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lead_r(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_r + se_v_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_r - se_v_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lag_r(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_r + se_v_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_r - se_v_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
xlim([-0.17 0.17])
ylabel('Velocity (mm/s)')
xlabel('Lick dmax (ms)')
title('Lick trace | R licks')

subplot(3, 4, 9)
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all + se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all - se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_all + se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_all - se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all + se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all - se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all + se_SS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all -se_SS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
%     ylabel('Change in SS firing (spks/s)')
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick dmax (ms)')
title('SS firing | All licks')

subplot(3, 4, 10)
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l + se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l - se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_l + se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_l - se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l + se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l - se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l + se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l - se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
%     ylabel('Change in SS firing (spks/s)')
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick dmax (ms)')
title('SS firing | L licks')

subplot(3, 4, 11)
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g + se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g - se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_g(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_g + se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_g - se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g + se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g - se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g + se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g - se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
%     ylabel('Change in SS firing (spks/s)')
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick dmax (ms)')
title('SS firing | G licks')

subplot(3, 4, 12)
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r + se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r - se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_r + se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_r - se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r + se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r - se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r + se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r - se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
%     ylabel('Change in SS firing (spks/s)')
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick dmax (ms)')
title('SS firing | R licks')
%% PLOT LABELED TRACES - Kinematics
figure
subplot(2, 4, 1)
hold on
yyaxis left
plot(inds_span_corr, mean_d_tip_all, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_all + se_d_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_all - se_d_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_all(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_all + se_d_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_all - se_d_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_all(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_all + se_d_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_all - se_d_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
ylabel('Displacement (mm)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all + se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all - se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all + se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all - se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all + se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all - se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all + se_SS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all - se_SS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick dmax (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | All licks')

subplot(2, 4, 2)
hold on
yyaxis left
plot(inds_span_corr, mean_d_tip_l, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_l + se_d_tip_rhythmic_none_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_l - se_d_tip_rhythmic_none_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_l(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_l + se_d_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_l - se_d_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_l(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_l + se_d_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_l - se_d_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
ylabel('Displacement (mm)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l + se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l - se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l + se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l - se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l + se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l - se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l + se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l - se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick dmax (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | L licks')


subplot(2, 4, 3)
hold on
yyaxis left
plot(inds_span_corr, mean_d_tip_g, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_g + se_d_tip_rhythmic_none_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_g - se_d_tip_rhythmic_none_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_g(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_g + se_d_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_g - se_d_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_g(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_g + se_d_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_g - se_d_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
ylabel('Displacement (mm)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g + se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g - se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_g(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_g + se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_g - se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g + se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g - se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g + se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g - se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick dmax (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | G licks')


subplot(2, 4, 4)
hold on
yyaxis left
plot(inds_span_corr, mean_d_tip_r, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_r + se_d_tip_rhythmic_none_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_r - se_d_tip_rhythmic_none_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_r(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_r + se_d_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_r - se_d_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_r(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_r + se_d_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_r - se_d_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
ylabel('Displacement (mm)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r + se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r - se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r + se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r - se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r + se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r - se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r + se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r - se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick dmax (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | R licks')

subplot(2, 4, 5)
hold on
yyaxis left
plot(inds_span_corr, mean_v_tip_all, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_all + se_v_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_all - se_v_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_all(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_all + se_v_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_all - se_v_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_all(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_all + se_v_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_all - se_v_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
ylabel('Velocity (mm/s)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all + se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all - se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all + se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all - se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all + se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all - se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all + seSS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all - se_SS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick dmax (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | All licks')

subplot(2, 4, 6)
hold on
yyaxis left
plot(inds_span_corr, mean_v_tip_l, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_l + se_v_tip_rhythmic_none_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_l - se_v_tip_rhythmic_none_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_l(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_l + se_v_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_l - se_v_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_l(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_l + se_v_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_l - se_v_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
ylabel('Velocity (mm/s)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l + se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l - se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l + se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l - se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l + se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l - se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l + se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l - se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick dmax (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | L licks')

subplot(2, 4, 7)
hold on
yyaxis left
plot(inds_span_corr, mean_v_tip_g, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_g + se_v_tip_rhythmic_none_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_g - se_v_tip_rhythmic_none_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_g(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_g + se_v_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_g - se_v_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_g(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_g + se_v_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_g - se_v_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
ylabel('Velocity (mm/s)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g + se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g - se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l + se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l - se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g + se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g - se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g + se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g - se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick dmax (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | G licks')

subplot(2, 4, 8)
hold on
yyaxis left
plot(inds_span_corr, mean_v_tip_r, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_r + se_v_tip_rhythmic_none_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_r - se_v_tip_rhythmic_none_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_r(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_r + se_v_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_r - se_v_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_r(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_r + se_v_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_r - se_v_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
ylabel('Velocity (mm/s)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r + se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r - se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r + se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r - se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r + se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r - se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r + se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r - se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick dmax (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | R licks')

%%  ANALYSIS PLOTS BASED ON LABELS - VMAX
clearvars -except ALL_PCELL_COMPRESSED_DATA POPULATION

number_pCell =  length(ALL_PCELL_COMPRESSED_DATA);

for counter_pCell = 1 : number_pCell
    id(counter_pCell, 1) = string(ALL_PCELL_COMPRESSED_DATA(counter_pCell).id{1, 1});
    
    %%% d_tip %%%
    d_tip_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmax.VID_d_tip_100_all_corr);
    d_tip_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmax.VID_d_tip_100_grooming_corr);
    d_tip_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmax.VID_d_tip_100_r_corr);
    d_tip_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmax.VID_d_tip_100_l_corr);
    
    %%% v_tip %%%
    v_tip_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmax.VID_v_tip_100_all_corr);
    v_tip_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmax.VID_v_tip_100_grooming_corr);
    v_tip_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmax.VID_v_tip_100_r_corr);
    v_tip_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmax.VID_v_tip_100_l_corr);
    
    %%% SS_str_bout %%%
    SS_str_bout_all(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmax.train_data_logic_SS_str_bout * 100), 15, 'sgolay', 2);
    SS_str_bout_baseline_all(counter_pCell, :) = nanmean(SS_str_bout_all(counter_pCell, 1:50));
    SS_str_bout_r(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmax.train_data_logic_SS_str_bout_r * 100), 15, 'sgolay', 2);
    SS_str_bout_baseline_r(counter_pCell, :) = nanmean(SS_str_bout_r(counter_pCell, 1:50));
    SS_str_bout_l(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmax.train_data_logic_SS_str_bout_l * 100), 15, 'sgolay', 2);
    SS_str_bout_baseline_l(counter_pCell, :) = nanmean(SS_str_bout_l(counter_pCell, 1:50));
    
    %%% SS_lick %%%
    SS_lick_all(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmax.train_data_logic_SS_all_corr * 100), 15, 'sgolay', 2);
    change_SS_lick_all(counter_pCell,:) = SS_lick_all(counter_pCell,:) - SS_str_bout_baseline_all(counter_pCell, 1);
    mean_window_SS_lick_all(counter_pCell,1) = nanmean(SS_lick_all(counter_pCell, :));
    mean_change_window_SS_lick_all(counter_pCell,1) = nanmean(change_SS_lick_all(counter_pCell, :));
    SS_lick_all(counter_pCell,:) = SS_lick_all(counter_pCell,:)/mean_window_SS_lick_all(counter_pCell,1);
    
    SS_lick_r(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmax.train_data_logic_SS_r_corr * 100), 15, 'sgolay', 2);
    change_SS_lick_r(counter_pCell,:) = SS_lick_r(counter_pCell,:) - SS_str_bout_baseline_r(counter_pCell, 1);
    mean_window_SS_lick_r(counter_pCell,1) = nanmean(SS_lick_r(counter_pCell, :));
    mean_change_window_SS_lick_r(counter_pCell,1) = nanmean(change_SS_lick_r(counter_pCell, :));
    SS_lick_r(counter_pCell,:) = SS_lick_r(counter_pCell,:)/mean_window_SS_lick_r(counter_pCell,1);
    
    SS_lick_g(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmax.train_data_logic_SS_grooming_corr * 100), 15, 'sgolay', 2);
    change_SS_lick_g(counter_pCell,:) = SS_lick_g(counter_pCell,:) - SS_str_bout_baseline_all(counter_pCell, 1);
    mean_window_SS_lick_g(counter_pCell,1) = nanmean(SS_lick_g(counter_pCell, :));
    mean_change_window_SS_lick_g(counter_pCell,1) = nanmean(change_SS_lick_g(counter_pCell, :));
    SS_lick_g(counter_pCell,:) = SS_lick_g(counter_pCell,:)/mean_window_SS_lick_g(counter_pCell,1);
    
    SS_lick_l(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmax.train_data_logic_SS_l_corr * 100), 15, 'sgolay', 2);
    change_SS_lick_l(counter_pCell,:) = SS_lick_l(counter_pCell,:) - SS_str_bout_baseline_l(counter_pCell, 1);
    mean_window_SS_lick_l(counter_pCell,1) = nanmean(SS_lick_l(counter_pCell, :));
    mean_change_window_SS_lick_l(counter_pCell,1) = nanmean(change_SS_lick_l(counter_pCell, :));
    SS_lick_l(counter_pCell,:) = SS_lick_l(counter_pCell,:)/mean_window_SS_lick_l(counter_pCell,1);
    
end

%%% LOAD LABELS FROM POPULATION STRUCTURE %%%
label_large_increase_all = POPULATION.dmax.label_large_increase_all;
label_medium_increase_all  = POPULATION.dmax.label_medium_increase_all;
label_small_increase_all = POPULATION.dmax.label_small_increase_all;
label_phase_all = POPULATION.dmax.label_phase_all;
label_small_decrease_all = POPULATION.dmax.label_small_decrease_all;
label_medium_decrease_all = POPULATION.dmax.label_medium_decrease_all;
label_large_decrease_all = POPULATION.dmax.label_large_decrease_all;
label_large_r = POPULATION.dmax.label_large_r;
label_medium_r = POPULATION.dmax.label_medium_r;
label_small_r = POPULATION.dmax.label_small_r;
label_none_dir = POPULATION.dmax.label_none_dir;
label_small_l = POPULATION.dmax.label_small_l;
label_medium_l = POPULATION.dmax.label_medium_l;
label_large_l = POPULATION.dmax.label_large_l;

label_rhythmic_all = POPULATION.dmax.label_rhythmic_all;
label_weakrhythmic_all = POPULATION.dmax.label_weakrhythmic_all;
label_rhythmic_l = POPULATION.dmax.label_rhythmic_l;
label_weakrhythmic_l = POPULATION.dmax.label_weakrhythmic_l;
label_rhythmic_g = POPULATION.dmax.label_rhythmic_g;
label_weakrhythmic_g = POPULATION.dmax.label_weakrhythmic_g;
label_rhythmic_r = POPULATION.dmax.label_rhythmic_r;
label_weakrhythmic_r = POPULATION.dmax.label_weakrhythmic_r;

label_lead_all = POPULATION.dmax.label_lead_all;
label_lead30to90_all  = POPULATION.dmax.label_lead30to90_all;
label_lead90to150_all  = POPULATION.dmax.label_lead90to150_all;
label_lag_all = POPULATION.dmax.label_lag_all;
label_lag30to90_all = POPULATION.dmax.label_lag30to90_all;
label_lag90to150_all = POPULATION.dmax.label_lag90to150_all;
label_phase_all = POPULATION.dmax.label_phase_all;
label_antiphase_all = POPULATION.dmax.label_antiphase_all;

label_lead_l = POPULATION.dmax.label_lead_l;
label_lead30to90_l  = POPULATION.dmax.label_lead30to90_l;
label_lead90to150_l  = POPULATION.dmax.label_lead90to150_l;
label_lag_l = POPULATION.dmax.label_lag_l;
label_lag30to90_l = POPULATION.dmax.label_lag30to90_l;
label_lag90to150_l = POPULATION.dmax.label_lag90to150_l;
label_phase_l = POPULATION.dmax.label_phase_l;
label_antiphase_l = POPULATION.dmax.label_antiphase_l;

label_lead_g = POPULATION.dmax.label_lead_g;
label_lead30to90_g  = POPULATION.dmax.label_lead30to90_g;
label_lead90to150_g  = POPULATION.dmax.label_lead90to150_g;
label_lag_g = POPULATION.dmax.label_lag_g;
label_lag30to90_g = POPULATION.dmax.label_lag30to90_g;
label_lag90to150_g = POPULATION.dmax.label_lag90to150_g;
label_phase_g = POPULATION.dmax.label_phase_g;
label_antiphase_g = POPULATION.dmax.label_antiphase_g;

label_lead_r = POPULATION.dmax.label_lead_r;
label_lead30to90_r  = POPULATION.dmax.label_lead30to90_r;
label_lead90to150_r  = POPULATION.dmax.label_lead90to150_r;
label_lag_r = POPULATION.dmax.label_lag_r;
label_lag30to90_r = POPULATION.dmax.label_lag30to90_r;
label_lag90to150_r = POPULATION.dmax.label_lag90to150_r;
label_phase_r = POPULATION.dmax.label_phase_r;
label_antiphase_r = POPULATION.dmax.label_antiphase_r;

%%% Mean d_tip based on labels
mean_d_tip_all = nanmean(d_tip_all);
mean_d_tip_l = nanmean(d_tip_l);
mean_d_tip_g = nanmean(d_tip_g);
mean_d_tip_r = nanmean(d_tip_r);

mean_d_tip_rhythmic_phase_all = nanmean(d_tip_all(label_rhythmic_all & label_phase_all,:));
se_d_tip_rhythmic_phase_all = nanstd(d_tip_all(label_rhythmic_all & label_phase_all,:))/sqrt(sum(label_rhythmic_all & label_phase_all));

mean_d_tip_rhythmic_antiphase_all = nanmean(d_tip_all(label_rhythmic_all & label_antiphase_all,:));
se_d_tip_rhythmic_antiphase_all = nanstd(d_tip_all(label_rhythmic_all & label_antiphase_all,:))/sqrt(sum(label_rhythmic_all & label_antiphase_all));

mean_d_tip_rhythmic_lead_all = nanmean(d_tip_all(label_rhythmic_all & label_lead_all,:));
se_d_tip_rhythmic_lead_all = nanstd(d_tip_all(label_rhythmic_all & label_lead_all,:))/sqrt(sum(label_rhythmic_all & label_lead_all));
mean_d_tip_rhythmic_lead30to90_all = nanmean(d_tip_all(label_rhythmic_all & label_lead30to90_all,:));
se_d_tip_rhythmic_lead30to90_all = nanstd(d_tip_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_d_tip_rhythmic_lead90to150_all = nanmean(d_tip_all(label_rhythmic_all & label_lead90to150_all,:));
se_d_tip_rhythmic_lead90to150_all = nanstd(d_tip_all(label_rhythmic_all & label_lead90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lead90to150_all));

mean_d_tip_rhythmic_lag_all = nanmean(d_tip_all(label_rhythmic_all & label_lag_all,:));
se_d_tip_rhythmic_lag_all = nanstd(d_tip_all(label_rhythmic_all & label_lag_all,:))/sqrt(sum(label_rhythmic_all & label_lag_all));
mean_d_tip_rhythmic_lead30to90_all = nanmean(d_tip_all(label_rhythmic_all & label_lead30to90_all,:));
se_d_tip_rhythmic_lead30to90_all = nanstd(d_tip_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_d_tip_rhythmic_lag90to150_all = nanmean(d_tip_all(label_rhythmic_all & label_lag90to150_all,:));
se_d_tip_rhythmic_lag90to150_all = nanstd(d_tip_all(label_rhythmic_all & label_lag90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lag90to150_all));

mean_d_tip_rhythmic_phase_l = nanmean(d_tip_l(label_rhythmic_l & label_phase_l,:));
se_d_tip_rhythmic_phase_l = nanstd(d_tip_l(label_rhythmic_l & label_phase_l,:))/sqrt(sum(label_rhythmic_l & label_phase_l));

mean_d_tip_rhythmic_antiphase_l = nanmean(d_tip_l(label_rhythmic_l & label_antiphase_l,:));
se_d_tip_rhythmic_antiphase_l = nanstd(d_tip_l(label_rhythmic_l & label_antiphase_l,:))/sqrt(sum(label_rhythmic_l & label_antiphase_l));

mean_d_tip_rhythmic_lead_l = nanmean(d_tip_l(label_rhythmic_l & label_lead_l,:));
se_d_tip_rhythmic_lead_l = nanstd(d_tip_l(label_rhythmic_l & label_lead_l,:))/sqrt(sum(label_rhythmic_l & label_lead_l));
mean_d_tip_rhythmic_lead30to90_l = nanmean(d_tip_l(label_rhythmic_l & label_lead30to90_l,:));
se_d_tip_rhythmic_lead30to90_l = nanstd(d_tip_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_d_tip_rhythmic_lead90to150_l = nanmean(d_tip_l(label_rhythmic_l & label_lead90to150_l,:));
se_d_tip_rhythmic_lead90to150_l = nanstd(d_tip_l(label_rhythmic_l & label_lead90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lead90to150_l));

mean_d_tip_rhythmic_lag_l = nanmean(d_tip_l(label_rhythmic_l & label_lag_l,:));
se_d_tip_rhythmic_lag_l = nanstd(d_tip_l(label_rhythmic_l & label_lag_l,:))/sqrt(sum(label_rhythmic_l & label_lag_l));
mean_d_tip_rhythmic_lead30to90_l = nanmean(d_tip_l(label_rhythmic_l & label_lead30to90_l,:));
se_d_tip_rhythmic_lead30to90_l = nanstd(d_tip_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_d_tip_rhythmic_lag90to150_l = nanmean(d_tip_l(label_rhythmic_l & label_lag90to150_l,:));
se_d_tip_rhythmic_lag90to150_l = nanstd(d_tip_l(label_rhythmic_l & label_lag90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lag90to150_l));

mean_d_tip_rhythmic_phase_g = nanmean(d_tip_g(label_rhythmic_g & label_phase_g,:));
se_d_tip_rhythmic_phase_g = nanstd(d_tip_g(label_rhythmic_g & label_phase_g,:))/sqrt(sum(label_rhythmic_g & label_phase_g));

mean_d_tip_rhythmic_antiphase_g = nanmean(d_tip_g(label_rhythmic_g & label_antiphase_g,:));
se_d_tip_rhythmic_antiphase_g = nanstd(d_tip_g(label_rhythmic_g & label_antiphase_g,:))/sqrt(sum(label_rhythmic_g & label_antiphase_g));

mean_d_tip_rhythmic_lead_g = nanmean(d_tip_g(label_rhythmic_g & label_lead_g,:));
se_d_tip_rhythmic_lead_g = nanstd(d_tip_g(label_rhythmic_g & label_lead_g,:))/sqrt(sum(label_rhythmic_g & label_lead_g));
mean_d_tip_rhythmic_lead30to90_g = nanmean(d_tip_g(label_rhythmic_g & label_lead30to90_g,:));
se_d_tip_rhythmic_lead30to90_g = nanstd(d_tip_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_d_tip_rhythmic_lead90to150_g = nanmean(d_tip_g(label_rhythmic_g & label_lead90to150_g,:));
se_d_tip_rhythmic_lead90to150_g = nanstd(d_tip_g(label_rhythmic_g & label_lead90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lead90to150_g));

mean_d_tip_rhythmic_lag_g = nanmean(d_tip_g(label_rhythmic_g & label_lag_g,:));
se_d_tip_rhythmic_lag_g = nanstd(d_tip_g(label_rhythmic_g & label_lag_g,:))/sqrt(sum(label_rhythmic_g & label_lag_g));
mean_d_tip_rhythmic_lead30to90_g = nanmean(d_tip_g(label_rhythmic_g & label_lead30to90_g,:));
se_d_tip_rhythmic_lead30to90_g = nanstd(d_tip_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_d_tip_rhythmic_lag90to150_g = nanmean(d_tip_g(label_rhythmic_g & label_lag90to150_g,:));
se_d_tip_rhythmic_lag90to150_g = nanstd(d_tip_g(label_rhythmic_g & label_lag90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lag90to150_g));

mean_d_tip_rhythmic_phase_r = nanmean(d_tip_r(label_rhythmic_r & label_phase_r,:));
se_d_tip_rhythmic_phase_r = nanstd(d_tip_r(label_rhythmic_r & label_phase_r,:))/sqrt(sum(label_rhythmic_r & label_phase_r));

mean_d_tip_rhythmic_antiphase_r = nanmean(d_tip_g(label_rhythmic_r & label_antiphase_r,:));
se_d_tip_rhythmic_antiphase_r = nanstd(d_tip_g(label_rhythmic_r & label_antiphase_r,:))/sqrt(sum(label_rhythmic_r & label_antiphase_r));

mean_d_tip_rhythmic_lead_r = nanmean(d_tip_r(label_rhythmic_r & label_lead_r,:));
se_d_tip_rhythmic_lead_r = nanstd(d_tip_r(label_rhythmic_r & label_lead_r,:))/sqrt(sum(label_rhythmic_r & label_lead_r));
mean_d_tip_rhythmic_lead30to90_r = nanmean(d_tip_r(label_rhythmic_r & label_lead30to90_r,:));
se_d_tip_rhythmic_lead30to90_r = nanstd(d_tip_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_d_tip_rhythmic_lead90to150_r = nanmean(d_tip_r(label_rhythmic_r & label_lead90to150_r,:));
se_d_tip_rhythmic_lead90to150_r = nanstd(d_tip_r(label_rhythmic_r & label_lead90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lead90to150_r));

mean_d_tip_rhythmic_lag_r = nanmean(d_tip_r(label_rhythmic_r & label_lag_r,:));
se_d_tip_rhythmic_lag_r = nanstd(d_tip_r(label_rhythmic_r & label_lag_r,:))/sqrt(sum(label_rhythmic_r & label_lag_r));
mean_d_tip_rhythmic_lead30to90_r = nanmean(d_tip_r(label_rhythmic_r & label_lead30to90_r,:));
se_d_tip_rhythmic_lead30to90_r = nanstd(d_tip_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_d_tip_rhythmic_lag90to150_r = nanmean(d_tip_r(label_rhythmic_r & label_lag90to150_r,:));
se_d_tip_rhythmic_lag90to150_r = nanstd(d_tip_r(label_rhythmic_r & label_lag90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lag90to150_r));

%%% Mean v_tip based on labels
mean_v_tip_all = nanmean(v_tip_all);
mean_v_tip_l = nanmean(v_tip_l);
mean_v_tip_g = nanmean(v_tip_g);
mean_v_tip_r = nanmean(v_tip_r);

mean_v_tip_rhythmic_phase_all = nanmean(v_tip_all(label_rhythmic_all & label_phase_all,:));
se_v_tip_rhythmic_phase_all = nanstd(v_tip_all(label_rhythmic_all & label_phase_all,:))/sqrt(sum(label_rhythmic_all & label_phase_all));

mean_v_tip_rhythmic_antiphase_all = nanmean(v_tip_all(label_rhythmic_all & label_antiphase_all,:));
se_v_tip_rhythmic_antiphase_all = nanstd(v_tip_all(label_rhythmic_all & label_antiphase_all,:))/sqrt(sum(label_rhythmic_all & label_antiphase_all));

mean_v_tip_rhythmic_lead_all = nanmean(v_tip_all(label_rhythmic_all & label_lead_all,:));
se_v_tip_rhythmic_lead_all = nanstd(v_tip_all(label_rhythmic_all & label_lead_all,:))/sqrt(sum(label_rhythmic_all & label_lead_all));
mean_v_tip_rhythmic_lead30to90_all = nanmean(v_tip_all(label_rhythmic_all & label_lead30to90_all,:));
se_v_tip_rhythmic_lead30to90_all = nanstd(v_tip_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_v_tip_rhythmic_lead90to150_all = nanmean(v_tip_all(label_rhythmic_all & label_lead90to150_all,:));
se_v_tip_rhythmic_lead90to150_all = nanstd(v_tip_all(label_rhythmic_all & label_lead90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lead90to150_all));

mean_v_tip_rhythmic_lag_all = nanmean(v_tip_all(label_rhythmic_all & label_lag_all,:));
se_v_tip_rhythmic_lag_all = nanstd(v_tip_all(label_rhythmic_all & label_lag_all,:))/sqrt(sum(label_rhythmic_all & label_lag_all));
mean_v_tip_rhythmic_lead30to90_all = nanmean(v_tip_all(label_rhythmic_all & label_lead30to90_all,:));
se_v_tip_rhythmic_lead30to90_all = nanstd(v_tip_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_v_tip_rhythmic_lag90to150_all = nanmean(v_tip_all(label_rhythmic_all & label_lag90to150_all,:));
se_v_tip_rhythmic_lag90to150_all = nanstd(v_tip_all(label_rhythmic_all & label_lag90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lag90to150_all));

mean_v_tip_rhythmic_phase_l = nanmean(v_tip_l(label_rhythmic_l & label_phase_l,:));
se_v_tip_rhythmic_phasell = nanstd(v_tip_l(label_rhythmic_l & label_phase_l,:))/sqrt(sum(label_rhythmic_l & label_phase_l));

mean_v_tip_rhythmic_antiphase_l = nanmean(v_tip_l(label_rhythmic_l & label_antiphase_l,:));
se_v_tip_rhythmic_antiphase_l = nanstd(v_tip_l(label_rhythmic_l & label_antiphase_l,:))/sqrt(sum(label_rhythmic_l & label_antiphase_l));

mean_v_tip_rhythmic_lead_l = nanmean(v_tip_l(label_rhythmic_l & label_lead_l,:));
se_v_tip_rhythmic_lead_l = nanstd(v_tip_l(label_rhythmic_l & label_lead_l,:))/sqrt(sum(label_rhythmic_l & label_lead_l));
mean_v_tip_rhythmic_lead30to90_l = nanmean(v_tip_l(label_rhythmic_l & label_lead30to90_l,:));
se_v_tip_rhythmic_lead30to90_l = nanstd(v_tip_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_v_tip_rhythmic_lead90to150_l = nanmean(v_tip_l(label_rhythmic_l & label_lead90to150_l,:));
se_v_tip_rhythmic_lead90to150_l = nanstd(v_tip_l(label_rhythmic_l & label_lead90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lead90to150_l));

mean_v_tip_rhythmic_lag_l = nanmean(v_tip_l(label_rhythmic_l & label_lag_l,:));
se_v_tip_rhythmic_lag_l = nanstd(v_tip_l(label_rhythmic_l & label_lag_l,:))/sqrt(sum(label_rhythmic_l & label_lag_l));
mean_v_tip_rhythmic_lead30to90_l = nanmean(v_tip_l(label_rhythmic_l & label_lead30to90_l,:));
se_v_tip_rhythmic_lead30to90_l = nanstd(v_tip_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_v_tip_rhythmic_lag90to150_l = nanmean(v_tip_l(label_rhythmic_l & label_lag90to150_l,:));
se_v_tip_rhythmic_lag90to150_l = nanstd(v_tip_l(label_rhythmic_l & label_lag90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lag90to150_l));

mean_v_tip_rhythmic_phase_g = nanmean(v_tip_g(label_rhythmic_g & label_phase_g,:));
se_v_tip_rhythmic_phase_g = nanstd(v_tip_g(label_rhythmic_g & label_phase_g,:))/sqrt(sum(label_rhythmic_g & label_phase_g));

mean_v_tip_rhythmic_antiphase_g = nanmean(v_tip_g(label_rhythmic_g & label_antiphase_g,:));
se_v_tip_rhythmic_antiphase_g = nanstd(v_tip_g(label_rhythmic_g & label_antiphase_g,:))/sqrt(sum(label_rhythmic_g & label_antiphase_g));

mean_v_tip_rhythmic_lead_g = nanmean(v_tip_g(label_rhythmic_g & label_lead_g,:));
se_v_tip_rhythmic_lead_g = nanstd(v_tip_g(label_rhythmic_g & label_lead_g,:))/sqrt(sum(label_rhythmic_g & label_lead_g));
mean_v_tip_rhythmic_lead30to90_g = nanmean(v_tip_g(label_rhythmic_g & label_lead30to90_g,:));
se_v_tip_rhythmic_lead30to90_g = nanstd(v_tip_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_v_tip_rhythmic_lead90to150_g = nanmean(v_tip_g(label_rhythmic_g & label_lead90to150_g,:));
se_v_tip_rhythmic_lead90to150_g = nanstd(v_tip_g(label_rhythmic_g & label_lead90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lead90to150_g));

mean_v_tip_rhythmic_lag_g = nanmean(v_tip_g(label_rhythmic_g & label_lag_g,:));
se_v_tip_rhythmic_lag_g = nanstd(v_tip_g(label_rhythmic_g & label_lag_g,:))/sqrt(sum(label_rhythmic_g & label_lag_g));
mean_v_tip_rhythmic_lead30to90_g = nanmean(v_tip_g(label_rhythmic_g & label_lead30to90_g,:));
se_v_tip_rhythmic_lead30to90_g = nanstd(v_tip_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_v_tip_rhythmic_lag90to150_g = nanmean(v_tip_g(label_rhythmic_g & label_lag90to150_g,:));
se_v_tip_rhythmic_lag90to150_g = nanstd(v_tip_g(label_rhythmic_g & label_lag90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lag90to150_g));

mean_v_tip_rhythmic_phase_r = nanmean(v_tip_r(label_rhythmic_r & label_phase_r,:));
se_v_tip_rhythmic_phase_r = nanstd(v_tip_r(label_rhythmic_r & label_phase_r,:))/sqrt(sum(label_rhythmic_r & label_phase_r));

mean_v_tip_rhythmic_antiphase_r = nanmean(v_tip_r(label_rhythmic_r & label_antiphase_r,:));
se_v_tip_rhythmic_antiphaser = nanstd(v_tip_r(label_rhythmic_r & label_antiphase_r,:))/sqrt(sum(label_rhythmic_r & label_antiphase_r));

mean_v_tip_rhythmic_lead_r = nanmean(v_tip_r(label_rhythmic_r & label_lead_r,:));
se_v_tip_rhythmic_lead_r = nanstd(v_tip_r(label_rhythmic_r & label_lead_r,:))/sqrt(sum(label_rhythmic_r & label_lead_r));
mean_v_tip_rhythmic_lead30to90_r = nanmean(v_tip_r(label_rhythmic_r & label_lead30to90_r,:));
se_v_tip_rhythmic_lead30to90_r = nanstd(v_tip_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_v_tip_rhythmic_lead90to150_r = nanmean(v_tip_r(label_rhythmic_r & label_lead90to150_r,:));
se_v_tip_rhythmic_lead90to150_r = nanstd(v_tip_r(label_rhythmic_r & label_lead90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lead90to150_r));

mean_v_tip_rhythmic_lag_r = nanmean(v_tip_r(label_rhythmic_r & label_lag_r,:));
se_v_tip_rhythmic_lag_r = nanstd(v_tip_r(label_rhythmic_r & label_lag_r,:))/sqrt(sum(label_rhythmic_r & label_lag_r));
mean_v_tip_rhythmic_lead30to90_r = nanmean(v_tip_r(label_rhythmic_r & label_lead30to90_r,:));
se_v_tip_rhythmic_lead30to90_r = nanstd(v_tip_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_v_tip_rhythmic_lag90to150_r = nanmean(v_tip_r(label_rhythmic_r & label_lag90to150_r,:));
se_v_tip_rhythmic_lag90to150_r = nanstd(v_tip_r(label_rhythmic_r & label_lag90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lag90to150_r));

%%% Mean change_SS_lick based on labels
mean_SS_lick_rhythmic_phase_all = nanmean(SS_lick_all(label_rhythmic_all & label_phase_all,:));
se_SS_lick_rhythmic_phase_all = nanstd(SS_lick_all(label_rhythmic_all & label_phase_all,:))/sqrt(sum(label_rhythmic_all & label_phase_all));

mean_SS_lick_rhythmic_antiphase_all = nanmean(SS_lick_all(label_rhythmic_all & label_antiphase_all,:));
se_SS_lick_rhythmic_antiphase_all = nanstd(SS_lick_all(label_rhythmic_all & label_antiphase_all,:))/sqrt(sum(label_rhythmic_all & label_antiphase_all));

mean_SS_lick_rhythmic_lead_all = nanmean(SS_lick_all(label_rhythmic_all & label_lead_all,:));
se_SS_lick_rhythmic_lead_all = nanstd(SS_lick_all(label_rhythmic_all & label_lead_all,:))/sqrt(sum(label_rhythmic_all & label_lead_all));
mean_SS_lick_rhythmic_lead30to90_all = nanmean(SS_lick_all(label_rhythmic_all & label_lead30to90_all,:));
se_SS_lick_rhythmic_lead30to90_all = nanstd(SS_lick_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_SS_lick_rhythmic_lead90to150_all = nanmean(SS_lick_all(label_rhythmic_all & label_lead90to150_all,:));
se_SS_lick_rhythmic_lead90to150_all = nanstd(SS_lick_all(label_rhythmic_all & label_lead90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lead90to150_all));

mean_SS_lick_rhythmic_lag_all = nanmean(SS_lick_all(label_rhythmic_all & label_lag_all,:));
se_SS_lick_rhythmic_lag_all = nanstd(SS_lick_all(label_rhythmic_all & label_lag_all,:))/sqrt(sum(label_rhythmic_all & label_lag_all));
mean_SS_lick_rhythmic_lead30to90_all = nanmean(SS_lick_all(label_rhythmic_all & label_lead30to90_all,:));
se_SS_lick_rhythmic_lead30to90_all = nanstd(SS_lick_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_SS_lick_rhythmic_lag90to150_all = nanmean(SS_lick_all(label_rhythmic_all & label_lag90to150_all,:));
se_SS_lick_rhythmic_lag90to150_all = nanstd(SS_lick_all(label_rhythmic_all & label_lag90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lag90to150_all));

mean_SS_lick_rhythmic_phase_l = nanmean(SS_lick_l(label_rhythmic_l & label_phase_l,:));
se_SS_lick_rhythmic_phase_l = nanstd(SS_lick_l(label_rhythmic_l & label_phase_l,:))/sqrt(sum(label_rhythmic_l & label_phase_l));

mean_SS_lick_rhythmic_antiphase_l = nanmean(SS_lick_l(label_rhythmic_l & label_antiphase_l,:));
se_SS_lick_rhythmic_antiphase_l = nanstd(SS_lick_l(label_rhythmic_l & label_antiphase_l,:))/sqrt(sum(label_rhythmic_l & label_antiphase_l));

mean_SS_lick_rhythmic_lead_l = nanmean(SS_lick_l(label_rhythmic_l & label_lead_l,:));
se_SS_lick_rhythmic_lead_l = nanstd(SS_lick_l(label_rhythmic_l & label_lead_l,:))/sqrt(sum(label_rhythmic_l & label_lead_l));
mean_SS_lick_rhythmic_lead30to90_l = nanmean(SS_lick_l(label_rhythmic_l & label_lead30to90_l,:));
se_SS_lick_rhythmic_lead30to90_l = nanstd(SS_lick_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_SS_lick_rhythmic_lead90to150_l = nanmean(SS_lick_l(label_rhythmic_l & label_lead90to150_l,:));
se_SS_lick_rhythmic_lead90to150_l = nanstd(SS_lick_l(label_rhythmic_l & label_lead90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lead90to150_l));

mean_SS_lick_rhythmic_lag_l = nanmean(SS_lick_l(label_rhythmic_l & label_lag_l,:));
se_SS_lick_rhythmic_lag_l = nanstd(SS_lick_l(label_rhythmic_l & label_lag_l,:))/sqrt(sum(label_rhythmic_l & label_lag_l));
mean_SS_lick_rhythmic_lead30to90_l = nanmean(SS_lick_l(label_rhythmic_l & label_lead30to90_l,:));
se_SS_lick_rhythmic_lead30to90_l = nanstd(SS_lick_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_SS_lick_rhythmic_lag90to150_l = nanmean(SS_lick_l(label_rhythmic_l & label_lag90to150_l,:));
se_SS_lick_rhythmic_lag90to150_l = nanstd(SS_lick_l(label_rhythmic_l & label_lag90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lag90to150_l));

mean_SS_lick_rhythmic_phase_g = nanmean(SS_lick_g(label_rhythmic_g & label_phase_g,:));
se_SS_lick_rhythmic_phase_g = nanstd(SS_lick_g(label_rhythmic_g & label_phase_g,:))/sqrt(sum(label_rhythmic_g & label_phase_g));

mean_SS_lick_rhythmic_antiphase_g = nanmean(SS_lick_g(label_rhythmic_g & label_antiphase_g,:));
se_SS_lick_rhythmic_antiphase_g = nanstd(SS_lick_g(label_rhythmic_g & label_antiphase_g,:))/sqrt(sum(label_rhythmic_g & label_antiphase_g));

mean_SS_lick_rhythmic_lead_g = nanmean(SS_lick_g(label_rhythmic_g & label_lead_g,:));
se_SS_lick_rhythmic_lead_g = nanstd(SS_lick_g(label_rhythmic_g & label_lead_g,:))/sqrt(sum(label_rhythmic_g & label_lead_g));
mean_SS_lick_rhythmic_lead30to90_g = nanmean(SS_lick_g(label_rhythmic_g & label_lead30to90_g,:));
se_SS_lick_rhythmic_lead30to90_g = nanstd(SS_lick_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_SS_lick_rhythmic_lead90to150_g = nanmean(SS_lick_g(label_rhythmic_g & label_lead90to150_g,:));
se_SS_lick_rhythmic_lead90to150_g = nanstd(SS_lick_g(label_rhythmic_g & label_lead90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lead90to150_g));

mean_SS_lick_rhythmic_lag_g = nanmean(SS_lick_g(label_rhythmic_g & label_lag_g,:));
se_SS_lick_rhythmic_lag_g = nanstd(SS_lick_g(label_rhythmic_g & label_lag_g,:))/sqrt(sum(label_rhythmic_g & label_lag_g));
mean_SS_lick_rhythmic_lead30to90_g = nanmean(SS_lick_g(label_rhythmic_g & label_lead30to90_g,:));
se_SS_lick_rhythmic_lead30to90_g = nanstd(SS_lick_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_SS_lick_rhythmic_lag90to150_g = nanmean(SS_lick_g(label_rhythmic_g & label_lag90to150_g,:));
se_SS_lick_rhythmic_lag90to150_g = nanstd(SS_lick_g(label_rhythmic_g & label_lag90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lag90to150_g));

mean_SS_lick_rhythmic_phase_r = nanmean(SS_lick_r(label_rhythmic_r & label_phase_r,:));
se_SS_lick_rhythmic_phase_r = nanstd(SS_lick_r(label_rhythmic_r & label_phase_r,:))/sqrt(sum(label_rhythmic_r & label_phase_r));

mean_SS_lick_rhythmic_antiphase_r = nanmean(SS_lick_r(label_rhythmic_r & label_antiphase_r,:));
se_SS_lick_rhythmic_antiphase_r = nanstd(SS_lick_r(label_rhythmic_r & label_antiphase_r,:))/sqrt(sum(label_rhythmic_r & label_antiphase_r));

mean_SS_lick_rhythmic_lead_r = nanmean(SS_lick_r(label_rhythmic_r & label_lead_r,:));
se_SS_lick_rhythmic_lead_r = nanstd(SS_lick_r(label_rhythmic_r & label_lead_r,:))/sqrt(sum(label_rhythmic_r & label_lead_r));
mean_SS_lick_rhythmic_lead30to90_r = nanmean(SS_lick_r(label_rhythmic_r & label_lead30to90_r,:));
se_SS_lick_rhythmic_lead30to90_r = nanstd(SS_lick_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_SS_lick_rhythmic_lead90to150_r = nanmean(SS_lick_r(label_rhythmic_r & label_lead90to150_r,:));
se_SS_lick_rhythmic_lead90to150_r = nanstd(SS_lick_r(label_rhythmic_r & label_lead90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lead90to150_r));

mean_SS_lick_rhythmic_lag_r = nanmean(SS_lick_r(label_rhythmic_r & label_lag_r,:));
se_SS_lick_rhythmic_lag_r = nanstd(SS_lick_r(label_rhythmic_r & label_lag_r,:))/sqrt(sum(label_rhythmic_r & label_lag_r));
mean_SS_lick_rhythmic_lead30to90_r = nanmean(SS_lick_r(label_rhythmic_r & label_lead30to90_r,:));
se_SS_lick_rhythmic_lead30to90_r = nanstd(SS_lick_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_SS_lick_rhythmic_lag90to150_r = nanmean(SS_lick_r(label_rhythmic_r & label_lag90to150_r,:));
se_SS_lick_rhythmic_lag90to150_r = nanstd(SS_lick_r(label_rhythmic_r & label_lag90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lag90to150_r));

%%% inds_span_corr %%%
inds_span_corr = -0.99 : 0.01 : 1;
%% PLOT LABELED TRACES
figure
subplot(3, 4, 1)
hold on
plot(inds_span_corr, mean_d_tip_rhythmic_phase_all(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_all + se_d_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_all - se_d_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_all(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_all + se_d_tip_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_all - se_d_tip_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lead_all(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_all + se_d_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_all - se_d_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lag_all(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_all + se_d_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_all - se_d_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-0.17 0.17])
ylabel('Displacement (mm)')
xlabel('Lick vmax (ms)')
title('Lick trace | All licks')

subplot(3, 4, 2)
hold on
plot(inds_span_corr, mean_d_tip_rhythmic_phase_l(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_l + se_d_tip_rhythmic_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_l - se_d_tip_rhythmic_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_l(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_l + se_d_tip_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_l - se_d_tip_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lead_l(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_l + se_d_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_l - se_d_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lag_l(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_l + se_d_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_l - se_d_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-0.17 0.17])
ylabel('Displacement (mm)')
xlabel('Lick vmax (ms)')
title('Lick trace | L licks')

subplot(3, 4, 3)
hold on
plot(inds_span_corr, mean_d_tip_rhythmic_phase_g(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_g + se_d_tip_rhythmic_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_g - se_d_tip_rhythmic_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_g(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_g + se_d_tip_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_g - se_d_tip_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lead_g(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_g + se_d_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_g - se_d_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lag_g(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_g + se_d_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_g - se_d_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-0.17 0.17])
ylabel('Displacement (mm)')
xlabel('Lick vmax (ms)')
title('Lick trace | G licks')

subplot(3, 4, 4)
hold on
plot(inds_span_corr, mean_d_tip_rhythmic_phase_r(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_r + se_d_tip_rhythmic_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_r - se_d_tip_rhythmic_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_r(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_r + se_d_tip_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_r - se_d_tip_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lead_r(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_r + se_d_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_r - se_d_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lag_r(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_r + se_d_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_r - se_d_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-0.17 0.17])
ylabel('Displacement (mm)')
xlabel('Lick vmax (ms)')
title('Lick trace | R licks')


subplot(3, 4, 5)
hold on
plot(inds_span_corr, mean_v_tip_rhythmic_phase_all(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_all + se_v_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_all - se_v_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_all(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_all + se_v_tip_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_all - se_v_tip_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lead_all(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_all +se_v_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_all - se_v_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lag_all(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_all + se_v_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_all - se_v_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
xlim([-0.17 0.17])
ylabel('Velocity (mm/s)')
xlabel('Lick vmax (ms)')
title('Lick trace | All licks')

subplot(3, 4, 6)
hold on
plot(inds_span_corr, mean_v_tip_rhythmic_phase_l(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_l + se_v_tip_rhythmic_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_l - se_v_tip_rhythmic_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_l(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_l + se_v_tip_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_l - se_v_tip_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lead_l(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_l + se_v_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_l - se_v_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lag_l(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_l + se_v_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_l - se_v_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
xlim([-0.17 0.17])
ylabel('Velocity (mm/s)')
xlabel('Lick vmax (ms)')
title('Lick trace | L licks')

subplot(3, 4, 7)
hold on
plot(inds_span_corr, mean_v_tip_rhythmic_phase_g(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_g + se_v_tip_rhythmic_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_g - se_v_tip_rhythmic_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_g(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_g + se_v_tip_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_g - se_v_tip_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lead_g(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_g + se_v_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_g - se_v_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lag_g(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_g + se_v_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_g - se_v_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
xlim([-0.17 0.17])
ylabel('Velocity (mm/s)')
xlabel('Lick vmax (ms)')
title('Lick trace | G licks')

subplot(3, 4,8)
hold on
plot(inds_span_corr, mean_v_tip_rhythmic_phase_r(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_r + se_v_tip_rhythmic_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_r - se_v_tip_rhythmic_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_r(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_r + se_v_tip_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_r - se_v_tip_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lead_r(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_r + se_v_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_r - se_v_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lag_r(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_r + se_v_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_r - se_v_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
xlim([-0.17 0.17])
ylabel('Velocity (mm/s)')
xlabel('Lick vmax (ms)')
title('Lick trace | R licks')

subplot(3, 4, 9)
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all + se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all - se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_all + se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_all - se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all + se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all - se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all + se_SS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all -se_SS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
%     ylabel('Change in SS firing (spks/s)')
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmax (ms)')
title('SS firing | All licks')

subplot(3, 4, 10)
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l + se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l - se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_l + se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_l - se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l + se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l - se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l + se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l - se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
%     ylabel('Change in SS firing (spks/s)')
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmax (ms)')
title('SS firing | L licks')

subplot(3, 4, 11)
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g + se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g - se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_g(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_g + se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_g - se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g + se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g - se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g + se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g - se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
%     ylabel('Change in SS firing (spks/s)')
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmax (ms)')
title('SS firing | G licks')

subplot(3, 4, 12)
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r + se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r - se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_r + se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_r - se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r + se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r - se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r + se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r - se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
%     ylabel('Change in SS firing (spks/s)')
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmax (ms)')
title('SS firing | R licks')
%% PLOT LABELED TRACES - Kinematics
figure
subplot(2, 4, 1)
hold on
yyaxis left
plot(inds_span_corr, mean_d_tip_all, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_all + se_d_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_all - se_d_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_all(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_all + se_d_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_all - se_d_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_all(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_all + se_d_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_all - se_d_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
ylabel('Displacement (mm)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all + se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all - se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all + se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all - se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all + se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all - se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all + se_SS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all - se_SS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmax (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | All licks')

subplot(2, 4, 2)
hold on
yyaxis left
plot(inds_span_corr, mean_d_tip_l, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_l + se_d_tip_rhythmic_none_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_l - se_d_tip_rhythmic_none_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_l(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_l + se_d_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_l - se_d_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_l(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_l + se_d_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_l - se_d_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
ylabel('Displacement (mm)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l + se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l - se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l + se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l - se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l + se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l - se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l + se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l - se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmax (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | L licks')


subplot(2, 4, 3)
hold on
yyaxis left
plot(inds_span_corr, mean_d_tip_g, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_g + se_d_tip_rhythmic_none_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_g - se_d_tip_rhythmic_none_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_g(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_g + se_d_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_g - se_d_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_g(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_g + se_d_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_g - se_d_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
ylabel('Displacement (mm)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g + se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g - se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_g(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_g + se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_g - se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g + se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g - se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g + se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g - se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmax (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | G licks')


subplot(2, 4, 4)
hold on
yyaxis left
plot(inds_span_corr, mean_d_tip_r, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_r + se_d_tip_rhythmic_none_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_r - se_d_tip_rhythmic_none_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_r(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_r + se_d_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_r - se_d_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_r(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_r + se_d_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_r - se_d_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
ylabel('Displacement (mm)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r + se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r - se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r + se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r - se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r + se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r - se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r + se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r - se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmax (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | R licks')

subplot(2, 4, 5)
hold on
yyaxis left
plot(inds_span_corr, mean_v_tip_all, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_all + se_v_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_all - se_v_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_all(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_all + se_v_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_all - se_v_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_all(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_all + se_v_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_all - se_v_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
ylabel('Velocity (mm/s)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all + se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all - se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all + se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all - se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all + se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all - se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all + seSS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all - se_SS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmax (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | All licks')

subplot(2, 4, 6)
hold on
yyaxis left
plot(inds_span_corr, mean_v_tip_l, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_l + se_v_tip_rhythmic_none_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_l - se_v_tip_rhythmic_none_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_l(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_l + se_v_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_l - se_v_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_l(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_l + se_v_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_l - se_v_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
ylabel('Velocity (mm/s)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l + se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l - se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l + se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l - se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l + se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l - se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l + se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l - se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmax (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | L licks')

subplot(2, 4, 7)
hold on
yyaxis left
plot(inds_span_corr, mean_v_tip_g, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_g + se_v_tip_rhythmic_none_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_g - se_v_tip_rhythmic_none_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_g(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_g + se_v_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_g - se_v_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_g(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_g + se_v_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_g - se_v_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
ylabel('Velocity (mm/s)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g + se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g - se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l + se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l - se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g + se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g - se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g + se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g - se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmax (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | G licks')

subplot(2, 4, 8)
hold on
yyaxis left
plot(inds_span_corr, mean_v_tip_r, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_r + se_v_tip_rhythmic_none_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_r - se_v_tip_rhythmic_none_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_r(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_r + se_v_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_r - se_v_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_r(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_r + se_v_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_r - se_v_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
ylabel('Velocity (mm/s)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r + se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r - se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r + se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r - se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r + se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r - se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r + se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r - se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmax (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | R licks')

%%  ANALYSIS PLOTS BASED ON LABELS - VMIN
clearvars -except ALL_PCELL_COMPRESSED_DATA POPULATION

number_pCell =  length(ALL_PCELL_COMPRESSED_DATA);

for counter_pCell = 1 : number_pCell
    id(counter_pCell, 1) = string(ALL_PCELL_COMPRESSED_DATA(counter_pCell).id{1, 1});
    
    %%% d_tip %%%
    d_tip_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmin.VID_d_tip_100_all_corr);
    d_tip_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmin.VID_d_tip_100_grooming_corr);
    d_tip_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmin.VID_d_tip_100_r_corr);
    d_tip_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmin.VID_d_tip_100_l_corr);
    
    %%% v_tip %%%
    v_tip_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmin.VID_v_tip_100_all_corr);
    v_tip_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmin.VID_v_tip_100_grooming_corr);
    v_tip_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmin.VID_v_tip_100_r_corr);
    v_tip_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmin.VID_v_tip_100_l_corr);
    
    %%% SS_str_bout %%%
    SS_str_bout_all(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmin.train_data_logic_SS_str_bout * 100), 15, 'sgolay', 2);
    SS_str_bout_baseline_all(counter_pCell, :) = nanmean(SS_str_bout_all(counter_pCell, 1:50));
    SS_str_bout_r(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmin.train_data_logic_SS_str_bout_r * 100), 15, 'sgolay', 2);
    SS_str_bout_baseline_r(counter_pCell, :) = nanmean(SS_str_bout_r(counter_pCell, 1:50));
    SS_str_bout_l(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmin.train_data_logic_SS_str_bout_l * 100), 15, 'sgolay', 2);
    SS_str_bout_baseline_l(counter_pCell, :) = nanmean(SS_str_bout_l(counter_pCell, 1:50));
    
    %%% SS_lick %%%
    SS_lick_all(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmin.train_data_logic_SS_all_corr * 100), 15, 'sgolay', 2);
    change_SS_lick_all(counter_pCell,:) = SS_lick_all(counter_pCell,:) - SS_str_bout_baseline_all(counter_pCell, 1);
    mean_window_SS_lick_all(counter_pCell,1) = nanmean(SS_lick_all(counter_pCell, :));
    mean_change_window_SS_lick_all(counter_pCell,1) = nanmean(change_SS_lick_all(counter_pCell, :));
    SS_lick_all(counter_pCell,:) = SS_lick_all(counter_pCell,:)/mean_window_SS_lick_all(counter_pCell,1);
    
    SS_lick_r(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmin.train_data_logic_SS_r_corr * 100), 15, 'sgolay', 2);
    change_SS_lick_r(counter_pCell,:) = SS_lick_r(counter_pCell,:) - SS_str_bout_baseline_r(counter_pCell, 1);
    mean_window_SS_lick_r(counter_pCell,1) = nanmean(SS_lick_r(counter_pCell, :));
    mean_change_window_SS_lick_r(counter_pCell,1) = nanmean(change_SS_lick_r(counter_pCell, :));
    SS_lick_r(counter_pCell,:) = SS_lick_r(counter_pCell,:)/mean_window_SS_lick_r(counter_pCell,1);
    
    SS_lick_g(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmin.train_data_logic_SS_grooming_corr * 100), 15, 'sgolay', 2);
    change_SS_lick_g(counter_pCell,:) = SS_lick_g(counter_pCell,:) - SS_str_bout_baseline_all(counter_pCell, 1);
    mean_window_SS_lick_g(counter_pCell,1) = nanmean(SS_lick_g(counter_pCell, :));
    mean_change_window_SS_lick_g(counter_pCell,1) = nanmean(change_SS_lick_g(counter_pCell, :));
    SS_lick_g(counter_pCell,:) = SS_lick_g(counter_pCell,:)/mean_window_SS_lick_g(counter_pCell,1);
    
    SS_lick_l(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_vmin.train_data_logic_SS_l_corr * 100), 15, 'sgolay', 2);
    change_SS_lick_l(counter_pCell,:) = SS_lick_l(counter_pCell,:) - SS_str_bout_baseline_l(counter_pCell, 1);
    mean_window_SS_lick_l(counter_pCell,1) = nanmean(SS_lick_l(counter_pCell, :));
    mean_change_window_SS_lick_l(counter_pCell,1) = nanmean(change_SS_lick_l(counter_pCell, :));
    SS_lick_l(counter_pCell,:) = SS_lick_l(counter_pCell,:)/mean_window_SS_lick_l(counter_pCell,1);
    
end

%%% LOAD LABELS FROM POPULATION STRUCTURE %%%
label_large_increase_all = POPULATION.dmax.label_large_increase_all;
label_medium_increase_all  = POPULATION.dmax.label_medium_increase_all;
label_small_increase_all = POPULATION.dmax.label_small_increase_all;
label_phase_all = POPULATION.dmax.label_phase_all;
label_small_decrease_all = POPULATION.dmax.label_small_decrease_all;
label_medium_decrease_all = POPULATION.dmax.label_medium_decrease_all;
label_large_decrease_all = POPULATION.dmax.label_large_decrease_all;
label_large_r = POPULATION.dmax.label_large_r;
label_medium_r = POPULATION.dmax.label_medium_r;
label_small_r = POPULATION.dmax.label_small_r;
label_none_dir = POPULATION.dmax.label_none_dir;
label_small_l = POPULATION.dmax.label_small_l;
label_medium_l = POPULATION.dmax.label_medium_l;
label_large_l = POPULATION.dmax.label_large_l;

label_rhythmic_all = POPULATION.dmax.label_rhythmic_all;
label_weakrhythmic_all = POPULATION.dmax.label_weakrhythmic_all;
label_rhythmic_l = POPULATION.dmax.label_rhythmic_l;
label_weakrhythmic_l = POPULATION.dmax.label_weakrhythmic_l;
label_rhythmic_g = POPULATION.dmax.label_rhythmic_g;
label_weakrhythmic_g = POPULATION.dmax.label_weakrhythmic_g;
label_rhythmic_r = POPULATION.dmax.label_rhythmic_r;
label_weakrhythmic_r = POPULATION.dmax.label_weakrhythmic_r;

label_lead_all = POPULATION.dmax.label_lead_all;
label_lead30to90_all  = POPULATION.dmax.label_lead30to90_all;
label_lead90to150_all  = POPULATION.dmax.label_lead90to150_all;
label_lag_all = POPULATION.dmax.label_lag_all;
label_lag30to90_all = POPULATION.dmax.label_lag30to90_all;
label_lag90to150_all = POPULATION.dmax.label_lag90to150_all;
label_phase_all = POPULATION.dmax.label_phase_all;
label_antiphase_all = POPULATION.dmax.label_antiphase_all;

label_lead_l = POPULATION.dmax.label_lead_l;
label_lead30to90_l  = POPULATION.dmax.label_lead30to90_l;
label_lead90to150_l  = POPULATION.dmax.label_lead90to150_l;
label_lag_l = POPULATION.dmax.label_lag_l;
label_lag30to90_l = POPULATION.dmax.label_lag30to90_l;
label_lag90to150_l = POPULATION.dmax.label_lag90to150_l;
label_phase_l = POPULATION.dmax.label_phase_l;
label_antiphase_l = POPULATION.dmax.label_antiphase_l;

label_lead_g = POPULATION.dmax.label_lead_g;
label_lead30to90_g  = POPULATION.dmax.label_lead30to90_g;
label_lead90to150_g  = POPULATION.dmax.label_lead90to150_g;
label_lag_g = POPULATION.dmax.label_lag_g;
label_lag30to90_g = POPULATION.dmax.label_lag30to90_g;
label_lag90to150_g = POPULATION.dmax.label_lag90to150_g;
label_phase_g = POPULATION.dmax.label_phase_g;
label_antiphase_g = POPULATION.dmax.label_antiphase_g;

label_lead_r = POPULATION.dmax.label_lead_r;
label_lead30to90_r  = POPULATION.dmax.label_lead30to90_r;
label_lead90to150_r  = POPULATION.dmax.label_lead90to150_r;
label_lag_r = POPULATION.dmax.label_lag_r;
label_lag30to90_r = POPULATION.dmax.label_lag30to90_r;
label_lag90to150_r = POPULATION.dmax.label_lag90to150_r;
label_phase_r = POPULATION.dmax.label_phase_r;
label_antiphase_r = POPULATION.dmax.label_antiphase_r;

%%% Mean d_tip based on labels
mean_d_tip_all = nanmean(d_tip_all);
mean_d_tip_l = nanmean(d_tip_l);
mean_d_tip_g = nanmean(d_tip_g);
mean_d_tip_r = nanmean(d_tip_r);

mean_d_tip_rhythmic_phase_all = nanmean(d_tip_all(label_rhythmic_all & label_phase_all,:));
se_d_tip_rhythmic_phase_all = nanstd(d_tip_all(label_rhythmic_all & label_phase_all,:))/sqrt(sum(label_rhythmic_all & label_phase_all));

mean_d_tip_rhythmic_antiphase_all = nanmean(d_tip_all(label_rhythmic_all & label_antiphase_all,:));
se_d_tip_rhythmic_antiphase_all = nanstd(d_tip_all(label_rhythmic_all & label_antiphase_all,:))/sqrt(sum(label_rhythmic_all & label_antiphase_all));

mean_d_tip_rhythmic_lead_all = nanmean(d_tip_all(label_rhythmic_all & label_lead_all,:));
se_d_tip_rhythmic_lead_all = nanstd(d_tip_all(label_rhythmic_all & label_lead_all,:))/sqrt(sum(label_rhythmic_all & label_lead_all));
mean_d_tip_rhythmic_lead30to90_all = nanmean(d_tip_all(label_rhythmic_all & label_lead30to90_all,:));
se_d_tip_rhythmic_lead30to90_all = nanstd(d_tip_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_d_tip_rhythmic_lead90to150_all = nanmean(d_tip_all(label_rhythmic_all & label_lead90to150_all,:));
se_d_tip_rhythmic_lead90to150_all = nanstd(d_tip_all(label_rhythmic_all & label_lead90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lead90to150_all));

mean_d_tip_rhythmic_lag_all = nanmean(d_tip_all(label_rhythmic_all & label_lag_all,:));
se_d_tip_rhythmic_lag_all = nanstd(d_tip_all(label_rhythmic_all & label_lag_all,:))/sqrt(sum(label_rhythmic_all & label_lag_all));
mean_d_tip_rhythmic_lead30to90_all = nanmean(d_tip_all(label_rhythmic_all & label_lead30to90_all,:));
se_d_tip_rhythmic_lead30to90_all = nanstd(d_tip_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_d_tip_rhythmic_lag90to150_all = nanmean(d_tip_all(label_rhythmic_all & label_lag90to150_all,:));
se_d_tip_rhythmic_lag90to150_all = nanstd(d_tip_all(label_rhythmic_all & label_lag90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lag90to150_all));

mean_d_tip_rhythmic_phase_l = nanmean(d_tip_l(label_rhythmic_l & label_phase_l,:));
se_d_tip_rhythmic_phase_l = nanstd(d_tip_l(label_rhythmic_l & label_phase_l,:))/sqrt(sum(label_rhythmic_l & label_phase_l));

mean_d_tip_rhythmic_antiphase_l = nanmean(d_tip_l(label_rhythmic_l & label_antiphase_l,:));
se_d_tip_rhythmic_antiphase_l = nanstd(d_tip_l(label_rhythmic_l & label_antiphase_l,:))/sqrt(sum(label_rhythmic_l & label_antiphase_l));

mean_d_tip_rhythmic_lead_l = nanmean(d_tip_l(label_rhythmic_l & label_lead_l,:));
se_d_tip_rhythmic_lead_l = nanstd(d_tip_l(label_rhythmic_l & label_lead_l,:))/sqrt(sum(label_rhythmic_l & label_lead_l));
mean_d_tip_rhythmic_lead30to90_l = nanmean(d_tip_l(label_rhythmic_l & label_lead30to90_l,:));
se_d_tip_rhythmic_lead30to90_l = nanstd(d_tip_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_d_tip_rhythmic_lead90to150_l = nanmean(d_tip_l(label_rhythmic_l & label_lead90to150_l,:));
se_d_tip_rhythmic_lead90to150_l = nanstd(d_tip_l(label_rhythmic_l & label_lead90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lead90to150_l));

mean_d_tip_rhythmic_lag_l = nanmean(d_tip_l(label_rhythmic_l & label_lag_l,:));
se_d_tip_rhythmic_lag_l = nanstd(d_tip_l(label_rhythmic_l & label_lag_l,:))/sqrt(sum(label_rhythmic_l & label_lag_l));
mean_d_tip_rhythmic_lead30to90_l = nanmean(d_tip_l(label_rhythmic_l & label_lead30to90_l,:));
se_d_tip_rhythmic_lead30to90_l = nanstd(d_tip_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_d_tip_rhythmic_lag90to150_l = nanmean(d_tip_l(label_rhythmic_l & label_lag90to150_l,:));
se_d_tip_rhythmic_lag90to150_l = nanstd(d_tip_l(label_rhythmic_l & label_lag90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lag90to150_l));

mean_d_tip_rhythmic_phase_g = nanmean(d_tip_g(label_rhythmic_g & label_phase_g,:));
se_d_tip_rhythmic_phase_g = nanstd(d_tip_g(label_rhythmic_g & label_phase_g,:))/sqrt(sum(label_rhythmic_g & label_phase_g));

mean_d_tip_rhythmic_antiphase_g = nanmean(d_tip_g(label_rhythmic_g & label_antiphase_g,:));
se_d_tip_rhythmic_antiphase_g = nanstd(d_tip_g(label_rhythmic_g & label_antiphase_g,:))/sqrt(sum(label_rhythmic_g & label_antiphase_g));

mean_d_tip_rhythmic_lead_g = nanmean(d_tip_g(label_rhythmic_g & label_lead_g,:));
se_d_tip_rhythmic_lead_g = nanstd(d_tip_g(label_rhythmic_g & label_lead_g,:))/sqrt(sum(label_rhythmic_g & label_lead_g));
mean_d_tip_rhythmic_lead30to90_g = nanmean(d_tip_g(label_rhythmic_g & label_lead30to90_g,:));
se_d_tip_rhythmic_lead30to90_g = nanstd(d_tip_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_d_tip_rhythmic_lead90to150_g = nanmean(d_tip_g(label_rhythmic_g & label_lead90to150_g,:));
se_d_tip_rhythmic_lead90to150_g = nanstd(d_tip_g(label_rhythmic_g & label_lead90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lead90to150_g));

mean_d_tip_rhythmic_lag_g = nanmean(d_tip_g(label_rhythmic_g & label_lag_g,:));
se_d_tip_rhythmic_lag_g = nanstd(d_tip_g(label_rhythmic_g & label_lag_g,:))/sqrt(sum(label_rhythmic_g & label_lag_g));
mean_d_tip_rhythmic_lead30to90_g = nanmean(d_tip_g(label_rhythmic_g & label_lead30to90_g,:));
se_d_tip_rhythmic_lead30to90_g = nanstd(d_tip_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_d_tip_rhythmic_lag90to150_g = nanmean(d_tip_g(label_rhythmic_g & label_lag90to150_g,:));
se_d_tip_rhythmic_lag90to150_g = nanstd(d_tip_g(label_rhythmic_g & label_lag90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lag90to150_g));

mean_d_tip_rhythmic_phase_r = nanmean(d_tip_r(label_rhythmic_r & label_phase_r,:));
se_d_tip_rhythmic_phase_r = nanstd(d_tip_r(label_rhythmic_r & label_phase_r,:))/sqrt(sum(label_rhythmic_r & label_phase_r));

mean_d_tip_rhythmic_antiphase_r = nanmean(d_tip_g(label_rhythmic_r & label_antiphase_r,:));
se_d_tip_rhythmic_antiphase_r = nanstd(d_tip_g(label_rhythmic_r & label_antiphase_r,:))/sqrt(sum(label_rhythmic_r & label_antiphase_r));

mean_d_tip_rhythmic_lead_r = nanmean(d_tip_r(label_rhythmic_r & label_lead_r,:));
se_d_tip_rhythmic_lead_r = nanstd(d_tip_r(label_rhythmic_r & label_lead_r,:))/sqrt(sum(label_rhythmic_r & label_lead_r));
mean_d_tip_rhythmic_lead30to90_r = nanmean(d_tip_r(label_rhythmic_r & label_lead30to90_r,:));
se_d_tip_rhythmic_lead30to90_r = nanstd(d_tip_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_d_tip_rhythmic_lead90to150_r = nanmean(d_tip_r(label_rhythmic_r & label_lead90to150_r,:));
se_d_tip_rhythmic_lead90to150_r = nanstd(d_tip_r(label_rhythmic_r & label_lead90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lead90to150_r));

mean_d_tip_rhythmic_lag_r = nanmean(d_tip_r(label_rhythmic_r & label_lag_r,:));
se_d_tip_rhythmic_lag_r = nanstd(d_tip_r(label_rhythmic_r & label_lag_r,:))/sqrt(sum(label_rhythmic_r & label_lag_r));
mean_d_tip_rhythmic_lead30to90_r = nanmean(d_tip_r(label_rhythmic_r & label_lead30to90_r,:));
se_d_tip_rhythmic_lead30to90_r = nanstd(d_tip_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_d_tip_rhythmic_lag90to150_r = nanmean(d_tip_r(label_rhythmic_r & label_lag90to150_r,:));
se_d_tip_rhythmic_lag90to150_r = nanstd(d_tip_r(label_rhythmic_r & label_lag90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lag90to150_r));

%%% Mean v_tip based on labels
mean_v_tip_all = nanmean(v_tip_all);
mean_v_tip_l = nanmean(v_tip_l);
mean_v_tip_g = nanmean(v_tip_g);
mean_v_tip_r = nanmean(v_tip_r);

mean_v_tip_rhythmic_phase_all = nanmean(v_tip_all(label_rhythmic_all & label_phase_all,:));
se_v_tip_rhythmic_phase_all = nanstd(v_tip_all(label_rhythmic_all & label_phase_all,:))/sqrt(sum(label_rhythmic_all & label_phase_all));

mean_v_tip_rhythmic_antiphase_all = nanmean(v_tip_all(label_rhythmic_all & label_antiphase_all,:));
se_v_tip_rhythmic_antiphase_all = nanstd(v_tip_all(label_rhythmic_all & label_antiphase_all,:))/sqrt(sum(label_rhythmic_all & label_antiphase_all));

mean_v_tip_rhythmic_lead_all = nanmean(v_tip_all(label_rhythmic_all & label_lead_all,:));
se_v_tip_rhythmic_lead_all = nanstd(v_tip_all(label_rhythmic_all & label_lead_all,:))/sqrt(sum(label_rhythmic_all & label_lead_all));
mean_v_tip_rhythmic_lead30to90_all = nanmean(v_tip_all(label_rhythmic_all & label_lead30to90_all,:));
se_v_tip_rhythmic_lead30to90_all = nanstd(v_tip_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_v_tip_rhythmic_lead90to150_all = nanmean(v_tip_all(label_rhythmic_all & label_lead90to150_all,:));
se_v_tip_rhythmic_lead90to150_all = nanstd(v_tip_all(label_rhythmic_all & label_lead90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lead90to150_all));

mean_v_tip_rhythmic_lag_all = nanmean(v_tip_all(label_rhythmic_all & label_lag_all,:));
se_v_tip_rhythmic_lag_all = nanstd(v_tip_all(label_rhythmic_all & label_lag_all,:))/sqrt(sum(label_rhythmic_all & label_lag_all));
mean_v_tip_rhythmic_lead30to90_all = nanmean(v_tip_all(label_rhythmic_all & label_lead30to90_all,:));
se_v_tip_rhythmic_lead30to90_all = nanstd(v_tip_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_v_tip_rhythmic_lag90to150_all = nanmean(v_tip_all(label_rhythmic_all & label_lag90to150_all,:));
se_v_tip_rhythmic_lag90to150_all = nanstd(v_tip_all(label_rhythmic_all & label_lag90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lag90to150_all));

mean_v_tip_rhythmic_phase_l = nanmean(v_tip_l(label_rhythmic_l & label_phase_l,:));
se_v_tip_rhythmic_phasell = nanstd(v_tip_l(label_rhythmic_l & label_phase_l,:))/sqrt(sum(label_rhythmic_l & label_phase_l));

mean_v_tip_rhythmic_antiphase_l = nanmean(v_tip_l(label_rhythmic_l & label_antiphase_l,:));
se_v_tip_rhythmic_antiphase_l = nanstd(v_tip_l(label_rhythmic_l & label_antiphase_l,:))/sqrt(sum(label_rhythmic_l & label_antiphase_l));

mean_v_tip_rhythmic_lead_l = nanmean(v_tip_l(label_rhythmic_l & label_lead_l,:));
se_v_tip_rhythmic_lead_l = nanstd(v_tip_l(label_rhythmic_l & label_lead_l,:))/sqrt(sum(label_rhythmic_l & label_lead_l));
mean_v_tip_rhythmic_lead30to90_l = nanmean(v_tip_l(label_rhythmic_l & label_lead30to90_l,:));
se_v_tip_rhythmic_lead30to90_l = nanstd(v_tip_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_v_tip_rhythmic_lead90to150_l = nanmean(v_tip_l(label_rhythmic_l & label_lead90to150_l,:));
se_v_tip_rhythmic_lead90to150_l = nanstd(v_tip_l(label_rhythmic_l & label_lead90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lead90to150_l));

mean_v_tip_rhythmic_lag_l = nanmean(v_tip_l(label_rhythmic_l & label_lag_l,:));
se_v_tip_rhythmic_lag_l = nanstd(v_tip_l(label_rhythmic_l & label_lag_l,:))/sqrt(sum(label_rhythmic_l & label_lag_l));
mean_v_tip_rhythmic_lead30to90_l = nanmean(v_tip_l(label_rhythmic_l & label_lead30to90_l,:));
se_v_tip_rhythmic_lead30to90_l = nanstd(v_tip_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_v_tip_rhythmic_lag90to150_l = nanmean(v_tip_l(label_rhythmic_l & label_lag90to150_l,:));
se_v_tip_rhythmic_lag90to150_l = nanstd(v_tip_l(label_rhythmic_l & label_lag90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lag90to150_l));

mean_v_tip_rhythmic_phase_g = nanmean(v_tip_g(label_rhythmic_g & label_phase_g,:));
se_v_tip_rhythmic_phase_g = nanstd(v_tip_g(label_rhythmic_g & label_phase_g,:))/sqrt(sum(label_rhythmic_g & label_phase_g));

mean_v_tip_rhythmic_antiphase_g = nanmean(v_tip_g(label_rhythmic_g & label_antiphase_g,:));
se_v_tip_rhythmic_antiphase_g = nanstd(v_tip_g(label_rhythmic_g & label_antiphase_g,:))/sqrt(sum(label_rhythmic_g & label_antiphase_g));

mean_v_tip_rhythmic_lead_g = nanmean(v_tip_g(label_rhythmic_g & label_lead_g,:));
se_v_tip_rhythmic_lead_g = nanstd(v_tip_g(label_rhythmic_g & label_lead_g,:))/sqrt(sum(label_rhythmic_g & label_lead_g));
mean_v_tip_rhythmic_lead30to90_g = nanmean(v_tip_g(label_rhythmic_g & label_lead30to90_g,:));
se_v_tip_rhythmic_lead30to90_g = nanstd(v_tip_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_v_tip_rhythmic_lead90to150_g = nanmean(v_tip_g(label_rhythmic_g & label_lead90to150_g,:));
se_v_tip_rhythmic_lead90to150_g = nanstd(v_tip_g(label_rhythmic_g & label_lead90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lead90to150_g));

mean_v_tip_rhythmic_lag_g = nanmean(v_tip_g(label_rhythmic_g & label_lag_g,:));
se_v_tip_rhythmic_lag_g = nanstd(v_tip_g(label_rhythmic_g & label_lag_g,:))/sqrt(sum(label_rhythmic_g & label_lag_g));
mean_v_tip_rhythmic_lead30to90_g = nanmean(v_tip_g(label_rhythmic_g & label_lead30to90_g,:));
se_v_tip_rhythmic_lead30to90_g = nanstd(v_tip_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_v_tip_rhythmic_lag90to150_g = nanmean(v_tip_g(label_rhythmic_g & label_lag90to150_g,:));
se_v_tip_rhythmic_lag90to150_g = nanstd(v_tip_g(label_rhythmic_g & label_lag90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lag90to150_g));

mean_v_tip_rhythmic_phase_r = nanmean(v_tip_r(label_rhythmic_r & label_phase_r,:));
se_v_tip_rhythmic_phase_r = nanstd(v_tip_r(label_rhythmic_r & label_phase_r,:))/sqrt(sum(label_rhythmic_r & label_phase_r));

mean_v_tip_rhythmic_antiphase_r = nanmean(v_tip_r(label_rhythmic_r & label_antiphase_r,:));
se_v_tip_rhythmic_antiphaser = nanstd(v_tip_r(label_rhythmic_r & label_antiphase_r,:))/sqrt(sum(label_rhythmic_r & label_antiphase_r));

mean_v_tip_rhythmic_lead_r = nanmean(v_tip_r(label_rhythmic_r & label_lead_r,:));
se_v_tip_rhythmic_lead_r = nanstd(v_tip_r(label_rhythmic_r & label_lead_r,:))/sqrt(sum(label_rhythmic_r & label_lead_r));
mean_v_tip_rhythmic_lead30to90_r = nanmean(v_tip_r(label_rhythmic_r & label_lead30to90_r,:));
se_v_tip_rhythmic_lead30to90_r = nanstd(v_tip_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_v_tip_rhythmic_lead90to150_r = nanmean(v_tip_r(label_rhythmic_r & label_lead90to150_r,:));
se_v_tip_rhythmic_lead90to150_r = nanstd(v_tip_r(label_rhythmic_r & label_lead90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lead90to150_r));

mean_v_tip_rhythmic_lag_r = nanmean(v_tip_r(label_rhythmic_r & label_lag_r,:));
se_v_tip_rhythmic_lag_r = nanstd(v_tip_r(label_rhythmic_r & label_lag_r,:))/sqrt(sum(label_rhythmic_r & label_lag_r));
mean_v_tip_rhythmic_lead30to90_r = nanmean(v_tip_r(label_rhythmic_r & label_lead30to90_r,:));
se_v_tip_rhythmic_lead30to90_r = nanstd(v_tip_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_v_tip_rhythmic_lag90to150_r = nanmean(v_tip_r(label_rhythmic_r & label_lag90to150_r,:));
se_v_tip_rhythmic_lag90to150_r = nanstd(v_tip_r(label_rhythmic_r & label_lag90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lag90to150_r));

%%% Mean change_SS_lick based on labels
mean_SS_lick_rhythmic_phase_all = nanmean(SS_lick_all(label_rhythmic_all & label_phase_all,:));
se_SS_lick_rhythmic_phase_all = nanstd(SS_lick_all(label_rhythmic_all & label_phase_all,:))/sqrt(sum(label_rhythmic_all & label_phase_all));

mean_SS_lick_rhythmic_antiphase_all = nanmean(SS_lick_all(label_rhythmic_all & label_antiphase_all,:));
se_SS_lick_rhythmic_antiphase_all = nanstd(SS_lick_all(label_rhythmic_all & label_antiphase_all,:))/sqrt(sum(label_rhythmic_all & label_antiphase_all));

mean_SS_lick_rhythmic_lead_all = nanmean(SS_lick_all(label_rhythmic_all & label_lead_all,:));
se_SS_lick_rhythmic_lead_all = nanstd(SS_lick_all(label_rhythmic_all & label_lead_all,:))/sqrt(sum(label_rhythmic_all & label_lead_all));
mean_SS_lick_rhythmic_lead30to90_all = nanmean(SS_lick_all(label_rhythmic_all & label_lead30to90_all,:));
se_SS_lick_rhythmic_lead30to90_all = nanstd(SS_lick_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_SS_lick_rhythmic_lead90to150_all = nanmean(SS_lick_all(label_rhythmic_all & label_lead90to150_all,:));
se_SS_lick_rhythmic_lead90to150_all = nanstd(SS_lick_all(label_rhythmic_all & label_lead90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lead90to150_all));

mean_SS_lick_rhythmic_lag_all = nanmean(SS_lick_all(label_rhythmic_all & label_lag_all,:));
se_SS_lick_rhythmic_lag_all = nanstd(SS_lick_all(label_rhythmic_all & label_lag_all,:))/sqrt(sum(label_rhythmic_all & label_lag_all));
mean_SS_lick_rhythmic_lead30to90_all = nanmean(SS_lick_all(label_rhythmic_all & label_lead30to90_all,:));
se_SS_lick_rhythmic_lead30to90_all = nanstd(SS_lick_all(label_rhythmic_all & label_lead30to90_all,:))/sqrt(sum(label_rhythmic_all & label_lead30to90_all));
mean_SS_lick_rhythmic_lag90to150_all = nanmean(SS_lick_all(label_rhythmic_all & label_lag90to150_all,:));
se_SS_lick_rhythmic_lag90to150_all = nanstd(SS_lick_all(label_rhythmic_all & label_lag90to150_all,:))/sqrt(sum(label_rhythmic_all & label_lag90to150_all));

mean_SS_lick_rhythmic_phase_l = nanmean(SS_lick_l(label_rhythmic_l & label_phase_l,:));
se_SS_lick_rhythmic_phase_l = nanstd(SS_lick_l(label_rhythmic_l & label_phase_l,:))/sqrt(sum(label_rhythmic_l & label_phase_l));

mean_SS_lick_rhythmic_antiphase_l = nanmean(SS_lick_l(label_rhythmic_l & label_antiphase_l,:));
se_SS_lick_rhythmic_antiphase_l = nanstd(SS_lick_l(label_rhythmic_l & label_antiphase_l,:))/sqrt(sum(label_rhythmic_l & label_antiphase_l));

mean_SS_lick_rhythmic_lead_l = nanmean(SS_lick_l(label_rhythmic_l & label_lead_l,:));
se_SS_lick_rhythmic_lead_l = nanstd(SS_lick_l(label_rhythmic_l & label_lead_l,:))/sqrt(sum(label_rhythmic_l & label_lead_l));
mean_SS_lick_rhythmic_lead30to90_l = nanmean(SS_lick_l(label_rhythmic_l & label_lead30to90_l,:));
se_SS_lick_rhythmic_lead30to90_l = nanstd(SS_lick_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_SS_lick_rhythmic_lead90to150_l = nanmean(SS_lick_l(label_rhythmic_l & label_lead90to150_l,:));
se_SS_lick_rhythmic_lead90to150_l = nanstd(SS_lick_l(label_rhythmic_l & label_lead90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lead90to150_l));

mean_SS_lick_rhythmic_lag_l = nanmean(SS_lick_l(label_rhythmic_l & label_lag_l,:));
se_SS_lick_rhythmic_lag_l = nanstd(SS_lick_l(label_rhythmic_l & label_lag_l,:))/sqrt(sum(label_rhythmic_l & label_lag_l));
mean_SS_lick_rhythmic_lead30to90_l = nanmean(SS_lick_l(label_rhythmic_l & label_lead30to90_l,:));
se_SS_lick_rhythmic_lead30to90_l = nanstd(SS_lick_l(label_rhythmic_l & label_lead30to90_l,:))/sqrt(sum(label_rhythmic_l & label_lead30to90_l));
mean_SS_lick_rhythmic_lag90to150_l = nanmean(SS_lick_l(label_rhythmic_l & label_lag90to150_l,:));
se_SS_lick_rhythmic_lag90to150_l = nanstd(SS_lick_l(label_rhythmic_l & label_lag90to150_l,:))/sqrt(sum(label_rhythmic_l & label_lag90to150_l));

mean_SS_lick_rhythmic_phase_g = nanmean(SS_lick_g(label_rhythmic_g & label_phase_g,:));
se_SS_lick_rhythmic_phase_g = nanstd(SS_lick_g(label_rhythmic_g & label_phase_g,:))/sqrt(sum(label_rhythmic_g & label_phase_g));

mean_SS_lick_rhythmic_antiphase_g = nanmean(SS_lick_g(label_rhythmic_g & label_antiphase_g,:));
se_SS_lick_rhythmic_antiphase_g = nanstd(SS_lick_g(label_rhythmic_g & label_antiphase_g,:))/sqrt(sum(label_rhythmic_g & label_antiphase_g));

mean_SS_lick_rhythmic_lead_g = nanmean(SS_lick_g(label_rhythmic_g & label_lead_g,:));
se_SS_lick_rhythmic_lead_g = nanstd(SS_lick_g(label_rhythmic_g & label_lead_g,:))/sqrt(sum(label_rhythmic_g & label_lead_g));
mean_SS_lick_rhythmic_lead30to90_g = nanmean(SS_lick_g(label_rhythmic_g & label_lead30to90_g,:));
se_SS_lick_rhythmic_lead30to90_g = nanstd(SS_lick_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_SS_lick_rhythmic_lead90to150_g = nanmean(SS_lick_g(label_rhythmic_g & label_lead90to150_g,:));
se_SS_lick_rhythmic_lead90to150_g = nanstd(SS_lick_g(label_rhythmic_g & label_lead90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lead90to150_g));

mean_SS_lick_rhythmic_lag_g = nanmean(SS_lick_g(label_rhythmic_g & label_lag_g,:));
se_SS_lick_rhythmic_lag_g = nanstd(SS_lick_g(label_rhythmic_g & label_lag_g,:))/sqrt(sum(label_rhythmic_g & label_lag_g));
mean_SS_lick_rhythmic_lead30to90_g = nanmean(SS_lick_g(label_rhythmic_g & label_lead30to90_g,:));
se_SS_lick_rhythmic_lead30to90_g = nanstd(SS_lick_g(label_rhythmic_g & label_lead30to90_g,:))/sqrt(sum(label_rhythmic_g & label_lead30to90_g));
mean_SS_lick_rhythmic_lag90to150_g = nanmean(SS_lick_g(label_rhythmic_g & label_lag90to150_g,:));
se_SS_lick_rhythmic_lag90to150_g = nanstd(SS_lick_g(label_rhythmic_g & label_lag90to150_g,:))/sqrt(sum(label_rhythmic_g & label_lag90to150_g));

mean_SS_lick_rhythmic_phase_r = nanmean(SS_lick_r(label_rhythmic_r & label_phase_r,:));
se_SS_lick_rhythmic_phase_r = nanstd(SS_lick_r(label_rhythmic_r & label_phase_r,:))/sqrt(sum(label_rhythmic_r & label_phase_r));

mean_SS_lick_rhythmic_antiphase_r = nanmean(SS_lick_r(label_rhythmic_r & label_antiphase_r,:));
se_SS_lick_rhythmic_antiphase_r = nanstd(SS_lick_r(label_rhythmic_r & label_antiphase_r,:))/sqrt(sum(label_rhythmic_r & label_antiphase_r));

mean_SS_lick_rhythmic_lead_r = nanmean(SS_lick_r(label_rhythmic_r & label_lead_r,:));
se_SS_lick_rhythmic_lead_r = nanstd(SS_lick_r(label_rhythmic_r & label_lead_r,:))/sqrt(sum(label_rhythmic_r & label_lead_r));
mean_SS_lick_rhythmic_lead30to90_r = nanmean(SS_lick_r(label_rhythmic_r & label_lead30to90_r,:));
se_SS_lick_rhythmic_lead30to90_r = nanstd(SS_lick_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_SS_lick_rhythmic_lead90to150_r = nanmean(SS_lick_r(label_rhythmic_r & label_lead90to150_r,:));
se_SS_lick_rhythmic_lead90to150_r = nanstd(SS_lick_r(label_rhythmic_r & label_lead90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lead90to150_r));

mean_SS_lick_rhythmic_lag_r = nanmean(SS_lick_r(label_rhythmic_r & label_lag_r,:));
se_SS_lick_rhythmic_lag_r = nanstd(SS_lick_r(label_rhythmic_r & label_lag_r,:))/sqrt(sum(label_rhythmic_r & label_lag_r));
mean_SS_lick_rhythmic_lead30to90_r = nanmean(SS_lick_r(label_rhythmic_r & label_lead30to90_r,:));
se_SS_lick_rhythmic_lead30to90_r = nanstd(SS_lick_r(label_rhythmic_r & label_lead30to90_r,:))/sqrt(sum(label_rhythmic_r & label_lead30to90_r));
mean_SS_lick_rhythmic_lag90to150_r = nanmean(SS_lick_r(label_rhythmic_r & label_lag90to150_r,:));
se_SS_lick_rhythmic_lag90to150_r = nanstd(SS_lick_r(label_rhythmic_r & label_lag90to150_r,:))/sqrt(sum(label_rhythmic_r & label_lag90to150_r));

%%% inds_span_corr %%%
inds_span_corr = -0.99 : 0.01 : 1;
%% PLOT LABELED TRACES
figure
subplot(3, 4, 1)
hold on
plot(inds_span_corr, mean_d_tip_rhythmic_phase_all(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_all + se_d_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_all - se_d_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_all(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_all + se_d_tip_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_all - se_d_tip_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lead_all(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_all + se_d_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_all - se_d_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lag_all(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_all + se_d_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_all - se_d_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-0.17 0.17])
ylabel('Displacement (mm)')
xlabel('Lick vmin (ms)')
title('Lick trace | All licks')

subplot(3, 4, 2)
hold on
plot(inds_span_corr, mean_d_tip_rhythmic_phase_l(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_l + se_d_tip_rhythmic_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_l - se_d_tip_rhythmic_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_l(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_l + se_d_tip_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_l - se_d_tip_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lead_l(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_l + se_d_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_l - se_d_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lag_l(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_l + se_d_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_l - se_d_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-0.17 0.17])
ylabel('Displacement (mm)')
xlabel('Lick vmin (ms)')
title('Lick trace | L licks')

subplot(3, 4, 3)
hold on
plot(inds_span_corr, mean_d_tip_rhythmic_phase_g(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_g + se_d_tip_rhythmic_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_g - se_d_tip_rhythmic_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_g(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_g + se_d_tip_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_g - se_d_tip_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lead_g(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_g + se_d_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_g - se_d_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lag_g(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_g + se_d_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_g - se_d_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-0.17 0.17])
ylabel('Displacement (mm)')
xlabel('Lick vmin (ms)')
title('Lick trace | G licks')

subplot(3, 4, 4)
hold on
plot(inds_span_corr, mean_d_tip_rhythmic_phase_r(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_r + se_d_tip_rhythmic_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_r - se_d_tip_rhythmic_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_r(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_r + se_d_tip_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_antiphase_r - se_d_tip_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lead_r(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_r + se_d_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_r - se_d_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_d_tip_rhythmic_lag_r(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_r + se_d_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_r - se_d_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
xlim([-0.17 0.17])
ylabel('Displacement (mm)')
xlabel('Lick vmin (ms)')
title('Lick trace | R licks')


subplot(3, 4, 5)
hold on
plot(inds_span_corr, mean_v_tip_rhythmic_phase_all(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_all + se_v_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_all - se_v_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_all(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_all + se_v_tip_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_all - se_v_tip_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lead_all(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_all +se_v_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_all - se_v_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lag_all(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_all + se_v_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_all - se_v_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
xlim([-0.17 0.17])
ylabel('Velocity (mm/s)')
xlabel('Lick vmin (ms)')
title('Lick trace | All licks')

subplot(3, 4, 6)
hold on
plot(inds_span_corr, mean_v_tip_rhythmic_phase_l(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_l + se_v_tip_rhythmic_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_l - se_v_tip_rhythmic_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_l(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_l + se_v_tip_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_l - se_v_tip_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lead_l(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_l + se_v_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_l - se_v_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lag_l(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_l + se_v_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_l - se_v_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
xlim([-0.17 0.17])
ylabel('Velocity (mm/s)')
xlabel('Lick vmin (ms)')
title('Lick trace | L licks')

subplot(3, 4, 7)
hold on
plot(inds_span_corr, mean_v_tip_rhythmic_phase_g(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_g + se_v_tip_rhythmic_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_g - se_v_tip_rhythmic_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_g(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_g + se_v_tip_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_g - se_v_tip_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lead_g(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_g + se_v_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_g - se_v_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lag_g(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_g + se_v_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_g - se_v_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
xlim([-0.17 0.17])
ylabel('Velocity (mm/s)')
xlabel('Lick vmin (ms)')
title('Lick trace | G licks')

subplot(3, 4,8)
hold on
plot(inds_span_corr, mean_v_tip_rhythmic_phase_r(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_r + se_v_tip_rhythmic_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_r - se_v_tip_rhythmic_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_r(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_r + se_v_tip_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_antiphase_r - se_v_tip_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lead_r(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_r + se_v_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_r - se_v_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_v_tip_rhythmic_lag_r(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_r + se_v_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_r - se_v_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
xlim([-0.17 0.17])
ylabel('Velocity (mm/s)')
xlabel('Lick vmin (ms)')
title('Lick trace | R licks')

subplot(3, 4, 9)
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all + se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all - se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_all + se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_all - se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all + se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all - se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all + se_SS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all -se_SS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
%     ylabel('Change in SS firing (spks/s)')
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmin (ms)')
title('SS firing | All licks')

subplot(3, 4, 10)
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l + se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l - se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_l + se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_l - se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l + se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l - se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l + se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l - se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
%     ylabel('Change in SS firing (spks/s)')
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmin (ms)')
title('SS firing | L licks')

subplot(3, 4, 11)
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g + se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g - se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_g(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_g + se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_g - se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g + se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g - se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g + se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g - se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
%     ylabel('Change in SS firing (spks/s)')
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmin (ms)')
title('SS firing | G licks')

subplot(3, 4, 12)
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r(:), 'k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r + se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r - se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r(:), 'm', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_r + se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_antirhythmic_phase_r - se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r(:), 'r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r + se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r - se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r(:), 'b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r + se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r - se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
%     ylabel('Change in SS firing (spks/s)')
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmin (ms)')
title('SS firing | R licks')
%% PLOT LABELED TRACES - Kinematics
figure
subplot(2, 4, 1)
hold on
yyaxis left
plot(inds_span_corr, mean_d_tip_all, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_all + se_d_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_phase_all - se_d_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_all(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_all + se_d_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_all - se_d_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_all(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_all + se_d_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_all - se_d_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
ylabel('Displacement (mm)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all + se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all - se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all + se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all - se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all + se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all - se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all + se_SS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all - se_SS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmin (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | All licks')

subplot(2, 4, 2)
hold on
yyaxis left
plot(inds_span_corr, mean_d_tip_l, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_l + se_d_tip_rhythmic_none_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_l - se_d_tip_rhythmic_none_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_l(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_l + se_d_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_l - se_d_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_l(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_l + se_d_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_l - se_d_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
ylabel('Displacement (mm)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l + se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l - se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l + se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l - se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l + se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l - se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l + se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l - se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmin (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | L licks')


subplot(2, 4, 3)
hold on
yyaxis left
plot(inds_span_corr, mean_d_tip_g, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_g + se_d_tip_rhythmic_none_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_g - se_d_tip_rhythmic_none_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_g(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_g + se_d_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_g - se_d_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_g(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_g + se_d_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_g - se_d_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
ylabel('Displacement (mm)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g + se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g - se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_g(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_g + se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_g - se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g + se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g - se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g + se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g - se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmin (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | G licks')


subplot(2, 4, 4)
hold on
yyaxis left
plot(inds_span_corr, mean_d_tip_r, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_r + se_d_tip_rhythmic_none_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_none_r - se_d_tip_rhythmic_none_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_r(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_r + se_d_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lead_r - se_d_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_r(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_r + se_d_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_d_tip_rhythmic_lag_r - se_d_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-1 18])
ylabel('Displacement (mm)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r + se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r - se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r + se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r - se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r + se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r - se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r + se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r - se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmin (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | R licks')

subplot(2, 4, 5)
hold on
yyaxis left
plot(inds_span_corr, mean_v_tip_all, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_all + se_v_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_phase_all - se_v_tip_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_all(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_all + se_v_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_all - se_v_tip_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_all(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_all + se_v_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_all - se_v_tip_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
ylabel('Velocity (mm/s)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all + se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_all - se_SS_lick_rhythmic_phase_all, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all + se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_all - se_SS_lick_rhythmic_antiphase_all, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all + se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_all - se_SS_lick_rhythmic_lead_all, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all + seSS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_all - se_SS_lick_rhythmic_lag_all, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmin (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | All licks')

subplot(2, 4, 6)
hold on
yyaxis left
plot(inds_span_corr, mean_v_tip_l, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_l + se_v_tip_rhythmic_none_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_l - se_v_tip_rhythmic_none_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_l(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_l + se_v_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_l - se_v_tip_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_l(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_l + se_v_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_l - se_v_tip_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
ylabel('Velocity (mm/s)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l + se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_l - se_SS_lick_rhythmic_phase_l, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l + se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l - se_SS_lick_rhythmic_antiphase_l, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l + se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_l - se_SS_lick_rhythmic_lead_l, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l + se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_l - se_SS_lick_rhythmic_lag_l, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmin (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | L licks')

subplot(2, 4, 7)
hold on
yyaxis left
plot(inds_span_corr, mean_v_tip_g, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_g + se_v_tip_rhythmic_none_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_g - se_v_tip_rhythmic_none_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_g(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_g + se_v_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_g - se_v_tip_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_g(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_g + se_v_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_g - se_v_tip_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
ylabel('Velocity (mm/s)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g + se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_g - se_SS_lick_rhythmic_phase_g, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l + se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_l - se_SS_lick_rhythmic_antiphase_g, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g + se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_g - se_SS_lick_rhythmic_lead_g, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g + se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_g - se_SS_lick_rhythmic_lag_g, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmin (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | G licks')

subplot(2, 4, 8)
hold on
yyaxis left
plot(inds_span_corr, mean_v_tip_r, 'Color' , [0.7 0.7 0.7], 'LineWidth', 3)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_r + se_v_tip_rhythmic_none_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_none_r - se_v_tip_rhythmic_none_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_r(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_r + se_v_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lead_r - se_v_tip_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_r(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_r + se_v_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_v_tip_rhythmic_lag_r - se_v_tip_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([-400 400])
ylabel('Velocity (mm/s)')
set(gca, 'YColor', [0.7 0.7 0.7])
yyaxis right;
hold on
plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r(:), '-k', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r + se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_phase_r - se_SS_lick_rhythmic_phase_r, 'k', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r(:), '-m', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r + se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_antiphase_r - se_SS_lick_rhythmic_antiphase_r, 'm', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r(:), '-r', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r + se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lead_r - se_SS_lick_rhythmic_lead_r, 'r', 'LineWidth', 1)
plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r(:), '-b', 'LineWidth', 2)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r + se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
%     plot(inds_span_corr, mean_SS_lick_rhythmic_lag_r - se_SS_lick_rhythmic_lag_r, 'b', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2);
ylim([0.85 1.15])
xlim([-0.17 0.17])
ylabel('Normalized SS firing (spks/s)')
xlabel('Lick vmin (s)')
set(gca, 'YColor', [0 0 0])
title('Lick trace | R licks')




















































%% OLD CODE BELOW

%% PLOT LABELED Traces
clearvars -except ALL_PCELL_COMPRESSED_DATA POPULATION

number_pCell =  length(ALL_PCELL_COMPRESSED_DATA);

%%% Load d_tip & SS %%%
for counter_pCell = 1 : number_pCell
    %%% d_tip %%%
    d_tip_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_all_corr);
    d_tip_g(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_grooming_corr);
    d_tip_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_r_corr);
    d_tip_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_l_corr);
    
    %%% d_tip_str_bout %%%
    d_tip_str_bout_all(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_str_bout);
    d_tip_str_bout_r(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_str_bout_r);
    d_tip_str_bout_l(counter_pCell,:) = nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.VID_d_tip_100_str_bout_l);
    
    %%% SS_lick %%%
    SS_lick_all(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_all_corr * 100), 15, 'sgolay', 2);
    SS_lick_r(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_r_corr * 100), 15, 'sgolay', 2);
    SS_lick_g(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_grooming_corr * 100), 15, 'sgolay', 2);
    SS_lick_l(counter_pCell,:) = smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_l_corr * 100), 15, 'sgolay', 2);
    
    %%% SS_lick_str_bout %%%
    SS_str_bout_all(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_str_bout * 100), 15, 'sgolay', 2);
    SS_str_bout_baseline_all(counter_pCell, 1) = nanmean(SS_str_bout_all(counter_pCell, 1:50));
    change_SS_str_bout_all(counter_pCell,:) = SS_str_bout_all(counter_pCell,:) - SS_str_bout_baseline_all(counter_pCell,1);
    
    SS_str_bout_l(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_str_bout_l * 100), 15, 'sgolay', 2);
    SS_str_bout_baseline_l(counter_pCell, 1) = nanmean(SS_str_bout_l(counter_pCell, 1:50));
    change_SS_str_bout_l(counter_pCell,:) = SS_str_bout_l(counter_pCell,:) - SS_str_bout_baseline_l(counter_pCell,1);
    
    SS_str_bout_r(counter_pCell,:) =  smooth(nanmean(ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_str_bout_r * 100), 15, 'sgolay', 2);
    SS_str_bout_baseline_r(counter_pCell, 1) = nanmean(SS_str_bout_r(counter_pCell, 1:50));
    change_SS_str_bout_r(counter_pCell,:) = SS_str_bout_r(counter_pCell,:) - SS_str_bout_baseline_r(counter_pCell,1);
    
end

%%% Load labels %%%
label_large_increase_all = POPULATION.bout.label_large_increase_all;
label_medium_increase_all = POPULATION.bout.label_medium_increase_all;
label_small_increase_all = POPULATION.bout.label_small_increase_all;
label_phase_all = POPULATION.bout.label_phase_all;
label_large_decrease_all = POPULATION.bout.label_large_decrease_all;
label_medium_decrease_all = POPULATION.bout.label_medium_decrease_all;
label_small_decrease_all = POPULATION.bout.label_small_decrease_all;


%%% Means change_SS_str_bout based on labels
mean_change_SS_str_bout_large_increase_all = nanmean(change_SS_str_bout_all(label_large_increase_all,:));
se_change_SS_str_bout_large_increase_all = nanstd(change_SS_str_bout_all(label_large_increase_all,:))/sqrt(sum(label_large_increase_all));

mean_change_SS_str_bout_medium_increase_all = nanmean(change_SS_str_bout_all(label_medium_increase_all,:));
se_change_SS_str_bout_medium_increase_all = nanstd(change_SS_str_bout_all(label_medium_increase_all,:))/sqrt(sum(label_medium_increase_all));

mean_change_SS_str_bout_small_increase_all = nanmean(change_SS_str_bout_all(label_small_increase_all,:));
se_change_SS_str_bout_small_increase_all = nanstd(change_SS_str_bout_all(label_small_increase_all,:))/sqrt(sum(label_small_increase_all));

mean_change_SS_str_bout_phase_all = nanmean(change_SS_str_bout_all(label_phase_all,:));
se_change_SS_str_bout_phase_all = nanstd(change_SS_str_bout_all(label_phase_all,:))/sqrt(sum(label_phase_all));

mean_change_SS_str_bout_large_decrease_all = nanmean(change_SS_str_bout_all(label_large_decrease_all,:));
se_change_SS_str_bout_large_decrease_all = nanstd(change_SS_str_bout_all(label_large_decrease_all,:))/sqrt(sum(label_large_decrease_all));

mean_change_SS_str_bout_medium_decrease_all = nanmean(change_SS_str_bout_all(label_medium_decrease_all,:));
se_change_SS_str_bout_medium_decrease_all = nanstd(change_SS_str_bout_all(label_medium_decrease_all,:))/sqrt(sum(label_medium_decrease_all));

mean_change_SS_str_bout_small_decrease_all = nanmean(change_SS_str_bout_all(label_small_decrease_all,:));
se_change_SS_str_bout_small_decrease_all = nanstd(change_SS_str_bout_all(label_small_decrease_all,:))/sqrt(sum(label_small_decrease_all));

%%% inds_span_bout %%%
inds_span_bout = -0.99 : 0.01 : 1;

%% FIGURES
figure
sgtitle(['pCells: ' num2str(number_pCell)])

subplot(1,2,1)
hold on
plot(inds_span_bout, mean_change_SS_str_bout_large_increase_all(:), 'r', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_large_increase_all + se_change_SS_str_bout_large_increase_all, 'r', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_large_increase_all - se_change_SS_str_bout_large_increase_all, 'r', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_medium_increase_all(:), 'b', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_medium_increase_all + se_change_SS_str_bout_medium_increase_all, 'b', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_medium_increase_all - se_change_SS_str_bout_medium_increase_all, 'b', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_small_increase_all(:), 'g', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_small_increase_all + se_change_SS_str_bout_small_increase_all, 'g', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_small_increase_all - se_change_SS_str_bout_small_increase_all, 'g', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_none_all(:), 'k', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_none_all + se_change_SS_str_bout_none_all, 'k', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_none_all - se_change_SS_str_bout_none_all, 'k', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_small_decrease_all(:), 'm', 'LineWidth', 2)
plot(inds_span_bout, mean_change_SS_str_bout_small_decrease_all + se_change_SS_str_bout_small_decrease_all, 'm', 'LineWidth', 1)
plot(inds_span_bout, mean_change_SS_str_bout_small_decrease_all - se_change_SS_str_bout_small_decrease_all, 'm', 'LineWidth', 1)
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
yline(30,'--r', 'LineWidth', 1);
yline(-30,'--r', 'LineWidth', 1);
yline(15,'--b', 'LineWidth', 1);
yline(-15,'--b', 'LineWidth', 1);
yline(5,'--g', 'LineWidth', 1);
yline(-5,'--m', 'LineWidth', 1);
xlim([-1 1])
ylim([-30 100])
xlabel('Bout start (s)')
ylabel('Change in SS firing (spks/s)')
title(['Effect of bout on SS firing | All'])

%% CELL RESPONSE PROPERTIES
%%% bouts combined %%%
for counter_pCell = 1: length(ALL_PCELL_COMPRESSED_DATA)
    data_bout(counter_pCell,:) =ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_str_bout * 100;
    data_baseline_bout(counter_pCell,1) = mean(data_bout(counter_pCell,1:50), 2);
    data_diff_bout(counter_pCell,:) = data_bout(counter_pCell,:) - data_baseline_bout(counter_pCell,1);
    data_smooth_bout(counter_pCell,:) = smooth(data_diff_bout(counter_pCell, :), 5, 'sgolay', 2);
end
[~, pca_bout, ~] = pca(data_smooth_bout(:,50:200));
label_1_bout = pca_bout(:,1) > 0;
label_2_bout = pca_bout(:,1) < 0 ;

mean_label_1_bout = mean(data_smooth_bout(label_1_bout, :), 1);
se_label_1_bout = std(data_smooth_bout(label_1_bout, :)) / sqrt(sum(label_1_bout));
mean_label_2_bout = mean(data_smooth_bout(label_2_bout, :), 1);
se_label_2_bout = std(data_smooth_bout(label_2_bout, :)) / sqrt(sum(label_2_bout));


%%% licks combined %%%
for counter_pCell = 1: length(ALL_PCELL_COMPRESSED_DATA)
    data_lick_g(counter_pCell,:) =ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_grooming * 100;
    data_lick_l(counter_pCell,:) =ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_l * 100;
    data_lick_r(counter_pCell,:) =ALL_PCELL_COMPRESSED_DATA(counter_pCell).raster_data_lick_onset.train_data_logic_SS_r * 100;
    data_diff_lick_g(counter_pCell,:) = data_lick_g(counter_pCell,:) - data_baseline_bout(counter_pCell,1);
    data_diff_lick_l(counter_pCell,:) = data_lick_l(counter_pCell,:) - data_baseline_bout(counter_pCell,1);
    data_diff_lick_r(counter_pCell,:) = data_lick_r(counter_pCell,:) - data_baseline_bout(counter_pCell,1);
    data_diff_lick(counter_pCell,:) = data_diff_lick_g(counter_pCell,:) * ((num_g(counter_pCell)/(num_g(counter_pCell)+num_l(counter_pCell)+num_r(counter_pCell)))) + ...
        data_diff_lick_l(counter_pCell,:) * ((num_l(counter_pCell)/(num_g(counter_pCell)+num_l(counter_pCell)+num_r(counter_pCell)))) + ...
        data_diff_lick_r(counter_pCell,:) * ((num_r(counter_pCell)/(num_g(counter_pCell)+num_l(counter_pCell)+num_r(counter_pCell))));
    data_smooth_lick(counter_pCell,:) = smooth(data_diff_lick(counter_pCell, :), 5, 'sgolay', 2);
end
[~, pca_lick, ~] = pca(data_smooth_lick(:,25:60));

mean_label_1_lick = mean(data_smooth_lick(label_1_bout, :), 1);
se_label_1_lick = std(data_smooth_lick(label_1_bout, :)) / sqrt(sum(label_1_bout));
mean_label_2_lick = mean(data_smooth_lick(label_2_bout, :), 1);
se_label_2_lick = std(data_smooth_lick(label_2_bout, :)) / sqrt(sum(label_2_bout));

figure
subplot(2, 3, 1)
hold on
plot(-0.99:0.01:1, data_smooth_bout(label_1_bout,:), 'r', 'LineWidth', 1 )
plot(-0.99:0.01:1, data_smooth_bout(label_2_bout,:), 'b', 'LineWidth', 1 )
% for counter_pCell = 1: length(ALL_PCELL_COMPRESSED_DATA)
%     plot(-0.99:0.01:1, data_smooth_bout(counter_pCell, :), 'LineWidth', 2 )
% end
xline(0, 'LineWidth', 2);
ylim([-60 160])
xlabel('Bout onset (s)')
ylabel('Change in SS firing (spks/s)')
title(['pCells : ' num2str(size(data_smooth_bout,1))])

subplot(2, 3, 2)
hold on
plot(pca_bout(label_1_bout, 1), pca_bout(label_1_bout, 2), 'or', 'MarkerSize', 5)
plot(pca_bout(label_2_bout, 1), pca_bout(label_2_bout, 2), 'ob', 'MarkerSize', 5)

xlabel('PC 1')
ylabel('PC 2 ')
title(['PCA | pCells: ' num2str(size(pca_bout,1))])

subplot(2, 3, 3)
hold on
plot(-0.99:0.01:1, mean_label_1_bout, 'r', 'LineWidth', 2 )
plot(-0.99:0.01:1, mean_label_1_bout + se_label_1_bout, 'r', 'LineWidth', 1 )
plot(-0.99:0.01:1, mean_label_1_bout - se_label_1_bout , 'r', 'LineWidth', 1 )
plot(-0.99:0.01:1, mean_label_2_bout, 'b', 'LineWidth', 2 )
plot(-0.99:0.01:1, mean_label_2_bout + se_label_2_bout, 'b', 'LineWidth', 1 )
plot(-0.99:0.01:1, mean_label_2_bout - se_label_2_bout , 'b', 'LineWidth', 1 )
xline(0, 'LineWidth', 2);
ylim([-60 160])
xlabel('Bout onset (s)')
ylabel('Change in SS firing (spks/s)')
title('Mean of clusters')

subplot(2, 3, 4)
hold on
plot(-299:10:300, data_smooth_lick(label_1_bout,:), 'r', 'LineWidth', 1 )
plot(-299:10:300, data_smooth_lick(label_2_bout,:), 'b', 'LineWidth', 1 )

xline(0, 'LineWidth', 2);
ylim([-60 160])
xlabel('Lick onset (ms)')
ylabel('Change in SS firing (spks/s)')
title(['pCells : ' num2str(size(data_smooth_lick,1))])

subplot(2, 3, 5)
hold on
plot(pca_lick(label_1_bout, 1), pca_lick(label_1_bout, 2), 'or', 'MarkerSize', 5)
plot(pca_lick(label_2_bout, 1), pca_lick(label_2_bout, 2), 'ob', 'MarkerSize', 5)
xlabel('PC 1')
ylabel('PC 2 ')
title(['PCA | pCells: ' num2str(size(pca_lick,1))])

subplot(2, 3, 6)
hold on
plot(-299:10:300, mean_label_1_lick, 'r', 'LineWidth', 2 )
plot(-299:10:300, mean_label_1_lick + se_label_1_lick, 'r', 'LineWidth', 1 )
plot(-299:10:300, mean_label_1_lick - se_label_1_lick , 'r', 'LineWidth', 1 )
plot(-299:10:300, mean_label_2_lick, 'b', 'LineWidth', 2 )
plot(-299:10:300, mean_label_2_lick + se_label_2_lick, 'b', 'LineWidth', 1 )
plot(-299:10:300, mean_label_2_lick - se_label_2_lick , 'b', 'LineWidth', 1 )
xline(0, 'LineWidth', 2);
ylim([-60 160])
xlabel('Lick onset (s)')
ylabel('Change in SS firing (spks/s)')
title('Mean of clusters')

%% SS BOUT PROPERTIES
%%% L-vermis L bouts %%%
for counter_pCell = 1 : length(inds_leftVermis)
    if (sum(contains(fields(ALL_PCELL_COMPRESSED_DATA(inds_leftVermis(counter_pCell)).raster_data_lick_onset), 'onset_l'), 1) > 0 )
        data_ll_bout(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(inds_leftVermis(counter_pCell)).raster_data_lick_onset.train_data_logic_SS_str_bout_l * 100;
        data_baseline_ll_bout(counter_pCell,1) = mean(data_ll_bout(counter_pCell, 1:50), 2)';
        data_diff_ll_bout(counter_pCell,:) = data_ll_bout(counter_pCell,:) - data_baseline_ll_bout(counter_pCell,1);
        data_smooth_ll_bout(counter_pCell,:) = smooth(data_diff_ll_bout(counter_pCell,:), 5, 'sgolay', 2)';
    end
end
%%% L-vermis R bouts %%%
for counter_pCell = 1 : length(inds_leftVermis)
    if (sum(contains(fields(ALL_PCELL_COMPRESSED_DATA(inds_leftVermis(counter_pCell)).raster_data_lick_onset), 'onset_l'), 1) > 0 )
        data_lr_bout(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(inds_leftVermis(counter_pCell)).raster_data_lick_onset.train_data_logic_SS_str_bout_r * 100;
        data_baseline_lr_bout(counter_pCell,1) = mean(data_lr_bout(counter_pCell, 1:50), 2)';
        data_diff_lr_bout(counter_pCell,:) = data_lr_bout(counter_pCell,:) - data_baseline_lr_bout(counter_pCell,1);
        data_smooth_lr_bout(counter_pCell,:) = smooth(data_diff_lr_bout(counter_pCell,:), 5, 'sgolay', 2)';
    end
end
%%% diff in LL vs LR %%%
diff_ll_lr = data_smooth_ll_bout - data_smooth_lr_bout;
%%% PCA diff in LL vs LR %%%
[~, pca_diff_ll_lr, ~] = pca(diff_ll_lr(:,50:200));
label_1_l = pca_diff_ll_lr(:,1) > 50;
label_2_l = pca_diff_ll_lr(:,1) > -50 & pca_diff_ll_lr(:,1) < 50  ;
label_3_l = pca_diff_ll_lr(:,1) < -50;

%%% means of labels_l %%%
mean_label_1_l = mean(diff_ll_lr(label_1_l, :), 1);
se_label_1_l = std(diff_ll_lr(label_1_l, :)) / sqrt(sum(label_1_l));
mean_label_2_l = mean(diff_ll_lr(label_2_l, :), 1);
se_label_2_l = std(diff_ll_lr(label_2_l, :)) / sqrt(sum(label_2_l));
mean_label_3_l = mean(diff_ll_lr(label_3_l, :) , 1);
se_label_3_l = std(diff_ll_lr(label_3_l, :)) / sqrt(sum(label_3_l));

%%% R-vermis L bouts %%%
for counter_pCell = 1 : length(inds_rightVermis)
    if (sum(contains(fields(ALL_PCELL_COMPRESSED_DATA(inds_rightVermis(counter_pCell)).raster_data_lick_onset), 'onset_l'), 1) > 0 )
        data_rl_bout(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(inds_rightVermis(counter_pCell)).raster_data_lick_onset.train_data_logic_SS_str_bout_l * 100;
        data_baseline_rl_bout(counter_pCell,1) = mean(data_rl_bout(counter_pCell, 1:50), 2)';
        data_diff_rl_bout(counter_pCell,:) = data_rl_bout(counter_pCell,:) - data_baseline_rl_bout(counter_pCell,1);
        data_smooth_rl_bout(counter_pCell,:) = smooth(data_diff_rl_bout(counter_pCell,:), 5, 'sgolay', 2)';
    end
end
%%% R-vermis R bouts %%%
for counter_pCell = 1 : length(inds_rightVermis)
    if (sum(contains(fields(ALL_PCELL_COMPRESSED_DATA(inds_rightVermis(counter_pCell)).raster_data_lick_onset), 'onset_l'), 1) > 0 )
        data_rr_bout(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(inds_rightVermis(counter_pCell)).raster_data_lick_onset.train_data_logic_SS_str_bout_r * 100;
        data_baseline_rr_bout(counter_pCell,1) = mean(data_rr_bout(counter_pCell, 1:50), 2)';
        data_diff_rr_bout(counter_pCell,:) = data_rr_bout(counter_pCell,:) - data_baseline_rr_bout(counter_pCell,1);
        data_smooth_rr_bout(counter_pCell,:) = smooth(data_diff_rr_bout(counter_pCell,:), 5, 'sgolay', 2)';
    end
end
%%% diff in LL vs LR %%%
diff_rl_rr = data_smooth_rl_bout - data_smooth_rr_bout;
%%% PCA diff in RL vs RR %%%
[~, pca_diff_rl_rr, ~] = pca(diff_rl_rr(:,50:200));
label_1_r = pca_diff_rl_rr(:,1) > 50;
label_2_r = pca_diff_rl_rr(:,1) > -50 & pca_diff_rl_rr(:,1) < 50  ;
label_3_r = pca_diff_rl_rr(:,1) < -50;

%%% means of labels_r %%%
mean_label_1_r = mean(diff_rl_rr(label_1_r, :), 1);
se_label_1_r = std(diff_rl_rr(label_1_r, :)) / sqrt(sum(label_1_r));
mean_label_2_r = mean(diff_rl_rr(label_2_r, :), 1);
se_label_2_r = std(diff_rl_rr(label_2_r, :)) / sqrt(sum(label_2_r));
mean_label_3_r = mean(diff_rl_rr(label_3_r, :) , 1);
se_label_3_r = std(diff_rl_rr(label_3_r, :)) / sqrt(sum(label_3_r));

figure
subplot(3,4,1)
hold on
plot(-0.99:0.01:1,data_smooth_ll_bout(label_1_l,:), 'r', 'LineWidth', 2)
plot(-0.99:0.01:1,data_smooth_ll_bout(label_2_l,:), 'g', 'LineWidth', 2)
plot(-0.99:0.01:1,data_smooth_ll_bout(label_3_l,:), 'b', 'LineWidth', 2)
% plot(-0.99:0.01:1,data_smooth_ll_bout, 'LineWidth', 2)
xline(0, 'LineWidth', 2);
xlabel('Bout onset (s)')
ylabel('Change in SS firing (spks/s)')
ylim([-60 160])
title(['L-verm, L-bouts | pCells: ' num2str(length(inds_leftVermis))])

subplot(3,4,2)
hold on
plot(-0.99:0.01:1,data_smooth_lr_bout(label_1_l,:), 'r', 'LineWidth', 2)
plot(-0.99:0.01:1,data_smooth_lr_bout(label_2_l,:), 'g', 'LineWidth', 2)
plot(-0.99:0.01:1,data_smooth_lr_bout(label_3_l,:), 'b', 'LineWidth', 2)
% plot(-0.99:0.01:1,data_smooth_lr_bout, 'LineWidth', 2)
xline(0, 'LineWidth', 2);
xlabel('Bout onset (s)')
ylabel('Change in SS firing (spks/s)')
ylim([-60 160])
title(['L-verm, R-bouts | pCells: ' num2str(length(inds_leftVermis))])

subplot(3,4,3)
hold on
plot(-0.99:0.01:1,data_smooth_rl_bout(label_1_r,:), 'r', 'LineWidth', 2)
plot(-0.99:0.01:1,data_smooth_rl_bout(label_2_r,:), 'g', 'LineWidth', 2)
plot(-0.99:0.01:1,data_smooth_rl_bout(label_3_r,:), 'b', 'LineWidth', 2)
% plot(-0.99:0.01:1,data_smooth_rl_bout, 'LineWidth', 2)
xline(0, 'LineWidth', 2);
xlabel('Bout onset (s)')
ylabel('Change in SS firing (spks/s)')
ylim([-60 160])
title(['R-verm, L-bouts | pCells: ' num2str(length(inds_rightVermis))])

subplot(3,4,4)
hold on
plot(-0.99:0.01:1,data_smooth_rr_bout(label_1_r,:), 'r', 'LineWidth', 2)
plot(-0.99:0.01:1,data_smooth_rr_bout(label_2_r,:), 'g', 'LineWidth', 2)
plot(-0.99:0.01:1,data_smooth_rr_bout(label_3_r,:), 'b', 'LineWidth', 2)
% plot(-0.99:0.01:1,data_smooth_rr_bout, 'LineWidth', 2)
xline(0, 'LineWidth', 2);
xlabel('Bout onset (s)')
ylabel('Change in SS firing (spks/s)')
title(['R-verm, R-bouts | pCells: ' num2str(length(inds_rightVermis))])

subplot(3,4,5)
hold on
plot(-0.99:0.01:1,diff_ll_lr(label_1_l,:), 'r', 'LineWidth', 2);
plot(-0.99:0.01:1,diff_ll_lr(label_2_l,:), 'g', 'LineWidth', 2);
plot(-0.99:0.01:1,diff_ll_lr(label_3_l,:), 'b', 'LineWidth', 2);
% plot(-0.99:0.01:1,diff_ll_lr, 'LineWidth', 2);
xline(0, 'LineWidth', 2);
xlabel('Bout onset (s)')
ylabel('Diff in change in SS firing (spks/s)')
ylim([-60 160])
title('L-verm(L-bouts - R-bouts)')

subplot(3,4,6)
hold on
plot(pca_diff_ll_lr(label_1_l, 1), pca_diff_ll_lr(label_1_l, 2), 'or', 'MarkerSize', 5)
plot(pca_diff_ll_lr(label_2_l, 1), pca_diff_ll_lr(label_2_l, 2), 'og', 'MarkerSize', 5)
plot(pca_diff_ll_lr(label_3_l, 1), pca_diff_ll_lr(label_3_l, 2), 'ob', 'MarkerSize', 5)
% for counter_pCell = 1 : length(inds_leftVermis)
% plot(pca_diff_ll_lr(counter_pCell, 1), pca_diff_ll_lr(counter_pCell, 2), 'o', 'MarkerSize', 5)
% end
xlabel('PC 1')
ylabel('PC 2 ')
title('PCA: L-verm(L-bouts - R-bouts) [-0.5 1]s')


subplot(3,4,7)
hold on
plot(-0.99:0.01:1,diff_rl_rr(label_1_r,:), 'r', 'LineWidth', 2);
plot(-0.99:0.01:1,diff_rl_rr(label_2_r,:), 'g', 'LineWidth', 2);
plot(-0.99:0.01:1,diff_rl_rr(label_3_r,:), 'b', 'LineWidth', 2);
% plot(-0.99:0.01:1,diff_rl_rr, 'LineWidth', 2);
xline(0, 'LineWidth', 2);
xlabel('Bout onset (s)')
ylabel('Diff in change in SS firing (spks/s)')
ylim([-60 160])
title('R-verm(L-bouts - R-bouts)')


subplot(3,4,8)
hold on
plot(pca_diff_rl_rr(label_1_r, 1), pca_diff_rl_rr(label_1_r, 2), 'or', 'MarkerSize', 5)
plot(pca_diff_rl_rr(label_2_r, 1), pca_diff_rl_rr(label_2_r, 2), 'og', 'MarkerSize', 5)
plot(pca_diff_rl_rr(label_3_r, 1), pca_diff_rl_rr(label_3_r, 2), 'ob', 'MarkerSize', 5)
% for counter_pCell = 1 : length(inds_rightVermis)
% plot(pca_diff_rl_rr(counter_pCell, 1), pca_diff_rl_rr(counter_pCell, 2), 'o', 'MarkerSize', 5)
% end
xlabel('PC 1')
ylabel('PC 2 ')
title('PCA: R-verm(L-bouts - R-bouts) [-0.5 1]s')

subplot(3,4, [9 10])
hold on
plot(-0.99:0.01:1,mean_label_1_l, 'r', 'LineWidth', 2)
plot(-0.99:0.01:1,mean_label_1_l + se_label_1_l, 'r', 'LineWidth', 1)
plot(-0.99:0.01:1,mean_label_1_l - se_label_1_l, 'r', 'LineWidth', 1)
plot(-0.99:0.01:1,mean_label_2_l, 'g', 'LineWidth', 2)
plot(-0.99:0.01:1,mean_label_2_l + se_label_2_l, 'g', 'LineWidth', 1)
plot(-0.99:0.01:1,mean_label_2_l - se_label_2_l, 'g', 'LineWidth', 1)
plot(-0.99:0.01:1,mean_label_3_l, 'b', 'LineWidth', 2)
plot(-0.99:0.01:1,mean_label_3_l + se_label_3_l, 'b', 'LineWidth', 1)
plot(-0.99:0.01:1,mean_label_3_l - se_label_3_l, 'b', 'LineWidth', 1)
xline(0, 'LineWidth', 2);
xlabel('Bout onset (s)')
ylabel('Diff in change in SS firing (spks/s)')
ylim([-60 160])
title('Mean of L-verm(L-bouts - R-bouts) clusters')

subplot(3,4, [11 12])
hold on
plot(-0.99:0.01:1,mean_label_1_r, 'r', 'LineWidth', 2)
plot(-0.99:0.01:1,mean_label_1_r + se_label_1_r, 'r', 'LineWidth', 1)
plot(-0.99:0.01:1,mean_label_1_r - se_label_1_r, 'r', 'LineWidth', 1)
plot(-0.99:0.01:1,mean_label_2_r, 'g', 'LineWidth', 2)
plot(-0.99:0.01:1,mean_label_2_r + se_label_2_r, 'g', 'LineWidth', 1)
plot(-0.99:0.01:1,mean_label_2_r - se_label_2_r, 'g', 'LineWidth', 1)
plot(-0.99:0.01:1,mean_label_3_r, 'b', 'LineWidth', 2)
plot(-0.99:0.01:1,mean_label_3_r + se_label_3_r, 'b', 'LineWidth', 1)
plot(-0.99:0.01:1,mean_label_3_r - se_label_3_r, 'b', 'LineWidth', 1)
xline(0, 'LineWidth', 2);
xlabel('Bout onset (s)')
ylabel('Diff in change in SS firing (spks/s)')
ylim([-60 160])
title('Mean of R-verm(L-bouts - R-bouts) clusters')

%% SCRATCH PAD
[xcorr_value,xcorr_lag] = xcorr(SS_lick_all(1,:), d_tip_all(1,:)); % cross-correlate signals with each other
[~,ind_max_xcross] = max(abs(xcorr_value));
sample_diff = xcorr_lag(ind_max_xcross);

POPULATION.bout.change_SS_str_bout_all = change_SS_str_bout_all;
POPULATION.bout.change_SS_str_bout_l = change_SS_str_bout_l;
POPULATION.bout.change_SS_str_bout_r = change_SS_str_bout_r;
POPULATION.bout.SS_str_bout_baseline_all = SS_str_bout_baseline_all;
POPULATION.bout.SS_str_bout_baseline_l = SS_str_bout_baseline_l;
POPULATION.bout.SS_str_bout_baseline_r = SS_str_bout_baseline_r;
POPULATION.bout.SS_str_bout_post_all = SS_str_bout_post_all;
POPULATION.bout.SS_str_bout_post_l = SS_str_bout_post_l;
POPULATION.bout.SS_str_bout_post_r = SS_str_bout_post_r;
POPULATION.bout.mean_change_SS_str_bout_all = mean_change_SS_str_bout_all;
POPULATION.bout.se_change_SS_str_bout_all = se_change_SS_str_bout_all;
POPULATION.bout.mean_change_SS_str_bout_l = mean_change_SS_str_bout_l;
POPULATION.bout.se_change_SS_str_bout_l = se_change_SS_str_bout_l;
POPULATION.bout.mean_change_SS_str_bout_r = mean_change_SS_str_bout_r;
POPULATION.bout.se_change_SS_str_bout_r = se_change_SS_str_bout_r;

