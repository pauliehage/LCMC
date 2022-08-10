%% function PGH_freq_analysis
function frequency_data = PGH_freq_analysis(LICKS_ALL_DATA,params,funcs,tags_of_interest)
% if no tags specified then use all
if nargin<4
    tag_id = [1:9]; % default
elseif nargin==4
    tag_id = tags_of_interest;
end

Fs = 1000; %Hz

ang_edges = params.lick.ang_edges;
ang_values = params.lick.ang_values;
tag_name_list = params.lick.tag_name_list;

tag_bin = 1 : length(tag_name_list);
%% compute SS and dm freq
tongue_ang_bin = discretize(LICKS_ALL_DATA.tongue_ang_max, ang_edges);
num_ang_bin = length(ang_edges) - 1;
num_tag_bin = length(tag_bin);


for counter_tag = 1 : num_tag_bin
    for counter_ang = 1 : num_ang_bin
        idx_tag = (LICKS_ALL_DATA.tag == tag_bin(counter_tag));
        idx_ang = (tongue_ang_bin == counter_ang);
        idx_ = idx_tag&idx_ang;
        if sum(idx_) ~= 0
            % SS
            SS_onset = ESN_smooth(mean(LICKS_ALL_DATA.neuro_SS_onset(:,idx_),2)*1000,3);
            SS_onset = SS_onset - mean(SS_onset); % remove dc component of SS signal
            [power_SS,frequency_SS_] = pspectrum(SS_onset,Fs);
            [~, ind] = max(power_SS);
            frequency_SS(counter_tag, counter_ang) = frequency_SS_(ind);
           
            % CS
            CS_onset = ESN_smooth(mean(LICKS_ALL_DATA.neuro_CS_onset(:,idx_),2)*1000,3);
            CS_onset = CS_onset - mean(CS_onset); % remove dc component of SS signal
            [power_CS,frequency_CS_] = pspectrum(CS_onset,Fs);
            [~, ind] = max(power_CS);
            frequency_CS(counter_tag, counter_ang) = frequency_CS_(ind);

            % dm
            dm_onset= mean(LICKS_ALL_DATA.tongue_dm_onset(:,idx_),2);
            dm_onset = dm_onset - mean(dm_onset); % remove dc component of dm signal
            [power_dm,frequency_dm_] = pspectrum(dm_onset,Fs);
            [~, ind] = max(power_dm);
            frequency_dm(counter_tag, counter_ang) = frequency_dm_(ind);
            lick_count(counter_tag, counter_ang)  = sum(idx_);
        else
            frequency_SS(counter_tag, counter_ang) = nan;
            frequency_CS(counter_tag, counter_ang) = nan;
            frequency_dm(counter_tag, counter_ang) = nan;
        end
    end
end

%% Weighted average accross tags
frequency_SS_avg = nan(size(lick_count( 1, :)));
frequency_CS_avg = nan(size(lick_count( 1, :)));
frequency_dm_avg = nan(size(lick_count( 1, :)));

lick_count_avg = zeros(size(lick_count( 1, :)));

tags_ang_avg = tag_id; % 'prim_success' tag 1 % 'corr_success' tag 4 % 'back_center_irrelev' tag 8

for tag_ = tags_ang_avg
    frequency_SS_avg  = nansum([frequency_SS_avg ; frequency_SS( tag_, :).*lick_count(tag_,:)]);
    frequency_CS_avg  = nansum([frequency_CS_avg ; frequency_CS( tag_, :).*lick_count(tag_,:)]);
    frequency_dm_avg  = nansum([frequency_dm_avg ; frequency_dm( tag_, :).*lick_count(tag_,:)]);
    lick_count_avg = lick_count_avg + lick_count(tag_, :);
end
    frequency_SS_avg = frequency_SS_avg ./ lick_count_avg;
    frequency_CS_avg = frequency_CS_avg ./ lick_count_avg;
    frequency_dm_avg = frequency_dm_avg ./ lick_count_avg;

%% Build frequency_data
frequency_data.frequency_SS = frequency_SS;
frequency_data.frequency_CS = frequency_CS;
frequency_data.frequency_dm = frequency_dm;

frequency_data.frequency_SS_avg = frequency_SS_avg;
frequency_data.frequency_CS_avg = frequency_CS_avg;
frequency_data.frequency_dm_avg = frequency_dm_avg;


% figure
% x_axis = -90:45:90;
% hold on
% plot(x_axis,frequency_dm_avg, 'k')
% plot(x_axis,frequency_SS_avg, 'b')
% plot(x_axis,frequency_CS_avg, 'r')
% xlabel('Direction (deg)')
% ylabel('Frequency (Hz)')
% ylim([0 inf]);
% xlim([x_axis(1) x_axis(end)])
% xticks([x_axis])
% ESN_Beautify_Plot(gcf, [20 10])

end


