% AUTHORS:
%   Paul Hage and Lucas S. Mandacaru Guerra
function [pupil_area, validity] = filter_pupil(pupil_area_raw,measurement_time,counter_session,path)
pupil_area = pupil_area_raw;
% cap negative values at 0
pupil_area(pupil_area<0) = 0;
% find lower prom threshold using kmeans on all peaks
[peaks,~,~,~] = findpeaks(-pupil_area, 'Annotate', 'extent');
[~, cent] = kmeans(peaks', 2,"Start",[median(peaks) ; 0]);
threshold = abs(round((cent(1) + cent(2))/2));
% detect troughs
[~,ind_peak_neg,half_width,~] = findpeaks(-pupil_area, 'MinPeakHeight', -threshold,'Annotate', 'extent');
width = half_width; width(isoutlier(width)) = nanmean(width(~isoutlier(width))); width = round(width);
for counter_peak = 1 : length(ind_peak_neg)
    deletion_window = -width(counter_peak):1:width(counter_peak);
    inds_nan = ind_peak_neg(counter_peak) + deletion_window;
    inds_nan(inds_nan<=0 | inds_nan > length(pupil_area)) = [];
    pupil_area(inds_nan) = nan;
end
pupil_area(pupil_area < threshold) = nan;
pupil_area_paul = pupil_area(:);
% find upper prom threshold using kmeans on all peaks
threshold_up = mean(pupil_area,'omitnan')+4*std(pupil_area,'omitnan');
pupil_area(pupil_area > threshold_up) = nan;
% find dA/dt
pupil_velocity = [];
pupil_velocity(1) = abs(pupil_area(2)-pupil_area(1))/(measurement_time(2)-measurement_time(1));
for i = 2:length(pupil_area)-1
    pupil_velocity(i) = max([abs(pupil_area(i)-pupil_area(i-1))/(measurement_time(i)-measurement_time(i-1)),abs(pupil_area(i+1)-pupil_area(i))/(measurement_time(i+1)-measurement_time(i))]);
end
pupil_velocity(end+1) = abs(pupil_area(end)-pupil_area(end-1))/(measurement_time(end)-measurement_time(end-1));
% find dA/dt threshold using kmeans on all peaks
[peaks,~,~,~] = findpeaks(pupil_velocity, 'Annotate', 'extent');
[~, cent] = kmeans(peaks', 3,"Start",[median(peaks) ; 0 ; 0]);
mean_velocities = sort(cent);
threshold_velocity = abs(round((mean_velocities(1) + mean_velocities(2))/2));
% find new lower/upper prom thresholds using kmeans on all peaks
[peaks,~,~,~] = findpeaks(-pupil_area, 'Annotate', 'extent');
[~, cent] = kmeans(peaks', 3,"Start",[median(peaks) ; 0; 0]);
mean_neg_peaks = sort(cent);
lower_threshold = abs(round((mean_neg_peaks(2) + mean_neg_peaks(3))/2));
[peaks,~,~,~] = findpeaks(pupil_area, 'Annotate', 'extent');
[~, cent] = kmeans(peaks', 3,"Start",[median(peaks) ; 0; 0]);
mean_pos_peaks = sort(cent);
upper_threshold = abs(round((mean_pos_peaks(2) + mean_pos_peaks(3))/2));
% detect troughs
[~,ind_peak,half_width,~] = findpeaks(pupil_velocity, 'MinPeakHeight', threshold_velocity, 'Annotate', 'extent');
width = half_width; width(isoutlier(width)) = nanmean(width(~isoutlier(width))); width = round(width);
for counter_peak = 1 : length(ind_peak)
    if pupil_area(ind_peak(counter_peak)) < lower_threshold || pupil_area(ind_peak(counter_peak)) > upper_threshold
        deletion_window = -width(counter_peak):1:width(counter_peak);
        inds_nan = ind_peak(counter_peak) + deletion_window;
        inds_nan(inds_nan<=0 | inds_nan > length(pupil_area)) = [];
        pupil_area(inds_nan) = nan;
    end
end

% build validity
validity = ones(length(pupil_area)); validity(isnan(pupil_area)) = 0;

% compare raw vs filtered pupil signal
plot_check = 1;
if plot_check == 1
    figure
    hold on
    plot(pupil_area_raw,'r')
    plot(pupil_area_paul,'b')
    plot(pupil_area,'k')
    ylabel('pupil area')
    xlabel('time')
    ylim([0 inf])
    title(['session:' num2str(counter_session)])
    legend('raw', 'filtered - paul', 'filtered - lucas')
%     image_name = append(path, 'Pupil_', num2str(counter_session), '.fig')
%     saveas(gcf,image_name)
    ESN_Beautify_Plot(gcf,[20,8])
    saveas(gcf, [path.out_path 'PUPIL' filesep path.path_data_monkey_sorted '_pupil_' num2str(counter_session) ], 'pdf');
    close all
end
end
