function PGH_analyze_neural_properties(NP)
NP = load('Z:\video_10TB\Paul\FN\population_neural_properties.mat');
NP = NP.population_neural_properties;
%% Baseline FR
clearvars -except NP
FR = NP.SS_firing_rate;

mean_FR = mean(FR);
sd_FR = std(FR);
se_FR = sd_FR/sqrt(length(FR));

figure
histogram(FR, 'FaceColor', 'k');
xline(mean_FR, '-r', 'LineWidth', 2);
xline(mean_FR + sd_FR, '--r', 'LineWidth', 1);
xline(mean_FR - sd_FR, '--r', 'LineWidth', 1);

xlabel('Firing rate (spks/s)')
ylabel('Number of neurons')
title(['N = ' num2str(length(FR)) ' | ' num2str(mean_FR) '+/-' num2str(sd_FR) ' (spks/s; mean+/-sd)'])

ESN_Beautify_Plot

%% Autoprobability
clearvars -except NP 

time_window = -49:50;
autoprob = NP.Corr_data_SS_SSxSS_AUTO;
autoprob(:,50) = 0;
FR = NP.SS_firing_rate;

%normalize data
for counter_cell = 1 : length(autoprob)
    norm_autoprob(counter_cell,:) = autoprob(counter_cell,:)./FR(counter_cell,:)*1e3;
end

% do not normalize
norm_autoprob = autoprob;

%plot norm_autoprob
figure
hold on
for counter_cell = 1 : length(norm_autoprob)
    plot(time_window,norm_autoprob(counter_cell,:))
end
xline(0);
xlabel('Time (ms)')
ylabel('Normalized prob.')
title(['Normalized autoprob | N = ' num2str(length(norm_autoprob))])

window = 1:size(norm_autoprob,2);
[umap, ~, ~, ~]=run_umap(norm_autoprob(:,window));
umap = umap(:, 1:2);

% Kmeans cluster check
cluster_check = 1;
if cluster_check == 1
    for counter_k = 1:10
        [label, centroid_SS_waveform, sum_distances] = kmeans(umap, counter_k);
        WCSS(counter_k) = sum(sum_distances);
    end
    plot(WCSS)
end

k =6;

[label, centroid_SS_waveform, sum_distances] = kmeans(umap, k);

for counter_k = 1 : k
    mean_label(counter_k,:) = nanmean(norm_autoprob(label == counter_k, :));
end

colors1 = ['or' ;'ob';'og' ;'oc'; 'om';'oy';'*r'; '*b'; '*g'];
colors2 = ['-r' ;'-b';'-g' ;'-c';'-m';'-y'; '-r'; '-b'; '-g'];


check_scatter = 1;
if check_scatter == 1
figure
hold on
plot(umap(:,1), umap(:,2), 'ok')
for counter_k = 1 : k
    plot(umap(label==counter_k,1), umap(label==counter_k,2),colors1(counter_k,:))
    hold on
end
xlabel('UMAP 1')
ylabel('UMAP 2')
title('UMAP')
end

figure
for counter_k = 1 : k
    subplot(1,k,counter_k)
    hold on
    norm_autoprob_label = norm_autoprob(label == counter_k, :);
    for counter_num_in_k = 1 : sum(label == counter_k)
        plot(norm_autoprob_label(counter_num_in_k,:), '-k')
    end
    plot( mean_label(counter_k,:), colors2(counter_k,:), 'LineWidth', 2)
        ylim([0 inf])

end

figure
hold on
for counter_k = 1 : k
    plot(time_window, mean_label(counter_k,:), colors2(counter_k,:), 'LineWidth', 2)
end
xlabel('Time (ms)')
ylabel('Normalized prob.')
title('Cluster autoprobability averages')

DATA.label_autoprob = label;

%% Waveform
close all; clc
clearvars -except NP DATA

waveform = NP.waveform;
time_window = length(waveform)/30000;

%flip upward spiked
for counter_cell = 1 : length(waveform)
    if max(waveform(counter_cell,:))>abs(min(waveform(counter_cell,:)))
        waveform(counter_cell,:) = waveform(counter_cell,:)* -1;        
    end
end

%normalize data
for counter_cell = 1 : length(waveform)
    norm_waveform(counter_cell,:) = waveform(counter_cell,:)/max(abs(waveform(counter_cell,:)));
end

% plot all normalized waveforms
figure
hold on
for counter_cell = 1 : length(waveform)
    plot(waveform((counter_cell), :))
end
xline(0);
title('Normalized waveform')

% plot metronome normalized waveforms
% norm_waveform_metronome = norm_waveform(DATA.ind_metronome,:);
% plot all normalized waveforms
figure
hold on
for counter_cell = 1 : size(norm_waveform_metronome,1)
    plot(norm_waveform_metronome((counter_cell), :))
end
xline(0);
title('Normalized metronome waveform')


% disclude metronome
norm_waveform = norm_waveform(DATA.ind_non_metronome,:);

% cluster waveform discluding
window = 1:size(norm_waveform,2);
[umap, ~, ~, ~]=run_umap(norm_waveform(:,window));
umap = umap(:, 1:2);

% Kmeans cluster check
cluster_check = 0;
if cluster_check == 1
    for counter_k = 1:10
        [label, centroid_SS_waveform, sum_distances] = kmeans(umap, counter_k);
        WCSS(counter_k) = sum(sum_distances);
    end
    plot(WCSS)
end

k =6;

[label, centroid_SS_waveform, sum_distances] = kmeans(umap, k);

for counter_k = 1 : k
    mean_label(counter_k,:) = nanmean(norm_waveform(label == counter_k, :));
end

colors1 = ['or' ;'ob';'og' ;'oc'; 'om';'oy';'*r'; '*b'; '*g'];
colors2 = ['-r' ;'-b';'-g' ;'-c';'-m';'-y'; '-r'; '-b'; '-g'];

check_scatter = 0;
if check_scatter == 1
    figure
    hold on
    plot(umap(:,1), umap(:,2), 'ok')
    for counter_k = 1 : k
        plot(umap(label==counter_k,1), umap(label==counter_k,2),colors1(counter_k,:))
        hold on
    end
    xlabel('UMAP 1')
    ylabel('UMAP 2')
    title('UMAP')
end


figure
hold on
for counter_k = 1 : k
    plot( mean_label(counter_k,:), colors2(counter_k,:), 'LineWidth', 2)
end
xlabel('Time (ms)')
ylabel('Normalized prob.')
sgtitle('Cluster waveform averages')

figure
for counter_k = 1 : k
    subplot(1,k,counter_k)
    hold on
    norm_waveform_label = norm_waveform(label == counter_k, :);
    for counter_num_in_k = 1 : sum(label == counter_k)
        plot(norm_waveform_label(counter_num_in_k,:), '-k')
    end
    plot( mean_label(counter_k,:), colors2(counter_k,:), 'LineWidth', 2)
%     xlim([0 15])
end

DATA.label_waveform = label;

%% Confusion matrix: autoprob vs waveform
figure
[mat,ind] = (build_confusion_matrix(DATA.label_autoprob  , DATA.label_waveform  ));
imagesc(build_confusion_matrix(DATA.label_autoprob  , DATA.label_waveform  ))
xtick = ind(:,2);
% xtick = [1 3 2 5 6 4];
ytick = ind(:,1);
% ytick = [2 3 5 4 6 1];
set(gca, 'XTickLabel', xtick, 'YtickLabel', ytick);
xlabel('Waveform')
ylabel('Autoprobability')
colorbar

title('Transformed confusion matrix')

ESN_Beautify_Plot

%% function build_confusion_matrix
function [trans_conf_mat, correspondence_mat] = build_confusion_matrix(clustering_id_1, clustering_id_2)
num_classes_1 = max(clustering_id_1);
num_classes_2 = max(clustering_id_2);
confusion_matrix = zeros(num_classes_1, num_classes_2);
for counter_1 = 1 : num_classes_1
    idx_1_ = (clustering_id_1 == counter_1);
    for counter_2 = 1 : num_classes_2
        idx_2_ = (clustering_id_2 == counter_2);
        num_overlap = sum( idx_1_.*idx_2_ );
        prob_overlap = num_overlap ./ min([sum(idx_1_) sum(idx_2_)]);
        confusion_matrix(counter_1, counter_2) = prob_overlap;
    end
end
[trans_conf_mat, ind_mat] = transform_confusion_matrix(confusion_matrix);

ind_mat_diag = diag(ind_mat);
conf_mat_diag = diag(trans_conf_mat);
[row_,col_] = ind2sub(size(ind_mat),ind_mat_diag);
correspondence_mat = [row_, col_, conf_mat_diag];

end

%% function transform_confusion_matrix(conf_mat)
function [trans_conf_mat, ind_mat]= transform_confusion_matrix(conf_mat)
num_iteration = min([size(conf_mat,1) size(conf_mat,2)]);
ind_mat = reshape((1:numel(conf_mat))', size(conf_mat));
trans_conf_mat = conf_mat; % nan(size(conf_mat));
for counter_iter = 1 : num_iteration
    mat_1_ = trans_conf_mat;
    if counter_iter > 1
        mat_1_(1:counter_iter-1, :) = nan;
        mat_1_(:, 1:counter_iter-1) = nan;
    end
    [~, ind_max] = max(mat_1_,[],'all','linear','omitnan');
    [row_,col_] = ind2sub(size(mat_1_),ind_max);
    mat_1_ = trans_conf_mat;
    
    mat_2_ = mat_1_;
    mat_2_(counter_iter, :) = mat_1_(row_, :);
    mat_2_(row_, :) = mat_1_(counter_iter, :);
    
    ind_mat_ = ind_mat;
    ind_mat_(counter_iter, :) = ind_mat(row_, :);
    ind_mat_(row_, :) = ind_mat(counter_iter, :);
    
    mat_3_ = mat_2_;
    mat_3_(:, counter_iter) = mat_2_(:, col_);
    mat_3_(:, col_) = mat_2_(:, counter_iter);
    trans_conf_mat = mat_3_;
    
    ind_mat = ind_mat_;
    ind_mat_(:, counter_iter) = ind_mat(:, col_);
    ind_mat_(:, col_) = ind_mat(:, counter_iter);
    ind_mat = ind_mat_;
end

end
end