% clc; clear all; close all;
warning('off','all')

%% LOAD DATA
% [file_name, file_path] = uigetfile([pwd], 'Select ALL_CELL_COMPRESSED_DATA file');
% load([file_path file_name ] )

%% ANALYSIS OF NEURAL PROPERTIES
% clearvars -except ALL_CELL_COMPRESSED_DATA POPULATION

% Load data
num_DCN_cell = length(ALL_CELL_COMPRESSED_DATA);
time_SS_waveform = -1.9:0.1:2;
time_SSxSS = -19:20;

for counter_Cell = 1 : num_DCN_cell
    id(counter_Cell, 1) = string(ALL_CELL_COMPRESSED_DATA(counter_Cell).id{1, 1});
    SS_waveform(counter_Cell,:) = ALL_CELL_COMPRESSED_DATA(counter_Cell).Neural_Properties_data.SS_waveform(41:80);   
    SS_waveform_normalized(counter_Cell,:) = SS_waveform(counter_Cell,:)/norm(SS_waveform(counter_Cell,:));
    if max(SS_waveform_normalized(counter_Cell,:)) > abs(min(SS_waveform_normalized(counter_Cell,:)))
    SS_waveform_normalized(counter_Cell,:) = SS_waveform_normalized(counter_Cell,:)*-1;
    end
    SSxSS(counter_Cell, :) = ALL_CELL_COMPRESSED_DATA(counter_Cell).Neural_Properties_data.Corr_data.SS_SSxSS_AUTO(31:70);
    SSxSS(counter_Cell, 20) = 0;
    SSxSS_normalized(counter_Cell,:) = SSxSS(counter_Cell,:)/norm(SSxSS(counter_Cell,:));
    SS_FR(counter_Cell,1) = ALL_CELL_COMPRESSED_DATA(counter_Cell).Neural_Properties_data.SS_firing_rate;  
end

% Build SS_ISI
for counter_Cell = 1 : num_DCN_cell
    length_SS_ISI(counter_Cell,:) = length(ALL_CELL_COMPRESSED_DATA(counter_Cell).Neural_Properties_data.SS_time);
end
SS_ISI_ = nan(num_DCN_cell, max(length_SS_ISI));
for counter_Cell = 1 : num_DCN_cell
    SS_ISI_(counter_Cell,1:length_SS_ISI(counter_Cell)) = [0 diff(ALL_CELL_COMPRESSED_DATA(counter_Cell).Neural_Properties_data.SS_time')];
end
SS_ISI_(SS_ISI_ < 0) = 0;
for counter_Cell = 1 : num_DCN_cell
SS_ISI(counter_Cell,:) = histcounts(SS_ISI_(counter_Cell,:),0:0.001:0.1, 'Normalization', 'probability');
SS_ISI_normalized(counter_Cell,:) = SS_ISI(counter_Cell,:)/norm(SS_ISI(counter_Cell,:));
end

% Combine features: waveform, SSxSS, and ISI
combined_features = [SS_waveform_normalized SSxSS_normalized SS_ISI_normalized];

% Find SS_waveform characteristics
for counter_Cell = 1: num_DCN_cell
[peak_SS_waveform(counter_Cell,1), ind_peak_SS_waveform(counter_Cell,1), width_SS_waveform(counter_Cell,1), prom_SS_waveform(counter_Cell, 1)] = findpeaks(-1*SS_waveform_normalized(counter_Cell,:), time_SS_waveform, 'MinPeakProminence', 0.05);
peak_SS_waveform(counter_Cell,1) = peak_SS_waveform(counter_Cell,1) * -1;
%findpeaks(-1*SS_waveform_normalized(counter_Cell,:),time_SS_waveform, 'MinPeakProminence', 0.05, 'Annotate', 'extent')
end

% Find SSxSS characteristics
for counter_Cell = 1: num_DCN_cell
[peak_SSxSS(counter_Cell,1), ind_peak_SSxSS(counter_Cell,1), width_SSxSS(counter_Cell,1), prom_SSxSS(counter_Cell, 1)] = findpeaks(-1*SSxSS_normalized(counter_Cell,:), time_SSxSS, 'MinPeakProminence', 0.05);
peak_SSxSS(counter_Cell,1) = peak_SSxSS(counter_Cell,1) * -1;
%findpeaks(-1*SSxSS_normalized(counter_Cell,:),time_SSxSS, 'MinPeakProminence', 0.05, 'Annotate', 'extent')
end

% Gamma fit on SS ISI
for counter_Cell = 1 : num_DCN_cell
[gamma_param_SS_ISI(counter_Cell,:)] = gamfit(SS_ISI_(counter_Cell,1:length_SS_ISI(counter_Cell)));
end

% % PCA
% [~, pca_SS_waveform, ~] = pca(SS_waveform_normalized);
% pca_SS_waveform = pca_SS_waveform(:,1:2);

% UMAP
% addpath('C:\Users\Paul\OneDrive - Johns Hopkins\Shadmehr Lab\Code\umapFileExchange\umap')
[umap_SS_waveform, ~, ~, ~] = run_umap(SS_waveform_normalized);
[umap_SSxSS, ~, ~, ~] = run_umap(SSxSS_normalized);
[umap_SS_ISI, ~, ~, ~] = run_umap(SS_ISI_normalized);
[umap_combined_features, ~, ~, ~] = run_umap(combined_features);
umap_SS_waveform = umap_SS_waveform(:,1:2);
umap_SSxSS = umap_SSxSS(:,1:2);
umap_SS_ISI = umap_SS_ISI(:,1:2);
umap_combined_features = umap_combined_features(:,1:2);
close all;

% Kmeans cluster check
cluster_check = 0;
if cluster_check == 1
for counter_k = 1:10
    [label_SS_waveform, centroid_SS_waveform, sum_distances_SS_waveform] = kmeans(umap_SS_waveform, counter_k);
    WCSS_SS_waveform(counter_k) = sum(sum_distances_SS_waveform);
    [label_SSxSS, centroid_SSxSS, sum_distances_SSxSS] = kmeans(umap_SSxSS, counter_k);
    WCSS_SSxSS(counter_k) = sum(sum_distances_SSxSS);    
    [label_gamma_param_SS_ISI, centroid_gamma_param_SS_ISI,sum_distances_gamma_param_SS_ISI] = kmeans(gamma_param_SS_ISI,counter_k);
    WCSS_gamma_param_SS_ISI(counter_k) = sum(sum_distances_gamma_param_SS_ISI);   
    [label_SS_ISI, centroid_SS_ISI, sum_distances_SS_ISI] = kmeans(umap_SS_ISI, counter_k);
    WCSS_SS_ISI(counter_k) = sum(sum_distances_SS_ISI);
    [label_combined_features, centroid_combined_features, sum_distances_combined_features] = kmeans(umap_combined_features, counter_k);
    WCSS_combined_features(counter_k) = sum(sum_distances_combined_features); 
    
end
% Plot: WCSS x number of clusters
figure
subplot(5,1,1)
plot(1:10, WCSS_SS_waveform, '-k', 'LineWidth', 2)
xlabel('Number of clusters')
ylabel('WCSS')
title('SS waveform: WCSS vs number of clusters')
subplot(5,1,2)
plot(1:10, WCSS_SSxSS, '-k', 'LineWidth', 2)
xlabel('Number of clusters')
ylabel('WCSS')
title('SSxSS: WCSS vs number of clusters')
subplot(5,1,3)
plot(1:10, WCSS_SS_ISI, '-k', 'LineWidth', 2)
xlabel('Number of clusters')
ylabel('WCSS')
title('SS ISI: WCSS vs number of clusters')
subplot(5,1,4)
plot(1:10, WCSS_gamma_param_SS_ISI, '-k', 'LineWidth', 2)
xlabel('Number of clusters')
ylabel('WCSS')
title('SS ISI (gamma fit params): WCSS vs number of clusters')
subplot(5,1,5)
plot(1:10, WCSS_combined_features, '-k', 'LineWidth', 2)
xlabel('Number of clusters')
ylabel('WCSS')
title('Combined Features: WCSS vs number of clusters')
end

num_clusters = 5;

%Kmeans with selected number of clusters
[label_SS_waveform, centroid_SS_waveform, sum_distances_SS_waveform] = kmeans(umap_SS_waveform, 5);
[label_SSxSS, centroid_SSxSS, sum_distances_SSxSS] = kmeans(umap_SSxSS, num_clusters);
[label_gamma_param_SS_ISI, centroid_gamma_param_SS_ISI, sum_distances_gamma_param_SS_ISI] = kmeans(gamma_param_SS_ISI, num_clusters);
[label_SS_ISI, centroid_SS_ISI, sum_distances_SS_ISI] = kmeans(umap_SS_ISI, num_clusters);
[label_combined_features, centroid_combined_features, sum_distances_combined_features] = kmeans(umap_combined_features, num_clusters);

% Compute mean of SS_waveform clusters
for counter_cluster = 1:num_clusters
    eval(['mean_SS_waveform_' num2str(counter_cluster) ' = nanmean(SS_waveform_normalized(label_SS_waveform == ' num2str(counter_cluster) ',:));']);
    eval(['se_SS_waveform_' num2str(counter_cluster) ' = nanstd(SS_waveform_normalized(label_SS_waveform == ' num2str(counter_cluster) ',:)) / sqrt(sum(label_SS_waveform == ' num2str(counter_cluster) '));']);
    eval(['mean_SS_waveform_FR_' num2str(counter_cluster) ' = nanmean(SS_FR(label_SS_waveform == ' num2str(counter_cluster) '));']);
    eval(['se_SS_waveform_FR_' num2str(counter_cluster) ' = nanstd(SS_FR(label_SS_waveform == ' num2str(counter_cluster) ')) / sqrt(sum(label_SS_waveform ==' num2str(counter_cluster) '));']);
end

% Compute mean of SSxSS clusters
for counter_cluster = 1:num_clusters
    eval(['mean_SSxSS_' num2str(counter_cluster) ' = nanmean(SSxSS_normalized(label_SSxSS == ' num2str(counter_cluster) ',:));'])
    eval(['se_SSxSS_' num2str(counter_cluster) ' = nanstd(SSxSS_normalized(label_SSxSS == ' num2str(counter_cluster) ',:)) / sqrt(sum(label_SSxSS == ' num2str(counter_cluster) '));'])
    eval(['mean_SSxSS_FR_' num2str(counter_cluster) ' = nanmean(SS_FR(label_SSxSS == ' num2str(counter_cluster) '));']);
    eval(['se_SSxSS_FR_' num2str(counter_cluster) ' = nanstd(SS_FR(label_SSxSS == ' num2str(counter_cluster) ')) / sqrt(sum(label_SSxSS ==' num2str(counter_cluster) '));']);
end

% Compute mean of SS_ISI clusters
for counter_cluster = 1:num_clusters
    eval(['mean_SS_ISI_' num2str(counter_cluster) ' = nanmean(SS_ISI_normalized(label_SS_ISI == ' num2str(counter_cluster) ',:));'])
    eval(['se_SS_ISI_' num2str(counter_cluster) ' = nanstd(SS_ISI_normalized(label_SS_ISI == ' num2str(counter_cluster) ',:)) / sqrt(sum(label_SS_ISI == ' num2str(counter_cluster) '));'])
    eval(['mean_SS_ISI_FR_' num2str(counter_cluster) ' = nanmean(SS_FR(label_SS_ISI == ' num2str(counter_cluster) '));']);
    eval(['se_SS_ISI_FR_' num2str(counter_cluster) ' = nanstd(SS_FR(label_SS_ISI == ' num2str(counter_cluster) ')) / sqrt(sum(label_SS_ISI ==' num2str(counter_cluster) '));']);
end

% % Compute mean of SS_ISI_gamma clusters
% mean_SS_ISI_gamma_histogram_1 = nanmean(SS_ISI_normalized(label_gamma_param_SS_ISI ==1,:));
% se_SS_ISI_gamma_histogram_1 = nanstd(SS_ISI_normalized(label_gamma_param_SS_ISI ==1,:))/sqrt(sum(label_gamma_param_SS_ISI==1));

% Compute mean of combined_features
for counter_cluster = 1:num_clusters
    eval(['mean_combined_features_' num2str(counter_cluster) ' = nanmean(combined_features(label_combined_features == ' num2str(counter_cluster) ',:));'])
    eval(['se_combined_features_' num2str(counter_cluster) ' = nanstd(combined_features(label_combined_features == ' num2str(counter_cluster) ',:)) / sqrt(sum(label_combined_features == ' num2str(counter_cluster) '));'])
    eval(['mean_combined_features_FR_' num2str(counter_cluster) ' = nanmean(SS_FR(label_combined_features == ' num2str(counter_cluster) '));']);
    eval(['se_combined_features_FR_' num2str(counter_cluster) ' = nanstd(SS_FR(label_combined_features == ' num2str(counter_cluster) ')) / sqrt(sum(label_combined_features ==' num2str(counter_cluster) '));']);
end

% Compute xxx duration of mean SS_waveform clusters
ind_spike_onset = find(abs(diff(mean_SS_waveform_1)) > 0.01, 1, 'first') ;
[~,ind_spike_trough] = min(mean_SS_waveform_1);
%% Cell select
DCN_cell = [17];
%% Plot: data 
fig_handle = figure;
fig_handle.WindowState = 'maximized';

sgtitle([num2str(num_DCN_cell) ' DCN cells | Selected cell: ' num2str(DCN_cell)])

subplot(4,2,1)
hold on
plot(1:num_DCN_cell,SS_FR, '.k', 'MarkerSize', 11)
plot(DCN_cell, SS_FR(DCN_cell), '.r', 'MarkerSize', 11)
xlabel('Cell #')
ylabel('Firing rate (spks/s)')
title('All Cells Baseline Firing Rate ')

subplot(4,2,2)
hold on
histogram(SS_FR, 'FaceColor', 'black')
% xline(SS_FR(DCN_cell), 'r', 'LineWidth', 2);
xline(nanmean(SS_FR), 'r', 'LineWidth', 2);
xline(nanmean(SS_FR)+nanstd(SS_FR), 'r', 'LineWidth', 1);
xline(nanmean(SS_FR)-nanstd(SS_FR), 'r', 'LineWidth', 1);

xlabel('Firing Rate (spks/s)')
ylabel('Count')
title('Baseline Firing Rate Distribution ')

subplot(4,2,3)
xlabel_waveform = -1.9:0.1:2;
hold on
for counter_Cell = 1 : num_DCN_cell
plot(xlabel_waveform,SS_waveform_normalized(counter_Cell,:),'-k', 'LineWidth', 0.5)
% plot(xlabel_waveform,SS_waveform_normalized(DCN_cell,:),'r', 'LineWidth', 2)
xlabel('Time (ms)')
ylabel('Normalized uV')
end
title('Normalized SS Waveform')

subplot(4,2,4)
hold on
plot(width_SS_waveform, SS_FR, '.k', 'MarkerSize', 11)
% plot(width_SS_waveform(DCN_cell), SS_FR(DCN_cell), '.r', 'MarkerSize', 11)
xlabel('Spike waveform half-peak width (ms)')
ylabel('Baseline firing rate (spks/s)')
% hold on
% for counter_Cell = 1 : num_DCN_cell
% plot(umap_SS_waveform(counter_Cell,1),umap_SS_waveform(counter_Cell,2),'.k', 'MarkerSize', 11)
% plot(umap_SS_waveform(DCN_cell,1),umap_SS_waveform(DCN_cell,2),'.r', 'MarkerSize', 11)
% end
% xlabel('UMAP 1')
% ylabel('UMAP 2')
% title('UMAP Normalized SS Waveform')

subplot(4,2,5)
xlabel_SSxSS = -19:20;
hold on
for counter_Cell = 1 : num_DCN_cell
plot(xlabel_SSxSS,SSxSS_normalized(counter_Cell,:),'-k', 'LineWidth', 0.5)
% plot(xlabel_SSxSS,SSxSS_normalized(DCN_cell,:),'r', 'LineWidth', 2)
xlabel('Time (ms)')
ylabel('Normalized Probability')
end
title('Normalized SSxSS')

subplot(4,2,6)
hold on
plot(width_SSxSS, SS_FR, '.k', 'MarkerSize', 11)
% plot(width_SSxSS(DCN_cell), SS_FR(DCN_cell), '.r', 'MarkerSize', 11)
xlabel('SSxSS half-peak width (ms)')
ylabel('Baseline firing rate (spks/s)')
% hold on
% for counter_Cell = 1 : num_DCN_cell
% plot(umap_SSxSS(counter_Cell,1),umap_SSxSS(counter_Cell,2),'.k', 'MarkerSize', 11)
% plot(umap_SSxSS(DCN_cell,1),umap_SSxSS(DCN_cell,2),'.r', 'MarkerSize', 11)
% end
% xlabel('UMAP 1')
% ylabel('UMAP 2')
% title('UMAP Normalized SSxSS')

subplot(4,2,7)
xlabel_SS_ISI =(0:0.001:0.099);
hold on
for counter_Cell = 1 : num_DCN_cell
% histogram(SS_ISI(counter_Cell,:),0:0.001:0.1,'DisplayStyle', 'stairs', 'Normalization', 'probability', 'EdgeColor', 'black', 'LineWidth', 0.5)
plot(xlabel_SS_ISI,SS_ISI_normalized(counter_Cell, :), '-k','LineWidth', 0.5 );
end
% histogram(SS_ISI(DCN_cell,:),0:0.001:0.1,'DisplayStyle', 'stairs', 'Normalization', 'probability', 'EdgeColor', 'red', 'LineWidth', 2)
% plot(xlabel_SS_ISI,SS_ISI_normalized(DCN_cell, :), 'r','LineWidth', 2 );
xlim([0 0.1])
xlabel('Normalized ISI (s)')
ylabel('Probability')
title('Normalized SS ISI distribution')

subplot(4,2,8)
hold on
plot(nanmean(SS_ISI,2), SS_FR, '.k', 'MarkerSize', 10)
% plot(nanmean(SS_ISI(DCN_cell),2), SS_FR(DCN_cell), '.r', 'MarkerSize', 10)
xlabel('Mean ISI (s)')
ylabel('Baseline firing rate (spks/s)')
% hold on
% for counter_Cell = 1 : num_DCN_cell
% plot(umap_SS_ISI(counter_Cell,1), umap_SS_ISI(counter_Cell,2),'.k', 'MarkerSize', 11)
% end
% plot(umap_SS_ISI(DCN_cell,1), umap_SS_ISI(counter_Cell,2),'.r', 'MarkerSize', 11)
% xlabel('UMAP 1')
% ylabel('UMAP 2')
% title('UMAP Normalized SS ISI')


% subplot(4,2,8)
% hold on
% for counter_Cell = 1 : num_DCN_cell
% plot(gamma_param_SS_ISI(counter_Cell,1), gamma_param_SS_ISI(counter_Cell,2),'.k', 'MarkerSize', 11)
% end
% plot(gamma_param_SS_ISI(DCN_cell,1), gamma_param_SS_ISI(counter_Cell,2),'.r', 'MarkerSize', 11)
% xlabel('Shape param (k)')
% ylabel('Scale param (\theta)')
% title('SS ISI Gamma Parameters')

ESN_Beautify_Plot
%% Plot: clustered data 
fig_handle = figure;
% fig_handle.WindowState = 'maximized';
sgtitle('Kmeans Clustering: Waveforms, SSxSS, and ISI')
subplot(3,4,1)
xlabel_waveform = -1.9:0.1:2;
hold on
if length(se_SS_waveform_1) > 1 
shade(xlabel_waveform, mean_SS_waveform_1 + se_SS_waveform_1, 'r', xlabel_waveform, mean_SS_waveform_1 - se_SS_waveform_1, 'r',  'LineWidth', 0.1, 'FillColor', 'red', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(xlabel_waveform, mean_SS_waveform_1, 'r')
if length(se_SS_waveform_2) > 1 
shade(xlabel_waveform, mean_SS_waveform_2 + se_SS_waveform_2, '-g', xlabel_waveform, mean_SS_waveform_2 - se_SS_waveform_2, '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(xlabel_waveform, mean_SS_waveform_2, '-g')
if length(se_SS_waveform_3) > 1 
shade(xlabel_waveform, mean_SS_waveform_3 + se_SS_waveform_3, 'b', xlabel_waveform, mean_SS_waveform_3 - se_SS_waveform_3, 'b',  'LineWidth', 0.1, 'FillColor', 'blue', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(xlabel_waveform, mean_SS_waveform_3, 'b')
if length(se_SS_waveform_4) > 1 
shade(xlabel_waveform, mean_SS_waveform_4 + se_SS_waveform_4, 'c', xlabel_waveform, mean_SS_waveform_4 - se_SS_waveform_4, 'c',  'LineWidth', 0.1, 'FillColor', 'cyan', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(xlabel_waveform, mean_SS_waveform_4, 'c')
if length(se_SS_waveform_5) > 1 
shade(xlabel_waveform, mean_SS_waveform_5 + se_SS_waveform_5, 'm', xlabel_waveform, mean_SS_waveform_5 - se_SS_waveform_5, 'm',  'LineWidth', 0.1, 'FillColor', 'magenta', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(xlabel_waveform, mean_SS_waveform_5, 'm')
xlabel('Time (ms)')
ylabel('Normalized uV')
title('Clustered Waveforms')

subplot(3,4,2)
hold on
plot(umap_SS_waveform(label_SS_waveform==1,1), umap_SS_waveform(label_SS_waveform==1,2), '.r', 'MarkerSize', 10)
plot(umap_SS_waveform(label_SS_waveform==2,1), umap_SS_waveform(label_SS_waveform==2,2), '.g', 'MarkerSize', 10 )
plot(umap_SS_waveform(label_SS_waveform==3,1), umap_SS_waveform(label_SS_waveform==3,2), '.b', 'MarkerSize', 10)
plot(umap_SS_waveform(label_SS_waveform==4,1), umap_SS_waveform(label_SS_waveform==4,2), '.c', 'MarkerSize', 10)
plot(umap_SS_waveform(label_SS_waveform==5,1), umap_SS_waveform(label_SS_waveform==5,2), '.m', 'MarkerSize', 10)
xlabel('UMAP 1')
ylabel('UMAP 2')
title('Clustered UMAP: Waveform')

subplot(3,4,3)
hold on
bar_1 = bar(1, mean_SS_waveform_FR_1, 'FaceColor', 'red');
bar(1, SS_FR(label_SS_waveform==1),'FaceColor', 'red' );
bar_1.FaceAlpha = 0.3;
bar_2 = bar(2, mean_SS_waveform_FR_2, 'FaceColor', 'green');
bar(2, SS_FR(label_SS_waveform==2),'FaceColor', 'green' );
bar_2.FaceAlpha = 0.3;
bar_3 = bar(3, mean_SS_waveform_FR_3, 'FaceColor', 'blue');
bar(3, SS_FR(label_SS_waveform==3),'FaceColor', 'blue' );
bar_3.FaceAlpha = 0.3;
bar_4 = bar(4, mean_SS_waveform_FR_4, 'FaceColor', 'cyan');
bar(4, SS_FR(label_SS_waveform==4),'FaceColor', 'cyan' );
bar_4.FaceAlpha = 0.3;
bar_5 = bar(5, mean_SS_waveform_FR_5, 'FaceColor', 'magenta');
bar(5, SS_FR(label_SS_waveform==5),'FaceColor', 'magenta' );
bar_5.FaceAlpha = 0.3;
xticks([1 2 3 4 5])
xlabel('Cluster #')
ylabel('Firing Rate (spks/s)')
ylim([0 90])
title('Mean Cluster Firing Rate')

pie_1 = subplot(3,4,4);
pie_values = [sum(label_SS_waveform==1) sum(label_SS_waveform==2) sum(label_SS_waveform==3) sum(label_SS_waveform==4) sum(label_SS_waveform==5)];
pie_labels = {num2str(sum(label_SS_waveform==1)) num2str(sum(label_SS_waveform==2)) num2str(sum(label_SS_waveform==3)) num2str(sum(label_SS_waveform==4)) num2str(sum(label_SS_waveform==5))};
pie(pie_values, pie_labels)
colormap(pie_1, [1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1]);
title('# Cells Per Cluster')

subplot(3,4,5)
xlabel_SSxSS = -19:20;
hold on
if length(se_SSxSS_1) > 1 
shade(xlabel_SSxSS, mean_SSxSS_1 + se_SSxSS_1, 'r', xlabel_SSxSS, mean_SSxSS_1 - se_SSxSS_1, 'r',  'LineWidth', 0.1, 'FillColor', 'red', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(xlabel_SSxSS, mean_SSxSS_1, 'r')
if length(se_SSxSS_2) > 1 
shade(xlabel_SSxSS, mean_SSxSS_2 + se_SSxSS_2, '-g', xlabel_SSxSS, mean_SSxSS_2 - se_SSxSS_2, '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(xlabel_SSxSS, mean_SSxSS_2, '-g')
if length(se_SSxSS_3) > 1 
shade(xlabel_SSxSS, mean_SSxSS_3 + se_SSxSS_3, 'b', xlabel_SSxSS, mean_SSxSS_3 - se_SSxSS_3, 'b',  'LineWidth', 0.1, 'FillColor', 'blue', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(xlabel_SSxSS, mean_SSxSS_3, 'b')
if length(se_SSxSS_4) > 1 
shade(xlabel_SSxSS, mean_SSxSS_4 + se_SSxSS_4, 'c', xlabel_SSxSS, mean_SSxSS_4 - se_SSxSS_4, 'c',  'LineWidth', 0.1, 'FillColor', 'cyan', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(xlabel_SSxSS, mean_SSxSS_4, 'c')
if length(se_SSxSS_5) > 1 
shade(xlabel_SSxSS, mean_SSxSS_5 + se_SSxSS_5, 'm', xlabel_SSxSS, mean_SSxSS_5 - se_SSxSS_5, 'm',  'LineWidth', 0.1, 'FillColor', 'magenta', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(xlabel_SSxSS, mean_SSxSS_5, 'm')
xlabel('Time (ms)')
ylabel('Normalized Probability')
title('Clustered SSxSS')

subplot(3,4,6)
hold on
plot(umap_SSxSS(label_SSxSS==1,1), umap_SSxSS(label_SSxSS==1,2), '.r', 'MarkerSize', 10)
plot(umap_SSxSS(label_SSxSS==2,1), umap_SSxSS(label_SSxSS==2,2), '.g', 'MarkerSize', 10)
plot(umap_SSxSS(label_SSxSS==3,1), umap_SSxSS(label_SSxSS==3,2), '.b', 'MarkerSize', 10)
plot(umap_SSxSS(label_SSxSS==4,1), umap_SSxSS(label_SSxSS==4,2), '.c', 'MarkerSize', 10)
plot(umap_SSxSS(label_SSxSS==5,1), umap_SSxSS(label_SSxSS==5,2), '.m', 'MarkerSize', 10)
xlabel('UMAP 1')
ylabel('UMAP 2')
title('Clustered UMAP: SSxSS')

subplot(3,4,7)
hold on
bar_1 = bar(1, mean_SSxSS_FR_1, 'FaceColor', 'red');
bar(1, SS_FR(label_SSxSS == 1), 'FaceColor', 'red' );
bar_1.FaceAlpha = 0.3;
bar_2 = bar(2, mean_SSxSS_FR_2, 'FaceColor', 'green');
bar(2, SS_FR(label_SSxSS == 2), 'FaceColor', 'green' );
bar_2.FaceAlpha = 0.3;
bar_3 = bar(3, mean_SSxSS_FR_3, 'FaceColor', 'blue');
bar(3, SS_FR(label_SSxSS == 3), 'FaceColor', 'blue' );
bar_3.FaceAlpha = 0.3;
bar_4 = bar(4, mean_SSxSS_FR_4, 'FaceColor', 'cyan');
bar(4, SS_FR(label_SSxSS == 4), 'FaceColor', 'cyan' );
bar_4.FaceAlpha = 0.3;
bar_5 = bar(5, mean_SSxSS_FR_5, 'FaceColor', 'magenta');
bar(5, SS_FR(label_SSxSS == 5), 'FaceColor', 'magenta' );
bar_5.FaceAlpha = 0.3;
xticks([1 2 3 4 5])
xlabel('Cluster #')
ylabel('Firing Rate (spks/s)')
ylim([0 90])

pie_2 = subplot(3,4,8);
pie_values = [sum(label_SSxSS==1) sum(label_SSxSS==2) sum(label_SSxSS==3) sum(label_SSxSS==4) sum(label_SSxSS==5)];
pie_labels = {num2str(sum(label_SSxSS==1)) num2str(sum(label_SSxSS==2)) num2str(sum(label_SSxSS==3)) num2str(sum(label_SSxSS==4)) num2str(sum(label_SSxSS==5))};
pie(pie_values, pie_labels)
colormap(pie_2, [1 0 0; 0 1 0; 0 0 1; 0 1 1 ; 1 0 1]);

subplot(3,4,9)
xlabel_SS_ISI = (0:0.001:0.099);
hold on
if length(mean_SS_ISI_1) > 1 
shade(xlabel_SS_ISI, (mean_SS_ISI_1 + se_SS_ISI_1), 'r', xlabel_SS_ISI, (mean_SS_ISI_1 - se_SS_ISI_1), 'r',  'LineWidth', 0.1, 'FillColor', 'red', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(xlabel_SS_ISI,(mean_SS_ISI_1), 'r')
if length(mean_SS_ISI_2) > 1 
shade(xlabel_SS_ISI, (mean_SS_ISI_2 + se_SS_ISI_2), '-g', xlabel_SS_ISI, (mean_SS_ISI_2 - se_SS_ISI_2), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(xlabel_SS_ISI, (mean_SS_ISI_2), '-g')
if length(mean_SS_ISI_3) > 1 
shade(xlabel_SS_ISI, (mean_SS_ISI_3 + se_SS_ISI_3), 'b', xlabel_SS_ISI, (mean_SS_ISI_3 - se_SS_ISI_3), 'b',  'LineWidth', 0.1, 'FillColor', 'blue', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(xlabel_SS_ISI, (mean_SS_ISI_3), 'b')
if length(mean_SS_ISI_4) > 1 
shade(xlabel_SS_ISI, (mean_SS_ISI_4 + se_SS_ISI_4), 'c', xlabel_SS_ISI, (mean_SS_ISI_4 - se_SS_ISI_4), 'c',  'LineWidth', 0.1, 'FillColor', 'cyan', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(xlabel_SS_ISI, (mean_SS_ISI_4), 'c')
if length(mean_SS_ISI_5) > 1 
shade(xlabel_SS_ISI, (mean_SS_ISI_5 + se_SS_ISI_5), 'm', xlabel_SS_ISI, (mean_SS_ISI_5 - se_SS_ISI_5), 'm',  'LineWidth', 0.1, 'FillColor', 'magenta', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(xlabel_SS_ISI, (mean_SS_ISI_5), 'm')
xlabel('Normalized ISI (s)')
ylabel('Probability')
title('Clustered ISI')

subplot(3,4,10)
hold on
plot(umap_SS_ISI(label_SS_ISI==1,1), umap_SS_ISI(label_SS_ISI==1,2), '.r', 'MarkerSize', 10)
plot(umap_SS_ISI(label_SS_ISI==2,1), umap_SS_ISI(label_SS_ISI==2,2), '.g', 'MarkerSize', 10)
plot(umap_SS_ISI(label_SS_ISI==3,1), umap_SS_ISI(label_SS_ISI==3,2), '.b', 'MarkerSize', 10)
plot(umap_SS_ISI(label_SS_ISI==4,1), umap_SS_ISI(label_SS_ISI==4,2), '.c', 'MarkerSize', 10)
plot(umap_SS_ISI(label_SS_ISI==5,1), umap_SS_ISI(label_SS_ISI==5,2), '.m', 'MarkerSize', 10)
xlabel('UMAP 1')
ylabel('UMAP 2')
title('Clustered UMAP: ISI')

subplot(3,4,11)
hold on
bar_1 = bar(1, mean_SS_ISI_FR_1, 'FaceColor', 'red');
bar(1, SS_FR(label_SS_ISI == 1), 'FaceColor', 'red' );
bar_1.FaceAlpha = 0.3;
bar_2 = bar(2, mean_SS_ISI_FR_2, 'FaceColor', 'green');
bar(2, SS_FR(label_SS_ISI == 2), 'FaceColor', 'green' );
bar_2.FaceAlpha = 0.3;
bar_3 = bar(3, mean_SS_ISI_FR_3, 'FaceColor', 'blue');
bar(3, SS_FR(label_SS_ISI ==3), 'FaceColor', 'blue' );
bar_3.FaceAlpha = 0.3;
bar_4 = bar(4, mean_SS_ISI_FR_4, 'FaceColor', 'cyan');
bar(4, SS_FR(label_SS_ISI ==4), 'FaceColor', 'cyan' );
bar_4.FaceAlpha = 0.3;
bar_5 = bar(5, mean_SS_ISI_FR_5, 'FaceColor', 'magenta');
bar(5, SS_FR(label_SS_ISI ==5), 'FaceColor', 'magenta' );
bar_5.FaceAlpha = 0.3;
xticks([1 2 3 4 5])
xlabel('Cluster #')
ylabel('Firing Rate (spks/s)')
ylim([0 90])

pie_3 = subplot(3,4,12);
pie_values = [sum(label_SS_ISI==1) sum(label_SS_ISI==2) sum(label_SS_ISI==3) sum(label_SS_ISI==4) sum(label_SS_ISI==5)];
pie_labels = {num2str(sum(label_SS_ISI==1)) num2str(sum(label_SS_ISI==2)) num2str(sum(label_SS_ISI==3)) num2str(sum(label_SS_ISI==4)) num2str(sum(label_SS_ISI==5))};
pie(pie_values, pie_labels)
colormap(pie_3, [1 0 0; 0 1 0; 0 0 1; 0 1 1 ; 1 0 1]);
ESN_Beautify_Plot

%% Plot: clustered data - combined features 
fig_handle = figure;
fig_handle.WindowState = 'maximized';
sgtitle('Kmeans Clustering: Waveforms, SSxSS, and ISI')
% xlabel_waveform = -1.9:0.1:2;
% xlabel_SSxSS = -19:20;
% xlabel_SS_ISI = (0:0.001:0.099);
% xlabel_combined = [xlabel_waveform xlabel_SSxSS xlabel_SS_ISI];
xlabel_combined = 1:180;
subplot(3,2,1)
hold on
if length(se_combined_features_1) > 1 
shade(xlabel_combined, mean_combined_features_1 + se_combined_features_1, 'r', xlabel_combined, mean_combined_features_1 - se_combined_features_1, 'r',  'LineWidth', 0.1, 'FillColor', 'red', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(xlabel_combined, mean_combined_features_1, 'r')
if length(se_combined_features_2) > 1 
shade(xlabel_combined, mean_combined_features_2 + se_combined_features_2, '-g', xlabel_combined, mean_combined_features_2 - se_combined_features_2, '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(xlabel_combined, mean_combined_features_2, '-g')
if length(se_combined_features_3) > 1 
shade(xlabel_combined, mean_combined_features_3 + se_combined_features_3, 'b', xlabel_combined, mean_combined_features_3 - se_combined_features_3, 'b',  'LineWidth', 0.1, 'FillColor', 'blue', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(xlabel_combined, mean_combined_features_3, 'b')
if length(se_combined_features_4) > 1 
shade(xlabel_combined, mean_combined_features_4 + se_combined_features_4, 'c', xlabel_combined, mean_combined_features_4 - se_combined_features_4, 'c',  'LineWidth', 0.1, 'FillColor', 'cyan', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(xlabel_combined, mean_combined_features_4, 'c')
if length(se_combined_features_5) > 1 
shade(xlabel_combined, mean_combined_features_5 + se_combined_features_5, 'm', xlabel_combined, mean_combined_features_5 - se_combined_features_5, 'm',  'LineWidth', 0.1, 'FillColor', 'magenta', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
end
plot(xlabel_combined, mean_combined_features_5, 'm')
xlim([0 180])
xlabel('Sample #')
ylabel('Normalized unit')
title('Clustered Combined Features')

subplot(3,2,2)
hold on
plot(umap_combined_features(label_combined_features==1,1), umap_combined_features(label_combined_features==1,2), '.r', 'MarkerSize', 10)
plot(umap_combined_features(label_combined_features==2,1), umap_combined_features(label_combined_features==2,2), '.g', 'MarkerSize', 10 )
plot(umap_combined_features(label_combined_features==3,1), umap_combined_features(label_combined_features==3,2), '.b', 'MarkerSize', 10)
plot(umap_combined_features(label_combined_features==4,1), umap_combined_features(label_combined_features==4,2), '.c', 'MarkerSize', 10)
plot(umap_combined_features(label_combined_features==5,1), umap_combined_features(label_combined_features==5,2), '.m', 'MarkerSize', 10)
xlabel('UMAP 1')
ylabel('UMAP 2')
title('Clustered UMAP: Combined Features')

subplot(3,2,3)
hold on
bar_1 = bar(1, mean_combined_features_FR_1, 'FaceColor', 'red');
bar(1, SS_FR(label_combined_features==1),'FaceColor', 'red' );
bar_1.FaceAlpha = 0.3;
bar_2 = bar(2, mean_combined_features_FR_2, 'FaceColor', 'green');
bar(2, SS_FR(label_combined_features==2),'FaceColor', 'green' );
bar_2.FaceAlpha = 0.3;
bar_3 = bar(3, mean_combined_features_FR_3, 'FaceColor', 'blue');
bar(3, SS_FR(label_combined_features==3),'FaceColor', 'blue' );
bar_3.FaceAlpha = 0.3;
bar_4 = bar(4, mean_combined_features_FR_4, 'FaceColor', 'cyan');
bar(4, SS_FR(label_combined_features==4),'FaceColor', 'cyan' );
bar_4.FaceAlpha = 0.3;
bar_5 = bar(5, mean_combined_features_FR_5, 'FaceColor', 'magenta');
bar(5, SS_FR(label_combined_features==5),'FaceColor', 'magenta' );
bar_5.FaceAlpha = 0.3;
xticks([1 2 3 4 5])
xlabel('Cluster #')
ylabel('Firing Rate (spks/s)')
ylim([0 90])
title('Cluster Firing Rate')

pie_1 = subplot(3,2,4);
pie_values = [sum(label_combined_features==1) sum(label_combined_features==2) sum(label_combined_features==3) sum(label_combined_features==4) sum(label_combined_features==5)];
pie_labels = {num2str(sum(label_combined_features==1)) num2str(sum(label_combined_features==2)) num2str(sum(label_combined_features==3)) num2str(sum(label_combined_features==4)) num2str(sum(label_combined_features==5))};
pie(pie_values, pie_labels)
colormap(pie_1, [1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1]);
title('# Cells Per Cluster')

subplot(3,2,5)
hold on
bar_1 = bar(1, mean_combined_features_FR_1, 'FaceColor', 'red');
bar(1, SS_FR(label_combined_features==1),'FaceColor', 'red' );
bar_1.FaceAlpha = 0.3;
bar_2 = bar(2, mean_combined_features_FR_2, 'FaceColor', 'green');
bar(2, SS_FR(label_combined_features==2),'FaceColor', 'green' );
bar_2.FaceAlpha = 0.3;
bar_3 = bar(3, mean_combined_features_FR_3, 'FaceColor', 'blue');
bar(3, SS_FR(label_combined_features==3),'FaceColor', 'blue' );
bar_3.FaceAlpha = 0.3;
bar_4 = bar(4, mean_combined_features_FR_4, 'FaceColor', 'cyan');
bar(4, SS_FR(label_combined_features==4),'FaceColor', 'cyan' );
bar_4.FaceAlpha = 0.3;
bar_5 = bar(5, mean_combined_features_FR_5, 'FaceColor', 'magenta');
bar(5, SS_FR(label_combined_features==5),'FaceColor', 'magenta' );
bar_5.FaceAlpha = 0.3;
xticks([1 2 3 4 5])
xlabel('Cluster #')
ylabel('Duration (ms)')
title('Cluster Waveform Duration')
ESN_Beautify_Plot

%% ANALYSIS OF NEURAL MODULATION: SACCADES
% clearvars -except ALL_CELL_COMPRESSED_DATA POPULATION
alignment_ = ["cue_present" "primSac_onset" "prim_Sac_offset" "corrSac_onset" "corrSac_offset"];
alignment = alignment_(2);
direction = ["000" "045" "090" "135" "180" "225" "270" "315"];

% Load data
num_DCN_cell = length(ALL_CELL_COMPRESSED_DATA);
time = -299:300;

for counter_Cell = 1 : num_DCN_cell
    for counter_direction = 1:length(direction)
        SS_FR(counter_Cell,1) = ALL_CELL_COMPRESSED_DATA(counter_Cell).Neural_Properties_data.SS_firing_rate;  
        eval(['mean_SS_' convertStringsToChars(direction(counter_direction)) '(counter_Cell,:) = ESN_smooth(nanmean(ALL_CELL_COMPRESSED_DATA(counter_Cell).raster_data_' convertStringsToChars(alignment) '.train_data_logic_SS_' convertStringsToChars(direction(counter_direction)) ')*1000, 10);']);
        eval(['se_SS_' convertStringsToChars(direction(counter_direction)) '(counter_Cell,:) = ESN_smooth((nanstd(ALL_CELL_COMPRESSED_DATA(counter_Cell).raster_data_' convertStringsToChars(alignment) '.train_data_logic_SS_' convertStringsToChars(direction(counter_direction)) ')) / sqrt(size(ALL_CELL_COMPRESSED_DATA(counter_Cell).raster_data_' convertStringsToChars(alignment) '.train_data_logic_SS_' convertStringsToChars(direction(counter_direction)) ',2))*1000, 10);']);
        eval(['change_mean_SS_' convertStringsToChars(direction(counter_direction)) '(counter_Cell,:) = mean_SS_' convertStringsToChars(direction(counter_direction)) '(counter_Cell,:) - nanmean(mean_SS_' convertStringsToChars(direction(counter_direction)) '(counter_Cell,1:200));']);
        eval(['mean_velocity_' convertStringsToChars(direction(counter_direction)) '(counter_Cell,:) = nanmean(ALL_CELL_COMPRESSED_DATA(counter_Cell).raster_data_' convertStringsToChars(alignment) '.velocity_data_' convertStringsToChars(direction(counter_direction)) ');']);
        eval(['se_velocity_' convertStringsToChars(direction(counter_direction)) '(counter_Cell,:) = (nanstd(ALL_CELL_COMPRESSED_DATA(counter_Cell).raster_data_' convertStringsToChars(alignment) '.velocity_data_' convertStringsToChars(direction(counter_direction)) ')) / sqrt(size(ALL_CELL_COMPRESSED_DATA(counter_Cell).raster_data_' convertStringsToChars(alignment) '.velocity_data_' convertStringsToChars(direction(counter_direction)) ',2));']);
    end
end

% Build tuning curve
for counter_Cell = 1 : num_DCN_cell
    for counter_direction = 1:length(direction)
        eval(['[max_change_mean_SS_' convertStringsToChars(direction(counter_direction)) '(counter_Cell, 1), ind_max_change_mean_SS_' convertStringsToChars(direction(counter_direction)) '(counter_Cell, 1)] = max(change_mean_SS_' convertStringsToChars(direction(counter_direction)) '(counter_Cell,250:350));']);
        eval(['[min_change_mean_SS_' convertStringsToChars(direction(counter_direction)) '(counter_Cell, 1), ind_min_change_mean_SS_' convertStringsToChars(direction(counter_direction)) '(counter_Cell, 1)] = min(change_mean_SS_' convertStringsToChars(direction(counter_direction)) '(counter_Cell,250:350));']);
        if eval(['max_change_mean_SS_' convertStringsToChars(direction(counter_direction)) '(counter_Cell, 1) > abs(min_change_mean_SS_'  convertStringsToChars(direction(counter_direction)) '(counter_Cell, 1));' ])
            eval(['tuning_curve_' convertStringsToChars(direction(counter_direction)) '(counter_Cell, 1) = max_change_mean_SS_' convertStringsToChars(direction(counter_direction)) '(counter_Cell, 1);']);
            eval(['ind_tuning_curve_' convertStringsToChars(direction(counter_direction)) '(counter_Cell, 1) = ind_max_change_mean_SS_' convertStringsToChars(direction(counter_direction)) '(counter_Cell, 1);']);
        else
            eval(['tuning_curve_' convertStringsToChars(direction(counter_direction)) '(counter_Cell, 1) = min_change_mean_SS_' convertStringsToChars(direction(counter_direction)) '(counter_Cell, 1);']);
            eval(['ind_tuning_curve_' convertStringsToChars(direction(counter_direction)) '(counter_Cell, 1) = ind_min_change_mean_SS_' convertStringsToChars(direction(counter_direction)) '(counter_Cell, 1);']) ;          
        end
    end
end
tuning_curve = [tuning_curve_000 tuning_curve_045 tuning_curve_090 tuning_curve_135 tuning_curve_180 tuning_curve_225 tuning_curve_270 tuning_curve_315];
ind_tuning_curve = [ind_tuning_curve_000 ind_tuning_curve_045 ind_tuning_curve_090 ind_tuning_curve_135 ind_tuning_curve_180 ind_tuning_curve_225 ind_tuning_curve_270 ind_tuning_curve_315];
time_tuning_curve = time(ind_tuning_curve) + 250;
% UMAP
% addpath('C:\Users\Paul\OneDrive - Johns Hopkins\Shadmehr Lab\Code\umapFileExchange\umap')
[umap_tuning_curve, ~, ~, ~] = run_umap(tuning_curve);
umap_tuning_curve = umap_tuning_curve(:,1:2);
close all;

% Kmeans cluster check
cluster_check = 0;
if cluster_check == 1
    for counter_k = 1:10
        [label_tuning_curve, centroid_tuning_curve, sum_distances_tuning_curve] = kmeans(umap_tuning_curve, counter_k);
        WCSS_tuning_curve(counter_k) = sum(sum_distances_tuning_curve);
    end  
    % Plot: WCSS x number of clusters
    figure
    plot(1:10, WCSS_tuning_curve, '-k', 'LineWidth', 2)
    xlabel('Number of clusters')
    ylabel('WCSS')
    title('Change in Firing Rate Tuning Curve: WCSS vs number of clusters')
end

%Kmeans with selected number of clusters
num_clusters = 3;
[label_tuning_curve, centroid_tuning_curve, sum_distances_tuning_curve] = kmeans(umap_tuning_curve, num_clusters);

% Compute mean of tuning_curve clusters
for counter_cluster = 1:num_clusters
    eval(['mean_tuning_curve_' num2str(counter_cluster) ' = nanmean(tuning_curve(label_tuning_curve == ' num2str(counter_cluster) ',:));'])
    eval(['se_tuning_curve_' num2str(counter_cluster) ' = nanstd(tuning_curve(label_tuning_curve == ' num2str(counter_cluster) ',:)) / sqrt(sum(label_tuning_curve == ' num2str(counter_cluster) '));'])
end

% Compute tuning curve mean of clusters
for counter_cluster = 1:num_clusters
    eval(['mean_tuning_curve_' num2str(counter_cluster) ' = nanmean(tuning_curve(label_tuning_curve == ' num2str(counter_cluster) ',:));'])
    eval(['se_tuning_curve_' num2str(counter_cluster) ' = nanstd(tuning_curve(label_tuning_curve == ' num2str(counter_cluster) ',:)) / sqrt(sum(label_tuning_curve == ' num2str(counter_cluster) '));'])
end

% Compute modulation mean of clusters
for counter_Cell = 1 : num_DCN_cell
    for counter_direction = 1:length(direction)
        for counter_cluster = 1:num_clusters
        eval(['change_mean_SS_' convertStringsToChars(direction(counter_direction)) '_' num2str(counter_cluster) '(counter_Cell,:) = nanmean(change_mean_SS_' convertStringsToChars(direction(counter_direction))  '(label_tuning_curve ==' num2str(counter_cluster) ',:));']);
        end
    end
end
%% Cell select
DCN_cell = [6];
%% Plot: data
fig_handle = figure;
fig_handle.WindowState = 'maximized';
set(0, 'DefaultTextInterpreter', 'none')
sgtitle(['Single cell modulation | Selected cell: ' num2str(DCN_cell)])
xaxis = -299:300;

subplot(3,3,1)
hold on
shade(xaxis, (mean_SS_135(DCN_cell,:) + se_SS_135(DCN_cell,:)), '-k', xaxis,  (mean_SS_135(DCN_cell,:) - se_SS_135(DCN_cell,:)), '-k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, (mean_SS_135(DCN_cell,:)), '-k')
xline(0, '-k');
xlabel(char(['Time from' alignment ' (ms)']))
ylabel('Firing rate (spks/s)')
ylim([0 300])
yyaxis right
shade(xaxis, mean_velocity_135(DCN_cell,:) + se_velocity_135(DCN_cell,:), '-g', xaxis,  mean_velocity_135(DCN_cell,:) - se_velocity_135(DCN_cell,:), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, mean_velocity_135(DCN_cell,:), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

subplot(3,3,2)
hold on
shade(xaxis, (mean_SS_090(DCN_cell,:) + se_SS_090(DCN_cell,:)), '-k', xaxis,  (mean_SS_090(DCN_cell,:) - se_SS_090(DCN_cell,:)), '-k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, (mean_SS_090(DCN_cell,:)), '-k')
xline(0, '-k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Firing rate (spks/s)')
ylim([0 300])
yyaxis right
shade(xaxis, mean_velocity_090(DCN_cell,:) + se_velocity_090(DCN_cell,:), '-g', xaxis,  mean_velocity_090(DCN_cell,:) - se_velocity_090(DCN_cell,:), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, mean_velocity_090(DCN_cell,:), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

subplot(3,3,3)
hold on
shade(xaxis, (mean_SS_045(DCN_cell,:) + se_SS_045(DCN_cell,:)), '-k', xaxis,  (mean_SS_045(DCN_cell,:) - se_SS_045(DCN_cell,:)), '-k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, (mean_SS_045(DCN_cell,:)), '-k')
xline(0, '-k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Firing rate (spks/s)')
ylim([0 300])
yyaxis right
shade(xaxis, mean_velocity_045(DCN_cell,:) + se_velocity_045(DCN_cell,:), '-g', xaxis,  mean_velocity_045(DCN_cell,:) - se_velocity_045(DCN_cell,:), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, mean_velocity_045(DCN_cell,:), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

subplot(3,3,4)
hold on
shade(xaxis, (mean_SS_180(DCN_cell,:) + se_SS_180(DCN_cell,:)), '-k', xaxis,  (mean_SS_180(DCN_cell,:) - se_SS_180(DCN_cell,:)), '-k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, (mean_SS_180(DCN_cell,:)), '-k')
xline(0, '-k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Firing rate (spks/s)')
ylim([0 300])
yyaxis right
shade(xaxis, mean_velocity_180(DCN_cell,:) + se_velocity_180(DCN_cell,:), '-g', xaxis,  mean_velocity_180(DCN_cell,:) - se_velocity_180(DCN_cell,:), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, mean_velocity_180(DCN_cell,:), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

subplot(3,3,5)
hold on
plot(tuning_curve(DCN_cell,:), 'k', 'LineWidth', 2)
xticklabels({'0', '45', '90', '135',  '180', '225', '270', '315'});
xticks([1 2 3 4 5 6 7 8])
yline(0, 'k', 'LineWidth', 1);
xlabel('Saccade direction (deg)')
ylabel('Change in firing rate (spks/s)')
yyaxis right 
plot(time_tuning_curve(DCN_cell,:),'c', 'LineWidth', 2)
ylabel('Time of peak modulation from alignment (ms)')
yline(0, 'c', 'LineWidth', 1);
title('Change in Firing Rate Tuning Curve')
set(gca, 'ycolor', 'c')

subplot(3,3,6)
hold on
shade(xaxis, (mean_SS_000(DCN_cell,:) + se_SS_000(DCN_cell,:)), '-k', xaxis,  (mean_SS_000(DCN_cell,:) - se_SS_000(DCN_cell,:)), '-k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, (mean_SS_000(DCN_cell,:)), '-k')
xline(0, '-k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Firing rate (spks/s)')
ylim([0 300])
yyaxis right
shade(xaxis, mean_velocity_000(DCN_cell,:) + se_velocity_000(DCN_cell,:), '-g', xaxis,  mean_velocity_000(DCN_cell,:) - se_velocity_000(DCN_cell,:), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, mean_velocity_000(DCN_cell,:), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

subplot(3,3,7)
hold on
shade(xaxis, (mean_SS_225(DCN_cell,:) + se_SS_225(DCN_cell,:)), '-k', xaxis,  (mean_SS_225(DCN_cell,:) - se_SS_225(DCN_cell,:)), '-k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, (mean_SS_225(DCN_cell,:)), '-k')
xline(0, '-k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Firing rate (spks/s)')
ylim([0 300])
yyaxis right
shade(xaxis, mean_velocity_225(DCN_cell,:) + se_velocity_225(DCN_cell,:), '-g', xaxis,  mean_velocity_225(DCN_cell,:) - se_velocity_225(DCN_cell,:), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, mean_velocity_225(DCN_cell,:), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

subplot(3,3,8)
hold on
shade(xaxis, (mean_SS_270(DCN_cell,:) + se_SS_270(DCN_cell,:)), '-k', xaxis,  (mean_SS_270(DCN_cell,:) - se_SS_270(DCN_cell,:)), '-k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis,( mean_SS_270(DCN_cell,:)), '-k')
xline(0, '-k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Firing rate (spks/s)')
ylim([0 300])
yyaxis right
shade(xaxis, mean_velocity_270(DCN_cell,:) + se_velocity_270(DCN_cell,:), '-g', xaxis,  mean_velocity_270(DCN_cell,:) - se_velocity_270(DCN_cell,:), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, mean_velocity_270(DCN_cell,:), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

subplot(3,3,9)
hold on
shade(xaxis, (mean_SS_315(DCN_cell,:) + se_SS_315(DCN_cell,:)), '-k', xaxis,  (mean_SS_315(DCN_cell,:) - se_SS_315(DCN_cell,:)), '-k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, (mean_SS_315(DCN_cell,:)), '-k')
xline(0, '-k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Firing rate (spks/s)')
ylim([0 300])
yyaxis right
shade(xaxis, mean_velocity_315(DCN_cell,:) + se_velocity_315(DCN_cell,:), '-g', xaxis,  mean_velocity_315(DCN_cell,:) - se_velocity_315(DCN_cell,:), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, mean_velocity_315(DCN_cell,:), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

%  ESN_Beautify_Plot

%% Plot: clustered data based on tuning curves
fig_handle = figure;
fig_handle.WindowState = 'maximized';
sgtitle('Clustering based on tuning curve')
xaxis = -299:300;

subplot(4,3,1)
hold on
plot(umap_tuning_curve(label_tuning_curve==1,1), umap_tuning_curve(label_tuning_curve==1,2), '.r', 'MarkerSize', 10)
plot(umap_tuning_curve(label_tuning_curve==2,1), umap_tuning_curve(label_tuning_curve==2,2), '.g', 'MarkerSize', 10 )
plot(umap_tuning_curve(label_tuning_curve==3,1), umap_tuning_curve(label_tuning_curve==3,2), '.b', 'MarkerSize', 10 )
xlabel('UMAP 1')
ylabel('UMAP 2')
title('Clustered UMAP: tuning curve')

subplot(4,3,2)
hold on
tuning_curve_1 = tuning_curve(label_tuning_curve == 1,:);
for counter_cluster = 1 : sum(label_tuning_curve == 1)
plot(tuning_curve_1(counter_cluster,:),'-r', 'LineWidth', 0.5)
end
tuning_curve_2 = tuning_curve(label_tuning_curve == 2,:);
for counter_cluster = 1 : sum(label_tuning_curve == 2)
plot(tuning_curve_2(counter_cluster,:),'-g', 'LineWidth', 0.5)
end
tuning_curve_3 = tuning_curve(label_tuning_curve == 3,:);
for counter_cluster = 1 : sum(label_tuning_curve == 3)
plot(tuning_curve_3(counter_cluster,:),'-b', 'LineWidth', 0.5)
end
xticklabels({'0', '45', '90', '135',  '180', '225', '270', '315'});
xticks([1 2 3 4 5 6 7 8])
yline(0, 'k', 'LineWidth', 1);
xlabel('Saccade direction (deg)')
ylabel('Change in firing rate (spks/s)')
title('Change in Firing Rate Tuning Curve')

pie_1 = subplot(4,3,3);
pie_values = [sum(label_tuning_curve==1) sum(label_tuning_curve==2) sum(label_tuning_curve==3)];
pie_labels = {num2str(sum(label_tuning_curve==1)) num2str(sum(label_tuning_curve==2)) num2str(sum(label_tuning_curve==3))};
pie(pie_values, pie_labels)
colormap(pie_1, [1 0 0; 0 1 0; 0 0 1]);
title('Number of Cells Per Cluster')

subplot(4,3,4)
hold on
plot(xaxis, (change_mean_SS_135_1), '-r', 'LineWidth', 2)
plot(xaxis, (change_mean_SS_135_2), '-g', 'LineWidth', 2)
plot(xaxis, (change_mean_SS_135_3), '-b', 'LineWidth', 2)
xline(0, '-k');
yline(0, '-k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Change in firing rate (spks/s)')
ylim([-20 100])
yyaxis right
% shade(xaxis, mean_velocity_135(DCN_cell,:) + se_velocity_135(DCN_cell,:), '-g', xaxis,  mean_velocity_135(DCN_cell,:) - se_velocity_135(DCN_cell,:), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
% plot(xaxis, mean_velocity_135(DCN_cell,:), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

subplot(4,3,5)
hold on
plot(xaxis, (change_mean_SS_090_1), '-r', 'LineWidth', 2)
plot(xaxis, (change_mean_SS_090_2), '-g', 'LineWidth', 2)
plot(xaxis, (change_mean_SS_090_3), '-b', 'LineWidth', 2)
xline(0, '-k');
yline(0, '-k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Change in firing rate (spks/s)')
ylim([-20 100])
yyaxis right
% shade(xaxis, mean_velocity_135(DCN_cell,:) + se_velocity_135(DCN_cell,:), '-g', xaxis,  mean_velocity_135(DCN_cell,:) - se_velocity_135(DCN_cell,:), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
% plot(xaxis, mean_velocity_135(DCN_cell,:), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

subplot(4,3,6)
hold on
plot(xaxis, (change_mean_SS_045_1), '-r', 'LineWidth', 2)
plot(xaxis, (change_mean_SS_045_2), '-g', 'LineWidth', 2)
plot(xaxis, (change_mean_SS_045_3), '-b', 'LineWidth', 2)
xline(0, '-k');
yline(0, '-k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Change in firing rate (spks/s)')
ylim([-20 100])
yyaxis right
% shade(xaxis, mean_velocity_135(DCN_cell,:) + se_velocity_135(DCN_cell,:), '-g', xaxis,  mean_velocity_135(DCN_cell,:) - se_velocity_135(DCN_cell,:), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
% plot(xaxis, mean_velocity_135(DCN_cell,:), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

subplot(4,3,7)
hold on
plot(xaxis, (change_mean_SS_180_1), '-r', 'LineWidth', 2)
plot(xaxis, (change_mean_SS_180_2), '-g', 'LineWidth', 2)
plot(xaxis, (change_mean_SS_180_3), '-b', 'LineWidth', 2)
xline(0, '-k');
yline(0, '-k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Change in firing rate (spks/s)')
ylim([-20 100])
yyaxis right
% shade(xaxis, mean_velocity_135(DCN_cell,:) + se_velocity_135(DCN_cell,:), '-g', xaxis,  mean_velocity_135(DCN_cell,:) - se_velocity_135(DCN_cell,:), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
% plot(xaxis, mean_velocity_135(DCN_cell,:), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

subplot(4,3,8)
hold on
plot(mean_tuning_curve_1,'-r', 'LineWidth', 2)
plot(mean_tuning_curve_2,'-g', 'LineWidth', 2)
plot(mean_tuning_curve_3,'-b', 'LineWidth', 2)
xticklabels({'0', '45', '90', '135',  '180', '225', '270', '315'});
xticks([1 2 3 4 5 6 7 8])
yline(0, 'k', 'LineWidth', 1);
ylim([-50 100])
xlabel('Saccade direction (deg)')
ylabel('Change in firing rate (spks/s)')
title('Change in Firing Rate Tuning Curve')

subplot(4,3,9)
hold on
plot(xaxis, (change_mean_SS_000_1), '-r', 'LineWidth', 2)
plot(xaxis, (change_mean_SS_000_2), '-g', 'LineWidth', 2)
plot(xaxis, (change_mean_SS_000_3), '-b', 'LineWidth', 2)
xline(0, '-k');
yline(0, '-k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Change in firing rate (spks/s)')
ylim([-20 100])
yyaxis right
% shade(xaxis, mean_velocity_135(DCN_cell,:) + se_velocity_135(DCN_cell,:), '-g', xaxis,  mean_velocity_135(DCN_cell,:) - se_velocity_135(DCN_cell,:), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
% plot(xaxis, mean_velocity_135(DCN_cell,:), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

subplot(4,3,10)
hold on
plot(xaxis, (change_mean_SS_225_1), '-r', 'LineWidth', 2)
plot(xaxis, (change_mean_SS_225_2), '-g', 'LineWidth', 2)
plot(xaxis, (change_mean_SS_225_3), '-b', 'LineWidth', 2)
xline(0, '-k');
yline(0, '-k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Change in firing rate (spks/s)')
ylim([-20 100])
yyaxis right
% shade(xaxis, mean_velocity_135(DCN_cell,:) + se_velocity_135(DCN_cell,:), '-g', xaxis,  mean_velocity_135(DCN_cell,:) - se_velocity_135(DCN_cell,:), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
% plot(xaxis, mean_velocity_135(DCN_cell,:), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

subplot(4,3,11)
hold on
plot(xaxis, (change_mean_SS_270_1), '-r', 'LineWidth', 2)
plot(xaxis, (change_mean_SS_270_2), '-g', 'LineWidth', 2)
plot(xaxis, (change_mean_SS_270_3), '-b', 'LineWidth', 2)
xline(0, '-k');
yline(0, '-k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Change in firing rate (spks/s)')
ylim([-20 100])
yyaxis right
% shade(xaxis, mean_velocity_135(DCN_cell,:) + se_velocity_135(DCN_cell,:), '-g', xaxis,  mean_velocity_135(DCN_cell,:) - se_velocity_135(DCN_cell,:), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
% plot(xaxis, mean_velocity_135(DCN_cell,:), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

subplot(4,3,12)
hold on
plot(xaxis, (change_mean_SS_315_1), '-r', 'LineWidth', 2)
plot(xaxis, (change_mean_SS_315_2), '-g', 'LineWidth', 2)
plot(xaxis, (change_mean_SS_315_3), '-b', 'LineWidth', 2)
xline(0, '-k');
yline(0, '-k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Change in firing rate (spks/s)')
ylim([-20 100])
yyaxis right
% shade(xaxis, mean_velocity_135(DCN_cell,:) + se_velocity_135(DCN_cell,:), '-g', xaxis,  mean_velocity_135(DCN_cell,:) - se_velocity_135(DCN_cell,:), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
% plot(xaxis, mean_velocity_135(DCN_cell,:), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')
%% Plot: data averaged
fig_handle = figure;
fig_handle.WindowState = 'maximized';
sgtitle(['Population Response | ' num2str(num_DCN_cell) ' DCN cells'])

xaxis = -299:300;

subplot(3,3,1)
hold on
shade(xaxis, ESN_smooth(nanmean(change_mean_SS_135) + nanstd(change_mean_SS_135) / sqrt(size(change_mean_SS_135,1)),5), '-k', xaxis, ESN_smooth(nanmean(change_mean_SS_135) - nanstd(mean_SS_135) / sqrt(size(change_mean_SS_135,1)),5), '-k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, ESN_smooth(nanmean(change_mean_SS_135),5), '-k')
xline(0, '-k');yline(0,'k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Change in Firing rate (spks/s)')
ylim([-10 30])
yyaxis right
shade(xaxis, nanmean(mean_velocity_135) + nanstd(mean_velocity_135) / sqrt(size(mean_velocity_135,1)), '-g', xaxis, nanmean(mean_velocity_135) - nanstd(mean_velocity_135) / sqrt(size(mean_velocity_135,1)), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, nanmean(mean_velocity_135), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

subplot(3,3,2)
hold on
shade(xaxis, ESN_smooth(nanmean(change_mean_SS_090) + nanstd(change_mean_SS_090) / sqrt(size(change_mean_SS_090,1)),5), '-k', xaxis, ESN_smooth(nanmean(change_mean_SS_090) - nanstd(mean_SS_090) / sqrt(size(change_mean_SS_090,1)),5), '-k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, ESN_smooth(nanmean(change_mean_SS_090),5), '-k')
xline(0, '-k');yline(0,'k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Change in Firing rate (spks/s)')
ylim([-10 30])
yyaxis right
shade(xaxis, nanmean(mean_velocity_090) + nanstd(mean_velocity_090) / sqrt(size(mean_velocity_090,1)), '-g', xaxis, nanmean(mean_velocity_090) - nanstd(mean_velocity_090) / sqrt(size(mean_velocity_090,1)), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, nanmean(mean_velocity_090), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

subplot(3,3,3)
hold on
shade(xaxis, ESN_smooth(nanmean(change_mean_SS_045) + nanstd(change_mean_SS_045) / sqrt(size(change_mean_SS_045,1)),5), '-k', xaxis, ESN_smooth(nanmean(change_mean_SS_045) - nanstd(change_mean_SS_045) / sqrt(size(change_mean_SS_045,1)),5), '-k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, ESN_smooth(nanmean(change_mean_SS_045),5), '-k')
xline(0, '-k');yline(0,'k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Change in Firing rate (spks/s)')
ylim([-10 30])
yyaxis right
shade(xaxis, nanmean(mean_velocity_045) + nanstd(mean_velocity_045) / sqrt(size(mean_velocity_045,1)), '-g', xaxis, nanmean(mean_velocity_045) - nanstd(mean_velocity_045) / sqrt(size(mean_velocity_045,1)), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, nanmean(mean_velocity_045), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

subplot(3,3,4)
hold on
shade(xaxis, ESN_smooth(nanmean(change_mean_SS_180) + nanstd(change_mean_SS_180) / sqrt(size(change_mean_SS_180,1)),5), '-k', xaxis, ESN_smooth(nanmean(change_mean_SS_180) - nanstd(change_mean_SS_180) / sqrt(size(change_mean_SS_180,1)),5), '-k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, ESN_smooth(nanmean(change_mean_SS_180),5), '-k')
xline(0, '-k');yline(0,'k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Change in Firing rate (spks/s)')
ylim([-10 30])
yyaxis right
shade(xaxis, nanmean(mean_velocity_180) + nanstd(mean_velocity_180) / sqrt(size(mean_velocity_180,1)), '-g', xaxis, nanmean(mean_velocity_180) - nanstd(mean_velocity_180) / sqrt(size(mean_velocity_180,1)), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, nanmean(mean_velocity_180), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

subplot(3,3,6)
hold on
shade(xaxis, ESN_smooth(nanmean(change_mean_SS_000) + nanstd(change_mean_SS_000) / sqrt(size(change_mean_SS_000,1)),5), '-k', xaxis, ESN_smooth(nanmean(change_mean_SS_000) - nanstd(change_mean_SS_000) / sqrt(size(change_mean_SS_000,1)),5), '-k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, ESN_smooth(nanmean(change_mean_SS_000),5), '-k')
xline(0, '-k');yline(0,'k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Change in Firing rate (spks/s)')
ylim([-10 30])
yyaxis right
shade(xaxis, nanmean(mean_velocity_000) + nanstd(mean_velocity_000) / sqrt(size(mean_velocity_000,1)), '-g', xaxis, nanmean(mean_velocity_000) - nanstd(mean_velocity_000) / sqrt(size(mean_velocity_000,1)), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, nanmean(mean_velocity_000), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

subplot(3,3,5)
hold on
plot(nanmean(tuning_curve), 'k', 'LineWidth', 2)
xticklabels({'0', '45', '90', '135',  '180', '225', '270', '315'});
xticks([1 2 3 4 5 6 7 8])
yline(0, 'k', 'LineWidth', 1);
xlabel('Saccade direction (deg)')
ylabel('Change in firing rate (spks/s)')
title('Change in Firing Rate Tuning Curve')

subplot(3,3,7)
hold on
shade(xaxis, ESN_smooth(nanmean(change_mean_SS_225) + nanstd(change_mean_SS_225) / sqrt(size(change_mean_SS_225,1)),5), '-k', xaxis, ESN_smooth(nanmean(change_mean_SS_225) - nanstd(change_mean_SS_225) / sqrt(size(change_mean_SS_225,1)),5), '-k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, ESN_smooth(nanmean(change_mean_SS_225),5), '-k')
xline(0, '-k');yline(0,'k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Change in Firing rate (spks/s)')
ylim([-10 30])
yyaxis right
shade(xaxis, nanmean(mean_velocity_225) + nanstd(mean_velocity_225) / sqrt(size(mean_velocity_225,1)), '-g', xaxis, nanmean(mean_velocity_225) - nanstd(mean_velocity_225) / sqrt(size(mean_velocity_225,1)), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, nanmean(mean_velocity_225), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

subplot(3,3,8)
hold on
shade(xaxis, ESN_smooth(nanmean(change_mean_SS_270) + nanstd(change_mean_SS_270) / sqrt(size(change_mean_SS_270,1)),5), '-k', xaxis, ESN_smooth(nanmean(change_mean_SS_270) - nanstd(change_mean_SS_270) / sqrt(size(change_mean_SS_270,1)),5), '-k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, ESN_smooth(nanmean(change_mean_SS_270),5), '-k')
xline(0, '-k');yline(0,'k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Change in Firing rate (spks/s)')
ylim([-10 30])
yyaxis right
shade(xaxis, nanmean(mean_velocity_270) + nanstd(mean_velocity_270) / sqrt(size(mean_velocity_270,1)), '-g', xaxis, nanmean(mean_velocity_270) - nanstd(mean_velocity_270) / sqrt(size(mean_velocity_270,1)), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, nanmean(mean_velocity_270), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')

subplot(3,3,9)
hold on
shade(xaxis, ESN_smooth(nanmean(change_mean_SS_315) + nanstd(change_mean_SS_315) / sqrt(size(change_mean_SS_315,1)),5), '-k', xaxis, ESN_smooth(nanmean(change_mean_SS_315) - nanstd(change_mean_SS_315) / sqrt(size(change_mean_SS_315,1)),5), '-k',  'LineWidth', 0.1, 'FillColor', 'black', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, ESN_smooth(nanmean(change_mean_SS_315),5), '-k')
xline(0, '-k');yline(0,'k');
xlabel(char(['Time from' alignment ' (ms)'])) 
ylabel('Change in Firing rate (spks/s)')
ylim([-10 30])
yyaxis right
shade(xaxis, nanmean(mean_velocity_315) + nanstd(mean_velocity_315) / sqrt(size(mean_velocity_315,1)), '-g', xaxis, nanmean(mean_velocity_315) - nanstd(mean_velocity_315) / sqrt(size(mean_velocity_315,1)), '-g',  'LineWidth', 0.1, 'FillColor', 'green', 'Filltype', [1 2; 2 1], 'FillAlpha' , 0.1)
plot(xaxis, nanmean(mean_velocity_315), '-g')
ylabel('Velocity (deg/s)')
ylim([0 700])
set(gca, 'ycolor', 'g')


