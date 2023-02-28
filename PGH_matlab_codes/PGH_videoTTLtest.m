clear all; close all; 
%% load 
dir_events = dir('*.events');
name_events = dir_events(1).name;
[ch_data, ch_time, ch_header] = load_open_ephys_data(name_events);
data = [ch_time ch_data ch_header.eventId(:)];

dir_vid = dir('*.avi');
for counter_vid = 1 : size(dir_vid,1)
    name_vid = dir_vid(counter_vid).name;
    vid{counter_vid,1} = VideoReader(name_vid);
end


%% plot
% initialize
ttl1 = zeros(1,length(ch_time));
ttl2 = zeros(1,length(ch_time));
ttl_time = data(:,1);

% find rising and falling edges
ind_ttl1_rise = find(data(:,2) == 5 & data(:,3) == 1);
ind_ttl1_fall = find(data(:,2) == 5 & data(:,3) == 0);
duration_ttl1 = ttl_time(ind_ttl1_rise(end)) - ttl_time(ind_ttl1_fall(1));

ind_ttl2_rise = find(data(:,2) == 6 & data(:,3) == 1);
ind_ttl2_fall = find(data(:,2) == 6 & data(:,3) == 0);
duration_ttl2 = ttl_time(ind_ttl2_rise(end)) - ttl_time(ind_ttl2_fall(1));


% build pulse train
for counter_pulse = 1 : min([length(ind_ttl1_rise) length(ind_ttl1_fall)])
    if ind_ttl1_fall(counter_pulse) - ind_ttl1_rise(counter_pulse) > 0
        ttl1(ind_ttl1_rise(counter_pulse):ind_ttl1_fall(counter_pulse)) = 1;
    else
        error('check rise/fall order')
    end
end

for counter_pulse = 1 : min([length(ind_ttl2_rise) length(ind_ttl2_fall)])
    if ind_ttl2_fall(counter_pulse) - ind_ttl2_rise(counter_pulse) > 0
        ttl2(ind_ttl2_rise(counter_pulse):ind_ttl2_fall(counter_pulse)) = 1;
    else
        error('check rise/fall order')
    end
end

figure
subplot(3,1,1)
hold on
plot(diff(ind_ttl1_rise),'g')
plot(diff(ind_ttl2_rise), 'r')
xlabel('sample #')
ylabel('diff_ttl_rise')

subplot(3,1,2)
hold on
plot(ttl1,'g')
plot(ttl2,'r')
xlabel('sample #')
ylabel('on/off')
ylim([-0.3 1.3])

subplot(3,1,3)
hold on
plot(ttl_time,ttl1,'g')
plot(ttl_time,ttl2,'r')
xlabel('time (s)')
ylabel('on/off')
ylim([-0.3 1.3])

sgtitle(['# frames(ttl vs vid): 1) ' num2str(max([length(ind_ttl1_rise) length(ind_ttl1_fall)])) ' vs ' num2str(vid{1, 1}.NumFrames) ', 2) ' num2str(max([length(ind_ttl2_rise) length(ind_ttl2_fall)]))  ' vs ' num2str(vid{2, 1}.NumFrames) newline ... 
    'duration(s; ttl vs vid): 1) ' num2str(duration_ttl1) ' vs '  num2str(vid{1, 1}.Duration)  ', 2) ' num2str(duration_ttl2)  ' vs '  num2str(vid{2, 1}.Duration)], 'Interpreter', 'none');

ESN_Beautify_Plot(gcf, [20 10])