%% function PGH_CS_on_analysis
function CS_on_data = PGH_CS_on_analysis(LICKS_ALL_DATA,params,funcs,tags_of_interest)
if nargin<4
    tag_id = [1:9]; % default
elseif nargin==4
    tag_id = tags_of_interest;
end
ang_edges = params.lick.ang_edges;
ang_values = params.lick.ang_values;
tag_name_list = params.lick.tag_name_list;

tag_bin = 1 : length(tag_name_list);

%% compute CS-on
tongue_ang_bin = discretize(LICKS_ALL_DATA.tongue_ang_max, ang_edges);
neuro_CS_count = LICKS_ALL_DATA.neuro_CS_count_onset_offset;
neuro_CS_count_period = LICKS_ALL_DATA.neuro_CS_count_onset_offset_period;
% 1: -90deg % 2: -45deg % 3: 0deg % 4: 45deg % 5: 90deg 

num_ang_bin = length(ang_edges) - 1;
num_tag_bin = length(tag_bin);
CS_count  = zeros(num_tag_bin, num_ang_bin);
lick_count = zeros(num_tag_bin, num_ang_bin);
% CS_prob   = zeros(num_tag_bin, num_ang_bin);
CS_fr     = zeros(num_tag_bin, num_ang_bin);
for counter_tag = 1 : num_tag_bin
    for counter_ang = 1 : num_ang_bin
        idx_tag = (LICKS_ALL_DATA.tag == tag_bin(counter_tag));
        idx_ang = (tongue_ang_bin == counter_ang);
        idx_ = idx_tag&idx_ang;
        CS_count(counter_tag, counter_ang) = nansum(neuro_CS_count(1, idx_));
        lick_count(counter_tag, counter_ang) = nansum(idx_);
        %         CS_prob(counter_tag, counter_ang) = CS_count(counter_tag, counter_ang) ./ lick_count(counter_tag, counter_ang);
        CS_fr(counter_tag, counter_ang) = sum(neuro_CS_count(1, idx_)./neuro_CS_count_period(1,idx_))./lick_count(counter_tag, counter_ang).*1000; % avg. firing rate in Hz
    end
end

% r = nansum(CS_prob.* repmat(exp(1i*deg2rad(ang_values)), num_tag_bin, 1) , 2); % compute weighted sum of cos and sin of angles
r = nansum(CS_fr.* repmat(exp(1i*deg2rad(ang_values)), num_tag_bin, 1) , 2); % compute weighted sum of cos and sin of angles
CS_ang =  (rad2deg(angle(r))); % Computes the mean direction for circular data.
% CS_prob_sum = nansum(CS_prob,2); % sum of weights
% CS_rho = abs(r) ./ CS_prob_sum; % Computes mean resultant vector length for circular data.
CS_fr_sum = nansum(CS_fr,2); % sum of weights
CS_rho = abs(r) ./ CS_fr_sum; % Computes mean resultant vector length for circular data.

% CS_count_avg = zeros(size(CS_count( 1, :)));
CS_fr_avg = zeros(size(CS_count( 1, :)));
lick_count_avg = zeros(size(lick_count( 1, :)));

tags_CS_ang_avg = tag_id; % 'prim_success' tag 1 % 'corr_success' tag 4 % 'back_center_irrelev' tag 8

for tag_ = tags_CS_ang_avg
    %     CS_count_avg  = CS_count_avg  + CS_count( tag_, :);
    CS_fr_avg  = nansum([CS_fr_avg ; CS_fr( tag_, :).*lick_count(tag_,:)]);
    lick_count_avg = lick_count_avg + lick_count(tag_, :);
end
% CS_prob_avg = CS_count_avg ./ lick_count_avg;
    CS_fr_avg = CS_fr_avg ./ lick_count_avg;

% r_avg = nansum(CS_prob_avg.* exp(1i*deg2rad(ang_values)) , 2); % compute weighted sum of cos and sin of angles
    r_avg = nansum(CS_fr_avg.* exp(1i*deg2rad(ang_values)) , 2); % compute weighted sum of cos and sin of angles
CS_ang_avg =  (rad2deg(angle(r_avg)));

% CS_rho_avg = abs(r_avg) ./ nansum(CS_prob_avg,2);
    CS_rho_avg = abs(r_avg) ./ nansum(CS_fr_avg,2);

idx_CS_on_dir = discretize(CS_ang_avg, ang_edges);

idx_ = idx_CS_on_dir - 1; % make it 0-index format
idx_CS_tuned = mod((idx_ : 1 : idx_+4), 5) + 1;
% CS_prob_avg_tuned = CS_prob_avg(idx_CS_tuned);
    CS_fr_avg_tuned = CS_fr_avg(idx_CS_tuned);

%% Von Mises
if CS_rho_avg < 0.53
    vonMises_kappa = 2*CS_rho_avg + CS_rho_avg^3 + 5*CS_rho_avg^5/6;
elseif CS_rho_avg>=0.53 && CS_rho_avg<0.85
    vonMises_kappa = -.4 + 1.39*CS_rho_avg + 0.43/(1-CS_rho_avg);
else
    vonMises_kappa = 1/(CS_rho_avg^3 - 4*CS_rho_avg^2 + 3*CS_rho_avg);
end
% evaluate pdf
vonMises_C = 1/(2*pi*besseli(0,vonMises_kappa));
vonMises_pdf = vonMises_C * exp(vonMises_kappa*cosd(ang_values-CS_ang_avg));
vonMises_var = 1 - (besseli(1,vonMises_kappa) / besseli(0,vonMises_kappa));
vonMises_std = (rad2deg(sqrt(vonMises_var)));

%% Permutation for CS-on
idx_tag_perm = false(size(LICKS_ALL_DATA.tag));
for tag_ = tags_CS_ang_avg
    idx_tag_perm = idx_tag_perm | (LICKS_ALL_DATA.tag == tag_);
end
neuro_CS_count_perm = neuro_CS_count(idx_tag_perm);
neuro_CS_count_period_perm = neuro_CS_count_period(idx_tag_perm);
tongue_ang_perm = tongue_ang_bin(idx_tag_perm);
idx_tag_perm = find(idx_tag_perm);

num_perm = 1000;
CS_ang_perm = nan(num_perm, 1);
for counter_perm = 1 : num_perm
    idx_perm_ = randi(length(idx_tag_perm), 1, length(idx_tag_perm));
    CS_count_perm_ = neuro_CS_count_perm(idx_perm_);
   CS_count_period_perm = neuro_CS_count_period_perm(idx_perm_);
    tongue_ang_perm = tongue_ang_perm(idx_perm_);
%     CS_prob_perm__   = zeros(1, num_ang_bin);
        CS_fr_perm__   = zeros(1, num_ang_bin);
    for counter_ang = 1 : num_ang_bin
        idx_ang_ = (tongue_ang_perm == counter_ang);
                    CS_fr_perm__(1, counter_ang) = nansum(CS_count_perm_(idx_ang_)./CS_count_period_perm(idx_ang_))./nansum(idx_ang_).*1000;
%         CS_count_perm__ = nansum(CS_count_perm_(idx_ang_));
%         lick_count_perm__ = nansum(idx_ang_);
%         CS_prob_perm__(1, counter_ang) = CS_count_perm__ ./ lick_count_perm__;
    end
%     r_perm_ = nansum(CS_prob_perm__.* exp(1i*deg2rad(ang_values)) , 2);
        r_perm_ = nansum(CS_fr_perm__.* exp(1i*deg2rad(ang_values)) , 2);
   
CS_ang_perm(counter_perm, 1) =  (rad2deg(angle(r_perm_)));
end
CS_ang_avg_perm = mean(CS_ang_perm);
CS_ang_std_perm = std(CS_ang_perm);

%% Build CS_on_data
CS_on_data.lick_count = lick_count;
CS_on_data.CS_count  = CS_count;
% CS_on_data.CS_prob = CS_prob;
CS_on_data.CS_fr = CS_fr;
CS_on_data.CS_ang  = CS_ang;
CS_on_data.CS_rho  = CS_rho;
% CS_on_data.CS_prob_avg = CS_prob_avg;
% CS_on_data.CS_prob_avg_tuned = CS_prob_avg_tuned;
CS_on_data.CS_fr_avg = CS_fr_avg;
CS_on_data.CS_fr_avg_tuned = CS_fr_avg_tuned;
CS_on_data.CS_ang_avg  = CS_ang_avg;
CS_on_data.CS_rho_avg  = CS_rho_avg;
CS_on_data.vonMises_kappa  = vonMises_kappa;
CS_on_data.vonMises_var    = vonMises_var;
CS_on_data.vonMises_std    = vonMises_std;
CS_on_data.vonMises_pdf    = vonMises_pdf;
CS_on_data.idx_CS_on_dir   = idx_CS_on_dir;
CS_on_data.idx_CS_tuned    = idx_CS_tuned;
CS_on_data.tongue_ang_bin  = tongue_ang_bin;
CS_on_data.CS_ang_avg_perm = CS_ang_avg_perm;
CS_on_data.CS_ang_std_perm = CS_ang_std_perm;


end