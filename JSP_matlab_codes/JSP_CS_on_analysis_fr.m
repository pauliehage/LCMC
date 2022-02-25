function [CS_on_data] = JSP_CS_on_analysis_fr(SACS_ALL_DATA, params, funcs, tags_CS_ang_avg)

    if nargin == 3
        tags_CS_ang_avg = params.sac.tags_CS_ang_avg;
    end

    ang_edges       = params.sac.ang_edges;
    ang_values      = params.sac.ang_values;
    tag_name_list   = params.sac.tag_name_list;

    tag_bin = 1 : length(tag_name_list);
    % compute CS-on
    visual_px_offset = SACS_ALL_DATA.visual_px_offset;
    visual_py_offset = SACS_ALL_DATA.visual_py_offset;
    eye_r_px_onset   = SACS_ALL_DATA.eye_r_px_onset;
    eye_r_py_onset   = SACS_ALL_DATA.eye_r_py_onset;
    delta_x = visual_px_offset - eye_r_px_onset;
    delta_y = visual_py_offset - eye_r_py_onset;
    visual_ang = wrapTo360(atan2d(delta_y, delta_x));
%     visual_ang = wrapTo360(SACS_ALL_DATA.eye_r_ang);
    neuro_CS_count        = SACS_ALL_DATA.neuro_CS_count_visual;
    neuro_CS_count_period = SACS_ALL_DATA.neuro_CS_count_visual_period;
%     if (counter_pCell >= range_cell_with_4dir_behave(1)) && (counter_pCell <= range_cell_with_4dir_behave(2))
%         visual_ang_bin = discretize(ESN_Round(visual_ang, 90.0, 'round'), ang_edges);
%     else
%         
%     end
    visual_ang_bin = discretize(visual_ang, ang_edges);
    last_bin_id = length(ang_edges) - 1;
    visual_ang_bin(visual_ang_bin == last_bin_id) = 1; % wrap the circle around
    % 1: 0deg % 2: 45deg % 3: 90deg % 4: 135deg % 5: 180deg % 6: 225deg % 7: 270deg % 8: 315deg
    num_ang_bin = length(ang_edges) - 2;
    num_tag_bin = length(tag_bin);
    CS_count  = zeros(num_tag_bin, num_ang_bin);
    sac_count = zeros(num_tag_bin, num_ang_bin);
    CS_fr     = zeros(num_tag_bin, num_ang_bin);
    idx_validity = SACS_ALL_DATA.validity; 
    for counter_tag = 1 : num_tag_bin
        for counter_ang = 1 : num_ang_bin
            idx_tag      = (SACS_ALL_DATA.tag == tag_bin(counter_tag));
            idx_ang      = (visual_ang_bin == counter_ang);
%             flag_last_cue = SACS_ALL_DATA.flag_last_cue;
            idx_         = idx_tag & idx_ang & idx_validity;%&flag_last_cue;
            CS_count(counter_tag, counter_ang) = nansum(neuro_CS_count(1, idx_));
            sac_count(counter_tag, counter_ang) = nansum(idx_);            
            CS_fr(counter_tag, counter_ang) = nansum(neuro_CS_count(1, idx_)./neuro_CS_count_period(1,idx_))./sac_count(counter_tag, counter_ang).*1000; % avg. firing rate in Hz
        end
    end
    
    r = nansum(CS_fr.* repmat(exp(1i*deg2rad(ang_values)), num_tag_bin, 1) , 2); % compute weighted sum of cos and sin of angles
    CS_ang =  wrapTo360(rad2deg(angle(r))); % Computes the mean direction for circular data.
    CS_fr_sum = nansum(CS_fr,2); % sum of weights
    CS_rho = abs(r) ./ CS_fr_sum; % Computes mean resultant vector length for circular data.
    
    CS_fr_avg = zeros(size(CS_count( 1, :)));
    sac_count_avg = zeros(size(sac_count( 1, :)));

    for tag_ = tags_CS_ang_avg
        CS_fr_avg  = CS_fr_avg  + CS_fr( tag_, :).*sac_count(tag_,:);
        sac_count_avg = sac_count_avg + sac_count(tag_, :);
    end
    CS_fr_avg = CS_fr_avg ./ sac_count_avg;
    
    r_avg = nansum(CS_fr_avg.* exp(1i*deg2rad(ang_values)) , 2); % compute weighted sum of cos and sin of angles
    CS_ang_avg =  wrapTo360(rad2deg(angle(r_avg)));
    CS_rho_avg = abs(r_avg) ./ nansum(CS_fr_avg,2);
%     if (counter_pCell >= range_cell_with_4dir_behave(1)) && (counter_pCell <= range_cell_with_4dir_behave(2))
%         idx_CS_on_dir = discretize( ESN_Round(CS_ang_avg, 90.0, 'round'), ang_edges);
%     else
%         idx_CS_on_dir = discretize(CS_ang_avg, ang_edges);
%     end
    idx_CS_on_dir = discretize(CS_ang_avg, ang_edges);
    if idx_CS_on_dir == last_bin_id; idx_CS_on_dir = 1; end
    
    idx_ = idx_CS_on_dir - 1; % make it 0-index format
    if (idx_ == 8); idx_ = 0; end
    idx_CS_tuned = mod((idx_ : 1 : idx_+7), 8) + 1;
    CS_fr_avg_tuned = CS_fr_avg(idx_CS_tuned);
    
    % Von Mises
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
    vonMises_std = wrapTo360(rad2deg(sqrt(vonMises_var)));
    
    % Permutation for CS-on
    idx_tag_perm = false(size(SACS_ALL_DATA.tag));
    for tag_ = tags_CS_ang_avg
        idx_tag_perm = idx_tag_perm | (SACS_ALL_DATA.tag == tag_);
    end
    neuro_CS_count_perm = neuro_CS_count(idx_tag_perm);
    neuro_CS_count_period_perm = neuro_CS_count_period(idx_tag_perm);
    visual_ang_perm = visual_ang_bin(idx_tag_perm);
    idx_tag_perm = find(idx_tag_perm);
    
    num_perm = 1000;
    CS_ang_perm = nan(num_perm, 1);
    for counter_perm = 1 : num_perm
        idx_perm_ = randi(length(idx_tag_perm), 1, length(idx_tag_perm));
        CS_count_perm_ = neuro_CS_count_perm(idx_perm_);
        CS_count_period_perm = neuro_CS_count_period_perm(idx_perm_);
        visual_ang_perm_ = visual_ang_perm(idx_perm_);
        CS_fr_perm__   = zeros(1, num_ang_bin);
        for counter_ang = 1 : num_ang_bin
            idx_ang_ = (visual_ang_perm_ == counter_ang);
            CS_fr_perm__(1, counter_ang) = nansum(CS_count_perm_(idx_ang_)./CS_count_period_perm(idx_ang_))./nansum(idx_ang_).*1000;
        end
        r_perm_ = nansum(CS_fr_perm__.* exp(1i*deg2rad(ang_values)) , 2);
        CS_ang_perm(counter_perm, 1) =  wrapTo360(rad2deg(angle(r_perm_)));
    end
    CS_ang_avg_perm = mean(CS_ang_perm);
    CS_ang_std_perm = std(CS_ang_perm);
    
    % Build CS_on_data
    CS_on_data.sac_count = sac_count;
    CS_on_data.CS_count  = CS_count;
    CS_on_data.CS_fr = CS_fr;
    CS_on_data.CS_ang  = CS_ang;
    CS_on_data.CS_rho  = CS_rho;
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
    CS_on_data.visual_ang_bin  = visual_ang_bin;
    CS_on_data.CS_ang_avg_perm = CS_ang_avg_perm;
    CS_on_data.CS_ang_std_perm = CS_ang_std_perm;
end