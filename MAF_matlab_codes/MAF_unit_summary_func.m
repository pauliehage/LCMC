%% function units_summary
function unit_summary = MAF_unit_summary_func(SACS_ALL_DATA, LICKS_ALL_DATA, Neural_Properties, EXPERIMENT_PARAMS, params, funcs)

% type
unit_summary.type = Neural_Properties.type;

% waveform
unit_summary.ss.waveform.wave = Neural_Properties.waveform.ss_wave;
unit_summary.ss.waveform.ch   = Neural_Properties.waveform.ss_wave_ch;
unit_summary.cs.waveform.wave = Neural_Properties.waveform.cs_wave;
unit_summary.cs.waveform.ch   = Neural_Properties.waveform.cs_wave_ch;
unit_summary.ch_map           = Neural_Properties.waveform.ch_map;

% Cross/Auto-Corr
unit_summary.ss.corr          = Neural_Properties.Corr_data_SS_SSxSS_AUTO;
unit_summary.ss.corr(end/2)   = nan;
unit_summary.ss.corr_span     = Neural_Properties.Corr_data_SS_inds_span;
unit_summary.cs.corr          = Neural_Properties.Corr_data_CS_CSxSS_AUTO;
unit_summary.cs.corr_span     = Neural_Properties.Corr_data_CS_inds_span;

% isi
ss_isi = diff(Neural_Properties.SS_time);
ss_isi_min = 0;
ss_isi_max = .1;
ss_isi_binNum = 50;
ss_isi_edges = linspace(ss_isi_min, ss_isi_max, ss_isi_binNum);
unit_summary.ss.isi      = histcounts(ss_isi,ss_isi_edges);
unit_summary.ss.isi_edge = ss_isi_edges;

cs_isi = diff(Neural_Properties.CS_time);
cs_isi_min = 0;
cs_isi_max = 2;
cs_isi_binNum = 25;
cs_isi_edges = linspace(cs_isi_min, cs_isi_max, cs_isi_binNum);
unit_summary.cs.isi      = histcounts(cs_isi,cs_isi_edges);
unit_summary.cs.isi_edge = cs_isi_edges;

% Eye Modulation
% CS
tags_of_interest = params.lick.tags_CS_bundle_app;
alignment = params.lick.align_CS_bundle_app;
sac_data_dir = funcs.buildSacData(SACS_ALL_DATA,tags_of_interest,alignment,params,funcs);
mod_span = (0:200) + 250;
for counter_dir = 1:8
    unit_summary.cs.sac_raster{counter_dir} = ...
        logical(sac_data_dir(counter_dir).CS(mod_span,:));
end
unit_summary.cs.sac_raster_span = mod_span-250;

unit_summary.cs.sac_cs_on_data = funcs.sac_cs_on_analysis(SACS_ALL_DATA,params,funcs);

% SS
tags_of_interest = params.lick.tags_SS_bundle_app;
alignment = params.lick.align_SS_bundle_app;
sac_data_dir = funcs.buildSacData(SACS_ALL_DATA,tags_of_interest,alignment,params,funcs);
mod_span = (-100:100) + 250;
for counter_dir = 1:8
    unit_summary.ss.sac_raster{counter_dir} = ...
        logical(sac_data_dir(counter_dir).SS(mod_span,:));
end
unit_summary.ss.sac_raster_span = mod_span-250;

% Tongue Modulation
if not(isempty(LICKS_ALL_DATA))
    % CS
    tags_of_interest = params.lick.tags_CS_bundle_app;
    alignment = params.lick.align_CS_bundle_app;
    lick_data_dir = funcs.buildLickData(LICKS_ALL_DATA,tags_of_interest,alignment,params,funcs);
    mod_span = (-50:250) + 250;
    for counter_dir = 1:5
        unit_summary.cs.lick_raster{counter_dir} = ...
            logical(lick_data_dir(counter_dir).CS(mod_span,:));
    end
    unit_summary.cs.lick_raster_span = mod_span-250;
    unit_summary.cs.lick_cs_on_data = funcs.lick_cs_on_analysis(LICKS_ALL_DATA,params,funcs,tags_of_interest);
    
    % 6:8 bout
    if strcmp(params.lick.align_CS_bundle_app,'offset')
        alignment = 'offset';
    else
        alignment = 'onset';
    end
    bout_data_dir = buildBoutData(LICKS_ALL_DATA, alignment, params, funcs);
    % left, groom, right
    bout_order = [2,3,1];
    for counter_dir = 1:3
        unit_summary.cs.lick_raster{counter_dir + size(lick_data_dir,2)} = ...
            logical(bout_data_dir(bout_order(counter_dir)).CS(mod_span,:));
    end

    % SS
    tags_of_interest = params.lick.tags_SS_bundle_app;
    alignment = params.lick.align_SS_bundle_app;
    lick_data_dir = funcs.buildLickData(LICKS_ALL_DATA,tags_of_interest,alignment,params,funcs);
    mod_span = (-249:250) + 250;
    for counter_dir = 1:5
        unit_summary.ss.lick_raster{counter_dir} = ...
            logical(lick_data_dir(counter_dir).SS(mod_span,:));
    end
    unit_summary.ss.lick_raster_span = mod_span-250;

    % 6:8 bout
    if strcmp(params.lick.align_SS_bundle_app,'offset')
        alignment = 'offset';
    else
        alignment = 'onset';
    end
    bout_data_dir = buildBoutData(LICKS_ALL_DATA, alignment, params, funcs);
    % left, groom, right
    bout_order = [2,3,1];
    for counter_dir = 1:3
        unit_summary.ss.lick_raster{counter_dir + size(lick_data_dir,2)} = ...
            logical(bout_data_dir(bout_order(counter_dir)).SS(mod_span,:));
    end
    
end

end
