%% params
params.sac.ang_step  = 45;
params.sac.length_trace = 500;
params.sac.trial_thresh = 100;
params.sac.inds_span        = ((-(params.sac.length_trace/2)+1) : 1 : (params.sac.length_trace/2))';
params.sac.ang_step         = 45;
params.sac.ang_edges        = (0 - (params.sac.ang_step/2)) : params.sac.ang_step : (360 + (params.sac.ang_step/2));
params.sac.ang_values       = (0) : params.sac.ang_step : (360 - params.sac.ang_step);
params.sac.tags_CS_ang_avg  = [1,4,6];
params.sac.tags_SS_ang_avg  = [1,4,6];
params.sac.tag_name_list    = { ...
    'prim_success', ... % tag 1
    'prim_attempt', ... % tag 2
    'prim_fail', ... % tag 3
    'corr_success', ... % tag 4
    'corr_fail', ... % tag 5
    'back_center_success', ... % tag 6
    'back_center_prim', ... % tag 7
    'back_center_irrelev', ... % tag 8
    'target_irrelev', ... % tag 9
    'other_irrelev', ... % tag 10
    'prim_no_corr',... % tag 11; prim. sac. that is not followed by corr. sac.
    'db_corr_success',... % tag 12; sac. that follows first corr. sac., back to 2nd jumped cue
    'corr_no_db_corr',... % tag 13; corr. sac. that is not followed by another corr. sac.
    'other_irrelev_visual',... % tag 14; like tag 10, but visual ang. based on the pursuit target present on the screen
    'back_center_irrelev_visual',... % tag 15; like tag 10, but visual start position based on offset of previous saccade
    ... % Add more tags here, do not reorder or change the tags defined above.
    };
params.sac.tags_CS_bundle_app = [1,4,6];
params.sac.tags_SS_bundle_app = [1,6,7];
params.sac.align_CS_bundle_app = 'visual';
params.sac.align_SS_bundle_app = 'onset';

params.lick.length_trace    = 2000;
params.lick.inds_span       = ((-(params.lick.length_trace/2)+1) : 1 : (params.lick.length_trace/2))';
params.lick.ang_step        = 45;
params.lick.ang_edges       = -90 - params.lick.ang_step /2 : params.lick.ang_step  : 90 + params.lick.ang_step /2 ;
params.lick.ang_values      = -90 : params.lick.ang_step : 90;
params.lick.amp_step    = 5;
params.lick.amp_edges    = 0 : params.lick.amp_step : 25;
params.lick.vel_step    = 150;
params.lick.vel_edges    = 0 : params.lick.vel_step : 750;
params.lick.dur_step = 50;
params.lick.dur_edges = 0 : params.lick.dur_step : 500;
params.lick.tags_CS_ang_avg  = [1:7];
params.lick.tag_name_list = {  ...
    'groom', ... % tag 1
    'inner_tube_success', ... % tag 2
    'inner_tube_fail', ... % tag 3
    'outer_edge_success', ... % tag 4
    'outer_edge_fail', ... % tag 5
    'under_tube_success', ... % tag 6
    'under_tube_fail', ... % tag 7
    };
params.lick.tag_bout_name_list = {  ...
    'bout_start', ... % tag 1
    'bout_end', ... % tag 2
    };
params.lick.tag_harvest_name_list = {  ...
    'harvest_start',  ...% tag 1
    'harvest_end', ... % tag 2
    };
params.lick.line_colors = [0,0,0; pink(round(1.5*length(params.lick.amp_edges)))];
params.lick.tags_CS_bundle_app = [1:7];
params.lick.tags_SS_bundle_app = [1:7];
params.lick.align_CS_bundle_app = 'dmax';
params.lick.align_SS_bundle_app = 'onset';
%% funcs
funcs.file_sorter           = @PGH_reorganize_files;
funcs.spike_sorter          = @KS;
funcs.monkey_behavior_sac   = @MAF_monkey_behavior_sac;
funcs.Sac_Sorter            = @JSP_Sac_Sorter;
funcs.Sac_Recal             = @JSP_recalibrate_all_sac;
funcs.add_ephys_sac_sorter  = @MAF_add_ephys_sac_sorter;
funcs.sac_cs_on_analysis    = @JSP_CS_on_analysis_fr;
funcs.sac_ss_on_analysis    = @MAF_SS_on_analysis_fr;
funcs.buildSacData          = @ESN_buildSacData;

funcs.monkey_behavior_lick  = @PGH_monkey_behavior_lick;
funcs.add_ephys_lick_sorter = @PGH_add_ephys_lick_sorter;
funcs.lick_cs_on_analysis   = @PGH_CS_on_analysis;
funcs.buildLickData         = @PGH_buildLickData;
funcs.buildBoutData         = @PGH_buildBoutData;

funcs.unit_summary_func     = @MAF_unit_summary_func;
