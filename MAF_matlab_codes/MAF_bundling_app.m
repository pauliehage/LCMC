classdef MAF_bundling_app < handle
%MAF_BUNDLING_APP
%path_data_monkey_sorted, session, params, funcs
%ARGS:
%   path_data_monkey_sorted - string of the main directory of subject
%   session                 - session string in the format of YYYY-MM-DD
%   params                  - uses ang_values for lick and sac 
%   funcs                   - not used
%OUTPUTS:
%   rec meta data           - writes/modifies the rec-metadata with
%                             assigned type and bundles
%NOTE: 
%   *_units_summary.mat and session meta data for all recordings is
%   required. The app read from the rec-meta data (if saved before) and
%   overwrite it.
    properties
        elec                    % electrode properties
        recs                    % units_summary of recordings
        bndl                    % bundling result
        type                    % type result
        mpm                     % mpm values of recordings
        sess_name               % session name
        sess_path               % session path
        rec_path                % recs path
        sess_meta_data          % session metadata
        rec_meta_data           % recs metadata
    end
    properties (Access = private)
        colors                  % different colors for each rec
        mod_dd                  % SS/CS switch
        % Axes
        ax_wav                  % waveform plot
        ax_isi                  % isi plot
        ax_cor                  % corr plot
        mod_axes                % 1x8 cell of all mod axes
        ax_mod1                 % mod plot 0
        ax_mod2                 % mod plot 45
        ax_mod3                 % mod plot 90
        ax_mod4                 % mod plot 135
        ax_mod5                 % mod plot 180
        ax_mod6                 % mod plot 225
        ax_mod7                 % mod plot 270
        ax_mod8                 % mod plot 315
        panel_tune
        ax_tune
        % drop down
        dd                      % drop downs
        % check boxes
        chbx                    % checkbox 
        lbl_rec_name            % recording name label
        lbl_mpm                 % label mpm
        % radio buttons
        rb                      % radio buttons
        bg                      % button group for radio buttons
        % bundling labels
        bndl_lbl                % bundling label   
        lbl_max_lbl             % maximum bundling label
        % text inputs
        type_txa                % neuron type text input
        bndl_txa                % bundle text
        % buttons
        auto_btn                % auto button
        bndl_btn                % bundle button
        delt_btn                % delete button
        reset_btn               % reset button
        save_btn                % save button
        prev_btn                % reset button
        next_btn                % save button
        % slider thresh
        lbl_thresh              % threshold label
        sldr_thresh             % threshold slider
        thresh                  % threshold
        % slider KNN
        lbl_KNN                 % threshold label
        sldr_KNN                % threshold slider
        KNN                     % KNN
        % slider zoom
        lbl_zoom                % zoom label
        sldr_zoom               % zoom slider
        zm                      % zoom value
        % sac/lick radio buttons
        sac_lick_bg
        rb_lick
        rb_sac
        % params and funcs
        params
        funcs
    end
    methods
        function obj = MAF_bundling_app(path_data_monkey_sorted, session, params, funcs)

            obj.params = params;
            obj.funcs  = funcs;
            
            if ~strcmp(path_data_monkey_sorted(end), filesep)
                path_data_monkey_sorted = [path_data_monkey_sorted filesep];
            end
            
            obj.sess_path = [path_data_monkey_sorted, session(1:7), filesep,...
                session, filesep];
            
            obj.sess_name = regexprep(session(3:end),'-','');
            
            load_data(obj)
            
            n_recs = length(obj.recs);
            
            obj.colors = turbo(n_recs);
            obj.KNN = 3;
            obj.thresh = .85;
            obj.zm = 1;
            
            init_gui(obj);
            update_(obj);
        end

        function load_data(obj)
            obj.sess_meta_data = readtable([obj.sess_path obj.sess_name '.xls']);

            ind_rec_used = logical...
                (obj.sess_meta_data.ephys .* obj.sess_meta_data.eye);

            elec_id_ = obj.sess_meta_data.elec(ind_rec_used);
            elec_id = elec_id_(1);
            
            rec_list = obj.sess_meta_data.folder_name(ind_rec_used);

            obj.mpm = obj.sess_meta_data.MPM(ind_rec_used);

            n_recs = length(rec_list);

            obj.rec_path = cell(n_recs,1);

            flag_rec_meta_data = 1;
            % check rec_meta_data
            for counter_recs = 1:n_recs
                current_rec = rec_list{counter_recs};
                rec_name = regexprep(current_rec(3:end),'-','');
                obj.rec_path{counter_recs} = ...
                    [obj.sess_path current_rec filesep];
                % read rec_meta_data
                if not(isfile([obj.rec_path{counter_recs} rec_name '.xls']))
                    flag_rec_meta_data = 0;
                end
            end

            for counter_recs = 1:n_recs
                current_rec = rec_list{counter_recs};
                rec_name = regexprep(current_rec(3:end),'-','');
                obj.rec_path{counter_recs} = ...
                    [obj.sess_path current_rec filesep];
                obj.recs{counter_recs} = load([obj.rec_path{counter_recs} ...
                    rec_name '_units_summary.mat'],'units_summary');
                
                units_summary = obj.recs{counter_recs}.units_summary;

                n_units = length(units_summary);
                
                % read rec_meta_data
                if flag_rec_meta_data
                    obj.rec_meta_data = ...
                        readtable([obj.rec_path{counter_recs} rec_name '.xls']);
                    bndl_ = obj.rec_meta_data.cell;
                    bndl_(isnan(bndl_)) = 0;
                    if iscell(obj.rec_meta_data.type)
                        type_ = obj.rec_meta_data.type;
                    else
                        type_ = repmat({''},numel(obj.rec_meta_data.type),1);
                    end
                    obj.bndl{counter_recs} = bndl_;
                    obj.type{counter_recs} = type_;
    
                    if  not(n_units == length(bndl_))
                        q_ans = questdlg('Inconsistent Recording meta-data! Delete all recs meta data?',...
                            'Delete','Yes','No');
                        if strcmp(q_ans,'Yes')
                            for counter_recs2 = 1:n_recs
                                current_rec = rec_list{counter_recs2};
                                rec_name = regexprep(current_rec(3:end),'-','');
                                if isfile([obj.rec_path{counter_recs2} rec_name '.xls'])
                                    delete([obj.rec_path{counter_recs2} rec_name '.xls'])
                                end
                            end
                        end
                        error('Inconsistent Recording meta-data!');
                    end
                else
                    obj.bndl{counter_recs} = zeros(n_units,1);
                    obj.type{counter_recs} = ...
                        {obj.recs{counter_recs}.units_summary.type}';
                end
            end

            switch(elec_id)
                case(1)
                    obj.elec.x = 16*repmat([0,2,1,3],1,16);
                    obj.elec.y = 20*repelem(1:32,2);
                case(2)
                    obj.elec.x = zeros(1,64);
                    obj.elec.y = 20*(1:64);
                case(3)
                    obj.elec.x = 3*[sqrt(3),sqrt(3),sqrt(3)/2,3*sqrt(3)/2];
                    obj.elec.y = 6*[3,2,3/2,3/2];
                case(4)
                    obj.elec.x = 3*[0, sqrt(3), sqrt(3)/2, sqrt(3),...
                        3*sqrt(3)/2, sqrt(3), 2*sqrt(3)];
                    obj.elec.y = 6*[3, 0, 3/2, 2, 3/2, 3, 3];
            end
        end

        function init_gui(obj)
            
            set(0,'units','pixels')
            sz = get(0,'ScreenSize');

            x = sz(3);
            y = sz(4);
            br = 10;
            dist = 50;
            
            h = 30;
            fs = 14;
            w = 150;
            
            fig = uifigure('Position',[0,0,x,y],'Color','White');
            
            n_recs = numel(obj.recs);

            % Properties panel
            w_prop = round(x/4);
            
            % create Axes
            % wave axes
            l_wave = y/2;
            obj.ax_wav = uiaxes('Parent',fig,...
                'Position',[0 y-l_wave w_prop l_wave]+[1,1,-1,-1]*10*br);
            title(obj.ax_wav,'Waveform');
            set(obj.ax_wav, 'xcolor', 'w');
            set(obj.ax_wav, 'ycolor', 'w');
            set(obj.ax_wav, 'xtick', []);
            set(obj.ax_wav, 'ytick', []);
            MAF_Beautify_Ax(obj.ax_wav)

            % cor axes
            l_cor = y/3;
            obj.ax_cor = uiaxes('Parent',fig,...
                'Position',[0 y-l_cor-l_wave w_prop l_cor]+[1,1,-1,-1]*4*br);
            title(obj.ax_cor,'Auto Prob.');
            ylabel(obj.ax_cor,'prob.');
            xlabel(obj.ax_cor,'time (ms)');
            MAF_Beautify_Ax(obj.ax_cor)
            
            % isi axes
            l_isi = y/6;
            obj.ax_isi = uiaxes('Parent',fig,...
                'Position',[0 0 w_prop l_isi]+[1,1,-1,-1]*4*br);
            title(obj.ax_isi,'ISI dist');
            ylabel(obj.ax_isi,'prob.');
            xlabel(obj.ax_isi,'ISI (ms)');
            MAF_Beautify_Ax(obj.ax_isi);

            
            % mod axes
            w_mod = round(x/2);
            obj.mod_axes = {obj.ax_mod1, obj.ax_mod2, obj.ax_mod3,...
                obj.ax_mod4, obj.ax_mod5, obj.ax_mod6,...
                obj.ax_mod7, obj.ax_mod8};
            l_mod_ax = y/3;
            w_mod_ax = w_mod/3;
            x_mod_axes = [2*w_mod_ax, 2*w_mod_ax, w_mod_ax,...
                0, 0, 0,...
                w_mod_ax, 2*w_mod_ax];
            y_mod_axes = [l_mod_ax, 2*l_mod_ax, 2*l_mod_ax,...
                2*l_mod_ax, l_mod_ax, 0,...
                0,0,0];
            for counter_axes = 1:8
                obj.mod_axes{counter_axes} = uiaxes('Parent',fig,...
                    'Position',[w_prop+x_mod_axes(counter_axes)...
                    y_mod_axes(counter_axes) w_mod_ax l_mod_ax]+...
                    [1,1,-1,-1]*4*br);
                MAF_Beautify_Ax(obj.mod_axes{counter_axes});
                if ismember(counter_axes,[4,5,6])
                    ylabel(obj.mod_axes{counter_axes},'trials');
                end
                if ismember(counter_axes,[6,7,8])
                    xlabel(obj.mod_axes{counter_axes},'time (ms)');
                end
            end

            obj.panel_tune = uipanel(fig,'Position',[w_prop+w_mod_ax...
                    l_mod_ax w_mod_ax l_mod_ax]+...
                    [1,1,-1,-1]*4*br,...
                    'Title','CS rate following target onset',...
                    'BackgroundColor','white',...
                    'TitlePosition','centertop','BorderType','none');
            obj.ax_tune = polaraxes(obj.panel_tune);
            axp = obj.ax_tune;
            hold(axp,'on')
            MAF_Beautify_Ax(axp)
            
            % lick sac radio button
            obj.sac_lick_bg = uibuttongroup(fig,'Title','Modulation',...
                'SelectionChangedFcn', @(source,eventdata) update_(obj),'Position',...
                [round(3*x/4)+br, y-w/2,...
                w/2,w/2]);
            obj.rb_lick = uiradiobutton(obj.sac_lick_bg,'Position',[br/2,br/2,w/2,w/4-br],'Text','lick','Value',0);
            obj.rb_sac = uiradiobutton(obj.sac_lick_bg,'Position',[br/2,2*br,w/2,w/4-br],'Text','sac','Value',1);

            ind_rec_used = logical...
                (obj.sess_meta_data.ephys .* obj.sess_meta_data.eye);
            if sum(obj.sess_meta_data.tongue(ind_rec_used)) == 0
                obj.sac_lick_bg.Enable = 'off';
            end

            % recording selection drop downs and text inputs
            obj.dd = cell(n_recs,1);
            obj.bndl_txa = cell(n_recs,1);
            obj.type_txa = cell(n_recs,1);
            obj.chbx = cell(n_recs,1);
            obj.rb = cell(n_recs,1);
            
            obj.bg = uibuttongroup(fig,'BackgroundColor','w','BorderType','none',...
                'SelectionChangedFcn', @(source,eventdata) rb_changed(obj),'Position',...
                [round(3*x/4)+3*w/2+1.5*br, y-9.5*br-n_recs*(dist),...
                w/6,n_recs*(dist)]);

            MPM = arrayfun(@num2str, obj.mpm,'UniformOutput', 0);
            
            for counter_recs = 1:n_recs
                rec_ = obj.recs{counter_recs}.units_summary;
                names = {rec_.name};
                rec_name = names{1}(1:13);
                unit_names = cellfun(@(x) x(15:end), names,...
                    'UniformOutput',false);
                
                obj.chbx{counter_recs} = uicheckbox(fig,'Text',...
                    '▮','FontSize',fs,'Value',1,...
                    'FontColor',obj.colors(counter_recs,:));
                obj.chbx{counter_recs}.Position = [round(3*x/4)+br/2,...
                    y-br*10-counter_recs*dist, w/5, h];
                
                obj.lbl_rec_name = uilabel(fig,'Text',rec_name,...
                    'FontSize',fs);
                obj.lbl_rec_name.Position = [round(3*x/4)+4*br,...
                    y-br*10-counter_recs*dist, w, h];
                
                obj.dd{counter_recs} = uidropdown(fig,...
                            'Position',[round(3*x/4)+w+br/4,...
                            y-br*10-counter_recs*dist, w/2+br, h],...
                            'Items',unit_names,...
                            'ItemsData',1:numel(unit_names),...
                            'FontSize',fs);
                        
                obj.rb{counter_recs} = uiradiobutton(obj.bg,'Position',...
                    [br/2,(counter_recs-1)*dist+br/2,br,br],'Text','');
                
                obj.bndl_txa{counter_recs} = uieditfield(fig,'numeric',...
                    'Editable','on','Position',...
                    [round(3*x/4)+3*w/2+7*br/2, y-br*10-counter_recs*dist,...
                    w/4, h],'FontSize',fs);

                obj.type_txa{counter_recs} = uieditfield(fig,'text',...
                    'Editable','on','Position',...
                    [round(3*x/4)+2*w+br/2,...
                    y-br*10-counter_recs*dist,...
                    w/4, h],'FontSize',fs);

                obj.lbl_mpm = uilabel(fig,'Text',MPM{counter_recs},...
                    'FontSize',fs);
                obj.lbl_mpm.Position = [round(3*x/4)+2*w+w/4+br,...
                    y-br*10-counter_recs*dist, w, h];
             
            end
            
            uilabel(fig,'Text','Max Bundle: ',...
                'FontSize',fs,'Position',...
                [round(3*x/4)+br, y-br*10, w, h]);
            obj.bndl_lbl = max(cellfun(@max,obj.bndl));
            obj.lbl_max_lbl = uilabel(fig,'Text',num2str(obj.bndl_lbl),...
                'FontSize',fs);
            obj.lbl_max_lbl.Position = [round(3*x/4)+br+w, y-br*10,...
                w/5, h];
            
            for counter_recs = 1:n_recs
                obj.type_txa{counter_recs}.ValueChangedFcn = ...
                    @(textarea,event) type_textEntered(obj);
                obj.bndl_txa{counter_recs}.ValueChangedFcn = ...
                    @(textarea,event) bndl_textEntered(obj);
	            obj.dd{counter_recs}.ValueChangedFcn = ...
                    @(dd1,event) update_(obj);
                obj.chbx{counter_recs}.ValueChangedFcn = ...
                    @(cbx,event) update_(obj);
            end

            obj.mod_dd = uidropdown(fig,...
                            'Position',[round(3*x/4)+3*br+w+w/5,...
                            y-br*10, w/2, h],...
                            'Items',{'ss','cs'},...
                            'FontSize',fs,'ValueChangedFcn',...
                            @(dd2,event) update_(obj));
            
            obj.auto_btn = uibutton(fig, 'Text', 'Auto', 'Position',...
                [round(3*x/4)+3*br,y-br*6-(n_recs+3)*dist, 2*w/3, 2*h],...
                'FontSize',fs,'ButtonPushedFcn',...
                @(auto_btn,event) auto_bt_pressed_(obj));
            
            obj.bndl_btn = uibutton(fig, 'Text', 'Bundle', 'Position',...
                [round(3*x/4)+3.5*br+2*w/3,y-br*6-(n_recs+3)*dist, 2*w/3, 2*h],...
                'FontSize',fs,'ButtonPushedFcn',...
                @(bndl_btn,event) bndl_bt_pressed_(obj));

            obj.delt_btn = uibutton(fig, 'Text', 'Delete', 'Position',...
                [round(3*x/4)+4*br+4*w/3,y-br*6-(n_recs+3)*dist, 2*w/3, 2*h],...
                'FontSize',fs,'ButtonPushedFcn',...
                @(delt_btn,event) delt_bt_pressed_(obj));

            obj.prev_btn = uibutton(fig, 'Text', 'prev', 'Position',...
                [round(3*x/4)+3*br,y-br*6-(n_recs+3)*dist+2*h+br, w, h],...
                'FontSize',fs,'ButtonPushedFcn',...
                @(prev_btn,event) prev_bt_pressed_(obj));
            
            obj.next_btn = uibutton(fig, 'Text', 'next', 'Position',...
                [round(3*x/4)+4*br+w,y-br*6-(n_recs+3)*dist+2*h+br, w, h],...
                'FontSize',fs,'ButtonPushedFcn',...
                @(next_btn,event) next_bt_pressed_(obj));
            
            obj.lbl_thresh = uilabel(fig,'Text','Threshold: ','FontSize',fs);
            obj.lbl_thresh.Position = [round(3*x/4)+4*br, y-br*6-...
                (n_recs+2)*dist-2*h-2*br, w, h/2];
            
            obj.sldr_thresh = uislider(fig,'Position',...
                [round(3*x/4)+br+w,y-br*6-(n_recs+2)*dist-2*h-br, w, 3],...
                'Limits',[.5,1],'Value', obj.thresh, ...
                'ValueChangedFcn',...
                @(sld,event) update_thresh(obj,sld));

            obj.lbl_KNN = uilabel(fig,'Text','KNN: ','FontSize',fs);
            obj.lbl_KNN.Position = [round(3*x/4)+4*br,...
                y-br*6-(n_recs+3)*dist-2*h-2*br, w, h/2];
            
            obj.sldr_KNN = uislider(fig,'Position',...
                [round(3*x/4)+br+w,y-br*6-(n_recs+3)*dist-2*h-br, w, 3],...
                'Limits',[1,min(length(obj.elec.x),32)],'Value', obj.KNN,...
                'ValueChangedFcn',...
                @(sld,event) update_KNN(obj,sld));

            obj.lbl_zoom = uilabel(fig,'Text','Zoom: ','FontSize',fs);
            obj.lbl_zoom.Position = [round(3*x/4)+4*br,...
                y-br*8-(n_recs+3)*dist-3*h-2*br, w, h/2];
            
            obj.sldr_zoom = uislider(fig,'Position',...
                [round(3*x/4)+br+w,y-br*8-(n_recs+3)*dist-3*h-br, w, 3],...
                'Limits',[0,2],'Value', obj.zm,...
                'ValueChangedFcn',...
                @(sld,event) update_zoom(obj,sld));


            obj.reset_btn = uibutton(fig, 'Text', 'Reset', 'Position',...
                [round(3*x/4)+3*br,y-br*12-(n_recs+3)*dist-7*h-br, w, 2*h],...
                'FontSize',fs,'ButtonPushedFcn',...
                @(reset_btn,event) reset_bt_pressed_(obj));
            
            obj.save_btn = uibutton(fig, 'Text', 'Save', 'Position',...
                [round(3*x/4)+4*br+w,y-br*12-(n_recs+3)*dist-7*h-br, w, 2*h],...
                'FontSize',fs,'ButtonPushedFcn',...
                @(save_btn,event) save_bt_pressed_(obj));

        end
        
        function rb_changed(obj)
            n_recs = numel(obj.rb);

            for counter_rec = 1:n_recs
                if obj.rb{counter_rec}.Value
                    rec_selected = n_recs - counter_rec + 1;
                end
            end

            obj.chbx{rec_selected}.Value = 1;
            update_(obj);
        end

        function update_thresh(obj,sld)
            obj.thresh = sld.Value;
            update_(obj);
        end
        
        function update_KNN(obj,sld)
            roundedValue = round(sld.Value);
            obj.sldr_KNN.Value = roundedValue;
            obj.KNN = roundedValue;
            update_(obj);
        end

        function update_zoom(obj,sld)
            obj.sldr_zoom.Value = sld.Value;
            obj.zm = sld.Value;
            update_(obj);
        end

        function prev_bt_pressed_(obj)
            n_recs = numel(obj.rb);
            
            for counter_rec = 1:n_recs
                if obj.rb{counter_rec}.Value
                    rec_selected = n_recs - counter_rec + 1;
                end
            end
            
            current_unit = obj.dd{rec_selected}.Value;

            if current_unit == 1
                prev = length(obj.recs{rec_selected}.units_summary);
            else
                prev = current_unit - 1;
            end

            obj.dd{rec_selected}.Value = prev;

            new_type = obj.type{rec_selected}{obj.dd{rec_selected}.Value};
            if strcmp(new_type, 'CS')
                obj.mod_dd.Value = 'cs';
            else
                obj.mod_dd.Value = 'ss';
            end

            update_(obj);
            
        end

        function next_bt_pressed_(obj)
            n_recs = numel(obj.rb);
            
            for counter_rec = 1:n_recs
                if obj.rb{counter_rec}.Value
                    rec_selected = n_recs - counter_rec + 1;
                end
            end
            
            current_unit = obj.dd{rec_selected}.Value;

            if current_unit == length(obj.recs{rec_selected}.units_summary)
                next = 1;
            else
                next = current_unit + 1;
            end

            obj.dd{rec_selected}.Value = next;

            new_type = obj.type{rec_selected}{obj.dd{rec_selected}.Value};
            if strcmp(new_type, 'CS')
                obj.mod_dd.Value = 'cs';
            else
                obj.mod_dd.Value = 'ss';
            end

            update_(obj);
        end

        function save_bt_pressed_(obj)
            n_recs = numel(obj.recs);
            for counter_recs = 1:n_recs
                cell_idx_max = max(cellfun(@max,obj.bndl));
                rec_ = obj.recs{counter_recs}.units_summary;
                names = {rec_.name};
                rec_name = names{1}(1:13);
                unit_names = cellfun(@(x) x(15:end), names, 'UniformOutput',false);
                bndl_ = obj.bndl{counter_recs};
                bndl_(bndl_ == 0) = ...
                    (1:sum(bndl_ == 0))+cell_idx_max;
                type_ = obj.type{counter_recs};
                obj.bndl{counter_recs} = bndl_;
                obj.rec_meta_data = table(unit_names',bndl_,type_);
                obj.rec_meta_data.Properties.VariableNames = ...
                    {'unit_name','cell','type'};
                writetable(obj.rec_meta_data,...
                    [obj.rec_path{counter_recs} rec_name],...
                    'FileType','spreadsheet');
            end
            update_(obj);
        end
        
        function reset_bt_pressed_(obj)
            n_recs = numel(obj.recs);
            for counter_recs = 1:n_recs
                rec_ = obj.recs{counter_recs}.units_summary;
                names = {rec_.name};
                unit_names = cellfun(@(x) x(15:end), names, 'UniformOutput',false);
                obj.bndl{counter_recs} = zeros(numel(unit_names),1);
            end
            update_(obj);
            
        end
        
        function auto_bt_pressed_(obj)
        
            n_recs = numel(obj.rb);
            
            for counter_rec = 1:n_recs
                if obj.rb{counter_rec}.Value
                    rec_selected = n_recs - counter_rec + 1;
                end
            end

            mod = (obj.mod_dd.Value);
            
            unit_selected = obj.recs{rec_selected}.units_summary...
                (obj.dd{rec_selected}.Value);
            wave_selected = unit_selected.(mod).waveform.wave;
            ch_idx_selected  = unit_selected.(mod).waveform.ch+1;
            if isnan(ch_idx_selected)
                errordlg(['No' mod 'shape for this unit!','Error']);
                return;
            end
            wave_ch_selected = unit_selected.ch_map.map(ch_idx_selected);
            auto_prob_selected = unit_selected.(mod).corr;
            auto_prob_selected(end/2) = [];
            ch_selected = str2double(unit_selected.name(15:16));
        
            elec_mat = [obj.elec.x' obj.elec.y'];
            idx_knn = knnsearch(elec_mat,...
                [obj.elec.x(ch_selected) obj.elec.y(ch_selected)],...
                'K',obj.KNN,'Distance','euclidean');
            
            for counter_rec = 1:n_recs
                if counter_rec == rec_selected
                    continue;
                end
                rec_ = obj.recs{counter_rec}.units_summary;
                rec_chs = cell2mat(cellfun(@(x) str2double(x(15:16)),...
                    {rec_.name},'UniformOutput',false));
                rec_KNN = find(ismember(rec_chs,idx_knn));
                if numel(rec_KNN) < 1
                    obj.chbx{counter_rec}.Value = 0;
                    update_(obj);
                    continue;
                end
        
                score_max = 0;
                best_match_id = rec_KNN(1);
                
                for counter_KNN = 1:numel(rec_KNN)
                    
                    id_ = rec_KNN(counter_KNN);
                    wave_ = rec_(id_).(mod).waveform.wave;
                    if isnan(wave_)
                        continue;
                    end
                    wave_ch_ = rec_(id_).ch_map.map...
                        (rec_(id_).(mod).waveform.ch+1);
        
                    [~,idx_ch_selected,idx_ch_knn] = ...
                        intersect(wave_ch_selected,wave_ch_);
        
                    auto_prob_ = rec_(id_).(mod).corr;
                    auto_prob_(50) = [];
                    
                    wv_R = corr2(wave_selected(idx_ch_selected,:),...
                        wave_(idx_ch_knn,:));
                    ap_R = corrcoef(auto_prob_selected,auto_prob_);
                    
                    wv_cc = abs(wv_R);
                    ap_cc = abs(ap_R(1,2));
                    
                    if strcmp(mod,'cs')
                        score = wv_cc;
                    else
                        score = (wv_cc+ap_cc)/2;
                    end
                    if score >= score_max
                        best_match_id = id_;
                        score_max = score;
                    end
                    
                end
        
                obj.dd{counter_rec}.Value = best_match_id;
                if score_max < obj.thresh
                    obj.chbx{counter_rec}.Value = 0;
                else
                    obj.chbx{counter_rec}.Value = 1;
                end
                
            end
            
            update_(obj);
        
        end
        
        function bndl_bt_pressed_(obj)
            
            nxt_bndl_lbl = max(cellfun(@max,obj.bndl)) + 1;
            n_recs = numel(obj.recs);
            flag_selected = zeros(n_recs,1);
            counter_selected = 1;
            for counter_rec = 1:n_recs
                if obj.chbx{counter_rec}.Value
                    flag_selected(counter_rec) = 1;
                    type_selected{counter_selected} =...
                        obj.type_txa{counter_rec}.Value;
                    counter_selected = counter_selected + 1;
                end
            end
            ind_selected = find(flag_selected == 1);
            for counter_ind = 1:numel(ind_selected)
                current_ind = ind_selected(counter_ind);
                obj.type_txa{current_ind}.Value = type_selected{1};
                obj.bndl_txa{current_ind}.Value = nxt_bndl_lbl;
            end
            obj.lbl_max_lbl.Text = num2str(nxt_bndl_lbl);
            bndl_textEntered(obj)
            type_textEntered(obj)
        end

        function delt_bt_pressed_(obj)
            
            delete_bndl_lbl = -1;
            n_recs = numel(obj.recs);
            flag_selected = zeros(n_recs,1);
            counter_selected = 1;
            for counter_rec = 1:n_recs
                if obj.chbx{counter_rec}.Value
                    flag_selected(counter_rec) = 1;
                    counter_selected = counter_selected + 1;
                end
            end
            ind_selected = find(flag_selected == 1);
            for counter_ind = 1:numel(ind_selected)
                current_ind = ind_selected(counter_ind);
                obj.bndl_txa{current_ind}.Value = delete_bndl_lbl;
            end
            bndl_textEntered(obj);
            obj.bndl_lbl = max(cellfun(@max,obj.bndl));
            obj.lbl_max_lbl.Text = num2str(obj.bndl_lbl);
        end
        
        function bndl_textEntered(obj)
            for counter_rec = 1:numel(obj.recs)
                obj.bndl{counter_rec}(obj.dd{counter_rec}.Value) = ...
                    obj.bndl_txa{counter_rec}.Value;
            end
            obj.bndl_lbl = max(cellfun(@max,obj.bndl));
            obj.lbl_max_lbl.Text = num2str(obj.bndl_lbl);
        end
        
        function type_textEntered(obj)
            for counter_rec = 1:numel(obj.recs)
                obj.type{counter_rec}{obj.dd{counter_rec}.Value} = ...
                    obj.type_txa{counter_rec}.Value;
            end
        end

        function update_(obj)
        
            n_recs = numel(obj.recs);
            
            cla(obj.ax_wav);
            cla(obj.ax_isi);
            cla(obj.ax_cor);
            for counter_ax = 1:numel(obj.mod_axes)
                yyaxis(obj.mod_axes{counter_ax},'right');
                cla(obj.mod_axes{counter_ax});
                yyaxis(obj.mod_axes{counter_ax},'left');
                cla(obj.mod_axes{counter_ax});
            end
            cla(obj.ax_tune);

            if obj.rb_sac.Value
                behav = 'sac';
                obj.ax_tune.ThetaZeroLocation = 'right';
                obj.ax_tune.ThetaDir = 'counterclockwise';
                thetaticks(obj.ax_tune, 0:45:315);
                thetaticklabels(obj.ax_tune,{'0','','90','', '180','','270', ''});
                thetalim(obj.ax_tune,[0 360]);
            else
                behav = 'lick';
                obj.ax_tune.ThetaZeroLocation = 'top';
                obj.ax_tune.ThetaDir = 'clockwise';
                thetaticks(obj.ax_tune,[-90 -45 0 45 90])
                thetaticklabels(obj.ax_tune,{'-90', '-45', '0','45', '90'})
                thetalim(obj.ax_tune,[-90 90]);
            end

            ang_values   = obj.params.(behav).ang_values;

            mod = (obj.mod_dd.Value);
            
            num_trial_tot = zeros(8,1);

            for counter_rec = 1:n_recs
                
                obj.bndl_txa{counter_rec}.Value = ...
                    obj.bndl{counter_rec}(obj.dd{counter_rec}.Value);

                obj.type_txa{counter_rec}.Value = ...
                    obj.type{counter_rec}{obj.dd{counter_rec}.Value};
                
                if not(obj.chbx{counter_rec}.Value)
                    continue;
                end
                
                unit_selected = obj.recs{counter_rec}.units_summary...
                    (obj.dd{counter_rec}.Value);

                % waveform
                wave = unit_selected.(mod).waveform.wave * obj.zm;
                sample_rate = 3e4;
            
                if isnan(wave)
                    continue;
                end
            
                x_ = unit_selected.ch_map.x * 4;
                y_ = unit_selected.ch_map.y * 100;
            
                n_ch = size(wave,1);
                n_sig = length(wave);
                wave_ch = unit_selected.(mod).waveform.ch;
            
                x = x_(wave_ch+1);
                y = y_(wave_ch+1);
                ch_num = unit_selected.ch_map.map(wave_ch+1);
            
                span_ind = (0:n_sig-1)/sample_rate;
                span_group_ = repmat([span_ind,nan],n_ch,1);
                span_group = reshape((span_group_*1e3+x)',1,n_ch*(n_sig+1));
            
                if not(isnan(wave))
                    wave_ = [wave, nan(n_ch,1)];
                    wave_group = reshape((wave_+y)',1,n_ch*(n_sig+1));
                end

                ch_map = arrayfun(@num2str, ch_num+1,'UniformOutput', 0);
                text(obj.ax_wav,x-100*n_sig/sample_rate,y,ch_map);

                plot(obj.ax_wav,span_group,...
                    wave_group,'LineWidth',2,'Color',obj.colors(counter_rec,:));
                hold(obj.ax_wav,'on');
                
                % isi dist
                histogram(obj.ax_isi,'BinCounts',...
                    unit_selected.(mod).isi,'BinEdges', ...
                    unit_selected.(mod).isi_edge,...
                    'DisplayStyle', 'stairs', 'EdgeColor',...
                    obj.colors(counter_rec,:),'FaceColor',...
                    'none', 'linewidth', 2);
                hold(obj.ax_isi,'on');
                
                % cross prob
                plot(obj.ax_cor,unit_selected.(mod).corr_span,...
                    unit_selected.(mod).corr,...
                    'LineWidth',2,'Color',obj.colors(counter_rec,:));
                hold(obj.ax_cor,'on');
                ylim(obj.ax_cor,[0,inf]);
            
                % plot raster
                raster = unit_selected.(mod).([behav '_raster']);
                inds_span = unit_selected.(mod).([behav '_raster_span']);
                for counter_dir = 1:numel(obj.mod_axes)
                    raster_ = raster{counter_dir}';
                    [x_axis_raster_, y_axis_raster_] = ESN_raster_plot_axes( ...
                        raster_, inds_span, 0.5);
                    if strcmp(mod,'ss')
                        yyaxis(obj.mod_axes{counter_dir},'left');
                        plot(obj.mod_axes{counter_dir}, x_axis_raster_(:), ...
                            y_axis_raster_(:)+num_trial_tot(counter_dir),...
                            'Color',obj.colors(counter_rec,:),...
                            'LineWidth',1,'LineStyle','-','Marker','none');
                        firing_SS_ = mean(raster_) * 1000;
                        firing_SS_ = ESN_smooth(firing_SS_);
                        firing_SS_(firing_SS_<0)=0;
                        % plot SS rate
                        yyaxis(obj.mod_axes{counter_dir},'right');
                        plot(obj.mod_axes{counter_dir}, inds_span,...
                            firing_SS_,'Color',obj.colors(counter_rec,:),...
                            'LineWidth',2,'LineStyle','-','Marker','none');
                        ylim(obj.mod_axes{counter_dir},[0,200])
                        obj.mod_axes{counter_dir}.YAxis(1).Color = 'k';
                        obj.mod_axes{counter_dir}.YAxis(2).Color = 'k';
                        
                    else
                        plot(obj.mod_axes{counter_dir}, x_axis_raster_(:), ...
                            y_axis_raster_(:)+num_trial_tot(counter_dir),...
                            'Marker','.','MarkerSize',10,...
                            'LineWidth',5,'Color',obj.colors(counter_rec,:));
                        obj.mod_axes{counter_dir}.YAxis(1).Color = 'k';
                        obj.mod_axes{counter_dir}.YAxis(2).Color = 'k';
                    end
                    hold(obj.mod_axes{counter_dir},'on');
                    num_trial_tot(counter_dir) = num_trial_tot(counter_dir) + size(raster_,1);
                    if counter_rec == 1
                        xline(obj.mod_axes{counter_dir},0,'LineStyle',...
                        '--','LineWidth',1);
                    end
                end

                CS_on_data = unit_selected.cs.([behav '_cs_on_data']);
                fr_amplitude = CS_on_data.CS_fr_avg;

                vonMises_std = CS_on_data.vonMises_std;
                CS_ang_avg = CS_on_data.CS_ang_avg;
                CS_rho_avg = CS_on_data.CS_rho_avg;
                std_curv_ang = (CS_on_data.CS_ang_avg-vonMises_std) : 2 : (CS_on_data.CS_ang_avg+vonMises_std);
                std_curv_amp = repmat(CS_rho_avg, length(std_curv_ang), 1);
            
                plot_data_amp_mean = [fr_amplitude,...
                    fr_amplitude(1), nan]';
                plot_data_deg_mean = [ang_values, ang_values(1), nan]';

                polarplot(obj.ax_tune, deg2rad(plot_data_deg_mean),...
                    plot_data_amp_mean,...
                    'LineWidth', 1, 'Color', obj.colors(counter_rec,:),...
                    'LineStyle','-','Marker','none');
                polarplot(obj.ax_tune, deg2rad(std_curv_ang),...
                    std_curv_amp, '-', ...
                    'LineWidth', 1.5, 'Color', obj.colors(counter_rec,:),...
                    'LineStyle','-','Marker','none');
                polarplot(obj.ax_tune, [0 deg2rad(CS_ang_avg)],...
                    [0 CS_rho_avg], '-', ...
                    'LineWidth', 1.5, 'Color', obj.colors(counter_rec,:),...
                    'LineStyle','-','Marker','none');
                hold(obj.ax_tune,'on');
            
            end

            for counter_dir = 1:numel(obj.mod_axes)
                yyaxis(obj.mod_axes{counter_dir},'left');
                if num_trial_tot(counter_dir) > 0
                    ylim(obj.mod_axes{counter_dir},[.5,num_trial_tot(counter_dir)+.5]);
                end
            end

            if exist('ch_num','var')
                n_recs = numel(obj.rb);
                
                for counter_rec = 1:n_recs
                    if obj.rb{counter_rec}.Value
                        rec_selected = n_recs - counter_rec + 1;
                    end
                end
                
                unit_selected = obj.recs{rec_selected}.units_summary...
                    (obj.dd{rec_selected}.Value);
    
                ch_selected = str2double(unit_selected.name(15:16));
    
                elec_mat = [obj.elec.x' obj.elec.y'];
                idx_knn = knnsearch(elec_mat,...
                    [obj.elec.x(ch_selected) obj.elec.y(ch_selected)],...
                    'K',obj.KNN,'Distance','euclidean');
    
                idx_knn_ch = ismember(ch_num+1,idx_knn);
                ch_map_KNN = arrayfun(@num2str, ...
                    ch_num(idx_knn_ch)+1,...
                    'UniformOutput', 0);
                text(obj.ax_wav,x(idx_knn_ch)-100*n_sig/sample_rate,...
                    y(idx_knn_ch),ch_map_KNN,'Color','red');
    
                obj.bndl_lbl = max(cellfun(@max,obj.bndl));
                obj.lbl_max_lbl.Text = num2str(obj.bndl_lbl);
            end
            
            end
    end
end