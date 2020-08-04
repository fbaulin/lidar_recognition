classdef Crossvalidation
    %Crossvalidation ����� ��� ���������� ������������ ������ �� ������
    %���������
    %   � ������ ���������� ������ ��� ���������� �������������� �������
    %   ����������� ��������� � ������������ ����������� ����� �������
    
    properties(Constant=true, Hidden=true)
        marker_type = {  '+' 'x' 's' 'o' '^' 'v' 'd' 'p' 'h' '*' '.' '+' 'x' 's' 'o' '^' 'v' 'd' 'p' 'h' '.' '*' };
    end
    
    methods
        function obj = Crossvalidation()
            %Crossvalidation �����������
            %   ����� �����������
            
        end
    end
    
    methods(Static=true)
        
        function reduction(varargin)
        %reduction ���������� �������������� ������� ����������� ���������
        %   ����������� �������������� ������� ����������� ��������� �
        %   �������������� ������������ �������� � ������������� ��
        %   ��������������. ���� �������� ��������� ����������� �� ������� ������ �
        %   /dataset_folder_name/fvs, � ����������� � /dataset_folder_name/results
        %
        %   ���������: 
        %   ���, ��������_��-��������� - ��������
        %   'dset_folder_name', 'test' - ������� ��������� ������� ����� ����� �� ����� dset_folder_name/fvs 
        %   't_select_v', 1 - ������ � ������� �������� �������������� ��������������                    
        %   't_types', 'legacy" - ������ ����������� ��������������
        %   'reduction_method', 'prat' - ����� ��������:
        %                       prat - ������������ �������� �����������
        %                       buryi - ��������������� ����������� �������
        %                       nmin - ��������������� n-� �������
        %                       fisher - ����� ������� ������
        %   'n_metrics_fraction', 0.03 -  ���� ������, ��� ������� ��������� ����������
        %   'n_nearest', 10- ����� ��������� ������� �� ������� ����������� �����������
        %   'max_dim', 8   - ������������ ����������� ����������� 
        %   'step', 1      - ��� � ������� ����������� ������ �� ������������
        %   'imp_hwidth_ns', 1 - ������������ �� �� ����������� ������, ��
        %   'mse_noise', 0.019 - ��� ����������� ���� ��� �������� ������� ��
        %   't_bfuns', {{'x','x','x','dtf2','dtf2','sym5','x'}} - ������������ �������� �-���
        %   'cv_method', 'kfold' - ��� ��������:
        %                       kfold - ������������, 
        %                       holdout - ���������� ������� 
        %   'cv_param', 4 - ����� n ������� ��� �������� (��� holdout - 1/n ���� �����. �������)
        %   'classifier', 'pnet' - ��� ��������������: 
        %                       pnet - �� � �������������� �������� ����� (��. patternnet()) 
        %                       mlp - �� � �������� �������� ����� (��. feedforwardnet()) 
        %   'layers', [ 8 4 ] - �������������� ����
        %   'performance', 'default' - �������� ������ ��������:
        %                       default - �������� �� ���������
        %                       crossentropy - ������������ �������� (������ ��� pnet) 
        %                       mse - ��� ������ (������ ��� mlp)
        %   'use_gpu', 'no' - ������������ ������������ ��� ��������� ����������
            fprintf('(%s) %s\n', datetime('now','Format','HH:mm:ss'), ' ======= ������ crossval_check =======')
            kwargs = KeywordArguments(...
                'dset_folder_name', 'test', ...           
                't_select_v', 1, ...                     % ����� ������ ��������������
                't_types', 'legacy', ...
                'reduction_method', 'prat',...          % ����� �������� ������������ ���������
                'n_metrics_fraction', 0.03,...              % ����� ������, ��� ������� ��������� ����������
                'k_alien', 3,...
                'n_nearest', 10,...
                'max_dim', 8,...
                'step', 1,...
                'imp_hwidth_ns', 1,...          % ����������������� �������� �� ����������� ������, ��
                'mse_noise', 0.019,...
                't_bfuns', {{'x','x','x','dtf2','dtf2','sym5','x'}},...
                'cv_method', 'kfold',...
                'cv_param', 4 ,...
                'classifier', 'pnet',...
                'layers', [ 8 4 ],...
                'performance', 'default',...
                'use_gpu', 'no');
            [   dset_folder_name, t_select_v, t_types, reduction_method, n_metrics_fraction,...
                k_alien, n_nearest, max_dim, step, imp_hwidth_ns, mse_noise, t_bfuns, ...
                cv_method, cv_param, clf_type, clf_layers, clf_perf_score, use_gpu] ...
                = kwargs.parse_input_cell(varargin);

            if strcmp(t_types, 'legacy')
                t_types = {'none'   'fft'   'afft'  'cwt'   'acwt'   'wt'    'pca'};
                t_types = t_types(t_select_v);  t_bfuns = t_bfuns(t_select_v);
                warning('���� ����� ����� �������������� ���������� � ����� ����� ������')
            end

            % ������������ ��������� ������������
            switch reduction_method
                case 'nmin',        red_param = num2str(n_metrics_fraction);
                case {'buryi','fisher'},            red_param = '';
                case {'minalien','prat','wprat'},   red_param = num2str(n_nearest);
                case 'refmmin',                   red_param = [num2str(k_alien) 'of' num2str(n_nearest)];
                otherwise,          error('����������� ����� ��������')
            end
            indexing_string = ...
                [ cv_method ...             �������� �������
                '_t' num2str(imp_hwidth_ns,'%02.1f') ...        t ��������
                '_mse' num2str(max(mse_noise)) ...              ��� ����
                '_' reduction_method num2str(red_param)...      ������ � �������� ����
                '_dim' num2str(step) '-' num2str(max_dim) ...   ������������ �����������
                ]; 
            n_types = length(t_types);  % ����� ��������������
            dimensions_init = 0:step:max_dim;
            for i_type = 1:n_types
                % �������� ������
                t_type = t_types{i_type};   % ������� ��� ������
                t_bfun = t_bfuns{i_type};   % ������ ������� ������� ��������������
                in_filename = [ dset_folder_name '/fvs/' t_type '-' t_bfun ...
                        '_t' num2str(imp_hwidth_ns,'%02.1f') '_mse' num2str(max(mse_noise)) '.mat'];
                load(in_filename,'fvs','meta','feat_header');   % ��������� ������
                % ��������������  
                fprintf('(%s) %s\n', datetime('now','Format','HH:mm:ss'),[ '�������� ' t_type ]) 
                qcheck = QualityCheck('dimensions',dimensions_init, ... % ������ �������� ��������
                    'reduction_method', reduction_method,...
                    'n_metrics_fraction',n_metrics_fraction,...
                    'k_alien', k_alien, 'n_nearest', n_nearest );  
                qcheck.setup_classifier( 'classifier_type',clf_type, 'performance',clf_perf_score,...
                    'layers',clf_layers, 'gpu', use_gpu );
                [ fs_perf, fs_chars, fs_maps, conf_mx ] = ...
                    Crossvalidation.cv_check(fvs, meta, qcheck, 'cv_method', cv_method, 'cv_param', cv_param);

                obj_names = cell(1, length(meta));
                for i_o = 1:length(meta), obj_names{i_o} = meta{1,i_o}(1).name; end
                % ����������
                output_file = [ dset_folder_name '/results/' indexing_string '_' t_type '_' t_bfun '.mat'];
                save(output_file, ...
                    'fs_perf', 'fs_chars', 'fs_maps', 'feat_header', 'conf_mx', ...
                    'obj_names', 't_type', 't_bfun', 'mse_noise', 'imp_hwidth_ns', 'reduction_method')
                fprintf('(%s) sep char: [%s]\n', datetime('now','Format','HH:mm:ss'),num2str(fs_chars.', '%6.3f '));
                fprintf('(%s) clf perf: [%s]\n', datetime('now','Format','HH:mm:ss'),num2str(fs_perf.', '%6.3f '));
                fprintf('\t�������� ���� � ������������ %s\n',output_file)
            end
        end
                
        function crossval_list_vis(filelist,varargin)
        %crossval_list_vis - ������������ ����������� �������
        %   ������� �������� �������� ������, �������������� ������� reduction, �
        %   ���������� �������� �������
        % 
        %   ���������:
        %   1) filelist - ������ ������� ������ � ������������ �������
        %               {'.\file1', '.\file2'} - ������ � char ��������� - ��������
        %               [] - ������� ��������� ��� ������ ������
        %   �����������:
        %   'maps', false       - ���������� ����� ������ ���������
        %   'confusion', false  - ���������� ������� ������
        %   'clf_score', 'erroe_prob'   - �������������� �������� �������������:
        %               error_prob  - ���� ��������� �������
        %               clf         - ���������� �������������� ��� �������
        %               �������������� (crossentropy/mse)
            kwargs = KeywordArguments(...
                    'maps',false,...
                    'confusion',false,...
                    'clf_score', 'error_prob');
            [   maps_flag, confusion_mx_flag, clf_score ] = kwargs.parse_input_cell(varargin);
            if ~iscell(filelist), filelist = Crossvalidation.ui_open_files(); end
            n_files = length(filelist);
            compare_params =    {'t_type' 't_bfun' 'imp_hwidth_ns' 'mse_noise' 'reduction_method'};
            compare_params_tex ={''       ''       '\tau_n_s='         '\sigma_n='   'Q - '};
            datastruct = Crossvalidation.extract_data(filelist);

            figure('Name','comp red','color', 'white','WindowStyle','docked'); 
            if maps_flag
                n_columns = n_files;
                ax_fs_score = subplot(2,n_columns,1:ceil(n_columns/2));
                ax_clf_score = subplot(2,n_columns,ceil(n_columns/2)+1:n_columns);
                ax_fs_maps = arrayfun(@(i) subplot(2,n_columns,n_columns+i), 1:n_files);
            else
                ax_fs_score = subplot(1,2,1);
                ax_clf_score = subplot(1,2,2);
            end

            [leg_str, tit_str] = Crossvalidation.compile_legend(datastruct, compare_params, compare_params_tex);
            Crossvalidation.draw_fs_score(datastruct, leg_str, tit_str, ax_fs_score);
            switch clf_score 
                case 'clf', Crossvalidation.draw_clf_score(datastruct, leg_str, tit_str, ax_clf_score);
                case 'error_prob' , Crossvalidation.draw_p_err(datastruct, leg_str, tit_str, ax_clf_score);
                otherwise, error('�� �������������� ��� �������� �������������');
            end
            if maps_flag,   Crossvalidation.draw_fs_maps(datastruct, ax_fs_maps, leg_str); end
            if confusion_mx_flag,   Crossvalidation.confusion_mx(datastruct,leg_str); end
        end
        
        function compare_nmins(filelist,varargin)
        %compare_nmins - ������� ��� ������ ������������ �������� ���� ������
        %   ������� ������������ ��� ���������� �������� ��� ������ �������� ����
        %   ������ � ����������� ��������� ���� ������
            kwargs = KeywordArguments(...
                'dimensions', 2:6, ...
                'maps',0);
            [  dims, maps_flag ] = kwargs.parse_input_cell(varargin);
            if ~iscell(filelist), filelist = Crossvalidation.ui_open_files(); end
            n_files = length(filelist);
            compare_params =    {'t_type' 't_bfun' 'imp_hwidth_ns' 'mse_noise' 'reduction_method'};
            compare_params_tex ={''       ''       '\tau_n_s='         '\sigma_n='   'Q - '};
            datastruct = Crossvalidation.extract_data(filelist);

            figure('Name','comp red','color', 'white','WindowStyle','docked'); 
            if maps_flag
                n_columns = n_files;
                ax_clf_score = subplot(2,n_columns,1:ceil(n_columns/2));
                ax_clf_score_nminwise = subplot(2,n_columns,ceil(n_columns/2)+1:n_columns);
                ax_fs_maps = arrayfun(@(i) subplot(2,n_columns,n_columns+i), 1:n_files);
            else
                ax_clf_score = subplot(1,2,1);
                ax_clf_score_nminwise = subplot(1,2,2);
            end

            [leg_str, tit_str] = Crossvalidation.compile_legend(datastruct, compare_params, compare_params_tex);
            Crossvalidation.draw_clf_score(datastruct, leg_str, tit_str, ax_clf_score);
            Crossvalidation.draw_clf_score_filewise(datastruct, ax_clf_score_nminwise, dims);
            if maps_flag,   Crossvalidation.draw_fs_maps(datastruct, ax_fs_maps, leg_str); end
        end
        
        %TODO �������
        function [] = metric_histogram(varargin)
        %HIST_COUNTS Summary of this function goes here
        %   Detailed explanation goes here
            kwargs = KeywordArguments(...
                'edges', 'auto', ...
                'file', 'none', ...
                'show', 1 );
            [hist_edges, file, show_flag] =  ...
                kwargs.parse_input_cell(varargin);
            [fvs_mx, targets] = Crossvalidation.open_fv_base(file);
            
            rhos = Crossvalidation.calc_metrics(fvs_mx, targets);
            if strcmp(hist_edges, 'auto')
                hist_edges = linspace(0,max(rhos),100);
            elseif isnumeric(hist_edges)
                if length(hist_edges)==1
                    hist_edges = linspace(0,max(rhos),hist_edges);
                end
            else, error('������ ������� �����������')
            end
            [hist_counts, ~] = histcounts(rhos, hist_edges);
            if show_flag,   Crossvalidation.integral_hist_vis(hist_counts, hist_edges, file);    end
        end
        
        function [] = compare_fs_hists(varargin)
        %compare_fs_hists - ��������� ��������� ��������, ������������ ��� ��������
        %���������
        %
        %   ����������� ���������:
        %   'file', 'none',         - ������ ������
        %           {'.\file1', '.\file2'} - ������ � char ��������� - ��������
        %           none - ������� ��������� ��� ������ ������
        %   'edges', 'auto'         - ������ ��������� �����������
        %               auto    - 99 ���������� �� 0 �� ����. ����. �������
        %               n       - n ���������� �� 0 �� ����. ����. �������
        %               [array] - ������ ������� ����� �������
        %   'feature_sets', 'all'   - ������� ��������
        %               all             - ��� ��������
        %               {feat1, ...}    - ����� ���������
        %               [fnum1, ...]    - ���������� ������ ���������
        %   'axis', 'new'   - ������ ��� ���������
        %               'new'  - ����� ������
        %               ax_obj - ���������� ������ �� �������� ���
            kwargs = KeywordArguments(...
                'file', 'none', ...
                'edges', 'auto', ...
                'feature_sets', 'all',...
                'axis', 'new');
            [file, hist_edges, feature_query, ax] =  ...
                kwargs.parse_input_cell(varargin);

            [fvs_mx, targets, file_feature_names] = Crossvalidation.open_fv_base(file);
            query_maps = Crossvalidation.process_feature_query(feature_query, file_feature_names);
            n_fs = size(query_maps,1);
            rhos = cell(n_fs,1);
            for i_fs = 1:n_fs
                rhos{i_fs} = Crossvalidation.calc_metrics(fvs_mx(:,query_maps(i_fs,:)), targets).';
            end
            rhos = cell2mat(rhos);
            if strcmp(hist_edges, 'auto')
                hist_edges = linspace(0,max(rhos,[],'all'),101);
            elseif isnumeric(hist_edges)
                if length(hist_edges)==1
                    hist_edges = linspace(0,max(rhos,[],'all'),hist_edges+1);
                else	% hist_edges = hist_edges ������������ �������� ������
                end
            else, error('������ ������� �����������')
            end
            hist_counts = zeros(n_fs, length(hist_edges)-1);
            for i_fs = 1:n_fs
                [hist_counts(i_fs,:), ~] = histcounts(rhos(i_fs,:), hist_edges);
            end
            Crossvalidation.metric_prob_int(ax, hist_counts, hist_edges, feature_query);

        end
    end
    
        
    methods(Static=true, Hidden=true)
        
        function output_args = ui_open_files()
        %FGETFILELIST �������� ������ ������
        %   ������� ���������� ����������� ��������� ��� ������ ������ �
        %   ��������� �� ������ ������ cell � ������� ������� ������ ���
        %   cell � �����
            [filenames,path] = ...
                uigetfile(...
                ['*.', 'mat'],'MultiSelect','on');  % ������� ����� ����� ���
            if ~iscell(filenames)                       % ������� <2 ������
                if (size(filenames,2)==1);filenames={0};% 0 ������
                else; filenames = {[path,filenames]};   % 1 ���� - ��������� ������ ���
                end
            else                                        % ��������� ������
                for i=1:(size(filenames,2))
                    filenames{i} = [path,filenames{i}];
                end
            end
            output_args = filenames; % ������� ������
        end
        
        % ���������� �������� �� ������������ �������
        function [ fs_perf, fs_chars, varargout ] = cv_check(fvs, meta, qcheck, varargin)
            kwargs = KeywordArguments('cv_method','kfold','cv_param',4);
            [cv_method, cv_param] = kwargs.parse_input_cell(varargin);
            if strcmp(cv_method,'holdout'), cv_param=1/cv_param; end    % ���� ���������� �������, �� ����������� � ����
            [x,y] = qcheck.unwrap_cell_data(fvs, meta, 'matrix_output', false); % ������������ �������
            cv_results = crossval(...
                @(x_tn,y_tn,x_ts,y_ts)  Crossvalidation.eval_wrapper(qcheck,x_tn, y_tn, x_ts, y_ts),...
                x,y,cv_method,cv_param);%'kfold',4);  % ����������� �������� ��������
            switch nargout
                case 3  
                    [ fs_perf, fs_chars, varargout{1} ] = ...
                        Crossvalidation.process_cv_results(cv_results); 
                case 4  
                    [ fs_perf, fs_chars, varargout{1}, varargout{2} ] = ...
                        Crossvalidation.process_cv_results(cv_results); 
            end
        end

        % ���������� ����������� �� �������� 
        function [clf_perf_score, rho_min, varargout] = process_cv_results(cv_results)
            clf_perf_score_mx = horzcat(cv_results{:,1});
            fs_score_mx = horzcat(cv_results{:,2});
            clf_perf_score = mean(clf_perf_score_mx, 2);            % �������� ������������� �����������

            % ���������� ���������� � ��������� ��������� �������� ������
            [fs_score_max, best_fold] = max(fs_score_mx,[],2);   % ��������� ������ ����������
            if nargout>=3   % ���� ��������� ������� ��������� ����������� ���������
                fs_maps = zeros(size(cv_results{1,3}));
                for i_dim = 1:length(best_fold)
                    fs_maps(i_dim,:) = cv_results{best_fold(i_dim),3}(i_dim,:);
                end
                varargout{1} = fs_maps;
            end
            if nargout>=4   % ���� ��������� ������� ������
                varargout{2} = sum(cat(4,cv_results{:,4}),4); % ����� ������� ������
                
%                 [~, ind_dim] = max(fs_score_max); % ������ ����������� ����������� (�.�. ������ ������ ������)
%                 varargout{2} = cv_results{ best_fold(ind_dim),4};
                
            end

            rho_min = fs_score_max;

        end

        % ������� ��� ���������� �������� � ���� ���������� �������� �� ������
        function output = eval_wrapper(obj,x_tn,y_tn,x_ts,y_ts)
        %eval_wrapper ������� ��� ��������� ���������� �����������
            [score, rho_min, fs_maps, conf_mx ] = obj.assess_quality(x_tn,y_tn,x_ts,y_ts);
            output = { score, rho_min, fs_maps, conf_mx };   % ���� ������� 
        end
        
        %������� ���� � ��������� ���������
        function [varargout] = open_fv_base(file)
            
            if ~strcmp(file,'none')
                switch file(end-3:end)
                    case '.mat'
                        in_data = load(file, 'fvs', 'meta','feat_header');
                        [fvs_mx, targets] = QualityCheck.unwrap_cell_data(...
                            in_data.fvs, in_data.meta);
                        file_features = in_data.feat_header;
                    case '.csv'
                        in_data = readtable(file);
                        fvs_mx = table2array(in_data(:,1:end-3));
                        targets = table2array(in_data(:,end));
                        file_features = in_data(:,1:end-3).Properties.VariableNames;
                    otherwise
                        error('�� �������������� ��� �����');
                end
            else
                filelist = Crossvalidation.ui_open_files();
                file = filelist{1};
                in_data = load(file, 'fvs', 'meta','feat_header');
                        [fvs_mx, targets] = QualityCheck.unwrap_cell_data(...
                            in_data.fvs, in_data.meta);
                file_features = in_data.feat_header;
            end
            if ischar(file_features)
                file_features = strsplit(file_features,','); 
                file_features(end) = [];
            end
            varargout{1} = fvs_mx;
            if nargout>=2, varargout{2} = targets; end
            if nargout>=3, varargout{3} = file_features; end
        end
                
        %% ��������� ������� ������������
        % Visualize fs score
        function draw_fs_score(datastruct, leg_str, tit_str, varargin)
            if nargin>=4, axes(varargin{1}); end
            n_files = size(datastruct,2); 
            hold on;
            curve_h = zeros(1,n_files);
            for i_file = 1:n_files
                curve_h(i_file) = stem(datastruct(i_file).('dimensions'),datastruct(i_file).('fs_chars'),...
                    [ ':' Crossvalidation.marker_type{i_file} 'k']);
                plot(datastruct(i_file).('dimensions'),datastruct(i_file).('fs_chars'), ':');
            end
            xticks(datastruct(i_file).('dimensions'))
            ylabel(datastruct(i_file).('reduction_method'))
            legend(curve_h, leg_str) 
            title(['�������� ������ Q (' tit_str ')'])
        end

        % Visualize clf score
        function draw_clf_score(datastruct, leg_str, tit_str, varargin)
            if nargin>=4, axes(varargin{1}); end
            hold on;
            n_files = size(datastruct,2); 
            curve_h = zeros(1,n_files);

            for i_file = 1:n_files
                curve_h(i_file) = stem(...
                    datastruct(i_file).('dimensions'),datastruct(i_file).('fs_perf'),...
                    [ ':' Crossvalidation.marker_type{i_file} 'k']);
                plot(datastruct(i_file).('dimensions'),datastruct(i_file).('fs_perf'),':');
            end
            xticks(datastruct(i_file).('dimensions'))
            legend(curve_h,leg_str)
            title(['������ (' tit_str ')'])
        end

        % Visualize clf score
        function draw_clf_score_filewise(datastruct,  varargin)
            if nargin>=2, axes(varargin{1}); end
            if nargin>=3, dims = varargin{2}; else, dims = 2:6; end
            if nargin>=4, tit_str = varargin{3}; else, tit_str = 'nmin'; end
            hold on;
            n_files = size(datastruct,2); 

            
            files_sel = ones(1,n_files);
            param_values = zeros(1,n_files);
            curve_values = zeros(length(dims),n_files);
            
            for i=1:n_files 
                if isfield(datastruct(i),'reduction_param')
                    param_values(i) = datastruct(i).('reduction_param');
                    curve_values(:,i) = datastruct(i).('fs_perf')(dims+1).';
                else
                    warning(['���� ����������� ������ �������(' datastruct(i).('reduction_method') ')'])
                    files_sel(i)=0;
                end
            end
            param_values(files_sel==0)=[];
            curve_values(repmat(files_sel==0,length(dims),1))=[];
            [param_values, sort_ids] = sort(param_values);
            
            curve_values = curve_values(:,sort_ids);
            mesh(param_values, dims, curve_values)
                        
%             curve_h = zeros(1,length(dims));
%             for i_dim = 1:length(dims)
%                 curve_h(i_dim) = stem(...
%                     param_values,curve_values(i_dim,:).',...
%                     [ ':' Crossvalidation.marker_type{i_dim} 'k']);
%                 plot(param_values,curve_values(i_dim,:).',':');
%             end
            
%             if min(param_values)>0, set(gca, 'XScale', 'log'); end
%             xlabel('\rho_min'); ylabel('performance')
%             legend(curve_h,arrayfun(@ (d) num2str(d,'d=%u'),dims,'UniformOutput',false));

            xticks(param_values)
            
            xlabel('\rho_min'); ylabel('dims'); zlabel('performance')
            zlim([0 max(curve_values,[],'all')])
            axis vis3d; view(45,45);
            title(['���������� �. (' tit_str ')'])
        end 
        
        % Visualize fs maps
        function draw_fs_maps(datastruct, ax, leg_str)
            n_files = size(datastruct,2); 
            for i_file = 1:n_files
                feature_names = strsplit(strrep(datastruct(i_file).('feat_header'),'_',' '),',');
                axes(ax(i_file));  %#ok<LAXES>
                hold on;
                spy(datastruct(i_file).('fs_maps').',[Crossvalidation.marker_type{i_file} 'k']);
                xticks(ax(i_file),1:length(datastruct(i_file).('dimensions')))
                xticklabels(ax(i_file),datastruct(i_file).('dimensions'))
                xlabel(ax(i_file),'dimensionality, n')
                ylabel(ax(i_file),[datastruct(i_file).('t_type') ' features'])
                title(ax(i_file),['�������� ��� ' leg_str(i_file)])
                ylim(ax(i_file),[0 size(datastruct(i_file).('fs_maps'),2)+1])
                yticks(ax(i_file),1:size(datastruct(i_file).('fs_maps'),2))
                yticklabels(ax(i_file),feature_names)
            end
        end
        
        function draw_p_err(datastruct,leg_str, tit_str, varargin)
            if nargin>=4, axes(varargin{1}); end
            hold on;
            n_files = size(datastruct,2); 
            curve_h = zeros(1,n_files);

            for i_file = 1:n_files
                conf_mx = datastruct(i_file).('conf_mx');
                
                [n_dim, n_obj, ~] = size(conf_mx);
                mask = repmat(permute(eye(n_obj,'logical'),[3,1,2]),n_dim,1,1);
                total_samp = sum(conf_mx,[2 3]);
                conf_mx(mask)=0;
                err_samp = sum(conf_mx,[2 3]);
                p_err = err_samp./total_samp;
                
                curve_h(i_file) = stem(...
                    datastruct(i_file).('dimensions'),p_err,...
                    [ ':' Crossvalidation.marker_type{i_file} 'k']);
                plot(datastruct(i_file).('dimensions'),p_err,':');
            end
            xticks(datastruct(i_file).('dimensions'))
            legend(curve_h,leg_str)
            title(['P_{err} (' tit_str ')'])
        end
        
        % Visualize confusion matricies
        function confusion_mx(datastruct,leg_str)
            n_files = size(datastruct,2);
            n_cols = 2;
            n_rows = ceil(n_files/2);
            figure('Name','','color', 'white','WindowStyle','docked');
            for i_type = 1:n_files
                subplot(n_rows, n_cols, i_type)
                
                fs_perf = datastruct(i_type).('fs_perf');
                [~,i_best_dim] = min(fs_perf);
                cmx = datastruct(i_type).('conf_mx');
                cmx = permute(cmx(i_best_dim,:,:), [2 3 1]);
                cm = confusionchart(cmx,...
                    datastruct(i_type).('obj_names'));
                cm.RowSummary = 'row-normalized';
                cm.ColumnSummary = 'column-normalized';
                n_obj = size(cmx,1);
                perr = 1- sum(cmx(eye(n_obj,'logical')),'all')/sum(cmx,'all');
                title(horzcat(...
                    leg_str{i_type}, ', p_{err}(n=', ...
                    num2str(i_best_dim,'%d'), ')=', ...
                    num2str(perr,'%0.3f')))
            end
        end

        % Extraxt data
        function datastruct = extract_data(filelist)
            n_files = length(filelist);
            datastruct = cell(1,n_files);
            
            for i_file = 1:n_files
                datastruct{i_file} = load(filelist{i_file});
                datastruct{i_file}.('dimensions') = sum(datastruct{i_file}.('fs_maps'),2);
                file_name = filelist{i_file};
                crop_ind = strfind(file_name,'\');
                
                datastruct{i_file}.('filename') = file_name(crop_ind(end)+1:end-4);
                
                
                red_method = datastruct{i_file}.('reduction_method');
                switch red_method
                    case {'nmin','minalien','prat'}
                        crop_ind_start = strfind(filelist{i_file},['_' red_method])+1+length(red_method);
                        crop_ind_end = crop_ind_start+strfind(filelist{i_file}(crop_ind_start:end),'_')-2;
                        crop_ind_end = crop_ind_end(1);
                        file_name = filelist{i_file};
                        datastruct{i_file}.('reduction_param') = str2double(file_name(crop_ind_start:crop_ind_end));
                    otherwise
                        datastruct{i_file}.('reduction_param') = 0;
                end
                switch datastruct{i_file}.('reduction_method')
                    case {'nmin', 'buryi'}
                        datastruct{i_file}.('fs_chars')(2:end) = datastruct{i_file}.('fs_chars')(2:end)./sqrt(datastruct{i_file}.('dimensions')(2:end));
                    case {'minalien', 'prat'}
                        datastruct{i_file}.('fs_chars') = 1-datastruct{i_file}.('fs_chars');
                    case 'fisher'
                        datastruct{i_file}.('reduction_param') = 0;
                    otherwise, warning('�� ����� ����� ������������');
                end
                
                
                
            end
            
            datastruct = [datastruct{:}];
        end

        % COmpile titles and legends
        function [legend_cell, title_str] = compile_legend(dataset, params, varargin)
            if nargin==3,   tex_names =  varargin{1}; 
            else,           tex_names = params; 
            end
            n_expers = length(dataset);
            legend_cell = cell(1, n_expers);
            title_str = '';
            for i_param=1:length(params)
                param=params{i_param};
                param_vals = cellfun(@ (a) num2str(a), {dataset.(param)},'UniformOutput', false);
                c = unique(param_vals);
                if length(c)>1
                    for i_exper = 1:n_expers
                        legend_cell{i_exper} = [legend_cell{i_exper} tex_names{i_param} param_vals{i_exper} ','];
                    end
                else
                    title_str = [title_str,tex_names{i_param},c{:},',']; %#ok<AGROW>
                end

            end
            title_str(end)=[];
            
            if ~isempty(horzcat(legend_cell{:}))
                if length(unique(legend_cell))==n_expers
                    for i=1:length(legend_cell), legend_cell{i}(end)=[]; end
                else
                    legend_cell = cellfun(@(s) strrep(s,'_',' '),...
                        {dataset.('filename')},'UniformOutput',false);
                end
            else
                legend_cell = cellfun(@(s) strrep(s,'_',' '),...
                    {dataset.('filename')},'UniformOutput',false);
            end
            
        end
        
        %% ��������� ������� ��� ���������� ������

        function integral_hist_vis(hist_counts, hist_edges, filename)
            crop_ind_start = strfind(filename,'\');
            title_str = filename(crop_ind_start(end)+1:end-4);
            hist_bins = movmean(hist_edges,2,'Endpoints','discard');
            int_counts = cumsum([0 hist_counts]);
            int_counts = int_counts/int_counts(end);

            figure('Name','Hist analysis',...
                'color', 'white','WindowStyle','docked'); 
            hold on;
            marker_type = { '*' 's' 'd' '+' 'x' 'o'};
            line_type = {'-','--',':','-.', '-','--'};

            i_band =1; 
            plot(hist_edges,int_counts...
               ,[':' marker_type{i_band} 'k'],'LineWidth',1.2)

            Slozhno.setup_axis_labels('\rho_{th}', 'N(\rho>\rho_{th})')
            set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
            title(title_str);            
        end

        % ����� ��������� �����������
        function metric_prob_int(ax, bms_counts, hist_edges, feature_query )    
            if ischar(ax)
                figure('Name','Hist analysis',...
                    'color', 'white','WindowStyle','docked'); 
            else
                axes(ax);
            end
            hold on;
            
            line_type = {'-','--',':','-.'};
            
            n_fs = size(bms_counts,1);
            leg_str = cell(1,n_fs);
            curves_h = zeros(n_fs,1);
            for i_fs = 1:n_fs
                int_counts = cumsum([0 bms_counts(i_fs,:)]);
                int_counts = int_counts/int_counts(end);
                curves_h(i_fs) = plot(hist_edges(1:end),int_counts...
                   ,[line_type{mod(i_fs,4)+1} Crossvalidation.marker_type{i_fs} 'k'],'LineWidth',1.2);
                leg_str(i_fs) = join(feature_query{i_fs},',');
            end
            legend(curves_h, leg_str,'Interpreter','none')
            Slozhno.setup_axis_labels('\rho_{th}', 'R(\rho>\rho_{th})')
            
            
            set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
        end

        function rhos = calc_metrics(fvs_mx,trgs)
            [n_samples, n_dim] = size(fvs_mx);
            [~,~,trgs_id] = unique(trgs);
            alien_mask = and(...    % ������ �������
                trgs_id~=trgs_id.',...
                triu(ones(n_samples,'logical'),1));
            rhos = zeros(nnz(alien_mask),1,'like',fvs_mx);

            for i=1:n_dim
                rhos_cur = fvs_mx(:,i) - fvs_mx(:,i).';
                if any(diag(rhos_cur)>0)
                    warning('>0'); 
                end
                rhos = rhos+(rhos_cur(alien_mask)).^2;
            end
            rhos = sqrt(rhos);
        end
        
        function query_maps = process_feature_query(feature_query, feature_names)
            
            dim_fs = length(feature_names);
            if ischar(feature_query)
                if strcmp(feature_query,'all'), query_maps = ones(1,dim_fs,'logical');
                else, error('����������� ������� ��������� ������ ���������')
                end
            elseif iscell(feature_query)
                n_fs = length(feature_query);
                query_maps = zeros(n_fs,dim_fs,'logical');
                for i_fs = 1:n_fs
                    if iscell(feature_query{i_fs})
                        for j = 1:length(feature_query{i_fs})
                            query_maps(i_fs,strcmp(feature_names, feature_query{i_fs}{j}))=true;
                        end
                    elseif isnumeric(feature_query{i_fs})
                        query_maps(i_fs, feature_query{i_fs})=true;
                    elseif ischar(feature_query{i_fs})
                        query_maps(i_fs,strcmp(feature_names, feature_query{i_fs}))=true;
                    end
                end
            else, error('����������� ������� ��������� ������ ���������')
            end
        end
    end
end

