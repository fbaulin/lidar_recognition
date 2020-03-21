function reduction_assessment_cv(varargin)
    fprintf('(%s) %s\n', datetime('now','Format','HH:mm:ss'), ' ======= ������ crossval_check =======')
    kwargs = KeywordArguments(...
        'dset_folder_name', 'test', ...           
        't_select_v', 1, ...                     % ����� ������ ��������������
        't_types', 'legacy', ...
        'reduction_method', 'prat',...          % ����� �������� ������������ ���������
        'n_metrics_fraction', 0.03,...              % ����� ������, ��� ������� ��������� ����������
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
        n_nearest, max_dim, step, imp_hwidth_ns, mse_noise, t_bfuns, ...
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
            'reduction_method', reduction_method,'n_nearest', n_nearest );  
        qcheck.setup_classifier( 'classifier_type',clf_type, 'performance',clf_perf_score,...
            'layers',clf_layers, 'gpu', use_gpu );
        [ fs_perf, fs_chars, fs_maps, conf_mx ] = ...
            cv_check(fvs, meta, qcheck, 'cv_method', cv_method, 'cv_param', cv_param);

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

% ���������� �������� �� ������������ �������
function [ fs_perf, fs_chars, varargout ] = cv_check(fvs, meta, qcheck, varargin)
    kwargs = KeywordArguments('cv_method','kfold','cv_param',4);
    [cv_method, cv_param] = kwargs.parse_input_cell(varargin);
    if strcmp(cv_method,'holdout'), cv_param=1/cv_param; end    % ���� ���������� �������, �� ����������� � ����
    [x,y] = qcheck.unwrap_cell_data(fvs, meta, 'matrix_output', false); % ������������ �������
    cv_results = crossval(...
        @(x_tn,y_tn,x_ts,y_ts)  eval_wrapper(qcheck,x_tn, y_tn, x_ts, y_ts),...
        x,y,cv_method,cv_param);%'kfold',4);  % ����������� �������� ��������
    switch nargout
        case 3  
            [ fs_perf, fs_chars, varargout{1} ] = ...
                process_cv_results(cv_results); 
        case 4  
            [ fs_perf, fs_chars, varargout{1}, varargout{2} ] = ...
                process_cv_results(cv_results); 
    end
    if nargout==2

    end
end

% ���������� ����������� �� �������� 
function [perf_score, rho_min, varargout] = process_cv_results(cv_results)
    perf_score_mx = horzcat(cv_results{:,1});
    rho_max_mx = horzcat(cv_results{:,2});
    perf_score = mean(perf_score_mx, 2);            % �������� ������������� �����������

    % ���������� ���������� � ��������� ��������� �������� ������
    [rhos_max, best_fold] = max(rho_max_mx,[],2);   % ��������� ������ ����������
    if nargout>=3   % ���� ��������� ������� ��������� ����������� ���������
        fs_maps = zeros(size(cv_results{1,3}));
        for i_dim = 1:length(best_fold)
            fs_maps(i_dim,:) = cv_results{best_fold(i_dim),3}(i_dim,:);
        end
        varargout{1} = fs_maps;
    end
    if nargout>=4   % ���� ��������� ������� ������
        [~, ind_dim] = max(rhos_max); % ������ ����������� ����������� (�.�. ������ ������ ������)
        varargout{2} = cv_results{ best_fold(ind_dim),4};
    end

    rho_min = rhos_max;

end

% ������� ��� ���������� �������� � ���� ���������� �������� �� ������
function output = eval_wrapper(obj,x_tn,y_tn,x_ts,y_ts)
%eval_wrapper ������� ��� ��������� ���������� �����������
    [score, rho_min, fs_maps, conf_mx ] = obj.assess_quality(x_tn,y_tn,x_ts,y_ts);
    output = { score, rho_min, fs_maps, conf_mx };   % ���� ������� 
end



