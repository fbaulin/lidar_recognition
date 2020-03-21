%% �������� �������� �� �����-��������� �������
fprintf('(%s) %s\n', datetime('now','Format','HH:mm:ss'), ' ======= ������ crossval_check =======')
%% ������������� 


if ~exist('external_input','var')
    
    fprintf('(%s) %s\n', datetime('now','Format','HH:mm:ss'), '��������� �� ���������� �����')
    % ������ ���������
    clearvars;
    %%% �������� ������
    dset_folder_name = 'az_scan_30';
    t_select_v = [ 3 ];                     % ����� ������ ��������������
   

    % for rm = { 'minalien' 'mhist' }%'nmin' 'buryi' }

    %%% ��������� ��������� ������
    reduction_method = 'minalien';          % ����� �������� ������������ ���������
    n_metrics_fraction = 0.03;              % ����� ������, ��� ������� ��������� ����������
    hist_max_value = 1;
    n_hist_edges = 51;
    n_nearest = 10;
    k_alien_max = 4;
    
    max_dim = 8;
    step = 1;
    %%% ��������� �������
    imp_hwidth_ns = 5;          % ����������������� �������� �� ����������� ������, ��
    mse_noise = 0.0001;
    show_graphics=1;
    t_bfuns = {'x'      'x'     'x'     'dtf2'  'dtf2' 'sym5'   'x'};
    cv_method = 'kfold';
    cv_param = 4;
end
t_types = {'none'   'fft'   'afft'  'cwt'   'acwt' 'wt'     'pca'};
t_types = t_types(t_select_v);
t_bfuns = t_bfuns(t_select_v); 

dimensions_init = 0:step:max_dim;
hist_edges = single(linspace(0,hist_max_value,n_hist_edges));    % �������� ������ ��� �����. �����������
%% ������������ ��������� ������������
switch reduction_method
    case 'mhist',       red_param = [num2str(n_metrics_fraction) '-' num2str(hist_max_value)];
    case 'nmin',        red_param = num2str(n_metrics_fraction);
    case {'buryi','fisher'},            red_param = '';
    case {'minalien','prat','wprat'},   red_param = num2str(n_nearest);
    otherwise,          error('����������� ����� ��������')
end

indexing_string = ...
    [ cv_method ...             �������� �������
    '_t' num2str(imp_hwidth_ns,'%02.1f') ...    t ��������
    '_mse' num2str(max(mse_noise)) ...          ��� ����
    '_' reduction_method num2str(red_param)...       ������ � �������� ����
    '_dim' num2str(step) '-' num2str(max_dim) ...         ������������ �����������
    ]; 


%% ������
n_types = length(t_types);  % ����� ��������������
fs_perf_cell = cell(1, n_types);    fs_chars_cell = cell(1, n_types);
fs_maps_cell = cell(1, n_types);    feat_header_cell = cell(1, n_types);
conf_mx_cell = cell(1, n_types); 

for i_type = 1:n_types
  % �������� ������
    t_type = t_types{i_type};   % ������� ��� ������
    t_bfun = t_bfuns{i_type};   % ������ ������� ������� ��������������
    in_filename = [ dset_folder_name '/fvs/' t_type '-' t_bfun ...
            '_t' num2str(imp_hwidth_ns,'%02.1f') '_mse' num2str(max(mse_noise)) '.mat'];
    load(in_filename)           % ��������� ������
  % ��������������  
    fprintf('(%s) %s\n', datetime('now','Format','HH:mm:ss'),[ '�������� ' t_type ]) 
    qcheck = QualityCheck(...                       % ������ �������� ��������
        'dimensions',dimensions_init,...
        'reduction_method', reduction_method,...
        'k_alien', k_alien_max, 'n_nearest', n_nearest, ...
        'hist_edges', hist_edges, 'n_metrics_fraction', n_metrics_fraction      );  
    qcheck.setup_classifier( ...                    % ��������� ���������������
        'classifier_type','pnet', ...
        'performance','default',...
        'layers', [15 7], ...
        'gpu', 'no'     );
    [ fs_perf, fs_chars, fs_maps, conf_mx ] = ...
        cv_check(fvs, meta, qcheck, '');
    
    fs_perf_cell{i_type} = fs_perf;     fs_chars_cell{i_type} = fs_chars;    
    fs_maps_cell{i_type} = fs_maps;     conf_mx_cell{i_type} = conf_mx;
    feat_header_cell{i_type} = feat_header; 
    
    obj_names = cell(1, length(meta));
    for i_o = 1:length(meta), obj_names{i_o} = meta{1,i_o}(1).name; end
    
%     fprintf('(%s) %s\n', datetime('now','Format','HH:mm:ss'),...
%     ['�������� ' t_type ' ���������'])
  % ����������
    output_file = [ dset_folder_name '/results/' indexing_string '_' t_type '_' t_bfun '.mat'];
    save(output_file, ...
        'fs_perf', 'fs_chars', 'fs_maps', 'feat_header', 'conf_mx', ...
        'obj_names', 't_type', 't_bfun', 'mse_noise', 'imp_hwidth_ns', 'reduction_method')
    fprintf('(%s) sep char: [%s]\n', datetime('now','Format','HH:mm:ss'),num2str(fs_chars.', '%6.3f '));
    fprintf('(%s) clf perf: [%s]\n', datetime('now','Format','HH:mm:ss'),num2str(fs_perf.', '%6.3f '));
    fprintf('\t�������� ���� � ������������ %s\n',output_file)
        
end
clear n_filename fs_perf fs_chars fs_maps i_type t_type

if show_graphics    % ������ ������������
    extern_input = 1; 
    crossval_check_vis; 
    clear external_input
end 

%% ��������� ������� 

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
    if nargout>=3
        fs_maps = zeros(size(cv_results{1,3}));
        for i_dim = 1:length(best_fold)
            fs_maps(i_dim,:) = cv_results{best_fold(i_dim),3}(i_dim,:);
        end
        varargout{1} = fs_maps;
    end
    if nargout>=4 
        [~, ind_dim] = max(rhos_max); % ������ ������ ����������� (�.�. ������ ������ ������)
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
