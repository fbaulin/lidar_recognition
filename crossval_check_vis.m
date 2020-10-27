%mhist_vis
%Сценарий визуализации результатов вычиселния mhist
%% Начальные установки
if ~exist('external_input','var') % ручной запуск скрипта
    clearvars
    dset_folder_name = 'az_scan_30';
    
    t_select_v = [ 3 ];                     % маска выбора объектов

    max_dim = 6;
    step=2;
    imp_hwidth_ns = 5;        % продолжительность импульса по половинному уровню, нс
    mse_noise = 0.0001;
    n_metrics_fraction = 0.03;
    cv_method = 'kfold'
    
end
t_types = {'none' 'fft' 'afft' 'cwt' 'acwt' 'wt' 'pca' 'dca'};
t_types = t_types(t_select_v);

t_names = {'no' 'Fourier' 'amplitude Fourier' 'comlex wavelet' 'complex wavelet' 'wavelet', 'PCA', 'DCA'};
t_names = t_names(t_select_v);

% for rm =  rm% 'buryi' 'nmin' }
%  reduction_method = 'mhist';

switch reduction_method
    case 'mhist',       red_param = [num2str(n_metrics_fraction) '-' num2str(hist_max_value)];
    case {'nmin','integrnmin'},        red_param = num2str(n_metrics_fraction);
    case {'buryi', 'fisher', 'auhist'},            red_param = '';
    case {'minalien','prat','wprat'},    red_param = num2str(n_nearest);
    otherwise,          error('Неизвестный метод редукции')
end

indexing_string = ...
    [ cv_method ...             исходная выборка
    '_t' num2str(imp_hwidth_ns,'%02.1f') ...    t импульса
    '_mse' num2str(max(mse_noise)) ...          СКО шума
    '_' reduction_method num2str(red_param)...       метрик в штрафной зоне
    '_dim' num2str(step) '-' num2str(max_dim) ...         максимальная размерность
    ]; 

%%
n_types = length(t_types);

%% Load data
fs_chars_cell = cell(1, n_types);
fs_maps_cell = cell(1, n_types);
feat_header_cell = cell(1, n_types);
dimensions_cell = cell(1, n_types);
fs_perf_cell = cell(1, n_types);
conf_mx_cell = cell(1, n_types);
obj_names_cell = cell(1, n_types);
for i_type = 1:n_types
    t_type = t_types{i_type};
    t_bfun = t_bfuns{i_type};
    in_filename = [ dset_folder_name '/results/' indexing_string '_' t_type '_' t_bfun '.mat'];
    in_data = load(in_filename);
    fs_perf_cell{i_type} = in_data.fs_perf;
    
    fs_maps_cell{i_type} = in_data.fs_maps;
    feat_header_cell{i_type} = in_data.feat_header;
    dimensions_cell{i_type} = sum(in_data.fs_maps,2);
    fs_chars_cell{i_type} = in_data.fs_chars;
    switch reduction_method
    case {'nmin', 'buryi','auhist'}
        fs_chars_cell{i_type}(2:end) = fs_chars_cell{i_type}(2:end)./sqrt(dimensions_cell{i_type}(2:end));
        y_scale = [0 0.03];
    case {'minalien', 'prat'}
        fs_chars_cell{i_type} = 1-in_data.fs_chars;
        y_scale = [0 0.2];
    otherwise 
        warning('Не задан метод визуализации');
    end
    conf_mx_cell{i_type} = in_data.conf_mx;
    obj_names_cell{i_type} = in_data.obj_names;
end, clear i_type t_type in_filename in_data

%% Visualize characterstics
marker_type = {  '*' '+' 'x' 's' 'd' 'o' '.'};

figure('Name',reduction_method,'color', 'white','WindowStyle','docked'); 
horz_grid = n_types+mod(n_types,2);
subplot(2,horz_grid,1:horz_grid/2); hold on;

curve_h = zeros(1,n_types);


for i_type = 1:n_types
    curve_h(i_type) = stem(dimensions_cell{i_type},fs_chars_cell{i_type},...
        [ ':' marker_type{i_type} 'k']);
    plot(dimensions_cell{i_type},fs_chars_cell{i_type},':');
end
xticks(dimensions_cell{i_type})
% set(gca, 'YScale', 'log');
% ylim(y_scale)
ylabel(reduction_method)
legend(curve_h,t_types)
title(['Граничные метрики (выборка ' reduction_method ...
    ', \tau_и=' num2str(imp_hwidth_ns) ...
    'нс, СКО=' num2str(mse_noise) ')'])

subplot(2,horz_grid,horz_grid/2+1:horz_grid); hold on;
for i_type = 1:n_types
    curve_h(i_type) = stem(dimensions_cell{i_type},fs_perf_cell{i_type},...
        [ ':' marker_type{i_type} 'k']);
    plot(dimensions_cell{i_type},fs_perf_cell{i_type},':');
end
xticks(dimensions_cell{i_type})
% set(gca, 'YScale', 'log');
ylim([0 0.2])
legend(curve_h,t_types)
title(['Ошибки (выборка ' reduction_method ...
    ', \tau_и=' num2str(imp_hwidth_ns) ...
    'нс, СКО=' num2str(mse_noise) ')'])

%% Visualize sparse matrices
for i_type = 1:n_types
%     dim_vect = 1:min(max_dim,length(fs_chars_cell{i_type}));
    feature_names = strsplit(strrep(feat_header_cell{i_type},'_',' '),',');
    
    subplot(2,horz_grid,horz_grid+i_type); hold on;
    spy(fs_maps_cell{i_type}.',[marker_type{i_type} 'k']);
    xticks(1:length(dimensions_cell{i_type}))
    xticklabels(dimensions_cell{i_type})
    xlabel('dimensionality, n')
    ylabel([t_types{i_type} ' features'])
    title(['Feature selection for case ' t_names{i_type} ' transform'])
    ylim([0 size(fs_maps_cell{i_type},2)+1])
    yticks(1:size(fs_maps_cell{i_type},2))
    yticklabels(feature_names)
end

%% Visualize confusion matricies
% figure('Name','','color', 'white','WindowStyle','docked');
% for i_type = 1:n_types
%    
%     subplot(n_types, 1, i_type)
%     cm = confusionchart(conf_mx_cell{1},obj_names_cell{i_type});
%     cm.RowSummary = 'row-normalized';
%     cm.ColumnSummary = 'column-normalized';
% end

% end