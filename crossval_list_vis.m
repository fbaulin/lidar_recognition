function crossval_list_vis(filelist,varargin)
%Сценарий визуализации результатов вычиселния mhist
    kwargs = KeywordArguments(...
            'maps',1,...
            'confusion',1 );
    [   maps_flag, confusion_mx_flag ] = kwargs.parse_input_cell(varargin);
    n_files = length(filelist);
    compare_params =    {'t_type' 't_bfun' 'imp_hwidth_ns' 'mse_noise' 'reduction_method'};
    compare_params_tex ={''       ''       '\tau_n_s='         '\sigma_n='   'Q - '};
    datastruct = extract_data(filelist);

    marker_type = { '*' '+' 'x' 's' 'o' '^' 'v' 'd' 'p' 'h' '.'};
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

    [leg_str, tit_str] = compile_legend(datastruct, compare_params, compare_params_tex);
    draw_fs_score(datastruct, leg_str, tit_str, ax_fs_score, marker_type);
    draw_clf_score(datastruct, leg_str, tit_str, ax_clf_score, marker_type);
    if maps_flag,   draw_fs_maps(datastruct, ax_fs_maps, leg_str, marker_type); end
    if confusion_mx_flag,   confusion_mx(datastruct,leg_str); end
end
%% Служебные функции
% Visualize fs score
function draw_fs_score(datastruct, leg_str, tit_str, varargin)
    if nargin>=4, axes(varargin{1}); end
    if nargin>=5, marker_type = varargin{4};
    else, marker_type = { '*' '+' 'x' 's' 'o' '^' 'v' 'd' 'p' 'h' '.'};
    end
    n_files = size(datastruct,2); 
    hold on;
    curve_h = zeros(1,n_files);
    for i_file = 1:n_files
        curve_h(i_file) = stem(datastruct(i_file).('dimensions'),datastruct(i_file).('fs_chars'),...
            [ ':' marker_type{i_file} 'k']);
        plot(datastruct(i_file).('dimensions'),datastruct(i_file).('fs_chars'), ':');
    end
    xticks(datastruct(i_file).('dimensions'))
    ylabel(datastruct(i_file).('reduction_method'))
    legend(curve_h, leg_str) 
    title(['Критерий отбора Q (' tit_str ')'])
end

% Visualize clf score
function draw_clf_score(datastruct, leg_str, tit_str, varargin)
    if nargin>=4, axes(varargin{1}); end
    if nargin>=5, marker_type = varargin{4};
    else, marker_type = { '*' '+' 'x' 's' 'o' '^' 'v' 'd' 'p' 'h' '.'};
    end
    hold on;
    n_files = size(datastruct,2); 
    curve_h = zeros(1,n_files);
    
    for i_file = 1:n_files
        curve_h(i_file) = stem(...
            datastruct(i_file).('dimensions'),datastruct(i_file).('fs_perf'),...
            [ ':' marker_type{i_file} 'k']);
        plot(datastruct(i_file).('dimensions'),datastruct(i_file).('fs_perf'),':');
    end
    xticks(datastruct(i_file).('dimensions'))
    legend(curve_h,leg_str)
    title(['Ошибки (' tit_str ')'])
end

% Visualize fs maps
function draw_fs_maps(datastruct, ax, leg_str, marker_type)
    n_files = size(datastruct,2); 
    for i_file = 1:n_files
        feature_names = strsplit(strrep(datastruct(i_file).('feat_header'),'_',' '),',');
        axes(ax(i_file)); 
        hold on;
        spy(datastruct(i_file).('fs_maps').',[marker_type{i_file} 'k']);
        xticks(1:length(datastruct(i_file).('dimensions')))
        xticklabels(datastruct(i_file).('dimensions'))
        xlabel('dimensionality, n')
        ylabel([datastruct(i_file).('t_type') ' features'])
        title(['Признаки для ' leg_str(i_file)])
        ylim([0 size(datastruct(i_file).('fs_maps'),2)+1])
        yticks(1:size(datastruct(i_file).('fs_maps'),2))
        yticklabels(feature_names)
    end
end

% Visualize confusion matricies
function confusion_mx(datastruct,leg_str)
    n_files = size(datastruct,2);
    figure('Name','','color', 'white','WindowStyle','docked');
    for i_type = 1:n_files
        subplot(n_files, 1, i_type)
        cm = confusionchart(datastruct(i_type).('conf_mx'),...
            datastruct(i_type).('obj_names'));
        cm.RowSummary = 'row-normalized';
        cm.ColumnSummary = 'column-normalized';
        title(leg_str(i_type))
    end
end

% Extraxt data
function datastruct = extract_data(filelist)
    
    n_files = length(filelist);
    datastruct = cell(1,n_files);
    
    for i_file = 1:n_files
        datastruct{i_file} = load(filelist{i_file});
        datastruct{i_file}.('dimensions') = sum(datastruct{i_file}.('fs_maps'),2);
        
        crop_ind = strfind(filelist{i_file},'_t');
        file_name = filelist{i_file};
        datastruct{i_file}.('filename') = file_name(crop_ind(end)+1:end-4);
        
        switch datastruct{i_file}.('reduction_method')
            case {'nmin', 'buryi'}
                datastruct{i_file}.('fs_chars')(2:end) = datastruct{i_file}.('fs_chars')(2:end)./sqrt(datastruct{i_file}.('dimensions')(2:end));
            case {'minalien', 'prat'}
                datastruct{i_file}.('fs_chars') = 1-datastruct{i_file}.('fs_chars');
            case 'fisher'
            otherwise, warning('Не задан метод визуализации');
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
            title_str = [title_str,tex_names{i_param},c{:},','];
        end
        
    end
    title_str(end)=[];
    for i=1:length(legend_cell), legend_cell{i}(end)=[]; end
    if isempty(horzcat(legend_cell{:}))
        legend_cell = cellfun(@(s) strrep(s,'_',','),...
            {dataset.('filename')},'UniformOutput',false);
    end
    
end