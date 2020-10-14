classdef RecursiveReduction < handle
    %RecursiveReduction ������������� ����������� ������ ������ ���������
    %   ������ ������ ��������� ������������� �� ������ ������������ ������� ����������� ��������� 
    
    properties
        rhos_sqr        % ������� �������
        feat_map        % ������ ��������� ������������ ���������
        fvs             % ������� ���������
        targets_ids     % ����������������� ������, ������������� ������ ����� ������
        alien_mask
        weight_fun      % ������������ �������
        
        use_gpu         % ������������ gpu
        reduction       % ����� ��������
        heurist         % ������������� ��������
        estimate_q      % ������� ������ �������� ������������
        
    end
    
    properties(Hidden=true)
        upd=true    % �������� ������� �� ������� �������
        gpu_supported_methods = { 'test', 'nmin', 'buryi', 'buryisum', 'refmmin', 'integrnmin', 'integrnmax', 'integr', 'auhist', 'integr_min', 'minalien', 'prat'}
    end
    
    methods
        function obj = RecursiveReduction(varargin)
        %RecursiveReduction ������������ ��� ������� ���������� ��������
        %   Detailed explanation goes here
        %   ����������� ���������:
        %   'gpu_on' - 'True'/'False' - ������������ ����������� �������
        %   'reduction' - ��� ��������
        %       'buryi'
        %       'nmin'
        %       'mhist'
        %       'minalien'
        %       'prat'
        %       'wprat'
        %       'fisher'
        %       'refmin'
        %   'heuristic' - ������������� ����� ���������� ����������
        %       'adddel' - ����� ������������ ���������� � �������� ��������� (����������� �� ���
        %       ������������ ���������)
        %       'gray'  - ����������� ����������� ����� ����
        %   'dimensions' - ����������� ������������ ��������� ��� ������� ����� ����� �������
        %   �������������� �������:
        %   'n_metrics_fraction'
        %   'm_technique'
        %   'k_alien'
        %   'n_nearest'
        %   'hist_edges'
            
            kw = KeywordArguments();
            obj.use_gpu = kw.get_value(varargin, 'gpu_on', gpuDeviceCount());
            obj.heurist = kw.get_value(varargin, 'heuristic', 'adddel');          % ��������� ������������
            obj.reduction = kw.get_value(varargin, 'reduction_type', 'buryi');   % ����������� �����������
            if or(strcmp(obj.heurist, 'fisher'), strcmp(obj.reduction, 'fisher'))
                obj.reduction = 'fisher'; obj.heurist = 'fisher';
            end
            obj.setup_reduction(varargin);
    
        end
        
        function [fspace_qs, fspace_map, varargout] = reduce(obj, fvs, targets, varargin)
            kw = KeywordArguments();
            obj.use_gpu = and(...
                any(strcmp(kw.get_value(varargin, 'reduction_type'), obj.gpu_supported_methods)), ...
                kw.get_value(varargin, 'gpu_on', gpuDeviceCount()) );    % �������� ���� �� ���� ������� ����� �����
            obj.heurist = kw.get_value(varargin, 'heuristic', obj.heurist);    
            obj.write_data(fvs, targets);    % ������ ������
            obj.setup_reduction(varargin{:});   % ��������� ��������
            
            dimensions = kw.get_value(varargin,'dimensions',1:size(fvs,2));
            dimensions = dimensions(dimensions<=size(fvs,2)); % �������� ������ ������������
            
            fprintf('(%s) H.=%s, ',datetime('now','Format','HH:mm:ss'), obj.heurist);
            switch obj.heurist
                case 'adddel',              [fspace_qs, fspace_map] = ...
                        obj.add_del_reduction(obj.estimate_q, dimensions);
                case 'gray',                [fspace_qs, fspace_map] = ...
                        obj.gray_reduction(obj.estimate_q, dimensions);
                case 'fisher',              [fspace_qs, fspace_map] = ...
                        obj.fisher_selection(dimensions);
                otherwise, error('������ ������ �������������� ��������� ��������')
            end
            fprintf('%s: [ %s ]\n',obj.reduction,num2str(fspace_qs.', '%6.3f '));
            if nargout==3, varargout{1} = dimensions; end
            
        end
        
        function write_data(obj, fvs, targets)
            dim_fv = size(fvs,2);            % ����������� ��
            obj.feat_map = zeros(1,dim_fv,'single');
            [~, ~, obj.targets_ids] =  unique(targets); 
            if obj.use_gpu
                obj.fvs = gpuArray(fvs);
                obj.targets_ids = gpuArray(uint8(obj.targets_ids));
                obj.alien_mask = gpuArray(and(...
                    obj.targets_ids~=obj.targets_ids.',...
                    triu(ones(size(obj.targets_ids,1),'like',obj.targets_ids),1)));     % ����� ������� ������                        
            else
                obj.fvs = fvs;
                obj.targets_ids = uint8(obj.targets_ids);
                obj.alien_mask = and(...
                    obj.targets_ids~=obj.targets_ids.',...
                    triu(ones(size(obj.targets_ids,1),'logical'),1));   % ����� ������� ������
            end            
            obj.rhos_sqr = zeros(size(fvs,1),'like',fvs);
        end
        
        function setup_reduction(obj, varargin)
            kw = KeywordArguments();
            obj.reduction = kw.get_value(varargin, 'reduction_type', obj.reduction);      % ����������� �����������
            n_metrics_fraction = kw.get_value(varargin, 'n_metrics_fraction', 0.05);      % ����� ������
%             hist_edges =    kw.get_value(varargin, 'hist_edges', linspace(0,2,51));  % ������� ���������� �����������
            n_nearest =     kw.get_value(varargin, 'n_nearest', 5);               % ����� ����������� �������
            m_technique =   kw.get_value(varargin, 'm_technique','mean');            % ����� ������
            k_alien =       kw.get_value(varargin, 'k_alien', 2);                   % ����� �����
            switch obj.reduction
                case 'test',        rho_estimate = @(t_map) obj.local_metric_sum(t_map, n_nearest);
                case 'nmin',        rho_estimate = @(t_map) obj.nmin_metric_balanced(t_map, n_metrics_fraction);
                case 'integrnmin',  rho_estimate = @(t_map) obj.integr_nmin(t_map, n_metrics_fraction);
                case 'integrnmax',  rho_estimate = @(t_map) obj.integr_nmax(t_map, n_metrics_fraction);
                case 'integr',      rho_estimate = @(t_map) obj.integr(t_map);
                case 'integr_min',  rho_estimate = @(t_map) obj.integr_min(t_map);
                case 'buryi',       rho_estimate = @(t_map) obj.buryi(t_map);
                case 'buryisum',       rho_estimate = @(t_map) obj.buryi_sum(t_map);
%                 case 'mhist',   rho_estimate = @(t_map) obj.hist_acc(t_map, hist_edges, n_metrics_fraction);
                case 'minalien',    rho_estimate = @(t_map) obj.minalien(t_map, n_nearest, m_technique, k_alien);
                case 'auhist',      rho_estimate = @(t_map) obj.auhist(t_map);
                case 'prat',        rho_estimate = @(t_map) obj.max_p_relative(t_map, n_nearest);
                case 'wprat'   
                    rho_estimate = @(t_map) obj.w_p_relative(t_map, n_nearest);
                    obj.weight_fun = ...
                        linspace(1,0,n_nearest+1 );
                case 'refmmin'
                    obj.refine_aliens(k_alien, n_nearest)
                    rho_estimate = @(t_map) obj.buryi(t_map);
                case 'fisher',      obj.heurist = 'fisher';
                
                otherwise,      error('����������� ����� ����� ��������. ��������������: nmin, buryi, mhist, minalien')
            end
            if obj.use_gpu
                if any(strcmp(obj.reduction,obj.gpu_supported_methods))
                    obj.estimate_q = @(t_map) gather(rho_estimate(t_map));
                else
                    obj.use_gpu = false;
                    obj.estimate_q = rho_estimate;
                end
            end
        end
        
        %% ������������� ��������� ��������
        
        % �������� �� ������ ����������� ������� ��������� ����� ����
        function [fspace_qs, fspace_map] = gray_reduction(obj, rho_estimate, dimensions)
            dim_fv = size(obj.fvs,2);
            n_spaces = length(dimensions);
            fspace_qs = zeros(n_spaces,1);          % �������� ��������� ������� ��� ������������
            fspace_map = zeros(n_spaces, dim_fv);   % ������� ����������� ����������� - �� ������� �� �����������
            n_spaces = 2^dim_fv-1;
            
            line_len = fprintf('�����������: %d. \t���������: 0.00%%', n_spaces);
            cur_percent = 0;
            obj.upd = true;
            for feat_val = 1:n_spaces               % �� ���� �������� ������������� ���������
                t_feat_map = de2bi(feat_val,dim_fv,'left-msb'); % ��������� ����� �������� ������������
                t_feat_map_sh = circshift(t_feat_map,1);
                t_feat_map_sh(1) = 0;
                t_feat_map = xor(t_feat_map,t_feat_map_sh);     % ������������� ���� ���� ��� ����������� ��������
                i_fmap = find(sum(t_feat_map)==dimensions); % ������ ��������������� �����������
                if i_fmap
                    [cur_q_val] = rho_estimate(t_feat_map);
                    if cur_q_val>=fspace_qs(i_fmap)     % ���� �������� ������� ������ �������� ��� ��. ��������� ����� �����������
                        fspace_qs(i_fmap) = cur_q_val;
                        fspace_map(i_fmap,:) = t_feat_map(:);
                    end
                    if feat_val/n_spaces>=cur_percent
                        cur_percent = cur_percent+0.001;
                        fprintf('\b\b\b\b\b\b%05.1f%%',cur_percent*100);
                    end
                end
            end
            fprintf(repmat('\b', 1, line_len));
        end
        
        % �������� �� ������ ��������� ������� ������
        function [fspace_qs, fspace_map] = fisher_selection(obj, dimensions)
        % ������������ ����� ������� � �������� ���������
            dim_fv = size(obj.fvs,2);       % ����������� ���. �������
            n_spaces = length(dimensions);  % ����� ��������������
            fspace_qs = zeros(n_spaces,1);  % �������� ��������� ������� ��� ������������
            fspace_map = zeros(n_spaces, length(obj.feat_map));   % ������� ����������� ����������� - �� ������� �� �����������
            
            m_exp_i     = mean(obj.fvs,1);          % ����������� ���� �������
            classes_ids = unique(obj.targets_ids);    % �������������� �������
            n_classes   = length(classes_ids);                % ����� �������
            m_exp_ij    = zeros(dim_fv,n_classes);  % ������� ����������� �������
            var_ij      = zeros(dim_fv,n_classes);  % ������� ��� �������
            fisher_score = zeros(dim_fv,1);         % ������ ������ ������
            n_j         = zeros(1, n_classes);      % ����� �������� � �������
            % ������ ������ ������
            for i_class=1:n_classes
                n_j(i_class) = sum(obj.targets_ids==classes_ids(i_class));     % ����� �������� ������
                m_exp_ij(:,i_class) = mean(obj.fvs(obj.targets_ids==classes_ids(i_class),:),1);   % ����������� �������
                var_ij(:,i_class) = var(obj.fvs(obj.targets_ids==classes_ids(i_class),:),0,1);    % ��������� � ���������� �� N-1
            end
            for i_feat = 1:dim_fv             % ���� �� ���������
                fisher_score(i_feat) = sum(n_j.*(m_exp_ij(i_feat,:)-m_exp_i(i_feat)).^2)/sum(n_j.*var_ij(i_feat));
            end
            % ��������� ��������
            [fisher_score, fisher_ids] = sort(fisher_score,'descend');
            for i_dim = 1:n_spaces
                curr_dmty = dimensions(i_dim); % ������� �����������
                fspace_qs(i_dim) = sum(fisher_score(1:curr_dmty));
                fspace_map(i_dim, fisher_ids(1:curr_dmty))=1;    % ������������� �������� �� ����                
            end
        end
        
        % �������� ������������ �������� � ���������� ���������
        function [fspace_qs, fspace_map] = add_del_reduction(obj, rho_estimate, dimensions)
            
            n_spaces = length(dimensions);
            fspace_qs = zeros(n_spaces,1);          % �������� ��������� ������� ��� ������������
            fspace_map = zeros(n_spaces, length(obj.feat_map));   % ������� ����������� ����������� - �� ������� �� �����������
            
            obj.upd = false;
            cur_feat_map = obj.feat_map;
            alg_steps = [1 1 1 -1 -1];  % ��� ��������� 1 - ��� ������, -1 - ��� �����
            
            line_len = fprintf('������(%02d): 00',max(dimensions));
            for i_dim = dimensions(dimensions>0) %TODO: ���� �� ���������������� ����������� - ����� ������
                % ������� 3 ���� ������ � ��� - �����
                fprintf('\b\b%02d',i_dim);
                for step=(alg_steps==-1)
                    % ������� ���� ����������� i_dim+1
                    z_idxs = find(cur_feat_map==step); % �������� ������ � ��������� ������� ��������� (m ����)
                    if or(...
                            isempty(z_idxs), ...    % ���� ���������� ������ ���
                            and(sum(cur_feat_map)<(i_dim+1),step)...  % ������� ������������ ��� ��������
                            )
                        continue;  
                    end
                    q_sep_vect = zeros(length(z_idxs),1);
                    for i_idx = 1:length(z_idxs)            % ��������� m ���������
                        t_feat_map = cur_feat_map;          % ����������� �����
                        t_feat_map(z_idxs(i_idx)) = ~step;
                        q_sep_vect(i_idx) = rho_estimate(t_feat_map);   %
                    end % ���������
                    [~, best_idx] = max(q_sep_vect);         % ����� id �������� ��������
                    cur_feat_map(z_idxs(best_idx)) = ~step; % �������� ������ ������������
                end
                obj.upd=true; 
                dim_idx = find(dimensions==i_dim);
                fspace_qs(dim_idx) = rho_estimate(cur_feat_map); 
                fspace_map(dim_idx,:) = cur_feat_map;
                obj.upd=false; % �������� ������������
            end
            fprintf(repmat('\b', 1, line_len));
            
        end
           
        
        
        %% ���� �������� ������������
        
        % ������ ���������������� ����� �����
        function space_q = minalien(obj, t_map, n_nearest, m_technique, k_alien)            
            
            n_samples = size(obj.fvs, 1);
            [~, sort_indexes] = mink(obj.get_metrics(t_map), n_nearest+1, 2);  % �������� ������� ������� ��� ������� �������
            %FIXME: � ��������� ������� ��������� ����� �� ���� ����� �������� (����� �� ������ �������� ���� 0)
            neq_matrix = zeros(n_samples, n_nearest); % ������������ �-�� �����-��� 
            
            for i=1:n_nearest 
                 neq_matrix(:, i) = ...
                     obj.targets_ids ~= obj.targets_ids(sort_indexes(:, i+1)); % i+1, �.�. 1� - ��� ������ ���������
            end
            switch m_technique
                case 'mean'     % ������� ����� ����� ������� ������ ����� �������
                    q_sep = mean(sum(neq_matrix, 2)/n_nearest);
                case 'thresh'   % ���� �������� ��� ������� ����� ��������� ������� ������ ������
                    q_sep = sum(sum(neq_matrix, 2)>(k_alien))/n_samples;
            end
            space_q = 1 - q_sep;
        end
        
        function space_q = local_metric_sum(obj, t_map, n_nearest)            
        %local_metric_sum ��������� ������ ����� ����������� ������ � ��������� ���������� ��������
        %NOTE: ������������� ������������
            n_samples = size(obj.fvs, 1);
            rhos = obj.get_metrics(t_map);                                  % ������� ������
            rhos(1:n_samples+1:end)=nan;
            least_rhos = zeros(n_samples, n_nearest,'like', obj.fvs);       % ������� ���������� ������
            neq_matrix = false(n_samples, n_nearest,'like', obj.alien_mask);
            for i=1:n_nearest 
                [least_rhos(:,i), sort_indexes] = min(rhos,[], 2);                  % ����������� �������
                rhos(sub2ind(size(rhos), (1:n_samples).', sort_indexes)) = nan;
                neq_matrix(:,i) = ...
                    obj.targets_ids == obj.targets_ids(sort_indexes); % i+1, �.�. 1� - ��� ������ ���������
            end
            
            space_q = sum(least_rhos(neq_matrix),'all')-sum(least_rhos(~neq_matrix),'all');
        end
        
        function space_q = nmin_metric_balanced(obj, t_map, n_min_ratio)
        %nmin_metric ��������� ������ ������������ ������������ ��������� �� �������� �����������
        %   fvs_red     - ������ ����� � ��������� ��������� 
        %   n_min_ratio - ���� ������, �� ������� ��������� �������� ������� (<1)
        %   w_rep       - ������������ ����������
            rhos = obj.get_metrics(t_map);              % ������ ������
            rhos = sort( rhos(obj.alien_mask), 'ascend');  % ��������� ������� ������ � �������������� ����������� �������
            n_rho_index = floor(length(rhos)*n_min_ratio);
            space_q = rhos(n_rho_index);        % ����� n-��� �������
        end
        
        function space_q = integr_nmin(obj, t_map, n_min_ratio)
            rhos = obj.get_metrics(t_map);
            rhos = sort(rhos(obj.alien_mask), 'ascend');
%             if obj.upd, Loggers.log('integr_16db.log',{sum(t_map),rhos}, 'sep', ','); end
            n_rho_index = floor(length(rhos)*n_min_ratio);
            space_q = sum(rhos(1:n_rho_index));        % ����� n-��� �������
        end
        
        function space_q = integr_min(obj, t_map)
        % ������� ����� ������ ����� ���������� �������� ������� ������
        %NOTE: ������������� ������������
            rhos = obj.get_metrics(t_map);
            rhos(~or(obj.alien_mask,obj.alien_mask.')) = NaN;
            rhos = min(rhos,[],2);
            space_q = sum(rhos);                % 
        end

        function space_q = integr_nmax(obj, t_map, n_min_ratio)
            rhos = obj.get_metrics(t_map);
            rhos = sort(rhos(obj.alien_mask), 'ascend');
%             if obj.upd, Loggers.log('integr_16db.log',{sum(t_map),rhos}, 'sep', ','); end
            n_rho_index = floor(length(rhos)*n_min_ratio);
            space_q = sum(rhos(n_rho_index:end));        % ����� n-��� �������
        end
        
        function space_q = integr(obj, t_map)
            rhos = obj.get_metrics(t_map);
            space_q = sum(rhos(obj.alien_mask));
        end
        
        function space_q = auhist(obj, t_map)
            % ���� ���������� �� ���������� �������� ����� ������ � ���� AUC ROC
            rhos = obj.get_metrics(t_map);
            rhos = rhos(obj.alien_mask);
            sum_rhos = sum(rhos);
            max_rhos = max(rhos);
            space_q = sum_rhos*sum_rhos/(length(rhos)*max_rhos);
        end 
        
        function space_q = buryi_sum(obj, t_map)
            rhos = obj.get_metrics(t_map);      % ������ ������
%             n_samples = size(obj.fvs, 1);       % ����� ��������
%             if obj.upd, Loggers.log('.\metric_log\buryi_16db.log',{sum(t_map),gather(sort(rhos(obj.alien_mask)))}, 'sep', ',','header',''); end
            space_q = sum(rhos, 'all');
            % ��������� ������� ������ � �������������� ����������� �������
        end
        
        function space_q = buryi(obj, t_map)
            rhos = obj.get_metrics(t_map);      % ������ ������
%             n_samples = size(obj.fvs, 1);       % ����� ��������
%             if obj.upd, Loggers.log('.\metric_log\buryi_16db.log',{sum(t_map),gather(sort(rhos(obj.alien_mask)))}, 'sep', ',','header',''); end
            space_q = min(rhos(obj.alien_mask));
            % ��������� ������� ������ � �������������� ����������� �������
        end
        
        % �������� ������� ��������� �����������
        function space_q = max_p_relative(obj, t_map, n_nearest)
            [~, local_groups] = mink(obj.get_metrics(t_map),n_nearest, 2);  % �������� ������� ������� ��� ������� �������
            local_groups = obj.targets_ids(local_groups);   % � ������� -> ������
            [~, occur] = mode(local_groups,2);              % ������������� ����� � ������
            space_q = mean(occur)/n_nearest;                % 
%             group_modes = mode(local_groups,2);                 % ������������� ����� � ������
%             space_q = mean(sum(group_modes==local_groups,2))/n_nearest; % ������� ����� ������.   
        end
        
       
        % �������� ������� ��������� �����������
        function space_q = w_p_relative(obj, t_map, n_nearest)
            [~, local_groups] = mink(obj.get_metrics(t_map),n_nearest, 2);  % �������� ������� ������� ��� ������� �������
            local_groups = obj.targets_ids(local_groups);       % � ������� -> ������
            group_modes = mode(local_groups,2);                 % ������������� ����� � ������
            space_q = 1-mean(sum(obj.weight_fun.*(group_modes~=local_groups),2))/n_nearest; % ���������� ����� �����   
        end
        
        %% ������ ��������� �������
        
        % ������ ������
        function [ rhos ] = get_metrics( obj, t_map )
        %get_metrics ������� ������� ������� ����� ��������� ���� �������
        %   ������� ��������� ��� ������ � m(C1Data) � n(C2Data) ���������.
        %   ������������ �������� ����� ����� ��������
        %   �� ������ mxn ������� ���������� ������ 
            cor_map = t_map - obj.feat_map; % ������ �������: ������-1, �����+1
            if and(nnz(cor_map) >= nnz(t_map), obj.upd)   % ���� ����� ��������� �� ������ ����� ��������� �� ��������� ������
                obj.rhos_sqr(:) = 0;        % �������� �������
                obj.feat_map(:) = 0;        % �������� �����
                cor_map = t_map;            % ����� ������������� = ����� �������� ������������ ���������
            end
            fvs_red = obj.fvs(:,cor_map~=0);    % �������� ���������� ��, ������������ ��� ���������
            metr_dif = RecursiveReduction.get_correction_components(fvs_red,fvs_red).^2; % ��������� �������� �������������
            cor_map(cor_map==0) = [];       % ������ ���������� ��������
            metr_dif(cor_map.'<0,:,:) = -metr_dif(cor_map.'<0,:,:);
            metr_dif = squeeze(sum(metr_dif, 1));

            if obj.upd  % ���� ����� ������� ����������, �� ��������
                obj.rhos_sqr = obj.rhos_sqr + metr_dif;
                obj.feat_map = t_map;
                rhos = obj.rhos_sqr;
%                 rhos = sqrt(obj.rhos_sqr);
            else
%                 rhos = sqrt(obj.rhos_sqr + metr_dif);
                rhos = obj.rhos_sqr + metr_dif;
            end
        end
        
        
        % ������� ������� �������� �����������
        function [] = refine_aliens(obj, k_alien, n_nearest)
            
            % ���������� �������
            % �����������
            % �������            
            update_status = obj.upd;
            obj.upd = false;
            [n_samples, dim_fv] = size(obj.fvs);
            t_map = ones(1,dim_fv,'single');
            [~, sort_indexes] = mink(obj.get_metrics(t_map), n_nearest+1, 2);  % �������� ������� ������� ��� ������� �������
            neq_matrix = zeros(n_samples, n_nearest, 'like', obj.targets_ids); % ������������ �-�� �����-���
            
            for i=1:n_nearest 
                 neq_matrix(:, i) = ...
                     obj.targets_ids ~= obj.targets_ids(sort_indexes(:, i+1)); % i+1, �.�. 1� - ��� ������ ���������
            end
            % �������� ������� ��������������� n_samples x n_nearest
            del_idx = sum(neq_matrix, 2)>(k_alien);
            
            obj.fvs(del_idx,:)=[];
            obj.rhos_sqr = zeros(size(obj.fvs,1),'like',obj.fvs);   % ������� �������
            obj.targets_ids(del_idx,:)=[];
            obj.alien_mask = and(...
                obj.targets_ids~=obj.targets_ids.',...
                triu(ones(size(obj.targets_ids,1),'like',obj.targets_ids),1)); % ����� ������� ������
            obj.upd = update_status;
        end
       
        
        %% legacy
        % �������� ������������ ���������
        function [varargout] = fs_reduction(obj, fvs, targets, varargin)
            varargout=cell(1,nargout);
            
            reduction_name = KeywordArguments.get_value(varargin, 'reduction_type');
            supported = any(strcmp(reduction_name,obj.gpu_supported_methods));
            if and(gpuDeviceCount() >=1, supported)
                [varargout{:}] = fs_reduction_gpu(obj, fvs, targets, varargin{:});
            else
                [varargout{:}] = fs_reduction_cpu(obj, fvs, targets, varargin{:});
            end
           
        end
        % �������� ������������ ��������� � ������ ��
        function [fspace_qs, fspace_map, varargout] = fs_reduction_cpu(obj, fvs, targets, varargin)
        %fs_reduction �������� � ������ ����������� �������
        % �� ������
        % fspace_rhos - ������-������� � ���������
        % fspace_maps - ������� �� ������� ������� ����������� ��������� ������������
        % ���������: 
        % fvs - ������� ��������� (�� �������)
        % ����������� ���������:
        % dimensions, reduction_type, n_metrics_fraction, hist_edges,
        % obj_weight
            obj.fvs = fvs;
            dim_fv = size(fvs,2);            % ����������� ��
            obj.feat_map = zeros(1,dim_fv,'single');
            obj.rhos_sqr = zeros(size(fvs,1),'like',fvs);
            [~, ~, obj.targets_ids] = unique(targets);        % ������������� ����� � id-������
            obj.alien_mask = and(...
                obj.targets_ids~=obj.targets_ids.',...
                triu(ones(size(obj.targets_ids,1),'logical'),1)); % ����� ������� ������             
            kwargs = KeywordArguments(...
                'dimensions',1:dim_fv, ...          % ����������� �����������
                'reduction_type', 'fisher',...      % ����������� �����������
                'n_metrics_fraction', 0.05,...      % ����� ������
                'hist_edges', linspace(0,2,51),...  % ������� ���������� �����������
                'n_nearest', 5, ...                 % ����� ����������� �������
                'm_technique','mean',...            % ����� ������
                'k_alien', 2, ...                   % ����� �����
                'heuristic', 'adddel');               % ��������� ������������
            [dimensions, reduction_type, n_metrics_fraction, hist_edges, ...
                n_nearest, m_technique, k_alien, heuristic] =  ...
                kwargs.parse_input_cell(varargin);
            % ��������� ������������� ��� ������
            switch reduction_type
                case 'nmin',    rho_estimate = @(t_map) obj.nmin_metric_balanced(t_map, n_metrics_fraction);
                case 'buryi',   rho_estimate = @(t_map) obj.buryi(t_map);
                case 'mhist',   rho_estimate = @(t_map) obj.hist_acc(t_map, hist_edges, n_metrics_fraction);
                case 'minalien',rho_estimate = @(t_map) obj.minalien(t_map, n_nearest, m_technique, k_alien);
                case 'prat',    rho_estimate = @(t_map) obj.max_p_relative(t_map, n_nearest);
                case 'wprat'   
                    rho_estimate = @(t_map) obj.w_p_relative(t_map, n_nearest);
                    obj.weight_fun = ...
                        linspace(1,0,n_nearest+1 );
                case 'fisher',  heuristic = 'fisher';
                case 'refmmin'
                    obj.refine_aliens(k_alien, n_nearest)
                    rho_estimate = @(t_map) obj.buryi(t_map);
                otherwise,      error('����������� ����� ����� ��������. ��������������: nmin, buryi, mhist, minalien')
            end
            dimensions = dimensions(dimensions<=dim_fv); % �������� ������ ������������
            fprintf('(%s) H.=%s, ',datetime('now','Format','HH:mm:ss'), heuristic);
            switch heuristic
                case 'adddel',              [fspace_qs, fspace_map] = ...
                        obj.add_del_reduction(rho_estimate, dimensions);
                case 'gray',                [fspace_qs, fspace_map] = ...
                        obj.gray_reduction(rho_estimate, dimensions);
                case 'fisher',              [fspace_qs, fspace_map] = ...
                        obj.fisher_selection(dimensions);
                otherwise, error('������ ������ �������������� ��������� ��������')
            end
            fprintf('%s: [ %s ]\n',reduction_type,num2str(fspace_qs.', '%6.3f '));
            if nargout==3, varargout{1} = dimensions; end
        end
        
        % �������� ������������ ��������� � ������ ��
        function [fspace_qs, fspace_map, varargout] = fs_reduction_gpu(obj, fvs, targets, varargin)
        %fs_reduction �������� � ������ ����������� �������
        % �� ������
        % fspace_rhos - ������-������� � ���������
        % fspace_maps - ������� �� ������� ������� ����������� ��������� ������������
        % ���������: 
        % fvs - ������� ��������� (�� �������)
        % ����������� ���������:
        % dimensions, reduction_type, n_metrics_fraction, hist_edges,
        % obj_weight
            obj.fvs = gpuArray(fvs);
            dim_fv = size(fvs,2);            % ����������� ��
            obj.feat_map = zeros(1,dim_fv,'single');
            obj.rhos_sqr = zeros(size(fvs,1),'single','gpuArray');
            [~, ~, obj.targets_ids] =  unique(targets);        % ������������� ����� � id-������
            obj.targets_ids = gpuArray(uint8(obj.targets_ids));
            obj.alien_mask = gpuArray(and(...
                obj.targets_ids~=obj.targets_ids.',...
                triu(ones(size(obj.targets_ids,1),'like',obj.targets_ids),1))); % ����� ������� ������
                        
            kwargs = KeywordArguments(...
                'dimensions',1:dim_fv, ...          % ����������� �����������
                'reduction_type', 'fisher',...      % ����������� �����������
                'n_metrics_fraction', 0.05,...      % ����� ������
                'hist_edges', linspace(0,2,51),...  % ������� ���������� �����������
                'n_nearest', 5, ...                 % ����� ����������� �������
                'm_technique','mean',...            % ����� ������
                'k_alien', 2, ...                   % ����� �����
                'heuristic', 'adddel');               % ��������� ������������
            [dimensions, reduction_type, n_metrics_fraction, hist_edges, ...
                n_nearest, m_technique, k_alien, heuristic] =  ...
                kwargs.parse_input_cell(varargin);
            % ��������� ������������� ��� ������
            switch reduction_type
                case 'nmin',    rho_estimate = @(t_map) gather(obj.nmin_metric_balanced(t_map, n_metrics_fraction));
                case 'buryi',   rho_estimate = @(t_map) gather(obj.buryi(t_map));
                case 'minalien',rho_estimate = @(t_map) gather(obj.minalien(t_map, n_nearest, m_technique, k_alien));
                case 'mhist',   rho_estimate = @(t_map) gather(obj.hist_acc(t_map, hist_edges, n_metrics_fraction));
                case 'prat',    rho_estimate = @(t_map) gather(obj.max_p_relative(t_map, n_nearest));
                case 'fisher',  heuristic = 'fisher';
                case 'refmmin'
                    obj.refine_aliens(k_alien, n_nearest)
                    rho_estimate = @(t_map) gather(obj.buryi(t_map));
                otherwise,      error('����������� ����� ����� ��������. ��������������: nmin, buryi, mhist, minalien')
            end
            dimensions = dimensions(dimensions<=dim_fv); % �������� ������ ������������
            fprintf('(%s) H.=%s, ',datetime('now','Format','HH:mm:ss'), heuristic);
            switch heuristic
                case 'adddel',              [fspace_qs, fspace_map] = ...
                        obj.add_del_reduction(rho_estimate, dimensions);
                case 'gray',                [fspace_qs, fspace_map] = ...
                        obj.gray_reduction(rho_estimate, dimensions);
                case 'fisher',              [fspace_qs, fspace_map] = ...
                        obj.fisher_selection(dimensions);
                otherwise, error('������ ������ �������������� ��������� ��������')
            end
            fprintf('%s: [ %s ]\n',reduction_type,num2str(fspace_qs.', '%6.3f '));
            if nargout==3, varargout{1} = dimensions; end
        end
        
        
    end
    
    methods(Static=true)
        
        % ������ ������������� ������ �� ����������
        function [ c1 ] = get_correction_components( c1, c2 )
        %fRhoCalc ������� ������� �������� ������� ����� ��������� ���� �������
        %   ������� ��������� ��� ������ � m(C1Data) � n(C2Data) ���������.
        %   [mem_limit_MB] 
        %   ������������ �������� ����� ����� ��������
        %   �� ������ (v_dim x m x n) ������� ��������� 
            n_obj = [size(c1,1) size(c2,1)];     % ����������� ����� �������� � �������
            c1 = repmat(c1.',1,1,n_obj(2)) - repmat(permute(c2,[2 3 1]), 1, n_obj(1), 1);
            
        end
        
    end
    
end

