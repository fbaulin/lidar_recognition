classdef RecursiveReduction < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        rhos_sqr    % текущие метрики
        feat_map    % текщая настройка пространства признаков
        fvs         % векторы признаков
        targets     % разметка векторов признаков
        targets_ids  % 
        weight_fun   % взвешивающая функция
    end
    
    properties(Hidden=true)
        upd = true  % обновить метрики на текущем расчете
        alien_mask
    end
    
    methods
        function obj = RecursiveReduction()
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            
        end
        
        % редукция пространства признаков
        function [fspace_qs, fspace_map, varargout] = fs_reduction(obj, fvs, targets, varargin)
        %fs_reduction редукция с учетом размерности выборки
        % на выходе
        % fspace_rhos - вектор-столбец с метриками
        % fspace_maps - матрица по строкам которой расположены настройки пространства
        % Аргументы: 
        % fvs - вектора признаков (по строкам)
        % Именованные параметры:
        % dimensions, reduction_type, n_metrics_fraction, hist_edges,
        % obj_weight
            obj.fvs = fvs;
            dim_fv = size(fvs,2);            % размерность ВП
            obj.feat_map = zeros(1,dim_fv,'single');
            obj.rhos_sqr = zeros(size(fvs,1),'like',fvs);
            obj.targets = targets;
            [~, ~, obj.targets_ids] = unique(targets);        % преобразовать имена в id-номера
            obj.alien_mask = and(...
                obj.targets_ids~=obj.targets_ids.',...
                triu(ones(size(obj.targets_ids,1),'logical'),1)); % маска расчета метрик
            
            kwargs = KeywordArguments(...
                'dimensions',1:dim_fv, ...          % исследуемые размерности
                'reduction_type', 'none',...        % исследуемые размерности
                'n_metrics_fraction', 0.05,...      % число метрик
                'hist_edges', linspace(0,2,51),...  % границы интервалов гистограммы
                'n_nearest', 5, ...                 % число проверяемых соседей
                'm_technique','mean',...            % метод расчет
                'k_alien', 2, ...                   % число своих
                'heuristic', 'adddel');               % эвристика сканирования
            [dimensions, reduction_type, n_metrics_fraction, hist_edges, ...
                n_nearest, m_technique, k_alien, heuristic] =  ...
                kwargs.parse_input_cell(varargin);
            % настройка специфическая для метода
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
                otherwise,      error('Неправильно задан метод редукции. Поддерживаемые: nmin, buryi, mhist, minalien')
            end
            dimensions = dimensions(dimensions<=dim_fv); % обновить вектор размерностей
            fprintf('(%s) H.=%s, ',datetime('now','Format','HH:mm:ss'), heuristic);
            switch heuristic
                case 'adddel',              [fspace_qs, fspace_map] = ...
                        obj.add_del_reduction(rho_estimate, dimensions);
                case 'gray',                [fspace_qs, fspace_map] = ...
                        obj.gray_reduction(rho_estimate, dimensions);
                case 'fisher',              [fspace_qs, fspace_map] = ...
                        obj.fisher_selection(dimensions);
                otherwise, error('Ошибка выбора эвристического алгоритма редукции')
            end
            fprintf('%s: [ %s ]\n',reduction_type,num2str(fspace_qs.', '%6.3f '));
            if nargout==3, varargout{1} = dimensions; end
        end
        
        
        function [fspace_qs, fspace_map] = gray_reduction(obj, rho_estimate, dimensions)
            dim_fv = size(obj.fvs,2);
            n_spaces = length(dimensions);
            fspace_qs = zeros(n_spaces,1);          % значения пороговой метрики для размерностей
            fspace_map = zeros(n_spaces, dim_fv);   % матрица оптимальных пространств - по строкам от размерности
            n_spaces = 2^dim_fv-1;
            
            line_len = fprintf('Пространств: %d. \tВыполнено: 0.00%%', n_spaces);
            cur_percent = 0;
            obj.upd = true;
            for feat_val = 1:n_spaces               % по всем возможны пространствам признаков
                t_feat_map = de2bi(feat_val,dim_fv,'left-msb'); % получение карты текущего пространства
                t_feat_map_sh = circshift(t_feat_map,1);
                t_feat_map_sh(1) = 0;
                t_feat_map = xor(t_feat_map,t_feat_map_sh);     % использвоание кода грея для оптимизации рекурсии
                i_fmap = find(sum(t_feat_map)==dimensions); % индекс соответствующий размерности
                if i_fmap
                    [cur_q_val] = rho_estimate(t_feat_map);
                    if cur_q_val>=fspace_qs(i_fmap)     % если значение метрики больще значения для пр. признаков такой размерности
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
        
        function [fspace_qs, fspace_map] = fisher_selection(obj, dimensions)
        % сформировать набор классов и векторов признаков
            dim_fv = size(obj.fvs,2);       % размерность исх. выборки
            n_spaces = length(dimensions);  % число подпространств
            fspace_qs = zeros(n_spaces,1);  % значения пороговой метрики для размерностей
            fspace_map = zeros(n_spaces, length(obj.feat_map));   % матрица оптимальных пространств - по строкам от размерности
            
            m_exp_i     = mean(obj.fvs,1);          % матожидание всей выборки
            classes_ids = unique(obj.targets_ids);    % идентификаторы классов
            n_classes   = length(classes_ids);                % число классов
            m_exp_ij    = zeros(dim_fv,n_classes);  % матрица матожиданий классов
            var_ij      = zeros(dim_fv,n_classes);  % матрица СКО классов
            fisher_score = zeros(dim_fv,1);         % вектор метрик фишера
            n_j         = zeros(1, n_classes);      % число объектов в классах
            % расчет метрик фишера
            for i_class=1:n_classes
                n_j(i_class) = sum(obj.targets_ids==classes_ids(i_class));     % число векторов класса
                m_exp_ij(:,i_class) = mean(obj.fvs(obj.targets_ids==classes_ids(i_class),:),1);   % матожидание классов
                var_ij(:,i_class) = var(obj.fvs(obj.targets_ids==classes_ids(i_class),:),0,1);    % дисперсия с нормирокой на N-1
            end
            for i_feat = 1:dim_fv             % цикл по признакам
                fisher_score(i_feat) = sum(n_j.*(m_exp_ij(i_feat,:)-m_exp_i(i_feat)).^2)/sum(n_j.*var_ij(i_feat));
            end
            % расчитать критерии
            [fisher_score, fisher_ids] = sort(fisher_score,'descend');
            for i_dim = 1:n_spaces
                curr_dmty = dimensions(i_dim); % текущая размерность
                fspace_qs(i_dim) = sum(fisher_score(1:curr_dmty));
                fspace_map(i_dim, fisher_ids(1:curr_dmty))=1;    % задействовать признаки по мере                
            end
        end
        
        function [fspace_qs, fspace_map] = add_del_reduction(obj, rho_estimate, dimensions)
            
            n_spaces = length(dimensions);
            fspace_qs = zeros(n_spaces,1);          % значения пороговой метрики для размерностей
            fspace_map = zeros(n_spaces, length(obj.feat_map));   % матрица оптимальных пространств - по строкам от размерности
            
            obj.upd = false;
            cur_feat_map = obj.feat_map;
            alg_steps = [1 1 1 -1 -1];  % код алгоритма 1 - шаг вперед, -1 - шаг назад
            
            line_len = fprintf('Анализ(%02d): 00',max(dimensions));
            for i_dim = dimensions(dimensions>0) %TODO: если не последовательные размерности - будет ошибка
                % Сделать 3 шага вперед и два - назад
                fprintf('\b\b%02d',i_dim);
                for step=(alg_steps==-1)
                    % перебор всех пространств i_dim+1
                    z_idxs = find(cur_feat_map==step); % вытащить вектор с индексами нулевых элементов (m штук)
                    if or(...
                            isempty(z_idxs), ...    % если прибавлять некуда или
                            and(sum(cur_feat_map)<(i_dim+1),step)...  % целевое простраснтво уже пройдено
                            )
                        continue;  
                    end
                    q_sep_vect = zeros(length(z_idxs),1);
                    for i_idx = 1:length(z_idxs)            % перебрать m вариантов
                        t_feat_map = cur_feat_map;          % проверяемая карта
                        t_feat_map(z_idxs(i_idx)) = ~step;
                        q_sep_vect(i_idx) = rho_estimate(t_feat_map);   %
                    end % повторить
                    [~, best_idx] = max(q_sep_vect);         % найти id элемента оптимума
                    cur_feat_map(z_idxs(best_idx)) = ~step; % заменить текщее пространство
                end
                obj.upd=true; 
                dim_idx = find(dimensions==i_dim);
                fspace_qs(dim_idx) = rho_estimate(cur_feat_map); 
                fspace_map(dim_idx,:) = cur_feat_map;
                obj.upd=false; % обновить пространство
            end
            fprintf(repmat('\b', 1, line_len));
            
        end
            
        
        % Расчет максимизируемого числа своих
        function space_q = minalien(obj, t_map, n_nearest, m_technique, k_alien)            
            
            n_samples = size(obj.fvs, 1);
            [~, sort_indexes] = mink(obj.get_metrics(t_map), n_nearest+1, 2);  % получить порядок соседей для каждого вектора
            neq_matrix = zeros(n_samples, n_nearest); % предсоздание м-цы эквив-сти
            
            for i=1:n_nearest 
                 neq_matrix(:, i) = ...
                     obj.targets_ids ~= obj.targets_ids(sort_indexes(:, i+1)); % i+1, т.к. 1й - сам вектор признаков
            end
            switch m_technique
                case 'mean'     % среднее число точек другого класса среди соседей
                    q_sep = mean(sum(neq_matrix, 2)/n_nearest);
                case 'thresh'   % доля объектов для которых число ближайших соседей больше порога
                    q_sep = sum(sum(neq_matrix, 2)>(k_alien))/n_samples;
            end
            space_q = 1 - q_sep;
        end
        
        function space_q = nmin_metric_balanced(obj, t_map, n_min_ratio)
        %nmin_metric выполнить оценку устойчивости пространства признаков по заданной вероятности
        %   fvs_red     - массив ячеек с векторами признаков 
        %   n_min_ratio - доля метрик, по которым считается величина метрики (<1)
        %   w_rep       - коэффициенты репликации
            rhos = obj.get_metrics(t_map);              % апдейт метрик
%             n_samples = size(obj.fvs, 1);               % число объектов
%             alien_mask = and(...
%                 obj.targets_ids~=obj.targets_ids.',...
%                 triu(ones(n_samples,'logical'),1));
            rhos = sort( rhos(obj.alien_mask), 'ascend');  % получение массива метрик с использованием треугольной матрицы
%             n_rho_index = floor(numel(alien_mask)*n_min_ratio); % ПЕРЕСЧЕТ ИЗ ДОЛИ ОБЩЕГО ЧИСЛА МЕТРИК В ДОЛЮ ЧУЖИХ
            n_rho_index = floor(length(rhos)*n_min_ratio);
%             space_q = mean(rhos(1:n_rho_index));        % поиск n-ной метрики
            space_q = rhos(n_rho_index);        % поиск n-ной метрики
        end
        
        function space_q = buryi(obj, t_map)
            rhos = obj.get_metrics(t_map);      % апдейт метрик
%             n_samples = size(obj.fvs, 1);       % число объектов
            space_q = min(rhos(obj.alien_mask));
            % получение массива метрик с использованием треугольной матрицы
        end
        
        % Критерий пиковой плотности вероятности
        function space_q = max_p_relative(obj, t_map, n_nearest)
            [~, local_groups] = mink(obj.get_metrics(t_map),n_nearest, 2);  % получить порядок соседей для каждого вектора
            local_groups = obj.targets_ids(local_groups);   % № объекта -> классы
            [~, occur] = mode(local_groups,2);              % превалирующий класс в группе
            space_q = mean(occur)/n_nearest;                % 
%             group_modes = mode(local_groups,2);                 % превалирующий класс в группе
%             space_q = mean(sum(group_modes==local_groups,2))/n_nearest; % среднее число привал.   
        end
        
        % Критерий пиковой плотности вероятности
        function space_q = w_p_relative(obj, t_map, n_nearest)
            [~, local_groups] = mink(obj.get_metrics(t_map),n_nearest+1, 2);  % получить порядок соседей для каждого вектора
            local_groups = obj.targets_ids(local_groups);       % № объекта -> классы
            group_modes = mode(local_groups,2);                 % превалирующий класс в группе
            space_q = 1-mean(sum(obj.weight_fun.*(group_modes~=local_groups),2))/n_nearest; % взвешенная сумма чужих   
        end
        
        % Расчет метрик
        function [ rhos ] = get_metrics( obj, t_map )
        %fRhoCalc Функция расчета метрики между объектами двух классов
        %   Функция принимает два масива с m(C1Data) и n(C2Data) векторами.
        %   Рассчитывает разности между всеми строками
        %   На выходе mxn матрица евклидовых метрик 
            cor_map = t_map - obj.feat_map; % вектор отличия: убрать-1, добав+1
            if and(nnz(cor_map) >= nnz(t_map), obj.upd)   % если число коррекций не меньше числа ненулевых то расчитать заново
                obj.rhos_sqr(:) = 0;        % сбросить метрики
                obj.feat_map(:) = 0;        % сбросить карту
                cor_map = t_map;            % карта корректировки = карте целевого пространства признаков
            end
            fvs_red = obj.fvs(:,cor_map~=0);    % отобрать компоненты ВП, используемые для коррекции
            metr_dif = RecursiveReduction.get_correction_components(fvs_red,fvs_red).^2; % расчитать элементы корректировки
            
            cor_map(cor_map==0) = [];       % убрать неизменные сегменты
            
            metr_dif(cor_map.'<0,:,:) = -metr_dif(cor_map.'<0,:,:);
            metr_dif = squeeze(sum(metr_dif, 1));
            
%             metr_dif = squeeze(sum(...      % добавка вычисляется как линейная комбинация
%                 metr_dif.*...               % добавок с множителями удал/добавл
%                 repmat( reshape(cor_map,[],1),...   % репликации вектора множителей
%                 [1,size(metr_dif,2),size(metr_dif,3)] ), ...
%                 1));
            
            if obj.upd  % если стоит признак обновления, то обновить
                obj.rhos_sqr = obj.rhos_sqr + metr_dif;
                obj.feat_map = t_map;
                rhos = obj.rhos_sqr;
%                 rhos = sqrt(obj.rhos_sqr);
            else
%                 rhos = sqrt(obj.rhos_sqr + metr_dif);
                rhos = obj.rhos_sqr + metr_dif;
            end
        end
        
       
    end
    
    methods(Static=true)
        
        % Расчет распределения метрик по измерениям
        function [ c1 ] = get_correction_components( c1, c2 )
        %fRhoCalc Функция расчета квадрата метрики между объектами двух классов
        %   Функция принимает два масива с m(C1Data) и n(C2Data) векторами.
        %   [mem_limit_MB] 
        %   Рассчитывает разности между всеми строками
        %   На выходе (v_dim x m x n) матрица разностей 
            n_obj = [size(c1,1) size(c2,1)];     % Определение числа объектов в классах
            c1 = repmat(c1.',1,1,n_obj(2)) - repmat(permute(c2,[2 3 1]), 1, n_obj(1), 1);
            
        end
        
    end
end

