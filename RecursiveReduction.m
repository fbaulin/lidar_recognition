classdef RecursiveReduction < handle
    %RecursiveReduction Вычислительно эффективные методы отбора признаков
    %   Методы отбора признаков реализованные на основе метрического анализа пространств признаков 
    
    properties
        rhos_sqr        % текущие метрики
        feat_map        % текщая настройка пространства признаков
        fvs             % векторы признаков
        targets_ids     % идентификационные номера, присваиваемые внотри этого класса
        alien_mask
        weight_fun      % взвешивающая функция
        
        use_gpu         % импользовать gpu
        reduction       % метод редукции
        heurist         % эвристический алгоритм
        estimate_q      % функция оценки качества пространства
        
    end
    
    properties(Hidden=true)
        upd=true    % обновить метрики на текущем расчете
        gpu_supported_methods = { 'test', 'nmin', 'buryi', 'buryisum', 'refmmin', 'integrnmin', 'integrnmax', 'integr', 'auhist', 'integr_min', 'minalien', 'prat'}
    end
    
    methods
        function obj = RecursiveReduction(varargin)
        %RecursiveReduction Используется для задания параметров редукции
        %   Detailed explanation goes here
        %   Именованные параметры:
        %   'gpu_on' - 'True'/'False' - использовать графический адаптер
        %   'reduction' - тип редукции
        %       'buryi'
        %       'nmin'
        %       'mhist'
        %       'minalien'
        %       'prat'
        %       'wprat'
        %       'fisher'
        %       'refmin'
        %   'heuristic' - эвристический метод сокращения вычислений
        %       'adddel' - метод поочередного добавления и удаления компонент (сканируется не все
        %       пространство признаков)
        %       'gray'  - кодирование пространств кодом грея
        %   'dimensions' - размерности пространства признаков для которых нужно найти решения
        %   Гиперпараметры методов:
        %   'n_metrics_fraction'
        %   'm_technique'
        %   'k_alien'
        %   'n_nearest'
        %   'hist_edges'
            
            kw = KeywordArguments();
            obj.use_gpu = kw.get_value(varargin, 'gpu_on', gpuDeviceCount());
            obj.heurist = kw.get_value(varargin, 'heuristic', 'adddel');          % эвристика сканирования
            obj.reduction = kw.get_value(varargin, 'reduction_type', 'buryi');   % исследуемые размерности
            if or(strcmp(obj.heurist, 'fisher'), strcmp(obj.reduction, 'fisher'))
                obj.reduction = 'fisher'; obj.heurist = 'fisher';
            end
            obj.setup_reduction(varargin);
    
        end
        
        function [fspace_qs, fspace_map, varargout] = reduce(obj, fvs, targets, varargin)
            kw = KeywordArguments();
            obj.use_gpu = and(...
                any(strcmp(kw.get_value(varargin, 'reduction_type'), obj.gpu_supported_methods)), ...
                kw.get_value(varargin, 'gpu_on', gpuDeviceCount()) );    % обновить если на вход передан новый парам
            obj.heurist = kw.get_value(varargin, 'heuristic', obj.heurist);    
            obj.write_data(fvs, targets);    % запись данных
            obj.setup_reduction(varargin{:});   % настройка редукции
            
            dimensions = kw.get_value(varargin,'dimensions',1:size(fvs,2));
            dimensions = dimensions(dimensions<=size(fvs,2)); % обновить вектор размерностей
            
            fprintf('(%s) H.=%s, ',datetime('now','Format','HH:mm:ss'), obj.heurist);
            switch obj.heurist
                case 'adddel',              [fspace_qs, fspace_map] = ...
                        obj.add_del_reduction(obj.estimate_q, dimensions);
                case 'gray',                [fspace_qs, fspace_map] = ...
                        obj.gray_reduction(obj.estimate_q, dimensions);
                case 'fisher',              [fspace_qs, fspace_map] = ...
                        obj.fisher_selection(dimensions);
                otherwise, error('Ошибка выбора эвристического алгоритма редукции')
            end
            fprintf('%s: [ %s ]\n',obj.reduction,num2str(fspace_qs.', '%6.3f '));
            if nargout==3, varargout{1} = dimensions; end
            
        end
        
        function write_data(obj, fvs, targets)
            dim_fv = size(fvs,2);            % размерность ВП
            obj.feat_map = zeros(1,dim_fv,'single');
            [~, ~, obj.targets_ids] =  unique(targets); 
            if obj.use_gpu
                obj.fvs = gpuArray(fvs);
                obj.targets_ids = gpuArray(uint8(obj.targets_ids));
                obj.alien_mask = gpuArray(and(...
                    obj.targets_ids~=obj.targets_ids.',...
                    triu(ones(size(obj.targets_ids,1),'like',obj.targets_ids),1)));     % маска расчета метрик                        
            else
                obj.fvs = fvs;
                obj.targets_ids = uint8(obj.targets_ids);
                obj.alien_mask = and(...
                    obj.targets_ids~=obj.targets_ids.',...
                    triu(ones(size(obj.targets_ids,1),'logical'),1));   % маска расчета метрик
            end            
            obj.rhos_sqr = zeros(size(fvs,1),'like',fvs);
        end
        
        function setup_reduction(obj, varargin)
            kw = KeywordArguments();
            obj.reduction = kw.get_value(varargin, 'reduction_type', obj.reduction);      % исследуемые размерности
            n_metrics_fraction = kw.get_value(varargin, 'n_metrics_fraction', 0.05);      % число метрик
%             hist_edges =    kw.get_value(varargin, 'hist_edges', linspace(0,2,51));  % границы интервалов гистограммы
            n_nearest =     kw.get_value(varargin, 'n_nearest', 5);               % число проверяемых соседей
            m_technique =   kw.get_value(varargin, 'm_technique','mean');            % метод расчет
            k_alien =       kw.get_value(varargin, 'k_alien', 2);                   % число чужих
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
                
                otherwise,      error('Неправильно задан метод редукции. Поддерживаемые: nmin, buryi, mhist, minalien')
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
        
        %% Эвристические алгоритмы редукции
        
        % Редукция на основе кодирования вектора компонент кодом Грея
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
        
        % Редукция на основе сравнения метрики Фишера
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
        
        % Алгоритм поочередного удаления и добавления компонент
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
           
        
        
        %% Меры качества пространства
        
        % Расчет максимизируемого числа своих
        function space_q = minalien(obj, t_map, n_nearest, m_technique, k_alien)            
            
            n_samples = size(obj.fvs, 1);
            [~, sort_indexes] = mink(obj.get_metrics(t_map), n_nearest+1, 2);  % получить порядок соседей для каждого вектора
            %FIXME: В некоторых случаях ближайший может не быть самим вектором (когда до других векторов тоже 0)
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
        
        function space_q = local_metric_sum(obj, t_map, n_nearest)            
        %local_metric_sum выполняет оценку суммы минимальных метрик с заданными ближайшими соседями
        %NOTE: ИСПОЛЬЗОВАНИЕ НЕЭФФЕКТИВНО
            n_samples = size(obj.fvs, 1);
            rhos = obj.get_metrics(t_map);                                  % матрица метрик
            rhos(1:n_samples+1:end)=nan;
            least_rhos = zeros(n_samples, n_nearest,'like', obj.fvs);       % матрица наименьших метрик
            neq_matrix = false(n_samples, n_nearest,'like', obj.alien_mask);
            for i=1:n_nearest 
                [least_rhos(:,i), sort_indexes] = min(rhos,[], 2);                  % минимальная метрика
                rhos(sub2ind(size(rhos), (1:n_samples).', sort_indexes)) = nan;
                neq_matrix(:,i) = ...
                    obj.targets_ids == obj.targets_ids(sort_indexes); % i+1, т.к. 1й - сам вектор признаков
            end
            
            space_q = sum(least_rhos(neq_matrix),'all')-sum(least_rhos(~neq_matrix),'all');
        end
        
        function space_q = nmin_metric_balanced(obj, t_map, n_min_ratio)
        %nmin_metric выполнить оценку устойчивости пространства признаков по заданной вероятности
        %   fvs_red     - массив ячеек с векторами признаков 
        %   n_min_ratio - доля метрик, по которым считается величина метрики (<1)
        %   w_rep       - коэффициенты репликации
            rhos = obj.get_metrics(t_map);              % апдейт метрик
            rhos = sort( rhos(obj.alien_mask), 'ascend');  % получение массива метрик с использованием треугольной матрицы
            n_rho_index = floor(length(rhos)*n_min_ratio);
            space_q = rhos(n_rho_index);        % поиск n-ной метрики
        end
        
        function space_q = integr_nmin(obj, t_map, n_min_ratio)
            rhos = obj.get_metrics(t_map);
            rhos = sort(rhos(obj.alien_mask), 'ascend');
%             if obj.upd, Loggers.log('integr_16db.log',{sum(t_map),rhos}, 'sep', ','); end
            n_rho_index = floor(length(rhos)*n_min_ratio);
            space_q = sum(rhos(1:n_rho_index));        % поиск n-ной метрики
        end
        
        function space_q = integr_min(obj, t_map)
        % оценить сумму метрик между ближайшими соседями другого класса
        %NOTE: ИСПОЛЬЗОВАНИЕ НЕЭФФЕКТИВНО
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
            space_q = sum(rhos(n_rho_index:end));        % поиск n-ной метрики
        end
        
        function space_q = integr(obj, t_map)
            rhos = obj.get_metrics(t_map);
            space_q = sum(rhos(obj.alien_mask));
        end
        
        function space_q = auhist(obj, t_map)
            % Мера основанная на комбинации значения суммы метрик и меры AUC ROC
            rhos = obj.get_metrics(t_map);
            rhos = rhos(obj.alien_mask);
            sum_rhos = sum(rhos);
            max_rhos = max(rhos);
            space_q = sum_rhos*sum_rhos/(length(rhos)*max_rhos);
        end 
        
        function space_q = buryi_sum(obj, t_map)
            rhos = obj.get_metrics(t_map);      % апдейт метрик
%             n_samples = size(obj.fvs, 1);       % число объектов
%             if obj.upd, Loggers.log('.\metric_log\buryi_16db.log',{sum(t_map),gather(sort(rhos(obj.alien_mask)))}, 'sep', ',','header',''); end
            space_q = sum(rhos, 'all');
            % получение массива метрик с использованием треугольной матрицы
        end
        
        function space_q = buryi(obj, t_map)
            rhos = obj.get_metrics(t_map);      % апдейт метрик
%             n_samples = size(obj.fvs, 1);       % число объектов
%             if obj.upd, Loggers.log('.\metric_log\buryi_16db.log',{sum(t_map),gather(sort(rhos(obj.alien_mask)))}, 'sep', ',','header',''); end
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
            [~, local_groups] = mink(obj.get_metrics(t_map),n_nearest, 2);  % получить порядок соседей для каждого вектора
            local_groups = obj.targets_ids(local_groups);       % № объекта -> классы
            group_modes = mode(local_groups,2);                 % превалирующий класс в группе
            space_q = 1-mean(sum(obj.weight_fun.*(group_modes~=local_groups),2))/n_nearest; % взвешенная сумма чужих   
        end
        
        %% Прочие служебные функции
        
        % Расчет метрик
        function [ rhos ] = get_metrics( obj, t_map )
        %get_metrics Функция расчета метрики между объектами двух классов
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
        
        
        % Функция очистки областей пересечения
        function [] = refine_aliens(obj, k_alien, n_nearest)
            
            % рассчитать метрики
            % сортировать
            % удалить            
            update_status = obj.upd;
            obj.upd = false;
            [n_samples, dim_fv] = size(obj.fvs);
            t_map = ones(1,dim_fv,'single');
            [~, sort_indexes] = mink(obj.get_metrics(t_map), n_nearest+1, 2);  % получить порядок соседей для каждого вектора
            neq_matrix = zeros(n_samples, n_nearest, 'like', obj.targets_ids); % предсоздание м-цы эквив-сти
            
            for i=1:n_nearest 
                 neq_matrix(:, i) = ...
                     obj.targets_ids ~= obj.targets_ids(sort_indexes(:, i+1)); % i+1, т.к. 1й - сам вектор признаков
            end
            % получена матрица эквивалентности n_samples x n_nearest
            del_idx = sum(neq_matrix, 2)>(k_alien);
            
            obj.fvs(del_idx,:)=[];
            obj.rhos_sqr = zeros(size(obj.fvs,1),'like',obj.fvs);   % текущие метрики
            obj.targets_ids(del_idx,:)=[];
            obj.alien_mask = and(...
                obj.targets_ids~=obj.targets_ids.',...
                triu(ones(size(obj.targets_ids,1),'like',obj.targets_ids),1)); % маска расчета метрик
            obj.upd = update_status;
        end
       
        
        %% legacy
        % рудукция пространства признаков
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
        % редукция пространства признаков в режиме ЦП
        function [fspace_qs, fspace_map, varargout] = fs_reduction_cpu(obj, fvs, targets, varargin)
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
            [~, ~, obj.targets_ids] = unique(targets);        % преобразовать имена в id-номера
            obj.alien_mask = and(...
                obj.targets_ids~=obj.targets_ids.',...
                triu(ones(size(obj.targets_ids,1),'logical'),1)); % маска расчета метрик             
            kwargs = KeywordArguments(...
                'dimensions',1:dim_fv, ...          % исследуемые размерности
                'reduction_type', 'fisher',...      % исследуемые размерности
                'n_metrics_fraction', 0.05,...      % число метрик
                'hist_edges', linspace(0,2,51),...  % границы интервалов гистограммы
                'n_nearest', 5, ...                 % число проверяемых соседей
                'm_technique','mean',...            % метод расчет
                'k_alien', 2, ...                   % число чужих
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
                case 'refmmin'
                    obj.refine_aliens(k_alien, n_nearest)
                    rho_estimate = @(t_map) obj.buryi(t_map);
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
        
        % редукция пространства признаков в режиме ГП
        function [fspace_qs, fspace_map, varargout] = fs_reduction_gpu(obj, fvs, targets, varargin)
        %fs_reduction редукция с учетом размерности выборки
        % на выходе
        % fspace_rhos - вектор-столбец с метриками
        % fspace_maps - матрица по строкам которой расположены настройки пространства
        % Аргументы: 
        % fvs - вектора признаков (по строкам)
        % Именованные параметры:
        % dimensions, reduction_type, n_metrics_fraction, hist_edges,
        % obj_weight
            obj.fvs = gpuArray(fvs);
            dim_fv = size(fvs,2);            % размерность ВП
            obj.feat_map = zeros(1,dim_fv,'single');
            obj.rhos_sqr = zeros(size(fvs,1),'single','gpuArray');
            [~, ~, obj.targets_ids] =  unique(targets);        % преобразовать имена в id-номера
            obj.targets_ids = gpuArray(uint8(obj.targets_ids));
            obj.alien_mask = gpuArray(and(...
                obj.targets_ids~=obj.targets_ids.',...
                triu(ones(size(obj.targets_ids,1),'like',obj.targets_ids),1))); % маска расчета метрик
                        
            kwargs = KeywordArguments(...
                'dimensions',1:dim_fv, ...          % исследуемые размерности
                'reduction_type', 'fisher',...      % исследуемые размерности
                'n_metrics_fraction', 0.05,...      % число метрик
                'hist_edges', linspace(0,2,51),...  % границы интервалов гистограммы
                'n_nearest', 5, ...                 % число проверяемых соседей
                'm_technique','mean',...            % метод расчет
                'k_alien', 2, ...                   % число чужих
                'heuristic', 'adddel');               % эвристика сканирования
            [dimensions, reduction_type, n_metrics_fraction, hist_edges, ...
                n_nearest, m_technique, k_alien, heuristic] =  ...
                kwargs.parse_input_cell(varargin);
            % настройка специфическая для метода
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

