classdef QualityCheck < handle
    %QualityCheck Класс для оценки качества распознавания
    %   Выполняет оценку качества распознавания
    
    properties
        reduction_method    % метод редукции пространства признаков (метрики, гистограммы, etc)
        heuristic           % метод сокращения числа вычислений
        classifier_type     % название классификатора 
        dimensions          % размерность пространства признаков
        clf                 % классификатор, поддерживающий configure и train
        decision_mode       % принцип работы решающего устройства
        use_GPU = 'no'      % 
        fs_map              % карта признаков
        score_type          % тип меры качества распознавания
        
        maps_filename       % имя файла с картами пространств признаков
        
        batch_filename      % имя файла для сохранения 
        feature_header      % заголовок для файла
        
        n_nearest           % число ближайших соседей
        m_technique         % тип метрики ближйших соседей
        k_alien             % граничное число чужих
        n_metrics_fraction  % доля метрик для которых допустимо быть меньше заданного значения
        hist_edges          % границы гистограмм для метода оценки по гистограммам
    end

   
    methods
        
        % Конструктор класса
        function obj = QualityCheck(varargin)
        %QualityCheck Сформировать реализацию класса
        %   Сформировать реализацию класса
        %   Args
        %       'reduction_method'  - метод редукции [nmin, buryi, hist].
        %       'classifier_type'   - тип классификатора, 'mlp' - многослойный персептрон.
        %       'layers',           - задать структуру скрытых слоев, [160 80 40].
        %       'dimensions'        - размерности, 1.
        %       'batch_filename'    - имя файла обучающей выборки, false.
        %       'maps_filename'     - имя файла с картами признаков,  false.
            kwargs = KeywordArguments(...
                'reduction_method','nmin',...
                'heuristic', 'adddel',...   % метод перебора
                'classifier_type', 'none',...
                'dimensions', 1, ...
                'batch_filename', false,...
                'maps_filename', false,... 
                'n_nearest',5, ...          % число проверяемых соседей
                'm_technique','mean',...    % метод расчета
                'k_alien',2, ...               % число своих
                'n_metrics_fraction', 0.05,...
                'hist_edges', linspace(0,2,51)...  
                );   % значения по-умолчанию
            [   obj.reduction_method, ...
                obj.heuristic,...
                obj.classifier_type, ...
                obj.dimensions,...
                obj.batch_filename,...
                obj.maps_filename,...
                obj.n_nearest,...
                obj.m_technique,...
                obj.k_alien,...
                obj.n_metrics_fraction,... 
                obj.hist_edges] = kwargs.parse_input_cell(varargin); 
            if strcmp(obj.classifier_type,'mlp')
                obj.setup_classifier(...
                    'classifier_type',obj.classifier_type,...
                    'layers',[40 20 10])
            end
        end
                
        % Настройка классификатора
        function setup_classifier(obj, varargin)
        %SETUP_CLASSIFIER Настройка классификатора
        %   Args, default:
        %       'classifier_type', 'mlp'    - тип классификатора    
        %       'layers', [40 20 10]        - число слоёв
        %       'performance', 'mse'        - функция по которой расчитывается качество работы
        %       'train_fun', 'default'      - метод обучения
        %       'gpu', obj.use_GPU          - использовать ГПУ для обучения
            kwargs = KeywordArguments(...
                'classifier_type', 'mlp',...    
                'layers', [40 20 10],...
                'performance', 'default',...
                'decision_mode', 'metric',...   % если 'metric', то решение принимается по метрике, иначе - по 
                'train_fun', 'default',...
                'gpu', obj.use_GPU);
            [   obj.classifier_type, ...
                layers, ...
                performance,...
                obj.decision_mode, ...
                train_fun,...
                obj.use_GPU     ] = kwargs.parse_input_cell(varargin);
                
            switch obj.classifier_type
                case 'mlp'   
                    if strcmp(train_fun, 'default') % если значение по-умолчанию, то
                        if strcmp(obj.use_GPU, 'yes')
                            train_fun = 'trainscg'; % для ГПУ - градиентная
                        else
                            train_fun = 'trainlm';  % для ЦПУ - Левенберг-Марк.
                        end
                    end
                    obj.clf = feedforwardnet(layers,train_fun); % создать персептрон
                    if strcmp(performance,'default'), performance='mse'; end
                    obj.clf.performFcn = performance;           % функция ошибок
                    fprintf('\tНейросеть: [%s], качество: %s\n',...
                        num2str(layers,'%d '), performance);
                case 'pnet'
                    if strcmp(performance,'default'), performance='crossentropy'; end
                    if strcmp(train_fun, 'default'), train_fun = 'trainscg'; end
                    obj.clf = patternnet(layers,train_fun,performance);       % создать персептрон
                    fprintf('\tНейросетевой классификатор: [%s], качество: %s\n',...
                        num2str(layers,'%d '), performance);
                otherwise               % если заданный классификатор не реализован
                    error('Тип классификатора не поддерживается')
            end
            obj.clf.trainParam.showWindow = false;      % отключить окно со стат. обучения
            obj.clf.trainParam.showCommandLine = false; % отключить вывод в консоль
            
        end
        
        % Подготовка обучающих выборок
        function [ x_train, y_train, x_test, y_test ] = prepare_dataset(obj, x_train, y_train_str, x_test, y_test_str)
            %NOTE: при выполненнии отбора признаков после нормировки данных
            % ошибки на распознавании оказываются выше, чем ошибки в случае
            % выполнения нормировки после отбора данных. Это
            % свидетельствует о том, что при применении отбора признаков по
            % алгоритму nmin важно, чтобы СКО аддитивного шума было равным
            % для всех признаков (не путать с СКО самих признаков)
            n_test = size(y_test_str,1);
            n_train = size(y_train_str,1);
            y = vertcat(y_train_str, y_test_str);       % создать общую выборку для правильного формирования целевого вектора
            y = obj.string2vector(y);                   % преобразовать названия классов в векторы
            y_train = y(1:size(y_train_str, 1),:);      % вытащить из полученного вектора обучающую выборку
            y_test =  y(size(y_train_str, 1)+1:end,:);      % вытащить тестовую выборку
            
            data_means = mean(x_train,1);       % среднее значение
            data_stds = std(x_train,[],1);      % расчет СКО
            
            x_train =  (x_train - data_means)...
                ./repmat(data_stds,n_train,1);         % нормировать тестовые вектора признаков;
            x_test = (x_test - data_means)...
                ./repmat(data_stds,n_test,1);          % нормировать тестовые вектора признаков
            
        end
        
        % Оценить качество распознавания
        function [score, varargout] = assess_quality(obj, x_train, y_train_str, x_test, y_test_str)
        %assess_quality Оценить качество выделения признаков и распознавания
        %   Выполнить выделение признаков и оценку качества распознавания для размерностей
        %   заданных в реализации класса.
        %   Args:
        %       x_train     - вектора признаков обучающей выборки.
        %       y_train     - решения (названия объектов) обучающей выборки.
        %       x_test      - вектора признаков тестовой выборки.
        %       y_test      - решения (названия объектов) тестовой выборки.
        %   Returns: 
        %       score       - мера качества распознавания, обеспечиваемого при классификации.
        %       [rho_min]   - мера качества сформированного пространства признаков. 
        %       [fs_maps]   - карта соответствующего пространства признаков <n_dim x fv_dim>.
        %       [conf_mx]   - матрица ошибок для оптимальной размерности.

%             fvs = obj.wrap_data(x_train,y_train_str);   % преобразовать данные в массив ячеек
%             [rho_min, fs_maps,obj.dimensions] = SystemModel.fs_reduction(fvs,...
%                 'dimensions',obj.dimensions, 'reduction_type',obj.reduction_method...
%                 );   % подобрать оптимальное пространство признаков для заданной размерности

            reduction = RecursiveReduction();
            [rho_min, fs_maps,obj.dimensions] = reduction.reduce(x_train, y_train_str,...
                'dimensions',obj.dimensions, 'reduction_type',obj.reduction_method,...
                'heuristic',obj.heuristic,...
                'n_metrics_fraction',obj.n_metrics_fraction,...
                'hist_edges', obj.hist_edges,...
                'n_nearest', obj.n_nearest,...
                'm_technique', obj.m_technique,...
                'k_alien',obj.k_alien);   % подобрать оптимальное пространство признаков для заданной размерности
            
            obj.fs_map = logical(fs_maps);      % преобразовать в логические для индексации
            n_dims = length(rho_min);           % обновить число исследованных размерностей
            score = zeros(n_dims,1);            % вектор-столбец с оценкой рез-тата распознавания
            if nargout>=4, conf_mx = cell(n_dims,1); end        % матрица ошибок
%             conf_mx = cell(length(fvs), length(fvs));
            % подготовка данных: в т.ч. нормализация
            [ x_train, y_train, x_test, y_test ] = obj.prepare_dataset(x_train, y_train_str, x_test, y_test_str);
            
            if strcmp(obj.use_GPU,'yes') % обучение на ГП не работает с single признаками
                x_train = double(x_train);
                x_test = double(x_test);
            end
            n_arg_out = nargout;
            parfor i_dim = 1:n_dims
                x_train_red = x_train(:, obj.fs_map(i_dim,:));   %#ok<PFBNS> % редуцировать обучающую выборку
                x_test_red = x_test(:, obj.fs_map(i_dim,:));     %#ok<PFBNS> % редуцировать тестовую выборку
               % обучить классификатор
                inst_clf = configure(obj.clf,x_train_red.',y_train.');   % конф. входного и выходного слоёв
                inst_clf = train(inst_clf,x_train_red.',y_train.','useGPU',obj.use_GPU);       % обучение
                y_pred = inst_clf(x_test_red.').';                       % формирование ответов для тестовой выборки
                score(i_dim) = perform(inst_clf, y_test.', y_pred.');    % расчет качетсва распознавания
                if n_arg_out>=4
                    switch obj.decision_mode
                        case 'metric', y_pred = QualityCheck.metric_decision(y_pred);
                        case 'max'
                        otherwise, error(['Задан неизвестный метод принятия решения (', obj.decision_mode, ')'])
                    end
                    [~,i_tst] = max(y_test,[],2); [~,i_pred] = max(y_pred,[],2);    % one-hot -> индексы классов
                    conf_mx{i_dim} = confusionmat(i_tst,i_pred);
                end    % расчитать матрицу
            end
            if nargout>=2, varargout{1} = rho_min; end      % если более одного арг - выдать rho_min
            if nargout>=3, varargout{2} = fs_maps; end      % если 3 аргумента выдать карту признаков
            if nargout>=4, varargout{3} = permute(cat(3,conf_mx{:}),[3 1 2]); end      % выдать матрицу ошибок
        end
        
        
        
    end
    
    methods(Access = private, Hidden = true)
        
        % Сохранить выборку в csv
        function [] = save_batch(obj, features, meta)
        %SAVE_CSV Сохраняет в csv признаки
        %   В первой строке содрежатся названия признаков, а t обозначает
        %   объект. Так, для признаков. 
        %   Args:
        %       features    - матрица признаков.
        %       meta        - строковые обозначения элементов (имя объекта [ракурс]).
        %       f_header    - признаков.
        %       filename    - имя файла.
            f_header = obj.feature_header;
            filename = obj.batch_filename;
            [ n_obj, ~ ] = size(features);
            %meta = RPTools.meta2str(meta); % восстановить если понадобится вывод с
            %ракурсами
            %meta = {meta.name}; % 
            fid = fopen(filename,'w');
            fprintf(fid,'%s,',f_header);
            fprintf(fid,'%s\n','t');
            for i_line = 1:n_obj
                fprintf(fid,'%f,', features(i_line,:));
                fprintf(fid,'%s\n', meta(i_line).name);
            end
            fclose(fid);
            
            disp(['Файл ' filename ' сохранен'])
            
        end
        
    end
    
    methods(Static=true)
        
        function y_dec = metric_decision(y_pred, varargin)
        %   
        %   Если расстояние до целевого вектора (one hot векторы) меньше порогового (по умолчанию (sqrt(2)/2), то
        %   вектор заменяется на соответствующий one-hot вектор, иначе - на нулевой вектор.
            if nargin==2, thresh_sqr = sqrt(varargin{1}); else thresh_sqr = 1/2; end
            [n_samples, n_class] = size(y_pred);    % определить число классов
            y_dec = zeros(n_samples, n_class, 'like', y_pred);
            % цикл по классам 
            for i=1:n_class
                y_h = zeros(n_samples, n_class, 'like', y_pred);
                y_h(:,i) = 1;
                dists = sum((y_h-y_pred).^2, 2);    % рассчитать квадраты метрик до целевого вектора класса
                y_dec(dists<thresh_sqr, i) = 1;            % забить все векторы для которых расстояние меньше порога
            end
            y_dec(sum(y_dec,2)>1,:) = 0; % выбить векторы, для которых более одного рещения 
        end
        
        % Преобразовать выборку из набора ячеек в единый массив
        function [x,y] = unwrap_cell_data(fvs,meta,varargin)
        %unwrap_cell_data Преобразовать массив ячеек в данные
        %   Переводит данные, разбитые по ячейкам в матрицу признаков и массив ячеек с именем
        %   класса для каждого элемента выборки.
        %   Args:
        %       fvs     - массив ячеек с векторами признаков.
        %       meta    - массив ячеек со структурой с именами классов.
        %       'matrix_output' - вместо массива ячеек матрица,  false. 
            kwargs = KeywordArguments('matrix_output', false );
            [ matrix_output ] = kwargs.parse_input_cell(varargin);
            n_obj = length(fvs);    % число объектов
            x = vertcat(fvs{:});    % свернуть все вектора признаков в одну матрицу
            y = arrayfun(@(i) {meta{i}.name}.', 1:n_obj, 'UniformOutput', false); % вывести названия классов
            y = vertcat(y{:});  % свернуть массивы (по объектам) массивов ячеек, в один массив ячеек
            if matrix_output,   y = cellfun(@(name) string(name), y);    end     % свернуть массив ячеек в матрицу
        end
        
        % Преобразовать данные из массива в массив ячеек по объектам
        function [varargout] = wrap_data(x,y)
        %wrap_data преобразовать массив с разметкой в ячейки.
        %   Переводит данные x в массивы ячеек, соответствующие объектом
        %   информация о разбивке берется из y.
        %   Args:
        %       x - марица с векторами признаков.
        %       y - массив ячеек с именами объектов.
        %   Returns:
        %       массив ячеек, содержащий матрицы векторов признаков, число ячеек в котором равно числу объектов.
        %       если возвращается 2 результата, то второй - имена объектов.
            [object_names, ~, ci] = unique(y);  % собрать имена объектов
            % здесь ci - номера объекта, соответствует i-му имени в массиве object_names
            n_objects = length(object_names);   % число объектов
            fvs = cell(n_objects, 1);           % сформировать массив с ячейкой под каждый объект
            for i_obj = 1:n_objects             % цикл по объектам
                fvs{i_obj} = x(i_obj==ci,:);    % найти вектора признаков i-го объекта и упаковать их
            end
            varargout{1} = fvs;
            if nargout>1, varargout{2} = object_names; end
                
        end
        
        % Сформировать целевые one-hot векторы
        function [y_vect] = string2vector(y_string)
        %string2vector  Сформировать целевые векторы на основании массива
            [object_names, ~, ci] = unique(y_string);           % выделить имена объектов
            y_vect = zeros(length(ci), length(object_names));   % сформировать матрицу под целевые вектора
            for i_obj = 1:length(object_names)  % цикл по объектам
                y_vect(i_obj==ci,i_obj) = 1;    % в столбцах нужных объектов поставить единицу
            end 
        end
        
    end
    
end

