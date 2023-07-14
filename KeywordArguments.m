classdef KeywordArguments < handle 
    %KEYWORDARGUMENTS V0.0 Обработка входный именованных аргументов
    %TODO: Добавить метод extract_value   
    %   Класс предназначен для обработки входных именованных параметров.
    %   При разработке функции, в которой используются пары именованных параметров, значения параметров
    %   по умолчанию могут быть заданы заранее. Использования в наиболее простом варианте предусматривает
    %   два этапа: 
    %       1) создание объекта с заданием параметров и их значений по-умолчанию 
    %           порядок следования параметров не имеет значения - главное, чтобы первым шло название параметра,
    %           а за ним - его значение.
    %       kwargs = KeywordArguments(...                               % в качестве занчений можно использовать
    %           'Param_name', 'some default value',...                  % строки
    %           'Other_name', 2, ...                                    % числа, векторы и матрицы
    %           'cell parameter', {'a', 2, some_value_calced earlier},  % ячейки или значения заданные выше
    %           'fun', @(a) a^2);                                        % лямбда-функции
    %       2) обработка входных переменных varargin с использованием сгенерированного объекта:
    %       [param_internal_name, other_par, cell_par, fun] = kwargs.parse(varargin);
    %       
    %       Если функция не спользует часть параметров, а передает их в одну из вызываемых ей функций, то следует 
    %       воспользоваться функцией сохранения остатка. Эта функция позволяет сохранить параметры, не заданные
    %       при инициализации объекта. Для использования этой функции следует включить параметр keep_residue.
    %       После включения параметра keep_residue "невостребованные" параметры из varargin будут сохранены в объекте.
    %       Доступ к ним можно получить, обратившись к полю residue или с использованием функции get_residue.
    %       Этот параметр включается автоматически при вызове функции parse_with_residue(). Параметры, полученные 
    %       в результате, представляют собой массив ячеек и могут быть переданы в дочернюю функцию, напр.:    
    %       kwargs = KeywordArguments(__);
    %       kwargs.keep_residue = true;
    %       [__] = kwargs.parse(varargin);
    %       some_subfun(kwargs.residue{:});
    %                                   
    %%
    properties
        vars
        check_type = false      %TODO: Добавить элемент с описанием для каждого аргумента - повысит гибкость
        types                   % типы переменных
        option_names            % имена параметров в порядке заданном при инициализации объекта
        keep_residue = false    % сохранять элементы, не заданные при инициализации объекта
                                % если false, то при передаче непредусмотренных параметров выдается ошибка
        residue                 % здесь хранятся переданные в varargin не предусмотренные параметры
    end
    
    %%
    methods
        
        % Конструктор, инициализирующий объект
        function obj = KeywordArguments(varargin)
        %KEYWORDARGUMENTS Construct an instance of this class
        %   Передается набор параметров и значений по умолчанию
        %   
            if mod(length(varargin),2)~=0
                error('Число имен должно быть равно числу значений')
            end
            obj.vars = struct(varargin{:});         % создать структуру
            obj.option_names = varargin(1:2:end);
            
        end
               
        % Обработать varargin,
        function varargout = parse(obj, arg_in)
        %PARSE_INPUT_CELL Обработать входные параметры и значения в
        %порядке следования по-умолчанию
            obj.residue = {};
            n_options = length(obj.option_names);   % общее число опций
            if nargout ~= n_options
                error('Число входных и выходных параметров должно совпадать')
            end
            obj.parse_input_struct(arg_in);         % записать аргументы в объект
            varargout = cell(1,n_options);          % сформировать массив
            for i = 1:n_options                     
                varargout{i} = obj.vars.(obj.option_names{i}); % 
            end
        end
        
        % Обработать varargin, вернув непредусмотренные параметры
        function varargout = parse_with_residue(obj, arg_in)
        %PARSE_INPUT_CELL Обработать входные параметры и значения в
        %порядке следования по-умолчанию
        %   Функция принимает на вход массив ячеек, содержащих имена параметров и значения параметров.
        %   Все параметры, не входящие в заранее заданный список сохраняются в переменную объекта residue
        %   и возвращаются этой функцией последним элементом, являющимся массивом ячеек.
            obj.keep_residue = 1;
            obj.residue = {};
            n_options = length(obj.option_names);   % общее число опций
            if nargout ~= n_options
                error('Число входных и выходных параметров должно совпадать')
            end
            obj.parse_input_struct(arg_in);         % записать аргументы в объект
            varargout = cell(1,n_options+1);          % сформировать массив
            for i = 1:n_options                     
                varargout{i} = obj.vars.(obj.option_names{i}); % 
            end
            varargout{end} = obj.residue;
        end
        
        % Обработать varargin, вернув непредусмотренные параметры
        function varargout = get_residue(obj)
        %get_residue Возвратить параметры из varargin, не входящие в список заранее заданных.
            varargout = obj.residue;
        end
        
        % Повторить пары имя-значение, но в порядке следования, установленном при инициализации объекта
        function varargout = produce_input(obj)
        %reproduce_input возвратить набор параметров в виде массива ячеек
        %   Возвращает набор параметров в виде массива ячеек, содержащего на нечетных местах
        %   названия параметров, а на четных - их значения.
            varargout = cell(1,2*length(obj.option_names));     % создать массив ячеек
            varargout(1:2:end) = obj.option_names(:);           % заполнить названия параметров
            varargout(2:2:end) = cellfun(@(fn) obj.vars.(fn),...
                obj.option_names,'UniformOutput',false);        % заполнить значения параметров
        end
                
        % Задать проверку 
        function set_types(obj,varargin)
            if length(varargin) == 1
                if strcmp(varargin{1},'default')
                    for name = fieldnames(obj.vars)
                        obj.types.(name) = class(obj.vars.(name));
                    end
                    obj.check_type = true;  %TODO: сделать параметр векторным
                else
                    warning('Неподдерживаемый формат ввода')
                end
            elseif mod(length(varargin),2)==0 
                for pair = reshape(varargin,2,[])
                    obj.types.(pair{1}) = pair{2};
                end
            else
                warning('Нечетное число элементов. В функцию set_types следует передавать пары.')
            end
            
        end        
        
        %% legacy v0.0
        function varargout = parse_input_cell(obj, arg_in)
        %PARSE_INPUT_CELL Обработать входные параметры и значения в
        %порядке следования по-умолчанию
            obj.residue = {};
            n_options = length(obj.option_names);   % общее число опций
            if nargout ~= n_options
                error('Число входных и выходных параметров должно совпадать')
            end
            obj.parse_input_struct(arg_in);         % записать аргументы в объект
            varargout = cell(1,n_options);          % сформировать массив
            for i = 1:n_options                     
                varargout{i} = obj.vars.(obj.option_names{i}); % 
            end
        end
        
        function varargout = parse_input_cell_with_residue(obj, arg_in)
        %PARSE_INPUT_CELL Обработать входные параметры и значения в
        %порядке следования по-умолчанию
        %   Функция принимает на вход массив ячеек, содержащих имена параметров и значения параметров.
        %   Все параметры, не входящие в заранее заданный список сохраняются в переменную объекта residue.
        %   Получить доступ к этим параметрам можно напрямую, либо с помощью функции get_residue.
            obj.keep_residue = 1;
            obj.residue = {};
            n_options = length(obj.option_names);   % общее число опций
            if nargout ~= n_options
                error('Число входных и выходных параметров должно совпадать')
            end
            obj.parse_input_struct(arg_in);         % записать аргументы в объект
            varargout = cell(1,n_options+1);          % сформировать массив
            for i = 1:n_options                     
                varargout{i} = obj.vars.(obj.option_names{i}); % 
            end
            varargout{end} = obj.residue;
        end
        
    end
    
    %%
    methods(Hidden=true)
        function vars = parse_input_struct(obj, arg_in)
        %PARSE_INPUT_STRUCT Обновить объект и сформировать структуру
            n_args = length(arg_in);
            if mod(n_args,2)~=0
                error('Число имен и значений не совпадает')
            end
            for pair = reshape(arg_in,2,[])                 % сделать пары имя-значение по строкам
                if any(strcmp(pair{1},obj.option_names))    % если такая опция предусмотрена
                    if obj.check_type                       % если нужно проверить тип
                        if any(strcmp(class(obj.vars.(pair{1})),obj.types.(pair{1})),'all')
                        else %TODO: проверить правильность работы.
                            warning('Тип аргумента "%s" не совпадает',pair{1}); % дать предупрежение
                            continue;                       % пропустить
                        end
                    end
                    obj.vars.(pair{1}) = pair{2};           % обновить значение опции
                elseif obj.keep_residue
                    obj.residue = cat(2, obj.residue, pair.');
                else 
                    warning(...
                        '"%s" неизвестный параметр. Чтобы подавить сообщение, установите свойство "keep_residue"',...
                        pair{1})
                end
            end
            vars = obj.vars;
        end
    end
    
    methods(Static=true)
        
        % поиск значения именнованного параметра по его имени
        function parameter_value = get_value(arg_in, parameter_name, varargin)
        %get_value Функция выполняет поиск имени параметра и выводит его значение
        %   Аргументы:
        %   arg_in - массив ячеек с последовательно идущими именами и значениями параметров - можно вводить varargin
        %   parameter_name - строковое значение имени параметра, значение которого требуется найти
        %   varargin - третьим элементом может быть задано значение параметра по умолчанию.
        %   KeywordArguments(varargin, 'Parameter_to_find', 69) % 69 будет присвоено, если параметра нет в varargin
            i = find(strcmpi(parameter_name,arg_in),1);
            if isempty(i)
                if nargin == 3
                    parameter_value = varargin{1};
                else 
                    error(['Параметр ' parameter_name ' не найден в вводе'])
                end
            elseif length(i)==1
                parameter_value = arg_in{i+1};    %
            else
                error(['Имя параметра ' parameter_name ' найдено несколько раз.'])
            end
        end
        
    end

end
