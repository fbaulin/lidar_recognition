classdef KeywordArguments < handle 
    %KEYWORDARGUMENTS Обработка входный именованных аргументов
    %   
    %   Detailed explanation goes her
    %   Базовое использование класса
    properties
        vars
        check_type = false
        types
        option_names
    end
    
    methods
        function obj = KeywordArguments(varargin)
            %KEYWORDARGUMENTS Construct an instance of this class
            %   передается набор параметров и значений по умолчанию
            %   
            if mod(length(varargin),2)~=0
                error('Число имен должно быть равно числу значений')
            end
            obj.vars = struct(varargin{:}); % создать структуру
            obj.option_names = varargin(1:2:end);
            
        end
        
        function vars = parse_input_struct(obj, arg_in)
            %PARSE_INPUT_STRUCT Обновить объект и сформировать структуру
            n_args = length(arg_in);
            if mod(n_args,2)~=0
                error('Число имен и значений не совпадает')
            end
            
            for pair = reshape(arg_in,2,[])                 % сделать пары имя-значение по строкам
                if any(strcmp(pair{1},obj.option_names))    % если такая опция предусмотрена
                    if obj.check_type                       % если нужно проверить тип
                        if ~strcmp(class(obj.vars.(pair{1})),obj.types.(pair{1}))
                        else 
                            warning('Тип аргумента %s не совпадает',pair{1}); % дать предупрежение
                            continue;                       % пропустить
                        end
                    end
                    obj.vars.(pair{1}) = pair{2};           % обновить значение опции
                else,   warning('%s неизвестный параметр',pair{1})
                end
            end
            vars = obj.vars;
        end
        
        function varargout = parse_input_cell(obj, arg_in)
            %PARSE_INPUT_CELL Обработать входные параметры и значения в
            %порядке следования по-умолчанию
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
        
        function varargout = produce_input(obj)
            varargout = cell(1,2*length(obj.option_names));
            varargout(1:2:end) = obj.option_names(:);
            varargout(2:2:end) = cellfun(@(fn) obj.vars.(fn),obj.option_names,'UniformOutput',false);
        end
        
        function set_types(obj,varargin)
            if length(varargin) == 1
                if strcmp(varargin{1},'as_default')
                    for name = fieldnames(obj.vars)
                        obj.types.(name) = class(obj.vars.(name));
                        obj.check_type = true;
                    end
                end
            end
            warning('Неподдерживаемый формат ввода')
        end
       
    end

end
