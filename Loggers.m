classdef Loggers
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    % saves data to file. If used as object - keeps context in file saved
    %%
    properties
        data_output_file
        context_file
    end
    
    properties(Constant=true)
        
    end
    
    methods
        %%
        
        function obj = untitled2(inputArg1,inputArg2)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
    
    methods(Static=true)
       
        function log_init(filename, varargin)
            % create file
            % check if file csv, mat or log
            header = Loggers.get_value(varargin, 'header',...
                strcat( string(datetime(datetime, 'Format','yyyy-MM-dd HH:mm:ss')), ", file: ", filename));
            idx = strfind(filename,'.');
            extension = filename(idx(end):end);
            switch extension
                case '.mat'
                case '.csv'
                otherwise
                    fileid = fopen(filename, 'w');
                    fprintf(fileid, '%s\n', header);
            end
            fclose(fileid);
        end
        
        function reproduced_values = log(filename, values, varargin)
        %log_log - Записать значения values и names в filename
        % line_ids is vertical cell array (Nx1), num array or empty
        % values - num arry
        % 'header' - заголовок файла. Может быть передана инфомация о значениях
        % 'sep' - разделительный знак
        % 'break' - поставить знак разделения
        % returns values
        % 
            header = Loggers.get_value(varargin, 'header', strcat(...
                string(datetime(datetime, 'Format','yyyy-MM-dd HH:mm:ss')), ", file: ", filename));
            sep = char(Loggers.get_value(varargin, 'sep', ' '));
            reproduced_values = values;
            if isnumeric(header), header = num2str(header,['%d' sep]); end
            if isnumeric(values) 
                char_values = num2str(values, ['%d' sep] );
                char_values = char_values(1:end-1);
            elseif iscell(values)
                char_values = '';
                for i = 1:length(values)
                    if isnumeric(values{i})
                        if size(values{i},1)>1
                            if size(values{i},2)==1,    values{i} = values{i}.';
                            else, error('log: Не поддерживаемый ввод. В ячейке может содержаться вектор строка или столбец.')
                            end
                        end
                        char_values = [char_values num2str(values{i},['%d' sep])];
                    elseif isstring(val), char_values = append(char_values, char(val));
                    end
                end
                char_values=char_values(1:end-length(sep));
            else
                error('log: Не поддерживаемый ввода. Только векторы и массивы ячеек')
            end
            try
                if exist(filename, 'file')             % если файл не существует, создать и записать шапку 
                    fileid = fopen(filename,'a+');
                else
                    fileid = fopen(filename,'w+');
                    fprintf(fileid, '%s\n', header);
                end
                fprintf(fileid, '%s\n', char_values);
                fclose(fileid);
            catch
                fclose(fileid);
            end
        end
        
        function filename_date = stamp_filename(filename)
        %stamp_filename Добавить к имени файла дату и время UTC
            
        end
        
                % поиск значения именнованного параметра по его имени
        function parameter_value = get_value(arg_in, parameter_name, varargin)
        %get_value Функция выполняет поиск имени параметра и выводит его значение
        %   Аргументы:
        %   arg_in - массив ячеек с последовательно идущими именами и значениями параметров - можно вводить varargin
        %   parameter_name - строковое значение имени параметра, значение которого требуется найти
        %   varargin - третьим элементом может быть задано значение параметра по умолчанию.
        %   KeywordArguments(varargin, 'Parameter_to_find', 69) % 69 будет присвоено, если параметра нет в varargin
            i = find(strcmp(parameter_name,arg_in),1);
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

