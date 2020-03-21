classdef KeywordArguments < handle 
    %KEYWORDARGUMENTS ��������� ������� ����������� ����������
    %   
    %   Detailed explanation goes her
    %   ������� ������������� ������
    properties
        vars
        check_type = false
        types
        option_names
    end
    
    methods
        function obj = KeywordArguments(varargin)
            %KEYWORDARGUMENTS Construct an instance of this class
            %   ���������� ����� ���������� � �������� �� ���������
            %   
            if mod(length(varargin),2)~=0
                error('����� ���� ������ ���� ����� ����� ��������')
            end
            obj.vars = struct(varargin{:}); % ������� ���������
            obj.option_names = varargin(1:2:end);
            
        end
        
        function vars = parse_input_struct(obj, arg_in)
            %PARSE_INPUT_STRUCT �������� ������ � ������������ ���������
            n_args = length(arg_in);
            if mod(n_args,2)~=0
                error('����� ���� � �������� �� ���������')
            end
            
            for pair = reshape(arg_in,2,[])                 % ������� ���� ���-�������� �� �������
                if any(strcmp(pair{1},obj.option_names))    % ���� ����� ����� �������������
                    if obj.check_type                       % ���� ����� ��������� ���
                        if ~strcmp(class(obj.vars.(pair{1})),obj.types.(pair{1}))
                        else 
                            warning('��� ��������� %s �� ���������',pair{1}); % ���� �������������
                            continue;                       % ����������
                        end
                    end
                    obj.vars.(pair{1}) = pair{2};           % �������� �������� �����
                else,   warning('%s ����������� ��������',pair{1})
                end
            end
            vars = obj.vars;
        end
        
        function varargout = parse_input_cell(obj, arg_in)
            %PARSE_INPUT_CELL ���������� ������� ��������� � �������� �
            %������� ���������� ��-���������
            n_options = length(obj.option_names);   % ����� ����� �����
            if nargout ~= n_options
                error('����� ������� � �������� ���������� ������ ���������')
            end
            obj.parse_input_struct(arg_in);         % �������� ��������� � ������
            varargout = cell(1,n_options);          % ������������ ������
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
            warning('���������������� ������ �����')
        end
       
    end

end
