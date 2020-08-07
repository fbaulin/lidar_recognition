classdef KeywordArguments < handle 
    %KEYWORDARGUMENTS V0.0 ��������� ������� ����������� ����������
    %TODO: �������� ����� extract_value   
    %   ����� ������������ ��� ��������� ������� ����������� ����������.
    %   ��� ���������� �������, � ������� ������������ ���� ����������� ����������, �������� ����������
    %   �� ��������� ����� ���� ������ �������. ������������� � �������� ������� �������� ���������������
    %   ��� �����: 
    %       1) �������� ������� � �������� ���������� � �� �������� ��-��������� 
    %           ������� ���������� ���������� �� ����� �������� - �������, ����� ������ ��� �������� ���������,
    %           � �� ��� - ��� ��������.
    %       kwargs = KeywordArguments(...                               % � �������� �������� ����� ������������
    %           'Param_name', 'some default value',...                  % ������
    %           'Other_name', 2, ...                                    % �����, ������� � �������
    %           'cell parameter', {'a', 2, some_value_calced earlier},  % ������ ��� �������� �������� ����
    %           'fun', @(a) a^2);                                        % ������-�������
    %       2) ��������� ������� ���������� varargin � �������������� ���������������� �������:
    %       [param_internal_name, other_par, cell_par, fun] = kwargs.parse(varargin);
    %       
    %       ���� ������� �� ��������� ����� ����������, � �������� �� � ���� �� ���������� �� �������, �� ������� 
    %       ��������������� �������� ���������� �������. ��� ������� ��������� ��������� ���������, �� ��������
    %       ��� ������������� �������. ��� ������������� ���� ������� ������� �������� �������� keep_residue.
    %       ����� ��������� ��������� keep_residue "����������������" ��������� �� varargin ����� ��������� � �������.
    %       ������ � ��� ����� ��������, ����������� � ���� residue ��� � �������������� ������� get_residue.
    %       ���� �������� ���������� ������������� ��� ������ ������� parse_with_residue(). ���������, ���������� 
    %       � ����������, ������������ ����� ������ ����� � ����� ���� �������� � �������� �������, ����.:    
    %       kwargs = KeywordArguments(__);
    %       kwargs.keep_residue = true;
    %       [__] = kwargs.parse(varargin);
    %       some_subfun(kwargs.residue{:});
    %                                   
    %%
    properties
        vars
        check_type = false      %TODO: �������� ������� � ��������� ��� ������� ��������� - ������� ��������
        types                   % ���� ����������
        option_names            % ����� ���������� � ������� �������� ��� ������������� �������
        keep_residue = false    % ��������� ��������, �� �������� ��� ������������� �������
                                % ���� false, �� ��� �������� ����������������� ���������� �������� ������
        residue                 % ����� �������� ���������� � varargin �� ��������������� ���������
    end
    
    %%
    methods
        
        % �����������, ���������������� ������
        function obj = KeywordArguments(varargin)
        %KEYWORDARGUMENTS Construct an instance of this class
        %   ���������� ����� ���������� � �������� �� ���������
        %   
            if mod(length(varargin),2)~=0
                error('����� ���� ������ ���� ����� ����� ��������')
            end
            obj.vars = struct(varargin{:});         % ������� ���������
            obj.option_names = varargin(1:2:end);
            
        end
               
        % ���������� varargin,
        function varargout = parse(obj, arg_in)
        %PARSE_INPUT_CELL ���������� ������� ��������� � �������� �
        %������� ���������� ��-���������
            obj.residue = {};
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
        
        % ���������� varargin, ������ ����������������� ���������
        function varargout = parse_with_residue(obj, arg_in)
        %PARSE_INPUT_CELL ���������� ������� ��������� � �������� �
        %������� ���������� ��-���������
        %   ������� ��������� �� ���� ������ �����, ���������� ����� ���������� � �������� ����������.
        %   ��� ���������, �� �������� � ������� �������� ������ ����������� � ���������� ������� residue
        %   � ������������ ���� �������� ��������� ���������, ���������� �������� �����.
            obj.keep_residue = 1;
            obj.residue = {};
            n_options = length(obj.option_names);   % ����� ����� �����
            if nargout ~= n_options
                error('����� ������� � �������� ���������� ������ ���������')
            end
            obj.parse_input_struct(arg_in);         % �������� ��������� � ������
            varargout = cell(1,n_options+1);          % ������������ ������
            for i = 1:n_options                     
                varargout{i} = obj.vars.(obj.option_names{i}); % 
            end
            varargout{end} = obj.residue;
        end
        
        % ���������� varargin, ������ ����������������� ���������
        function varargout = get_residue(obj)
        %get_residue ���������� ��������� �� varargin, �� �������� � ������ ������� ��������.
            varargout = obj.residue;
        end
        
        % ��������� ���� ���-��������, �� � ������� ����������, ������������� ��� ������������� �������
        function varargout = produce_input(obj)
        %reproduce_input ���������� ����� ���������� � ���� ������� �����
        %   ���������� ����� ���������� � ���� ������� �����, ����������� �� �������� ������
        %   �������� ����������, � �� ������ - �� ��������.
            varargout = cell(1,2*length(obj.option_names));     % ������� ������ �����
            varargout(1:2:end) = obj.option_names(:);           % ��������� �������� ����������
            varargout(2:2:end) = cellfun(@(fn) obj.vars.(fn),...
                obj.option_names,'UniformOutput',false);        % ��������� �������� ����������
        end
                
        % ������ �������� 
        function set_types(obj,varargin)
            if length(varargin) == 1
                if strcmp(varargin{1},'default')
                    for name = fieldnames(obj.vars)
                        obj.types.(name) = class(obj.vars.(name));
                    end
                    obj.check_type = true;  %TODO: ������� �������� ���������
                else
                    warning('���������������� ������ �����')
                end
            elseif mod(length(varargin),2)==0 
                for pair = reshape(varargin,2,[])
                    obj.types.(pair{1}) = pair{2};
                end
            else
                warning('�������� ����� ���������. � ������� set_types ������� ���������� ����.')
            end
            
        end        
        
        %% legacy v0.0
        function varargout = parse_input_cell(obj, arg_in)
        %PARSE_INPUT_CELL ���������� ������� ��������� � �������� �
        %������� ���������� ��-���������
            obj.residue = {};
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
        
        function varargout = parse_input_cell_with_residue(obj, arg_in)
        %PARSE_INPUT_CELL ���������� ������� ��������� � �������� �
        %������� ���������� ��-���������
        %   ������� ��������� �� ���� ������ �����, ���������� ����� ���������� � �������� ����������.
        %   ��� ���������, �� �������� � ������� �������� ������ ����������� � ���������� ������� residue.
        %   �������� ������ � ���� ���������� ����� ��������, ���� � ������� ������� get_residue.
            obj.keep_residue = 1;
            obj.residue = {};
            n_options = length(obj.option_names);   % ����� ����� �����
            if nargout ~= n_options
                error('����� ������� � �������� ���������� ������ ���������')
            end
            obj.parse_input_struct(arg_in);         % �������� ��������� � ������
            varargout = cell(1,n_options+1);          % ������������ ������
            for i = 1:n_options                     
                varargout{i} = obj.vars.(obj.option_names{i}); % 
            end
            varargout{end} = obj.residue;
        end
        
    end
    
    %%
    methods(Hidden=true)
        function vars = parse_input_struct(obj, arg_in)
        %PARSE_INPUT_STRUCT �������� ������ � ������������ ���������
            n_args = length(arg_in);
            if mod(n_args,2)~=0
                error('����� ���� � �������� �� ���������')
            end
            for pair = reshape(arg_in,2,[])                 % ������� ���� ���-�������� �� �������
                if any(strcmp(pair{1},obj.option_names))    % ���� ����� ����� �������������
                    if obj.check_type                       % ���� ����� ��������� ���
                        if any(strcmp(class(obj.vars.(pair{1})),obj.types.(pair{1})),'all')
                        else %TODO: ��������� ������������ ������.
                            warning('��� ��������� "%s" �� ���������',pair{1}); % ���� �������������
                            continue;                       % ����������
                        end
                    end
                    obj.vars.(pair{1}) = pair{2};           % �������� �������� �����
                elseif obj.keep_residue
                    obj.residue = cat(2, obj.residue, pair.');
                else 
                    warning(...
                        '"%s" ����������� ��������. ����� �������� ���������, ���������� �������� "keep_residue"',...
                        pair{1})
                end
            end
            vars = obj.vars;
        end
    end
    
    methods(Static=true)
        
        % ����� �������� ������������� ��������� �� ��� �����
        function parameter_value = get_value(arg_in, parameter_name, varargin)
        %get_value ������� ��������� ����� ����� ��������� � ������� ��� ��������
        %   ���������:
        %   arg_in - ������ ����� � ��������������� ������� ������� � ���������� ���������� - ����� ������� varargin
        %   parameter_name - ��������� �������� ����� ���������, �������� �������� ��������� �����
        %   varargin - ������� ��������� ����� ���� ������ �������� ��������� �� ���������.
        %   KeywordArguments(varargin, 'Parameter_to_find', 69) % 69 ����� ���������, ���� ��������� ��� � varargin
            i = find(strcmp(parameter_name,arg_in),1);
            if isempty(i)
                if nargin == 3
                    parameter_value = varargin{1};
                else 
                    error(['�������� ' parameter_name ' �� ������ � �����'])
                end
            elseif length(i)==1
                parameter_value = arg_in{i+1};    %
            else
                error(['��� ��������� ' parameter_name ' ������� ��������� ���.'])
            end
        end
        
    end

end
