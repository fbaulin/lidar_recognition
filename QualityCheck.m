classdef QualityCheck < handle
    %QualityCheck ����� ��� ������ �������� �������������
    %   ��������� ������ �������� �������������
    
    properties
        reduction_method    % ����� �������� ������������ ��������� (�������, �����������, etc)
        heuristic           % ����� ���������� ����� ����������
        classifier_type     % �������� �������������� 
        dimensions          % ����������� ������������ ���������
        clf                 % �������������, �������������� configure � train
        decision_mode       % ������� ������ ��������� ����������
        use_GPU = 'no'      % 
        fs_map              % ����� ���������
        score_type          % ��� ���� �������� �������������
        
        maps_filename       % ��� ����� � ������� ����������� ���������
        
        batch_filename      % ��� ����� ��� ���������� 
        feature_header      % ��������� ��� �����
        
        n_nearest           % ����� ��������� �������
        m_technique         % ��� ������� �������� �������
        k_alien             % ��������� ����� �����
        n_metrics_fraction  % ���� ������ ��� ������� ��������� ���� ������ ��������� ��������
        hist_edges          % ������� ���������� ��� ������ ������ �� ������������
    end

   
    methods
        
        % ����������� ������
        function obj = QualityCheck(varargin)
        %QualityCheck ������������ ���������� ������
        %   ������������ ���������� ������
        %   Args
        %       'reduction_method'  - ����� �������� [nmin, buryi, hist].
        %       'classifier_type'   - ��� ��������������, 'mlp' - ������������ ����������.
        %       'layers',           - ������ ��������� ������� �����, [160 80 40].
        %       'dimensions'        - �����������, 1.
        %       'batch_filename'    - ��� ����� ��������� �������, false.
        %       'maps_filename'     - ��� ����� � ������� ���������,  false.
            kwargs = KeywordArguments(...
                'reduction_method','nmin',...
                'heuristic', 'adddel',...   % ����� ��������
                'classifier_type', 'none',...
                'dimensions', 1, ...
                'batch_filename', false,...
                'maps_filename', false,... 
                'n_nearest',5, ...          % ����� ����������� �������
                'm_technique','mean',...    % ����� �������
                'k_alien',2, ...               % ����� �����
                'n_metrics_fraction', 0.05,...
                'hist_edges', linspace(0,2,51)...  
                );   % �������� ��-���������
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
                
        % ��������� ��������������
        function setup_classifier(obj, varargin)
        %SETUP_CLASSIFIER ��������� ��������������
        %   Args, default:
        %       'classifier_type', 'mlp'    - ��� ��������������    
        %       'layers', [40 20 10]        - ����� ����
        %       'performance', 'mse'        - ������� �� ������� ������������� �������� ������
        %       'train_fun', 'default'      - ����� ��������
        %       'gpu', obj.use_GPU          - ������������ ��� ��� ��������
            kwargs = KeywordArguments(...
                'classifier_type', 'mlp',...    
                'layers', [40 20 10],...
                'performance', 'default',...
                'decision_mode', 'metric',...   % ���� 'metric', �� ������� ����������� �� �������, ����� - �� 
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
                    if strcmp(train_fun, 'default') % ���� �������� ��-���������, ��
                        if strcmp(obj.use_GPU, 'yes')
                            train_fun = 'trainscg'; % ��� ��� - �����������
                        else
                            train_fun = 'trainlm';  % ��� ��� - ���������-����.
                        end
                    end
                    obj.clf = feedforwardnet(layers,train_fun); % ������� ����������
                    if strcmp(performance,'default'), performance='mse'; end
                    obj.clf.performFcn = performance;           % ������� ������
                    fprintf('\t���������: [%s], ��������: %s\n',...
                        num2str(layers,'%d '), performance);
                case 'pnet'
                    if strcmp(performance,'default'), performance='crossentropy'; end
                    if strcmp(train_fun, 'default'), train_fun = 'trainscg'; end
                    obj.clf = patternnet(layers,train_fun,performance);       % ������� ����������
                    fprintf('\t������������ �������������: [%s], ��������: %s\n',...
                        num2str(layers,'%d '), performance);
                otherwise               % ���� �������� ������������� �� ����������
                    error('��� �������������� �� ��������������')
            end
            obj.clf.trainParam.showWindow = false;      % ��������� ���� �� ����. ��������
            obj.clf.trainParam.showCommandLine = false; % ��������� ����� � �������
            
        end
        
        % ���������� ��������� �������
        function [ x_train, y_train, x_test, y_test ] = prepare_dataset(obj, x_train, y_train_str, x_test, y_test_str)
            %NOTE: ��� ����������� ������ ��������� ����� ���������� ������
            % ������ �� ������������� ����������� ����, ��� ������ � ������
            % ���������� ���������� ����� ������ ������. ���
            % ��������������� � ���, ��� ��� ���������� ������ ��������� ��
            % ��������� nmin �����, ����� ��� ����������� ���� ���� ������
            % ��� ���� ��������� (�� ������ � ��� ����� ���������)
            n_test = size(y_test_str,1);
            n_train = size(y_train_str,1);
            y = vertcat(y_train_str, y_test_str);       % ������� ����� ������� ��� ����������� ������������ �������� �������
            y = obj.string2vector(y);                   % ������������� �������� ������� � �������
            y_train = y(1:size(y_train_str, 1),:);      % �������� �� ����������� ������� ��������� �������
            y_test =  y(size(y_train_str, 1)+1:end,:);      % �������� �������� �������
            
            data_means = mean(x_train,1);       % ������� ��������
            data_stds = std(x_train,[],1);      % ������ ���
            
            x_train =  (x_train - data_means)...
                ./repmat(data_stds,n_train,1);         % ����������� �������� ������� ���������;
            x_test = (x_test - data_means)...
                ./repmat(data_stds,n_test,1);          % ����������� �������� ������� ���������
            
        end
        
        % ������� �������� �������������
        function [score, varargout] = assess_quality(obj, x_train, y_train_str, x_test, y_test_str)
        %assess_quality ������� �������� ��������� ��������� � �������������
        %   ��������� ��������� ��������� � ������ �������� ������������� ��� ������������
        %   �������� � ���������� ������.
        %   Args:
        %       x_train     - ������� ��������� ��������� �������.
        %       y_train     - ������� (�������� ��������) ��������� �������.
        %       x_test      - ������� ��������� �������� �������.
        %       y_test      - ������� (�������� ��������) �������� �������.
        %   Returns: 
        %       score       - ���� �������� �������������, ��������������� ��� �������������.
        %       [rho_min]   - ���� �������� ��������������� ������������ ���������. 
        %       [fs_maps]   - ����� ���������������� ������������ ��������� <n_dim x fv_dim>.
        %       [conf_mx]   - ������� ������ ��� ����������� �����������.

%             fvs = obj.wrap_data(x_train,y_train_str);   % ������������� ������ � ������ �����
%             [rho_min, fs_maps,obj.dimensions] = SystemModel.fs_reduction(fvs,...
%                 'dimensions',obj.dimensions, 'reduction_type',obj.reduction_method...
%                 );   % ��������� ����������� ������������ ��������� ��� �������� �����������

            reduction = RecursiveReduction();
            [rho_min, fs_maps,obj.dimensions] = reduction.reduce(x_train, y_train_str,...
                'dimensions',obj.dimensions, 'reduction_type',obj.reduction_method,...
                'heuristic',obj.heuristic,...
                'n_metrics_fraction',obj.n_metrics_fraction,...
                'hist_edges', obj.hist_edges,...
                'n_nearest', obj.n_nearest,...
                'm_technique', obj.m_technique,...
                'k_alien',obj.k_alien);   % ��������� ����������� ������������ ��������� ��� �������� �����������
            
            obj.fs_map = logical(fs_maps);      % ������������� � ���������� ��� ����������
            n_dims = length(rho_min);           % �������� ����� ������������� ������������
            score = zeros(n_dims,1);            % ������-������� � ������� ���-���� �������������
            if nargout>=4, conf_mx = cell(n_dims,1); end        % ������� ������
%             conf_mx = cell(length(fvs), length(fvs));
            % ���������� ������: � �.�. ������������
            [ x_train, y_train, x_test, y_test ] = obj.prepare_dataset(x_train, y_train_str, x_test, y_test_str);
            
            if strcmp(obj.use_GPU,'yes') % �������� �� �� �� �������� � single ����������
                x_train = double(x_train);
                x_test = double(x_test);
            end
            n_arg_out = nargout;
            parfor i_dim = 1:n_dims
                x_train_red = x_train(:, obj.fs_map(i_dim,:));   %#ok<PFBNS> % ������������ ��������� �������
                x_test_red = x_test(:, obj.fs_map(i_dim,:));     %#ok<PFBNS> % ������������ �������� �������
               % ������� �������������
                inst_clf = configure(obj.clf,x_train_red.',y_train.');   % ����. �������� � ��������� ����
                inst_clf = train(inst_clf,x_train_red.',y_train.','useGPU',obj.use_GPU);       % ��������
                y_pred = inst_clf(x_test_red.').';                       % ������������ ������� ��� �������� �������
                score(i_dim) = perform(inst_clf, y_test.', y_pred.');    % ������ �������� �������������
                if n_arg_out>=4
                    switch obj.decision_mode
                        case 'metric', y_pred = QualityCheck.metric_decision(y_pred);
                        case 'max'
                        otherwise, error(['����� ����������� ����� �������� ������� (', obj.decision_mode, ')'])
                    end
                    [~,i_tst] = max(y_test,[],2); [~,i_pred] = max(y_pred,[],2);    % one-hot -> ������� �������
                    conf_mx{i_dim} = confusionmat(i_tst,i_pred);
                end    % ��������� �������
            end
            if nargout>=2, varargout{1} = rho_min; end      % ���� ����� ������ ��� - ������ rho_min
            if nargout>=3, varargout{2} = fs_maps; end      % ���� 3 ��������� ������ ����� ���������
            if nargout>=4, varargout{3} = permute(cat(3,conf_mx{:}),[3 1 2]); end      % ������ ������� ������
        end
        
        
        
    end
    
    methods(Access = private, Hidden = true)
        
        % ��������� ������� � csv
        function [] = save_batch(obj, features, meta)
        %SAVE_CSV ��������� � csv ��������
        %   � ������ ������ ���������� �������� ���������, � t ����������
        %   ������. ���, ��� ���������. 
        %   Args:
        %       features    - ������� ���������.
        %       meta        - ��������� ����������� ��������� (��� ������� [������]).
        %       f_header    - ���������.
        %       filename    - ��� �����.
            f_header = obj.feature_header;
            filename = obj.batch_filename;
            [ n_obj, ~ ] = size(features);
            %meta = RPTools.meta2str(meta); % ������������ ���� ����������� ����� �
            %���������
            %meta = {meta.name}; % 
            fid = fopen(filename,'w');
            fprintf(fid,'%s,',f_header);
            fprintf(fid,'%s\n','t');
            for i_line = 1:n_obj
                fprintf(fid,'%f,', features(i_line,:));
                fprintf(fid,'%s\n', meta(i_line).name);
            end
            fclose(fid);
            
            disp(['���� ' filename ' ��������'])
            
        end
        
    end
    
    methods(Static=true)
        
        function y_dec = metric_decision(y_pred, varargin)
        %   
        %   ���� ���������� �� �������� ������� (one hot �������) ������ ���������� (�� ��������� (sqrt(2)/2), ��
        %   ������ ���������� �� ��������������� one-hot ������, ����� - �� ������� ������.
            if nargin==2, thresh_sqr = sqrt(varargin{1}); else thresh_sqr = 1/2; end
            [n_samples, n_class] = size(y_pred);    % ���������� ����� �������
            y_dec = zeros(n_samples, n_class, 'like', y_pred);
            % ���� �� ������� 
            for i=1:n_class
                y_h = zeros(n_samples, n_class, 'like', y_pred);
                y_h(:,i) = 1;
                dists = sum((y_h-y_pred).^2, 2);    % ���������� �������� ������ �� �������� ������� ������
                y_dec(dists<thresh_sqr, i) = 1;            % ������ ��� ������� ��� ������� ���������� ������ ������
            end
            y_dec(sum(y_dec,2)>1,:) = 0; % ������ �������, ��� ������� ����� ������ ������� 
        end
        
        % ������������� ������� �� ������ ����� � ������ ������
        function [x,y] = unwrap_cell_data(fvs,meta,varargin)
        %unwrap_cell_data ������������� ������ ����� � ������
        %   ��������� ������, �������� �� ������� � ������� ��������� � ������ ����� � ������
        %   ������ ��� ������� �������� �������.
        %   Args:
        %       fvs     - ������ ����� � ��������� ���������.
        %       meta    - ������ ����� �� ���������� � ������� �������.
        %       'matrix_output' - ������ ������� ����� �������,  false. 
            kwargs = KeywordArguments('matrix_output', false );
            [ matrix_output ] = kwargs.parse_input_cell(varargin);
            n_obj = length(fvs);    % ����� ��������
            x = vertcat(fvs{:});    % �������� ��� ������� ��������� � ���� �������
            y = arrayfun(@(i) {meta{i}.name}.', 1:n_obj, 'UniformOutput', false); % ������� �������� �������
            y = vertcat(y{:});  % �������� ������� (�� ��������) �������� �����, � ���� ������ �����
            if matrix_output,   y = cellfun(@(name) string(name), y);    end     % �������� ������ ����� � �������
        end
        
        % ������������� ������ �� ������� � ������ ����� �� ��������
        function [varargout] = wrap_data(x,y)
        %wrap_data ������������� ������ � ��������� � ������.
        %   ��������� ������ x � ������� �����, ��������������� ��������
        %   ���������� � �������� ������� �� y.
        %   Args:
        %       x - ������ � ��������� ���������.
        %       y - ������ ����� � ������� ��������.
        %   Returns:
        %       ������ �����, ���������� ������� �������� ���������, ����� ����� � ������� ����� ����� ��������.
        %       ���� ������������ 2 ����������, �� ������ - ����� ��������.
            [object_names, ~, ci] = unique(y);  % ������� ����� ��������
            % ����� ci - ������ �������, ������������� i-�� ����� � ������� object_names
            n_objects = length(object_names);   % ����� ��������
            fvs = cell(n_objects, 1);           % ������������ ������ � ������� ��� ������ ������
            for i_obj = 1:n_objects             % ���� �� ��������
                fvs{i_obj} = x(i_obj==ci,:);    % ����� ������� ��������� i-�� ������� � ��������� ��
            end
            varargout{1} = fvs;
            if nargout>1, varargout{2} = object_names; end
                
        end
        
        % ������������ ������� one-hot �������
        function [y_vect] = string2vector(y_string)
        %string2vector  ������������ ������� ������� �� ��������� �������
            [object_names, ~, ci] = unique(y_string);           % �������� ����� ��������
            y_vect = zeros(length(ci), length(object_names));   % ������������ ������� ��� ������� �������
            for i_obj = 1:length(object_names)  % ���� �� ��������
                y_vect(i_obj==ci,i_obj) = 1;    % � �������� ������ �������� ��������� �������
            end 
        end
        
    end
    
end

