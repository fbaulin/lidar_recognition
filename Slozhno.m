classdef Slozhno < handle
    %Slozhno Support functions for slozhno experiment
    %   Slozhno experiment
    %#ok<*PROPLC>
    properties
                
        window_m = 3;           % ������������ ������������� �������
        icr_dim = 200;          % ����� �������� � ���
        hwidth2sig = 2.35;      % ��������� ��� �������� ����� � ����� �� ����������� ������
        resol_ns                % �������� ���������� ���
        rp_dim
        
        scan_type               % ��� ������������ ��������������� �������������
        obj_lst                 % ������ ��������
        imp_hwidths_ns          % ������ ������������� �� �� ���������� ������
        t_type                  % ��� ������������� ��������������
        wave_name               % ��� ��������
        cwave_name
        band_factors            % ��������� fd/�
        
        f_sampling_mhz          % ������� ������������, ���
        T_sampling_ns           % �������� �������������, ��
        f_imp_mhz               % ������ �������, ���
        f_lpfilter_mhz          % ������ �������, ���
    end
        
    methods
        
        % �����������
        function obj = Slozhno() 
            
            obj.resol_ns = obj.window_m * 2 / obj.icr_dim / 3e8 * 1e9;      % �������� ��������� ���������� ���, ��
            obj.obj_lst ={ 'brk3x1' 'con2x1','con3x1','cyl3x1','sph3x1' };
            obj.scan_type = 'default';      % �� ��������� 
            obj.t_type = 'none';            % �� ��������� ��� ��������������
            obj.wave_name = 'dtf2';         % �� ��������� ����������� �������
        end
        
        % ��������� ��� ������ ������� �� icr ������
        function [icrs, meta] = get_icrs(obj, icr_folder) 
        %GET_ICRS ��������� ��� � ���������� (���� � ��.) �� �����
            %   ���������:
            %   icr_folder -    - ����� � ���
            %   obj_id          - ��� ��� ��
            % �������� ���
            mat_files = dir([icr_folder '/*.mat']);     % �������� ����� ���
            if size(mat_files,1)~=0     % ���� ���-�����              
                load([mat_files(1).folder '\' mat_files(1).name],'icrs','meta');   % ��������� ��� �� mat �����
            else
                filelist = dir([icr_folder '/*.icr']);      % �������� ��������� ������ ������ �� �������� �����
                for i = 1:length(filelist)                  % ������� ������ ������ ���� ������ �� ���������
                    filelist(i).name = horzcat(filelist(i).folder, '\', filelist(i).name); 
                end
                [icrs, meta] = RPTools.open_icr_file(...
                    {filelist(1:end).name}.');              % �������� ��� � ����������
                if size(icrs,2)~=obj.icr_dim
                    warning(['����� �������� � ��� �� ����� ' num2str(obj.icr_dim)]);
                end
                icrs_db_filename = [filelist(i).folder '\' meta(1).name '.mat'];
                save(icrs_db_filename,'icrs','meta');
                disp(['���������� ���� ��� � ' icrs_db_filename]);
            end
        end
        
        % �������� ��� ���������� �������� �� mat ������
        function [icrs, meta] = get_all_icrs(obj)
            n_obj = length(obj.obj_lst);    
            icrs = cell(1,n_obj);   % ������� ��� ���
            meta = cell(1, n_obj);  % ������� ��� ����������
            for i_obj = 1:n_obj     % �������� ��� 
                icr_folder = [ 'icrs/' obj.scan_type '/' obj.obj_lst{i_obj}];       % ������ ���� ����� � ���
                %icrs_data = load([icr_folder '/' obj.obj_lst{i_obj} '.mat']);   % ��������� ��� �� mat �����
                [icrs{i_obj}, meta{i_obj}] = obj.get_icrs(icr_folder);
                
            end
        end
        
        % ������ ��������� ������
        function rps = lls_model(obj, icrs, imp_hwidth_ns, lp_band_mhz, fd_mhz_q) 
        %LLS_MODEL ������������� ������ ���
        % icrs - ���, ������������� �� �����������
        % imp_hwidth_ns - ������������ �������� �� �������� ������
        % lp_band_mhz - ������ �� �������
        % fd_mhz_q - ������ ����� {(fd, ���),(��������)}, 
        % 'n' ��������� ������������� � �����������
            % ������������ ��
            
            fprintf('\tllis_model: ')
            imp_hwidth = ceil(imp_hwidth_ns/obj.resol_ns);  % ������������ ��������
            if ~strcmp(imp_hwidth_ns,'n')
                if imp_hwidth_ns>0.2                            % ���� �� ���������� ������������
                    rps = RPTools.make_rps(icrs, imp_hwidth);   % ������ ��
                    fprintf('[����. �� � t=%0.2f��]',imp_hwidth_ns)
                else                    % ����� ������������ �� ������� ����:
                    rps = icrs;         % ������� ��� �� ��
                    fprintf('[�� - ���]');
                end
            else
                rps = icrs;
                fprintf('[����. ��]');
            end
            [n_rps, obj.rp_dim] = size(rps);
            % ����������
            if ~strcmp(lp_band_mhz,'n')                         % ���� �������� 
                n_filter = 2;
                rps = Slozhno.bworth_filter(...
                    rps, obj.resol_ns, lp_band_mhz, n_filter);  % ����������
                fprintf('[���(%d�) F�=%0.1f���]', ...
                    n_filter, lp_band_mhz);
                obj.f_lpfilter_mhz = lp_band_mhz;
            end
            % ���
            if iscell(fd_mhz_q)     % ������������� � �����������
                fd_mhz  = fd_mhz_q{1};  % ���������� ������� �������������
                q_bit   = fd_mhz_q{2};  % ���������� ����� ������� �����������
                if ~strcmp(fd_mhz,'n')  % ������������� �������
                % ����������� ����������
                    fd_src_mhz = 1/obj.resol_ns*1e3;    % �������� ������� �������. (���)
                    rps = RPTools.spline_interpol(rps, fd_src_mhz, fd_mhz);
                    fprintf('[��� F�=%d���', floor(fd_mhz));
                    obj.f_sampling_mhz = fd_mhz;        % �������� ������� �������������
                    obj.T_sampling_ns = 1/fd_mhz*1e3;   % �������� ���������� �� �������
                    
                end
                if ~strcmp(q_bit,'n')   % ����������� 
                    fprintf(',q=%d', q_bit);
                    rps = RPTools.adc(rps,q_bit);                    
                end
                fprintf(']');
            else, warning('    ������� �������� ���');    
            end
            if mod(size(rps,2),2)==1;rps=[rps, zeros(n_rps,1)];end
            fprintf('\n');
            %rps = RPTools.norm_rp(rps,'mod');  % ���������� �� ������
        end
        
        % ����������� ������
        function [bms_counts] = gpu_histcounts(~, fvs, hist_bins) 
        %gpu_histcounts
            n_obj = length(fvs);
            bms_mod = cell(n_obj, n_obj);       % ������� ��� �������
            bms_counts = cell(n_obj, n_obj);    % ������� �������� ������� ������
            gpu_hist_bins = gpuArray(hist_bins);
            for i_obj = 1:(n_obj-1)             % ������ ������
                fvs_1 = gpuArray(fvs{i_obj});   % �������� ������� ��������� � ��
                %dim_fv = size(fvs_1,2);
                fv_mods = sqrt(sum(fvs_1.^2,2)); %repmat(sqrt(sum(fvs_1.^2,2)), 1, dim_fv);     % ��������� ���������� �������� �� �� �����
                fvs_1 = fvs_1./fv_mods;                   % ���������� ��� ������� ������� ������
                for j_obj = (i_obj+1):n_obj
                    fvs_2 = gpuArray(fvs{j_obj});
                    fv_mods = sqrt(sum(fvs_2.^2,2));%repmat(sqrt(sum(fvs_2.^2,2)), 1, dim_fv);
                    fvs_2 = fvs_2./fv_mods; % ���������� ��� ������� ������� ������
                    clear fv_mods
                    bms_mod{i_obj,j_obj} = 2*asin(...
                        Slozhno.calc_b_metrics(fvs_1,fvs_2,8)/2); % ������ ������� ������
                    bms_counts{i_obj,j_obj} = ...
                        histcounts(bms_mod{i_obj,j_obj}, gpu_hist_bins );
                    bms_counts{i_obj,j_obj}=gather(bms_counts{i_obj,j_obj});
                    
                end
            end
        end
        
        % ����������� ������
        function [bms_counts] = gpu_linear_histcounts(~, fvs, hist_bins) 
        %gpu_histcounts
            n_obj = length(fvs);
            bms_mod = cell(n_obj, n_obj);       % ������� ��� �������
            bms_counts = cell(n_obj, n_obj);    % ������� �������� ������� ������
            gpu_hist_bins = gpuArray(hist_bins);
            for i_obj = 1:(n_obj-1)             % ������ ������
                fvs_1 = gpuArray(fvs{i_obj});   % �������� ������� ��������� � ��
                for j_obj = (i_obj+1):n_obj
                    fvs_2 = gpuArray(fvs{j_obj});
                    bms_mod{i_obj,j_obj} = Slozhno.calc_b_metrics(fvs_1,fvs_2); % ������ ������� ������
                    bms_counts{i_obj,j_obj} = ...
                        histcounts(bms_mod{i_obj,j_obj}, gpu_hist_bins );
                    bms_counts{i_obj,j_obj}=gather(bms_counts{i_obj,j_obj});
                    
                end
            end
        end
        
        % ������� �� �����������
        function [] = gpu_bms_vectors(obj, fvs, meta, imp_hwidth_ns) 
        %GPU_BMS_VECTORS ������������ � ���������� ���������� �������� ��
        %��
        %   gpu_bms_vectors(fvs, meta)
        %   fvs - ������ � ����������� ��
        %   meta - ������ � ����������� ��� ����������
            n_obj = length(obj.obj_lst);
            for i_obj = 1:(n_obj-1)
                fvs_1 = gpuArray(fvs{i_obj});   % �������� ������� ��������� � ��
                for j_obj = i_obj+1:n_obj
                    fvs_2 = gpuArray(fvs{j_obj});
                    obj_names = {obj.obj_lst{i_obj} obj.obj_lst{j_obj}};
                    bms = gather(Slozhno.calc_b_vectors(fvs_1, fvs_2));
                    save([...
                        '.\bms_vectors\' obj.scan_type ' '...
                        obj.obj_lst{i_obj} '-' obj.obj_lst{j_obj} ' ',...
                        't=' num2str(imp_hwidth_ns) 'ns '...
                        datestr(fix(clock),'yyyymmdd HH-MM-SS'),...
                        '.mat'], 'obj_names', 'meta' , 'fvs' , 'bms');
                end
            end
        end
       
        % ��������������
        function [ fvs ] = transform(obj, rps, varargin)
        %TRANSFORM ������������ ������� ���������
        %   ��������� �������������� � �������������� ������� ������ RPTools
        %   ���������:
        %   .obj        
        %   rps         - ������� ����� ��
        %   [t_type]    - ��� ��������������
        %   [wave_name] - ����������� ��������
            if nargin>=3;   t_type = varargin{1};   % ��� ��������������
            else; t_type = obj.t_type; end 
            if nargin>=4; wave_name = varargin{2};  % ��� ��������
            else; wave_name = obj.wave_name; end
            switch t_type
                case 'none'
                    tr = @(rp) rp;
                case 'afft' 
                    tr = @(rp) RPTools.afft(rp);
                case 'fft'
                    tr = @(rp) RPTools.fft(rp);
                case 'cwt'  
                    tr = @(rp) RPTools.cwt(rp,wave_name);
                case 'acwt' 
                    tr = @(rp) RPTools.acwt(rp,wave_name);
                case 'wt'   
                    tr = @(rp) RPTools.wt(rp,wave_name);
                case 'pca'   
                    tr = @(rp) RPTools.pca(rp);
            end
            if iscell(rps); fvs = cellfun(tr,rps,'UniformOutput',0);
            else;           fvs = tr(rps); 
            end
        end
        
        % �������� �������� ��������� ��� ��������� ��������������
        function [ fvs, coef_map ] = transform_wmap(obj, rps, varargin)
        %TRANSFORM ������������ ������� ���������
        %   ��������� �������������� � �������������� ������� ������ RPTools
        %   ���������:
        %   .obj        
        %   rps         - ������� ����� ��
        %   [t_type]    - ��� ��������������
        %   [wave_name] - ����������� ��������
            if nargin>=3;   t_type = varargin{1};   % ��� ��������������
            else; t_type = obj.t_type; end 
            if nargin>=4; wave_name = varargin{2};  % ��� ��������
            else; wave_name = obj.wave_name; end
            switch t_type
                case 'none'
                    fvs =  rps;
                    coef_map = size(fvs,2);
                case 'afft' 
                    fvs = RPTools.afft(rps);
                    coef_map = size(fvs,2);
                case 'fft'
                    fvs = RPTools.fft(rps);
                    coef_map = size(fvs,2);
                case 'cwt'   
                    [fvs, coef_map] = RPTools.cwt(rps,wave_name);
                case 'acwt' 
                    [fvs, coef_map] = RPTools.acwt(rps,wave_name);
                case 'wt'   
                    [fvs, coef_map] = RPTools.wt(rps,wave_name);
                case 'pca'
                    [fvs] = RPTools.pca(rps);
                    coef_map = size(fvs,2);
            end
            
        end
        
        % ����������� ���������� ������
        function [wms_counts] = wms_hcounts(~, fvs, hist_bins)
            n_obj = length(fvs);
            wms_counts = cell(1, n_obj);
            hist_bins = gpuArray(hist_bins);
            for i_obj = 1:n_obj
                gpu_fvs = (fvs{i_obj});   
                gpu_fvs = gpu_fvs./sum(gpu_fvs,2);  % ���������� ��� ���������� ������� ����
                wms = asin(Slozhno.calc_w_metrics(gpu_fvs)/2)*2;  % ������ ������� �������
                clear gpu_fvs;
                wms_counts{1,i_obj} = (...
                        histcounts(wms, hist_bins )); % ���� ������
            end
            clear hist_bins
        end
            
        % ����� ���������
        function plot_metric_hist(obj, hist_bins, bms_counts, legend_entries)

            n_band = length(bms_counts);

            figure('Name',['Hist ' obj.scan_type],...
                'color', 'white','WindowStyle','docked'); 
            hold on;
            marker_type = { '*' 's' 'd' '+' 'x' 'o'};
            line_width = linspace(1, 2.8, n_band);
            line_color = linspace(0.2, 0.8, n_band);
            fig_h = zeros(1,n_band);
            %stem_offset = (hist_bins(2)-hist_bins(1)).*(1:n_band)./(n_band+2);
            stem_offset = linspace(-0.3,0.3, n_band)*(hist_bins(2)-hist_bins(1));
            max_count_val = 0;
            for i_band = 1:n_band
                fig_h(i_band) = stem(hist_bins(2:end)+stem_offset(i_band),...
                   bms_counts{i_band},['-' marker_type{i_band}],'filled',...
                   'LineWidth',line_width(i_band), 'Color', [1 1 1].*line_color(i_band));
                plot(hist_bins(2:end)+stem_offset(i_band),bms_counts{i_band}...
                   ,[':' marker_type{i_band} ''],'LineWidth',1)
                if max(bms_counts{i_band})>max_count_val, max_count_val =  max(bms_counts{i_band});end
            end
            %xlabel('\rho, rad','FontSize',10);ylabel('N(\rho)','FontSize',10);
            obj.setup_axis_labels('\rho','N(\rho)')
            xticks(hist_bins(2:end))
            xticklabels([num2str(hist_bins(1:end-1).','%3.2f-'), ...
                num2str(hist_bins(2:end).','%3.2f')])
            xtickangle(45)
            bar( hist_bins(2:end), zeros(1,length(hist_bins)-1) + max_count_val/200, 'k' )
            legend(fig_h,legend_entries)
            set(gca, 'XMinorTick', 'off', 'YMinorTick', 'on')
            xlim([hist_bins(1) hist_bins(end)])
            ylim([0 max_count_val])
        end
                
    end
    
    methods(Access = public, Static = true)
        
        % ������������ ��
        function [rps, meta] = generate_rps(icr_folder, gauss_hwidth)
        %GENERATE_RPS ������������ �� �� ��������� ���
            %   ���������:
            %   icr_folder -    - ����� � ���
            %   obj_id          - ��� ��� ��
            
            dim_icr = 200;      % ����� �������� � ���
            window_m = 3;       % ����� ��������� ���������, �.
            resol_ns = window_m * 2 / dim_icr / 3e8 * 1e9;  % �������� ��������� ���������� ���, ��
            % �������� ���
            filelist = dir([icr_folder '/*.icr']);        % �������� ��������� ������ ������ �� �������� �����
            for i = 1:length(filelist)      % ������� ������ ������ ���� ������ �� ���������
                filelist(i).name = horzcat(filelist(i).folder, '\', filelist(i).name); 
            end
            [icrs, meta] = RPTools.open_icr_file(...
                {filelist(1:end).name}.');  % �������� �� � ����������
            clear filelist
            % ������������ ��
            gauss_width_ns = gauss_hwidth * resol_ns;    % ������������ �������� � �c
            rps = RPTools.make_rps(icrs, gauss_hwidth);  % ������ ��
            disp(['������������ �� �� ����������: ' num2str(gauss_width_ns) '���']);
        end
            
        % ������������ ������� ���������
        function [rp_afft,meta] = generate_fvs(icr_folder, gauss_hwidth)
            %GENERATE_FVS ������������ ������� ��������� �� ��������� ���
            %   ���������:
            %   icr_folder -    - ����� � ���
            %   obj_id          - ��� ��� ��
            
            dim_icr = 200;      % ����� �������� � ���
            window_m = 3;       % ����� ��������� ���������, �.
            resol_ns = window_m * 2 / dim_icr / 3e8 * 1e9;  % �������� ��������� ���������� ���, ��
            % �������� ���
            filelist = dir([icr_folder '/*.icr']);        % �������� ��������� ������ ������ �� �������� �����
            for i = 1:length(filelist)      % ������� ������ ������ ���� ������ �� ���������
                filelist(i).name = horzcat(filelist(i).folder, '\', filelist(i).name); 
            end
            [icrs, meta] = RPTools.open_icr_file(...
                {filelist(1:end).name}.');  % �������� �� � ����������
            clear filelist
            % ������������ ��
            gauss_width_ns = gauss_hwidth * resol_ns;    % ������������ �������� � �c
            rps = RPTools.make_rps(icrs, gauss_hwidth);  % ������ ��
            disp(['������������ �� �� ����������: ' num2str(gauss_width_ns) '��']);
            %mesh(rps); axis vis3d;  % ������������ ��

            rp_afft = single(RPTools.afft(rps));                % AFFT
            rp_afft = RPTools.norm_rp(rp_afft, 'energy');
            disp('��������� ����������� ��������, ������������� �� �������')
            %mesh(rp_afft); axis vis3d;                 % ������������ ��������
%             ���������� ����� ��
%             fvdb_fname = ['fvdb_' obj_id '.csv'];       % ��� ����� ��
%             f_header = RPTools.get_feature_map...
%                 ('afft', size(rp_afft,2));              % ������� ��������� ��� �����
%             RPTools.save_csv...
%                 (rp_afft, meta, f_header, fvdb_fname);  % ���������� ���� �������� ���������
%             disp('���������� ���� �������� ���������')
        end
        
        % ������ ������ ��� ��������� �������� ���������
        function w_met = calc_w_metrics(feature_values)
        %FENORMSW ������������ ������� ������ ��������� �������� ���������
        %   ������������ ������� ����� ��������� ���������, 
        %   feature_values  - ������� ��������� �� ������� �������
            [n_aspect,nCounts] = size(feature_values); % ����� �� � �� �����������
            n_metric = nchoosek(n_aspect,2);  % ����� ������
            w_met=zeros(n_metric,nCounts...
                ,'like',feature_values);    % ���� ������ � �� ���������
            i_met = 0;
            disp(['�� ' num2str(n_metric) ' ������ ����������:     ']);
            for i=1:n_aspect
                for j=(i+1):n_aspect
                    i_met=i_met+1;
                    w_met(i_met,:) = feature_values(i,:) - feature_values(j,:); % �������� ���������
                end
                if i_met<1000000 ;fprintf('\b\b\b\b\b\b%6.0d',i_met); end
            end
            fprintf('\r');
            w_met = sqrt(sum(w_met.^2,2)); % ����������� ������������ �� �������
            % ����� ������������ �������� ����������� �������-�������
            % ���������� � �������
            
        end
                
        % ������ ������
        function [ rhos ] = calc_b_metrics( c1, c2, varargin )
        %fRhoCalc ������� ������� ������� ����� ��������� ���� �������
        %   ������� ��������� ��� ������ � m(C1Data) � n(C2Data) ���������.
        %   ������������ �������� ����� ����� ��������
        %   �� ������ mxn ������� ���������� ������ 
            n_obj = [size(c1,1) size(c2,1)];     % ����������� ����� �������� � �������
            v_dim = size(c1,2);
            if nargin == 3
                batch_dim = min([varargin{1}, v_dim]);
            else
                memory_info = memory; 
                memory_info = memory_info.MaxPossibleArrayBytes;
                memory_ratio = prod(n_obj)*v_dim*8/memory_info*(1.1);
                if memory_ratio>1
                    batch_dim = floor(v_dim/memory_ratio);
                else
                    batch_dim = v_dim;
                end
            end
            batch_edges = [1:batch_dim:v_dim v_dim+1];
            rhos = zeros(n_obj(1), n_obj(2), 'like', c1);
            n_batch = length(batch_edges)-1;
            for i_batch=1:n_batch
                rhos = rhos+sum( Slozhno.calc_b_vectors( ...
                    c1(:,batch_edges(i_batch):batch_edges(i_batch+1)-1),...
                    c2(:,batch_edges(i_batch):batch_edges(i_batch+1)-1)...
                    ).^2 , 3);
            end
            rhos = sqrt(rhos);
                
        end
        
        % ������ ������������� ������ �� ����������
        function [ c1 ] = calc_b_vectors( c1, c2 )
        %fRhoCalc ������� ������� ������� ����� ��������� ���� �������
        %   ������� ��������� ��� ������ � m(C1Data) � n(C2Data) ���������.
        %   [mem_limit_MB] 
        %   ������������ �������� ����� ����� ��������
        %   �� ������ (m x n x v_dim) ������� ��������� 
            n_obj = [size(c1,1) size(c2,1)];     % ����������� ����� �������� � �������
            v_dim = size(c1,2);
            c1 = repmat(c1.',1,1,n_obj(2));      % ���������� ������� -> v_dim x m x n
            c2 = permute(c2,[2 3 1]);            % ���������������� ������� -> v_dim x 1 x n
            for i_obj = 1:n_obj(1)
                c1(1:v_dim,i_obj,:) = (c1(1:v_dim,i_obj,:) - c2);   % ������ ������
            end
            c1 = permute(c1,[2 3 1]);           % ���������������� -> m x n x v_dim
        end
        
        % �������������� ������� ������� � ������ �����
        function [tr_mx] = met2rp(n_rp)
        %met2icr �������� ������� �������� �� ������� � i,j ��
        %   ���������� ������� ��������� <2 x n_wms>
        %   n_rp  - ����� ��
            n_wms = (n_rp^2-n_rp)/2; 
            tr_mx = zeros(n_wms, 2, 'uint16');
            i_met = 1;
            for i=1:n_rp
                for j=(i+1):n_rp
                    tr_mx(i_met,:)=[i j];
                    i_met= i_met+1;
                end
            end
        end
        
        % �������������� ������� ������� � ������ �����
        function [met_id] = rp2met(varargin)
        %met2icr ����������� ������ ������� ������ � ������� ��
        %   1 - ����� ��
        %   2 - 1� ������ ��
        %   [3 - 2� ������ ��]
        %   ������� �������� � ������, �������������� �� ���� ������� 
        %   ����������� ������� ������, ��������� �� �������.
            % ���� ������� ���-� ������ ��, �� ������ ���-� ���� �����. ������
            if (nargin==2)
                n_rp = uint32(varargin{1});     % �������������� � ����� ������
                i = uint32(varargin{2});        % �������������� � ����� ������
                if i>1 % ���� ������ ��, �� ������� �� �����
                    k_c = zeros(1,i-1,'uint32');    % ������ ������ �������
                    for j=1:i-1
                        k_c(j) = n_rp * j - sum(1:j)+i-n_rp;    % ������ �������� ������ � �������
                    end
                else
                    k_c = []; % ����
                end
                if i<n_rp
                    k_r = (i*n_rp-n_rp+i+1-sum(1:i)):1:(n_rp*i-sum(1:i));
                else 
                    k_r = [];
                end
                met_id = [k_c k_r];
            % ���� ������� ���-� ���� ��, �� ������ ���-� ������� ���. ����
            elseif nargin ==3
                n_rp = varargin{1};
                i = max([varargin{2} varargin{3}]); 
                j = min([varargin{2} varargin{3}]); 
                met_id = n_rp * j - sum(1:j)+i-n_rp; % ������ ���-�� �������
            % ����� ����� ���������� �� ������������� ���������
            else 
                error('Argument input error')
            end
            
        end
        
        % ������������ ���������� ������
        function [metrics_hist] = reduce_ms(wms,thr_ms)
        %REDUCE_MS - ������� ����� ������, ����������� ������ � �������
        %thr_ms
        %   wms          - ������ � ���������
        %   thr_ms       - ������ � ���������� ����������
        %   metrics_hist - ������ � ������ ������, ����������� ��������� ��������
            wms = sort(wms);            % ����������� ������� �� �����������
            n_wms = length(wms);        % ����� ������
            n_thr_ms = length(thr_ms);  % ����� ��������� ��������
            metrics_hist = zeros(size(thr_ms,1),size(thr_ms,2)); 
            for i = 1:n_thr_ms
                ind = find(wms>thr_ms(i),1);            % ����� ������ ������ ������� ������� ������ = ����� ������ ������ ������
                if isempty(ind); metrics_hist(i) = 0;   % ���� ���, �� 0 ������ ������ ������
                else; metrics_hist(i) = n_wms-ind+1;    % ����� ����� ������ ������ ������ ����� ������ ����� ����� ����� ������ ������� ������� ������
                end
            end
        end
        
        % ������������ ������������� ���� ������
        function [en_rp] = reduce_rpdb(wms,thr_ms)
        %REDUCE_RPDB �������� �� �� �������� ������
            wms = single(wms);
            thr_ms = single(thr_ms);

            n_thr = length(thr_ms);     % ����� �������/��������
            n_wms = length(wms);        % ����� ������
            n_rp = (1+sqrt(1+8*n_wms))/2; % ����� ��

            [~,tr_s_wms] = sort(wms);   % ������� ��������� � �������������� ������
            tr_s_wms = uint32(tr_s_wms);
            tr_wms2rp = Slozhno.met2rp(n_rp);       % ������� �������� �� �������� ������� ������ � �������� ��
            en_rp = ones(n_thr, n_rp,'logical');    % ������� ������ ������� �� (�������)
            en_wms = ones(1,n_wms,'logical');       % ������� ������� ������
            i_thr = 1;                  % ������� ������������ ��������� ������
            i_left_rp = n_rp;
            for i_s_wms=1:n_wms         % ���� �� ������� ������������� �������� ������� ������
                i_wms = tr_s_wms(i_s_wms);          % ��������� ������� �������
                if en_wms(i_wms)                    % �� ��������� �� ������� �����
                    while wms(i_wms)>thr_ms(i_thr)  % �� ����� �� ������� �� ���. ��������
                        i_thr=i_thr+1;              % ���� ��, ������� � ����. ���������
                        if i_thr>n_thr; return; end % ����� �� ������������ ���� ��������� �����
                    end
                    i_rp = tr_wms2rp(i_wms,:);  % ��������� ���� �� ����� ���. ����. �������
                    if i_left_rp==2             % ���� ��������� �������
                        en_rp(i_thr:end, i_rp(2))=false;    % ������� ���� � ����� ��
                        break;
                    end
                    min_met = zeros(1,2);       % ����� ��� ���. �������
                    i_ex_wms = cell(2,1);       % ����� ��� ������� ������ ������
                    for i = 1:2                 
                        tmp = Slozhno.rp2met(n_rp, i_rp(i));% ��������� �������� ���� ������ ����� ��������
                        i_ex_wms{i} = tmp(en_wms(tmp));     % ����-�� ���-�� ����� ����-� ������
                        tmp = sort(wms(i_ex_wms{i}));
                        min_met(i) = tmp(2);    % ��������� ����������� �������
                    end
                    if min_met(1)<min_met(2)    % ���� ���. ������� 1��� �� ������
                        en_wms(i_ex_wms{1})=false;          % ������� ���� � ����� ������
                        en_rp(i_thr:end, i_rp(1))=false;    % ������� ���� � ����� ��
                    else
                        en_wms(i_ex_wms{2})=false;          % ������� ���� � ����� ������
                        en_rp(i_thr:end, i_rp(2))=false;    % ������� ���� � ����� ��
                    end
                    i_left_rp = i_left_rp-1;                % ��������� ����� ���. ��
                end
            end 
        end
        
        % ������ ����������� n-��� �������
        function signal = bworth_filter(signal,t_step_ns,filt_band_mhz, varargin)
        %filter_2order - ������� ���������� ������� �������� � ��� �����������
        %���������:   
        %   signal          - ������� ������� ������ ����� ������ ����� ��������
        %   t_step_ns       - ��� �������������, ��
        %   filt_band_mhz   - ��������� ������� �� -3��, ���
        %   [n_filter]      - ������� �������
            if nargin == 4,  n_filter = varargin{1};    % ������ �������, ��������� �� �����
            else,            n_filter = 2;              % �� ��������� ������ 2-�� �������
            end
            f_sampling_mhz = 1e3/t_step_ns;             % ������� �������������
            w_n = filt_band_mhz/(f_sampling_mhz/2);                 % ������������� ������� �����
            [b,a] = butter(n_filter, w_n, 'low');  % ������������ ������������� ������� �����������
            signal = filter(b,a,signal.').';         % ���������� �������
        end
        
        % ����������� ��������� ���� ��������
        function setup_axis_labels( x_label,y_label )
        %SETUP_AXIS_LABEL( xLabelH,yLabelH ) ��������� ����������������� ������� ����
        % ������� �������
        %   xLabelH - ������������� ������� ��� x - ���������� �� ���. xlabel(_)
        %   yLabelH - ������������� ������� ��� y - ���������� �� ���. ylabel(_)
        %   ��������� ��������� ������� ��� x � ������ ������ ����, 
        %   ������� ��� y - � ������� ����� ���� � ������������ � �����. ���������
        
        x_label_h = xlabel(x_label,'FontSize', 12);
        y_label_h = ylabel(y_label,'FontSize', 12);
        
        set(x_label_h,'Units','normalized')
        set(x_label_h...
            ,'Position', [1 0 0]...
            ,'VerticalAlignment','top'...
            ,'HorizontalAlignment', 'center'...
            );
        
        set(y_label_h,'Units','normalized')
        set(y_label_h...
            ,'Position', [0 1.01 0]...
            ,'Rotation',0 ...
            ,'VerticalAlignment','bottom'...
            ,'HorizontalAlignment', 'right'...
            );
        end

        % ������ ����������� �������
        function cur_rho = buryi(fvs_red)
            n_obj = length(fvs_red);
            cur_rho = ones(n_obj-1, n_obj);   % ��������� ������
            for i_obj = 1:n_obj-1
                for j_obj = (i_obj+1):n_obj
                    cur_rho(i_obj,j_obj) = min(Slozhno.calc_b_metrics(...
                        fvs_red{i_obj},fvs_red{j_obj}),[],'all');  % ������ ����������� ������� ��� i-j ���������
                end
            end
            cur_rho = min(cur_rho,[],'all');   % �����
        end
        
        % ������ �� n-� ����������� �������
        function cur_rho = nmin_metric(fvs_red, n_min_ratio)
        %nmin_metric ��������� ������ ������������ ������������ ��������� �� �������� �����������
        %   fvs_red     - ������ ����� � ��������� ��������� 
        %   n_min_ratio - ���� ������, �� ������� ��������� �������� ������� (<1)
        %   w_rep       - ������������ ����������
        %   TODO: ������������� ����� ������������ ��������� ��������
            n_obj = length(fvs_red);
            n_asps = cellfun(@(fvs) size(fvs,1), fvs_red);          % ����� �������� � ������
            mat_ij = repmat((1:n_obj),n_obj,1);                 % 
            asp_i = mat_ij(triu(ones(n_obj,'logical'),1).');    % ������ ������� ���������
            mat_ij = mat_ij.';                                  %
            asp_j = mat_ij(triu(ones(n_obj,'logical'),1).');    % ������ ������� ���������
            n_met_total = sum(n_asps(asp_i).*n_asps(asp_j));     % ����� ��
            % ���������-�����������
            elem_rho = prod(n_asps, 'all');
            m_ij = arrayfun(@(i,j) elem_rho/(n_asps(i)*n_asps(j)), ...
                asp_i, asp_j);              % ��������� ����������
            m_ij = floor(m_ij/min(m_ij));
            % ��� ����� ������ * ���� ���������� ������ * ����� �� ������� ����� �������������
            cur_rho = cell(length(m_ij),1);
            n_met_thresh = ceil((n_met_total * n_min_ratio)./m_ij);  % 
            for i_case = 1:length(m_ij)
                cur_rho{i_case} = sort( reshape(...
                    Slozhno.calc_b_metrics(fvs_red{asp_i(i_case)}, fvs_red{asp_j(i_case)}), ...
                    1,[]), 'ascend');  % ������ ������ ��� �����., ������� � ������ � ����������

                cur_rho{i_case} = repmat( ...
                    cur_rho{i_case}(1:n_met_thresh(i_case))... % �������� �� ������� n/ ����. ����������
                    , 1, m_ij(i_case));          % ���������� (�������� nthr/���� ����������)
                                 
            end
            cur_rho = sort(horzcat(cur_rho{:}),'ascend');   % ���������� ����� �������
            cur_rho = cur_rho(ceil(length(cur_rho)*n_min_ratio)); 
        end
        
        % ������ �� n-� ����������� ������� ��� ���������������� �������
        function cur_rho = nmin_metric_balanced(fvs_red, n_min_ratio)
        %nmin_metric ��������� ������ ������������ ������������ ��������� �� �������� �����������
        %   fvs_red     - ������ ����� � ��������� ��������� 
        %   n_min_ratio - ���� ������, �� ������� ��������� �������� ������� (<1)
        %   w_rep       - ������������ ����������
            n_obj = length(fvs_red);
            n_asps = cellfun(@(fvs) size(fvs,1), fvs_red);      % ����� �������� � ������
            mat_ij = repmat((1:n_obj),n_obj,1);                 % 
            asp_i = mat_ij(triu(ones(n_obj,'logical'),1).');    % ������ ������� ���������
            mat_ij = mat_ij.';                                  %
            asp_j = mat_ij(triu(ones(n_obj,'logical'),1).');    % ������ ������� ���������
            n_met_total = sum(n_asps(asp_i).*n_asps(asp_j));    % ����� ������
            % ���������-�����������
            elem_rho = prod(n_asps, 'all');
            m_ij = arrayfun(@(i,j) elem_rho/(n_asps(i)*n_asps(j)), ...
                asp_i, asp_j);              % ��������� ����������
            m_ij = floor(m_ij/min(m_ij));
            % ��� ����� ������ * ���� ���������� ������ * ����� �� ������� ����� �������������
            cur_rho = cell(length(m_ij),1);
            n_met_thresh = ceil((n_met_total * n_min_ratio)./m_ij);  % 
            for i_case = 1:length(m_ij)
                cur_rho{i_case} = sort( reshape(...
                    Slozhno.calc_b_metrics(fvs_red{asp_i(i_case)}, fvs_red{asp_j(i_case)}), ...
                    1,[]), 'ascend');  % ������ ������ ��� �����., ������� � ������ � ����������

                cur_rho{i_case} = repmat( ...
                    cur_rho{i_case}(1:n_met_thresh(i_case))... % �������� �� ������� n/ ����. ����������
                    , 1, m_ij(i_case));          % ���������� (�������� nthr/���� ����������)
                                 
            end
            cur_rho = sort(horzcat(cur_rho{:}),'ascend');   % ���������� ����� �������
            cur_rho = cur_rho(ceil(length(cur_rho)*n_min_ratio)); 
        end
        
        % ������ ������������ ��������� �� ������ ����������
        function cur_rho = hist_acc(fvs_red, hist_bins, n_metrics_fraction)
        % ������ ����������� ������
                    hist_counts = Slozhno.cpu_linear_histcounts(fvs_red,hist_bins);         % �����������
                    hist_counts = vertcat(hist_counts{triu(ones(length(fvs_red),'logical'),1)});  % �������� �������� �� ������� ����������� ������� ��� ���������
                    hist_counts = sum(hist_counts,1);                   % ��������������
                    acc_hist = cumsum(hist_counts);                     % �������� �� ����
                    thrld = acc_hist(end)*n_metrics_fraction;           % ����� �� �����������
                    cur_rho = hist_bins(find(acc_hist>thrld,1));        % �������� �������, ��������������� ������ ��� ���. ������������ ���������
        end
        
        % �������� ������������ ���������
        function [fspace_rhos, fspace_map, varargout] = fs_reduction(fvs, varargin)
        %fs_reduction �������� � ������ ����������� �������
        % �� ������
        % fspace_rhos - ������-������� � ���������
        % fspace_maps - ������� �� ������� ������� ����������� ��������� ������������
        % ���������: 
        % fvs - ������� ��������� (�� �������)
        % ����������� ���������:
        % dimensions, reduction_type, n_metrics_fraction, hist_edges,
        % obj_weight
        warning('����� ������� �������. ������������ ������������� ������������ ����� RecursiveReduction.')
            dim_fv = size(fvs{1},2);            % ����������� ��
            n_obj = length(fvs);
            kwargs = KeywordArguments(...
                'dimensions',1:dim_fv, ...          % ����������� �����������
                'reduction_type', 'none',...        % ����������� �����������
                'n_metrics_fraction', 0.05,...      % ����� ������
                'hist_edges', linspace(0,2,51),...  % ������� ���������� �����������
                'obj_weights',ones(n_obj,1)/n_obj,... % ���� �������� (��� ����� ��������� ���-���)
                'n_nearest',5, ...          % ����� ����������� �������
                'm_technique','mean',...    % ����� ������
                'k_alien',2 ...               % ����� �����
                );
            [dimensions, reduction_type, n_metrics_fraction, hist_edges, obj_weights, ...
                n_nearest, m_technique, k_alien] =  ...
                kwargs.parse_input_cell(varargin);
            % ��������� ������������� ��� ������
            switch reduction_type
                case 'nmin'
                    if range(arrayfun(@(i) size(fvs{i},1),1:n_obj))==0  % ���� ��� ����������, �� ������� ������� ����������������
                        rho_estimate = @(fvs) Slozhno.nmin_metric_balanced(fvs, n_metrics_fraction);
                    else
                        rho_estimate = @(fvs) Slozhno.nmin_metric(fvs, n_metrics_fraction);
                    end
                case 'buryi'
                    rho_estimate = @(fvs) Slozhno.buryi(fvs);
                case 'mhist'
                    rho_estimate = @(fvs) Slozhno.hist_acc(fvs, hist_edges, n_metrics_fraction);
                case 'minalien'
                    rho_estimate = @(fvs) Slozhno.minalien(fvs, n_nearest, m_technique, k_alien);
                otherwise
                    error('����������� ����� ����� ��������. ��������������: nmin, buryi, mhist, minalien')
            end
            dimensions = dimensions(dimensions<=dim_fv);
            n_spaces = length(dimensions);
            
            fspace_rhos = zeros(n_spaces,1);        % �������� ��������� ������� ��� ������������
            fspace_map = zeros(n_spaces, dim_fv);   % ������� ����������� ����������� - �� ������� �� �����������
            fprintf('\t��������� �����������: %d.',2^dim_fv);
            fprintf('\t���������: 00.00%%');
            for feat_val = 1:2^dim_fv-1             % �� ���� �������� ������������� ���������
                cur_feat_map = de2bi(feat_val,dim_fv,'left-msb'); % ��������� ����� �������� ������������
                i_fmap = find(sum(cur_feat_map)==dimensions);  % ������ ��������������� �����������
                if i_fmap
                    cur_feat_map = logical(cur_feat_map);
                    fvs_red = cellfun(@(fvs) fvs(:,cur_feat_map), fvs, 'UniformOutput',false); % ��������
                    
                    [cur_rho] = rho_estimate(fvs_red);
                    
                    if cur_rho>=fspace_rhos(i_fmap)     % ���� �������� ������� ������ �������� ��� ��. ��������� ����� �����������
                        fspace_rhos(i_fmap) = cur_rho;
                        fspace_map(i_fmap,:) = cur_feat_map(:);
                    end
                    fprintf('\b\b\b\b\b\b%05.2f%%',feat_val/2^dim_fv*99.99);
                end
            end
            fprintf('\b\b\b\b\b\b\b 100%%\n');
            if nargout==3, varargout{1} = dimensions; end
                        
        end
        
        % ����������� ������
        function [bms_counts] = cpu_linear_histcounts(fvs, hist_bins) 
        %cpu_histcounts
        % 
            n_obj = length(fvs);
            bms_mod = cell(n_obj, n_obj);       % ������� ��� �������
            bms_counts = cell(n_obj, n_obj);    % ������� �������� ������� ������
            for i_obj = 1:(n_obj-1)             % ������ ������
                for j_obj = (i_obj+1):n_obj
                    bms_mod{i_obj,j_obj} = Slozhno.calc_b_metrics(fvs{i_obj},fvs{j_obj}); % ������ ������� ������
                    bms_counts{i_obj,j_obj} = ...
                        histcounts(bms_mod{i_obj,j_obj}, hist_bins );
                end
            end
        end
        
        % ������� ������� �������� ������������ ��������� tsne
        function q_sep = tsne_separability(data, varargin)
        % ��������� ���������� - ������ tsne � ������� �������    
            kwargs = KeywordArguments(...
                'n_nearest',false, ...          % ����� ����������� �������
                'm_technique','mean',...    % ����� ������
                'k_alien',3 ...               % ����� �����
                );
            [n_nearest, m_technique, k_own] =  ...
                kwargs.parse_input_cell(varargin);
            if ischar(data)
                in_data = load(data,'tsne_data');   % ��������� ����������
                features = vertcat(in_data.tsne_data.tsne_scatter);
                targets = {in_data.tsne_data.name};
                if n_nearest
                else, n_nearest = in_data.perplex;  % ��������� ����� ����. �����. ���-�� tsne
                end
            elseif(iscell(data))
                if nargin==3
                    features = vertcat(data{:});
                    targets = varargin{1};
                else
                    error('Missing argument');
                end
            else
                error('������ � ������� ����������')
            end
            q_sep = Slozhno.nn_separability(features, targets, ...
                'n_nearest',n_nearest, 'm_technique',m_technique, 'k_alien', k_own);
        end
        
        % ������� ������� ������������ �� �������� ���������
        function q_sep = nn_separability(features, targets_str, varargin)
        % ��������� ���������� - ������ ������� ��������� � ������� �������    
            kwargs = KeywordArguments(...
                'n_nearest',5, ...          % ����� ����������� �������
                'm_technique','mean',...    % ����� ������
                'k_alien',2 ...             % ����� �����
                );
            [n_nearest, m_technique, k_alien] =  ...
                kwargs.parse_input_cell(varargin);
            n_samples = size(features, 1);
            rhos = Slozhno.calc_b_metrics(features, features);  % ��������� �������
            [~, sort_indexes] = mink(rhos, n_nearest+1, 2);                  % �����������
            clear('rhos')
            [~, ~, targets] = unique(targets_str);        % ������������� ����� � id-������
            
            switch m_technique
                case {'mean' 'thresh'}
                    neq_matrix = zeros(n_samples, n_nearest); % ������������ �-�� �����-���
                    for i=1:n_nearest
                        neq_matrix(:, i) = targets ~= targets(sort_indexes(:, i+1));
                    end
                    switch m_technique
                        case 'mean'     % ������� ����� ����� ������� ������ ����� �������
                            q_sep = mean(sum(neq_matrix, 2)/n_nearest);
                        case 'thresh'   % ���� �������� ��� ������� ����� ��������� ������� ������ ������
                            q_sep = sum((sum(neq_matrix, 2)>(k_alien)))/n_samples;
                    end
                    
                case 'prat'
                    local_groups = targets(sort_indexes(:,1:n_nearest));   % � ������� -> ������
                    group_modes = mode(local_groups,2);              % ������������� ����� � ������
                    q_sep = 1-mean(sum(group_modes==local_groups,2))/n_nearest; % ������� ����� ������.
                    
            end
            
        end
        
        % ������ ��������������� ����� �����
        function cur_rho = minalien(fvs, n_nearest, m_technique, k_alien)            
            n_obj = length(fvs);
            n_fvs = cellfun(@length, fvs);
            targets = zeros(sum(n_fvs),1);
            idx=1;
            for i_obj=1:n_obj
                targets(idx:idx+n_fvs(i_obj)-1)=i_obj;
                idx = idx+n_fvs(i_obj);
            end

            fvs = vertcat(fvs{:});
            
            n_samples = size(fvs, 1);
            rhos = Slozhno.calc_b_metrics(fvs, fvs);  % ��������� �������
            [~, sort_indexes] = sort(rhos, 2);                  % �����������
            clear('rhos')
            neq_matrix = zeros(n_samples, n_nearest); % ������������ �-�� �����-���
%             [~, ~, targets] = unique(targets_str);        % ������������� ����� � id-������
            for i=1:n_nearest
                 neq_matrix(:, i) = targets ~= targets(sort_indexes(:, i+1)); % i+1 -���� ����, ��� 1� ����� - ��� �����
            end
            switch m_technique
                case 'mean'     % ������� ����� ����� ������� ������ ����� �������
                    q_sep = mean(sum(neq_matrix, 2)/n_nearest);
                case 'thresh'   % ���� �������� ��� ������� ����� ��������� ������� ������ ������
                    q_sep = sum((sum(neq_matrix, 2)>(k_alien)))/n_samples;
            end
            cur_rho = 1 - q_sep;
        end
        
    end
        
end

