classdef RPTools
    %RPTools Ver.01.02 ����� �������� ������� ���  ���������� ������������
    %   ��� ������� ��������� ����������
    
    % �������� �������
    methods(Access = public, Static = true)
        
        % ������������ ��
        function [rps] = make_rps(icrs, gauss_hwidth)
            % ����������� ������� �� ���� ��� � ������� gauss_width ��
            % ����������
            hw2sig = 2.35;              % ��������� ������ �� �������� � �����
            sig = gauss_hwidth/hw2sig;  % �������� � ����� �����
            t_axis = linspace(-4,4, sig*8);
            imp = exp(-t_axis.^2/2);    % ������� ��
            rps = conv2(icrs, imp);     % ������� ��� � ��������� ��������
            if mod(size(rps,2),2)
                rps = rps(:,1:end-1);   % ���� �������� ����� ��������, �� �������� ���������
            end
        end
        
        % ���������� ���
        function [rp_samples, meta] = add_noise(rps, meta, sigma, n_real, multiplier)
        % ������������ ��������� ���������� � ���������� ����������� �����
        %   ����: 
        %   rps � �������, ���������� ������� ��
        %   meta � cell array � �������������� �������� � ������ �������� ? � ?
        %   sigma � �������� ��� ��� ���
        %   n_real � ����� ����������
        %   �����:
        %   rps � �������, ���������� ��������� ���������� ��
        %   meta � cell array, ���������� ���������� ��� ������������ ������� ���������� ��
            signal_power = mean(rps.^2,2);  % �������� �� �� ���������, �������. 3� ������
            [n_rps, dim_rp] = size(rps);
            total_samples = n_rps*n_real;   % ����� ����������
            fprintf('\t���. �/�: ')
            if length(sigma)==1             % ���� ������ ���� �������� �/�
                signal_power = num2cell(signal_power/sigma^2);  % �/� = ��.����.��(3�����)/��.����.����(3�����)
            else                            % ���� ����� ����� ��������� �/�
                %TODO: �������� ����� ������ ������ �/�
                signal_power = num2cell(signal_power/max(sigma)^2); % ������� ������ �/�
            end
            if ~isempty(meta)
                [meta(:).snr] = deal(signal_power{:});
                [~,i_min] = min([meta.snr]); [~,i_max] = max([meta.snr]);
                fprintf('�� %01.1f(%s(%+03d,%+03d)\x00B0) �� %01.1f(%s(%+03d,%+03d)\x00B0)\n',...
                    meta(i_min).snr,meta(i_min).name, meta(i_min).asp_a, meta(i_min).asp_b,...
                    meta(i_max).snr,meta(i_max).name, meta(i_max).asp_a, meta(i_max).asp_b);
            else, warning('No metadata')
            end            
            % ��������� ������
            max_shift = dim_rp*(multiplier-1);      % ������������ ��������� �����
            if(multiplier>1)
                fprintf('\t���� ���� ������: %d ���� ���\n', multiplier);
                r_start  = 1 + ...
                    (floor(max_shift * rand(total_samples,1))); % ������������ ������� �� ���������� ��������� ������� 0.1*ones(total_samples,1))); 
            else
                r_start = ones(total_samples,1);    % ���� ���������������
            end
            r_stop   = r_start + dim_rp - 1;        % �������� �������
            sigma_mix = (max(sigma)-min(sigma))*rand(total_samples,1,'like', rps)+...
                min(sigma); % ��������� ���������� ������������ �� ����� ��� �� ����� ����
                        
            rp_samples = randn( total_samples, dim_rp*multiplier, 'like', rps );    % ������� � �������
            rp_samples = rp_samples .* sigma_mix;
            % ���������� � ������� �������� �������� ��������� �������
            for i = 1 : total_samples               % ���� �� ���� �����������
                rp_samples(i, r_start(i):r_stop(i)) = ...
                        rp_samples( i , r_start(i):r_stop(i) ) + ...
                        rps(mod(i-1, n_rps)+1,:);
            end
            meta = repmat(meta,n_real,1);
            
        end
        
        % �������������
        function rps = decimate(rps, frac)
            rps = movmean(rps, frac, 2);
            rps = rps(:, 1:frac:end);
        end
                
        % ����������� ���
        function dig_rp = adc(rp,varargin)
        %adc ��������� �������������� ������������ ������ � ���������
            if nargin == 1      % ���� �������� �� ��������, �� ����������� �� ���������
                dtype = 'uint8';    % �� ��������� 8 ���
                nMax = 2^8;         % �� ��������� 8 ���
            elseif nargin == 2  % ���� �������� ��������
                bit = floor(varargin{1});
                nMax = 2^bit;   % ������������ ���������� ��������
                if      and(bit >3, bit<=8);    dtype = 'uint8' ;   % 8 ���
                elseif  and(bit > 8, bit<=16);  dtype = 'uint16';   % 16 ���
                %elseif  and(bit > 16, bit<=32); dtype = 'uint32';   % 32 ����
                else; error('�������� ��� ������ ���� � �������� �� 4 �� 16');
                end
            else 
                error('�������� ����� ������� ����������')
            end
            rp_peak_value = max(rp,[],'all'); 
            if rp_peak_value>1
                warning('������������ �������� ���: �������.'); 
                rp(rp > rp_peak_value) = 1; % �������� ���������
            end
            dig_rp = cast(floor(rp * nMax),dtype);
        end
        
        % ����������
        function rps = norm_rp( rps, norm_mode )
            [n_rps, ~] = size(rps);
            switch norm_mode
                case 'mean'     % ���������� �� �������� ��������
                    norm_factors = 1./mean(rps,2);
                    norm_factors = norm_factors/max(norm_factors);
                case 'max_peak' % ���������� �� ���� ���������� ��
                    norm_factors = 1./max(rps,[],2)*.7; % ��� ���������� ������� ����������� ������ 0.7
                case 'energy'
                    norm_factors = 1./sqrt((sum(rps.^2, 2))); % ������������ ������� � ������������� �� N                    
                case 'mod'
                    norm_factors = 1./sqrt((sum(rps.^2, 2))); % ������������ ������� � ������������� �� N                    
                otherwise
                    if ~strcmp(norm_mode,'none') % ������� �������������� ���� ���� �� �������
                        warning('���������� �� �����������')
                    end
                    return
            end
            norm_factors(norm_factors==inf) = 0;
            for i = 1:n_rps
                rps(i,:) = rps(i,:)*norm_factors(i);
            end
                        
        end
        
        % ������� ���c�� ������ �������
        function ind = find_center(randrp)
        %FINDCENTER     ����� ������ ���� � ������� �� �������
            [nrp, nDim] = size(randrp);
            halfMass = sum(randrp,2)/2;         % ����� �������� 
            ind = zeros(nrp,1,'uint32');        % ������� ������� ��������
            sbuf = zeros(nrp,1);                % ������� ������� �������� ����
            for i = 1:nDim 
                sbuf = sbuf + cast(randrp(:,i),'like',sbuf);      % ��������� ������������� ��������
                mask = sbuf <= halfMass;
                ind = ind + cast(mask,'uint32');% ��������� ������
            end
        end
        
        % ��������� ��������� ����������������
        function rps = rp_time_snap( rand_rps, snap_mode, nc_offset )
        %RP_TIME_SNAP
        %   ��������� ��������� ���������������
        %   rand_rps    - ������� �� �� �������
        %   snap_mode   - ������ ���������� ��������� ���������������: 
        %       energy      - ����� �������
        %       max_peak    - �� ����
        %   nc_offset   - ������ ��������� ��. ���������� ����� ������� ��
        %   � ������ ������� �� ����� �������� ����� ����.
            [n_rps, ~] = size(rand_rps);
            % ����� ����������� ����� �� ��� ��������
            switch snap_mode
                case 'energy'       % ����������� ����� - ����� �������
                    snap_index = RPTools.find_center(rand_rps);
                case 'max_peak'     % ����������� ����� - ���������� ��������
                    [~, snap_index] = max(rand_rps,[],2);
                case 'none'         % �� ������������ ��� ����
                    rps = rand_rps;
                    return
                otherwise
                    error(['����� ' snap_mode ' �� ��������������']);
            end
            
            rps = zeros(n_rps, 2*nc_offset, 'like',rand_rps);   % ����������� ������� ��� ��
            rand_rps = [...
                zeros(n_rps,nc_offset,'like',rand_rps), ...
                rand_rps, ...
                zeros(n_rps,nc_offset,'like',rand_rps)];        % ��������� ������� �� ������ ������ �� �������
            start = snap_index + 1;         %
            stop = snap_index + 2*nc_offset;
            for i = 1:n_rps
                rps(i,:) = rand_rps( i, start(i):stop(i) );
            end
            
        end

        % �������� ��������� ��������� ���
        function [ features, coef_index ] = cwt_reduce( features_src, coef_map ,n_features)
        %CWT_REDUCE
        % ����������� ����� ���������, ��� ��� ����� ���������� �� ������
        %TODO: ��������� ������������� �� ����� ������ ��������� �������,
        %����������� ��� ���
            [n_rps, ~] = size(features_src);
            features = zeros(n_rps, n_features);
            coef_index = zeros(2 ,n_features); % ������������ ������� � ����� ������������
            
            n_coef = coef_map(end); % ���������� ����� �������������
            cntr = 100;
            i_lev = length(coef_map)+1;
            for i = 1:n_features
                if cntr > n_coef
                    i_lev = i_lev - 1;          % ��������� �������
                    n_coef = coef_map(i_lev);   % ����������� ������ ����������
                    cntr = 1;                   % ������� ���������
                    lev_premult = + 1 - 2 * mod( n_coef, 2 ); % ��������� ��������� �� �������� ����� �������������
                    % odd = + ; even = - 
                end
                offset = floor(cntr / 2);
                sign_premult = 1 - 2 * mod( cntr, 2 ); % ���������� � + 0

                i_coef =  ceil(n_coef/2) +...
                   lev_premult * sign_premult * offset;  % ��������� �� ������ ���� acwt (5/2 = 2.5 -> 3 -> 2 -> 4...)
                ind = sum( coef_map(1:i_lev) ) - coef_map(i_lev)  + i_coef;
                features(:,i) = features_src(:,ind);
                coef_index(:,i) = [i_lev; i_coef];
                cntr = cntr + 1;     % ��������� �������
            end
        end
        
        % ������������ ���������� �� �� ��������� ������� � ���������� �����
        function samples = generate_samples(rp,nSamples,nExtraDim,dSigm)
        %FRANDOMSHIFT ������������ ���������� ������ ��
        %   ������������ �� ���������� �� �� ������ ���������� ����������
        %   ������� �������
        %   Args: 
        %       rp - �������, � ������� �� ����������� �� �������.
        %       nSamples - ����� ����� ���������� �� ������ �������� ��.
        %       nExtraDim - �� ������� ��� ������������� �������� ����������.
        %       dSigm - ��� ���.
            [nIcrs, dimrp] = size(rp);      % ����� �� � �� �����������
            totalSamples = nIcrs*nSamples;  % ����� ����������
            
            samples = dSigm * ... 
                randn( totalSamples, dimrp*nExtraDim, 'like',rp );    % ������� � �������
            maxTShift = dimrp*(nExtraDim-1);% ������������ ��������� �����
            rStart  = 1 + (floor(maxTShift * rand(totalSamples,1)));% ������������ ������� �� ���������� ��������� �������
            rStop   = rStart + dimrp - 1;   % �������� �������
            i = 0;
            for iIcr = 1 : nIcrs            % ���� �� ���� �������� (���)
                for iSample = 1 : nSamples  % ���� �� ���� �����������
                    i = i + 1;
                    samples( i , rStart(i):rStop(i) ) =...
                        samples( i , rStart(i):rStop(i) ) + rp(iIcr,:);
                end
            end
        end
    
    % ������������ 
        
        % ����������� ������� ������������
        function signal = interpol(signal, fs, fq, interp_method)
        %interol �������� ������� ������������� � �������������� ������������
        %   ������� ������������ ������� �� �������� ������� ������������ fs 
        %   � ������� ������� fq �� ��������� ������ ������������
        %   Args:
        %       signal - �������, � ������� �� ����������� �� �������.
        %       fs - �������� ������� ������������� �������.
        %       fq - ������� ������� ������������� �������.
        %       interp_method - ����� ������������. �������������� ������:
        %           none     �� ��������� ��������� �������;
        %           freq     ��������� ������������;
        %           poly     �������������� ������������;
        %           Lagrange ���������� ������������;
        %           spline   ������������ ���������.
            switch interp_method
                case 'none'        
                case 'freq',        signal = RPTools.freq_interpol(signal, fs, fq);     % ��������� ������������
                case 'poly',        signal = RPTools.poly_interpol(signal, fs, fq);     % �������������� ������������
                case 'Lagrange',    signal = RPTools.lgrnge_interpol(signal, fs, fq);   % ���������� ������������
                case 'spline',      signal = RPTools.spline_interpol(signal, fs, fq);   % ������������ ���������
                otherwise, warning('����������� ����� ������������! ������������ �� ���������')
            end
        end
        
        % ������������ ����������
        function signal = poly_interpol(signal, fs, fq)
            t_axis = 0:1/fs:1/fs*(size(signal,2)-1);
            t_axis_est = 0:1/fq:t_axis(end);
            signal = interp1(t_axis,signal.',t_axis_est, 'pchip').'; % ������������ �����. 3� ����
            
        end
        
        % ������������ � �������������� ���������� ������
        function signal = freq_interpol(signal, fs, fq, varargin)
        % 
            if nargin>3, err = varargin{1}; else, err = 0.01; end
            [p, q]=rat(fq/fs, err);                 % ��������� ���� � ���� ��� �����  
            signal = resample(double(signal.'), p,q).';
        end
        
        % ������������ ���������
        function signal = spline_interpol(signal, fs, fq)

            t_axis = 0:1/fs:1/fs*(size(signal,2)-1);
            t_axis_est = 0:1/fq:t_axis(end);
            signal = spline(t_axis, signal, t_axis_est);
        end
        
        % ���������� ������������
        function signal = lgrnge_interpol(signal, fs, fq, varargin)
            if nargin>3,    n_lagrange = varargin{1};   % ������� ���������
            else,           n_lagrange = 5;             % ������� ���������
            end

            [int_factor, q]=rat(fq/fs, 0.01);         % ��������� ���� � ���� ��� �����  
            lgr_filter = intfilt(int_factor,n_lagrange,'Lagrange');   % ������ �/������������

            signal = upsample(signal.',int_factor);            % ���������� �����
            signal = filter(lgr_filter,1,signal);       % ���������� �������
            signal(1:floor(mean(grpdelay(lgr_filter))),:) = []; % ������� ������� ��������
            signal = signal(1:q:end,:).';
        end
        
    % ������ � �������
        
        % ��������� csv
        function [] = save_csv( features, meta, f_header, filename )
        %SAVE_CSV ��������� � csv ��������
        %   ����:
        %   features    - ������� ���������
        %   meta        - ��������� ����������� ��������� (��� ������� [������])
        %   f_header    - ���������
        %   filename    - ��� �����
        %   � ������ ������ ���������� �������� ���������, � t ����������
        %   ������. ���, ��� ��������� 
                       
            [ n_obj, ~ ] = size(features);
            %meta = RPTools.meta2str(meta); % ������������ ���� ����������� ����� �
            %���������
            %meta = {meta.name}; % 
            fid = fopen(filename,'w');
            fprintf(fid,'%s',f_header);
            fprintf(fid,'%s\n','alfa,beta,t');
            for i_line = 1:n_obj
                fprintf(fid,'%f,',features(i_line,:));
                fprintf(fid,'%d,%d,%s\n',...
                    meta(i_line).asp_a,...
                    meta(i_line).asp_b,...
                    meta(i_line).name);
            end
            fclose(fid);
            
            disp(['���� ' filename ' ��������'])
            
        end
            
        % ������� ����� ���
        function [icrs, meta] = open_icr_file(varargin)
        %open_icr_file ��������� ����, ��� ��������� ������
        %   ����:
        %   filename � ������ ��� cell array ����� (������������)
        %   �����:
        %   icrs � ������� ����������� n_objects x n_counts, ��� n_objects � ����� ��������
            obj_names = {'sph' 'con' 'cyl' 'brk' 'dsk' 'plt'};
            % ----- ��������� ��������� ����� 
            if nargin == 0      % ������ ����
                filelist = RPTools.fGetFileList('icr').';
            else 
                if iscell(varargin)
                    filelist = varargin{1};
                else; error('����� ���������� ������ cellarray')
                end
            end
            % ----- �������� ������
            nFiles = size(filelist,1);     % ��������� ����� ������
            meta = cell(nFiles,2);  % ������ � ���������������
            icrs = cell(nFiles,1);  % ������ � ��������� ���
            for i_file = 1:nFiles
                fileID = fopen(filelist{i_file});   % ��������� ����� ���������� ID
                icrs{i_file} = textscan(fileID, '%f32');          % ������� �� ����� ������� ���(float32)
                icrs{i_file} = icrs{i_file}{1}.';
                textscan(fileID, '%s',6);           % ���������� ��������� ��������������
                f_meta = textscan(fileID, '%u %f32 %f32 %f32 %f32 %f32');   % ��������� ��������������
                fclose(fileID);                         % �������� ID
                % 1) N_obj; 2) H_obj; 3) R_obj; 4) L_elm; 5) Fi_az; 6) Fi_um
                meta{i_file, 1} = [ obj_names{f_meta{1}} '('... 
                    num2str(f_meta{2}) 'x' num2str(f_meta{3}) ')' ];                % ����� �����
                meta{i_file, 2} = f_meta{5}; meta{i_file, 3} = f_meta{6};   % ���� �������
            end
            icrs = vertcat(icrs{:});
            meta = struct('name', meta(:,1), 'asp_a', meta(:,2),'asp_b', meta(:,3));
            icrs = fillmissing(icrs,'movmean',5,2);     % ��������� NaN ��������
            disp(['������� ' num2str(nFiles) ' ������'])
        end
        
        
        % �������� ����� ��
        function [rps, meta] = open_dat_file(varargin)
        %open_dat_file ������� ���� dat � ��
        %   ����:
        %   filename � ������ ��� cell array ����� (������������)
        %   �����:
        %   rps � �� ��� ��������� ��
        %   n_counts � ����� �������� � ��
            % ----- ��������� ��������� ����� 
            if nargin == 0      % ������ ����
                filelist = RPTools.fGetFileList('dat').';
            else 
                if iscell(varargin)
                    filelist = varargin{1};
                else; error('����� ���������� ������ cellarray')
                end
            end
            % ----- �������� ������
            nFiles = size(filelist,1);     % ��������� ����� ������
            meta = cell(nFiles,2);  % ������ � ���������������
            for i_file = 1:nFiles
                filename = filelist{i_file}; % �������� ��� �����
                fileID = fopen(filename);   % ��������� ����� ���������� ID
                rp = textscan(fileID, '%u %f32 %f32');          % ������� �� ����� ������� ���(float32)
                if (i_file==1); rps = zeros(nFiles, length(rp{2})); end
                rps(i_file,:) = rp{2}.';
                fclose(fileID);                         % �������� ID
                filename = filename(find(filename=='\', 1, 'last')+1:end);
                meta{i_file, 1} = filename(1:2);   % ����� �����
                meta{i_file, 2} = str2double(filename(4:5)); 
                meta{i_file, 3} = str2double(filename(7:8));   % ���� �������
            end
            meta = struct('name', meta(:,1), 'asp_a', meta(:,2),'asp_b', meta(:,3));
            disp(['������� ' num2str(nFiles) ' ������'])
        end
        
        % �������� ������ ������ � cellarray
        function [output_args] = fGetFileList(ext)
        %FGETFILELIST �������� ������ ������
        %   ������� ���������� ����������� ��������� ��� ������ ������ �
        %   ��������� �� ������ ������ cell � ������� ������� ������ ���
        %   cell � �����
            [filenames,path] = ...
                uigetfile(...
                ['*.', ext],'MultiSelect','on');  % ������� ����� ����� ���
            if ~iscell(filenames)                       % ������� <2 ������
                if (size(filenames,2)==1);filenames={0};% 0 ������
                else; filenames = {[path,filenames]};   % 1 ���� - ��������� ������ ���
                end
            else                                        % ��������� ������
                for i=1:(size(filenames,2))
                    filenames{i} = [path,filenames{i}];
                end
            end
            output_args = filenames; % ������� ������
        end
        
        % ��������� ��� + ���������  �� �����
        function [ output_args ] = fOpenICR( filename )
        %FOPENICR ��������� ��� �������
        %   ��������� ��� ������� �� ����� filename
        %   output_args{1:7} ��������� � �������: 
        %   1       2       3       4       5       6
        %   N_obj   H_obj   R_obj   L_elm   N_cnts  [Fi_az Fi_um] 
        %   7
        %   icrValues
            fileID = fopen(filename);               % ��������� ����� ���������� ID
            icrCell = textscan(fileID, '%f32');     % ������� �� ����� ������ ����� float32
            textscan(fileID, '%s',6); 
            Parameters = textscan(fileID, '%u %f32 %f32 %f32 %f32 %f32');    % ��������� ����� �����
            fclose(fileID);                         % �������� ID
            output_args(1:4) = Parameters(1:4);     % ����� ����������
            output_args{5} = size(icrCell{1,1},1);  % ����� ����� �������� � ���
            output_args{6} = horzcat(Parameters{5:6});     % ����� ���� �������
            output_args{7} = icrCell{1,1}.';        % ����� ��������  
        end
    
    % ������������ ��������������
        
        % ����������� ��
        function [ FeatDataF ] = afft( inData, varargin )
        %AFFT �������������� � ������������ ��������� �����
        %   ���������� ���������� �������������� � ������������ ��������� �����
        %   �� ����� ������� inData: ������ ������ - ��
        %   �������� ��������� � ������� �� n-� - ����� ���������������;
        %   ���� ����� ��������������� ���������� 
        %   ����� ��������� ������������ ����� ������� ��������� ����������
        %   ��� ������� ������������ ������������ ������� �� sqrt2.
            [ ~ , nDim] = size(inData);
%             nDim = 2*nDim;  % �������� ��� ����������� ����� �������� �� � ��
            if nargin==2
                nDim = varargin{1};
            end
            FeatDataF = abs(fft(inData,nDim,2))/sqrt(nDim);                 % ��� ����� ������� ���������
            % ��������� ��������� ���������, �.�. f� ������������
            % ������. ������� ���������, �.�. ��-�� ����������
            % ��� ����� �������� ���������� (��������� ���� � ��������)
            FeatDataF = FeatDataF(:,1:(floor(nDim/2)+1));
            FeatDataF(:,2:ceil(nDim/2)) = sqrt(2)*FeatDataF(:,2:ceil(nDim/2)); % ���������� �� sqrt ��-�� ���������� ��������
        end
        
        % �������������� ����� (������������� ������)
        function [ FeatDataF ] = fft( inData )
        %FFFT �������������� � ������������ ��������� �����
        %   ���������� ���������� �������������� � ������������ ��������� �����
        %   �� ����� ������� inData: ������ ������ - ��
        %   �� ����� ������ 
        %   1 [0� re ����]  2 [1� re ����] 3 [1� im ����] ... 2n [n-� ����(��������������)]
            [nObj,nDim] = size(inData);
            if mod(nDim,2); warning('�������� ����� ��������'); end
            FeatDataFC = fft(inData, nDim, 2)/sqrt(nDim);           % ��� ����� ������� ���������
            FeatDataF = zeros(nObj,nDim,'like',inData);
            FeatDataF(:,1) = FeatDataFC(:,1);           % 0-� ���������
            FeatDataF(:,end) = FeatDataFC(:,nDim/2+1);  % ��������� ���������
            FeatDataF(:,2:2:end-2) = real(FeatDataFC(:,2:nDim/2));  % �������������� ����������
            FeatDataF(:,3:2:end-1) = imag(FeatDataFC(:,2:nDim/2));  % ������ ����������
            FeatDataF(:,2:(end-1)) = sqrt(2)*FeatDataF(:,2:(end-1)); % ������������� ������������ �� sqrt 2
        end
        
        % ����������� ��
        function [ features, coef_map ] = cwt( rps, varargin )
        %CWT ��������� ���
        %   ���������� ������� � ��, ����� �������� �������� ��������
        %   ����� �������� �������������:
        %   [�� 1 1i 2 2i...] ... [�� 1 1i 2 2i... ]
            if nargin == 1
                wave_name = 'dtf2';
            elseif nargin == 2
                wave_name = varargin{1};
            else, error('������� ����� ����������')
            end
            
            [rps, coef_map] = RPTools.fDTCWT(rps, wave_name);  % ������������ ����������� ������������
            [n_rps, wt_dim] =size(rps);   % ���������� ����������� ������� ��������� cwt
            features = zeros(n_rps, wt_dim);
            features(:,1:2:end) = rps(:,1:wt_dim/2); 
            features(:,2:2:end) = rps(:,wt_dim/2+1:wt_dim);
            coef_map = coef_map * 2;      % ������������� �������� ��-�� �������������
        end
        
        % ����������� ���
        function [ FeatDataW, coef_map ] = acwt( InData, varargin )
        %FACWT ������������ ����������� ������������ ������������ �������
        %��������������
            if nargin == 1
                wave_name = 'dtf2';
            elseif nargin == 2
                wave_name = varargin{1};
            else, error('������� ����� ����������')
            end
            [FeatDataWC, coef_map] = RPTools.fDTCWT(InData,wave_name);  % ������������ ����������� ������������
            nDim =size(FeatDataWC,2);   % ���������� ����������� ������� ��������� cwt
            FeatDataW = abs( FeatDataWC(:,1:nDim/2)+...
                1i*FeatDataWC(:,nDim/2+1:nDim));
            %coef_map = flip(coef_map, 2);
            
        end
        
        % ��� - ������� ������� ������
        function [ FeatDataW, coef_map ] = fDTCWT( InData, wavenm )
        %fDTCWT ����������� ������� �������������� �� ������ ���������
        %   ������� ��������� ����������� �������-�������������� ������� ������ �� 
        %   ����� ����� ��������� ����� �������� 
        %   �������������� �������� ���������� � ��������� �� ������ �����
        %   �������� � 2^m ��� ����������� ����� �������� �� �����
        %   ��� ���� ����������� ���������� ������ �� ����������� ������
        %   ��������� ������� ������.
        %   �������� 1-������ ������ ����� ������
        %      1[re �� ���� �� 1] ... [re �� ���� �� n]m,    m+1[re �� ����]n
        %    n+1[im �� ���� �� 1] ... [im �� ���� �� n]n+m,n+m+1[im �� ����]2n
            [ nObj,nDim ] = size(InData);   % �����������
            decfilt = dtfilters(wavenm);    % ��������� ���� ������� �� ��������
            filtMaxLen = max([...
                size(decfilt{1}{1},1) size(decfilt{1}{2},1)...
                size(decfilt{2}{1},1) size(decfilt{2}{2},1)]);  % ����������� ������������ ����� ��������
            % ����������� ������ ����������
            declvl = floor(log2(nDim/filtMaxLen)+1);    % ������� ���������� �� ������� x_len>=flen*2^(L-1)
            newDim = 2^declvl * ceil(nDim/2^declvl);    % ����� ����������� 
            
            dimAdd = newDim-nDim;           % ��������� �������
            if (dimAdd~=0)                  % ���� ������� �� �������
                InData = [InData , zeros(nObj,dimAdd)]; % ��������� ������� ������ ��� ���������� �������
                nDim = newDim;                      % � �������� �����������
            end
            % ���������� ������������� (�������: ����� ������������� = nDim �2
            FeatDataW = zeros(nObj,2*nDim);         % ��������� �������
            for iObj = 1:nObj                       % ��� ���� ��������
                cwtData = dddtree(...
                    'cplxdt',InData(iObj,:),declvl,wavenm);         % ��������� �������������� 
                cfs = cellfun(@squeeze, cwtData.cfs,'UniformOutput',0); % ����� ���������
                cfs = vertcat(cfs{:});
                FeatDataW(iObj,1:2*nDim) = [cfs(:,1); cfs(:,2)].';  % ��������� ��. � ��. ������������
            end
            coef_map = cellfun(@length,cwtData.cfs);
        end
        
        % ��
        function [ cfs, coef_map ] = wt( rps, varargin )
            [n_rps, dim_rp] = size(rps);
            if nargin == 1
                wave_name = 'sym4';
            elseif nargin == 2
                wave_name = varargin{1};
            else, error('������� ����� ����������')
            end
            dec_lev = wmaxlev(dim_rp, wave_name);   % ������������ ������� ������������
            [wt_fset, coef_map] = wavedec(rps(1,:), dec_lev, wave_name); % ���������� ������� ��
            coef_map = coef_map(1:end - 1);         % ��������� ��������� ������� (����� ���. �������)
            cfs = zeros(n_rps, length(wt_fset));    % ������� ������� �������������
            cfs(1,:) = wt_fset; % ����������� ������������ ������� ������������ ��
            for i = 2:n_rps     % ���������� ������������
                cfs(i,:) = wavedec(rps(i,:), dec_lev, wave_name);
            end
        end
        
        % ����
        function [ signal ] = ifft(fft_fvs)
            
            n_dim = size(fft_fvs,2);
            rec_spectre = sqrt(n_dim)*[ ...
                fft_fvs(:,1)    (fft_fvs(:,2:2:end-2)+1i*fft_fvs(:,3:2:end-1))/sqrt(2) ...
                fft_fvs(:,end)  (fft_fvs(:,end-2:-2:2)-1i*fft_fvs(:,end-1:-2:3))/sqrt(2) ];
            signal = ifft(rec_spectre, n_dim, 2);
        end
        
        % PCA
        function [varargout] = pca(rps)
        %pca �������������� ��������-�����
            coeff = pca(rps);       % ������� ��������������
            pca_rps = rps*coeff;    % ��������������
            varargout{1} = pca_rps;
            if nargout>1, varargout{2} = coeff; end
        end
        
        function [dca_rps] = dca_total(rps,trgs)
        %dca ����� ������� ���������� ���������
        % ����������� ������ ���������� �������� ����� ��������� ���������, �������������� ������ ��������
        % ��� ��������� �������� ����������� ��������� ���������. 
        % ���������� ��������� ���������� � ���, � ���������� ���� ����������� �������������� ���������
        % ���������:
        %   rps - ������� ��������� (��)
        %   trgs - ������� �������
        % �����:
        %   ��������� ��������������[, ������������ ��������������: x*coeff]
            [object_names, ~, ci] = unique(trgs,'rows');
            [~, dim_rps] = size(rps);
            n_obj = size(object_names,1);
            rho_cell = cell(1,(n_obj^2));
            for i_obj = 1:(n_obj-1)
                for j_obj = (i_obj+1):n_obj
                    
                    c1 = rps(ci==i_obj, :);    % ������������ ������ �������
                    c2 = rps(ci==j_obj, :);    % ������������ ������� �������
                    n_fvs = [size(c1,1) size(c2,1)];    % ����������� ����� ����������
                    vrhos = repmat(c1.',1,1,n_fvs(2)) - repmat(...
                        permute(c2,[2 3 1]), 1, n_fvs(1), 1);   % ������ ���������� ��������
                    vrhos = reshape(vrhos, dim_rps, []);
                    rho_cell{i_obj,j_obj} = vrhos.';
                end
            end
            clear vrhos;
            rho_cell = vertcat(rho_cell{:});
            rho_cell = vertcat(rho_cell, -rho_cell);
            coeff = pca(rho_cell);  % ������� ��������������
            dca_rps = rps*coeff;    % ��������������
            
        end
        
        % DCA
        function [varargout] = dca(rps, trgs)
        %dca ����� ������� ��������� ���������� ���������
        % ����������� ������ ���������� �������� ����� ���������� ��������� ���������, �������������� ������ ��������
        % ��� ��������� �������� ����������� ��������� ���������. 
        % ���������� ��������� ���������� � ���, � ���������� ���� ����������� �������������� ���������
        % ���������:
        %   rps - ������� ��������� (��)
        %   trgs - ������� �������
        % �����:
        %   ��������� ��������������[, ������������ ��������������: x*coeff]
        
            [object_names, ~, ci] = unique(trgs,'rows');
            [~, dim_rps] = size(rps);
            n_obj = size(object_names,1);
            minrho = cell(1,(n_obj^2));
            for i_obj = 1:(n_obj-1)
                for j_obj = (i_obj+1):n_obj
                    
                    c1 = rps(ci==i_obj, :);    % ������������ ������ �������
                    c2 = rps(ci==j_obj, :);    % ������������ ������� �������
                    n_fvs = [size(c1,1) size(c2,1)];    % ����������� ����� ����������
                    vrhos = repmat(c1.',1,1,n_fvs(2)) - repmat(...
                        permute(c2,[2 3 1]), 1, n_fvs(1), 1);   % ������ ���������� ��������
                    
                    rhos = sum((vrhos).^2, 1);              % ������ ������ ( ���� 1 � n_fvs(1) x n_fvs(2) )
                    [~, i_obj_mins] = min(rhos,[],3);
                    [~, j_obj_mins] = min(rhos,[],2);
                    
                    vrhos = reshape(vrhos, dim_rps, []);
                    i_obj_mins = sub2ind(size(squeeze(rhos)), (1:n_fvs(1)).', i_obj_mins(:));
                    j_obj_mins = sub2ind(size(squeeze(rhos)), j_obj_mins(:), (1:n_fvs(2)).');
                    minrho{i_obj,j_obj} = vrhos(:, i_obj_mins).';
                    minrho{j_obj,i_obj} = vrhos(:, j_obj_mins).';
                    
                end
            end
            clear vrhos;
            minrho = vertcat(minrho{:});
            minrho = vertcat(minrho, -minrho);  % �������� �������
            coeff = pca(minrho);    % ������� ��������������
            dca_rps = rps*coeff;    % ��������������
            varargout{1} = dca_rps;
            if nargout>1, varargout{2} = coeff; end
        end
        
    % ����������� � ������

        % ��������� ����� ��������������/�������� ���������
        function [ feature_map ] = get_feature_map( transform_type, coef_map )
        %GET_FEATURE_MAP ������������ �������� �������������
        %   ����:
        %   transform_type  - �������� ��������������
        %   coef_map        - ����� �������������
            switch transform_type
                case 'none'
                    feature_map = num2str(0:coef_map-1, 't_%03d,');
                case 'afft' % ������ ���� �������� ����������� ��������
                    feature_map = num2str(0:coef_map-1, 'f_%03d,');
                
                case 'fft' 
                    feature_map = [ 'f_000r,'...
                        num2str(reshape([1:coef_map/2-1;1:coef_map/2-1],1,[]), 'f_%03dr,f_%03di,') ...
                        num2str(coef_map/2, 'f_%03dr,')];
                    
                case 'cwt'  % ������ ���� �������� ����� �������������
                    % ��������� ���������� ������� ������������
                    coef_map = coef_map/2;
                    max_lev = length(coef_map)-1;
                    coefs = cell(max_lev+1,1);
                    for m = 1:max_lev+1
                        cur_lev = mod(m,max_lev+1);
                        m_vect = cur_lev+zeros(coef_map(m), 1);    % ������ ���������� ���������
                        n_vect = (0:(coef_map(m)-1)).';     % ������ ��������� �������
                        if coef_map(m)<100, formatSpecs = 'c_%02d_%02d';
                        elseif coef_map<1000, formatSpecs = 'c_%03d_%03d';
                        else, error('������� ����� �������������')
                        end
                        coefs{m,1} = [ ... 
                            num2str([m_vect n_vect], [formatSpecs '_r,']) ...
                            num2str([m_vect n_vect], [formatSpecs '_i,']) ];
                        coefs{m,1} = reshape(coefs{m,1}.',1,[]);
                    end
                    feature_map = horzcat(coefs{:,1});
                    
                case 'acwt' %
                    max_lev = length(coef_map)-1; % ������������ ������� ����������
                    coefs = cell(max_lev+1,1);
                    for m = 1:max_lev+1
                        cur_lev = mod(m,max_lev+1);
                        m_vect = cur_lev+zeros(coef_map(m), 1);    % ������ ���������� ���������
                        n_vect = (0:(coef_map(m)-1)).';     % ������ ��������� �������
                        if coef_map(m)<100, formatSpecs = 'c_%02d_%02d';
                        elseif coef_map<1000, formatSpecs = 'c_%03d_%03d';
                        else, error('������� ����� �������������')
                        end
                        coefs{m,1} = num2str([m_vect n_vect], [formatSpecs ',']);
                        coefs{m,1} = reshape(coefs{m,1}.',1,[]);
                    end
                    feature_map = horzcat(horzcat(coefs{:,1}));
                    
                case 'wt'   % ��� ������� ��������������
                    max_lev = length(coef_map)-1; % ������������ ������� ����������
                    coefs = cell(max_lev+1,1);
                    for m = 1:max_lev+1
                        cur_lev = mod(m,max_lev+1);
                        dim_cur_lev = coef_map(end-m+1);
                        m_vect = cur_lev+zeros(dim_cur_lev, 1);    % ������ ���������� ���������
                        n_vect = ( 0:(dim_cur_lev -1) ).';     % ������ ��������� �������
                        if coef_map(m)<100, formatSpecs = 'c_%02d_%02d';
                        elseif coef_map<1000, formatSpecs = 'c_%03d_%03d';
                        else, error('������� ����� �������������')
                        end
                        coefs{end-m+1,1} = num2str([m_vect n_vect], [formatSpecs ',']);
                        coefs{end-m+1,1} = reshape(coefs{end-m+1,1}.',1,[]);
                    end
                    feature_map = horzcat(horzcat(coefs{:,1}));
                case 'pca'
                    feature_map = num2str(0:coef_map-1, 'pca_%03d,');
                case 'dca'
                    feature_map = num2str(0:coef_map-1, 'dca_%03d,');
                otherwise
                    warning('����������� ��������������');
            end
            
        end
          
        % ������������� ���������� � ������
        function [meta_str] = meta2str( meta )
            %META2STR 
            n_rps = size(meta, 1);
            meta_str = cell(n_rps,1);
            for i = 1:n_rps
                meta_str{i} = horzcat(meta(i).name, ...
                        num2str(meta(i).asp_a,'%02u'), ' ',...
                        num2str(meta(i).asp_b,'%02u') );
            end
        end
        
        % ������� ��������� ���������� �������
        function plot_rps(h, rps, meta, list)
            [~, dim_rp] = size(rps);
            obj_names = unique({meta.name});
            obj_marks = {'','','s','^','v','*'};
            asp_colors = repmat(linspace(0, 0.4, 10), 3, 1).';
            asp_widths = linspace(0.8, 1.4, 10);
            asp_lines =  {'--', '-', '-.', ':', '--', '-.', '-', ':', '--', '-.'};
            hold on
            legend_str = cell(length(list),1);
            for i = 1:length(list)
                i_rp = list(i);
                name = meta(i_rp).name; asp_a = meta(i_rp).asp_a; asp_b = meta(i_rp).asp_b;
                curve = plot(h,...
                    1:dim_rp, rps(i_rp,:),...
                    [':' obj_marks{(contains(obj_names, name))} 'k']...
                    );
                asp_ind = i;%floor(mod(ceil(max([asp_a, asp_b]))/15 -1, 10) +1);
                set(curve, ...
                    'MarkerSize', 5,...
                    'Color', asp_colors(asp_ind, :),...
                    'LineWidth', asp_widths(asp_ind),...
                    'LineStyle', asp_lines{asp_ind});
                legend_str{i} = [ name '^{'...
                    num2str(asp_a) '^\circ,'...
                    num2str(asp_b) '^\circ}'];  % ������ �������� ��
%                 if i <= 6     % ��������� ����� �� ���. ���� �� ����� �� 6
%                     [y, x] = max(rps(i_rp, :));    % ����� ��� ���������� ������ (���� ��������)
%                     text(cast(x, 'double'), cast(y, 'double'), [ legend_str{i}],'VerticalAlignment','bottom');
%                 end
            end
            set(gca,'FontName','ISOCPEUR','FontAngle','italic',...
                'xGrid','on','yGrid','on',...
                'XMinorTick','on')
            if length(list) < 7; legend(legend_str); end
            
        end
        
        % ����� ������������� ������������ ������� ��������������
        function plot_cwt(h, cfs, coef_map, meta, varargin)
            if nargin == 4
                list = 1:size(cfs, 1);
            elseif nargin == 5
                list = varargin{1};
            else, error('������������ ����� ����������')
            end
            colors = 'rgbmcykw';
            
            n_samples = length(list);
            legend_list = cell(n_samples, 1);
            hold on;
            max_lvl = length(coef_map);
            h_vector = zeros(n_samples, 1);
            for i_sample = list
                start = 1;
                for i = 1:max_lvl
                    stop = start + coef_map(i) - 1;
                    z_axis = cfs(i_sample, start:stop);
                    x_axis = linspace(1,coef_map(end),coef_map(i)); 
                    y_axis = i * ones(coef_map(i),1) - 1 +...
                        0.5 * i_sample / length(list);
                    h_line = plot3(h, x_axis, y_axis, z_axis,['-+' colors(i_sample)]);
                    %[z, max_i] = max(z_axis);    % ����� ��� ���������� ������ (���� ��������)
                    %y = y_axis(1);
                    %text(x_axis(max_i), y, z, meta(i_sample).name); %
                    start = stop + 1;
                end
                h_vector(i_sample) = h_line;
                legend_list{i_sample} = horzcat(meta(i_sample).name, ...
                    num2str(meta(i_sample).asp_a), ',',...
                    num2str(meta(i_sample).asp_b) );
            end
            xlabel('��� ��������'); ylabel('������� ����������');   zlabel('�������� ���������');
            legend(h_vector, legend_list);
            zlim([min(min(cfs)), max(max(cfs))])
            objAxe = gca;
            objAxe.YTick = 0:max_lvl;
            %objAxe.YColor = [0.6 0.6 0.6];
            objAxe.GridAlpha = 1;
            grid on; axis vis3d; view(-45, 30);
            camzoom(0.6);
            objAxe.XGrid = 'off'; objAxe.ZGrid = 'off';
%             objAxe.XDir = 'reverse';
            objAxe.GridLineStyle = ':';
        end
                
    end
    
end
