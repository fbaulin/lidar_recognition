classdef Slozhno < handle
    %Slozhno Support functions for slozhno experiment
    %   Slozhno experiment
    %#ok<*PROPLC>
    properties
                
        window_m = 3;           % максимальная протяженность объекта
        icr_dim = 200;          % число отсчетов в ИХР
        hwidth2sig = 2.35;      % множитель для перевода сигмы в длину по половинному уровню
        resol_ns                % условное разрешение ИХР
        rp_dim
        
        scan_type               % тип сканирования устанавливается пользователем
        obj_lst                 % список объектов
        imp_hwidths_ns          % вектор длительностей ЗИ по половинной высоте
        t_type                  % тип используемого преобразования
        wave_name               % тип вейвлета
        cwave_name
        band_factors            % отношения fd/П
        
        f_sampling_mhz          % частота дисретизации, МГц
        T_sampling_ns           % интервал дискретизации, нс
        f_imp_mhz               % полоса сигнала, МГц
        f_lpfilter_mhz          % полоса фильтра, МГц
    end
        
    methods
        
        % конструктор
        function obj = Slozhno() 
            
            obj.resol_ns = obj.window_m * 2 / obj.icr_dim / 3e8 * 1e9;      % условное временное разрешение ИХР, нс
            obj.obj_lst ={ 'brk3x1' 'con2x1','con3x1','cyl3x1','sph3x1' };
            obj.scan_type = 'default';      % по умолчанию 
            obj.t_type = 'none';            % по умолчанию нет преобразования
            obj.wave_name = 'dtf2';         % по умолчанию комплексный вейвлет
        end
        
        % загрузить ИХР одного объекта из icr файлов
        function [icrs, meta] = get_icrs(obj, icr_folder) 
        %GET_ICRS Загрузить ИХР и метаданные (углы и пр.) из папки
            %   Аргументы:
            %   icr_folder -    - Папка с ИХР
            %   obj_id          - Имя для БД
            % Загрузка ИХР
            mat_files = dir([icr_folder '/*.mat']);     % получить набор ИХР
            if size(mat_files,1)~=0     % есть мат-файлы              
                load([mat_files(1).folder '\' mat_files(1).name],'icrs','meta');   % загрузить ИХР из mat файла
            else
                filelist = dir([icr_folder '/*.icr']);      % получить структуру списка файлов из заданной папки
                for i = 1:length(filelist)                  % собрать список полных имен файлов из структуры
                    filelist(i).name = horzcat(filelist(i).folder, '\', filelist(i).name); 
                end
                [icrs, meta] = RPTools.open_icr_file(...
                    {filelist(1:end).name}.');              % получить ИХР и метаданные
                if size(icrs,2)~=obj.icr_dim
                    warning(['Число отсчетов в ИХР не равно ' num2str(obj.icr_dim)]);
                end
                icrs_db_filename = [filelist(i).folder '\' meta(1).name '.mat'];
                save(icrs_db_filename,'icrs','meta');
                disp(['Сохранение базы ИХР в ' icrs_db_filename]);
            end
        end
        
        % загрузка ИХР нескольких объектов из mat файлов
        function [icrs, meta] = get_all_icrs(obj)
            n_obj = length(obj.obj_lst);    
            icrs = cell(1,n_obj);   % матрица под ИХР
            meta = cell(1, n_obj);  % матрица под метаданные
            for i_obj = 1:n_obj     % загрузка ИХР 
                icr_folder = [ 'icrs/' obj.scan_type '/' obj.obj_lst{i_obj}];       % задать путь папки с ИХР
                %icrs_data = load([icr_folder '/' obj.obj_lst{i_obj} '.mat']);   % загрузить ИХР из mat файла
                [icrs{i_obj}, meta{i_obj}] = obj.get_icrs(icr_folder);
                
            end
        end
        
        % модель приемного тракта
        function rps = lls_model(obj, icrs, imp_hwidth_ns, lp_band_mhz, fd_mhz_q) 
        %LLS_MODEL Смоделировать работу ЛЛС
        % icrs - ИХР, расположенные по горизонтали
        % imp_hwidth_ns - длительность импульса по половине высоты
        % lp_band_mhz - полоса НЧ фильтра
        % fd_mhz_q - массив ячеек {(fd, МГц),(битность)}, 
        % 'n' отключает дискретизацию и квантование
            % формирование ДП
            
            fprintf('\tllis_model: ')
            imp_hwidth = ceil(imp_hwidth_ns/obj.resol_ns);  % длительность импульса
            if ~strcmp(imp_hwidth_ns,'n')
                if imp_hwidth_ns>0.2                            % если ЗИ приемлемой длительности
                    rps = RPTools.make_rps(icrs, imp_hwidth);   % синтез ДП
                    fprintf('[расс. ЗИ с t=%0.2fнс]',imp_hwidth_ns)
                else                    % иначе длительность ЗИ слишком мала:
                    rps = icrs;         % принять ИХР за ДП
                    fprintf('[ДП - ИХР]');
                end
            else
                rps = icrs;
                fprintf('[загр. ДП]');
            end
            [n_rps, obj.rp_dim] = size(rps);
            % фильтрация
            if ~strcmp(lp_band_mhz,'n')                         % если передано 
                n_filter = 2;
                rps = Slozhno.bworth_filter(...
                    rps, obj.resol_ns, lp_band_mhz, n_filter);  % фильтрация
                fprintf('[ФНЧ(%dп) Fф=%0.1fМГц]', ...
                    n_filter, lp_band_mhz);
                obj.f_lpfilter_mhz = lp_band_mhz;
            end
            % АЦП
            if iscell(fd_mhz_q)     % дискретизация и квантование
                fd_mhz  = fd_mhz_q{1};  % извлечение частоты дискретизации
                q_bit   = fd_mhz_q{2};  % извлечение числа уровней квантования
                if ~strcmp(fd_mhz,'n')  % дискретизация сигнала
                % Реализовать ресэмплинг
                    fd_src_mhz = 1/obj.resol_ns*1e3;    % исходная частота дискрет. (ИХР)
                    rps = RPTools.spline_interpol(rps, fd_src_mhz, fd_mhz);
                    fprintf('[АЦП Fд=%dМГц', floor(fd_mhz));
                    obj.f_sampling_mhz = fd_mhz;        % обновить частоту дискретизации
                    obj.T_sampling_ns = 1/fd_mhz*1e3;   % обновить разрешение по времени
                    
                end
                if ~strcmp(q_bit,'n')   % квантование 
                    fprintf(',q=%d', q_bit);
                    rps = RPTools.adc(rps,q_bit);                    
                end
                fprintf(']');
            else, warning('    Пропуск операции АЦП');    
            end
            if mod(size(rps,2),2)==1;rps=[rps, zeros(n_rps,1)];end
            fprintf('\n');
            %rps = RPTools.norm_rp(rps,'mod');  % нормировка по модулю
        end
        
        % гистограмма метрик
        function [bms_counts] = gpu_histcounts(~, fvs, hist_bins) 
        %gpu_histcounts
            n_obj = length(fvs);
            bms_mod = cell(n_obj, n_obj);       % матрица под метрики
            bms_counts = cell(n_obj, n_obj);    % матрица отсчетов внешних метрик
            gpu_hist_bins = gpuArray(hist_bins);
            for i_obj = 1:(n_obj-1)             % расчет метрик
                fvs_1 = gpuArray(fvs{i_obj});   % загрузка вектора признаков в ГП
                %dim_fv = size(fvs_1,2);
                fv_mods = sqrt(sum(fvs_1.^2,2)); %repmat(sqrt(sum(fvs_1.^2,2)), 1, dim_fv);     % выполнить нормировку векторов на их длину
                fvs_1 = fvs_1./fv_mods;                   % нормировка для расчета угловых метрик
                for j_obj = (i_obj+1):n_obj
                    fvs_2 = gpuArray(fvs{j_obj});
                    fv_mods = sqrt(sum(fvs_2.^2,2));%repmat(sqrt(sum(fvs_2.^2,2)), 1, dim_fv);
                    fvs_2 = fvs_2./fv_mods; % нормировка для расчета угловых метрик
                    clear fv_mods
                    bms_mod{i_obj,j_obj} = 2*asin(...
                        Slozhno.calc_b_metrics(fvs_1,fvs_2,8)/2); % расчет угловых метрик
                    bms_counts{i_obj,j_obj} = ...
                        histcounts(bms_mod{i_obj,j_obj}, gpu_hist_bins );
                    bms_counts{i_obj,j_obj}=gather(bms_counts{i_obj,j_obj});
                    
                end
            end
        end
        
        % гистограмма метрик
        function [bms_counts] = gpu_linear_histcounts(~, fvs, hist_bins) 
        %gpu_histcounts
            n_obj = length(fvs);
            bms_mod = cell(n_obj, n_obj);       % матрица под метрики
            bms_counts = cell(n_obj, n_obj);    % матрица отсчетов внешних метрик
            gpu_hist_bins = gpuArray(hist_bins);
            for i_obj = 1:(n_obj-1)             % расчет метрик
                fvs_1 = gpuArray(fvs{i_obj});   % загрузка вектора признаков в ГП
                for j_obj = (i_obj+1):n_obj
                    fvs_2 = gpuArray(fvs{j_obj});
                    bms_mod{i_obj,j_obj} = Slozhno.calc_b_metrics(fvs_1,fvs_2); % расчет угловых метрик
                    bms_counts{i_obj,j_obj} = ...
                        histcounts(bms_mod{i_obj,j_obj}, gpu_hist_bins );
                    bms_counts{i_obj,j_obj}=gather(bms_counts{i_obj,j_obj});
                    
                end
            end
        end
        
        % метрики по компонентам
        function [] = gpu_bms_vectors(obj, fvs, meta, imp_hwidth_ns) 
        %GPU_BMS_VECTORS Формирование и созданение разностных векторов на
        %ГП
        %   gpu_bms_vectors(fvs, meta)
        %   fvs - ячейки с множествами ДП
        %   meta - ячейки с метаданными для сохранения
            n_obj = length(obj.obj_lst);
            for i_obj = 1:(n_obj-1)
                fvs_1 = gpuArray(fvs{i_obj});   % загрузка вектора признаков в ГП
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
       
        % преобразование
        function [ fvs ] = transform(obj, rps, varargin)
        %TRANSFORM Сформировать вектора признаков
        %   Выполняет преобразование с использованием методов класса RPTools
        %   Аргументы:
        %   .obj        
        %   rps         - матрица ячеек ДП
        %   [t_type]    - тип преобразования
        %   [wave_name] - обозначение вейвлета
            if nargin>=3;   t_type = varargin{1};   % тип преобразования
            else; t_type = obj.t_type; end 
            if nargin>=4; wave_name = varargin{2};  % тип вейвлета
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
        
        % получить названия признаков для заданного преобразования
        function [ fvs, coef_map ] = transform_wmap(obj, rps, varargin)
        %TRANSFORM Сформировать вектора признаков
        %   Выполняет преобразование с использованием методов класса RPTools
        %   Аргументы:
        %   .obj        
        %   rps         - матрица ячеек ДП
        %   [t_type]    - тип преобразования
        %   [wave_name] - обозначение вейвлета
            if nargin>=3;   t_type = varargin{1};   % тип преобразования
            else; t_type = obj.t_type; end 
            if nargin>=4; wave_name = varargin{2};  % тип вейвлета
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
        
        % гистограммы внутренних метрик
        function [wms_counts] = wms_hcounts(~, fvs, hist_bins)
            n_obj = length(fvs);
            wms_counts = cell(1, n_obj);
            hist_bins = gpuArray(hist_bins);
            for i_obj = 1:n_obj
                gpu_fvs = (fvs{i_obj});   
                gpu_fvs = gpu_fvs./sum(gpu_fvs,2);  % нормировка для вычисления угловой меры
                wms = asin(Slozhno.calc_w_metrics(gpu_fvs)/2)*2;  % расчет угловой метрики
                clear gpu_fvs;
                wms_counts{1,i_obj} = (...
                        histcounts(wms, hist_bins )); % сбор данных
            end
            clear hist_bins
        end
            
        % вывод гистограм
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
        
        % Сформировать ДП
        function [rps, meta] = generate_rps(icr_folder, gauss_hwidth)
        %GENERATE_RPS Сформировать ДП на основании ИХР
            %   Аргументы:
            %   icr_folder -    - Папка с ИХР
            %   obj_id          - Имя для БД
            
            dim_icr = 200;      % число отсчетов в ИХР
            window_m = 3;       % длина интервала измерений, м.
            resol_ns = window_m * 2 / dim_icr / 3e8 * 1e9;  % условное временное разрешение ИХР, нс
            % Загрузка ИХР
            filelist = dir([icr_folder '/*.icr']);        % получить структуру списка файлов из заданной папки
            for i = 1:length(filelist)      % собрать список полных имен файлов из структуры
                filelist(i).name = horzcat(filelist(i).folder, '\', filelist(i).name); 
            end
            [icrs, meta] = RPTools.open_icr_file(...
                {filelist(1:end).name}.');  % получить ДП и метаданные
            clear filelist
            % Формирование ДП
            gauss_width_ns = gauss_hwidth * resol_ns;    % длительность импульса в нc
            rps = RPTools.make_rps(icrs, gauss_hwidth);  % синтез ДП
            disp(['Длительность ЗИ по полувысоте: ' num2str(gauss_width_ns) 'МГц']);
        end
            
        % Сформировать вектора признаков
        function [rp_afft,meta] = generate_fvs(icr_folder, gauss_hwidth)
            %GENERATE_FVS Сформировать вектора признаков на основании ИХР
            %   Аргументы:
            %   icr_folder -    - Папка с ИХР
            %   obj_id          - Имя для БД
            
            dim_icr = 200;      % число отсчетов в ИХР
            window_m = 3;       % длина интервала измерений, м.
            resol_ns = window_m * 2 / dim_icr / 3e8 * 1e9;  % условное временное разрешение ИХР, нс
            % Загрузка ИХР
            filelist = dir([icr_folder '/*.icr']);        % получить структуру списка файлов из заданной папки
            for i = 1:length(filelist)      % собрать список полных имен файлов из структуры
                filelist(i).name = horzcat(filelist(i).folder, '\', filelist(i).name); 
            end
            [icrs, meta] = RPTools.open_icr_file(...
                {filelist(1:end).name}.');  % получить ДП и метаданные
            clear filelist
            % Формирование ДП
            gauss_width_ns = gauss_hwidth * resol_ns;    % длительность импульса в нc
            rps = RPTools.make_rps(icrs, gauss_hwidth);  % синтез ДП
            disp(['Длительность ЗИ по полувысоте: ' num2str(gauss_width_ns) 'нс']);
            %mesh(rps); axis vis3d;  % Визуализация ДП

            rp_afft = single(RPTools.afft(rps));                % AFFT
            rp_afft = RPTools.norm_rp(rp_afft, 'energy');
            disp('Получение амплитудных спектров, нормированных по энергии')
            %mesh(rp_afft); axis vis3d;                 % Визуализация спектров
%             Сохранение файла БД
%             fvdb_fname = ['fvdb_' obj_id '.csv'];       % имя файла БД
%             f_header = RPTools.get_feature_map...
%                 ('afft', size(rp_afft,2));              % полчить заголовок для файла
%             RPTools.save_csv...
%                 (rp_afft, meta, f_header, fvdb_fname);  % сохранение базы векторов признаков
%             disp('Сохранение базы векторов признаков')
        end
        
        % Расчет метрик для множества векторов признаков
        function w_met = calc_w_metrics(feature_values)
        %FENORMSW Рассчитывает метрики внутри множества векторов признаков
        %   Рассчитывает метрики между веткорами признаков, 
        %   feature_values  - вектора признаков по строкам матрицы
            [n_aspect,nCounts] = size(feature_values); % число ДП и их размерность
            n_metric = nchoosek(n_aspect,2);  % число метрик
            w_met=zeros(n_metric,nCounts...
                ,'like',feature_values);    % иниц метрик в пр признаков
            i_met = 0;
            disp(['Из ' num2str(n_metric) ' метрик рассчитано:     ']);
            for i=1:n_aspect
                for j=(i+1):n_aspect
                    i_met=i_met+1;
                    w_met(i_met,:) = feature_values(i,:) - feature_values(j,:); % разности признаков
                end
                if i_met<1000000 ;fprintf('\b\b\b\b\b\b%6.0d',i_met); end
            end
            fprintf('\r');
            w_met = sqrt(sum(w_met.^2,2)); % выполняется суммирование по строкам
            % после суммирования элементы полученного вектора-столбца
            % возводятся в квадрат
            
        end
                
        % Расчет метрик
        function [ rhos ] = calc_b_metrics( c1, c2, varargin )
        %fRhoCalc Функция расчета метрики между объектами двух классов
        %   Функция принимает два масива с m(C1Data) и n(C2Data) векторами.
        %   Рассчитывает разности между всеми строками
        %   На выходе mxn матрица евклидовых метрик 
            n_obj = [size(c1,1) size(c2,1)];     % Определение числа объектов в классах
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
        
        % Расчет распределения метрик по измерениям
        function [ c1 ] = calc_b_vectors( c1, c2 )
        %fRhoCalc Функция расчета метрики между объектами двух классов
        %   Функция принимает два масива с m(C1Data) и n(C2Data) векторами.
        %   [mem_limit_MB] 
        %   Рассчитывает разности между всеми строками
        %   На выходе (m x n x v_dim) матрица разностей 
            n_obj = [size(c1,1) size(c2,1)];     % Определение числа объектов в классах
            v_dim = size(c1,2);
            c1 = repmat(c1.',1,1,n_obj(2));      % Репликация матрицы -> v_dim x m x n
            c2 = permute(c2,[2 3 1]);            % транспонирование матрицы -> v_dim x 1 x n
            for i_obj = 1:n_obj(1)
                c1(1:v_dim,i_obj,:) = (c1(1:v_dim,i_obj,:) - c2);   % Расчет метрик
            end
            c1 = permute(c1,[2 3 1]);           % Транспонирование -> m x n x v_dim
        end
        
        % Преобразование индекса метрики в индекс ИХРов
        function [tr_mx] = met2rp(n_rp)
        %met2icr Получить матрицу перехода от метрики к i,j ДП
        %   Получаемая матрица переходов <2 x n_wms>
        %   n_rp  - число ДП
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
        
        % Преобразование индекса метрики в индекс ИХРов
        function [met_id] = rp2met(varargin)
        %met2icr Преобразует индекс вектора метрик в индексы ДП
        %   1 - число ДП
        %   2 - 1й индекс ДП
        %   [3 - 2й индекс ДП]
        %   Метрики занесены в вектор, представляющий из себя верхнюю 
        %   треугольную матрицу метрик, вытянутую по строкам.
            % Если передан инд-с одного ДП, то выдать инд-ы всех соотв. метрик
            if (nargin==2)
                n_rp = uint32(varargin{1});     % преобразования к целым числам
                i = uint32(varargin{2});        % преобразования к целым числам
                if i>1 % если первый ДП, то колонки не будет
                    k_c = zeros(1,i-1,'uint32');    % задать вектор колонки
                    for j=1:i-1
                        k_c(j) = n_rp * j - sum(1:j)+i-n_rp;    % расчет индексов метрик в колонке
                    end
                else
                    k_c = []; % если
                end
                if i<n_rp
                    k_r = (i*n_rp-n_rp+i+1-sum(1:i)):1:(n_rp*i-sum(1:i));
                else 
                    k_r = [];
                end
                met_id = [k_c k_r];
            % Если передан инд-с двух ДП, то выдать инд-с метрики меж. ними
            elseif nargin ==3
                n_rp = varargin{1};
                i = max([varargin{2} varargin{3}]); 
                j = min([varargin{2} varargin{3}]); 
                met_id = n_rp * j - sum(1:j)+i-n_rp; % расчет инд-са метрики
            % Иначе число аргументов не соответствует заданному
            else 
                error('Argument input error')
            end
            
        end
        
        % Итерационное исключение метрик
        function [metrics_hist] = reduce_ms(wms,thr_ms)
        %REDUCE_MS - считает число метрик, превышающих пороги в векторе
        %thr_ms
        %   wms          - вектор с метриками
        %   thr_ms       - вектор с пороговыми значениями
        %   metrics_hist - вектор с числом метрик, превышающих пороговые значения
            wms = sort(wms);            % упорядочить метрики по возрастанию
            n_wms = length(wms);        % число метрик
            n_thr_ms = length(thr_ms);  % число пороговых значений
            metrics_hist = zeros(size(thr_ms,1),size(thr_ms,2)); 
            for i = 1:n_thr_ms
                ind = find(wms>thr_ms(i),1);            % найти индекс первой метрики большей порога = число метрик меньше порога
                if isempty(ind); metrics_hist(i) = 0;   % если нет, то 0 метрик больше порога
                else; metrics_hist(i) = n_wms-ind+1;    % иначе число метрик больше порога равно общему числу минус номер первой метрики большей порога
                end
            end
        end
        
        % Итерационное прореживаение базы данных
        function [en_rp] = reduce_rpdb(wms,thr_ms)
        %REDUCE_RPDB Удаление ДП по величине метрик
            wms = single(wms);
            thr_ms = single(thr_ms);

            n_thr = length(thr_ms);     % число порогов/итераций
            n_wms = length(wms);        % число метрик
            n_rp = (1+sqrt(1+8*n_wms))/2; % число дп

            [~,tr_s_wms] = sort(wms);   % матрица пререхода к сортированному списку
            tr_s_wms = uint32(tr_s_wms);
            tr_wms2rp = Slozhno.met2rp(n_rp);       % матрица перехода от индексов вектора метрик к индексам ДП
            en_rp = ones(n_thr, n_rp,'logical');    % создать вектор наличия ДП (матрица)
            en_wms = ones(1,n_wms,'logical');       % матрица наличия метрик
            i_thr = 1;                  % счетчик исследуемого диапазона метрик
            i_left_rp = n_rp;
            for i_s_wms=1:n_wms         % шаги по вектору сортированных индексов матрицы метрик
                i_wms = tr_s_wms(i_s_wms);          % получение индекса метрики
                if en_wms(i_wms)                    % не исключена ли метрика ранее
                    while wms(i_wms)>thr_ms(i_thr)  % не вышла ли метрика за тек. диапазон
                        i_thr=i_thr+1;              % если да, перейти к след. диапазону
                        if i_thr>n_thr; return; end % выйти из подпрограммы если достигнут конец
                    end
                    i_rp = tr_wms2rp(i_wms,:);  % получение двух ДП между кот. расч. метрика
                    if i_left_rp==2             % если последния метрика
                        en_rp(i_thr:end, i_rp(2))=false;    % занести нули в маску ДП
                        break;
                    end
                    min_met = zeros(1,2);       % буфер под мин. метрику
                    i_ex_wms = cell(2,1);       % буфер под индексы исключ метрик
                    for i = 1:2                 
                        tmp = Slozhno.rp2met(n_rp, i_rp(i));% получение индексов всех метрик этого элемента
                        i_ex_wms{i} = tmp(en_wms(tmp));     % искл-ие инд-ов ранее искл-х метрик
                        tmp = sort(wms(i_ex_wms{i}));
                        min_met(i) = tmp(2);    % запомнить минимальную метрику
                    end
                    if min_met(1)<min_met(2)    % если мин. метрика 1ого ДП меньше
                        en_wms(i_ex_wms{1})=false;          % занести нули в маску метрик
                        en_rp(i_thr:end, i_rp(1))=false;    % занести нули в маску ДП
                    else
                        en_wms(i_ex_wms{2})=false;          % занести нули в маску метрик
                        en_rp(i_thr:end, i_rp(2))=false;    % занести нули в маску ДП
                    end
                    i_left_rp = i_left_rp-1;                % уменьшить число ост. ДП
                end
            end 
        end
        
        % фильтр баттерворта n-ого порядка
        function signal = bworth_filter(signal,t_step_ns,filt_band_mhz, varargin)
        %filter_2order - Функция фильтрации сигнала фильтром с АЧХ Баттерфорта
        %Аргументы:   
        %   signal          - выборки сигнала должны иметь четное число отсчетов
        %   t_step_ns       - шаг дискретизации, нс
        %   filt_band_mhz   - граничная частота по -3дБ, МГц
        %   [n_filter]      - порядок фильтра
            if nargin == 4,  n_filter = varargin{1};    % задать порядок, указанный на входе
            else,            n_filter = 2;              % по умолчанию фильтр 2-го порядка
            end
            f_sampling_mhz = 1e3/t_step_ns;             % частота дискретизации
            w_n = filt_band_mhz/(f_sampling_mhz/2);                 % нормированная частота среза
            [b,a] = butter(n_filter, w_n, 'low');  % формирование коэффициентов фильтра Баттерворта
            signal = filter(b,a,signal.').';         % фильтрация сигнала
        end
        
        % Стандартная установка осей графиков
        function setup_axis_labels( x_label,y_label )
        %SETUP_AXIS_LABEL( xLabelH,yLabelH ) Программа репозиционирующая подписи осей
        % должным образом
        %   xLabelH - идентификатор подписи оси x - получается на вых. xlabel(_)
        %   yLabelH - идентификатор подписи оси y - получается на вых. ylabel(_)
        %   Программа переносит подпись оси x в правый нижний угол, 
        %   подпись оси y - в верхний левый угол и поворачивает в гориз. положение
        
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

        % Оценка минимальной метрики
        function cur_rho = buryi(fvs_red)
            n_obj = length(fvs_red);
            cur_rho = ones(n_obj-1, n_obj);   % контейнер метрик
            for i_obj = 1:n_obj-1
                for j_obj = (i_obj+1):n_obj
                    cur_rho(i_obj,j_obj) = min(Slozhno.calc_b_metrics(...
                        fvs_red{i_obj},fvs_red{j_obj}),[],'all');  % расчет минимальной метрики для i-j сочетания
                end
            end
            cur_rho = min(cur_rho,[],'all');   % поиск
        end
        
        % Оценка по n-й минимальной метрике
        function cur_rho = nmin_metric(fvs_red, n_min_ratio)
        %nmin_metric выполнить оценку устойчивости пространства признаков по заданной вероятности
        %   fvs_red     - массив ячеек с векторами признаков 
        %   n_min_ratio - доля метрик, по которым считается величина метрики (<1)
        %   w_rep       - коэффициенты репликации
        %   TODO: предусмотреть карту вероятностей отдельных ракурсов
            n_obj = length(fvs_red);
            n_asps = cellfun(@(fvs) size(fvs,1), fvs_red);          % число ракурсов в каждом
            mat_ij = repmat((1:n_obj),n_obj,1);                 % 
            asp_i = mat_ij(triu(ones(n_obj,'logical'),1).');    % первый элемент сочитания
            mat_ij = mat_ij.';                                  %
            asp_j = mat_ij(triu(ones(n_obj,'logical'),1).');    % второй элемент сочетания
            n_met_total = sum(n_asps(asp_i).*n_asps(asp_j));     % всего ДП
            % множитель-дополнитель
            elem_rho = prod(n_asps, 'all');
            m_ij = arrayfun(@(i,j) elem_rho/(n_asps(i)*n_asps(j)), ...
                asp_i, asp_j);              % множители репликации
            m_ij = floor(m_ij/min(m_ij));
            % общ число метрик * доля допустимых ощибок * число на которое будем реплицировать
            cur_rho = cell(length(m_ij),1);
            n_met_thresh = ceil((n_met_total * n_min_ratio)./m_ij);  % 
            for i_case = 1:length(m_ij)
                cur_rho{i_case} = sort( reshape(...
                    Slozhno.calc_b_metrics(fvs_red{asp_i(i_case)}, fvs_red{asp_j(i_case)}), ...
                    1,[]), 'ascend');  % расчет метрик для сочет., перевод в строку и сортировка

                cur_rho{i_case} = repmat( ...
                    cur_rho{i_case}(1:n_met_thresh(i_case))... % редукция до метрики n/ коэф. репликации
                    , 1, m_ij(i_case));          % репликация (вырезать nthr/коэф репликации)
                                 
            end
            cur_rho = sort(horzcat(cur_rho{:}),'ascend');   % сортировка общей матрицы
            cur_rho = cur_rho(ceil(length(cur_rho)*n_min_ratio)); 
        end
        
        % Оценка по n-й минимальной метрике при сбалансированной выборке
        function cur_rho = nmin_metric_balanced(fvs_red, n_min_ratio)
        %nmin_metric выполнить оценку устойчивости пространства признаков по заданной вероятности
        %   fvs_red     - массив ячеек с векторами признаков 
        %   n_min_ratio - доля метрик, по которым считается величина метрики (<1)
        %   w_rep       - коэффициенты репликации
            n_obj = length(fvs_red);
            n_asps = cellfun(@(fvs) size(fvs,1), fvs_red);      % число ракурсов в каждом
            mat_ij = repmat((1:n_obj),n_obj,1);                 % 
            asp_i = mat_ij(triu(ones(n_obj,'logical'),1).');    % первый элемент сочитания
            mat_ij = mat_ij.';                                  %
            asp_j = mat_ij(triu(ones(n_obj,'logical'),1).');    % второй элемент сочетания
            n_met_total = sum(n_asps(asp_i).*n_asps(asp_j));    % всего метрик
            % множитель-дополнитель
            elem_rho = prod(n_asps, 'all');
            m_ij = arrayfun(@(i,j) elem_rho/(n_asps(i)*n_asps(j)), ...
                asp_i, asp_j);              % множители репликации
            m_ij = floor(m_ij/min(m_ij));
            % общ число метрик * доля допустимых ощибок * число на которое будем реплицировать
            cur_rho = cell(length(m_ij),1);
            n_met_thresh = ceil((n_met_total * n_min_ratio)./m_ij);  % 
            for i_case = 1:length(m_ij)
                cur_rho{i_case} = sort( reshape(...
                    Slozhno.calc_b_metrics(fvs_red{asp_i(i_case)}, fvs_red{asp_j(i_case)}), ...
                    1,[]), 'ascend');  % расчет метрик для сочет., перевод в строку и сортировка

                cur_rho{i_case} = repmat( ...
                    cur_rho{i_case}(1:n_met_thresh(i_case))... % редукция до метрики n/ коэф. репликации
                    , 1, m_ij(i_case));          % репликация (вырезать nthr/коэф репликации)
                                 
            end
            cur_rho = sort(horzcat(cur_rho{:}),'ascend');   % сортировка общей матрицы
            cur_rho = cur_rho(ceil(length(cur_rho)*n_min_ratio)); 
        end
        
        % Оценка пространства признаков по методу гистограмм
        function cur_rho = hist_acc(fvs_red, hist_bins, n_metrics_fraction)
        % расчет гистограммы метрик
                    hist_counts = Slozhno.cpu_linear_histcounts(fvs_red,hist_bins);         % гистограмма
                    hist_counts = vertcat(hist_counts{triu(ones(length(fvs_red),'logical'),1)});  % вытащить значения по верхней треугольной матрице без диагонали
                    hist_counts = sum(hist_counts,1);                   % просуммировать
                    acc_hist = cumsum(hist_counts);                     % интеграл от гист
                    thrld = acc_hist(end)*n_metrics_fraction;           % порог по вероятности
                    cur_rho = hist_bins(find(acc_hist>thrld,1));        % значение метрики, соответствующей порогу для тек. пространства признаков
        end
        
        % редукция пространства признаков
        function [fspace_rhos, fspace_map, varargout] = fs_reduction(fvs, varargin)
        %fs_reduction редукция с учетом размерности выборки
        % на выходе
        % fspace_rhos - вектор-столбец с метриками
        % fspace_maps - матрица по строкам которой расположены настройки пространства
        % Аргументы: 
        % fvs - вектора признаков (по строкам)
        % Именованные параметры:
        % dimensions, reduction_type, n_metrics_fraction, hist_edges,
        % obj_weight
        warning('Набор функций устарел. Настоятельно рекомендуется использовать класс RecursiveReduction.')
            dim_fv = size(fvs{1},2);            % размерность ВП
            n_obj = length(fvs);
            kwargs = KeywordArguments(...
                'dimensions',1:dim_fv, ...          % исследуемые размерности
                'reduction_type', 'none',...        % исследуемые размерности
                'n_metrics_fraction', 0.05,...      % число метрик
                'hist_edges', linspace(0,2,51),...  % границы интервалов гистограммы
                'obj_weights',ones(n_obj,1)/n_obj,... % веса объектов (для учета априорной вер-сти)
                'n_nearest',5, ...          % число проверяемых соседей
                'm_technique','mean',...    % метод расчет
                'k_alien',2 ...               % число своих
                );
            [dimensions, reduction_type, n_metrics_fraction, hist_edges, obj_weights, ...
                n_nearest, m_technique, k_alien] =  ...
                kwargs.parse_input_cell(varargin);
            % настройка специфическая для метода
            switch reduction_type
                case 'nmin'
                    if range(arrayfun(@(i) size(fvs{i},1),1:n_obj))==0  % если все одинаковые, то считаем выборку сбалансированной
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
                    error('Неправильно задан метод редукции. Поддерживаемые: nmin, buryi, mhist, minalien')
            end
            dimensions = dimensions(dimensions<=dim_fv);
            n_spaces = length(dimensions);
            
            fspace_rhos = zeros(n_spaces,1);        % значения пороговой метрики для размерностей
            fspace_map = zeros(n_spaces, dim_fv);   % матрица оптимальных пространств - по строкам от размерности
            fprintf('\tВозможных пространств: %d.',2^dim_fv);
            fprintf('\tВыполнено: 00.00%%');
            for feat_val = 1:2^dim_fv-1             % по всем возможны пространствам признаков
                cur_feat_map = de2bi(feat_val,dim_fv,'left-msb'); % получение карты текущего пространства
                i_fmap = find(sum(cur_feat_map)==dimensions);  % индекс соответствующий размерности
                if i_fmap
                    cur_feat_map = logical(cur_feat_map);
                    fvs_red = cellfun(@(fvs) fvs(:,cur_feat_map), fvs, 'UniformOutput',false); % редукция
                    
                    [cur_rho] = rho_estimate(fvs_red);
                    
                    if cur_rho>=fspace_rhos(i_fmap)     % если значение метрики больще значения для пр. признаков такой размерности
                        fspace_rhos(i_fmap) = cur_rho;
                        fspace_map(i_fmap,:) = cur_feat_map(:);
                    end
                    fprintf('\b\b\b\b\b\b%05.2f%%',feat_val/2^dim_fv*99.99);
                end
            end
            fprintf('\b\b\b\b\b\b\b 100%%\n');
            if nargout==3, varargout{1} = dimensions; end
                        
        end
        
        % гистограмма метрик
        function [bms_counts] = cpu_linear_histcounts(fvs, hist_bins) 
        %cpu_histcounts
        % 
            n_obj = length(fvs);
            bms_mod = cell(n_obj, n_obj);       % матрица под метрики
            bms_counts = cell(n_obj, n_obj);    % матрица отсчетов внешних метрик
            for i_obj = 1:(n_obj-1)             % расчет метрик
                for j_obj = (i_obj+1):n_obj
                    bms_mod{i_obj,j_obj} = Slozhno.calc_b_metrics(fvs{i_obj},fvs{j_obj}); % расчет угловых метрик
                    bms_counts{i_obj,j_obj} = ...
                        histcounts(bms_mod{i_obj,j_obj}, hist_bins );
                end
            end
        end
        
        % Функция расчета величины разделимости множества tsne
        function q_sep = tsne_separability(data, varargin)
        % загрузить информацию - данные tsne и целевые вектора    
            kwargs = KeywordArguments(...
                'n_nearest',false, ...          % число проверяемых соседей
                'm_technique','mean',...    % метод расчет
                'k_alien',3 ...               % число своих
                );
            [n_nearest, m_technique, k_own] =  ...
                kwargs.parse_input_cell(varargin);
            if ischar(data)
                in_data = load(data,'tsne_data');   % загрузить информацию
                features = vertcat(in_data.tsne_data.tsne_scatter);
                targets = {in_data.tsne_data.name};
                if n_nearest
                else, n_nearest = in_data.perplex;  % прочитать число ближ. сосед. алг-ма tsne
                end
            elseif(iscell(data))
                if nargin==3
                    features = vertcat(data{:});
                    targets = varargin{1};
                else
                    error('Missing argument');
                end
            else
                error('Ошибка в формате аргументов')
            end
            q_sep = Slozhno.nn_separability(features, targets, ...
                'n_nearest',n_nearest, 'm_technique',m_technique, 'k_alien', k_own);
        end
        
        % Функция расчета разделимости по векторам признаков
        function q_sep = nn_separability(features, targets_str, varargin)
        % загрузить информацию - данные вектора признаков и целевые вектора    
            kwargs = KeywordArguments(...
                'n_nearest',5, ...          % число проверяемых соседей
                'm_technique','mean',...    % метод расчет
                'k_alien',2 ...             % число своих
                );
            [n_nearest, m_technique, k_alien] =  ...
                kwargs.parse_input_cell(varargin);
            n_samples = size(features, 1);
            rhos = Slozhno.calc_b_metrics(features, features);  % расчитать метрики
            [~, sort_indexes] = mink(rhos, n_nearest+1, 2);                  % сортировать
            clear('rhos')
            [~, ~, targets] = unique(targets_str);        % преобразовать имена в id-номера
            
            switch m_technique
                case {'mean' 'thresh'}
                    neq_matrix = zeros(n_samples, n_nearest); % предсоздание м-цы эквив-сти
                    for i=1:n_nearest
                        neq_matrix(:, i) = targets ~= targets(sort_indexes(:, i+1));
                    end
                    switch m_technique
                        case 'mean'     % среднее число точек другого класса среди соседей
                            q_sep = mean(sum(neq_matrix, 2)/n_nearest);
                        case 'thresh'   % доля объектов для которых число ближайших соседей больше порога
                            q_sep = sum((sum(neq_matrix, 2)>(k_alien)))/n_samples;
                    end
                    
                case 'prat'
                    local_groups = targets(sort_indexes(:,1:n_nearest));   % № объекта -> классы
                    group_modes = mode(local_groups,2);              % превалирующий класс в группе
                    q_sep = 1-mean(sum(group_modes==local_groups,2))/n_nearest; % среднее число привал.
                    
            end
            
        end
        
        % Расчет махимизируемого числа своих
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
            rhos = Slozhno.calc_b_metrics(fvs, fvs);  % расчитать метрики
            [~, sort_indexes] = sort(rhos, 2);                  % сортировать
            clear('rhos')
            neq_matrix = zeros(n_samples, n_nearest); % предсоздание м-цы эквив-сти
%             [~, ~, targets] = unique(targets_str);        % преобразовать имена в id-номера
            for i=1:n_nearest
                 neq_matrix(:, i) = targets ~= targets(sort_indexes(:, i+1)); % i+1 -учет того, что 1й сосед - сам ветор
            end
            switch m_technique
                case 'mean'     % среднее число точек другого класса среди соседей
                    q_sep = mean(sum(neq_matrix, 2)/n_nearest);
                case 'thresh'   % доля объектов для которых число ближайших соседей больше порога
                    q_sep = sum((sum(neq_matrix, 2)>(k_alien)))/n_samples;
            end
            cur_rho = 1 - q_sep;
        end
        
    end
        
end

