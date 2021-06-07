classdef RPTools
    %RPTools Ver.01.02 Здесь хранятся функции для  выполнения эксперимента
    %   Все функции выполнены статически
    
    % Основные функции
    methods(Access = public, Static = true)
        
        % Сформировать ДП
        function [rps] = make_rps(icrs, gauss_hwidth)
            % гауссовская функция по трем СКО с шириной gauss_width по
            % полувысоте
            hw2sig = 2.35;              % отношение ширины по половине к сигма
            sig = gauss_hwidth/hw2sig;  % отсчетов в одной сигма
            t_axis = linspace(-4,4, sig*8);
            imp = exp(-t_axis.^2/2);    % импульс ЗИ
            rps = conv2(icrs, imp);     % свертка ИХР с огибающей импульса
            if mod(size(rps,2),2)
                rps = rps(:,1:end-1);   % если нечетное число отсчетов, то вырезать последний
            end
        end
        
        % Аддитивный шум
        function [rp_samples, meta] = add_noise(rps, meta, sigma, n_real, multiplier)
        % Сформировать множество реализаций с аддитивным гауссовским шумом
        %   Ввод: 
        %   rps – матрица, содержащая отсчеты ДП
        %   meta – cell array с наименованиями объектов и углами ракурсов ? и ?
        %   sigma – значение СКО для БГШ
        %   n_real – число реализаций
        %   Вывод:
        %   rps – матрица, содержащая множество реализаций ДП
        %   meta – cell array, содержащий информацию для обновленного массива реализаций ДП
            signal_power = mean(rps.^2,2);  % мощность ДП на интервале, соответ. 3м метрам
            [n_rps, dim_rp] = size(rps);
            total_samples = n_rps*n_real;   % всего реализаций
            fprintf('\tОтн. С/Ш: ')
            if length(sigma)==1             % если задано одно значение с/ш
                signal_power = num2cell(signal_power/sigma^2);  % с/ш = ср.мощн.ДП(3метра)/ср.мощн.шума(3метра)
            else                            % если задан набор отношений с/ш
                %TODO: Добавить более точный расчет С/Ш
                signal_power = num2cell(signal_power/max(sigma)^2); % считать худшее с/ш
            end
            if ~isempty(meta)
                [meta(:).snr] = deal(signal_power{:});
                [~,i_min] = min([meta.snr]); [~,i_max] = max([meta.snr]);
                fprintf('от %01.1f(%s(%+03d,%+03d)\x00B0) до %01.1f(%s(%+03d,%+03d)\x00B0)\n',...
                    meta(i_min).snr,meta(i_min).name, meta(i_min).asp_a, meta(i_min).asp_b,...
                    meta(i_max).snr,meta(i_max).name, meta(i_max).asp_a, meta(i_max).asp_b);
            else, warning('No metadata')
            end            
            % временные сдвиги
            max_shift = dim_rp*(multiplier-1);      % максимальный временной сдвиг
            if(multiplier>1)
                fprintf('\tОкно случ сдвига: %d длин исх\n', multiplier);
                r_start  = 1 + ...
                    (floor(max_shift * rand(total_samples,1))); % формирование матрицы со значениями случайных сдвигов 0.1*ones(total_samples,1))); 
            else
                r_start = ones(total_samples,1);    % если неоднозначности
            end
            r_stop   = r_start + dim_rp - 1;        % конечные индексы
            sigma_mix = (max(sigma)-min(sigma))*rand(total_samples,1,'like', rps)+...
                min(sigma); % дисперсии равномерно распределены от сигма мин до сигма макс
                        
            rp_samples = randn( total_samples, dim_rp*multiplier, 'like', rps );    % шумовые в выборки
            rp_samples = rp_samples .* sigma_mix;
            % Прибавляем к шумовым выборкам случайно сдвинутый портрет
            for i = 1 : total_samples               % цикл по всем реализациям
                rp_samples(i, r_start(i):r_stop(i)) = ...
                        rp_samples( i , r_start(i):r_stop(i) ) + ...
                        rps(mod(i-1, n_rps)+1,:);
            end
            meta = repmat(meta,n_real,1);
            
        end
        
        % Дискретизация
        function rps = decimate(rps, frac)
            rps = movmean(rps, frac, 2);
            rps = rps(:, 1:frac:end);
        end
                
        % Квантование АЦП
        function dig_rp = adc(rp,varargin)
        %adc выполняет преобразование передаваемых данных в множество
            if nargin == 1      % если битность не передана, то выполняется по умолчанию
                dtype = 'uint8';    % по умолчанию 8 бит
                nMax = 2^8;         % по умолчанию 8 бит
            elseif nargin == 2  % если передана битность
                bit = floor(varargin{1});
                nMax = 2^bit;   % максимальное допустимое значение
                if      and(bit >3, bit<=8);    dtype = 'uint8' ;   % 8 бит
                elseif  and(bit > 8, bit<=16);  dtype = 'uint16';   % 16 бит
                %elseif  and(bit > 16, bit<=32); dtype = 'uint32';   % 32 бита
                else; error('Битность АЦП должна быть в пределах от 4 до 16');
                end
            else 
                error('Неверное число входных аргументов')
            end
            rp_peak_value = max(rp,[],'all'); 
            if rp_peak_value>1
                warning('Динамический интервал АЦП: Отсечка.'); 
                rp(rp > rp_peak_value) = 1; % отсекаем выходящие
            end
            dig_rp = cast(floor(rp * nMax),dtype);
        end
        
        % Нормировка
        function rps = norm_rp( rps, norm_mode )
            [n_rps, ~] = size(rps);
            switch norm_mode
                case 'mean'     % нормировка по среднему значению
                    norm_factors = 1./mean(rps,2);
                    norm_factors = norm_factors/max(norm_factors);
                case 'max_peak' % нормировка по глоб мамксимуму ДП
                    norm_factors = 1./max(rps,[],2)*.7; % при нормировке берется коэффициент запаса 0.7
                case 'energy'
                    norm_factors = 1./sqrt((sum(rps.^2, 2))); % рассчитываем энергию с нормированием на N                    
                case 'mod'
                    norm_factors = 1./sqrt((sum(rps.^2, 2))); % рассчитываем энергию с нормированием на N                    
                otherwise
                    if ~strcmp(norm_mode,'none') % вывести предупреждение если явно не указано
                        warning('Нормировка не выполняется')
                    end
                    return
            end
            norm_factors(norm_factors==inf) = 0;
            for i = 1:n_rps
                rps(i,:) = rps(i,:)*norm_factors(i);
            end
                        
        end
        
        % Функция поиcка центра тяжести
        function ind = find_center(randrp)
        %FINDCENTER     Поиск центра масс в массиве по строкам
            [nrp, nDim] = size(randrp);
            halfMass = sum(randrp,2)/2;         % найти половину 
            ind = zeros(nrp,1,'uint32');        % создать матрицу индексов
            sbuf = zeros(nrp,1);                % создать матрицу значений сумм
            for i = 1:nDim 
                sbuf = sbuf + cast(randrp(:,i),'like',sbuf);      % вычислить промежуточное значение
                mask = sbuf <= halfMass;
                ind = ind + cast(mask,'uint32');% прибавить индекс
            end
        end
        
        % Устранить временную неопределенность
        function rps = rp_time_snap( rand_rps, snap_mode, nc_offset )
        %RP_TIME_SNAP
        %   Устранить временную неоднозначность
        %   rand_rps    - матрица ДП по строкам
        %   snap_mode   - способ устранения временной неоднозначности: 
        %       energy      - центр тяжести
        %       max_peak    - по пику
        %   nc_offset   - размер исходного ДП. Определяет какой участок ДП
        %   в каждую сторону от точки привязки будет взят.
            [n_rps, ~] = size(rand_rps);
            % Поиск характерной точки ДП для привязки
            switch snap_mode
                case 'energy'       % характерная точка - центр тяжести
                    snap_index = RPTools.find_center(rand_rps);
                case 'max_peak'     % характерная точка - глобальный максимум
                    [~, snap_index] = max(rand_rps,[],2);
                case 'none'         % ДП используются как есть
                    rps = rand_rps;
                    return
                otherwise
                    error(['Метод ' snap_mode ' не поддерживается']);
            end
            
            rps = zeros(n_rps, 2*nc_offset, 'like',rand_rps);   % подготовить матрицу под ДП
            rand_rps = [...
                zeros(n_rps,nc_offset,'like',rand_rps), ...
                rand_rps, ...
                zeros(n_rps,nc_offset,'like',rand_rps)];        % дополнить выборку на случай выхода за пределы
            start = snap_index + 1;         %
            stop = snap_index + 2*nc_offset;
            for i = 1:n_rps
                rps(i,:) = rand_rps( i, start(i):stop(i) );
            end
            
        end

        % редукция множества признаков КВП
        function [ features, coef_index ] = cwt_reduce( features_src, coef_map ,n_features)
        %CWT_REDUCE
        % Выполняется отбор признаков, при том отбор начинается из центра
        %TODO: проверить соответствует ли схема отбора признаков матрице,
        %формируемой при КВП
            [n_rps, ~] = size(features_src);
            features = zeros(n_rps, n_features);
            coef_index = zeros(2 ,n_features); % записывается уровень и номер коэффициента
            
            n_coef = coef_map(end); % установить число коэффициентов
            cntr = 100;
            i_lev = length(coef_map)+1;
            for i = 1:n_features
                if cntr > n_coef
                    i_lev = i_lev - 1;          % следующий уровень
                    n_coef = coef_map(i_lev);   % размерность уровня разложения
                    cntr = 1;                   % нулевая отстройка
                    lev_premult = + 1 - 2 * mod( n_coef, 2 ); % множетель зависящий от четности числа коэффициентов
                    % odd = + ; even = - 
                end
                offset = floor(cntr / 2);
                sign_premult = 1 - 2 * mod( cntr, 2 ); % начинается с + 0

                i_coef =  ceil(n_coef/2) +...
                   lev_premult * sign_premult * offset;  % округляем на случай если acwt (5/2 = 2.5 -> 3 -> 2 -> 4...)
                ind = sum( coef_map(1:i_lev) ) - coef_map(i_lev)  + i_coef;
                features(:,i) = features_src(:,ind);
                coef_index(:,i) = [i_lev; i_coef];
                cntr = cntr + 1;     % увеличить счетчик
            end
        end
        
        % Формирование реализаций ДП со случайным сдвигом и аддитивным шумом
        function samples = generate_samples(rp,nSamples,nExtraDim,dSigm)
        %FRANDOMSHIFT Формирование случайного сдвига ДП
        %   Сформировать из переданных ДП на основе переданных параметров
        %   шумовые выборки
        %   Args: 
        %       rp - матрица, в которой ДП расположены по строкам.
        %       nSamples - целое число реализаций на каждый исходный ДП.
        %       nExtraDim - во сколько раз увеличивается интервал наблюдения.
        %       dSigm - СКО БГШ.
            [nIcrs, dimrp] = size(rp);      % число ДП и их размерность
            totalSamples = nIcrs*nSamples;  % всего реализаций
            
            samples = dSigm * ... 
                randn( totalSamples, dimrp*nExtraDim, 'like',rp );    % шумовые в выборки
            maxTShift = dimrp*(nExtraDim-1);% максимальный временной сдвиг
            rStart  = 1 + (floor(maxTShift * rand(totalSamples,1)));% формирование матрицы со значениями случайных сдвигов
            rStop   = rStart + dimrp - 1;   % конечные индексы
            i = 0;
            for iIcr = 1 : nIcrs            % цикл по всем ракурсам (ихр)
                for iSample = 1 : nSamples  % цикл по всем реализациям
                    i = i + 1;
                    samples( i , rStart(i):rStop(i) ) =...
                        samples( i , rStart(i):rStop(i) ) + rp(iIcr,:);
                end
            end
        end
    
    % Интерполяция 
        
        % Управляющая функция интерполяции
        function signal = interpol(signal, fs, fq, interp_method)
        %interol Изменить частоту дискретизации с использованием интерполяции
        %   Функция осуществляет переход от исходной частоты дискретизаци fs 
        %   к целевой частоте fq по заданному методу интерполяции
        %   Args:
        %       signal - матрица, в которой ДП расположены по строкам.
        %       fs - исходная частота дискретизации сигнала.
        %       fq - целевая частота дискретизации сигнала.
        %       interp_method - метод интерполяции. Поддерживаемые методы:
        %           none     не выполнять изменение частоты;
        %           freq     частотная интерполяция;
        %           poly     полиномиальная интерполяция;
        %           Lagrange Лагранжева интерполяция;
        %           spline   интерполяция сплайнами.
            switch interp_method
                case 'none'        
                case 'freq',        signal = RPTools.freq_interpol(signal, fs, fq);     % частотная интерполяция
                case 'poly',        signal = RPTools.poly_interpol(signal, fs, fq);     % полиномиальная интерполяция
                case 'Lagrange',    signal = RPTools.lgrnge_interpol(signal, fs, fq);   % Лагранжева интерполяция
                case 'spline',      signal = RPTools.spline_interpol(signal, fs, fq);   % интерполяция сплайнами
                otherwise, warning('Неизвестный метод интерполяции! Интерполяция не выполнена')
            end
        end
        
        % Интерполяция полиномами
        function signal = poly_interpol(signal, fs, fq)
            t_axis = 0:1/fs:1/fs*(size(signal,2)-1);
            t_axis_est = 0:1/fq:t_axis(end);
            signal = interp1(t_axis,signal.',t_axis_est, 'pchip').'; % интерполяция полин. 3й степ
            
        end
        
        % Интерполяция с использованием частотного метода
        function signal = freq_interpol(signal, fs, fq, varargin)
        % 
            if nargin>3, err = varargin{1}; else, err = 0.01; end
            [p, q]=rat(fq/fs, err);                 % получение числ и знам рац дроби  
            signal = resample(double(signal.'), p,q).';
        end
        
        % Интерполяция сплайнами
        function signal = spline_interpol(signal, fs, fq)

            t_axis = 0:1/fs:1/fs*(size(signal,2)-1);
            t_axis_est = 0:1/fq:t_axis(end);
            signal = spline(t_axis, signal, t_axis_est);
        end
        
        % Лагранжева интерполяция
        function signal = lgrnge_interpol(signal, fs, fq, varargin)
            if nargin>3,    n_lagrange = varargin{1};   % порядок полиномов
            else,           n_lagrange = 5;             % порядок полиномов
            end

            [int_factor, q]=rat(fq/fs, 0.01);         % получение числ и знам рац дроби  
            lgr_filter = intfilt(int_factor,n_lagrange,'Lagrange');   % фильтр д/интерполяции

            signal = upsample(signal.',int_factor);            % добавление нулей
            signal = filter(lgr_filter,1,signal);       % применение фильтра
            signal(1:floor(mean(grpdelay(lgr_filter))),:) = []; % срезать отсчеты задержки
            signal = signal(1:q:end,:).';
        end
        
    % Работа с файлами
        
        % Сохранить csv
        function [] = save_csv( features, meta, f_header, filename )
        %SAVE_CSV Сохраняет в csv признаки
        %   Ввод:
        %   features    - матрица признаков
        %   meta        - строковые обозначения элементов (имя объекта [ракурс])
        %   f_header    - признаков
        %   filename    - имя файла
        %   В первой строке содрежатся названия признаков, а t обозначает
        %   объект. Так, для признаков 
                       
            [ n_obj, ~ ] = size(features);
            %meta = RPTools.meta2str(meta); % восстановить если понадобится вывод с
            %ракурсами
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
            
            disp(['Файл ' filename ' сохранен'])
            
        end
            
        % Открыть файлы ИХР
        function [icrs, meta] = open_icr_file(varargin)
        %open_icr_file Открывает файл, или несколько файлов
        %   Ввод:
        %   filename – строка или cell array строк (опциональный)
        %   Вывод:
        %   icrs – матрица размерности n_objects x n_counts, где n_objects – число объектов
            obj_names = {'sph' 'con' 'cyl' 'brk' 'dsk' 'plt'};
            % ----- обработка вариантов ввода 
            if nargin == 0      % пустой ввод
                filelist = RPTools.fGetFileList('icr').';
            else 
                if iscell(varargin)
                    filelist = varargin{1};
                else; error('Можно передавать только cellarray')
                end
            end
            % ----- открытие файлов
            nFiles = size(filelist,1);     % прочитать число файлов
            meta = cell(nFiles,2);  % массив с метаинформацией
            icrs = cell(nFiles,1);  % массив с отсчетами ИХР
            for i_file = 1:nFiles
                fileID = fopen(filelist{i_file});   % присвоить файлу внутренний ID
                icrs{i_file} = textscan(fileID, '%f32');          % считать из файла отсчеты ИХР(float32)
                icrs{i_file} = icrs{i_file}{1}.';
                textscan(fileID, '%s',6);           % считывание заголовка метаинформации
                f_meta = textscan(fileID, '%u %f32 %f32 %f32 %f32 %f32');   % прочитать метаинформацию
                fclose(fileID);                         % отвязать ID
                % 1) N_obj; 2) H_obj; 3) R_obj; 4) L_elm; 5) Fi_az; 6) Fi_um
                meta{i_file, 1} = [ obj_names{f_meta{1}} '('... 
                    num2str(f_meta{2}) 'x' num2str(f_meta{3}) ')' ];                % вывод имени
                meta{i_file, 2} = f_meta{5}; meta{i_file, 3} = f_meta{6};   % углы ракурса
            end
            icrs = vertcat(icrs{:});
            meta = struct('name', meta(:,1), 'asp_a', meta(:,2),'asp_b', meta(:,3));
            icrs = fillmissing(icrs,'movmean',5,2);     % заполнить NaN значения
            disp(['Открыто ' num2str(nFiles) ' файлов'])
        end
        
        
        % Откртыть файлы ДП
        function [rps, meta] = open_dat_file(varargin)
        %open_dat_file Открыть файл dat с ДП
        %   Ввод:
        %   filename – строка или cell array строк (опциональный)
        %   Вывод:
        %   rps – ДП или множество ДП
        %   n_counts – число отсчетов в ДП
            % ----- обработка вариантов ввода 
            if nargin == 0      % пустой ввод
                filelist = RPTools.fGetFileList('dat').';
            else 
                if iscell(varargin)
                    filelist = varargin{1};
                else; error('Можно передавать только cellarray')
                end
            end
            % ----- открытие файлов
            nFiles = size(filelist,1);     % прочитать число файлов
            meta = cell(nFiles,2);  % массив с метаинформацией
            for i_file = 1:nFiles
                filename = filelist{i_file}; % выделить имя файла
                fileID = fopen(filename);   % присвоить файлу внутренний ID
                rp = textscan(fileID, '%u %f32 %f32');          % считать из файла отсчеты ИХР(float32)
                if (i_file==1); rps = zeros(nFiles, length(rp{2})); end
                rps(i_file,:) = rp{2}.';
                fclose(fileID);                         % отвязать ID
                filename = filename(find(filename=='\', 1, 'last')+1:end);
                meta{i_file, 1} = filename(1:2);   % вывод имени
                meta{i_file, 2} = str2double(filename(4:5)); 
                meta{i_file, 3} = str2double(filename(7:8));   % углы ракурса
            end
            meta = struct('name', meta(:,1), 'asp_a', meta(:,2),'asp_b', meta(:,3));
            disp(['Открыто ' num2str(nFiles) ' файлов'])
        end
        
        % Получить список файлов в cellarray
        function [output_args] = fGetFileList(ext)
        %FGETFILELIST Получить список файлов
        %   Функция предлагает графический интерфейс для выбора файлов и
        %   формирует на выходе список cell с полными именами файлов или
        %   cell с нулем
            [filenames,path] = ...
                uigetfile(...
                ['*.', ext],'MultiSelect','on');  % открыть файлы через ГИП
            if ~iscell(filenames)                       % Открыто <2 файлов
                if (size(filenames,2)==1);filenames={0};% 0 файлов
                else; filenames = {[path,filenames]};   % 1 файл - составить полное имя
                end
            else                                        % несколько файлов
                for i=1:(size(filenames,2))
                    filenames{i} = [path,filenames{i}];
                end
            end
            output_args = filenames; % вывести данные
        end
        
        % Прочитать ИХР + параметры  из файла
        function [ output_args ] = fOpenICR( filename )
        %FOPENICR Загрузить ИХР объекта
        %   Загружает ИХР объекта из файла filename
        %   output_args{1:7} параметры в формате: 
        %   1       2       3       4       5       6
        %   N_obj   H_obj   R_obj   L_elm   N_cnts  [Fi_az Fi_um] 
        %   7
        %   icrValues
            fileID = fopen(filename);               % присвоить файлу внутренний ID
            icrCell = textscan(fileID, '%f32');     % считать из файла массив числе float32
            textscan(fileID, '%s',6); 
            Parameters = textscan(fileID, '%u %f32 %f32 %f32 %f32 %f32');    % прочитать конец файла
            fclose(fileID);                         % отвязать ID
            output_args(1:4) = Parameters(1:4);     % вывод параметров
            output_args{5} = size(icrCell{1,1},1);  % вывод числа отсчетов в ИХР
            output_args{6} = horzcat(Parameters{5:6});     % вывод угла ракурса
            output_args{7} = icrCell{1,1}.';        % вывод отсчетов  
        end
    
    % Интегральные преобразования
        
        % Амплитудное ПФ
        function [ FeatDataF ] = afft( inData, varargin )
        %AFFT Преобразование в пространство признаков Фурье
        %   Производит построчное преобразование в пространство признаков Фурье
        %   На входе матрица inData: каждая строка - ДП
        %   Выделяет гармоники с нулевой до n-й - самой высокочастотной;
        %   если самая высокочастотная компонента 
        %   После вырезания симметричной части спектра домножает компоненты
        %   для которых симметричная составляющая усечена на sqrt2.
            [ ~ , nDim] = size(inData);
%             nDim = 2*nDim;  % удвоение для уравнивания числа отсчётоы ПФ и ВП
            if nargin==2
                nDim = varargin{1};
            end
            FeatDataF = abs(fft(inData,nDim,2))/sqrt(nDim);                 % БПФ вдоль второго измерения
            % Последнюю гармонику исключаем, т.к. fд относительно
            % высока. Нулевую оставляем, т.к. из-за нормировки
            % она несет полезную информацию (отношение пика к среднему)
            FeatDataF = FeatDataF(:,1:(floor(nDim/2)+1));
            FeatDataF(:,2:ceil(nDim/2)) = sqrt(2)*FeatDataF(:,2:ceil(nDim/2)); % домножение на sqrt из-за элиминации гармоник
        end
        
        % Преобразование Фурье (односторонний спектр)
        function [ FeatDataF ] = fft( inData )
        %FFFT Преобразование в пространство признаков Фурье
        %   Производит построчное преобразование в пространство признаков Фурье
        %   На входе матрица inData: каждая строка - ДП
        %   На выход вектор 
        %   1 [0я re гарм]  2 [1я re гарм] 3 [1я im гарм] ... 2n [n-я гарм(действительная)]
            [nObj,nDim] = size(inData);
            if mod(nDim,2); warning('Нечетное число отсчетов'); end
            FeatDataFC = fft(inData, nDim, 2)/sqrt(nDim);           % БПФ вдоль второго измерения
            FeatDataF = zeros(nObj,nDim,'like',inData);
            FeatDataF(:,1) = FeatDataFC(:,1);           % 0-я гармоника
            FeatDataF(:,end) = FeatDataFC(:,nDim/2+1);  % последняя гармоника
            FeatDataF(:,2:2:end-2) = real(FeatDataFC(:,2:nDim/2));  % действительные компоненты
            FeatDataF(:,3:2:end-1) = imag(FeatDataFC(:,2:nDim/2));  % мнимые компоненты
            FeatDataF(:,2:(end-1)) = sqrt(2)*FeatDataF(:,2:(end-1)); % корректировка коэфициентов на sqrt 2
        end
        
        % Комплексное ВП
        function [ features, coef_map ] = cwt( rps, varargin )
        %CWT выполняет КВП
        %   Передается матрица с ДП, можно передать название вейвлета
        %   Карта выходных коэффициентов:
        %   [ВЧ 1 1i 2 2i...] ... [НЧ 1 1i 2 2i... ]
            if nargin == 1
                wave_name = 'dtf2';
            elseif nargin == 2
                wave_name = varargin{1};
            else, error('Слишком много аргументов')
            end
            
            [rps, coef_map] = RPTools.fDTCWT(rps, wave_name);  % сформировать комплексные составляющие
            [n_rps, wt_dim] =size(rps);   % рассчитать размерность вектора признаков cwt
            features = zeros(n_rps, wt_dim);
            features(:,1:2:end) = rps(:,1:wt_dim/2); 
            features(:,2:2:end) = rps(:,wt_dim/2+1:wt_dim);
            coef_map = coef_map * 2;      % корректировка размеров из-за комплексности
        end
        
        % Амплитудное КВП
        function [ FeatDataW, coef_map ] = acwt( InData, varargin )
        %FACWT Формирование амплитудных составляющих комплексного вейвлет
        %преобразования
            if nargin == 1
                wave_name = 'dtf2';
            elseif nargin == 2
                wave_name = varargin{1};
            else, error('Слишком много аргументов')
            end
            [FeatDataWC, coef_map] = RPTools.fDTCWT(InData,wave_name);  % сформировать комплексные составляющие
            nDim =size(FeatDataWC,2);   % рассчитать размерность вектора признаков cwt
            FeatDataW = abs( FeatDataWC(:,1:nDim/2)+...
                1i*FeatDataWC(:,nDim/2+1:nDim));
            %coef_map = flip(coef_map, 2);
            
        end
        
        % КВП - функция нижнего уровня
        function [ FeatDataW, coef_map ] = fDTCWT( InData, wavenm )
        %fDTCWT Комплексное вейвлет преобразование по методу кингсбери
        %   Функция выполняет комплексное вейвлет-преобразование массива данных на 
        %   своем входе используя метод Кингсбёри 
        %   Преобразование является избыточным и формирует на выходе число
        %   отсчетов в 2^m раз превышающее число отсчетов на входе
        %   При этом выполняется дополнение нулями до размерности равной
        %   ближайшей степени двойки.
        %   Выходной 1-мерный вектор имеет формат
        %      1[re ВЧ коэф ур 1] ... [re ВЧ коэф ур n]m,    m+1[re НЧ коэф]n
        %    n+1[im ВЧ коэф ур 1] ... [im ВЧ коэф ур n]n+m,n+m+1[im НЧ коэф]2n
            [ nObj,nDim ] = size(InData);   % размерность
            decfilt = dtfilters(wavenm);    % получение коэф фильтра по названию
            filtMaxLen = max([...
                size(decfilt{1}{1},1) size(decfilt{1}{2},1)...
                size(decfilt{2}{1},1) size(decfilt{2}{2},1)]);  % определение максимальной длины фильтров
            % определение уровня разложения
            declvl = floor(log2(nDim/filtMaxLen)+1);    % уровень разложения по условию x_len>=flen*2^(L-1)
            newDim = 2^declvl * ceil(nDim/2^declvl);    % новая размерность 
            
            dimAdd = newDim-nDim;           % вычислить добавку
            if (dimAdd~=0)                  % если добавка не нулевая
                InData = [InData , zeros(nObj,dimAdd)]; % дополнить выборку нулями для выполнения условия
                nDim = newDim;                      % и обновить размерность
            end
            % вычисление коэффициентов (спойлер: число коэффициентов = nDim х2
            FeatDataW = zeros(nObj,2*nDim);         % выделение массива
            for iObj = 1:nObj                       % для всех объектов
                cwtData = dddtree(...
                    'cplxdt',InData(iObj,:),declvl,wavenm);         % выполнить преобразование 
                cfs = cellfun(@squeeze, cwtData.cfs,'UniformOutput',0); % ужать измерения
                cfs = vertcat(cfs{:});
                FeatDataW(iObj,1:2*nDim) = [cfs(:,1); cfs(:,2)].';  % присвоить де. и мн. коэффициенты
            end
            coef_map = cellfun(@length,cwtData.cfs);
        end
        
        % ВП
        function [ cfs, coef_map ] = wt( rps, varargin )
            [n_rps, dim_rp] = size(rps);
            if nargin == 1
                wave_name = 'sym4';
            elseif nargin == 2
                wave_name = varargin{1};
            else, error('Слишком много аргументов')
            end
            dec_lev = wmaxlev(dim_rp, wave_name);   % максимальный уровень декомпозиции
            [wt_fset, coef_map] = wavedec(rps(1,:), dec_lev, wave_name); % разложение первого ДП
            coef_map = coef_map(1:end - 1);         % отбросить последний элемент (длина исх. выборки)
            cfs = zeros(n_rps, length(wt_fset));    % создать матрицу коэффициентов
            cfs(1,:) = wt_fset; % скопировать коэффициенты первого разложенного ДП
            for i = 2:n_rps     % продолжить декомпозицию
                cfs(i,:) = wavedec(rps(i,:), dec_lev, wave_name);
            end
        end
        
        % ОБПФ
        function [ signal ] = ifft(fft_fvs)
            
            n_dim = size(fft_fvs,2);
            rec_spectre = sqrt(n_dim)*[ ...
                fft_fvs(:,1)    (fft_fvs(:,2:2:end-2)+1i*fft_fvs(:,3:2:end-1))/sqrt(2) ...
                fft_fvs(:,end)  (fft_fvs(:,end-2:-2:2)-1i*fft_fvs(:,end-1:-2:3))/sqrt(2) ];
            signal = ifft(rec_spectre, n_dim, 2);
        end
        
        % PCA
        function [varargout] = pca(rps)
        %pca Преобразование Карунена-Лоэва
            coeff = pca(rps);       % матрица преобразования
            pca_rps = rps*coeff;    % преобразование
            varargout{1} = pca_rps;
            if nargout>1, varargout{2} = coeff; end
        end
        
        function [dca_rps] = dca_total(rps,trgs)
        %dca Метод анализа разностных компонент
        % Выполняется расчет разностных векторов между векторами признаков, принадлежащими разным объектам
        % Это множество векторов дополняется обратными векторами. 
        % Полученное множество загоняется в АГК, в результате чего формируется преобразовнное множество
        % Аргументы:
        %   rps - входное множество (ДП)
        %   trgs - целевые векторы
        % Вывод:
        %   Результат преобразования[, коэффициенты преобразования: x*coeff]
            [object_names, ~, ci] = unique(trgs,'rows');
            [~, dim_rps] = size(rps);
            n_obj = size(object_names,1);
            rho_cell = cell(1,(n_obj^2));
            for i_obj = 1:(n_obj-1)
                for j_obj = (i_obj+1):n_obj
                    
                    c1 = rps(ci==i_obj, :);    % формирование первой матрицы
                    c2 = rps(ci==j_obj, :);    % формирование воторой матрицы
                    n_fvs = [size(c1,1) size(c2,1)];    % определение числа реализаций
                    vrhos = repmat(c1.',1,1,n_fvs(2)) - repmat(...
                        permute(c2,[2 3 1]), 1, n_fvs(1), 1);   % расчет разностных векторов
                    vrhos = reshape(vrhos, dim_rps, []);
                    rho_cell{i_obj,j_obj} = vrhos.';
                end
            end
            clear vrhos;
            rho_cell = vertcat(rho_cell{:});
            rho_cell = vertcat(rho_cell, -rho_cell);
            coeff = pca(rho_cell);  % матрица преобразования
            dca_rps = rps*coeff;    % преобразование
            
        end
        
        % DCA
        function [varargout] = dca(rps, trgs)
        %dca Метод анализа локальных разностных компонент
        % Выполняется расчет разностных векторов между ближайшими векторами признаков, принадлежащими разным объектам
        % Это множество векторов дополняется обратными векторами. 
        % Полученное множество загоняется в АГК, в результате чего формируется преобразовнное множество
        % Аргументы:
        %   rps - входное множество (ДП)
        %   trgs - целевые векторы
        % Вывод:
        %   Результат преобразования[, коэффициенты преобразования: x*coeff]
        
            [object_names, ~, ci] = unique(trgs,'rows');
            [~, dim_rps] = size(rps);
            n_obj = size(object_names,1);
            minrho = cell(1,(n_obj^2));
            for i_obj = 1:(n_obj-1)
                for j_obj = (i_obj+1):n_obj
                    
                    c1 = rps(ci==i_obj, :);    % формирование первой матрицы
                    c2 = rps(ci==j_obj, :);    % формирование воторой матрицы
                    n_fvs = [size(c1,1) size(c2,1)];    % определение числа реализаций
                    vrhos = repmat(c1.',1,1,n_fvs(2)) - repmat(...
                        permute(c2,[2 3 1]), 1, n_fvs(1), 1);   % расчет разностных векторов
                    
                    rhos = sum((vrhos).^2, 1);              % расчет метрик ( матр 1 х n_fvs(1) x n_fvs(2) )
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
            minrho = vertcat(minrho, -minrho);  % обратные векторы
            coeff = pca(minrho);    % матрица преобразования
            dca_rps = rps*coeff;    % преобразование
            varargout{1} = dca_rps;
            if nargout>1, varargout{2} = coeff; end
        end
        
    % Отображение и прочее

        % Получение карты преобразования/названий признаков
        function [ feature_map ] = get_feature_map( transform_type, coef_map )
        %GET_FEATURE_MAP Сформировать названия коэффициентов
        %   Ввод:
        %   transform_type  - название преобразования
        %   coef_map        - карта коэффициентов
            switch transform_type
                case 'none'
                    feature_map = num2str(0:coef_map-1, 't_%03d,');
                case 'afft' % должно быть передано амплитудных значений
                    feature_map = num2str(0:coef_map-1, 'f_%03d,');
                
                case 'fft' 
                    feature_map = [ 'f_000r,'...
                        num2str(reshape([1:coef_map/2-1;1:coef_map/2-1],1,[]), 'f_%03dr,f_%03di,') ...
                        num2str(coef_map/2, 'f_%03dr,')];
                    
                case 'cwt'  % должна быть передана карта коэффициентов
                    % частотные компоненты вейвлет составляющей
                    coef_map = coef_map/2;
                    max_lev = length(coef_map)-1;
                    coefs = cell(max_lev+1,1);
                    for m = 1:max_lev+1
                        cur_lev = mod(m,max_lev+1);
                        m_vect = cur_lev+zeros(coef_map(m), 1);    % вектор частотного положения
                        n_vect = (0:(coef_map(m)-1)).';     % вектор временных сдвигов
                        if coef_map(m)<100, formatSpecs = 'c_%02d_%02d';
                        elseif coef_map<1000, formatSpecs = 'c_%03d_%03d';
                        else, error('Слишком много коэффициентов')
                        end
                        coefs{m,1} = [ ... 
                            num2str([m_vect n_vect], [formatSpecs '_r,']) ...
                            num2str([m_vect n_vect], [formatSpecs '_i,']) ];
                        coefs{m,1} = reshape(coefs{m,1}.',1,[]);
                    end
                    feature_map = horzcat(coefs{:,1});
                    
                case 'acwt' %
                    max_lev = length(coef_map)-1; % максимальный уровень разложение
                    coefs = cell(max_lev+1,1);
                    for m = 1:max_lev+1
                        cur_lev = mod(m,max_lev+1);
                        m_vect = cur_lev+zeros(coef_map(m), 1);    % вектор частотного положения
                        n_vect = (0:(coef_map(m)-1)).';     % вектор временных сдвигов
                        if coef_map(m)<100, formatSpecs = 'c_%02d_%02d';
                        elseif coef_map<1000, formatSpecs = 'c_%03d_%03d';
                        else, error('Слишком много коэффициентов')
                        end
                        coefs{m,1} = num2str([m_vect n_vect], [formatSpecs ',']);
                        coefs{m,1} = reshape(coefs{m,1}.',1,[]);
                    end
                    feature_map = horzcat(horzcat(coefs{:,1}));
                    
                case 'wt'   % для вейвлет преобразования
                    max_lev = length(coef_map)-1; % максимальный уровень разложение
                    coefs = cell(max_lev+1,1);
                    for m = 1:max_lev+1
                        cur_lev = mod(m,max_lev+1);
                        dim_cur_lev = coef_map(end-m+1);
                        m_vect = cur_lev+zeros(dim_cur_lev, 1);    % вектор частотного положения
                        n_vect = ( 0:(dim_cur_lev -1) ).';     % вектор временных сдвигов
                        if coef_map(m)<100, formatSpecs = 'c_%02d_%02d';
                        elseif coef_map<1000, formatSpecs = 'c_%03d_%03d';
                        else, error('Слишком много коэффициентов')
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
                    warning('Неизвестное преобразование');
            end
            
        end
          
        % Преобразовать метаданные в строку
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
        
        % Вывести множество одномерных функций
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
                    num2str(asp_b) '^\circ}'];  % задать описание ДП
%                 if i <= 6     % добавлять текст на рис. если не более ДП 6
%                     [y, x] = max(rps(i_rp, :));    % коорд для размещения текста (глоб максимум)
%                     text(cast(x, 'double'), cast(y, 'double'), [ legend_str{i}],'VerticalAlignment','bottom');
%                 end
            end
            set(gca,'FontName','ISOCPEUR','FontAngle','italic',...
                'xGrid','on','yGrid','on',...
                'XMinorTick','on')
            if length(list) < 7; legend(legend_str); end
            
        end
        
        % вывод коэффициентов комплексного вейвлет преобразования
        function plot_cwt(h, cfs, coef_map, meta, varargin)
            if nargin == 4
                list = 1:size(cfs, 1);
            elseif nargin == 5
                list = varargin{1};
            else, error('Неправильное число аргументов')
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
                    %[z, max_i] = max(z_axis);    % коорд для размещения текста (глоб максимум)
                    %y = y_axis(1);
                    %text(x_axis(max_i), y, z, meta(i_sample).name); %
                    start = stop + 1;
                end
                h_vector(i_sample) = h_line;
                legend_list{i_sample} = horzcat(meta(i_sample).name, ...
                    num2str(meta(i_sample).asp_a), ',',...
                    num2str(meta(i_sample).asp_b) );
            end
            xlabel('Ось отсчетов'); ylabel('Уровень разложения');   zlabel('Значение амплитуды');
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
