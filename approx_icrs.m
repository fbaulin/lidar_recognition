%%

path = ("az_scan_60\icrs\");
process_icrs(path)



%% ДП Сфер
function data_est = approximate_sphere(data, meta, R_m)
    r = round(200*R_m/3);
    [n, dim] = size(data);
    x = 1:r;
    peaks_est = mean(data(:,1:2),'all');
    rp = [((r-x)/r).^2, zeros(1, dim-r)]*peaks_est;
    data_est = repmat(rp, n, 1);
end

%% ДП конусов
function data_est = approximate_cone(data, meta, H_m, R_m)
    h = 200*H_m/3; r = 200*R_m/3;   % пересчёт значений в мемтрах в отсчёты
    angles = deg2rad([meta.asp_a].');
    c = norm([h,r]);        % длина образующей
    alpha_c = acos(h/c);    % угол при вершине
    front_mask = (angles+alpha_c)<pi/2;     % разделение ракурсов
    % Вид спереди
    if any(front_mask)
        front_data = data(front_mask,:);
        front_data = cone_front_view(front_data, h, r, angles(front_mask));
    else, front_data = [];
    end
    % Вид сбоку
    if any(~front_mask)
        back_data = data(~front_mask,:);
        back_data = cone_side_view(back_data, h, r, angles(~front_mask));
    else, back_data = [];
    end
    data_est = [front_data; back_data];
%     mesh_plot(data_est, 'title', 'model');
end


function data_est = cone_side_view(data, h, r, alphas)
    [n, ~] = size(data);
    data_est = zeros(size(data),'like', data);
    peak_model = cone_peak_model(h, r, alphas);
    ending_est = estimate_endpoint(data); % конец ДП
    peaks_est = estimate_peakpoint(data, alphas);
    for i=1:n
        data_est(i,1:ending_est) = ...
            gaussmf(1:ending_est, [ending_est(i)/9, peak_model(i)]);
    end
    data_est=data_est.*peaks_est;
end


function data_est = cone_front_view(data, h, r, alphas)
    data_est = zeros(size(data), 'like', data); % подготовка массива
    peak_model = cone_peak_model(h, r, alphas);
    ending_est = estimate_endpoint(data); % конец ДП
    peaks_est = estimate_peakpoint(data, alphas);
    % геометрия задней части и геометрия параболы
    a_par = peaks_est./(ending_est-peak_model).^4;     % множитель параболы
    for i = 1:length(alphas)
        linend = peak_model(i); rpend = ending_est(i);
        if linend>=2
            p = polyfit(1:linend, data(i,1:linend), 1);
            data_est(i,1:linend) = polyval(p, 1:linend); % формировние линейной части
        end
        if (rpend-linend)>=2
            data_est(i,linend+1:rpend) = a_par(i)*((linend+1:rpend)-rpend).^4;
        end
    end
end


function [varargout] = cone_peak_model(h, r, alphas)
    c = norm([h,r]);        % длина образующей
    alpha_c = acos(h/c);    % угол при вершине
    peak_values = round(c*(-0.5).^((alpha_c+alphas)>pi/2).*cos(alphas+alpha_c)); % конец линейного участка
    varargout{1} = peak_values;
    if nargout>1
        varargout{2} = alpha_c;
    end
end


function peaks = estimate_peakpoint(data, varargin)
    % аппроксимация значений в максимуме
    if nargin==2, x_axis = varargin{1}; 
    else, x_axis = 1:size(data,2);
    end
    peaks_v = max(data, [], 2); % величины максимумов
    p = polyfit(x_axis, peaks_v, 2);
    peaks = polyval(p, x_axis);
end 


function ending = estimate_endpoint(data)
    [n, dim] = size(data);
    peak_v = max(data, [], 2);
    ending = arrayfun(@(x)find(data(x,:)>(peak_v(x)/50),1,'last'),1:n).';
    % два случая: ending - последняя точка и не последняя точка
    ap_mask = ending<dim;       % для аппрокс используются только не последние точки
    angle_axis = (1:n).';
    p = polyfit(angle_axis(ap_mask), ending(ap_mask), 2);
    ending(ap_mask) = round(polyval(p, angle_axis(ap_mask)));
end


function data = avrg_mean(data, varargin)
    % Find the moving average of the data and plot it against the original data.
    window_size = get_value(varargin, 'window_size', 5);
    b = (1/window_size)*ones(1,window_size);
    a = 1;
    data = filter(b,a,data);

end


%% Служебные функции
function ax = mesh_plot(data, varargin)
    plot_title = get_value(varargin, 'title', 'ICRS');
    angle_limits = get_value(varargin, 'angles', [1, size(data,1)]);
    time_limit = get_value(varargin, 'time', size(data,2)-1);
    figure('Name',plot_title, 'color', 'white','WindowStyle','docked');
    hold on;
    [X, Y] = meshgrid(...
        linspace(0, time_limit, size(data, 2)),...
        linspace(angle_limits(1), angle_limits(2), size(data, 1))...
        );
    mesh(X, Y, data); ax = gca;
    axis vis3d; view(45,20);
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

function process_icrs(icr_path)
    
    load(icr_path + "con1x1.mat", 'icrs', 'meta');
    icrs_int = avrg_mean(icrs, meta, 'window_size',3);
    icrs = approximate_cone(icrs_int, meta, 1, 1);
    mesh_plot(icrs, 'angles', [0 30], 'time', 3, 'title', 'icrs');
    save(icr_path +"m_" +"con1x1.mat",'icrs','meta')
    
    load(icr_path + "con3x1.mat", 'icrs', 'meta');
    icrs_int = avrg_mean(icrs, meta, 'window_size',3);
    icrs = approximate_cone(icrs_int, meta, 3, 1);
    mesh_plot(icrs, 'angles', [0 30], 'time', 3, 'title', 'icrs');
    save(icr_path +"m_" +"con3x1.mat",'icrs','meta')
    
    load(icr_path + "sphx15.mat", 'icrs', 'meta');
    icrs_int = avrg_mean(icrs, meta, 'window_size',3);
    icrs = approximate_sphere(icrs_int, meta, 1.5);
    mesh_plot(icrs, 'angles', [0 30], 'time', 3, 'title', 'icrs');
    save(icr_path +"m_" +"sphx15.mat",'icrs','meta')

end