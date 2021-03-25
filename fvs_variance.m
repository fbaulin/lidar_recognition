% fvs variance estimation script
% to be placed in lidar_recognition project path

clear;
clc;

% project abs path (to be changed for win)
default_abs_dir = '/Users/ivnvalex/Documents/BMSTU/nirs2021/studvesna/lidar_recognition'; 
% fvs path
fvs_dir = 'test/fvs/';  % k = 0 fvs should be in 'test/fvs/0', etc

transforms = {'none', 'afft', 'acwt', 'pca'};  % transfortms cell array

for k = 0:10
    
    kdir = fvs_dir + string(k);  % dir of the k-th fvs
    cd(default_abs_dir);  % to make sure we're in the project dir
    cd(kdir);  % k-th fvs dir
    
    for t = 1:length(transforms)  % for every transform
        
        transform_str = transforms{t};  % current transform str
        formatSpec = '*%s*.mat';  % file format 
        str = sprintf(formatSpec, transform_str);  % file format pattern *transform_str*.mat
        file = dir(str);  % fvs filename with str pattern
        transform = load(file.name);  % load current fvs

        a = size(transform.fvs);  % group fvs size
        n_objs = a(2);  % number of objects
        n_angs = 36;  % number of angles
        n_reps = length(transform.fvs{1,1}) / n_angs;  % number of repeats for each angle
        b = size(transform.fvs{1,1});  % each fvs size
        n_points = b(2);  % number of points in the fvs
        reps = zeros(n_reps, n_points);  % zero matrix for every repeat of an angle
        angs = zeros(n_angs, n_objs);  % zero matrix for every angle

        for obj = 1:n_objs  % for every object
            for ang = 1:n_angs  % for every angle
                for rep = 1:n_reps  % for every repeat of an angle
                    reps(rep, :) = transform.fvs{1, obj}(ang + (rep-1)*n_angs, :);  % fvs for every repeat of an angle
                end
                angs(ang, obj) = mean(var(reps));  % mean variance for every angle and object
            end
        end
        
        variance = mean(mean(angs));  % mean variance for the angles and the objects
        formatSpec = 'k = %d: %s variance = ';  % output pattern
        str = sprintf(formatSpec, k, transform_str);  % output string
        disp(str + string(variance));  % output
       
    end
end