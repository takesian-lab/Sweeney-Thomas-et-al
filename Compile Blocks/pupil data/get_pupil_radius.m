function [Radius, Eye_Area] = get_pupil_radius(csvfiles, DLC_training_name, plot_figure)
%% Load csv files made with DeepLabCut to extract pupil measurements
% Written by Selorm Quarshie, August 2022
%
% Argument(s):
%   csvfiles (cell array) list of csv filenames to analyze
%   DLC_training_name (string) name of DLC training file used to create csvs
%
% Returns:
%   Radius (cell) nFiles x 1 cell with pupil radius for each frame
%   Eye_Area (cell) nFiles x 1 cell with eyelid area for each frame
%
% Version history:
% - V1 = current version

% Notes:
% - likelihood thresholding works very well for pupil but not for eyelid

%% Setup

pixel_size_multiplier = 1; %We could use this if we ever want to try converting from pixels to cm
likelihood_threshold = 0.99; %Cutoff for determining DLC points to exclude

switch DLC_training_name
    
    case 'DLC_resnet101_Pupillometry_01_15_23Jan15shuffle1_500000' 
        %'DLC_resnet101_Pupillometry 07-21-22Jul21shuffle1_500000' <- first iteration of model
        
        %Columns are:
        %2:4    upper palebra maxima
        %5:7    upper palebra medial
        %8:10   upper palebra lateral
        %11:13  lower palebra minima
        %14:16  lower palebra medial
        %17:19  lower palebra lateral
        %20:22  medial canthus
        %23:25  lateral canthus
        %26:28  pupil maxima
        %29:31  pupil minima
        %32:34  pupil medial
        %35:37  pupil lateral
        %38:40  pupil maxima x medial
        %41:43  pupil maxima x lateral
        %44:46  pupil minima x medial
        %47:49  pupil minima x lateral
        
        %Columns in .csv corresponding to pupil
        pupil_x_ind = 26:3:49;
        pupil_y_ind = pupil_x_ind + 1;
        pupil_likelihood_ind = pupil_x_ind + 2;
        
        %Columns in .csv corresponding to eyelid
        upper_palebra_x_ind = [2:3:8 20 23];
        upper_palebra_y_ind = upper_palebra_x_ind + 1;
        upper_eyelid_likelihood_ind = upper_palebra_x_ind + 2;
        lower_palebra_x_ind = 11:3:23;
        lower_palebra_y_ind = lower_palebra_x_ind + 1;
        lower_eyelid_likelihood_ind = lower_palebra_x_ind + 2;
        
    otherwise
        error('DLC training name does not match known training datasets')
end

%% Iterate over csv files
%iterate over cells in csv for every pupil radius

Radius = cell(size(csvfiles));
Eye_Area = cell(size(csvfiles));

for i = 1:length(csvfiles)
    PF = csvfiles{i};
    M = readmatrix(PF);
    
    pupil_x = M(: , pupil_x_ind);        %all the x coordinates for the pupil
    pupil_y = M(: , pupil_y_ind);        %all the y coordinates for the pupil

    upper_palebra_x = M(: , upper_palebra_x_ind);           %all the xy coordinates for the upper eyelid
    upper_palebra_y = M(: , upper_palebra_y_ind);
    lower_palebra_x = M(: , lower_palebra_x_ind);            %all the xy coordinates for the lower eyelid
    lower_palebra_y = M(: , lower_palebra_y_ind);

    frames = M(:,1);
        
    %%Removing DLC points with low likelihood values

    pupil_likelihood = M(: , pupil_likelihood_ind);                         %pupil likelihood values from DLC
    pupil_likelihood_filter = pupil_likelihood >= likelihood_threshold;     %creates threshold for omitting pupil values based on likelihood
    
    upper_eyelid_likelihood = M(: , upper_eyelid_likelihood_ind);                        %eyelid likelihood values from DLC
    upper_eyelid_likelihood_filter = upper_eyelid_likelihood >= likelihood_threshold;    %creates threshold for omitting eyelid values based on likelihood
    lower_eyelid_likelihood = M(: , lower_eyelid_likelihood_ind);
    lower_eyelid_likelihood_filter = lower_eyelid_likelihood >= likelihood_threshold;
    
    %Threshold pupil points
    pupil_x(pupil_likelihood_filter == 0) = nan;
    pupil_y(pupil_likelihood_filter == 0) = nan;

    %Threshold eyelid points
    upper_palebra_x(upper_eyelid_likelihood_filter == 0) = nan;
    upper_palebra_y(upper_eyelid_likelihood_filter == 0) = nan;
    lower_palebra_x(lower_eyelid_likelihood_filter == 0) = nan;
    lower_palebra_y(lower_eyelid_likelihood_filter == 0) = nan;

    %%concatenates upper and lower eyelid values
    eyelid_x = horzcat(upper_palebra_x, lower_palebra_x);       
    eyelid_y = horzcat(upper_palebra_y, lower_palebra_y);
    
    %%Compute pupil radius and eyelid area for each frame
    
    %%Preallocate variables
    xcenter = nan(size(frames));
    ycenter = nan(size(frames));
    radius = nan(size(frames));
    eye_area = nan(size(frames));
    
    for j = 1:length(frames)        
       
        %Find good points for the pupil
        temp_pupil_xy = [pupil_x(j,:); pupil_y(j,:)];       %pair pupil data
        pupil_nanColumns = sum(isnan(temp_pupil_xy)) > 0;   %find points where either x or y are nan
        temp_pupil_xy(:, pupil_nanColumns) = [];            %remove those points
        
        %Only compute circle from pupil if we have at least 3 points
        if sum(pupil_nanColumns) < 3
        
            [xc, yc , R, a] = circfit(temp_pupil_xy(1,:), temp_pupil_xy(2,:));

            xcenter(j) = xc;
            ycenter(j) = yc;
            radius(j) = R * pixel_size_multiplier;
        end
        
        %Find good points for the eyelids
        temp_eyelid_xy = [eyelid_x(j,:); eyelid_y(j,:)];    %pair eyelid data
        eyelid_nanColumns = sum(isnan(temp_eyelid_xy)) > 0; %find points where either x or y are nan
        temp_eyelid_xy(:, eyelid_nanColumns) = [];          %remove those points
        
        %Only compute eyelid area if we have at least 3 points
        if sum(eyelid_nanColumns) < 3
            %NOTE: We are still working on the best computation for eyelid area
            [k,v] = boundary(temp_eyelid_xy(1,:)', temp_eyelid_xy(2,:)'); 
            eye_area(j) = v * pixel_size_multiplier;
        end
        
        %Plot figure to troubleshoot
        plot_test_figure = 0;
        if plot_test_figure
            figure; hold on
            scatter(temp_pupil_xy(1,:), temp_pupil_xy(2,:), 'red') %Pupil points
            scatter(xc, yc, 'red', 'x') %Center of the pupil
            line([xc xc+R], [yc yc], 'Color', 'red') %Radius of the pupil
            scatter(temp_eyelid_xy(1,:), temp_eyelid_xy(2,:), 'blue')
            scatter(upper_palebra_x(j,:), upper_palebra_y(j,:), 'cyan')
            scatter(lower_palebra_x(j,:), lower_palebra_y(j,:), 'blue')
        end
    end
    
    Radius{i} = radius;
    Eye_Area{i} = eye_area; 
end

%%  Plot concatenated trials

if plot_figure
    
    %Concatenate data from all trials
    concat_radius = [];
    concat_eyelid = [];
    trial_start = nan(size(Radius));
    for i = 1:length(Radius)
        trial_start(i) = length(concat_radius);
        concat_radius = [concat_radius, Radius{i}'];
        concat_eyelid = [concat_eyelid, Eye_Area{i}'];
    end
    trial_start = trial_start + 1;

    %Plot figure
    figure
    subplot(2,1,1); hold on
    title('Pupil')
    plot(concat_radius)
    vline(trial_start)
    for i = 1:length(trial_start)
        text(trial_start(i), nanmedian(concat_radius)*1.15, num2str(i))
    end
    h = gca;

    subplot(2,1,2); hold on

    title('Eyelid')
    plot(concat_eyelid)
    vline(trial_start)
    for i = 1:length(trial_start)
        text(trial_start(i), nanmedian(concat_eyelid)*1.5, num2str(i))
    end
    h(2) = gca;
    
    linkaxes(h, 'x');
    
end
