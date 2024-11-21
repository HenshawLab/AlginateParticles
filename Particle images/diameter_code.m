clear all;
clc;

% Set up directories
DataDir = 'C:\Users\romanelli\Desktop\MATLAB_code_particle'; % Update this to your actual folder path
nd2FileName = fullfile(DataDir, 'IBIDI_102024_time2.nd2');
OutputMainDir = 'C:\Users\romanelli\Desktop\MATLAB_code_particle\AnalysisOutput';
mkdir(OutputMainDir);

% Add Bio-Formats toolbox to MATLAB path
addpath('C:\Users\romanelli\Desktop\bfmatlab'); % Update this path

% Verify if BioformatsImage is accessible
if exist('BioformatsImage', 'class') ~= 8
    error('BioformatsImage is not accessible. Make sure the bfmatlab folder path is correct.');
end

% Parameters
NumChannels = 2;
NumZStacks = 11;
TimeInterval = 20; % minutes per time point
TotalDuration = 20 * 69; % duration in minutes for 20+ hours
NumTimePoints = ceil(TotalDuration / TimeInterval);
    
% Initialize tables to store radius and fluorescence data
radius_data = table();
fluorescence_data = table();

% Check if the .nd2 file exists
if ~isfile(nd2FileName)
    error(['ND2 file does not exist: ' nd2FileName]);
end

% Open the ND2 file
bfr = BioformatsImage(nd2FileName);

% Set output directory for Series 1
OutDir = fullfile(OutputMainDir, 'Series1');
mkdir(OutDir);

% Initialize arrays to store radius and fluorescence measurements
radii_series = nan(NumTimePoints, 1); % Radius for each time point
fluorescence_series = nan(NumTimePoints, 1); % Fluorescence intensity for each time point

% Loop through each timepoint for Series 1
for t = 1:NumTimePoints
    % Initialize mask as empty at the start of each timepoint
    mask = []; 

    % Loop through each channel
    for ch = 1:NumChannels
        % Preallocate stack for Z-stacks
        img_stack = zeros(bfr.height, bfr.width, NumZStacks);

        % Load each Z-plane for the current series, channel, and timepoint
        for z = 1:NumZStacks
            img = getPlane(bfr, 1, ch, t, z); % Load Series 1
            img_stack(:,:,z) = mat2gray(img); % Normalize intensity
        end

        % Process the second channel for particle detection and radius measurement
        if ch == 2
            % Use max projection across Z-stacks
            max_proj = max(img_stack, [], 3);

            % Initialize radii and centers as empty for safety
            radii = []; 
            centers = []; 

            % Detect particle and measure radius in Channel 2
            [centers, radii] = imfindcircles(img_stack(:,:,z), [799, 801], 'ObjectPolarity', 'bright', 'Sensitivity', 0.9);

            if ~isempty(radii)
                % Assume the largest detected circle is the particle
                radii_series(t) = max(radii); % Keep only the largest radius detected
                particle_center = centers(radii == max(radii), :);

                % Create a circular mask based on the detected radius and center
                [X, Y] = meshgrid(1:size(max_proj, 2), 1:size(max_proj, 1));
                mask = ((X - particle_center(1)).^2 + (Y - particle_center(2)).^2) <= radii_series(t)^2;

                % Visualization: Display each frame with detected circle
                figure;
                imshow(max_proj, []); % Display the max-projected image
                hold on;
                viscircles(particle_center, radii_series(t), 'EdgeColor', 'r'); % Overlay detected circle
                title(['Timepoint: ' num2str(t) ', Radius: ' num2str(radii_series(t)) ' pixels']);
                hold off;
                pause(0.5); % Pause to allow viewing each frame
            else
                disp(['No particle detected in Series 1, Timepoint ' num2str(t)]);
                mask = []; % No particle detected, mask is empty
            end
        end

        % Process the first channel for fluorescence measurement using the mask
        if ch == 1 && ~isempty(mask)
            % Use max projection across Z-stacks for fluorescence measurement
            max_proj = max(img_stack, [], 3);

            % Apply mask to extract fluorescence within detected particle area
            masked_image = max_proj .* mask;
            fluorescence_series(t) = sum(masked_image(:)); % Sum of intensity inside the mask
        end
    end
end

% Save Series 1 data to tables for radius and fluorescence
radius_data = table((1:NumTimePoints)' * TimeInterval, radii_series, ...
                    'VariableNames', {'Time_Minutes', 'Radius_Series1'});
fluorescence_data = table((1:NumTimePoints)' * TimeInterval, fluorescence_series, ...
                          'VariableNames', {'Time_Minutes', 'Fluorescence_Series1'});

% Export data to CSV files
writetable(radius_data, fullfile(OutputMainDir, 'radius_data_Series1.csv'));
writetable(fluorescence_data, fullfile(OutputMainDir, 'fluorescence_data_Series1.csv'));

disp('Analysis complete for Series 1. Radius and fluorescence data saved to CSV files.');
