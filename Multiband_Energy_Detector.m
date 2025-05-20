%% ------------------ Multiband Energy Detector by Theresia Muhr ---------------------------
clc; clear; close all;

%% Define recordings
recording = {"your_file1.wav"};
% recording = {"your_file1.wav", "your_file2.wav"};

%% Define Detection Parameters
flow = 100;        % Lower frequency bound
fhigh = 4000;      % Higher frequency bound
bands_per_decade = 10; % Number of frequency bands
chunk_length = 1;  % Length of chunks in seconds
chunk_overlap = 0.5; % Overlap between chunks (50%)
sta_length = 3;    % STA window length in minutes
threshold = 0.5;   % Detection threshold
min_gap = 1;       % Minimum gap between events in windows

%% Logarithmic Band spacing
num_bands = round(bands_per_decade * log10(fhigh / flow)); % Total number of bands
f_center = flow * 10.^( (0:num_bands - 1) / bands_per_decade);
f_lower = f_center / 10^(0.5 / bands_per_decade); %Lower edge
f_upper = f_center * 10^(0.5 / bands_per_decade); %Upper edge
f_upper = [f_upper, f_upper(end) * 10^(0.5 / bands_per_decade)]; %One extra element for proper indexing

%% Process each recording
for rec_idx = 1:length(recording)
    fprintf('Processing recording %d/%d\n', rec_idx, length(recording));
    
    %% Extract audio file information
    current_rec = recording{rec_idx};
    audio_file = current_rec{1};
    
    %% Get audio file info
    audio_info = audioinfo(audio_file);
    total_duration = audio_info.Duration; % Length of the recording in seconds

    lta_length = total_duration / 60;   % LTA window length in minutes = length of the recording in minutes 

    fs = audio_info.SampleRate; %Sampling frequency
    
    %% Prepare processing parameters
    sta_samples = round(sta_length * 60 * fs); % Amount of samples per STA window
    sta_duration = sta_length * 60;  % STA duration in seconds
   
    %% Preallocate storage for all events
    all_events = [];
    
    %% Prepare start and end times for STA windows
    start_time = 0 : sta_duration : total_duration;  % STA Window start time in seconds
    end_time = min(start_time + sta_duration, total_duration); % STA Window end time in seconds

    %% Preallocate storage for energy per band per STA window
    energy_per_band = nan(length(start_time), num_bands);

    %% Process file in STA-length chunks
    for z = 1 : length(start_time) % Loop through each STA-sized chunk
        current_start_time = start_time(z); % Get the current start time from the array
        current_end_time = end_time(z); % Get the current end time from the array

        start_sample = round(current_start_time * fs) + 1; %+1 because MATLAB indexing starts at 1
        end_sample = round(current_end_time * fs); % End time in samples

        %% Processing info
        start_timestamp = sec_to_ddhhmmss(current_start_time); % Self-written function, converts sec into dd:hh:mm:ss format
        end_timestamp = sec_to_ddhhmmss(current_end_time); % Self-written function, converts sec into dd:hh:mm:ss format
        
        %% Extract audio chunk using audioread
        [signal_v, ~] = audioread(audio_file, [start_sample, end_sample]); % Extract only one STA window using sample values

        %% Convert Voltage to Pascal
        sens_linear = 10^(-187.2 / 20);
        signal = signal_v / sens_linear;
        
        %% Prepare processing parameters
        chunk_samples = round(chunk_length * fs); % Number of samples per chunk (chunk is 1 second)
        chunk_overlap_samples = round(chunk_samples * chunk_overlap); % Number of samples for overlap
        
        %% Compute chunks with overlap
        num_chunks = floor((length(signal) - chunk_overlap_samples) / (chunk_samples - chunk_overlap_samples)); % Number of chunks per STA window
        band_energies = nan(num_chunks, num_bands); % Preallocation 
        
        %% Process each chunk in the STA window
        for i = 1:num_chunks
            % Calculate chunk indices
            start_idx = 1 + (i-1) * (chunk_samples - chunk_overlap_samples);
            end_idx = start_idx + chunk_samples - 1;
            
            % Ensure signal length doesn't exceed the available length
            if end_idx > length(signal)
                break;
            end
            
            % Extract and window chunk
            chunk = signal(start_idx:end_idx); % Get the chunk 
            window = blackmanharris(length(chunk)); % Create the Blackman-Harris window
            chunk_windowed = chunk .* window; % Apply Blackman-Harris window
                    
            % Compute FFT of chunk
            fft_chunk = fft(chunk_windowed);
            fft_mag = abs(fft_chunk(1:floor(length(fft_chunk)/2+1))).^2; % Power spectrum
            freqs = linspace(0, fs/2, length(fft_mag)); % Frequency axis
            
            % Compute energy for each frequency band
            for j = 1:num_bands % Loop trough each frequency band
                if j == num_bands  % Last frequency band (upper bound inclusive)
                    band_idx = (freqs >= f_lower(j)) & (freqs <= f_upper(j)); % Select which bins from the FFT belong to the current frequency band j
                else
                    band_idx = (freqs >= f_lower(j)) & (freqs < f_upper(j+1)); % Select which bins from the FFT belong to the current frequency band j
                end

                % Calculate the energy for the current band by averaging the magnitudes
                band_energies(i, j) = mean(fft_mag(band_idx));  % Average power per band (i = chunk ; j = freq band)
            end
        end
        
        %% Store mean energy for each band in this STA window
        energy_per_band(z, :) = mean(band_energies, 1);  % Mean energy across chunks for each band
    end
    
    %% Compute LTA (background noise level) for each frequency band
    num_short_in_long = round(lta_length * 60 / (sta_length * 60)); % How many short windows are used for the long window
    mean_bn = nan(length(energy_per_band), num_bands); % Preallocate
    
    for q = 1:length(energy_per_band)
        start_long_idx = max(1, q - floor(num_short_in_long / 2));
        end_long_idx = min(length(energy_per_band), q + floor(num_short_in_long / 2));
        
        % For each frequency band, calculate the LTA
        for j = 1:num_bands
            mean_bn(q, j) = mean(energy_per_band(start_long_idx:end_long_idx, j));
        end
    end
    
    %% Compute STA/LTA ratio for each band
    energy_ratio = energy_per_band ./ mean_bn;  % STA / LTA
    
    %% Apply detection threshold to max of the ratios
    energy_sig_fin = max(energy_ratio, [], 2);  % Maximum ratio across all bands
    detection = energy_sig_fin >= threshold;
    
    %% Group detections that are close together into single events
    chunk_events = [];
    in_event = false;  % To track if in an event
    event_start = 0;   % Initialize event_start
    event_end = 0;     % Initialize event_end
    
    for i = 1:length(detection)
        if detection(i)  % If detection is true (1)
            if ~in_event  % If we are not already in an event, start a new event
                event_start = i;  % Set the start of the event to the current index
                in_event = true;  % Mark that we're in an event
            end
            % Always update the end when we have a detection
            event_end = i;
        else  % If detection is false (0)
            if in_event  % If we were in an event and now we're not
                % Add the completed event to our list
                chunk_events = [chunk_events; event_start event_end];
                in_event = false;  % Mark that we're no longer in an event
            end
        end
    end
    
    % If we're still in an event at the end of the loop, close the event
    if in_event
        chunk_events = [chunk_events; event_start event_end];
    end
    
    %% Convert detected events to global time in seconds
    all_events = [];
    if ~isempty(chunk_events)
        % Convert STA window indices to global time in seconds
        % Each index 'i' corresponds to the i-th STA window of duration sta_duration
        all_events = [all_events; 
                     (chunk_events(:,1) - 1) * sta_duration + start_time(1), ... % Start time in seconds
                     chunk_events(:,2) * sta_duration + start_time(1)];          % End time in seconds
    end
    
    %% Convert detected events to timestamps in dd:hh:mm:ss format
    num_events = size(all_events, 1);
    event_timestamps = cell(num_events, 2);
    
    for i = 1:num_events
        % Convert seconds to dd:hh:mm:ss format using sec_to_ddhhmmss
        event_timestamps{i,1} = sec_to_ddhhmmss(all_events(i,1));
        event_timestamps{i,2} = sec_to_ddhhmmss(all_events(i,2));
    end
    
    % Show the detection table
    event_table = table((1:num_events)', ...
        event_timestamps(:,1), ...
        event_timestamps(:,2), ...
        'VariableNames', {'EventIndex', 'StartTime', 'EndTime'});
    
    disp('Event Detection Table for Recording ' + recording{rec_idx} + ' with Threshold: ' + threshold);
    disp(event_table);

    %% Plot 
    figure('Name', 'Detection Timeline', 'Position', [400 400 700 400], 'color', 'w');
    
    % Create time axis for STA windows
    time_axis = (0:size(energy_per_band, 1)-1) * sta_duration;
    
    % Plot the timeline with threshold
    plot(time_axis, energy_sig_fin, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(time_axis, threshold*ones(size(energy_sig_fin)), 'r--', 'LineWidth', 1.5);
    xlabel('Time (s)', 'interpreter', 'latex', 'Fontsize', 12);
    xlim([min(time_axis) max(time_axis)]);
    ylabel('Max STA/LTA Ratio', 'interpreter', 'latex', 'Fontsize', 12);
    title('Recording: ' + recording{rec_idx}, 'interpreter', 'latex', 'Fontsize', 12);
    grid on;
    
    % Highlight detection periods
    for i = 1:size(all_events, 1)
        event_start = all_events(i, 1);
        event_end = all_events(i, 2);
        
        % Fill area with semi-transparent color
        y_min = 0;
        y_max = max(energy_sig_fin) * 1.1;
        fill([event_start, event_end, event_end, event_start], [y_min, y_min, y_max, y_max], ...
             'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        
        % % Uncomment the two lines below to add event number label in Plot
        % text((event_start + event_end)/2, y_max*0.9, ['Event ', num2str(i)], ...
        %      'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
    legend('STA/LTA Ratio', sprintf('Detection Threshold: %.2f', threshold), 'Detected Events', 'interpreter', 'latex', 'Fontsize', 12);
    ylim([0, max(max(energy_sig_fin)*1.1, threshold*1.5)]);

end
