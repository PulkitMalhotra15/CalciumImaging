function [  ] = analysis(  )
%ANALYSIS Summary of this function goes here
%   Detailed explanation goes here
root_fp = '\Users\mshtrahman\Desktop\Pilo Stuff\Shannon Final Pilo Analysis Files';
test_fp = [root_fp,'Middle Aged Mice\control'];

%% Load the imported data
imported_data = load( [ test_fp, '20171108_PiloMouse0_001 intensity data.mat'] );

% Loads the intensity data, and the firing cells' intensity data converted
% into binary.
intensityData = imported_data.intensityData; % [totalCells×frames]
binaryFiring = imported_data.binaryFiring; % [totalActiveCells×frames]
[totalActiveCells, num_frames] = size(binaryFiring);
totalCells = imported_data.totalCells;
% Less important data
firedNeurons = imported_data.firedNeurons; % list `totalActiveCells` long
threshold = imported_data.threshold;percentActiveCells = totalActiveCells/totalCells;
percentActiveCells = totalActiveCells/totalCells;
% Note, the following are not loaded nor necessary: amplitudes, number_of_events, firingTimes



%% Quantify activity
    function[num_active_cells_in_window] = get_num_active_cells_window( window_start, window_size )
        % Returns the number of active cells in a given window
        binaryFiring_window = binaryFiring(:,window_start:window_start+window_size);
        num_active_cells_in_window = 0;
        for cell_i=1:totalActiveCells
            binaryFiring_window_single_cell = binaryFiring_window(cell_i,:);
            % cell_fires=1 if binaryFiring_window_single_cell contains a 1
            % cell_fires=0 if binaryFiring_window_single_cell does not contain 1
            cell_fires = any(binaryFiring_window_single_cell);
            if cell_fires
                num_active_cells_in_window = num_active_cells_in_window + 1;
            end
        end
    end
    function[segmented_active_cells] = get_segmented_active_cells( ...
            window_size, as_percent, ignore_partial_window )
        % Returns list of the activity windowed
        % if as_percent=true: returns percent active in each window bin
        % if as_percent=false: returns number active in each window bin
        % if ignore_partial_window=true: only returns activity for each
        %        FULL window (typically ignores last bin)
        segmented_active_cells = [];
        
        for window_start=1:window_size:num_frames
            % If this is the last window and we are ignoring the last
            % window, skip
            if window_start+window_size>num_frames && ignore_partial_window
                fprintf('Incomplete window after frame %d (out of %d), ignoring remaining partial window\n\n'...
                    ,window_start, num_frames)
                break
            end
            
            num_active_in_window = get_num_active_cells_window( window_start, window_size );
            segmented_active_cells = [segmented_active_cells,num_active_in_window];
        end
        
        if as_percent
            segmented_active_cells = 100*segmented_active_cells/totalActiveCells;
        end
        
    end

% Make histogram of binned cell activity
window_size_in_frames = 25;
return_as_percent = true;
ignoring_partial_window = true;
segmented_active_cells = get_segmented_active_cells( ...
    window_size_in_frames, return_as_percent, ignoring_partial_window );

hist_step = 5;
hist_largest_val = 100;
edges = (0:hist_step:hist_largest_val);
nbins = hist_largest_val/hist_step;
hist(segmented_active_cells,edges,nbins)

title(['Histogram of windowed cell activity, window size: ',string(window_size_in_frames)])
xlabel('Percent active cells') 
ylabel('Number of bins with this percent of active cells') 
percentActive_mean = mean(segmented_active_cells);
percentActive_std = std(segmented_active_cells);
fprintf('__Windowed Percent Active Cells__\n')
fprintf('  -  mean:%d    std:%d   \n\n',percentActive_mean,percentActive_std)
pause();




%% Quantify events
    function[ spike_indices ] = get_indices_of_initial_spike( binaryFiring_single_cell )
        % Return binaryFiring array, ignoring indices after first spike
        % frame. Example: [0,1,1,1,0,1,1] converted to [0,1,0,0,0,1,0]
        spike_indices = [];
        for frame=2:num_frames
            val_at_curr_frame = binaryFiring_single_cell(frame);
            val_at_prev_frame = binaryFiring_single_cell(frame-1);
            if val_at_curr_frame==1 && val_at_prev_frame==0
                spike_indices = [spike_indices, frame];
            end
        end
    end
    function[amplitude_array] = get_amplitude_array( as_list )
        % Returns an array of ALL amplitudes for this movie
        % as_list will determine whether the amplitudes are returned as a
        % list or as a cell_array
        if as_list
                amplitude_array = [];
            else
                amplitude_array = {};
        end
        
        for cell_i=1:totalActiveCells
            cell_num = firedNeurons(cell_i);
            binaryFiring_single_cell = binaryFiring(cell_i,:);
            intensityData_single_cell = intensityData(cell_num,:);
            
            initial_spike_indices = get_indices_of_initial_spike( binaryFiring_single_cell );
            amplitudes = intensityData_single_cell( initial_spike_indices );
            if as_list
                amplitude_array = [amplitude_array,amplitudes];
            else
                amplitude_array{cell_i} = amplitudes;
            end
        end
        
    end
% Make histogram of event amplitudes
amplitude_list = get_amplitude_array( true );
hist( amplitude_list, 50 )

title('Histogram of windowed event amplitudes')
xlabel('Firing amplitude') 
ylabel('Number of times a cell has fired with this amplitude') 
amplitude_mean = mean(amplitude_list);
amplitude_std = std(amplitude_list);
fprintf('__Amplitude__\n')
fprintf('  -  mean:%d    std:%d   \n\n',amplitude_mean,amplitude_std)
pause();


    function[ firing_frames ] = get_start_and_end_firing_frames_single_cell( binaryFiring_single_cell )
        firing_frames = [];
        event_number = 1; % The event number is the index of firing_times
        for frame=2:num_frames
            val_at_curr_frame = binaryFiring_single_cell(frame);
            val_at_prev_frame = binaryFiring_single_cell(frame-1);
            % If the current frame is a start frame 
            if val_at_curr_frame==1 && val_at_prev_frame==0
                start_frame = frame;
                for subframe=frame:num_frames
                    subVal_at_curr_frame = binaryFiring_single_cell(subframe);
                    subVal_at_prev_frame = binaryFiring_single_cell(subframe-1);
                    if (subVal_at_curr_frame==0 && subVal_at_prev_frame==1) || subframe==num_frames
                        end_frame = subframe-1;
                        break
                    end
                end
                firing_frames{event_number} = [start_frame,end_frame];
                event_number = event_number+1;
            end
        end
    end
    function[ area_list ] = get_area_list( )
        
        area_list = [];
        
        for cell_i=1:totalActiveCells
            cell_num = firedNeurons(cell_i);
            binaryFiring_single_cell = binaryFiring(cell_i,:);
            intensityData_single_cell = intensityData(cell_num,:);
            % firing_times encoded as: 
            %  start_time = firing_times{firing_index}(1)
            %  end_time = firing_times{firing_index}(2)
            firing_times = get_start_and_end_firing_frames_single_cell( binaryFiring_single_cell );

            for firing_i=1:length(firing_times)
                start_time = firing_times{firing_i}(1);
                end_time = firing_times{firing_i}(2);
                area_under_peak = 0;
                for firing_time=start_time:end_time
                    area_under_peak = area_under_peak + intensityData_single_cell(firing_time);
                end
                area_list = [area_list, area_under_peak];
            end
            
        end 
    end

% Make histogram of event areas
area_list = get_area_list( );
hist( area_list, 35 )

title('Histogram of areas under firing curve')
xlabel('Area under firing curve') 
ylabel('Number of events with this area') 
area_mean = mean(area_list);
area_std = std(area_list);
fprintf('__Area__\n')
fprintf('  -  mean:%d    std:%d   \n\n',area_mean,area_std)
pause();

% Can check first value with:
% id.intensityData( 2,id.binaryFiring(1,:)==1 ) == [0.2694,    0.0992,    0.0822]




%% Quantify synchrony
    function[synchrony_array] = get_synchronous_firing_bins( )
        % Return a NON-normalized cell array of synchronous firing
        % synchrony_array{2}=number of times two cells fire on the same frame
        synchrony_array = {};
        for i=1:10
            synchrony_array{i}=0;
        end
        firing_times = {};
        
        for frame=1:num_frames
            binaryFiring_all_cells_one_frame = binaryFiring(:,frame);
            % Summing all binaryFiring data over 1 frame gives number of
            % cells firing at the same time
            number_of_synchronous_cells = sum(binaryFiring_all_cells_one_frame(:) == 1);
            if number_of_synchronous_cells~=0
                if length(synchrony_array{number_of_synchronous_cells})==0
                    synchrony_array{number_of_synchronous_cells} = 1;
                else
                    synchrony_array{number_of_synchronous_cells} = ...
                        synchrony_array{number_of_synchronous_cells} + 1;
                end
            end
        end
    end

    function[ firing_frames ] = get_start_firing_frames_single_cell( binaryFiring_single_cell )
        firing_frames = [];
        event_number = 1; % The event number is the index of firing_times
        for frame=2:num_frames
            val_at_curr_frame = binaryFiring_single_cell(frame);
            val_at_prev_frame = binaryFiring_single_cell(frame-1);
            % If the current frame is a start frame 
            if val_at_curr_frame==1 && val_at_prev_frame==0
                start_frame = frame;
                firing_frames = [firing_frames, start_frame];
                event_number = event_number+1;
            end
        end
    end

    function[num_events_in_window_per_frame] = get_synchronous_firing_bins_v2( synchrony_window_frames )
        % Return a NON-normalized cell array of synchronous firing
        % synchrony_array{2}=number of times two cells fire on the same frame
        
        firing_times = {};
        for cell_i=1:totalActiveCells
            binaryFiring_single_cell = binaryFiring(cell_i,:);
            % firing_times encoded as: 
            %  start_time_1 = firing_times(1)
            firing_times{cell_i} = get_start_firing_frames_single_cell( binaryFiring_single_cell );
        end
        
        num_events_in_window_per_frame = zeros(num_frames,1);
        for frame=2:num_frames
            num_events_in_window = 0;
            for cell_i=1:synchrony_window_frames:totalActiveCells
                firing_times_single_cell = firing_times{ cell_i };
                
                for firing_time_single_cell_i=1:length(firing_times_single_cell)
                    % Find number of frames seperating `frame` and the current firing time
                    firing_time_single_cell = firing_times_single_cell(firing_time_single_cell_i);
                    frame_seperation = abs( firing_time_single_cell-frame );
                    % If the firing time is close enough to the frame, it
                    %     counts as firing in this window
                    if frame_seperation <= synchrony_window_frames
                        num_events_in_window = num_events_in_window + 1;
                    end
                end
            end
            num_events_in_window_per_frame(frame) = num_events_in_window;
        end
    end

    function[num_events_in_window_per_frame] = get_synchronous_firing_bins_v3( synchrony_window_frames )
        % Return a NON-normalized cell array of synchronous firing
        % synchrony_array{2}=number of times two cells fire on the same frame
        
        firing_times = {};
        % firing_synchrony_record contains a correspondance of firing time to number of synchronous cells
        % Example: `firing_synchrony_record{3}==[1,2]` if the third active
        % cell fires twice, once alone, and once synchronized with one other cell
        firing_synchrony_record = {};
        for cell_i=1:totalActiveCells
            binaryFiring_single_cell = binaryFiring(cell_i,:);
            % firing_times encoded as: 
            %  start_time_1 = firing_times(1)
            firing_start_times = get_start_firing_frames_single_cell( binaryFiring_single_cell );
            firing_times{cell_i} = firing_start_times;
            firing_synchrony_record{cell_i} = zeros(length(firing_start_times),1);
        end
        
        num_events_in_window_per_frame = zeros(num_frames,1);
        for frame=2:num_frames
            % First need to find number of synchronous events in the window
            num_events_in_window = 0;
            for cell_i=1:synchrony_window_frames:totalActiveCells
                firing_times_single_cell = firing_times{ cell_i };
                
                for firing_time_single_cell_i=1:length(firing_times_single_cell)
                    % Find number of frames seperating `frame` and the current firing time
                    firing_time_single_cell = firing_times_single_cell(firing_time_single_cell_i);
                    frame_seperation = abs( firing_time_single_cell-frame );
                    % If the firing time is close enough to the frame, it
                    %     counts as firing in this window
                    if frame_seperation <= synchrony_window_frames
                        num_events_in_window = num_events_in_window + 1;
                    end
                end
            end
            num_events_in_window_per_frame(frame) = num_events_in_window;
            
            % Now, record the number of synchronous events in the window in
            % `firing_synchrony_record`
            for cell_i=1:synchrony_window_frames:totalActiveCells
                firing_synchrony_record_single_cell = firing_synchrony_record{ cell_i };
                
                for firing_time_single_cell_i=1:length(firing_times_single_cell)
                    % Find number of frames seperating `frame` and the current firing time
                    firing_time_single_cell = firing_times_single_cell(firing_time_single_cell_i);
                    frame_seperation = abs( firing_time_single_cell-frame );
                    % If the firing time is close enough to the frame, it
                    %     counts as firing in this window
                    if frame_seperation <= synchrony_window_frames
                        current_stored_val = firing_synchrony_record_single_cell{cell_i}(firing_time_single_cell_i);
                        % Store the largest number of synchronous cells found in the sliding window. 
                        % If the cell firing alone in one window, and firing with 2 other cells in another window, then
                        % record as firing with 2 other cells.
                        firing_synchrony_record_single_cell{cell_i}(firing_time_single_cell_i) ...
                            = max( current_stored_val, num_events_in_window);
                    end
                end
            end
        end
    end

% Naive implementation, counts way too many
% synchrony_array = get_synchronous_firing_bins_v2( 6 ); % NOT NORMALIZED
% synchrony_list = [];
% for num_synchronous_cells=1:length(synchrony_array)
%     num_times_occurring = synchrony_array{num_synchronous_cells};
%     for num_times_occurring=1:num_times_occurring
%         synchrony_list = [synchrony_list,num_synchronous_cells];
%     end
% end
% hist( synchrony_list, 5 )


synchrony_list = get_synchronous_firing_bins_v2( 6 ); % NOT NORMALIZED
hist( synchrony_list, 5 )

title('Histogram of synchronized cell activity')
xlabel('Simultaneous event count') 
ylabel('Number of occurences of this many cells simultaneously firing') 
synchrony_mean = mean(synchrony_list);
synchrony_std = std(synchrony_list);
fprintf('__Synchrony__\n')
fprintf('  -  mean:%d    std:%d   \n\n',synchrony_mean,synchrony_std)
pause();
pause();

end