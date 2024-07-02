% Define the selected trials vector
%selected_trials = [1, 0, 1, ..., 0];  % example vector

% Extract selected trials
%selected_idx = find(selected_trials == 1);
%EEG_selected = pop_select(EEG, 'trial', selected_idx);

% Compute the ERP
ERP = mean(EEG.data, 3);  % averaging across the third dimension (trials)
%ERP = mean(EEG_selected.data, 3);  % averaging across the third dimension (trials)

% Plot the ERP for a specific channel
time_vector = EEG.times;  % time points from the EEG structure
channel_to_plot = 19;  % channel number to plot

figure;
plot(time_vector, ERP(channel_to_plot, :));
xlabel('Time (ms)');
ylabel('Amplitude (µV)');
title(['ERP for Channel ', num2str(channel_to_plot)]);
% Add vertical lines at specified time points
hold on;
xline([-1300 -1200, -600, 0, 800, 1500,2000], '--r', {'fix point', 'baseline on', 'cue change', 'evidence on', 'response allowed', 'DDL-1.5s','DDL-2s'});
hold off;


% % Plot scalp topography at 300 ms
% [~, time_300ms_idx] = min(abs(time_vector - 300));
% 
% figure;
% topoplot(ERP(:, time_300ms_idx), EEG.chanlocs, 'maplimits', [-max(abs(ERP(:))), max(abs(ERP(:)))]);
% colorbar;
% title('Scalp Topography at 300 ms');
% Define time points for scalp topographies (from -1.5 to 3 seconds)

num_plots = 10;
time_points = linspace(-1200, 2000, num_plots);  % in milliseconds

% Find the indices of the defined time points
time_indices = arrayfun(@(t) find(EEG.times >= t, 1), time_points);

% Plot scalp topographies
figure;
for i = 1:num_plots
    subplot(2, 5, i);  % create a 2x5 subplot
    topoplot(ERP(:, time_indices(i)), EEG.chanlocs, 'maplimits', [-max(abs(ERP(:))), max(abs(ERP(:)))]);
    title([num2str(time_points(i)), ' ms']);
    colorbar;
end