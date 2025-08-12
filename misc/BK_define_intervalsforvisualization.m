function [interval_start,interval_end]=BK_define_intervalsforvisualization(data,interval_size,numberoftopintervals)
% Step 1: Define the array (example data)
% data = rand(100, 1) * 100; % 100 random numbers between 0 and 100
% 
% Step 2: Define the interval size
% interval_size = 10; % Size of each interval

% Step 3: Calculate the number of intervals
min_val = min(data); % Minimum value in the data
max_val = max(data); % Maximum value in the data

% Determine the range of values (the intervals should cover this range)
num_intervals = floor((max_val - min_val) / interval_size);

% Step 4: Initialize a vector to store the count of elements in each interval
interval_counts = zeros(num_intervals, 1);

% Step 5: Count how many values fall into each interval
for i = 1:length(data)
    % Find which interval the current value belongs to
    interval_index = floor((data(i) - min_val) / interval_size) + 1;
    
    % Make sure the index is within the bounds
    if interval_index >= 1 && interval_index <= num_intervals
        interval_counts(interval_index) = interval_counts(interval_index) + 1;
    end
end

% Step 6: Find the top 4 intervals with the most elements
[~, sorted_indices] = sort(interval_counts, 'descend');
top_intervals = sorted_indices(1:numberoftopintervals); % Indices of the top 4 intervals

% Step 7: Display the results
fprintf('Top 4 intervals with the most elements:\n');
for i = 1:4
    interval_start(i) = min_val + (top_intervals(i) - 1) * interval_size;
    interval_end(i) = interval_start(i) + interval_size;
    fprintf('Interval %d: %.2f to %.2f, Count: %d\n', i, interval_start(i), interval_end(i), interval_counts(top_intervals(i)));
end
