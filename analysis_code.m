%% load
% clear all;
load ('SpikesX10U12D.mat');
SpikesFullData = SpikesX10U12D; clear SpikesX10U12D; 

%% test
TrialDuration = 1280/1000; % s;
TimeBin = 20/1000; % s

[N_units, N_dirct, N_rept] = size (SpikesFullData);
N_bins = floor(TrialDuration / TimeBin);
V_bins = 0:TimeBin:TrialDuration;

Directions = pi/6*(0:N_dirct-1);
DirDegrees = rad2deg(Directions);

SpikesImages = zeros (N_units, N_bins, N_dirct, N_rept); 
SpikesCounts = zeros (N_units, N_dirct, N_rept); % This variable will be used in part 2

% Run over each unit, direction and repetition.
for i_unit = 1:N_units
	for i_dirc = 1:N_dirct
		for i_rept = 1:N_rept
            % Split the spikes of the trial into the time bins
			SpikesTimes = SpikesFullData(i_unit, i_dirc, i_rept).TimeList;
            psth = histcounts(SpikesTimes.', V_bins); 
            for i_bins = 1:N_bins
                SpikesImages(i_unit, i_bins, i_dirc, i_rept) = psth(i_bins);
            end
            % Count the total spikes of the repetition
            SpikesCounts(i_unit, i_dirc, i_rept) = sum(psth);
        end
    end
end

DisplayedUnit = 4;

% extract the rates of each direction and time-bin for the displayed unit.
Rate = sum(squeeze(SpikesImages(DisplayedUnit,:,:,:)),3);
% Fix the unit - spikes/sec
Rate = Rate/TimeBin/N_rept; 
% Plot the selected unit's rate
waterfall(V_bins(1:64),DirDegrees,Rate');
title("PSTH Signals Per Direction - Unit 4");
xlabel('Time(sec)'); ylabel('Direction (deg)'); zlabel('Rate (Spikes/Sec)');
xlim([V_bins(1), V_bins(end)]);
ylim([DirDegrees(1), DirDegrees(end)]);


%% Part 2

% The unit numbers of the direction specific units
DirectionUnits = [4, 5];

figure('Color', 'w', 'Units', 'centimeters', 'Position', [0 0 25 10]); hold on;

ResponseM = zeros(1,N_dirct);
ResponseSD = zeros(1,N_dirct);

% Calculate the mean and std of total spikes, for each unit and direction
for unit_idx = 1:10
    for i_dirc = 1:N_dirct
        ResponseM(i_dirc) = mean(SpikesCounts(unit_idx, i_dirc,:));
        ResponseSD(i_dirc) = std(SpikesCounts(unit_idx, i_dirc,:));
    end
    
    subplot (3,5, unit_idx); hold on; % create subplot on the main figure
    xlabel('Direction (deg)','fontsize',12);
    ylabel('Mean Spikes (Spikes/Sec)','fontsize',12);
    title ("Unit " + unit_idx);
    
    % plot the experimental data:
    ExData = errorbar(DirDegrees, ResponseM, ResponseSD, 'o');
    % fit the axis to the data dimensions
    DataMin = min(ResponseM - ResponseSD) - 1;
    DataMax = max(ResponseM + ResponseSD) + 1;
    ylim([DataMin,DataMax]);
    
    % plot the fitted curve :
    if ismember(unit_idx,DirectionUnits)
        % Use direction selective von mises
        VM = fit (Directions', ResponseM', fittype ('A * exp (k * cos (x - PO) )', ...
        'coefficients', {'A', 'k', 'PO'}, 'independent', 'x'));
    else
        % Use orientation selective von mises
        VM = fit (Directions', ResponseM', fittype ('A * exp (k * cos (2*(x - PO)) )', ...
        'coefficients', {'A', 'k', 'PO'}, 'independent', 'x'));
    end
    % Add the fitted von mises to the graph
    VonMises =  plot (DirDegrees, VM(Directions), 'r');
    set (gca, 'FontSize', 13);
    xlim([DirDegrees(1), DirDegrees(end)]);
    
    leg = subplot (3,5,13);  
    legend(leg,[ExData,VonMises], "Mean Spikes","Von Mises");
    axis (leg, 'off');
    

end
