
% Remove residuals
clear all;

% Path vars
PATH_EEGLAB = '/home/plkn/eeglab2022.1/';
PATH_DATA = '/mnt/data_fast/fatigued_blinking/0_blink_events_continuous/';


% The subject list
subject_list =    {'Vp01', 'Vp02', 'Vp03',...
                   'Vp04', 'Vp05', 'Vp06',...
                   'Vp08', 'Vp09', 'Vp10',...
                   'Vp11', 'Vp12', 'Vp13',...
                   'Vp14', 'Vp40', 'Vp41',...
                   'Vp43', 'Vp45'};
% Exclude subjects
subject_list = setdiff(subject_list, {'Vp12', 'Vp43'});

% Init eeglab
addpath(PATH_EEGLAB);
eeglab;

% Init fieldtrip
ft_path = '/home/plkn/fieldtrip-master/';
addpath(ft_path);
ft_defaults;

% Subject 
remaining_blinks = [];
erps = [];
blink_distributions = [];
coefs_true = [];
coefs_fake = [];

% Iterating subject list
for s = 1 : length(subject_list)

    % Current iteration subject
    subject = subject_list{s};
    id = str2num(subject(3 : 4));

    % Load data
    EEG = pop_loadset('filename', ['BR_', subject, '.set'], 'filepath', PATH_DATA, 'loadmode', 'all');

    % The order of things
    new_order_labels = {...
    'Fp1',...
    'Fp2',...
    'AFz',...
    'F7',...
    'F3',...
    'Fz',...
    'F4',...
    'F8',...
    'FT9',...
    'FT7',...
    'FC5',...
    'FC3',...
    'FC1',...
    'FCz',...
    'FC2',...
    'FC4',...
    'FC6',...
    'FT8',...
    'FT10',...
    'T7',...
    'C5',...
    'C3',...
    'C1',...
    'Cz',...
    'C2',...
    'C4',...
    'C6',...
    'T8',...
    'TP9',...
    'TP7',...
    'CP5',...
    'CP3',...
    'CP1',...
    'CPz',...
    'CP2',...
    'CP4',...
    'CP6',...
    'TP8',...
    'TP10',...
    'P7',...
    'P5',...
    'P3',...
    'P1',...
    'Pz',...
    'P2',...
    'P4',...
    'P6',...
    'P8',...
    'PO9',...
    'PO7',...
    'PO5',...
    'PO1',...
    'POz',...
    'PO2',...
    'PO6',...
    'PO8',...
    'PO10',...
    'O1',...
    'Oz',...
    'O2',...
    };

    % Get original chanloc labels
    chanlocs_labels = {};
    for ch = 1 : length(EEG.chanlocs)
        chanlocs_labels{end + 1} = EEG.chanlocs(ch).labels;
    end

    % Get new order indices
    new_order_idx = [];
    for ch = 1 : length(EEG.chanlocs)
        new_order_idx(end + 1) = find(strcmpi(new_order_labels{ch}, chanlocs_labels));
    end

    % Re-order chanlocs
    EEG.chanlocs = EEG.chanlocs(new_order_idx);

    % Re-order data
    EEG.data = EEG.data(new_order_idx, :, :);

    % Loop events
    block = 0;
    subblock = 0;

    % Loop events to select blinks
    blink_count = 0;
    stim_count = 0;
    EEG.trialinfo = [];
    for e = 1 : length(EEG.event)

        % If stim event
        if strcmpi(EEG.event(e).type(1), 'S') & ismember(str2num(EEG.event(e).type(2 : end)), [1 : 108])

            % Count stimuli
            stim_count = stim_count + 1;

            % Get event number
            enum = str2num(EEG.event(e).type(2 : end));

            % If block changes
            if ceil(enum / 36) > block
                
                % change block
                block = ceil(enum / 36);

                % Set block offset latency
                block_lat_offset = EEG.event(e).latency;

            end

            % If subblock changes
            if ceil(enum / 12) > subblock

                % change subblock
                subblock = ceil(enum / 12);

                % Set subblock offset latency
                subblock_lat_offset = EEG.event(e).latency;

            end

            % Set event offset latency
            event_lat_offset = EEG.event(e).latency;

        end % End stim-event check

        % If it is a blink in a subblock
        if strcmpi(EEG.event(e).type, 'max') & subblock > 0

            % Set latencies
            lat_all = EEG.event(e).latency;
            lat_block = EEG.event(e).latency - block_lat_offset;
            lat_subblock = EEG.event(e).latency - subblock_lat_offset;
            lat_event = EEG.event(e).latency - event_lat_offset;

            % Select blink
            blink_count = blink_count + 1;
            EEG.event(e).type = 'selected';

            % Collect time on task predictors for blink events
            EEG.trialinfo(blink_count, :) = [id, blink_count, block, subblock, lat_all, lat_block, lat_subblock, lat_event];

        end

    end % End event loop

    % Epoch data
    [EEG, epoch_idx] = pop_epoch(EEG, {'selected'}, [-0.5, 1.2], 'newname', 'blink', 'epochinfo', 'yes');
    EEG = pop_rmbase(EEG, [-300, -100]);
    EEG.trialinfo = EEG.trialinfo(epoch_idx, :);

    % Exclude blinks that are too close to stimulus (response locked...)
    to_keep = find(EEG.trialinfo(:, 8) >= 250 & EEG.trialinfo(:, 8) <= 500);
    EEG = pop_select(EEG, 'trial', to_keep);
    EEG.trialinfo = EEG.trialinfo(to_keep, :);

    % Compute the kernel density estimation for blink latencies
    [blink_distributions(s, :), blink_distribution_x] = ksdensity(EEG.trialinfo(:, 8));

    % Get number of remaining blinks
    remaining_blinks(s) = length(to_keep);

    % Calculate erps
    erps(s, :, :) = squeeze(mean(EEG.data, 3));

    % Regression design matrix
    desmat = [ones(size(EEG.trialinfo, 1), 1), EEG.trialinfo(:, 5), EEG.trialinfo(:, 6), EEG.trialinfo(:, 7)];

    % Scale predictors
    desmat(:, 2) = desmat(:, 2) / max(abs(desmat(:, 2)));
    desmat(:, 3) = desmat(:, 3) / max(abs(desmat(:, 3)));
    desmat(:, 4) = desmat(:, 4) / max(abs(desmat(:, 4)));

    % regression
    for ch = 1 : EEG.nbchan

        % Get channel data with trials as rows
        d = squeeze(EEG.data(ch, :, :))';

        % OLS fit
        tmp = (desmat' * desmat) \ desmat' * d;
        coefs_1(s, ch, :) = squeeze(tmp(2, :));
        coefs_2(s, ch, :) = squeeze(tmp(3, :));
        coefs_3(s, ch, :) = squeeze(tmp(4, :));


    end % End chanit

end % End subjecct loop

% Restructure coordinates
chanlocs = EEG.chanlocs;
chanlabs = {};
coords = [];
for c = 1 : numel(chanlocs)
    chanlabs{c} = chanlocs(c).labels;
    coords(c, :) = [chanlocs(c).X, chanlocs(c).Y, chanlocs(c).Z];
end

% A sensor struct
sensors = struct();
sensors.label = chanlabs;
sensors.chanpos = coords;
sensors.elecpos = coords;

% Prepare neighbor struct
cfg                 = [];
cfg.elec            = sensors;
cfg.feedback        = 'no';
cfg.method          = 'triangulation'; 
neighbours          = ft_prepare_neighbours(cfg);

% A template for GA structs
cfg = [];
cfg.keepindividual = 'yes';
ga_template = [];
ga_template.dimord = 'chan_time';
ga_template.label = chanlabs;
ga_template.time = EEG.times;

% GA struct true
GA = {};
for s = 1 : length(subject_list)
    ga_template.avg = squeeze(coefs_1(s, :, :));
    GA{s} = ga_template;
end 
GA_1 = ft_timelockgrandaverage(cfg, GA{1, :});
GA = {};
for s = 1 : length(subject_list)
    ga_template.avg = squeeze(coefs_2(s, :, :));
    GA{s} = ga_template;
end 
GA_2 = ft_timelockgrandaverage(cfg, GA{1, :});
GA = {};
for s = 1 : length(subject_list)
    ga_template.avg = squeeze(coefs_3(s, :, :));
    GA{s} = ga_template;
end 
GA_3 = ft_timelockgrandaverage(cfg, GA{1, :});

% GA struct fake
GA = {};
for s = 1 : length(subject_list)
    ga_template.avg = zeros(size(squeeze(coefs_1(s, :, :))));
    GA{s} = ga_template;
end 
GA_fake = ft_timelockgrandaverage(cfg, GA{1, :});

% Plot coefs
figure()

subplot(1, 3, 1)
pd = squeeze(mean(coefs_1, 1));
contourf(EEG.times,[1 : 60], pd, 50, 'linecolor','none')
%hold on
%contour(EEG.times,[1 : 60], stat.mask, 1, 'linecolor', 'k', 'LineWidth', 2)
colormap('jet')
clim([-2, 2])

subplot(1, 3, 2)
pd = squeeze(mean(coefs_2, 1));
contourf(EEG.times,[1 : 60], pd, 50, 'linecolor','none')
%hold on
%contour(EEG.times,[1 : 60], stat.mask, 1, 'linecolor', 'k', 'LineWidth', 2)
colormap('jet')
clim([-2, 2])

subplot(1, 3, 3)
pd = squeeze(mean(coefs_3, 1));
contourf(EEG.times,[1 : 60], pd, 50, 'linecolor','none')
%hold on
%contour(EEG.times,[1 : 60], stat.mask, 1, 'linecolor', 'k', 'LineWidth', 2)
colormap('jet')
clim([-2, 2])

aa = bb

% Testparams
testalpha  = 0.05;
voxelalpha  = 0.01;
nperm = 1000;

% Set config
cfg = [];
cfg.tail             = 1;
cfg.statistic        = 'depsamplesFmultivariate';
cfg.alpha            = testalpha;
cfg.neighbours       = neighbours;
cfg.minnbchan        = 2;
cfg.method           = 'montecarlo';
cfg.correctm         = 'cluster';
cfg.clustertail      = 1;
cfg.clusteralpha     = voxelalpha;
cfg.clusterstatistic = 'maxsum';
cfg.numrandomization = nperm;
cfg.computecritval   = 'yes'; 
cfg.ivar             = 1;
cfg.uvar             = 2;

% Set up design
n_subjects = length(subject_list);
design = zeros(2, n_subjects * 2);
design(1, :) = [ones(1, n_subjects), 2 * ones(1, n_subjects)];
design(2, :) = [1 : n_subjects, 1 : n_subjects];
cfg.design = design;

% The tests
[stat]  = ft_timelockstatistics(cfg, GA_1, GA_fake);




% Plot erps at Pz
figure();
subplot(2,1,1)
pd = squeeze(erps(:, 44, :))';
plot(EEG.times, pd);
subplot(2,1,2)
plot(blink_distribution_x, blink_distributions);
