
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
subject_list = setdiff(subject_list, {'Vp01'});

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
    EEG.trialinfo = [];
    for e = 1 : length(EEG.event)

        % If stim event
        if strcmpi(EEG.event(e).type(1), 'S') & ismember(str2num(EEG.event(e).type(2 : end)), [1 : 108])

            % Reset post-stim blink counter
            post_stim_blink_nr = 0;

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
        if strcmpi(EEG.event(e).type, 'selected') & subblock > 0

            % Set latencies
            lat_all = EEG.event(e).latency;
            lat_block = EEG.event(e).latency - block_lat_offset;
            lat_subblock = EEG.event(e).latency - subblock_lat_offset;
            lat_event = EEG.event(e).latency - event_lat_offset;

            % Select blink
            blink_count = blink_count + 1;
            EEG.event(e).type = 'selected_again';

            % Increase post-stim blink counter
            post_stim_blink_nr = post_stim_blink_nr + 1;

            % Collect time on task predictors for blink events
            EEG.trialinfo(blink_count, :) = [id, blink_count, block, subblock, lat_all, lat_block, lat_subblock, lat_event, post_stim_blink_nr];

        end

    end % End event loop

    % Epoch data
    [EEG, epoch_idx] = pop_epoch(EEG, {'selected_again'}, [-0.5, 1.2], 'newname', 'blink', 'epochinfo', 'yes');
    EEG = pop_rmbase(EEG, [-300, -100]);
    EEG.trialinfo = EEG.trialinfo(epoch_idx, :);

    % Compute the kernel density estimation for blink latencies
    [blink_distributions(s, :), blink_distribution_x] = ksdensity(EEG.trialinfo(:, 8));
    
    % Exclude blinks that are too close to stimulus (response locked...)
    %to_keep = find(EEG.trialinfo(:, 8) >= 100 & EEG.trialinfo(:, 8) <= 500 & EEG.trialinfo(:, 9) > 1);
    to_keep = find(EEG.trialinfo(:, 9) > 1);
    EEG = pop_select(EEG, 'trial', to_keep);
    EEG.trialinfo = EEG.trialinfo(to_keep, :);

    % Get number of remaining blinks
    remaining_blinks(s) = length(to_keep);

    % Calculate erps
    for sbl = 1 : 9
        idx_sbl = EEG.trialinfo(:, 4) == sbl;
        erps(s, sbl, :, :) = squeeze(mean(EEG.data(:, :, idx_sbl), 3));

        % Get number of remaining blinks
        remaining_blinks(s, sbl) = sum(idx_sbl);

    end

    % Regression design matrix
    desmat = [ones(size(EEG.trialinfo, 1), 1), EEG.trialinfo(:, 5), EEG.trialinfo(:, 6)];

    % Scale predictors
    desmat(:, 2) = desmat(:, 2) / max(abs(desmat(:, 2)));
    desmat(:, 3) = desmat(:, 3) / max(abs(desmat(:, 3)));

    % regression
    for ch = 1 : EEG.nbchan

        % Get channel data with trials as rows
        d = squeeze(EEG.data(ch, :, :))';

        % OLS fit
        tmp = (desmat' * desmat) \ desmat' * d;
        coefs_1(s, ch, :) = squeeze(tmp(2, :));
        coefs_2(s, ch, :) = squeeze(tmp(3, :));

    end % End chanit

end % End subjecct loop

% Plot erps
figure()
idx_chan = 14;
for s = 1 : length(subject_list)

    subplot(4, 5, s)
    for sbl = 1 : 9
        pd = squeeze(erps(s, sbl, idx_chan, :));
        plot(EEG.times, pd);
        title(num2str(s))
        if sbl ~= 9
            hold on
        end
        if s == length(subject_list)
            legend({'1','2','3','4','5','6','7','8','9'})
        end
    end
end

figure()
idx_chans = [6, 14, 24, 34, 44, 59];
for ch = 1 : length(idx_chans)
    subplot(2, 3, ch)
    for sbl = 1 : 9
        pd = mean(squeeze(erps(:, sbl, idx_chans(ch), :)), 1);
        plot(EEG.times, pd);
        title(new_order_labels(idx_chans(ch)))
        if sbl ~= 9
            hold on
        end
    end
    if ch == length(idx_chans)
        legend({'1','2','3','4','5','6','7','8','9'})
    end
end




% Plot distributions
figure()
idx_chan = 14;
for s = 1 : length(subject_list)
    subplot(4, 5, s)
    pd = blink_distributions(s, :);
    plot(blink_distribution_x, pd);
    title(num2str(s))
    if s == length(subject_list)
        legend({'1','2','3','4','5','6','7','8','9'})
    end
end

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

% GA struct coefs
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

% GA struct null
GA = {};
for s = 1 : length(subject_list)
    ga_template.avg = zeros(size(squeeze(coefs_1(s, :, :))));
    GA{s} = ga_template;
end 
GA_null = ft_timelockgrandaverage(cfg, GA{1, :});

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
[stat_all]  = ft_timelockstatistics(cfg, GA_1, GA_null);
[stat_blk]  = ft_timelockstatistics(cfg, GA_2, GA_null);

% Plot coefs
figure()

subplot(2, 2, 1)
pd = squeeze(mean(coefs_1, 1));
contourf(EEG.times,[1 : 60], pd, 50, 'linecolor','none')
hold on
contour(EEG.times,[1 : 60], stat_all.mask, 1, 'linecolor', 'k', 'LineWidth', 2)
colormap('jet')
clim([-2, 2])
title('whole experiment')

subplot(2, 2, 2)
pd = squeeze(mean(coefs_2, 1));
contourf(EEG.times,[1 : 60], pd, 50, 'linecolor','none')
hold on
contour(EEG.times,[1 : 60], stat_blk.mask, 1, 'linecolor', 'k', 'LineWidth', 2)
colormap('jet')
clim([-2, 2])
title('block')


% Plot coef topos at selected time points
figure()
tpoints = [100, 225, 300, 400, 500];
for t = 1 : length(tpoints)
    subplot(1, length(tpoints), t)
    tidx = EEG.times >= tpoints(t) - 25 & EEG.times <= tpoints(t) + 25;
    tmp = squeeze(mean(coefs_2, 1));
    pd = mean(tmp(:, tidx), 2);
    topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'on');
    title(num2str(tpoints(t)))
    colormap(jet);
end







