
% Remove residuals
clear all;

% ========================= PATH VARIABLES ===============================================================================================

% Path vars (there are some relative paths used in the script):
PATH_EEGLAB           = '/home/plkn/eeglab2022.1/';
PATH_RAW_DATA         = '/mnt/data_dump/fatigued_blinking/0_raw/';
PATH_COMBINED         = '/mnt/data_dump/fatigued_blinking/1_combined/';
PATH_ICSET            = '/mnt/data_dump/fatigued_blinking/2_icset/';
PATH_EYE_CATCHER      = '/mnt/data_dump/fatigued_blinking/eye-catch-master/';

% The subject list
subject_list =    {'Vp01', 'Vp02', 'Vp03',...
                   'Vp04', 'Vp05', 'Vp06',...
                   'Vp07', 'Vp08', 'Vp09',...
                   'Vp10', 'Vp11', 'Vp12',...
                   'Vp13', 'Vp14', 'Vp15',...
                   'Vp40', 'Vp41', 'Vp42',...
                   'Vp43', 'Vp44', 'Vp45'};

%subject_list =    {'Vp01'};

% Init EEGlab
addpath(PATH_EEGLAB);
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;


addpath(PATH_EYE_CATCHER);

% SWITCH: Switch parts of script on/off
to_execute = {'part3'};

% The new order
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

% ========================= PART 1: Combine blocks  =========================================================================================================
if ismember('part1', to_execute)

    % Find chanlocfile
    channel_location_file = [PATH_RAW_DATA, 'standard-10-5-cap385.elp'];

    % Iterating subject list
    for s = 1 : length(subject_list)

        % Current iteration subject
        subject = subject_list{s};

        % Load blocks
        EEG1 = pop_loadbv(PATH_RAW_DATA, [subject, '_Block1.vhdr'], [], [1 : 60]);
        EEG2 = pop_loadbv(PATH_RAW_DATA, [subject, '_Block2.vhdr'], [], [1 : 60]);
        EEG3 = pop_loadbv(PATH_RAW_DATA, [subject, '_Block3.vhdr'], [], [1 : 60]);

        % Merge blocks
        EEG = pop_mergeset(EEG1, EEG2);
        EEG = pop_mergeset(EEG, EEG3);

        % Resample data
        EEG = pop_resample(EEG, 250);

        % Add chanlocs
        EEG = pop_chanedit(EEG, 'lookup', channel_location_file);
        EEG.chanlocs_original = EEG.chanlocs;

        % Save combined data
        EEG = pop_saveset(EEG, 'filename', [subject '_icset.set'], 'filepath', PATH_COMBINED, 'check', 'on');

    end
end

% ========================= PART 2: Calculate blink ICs =========================================================================================================
if ismember('part2', to_execute)

    % Iterating subject list
    for s = 1 : length(subject_list)

        % Current iteration subject
        subject = subject_list{s};
        id = str2num(subject(3 : 4));

        % Load data
        EEG = pop_loadset('filename', [subject '_icset.set'], 'filepath', PATH_COMBINED, 'loadmode', 'all');

        % Init event struct
        new_events = struct('latency', {},...
                            'oldtype', {},...
                            'type', {},...
                            'code', {},...
                            'id', {},...
                            'stimloc', {},... 
                            'corresp', {},...   
                            'saliency', {},...
                            'bl', {},...  
                            'sbl', {},...      
                            'sbl_total', {},...    
                            'stimnum', {},...     
                            'stimnum_sbl', {},...
                            'resploc', {},...    
                            'acc', {},...     
                            'rt', {},...         
                            'urevent', {},...
                            'duration', {}...
                            );

        % Sort events by latency
        [~, idx] = sort(cell2mat({EEG.event.latency}));
        EEG.event = EEG.event(idx);

        % Code experimental conditions
        stimnum = 0;
        stimnum_sbl = 0;
        sblnum = 0;
        for e = 1 : length(EEG.event)

            % If stim event
            if strcmpi(EEG.event(e).type(1), 'S') & ismember(str2num(EEG.event(e).type(2 : end)), [1 : 108])

                % Get eventcode and increase eventcount 
                enum = str2num(EEG.event(e).type(2 : end));
                stimnum = stimnum + 1;

                % Track change of sbl and increase sbl stimcount
                if ceil(enum / 12) > sblnum
                    sblnum = sblnum + 1;
                    stimnum_sbl = 0;
                end
                stimnum_sbl = stimnum_sbl + 1;

                % Code all the things
                new_events(stimnum).latency = EEG.event(e).latency;
                new_events(stimnum).oldtype = EEG.event(e).type;
                new_events(stimnum).type = 'stim';
                new_events(stimnum).code = 'stim';
                new_events(stimnum).id = id;
                if ismember(mod(enum, 4), [1, 0])
                    new_events(stimnum).stimloc = 'left';
                else
                    new_events(stimnum).stimloc = 'right';
                end
                if ismember(mod(enum, 4), [1, 2])
                    new_events(stimnum).corresp = 1; % corresponding
                else
                    new_events(stimnum).corresp = 2; % non-corresponding
                end
                if ismember(mod(enum, 12), [1, 2, 3, 4])
                    new_events(stimnum).saliency = 'low';
                elseif ismember(mod(enum, 12), [5, 6, 7, 8])
                    new_events(stimnum).saliency = 'mid';
                else
                    new_events(stimnum).saliency = 'high';
                end
                new_events(stimnum).bl = ceil(enum / 36);
                new_events(stimnum).sbl = ceil(ceil(enum / 12) / 3);
                new_events(stimnum).sbl_total = ceil(enum / 12);
                new_events(stimnum).stimnum = stimnum;
                new_events(stimnum).stimnum_sbl = stimnum_sbl;

                % Loop for response
                f = e + 1;
                stimcodes = {};
                for n = 1 : 108
                    if n < 10
                        stimcodes{end + 1} = ['S  ' num2str(n)];
                    elseif n < 100
                        stimcodes{end + 1} = ['S ' num2str(n)];
                    else
                        stimcodes{end + 1} = ['S' num2str(n)];
                    end
                end
                while ~strcmpi(EEG.event(f).type(1), 'R') &...
                      ~ismember(EEG.event(f).type, stimcodes) &...
                      f < length(EEG.event)
                    f = f + 1;
                end
                rnum = 0;
                if strcmpi(EEG.event(f).type(1), 'R')
                    rnum = str2num(EEG.event(f).type(2 : end));
                end

                if ismember(rnum, [3, 13])
                    new_events(stimnum).resploc = 'left';
                    new_events(stimnum).acc = 1;
                    new_events(stimnum).rt = (EEG.event(f).latency - EEG.event(e).latency) * 1000 / EEG.srate;
                elseif ismember(rnum, [4, 14])
                    new_events(stimnum).resploc = 'right';
                    new_events(stimnum).acc = 1;
                    new_events(stimnum).rt = (EEG.event(f).latency - EEG.event(e).latency) * 1000 / EEG.srate;
                elseif ismember(rnum, [5, 6, 9, 10, 11, 12, 15, 16])
                    new_events(stimnum).resploc = 'none';
                    new_events(stimnum).acc = 0;
                    new_events(stimnum).rt = NaN;
                elseif ismember(rnum, [7, 8])
                    new_events(stimnum).resploc = 'none';
                    new_events(stimnum).acc = 2;
                    new_events(stimnum).rt = NaN;
                else
                    new_events(stimnum).resploc = 'none';
                    new_events(stimnum).acc = NaN;
                    new_events(stimnum).rt = NaN;
                end

                new_events(stimnum).urevent = e;
                new_events(stimnum).duration = 1;

            end % End event check
        end % End eventit

        % Collect boundaries
        for e = 1 : length(EEG.event)
            if strcmpi(EEG.event(e).type, 'boundary')
                new_events(end + 1).latency = EEG.event(e).latency;
                new_events(end).oldtype = 'boundary';
                new_events(end).type = 'boundary';
                new_events(end).code = 'boundary';
                new_events(end).id = 0;
                new_events(end).stimloc = 'none';
                new_events(end).corresp = 0;
                new_events(end).saliency = 'none';
                new_events(end).bl = NaN;
                new_events(end).sbl = NaN;
                new_events(end).sbl_total = NaN;
                new_events(end).stimnum = NaN;
                new_events(end).stimnum_sbl = NaN;
                new_events(end).resploc = 'none';
                new_events(end).acc = NaN;
                new_events(end).rt = NaN;
                new_events(end).urevent = length(new_events);
                new_events(end).duration = 1;
            end
        end

        % Sort new events by latency
        [~, idx] = sort(cell2mat({new_events.latency}));
        new_events = new_events(idx);

        % Replace events by new events
        EEG.event = new_events;
        EEG = eeg_checkset(EEG, 'eventconsistency');

        % Save a backup before filtering
        BACKUP = EEG;

        % Bandpass filter data (ERPlab toolbox function) 
        EEG = pop_basicfilter(EEG, [1 : EEG.nbchan], 'Cutoff', [2, 20], 'Design', 'butter', 'Filter', 'bandpass', 'Order', 4, 'RemoveDC', 'on', 'Boundary', 'boundary');

        % Bad channel detection
        [EEG, i1] = pop_rejchan(EEG, 'elec', [1 : EEG.nbchan], 'threshold', 10, 'norm', 'on', 'measure', 'kurt');
        [EEG, i2] = pop_rejchan(EEG, 'elec', [1 : EEG.nbchan], 'threshold', 5, 'norm', 'on', 'measure', 'prob');
        EEG.chans_rejected = horzcat(i1, i2);
        EEG.n_chans_rejected = length(horzcat(i1, i2));

        % Interpolate channels
        EEG = pop_interp(EEG, EEG.chanlocs_original, 'spherical');

        % Reref common average
        EEG = pop_reref(EEG, []);

        % Determine rank of data
        dataRank = sum(eig(cov(double(EEG.data'))) > 1e-6);

        % Resample data
        EEG = pop_resample(EEG, 125);

        % Epoch data 150 ms stim + jittered 2800ms ISI
        EEG = pop_epoch(EEG, {'stim'}, [-2, 2], 'newname', [num2str(id) '_seg'], 'epochinfo', 'yes');
        
        % Autoreject data before ICA
        EEG.segs_original_n = size(EEG.data, 3);
        [EEG, rejsegs] = pop_autorej(EEG, 'nogui', 'on', 'threshold', 1000, 'startprob', 5, 'maxrej', 5, 'eegplot', 'off');
        EEG.segs_rejected_before_ica = length(rejsegs);

        % Run ICA on randsample of 666 trials
        idx = randsample([1 : size(EEG.data, 3)], 600);
        EEG = pop_selectevent(EEG, 'epoch', idx, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');

        % Runica & ICLabel
        EEG = pop_runica(EEG, 'extended', 1, 'interrupt', 'on', 'PCA', dataRank);
        EEG = iclabel(EEG);

        % Copy ICs to backup set
        BACKUP = pop_editset(BACKUP, 'icachansind', 'EEG.icachansind', 'icaweights', 'EEG.icaweights', 'icasphere', 'EEG.icasphere');
        BACKUP.etc = EEG.etc;

        % Save IC set
        BACKUP = pop_saveset(BACKUP, 'filename', [subject '_icset.set'], 'filepath', PATH_ICSET, 'check', 'on');

    end % End subject loop

end% End part2

% ========================= PART 3: Create blink events =========================================================================================================
if ismember('part3', to_execute)

    % Iterating subject list
    for s = 1 : length(subject_list)

        % Current iteration subject
        subject = subject_list{s};
        id = str2num(subject(3 : 4));

        % Load IC-data
        EEG = pop_loadset('filename', [subject '_icset.set'], 'filepath', PATH_ICSET, 'loadmode', 'all');

        % Bandpass filter data, mainly to get rid of slow drifts
        EEG = pop_basicfilter(EEG, [1 : EEG.nbchan], 'Cutoff', [0.1, 30], 'Design', 'butter', 'Filter', 'bandpass', 'Order', 4, 'RemoveDC', 'on', 'Boundary', 'boundary');

        % Get eye-IC indices
        eye_ICs_idx = find(EEG.etc.ic_classification.ICLabel.classifications(:, 3) >= 0.8);

        % Get frontal channel indices
        frontal_channel_idx_left  = ismember({EEG.chanlocs.labels}, {'Fp1','F3','F7'});
        frontal_channel_idx_right = ismember({EEG.chanlocs.labels}, {'Fp2','F4','F8'});

        % Correlate left and right projections of eye-ICs in sensor space
        mask_positive_correlation = [];
        for ic = 1 : length(eye_ICs_idx)

            % Get IC activation
            ICA_ACT = pop_subcomp(EEG, eye_ICs_idx(ic), 0, 1);

            % Correlate the average left and right sensor projections
            rr = corrcoef(mean(ICA_ACT.data(frontal_channel_idx_left, :))', mean(ICA_ACT.data(frontal_channel_idx_right, :))');

            % Set mask to separate blink-ICs from saccade ICs
            if rr > 0
                mask_positive_correlation(ic) = 1;
            else
                mask_positive_correlation(ic) = 0;
            end

        end

        % Project blink activation to sensor space
        EEG = pop_subcomp(EEG, eye_ICs_idx(find(mask_positive_correlation)), 0, 1);

        % Run blinker
        params = [];
        params.verbose = false;
        params.showMaxDistribution = false;
        [~, ~, blinks] = pop_blinker(EEG, params);

        % Get blink latencies at afz
        idx_afz = ismember({blinks.signalData.signalLabel}, {'afz'});
        blink_latencies = blinks.signalData(idx_afz).blinkPositions';

        % Create blink events at maximum
        for e = 1 : size(blink_latencies, 1)

            % Get blink peak latency at afz
            tmp = blinks.signalData(idx_afz).signal;
            tmp = tmp(blink_latencies(e, 1) : blink_latencies(e, 2));
            [~, idx_max]= max(tmp);
            idx_max = blink_latencies(e, 1) + idx_max;

            % Create event
            event_idx = length(EEG.event) + 1;
            EEG.event(event_idx).latency = idx_max;
            EEG.event(event_idx).type = 'blink';
            EEG.event(event_idx).code = 'blink';

        end

        % Fix latencies
        EEG = eeg_checkset(EEG, 'eventconsistency');

        % Save events
        blink_events = EEG.event;

        % Re-load continuous data
        EEG = pop_loadset('filename', [subject '_icset.set'], 'filepath', PATH_ICSET, 'loadmode', 'all');

        % Replace events
        EEG.event = blink_events;

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
            if strcmpi(EEG.event(e).type, 'blink') & subblock > 0

                % Set latencies
                lat_all = EEG.event(e).latency;
                lat_block = EEG.event(e).latency - block_lat_offset;
                lat_subblock = EEG.event(e).latency - subblock_lat_offset;
                lat_event = EEG.event(e).latency - event_lat_offset;

                % Select blink
                blink_count = blink_count + 1;
                EEG.event(e).type = 'selected';

                % Increase post-stim blink counter
                post_stim_blink_nr = post_stim_blink_nr + 1;

                % Collect time on task predictors for blink events
                EEG.trialinfo(blink_count, :) = [id, blink_count, block, subblock, lat_all, lat_block, lat_subblock, lat_event, post_stim_blink_nr];

            end

        end % End event loop

        % Epoch data
        [EEG, epoch_idx] = pop_epoch(EEG, {'selected'}, [-0.5, 1.2], 'newname', 'blink', 'epochinfo', 'yes');
        EEG = pop_rmbase(EEG, [-300, -100]);
        EEG.trialinfo = EEG.trialinfo(epoch_idx, :);

        aa = bb;




    end % End subject loop

end% End part3