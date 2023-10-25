
% Remove residuals
clear all;

% ========================= PATH VARIABLES ===============================================================================================

% Path vars (there are some relative paths used in the script):
PATH_EEGLAB           = '/home/plkn/eeglab2022.1/';
PATH_RAW_DATA         = '/mnt/data_dump/fatigued_blinking/0_raw/';
PATH_COMBINED         = '/mnt/data_dump/fatigued_blinking/1_combined/';
PATH_ICSET            = '/mnt/data_dump/fatigued_blinking/2_icset/';

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

% SWITCH: Switch parts of script on/off
to_execute = {'part2'};

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

% ========================= PART 2: Calculate ICs =========================================================================================================
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

end% End part1
