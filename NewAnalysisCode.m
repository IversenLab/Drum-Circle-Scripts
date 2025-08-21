%analyze rhythm network experiments conducted at McMaster March 2025

close all

%location of data
fpath = '/Users/jonathankirsh/Documents/DrumCircle/Aug15Pilot/07_ring_trios/';
fextension = '.experiment';

%% Experiment Mac 3 (08/12/25)
if 1
    snames = {'DCP21','DCP22','DCP23','DCP24','DCP25','DCP26'};
    exppath = '/Users/jonathankirsh/Documents/DrumCircle/Aug15Pilot/FIGURES/';
    fnames = {'07a_ring_nodelay_trios','07b_ring_135ft_trios','07c_ring_270ft_trios'};
    nTappers = 6;
    fnames = fnames(1:3)
end

%% Experiment Mac 2 (03/28/25)
%if 0
%    snames = {'A1','A2','A3','A4','A5','A6'};
%    exppath = 'Experiment Mac2 032825';
%    fnames = {'03a2_alltoall,free_DS2', '03a3_alltoall,free_DS3','04_alltoall,unsynchronized', ...
%        '07a_ring,instant_DS1', '07a_ring,instant_DS2', '08_whisperRing_DS1', '09_subnetcollideA', '09_subnetcollideB', '10_alltoall,freeJam_DS3'};
%    nTappers = 6;
%end

longgap = 2500;

%% Set up for R

outdir = fullfile(exppath, 'exports_R');
if ~exist(outdir, 'dir'), mkdir(outdir); end

combined_csv = fullfile(outdir, 'all_to_all_trials_long.csv');

%% iterate over the experiment trials
ftodo = 1:length(fnames);
for fi = ftodo
    
    fbase = [fnames{fi} fextension];
    loadfname = fullfile(fpath, fbase);
    S = loadXMLPlist(loadfname);

    %% process experiment events
    %it's a long string, delimited by line breaks
    % lines have the format:
    %time (ns)  channel note velocity
    %389053395	16	64	100

    events = eval([ '[' S.recordedEvents ']' ]);     

    % there are repeated events stored in the recording.
    % Identify unique events based on columns 2 to 4, and keep only first
    % Initialize output with the first row
    filtered = events(1, :);
    % Loop through from the second row
    for i = 2:size(events, 1)
        prev = events(i-1, 2:4);  % previous event (no timestamp)
        curr = events(i, 2:4);    % current event (no timestamp)

        if ~isequal(curr, prev)
            filtered = [filtered; events(i, :)]; %#ok<AGROW>
        end
    end
    events = filtered;

    %% initialize the figure
    figure
    jisubplot(1,1,1,'landscape')
    hold on
    timestamp

    %% make timeseries of ITI for tappers
    % various possible color gradients
    % colors = 'rgbmyk';
    % colors = [
    %     0.00, 0.45, 0.74;  % blue
    %     0.00, 0.75, 0.75;  % cyan
    %     0.00, 0.60, 0.00;  % green
    %     0.80, 0.80, 0.00;  % yellow
    %     0.85, 0.33, 0.10;  % orange
    %     0.93, 0.00, 0.00;  % red
    %     ];
    colors_all = lines(6);  % Using the 'lines' colormap for 6 distinguishable colors
   
    lh = []; %line handles
    allevents = [];
    clear xtap

    %chans = 1:nTappers;
    %numChannels = 6; %trying this out

    % Define groups of tappers
    groups = {
        [1, 2, 3],  ... Group 1: define channels
        [4, 5, 6]   % Group 2: define channels
    };
    
    % Set the group index
    iG = 1;  % Change this to switch between groups
    nTappers = length(groups{iG});
    chans = groups{iG};
    numChannels = numel(chans); % Number of channels in the selected group

    colors_local = zeros(numChannels,3);
    for k = 1:numChannels
        colors_local(k,:) = colors_all(chans(k),:);
    end

    times = cell(numChannels, 1); % Initialize cell array for tap times
    ITI = cell(numChannels, 1); % Initialize cell array for inter-tap intervals
    velocities = cell(numChannels, 1); % Initialize cell array for velocities
    lh = gobjects(numChannels, 1); % Initialize line handles for plotting
    xtap = cell(numChannels,1);

    for k = 1:numChannels
        chan = chans(k);
        iChan = find(events(:,2)==chan-1);

        if ~isempty(iChan)
            % get tap times, compute ITI
            chanEvents = events(iChan,:);
            times{k} = chanEvents(:,1) / 1e6; %ns to ms
            ITI{k} = diff(times{k});
            ITI{k}(ITI{k}>longgap) = nan; %blank any long gaps

            % remove missing taps
            % Set the window size (e.g., 5 samples before and after)
            windowSize = 11;  % total window size, must be odd for symmetry
            % Compute the local median using a sliding window
            localMedian = movmedian(ITI{k}, windowSize, 'omitnan');
            % Mark values that are > 1.8 × local median as NaN
            threshold = 1.8;
            ITI_filt = ITI{k};  % copy original
            ITI_filt(ITI{k} > threshold * localMedian) = NaN;
            ITI{k} = ITI_filt;

            velocities{k} = chanEvents(:,4);
            times{k}(1)=[]; %omit first tap to use as time axis for ITI timeseries
            velocities{k}(1)=[];

            % plot ITI timeseries
            hold on
            lh(k) = plot(times{k}/1e3, ITI{k}, 'color', colors_local(k,:));
            % plot actual taps as a raster
            %plot_grid(times{chan}/1e3, [0 100],[],[],colors(chan));
            xtap{k} = times{k}/1e3;
            ytap = (k - (numChannels+1)/2)*20 * ones(size(xtap{k}));
            plot(xtap{k}, ytap, '.','MarkerEdgeColor',colors_local(k,:), 'MarkerFaceColor', colors_local(k,:), 'markersize',5)
            allevents = [allevents; xtap{k}];
            
            % optional: plot tap velocities
            % jisubplot(2,1,2)
            %hold on
            %plot(times{chan} / 1e3,velocities{chan},[colors(chan) '-'])

        end %if we have tap data
    end %loop over participants

    xlabel('time (s)')
    ylabel('Inter-tap Interval (ms)')

    %need to parse experiment details to get this kind of info
    %jisubplot(2,1,1)
    %gridy(1000,[.5 .5 1])
    %gridy(800,[1 .5 .5])
    %gridy(500,[.5 1 .5])

    %plot_grid_x([1 11 31]*1000)
    ylim([0 1000])
    title(protect_underscore([exppath ': ' fnames{fi}]))

    %jisubplot(2,1,2)
    %ylim([0 128])
    %title('velocity')
    
    % optional: event density plot
    if 0
        jisubplot(2,1,2)
        title('event density')
        allevents = sort(allevents);
        maxt=max(allevents);
        bin = 0.05;
        step = .01;
        alln = [];
        allt = [];
        for offset = (0:step:(bin-step))
            t = (0:bin:maxt) + offset;
            n=histc(allevents,t);
            alln = [alln; n];
            allt = [allt; t'];
        end
        [t, idx] = sort(allt);
        n = alln(idx);
        tplot = t + bin/2;
        plot(tplot,n,'color',[1 .5 .5])

        %variance measure of synchrony
        %   (problem, scales w/ rate, too), smple try: normalize by sum
        %       sum must be proportional to number of taps
        z=[];
        s=[];
        stride = 10;
        tidx = 1:stride:length(tplot)-stride;
        for i = tidx
            span = i:(i+stride-1);
            z = [z std(n(span))];
            s = [s sum(n(span))];
        end
        z=z.*z; %variance

        tidx = tidx+stride/2;
        hold on
        plot(tplot(tidx),z*(stride/100),'g')
        plot(tplot(tidx),z./s*(stride*2),'b','linewidth',1)
        plot(tplot(tidx),s/(stride/5),'k')
        %ylim([0 8])
        xl=xlim;
        xlim([0 xl(2)])
        delete(gca)
    end
    
    %jisubplot(2,1,1)
    %plot(tplot,n*10,'color',[1 .5 .5])
    xl=xlim;
    xlim([0 xl(2)])
    yl=ylim;
    ylim([-60 yl(2)])

    %% global measure of synchrony
    % instead of pairwise relative phase, compute a
    % continuous phase measure for each event series and use a kuromoto order parameter on
    % these:

    % events: 1xN cell array of event time vectors
    % t_grid: vector of times at which to evaluate synchrony
    % relphase: your existing function
    tt = events(:,1) / 1e9; %ns to s
    t0=tt(1); %adjust first event time to 0
    tt = tt-t0;
    tend = tt(end);

    t_grid = 0:.001:tend; %calculate phase at 1ms intervals

    
    phi = NaN(numChannels, numel(t_grid));  % preallocate phases for each participant at each time point
    
    for k = 1:numChannels
        if ~isempty(times{k})
        rp = relphase([], times{k}/1e3, t_grid);  % relative phase [0, 1)
        phi(k, :) = 2 * pi * rp;                  % convert to radians
        end
    end

    % Compute Kuramoto order parameter over time
    R = abs(nanmean(exp(1i * phi), 1));  % 1 x length(t_grid)

    Rs = jsmooth(2001, R); %3s window (smooth over multiple taps)

    % Savitsky-Golay filter
    dt = mean(diff(t_grid));
    windowSize = round(2 / dt);
    if mod(windowSize,2)==0, windowSize=windowSize+1; end  % must be odd
    polyOrder = 2;  % quadratic fit
    R_sgolay = sgolayfilt(R, polyOrder, windowSize);

    plot(t_grid, R*100+100,'r','LineWidth',2) %raw
    plot(t_grid, Rs*100+100,'b:') %smoothed (for comparison)
    plot(t_grid, R_sgolay*100+100,'k','LineWidth',2) %filtered

    gridy([100 200])

        %legend
    legend(lh, snames(chans),'location','northeast')

    % ========= LONG-FORMAT EXPORT TO CSV (per trial + combined) =========

    % We'll create one row per ITI sample per drummer, including metadata and synchrony.
    % Interpolate synchrony (R) onto each drummer's ITI time stamps for convenience
    % (tap times are in seconds: times{k}/1e3 gave seconds earlier; here we already used /1e3 when plotting)
    trial_id   = string(fnames{fi});
    group_id   = iG;        % which subgroup (your selector)
    trial_csv  = fullfile(outdir, [fnames{fi} '_long.csv']);
    T_long = table;  % accumulator for this trial
 
    for k = 1:numChannels
    % Guard against empty channels (e.g., dropouts)
        if isempty(times{k}) || isempty(ITI{k})
            continue
        end
 
        tap_time_s = times{k} / 1e3;   % you plotted with /1e3 earlier; keep consistent units here
        iti_ms     = ITI{k};
        vel        = velocities{k};
        nRows      = numel(iti_ms);
 
        % Synchrony at those times

        R_at_t     = interp1(t_grid, R,        tap_time_s, 'linear', 'extrap');
        Rs_at_t    = interp1(t_grid, R_sgolay, tap_time_s, 'linear', 'extrap');
 
        drummer_id   = chans(k);                 % original channel number (1..6)
        drummer_name = string(snames{drummer_id});
 
        % Assemble trial-long rows

        T_k = table( ...
            repmat(trial_id, nRows, 1), ...
            repmat(group_id, nRows, 1), ...
            repmat(fi, nRows, 1), ...
            repmat(drummer_id, nRows, 1), ...
            repmat(drummer_name, nRows, 1), ...
            tap_time_s(:), ...
            iti_ms(:), ...
            vel(:), ...
            R_at_t(:), ...
            Rs_at_t(:), ...
            'VariableNames', {'Trial','Group','TrialIndex','Drummer','DrummerName', ...
                              'TapTime_s','ITI_ms','Velocity','R_raw','R_sgolay'});
 
        T_long = [T_long; T_k]; %#ok<AGROW>
    end
 
    % Write per-trial CSV (new file each time)
    writetable(T_long, trial_csv, 'FileType','text','Delimiter',',','QuoteStrings',false,'WriteRowNames',false);
 
    % Append to combined CSV (create if it doesn't exist)
    if ~isfile(combined_csv)
        writetable(T_long, combined_csv, 'FileType','text','Delimiter',',','QuoteStrings',false,'WriteRowNames',false);
    else
            % R2020b+ supports WriteMode
        writetable(T_long, combined_csv, 'FileType','text','Delimiter',',','QuoteStrings',false,'WriteRowNames',false, 'WriteMode','append');
    end

% ====================================================================

    % save figure
        expmulti(gcf, 'pdf', fullfile(exppath, fnames{fi}), 600);
        close

end %loop over experiments

%%analyze driving stimulus timing
if 0
    chan = 15;
    iChan = find(events(:,2)==chan);
    chanEvents = events(iChan,:);
    times = chanEvents(:,1) / 1e6; %to milliseconds

    %now pull out scheduled times
    scheduledTimes = eval([ '[' S.partTiming{3}.subEventTimes ']' ] ) / 1e6;

    %compare absolute error--in 1e-5 range (that's ms, so on order of 10 ns)
    st = [scheduledTimes times];
    error = diff(st,1,2)
    mean(error)

    %cumulative IOI error--do they drift apart? (NO: sum of errors is insignificant: 1e-13)
    d = diff(st);
    error = diff(d,1,2)
    sum(error)
end


