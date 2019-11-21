function sd = rhythmGLMfit_shufPRE(cfg_in)
% GLM for spike prediction with LFP features
%
% as rhythmGLMfit, but add exclusion of duplicate cells
%
% This top-level function fits a number of GLMs to single session spike train data.
%
% The overall idea is to test whether addition of LFP features, such as
% theta power and/or gamma phase, improve cross-validated model prediction.
%
% The analysis proceeds as follows:
% - load data
% - define models (in sd.m.modelspec), MUST include a baseline model
% - prepare session-wide variables (linearized position, LFP features,
%   speed etc) on common timebase ('TVECc')
% - for each cell:
%   * prepare regressors for this ell
%   - for each cross-validation run ("pleat"; a set of folds) and fold:
%     + fit models on training data
%     + test models
%     + store error relative to baseline model
% - for each model:
%   * plot error across cells
%   * plot error across time to reward (tuning curve)
%
% run this with data session to analyze as the working folder.
%
% required path: striatal-spike-rhythms github repo

%% params
cfg_master = []; % overall params
cfg_master.dt = 0.001;
cfg_master.f.delta = [2 5];
cfg_master.f.theta = [6.5 9.5];
cfg_master.f.beta = [14 25];
cfg_master.f.lowGamma = [40 65];
cfg_master.f.highGamma = [70 100];
cfg_master.nPleats = 1;
cfg_master.kFold = 2;
cfg_master.plotOutput = 1;
cfg_master.writeOutput = 0;
cfg_master.writeFullError = 0; % write full nSamples x nCells x nModels error for each model (NOTE: takes up to 1GB per session)
cfg_master.smooth = 501; % smoothing window (in samples) for error
cfg_master.output_prefix = 'S0_';
cfg_master.output_dir = 'C:\temp';
cfg_master.linposBins = 101; % number of position bins (is autoscaled for each session)
cfg_master.nMinSpikes = 100; % minimum number of spikes needed to include cell
cfg_master.ccMethod = 'MvdM'; % cell type classification method
cfg_master.maxPrevCorr = 0.99; % if wv correlation with previous day is bigger than this, cell is possible duplicate
cfg_master.maxPeakn = 0.2; % if peak wv difference (normalized) with previous day is smaller than this, cell is possible duplicate
cfg_master.iS = []; % current session number out of fd list, get this from input cfg
cfg_master.fd = []; % full list of session fd's, get this from input cfg
cfg_master.fd_extra = []; % get this from input cfg
cfg_master.nShuf = 1;

cfg_master = ProcessConfig(cfg_master,cfg_in);

%% loading
% load data
LoadExpKeys;

% spikes
sd.S = LoadSpikesTarget(cfg_master);
nSpikes = cellfun(@length, sd.S.t); keep = nSpikes >= cfg_master.nMinSpikes;
sd.S = SelectTS([], sd.S, keep);
nCells = length(sd.S.t); if nCells == 0, sd = []; return; end

%% Categorize cells and add tetrode depths
cfg_wv = []; cfg_wv.cMethod = cfg_master.ccMethod;
s_out = CategorizeStriatumWave(cfg_wv, sd.S);

s_out.unit = [s_out.other s_out.msn s_out.fsi];
s_out.ident = [zeros(1, length(s_out.other)) ones(1, length(s_out.msn)) repmat(2, 1, length(s_out.fsi))];

cfg_tt = []; cfg_tt.verbose = 1;
cfg_tt.this_rat_ID = cfg_master.fd_extra.ratID_num(cfg_master.iS);
cfg_tt.this_date = cfg_master.fd_extra.fd_date_num(cfg_master.iS);

for iC = 1:length(sd.S.t)
    sd.S.usr.cell_type(iC) = s_out.ident(find(s_out.unit == iC));
    sd.S.usr.tetrodeDepths(iC) = ExpKeys.TetrodeDepths(sd.S.usr.tt_num(iC));
    
    cfg_tt.ttno = sd.S.usr.tt_num(iC);
    [sd.S.usr.distanceTurned(iC), prev_fd] = DistanceTurned(cfg_tt, cfg_master.fd, cfg_master.fd_extra);
    cfg_tt.verbose = 0;
end

% correlate with previous session waveforms if available
if isempty(prev_fd) % no previous day available
    sd.S.usr.duplicate = zeros(size(sd.S.usr.tt_num));
else
    pushdir(prev_fd);
    S2 = LoadSpikes([]);
    nSpikes = cellfun(@length, S2.t); keep = nSpikes >= cfg_master.nMinSpikes;
    S2 = SelectTS([], S2, keep);
    
    s_out2 = CategorizeStriatumWave(cfg_wv, S2);
    s_out = CalcWVDistances([], s_out, s_out2); % add comparison with previous day's waveforms
    
    popdir;
    
    % for each cell in current session, decide if duplicate
    for iC = 1:length(sd.S.t)
        
        this_tt_no = sd.S.usr.tt_num(iC);
        prev_day_cells = find(S2.usr.tt_num == this_tt_no);
        
        if isempty(prev_day_cells) % no cells recorded fron this tt in previous session
            sd.S.usr.duplicate(iC) = 0;
        else % previous day cells found
            temp_corr = s_out.corr(iC, prev_day_cells);
            temp_peakn = s_out.peakdiffn(iC, prev_day_cells);
            
            if temp_corr > cfg_master.maxPrevCorr & abs(temp_peakn) < cfg_master.maxPeakn % wv correlation big, peak difference small
                sd.S.usr.duplicate(iC) = 1;
            else
                sd.S.usr.duplicate(iC) = 0;
            end
        end
    end
end % of previous day available checks

% LFP - vStr
if isfield(ExpKeys,'goodGamma_vStr')
    cfg = []; cfg.fc = ExpKeys.goodGamma_vStr(1);
elseif isfield(ExpKeys,'goodGamma')
     cfg = []; cfg.fc = ExpKeys.goodGamma(1);
else
    error('Don''t know what LFP to load.');
end
csc = LoadCSC(cfg);

%%finding the first large timeskip (greatert than 1 sec)
timeDiff = diff(csc.tvec);
timeSkipIdx = find(timeDiff > 1);
timeSkip = csc.tvec(timeSkipIdx(1)); %the first "big" jump

lfp_tt = regexp(cfg.fc, 'CSC\d+', 'match');
lfp_tt = str2double(lfp_tt{1}{1}(4:end)); % need this to skip cells from same tt (could make into function)
fprintf('LFP ttno is %d\n', lfp_tt);

cfg_phi = []; % LFP features
cfg_phi.dt = median(diff(csc.tvec));
cfg_phi.ord = 100;
cfg_phi.bins = -pi:pi/18:pi;
cfg_phi.interp = 'nearest';

% position
if isempty(ls('*.nvt')) % hmm
    fn = ls('*1.zip');
    system(cat(2, '7z x ', fn));
end
pos = LoadPos([]);

% reward deliveries
reward_t = getRewardTimes;
fprintf('%d trials found.\n', length(reward_t));

%% initialize variables
% overall plan:
%
% for each cell, p table contains regressors. these are NOT stored across cells (into session data) for memory reasons
% but maybe session-wide variables should be stored (time_to_reward,
% linpos, LFP features) for later analysis
%
% across cells, sd (session data) struct tracks:
% - TVECc: centers of time bins for data to be fitted [1 x nTimeBins]
% for each model:
% - .m.(modelName).err (nCells x nTimeBins)
% - .m.(modelName).modelspec

p = table;

% common timebase for this session
%TVECe = ExpKeys.TimeOnTrack:cfg_master.dt:ExpKeys.TimeOffTrack; % edges

% Changing Timebase to include only prebehabiour resting period
TVECe = 0:cfg_master.dt:min(timeSkip,(ExpKeys.TimeOnTrack-10)); %edges
%%Restrict timebase to start fro beginning of csc
startIdx = find(TVECe >= csc.tvec(1),1);
TVECe = TVECe(startIdx:end);

sd.TVECc = TVECe(1:end-1)+cfg_master.dt/2; % centers


nMaxVars = 15; % only used for initializing t-stat matrix
% baseline model MUST be defined first or things will break!
sd.m.baseline.modelspec = 'spk ~ 1 + cif';
sd.m.tphi.modelspec = 'spk ~ 1 + cif + theta_phase';
sd.m.bphi.modelspec = 'spk ~ 1 + cif + beta_phase';
sd.m.dphi.modelspec = 'spk ~ 1 + cif + delta_phase';
sd.m.lgphi.modelspec = 'spk ~ 1 + cif + lowGamma_phase';
sd.m.hgphi.modelspec = 'spk ~ 1 + cif + highGamma_phase';
sd.m.allphi.modelspec = 'spk ~ 1 + cif + delta_phase + beta_phase + theta_phase + lowGamma_phase + highGamma_phase';

% init error vars -- now done for each cell
%mn = fieldnames(sd.m);
%for iM = 1:length(mn)
%   sd.m.(mn{iM}).err = zeros(nCells, length(sd.TVECc)); % needs to be zeros because error output will be added to this
%   sd.m.(mn{iM}).tstat = nan(nCells, nMaxVars);
%end

% define training and testing sets
for iPleat = cfg_master.nPleats:-1:1
    C{iPleat} = cvpartition(length(sd.TVECc), 'KFold', cfg_master.kFold);
end

% LFP features
disp('Computing session-wide LFP features...');
fb_names = fieldnames(cfg_master.f);
cfg_phi.debug = 0;
for iF = 1:length(fb_names)
   
    cfg_phi.fpass = cfg_master.f.(fb_names{iF});
    fprintf('%s: [%d %d]\n', fb_names{iF}, cfg_phi.fpass(1), cfg_phi.fpass(2));
    
    [FF.(fb_names{iF}).phase, FF.(fb_names{iF}).env] = ComputePhase(cfg_phi, csc);
    
end


%%% time / nTrials predictor
p.time = sd.TVECc';

%% loop over all cells 
cc = 1;
sd.duplicateCells = [];
sd.fewSpikes = [];
sd.sameTT = [];
for iC = nCells:-1:1

    fprintf('Cell %d/%d...\n',iC,nCells);
    % skip if not turned & correlated with prev session
    if sd.S.usr.duplicate(iC) & sd.S.usr.distanceTurned(iC) < 80
        fprintf('Cell skipped - likely duplicate.\n');
        sd.duplicateCells = [sd.duplicateCells; sd.S.label(iC), sd.S.usr.cell_type(iC)];
        continue;
    end
    
    % skip if on same tt as LFP
    if sd.S.usr.tt_num(iC) == lfp_tt
        fprintf('Cell skipped - same TT as LFP.\n');
        sd.sameTT = [sd.sameTT; sd.S.label(iC), sd.S.usr.cell_type(iC)];
        continue;
    end
    
    
    % dependent variable: binned spike train
    spk_binned = histc(sd.S.t{iC}, TVECe); spk_binned = spk_binned(1:end - 1);
    
    spk_binned = logical(spk_binned > 0); % binarize
    %k = gausskernel(50,2); % smoothing works too, should make into option
    %spk_binned = conv2(spk_binned,k,'same');
    
    %%% PREDICTOR: conditional intensity function
    % note we need to do this before restricting the spike train to avoid weird breaks
    cif = conv(spk_binned,[0 0 0 1 1],'same'); % 2 ms refractory period
    cif = (cif > 0);
    
    %%% PREDICTOR: better CIF based on acorr
    % why does this work so well?
    % what happens for a symmetric acorr?
    % what happens when x-val is done on "blocked" rather than random folds?
    % how do results depend on this thing being present or not?
    cfg_acf = []; cfg_acf.maxlag = 500; cfg_acf.binsize = cfg_master.dt;
    cfg_acf.sided = 'onezero';
    [acf, acf_tvec] = ComputeACF(cfg_acf, spk_binned);
    
    cif_full = conv(spk_binned, acf ./ sum(acf), 'same');
    cif_full = cif_full-nanmean(cif_full);
    cif_full(cif) = 0; % absolute refractory period
    
    p.cif = cif_full;
    
    %%% SKIP CELL IF NOT ENOUGH SPIKES
    if sum(spk_binned) <= cfg_master.nMinSpikes
       fprintf('\n\n*** CELL SKIPPED - INSUFFICIENT SPIKES ***\n\n');
       sd.fewSpikes = [sd.fewSpikes; sd.S.label(iC), sd.S.usr.cell_type(iC)];
       continue;
    end
    
    %%% INCLUDE SOME INFO ABOUT THIS CELL
    sd.cellType(cc) = sd.S.usr.cell_type(iC);
    sd.cellLabel{cc} = sd.S.label(iC);
    sd.cellDepth(cc) = sd.S.usr.tetrodeDepths(iC);
    sd.cellID(cc) = s_out.cq.id(iC); sd.cellLr(cc) = s_out.cq.lr(iC); sd.cellAmpl(cc) = s_out.cq.ampl(iC);


    %%% PREDICTOR: LFP phase and envelope for all defined bands %%%
    what = {'phase','env'};
    for iF = 1:length(fb_names) % loop over frequency bands
        for iW = 1:length(what)
            this_name = cat(2,fb_names{iF},'_',what{iW});
            switch iW
                case 1 % use phase bins
                    cfg_temp = cfg_phi;
                case 2 % envelope bins depend on data
                    cfg_temp = cfg_phi; cfg_temp.bins = linspace(min(FF.(fb_names{iF}).(what{iW})), 0.5*max(FF.(fb_names{iF}).(what{iW})), 100);
            end
            p.(this_name) = MakeTC_1D(cfg_temp, csc.tvec, FF.(fb_names{iF}).(what{iW}), sd.TVECc, spk_binned)';
        end
    end
    
    % init model vars
    mn = fieldnames(sd.m);
    for iM = 1:length(mn)
       sd.m.(mn{iM}).err(cc, :) = zeros(1, length(sd.TVECc)); % needs to be zeros because error output will be added to this
       sd.m.(mn{iM}).shufErr(cc, 1:cfg_master.nShuf, 1) = 0; % can't save all shuffles, so need to just save the average value
       sd.m.(mn{iM}).tstat(cc, :) = nan(1, nMaxVars);
    end

    %% x-val loop
    p.spk = spk_binned;
    clear cif cif_full;
    
    for iPleat = 1:cfg_master.nPleats
        
        for iFold = 1:C{iPleat}.NumTestSets
            
            fprintf('Pleat %d fold %d...\n', iPleat, iFold);
            
            % get idxs for training and testing set
            tr_idx = C{iPleat}.training(iFold); te_idx = C{iPleat}.test(iFold);
            
            for iModel = 1:length(mn)
                
                % train initial model
                this_m = fitglm(p(tr_idx,:), sd.m.(mn{iModel}).modelspec, 'Distribution', 'binomial');
                sd.m.(mn{iModel}).tstat(cc, 1:length(this_m.Coefficients.tStat)) = this_m.Coefficients.tStat; % should initialize this, but that's a pain
                sd.m.(mn{iModel}).varnames = this_m.PredictorNames;
                
                % refine model by throwing out ineffective predictors
                %toss = this_m.Coefficients.Row(abs(this_m.Coefficients.tStat) < 2);
                %for iToss = 1:length(toss)
                %    if ~strcmp('(Intercept)',toss{iToss}) % can't remove intercept
                %        this_m.removeTerms(toss{iToss});
                %   end
                %end
                
                % test it and add resulting error to running total
                this_err = this_m.predict(p(te_idx,:));
                this_err = (this_err - spk_binned(te_idx)).^2;
                sd.m.(mn{iModel}).err(cc, te_idx) = sd.m.(mn{iModel}).err(cc, te_idx) + (this_err ./ cfg_master.nPleats)';
                
            end % across models
            
        end % over folds
        
    end % over pleats
    
    % run some shuffles on final model
    %iModel = length(mn);
    for iShuf = cfg_master.nShuf:-1:1
        
        fprintf('Shuffle %d...\n', iShuf);
        
        % shift LFP data relative to spike data & fit model again
        p_shuf = shiftLFPdata2(p);
           
        for iModel = 1:length(mn)
            %for each model, the placeholder variable needs to be set to
            %zero at the beginning of the loop
            this_shufErr = zeros(1, length(sd.TVECc)); 
            for iPleat = 1:cfg_master.nPleats
                
                for iFold = 1:C{iPleat}.NumTestSets
                    
                    tr_idx = C{iPleat}.training(iFold); te_idx = C{iPleat}.test(iFold);
                    
                    this_m = fitglm(p_shuf(tr_idx,:), sd.m.(mn{iModel}).modelspec, 'Distribution', 'binomial');
                    
                    % refine model by throwing out ineffective predictors
                    %toss = this_m.Coefficients.Row(abs(this_m.Coefficients.tStat) < 2);
                    %for iToss = 1:length(toss)
                    %    if ~strcmp('(Intercept)',toss{iToss}) % can't remove intercept
                    %        this_m.removeTerms(toss{iToss});
                    %    end
                    %end
                    
                    % test it and add resulting error to running total
                    this_err = this_m.predict(p_shuf(te_idx,:));
                    this_err = (this_err - spk_binned(te_idx)).^2;
                    
                    
                    
                    
                    % NOTE THIS ONLY WORKS WHEN DOING 1 SHUFFLE! OTHERWISE NEED TO
                    % INDEX BY SHUFFLE
%                     sd.m.(mn{iModel}).shufErr(cc, te_idx) = sd.m.(mn{iModel}).shufErr(cc, te_idx) + ((this_err ./ cfg_master.nShuf) ./ cfg_master.nPleats)';
                    this_shufErr(te_idx) = this_shufErr(te_idx) + (this_err ./ cfg_master.nPleats)';
                end % of shuffle folds
                
            end % of pleats
        this_shuf_Err = sd.m.baseline.err(cc,:) - this_shufErr;
        sd.m.(mn{iModel}).shufErr(cc, iShuf, 1) = nanmean(sd.m.baseline.err(cc,:) - this_shufErr);
        end % of models
        
    end % of shuffles
    
    % average across folds
    %shuf_err = nansum(shuf_err, 2); % first sum across folds
        
    % average across shuffles if needed
    %if cfg_master.nShuf > 1
    %    sd.m.(mn{iModel}).shufErrMean = nanmean(shuf_err);
    %    sd.m.(mn{iModel}).shufErrSD = nanstd(shuf_err);
    %else
    %    sd.m.(mn{iModel}).shufErrMean = shuf_err;
    %    sd.m.(mn{iModel}).shufErrSD = [];
    %end
    
    cc = cc + 1;
end % over cells
cc = cc - 1;

if cc == 0
    return;
end

%% analyze models
this_f = ones(cfg_master.smooth, 1) ./ cfg_master.smooth;

if cfg_master.plotOutput
    
    
    for iM = 1:length(mn)
        
        if strcmp(mn{iM},'baseline')
            continue;
        end
        
        figure;
        
        this_err = sd.m.baseline.err - sd.m.(mn{iM}).err;
        
        % mean error over cells
        celldiffmean = nanmean(this_err,2);
        
        subplot(311);
        plot(celldiffmean,'LineWidth',1); hold on;
        keep = find(celldiffmean > 0); plot(keep,celldiffmean(keep),'.g','MarkerSize',20);
        keep = find(celldiffmean < 0); plot(keep,celldiffmean(keep),'.r','MarkerSize',20);
        set(gca,'LineWidth',1,'FontSize',18); box off;
        ylabel('Prediction improvement'); xlabel('cell #');
        title(mn{iM});
        
        % t-stat across neurons -- HOW TO GET THE LABELS FROM THE MODELSPEC?
        subplot(312);
        imagesc(sd.m.(mn{iM}).tstat); caxis([0 10]);
        
        subplot(313);
        imagesc(corrcoef(sd.m.(mn{iM}).tstat));
        
        % could show each cell's model improvement across time
        
    end
end

%% prepare and write output
if cfg_master.writeOutput

    sd.cfg = cfg_master;
    
    % compute difference TCs for each model and each cell
    for iM = 1:length(mn)
    
        % get error diff for this model
        this_err = sd.m.baseline.err - sd.m.(mn{iM}).err;
        this_errShuf = sd.m.baseline.err - sd.m.(mn{iM}).shufErr;
        
        for iC = cc:-1:1
            
            % smooth
            this_cell_err = this_err(iC,:); this_cell_errShuf = this_errShuf(iC,:);
            this_cell_err = conv(this_cell_err, this_f, 'same'); this_cell_errShuf = conv(this_cell_errShuf, this_f, 'same');
            
        end % loop over cells
    
    end % loop over models
    
    % if requested, save full error variables
    if cfg_master.writeFullError
        fprintf('\n*** WARNING: saving full error variables in output! (may take up GB''s of space!\n\n');
    else % set error fields to average across samples
        sd.m.baseline.err = nanmean(sd.m.baseline.err, 2);
        for iM = 1:length(mn)
            sd.m.(mn{iM}).err = nanmean(sd.m.(mn{iM}).err, 2);
            sd.m.(mn{iM}).shufErr = nanmean(sd.m.(mn{iM}).shufErr, 2);
        end
    end
    
    % write
    [~, fp, ~] = fileparts(pwd);
    
    pushdir(cfg_master.output_dir);
    fn_out = cat(2,cfg_master.output_prefix, fp, '_sd.mat');
    save(fn_out,'sd'); % should add option to save in specified output dir
    popdir;
    
end
  
end % of top-level function

%% TODO: correlate GLM results with things like spike spectrum, STA spectrum (across cells)

%% TODO: correlate tuning curves with GLM results
% why are the TTR GLM fits so "streaky"?
% is there some correlation between TTR, space, movement onset TCs

%%
function out = zmat(in,cfg)
% z-score each row of input matrix
out = nan(size(in));
for iRow = size(in, 1):-1:1
    this_row = in(iRow,:);
    switch(cfg)
        case 'z'
            out(iRow,:) = (this_row - nanmean(this_row)) ./ nanstd(this_row);
        case 'max'
            out(iRow,:) = this_row ./ max(this_row);
    end
end
end
% could speed up by repmatting

%%
function S = KeepCells(S)
thr = 500;

nCells = length(S.t);
for iC = nCells:-1:1
   l(iC) = length(S.t{iC}); 
end
keep = l >= thr;

S.label = S.label(keep);
S.t = S.t(keep);
end

%%
function [phase, env] = ComputePhase(cfg, csc)
% compute csc.data phase and envelope in frequency band cfg.fpass

Fs = (1 ./ cfg.dt);

d = fdesign.bandpass('N,Fp1,Fp2,Ap', 10, cfg.fpass(1), cfg.fpass(2), .5, Fs);
Hd = design(d, 'cheby1');

if cfg.debug
   figure;
   fvtool(Hd);
end

data = filtfilt(Hd.sosMatrix, Hd.ScaleValues, csc.data);
h = hilbert(data);
phase = angle(h);
env = abs(h);

end

%%
function S = LoadSpikesTarget(cfg_in)

if ~isfield(cfg_in, 'Target') % no target specified, load them all
    S = LoadSpikes([]);
    return;
end
LoadExpKeys;

% target specified, need to do some work
% first see if this session has more than one target

nTargets = length(ExpKeys.Target);
if ~iscell(ExpKeys.Target) | (iscell(ExpKeys.Target) && length(ExpKeys.Target) == 1) % one target
    
    target_idx = strmatch(cfg_in.Target, ExpKeys.Target);
    if isempty(target_idx)
        S = ts;
    else
        S = LoadSpikes([]);
    end
    
else % multiple targets, assume TetrodeTargets exists
    
    please = []; please.getTTnumbers = 1;
    S = LoadSpikes(please);
    
    target_idx = strmatch(cfg_in.Target, ExpKeys.Target);
    tt_num_keep = find(ExpKeys.TetrodeTargets == target_idx);
    
    keep = ismember(S.usr.tt_num, tt_num_keep);
    S = SelectTS([], S, keep);
    
end

end

%%
function lfp_shifted_data = shiftLFPdata(input_data)
% circularly shift phase rows of data table by the same amount for each
% trial

trial_len = 10000; % 10 seconds @ 1 ms bin size
this_shift = round(trial_len * rand(1));
nTrials = size(input_data, 1) ./ trial_len;

rows_to_shift = find(cellfun(@isempty, (strfind(fieldnames(input_data),'_phase'))) == 0);

lfp_shifted_data = input_data.Variables;

for iT = 1:nTrials
   
    start_idx = (iT-1)*trial_len + 1;
    end_idx = iT*trial_len;
    
    lfp_shifted_data(start_idx:end_idx,rows_to_shift) = circshift(lfp_shifted_data(start_idx:end_idx,rows_to_shift), this_shift, 1);
    
end

lfp_shifted_data = array2table(lfp_shifted_data, 'VariableNames', input_data.Properties.VariableNames);

end

function lfp_shifted_data = shiftLFPdata2(input_data)
% circularly shift phase rows of data table by the same amount (no concept
% of trial in pre-trial data)

this_shift = round(size(input_data, 1) * rand(1));

rows_to_shift = find(cellfun(@isempty, (strfind(fieldnames(input_data),'_phase'))) == 0);

lfp_shifted_data = input_data.Variables;

lfp_shifted_data(:,rows_to_shift) = circshift(lfp_shifted_data(:,rows_to_shift), this_shift, 1);

lfp_shifted_data = array2table(lfp_shifted_data, 'VariableNames', input_data.Properties.VariableNames);

end
