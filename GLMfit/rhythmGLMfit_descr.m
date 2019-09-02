function sd = rhythmGLMfit_descr(cfg_in)
% Descriptive statistics for GLM fit procedure
%
% doesn't actually fit GLMs, just gets tuning curves, cell counts etc.
%
% This top-level function fits a number of GLMs to single session spike train data.
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
cfg_master.ttr_bins = [-5:0.1:5]; % time bin edges for time-to-reward tuning curves
cfg_master.linposBins = 101; % number of position bins (is autoscaled for each session)
cfg_master.nMinSpikes = 100; % minimum number of spikes needed to include cell
cfg_master.ccMethod = 'MvdM'; % cell type classification method
cfg_master.maxPrevCorr = 0.99; % if wv correlation with previous day is bigger than this, cell is possible duplicate
cfg_master.maxPeakn = 0.2; % if peak wv difference (normalized) with previous day is smaller than this, cell is possible duplicate
cfg_master.iS = []; % current session number out of fd list, get this from input cfg
cfg_master.fc = []; % full list of session fd's, get this from input cfg
cfg_master.fc_extra = []; % get this from input cfg
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

% LFP - vStr
if isfield(ExpKeys,'goodGamma_vStr')
    cfg = []; cfg.fc = ExpKeys.goodGamma_vStr(1);
elseif isfield(ExpKeys,'goodGamma')
     cfg = []; cfg.fc = ExpKeys.goodGamma(1);
else
    error('Don''t know what LFP to load.');
end
lfp_tt = regexp(cfg.fc, 'CSC\d+', 'match');
lfp_tt = str2double(lfp_tt{1}{1}(4:end)); % need this to skip cells from same tt (could make into function)

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
TVECe = ExpKeys.TimeOnTrack:cfg_master.dt:ExpKeys.TimeOffTrack; % edges
sd.TVECc = TVECe(1:end-1)+cfg_master.dt/2; % centers

% time-to-reward variable is used to restrict samples
t_to_reward = ComputeTTR(reward_t,sd.TVECc);
cfg_ttr = []; cfg_ttr.interp = 'nearest'; cfg_ttr.bins = cfg_master.ttr_bins;
spk_binned = logical(zeros(size(sd.TVECc))); % dummy spike train just to get binned TTR
[ttr, lambda_ttr, ttr_binned] = MakeTC_1D(cfg_ttr, sd.TVECc, t_to_reward, sd.TVECc, spk_binned);

MASTER_keep = ~isnan(ttr); % NEED THIS TO RESTRICT BINNED SPIKE TRAIN LATER!

t_to_reward = t_to_reward(MASTER_keep);
sd.TVECc = sd.TVECc(MASTER_keep);

% some counters
sd.descr.duplicate = nan(nCells, 1);
sd.descr.sameTT = nan(nCells, 1);
sd.descr.insuffSpk = nan(nCells, 1);
sd.descr.accepted = ones(nCells, 1); % accept by default, but set to 0 if fail

%% loop over all cells 
for iC = nCells:-1:1

    fprintf('Cell %d/%d...\n',iC,nCells);
 
    % skip if not turned & correlated with prev session
    if sd.S.usr.duplicate(iC) & sd.S.usr.distanceTurned(iC) < 80
        sd.descr.duplicate(iC) = 1;
        sd.descr.accepted(iC) = 0;
    else
        sd.descr.duplicate(iC) = 0;
    end
    
    % skip if on same tt as LFP
    if sd.S.usr.tt_num(iC) == lfp_tt
        sd.descr.sameTT(iC) = 1;
        sd.descr.accepted(iC) = 0;
    else
        sd.descr.sameTT(iC) = 0;
    end
    
    % dependent variable: binned spike train
    spk_binned = histc(sd.S.t{iC}, TVECe); spk_binned = spk_binned(1:end - 1);    
    spk_binned = logical(spk_binned > 0); % binarize
    
    %%% now can restrict spike train
    spk_binned = spk_binned(MASTER_keep);
    
    %%% SKIP CELL IF NOT ENOUGH NEURONS
    if sum(spk_binned) <= cfg_master.nMinSpikes
        sd.descr.insuffSpk(iC) = 1;
        sd.descr.accepted(iC) = 0;
    else
        sd.descr.insuffSpk(iC) = 0;
    end
    
    %%% INCLUDE SOME INFO ABOUT THIS CELL
    sd.cellType(iC) = sd.S.usr.cell_type(iC);
    sd.cellLabel{iC} = sd.S.label(iC);
    sd.cellDepth(iC) = sd.S.usr.tetrodeDepths(iC);
    sd.cellID(iC) = s_out.cq.id(iC); sd.cellLr(iC) = s_out.cq.lr(iC); sd.cellAmpl(iC) = s_out.cq.ampl(iC);
          
    %%% PREDICTOR: time to reward
    [ttr, lambda_ttr, ttr_binned] = MakeTC_1D(cfg_ttr, sd.TVECc, t_to_reward, sd.TVECc, spk_binned);
    sd.descr.ttr(iC, :) = lambda_ttr ./ cfg_master.dt;

end % over cells

%% prepare and write output
if cfg_master.writeOutput

    sd.cfg = cfg_master;
    
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
function linpos_tsd = ComputeLinPos(pos)

f = FindFiles('*Coord*');
if ~isempty(f)
    load(f{1});
else
    LoadMetadata;
    Coord = metadata.Linpath.Coord8;
end

Coord_out = [];
Coord_out.coord = Coord;
Coord_out.nPoints = size(Coord,2);
Coord_out.pointDist = pdist(Coord(:,1:2)','euclidean');
Coord_out.standardized = 0;
Coord_out.units = 'px';

linpos_tsd = LinearizePos([],pos,Coord_out);

end

%%
function t_to_reward = ComputeTTR(reward_t,TVECc)
% find time to closest reward at each time point in TVECc

t_to_reward = TVECc-reward_t(1); % running total
for iR = 2:length(reward_t)
    
    this_dt = TVECc-reward_t(iR);
    swap_idx = find(abs(this_dt) < abs(t_to_reward)); % samples that are closer than running total
    t_to_reward(swap_idx) = this_dt(swap_idx);
end

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
