% COLLECT_rhythmGLMfit.m
%
% Collector script for output of ALL_rhythmGLMfit

clear;
cfg.inputDir = 'C:\temp\GLMfit'; % where the files to load are
%cfg.inputDir = 'C:\temp\GLMfit\HC'; % where the files to load are
cfg.input_prefix = 'D0_'; 
%cfg.models = {'allphi'};
cfg.nMaxCells = 1000;
cfg.nTimeBins = 100;
cfg.nSpaceBins = 100;

%%
pushdir(cfg.inputDir);

fd = FindFiles(cat(2,cfg.input_prefix,'*.mat'), 'CheckSubdirs', 0);
nSessions = length(fd);

popdir;

%% initialize variables
% for all non-baseline models:
nBins = 100;

ALL_celltype = nan(1, cfg.nMaxCells); % redundant across models but will be easier to deal with later
ALL_ttr = nan(nBins, cfg.nMaxCells);
ALL_duplicate = nan(1, cfg.nMaxCells);
ALL_sameTT = nan(1, cfg.nMaxCells);
ALL_insuffSpk = nan(1, cfg.nMaxCells);
ALL_accepted = nan(1, cfg.nMaxCells);

% counters
cellCount = 1; % MATLAB indexing starts at 1


%% loop across sessions to build variables

for iS = 1:nSessions

    load(fd{iS});
    
    nCells = length(sd.cellLabel);
    start_idx = cellCount; end_idx = cellCount + nCells - 1;

    ALL_cellType(start_idx:end_idx) = sd.cellType;
    ALL_ttr(:, start_idx:end_idx) = sd.descr.ttr';
    ALL_duplicate(start_idx:end_idx) = sd.descr.duplicate;
    ALL_sameTT(start_idx:end_idx) = sd.descr.sameTT;
    ALL_insuffSpk(start_idx:end_idx) = sd.descr.insuffSpk;
    ALL_accepted(start_idx:end_idx) = sd.descr.accepted;
    
    cellCount = cellCount + nCells;

end
cellCount = cellCount - 1;

%% plot
%cellTypes = {'MSNs', 'FSIs', 'all', 'both'};
cellTypes = {'both'};

for iC = 1:length(cellTypes)

    switch cellTypes{iC}
       
        case 'MSNs'
            cell_keep = find(ALL_cellType(1,:) == 1);
        case 'FSIs'
            cell_keep = find(ALL_cellType(1,:) == 2);
        case 'both'
            cell_keep = find(ALL_cellType(1,:) == 1 | ALL_cellType(1,:) == 2);
        otherwise
            cell_keep = find(~isnan(ALL_cellType(1,:)));
    end
    
    this_keep = intersect(cell_keep, find(ALL_accepted));
    %cellCount = length(cell_keep);
    
    % TCs - accept
    this_tc = ALL_ttr(:, this_keep); 
    
    subplot(121);
    imagesc(this_tc');
    
    % TCs - reject
    this_keep = intersect(cell_keep, find(ALL_insuffSpk));
    this_tc = ALL_ttr(:, this_keep); 
    
    subplot(122);
    imagesc(this_tc');
    
    
end % of cell types