cd('/Users/manishm/Work/vanDerMeerLab/ADRLabData/R117/R117-2007-06-02');
cfg.output_dir = '/Users/manishm/Work/vanDerMeerLab/GLM/temp';
cfg.writeOutput = 1;
cfg.plotOutput = 0;
cfg.output_prefix = 'S2_Task_'; % prefix filenames with this (identify runs)
please.rats = {'R117'};
[cfg.fd,cfg.fd_extra] = getDataPath(please);
cfg.iS = 1;
cfg.nMinSpikes = 100;
% cfg.nShuf = 10;
sd = rhythmGLM_TC_Task(cfg);