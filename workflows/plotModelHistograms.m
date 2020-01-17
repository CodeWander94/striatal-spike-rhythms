%%Plot the prediction improvement of models across days%%
cd('/Users/manishm/Work/vanDerMeerLab/GLM/PreTrial')
rats = {'R117','R119','R131','R132'};
%rats = {'R132'};
models = {'dphi', 'tphi', 'bphi', 'lgphi', 'hgphi', 'allphi'};
mod_pred = {};
for idx = 1:length(models)
    mod_pred.(models{idx}) = [];
end

for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat('S2Pre_',curRat,'*sd.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name);
%        N = find(sd.cellType == 0);
        N = length(sd.cellType);
        for iC = 1:length(N)
            for iM = 1:length(models)
                mod_pred.(models{iM}) = [mod_pred.(models{iM});str2num(curRat(2:end)),sd.cellType(iC),nanmean(sd.m.baseline.err(N(iC),:) - sd.m.(models{iM}).err(N(iC),:))];
            end
        end
    end
end
%% For dphi
c0 = find(mod_pred.dphi(:,2) == 0);
c1 = find(mod_pred.dphi(:,2) == 1);
c2 = find(mod_pred.dphi(:,2) == 2);
figure
graphTitle = 'dphi';
g = subplot(5,1,1);
t = title(graphTitle,'FontSize',40);
g.Visible = 'off';
t.Visible = 'on';
subplot(5,1,2);
hist(mod_pred.dphi(:,3), -0.000005:0.0000005:0.0001)
title('All cells', 'fontsize',20)
subplot(5,1,3);
hist(mod_pred.dphi(c0,3), -0.000005:0.0000005:0.0001)
title('Type0', 'fontsize',20)
subplot(5,1,4);
hist(mod_pred.dphi(c1,3), -0.000005:0.0000005:0.0001)
title('Type1', 'fontsize',20)
subplot(5,1,5);
hist(mod_pred.dphi(c2,3), -0.000005:0.0000005:0.0001)
title('Type2', 'fontsize',20)

%% For tphi
c0 = find(mod_pred.tphi(:,2) == 0);
c1 = find(mod_pred.tphi(:,2) == 1);
c2 = find(mod_pred.tphi(:,2) == 2);
figure
graphTitle = 'tphi';
g = subplot(5,1,1);
t = title(graphTitle,'FontSize',40);
g.Visible = 'off';
t.Visible = 'on';
subplot(5,1,2);
hist(mod_pred.tphi(:,3), -0.000005:0.0000005:0.0001)
title('All cells', 'fontsize',20)
subplot(5,1,3);
hist(mod_pred.tphi(c0,3), -0.000005:0.0000005:0.0001)
title('Type0', 'fontsize',20)
subplot(5,1,4);
hist(mod_pred.tphi(c1,3), -0.000005:0.0000005:0.0001)
title('Type1', 'fontsize',20)
subplot(5,1,5);
hist(mod_pred.tphi(c2,3), -0.000005:0.0000005:0.0001)
title('Type2', 'fontsize',20)

%% For bphi
c0 = find(mod_pred.bphi(:,2) == 0);
c1 = find(mod_pred.bphi(:,2) == 1);
c2 = find(mod_pred.bphi(:,2) == 2);
figure
graphTitle = 'bphi';
g = subplot(5,1,1);
t = title(graphTitle,'FontSize',40);
g.Visible = 'off';
t.Visible = 'on';
subplot(5,1,2);
hist(mod_pred.bphi(:,3), -0.000005:0.0000005:0.0001)
title('All cells', 'fontsize',20)
subplot(5,1,3);
hist(mod_pred.bphi(c0,3), -0.000005:0.0000005:0.0001)
title('Type0', 'fontsize',20)
subplot(5,1,4);
hist(mod_pred.bphi(c1,3), -0.000005:0.0000005:0.0001)
title('Type1', 'fontsize',20)
subplot(5,1,5);
hist(mod_pred.bphi(c2,3), -0.000005:0.0000005:0.0001)
title('Type2', 'fontsize',20)

%% For lgphi
c0 = find(mod_pred.lgphi(:,2) == 0);
c1 = find(mod_pred.lgphi(:,2) == 1);
c2 = find(mod_pred.lgphi(:,2) == 2);
figure
graphTitle = 'lgphi';
g = subplot(5,1,1);
t = title(graphTitle,'FontSize',40);
g.Visible = 'off';
t.Visible = 'on';
subplot(5,1,2);
hist(mod_pred.lgphi(:,3), -0.000005:0.0000005:0.0001)
title('All cells', 'fontsize',20)
subplot(5,1,3);
hist(mod_pred.lgphi(c0,3), -0.000005:0.0000005:0.0001)
title('Type0', 'fontsize',20)
subplot(5,1,4);
hist(mod_pred.lgphi(c1,3), -0.000005:0.0000005:0.0001)
title('Type1', 'fontsize',20)
subplot(5,1,5);
hist(mod_pred.lgphi(c2,3), -0.000005:0.0000005:0.0001)
title('Type2', 'fontsize',20)

%% For hgphi
c0 = find(mod_pred.hgphi(:,2) == 0);
c1 = find(mod_pred.hgphi(:,2) == 1);
c2 = find(mod_pred.hgphi(:,2) == 2);
figure
graphTitle = 'hgphi';
g = subplot(5,1,1);
t = title(graphTitle,'FontSize',40);
g.Visible = 'off';
t.Visible = 'on';
subplot(5,1,2);
hist(mod_pred.hgphi(:,3), -0.000005:0.0000005:0.0001)
title('All cells', 'fontsize',20)
subplot(5,1,3);
hist(mod_pred.hgphi(c0,3), -0.000005:0.0000005:0.0001)
title('Type0', 'fontsize',20)
subplot(5,1,4);
hist(mod_pred.hgphi(c1,3), -0.000005:0.0000005:0.0001)
title('Type1', 'fontsize',20)
subplot(5,1,5);
hist(mod_pred.hgphi(c2,3), -0.000005:0.0000005:0.0001)
title('Type2', 'fontsize',20)

%% Box Plots
C3 = mod_pred.dphi(:,3)';
C0 = mod_pred.dphi(c0,3)';
C1 = mod_pred.dphi(c1,3)';
C2 = mod_pred.dphi(c2,3)';
C = [C3 C0 C1 C2];
grp = [repmat("All cells", 1, length(C3)),repmat("Type0", 1, length(C0)),repmat("Type1", 1, length(C1)),repmat("Type2", 1, length(C2))];
figure
boxplot(C,grp);
title('dphi', 'Fontsize', 40)

C3 = mod_pred.tphi(:,3)';
C0 = mod_pred.tphi(c0,3)';
C1 = mod_pred.tphi(c1,3)';
C2 = mod_pred.tphi(c2,3)';
C = [C3 C0 C1 C2];
grp = [repmat("All cells", 1, length(C3)),repmat("Type0", 1, length(C0)),repmat("Type1", 1, length(C1)),repmat("Type2", 1, length(C2))];
figure
boxplot(C,grp);
title('tphi', 'Fontsize', 40)

C3 = mod_pred.bphi(:,3)';
C0 = mod_pred.bphi(c0,3)';
C1 = mod_pred.bphi(c1,3)';
C2 = mod_pred.bphi(c2,3)';
C = [C3 C0 C1 C2];
grp = [repmat("All cells", 1, length(C3)),repmat("Type0", 1, length(C0)),repmat("Type1", 1, length(C1)),repmat("Type2", 1, length(C2))];
figure
boxplot(C,grp);
title('bphi', 'Fontsize', 40)

C3 = mod_pred.lgphi(:,3)';
C0 = mod_pred.lgphi(c0,3)';
C1 = mod_pred.lgphi(c1,3)';
C2 = mod_pred.lgphi(c2,3)';
C = [C3 C0 C1 C2];
grp = [repmat("All cells", 1, length(C3)),repmat("Type0", 1, length(C0)),repmat("Type1", 1, length(C1)),repmat("Type2", 1, length(C2))];
figure
boxplot(C,grp);
title('lgphi', 'Fontsize', 40)

C3 = mod_pred.hgphi(:,3)';
C0 = mod_pred.hgphi(c0,3)';
C1 = mod_pred.hgphi(c1,3)';
C2 = mod_pred.hgphi(c2,3)';
C = [C3 C0 C1 C2];
grp = [repmat("All cells", 1, length(C3)),repmat("Type0", 1, length(C0)),repmat("Type1", 1, length(C1)),repmat("Type2", 1, length(C2))];
figure
boxplot(C,grp);
title('hgphi', 'Fontsize', 40)