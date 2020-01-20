%%Plot the prediction improvement of models across days%%
cd('/Users/manishm/Work/vanDerMeerLab/GLM/temp/Shuffle')
rats = {'R117','R119','R131','R132'};
% models = {'dphi', 'tphi', 'bphi', 'lgphi', 'hgphi', 'allphi'};
models = {'allphi'};
mod_pred = {};
for idx = 1:length(models)
    mod_pred.(models{idx}) = [];
end

for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat('S2Post_',curRat,'*sd.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name);
%        N = find(sd.cellType == 0);
        N = length(sd.cellType);
        for iC = 1:N
            for iM = 1:length(models)
                mod_pred.(models{iM}) = [mod_pred.(models{iM});str2double(curRat(2:end)),sd.cellType(iC),sd.m.(models{iM}).err(iC), mean(sd.m.(models{iM}).shufErr(iC,:))];
            end
        end
    end
end

%% For all models
for iM = 1:length(models)
    c1 = mod_pred.(models{iM})(:,2) == 1;
    c_1 = mod_pred.(models{iM})(c1,3);
    c_1s = mod_pred.(models{iM})(c1,4);
    c_1 = c_1(:);
    c_1s = c_1s(:);
    c_1d = c_1 - c_1s;
    c2 = mod_pred.(models{iM})(:,2) == 2;
    c_2 = mod_pred.(models{iM})(c2,3);
    c_2s = mod_pred.(models{iM})(c2,4);
    c_2 = c_2(:);
    c_2s = c_2s(:);
    c_2d = c_2 - c_2s;
    c_12 = [c_1;c_2];
    c_12s = [c_1s;c_2s];
    c_12d = c_12 - c_12s;
    figure;
    bin_edges = -0.00001:0.0000005:0.00005;
    histogram(c_1s,bin_edges,'FaceColor','red', 'FaceAlpha',0.4);
    hold on
    histogram(c_1,bin_edges,'FaceColor','green', 'FaceAlpha',0.4);
    ax = gca;
    ax.Box = 'off';
    ax.FontSize = 15;
    ax.TickDir = 'out';
    ax.YLim = [0 200];
    output_name = strcat('Post task', {' '}, (models{iM}),' - type 1');
    title(output_name, 'fontsize',30);
    xlabel('Model Improvement', 'FontSize', 15);
    WriteFig(output_name{1})
    figure;
    bin_edges = -0.00001:0.0000005:0.00005;
    histogram(c_2s,bin_edges,'FaceColor','red', 'FaceAlpha',0.4);
    hold on
    histogram(c_2,bin_edges,'FaceColor','green', 'FaceAlpha',0.4);
    ax = gca;
    ax.Box = 'off';
    ax.FontSize = 15;
    ax.TickDir = 'out';
    ax.YLim = [0 80];
    output_name = strcat('Post task', {' '}, (models{iM}),' - type 2');
    title(output_name, 'fontsize',30);
    xlabel('Model Improvement', 'FontSize', 15);
    WriteFig(output_name{1})
    figure;
    bin_edges = -0.00001:0.0000005:0.00005;
    histogram(c_12s,bin_edges,'FaceColor','red', 'FaceAlpha',0.4);
    hold on
    histogram(c_12,bin_edges,'FaceColor','green', 'FaceAlpha',0.4);
    ax = gca;
    ax.Box = 'off';
    ax.FontSize = 15;
    ax.TickDir = 'out';
    ax.YLim = [0 250];
    output_name = strcat('Post task', {' '}, (models{iM}),' - type 1 and 2');
    title(output_name, 'fontsize',30);
    xlabel('Model Improvement', 'FontSize', 15);
    WriteFig(output_name{1})
    figure;
    mask = c_1d > 0;
    bin_edges = -0.00001:0.0000005:0.00005;
    histogram(c_1d(mask),bin_edges,'FaceColor','green');
    hold on;
    histogram(c_1d(~mask),bin_edges,'FaceColor','red');
    hold on;
    pc = 100 * sum(mask)/length(mask);
    text(0.00002,25,strcat('Percentage of improved models :',{' '}, num2str(round(pc,2))),'fontsize',15);
    ax = gca;
    ax.Box = 'off';
    ax.FontSize = 15;
    ax.TickDir = 'out';
    ax.YLim = [0 30];
    output_name = strcat('Post task (Real - shuffle) ', {' '}, (models{iM}),' - type 1');
    title(output_name, 'fontsize',25);
    xlabel('Model Improvement', 'FontSize', 15);
    hold off;
    WriteFig(output_name{1})
    figure;
    mask = c_2d > 0;
    histogram(c_2d(mask),bin_edges,'FaceColor','green');
    hold on;
    histogram(c_2d(~mask),bin_edges,'FaceColor','red');
    hold on;
    pc = 100 * sum(mask)/length(mask);
    text(0.00002,8,strcat('Percentage of improved models :',{' '}, num2str(round(pc,2))),'fontsize',15);
    ax = gca;
    ax.Box = 'off';
    ax.FontSize = 15;
    ax.TickDir = 'out';
    ax.YLim = [0 10];
    output_name = strcat('Post task (Real - shuffle) ', {' '}, (models{iM}),' - type 2');
    title(output_name, 'fontsize',25);
    xlabel('Model Improvement', 'FontSize', 15);
    hold off;
    WriteFig(output_name{1})
    figure;
    mask = c_12d > 0;
    histogram(c_12d(mask),bin_edges,'FaceColor','green');
    hold on;
    histogram(c_12d(~mask),bin_edges,'FaceColor','red');
    hold on;
    pc = 100 * sum(mask)/length(mask);
    text(0.00002,25,strcat('Percentage of improved models :',{' '}, num2str(round(pc,2))),'fontsize',15);
    ax = gca;
    ax.Box = 'off';
    ax.FontSize = 15;
    ax.TickDir = 'out';
    ax.YLim = [0 35];
    output_name = strcat('Post task (Real - shuffle) ', {' '}, (models{iM}),' - type 1 and 2');
    title(output_name, 'fontsize',25);
    xlabel('Model Improvement', 'FontSize', 15);
    hold off;
    WriteFig(output_name{1})
end

% %% Box Plots
% C3 = mod_pred.dphi(:,3)';
% C0 = mod_pred.dphi(c0,3)';
% C1 = mod_pred.dphi(c1,3)';
% C2 = mod_pred.dphi(c2,3)';
% C = [C3 C0 C1 C2];
% grp = [repmat("All cells", 1, length(C3)),repmat("Type0", 1, length(C0)),repmat("Type1", 1, length(C1)),repmat("Type2", 1, length(C2))];
% figure
% boxplot(C,grp);
% title('dphi', 'Fontsize', 40)
% 
% C3 = mod_pred.tphi(:,3)';
% C0 = mod_pred.tphi(c0,3)';
% C1 = mod_pred.tphi(c1,3)';
% C2 = mod_pred.tphi(c2,3)';
% C = [C3 C0 C1 C2];
% grp = [repmat("All cells", 1, length(C3)),repmat("Type0", 1, length(C0)),repmat("Type1", 1, length(C1)),repmat("Type2", 1, length(C2))];
% figure
% boxplot(C,grp);
% title('tphi', 'Fontsize', 40)
% 
% C3 = mod_pred.bphi(:,3)';
% C0 = mod_pred.bphi(c0,3)';
% C1 = mod_pred.bphi(c1,3)';
% C2 = mod_pred.bphi(c2,3)';
% C = [C3 C0 C1 C2];
% grp = [repmat("All cells", 1, length(C3)),repmat("Type0", 1, length(C0)),repmat("Type1", 1, length(C1)),repmat("Type2", 1, length(C2))];
% figure
% boxplot(C,grp);
% title('bphi', 'Fontsize', 40)
% 
% C3 = mod_pred.lgphi(:,3)';
% C0 = mod_pred.lgphi(c0,3)';
% C1 = mod_pred.lgphi(c1,3)';
% C2 = mod_pred.lgphi(c2,3)';
% C = [C3 C0 C1 C2];
% grp = [repmat("All cells", 1, length(C3)),repmat("Type0", 1, length(C0)),repmat("Type1", 1, length(C1)),repmat("Type2", 1, length(C2))];
% figure
% boxplot(C,grp);
% title('lgphi', 'Fontsize', 40)
% 
% C3 = mod_pred.hgphi(:,3)';
% C0 = mod_pred.hgphi(c0,3)';
% C1 = mod_pred.hgphi(c1,3)';
% C2 = mod_pred.hgphi(c2,3)';
% C = [C3 C0 C1 C2];
% grp = [repmat("All cells", 1, length(C3)),repmat("Type0", 1, length(C0)),repmat("Type1", 1, length(C1)),repmat("Type2", 1, length(C2))];
% figure
% boxplot(C,grp);
% title('hgphi', 'Fontsize', 40)