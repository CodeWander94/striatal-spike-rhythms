%%Plot the prediction improvement of models across days%%
cd('/Users/manishm/Work/vanDerMeerLab/GLM/temp/Shuffle')
rats = {'R117','R119','R131','R132'};

fid = fopen('Post_results.csv', 'a+' );
fprintf(fid, 'Rat,session,Total Cells,Duplicate Cells,Few Spikes,Same TT,Type 0,Type 1,Type 2,Type 1 and 2\n');
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat('S2Post_',curRat,'*sd.mat');
    ofiles = dir(searchString);  
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name);
        t0 = sum(find(sd.cellType == 0));
        t1 = sum(find(sd.cellType == 1));
        t2 = sum(find(sd.cellType == 2));
        fS = length(sd.fewSpikes);
        sT = length(sd.sameTT);
        dC = length(sd.duplicateCells);
        q = split(ofiles(jdx).name, '_');
        q = split(q(2), '-');
        session = join(q(2:end),'/');
%         res = [curRat,session,num2str(t0+t1+t2+fS+sT+dC),num2str(dC),num2str(fS), num2str(sT),num2str(t0),num2str(t1),num2str(t2),num2str(t1+t2)];
        fprintf(fid, '%s,%s,%d,%d,%d,%d,%d,%d,%d,%d\n',curRat,cell2mat(session(1)),t0+t1+t2+fS+sT+dC,dC,fS,sT,t0,t1,t2,t1+t2);
    end
end
fclose(fid);