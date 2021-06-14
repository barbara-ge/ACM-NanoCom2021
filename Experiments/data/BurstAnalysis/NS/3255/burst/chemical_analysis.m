%%4AP && gabazine

%control 28DIV
%NS
list_mea = [3255, 18136, 18142, 18147, 18148, 18154, 38785];

for m = 1: length(list_mea)
    
    mea_num = double(list_mea(m));
    path_cd = ['/media/barbara/LaCie/MEA_data/winter2021/17022021_week4/CSV/NS/analysed/' num2str(mea_num) '/burst/'];
    cd(path_cd)
    load('burstAnalysis.mat')
    
    BurstInfoNArray = table2cell(BurstInfoN);
    lengthST = DataCell{k, 4} * 1000; %from s to ms
    
    SR_mat = cell2mat(BurstInfoNArray(:,2))/(lengthST/1000/60);
    BR_mat = cell2mat(BurstInfoNArray(:,4))/(lengthST/1000/60);
    BD_mat = cell2mat(BurstInfoNArray(:,6));
    SB_mat = cell2mat(BurstInfoNArray(:,7));
    SpikeType = cell2mat(BurstInfoNArray(m,13));
    BurstSpike = SpikeType(SpikeType==1);
    SpikeCount = cell2mat(BurstInfoNArray(m,2));
    
    NS_SR(m,1) = mean(SR_mat(SR_mat>3));
    NS_BR(m,1) = mean(BR_mat(BR_mat>0));
    NS_BD(m,1) = mean(BD_mat(BD_mat>0));
    NS_SB(m,1) = mean(SB_mat(SB_mat>0));
    NS_SBPerc(m,1) = length(BurstSpike)/SpikeCount;
    
    
end

%9010
list_mea = [3247, 18139, 18140, 18141, 18144, 18151, 38776, 38780];

for m = 1: length(list_mea)
    
    mea_num = double(list_mea(m));
    path_cd = ['/media/barbara/LaCie/MEA_data/winter2021/17022021_week4/CSV/9010/analysed/' num2str(mea_num) '/burst/'];
    cd(path_cd)
    load('burstAnalysis.mat')
    
    BurstInfoNArray = table2cell(BurstInfoN);
        lengthST = DataCell{k, 4} * 1000; %from s to ms

    SR_mat = cell2mat(BurstInfoNArray(:,2))/(lengthST/1000/60);
    BR_mat = cell2mat(BurstInfoNArray(:,4))/(lengthST/1000/60);
    BD_mat = cell2mat(BurstInfoNArray(:,6));
    SB_mat = cell2mat(BurstInfoNArray(:,7));
    SpikeType = cell2mat(BurstInfoNArray(m,13));
    BurstSpike = SpikeType(SpikeType==1);
    SpikeCount = cell2mat(BurstInfoNArray(m,2));
    
    SR_9010(m,1) = mean(SR_mat(SR_mat>3));
    BR_9010(m,1) = mean(BR_mat(BR_mat>0));
    BD_9010(m,1) = mean(BD_mat(BD_mat>0));
    SB_9010(m,1) = mean(SB_mat(SB_mat>0));
    SBPerc_9010(m,1) = length(BurstSpike)/SpikeCount;
    
    
end

%%8020
list_mea = [18126, 18127, 18134, 18135, 18137, 18143, 18145, 38788];

for m = 1: length(list_mea)
    
    mea_num = double(list_mea(m));
    path_cd = ['/media/barbara/LaCie/MEA_data/winter2021/17022021_week4/CSV/8020/analysed/' num2str(mea_num) '/burst/'];
    cd(path_cd)
    load('burstAnalysis.mat')
    
    BurstInfoNArray = table2cell(BurstInfoN);
        lengthST = DataCell{k, 4} * 1000; %from s to ms

    SR_mat = cell2mat(BurstInfoNArray(:,2))/(lengthST/1000/60);
    BR_mat = cell2mat(BurstInfoNArray(:,4))/(lengthST/1000/60);
    BD_mat = cell2mat(BurstInfoNArray(:,6));
    SB_mat = cell2mat(BurstInfoNArray(:,7));
    SpikeType = cell2mat(BurstInfoNArray(m,13));
    BurstSpike = SpikeType(SpikeType==1);
    SpikeCount = cell2mat(BurstInfoNArray(m,2));
    
    SR_8020(m,1) = mean(SR_mat(SR_mat>3));
    BR_8020(m,1) = mean(BR_mat(BR_mat>0));
    BD_8020(m,1) = mean(BD_mat(BD_mat>0));
    SB_8020(m,1) = mean(SB_mat(SB_mat>0));
    SBPerc_8020(m,1) = length(BurstSpike)/SpikeCount;
    
    
end

%%7030
list_mea = [18125, 18128, 18129, 18130, 18133, 18138, 18152, 18155];

for m = 1: length(list_mea)
    
    mea_num = double(list_mea(m));
    path_cd = ['/media/barbara/LaCie/MEA_data/winter2021/17022021_week4/CSV/7030/analysed/' num2str(mea_num) '/burst/'];
    cd(path_cd)
    load('burstAnalysis.mat')
    
    BurstInfoNArray = table2cell(BurstInfoN);
        lengthST = DataCell{k, 4} * 1000; %from s to ms

    SR_mat = cell2mat(BurstInfoNArray(:,2))/(lengthST/1000/60);
    BR_mat = cell2mat(BurstInfoNArray(:,4))/(lengthST/1000/60);
    BD_mat = cell2mat(BurstInfoNArray(:,6));
    SB_mat = cell2mat(BurstInfoNArray(:,7));
    SpikeType = cell2mat(BurstInfoNArray(m,13));
    BurstSpike = SpikeType(SpikeType==1);
    SpikeCount = cell2mat(BurstInfoNArray(m,2));
    
    SR_7030(m,1) = mean(SR_mat(SR_mat>3));
    BR_7030(m,1) = mean(BR_mat(BR_mat>0));
    BD_7030(m,1) = mean(BD_mat(BD_mat>0));
    SB_7030(m,1) = mean(SB_mat(SB_mat>0));
    SBPerc_7030(m,1) = length(BurstSpike)/SpikeCount;
    
    
end

%%5050
list_mea = [18131, 18146, 18150, 38777, 38781, 38782, 38784, 38786];

for m = 1: length(list_mea)
    
    mea_num = double(list_mea(m));
    path_cd = ['/media/barbara/LaCie/MEA_data/winter2021/17022021_week4/CSV/5050/analysed/' num2str(mea_num) '/burst/'];
    cd(path_cd)
    load('burstAnalysis.mat')
    
    BurstInfoNArray = table2cell(BurstInfoN);
        lengthST = DataCell{k, 4} * 1000; %from s to ms

    SR_mat = cell2mat(BurstInfoNArray(:,2))/(lengthST/1000/60);
    BR_mat = cell2mat(BurstInfoNArray(:,4))/(lengthST/1000/60);
    BD_mat = cell2mat(BurstInfoNArray(:,6));
    SB_mat = cell2mat(BurstInfoNArray(:,7));
    SpikeType = cell2mat(BurstInfoNArray(m,13));
    BurstSpike = SpikeType(SpikeType==1);
    SpikeCount = cell2mat(BurstInfoNArray(m,2));
    
    SR_5050(m,1) = mean(SR_mat(SR_mat>3));
    BR_5050(m,1) = mean(BR_mat(BR_mat>0));
    BD_5050(m,1) = mean(BD_mat(BD_mat>0));
    SB_5050(m,1) = mean(SB_mat(SB_mat>0));
    SBPerc_5050(m,1) = length(BurstSpike)/SpikeCount;
    
    
end

fig1 = figure(1);
columnParams = [NS_SR(:,1); SR_9010(:,1); SR_8020(:,1); SR_7030(:,1); SR_5050(:,1)];
columnLabels = [ones(size(NS_SR(:,1))); 2*ones(size(SR_9010(:,1)));3*ones(size(SR_8020(:,1)));4*ones(size(SR_7030(:,1)));5*ones(size(SR_7030(:,1)))];
boxplot(columnParams, columnLabels)
xticklabels({'NS', '9010', '8020', '7030', '5050'})
xlabel('Population')
ylabel('Spike rate (spikes/min)')

fig2 = figure(2);
columnParams = [NS_BR(:,1); BR_9010(:,1); BR_8020(:,1); BR_7030(:,1); BR_5050(:,1)];
columnLabels = [ones(size(NS_BR(:,1))); 2*ones(size(BR_9010(:,1)));3*ones(size(BR_8020(:,1)));4*ones(size(BR_7030(:,1)));5*ones(size(BR_7030(:,1)))];
boxplot(columnParams, columnLabels)
xticklabels({'NS', '9010', '8020', '7030', '5050'})
xlabel('Population')
ylabel('Burst rate (burst/min)')