%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Analysis and plotting for the experiments results 
% (see publication with doi: "add doi here")
% author: Barbara Genocchi 
% date: 06/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

cd('/Volumes/LaCie/NanoCom_github/Experiments/analysis/')
%control 28DIV
%NS
list_mea =  [18154, 18148, 18136, 3255, 18142, 38785];

%pre
SR_list_NS = [];
BR_list_NS = [];
centrality_degree_NS = [];
L_MEA_path = [];
num_CH = 60;
weight_NS = {};
node_NS = {};
SR_NS = {};
BR_NS = {};
weight_NS_list = [];
degree_NS_list = [];

for mm = 1:length(list_mea)
    load_file = ['graphAnalysis_' num2str(list_mea(mm)) '.mat'];
    
    if exist(load_file) == 0
        continue
    else
        load(load_file)
        if isempty(conn_matrix_NS{mm,2}) == 1
            continue
        else
            load_file_2 = ['/Volumes/LaCie/NanoCom_github/Experiments/data/BurstAnalysis/NS/' num2str(list_mea(mm)) '/burst/burstAnalysis.mat'];
            load(load_file_2)
            BR_list_NS = [BR_list_NS; cell2mat(BR)];
            
            SR_list_NS = [SR_list_NS; cell2mat(SR)];
            
            BR_mat = cell2mat(BR);
            BR_mean = mean(BR_mat(BR_mat > 0));
            BR_NS{mm,1} = list_mea(mm);
            BR_NS{mm,2} = BR_mean;
            
            SR_mat = cell2mat(SR);
            SR_mean = mean(SR_mat(SR_mat > 0));
            SR_NS{mm,1} = list_mea(mm);
            SR_NS{mm,2} = SR_mean;
            weight_NS_for = [];
            for ii = 2: num_CH
                for jj = 1: ii-1
                    weight_NS_for = [weight_NS_for; conn_matrix_NS{mm,2}(ii,jj)];
                    
                end
            end
            weight_NS{mm,1} = list_mea(mm);
            weight_NS{mm,2} = mean(weight_NS_for);
            weight_NS_list = [weight_NS_list; weight_NS_for];
            centrality_degree = [];
            centrality_degree = centrality(graph_NS{mm,2}, 'degree', 'Importance',graph_NS{mm,2}.Edges.Weight);
            centrality_degree_NS = [centrality_degree_NS; centrality_degree];
            degree_NS_mat = cell2mat(degree_NS(mm,2));
            degree_NS_list = [degree_NS_list; degree_NS_mat];
            node_NS{mm,1} = list_mea(mm);
            node_NS{mm,2} = length(degree_NS_mat(degree_NS_mat > 0));
        end
    end
end

load_file = ['graphAnalysis_' num2str(list_mea(6)) '.mat'];
load(load_file)


degree_NS = [cell2mat(degree_NS(1, 3)); cell2mat(degree_NS(2, 3)); cell2mat(degree_NS(3, 3)); cell2mat(degree_NS(4, 3)); cell2mat(degree_NS(5, 3)); cell2mat(degree_NS(6, 3))];
edges_NS = [cell2mat(edges_NS(1, 2)); cell2mat(edges_NS(2, 2)); cell2mat(edges_NS(3,2)); cell2mat(edges_NS(4,2)); cell2mat(edges_NS(5, 2)); cell2mat(edges_NS(6,2))];
L_min_NS = [cell2mat(L_min_NS(1,3)); cell2mat(L_min_NS(2,3)); cell2mat(L_min_NS(3, 3)); cell2mat(L_min_NS(4, 3)); cell2mat(L_min_NS(5,3)); cell2mat(L_min_NS(6, 3))];
centrality_degree_NS_list = centrality_degree_NS;
nodes_NS = [cell2mat(node_NS(1, 2)); cell2mat(node_NS(2, 2)); cell2mat(node_NS(3, 2)); cell2mat(node_NS(4, 2)); cell2mat(node_NS(5, 2)); cell2mat(node_NS(6, 2))];
weights_NS = [cell2mat(weight_NS(1, 2)); cell2mat(weight_NS(2, 2)); cell2mat(weight_NS(3, 2)); cell2mat(weight_NS(4, 2)); cell2mat(weight_NS(5, 2)); cell2mat(weight_NS(6, 2))];
SR_mean_NS = [cell2mat(SR_NS(1, 2)); cell2mat(SR_NS(2, 2)); cell2mat(SR_NS(3, 2)); cell2mat(SR_NS(4, 2)); cell2mat(SR_NS(5, 2)); cell2mat(SR_NS(6, 2))];
BR_mean_NS = [cell2mat(BR_NS(1, 2)); cell2mat(BR_NS(2, 2)); cell2mat(BR_NS(3, 2)); cell2mat(BR_NS(4, 2)); cell2mat(BR_NS(5, 2)); cell2mat(BR_NS(6, 2))];

%control 28DIV
%9010
list_mea =  [3247, 18139, 18144, 38780, 18140, 18151];

%pre
SR_list_9010 = [];
BR_list_9010 = [];
centrality_degree_9010 = [];
L_MEA_path = [];
num_CH = 60;
weight_9010 = {};
node_9010 = {};
SR_9010 = {};
BR_9010 = {};
weight_9010_list = [];
degree_9010_list = [];

for mm = 1:length(list_mea)
    load_file = ['graphAnalysis_' num2str(list_mea(mm)) '.mat'];
    if exist(load_file) == 0
        continue
    else
        load(load_file)
        if isempty(conn_matrix_9010{mm,2}) == 1
            continue
        else    load_file_2 = ['/Volumes/LaCie/NanoCom_github/Experiments/data/BurstAnalysis/9010/' num2str(list_mea(mm)) '/burst/burstAnalysis.mat'];
            load(load_file_2)
            BR_list_9010 = [BR_list_9010; cell2mat(BR)];
            SR_list_9010 = [SR_list_9010; cell2mat(SR)];
            BR_mat = cell2mat(BR);
            BR_mean = mean(BR_mat(BR_mat > 0));
            BR_9010{mm,1} = list_mea(mm);
            BR_9010{mm,2} = BR_mean;
            
            SR_mat = cell2mat(SR);
            SR_mean = mean(SR_mat(SR_mat > 0));
            SR_9010{mm,1} = list_mea(mm);
            SR_9010{mm,2} = SR_mean;
            weight_9010_for = [];
            for ii = 2: num_CH
                for jj = 1: ii-1
                    weight_9010_for = [weight_9010_for; conn_matrix_9010{mm,2}(ii,jj)];
                    
                end
            end
            weight_9010{mm,1} = list_mea(mm);
            weight_9010{mm,2} = mean(weight_9010_for);
            weight_9010_list = [weight_9010_list; weight_9010_for];
            centrality_degree = [];
            centrality_degree = centrality(graph_9010{mm,2}, 'degree', 'Importance',graph_9010{mm,2}.Edges.Weight);
            centrality_degree_9010 = [centrality_degree_9010; centrality_degree];
            degree_9010_mat = cell2mat(degree_9010(mm,2));
            degree_9010_list = [degree_9010_list; degree_9010_mat];
            node_9010{mm,1} = list_mea(mm);
            node_9010{mm,2} = length(degree_9010_mat(degree_9010_mat > 0));
        end
    end
end

load_file = ['graphAnalysis_' num2str(list_mea(6)) '.mat'];
load(load_file)


degree_9010 = [cell2mat(degree_9010(1, 3)); cell2mat(degree_9010(2, 3)); cell2mat(degree_9010(3, 3)); cell2mat(degree_9010(4, 3)); cell2mat(degree_9010(5, 3)); cell2mat(degree_9010(6, 3))];
edges_9010 = [cell2mat(edges_9010(1, 2)); cell2mat(edges_9010(2,2)); cell2mat(edges_9010(3,2)); cell2mat(edges_9010(4,2)); cell2mat(edges_9010(5,2)); cell2mat(edges_9010(6,2))];
L_min_9010 = [cell2mat(L_min_9010(1,3)); cell2mat(L_min_9010(2, 3)); cell2mat(L_min_9010(3, 3)); cell2mat(L_min_9010(4, 3)); cell2mat(L_min_9010(5, 3)); cell2mat(L_min_9010(6, 3))];
centrality_degree_9010_list = centrality_degree_9010;
nodes_9010 = [cell2mat(node_9010(1, 2)); cell2mat(node_9010(2, 2)); cell2mat(node_9010(3, 2)); cell2mat(node_9010(4, 2)); cell2mat(node_9010(5, 2)); cell2mat(node_9010(6, 2))];
weights_9010 = [cell2mat(weight_9010(1, 2)); cell2mat(weight_9010(2, 2)); cell2mat(weight_9010(3, 2)); cell2mat(weight_9010(4, 2)); cell2mat(weight_9010(5, 2)); cell2mat(weight_9010(6, 2))];
SR_mean_9010 = [cell2mat(SR_9010(1, 2)); cell2mat(SR_9010(2, 2)); cell2mat(SR_9010(3, 2)); cell2mat(SR_9010(4, 2)); cell2mat(SR_9010(5, 2)); cell2mat(SR_9010(6, 2))];
BR_mean_9010 = [cell2mat(BR_9010(1, 2)); cell2mat(BR_9010(2, 2)); cell2mat(BR_9010(3, 2)); cell2mat(BR_9010(4, 2)); cell2mat(BR_9010(5, 2)); cell2mat(BR_9010(6, 2))];

%control 28DIV
%8020
list_mea =  [18134, 18135, 18137, 18143, 18127, 18126, 38788, 18145];

%pre
SR_list_8020 = [];
BR_list_8020 = [];
centrality_degree_8020 = [];
L_MEA_path = [];
num_CH = 60;
weight_8020 = {};
node_8020 = {};
SR_8020 = {};
BR_8020 = {};
weight_8020_list = [];
degree_8020_list = [];

for mm = 1:length(list_mea)
    load_file = ['graphAnalysis_' num2str(list_mea(mm)) '.mat'];
    if exist(load_file) == 0
        continue
    else
        load(load_file)
        if isempty(conn_matrix_8020{mm,2}) == 1
            continue
        else    load_file_2 = ['/Volumes/LaCie/NanoCom_github/Experiments/data/BurstAnalysis/8020/' num2str(list_mea(mm)) '/burst/burstAnalysis.mat'];
            load(load_file_2)
            BR_list_8020 = [BR_list_8020; cell2mat(BR)];
            SR_list_8020 = [SR_list_8020; cell2mat(SR)];
            BR_8020{mm,1} = list_mea(mm);
            BR_8020{mm,2} = BR_mean;
            
            SR_mat = cell2mat(SR);
            SR_mean = mean(SR_mat(SR_mat > 0));
            SR_8020{mm,1} = list_mea(mm);
            SR_8020{mm,2} = SR_mean;
            weight_8020_for = [];
            for ii = 2: num_CH
                for jj = 1: ii-1
                    weight_8020_for = [weight_8020_for; conn_matrix_8020{mm,2}(ii,jj)];
                    
                end
            end
            weight_8020{mm,1} = list_mea(mm);
            weight_8020{mm,2} = mean(weight_8020_for);
            weight_8020_list = [weight_8020_list; weight_8020_for];
            centrality_degree = [];
            centrality_degree = centrality(graph_8020{mm,2}, 'degree', 'Importance',graph_8020{mm,2}.Edges.Weight);
            centrality_degree_8020 = [centrality_degree_8020; centrality_degree];
            degree_8020_mat = cell2mat(degree_8020(mm,2));
            degree_8020_list = [degree_8020_list; degree_8020_mat];
            node_8020{mm,1} = list_mea(mm);
            node_8020{mm,2} = length(degree_8020_mat(degree_8020_mat > 0));
        end
    end
end

load_file = ['graphAnalysis_' num2str(list_mea(8)) '.mat'];
load(load_file)


degree_8020 = [cell2mat(degree_8020(1, 3)); cell2mat(degree_8020(2, 3)); cell2mat(degree_8020(3, 3)); cell2mat(degree_8020(4, 3)); cell2mat(degree_8020(5, 3)); cell2mat(degree_8020(6, 3)); cell2mat(degree_8020(7, 3)); cell2mat(degree_8020(8, 3))];
edges_8020 = [cell2mat(edges_8020(1, 2)); cell2mat(edges_8020(2,2)); cell2mat(edges_8020(3,2)); cell2mat(edges_8020(4,2)); cell2mat(edges_8020(5,2)); cell2mat(edges_8020(6,2)); cell2mat(edges_8020(7,2)); cell2mat(edges_8020(8,2))];
L_min_8020 = [cell2mat(L_min_8020(1,3)); cell2mat(L_min_8020(2, 3)); cell2mat(L_min_8020(3, 3)); cell2mat(L_min_8020(4, 3)); cell2mat(L_min_8020(5, 3)); cell2mat(L_min_8020(6, 3)); cell2mat(L_min_8020(7, 3)); cell2mat(L_min_8020(8, 3))];
centrality_degree_8020_list = centrality_degree_8020;
nodes_8020 = [cell2mat(node_8020(1, 2)); cell2mat(node_8020(2, 2)); cell2mat(node_8020(3, 2)); cell2mat(node_8020(4, 2)); cell2mat(node_8020(5, 2)); cell2mat(node_8020(6, 2)); cell2mat(node_8020(7, 2)); cell2mat(node_8020(8, 2))];
weights_8020 = [cell2mat(weight_8020(1, 2)); cell2mat(weight_8020(2, 2)); cell2mat(weight_8020(3, 2)); cell2mat(weight_8020(4, 2)); cell2mat(weight_8020(5, 2)); cell2mat(weight_8020(6, 2)); cell2mat(weight_8020(7, 2)); cell2mat(weight_8020(8, 2))];
SR_mean_8020 = [cell2mat(SR_8020(1, 2)); cell2mat(SR_8020(2, 2)); cell2mat(SR_8020(3, 2)); cell2mat(SR_8020(4, 2)); cell2mat(SR_8020(5, 2)); cell2mat(SR_8020(6, 2)); cell2mat(SR_8020(7, 2)); cell2mat(SR_8020(8, 2))];
BR_mean_8020 = [cell2mat(BR_8020(1, 2)); cell2mat(BR_8020(2, 2)); cell2mat(BR_8020(3, 2)); cell2mat(BR_8020(4, 2)); cell2mat(BR_8020(5, 2)); cell2mat(BR_8020(6, 2)); cell2mat(BR_8020(7, 2)); cell2mat(BR_8020(8, 2))];

%control 28DIV
%7030
list_mea =  [18125, 18152, 18130, 18129, 18133, 18128, 18138, 18155];

%pre
SR_list_7030 = [];
BR_list_7030 = [];
centrality_degree_7030 = [];
L_MEA_path = [];
num_CH = 60;
weight_7030 = {};
node_7030 = {};
SR_7030 = {};
BR_7030 = {};
weight_7030_list = [];
degree_7030_list = [];

for mm = 1:length(list_mea)
    load_file = ['graphAnalysis_' num2str(list_mea(mm)) '.mat'];
    if exist(load_file) == 0
        continue
    else
        load(load_file)
        if isempty(conn_matrix_7030{mm,2}) == 1
            continue
        else    load_file_2 = ['/Volumes/LaCie/NanoCom_github/Experiments/data/BurstAnalysis/7030/' num2str(list_mea(mm)) '/burst/burstAnalysis.mat'];
            load(load_file_2)
            BR_list_7030 = [BR_list_7030; cell2mat(BR)];
            SR_list_7030 = [SR_list_7030; cell2mat(SR)];
            BR_mat = cell2mat(BR);
            BR_mean = mean(BR_mat(BR_mat > 0));
            BR_7030{mm,1} = list_mea(mm);
            BR_7030{mm,2} = BR_mean;
            
            SR_mat = cell2mat(SR);
            SR_mean = mean(SR_mat(SR_mat > 0));
            SR_7030{mm,1} = list_mea(mm);
            SR_7030{mm,2} = SR_mean;
            weight_7030_for = [];
            for ii = 2: num_CH
                for jj = 1: ii-1
                    weight_7030_for = [weight_7030_for; conn_matrix_7030{mm,2}(ii,jj)];
                    
                end
            end
            weight_7030{mm,1} = list_mea(mm);
            weight_7030{mm,2} = mean(weight_7030_for);
            weight_7030_list = [weight_7030_list; weight_7030_for];
            centrality_degree = [];
            centrality_degree = centrality(graph_7030{mm,2}, 'degree', 'Importance',graph_7030{mm,2}.Edges.Weight);
            centrality_degree_7030 = [centrality_degree_7030; centrality_degree];
            degree_7030_mat = cell2mat(degree_7030(mm,2));
            degree_7030_list = [degree_7030_list; degree_7030_mat];
            node_7030{mm,1} = list_mea(mm);
            node_7030{mm,2} = length(degree_7030_mat(degree_7030_mat > 0));
        end
    end
end

load_file = ['graphAnalysis_' num2str(list_mea(8)) '.mat'];
load(load_file)


degree_7030 = [cell2mat(degree_7030(1, 3)); cell2mat(degree_7030(2, 3)); cell2mat(degree_7030(3, 3)); cell2mat(degree_7030(4, 3)); cell2mat(degree_7030(5, 3)); cell2mat(degree_7030(6, 3)); cell2mat(degree_7030(7, 3)); cell2mat(degree_7030(8, 3))];
edges_7030 = [cell2mat(edges_7030(1, 2)); cell2mat(edges_7030(2,2)); cell2mat(edges_7030(3,2)); cell2mat(edges_7030(4,2)); cell2mat(edges_7030(5,2)); cell2mat(edges_7030(6,2)); cell2mat(edges_7030(7,2)); cell2mat(edges_7030(8,2))];
L_min_7030 = [cell2mat(L_min_7030(1,3)); cell2mat(L_min_7030(2, 3)); cell2mat(L_min_7030(3, 3)); cell2mat(L_min_7030(4, 3)); cell2mat(L_min_7030(5, 3)); cell2mat(L_min_7030(6, 3)); cell2mat(L_min_7030(7, 3)); cell2mat(L_min_7030(8, 3))];
centrality_degree_7030_list = centrality_degree_7030;
nodes_7030 = [cell2mat(node_7030(1, 2)); cell2mat(node_7030(2, 2)); cell2mat(node_7030(3, 2)); cell2mat(node_7030(4, 2)); cell2mat(node_7030(5, 2)); cell2mat(node_7030(6, 2)); cell2mat(node_7030(7, 2)); cell2mat(node_7030(8, 2))];
weights_7030 = [cell2mat(weight_7030(1, 2)); cell2mat(weight_7030(2, 2)); cell2mat(weight_7030(3, 2)); cell2mat(weight_7030(4, 2)); cell2mat(weight_7030(5, 2)); cell2mat(weight_7030(6, 2)); cell2mat(weight_7030(7, 2)); cell2mat(weight_7030(8, 2))];
SR_mean_7030 = [cell2mat(SR_7030(1, 2)); cell2mat(SR_7030(2, 2)); cell2mat(SR_7030(3, 2)); cell2mat(SR_7030(4, 2)); cell2mat(SR_7030(5, 2)); cell2mat(SR_7030(6, 2)); cell2mat(SR_7030(7, 2)); cell2mat(SR_7030(8, 2))];
BR_mean_7030 = [cell2mat(BR_7030(1, 2)); cell2mat(BR_7030(2, 2)); cell2mat(BR_7030(3, 2)); cell2mat(BR_7030(4, 2)); cell2mat(BR_7030(5, 2)); cell2mat(BR_7030(6, 2)); cell2mat(BR_7030(7, 2)); cell2mat(BR_7030(8, 2))];



L_min_NS = L_min_NS(isnan(L_min_NS) == 0);
L_min_NS = L_min_NS(isinf(L_min_NS) == 0);


L_min_9010 = L_min_9010(isnan(L_min_9010) == 0);
L_min_9010 = L_min_9010(isinf(L_min_9010) == 0);

L_min_8020 = L_min_8020(isnan(L_min_8020) == 0);
L_min_8020 = L_min_8020(isinf(L_min_8020) == 0);


L_min_7030 = L_min_7030(isnan(L_min_7030) == 0);
L_min_7030 = L_min_7030(isinf(L_min_7030) == 0);



[R_NS,P_NS] = corrcoef(SR_list_NS(degree_NS_list>0), degree_NS_list(degree_NS_list>0));
[R_9010,P_9010] = corrcoef(SR_list_9010(degree_9010_list>0), degree_9010_list(degree_9010_list>0));
[R_8020,P_8020] = corrcoef(SR_list_8020(degree_8020_list>0), degree_8020_list(degree_8020_list>0));
[R_7030,P_7030] = corrcoef(SR_list_7030(degree_7030_list>0), degree_7030_list(degree_7030_list>0));



line1 = fitlm(degree_NS_list(degree_NS_list>0),SR_list_NS(degree_NS_list>0));
line2 = fitlm(degree_9010_list(degree_9010_list>0),SR_list_9010(degree_9010_list>0));
line3 = fitlm(degree_8020_list(degree_8020_list>0),SR_list_8020(degree_8020_list>0));
line4 = fitlm(degree_7030_list(degree_7030_list>0),SR_list_7030(degree_7030_list>0));

cd ../../..
save('Workspace_graph_analysis_experiments.mat')

int1 = line1.Coefficients{1,1};
reg1 = line1.Coefficients{2,1};

int2 = line2.Coefficients{1,1};
reg2 = line2.Coefficients{2,1};

int3 = line3.Coefficients{1,1};
reg3 = line3.Coefficients{2,1};

int4 = line4.Coefficients{1,1};
reg4 = line4.Coefficients{2,1};



label_NS = ['NS - R = ' num2str(R_NS(1,2)) ''];
label_9010 = ['90/10 - R = ' num2str(R_9010(1,2)) ''];
label_8020 = ['80/20 - R = ' num2str(R_8020(1,2)) ''];
label_7030 = ['70/30 - R = ' num2str(R_7030(1,2)) ''];
%plot corr

lineNS = int1 + reg1.*[0:40];
line9010 = int2 + reg2.*[0:40];
line8020 = int3 + reg3.*[0:40];
line7030 = int4 + reg4.*[0:40];

f1 = figure(1);

bubblechart( degree_NS_list(:,1), SR_list_NS(:,1),centrality_degree_NS(:,1),'b','DisplayName', label_NS)
hold on
plot(lineNS, '-b', 'LineWidth', 4)
hold on
bubblechart( degree_9010_list(:,1),SR_list_9010(:,1),centrality_degree_9010(:,1),'c','DisplayName', label_9010)
hold on
plot(line9010,'-c', 'LineWidth', 4)
hold on
bubblechart( degree_8020_list(:,1),SR_list_8020(:,1),centrality_degree_8020(:,1),'r','DisplayName', label_8020)
hold on
plot(line8020, '-r', 'LineWidth', 4)
hold on
bubblechart( degree_7030_list(:,1),SR_list_7030(:,1), centrality_degree_7030(:,1),'y','DisplayName', label_7030)
hold on
plot(line7030, '-y', 'LineWidth', 4)
xlim([0, 40])
ylim([0, 800])
xlabel('Degree')
ylabel('Spike rate (spikes/min)')
set(gca,'FontSize',30, 'FontWeight', 'bold')
legend('Orientation', 'horizontal', 'Location', 'northoutside')

%br
% [R_BR_NS,P_BR_NS] = corrcoef(BR_list_NS, degree_NS_list);
% [R_BR_9010,P_BR_9010] = corrcoef(BR_list_9010, degree_9010_list);
% [R_BR_8020,P_BR_8020] = corrcoef(BR_list_8020, degree_8020_list);
% [R_BR_7030,P_BR_7030] = corrcoef(BR_list_7030, degree_7030_list);
% 
% line1 = fitlm(degree_NS_list(degree_NS_list>0),BR_list_NS(degree_NS_list>0));
% line2 = fitlm(degree_9010_list(degree_9010_list>0),BR_list_9010(degree_9010_list>0));
% line3 = fitlm(degree_8020_list(degree_8020_list>0),BR_list_8020(degree_8020_list>0));
% line4 = fitlm(degree_7030_list(degree_7030_list>0),BR_list_7030(degree_7030_list>0));
% 
% int1 = line1.Coefficients{1,1};
% reg1 = line1.Coefficients{2,1};
% 
% int2 = line2.Coefficients{1,1};
% reg2 = line2.Coefficients{2,1};
% 
% int3 = line3.Coefficients{1,1};
% reg3 = line3.Coefficients{2,1};
% 
% int4 = line4.Coefficients{1,1};
% reg4 = line4.Coefficients{2,1};
% 
% lineNS = int1 + reg1.*[0:40];
% line9010 = int2 + reg2.*[0:40];
% line8020 = int3 + reg3.*[0:40];
% line7030 = int4 + reg4.*[0:40];
% 
% label_NS = ['NS - R = ' num2str(R_BR_NS(1,2)) ''];
% label_9010 = ['90/10 - R = ' num2str(R_BR_9010(1,2)) ''];
% label_8020 = ['80/20 - R = ' num2str(R_BR_8020(1,2)) ''];
% label_7030 = ['70/30 - R = ' num2str(R_BR_7030(1,2)) ''];
% %plot corr
% f2 = figure(2);
% 
% bubblechart( degree_NS_list(:,1), BR_list_NS(:,1),centrality_degree_NS(:,1),'b','DisplayName', label_NS)
% hold on
% plot(lineNS, '-b', 'LineWidth', 4)
% hold on
% bubblechart( degree_9010_list(:,1),BR_list_9010(:,1),centrality_degree_9010(:,1),'c','DisplayName', label_9010)
% hold on
% plot(line9010,'-c', 'LineWidth', 4)
% hold on
% bubblechart( degree_8020_list(:,1),BR_list_8020(:,1),centrality_degree_8020(:,1),'r','DisplayName', label_8020)
% hold on
% plot(line8020, '-r', 'LineWidth', 4)
% hold on
% bubblechart( degree_7030_list(:,1),BR_list_7030(:,1), centrality_degree_7030(:,1),'y','DisplayName', label_7030)
% hold on
% plot(line7030, '-y', 'LineWidth', 4)
% xlabel('Burst rate (bursts/min)')
% ylabel('Degree')
% set(gca,'FontSize',30, 'FontWeight', 'bold')
% legend('Orientation', 'horizontal', 'Location', 'northoutside')


%pairwise Mann-Whitney u tests
[p_SR_9010] = ranksum(SR_mean_NS, SR_mean_9010);
[p_SR_8020] = ranksum(SR_mean_NS, SR_mean_8020);
[p_SR_7030] = ranksum(SR_mean_NS, SR_mean_7030);

[p_BR_9010] = ranksum(BR_mean_NS, BR_mean_9010);
[p_BR_8020] = ranksum(BR_mean_NS, BR_mean_8020);
[p_BR_7030] = ranksum(BR_mean_NS, BR_mean_7030);

[p_degree_9010] = ranksum(degree_NS, degree_9010);
[p_degree_8020] = ranksum(degree_NS, degree_8020);
[p_degree_7030] = ranksum(degree_NS, degree_7030);

[p_weights_9010] = ranksum(weights_NS, weights_9010);
[p_weights_8020] = ranksum(weights_NS, weights_8020);
[p_weights_7030] = ranksum(weights_NS, weights_7030);

[p_node_9010] = ranksum(nodes_NS, nodes_9010);
[p_node_8020] = ranksum(nodes_NS, nodes_8020);
[p_node_7030] = ranksum(nodes_NS, nodes_7030);

[p_weight_9010] = ranksum(weight_NS_list, weight_9010_list);
[p_weight_8020] = ranksum(weight_NS_list, weight_8020_list);
[p_weight_7030] = ranksum(weight_NS_list, weight_7030_list);

[p_L_min_9010] = ranksum(L_min_NS, L_min_9010);
[p_L_min_8020] = ranksum(L_min_NS, L_min_8020);
[p_L_min_7030] = ranksum(L_min_NS, L_min_7030);

[p_edges_9010] = ranksum(edges_NS, edges_9010);
[p_edges_8020] = ranksum(edges_NS, edges_8020);
[p_edges_7030] = ranksum(edges_NS, edges_7030);

[p_centrality_degree_9010] = ranksum(centrality_degree_NS, centrality_degree_9010);
[p_centrality_degree_8020] = ranksum(centrality_degree_NS, centrality_degree_8020);
[p_centrality_degree_7030] = ranksum(centrality_degree_NS, centrality_degree_7030);

[p_SR_list_9010] = ranksum(SR_list_NS, SR_list_9010);
[p_SR_list_8020] = ranksum(SR_list_NS, SR_list_8020);
[p_SR_list_7030] = ranksum(SR_list_NS, SR_list_7030);

[p_BR_list_9010] = ranksum(BR_list_NS, BR_list_9010);
[p_BR_list_8020] = ranksum(BR_list_NS, BR_list_8020);
[p_BR_list_7030] = ranksum(BR_list_NS, BR_list_7030);

[p_degree_list_9010] = ranksum(degree_NS_list, degree_9010_list);
[p_degree_list_8020] = ranksum(degree_NS_list, degree_8020_list);
[p_degree_list_7030] = ranksum(degree_NS_list, degree_7030_list);
