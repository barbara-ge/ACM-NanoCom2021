%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Graph theory analysis for the experiments results 
% (see publication with doi: "add doi here")
% author: Barbara Genocchi 
% date: 06/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all
addpath '/Volumes/LaCie/NanoCom_github/Experiments/code' %adapt to your path here

%NS
list_mea = [18154, 18148, 18136, 3255, 18142, 38785];

for mm = 1: length(list_mea)
    
    mm
    
    mea_num = double(list_mea(mm));
    path_cd = ['/Volumes/LaCie/NanoCom_github/Experiments/data/']; %change path here to where you saved the data
    cd(path_cd)
    file_load = ['DataCell_NS_' num2str(mea_num) '.mat'];
    
    load(file_load)
    DataCell = DataCell_2;
    num_CH = 60;
    lengthST = cell2mat(DataCell(1,4)); %length of spike train. Do not manually change it
    
    [conn_matrix] = correlation_MEA(DataCell, num_CH, lengthST);
    
    
    IDs = [];
    for hh = 1: num_CH
        if hh == 15
            IDs(hh,1) = 15;
        else
            IDs(hh,1) = str2num(cell2mat(DataCell(hh,1)));
        end
    end
    
    %clear noisy channels
    artefacted_channels = [15]; %change if your channels have artefact, removed only ground here
    for u = 1 : length(artefacted_channels)
        uu = find(IDs == artefacted_channels(1,u));
        conn_matrix(uu , :)=0;
        conn_matrix(: , uu)=0;
    end
    
    max_corr = max(conn_matrix, [],'all');
    conn_matrix_norm = conn_matrix./max_corr;
    conn_matrix_norm_rounded = round(conn_matrix_norm, 2); %if too many decimal places matrix becomes not symmetric due to internal rounding
    
    if issymmetric(conn_matrix_norm_rounded)==0
        continue
    else
        %graph creation from cleaned and normalized correlation matrix
        graph_creation;
        
        conn_matrix_NS{mm,1} = mea_num;
        conn_matrix_NS{mm,2} = conn_matrix_norm_rounded;
        
        graph_NS{mm,1} = mea_num;
        graph_NS{mm,2} = graph_M_th;
        
        %Analysis
        
        %edges
        
        edges_MEA = numedges(graph_M_th);
        
        edges_NS{mm,1} = mea_num;
        edges_NS{mm,2} = edges_MEA;
        
        %degree_MEAree
        degree_MEA = degree(graph_M_th);
        mean_degree_MEA = mean(degree_MEA);
        std_degree_MEA = std(degree_MEA);
        
        degree_NS{mm,1} = mea_num;
        degree_NS{mm,2} = degree_MEA;
        degree_NS{mm,3} = mean_degree_MEA;
        degree_NS{mm,4} = std_degree_MEA;
        
        %shortest path
        L_MEA_min = [];
        
        for bb = 1: num_CH
            for cc = 1: num_CH
                
                if bb~=cc %no self connections
                    [L_path, min_path] = shortestpath(graph_M_th,bb,cc);
                    
                    L_MEA_min = [L_MEA_min, min_path];
                end
            end
        end
        L_MEA_mean = mean(L_MEA_min(find(L_MEA_min ~= Inf)));
        L_MEA_std = std(L_MEA_min(find(L_MEA_min ~= Inf)));
        
        L_min_NS{mm,1} = mea_num;
        L_min_NS{mm,2} = L_MEA_min;
        L_min_NS{mm,3} = L_MEA_mean;
        L_min_NS{mm,4} = L_MEA_std;
        
        
        for i = 1: num_CH
            nodeID(i,1) = str2num(ids{i});
        end
        
        %plot of MEA graphs
        %plot_graph; %plot graph if necessary
        
        %save results
        cd('/Volumes/LaCie/NanoCom_github/Experiments/analysis/')  %adapt to your path here
        
        save_name = ['graphAnalysis_' num2str(mea_num) ''];
        save(save_name, 'SR', 'conn_matrix_NS', 'graph_NS', 'degree_NS', 'edges_NS', 'L_min_NS')
        
     end
end
%9010
list_mea = [3247, 18139, 18144, 38780, 18140, 18151];

for mm = 1: length(list_mea)
    mm
    mea_num = double(list_mea(mm));
    path_cd = ['/Volumes/LaCie/NanoCom_github/Experiments/data/']; %adapt to your path here
    cd(path_cd)
    file_load = ['DataCell_9010_' num2str(mea_num) '.mat'];
    
    load(file_load)
    DataCell = DataCell_2;
    num_CH = 60;
    lengthST = cell2mat(DataCell(1,4));
    
    [conn_matrix] = correlation_MEA(DataCell, num_CH, lengthST);
    
    
    IDs = [];
    for hh = 1: num_CH
        if hh == 15
            IDs(hh,1) = 15;
        else
            IDs(hh,1) = str2num(cell2mat(DataCell(hh,1)));
        end
    end
    
    %clear noisy channels
    artefacted_channels = [15];
    for u = 1 : length(artefacted_channels)
        uu = find(IDs == artefacted_channels(1,u));
        conn_matrix(uu , :)=0;
        conn_matrix(: , uu)=0;
    end
    
    max_corr = max(conn_matrix, [],'all');
    conn_matrix_norm = conn_matrix./max_corr;
    conn_matrix_norm_rounded = round(conn_matrix_norm, 2); %if too many decimal places matrix becomes not symmetric due to internal rounding
    
    if issymmetric(conn_matrix_norm_rounded)==0
        continue
    else
        %graph creation from cleaned and normalized correlation matrix
        graph_creation;
        
        conn_matrix_9010{mm,1} = mea_num;
        conn_matrix_9010{mm,2} = conn_matrix_norm_rounded;
        
        graph_9010{mm,1} = mea_num;
        graph_9010{mm,2} = graph_M_th;
        
        %Analysis
        
        %edges
        
        edges_MEA = numedges(graph_M_th);
        
        edges_9010{mm,1} = mea_num;
        edges_9010{mm,2} = edges_MEA;
        
        %degree_MEAree
        degree_MEA = degree(graph_M_th);
        mean_degree_MEA = mean(degree_MEA);
        std_degree_MEA = std(degree_MEA);
        
        degree_9010{mm,1} = mea_num;
        degree_9010{mm,2} = degree_MEA;
        degree_9010{mm,3} = mean_degree_MEA;
        degree_9010{mm,4} = std_degree_MEA;
        
        %shortest path
        L_MEA_min = [];
        
        for bb = 1: num_CH
            for cc = 1: num_CH
                
                if bb~=cc %no self connections
                    [L_path, min_path] = shortestpath(graph_M_th,bb,cc);
                    
                    L_MEA_min = [L_MEA_min, min_path];
                end
            end
        end
        L_MEA_mean = mean(L_MEA_min(find(L_MEA_min ~= Inf)));
        L_MEA_std = std(L_MEA_min(find(L_MEA_min ~= Inf)));
        
        L_min_9010{mm,1} = mea_num;
        L_min_9010{mm,2} = L_MEA_min;
        L_min_9010{mm,3} = L_MEA_mean;
        L_min_9010{mm,4} = L_MEA_std;
        
                
        for i = 1: num_CH
            nodeID(i,1) = str2num(ids{i});
        end
        
        %plot of MEA graphs
        %plot_graph;
        
        %save results
        cd('/Volumes/LaCie/NanoCom_github/Experiments/analysis/') %adapt to your path here
        
        save_name = ['graphAnalysis_' num2str(mea_num) ''];
        save(save_name, 'SR', 'conn_matrix_9010','graph_9010', 'degree_9010', 'edges_9010', 'L_min_9010')
        
        
    end
end
%8020
list_mea = [18134, 18135, 18137, 18143, 18127, 18126, 38788, 18145];

for mm = 1: length(list_mea)
    mm
    mea_num = double(list_mea(mm));
    path_cd = ['/Volumes/LaCie/NanoCom_github/Experiments/data/']; %adapt to your path here
    cd(path_cd)
    file_load = ['DataCell_8020_' num2str(mea_num) '.mat'];
    
    load(file_load)
    DataCell = DataCell_2;
    num_CH = 60;
    lengthST = cell2mat(DataCell(1,4));
    
    [conn_matrix] = correlation_MEA(DataCell, num_CH, lengthST);
    
    
    IDs = [];
    for hh = 1: num_CH
        if hh == 15
            IDs(hh,1) = 15;
        else
            IDs(hh,1) = str2num(cell2mat(DataCell(hh,1)));
        end
    end
    
    %clear noisy channels
    artefacted_channels = [15];
    for u = 1 : length(artefacted_channels)
        uu = find(IDs == artefacted_channels(1,u));
        conn_matrix(uu , :)=0;
        conn_matrix(: , uu)=0;
    end
    
    max_corr = max(conn_matrix, [],'all');
    conn_matrix_norm = conn_matrix./max_corr;
    conn_matrix_norm_rounded = round(conn_matrix_norm, 2); %if too many decimal places matrix becomes not symmetric due to internal rounding
    
    if issymmetric(conn_matrix_norm_rounded)==0
        continue
    else
        %graph creation from cleaned and normalized correlation matrix
        graph_creation;
        
        conn_matrix_8020{mm,1} = mea_num;
        conn_matrix_8020{mm,2} = conn_matrix_norm_rounded;
        
        graph_8020{mm,1} = mea_num;
        graph_8020{mm,2} = graph_M_th;
        
        %Analysis
        
        %edges
        
        edges_MEA = numedges(graph_M_th);
        
        edges_8020{mm,1} = mea_num;
        edges_8020{mm,2} = edges_MEA;
        
        %degree_MEAree
        degree_MEA = degree(graph_M_th);
        mean_degree_MEA = mean(degree_MEA);
        std_degree_MEA = std(degree_MEA);
        
        degree_8020{mm,1} = mea_num;
        degree_8020{mm,2} = degree_MEA;
        degree_8020{mm,3} = mean_degree_MEA;
        degree_8020{mm,4} = std_degree_MEA;
        
        %shortest path
        L_MEA_min = [];
        
        for bb = 1: num_CH
            for cc = 1: num_CH
                
                if bb~=cc %no self connections
                    [L_path, min_path] = shortestpath(graph_M_th,bb,cc);
                    
                    L_MEA_min = [L_MEA_min, min_path];
                end
            end
        end
        L_MEA_mean = mean(L_MEA_min(find(L_MEA_min ~= Inf)));
        L_MEA_std = std(L_MEA_min(find(L_MEA_min ~= Inf)));
        
        L_min_8020{mm,1} = mea_num;
        L_min_8020{mm,2} = L_MEA_min;
        L_min_8020{mm,3} = L_MEA_mean;
        L_min_8020{mm,4} = L_MEA_std;
        
        
        
        for i = 1: num_CH
            nodeID(i,1) = str2num(ids{i});
        end
        
        %plot of MEA graphs
        %plot_graph;
        
        %save results
        cd('/Volumes/LaCie/NanoCom_github/Experiments/analysis/') %adapt to your path here
        
        save_name = ['graphAnalysis_' num2str(mea_num) ''];
        save(save_name, 'SR','conn_matrix_8020','graph_8020', 'degree_8020', 'edges_8020', 'L_min_8020')
        
        
    end
end
%7030
list_mea = [18125, 18152, 18130, 18129, 18133, 18128, 18138, 18155];

for mm = 1: length(list_mea)
    mm
    mea_num = double(list_mea(mm));
    path_cd = ['/Volumes/LaCie/NanoCom_github/Experiments/data/']; %adapt to your path here
    cd(path_cd)
    file_load = ['DataCell_7030_' num2str(mea_num) '.mat'];
    
    load(file_load)
    DataCell = DataCell_2;
    num_CH = 60;
    lengthST = cell2mat(DataCell(1,4));
    
    [conn_matrix] = correlation_MEA(DataCell, num_CH, lengthST);
    
    
    IDs = [];
    for hh = 1: num_CH
        if hh == 15
            IDs(hh,1) = 15;
        else
            IDs(hh,1) = str2num(cell2mat(DataCell(hh,1)));
        end
    end
    
    %clear noisy channels
    artefacted_channels = [15];
    for u = 1 : length(artefacted_channels)
        uu = find(IDs == artefacted_channels(1,u));
        conn_matrix(uu , :)=0;
        conn_matrix(: , uu)=0;
    end
    
    max_corr = max(conn_matrix, [],'all');
    conn_matrix_norm = conn_matrix./max_corr;
    conn_matrix_norm_rounded = round(conn_matrix_norm, 2); %if too many decimal places matrix becomes not symmetric due to internal rounding
    
    if issymmetric(conn_matrix_norm_rounded)==0
        continue
    else
        %graph creation from cleaned and normalized correlation matrix
        graph_creation;
        
        conn_matrix_7030{mm,1} = mea_num;
        conn_matrix_7030{mm,2} = conn_matrix_norm_rounded;
        
        graph_7030{mm,1} = mea_num;
        graph_7030{mm,2} = graph_M_th;
        
        %Analysis
        
        %edges
        
        edges_MEA = numedges(graph_M_th);
        
        edges_7030{mm,1} = mea_num;
        edges_7030{mm,2} = edges_MEA;
        
        %degree_MEAree
        degree_MEA = degree(graph_M_th);
        mean_degree_MEA = mean(degree_MEA);
        std_degree_MEA = std(degree_MEA);
        
        degree_7030{mm,1} = mea_num;
        degree_7030{mm,2} = degree_MEA;
        degree_7030{mm,3} = mean_degree_MEA;
        degree_7030{mm,4} = std_degree_MEA;
        
        %shortest path
        L_MEA_min = [];
        
        for bb = 1: num_CH
            for cc = 1: num_CH
                
                if bb~=cc %no self connections
                    [L_path, min_path] = shortestpath(graph_M_th,bb,cc);
                    
                    L_MEA_min = [L_MEA_min, min_path];
                end
            end
        end
        L_MEA_mean = mean(L_MEA_min(find(L_MEA_min ~= Inf)));
        L_MEA_std = std(L_MEA_min(find(L_MEA_min ~= Inf)));
        
        L_min_7030{mm,1} = mea_num;
        L_min_7030{mm,2} = L_MEA_min;
        L_min_7030{mm,3} = L_MEA_mean;
        L_min_7030{mm,4} = L_MEA_std;
        
        
        
        for i = 1: num_CH
            nodeID(i,1) = str2num(ids{i});
        end
        
        %plot of MEA graphs
        %plot_graph;
        
        %save results
        cd('/Volumes/LaCie/NanoCom_github/Experiments/analysis/') %adapt to your path here
        
        save_name = ['graphAnalysis_' num2str(mea_num) ''];
        save(save_name, 'SR','conn_matrix_7030','graph_7030', 'degree_7030', 'edges_7030', 'L_min_7030')
        
        
    end
end
close all;

