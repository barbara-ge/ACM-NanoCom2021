%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graph theory analysis and plot for the simulations results from INEXA
% (see publication with doi: "add doi here")
% author: Barbara Genocchi 
% date: 06/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%NS
NumberOfNeurons = 250;

num_CH = 64;
lengthST = 300;

fprintf('NS')

for ii = 1:10
    
    ii
    cd(['/Volumes/LaCie/NanoCom_github/Simulations/data/NS/NS_results_' num2str(ii) '/'])
    
    
    NeuronInfo = readtable('NeuronNetworkTopology.csv');
    NeuronInfo = NeuronInfo(1:end,1:end);
    NeuronInfo = table2array(NeuronInfo);
    
    false_MEA = zeros(8,8);
    pos_MEA_x = [];
    pos_MEA_y = [];
    
    for aa = 0: 93.75: 656.25
        for bb = 0: 93.75: 656.25
            
            Lx = find(NeuronInfo(1,:) < (aa + 93.75) & NeuronInfo(1,:) > aa);
            Ly = find(NeuronInfo(2,:) < (bb + 93.75) & NeuronInfo(2,:) > bb);
            
            Loc = intersect(Lx, Ly);
            
            
            
            if isempty(Loc)==1
                continue
            else
                
                pos = randi(length(Loc));
                sel_neu = Loc(pos);
                
                grid_x = (aa + 93.75)/93.75;
                grid_y = (bb + 93.75)/93.75;
                false_MEA(grid_x, grid_y) = sel_neu;
                
            end
        end
    end
    load(['/Volumes/LaCie/NanoCom_github/Simulations/data/NS/NS_results_' num2str(ii) '/burst/BurstInfoN_0.2.mat'])
    
    SR_MEA = zeros(8,8);
    BR_MEA = zeros(8,8);
    ns = [];
    
    SR_graph = [];
    for aa = 1:8
        for bb = 1:8
            neu = false_MEA(aa, bb);
            if neu == 0
                SR_MEA(aa,bb) = 0;
                BR_MEA(aa,bb) = 0;
                SR_graph = [SR_graph; 0];
                pos_MEA_x = [pos_MEA_x; 0];
                pos_MEA_y = [pos_MEA_y; 0
                    ];
            else
                SR_MEA(aa,bb) = cell2mat(SR(neu));
                BR_MEA(aa,bb) = cell2mat(BR(neu));
                SR_graph = [SR_graph; cell2mat(SR(neu))];
                ns = [ns; neu];
                pos_MEA_x = [pos_MEA_x; ((aa*93.75)-46.87)];
                pos_MEA_y = [pos_MEA_y; ((bb*93.75)-46.87)];
            end
        end
    end
    
    SR_MEA_NS_run{ii,1} = ii;
    SR_MEA_NS_run{ii,2} = SR_MEA;
    
    BR_MEA_NS_run{ii,1} = ii;
    BR_MEA_NS_run{ii,2} = BR_MEA;
    
    for row = 1: length(ns)
        for col = 1:length(ns)
            if row == col
                continue
            else
                
                timest_x =cell2mat(DataCell(ns(row), 3));
                timest_y =cell2mat(DataCell(ns(col), 3));
                
                series_x = zeros(1, lengthST*1000);
                
                for i = 1: length(timest_x)
                    a = round(timest_x(i));
                    series_x(1,a) = 1;
                end
                
                series_y = zeros(1, lengthST*1000);
                
                for i = 1: length(timest_y)
                    b = round(timest_y(i));
                    series_y(1,b) = 1;
                end
                
                [r,lags] = xcorr(series_x,series_y);
                
                
                
                zeroLag_index = find(lags==0);
                
                
                maxcorr=max(r(zeroLag_index : zeroLag_index));
                
                conn_matrix(row, col) = maxcorr;
            end
        end
    end
    
    max_corr = max(conn_matrix, [],'all');
    conn_matrix_norm = conn_matrix./max_corr;
    conn_matrix_norm_rounded = round(conn_matrix_norm, 3, 'significant'); %if too many decimal places matrix becomes not symmetric due to internal rounding
    
    if issymmetric(conn_matrix_norm_rounded)==0
        conn_matrix_norm_rounded_2 = zeros(num_CH, num_CH);
        for lll = 2: num_CH
            for mmm = 1: (lll-1) 
                a=conn_matrix_norm_rounded(mmm,lll);
                conn_matrix_norm_rounded_2(mmm,lll) = a;
                conn_matrix_norm_rounded_2(lll,mmm) = a;
            end
        end
        conn_matrix_norm_rounded = conn_matrix_norm_rounded_2;
    end
        %graph creation from cleaned and normalized correlation matrix
        graph_M = graph(conn_matrix_norm_rounded);
        
        conn_matrix_NS{ii,1} = ii;
        conn_matrix_NS{ii,2} = conn_matrix_norm_rounded;
        
        graph_NS{ii,1} = ii;
        graph_NS{ii,2} = graph_M;
        
        conn_matrix_th = zeros(num_CH,num_CH);
        th = 0.65;
        
        for i = 1:num_CH
            for j = 1:num_CH
                if conn_matrix_norm_rounded(i, j) > th
                    conn_matrix_th(i,j) = conn_matrix_norm_rounded(i, j);
                end
            end
        end
        
        rows = ones(num_CH,1);
        names = [1:num_CH]';
        names = mat2cell(num2str(names), rows);
        
        graph_M_th = graph(conn_matrix_th, names);        
        graph_NS_th{ii,1} = ii;
        graph_NS_th{ii,1} = graph_M_th;
        
        %Analysis
        
        %edges
        
        edges_MEA = numedges(graph_M_th);
        
        edges_NS{ii,1} = ii;
        edges_NS{ii,2} = edges_MEA;
        
        %degree_MEAree
        degree_MEA = degree(graph_M_th);
        mean_degree_MEA = mean(degree_MEA);
        std_degree_MEA = std(degree_MEA);
        
        degree_NS{ii,1} = ii;
        degree_NS{ii,2} = degree_MEA;
        degree_NS{ii,3} = mean_degree_MEA;
        degree_NS{ii,4} = std_degree_MEA;
        
        %shortest path
        L_MEA_min = [];
        
        for dd = 1: num_CH
            for cc = 1: num_CH
                
                if dd~=cc %no self connections
                    [L_path, min_path] = shortestpath(graph_M_th,dd,cc);
                    
                    L_MEA_min = [L_MEA_min, min_path];
                end
            end
        end
        L_MEA_mean = mean(L_MEA_min(find(L_MEA_min ~= Inf)));
        L_MEA_std = std(L_MEA_min(find(L_MEA_min ~= Inf)));
        
        L_min_NS{ii,1} = ii;
        L_min_NS{ii,2} = L_MEA_min;
        L_min_NS{ii,3} = L_MEA_mean;
        L_min_NS{ii,4} = L_MEA_std;

        node_NS{ii,1} = length(degree_MEA(degree_MEA > 0));
        
                
        % degree importance
        importance = centrality(graph_M_th, 'degree', 'Importance', graph_M_th.Edges.Weight);
        importance_list_NS{ii,1} = importance;
%         fig1 = figure(1);
%         ax1 = axes;
%         plot(ax1, graph_M_th, 'XData', pos_MEA_x, 'YData', pos_MEA_y, 'NodeCData',SR_graph, 'MarkerSize',graph_M_th.degree+1, 'EdgeAlpha',0.5, 'NodeLabel',{})
%          colormap(ax1, winter(256));
%          %%Then add colorbars and get everything lined up
%         %set(ax1,'Position',[.17 .11 .685 .815]);
%         cb1 = colorbar(ax1,'Position',[.88 .11 .0675 .815]);
%         cb1.Label.String = 'Spike rate (spikes/min)';
%         savefig(fig1, ['network_NS_' num2str(ii) '.fig'])
%         close(fig1)
%         
%        
%         
%         
    
    
    SR_mean = mean(SR_MEA_NS_run{ii,2}, 'all');
    BR_mean = mean(BR_MEA_NS_run{ii,2}, 'all');
    
    SR_NS(ii,1) = SR_mean;
    BR_NS(ii,1) = BR_mean;
    
   
end

%9010


NumberOfAstrocytes = 28;
fprintf('9010')

for ii = 1:10
    ii
    cd(['/Volumes/LaCie/NanoCom_github/Simulations/data/9010/9010_results_' num2str(ii) '/'])
    
       
    NeuronInfo = readtable('NeuronNetworkTopology.csv');
    NeuronInfo = NeuronInfo(1:end,1:end);
    NeuronInfo = table2array(NeuronInfo);
    
    
    
    false_MEA = zeros(8,8);
    pos_MEA_x = [];
    pos_MEA_y = [];
    
    
    for aa = 0: 93.75: 656.25
        for bb = 0: 93.75: 656.25
            
            Lx = find(NeuronInfo(1,:) < (aa + 93.75) & NeuronInfo(1,:) > aa);
            Ly = find(NeuronInfo(2,:) < (bb + 93.75) & NeuronInfo(2,:) > bb);
            Loc = intersect(Lx, Ly);
            
            
            
            if isempty(Loc)==1
                continue
            else
                
                pos = randi(length(Loc));
                sel_neu = Loc(pos);
                
                grid_x = (aa + 93.75)/93.75;
                grid_y = (bb + 93.75)/93.75;
                false_MEA(grid_x, grid_y) = sel_neu;
                
            end
        end
    end
    load(['/Volumes/LaCie/NanoCom_github/Simulations/data/9010/9010_results_' num2str(ii) '/burst/BurstInfoN_0.2.mat'])
    
    SR_MEA = zeros(8,8);
    BR_MEA = zeros(8,8);
    ns = [];
    
    SR_graph = [];
    for aa = 1:8
        for bb = 1:8
            neu = false_MEA(aa, bb);
            if neu == 0
                SR_MEA(aa,bb) = 0;
                BR_MEA(aa,bb) = 0;
                SR_graph = [SR_graph; 0];
                pos_MEA_x = [pos_MEA_x; 0];
                pos_MEA_y = [pos_MEA_y; 0];
            else
                SR_MEA(aa,bb) = cell2mat(SR(neu));
                BR_MEA(aa,bb) = cell2mat(BR(neu));
                SR_graph = [SR_graph; cell2mat(SR(neu))];
                ns = [ns; neu];
                pos_MEA_x = [pos_MEA_x; ((aa*93.75)-46.87)];
                pos_MEA_y = [pos_MEA_y; ((bb*93.75)-46.87)];
            end
        end
    end
    
    SR_MEA_9010_run{ii,1} = ii;
    SR_MEA_9010_run{ii,2} = SR_MEA;
    
    BR_MEA_9010_run{ii,1} = ii;
    BR_MEA_9010_run{ii,2} = BR_MEA;
    
    timest_x = [];
    timest_y = [];
    conn_matrix = zeros(num_CH, num_CH);
    for row = 1: length(ns)
        for col = 1:length(ns)
            if row == col
                continue
            else
                
                timest_x =cell2mat(DataCell(ns(row), 3));
                timest_y =cell2mat(DataCell(ns(col), 3));
                
                series_x = zeros(1, lengthST*1000);
                
                for i = 1: length(timest_x)
                    a = round(timest_x(i));
                    series_x(1,a) = 1;
                end
                
                series_y = zeros(1, lengthST*1000);
                
                for i = 1: length(timest_y)
                    b = round(timest_y(i));
                    series_y(1,b) = 1;
                end
                
                [r,lags] = xcorr(series_x,series_y, round(2*(lengthST-1)));
                
                
                
                zeroLag_index = find(lags==0);
                
                
                maxcorr=max(r(zeroLag_index : zeroLag_index));
                
                conn_matrix(row, col) = maxcorr;
            end
        end
    end
    
    max_corr = max(conn_matrix, [],'all');
    conn_matrix_norm_rounded = zeros(num_CH, num_CH);
    conn_matrix_norm = zeros(num_CH, num_CH);
    
    conn_matrix_norm = conn_matrix./max_corr;
    conn_matrix_norm_rounded = round(conn_matrix_norm, 3, 'significant'); %if too many decimal places matrix becomes not symmetric due to internal rounding
    
    if issymmetric(conn_matrix_norm_rounded)==0
        conn_matrix_norm_rounded_2 = zeros(num_CH, num_CH);
        for lll = 2: num_CH
            for mmm = 1: (lll-1) 
                a=conn_matrix_norm_rounded(mmm,lll);
                conn_matrix_norm_rounded_2(mmm,lll) = a;
                conn_matrix_norm_rounded_2(lll,mmm) = a;
            end
        end
        conn_matrix_norm_rounded = conn_matrix_norm_rounded_2;
    end
        %graph creation from cleaned and normalized correlation matrix
        graph_M = graph(conn_matrix_norm_rounded);
        
        conn_matrix_9010{ii,1} = ii;
        conn_matrix_9010{ii,2} = conn_matrix_norm_rounded;
        
        graph_9010{ii,1} = ii;
        graph_9010{ii,2} = graph_M;
        
        conn_matrix_th = zeros(num_CH,num_CH);
        th = 0.65;
        
        for i = 1:num_CH
            for j = 1:num_CH
                if conn_matrix_norm_rounded(i, j) > th
                    conn_matrix_th(i,j) = conn_matrix_norm_rounded(i, j);
                end
            end
        end
        
        rows = ones(num_CH,1);
        names = [1:num_CH]';
        names = mat2cell(num2str(names), rows);
        
        graph_M_th = graph(conn_matrix_th, names);        
        graph_9010_th{ii,1} = ii;
        graph_9010_th{ii,1} = graph_M_th;
        
        %Analysis
        
        %edges
        
        edges_MEA = numedges(graph_M_th);
        
        edges_9010{ii,1} = ii;
        edges_9010{ii,2} = edges_MEA;
        
        %degree_MEA
        degree_MEA = degree(graph_M_th);
        mean_degree_MEA = mean(degree_MEA);
        std_degree_MEA = std(degree_MEA);
        
        degree_9010{ii,1} = ii;
        degree_9010{ii,2} = degree_MEA;
        degree_9010{ii,3} = mean_degree_MEA;
        degree_9010{ii,4} = std_degree_MEA;
        
        %shortest path
        L_MEA_min = [];
        
        for dd = 1: num_CH
            for cc = 1: num_CH
                
                if dd~=cc %no self connections
                    [L_path, min_path] = shortestpath(graph_M_th,dd,cc);
                    
                    L_MEA_min = [L_MEA_min, min_path];
                end
            end
        end
        L_MEA_mean = mean(L_MEA_min(find(L_MEA_min ~= Inf)));
        L_MEA_std = std(L_MEA_min(find(L_MEA_min ~= Inf)));
        
        L_min_9010{ii,1} = ii;
        L_min_9010{ii,2} = L_MEA_min;
        L_min_9010{ii,3} = L_MEA_mean;
        L_min_9010{ii,4} = L_MEA_std;

        node_9010{ii,1} = length(degree_MEA(degree_MEA > 0));
        
        
         % degree importance
        importance = centrality(graph_M_th, 'degree', 'Importance', graph_M_th.Edges.Weight);
        importance_list_9010{ii,1} = importance;
        
        AstroConnections = readtable('AstrocyteConnections.csv');
        AstroConnections = AstroConnections(1:end,1:end);
        AstroConnections = table2array(AstroConnections);
        AstroData = importdata('AstroData_0.0200_0.7000_0.0000_0.0000_0.0000_0.7000_0.0000_0.0000_0.0000_0.0000_0.0000.csv');
        AstroInfo = readtable('AstrocyteNetworkTopology.csv');
        AstroInfo = AstroInfo(1:end,1:end);
        AstroInfo = table2array(AstroInfo);
        activity = AstroData(14:14+NumberOfAstrocytes, 1:end-1);
        activity_color = zeros(NumberOfAstrocytes,1);
        
        rows = ones(NumberOfAstrocytes,1);
        names = [1:NumberOfAstrocytes]';
        names = mat2cell(num2str(names), rows);
        
        Adj_mat_ast = AstroConnections;
        
        Astro_graph = graph(Adj_mat_ast, names);
        
                L_astro_min = [];
        
        for bb = 1: length(AstroConnections)
            for cc = 1: length(AstroConnections)
                
                if bb~=cc %no self connections
                
                [L_path, min_path] = shortestpath(Astro_graph,bb,cc);
                
                L_astro_min = [L_astro_min,min_path];
                end
            end
        end
        
        L_astro_min = L_astro_min(isinf(L_astro_min) == 0);
        L_astro_min = L_astro_min(isnan(L_astro_min) == 0);
        L_astro_mean = mean(L_astro_min);
        L_astro_std = std(L_astro_min);
        
        L_astro_run_9010{ii,1} = L_astro_mean;
         L_astro_run_9010{ii,2} = L_astro_std;
        %% degree from graph
        
        k_ast = degree(Astro_graph);
        k_ast_mean = mean(k_ast);
        k_ast_std = std(k_ast);
        
        astro_degree_9010{ii, 1} = k_ast_mean;
        astro_degree_9010{ii, 2} = k_ast_std;
         
        for kk = 1:NumberOfAstrocytes
            count = 0;
            for tt = 1:(length(activity)-1)
                if activity(kk, tt) == 1 && activity(kk,tt+1) == -1
                    count = count+1;
                end
            end
            activity_color(kk,1) = count;
        end
        
        activity_size = ones(NumberOfAstrocytes,1);
        
        
        for kk = 1: NumberOfAstrocytes
            if Astro_graph.degree(kk)>0
                activity_size(kk,1) = Astro_graph.degree(kk)*3;
            end
        end
        
        activity_astro_9010{ii,1} = mean(activity_color);
        
%         fig1 = figure(1);
%         ax1 = axes;
%         plot(ax1, Astro_graph, 'XData', AstroInfo(1,:), 'YData', AstroInfo(2,:), 'NodeCData',activity_color, 'Marker', 'diamond', 'MarkerSize',activity_size+1, 'EdgeAlpha',0.5, 'LineStyle', '-.','EdgeColor', 'r', 'NodeLabel',{})
%          colormap(ax1, autumn(256));
%          ax2 = axes;
%         plot(ax2, graph_M_th, 'XData', pos_MEA_x, 'YData', pos_MEA_y, 'NodeCData',SR_graph, 'MarkerSize',graph_M_th.degree+1, 'EdgeAlpha',0.5, 'EdgeColor','b', 'NodeLabel',{})
%          colormap(ax2, winter(256));
%          ax2.Visible = 'off';
%          %%Then add colorbars and get everything lined up
%         %set([ax1,ax2],'Position',[.17 .11 .685 .815]);
%         cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]);
%         cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
%         cb1.Label.String = 'Astro active time (s)';
%         cb2.Label.String = 'Spike rate (spikes/min)';
%         savefig(fig1, ['network_combined_9010_' num2str(ii) '.fig'])
%         close(fig1)
        

%         

    
    
    SR_mean = mean(SR_MEA_9010_run{ii,2}, 'all');
    BR_mean = mean(BR_MEA_9010_run{ii,2}, 'all');
    
    SR_9010(ii,1) = SR_mean;
    BR_9010(ii,1) = BR_mean;
    
    
end


NumberOfAstrocytes = 62;
fprintf('8020')

for ii = 1:10
    ii
    cd(['/Volumes/LaCie/NanoCom_github/Simulations/data/8020/8020_results_' num2str(ii) '/'])
    
       
    NeuronInfo = readtable('NeuronNetworkTopology.csv');
    NeuronInfo = NeuronInfo(1:end,1:end);
    NeuronInfo = table2array(NeuronInfo);
    
    
    
    false_MEA = zeros(8,8);
    pos_MEA_x = [];
    pos_MEA_y = [];
    
    
    for aa = 0: 93.75: 656.25
        for bb = 0: 93.75: 656.25
            
            Lx = find(NeuronInfo(1,:) < (aa + 93.75) & NeuronInfo(1,:) > aa);
            Ly = find(NeuronInfo(2,:) < (bb + 93.75) & NeuronInfo(2,:) > bb);
            Loc = intersect(Lx, Ly);
            
            
            
            if isempty(Loc)==1
                continue
            else
                
                pos = randi(length(Loc));
                sel_neu = Loc(pos);
                
                grid_x = (aa + 93.75)/93.75;
                grid_y = (bb + 93.75)/93.75;
                false_MEA(grid_x, grid_y) = sel_neu;
                
            end
        end
    end
    load(['/Volumes/LaCie/NanoCom_github/Simulations/data/8020/8020_results_' num2str(ii) '/burst/BurstInfoN_0.2.mat'])
    
    SR_MEA = zeros(8,8);
    BR_MEA = zeros(8,8);
    ns = [];
    
    SR_graph = [];
    for aa = 1:8
        for bb = 1:8
            neu = false_MEA(aa, bb);
            if neu == 0
                SR_MEA(aa,bb) = 0;
                BR_MEA(aa,bb) = 0;
                SR_graph = [SR_graph; 0];
                pos_MEA_x = [pos_MEA_x; 0];
                pos_MEA_y = [pos_MEA_y; 0];
            else
                SR_MEA(aa,bb) = cell2mat(SR(neu));
                BR_MEA(aa,bb) = cell2mat(BR(neu));
                SR_graph = [SR_graph; cell2mat(SR(neu))];
                ns = [ns; neu];
                pos_MEA_x = [pos_MEA_x; ((aa*93.75)-46.87)];
                pos_MEA_y = [pos_MEA_y; ((bb*93.75)-46.87)];
            end
        end
    end
    
    SR_MEA_8020_run{ii,1} = ii;
    SR_MEA_8020_run{ii,2} = SR_MEA;
    
    BR_MEA_8020_run{ii,1} = ii;
    BR_MEA_8020_run{ii,2} = BR_MEA;
    
    timest_x = [];
    timest_y = [];
    conn_matrix = zeros(num_CH, num_CH);
    for row = 1: length(ns)
        for col = 1:length(ns)
            if row == col
                continue
            else
                
                timest_x =cell2mat(DataCell(ns(row), 3));
                timest_y =cell2mat(DataCell(ns(col), 3));
                
                series_x = zeros(1, lengthST*1000);
                
                for i = 1: length(timest_x)
                    a = round(timest_x(i));
                    series_x(1,a) = 1;
                end
                
                series_y = zeros(1, lengthST*1000);
                
                for i = 1: length(timest_y)
                    b = round(timest_y(i));
                    series_y(1,b) = 1;
                end
                
                [r,lags] = xcorr(series_x,series_y, round(2*(lengthST-1)));
                
                
                
                zeroLag_index = find(lags==0);
                
                
                maxcorr=max(r(zeroLag_index : zeroLag_index));
                
                conn_matrix(row, col) = maxcorr;
            end
        end
    end
    
    max_corr = max(conn_matrix, [],'all');
    conn_matrix_norm_rounded = zeros(num_CH, num_CH);
    conn_matrix_norm = zeros(num_CH, num_CH);
    
    conn_matrix_norm = conn_matrix./max_corr;
    conn_matrix_norm_rounded = round(conn_matrix_norm, 3, 'significant'); %if too many decimal places matrix becomes not symmetric due to internal rounding
    
    if issymmetric(conn_matrix_norm_rounded)==0
        conn_matrix_norm_rounded_2 = zeros(num_CH, num_CH);
        for lll = 2: num_CH
            for mmm = 1: (lll-1) 
                a=conn_matrix_norm_rounded(mmm,lll);
                conn_matrix_norm_rounded_2(mmm,lll) = a;
                conn_matrix_norm_rounded_2(lll,mmm) = a;
            end
        end
        conn_matrix_norm_rounded = conn_matrix_norm_rounded_2;
    end
        %graph creation from cleaned and normalized correlation matrix
        graph_M = graph(conn_matrix_norm_rounded);
        
        conn_matrix_8020{ii,1} = ii;
        conn_matrix_8020{ii,2} = conn_matrix_norm_rounded;
        
        graph_8020{ii,1} = ii;
        graph_8020{ii,2} = graph_M;
        
        conn_matrix_th = zeros(num_CH,num_CH);
        th = 0.65;
        
        for i = 1:num_CH
            for j = 1:num_CH
                if conn_matrix_norm_rounded(i, j) > th
                    conn_matrix_th(i,j) = conn_matrix_norm_rounded(i, j);
                end
            end
        end
        
        rows = ones(num_CH,1);
        names = [1:num_CH]';
        names = mat2cell(num2str(names), rows);
        
        graph_M_th = graph(conn_matrix_th, names);    
        
        graph_8020_th{ii,1} = ii;
        graph_8020_th{ii,1} = graph_M_th;
        
        %Analysis
        
        %edges
        
        edges_MEA = numedges(graph_M_th);
        
        edges_8020{ii,1} = ii;
        edges_8020{ii,2} = edges_MEA;
        
        %degree_MEAree
        degree_MEA = degree(graph_M_th);
        mean_degree_MEA = mean(degree_MEA);
        std_degree_MEA = std(degree_MEA);
        
        degree_8020{ii,1} = ii;
        degree_8020{ii,2} = degree_MEA;
        degree_8020{ii,3} = mean_degree_MEA;
        degree_8020{ii,4} = std_degree_MEA;
        
        %shortest path
        L_MEA_min = [];
        
        for dd = 1: num_CH
            for cc = 1: num_CH
                
                if dd~=cc %no self connections
                    [L_path, min_path] = shortestpath(graph_M_th,dd,cc);
                    
                    L_MEA_min = [L_MEA_min, min_path];
                end
            end
        end
        L_MEA_mean = mean(L_MEA_min(find(L_MEA_min ~= Inf)));
        L_MEA_std = std(L_MEA_min(find(L_MEA_min ~= Inf)));
        
        L_min_8020{ii,1} = ii;
        L_min_8020{ii,2} = L_MEA_min;
        L_min_8020{ii,3} = L_MEA_mean;
        L_min_8020{ii,4} = L_MEA_std;

        node_8020{ii,1} = length(degree_MEA(degree_MEA > 0));
        
        
         % degree importance
        importance = centrality(graph_M_th, 'degree', 'Importance', graph_M_th.Edges.Weight);
        importance_list_8020{ii,1} = importance;
        
        AstroConnections = readtable('AstrocyteConnections.csv');
        AstroConnections = AstroConnections(1:end,1:end);
        AstroConnections = table2array(AstroConnections);
        AstroData = importdata('AstroData_0.0200_0.7000_0.0000_0.0000_0.0000_0.7000_0.0000_0.0000_0.0000_0.0000_0.0000.csv');
        AstroInfo = readtable('AstrocyteNetworkTopology.csv');
        AstroInfo = AstroInfo(1:end,1:end);
        AstroInfo = table2array(AstroInfo);
        activity = AstroData(14:14+NumberOfAstrocytes, 1:end-1);
        activity_color = zeros(NumberOfAstrocytes,1);
        
        rows = ones(NumberOfAstrocytes,1);
        names = [1:NumberOfAstrocytes]';
        names = mat2cell(num2str(names), rows);
        
        Adj_mat_ast = AstroConnections;
        
        Astro_graph = graph(Adj_mat_ast, names);
        
                L_astro_min = [];
        
        for bb = 1: length(AstroConnections)
            for cc = 1: length(AstroConnections)
                
                if bb~=cc %no self connections
                
                [L_path, min_path] = shortestpath(Astro_graph,bb,cc);
                
                L_astro_min = [L_astro_min,min_path];
                end
            end
        end
        
        L_astro_min = L_astro_min(isinf(L_astro_min) == 0);
        L_astro_min = L_astro_min(isnan(L_astro_min) == 0);
        L_astro_mean = mean(L_astro_min);
        L_astro_std = std(L_astro_min);
        
        L_astro_run{ii,1} = L_astro_mean;
         L_astro_run{ii,2} = L_astro_std;
        %% degree from graph
        
        k_ast = degree(Astro_graph);
        k_ast_mean = mean(k_ast);
        k_ast_std = std(k_ast);
        
        astro_degree_8020{ii, 1} = k_ast_mean;
         astro_degree_8020{ii, 2} = k_ast_std;
         
        for kk = 1:NumberOfAstrocytes
            count = 0;
            for tt = 1:(length(activity)-1)
                if activity(kk, tt) == 1 && activity(kk,tt+1) == -1
                    count = count+1;
                end
            end
            activity_color(kk,1) = count;
        end
        
        activity_size = ones(NumberOfAstrocytes,1);
        
        
        for kk = 1: NumberOfAstrocytes
            if Astro_graph.degree(kk)>0
                activity_size(kk,1) = Astro_graph.degree(kk)*3;
            end
        end
        
        activity_astro_8020{ii,1} = mean(activity_color);
        
%         fig1 = figure(1);
%         ax1 = axes;
%         plot(ax1, Astro_graph, 'XData', AstroInfo(1,:), 'YData', AstroInfo(2,:), 'NodeCData',activity_color, 'Marker', 'diamond', 'MarkerSize',activity_size+1, 'EdgeAlpha',0.5, 'LineStyle', '-.','EdgeColor', 'r', 'NodeLabel',{})
%          colormap(ax1, autumn(256));
%          ax2 = axes;
%         plot(ax2, graph_M_th, 'XData', pos_MEA_x, 'YData', pos_MEA_y, 'NodeCData',SR_graph, 'MarkerSize',graph_M_th.degree+1, 'EdgeAlpha',0.5, 'EdgeColor','b', 'NodeLabel',{})
%          colormap(ax2, winter(256));
%          ax2.Visible = 'off';
%          %%Then add colorbars and get everything lined up
%         %set([ax1,ax2],'Position',[.17 .11 .685 .815]);
%         cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]);
%         cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
%         cb1.Label.String = 'Astro active time (s)';
%         cb2.Label.String = 'Spike rate (spikes/min)';
%         savefig(fig1, ['network_combined_8020_' num2str(ii) '.fig'])
%         close(fig1)
%         

    
    SR_mean = mean(SR_MEA_8020_run{ii,2}, 'all');
    BR_mean = mean(BR_MEA_8020_run{ii,2}, 'all');
    
    SR_8020(ii,1) = SR_mean;
    BR_8020(ii,1) = BR_mean;
    
    
end


NumberOfAstrocytes = 107;
fprintf('7030')

for ii = 1:10
    ii
   cd(['/Volumes/LaCie/NanoCom_github/Simulations/data/7030/7030_results_' num2str(ii) '/'])
    
       
    NeuronInfo = readtable('NeuronNetworkTopology.csv');
    NeuronInfo = NeuronInfo(1:end,1:end);
    NeuronInfo = table2array(NeuronInfo);
    
    
    
    false_MEA = zeros(8,8);
    pos_MEA_x = [];
    pos_MEA_y = [];
    
    
    for aa = 0: 93.75: 656.25
        for bb = 0: 93.75: 656.25
            
            Lx = find(NeuronInfo(1,:) < (aa + 93.75) & NeuronInfo(1,:) > aa);
            Ly = find(NeuronInfo(2,:) < (bb + 93.75) & NeuronInfo(2,:) > bb);
            Loc = intersect(Lx, Ly);
            
            
            
            if isempty(Loc)==1
                continue
            else
                
                pos = randi(length(Loc));
                sel_neu = Loc(pos);
                
                grid_x = (aa + 93.75)/93.75;
                grid_y = (bb + 93.75)/93.75;
                false_MEA(grid_x, grid_y) = sel_neu;
                
            end
        end
    end
    load(['/Volumes/LaCie/NanoCom_github/Simulations/data/7030/7030_results_' num2str(ii) '/burst/BurstInfoN_0.2.mat'])
    
    SR_MEA = zeros(8,8);
    BR_MEA = zeros(8,8);
    ns = [];
    
    SR_graph = [];
    for aa = 1:8
        for bb = 1:8
            neu = false_MEA(aa, bb);
            if neu == 0
                SR_MEA(aa,bb) = 0;
                BR_MEA(aa,bb) = 0;
                SR_graph = [SR_graph; 0];
                pos_MEA_x = [pos_MEA_x; 0];
                pos_MEA_y = [pos_MEA_y; 0];
            else
                SR_MEA(aa,bb) = cell2mat(SR(neu));
                BR_MEA(aa,bb) = cell2mat(BR(neu));
                SR_graph = [SR_graph; cell2mat(SR(neu))];
                ns = [ns; neu];
                pos_MEA_x = [pos_MEA_x; ((aa*93.75)-46.87)];
                pos_MEA_y = [pos_MEA_y; ((bb*93.75)-46.87)];
            end
        end
    end
    
    SR_MEA_7030_run{ii,1} = ii;
    SR_MEA_7030_run{ii,2} = SR_MEA;
    
    BR_MEA_7030_run{ii,1} = ii;
    BR_MEA_7030_run{ii,2} = BR_MEA;
    
    timest_x = [];
    timest_y = [];
    conn_matrix = zeros(num_CH, num_CH);
    for row = 1: length(ns)
        for col = 1:length(ns)
            if row == col
                continue
            else
                
                timest_x =cell2mat(DataCell(ns(row), 3));
                timest_y =cell2mat(DataCell(ns(col), 3));
                
                series_x = zeros(1, lengthST*1000);
                
                for i = 1: length(timest_x)
                    a = round(timest_x(i));
                    series_x(1,a) = 1;
                end
                
                series_y = zeros(1, lengthST*1000);
                
                for i = 1: length(timest_y)
                    b = round(timest_y(i));
                    series_y(1,b) = 1;
                end
                
                [r,lags] = xcorr(series_x,series_y, round(2*(lengthST-1)));
                
                
                
                zeroLag_index = find(lags==0);
                
                
                maxcorr=max(r(zeroLag_index : zeroLag_index));
                
                conn_matrix(row, col) = maxcorr;
            end
        end
    end
    
    max_corr = max(conn_matrix, [],'all');
    conn_matrix_norm_rounded = zeros(num_CH, num_CH);
    conn_matrix_norm = zeros(num_CH, num_CH);
    
    conn_matrix_norm = conn_matrix./max_corr;
    conn_matrix_norm_rounded = round(conn_matrix_norm, 3, 'significant'); %if too many decimal places matrix becomes not symmetric due to internal rounding
    
    if issymmetric(conn_matrix_norm_rounded)==0
        conn_matrix_norm_rounded_2 = zeros(num_CH, num_CH);
        for lll = 2: num_CH
            for mmm = 1: (lll-1) 
                a=conn_matrix_norm_rounded(mmm,lll);
                conn_matrix_norm_rounded_2(mmm,lll) = a;
                conn_matrix_norm_rounded_2(lll,mmm) = a;
            end
        end
        conn_matrix_norm_rounded = conn_matrix_norm_rounded_2;
    end
        %graph creation from cleaned and normalized correlation matrix
        graph_M = graph(conn_matrix_norm_rounded);
        
        conn_matrix_7030{ii,1} = ii;
        conn_matrix_7030{ii,2} = conn_matrix_norm_rounded;
        
        graph_7030{ii,1} = ii;
        graph_7030{ii,2} = graph_M;
        
        conn_matrix_th = zeros(num_CH,num_CH);
        th = 0.65;
        
        for i = 1:num_CH
            for j = 1:num_CH
                if conn_matrix_norm_rounded(i, j) > th
                    conn_matrix_th(i,j) = conn_matrix_norm_rounded(i, j);
                end
            end
        end
        
        rows = ones(num_CH,1);
        names = [1:num_CH]';
        names = mat2cell(num2str(names), rows);
        
        graph_M_th = graph(conn_matrix_th, names);        
        graph_7030_th{ii,1} = ii;
        graph_7030_th{ii,1} = graph_M_th;
        
        %Analysis
        
        %edges
        
        edges_MEA = numedges(graph_M_th);
        
        edges_7030{ii,1} = ii;
        edges_7030{ii,2} = edges_MEA;
        
        %degree_MEAree
        degree_MEA = degree(graph_M_th);
        mean_degree_MEA = mean(degree_MEA);
        std_degree_MEA = std(degree_MEA);
        
        degree_7030{ii,1} = ii;
        degree_7030{ii,2} = degree_MEA;
        degree_7030{ii,3} = mean_degree_MEA;
        degree_7030{ii,4} = std_degree_MEA;
        
        %shortest path
        L_MEA_min = [];
        
        for dd = 1: num_CH
            for cc = 1: num_CH
                
                if dd~=cc %no self connections
                    [L_path, min_path] = shortestpath(graph_M_th,dd,cc);
                    
                    L_MEA_min = [L_MEA_min, min_path];
                end
            end
        end
        L_MEA_mean = mean(L_MEA_min(find(L_MEA_min ~= Inf)));
        L_MEA_std = std(L_MEA_min(find(L_MEA_min ~= Inf)));
        
        L_min_7030{ii,1} = ii;
        L_min_7030{ii,2} = L_MEA_min;
        L_min_7030{ii,3} = L_MEA_mean;
        L_min_7030{ii,4} = L_MEA_std;

        node_7030{ii,1} = length(degree_MEA(degree_MEA > 0));
        
        
        
         % degree importance
        importance = centrality(graph_M_th, 'degree', 'Importance', graph_M_th.Edges.Weight);
        importance_list_7030{ii,1} = importance;
        
        AstroConnections = readtable('AstrocyteConnections.csv');
        AstroConnections = AstroConnections(1:end,1:end);
        AstroConnections = table2array(AstroConnections);
        AstroData = importdata('AstroData_0.0200_0.7000_0.0000_0.0000_0.0000_0.7000_0.0000_0.0000_0.0000_0.0000_0.0000.csv');
        AstroInfo = readtable('AstrocyteNetworkTopology.csv');
        AstroInfo = AstroInfo(1:end,1:end);
        AstroInfo = table2array(AstroInfo);
        activity = AstroData(14:14+NumberOfAstrocytes, 1:end-1);
        activity_color = zeros(NumberOfAstrocytes,1);
        
        rows = ones(NumberOfAstrocytes,1);
        names = [1:NumberOfAstrocytes]';
        names = mat2cell(num2str(names), rows);
        
        Adj_mat_ast = AstroConnections;
        
        Astro_graph = graph(Adj_mat_ast, names);
        
                 L_astro_min = [];
        
        for bb = 1: length(AstroConnections)
            for cc = 1: length(AstroConnections)
                
                if bb~=cc %no self connections
                
                [L_path, min_path] = shortestpath(Astro_graph,bb,cc);
                
                L_astro_min = [L_astro_min,min_path];
                end
            end
        end
        
        L_astro_min = L_astro_min(isinf(L_astro_min) == 0);
        L_astro_min = L_astro_min(isnan(L_astro_min) == 0);
        L_astro_mean = mean(L_astro_min);
        L_astro_std = std(L_astro_min);
        
        L_astro_run{ii,1} = L_astro_mean;
         L_astro_run{ii,2} = L_astro_std;
        %% degree from graph
        
        k_ast = degree(Astro_graph);
        k_ast_mean = mean(k_ast);
        k_ast_std = std(k_ast);
        
        astro_degree_7030{ii, 1} = k_ast_mean;
         astro_degree_7030{ii, 2} = k_ast_std;
         
         
        for kk = 1:NumberOfAstrocytes
            count = 0;
            for tt = 1:(length(activity)-1)
                if activity(kk, tt) == 1 && activity(kk,tt+1) == -1
                    count = count+1;
                end
            end
            activity_color(kk,1) = count;
        end
        
        activity_size = ones(NumberOfAstrocytes,1);
        
        
        for kk = 1: NumberOfAstrocytes
            if Astro_graph.degree(kk)>0
                activity_size(kk,1) = Astro_graph.degree(kk)*3;
            end
        end
        
        activity_astro_7030{ii,1} = mean(activity_color);
        
%         fig1 = figure(1);
%         ax1 = axes;
%         plot(ax1, Astro_graph, 'XData', AstroInfo(1,:), 'YData', AstroInfo(2,:), 'NodeCData',activity_color, 'Marker', 'diamond', 'MarkerSize',activity_size+1, 'EdgeAlpha',0.5, 'LineStyle', '-.','EdgeColor', 'r', 'NodeLabel',{})
%          colormap(ax1, autumn(256));
%          ax2 = axes;
%         plot(ax2, graph_M_th, 'XData', pos_MEA_x, 'YData', pos_MEA_y, 'NodeCData',SR_graph, 'MarkerSize',graph_M_th.degree+1, 'EdgeAlpha',0.5, 'EdgeColor','b', 'NodeLabel',{})
%          colormap(ax2, winter(256));
%          ax2.Visible = 'off';
%          %%Then add colorbars and get everything lined up
%         %set([ax1,ax2],'Position',[.17 .11 .685 .815]);
%         cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]);
%         cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
%         cb1.Label.String = 'Astro active time (s)';
%         cb2.Label.String = 'Spike rate (spikes/min)';
%         savefig(fig1, ['network_combined_7030_' num2str(ii) '.fig'])
%         close(fig1)
%         

    
    SR_mean = mean(SR_MEA_7030_run{ii,2}, 'all');
    BR_mean = mean(BR_MEA_7030_run{ii,2}, 'all');
    
    SR_7030(ii,1) = SR_mean;
    BR_7030(ii,1) = BR_mean;
    
    
end
cd ../../..
save('Workspace_graph_analysis.mat')

SR_list_NS = [];
degree_list_NS = [];
degree_importance_NS = [];

SR_list_9010 = [];
degree_list_9010 = [];
degree_importance_9010 = [];

SR_list_8020 = [];
degree_list_8020 = [];
degree_importance_8020 = [];

SR_list_7030 = [];
degree_list_7030 = [];
degree_importance_7030 = [];

for ii = 1:10
    
    a = reshape(SR_MEA_NS_run{ii,2}, 1, []);
    SR_list_NS = [SR_list_NS; a'];
    degree_list_NS = [degree_list_NS; degree_NS{ii, 2}];
    degree_importance_NS = [degree_importance_NS; importance_list_NS{ii,1}];

    b = reshape(SR_MEA_9010_run{ii,2}, 1, []);
    SR_list_9010 = [SR_list_9010; b'];
    degree_list_9010 = [degree_list_9010; degree_9010{ii, 2}];
    degree_importance_9010 = [degree_importance_9010; importance_list_9010{ii,1}];
    
    c = reshape(SR_MEA_8020_run{ii,2}, 1, []);
    SR_list_8020 = [SR_list_8020; c'];
    degree_list_8020 = [degree_list_8020; degree_8020{ii, 2}];
    degree_importance_8020 = [degree_importance_8020; importance_list_8020{ii,1}];
    
    d = reshape(SR_MEA_7030_run{ii,2}, 1, []);
    SR_list_7030 = [SR_list_7030; d'];
    degree_list_7030 = [degree_list_7030; degree_7030{ii, 2}];
    degree_importance_7030 = [degree_importance_7030; importance_list_7030{ii,1}];
end

[R_NS,P_NS] = corrcoef(SR_list_NS(degree_list_NS>0), degree_list_NS(degree_list_NS>0));
[R_9010,P_9010] = corrcoef(SR_list_9010(degree_list_9010>0), degree_list_9010(degree_list_9010>0));
[R_8020,P_8020] = corrcoef(SR_list_8020(degree_list_8020>0), degree_list_8020(degree_list_8020>0));
[R_7030,P_7030] = corrcoef(SR_list_7030(degree_list_7030>0), degree_list_7030(degree_list_7030>0));


line1 = fitlm(degree_list_NS(degree_list_NS>0),SR_list_NS(degree_list_NS>0));
line2 = fitlm(degree_list_9010(degree_list_9010>0),SR_list_9010(degree_list_9010>0));
line3 = fitlm(degree_list_8020(degree_list_8020>0),SR_list_8020(degree_list_8020>0));
line4 = fitlm(degree_list_7030(degree_list_7030>0),SR_list_7030(degree_list_7030>0));

int1 = line1.Coefficients{1,1};
reg1 = line1.Coefficients{2,1};

int2 = line2.Coefficients{1,1};
reg2 = line2.Coefficients{2,1};

int3 = line3.Coefficients{1,1};
reg3 = line3.Coefficients{2,1};

int4 = line4.Coefficients{1,1};
reg4 = line4.Coefficients{2,1};

lineNS = int1 + reg1.*[0:40];
line9010 = int2 + reg2.*[0:40];
line8020 = int3 + reg3.*[0:40];
line7030 = int4 + reg4.*[0:40];

label_NS = ['NS - R = ' num2str(R_NS(1,2)) ''];
label_9010 = ['90/10 - R = ' num2str(R_9010(1,2)) ''];
label_8020 = ['80/20 - R = ' num2str(R_8020(1,2)) ''];
label_7030 = ['70/30 - R = ' num2str(R_7030(1,2)) ''];

f1 = figure(1);

bubblechart( degree_list_NS(:,1), SR_list_NS(:,1),degree_importance_NS(:,1),'b','DisplayName', label_NS)
hold on
plot(lineNS, '-b', 'LineWidth', 4)
hold on
bubblechart( degree_list_9010(:,1),SR_list_9010(:,1),degree_importance_9010(:,1),'c','DisplayName', label_9010)
hold on
plot(line9010,'-c', 'LineWidth', 4)
hold on
bubblechart( degree_list_8020(:,1),SR_list_8020(:,1),degree_importance_8020(:,1),'r','DisplayName', label_8020)
hold on
plot(line8020, '-r', 'LineWidth', 4)
hold on
bubblechart( degree_list_7030(:,1),SR_list_7030(:,1), degree_importance_7030(:,1),'y','DisplayName', label_7030)
hold on
plot(line7030, '-y', 'LineWidth', 4)
xlim([0, 30])
ylim([0, 400])
xlabel('Degree')
ylabel('Spike rate (spikes/min)')
set(gca,'FontSize',30, 'FontWeight', 'bold')
legend('Orientation', 'horizontal', 'Location', 'northoutside')


%
% [p_SR_9010] = ranksum(SR_NS, SR_9010);
% [p_SR_8020] = ranksum(SR_NS, SR_8020);
% [p_SR_7030] = ranksum(SR_NS, SR_7030);
%
% [p_BR_9010] = ranksum(BR_NS, BR_9010);
% [p_BR_8020] = ranksum(BR_NS, BR_8020);
% [p_BR_7030] = ranksum(BR_NS, BR_7030);
%
% [p_degree_9010] = ranksum(degree_NS, degree_NS_9010);
% [p_degree_8020] = ranksum(degree_NS, degree_NS_8020);
% [p_degree_7030] = ranksum(degree_NS, degree_NS_7030);
%
% [p_L_9010] = ranksum(L_NS, L_NS_9010);
% [p_L_8020] = ranksum(L_NS, L_NS_8020);
% [p_L_7030] = ranksum(L_NS, L_NS_7030);
%
