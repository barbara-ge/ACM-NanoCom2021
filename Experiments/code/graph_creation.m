%graph creation

IDs = [];
    for hh = 1: num_CH
        if hh == 15
            IDs(hh,1) = 15;
        else
            IDs(hh,1) = str2num(cell2mat(DataCell(hh,1)));
        end
    end
    
    %clear noisy channels
    artefacted_channels = [15, 35, 58, 24, 47];
    for u = 1 : length(artefacted_channels)
        uu = find(IDs == artefacted_channels(1,u));
        conn_matrix(uu , :)=0;
        conn_matrix(: , uu)=0;
    end
    
    max_corr = max(conn_matrix, [],'all');
    conn_matrix_norm = conn_matrix./max_corr;
    conn_matrix_norm_rounded = round(conn_matrix_norm, 2); %if too many decimal places matrix becomes not symmetric due to internal rounding
    graph_M = graph(conn_matrix_norm_rounded);
    
    %check symmetric
%     sym_mat = zeros(num_CH, num_CH);
%     for ii = 1: num_CH
%         for jj = 1: num_CH
%             if conn_matrix_norm(ii,jj) == conn_matrix_norm(jj,ii)
%                 sym_mat(ii,jj) = 1;
%                 sym_mat(jj,ii) = 1;
%             end
%         end
%     end
%     
%     heatmap(sym_mat)
    
    for k=1:num_CH
        if k == 15
            ChannelName{k} = 15;
            SR{k} = 0;
        else
            %ChannelName{k} = DataCell{channels(k),1};
            ChannelName{k} = str2num(cell2mat(DataCell(k,1)));
            SR{k} = cell2mat(DataCell(k,2))./(lengthST/60);
        end
        ids{k}=['' num2str(ChannelName{1,k}) ''];
        
        
    end
    
    page_Size = [400 400];
    
    XpageSize = page_Size(1);
    YpageSize = page_Size(2);
    
    for i = 1:length(ids)
        nodeID = str2num(ids{i});
        tens = floor(nodeID / 10);
        ones = mod(nodeID , 10);
        xposition = floor(XpageSize / 8) * tens ;
        yposition = floor(YpageSize / 8) * (9-ones); %(9-ones) for MCS setup; (ones) for Axion setup
        
        NodePosition(1, i) = xposition;
        NodePosition(2, i) = yposition;
        
    end
    
    
    conn_matrix_th = zeros(num_CH,num_CH);
    th = 0.65;
    
    for i = 1:num_CH
        for j = 1:num_CH
            if conn_matrix_norm_rounded(i, j) > th
                conn_matrix_th(i,j) = conn_matrix_norm_rounded(i, j);
            end
        end
    end
    
    graph_M_th = graph(conn_matrix_th);