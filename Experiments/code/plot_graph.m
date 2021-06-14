%plot graphs

colormap autumn
    nSizes = 2*sqrt(degree_MEA-min(degree_MEA)+0.2);
    nColors = degree_MEA;
    f1 = figure(1);
    plot(graph_M_th,'XData',NodePosition(1,:),'YData',NodePosition(2,:),'MarkerSize',nSizes,'NodeCData',nColors,'EdgeAlpha',0.2,'NodeLabel',nodeID,'LineWidth',2)
    colorbar
    title_MEA_col = ['Graph-coloured-degree-map-MEA-' num2str(mea_num) ''];
    title(title_MEA_col)
    title_MEA_fig_col = [title_MEA_col, '.fig'];
    savefig(f1, title_MEA_fig_col)
    
    f2 = figure(2);
    colormap autumn
    p = plot(graph_M_th,'XData',NodePosition(1,:),'YData',NodePosition(2,:),'MarkerSize',nSizes, 'NodeLabel',nodeID);
    wcc = centrality(graph_M_th,'degree','Importance',graph_M_th.Edges.Weight);
    p.NodeCData = wcc;
    colorbar
    title_MEA_col = ['Graph-coloured-degree-importance-weight-map-MEA-' num2str(mea_num) ''];
    title(title_MEA_col)
    title_MEA_fig_col = [title_MEA_col, '.fig'];
    savefig(f2, title_MEA_fig_col)
    
    f3 = figure(3);
    colormap autumn
    p = plot(graph_M_th,'XData',NodePosition(1,:),'YData',NodePosition(2,:),'MarkerSize',nSizes, 'NodeLabel',nodeID);
    bcc = centrality(graph_M_th,'betweenness','Cost',graph_M_th.Edges.Weight);
    p.NodeCData = bcc;
    colorbar
    title_MEA_col = ['Graph-coloured-weight-betweenness-map-MEA' num2str(mea_num) ''];
    title(title_MEA_col)
    title_MEA_fig_col = [title_MEA_col, '.fig'];
    savefig(f3, title_MEA_fig_col)
    
    SR_mat = cell2mat(SR);
    
    f4 = figure(4);
    colormap autumn
    p = plot(graph_M_th,'XData',NodePosition(1,:),'YData',NodePosition(2,:),'MarkerSize',nSizes, 'NodeCData', SR_mat, 'NodeLabel',nodeID);
    colorbar
    title_MEA_col = ['Graph-coloured-SR-map-MEA-' num2str(mea_num) ''];
    title(title_MEA_col)
    title_MEA_fig_col = [title_MEA_col, '.fig'];
    savefig(f4, title_MEA_fig_col)
    
    f5 = figure(5);
    colormap autumn
    p = plot(graph_M_th,'XData',NodePosition(1,:),'YData',NodePosition(2,:),'MarkerSize',nSizes, 'NodeLabel',nodeID);
    wcc = centrality(graph_M_th,'closeness','Cost',graph_M_th.Edges.Weight);
    p.NodeCData = wcc;
    colorbar
    title_MEA_col = ['Graph-coloured-centrality-weight-weight-map-MEA-' num2str(mea_num) ''];
    title(title_MEA_col)
    title_MEA_fig_col = [title_MEA_col, '.fig'];
    savefig(f5, title_MEA_fig_col)
    