function [conn_matrix] = correlation_MEA(DataCell, chNUM, lengthST)

%author: Barbara Genocchi
%This function load the timestamp of firing of the two channels to be correlated (timest_x & timest_y), creates a binary signal of
%spiking (series_x & series y). And makes the cross correlation between all
%channels of the MEA.
%INPUT:  DataCell. For example:_ 
%      
%       chNUM = number of channels in MEA
%       lengthST = length of spike train
%OUTPUT: conn_matrix = correlation matrix between MEA channels (not normalized)
        
conn_matrix = zeros(chNUM,chNUM);


for row = 1:chNUM
    for col = 1:chNUM
        
        if row == col
            continue
        else
            
            timest_x =cell2mat(DataCell(row, 3));
            timest_y =cell2mat(DataCell(col, 3));
            
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

end