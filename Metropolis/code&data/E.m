function err = E(theta, data)
%E Summary of this function goes here
%   Detailed explanation goes here
    %[i1, j1] = size(data);
    global ReprData1718 stdev
    if nargin == 1
        data = ReprData1718(:,2:end);
    end

    checkpoints = M(theta);
    
    err = (1/2)*sum(sum((data(:,:)-checkpoints(:,:)).^2))/(stdev^2);
end

