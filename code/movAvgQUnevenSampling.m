% movAvgQUnevenSampling Averages and interpolates unregularly spaced orientation quaternions
%    [newTime, qNew] = movAvgQUnevenSampling(time, qOld, newTime, windowSize)
%    resamples and smoothes the orientation quaternions given in the N-by-4
%    matrix qOld, sampled at time instants time. The quaternions are
%    averaged over windows of size windowSize, centered at the time
%    instants given in the M-by-1 vector newTime. 
%    The function supports a different window size for each time sample of
%    newTime. In this case newTime is a M-by-1 vector of window sizes. The
%    window size is given in number of samples.
%    Returned is a copy of newTime and qNew, the M-by-4 quaternion matrix
%    of the M averaged and interpolated orientation quaternions.

function [newTime, qNew] = movAvgQUnevenSampling(time, qOld, newTime, windowSize)
    qNew = zeros(size(newTime,1), 4);

    % Always same windowSize
    if numel(windowSize)==1
        halfWindowSize = floor(windowSize/2);
    
        for i=1:numel(newTime)
            [~, startTimeInd] = min(abs(time-(newTime(i)-halfWindowSize)));
            [~, stopTimeInd] = min(abs(time-(newTime(i)+halfWindowSize)));

            qNew(i,:) = quat_mean(qOld(startTimeInd:stopTimeInd,:));
        end

    % Different window size for each time sample
    elseif numel(windowSize)==numel(newTime)
        for i=1:numel(newTime)
            halfWindowSize = floor(windowSize(i)/2);
            [~, startTimeInd] = min(abs(time-(newTime(i)-halfWindowSize)));
            [~, stopTimeInd] = min(abs(time-(newTime(i)+halfWindowSize)));

            qNew(i,:) = quat_mean(qOld(startTimeInd:stopTimeInd,:));
        end
    else
        error('movAvgQUnevenSampling:timeWindows:wrong', 'Wrong time window vector length. WindowSize must be either 1 or match exactly the new timestamps');
    end
end