% estimateJointDrift Estimates the joint drift for one joint
%    jointDriftQInterp = estimateJointDrift(proximalAcc, proximalGyr, proximalOrientation, distalAcc, distalGyr, distalOrientation, gravityValue, fs)
%    estimates the joint drift for two sensors connected by a joint. The
%    "virtual" acceleration at the joint center is given in the variables
%    proximalAcc and distalAcc. 
%    Returned is a N-by-4 quaternion with the drift estimation of the
%    distal sensor relative to the proximal sensor for each of the N
%    samples.

function jointDriftQInterp = estimateJointDrift(proximalAcc, proximalGyr, proximalOrientation, distalAcc, distalGyr, distalOrientation, gravityValue, fs)
    
    % Convert units to m/s2 and store in structure so that we can use the
    % convertFrames function to convert the data into the global frame
    proximalLF.acc = proximalAcc.*gravityValue;
    proximalLF.gyr = proximalGyr;
    distalLF.acc = distalAcc.*gravityValue;
    distalLF.gyr = distalGyr;
    
    % Convert to global frame
    proximalGF = convertFrames(proximalLF, proximalOrientation);
    distalGF = convertFrames(distalLF, distalOrientation);
    
    
    %% Determine total, 3D drift
    
    % Estimate drift only where we have a reliable estimate
    proximalNorm = sqrt(sum(proximalGF.acc.^2,2));
    distalNorm = sqrt(sum(distalGF.acc.^2,2));
    
    % There must be some acceleration in order to have a high enough signal
    % to noise ratio. For skiing with a lot of vibrations the optimal value
    % seems to be around 6m/sec^2. For other applications, this value could
    % also be lowered.
    normTooSmall = proximalNorm<6 | distalNorm<6;
    
    % Difference in norm must be at most 20%
    normDifferenceTooLarge = abs((proximalNorm-distalNorm)./(0.5*(proximalNorm + distalNorm)))>0.2;
    
    % Remove instants with low angular velocity
    [b,a] = butter(2, 1.2465*1/fs*2, 'low');
    proximalGyrNorm = filtfilt(b,a,sqrt(sum(proximalGyr.^2,2)));
    distalGyrNorm = filtfilt(b,a,sqrt(sum(distalGyr.^2,2)));
    gyroMotionless = proximalGyrNorm<10 | distalGyrNorm<10;
    
    exclude = normDifferenceTooLarge | normTooSmall | gyroMotionless;

    % Only continue we we have at least some samples with drift correction
    if sum(double(~exclude))<50
        jointDriftQInterp = [];
        warning('Not enough reliable samples available for joint drift correction.');
    else
        
% %         % Only for visualization
% %         totalDriftAngleVisualization = zeros(size(proximalGF.acc,1),1);
% %         totalDriftAxisVisualization = zeros(size(proximalGF.acc,1),3);
        totalDriftQ = zeros(size(proximalGF.acc,1),4);
        totalDriftQ(1,:) = 1;
        for i=1:size(totalDriftQ,1)
            if ~exclude(i)
                [totalDriftAxis, totalDriftAngle] = vec2helic(distalGF.acc(i,:), proximalGF.acc(i,:));
                totalDriftQ(i,:) = helic2quat(totalDriftAxis, totalDriftAngle);

% %                 % Only for visualization
% %                 [totalDriftAxisVisualization(i,:), totalDriftAngleVisualization(i)] = vec2helic(distalGF.acc(i,:), proximalGF.acc(i,:));
% %                 if totalDriftAxisVisualization(i,2)<0
% %                     totalDriftAngleVisualization(i) = -totalDriftAngleVisualization(i);
% %                 end
            end
        end

        time = 1:numel(normDifferenceTooLarge);
        oldTime = time(~exclude);
        newTimeTimestep = 0.04;
        newTime = (1:newTimeTimestep*fs:size(proximalLF.acc,1)-1)';   

        % Different window size for each interpolated datapoint. Window size
        % depends on signal period: take one period to avoid correcting wrongly
        % for soft tissue artefact. Here we only want to estimate overall
        % sensor drift.

        % Get principal axis of movement
        [b,a] = butter(2, 1.2465*5/fs*2, 'low');
        [~, angVelPca] = pca(filtfilt3D(b,a,proximalLF.gyr));

        % Get approximate signal periodicity -> get all ang vel peaks
        [b,a] = butter(2, 1.2465*0.5/fs*2, 'low');
        angVelLow = filtfilt(b,a,angVelPca(:,1));
        minPeakHeight = 0.5.*std(angVelLow);
        [peakVals, locs] = findpeaks(angVelLow, 'minPeakHeight', minPeakHeight);

        % Control plot
        % figure; plot(angVelLow); hold on; plot(locs, peakVals, '.', 'markerSize', 16'); title('Movement cycle detection');

        % Get period
        cadenceOrig = diff(locs);
        cadenceVals = [cadenceOrig(1); cadenceOrig; cadenceOrig(end)];
        avgCadence = median(cadenceVals);
        
        % Remove outliers
        stdCadence = std(cadenceVals);
        cadenceVals(cadenceVals>(avgCadence+2*stdCadence) | cadenceVals<(avgCadence-2*stdCadence)) = avgCadence;
        
        cadenceTime = [1; locs(1:end-1)+0.5*cadenceOrig; size(angVelLow,1)];

        % Compute windowsize --> window size should approximately match period
        windowSize = interp1(cadenceTime, cadenceVals, newTime, 'linear', mean(cadenceVals));

        % We only want to estimate the drift if enough samples are available.
        % Want at least 25% of valid estimates
        validSamples = movAvg(double(exclude), round(3*avgCadence))<0.75;
        validSamples(1:round(3*avgCadence)) = 0;
        validSamples(end-round(3*avgCadence):end) = 0;

        % Downsample to match newTime and windowSize vectors
        validSamples = logical(interp1(time, double(validSamples), newTime, 'nearest', 0));

        newTimeSub = newTime(validSamples);
        windowSizeSub = windowSize(validSamples);

        % Control plot
        % figure; plot(cadenceTime, cadenceVals, '.-'); hold on; plot(newTime, windowSize); title('Estimated cycle duration');

        % Estimate only during valid estimate
        % windows where at least 50% of samples are available. Return empty
        % variable when no drift correction available
        if isempty(newTimeSub) || numel(newTimeSub)<10
            jointDriftQInterp = [];
            display('Warning: No samples available for joint drift correction.');
        else

            % Get drift estimation during valid samples
            [~, driftEstimated] = movAvgQUnevenSampling(oldTime, totalDriftQ(~exclude,:), newTimeSub, 2.*windowSizeSub);

%             % Control plot
%             figure; hold on; plot(oldTime, totalDriftQ(~exclude,:)); 
%             plot(newTimeSub, driftEstimated); 

            % Set beginning and end to "default" values and interpolate
            % linearly
            driftEstimated = [driftEstimated(1,:); driftEstimated; driftEstimated(end,:)];
            newTimeSub = [0; newTimeSub; size(distalOrientation,1)];

            jointDriftQ = interp1(newTimeSub, driftEstimated, newTime, 'linear');

            % Smooth. Design a very simple filter
            b = 1; a = [1 -0.975];
            scale = 0.025;  % h = freqz(b,a, 1000, 1/newTimeTimestep);  scale = 1/max(abs(h));
            b = b*scale;
            for i=1:4
                jointDriftQ(:,i) = filtfilt(b,a,jointDriftQ(:,i));
            end

%             plot(newTime, jointDriftQ);


            % Interpolate / upsample to 500Hz
            jointDriftQInterp(:,1) = interp1(newTime, jointDriftQ(:,1), time, 'linear', 1);
            jointDriftQInterp(:,2) = interp1(newTime, jointDriftQ(:,2), time, 'linear', 0);
            jointDriftQInterp(:,3) = interp1(newTime, jointDriftQ(:,3), time, 'linear', 0);
            jointDriftQInterp(:,4) = interp1(newTime, jointDriftQ(:,4), time, 'linear', 0);

            % Can happen that we don't have any info at end: replace by last known
            % drift value
            lastKnownDrift = find(jointDriftQInterp(:,1)~=1, 1, 'last');
            nbMissingSamples = size(jointDriftQInterp,1)-lastKnownDrift;
            if nbMissingSamples>0
                jointDriftQInterp(lastKnownDrift+1:end,:) = repmat(jointDriftQInterp(lastKnownDrift,:), nbMissingSamples, 1);
            end

            % Normalize to make again unit quaternions representing rotations
            for i=1:numel(time)
                jointDriftQInterp(i,:) = quat_normalize(jointDriftQInterp(i,:));
            end
        end
    end
end