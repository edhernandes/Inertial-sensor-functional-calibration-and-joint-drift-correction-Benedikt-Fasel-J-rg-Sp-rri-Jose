% determineAzimuthOffset Determines the azimuth offset between two acceleration streams
%    azimutOffsetQuat = determineAzimutOffset(proximalAcc, proximalGyr, proximalOrientation, distalAcc, distalGyr, distalOrientation, gravityValue, fs)
%    determines the azimuth offset such that the acceleration from the
%    distal sensor (distalAcc) would be aligned with the acceleration from
%    the proximal sensor (proximalAcc).
%    Returned is the quaterion representing this rotation.

function azimutOffsetQuat = determineAzimutOffset(proximalAcc, proximalGyr, proximalOrientation, distalAcc, distalGyr, distalOrientation, gravityValue, fs)

    % Convert units to m/s2 and store in structure so that we can use the
    % convertFrames function to convert the data into the global frame
    proximalLF.acc = proximalAcc.*gravityValue;
    proximalLF.gyr = proximalGyr;
    distalLF.acc = distalAcc.*gravityValue;
    distalLF.gyr = distalGyr;
    
    % Convert to global frame
    proximalGF = convertFrames(proximalLF, proximalOrientation);
    distalGF = convertFrames(distalLF, distalOrientation);
    
    % Set y-axis to 0 since we only want azimuth offset
    distalGF.acc(:,2) = 0;
    proximalGF.acc(:,2) = 0;
    
    
    %% Determine total, 3D drift
    
    % Estimate drift only where we have a reliable estimate
    proximalNorm = sqrt(sum(proximalGF.acc.^2,2));
    distalNorm = sqrt(sum(distalGF.acc.^2,2));
    
    % There must be some acceleration
    normTooSmall = proximalNorm<0.6 | distalNorm<0.6;
    
    % Difference must be at most 35%
    normDifferenceTooLarge = abs((proximalNorm-distalNorm)./(0.5*(proximalNorm + distalNorm)))>0.25;
    
    % Remove instants with low angular velocity
    [b,a] = butter(2, 1.2465*1/fs*2, 'low');
    proximalGyrNorm = filtfilt(b,a,sqrt(sum(proximalGyr.^2,2)));
    distalGyrNorm = filtfilt(b,a,sqrt(sum(distalGyr.^2,2)));
    gyroMotionless = proximalGyrNorm<10 | distalGyrNorm<10;
    
    exclude = normDifferenceTooLarge | normTooSmall | gyroMotionless;
    
    
    % Only continue we we have at least some samples with drift correction
    if sum(double(~exclude))<50
        azimutOffsetQuat = [1 0 0 0];
        warning('Not enough reliable samples available for azimut offset estimation.');
        
    else
        totalDriftQ = zeros(size(proximalGF.acc,1),4);
        totalDriftQ(1,:) = 1;
        for i=1:size(totalDriftQ,1)
            if ~exclude(i)
                [totalDriftAxis, totalDriftAngle] = vec2helic(distalGF.acc(i,:), proximalGF.acc(i,:));
                totalDriftQ(i,:) = helic2quat(totalDriftAxis, totalDriftAngle);
            end
        end
        
        % Control plot. Uncomment to show
        % figure; plot(totalDriftQ(~exclude,:)); title('Azimut drift estimation quaternion');
        
        % Compute average rotation offset quaternion
        azimutOffsetQuat = quat_mean(totalDriftQ(~exclude,:));
    end
end