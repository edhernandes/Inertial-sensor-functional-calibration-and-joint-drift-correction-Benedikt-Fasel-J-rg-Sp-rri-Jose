% determineFlexionOffset Determines a flexion ofset of a segment
%    flexionOffsetQuat = determineFlexionOffset(proximalLF, proximalOrientation, distalLF, distalOrientation, gravityValue, fs)
%    determines the flexion offset of one segment relative to the other by
%    comparing the recorded acceleration at their common joint centre.

function flexionOffsetQuat = determineFlexionOffset(proximalLF, proximalOrientation, distalLF, distalOrientation, gravityValue, fs)

    % Convert units to m/s^2. Leave angular velocity in deg/sec
    proximalLF.acc = proximalLF.acc .* gravityValue;
    distalLF.acc = distalLF.acc .* gravityValue;
    
    % Convert to global frame
    proximalGF = convertFrames(proximalLF, proximalOrientation);
    distalGF = convertFrames(distalLF, distalOrientation);
    
    % Set z-axis to 0 since we only want to work in sagittal plane
    distalGF.acc(:,3) = 0;
    proximalGF.acc(:,3) = 0;
    
    
    %% Determine total, 3D drift
    
    % Estimate drift only where we have a reliable estimate
    proximalNorm = sqrt(sum(proximalGF.acc.^2,2));
    distalNorm = sqrt(sum(distalGF.acc.^2,2));
    
    % There must be some acceleration
    normTooSmall = proximalNorm<0.6 | distalNorm<0.6;
    
    % Difference must be at most 35%
    normDifferenceTooLarge = abs((proximalNorm-distalNorm)./(0.5*(proximalNorm + distalNorm)))>0.35;
    
    % Remove instants with low angular velocity
    [b,a] = butter(2, 1.2465*1/fs*2, 'low');
    proximalGyrNorm = filtfilt(b,a,sqrt(sum(proximalLF.gyr.^2,2)));
    distalGyrNorm = filtfilt(b,a,sqrt(sum(distalLF.gyr.^2,2)));
    gyroMotionless = proximalGyrNorm<10 | distalGyrNorm<10;
    
    exclude = normDifferenceTooLarge | normTooSmall | gyroMotionless;
    
    
    % Only continue if we have at least some samples with drift correction
    if sum(double(~exclude))<50
        flexionOffsetQuat = [1 0 0 0];
        warning('functionalCalibration:thigh:flexionCorrection', 'Not enough reliable samples available for azimut joint drift correction.');
        
    else
        orientationDifferenceQ = zeros(size(proximalGF.acc,1),4);
        orientationDifferenceQ(1,:) = 1;
        for i=1:size(orientationDifferenceQ,1)
            if ~exclude(i)
                [orientationDifferenceAxis, orientationDifferenceAngle] = vec2helic(distalGF.acc(i,:), proximalGF.acc(i,:));
                orientationDifferenceQ(i,:) = helic2quat(orientationDifferenceAxis, orientationDifferenceAngle);
            end
        end
        
        % Compute average rotation offset quaternion
        flexionOffsetQuat = quat_mean(orientationDifferenceQ(~exclude,:));
    end
end