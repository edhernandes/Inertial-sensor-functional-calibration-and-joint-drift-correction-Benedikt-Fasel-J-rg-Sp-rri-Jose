% armCalibration Calibrate the arm IMUs
%     rotMatrix = armCalibration(imuData, calibInfo, sensorName, gravityVector, fs)
%     finds the rotation matrix that orients the arm sensor to the 
%     arm's anatomical axes. 
%      
%     Assumptions: Arm swing takes place in the saggital plane where the main
%                  rotation happens around the medio-lateral axis. 
%     Algorithm: First the medio-lateral axis is estimated and second the vertical
%                axis using upright standing

function armRotMatrix = armCalibration(imuData, calibInfo, sensorName, gravityVector, fs)

           
    % Lowpass filter parameters
    [bNormal,aNormal] = butter(2, 1.2465*5/fs*2, 'low');  % cutoff freq. 5 Hz
    [bLow,aLow] = butter(2, 1.2465*1/fs*2, 'low');        % cutoff freq. 1 Hz 
    
    
    %% Find azimut alignment rotation matrix
    
    % Hypothesis:
    % Principal axis of rotation during the armswing corresponds to the
    % medio-lateral Z-axis [0 0 1].
    
    % Get acceleration and angular velocity for the armswing
    armswing.acc = imuData.(sensorName).acc(calibInfo.armswing_start:calibInfo.armswing_stop,:);
    armswing.gyr = imuData.(sensorName).gyr(calibInfo.armswing_start:calibInfo.armswing_stop,:);
    
    % Low pass filter angular velocity to remove noise
    armswingAngVel = filtfilt3D(bNormal, aNormal, armswing.gyr);
    armswingAngVelNorm = sqrt(sum(armswingAngVel.^2,2));
    
    % To improve estimation accuracy we want to remove all instants of very
    % low rotation speeds (since this low rotation speed could by in any direction)
    % The angular velocity threshold was fixed at 30 deg/sec
    angVelThres = 30;
    
    % Compute the PCA to find the principal axis of rotation
    rotationMatrix = pca(armswingAngVel(armswingAngVelNorm>angVelThres,:));

    % Find the rotation that aligns the principal axis to the Z-axis
    [rotationAxis, rotationAngle] = vec2helic(rotationMatrix(:,1), [0 0 1]);
    armRotMatrix = quat2matrix(helic2quat(rotationAxis, rotationAngle));

    
    % Check orientation of the z axis (must point to the right)
    % Hypothesis: first peak of the angular veloctiy in flexion must be 
    %             positive (start of moving/rotating arm upwards)
    
    % Align the angular velocity according to the orientation found above
    armswingAligned = armswingAngVel * armRotMatrix';
    
    % Massively lowpass filter because we only want global oscillation pattern
    armswingLow = filtfilt(bLow,aLow,armswingAligned(:,3));
    
    % Find all the angular velocity peaks (positive and negative)
    % To be insensitive to the movement intensity: chose a variable peak
    % height (50% of the standard deviation) for the peak detection.
    [~, locs] = findpeaks(abs(armswingLow), 'minpeakheight', 0.5*std(armswingLow));
    
    % First peak negative ==> Axis is inverted (points to the left instead
    % of to the right). Rotate by 180° around the "vertical" axis
    if armswingLow(locs(1))<0
        armRotMatrix = yRotation(pi) * armRotMatrix;
    end
    
    
    %% Find the vertical axis
     
    % Get the data for the upright posture (LF = local frame)
    upright_LF.acc = imuData.(sensorName).acc(calibInfo.upright_start:calibInfo.upright_stop,:);
    upright_LF.gyr = imuData.(sensorName).gyr(calibInfo.upright_start:calibInfo.upright_stop,:);

    % Transform to the frame found above (AF = anatomical frame)
    upright_AF.acc = upright_LF.acc * armRotMatrix';
    upright_AF.gyr = upright_LF.gyr * armRotMatrix';

    % Get the average acceleration during the upright posture
    avgAcceleration = mean(upright_AF.acc);
    
    % Find the rotation difference to the true upright [0 1 0]. Project
    % onto the sagittal plane since we already aligned the medio-lateral
    % axis in step 1.
    [rotationAxis, rotationAngle] = vec2helic([avgAcceleration(1:2) 0], gravityVector);
    uprightRotation = quat2matrix(helic2quat(rotationAxis, rotationAngle));
    
    % Apply the rotation to the previously found rotation matrix.
    armRotMatrix  = uprightRotation * armRotMatrix;
end