% thighCalibration Calibrate the thigh IMU
%     thighRotMatrix = thighCalibration(imuData, calibInfo, side, sacrumRotMatrix, gravityVector, gravityValue, fs)
%     computes the rotation matrix required for the functional calibration
%     of the thigh.

function thighRotMatrix = thighCalibration(imuData, calibInfo, side, sacrumRotMatrix, gravityVector, gravityValue, fs)
   

    % Lowpass filter parameters
    [bNormal,aNormal] = butter(2, 1.2465*5/fs*2, 'low');  % cutoff freq. 5 Hz
    [bLow,aLow] = butter(2, 1.2465*1/fs*2, 'low');        % cutoff freq. 1 Hz 
    
    
    % Get acceleration and angular velocity for the squats
    squatSacrum_LF.acc = imuData.sacrum.acc(calibInfo.squat_start:calibInfo.squat_stop,:);
    squatSacrum_LF.gyr = imuData.sacrum.gyr(calibInfo.squat_start:calibInfo.squat_stop,:);
    if strcmpi(side, 'left')
        squatThigh_LF.acc = imuData.leftThigh.acc(calibInfo.squat_start:calibInfo.squat_stop,:);
        squatThigh_LF.gyr = imuData.leftThigh.gyr(calibInfo.squat_start:calibInfo.squat_stop,:);
    else
        squatThigh_LF.acc = imuData.rightThigh.acc(calibInfo.squat_start:calibInfo.squat_stop,:);
        squatThigh_LF.gyr = imuData.rightThigh.gyr(calibInfo.squat_start:calibInfo.squat_stop,:);
    end
       
       
    %% Find medio-lateral alignment rotation matrix
    
    % Hypothesis:
    % Principal axis of rotation during the squats corresponds to the
    % medio-lateral Z-axis [0 0 1].
    
    % Note: For more details please also refer to the comments in the 
    % function spineCalibration
    
    % Extract angular velocity
    thighAngVel = filtfilt3D(bNormal, aNormal, squatThigh_LF.gyr);

    % Compute the norm
    thighAngVelNorm = sqrt(sum(thighAngVel.^2,2));
    
    % Only consider moments with some rotation to avoid considering also movement noise
    angVelThres = 25; 
    takeValues = thighAngVelNorm>angVelThres;

    % Ignore this calibration if no movement present
    if sum(double(takeValues))<10
        thighRotMatrix = eye(3);
        warning('functionalCalibration:thigh:noMovement', 'Squat movement not present for calibration of thigh. Ignoring...');
    else

        % Compute PCA
        thighCoeff = pca(thighAngVel(takeValues, :));

        % Get first rotation axis and align with medio-lateral (z) axis
        [thighRotationAxis, thighRotationAngle] = vec2helic(thighCoeff(:,1), [0 0 1]);

        thighRotMatrix = quat2matrix(helic2quat(thighRotationAxis, thighRotationAngle));

        % Check that z-axis is pointing to the right
        % --> want to have first peak in ang vel positive
        thighAngVelCalib = thighAngVel * thighRotMatrix';

        % Massively lowpass filter because we only want global oscillation pattern
        squatLow = filtfilt(bLow,aLow, thighAngVelCalib(:,3));

        % Find all the angular velocity peaks
        [~, locs] = findpeaks(abs(squatLow), 'minpeakheight', std(squatLow)/2);

        % Correct orientation
        if squatLow(locs(1))<0 
            thighRotMatrix = yRotation(pi) * thighRotMatrix;
        end 
    end
    
    
    %% Approximate thigh flexion 
    
    % Vertically orient thigh during motionless at beginning of squats
    lastValidInd = max([50 round(0.8*find(thighAngVelNorm>20, 1, 'first'))]);
    
    % Get acceleration in AF and compute mean during motionless
    squatThigh_AF.acc = squatThigh_LF.acc * thighRotMatrix';
    squatThigh_AF.gyr = squatThigh_LF.gyr * thighRotMatrix';
    motionlessAvgAcc = mean(squatThigh_AF.acc(1:lastValidInd,:));
    
    % Find flexion offset
    [thighRotationAxis, thighRotationAngle] = vec2helic([motionlessAvgAcc(1) motionlessAvgAcc(2) 0], gravityVector);
    rotMatrix = quat2matrix(helic2quat(thighRotationAxis, thighRotationAngle));
    thighRotMatrix = rotMatrix * thighRotMatrix;


    
    %% Get true flexion offset
    
    % Compute sacrum orientation during squats
    squatSacrum_AF.acc = squatSacrum_LF.acc * sacrumRotMatrix';
    squatSacrum_AF.gyr = squatSacrum_LF.gyr * sacrumRotMatrix';
    sacrumOrientation = strapdownStaticDriftCorrection(squatSacrum_AF, gravityVector, fs);

    % Compute thigh orientation during squats
    squatThigh_AF.acc = squatThigh_LF.acc * thighRotMatrix';
    squatThigh_AF.gyr = squatThigh_LF.gyr * thighRotMatrix';
    thighOrientation = strapdownStaticDriftCorrection(squatThigh_AF, gravityVector, fs);

    % Translate sacrum acc to hip joint
    
    % Sensor-to-hip-joint-centre vectors
    if strcmp(side, 'left')
        rSacrumHip = [0.05 -0.10 0]; % z-axis not really important since we will only work in sagittal plane
        rThighHip = [-0.05 0.3 0];
    else
        rSacrumHip = [0.05 -0.10 0];
        rThighHip = [-0.05 0.3 0];
    end
    
    squatSacrumAtHip.acc = translateAcc(squatSacrum_AF, rSacrumHip, gravityValue, fs);
    squatSacrumAtHip.gyr = squatSacrum_AF.gyr; % Solid body: angular velocity does not change
    
    squatThighAtHip.acc = translateAcc(squatThigh_AF, rThighHip, gravityValue, fs);
    squatThighAtHip.gyr = squatThigh_AF.gyr; 
    
    % Determine the segment flexion offset that minimizes the acceleration
    % difference at the hip joint centre
    flexionOffsetQuat = determineFlexionOffset(squatSacrumAtHip, sacrumOrientation, squatThighAtHip, thighOrientation, gravityValue, fs);

    % Print for debugging
    fprintf('    Thigh flexion optimization: %1.2fdeg offset\n', 2*asind(flexionOffsetQuat(4)));

    % Correct
    thighRotMatrix = quat2matrix(flexionOffsetQuat) * thighRotMatrix;
end