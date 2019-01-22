% spineCalibration Calibrate any spine/trunk-fixed IMU
%     rotMatrix = spineCalibration(imuData, calibInfo, sensorName, gravityVector, fs)
%     finds the rotation matrix that orients a spine/trunk-fixed sensor to 
%     its corresponding functional frame.
%
%     Assumptions: Squats take place in the saggital plane where the main
%                  rotation happens around the medio-lateral axis. 
%     Algorithm: First the medio-lateral axis is estimated and second the 
%                vertical axis. Finally any remaining abduction and flexion 
%                offsets are corrected based on the upright posture.

function spineRotMatrix = spineCalibration(imuData, calibInfo, sensorName, gravityVector, fs)
   
       
    % Lowpass filter parameters
    [bNormal,aNormal] = butter(2, 1.2465*5/fs*2, 'low');  % cutoff freq. 5 Hz
    [bLow,aLow] = butter(2, 1.2465*1/fs*2, 'low');        % cutoff freq. 1 Hz 
    
    
    %% Find azimut alignment rotation matrix
    
    % Hypothesis:
    % Principal axis of rotation during the squats corresponds to the
    % medio-lateral Z-axis [0 0 1].
    
    % Get acceleration and angular velocity for the squats
    squats.acc = imuData.(sensorName).acc(calibInfo.squat_start:calibInfo.squat_stop,:);
    squats.gyr = imuData.(sensorName).gyr(calibInfo.squat_start:calibInfo.squat_stop,:);

    % Low pass filter angular velocity to remove noise
    squatAngVel = filtfilt3D(bNormal, aNormal, squats.gyr);
    squatAngVelNorm = sqrt(sum(squatAngVel.^2,2));
    
    % To improve estimation accuracy we want to remove all instants of very
    % low rotation speeds (since this low rotation speed could by in any direction)
    % The angular velocity threshold was fixed at 30 deg/sec
    angVelThres = 30;
    
    % Compute the PCA to find the principal axis of rotation
    principalAxes = pca(squatAngVel(squatAngVelNorm>angVelThres,:));

    % Find the rotation that aligns the principal axis to the Z-axis
    [rotationAxis, rotationAngle] = vec2helic(principalAxes(:,1), [0 0 1]);
    spineRotMatrix = quat2matrix(helic2quat(rotationAxis, rotationAngle));


    % Check orientation of the z axis (must point to the right)
    % Hypothesis: first peak of the angular velocity in flexion must be
    %             negative (trunk forwards bending)
    
    % Align the angular velocity according to the orientation found above
    squatAligned = squatAngVel * spineRotMatrix';
    
    % Massively lowpass filter because we only want global oscillation pattern
    flexionLow = filtfilt(bLow,aLow,squatAligned(:,3));
    
    % Find all the angular velocity peaks (positive and negative)
    % To be insensitive to the movement intensity: chose a variable peak
    % height (75% of the standard deviation) for the peak detection.
    [~, locs] = findpeaks(abs(flexionLow), 'minpeakheight', 0.75*std(flexionLow));
    
    % First peak positive ==> Axis is inverted (points to the left instead
    % of to the right). Rotate by 180° around the "vertical" axis
    if flexionLow(locs(1))>0
        spineRotMatrix = yRotation(pi) * spineRotMatrix;
    end
    
    
    %% Find the vertical axis
    
    % Get the data for the trunk rotation (LF = local frame)
    rotation_LF.acc = imuData.(sensorName).acc(calibInfo.rotation_start:calibInfo.rotation_stop,:);
    rotation_LF.gyr = imuData.(sensorName).gyr(calibInfo.rotation_start:calibInfo.rotation_stop,:);
    
    % Transform to the frame found above (AF = anatomical frame)
    rotation_AF.acc = rotation_LF.acc * spineRotMatrix';
    rotation_AF.gyr = rotation_LF.gyr * spineRotMatrix';
    
    % Lowpass filter
    rotationAngVel = filtfilt3D(bNormal, aNormal, rotation_AF.gyr);
    rotationAngVelNorm = sqrt(sum(rotationAngVel.^2,2));
    
    % Same as for the squats: we only want to consider moments with "large"
    % rotational speeds. Since the vertical trunk movement is slower, lower
    % the threshold to 20 deg/sec
    angVelThres = 20;
    
    % Find the principal axes. Only consider the X- and Y-axis since the
    % Z-axis is already aligned.
    principalAxes = pca(rotationAngVel(rotationAngVelNorm>angVelThres,1:2));
    
    % Find rotation matrix to rotate the princiapal axis to the Y-axis
    [rotationAxis, rotationAngle] = vec2helic([principalAxes(:,1); 0], [0 1 0]);
    spineRotMatrix = quat2matrix(helic2quat(rotationAxis, rotationAngle)) * spineRotMatrix;
    
    % Check orientation of the y axis
    % Hypothesis: must be approximately aligned with gravity during the 
    %             entire trunk rotation movement
    
    % Align the acceleration to the updated anatomical frame
    rotation_AF.acc = rotation_LF.acc * spineRotMatrix'; 
    
    % If mean acceleration is <0 ==> Axis is inverted (points downwards
    % instead of upwards). Rotate by 180° around the medio-lateral axis
    if mean(rotation_AF.acc(:,2))<0
        spineRotMatrix = zRotation(pi) * spineRotMatrix;
    end
    
    
    %% Correction during upright
    
    % The movements may not have been perfectly executed: refine
    % calibration based on the upright posture
    
    
    % Constraint 1: zero abuction
    % ---------------------------
    
    % Get the data for the upright posture
    upright_LF.acc = imuData.(sensorName).acc(calibInfo.upright_start:calibInfo.upright_stop,:);
    upright_LF.gyr = imuData.(sensorName).gyr(calibInfo.upright_start:calibInfo.upright_stop,:);
    
    % Transform to the frame found above 
    upright_AF.acc = upright_LF.acc * spineRotMatrix';
    upright_AF.gyr = upright_LF.gyr * spineRotMatrix';

    
    % Get the abduction angle (rotation around x axis)
    tmp = mean(upright_AF.acc);
    [axis, abduction] = vec2helic([0 tmp(2:3)], gravityVector);
    abductionCorrection = quat2matrix(helic2quat(axis, abduction));
       
    fprintf('    %s abduction correction %1.2fdeg\n', sensorName, abduction/pi*180);
    
    % Update the calibration matrix
    spineRotMatrix  = abductionCorrection * spineRotMatrix;
    
    
    
    % Constraint 2: zero flexion
    % --------------------------
    
    % Update the "anatomical" frame
    upright_AF.acc = upright_LF.acc * spineRotMatrix';
    upright_AF.gyr = upright_LF.gyr * spineRotMatrix';
    
    % Get the flexion angle (rotation around z axis)
    tmp = mean(upright_AF.acc);
    [axis, flexion] = vec2helic([tmp(1:2) 0], gravityVector);
    flexionCorrection = quat2matrix(helic2quat(axis, flexion));
       
    fprintf('    %s flexion correction %1.2fdeg\n', sensorName, flexion/pi*180);
    
    spineRotMatrix  = flexionCorrection * spineRotMatrix;
    
end