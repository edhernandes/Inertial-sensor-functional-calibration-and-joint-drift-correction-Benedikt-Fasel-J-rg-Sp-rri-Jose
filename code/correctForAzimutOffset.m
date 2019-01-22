% correctForAzimuthOffset Corrects initial azimuth offset
%    imu = correctForAzimutOffset(imu, gravityValue, fs) corrects an
%    initial azimuth offset that could be present due to a non-perfect
%    upright posture at initialization of strap-down. Only the upper limb
%    sensors are corrected
% 
%    Returned is the corrected orientation for each sensor

function imu = correctForAzimutOffset(imu, gravityValue, fs)

    % Take approximative default values
    dists.sternumLeftShoulder = [-0.05 0.15 -0.15];
    dists.sternumRightShoulder = [-0.05 0.15 0.15];
    dists.leftArmShoulder = [0 0.15 0.03];
    dists.rightArmShoulder = [0 0.15 -0.03];
    
    dists.leftArmElbow = [0 -0.1 0.03];
    dists.rightArmElbow = [0 -0.1 -0.03];
    dists.leftWristElbow = [0 0.15 0];
    dists.rightWristElbow = [0 0.15 0];
    
    
    dists.headNeck = [-0.05 -0.25 0];
    dists.sternumNeck = [-0.08 0.2 0];
    
    % Correct the right arm sensor. Assume that sternum is correct.
    if isfield(imu, 'sternum') && isfield(imu, 'rightArm')
        fprintf('Estimating azimut offset for right arm sensor... ');

        % Compute virtual acceleration at shoulder based on the sternum sensor
        accSternum_atShoulder_wGrav = translateAcc(imu.sternum, dists.sternumRightShoulder, gravityValue, fs);
        
        % Compute virtual acceleration at shoulder based on the arm sensor
        accRightArm_atShoulder_wGrav = translateAcc(imu.rightArm, dists.rightArmShoulder, gravityValue, fs);
        
        % Find average azimuth orientation difference such that both
        % accelerations would be aligned
        rightArmOffsetQuat = determineAzimutOffset(accSternum_atShoulder_wGrav, imu.sternum.gyr, imu.sternum.orientation, accRightArm_atShoulder_wGrav, imu.rightArm.gyr, imu.rightArm.orientation, gravityValue, fs);

        % Print
        fprintf('%1.3f deg offset.\n', 2*asind(rightArmOffsetQuat(3)));
        
        % Correct
        for i=1:size(imu.rightArm.orientation,1)
            imu.rightArm.orientation(i,:) = quat_multiply(rightArmOffsetQuat, imu.rightArm.orientation(i,:));
        end
    end
    
    % Same for the left arm sensor
    if isfield(imu, 'sternum') && isfield(imu, 'leftArm')
        fprintf('Estimating azimut offset for left arm sensor... ');

        accSternum_atShoulder_wGrav = translateAcc(imu.sternum, dists.sternumLeftShoulder, gravityValue, fs);
        accLeftArm_atShoulder_wGrav = translateAcc(imu.leftArm, dists.leftArmShoulder, gravityValue, fs);
        leftArmOffsetQuat = determineAzimutOffset(accSternum_atShoulder_wGrav, imu.sternum.gyr, imu.sternum.orientation, accLeftArm_atShoulder_wGrav, imu.leftArm.gyr, imu.leftArm.orientation, gravityValue, fs);

        % Print
        fprintf('%1.3f deg offset.\n', 2*asind(leftArmOffsetQuat(3)));
        
        % Correct
        for i=1:size(imu.leftArm.orientation,1)
            imu.leftArm.orientation(i,:) = quat_multiply(leftArmOffsetQuat, imu.leftArm.orientation(i,:));
        end
    end
       
    % Now assume that right arm is correct and correct the forearm/wrist
    % sensor relative to the arm sensor.
    if isfield(imu, 'rightArm') && isfield(imu, 'rightWrist')
        fprintf('Estimating azimut offset for right wrist sensor... ');

        accRightArm_atElbow_wGrav = translateAcc(imu.rightArm, dists.rightArmElbow, gravityValue, fs);
        accRightWrist_atElbow_wGrav = translateAcc(imu.rightWrist, dists.rightWristElbow, gravityValue, fs);
        rightForearmOffsetQuat = determineAzimutOffset(accRightArm_atElbow_wGrav, imu.rightArm.gyr, imu.rightArm.orientation, accRightWrist_atElbow_wGrav, imu.rightWrist.gyr, imu.rightWrist.orientation, gravityValue, fs);

        % Print
        fprintf('%1.3f deg offset.\n', 2*asind(rightForearmOffsetQuat(3)));
        
        % Correct
        for i=1:size(imu.rightWrist.orientation,1)
            imu.rightWrist.orientation(i,:) = quat_multiply(rightForearmOffsetQuat, imu.rightWrist.orientation(i,:));
        end
    end
    
    % Same for left side
    if isfield(imu, 'leftArm') && isfield(imu, 'leftWrist')
        fprintf('Estimating azimut offset for left wrist sensor... ');

        accLeftArm_atElbow_wGrav = translateAcc(imu.leftArm, dists.leftArmElbow, gravityValue, fs);
        accLeftWrist_atElbow_wGrav = translateAcc(imu.leftWrist, dists.leftWristElbow, gravityValue, fs);
        leftForearmOffsetQuat = determineAzimutOffset(accLeftArm_atElbow_wGrav, imu.leftArm.gyr, imu.leftArm.orientation, accLeftWrist_atElbow_wGrav, imu.leftWrist.gyr, imu.leftWrist.orientation, gravityValue, fs);

        % Print
        fprintf('%1.3f deg offset.\n', 2*asind(leftForearmOffsetQuat(3)));
        
        % Correct
        for i=1:size(imu.leftWrist.orientation,1)
            imu.leftWrist.orientation(i,:) = quat_multiply(leftForearmOffsetQuat, imu.leftWrist.orientation(i,:));
        end
    end
    
    % Head
    if isfield(imu, 'head') && isfield(imu, 'sternum')
        fprintf('Estimating azimut offset for head sensor... ');
        
        accHead_atNeck_wGrav = translateAcc(imu.head, dists.headNeck, gravityValue, fs);
        accSternum_atNeck_wGrav = translateAcc(imu.sternum, dists.sternumNeck, gravityValue, fs);
        headOffsetQuat = determineAzimutOffset(accHead_atNeck_wGrav, imu.head.gyr, imu.head.orientation, accSternum_atNeck_wGrav, imu.sternum.gyr, imu.sternum.orientation, gravityValue, fs);
        
        % Print
        fprintf('%1.3f deg offset.\n', 2*asind(headOffsetQuat(3)));
        
        % Correct
        for i=1:size(imu.head.orientation,1)
            imu.head.orientation(i,:) = quat_multiply(headOffsetQuat, imu.head.orientation(i,:));
        end
    end
end