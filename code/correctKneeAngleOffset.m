% correctKneeAngleOffset Corrects knee flexion and abduction offsets
%    [shankRotMatrix, thighRotMatrix] = correctKneeAngleOffset(imuData, calibInfo, side, shankRotMatrix, thighRotMatrix, gravityVector, fs)
%    corrects the knee flexion and abduction offsets. The assumption is
%    zero average knee flexion during the abduction movements and zero knee
%    abduction during the upright standing.
%    Returned are the updated calibration matrices for both the shank and
%    thigh.

function [shankRotMatrix, thighRotMatrix] = correctKneeAngleOffset(imuData, calibInfo, side, shankRotMatrix, thighRotMatrix, gravityVector, fs)

    % Load data
    if strcmpi(side, 'left')
        % Abduction
        abductionThigh_LF.acc = imuData.leftThigh.acc(calibInfo.left_abduction_start:calibInfo.left_abduction_stop,:);
        abductionThigh_LF.gyr = imuData.leftThigh.gyr(calibInfo.left_abduction_start:calibInfo.left_abduction_stop,:);
        abductionShank_LF.acc = imuData.leftShank.acc(calibInfo.left_abduction_start:calibInfo.left_abduction_stop,:);
        abductionShank_LF.gyr = imuData.leftShank.gyr(calibInfo.left_abduction_start:calibInfo.left_abduction_stop,:);
        
        % Upright
        uprightThigh_LF.acc = imuData.leftThigh.acc(calibInfo.upright_start:calibInfo.upright_stop,:);
        uprightThigh_LF.gyr = imuData.leftThigh.gyr(calibInfo.upright_start:calibInfo.upright_stop,:);
        uprightShank_LF.acc = imuData.leftShank.acc(calibInfo.upright_start:calibInfo.upright_stop,:);
        uprightShank_LF.gyr = imuData.leftShank.gyr(calibInfo.upright_start:calibInfo.upright_stop,:);
    else
        % Abduction
        abductionThigh_LF.acc = imuData.rightThigh.acc(calibInfo.right_abduction_start:calibInfo.right_abduction_stop,:);
        abductionThigh_LF.gyr = imuData.rightThigh.gyr(calibInfo.right_abduction_start:calibInfo.right_abduction_stop,:);
        abductionShank_LF.acc = imuData.rightShank.acc(calibInfo.right_abduction_start:calibInfo.right_abduction_stop,:);
        abductionShank_LF.gyr = imuData.rightShank.gyr(calibInfo.right_abduction_start:calibInfo.right_abduction_stop,:);
        
        % Upright
        uprightThigh_LF.acc = imuData.rightThigh.acc(calibInfo.upright_start:calibInfo.upright_stop,:);
        uprightThigh_LF.gyr = imuData.rightThigh.gyr(calibInfo.upright_start:calibInfo.upright_stop,:);
        uprightShank_LF.acc = imuData.rightShank.acc(calibInfo.upright_start:calibInfo.upright_stop,:);
        uprightShank_LF.gyr = imuData.rightShank.gyr(calibInfo.upright_start:calibInfo.upright_stop,:);
    end

    
    
    %% Impose 0 deg average knee flexion during entire abduction movement
    
    % Convert to "anatomical" frame
    abductionShank_AF.acc = abductionShank_LF.acc * shankRotMatrix';
    abductionShank_AF.gyr = abductionShank_LF.gyr * shankRotMatrix';
    abductionThigh_AF.acc = abductionThigh_LF.acc * thighRotMatrix';
    abductionThigh_AF.gyr = abductionThigh_LF.gyr * thighRotMatrix';
    
    % Compute orientation
    shankOrientation = strapdownStaticDriftCorrection(abductionShank_AF, gravityVector, fs);
    thighOrientation = strapdownStaticDriftCorrection(abductionThigh_AF, gravityVector, fs);
    
    % Compute 3D knee angles
    kneeAngles = computeKneeAngles(shankOrientation, thighOrientation, side);
    
    % Get average flexion angles
    averageFlexion = mean(kneeAngles(:,1));
    
    % Correct flexion offset half on thigh, half on shank
    kneeFlexionOffset = averageFlexion/2;
    shankRotMatrix = zRotation(-kneeFlexionOffset/180*pi) * shankRotMatrix;
    thighRotMatrix = zRotation(kneeFlexionOffset/180*pi) * thighRotMatrix;
    
    
    
    %% Impose 0 deg average knee abduction during entire upright posture
    
    % Convert to "anatomical" frame
    uprightShank_AF.acc = uprightShank_LF.acc * shankRotMatrix';
    uprightShank_AF.gyr = uprightShank_LF.gyr * shankRotMatrix';
    uprightThigh_AF.acc = uprightThigh_LF.acc * thighRotMatrix';
    uprightThigh_AF.gyr = uprightThigh_LF.gyr * thighRotMatrix';
    
    % Get the abduction angle (rotation around x axis) based on measured gravity
    tmp = mean(uprightShank_AF.acc);
    [shankAxis, shankAbduction] = vec2helic([0 tmp(2:3)], gravityVector);
    shankCorrection = quat2matrix(helic2quat(shankAxis, shankAbduction));
    
    tmp = mean(uprightThigh_AF.acc);
    [thighAxis, thighAbduction] = vec2helic([0 tmp(2:3)], gravityVector);
    thighCorrection = quat2matrix(helic2quat(thighAxis, thighAbduction));
    
    fprintf('    Corrected %s shank abduction %1.2fdeg; %s thigh abduction %1.2fdeg\n', side, shankAbduction/pi*180, side, thighAbduction/pi*180);
    
    shankRotMatrix  = shankCorrection * shankRotMatrix;
    thighRotMatrix  = thighCorrection * thighRotMatrix;
end