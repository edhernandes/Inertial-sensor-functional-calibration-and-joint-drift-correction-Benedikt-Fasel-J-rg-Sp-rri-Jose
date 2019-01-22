% symmetrizeLowerlimbs Ensures symmetric orientation of the lower limbs
%    calibMatrix = symmetrizeLowerlimbs(imuData, calibInfo, calibMatrix, gravityVector, fs)
%    ensures that the shanks' and thighs' inclination (flexion) are
%    identical during the static period at the beginning of the squat
%    movements.

function calibMatrix = symmetrizeLowerlimbs(imuData, calibInfo, calibMatrix, gravityVector, fs)

    % Load data
    squats_LF.leftThigh.acc = imuData.leftThigh.acc(calibInfo.squat_start:calibInfo.squat_stop,:);
    squats_LF.leftThigh.gyr = imuData.leftThigh.gyr(calibInfo.squat_start:calibInfo.squat_stop,:);
    squats_LF.leftShank.acc = imuData.leftShank.acc(calibInfo.squat_start:calibInfo.squat_stop,:);
    squats_LF.leftShank.gyr = imuData.leftShank.gyr(calibInfo.squat_start:calibInfo.squat_stop,:);
    
    squats_LF.rightThigh.acc = imuData.rightThigh.acc(calibInfo.squat_start:calibInfo.squat_stop,:);
    squats_LF.rightThigh.gyr = imuData.rightThigh.gyr(calibInfo.squat_start:calibInfo.squat_stop,:);
    squats_LF.rightShank.acc = imuData.rightShank.acc(calibInfo.squat_start:calibInfo.squat_stop,:);
    squats_LF.rightShank.gyr = imuData.rightShank.gyr(calibInfo.squat_start:calibInfo.squat_stop,:);
    

    % Convert into anatomical frame
    squats_AF.leftThigh.acc = squats_LF.leftThigh.acc * calibMatrix.leftThigh';
    squats_AF.leftThigh.gyr = squats_LF.leftThigh.gyr * calibMatrix.leftThigh';
    squats_AF.leftShank.acc = squats_LF.leftShank.acc * calibMatrix.leftShank';
    squats_AF.leftShank.gyr = squats_LF.leftShank.gyr * calibMatrix.leftShank';
    
    squats_AF.rightThigh.acc = squats_LF.rightThigh.acc * calibMatrix.rightThigh';
    squats_AF.rightThigh.gyr = squats_LF.rightThigh.gyr * calibMatrix.rightThigh';
    squats_AF.rightShank.acc = squats_LF.rightShank.acc * calibMatrix.rightShank';
    squats_AF.rightShank.gyr = squats_LF.rightShank.gyr * calibMatrix.rightShank';
    

    % Get the flexion angle (rotation around z axis) at the static period
    % at the beginning of the squats. Take the first 0.3 seconds
    
    % Shanks
    tmp = mean(squats_AF.leftShank.acc(1:round(0.3*fs),:));
    [leftAxis, leftShankFlexion] = vec2helic([tmp(1:2) 0], gravityVector);
    if leftAxis(3)<0
        leftShankFlexion = -leftShankFlexion;
    end
    
    tmp = mean(squats_AF.rightShank.acc(1:round(0.3*fs),:));
    [rightAxis, rightShankFlexion] = vec2helic([tmp(1:2) 0], gravityVector);
    if rightAxis(3)<0
        rightShankFlexion = -rightShankFlexion;
    end
    
    fprintf('    Left shank flexion %1.2fdeg; right shank flexion %1.2fdeg\n', leftShankFlexion/pi*180, rightShankFlexion/pi*180);
    
    % Thighs
    tmp = mean(squats_AF.leftThigh.acc(1:round(0.3*fs),:));
    [leftAxis, leftThighFlexion] = vec2helic([tmp(1:2) 0], gravityVector);
    if leftAxis(3)<0
        leftThighFlexion = -leftThighFlexion;
    end
    
    tmp = mean(squats_AF.rightThigh.acc(1:round(0.3*fs),:));
    [rightAxis, rightThighFlexion] = vec2helic([tmp(1:2) 0], gravityVector);
    if rightAxis(3)<0
        rightThighFlexion = -rightThighFlexion;
    end
    
    fprintf('    Left thigh flexion %1.2fdeg; right thigh flexion %1.2fdeg\n', leftThighFlexion/pi*180, rightThighFlexion/pi*180);

    
    % Correct flexion offset: set to mean angle
    meanShankFlexion = 0.5*(leftShankFlexion+rightShankFlexion);
    
    leftShankCorrection = meanShankFlexion - leftShankFlexion;
    rightShankCorrection = meanShankFlexion - rightShankFlexion;
    
    meanThighFlexion = 0.5*(leftThighFlexion+rightThighFlexion);
    
    leftThighCorrection = meanThighFlexion - leftThighFlexion;
    rightThighCorrection = meanThighFlexion - rightThighFlexion;
    
    
    % Apply the correction
    calibMatrix.leftShank  = zRotation(-leftShankCorrection) * calibMatrix.leftShank;
    calibMatrix.leftThigh  = zRotation(-leftThighCorrection) * calibMatrix.leftThigh;
    calibMatrix.rightShank  = zRotation(-rightShankCorrection) * calibMatrix.rightShank;
    calibMatrix.rightThigh  = zRotation(-rightThighCorrection) * calibMatrix.rightThigh;
end