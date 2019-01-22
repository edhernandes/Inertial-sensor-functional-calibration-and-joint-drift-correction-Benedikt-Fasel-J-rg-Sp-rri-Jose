% computeDistances Computes the distances between all sensors and joint centers
%    dists = computeDistances(imu, gravityVector, gravityValue, fs)
%    computes all distances between the sensors present in the structure
%    imu and their respective joint centers (e.g. for the thigh sensor the
%    distances to the knee and hip joint are computes).
%    Returned are the distances in m.

function dists = computeDistances(imu, gravityVector, gravityValue, fs)

    if isfield(imu, 'leftThigh') && isfield(imu, 'leftShank')
        [dists.leftThighKnee, dists.leftShankKnee] = findJointDistances(imu.leftThigh, [-0.03 -0.15 0.03], imu.leftShank, [-0.01 0.1 -0.02], 'left knee', gravityVector, gravityValue, fs);
    end

    if isfield(imu, 'rightThigh') && isfield(imu, 'rightShank')
        [dists.rightThighKnee, dists.rightShankKnee] = findJointDistances(imu.rightThigh, [-0.03 -0.15 -0.03], imu.rightShank, [-0.01 0.1 0.02], 'right knee', gravityVector, gravityValue, fs);
    end

    if isfield(imu, 'sacrum') && isfield(imu, 'leftThigh')
        [dists.sacrumLeftHip, dists.leftThighHip] = findJointDistances(imu.sacrum, [0.05 -0.05 -0.12], imu.leftThigh, [-0.03 0.3 0.03], 'left hip', gravityVector, gravityValue, fs);
    end

    if isfield(imu, 'sacrum') && isfield(imu, 'rightThigh')
        [dists.sacrumRightHip, dists.rightThighHip] = findJointDistances(imu.sacrum, [0.05 -0.05 0.12], imu.rightThigh, [-0.03 0.3 -0.03], 'right hip', gravityVector, gravityValue, fs);
    end

    
    % Trunk
    if isfield(imu, 'sacrum') && isfield(imu, 'sternum')
        [dists.sacrumTrunkCenterSt, dists.sternumTrunkCenterSt] = findJointDistances(imu.sacrum, [0.02 0.2 0.0], imu.sternum, [-0.1 -0.2 0.0], 'trunk center (sacr / strn)', gravityVector, gravityValue, fs);
    end

    if isfield(imu, 'sternum') && isfield(imu, 'head')
        [dists.sternumNeck, dists.headNeckSt] = findJointDistances(imu.sternum, [-0.1 0.3 0.0], imu.head, [0.0 -0.17 0.0], 'neck (strn / head)', gravityVector, gravityValue, fs);
    end    


    % Upper limbs
    if isfield(imu, 'sternum') && isfield(imu, 'leftArm')
        [dists.sternumLeftShoulder, dists.leftArmShoulder] = findJointDistances(imu.sternum, [-0.05 0.15 -0.15], imu.leftArm, [-0.02 0.2 0.0], 'left shoulder', gravityVector, gravityValue, fs);
    end

    if isfield(imu, 'sternum') && isfield(imu, 'rightArm')
        [dists.sternumRightShoulder, dists.rightArmShoulder] = findJointDistances(imu.sternum, [-0.05 0.15 0.15], imu.rightArm, [-0.02 0.2 0.0], 'right shoulder', gravityVector, gravityValue, fs);
    end
    if isfield(imu, 'leftArm') && isfield(imu, 'leftWrist')
        [dists.leftArmElbow, dists.leftWristElbow] = findJointDistances(imu.leftArm, [0.02 -0.15 0.0], imu.leftWrist, [0.0 0.2 0.0], 'left elbow', gravityVector, gravityValue, fs);
    end

    if isfield(imu, 'rightArm') && isfield(imu, 'rightWrist')
        [dists.rightArmElbow, dists.rightWristElbow] = findJointDistances(imu.rightArm, [0.02 -0.15 0.0], imu.rightWrist, [0.0 0.2 0.0], 'right elbow', gravityVector, gravityValue, fs);
    end

    
    fprintf('\nFinished computing distances\n\n');
end