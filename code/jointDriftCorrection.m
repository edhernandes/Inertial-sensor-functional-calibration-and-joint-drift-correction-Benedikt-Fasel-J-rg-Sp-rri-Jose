function imu = jointDriftCorrection(imu, dists, gravityValue, fs)

    % Set this variable to false if no second azimut correction shall be
    % done
    azimutCorrection = true;

    % Use sacrum as reference sensor. Its drift is given by the average
    % drift from left and right hip and trunk center joint


    %% Drift sacrum sensor


    % In order for the code to work well all three sensors must
    % be available
    if ~isfield(imu, 'sacrum') || ~isfield(imu, 'leftThigh') || ~isfield(imu, 'rightThigh') || ~isfield(imu, 'sternum')
        error('For the joint drift correction method to work properly all the four sensors sacrum, sternum, left thigh, right thigh must be present. One or multiple sensors are missing.');
    end

    
    fprintf('Estimating joint drift for sacrum sensor... ');

    % Left hip drift
    if isfield(imu, 'leftThigh')
        
        % Translate the accelerations to the common joint
        accSacrum_atHip_wGrav = translateAcc(imu.sacrum, dists.sacrumLeftHip, gravityValue, fs); % return unit: g
        accThigh_atHip_wGrav = translateAcc(imu.leftThigh, dists.leftThighHip, gravityValue, fs); % return unit: g
        
        % Estimate the joint drift (relative orientation difference between
        % the two translated accelerations)
        leftHipDriftQ = estimateJointDrift(accSacrum_atHip_wGrav, imu.sacrum.gyr, imu.sacrum.orientation, accThigh_atHip_wGrav, imu.leftThigh.gyr, imu.leftThigh.orientation, gravityValue, fs);
    end
    
    % Right hip drift
    if isfield(imu, 'rightThigh')
        accSacrum_atHip_wGrav = translateAcc(imu.sacrum, dists.sacrumRightHip, gravityValue, fs);
        accThigh_atHip_wGrav = translateAcc(imu.rightThigh, dists.rightThighHip, gravityValue, fs);
        rightHipDriftQ = estimateJointDrift(accSacrum_atHip_wGrav, imu.sacrum.gyr, imu.sacrum.orientation, accThigh_atHip_wGrav, imu.rightThigh.gyr, imu.rightThigh.orientation, gravityValue, fs);
    end
    
    % Trunk center with sternum
    if isfield(imu, 'sternum')
        accSacrum_atTrunk_wGrav = translateAcc(imu.sacrum, dists.sacrumTrunkCenterSt, gravityValue, fs);
        accSternum_atTrunk_wGrav = translateAcc(imu.sternum, dists.sternumTrunkCenterSt, gravityValue, fs);
        trunkCenterDriftQ_st = estimateJointDrift(accSacrum_atTrunk_wGrav, imu.sacrum.gyr, imu.sacrum.orientation, accSternum_atTrunk_wGrav, imu.sternum.gyr, imu.sternum.orientation, gravityValue, fs);
    end
    

    
    % Compute average drift 
    sacrumDrift = zeros(size(imu.sacrum.gyr,1),4);
    sacrumDrift(:,1) = 1;
    for i=1:size(sacrumDrift,1)
        allDrifts = zeros(3, 4);
        
        if isfield(imu, 'leftThigh')
            [axis, angle] = quat2helic(leftHipDriftQ(i,:));
            allDrifts(1,:) = helic2quat(axis, -angle);
        end
        
        if isfield(imu, 'rightThigh')
            [axis, angle] = quat2helic(rightHipDriftQ(i,:));
            allDrifts(2,:) = helic2quat(axis, -angle);
        end
        
        if isfield(imu, 'sternum')
            [axis, angle] = quat2helic(trunkCenterDriftQ_st(i,:));
            allDrifts(3,:) = helic2quat(axis, -angle);
        end
        
        if ~any(isnan(allDrifts(:,1)))
            sacrumDrift(i,:) = quat_mean(allDrifts);
        end
    end
    

    % Lowpass filter
    [b,a] = butter(2, 1.2465*2/fs*2, 'low');
    for i=1:4
        sacrumDrift(:,i) = filtfilt(b,a,sacrumDrift(:,i));
    end

    % Normalize
    for i=1:size(sacrumDrift,1)
        sacrumDrift(i,:) = quat_normalize(sacrumDrift(i,:));
    end

    
    % Correct drift and update angular velocity
    imu.sacrum.orientation = correctSegmentDrift(imu.sacrum.orientation, sacrumDrift);
    sacrumAngVel = inverseStrapdown(fs, imu.sacrum.orientation)/pi*180;
    imu.sacrum.gyr(~isnan(sacrumAngVel(:,1)),:) = sacrumAngVel(~isnan(sacrumAngVel(:,1)),:);
    

    fprintf('[done]\n');


    %% Joint drift correction along kinematic chain

    % Drift of left and right thigh
    % -----------------------------

    fprintf('Estimating joint drift for left and right thigh sensor... ');

    % Left hip drift
    if isfield(imu, 'leftThigh')
        accSacrum_atHip_wGrav = translateAcc(imu.sacrum, dists.sacrumLeftHip, gravityValue, fs);
        accThigh_atHip_wGrav = translateAcc(imu.leftThigh, dists.leftThighHip, gravityValue, fs);
        leftHipDriftQ = estimateJointDrift(accSacrum_atHip_wGrav, imu.sacrum.gyr, imu.sacrum.orientation, accThigh_atHip_wGrav, imu.leftThigh.gyr, imu.leftThigh.orientation, gravityValue, fs);
    
        if ~isempty(leftHipDriftQ)
            imu.leftThigh.orientation = correctSegmentDrift(imu.leftThigh.orientation, leftHipDriftQ);
            leftThighAngVel = inverseStrapdown(fs, imu.leftThigh.orientation)/pi*180;
            imu.leftThigh.gyr(~isnan(leftThighAngVel(:,1)),:) = leftThighAngVel(~isnan(leftThighAngVel(:,1)),:);
        end
        
        if azimutCorrection
            accSacrum_atHip_wGrav = translateAcc(imu.sacrum, dists.sacrumLeftHip, gravityValue, fs);
            accThigh_atHip_wGrav = translateAcc(imu.leftThigh, dists.leftThighHip, gravityValue, fs);
            leftHipDriftQ = estimateJointDrift_AzimutOnly(accSacrum_atHip_wGrav, imu.sacrum.gyr, imu.sacrum.orientation, accThigh_atHip_wGrav, imu.leftThigh.gyr, imu.leftThigh.orientation, gravityValue, fs);

            if ~isempty(leftHipDriftQ)
                imu.leftThigh.orientation = correctSegmentDrift(imu.leftThigh.orientation, leftHipDriftQ);
                leftThighAngVel = inverseStrapdown(fs, imu.leftThigh.orientation)/pi*180;
                imu.leftThigh.gyr(~isnan(leftThighAngVel(:,1)),:) = leftThighAngVel(~isnan(leftThighAngVel(:,1)),:);
            end
        end
    end
    
    % Right hip drift
    if isfield(imu, 'rightThigh')
        accSacrum_atHip_wGrav = translateAcc(imu.sacrum, dists.sacrumRightHip, gravityValue, fs);
        accThigh_atHip_wGrav = translateAcc(imu.rightThigh, dists.rightThighHip, gravityValue, fs);
        rightHipDriftQ = estimateJointDrift(accSacrum_atHip_wGrav, imu.sacrum.gyr, imu.sacrum.orientation, accThigh_atHip_wGrav, imu.rightThigh.gyr, imu.rightThigh.orientation, gravityValue, fs);
    
        if ~isempty(rightHipDriftQ)
            imu.rightThigh.orientation = correctSegmentDrift(imu.rightThigh.orientation, rightHipDriftQ);
            rightThighAngVel = inverseStrapdown(fs, imu.rightThigh.orientation)/pi*180;
            imu.rightThigh.gyr(~isnan(rightThighAngVel(:,1)),:) = rightThighAngVel(~isnan(rightThighAngVel(:,1)),:);
        end
        
        if azimutCorrection
            accSacrum_atHip_wGrav = translateAcc(imu.sacrum, dists.sacrumRightHip, gravityValue, fs);
            accThigh_atHip_wGrav = translateAcc(imu.rightThigh, dists.rightThighHip, gravityValue, fs);
            rightHipDriftQ = estimateJointDrift_AzimutOnly(accSacrum_atHip_wGrav, imu.sacrum.gyr, imu.sacrum.orientation, accThigh_atHip_wGrav, imu.rightThigh.gyr, imu.rightThigh.orientation, gravityValue, fs);

            if ~isempty(rightHipDriftQ)
                imu.rightThigh.orientation = correctSegmentDrift(imu.rightThigh.orientation, rightHipDriftQ);
                rightThighAngVel = inverseStrapdown(fs, imu.rightThigh.orientation)/pi*180;
                imu.rightThigh.gyr(~isnan(rightThighAngVel(:,1)),:) = rightThighAngVel(~isnan(rightThighAngVel(:,1)),:);
            end
        end
    end
    
    fprintf('[done]\n');
    clear accSacrum_atHip_wGrav accThigh_atHip_wGrav rightHipDriftQ leftHipDriftQ
    clear rightThighAngVel leftThighAngVel

    


    % Drift of left and right shank
    % -----------------------------

    fprintf('Estimating joint drift for left and right shank sensor... ');

    % Left knee drift
    if isfield(imu, 'leftShank') && isfield(imu, 'leftThigh')
        accThigh_atKnee_wGrav = translateAcc(imu.leftThigh, dists.leftThighKnee, gravityValue, fs);
        accShank_atKnee_wGrav = translateAcc(imu.leftShank, dists.leftShankKnee, gravityValue, fs);
        leftKneeDriftQ = estimateJointDrift(accThigh_atKnee_wGrav, imu.leftThigh.gyr, imu.leftThigh.orientation, accShank_atKnee_wGrav, imu.leftShank.gyr, imu.leftShank.orientation, gravityValue, fs);

        if ~isempty(leftKneeDriftQ)
            imu.leftShank.orientation = correctSegmentDrift(imu.leftShank.orientation, leftKneeDriftQ);
            leftShankAngVel = inverseStrapdown(fs, imu.leftShank.orientation)/pi*180;
            imu.leftShank.gyr(~isnan(leftShankAngVel(:,1)),:) = leftShankAngVel(~isnan(leftShankAngVel(:,1)),:);
        end
        
        if azimutCorrection
            accThigh_atKnee_wGrav = translateAcc(imu.leftThigh, dists.leftThighKnee, gravityValue, fs);
            accShank_atKnee_wGrav = translateAcc(imu.leftShank, dists.leftShankKnee, gravityValue, fs);
            leftKneeDriftQ = estimateJointDrift_AzimutOnly(accThigh_atKnee_wGrav, imu.leftThigh.gyr, imu.leftThigh.orientation, accShank_atKnee_wGrav, imu.leftShank.gyr, imu.leftShank.orientation, gravityValue, fs);

            if ~isempty(leftKneeDriftQ)
                imu.leftShank.orientation = correctSegmentDrift(imu.leftShank.orientation, leftKneeDriftQ);
                leftShankAngVel = inverseStrapdown(fs, imu.leftShank.orientation)/pi*180;
                imu.leftShank.gyr(~isnan(leftShankAngVel(:,1)),:) = leftShankAngVel(~isnan(leftShankAngVel(:,1)),:);
            end
        end
    end
    
    % Right knee drift
    if isfield(imu, 'rightShank') && isfield(imu, 'rightThigh')
        accThigh_atKnee_wGrav = translateAcc(imu.rightThigh, dists.rightThighKnee, gravityValue, fs);
        accShank_atKnee_wGrav = translateAcc(imu.rightShank, dists.rightShankKnee, gravityValue, fs);
        rightKneeDriftQ = estimateJointDrift(accThigh_atKnee_wGrav, imu.rightThigh.gyr, imu.rightThigh.orientation, accShank_atKnee_wGrav, imu.rightShank.gyr, imu.rightShank.orientation, gravityValue, fs);

        if ~isempty(rightKneeDriftQ)
            imu.rightShank.orientation = correctSegmentDrift(imu.rightShank.orientation, rightKneeDriftQ);
            rightShankAngVel = inverseStrapdown(fs, imu.rightShank.orientation)/pi*180;
            imu.rightShank.gyr(~isnan(rightShankAngVel(:,1)),:) = rightShankAngVel(~isnan(rightShankAngVel(:,1)),:);
        end
        
        if azimutCorrection
            accThigh_atKnee_wGrav = translateAcc(imu.rightThigh, dists.rightThighKnee, gravityValue, fs);
            accShank_atKnee_wGrav = translateAcc(imu.rightShank, dists.rightShankKnee, gravityValue, fs);
            rightKneeDriftQ = estimateJointDrift_AzimutOnly(accThigh_atKnee_wGrav, imu.rightThigh.gyr, imu.rightThigh.orientation, accShank_atKnee_wGrav, imu.rightShank.gyr, imu.rightShank.orientation, gravityValue, fs);

            if ~isempty(rightKneeDriftQ)
                imu.rightShank.orientation = correctSegmentDrift(imu.rightShank.orientation, rightKneeDriftQ);
                rightShankAngVel = inverseStrapdown(fs, imu.rightShank.orientation)/pi*180;
                imu.rightShank.gyr(~isnan(rightShankAngVel(:,1)),:) = rightShankAngVel(~isnan(rightShankAngVel(:,1)),:);
            end
        end
    end
    
    fprintf('[done]\n');
    clear accThigh_atKnee_wGrav accShank_atKnee_wGrav rightKneeDriftQ leftKneeDriftQ
    clear rightShankAngVel leftShankAngVel

    

    % Drift of sternum
    % ----------------

    if isfield(imu, 'sternum')
        fprintf('Estimating joint drift for sternum sensor... ');
        
        accSacrum_atTrunk_wGrav = translateAcc(imu.sacrum, dists.sacrumTrunkCenterSt, gravityValue, fs);
        accSternum_atTrunk_wGrav = translateAcc(imu.sternum, dists.sternumTrunkCenterSt, gravityValue, fs);
        trunkCenterDriftQ = estimateJointDrift(accSacrum_atTrunk_wGrav, imu.sacrum.gyr, imu.sacrum.orientation, accSternum_atTrunk_wGrav, imu.sternum.gyr, imu.sternum.orientation, gravityValue, fs);

        if ~isempty(trunkCenterDriftQ)
            imu.sternum.orientation = correctSegmentDrift(imu.sternum.orientation, trunkCenterDriftQ);
            sternumAngVel = inverseStrapdown(fs, imu.sternum.orientation)/pi*180;
            imu.sternum.gyr(~isnan(sternumAngVel(:,1)),:) = sternumAngVel(~isnan(sternumAngVel(:,1)),:);
        end
        
        if azimutCorrection
            accSacrum_atTrunk_wGrav = translateAcc(imu.sacrum, dists.sacrumTrunkCenterSt, gravityValue, fs);
            accSternum_atTrunk_wGrav = translateAcc(imu.sternum, dists.sternumTrunkCenterSt, gravityValue, fs);
            trunkCenterDriftQ = estimateJointDrift_AzimutOnly(accSacrum_atTrunk_wGrav, imu.sacrum.gyr, imu.sacrum.orientation, accSternum_atTrunk_wGrav, imu.sternum.gyr, imu.sternum.orientation, gravityValue, fs);

            if ~isempty(trunkCenterDriftQ)
                imu.sternum.orientation = correctSegmentDrift(imu.sternum.orientation, trunkCenterDriftQ);
                sternumAngVel = inverseStrapdown(fs, imu.sternum.orientation)/pi*180;
                imu.sternum.gyr(~isnan(sternumAngVel(:,1)),:) = sternumAngVel(~isnan(sternumAngVel(:,1)),:);
            end
        end

        fprintf('[done]\n');
        clear accSacrum_atTrunk_wGrav accSternum_atTrunk_wGrav trunkCenterDriftQ sternumAngVel
    end
    

    % Drift of head
    % -------------

    % Sternum sensor available but not spine C7
    if isfield(imu, 'head') && isfield(imu, 'sternum') && ~isfield(imu, 'spineC7')
        fprintf('Estimating joint drift for head sensor... ');

        accSternum_atNeck_wGrav = translateAcc(imu.sternum, dists.sternumNeck, gravityValue, fs);
        accHead_atNeck_wGrav = translateAcc(imu.head, dists.headNeckSt, gravityValue, fs);
        neckDriftQ = estimateJointDrift(accSternum_atNeck_wGrav, imu.sternum.gyr, imu.sternum.orientation, accHead_atNeck_wGrav, imu.head.gyr, imu.head.orientation, gravityValue, fs);

        if ~isempty(neckDriftQ)
            imu.head.orientation = correctSegmentDrift(imu.head.orientation, neckDriftQ);
            headAngVel = inverseStrapdown(fs, imu.head.orientation)/pi*180;
            imu.head.gyr(~isnan(headAngVel(:,1)),:) = headAngVel(~isnan(headAngVel(:,1)),:);  
        end
        
        if azimutCorrection
            accSternum_atNeck_wGrav = translateAcc(imu.sternum, dists.sternumNeck, gravityValue, fs);
            accHead_atNeck_wGrav = translateAcc(imu.head, dists.headNeckSt, gravityValue, fs);
            neckDriftQ = estimateJointDrift_AzimutOnly(accSternum_atNeck_wGrav, imu.sternum.gyr, imu.sternum.orientation, accHead_atNeck_wGrav, imu.head.gyr, imu.head.orientation, gravityValue, fs);

            if ~isempty(neckDriftQ)
                imu.head.orientation = correctSegmentDrift(imu.head.orientation, neckDriftQ);
                headAngVel = inverseStrapdown(fs, imu.head.orientation)/pi*180;
                imu.head.gyr(~isnan(headAngVel(:,1)),:) = headAngVel(~isnan(headAngVel(:,1)),:);  
            end
        end

        fprintf('[done]\n');
        clear accSternum_atNeck_wGrav accHead_atNeck_wGrav neckDriftQ headAngVel
    end
    
    

    % Drift of right arm
    % -------------------

    if isfield(imu, 'sternum') && isfield(imu, 'rightArm')
        fprintf('Estimating joint drift for right arm sensor... ');

        accSternum_atShoulder_wGrav = translateAcc(imu.sternum, dists.sternumRightShoulder, gravityValue, fs);
        accRightArm_atShoulder_wGrav = translateAcc(imu.rightArm, dists.rightArmShoulder, gravityValue, fs);
        shoulderDriftQ = estimateJointDrift(accSternum_atShoulder_wGrav, imu.sternum.gyr, imu.sternum.orientation, accRightArm_atShoulder_wGrav, imu.rightArm.gyr, imu.rightArm.orientation, gravityValue, fs);

        if ~isempty(shoulderDriftQ)
            imu.rightArm.orientation = correctSegmentDrift(imu.rightArm.orientation, shoulderDriftQ);
            rightArmAngVel = inverseStrapdown(fs, imu.rightArm.orientation)/pi*180;
            imu.rightArm.gyr(~isnan(rightArmAngVel(:,1)),:) = rightArmAngVel(~isnan(rightArmAngVel(:,1)),:);
        end
        
        if azimutCorrection
            accSternum_atShoulder_wGrav = translateAcc(imu.sternum, dists.sternumRightShoulder, gravityValue, fs);
            accRightArm_atShoulder_wGrav = translateAcc(imu.rightArm, dists.rightArmShoulder, gravityValue, fs);
            shoulderDriftQ = estimateJointDrift_AzimutOnly(accSternum_atShoulder_wGrav, imu.sternum.gyr, imu.sternum.orientation, accRightArm_atShoulder_wGrav, imu.rightArm.gyr, imu.rightArm.orientation, gravityValue, fs, false);

            if ~isempty(shoulderDriftQ)
                imu.rightArm.orientation = correctSegmentDrift(imu.rightArm.orientation, shoulderDriftQ);
                rightArmAngVel = inverseStrapdown(fs, imu.rightArm.orientation)/pi*180;
                imu.rightArm.gyr(~isnan(rightArmAngVel(:,1)),:) = rightArmAngVel(~isnan(rightArmAngVel(:,1)),:);
            end
        end
        
        fprintf('[done]\n');
        clear accSternum_atShoulder_wGrav accRightArm_atShoulder_wGrav shoulderDriftQ rightArmAngVel
    end


    % Drift of left arm
    % -------------------

    if isfield(imu, 'sternum') && isfield(imu, 'leftArm')
        fprintf('Estimating joint drift for left arm sensor... ');

        accSternum_atShoulder_wGrav = translateAcc(imu.sternum, dists.sternumLeftShoulder, gravityValue, fs);
        accLeftArm_atShoulder_wGrav = translateAcc(imu.leftArm, dists.leftArmShoulder, gravityValue, fs);
        shoulderDriftQ = estimateJointDrift(accSternum_atShoulder_wGrav, imu.sternum.gyr, imu.sternum.orientation, accLeftArm_atShoulder_wGrav, imu.leftArm.gyr, imu.leftArm.orientation, gravityValue, fs);

        if ~isempty(shoulderDriftQ)
            imu.leftArm.orientation = correctSegmentDrift(imu.leftArm.orientation, shoulderDriftQ);
            leftArmAngVel = inverseStrapdown(fs, imu.leftArm.orientation)/pi*180;
            imu.leftArm.gyr(~isnan(leftArmAngVel(:,1)),:) = leftArmAngVel(~isnan(leftArmAngVel(:,1)),:);
        end
        
        if azimutCorrection
            accSternum_atShoulder_wGrav = translateAcc(imu.sternum, dists.sternumLeftShoulder, gravityValue, fs);
            accLeftArm_atShoulder_wGrav = translateAcc(imu.leftArm, dists.leftArmShoulder, gravityValue, fs);
            shoulderDriftQ = estimateJointDrift_AzimutOnly(accSternum_atShoulder_wGrav, imu.sternum.gyr, imu.sternum.orientation, accLeftArm_atShoulder_wGrav, imu.leftArm.gyr, imu.leftArm.orientation, gravityValue, fs, false);

            if ~isempty(shoulderDriftQ)
                imu.leftArm.orientation = correctSegmentDrift(imu.leftArm.orientation, shoulderDriftQ);
                leftArmAngVel = inverseStrapdown(fs, imu.leftArm.orientation)/pi*180;
                imu.leftArm.gyr(~isnan(leftArmAngVel(:,1)),:) = leftArmAngVel(~isnan(leftArmAngVel(:,1)),:);
            end
        end
        

        fprintf('[done]\n');
        clear accC7_atShoulder_wGrav accLeftArm_atShoulder_wGrav shoulderDriftQ leftArmAngVel
    end

    
    % Drift of right forearm
    % ----------------------
    
    if isfield(imu, 'rightArm') && isfield(imu, 'rightWrist')
        fprintf('Estimating joint drift for right forearm sensor... ');

        accRightArm_atElbow_wGrav = translateAcc(imu.rightArm, dists.rightArmElbow, gravityValue, fs);
        accRightWrist_atElbow_wGrav = translateAcc(imu.rightWrist, dists.rightWristElbow, gravityValue, fs);
        elbowDriftQ = estimateJointDrift(accRightArm_atElbow_wGrav, imu.rightArm.gyr, imu.rightArm.orientation, accRightWrist_atElbow_wGrav, imu.rightWrist.gyr, imu.rightWrist.orientation, gravityValue, fs);

        if ~isempty(elbowDriftQ)
            imu.rightWrist.orientation = correctSegmentDrift(imu.rightWrist.orientation, elbowDriftQ);
            rightWristAngVel = inverseStrapdown(fs, imu.rightWrist.orientation)/pi*180;
            imu.rightWrist.gyr(~isnan(rightWristAngVel(:,1)),:) = rightWristAngVel(~isnan(rightWristAngVel(:,1)),:);
        end
        
        if azimutCorrection
            accRightArm_atElbow_wGrav = translateAcc(imu.rightArm, dists.rightArmElbow, gravityValue, fs);
            accRightWrist_atElbow_wGrav = translateAcc(imu.rightWrist, dists.rightWristElbow, gravityValue, fs);
            elbowDriftQ = estimateJointDrift_AzimutOnly(accRightArm_atElbow_wGrav, imu.rightArm.gyr, imu.rightArm.orientation, accRightWrist_atElbow_wGrav, imu.rightWrist.gyr, imu.rightWrist.orientation, gravityValue, fs, false);

            if ~isempty(elbowDriftQ)
                imu.rightWrist.orientation = correctSegmentDrift(imu.rightWrist.orientation, elbowDriftQ);
                rightWristAngVel = inverseStrapdown(fs, imu.rightWrist.orientation)/pi*180;
                imu.rightWrist.gyr(~isnan(rightWristAngVel(:,1)),:) = rightWristAngVel(~isnan(rightWristAngVel(:,1)),:);
            end
        end

        fprintf('[done]\n');
        clear accRightArm_atElbow_wGrav accRightWrist_atElbow_wGrav elbowDriftQ rightWristAngVel
    end


    % Drift of left forearm
    % ---------------------

    if isfield(imu, 'leftArm') && isfield(imu, 'leftWrist')
        fprintf('Estimating joint drift for left forearm sensor... ');

        accLeftArm_atElbow_wGrav = translateAcc(imu.leftArm, dists.leftArmElbow, gravityValue, fs);
        accLeftWrist_atElbow_wGrav = translateAcc(imu.leftWrist, dists.leftWristElbow, gravityValue, fs);
        elbowDriftQ = estimateJointDrift(accLeftArm_atElbow_wGrav, imu.leftArm.gyr, imu.leftArm.orientation, accLeftWrist_atElbow_wGrav, imu.leftWrist.gyr, imu.leftWrist.orientation, gravityValue, fs);

        if ~isempty(elbowDriftQ)
            imu.leftWrist.orientation = correctSegmentDrift(imu.leftWrist.orientation, elbowDriftQ);
            leftWristAngVel = inverseStrapdown(fs, imu.leftWrist.orientation)/pi*180;
            imu.leftWrist.gyr(~isnan(leftWristAngVel(:,1)),:) = leftWristAngVel(~isnan(leftWristAngVel(:,1)),:);
        end
        
        if azimutCorrection
            accLeftArm_atElbow_wGrav = translateAcc(imu.leftArm, dists.leftArmElbow, gravityValue, fs);
            accLeftWrist_atElbow_wGrav = translateAcc(imu.leftWrist, dists.leftWristElbow, gravityValue, fs);
            elbowDriftQ = estimateJointDrift_AzimutOnly(accLeftArm_atElbow_wGrav, imu.leftArm.gyr, imu.leftArm.orientation, accLeftWrist_atElbow_wGrav, imu.leftWrist.gyr, imu.leftWrist.orientation, gravityValue, fs, false);

            if ~isempty(elbowDriftQ)
                imu.leftWrist.orientation = correctSegmentDrift(imu.leftWrist.orientation, elbowDriftQ);
                leftWristAngVel = inverseStrapdown(fs, imu.leftWrist.orientation)/pi*180;
                imu.leftWrist.gyr(~isnan(leftWristAngVel(:,1)),:) = leftWristAngVel(~isnan(leftWristAngVel(:,1)),:);
            end
        end
        

        fprintf('[done]\n');
        clear accLeftArm_atElbow_wGrav accLeftWrist_atElbow_wGrav elbowDriftQ leftWristAngVel
    end

end