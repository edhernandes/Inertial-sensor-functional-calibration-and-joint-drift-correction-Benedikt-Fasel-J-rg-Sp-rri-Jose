% shankCalibration Calibrate the shank IMU
%     shankRotMatrix = shankCalibration(imuData, calibInfo, side, gravityVector, fs)
%     computes the rotation matrix required for the functional calibration
%     of the shank.

function shankRotMatrix = shankCalibration(imuData, calibInfo, side, gravityVector, fs)
    
    
    % Lowpass filter parameters
    [bLow,aLow] = butter(2, 1.2465*1/fs*2, 'low');        % cutoff freq. 1 Hz 
    
    % Load the data
    if strcmpi(side, 'left')
        abductionShank_LF.acc = imuData.leftShank.acc(calibInfo.left_abduction_start:calibInfo.left_abduction_stop,:);
        abductionShank_LF.gyr = imuData.leftShank.gyr(calibInfo.left_abduction_start:calibInfo.left_abduction_stop,:);
    else
        abductionShank_LF.acc = imuData.rightShank.acc(calibInfo.right_abduction_start:calibInfo.right_abduction_stop,:);
        abductionShank_LF.gyr = imuData.rightShank.gyr(calibInfo.right_abduction_start:calibInfo.right_abduction_stop,:);
    end
    
    
    %% Find anterior-posterior alignment rotation matrix
    
    % Hypothesis:
    % Principal axis of rotation during the abduction movement corresponds to the
    % posterior-anterior X-axis [1 0 0].
    
    % Note: For more details please also refer to the comments in the 
    % function spineCalibration
    
    
    % Extract angular velocity
    shankAngVel = filtfilt3D(bLow, aLow, abductionShank_LF.gyr);
    
    % Compute the norm
    shankAngVelNorm = sqrt(sum(shankAngVel.^2,2));
    
    % Only consider moments with some rotation to avoid considering also movement noise
    angVelThres = 30; 
    takeValues = shankAngVelNorm>angVelThres;
    
	% Ignore this calibration if no movement present
	if sum(double(takeValues))<10
		shankRotMatrix = eye(3);
		warning('functionalCalibration:shank:noMovement', 'Abduction movement not present for calibration. Ignoring...');
	else

		shankCoeff = pca(shankAngVel(takeValues, :));

		% Get first rotation axis
		[shankRotationAxis, shankRotationAngle] = vec2helic(shankCoeff(:,1), [1 0 0]);
		
		shankRotMatrix = quat2matrix(helic2quat(shankRotationAxis, shankRotationAngle));
		
		shankAngVelCalib = shankAngVel * shankRotMatrix';
		
	  
		% Check orientation of the x axis (must point forwards)
		% Hypothesis: first peak of angular velocity must be
		% <0 for the right side
		% >0 for the left side
			
		% Massively lowpass filter because we only want global oscillation pattern
		abductionLow = filtfilt(bLow,aLow, shankAngVelCalib(:,1));
		
		% Find all the angular velocity peaks
		[~, locs] = findpeaks(abs(abductionLow), 'minpeakheight', std(abductionLow)/2);
		        
        if (abductionLow(locs(1))>0 && strcmp(side, 'right')) || (abductionLow(locs(1))<0 && strcmp(side, 'left'))
			shankRotMatrix = yRotation(pi) * shankRotMatrix;
        end 
	end
    
    
    
    %% Find medio-lateral/vertical axis
    % Hypothesis: Gravity must be 0 along medio-lateral axis and maximum along vertical axis
    
    % Get static acceleration (suppose no movement during the first 0.5 seconds)
    staticShank = median(abductionShank_LF.acc(1:round(0.5*fs),:)) * shankRotMatrix';
       
    % Orient gravity on shank
    [rotationAxis, rotationAngle] = vec2helic([0 staticShank(2) staticShank(3)], gravityVector);
    gravityRotMatrix = quat2matrix(helic2quat(rotationAxis, rotationAngle));
    
    shankRotMatrix = gravityRotMatrix * shankRotMatrix;   
end