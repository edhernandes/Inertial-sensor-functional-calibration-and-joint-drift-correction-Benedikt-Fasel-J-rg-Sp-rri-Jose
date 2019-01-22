% strapdownStaticDriftCorrection Performs strapdown integration on angular velocities
%    orientation = strapdownStaticDriftCorrection(imuData, gravityVector, fs)
%    returns the drift corrected orientation of the IMU in the global frame 
%    by 3D integration of the angular velocity. Each motionless instant is
%    used to correct the drift. The drift for the azimut is NOT corrected.
%    For correct operation, the measurement should start with a motionless
%    period so that the initial conditions of the strapdown can be
%    estimated reliably.
% 
%    Choice of initial conditions:
%    - Segment inclination is estimated based on the measured gravity 
%    - Segment azimuth / heading is set to 0
%
%    The implemented drift correction is based on the work (adapted for
%    numerical stability) from:
%    - Favre, J., Jolles, B. M., Siegrist, O., & Aminian, K. (2006). 
%        Quaternion-based fusion of gyroscopes and accelerometers to improve 
%        3D angle measurement. Electronics Letters, 42(11), 612.
%    - Sabatini, A. M. (2005). Quaternion-based strap-down integration method 
%        for applications of inertial sensing to gait analysis. Medical & 
%        Biological Engineering & Computing, 43(1), 94–101.

function orientation = strapdownStaticDriftCorrection(imuData, gravityVector, fs)

    % Number of samples required for 0.01 seconds
    oneHundreth = round(0.01*fs);
    
    % Convert the angular velocity into rad/sec
    imuData.gyr = imuData.gyr / 180 * pi;
    
    % Estimate all motionless instants based on the constraints that
    % gravity norm is close to 1 and not changing fast + low ang velocity
    [accLowpass(:,1), accStd(:,1)] = movAvg(imuData.acc(:,1), oneHundreth);
    [accLowpass(:,2), accStd(:,2)] = movAvg(imuData.acc(:,2), oneHundreth);
    [accLowpass(:,3), accStd(:,3)] = movAvg(imuData.acc(:,3), oneHundreth);
    [angLowpass(:,1), angStd(:,1)] = movAvg(imuData.gyr(:,1), oneHundreth);
    [angLowpass(:,2), angStd(:,2)] = movAvg(imuData.gyr(:,2), oneHundreth);
    [angLowpass(:,3), angStd(:,3)] = movAvg(imuData.gyr(:,3), oneHundreth);
    
    % Fill in first part of zeros with original values
    angLowpass(1:round(oneHundreth/2),:) = imuData.gyr(1:round(oneHundreth/2),:);
    accLowpass(1:round(oneHundreth/2),:) = imuData.acc(1:round(oneHundreth/2),:);

    % Compute norm of lowpass filtered acceleration and angular velocity
    accNorm = sqrt(sum(accLowpass.^2, 2));
    angNorm = sqrt(sum(angLowpass.^2, 2));
    
    % Smooth the moving standard deviations (cutoff freq. 10Hz)
    [b,a] = butter(2, 1.2465*10/fs*2, 'low');
    accStd = filtfilt(b,a,sqrt(sum(accStd.^2, 2)));
    angStd = filtfilt(b,a, sqrt(sum(angStd.^2, 2)));
    
    % Due to calibration inaccuracies or a change of Earth's gravity field
    % with respect to its initial calibration it may happen that acceleration 
    % norm at rest / during motionless is not exactly 1g. Thus, find the 
    % "true" measured gravity at rest. Suppose that rest periods are
    % moments with very little change in acceleration over time --> very
    % low moving standard deviation   
    accMag = median(accNorm(accStd<0.05));
    if isnan(accMag)
        accMag = 1;
    end
    
    % Determine all motionless moments. For a sample to be motionless
    % multiple conditions need to be fulfilled. The proposed thresholds
    % below are working for a large number of situations. Reducing the
    % thresholds results in fewer motionless moments but a more precise
    % estimation of true orientation for the static drift correction.
    motionless = abs(accNorm-accMag)<0.06 & angNorm<(15/180*pi) & accStd<0.06 & angStd<0.04;

    % Smooth the detected motionless periods
    motionless = movAvg(double(motionless), round(0.15*fs));
    motionless(motionless<0.5) = 0;
    motionless(motionless>0) = 1;
    motionless = logical(motionless);
    motionless = [motionless; false];
    
    % Find the very first instant of motionless
    firstMotionless = find(motionless, 1, 'first');
    
    % If no motionless present, define first sample after 0.025 second as
    % the first instant of motionless. Initial conditions will be wrongly
    % estimated but this way allows still to obtain a segment orientation.
    if isempty(firstMotionless)
        firstMotionless = max(11, ceil(0.025*fs));
        warning('strapdown:motionless:noMotionlessFound', 'No motionless period was found in the measurement. The computed orientation will most likely be wrong.');
    end
        
    % If there is a lot of movement during the first 1.5 seconds of
    % recording
    if firstMotionless>1.5*fs
        firstMotionless = round(0.2*fs);
        warning('strapdown:motionless:motionlessTooLate', 'No motionless period was found in the beginning of the measurement. Please try to resegment your measurement to obtain accurate results.');
    end
    

    % Initialise orientation matrix to store the orientation quaternions
    orientation = zeros(size(imuData.gyr,1), 4);  
    
    % Correct angular velocity by sampling frequency for integration
    angVel = imuData.gyr ./ fs; 
    
    % Find initial orientation at time instant motionlessStart(1)
    % Compute initial orientation at strapdown start. Take average
    % acceleration during 21 samples (independent of sampling frequency)
    earthAcc_imu = mean(imuData.acc(max(1,firstMotionless-10):firstMotionless+10,:));
    
    % Normalize to one
    earthAcc_imu = earthAcc_imu/norm(earthAcc_imu);

    % Find the required rotation to align the segment with gravity
    % --> this is the initial orientation
    [rotationAxis, rotationAngle] = vec2helic(earthAcc_imu, gravityVector);
    
    % Numerical stability: if segment was already aligned at the beginning 
    % reset rotation axis to [1 0 0].
    if rotationAngle == 0;
        rotationAxis = [1 0 0];
    end
    
    % Transform into quaternion representation
    qInit = helic2quat(rotationAxis, rotationAngle);
    
    
    % Reset azimut to 0
    initOrientation = quat2matrix(qInit);
    initOrientation = initOrientation(:,1)';
    initOrientation(2) = 0; % Project onto x-z plane
    initOrientation = initOrientation ./ norm(initOrientation);
    azimutRot = yRotation(0);
    [rotationAxis, rotationAngle] = vec2helic(initOrientation, azimutRot(:,1)'); % Find azimut angle
    
    % Numerical stability: if segment was already aligned at the beginning 
    % reset rotation axis to [1 0 0].
    if rotationAngle == 0;
        rotationAxis = [1 0 0];
    end
    
    % Azimuth "offset"
    qAzimutCorrection = helic2quat(rotationAxis, rotationAngle);

    % Store final initial orientation
    orientation(firstMotionless,:) = quat_normalize(quat_multiply(qAzimutCorrection, qInit));
    
    % Set orientation previous to this first motionless constant
    % Setting it constant avoids a discontinuity at the first motionless
    % instant. If one would like to be very correct, then a time-reversed
    % strapdown could be used to estimate the real orientation prior to the
    % first motionless.
    orientation(1:firstMotionless-1,:) = repmat(orientation(firstMotionless,:), firstMotionless-1, 1);
    
    % Store the index of the last motionless 
    lastCorrectedMotionlessInd = 1; 
    
    % Strapdown
    for i=firstMotionless:size(orientation,1)-1
        
        % Find axis-angle representation of angular displacement
        % based on the measured angular velocity
        angVelGlobal = angVel(i,:)*quat2matrix(orientation(i,:))';
        angle = norm(angVelGlobal);
        dq = helic2quat(angVelGlobal, angle);

        % Update
        orientation(i+1, :) = quat_normalize(quat_multiply(dq, orientation(i, :)));
        
        % If we are at motionless moment, correct drift
        if motionless(i+1)
            
            % Get measured acceleration in GF
            measuredGravityQ = quat_multiply(orientation(i+1,:), quat_multiply([0 imuData.acc(i+1,:)], quat_inv(orientation(i+1,:))));
            
            % Find the drift
            [axis, angle] = vec2helic(measuredGravityQ(2:4), gravityVector);
            
            % For numerical stability of no drift was found
            if angle==0
                axis = [1 0 0];
            end
            
            % Hypothesis: drift is constant over time (between initialOrientation
            % and finalOrientation)
            correctionStart = lastCorrectedMotionlessInd;
            correctionStop = i+1;
            correctionLength = correctionStop-correctionStart;
            
            driftPerSecond = abs(angle)/pi*180/(correctionLength/fs);
           
            % If drift is larger than 2.5 degrees per second there must be an
            % error in its estimation (noise) -> do not correct
            if driftPerSecond<2.5
                angleCorrection = angle / correctionLength;

                % Correct orientation
                for ii=correctionStart+1:i+1
                    driftCorrection = helic2quat(axis, angleCorrection*(ii-correctionStart));
                    orientation(ii,:) = quat_multiply(driftCorrection, orientation(ii,:));
                end
                
                lastCorrectedMotionlessInd = i+1;
            end
        end
    end
end