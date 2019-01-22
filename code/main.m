% main Functional calibration and joint drift correction for lower & upper limbs and trunk sensors

function main
    
    % Add the library folder to the matlab search path
    addpath('lib');
    

    % ----------------
    % Define constants
    % ----------------
    
    % Note: the script has only been tested with the below configuration.
    %     X-axis: posterior-anterior axis, pointing forwards
    %     Y-axis: inferior-superior axis, pointing upwards (parallel to
    %             gravity during perfect upright standing)
    %     Z-axis: medio-lateral axis, pointing towards the right
       
    % How gravity is measured in the global frame. Must be [0 1 0].
    gravityVector = [0 1 0];
    

    % The value of 1g in m/s^2
    %     Note: this value may change depending on your location and how you
    %     calibrated your sensors. The data of this example has been calibrated
    %     according to Ferraris (1995) where the measured acceleration during
    %     a static period is equal to 1g.
    gravityValue = 9.80604719; % at EPFL, from map.geo.admin.ch (epfl reference station)
    
    
    
    % -------------
    % Load the data
    % -------------
    
    % Segmentation information
    %     Note: This .mat file contains the sample numbers of the start and
    %     end points of each calibration movement and of an example
    %     treadmill skiing movement. It is best to choose these moments manually.
    calibInfo = load('../data/calibInfo.mat');
    
    % Inertial sensor data
    %     Note: All sensor data has been synchronized during a
    %     pre-processing step. Accelerometer offset and sensitivity as well
    %     as gyroscope offset were corrected also during a pre-processing
    %     step. These steps are independent of the functional calibration
    %     and may change between different sensor types and manufacturers.
    %     Each sensor is represented as a structure with the fields .acc
    %     (3D acceleraton, in g), .gyr (3D angular velocity, in deg/sec), 
    %     .fs (sampling frequency, in Hz. Equal for all sensors).
    imuData = load('../data/imuData.mat');

    % Extract the sampling frequency
    fs = imuData.fs;
    
    
    
    % --------------------------------
    % Calibrate the different segments
    % --------------------------------
    
    % The calibration procedure is described on protocols.io:
    % https://www.protocols.io/view/functional-calibration-for-trunk-and-lower-limb-fi-itrcem6
    % https://www.protocols.io/view/functional-calibration-for-trunk-lower-and-upper-l-jzncp5e
    % 
    % The movements are validated in the following publication:
    % B. Fasel, J. Spörri, P. Schütz, S. Lorenzetti and K. Aminian. Validation of functional 
    % calibration and strap-down joint drift correction for computing 3D joint angles of knee,
    % hip, and trunk in alpine skiing, in PLOS ONE, vol. 12, num. 7, 2017. 
    % DOI: 10.1371/journal.pone.0181446
    
    % Store the rotation matrices from the calibration 
    calibMatrix = struct();
    
    display('Start calibration');
    
    % Sacrum
    display('Calibrate sacrum');
    calibMatrix.sacrum = spineCalibration(imuData, calibInfo, 'sacrum', gravityVector, fs);
    display('Done');
    display(' ');
    
    % Sternum
    display('Calibrate sternum');
    calibMatrix.sternum = spineCalibration(imuData, calibInfo, 'sternum', gravityVector, fs);
    display('Done');
    display(' ');
    
    % Left thigh
    display('Calibrate left thigh');
    calibMatrix.leftThigh = thighCalibration(imuData, calibInfo, 'left', calibMatrix.sacrum, gravityVector, gravityValue, fs);
    display('Done');
    display(' ');
    
    % Right thigh
    display('Calibrate right thigh');
    calibMatrix.rightThigh = thighCalibration(imuData, calibInfo, 'right', calibMatrix.sacrum, gravityVector, gravityValue, fs);
    display('Done');
    display(' ');

    % Left shank
    display('Calibrate left shank');
    calibMatrix.leftShank = shankCalibration(imuData, calibInfo, 'left', gravityVector, fs);
    display('Done');
    display(' ');
    
    % Right shank
    display('Calibrate right shank');
    calibMatrix.rightShank = shankCalibration(imuData, calibInfo, 'right', gravityVector, fs);
    display('Done');
    display(' ');
    
    % Head (use same method as for trunk sensors since movement is the same)
    display('Calibrate head');
    calibMatrix.head = spineCalibration(imuData, calibInfo, 'head', gravityVector, fs);
    display('Done');
    display(' ');
    
    % Left upper arm
    display('Calibrate left arm');
    calibMatrix.leftArm = armCalibration(imuData, calibInfo, 'leftArm', gravityVector, fs);
    display('Done');
    display(' ');
    
    % Right upper arm
    display('Calibrate right arm');
    calibMatrix.rightArm = armCalibration(imuData, calibInfo, 'rightArm', gravityVector, fs);
    display('Done');
    display(' ');
    
    % Left forearm / wrist
    display('Calibrate left forearm');
    calibMatrix.leftWrist = forearmCalibration(imuData, calibInfo, 'leftWrist', 'left', gravityVector, fs);
    display('Done');
    display(' ');
    
    % Right forearm / wrist
    display('Calibrate right forearm');
    calibMatrix.rightWrist = forearmCalibration(imuData, calibInfo, 'rightWrist', 'right', gravityVector, fs);
    display('Done');
    display(' ');
    
    
    % -----------------------------------------
    % Improve lower limb functional calibration
    % -----------------------------------------
    
    
    % Correct rotation and abduction offsets: left side
    display('Correct left limb orientations');
    [calibMatrix.leftShank, calibMatrix.leftThigh] = correctKneeAngleOffset(imuData, calibInfo, 'left', calibMatrix.leftShank, calibMatrix.leftThigh, gravityVector, fs);
    display('Done');
    display(' ');
    
    % Correct rotation and abduction offsets: right side
    display('Correct right limb orientations');
    [calibMatrix.rightShank, calibMatrix.rightThigh] = correctKneeAngleOffset(imuData, calibInfo, 'right', calibMatrix.rightShank, calibMatrix.rightThigh, gravityVector, fs);
    display('Done');
    display(' ');
     
    % Make sure that left and right shank & thigh inclination (e.g. flexion) 
    % are the same during the upright period at the beginning of the squats
    display('Symmetrize lower limbs');
    calibMatrix = symmetrizeLowerlimbs(imuData, calibInfo, calibMatrix, gravityVector, fs);
    display('Done');
    display(' ');
    
    
    %% Compute segment orientations and perform drift correction
    
    % Final goal: plot 3D knee angles for left and right knee during 
    %             carpet / treadmill skiing.
    
    % The joint drift correction is explained in the following publication
    % B. Fasel, J. Spörri, P. Schütz, S. Lorenzetti and K. Aminian. An inertial 
    % sensor-based method for estimating the athlete’s relative joint center 
    % positions and center of mass kinematics in alpine ski racing, in Frontiers
    % in Physiology, 2017. 
    % DOI: 10.3389/fphys.2017.00850   
    
    % Step 1: 
    % Extract relevant data and find orientation (without joint
    % drift correction)
    
    display('Compute all segment orientations during skiing');
        
    % Apply calibration for the skiing movement
    for currentSegment = fieldnames(imuData)'
        if ~strcmp(currentSegment{1}, 'fs')
            imuSkiing.(currentSegment{1}).acc = imuData.(currentSegment{1}).acc(calibInfo.skiing_start:calibInfo.skiing_stop,:) * calibMatrix.(currentSegment{1})';
            imuSkiing.(currentSegment{1}).gyr = imuData.(currentSegment{1}).gyr(calibInfo.skiing_start:calibInfo.skiing_stop,:) * calibMatrix.(currentSegment{1})';
        end
    end
    
    % Strapdown with static drift correction for each segment
    for currentSegment = fieldnames(imuData)'
        if ~strcmp(currentSegment{1}, 'fs')
            imuSkiing.(currentSegment{1}).orientation = strapdownStaticDriftCorrection(imuSkiing.(currentSegment{1}), gravityVector, fs);
        end
    end
    
    display('Done');
    
    
    % Step 2:
    % Correct for initial azimuth offset (especially for upper limbs)
    
    display('Correct initial azimuth offsets');
    imuSkiing = correctForAzimutOffset(imuSkiing, gravityValue, fs);
    display('Done');
    
    % Step 3:
    % Joint drift correction 

    % Estimate distances between sensors and joint centers
    fprintf('Finding sensor distance to joint centres... \n');
    dists = computeDistances(imuSkiing, gravityVector, gravityValue, fs);
    fprintf('\nFinished computing distances\n\n');
        
    % Joint drift correction
    display('Joint drift correction');
    imuSkiing = jointDriftCorrection(imuSkiing, dists, gravityValue, fs);
    display('Done');
    
    
    % Step 4:
    % Compute knee angles and plot data
    
    % Compute the knee angles
    display('Computing the knee angles');
    leftKneeAngles = computeKneeAngles(imuSkiing.leftShank.orientation, imuSkiing.leftThigh.orientation, 'left');
    rightKneeAngles = computeKneeAngles(imuSkiing.rightShank.orientation, imuSkiing.rightThigh.orientation, 'right');
    display('Done');
    display(' ');
    
    % Plot acceleration and angular velocity for left leg and knee angles
    display('Preparing figures for left leg inertial data and left knee angle');
    
    % Timestamp to get x-axis in seconds
    timestamps = (1:size(imuSkiing.leftThigh.acc,1))'./fs;
    
    % Figure for left thigh inertial data
    fLeftThigh = figure('name', 'Left Thigh Inertial Data');
    sA(1) = subplot(2,1,1);
    plot(timestamps, imuSkiing.leftThigh.acc);
    title('AF Left Thigh Acceleration During Skiing');
    ylabel('Acceleration, g');
    xlabel('Time, sec');
    
    sA(2) = subplot(2,1,2);
    plot(timestamps, imuSkiing.leftThigh.gyr);
    title('AF Left Thigh Angular Velocity During Skiing');
    ylabel('Angular Velocity, deg/sec');
    xlabel('Time, sec');
    legend('Post-Ant', 'Inf-Sup', 'Med-Lat');
    linkaxes(sA, 'x');
    xlim([0 30]);
    
    % Figure for left shank inertial data
    fLeftShank = figure('name', 'Left Shank Inertial Data');
    sB(1) = subplot(2,1,1);
    plot(timestamps, imuSkiing.leftShank.acc);
    title('AF Left Shank Acceleration During Skiing');
    ylabel('Acceleration, g');
    xlabel('Time, sec');
    
    sB(2) = subplot(2,1,2);
    plot(timestamps, imuSkiing.leftShank.gyr);
    title('AF Left Shank Angular Velocity During Skiing');
    ylabel('Angular Velocity, deg/sec');
    xlabel('Time, sec');
    legend('Post-Ant', 'Inf-Sup', 'Med-Lat');
    linkaxes(sB, 'x');
    xlim([0 30]);
    
    % Figure for left knee angles
    fKneeAngles = figure('name', 'Left and Right Knee Angles'); 
    hold on;
    plot(timestamps, leftKneeAngles, ':', 'lineWidth', 1);
    
    % Restart color order
    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(timestamps, rightKneeAngles, '-', 'lineWidth', 1);
    title('Knee Angles During Skiing');
    ylabel('Angle, deg');
    xlabel('Time, sec');
    legend('Left Flex.', 'Left Abd.', 'Left Rot.', 'Right Flex.', 'Right Abd.', 'Right Rot.');
    xlim([0 30]);
    
    display('Done');
    display(' ');
    
    
    % Store figures
    display('Saving figures in the results folder');
    saveas(fLeftShank, '../results/leftShank.png');
    saveas(fLeftThigh, '../results/leftThigh.png');
    saveas(fKneeAngles, '../results/kneeAngles.png');
    display('Done');
    display(' ');
    
    display('============================================');
    display('            Processing finished             ');
    display('============================================');
end