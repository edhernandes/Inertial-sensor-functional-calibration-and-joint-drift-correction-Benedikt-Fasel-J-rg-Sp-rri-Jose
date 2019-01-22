% findJointDistances Estimates the sensor distances to common joint
%     [rProximal, rDistal] = findJointDistances(proximal, rProximalInit, distal, rDistalInit, jointName, gravityVector, gravityValue, fs)
%     estimates the sensor distances from two adacent sensors to the
%     common joint center. proximal and distal are the structures containing the
%     calibrated sensor readings and orientations in the global frame.
%     Returned are two 1-by-3 vectors rProximal and rDistal containing the
%     distances for sensors to the joint center in m. 

function [rProximal, rDistal] = findJointDistances(proximal, rProximalInit, distal, rDistalInit, jointName, gravityVector, gravityValue, fs)
    
    % Check that initial sampling rate was 500Hz (other sampling rates not
    % supported as downsampling to 50Hz is hard-coded)
    if fs~=500
        error('SegmentOrientations:jointDistance:wrongSamplingRate', 'The sampling rate must be 500Hz, other sampling rates are not supported');
    end

    % Lowpass filter at 10Hz
    [b,a] = butter(2, 1.2465*10/fs*2, 'low');
    
    proximalLow.acc = filtfilt3D(b,a, proximal.acc);
    proximalLow.gyr = filtfilt3D(b,a, proximal.gyr);
    distalLow.acc = filtfilt3D(b,a, distal.acc);
    distalLow.gyr = filtfilt3D(b,a, distal.gyr);
    
    
    % Downsample to 50Hz (take only every 10th sample)
    fs_down = 50;
    proximalLow.orientation = proximal.orientation(1:10:end-1, :);
    distalLow.orientation = distal.orientation(1:10:end-1, :);
    
    proximalLow.acc = proximalLow.acc(1:10:end,:);
    proximalLow.gyr = proximalLow.gyr(1:10:end,:);
    distalLow.acc = distalLow.acc(1:10:end,:);
    distalLow.gyr = distalLow.gyr(1:10:end,:);
    
    
    % Remove gravity
    proximalGravity = zeros(size(proximalLow.orientation,1),3);
    distalGravity = zeros(size(distalLow.orientation,1),3);
    
    for i=1:size(distalLow.orientation,1)
        proximalGravity(i,:) = gravityVector*quat2matrix(proximalLow.orientation(i,:));
        distalGravity(i,:) = gravityVector*quat2matrix(distalLow.orientation(i,:)); 
    end
    
    % Convert to m/s2
    proximalAccCorrected = (proximalLow.acc - proximalGravity)*gravityValue;
    distalAccCorrected = (distalLow.acc - distalGravity)*gravityValue;
    
    
%     % Control plots for debugging
%     plotSensorData([distalGravity proximalGravity], 'Gravity, Distal and Proximal');
%     plotSensorData([distalAccCorrected proximalAccCorrected], 'Gravity Corrected Acceleration, Distal and Proximal, Newton');


    % Compute angular acceleration
    wProximal = proximalLow.gyr/180*pi;
    aProximal = [0 0 0; (wProximal(3:end,:)-wProximal(1:end-2,:)).*(fs_down/2); 0 0 0];
    wDistal = distalLow.gyr/180*pi;
    aDistal = [0 0 0; (wDistal(3:end,:)-wDistal(1:end-2,:)).*(fs_down/2); 0 0 0];

    

    %% Estimation of Sensor Distance
   
    % Take approximately 2500 randomly sampled points
    maxNbSamples = min([2500 size(proximalLow.orientation,1)]);
    mask = rand(size(proximalLow.orientation,1),1)<maxNbSamples/size(proximalLow.orientation,1);
    
    % Remove angles that are NaN or [1 0 0 0]
    mask(isnan(proximalLow.orientation(:,1)) | proximalLow.orientation(:,1)==1 | proximalLow.orientation(:,1)==0) = false;
    mask(isnan(distalLow.orientation(:,1)) | distalLow.orientation(:,1)==1 | distalLow.orientation(:,1)==0) = false;
    
    
    % Get initial estimate of distance
    initEstimate = [rProximalInit rDistalInit];
    
    fprintf('%s: computing sensor distances (N=%d)... ', upper(jointName), sum(double(mask)));

    % Do the optimization
    [x, fVal, ~, output] = fminsearch(@(r)findAccError3DGlobalFrame(wProximal(mask,:), aProximal(mask,:), wDistal(mask,:), aDistal(mask,:), proximalAccCorrected(mask,:), distalAccCorrected(mask,:), r, proximalLow.orientation(mask,:), distalLow.orientation(mask,:)), initEstimate, optimset('Display', 'off', 'MaxIter', 300, 'TolX', 0.002, 'TolFun', 5));
    rProximal = x(1:3);
    rDistal = x(4:6);
        
    fprintf('done\n');
    
    % Show statistics
    fprintf('    Number of Iterations: %d \t Error: %1.2f m/s^2/sample\n', output.iterations, fVal/sum(double(mask)));
    fprintf('    Proximal Distance: (%1.2f, %1.2f, %1.2f)cm \t Distal Distance: (%1.2f, %1.2f, %1.2f)cm\n', rProximal(1)*100, rProximal(2)*100, rProximal(3)*100, rDistal(1)*100, rDistal(2)*100, rDistal(3)*100);
end