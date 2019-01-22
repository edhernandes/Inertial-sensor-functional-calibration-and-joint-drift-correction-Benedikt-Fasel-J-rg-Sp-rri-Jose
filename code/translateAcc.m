% translateAcc Translates acceleration to a virtual point
%     acc_translated = translateAcc(imuData, translation, gravityValue, fs)
%     translates the acceleration measured in imuData by translation meters
%     relative to the sensor position (i.e. from the physical sensor position).
%     The translation is performed in the local frame in which imuData is
%     represented.
%     Returned is the acceleration at the new position, in g

function acc_translated = translateAcc(imuData, translation, gravityValue, fs)
    
    % Convert angular velocity into rad/sec   
    imuAngVel = imuData.gyr ./ 180 * pi;
    
    % Convert acceleration into m/sec^2
    imuAcc = imuData.acc*gravityValue;
    
    % Compute the angular acceleration and retrieve the acceleration
    imuAngAcc = [0 0 0; (imuAngVel(3:end,:)-imuAngVel(1:end-2,:)).*(fs/2); 0 0 0];
    
    % Lowpass filter slightly (20 Hz)
    [b,a] = butter(2, 1.2465*20/fs*2, 'low');
    imuAngAcc = filtfilt3D(b,a,imuAngAcc);

    
    % Translate to new point
    % ----------------------
    
    % Initialize variables
    analysisLength = size(imuData.acc,1);
    acc_translated = zeros(analysisLength,3);

    % Sample-by-sample processing
    for i=1:analysisLength
        % Formula:
        % acc_new = acc_old + cross(angAcc, r) + cross(angVel, cross(angVel, r)) 
        % To increase speed: write out individual multiplications and sums
        
        % compute cross(angAcc(i,:), r)
        a = imuAngAcc(i,:); b = translation;
        c1 = [a(2).*b(3)-a(3).*b(2) a(3).*b(1)-a(1).*b(3) a(1).*b(2)-a(2).*b(1)];
    
        % Compute cross(angVel(i,:), r)
        a = imuAngVel(i,:); b = translation;
        c2 = [a(2).*b(3)-a(3).*b(2) a(3).*b(1)-a(1).*b(3) a(1).*b(2)-a(2).*b(1)];
        
        % Compute cross(wThigh(i,:), c2)
        a = imuAngVel(i,:); b = c2;
        c3 = [a(2).*b(3)-a(3).*b(2) a(3).*b(1)-a(1).*b(3) a(1).*b(2)-a(2).*b(1)];
    
        acc_translated(i,:) = imuAcc(i,:) + c1 + c3;
    end
   
    % Convert acceleration units back to g
    acc_translated = acc_translated./gravityValue;
end