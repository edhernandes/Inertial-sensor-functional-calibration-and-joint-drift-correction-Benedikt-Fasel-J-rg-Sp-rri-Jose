% findAccError3DGlobalFrame Computes sum of norm of translated acceleration difference
%    totalErr = findAccError3DGlobalFrame(wProximal, aProximal, wDistal, aDistal, proximalAcc, distalAcc, r, orientationProximal, orientationDistal)
%    returns the sum of norm of the acceleration difference at distance r
%    from the proximal and distal sensors. This function is minimized for
%    finding the sensor-joint center distances. 
%    r is a 1-by-6 vector where the first three elements refer to the 3D
%    distance vector for the proximal sensor and the last three element
%    refer to the 3D distance vector for the distal sensor.

function totalErr = findAccError3DGlobalFrame(wProximal, aProximal, wDistal, aDistal, proximalAcc, distalAcc, r, orientationProximal, orientationDistal)
    rProximal = r(1:3);
    rDistal = r(4:6);

    totalErr = 0;
    for i=1:size(wProximal,1)
        
        % Computations for thigh
        % ----------------------
        
        % compute cross(aThigh(i,:), rThigh)
        a = aProximal(i,:); b = rProximal;
        c1 = [a(2).*b(3)-a(3).*b(2) a(3).*b(1)-a(1).*b(3) a(1).*b(2)-a(2).*b(1)];
    
        % Compute cross(wThigh(i,:), rThigh)
        a = wProximal(i,:); b = rProximal;
        c2 = [a(2).*b(3)-a(3).*b(2) a(3).*b(1)-a(1).*b(3) a(1).*b(2)-a(2).*b(1)];
        
        % Compute cross(wThigh(i,:), c2)
        a = wProximal(i,:); b = c2;
        c3 = [a(2).*b(3)-a(3).*b(2) a(3).*b(1)-a(1).*b(3) a(1).*b(2)-a(2).*b(1)];
    
        accThigh_atKnee = proximalAcc(i,:) + c1 + c3;
        
        
        % Computations for shank
        % ----------------------
        
        % compute cross(aShank(i,:), rShank)
        a = aDistal(i,:); b = rDistal;
        c1 = [a(2).*b(3)-a(3).*b(2) a(3).*b(1)-a(1).*b(3) a(1).*b(2)-a(2).*b(1)];
    
        % Compute cross(wShank(i,:), rShank)
        a = wDistal(i,:); b = rDistal;
        c2 = [a(2).*b(3)-a(3).*b(2) a(3).*b(1)-a(1).*b(3) a(1).*b(2)-a(2).*b(1)];
        
        % Compute cross(wShank(i,:), c2)
        a = wDistal(i,:); b = c2;
        c3 = [a(2).*b(3)-a(3).*b(2) a(3).*b(1)-a(1).*b(3) a(1).*b(2)-a(2).*b(1)];
    
        accShank_atKnee = distalAcc(i,:) + c1 + c3;
        
        

        % Change into global frame
        accThigh_atKnee_globalFrame = accThigh_atKnee * (quat2matrix(orientationProximal(i,:))');
        accShank_atKnee_globalFrame = accShank_atKnee * (quat2matrix(orientationDistal(i,:))');
        
        difference = accShank_atKnee_globalFrame-accThigh_atKnee_globalFrame;
        totalErr = totalErr + sqrt(difference(1)*difference(1) + difference(2)*difference(2) + difference(3)*difference(3));
    end
    %fprintf('Total error:\t%1.2f\n', totalErr);
end