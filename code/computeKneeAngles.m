% computeKneeAngles Computes the 3D knee joint angles
%     angles3D = computeKneeAngles(qShank, qThigh, side) computes the 3D 
%     knee joint angles where the orientation of the shank is given in
%     qShank, the orientation of the thigh in qThigh and the attribute
%     side is either 'left' or 'right', indicating the side. qShank and
%     qThigh are N-by-4 matrices where each line represents a sample.
%     angles3D is a N-by-3 vector of the knee joint angles in degrees where
%     the first column is flexion, second colum abduction and third column
%     internal-external rotation

function angles3D = computeKneeAngles(qShank, qThigh, side)

    angles3D = zeros(size(qThigh,1),3);

    for ind=1:size(qThigh,1)    

        % Get orientation matrix of thigh and shank
        thighO = quat2matrix(qThigh(ind,:));
        shankO = quat2matrix(qShank(ind,:));


        % Grood & Suntay 1983
        % -------------------

        % Coordinate transform from my system into Grood&Suntay
        % shankGrood = [shankO(:,3) shankO(:,1) shankO(:,2)];
        % thighGrood = [thighO(:,3) thighO(:,1) thighO(:,2)];

        % Axes of rotation that we need (includes the coordinate system
        % transform)
        I = thighO(:,3);
        J = thighO(:,1);
        K = thighO(:,2);
        i = shankO(:,3);
        j = shankO(:,1);
        k = shankO(:,2);

        e1 = thighO(:,3);
        e3 = shankO(:,2);

        % Compute the floating axis e2 based on e3 and e1
        e2 = cross(e3, e1);
        e2 = e2 ./ norm(e2);

        % Compute the angles
        flexion = real(acos(J'*e2));
        extRot = real(acos(j'*e2));
        
        if strcmpi(side, 'left')
            add = pi/2 - acos(I'*k); % For left knee
        else
            add = -pi/2 + acos(I'*k); % For right knee
        end
        
        
        % Determine the sign for flexion and extRot
        flexSign = sign(asin(-e2'*K)); flexion = flexion * flexSign;
        extSign = sign(asin(e2'*i)); extRot = extRot * extSign;

        angles3D(ind,:) = [flexion add extRot];
    end

    angles3D = angles3D / pi*180;
end