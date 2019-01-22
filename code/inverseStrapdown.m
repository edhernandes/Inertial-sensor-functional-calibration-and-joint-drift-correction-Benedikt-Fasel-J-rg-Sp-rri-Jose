% inverseStrapdown Computes the inverse strapdown operation
%    angularVelocity = inverseStrapdown(fs, orientation) performes the
%    inverse strapdown in order to find the angular velocities that were
%    needed to obtain orientation. fs is the sampling frequency in Hz
%    and the angular velocity is radians per second

function  angularVelocity = inverseStrapdown(fs, orientation)
    nbSamples = size(orientation,1);
    
    angularVelocity = zeros(nbSamples-1, 3);
    
    for i=1:nbSamples-1
        dq = quat_multiply(orientation(i+1,:), quat_inv(orientation(i,:)));
        [rotationAxis, rotationAngle] = quat2helic(dq);
        
        angVelGlobalFrame = rotationAxis .* rotationAngle;
        angularVelocity(i,:) = angVelGlobalFrame * quat2matrix(quat_inv(orientation(i,:)))';
    end
    
    angularVelocity = angularVelocity .* fs;
end