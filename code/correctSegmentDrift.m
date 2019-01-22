% correctSegmentDrift Removes previously computed drift from a segment
%    orientationCorrected = correctSegmentDrift(orientation, drift) removes
%    the drift specified in the N-by-4 quaternion drift from the segment
%    orientation provided in the N-by-4 quaternion orientation.
%    Returned is the N-by-4 segment orientation quaternion without drift.

function orientationCorrected = correctSegmentDrift(orientation, drift)

    orientationCorrected = zeros(size(orientation));
           
    % Correct orientation
    for i=1:size(drift,1)
        orientationCorrected(i,:) = quat_multiply(drift(i,:), orientation(i,:));
    end
end