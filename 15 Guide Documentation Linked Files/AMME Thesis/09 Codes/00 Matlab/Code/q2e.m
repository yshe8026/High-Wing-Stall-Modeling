%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert quaternions to Euler angles
% Input: 
%   - quat: [4 x n] vector of quaternions
% Output: 
%   - euler: [3 x n] vector of Euler angles in radians
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function euler =  q2e(quat)
    %  Shorthand
    q0 = quat(1, :);
    q1 = quat(2, :);
    q2 = quat(3, :);
    q3 = quat(4, :);
    
    %  Caculate Euler angles (result in radian)
    theta = atan2(q0.*q2 - q1.*q3, ((q0.^2+q1.^2-0.5).^2 + (q1.*q2+q0.*q3).^2).^0.5);
    phi = atan2(q2.*q3 + q0.*q1, q0.^2 + q3.^2 - 0.5);
    psi = atan2(q1.*q2 + q0.*q3, q0.^2 + q1.^2 - 0.5);
    
    %  Assign Euler angles to the [3 x n] vector euler and convert them to degree
    euler = [phi; theta; psi];
end