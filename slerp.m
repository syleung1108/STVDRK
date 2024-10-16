% This program is provided "as is" without warranty of any kind. 
% Use at your own risk. If you use this program in a publication, 
% please cite
%
% Shingyu Leung, Wai Ming Chau, Young Kyu Lee.
% SLERP-TVDRK (STVDRK) Methods for Ordinary Differential Equations on Spheres. 
% J. Sci. Comput. (arXiv:2410.10420), 2024.

function qm = slerp(qi, qn, t)

qi = qi(:);
qn = qn(:);

c = dot(qi, qn);
theta = acos(c);

if t == 0
    qm = qi;
elseif t == 1
    qm = qn;
else
    qm = qi * (sin((1 - t) * theta) / sin(theta)) + qn * (sin(t * theta) / sin(theta));
end
end