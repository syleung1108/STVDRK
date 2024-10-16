% This program is provided "as is" without warranty of any kind. 
% Use at your own risk. If you use this program in a publication, 
% please cite
%
% Shingyu Leung, Wai Ming Chau, Young Kyu Lee.
% SLERP-TVDRK (STVDRK) Methods for Ordinary Differential Equations on Spheres. 
% J. Sci. Comput. (arXiv:2410.10420), 2024.

function [tout, yout] = tvdrk2_proj_single(FunFcn, tspan, y0, ssize)

% Initialization

t0=tspan(1);
tfinal=tspan(2);
pm = sign(tfinal - t0);  % Which way are we computing?
if nargin < 4, ssize = (tfinal - t0)/100; end
if ssize < 0, ssize = -ssize; end
h = pm*ssize;
t = t0;
y = y0(:);

% We need to compute the number of steps.

dt = abs(tfinal - t0);
N = floor(dt/ssize) + 1;
if (N-1)*ssize < dt
  N = N + 1;
end

% Initialize the output.

tout = zeros(N,1);
tout(1) = t;
yout = zeros(N,size(y,1));
yout(1,:) = y.';
k = 1;

% The main loop
while (k < N)
  if pm*(t + h - tfinal) > 0 
    h = tfinal - t; 
    tout(k+1) = tfinal;
  else
    tout(k+1) = t0 +k*h;
  end
  k = k + 1;

  y1 = y + h*feval(FunFcn, t, y); 
  y2 = y1 + h*feval(FunFcn, t+h, y1); 
  y = (y+y2)/2;
  y=y/norm(y);
  y = y';

  t = tout(k);
  yout(k,:) = y.';
end

