function dy=functiondy(t,y)

v1=[1 -1 1]';
v2=[1 -1 -1]';
v3=[-2 1 0]';
v4=[-1 -1 1]';
v1=v1/norm(v1);
v2=v2/norm(v2);
v3=v3/norm(v3);
v4=v4/norm(v4);

dy=cross(v1,y)/2/(1-dot(v1,y))+cross(v2,y)/2/(1-dot(v2,y))+...
    cross(v3,y)/2/(1-dot(v3,y))+cross(v4,y)/2/(1-dot(v4,y));
dy=dy/2/pi;

return