% This program is provided "as is" without warranty of any kind. 
% Use at your own risk. If you use this program in a publication, 
% please cite
%
% Shingyu Leung, Wai Ming Chau, Young Kyu Lee.
% SLERP-TVDRK (STVDRK) Methods for Ordinary Differential Equations on Spheres. 
% J. Sci. Comput. (arXiv:2410.10420), 2024.

clear all
close all

tf=2; 
y0=[1,0,0];
t0=0;
max_I=12;

dt_error=zeros(max_I,1);
error=zeros(max_I,13);
error_length=zeros(max_I,13);
casecount=0;
load Exact_Reference.mat

casecount=casecount+1;
for i=1:max_I
    m=2^(i+1); dt=(tf-t0)/(m-1);
    dt_error(i)=dt;    
    [t,tmp]=tvdrk2(@functiondy,[t0 tf],y0,dt);
    error(i,casecount)=norm(tmp(end,:)-X_exact);
    error_length(i,casecount)=abs(norm(tmp(end,:))-1);
end

casecount=casecount+1;
for i=1:max_I
    m=2^(i+1); dt=(tf-t0)/(m-1);
    [t,tmp]=tvdrk3(@functiondy,[t0 tf],y0,dt);
    error(i,casecount)=norm(tmp(end,:)-X_exact);
    error_length(i,casecount)=abs(norm(tmp(end,:))-1);
end

casecount=casecount+1;
for i=1:max_I
    m=2^(i+1); dt=(tf-t0)/(m-1);
    [t,tmp]=tvdrk2_proj(@functiondy,[t0 tf],y0,dt);
    error(i,casecount)=norm(tmp(end,:)-X_exact);
    error_length(i,casecount)=abs(norm(tmp(end,:))-1);
end

casecount=casecount+1;
for i=1:max_I
    m=2^(i+1); dt=(tf-t0)/(m-1);
    [t,tmp]=tvdrk2_proj_single(@functiondy,[t0 tf],y0,dt);
    error(i,casecount)=norm(tmp(end,:)-X_exact);
    error_length(i,casecount)=abs(norm(tmp(end,:))-1);
end

casecount=casecount+1;
for i=1:max_I
    m=2^(i+1); dt=(tf-t0)/(m-1);
    [t,tmp]=tvdrk3_proj(@functiondy,[t0 tf],y0,dt);
    error(i,casecount)=norm(tmp(end,:)-X_exact);
    error_length(i,casecount)=abs(norm(tmp(end,:))-1);
end

casecount=casecount+1;
for i=1:max_I
    m=2^(i+1); dt=(tf-t0)/(m-1);
    [t,tmp]=tvdrk3_proj_single(@functiondy,[t0 tf],y0,dt);
    error(i,casecount)=norm(tmp(end,:)-X_exact);
    error_length(i,casecount)=abs(norm(tmp(end,:))-1);
end

casecount=casecount+1;
for i=1:max_I
    m=2^(i+1); dt=(tf-t0)/(m-1);
    [t,tmp]=slerp_fe(@functiondy,[t0 tf],y0,dt);
    error(i,casecount)=norm(tmp(end,:)-X_exact);
    error_length(i,casecount)=abs(norm(tmp(end,:))-1);
end

casecount=casecount+1;
for i=1:max_I
    m=2^(i+1); dt=(tf-t0)/(m-1);
    [t,tmp]=slerp_tvdrk2(@functiondy,[t0 tf],y0,dt);
    error(i,casecount)=norm(tmp(end,:)-X_exact);
    error_length(i,casecount)=abs(norm(tmp(end,:))-1);
end

casecount=casecount+1;
for i=1:max_I
    m=2^(i+1); dt=(tf-t0)/(m-1);
    [t,tmp]=slerp_tvdrk3(@functiondy,[t0 tf],y0,dt);
    error(i,casecount)=norm(tmp(end,:)-X_exact);
    error_length(i,casecount)=abs(norm(tmp(end,:))-1);
end

%%%%%%
%   plot graph
%%%%%%

range=1:max_I;
p=zeros(9,2);
for i=1:9
    p(i,:)=polyfit(log(dt_error(range)),log(error(range,i)),1);
    q=exp(p(i,2))*dt_error(range).^p(i,1);
    loglog(dt_error(range),q,'--')
    hold on
end
loglog(dt_error,error,'o')
legend( ...
    strcat('TVDRK2: ',num2str(p(1,1))), ...
    strcat('TVDRK3: ',num2str(p(2,1))),...
    strcat('PTVDRK2^\prime: ',num2str(p(3,1))),...
    strcat('PTVDRK2: ',num2str(p(4,1))),...
    strcat('PTVDRK3^\prime: ',num2str(p(5,1))),...
    strcat('PTVDRK3: ',num2str(p(6,1))),...
    strcat('Spherical FE: ',num2str(p(7,1))),...
    strcat('STVDRK2: ',num2str(p(8,1))), ...
    strcat('STVDRK3: ',num2str(p(9,1))), ...
    'Location','northeast')
axis([1e-4 1e5 1e-16 1e-1])
print -djpeg Convergence1.jpg

figure
for i=[1 3 4 5 7 8]
    p(i,:)=polyfit(log(dt_error(range)),log(error(range,i)),1);
    q=exp(p(i,2))*dt_error(range).^p(i,1);
    loglog(dt_error(range),q)
    hold on
end
loglog(dt_error,error(:,[1 3 4 5 7 8]),'o')
legend( ...
    strcat('TVDRK2: ',num2str(p(1,1))), ...
    strcat('PTVDRK2^\prime: ',num2str(p(3,1))),...
    strcat('PTVDRK2: ',num2str(p(4,1))),...
    strcat('PTVDRK3^\prime: ',num2str(p(5,1))),...
    strcat('SFE: ',num2str(p(8,1))),...
    strcat('STVDRK2: ',num2str(p(9,1))), ...
    'Location','southeast')
axis([1e-4 1e0 1e-16 1e-1])
print -djpeg Convergence2.jpg

figure
for i=[2 6 9]
    p(i,:)=polyfit(log(dt_error(range)),log(error(range,i)),1);
    q=exp(p(i,2))*dt_error(range).^p(i,1);
    loglog(dt_error(range),q)
    hold on
end
loglog(dt_error,error(:,[2 6 9]),'o')
legend( ...
    strcat('TVDRK3: ',num2str(p(2,1))),...
    strcat('PTVDRK3: ',num2str(p(6,1))),...
    strcat('STVDRK3: ',num2str(p(9,1))), ...
    'Location','northwest')
axis([1e-4 1e0 1e-16 1e-1])
print -djpeg Convergence3.jpg

