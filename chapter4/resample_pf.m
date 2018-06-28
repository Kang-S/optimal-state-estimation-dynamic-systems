function [x_resamp,w_resamp,index]=resample_pf(x_particle,w_particle);

% Get Length of Particles
n=length(x_particle);

% Cumulative Sum of Particles
w_particle=w_particle(:);
c=cumsum(w_particle);

% Compute u Vector
u=zeros(n,1);
u(1)=inv(n)*rand(1);
u(2:n)=u(1)+inv(n)*(1:n-1)';

% Pre-allocate Index
index=zeros(n,1);

% Compute Index for Resampling 
i=1;
for j=1:n
    while u(j)>c(i)
        i=i+1;
    end
    index(j)=i;
end

% Resampled Data
x_resamp=x_particle(index,:);
w_resamp=inv(n)*ones(n,1);
