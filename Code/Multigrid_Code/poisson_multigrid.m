close all
n = 1024; %Number of grid points (make a power of 2 for simplicity)
T = 3; %Number of iterations at each step
dx = 1/n; %Grid resolution
[X,Y] = meshgrid(0:dx:1);

f = ones(n+1,n+1); %Right hand side
u = zeros(n+1,n+1);

%tic; u = poisson_solver(u,f,1e-2,false);toc;

tic;
err = 1;
while err > dx^2
   [u,r] = poisson_vcycle(u,f,T,false);
   err = max(max(abs(r(2:end-1,2:end-1))))
end
toc;
