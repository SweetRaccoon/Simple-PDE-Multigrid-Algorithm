close all
n = 128; %Number of grid points (make a power of 2 for simplicity)
T = 3; %Number of iterations at each step

f = ones(n+1,n+1); %Right hand side
u = zeros(n+1,n+1);

tic;poisson_solver(u,f,1e-2,false);toc;

tic;
err = 1;
while err > 1e-2

   %Iterate at fine scale
   [u,r] = poisson_iterations(u,f,T,false); 
   err = max(max(abs(r(2:end-1,2:end-1))));

   %Subsample residual
   rs = r(1:2:end,1:2:end);

   %Iterate at coarse scale
   vs = poisson_iterations(zeros(size(rs)),rs,T*20,false); 

   %This is where you insert 1 additional coarsening
   v = interp2(vs,1,'linear');

   %add residual back to u
   u = u + v; 

end
toc;
