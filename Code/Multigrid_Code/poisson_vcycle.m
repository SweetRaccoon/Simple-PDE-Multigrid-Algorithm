function [u,r] = poisson_vcycle(u,f,T,pl)
%Poisson vcycle
%f should be square nxn matrix
%Runs for T iterations 

   n = max(size(f)); %number of grid points
   dx = 1/(n-1); %Grid resolution
   if pl
      [X,Y] = meshgrid(0:dx:1);
   end
   
   %Time step
   dt = 0.9*dx^2/4;
   
   if n < 5
      u = poisson_solver(u,f,1e-1,pl);
   else
      for i = 1:T

         %Compute partial derivatives
         uxx = (u([2:n,n],:) - 2*u + u([1,1:n-1],:))/dx^2;
         uyy = (u(:,[2:n,n]) - 2*u + u(:,[1,1:n-1]))/dx^2;
         Laplacian = uxx + uyy;
         r = Laplacian + f;
         
         %Update
         u = u + dt*(Laplacian + f);

         %Boundary condition
         u(1,:) = 0; u(n,:) = 0; u(:,1) = 0; u(:,n) = 0;

         if pl
            surf(X,Y,u);
            drawnow
            %pause
         end
      end

      %Subsample residual and f
      rs = r(1:2:end,1:2:end);
      
      %Recurse
      vs = poisson_vcycle(zeros(size(rs)),rs,T,pl);

      %Interpolate back to this grid
      v = interp2(vs,1,'linear');

      %add residual back to u
      u = u + v; 

      %Iterate T more times
      for i = 1:T

         %Compute partial derivatives
         uxx = (u([2:n,n],:) - 2*u + u([1,1:n-1],:))/dx^2;
         uyy = (u(:,[2:n,n]) - 2*u + u(:,[1,1:n-1]))/dx^2;
         Laplacian = uxx + uyy;
         r = Laplacian + f;
         
         %Update
         u = u + dt*(Laplacian + f);

         %Boundary condition
         u(1,:) = 0; u(n,:) = 0; u(:,1) = 0; u(:,n) = 0;

         if pl
            surf(X,Y,u);
            drawnow
            %pause
         end
      end
   end
end
