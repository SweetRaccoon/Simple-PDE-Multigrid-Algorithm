function [u,res] = poisson_iterations(u,f,T,pl)
%Poisson solver -Delta u = f, u=0 on boundary of box
%f should be square nxn matrix
%Runs for T iterations starting from u

   n = max(size(f)); %number of grid points
   dx = 1/(n-1); %Grid resolution
   [X,Y] = meshgrid(0:dx:1);

   %Time step
   dt = 0.8*dx^2/4;
   
   for i = 1:T

      %Compute partial derivatives
      uxx = (u([2:n,n],:) - 2*u + u([1,1:n-1],:))/dx^2;
      uyy = (u(:,[2:n,n]) - 2*u + u(:,[1,1:n-1]))/dx^2;
      Laplacian = uxx + uyy;
      res = Laplacian + f;
      err = max(max(abs(res(2:n-1,2:n-1))));
      
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
