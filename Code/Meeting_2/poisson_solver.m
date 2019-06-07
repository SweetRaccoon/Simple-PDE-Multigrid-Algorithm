function u = poisson_solver(u,f,epsilon,pl)
%Poisson solver -Delta u = f, u=0 on boundary of box
%f should be square nxn matrix

   n = max(size(f)); %number of grid points
   dx = 1/(n-1); %Grid resolution
   [X,Y] = meshgrid(0:dx:1);

   %Time step
   dt = 0.8*dx^2/4;
   
   iteration_count = 0;
   err = epsilon + 1;
   while err > epsilon

      iteration_count = iteration_count + 1;

      %Compute partial derivatives
      uxx = (u([2:n,n],:) - 2*u + u([1,1:n-1],:))/dx^2;
      uyy = (u(:,[2:n,n]) - 2*u + u(:,[1,1:n-1]))/dx^2;
      Laplacian = uxx + uyy;
      err = abs(Laplacian + f);
      err = max(max(err(2:n-1,2:n-1)));
      
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
