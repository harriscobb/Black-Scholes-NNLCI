function U = MCS2D(sigx,sigy,rho,r,Nx)

    %========== 1. Setup tilde{x}, tilde{y} grids =======================
    T=1;
    K=0.05;
    Ny=Nx;
    dx = 1/Nx;      % spacing in tilde{x}
    dy = 1/Ny;      % spacing in tilde{y}
    tildeX = linspace(0,1,Nx+1);   % 0..1
    tildeY = linspace(0,1,Ny+1);   % 0..1

    %========== 2. Time-stepping and data structures ====================
    Nt = Nx*2;
    dt = T/Nt;
    U = zeros(Nx+1,Ny+1);

    %========== 3. Terminal/Initial condition: worst-of put payoff ======
    for i=1:Nx + 1
        for j=1:Ny + 1
            s1 = - log(tildeX(i));
            s2 = - log(tildeY(j));
            U(i,j) = max(K - min(s1, s2), 0);
        end
    end


    %========== 4. PDE splitting:  A = A0 + A1 + A2 =====================
    % A1: x-direction operator  (includes -r/2 * I)
    % A2: y-direction operator  (includes -r/2 * I)
    % A0: cross-derivative operator + any leftover pieces


    for timesteps = 1:Nt
        %Step 1: Y0 = U_{n-1}+dt*A*U_{n-1}
        Y0=U;
        for i=2:Nx
            x=tildeX(i);
            for j=2:Ny
                y=tildeY(j);
                Y0(i,j)=U(i,j)+dt*0.5*(sigx^2*log(x)^2*(x*(U(i+1,j)-U(i-1,j))/(2*dx)+x^2*(U(i+1,j)-2*U(i,j)+U(i-1,j))/(dx^2))+...
                                   sigy^2*log(y)^2*(y*(U(i,j+1)-U(i,j-1))/(2*dy)+y^2*(U(i,j+1)-2*U(i,j)+U(i,j-1))/(dy^2))+...
                                   2*rho*sigx*sigy*log(x)*log(y)*x*y*(U(i+1,j+1)+U(i-1,j-1)-U(i+1,j-1)-U(i-1,j+1))/(4*dx*dy))+...
                                   dt*r*(log(x)*x*(U(i+1,j)-U(i-1,j))/(2*dx)+log(y)*y*(U(i,j+1)-U(i,j-1))/(2*dy))-dt*r*U(i,j);
        
            end
        end
    
    
    
    
        %Step 2: (I-theta*dt*A1)*Y1=Y0-theta*dt*A1*U_{n-1}
        %A1 has 1/2*sigx^2*log(x)^2*(x*du/dx+x^2*d^2u/dx^2)+r*ln(x)*x*du/dx-r*u/2
        Y1=U;
        a1 = zeros(Nx+1,1);
        b1 = zeros(Nx+1,1);
        c1 = zeros(Nx+1,1);
        f1 = zeros(Nx+1,1);
        for i=1:Nx+1
            x=tildeX(i);
            a1(i) = 0.5*sigx^2*log(x)^2*(-x/(2*dx)+x^2/(dx^2))-r*log(x)*x/(2*dx);
            b1(i) = 0.5*sigx^2*log(x)^2*(-2*x^2/(dx^2))-r/2;
            c1(i) = 0.5*sigx^2*log(x)^2*(x/(2*dx)+x^2/(dx^2))+r*log(x)*x/(2*dx);
        end
        %b1(Nx)=b1(Nx)+c1(Nx);
        for j=2:Ny
            for i=2:Nx
                %i0=min(max(2,i),Nx);
                f1(i) = Y0(i,j)-0.5*dt*(a1(i)*U(i-1,j)+b1(i)*U(i,j)+c1(i)*U(i+1,j));
            end
            f1(Nx)=f1(Nx)+0.5*dt*c1(Nx)*U(Nx+1,j);
            Y1(2:Nx,j) = thomas_algorithm(-0.5*dt*a1(3:Nx),ones(Nx-1,1)-0.5*dt*b1(2:Nx),-0.5*dt*c1(2:Nx-1),f1(2:Nx));
        end
    
                
        %Step 3: (I-dt/2*A1)*Y1=Y0-dt/2*A1*U_{n-1}
        %A2 has 1/2*sigy^2(log(y)^2*y*du/dy+y^2*d^2u/dy^2)+r*ln(y)*y*du/dy-r*u/2
        Y2=U;
        a2 = zeros(Ny+1,1);
        b2 = zeros(Ny+1,1);
        c2 = zeros(Ny+1,1);
        f2 = zeros(Ny+1,1);
        for j=2:Ny
            a2(j)=0.5*sigy^2*log(tildeY(j))^2*(-tildeY(j)/(2*dy)+tildeY(j)^2/(dy^2))-r*log(tildeY(j))*tildeY(j)/(2*dy);
            b2(j)=0.5*sigy^2*log(tildeY(j))^2*(-2*tildeY(j)^2/(dy^2))-r/2;
            c2(j)=0.5*sigy^2*log(tildeY(j))^2*(tildeY(j)/(2*dy)+tildeY(j)^2/(dy^2))+r*log(tildeY(j))*tildeY(j)/(2*dy);
        end
        %b2(Ny)=b2(Ny)+c2(Ny);
        for i=2:Nx
            for j=2:Ny
                %j0=min(max(2,j),Ny);
                f2(j) = Y1(i,j)-dt/2*(a2(j)*U(i,j-1)+b2(j)*U(i,j)+c2(j)*U(i,j+1));
            end
            f2(Ny)=f2(Ny)+0.5*dt*c2(Ny)*U(i,Ny+1);
            Y2(i,2:Nx) = thomas_algorithm(-0.5*dt*a2(3:Ny),ones(Ny-1,1)-0.5*dt*b2(2:Ny),-0.5*dt*c2(2:Ny-1),f2(2:Ny));
        end
    
        %Step 4: Y0hat = Y0 +1/2 dt A0(Y2-U_{n-1})
        Y0hat = U;
        for i = 2:Nx
            %i0=min(max(2,i),Nx);
            for j = 2:Ny
                %j0=min(max(2,j),Ny);
                Y0hat(i,j) = Y0(i,j)+dt/2*rho*sigx*sigy*log(tildeX(i))*log(tildeY(j))*tildeX(i)*tildeY(j)*(Y2(i+1,j+1)-U(i+1,j+1)+Y2(i-1,j-1)-U(i-1,j-1)-Y2(i+1,j-1)+U(i+1,j-1)-Y2(i-1,j+1)+U(i-1,j+1))/(4*dx*dy);
            end
        end
    
        %Step 6: (skip 5 since theta = 1/2)
        % (I-dt/2*A1)*Y1hat = Y0hat-dt/2*A1*U_{n-1}
        Y1hat=U;
        for j=2:Ny
            for i=2:Nx
                %i0=min(max(2,i),Nx);
                f1(i) = Y0hat(i,j)-dt/2*(a1(i)*U(i-1,j)+b1(i)*U(i,j)+c1(i)*U(i+1,j));
            end
            f1(Nx)=f1(Nx)+0.5*dt*c1(Nx)*U(Nx+1,j);
            Y1hat(2:Nx,j) = thomas_algorithm(-0.5*dt*a1(3:Nx),ones(Nx-1,1)-0.5*dt*b1(2:Nx),-0.5*dt*c1(2:Nx-1),f1(2:Nx));
        end
    
        %Step 7: (I-dt/2*A2)*Y2hat = Y1hat-dt/2*A2*U_{n-1}
        Y2hat=U;
        for i=2:Nx
            for j=2:Ny
                %j0=min(max(2,j),Ny);
                f2(j) = Y1hat(i,j)-dt/2*(a2(j)*U(i,j-1)+b2(j)*U(i,j)+c2(j)*U(i,j+1));
            end
            f2(Ny)=f2(Ny)+0.5*dt*c2(Ny)*U(i,Ny+1);
            Y2hat(i,2:Ny) = thomas_algorithm(-0.5*dt*a2(3:Ny),ones(Ny-1,1)-0.5*dt*b2(2:Ny),-0.5*dt*c2(2:Ny-1),f2(2:Ny));
        end
    
        %Step 8
        U(2:Nx,2:Ny)=Y2hat(2:Nx,2:Ny);
        for i=1:Nx+1
            U(i,Ny+1)=K*exp(-r*timesteps*dt);
        end
        for j=1:Ny+1
            U(Nx+1,j)=K*exp(-r*timesteps*dt);
        end
    end

    % Plot
    figure;
    x=linspace(0,3*K,Nx+1);
    y=linspace(0,3*K,Ny+1);

    % Create the new transformed grid points
    newX = exp(-x);  % Transform x to exp(-x)
    newY = exp(-y);  % Transform y to exp(-y)
    
    % Initialize the new matrix u
    %u = zeros(Nx + 1, Ny + 1);
    
    % Perform interpolation using interp2
    [X, Y] = meshgrid(tildeX, tildeY);       % Original grid
    [X_new, Y_new] = meshgrid(newX, newY);   % New grid for u
    u = interp2(X, Y, U, X_new, Y_new, 'linear');  % Interpolation
    [X0,Y0]=meshgrid(x,y);
    %u_exact=zeros(Nx+1,Ny+1);
    tau=T;

    % Initialize u_exact
    u_exact = zeros(Nx+1, Ny+1);
    
    for i = 1:Nx+1
        for j = 1:Ny+1
            
            % Extract asset prices from the mesh
            s1 = x(i);  % s1
            s2 = y(j);  % s2
            u_exact(i,j) = Worst_of_Exact(sigx,sigy,rho,r,s1,s2,K,tau);
            
        end
    end


    %[X, Y] = meshgrid(tildeX, tildeY);  
    surf(X0, Y0, u_exact-u);
    xlabel('x'); ylabel('y'); zlabel('u_exact-u');
    title('Difference Between Numerical and Exact Solutions');
    %Infinity norm relative to K
    max(abs(u-u_exact),[],"all")/K

end