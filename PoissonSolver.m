%% PoissonSolver
%  Solve (approximately) the Poisson equation for the pressure correction.
%  Find a field p such that div(grad(p)) = div(u)

function [p,nIt,eps]=PoissonSolver(div,rhoPhi,rho,dt,dx,omega,nItMax,epsMax) 
    % Initialisation
    [Ima,Jma]=size(rhoPhi.rhoU);
    p = zeros(Ima,Jma);

    % Loop over iteration steps. 
    eps = 1.0E+9; nIt = 0;
    while (nIt < nItMax) && (eps > epsMax)
          % Increment iteration counter
            nIt = nIt + 1;
            
            % Store pressure field of last iteration step
            pold = p;

            % Compute new pressure with Jacobi point iteration
            p(2:Ima-1,2:Jma-1) = 0.25 * (pold(3:Ima,2:Jma-1) + pold(1:Ima-2,2:Jma-1) ...
                                       + pold(2:Ima-1,3:Jma) + pold(2:Ima-1,1:Jma-2) ...
                                      - omega.*(div(2:Ima-1,2:Jma-1)*rho*dx./dt) );

            % Set zero gradient BC at the inlet and fixed value of zero at all            
            p(1,:) = p(2,:); p(Ima,:) = -p(Ima-1,:);
            p(:,1) = p(:,2); p(:,Jma) = -p(:,Jma-1);
                  
            % Compute maximum error
            eps = max (max (abs (p-pold)));  
    end
end   





























