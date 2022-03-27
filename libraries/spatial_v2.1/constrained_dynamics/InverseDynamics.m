%INVERSEDYNAMICS   Inverse dynamics in a constrained system.
%   u = InverseDynamics(x,qdotdot,model,method) returns the input force u
%   needed at state x to achieve an acceleration qdotdot in a
%   constrianed systems.
%
%   OUTPUTS:
%       - u is the input force in the inverse dynamics problem.  
%       In constrained systems, size(u) is not equal to size(q) and use the
%       projection matrix B to obtain tau=Bu.
%   INPUTS:
%       - x: state.
%       - qdotdot: desired acceleration. It should be a feasible
%       acceleration for a constrained system (Phidotdot = 0).
%       - model: model function from Featherstone Toolbox.
%       - method: determines the Inverse Dynamics method to be used.
%           1. Projected Inverse Dynamics [Aghili ]
%           2. IRI approach
%           3. Featherstone approach with local coordinates map
%           4. Featherstone approach with coordinate partitioning
%           5. Minimum norm projected inverse dynamics [icra2009-Aghili]
%
%   Ricard Bordalba, IRI-UPC, rbordalba@iri.upc.edu


function u = InverseDynamics(x,qdotdot,model,method)
  
n_q = length(x)/2;
q = x(1:n_q); %position
v = x(n_q+1:2*n_q); %velocity


if isfield(model,'loop')  
    [M, f, Phi, Phi_q, xi] = HandC(model,q,v);
else
    % Openloop dynamic terms
    [M, f] = HandC(model,q,v,[]) ;

    % constraints equations to close the loop
    [~, Phi_q, xi] = model.constraint(model,q,v);
end
    
n_lambda = length(xi);
n_q = length(q);
B = model.B;   
BB = B*B';
n_u = size(B,2);


tau_ID = M*qdotdot + f;
switch method  
    case 1 % Projected inverse
        P = null(Phi_q)*null(Phi_q)';%eye(length(Phi_q))-pinv(Phi_q)*Phi_q;
        tau = P*tau_ID;  % potent input force in full dimension
        %tau = PB u  and we want to find u
        P1=(P*B);
        %u = ((P1'*P1)^-1*P1')*tau;
        u=(P1)\tau;
        %u=pinv(P1)*tau;
    case 2 % Our approach
        %temp = [B Phi_q']\tau_ID;
        temp = pinv([B Phi_q'])*tau_ID;
        u = temp(1:n_u);
   case 3 % Featherstone & local coordinates   
    %Dpsi'*(tau_ID-tau)=0 and tau=B*u ---> (Dpsi'*B)*u = DPsi'*tau_ID
        U0 = null(Phi_q);
        DPsi = [Phi_q;U0']\[zeros(size(U0,1)-size(U0,2),size(U0,2));eye(size(U0,2))];%U0 in case the tangent space center is the current state 
        %u = (DPsi'*B)\(DPsi'*tau_ID);
        u = pinv(DPsi'*B)*(DPsi'*tau_ID);
   case 4 % Featherstone & coordinate partitioning
       % Choosing actuated joints as independent
        dep = find(~any(B, 2));    
        Phi_qi = Phi_q*B; 
        Phi_qd = Phi_q(:,dep);
        A = B*eye(n_u);
        %A(dep,:)=-Phi_qd\Phi_qi;         
        A(dep,:)=-pinv(Phi_qd)*Phi_qi;       
        %u = (A'*B)\(A'*tau_ID);
        u = pinv(A'*B)*(A'*tau_ID);
   case 5 % Projected minimum torque
       % Choosing actuated joints as independent
        P = null(Phi_q)*null(Phi_q)';%eye(length(Phi_q))-pinv(Phi_q)*Phi_q;
        tau = P*tau_ID;  % potent input force in full dimension
        umax = 1000; % big number 
        options = optimoptions('fmincon','Display','off');
        u = fmincon(@(u)norm(u),pinv(P*B)*tau,[eye(n_u);-eye(n_u)],umax*ones(2*n_u,1),P*B,tau,[],[],[],options);       
    case 6
        %lambda = -(B_bar'*Phi_q')\(B_bar'*tau_ID);
        BB_bar = eye(n_q)-BB;
        B_bar = BB_bar(:,any(BB_bar));
        
        Phi_qa = Phi_q*B; 
        Phi_qr = Phi_q*B_bar;
        
        tmp = pinv([eye(n_u)           Phi_qa';
                    zeros(n_q-n_u,n_u) Phi_qr'])*[B';B_bar']*tau_ID;
         
        u = tmp(1:n_u);

        %lambda = -pinv(B_bar'*Phi_q')*(B_bar'*tau_ID);
        %u = (B'*tau_ID)+(Phi_q*B)'*lambda;
        
end    
end

