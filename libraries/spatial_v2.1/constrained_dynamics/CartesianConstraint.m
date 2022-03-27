function  [r, R, J, Jdot] = CartesianConstraint( model, q, qd )

J = zeros(6,model.NB);
Jdot = zeros(6,model.NB);
% Calculates cartesian constraint (end effector)
for i = 1:max(find(model.ee==1)) %model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
  vJ = S{i}*qd(i);
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    avp{i} = Xup{i} * zeros(6,1);  
    %Compute up transform i^X_i-1
    Xa{i} = Xup{i}; 
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    avp{i} = Xup{i}*avp{model.parent(i)} + crm(v{i})*vJ;   
    %Compute up transform from base i^X_0
    Xa{i} = Xup{i} * Xa{model.parent(i)};
  end  
  %transform Si to 0 coordinates    0^X_i * S{i}
  %J(:,i) = Xa{i}\S{i};
  %Jdot(:,i) = Xa{i}\(crm(v{i})*S{i});
end

ee = max(find(model.ee==1)); % find end effector
%dr = Xa{ee}\v{ee};
[r,R] = X2vec(Xa{ee}^-1*(model.Xee)^-1);

% remove dependent end-effector coordinates
if isfield(model,'planar')
   %J = J(model.cartesian,:);
   %Jdot = Jdot(model.cartesian,:);
   r = r(model.cartesian-3);    
end
end



function [r, R] = X2vec(X)

   R = X(1:3,1:3); rl = X(4:6,1:3);
   Sr = rl*R^-1; %obtain screw(r) in global coordinates

   r = 0.5*[Sr(3,2)-Sr(2,3);
            Sr(1,3)-Sr(3,1);
            Sr(2,1)-Sr(1,2)];
end