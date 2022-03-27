%LINEARIZE   Linearize nonlinear system.
%   [A,B] = LINEARIZE(model,x0,u0,tangent) returns A,B matrices
%   of the system linearization defined in 'model' around x0,u0. If tangent
%   is 1, then the linearization is done in the the parameter space (tangent).
%   
%   OUTPUTS:
%       - A and B matrices corresponding to the linear system x_dot = Ax+Bu
%       if tangent=0 or y_dot = Ay+Bu if tangent=1. A =partial f/ partial x
%       B = partial f/ partial u, where f(x,u)=x_dot.
%   INPUTS:
%       - model: model function from Featherstone Toolbox.
%       - x0: state where the linearization should be done.
%       - u0: input where the linearization should be done.
%       - tangent: linearization is done in tangent if 1, linearization is
%       done in ambient otherwise.
%
%   Ricard Bordalba, IRI-UPC, rbordalba@iri.upc.edu
function [A,B] = Linearize(model,y0,x0,u0,tangent)

epsilon = 1e-7;

switch tangent
    case 1
        integration  = 1; % normal
        f0 = ydot([],y0,x0,u0,model,integration);
        A = [];
        B = [];
        for i=1:length(y0)
            y = y0;
            temp = y(i) + epsilon;
            epsilon1 = temp - y(i);
            y(i)=y(i)+epsilon1;
            
            f1 = ydot([],y,x0,u0,model,integration);
            
            A = [A (f1-f0)/epsilon1];
        end
                
        for i=1:length(u0)
            u = u0;
            temp = u(i) + epsilon;
            epsilon1 = temp - u(i);
            u(i)=u(i)+epsilon1;
           
            f1 = ydot([],y0,x0,u,model,integration);
            
            B = [B (f1-f0)/epsilon1];   
        end       
        %A = gradest(@(y)ydot([],y,x0,u0,model,integration),y0);
        %B = gradest(@(u)ydot([],y0,x0,u,model,integration),y0);

        %[A,err] = jacobianest(@(y)ydot([],y,x0,u0,model,integration),y0);
        %[B,err] = jacobianest(@(u)ydot([],y0,x0,u,model,integration),u0);
    otherwise        
        integration  = 3; % normal ODE
        t = 1e-100;
        
        f0 = xdot([],x0,u0,model,integration);
        A = [];
        B = [];
        for i=1:length(x0)
            x = x0;
            x(i)=x(i)+j*t;
            A(:,i) = imag(xdot([],x,u0,model,integration))/t;
                        
        end
        
        for i=1:length(u0)
            u = u0;
            u(i)=u(i)+j*t;
            B(:,i) = imag(xdot([],x0,u,model,integration))/t;
        end  
        
%         [A,err] = jacobianest(@(x)xdot([],x,u0,model,integration),x0);
%         [B,err] = jacobianest(@(u)xdot([],x0,u,model,integration),u0);
end
end




% function [A,B] = Linearize(model,y0,x0,u0,tangent)
% 
% epsilon = 1e-7;
% 
% switch tangent
%     case 1
%         integration  = 1; % normal
%         f0 = ydot([],y0,x0,u0,model,integration);
%         A = [];
%         B = [];
%         for i=1:length(y0)
%             y = y0;
%             temp = y(i) + epsilon;
%             epsilon1 = temp - y(i);
%             y(i)=y(i)+epsilon1;
%             
%             f1 = ydot([],y,x0,u0,model,integration);
%             
%             A = [A (f1-f0)/epsilon1];
%         end
%         
%         for i=1:length(u0)
%             u = u0;
%             temp = u(i) + epsilon;
%             epsilon1 = temp - u(i);
%             u(i)=u(i)+epsilon1;
%            
%             f1 = ydot([],y0,x0,u,model,integration);
%             
%             B = [B (f1-f0)/epsilon1];   
%         end       
%         %A = gradest(@(y)ydot([],y,x0,u0,model,integration),y0);
%         %B = gradest(@(u)ydot([],y0,x0,u,model,integration),y0);
% 
%         %[A,err] = jacobianest(@(y)ydot([],y,x0,u0,model,integration),y0);
%         %[B,err] = jacobianest(@(u)ydot([],y0,x0,u,model,integration),u0);
%     otherwise        
%         integration  = 3; % normal ODE
%         
%         f0 = xdot([],x0,u0,model,integration);
%         A = [];
%         B = [];
%         for i=1:length(x0)
%             x = x0;
%             temp = x(i) + epsilon;
%             epsilon1 = temp - x(i);
%             x(i)=x(i)+epsilon1;
%             
%             f1 = xdot([],x,u0,model,integration);
%             
%             A = [A (f1-f0)/epsilon1];
%         end
%                 
%         for i=1:length(u0)
%             u = u0;
%             temp = u(i) + epsilon;
%             epsilon1 = temp - u(i);
%             u(i)=u(i)+epsilon1;
%             
%             f1 = xdot([],x0,u,model,integration);
%             
%             B = [B (f1-f0)/epsilon1];   
%         end  
%         
% %         [A,err] = jacobianest(@(x)xdot([],x,u0,model,integration),x0);
% %         [B,err] = jacobianest(@(u)xdot([],x0,u,model,integration),u0);
% end
% end
