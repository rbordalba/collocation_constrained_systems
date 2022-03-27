%orthonormalize   Orthonormalize C with respect to F.
%   [F1,C1]=orthonormalize(F,C) returns an orthonormal basis C1 wrt F
%
%   OUTPUTS:
%       - F1: New basis of F given by rows that is normalized 
%       (F*F'~I  while F1*F1'=I)
%       - C1: New basis of C given by columns that is orthogonal to F
%       (F1*C1=0) and that is normalized (C1'*C1=0)
%   INPUTS:
%       - F: basis given by rows
%       - C: basis given by columns
%
%   EXAMPLE: 
%       Given DF(k) and the previos tg space basis U(k-1), we want
%       to obtain a basis U(k) that is orthogonal wrt to DF (DF*U=0) and 
%       that is also normalized (F1*F1'=I).

%       [~,Unew]=orthonormalize(DF,Uprev);
%
%   Ricard Bordalba, IRI-UPC, rbordalba@iri.upc.edu

% F is given by rows
% C is given by columns
% At the end of this process
%   C1'*C1=I
%   F1*F1'=I
%   F1*C1=0

function [F1,C1]=orthonormalize(F,C)

  nf=size(F,1); % in rows
  nc=size(C,2); % in columns

  if size(F,2)~=size(C,1)
    error('wrong size');
  end
  
  F1=F;
  if ~isnumeric(C)
    import casadi.*
    C1 = SX.zeros(size(C,1),size(C,2)); % SX or MX depending on prev. defined
  else
    C1=C;
  end

  for i=1:nf
    for j=1:i-1  
      F1(i,:)=F1(i,:)-F1(j,:)*F1(i,:)'*F1(j,:);
    end
    F1(i,:)=F1(i,:)/norm(F1(i,:));
  end
  
  for i=1:nc
    j=1;
    C1(:,i)=C(:,i)-F1(j,:)'*(C(:,i)'*F1(j,:)')/(F1(j,:)*F1(j,:)');
    for j=2:nf
      C1(:,i)=C1(:,i)-F1(j,:)'*(C1(:,i)'*F1(j,:)')/(F1(j,:)*F1(j,:)');
    end
    %C1=C-F1'*F1*C; ------>> leads to C1'*C1~I
    for j=1:i-1  
      C1(:,i)=C1(:,i)-C1(:,j)*C1(:,i)'*C1(:,j);
    end
    C1(:,i)=C1(:,i)/norm(C1(:,i));
  end