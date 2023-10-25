function Aeff = tenseur_homo(K,wk,wj,yk,yj)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aeff :
% Calcul de la matrice Aeff .
%          
% INPUT * A : le tenseur caracteristique
%         K : matrice de rigidite
%     wk,wj : solution du probleme de cellule
%     yk,yj : coordonnee des sommets
%         
%
% OUTPUT - Aeff: valeur de la matrice homogeneise.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Aeff = zeros(2,2);
Aeff(1,1) = (K*yk)'*yk+2*(K*yk)'*wk+(K*wk)'*wk;
Aeff(1,2) = (K*yk)'*yj+(K*yk)'*wj+(K*wk)'*yj+(K*wk)'*wj;
Aeff(2,1) = Aeff(1,2);
Aeff(2,2) = (K*yj)'*yj+2*(K*yj)'*wj+(K*wj)'*wj;


endfunction
