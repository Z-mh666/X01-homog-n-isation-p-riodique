% =====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions de
% Neumann sur le maillage nom_maillage.msh
%
% | -Delta  u + u= f,   dans \Omega
% |         du/dn = 0,   sur le bord
%
% =====================================================


% lecture du maillage et affichage
% ---------------------------------
nom_maillage = 'geomCarre.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);


% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
LL = zeros(Nbpt,1);     % vecteur second membre

% boucle sur les triangles
% ------------------------
for i=1:Nbtri
  % Coordonnees des sommets du triangles
  % A COMPLETER
  k = Numtri(i,1);
  l = Numtri(i,2);
  m = Numtri(i,3);
  
  S1 = Coorneu(k,:);
  S2 = Coorneu(l,:);
  S3 = Coorneu(m,:);
  % calcul des matrices elementaires du triangle l 
  
  %Kel=matK_elem(S1, S2, S3);
  Kel = matK2_elem(S1,S2,S3);
           
  Mel=matM_elem(S1, S2, S3);
    
  % On fait l'assemmblage de la matrice globale et du second membre
  % A COMPLETER
  MM(k,k) += Mel(1,1);
  MM(l,l) += Mel(2,2);
  MM(m,m) += Mel(3,3);
  MM(k,l) += Mel(1,2);
  MM(k,m) += Mel(1,3);
  MM(l,m) += Mel(2,3);
  MM(l,k) += Mel(2,1);
  MM(m,k) += Mel(3,1);
  MM(m,l) += Mel(3,2);
  
  KK(k,k) += Kel(1,1);
  KK(l,l) += Kel(2,2);
  KK(m,m) += Kel(3,3);
  KK(k,l) += Kel(1,2);
  KK(k,m) += Kel(1,3);
  KK(l,m) += Kel(2,3);
  KK(l,k) += Kel(2,1);
  KK(m,k) += Kel(3,1);
  KK(m,l) += Kel(3,2);

end % for l

% Calcul du second membre L
% -------------------------
	% A COMPLETER
	% utiliser la routine f.m
FF = f(Coorneu(:,1),Coorneu(:,2));
LL = MM*FF;

% inversion
% ----------
UU = (MM+KK)\LL;

% visualisation
% -------------
affiche(UU, Numtri, Coorneu, sprintf('Neumann - %s', nom_maillage));

validation = 'non';
% validation
% ----------
if strcmp(validation,'non')
UU_exact = cos(pi*Coorneu(:,1)).*cos(2*pi*Coorneu(:,2));
% Calcul de l erreur L2
% A COMPLETER
errL2 = log(sqrt(UU'*MM*UU+UU_exact'*MM*UU_exact-2*UU'*MM*UU_exact)/sqrt(UU_exact'*MM*UU_exact))
%err_L2 = [-1.5731,-3.2949,-4.8630,-6.7224];%stocker les erreurs L2 pour A = 1
err_L2 = [-0.50531,-1.3261,-1.4388,-1.5049];%stocker les erreurs L2 pour A = sin(2*pi*x)*sin(2*pi*y)+2

figure(2)
plot([log(1/0.2),log(1/0.1),log(1/0.05),log(1/0.02)],err_L2);
xlabel('log(1/h)');
ylabel('erreur L2');

% Calcul de l erreur H1
% A COMPLETER
errH1 = log(sqrt(UU'*KK*UU+UU_exact'*KK*UU_exact-2*UU'*KK*UU_exact)/sqrt(UU_exact'*KK*UU_exact))
%err_H1 = [-2.1761,-3.4180,-4.7065,-6.3735];  %stocker les erreurs H1 pour A = 1
err_H1 = [-0.76800,-0.94492,-1.0757,-1.1548];%stocker les erreurs H1 pour A = sin(2*pi*x)*sin(2*pi*y)+2

figure(3)
plot([log(1/0.2),log(1/0.1),log(1/0.05),log(1/0.02)],err_L2);
xlabel('log(1/h)');
ylabel('erreur H1');
endif



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

