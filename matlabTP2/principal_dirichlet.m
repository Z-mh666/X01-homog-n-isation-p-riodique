% =====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions de
% Dirichlet sur le maillage nom_maillage.msh
%
% | -\Delta u + u= f,   dans \Omega
% |         u = 0,   sur le bord
%
% =====================================================


% lecture du maillage et affichage
% ---------------------------------

%nom_maillage = 'geomCarre_per.msh';  %(0,2)^2
nom_maillage = 'geomCarre_per2.msh'; %(0,1)^2
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
LL = zeros(Nbpt,1);     % vecteur second membre
e = 0.001;

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
  Kel = matK2_elem(S1,S2,S3,e);
           
  Mel = matM_elem(S1,S2,S3);
    
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
  %
FF = f(Coorneu(:,1),Coorneu(:,2));
LL = MM*FF;

% Projection sur l espace V_0
% ——————————————————�?
% matrice de projection 
num = 1;
PP = zeros(Nbpt-Nbaretes,Nbpt);
for i=1:Nbpt
  if Refneu(i)==0
    PP(num,i) = 1;
    num += 1;
  endif
endfor

AA = KK;
AA0 = PP*AA*PP';
LL0 = PP*LL;

% inversion
% ----------
UU0 = AA0\LL0;

% Expression de la solution dans toute la base
% ——————�?
UU = PP'*UU0;
UU_exact = sin(pi*Coorneu(:,1)).*sin(pi*Coorneu(:,2));

% visualisation
% -------------
affiche(UU, Numtri, Coorneu, sprintf('Dirichlet - %s', nom_maillage));

%affichage de la solution exacte
affiche(UU_exact, Numtri, Coorneu, sprintf('Dirichlet - %s', nom_maillage));

validation = 'non';
% validation
% ----------
if strcmp(validation,'oui')
UU_exact = sin(pi*Coorneu(:,1)).*sin(pi*Coorneu(:,2));
% Calcul de l erreur L2
% A COMPLETER
errL2 = sqrt(UU'*MM*UU+UU_exact'*MM*UU_exact-2*UU'*MM*UU_exact)


% Calcul de l erreur H1
% A COMPLETER
errH1 = sqrt(UU'*KK*UU+UU_exact'*KK*UU_exact-2*UU'*KK*UU_exact)

% attention de bien changer le terme source (dans FF)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

