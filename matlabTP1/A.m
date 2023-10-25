function val = A(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_A :
% Evaluation de la matrice A .
%
% SYNOPSIS val = mat_A(x,y)
%          
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%         
%
% OUTPUT - val: valeur de la matrice sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A COMPLETER


%val = [1 0;0 1];

%val = [2+sin(2*pi*x/10) 0;0 2+sin(2*pi*x/10)];

%val = [2+sin(2*pi*x/100) 0;0 2+sin(2*pi*x/100)];


val = sin(2*pi*x)*sin(2*pi*y)+2;
    
%val = sin(4*pi*x)*sin(4*pi*y)+2
    
%val = sin(8*pi*x)*sin(8*pi*y)+2
    
%val = sin(16*pi*x)*sin(16*pi*y)+2
    
%val = sin(32*pi*x)*sin(32*pi*y)+2

%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
