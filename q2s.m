function val_s= q2s(x)
% transform q={1,0} -> s={-1,+1}
% usage: (1) define symbolic variable x
%             syms x;
%        (2) fq=q2s(x); 
val_s=(1/2)*(1-x);
