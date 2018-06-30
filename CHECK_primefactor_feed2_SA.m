% ------------------------------------------------------------
% symbolic compution of factorization problem
% Created by: Suksmono
% problem adopted from:
% https://qiita.com/YuichiroMinato/items/9c9fba2b5bc978ec9e73
% the idea: 
% (1) a user enter his/her problem
% (2) server do symbolic calculation and simulation: 
%     (a)P0 inc k-body -> P1 with only 2-bdy in s-domain
%     (b) extract coefficients {b, hi, Jij}
%           -> fed to simulator
%           -> deliver to user for D-Wave input
% (3) server send back simulation results and parameter to user
% (4) user run his problem (b, hi, Jij) in D-wave, compare the
%     result from simulator
% ------------------------------------------------------------
clear all;
addpath functions;

% define symbolic binary variables 
syms q0 q1 q2;
% we will calculate the Hamiltonian
% H=(15 - p*q)^2, where p=(1+2q0+41q), q=(1+2q2)
% ------------------------------------------------------------
%  USER WRITES THE PROBLEM
% ------------------------------------------------------------
p=1+2*q0+4*q1;
q=1+2*q2;
H=(15-p*q)^2;

% ------------------------------------------------------------
%  SYMBOLIC COMPUTATION
% ------------------------------------------------------------

EQ0=expand(H);
%display result for checking
disp(sprintf('EQ0= %s',char(EQ0)));

%substitute qi^2->1
EQ1=simplify( ...
        subs(EQ0, {q0^2, q1^2, q2^2},{q0, q1, q2}) ...
        );

disp(sprintf('EQ1= %s',char(EQ1)));

% ------------------------------------------------------------
% k-body to 2-body transform
% ------------------------------------------------------------
syms q3 d;
H2b=simplify( ...
        subs(EQ1, {q1*q2},{q3}) ...
        + H2body(q1,q2,q3,d) ...
        );
disp(sprintf('H2b= %s',char(H2b)));

% ------------------------------------------------------------
% feed the result to be simulated with Ising System
% assume SA-simulator working on s=[-1,1] domain
% transform q={0,1} to s{-1,1}
% ------------------------------------------------------------
syms s0 s1 s2 s3;

% ------------------------------------------------------------
% considering d>max(H), letting q0=q1=q2=0, 
% H=(15-1)^2 => d=15^2
% ------------------------------------------------------------

H2b1=simplify(subs(H2b,{d},{255}));
disp(sprintf('H2b1= %s',char(H2b1)));

H2s=simplify( ...
        subs(H2b1, {q0, q1, q2, q3}, ...
                  {q2s(s0), q2s(s1), q2s(s2), q2s(s3)}) ...
        );
disp(sprintf('H2s= %s',char(H2s)));
[cx,tx]=coeffs(H2s,[s0 s1 s2 s3]);
% ------------------------------------------------------------
% terms and coeffs
% cx =[ 4, -14, 32, 4, 255/4, -255/2, 431/4, -255/2, 415/4, -271/2, 1141/4]
% tx =[ s0*s1, s0*s2, s0*s3, s0, s1*s2, s1*s3, s1, s2*s3, s2, s3, 1]
% copy to hi, Jij, b

% ------------------------------------------------------------
% EXTRACT THE RESULTS FROM SYMBOLIC SOLUTION
% ------------------------------------------------------------
% number of  qubits: q0, q1, q2, q3, => 4
NQ=4;
Jij=zeros(NQ,NQ); % init Jij: coupling coefficients
hi=zeros(NQ,1);
b=0; % bias

NT=length(cx);
for m=1:NT
    tTerm=tx(m);
    cTerm=cx(m);
    %% strange case
    if double(tTerm == s0)==1
        hi(1)=cTerm;
    end;
    if double(tTerm == s1)==1
        hi(2)=cTerm;
    end;
    if double(tTerm == s2)==1
        hi(3)=cTerm;
    end;
    if double(tTerm == s3)==1
        hi(4)=cTerm;
    end;
    %%
    if double(tTerm == s0*s1)==1
        Jij(1,2)=cTerm;
    end;
    if double(tTerm == s0*s2)==1
        Jij(1,3)=cTerm;
    end;
    if double(tTerm == s0*s3)==1
        Jij(1,4)=cTerm;
    end;
    if double(tTerm == s1*s2)==1
        Jij(2,3)=cTerm;
    end;
    if double(tTerm == s1*s3)==1
        Jij(2,4)=cTerm;
    end;
    if double(tTerm == s2*s3)==1
        Jij(3,4)=cTerm;
    end;
    %
    if double(tTerm == 1)==1
        b=cTerm;
    end;
end;

% ------------------------------------------------------------
% SEND THE RESULTS (Jij, hi, b) TO SIMULATOR (S.A.)
% ------------------------------------------------------------
ITR1=100;%2*10*1000; %set 1 for random init
vStart=0;
%% create a vector of [(N-1)x(N-1),1] spin with random {-,+] values
vSpin0=vRandSpin(NQ);

% Etresh= 0;%-2*(4*K-1)*(2*K-1);%N/4; % target energy
% annealing schedule
maxIter=1*10*100;
Ptresh=pMCSchedule(maxIter,0.5);
%-- Ptresh is the schedule
k=0;
vSpin=vSpin0;
% % init energy setting
Ecurr = vEnergyBx(b, hi, Jij,vSpin);
% NFLIP=N-NO;
Etresh=0;
while k<maxIter && Ecurr>Etresh
    k=k+1;
    Eprev=Ecurr;
    % try Monte Carlo
    % --flip a row randomly
    RD=randperm(NQ);
    % copy current spin to a template
    tvSpin=vSpin;
    tvSpin(RD(1)) = -1*tvSpin(RD(1));
    Ecurr=vEnergyBx(b, hi, Jij,tvSpin);
    if (Ecurr<Eprev)||rand>Ptresh(k) %accept if lower
        vSpin=tvSpin;
    end;
    Esys=vEnergyBx(b, hi, Jij,vSpin);
    Ecurr=Esys;
    vE(k)=Esys;
    disp(sprintf('LQ%d: Iter->%d, Tresh=%1.4f, E=%4.2f',...
        NQ,k,Ptresh(k),double(Esys)));
%     disp(sprintf('H%d: Iter->%d',N,k));
end;
% ------------------------------------------------------------
% DISPLAY THE RESULTS
% ------------------------------------------------------------
qSpin=vs2q(vSpin);
disp(sprintf('\nq0= %g,q1= %g,q2= %g,q3= %g',qSpin));
% the factors:  p=1+2*q0+4*q1; q=1+2*q2;

vp=subs(p,{q0,q1},{qSpin(1),qSpin(2)});
vq=subs(q,{q2},{qSpin(3)});
disp(sprintf('Factors of 15 are: p=%g and q= %g',double(vp),double(vq)));
