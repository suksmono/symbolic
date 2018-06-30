function tSched=pMCSchedule(N,pStart)
% create probabilistic schedule
% 1/4 part: pStart-0.99
% 1/4 part: 0.999
% 1/4 part: 0.9999
% 1/4 part: 1.0
N025=floor(N/4);
tt1=1:N025;
tt2=N025+1:2*N025;
tt3=2*N025+1:3*N025;
tt4=3*N025+1:N;
tSched=ones(1,N);
tSched(tt4)=1.;
tSched(tt3)=0.9999;
tSched(tt2)=0.999;
%
dGap=1.0-pStart;
xt1=tt1/N025;
st1=exp(-xt1);
st1=pStart+dGap*(1-st1)/max(1-st1);
tSched(1:N025)=st1;
%st1=st1/max(st1);
%st1=1-dGap*st1/(max(st1)-min(st1));