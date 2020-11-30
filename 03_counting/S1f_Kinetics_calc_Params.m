% Calculation of the photophysical parameters
% --Andreas Hartmann

clear all;

% parameter
tauC=0.15; % (s)

%% Values
%LacY
% kr1=12.913;errkr1=1.1945;
% kr2=1.8489;errkr2=0.1227;
% p1=0.5295;errp1=0.0496;
% p2=0.2217;errp2=0.0331;
% n=0.3047;errn=0.0213;
% kdkb=28.9248;errkdkb=0.3596;

% %LexAC3
% kr1=13.1443;errkr1=1.2048;
% kr2=1.6985;errkr2=0.1128;
% p1=0.6166;errp1=0.0535;
% p2=0.1926;errp2=0.0325;
% n=0.2195;errn=0.0201;
% kdkb=27.6361;errkdkb=0.3284;

% %LexAUV
% kr1=12.0378;errkr1=1.4132;
% kr2=1.7332;errkr2=0.1502;
% p1=0.6244;errp1=0.0711;
% p2=0.1843;errp2=0.0437;
% n=0.2427;errn=0.0208;
% kdkb=27.6554;errkdkb=0.3072;
% 
% % %All
% kr1=12.9816;errkr1=0.74146;
% kr2=1.7751;errkr2=0.07384;
% p1=0.58299;errp1=0.03104;
% p2=0.20173;errp2=0.020142;
% n=0.25;errn=0.0210;
% kdkb=27.5323;errkdkb=0.3722;
% 
% %LexA C3 + UV
% kr1=12.8194;errkr1=0.9111;
% kr2=1.7197;errkr2=0.0892;
% p1=0.6151;errp1=0.0412;
% p2=0.1895;errp2=0.0262;
% n=0.22714;errn=0.0203;
% kdkb=27.4024;errkdkb=0.38895;
%
% % %Dendra
% kr1=9.2857;errkr1=0.7896;
% kr2=1.3904;errkr2=0.0821;
% p1=0.3449;errp1=0.0362;
% p2=0.2647;errp2=0.0273;
% n=0.1126;errn=0.0089;
% kdkb=42.9576;errkdkb=0.3863;

%% Calulate errors
% ptc
rand_p1=errp1.*randn(10000,1)+p1;
rand_p2=1-rand_p1;
rand_kr1=errkr1.*randn(10000,1)+kr1;
rand_kr2=errkr2.*randn(10000,1)+kr2;

ptc=p1/(p1+p2)*exp(-kr1*tauC)+p2/(p1+p2)*exp(-kr2*tauC);
rand_ptc=rand_p1./(rand_p1+rand_p2).*exp(-rand_kr1.*tauC)+rand_p2./(rand_p1+rand_p2).*exp(-rand_kr2.*tauC);
errptc=std(rand_ptc);

% corrected n
rand_n=errn.*randn(10000,1)+n;
rand_ptc=errptc.*randn(10000,1)+ptc;

corrn=n/(ptc+n*(1-ptc));
rand_corrn=rand_n./(rand_ptc+rand_n.*(1-rand_ptc));
errcorrn=std(rand_corrn);

% rates
rand_kdkb=errkdkb.*randn(10000,1)+kdkb;

kd=kdkb*corrn;
rand_kd=rand_kdkb.*rand_corrn;
errkd=std(rand_kd);

kb=kdkb*(1-corrn);
rand_kb=rand_kdkb.*(1-rand_corrn);
errkb=std(rand_kb);

disp(['p_tc=(' num2str(ptc) '+-' num2str(errptc) ')']);
disp(['corr. eta=(' num2str(corrn) '+-' num2str(errcorrn) ')']);
disp(['kb=(' num2str(kb) '+-' num2str(errkb) ')s^-1']);
disp(['kd=(' num2str(kd) '+-' num2str(errkd) ')s^-1']);



