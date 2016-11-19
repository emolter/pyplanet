% This is Karpowicz's model from:
%   http://users.ece.gatech.edu/~psteffes/palpapers/karpowicz_data/water_model/karpowicz_h2o_model.m
function alpha_h2o=karpowicz_h2o_model(f,density_h2o,density_h2,density_he,T)
% This function also requires the vvwlineshape function originally written by Jim Hoffman (with removal of df factor). Also added shift term
% as us done in Rosenkranz's work.

mbars_to_bars=0.001;
inv_km_to_dB=4.342945;
convert_to_km=1e-4;
To=300;
Theta=300./T;
n=size(f,2);  %returns the number of columns in f
nones=ones(1,n);%Vector to get 2d arrays for multiplication

NA=6.0221415e23;%molecules/mol
M_amu=8.314472/0.46151805;%g/mol
isotope_partition=0.997317;
M_amu_he=4.0026;
M_amu_h2=2.01594;
        
P_h2o=(1.0/M_amu)*density_h2o*8.314472e-5*T;
P_He=(1.0/M_amu_he)*density_he*8.314310e-5*T;
P_H2=(1.0/M_amu_h2)*density_h2*8.314472e-5*T;
        
density_h2o=isotope_partition*(density_h2o/M_amu)*NA*(1.0/1e6); %need in molecules/cc     

% Center Frequencies in GHZ
f_o=[22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508, 443.0183, 448.0011, 470.8890, 474.6891, 488.4911, 556.9360, 620.7008, 752.0332, 916.1712];

% Line intensities
I_o=[.1314E-13, .2279E-11, .8058E-13, .2701E-11, .2444E-10,.2185E-11, .4637E-12, .2568E-10, .8392E-12, .3272E-11, .6676E-12, .1535E-08, .1711E-10, .1014E-08, .4238E-10];
% Temperature coefficients
E_o=[2.144, .668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405,3.597, 2.379, 2.852, .159, 2.391, .396, 1.441];


%self broadening parameters converted to bars
w_s=1/mbars_to_bars*[.01349, .01466, .01057, .01381, .01454, .009715, .00788, .01275, .00983, .01095, .01313, .01405, .011836, .01253, .01275];
x_s=[.61, .85, .54, .74, .89, .62, .50, .67, .65, .64, .72,1.0, .68, .84, .78];

w_H2=[2.395, 2.4000, 2.395, 2.395, 2.390, 2.395, 2.395, 2.395, 2.395, 2.395,2.395,2.395,2.395,2.395,2.395];
        
w_He=[0.67, 0.71, 0.67,0.67, 0.63, 0.67, 0.67, 0.67, 0.67,0.67,0.67, 0.67, 0.67,0.67, 0.67];
x_H2=[0.900, 0.950, 0.900, 0.900, 0.850, 0.900, 0.900, 0.900, 0.900, 0.900,0.900,0.900,0.900,0.900,0.900];
x_He=[0.515, 0.490, 0.515, 0.490, 0.540, 0.515, 0.515,0.515, 0.515,0.515,0.515,0.515,0.515,0.515,0.515];
SR=[ 0, 0, 0, 0, 0, 0.0,0.0,0.0,0.0,0.0,0.0, 0, 0.0,0.0,0.0];

expo=E_o*(1-Theta);
S=I_o.*Theta^(2.5).*exp(expo);

alpha_noshape_mat=transpose(S)*nones;


%df aka gamma delta-nu change aka pressure broadening term
df=w_s.*P_h2o.*Theta.^(x_s)+w_H2.*P_H2.*Theta.^(x_H2)+w_He.*P_He.*Theta.^(x_He);


% calculate van-vleck sum it up...
F = vvwlineshape_modified(f,transpose(f_o),transpose(df),transpose(SR));
vvw_alpha_matrix=alpha_noshape_mat.*F;
line_contribution=inv_km_to_dB*convert_to_km*density_h2o*sum(vvw_alpha_matrix,1);
sum(vvw_alpha_matrix,1)

%Continuum Terms Foreign, and self 
%Cf=((1/mbars_to_bars)^2) * 5.43e-10; %  (dB/km)/((GHz x kPa)^2)->db/km((GHz bars)^2)  %eqn 10 Rosenkranz,1998
%Cs=((1/mbars_to_bars)^2)* 1.8e-8*Theta.^(4.5);    %equation 5 (correction applied from ^+6,^-6) Rosenkranz,1998

%Foreign_Continuum=Cf.*P_f.*P_h2o.*f.^2*Theta^3; %foreign continuum term from Eqn 6, Rosenkranz
%Self_Continuum=Cs*P_h2o^2.*f.^2*Theta^3; %self continuum term from Eqn 6, Rosenkranz

Cf_He=((1.0/mbars_to_bars)^2)*1.03562010226e-10;%  (dB/km)/((GHz x kPa)^2)->db/km((GHz bars)^2)  #eqn 10 Rosenkranz,1998
Cf_H2=((1.0/mbars_to_bars)^2) *5.07722009423e-11;
           
Cs1=  4.36510480961e-07*power(Theta,13.3619799812);%equation 5 (correction applied from ^+6,^-6) Rosenkranz,1998
Cs2=  2.10003048186e-26 *((1.0/mbars_to_bars)*P_h2o)^(6.76418487001)*power(Theta,0.0435525417274);  
P_f=P_He+P_H2;

Foreign_Continuum_He=Cf_He*P_He*P_h2o*power(f,2)*power(Theta,3); %foreign continuum term from Eqn 6, Rosenkranz
Foreign_Continuum_H2=Cf_H2*P_H2*P_h2o*power(f,2)*power(Theta,3);
Foreign_Continuum=Foreign_Continuum_He+Foreign_Continuum_H2;

Self_Continuum=Cs1*((P_h2o*(1.0/mbars_to_bars))^2)*power(f,2.0)+Cs2*power(f,2.0);%*power(Theta,3) #self continuum term from Eqn 6, Rosenkranz

alpha_h2o=line_contribution+inv_km_to_dB*Foreign_Continuum+inv_km_to_dB*Self_Continuum;
%alpha_h2o=alpha_h2o/434294.5;

