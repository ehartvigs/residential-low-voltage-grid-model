% Function definition when calling the model
% from a script (to do multiple runs)

%function [FuseLimit Limit Likelihood Pop_Density HHNumber AptNumber xcord
%ycord] = MainModelEV_SWE()   

clear all
clc
close all

tic

% Import GIS data file
filename = 'Inspire_TotPop_Sweref_regionEVPerKom.shp';

AvgLLVMAX = [];
AvgTrCap = [];
AvgTr = [];
PopData=[272, 365, 495, 1242];
MeanTrCap = [];
LLVMAX_save = [];

% Import GIS data in given population density range (v1)
S = shaperead(filename,'Selector',{@(v1) (v1<(80000)) && (v1>(0000) ), 'Pop'});

% Load demand data for households.
Demand = load('ResidentialLP.mat');
DemandHH = Demand.HHLoadProfile;
DemandAP = Demand.APTLoadProfile;

% Randomly assign load profile.
HHNumber = randi(20)
AptNumber = randi(15)
Load = DemandHH(HHNumber,:);
LoadAP = DemandAP(AptNumber,:);


%% Load coincidence data
Obj = load('coincidence50100kWh.mat');

% For MATLAB to access the coincidence data in the network_model file, it
% needs to be loaded as a global file. The downside is that this creates a copy of the file for
% the network_model, and therefore takes a lot of RAM usage.
global CoincidenceEVLine
global coincidenceEVTR

CoincidenceEVLine = Obj.coincidence;
coincidenceEVTR =  Obj.coincidence(398,:,:);



x=extractfield(S,'X');
x=vec2mat(x,6);
y=extractfield(S,'Y');
y=vec2mat(y,6);
[g tmp] = size(S);

PeoplePerHH=2.22;
EnergyMatrix = zeros(3,3,3);
CapMatrix = zeros(3,3,3);
LimitMatrixVoltage = zeros(1,3);
LimitMatrixTransformer  = zeros(1,3);
CapPerCustMatrix = zeros(1,3);
PVProduction = zeros(3,52414);
ZLoop = zeros(3,209989);
CustomerPerAreaFile = zeros(3,103890);
AreaType = zeros(3,103890);



pp = 2;     % Voltage
mm = 3;     % Transformer margin
nn = 2;     % Load

%%%%% Control parameters
Zgrid = [0.3 0.4 0.55];
RXratio = [0.5 2 5];
voltageLimitMatrix = [1.1 0.9; 1.06 0.92; 1.05 0.95];
thermallimit_out = [0.8 1 1.2];

% Transformer sizing
alphavec = [1.2 1.5 1.8 2.1 2.4];
noload = 1;         % Set to 0 to have no residential load profile (only EV charging) when looking for violations
EnergyScaleFactor = [1 1 1.4];

datastore = [];


% Create variables for saving data from the runs.
deltaV = [];
llvmax = [];            % Length of longest feeder
Likelihood = [];        % Share of lambda that produce no violations
cap = [];
CustPerKm = [];
Transformers = [];      % Number of transformers
TransformersCap = [];   % Transformer capacity
Pop_Density = [];
Limit = [];
Z_earth = [];
FuseLimit = [];
PowerDemand = [];
CustomersCalc = [];
CustomersInitial = [];
CustomersPerTr = [];
CustEnergyUse = [];
CustomersPerArea = [];
Cable = [];
Pop_area = [];
PV_efficiency = [];
VoltLimit = [];
FuseLimit = [];
xcord = [];
ycord = [];
LikelihoodTr = [];

Production = zeros(1,52414);
Xgrid = Zgrid(2)/sqrt(1+RXratio(2)^2);
Rgrid = Xgrid*RXratio(2);
factor = EnergyScaleFactor(2);

thermallimit = thermallimit_out(2);
LoadProfileHH = Load';
LoadProfileAP1 = LoadAP';
voltageLimit = voltageLimitMatrix(pp,:);
alpha = alphavec(mm);

toc
for k=1:g     % number of km^2 with data, 1:l
    
    % Extract population density data and number of cars per household.
    PopDensity = S(k).Pop;
    CarsPerHH = S(k).CarsPerHH; % If adding marketshare, add here
    
    % Conversion mixup with X and Y coordinates
    x_cord = mean(y(k,1:5));            
    y_cord = mean(x(k,1:5));
    
    % Call the reference network model
    [Vlimit CustomersPerAreaOut fuselimit type fuseout CableSize z_loop AVG_LoadProfile PDemand CustEnergyUsetmp CustomersPerTransformer...
        TrCap CustomersCalcout CustomersInitialout voltage LV ll CustomersPerKm Trans ConnectionDensity Likeli Limiter LikTr]...
        =network_model_SWE_EV_2022(factor,PopDensity,PeoplePerHH,LoadProfileHH,LoadProfileAP1, Rgrid, Xgrid, noload, thermallimit, ...
        alpha, voltageLimit,CarsPerHH);
    
    
    % Save data from each run
    xcord = [xcord x_cord];
    ycord = [ycord y_cord];
    deltaV = [deltaV voltage];
    llvmax = [llvmax LV];
    cap = [cap ll'];
    Transformers = [Transformers Trans];
    Pop_Density = [Pop_Density PopDensity];
    CustPerKm = [CustPerKm CustomersPerKm];
    Likelihood = [Likelihood Likeli];
    Limit = [Limit Limiter'];
    FuseLimit = [FuseLimit fuselimit];
    CustomersCalc = [CustomersCalc CustomersCalcout];
    CustomersInitial = [CustomersInitial CustomersInitialout];
    TransformersCap = [TransformersCap TrCap];
    CustomersPerTr = [CustomersPerTr CustomersPerTransformer];
    CustEnergyUse = [CustEnergyUse CustEnergyUsetmp];
    PowerDemand = [PowerDemand PDemand];
    Z_earth = [Z_earth z_loop];
    Cable = [Cable CableSize];
    LikelihoodTr = [LikelihoodTr LikTr];
    Pop_area = [Pop_area type];
    VoltLimit = [VoltLimit Vlimit];
    CustomersPerArea = [CustomersPerArea CustomersPerAreaOut];

end

toc

output = sum(Likelihood==0);


%end


