%function [CapPerCustomer Energy HHNumber AptNumber cap] = Main_Model_SWE
%index = 1;
 clear all
 clc
 close all
tic
% SWEDEN
filename = 'Inspire_TotPop_Sweref_regionAptShare.shp';

  AvgLLVMAX = [];
    AvgTrCap = [];
    AvgTr = [];
PopData=[272, 365, 495, 1242];

%for aa = 1:length(PopData)
MeanTrCap = [];
LLVMAX_save = [];
% for mn = 1:length(PopData)
% S = shaperead(filename,'Selector',{@(v1) (v1<(PopData(mn)+1)) && (v1>(PopData(mn)-1) ), 'Pop'});
S = shaperead(filename,'Selector',{@(v1) (v1<(500000)) && (v1>(10000) ), 'Pop'});
%S = shaperead(filename);

% Sweden
% DemandHH=load('houseData.mat');
% DemandAP=load('aptData.mat');
Demand = load('ResidentialLP.mat');
DemandHH = Demand.HHLoadProfile;
DemandAP = Demand.APTLoadProfile;


% Randomly assign load profile
HHNumber = randi(20);
AptNumber = randi(15);
Load = DemandHH(HHNumber,:);
LoadAP = DemandAP(AptNumber,:);

% DemandHH=load('HouseLoadProfile.mat');
% DemandAP=load('AptLoadProfile.mat');
% 
% % For usign hourly load profiles
% load('HourlyLP_SWE.mat');
% House = repelem(House,[6]);
% Apt = repelem(Apt,[6]);
% Apt(find(Apt==max(Apt)))=0;
% House(find(House==max(House)))=0;
% House = House;
% 
% % Sweden
% %LoadProfileHH = (DemandHH.data4 + DemandHH.data5 + DemandHH.data2 + DemandHH.data3).*6.*1000./4;
% % data2 has highest peak load
% 
% LoadProfileHHMedium = DemandHH.house1tmp';
% max(LoadProfileHHMedium);
% LoadProfileHH0 = 0.*DemandHH.house5tmp'; %zeros(length(LoadProfileHHMedium),1);
% LoadProfileHHHigh = DemandHH.house20tmp'; %data2 high power
% max(LoadProfileHHHigh);
% Load = [LoadProfileHH0 LoadProfileHHMedium LoadProfileHHHigh];
% %Load = [LoadProfileHH0 House(1:52560)' 1.3.*House(1:52560)'];   % Hourly resultion
% %LoadProfileAG = LoadProfileHHelec;
% %LoadProfileAP = (DemandAP.data1+DemandAP.data3+DemandAP.data4+DemandAP.data5).*6.*1000/4;    % har enheten kwh/10min, d�rf�r mult. med 6
% % data3 has highest peak load
% % LoadProfileAP12 = DemandAP.data4*6.*1000;
% % max(LoadProfileAP12)
% LoadProfileAPMedium = DemandAP.apt8tmp';
% max(LoadProfileAPMedium);
% LoadProfileAPHigh = DemandAP.apt3tmp';
% max(LoadProfileAPHigh);
% LoadProfileAP0 = 0.*DemandAP.apt3tmp';
% LoadAP = [LoadProfileAP0 LoadProfileAPMedium LoadProfileAPHigh];
% %LoadAP = [LoadProfileAP0 Apt(1:52560)' 1.3.*Apt(1:52560)'];    % Hourly resultion



%kk=load('Solar_SWE.mat');
kk=load('PV_factor_XYcoordinatesSWE.mat');
PV_factor = kk.PV_factor;
PV_X = kk.X;
PV_Y = kk.Y;
%Sol = kk.Solar_Insolation;
%Sol = Sol(7:52420)./1000;       % i kW

% Price = load('SpotPrice.mat');

% SpotPrice = Price.SpotPrice;
% SpotPrice = SpotPrice(1:8760)./1000;
% SpotPrice = SpotPrice';

%shape = geoshape(shaperead(filename,'Selector',{@(v1) v1>500, 'Pop'}));
%mapshow(S,'FaceColor',[0.1 0.1 0.8],'edgecolor','none');
%hold on
%mapshow(S,'FaceColor',[0.1 0.1 0.8],'edgecolor','none');
%hold on

% X=extractfield(S2,'X');
% X=vec2mat(X,length(S2(103).X));
% Y=extractfield(S2,'Y');
% Y=vec2mat(Y,length(S2(103).Y));
x=extractfield(S,'X');
x=vec2mat(x,6);
y=extractfield(S,'Y');
y=vec2mat(y,6);
[l tmp] = size(S);

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




%Zgrid = [0.3 0.55];
%RXratio = [0.5 5];

pp = 3;     % Voltage
%hh = 3;
mm = 3;     % Transformer margin
nn = 2;     % Load


%%%%% Control parameters
Zgrid = [0.3 0.4 0.55];
RXratio = [0.5 2 5];

voltageLimitMatrix = [1.1 0.9; 1.06 0.92; 1.05 0.95];
thermallimit_out = [0.8 1 1.2];


% Transformer sizing
alphavec = [1.2 1.5 1.8 2.1 2.4];
% Load profile (see above)
% Gamma
%gamma = [1 1.1 1.2; 1 1 1; 1 1.2 1.3];
%fuse = [10 20; 6 16; 16 25];
noload = 1;
EnergyScaleFactor = [1 1 1.4];

datastore = [];

 %for pp=1:3
      %for nn=1:3      %Load
   %for mm=1:5      % Transformer sizing
    %pp = 1;
    %for hh=1:3
    %hh=pp;
    deltaV = [];
    llvmax = [];
    CapPerCustomer = [];
    cap = [];
    CustPerKm = [];
    Transformers = [];
    TransformersCap = [];
    Pop_Density = [];
    EnergyPerkW = [];
    Energy = [];
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
    PeakLoad = [];
    CustomersPerFeeder = [];
    ShareAPT = [];
  
    Production = zeros(1,52414);
    Xgrid = Zgrid(2)/sqrt(1+RXratio(2)^2);
    Rgrid = Xgrid*RXratio(2);
    factor = EnergyScaleFactor(2);
    
    thermallimit = thermallimit_out(2);
    LoadProfileHH = Load(:,nn);
    LoadProfileAP1 = LoadAP(:,nn);
    voltageLimit = voltageLimitMatrix(pp,:);
    alpha = alphavec(mm);
    
    %cableLengths = load2length(Zgrid(pp),RXratio(hh));
    % [268 476 1714 2239];  [362 537 811 838]; [3 3 5 5 9]; [10 25 26 42 78];
    %popd = [44 51 51 52 58];
   % [33941455.7699153;33941455.7699153;33941455.7699153]

 for k=1:l       % number of km^2 with data, 1:l
    % Sweden
    PopDensity = S(k).Pop;
    AptShare = S(k).ApartmentS;
    
    
    % UK
    %PopDensity = S(k).Pop;
    %PopDensity = k;
    %for g=1:290 % number of muncipality, sthlm 103
        x_cord = mean(y(k,1:5));            % Conversion mixup with X and Y coordinates
        y_cord = mean(x(k,1:5));
        [value ii] = min(abs(PV_X(:)-x_cord) + abs(PV_Y(:)-y_cord));
        [row col] = ind2sub(size(PV_X),ii);
        Sol_factor = PV_factor(row,col,:);
        Sol_factor = squeeze(Sol_factor);
        Sol_factor = repelem(Sol_factor,6);     % Six fold each value to get 10 min "resolution" of solar irradiation
        Sol_factor = Sol_factor(1:52560);   
        
        

            CustomersConmectionAPLarge = 40;
            CustomersConnectionAPSmall = 3;
            CustomersConnectionHH = 2;
            CustomersConnectionAG = 1;
            
            %CustomersConnectionTot = (HH*CustomersConnectionHH+AP*(CustomersConmectionAPLarge*APLarge+CustomersConnectionAPSmall*APSmall)+AG*CustomersConnectionAG)
            [Shareapt customersperfeeder MaxLoad Vlimit CustomersPerAreaOut fuselimit type fuseout CableSize z_loop solargeneration AVG_LoadProfile PDemand CustEnergyUsetmp CustomersPerTransformer...
                TrCap CustomersCalcout CustomersInitialout voltage LV ll CustomersPerKm Trans ConnectionDensity CapPerCust EnergykW EnergyPerKM2 Limiter]...
                =network_model_SWE_Final_2020(factor,PopDensity,PeoplePerHH,LoadProfileHH,LoadProfileAP1, Sol_factor, Rgrid, Xgrid, noload, thermallimit, ...
                alpha, voltageLimit,AptShare);
            % [solargeneration AVG_LoadProfile PowerDemand AVGEnergy CustomersPerTransformer TrCap CustomersCalc CustomersInitial voltage LLVMax solcap CustomersPerKm NumberOfTransformers NumberOfCustomers CapPerCustomer EnergyPerkW EnergyPerKM2 Limiter]
            % density
            % break out of for loop
            deltaV = [deltaV voltage];
            llvmax = [llvmax LV];
            cap = [cap ll'];
            Transformers = [Transformers Trans];
            Pop_Density = [Pop_Density PopDensity];
            CustPerKm = [CustPerKm CustomersPerKm];
            CapPerCustomer = [CapPerCustomer CapPerCust];
            EnergyPerkW = [EnergyPerkW EnergykW];
            Energy = [Energy EnergyPerKM2];
            Limit = [Limit Limiter'];
            FuseLimit = [FuseLimit fuselimit];
            CustomersCalc = [CustomersCalc CustomersCalcout];
            CustomersInitial = [CustomersInitial CustomersInitialout];
            TransformersCap = [TransformersCap TrCap];
            CustomersPerFeeder = [CustomersPerFeeder customersperfeeder];
            CustomersPerTr = [CustomersPerTr CustomersPerTransformer];
            CustEnergyUse = [CustEnergyUse CustEnergyUsetmp];
            PowerDemand = [PowerDemand PDemand];
            Z_earth = [Z_earth z_loop];
            Cable = [Cable CableSize];
            PeakLoad = [PeakLoad MaxLoad];
            %Fuse = [Fuse fuseout];
            Production = Production+solargeneration;
            Pop_area = [Pop_area type];
            VoltLimit = [VoltLimit Vlimit];
            ShareAPT = [ShareAPT Shareapt];
            CustomersPerArea = [CustomersPerArea CustomersPerAreaOut];
            PV_efficiency = [PV_efficiency max(Sol_factor)];
                        %S_voltage(k).Pop=max(voltage);
            %break
        %end
        
    %end
 end


 

% AvgLLVMAX = [AvgLLVMAX mean(llvmax)];
% AvgTrCap = [AvgTrCap mean(Transformers.*TransformersCap)];
% AvgTr = [AvgTr mean(Transformers)];

% nn
% CapacitySWE = sum(cap)
% EnergySWE = sum(Energy)
% MeanCap = mean(repelem(CapPerCustomer,CustomersPerArea))
% 
% datastore = [datastore CapacitySWE EnergySWE MeanCap];
   %end

% CapacitySWE = sum(cap)
% EnergySWE = sum(Energy)
% CustomersPerArea = round(CustomersPerArea);
% CapacityPerCustomerSWE = mean(repelem(CapPerCustomer,CustomersPerArea))
% VoltageLimitedSWE = sum(Limit==10)/length(Limit)
% TransformerLimitedSWE = sum(Limit==20)/length(Limit)
% FeederLimitedSWE = sum(Limit==25)/length(Limit)
% Energy = double(Energy);

% Energy = double(Energy);
% ShareApt = double(ShareAPT);
% S_data = S;
% 
% for h = 1:length(cap)
%     S_data(h).Cap = cap(h);
%     S_data(h).CapPerCust = CapPerCustomer(h);
%     S_data(h).Energy = Energy(h);
%     S_data(h).Feeder = llvmax(h);
%     S_data(h).TrCap = TransformersCap(h);
%     S_data(h).TrNumbber = Transformers(h);
%     S_data(h).CustPerFeeder = CustomersPerFeeder(h);
%     S_data(h).CustPerTr = CustomersPerTr(h);
%     S_data(h).Demand = PeakLoad(h);
%     S_data(h).ShareAPt = ShareApt(h);
% end
% 
% 
% %filename = sprintf('%s_%d','SWE_PVLimits',index);
% filename = sprintf('SWE_DataInBrief2');
% shapewrite(S_data,filename);

%  CapPerCustomer = mean(CapPerCustomer);
%  cap = sum(cap);
%  Energy = sum(Energy);

%save('LoopImpedanceSWE.mat','Z_earth','Pop_area')

% save('Data_SWE_V5','Pop_Density','CapPerCustomer','cap')

% 
% 
% EnergyMatrix(pp) = sum(Energy);
% CapMatrix(pp) = sum(cap);
% 
% %PVProduction(pp,:) = Production;
% %ZLoop(pp,:) = Z_earth;
% %CustomersPerAreaFile(pp,:) = CustomersPerArea;
% %AreaType(pp,:) = Pop_area;
% 
% end
%     
% 
% save('EnergyMatrixData.mat','EnergyMatrix')
% save('CapMatrixData.mat','CapMatrix')
% %save('Data','AreaType','CustomersPerAreaFile','ZLoop','PVProduction')


toc


%end


