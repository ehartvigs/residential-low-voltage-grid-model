function [VoltLimit CustomersPerArea FuseLimit type fuse CableSize z_loop AVG_LoadProfile PowerDemand AVGEnergy CustomersPerTransformer TrCap...
    CustomersCalc CustomersInitial VoltageLower LLVMax solcap CustomersPerKm NumberOfTransformers NumberOfCustomers Likelihood Limiter LikelihoodTr LikelihoodVL LikelihoodVL2 LikelihoodVU LikelihoodCL Likelihood11 voltage]...
    = network_model_SWE_EV_2022(factor, Pop_density, PeoplePerHouse,HH_LoadProfile, AP_LoadProfile1, Rgrid, Xgrid, no_load, thermal_limit, alpha, voltageLimit, CarsPerHH)
% -------------------------------------------------------------------------
% This is the modelling file that generates the low-voltage grids, and then
% populate them with data and check for violations. The bulk of the file
% concerns the generation of LV-grids, and only from about row 520 and
% onwards is about population and simulating the grids.
% The file is a bit messy when it comes to units for impedance (R, X and Z)
% as these are normally expressed per km, but the model us m as unit. This
% is seen by the occasional 1000 multiplier or divider.
% The file is divided into four sections:
% Section 1: Demand estimation
% Section 2: Feeder length calculation
% Section 3: Feeder sizing, capacity, voltage and tripping criteria
% Section 4: Adding new loads (EV, solar PV etc.)
% -------------------------------------------------------------------------









% ------------------------------------------------------------------------
%                                  SECTION 1
% ------------------------------------------------------------------------

%% ---------------------------------------------------------
% The use of global coincidence file is specifically for EV
global CoincidenceEVLine
global coincidenceEVTR

% SET VALUE OF LAMBDA MATRIX
% 50 for small coincidence matrix (to run on Therese laptop, otherwise 100)
lambda = 50;

%% ---------------------------------------------------------

VoltLimit = 10;
FuseLimit = 0;

Design_Z = 0.65;    % Maximum earth loop impedance in Swedish distribution grids
Design_voltage = [1.06 0.92];   % Design voltage

% Set feeder parameters based on demographic type
if Pop_density<=200
    LLAF_LV = 1;
    p = 1;
    type = 1;
elseif Pop_density>200 && Pop_density<=1000
    LLAF_LV = 1.1;
    p = 0;
    type = 2;
else
    LLAF_LV = 1.2;
    p = 1;
    type = 3;
end


Limiter = 0;
margin=alpha;
AptBuilding_share=1.25;     % https://www.sciencedirect.com/science/article/pii/S0378778811005019#fig0010
AVGEnergyHouse = factor*18.5*1000;     % Average annual electricity consumption for house in Sweden, Zimmermann
%AVGEnergyHouse = factor*8.5*1000;     % Average annual electricity consumption for house with district heating
AVGEnergyApt = factor*3.71*1000;       % Zimmermann, 6 high
solcap = 0;


AVGEnergyHouseFuture = AVGEnergyHouse;
AVGEnergyAptFuture = AVGEnergyApt*AptBuilding_share;

% Velander coefficients
k2_House = (0.025+0.05)/2;
k1_House = 0.0003;
k2_Apt = 0.014;
k1_Apt = 0.000264;

% Alternative parameter values for velander coefficients, first electric 
% heating and then for district heating.
% k2_House = 0.05;
% k1_House = 0.0003;
% k2_Apt = 0.04383;
% k1_Apt = 0.00026;

% k1_House = 0.0002;
% k2_House = 0.07;
% k1_Apt = 0.0002;
% k2_Apt = 0.07;

TransformerCap = [1250 800 500 315 200 100 50];
TransformerCost = [195272 134751 101565 70501 53509 38446 32140];

% Set power factor to 0.9 and calculate the corresponding active power
% share (also named PowerFactor).
PowerFactor = 0.9;
PowerFactor = sind(acosd(PowerFactor));

% Define different demographic areas and cable costs for each area. Costs
% are taken from normvardeslistan from EI.se
if Pop_density>1000    % Stadsomrades
    CostPerKmLineLV = 827000;
    CostPerKmLineMV = 1140746;
elseif Pop_density>200  % Tatort
    CostPerKmLineMV = 887790;
    CostPerKmLineLV = 540000;
else        % Landsbygd
    CostPerKmLineMV = 512337;
    CostPerKmLineLV = 177000;
end

% Estimate the share of apartments in each cell based on population
% density.
if Pop_density <= 4000
    frac_AP = min([(-2.5797*10^(-8)*Pop_density^2 + 0.000257*Pop_density + 0.3782) 1]);
    frac_AP = max([frac_AP 0]);
    frac_HH = 1-frac_AP;
else
    frac_AP = 1;
    frac_HH = 0;
end

fuse = 10*frac_AP+20*frac_HH;

CustomersConmectionAP = 1;
CustomersConnectionHH = 1;
CustomersConnectionTot = CustomersConmectionAP*frac_AP + CustomersConnectionHH*frac_HH;

% Persons per dwelling.
PeoplePerHouse = 2.7;
PeoplePerApt = 2;
PeoplePerHousehold = PeoplePerHouse*frac_HH + PeoplePerApt*frac_AP;

if Pop_density<4
    NumberOfCustomers = Pop_density;
    frac_AP = 0;
    frac_HH = 1;
else
    NumberOfCustomers = round(Pop_density/PeoplePerHousehold);
end


% Construct average loadprofile.
AVG_LoadProfile = (frac_HH.*HH_LoadProfile+frac_AP.*AP_LoadProfile1);  % Per connection
PowerDemand = max(AVG_LoadProfile);

AVGEnergy = (frac_AP*AVGEnergyAptFuture + frac_HH*AVGEnergyHouseFuture);

CustomersPerTransformer=NumberOfCustomers;
PeakLoad =margin*(k1_House*CustomersPerTransformer*frac_HH*AVGEnergyHouseFuture+k2_House*sqrt(CustomersPerTransformer*frac_HH*AVGEnergyHouseFuture)) + (k1_Apt*CustomersPerTransformer*frac_AP*AVGEnergyAptFuture+k2_Apt*sqrt(CustomersPerTransformer*frac_AP*AVGEnergyAptFuture));
n0 = ceil(PeakLoad/TransformerCap(1));











% ------------------------------------------------------------------------
%                                  SECTION 2
% ------------------------------------------------------------------------

% For loop iterates through all possible number of transformers that
% supply each cell and calculates the cost (power lines and transformers)
% in each case. Then find the cheapest option which is finally used.
ToTCost = zeros((200),length(TransformerCap));

for gg = n0:round(NumberOfCustomers)
    
    CustomersPerTransformer = round(NumberOfCustomers/gg);
    coincidenceTR = (k1_House*CustomersPerTransformer*AVGEnergy+k2_House*sqrt(CustomersPerTransformer*AVGEnergy))/(CustomersPerTransformer*(k1_House*AVGEnergy+k2_House*sqrt(AVGEnergy)));
    PeakLoadTR = margin*((k1_House*CustomersPerTransformer*frac_HH*AVGEnergyHouseFuture+k2_House*sqrt(CustomersPerTransformer*frac_HH*AVGEnergyHouseFuture)) + (k1_Apt*CustomersPerTransformer*frac_AP*AVGEnergyAptFuture+k2_Apt*sqrt(CustomersPerTransformer*frac_AP*AVGEnergyAptFuture)));
    
    
    indexTR = find((PeakLoadTR./(TransformerCap)<1)~=0,1, 'last');
    TransformerCostArea = TransformerCost(indexTR)*gg;
    
    
    A_TS = 1/gg; % in km^2
    NubmerOfConnectionsTransformer = round(NumberOfCustomers/gg);
    
    ConnectionPointDensity = NubmerOfConnectionsTransformer/A_TS;
    
    % Calcualte distance between low (d) and medium (dMV)
    % voltage supply points
    d = sqrt(A_TS)/(sqrt(NubmerOfConnectionsTransformer)+1)+0.02;
    
    dMV = 1/(sqrt(gg)+1);       % Length of medium voltage cables
    if mod(round(sqrt(NubmerOfConnectionsTransformer)),2)
        n = NubmerOfConnectionsTransformer-1;
    else
        n = NubmerOfConnectionsTransformer + sqrt(NubmerOfConnectionsTransformer) -2;
    end
    
    if round(sqrt(gg)) == 1
        nMV = 0.5;
    elseif mod(round(sqrt(gg)),2)
        nMV = gg-1;
    else
        nMV = gg + sqrt(gg) -2;
    end
 
    LengthLVPerTransformer = n*d;
    LengthMV = nMV*dMV;
    MV_cost = LengthMV*CostPerKmLineMV;         % Medium-voltave costs
    LV_cost = gg*LengthLVPerTransformer*CostPerKmLineLV;        % Low-voltage costs
    
    LineCost = LV_cost + MV_cost;
    ToTCost((gg-n0+1),indexTR) = TransformerCostArea + LineCost;

end

% Select option with lowest cost.
ToTCost(ToTCost==0)=NaN;
MinCost=min(min(ToTCost));
[xt,TrSize]=find(ToTCost==MinCost);
xt=xt(1);
TrSize=TrSize(1);
NumberOfTransformers = n0+xt-1;

TrCap = TransformerCap(TrSize);
A_TS = 1/NumberOfTransformers; % in km^2
w = sqrt(A_TS);

CustomersPerTransformer= round(NumberOfCustomers/NumberOfTransformers);

% Coincidence for each transformer.
coincidenceTR = (k1_House*CustomersPerTransformer*AVGEnergy+k2_House*sqrt(CustomersPerTransformer*AVGEnergy))/(CustomersPerTransformer*(k1_House*AVGEnergy+k2_House*sqrt(AVGEnergy)));


d = LLAF_LV*sqrt(A_TS)/(sqrt(CustomersPerTransformer));
LengthLVPerTransformer = CustomersPerTransformer*d;


% Calculate the length of the longest feeder. The first two (customers = 1
% or 2) are special cases.
if CustomersPerTransformer == 1
    LLVMax = 0.1;
    d = LLVMax;
    n = 1;
    nc = 1;
elseif CustomersPerTransformer == 2
    LLVMax =  LLAF_LV*sqrt(A_TS)/4 + 0.02;
    d = LLVMax;
    n = 2;
    nc = 2;
else
    n = round(sqrt(CustomersPerTransformer));
    nc = round(CustomersPerTransformer/4);
    if n == 2
        d = LLAF_LV*sqrt(A_TS)/(n);
        LLVMax = d;
        nc = 2;
    else
        LLVMax = LLAF_LV*sqrt(A_TS)*(n-1)/(n)+0.02;
    end
end








% ------------------------------------------------------------------------
%                                  SECTION 3
% ------------------------------------------------------------------------

c = 0.85;
Ufn = 230;

% Set data for fuse sizes, cable impedance and transformer impedance.
TransformerType = TrCap;
Transformer_Impedance = [6.5 10 13 20 32 65 130]/1000;
CableCapacity = [52 67 80 94 138 178 230];
fuses_tripping_size = [6 10 16 20 25 32 40 50 63 80 100 125 160 200 230;... % 250 315;...     % Table 3 in EN60269, for fuses less than 16A from BS3036 (UK standard)
    18 32 65 85 110 150 190 250 320 425 580 715 950 1250 1490]; % 1650 2200];               % 1490 through linear interpolation
targetFuseRatings = [6 10 16 20 25 32 40 50 63 80 100 125 160 200 230]; % 250 315];
Z_lineR_list = [1.83 1.15 0.76 0.641 0.32 0.206 0.125];
Z_lineX_list = [0.0861 0.0817 0.0783 0.0779 0.0762 0.0745 0.0752];

Z_TransformerR_list = [0.000848114120254821 0.00148841681507050 0.00251028652976876 0.00492770793724613 0.00776114000116266 0.0202385770250776 0.0404771540501553];
Z_TransformerX_list = [0.00763302708229339 0.0119073345205640 0.0125514326488438 0.0197108317489845 0.0310445600046506 0.0607157310752329 0.121431462150466];
Transformer_R = Z_TransformerR_list(TrSize);
Transformer_X = Z_TransformerX_list(TrSize);

Z_transformer = Transformer_Impedance(TrSize);

% Set number of brances based on number of customers per low-voltage grid.
if n > 10
    m = zeros(1,5);
    mc = zeros(1,5);
    m(1)= round((n-1)/5);
    m(2) = m(1);
    m(3) = m(1);
    m(4) = m(1);
    m(5) = n-m(1)-m(2)-m(3)-m(4);
    mc(1)= round((nc-1)/5);
    mc(2) = mc(1);
    mc(3) = mc(1);
    mc(4) = mc(1);
    mc(5) = nc-mc(1)-mc(2)-mc(3)-mc(4);
    p = 5;
elseif n > 8
    m = zeros(1,4);
    mc = zeros(1,4);
    m(1)= round((n-1)/4);
    m(2) = m(1);
    m(3) = m(1);
    m(4) = n-m(1)-m(2)-m(3);
    mc(1)= round((nc-1)/4);
    mc(2) = mc(1);
    mc(3) = mc(1);
    mc(4) = nc-mc(1)-mc(2)-mc(3);
    p = 4;
elseif n>6
    m = zeros(1,3);
    mc = zeros(1,3);
    m(1)= round((n-1)/3);
    m(2) = m(1);
    m(3) = n-m(1)-m(2);
    mc(1)= round((nc-1)/3);
    mc(2) = mc(1);
    mc(3) = nc-mc(1)-mc(2);
    p = 3;
elseif n == 1 || n == 2
    m = zeros(1,1);
    mc = zeros(1,1);
    m(1) = n;
    mc(1) = nc;
    p = 1;
else
    m = zeros(1,2);
    mc = zeros(1,2);
    m(1)= round((n-1)/2);
    m(2) = n-m(1);
    mc(1)= round((nc-1)/2);
    mc(2) = nc-mc(1);
    p = 2;
end

% Flip branch vector (so first element represents the branch closest to the
% transformer.
mc = flip(mc);
d_long = zeros(1,p);

for uu = 1 :p
    if n == 1
        d_long(uu) = (m(uu))*d*1000;
    else
        d_long(uu) = (m(uu)-1/p)*d*1000;
    end
end
d_long = flip(d_long);

% Create variables that are used for cable sizing.
Cable = zeros(1,p);
ixCable = [];
R = zeros(1,p);
X = zeros(1,p);
ZmaxThick = zeros(1,p);
coincidenceLine = zeros(1,p);
P_demand = zeros(1,p);
Z = zeros(1,p);
mp = zeros(1,p);
z_earth = zeros(1,p+1);
RX_multiplier = ones(1,p);
L_max = zeros(7,p);
I_line_fuse = zeros(1,p);
targetLineImpedance = [4.18 2.63 1.91 1.47 0.746 0.495 0.324];
Z_Line = targetLineImpedance;
Z_loop = [Z_transformer*1000];


%% ----------------------------------------------------------------------------------------
% Below code size cables according to tripping critera, voltage limits and thermal limits.
%% ----------------------------------------------------------------------------------------

% Loop over each branch (p) to check fuse size, and then if tripping
% critera is fullfilled.
for rt=1:p
    rr = sum(mc(rt:end));
    mp(rt) = rr;  % store cumulative customers in mp (and don't update it with rr anymore)
    coincidenceLine(rt) = (k1_House.*rr.*AVGEnergy+k2_House.*sqrt(rr.*AVGEnergy))./(rr.*(k1_House*AVGEnergy+k2_House*sqrt(AVGEnergy))); % Velanders formula
    P_demand(rt) = max(AVG_LoadProfile)*rr*coincidenceLine(rt);
    I_line = coincidenceLine(rt).*fuse*rr;
    
    % Check if cable current is compatible with fuse size. Otherwise
    % reduce demand through RX_multiplier.
    while I_line>max(targetFuseRatings)
        rr = round(rr/2);
        % RX_multiplier is a factor that consider multiple cables laid in parallell if needed.
        % If it is 0.5, then two cables are used, 0.25 four cables etc. 
        RX_multiplier(rt) =  RX_multiplier(rt)*0.5;     
        coincidenceLine(rt) = (k1_House.*rr.*AVGEnergy+k2_House.*sqrt(rr.*AVGEnergy))./(rr.*(k1_House*AVGEnergy+k2_House*sqrt(AVGEnergy))); % Velanders formula
        P_demand(rt) = max(AVG_LoadProfile)*rr*coincidenceLine(rt);
        I_line = coincidenceLine(rt).*fuse*rr;
        
    end
    
    I_line_fuse(rt) = targetFuseRatings(find((ceil(I_line)./targetFuseRatings<=1)~=0,1, 'first'));
    index = find((I_line_fuse(rt)./targetFuseRatings<=1)~=0,1, 'first') ;
    Iu = fuses_tripping_size(2,(index));
    ZmaxThick(rt) = 1000*c*Ufn/Iu;
    
    % Check if the tripping critera is fullfilled given length and capacity
    % of cables and fuse ratings. If needed use multiple parallell cables
    % (e.g. change RX_multiplier)
    L_max(:,rt) = (ZmaxThick(rt)-sum(Z_loop))./targetLineImpedance';
    while sum(L_max(:,rt)>d_long(rt)) == 0
        rr = round(rr/2);
        RX_multiplier(rt) = RX_multiplier(rt) * 0.5;
        coincidenceLine(rt) = (k1_House.*rr.*AVGEnergy+k2_House.*sqrt(rr.*AVGEnergy))./(rr.*(k1_House*AVGEnergy+k2_House*sqrt(AVGEnergy))); % Velanders formula
        P_demand(rt) = max(AVG_LoadProfile)*rr*coincidenceLine(rt);
        I_line = coincidenceLine(rt).*fuse*rr;
        I_line_fuse(rt) = targetFuseRatings(find((ceil(I_line)./targetFuseRatings<1)~=0,1, 'first'));
        index = find((I_line_fuse(rt)./targetFuseRatings<=1)~=0,1, 'first') ;
        Iu = fuses_tripping_size(2,(index));
        ZmaxThick(rt) = 1000*c*Ufn/Iu;
        L_max(:,rt) = (ZmaxThick(rt)-sum(Z_loop))./targetLineImpedance';
    end
    
    
    
    L_max(L_max<d_long(rt))=NaN;
    L_max(CableCapacity<I_line_fuse(rt),rt) = NaN;
    [void ixZloop] = min(L_max(:,rt));
    Z_loopNew = targetLineImpedance(ixZloop)*d_long(rt); % Z_loop(ixZloop,rt);
    Z_loop = [Z_loop Z_loopNew];
end



L_max(isnan(L_max))=0;
z_earth(1) = Z_transformer*1000;

% Check if cable capacity is sufficient and calcualte earth impedance.
% Earth impedance is used for validation purposes.
for gg=1:p
    ixCable(gg) = find(L_max(:,gg),1);
    Cable(gg) = CableCapacity(ixCable(gg));
    if CableCapacity(ixCable(gg))<I_line_fuse(gg)
        if I_line_fuse(gg) == 250 || I_line_fuse(gg) == 315
            break
        else
            ixCable(gg) = find((I_line_fuse(gg)./CableCapacity<=1)~=0,1, 'first');
            Cable(gg) = CableCapacity(ixCable(gg));
        end
    end
    R(gg) = Z_lineR_list(ixCable(gg))*d_long(gg)*RX_multiplier(gg);
    X(gg) = Z_lineX_list(ixCable(gg))*d_long(gg)*RX_multiplier(gg);
    z_earth(gg+1) = z_earth(gg) + targetLineImpedance(ixCable(gg))*d_long(gg)*RX_multiplier(gg);
    
end


% Update cable resistance and reactance if cable capacity has been changed.
% ELIAS: THIS WHOLE LOOP CAN BE COMMENTED OUT
%id = find(ixCable==max(ixCable),1);
%for hh=1:id
%    Cable(hh) = CableCapacity(ixCable(id));
%    R(hh) = Z_lineR_list(ixCable(id))*d_long(hh);
%    X(hh) = Z_lineX_list(ixCable(id))*d_long(hh);
%    z_earth(hh+1) = z_earth(hh) + targetLineImpedance(ixCable(hh))*d_long(hh)*RX_multiplier(hh);
%end


CustomersPerArea = NumberOfCustomers;

% Set correct unit for impedance values (per m, not per km)
R = R./1000;
X = X./1000;

N = sum(n)/sum(mc);

% Nomincal voltage is 400V (low-voltage standard)
Vn = 400;
V = zeros(1,p+1);

% Calcualte the initial voltage drop at the transformer.
V(1) = 400 -coincidenceTR*CustomersPerTransformer*(Transformer_R.*(fuse*230*3)+Transformer_X.*(fuse*230*3)*PowerFactor)./(Vn);
% Calculate voltage profile along the cables.
for kk=1:p
    V(kk+1) = V(kk)-(R(kk).*(P_demand(kk))+(X(kk)).*P_demand(kk)*PowerFactor)./(Vn);
end

% Update the cable capacity if voltage profile is outside limits.
while  min(V./Vn)<Design_voltage(2)
    [val pos] = min(ixCable);
    if val == 7
        [val2 pos] = max(RX_multiplier);
        RX_multiplier(pos) = RX_multiplier(pos)*0.5;
        R(pos) = Z_lineR_list(val)*d_long(pos)*RX_multiplier(pos)/1000;
        X(pos) = Z_lineX_list(val)*d_long(pos)*RX_multiplier(pos)/1000;
        ixCable(pos) = val;
        Cable(pos) = CableCapacity(val);
        z_earth(pos+1) = z_earth(pos) + targetLineImpedance((val))*d_long(pos)*RX_multiplier(pos);
    else
        R(pos) = Z_lineR_list(val+1)*d_long(pos)*RX_multiplier(pos)/1000;
        X(pos) = Z_lineX_list(val+1)*d_long(pos)*RX_multiplier(pos)/1000;
        ixCable(pos) = val +1;
        Cable(pos) = CableCapacity(val+1);
        z_earth(pos+1) = z_earth(pos) + targetLineImpedance((val+1))*d_long(pos)*RX_multiplier(pos);
    end
    
    
    V = zeros(1,p+1);
    V(1) = 400 -min((coincidenceTR*CustomersPerTransformer*(Transformer_R.*(fuse*230*3)+Transformer_X.*(fuse*230*3).*PowerFactor)./Vn));
    for kk=1:p
        V(kk+1) = V(kk)-(R(kk).*(P_demand(kk))+(X(kk)).*P_demand(kk).*PowerFactor)./(Vn);
    end

end


% Checking that loop impedance is lower than maximum value. Max value is
% taken from a CIRED paper written by Per Norborg at Vattenfall.
z_earth = z_earth/1000;     % Setting correct unit for impedance (per m)
while  max(z_earth)>Design_Z
    
    [val pos] = min(ixCable);
    if val == 7
        [val2 pos] = max(RX_multiplier);
        RX_multiplier(pos) = RX_multiplier(pos)*0.5;
        
        for gg = pos:length(ixCable)
            z_earth(gg+1) = z_earth(gg) + Z_Line(val)*d_long(gg)*RX_multiplier(gg)/1000;
        end
        R(pos) = Z_lineR_list(val)*d_long(pos)*RX_multiplier(pos)/1000;
        X(pos) = Z_lineX_list(val)*d_long(pos)*RX_multiplier(pos)/1000;
        ixCable(pos) = val;
        Cable(pos) = CableCapacity(val)*(1/RX_multiplier(pos));
    else
        for gg = pos:length(ixCable)
            z_earth(gg+1) = z_earth(gg) + Z_Line(val+1)*d_long(gg)*RX_multiplier(gg)/1000;
        end
        R(pos) = Z_lineR_list(val+1)*d_long(pos)*RX_multiplier(pos)/1000;
        X(pos) = Z_lineX_list(val+1)*d_long(pos)*RX_multiplier(pos)/1000;
        ixCable(pos) = val+1;
        Cable(pos) = CableCapacity(val+1);
    end
    
    
end
z_loop = repelem(z_earth,NumberOfTransformers);

CableSize = mean(Cable);

% If no load profile is used, e.g. only new technology (EV, PV..)
AVG_LoadProfile = AVG_LoadProfile*no_load;






% ------------------------------------------------------------------------
%                                  SECTION 4
% ------------------------------------------------------------------------

%% ---------- ONLY EDIT BELOW HERE IF ADDING NEW TECHNOLOGIES ---------- %%
% This part of the code concern the implementation/analysis of
% technologies. Currently it is written for EVs. 
% If a new technology is added (solar PV, heat pumps etc., then it is only
% the code below here that needs to be changed. The code contains three
% aspects.
% 1. Create new net demand profile.
% 2. Use the net demand profile to calculate voltage drop and power demand.
% 3. Check if any limits have been reached.
%% ----------------------------------------------------------------------%%


mp2 = round(mp.*RX_multiplier);
ChargePower = 6900;     % Charge power in watt
AdditionalVoltageLimit = 0.92;

  % Initial voltage drop (at the transformer)

    V1 = 400 - (CustomersPerTransformer*(Transformer_R.*(coincidenceTR.*AVG_LoadProfile+...
        CarsPerHH.*ChargePower.*reshape(coincidenceEVTR(1,:,:),[lambda 52560])'))...
        +Transformer_X.*(coincidenceTR.*AVG_LoadProfile))./Vn;
    
    CEVL = reshape(CoincidenceEVLine([ceil(CarsPerHH.*mp2)],:,:),[length(mp) lambda 52560]);
    CEVL = permute(CEVL,[3 1 2]);
    

    voltage1 = V1 - squeeze(sum((mp.*(R.*(coincidenceLine.*AVG_LoadProfile + CarsPerHH.*CEVL.*ChargePower)+...
                X.*(coincidenceLine.*AVG_LoadProfile*PowerFactor))./Vn),2));
    
    deltaCurrentCable1 = (mp2.*abs(repelem(coincidenceLine.*AVG_LoadProfile,1,1,lambda) + ...
        CarsPerHH.*CEVL.*ChargePower)/400/sqrt(3));
  
    deltaPower1 = ((coincidenceTR.*AVG_LoadProfile + ...
            CarsPerHH.*(reshape(coincidenceEVTR(1,:,:),[lambda 52560])').*ChargePower)).*CustomersPerTransformer;
    Cable1=repmat(Cable,52560,1,lambda);

        
    VoltageLowerLimit1 = (voltage1./400)<voltageLimit(2);
    VoltageLowerLimit2 = (voltage1./400)<AdditionalVoltageLimit;
    VoltageUpperLimit1 = (voltage1./400)>voltageLimit(1);
    TransformerLimit1 = deltaPower1>(TransformerType*1000)*thermal_limit;  
    CableLimit1 = deltaCurrentCable1>Cable1;
    CableLimit1 = squeeze(max(CableLimit1,[],2));


    VoltageRange = max(voltage1,[],'all') - min(voltage1,[],'all');
    VoltageLower = min(voltage1,[],'all');

    

    % The block of code below calculates how much extra power the violation 
    % require (e.g. what is needed to avoid them), and how long duration this
    % occurs over (in 10 min blocks). This can be used to identify
    % reinforcements. Cable demand is an vector while transformer demand is a
    % value. If required capacity to avoid violations is not wanted, the code
    % can be commented out. NOTE!! VOLTAGE VIOLATIONS ARE NOT COVERED.
%     TransformerDemand = max([TransformerType*(deltaPower./max((TransformerType*1000)*thermal_limit)-1) 0]);      % Unit, power (kVA)
%         TransformerDemand2 = ([TransformerType*(deltaPower./max((TransformerType*1000)*thermal_limit)-1) 0]);      % Unit, power (kVA)
% 
%     deltaPowerSum = max(sum(((coincidenceTR.*AVG_LoadProfile + ...
%          CarsPerHH.*(reshape(coincidenceEVTR(1,:,:),[lambda 52560])').*ChargePower).*CustomersPerTransformer)>(TransformerType*1000*thermal_limit)));     % Unit, time (number of 10 min blocks in one year)
%     CableDemand = sqrt(3)*400*Cable.*(max(deltaCurrentCable'./Cable) - 1);
%     CableDemand = CableDemand.*(CableDemand>0);          % Unit, power (W)
%     CableDemandSum = sum(deltaCurrentCable'>Cable);      % Unit, time (number of 10 min blocks in one year)


    ViolationMatrix1 = cat(3,VoltageLowerLimit1, VoltageUpperLimit1, TransformerLimit1, CableLimit1);
    ViolationMatrix1 = max(ViolationMatrix1,[],3);

    LikelihoodAll = sum(ViolationMatrix1~=0,2)./lambda;
    LikelihoodVLt = sum(VoltageLowerLimit1~=0,2)./lambda;
    LikelihoodVL2t = sum(VoltageLowerLimit2~=0,2)./lambda;
    LikelihoodVUt = sum(VoltageUpperLimit1~=0,2)./lambda;
    LikelihoodTLt = sum(TransformerLimit1~=0,2)./lambda;
    LikelihoodCLt = sum(CableLimit1~=0,2)./lambda;

    Likelihood11=sum(LikelihoodAll);

    LikelihoodVL = sum(VoltageLowerLimit1,"all")/lambda;
    LikelihoodVL2 = sum(VoltageLowerLimit2,"all")/lambda;
    LikelihoodVU = sum(VoltageUpperLimit1,"all")/lambda;
    LikelihoodCL = sum(CableLimit1,"all")/lambda;
    LikelihoodTr = sum(TransformerLimit1,"all")/lambda;

    LikelihoodTime = sum(ViolationMatrix1~=0,2)./lambda;
    Likelihood = sum(sum(ViolationMatrix1~=0,1)~=0)/lambda;  %kolla riktning på denna!!

    % Limiter: 1 for no violation, 2 or 3 for voltage, 4 for transformer
    % and 5 for cable. 
%    [tmp Limiter] = max([0 idxA idxB idxC idxD]);

    
    CustomersPerKm = CustomersPerTransformer/LengthLVPerTransformer;
    CustomersCalc = CustomersPerTransformer*NumberOfTransformers;
    CustomersInitial = Pop_density/PeoplePerHousehold;
    
    voltage = 0;
%    voltage = voltage1;
    



end


