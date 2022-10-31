function [frac_AP CustomersPerFeeder MaxLoad VoltLimit CustomersPerArea FuseLimit type fuse CableSize z_loop solargeneration AVG_LoadProfile PowerDemand AVGEnergy CustomersPerTransformer TrCap...
    CustomersCalc CustomersInitial voltage LLVMax solcap CustomersPerKm NumberOfTransformers NumberOfCustomers CapPerCustomer EnergyPerkW EnergyPerKM2 Limiter]...
    = network_model_SWE_v2(factor, Pop_density, PeoplePerHouse,HH_LoadProfile, AP_LoadProfile1, Sol, Rgrid, Xgrid, no_load, thermal_limit, alpha, voltageLimit,AptShare)


RXratio = [0.5 2 5];
FuseLimit = 0;
Design_Z = 0.65;    % Maximum earth loop impedance in Swedish distribution grids
Design_voltage = [1.06 0.92];   % Design voltage

% Set feeder parameters based on demographic type
if Pop_density<=200
    LLAF_LV = 1;
    p = 1;
    type = 1;
elseif Pop_density>200 && Pop_density<=1000
    LLAF_LV = 1.1;    % Should be 1.1
    p = 0;
    type = 2;
else
    LLAF_LV = 1.2;  % Should be 1.2
    p = 1;
    type = 3;
end


Limiter = 0;
%p = 0;
%margin=1.028^20;     % So customer can upgrade one step on the fuse, before 1.77
margin=alpha; 
AptBuilding_share=1.25;     % https://www.sciencedirect.com/science/article/pii/S0378778811005019#fig0010
AVGEnergyHouse = factor*18.5*1000;     % Average annual electricity consumption for house in Sweden, Zimmermann
%AVGEnergyHouse = factor*8.5*1000;     % Average annual electricity consumption for house with district heating
AVGEnergyApt = factor*3.71*1000;       % Zimmermann, 6 high
%AVGEnergyApt = 6.5*1000;
solcap = 0;

% Comparision with case studies
% Widen


AVGEnergyHouseFuture = AVGEnergyHouse;
AVGEnergyAptFuture = AVGEnergyApt*AptBuilding_share;
% k2_House = 0.025;     % Electric heating
% k2_House = 0.05;        % District heating
% k1_House = 0.0003;
% k2_Apt = k2_House;
% k1_Apt = k1_House;

k2_House = (0.025+0.05)/2;
k1_House = 0.0003;


k2_Apt = 0.014;
k1_Apt = 0.000264;

% k2_House = 0.05;
% k1_House = 0.0003;
% k2_Apt = 0.04383;
% k1_Apt = 0.00026;

TransformerCap = [1250 800 500 315 200 100 50];
TransformerCost = [195272 134751 101565 70501 53509 38446 32140];

% k1_House = 0.0002;
% k2_House = 0.07;
% k1_Apt = 0.0002;
% k2_Apt = 0.07;

%     if Pop_density>1000    % City/Urban
%         CostPerKmLineLV = 827000;
% 
%     elseif Pop_density>200  % Sub-urban
%         CostPerKmLineLV = 540000;
% 
%     else        % Rural
%         CostPerKmLineLV = 177000;

    %end
    % Medium voltage is 24kV, assume 24
V = 400;    % voltage level

PowerFactor = 0.9;
PowerFactor = sind(acosd(PowerFactor));

    if Pop_density>1000    % Stadsomr�de
        CostPerKmLineLV = 827000;
         CostPerKmLineMV = 1140746;
    elseif Pop_density>200  % T�tort
        CostPerKmLineMV = 887790;
        CostPerKmLineLV = 540000;
    else        % Landsbygd
        CostPerKmLineMV = 512337;
        CostPerKmLineLV = 177000;
    end
    
%     if Pop_density <= 4000
%         frac_AP = min([(-2.5797*10^(-8)*Pop_density^2 + 0.000257*Pop_density + 0.3782) 1]);
%         frac_AP = max([frac_AP 0]);
%         frac_HH = 1-frac_AP;
%     else
%         frac_AP = 1;
%         frac_HH = 0;
%     end
    if isnan(AptShare)
        AptShare = 0;
    end
    frac_AP = AptShare;
    frac_HH = 1-frac_AP;
    

fuse = 10*frac_AP+20*frac_HH;

%CustomersConmectionAP = 13/3; % Get from data on average from each country, EU data. Divide with 3 due to each customer being single phase
CustomersConmectionAP = 1; %  Divide with 3 due to each customer being single phase
CustomersConnectionHH = 1;
%CustomersConnectionAG = 1;
CustomersConnectionTot = CustomersConmectionAP*frac_AP + CustomersConnectionHH*frac_HH;
PeoplePerHouse = 2.7;
PeoplePerApt = 2;

PeoplePerHousehold = PeoplePerHouse*frac_HH + PeoplePerApt*frac_AP;
    %PeoplePerHousehold = 2.9;
    
    if Pop_density<4
        NumberOfCustomers = Pop_density;
        %ConnectionDensity = (NumberOfCustomers); % per km^2
        frac_AP = 0;
        frac_HH = 1;
    else
       % NumberOfUnits = Pop_density/PeoplePerHouse;
        %NumberOfAptUnits = NumberOfUnits*frac_AP_old;
        %NumberOfHouseUnits = NumberOfUnits*frac_HH_old;
        %NumberOfCustomers = NumberOfAptUnits + NumberOfHouseUnits;
       % frac_AP = NumberOfAptUnits/NumberOfCustomers;
       % frac_HH = NumberOfHouseUnits/NumberOfCustomers;
       % ConnectionDensity = (NumberOfCustomers)/CustomersConnectionTot; % per km^2
        NumberOfCustomers = round(Pop_density/PeoplePerHousehold);
    end
    
    
Pop_density;
    %Velander_constant_k1 = [0.3 0.27 0.3 0.27 0.27 0.3 0.23 0.20 0.18 0.20 0.18];   % Finnish phd thesis
    % Different consumers (households, single family, two family, row house etc
    %Velander_constant_k2 = [0 0.02 0 0.02 0 0.02 0 0.04 0.04 0.11 0.04 0.11];

    % CONSTRUCT AVG LOADPROFILE
    AP_LoadProfileAVG = (AP_LoadProfile1);      % + AP_LoadProfile2 + AP_LoadProfile3); 
 
    AVG_LoadProfile = (frac_HH.*HH_LoadProfile+frac_AP.*AP_LoadProfileAVG);  % Per connection
    PowerDemand = max(AVG_LoadProfile);
    
   %AVG_Energy = sum(AVG_LoadProfile)/(1000*6);
    AVGEnergy = (frac_AP*AVGEnergyAptFuture + frac_HH*AVGEnergyHouseFuture);
    %Tot_LoadProfile = AVG_LoadProfile.*NumberOfCustomers;   % Add coincidence?
    NumberOfCustomers;
    %coincidenceTR = (k1_House*NumberOfCustomers*AVGEnergy+k2_House*sqrt(NumberOfCustomers*AVGEnergy))/(NumberOfCustomers*(k1_House*AVGEnergy+k2_House*sqrt(AVGEnergy)));
    %if coincidenceTR<0.3
        %coincidenceTR = 0.3;     
    %end
    %(k1_House*AVGEnergyHouseFuture+k2_House*sqrt(AVGEnergyHouseFuture));
    %(k1_Apt*AVGEnergyAptFuture+k2_Apt*sqrt(AVGEnergyAptFuture));
    %PeakLoad = margin*coincidenceTR*max(Tot_LoadProfile)/(margin*coincidenceTR);
        CustomersPerTransformer=NumberOfCustomers;
                %PeakLoad =(k1_House*1*frac_HH*AVGEnergyHouseFuture+k2_House*sqrt(1*frac_HH*AVGEnergyHouseFuture)) + (k1_Apt*1*frac_AP*AVGEnergyAptFuture+k2_Apt*sqrt(1*frac_AP*AVGEnergyAptFuture));

        %coincidenceTR = (k1_House*CustomersPerTransformer*AVGEnergy+k2_House*sqrt(CustomersPerTransformer*AVGEnergy))/(CustomersPerTransformer*(k1_House*AVGEnergy+k2_House*sqrt(AVGEnergy)));
        PeakLoad =margin*(k1_House*CustomersPerTransformer*frac_HH*AVGEnergyHouseFuture+k2_House*sqrt(CustomersPerTransformer*frac_HH*AVGEnergyHouseFuture)) + (k1_Apt*CustomersPerTransformer*frac_AP*AVGEnergyAptFuture+k2_Apt*sqrt(CustomersPerTransformer*frac_AP*AVGEnergyAptFuture));
        n0 = ceil(PeakLoad/TransformerCap(1));
        ToTCost = zeros((200),length(TransformerCap));
        MaxLoad = PeakLoad/margin;
        
   for gg = n0:round(NumberOfCustomers)
       
        CustomersPerTransformer = round(NumberOfCustomers/gg);
        coincidenceTR = (k1_House*CustomersPerTransformer*AVGEnergy+k2_House*sqrt(CustomersPerTransformer*AVGEnergy))/(CustomersPerTransformer*(k1_House*AVGEnergy+k2_House*sqrt(AVGEnergy)));
        PeakLoadTR = margin*((k1_House*CustomersPerTransformer*frac_HH*AVGEnergyHouseFuture+k2_House*sqrt(CustomersPerTransformer*frac_HH*AVGEnergyHouseFuture)) + (k1_Apt*CustomersPerTransformer*frac_AP*AVGEnergyAptFuture+k2_Apt*sqrt(CustomersPerTransformer*frac_AP*AVGEnergyAptFuture)));

        
        indexTR = find((PeakLoadTR./(TransformerCap)<1)~=0,1, 'last');
        TransformerCostArea = TransformerCost(indexTR)*gg;


            A_TS = 1/gg; % in km^2
            NubmerOfConnectionsTransformer = round(NumberOfCustomers/gg); % Number of connections
            %coincidenceTR = (0.0002*NubmerOfConnectionsTransformer*AVGEnergy+0.07*sqrt(NubmerOfConnectionsTransformer*AVGEnergy))/(NubmerOfConnectionsTransformer*(0.0002*AVGEnergy+0.07*sqrt(AVGEnergy)));
            %CustomersPerTransformer= NumberOfCustomers/NumberOfTransformers;

            ConnectionPointDensity = NubmerOfConnectionsTransformer/A_TS;

            % Calcualte distance between low (d) and medium (dMV) 
            % voltage supply points

            d = sqrt(A_TS)/(sqrt(NubmerOfConnectionsTransformer)+1)+0.02;
            dMV = 1/(sqrt(gg)+1);
            %LengthLVPerTransformer = (NubmerOfConnectionsTransformer-1)*DistanceBetweenPoint*LLAF_LV+NubmerOfConnectionsTransformer*L_CP;
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
            gg;
            LengthLVPerTransformer = n*d;
            LengthMV = nMV*dMV;
            MV_cost = LengthMV*CostPerKmLineMV;
            LV_cost = gg*LengthLVPerTransformer*CostPerKmLineLV;
            
            LineCost = LV_cost + MV_cost;
            TransformerCostArea;
TransformerCostArea + LineCost;
    %LineCost = 0;
            ToTCost((gg-n0+1),indexTR) = TransformerCostArea + LineCost; % 0.17 0.08 for city centre, 0.03 for villa areas (215 people/km2)
        %end
        %if CustomersPerTransformer==round(NumberOfCustomers/NumberOfTransformers)
            %break
        %end
        %CustomersPerTransformer= NumberOfCustomers/NumberOfTransformers;
   end
   ToTCost;
   ToTCost(ToTCost==0)=NaN;

    MinCost=min(min(ToTCost));
    size(ToTCost);
    [xt,TrSize]=find(ToTCost==MinCost);
    xt=xt(1);
    TrSize=TrSize(1);
    NumberOfTransformers = n0+xt-1;
    
%     NumberOfTransformers=round(NumberOfTransformers/2);
%     CustomersPerTransformer= round(NumberOfCustomers/NumberOfTransformers);
%     coincidenceTR = (k1_House*CustomersPerTransformer*AVGEnergy+k2_House*sqrt(CustomersPerTransformer*AVGEnergy))/(CustomersPerTransformer*(k1_House*AVGEnergy+k2_House*sqrt(AVGEnergy)));
%     PeakLoadTR = 1*((k1_House*CustomersPerTransformer*frac_HH*AVGEnergyHouseFuture+k2_House*sqrt(CustomersPerTransformer*frac_HH*AVGEnergyHouseFuture)) + (k1_Apt*CustomersPerTransformer*frac_AP*AVGEnergyAptFuture+k2_Apt*sqrt(CustomersPerTransformer*frac_AP*AVGEnergyAptFuture)));
%     indexTR = find((PeakLoadTR./(TransformerCap)<1)~=0,1, 'last');
%     TrCap=TransformerCap(indexTR);
    
    TrCap = TransformerCap(TrSize);
    A_TS = 1/NumberOfTransformers; % in km^2
    w = sqrt(A_TS);
    %NubmerOfConnectionsTransformer = NumberOfCustomers/NumberOfTransformers; % Number of connections
    %NumberOfTransformers;
    %NubmerOfCustomersTransformer = NubmerOfConnectionsTransformer*CustPerConnection
    CustomersPerTransformer= round(NumberOfCustomers/NumberOfTransformers);
    %PeakLoadTR = ((k1_House*CustomersPerTransformer*frac_HH*AVGEnergyHouseFuture+k2_House*sqrt(CustomersPerTransformer*frac_HH*AVGEnergyHouseFuture)) + (k1_Apt*CustomersPerTransformer*frac_AP*AVGEnergyAptFuture+k2_Apt*sqrt(CustomersPerTransformer*frac_AP*AVGEnergyAptFuture)))
    margin;
    %coincidenceTR = PeakLoadTR/(CustomersPerTransformer*230*3*fuse/1000);
    coincidenceTR = (k1_House*CustomersPerTransformer*AVGEnergy+k2_House*sqrt(CustomersPerTransformer*AVGEnergy))/(CustomersPerTransformer*(k1_House*AVGEnergy+k2_House*sqrt(AVGEnergy)));
    
    d = LLAF_LV*sqrt(A_TS)/(sqrt(CustomersPerTransformer));
    LengthLVPerTransformer = CustomersPerTransformer*d;
    %ConnectionPointDensity = NubmerOfConnectionsTransformer/A_TS;

    %DistanceBetweenPoint = sqrt(1/ConnectionPointDensity);

    %LengthLVPerTransformer = (NubmerOfConnectionsTransformer-1)*DistanceBetweenPoint*LLAF_LV+NubmerOfConnectionsTransformer*L_CP;
    %LengthLVPerTransformer = NubmerOfConnectionsTransformer*sqrt(2)*w/(sqrt(NubmerOfConnectionsTransformer)-1);
    %LLVMax = (sqrt(A_TS)-DistanceBetweenPoint)*LLAF_LV+L_CP;
    
    if CustomersPerTransformer == 1
        %LLVMax =  LLAF_LV*sqrt(A_TS)/2 + 0.02;
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
        LLVMax = LLAF_LV*sqrt(A_TS)*(n-1)/(n)+0.02;     % 1.55 is line adjusting factors and is similar to finish guy
        %LLVMax = sqrt(A_TS)*(n-1)/(n+1)+0.02;
        %LLVMax = sqrt(A_TS)*(n)/(n+1)+0.01;
        %ll = LengthLVPerTransformer/NubmerOfConnectionsTransformer;
        end
    end
    TransformerType = TrCap;
    Transformer_Impedance = [6.5 10 13 20 32 65 130]/1000;
    CableCapacity = [52 67 80 94 138 178 230];
    
    Z_transformer = Transformer_Impedance(TrSize);
    Zgrid = sqrt(Rgrid^2 + Xgrid^2);
    Energy = 10000;
    c = 0.85;
    Ufn = 230;
    Max_Power_single = max(AVG_LoadProfile);
    AVG_LoadProfile = AVG_LoadProfile;%./sqrt(3);
    fuses_tripping_size = [6 10 16 20 25 32 40 50 63 80 100 125 160 200 230;... % 250 315;...     % Table 3 in EN60269, for fuses less than 16A from BS3036
          18 32 65 85 110 150 190 250 320 425 580 715 950 1250 1490]; % 1650 2200];               % 1490 through lienar interpolation
    targetFuseRatings = [6 10 16 20 25 32 40 50 63 80 100 125 160 200 230]; % 250 315];
    Z_lineR_list = [1.83 1.15 0.76 0.641 0.32 0.206 0.125];
    Z_lineX_list = [0.0861 0.0817 0.0783 0.0779 0.0762 0.0745 0.0752];
   
    Z_TransformerR_list = [0.000848114120254821 0.00148841681507050 0.00251028652976876 0.00492770793724613 0.00776114000116266 0.0202385770250776 0.0404771540501553];
    Z_TransformerX_list = [0.00763302708229339 0.0119073345205640 0.0125514326488438 0.0197108317489845 0.0310445600046506 0.0607157310752329 0.121431462150466];
    Transformer_R = Z_TransformerR_list(TrSize);
    Transformer_X = Z_TransformerX_list(TrSize);
    
    d;
    LLVMax;
    
    
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
        %Z_loop = zeros(7,p);
        L_max = zeros(7,p);
        I_line_fuse = zeros(1,p);
        targetLineImpedance = [4.18 2.63 1.91 1.47 0.746 0.495 0.324];
         Z_Line = targetLineImpedance; %sqrt(Z_lineR_list.^2+Z_lineX_list.^2);
        Z_transformer;
        status = 0;
        CustomersPerFeeder = sum(mc);

        %Z(1) = 0;
        Z_loop = [Z_transformer*1000];
        
            for rt=1:p
                rr = sum(mc(rt:end));
                mp(rt) = rr;
                coincidenceLine(rt) = (k1_House.*rr.*AVGEnergy+k2_House.*sqrt(rr.*AVGEnergy))./(rr.*(k1_House*AVGEnergy+k2_House*sqrt(AVGEnergy))); % Velanders formula
                P_demand(rt) = max(AVG_LoadProfile)*rr*coincidenceLine(rt);
                %I_line = ((k1_House*rr*frac_HH*AVGEnergyHouseFuture+k2_House*sqrt(rr*frac_HH*AVGEnergyHouseFuture)) + (k1_Apt*rr*frac_AP*AVGEnergyAptFuture+k2_Apt*sqrt(rr*frac_AP*AVGEnergyAptFuture)))*1000/230/3
                I_line = coincidenceLine(rt).*fuse*rr;
                while I_line>max(targetFuseRatings)
                    rr = round(rr/2);
                    RX_multiplier(rt) =  RX_multiplier(rt)*0.5;
                    %mp(rt) = rr;
                    coincidenceLine(rt) = (k1_House.*rr.*AVGEnergy+k2_House.*sqrt(rr.*AVGEnergy))./(rr.*(k1_House*AVGEnergy+k2_House*sqrt(AVGEnergy))); % Velanders formula
                    P_demand(rt) = max(AVG_LoadProfile)*rr*coincidenceLine(rt);
                    I_line = coincidenceLine(rt).*fuse*rr;
                end

                %targetFuseRatings(find((ceil(I_line)./targetFuseRatings<1)~=0,1, 'first'));
                I_line;
                I_line_fuse(rt) = targetFuseRatings(find((ceil(I_line)./targetFuseRatings<=1)~=0,1, 'first'));
                index = find((I_line_fuse(rt)./targetFuseRatings<=1)~=0,1, 'first') ;
                Iu = fuses_tripping_size(2,(index));     
                ZmaxThick(rt) = 1000*c*Ufn/Iu;
                %Z_loop(:,rt) = targetLineImpedance.*(rt-1)*d_long(rt)+Z_transformer*1000
                
                L_max(:,rt) = (ZmaxThick(rt)-sum(Z_loop))./targetLineImpedance';
                
                while sum(L_max(:,rt)>d_long(rt)) == 0
                    rr = round(rr/2);
                    RX_multiplier(rt) = 0.5;
                    %mp(rt) = rr;
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
        
            rr = n;
            ZmaxThick;
            Z_loop;
            I_line_fuse;
            coincidenceLine;
            I_line;
            L_max;
            ixCable;
            %unity=(L_max>d_long);
            %L_max = unity.*L_max
            L_max(isnan(L_max))=0;
            z_earth(1) = +Z_transformer*1000;
            
            
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

            
            id = find(ixCable==max(ixCable),1);
            
            for hh=1:id
                Cable(hh) = CableCapacity(ixCable(id));
                R(hh) = Z_lineR_list(ixCable(id))*d_long(id);
                X(hh) = Z_lineX_list(ixCable(id))*d_long(id);
                z_earth(hh+1) = z_earth(hh) + targetLineImpedance(ixCable(hh))*d_long(hh)*RX_multiplier(hh);
            end
            
            t = 0;
Cable;
            


%             while max(z_earth) > Design_Z
%                 for tt = 1:p
%                     for hh = p:-1:1
%                         gg = ixCable(hh)
%                         
%                         R(hh) = Z_lineR_list(gg+1)*d_long(hh);
%                         X(hh) = Z_lineX_list(gg+1)*d_long(hh);
%                         z_earth(hh) = z_earth(hh-1) + targetLineImpedance(gg+1)*d_long(hh)*RX_multiplier(hh)
%                         
%                     end
                
%                 while status==0
%                     t = t +1;
%                     for hh=1:p
%                         gg=ixCable(hh)+t;
%                         if gg >= 8
%                            status = 1;
%                            break
%                         end
%                         Cable(hh) = CableCapacity(gg);
%                         R(hh) = Z_lineR_list(gg)*d_long(hh);
%                         X(hh) = Z_lineX_list(gg)*d_long(hh);
%                         %z_earth(hh) = targetLineImpedance(gg)*d_long(hh)+Z_transformer*1000;
%                         z_earth(hh+1) = z_earth(gg) + targetLineImpedance(ixCable(gg))*d_long(hh)*RX_multiplier(hh);
%                         if max(z_earth) <= Design_Z
%                             status = 1;
%                             break
%                         end
% 
%                     end
%                     
%                 end

%                 end
%                 pause
%             end
            
            
            
            
            %z_loop = repelem(z_loop,round(NumberOfCustomers/p));
            CustomersPerArea = NumberOfCustomers;
%NumberOfTransformers = repelem(NumberOfTransformers,length(z_earth));
            R = R./1000;
            X = X./1000;
            
Cable;
P_demand;
        N = sum(n)/sum(mc);
        Vn = 400;
        V = zeros(1,p+1);
        V(1) = 400 -coincidenceTR*CustomersPerTransformer*(Transformer_R.*(fuse*230*3)+Transformer_X.*(fuse*230*3)*PowerFactor)./(Vn);

        for kk=1:p
            V(kk+1) = V(kk)-(R(kk).*(P_demand(kk))+(X(kk)).*P_demand(kk)*PowerFactor)./(Vn);
        end
%         V;
%         P_demand;
%         coincidenceLine;
%         for pl=p:-1:1   
%             if min(V./Vn)>Design_voltage(2)
%                 break
%             else
%                 if find((Z_lineR_list./(R(pl)/d_long(pl))<=1)~=0,1, 'first')<length(Z_lineR_list)-1
%                     nRX = find((Z_lineR_list./(R(pl)/d_long(pl))<=1)~=0,1, 'first');
%                     %nX = find((Z_lineX_list./(X(pl)/d)<1)~=0,1, 'first')
%                     R(pl) = Z_lineR_list(nRX+1)*d_long(pl)*RX_multiplier(p1);
%                     X(pl) = Z_lineX_list(nRX+1)*d_long(pl)*RX_multiplier(p1);
%                 else
% 
%                 end
%             end
%         Vn = 230;
%         V = zeros(1,p+1);
%         V(1) = 230; %-(Rgrid.*(P_demand(1))+(Xgrid).*P_demand(1)*PowerFactor)./(Vn);
% 
%         for kk=1:p
%             V(kk+1) = V(kk)-((R(kk)).*(P_demand(kk))+(X(kk)).*P_demand(kk)*PowerFactor)./(Vn);
%         end
%         
%         end
    Cable;
    
          while  min(V./Vn)<Design_voltage(2)
              [val pos] = min(ixCable);
                if min(ixCable) == 7
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
                    V(kk+1) = V(kk)-(R(kk).*(P_demand(kk).*1000)+(X(kk)).*P_demand(kk).*1000*PowerFactor)./(Vn);
                end
                V;
          end
            
          
          
          
          % Checking that loop impedance is lower than maximum value
          z_earth = z_earth/1000;

          while  max(z_earth)>Design_Z
              
              [val pos] = min(ixCable);
                if min(ixCable) == 7
                    [val2 pos] = max(RX_multiplier);
                    RX_multiplier(pos) = RX_multiplier(pos)*0.5;
                    
                    for gg = pos:length(ixCable)
                        z_earth(gg+1) = z_earth(gg) + Z_Line(val+1)*d_long(gg)*RX_multiplier(gg)/1000;
                    end
                    R(pos) = Z_lineR_list(val+1)*d_long(pos)*RX_multiplier(pos)/1000;
                    X(pos) = Z_lineX_list(val+1)*d_long(pos)*RX_multiplier(pos)/1000;
                    %z_earth(pos+1) =  Z_Line(val)*d_long(pos)*RX_multiplier(pos);
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
    V;
    
    
N;
Z;

%A = tril(ones(length(Z)),0)';
%z_loop = mean(Z*A);
Z;
AVG_LoadProfile=AVG_LoadProfile*no_load;
coincidenceLine;
size(AVG_LoadProfile);
% mp
% size(mp)
% mc
% R
% X
% N
% PowerFactor
% Vn
% size(AVG_LoadProfile)
% size(Sol)
% Transformer_X
% Transformer_R
% coincidenceLine
% coincidenceTR
% CustomersPerTransformer
% voltageLimit

CustomersPerTransformer;
V0 = min(400 - max((n*(Transformer_R.*(coincidenceLine(1)*AVG_LoadProfile)+Transformer_X.*(coincidenceLine(1)*AVG_LoadProfile*PowerFactor))./Vn))); %-min(N(1).*(Rgrid.*(coincidenceLine(1).*AVG_LoadProfile-(ll).*Sol.*1000)+Xgrid.*(coincidenceLine(1).*AVG_LoadProfile*PowerFactor))./Vn);
Cable;

voltage = min(V0- sum(max(m.*(R.*(coincidenceLine.*AVG_LoadProfile)+X.*(coincidenceLine.*AVG_LoadProfile*PowerFactor))./Vn)));
thermal_limit;
TransformerType;
mp;
%ll=0;
    for ll=0.5:0.5:100
        % Just add voltage calcuation for pv
        ll;
        %N(1)*min(Rgrid.*(coincidenceLine(1).*AVG_LoadProfile-(ll).*Sol.*1000))/230;
        V0 = 400 - min((CustomersPerTransformer*(Transformer_R.*(coincidenceTR*fuse*230*3-(ll).*Sol.*1000)+Transformer_X.*(coincidenceTR*fuse*230*3))./Vn)); %-min(N(1).*(Rgrid.*(coincidenceLine(1).*AVG_LoadProfile-(ll).*Sol.*1000)+Xgrid.*(coincidenceLine(1).*AVG_LoadProfile*PowerFactor))./Vn);
    if length(mp) == 1
       voltage = V0- min((mp.*(R.*(coincidenceLine.*AVG_LoadProfile-(ll).*Sol.*1000)+X.*(coincidenceLine.*AVG_LoadProfile*PowerFactor))./Vn)');
    else
        voltage = V0- min(sum((mp.*(R.*(coincidenceLine.*AVG_LoadProfile-(ll).*Sol.*1000)+X.*(coincidenceLine.*AVG_LoadProfile*PowerFactor))./Vn)'));
    end
        max(AVG_LoadProfile-ll.*Sol.*1000);
        min(AVG_LoadProfile-ll.*Sol.*1000);
        deltaCurrentPerCustomerDouble = mc.*max(abs(coincidenceLine.*AVG_LoadProfile-ll.*Sol.*1000))/400/sqrt(3);

       
        %max(coincidenceTR*AVG_LoadProfile);
        %max(Sol);
        deltaPowerPerCustomer = max(abs(coincidenceTR*AVG_LoadProfile-ll.*Sol.*1000));
        deltaPowerPerCustomerSingle = max(abs(AVG_LoadProfile-ll.*Sol.*1000));
        %deltaPowerPerCustomerDouble = n*max(abs(coincidenceLine(1).*AVG_LoadProfile-ll.*Sol.*1000))
        deltaPower = deltaPowerPerCustomer.*CustomersPerTransformer;

        
        
        if (max(max(voltage))/400)>voltageLimit(1) %|| sum(voltage<voltageLimit(2))
            solcap = (ll-0.5)*NumberOfCustomers;
            VoltLimit = (max(max(voltage))/400);
            Limiter = 10;
            break
            
        elseif deltaCurrentPerCustomerDouble>Cable %|| deltaPowerPerCustomerDouble*nk/(230)>CableThick
                VoltLimit = (max(max(voltage))/400);
                solcap = (ll-0.5)*NumberOfCustomers;
                Limiter = 25;
                break
        elseif deltaPower>(TransformerType*1000)*thermal_limit
            VoltLimit = (max(max(voltage))/400);
            solcap = (ll-0.5)*NumberOfCustomers;
            Limiter = 20;
            break 
        end
    end
    

Limiter ;
deltaPower;
TransformerType;
    EnergyPerkW = sum(Sol)/6;                               % i kWh
    EnergyPerKM2 = (NumberOfCustomers)*(ll-0.5)*sum(Sol)/6;       % i kW
    CapPerCustomer = (ll-0.5);
    CustomersPerKm = CustomersPerTransformer/LengthLVPerTransformer;
    CustomersCalc = CustomersPerTransformer*NumberOfTransformers;
    CustomersInitial = Pop_density/PeoplePerHousehold;
    %LengthLVPerKm2 = LengthLVPerTransformer*NubmerOfTransformers;
    
    solargeneration = (ll-0.5).*Sol.*NumberOfCustomers;
  
end










