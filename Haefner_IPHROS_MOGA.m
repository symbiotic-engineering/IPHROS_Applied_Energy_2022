% Matthew Haefner
% Applied Energy
% IPHROS MOGA

% Note: all V and E values are m^3/day and kWh/day


%% MOGA

% x = [E_r,gamma,gamma_RO,h,N_e1,N_e2,N_pv1,N_pv2]
fun = @objective;
nvars = 8;
A = [0 0 0 0 0 0 -1 1];
b = 0;
Aeq = [];
beq = [];
lb = [1 0.01 0.01 240 1 0 1 0];
ub = [100000e3 0.99 0.99 821 8 8 1.86e6 1.86e6];
nonlcon = @constraints;
intcon = [5,6,7,8];

options = optimoptions('gamultiobj','Display','iter','MaxGenerations',500,'UseParallel', true, 'UseVectorized', false);
rng default % For reproducibility
[x,fval,exitflag,output,population,scores] = gamultiobj(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,intcon,options);

%% Post Processing (minor)

E_er_array = -fval(:,1); % kWh/day
V_fwRO_array = -fval(:,2); % m^3/day
C_array = fval(:,3); % $ (proxy)

results = [E_er_array,V_fwRO_array,C_array];

%% Functions
% Objective Function
function J = objective(x)
    % General Parameters
    g = 9.81;             % gravitational acceleration [m/s^2]
    rho_sw = 1023.6;      % generic density of saltwater at 1 atm and 25 degrees celsius [kg/m^3]
    rho_fw = 996.9;       % generic density of freshwater at 1 atm and 25 degrees celsius [kg/m^3]
    eta_hp = 0.894;       % pump-side efficiency
    eta_ht = 0.894;       % pump-side efficiency [-]
    T = 25.0;             % temperature [degrees celsius]
    S_sw = 35.0;          % generic seawater salinity [g salt/kg seawater]
    M_salt = 58.44;       % molar mass of NaCl [g/mol]
    P_p = 14.696;         % permeate pressure, atmospheric for now [psi]
    R = 8.3145;           % universal gas constant [J/(mol*K)]
    MW_salt_avg = 30.548; % average molar mass of NaCl [g/mol]
    Am = 440;             % membrane surface area [ft^2]
    M_cost_e1 = 1.7;      % Proxy for membrane cost of 1st element in RO unit [$/unit]
    M_cost_e2 = 1;        % Proxy for membrane cost of other elements in RO unit [$/unit] 
    
    % Unpack x
    E_r = x(1); gamma = x(2); gamma_RO = x(3); h = x(4);
    N_e1 = x(5); N_e2 = x(6); N_pv1 = x(7); N_pv2 = x(8);
    
    % Go through the relevant modules
    [E_rp,E_rd] = initial_energy_division(E_r,gamma);
    V_wp = pump(E_rp,h,eta_hp,rho_sw,g);
    [V_wRO, V_swht, Q_f] = mountaintop_reservoir(V_wp,gamma_RO,N_pv1);
    [V_fwRO,V_oRO,S_oRO,rho_oRO,eta_ROio,eta_RO,rr_max,rr_min,Q_f_max,Q_f_min,Q_c_max,Q_c_min,Q_p_max,Q_p_min,LHS_SVM] = ro_system(h,Q_f,N_e1,N_e2,N_pv1,N_pv2,rho_sw,rho_fw,T,S_sw,M_salt,P_p,R,MW_salt_avg,Am);

    E_swht = turbine(V_swht,h,rho_sw,g,eta_ht);
    E_htRO = energy_recapture(V_oRO,h,eta_ROio,rho_oRO,g,eta_ht);
    E_er = energy_summation(E_rd,E_swht,E_htRO);
    C = RO_cost(N_e1,N_e2,N_pv1,N_pv2,M_cost_e1,M_cost_e2);

    E_er_min = -E_er; % maximizing
    V_fwRO_min = -V_fwRO; % maximizing
    C_min = C; % minimizing
    
    % Compiling
    J = [E_er_min,V_fwRO_min,C_min];
end

% Constraint Function
function [c,ceq] = constraints(x)
    % General Parameters
    g = 9.81;             % gravitational acceleration [m/s^2]
    rho_sw = 1023.6;      % generic density of saltwater at 1 atm and 25 degrees celsius [kg/m^3]
    rho_fw = 996.9;       % generic density of freshwater at 1 atm and 25 degrees celsius [kg/m^3]
    eta_hp = 0.894;       % pump-side efficiency
    eta_ht = 0.894;       % pump-side efficiency [-]
    T = 25.0;             % temperature [degrees celsius]
    S_sw = 35.0;          % generic seawater salinity [g salt/kg seawater]
    M_salt = 58.44;       % molar mass of NaCl [g/mol]
    P_p = 14.696;         % permeate pressure, atmospheric for now [psi]
    R = 8.3145;           % universal gas constant [J/(mol*K)]
    MW_salt_avg = 30.548; % average molar mass of NaCl [g/mol]
    Am = 440;             % membrane surface area [ft^2]
    M_cost_e1 = 1.7;      % Proxy for membrane cost of 1st element in RO unit [$/unit]
    M_cost_e2 = 1;        % Proxy for membrane cost of other elements in RO unit [$/unit] 
    
    % Unpack x
    E_r = x(1); gamma = x(2); gamma_RO = x(3); h = x(4);
    N_e1 = x(5); N_e2 = x(6); N_pv1 = x(7); N_pv2 = x(8);
    
    % Go through the relevant modules
    [E_rp,E_rd] = initial_energy_division(E_r,gamma);
    V_wp = pump(E_rp,h,eta_hp,rho_sw,g);
    [V_wRO, V_swht, Q_f] = mountaintop_reservoir(V_wp,gamma_RO,N_pv1);
    [V_fwRO,V_oRO,S_oRO,rho_oRO,eta_ROio,eta_RO,rr_max,rr_min,Q_f_max,Q_f_min,Q_c_max,Q_c_min,Q_p_max,Q_p_min,LHS_SVM] = ro_system(h,Q_f,N_e1,N_e2,N_pv1,N_pv2,rho_sw,rho_fw,T,S_sw,M_salt,P_p,R,MW_salt_avg,Am);

    [rho_ht,S_ht,V_ht] = mixing(V_oRO,V_swht,rho_oRO,S_oRO,gamma_RO,eta_RO,rho_sw,S_sw);

    % Salinity constraint
    S_ht_con = S_ht-40;

    % Feed flow rate to pressure vessel constraint (gpd constraint)
    Q_f_min_con = 21620-Q_f_min;
    Q_f_max_con = Q_f_max-98272;

    % Concentrate flow rate constraint (gpd constraint)
    Q_c_min_con = 21620-Q_c_min;
    Q_c_max_con = Q_c_max-98272;

    % Permeate flow rate constraint (gpd constraint)
    Q_p_min_con = -Q_p_min;
    Q_p_max_con = Q_p_max-8369;

    % Recovery rate for element constraints
    rr_min_con = -rr_min;
    rr_max_con = rr_max-0.13;

    % If N_e2 = 0, then N_pv2 = 0
    e2_pv2_con = abs(min([N_e2,1]) - min([N_pv2,1]));

    % Equation for SVM plane
    SVM_con = LHS_SVM;

    % Compiling
    c = [S_ht_con;Q_f_min_con;Q_f_max_con;Q_c_min_con;Q_c_max_con;Q_p_min_con;Q_p_max_con;rr_min_con;rr_max_con;e2_pv2_con;SVM_con];
    ceq = [];
end

% Modules
function [E_rp,E_rd] = initial_energy_division(E_r,gamma)
    E_rp = gamma*E_r; % kWh/day
    E_rd = (1-gamma)*E_r; % kWh/day
end

function V_wp = pump(E_rp,h,eta_hp,rho_sw,g)
    E_rp = E_rp*(3.6*10^6); % joules/day (from kWh/day)
    V_wp = E_rp*eta_hp/(rho_sw*g*h); % m^3/day
end

function [V_wRO,V_swht,Q_f] = mountaintop_reservoir(V_wp,gamma_RO,N_pv1)
    V_wRO = gamma_RO*V_wp; % m^3/day
    V_swht = (1-gamma_RO)*V_wp; % m^3/day
    Q_f = (V_wRO*264.172)/N_pv1; % gpd (from m^3/day)
end

function [V_fwRO,V_oRO,S_oRO,rho_oRO,eta_ROio,eta_RO,rr_max,rr_min,Q_f_max,Q_f_min,Q_c_max,Q_c_min,Q_p_max,Q_p_min,LHS_SVM] = ro_system(h,Q_f_input,N_e1,N_e2,N_pv1,N_pv2,rho_sw,rho_fw,T,S_sw,M_salt,P_p,R,MW_salt_avg,Am)

    % Initialize oopsies relative to Q_f_input
    if Q_f_input <= 98272 && Q_f_input >= 21620 % gpd
        oopsies = 0;
        oopsies1 = 0;
        oopsies2 = 0;
    else
        oopsies = 1;
        oopsies1 = 1;
        oopsies2 = 1;
    end

    % Collection Lists
    Q_p_collection = []; % gpd
    P_c_collection = []; % psi
    c_c_collection = []; % mol/L
    Q_c_collection = []; % gpd
    rr_collection = [];  % -
    P_f_collection = []; % psi
    S_f_collection = []; % g/kg
    Q_f_collection = []; % gpd
    
    % Initialize 1st stage element counter
    i1 = 1;
    
    % First stage
    while i1 <= N_e1 && oopsies ~= 1
        % Beginning with feed conditions into an element
        if i1 == 1
            % First element
            P_f = 1.4228424*h*(rho_sw/rho_fw); % psi
            c_f = (S_sw*rho_sw)/(1000*M_salt); % mol/L
            Q_f = Q_f_input; % gpd
        else
            % Feed conditions from previous element
            P_f = P_c_collection(i1-1); % psi
            Q_f = Q_c_collection(i1-1); % gpd
            c_f = c_c_collection(i1-1); % mol/L
        end

        pi_f = 1.12*(273+T)*2*c_f; % psi
        S_f = salinity_from_osmotic(pi_f,rho_fw,T,R,MW_salt_avg); % g/kg

        P_f_collection = [P_f_collection;P_f];
        S_f_collection = [S_f_collection;S_f];
        Q_f_collection = [Q_f_collection;Q_f];

        % Solving for Q_p and other numeric values %%%%%%%%%%%%%%%%%%%%%%%
        [Q_p_num,P_c_num,c_c_num,Q_c_num,rr_num,oopsies1] = symbolic_Q_p_stuff(c_f,T,Q_f,Am,P_f,P_p);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if oopsies1 == 1
            i1 = N_e1 + 1;
            oopsies2 = 1;
        else
            % Storing element outputs for further calculations
            Q_p_collection = [Q_p_collection;Q_p_num]; % gpd
            P_c_collection = [P_c_collection;P_c_num]; % psi
            c_c_collection = [c_c_collection;c_c_num]; % mol/L
            Q_c_collection = [Q_c_collection;Q_c_num]; % gpd
            rr_collection = [rr_collection;rr_num];    % -
            
            % Increasing counter
            i1 = i1 + 1;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up storage variables
    Q_p_collection2 = []; % gpd
    P_c_collection2 = []; % psi
    c_c_collection2 = []; % mol/L
    Q_c_collection2 = []; % gpd
    rr_collection2 = [];  % - 
    P_f_collection2 = []; % psi
    S_f_collection2 = []; % g/kg
    Q_f_collection2 = []; % gpd

    % Second stage
    if N_e2 == 0 || N_pv2 == 0 || oopsies == 1 || oopsies1 == 1
        i2 = N_e2+1;
    else        
        % initialize second stage element counter
        i2 = 1;
    end

    while i2 <= N_e2
        if i2 == 1
            % Feed conditions from previous stage (feed flowrate reconfigured)
            P_f = P_c_collection(end); % psi
            c_f = c_c_collection(end); % mol/L
            
            V_intermediate = Q_c_collection(end)*N_pv1; % gpd
            Q_f = V_intermediate/N_pv2; % gpd
        else
            % Feed conditions from previous element
            P_f = P_c_collection2(end); % psi
            Q_f = Q_c_collection2(end); % gpd
            c_f = c_c_collection2(end); % mol/L
        end

        pi_f = 1.12*(273+T)*2*c_f; % psi
        S_f = salinity_from_osmotic(pi_f,rho_fw,T,R,MW_salt_avg); % g/kg

        P_f_collection = [P_f_collection;P_f];
        S_f_collection = [S_f_collection;S_f];
        Q_f_collection = [Q_f_collection;Q_f];
        
        % Solving for Q_p and other numeric values %%%%%%%%%%%%%%%%%%%%%%%
        [Q_p_num,P_c_num,c_c_num,Q_c_num,rr_num,oopsies2] = symbolic_Q_p_stuff(c_f,T,Q_f,Am,P_f,P_p);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if oopsies2 == 1
            i2 = N_e2 + 1;
        else
            % Storing element outputs for further calculations
            Q_p_collection2 = [Q_p_collection2;Q_p_num]; % gpd
            P_c_collection2 = [P_c_collection2;P_c_num]; % psi
            c_c_collection2 = [c_c_collection2;c_c_num]; % mol/L
            Q_c_collection2 = [Q_c_collection2;Q_c_num]; % gpd
            rr_collection2 = [rr_collection2;rr_num];    % -
            
            % Increasing counter
            i2 = i2 + 1;
        end
    end
    if oopsies == 1 || oopsies1 == 1 || oopsies2 == 1
        % Dependent vars
        V_fwRO = 0; V_oRO = 0; S_oRO = 1e5; rho_oRO = 1e5; eta_ROio = 0; eta_RO = 0;

        % Constraint vars
        rr_max = 1e6; rr_min = -1e6; Q_f_max = 1e6; Q_f_min = -1e6; 
        Q_c_max = 1e6; Q_c_min = -1e6; Q_p_max = 1e6; Q_p_min = -1e6;
        LHS_SVM = 1e6;
        
    else
        % Processing for RO System Output Values
        V_fwRO = N_pv1*(sum(Q_p_collection)/264.172) + N_pv2*(sum(Q_p_collection2)/264.172); % m^3/day
        
        brine_collective = [Q_c_collection;Q_c_collection2];
        pressure_vessel_collective = [N_pv1;N_pv2];
        V_oRO = (brine_collective(end)/264.172)*pressure_vessel_collective(end);  % m^3/day (from gpd)
        
        c_c_collective = [c_c_collection;c_c_collection2];
        pi_c = 1.12*(273+T)*2*c_c_collective(end); % psi
        S_oRO = salinity_from_osmotic(pi_c,rho_fw,T,R,MW_salt_avg); % g/kg
    
        P_c_collective = [P_c_collection;P_c_collection2];
        rho_oRO = density_from_salinity(S_oRO,P_c_collective(end),rho_fw,T,P_p); % kg/m^3
        
        P_f_1st_element = 1.4228424*h*(rho_sw/rho_fw);                        % psi
        eta_ROio = 1-((P_f_1st_element-P_c_collective(end))/P_f_1st_element); % -
    
        rr_collective = [rr_collection;rr_collection2];
        eta_RO_product_term = 1;
        i_eta_RO = 1;
        N_e = N_e1 + N_e2;
        while i_eta_RO <= N_e
            eta_RO_product_term = eta_RO_product_term*(1-rr_collective(i_eta_RO));
            i_eta_RO = i_eta_RO + 1;
        end
        eta_RO = 1-eta_RO_product_term; % -
    
        %%%%%%%%%%%%%% (hyperplane constraint variable)
        P_f_array = [P_f_collection;P_f_collection2];
        S_f_array = [S_f_collection;S_f_collection2];
        Q_f_array = [Q_f_collection;Q_f_collection2];
        LHS_SVM = -1e6;
        for i = 1:length(P_f_array)
            var = -1334.934 + (18.222*P_f_array(i)) - (282.748*S_f_array(i)) - (200*(Q_f_array(i)/(24*264.172)));
            if var > LHS_SVM
                LHS_SVM = var;
            end
        end
        %%%%%%%%%%%%%%
    
        % Constraint Variables
        rr_max = max(rr_collective); 
        rr_min = min(rr_collective); 
        Q_f_max = Q_f_input; 
        Q_f_min = Q_f; 
        Q_c_max = max([Q_c_collection;Q_c_collection2]); 
        Q_c_min = min([Q_c_collection;Q_c_collection2]); 
        Q_p_max = max([Q_p_collection;Q_p_collection2]); 
        Q_p_min = max([Q_p_collection;Q_p_collection2]);
    end
end

function E_swht = turbine(V_swht,h,rho_sw,g,eta_ht)    
    E_swht_joules = rho_sw*V_swht*g*h*eta_ht; % joules/day
    E_swht = E_swht_joules/(3.6e6); % kWh/day
end

function E_htRO = energy_recapture(V_oRO,h,eta_ROio,rho_oRO,g,eta_ht)    
    E_htRO_joules = rho_oRO*V_oRO*g*h*eta_ht*eta_ROio; % joules/day
    E_htRO = E_htRO_joules/(3.6e6); % kWh/day
end

function [rho_ht,S_ht,V_ht] = mixing(V_oRO,V_swht,rho_oRO,S_oRO,gamma_RO,eta_RO,rho_sw,S_sw)
    rho_ht = ((rho_sw*(1-gamma_RO)) + (gamma_RO*rho_oRO*(1-eta_RO)))/(1-(eta_RO*gamma_RO)); % kg/m^3
    S_ht = ((S_sw*rho_sw*(1-gamma_RO)) + (gamma_RO*S_oRO*rho_oRO*(1-eta_RO)))/((rho_sw*(1-gamma_RO)) + (gamma_RO*rho_oRO*(1-eta_RO))); % g/kg
    V_ht = V_swht + V_oRO; % m^3/day
end

function E_er = energy_summation(E_rd,E_swht,E_htRO)
    E_er = E_rd+E_swht+E_htRO; % kWh/day
end

function C = RO_cost(N_e1,N_e2,N_pv1,N_pv2,M_cost_e1,M_cost_e2)
    C1 = N_pv1*(M_cost_e1 + ((N_e1-1)*M_cost_e2)); % $ (proxy)
    C2 = N_pv2*(M_cost_e1 + ((N_e2-1)*M_cost_e2)); % $ (proxy)
        
    C = C1 + C2; % $ (proxy)
end

% Miscellaneous Functions
function S_brine = salinity_from_osmotic(pi,rho_fw,T,R,MW_salt_avg)
    % Get Temperature and freshwater density
    t = T; % celsius
    T = (T)+273; % kelvin (from celsius)
    
    % Osmotic Coefficient Correlation (good for S>10)
    a1 = 8.9453233003*10^-1;
    a2 = 4.1560737424*10^-4;
    a3 = -4.6262121398*10^-6;
    a4 = 2.2211195897*10^-11;
    a5 = -1.1445456438*10^-4;
    a6 = -1.4783462366*10^-6;
    a7 = -1.3526263499*10^-11;
    a8 = 7.0132355546*10^-6;
    a9 = 5.6960486681*10^-8;
    a10 = -2.8624032584*10^-10;
    
    syms S real % g/kg
    
    phi = a1+(a2*t)+(a3*(t^2))+(a4*(t^4))+(a5*S)+(a6*S*t)+(a7*S*(t^3))+...
        (a8*(S^2))+(a9*t*(S^2))+(a10*(S^2)*(t^2));
    
    % Converting osmotic pressure units
    pi = pi/145.038; % MPa (from psi)
    
    % Getting numerical value for S
    eqn = pi == ((phi*R*T*rho_fw)/(10^6))*(1000/MW_salt_avg)*...
        (S/(1000-S));
    S_num = double(solve(eqn,S));
    
    S_brine = S_num; % g/kg
end

function rho_brine = density_from_salinity(S,P,rho_fw,T,P_p)
    t = T; % celsius
    P_o = P_p/145.038; % MPa (from psi)
    
    % Density at Atmospheric
    b1 = 8.020*10^2;
    b2 = -2.001;
    b3 = 1.677*10^-2;
    b4 = -3.060*10^-5;
    b5 = -1.613*10^-5;
    
    S_kg_kg = S/1000; % kg/kg
    
    rho_sw_atm = rho_fw+((b1*S_kg_kg)+(b2*S_kg_kg*t)+(b3*S_kg_kg*(t^2))+...
        (b4*S_kg_kg*(t^3))+(b5*(S_kg_kg^2)*(t^2))); % kg/m^3
    
    % Other Correction Term
    c1 = 5.0792*10^-4;
    c2 = -3.4168*10^-6;
    c3 = 5.6931*10^-8;
    c4 = -3.7263*10^-10;
    c5 = 1.4465*10^-12;
    c6 = -1.7058*10^-15;
    c7 = -1.3389*10^-6;
    c8 = 4.8603*10^-9;
    c9 = -6.8039*10^-13;
    
    d1 = -1.1077*10^-6;
    d2 = 5.5584*10^-9;
    d3 = -4.2539*10^-11;
    d4 = 8.3702*10^-9;
    
    P = P/145.038; % MPa (from psi)
    
    Fp = exp(((P-P_o)*(c1+(c2*t)+(c3*(t^2))+(c4*(t^3))+(c5*(t^4))+...
        (c6*(t^5))+(S*(d1+(d2*t)+(d3*(t^2))))))+((0.5*((P^2)-(P_o^2)))*...
        (c7+(c8*t)+(c9*(t^3))+(d4*S))));
    
    rho_brine = rho_sw_atm*Fp; % kg/m^3
end

% RO helper function (handles symbolic Q_p stuff)
function [Q_p_num,P_c_num,c_c_num,Q_c_num,rr_num,oopsies] = symbolic_Q_p_stuff(c_f,T,Q_f,Am,P_f,P_p)
    % Set up symbolic Q_p
    syms Q_p real
        
    % Getting Permeate Osmotic Pressure
    R = R_value(Q_p);            % -
    c_p = c_f*(1-R);             % mol/L
    pi_p = 1.12*(273+T)*(2*c_p); % psi
    
    % Getting Osmotic Pressure along Membrane
    rr = Q_p/Q_f;                    % -
    pf = exp(0.7*rr);                % -
    Q_c = Q_f-Q_p;                   % gpd
    c_c = ((Q_f*c_f)-(Q_p*c_p))/Q_c; % mol/L
    c_fc = (c_f+c_c)/2;              % mol/L
    c_m = (pf*(c_fc-c_p))+c_p;       % mol/L
    pi_m = 1.12*(273+T)*(2*c_m);     % psi
    
    % Getting Pressure Drop along Element
    Q_fc = 0.5*(Q_f+Q_c);           % gpd
    dP_fc = 0.01*((Q_fc/1440)^1.7); % psi
    
    % Getting Aw
    Aw = Aw_value(Q_p,c_p); % gfd/psi

    % Getting Numerical Value for Q_p
    eqn = Q_p == Aw*Am*(P_f-(0.5*dP_fc)-P_p-pi_m+pi_p);
    Q_p_num = double(vpasolve(eqn,Q_p,[60 8370])); % gpd

    if isempty(Q_p_num) == 1
        oopsies = 1;
        Q_p_num = 0;
        P_c_num = 0;
        c_c_num = 0;
        Q_c_num = 0;
        rr_num = 0;
    else
        oopsies = 0;
        % Subbing in Q_p to element outputs
        P_c_num = P_f-double(subs(dP_fc,Q_p,Q_p_num)); % psi
        c_c_num = double(subs(c_c,Q_p,Q_p_num));       % mol/L
        Q_c_num = double(subs(Q_c,Q_p,Q_p_num));       % gpd
        rr_num = double(subs(rr,Q_p,Q_p_num));         % -
    end
end

% Function for calculating R (symbolic)
function R = R_value(Q_p)
    R = 1.0034 - (12.2124*(Q_p^(-0.8122))); % -
end

% Function for calculating Aw (symbolic)
function Aw = Aw_value(Q_p,c_p)
    term1 = -85.71804787;
    term2 = 3.3699665311*log(c_p);
    term3 = 53.163612047*log(Q_p);
    term4 = 3.2938449078*log(c_p)*log(c_p);
    term5 = 0.726912142*log(c_p)*log(c_p)*log(c_p);
    term6 = 0.0516395755*log(c_p)*log(c_p)*log(c_p)*log(c_p);
    term7 = -12.22459358*log(Q_p)*log(Q_p);
    term8 = 1.1774872012*log(Q_p)*log(Q_p)*log(Q_p);
    term9 = -0.041266007*log(Q_p)*log(Q_p)*log(Q_p)*log(Q_p);
    
    term = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9;
    Aw = exp(term); % gfd/psi
end
