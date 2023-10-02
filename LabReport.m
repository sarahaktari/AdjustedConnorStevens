%% 
% *SARAH AKTARI* 
% 
% *The effect of injecting external currents on action potention firing rate 
% in an adjusted Connor Stevens neuron model* 
% 
% *ABSTRACT*
% 
% The objective of this study was to determine how applying varying levels of 
% current affects the firing rate of action potentials in a Connor Stevens model 
% of a neuron. In this study, we recorded the number of action potentials fired 
% with no applied current and compared that control value to the number of action 
% potentials fired after increasing the applied current by intervals of 100 picoAmps 
% until a max applied current of 3000 picoAmps. 30 data points of the firing rate 
% (action potentials/second) of the neuron model were collected, each with a different 
% applied current. The firing rate was plotted against the range for applied current, 
% and we found that injecting larger applied current results in an increased firing 
% rate and there is a jump in firing rate when the injected current is 1100pA. 
% In conclusion, injection of applied current in the Connor Stevens neuron model 
% leads to an increased firing of action potentials. 
% 
% *HYPOTHESIS*
% 
% We hypothesized that applying external currents to the neuron model will lead 
% to increased firing rates in the neuron. This would occur because applying external 
% depolarizing currents leads to a change in the conductance of Na+ channels which 
% further depolarizes the neuronal membrane due to Na+ influx. In these cases, 
% the increase in membrane potential causes the neuron to reach the threshold 
% voltage of -55mV, leading to the opening of both the m and h gates of voltage-gated 
% Na+ channels which results in the opening of the n gates of K+ channels and 
% a and b gates of the A-type K+ channels and hyperpolarizing currents. As a result, 
% we hypothesized that larger depolarizing currents would lead to increased firing 
% rates. 

% parameters
C_m = 100; %pF membrane capcitance 

E_leak = -70; %mV leak reversal potential that was modified to better model neuron at rest 
E_Na = 55; %mV sodium reversal potential
E_A = -75; %mV A-type potassium channel reversal potential 
E_K = -85; %mV potassium reversal potential 

G_leak = 30; %nS, leak conductance
g_Na_max = 12e3; %nS, maximum Na+ conductance
g_K_max = 2e3; %nS, maximun K+ conductance 
g_A_max = 4.77e3; %nS, maximum A-type conductance 

%% 
% *CONTROL AND EXPERIMENTAL CONDITIONS*
% 
% For our experimental conditions, we increased the degree of injected current 
% for each trial. We started by injecting 100 picoAmps, and recording how the 
% firing rate at this level of applied current compared to the firing rate of 
% the neuron with no applied current. We applied positive currents as this would 
% cause the membrane potential to depolarize, which will lead to a Na+ influx. 
% Once the neuron reaches the threshold, it would fire an action potential. We 
% chose to increase the applied current by intervals of 100pA until the maximum 
% applied current of 3000pA as it would provide sufficient enough data to create 
% an f-I curve displaying the trend of firing rates against applied current.
% 
% For our control condition, we recorded the behavior of the neuron model with 
% no external applied current. With no applied current, the neuron remained at 
% rest and fired no action potentials, and stayed at its resting leak potential 
% of -70mV. This means at rest, the firing rate of the neuron was zero and this 
% served as a baseline for comparing firing rates as a result of applied currents. 
% 
% *FIGURES*
% 
% _*Figure 1. Connor Stevens model neuron (with applied current of 3000pA at 
% 100ms)*_
% 
% This is a figure of voltage over time of the Connor Stevens model when there 
% is a max applied current of 3000pA. The peaks represent an action potential, 
% in which the membrane potential rises until it reaches the threshold at roughly 
% -55mV and then dips back and hyperpolarizes and the cycle continues as long 
% as there is an applied current to trigger the action potentials. The time frame 
% within this figure is 100 milliseconds. 
% 
% _*Figure 2. f-I curve of firing rate*_ 
% 
% This is a figure of the firing rates observed in the Connor Stevens model 
% at each level of applied current. When there is no applied current, the firing 
% rate is zero. The firing rate steadily in a linear manner with applied current 
% until a current of 1100pA is injected, where the firing rate drastically jumps 
% and then increases steadily again. 
% 
% *CONCLUSIONS*
% 
% From our experimental results, we can conclude that there exists a positive 
% correlation between injected current and the firing rate of the Connor Stevens 
% neuron model. From the results recorded in figure 2, we can conclude that this 
% neuron is a type II neuron, similar to a Hodgkin and Huxley model as there is 
% a drastic jump in firing rate when 1100pA of current is injected. The firing 
% rate of the neuron appears to increase steadily at all other points besides 
% at 1100 pA as the applied current increases, meaning there is a positive relationship 
% between the firing rate and applied current. This further supports our initial 
% hypothesis in which we expected the neuron to fire more action potentials as 
% more depolarizing current was injected into it. We had not expected this model 
% to behave as a type II neuron since the addition of the A-type potassium channel 
% was meant to slow down the rate of firing action potentials since it is an additional 
% hyperpolarizing current that would result in a longer refractory period than 
% seen in the Hodgkin and Huxley model. This would lead to smaller numbers of 
% action potentials fired over time as we would expect the neuron to need more 
% time to approach its resting membrane potential after hyperpolarization. However, 
% the parameters of the Connor Stevens model were modified in this experiment 
% to better mimic a neuron at rest, which can be the reason why is model behaves 
% as a type II neuron as opposed to a type I neuron. 
% 
% 
% 
% _Initialization of time, voltage, and gating variable vectors_

% time vector
dt = 0.01; %ms
tmax = 2000; %ms
time = (0:dt:tmax); %ms

% initialize vectors for variables
Vm = zeros(size(time)); %voltage vector
m = zeros(size(time)); % m gate vector variable in Na+ channel
n = zeros(size(time)); % h gate vector variable in Na+ channel
h = zeros(size(time)); % n gate vector variable in K+ channel
a = zeros(size(time)); % n gate vector variable in A-type potassium channel
b = zeros(size(time));% n gate vector variable in A-type potassium channel

Vm(1) = E_leak; %initialize voltage so at first index resting voltage is equal to leak conductance 
%% 
% _Currents for each channel with differential equations_ 

% applied current vector 
I_app = zeros(size(time)); 

I_step = (0:100:3000); %vector of applied current intervals 

count_AP = zeros(size(I_step)); % vector to store number of action potentials fired for each applied current interval 

for i = 1 : length(I_step) % outer for loop going through each interval of applied current 
    I_app(:) = I_step(i);

% Currents for each channel with differential equations 
  I_leak = zeros(size(time)); %current vector for leak channels 
  I_Na = zeros(size(time)); %current vector for Na+ channels
  I_K = zeros(size(time)); %current vector for K+ channels 
  I_A = zeros(size(time)); %current vector for A-type channels

for t = 1:length(time)-1 % inner for loop to find values of voltage over time 

    %rate constants 

    % m gate variable : equations from textbook
    alpha_m = (0.38.*(Vm(t)+29.7))./(1-exp(-0.1.*(Vm(t)+29.7)));
    beta_m = 15.2.*exp(-0.0556.*(Vm(t)+54.7));

    tau_m = 1./(alpha_m+beta_m); %time constant for reference (used in last section)
    m_inf = alpha_m./(alpha_m+beta_m); %steady state equation for reference (used in last section)

    % h variable 
    alpha_h = 0.266.*exp(-0.05.*(Vm(t)+48));
    beta_h = 3.8./(1+exp(-0.1.*(Vm(t)+18)));

    tau_h = 1./(alpha_h+beta_h); %time constant 
    h_inf = alpha_h./(alpha_h+beta_h); %steady state equation 

    % n variable 
    alpha_n = (0.02.*(Vm(t)+45.7))./(1-exp(-0.1.*(Vm(t)+45.7)));
    beta_n = 0.25.*exp(-0.0125.*(Vm(t)+55.7));

    tau_n = 1./(alpha_n+beta_n); %time constant for reference (used in last section)
    n_inf = alpha_n./(alpha_n+beta_n); %steady state equation 

    % a variable 
    a_inf = ((0.0761.*exp(0.0314.*(Vm(t))))./(1+exp(0.0346.*(Vm(t)+1.17))))^(1/3); %steady state equation 
    tau_a = 0.3632 + (1.158./(1+exp(0.0497.*(Vm(t)+55.96)))); %time constant for reference (used in last section)

    % b variable 
    b_inf = (1./(1+exp(0.0497.*(Vm(t)+53.3))))^4; %steady state equation
    tau_b = 1.24 + (2.678./(1+exp(0.0624.*(Vm(t)+50)))); %time constant for reference (used in last section)

 % initial condition for the gating variables for first index
 %  set first index of each gating variable to steady state value
    if t==1
        m(1) = m_inf;
        n(1) = n_inf;
        h(1) = h_inf;
        a(1) = a_inf;
        b(1) = b_inf;
    end

    % caluclate differential equations for each variable  
    
    % dm/dt
    dm = (alpha_m*(1-m(t)) - beta_m*(m(t)))*dt;
    m(t+1) = m(t) + dm;

    % dn/dt
    dn = (alpha_n*(1-n(t)) - beta_n*(n(t)))*dt;
    n(t+1) = n(t) + dn;

    % dh/dt
    dh = (alpha_h*(1-h(t)) - beta_h*(h(t)))*dt;
    h(t+1) = h(t) + dh;

    % da/dt
    da = ((a_inf-a(t))/tau_a)*dt; 
    a(t+1) = a(t) + da;

    % db/dt
    db = ((b_inf-b(t))/tau_b)*dt; 
    b(t+1) = b(t) + db;

    % each chunk of differential equation 
    I_leak(t) = -G_leak*(E_leak-Vm(t)); %l eak current
    I_Na(t) = -g_Na_max*m(t)*m(t)*m(t)*h(t)*(E_Na - Vm(t)); % Na current
    I_K(t) = -g_K_max*n(t)*n(t)*n(t)*n(t)*(E_K - Vm(t)); % K current
    I_A(t) = -g_A_max*a(t)*a(t)*a(t)*b(t)*(E_A - Vm(t)); % A-type current

    % the complete differential eqn for the membrane voltage
    dV = ((-I_leak(t) - I_Na(t) - I_K(t) - I_A(t) + I_app(t))/C_m)*dt; % negate to account for sign 
    Vm(t+1) = Vm(t) + dV;
end

    count_AP(i) = numel(findpeaks(Vm)); %counts peaks of action potentials after each applied current 

end 

rate_AP = count_AP ./ (tmax * 1e-3); %records firing rate in spikes per second

%% 
% _Voltage against time with applied current_

% plot of voltage against time in the first 100 milliseconds since there
% are so many spikes recorded 

figure;

plot(time,Vm,'LineWidth',3); grid on

xlabel('time (ms)')
ylabel('voltage (mV)')
title('figure 1. connor stevens model neuron')

xlim([0 100])

%% 
% _Calculating firing rate for each applied current_ 

% equation used to count number of action potentials for each I_app
%count_AP = numel(findpeaks(Vm)

figure % plotting f-I curve: this is a type I neuron
plot (I_step, rate_AP)
xlabel ('Applied Current (pA)')
ylabel ('Action Potential Firing Rate (# spikes/s)')
title ('figure 3. f-I curve')

%% 
%