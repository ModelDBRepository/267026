function ipRGCephys()

    Iapp = 100*1e-3; % applied current
    V0 = -30; % initial membrane voltage
    tf = 500; % final time
    
    figure(); hold all;
    
    for m1m4 = [1,4] % m1m4 - ipRGC type, can be 1 or 4
        K = getParameters(m1m4);
        G = getGatingFunctions(m1m4);
        
        ic = [V0,G.minf(V0),G.hinf(V0),G.ninf(V0),G.rinf(V0),G.finf(V0)]; % gating variables are at their steady values for V0
        [t,x] = ode23s(@(t,x) ipRGC_rhs(t,x,K,G,Iapp), [0,tf], ic);
        
        plot(t,x(:,1))
    end
    
    legend({'M1','M4'})
    xlabel('time (ms)')
    ylabel('voltage (mV)')
    title(['I_{app} = ' num2str(Iapp*1e3) 'pA'])
    
end


function xdot = ipRGC_rhs(t,x,K,G,Iapp)
    % state (t,x); parameters K; gating functions G; applied current Iapp

    V = x(1);
    m = x(2);
    h = x(3);
    n = x(4);
    r = x(5);
    f = x(6);

    xdot = NaN(size(x));
    xdot(1) = ( K.gNa*m^3*h*(K.eNa - V) + (K.gK*n^4)*(K.eK - V) + K.gCa*r*f*(K.eCa - V) + K.gL*(K.eL - V) + Iapp )/K.cM;
    xdot(2) = ( G.minf(V) - m )/(G.taum(V));
    xdot(3) = ( G.hinf(V) - h )/(G.tauh(V));
    xdot(4) = ( G.ninf(V) - n )/(G.taun(V));
    xdot(5) = ( G.rinf(V) - r )/(G.taur(V));
    xdot(6) = ( G.finf(V) - f )/(G.tauf(V));

end

function K = getParameters(m1m4)

	if m1m4==1 % M1
        K = struct('eCa',33.5,'eK',-50,'eL',-14.5,'eNa',45,'gCa',22.47,'gK',4.84,'gL',0.031,'gNa',79.18,'cM',1);
	elseif m1m4==4 % M4
        K = struct('eCa',70,'eK',-90.5722,'eL',-21.0611,'eNa',54,'gCa',21.196,'gK',4.565,'gL',.0296,'gNa',74.674,'cM',.15);
	end
    
end

function G = getGatingFunctions(m1m4)

	if m1m4==1 % M1
            
        Gparam.minfa = -0.254; Gparam.minfb = -4.4704;
        Gparam.hinfa = 0.364;  Gparam.hinfb = 11.2727;
        Gparam.ninfa = -2/17;  Gparam.ninfb = 14/17;
        Gparam.rinfa = -0.267; Gparam.rinfb = -3.3333;
        Gparam.finfa = 2/65;   Gparam.finfb = 4;

        Gparam.taumc = 0;       Gparam.tauma = -1/80;   Gparam.taumb = -143/80+log(5/21);
        Gparam.tauhc = 0.12;    Gparam.tauha = -0.28;   Gparam.tauhb = log(0.00561);
        Gparam.taunc = 0;       Gparam.tauna = -2/68;   Gparam.taunb = 67/68+log(5/21);
        Gparam.taurc = 0.738;   Gparam.taura = 0;       Gparam.taurb = -Inf;
        Gparam.taufc = 0;       Gparam.taufa = -1/110;  Gparam.taufb = log(1.79);
            
	elseif m1m4==4 % M4
            
        Gparam.minfa = -0.12383;  Gparam.minfb = -3.39282;
        Gparam.hinfa = 0.177865;  Gparam.hinfb = 9.74604;
        Gparam.ninfa = -180/3128; Gparam.ninfb = 4121/3128;
        Gparam.rinfa = -0.130435; Gparam.rinfb = -2.21376;
        Gparam.finfa = 36/2392;   Gparam.finfb = 9259/2392;

        Gparam.taumc = 0;          Gparam.tauma = -180/29440;  Gparam.taumb = -51079/29440+log(7/110);
        Gparam.tauhc = 0.0324545;  Gparam.tauha = -0.137783;   Gparam.tauhb = log(0.00490052);
        Gparam.taunc = 0;          Gparam.tauna = -45/3128;    Gparam.taunb = log(0.4821);
        Gparam.taurc = 0.1973;     Gparam.taura = 0;           Gparam.taurb = -Inf;
        Gparam.taufc = 0;          Gparam.taufa = -9/2024;     Gparam.taufb = log(0.547216);
                
	end
    
    G.minf = @(V) 1./(1 + exp(Gparam.minfa*V+Gparam.minfb));
    G.hinf = @(V) 1./(1 + exp(Gparam.hinfa*V+Gparam.hinfb));
    G.ninf = @(V) 1./(1 + exp(Gparam.ninfa*V+Gparam.ninfb)).^(1/4); % note the different form
    G.rinf = @(V) 1./(1 + exp(Gparam.rinfa*V+Gparam.rinfb));
    G.finf = @(V) 1./(1 + exp(Gparam.finfa*V+Gparam.finfb));

    G.taum = @(V) Gparam.taumc + exp(Gparam.tauma*V+Gparam.taumb);
    G.tauh = @(V) Gparam.tauhc + exp(Gparam.tauha*V+Gparam.tauhb);
    G.taun = @(V) Gparam.taunc + exp(Gparam.tauna*V+Gparam.taunb);
    G.taur = @(V) Gparam.taurc + exp(Gparam.taura*V+Gparam.taurb);
    G.tauf = @(V) Gparam.taufc + exp(Gparam.taufa*V+Gparam.taufb);
    
end