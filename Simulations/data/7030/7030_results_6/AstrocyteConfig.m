%% Astrocyte and presynapse config
% Eero Räisänen
% 24.2.2014
% v.1.4

% ***********************************************************************
% If De Pitta simulator for presynapse is turned on the values in matrix a
% are modified during computation and the original values are stored into
% matrix w.
%************************************************************************


% micrometers to each direction(x,y,z)
% 10% conn
% 250N : [750 750 10]
% 1000N : [1400 1400 10]
CultureSpace = [750 750 10];
MinimumNeuronDistance = 10;% micrometers
MinimumAstrocyteDistance = 30;% micrometers
AVGNeuro = 0; % average length of connections
NeuroSTD = 200; % standard deviation of connection distance
ANconnectivitySTD = 150; % Standard deviation of astrocyte connectivity to neuron
AVGAstro = 0; % barbara 11.6.2018
STDAstro = 150; % barbara 11.6.2018

% Absolute maximum distance after which the astrocyte will not connect to synapses
MaxAstrocyteReachDistance = 70;

% If a new (NN and AN) topology is not made an existing topology is loaded
MakeNewTopology = 1;

if(~MakeNewTopology )
    % Path for the Astrocyte Network.
    TopologyLoadPath = '/home/genocchi/Matlab/Astro/Networks_used/ANN_30/';
    % Path for the Neuron Network.
    PremadeNNTpath = '/home/genocchi/Matlab/Astro/Networks_used/ANN_30/';
end

% If new AN topology is made it is still possible to use an existing neuronal
% network defined here.
if(MakeNewTopology)
    % Is there an old neuron network made?
    UsePremadeNeuronNetwork = 0;
    % Path for the Neuron Network.
    PremadeNNTpath = '/home/genocchi/Matlab/Astro/Networks_used/ANN_30/';
end

% is De Pitta presynapse simulator on? boolean
preSynapse = 1;

% Is astrocyte simulator on? boolean
Astrocyte = 1;

%Is astrocyte network simulator on? boolean
AstrocyteNetwork = 1;

% If astrocyte network is on, forcing astrocytes on.
if (AstrocyteNetwork == 1)
    Astrocyte = 1;
end

% If astrocyte simulator is on, forcing presynapse simulator on.
if (Astrocyte == 1)
    preSynapse = 1;
end

if (AstrocyteNetwork == 1)
    SpatialCoordinateSystem = 1;
else
    NumberOfAstrocytes = 0;
    AstrocyteInhibition = 0;
end

if (preSynapse)
    
    % Calcium at terminal bound to sensors at the beginning between 0 and 1.
    Calcium = 0;

    % Resources at each synapse at the beginning between 0 and 1.
    Resources = 1;

	% Glutamate amount at the beginning between 0 and 1. With presynapse
	% simulator only this is a constant. With Astrocyte it becomes
	% variable.
    AstrocyteGlutamate = 0.0;
    
    alpha = 0.7;

    % Regenerate resources. Every ms proportion RegenRes of used resources is
    % added to resources. 0.01 meaning that 1% of used resources is regenerated
    % each cycle. Between 0 and 1.
    RegenRes = 0.02;

    % Regenerate Calcium equilibrium in presynapse. Each ms calcium bound to
    % sensors drops to proportion CaRegen of the previous value. 0.99 meaning
    % that 99% of calcium is left compared to previous amount after 1 ms.
    % Between 0 and 1. 
    CaRegen = 0.998;% 0.998
    
    % reblock rate of NMDAR magnesium blockade. Mg left after every cycle.
    % Between 0 and 1.
    ReBlock = 0.0;
    
end

if (Astrocyte)
    
    % percentage of glutamate left after every cycle. Between 0 and 1.
    GlutamateRemoval = 0.999923;
    %GlutamateRemovalvector = [0.999923; 0.00001;0.999923];
    % Astrocyte

    % release amount by astrocyte.
    % AstrocyteReleaseAmountVector = [0.3;0.1;0.3];
    AstrocyteReleaseAmount = 0.3;
    
    % level of Ca initializing Glutamate release
    % CaGlutamateReleaseLevelVector = [0.1;0.02;0.1];
    CaGlutamateReleaseLevel = 0.1;
    
    % cumulation factor for conversion from levels IP3 to calcium. Between
    % 0 and 1.
    IP3accumulation = 0.05;
    
    % removal of IP3 from single synapse. IP3 left after each ms.
    IP3degrading = 0.85873;
    
    
end

if (AstrocyteNetwork)    
    
    NumberOfAstrocytes = 107;
    
    %"constant" variables. Can be changed but don't need to be.
    
    TimeSliceLengthInMS = 1000*t;% taking variable from config
    AnimationLengthInMS = lengthST; % taking variable from config
    VideoDirectory = directory;  % taking variable from config
    
    %Altering variables
    
    % How strong an effect synapses Ca rise has to astrocyte Ca? Higher
    % number is stonger.
    SynapseCaEffect = 5;
    %SynapseCaEffectVector = [2.0;0.2;2.0];

    %AstrocyteMatrixdimensions = [10,10,1];% x,y,z dimensions of astrocytes
    %AstrocyteAverageDistance = 70;% micrometers
    
    % are astrocyte coordinates jittered? 1 yes, 0 no
    %Jitter = 1;
    %JitterStrength = 15;%the smaller the stronger
    connectionDistance = 100;% micrometers
    
    % Astrocyte parameters 
    SpontaneousActivation = 0;%TimeSliceLengthInMS/3000;
    ActivationProbability = TimeSliceLengthInMS/1500;
    ActivationTime = TimeSliceLengthInMS/7000;
    RefractoryTime = TimeSliceLengthInMS/5000;
    slope = 0.02;
    intercept = 0.205;
    % strength of GABA inhibition(negative number)
     AstrocyteInhibitionStrength = -0.01;
    %AstrocyteInhibitionStrengthVector = [-0.3;-0.05;-0.3];
   
end

if(Astrocyte)
    AstrocyteCyclingParameter1 =  [0.0;0.0;0.0];
    AstrocyteCyclingParameter2 =  [0.0;0.0;0.0];
else
    AstrocyteCyclingParameter1 = [0.0;0.0;0.0];
    AstrocyteCyclingParameter2 = [0.0;0.0;0.0];
end
