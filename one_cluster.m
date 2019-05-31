
%%%%%%% This is a STDP model as used in the paper 'Self-influencing synaptic plasticity: 
% Recurrent changes of synaptic weights can lead to specific functional properties'
% in JCNS, by M. Tamosiunaite, B. Porr, and F. Wörgötter
%
% Synaptic weigt develoment is calculated in a cluster of synapses, when the cluster has 
% three groups of inputs: correlated (svm), less correlated (svv) and uncorrelated (svd). 
% 
% First learning is driven by local dendritic spikes, later on (from time moment 'nuo_kur')
% a sharper and bigger amplitude back-propagating spike (BP) is added with
% a times-shift 'kiek_paslinkta'.
%
% Differential hebbian learning is used, with weights saturating to the interval [0,1]
% using hysteresis type function, with gradual approach and rapid go-off the zero
% or one (using function 'histereze').
% 
% Number of synapses in groups are in constants: kiek_mazu, kiek_vidutiniu, kiek_dideliu. 
%
% Intervals of dispersion of correlated, less correlated and uncorrelated input groups 
% are in constants: disp_maza, disp_didele, disp_nekoreliuotu. 
%
% Width of the BP-spike is in constant 'plotis'.
% 
% Threshold for dendritic spike generation is in constant 'thr'.
%
% Learning rate is in constant 'miu'
%%%%%%%


clear svm svv svd maxampa

%%%%%%% Initialization

kiek_mazu=3;        % Number of corelated inputs
kiek_vidutiniu=3;   % Number of less correlated inputs
kiek_dideliu=0;     % Number of uncorrelated inputs

disp_maza=5;        % Interval of distribution for correlated inputs
disp_didele=30;     % Interval of distribution of uncorrelated inputs
disp_nekoreliuotu=150; % Interval of distribution of uncorrelated inputs
plotis=0.1;         % Width parameter of the BP spike
kiek_paslinkta=10;  % Shift of the BP spike

impulsu_skaicius=600; % Number of timesteps in a simulation
nuo_kur=200;        % BP-spike emergence time moment
thr=0.14;           % Threshold for dendritic spike generation
miu=0.09;           % Learning rate

ilgis=400;          % Length of the signal we will be working with
centras=151;        % Time instance the input is centered around in the signal   

% Initializing synapses with average (0.5) weights
for j=1:kiek_mazu  
  svm(1,j)=0.5;     % svm - weights for correlated input synapses
end

for j=1:kiek_vidutiniu
  svv(1,j)=0.5;     % svv - weights for less correlated input synapses
end

for j=1:kiek_dideliu
  svd(1,j)=0.5;     % svd - weights for uncorrelated input synapses
end

%%%%%%% Learning cycle

for imp=2:impulsu_skaicius

signalas_susv=zeros(1,ilgis); % Initializing empty signal which we will fill with pulses further 

for j=1:kiek_mazu
  atsim(j)=centras+round((rand-0.5)*disp_maza);  % Generating the time moment of a correlated group input to emerge
  signalas_susv(atsim(j))=signalas_susv(atsim(j))+funh(svm(imp-1,j)); % Adding a pulse with appropriate weight at the generated time moment
                                                                      % funh preserves from weights exceeding 0 or 1 due to roundoff errors                                                               
end

for j=1:kiek_vidutiniu              
  atsiv(j)=centras+round((rand-0.5)*disp_didele); % Generating the time moment of less correlated group input to emerge
  signalas_susv(atsiv(j))=signalas_susv(atsiv(j))+funh(svv(imp-1,j)); % Adding a pulse with appropriate weight
end

for j=1:kiek_dideliu
  atsid(j)=centras+round((rand-0.5)*disp_nekoreliuotu); % Generating the time moment for uncorrelated group input to emerge
  signalas_susv(atsid(j))=signalas_susv(atsid(j))+funh(svd(imp-1,j)); % Adding a pulse with appropriate weight
end

% generating AMPA signal through filtering of the input pulses added to the signal just before 
xampa=filtras100(0.3,0.4,signalas_susv);


z=find(xampa>thr);  % time moments where AMPA signal excedes threshold
% Output spike is considered to be produced at the time moment the AMPA signal first excedes the threshold 'thr', that is at z(1) 

maxampa(imp)=max(xampa); % collecting maximas of AMPA signal for information how learning influences output development

% zero initialization of output signals
y=zeros(1,351);     % for a dendritic spike
ybp=zeros(1,501);   % for the back-propagating spike

if length(z)~=0     % If AMPA signal exceeded threshold
  centrasisej=z(1); % Adding a denfritic spike initialization pulse at z(1) time moment
  y(centrasisej)=1;
if (imp>=nuo_kur)   % If it is a time moment after BP spikes have strated emerging
  ybp(centrasisej+kiek_paslinkta)=1;  % Adding a BP spike with an appropriate shift in time
end
end

% Obtaining dendritic spike through filtering of the initiation signal
fis=0.0085;  % Filter parameters determining spike width and form
Qis=0.4;
yfd=filtras250(fis,Qis,y);


% Obtaining the BP spike through filtering of the BP initiation signal 
fis=plotis;  % Using BP width parameter here
Qis=0.4;
yfbp=filtras100(fis,Qis,ybp);
if sum(yfbp)>0
  yfbp=yfbp/sum(yfbp);  % Normalizing BP spike to have the same area independent of width
end

yf=yfd/6+yfbp*25;  % Adding the dendritic and the BP spikes, with coefficients on amplitudes providing a  
                   % BP spike of several times bigger amplitude than a dendritic spike
                
yfi=diff(yf);       % differentiation of an output for differential hebbian learning rule

%%%%%% weight update

for j=1:kiek_mazu  % correlated weight update one by one in the cycle
  signalas=zeros(1,ilgis); % generating signal with an input for the current input line
  signalas(atsim(j))=1;
  x_nmda=filtras200(0.017,0.4,signalas)/50; % filtering to provide the NMDA-signal-like shape
  integralas=sum(yfi.*x_nmda); % differential hebbian learning, where we integrate the product of input and output derivative
  svm(imp,j)=svm(imp-1,j)+histereze(svm(imp-1,j),integralas*miu); % weight update modified through hysteresis type saturation function
end

for j=1:kiek_vidutiniu % the same for less correlated weight update
  signalas=zeros(1,ilgis);
  signalas(atsiv(j))=1;
  x_nmda=filtras200(0.017,0.4,signalas)/50;
  integralas=sum(yfi.*x_nmda);
 svv(imp,j)=svv(imp-1,j)+histereze(svv(imp-1,j),integralas*miu);
end

for j=1:kiek_dideliu % the same for uncorrelated weight update
  signalas=zeros(1,ilgis);
  signalas(atsid(j))=1;
  x_nmda=filtras200(0.017,0.4,signalas)/50;
  integralas=sum(yfi.*x_nmda);
  svd(imp,j)=svd(imp-1,j)+histereze(svd(imp-1,j),integralas*miu);
end

end

%%%%%%% parameters of the model run needed for statistics like provided in appendix of the paper

vid_koreliuotu=mean(svm(impulsu_skaicius,:));   % calculating mean weight in correlated input group
vid_nekoreliuotu=mean(svv(impulsu_skaicius,:)); % calculating mean weight in less correlated group

%%%%%%% plotting weight development in the three groups, and maximum of MPA signal in each time step

figure(7)
hold off

for j=1:1:kiek_mazu
plot(svm(:,j),'k');
svm(impulsu_skaicius,j)
hold on
end

for j=1:1:kiek_vidutiniu
plot(svv(:,j),'k-.');
svv(impulsu_skaicius,j)
end

for j=1:1:kiek_dideliu
plot(svd(:,j),'k--');
svd(impulsu_skaicius,j)
end

figure(8)
hold off
plot(maxampa)



