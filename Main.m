%% This Generates the data for heat and current for the paper.(You will need mathematica to plot the final result.)
%% Loading the data and setting up the stuf!!!!
%%Final code
clc ;close all; clear all;
%% Global variables to use for the temperature:
global des act R Tmax Old_T sigma dtt R_d m_magnet ii I_power
diir = 1; % Direction Dipole
force_dir = 6; %Direction desiered Torque/force
load Tennis.mat; %The coil
load location_20100.mat; %Location of the dipole
% load Berk_location_10.mat;
[~,c_num]= size(B_tem_G);
Time_red = 0.1; % Sampeling rate for time (spacing samples, sec)
Time_diuration = 8; % Diuration of the sim (Min)
% Time_diuration = 0.1; % Diuration of the sim (Min)
Room_temperature = 27; % celsius
R = 3 * eye(c_num); % R = 3 of coils in Ohm (assuming they are identical)
I_max = 8*ones(c_num,1); % Max cont current amp (pos)
I_min = -8*ones(c_num,1); % Max cont current amp (neg)
I_power = 30; % Max count power supply  
YLIM = 250;
TempLIM = 150;
%% BOOK KIPPING
m = [0; 0; 0]; %Unit dipole
m(diir) = 1;
des = [0; 0; 0; 0; 0; 0];
des(force_dir) = 0.3;
mx = m(1);
my = m(2);
mz = m(3);
Time = Time_diuration * 60/Time_red;
field = B_tem_G;
dipol = [0 -mz my 0 0 0 0 0  ; mz  0 -mx 0 0 0 0 0 ; -my mx 0 0 0 0 0 0 ; 0 0 0 mx my mz 0 0 ; 0 0 0 0 mx 0 my  mz ; 0 0 0 -mz 0 mx -mz my ];
act = dipol * field;
% dtt= 0.05;
dtt =1;
m_magnet = magnet;
boxx = [0,0,Time_diuration,Time_diuration];
boxy = [TempLIM YLIM YLIM TempLIM];

%% PINV
Old_T = Room_temperature*ones(c_num,44); %Initializing the temperature of the coils to the room temp %Find what was 44! (I guess the numbers of layers fo wire for thermal models)
PINV_I = [];
Tmax = [];
PINV_time = [];
ii = 1;
T1 = [];
for dt = 0.0:dtt:Time
    PINV_time = [PINV_time ,dt]; 
    tic
    PINV_I = [PINV_I, PINV(des,act,3,I_max,I_min,I_power)];
    T1 = [T1,toc];
%     pinv_I(:,ii) = PINV(des,act,3,I_max,I_min,I_power);
     Old_T  = termal_model_T (PINV_I(:,ii), Room_temperature);
    Tmax = [Tmax , max(Old_T')'];
    ii = ii + 1; 
end
figure
subplot(1,5,1)
PINV_temp = Tmax;
plot(PINV_time/(60/Time_red),PINV_temp')
title('PINV','FontName', 'Times')
ylim([0 YLIM]);
patch(boxx,boxy,[1 0 0],'FaceAlpha',0.2,'LineStyle','none' )
xlabel('Time (min)','FontName', 'Times')
ylabel('Temperature (^{\circ}C)','FontName', 'Times')

%% RPINV
Old_T = Room_temperature*ones(c_num,44); %Initializing the temperature of the coils to the room temp %Find what was 44! (I guess the numbers of layers fo wire for thermal models)
RPINV_I = [];
Tmax = [];
RPINV_time = [];
ii = 1;
T2 = [];
for dt = 0.0:dtt:Time
    RPINV_time = [RPINV_time ,dt];
    tic
    RPINV_I = [RPINV_I, RPINV(des, act, 3, I_max,I_min, I_power,zeros(c_num,1), (1:c_num)')];
    T2 = [T2,toc];
    %     pinv_I(:,ii) = PINV(des,act,3,I_max,I_min,I_power);
     Old_T  = termal_model_T (RPINV_I(:,ii), Room_temperature);
    Tmax = [Tmax , max(Old_T')'];
    ii = ii + 1; 
end
subplot(1,5,2)
RPINV_temp = Tmax;
% figure
plot(RPINV_time/((60/Time_red)),RPINV_temp')
title ('RPINV','FontName', 'Times')
ylim([0 YLIM]);
patch(boxx,boxy,[1 0 0],'FaceAlpha',0.2,'LineStyle','none' )
xlabel('Time (min)','FontName', 'Times')
ylabel('Temperature (^{\circ}C)','FontName', 'Times')

%% RWPINV
Old_T = Room_temperature*ones(c_num,44); %Initializing the temperature of the coils to the room temp %Find what was 44! (I guess the numbers of layers fo wire for thermal models)
RWPINV_I = [];
Tmax = [];
RWPINV_time = [];
ii = 1;
W = eye(c_num);
T3 = [];
for dt = 0.0:dtt:Time
    RWPINV_time = [RWPINV_time ,dt]; 
    tic
    RWPINV_I = [RWPINV_I, RWPINV(des, act, 3, I_max,I_min, I_power,zeros(c_num,1), (1:c_num)',W)];
    T3 = [T3 ,toc];
%     disp('this is RWPIN')
    
    %     pinv_I(:,ii) = PINV(des,act,3,I_max,I_min,I_power);
    Old_T  = termal_model_T (RWPINV_I(:,ii), Room_temperature);
    Tmax = [Tmax , max(Old_T')'];
    for ll = 1:c_num
        W(ll,ll) = costT(Tmax(ll,ii));
    end
    ii = ii + 1; 
end
subplot(1,5,3)
RWPINV_temp = Tmax;
hold on
colorVec = {'b', 'r', 'c', 'm', 'g', 'c','g', 'b','r','c'};
linVec = {'-', '-.','--', '-.','--', '-.','-.','-','-','--'};
for i = 1:c_num
    plot(RWPINV_time/((60/Time_red)),RWPINV_temp(i,:),'Color',colorVec{i},'LineStyle',linVec{i},'LineWidth',1,'MarkerSize',1)
end
title ('RWPINV','FontName', 'Times')
ylim([0 YLIM]);
patch(boxx,boxy,[1 0 0],'FaceAlpha',0.2,'LineStyle','none' )
xlabel('Time (min)','FontName', 'Times')
ylabel('Temperature (^{\circ}C) ','FontName', 'Times')
% legend({'Coil 1','Coil 2','Coil 3','Coil 4','Coil 5','Coil 6','Coil 7','Coil 8'},'Location','northwest');



%% MIN SOLVER
Old_T = Room_temperature*ones(c_num,44); %Initializing the temperature of the coils to the room temp %Find what was 44! (I guess the numbers of layers fo wire for thermal models)
I = ones(c_num,1);
nonlcon = @power_suppy_eq_con;
obj = @Min_cost;
A = [];
b = [];
Aeq = [];
beq = [];
tic
% options=optimset('TolFun',.0001);
disp('=================Case Min_cost======================')
% objectiv = @termal_model_T_M ;
% objectiv = @curmin;
ii = 1;
T4 = [];
% MIN_time = [];
Tmax = [];
% options = optimoptions('fmincon','Algorithm','sqp','Display','off');
options = optimoptions('fmincon','Algorithm','sqp','Display','off');
% options.MaxFunctionEvaluations = 10000;
% options.FunctionTolerance = 0.00001;
% options.OptimalityTolerance = 0.1;
options.StepTolerance = 1e-04;
% options.MaxIterations= 400;

% options.StepTolerance = 1.0e-2;
for dt = 0.0:dtt:Time
%     Time
%     dt
    MIN_time(ii) = dt; 
    for i = 1:c_num
        if max(Old_T(i,:))>150
            Temp_I_max(i) = I_max(i)/(1+(exp(max(Old_T(i,:))-150))^6);
            Temp_I_min(i) = I_min(i)/(1+(exp(max(Old_T(i,:))-150))^6);
        else
            Temp_I_max(i) = I_max(i);
            Temp_I_min(i) = I_min(i);
        end
    end
    tic
    MIN_I(:,ii) = fmincon(obj,I,A,b,Aeq,beq,Temp_I_min,Temp_I_max,nonlcon,options);
    T4 = [T4,toc];
    Old_T  = termal_model_T (MIN_I(:,ii), Room_temperature);
    Tmax = [Tmax , max(Old_T')'];
    for ll = 1:c_num
        W(ll,ll) = costT(Tmax(ll,ii));
    end
    ii = ii + 1;
end
r_MIN_temp = Tmax;

subplot(1,5,4)
MIN_temp = Tmax;
hold on
colorVec = {'b', 'r', 'c', 'm', 'g', 'c','g', 'b','r','c'};
linVec = {'-', '-.','--', '-.','--', '-.','-.','-','-','--'};
for i = 1:c_num
    plot(MIN_time/((60/Time_red)),MIN_temp(i,:),'Color',colorVec{i},'LineStyle',linVec{i},'LineWidth',1,'MarkerSize',1)
end
title ('MIN','FontName', 'Times')
ylim([0 YLIM]);
patch(boxx,boxy,[1 0 0],'FaceAlpha',0.2,'LineStyle','none' )
xlabel('Time (min)','FontName', 'Times')
ylabel('Temperature (^{\circ}C) ','FontName', 'Times')

%% Comparing the performance 

subplot(1,5,5)
PINV_P = act*PINV_I;
RPINV_P = act*RPINV_I;
RWPINV_P = act*RWPINV_I;
MIN_P = act*MIN_I;
LL = Low_filter_force(MIN_P(force_dir,:));
disp('The current used')
curmin(MIN_I(:,1))
plot(PINV_time/((60/Time_red)),PINV_P(force_dir,:),'r*'); 
hold on
plot(RPINV_time/((60/Time_red)),RPINV_P(force_dir,:),'bo'); 
plot(RWPINV_time/((60/Time_red)),RWPINV_P(force_dir,:),'.c');
plot(MIN_time/((60/Time_red)),MIN_P(force_dir,:));
plot(MIN_time/((60/Time_red)),LL,'r');


title (strcat('The requested output was: ',num2str(des(force_dir))),'FontName', 'Times')
% ylim([0.5*des(force_dir) , 1.2*des(force_dir)]);
xlabel('Time (min)','FontName', 'Times')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Saving the .CSV file for mathematica plots (PINV)
% D = PINV_temp;
% t = PINV_time;
% [N_line,N_data] = size(D); 
% C = [];
% for l = 1:N_line
%     for d = 1:N_data
%         C = [C,l,t(d),D(l,d)];
%     end
% end
% csvwrite('PINV_temp.csv',C);
% D = PINV_I;
% [N_line,N_data] = size(D) ;
% C = [];
% for l = 1:N_line
%     for d = 1:N_data
%         C = [C,l,t(d),D(l,d)];
%     end
% end
% csvwrite('PINV_I.csv',C);
% 
% %% Saving the .CSV file for mathematica plots (RPINV)
% D = RPINV_temp;
% t = RPINV_time;
% [N_line,N_data] = size(D); 
% C = [];
% for l = 1:N_line
%     for d = 1:N_data
%         C = [C,l,t(d),D(l,d)];
%     end
% end
% csvwrite('RPINV_temp.csv',C);
% D = RPINV_I;
% [N_line,N_data] = size(D); 
% C = [];
% for l = 1:N_line
%     for d = 1:N_data
%         C = [C,l,t(d),D(l,d)];
%     end
% end
% csvwrite('RPINV_I.csv',C);
% 
% %% Saving the .CSV file for mathematica plots (RWPINV)
% D = RWPINV_temp;
% t = RWPINV_time;
% [N_line,N_data] = size(D); 
% C = [];
% for l = 1:N_line
%     for d = 1:N_data
%         C = [C,l,t(d),D(l,d)];
%     end
% end
% csvwrite('RWPINV_temp.csv',C);
% D = RWPINV_I;
% [N_line,N_data] = size(D); 
% C = [];
% for l = 1:N_line
%     for d = 1:N_data
%         C = [C,l,t(d),D(l,d)];
%     end
% end
% csvwrite('RWPINV_I.csv',C);
%%
RPINVF = RPINV_P(6,1);
PINVF = PINV_P(6,1);
(RPINVF-PINVF)/PINVF
PINV_I(:,1);
mean(T1)
mean(T3)
mean(T4)
% pinv(act)*(RPINV_P(:,1))*2.5/8
% PINV_I(:,1)*2.5/8
%%
disp("PINV:")
sum(abs(PINV_P(1:3,end)-[0,0,0]'))
rad2deg(subspace(PINV_P(4:6,end),[0,0,6]'))
disp("RWPINV:")
sum(abs(RWPINV_P(1:3,end)-[0,0,0]'))
rad2deg(subspace(RWPINV_P(4:6,end),[0,0,6]'))
disp("sqo:")
sum(abs(MIN_P(1:3,end)-[0,0,0]'))
rad2deg(subspace(MIN_P(4:6,end),[0,0,6]'))
